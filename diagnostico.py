"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: diagnostico.py — Fase 1 del orden de trabajo de calibración

Consume los volcados de telemetría del Hilo Rápido y emite un reporte de
consistencia del filtro. Se corre OFFLINE, después de cada sesión de Testnet;
nunca dentro del lazo de control.

    python diagnostico.py [directorio_telemetria]

CONTENIDO (Fase 1 de ORDEN_TRABAJO_CALIBRACION_1.1.md)
------------------------------------------------------
  1.1  NIS      — media y banda de confianza chi-cuadrado.
  1.3  Ljung-Box — test de blancura. ES LA COMPUERTA DE ENTRADA A LA FASE 2.
  1.4  Shapiro-Wilk — normalidad de los residuales; decide si ALS-IRLS vale la pena.

POR QUÉ EL ORDEN IMPORTA
------------------------
Diagnosticar -> calibrar -> meta-optimizar, sin saltos. ALS (Fase 2) es
matemáticamente válido SÓLO si A y H son correctas y las covarianzas de ruido son
las únicas incógnitas. Si A tiene estructura no modelada, ALS la absorberá dentro
de Q y devolverá una respuesta confiadamente equivocada. Ljung-Box es
precisamente lo que detecta ese caso, y por eso se ejecuta primero.

Los mensajes de este módulo van en ASCII: la consola de Windows es cp1252 y no
puede codificar letras griegas.
"""

import os
import sys
import glob

import numpy as np

DIR_TELEMETRIA_POR_DEFECTO = "telemetria"
DIM_MEDICION = 2  # m = rango del vector de medicion z_k (Sec. 7.3.1)

# A.2 de la v1.2: umbral de rho_1 por debajo del cual un rechazo de Ljung-Box se
# considera estadisticamente cierto pero SIN RELEVANCIA PRACTICA. La tabla de
# decision de la Seccion D usa 0.2 como frontera de "rho grande" (model mismatch)
# y 0.05 como "rho pequeno" (ruido mal calibrado).
UMBRAL_RHO_MATERIAL = 0.20
UMBRAL_RHO_DESPRECIABLE = 0.05


# ==============================================================================
# CARGA DE VOLCADOS
# ==============================================================================
def cargar_telemetria(directorio: str = DIR_TELEMETRIA_POR_DEFECTO) -> np.ndarray:
    """Concatena todos los bloques volcados y los ordena cronologicamente."""
    rutas = sorted(glob.glob(os.path.join(directorio, "telem_*.npy")))
    rutas_pq = sorted(glob.glob(os.path.join(directorio, "telem_*.parquet")))

    bloques = [np.load(r) for r in rutas]
    if rutas_pq:
        try:
            import pyarrow.parquet as pq

            for r in rutas_pq:
                tabla = pq.read_table(r)
                cols = {n: np.asarray(tabla[n]) for n in tabla.column_names}
                rec = np.zeros(
                    len(next(iter(cols.values()))),
                    dtype=[(n, np.float64) for n in cols],
                )
                for n, v in cols.items():
                    rec[n] = v
                bloques.append(rec)
        except ImportError:
            print("[AVISO] hay .parquet pero pyarrow no esta instalado; se omiten.")

    if not bloques:
        raise FileNotFoundError(
            f"No hay volcados de telemetria en '{directorio}'. "
            f"Corre el bot el tiempo suficiente para llenar al menos un bloque, "
            f"o cierralo con Ctrl-C (el cierre volca el bloque parcial)."
        )

    # El directorio acumula corridas de distintas versiones del TELEM_DTYPE, y
    # `np.concatenate` sobre dtypes estructurados distintos revienta con
    # DTypePromotionError. Se conserva el esquema MAS COMPLETO y se descartan los
    # bloques viejos, avisando: mezclarlos rellenando campos ausentes con ceros
    # seria peor -- un `hay_medicion` inventado a cero borraria todas las
    # observaciones de esos bloques sin que nadie se enterara.
    if len({b.dtype.names for b in bloques}) > 1:
        campos_max = max((b.dtype.names for b in bloques), key=len)
        completos = [b for b in bloques if b.dtype.names == campos_max]
        print(
            f"[AVISO] {len(bloques) - len(completos)} bloque(s) de telemetria con "
            f"un esquema anterior; se omiten. Se analizan {len(completos)} bloque(s) "
            f"con los {len(campos_max)} campos actuales."
        )
        bloques = completos

    datos = np.concatenate(bloques)
    # Se descartan las filas nunca escritas (t_wall == 0) del ultimo bloque parcial.
    datos = datos[datos["t_wall"] > 0]
    return datos[np.argsort(datos["t_wall"])]


def solo_observaciones(datos: np.ndarray) -> np.ndarray:
    """Filas en las que REALMENTE hubo una medicion. v1.3.

    Desde que el filtro corrige una vez por PAQUETE NUEVO y no una vez por ciclo
    de control (Sec. 7.3, divergencia documentada), la telemetria sigue teniendo
    una fila por ciclo (~90 Hz) pero la innovacion solo existe cuando llego un
    paquete. En los demas ciclos y0 e y1 valen cero POR RELLENO.

    ⚠ TODO analisis basado en la innovacion —NIS, Ljung-Box, Shapiro-Wilk, el A/B
    de la Seccion E— debe correr sobre este subconjunto. Sin filtrar, la serie es
    ~99 % ceros: rho_1 sale 0.0000, Shapiro rechaza normalidad por la masa
    puntual en cero, y el veredicto tiene apariencia de dato sin serlo.
    Medido en la corrida del 2026-08-04: 178 observaciones reales entre 20 000
    filas de telemetria.

    Compatible con volcados anteriores a la v1.3, que no llevan la bandera: en
    ese caso todas las filas tenian medicion por construccion.
    """
    if "hay_medicion" in datos.dtype.names:
        return datos[datos["hay_medicion"].astype(int) == 1]
    return datos


# ==============================================================================
# TRANSITORIO DE ARRANQUE
# ==============================================================================
def detectar_fin_transitorio(
    datos: np.ndarray, m: int = DIM_MEDICION, corridas: int = 10
) -> int:
    """Primer indice a partir del cual el filtro ya salio del arranque en frio.

    POR QUE HACE FALTA: el EAKF arranca con P_v, P_Rn ~ 1e6 y con x2 = 0 mientras
    la medicion del residuo macro vale ~45 000 (Sec. 7.3.4: incertidumbre inicial
    masiva para forzar la asimilacion). La innovacion del primer ciclo es por
    tanto enorme y perfectamente legitima. Si no se descarta:
      - la media del NIS queda dominada por unas pocas muestras (se observaron
        medias de ~3300 contra un teorico de 2),
      - la curtosis se dispara a ~3500 y Shapiro-Wilk rechaza normalidad por un
        unico outlier,
      - y Ljung-Box detecta "autocorrelacion" que en realidad es el escalon del
        arranque, no estructura no modelada.
    Es decir: sin este recorte los tres tests de la Fase 1 dan falsos positivos y
    la compuerta de la Fase 2 se cierra por la razon equivocada.

    Criterio: primer indice desde el cual `corridas` muestras consecutivas del NIS
    caen dentro de la banda de una sola muestra al 97.5%.
    """
    from scipy.stats import chi2

    if "nis" not in datos.dtype.names:
        return 0
    nis = datos["nis"]
    hi1 = float(chi2.ppf(0.975, m))
    ok = np.isfinite(nis) & (nis <= hi1)
    # Primer arranque de `corridas` valores consecutivos dentro de banda.
    racha = 0
    for i, v in enumerate(ok):
        racha = racha + 1 if v else 0
        if racha >= corridas:
            return i - corridas + 1
    return 0  # Nunca se estabilizo: se reporta sobre todo el registro


# ==============================================================================
# 1.1 NIS — Secuencia de Innovacion Normalizada al Cuadrado
# ==============================================================================
def reporte_nis(nis: np.ndarray, m: int = DIM_MEDICION, confianza: float = 0.95):
    """Media del NIS contra su banda chi-cuadrado.

    Si el filtro es consistente, eps_k ~ chi2 con m g.l., luego E[eps_k] = m.
      - Media persistentemente POR DEBAJO de m -> filtro conservador: sobreestima
        su propia incertidumbre y desaprovecha informacion.
      - POR ENCIMA de m -> sobreconfiado: divergencia inminente.
    """
    from scipy.stats import chi2

    v = nis[np.isfinite(nis)]
    if v.size < 2:
        return {"n": int(v.size), "media": float("nan"), "veredicto": "SIN DATOS"}

    n = v.size
    alfa = 1.0 - confianza
    lo = float(chi2.ppf(alfa / 2.0, n * m) / n)
    hi = float(chi2.ppf(1.0 - alfa / 2.0, n * m) / n)
    media = float(v.mean())

    # Fraccion de muestras individuales dentro de la banda de una sola muestra.
    lo1 = float(chi2.ppf(alfa / 2.0, m))
    hi1 = float(chi2.ppf(1.0 - alfa / 2.0, m))
    dentro = float(np.mean((v >= lo1) & (v <= hi1)))

    if media < lo:
        veredicto = "CONSERVADOR (sobreestima su incertidumbre)"
    elif media > hi:
        veredicto = "SOBRECONFIADO (riesgo de divergencia)"
    else:
        veredicto = "CONSISTENTE"

    return {
        "n": n, "media": media, "teorico": float(m), "banda": (lo, hi),
        "banda_muestral": (lo1, hi1), "fraccion_dentro": dentro,
        "veredicto": veredicto,
    }


# ==============================================================================
# 1.3 LJUNG-BOX — compuerta de entrada a la Fase 2
# ==============================================================================
def _rho_minimo_detectable(n: int, h: int = 20, alfa: float = 0.05) -> float:
    """rho_1 minimo que basta para rechazar Ljung-Box con este n (A.2 de la v1.2).

    Suponiendo que toda la potencia del estadistico viene del primer rezago:
        Q ~= n(n+2)·rho_1²/(n-1) >= chi2_critico   =>   rho_1 >= sqrt(...)
    Sirve para leer el veredicto en contexto: el mismo p-valor significa cosas
    muy distintas con n=500 que con n=50 000.
    """
    from scipy.stats import chi2

    if n <= h + 1:
        return float("nan")
    critico = float(chi2.ppf(1.0 - alfa, h))
    return float(np.sqrt(critico * (n - 1) / (n * (n + 2))))


def ljung_box(serie: np.ndarray, h: int = 20):
    """Ljung-Box univariante sobre una componente de la innovacion.

        Q_LB = n(n+2) * sum_{j=1..h} rho_j^2 / (n-j)      ~  chi2 con h g.l.

    rho_j es la autocorrelacion muestral al rezago j. Es INVARIANTE DE ESCALA, asi
    que no hace falta haber logeado S_k para normalizar la innovacion: basta con
    la serie cruda de ỹ_k.
    """
    from scipy.stats import chi2

    x = serie[np.isfinite(serie)]
    n = x.size
    if n <= h + 1:
        return {"n": int(n), "Q": float("nan"), "veredicto": "SIN DATOS"}

    x = x - x.mean()
    denom = float(np.dot(x, x))
    if denom <= 0.0:
        return {"n": int(n), "Q": float("nan"), "veredicto": "SERIE CONSTANTE"}

    Q = 0.0
    rhos = []
    for j in range(1, h + 1):
        rho_j = float(np.dot(x[j:], x[:-j]) / denom)
        rhos.append(rho_j)
        Q += rho_j * rho_j / (n - j)
    Q *= n * (n + 2)

    p = float(chi2.sf(Q, h))
    critico = float(chi2.ppf(0.95, h))

    # A.2 de la v1.2: EL VEREDICTO SE LEE SOBRE LA MAGNITUD, NO SOBRE EL P-VALOR.
    # Ljung-Box es extremadamente potente con n grande y rechaza por correlaciones
    # sin ninguna relevancia práctica. Con h=20 y alfa=0.05, el rho_1 minimo que
    # basta para rechazar es:
    #     n=500 -> 0.250 | n=2 000 -> 0.125 | n=11 500 -> 0.052 | n=50 000 -> 0.025
    # Con rho_1 ~ 0.8 no hay ambiguedad; con rho_1 ~ 0.03 el p-valor diria lo mismo
    # y significaria algo completamente distinto.
    rho1 = abs(rhos[0]) if rhos else 0.0
    if p >= 0.05:
        veredicto = "BLANCA"
    elif rho1 >= UMBRAL_RHO_MATERIAL:
        veredicto = "RECHAZA BLANCURA"
    else:
        veredicto = "rechazo NO MATERIAL"

    return {
        "n": int(n), "h": h, "Q": float(Q), "critico_95": critico, "p": p,
        "rho": rhos[:3], "rho1": rho1, "veredicto": veredicto,
    }


def ljung_box_multivariante(Y: np.ndarray, h: int = 20):
    """Portmanteau multivariante de Hosking sobre el vector de innovacion.

        Q = n^2 * sum_{j=1..h} tr(C_j' C_0^-1 C_j C_0^-1) / (n-j)   ~ chi2, h*m^2 g.l.

    Es la generalizacion natural del estadistico univariante de la Sec. 1.3 al
    caso m > 1: contempla la correlacion CRUZADA entre las componentes (precio y
    residuo macro), que el test por componente no ve.
    """
    from scipy.stats import chi2

    Y = Y[np.all(np.isfinite(Y), axis=1)]
    n, m = Y.shape
    if n <= h + 1:
        return {"n": int(n), "Q": float("nan"), "veredicto": "SIN DATOS"}

    Yc = Y - Y.mean(axis=0)
    C0 = (Yc.T @ Yc) / n
    try:
        C0_inv = np.linalg.inv(C0)
    except np.linalg.LinAlgError:
        return {"n": int(n), "Q": float("nan"), "veredicto": "C0 SINGULAR"}

    Q = 0.0
    for j in range(1, h + 1):
        Cj = (Yc[j:].T @ Yc[:-j]) / n
        Q += float(np.trace(Cj.T @ C0_inv @ Cj @ C0_inv)) / (n - j)
    Q *= n * n

    gl = h * m * m
    p = float(chi2.sf(Q, gl))
    return {
        "n": int(n), "h": h, "gl": gl, "Q": float(Q),
        "critico_95": float(chi2.ppf(0.95, gl)), "p": p,
        "veredicto": "BLANCA" if p >= 0.05 else "RECHAZA BLANCURA",
    }


# ==============================================================================
# 1.4 SHAPIRO-WILK
# ==============================================================================
def shapiro_normalidad(serie: np.ndarray, max_muestras: int = 5000, semilla: int = 0):
    """Normalidad de los residuales. Decide si ALS-IRLS (Fase 2.3) vale la pena.

    Colas pesadas o asimetria justifican sustituir el criterio cuadratico de ALS
    por la funcion de costo de Huber. No es un fin en si mismo.
    """
    from scipy.stats import shapiro, kurtosis, skew

    x = serie[np.isfinite(serie)]
    if x.size < 20:
        return {"n": int(x.size), "veredicto": "SIN DATOS"}

    # scipy.stats.shapiro pierde fiabilidad por encima de ~5000 muestras.
    if x.size > max_muestras:
        rng = np.random.default_rng(semilla)
        x_test = rng.choice(x, size=max_muestras, replace=False)
    else:
        x_test = x

    W, p = shapiro(x_test)
    return {
        "n": int(x.size), "n_test": int(x_test.size), "W": float(W), "p": float(p),
        "curtosis_exceso": float(kurtosis(x, fisher=True)),
        "asimetria": float(skew(x)),
        "veredicto": "NORMAL" if p >= 0.05 else "RECHAZA NORMALIDAD",
    }


# ==============================================================================
# REPORTE
# ==============================================================================
def generar_reporte(datos: np.ndarray, h: int = 20, descartar: int = -1) -> str:
    L = []
    a = L.append
    ciclos = len(datos)
    dur = float(datos["t_wall"].max() - datos["t_wall"].min())

    # v1.3: TODO lo que sigue se calcula sobre las filas con medicion real. Ver
    # `solo_observaciones`: con el filtro corrigiendo por paquete nuevo, la serie
    # sin filtrar es mayoritariamente ceros de relleno y los tres tests de la
    # Fase 1 dan resultados con apariencia de dato que no lo son.
    datos = solo_observaciones(datos)
    bruto = len(datos)
    if bruto < 50:
        return (
            "=" * 78
            + f"\nREPORTE ABORTADO: solo {bruto} observaciones reales en {ciclos} "
            f"ciclos de control.\nEl feed no esta entregando datos, o esta "
            f"entregando muy por debajo de la cadencia\ndel Hilo Rapido. Revisar "
            f"el estado del WebSocket antes de sacar conclusiones.\n" + "=" * 78
        )

    if descartar < 0:
        descartar = detectar_fin_transitorio(datos)
    if descartar > 0 and descartar < len(datos) - 100:
        datos = datos[descartar:]
    elif descartar >= len(datos) - 100:
        descartar = 0  # No queda registro util: se reporta sobre todo

    a("=" * 78)
    a("REPORTE DE CONSISTENCIA DEL FILTRO - Fase 1 del orden de trabajo")
    a("=" * 78)
    # v2.0: desde que el filtro corre bajo Delta_n = 1, cada fila de telemetria
    # es un TICK (una transaccion), no un ciclo de control. La etiqueta anterior
    # decia "ciclos" y era enganosa: sugeria que el 100 % de observaciones era
    # sospechoso, cuando es justo lo que se persigue.
    a(f"TICKS registrados: {ciclos}   Duracion: {dur:.1f} s   "
      f"Tasa: {ciclos/max(1e-9,dur):.2f} ticks/s")
    a(f"OBSERVACIONES REALES: {bruto} ({bruto/max(1,ciclos):.2%} de los ticks), "
      f"o sea {bruto/max(1e-9,dur):.2f} mediciones/s.")
    if bruto / max(1, ciclos) < 0.05:
        a("  [!] Menos del 5% de los ticks traen medicion. El filtro PREDICE mucho")
        a("      y CORRIGE poco, asi que P crece bastante entre observaciones y el")
        a("      NIS saldra bajo por construccion. No es un filtro conservador: es")
        a("      un feed lento. Mirar el estado del WebSocket antes que las matrices.")
    if descartar:
        a(f"Transitorio de arranque descartado: {descartar} muestras "
          f"({descartar/bruto:.2%}). Analizando {len(datos)}.")
        a("  (arranque en frio de la Sec. 7.3.4: P_v,P_Rn ~ 1e6 y x2=0 contra una")
        a("   medicion de ~45 000 producen una innovacion inicial legitima pero")
        a("   enorme, que sin recortar falsea los tres tests.)")
        if descartar / bruto > 0.10:
            a("  [!] El recorte supera el 10% del registro. Un transitorio real dura")
            a("      unas pocas decenas de ciclos; que haga falta descartar tanto")
            a("      significa que el NIS rara vez entra en banda, lo cual YA ES el")
            a("      hallazgo: el filtro no es consistente en regimen, no solo al")
            a("      arrancar. Contrastar con el veredicto del NIS de abajo.")
    else:
        a("Sin recorte de transitorio.")
    a("")

    # ---- 1.1 NIS ----
    a("-" * 78)
    a("1.1  NIS (Secuencia de Innovacion Normalizada al Cuadrado)")
    a("-" * 78)
    if "nis" in datos.dtype.names:
        r = reporte_nis(datos["nis"])
        if r["veredicto"] == "SIN DATOS":
            a("  SIN DATOS (ningun ciclo con medicion valida).")
        else:
            a(f"  media = {r['media']:.4f}   teorico = {r['teorico']:.1f}   "
              f"banda 95% de la media = [{r['banda'][0]:.4f}, {r['banda'][1]:.4f}]")
            a(f"  muestras individuales dentro de [{r['banda_muestral'][0]:.3f}, "
              f"{r['banda_muestral'][1]:.3f}]: {r['fraccion_dentro']:.1%} "
              f"(esperado ~95%)")
            a(f"  VEREDICTO: {r['veredicto']}")
    else:
        a("  El volcado no tiene columna 'nis' (telemetria anterior a la Fase 1).")

    # ---- Contraste con el criterio legacy ----
    if "dtr_dt" in datos.dtype.names:
        d = datos["dtr_dt"][np.isfinite(datos["dtr_dt"])]
        if d.size:
            a("")
            a(f"  [contraste] criterio legacy dTr/dt (Sec. 7.1): mediana="
              f"{np.median(d):.4f}  p99={np.percentile(d,99):.4f}  max={d.max():.4f}")
            a("  Recordatorio: dTr/dt depende de dt, de rho_k y del regimen; el NIS")
            a("  no. Por eso el NIS es el criterio primario (Fase 1.2).")

    # ---- 1.3 Ljung-Box ----
    a("")
    a("-" * 78)
    a("1.3  LJUNG-BOX (blancura de la innovacion) -- COMPUERTA A LA FASE 2")
    a("-" * 78)
    a(f"  n efectivo tras el recorte: {len(datos)}   rezagos h={h}")
    a("  (con este n, el rho_1 minimo que basta para rechazar es "
      f"{_rho_minimo_detectable(len(datos), h):.3f})")
    a("")
    Y = np.column_stack([datos["y0"], datos["y1"]])
    resultados = {}
    for etiqueta, serie in (("y0 (precio)", datos["y0"]), ("y1 (residuo)", datos["y1"])):
        r = ljung_box(serie, h)
        resultados[etiqueta] = r
        if np.isfinite(r["Q"]):
            rr = r.get("rho", [])
            a(f"  {etiqueta:<14} Q={r['Q']:11.2f}  p={r['p']:.2e}   "
              f"rho1={rr[0]:+.4f} rho2={rr[1]:+.4f} rho3={rr[2]:+.4f}"
              f"   -> {r['veredicto']}")
        else:
            a(f"  {etiqueta:<14} {r['veredicto']}")

    rm = ljung_box_multivariante(Y, h)
    resultados["multivariante"] = rm
    if np.isfinite(rm["Q"]):
        a(f"  {'multivariante':<14} Q={rm['Q']:11.2f}  p={rm['p']:.2e}   "
          f"critico_95(gl={rm['gl']})={rm['critico_95']:.1f}"
          f"   -> {rm['veredicto']}")
    else:
        a(f"  multivariante  {rm['veredicto']}")

    # Veredicto por MAGNITUD (A.2), no por p-valor.
    rho1_max = max(
        (r.get("rho1", 0.0) for r in resultados.values() if "rho1" in r), default=0.0
    )
    rechaza_material = rho1_max >= UMBRAL_RHO_MATERIAL
    rechaza_p = any(
        r.get("p", 1.0) < 0.05 for r in resultados.values() if "p" in r
    )
    # --- Multi-tasa: distinguir retención de orden cero de mismatch real -----
    a("")
    a("  Estructura de la autocorrelacion (corta vs persistente):")
    for etiqueta, serie in (("y0 (precio)", datos["y0"]), ("y1 (residuo)", datos["y1"])):
        v = serie[np.isfinite(serie)]
        if v.size < 100:
            continue
        frac = float(np.mean(np.diff(v) != 0.0))
        repite = 1.0 / frac if frac > 0 else float("inf")
        x = v - v.mean()
        den = float(np.dot(x, x))
        rr = []
        for j in (1, 10, 45):
            rr.append(float(np.dot(x[j:], x[:-j]) / den) if den > 0 else 0.0)
        a(f"    {etiqueta:<14} cambia en {frac:5.1%} de los pasos "
          f"(se repite ~{repite:.1f}x)   rho1={rr[0]:+.3f} rho10={rr[1]:+.3f} "
          f"rho45={rr[2]:+.3f}")
        if repite > 3.0:
            a(f"      [!] RETENCION DE ORDEN CERO. Esta componente llega a una tasa "
              f"~{repite:.0f}x menor")
            a("          que el lazo de control, y el filtro la asimila como si fuera una")
            a("          medicion NUEVA en cada ciclo. Eso fabrica autocorrelacion por si")
            a("          solo: su rho_1 NO es evidencia sobre la matriz A. Lo correcto es")
            a("          actualizar esa componente solo cuando llega un dato fresco")
            a("          (actualizacion secuencial multi-tasa), no en cada tick.")
    a("")
    a("  Lectura: una autocorrelacion que decae a ~0 en pocas decenas de rezagos es de")
    a("  CORTO ALCANCE (retencion o error de modelo suave). Un mismatch estructural de")
    a("  la Sec. E -- A asume velocidad constante y el mercado es oscilatorio -- deja")
    a("  correlacion a rezagos del orden del ciclo estructural, no solo en rho_1.")
    a("")
    a(f"  max|rho_1| sobre las componentes = {rho1_max:.4f}")
    if rechaza_material:
        a("  >> FASE 2 BLOQUEADA. rho_1 es GRANDE (>= "
          f"{UMBRAL_RHO_MATERIAL:.2f}): hay predictibilidad remanente que la")
        a("     matriz cinematica A no capturo. Eso es MODEL MISMATCH, no mala")
        a("     calibracion de ruido. Recalibrar Q y R aqui seria tapar el sintoma:")
        a("     ALS absorberia el error de modelo dentro de Q y devolveria una")
        a("     respuesta confiadamente equivocada.")
        a("     Respuesta correcta: AUMENTO DE ESTADO (Sec. 3.6 del PDF).")
        a("     PROHIBIDO: inflar Q o R con una constante ad-hoc. Medido sobre este")
        a("     sistema, Q x7 lleva el NIS de 19.35 a 3.22 pero solo baja rho_1 de")
        a("     +0.79 a +0.57, y Q x50 lo deja en +0.44. Un error determinista no lo")
        a("     describe ninguna matriz de covarianza: la inflacion no corrige, oculta.")
    elif rechaza_p:
        a(f"  >> Rechazo estadistico pero NO MATERIAL (rho_1 < {UMBRAL_RHO_MATERIAL:.2f}).")
        a("     Con n grande Ljung-Box rechaza por correlaciones irrelevantes. Si el")
        a("     NIS esta fuera de banda, la causa es ruido mal calibrado y el")
        a("     estimador correcto es ALS (Fase 2), no una constante a ojo.")
    else:
        a("  >> Compuerta superada: la innovacion pasa por blanca. Fase 2 habilitada.")

    # ---- 1.4 Shapiro-Wilk ----
    a("")
    a("-" * 78)
    a("1.4  SHAPIRO-WILK (normalidad de los residuales)")
    a("-" * 78)
    pesadas = False
    for etiqueta, serie in (("y0 (precio)", datos["y0"]), ("y1 (residuo)", datos["y1"])):
        r = shapiro_normalidad(serie)
        if r["veredicto"] == "SIN DATOS":
            a(f"  {etiqueta:<14} SIN DATOS")
            continue
        a(f"  {etiqueta:<14} W={r['W']:.5f}  p={r['p']:.3e}  "
          f"curtosis_exceso={r['curtosis_exceso']:+.3f}  "
          f"asimetria={r['asimetria']:+.3f}  -> {r['veredicto']}")
        if r["veredicto"] != "NORMAL" and r["curtosis_exceso"] > 1.0:
            pesadas = True
    a("")
    if pesadas:
        a("  >> Colas pesadas detectadas. Si se llega a la Fase 2, considerar")
        a("     ALS-IRLS (Huber) -- PERO correr ALS estandar primero: sin baseline")
        a("     no se puede medir la mejora.")
    else:
        a("  >> Sin evidencia fuerte de colas pesadas: ALS estandar deberia bastar.")

    # ---- SECCIÓN E de la v1.3: reporte comparativo del A/B ----
    a("")
    a(reporte_ab(datos, h=h))

    a("")
    a("=" * 78)
    return "\n".join(L)


# ==============================================================================
# SECCIÓN E de ORDEN_TRABAJO_RIESGO_1_3 — protocolo A/B y regla de decisión
# ==============================================================================
def reporte_ab(datos: np.ndarray, h: int = 20) -> str:
    """Compara el EAKF de control contra el SOMBRA sobre el mismo flujo.

    Regla de decisión de la Sec. E.3, sobre las innovaciones del canal `y0`, con
    recorte de transitorio y leyendo la MEDIANA del NIS:

      | rho_1(y0) del armonico < 0.20 Y menor que el de velocidad constante
      |      -> ADOPTAR A_arm
      | ambos con rho_1 >= 0.20
      |      -> no adoptar; el mismatch no es el que creiamos
      | armonico peor
      |      -> no adoptar, y VERIFICAR ANTES D.1 y D.2

    ⚠ NO SE USA `y1` PARA ESTA COMPUERTA. Su rho_1 = 0.87 es artefacto del
    solapamiento de ventanas del EMD (retencion de orden cero: el Hilo Lento
    publica a 2 Hz y el Rapido corre a ~89 Hz) y no es evidencia sobre `A`. Esto
    resuelve la pregunta que la v1.2 dejo abierta.

    ⚠ SE LEE LA MEDIANA DEL NIS, NO LA MEDIA. Medido en la v1.2: media 147 contra
    mediana 4.69, porque la cola es pesada.
    """
    L = []
    a = L.append
    a("=" * 78)
    a("SECCION E (v1.3)  PROTOCOLO A/B: A_arm CONTRA VELOCIDAD CONSTANTE")
    a("=" * 78)

    if "rama_A" not in datos.dtype.names or "y0_sombra" not in datos.dtype.names:
        a("  Telemetria anterior a la v1.3: sin campos rama_A / *_sombra.")
        a("  Volver a correr el bot para poder decidir sobre A_arm.")
        return "\n".join(L)

    # Solo filas con medicion: la innovacion del sombra tambien vale cero por
    # relleno en los ciclos sin paquete nuevo.
    datos = solo_observaciones(datos)
    if len(datos) < 50:
        a(f"  Solo {len(datos)} observaciones reales: insuficiente para decidir.")
        return "\n".join(L)

    rama = datos["rama_A"].astype(int)
    n_arm = int((rama == 1).sum())
    a(f"  Muestras: {len(datos)} | rama de control armonica en {n_arm} "
      f"({n_arm/max(1,len(datos)):.1%}), velocidad constante en {len(datos)-n_arm}")
    a("")

    # El sombra usa SIEMPRE la rama contraria a la de control, asi que las
    # innovaciones de cada modelo hay que reensamblarlas cruzando por `rama_A`.
    y0_ctrl, y0_som = datos["y0"], datos["y0_sombra"]
    y0_arm = np.where(rama == 1, y0_ctrl, y0_som)
    y0_vc = np.where(rama == 1, y0_som, y0_ctrl)
    nis_ctrl, nis_som = datos["nis"], datos["nis_sombra"]
    nis_arm = np.where(rama == 1, nis_ctrl, nis_som)
    nis_vc = np.where(rama == 1, nis_som, nis_ctrl)

    filas = []
    for etiqueta, y0s, niss in (
        ("velocidad constante", y0_vc, nis_vc),
        ("oscilador armonico", y0_arm, nis_arm),
    ):
        lb = ljung_box(y0s, h=h)
        if "rho1" not in lb:
            # `ljung_box` devuelve un dict reducido cuando n <= h+1 o la serie es
            # constante. Sin esta guarda el reporte reventaba con KeyError justo
            # en el caso que mas interesa diagnosticar: pocas observaciones.
            a(f"  {etiqueta:<22} {lb.get('veredicto', 'SIN DATOS')} (n={lb['n']})")
            a("  Insuficiente para decidir sobre A_arm.")
            return "\n".join(L)
        finitos = niss[np.isfinite(niss)]
        mediana = float(np.median(finitos)) if finitos.size else float("nan")
        media = float(np.mean(finitos)) if finitos.size else float("nan")
        filas.append((etiqueta, lb, mediana, media))
        a(f"  {etiqueta:<22} NIS mediana={mediana:8.3f}  (media={media:10.3f}, "
          f"inflada por cola pesada)")
        a(f"  {'':<22} rho_1={lb['rho1']:+.4f}  rho_2={lb['rho'][1]:+.4f}  "
          f"rho_3={lb['rho'][2]:+.4f}  p={lb['p']:.3e}")
        a(f"  {'':<22} n efectivo={lb['n']}  rho_1 minimo detectable="
          f"{_rho_minimo_detectable(lb['n'], h):.4f}")
        a("")

    rho_vc = abs(filas[0][1]["rho1"])
    rho_arm = abs(filas[1][1]["rho1"])
    n_obs = filas[0][1]["n"]
    rho_min = _rho_minimo_detectable(n_obs, h)
    a("-" * 78)
    a("REGLA DE DECISION (Sec. E.3), sobre y0 -- NUNCA sobre y1")
    a("-" * 78)

    # --- Condiciones bajo las que el A/B NO es decidible --------------------
    # Se comprueban ANTES de aplicar la regla, porque un veredicto emitido sobre
    # datos que no lo sostienen es peor que no emitir ninguno: se lee igual.
    duracion_h = (
        float(datos["t_wall"].max() - datos["t_wall"].min()) / 3600.0
        if len(datos) > 1
        else 0.0
    )
    tasa = n_obs / max(1e-9, duracion_h * 3600.0)
    no_decidible = []
    if duracion_h < 24.0:
        no_decidible.append(
            f"la corrida dura {duracion_h:.2f} h y la Sec. E.2 exige >= 24 h "
            f"continuas, cubriendo las sesiones asiatica, europea y americana"
        )
    if tasa < 5.0:
        no_decidible.append(
            f"solo {tasa:.2f} mediciones/s contra ~90 Hz de ciclo de control. "
            f"A_arm propaga el oscilador durante ~{1.0/max(tasa,1e-9):.1f} s entre "
            f"correcciones, asi que un error del 10 % en omega se acumula mucho "
            f"mas que en velocidad constante: el A/B mide el FEED, no el modelo"
        )
    if max(rho_vc, rho_arm) < rho_min:
        no_decidible.append(
            f"ambos rho_1 quedan por debajo del minimo detectable con n={n_obs} "
            f"({rho_min:.4f}): el rechazo no es distinguible del ruido muestral"
        )
    if no_decidible:
        a("  >> A/B NO DECIDIBLE TODAVIA. Lo que sigue es indicativo, no un veredicto:")
        for motivo in no_decidible:
            a(f"     - {motivo}")
        a("")

    if rho_arm < UMBRAL_RHO_MATERIAL and rho_arm < rho_vc:
        a(f"  >> ADOPTAR A_arm. |rho_1| baja de {rho_vc:.4f} a {rho_arm:.4f}, por")
        a(f"     debajo del umbral material {UMBRAL_RHO_MATERIAL:.2f}.")
        a("     CONSECUENCIA (Sec. E.4): la compuerta de Ljung-Box queda ABIERTA y")
        a("     la FASE 2 (ALS) se desbloquea. Es la via prevista para atacar")
        a("     r_S,base y r_EMD, medidos mal por 8x y 161x en la v1.2 y candidatos")
        a("     numero uno a explicar el NIS residual.")
    elif rho_vc >= UMBRAL_RHO_MATERIAL and rho_arm >= UMBRAL_RHO_MATERIAL:
        a(f"  >> NO ADOPTAR. Ambos modelos con |rho_1| >= {UMBRAL_RHO_MATERIAL:.2f} "
          f"({rho_vc:.4f} y {rho_arm:.4f}).")
        a("     El mismatch NO es el que creiamos: reabrir el diagnostico antes de")
        a("     seguir cambiando A.")
    else:
        a(f"  >> NO ADOPTAR: el armonico es PEOR ({rho_arm:.4f} contra {rho_vc:.4f}).")
        a("     >> VERIFICAR ANTES D.1 Y D.2. Un factor 2pi o un factor 125 se ven")
        a("       EXACTAMENTE asi, y ninguno de los dos da sintoma propio. Correr")
        a("       `python tests_v13.py` : test_D_trampa_2pi_periodo_implicito es la")
        a("       unica defensa contra ambos.")
    a("")
    a("  Nota: y1 NO entra en esta compuerta. Su rho_1 alto es artefacto del")
    a("  solapamiento de ventanas del EMD mas retencion de orden cero (Hilo Lento")
    a("  a 2 Hz contra Hilo Rapido a ~89 Hz), no evidencia sobre A.")
    return "\n".join(L)


def main(argv):
    directorio = DIR_TELEMETRIA_POR_DEFECTO
    descartar = -1  # -1 = deteccion automatica del transitorio
    resto = [x for x in argv[1:] if not x.startswith("--")]
    if resto:
        directorio = resto[0]
    episodio = None
    for x in argv[1:]:
        if x.startswith("--descartar="):
            descartar = int(x.split("=", 1)[1])
        elif x.startswith("--episodio="):
            episodio = int(x.split("=", 1)[1])

    try:
        datos = cargar_telemetria(directorio)
    except FileNotFoundError as err:
        print(f"[ERROR] {err}")
        return 1

    # Sec. C.5 de la v1.3: los ~30 episodios deben ser separables en el analisis.
    # Mezclarlos falsearia tanto Ljung-Box como el NIS, porque cada episodio
    # arranca con su propio transitorio de burn-in.
    if episodio is not None:
        if "id_episodio" not in datos.dtype.names:
            print("[ERROR] telemetria anterior a la v1.3: no lleva id_episodio")
            return 1
        mascara = datos["id_episodio"].astype(int) == episodio
        if not mascara.any():
            presentes = sorted(set(datos["id_episodio"].astype(int).tolist()))
            print(f"[ERROR] no hay muestras del episodio {episodio}. Presentes: {presentes}")
            return 1
        datos = datos[mascara]
        print(f"[FILTRO] Episodio {episodio}: {len(datos)} muestras.\n")

    print(generar_reporte(datos, descartar=descartar))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
