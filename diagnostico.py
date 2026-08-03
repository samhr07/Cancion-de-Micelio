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

    datos = np.concatenate(bloques)
    # Se descartan las filas nunca escritas (t_wall == 0) del ultimo bloque parcial.
    datos = datos[datos["t_wall"] > 0]
    return datos[np.argsort(datos["t_wall"])]


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
    for j in range(1, h + 1):
        rho_j = float(np.dot(x[j:], x[:-j]) / denom)
        Q += rho_j * rho_j / (n - j)
    Q *= n * (n + 2)

    p = float(chi2.sf(Q, h))
    critico = float(chi2.ppf(0.95, h))
    return {
        "n": int(n), "h": h, "Q": float(Q), "critico_95": critico, "p": p,
        "veredicto": "BLANCA" if p >= 0.05 else "RECHAZA BLANCURA",
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
    bruto = len(datos)
    dur = float(datos["t_wall"].max() - datos["t_wall"].min())

    if descartar < 0:
        descartar = detectar_fin_transitorio(datos)
    if descartar > 0 and descartar < len(datos) - 100:
        datos = datos[descartar:]
    elif descartar >= len(datos) - 100:
        descartar = 0  # No queda registro util: se reporta sobre todo

    a("=" * 78)
    a("REPORTE DE CONSISTENCIA DEL FILTRO - Fase 1 del orden de trabajo")
    a("=" * 78)
    a(f"Muestras: {bruto}   Duracion: {dur:.1f} s   "
      f"Cadencia media: {dur/max(1,bruto)*1e3:.2f} ms/ciclo")
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
    Y = np.column_stack([datos["y0"], datos["y1"]])
    resultados = {}
    for etiqueta, serie in (("y0 (precio)", datos["y0"]), ("y1 (residuo)", datos["y1"])):
        r = ljung_box(serie, h)
        resultados[etiqueta] = r
        if np.isfinite(r["Q"]):
            a(f"  {etiqueta:<14} Q={r['Q']:12.3f}  critico_95(h={h})="
              f"{r['critico_95']:8.3f}  p={r['p']:.3e}  -> {r['veredicto']}")
        else:
            a(f"  {etiqueta:<14} {r['veredicto']}")

    rm = ljung_box_multivariante(Y, h)
    resultados["multivariante"] = rm
    if np.isfinite(rm["Q"]):
        a(f"  {'multivariante':<14} Q={rm['Q']:12.3f}  critico_95(gl={rm['gl']})="
          f"{rm['critico_95']:8.3f}  p={rm['p']:.3e}  -> {rm['veredicto']}")
    else:
        a(f"  multivariante  {rm['veredicto']}")

    rechaza = any(
        r.get("veredicto") == "RECHAZA BLANCURA" for r in resultados.values()
    )
    a("")
    if rechaza:
        a("  >> FASE 2 BLOQUEADA. Hay predictibilidad remanente que la matriz")
        a("     cinematica A no capturo. Eso es MODEL MISMATCH, no mala calibracion")
        a("     de ruido. Recalibrar Q y R aqui seria tapar el sintoma: ALS")
        a("     absorberia el error de modelo dentro de Q y devolveria una respuesta")
        a("     confiadamente equivocada.")
        a("     Respuesta correcta: AUMENTO DE ESTADO (Sec. 3.6 del PDF) -- derivadas")
        a("     de orden superior o sesgos de friccion como nuevas incognitas.")
    else:
        a("  >> Compuerta superada: la innovacion pasa por blanca. La Fase 2 (ALS)")
        a("     esta habilitada.")

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

    a("")
    a("=" * 78)
    return "\n".join(L)


def main(argv):
    directorio = DIR_TELEMETRIA_POR_DEFECTO
    descartar = -1  # -1 = deteccion automatica del transitorio
    resto = [x for x in argv[1:] if not x.startswith("--")]
    if resto:
        directorio = resto[0]
    for x in argv[1:]:
        if x.startswith("--descartar="):
            descartar = int(x.split("=", 1)[1])

    try:
        datos = cargar_telemetria(directorio)
    except FileNotFoundError as err:
        print(f"[ERROR] {err}")
        return 1
    print(generar_reporte(datos, descartar=descartar))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
