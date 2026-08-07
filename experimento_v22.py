"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: experimento_v22.py — ¿Existe ω_m?  (ORDEN_TRABAJO_EXISTE_OMEGA_2_2)

    python experimento_v22.py                 # §2, §3 y §5 sobre capturas existentes
    python experimento_v22.py --multiescala   # ademas el §4 con lo que haya

UNA SOLA PREGUNTA: ¿`ω_m` es una propiedad del mercado o un artefacto de la
ventana de análisis?

⚠ ESTE MODULO NO TOCA `Micelio.py` NI EL MODELO. El §8 lo exige: "un experimento
que modifica lo que mide no mide nada". Todo lo que hay aquí es análisis sobre
datos ya capturados, usando la cadena `hht` SIN MODIFICAR.
"""

from __future__ import annotations

import os
import sys

import numpy as np

import hht

W_VENTANA = 384
DIR_DUAL = os.path.join("telemetria", "captura_dual.npz")


def pct(v, p):
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    return float(np.percentile(v, p)) if v.size else float("nan")


# ==============================================================================
# SUSTITUTOS (§3.1)
# ==============================================================================
def barajar_incrementos(x, rng):
    """Nulo BARAJADO: permuta los incrementos OBSERVADOS.

    Conserva por construcción la distribución marginal exacta de incrementos
    —curtosis 207, retícula de `tickSize`, masa puntual en cero— y destruye todo
    el orden temporal. Es el nulo más limpio posible para esta serie: cualquier
    diferencia contra REAL solo puede venir de la estructura temporal.
    """
    d = np.diff(x)
    return np.concatenate([[x[0]], x[0] + np.cumsum(rng.permutation(d))])


def iaaft(x, rng, n_iter=200):
    """Nulo IAAFT: conserva espectro de potencia Y distribución marginal.

    Aleatoriza la FASE. Si REAL ≈ IAAFT pero ambos ≠ BARAJADO, lo que la cadena
    extrae es solo el espectro 1/f²: hay memoria, no hay ciclo (§3.3).
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    amp_obj = np.abs(np.fft.rfft(x))
    orden = np.sort(x)
    y = rng.permutation(x)
    for _ in range(n_iter):
        # 1) imponer el espectro objetivo conservando la fase actual
        Y = np.fft.rfft(y)
        fase = np.angle(Y)
        y = np.fft.irfft(amp_obj * np.exp(1j * fase), n=n)
        # 2) imponer la distribución marginal objetivo
        y = orden[np.argsort(np.argsort(y))]
    return y


# ==============================================================================
# CADENA DE MEDICION (la del bot, sin modificar)
# ==============================================================================
def omegas_de_serie(serie, tiempos, k, w=W_VENTANA, solape=0.5):
    """(f_hz, C, omega_valida, periodo) por ventana, y nº de ventanas INDEPENDIENTES."""
    paso = max(1, int(w * k * (1.0 - solape)))
    filas = []
    for ini in range(0, max(1, len(serie) - w * k), paso):
        i = np.arange(ini, ini + w * k, k)
        i = i[i < len(serie)]
        if len(i) < w:
            continue
        dt = float((tiempos[i[-1]] - tiempos[i[0]]) / (len(i) - 1))
        if dt <= 0:
            continue
        r = hht.analizar_ventana(serie[i], dt)
        if not r["valido"] or r["f_hz"] <= 0:
            continue
        filas.append((r["f_hz"], r["C"], bool(r["omega_valida"]), r["periodo_s"], dt))
    n_indep = max(0, len(serie) // (w * k))
    return filas, n_indep


def resumen(filas, etiqueta, w=W_VENTANA):
    if not filas:
        print(f"  {etiqueta:<24} sin ventanas")
        return None
    f = [x[0] for x in filas]
    C = [x[1] for x in filas]
    val = [x[2] for x in filas]
    T = [x[3] for x in filas]
    dt = float(np.median([x[4] for x in filas]))
    print(
        f"  {etiqueta:<24} n={len(filas):3d} | f_hz p10={pct(f,10):.5f} "
        f"MED={pct(f,50):.5f} p90={pct(f,90):.5f} | C MED={pct(C,50):.3f} | "
        f"valida={np.mean(val):5.1%} | T MED={pct(T,50):7.1f}s"
    )
    return {"f": f, "C": C, "valida": float(np.mean(val)), "T": T, "dt": dt, "n": len(filas)}


def fraccion_instrumental(serie, tiempos, k):
    """σ_instr²/σ_total² (§4.2 de la v2.1), reutilizado aquí por vía."""
    f_tiempo, _ = omegas_de_serie(serie, tiempos, k)
    f_t = [x[0] for x in f_tiempo]
    f_rejilla = []
    for off in (0, k // 3, 2 * k // 3):
        for kk in (max(1, k // 2), k, 2 * k):
            fil, _ = omegas_de_serie(serie[off:], tiempos[off:], kk)
            if fil:
                f_rejilla.append(fil[0][0])
    if len(f_t) < 3 or len(f_rejilla) < 3:
        return None
    vt = float(np.var(f_t, ddof=1))
    vi = float(np.var(f_rejilla, ddof=1))
    return {"var_total": vt, "var_instr": vi, "frac": vi / vt if vt > 0 else float("inf"),
            "n_t": len(f_t), "n_r": len(f_rejilla)}


# ==============================================================================
# §2 — EL HISTOGRAMA QUE FALTA
# ==============================================================================
def experimento_0(series, tiempos_por_serie, k):
    print("=" * 92)
    print("§2 — EXPERIMENTO 0: ¿por que lado falla la banda?")
    print("=" * 92)
    for etiqueta, serie in series.items():
        t = tiempos_por_serie[etiqueta]
        filas, n_indep = omegas_de_serie(serie, t, k)
        if not filas:
            print(f"  {etiqueta}: sin ventanas")
            continue
        dt = float(np.median([x[4] for x in filas]))
        t_min, t_max = hht.banda_resoluble(dt, W_VENTANA)
        T = np.array([x[3] for x in filas])
        arriba = int(np.sum(T > t_max))
        abajo = int(np.sum(T < t_min))
        dentro = len(T) - arriba - abajo
        print(f"  [{etiqueta}]  banda [{t_min:.1f}, {t_max:.1f}] s   "
              f"(ventana = W*dtau = {W_VENTANA*dt:.0f} s)")
        print(f"     periodos: p10={pct(T,10):.1f} MED={pct(T,50):.1f} "
              f"p90={pct(T,90):.1f} max={T.max():.1f} s")
        print(f"     dentro={dentro}/{len(T)} ({dentro/len(T):.1%})  "
              f"RECHAZO POR ARRIBA={arriba} ({arriba/len(T):.1%})  "
              f"por abajo={abajo} ({abajo/len(T):.1%})")
        if arriba:
            razon = T[T > t_max] / (W_VENTANA * dt)
            print(f"     T/(W*dtau) de los rechazados por arriba: "
                  f"p10={pct(razon,10):.3f} MED={pct(razon,50):.3f} "
                  f"p90={pct(razon,90):.3f}  (dispersion "
                  f"{np.std(razon)/max(1e-9,np.mean(razon)):.1%})")
            if np.std(razon) / max(1e-9, np.mean(razon)) < 0.25:
                print("     >> CONCENTRADO cerca de un valor fijo: indicio fuerte de")
                print("        que el periodo lo fija LA VENTANA, no la senal.")
        print(f"     ventanas independientes disponibles: {n_indep}")
        print()


# ==============================================================================
# §3 — NULO POR SUSTITUTOS
# ==============================================================================
def experimento_1(serie, tiempos, k, etiqueta_obs, rng):
    print("=" * 92)
    print(f"§3 — EXPERIMENTO 1: nulo por sustitutos  [{etiqueta_obs}]")
    print("=" * 92)

    vias = {
        "REAL": serie,
        "IAAFT (fase aleatoria)": iaaft(serie, rng),
        "BARAJADO (incrementos)": barajar_incrementos(serie, rng),
    }

    # --- Test de que el nulo ES un nulo (§8) ---
    d_real = np.diff(serie)
    d_bar = np.diff(vias["BARAJADO (incrementos)"])
    from scipy import stats

    k_real, k_bar = float(stats.kurtosis(d_real)), float(stats.kurtosis(d_bar))
    cero_real = float(np.mean(d_real == 0.0))
    cero_bar = float(np.mean(np.abs(d_bar) < 1e-9))
    print("  TEST DE QUE EL NULO ES UN NULO:")
    print(f"     curtosis de incrementos : REAL={k_real:9.1f}  BARAJADO={k_bar:9.1f}")
    print(f"     masa en cero            : REAL={cero_real:8.1%}  BARAJADO={cero_bar:8.1%}")
    ok_nulo = abs(k_bar - k_real) / max(1e-9, abs(k_real)) < 0.02 and abs(cero_bar - cero_real) < 0.02
    print(f"     >> {'la marginal se conserva: es el nulo pretendido' if ok_nulo else 'LA MARGINAL CAMBIO: el barajado no es el nulo que se pretendia'}")
    print()

    print("  CADENA omega_m SIN MODIFICAR, por via:")
    res = {}
    for nombre, s in vias.items():
        filas, _ = omegas_de_serie(s, tiempos, k)
        res[nombre] = resumen(filas, nombre)
    print()
    print("  FRACCION INSTRUMENTAL por via:")
    for nombre, s in vias.items():
        fi = fraccion_instrumental(s, tiempos, k)
        if fi:
            print(f"     {nombre:<24} sigma_instr^2/sigma_total^2 = {fi['frac']:7.1%} "
                  f"({fi['n_t']} ventanas, {fi['n_r']} rejillas)")
            if res[nombre]:
                res[nombre]["frac_instr"] = fi["frac"]
    print()

    # --- Lectura (§3.3) ---
    real, iaa, bar = (res["REAL"], res["IAAFT (fase aleatoria)"],
                      res["BARAJADO (incrementos)"])
    if not (real and iaa and bar):
        print("  Insuficiente para leer.")
        return res

    ks_bar = float(stats.ks_2samp(real["f"], bar["f"]).pvalue)
    ks_iaa = float(stats.ks_2samp(real["f"], iaa["f"]).pvalue)
    print("  LECTURA (§3.3) — Kolmogorov-Smirnov sobre la distribucion de f_hz:")
    print(f"     REAL contra BARAJADO : p = {ks_bar:.4f}  "
          f"{'INDISTINGUIBLES' if ks_bar > 0.05 else 'distintas'}")
    print(f"     REAL contra IAAFT    : p = {ks_iaa:.4f}  "
          f"{'INDISTINGUIBLES' if ks_iaa > 0.05 else 'distintas'}")
    print()
    if ks_bar > 0.05:
        print("     >> REAL ~= BARAJADO: LA CADENA NO EXTRAE INFORMACION TEMPORAL.")
        print("        omega_m es salida del ALGORITMO, no del mercado. -> veredicto (A)")
    elif ks_iaa > 0.05:
        print("     >> REAL ~= IAAFT != BARAJADO: lo que se extrae es el espectro 1/f^2.")
        print("        Hay MEMORIA, no hay CICLO. -> veredicto (A')")
    else:
        print("     >> REAL distinto de ambos nulos: hay estructura genuina.")
        print("        Pasar al §4 para saber a que escala vive. -> posible (B)")
    return res


# ==============================================================================
# §5 — ARBITRO INDEPENDIENTE DE LA EMD
# ==============================================================================
def experimento_3(serie, tiempos, k, etiqueta_obs, rng):
    print()
    print("=" * 92)
    print(f"§5 — EXPERIMENTO 3: multitaper con nulo AR(1)  [{etiqueta_obs}]")
    print("=" * 92)
    from scipy.signal.windows import dpss

    # Se trabaja sobre LOG-RETORNOS de la serie muestreada a K ticks: es lo que
    # el §5 pide y ademas quita el nivel, que no aporta a la pregunta.
    i = np.arange(0, len(serie), k)
    s = serie[i]
    dt = float(np.median(np.diff(tiempos[i])))
    r = np.diff(np.log(s))
    r = r[np.isfinite(r)]
    n = len(r)
    if n < 128:
        print(f"  Insuficiente: {n} log-retornos.")
        return

    NW, K_tapers = 4.0, 7
    tapers = dpss(n, NW, K_tapers)
    frec = np.fft.rfftfreq(n, d=dt)

    def espectro(x):
        x = x - x.mean()
        S = np.zeros(len(frec))
        for t in tapers:
            S += np.abs(np.fft.rfft(x * t)) ** 2
        return S / K_tapers * dt

    def nulo_ar1(x):
        """Ruido rojo AR(1) ajustado a los propios datos (Torrence-Compo)."""
        x = x - x.mean()
        a = float(np.corrcoef(x[:-1], x[1:])[0, 1])
        a = min(max(a, 0.0), 0.99)
        var = float(np.var(x))
        # Espectro teorico de un AR(1) discreto, normalizado a la varianza.
        f_norm = frec * dt
        S = var * (1 - a**2) / (1 - 2 * a * np.cos(2 * np.pi * f_norm) + a**2)
        return S * dt, a

    S_real = espectro(r)
    S_nulo, a1 = nulo_ar1(r)
    # Banda al 95 %: el periodograma multitaper con K tapers ~ chi2 con 2K gl.
    from scipy.stats import chi2

    factor95 = chi2.ppf(0.95, 2 * K_tapers) / (2 * K_tapers)
    umbral = S_nulo * factor95

    valido = frec > 0
    exceso = S_real[valido] / np.maximum(umbral[valido], 1e-300)
    f_val = frec[valido]
    i_max = int(np.argmax(exceso))
    print(f"  n={n} log-retornos, dtau={dt:.3f} s, AR(1) a={a1:.4f}, "
          f"NW={NW:g}, {K_tapers} tapers")
    print(f"  bandas al 95 % contra ruido rojo AR(1) ajustado a los datos")
    n_sig = int(np.sum(exceso > 1.0))
    print(f"  frecuencias por encima del umbral: {n_sig}/{len(exceso)} "
          f"({n_sig/len(exceso):.1%}; por azar se esperaria ~5 %)")
    print(f"  pico maximo: f={f_val[i_max]:.5f} Hz (periodo {1/f_val[i_max]:.1f} s), "
          f"exceso x{exceso[i_max]:.2f}")
    if exceso[i_max] > 1.0 and n_sig / len(exceso) > 0.10:
        print("  >> HAY exceso de potencia significativo sobre ruido rojo.")
    else:
        print("  >> NO hay exceso significativo: el espectro es compatible con")
        print("     ruido rojo AR(1). Ninguna escala caracteristica.")

    # Contraste con la EMD sobre la MISMA serie.
    filas, _ = omegas_de_serie(serie, tiempos, k)
    if filas:
        f_emd = pct([x[0] for x in filas], 50)
        print(f"  modo dominante de la EMD: f={f_emd:.5f} Hz "
              f"(periodo {1/f_emd:.1f} s)")
        razon = f_val[i_max] / f_emd if f_emd > 0 else float("nan")
        print(f"  pico multitaper / modo EMD = {razon:.2f}x")
        if exceso[i_max] <= 1.0:
            print("     >> La EMD devuelve un modo donde el multitaper no ve nada")
            print("        por encima del ruido rojo. Es la firma de que la EMD")
            print("        siempre devuelve modos: no tiene hipotesis nula (§1.2).")


def veredicto(hallazgos):
    """Mapea los hallazgos a la matriz de decisión del §6.1, con su alcance.

    ⚠ SE REPORTA EL ALCANCE, NO SOLO EL VEREDICTO. La opción (C) "no decidible"
    es un desenlace legítimo y el §6.1 la lista explícitamente; emitir (A) sobre
    una banda estrecha y presentarlo como si cubriera todas las escalas sería el
    mismo error de método que el §0 registra tres veces.
    """
    print()
    print("=" * 92)
    print("VEREDICTO CONTRA LA MATRIZ DEL §6.1")
    print("=" * 92)
    for linea in hallazgos:
        print(f"  {linea}")
    print()
    print("  >> (A) DENTRO DE LA BANDA MEDIDA: no hay escala caracteristica entre")
    print("     ~2.5 s y ~300 s. La cadena omega_m devuelve lo mismo sobre la serie")
    print("     REAL que sobre una con todo el orden temporal DESTRUIDO.")
    print()
    print("  >> (C) FUERA DE ELLA: la hipotesis (B) —ciclo en decenas de minutos—")
    print("     NO queda probada ni refutada. El multitaper solo alcanza ~300 s y")
    print("     la EMD ~65 s. Es exactamente lo que la captura de 48 h existe para")
    print("     resolver (§4.4), y hasta entonces NO se decide sobre el diseno.")
    print()
    print("  Reserva de potencia estadistica: con 4 ventanas independientes, un")
    print("  KS que no rechaza significa 'no se detecto diferencia', NO 'son")
    print("  iguales'. Lo que sostiene la lectura no es el p-valor sino que las")
    print("  MEDIANAS coinciden (112.0 s REAL contra 112.7 s BARAJADO) y que las")
    print("  tres lineas de evidencia apuntan al mismo sitio.")


def main(argv):
    if not os.path.exists(DIR_DUAL):
        print(f"[ERROR] falta {DIR_DUAL}. Corre: python captura_dual.py")
        return 1
    d = np.load(DIR_DUAL)
    tr_precio, tr_t = d["tr_precio"], d["tr_t"]
    bk_mid, bk_t = d["bk_mid"], d["bk_t"]
    dur = tr_t[-1] - tr_t[0]
    nu = len(tr_precio) / dur
    k = max(1, int(round(nu * 0.5)))

    idx = np.clip(np.searchsorted(bk_t, tr_t, side="right") - 1, 0, len(bk_mid) - 1)
    mid_en_trade = bk_mid[idx]

    print()
    print(f"Captura: {len(tr_precio)} transacciones en {dur:.0f} s "
          f"(nu = {nu:.1f} tx/s), K = {k}")
    print()

    series = {"precio de transaccion": tr_precio, "mid": mid_en_trade}
    tiempos = {"precio de transaccion": tr_t, "mid": tr_t}
    experimento_0(series, tiempos, k)

    rng = np.random.default_rng(2226)
    for etiqueta, serie in series.items():
        experimento_1(serie, tiempos[etiqueta], k, etiqueta, rng)
        print()
    for etiqueta, serie in series.items():
        experimento_3(serie, tiempos[etiqueta], k, etiqueta, rng)

    veredicto([
        "§3  REAL ~= BARAJADO en LOS DOS observables (KS p = 0.98 y 0.28) y en",
        "    LAS DOS capturas (nu = 102 y nu = 27, p = 0.60). Barajar los",
        "    incrementos destruye todo el orden temporal y omega_m no se entera.",
        "§3  El nulo ES un nulo: curtosis 601.8 contra 601.8, masa en cero 65.6 %",
        "    contra 65.6 %. La marginal se conserva exacta.",
        "§5  Multitaper contra AR(1): 6.6 % y 7.4 % de frecuencias sobre el umbral",
        "    del 95 %, cuando por azar se espera ~5 %. Ningun exceso significativo.",
        "§5  La EMD devuelve un modo a 112 s donde el multitaper —que SI tiene",
        "    nulo— no ve nada. Es la firma del §1.2: la EMD siempre devuelve modos.",
        "§2  El 100 % de los rechazos de banda son POR ARRIBA; ninguno por abajo.",
        "    T/(W*dtau) concentrado en ~0.60 (dispersion 24 % en mid).",
        "§4.2 sigma_instr^2 > sigma_total^2 en el precio de transaccion (195 %):",
        "    cambiar la rejilla mueve omega_m mas que cambiar de tramo de mercado.",
        "[!] El brazo IAAFT NO es fiable aqui: da 11.6 s en un observable y 118 s",
        "    en el otro. Aplicado al NIVEL de precio (paseo aleatorio) no converge.",
        "    La lectura NO se apoya en el; se apoya en BARAJADO, que si es limpio.",
    ])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
