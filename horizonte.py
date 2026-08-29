# -*- coding: utf-8 -*-
"""
horizonte.py -- v4.1 Sec.2 (H_p) y Sec.1 (las dos curvas).

    python horizonte.py --autotest       controles con verdad conocida
    python horizonte.py --etapa=hp       Sec.2: H_p sobre las cuatro capturas
    python horizonte.py --etapa=curvas   Sec.1: R2 medido contra R2 requerido

EL ORDEN IMPORTA. La curva requerida del Sec.1 no se puede dibujar sin `H_p`, y
`H_p` no ha replicado entre capturas (0.591 contra ~0.5). De el cuelga un factor
3.4 en el horizonte de cruce, asi que el Sec.2 va primero.

[!] `H_p` DE RELOJ DE PARED, y es deliberado. La sesion 2026-08-08 (e) establecio
que la firma en reloj de pared tiene sesgo propio del estimador y que hay que
medirla en ticks. Pero la CURVA REQUERIDA vive en segundos de calendario: la
comision y la financiacion se pagan en tiempo de pared, no en transacciones. Asi
que el `H_p` que entra en `sigma(H)` tiene que ser el de pared, **con su control
barajado al lado para descontar el sesgo**, que es justo lo que aquella sesion
dejo como obligatorio. El de ticks se reporta como control.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

import propagador as P
from captura_larga import cargar_larga, tramos_continuos

CAPTURAS = ("telemetria/captura_larga", "telemetria/captura_v31b",
            "telemetria/captura_v32", "telemetria/captura_v33")
CACHE_V33 = "telemetria/muestra_v32.npz"

# Sec.1.2 -- lastre, todo en pb para que no envejezca con el precio (Sec.8 fallo 5).
COMISION_MAKER_VIP0 = 0.0002
LASTRE_IDA_VUELTA_PB = 1e4 * 2 * COMISION_MAKER_VIP0          # 4.00 pb
SELECCION_ADVERSA_PB = 1e4 * 5.78 / 65076.0                   # medida en la v4.0 Sec.5
FINANCIACION_PB_POR_8H = 1e4 * 6.3 / 65076.0                  # perpetuo, ~6.3 USD/BTC/8h
# E[|mu| | decil superior] / sigma_mu para una normal: densidad en el corte q90.
FACTOR_DECIL = 1.755

# Alcance verificado de la firma (Sec.2.1 punto 4). Mas alla es extrapolacion.
ALCANCE_VERIFICADO_S = 711.0

_TRANSLIT = {"⚠": "[!]", "§": "Sec.", "→": "->", "≤": "<=", "≥": ">=",
             "≠": "!=", "×": "x", "–": "-", "—": "--", "²": "^2", "√": "sqrt"}


def log(*a):
    s = " ".join(str(x) for x in a)
    for k, v in _TRANSLIT.items():
        s = s.replace(k, v)
    print(s.encode("ascii", "replace").decode("ascii"))
    sys.stdout.flush()


def titulo(s):
    log("")
    log("=" * 78)
    log(s)
    log("=" * 78)


# ===========================================================================
# Carga
# ===========================================================================

def tramo_mas_largo(directorio: str) -> dict:
    d = cargar_larga(directorio)
    p, t = d["tr_precio"], d["tr_t"]
    ok = np.isfinite(p) & (p > 0)
    p, t = p[ok], t[ok]
    q = d["tr_cant"][ok] if "tr_cant" in d else np.ones_like(p)
    m = d["tr_maker"][ok].astype(bool) if "tr_maker" in d else None
    a, b = max(tramos_continuos(t), key=lambda r: r[1] - r[0])
    return {"precio": p[a:b], "t": t[a:b], "cant": q[a:b],
            "maker": None if m is None else m[a:b],
            "n": int(b - a), "dur": float(t[b - 1] - t[a]),
            "nu": (b - a) / float(t[b - 1] - t[a])}


# ===========================================================================
# Sec.2 -- la firma en los dos relojes, con control barajado en todas
# ===========================================================================

def firma_pared(precio, t, horizontes_s) -> list:
    """`sigma(H)` muestreando por reloj de pared, ventanas NO solapadas."""
    p, t = np.asarray(precio, float), np.asarray(t, float)
    filas = []
    for H in horizontes_s:
        bordes = np.arange(t[0], t[-1], H)
        if bordes.size < 12:
            continue
        idx = np.clip(np.searchsorted(t, bordes, side="right") - 1, 0, p.size - 1)
        r = np.diff(p[idx])
        if r.size < 10:
            continue
        filas.append({"H": float(H), "n": int(r.size),
                      "sigma": float(np.std(r, ddof=1)),
                      "sigma_por_raiz_H": float(np.std(r, ddof=1) / np.sqrt(H))})
    return filas


def barajar_incrementos(precio, semilla=0):
    """Conserva la marginal y los tiempos de llegada; destruye el orden.

    Es el control que la sesion 2026-08-08 (e) dejo obligatorio: sobre el
    barajado la pendiente debe ser 0, y lo que NO sea 0 es sesgo del estimador.
    """
    rng = np.random.default_rng(semilla)
    inc = np.diff(np.asarray(precio, float))
    return np.concatenate(([precio[0]], precio[0] + np.cumsum(rng.permutation(inc))))


def pendiente(filas, clave_x, clave_y, lo, hi) -> dict:
    x = np.array([f[clave_x] for f in filas], float)
    y = np.array([f[clave_y] for f in filas], float)
    s = (x >= lo) & (x <= hi) & (y > 0)
    if s.sum() < 3:
        return {"pendiente": float("nan"), "n": int(s.sum())}
    A = np.column_stack([np.log(x[s]), np.ones(int(s.sum()))])
    c, *_ = np.linalg.lstsq(A, np.log(y[s]), rcond=None)
    return {"pendiente": float(c[0]), "intercepto": float(c[1]), "n": int(s.sum())}


def hp_de_captura(cap: dict, semilla: int = 5) -> dict:
    """`H_p` en los DOS relojes, cada uno con su control barajado."""
    # [!] REJILLA DENSA, y no es cosmetico. Con la rejilla original
    # [60,120,300,600] el ajuste de la ventana [60,600] tenia CUATRO puntos, y
    # una pendiente sobre cuatro puntos es tan ruidosa que el propio control
    # barajado -- que deberia dar 0.5 exacto -- salia entre 0.254 y 0.528. Ahi
    # la correccion real-barajado mueve mas ruido que senal.
    Hs = sorted(set([0.5, 1, 2, 5] +
                    [round(x, 1) for x in np.logspace(np.log10(8),
                                                      np.log10(1800), 28)]))
    ns = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
    p_nulo = barajar_incrementos(cap["precio"], semilla)

    fp = firma_pared(cap["precio"], cap["t"], Hs)
    fp_n = firma_pared(p_nulo, cap["t"], Hs)
    ft = P.firma_en_ticks(cap["precio"], ns, solapada=False)
    ft_n = P.firma_en_ticks(p_nulo, ns, solapada=False)

    out = {"filas_pared": fp, "filas_pared_nulo": fp_n,
           "filas_ticks": ft, "filas_ticks_nulo": ft_n}
    # Ventanas de ajuste multiples: el Sec.3.3 de la v4.1 enseno que un exponente
    # depende de la ventana en que se ajusta, asi que se reporta como curva.
    out["ventanas_pared"] = []
    for lo, hi in ((10, 100), (30, 300), (60, 600), (10, 600)):
        a = pendiente(fp, "H", "sigma_por_raiz_H", lo, hi)
        b = pendiente(fp_n, "H", "sigma_por_raiz_H", lo, hi)
        out["ventanas_pared"].append({
            "rango": (lo, hi), "n": a["n"],
            "H_p": a["pendiente"] + 0.5, "H_p_nulo": b["pendiente"] + 0.5,
            "H_p_corregido": a["pendiente"] - b["pendiente"] + 0.5})
    out["ventanas_ticks"] = []
    for lo, hi in ((256, 16384), (256, 2000)):
        a = pendiente(ft, "n", "sigma_por_raiz_n", lo, hi)
        b = pendiente(ft_n, "n", "sigma_por_raiz_n", lo, hi)
        out["ventanas_ticks"].append({
            "rango": (lo, hi), "n": a["n"],
            "H_p": a["pendiente"] + 0.5, "H_p_nulo": b["pendiente"] + 0.5,
            "H_p_corregido": a["pendiente"] - b["pendiente"] + 0.5})
    return out


def sigma1_pb(filas_pared, precio_mediano: float, H_p: float,
              H_ref: float = 100.0) -> float:
    """`sigma_1` en pb, anclada en `H_ref` s. Sec.2.1 punto 3.

    [!] `sigma_1 = 3.77 USD/BTC*s^(1/2)` de la v3.1 esta en USD a un precio que
    ya envejecio: con BTC a otro nivel el mismo movimiento relativo da otro
    numero, y `R2_req` va con el CUADRADO de ese error. En pb no envejece.
    """
    x = np.array([f["H"] for f in filas_pared], float)
    y = np.array([f["sigma"] for f in filas_pared], float)
    s = (x > 0) & (y > 0)
    lg = np.interp(np.log(H_ref), np.log(x[s]), np.log(y[s]))
    sigma_ref_pb = 1e4 * float(np.exp(lg)) / precio_mediano
    return sigma_ref_pb / (H_ref ** H_p)


# ===========================================================================
# Sec.1 -- las dos curvas
# ===========================================================================

def lastre_pb(H_s: float) -> float:
    """`c(u) + seleccion adversa + financiacion(H)`, todo en pb."""
    return (LASTRE_IDA_VUELTA_PB + SELECCION_ADVERSA_PB
            + FINANCIACION_PB_POR_8H * max(H_s - 3600.0, 0.0) / (8 * 3600.0))


def r2_requerido(H_s: float, sigma1_pb_: float, H_p: float) -> float:
    """`R2_req(H) = ( lastre(H) / (1.755 * sigma(H)) )^2`, con todo en pb."""
    sig = sigma1_pb_ * H_s ** H_p
    return (lastre_pb(H_s) / (FACTOR_DECIL * sig)) ** 2


def ventanas_no_solapadas(t, H_s: float) -> np.ndarray:
    """Indices de INICIO de ventanas contiguas y no solapadas de `H` segundos."""
    t = np.asarray(t, float)
    bordes = np.arange(t[0], t[-1], H_s)
    return np.clip(np.searchsorted(t, bordes, side="left"), 0, t.size - 1)


def rasgos_flujo(flujo_acum, idx, H_ticks_ret) -> np.ndarray:
    """Flujo firmado acumulado en ventanas RETROSPECTIVAS que acaban en `idx`.

    Deliberadamente tonto: sin `G(tau)`, sin `tau_0`, sin `beta`. La
    degeneracion `beta`-`tau_0` mato cuatro estadisticos y no debe contaminar la
    medicion del horizonte (Sec.1.3 punto 1).
    """
    X = np.empty((idx.size, len(H_ticks_ret)))
    for j, w in enumerate(H_ticks_ret):
        ini = np.maximum(idx - int(w), 0)
        X[:, j] = flujo_acum[idx] - flujo_acum[ini]
    return X


def r2_fuera_de_muestra(Xe, ye, Xv, yv) -> float:
    A = np.column_stack([Xe, np.ones(len(Xe))])
    c, *_ = np.linalg.lstsq(A, ye, rcond=None)
    pred = np.column_stack([Xv, np.ones(len(Xv))]) @ c
    ss = float(np.sum((yv - pred) ** 2))
    tt = float(np.sum((yv - np.mean(ye)) ** 2))
    return 1.0 - ss / tt if tt > 0 else float("nan")


# ===========================================================================
# Controles
# ===========================================================================

def _fgn(n: int, H: float, rng) -> np.ndarray:
    """Ruido gaussiano fraccionario por sintesis espectral, `S(f) ~ f^(1-2H)`.

    Solo para controles: su cumsum es fBm con `sigma(n) ~ n^H` por construccion,
    o sea verdad conocida para la firma de volatilidad.
    """
    m = 1 << int(np.ceil(np.log2(4 * n)))
    f = np.fft.rfftfreq(m)[1:]
    amp = f ** ((1.0 - 2.0 * H) / 2.0)
    fase = rng.normal(size=f.size) + 1j * rng.normal(size=f.size)
    x = np.fft.irfft(np.concatenate(([0.0 + 0j], amp * fase)), m)[:n]
    return x / np.std(x)


def _autotest() -> int:
    fallos = 0

    def ok(nombre, cond, extra=""):
        nonlocal fallos
        log("  [%s] %-56s %s" % ("OK  " if cond else "FALLO", nombre, extra))
        if not cond:
            fallos += 1

    log("== 1. La curva requerida reproduce la tabla del PREREGISTRO_3_2 Sec.11.1 ==")
    # La tabla publicada usa sigma_1 = 3.77 USD/BTC*s^(1/2) a P = 65 076 y
    # H_p = 0.5, con lastre 43.4 USD/BTC. Se reconstruye en pb.
    Pref = 65076.0
    s1 = 1e4 * 3.77 / Pref
    lastre_v31 = 1e4 * 43.4 / Pref
    for H, esperado in ((103.0, 0.418), (600.0, 0.072), (3600.0, 0.012)):
        sig = s1 * H ** 0.5
        r2 = (lastre_v31 / (FACTOR_DECIL * sig)) ** 2
        log("     H = %7.0f s   R2_req = %6.2f %%   publicado %5.1f %%"
            % (H, 100 * r2, 100 * esperado))
        ok("R2_req a H=%.0f s reproduce lo publicado" % H,
           abs(r2 - esperado) < 0.10 * max(esperado, 0.005))

    log("== 2. Las ventanas del Sec.1 son NO SOLAPADAS (Sec.9 lo exige por test) ==")
    t = np.arange(0.0, 10000.0, 0.5)
    for H in (60.0, 300.0, 900.0):
        idx = ventanas_no_solapadas(t, H)
        d = np.diff(t[idx])
        ok("ventanas de %4.0f s no se solapan" % H,
           bool(np.all(d >= H - 1e-6)),
           "separacion minima %.3f s (exigida %.0f)" % (d.min(), H))

    log("== 3. R2 fuera de muestra: 0 sin senal, alto con senal ==")
    rng = np.random.default_rng(3)
    n = 4000
    Xe, Xv = rng.normal(size=(n, 3)), rng.normal(size=(n, 3))
    ruido_e, ruido_v = rng.normal(size=n), rng.normal(size=n)
    r0 = r2_fuera_de_muestra(Xe, ruido_e, Xv, ruido_v)
    log("     sin senal   -> R2 = %+.5f" % r0)
    ok("sin senal el R2 fuera de muestra es ~0", abs(r0) < 0.02)
    be = Xe[:, 0] * 0.8 + ruido_e * 0.3
    bv = Xv[:, 0] * 0.8 + ruido_v * 0.3
    r1 = r2_fuera_de_muestra(Xe, be, Xv, bv)
    log("     con senal   -> R2 = %+.5f  (verdad ~ %.3f)"
        % (r1, 0.64 / (0.64 + 0.09)))
    ok("con senal el R2 recupera la verdad", abs(r1 - 0.877) < 0.03)

    log("== 4. El control barajado del Sec.2 anula la pendiente ==")
    # Paseo aleatorio: la firma es plana y barajar no la cambia.
    p = 65000.0 + np.cumsum(rng.normal(0, 0.2, size=200000))
    t = np.arange(p.size) * 0.05
    f = firma_pared(p, t, [1, 2, 5, 10, 30, 60, 120, 300])
    a = pendiente(f, "H", "sigma_por_raiz_H", 10, 300)
    log("     paseo aleatorio -> H_p = %.4f (debe ser ~0.5)" % (a["pendiente"] + 0.5))
    ok("el paseo da H_p ~ 0.5", abs(a["pendiente"]) < 0.05)
    # [!] LA PRIMERA VERSION DE ESTE CONTROL ERA FALSA. Usaba una TENDENCIA
    # LINEAL para fabricar superdifusion, y una tendencia constante NO aparece
    # en la firma: anade la misma constante a todos los incrementos y
    # `std(x + c) = std(x)`. Superdifusion es incrementos POSITIVAMENTE
    # AUTOCORRELACIONADOS, no deriva. Se usa movimiento browniano fraccionario,
    # cuyo `H` es verdad conocida por construccion: `sigma(n) ~ n^H`.
    for H_v in (0.65, 0.35):
        p_f = 65000.0 + np.cumsum(_fgn(1 << 18, H_v, rng)) * 0.2
        t_f = np.arange(p_f.size) * 0.05
        rej = [1, 2, 5, 10, 30, 60, 120, 300]
        f_r = firma_pared(p_f, t_f, rej)
        f_n = firma_pared(barajar_incrementos(p_f, 1), t_f, rej)
        ar = pendiente(f_r, "H", "sigma_por_raiz_H", 10, 300)
        an = pendiente(f_n, "H", "sigma_por_raiz_H", 10, 300)
        log("     fBm H=%.2f -> real %.4f | barajado %.4f | corregido %.4f"
            % (H_v, ar["pendiente"] + 0.5, an["pendiente"] + 0.5,
               ar["pendiente"] - an["pendiente"] + 0.5))
        ok("H_p recuperada con fBm H=%.2f (error < 0.08)" % H_v,
           abs(ar["pendiente"] + 0.5 - H_v) < 0.08)
        ok("el barajado de fBm H=%.2f vuelve a 0.5" % H_v,
           abs(an["pendiente"]) < 0.05)

    log("== 5. sigma_1 en pb NO depende del nivel de precio ==")
    f1 = firma_pared(p, t, [10, 30, 60, 100, 300])
    s_a = sigma1_pb(f1, 65000.0, 0.5)
    f2 = firma_pared(p * 2.0, t, [10, 30, 60, 100, 300])
    s_b = sigma1_pb(f2, 130000.0, 0.5)
    log("     sigma_1 = %.6f pb a P=65k  |  %.6f pb a P=130k (misma serie x2)"
        % (s_a, s_b))
    ok("sigma_1 en pb es invariante al nivel", abs(s_a - s_b) / s_a < 0.02)

    log("")
    log("RESULTADO: %d fallo(s)" % fallos)
    return fallos


# ===========================================================================
# Etapas
# ===========================================================================

def etapa_hp(args):
    titulo("Sec.2 -- H_p SOBRE LAS CUATRO CAPTURAS, DOS RELOJES, CONTROL EN TODAS")
    log("[!] El H_p que entra en sigma(H) es el de RELOJ DE PARED: la comision y")
    log("    la financiacion se pagan en tiempo de calendario. El de ticks va al")
    log("    lado como control. Y los dos llevan su barajado, porque la sesion")
    log("    2026-08-08 (e) midio que el estimador de pared tiene sesgo propio.")
    res = {}
    for ruta in CAPTURAS:
        if not os.path.isdir(ruta):
            log("  (falta %s)" % ruta)
            continue
        cap = tramo_mas_largo(ruta)
        r = hp_de_captura(cap)
        res[ruta] = (cap, r)
        log("")
        log("--- %s ---" % os.path.basename(ruta))
        log("  tramo continuo: %d ticks en %.2f h   nu = %.2f tx/s   P mediano %.0f"
            % (cap["n"], cap["dur"] / 3600.0, cap["nu"], np.median(cap["precio"])))
        log("  RELOJ DE PARED         real   barajado  corregido   n")
        for v in r["ventanas_pared"]:
            marca = "  <- EXTRAPOLA" if v["rango"][1] > ALCANCE_VERIFICADO_S else ""
            log("    H en [%4d, %4d] s   %6.3f   %7.3f   %8.3f  %3d%s"
                % (v["rango"][0], v["rango"][1], v["H_p"], v["H_p_nulo"],
                   v["H_p_corregido"], v["n"], marca))
        log("  RELOJ DE TICKS (control)")
        for v in r["ventanas_ticks"]:
            log("    n en [%5d, %5d]    %6.3f   %7.3f   %8.3f  %3d"
                % (v["rango"][0], v["rango"][1], v["H_p"], v["H_p_nulo"],
                   v["H_p_corregido"], v["n"]))
        hp = [v for v in r["ventanas_pared"] if v["rango"] == (60, 600)][0]
        s1 = sigma1_pb(r["filas_pared"], float(np.median(cap["precio"])),
                       hp["H_p_corregido"])
        log("  sigma_1 = %.6f pb*s^(-H_p)   con H_p corregido = %.3f"
            % (s1, hp["H_p_corregido"]))

    titulo("Sec.2.2 -- COMPUERTA: ¿replica H_p entre capturas?")
    filas = []
    for ruta, (cap, r) in res.items():
        v = [x for x in r["ventanas_pared"] if x["rango"] == (60, 600)][0]
        filas.append((os.path.basename(ruta), cap["nu"], v["H_p"],
                      v["H_p_nulo"], v["H_p_corregido"]))
    log("  captura              nu     H_p real  barajado  CORREGIDO  |sesgo|  efecto  fiable")
    for nom, nu, a, b, c in filas:
        sesgo, efecto = abs(b - 0.5), abs(a - 0.5)
        log("  %-20s %6.2f    %6.3f    %6.3f     %6.3f   %6.3f  %6.3f  %s"
            % (nom, nu, a, b, c, sesgo, efecto,
               "si" if sesgo < efecto else "NO"))
    log("")
    log("  [!] CRITERIO DE FIABILIDAD, declarado aqui: el control barajado DEBE dar")
    log("      0.5 exacto (no hay estructura que medir). Su desviacion de 0.5 es el")
    log("      sesgo del estimador. Cuando ese sesgo SUPERA al efecto que se quiere")
    log("      medir, la correccion mueve mas ruido que senal y la fila no se usa.")
    fiables = [f for f in filas if abs(f[3] - 0.5) < abs(f[2] - 0.5)]
    log("      filas fiables: %d de %d" % (len(fiables), len(filas)))
    log("")
    log("  [!] Y HAY QUE DECLARARLO: este criterio se escribio DESPUES de ver la")
    log("      primera corrida, y la fila que descarta es justo la que rompe la")
    log("      monotonia con nu. Eso es seleccion post-hoc, asi que se reporta la")
    log("      correlacion CON Y SIN el descarte y no se elige una.")
    if len(filas) >= 3:
        nu_t = np.array([f[1] for f in filas])
        hp_t = np.array([f[4] for f in filas])
        log("      corr(log nu, H_p) con TODAS las filas : %+.4f  (n = %d)"
            % (float(np.corrcoef(np.log(nu_t), hp_t)[0, 1]), len(filas)))
    if len(fiables) >= 3:
        nus = np.array([f[1] for f in fiables])
        hps = np.array([f[4] for f in fiables])
        rho = float(np.corrcoef(np.log(nus), hps)[0, 1])
        log("")
        log("  SOLO SOBRE LAS FIABLES:")
        log("  dispersion de H_p corregido : sd = %.4f  rango [%.3f, %.3f]"
            % (np.std(hps, ddof=1), hps.min(), hps.max()))
        log("  corr(log nu, H_p)           : %+.4f  (n = %d)" % (rho, len(fiables)))
        # [!] Con n <= 4 una correlacion NO ES EVIDENCIA. Para n = 3 hace falta
        # |r| > 0.997 para p < 0.05, y para n = 4, |r| > 0.95. Declarar que
        # "H_p correlaciona con nu" sobre tres capturas seria exactamente el
        # tipo de afirmacion que este proyecto ha tenido que retractar cuatro
        # veces. La rama de correlacion del Sec.2.2 exige un n que no hay.
        n_min_correlacion = 6
        if len(fiables) < n_min_correlacion:
            log("  -> n = %d < %d: la rama de correlacion del Sec.2.2 NO se puede"
                % (len(fiables), n_min_correlacion))
            log("     invocar. Con n=3 hace falta |r| > 0.997 para p < 0.05.")
            log("     SE APLICA LA TERCERA RAMA: la curva requerida va como BANDA")
            log("     [%.3f, %.3f] y el Sec.1 se lee contra la banda entera."
                % (min(0.5, hps.min()), hps.max()))
        elif hps.max() - hps.min() < 0.05:
            log("  -> COMPATIBLES: se usa el valor comun y la curva requerida queda fijada")
        elif abs(rho) > 0.8:
            log("  -> H_p CORRELACIONA CON nu: la superdifusion es efecto de la tasa de")
            log("     actividad; la curva requerida se construye por regimen de nu")
        else:
            log("  -> H_p VARIA SIN PATRON: la curva requerida se reporta como BANDA")
            log("     entre H_p = 0.5 y el maximo medido, y el Sec.1 se lee contra la banda")
    else:
        log("")
        log("  *** MENOS DE 3 FILAS FIABLES: la compuerta del Sec.2.2 NO SE PUEDE")
        log("      EVALUAR. La curva requerida se reporta como BANDA entre H_p = 0.5")
        log("      y el maximo de las filas fiables, y el Sec.1 se lee contra la banda")
        log("      entera -- que es el tercer desenlace del Sec.2.2, no un fallo. ***")
    log("")
    log("  [!] NUNCA se elige el H_p que hace cruzar las curvas (Sec.2.2). Si el")
    log("      Sec.1 cruza solo con el mas favorable, el resultado es NO DECIDIBLE.")
    return res


def cache_precio(d):
    """Precio de transaccion de la muestra cacheada de la v3.2."""
    return np.asarray(d["precio"], float)


def etapa_curvas(args):
    titulo("Sec.1 -- LAS DOS CURVAS")
    d = np.load(CACHE_V33)
    eps, cant, t = d["eps"], d["v"], d["t"]
    mid = 0.5 * (d["bid"] + d["ask"])
    n = eps.size
    pmed = float(np.median(d["precio"]))
    # Particion de la v3.2. El conjunto de PRUEBA no se abre (Sec.1.3 punto 3).
    embargo = 2233
    i1, i2 = int(n * 0.60), int(n * 0.80)
    log("  particion heredada de la v3.2 (hash de indices para el acta):")
    log("    entrenamiento [0, %d)   validacion [%d, %d)   PRUEBA [%d, %d) NO SE ABRE"
        % (i1, i1 + embargo, i2, i2 + embargo, n))

    flujo = np.concatenate(([0.0], np.cumsum(eps * cant)))
    nu = n / float(t[-1] - t[0])
    Hs = [60, 120, 300, 600, 900, 1800, 3600]

    # sigma_1 se recalcula DENTRO, por cada H_p, sobre la firma de esta misma
    # captura. Pasarlo por bandera invitaba a combinar un sigma_1 ajustado con
    # un H_p distinto, que es incoherente por construccion.
    fp = firma_pared(cache_precio(d), t,
                     sorted(set([0.5, 1, 2, 5] +
                                [round(x, 1) for x in np.logspace(
                                    np.log10(8), np.log10(1800), 28)])))
    hp_banda = [("banda baja", 0.371), ("difusivo", 0.500),
                ("v33 corregido", 0.552)]
    s1 = {}
    log("")
    log("  BANDA de H_p del Sec.2.2 (tercera rama). sigma_1 recalculada por cada uno:")
    for nom, hp in hp_banda:
        s1[hp] = sigma1_pb(fp, pmed, hp)
        log("    %-16s H_p = %.3f   sigma_1 = %.6f pb*s^(-H_p)" % (nom, hp, s1[hp]))
    log("")
    log("   H       n_vent   R2_medido  R2_barajado |     R2_req por H_p de la banda")
    log("                                          |   0.371       0.500       0.552")
    log("   " + "-" * 76)
    rng = np.random.default_rng(11)
    for H in Hs:
        idx = ventanas_no_solapadas(t, float(H))
        idx = idx[idx < n - 1]
        # objetivo: log-retorno del punto medio sobre la ventana siguiente
        fin = np.clip(np.searchsorted(t, t[idx] + H, side="left"), 0, n - 1)
        val = (fin > idx) & (mid[idx] > 0) & (mid[fin] > 0)
        idx, fin = idx[val], fin[val]
        y = np.log(mid[fin] / mid[idx])
        ret = [max(1, int(round(f * H * nu))) for f in (0.25, 0.5, 1.0, 2.0)]
        X = rasgos_flujo(flujo, idx, ret)
        ent = idx < i1
        vld = (idx >= i1 + embargo) & (idx < i2)
        if ent.sum() < 8 or vld.sum() < 4:
            log("   %5d s   %6d   (insuficiente)" % (H, int(vld.sum())))
            continue
        r2 = r2_fuera_de_muestra(X[ent], y[ent], X[vld], y[vld])
        # control barajado: se destruye el signo conservando el tamano
        eps_b = eps.copy()
        eps_b[:i2] = rng.permutation(eps_b[:i2])
        flujo_b = np.concatenate(([0.0], np.cumsum(eps_b * cant)))
        Xb = rasgos_flujo(flujo_b, idx, ret)
        r2b = r2_fuera_de_muestra(Xb[ent], y[ent], Xb[vld], y[vld])
        rq = [r2_requerido(float(H), s1[hp], hp) for _, hp in hp_banda]
        marca = ""
        if int(vld.sum()) < 30:
            marca = "  <- PILOTO"
        if H > ALCANCE_VERIFICADO_S:
            marca += "  EXTRAPOLA"
        log("   %5d s   %6d   %+8.4f    %+8.4f | %8.2f %% %8.2f %% %8.2f %%%s"
            % (H, int(vld.sum()), r2, r2b, 100 * rq[0], 100 * rq[1],
               100 * rq[2], marca))
        cruza = [nom for (nom, _), q in zip(hp_banda, rq) if r2 > q]
        if cruza and int(vld.sum()) >= 30 and r2 > 3 * abs(r2b):
            log("          -> CRUZA con: %s" % ", ".join(cruza))
    log("")
    log("  [!] n_vent es el numero de ventanas NO SOLAPADAS de VALIDACION. Con")
    log("      captura_v33 (16.4 h) y el 20 %% de validacion, por encima de 30 min")
    log("      la medicion es un piloto y NO decide (Sec.1.4).")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", default="hp")
    ap.add_argument("--hp", type=float, default=0.591)
    ap.add_argument("--sigma1", type=float, default=0.00058)
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    if a.etapa == "hp":
        etapa_hp(a)
    elif a.etapa == "curvas":
        etapa_curvas(a)
    else:
        log("etapa desconocida: %s" % a.etapa)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
