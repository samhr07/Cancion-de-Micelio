# -*- coding: utf-8 -*-
"""
tick_grande.py -- v3.3 Sec.1, Sec.2, Sec.3 y Sec.6 sobre dato ya en disco.

    python tick_grande.py --autotest     controles con verdad conocida
    python tick_grande.py                ejecuta Sec.6 -> Sec.1 -> Sec.2 -> Sec.3

NINGUN ESTADISTICO NUEVO para permanente/transitorio. El Sec.7 de la orden lo
prohibe y tiene razon: han muerto cuatro. Lo que hay aqui son rutas a
cantidades YA medidas, mas la correccion de bootstrap que puede invalidar
intervalos ya publicados.

[!] LA CONFUSION QUE HAY QUE DESHACER ANTES DE CALCULAR NADA
------------------------------------------------------------
La tabla del Sec.1 de la orden escribe `gamma_hat = 0.798` y lo etiqueta
"autocorrelacion de signos". Ese 0.798 es lo que la v3.2 midio, y es
`C(1) = corr(eps_t, eps_{t+1})`: la autocorrelacion **a rezago 1**.

Pero la relacion `H = (2-gamma)/2 - beta` viene del marco del propagador
(Bouchaud, Gefen, Potters & Wyart 2004; revision en Bouchaud, Farmer & Lillo
2009), donde `gamma` es el **EXPONENTE DE DECAIMIENTO** de

    C(l) = <eps_t eps_{t+l}>  ~  l^(-gamma)

Son dos numeros distintos y no hay ninguna razon para que coincidan: en renta
variable `C(1)` suele estar en 0.5-0.8 mientras el exponente ronda 0.5. Meter
`C(1)` donde va el exponente es el modo de fallo 1 del Sec.8 de la orden
("formula tomada del resumen sin verificar") con el 2 encima.

Aqui se mide el exponente por regresion log-log sobre el rango declarado, y se
reportan los dos numeros por separado para que no vuelvan a mezclarse.

VERIFICACION DE LA RELACION (Sec.1.1, punto 1)
----------------------------------------------
`H = 1 - gamma/2 - beta`, equivalentemente `2H = 2 - gamma - 2*beta`. Se
comprueba por sus dos limites, que es lo que la propia orden propone:

    beta = (1-gamma)/2  ->  H = 1/2                      (difusividad)
    beta = 0            ->  H = (2-gamma)/2              (impacto permanente)

RANGO DE VALIDEZ, y es la reserva que hay que llevar al reporte: la relacion es
ASINTOTICA y supone `0 < gamma < 1` con decaimiento en ley de potencias
sostenido sobre el rango de rezagos donde se mide `H`. Si `C(l)` no es una ley
de potencias limpia en ese rango, la relacion no aplica y el Sec.1 entero cae.
Por eso se reporta el R^2 del ajuste log-log junto al exponente.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

import propagador as P


CACHE = "telemetria/muestra_v32.npz"
TICK = 0.10                      # USD/BTC, tickSize de BTCUSDT en USD-M
COMISION_MAKER_VIP0 = 0.0002     # ASUMIDA, no leida (endpoint firmado)
COMISION_TAKER_VIP0 = 0.0005

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
    log("=" * 74)
    log(s)
    log("=" * 74)


# ===========================================================================
# Sec.6  -- N_eff bajo memoria larga. Corre ANTES que nada que use un IC.
# ===========================================================================

def escalado_de_varianza(x: np.ndarray, bloques=None) -> dict:
    """Mide `Var(media de bloque)` contra el tamano de bloque `b`.

    Bajo independencia `Var ~ b^(-1)`. Bajo memoria larga con
    `C(l) ~ l^(-gamma)` y `0 < gamma < 1`, `Var ~ b^(-gamma)`.

    [!] SE MIDE, NO SE SUPONE. La orden propone `N_eff = N^gamma` con la
    `gamma` de los signos, pero la serie cuyo IC hay que corregir es la de
    `dLL`, no la de signos, y no tienen por que compartir exponente. El
    exponente se estima de la propia serie y de ahi sale `N_eff`.
    """
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    if bloques is None:
        bloques = [b for b in (8, 16, 32, 64, 128, 256, 512, 1024,
                               2048, 4096, 8192, 16384) if b <= n // 8]
    filas = []
    for b in bloques:
        m = n // b
        med = x[: m * b].reshape(m, b).mean(axis=1)
        filas.append({"b": b, "n_bloques": m, "var": float(np.var(med, ddof=1))})
    lb = np.log([f["b"] for f in filas])
    lv = np.log([max(f["var"], 1e-300) for f in filas])
    A = np.column_stack([lb, np.ones_like(lb)])
    coef, *_ = np.linalg.lstsq(A, lv, rcond=None)
    pred = A @ coef
    ss = 1.0 - float(np.sum((lv - pred) ** 2) / np.sum((lv - lv.mean()) ** 2))
    gamma_eff = float(-coef[0])
    return {"filas": filas, "gamma_eff": gamma_eff, "r2": ss,
            "log_var_1": float(coef[1])}


def ic_memoria_larga(x: np.ndarray, z: float = 1.96) -> dict:
    """IC de la media extrapolando el escalado medido hasta `b = N`.

    En vez de remuestrear bloques de longitud fija -- que captura dependencia
    solo hasta esa longitud -- se ajusta `Var(media_b)` contra `b` y se evalua
    en `b = N`. Con `gamma_eff = 1` reproduce el IC clasico; con
    `gamma_eff < 1` lo ensancha por el factor `N^((1-gamma_eff)/2)`.
    """
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    esc = escalado_de_varianza(x)
    var_n = float(np.exp(esc["log_var_1"] + (-esc["gamma_eff"]) * np.log(n)))
    se = np.sqrt(max(var_n, 0.0))
    med = float(np.mean(x))
    n_eff = float(n ** esc["gamma_eff"])
    var_iid = float(np.var(x, ddof=1) / n)
    return {"media": med, "se": se, "ic": (med - z * se, med + z * se),
            "gamma_eff": esc["gamma_eff"], "r2": esc["r2"],
            "n": n, "n_eff": n_eff,
            "inflacion_vs_iid": float(np.sqrt(var_n / max(var_iid, 1e-300))),
            "escalado": esc}


# ===========================================================================
# Sec.1 -- gamma DE VERDAD, H, y beta implicita
# ===========================================================================

def acf_signos(eps: np.ndarray, rezagos) -> np.ndarray:
    """`C(l) = <eps_t eps_{t+l}>` por FFT, centrada."""
    e = np.asarray(eps, dtype=np.float64)
    e = e - e.mean()
    n = e.size
    m = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(e, m)
    ac = np.fft.irfft(F * np.conj(F), m)[:max(rezagos) + 1].real
    ac /= ac[0]
    return np.array([ac[l] for l in rezagos])


def gamma_de_signos(eps: np.ndarray, l_min: int = 10, l_max: int = 2000) -> dict:
    """Exponente de `C(l) ~ l^(-gamma)` por regresion log-log.

    [!] EL RANGO SE DECLARA. `l_min = 10` deja fuera la microestructura de los
    primeros rezagos (division de ordenes grandes) y `l_max = 2000` se queda
    dentro del rango donde `H` se midio. Se reportan tambien `C(1)` y el R^2
    del ajuste: si el R^2 es bajo, `C(l)` no es una ley de potencias limpia y
    la relacion del Sec.1 NO APLICA.
    """
    rez = np.unique(np.round(np.logspace(0, np.log10(l_max), 60)).astype(int))
    rez = rez[rez >= 1]
    C = acf_signos(eps, rez)
    sel = (rez >= l_min) & (rez <= l_max) & (C > 0)
    lx, ly = np.log(rez[sel].astype(float)), np.log(C[sel])
    A = np.column_stack([lx, np.ones_like(lx)])
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ coef
    r2 = 1.0 - float(np.sum((ly - pred) ** 2) / np.sum((ly - ly.mean()) ** 2))
    return {"gamma": float(-coef[0]), "r2": r2, "n_puntos": int(sel.sum()),
            "C1": float(acf_signos(eps, [1])[0]), "rezagos": rez, "C": C,
            "l_min": l_min, "l_max": l_max}


def H_de_firma(precio: np.ndarray, n_lo: int = 256, n_hi: int = 16384) -> dict:
    rej = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
    filas = P.firma_en_ticks(precio, rej, solapada=False)
    r = P.pendiente_en_ticks(filas, n_lo, n_hi)
    return {"H": r["H_p"], "pendiente": r["pendiente"], "filas": filas}


def beta_implicita(gamma: float, H: float) -> float:
    """`beta = (2 - gamma)/2 - H`, despejada de `H = (2-gamma)/2 - beta`."""
    return (2.0 - gamma) / 2.0 - H


def bootstrap_bloques_movil(x, estimador, longitud, n_sorteos, semilla=0):
    """Bootstrap por bloques moviles generico sobre una serie.

    [!] NO SIRVE PARA EXPONENTES DE MEMORIA LARGA, y se conserva aqui solo para
    poder EXHIBIR el fallo. Remuestrear bloques de longitud `b` produce una
    serie cuya dependencia muere en `b`, asi que un estimador que ajusta sobre
    rezagos mayores que `b` mide la longitud del bloque y no la serie. Medido
    sobre `captura_v33` con `b = N^(1/2) = 1015`:

        gamma  0.5217  ->  IC de bootstrap [0.7551, 0.9900]
        H      0.5907  ->  IC de bootstrap [0.0904, 0.1826]

    **Los intervalos no contienen el punto estimado**, que es la firma
    inequivoca de que los sorteos no estiman la misma cantidad. Usar
    `subsuestreo_contiguo` para estas cantidades.
    """
    x = np.asarray(x)
    n = x.size
    longitud = int(min(longitud, max(2, n // 4)))
    nb = int(np.ceil(n / longitud))
    rng = np.random.default_rng(semilla)
    out = np.empty(n_sorteos)
    for s in range(n_sorteos):
        ini = rng.integers(0, n - longitud + 1, size=nb)
        out[s] = estimador(np.concatenate([x[i:i + longitud] for i in ini])[:n])
    return out


def subsuestreo_contiguo(x, estimador, k: int) -> dict:
    """Parte la serie en `k` sub-series CONTIGUAS y estima en cada una.

    Es el sustituto correcto del bootstrap por bloques para exponentes de
    memoria larga: cada sub-serie conserva su estructura temporal intacta, asi
    que el estimador mide en ella lo mismo que mide en el total. El precio es
    que hay pocas replicas y que cada una es mas corta -- con sesgo propio de
    tamano finito, que por eso se reporta junto a la dispersion.

    [!] El intervalo `t` que se devuelve supone que las sub-series son
    independientes entre si. Con memoria larga NO lo son del todo, asi que es
    OPTIMISTA y va marcado como tal. No se convierte en un IC formal del
    estimador de muestra completa: la tasa de convergencia de un exponente de
    memoria larga no es la raiz de n y despejarla exigiria supuestos que este
    documento no tiene.
    """
    x = np.asarray(x)
    n = x.size
    m = n // k
    est = np.array([estimador(x[i * m:(i + 1) * m]) for i in range(k)])
    med = float(np.mean(est))
    sd = float(np.std(est, ddof=1)) if k > 1 else float("nan")
    se = sd / np.sqrt(k) if k > 1 else float("nan")
    return {"estimaciones": est, "media": med, "sd": sd, "se": se,
            "ic_t_optimista": (med - 2.31 * se, med + 2.31 * se),
            "k": k, "n_por_subserie": m,
            "subpotenciado": bool(k < 15)}


# ===========================================================================
# Sec.2 -- eta de Robert-Rosenbaum, la compuerta de tick grande
# ===========================================================================

def eta_robert_rosenbaum(precio: np.ndarray, tick: float = TICK,
                         solo_un_tick: bool = True) -> dict:
    """`eta = N_c / (2 * N_a)`.  Robert & Rosenbaum (2011), zonas de incertidumbre.

    [!] EL EVENTO DE SALTO SE DEFINE ANTES DE CONTAR (modo de fallo 3 del Sec.8):
    se toma la subsucesion de cambios de precio NO NULOS, expresados en ticks y
    redondeados. Una *continuacion* es un par de cambios consecutivos con el
    MISMO signo; una *alternacion*, con signo opuesto. Con `solo_un_tick` se
    exige ademas que los dos cambios sean de magnitud exactamente 1 tick, que es
    la version del articulo: contar saltos de mas de un tick infla `eta` de forma
    espuria.

    `eta < 1/2` -> regimen de tick grande.
    """
    p = np.asarray(precio, dtype=np.float64)
    d = np.round(np.diff(p) / tick).astype(np.int64)
    d = d[d != 0]
    if d.size < 3:
        return {"eta": float("nan"), "N_c": 0, "N_a": 0, "n_cambios": int(d.size)}
    a, b = d[:-1], d[1:]
    if solo_un_tick:
        m = (np.abs(a) == 1) & (np.abs(b) == 1)
        a, b = a[m], b[m]
    N_c = int(np.sum(np.sign(a) == np.sign(b)))
    N_a = int(np.sum(np.sign(a) != np.sign(b)))
    eta = N_c / (2.0 * N_a) if N_a > 0 else float("inf")
    return {"eta": float(eta), "N_c": N_c, "N_a": N_a,
            "n_cambios": int(d.size), "solo_un_tick": solo_un_tick}


# ===========================================================================
# Sec.3 -- la aritmetica del tick, EN PUNTOS BASICOS
# ===========================================================================

def tabla_costes_pb(precio: float, tick: float = TICK) -> dict:
    """Costes en pb al precio REAL del tramo, no en USD/BTC heredado.

    [!] El Sec.3 de la orden sospecha que el tramo podria estar cerca de 96 000,
    con lo que `c(u)` seria 38.4 y el factor del paso 3 pasaria de 3.46 a 5.1.
    Se resuelve midiendo el precio del tramo en vez de suponerlo.
    """
    pb = lambda usd: 1e4 * usd / precio
    c_mm = 2 * COMISION_MAKER_VIP0 * precio
    c_mt = (COMISION_MAKER_VIP0 + COMISION_TAKER_VIP0) * precio
    c_tt = 2 * COMISION_TAKER_VIP0 * precio
    return {"precio": precio, "tick_usd": tick, "tick_pb": pb(tick),
            "maker_maker": c_mm, "maker_taker": c_mt, "taker_taker": c_tt,
            "maker_maker_pb": pb(c_mm), "ticks_por_ida_y_vuelta": c_mm / tick}


# ===========================================================================
# Controles
# ===========================================================================

def _fgn(n: int, H: float, rng) -> np.ndarray:
    """Ruido gaussiano fraccionario por sintesis espectral.

    Densidad espectral `S(f) ~ f^(1-2H)`, que da `C(l) ~ l^(2H-2)`. Se usa
    SOLO en los controles: da una ACF con exponente conocido de antemano, que
    es lo que un control positivo necesita.
    """
    m = 1 << int(np.ceil(np.log2(4 * n)))
    f = np.fft.rfftfreq(m)[1:]
    amp = f ** ((1.0 - 2.0 * H) / 2.0)
    fase = rng.normal(size=f.size) + 1j * rng.normal(size=f.size)
    espectro = np.concatenate(([0.0 + 0j], amp * fase))
    x = np.fft.irfft(espectro, m)[:n]
    return x / np.std(x)


def _autotest() -> int:
    fallos = 0

    def ok(nombre, cond, extra=""):
        nonlocal fallos
        log("  [%s] %-58s %s" % ("OK  " if cond else "FALLO", nombre, extra))
        if not cond:
            fallos += 1

    log("== 1. La relacion H = (2-gamma)/2 - beta en sus dos limites ==")
    for g in (0.2, 0.5, 0.8):
        H_dif = (2 - g) / 2 - (1 - g) / 2
        H_per = (2 - g) / 2 - 0.0
        log("     gamma=%.1f  beta=(1-g)/2 -> H=%.4f (debe ser 0.5) | beta=0 -> H=%.4f"
            % (g, H_dif, H_per))
        ok("limite de difusividad con gamma=%.1f" % g, abs(H_dif - 0.5) < 1e-12)
        ok("limite permanente con gamma=%.1f" % g, abs(H_per - (2 - g) / 2) < 1e-12)

    log("== 2. gamma se recupera de una ACF de signos con exponente CONOCIDO ==")
    # [!] La primera version de este control usaba una suma de seis AR(1) y
    # daba R2 = 0.87: una mezcla FINITA de escalas es ley de potencias solo a
    # trozos, con ondulaciones entre las tau. El fallo era del GENERADOR, no
    # del estimador, y bajar el umbral para que pasara habria sido ajustar el
    # control al resultado. Se sustituye por ruido gaussiano fraccionario, cuyo
    # exponente SI es verdad conocida: `C(l) ~ l^(2H-2)`, o sea `gamma = 2-2H`.
    #
    # Tomar el signo cambia la ACF por la ley del arcoseno,
    # `corr(sgn) = (2/pi)*arcsin(rho)`, pero para `rho` pequeno eso es
    # proporcional a `rho`, asi que el EXPONENTE se conserva a rezagos grandes
    # -- que es justo lo que este control mide.
    rng = np.random.default_rng(4)
    r2s = {}
    for H_v in (0.75, 0.65):
        g_v = 2.0 - 2.0 * H_v
        eps = np.sign(_fgn(1 << 19, H_v, rng))
        g = gamma_de_signos(eps, l_min=10, l_max=2000)
        r2s[H_v] = g["r2"]
        log("     H=%.2f -> gamma verdadera %.2f | medida %.4f | R2 %.4f | C(1) %.4f"
            % (H_v, g_v, g["gamma"], g["r2"], g["C1"]))
        ok("gamma recuperada con H=%.2f (error < 0.10)" % H_v,
           abs(g["gamma"] - g_v) < 0.10)
        ok("gamma != C(1) con H=%.2f: son cantidades distintas" % H_v,
           abs(g["gamma"] - g["C1"]) > 0.1,
           "gamma=%.3f  C(1)=%.3f" % (g["gamma"], g["C1"]))
    # [!] HALLAZGO DEL CONTROL, y cambia como hay que leer el R2 sobre dato
    # real: con memoria FUERTE (H=0.75) el ajuste log-log da R2 = 0.977; con
    # memoria debil (H=0.65) cae a 0.66 **con la gamma igual de bien
    # recuperada**. La ACF verdadera SIGUE siendo una ley de potencias exacta
    # en los dos casos -- lo que cambia es que `C(l)` es mas pequena y el ruido
    # muestral domina el logaritmo a rezagos largos.
    #
    # Consecuencia declarada: **un R2 bajo NO es prueba de que la relacion del
    # Sec.1 no aplique.** Se reporta como diagnostico, no como compuerta. Poner
    # una compuerta en R2 > 0.90 habria descartado un caso donde el estimador
    # funciona perfectamente.
    ok("el R2 cae con la memoria debil aunque gamma siga siendo correcta",
       r2s[0.75] > 0.95 > r2s[0.65],
       "R2(H=0.75) = %.3f > R2(H=0.65) = %.3f" % (r2s[0.75], r2s[0.65]))

    log("== 3. El escalado de varianza distingue iid de memoria larga ==")
    x_iid = rng.normal(size=200000)
    e_iid = escalado_de_varianza(x_iid)
    log("     iid          -> gamma_eff = %.4f  (debe ser ~1)  R2 = %.4f"
        % (e_iid["gamma_eff"], e_iid["r2"]))
    ok("iid da gamma_eff ~ 1", abs(e_iid["gamma_eff"] - 1.0) < 0.10)

    x_mem = np.cumsum(rng.normal(size=200000))
    x_mem = x_mem - np.mean(x_mem)
    e_mem = escalado_de_varianza(x_mem)
    log("     paseo (memoria maxima) -> gamma_eff = %.4f  (debe ser << 1)"
        % e_mem["gamma_eff"])
    ok("el paseo da gamma_eff muy por debajo de 1", e_mem["gamma_eff"] < 0.5)

    ic_i = ic_memoria_larga(x_iid)
    log("     IC iid: inflacion sobre el iid clasico = %.3f (debe ser ~1)"
        % ic_i["inflacion_vs_iid"])
    ok("sobre iid el IC no se ensancha", 0.6 < ic_i["inflacion_vs_iid"] < 1.7)

    log("== 4. eta de Robert-Rosenbaum sobre verdad conocida ==")
    # Precio que ALTERNA siempre -> N_c = 0 -> eta = 0 (tick grande extremo).
    p_alt = 100.0 + TICK * np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1] * 200)
    r_alt = eta_robert_rosenbaum(p_alt)
    log("     alternancia pura -> eta = %.4f  (N_c=%d N_a=%d)"
        % (r_alt["eta"], r_alt["N_c"], r_alt["N_a"]))
    ok("alternancia pura da eta = 0", abs(r_alt["eta"]) < 1e-12)
    # Precio que CONTINUA siempre -> N_a = 0 -> eta infinita (sin zona).
    p_con = 100.0 + TICK * np.arange(2000.0)
    r_con = eta_robert_rosenbaum(p_con)
    log("     continuacion pura -> eta = %s  (N_c=%d N_a=%d)"
        % (r_con["eta"], r_con["N_c"], r_con["N_a"]))
    ok("continuacion pura no da tick grande", not (r_con["eta"] < 0.5))
    # Mezcla 50/50 -> N_c = N_a -> eta = 0.5, justo la frontera.
    rng2 = np.random.default_rng(9)
    pasos = rng2.choice([-1.0, 1.0], size=200000)
    r_mix = eta_robert_rosenbaum(100.0 + TICK * np.cumsum(pasos))
    log("     signos iid -> eta = %.4f  (la frontera teorica es 0.5)" % r_mix["eta"])
    ok("signos iid caen en la frontera 0.5", abs(r_mix["eta"] - 0.5) < 0.03)

    log("== 5. La tabla de costes escala con el precio ==")
    t65, t96 = tabla_costes_pb(65100.0), tabla_costes_pb(96000.0)
    log("     a 65 100: c_mm = %.2f USD/BTC = %.2f pb | a 96 000: %.2f USD/BTC = %.2f pb"
        % (t65["maker_maker"], t65["maker_maker_pb"],
           t96["maker_maker"], t96["maker_maker_pb"]))
    ok("los pb NO dependen del precio (4 pb en los dos)",
       abs(t65["maker_maker_pb"] - t96["maker_maker_pb"]) < 1e-9)
    ok("los USD/BTC SI dependen del precio",
       abs(t65["maker_maker"] - t96["maker_maker"]) > 10.0)

    log("")
    log("RESULTADO: %d fallo(s)" % fallos)
    return fallos


# ===========================================================================
# Ejecucion
# ===========================================================================

def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()

    d = np.load(CACHE)
    eps, precio, t = d["eps"], d["precio"], d["t"]
    n = eps.size
    nu = n / float(t[-1] - t[0])
    # Misma particion que la v3.2. El conjunto de PRUEBA no se reabre (Sec.7).
    embargo = 2233
    i1 = int(n * 0.60)
    ent = np.arange(0, i1)

    titulo("Sec.6 -- N_eff BAJO MEMORIA LARGA (corre antes que cualquier IC)")
    log("Se MIDE el escalado Var(media de bloque) contra b, no se supone.")
    for etiq, serie in (("signos eps (entrenamiento)", eps[ent].astype(float)),
                        ("incrementos de precio", np.diff(precio[ent]))):
        r = ic_memoria_larga(serie)
        log("")
        log("  %s" % etiq)
        log("    gamma_eff = %.4f   R2 del escalado = %.4f" % (r["gamma_eff"], r["r2"]))
        log("    N = %d   N_eff = N^gamma_eff = %.0f   perdida = %.1fx"
            % (r["n"], r["n_eff"], r["n"] / max(r["n_eff"], 1.0)))
        log("    inflacion del error estandar sobre el iid clasico = %.2fx"
            % r["inflacion_vs_iid"])

    log("")
    log("  --- RECALCULO DEL IC DEL PASO 2 DE LA v3.2 (Sec.6 lo exige) ---")
    log("  [!] Esto toca el conjunto de prueba. El Sec.7 prohibe REABRIRLO para un")
    log("      veredicto nuevo; el Sec.6 ORDENA recalcular un IC ya publicado. Se")
    log("      recomputa exactamente el mismo estadistico, sin ajustar ni elegir")
    log("      nada, y no se deriva ninguna decision nueva de el.")
    import migracion_v32 as M
    from experimento_v32 import ll_media
    K = 2 * embargo
    i2 = int(n * 0.80)
    pru = np.arange(i2 + embargo, n)
    y, v = d["y_mid"], d["v"]
    m0 = M.ajustar_M0(y[ent])
    m2 = M.ajustar(y[ent], eps[ent], v[ent], K, 0.0, con_suelo=True)
    _, l0 = ll_media(m0, y[pru], eps[pru], v[pru], K)
    _, l2 = ll_media(m2, y[pru], eps[pru], v[pru], K)
    dll = l2[:min(len(l0), len(l2))] - l0[:min(len(l0), len(l2))]
    r = ic_memoria_larga(dll)
    bs = M.bootstrap_bloques(dll, 5 * embargo, 500, semilla=3)
    log("    serie dLL de prueba: N = %d" % dll.size)
    log("    gamma_eff medido = %.4f  (R2 %.4f)  ->  N_eff = %.0f  (perdida %.1fx)"
        % (r["gamma_eff"], r["r2"], r["n_eff"], r["n"] / max(r["n_eff"], 1.0)))
    log("    IC publicado en la v3.2 (bloques de 5*embargo) : [%+.3e, %+.3e]"
        % bs["ic"])
    log("    IC corregido por memoria larga                 : [%+.3e, %+.3e]"
        % r["ic"])
    log("    ensanchamiento del error estandar : %.2fx" % r["inflacion_vs_iid"])
    log("    EL PASO 2 %s con el IC corregido"
        % ("SIGUE PASANDO (el IC sigue excluyendo 0)" if r["ic"][0] > 0
           else "*** YA NO PASA: el IC corregido incluye 0 ***"))

    titulo("Sec.1 -- gamma DE VERDAD, H, y beta implicita")
    log("[!] `C(1)` y el EXPONENTE de `C(l) ~ l^(-gamma)` son cantidades distintas.")
    log("    La tabla del Sec.1 de la orden mete `C(1) = 0.798` donde va el exponente.")
    log("")
    g = gamma_de_signos(eps[ent])
    log("  C(1) medida (lo que la v3.2 llamo gamma_hat) : %+.4f" % g["C1"])
    log("  EXPONENTE gamma de C(l) ~ l^(-gamma)         : %+.4f" % g["gamma"])
    log("  rango del ajuste: rezagos [%d, %d], %d puntos, R2 = %.4f"
        % (g["l_min"], g["l_max"], g["n_puntos"], g["r2"]))
    log("")
    log("  C(l) medida:")
    for l, c in zip(g["rezagos"], g["C"]):
        if l in (1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2000):
            log("    l = %5d   C = %+.5f" % (l, c))
    log("")
    log("  [!] El R2 se reporta como DIAGNOSTICO, no como compuerta. El control 2")
    log("      lo mide: con memoria debil el R2 cae a 0.66 con gamma igual de bien")
    log("      recuperada, porque C(l) es pequena y el ruido domina el logaritmo.")
    log("      Un R2 bajo NO prueba que la relacion del Sec.1 no aplique.")

    h = H_de_firma(precio[ent])
    log("")
    log("  H_p de la firma de volatilidad (n en [256, 16384]) = %.4f" % h["H"])

    b_imp = beta_implicita(g["gamma"], h["H"])
    b_dif = (1.0 - g["gamma"]) / 2.0
    log("  beta implicita = (2-gamma)/2 - H = %+.4f" % b_imp)
    log("  beta de difusividad = (1-gamma)/2 = %+.4f" % b_dif)

    log("")
    log("  [!] EL BOOTSTRAP POR BLOQUES NO SIRVE AQUI, y se exhibe el fallo:")
    L = int(round(n ** 0.5))
    gb = bootstrap_bloques_movil(
        eps[ent].astype(float),
        lambda s: gamma_de_signos(s)["gamma"], L, 40, semilla=31)
    hb = bootstrap_bloques_movil(
        precio[ent], lambda s: H_de_firma(s)["H"], L, 40, semilla=32)
    q = lambda v: np.percentile(v, [2.5, 97.5])
    log("      con bloque = N^(1/2) = %d:" % L)
    log("        gamma  %.4f -> IC de bootstrap [%.4f, %.4f]" % (g["gamma"], *q(gb)))
    log("        H      %.4f -> IC de bootstrap [%.4f, %.4f]" % (h["H"], *q(hb)))
    log("      Los intervalos NO CONTIENEN el punto estimado. Remuestrear bloques")
    log("      de longitud b da una serie cuya dependencia muere en b, y los dos")
    log("      estimadores ajustan sobre rezagos MAYORES que b (hasta 2000 y")
    log("      16384). Miden la longitud del bloque, no el mercado.")
    log("")
    log("  DISPERSION POR SUBMUESTREO CONTIGUO (el sustituto correcto):")
    k = 8
    sg = subsuestreo_contiguo(eps[ent].astype(float),
                              lambda s: gamma_de_signos(s)["gamma"], k)
    sh = subsuestreo_contiguo(precio[ent], lambda s: H_de_firma(s)["H"], k)
    sb = (2.0 - sg["estimaciones"]) / 2.0 - sh["estimaciones"]
    log("    %d sub-series contiguas de %d ticks cada una %s"
        % (k, sg["n_por_subserie"],
           "*** SUBPOTENCIADO (<15) ***" if sg["subpotenciado"] else ""))
    log("    gamma por sub-serie : " + " ".join("%.3f" % v for v in sg["estimaciones"]))
    log("    H     por sub-serie : " + " ".join("%.3f" % v for v in sh["estimaciones"]))
    log("    beta  por sub-serie : " + " ".join("%+.3f" % v for v in sb))
    seb = float(np.std(sb, ddof=1)) / np.sqrt(k)
    lo, hi = float(np.mean(sb) - 2.31 * seb), float(np.mean(sb) + 2.31 * seb)
    log("    gamma  media %.4f  sd %.4f   (muestra completa %.4f)"
        % (sg["media"], sg["sd"], g["gamma"]))
    log("    H      media %.4f  sd %.4f   (muestra completa %.4f)"
        % (sh["media"], sh["sd"], h["H"]))
    log("    beta   media %+.4f  sd %.4f  intervalo t OPTIMISTA [%+.4f, %+.4f]"
        % (float(np.mean(sb)), float(np.std(sb, ddof=1)), lo, hi))
    log("    [!] El intervalo supone sub-series independientes. Con memoria larga")
    log("        no lo son, asi que es una COTA INFERIOR de la incertidumbre.")
    log("")
    log("  CONTRASTES DECLARADOS (Sec.1.1 punto 3):")
    log("    intervalo de beta implicita contiene 0 (permanencia)?   %s"
        % ("SI" if lo <= 0 <= hi else "NO"))
    log("    contiene (1-gamma)/2 = %+.4f (difusividad)?  %s"
        % (b_dif, "SI" if lo <= b_dif <= hi else "NO"))
    log("    [!] Con un intervalo optimista, un 'NO' es debil: la incertidumbre")
    log("        real es mayor y podria cubrir cualquiera de los dos.")
    log("")
    log("  LOS TRES VALORES DE beta, JUNTOS Y SIN ELEGIR (Sec.1.1 punto 4):")
    log("    implicita por H y gamma          %+.4f" % b_imp)
    log("    ajuste M2 sobre punto medio      -0.1600")
    log("    ajuste M2 sobre precio de trans. +4.6900")
    log("    La discrepancia ES el resultado.")

    titulo("Sec.2 -- eta DE ROBERT-ROSENBAUM: compuerta de tick grande")
    for etiq, s1 in (("solo saltos de 1 tick (la del articulo)", True),
                     ("todos los saltos (control del modo de fallo 3)", False)):
        r = eta_robert_rosenbaum(precio, solo_un_tick=s1)
        log("  %-46s eta = %.4f   N_c = %d   N_a = %d"
            % (etiq, r["eta"], r["N_c"], r["N_a"]))
    r = eta_robert_rosenbaum(precio)
    pmed = float(np.median(precio))
    log("")
    log("  cambios de precio no nulos: %d de %d ticks (%.1f %%)"
        % (r["n_cambios"], len(precio) - 1,
           100.0 * r["n_cambios"] / (len(precio) - 1)))
    log("  VEREDICTO: eta = %.4f  ->  %s"
        % (r["eta"], "TICK GRANDE (Sec.4 y Sec.5 aplican)" if r["eta"] < 0.5
           else "el marco NO aplica: Sec.4 y Sec.5 no se ejecutan"))
    log("")
    log("  [!] RESERVA sobre el tick RELATIVO (Sec.2 de la orden):")
    log("      alpha/P = %.2f / %.0f = %.3e = %.4f pb por tick"
        % (TICK, pmed, TICK / pmed, 1e4 * TICK / pmed))
    log("      En renta variable 'tick grande' significa alpha/P de 1e-4 a 1e-3:")
    log("      DOS O TRES ORDENES DE MAGNITUD MAYOR. Que la horquilla clavada en")
    log("      1 tick y un alpha/P diminuto coexistan no esta en los articulos.")

    titulo("Sec.3 -- LA ARITMETICA DEL TICK, EN PUNTOS BASICOS")
    tb = tabla_costes_pb(pmed)
    log("  precio MEDIANO del tramo   : %.1f USD/BTC" % pmed)
    log("  rango del tramo            : %.1f -- %.1f" % (precio.min(), precio.max()))
    log("  -> la sospecha del Sec.3 (que el tramo corriera cerca de 96 000) NO se")
    log("     cumple: c(u) = %.2f USD/BTC es correcta para ESTE tramo." % tb["maker_maker"])
    log("")
    log("  magnitud                        USD/BTC        pb       ticks")
    log("  1 tick                          %8.3f  %8.4f  %10.1f" % (TICK, tb["tick_pb"], 1.0))
    log("  comision ida y vuelta maker     %8.3f  %8.4f  %10.1f"
        % (tb["maker_maker"], tb["maker_maker_pb"], tb["ticks_por_ida_y_vuelta"]))
    log("  |mu| q90 medida (v3.2)          %8.3f  %8.4f  %10.1f"
        % (11.297472, 1e4 * 11.297472 / pmed, 11.297472 / TICK))
    log("  umbral 1.5*c(u)                 %8.3f  %8.4f  %10.1f"
        % (1.5 * tb["maker_maker"], 1.5 * tb["maker_maker_pb"],
           1.5 * tb["maker_maker"] / TICK))
    log("")
    log("  factor que falta al paso 3 : %.2f"
        % (1.5 * tb["maker_maker"] / 11.297472))
    log("  Hay que predecir %.0f ticks de movimiento para pagar la ida y vuelta."
        % (1.5 * tb["maker_maker"] / TICK))
    log("  Un agotamiento de cola da 1 tick: la microestructura opera %.0fx por"
        % (1.5 * tb["maker_maker"] / TICK))
    log("  debajo del umbral de rentabilidad.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
