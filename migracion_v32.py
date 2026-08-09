# -*- coding: utf-8 -*-
"""
migracion_v32.py -- estimador de los modelos anidados M0 / M1 / M1' / M2 / M2-osc.

    python migracion_v32.py --autotest        controles con verdad conocida
    python migracion_v32.py --compuerta       evalua la compuerta y NO mide nada mas

⚠ ESTE MODULO NO EJECUTA LA DECISION. `PREREGISTRO_3_2.md` (commit 190edda) prohibe
ejecutar el §3 en adelante hasta que la compuerta de datos pase: >= 375*embargo
ticks continuos y limpios (731 250 con el piso de 1950). Lo que hay aqui es el
estimador y sus controles positivos, validados contra verdad conocida y sin tocar
la captura -- llegar a los 731 k ticks con el estimador sin validar seria el
desperdicio.

EL ESPACIO DE MODELOS
---------------------
    Delta p_t = SUM_k h(k) * x_{t-k} + eta_t,    x_t = eps_t * f(v_t)
    G(tau) = G0 * (1 + tau/tau0)^(-beta)         f(v) = v^delta
    h(0) = G(0),   h(k) = G(k) - G(k-1)

| id      | restriccion                | libres |
|---------|----------------------------|--------|
| M0      | G0 = 0                     | 0      |
| M1      | beta = 0, delta = 0.5      | 1      |
| M2      | ninguna                    | 4      |
| M2-osc  | + cos(omega_G*tau + phi)   | 5      |

`beta = 0` deja `G(tau) = G0` constante: impacto **instantaneo y permanente**,
o sea `h = [G0, 0, 0, ...]`. Esa es exactamente la especificacion de MkII, y por
eso M1 es un punto del mismo espacio y no una arquitectura rival.

DECISIONES NUMERICAS
--------------------
`G0` entra LINEALMENTE en la prediccion, asi que se perfila por minimos cuadrados
en cada evaluacion y solo se optimiza sobre `(log tau0, beta)`. Eso convierte un
problema de 3 parametros mal escalado en uno de 2 bien condicionado, y de paso da
el error estandar de `G0` sin trabajo extra.

Se optimiza sobre `log tau0` y no sobre `tau0`: la verosimilitud es casi plana en
tau0 grande y un optimizador sobre la escala lineal se queda donde lo dejen.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
from scipy import optimize


# Piso del embargo declarado en el preregistro Sec. 2.2 (50 s x 39 tx/s).
EMBARGO_PISO = 1950
# Bloques efectivos minimos para que un IC de bootstrap sea reportable.
BLOQUES_MINIMOS = 15


def ticks_de_compuerta(embargo: int = EMBARGO_PISO, f_prueba: float = 0.20,
                       bloques_min: int = BLOQUES_MINIMOS) -> int:
    """Ticks continuos y limpios que exige la compuerta. Se DERIVA, no se escribe.

    La banda de embargo se DESCARTA del conjunto de prueba, no se reparte, asi que

        (f_prueba * N  -  embargo)  >=  bloques_min * 5 * embargo
        =>  N  >=  (bloques_min*5 + 1) * embargo / f_prueba

    Con los valores declarados: 380 * embargo = 741 000 ticks.

    ⚠ La primera version de este numero decia 375*embargo = 731 250 por olvidar
    el embargo descontado, y daba **14 bloques, no 15**. Lo caza el control 6 de
    este modulo. Es exactamente el defecto que el preregistro obliga a marcar:
    presentar 12 bloques como 15.
    """
    return int(np.ceil((bloques_min * 5 + 1) * embargo / f_prueba))
# Rejilla de delta del §5.1. La celda 0.5 es la de MkII y se reporta siempre.
REJILLA_DELTA = (0.0, 0.25, 0.5, 1.0)


# ===========================================================================
# 1. Nucleo
# ===========================================================================

def nucleo_G(tau, tau0: float, beta: float) -> np.ndarray:
    """G(tau)/G0 = (1 + tau/tau0)^(-beta). Normalizado a G(0) = 1."""
    return (1.0 + np.asarray(tau, dtype=np.float64) / tau0) ** (-beta)


def nucleo_h(K: int, tau0: float, beta: float,
             omega_G: float = 0.0, phi: float = 0.0) -> np.ndarray:
    """Respuesta al impulso de los INCREMENTOS: h(0)=G(0), h(k)=G(k)-G(k-1).

    Con `beta > 0` sale h(0) = 1 y h(k) < 0: impacto instantaneo seguido de
    reversion mientras el impacto transitorio se relaja. Con `beta = 0` sale
    h = [1, 0, 0, ...]: instantaneo y permanente, que es M1.
    """
    tau = np.arange(K + 1, dtype=np.float64)
    G = nucleo_G(tau, tau0, beta)
    if omega_G != 0.0:
        G = G * np.cos(omega_G * tau + phi)
    h = np.empty(K + 1, dtype=np.float64)
    h[0] = G[0]
    h[1:] = G[1:] - G[:-1]
    return h


def convolucion_causal(x: np.ndarray, h: np.ndarray) -> np.ndarray:
    """SUM_k h(k) x_{t-k}, misma longitud que x y SIN mirar al futuro.

    ⚠ `np.convolve(..., mode="same")` CENTRA el nucleo y mete el futuro en el
    presente. Sobre un test de predictibilidad eso fabrica la senal que se
    pretende medir; costo una sesion entera en la v3.1 §3.
    """
    if x.size * h.size > 2_000_000:
        from scipy import signal as _sig
        return _sig.fftconvolve(x, h)[: x.size]
    return np.convolve(x, h)[: x.size]


def decaimiento(mod: dict) -> float:
    """G(K)/G(0): cuanto queda del impacto al final del rango ajustado.

    ⚠ ESTA es la cantidad identificada, y `beta` NO lo es por si sola. Con
    `tau0` grande, `(1 + tau/tau0)^(-5)` vale ~1 en todo el rango observado y
    es indistinguible de `beta = 0`. Lo que los datos determinan es la CURVA
    `G(tau)` sobre `[0, K]`, no la pareja `(tau0, beta)` que la genera.

    Encontrado por el control positivo 2 de este mismo modulo: con verdad
    `beta = 0` el ajuste devolvia `beta = 5.0` (el borde) y acertaba igual la
    curva. Ver la nota del §5.2 en `PREREGISTRO_3_2.md`.

      decaimiento ~ 1  ->  impacto permanente  ->  M1 / MkII
      decaimiento < 1  ->  impacto transitorio ->  propagador
    """
    return float(nucleo_G(mod["K"], mod["tau0"], mod["beta"]))


# ===========================================================================
# 2. Ajuste
# ===========================================================================

def _perfilar_G0(y: np.ndarray, z: np.ndarray) -> tuple[float, np.ndarray]:
    """G0 optimo por minimos cuadrados dado el regresor convolucionado `z`."""
    zz = float(np.dot(z, z))
    if zz <= 0.0:
        return 0.0, np.zeros_like(y)
    G0 = float(np.dot(y, z) / zz)
    return G0, G0 * z


def forzamiento(eps: np.ndarray, v: np.ndarray, delta: float,
                v_mediana: float | None = None) -> np.ndarray:
    """x_t = eps_t * (v_t / mediana(v))^delta.

    La normalizacion por la mediana del tramo es del §5.1 del preregistro: sin
    ella `G0` y `delta` quedan confundidos por la escala y sus errores estandar
    no son interpretables.
    """
    if v_mediana is None:
        v_mediana = float(np.median(v[v > 0])) if np.any(v > 0) else 1.0
    if delta == 0.0:
        return np.asarray(eps, dtype=np.float64).copy()
    return np.asarray(eps, np.float64) * (np.asarray(v, np.float64) / v_mediana) ** delta


def ajustar(y: np.ndarray, eps: np.ndarray, v: np.ndarray, K: int,
            delta: float, beta_libre: bool = True,
            con_oscilacion: bool = False, v_mediana: float | None = None) -> dict:
    """Ajusta el modelo sobre (y, eps, v). Devuelve parametros y sigma."""
    x = forzamiento(eps, v, delta, v_mediana)
    y = np.asarray(y, dtype=np.float64)

    def sse(par):
        if con_oscilacion:
            ltau0, beta, om, ph = par
        elif beta_libre:
            ltau0, beta = par
            om = ph = 0.0
        else:
            (ltau0,) = par
            beta = 0.0
            om = ph = 0.0
        tau0 = float(np.exp(ltau0))
        if not np.isfinite(tau0) or tau0 <= 1e-3 or beta < -1.0 or beta > 5.0:
            return 1e30
        h = nucleo_h(K, tau0, beta, om, ph)
        z = convolucion_causal(x, h)
        _, pred = _perfilar_G0(y, z)
        return float(np.sum((y - pred) ** 2))

    if not beta_libre and not con_oscilacion:
        # M1: beta = 0 hace G constante, asi que tau0 es inobservable. Se fija.
        tau0, beta, om, ph = 1.0, 0.0, 0.0, 0.0
    else:
        p0 = [np.log(50.0), 0.4] + ([0.01, 0.0] if con_oscilacion else [])
        r = optimize.minimize(sse, p0, method="Nelder-Mead",
                              options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-10})
        if con_oscilacion:
            tau0, beta, om, ph = float(np.exp(r.x[0])), float(r.x[1]), float(r.x[2]), float(r.x[3])
        else:
            tau0, beta = float(np.exp(r.x[0])), float(r.x[1])
            om = ph = 0.0

    h = nucleo_h(K, tau0, beta, om, ph)
    z = convolucion_causal(x, h)
    G0, pred = _perfilar_G0(y, z)
    resid = y - pred
    return {"G0": G0, "tau0": tau0, "beta": beta, "delta": delta,
            "omega_G": om, "phi": ph, "K": K,
            "sigma": float(np.std(resid, ddof=1)),
            "v_mediana": float(np.median(v[v > 0])) if v_mediana is None else v_mediana,
            "media_y": float(np.mean(y))}


def ajustar_M0(y: np.ndarray) -> dict:
    """El nulo: sin capa predictiva. La prediccion es la media."""
    y = np.asarray(y, dtype=np.float64)
    return {"G0": 0.0, "tau0": np.nan, "beta": np.nan, "delta": np.nan,
            "omega_G": 0.0, "phi": 0.0, "K": 0,
            "sigma": float(np.std(y, ddof=1)), "v_mediana": 1.0,
            "media_y": float(np.mean(y))}


def predecir(mod: dict, eps: np.ndarray, v: np.ndarray) -> np.ndarray:
    if mod["G0"] == 0.0:
        return np.full(len(eps), mod["media_y"], dtype=np.float64)
    x = forzamiento(eps, v, mod["delta"], mod["v_mediana"])
    h = nucleo_h(mod["K"], mod["tau0"], mod["beta"], mod["omega_G"], mod["phi"])
    return mod["G0"] * convolucion_causal(x, h)


def ll_por_obs(mod: dict, y: np.ndarray, eps: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Log-verosimilitud gaussiana POR OBSERVACION.

    Se devuelve el vector, no la suma: el bootstrap por bloques necesita
    remuestrear observaciones contiguas, no un escalar ya agregado.
    """
    pred = predecir(mod, eps, v)
    s = max(mod["sigma"], 1e-12)
    r = np.asarray(y, np.float64) - pred
    return -0.5 * np.log(2.0 * np.pi * s * s) - 0.5 * (r / s) ** 2


# ===========================================================================
# 3. Particion con embargo y bootstrap por bloques
# ===========================================================================

def particionar(n: int, embargo: int, fracciones=(0.60, 0.20, 0.20)) -> dict:
    """Corte temporal 60/20/20 con banda de embargo entre conjuntos.

    El embargo se DESCARTA, no se reparte: las observaciones dentro de una banda
    no pertenecen a ninguno de los dos lados.
    """
    f_ent, f_val, _ = fracciones
    i1 = int(n * f_ent)
    i2 = int(n * (f_ent + f_val))
    ent = np.arange(0, i1)
    val = np.arange(i1 + embargo, i2)
    pru = np.arange(i2 + embargo, n)
    if val.size <= 0 or pru.size <= 0:
        raise ValueError("embargo %d demasiado grande para n=%d" % (embargo, n))
    return {"entrenamiento": ent, "validacion": val, "prueba": pru,
            "embargo": int(embargo)}


def verificar_embargo(part: dict) -> bool:
    """El §2.2 exige verificarlo por test, no por comentario."""
    e = part["embargo"]
    ok1 = part["validacion"][0] - part["entrenamiento"][-1] > e
    ok2 = part["prueba"][0] - part["validacion"][-1] > e
    return bool(ok1 and ok2)


def bootstrap_bloques(serie: np.ndarray, longitud: int, n_sorteos: int,
                      semilla: int = 0) -> dict:
    """Bootstrap por bloques MOVILES sobre una serie por observacion.

    Devuelve la media y el IC al 95 %, mas el numero de bloques efectivos, que
    el preregistro exige reportar explicitamente: por debajo de 15 los
    intervalos van marcados como subpotenciados.
    """
    x = np.asarray(serie, dtype=np.float64)
    n = x.size
    longitud = int(min(longitud, max(1, n // 2)))
    n_bloques = int(np.ceil(n / longitud))
    rng = np.random.default_rng(semilla)
    inicios_max = n - longitud
    if inicios_max <= 0:
        return {"media": float(np.mean(x)), "ic": (np.nan, np.nan),
                "n_bloques": 1, "subpotenciado": True}
    medias = np.empty(n_sorteos)
    for s in range(n_sorteos):
        ini = rng.integers(0, inicios_max + 1, size=n_bloques)
        muestra = np.concatenate([x[i:i + longitud] for i in ini])[:n]
        medias[s] = float(np.mean(muestra))
    lo, hi = np.percentile(medias, [2.5, 97.5])
    return {"media": float(np.mean(x)), "ic": (float(lo), float(hi)),
            "n_bloques": n_bloques, "subpotenciado": bool(n_bloques < 15)}


def margen_bic(k: int, n: int, n_eff: int) -> float:
    """k * ln(N_eff) / (2N) -- el BIC con tamano muestral EFECTIVO, por observacion.

    Declarado en el §10.2 del preregistro. La orden escribia `k/(2N)`, que es
    medio AIC, y mezclaba el N ingenuo de ticks con la incertidumbre de N_eff de
    bloques: las opciones plausibles diferian en cinco ordenes de magnitud.
    """
    return k * np.log(max(n_eff, 2)) / (2.0 * max(n, 1))


# ===========================================================================
# 4. Contraste de omega_G contra NULO SIMULADO (nunca chi2)
# ===========================================================================

def contraste_omega_G(y, eps, v, K, delta, n_sorteos: int = 40,
                      semilla: int = 0) -> dict:
    """Razon de verosimilitud M2 contra M2-osc, con la distribucion del
    estadistico generada por SIMULACION bajo el nulo omega_G = 0.

    ⚠ Prohibido comparar 2*dLL contra la tabla chi2: con observaciones
    dependientes esa distribucion no se cumple, y seria la TERCERA vez que este
    proyecto usa una asintotica donde la finita es otra.
    """
    m2 = ajustar(y, eps, v, K, delta, beta_libre=True)
    mo = ajustar(y, eps, v, K, delta, con_oscilacion=True)
    ll2 = float(np.sum(ll_por_obs(m2, y, eps, v)))
    llo = float(np.sum(ll_por_obs(mo, y, eps, v)))
    obs = 2.0 * (llo - ll2)

    rng = np.random.default_rng(semilla)
    x = forzamiento(eps, v, delta, m2["v_mediana"])
    h = nucleo_h(K, m2["tau0"], m2["beta"])
    base = m2["G0"] * convolucion_causal(x, h)
    nulo = np.empty(n_sorteos)
    for s in range(n_sorteos):
        y_s = base + rng.normal(0.0, m2["sigma"], size=len(y))
        a = ajustar(y_s, eps, v, K, delta, beta_libre=True)
        b = ajustar(y_s, eps, v, K, delta, con_oscilacion=True)
        nulo[s] = 2.0 * (float(np.sum(ll_por_obs(b, y_s, eps, v)))
                         - float(np.sum(ll_por_obs(a, y_s, eps, v))))
    p = float(np.mean(nulo >= obs))
    return {"estadistico": obs, "p_simulado": p, "omega_G": mo["omega_G"],
            "nulo_p95": float(np.percentile(nulo, 95)), "n_sorteos": n_sorteos}


# ===========================================================================
# 5. Generador sintetico -- la verdad conocida de los controles
# ===========================================================================

def generar(n: int, G0: float, tau0: float, beta: float, delta: float,
            sigma: float, K: int = 200, gamma_signos: float = 0.0,
            semilla: int = 0) -> dict:
    """Serie con propagador conocido. `gamma_signos` da memoria al flujo de
    ordenes, que es la propiedad que el nulo de signos iid NO reproduce (leccion
    de la sesion 2026-08-08 e)."""
    rng = np.random.default_rng(semilla)
    if gamma_signos > 0.0:
        # Signos con autocorrelacion: cadena de Markov de dos estados.
        p = 0.5 * (1.0 + gamma_signos)
        eps = np.empty(n)
        eps[0] = rng.choice([-1.0, 1.0])
        cambia = rng.random(n) > p
        for t in range(1, n):
            eps[t] = -eps[t - 1] if cambia[t] else eps[t - 1]
    else:
        eps = rng.choice([-1.0, 1.0], size=n).astype(float)
    v = np.exp(rng.normal(0.0, 1.0, size=n))       # colas pesadas, como el volumen
    x = forzamiento(eps, v, delta)
    h = nucleo_h(K, tau0, beta)
    y = G0 * convolucion_causal(x, h) + rng.normal(0.0, sigma, size=n)
    return {"y": y, "eps": eps, "v": v,
            "verdad": {"G0": G0, "tau0": tau0, "beta": beta, "delta": delta,
                       "sigma": sigma}}


# ===========================================================================
# 6. Autotest -- controles POSITIVOS
# ===========================================================================

def _autotest() -> int:
    fallos = 0

    def ok(nombre, cond, detalle=""):
        nonlocal fallos
        if not cond:
            fallos += 1
        print("  [%s] %s %s" % ("OK  " if cond else "FALLA", nombre, detalle))

    K = 120

    print("== 1. CONTROL POSITIVO: recupera (G0, tau0, beta) conocidos ==")
    d = generar(40000, G0=0.5, tau0=60.0, beta=0.35, delta=0.5, sigma=0.20,
                K=K, semilla=1)
    m = ajustar(d["y"], d["eps"], d["v"], K, delta=0.5)
    dec_verdad = float(nucleo_G(K, 60.0, 0.35))
    print("     verdad  G0=0.500  tau0=60.0  beta=0.350  ->  G(K)/G(0)=%.4f" % dec_verdad)
    print("     medido  G0=%.3f  tau0=%.1f  beta=%.3f  ->  G(K)/G(0)=%.4f"
          % (m["G0"], m["tau0"], m["beta"], decaimiento(m)))
    ok("G0 al 10 %", abs(m["G0"] - 0.5) / 0.5 < 0.10)
    ok("decaimiento al 10 %", abs(decaimiento(m) - dec_verdad) / dec_verdad < 0.10,
       "(beta y tau0 por separado NO estan identificados)")

    print("== 2. CONTROL POSITIVO: impacto PERMANENTE verdadero (beta = 0) ==")
    d0 = generar(40000, G0=0.5, tau0=60.0, beta=0.0, delta=0.5, sigma=0.20,
                 K=K, semilla=2)
    m0 = ajustar(d0["y"], d0["eps"], d0["v"], K, delta=0.5)
    print("     beta medido = %.4f (verdad 0)  <- NO interpretable" % m0["beta"])
    print("     G(K)/G(0)   = %.4f (verdad 1.0000)  <- lo identificado"
          % decaimiento(m0))
    ok("recupera impacto permanente", abs(decaimiento(m0) - 1.0) < 0.05)

    print("== 3. M2 bate a M1 FUERA DE MUESTRA solo si el impacto DECAE ==")
    # Series largas a proposito: con bloques de 5*embargo = 9750 ticks hacen
    # falta >= 15 bloques efectivos, y el preregistro obliga a marcar como
    # subpotenciado cualquier intervalo que no llegue.
    dL = generar(250000, G0=0.5, tau0=60.0, beta=0.35, delta=0.5, sigma=0.20,
                 K=K, semilla=11)
    dL0 = generar(250000, G0=0.5, tau0=60.0, beta=0.0, delta=0.5, sigma=0.20,
                  K=K, semilla=12)
    for etiq, dd, espera in (("decae", dL, True), ("permanente", dL0, False)):
        part = particionar(len(dd["y"]), EMBARGO_PISO)
        ent, pru = part["entrenamiento"], part["prueba"]
        a1 = ajustar(dd["y"][ent], dd["eps"][ent], dd["v"][ent], K, 0.5,
                     beta_libre=False)
        a2 = ajustar(dd["y"][ent], dd["eps"][ent], dd["v"][ent], K, 0.5)
        l1 = ll_por_obs(a1, dd["y"][pru], dd["eps"][pru], dd["v"][pru])
        l2 = ll_por_obs(a2, dd["y"][pru], dd["eps"][pru], dd["v"][pru])
        b = bootstrap_bloques(l2 - l1, 5 * EMBARGO_PISO, 200, semilla=3)
        marg = margen_bic(3, pru.size, b["n_bloques"])
        gana = b["ic"][0] > marg
        print("     %s | dLL/N=%+.5f IC=[%+.5f, %+.5f] margen=%.2e bloques=%d"
              % (etiq, b["media"], b["ic"][0], b["ic"][1], marg, b["n_bloques"]))
        ok("M2 %s a M1 con %s" % ("gana" if espera else "NO gana", etiq), gana == espera)

    print("== 4. CONTROL NEGATIVO: sin forzamiento, M0 no pierde ==")
    dn = generar(40000, G0=0.0, tau0=60.0, beta=0.35, delta=0.5, sigma=0.20,
                 K=K, semilla=4)
    part = particionar(len(dn["y"]), EMBARGO_PISO)
    ent, pru = part["entrenamiento"], part["prueba"]
    a0 = ajustar_M0(dn["y"][ent])
    a2 = ajustar(dn["y"][ent], dn["eps"][ent], dn["v"][ent], K, 0.5)
    l0 = ll_por_obs(a0, dn["y"][pru], dn["eps"][pru], dn["v"][pru])
    l2 = ll_por_obs(a2, dn["y"][pru], dn["eps"][pru], dn["v"][pru])
    b = bootstrap_bloques(l2 - l0, 5 * EMBARGO_PISO, 200, semilla=5)
    marg = margen_bic(4, pru.size, b["n_bloques"])
    print("     dLL/N=%+.6f IC=[%+.6f, %+.6f] margen=%.2e"
          % (b["media"], b["ic"][0], b["ic"][1], marg))
    ok("M2 NO bate a M0 sin forzamiento", not (b["ic"][0] > marg))

    print("== 5. delta se recupera de la rejilla en VALIDACION ==")
    dd = generar(40000, G0=0.5, tau0=60.0, beta=0.35, delta=1.0, sigma=0.20,
                 K=K, semilla=6)
    part = particionar(len(dd["y"]), EMBARGO_PISO)
    ent, val = part["entrenamiento"], part["validacion"]
    mejor, mejor_ll = None, -np.inf
    for de in REJILLA_DELTA:
        a = ajustar(dd["y"][ent], dd["eps"][ent], dd["v"][ent], K, de)
        llv = float(np.mean(ll_por_obs(a, dd["y"][val], dd["eps"][val], dd["v"][val])))
        print("     delta=%.2f -> LL/N validacion = %+.5f" % (de, llv))
        if llv > mejor_ll:
            mejor, mejor_ll = de, llv
    ok("elige delta = 1.0 (verdad)", mejor == 1.0, "elegida %.2f" % mejor)

    print("== 6. El embargo se respeta y se verifica por TEST ==")
    p = particionar(100000, EMBARGO_PISO)
    print("     ent=%d val=%d pru=%d | huecos %d y %d ticks"
          % (p["entrenamiento"].size, p["validacion"].size, p["prueba"].size,
             p["validacion"][0] - p["entrenamiento"][-1],
             p["prueba"][0] - p["validacion"][-1]))
    ok("embargo verificado", verificar_embargo(p))
    # La compuerta del preregistro sale de encadenar los propios requisitos:
    # >= 15 bloques de longitud >= 5*embargo, sobre una prueba que es el 20 %.
    n_compuerta = ticks_de_compuerta()
    pc = particionar(n_compuerta, EMBARGO_PISO)
    bloques_en_prueba = pc["prueba"].size // (5 * EMBARGO_PISO)
    print("     compuerta = %d ticks (%.0f*embargo): prueba=%d -> %d bloques de %d"
          % (n_compuerta, n_compuerta / EMBARGO_PISO, pc["prueba"].size,
             bloques_en_prueba, 5 * EMBARGO_PISO))
    ok("la compuerta derivada da >= %d bloques" % BLOQUES_MINIMOS,
       bloques_en_prueba >= BLOQUES_MINIMOS)
    menos = particionar(375 * EMBARGO_PISO, EMBARGO_PISO)
    print("     contraprueba: 375*embargo da %d bloques -> insuficiente"
          % (menos["prueba"].size // (5 * EMBARGO_PISO)))
    ok("375*embargo se queda corto (por eso no es el umbral)",
       menos["prueba"].size // (5 * EMBARGO_PISO) < BLOQUES_MINIMOS)

    print("== 7. La convolucion es CAUSAL (no mete el futuro en el presente) ==")
    imp = np.zeros(50); imp[10] = 1.0
    h = nucleo_h(20, 5.0, 0.5)
    z = convolucion_causal(imp, h)
    ok("cero antes del impulso", np.allclose(z[:10], 0.0),
       "max|z[:10]| = %.2e" % float(np.max(np.abs(z[:10]))))

    print("== 8. omega_G contra NULO SIMULADO, sobre datos SIN oscilacion ==")
    dcorto = generar(6000, G0=0.5, tau0=60.0, beta=0.35, delta=0.5, sigma=0.20,
                     K=60, semilla=7)
    c = contraste_omega_G(dcorto["y"], dcorto["eps"], dcorto["v"], 60, 0.5,
                          n_sorteos=25, semilla=8)
    print("     2dLL observado = %.3f | p95 del nulo = %.3f | p simulado = %.3f"
          % (c["estadistico"], c["nulo_p95"], c["p_simulado"]))
    ok("no rechaza omega_G = 0 sobre datos sin oscilacion", c["p_simulado"] > 0.05)
    print("     (el estadistico NO se compara contra chi2; ese es el punto)")

    print("== 9. El margen de la Sec. 10.2 es el BIC con N_eff, no k/(2N) ==")
    for n, ne in ((150000, 15), (1000000, 138)):
        print("     N=%7d N_eff=%3d | k/(2N)=%.2e | declarado=%.2e"
              % (n, ne, 3 / (2 * n), margen_bic(3, n, ne)))
    ok("el margen declarado supera a k/(2N)",
       margen_bic(3, 150000, 15) > 3 / (2 * 150000))

    print("")
    print("RESULTADO: %d fallo(s)" % fallos)
    return fallos


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--compuerta", action="store_true")
    ap.add_argument("--dir", default="telemetria/captura_v32")
    a = ap.parse_args()
    if a.autotest:
        return 1 if _autotest() else 0
    if a.compuerta:
        import captura_larga
        return captura_larga.resumen(a.dir)
    print("uso: python migracion_v32.py --autotest | --compuerta")
    return 0


if __name__ == "__main__":
    sys.exit(main())
