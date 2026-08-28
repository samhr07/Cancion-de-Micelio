# -*- coding: utf-8 -*-
"""
actividad.py -- las tres preguntas del operador sobre `nu` (2026-08-27).

    python actividad.py --autotest      controles con verdad conocida
    python actividad.py --etapa=ticks   1. el horizonte en TIEMPO DE TRANSACCIONES
    python actividad.py --etapa=ucurva  2. la curva en U: ?volatilidad o actividad?
    python actividad.py --etapa=espectro 3. nulo espectral sobre la serie de nu

CONTEXTO: si el lazo actividad -> volatilidad -> retirada de liquidez ->
actividad es un OSCILADOR ENDOGENO. Tres partes, en orden de prioridad.

=========================================================================
1. EL HORIZONTE EN TIEMPO DE TRANSACCIONES
=========================================================================

La curva requerida vive hoy en segundos. El coste es fijo en pb y la volatilidad
crece con el horizonte, asi que `R2_req = (lastre / (1.755*sigma(H)))^2`. Si el
reloj natural es el de transacciones -- que es lo que dice Clark (1973) y lo que
H2 confirmo -- la curva hay que reexpresarla con `n` = numero de transacciones.

Linea base de Clark: `sigma^2(n) ~ n`, o sea `H_p_ticks = 0.5` exacto. La
correccion es la `H_p` medida en reloj de ticks.

⚠ **DECLARADO ANTES DE MIRAR** (lo pide el propio encargo): condicionar al decil
de `nu` sube `sigma` y por tanto BAJA `R2_req`. Pero **si `R2_medido` escala con
`nu` en la misma proporcion, la RAZON no mejora y el condicionamiento no compra
nada**. Lo que decide es la razon `R2_medido / R2_req` por decil, no `R2_req`
sola. Y la perdida de potencia por condicionar (ventanas por decil) va **en la
tabla**, no en una nota.

=========================================================================
2. LA CURVA EN U: ?VOLATILIDAD O ACTIVIDAD?
=========================================================================

`TRASPASO_SESION_2026-08-10.md` §4.3 dejo abierta la pregunta. Se contrasta
`sigma_t` contra `nu_t` como predictores de `sigma_{t+1}` fuera de muestra. Si
`nu` domina, la pregunta queda resuelta a favor de ACTIVIDAD y se retira.

⚠ Marco, no hallazgo: Clark (1973) y Tauchen & Pitts (1983).

=========================================================================
3. NULO ESPECTRAL SOBRE LA SERIE DE ACTIVIDAD
=========================================================================

El multitaper de la v2.2 (x1.78 sobre ruido rojo) se midio sobre el PRECIO y
sobre dato contaminado por transacciones con `p = 0`. La hipotesis de ciclo
endogeno vive en la serie de `nu` y **nunca se ha contrastado ahi**.

Tres correcciones sobre aquel procedimiento:
  a) sobre `nu`, no sobre el precio, y con dato limpio;
  b) nulo de ruido rojo **SIMULADO**, no tabla chi2 asintotica. La tabla ignora
     que `a` del AR(1) se estima del propio dato y no corrige comparaciones
     multiples; se usan sustitutos y ademas el estadistico del MAXIMO;
  c) **el ciclo diurno se retira primero.** Su periodo es 24 h, exogeno y ya
     medido: dejarlo dentro garantiza un pico que no responde a la pregunta.
     La hipotesis endogena vive en minutos.

⚠ **EXPECTATIVA DECLARADA ANTES DE CORRER: no habra pico.** El mecanismo
propuesto (1) carece de punto de consigna fijo -- cada participante calibra su
umbral con la volatilidad reciente --, (2) tiene retardos dispersos sobre seis
ordenes de magnitud y (3) es asimetrico: retirada en milisegundos, reposicion en
minutos. Un lazo con esas tres propiedades produce **memoria larga, no
periodicidad**. Si aparece pico es hallazgo; si no, cierra una pregunta abierta
desde la v2.2 y se reporta igual.
"""

from __future__ import annotations

import argparse
import gc
import os

import numpy as np

import estacionalidad as E
import horizonte as H
import identidad as I

log, titulo = H.log, H.titulo

BIN_DECIL_S = 900.0                       # 15 min para asignar decil de nu
NS_TICKS = [32, 128, 512, 2048, 8192]
ESQUINAS = I.ESQUINAS                     # (nombre, H_p, lastre_pb)
N_SUR = 200


# ===========================================================================
# 1. Horizonte en tiempo de transacciones
# ===========================================================================

def _bins_nu(f, bin_s=BIN_DECIL_S):
    """Indices de tick por casilla y `nu` de cada casilla."""
    t = f["t"]
    b0 = np.ceil(t[0] / bin_s) * bin_s
    idx = np.floor((t - b0) / bin_s).astype(np.int64)
    ok = idx >= 0
    nb = int(idx[ok].max()) + 1
    cnt = np.bincount(idx[ok], minlength=nb).astype(float)
    return idx, cnt / bin_s, nb


def etapa_ticks(args) -> int:
    titulo("1. EL HORIZONTE EN TIEMPO DE TRANSACCIONES, POR DECIL DE nu")
    log("")
    log("  Linea base de Clark: sigma^2(n) ~ n  =>  H_p_ticks = 0.5 exacto.")
    log("  DECLARADO ANTES DE MIRAR: condicionar sube sigma y baja R2_req, pero si")
    log("  R2_medido escala igual con nu la RAZON no mejora. Decide la razon.")

    # --- deciles de nu agrupando todos los fragmentos
    todas = []
    for f in I.fragmentos():
        _, nu, _ = _bins_nu(f)
        todas.append(nu[nu > 0])
    nu_all = np.concatenate(todas)
    cortes = np.percentile(nu_all, [10, 30, 50, 70, 90])
    log("")
    log("  cortes de nu [tx/s]: %s" % np.array2string(cortes, precision=2))
    grupos = ["D1 (mas bajo)", "D2", "D3", "D4", "D5", "D6 (mas alto)"]

    acc = {g: {n: {"r": [], "e": []} for n in NS_TICKS} for g in grupos}
    for f in I.fragmentos():
        idx, nu, nb = _bins_nu(f)
        gi = np.digitize(nu, cortes)                 # 0..5
        mid, eps = f["mid"], f["eps"].astype(float)
        # [!] BORDES POR BUSQUEDA BINARIA, no escaneando. `idx` es monotono
        # porque `t` lo es, asi que `searchsorted` da el rango de cada casilla en
        # log(n). La version con `flatnonzero(idx == b)` dentro del bucle era
        # O(n_casillas * n_ticks): sobre 20 M de ticks, miles de millones de
        # operaciones por horizonte.
        lo = np.searchsorted(idx, np.arange(nb), side="left")
        hi = np.searchsorted(idx, np.arange(nb), side="right")
        for n in NS_TICKS:
            for b in range(nb):
                if hi[b] - lo[b] < 2 * n:
                    continue
                k = np.arange(lo[b], hi[b] - n, n)
                if k.size < 2:
                    continue
                r = np.log(mid[k + n] / mid[k]) * 1e4
                ok = np.isfinite(r)
                acc[grupos[gi[b]]][n]["r"].append(r[ok])
                acc[grupos[gi[b]]][n]["e"].append(eps[k][ok])
        del idx, nu, gi, mid, eps
        gc.collect()

    log("")
    log("  %-14s %7s | %9s %9s | %10s %10s | %9s"
        % ("decil de nu", "n_ticks", "n_vent", "sigma[pb]", "R2_req fav", "R2_req adv", "R2_medido"))
    log("  " + "-" * 88)
    razones = {}
    for g in grupos:
        for n in NS_TICKS:
            r = np.concatenate(acc[g][n]["r"]) if acc[g][n]["r"] else np.array([])
            e = np.concatenate(acc[g][n]["e"]) if acc[g][n]["e"] else np.array([])
            if r.size < 60:
                continue
            sg = float(np.std(r))
            reqs = [(la / (H.FACTOR_DECIL * sg)) ** 2 for _, _, la in ESQUINAS]
            Rc = float(np.mean(e * r) - np.mean(e) * np.mean(r))
            r2m = Rc * Rc / (np.var(e) * np.var(r)) if np.var(e) > 0 else float("nan")
            razones.setdefault(n, []).append((g, r2m / max(reqs[0], 1e-30), r2m, reqs[0], r.size))
            log("  %-14s %7d | %9d %9.3f | %9.2f %% %9.2f %% | %9.6f"
                % (g, n, r.size, sg, 100 * reqs[0], 100 * reqs[1], r2m))
    log("")
    log("  --- LO QUE DECIDE: la razon R2_medido / R2_req(favorable) por decil ---")
    for n in NS_TICKS:
        if n not in razones:
            continue
        fil = razones[n]
        log("    n = %5d ticks : %s" % (n, "  ".join("%s %.2e" % (g.split()[0], z)
                                                     for g, z, _, _, _ in fil)))
        zs = [z for _, z, _, _, _ in fil]
        log("                      mejor %.2e en %s   |  razon mejor/peor = %.1fx"
            % (max(zs), fil[int(np.argmax(zs))][0].split()[0],
               max(zs) / max(min(zs), 1e-300)))
    log("")
    log("  [!] La perdida de potencia por condicionar esta en la columna n_vent, y")
    log("      es parte del resultado: cada decil se queda con ~1/6 de las ventanas.")
    return 0


# ===========================================================================
# 2. La curva en U
# ===========================================================================

def _bloques(n, k=4):
    """Pliegues de bloques CONTIGUOS: no hay fuga entre casillas vecinas."""
    b = np.array_split(np.arange(n), k)
    return b


def _r2_oos_bloques(X, y, k=4):
    n = y.size
    num = den = 0.0
    for te in _bloques(n, k):
        tr = np.setdiff1d(np.arange(n), te)
        if tr.size < X.shape[1] + 5 or te.size < 5:
            continue
        c, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None)
        p = X[te] @ c
        num += float(np.sum((y[te] - p) ** 2))
        den += float(np.sum((y[te] - y[tr].mean()) ** 2))
    return 1.0 - num / den if den > 0 else float("nan")


def _serie_frag(f, bin_s=300.0):
    """`log nu` y `log sigma` por casilla, de UN fragmento cacheado."""
    t, mid = f["t"], f["mid"]
    b0 = np.ceil(t[0] / bin_s) * bin_s
    b1 = np.floor(t[-1] / bin_s) * bin_s
    nb = int((b1 - b0) / bin_s)
    if nb < 12:
        return None
    idx = np.floor((t - b0) / bin_s).astype(np.int64)
    ok = (idx >= 0) & (idx < nb)
    cnt = np.bincount(idx[ok], minlength=nb).astype(float)
    g = np.arange(b0, b1 + 1e-9, 30.0)
    j = np.clip(np.searchsorted(t, g, side="left"), 0, t.size - 1)
    r = np.diff(np.log(mid[j])) * 1e4
    gi = np.floor((g[:-1] - b0) / bin_s).astype(np.int64)
    s2 = np.bincount(gi, weights=r ** 2, minlength=nb)
    nn = np.bincount(gi, minlength=nb)
    val = (nn >= int(bin_s / 30.0) - 1) & (cnt > 0) & (s2 > 0)
    return {"log_nu": np.log(cnt[val] / bin_s),
            "log_sigma": 0.5 * np.log(s2[val] / nn[val]), "val": val}


def etapa_ucurva(args) -> int:
    titulo("2. LA CURVA EN U: sigma_t contra nu_t prediciendo sigma_{t+1}")
    log("")
    log("  Marco: Clark (1973) y Tauchen & Pitts (1983). NO es hallazgo.")
    log("")
    log("  %-16s %6s | %9s %9s %9s %9s | %9s"
        % ("fragmento", "n", "constante", "sigma(t)", "nu(t)", "ambos", "nulo q95"))
    log("  " + "-" * 82)
    gana_nu = gana_si = 0
    rng = np.random.default_rng(9)
    for f in I.fragmentos():
        d = _serie_frag(f)
        if d is None:
            continue
        ls, ln = d["log_sigma"], d["log_nu"]
        if ls.size < 40:
            continue
        y = ls[1:]
        uno = np.ones(y.size)
        Xs = np.column_stack([uno, ls[:-1]])
        Xn = np.column_stack([uno, ln[:-1]])
        Xb = np.column_stack([uno, ls[:-1], ln[:-1]])
        r_c = _r2_oos_bloques(uno.reshape(-1, 1), y)
        r_s = _r2_oos_bloques(Xs, y)
        r_n = _r2_oos_bloques(Xn, y)
        r_b = _r2_oos_bloques(Xb, y)
        nul = [_r2_oos_bloques(Xn, np.roll(y, int(rng.integers(5, y.size - 5))))
               for _ in range(N_SUR)]
        q95 = float(np.percentile(np.abs(np.array(nul)), 95))
        log("  %-16s %6d | %+9.4f %+9.4f %+9.4f %+9.4f | %9.4f"
            % (f["nombre"], y.size, r_c, r_s, r_n, r_b, q95))
        if r_n > r_s:
            gana_nu += 1
        else:
            gana_si += 1
        del d
        gc.collect()
    log("")
    log("  nu gana en %d fragmentos, sigma en %d" % (gana_nu, gana_si))
    if gana_nu > gana_si:
        log("  *** RESUELTA A FAVOR DE ACTIVIDAD. La pregunta 'curva en U: volatilidad")
        log("      o actividad' del TRASPASO §4.3 se retira de la lista. ***")
    elif gana_si > gana_nu:
        log("  *** sigma domina: la pregunta NO se resuelve a favor de actividad. ***")
    else:
        log("  *** empate: no se retira nada. ***")
    return 0


# ===========================================================================
# 3. Nulo espectral sobre nu
# ===========================================================================

def _multitaper(x, dt, NW=4.0, K=7):
    from scipy.signal.windows import dpss
    n = x.size
    tap = dpss(n, NW, K)
    frec = np.fft.rfftfreq(n, d=dt)
    xc = x - x.mean()
    S = np.zeros(frec.size)
    for t in tap:
        S += np.abs(np.fft.rfft(xc * t)) ** 2
    return frec, S / K * dt


def _ar1_sur(x, rng):
    a = float(np.corrcoef(x[:-1], x[1:])[0, 1])
    a = min(max(a, 0.0), 0.995)
    s = float(np.std(x)) * np.sqrt(max(1 - a * a, 1e-9))
    y = np.empty(x.size)
    y[0] = rng.normal(0, float(np.std(x)))
    e = rng.normal(0, s, x.size)
    for i in range(1, x.size):
        y[i] = a * y[i - 1] + e[i]
    return y, a


def etapa_espectro(args) -> int:
    titulo("3. NULO ESPECTRAL SOBRE LA SERIE DE ACTIVIDAD")
    log("")
    log("  EXPECTATIVA DECLARADA ANTES DE CORRER: no habra pico. El mecanismo no")
    log("  tiene consigna fija, sus retardos se dispersan sobre seis ordenes de")
    log("  magnitud y es asimetrico. Eso produce MEMORIA LARGA, no periodicidad.")
    d = E.serie_binada(bin_s=300.0)
    y = d["log_nu"]
    # (c) el ciclo diurno se retira primero: su periodo es 24 h y es exogeno.
    X = E.diseno(d, E.K_ARM, True)
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c
    log("")
    log("  serie: %d casillas de 300 s.  varianza explicada por el ciclo diurno"
        " dentro de muestra: %.1f %%" % (y.size, 100 * (1 - np.var(res) / np.var(y))))
    # tramo contiguo mas largo del residuo, que es lo que el multitaper exige
    t = d["t"]
    corte = np.flatnonzero(np.diff(t) > 1.5 * 300.0)
    ini = np.r_[0, corte + 1]
    fin = np.r_[corte, t.size - 1]
    k = int(np.argmax(fin - ini))
    x = res[ini[k]:fin[k] + 1]
    log("  tramo contiguo mas largo del residuo: %d casillas = %.2f h"
        % (x.size, x.size * 300.0 / 3600.0))
    if x.size < 128:
        log("  *** insuficiente para multitaper ***")
        return 0
    frec, S = _multitaper(x, 300.0)
    rng = np.random.default_rng(21)
    sur = np.empty((N_SUR, frec.size))
    a_ = None
    for i in range(N_SUR):
        z, a_ = _ar1_sur(x, rng)
        sur[i] = _multitaper(z, 300.0)[1]
    q95 = np.percentile(sur, 95, axis=0)
    v = frec > 0
    exceso = S[v] / np.maximum(q95[v], 1e-300)
    fv = frec[v]
    # estadistico del MAXIMO: corrige comparaciones multiples
    max_sur = np.array([np.max(sur[i][v] / np.maximum(q95[v], 1e-300))
                        for i in range(N_SUR)])
    obs_max = float(np.max(exceso))
    p_glob = float(np.mean(max_sur >= obs_max))
    im = int(np.argmax(exceso))
    per = 1.0 / fv[im]
    log("")
    log("  AR(1) ajustado al residuo: a = %.4f   %d sustitutos SIMULADOS" % (a_, N_SUR))
    log("  frecuencias sobre el q95 puntual : %d de %d (%.1f %%, por azar 5 %%)"
        % (int(np.sum(exceso > 1)), exceso.size, 100 * np.mean(exceso > 1)))
    log("  pico: periodo %.0f s = %.2f h   exceso x%.2f   p GLOBAL (max) = %.4f"
        % (per, per / 3600.0, obs_max, p_glob))
    # guarda de banda de la v2.1 §2
    rango = x.size * 300.0
    log("")
    log("  --- guarda de banda de la v2.1 §2 ---")
    log("  rango ajustado = %.0f s;  periodo del pico / rango = %.3f" % (rango, per / rango))
    if per > rango / 3.0:
        log("  *** el periodo del pico excede un tercio del rango: es una INCLINACION,")
        log("      no un ciclo. No se lee como periodicidad. ***")
    elif p_glob < 0.05:
        log("  *** PICO SIGNIFICATIVO Y DENTRO DE BANDA: es hallazgo. ***")
    else:
        log("  *** NO HAY PICO. La expectativa declarada se cumple y la pregunta de")
        log("      ciclo endogeno abierta desde la v2.2 queda CERRADA sobre nu. ***")
    return 0


# ===========================================================================
# Controles
# ===========================================================================

def _autotest() -> int:
    titulo("CONTROLES DE actividad.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-52s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    rng = np.random.default_rng(3)
    n = 2048
    dt = 300.0
    # --- positivo: seno enterrado en AR(1) debe salir
    per = 20 * dt
    x = np.sin(2 * np.pi * np.arange(n) * dt / per)
    ru, _ = _ar1_sur(rng.normal(0, 1, n), rng)
    s = 1.2 * x + 2.0 * ru / np.std(ru)
    fr, S = _multitaper(s, dt)
    sur = np.array([_multitaper(_ar1_sur(s, rng)[0], dt)[1] for _ in range(60)])
    q = np.percentile(sur, 95, axis=0)
    v = fr > 0
    im = int(np.argmax(S[v] / np.maximum(q[v], 1e-300)))
    p_det = 1.0 / fr[v][im]
    chk(abs(p_det - per) / per < 0.15, "el multitaper encuentra un ciclo conocido",
        "%.0f s contra %.0f" % (p_det, per))
    # --- negativo: sobre AR(1) puro la TASA de rechazo debe rondar el 5 %
    #
    # [!] LA PRIMERA VERSION DE ESTE CONTROL ESTABA MAL DISENADA. Exigia que UN
    # solo sorteo diera `p > 0.05`, y bajo el nulo `p` es uniforme en [0,1]: ese
    # test falla el 5 % de las veces POR CONSTRUCCION, y fallo. Lo que hay que
    # comprobar es la CALIBRACION -- que la tasa de rechazo sobre muchas series
    # ronde el nominal --, no el resultado de una moneda. Es la segunda vez hoy
    # que una asercion mia pide de un estimador mas certeza de la que declara.
    ps = []
    for _ in range(20):
        z, _ = _ar1_sur(rng.normal(0, 1, 512), rng)
        fr2, S2 = _multitaper(z, dt)
        sur2 = np.array([_multitaper(_ar1_sur(z, rng)[0], dt)[1] for _ in range(40)])
        q2 = np.percentile(sur2, 95, axis=0)
        v2 = fr2 > 0
        exc2 = S2[v2] / np.maximum(q2[v2], 1e-300)
        mx = np.array([np.max(sur2[i][v2] / np.maximum(q2[v2], 1e-300)) for i in range(40)])
        ps.append(float(np.mean(mx >= np.max(exc2))))
    ps = np.array(ps)
    tasa = float(np.mean(ps < 0.05))
    chk(tasa <= 0.25, "sobre AR(1) puro la tasa de rechazo global no se dispara",
        "%.0f %% de 20 series (nominal 5 %%)" % (100 * tasa))
    chk(0.2 < float(np.mean(ps)) < 0.8, "y el p medio no esta pegado a un extremo",
        "p medio %.3f" % float(np.mean(ps)))

    # --- OOS por bloques: sin senal da ~0, con senal la recupera
    yy = rng.normal(0, 1, 600)
    Xr = np.column_stack([np.ones(600), rng.normal(0, 1, 600)])
    chk(abs(_r2_oos_bloques(Xr, yy)) < 0.05, "sin senal el R2 por bloques es ~0",
        "%.4f" % _r2_oos_bloques(Xr, yy))
    xs = rng.normal(0, 1, 600)
    chk(_r2_oos_bloques(np.column_stack([np.ones(600), xs]), 2 * xs + 0.3 * rng.normal(0, 1, 600)) > 0.9,
        "con senal fuerte la recupera")
    log("")
    log("RESULTADO: %d fallo(s)" % fallos)
    return 1 if fallos else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", default="")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    return {"ticks": etapa_ticks, "ucurva": etapa_ucurva,
            "espectro": etapa_espectro}.get(
        a.etapa, lambda _: (log("usa --autotest o --etapa=ticks|ucurva|espectro"), 2)[1])(a)


if __name__ == "__main__":
    raise SystemExit(main())
