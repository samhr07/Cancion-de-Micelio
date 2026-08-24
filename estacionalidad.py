# -*- coding: utf-8 -*-
"""
estacionalidad.py -- el ciclo diurno de `nu` y `sigma` como MODELO, no como media.

    python estacionalidad.py --autotest     controles con verdad conocida
    python estacionalidad.py --etapa=perfil descriptivo por hora y dia
    python estacionalidad.py --etapa=modelo la escalera de modelos, fuera de muestra

POR QUE EXISTE. El Sec.1 de la v4.1 se midio agrupando las 24 horas y el
veredicto quedo retractado el 2026-08-23: `R2_req` va como `(lastre/sigma)^2` y
`sigma_60` recorre de 3.29 a 9.54 pb segun la hora UTC, o sea que **el requisito
varia 8.4x**. Ademas `nu` entraba en el predictor como ESCALAR mientras recorre
0.55x a 2.17x de su media. Una media no describe esto: hace falta una funcion.

[!] ESTO NO RESUCITA `omega_m`, Y LA DISTINCION ES EL PUNTO ENTERO.
    `omega_m` murio porque era una frecuencia ENDOGENA sin ancla y sin nulo: la
    EMD devolvia un "ciclo" de 118 s sobre un paseo aleatorio puro (v2.2), y el
    periodo escalaba con la ventana de analisis en vez de con la senal.
    Aqui el periodo NO SE ESTIMA: es 24 h y 168 h, conocido a priori, con causa
    exogena verificable (sesiones de mercado y horario laboral humano). Lo que
    se estima es AMPLITUD y FASE de una frecuencia dada, que es un problema
    lineal, con hipotesis nula propia (barajar etiquetas de dia completo) y con
    validacion fuera de muestra por dias enteros. Es la diferencia entre buscar
    una frecuencia y medir la respuesta a una que ya sabes cual es.

EL MODELO, en la forma que pidio el operador:

    log X(t) = a_0
             + SUM_k [ a_k cos(2*pi*k*h/24) + b_k sin(2*pi*k*h/24) ]     <- topes y valles
             + finde * ( c_0 + SUM_k [ c_k cos(...) + d_k sin(...) ] )   <- amplitud Y FASE
             + SUM_s  A_s * exp( -delta_s(t) / tau )                     <- pulsos de sesion
                                                                            con decaimiento

`delta_s(t)` es el tiempo transcurrido desde la apertura `s` mas reciente, con
envoltura circular de 24 h. Un par (cos, sin) con coeficientes libres ES una
sinusoide de amplitud y fase libres, asi que dejar que el bloque de fin de
semana tenga los suyos es exactamente el "desfase de fin de semana": no hace
falta parametrizar la fase aparte, y asi el ajuste sigue siendo LINEAL y no
tiene minimos locales -- que es lo que hundio el estadistico de omega_G en la
v3.2. El unico parametro no lineal es `tau`, y se barre en rejilla y se reporta
la curva entera en vez de elegir un punto.
"""

from __future__ import annotations

import argparse
import datetime as DT
import os

import numpy as np

import curvas_estacional as C
import horizonte as H

log, titulo = H.log, H.titulo

BIN_S = 300.0          # 5 min: 288 casillas por dia
SUB_S = 30.0           # submuestreo para la volatilidad realizada dentro del bin
K_ARM = 4              # armonicos de 24 h. El 3o y el 4o dan la ASIMETRIA
                       # (subida rapida en la apertura, caida lenta de noche)
                       # que una sinusoide sola no puede representar.
# Aperturas en hora UTC. Son horarios de calendario, no cantidades ajustadas.
SESIONES = {"Tokio": 0.0, "Londres": 7.0, "NuevaYork": 13.5, "CierreUS": 20.0}
TAUS = [0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 9.0]      # horas
EPOCA = 0.0


# ===========================================================================
# Serie binada
# ===========================================================================

def serie_binada(bin_s: float = BIN_S) -> dict:
    """`log nu` y `log RV` por casilla de `bin_s`, etiquetadas en UTC."""
    tb, lnu, lrv = [], [], []
    for k in C._tramos_cache():
        m = dict(np.load(C.CACHE % k))
        t = m["t"].astype(float)
        mid = 0.5 * (m["bid"] + m["ask"])
        del m
        b0 = np.ceil(t[0] / bin_s) * bin_s
        b1 = np.floor(t[-1] / bin_s) * bin_s
        if b1 - b0 < 4 * bin_s:
            continue
        # --- nu por casilla
        idx = np.floor((t - b0) / bin_s).astype(np.int64)
        ok = (idx >= 0) & (idx < int((b1 - b0) / bin_s))
        cnt = np.bincount(idx[ok], minlength=int((b1 - b0) / bin_s)).astype(float)
        # --- volatilidad realizada por casilla, sobre rejilla de SUB_S
        g = np.arange(b0, b1 + 1e-9, SUB_S)
        j = np.clip(np.searchsorted(t, g, side="left"), 0, t.size - 1)
        r = np.diff(np.log(mid[j])) * 1e4                    # pb
        gi = np.floor((g[:-1] - b0) / bin_s).astype(np.int64)
        nb = int((b1 - b0) / bin_s)
        s2 = np.bincount(gi, weights=r ** 2, minlength=nb)
        nn = np.bincount(gi, minlength=nb)
        esperado = int(round(bin_s / SUB_S))
        val = (nn >= esperado - 1) & (cnt > 0) & (s2 > 0)
        centros = b0 + (np.arange(nb) + 0.5) * bin_s
        tb.append(centros[val])
        lnu.append(np.log(cnt[val] / bin_s))
        lrv.append(0.5 * np.log(s2[val] / nn[val]))          # log sigma por sub-paso
    t = np.concatenate(tb)
    o = np.argsort(t)
    t = t[o]
    d = {"t": t, "log_nu": np.concatenate(lnu)[o],
         "log_sigma": np.concatenate(lrv)[o]}
    u = [DT.datetime.utcfromtimestamp(x) for x in t]
    d["hora"] = np.array([x.hour + x.minute / 60.0 for x in u])
    d["dow"] = np.array([x.weekday() for x in u])
    d["dia"] = np.array([x.toordinal() for x in u])
    d["finde"] = (d["dow"] >= 5).astype(float)
    return d


# ===========================================================================
# Diseno
# ===========================================================================

def _armonicos(hora, K):
    col = []
    for k in range(1, K + 1):
        w = 2.0 * np.pi * k * hora / 24.0
        col += [np.cos(w), np.sin(w)]
    return col


def _pulsos(hora, tau_h):
    """`exp(-delta/tau)` con `delta` = horas desde la apertura mas reciente."""
    col = []
    for _, h0 in sorted(SESIONES.items()):
        d = np.mod(hora - h0, 24.0)
        col.append(np.exp(-d / tau_h))
    return col


def diseno(d: dict, K=K_ARM, finde=True, pulsos=None) -> np.ndarray:
    col = [np.ones_like(d["hora"])]
    if K > 0:
        col += _armonicos(d["hora"], K)
    if finde:
        col.append(d["finde"])
        if K > 0:
            col += [d["finde"] * c for c in _armonicos(d["hora"], K)]
    if pulsos is not None:
        col += _pulsos(d["hora"], pulsos)
        if finde:
            col += [d["finde"] * c for c in _pulsos(d["hora"], pulsos)]
    return np.column_stack(col)


# ===========================================================================
# Ajuste fuera de muestra POR DIAS ENTEROS
# ===========================================================================

def r2_oos_por_dias(X, y, dias, n_pliegues=4, semilla=7) -> float:
    """`R2` fuera de muestra con pliegues de DIAS COMPLETOS.

    [!] Partir por casillas sueltas seria fuga: casillas vecinas del mismo dia
    comparten el mismo estado de mercado. El dia entero entra o sale.
    """
    u = np.unique(dias)
    rng = np.random.default_rng(semilla)
    orden = rng.permutation(u)
    pliegues = np.array_split(orden, n_pliegues)
    num = den = 0.0
    for p in pliegues:
        te = np.isin(dias, p)
        tr = ~te
        if tr.sum() < X.shape[1] + 5 or te.sum() < 5:
            continue
        c, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None)
        pred = X[te] @ c
        num += float(np.sum((y[te] - pred) ** 2))
        den += float(np.sum((y[te] - np.mean(y[tr])) ** 2))
    return 1.0 - num / den if den > 0 else float("nan")


def escalera(d: dict, clave: str) -> list:
    y = d[clave]
    out = []
    out.append(("M0  constante", r2_oos_por_dias(diseno(d, 0, False), y, d["dia"])))
    out.append(("M1  armonicos 24 h", r2_oos_por_dias(diseno(d, K_ARM, False), y, d["dia"])))
    out.append(("M2  M1 + nivel de finde",
                r2_oos_por_dias(np.column_stack([diseno(d, K_ARM, False), d["finde"]]),
                                y, d["dia"])))
    out.append(("M3  M2 + amplitud y FASE de finde",
                r2_oos_por_dias(diseno(d, K_ARM, True), y, d["dia"])))
    for tau in TAUS:
        out.append(("M4  M3 + pulsos, tau = %4.2f h" % tau,
                    r2_oos_por_dias(diseno(d, K_ARM, True, tau), y, d["dia"])))
    return out


# ===========================================================================
# Etapas
# ===========================================================================

def etapa_perfil(args) -> int:
    titulo("PERFIL DIURNO DESCRIPTIVO -- `nu` y `sigma` no son escalares")
    d = serie_binada()
    log("")
    log("  casillas de %.0f s : %d   dias distintos: %d"
        % (BIN_S, d["t"].size, np.unique(d["dia"]).size))
    log("")
    log("  hora UTC |   habil: nu[tx/s]  sigma[pb/30s] |   finde: nu[tx/s]  sigma[pb/30s]")
    log("  " + "-" * 76)
    for h in range(24):
        s = (d["hora"] >= h) & (d["hora"] < h + 1)
        a = s & (d["finde"] == 0)
        b = s & (d["finde"] == 1)
        f = lambda m, c: (np.exp(np.mean(d[c][m])) if m.sum() else float("nan"))
        log("     %2d    |        %7.2f       %7.3f    |       %7.2f       %7.3f"
            % (h, f(a, "log_nu"), f(a, "log_sigma"), f(b, "log_nu"), f(b, "log_sigma")))
    for c, nom in (("log_nu", "nu"), ("log_sigma", "sigma")):
        hb = d[c][d["finde"] == 0]
        fd = d[c][d["finde"] == 1]
        log("")
        log("  %-6s habil/finde = %.2fx   (medias geometricas)"
            % (nom, np.exp(hb.mean() - fd.mean())))
    return 0


def etapa_modelo(args) -> int:
    titulo("ESCALERA DE MODELOS, FUERA DE MUESTRA POR DIAS ENTEROS")
    d = serie_binada()
    log("")
    log("  %d casillas de %.0f s sobre %d dias UTC. Pliegues de DIAS COMPLETOS."
        % (d["t"].size, BIN_S, np.unique(d["dia"]).size))
    for clave, nom in (("log_sigma", "log sigma (lo que fija R2_req)"),
                       ("log_nu", "log nu (lo que dimensiona las ventanas)")):
        log("")
        log("--- objetivo: %s" % nom)
        log("    %-36s  R2 fuera de muestra" % "modelo")
        mejor = None
        for etiqueta, r2 in escalera(d, clave):
            log("    %-36s  %+8.4f" % (etiqueta, r2))
            if mejor is None or r2 > mejor[1]:
                mejor = (etiqueta, r2)
        log("    -> mejor: %s   R2 = %+.4f" % mejor)
        # NULO: barajar la etiqueta de dia destruye la relacion hora-valor
        rng = np.random.default_rng(3)
        dn = dict(d)
        perm = rng.permutation(d["hora"].size)
        dn["hora"] = d["hora"][perm]
        dn["finde"] = d["finde"][perm]
        r2n = r2_oos_por_dias(diseno(dn, K_ARM, True, 2.0), d[clave], d["dia"])
        log("    NULO (hora y finde barajados)          %+8.4f" % r2n)
        if mejor[1] > 10 * max(abs(r2n), 1e-4):
            log("    [OK] la estructura diurna es real y sale del nulo con holgura")
        else:
            log("    *** el nulo alcanza al modelo: no se lee ***")
    return 0


# ===========================================================================
# Controles
# ===========================================================================

def _autotest() -> int:
    titulo("CONTROLES DE estacionalidad.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-58s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    rng = np.random.default_rng(1)
    # --- verdad conocida: armonico + pulso + finde desfasado
    n_dias = 14
    t = np.arange(0, n_dias * 24 * 3600, BIN_S) + 1.7e9
    u = [DT.datetime.utcfromtimestamp(x) for x in t]
    d = {"t": t,
         "hora": np.array([x.hour + x.minute / 60.0 for x in u]),
         "dow": np.array([x.weekday() for x in u]),
         "dia": np.array([x.toordinal() for x in u])}
    d["finde"] = (d["dow"] >= 5).astype(float)
    fase_hab, fase_fin = 15.0, 20.0          # el finde llega 5 h mas tarde
    verdad = (1.0 * np.cos(2 * np.pi * (d["hora"] - fase_hab) / 24.0) * (1 - d["finde"])
              + 0.6 * np.cos(2 * np.pi * (d["hora"] - fase_fin) / 24.0) * d["finde"]
              + 0.8 * np.exp(-np.mod(d["hora"] - 13.5, 24.0) / 2.0))
    d["log_sigma"] = verdad + 0.25 * rng.normal(size=t.size)
    d["log_nu"] = d["log_sigma"]

    r0 = r2_oos_por_dias(diseno(d, 0, False), d["log_sigma"], d["dia"])
    r1 = r2_oos_por_dias(diseno(d, K_ARM, False), d["log_sigma"], d["dia"])
    r3 = r2_oos_por_dias(diseno(d, K_ARM, True), d["log_sigma"], d["dia"])
    r4 = r2_oos_por_dias(diseno(d, K_ARM, True, 2.0), d["log_sigma"], d["dia"])
    techo = 1.0 - 0.25 ** 2 / float(np.var(d["log_sigma"]))
    log("     techo por ruido = %.3f | M0 %.3f | M1 %.3f | M3 %.3f | M4 %.3f"
        % (techo, r0, r1, r3, r4))
    chk(abs(r0) < 0.02, "M0 constante da R2 ~ 0 fuera de muestra", "%.4f" % r0)
    chk(r3 > r1 + 0.02, "el bloque de finde recupera el DESFASE de 5 h",
        "M3 - M1 = %+.4f" % (r3 - r1))
    chk(r4 > 0.85 * techo, "M4 se acerca al techo impuesto por el ruido",
        "%.3f contra %.3f" % (r4, techo))

    # --- fase recuperada. DOS defectos propios corregidos aqui, los dos del test:
    # (1) el signo. Con `cos(2*pi*(h-phi)/24)` el coeficiente del coseno es
    #     `A*cos(2*pi*phi/24)` y el del seno `A*sin(2*pi*phi/24)`, asi que
    #     `phi = atan2(c_sin, c_cos)`. Estaba escrito `atan2(-c_sin, c_cos)` y
    #     devolvia 24 - phi. Es la cuarta vez que este proyecto se juega un
    #     resultado en un signo (el 2*pi de la v1.3, el factor 125, la
    #     convencion de eps del propagador).
    # (2) se contrastaba contra una verdad que MEZCLA habil, finde y pulso. Un
    #     control tiene que poner a prueba una sola cosa: se usa una verdad de
    #     un solo armonico y se comprueban las dos fases por separado.
    for etiqueta, fase_v, mascara in (("habil", fase_hab, d["finde"] == 0),
                                      ("finde", fase_fin, d["finde"] == 1)):
        puro = np.cos(2 * np.pi * (d["hora"] - fase_v) / 24.0)
        X = diseno({"hora": d["hora"], "finde": d["finde"]}, 1, False)[mascara]
        c, *_ = np.linalg.lstsq(X, puro[mascara], rcond=None)
        fase = np.mod(np.arctan2(c[2], c[1]) * 24.0 / (2 * np.pi), 24.0)
        err = min(abs(fase - fase_v), 24 - abs(fase - fase_v))
        chk(err < 0.5, "la fase del armonico se recupera (%s)" % etiqueta,
            "%.2f h contra %.2f" % (fase, fase_v))

    # --- negativo: ruido puro
    d2 = dict(d)
    d2["log_sigma"] = rng.normal(size=t.size)
    rn = r2_oos_por_dias(diseno(d2, K_ARM, True, 2.0), d2["log_sigma"], d2["dia"])
    chk(abs(rn) < 0.03, "sobre ruido blanco el modelo completo da R2 ~ 0", "%.4f" % rn)

    # --- negativo: barajar hora destruye la relacion
    perm = rng.permutation(t.size)
    d3 = dict(d)
    d3["hora"] = d["hora"][perm]
    d3["finde"] = d["finde"][perm]
    rb = r2_oos_por_dias(diseno(d3, K_ARM, True, 2.0), d["log_sigma"], d["dia"])
    chk(rb < 0.05 * r4, "barajar la hora destruye la habilidad", "%.4f contra %.4f" % (rb, r4))

    # --- fuga: los pliegues no comparten dia
    dias = d["dia"]
    u_ = np.unique(dias)
    rngf = np.random.default_rng(7)
    pl = np.array_split(rngf.permutation(u_), 4)
    solapan = any(set(a) & set(b) for i, a in enumerate(pl) for b in pl[i + 1:])
    chk(not solapan, "los pliegues no comparten ningun dia (sin fuga)")

    log("")
    log("RESULTADO: %d fallo(s)" % fallos)
    return 1 if fallos else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", default="perfil")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    if a.etapa == "perfil":
        return etapa_perfil(a)
    if a.etapa == "modelo":
        return etapa_modelo(a)
    log("etapa desconocida: %s" % a.etapa)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
