# -*- coding: utf-8 -*-
"""
flujo_omega.py -- la metrica `phi` del flujo de mercado y su derivada `Omega`.

    python flujo_omega.py --autotest                     controles con verdad conocida
    python flujo_omega.py --etapa=serie  [--fuente=...]  construye y cachea la rejilla
    python flujo_omega.py --etapa=dia                    los Omega por dia y por bloque
    python flujo_omega.py --etapa=relacion               relacion entre variables

LA DEFINICION (la que pidio el operador, 2026-09-05):

    phi(t)   = Q_neto(t) * P(t) / tau_0(t)              [USD/s]
    Omega(t) = dphi/dt                                  [USD/s^2]

con `Q_neto` el flujo neto de volumen de mercado firmado [BTC], `P` el precio
[USD/BTC] y `tau_0` el tiempo de restablecimiento del libro [s].

⚠⚠ COLISION DE NOMBRES, Y ES GRAVE SI SE PASA POR ALTO. Este proyecto YA tiene
un `Omega` y un `phi`, los dos RETIRADOS por falta de sustento empirico:

  - el `Omega` de la Sec.1.4 del PDF ([Omega] = BTC/Ticks^2), que alimenta
    `q_inv(Omega)` y `R_base + kappa*Omega^2` en el NMPC, y que depende de
    `omega_m` -- refutada en la v2.2, la v3.0 y la v3.2;
  - `phi'` = ticks/BTC, retirada en el commit `9b2267e` por ser `1/q_bar`.

**El `phi` y el `Omega` de este modulo NO son ninguno de esos dos.** Tienen otras
unidades (USD/s y USD/s^2), otra definicion y otro proposito, y este modulo NO
se importa desde `Micelio.py` ni toca `constantes_micelio.py`. En cualquier acta
se escriben `phi_F` y `Omega_F` (F de flujo) si hay riesgo de confusion.

LA "EDP" ES LA REGLA DE LA CADENA, Y SE MIDE EXACTA. Como las tres entradas son
variables, la derivada total se abre en tres canales aditivos:

    Omega = (P/tau_0)*dQ/dt  +  (Q/tau_0)*dP/dt  -  (Q*P/tau_0^2)*dtau_0/dt
            \____ FLUJO ____/    \___ PRECIO ___/    \______ LIBRO ________/

y, donde `phi != 0`, en forma logaritmica:  Omega/phi = Q'/Q + P'/P - tau_0'/tau_0.

⚠ La descomposicion NO se implementa como aproximacion. Con diferencias
CENTRADAS y medias de los extremos, la regla del producto es una **identidad
algebraica exacta** sobre la rejilla (control 1: residuo < 1e-12 relativo). Eso
importa porque el reparto de varianza entre los tres canales es el resultado que
se va a leer, y un reparto que no suma al total no es un reparto.

⚠ `tau_0` NO SE TOMA DEL AJUSTE DEL NUCLEO. El `tau0` de `migracion_v32.py` es
el parametro de escala del propagador de Bouchaud, y este proyecto lo tiene
**retractado**: sale pegado a su cota inferior y la v4.1 Sec.3.4 lo declaro
"ajuste no identificado en el regimen tau_0 en cota". Ademas es UN numero por
ajuste, y aqui hace falta una SERIE `tau_0(t)`. Se mide directamente del libro,
con tres estimadores independientes (ver `TAU0`), y se reporta la concordancia
entre ellos en vez de elegir uno a ciegas.

⚠ ESTO ES UNA METRICA, NO UNA SENAL. La linea de microestructura de este
proyecto esta CERRADA (v4.2 Sec.4: 78 de 79 celdas negativas; Adenda C: falta un
factor 43x-577x). Nada de este modulo reabre esa decision, y el bloque
predictivo del final va marcado EXPLORATORIO con su nulo al lado, por la
convencion CONTEMPORANEO/PREDICTIVO del proyecto.
"""

from __future__ import annotations

import argparse
import io
import json
import os

import numpy as np

# `horizonte` solo aporta `log`/`titulo` (ASCII) y arrastra numpy; nada mas.
import horizonte as H

log, titulo = H.log, H.titulo

DT_REJILLA = 10.0          # [s] casilla de phi. La misma que usa `cont2014.py`.
BLOQUE_S = 3600.0          # [s] bloque de reporte -> 24 Omega por dia
CACHE = "telemetria/rejilla_omega_%s.npz"
THETAS = (0.30, 0.50, 0.70, 0.90)   # barrido del umbral de agotamiento
TAU_CENSURA_S = 60.0       # por encima de esto el evento de recuperacion se censura
N_SORTEOS_NULO = 200

# Nombres de los tres canales de la regla de la cadena, en orden.
CANALES = ("FLUJO", "PRECIO", "LIBRO")


# ===========================================================================
# 1. La descomposicion exacta de Omega
# ===========================================================================

def omega_tres_factores(t, X, Y, Z):
    """`Omega = d(X*Y*Z)/dt` y sus tres canales, con residuo EXACTAMENTE nulo.

    Diferencias centradas. Para `f = g*h` vale la identidad algebraica

        g[k+1]h[k+1] - g[k-1]h[k-1] = mean(g)*diff(h) + mean(h)*diff(g)

    con `mean` y `diff` tomados sobre los extremos `k-1` y `k+1`. Aplicada dos
    veces (primero separando `X` de `W = Y*Z`, luego `Y` de `Z`) reparte la
    diferencia de `phi` en tres sumandos que suman al total **sin residuo**, no
    hasta O(dt^2). Ojo: `mean(W)` es la media de los PRODUCTOS `Y*Z`, que no es
    `mean(Y)*mean(Z)`; usar esta ultima rompe la identidad.

    Devuelve arrays de longitud n-2, alineados a los indices interiores 1..n-2.
    """
    t = np.asarray(t, float)
    X = np.asarray(X, float); Y = np.asarray(Y, float); Z = np.asarray(Z, float)
    if not (t.shape == X.shape == Y.shape == Z.shape) or t.size < 3:
        raise ValueError("omega_tres_factores: formas incompatibles o n < 3")
    dt = t[2:] - t[:-2]
    Xm, Xd = 0.5 * (X[2:] + X[:-2]), X[2:] - X[:-2]
    Ym, Yd = 0.5 * (Y[2:] + Y[:-2]), Y[2:] - Y[:-2]
    Zm, Zd = 0.5 * (Z[2:] + Z[:-2]), Z[2:] - Z[:-2]
    W = Y * Z
    Wm = 0.5 * (W[2:] + W[:-2])
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = 1.0 / dt
    t1 = Wm * Xd * inv            # canal X
    t2 = Xm * Zm * Yd * inv       # canal Y
    t3 = Xm * Ym * Zd * inv       # canal Z
    phi = X * Y * Z
    om = (phi[2:] - phi[:-2]) * inv
    return {"omega": om, "c1": t1, "c2": t2, "c3": t3,
            "residuo": om - (t1 + t2 + t3), "dt": dt}


def descomponer(t, q_neto, precio, tau0, forma="P/tau"):
    """`phi` y `Omega` con sus canales FLUJO / PRECIO / LIBRO.

    `forma="P/tau"`  ->  phi = Q * P / tau_0        (la pedida)
    `forma="tau/P"`  ->  phi = Q * tau_0 / P        (la variante invertida)

    Las dos se dejan disponibles a proposito: el operador enuncio primero una y
    luego la otra, y la unica diferencia observable entre ambas es el SIGNO del
    canal LIBRO, que es justo lo que la descomposicion mide. Cambiar de forma no
    cuesta un reanalisis; leerlas al reves, si.
    """
    t = np.asarray(t, float)
    q = np.asarray(q_neto, float)
    p = np.asarray(precio, float)
    tau = np.asarray(tau0, float)
    if forma == "P/tau":
        Y, Z = p, 1.0 / tau
    elif forma == "tau/P":
        Y, Z = tau, 1.0 / p
    else:
        raise ValueError("forma desconocida: %r" % (forma,))
    d = omega_tres_factores(t, q, Y, Z)
    return {"t": t[1:-1], "phi": (q * Y * Z)[1:-1], "omega": d["omega"],
            "FLUJO": d["c1"], "PRECIO": d["c2"] if forma == "P/tau" else d["c3"],
            "LIBRO": d["c3"] if forma == "P/tau" else d["c2"],
            "residuo": d["residuo"], "forma": forma}


def reparto_varianza(d):
    """Cuanto de `Var(Omega)` aporta cada canal, contando los cruzados.

    Con `Omega = A + B + C`, `Var(Omega) = sum Var + 2*sum Cov`. Se reporta la
    parte de cada canal como `Cov(canal, Omega)/Var(Omega)`, que **suma 1 por
    construccion** y reparte los terminos cruzados donde corresponde. Las
    fracciones pueden salir negativas: eso significa que el canal se opone al
    total, y es informacion, no un error.
    """
    om = d["omega"]
    ok = np.isfinite(om)
    for c in CANALES:
        ok &= np.isfinite(d[c])
    v = float(np.var(om[ok])) if ok.sum() > 2 else np.nan
    out = {"n": int(ok.sum()), "var_omega": v}
    # [!] LA NORMALIZACION TIENE QUE SER LA MISMA EN LOS DOS SITIOS. `np.cov`
    # divide por n-1 y `np.var` por n; mezclarlas hace que el reparto sume
    # 1 + 1/(n-1) en vez de 1, que con n = 4000 son 2.5e-4 -- invisible a ojo y
    # detectado solo porque el control 5 pide la suma a 1e-10. La covarianza se
    # calcula a mano con ddof = 0.
    om_c = om[ok] - np.mean(om[ok])
    for c in CANALES:
        x = d[c][ok]
        out[c] = (float(np.mean((x - np.mean(x)) * om_c) / v)
                  if ok.sum() > 2 and v > 0 else np.nan)
        out["sd_" + c] = float(np.std(x)) if ok.sum() > 2 else np.nan
    return out


# ===========================================================================
# 2. tau_0 -- tres estimadores del tiempo de restablecimiento del libro
# ===========================================================================
#
# TAU0. Ninguno de los tres es "el" tau_0; se miden los tres y se reporta su
# concordancia. La leccion de la v4.1 Sec.3.1 -- que `gamma` recorria de +0.69 a
# -0.26 segun la banda y que por eso "no existe la gamma de este mercado" -- se
# aplica aqui por adelantado: si los tres estimadores no se mueven juntos, la
# cantidad no esta bien definida y hay que decirlo antes de construir nada
# encima.
#
#   tau_agot  = profundidad L1 / caudal negociado      [s]   PRIMARIO
#               Cuanto dura la cola visible al ritmo de consumo actual. Sin
#               ningun parametro libre: sale de dos cantidades medidas. En
#               regimen estacionario el ritmo de reposicion iguala al de
#               consumo, asi que es la escala de renovacion del nivel 1.
#
#   tau_recup(theta) = tiempo hasta que la profundidad vuelve a su nivel previo
#               tras caer una fraccion >= theta                [s]   LITERAL
#               Es la definicion textual de "restablecimiento". Tiene un umbral
#               libre, asi que se BARRE y se reporta como curva -- este proyecto
#               lleva cinco umbrales puestos a ojo que hubo que corregir
#               midiendo, y este no va a ser el sexto.
#
#               ⚠ Y es DISPERSO: solo esta definido en las casillas que tuvieron
#               un evento de caida >= theta. Con theta = 0.50 sobre la captura de
#               verificacion quedan 314 casillas de 17 041, repartidas en 226
#               rachas -- o sea que como entrada de `phi` fragmenta la serie y la
#               derivada casi no tiene donde apoyarse. Segunda razon para que el
#               primario sea `tau_agot`, que esta definido en todas.
#
#   tau_upd   = 1 / tasa de actualizacion del libro      [s]   DIAGNOSTICO
#               El reloj mas crudo del libro. No es restablecimiento: la mayoria
#               de las actualizaciones de `@bookTicker` son parpadeo de cotizacion.
#               Se reporta para tener una cota inferior de escala, no para usarla.
#
# ⚠ LIMITACION HEREDADA, y no se arregla con analisis: `@bookTicker` da solo el
# NIVEL 1. La reposicion que ocurre en niveles mas profundos es invisible. Ya
# esta medido que la liquidez que gobierna no vive en L1 (H1, 2026-08-23:
# parcial +0.0425 contra techo de nulo +0.0418; y `lambda` = 0.12 en `cont2014`),
# asi que `tau_0` medido aqui es el de la cola visible, no el del libro entero.


def primer_indice_sobre(prof, k, nivel, filas_max, B=64, trozo=100_000):
    """Para cada evento `i`, el primer `j > k[i]` con `prof[j] >= nivel[i]`.

    -1 si no lo hay dentro de `filas_max` filas. Busqueda por MAXIMOS DE BLOQUE:
    se salta un bloque entero de `B` filas cuando su maximo no alcanza el nivel,
    y solo se mira fila a fila el primer bloque que si lo alcanza.

    ⚠ Por que no una pila monotona. La pila resuelve en O(n) el caso "volver al
    nivel EXACTO previo" (`nxt[k-1]`), y esa fue la primera version. Pero el
    nivel objetivo aqui es la SEMI-recuperacion, que es distinto para cada
    evento, y la pila no admite objetivos arbitrarios. Esto si, sigue siendo una
    pasada vectorizada, y ademas impone el tope de busqueda en vez de encontrar
    una recuperacion un millon de filas despues para luego descartarla.
    """
    prof = np.asarray(prof, float)
    k = np.asarray(k, np.int64); nivel = np.asarray(nivel, float)
    n = prof.size
    if k.size == 0:
        return np.full(0, -1, np.int64)
    nbq = int(np.ceil(n / B))
    pm = np.full(nbq * B, -np.inf)
    pm[:n] = prof
    pm = pm.reshape(nbq, B)
    bmax = pm.max(axis=1)
    col = np.arange(B, dtype=np.int64)
    res = np.full(k.size, -1, np.int64)
    for a in range(0, k.size, trozo):
        sl = slice(a, min(a + trozo, k.size))
        kk, lv = k[sl], nivel[sl]
        hi = np.minimum(kk + filas_max, n - 1)
        blk = (kk + 1) // B
        listo = (kk + 1) > hi
        r = np.full(kk.size, -1, np.int64)
        for _ in range(int(np.ceil(filas_max / B)) + 2):
            act = (~listo) & (blk < nbq) & (blk * B <= hi)
            if not act.any():
                break
            cand = act & (bmax[np.minimum(blk, nbq - 1)] >= lv)
            ic = np.flatnonzero(cand)
            if ic.size:
                base = blk[ic] * B
                pos = base[:, None] + col[None, :]
                ok = ((pm[blk[ic]] >= lv[ic][:, None])
                      & (pos > kk[ic][:, None]) & (pos <= hi[ic][:, None]))
                tiene = ok.any(axis=1)
                if tiene.any():
                    j = ic[tiene]
                    r[j] = base[tiene] + np.argmax(ok[tiene], axis=1)
                    listo[j] = True
            blk = blk + 1
            listo |= (blk * B > hi)
        res[sl] = r
    return res


def eventos_recuperacion(t, prof, theta, censura=TAU_CENSURA_S):
    """(t_evento, tau_relajacion, n_censurados) para caidas de profundidad >= `theta`.

    Un evento en `k` es `prof[k] <= (1-theta)*prof[k-1]`. Se mide la SEMI-
    recuperacion: el primer instante en que la profundidad vuelve al punto medio
    entre lo que quedo y lo que habia, `0.5*(prof[k-1] + prof[k])`.

    ⚠ SEMI-RECUPERACION Y NO RECUPERACION TOTAL, y el motivo decide el
    estimador. Bajo relajacion exponencial hacia un equilibrio `D_eq`, la
    profundidad se acerca a su nivel previo pero **no lo alcanza en tiempo
    finito**: el retorno al nivel exacto es un tiempo de primer paso que lo
    dispara el RUIDO de cotizacion, no la reposicion. Medido sobre la captura
    sintetica de verificacion, esa version daba 0.8 s con theta = 0.30 (o sea,
    parpadeo) contra una constante verdadera de 3.4 a 24 s. La semi-recuperacion
    en cambio vale `tau*ln 2` exactamente bajo relajacion exponencial, asi que
    tiene verdad conocida y se puede controlar (control 7).

    Se devuelve ya dividida por `ln 2`, o sea como CONSTANTE DE RELAJACION en
    segundos, para que sea comparable con `tau_agot` sin conversiones al leer.
    """
    t = np.asarray(t, float)
    prof = np.asarray(prof, float)
    if t.size < 3:
        return np.empty(0), np.empty(0), 0
    k = np.flatnonzero(prof[1:] <= (1.0 - theta) * prof[:-1]) + 1
    if k.size == 0:
        return np.empty(0), np.empty(0), 0
    objetivo = 0.5 * (prof[k - 1] + prof[k])
    dtm = float(np.median(np.diff(t))) if t.size > 1 else 1.0
    filas_max = int(min(max(np.ceil(censura / max(dtm, 1e-6)), 8), 20000))
    j = primer_indice_sobre(prof, k, objetivo, filas_max)
    vivo = j >= 0
    tau = np.full(k.size, np.nan)
    tau[vivo] = (t[j[vivo]] - t[k[vivo]]) / np.log(2.0)
    ok = vivo & np.isfinite(tau) & (tau > 0) & (tau <= censura)
    return t[k[ok]], tau[ok], int((~ok).sum())


# ===========================================================================
# 3. Rejilla uniforme: una casilla de `dt` segundos
# ===========================================================================

def _bins(t, t0, dt, nb):
    k = np.floor((np.asarray(t, float) - t0) / dt).astype(np.int64)
    ok = (k >= 0) & (k < nb)
    return k, ok


def _acc_vacio(nb):
    c = ("q_neto", "q_tot", "n_tx", "p_sum", "prof_sum", "n_libro", "rv",
         "mid_ult", "p_ult")
    a = {x: np.zeros(nb) for x in c}
    for th in THETAS:
        a["rec_log_%.2f" % th] = np.zeros(nb)
        a["rec_n_%.2f" % th] = np.zeros(nb)
    # [!] Los censurados son un ESCALAR por umbral, no una serie por casilla.
    # La primera version los sumaba con `+= cens` sobre el array entero, asi que
    # cada casilla recibia el total y el recuento final salia multiplicado por
    # `nb`: 41 869 440 censurados sobre 1.7 M de filas de libro. Lo delato que
    # el numero era imposible, no una asercion -- por eso ahora hay una
    # (control 13).
    a["cens"] = {th: 0 for th in THETAS}
    return a


def _ultimo_por_casilla(k, v, dest):
    """Escribe en `dest[k]` el ULTIMO `v` de cada casilla (cierre)."""
    if k.size == 0:
        return
    o = np.argsort(k, kind="stable")
    kk, vv = k[o], v[o]
    ult = np.r_[np.flatnonzero(np.diff(kk)), kk.size - 1]
    dest[kk[ult]] = vv[ult]


def acumular_libro(bt, bb, bB, ba, bA, t0, dt, nb, acc, prev=None):
    """Profundidad, mid de cierre, volatilidad realizada y eventos de recuperacion.

    `prev` es la ultima fila de la parte anterior, para que el primer incremento
    de cada parte no se pierda ni se invente.
    """
    if bt.size == 0:
        return prev
    if prev is not None:
        bt = np.r_[prev[0], bt]; bb = np.r_[prev[1], bb]; bB = np.r_[prev[2], bB]
        ba = np.r_[prev[3], ba]; bA = np.r_[prev[4], bA]
    mid = 0.5 * (bb + ba)
    prof = 0.5 * (bB + bA)
    k, ok = _bins(bt, t0, dt, nb)
    if ok.any():
        acc["prof_sum"] += np.bincount(k[ok], weights=prof[ok], minlength=nb)
        acc["n_libro"] += np.bincount(k[ok], minlength=nb)
        _ultimo_por_casilla(k[ok], mid[ok], acc["mid_ult"])
    # volatilidad realizada: se asigna el incremento a la casilla de su EXTREMO
    if bt.size > 1:
        r = np.diff(np.log(mid))
        kr, okr = k[1:], ok[1:]
        m = okr & np.isfinite(r)
        if m.any():
            acc["rv"] += np.bincount(kr[m], weights=r[m] ** 2, minlength=nb)
    for th in THETAS:
        te, tau, cens = eventos_recuperacion(bt, prof, th)
        if te.size:
            ke, oke = _bins(te, t0, dt, nb)
            if oke.any():
                acc["rec_log_%.2f" % th] += np.bincount(
                    ke[oke], weights=np.log(tau[oke]), minlength=nb)
                acc["rec_n_%.2f" % th] += np.bincount(ke[oke], minlength=nb)
        acc["cens"][th] += cens          # ESCALAR: ver la nota en `_acc_vacio`
    return (float(bt[-1]), float(bb[-1]), float(bB[-1]),
            float(ba[-1]), float(bA[-1]))


def acumular_trades(tt, pr, q, maker, t0, dt, nb, acc):
    """Flujo neto firmado, volumen total, recuento y precio de cierre."""
    if tt.size == 0:
        return
    eps = np.where(np.asarray(maker, bool), -1.0, 1.0)   # convencion de `propagador`
    k, ok = _bins(tt, t0, dt, nb)
    if not ok.any():
        return
    acc["q_neto"] += np.bincount(k[ok], weights=(eps * q)[ok], minlength=nb)
    acc["q_tot"] += np.bincount(k[ok], weights=q[ok], minlength=nb)
    acc["n_tx"] += np.bincount(k[ok], minlength=nb)
    acc["p_sum"] += np.bincount(k[ok], weights=pr[ok], minlength=nb)
    _ultimo_por_casilla(k[ok], pr[ok], acc["p_ult"])


def cerrar_rejilla(acc, t0, dt, nb):
    """De los acumuladores a las variables del modelo, con su mascara de validez."""
    n_lib = acc["n_libro"]; n_tx = acc["n_tx"]; q_tot = acc["q_tot"]
    with np.errstate(divide="ignore", invalid="ignore"):
        prof = np.where(n_lib > 0, acc["prof_sum"] / np.maximum(n_lib, 1), np.nan)
        caudal = q_tot / dt                                   # [BTC/s]
        tau_agot = np.where(caudal > 0, prof / np.maximum(caudal, 1e-30), np.nan)
        tau_upd = np.where(n_lib > 0, dt / np.maximum(n_lib, 1), np.nan)
    r = {"t": t0 + (np.arange(nb) + 0.5) * dt, "dt": dt,
         "q_neto": acc["q_neto"], "q_tot": q_tot, "n_tx": n_tx,
         "nu": n_tx / dt, "prof": prof, "n_libro": n_lib,
         "mid": np.where(acc["mid_ult"] > 0, acc["mid_ult"], np.nan),
         "precio_tx": np.where(acc["p_ult"] > 0, acc["p_ult"], np.nan),
         "rv_pb": np.sqrt(np.maximum(acc["rv"], 0.0)) * 1e4,
         "tau_agot": tau_agot, "tau_upd": tau_upd}
    for th in THETAS:
        n = acc["rec_n_%.2f" % th]
        with np.errstate(divide="ignore", invalid="ignore"):
            r["tau_recup_%.2f" % th] = np.where(
                n > 0, np.exp(acc["rec_log_%.2f" % th] / np.maximum(n, 1)), np.nan)
        r["rec_n_%.2f" % th] = n
        r["rec_cens_%.2f" % th] = np.array([float(acc["cens"][th])])
    # El precio del modelo es el MID de cierre; el de transaccion es el respaldo.
    r["precio"] = np.where(np.isfinite(r["mid"]), r["mid"], r["precio_tx"])
    r["valida"] = (n_tx > 0) & (n_lib > 0) & (q_tot > 0) \
        & np.isfinite(r["precio"]) & np.isfinite(tau_agot) & (tau_agot > 0)
    return r


# ===========================================================================
# 4. Fuentes de datos
# ===========================================================================
#
# Una SOLA pasada por parte, sobre una rejilla global de toda la captura. Con
# `dt` = 10 s, tres semanas son 181 k casillas y ~30 MB: lo que no cabe es el
# dato crudo (223 M de filas de libro), no la rejilla. Por eso se acumula
# parte a parte y nunca hay mas de una en memoria -- la misma disciplina que
# `curvas_estacional.limites_tramos` tuvo que adoptar tras morir tres procesos.


def _rango_estacional():
    import curvas_estacional as C
    idx = C._indice("trades_")
    if not idx:
        raise SystemExit("no hay partes de trades en %s" % C.DIR)
    return min(x[1] for x in idx), max(x[2] for x in idx)


def rejilla_estacional(dt=DT_REJILLA, t_ini=None, t_fin=None):
    import curvas_estacional as C
    import pyarrow.parquet as pq
    a, b = _rango_estacional()
    t0 = a if t_ini is None else max(a, t_ini)
    t1 = b if t_fin is None else min(b, t_fin)
    nb = int(np.ceil((t1 - t0) / dt))
    log("  ventana: %.2f dias   casillas de %.0f s: %d" % ((t1 - t0) / 86400.0, dt, nb))
    acc = _acc_vacio(nb)
    prev = None
    partes = sorted([x for x in C._indice("libro_") if not (x[2] < t0 or x[1] > t1)],
                    key=lambda z: z[1])
    log("  partes de libro que solapan: %d" % len(partes))
    for i, (f, _, _, _) in enumerate(partes):
        try:
            tb = pq.read_table(f, columns=["t", "b", "B", "a", "A"])
        except Exception:
            continue
        bt = tb["t"].to_numpy().astype(float)
        bb = tb["b"].to_numpy().astype(float); bB = tb["B"].to_numpy().astype(float)
        ba = tb["a"].to_numpy().astype(float); bA = tb["A"].to_numpy().astype(float)
        del tb
        m = (bb > 0) & (ba > 0) & (bB > 0) & (bA > 0) & (bt >= t0) & (bt <= t1)
        if m.sum() < 2:
            continue
        bt, bb, bB, ba, bA = bt[m], bb[m], bB[m], ba[m], bA[m]
        o = np.argsort(bt, kind="stable")
        prev = acumular_libro(bt[o], bb[o], bB[o], ba[o], bA[o], t0, dt, nb, acc, prev)
        del bt, bb, bB, ba, bA, o, m
        if (i + 1) % 50 == 0:
            log("    libro %d/%d" % (i + 1, len(partes)))
    partes = sorted([x for x in C._indice("trades_") if not (x[2] < t0 or x[1] > t1)],
                    key=lambda z: z[1])
    log("  partes de trades que solapan: %d" % len(partes))
    for f, _, _, _ in partes:
        try:
            tb = pq.read_table(f, columns=["t", "precio", "cant", "maker"])
        except Exception:
            continue
        tt = tb["t"].to_numpy().astype(float)
        pr = tb["precio"].to_numpy().astype(float)
        q = tb["cant"].to_numpy().astype(float)
        mk = tb["maker"].to_numpy().astype(bool)
        del tb
        m = (pr > 0) & (q > 0) & (tt >= t0) & (tt <= t1)
        if m.any():
            acumular_trades(tt[m], pr[m], q[m], mk[m], t0, dt, nb, acc)
        del tt, pr, q, mk, m
    return cerrar_rejilla(acc, t0, dt, nb)


def rejilla_v33(dt=DT_REJILLA):
    """`captura_v33` -- 16.4 h en un solo tramo continuo. Cabe entera en memoria."""
    from captura_larga import cargar_larga
    d = cargar_larga("telemetria/captura_v33")
    bt, bb, bB = d["bk_t"], d["bk_b"], d["bk_B"]
    ba, bA = d["bk_a"], d["bk_A"]
    m = (bb > 0) & (ba > 0) & (bB > 0) & (bA > 0)
    bt, bb, bB, ba, bA = bt[m], bb[m], bB[m], ba[m], bA[m]
    o = np.argsort(bt, kind="stable")
    bt, bb, bB, ba, bA = bt[o], bb[o], bB[o], ba[o], bA[o]
    tt, pr, q = d["tr_t"], d["tr_precio"], d["tr_cant"]
    mk = d["tr_maker"].astype(bool)
    okt = np.isfinite(pr) & (pr > 0) & (q > 0)
    t0, t1 = float(min(bt[0], tt[okt][0])), float(max(bt[-1], tt[okt][-1]))
    nb = int(np.ceil((t1 - t0) / dt))
    log("  ventana: %.2f h   casillas de %.0f s: %d" % ((t1 - t0) / 3600.0, dt, nb))
    acc = _acc_vacio(nb)
    acumular_libro(bt, bb, bB, ba, bA, t0, dt, nb, acc, None)
    acumular_trades(tt[okt], pr[okt], q[okt], mk[okt], t0, dt, nb, acc)
    return cerrar_rejilla(acc, t0, dt, nb)


def cargar_rejilla(fuente):
    r = CACHE % fuente
    if not os.path.exists(r):
        raise SystemExit("falta %s -- corre primero --etapa=serie --fuente=%s"
                         % (r, fuente))
    d = np.load(r)
    return {k: (float(d[k]) if d[k].ndim == 0 else d[k]) for k in d.files}


# ===========================================================================
# 5. Agregacion por dia y por bloque -- "varios Omega por dia"
# ===========================================================================

def _dia_hora(t):
    d = np.floor(np.asarray(t, float) / 86400.0).astype(np.int64)
    return d, (np.asarray(t, float) - d * 86400.0)


def _fecha(dia):
    import datetime as _dt
    return (_dt.datetime(1970, 1, 1) + _dt.timedelta(days=int(dia))).strftime("%Y-%m-%d")


def serie_omega(r, tau="tau_agot", forma="P/tau"):
    """`phi` y `Omega` sobre las casillas VALIDAS y CONTIGUAS de la rejilla.

    ⚠ La contiguidad importa y por eso se corta. La captura tiene huecos (el de
    24.26 h del 2026-08-21 entre otros) y casillas sin transacciones; una
    diferencia centrada que cruce un hueco no es una derivada, es la pendiente
    entre dos regimenes distintos. Se parte la rejilla en RACHAS de casillas
    validas consecutivas y se deriva dentro de cada racha.
    """
    v = r["valida"] & np.isfinite(r[tau]) & (r[tau] > 0)
    idx = np.flatnonzero(v)
    if idx.size < 3:
        return {k: np.empty(0) for k in
                ("t", "phi", "omega") + CANALES + ("residuo",)}
    corte = np.flatnonzero(np.diff(idx) != 1)
    grupos = np.split(idx, corte + 1)
    piezas = []
    for g in grupos:
        if g.size < 3:
            continue
        d = descomponer(r["t"][g], r["q_neto"][g], r["precio"][g], r[tau][g], forma)
        d["idx"] = g[1:-1]
        piezas.append(d)
    if not piezas:
        return {k: np.empty(0) for k in
                ("t", "phi", "omega") + CANALES + ("residuo",)}
    out = {k: np.concatenate([p[k] for p in piezas])
           for k in ("t", "phi", "omega", "residuo") + CANALES}
    out["idx"] = np.concatenate([p["idx"] for p in piezas])
    out["forma"] = forma
    out["tau"] = tau
    out["n_rachas"] = len(piezas)
    return out


def por_bloque(r, s, bloque=BLOQUE_S):
    """Una fila por (dia, bloque). Es la tabla que el operador pidio."""
    if s["t"].size == 0:
        return []
    dia, seg = _dia_hora(s["t"])
    blq = np.floor(seg / bloque).astype(np.int64)
    clave = dia * 10000 + blq
    orden = np.argsort(clave, kind="stable")
    c = clave[orden]
    cortes = np.r_[0, np.flatnonzero(np.diff(c)) + 1, c.size]
    filas = []
    for i in range(cortes.size - 1):
        sl = orden[cortes[i]:cortes[i + 1]]
        if sl.size < 3:
            continue
        om = s["omega"][sl]; ph = s["phi"][sl]
        j = s["idx"][sl]
        f = {"dia": int(dia[sl[0]]), "fecha": _fecha(dia[sl[0]]),
             "bloque": int(blq[sl[0]]), "n": int(sl.size),
             "t0": float(s["t"][sl].min()),
             # los tres Omega del bloque, y son cantidades DISTINTAS:
             "omega_neta": float((ph[-1] - ph[0]) / max(s["t"][sl][-1] - s["t"][sl][0], 1e-9)),
             "omega_rms": float(np.sqrt(np.mean(om ** 2))),
             "omega_med": float(np.median(om)),
             "phi_med": float(np.median(ph)), "phi_sd": float(np.std(ph)),
             # variables de estado del bloque
             "q_neto": float(np.sum(r["q_neto"][j])),
             "q_tot": float(np.sum(r["q_tot"][j])),
             "precio": float(np.median(r["precio"][j])),
             "tau0": float(np.median(r[s["tau"]][j])),
             "prof": float(np.median(r["prof"][j])),
             "nu": float(np.median(r["nu"][j])),
             "rv_pb": float(np.sqrt(np.mean(r["rv_pb"][j] ** 2)))}
        rep = reparto_varianza({"omega": om, **{c2: s[c2][sl] for c2 in CANALES}})
        for c2 in CANALES:
            f["f_" + c2] = rep[c2]
        filas.append(f)
    return sorted(filas, key=lambda z: (z["dia"], z["bloque"]))


# ===========================================================================
# 6. Relacion entre variables
# ===========================================================================

def _rango(x):
    o = np.argsort(x, kind="stable")
    rk = np.empty(x.size, float)
    rk[o] = np.arange(x.size, dtype=float)
    # empates al promedio: sin esto Spearman se sesga con muchas casillas iguales
    xs = x[o]
    i = 0
    while i < xs.size:
        j = i
        while j + 1 < xs.size and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            rk[o[i:j + 1]] = 0.5 * (i + j)
        i = j + 1
    return rk


def corr(x, y, metodo="spearman"):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4:
        return np.nan
    x, y = x[ok], y[ok]
    if metodo == "spearman":
        x, y = _rango(x), _rango(y)
    sx, sy = np.std(x), np.std(y)
    if sx <= 0 or sy <= 0:
        return np.nan
    return float(np.mean((x - np.mean(x)) * (y - np.mean(y))) / (sx * sy))


def suelo_rotacion(x, y, n=N_SORTEOS_NULO, semilla=17, metodo="spearman"):
    """q95 de |corr| bajo ROTACION CIRCULAR de `y` contra `x`.

    ⚠ Rotacion y no barajado, y la diferencia no es cosmetica. La Adenda C midio
    que barajar subestima el suelo entre 75x y 292x porque destruye la memoria
    de la serie, que es lo que ensancha al estadistico. Estas series (`nu`,
    `tau_0`, `rv`) son fuertemente persistentes -- `rho_1(log nu) = 0.8845`
    medido el 2026-08-27 -- asi que aqui el barajado seria igual de enganoso.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    m = x.size
    if m < 20:
        return np.nan
    rng = np.random.default_rng(semilla)
    v = np.empty(n)
    for i in range(n):
        k = int(rng.integers(max(3, m // 20), m - max(3, m // 20)))
        v[i] = abs(corr(x, np.roll(y, k), metodo))
    return float(np.nanpercentile(v, 95))


def quitar_ciclo(x, dia, blq, en_log=True, iters=3):
    """Residuo de `x` tras retirar el NIVEL DEL DIA y la FORMA INTRADIA.

    ⚠ ESTO NO ES COSMETICA, ES LA CONDICION PARA QUE EL NULO SIGNIFIQUE ALGO.
    Casi todas las variables de aqui (`nu`, `tau_0`, `prof`, `rv`) llevan el
    mismo ciclo diurno de 24 h, medido el 2026-08-23: `nu` recorre 3.9x entre
    UTC 13-15 y UTC 4, y `sigma` 2.7x. Dos series con el mismo periodo conocido
    correlacionan a +0.99 **por tener el mismo periodo**, y una rotacion
    circular las vuelve a alinear en la siguiente vuelta -- por eso el suelo de
    rotacion sale tambien en 0.99 y la razon en 1.00. La correlacion cruda entre
    dos series diurnas no informa de nada; el nulo lo dice correctamente y hay
    que hacerle caso.

    Aquel dia fallaron TRES nulos propios por exactamente esto, y la regla que
    quedo escrita es: separar nivel diario y forma intradia ANTES de construir
    cualquier nulo. Ademas, el 2026-08-23 midio que el nivel del DIA es la
    palanca grande (`log sigma` recorre 13.75x entre dias, 46 % de su varianza),
    asi que retirarlo no es quitar una molestia: es quitar el efecto dominante
    para poder ver el resto.

    Se retiran las dos medias marginales de forma alternada (Gauss-Seidel), en
    logaritmo cuando la serie es estrictamente positiva.
    """
    x = np.asarray(x, float)
    y = np.log(x) if (en_log and np.all(np.isfinite(x)) and np.all(x > 0)) else x.copy()
    ok = np.isfinite(y)
    if ok.sum() < 8:
        return np.full(x.size, np.nan)
    out = np.where(ok, y - np.mean(y[ok]), np.nan)
    for _ in range(iters):
        for clave in (np.asarray(dia), np.asarray(blq)):
            for v in np.unique(clave[ok]):
                m = ok & (clave == v)
                if m.sum():
                    out[m] -= np.mean(out[m])
    return out


def informe_tau0(r):
    titulo("tau_0 -- CONCORDANCIA ENTRE LOS TRES ESTIMADORES")
    v = r["valida"]
    log("  casillas validas: %d de %d (%.1f %%)" % (v.sum(), v.size, 100 * v.mean()))
    log("")
    log("  %-22s %10s %10s %10s %10s" % ("estimador", "p10", "MEDIANA", "p90", "n"))
    cols = [("tau_agot [s]", "tau_agot"), ("tau_upd  [s]", "tau_upd")]
    cols += [("tau_recup(%.2f)" % th, "tau_recup_%.2f" % th) for th in THETAS]
    for nom, c in cols:
        x = r[c][v]
        x = x[np.isfinite(x) & (x > 0)]
        if x.size < 10:
            log("  %-22s %10s %10s %10s %10d" % (nom, "-", "-", "-", x.size))
            continue
        log("  %-22s %10.3f %10.3f %10.3f %10d"
            % (nom, np.percentile(x, 10), np.median(x), np.percentile(x, 90), x.size))
    log("")
    log("  correlacion de Spearman entre estimadores (casillas validas):")
    base = r["tau_agot"]
    for nom, c in cols[1:]:
        log("    tau_agot vs %-18s  rho = %+.4f" % (nom, corr(base[v], r[c][v])))
    log("")
    for th in THETAS:
        n = r["rec_n_%.2f" % th][v].sum()
        cn = r["rec_cens_%.2f" % th].sum()
        log("    theta = %.2f : %9.0f eventos con recuperacion, %8.0f censurados (%.1f %%)"
            % (th, n, cn, 100 * cn / max(n + cn, 1)))
    log("")
    log("  [!] `tau_agot` es el PRIMARIO: no tiene parametro libre. `tau_recup`")
    log("      es la definicion literal y lleva umbral, por eso va barrido.")
    log("      `tau_upd` es diagnostico: mide parpadeo de cotizacion, no reposicion.")


def informe_relacion(filas, s, forma):
    titulo("RELACION ENTRE VARIABLES -- nivel de BLOQUE")
    if len(filas) < 20:
        log("  solo %d bloques: por debajo de 20 no se lee nada." % len(filas))
        return
    col = {k: np.array([f[k] for f in filas], float)
           for k in ("omega_rms", "omega_neta", "phi_med", "phi_sd", "q_neto",
                     "q_tot", "precio", "tau0", "prof", "nu", "rv_pb",
                     "f_FLUJO", "f_PRECIO", "f_LIBRO")}
    col["abs_q_neto"] = np.abs(col["q_neto"])
    col["abs_omega_neta"] = np.abs(col["omega_neta"])
    dia = np.array([f["dia"] for f in filas])
    blq = np.array([f["bloque"] for f in filas])
    SIGNADAS = ("q_neto", "omega_neta", "f_FLUJO", "f_PRECIO", "f_LIBRO")
    res = {k: quitar_ciclo(v, dia, blq, en_log=(k not in SIGNADAS))
           for k, v in col.items()}
    nd = int(np.unique(dia).size)
    nb_u = int(np.unique(blq).size)
    log("  n = %d bloques   %d dias   forma phi = Q*%s" % (len(filas), nd, forma))
    if len(filas) < 4 * (nd + nb_u):
        log("")
        log("  [!] AVISO: %d bloques contra %d parametros de ciclo (%d dias + %d"
            % (len(filas), nd + nb_u - 1, nd, nb_u))
        log("      casillas horarias). La columna RESIDUO esta SOBRE-restada: se")
        log("      esta quitando casi tanta varianza como hay. Hacen falta >= 4")
        log("      observaciones por parametro, o sea >= %d bloques."
            % (4 * (nd + nb_u)))
    log("")
    log("  --- CONTEMPORANEO: Spearman contra el suelo de ROTACION (q95 de |rho|) ---")
    log("      CRUDO  = las series tal cual, con su ciclo diurno dentro")
    log("      RESIDUO = sin el nivel del dia ni la forma intradia (ver `quitar_ciclo`)")
    log("")
    log("  %-30s %8s %8s %6s   %8s %8s %6s"
        % ("par", "rho", "suelo", "razon", "rho_res", "suelo", "razon"))
    pares = [("omega_rms", "nu"), ("omega_rms", "rv_pb"), ("omega_rms", "q_tot"),
             ("omega_rms", "tau0"), ("omega_rms", "prof"), ("omega_rms", "precio"),
             ("phi_sd", "nu"), ("phi_sd", "rv_pb"),
             ("tau0", "nu"), ("tau0", "rv_pb"), ("prof", "nu"),
             ("abs_q_neto", "rv_pb"), ("q_neto", "rv_pb"),
             ("f_LIBRO", "nu"), ("f_LIBRO", "tau0")]
    for a, b in pares:
        rho = corr(col[a], col[b]); su = suelo_rotacion(col[a], col[b])
        rz = abs(rho) / su if (np.isfinite(su) and su > 0) else np.nan
        rr = corr(res[a], res[b]); sr = suelo_rotacion(res[a], res[b])
        rzr = abs(rr) / sr if (np.isfinite(sr) and sr > 0) else np.nan
        m1 = "<<" if np.isfinite(rzr) and rzr >= 3.0 else "  "
        log("  %-30s %+8.4f %8.4f %6.2f   %+8.4f %8.4f %6.2f %s"
            % ("%s vs %s" % (a, b), rho, su, rz, rr, sr, rzr, m1))
    log("")
    log("  [!] Se lee una fila solo si `razon >= 3` en la columna RESIDUO (`<<`).")
    log("      Un `razon ~ 1` en la columna CRUDO con `rho ~ 0.99` NO es una")
    log("      asociacion fuerte: es que las dos series tienen el mismo periodo")
    log("      de 24 h, y una rotacion circular vuelve a alinearlas.")
    log("")
    log("  --- reparto de Var(Omega) por canal, entre bloques ---")
    log("  %-10s %10s %10s %10s" % ("canal", "p10", "MEDIANA", "p90"))
    for c in CANALES:
        x = col["f_" + c]
        x = x[np.isfinite(x)]
        if x.size:
            log("  %-10s %10.4f %10.4f %10.4f"
                % (c, np.percentile(x, 10), np.median(x), np.percentile(x, 90)))
    log("")
    log("  --- PREDICTIVO EXPLORATORIO: Omega del bloque contra el retorno del")
    log("      bloque SIGUIENTE. NO es un veredicto y no reabre nada. ---")
    pr = col["precio"]
    ret = np.full(pr.size, np.nan)
    ret[:-1] = np.diff(np.log(pr)) * 1e4          # pb
    cont = np.full(pr.size, np.nan)
    cont[1:] = np.diff(np.log(pr)) * 1e4
    log("  %-40s %9s %9s %7s" % ("par", "rho", "suelo", "razon"))
    for nom, a, b in (("CONTEMPORANEO omega_neta vs ret(t)", col["omega_neta"], cont),
                      ("PREDICTIVO    omega_neta vs ret(t+1)", col["omega_neta"], ret),
                      ("PREDICTIVO    q_neto     vs ret(t+1)", col["q_neto"], ret),
                      ("PREDICTIVO    phi_med    vs ret(t+1)", col["phi_med"], ret)):
        rho = corr(a, b)
        su = suelo_rotacion(a, b)
        raz = abs(rho) / su if (np.isfinite(su) and su > 0) else np.nan
        marca = "  <<" if np.isfinite(raz) and raz >= 3.0 else ""
        log("  %-40s %+9.4f %9.4f %7.2f%s" % (nom, rho, su, raz, marca))


# ===========================================================================
# 7. Etapas
# ===========================================================================

def etapa_serie(args) -> int:
    titulo("REJILLA -- una casilla de %.0f s sobre `%s`" % (args.dt, args.fuente))
    if args.fuente == "estacional":
        r = rejilla_estacional(args.dt)
    elif args.fuente == "v33":
        r = rejilla_v33(args.dt)
    else:
        raise SystemExit("fuente desconocida: %r" % (args.fuente,))
    os.makedirs("telemetria", exist_ok=True)
    np.savez_compressed(CACHE % args.fuente, **r)
    log("")
    log("  guardado en %s" % (CACHE % args.fuente))
    v = r["valida"]
    log("  casillas %d, validas %d (%.1f %%), transacciones %d"
        % (v.size, v.sum(), 100 * v.mean(), int(r["n_tx"].sum())))
    if v.any():
        log("  nu mediana %.2f tx/s   precio %.1f -- %.1f   profundidad L1 mediana %.3f BTC"
            % (np.median(r["nu"][v]), np.nanmin(r["precio"][v]),
               np.nanmax(r["precio"][v]), np.nanmedian(r["prof"][v])))
    informe_tau0(r)
    return 0


def etapa_dia(args) -> int:
    r = cargar_rejilla(args.fuente)
    informe_tau0(r)
    s = serie_omega(r, tau=args.tau, forma=args.forma)
    titulo("phi Y Omega -- fuente `%s`, tau = `%s`, forma phi = Q*%s"
           % (args.fuente, args.tau, args.forma))
    if s["t"].size == 0:
        log("  no hay rachas validas de 3 casillas. Nada que derivar.")
        return 2
    res = np.max(np.abs(s["residuo"])) / max(np.max(np.abs(s["omega"])), 1e-30)
    log("  casillas con Omega: %d   rachas contiguas: %d" % (s["t"].size, s["n_rachas"]))
    log("  residuo de la regla de la cadena: %.3e relativo   (control 1: < 1e-9)" % res)
    log("  phi   [USD/s]  : p10 %+.4e  MED %+.4e  p90 %+.4e"
        % tuple(np.percentile(s["phi"], [10, 50, 90])))
    log("  Omega [USD/s^2]: p10 %+.4e  MED %+.4e  p90 %+.4e"
        % tuple(np.percentile(s["omega"], [10, 50, 90])))
    rep = reparto_varianza(s)
    log("")
    log("  reparto GLOBAL de Var(Omega):  " + "   ".join(
        "%s %+.4f" % (c, rep[c]) for c in CANALES)
        + "   (suma %+.4f)" % sum(rep[c] for c in CANALES))
    filas = por_bloque(r, s, args.bloque)
    titulo("LOS Omega POR DIA -- %d bloques de %.0f min" % (len(filas), args.bloque / 60.0))
    log("  omega_neta = (phi_fin - phi_ini)/T del bloque  [USD/s^2]  -- deriva neta")
    log("  omega_rms  = raiz de la media de Omega^2       [USD/s^2]  -- actividad")
    log("  f_X        = parte de Var(Omega) del canal X (los tres suman 1)")
    log("")
    dias = sorted(set(f["dia"] for f in filas))
    for d in dias:
        fs = [f for f in filas if f["dia"] == d]
        log("")
        log("  --- %s : %d bloques ---" % (fs[0]["fecha"], len(fs)))
        log("  %3s %6s %12s %12s %10s %9s %8s %7s %7s %7s %7s"
            % ("blq", "n", "omega_neta", "omega_rms", "q_neto", "tau0", "nu",
               "rv_pb", "fFLU", "fPRE", "fLIB"))
        for f in fs:
            log("  %3d %6d %+12.4e %12.4e %+10.3f %9.4f %8.2f %7.3f %+7.3f %+7.3f %+7.3f"
                % (f["bloque"], f["n"], f["omega_neta"], f["omega_rms"], f["q_neto"],
                   f["tau0"], f["nu"], f["rv_pb"],
                   f["f_FLUJO"], f["f_PRECIO"], f["f_LIBRO"]))
    titulo("RESUMEN POR DIA -- el cambio de un dia a otro")
    log("  %-12s %5s %13s %13s %10s %9s %8s %7s %7s"
        % ("fecha", "blq", "omega_rms MED", "omega_rms p90", "q_neto", "tau0",
           "nu", "rv_pb", "fLIB"))
    for d in dias:
        fs = [f for f in filas if f["dia"] == d]
        g = lambda c: np.array([f[c] for f in fs], float)
        log("  %-12s %5d %13.4e %13.4e %+10.2f %9.4f %8.2f %7.3f %+7.3f"
            % (fs[0]["fecha"], len(fs), np.median(g("omega_rms")),
               np.percentile(g("omega_rms"), 90), np.sum(g("q_neto")),
               np.median(g("tau0")), np.median(g("nu")), np.median(g("rv_pb")),
               np.median(g("f_LIBRO"))))
    if len(dias) > 1:
        rr = np.array([np.median([f["omega_rms"] for f in filas if f["dia"] == d])
                       for d in dias], float)
        log("")
        log("  recorrido de omega_rms entre dias: %.2fx  (min %.3e, max %.3e)"
            % (np.max(rr) / max(np.min(rr), 1e-30), np.min(rr), np.max(rr)))

    csv = "telemetria/omega_bloques_%s.csv" % args.fuente
    cab = ["fecha", "dia", "bloque", "n", "t0", "omega_neta", "omega_rms", "omega_med",
           "phi_med", "phi_sd", "q_neto", "q_tot", "precio", "tau0", "prof", "nu",
           "rv_pb", "f_FLUJO", "f_PRECIO", "f_LIBRO"]
    with io.open(csv, "w", encoding="ascii") as fh:
        fh.write(",".join(cab) + "\n")
        for f in filas:
            fh.write(",".join(("%s" % f[c]) if c == "fecha" else ("%.10g" % f[c])
                              for c in cab) + "\n")
    log("")
    log("  tabla completa en %s" % csv)
    json.dump({"fuente": args.fuente, "tau": args.tau, "forma": args.forma,
               "dt": args.dt, "bloque": args.bloque, "n_bloques": len(filas),
               "residuo_rel": res, "reparto": {c: rep[c] for c in CANALES}},
              io.open("telemetria/omega_resumen_%s.json" % args.fuente, "w",
                      encoding="ascii"))
    return 0


def etapa_relacion(args) -> int:
    r = cargar_rejilla(args.fuente)
    s = serie_omega(r, tau=args.tau, forma=args.forma)
    if s["t"].size == 0:
        log("  no hay serie. Corre --etapa=serie primero.")
        return 2
    filas = por_bloque(r, s, args.bloque)
    informe_relacion(filas, s, args.forma)
    return 0


# ===========================================================================
# 8. Controles con verdad conocida
# ===========================================================================

def _autotest() -> int:
    titulo("flujo_omega.py -- CONTROLES")
    n_ok = n_tot = 0

    def chk(ok, msg, det=""):
        nonlocal n_ok, n_tot
        n_tot += 1
        n_ok += bool(ok)
        log("  [%s] %s%s" % ("OK " if ok else "FALLA", msg,
                             ("   %s" % det) if det else ""))

    rng = np.random.default_rng(20260905)

    # --- 1. La regla de la cadena es EXACTA, no aproximada -------------------
    n = 4000
    t = np.arange(n) * 7.0
    Q = rng.normal(0, 4, n)
    P = 60000 + np.cumsum(rng.normal(0, 8, n))
    T = 0.5 + np.abs(rng.normal(3, 1.2, n))
    for forma in ("P/tau", "tau/P"):
        d = descomponer(t, Q, P, T, forma)
        rel = np.max(np.abs(d["residuo"])) / np.max(np.abs(d["omega"]))
        chk(rel < 1e-12, "1 regla de la cadena exacta (forma %s)" % forma,
            "residuo relativo %.2e" % rel)

    # --- 2. Recupera una derivada analitica ---------------------------------
    dt = 0.01
    tt = np.arange(0, 20, dt)
    Qa = 2.0 + np.sin(0.7 * tt)
    Pa = 60000.0 * (1.0 + 0.01 * np.cos(0.3 * tt))
    Ta = 3.0 + 0.5 * np.sin(0.11 * tt)
    dQ = 0.7 * np.cos(0.7 * tt)
    dP = -60000.0 * 0.01 * 0.3 * np.sin(0.3 * tt)
    dT = 0.5 * 0.11 * np.cos(0.11 * tt)
    exact = (Pa / Ta) * dQ + (Qa / Ta) * dP - (Qa * Pa / Ta ** 2) * dT
    d = descomponer(tt, Qa, Pa, Ta)
    err = np.max(np.abs(d["omega"] - exact[1:-1])) / np.max(np.abs(exact))
    chk(err < 1e-4, "2 Omega recupera la derivada analitica",
        "error relativo maximo %.2e con dt = %.3f s" % (err, dt))
    # y cada canal por separado, que es lo que se va a leer
    e1 = np.max(np.abs(d["FLUJO"] - ((Pa / Ta) * dQ)[1:-1])) / np.max(np.abs(exact))
    e3 = np.max(np.abs(d["LIBRO"] + ((Qa * Pa / Ta ** 2) * dT)[1:-1])) / np.max(np.abs(exact))
    chk(e1 < 1e-3 and e3 < 1e-3, "2b cada canal por separado, no solo la suma",
        "FLUJO %.2e   LIBRO %.2e" % (e1, e3))

    # --- 3. Canales degenerados: si la variable no se mueve, su canal es 0 ---
    d = descomponer(t, Q, P, np.full(n, 2.5))
    chk(np.all(d["LIBRO"] == 0.0), "3 tau_0 constante -> canal LIBRO exactamente 0",
        "max|LIBRO| = %.2e" % np.max(np.abs(d["LIBRO"])))
    cte = descomponer(t, np.full(n, 3.0), np.full(n, 100.0), np.full(n, 2.0))
    chk(np.max(np.abs(cte["omega"])) == 0.0, "3b phi constante -> Omega exactamente 0")

    # --- 4. La forma invertida da el canal LIBRO con el signo opuesto --------
    Tm = 1.0 + 0.001 * np.arange(n)          # tau_0 monotona creciente
    Qp = 1.0 + np.abs(rng.normal(0, 0.1, n))  # Q > 0 para que el signo sea legible
    a = descomponer(t, Qp, P, Tm, "P/tau")
    b = descomponer(t, Qp, P, Tm, "tau/P")
    chk(np.all(a["LIBRO"] < 0) and np.all(b["LIBRO"] > 0),
        "4 las dos formas dan el canal LIBRO con signo opuesto",
        "P/tau %+.3e   tau/P %+.3e" % (np.median(a["LIBRO"]), np.median(b["LIBRO"])))

    # --- 5. El reparto de varianza suma 1 -----------------------------------
    rep = reparto_varianza(descomponer(t, Q, P, T))
    su = sum(rep[c] for c in CANALES)
    chk(abs(su - 1.0) < 1e-10, "5 el reparto de Var(Omega) suma 1",
        "suma %.12f" % su)

    # --- 6. `primer_indice_sobre` contra una referencia por fuerza bruta -----
    p = rng.random(1200) * 10.0
    kq = rng.integers(0, 1150, 400)
    lv = rng.random(400) * 10.0
    got = primer_indice_sobre(p, kq, lv, filas_max=90)
    ref = np.full(kq.size, -1, np.int64)
    for i in range(kq.size):
        hi = min(kq[i] + 90, p.size - 1)
        w = np.flatnonzero(p[kq[i] + 1:hi + 1] >= lv[i])
        if w.size:
            ref[i] = kq[i] + 1 + w[0]
    chk(np.array_equal(got, ref), "6 busqueda por bloques == fuerza bruta O(n^2)",
        "%d consultas, tope de %d filas respetado" % (kq.size, 90))
    # con B distinto tiene que dar lo mismo: el bloque es detalle de implementacion
    chk(np.array_equal(primer_indice_sobre(p, kq, lv, 90, B=7), ref),
        "6b el resultado no depende del tamano de bloque")

    # --- 7. tau_recup recupera una constante de relajacion CONOCIDA ---------
    dtb, tau_v = 0.02, 5.0
    tb2 = np.arange(0, 4000, dtb)
    D_eq, caida = 10.0, 0.80
    prof = np.full(tb2.size, D_eq)
    golpes = np.arange(200, tb2.size - 2000, 4000)
    for g in golpes:                      # cae y relaja con constante `tau_v`
        rel = np.arange(tb2.size - g) * dtb
        prof[g:] = D_eq - (D_eq * caida) * np.exp(-rel / tau_v)
    te, tau, cens = eventos_recuperacion(tb2, prof, 0.5)
    err = abs(np.median(tau) - tau_v) / tau_v if tau.size else np.inf
    chk(tau.size == golpes.size and err < 0.02,
        "7 tau_recup recupera la constante de relajacion conocida",
        "%d eventos, tau = %.4f s (verdad %.1f, error %.2f %%)"
        % (tau.size, np.median(tau) if tau.size else np.nan, tau_v, 100 * err))
    # censura: si nunca se recupera, el evento NO se cuenta como rapido
    prof2 = np.r_[np.full(50, 10.0), np.full(50, 1.0)]
    te2, tau2, cens2 = eventos_recuperacion(np.arange(100.0), prof2, 0.5)
    chk(te2.size == 0 and cens2 == 1, "7b evento sin recuperar -> CENSURADO, no rapido",
        "eventos %d, censurados %d" % (te2.size, cens2))

    # --- 8. tau_agot recupera una razon profundidad/caudal conocida ---------
    dtg, nb, t0 = 10.0, 30, 0.0
    acc = _acc_vacio(nb)
    bt = np.arange(0, nb * dtg, 0.5)
    prof_v, caudal_v = 8.0, 2.0               # BTC y BTC/s -> tau = 4 s
    acumular_libro(bt, np.full(bt.size, 100.0), np.full(bt.size, prof_v),
                   np.full(bt.size, 100.1), np.full(bt.size, prof_v),
                   t0, dtg, nb, acc, None)
    tt2 = np.arange(0, nb * dtg, 0.25)
    acumular_trades(tt2, np.full(tt2.size, 100.0),
                    np.full(tt2.size, caudal_v * 0.25),
                    rng.random(tt2.size) < 0.5, t0, dtg, nb, acc)
    g = cerrar_rejilla(acc, t0, dtg, nb)
    ta = np.nanmedian(g["tau_agot"][g["valida"]])
    chk(abs(ta - prof_v / caudal_v) < 1e-9, "8 tau_agot = profundidad / caudal",
        "medido %.6f s, verdad %.6f s" % (ta, prof_v / caudal_v))

    # --- 9. La rejilla suma el flujo firmado que se le metio -----------------
    nb2 = 5
    acc2 = _acc_vacio(nb2)
    tq = np.array([1.0, 2.0, 11.0, 12.0, 13.0])
    qq = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    mk = np.array([False, True, False, True, False])     # +1 -1 +1 -1 +1
    acumular_trades(tq, np.full(5, 100.0), qq, mk, 0.0, 10.0, nb2, acc2)
    g2 = cerrar_rejilla(acc2, 0.0, 10.0, nb2)
    chk(abs(g2["q_neto"][0] - (1 - 2)) < 1e-12
        and abs(g2["q_neto"][1] - (4 - 8 + 16)) < 1e-12
        and abs(g2["q_tot"][1] - 28.0) < 1e-12,
        "9 flujo neto firmado por casilla",
        "casilla0 %+.1f  casilla1 %+.1f" % (g2["q_neto"][0], g2["q_neto"][1]))

    # --- 10. No se deriva a traves de un hueco ------------------------------
    nb3 = 60
    r3 = {"t": np.arange(nb3) * 10.0, "q_neto": rng.normal(0, 1, nb3),
          "precio": np.full(nb3, 60000.0), "tau_agot": np.full(nb3, 3.0),
          "valida": np.ones(nb3, bool)}
    r3["valida"][25:35] = False                 # hueco de 100 s
    s3 = serie_omega(r3)
    chk(s3["n_rachas"] == 2 and not np.any((s3["t"] > 240) & (s3["t"] < 360)),
        "10 un hueco parte la serie en vez de derivar a traves de el",
        "rachas %d, casillas %d de %d" % (s3["n_rachas"], s3["t"].size, nb3))

    # --- 11. El nulo por rotacion destruye la asociacion y conserva la marginal
    x = np.cumsum(rng.normal(0, 1, 600))          # persistente, como `nu`
    y = x + rng.normal(0, 0.3, 600)
    rho = corr(x, y)
    su = suelo_rotacion(x, y, n=60)
    yr = np.roll(y, 137)
    chk(rho > 0.9 and su < 0.9 and np.allclose(np.sort(yr), np.sort(y)),
        "11 la rotacion rompe el emparejamiento y conserva la marginal exacta",
        "rho %.3f -> suelo q95 %.3f" % (rho, su))

    # --- 12. Spearman con empates, contra un caso de verdad conocida --------
    a1 = np.array([1.0, 2, 3, 4, 5]); b1 = np.array([5.0, 4, 3, 2, 1])
    e = np.array([1.0, 1, 2, 2, 3]); f = np.array([1.0, 1, 2, 2, 3])
    chk(abs(corr(a1, b1) + 1.0) < 1e-12 and abs(corr(e, f) - 1.0) < 1e-12,
        "12 Spearman: -1 exacto en orden inverso, +1 con empates promediados",
        "%.6f / %.6f" % (corr(a1, b1), corr(e, f)))

    # --- 13. Los censurados son un ESCALAR, no una serie por casilla --------
    nb4 = 40
    acc4 = _acc_vacio(nb4)
    tb4 = np.arange(0, nb4 * 10.0, 0.5)
    pf = np.full(tb4.size, 10.0)
    pf[100::200] = 0.5                    # caidas sin recuperacion posible al final
    acumular_libro(tb4, np.full(tb4.size, 100.0), pf, np.full(tb4.size, 100.1), pf,
                   0.0, 10.0, nb4, acc4, None)
    g4 = cerrar_rejilla(acc4, 0.0, 10.0, nb4)
    tot_ev = float(g4["rec_n_0.50"].sum()) + float(g4["rec_cens_0.50"].sum())
    chk(g4["rec_cens_0.50"].size == 1 and tot_ev <= tb4.size,
        "13 censurados = escalar acotado por el numero de filas",
        "eventos+censurados %.0f sobre %d filas" % (tot_ev, tb4.size))

    log("")
    log("  %d / %d" % (n_ok, n_tot))
    return 0 if n_ok == n_tot else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", choices=("serie", "dia", "relacion"))
    ap.add_argument("--fuente", default="estacional", choices=("estacional", "v33"))
    ap.add_argument("--tau", default="tau_agot",
                    help="tau_agot | tau_upd | tau_recup_0.50 ...")
    ap.add_argument("--forma", default="P/tau", choices=("P/tau", "tau/P"))
    ap.add_argument("--dt", type=float, default=DT_REJILLA)
    ap.add_argument("--bloque", type=float, default=BLOQUE_S)
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    if a.etapa == "serie":
        return etapa_serie(a)
    if a.etapa == "dia":
        return etapa_dia(a)
    if a.etapa == "relacion":
        return etapa_relacion(a)
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
