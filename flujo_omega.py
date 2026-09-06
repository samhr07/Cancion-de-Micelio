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
# ---------------------------------------------------------------------------
# RUTAS. No estan cableadas: las capturas del operador viven en una USB (D:) y
# el modulo tiene que poder leerlas desde donde esten sin editar codigo.
#
#   --datos=RUTA   directorio de `captura_estacional`   (env MICELIO_DATOS)
#   --v33=RUTA     directorio de `captura_v33`          (env MICELIO_V33)
#   --salida=RUTA  donde se escriben npz / csv / json   (env MICELIO_SALIDA)
#
# En Windows valen tanto "D:/Maxwell/telemetria/estacional" como
# "D:\Maxwell\telemetria\estacional"; `os.path.join` no toca lo que ya es
# absoluto. La SALIDA se separa de los DATOS a proposito: escribir la rejilla
# cacheada en la USB la ata al medio extraible, y esa rejilla es justo lo que
# conviene tener en el disco -- son decenas de MB y de ella comen todas las
# etapas.
# ---------------------------------------------------------------------------

DIR_SALIDA = os.environ.get("MICELIO_SALIDA", "telemetria")
DIR_V33 = os.environ.get("MICELIO_V33", "telemetria/captura_v33")


def ruta_cache(fuente):
    return os.path.join(DIR_SALIDA, "rejilla_omega_%s.npz" % fuente)


def ruta_salida(nombre):
    return os.path.join(DIR_SALIDA, nombre)


def configurar_rutas(args):
    """Aplica --datos / --v33 / --salida. Devuelve lo que quedo resuelto."""
    global DIR_SALIDA, DIR_V33
    import curvas_estacional as C
    if getattr(args, "datos", None):
        C.DIR = args.datos
    elif os.environ.get("MICELIO_DATOS"):
        C.DIR = os.environ["MICELIO_DATOS"]
    if getattr(args, "v33", None):
        DIR_V33 = args.v33
    if getattr(args, "salida", None):
        DIR_SALIDA = args.salida
    os.makedirs(DIR_SALIDA, exist_ok=True)
    return {"datos": C.DIR, "v33": DIR_V33, "salida": DIR_SALIDA}


def _anadir_rutas(ap):
    ap.add_argument("--datos", help="directorio de captura_estacional"
                                    " (env MICELIO_DATOS)")
    ap.add_argument("--v33", help="directorio de captura_v33 (env MICELIO_V33)")
    ap.add_argument("--salida", help="donde escribir npz/csv/json"
                                     " (env MICELIO_SALIDA)")
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


def piso_parpadeo(prof, vol):
    """Nivel de PARPADEO del libro, medido: |dprof|/prof donde NO hubo flujo.

    Es el umbral que separa "la profundidad bajo porque se la comieron" de "la
    profundidad bajo porque el creador de mercado movio su cotizacion". Se mide
    sobre los intervalos con `vol == 0`, donde por construccion ningun cambio de
    profundidad es atribuible a transacciones. **No es un parametro elegido**:
    sale del propio dato, y se reporta.
    """
    prof = np.asarray(prof, float); vol = np.asarray(vol, float)
    if prof.size < 3:
        return np.nan
    d = np.abs(prof[1:] - prof[:-1]) / np.maximum(prof[:-1], 1e-12)
    sin = (vol[1:] <= 0) & np.isfinite(d)
    if sin.sum() < 20:
        return np.nan
    return float(np.median(d[sin]))


def eventos_recuperacion_local(t, prof, vol, censura=TAU_CENSURA_S, piso=None):
    """Recuperacion con `theta` MEDIDA, no elegida. Retroalimenta `tau_agot`.

    ⚠ ESTE ES EL ESTIMADOR QUE PIDIO EL OPERADOR, y su objecion era correcta:
    fijar `theta` en una fraccion constante es vacio, porque no hay ninguna
    razon para que el libro se agote un 30 % o un 90 % en un intervalo dado.
    Aqui `theta` es **lo que el flujo llego a agotar**:

        theta_k = volumen negociado en (t[k-1], t[k]] / profundidad en t[k-1]

    que es exactamente `dt_k / tau_agot_k`, o sea la fraccion de la cola visible
    que el caudal consumio en ese intervalo. Un evento es que la profundidad
    caiga **al menos lo que el flujo se llevo**:

        prof[k] <= (1 - theta_k) * prof[k-1]

    o sea, que la reposicion no compensara al consumo. Con eso `tau_recup` y
    `tau_agot` dejan de ser dos medidas independientes del libro y pasan a ser
    las dos mitades de un mismo ciclo: cuanto tarda en vaciarse y cuanto en
    volver.

    EL PISO NO ES UN PARAMETRO LIBRE. Sin piso, un intervalo con `vol = 0` da
    `theta = 0` y entonces CUALQUIER bajada cuenta como evento -- incluido el
    parpadeo de cotizacion, que es justo lo que hundio la primera version del
    estimador. El piso es el parpadeo MEDIDO (`piso_parpadeo`): la mediana de
    |dprof|/prof sobre los intervalos SIN transacciones, donde ningun cambio es
    atribuible al flujo. Se mide, se reporta y no se elige.
    """
    t = np.asarray(t, float); prof = np.asarray(prof, float)
    vol = np.asarray(vol, float)
    if t.size < 3:
        return np.empty(0), np.empty(0), np.empty(0), 0, np.nan
    if piso is None:
        piso = piso_parpadeo(prof, vol)
    if not np.isfinite(piso):
        piso = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        theta = vol[1:] / np.maximum(prof[:-1], 1e-12)
    theta = np.clip(theta, 0.0, 0.99)
    k = np.flatnonzero((theta > piso) & (prof[1:] <= (1.0 - theta) * prof[:-1])) + 1
    if k.size == 0:
        return np.empty(0), np.empty(0), np.empty(0), 0, piso
    objetivo = 0.5 * (prof[k - 1] + prof[k])
    dtm = float(np.median(np.diff(t))) if t.size > 1 else 1.0
    filas_max = int(min(max(np.ceil(censura / max(dtm, 1e-6)), 8), 20000))
    j = primer_indice_sobre(prof, k, objetivo, filas_max)
    vivo = j >= 0
    tau = np.full(k.size, np.nan)
    tau[vivo] = (t[j[vivo]] - t[k[vivo]]) / np.log(2.0)
    ok = vivo & np.isfinite(tau) & (tau > 0) & (tau <= censura)
    return t[k[ok]], tau[ok], theta[k - 1][ok], int((~ok).sum()), piso


# ===========================================================================
# 3. Rejilla uniforme: una casilla de `dt` segundos
# ===========================================================================

def _bins(t, t0, dt, nb):
    k = np.floor((np.asarray(t, float) - t0) / dt).astype(np.int64)
    ok = (k >= 0) & (k < nb)
    return k, ok


def _acc_vacio(nb):
    c = ("q_neto", "eps_neto", "q_tot", "n_tx", "p_sum", "prof_sum", "n_libro",
         "rv", "mid_ult", "p_ult")
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
    # theta MEDIDA (retroalimentada de tau_agot), ver `eventos_recuperacion_local`
    a["loc_log"] = np.zeros(nb)      # sum log(tau_recup) por casilla
    a["loc_n"] = np.zeros(nb)
    a["loc_theta"] = np.zeros(nb)    # sum de la theta medida, para reportarla
    a["cens_loc"] = 0
    a["piso"] = []                   # parpadeo medido, una entrada por parte
    return a


def _ultimo_por_casilla(k, v, dest):
    """Escribe en `dest[k]` el ULTIMO `v` de cada casilla (cierre)."""
    if k.size == 0:
        return
    o = np.argsort(k, kind="stable")
    kk, vv = k[o], v[o]
    ult = np.r_[np.flatnonzero(np.diff(kk)), kk.size - 1]
    dest[kk[ult]] = vv[ult]


def acumular_libro(bt, bb, bB, ba, bA, t0, dt, nb, acc, prev=None, vol_fn=None):
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
    # theta MEDIDA: necesita el volumen negociado en cada intervalo del libro
    if vol_fn is not None:
        vol = vol_fn(bt)
        te, tau, thm, cens, piso = eventos_recuperacion_local(bt, prof, vol)
        if np.isfinite(piso):
            acc["piso"].append(float(piso))
        acc["cens_loc"] += cens
        if te.size:
            ke, oke = _bins(te, t0, dt, nb)
            if oke.any():
                acc["loc_log"] += np.bincount(ke[oke], weights=np.log(tau[oke]),
                                              minlength=nb)
                acc["loc_n"] += np.bincount(ke[oke], minlength=nb)
                acc["loc_theta"] += np.bincount(ke[oke], weights=thm[oke], minlength=nb)
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
    # Desequilibrio por CONTEO de transacciones (compras menos ventas), que es
    # otra cantidad que la de volumen: Jones, Kaul & Lipson (1994) y el propio
    # `cont2014.py` de este proyecto midieron que el NUMERO de operaciones lleva
    # informacion que su tamano no lleva -- alli el desequilibrio de
    # transacciones salio significativo en el 91 % de submuestras contra el 31 %
    # del articulo original.
    acc["eps_neto"] += np.bincount(k[ok], weights=eps[ok], minlength=nb)
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
         "q_neto": acc["q_neto"], "eps_neto": acc["eps_neto"],
         "q_tot": q_tot, "n_tx": n_tx,
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
    n = acc["loc_n"]
    with np.errstate(divide="ignore", invalid="ignore"):
        r["tau_recup_loc"] = np.where(n > 0, np.exp(acc["loc_log"] / np.maximum(n, 1)),
                                      np.nan)
        r["theta_medida"] = np.where(n > 0, acc["loc_theta"] / np.maximum(n, 1), np.nan)
        # RESILIENCIA: cuanto tarda en vaciarse contra cuanto tarda en volver.
        # >> 1 el libro repone mucho mas rapido de lo que el flujo lo consume;
        # << 1 el flujo lo vacia mas rapido de lo que el libro repone.
        r["resiliencia"] = tau_agot / np.where(r["tau_recup_loc"] > 0,
                                               r["tau_recup_loc"], np.nan)
    r["rec_loc_n"] = n
    r["rec_loc_cens"] = np.array([float(acc["cens_loc"])])
    r["piso_parpadeo"] = np.array([float(np.median(acc["piso"])) if acc["piso"]
                                   else np.nan])
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


def volumen_por_intervalo(tt, q):
    """Devuelve `vol_fn(bt)`: volumen negociado en cada intervalo `(bt[k-1], bt[k]]`.

    `vol[0] = 0` por definicion (no hay intervalo anterior al primer snapshot).
    Es una suma acumulada mas dos busquedas binarias, o sea O(n log n) sin
    materializar nada del tamano del producto.
    """
    tt = np.asarray(tt, float); q = np.asarray(q, float)
    o = np.argsort(tt, kind="stable")
    tt, q = tt[o], q[o]
    cum = np.concatenate([[0.0], np.cumsum(q)])

    def _fn(bt):
        idx = np.searchsorted(tt, np.asarray(bt, float), side="right")
        v = np.empty(idx.size)
        v[0] = 0.0
        v[1:] = cum[idx[1:]] - cum[idx[:-1]]
        return np.maximum(v, 0.0)
    return _fn


class _TradesPorVentana:
    """Carga solo las partes de transacciones que solapan con la ventana pedida.

    ⚠ MEMORIA. Las 32 M de transacciones del estacional son ~512 MB si se
    concatenan, y esta maquina tiene ~2 GB libres. Se cachea UNA sola ventana a
    la vez: como las partes de libro y de transacciones van las dos en orden
    temporal, cada parte de transacciones se lee un par de veces y nunca hay mas
    de una ventana viva.
    """

    def __init__(self, indice):
        self.idx = sorted(indice, key=lambda z: z[1])
        self.rango = (np.inf, -np.inf)
        self.fn = None

    def para(self, t0, t1):
        if self.fn is not None and t0 >= self.rango[0] and t1 <= self.rango[1]:
            return self.fn
        import pyarrow.parquet as pq
        m = 60.0                        # margen, para no recargar por un borde
        sel = [x for x in self.idx if not (x[2] < t0 - m or x[1] > t1 + m)]
        ts, qs = [], []
        for f, _, _, _ in sel:
            try:
                tb = pq.read_table(f, columns=["t", "cant"])
            except Exception:
                continue
            a_ = tb["t"].to_numpy().astype(float)
            b_ = tb["cant"].to_numpy().astype(float)
            del tb
            k = (a_ >= t0 - m) & (a_ <= t1 + m) & (b_ > 0)
            if k.any():
                ts.append(a_[k]); qs.append(b_[k])
        if not ts:
            self.rango = (t0 - m, t1 + m)
            self.fn = lambda bt: np.zeros(np.asarray(bt).size)
            return self.fn
        self.rango = (t0 - m, t1 + m)
        self.fn = volumen_por_intervalo(np.concatenate(ts), np.concatenate(qs))
        return self.fn


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
    trades = _TradesPorVentana(C._indice("trades_"))
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
        vf = trades.para(float(bt[o][0]) - 5.0, float(bt[o][-1]) + 5.0)
        prev = acumular_libro(bt[o], bb[o], bB[o], ba[o], bA[o], t0, dt, nb, acc,
                              prev, vol_fn=vf)
        del bt, bb, bB, ba, bA, o, m
        # cada 50 partes no vale para todos los casos: con 12 partes (una
        # captura de 12 dias descargada) no se imprimiria NUNCA y la corrida
        # parece colgada. El paso se adapta al numero de partes.
        paso = max(1, len(partes) // 20)
        if (i + 1) % paso == 0 or i + 1 == len(partes):
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
    d = cargar_larga(DIR_V33)
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
    acumular_libro(bt, bb, bB, ba, bA, t0, dt, nb, acc, None,
                   vol_fn=volumen_por_intervalo(tt[okt], q[okt]))
    acumular_trades(tt[okt], pr[okt], q[okt], mk[okt], t0, dt, nb, acc)
    return cerrar_rejilla(acc, t0, dt, nb)


def procedencia(dir_datos, dt, t0, t1, extra=None):
    """Sello de PROCEDENCIA de la rejilla: de que foto de la captura salio.

    ⚠ EXISTE PORQUE LA CAPTURA SIGUE ESCRIBIENDO. Aviso del operador
    (2026-09-06): la rejilla que se genere hoy es una foto de hoy, su ultimo dia
    puede estar PARCIAL, y regenerarla manana NO da lo mismo. Una cifra que se
    cite tiene que poder decir de que version salio, asi que la rejilla se lleva
    dentro el instante de generacion, la ventana, el recuento de partes y un
    hash de la lista de ficheros fuente (nombre, tamano, filas). Regenerar y
    comparar el hash dice en un segundo si es la misma foto o no.
    """
    import hashlib
    import datetime as _dt
    meta = {"generado_utc": _dt.datetime.now(_dt.timezone.utc)
            .strftime("%Y-%m-%dT%H:%M:%SZ"),
            "dir_datos": os.path.abspath(dir_datos) if dir_datos else None,
            "dt": dt, "t_ini": t0, "t_fin": t1,
            "modulo": "flujo_omega"}
    try:
        import curvas_estacional as C
        firma = []
        for pref in ("trades_", "libro_"):
            idx = C._indice(pref)
            meta["n_partes_" + pref.strip("_")] = len(idx)
            for f, a_, b_, n_ in idx:
                firma.append("%s|%d|%.3f|%.3f" % (os.path.basename(f), n_, a_, b_))
        meta["hash_fuente"] = hashlib.sha256(
            "\n".join(sorted(firma)).encode()).hexdigest()[:16]
    except Exception as e:
        meta["hash_fuente"] = None
        meta["nota"] = "sin indice de partes: %s" % e
    if extra:
        meta.update(extra)
    return meta


def cargar_rejilla(fuente):
    r = ruta_cache(fuente)
    if not os.path.exists(r):
        raise SystemExit(
            "falta %s\n"
            "  corre primero:  --etapa=serie --fuente=%s\n"
            "  y si las capturas no estan en ./telemetria, pasa --datos=RUTA"
            % (r, fuente))
    d = np.load(r, allow_pickle=False)
    out = {}
    for k in d.files:
        v = d[k]
        # el sello de procedencia es texto, no un escalar numerico
        out[k] = str(v) if v.dtype.kind in "USO" else (float(v) if v.ndim == 0 else v)
    return out


# ===========================================================================
# 4.bis  Calendario de NUEVA YORK -- el ciclo diurno que importa es el humano
# ===========================================================================
#
# ⚠ POR QUE NY Y NO UTC. El ciclo de 24 h que este proyecto midio el 2026-08-23
# tiene su pico en UTC 13-15 y su valle en UTC 4, y ese pico ES la apertura de
# Nueva York. Agrupar por hora UTC funciona por casualidad -- porque el offset
# es constante dentro de cada regimen de horario -- pero se rompe dos veces al
# ano, en los cambios de horario de verano, y ademas mezcla el sabado y el
# domingo de NY con el lunes de UTC. Fijar el reloj en NY convierte el ciclo de
# una PRESUNCION sobre el dato en una variable exogena con causa conocida.
#
# Se implementa la regla de EE.UU. a mano y no con `zoneinfo` a proposito: en
# Windows `zoneinfo` necesita el paquete `tzdata` instalado aparte, y este
# modulo tiene que correr en la maquina de la captura sin anadir dependencias.
# La regla vigente desde 2007 es: EDT (UTC-4) desde el 2do domingo de marzo a
# las 02:00 locales hasta el 1er domingo de noviembre a las 02:00 locales; EST
# (UTC-5) el resto.

# 1970-01-01 fue JUEVES, y con 0 = lunes eso es 3, no 4. La primera version
# puso 4 y el control lo caza al instante: 2026-03-08 es domingo y salia habil.
# Quinta vez que este proyecto se juega algo en un offset o un signo (tras el
# 2pi, el factor 125, la convencion de `eps` y la fase de `atan2`).
EPOCH_DOW = 3


def _domingo_n(ano, mes, n):
    """Instante UTC (en dias desde epoch) del n-esimo domingo de `ano`/`mes`."""
    import datetime as _dt
    d = _dt.date(ano, mes, 1)
    # weekday(): 0 = lunes ... 6 = domingo
    primer = 1 + (6 - d.weekday()) % 7
    return _dt.date(ano, mes, primer + 7 * (n - 1))


def offset_ny(t):
    """Offset de Nueva York en segundos (-5 h o -4 h) para cada instante UTC."""
    import datetime as _dt
    t = np.asarray(t, float)
    off = np.full(t.size, -5 * 3600.0)
    anos = np.unique((t // 31556952.0).astype(np.int64) + 1970)
    for a0 in range(int(anos.min()) - 1, int(anos.max()) + 2):
        try:
            ini = _domingo_n(a0, 3, 2)      # 2do domingo de marzo, 02:00 EST = 07:00 UTC
            fin = _domingo_n(a0, 11, 1)     # 1er domingo de noviembre, 02:00 EDT = 06:00 UTC
        except ValueError:
            continue
        t0 = (ini - _dt.date(1970, 1, 1)).days * 86400.0 + 7 * 3600.0
        t1 = (fin - _dt.date(1970, 1, 1)).days * 86400.0 + 6 * 3600.0
        off[(t >= t0) & (t < t1)] = -4 * 3600.0
    return off


def calendario_ny(t):
    """dia local, hora local, dia de la semana y bandera de fin de semana."""
    t = np.asarray(t, float)
    loc = t + offset_ny(t)
    dia = np.floor(loc / 86400.0).astype(np.int64)
    hora = (loc - dia * 86400.0) / 3600.0
    dow = (dia + EPOCH_DOW) % 7                 # 0 = lunes ... 6 = domingo
    return {"dia": dia, "hora": hora, "dow": dow, "finde": dow >= 5, "local": loc}


def _fecha_ny(dia):
    import datetime as _dt
    d = _dt.date(1970, 1, 1) + _dt.timedelta(days=int(dia))
    return d.strftime("%Y-%m-%d %a")


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


def por_bloque(r, s, bloque=BLOQUE_S, zona="ny"):
    """Una fila por (dia, bloque). Es la tabla que el operador pidio.

    ⚠ El dia y la hora son de NUEVA YORK por omision, no UTC, y cada fila lleva
    su bandera `finde`. El ciclo diurno de este mercado tiene su pico en la
    apertura de NY (medido el 2026-08-23: `nu` recorre 3.9x entre UTC 13-15 y
    UTC 4), asi que anclarlo al reloj de NY lo convierte en una variable exogena
    con causa conocida en vez de una presuncion sobre el dato. `--zona=utc`
    vuelve al comportamiento anterior.
    """
    if s["t"].size == 0:
        return []
    if zona == "ny":
        cal = calendario_ny(s["t"])
        dia, seg, finde = cal["dia"], cal["hora"] * 3600.0, cal["finde"]
    else:
        dia, seg = _dia_hora(s["t"])
        finde = np.zeros(dia.size, bool)
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

        def med(col):
            x = r[col][j] if col in r else np.full(j.size, np.nan)
            x = x[np.isfinite(x)]
            return float(np.median(x)) if x.size else np.nan

        f = {"dia": int(dia[sl[0]]),
             "fecha": (_fecha_ny(dia[sl[0]]) if zona == "ny"
                       else _fecha(dia[sl[0]])),
             "finde": int(bool(finde[sl[0]])),
             "bloque": int(blq[sl[0]]), "n": int(sl.size),
             "t0": float(s["t"][sl].min()),
             # los tres Omega del bloque, y son cantidades DISTINTAS:
             "omega_neta": float((ph[-1] - ph[0])
                                 / max(s["t"][sl][-1] - s["t"][sl][0], 1e-9)),
             "omega_rms": float(np.sqrt(np.mean(om ** 2))),
             "omega_med": float(np.median(om)),
             "phi_med": float(np.median(ph)), "phi_sd": float(np.std(ph)),
             # variables de estado del bloque
             "q_neto": float(np.sum(r["q_neto"][j])),
             "q_tot": float(np.sum(r["q_tot"][j])),
             "precio": med("precio"),
             "tau0": med(s["tau"]),
             "tau_agot": med("tau_agot"),
             "tau_recup_loc": med("tau_recup_loc"),
             "theta_medida": med("theta_medida"),
             "resiliencia": med("resiliencia"),
             "prof": med("prof"), "nu": med("nu"),
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
    titulo("tau_0 -- CONCORDANCIA ENTRE LOS ESTIMADORES")
    v = r["valida"]
    log("  casillas validas: %d de %d (%.1f %%)" % (v.sum(), v.size, 100 * v.mean()))
    log("")
    log("  %-24s %10s %10s %10s %10s" % ("estimador", "p10", "MEDIANA", "p90", "n"))
    cols = [("tau_agot [s]", "tau_agot"), ("tau_upd  [s]", "tau_upd"),
            ("tau_recup_LOCAL", "tau_recup_loc")]
    cols += [("tau_recup(%.2f) fijo" % th, "tau_recup_%.2f" % th) for th in THETAS]
    for nom, c in cols:
        x = r[c][v] if c in r else np.zeros(0)
        x = x[np.isfinite(x) & (x > 0)]
        if x.size < 10:
            log("  %-24s %10s %10s %10s %10d" % (nom, "-", "-", "-", x.size))
            continue
        log("  %-24s %10.3f %10.3f %10.3f %10d"
            % (nom, np.percentile(x, 10), np.median(x), np.percentile(x, 90), x.size))
    log("")
    log("  correlacion de Spearman contra `tau_agot` (casillas validas):")
    for nom, c in cols[1:]:
        if c in r:
            log("    tau_agot vs %-20s  rho = %+.4f" % (nom, corr(r["tau_agot"][v], r[c][v])))
    log("")
    log("  --- theta MEDIDA (retroalimentada de tau_agot), no elegida ---")
    if "piso_parpadeo" in r:
        log("    piso de parpadeo medido (|dprof|/prof sin flujo): %.4f"
            % float(np.ravel(r["piso_parpadeo"])[0]))
    if "theta_medida" in r:
        x = r["theta_medida"][v]
        x = x[np.isfinite(x)]
        if x.size:
            log("    theta medida en los eventos: p10 %.4f  MED %.4f  p90 %.4f"
                % tuple(np.percentile(x, [10, 50, 90])))
    if "rec_loc_n" in r:
        ne = float(r["rec_loc_n"][v].sum()); ce = float(np.ravel(r["rec_loc_cens"])[0])
        log("    eventos con recuperacion %.0f, censurados %.0f (%.1f %%),"
            " casillas con dato %d de %d"
            % (ne, ce, 100 * ce / max(ne + ce, 1),
               int((r["rec_loc_n"][v] > 0).sum()), int(v.sum())))
    log("")
    log("  --- RESILIENCIA = tau_agot / tau_recup_LOCAL  (adimensional) ---")
    log("      >> 1 : el libro repone mucho mas rapido de lo que el flujo lo consume")
    log("      << 1 : el flujo lo vacia mas rapido de lo que el libro repone")
    if "resiliencia" in r:
        x = r["resiliencia"][v]
        x = x[np.isfinite(x) & (x > 0)]
        if x.size:
            log("      p10 %.3f   MEDIANA %.3f   p90 %.3f   n = %d"
                % (np.percentile(x, 10), np.median(x), np.percentile(x, 90), x.size))
            log("      fraccion de casillas con resiliencia < 1 (libro perdiendo): %.2f %%"
                % (100 * np.mean(x < 1.0)))
    log("")
    log("  [!] `tau_agot` es el PRIMARIO: no tiene parametro libre.")
    log("      `tau_recup_LOCAL` tampoco: su theta la fija el propio flujo y su")
    log("      piso lo fija el parpadeo medido. Los `tau_recup(theta)` de umbral")
    log("      FIJO se conservan solo como contraste, y son los que el operador")
    log("      objeto con razon: una fraccion constante no describe nada.")
    log("      `tau_upd` es diagnostico: mide parpadeo de cotizacion, no reposicion.")


def informe_relacion(filas, s, forma):
    titulo("RELACION ENTRE VARIABLES -- nivel de BLOQUE")
    if len(filas) < 20:
        log("  solo %d bloques: por debajo de 20 no se lee nada." % len(filas))
        return
    col = {k: np.array([f.get(k, np.nan) for f in filas], float)
           for k in ("omega_rms", "omega_neta", "phi_med", "phi_sd", "q_neto",
                     "q_tot", "precio", "tau0", "tau_agot", "tau_recup_loc",
                     "theta_medida", "resiliencia", "prof", "nu", "rv_pb",
                     "f_FLUJO", "f_PRECIO", "f_LIBRO")}
    col["abs_q_neto"] = np.abs(col["q_neto"])
    col["abs_omega_neta"] = np.abs(col["omega_neta"])
    dia = np.array([f["dia"] for f in filas])
    finde = np.array([f.get("finde", 0) for f in filas])
    # ⚠ La celda intradia separa HABIL de FIN DE SEMANA. El 2026-08-23 se midio
    # que el fin de semana tiene amplitud Y FASE propias -- su pico llega ~7 h
    # mas tarde -- asi que una sola forma intradia para los siete dias mezcla dos
    # ciclos distintos y no retira ninguno de los dos.
    blq = np.array([f["bloque"] for f in filas]) + 1000 * finde
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
             ("f_LIBRO", "nu"), ("f_LIBRO", "tau0"),
             # el ciclo del libro: agotamiento contra reposicion, y el flujo
             ("tau_agot", "tau_recup_loc"), ("resiliencia", "nu"),
             ("resiliencia", "q_tot"), ("resiliencia", "abs_q_neto"),
             ("resiliencia", "rv_pb"), ("resiliencia", "omega_rms"),
             ("theta_medida", "nu")]
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
    hab = np.array([not f.get("finde", 0) for f in filas])
    if hab.any() and (~hab).any():
        log("  --- HABIL contra FIN DE SEMANA (reloj de NY), sobre el RESIDUO ---")
        log("  %-30s %10s %10s" % ("par", "habil", "finde"))
        for a2, b2 in (("resiliencia", "nu"), ("omega_rms", "nu"),
                       ("tau_agot", "tau_recup_loc"), ("omega_rms", "rv_pb")):
            log("  %-30s %+10.4f %+10.4f"
                % ("%s vs %s" % (a2, b2), corr(res[a2][hab], res[b2][hab]),
                   corr(res[a2][~hab], res[b2][~hab])))
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
                      ("PREDICTIVO    phi_med    vs ret(t+1)", col["phi_med"], ret),
                      ("PREDICTIVO    resiliencia vs |ret(t+1)|", col["resiliencia"],
                       np.abs(ret))):
        rho = corr(a, b)
        su = suelo_rotacion(a, b)
        raz = abs(rho) / su if (np.isfinite(su) and su > 0) else np.nan
        marca = "  <<" if np.isfinite(raz) and raz >= 3.0 else ""
        log("  %-40s %+9.4f %9.4f %7.2f%s" % (nom, rho, su, raz, marca))


# ===========================================================================
# 7. Etapas
# ===========================================================================

def etapa_rutas(args) -> int:
    """Comprobacion previa: que ve el modulo y desde donde. Sin tocar dato.

    Existe porque la primera corrida en la maquina del operador es sobre una
    USB y una ruta de Windows, y un `FileNotFoundError` a los veinte minutos de
    escanear parquet no dice cual de las tres rutas estaba mal.
    """
    import curvas_estacional as C
    titulo("RUTAS -- que ve el modulo y desde donde")
    log("  datos  (captura_estacional) : %s" % os.path.abspath(C.DIR))
    log("  datos  (captura_v33)        : %s" % os.path.abspath(DIR_V33))
    log("  salida (npz / csv / json)   : %s" % os.path.abspath(DIR_SALIDA))
    log("")
    ok = True
    try:
        import pyarrow  # noqa: F401
        import pyarrow.parquet  # noqa: F401
        log("  [OK ] pyarrow disponible")
    except Exception as e:
        ok = False
        log("  [FALLA] pyarrow NO disponible: %s" % e)
        log("          instalalo en el Python que corre esto:")
        log("          C:/Users/Usuario/miniconda3/python.exe -m pip install pyarrow")
    for nom, d in (("captura_estacional", C.DIR), ("captura_v33", DIR_V33)):
        if not os.path.isdir(d):
            log("  [ -- ] %s: no existe %s" % (nom, d))
            continue
        try:
            sub = sorted(os.listdir(d))
        except Exception as e:
            ok = False
            log("  [FALLA] %s: no se puede listar %s (%s)" % (nom, d, e))
            continue
        tr = [x for x in sub if x.startswith("trades_")]
        lb = [x for x in sub if x.startswith("libro_")]
        log("  [OK ] %s: %d entradas (%d trades_, %d libro_)"
            % (nom, len(sub), len(tr), len(lb)))
        if nom == "captura_estacional" and (not tr or not lb):
            log("         ⚠ se esperan subdirectorios `trades_*` y `libro_*` con")
            log("           parquet dentro. Si tu captura tiene otra forma, apunta")
            log("           --datos al directorio que los CONTIENE.")
    try:
        p = os.path.join(DIR_SALIDA, "._prueba_escritura")
        with io.open(p, "w", encoding="ascii") as fh:
            fh.write("ok")
        os.remove(p)
        log("  [OK ] la salida es escribible")
    except Exception as e:
        ok = False
        log("  [FALLA] no se puede escribir en la salida: %s" % e)
    log("")
    log("  rejilla cacheada esperada: %s  (%s)"
        % (ruta_cache("estacional"),
           "EXISTE" if os.path.exists(ruta_cache("estacional")) else "no existe aun"))
    return 0 if ok else 2


def _ventana(args):
    """(t_ini, t_fin) desde --desde / --hasta / --dias. None = sin limite."""
    import datetime as _dt

    def _p(x):
        if not x:
            return None
        for f in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d"):
            try:
                return _dt.datetime.strptime(x, f).replace(
                    tzinfo=_dt.timezone.utc).timestamp()
            except ValueError:
                continue
        raise SystemExit("fecha no reconocida: %r (usa YYYY-MM-DD)" % x)

    t0, t1 = _p(getattr(args, "desde", None)), _p(getattr(args, "hasta", None))
    if getattr(args, "dias", 0):
        a0, b0 = _rango_estacional()
        if t0 is None:
            t0 = a0
        t1 = min(b0 if t1 is None else t1, t0 + args.dias * 86400.0)
    return t0, t1


def etapa_serie(args) -> int:
    titulo("REJILLA -- una casilla de %.0f s sobre `%s`" % (args.dt, args.fuente))
    if args.fuente == "estacional":
        # ⚠ ACOTAR LA VENTANA ES LO PRIMERO QUE CONVIENE HACER en una captura
        # nueva. El estacional son ~223 M de filas de libro, y la primera pasada
        # sobre una USB no es rapida: `--dias=1` da una prueba de humo que
        # termina en minutos y produce una rejilla legible por las tres etapas.
        # Si esa sale bien, se relanza sin `--dias` para la captura entera.
        t0, t1 = _ventana(args)
        if t0 is not None or t1 is not None:
            log("  ventana acotada por bandera")
        r = rejilla_estacional(args.dt, t_ini=t0, t_fin=t1)
    elif args.fuente == "v33":
        r = rejilla_v33(args.dt)
    else:
        raise SystemExit("fuente desconocida: %r" % (args.fuente,))
    import curvas_estacional as C
    meta = procedencia(C.DIR if args.fuente == "estacional" else DIR_V33,
                       args.dt, float(r["t"][0]), float(r["t"][-1]),
                       {"fuente": args.fuente})
    r["_procedencia"] = np.array(json.dumps(meta, ensure_ascii=True))
    np.savez_compressed(ruta_cache(args.fuente), **r)
    log("")
    log("  guardado en %s" % ruta_cache(args.fuente))
    log("  PROCEDENCIA  generada %s   hash de fuente %s"
        % (meta["generado_utc"], meta.get("hash_fuente")))
    log("               partes: %s trades, %s libro"
        % (meta.get("n_partes_trades"), meta.get("n_partes_libro")))
    log("               ⚠ la captura sigue escribiendo: esta rejilla es una FOTO.")
    log("                 Cita la version por su hash; regenerarla no da lo mismo.")
    v = r["valida"]
    log("  casillas %d, validas %d (%.1f %%), transacciones %d"
        % (v.size, v.sum(), 100 * v.mean(), int(r["n_tx"].sum())))
    if v.any():
        log("  nu mediana %.2f tx/s   precio %.1f -- %.1f   profundidad L1 mediana %.3f BTC"
            % (np.median(r["nu"][v]), np.nanmin(r["precio"][v]),
               np.nanmax(r["precio"][v]), np.nanmedian(r["prof"][v])))
    cobertura_por_dia(r, args.zona)
    informe_tau0(r)
    return 0


def cobertura_por_dia(r, zona="ny"):
    """Fraccion de casillas validas por dia. Un dia PARCIAL se ve aqui.

    El ultimo dia de una captura viva casi siempre esta a medias, y una fila con
    el 30 % de cobertura no es comparable con una del 99 %. Se reporta antes que
    nada para que no se lea como si lo fuera.
    """
    titulo("COBERTURA POR DIA -- un dia parcial se ve aqui")
    cal = calendario_ny(r["t"]) if zona == "ny" else None
    dia = (cal["dia"] if cal is not None
           else np.floor(r["t"] / 86400.0).astype(np.int64))
    esperadas = 86400.0 / float(r["dt"])
    log("  %-16s %3s %10s %10s %9s %9s"
        % ("fecha", "fin", "validas", "esperadas", "cobertura", "nu MED"))
    for d in np.unique(dia):
        m = dia == d
        v = r["valida"] & m
        nu = r["nu"][v]
        cob = v.sum() / esperadas
        marca = "   <- PARCIAL" if cob < 0.90 else ""
        log("  %-16s %3s %10d %10.0f %8.1f %% %9.2f%s"
            % (_fecha_ny(d) if cal is not None else _fecha(d),
               "SI" if (cal is not None and cal["finde"][m][0]) else "-",
               int(v.sum()), esperadas, 100 * cob,
               float(np.median(nu)) if nu.size else np.nan, marca))


def etapa_tau(args) -> int:
    """Diagnostico de `tau_recup`: por que es tan rapida, y ESTRATIFICADO.

    ⚠ ESTA ETAPA EXISTE PORQUE LA LECTURA AGRUPADA ERA ENGANOSA. La primera
    corrida sobre dato real reporto `tau_recup` (MED 0.007 s) "pegada a
    `tau_upd`" (MED 0.006 s) y concluyo que el estimador estaba en el suelo de
    resolucion. Objecion del operador: eso se midio promediando los 12 dias, sin
    separar los extremos de la U de volatilidad. Tenia razon, y al estratificar
    la lectura CAMBIA -- ver la tabla que imprime esta etapa.
    """
    r = cargar_rejilla(args.fuente)
    v = r["valida"] & np.isfinite(r["tau_recup_loc"]) & (r["tau_recup_loc"] > 0)
    titulo("tau_recup -- DIAGNOSTICO ESTRATIFICADO")
    if v.sum() < 500:
        log("  casillas insuficientes."); return 2
    ta, tu, tr = r["tau_agot"][v], r["tau_upd"][v], r["tau_recup_loc"][v]
    raz = tr / tu
    log("  casillas con tau_recup: %d de %d validas" % (v.sum(), r["valida"].sum()))
    log("")
    log("  --- LAS DOS MITADES DE LA PREGUNTA ---")
    log("")
    log("  (a) De los que VUELVEN: cuantas actualizaciones de libro tardan")
    log("      tau_recup / tau_upd   (1.0 = una sola actualizacion)")
    log("      p10 %.2f  p25 %.2f  MED %.2f  p75 %.2f  p90 %.2f  p99 %.2f"
        % tuple(np.percentile(raz, [10, 25, 50, 75, 90, 99])))
    log("      fraccion en 2 actualizaciones o menos: %.1f %%" % (100 * np.mean(raz <= 2)))
    ne = float(r["rec_loc_n"][v].sum())
    nc = float(np.ravel(r["rec_loc_cens"])[0])
    log("")
    log("  (b) De los que NO vuelven: censurados a %.0f s" % TAU_CENSURA_S)
    log("      recuperados %.0f   censurados %.0f   -> **%.1f %% no vuelve**"
        % (ne, nc, 100 * nc / max(ne + nc, 1)))
    log("")
    log("  --- ESTRATIFICADO POR LA U DE VOLATILIDAD (quintiles de rv_pb) ---")
    rv = r["rv_pb"][v]
    q = np.percentile(rv, [20, 40, 60, 80])
    est = np.digitize(rv, q)
    nom = ("muy baja", "baja", "media", "alta", "MUY ALTA")
    log("  %-10s %8s %8s %8s %10s %10s %9s %9s"
        % ("estrato", "n", "rv MED", "nu MED", "tau_upd", "tau_rec", "rec/upd",
           "tau_agot"))
    for k in range(5):
        m = est == k
        if m.sum() < 50:
            continue
        log("  %-10s %8d %8.3f %8.2f %10.4f %10.4f %9.2f %9.3f"
            % (nom[k], m.sum(), np.median(rv[m]), np.median(r["nu"][v][m]),
               np.median(tu[m]), np.median(tr[m]), np.median(raz[m]),
               np.median(ta[m])))
    log("")
    log("  --- LA CONCORDANCIA, DENTRO de cada estrato ---")
    log("  (agrupando los 12 dias salia rho = -0.575; ver si eso era el regimen)")
    log("  %-10s %26s %26s" % ("estrato", "rho(tau_agot, tau_recup)",
                               "rho(tau_agot, tau_upd)"))
    for k in range(5):
        m = est == k
        if m.sum() < 50:
            continue
        log("  %-10s %26.4f %26.4f"
            % (nom[k], corr(ta[m], tr[m]), corr(ta[m], tu[m])))
    log("")
    log("  [!] Si la concordancia dentro de estrato es mucho mas debil que")
    log("      agrupada, la anticorrelacion agrupada era CONFUSION POR REGIMEN:")
    log("      `tau_agot` baja con la actividad y `tau_recup` sube con ella, asi")
    log("      que mezclar regimenes la fabrica. No son opuestas: son casi")
    log("      ortogonales dentro de un regimen, o sea que miden cosas")
    log("      DISTINTAS -- que no es lo mismo que medir cosas contradictorias.")
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
    filas = por_bloque(r, s, args.bloque, zona=args.zona)
    titulo("LOS Omega POR DIA -- %d bloques de %.0f min, reloj de %s"
           % (len(filas), args.bloque / 60.0, "NUEVA YORK" if args.zona == "ny" else "UTC"))
    log("  omega_neta = (phi_fin - phi_ini)/T del bloque  [USD/s^2]  -- deriva neta")
    log("  omega_rms  = raiz de la media de Omega^2       [USD/s^2]  -- actividad")
    log("  resil      = tau_agot / tau_recup_LOCAL (adimensional; < 1 = libro perdiendo)")
    log("  f_X        = parte de Var(Omega) del canal X (los tres suman 1)")
    log("")
    dias = sorted(set(f["dia"] for f in filas))
    for d in dias:
        fs = [f for f in filas if f["dia"] == d]
        log("")
        log("  --- %s : %d bloques ---" % (fs[0]["fecha"], len(fs)))
        log("  %3s %6s %12s %12s %10s %9s %8s %8s %8s %7s %7s %7s"
            % ("blq", "n", "omega_neta", "omega_rms", "q_neto", "tau0", "t_rec",
               "resil", "nu", "rv_pb", "fFLU", "fLIB"))
        for f in fs:
            log("  %3d %6d %+12.4e %12.4e %+10.3f %9.4f %8.3f %8.3f %8.2f %7.3f %+7.3f %+7.3f"
                % (f["bloque"], f["n"], f["omega_neta"], f["omega_rms"], f["q_neto"],
                   f["tau0"], f["tau_recup_loc"], f["resiliencia"], f["nu"],
                   f["rv_pb"], f["f_FLUJO"], f["f_LIBRO"]))
    def _nanmed(x):
        x = np.asarray(x, float)
        x = x[np.isfinite(x)]
        return float(np.median(x)) if x.size else np.nan

    titulo("RESUMEN POR DIA -- el cambio de un dia a otro")
    log("  %-16s %3s %5s %13s %10s %9s %8s %8s %8s %7s"
        % ("fecha", "fin", "blq", "omega_rms MED", "q_neto", "tau0", "t_rec",
           "resil", "nu", "rv_pb"))
    for d in dias:
        fs = [f for f in filas if f["dia"] == d]
        g = lambda c: np.array([f[c] for f in fs], float)
        log("  %-16s %3s %5d %13.4e %+10.2f %9.4f %8.3f %8.3f %8.2f %7.3f"
            % (fs[0]["fecha"], "SI" if fs[0]["finde"] else "-", len(fs),
               _nanmed(g("omega_rms")), np.nansum(g("q_neto")), _nanmed(g("tau0")),
               _nanmed(g("tau_recup_loc")), _nanmed(g("resiliencia")),
               _nanmed(g("nu")), _nanmed(g("rv_pb"))))
    if len(dias) > 1:
        rr = np.array([_nanmed([f["omega_rms"] for f in filas if f["dia"] == d])
                       for d in dias], float)
        rr = rr[np.isfinite(rr) & (rr > 0)]
        if rr.size > 1:
            log("")
            log("  recorrido de omega_rms entre dias: %.2fx  (min %.3e, max %.3e)"
                % (np.max(rr) / max(np.min(rr), 1e-30), np.min(rr), np.max(rr)))
    # --- habil contra fin de semana, hora de NUEVA YORK --------------------
    hab = [f for f in filas if not f["finde"]]
    fin_ = [f for f in filas if f["finde"]]
    if hab and fin_:
        titulo("HABIL contra FIN DE SEMANA -- reloj de NUEVA YORK")
        log("  %-10s %6s %13s %10s %9s %8s %8s %8s"
            % ("grupo", "blq", "omega_rms MED", "|q_neto|", "tau0", "t_rec",
               "resil", "nu"))
        for nom, g in (("habil", hab), ("finde", fin_)):
            a_ = lambda c: np.array([f[c] for f in g], float)
            log("  %-10s %6d %13.4e %10.2f %9.4f %8.3f %8.3f %8.2f"
                % (nom, len(g), _nanmed(a_("omega_rms")),
                   _nanmed(np.abs(a_("q_neto"))), _nanmed(a_("tau0")),
                   _nanmed(a_("tau_recup_loc")), _nanmed(a_("resiliencia")),
                   _nanmed(a_("nu"))))
        r1 = _nanmed([f["omega_rms"] for f in hab])
        r2 = _nanmed([f["omega_rms"] for f in fin_])
        if np.isfinite(r1) and np.isfinite(r2) and r2 > 0:
            log("")
            log("  razon habil/finde de omega_rms: %.3fx" % (r1 / r2))
        log("")
        log("  --- perfil horario de omega_rms (hora de NY) ---")
        log("  %5s %14s %14s" % ("hora", "habil", "finde"))
        hs = sorted(set(f["bloque"] for f in filas))
        for h in hs:
            v1 = _nanmed([f["omega_rms"] for f in hab if f["bloque"] == h])
            v2 = _nanmed([f["omega_rms"] for f in fin_ if f["bloque"] == h])
            log("  %5d %14s %14s"
                % (h, "%.4e" % v1 if np.isfinite(v1) else "-",
                   "%.4e" % v2 if np.isfinite(v2) else "-"))

    csv = ruta_salida("omega_bloques_%s.csv" % args.fuente)
    cab = ["fecha", "dia", "finde", "bloque", "n", "t0", "omega_neta", "omega_rms",
           "omega_med", "phi_med", "phi_sd", "q_neto", "q_tot", "precio", "tau0",
           "tau_agot", "tau_recup_loc", "theta_medida", "resiliencia", "prof", "nu",
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
              io.open(ruta_salida("omega_resumen_%s.json" % args.fuente), "w",
                      encoding="ascii"))
    return 0


def etapa_relacion(args) -> int:
    r = cargar_rejilla(args.fuente)
    s = serie_omega(r, tau=args.tau, forma=args.forma)
    if s["t"].size == 0:
        log("  no hay serie. Corre --etapa=serie primero.")
        return 2
    filas = por_bloque(r, s, args.bloque, zona=args.zona)
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
    # `eps_neto` cuenta TRANSACCIONES, no volumen: la casilla 1 tiene +1 -1 +1
    chk(abs(g2["eps_neto"][0] - 0.0) < 1e-12
        and abs(g2["eps_neto"][1] - 1.0) < 1e-12,
        "9b eps_neto es desequilibrio de CONTEO, no de volumen",
        "casilla0 %+.1f (1 compra, 1 venta)  casilla1 %+.1f (2 compras, 1 venta)"
        % (g2["eps_neto"][0], g2["eps_neto"][1]))

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

    # --- 14. Calendario de NY: fronteras de horario de verano EXACTAS -------
    import datetime as _dt

    def _ts(x):
        return _dt.datetime.strptime(x, "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=_dt.timezone.utc).timestamp()
    casos = [("2026-01-15T12:00:00Z", -5), ("2026-07-15T12:00:00Z", -4),
             ("2026-03-08T06:59:00Z", -5), ("2026-03-08T07:01:00Z", -4),
             ("2026-11-01T05:59:00Z", -4), ("2026-11-01T06:01:00Z", -5)]
    ts = np.array([_ts(c[0]) for c in casos])
    got = offset_ny(ts) / 3600.0
    chk(np.array_equal(got, np.array([c[1] for c in casos], float)),
        "14 offset de NY exacto en las dos fronteras de horario de verano",
        "medido %s" % np.array2string(got, precision=0))

    # --- 15. dia de la semana contra la biblioteca estandar -----------------
    ts2 = 1750000000.0 + np.arange(400) * 86400.0
    cal = calendario_ny(ts2)
    ref = np.array([(_dt.date(1970, 1, 1) + _dt.timedelta(days=int(d))).weekday()
                    for d in cal["dia"]])
    chk(np.array_equal(cal["dow"], ref) and np.array_equal(cal["finde"], ref >= 5),
        "15 dia de la semana == datetime.weekday() en 400 dias",
        "y `finde` = sabado o domingo")

    # --- 16. volumen por intervalo: suma exacta -----------------------------
    tt3 = np.array([0.5, 1.5, 2.5, 3.5, 9.5])
    q3 = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
    fn = volumen_por_intervalo(tt3, q3)
    v = fn(np.array([0.0, 2.0, 4.0, 10.0]))
    chk(np.allclose(v, [0.0, 3.0, 12.0, 16.0]),
        "16 volumen por intervalo del libro, suma exacta", "medido %s" % v)

    # --- 17/18. theta MEDIDA y tau_recup con verdad conocida ----------------
    dtb, tau_v, D_eq = 0.02, 5.0, 10.0
    nb5 = 60000
    tb5 = np.arange(nb5) * dtb
    prof5 = np.full(nb5, D_eq)
    vol5 = np.zeros(nb5)
    golpes = np.arange(2000, nb5 - 3000, 4000)
    v_golpe = 4.0                              # theta verdadera = 4/10 = 0.40
    for g in golpes:
        rel = np.arange(nb5 - g) * dtb
        prof5[g:] = D_eq - v_golpe * np.exp(-rel / tau_v)
        # [!] EL NIVEL PREVIO SE FIJA EN D_eq EXACTO. Sin esto queda la cola de
        # la relajacion del golpe anterior (~1e-3) y la condicion de evento
        # -- que la profundidad caiga AL MENOS lo que el flujo se llevo -- se
        # evalua justo en su frontera y falla por un pelo: 1 evento de 14. No es
        # un fallo del estimador, es que el caso de prueba estaba construido en
        # el punto exacto de igualdad. Fijarlo lo deja en el caso MAS ESTRICTO
        # posible (igualdad exacta), que es el que hay que probar.
        prof5[g - 1] = D_eq
        vol5[g] = v_golpe
    piso = piso_parpadeo(prof5, vol5)
    # Sin ruido de cotizacion el piso no es cero: es la COLA DE LA RELAJACION,
    # que tambien cambia la profundidad sin que haya transacciones. Lo que se
    # exige es que sea despreciable frente a cualquier theta real (~0.4).
    chk(np.isfinite(piso) and piso < 1e-4,
        "17 sin parpadeo el piso es despreciable frente a theta",
        "medido %.3e contra theta verdadera 0.40" % piso)
    te, tau, thm, cens, _p = eventos_recuperacion_local(tb5, prof5, vol5)
    err = abs(np.median(tau) - tau_v) / tau_v if tau.size else np.inf
    chk(tau.size == golpes.size and err < 0.02 and abs(np.median(thm) - 0.4) < 1e-9,
        "18 theta MEDIDA y tau_recup local recuperan su verdad conocida",
        "%d eventos, theta = %.4f (verdad 0.40), tau = %.4f s (verdad %.1f)"
        % (tau.size, np.median(thm) if thm.size else np.nan,
           np.median(tau) if tau.size else np.nan, tau_v))

    # --- 19. el piso de parpadeo SI descarta el ruido de cotizacion ---------
    rng2 = np.random.default_rng(5)
    prof6 = prof5 * np.exp(rng2.normal(0, 0.20, nb5))     # parpadeo fuerte
    piso6 = piso_parpadeo(prof6, vol5)
    te6, tau6, thm6, _c6, _p6 = eventos_recuperacion_local(tb5, prof6, vol5)
    chk(piso6 > 0.05 and te6.size <= golpes.size,
        "19 con parpadeo, el piso sube y NO se inventan eventos",
        "piso %.4f, eventos %d contra %d golpes reales"
        % (piso6, te6.size, golpes.size))

    # --- 20. resiliencia = tau_agot / tau_recup, con las dos conocidas ------
    nb7, dt7 = 20, 10.0
    acc7 = _acc_vacio(nb7)
    bt7 = np.arange(0, nb7 * dt7, 0.5)
    pv, cv = 8.0, 2.0                          # tau_agot = 8/2 = 4 s
    acumular_libro(bt7, np.full(bt7.size, 100.0), np.full(bt7.size, pv),
                   np.full(bt7.size, 100.1), np.full(bt7.size, pv), t0=0.0,
                   dt=dt7, nb=nb7, acc=acc7, prev=None)
    tt7 = np.arange(0, nb7 * dt7, 0.25)
    acumular_trades(tt7, np.full(tt7.size, 100.0), np.full(tt7.size, cv * 0.25),
                    rng.random(tt7.size) < 0.5, 0.0, dt7, nb7, acc7)
    g7 = cerrar_rejilla(acc7, 0.0, dt7, nb7)
    chk(abs(np.nanmedian(g7["tau_agot"][g7["valida"]]) - 4.0) < 1e-9
        and "resiliencia" in g7 and g7["resiliencia"].size == nb7,
        "20 la resiliencia se publica y su numerador sigue siendo tau_agot",
        "tau_agot %.6f s" % np.nanmedian(g7["tau_agot"][g7["valida"]]))

    # --- 21. Las rutas se pueden reapuntar sin editar codigo ---------------
    import curvas_estacional as _C
    guarda = (_C.DIR, DIR_V33, DIR_SALIDA)
    try:
        class _A:
            datos = os.path.join("D:", "Maxwell", "telemetria", "estacional")
            v33 = None
            salida = os.path.join("/tmp", "salida_prueba_micelio")
        r21 = configurar_rutas(_A())
        chk(r21["datos"] == _A.datos
            and ruta_cache("estacional").startswith(_A.salida)
            and ruta_salida("x.csv") == os.path.join(_A.salida, "x.csv"),
            "21 --datos y --salida reapuntan lectura y escritura",
            "datos %s   cache %s" % (r21["datos"], ruta_cache("estacional")))
    finally:
        _C.DIR, DIR_V33_, DIR_SALIDA_ = guarda
        globals()["DIR_V33"], globals()["DIR_SALIDA"] = DIR_V33_, DIR_SALIDA_
        import shutil
        shutil.rmtree("/tmp/salida_prueba_micelio", ignore_errors=True)
    chk(_C.DIR == guarda[0] and DIR_SALIDA == guarda[2],
        "21b el control restaura las rutas y no contamina las demas etapas")

    log("")
    log("  %d / %d" % (n_ok, n_tot))
    return 0 if n_ok == n_tot else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", choices=("rutas", "serie", "tau", "dia",
                                        "relacion"))
    _anadir_rutas(ap)
    ap.add_argument("--fuente", default="estacional", choices=("estacional", "v33"))
    ap.add_argument("--tau", default="tau_agot",
                    help="tau_agot | tau_upd | tau_recup_0.50 ...")
    ap.add_argument("--forma", default="P/tau", choices=("P/tau", "tau/P"))
    ap.add_argument("--dt", type=float, default=DT_REJILLA)
    ap.add_argument("--bloque", type=float, default=BLOQUE_S)
    ap.add_argument("--dias", type=float, default=0.0,
                    help="procesar solo los primeros N dias (prueba de humo)")
    ap.add_argument("--desde", help="inicio de ventana, YYYY-MM-DD[THH:MM:SS] UTC")
    ap.add_argument("--hasta", help="fin de ventana, YYYY-MM-DD[THH:MM:SS] UTC")
    ap.add_argument("--zona", default="ny", choices=("ny", "utc"),
                    help="reloj para dia/hora/finde. NY por omision (Sec.4.bis)")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    configurar_rutas(a)
    if a.etapa == "rutas":
        return etapa_rutas(a)
    if a.etapa == "serie":
        return etapa_serie(a)
    if a.etapa == "tau":
        return etapa_tau(a)
    if a.etapa == "dia":
        return etapa_dia(a)
    if a.etapa == "relacion":
        return etapa_relacion(a)
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
