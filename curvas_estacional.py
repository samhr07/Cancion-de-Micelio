# -*- coding: utf-8 -*-
"""
curvas_estacional.py -- v4.1 Sec.2 y Sec.1 sobre `captura_estacional`, POR TRAMO.

    python curvas_estacional.py --muestras [--tramo=K]   construye/cachea muestras
    python curvas_estacional.py --etapa=hp               Sec.2: H_p por tramo
    python curvas_estacional.py --etapa=curvas           Sec.1: las dos curvas

POR QUE UN MODULO NUEVO Y NO UNA BANDERA EN `horizonte.py`:
`horizonte.py` lleva 14 controles verdes y su Sec.1 esta atado a la particion y
al cache de la v3.2. Aqui cambian tres cosas -- la fuente, la particion y el
alcance de la rejilla -- y no conviene tocar un modulo probado para eso. Se
IMPORTAN sus funciones, asi que la aritmetica de las dos curvas es literalmente
la misma que produjo la tabla ya publicada.

QUE CAMBIA RESPECTO A LA CORRIDA DEL 2026-08-10, y todo esta declarado:

1. TRES TRAMOS CONTINUOS por separado, no uno. El proyecto lleva desde la v3.2
   encontrando que nada replica entre capturas con `nu` distinta; con tramos de
   `nu` 18.8 / 67.3 / 32.8 tx/s el regimen pasa a ser una variable MEDIDA en vez
   de una nota al pie. Un resultado que solo aparece en uno de los tres no es un
   resultado, y esa lectura hay que poder hacerla.

2. EMBARGO EN SEGUNDOS, = al mayor `H` probado (4 h), no los 2233 ticks de la
   v3.2. Aquel embargo se derivo para un objetivo a 128 s. Aqui el objetivo se
   mide a 4 h, asi que un embargo de 128 s dejaria la ultima ventana de
   entrenamiento SOLAPANDO con la primera de validacion. Es mas conservador que
   el heredado, no menos.

3. DOS CONTROLES BARAJADOS, no uno. El Sec.1.3 punto 4 pide barajar los signos
   "dentro de cada ventana"; la corrida publicada barajo GLOBALMENTE sobre
   entrenamiento+validacion. Se reportan los dos: el global para poder comparar
   con lo ya publicado, y el de ventana porque es el que el documento pide.

4. NINGUN CONJUNTO DE PRUEBA SE ABRE. Cada tramo se parte 60/20/20 y su ultimo
   20 % queda sin tocar, para que el contraste economico -- si el Sec.1 llega a
   habilitarlo -- tenga holdout propio. El de la v3.2 tampoco se abre.

[!] MEMORIA. La primera version cargaba las 36.5 M de transacciones enteras y
    ADEMAS reescaneaba las 1605 partes de libro una vez por tramo: el proceso
    murio antes de cerrar el primero. Ahora hay un indice de rangos temporales
    cacheado en disco y se procesa UN TRAMO POR INVOCACION, leyendo solo las
    partes que solapan con su ventana.
"""

from __future__ import annotations

import argparse
import gc
import glob
import io
import json
import os
import time

import numpy as np
import pyarrow.parquet as pq

import horizonte as H
import propagador as P

DIR = "telemetria/estacional"
CACHE = "telemetria/muestra_estacional_%02d.npz"
CORTE_TRAMO_S = 300.0        # el mismo criterio de tramo continuo de la v2.2
MIN_HORAS_TRAMO = 20.0       # por debajo no aporta ventanas a la banda que importa
HS = [60, 120, 300, 600, 900, 1800, 3600, 7200, 14400]
EMBARGO_S = float(max(HS))   # ver punto 2 de la cabecera
N_SORTEOS_POTENCIA = 25      # direcciones de senal del control positivo

log, titulo = H.log, H.titulo


# ===========================================================================
# Indice de partes y carga por ventana temporal
# ===========================================================================

def _partes(pref):
    out = []
    for d in sorted(os.listdir(DIR)):
        if d.startswith(pref):
            out += sorted(glob.glob(os.path.join(DIR, d, "*.parquet")))
    return out


def _indice(pref: str) -> list:
    """[(ruta, t_min, t_max, filas)] por parte, cacheado en disco."""
    cache = os.path.join(DIR, "_indice_%s.json" % pref.strip("_"))
    fs = _partes(pref)
    if os.path.exists(cache):
        try:
            d = json.load(io.open(cache, encoding="utf-8"))
            if len(d) == len(fs) and all(a[0] == b for a, b in zip(d, fs)):
                return [tuple(x) for x in d]
        except Exception:
            pass
    out = []
    for f in fs:
        try:
            t = pq.read_table(f, columns=["t"])["t"].to_numpy()
        except Exception:
            continue
        if t.size:
            out.append((f, float(np.min(t)), float(np.max(t)), int(t.size)))
    json.dump(out, io.open(cache, "w", encoding="utf-8"))
    return out


def limites_tramos() -> list:
    """[(t0, t1, n)] de los tramos continuos largos, el mas largo primero.

    [!] SE RECORRE PARTE A PARTE Y SE CACHEA. Esta maquina tiene 7.7 GB de RAM
    con ~2 GB libres: concatenar los 36.5 M de `t` para despues ordenarlos ya
    mataba el proceso. Aqui nunca hay mas de una parte en memoria, y el corte
    entre partes se evalua con la ultima `t` de la anterior.
    """
    cache = os.path.join(DIR, "_tramos.json")
    idx = _indice("trades_")
    if os.path.exists(cache):
        try:
            d = json.load(io.open(cache, encoding="utf-8"))
            if d.get("n_partes") == len(idx) and d.get("ultimo") == idx[-1][2]:
                return [tuple(x) for x in d["tramos"]]
        except Exception:
            pass
    cortes = []          # (t0, t1, n) provisional
    t0 = t_ant = None
    n = 0
    for f, _, _, _ in idx:
        try:
            t = np.sort(pq.read_table(f, columns=["t"])["t"].to_numpy().astype(np.float64))
        except Exception:
            continue
        if t.size == 0:
            continue
        if t_ant is not None and t[0] - t_ant > CORTE_TRAMO_S:
            cortes.append((t0, t_ant, n))
            t0, n = None, 0
        if t0 is None:
            t0 = float(t[0])
        d = np.diff(t)
        for j in np.flatnonzero(d > CORTE_TRAMO_S):
            cortes.append((t0, float(t[j]), n + int(j) + 1))
            t0, n = float(t[j + 1]), -int(j) - 1
        n += t.size
        t_ant = float(t[-1])
    if t0 is not None:
        cortes.append((t0, t_ant, n))
    r = [(a, b, c) for a, b, c in cortes if (b - a) / 3600.0 >= MIN_HORAS_TRAMO]
    r = sorted(r, key=lambda z: -(z[1] - z[0]))
    json.dump({"n_partes": len(idx), "ultimo": idx[-1][2], "tramos": r},
              io.open(cache, "w", encoding="utf-8"))
    return r


def transacciones_en(t0: float, t1: float) -> dict:
    """Transacciones de [t0, t1] en arrays PREASIGNADOS.

    `np.concatenate` sobre una lista de partes mantiene la lista y el resultado
    a la vez: dos copias del tramo. Con 20 M de ticks eso son 1.5 GB y no caben.
    Se cuenta primero y se rellena despues.
    """
    sel = [(f, a, b) for f, a, b, _ in _indice("trades_") if not (b < t0 or a > t1)]
    tot = 0
    trozos = []
    for f, _, _ in sel:
        try:
            tb = pq.read_table(f, columns=["t", "precio", "cant", "maker"])
        except Exception:
            continue
        tt = tb["t"].to_numpy().astype(np.float64)
        s = (tt >= t0) & (tt <= t1)
        k = int(s.sum())
        if k:
            trozos.append((f, s, k))
            tot += k
        del tb, tt
    d = {"t": np.empty(tot, np.float64), "precio": np.empty(tot, np.float64),
         "cant": np.empty(tot, np.float32), "maker": np.empty(tot, bool)}
    o = 0
    for f, s, k in trozos:
        tb = pq.read_table(f, columns=["t", "precio", "cant", "maker"])
        d["t"][o:o + k] = tb["t"].to_numpy().astype(np.float64)[s]
        d["precio"][o:o + k] = tb["precio"].to_numpy().astype(np.float64)[s]
        d["cant"][o:o + k] = tb["cant"].to_numpy().astype(np.float32)[s]
        d["maker"][o:o + k] = tb["maker"].to_numpy().astype(bool)[s]
        o += k
        del tb
    del trozos
    gc.collect()
    if not np.all(np.diff(d["t"]) >= 0):
        # No deberia pasar: la captura escribe en orden. Si pasa, se ordena.
        ordn = np.argsort(d["t"], kind="stable")
        for c in list(d):
            d[c] = d[c][ordn]
        del ordn
        gc.collect()
    # Ceros del feed (limitacion 3 del Sec.12 de la v3.2): fuera en origen.
    ok = np.isfinite(d["precio"]) & (d["precio"] > 0) & (d["cant"] > 0)
    if not ok.all():
        for c in list(d):
            d[c] = d[c][ok]
    del ok
    gc.collect()
    return d


def alinear_libro(tt: np.ndarray) -> dict:
    """Ultimo bookTicker ESTRICTAMENTE anterior a cada transaccion (Sec.4.1).

    Se recorre el libro PARTE A PARTE en orden temporal arrastrando la ultima
    fila de la anterior. De todo el libro solo se necesita su valor en los
    instantes de transaccion, asi que nunca hay mas de una parte en memoria --
    que con 223 M de filas y 2 GB libres es la unica forma de que quepa.
    """
    n = tt.size
    bid = np.full(n, np.nan)
    ask = np.full(n, np.nan)
    t0, t1 = float(tt[0]), float(tt[-1])
    ult = None                       # (t, b, a) arrastrada de la parte anterior
    for f, ra, rb, _ in _indice("libro_"):
        if rb < t0 - 3600.0 or ra > t1 + 60.0:
            continue
        try:
            tb = pq.read_table(f, columns=["t", "b", "a"])
        except Exception:
            continue
        bt = tb["t"].to_numpy().astype(np.float64)
        bb = tb["b"].to_numpy().astype(np.float64)
        ba = tb["a"].to_numpy().astype(np.float64)
        del tb
        ok = np.isfinite(bb) & (bb > 0) & np.isfinite(ba) & (ba > 0)
        bt, bb, ba = bt[ok], bb[ok], ba[ok]
        if bt.size == 0:
            continue
        o = np.argsort(bt, kind="stable")
        bt, bb, ba = bt[o], bb[o], ba[o]
        if ult is not None:
            bt = np.r_[ult[0], bt]
            bb = np.r_[ult[1], bb]
            ba = np.r_[ult[2], ba]
        lo = int(np.searchsorted(tt, bt[0], side="left"))
        hi = int(np.searchsorted(tt, bt[-1], side="right"))
        if hi > lo:
            j = np.searchsorted(bt, tt[lo:hi], side="left") - 1
            v = j >= 0
            idx = np.arange(lo, hi)[v]
            bid[idx] = bb[j[v]]
            ask[idx] = ba[j[v]]
            del j, v, idx
        ult = (float(bt[-1]), float(bb[-1]), float(ba[-1]))
        del bt, bb, ba
    gc.collect()
    return {"bid": bid, "ask": ask}


def construir(tr: dict) -> dict:
    lb = alinear_libro(tr["t"])
    bid, ask = lb["bid"], lb["ask"]
    del lb
    ok = np.isfinite(bid) & np.isfinite(ask) & (ask >= bid)
    cob = float(ok.mean())
    m = {}
    # Se vacia `tr` a medida que se copia: dos copias del tramo no caben.
    for orig, dest in (("t", "t"), ("precio", "precio"), ("cant", "cant"),
                       ("maker", "maker")):
        m[dest] = tr[orig][ok] if not ok.all() else tr[orig]
        tr[orig] = None
        gc.collect()
    m["bid"] = bid[ok] if not ok.all() else bid
    del bid
    m["ask"] = ask[ok] if not ok.all() else ask
    del ask
    gc.collect()
    m["eps"] = P.signo_transaccion(m.pop("maker")).astype(np.int8)
    m["cobertura_libro"] = np.array([cob])
    return m


def cachear_muestras(args) -> int:
    titulo("MUESTRAS POR TRAMO -- alineado libro-transaccion del Sec.4.1")
    lim = limites_tramos()
    log("tramos continuos de >= %.0f h: %d" % (MIN_HORAS_TRAMO, len(lim)))
    for k, (t0, t1, n) in enumerate(lim):
        log("   tramo %d: %8.2f h  %10d ticks  nu = %6.2f tx/s"
            % (k, (t1 - t0) / 3600.0, n, n / (t1 - t0)))
    objetivo = range(len(lim)) if args.tramo < 0 else [args.tramo]
    for k in objetivo:
        if k >= len(lim):
            continue
        if os.path.exists(CACHE % k) and not args.rehacer:
            log("")
            log("--- tramo %d: ya cacheado, se omite (--rehacer para forzar)" % k)
            continue
        t0, t1, n = lim[k]
        t_ini = time.time()
        log("")
        log("--- tramo %d: %.2f h, %d ticks" % (k, (t1 - t0) / 3600.0, n))
        tr = transacciones_en(t0, t1)
        log("    transacciones validas: %d" % tr["t"].size)
        m = construir(tr)
        del tr
        np.savez_compressed(CACHE % k, **m)
        sp = m["ask"] - m["bid"]
        mid = 0.5 * (m["bid"] + m["ask"])
        ym = np.diff(mid)
        g0 = float(np.mean(ym * m["eps"][1:]))
        log("    alineadas %d (%.2f %% de cobertura de libro)"
            % (m["t"].size, 100 * float(m["cobertura_libro"][0])))
        log("    precio %.1f - %.1f   spread mediano %.4f USD (p90 %.4f)"
            % (m["precio"].min(), m["precio"].max(),
               float(np.median(sp)), float(np.percentile(sp, 90))))
        log("    y_mid nulos %.2f %%   E[y_mid*eps] = %+.6f   %s"
            % (100 * float(np.mean(ym == 0)), g0,
               "[OK] G(0) > 0" if g0 > 0 else "*** G(0) <= 0: PARAR ***"))
        log("    cache: %s   (%.1f s)" % (CACHE % k, time.time() - t_ini))
        del m
    return 0


def _tramos_cache():
    return [k for k in range(12) if os.path.exists(CACHE % k)]


# ===========================================================================
# Sec.2 -- H_p por tramo
# ===========================================================================

def _rejilla_H(dur_s: float):
    """Rejilla densa; el tope es el mayor `H` con >= 30 ventanas del tramo."""
    tope = min(14400.0, dur_s / 30.0)
    return sorted(set([0.5, 1, 2, 5] +
                      [round(x, 1) for x in np.logspace(
                          np.log10(8), np.log10(max(tope, 60.0)), 34)]))


def hp_tramo(m: dict, semilla: int = 5) -> dict:
    p, t = np.asarray(m["precio"], float), np.asarray(m["t"], float)
    dur = float(t[-1] - t[0])
    Hs = _rejilla_H(dur)
    nulo = H.barajar_incrementos(p, semilla)
    fp = H.firma_pared(p, t, Hs)
    fn = H.firma_pared(nulo, t, Hs)
    out = {"filas": fp, "filas_nulo": fn,
           "alcance": max(f["H"] for f in fp), "ventanas": []}
    for lo, hi in ((10, 100), (30, 300), (60, 600), (10, 600),
                   (60, 3600), (300, 14400)):
        a = H.pendiente(fp, "H", "sigma_por_raiz_H", lo, hi)
        b = H.pendiente(fn, "H", "sigma_por_raiz_H", lo, hi)
        if a["n"] < 3:
            continue
        out["ventanas"].append({
            "rango": (lo, hi), "n": a["n"],
            "H_p": a["pendiente"] + 0.5, "H_p_nulo": b["pendiente"] + 0.5,
            "H_p_corr": a["pendiente"] - b["pendiente"] + 0.5,
            "sesgo": abs(b["pendiente"])})
    return out


def etapa_hp(args) -> int:
    titulo("Sec.2 -- H_p POR TRAMO (reloj de pared, control barajado en todos)")
    log("")
    log("  El Sec.2.2 resolvio por BANDA porque con n <= 4 capturas una correlacion")
    log("  no es evidencia. Estos tramos anaden replicas, no sustituyen la banda.")
    res = {}
    for k in _tramos_cache():
        m = dict(np.load(CACHE % k))
        dur = float(m["t"][-1] - m["t"][0])
        r = hp_tramo(m)
        res[k] = r
        log("")
        log("--- tramo %d:  %.2f h   nu = %.2f tx/s   precio %.0f - %.0f"
            % (k, dur / 3600.0, m["t"].size / dur,
               m["precio"].min(), m["precio"].max()))
        log("    alcance verificado de la firma: %.0f s" % r["alcance"])
        log("    %-14s %4s   %9s %9s %10s   %s"
            % ("ventana [s]", "n", "H_p", "barajado", "corregido", "|sesgo|"))
        for v in r["ventanas"]:
            log("    %-14s %4d   %+9.4f %+9.4f %+10.4f   %.4f"
                % ("[%d, %d]" % v["rango"], v["n"], v["H_p"], v["H_p_nulo"],
                   v["H_p_corr"], v["sesgo"]))
        del m
    log("")
    todos = [v["H_p_corr"] for k in res for v in res[k]["ventanas"]]
    if todos:
        log("  BANDA de H_p sobre los tramos y todas las ventanas de ajuste:")
        log("    min %.4f   mediana %.4f   max %.4f   (n = %d ajustes)"
            % (min(todos), float(np.median(todos)), max(todos), len(todos)))
    log("")
    log("  [!] NUNCA se elige el H_p que hace cruzar las curvas (Sec.2.2).")
    return 0


# ===========================================================================
# Sec.1 -- las dos curvas, por tramo
# ===========================================================================

def curvas_tramo(m: dict, hp_banda, semilla=11):
    t = np.asarray(m["t"], float)
    eps = np.asarray(m["eps"], float)
    cant = np.asarray(m["cant"], float)
    mid = 0.5 * (m["bid"] + m["ask"])
    n = t.size
    dur = float(t[-1] - t[0])
    nu = n / dur
    pmed = float(np.median(m["precio"]))

    fp = H.firma_pared(m["precio"], t, _rejilla_H(dur))
    s1 = {hp: H.sigma1_pb(fp, pmed, hp) for _, hp in hp_banda}

    # Particion 60/20/20 EN TIEMPO, con embargo de EMBARGO_S entre bloques.
    t_e = t[0] + 0.60 * dur
    t_v = t[0] + 0.80 * dur
    flujo = np.concatenate(([0.0], np.cumsum(eps * cant)))
    rng = np.random.default_rng(semilla)

    # control global: signos barajados sobre entrenamiento+validacion
    lim = int(np.searchsorted(t, t_v, side="left"))
    eps_g = eps.copy()
    eps_g[:lim] = rng.permutation(eps_g[:lim])
    flujo_g = np.concatenate(([0.0], np.cumsum(eps_g * cant)))

    filas = []
    for Hh in HS:
        idx = H.ventanas_no_solapadas(t, float(Hh))
        idx = idx[idx < n - 1]
        fin = np.clip(np.searchsorted(t, t[idx] + Hh, side="left"), 0, n - 1)
        val = (fin > idx) & (mid[idx] > 0) & (mid[fin] > 0)
        idx, fin = idx[val], fin[val]
        if idx.size < 12:
            continue
        y = np.log(mid[fin] / mid[idx])
        ret = [max(1, int(round(f * Hh * nu))) for f in (0.25, 0.5, 1.0, 2.0)]
        X = H.rasgos_flujo(flujo, idx, ret)
        # la ventana entera tiene que caber en su bloque: se juzga por su FIN
        ent = t[fin] <= t_e
        vld = (t[idx] >= t_e + EMBARGO_S) & (t[fin] <= t_v)
        if ent.sum() < 8 or vld.sum() < 4:
            filas.append({"H": Hh, "n": int(vld.sum()), "insuf": True})
            continue
        r2 = H.r2_fuera_de_muestra(X[ent], y[ent], X[vld], y[vld])
        Xg = H.rasgos_flujo(flujo_g, idx, ret)
        r2g = H.r2_fuera_de_muestra(Xg[ent], y[ent], Xg[vld], y[vld])
        # control por VENTANA: se baraja el signo dentro de cada ventana de H s
        eps_w = eps.copy()
        bordes = np.searchsorted(t, np.arange(t[0], t[-1], Hh), side="left")
        for u, v in zip(bordes[:-1], bordes[1:]):
            if v > u + 1:
                eps_w[u:v] = rng.permutation(eps_w[u:v])
        flujo_w = np.concatenate(([0.0], np.cumsum(eps_w * cant)))
        Xw = H.rasgos_flujo(flujo_w, idx, ret)
        r2w = H.r2_fuera_de_muestra(Xw[ent], y[ent], Xw[vld], y[vld])
        filas.append({"H": Hh, "n": int(vld.sum()), "n_ent": int(ent.sum()),
                      "r2": r2, "r2_glob": r2g, "r2_vent": r2w,
                      "req": [H.r2_requerido(float(Hh), s1[hp], hp)
                              for _, hp in hp_banda], "insuf": False})
    return filas, s1


def etapa_curvas(args) -> int:
    titulo("Sec.1 -- LAS DOS CURVAS, POR TRAMO")
    hp_banda = [("banda baja", 0.371), ("difusivo", 0.500), ("banda alta", 0.552)]
    log("")
    log("  particion 60/20/20 EN TIEMPO por tramo, embargo de %.0f s (= el mayor H"
        % EMBARGO_S)
    log("  probado). El ultimo 20 %% de cada tramo NO SE ABRE, ni el de la v3.2.")
    log("  Ventanas NO solapadas; una ventana pertenece a un bloque solo si CABE")
    log("  entera en el.")

    resumen = []
    for k in _tramos_cache():
        m = dict(np.load(CACHE % k))
        dur = float(m["t"][-1] - m["t"][0])
        filas, s1 = curvas_tramo(m, hp_banda)
        log("")
        log("=" * 78)
        log("TRAMO %d   %.2f h   nu = %.2f tx/s   precio %.0f - %.0f  (recorrido %.1f %%)"
            % (k, dur / 3600.0, m["t"].size / dur, m["precio"].min(),
               m["precio"].max(),
               100 * (m["precio"].max() - m["precio"].min()) / m["precio"].min()))
        log("=" * 78)
        for nom, hp in hp_banda:
            log("    sigma_1 (%-11s H_p = %.3f) = %.6f pb*s^(-H_p)"
                % (nom, hp, s1[hp]))
        log("")
        log("     H     n_ent  n_val | R2_medido  bar_global  bar_ventana |"
            "   R2_req  0.371    0.500    0.552")
        log("   " + "-" * 96)
        for f in filas:
            if f["insuf"]:
                log("   %5d s          %5d | (insuficiente)" % (f["H"], f["n"]))
                continue
            marca = "" if f["n"] >= 30 else "  <- PILOTO"
            log("   %5d s  %7d %6d | %+9.4f  %+9.4f   %+9.4f | %8.2f %% %7.2f %% %7.2f %%%s"
                % (f["H"], f["n_ent"], f["n"], f["r2"], f["r2_glob"],
                   f["r2_vent"], 100 * f["req"][0], 100 * f["req"][1],
                   100 * f["req"][2], marca))
            ctrl = max(abs(f["r2_glob"]), abs(f["r2_vent"]))
            cruza = [nom for (nom, _), q in zip(hp_banda, f["req"]) if f["r2"] > q]
            if cruza and f["n"] >= 30 and f["r2"] > 3 * ctrl:
                log("            -> CRUZA con: %s   (control limpio)"
                    % ", ".join(cruza))
            resumen.append((k, f))
        log("")
        del m

    # --- lectura del Sec.1.4, aplicada literalmente ---
    titulo("REGLA DE LECTURA DEL Sec.1.4, APLICADA")
    dec = [(k, f) for k, f in resumen if f["n"] >= 30]
    log("")
    log("  filas con n_val >= 30 (decidibles): %d de %d" % (len(dec), len(resumen)))
    if not dec:
        log("  -> NO DECIDIBLE en todo el rango: ninguna H alcanza 30 ventanas.")
        return 0
    cruces = [(k, f) for k, f in dec
              if any(f["r2"] > q for q in f["req"])
              and f["r2"] > 3 * max(abs(f["r2_glob"]), abs(f["r2_vent"]))]
    log("  filas decidibles que CRUZAN con control limpio  : %d" % len(cruces))
    for k, f in cruces:
        log("      tramo %d  H = %d s  R2 = %+.4f" % (k, f["H"], f["r2"]))
    log("")
    if not cruces:
        log("  *** DESENLACE: no hay banda viable con este predictor. En ninguna H")
        log("      con potencia el R2 medido alcanza al requerido, en ninguno de los")
        log("      regimenes de nu medidos. Es el segundo desenlace de la tabla del")
        log("      Sec.1.4, y el Sec.1.1 lo declaro admisible antes de medir. ***")
    else:
        log("  *** HAY BANDA VIABLE. Se preregistra el contraste economico. ***")
    log("")
    log("  margen mas favorable (R2_medido contra el R2_req mas benevolo):")
    mej = sorted(dec, key=lambda z: -(z[1]["r2"] / max(min(z[1]["req"]), 1e-12)))[:8]
    for k, f in mej:
        log("      tramo %d  H = %5d s  n = %4d  R2 = %+.4f  req_min = %.4f  razon = %.4f"
            % (k, f["H"], f["n"], f["r2"], min(f["req"]),
               f["r2"] / max(min(f["req"]), 1e-12)))
    return 0



# ===========================================================================
# Control POSITIVO de potencia -- sin esto, "R2 = 0" no es un resultado
# ===========================================================================

def potencia_tramo(m: dict, hp_banda, semilla=23):
    """Se INYECTA una senal de tamano exactamente `R2_req` y se comprueba que
    el mismo procedimiento la recupera.

    [!] ES OBLIGATORIO Y ES LA CUARTA REGLA DEL Sec.1.4. Un `R2` medido de cero
    admite dos lecturas -- "no hay senal" y "no hay potencia" -- y solo una de
    ellas decide el proyecto. La cuarta fila de la tabla del Sec.1.4 manda no
    leer una `H` donde el barajado es comparable al real; aqui los dos salen
    ~0, asi que hay que demostrar por separado que una senal del tamano
    REQUERIDO habria sido visible. Si el control recupera `R2_req`, "cero"
    significa cero.

    La senal inyectada es una combinacion lineal FIJA de los mismos rasgos, asi
    que es recuperable por construccion: lo que se pone a prueba no es el
    predictor sino el tamano muestral.
    """
    t = np.asarray(m["t"], float)
    eps = np.asarray(m["eps"], float)
    cant = np.asarray(m["cant"], float)
    mid = 0.5 * (m["bid"] + m["ask"])
    n = t.size
    dur = float(t[-1] - t[0])
    nu = n / dur
    pmed = float(np.median(m["precio"]))
    fp = H.firma_pared(m["precio"], t, _rejilla_H(dur))
    s1 = {hp: H.sigma1_pb(fp, pmed, hp) for _, hp in hp_banda}
    t_e = t[0] + 0.60 * dur
    t_v = t[0] + 0.80 * dur
    flujo = np.concatenate(([0.0], np.cumsum(eps * cant)))
    rng = np.random.default_rng(semilla)

    filas = []
    for Hh in HS:
        idx = H.ventanas_no_solapadas(t, float(Hh))
        idx = idx[idx < n - 1]
        fin = np.clip(np.searchsorted(t, t[idx] + Hh, side="left"), 0, n - 1)
        val = (fin > idx) & (mid[idx] > 0) & (mid[fin] > 0)
        idx, fin = idx[val], fin[val]
        if idx.size < 12:
            continue
        y = np.log(mid[fin] / mid[idx])
        ret = [max(1, int(round(f * Hh * nu))) for f in (0.25, 0.5, 1.0, 2.0)]
        X = H.rasgos_flujo(flujo, idx, ret)
        ent = t[fin] <= t_e
        vld = (t[idx] >= t_e + EMBARGO_S) & (t[fin] <= t_v)
        if ent.sum() < 8 or vld.sum() < 30:
            continue
        req = min(H.r2_requerido(float(Hh), s1[hp], hp) for _, hp in hp_banda)
        req = float(min(max(req, 1e-4), 0.95))
        # [!] NORMALIZACION POR BLOQUE, y las dos versiones anteriores estaban
        # mal por la misma razon con distinto disfraz. Con estadisticos de
        # ENTRENAMIENTO la senal llegaba encogida a validacion; con
        # estadisticos GLOBALES seguia contaminada porque la varianza del
        # retorno real cambia mucho de un bloque a otro (agrupamiento de
        # volatilidad), y entonces el control media heterocedasticidad y no
        # potencia. La pregunta correcta es: "si en el bloque de VALIDACION
        # una fraccion R2_req de la varianza fuera predecible, la veriamos?".
        # Eso exige que la fraccion inyectada sea req DENTRO de cada bloque.
        # Normalizar asi no es mirar el dato: el objetivo sintetico se
        # construye entero y la medicion real del Sec.1 no lo usa.
        # [!] SE PROMEDIA SOBRE VARIAS DIRECCIONES. Una sola `w` es un
        # estimador ruidoso de la potencia: los cuatro rasgos son sumas
        # acumuladas anidadas y estan muy correlacionados, asi que hay
        # direcciones casi degeneradas en las que cualquier metodo falla. La
        # potencia de la medicion es la TIPICA, no la de un sorteo.
        rec = []
        esc = []
        for _ in range(N_SORTEOS_POTENCIA):
            w = rng.normal(size=X.shape[1])
            y_syn = np.zeros_like(y)
            esc = []
            for blq in (ent, vld):
                Z = (X[blq] - X[blq].mean(0)) / (X[blq].std(0) + 1e-30)
                g = Z @ w
                g = (g - g.mean()) / (g.std() + 1e-30)
                yb = (y[blq] - y[blq].mean()) / (y[blq].std() + 1e-30)
                y_syn[blq] = np.sqrt(1.0 - req) * yb + np.sqrt(req) * g
                esc.append(float(np.std(y[blq])))
            rec.append(H.r2_fuera_de_muestra(X[ent], y_syn[ent], X[vld], y_syn[vld]))
        r2 = float(np.median(rec))
        r2_p10 = float(np.percentile(rec, 10))
        esc = esc[1] / max(esc[0], 1e-30)
        filas.append({"H": Hh, "n": int(vld.sum()), "req": req,
                      "r2_rec": r2, "r2_p10": r2_p10, "esc": esc})
    return filas


def etapa_potencia(args) -> int:
    titulo("CONTROL POSITIVO DE POTENCIA (cuarta regla del Sec.1.4)")
    hp_banda = [("banda baja", 0.371), ("difusivo", 0.500), ("banda alta", 0.552)]
    log("")
    log("  Se inyecta en el objetivo una senal de tamano igual al R2 REQUERIDO mas")
    log("  benevolo de la banda y se mide con EL MISMO procedimiento. Si vuelve,")
    log("  el 'R2 = 0' del Sec.1 significa ausencia de senal y no falta de potencia.")
    fallos = 0
    for k in _tramos_cache():
        m = dict(np.load(CACHE % k))
        filas = potencia_tramo(m, hp_banda)
        del m
        log("")
        log("--- tramo %d" % k)
        log("     H      n_val |  R2 inyectado   R2 recup MED    p10     razon   sd(y)val/ent")
        for f in filas:
            raz = f["r2_rec"] / max(f["req"], 1e-12)
            ok = "OK" if raz > 0.5 else "*** SIN POTENCIA ***"
            if raz <= 0.5:
                fallos += 1
            log("   %5d s %7d |    %8.4f      %+9.4f %+8.4f %6.3f  %8.3f   %s"
                % (f["H"], f["n"], f["req"], f["r2_rec"], f["r2_p10"],
                   raz, f["esc"], ok))
    log("")
    if fallos == 0:
        log("  *** LA MEDICION TIENE POTENCIA en todas las filas decidibles: una senal")
        log("      del tamano requerido se recupera. El cero del Sec.1 es un cero. ***")
    else:
        log("  *** %d fila(s) SIN POTENCIA: en esas H el Sec.1 no se lee. ***" % fallos)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--muestras", action="store_true")
    ap.add_argument("--tramo", type=int, default=-1)
    ap.add_argument("--rehacer", action="store_true")
    ap.add_argument("--etapa", default="")
    a = ap.parse_args(argv)
    if a.muestras:
        return cachear_muestras(a)
    if a.etapa == "hp":
        return etapa_hp(a)
    if a.etapa == "curvas":
        return etapa_curvas(a)
    if a.etapa == "potencia":
        return etapa_potencia(a)
    log("nada que hacer: usa --muestras, --etapa=hp, =curvas o =potencia")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
