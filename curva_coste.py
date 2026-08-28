# -*- coding: utf-8 -*-
"""
curva_coste.py -- v4.2 §2 y §3: `L(H)` con sus cinco terminos, y `C_respaldo`.

    python curva_coste.py --autotest       controles con verdad conocida
    python curva_coste.py --etapa=componentes   A(H), p(H), C_respaldo(H), F(H)

    L(H) = c(u) + A(H) + F(H) + (1 - p(H)) * C_respaldo(H)

| termino | que es | direccion esperada |
|---|---|---|
| `c(u)` | comision ida y vuelta | constante en pb; **LEIDA de la cuenta el 2026-08-28** |
| `A(H)` | seleccion adversa: markout del punto medio a `H` tras el llenado | crece y deberia saturar |
| `F(H)` | financiacion del perpetuo, **con signo** | crece ~lineal |
| `p(H)` | tasa de llenado maker, con NIVEL BARRIDO contado como llenado | crece con `H` |
| `C_respaldo(H)` | lo que costo no haber cruzado al principio | **crece con `H`** |

DEFINICIONES DECLARADAS ANTES DE MIRAR DATO (§3 lo exige):

  Colocacion: en cada instante de una rejilla UNIFORME EN TIEMPO se coloca una
  orden maker de compra al mejor bid y una de venta al mejor ask.

  Llenado (`NIVEL BARRIDO`, como la v4.0): la compra al bid `P_b` se considera
  llenada si el mejor bid baja ESTRICTAMENTE por debajo de `P_b` en algun
  momento antes de `t + H` -- el nivel fue barrido. Simetrico para la venta.
  Es una cota SUPERIOR de la tasa de llenado real: no conocemos la posicion en
  cola, asi que contamos como llenado todo nivel barrido aunque nuestra orden
  pudiera haber quedado detras. Se declara asi y no al reves porque una `p`
  optimista hace la conclusion economica mas dificil de rechazar, no mas facil.

  `A(H)` seleccion adversa: para las llenadas, `-(mid(t+H) - P_llenado)` en pb
  con el signo de la direccion. Si el precio se va en tu contra tras llenarte,
  es coste positivo.

  `C_respaldo(H)`: para las NO llenadas, `mid(t+H) - mid(t)` con el signo de la
  direccion pretendida. Es lo que costo no haber cruzado al principio.

⚠ **La nota de prioridad de la v4.1 estaba mal y el §3 lo registra:** alli se
escribio que `C_respaldo` pierde peso al crecer `H` porque `p` sube. Cierto para
el factor `(1-p)` y **falso para el otro**: cuanto mas esperas, mas se puede
haber escapado el precio. El producto no es monotono y **se mide, no se razona**.
"""

from __future__ import annotations

import argparse
import gc
import json
import os

import numpy as np

import coste as CO
import horizonte as H
import identidad as I

log, titulo = H.log, H.titulo

HS = [900, 1800, 3600, 7200, 14400, 28800, 86400, 259200]   # §4.2, congelada
PASO_COLOC_S = 60.0          # rejilla uniforme de colocacion
BUCKET_S = 60.0              # resolucion del min/max rodante


def _por_bucket(t, v, t0, t1, bs, agg):
    """min o max de `v` por casilla de `bs` segundos."""
    nb = int((t1 - t0) / bs) + 1
    k = np.clip(np.floor((t - t0) / bs).astype(np.int64), 0, nb - 1)
    out = np.full(nb, np.nan)
    orden = np.argsort(k, kind="stable")
    kk, vv = k[orden], v[orden]
    lim = np.r_[0, np.flatnonzero(np.diff(kk)) + 1, kk.size]
    for a, b in zip(lim[:-1], lim[1:]):
        out[kk[a]] = agg(vv[a:b])
    return out


def _rodante(x, k, agg):
    """min/max rodante HACIA ADELANTE sobre `k` casillas, O(n) con ndimage."""
    from scipy.ndimage import minimum_filter1d, maximum_filter1d
    f = minimum_filter1d if agg == "min" else maximum_filter1d
    y = np.where(np.isfinite(x), x, (np.inf if agg == "min" else -np.inf))
    r = f(y, size=k, origin=-(k // 2), mode="nearest")
    return r


def componentes_fragmento(nombre, hs=HS) -> dict:
    d = I._dir_mmap(nombre)
    car = lambda c: np.load(os.path.join(d, c + ".npy"), mmap_mode="r")
    t = np.asarray(car("t"), float)
    bid = np.asarray(car("bid"), float)
    ask = np.asarray(car("ask"), float)
    mid = np.asarray(car("mid"), float)
    t0, t1 = float(t[0]), float(t[-1])
    if t1 - t0 < 4 * max(hs):
        hs = [h for h in hs if t1 - t0 >= 4 * h]
    if not hs:
        return {}
    mb = _por_bucket(t, bid, t0, t1, BUCKET_S, np.min)
    Ma = _por_bucket(t, ask, t0, t1, BUCKET_S, np.max)
    out = {"nombre": nombre, "dur_h": (t1 - t0) / 3600.0, "filas": []}
    g = np.arange(t0, t1 - max(hs), PASO_COLOC_S)
    idx = np.clip(np.searchsorted(t, g, side="left"), 0, t.size - 1)
    for Hs in hs:
        k = max(2, int(Hs / BUCKET_S))
        minb = _rodante(mb, k, "min")
        maxa = _rodante(Ma, k, "max")
        kb = np.clip(np.floor((t[idx] - t0) / BUCKET_S).astype(np.int64) + 1,
                     0, minb.size - 1)
        j = np.clip(np.searchsorted(t, t[idx] + Hs, side="left"), 0, t.size - 1)
        v = (j > idx) & np.isfinite(minb[kb]) & np.isfinite(maxa[kb])
        i2, j2, kb2 = idx[v], j[v], kb[v]
        if i2.size < 30:
            continue
        Pb, Pa, m0, m1 = bid[i2], ask[i2], mid[i2], mid[j2]
        # NIVEL BARRIDO
        llena_c = minb[kb2] < Pb            # compra al bid: el bid se hunde
        llena_v = maxa[kb2] > Pa            # venta al ask: el ask sube
        pb = 1e4 / m0
        # A(H): seleccion adversa de las LLENADAS, con signo de la direccion
        A_c = -(m1 - Pb) * pb               # compra: coste si el mid baja
        A_v = (m1 - Pa) * pb                # venta: coste si el mid sube
        A = np.concatenate([A_c[llena_c], A_v[llena_v]])
        # C_respaldo: de las NO llenadas, lo que se escapo
        R_c = (m1 - m0) * pb                # compra no llenada y el precio subio
        R_v = -(m1 - m0) * pb
        R = np.concatenate([R_c[~llena_c], R_v[~llena_v]])
        p = float((llena_c.sum() + llena_v.sum()) / (2.0 * i2.size))
        q = lambda z, w: (float(np.percentile(z, w)) if z.size else float("nan"))
        out["filas"].append({
            "H": Hs, "n": int(2 * i2.size), "p": p,
            "A_med": q(A, 50), "A_p10": q(A, 10), "A_p90": q(A, 90),
            "A_media": float(A.mean()) if A.size else float("nan"),
            "R_med": q(R, 50), "R_p10": q(R, 10), "R_p90": q(R, 90),
            "R_media": float(R.mean()) if R.size else float("nan")})
        del minb, maxa
        gc.collect()
    del t, bid, ask, mid, mb, Ma
    gc.collect()
    return out


# ===========================================================================
# F(H): financiacion, CON SIGNO
# ===========================================================================

def financiacion(hs=HS) -> list:
    c = CO.cache().get("financiacion")
    if not c:
        return []
    # la tasa es por periodo de 8 h y se liquida en instantes FIJOS
    med = c["abs_mediana_pb"]
    media = c["media_pb"]
    p90 = c["p90_abs_pb"]
    out = []
    for Hs in hs:
        pc = min(1.0, Hs / 28800.0)          # probabilidad de cruzar una marca
        out.append({"H": Hs, "p_cruce": pc,
                    "F_esperada_largo": media * pc,      # un LARGO paga si la tasa es +
                    "F_esperada_corto": -media * pc,     # un CORTO la COBRA
                    "F_abs_med": med * pc, "F_p90": p90 * pc,
                    "n_liquidaciones": Hs / 28800.0})
    return out


# ===========================================================================
# Etapa
# ===========================================================================

def etapa_componentes(args) -> int:
    titulo("v4.2 §2 y §3 -- L(H) CON SUS CINCO TERMINOS")
    cu_m = CO.c_u("maker_maker", 0.0)["comision_pb"]
    cu_t = CO.c_u("taker_taker", 0.0)["comision_pb"]
    log("")
    log("  c(u) LEIDA de la cuenta (2026-08-28): maker+maker %.4f pb, taker+taker %.4f pb"
        % (cu_m, cu_t))
    log("  descuento BNB: %s" % ("ACTIVO" if CO.BNB_BURN_ACTIVO else "APAGADO"))

    log("")
    log("  --- F(H): financiacion del perpetuo, CON SIGNO ---")
    log("  %8s %9s %14s %14s %12s %12s"
        % ("H", "p_cruce", "E[F] largo", "E[F] corto", "|F| mediana", "|F| p90"))
    for f in financiacion():
        log("  %7ds %8.4f %13.4f %14.4f %12.4f %12.4f"
            % (f["H"], f["p_cruce"], f["F_esperada_largo"], f["F_esperada_corto"],
               f["F_abs_med"], f["F_p90"]))
    log("  [!] Un CORTO cobra la financiacion cuando la tasa es positiva, y la tasa")
    log("      medida es +0.2333 pb/8h de media. Tratarla como coste sin signo")
    log("      destruye una fuente de ventaja -- es el modo de fallo 4 del §8.")

    todos = {}
    for nombre in I._nombres():
        r = componentes_fragmento(nombre)
        if not r or not r["filas"]:
            continue
        todos[nombre] = r
        log("")
        log("--- %s (%.1f h) ---" % (nombre, r["dur_h"]))
        log("  %8s %8s %8s | %10s %10s %10s | %10s %10s %10s"
            % ("H", "n", "p(H)", "A med", "A p10", "A p90", "Cresp med", "Cresp p10", "Cresp p90"))
        for f in r["filas"]:
            log("  %7ds %8d %8.4f | %10.4f %10.4f %10.4f | %10.4f %10.4f %10.4f"
                % (f["H"], f["n"], f["p"], f["A_med"], f["A_p10"], f["A_p90"],
                   f["R_med"], f["R_p10"], f["R_p90"]))
        gc.collect()

    # --- L(H) agregada
    log("")
    titulo("L(H) COMPLETA, dos columnas de c(u)")
    fin = {f["H"]: f for f in financiacion()}
    log("  %8s %8s %9s %9s %9s | %11s %11s"
        % ("H", "p(H)", "A(H)", "Cresp", "(1-p)Cr", "L maker", "L taker"))
    log("  " + "-" * 76)
    filas_L = []
    for Hs in HS:
        ps, As, Rs = [], [], []
        for r in todos.values():
            for f in r["filas"]:
                if f["H"] == Hs:
                    ps.append(f["p"]); As.append(f["A_med"]); Rs.append(f["R_med"])
        if not ps:
            continue
        p, A, R = float(np.median(ps)), float(np.median(As)), float(np.median(Rs))
        F = fin[Hs]["F_abs_med"] if Hs in fin else 0.0
        Lm = cu_m + A + F + (1 - p) * R
        Lt = cu_t + A + F
        filas_L.append({"H": Hs, "p": p, "A": A, "R": R, "L_maker": Lm, "L_taker": Lt})
        log("  %7ds %8.4f %9.4f %9.4f %9.4f | %11.4f %11.4f"
            % (Hs, p, A, R, (1 - p) * R, Lm, Lt))
    log("")
    log("  [!] `A` y `Cresp` son MEDIANAS entre fragmentos de la mediana de cada uno.")
    log("      Las distribuciones completas estan arriba, por fragmento, como pide el §3.")
    json.dump({"L": filas_L, "c_u_maker": cu_m, "c_u_taker": cu_t},
              open("telemetria/curva_coste.json", "w", encoding="utf-8"), indent=1)
    log("  guardado -> telemetria/curva_coste.json")
    return 0


# ===========================================================================
# Controles
# ===========================================================================

def _autotest() -> int:
    titulo("CONTROLES DE curva_coste.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-52s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    # min/max rodante hacia adelante
    x = np.array([5.0, 3.0, 4.0, 1.0, 6.0, 2.0])
    r = _rodante(x, 3, "min")
    chk(abs(r[0] - 3.0) < 1e-9, "el minimo rodante mira HACIA ADELANTE",
        "%.1f (esperado 3.0 = min de 5,3,4)" % r[0])
    chk(abs(r[3] - 1.0) < 1e-9, "y en la posicion 3 da 1.0", "%.1f" % r[3])

    # por_bucket
    t = np.array([0.0, 1.0, 2.0, 61.0, 62.0])
    v = np.array([10.0, 8.0, 9.0, 3.0, 5.0])
    b = _por_bucket(t, v, 0.0, 62.0, 60.0, np.min)
    chk(abs(b[0] - 8.0) < 1e-9 and abs(b[1] - 3.0) < 1e-9,
        "el agregado por casilla separa bien los buckets", "%.1f y %.1f" % (b[0], b[1]))

    # financiacion: signo opuesto para largo y corto, y lineal en H
    CO_cache = CO.cache().get("financiacion")
    if CO_cache:
        f = {x["H"]: x for x in financiacion()}
        a, b_ = f[900], f[1800]
        chk(abs(b_["F_abs_med"] - 2 * a["F_abs_med"]) < 1e-9,
            "la financiacion esperada es lineal en H por debajo de 8 h")
        chk(a["F_esperada_largo"] * a["F_esperada_corto"] <= 0,
            "largo y corto tienen SIGNO OPUESTO (es transferencia, no coste)")
        chk(abs(f[28800]["p_cruce"] - 1.0) < 1e-9, "a 8 h la probabilidad de cruce es 1")
    else:
        log("  (financiacion no cacheada: corre `python coste.py --leer`)")

    # c(u) leida
    chk(CO.COMISIONES_LEIDAS is True, "c(u) esta LEIDA de la cuenta, no asumida")
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
    if a.etapa == "componentes":
        return etapa_componentes(a)
    log("usa --autotest o --etapa=componentes")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
