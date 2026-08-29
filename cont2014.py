# -*- coding: utf-8 -*-
"""
cont2014.py -- Tarea 1: calibracion contra Cont, Kukanov & Stoikov (2014).

    python cont2014.py --autotest    controles con verdad conocida
    python cont2014.py --etapa=tabla la Tabla 4 Panel A sobre nuestras capturas
    python cont2014.py --etapa=lambda  beta_i = c / AD_i^lambda

FUENTE: Cont, R., Kukanov, A. & Stoikov, S. (2014), *The Price Impact of Order
Book Events*, Journal of Financial Econometrics 12(1), 47-88. Sus cifras de
referencia, todas **CONTEMPORANEAS**, dt = 10 s, punto medio, 50 acciones NYSE,
abril 2010:

    R2 medio de  dP_k = alpha + beta*OFI_k + e_k        : 65 %
    idem excluyendo del OFI los eventos que mueven el precio : 35-60 %
    desequilibrio de transacciones solo                 : 32 %
    los dos regresores juntos                           : 67 %
      (el coeficiente de transacciones es significativo en solo el 31 % de submuestras)
    con termino cuadratico OFI*|OFI|                     : 65 % -> 68 %, insignificante
    beta_i = c / AD_i^lambda  con  lambda ~ 0.98, y lambda = 1 no rechazable en 35/50
    sobre precio de TRANSACCION en tiempo de operaciones : 14 % (L=2), 38 % (L=5), 51 % (L=10)

⚠ **INTERPRETACION DECLARADA ANTES DE CORRER.** Un `R2` sustancialmente por
debajo de 35-65 % **no es un hallazgo sobre BTCUSDT**: es una medida de la
degradacion de nuestro **OFI-L1 aproximado**. `@bookTicker` es un snapshot con
estrangulamiento, no el flujo diferencial de `@depth`, y la v3.2 §4.3 ya midio
que el residuo de reconciliacion entre variaciones de cantidad y transacciones
vistas es del **93.5 %**. Hay que reportar las dos cosas como UNA conclusion:
cuanto del deficit es instrumento y cuanto es mercado.

⚠ **CONTEMPORANEO, NO PREDICTIVO.** Todo lo de este modulo regresa `dP_k` contra
`OFI_k` del **mismo** intervalo. No es comparable con nada de la Adenda C ni con
el §1 de la v4.1, que predicen el futuro desde el pasado. Ver la convencion de
la Tarea 2 en `CLAUDE.md`.
"""

from __future__ import annotations

import argparse
import glob
import os

import numpy as np

import horizonte as H
from captura_larga import ofi_l1

log, titulo = H.log, H.titulo

DT_REJILLA = 10.0        # su dt
SUB_S = 1800.0           # submuestras de media hora, como el articulo
REF = {"ofi": 0.65, "ofi_sin_cambio": (0.35, 0.60), "ti": 0.32, "ambos": 0.67,
       "cuadratico": 0.68, "lambda": 0.98}


# ===========================================================================
# Ingesta: rejilla uniforme de 10 s
# ===========================================================================

def _acumular(bt, bb, bB, ba, bA, t0, dt, nb, acc):
    """OFI, OFI sin eventos que mueven el precio, profundidad y mid por casilla."""
    e = ofi_l1(bb, bB, ba, bA)                      # longitud L-1
    # nota 3 del articulo: control de tautologia. Se excluyen los eventos que
    # CAMBIAN el precio; quedan solo los que mueven cantidades a precio fijo.
    quieto = (bb[1:] == bb[:-1]) & (ba[1:] == ba[:-1])
    k = np.floor((bt[1:] - t0) / dt).astype(np.int64)
    ok = (k >= 0) & (k < nb)
    k, e, quieto = k[ok], e[ok], quieto[ok]
    acc["ofi"] += np.bincount(k, weights=e, minlength=nb)
    acc["ofi_q"] += np.bincount(k, weights=e * quieto, minlength=nb)
    prof = (bB + bA)[1:][ok]
    acc["prof"] += np.bincount(k, weights=prof, minlength=nb)
    acc["n"] += np.bincount(k, minlength=nb)
    # mid: se guarda el ULTIMO de cada casilla (cierre), que es lo que define dP
    mid = (0.5 * (bb + ba))[1:][ok]
    orden = np.argsort(k, kind="stable")
    kk, mm = k[orden], mid[orden]
    ult = np.r_[np.flatnonzero(np.diff(kk)), kk.size - 1] if kk.size else np.array([], int)
    acc["mid"][kk[ult]] = mm[ult]


def serie_10s_estacional(t0, t1) -> dict:
    import curvas_estacional as C
    nb = int((t1 - t0) / DT_REJILLA)
    acc = {c: np.zeros(nb) for c in ("ofi", "ofi_q", "prof", "n", "mid", "ti")}
    for f, ra, rb, _ in C._indice("libro_"):
        if rb < t0 or ra > t1:
            continue
        import pyarrow.parquet as pq
        tb = pq.read_table(f, columns=["t", "b", "B", "a", "A"])
        bt = tb["t"].to_numpy().astype(float)
        bb, bB = tb["b"].to_numpy().astype(float), tb["B"].to_numpy().astype(float)
        ba, bA = tb["a"].to_numpy().astype(float), tb["A"].to_numpy().astype(float)
        del tb
        m = (bb > 0) & (ba > 0) & (bB > 0) & (bA > 0) & (bt >= t0 - 1) & (bt <= t1 + 1)
        if m.sum() < 2:
            continue
        _acumular(bt[m], bb[m], bB[m], ba[m], bA[m], t0, DT_REJILLA, nb, acc)
    for f, ra, rb, _ in C._indice("trades_"):
        if rb < t0 or ra > t1:
            continue
        import pyarrow.parquet as pq
        tb = pq.read_table(f, columns=["t", "cant", "maker"])
        tt = tb["t"].to_numpy().astype(float)
        q = tb["cant"].to_numpy().astype(float)
        mk = tb["maker"].to_numpy().astype(bool)
        del tb
        k = np.floor((tt - t0) / DT_REJILLA).astype(np.int64)
        ok = (k >= 0) & (k < nb) & (q > 0)
        acc["ti"] += np.bincount(k[ok], weights=(np.where(mk[ok], -1.0, 1.0) * q[ok]),
                                 minlength=nb)
    return acc


def serie_10s_v33() -> dict:
    from captura_larga import cargar_larga, tramos_continuos
    d = cargar_larga("telemetria/captura_v33")
    bt, bb, bB = d["bk_t"], d["bk_b"], d["bk_B"]
    ba, bA = d["bk_a"], d["bk_A"]
    m = (bb > 0) & (ba > 0) & (bB > 0) & (bA > 0)
    bt, bb, bB, ba, bA = bt[m], bb[m], bB[m], ba[m], bA[m]
    o = np.argsort(bt, kind="stable")
    bt, bb, bB, ba, bA = bt[o], bb[o], bB[o], ba[o], bA[o]
    tt = d["tr_t"]
    q = d["tr_cant"]
    mk = d["tr_maker"].astype(bool)
    okt = np.isfinite(d["tr_precio"]) & (d["tr_precio"] > 0) & (q > 0)
    a, b = max(tramos_continuos(tt[okt]), key=lambda z: z[1] - z[0])
    t0, t1 = float(tt[okt][a]), float(tt[okt][b - 1])
    nb = int((t1 - t0) / DT_REJILLA)
    acc = {c: np.zeros(nb) for c in ("ofi", "ofi_q", "prof", "n", "mid", "ti")}
    sel = (bt >= t0 - 1) & (bt <= t1 + 1)
    _acumular(bt[sel], bb[sel], bB[sel], ba[sel], bA[sel], t0, DT_REJILLA, nb, acc)
    k = np.floor((tt[okt] - t0) / DT_REJILLA).astype(np.int64)
    ok = (k >= 0) & (k < nb)
    acc["ti"] += np.bincount(k[ok], weights=(np.where(mk[okt][ok], -1.0, 1.0) * q[okt][ok]),
                             minlength=nb)
    return acc


# ===========================================================================
# Regresion con errores estandar de White
# ===========================================================================

def ols_white(X, y):
    XtX = X.T @ X
    try:
        bi = np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        return None
    beta = bi @ (X.T @ y)
    res = y - X @ beta
    S = (X * res[:, None]).T @ (X * res[:, None])
    V = bi @ S @ bi
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - float(np.sum(res ** 2)) / sst if sst > 0 else float("nan")
    return {"beta": beta, "se": np.sqrt(np.maximum(np.diag(V), 0)), "r2": r2}


def por_submuestras(acc, etiqueta) -> dict:
    nb = acc["mid"].size
    mid = acc["mid"].copy()
    # casillas sin ninguna actualizacion de libro: se arrastra el ultimo mid
    vacia = mid <= 0
    idx = np.arange(nb)
    lleno = idx[~vacia]
    if lleno.size < 100:
        return {}
    mid = np.interp(idx, lleno, mid[lleno])
    dP = np.diff(mid)
    ofi, ofiq, ti = acc["ofi"][1:], acc["ofi_q"][1:], acc["ti"][1:]
    prof = np.where(acc["n"][1:] > 0, acc["prof"][1:] / np.maximum(acc["n"][1:], 1), np.nan)
    por = int(SUB_S / DT_REJILLA)
    res = {c: [] for c in ("ofi", "ofi_q", "ti", "ambos", "cuad")}
    sig_ti, betas, ADs = 0, [], []
    nsub = 0
    for a in range(0, dP.size - por, por):
        s = slice(a, a + por)
        y = dP[s]
        if np.std(y) <= 0 or np.std(ofi[s]) <= 0:
            continue
        uno = np.ones(por)
        r1 = ols_white(np.column_stack([uno, ofi[s]]), y)
        r2_ = ols_white(np.column_stack([uno, ofiq[s]]), y) if np.std(ofiq[s]) > 0 else None
        r3 = ols_white(np.column_stack([uno, ti[s]]), y) if np.std(ti[s]) > 0 else None
        r4 = ols_white(np.column_stack([uno, ofi[s], ti[s]]), y) if np.std(ti[s]) > 0 else None
        r5 = ols_white(np.column_stack([uno, ofi[s], ofi[s] * np.abs(ofi[s])]), y)
        if r1:
            res["ofi"].append(r1["r2"])
            betas.append(r1["beta"][1])
            ADs.append(np.nanmean(prof[s]))
        if r2_:
            res["ofi_q"].append(r2_["r2"])
        if r3:
            res["ti"].append(r3["r2"])
        if r4:
            res["ambos"].append(r4["r2"])
            if abs(r4["beta"][2]) > 1.96 * max(r4["se"][2], 1e-30):
                sig_ti += 1
        if r5:
            res["cuad"].append(r5["r2"])
        nsub += 1
    return {"etiqueta": etiqueta, "n_sub": nsub,
            "r2": {k: (float(np.mean(v)) if v else float("nan")) for k, v in res.items()},
            "frac_sig_ti": sig_ti / max(nsub, 1),
            "betas": np.array(betas), "AD": np.array(ADs)}


# ===========================================================================
# Etapas
# ===========================================================================

def _fuentes(cual=None):
    """[!] UNA FUENTE POR INVOCACION. Correr dos analisis pesados a la vez los
    mata entre ellos -- ya paso. El tramo de 131.83 h son ~130 M de snapshots de
    libro que hay que leer de parquet, y eso no convive con nada."""
    import curvas_estacional as C
    out = []
    if os.path.isdir("telemetria/captura_v33"):
        out.append(("captura_v33", serie_10s_v33))
    lim = C.limites_tramos()
    for k, (t0, t1, _) in enumerate(lim):
        out.append(("estacional_tramo%d" % k,
                    (lambda a, b: (lambda: serie_10s_estacional(a, b)))(t0, t1)))
    if cual:
        out = [x for x in out if x[0] == cual]
    return out


def etapa_tabla(args) -> int:
    titulo("TAREA 1 -- TABLA 4 PANEL A de Cont, Kukanov & Stoikov (2014)")
    log("")
    log("  ⚠ TODO ESTO ES CONTEMPORANEO: dP_k contra OFI_k del MISMO intervalo de")
    log("    10 s. No es comparable con la Adenda C ni con el §1, que son predictivos.")
    log("")
    log("  %-20s %6s | %8s %10s %8s %8s %8s | %9s"
        % ("serie", "n_sub", "OFI", "OFI sin dP", "transac", "ambos", "cuadrat", "sig(TI)"))
    log("  " + "-" * 92)
    RES = "telemetria/cont2014_resultados.json"
    import json
    guardado = {}
    prev = json.load(open(RES, encoding="utf-8")) if os.path.exists(RES) else {}
    for nom, r in prev.items():
        q = r["r2"]
        log("  %-20s %6d | %7.1f%% %9.1f%% %7.1f%% %7.1f%% %7.1f%% | %8.1f%%   (previo)"
            % (nom, r["n_sub"], 100*q["ofi"], 100*q["ofi_q"], 100*q["ti"],
               100*q["ambos"], 100*q["cuad"], 100*r["frac_sig_ti"]))
    for nom, fn in _fuentes(args.fuente):
        acc = fn()
        r = por_submuestras(acc, nom)
        if not r:
            log("  %-20s (insuficiente)" % nom)
            continue
        guardado[nom] = r
        prev[nom] = {"n_sub": r["n_sub"], "r2": r["r2"], "frac_sig_ti": r["frac_sig_ti"]}
        json.dump(prev, open(RES, "w", encoding="utf-8"), indent=1)
        q = r["r2"]
        log("  %-20s %6d | %7.1f%% %9.1f%% %7.1f%% %7.1f%% %7.1f%% | %8.1f%%"
            % (nom, r["n_sub"], 100 * q["ofi"], 100 * q["ofi_q"], 100 * q["ti"],
               100 * q["ambos"], 100 * q["cuad"], 100 * r["frac_sig_ti"]))
        del acc
    log("  %-20s %6s | %7.1f%% %9s %7.1f%% %7.1f%% %7.1f%% | %8.1f%%"
        % ("PUBLICADO (NYSE)", 50, 100 * REF["ofi"], "35-60%", 100 * REF["ti"],
           100 * REF["ambos"], 100 * REF["cuadratico"], 31.0))
    log("")
    log("  --- LA CONCLUSION, QUE ES UNA SOLA: instrumento contra mercado ---")
    for nom, r in prev.items():
        d = r["r2"]["ofi"]
        log("    %-20s R2(OFI) = %.1f %%  contra 65 %% publicado -> deficit %.0fx"
            % (nom, 100 * d, REF["ofi"] / max(d, 1e-9)))
    log("")
    log("    El OFI de este proyecto es **OFI-L1 APROXIMADO desde `bookTicker`**, no")
    log("    el flujo diferencial de `@depth`. La v3.2 §4.3 midio que el residuo de")
    log("    reconciliacion entre variaciones de cantidad y transacciones vistas es")
    log("    del 93.5 %: casi todo el movimiento del libro ocurre en actividad de")
    log("    limite INVISIBLE entre snapshots. Un deficit grande mide ESO, no el")
    log("    mercado, y asi hay que citarlo.")
    if guardado:
        viejo = dict(np.load("telemetria/cont2014.npz")) if os.path.exists(
            "telemetria/cont2014.npz") else {}
        viejo.update({("%s_%s" % (n, k)): v for n, r in guardado.items()
                      for k, v in (("betas", r["betas"]), ("AD", r["AD"]))})
        np.savez("telemetria/cont2014.npz", **viejo)
    return 0


def etapa_lambda(args) -> int:
    titulo("TAREA 1b -- beta_i = c / AD_i^lambda")
    log("")
    log("  El articulo: lambda ~ 0.98 de media, y lambda = 1 no rechazable en 35 de")
    log("  50 acciones. El coeficiente de impacto es inversamente proporcional a la")
    log("  profundidad, y de ahi sale su estacionalidad intradia.")
    log("")
    if not os.path.exists("telemetria/cont2014.npz"):
        log("  falta `--etapa=tabla` primero")
        return 2
    d = np.load("telemetria/cont2014.npz")
    for nom in sorted({k.rsplit("_", 1)[0] for k in d.files}):
        b = d["%s_betas" % nom]
        A = d["%s_AD" % nom]
        m = np.isfinite(b) & np.isfinite(A) & (b > 0) & (A > 0)
        if m.sum() < 10:
            log("  %-20s (insuficiente: %d submuestras utiles)" % (nom, int(m.sum())))
            continue
        X = np.column_stack([np.ones(int(m.sum())), np.log(A[m])])
        r = ols_white(X, np.log(b[m]))
        lam = -r["beta"][1]
        se = r["se"][1]
        z = (lam - 1.0) / max(se, 1e-30)
        log("  %-20s n=%4d   lambda = %+.4f +- %.4f   z(lambda=1) = %+.2f  -> %s"
            % (nom, int(m.sum()), lam, se, z,
               "compatible con 1" if abs(z) < 1.96 else "*** DIFIERE de 1 ***"))
    log("")
    log("  LECTURA DECLARADA: si lambda sale lejos de 1, la profundidad mas alla del")
    log("  primer nivel domina, y ESO SI seria argumento para ingerir `@depth`.")
    log("  Hasta ahora no habia ninguno.")
    return 0


# ===========================================================================
# Controles
# ===========================================================================

def _autotest() -> int:
    titulo("CONTROLES DE cont2014.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-52s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    rng = np.random.default_rng(7)
    n = 5000
    x = rng.normal(0, 1, n)
    y = 2.5 * x + 0.5 * rng.normal(0, 1, n)
    r = ols_white(np.column_stack([np.ones(n), x]), y)
    chk(abs(r["beta"][1] - 2.5) < 0.05, "OLS recupera la pendiente", "%.4f" % r["beta"][1])
    chk(abs(r["r2"] - (2.5 ** 2 / (2.5 ** 2 + 0.25))) < 0.02, "y el R2 teorico",
        "%.4f" % r["r2"])
    # White: con heterocedasticidad el se de White > el clasico
    y2 = 2.5 * x + np.abs(x) * rng.normal(0, 1, n)
    rw = ols_white(np.column_stack([np.ones(n), x]), y2)
    res = y2 - np.column_stack([np.ones(n), x]) @ rw["beta"]
    se_cl = float(np.sqrt(np.var(res) * np.linalg.inv(
        np.column_stack([np.ones(n), x]).T @ np.column_stack([np.ones(n), x]))[1, 1]))
    chk(rw["se"][1] > se_cl, "el se de White supera al clasico con heterocedasticidad",
        "%.5f contra %.5f" % (rw["se"][1], se_cl))
    # OFI: signo correcto con una subida pura del bid
    bb = np.array([100.0, 100.1]); bB = np.array([5.0, 7.0])
    ba = np.array([100.2, 100.2]); bA = np.array([5.0, 5.0])
    chk(ofi_l1(bb, bB, ba, bA)[0] > 0, "OFI positivo si el bid sube",
        "%.1f" % ofi_l1(bb, bB, ba, bA)[0])
    bb2 = np.array([100.0, 99.9])
    chk(ofi_l1(bb2, bB, ba, bA)[0] < 0, "OFI negativo si el bid baja",
        "%.1f" % ofi_l1(bb2, bB, ba, bA)[0])
    # lambda sintetica: beta = c/AD, se recupera 1
    AD = np.exp(rng.normal(0, 0.5, 300))
    be = 3.0 / AD ** 1.0 * np.exp(rng.normal(0, 0.05, 300))
    rr = ols_white(np.column_stack([np.ones(300), np.log(AD)]), np.log(be))
    chk(abs(-rr["beta"][1] - 1.0) < 0.05, "recupera lambda = 1 conocida",
        "%.4f" % (-rr["beta"][1]))
    log("")
    log("RESULTADO: %d fallo(s)" % fallos)
    return 1 if fallos else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", default="")
    ap.add_argument("--fuente", default="")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    return {"tabla": etapa_tabla, "lambda": etapa_lambda}.get(
        a.etapa, lambda _: (log("usa --autotest o --etapa=tabla|lambda"), 2)[1])(a)


if __name__ == "__main__":
    raise SystemExit(main())
