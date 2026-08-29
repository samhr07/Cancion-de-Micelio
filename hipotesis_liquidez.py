# -*- coding: utf-8 -*-
"""
hipotesis_liquidez.py -- H1 (presion de libro) y H2 (sistema perturbado).

    python hipotesis_liquidez.py --autotest    controles con verdad conocida
    python hipotesis_liquidez.py --etapa=libro cachea la serie binada del libro
    python hipotesis_liquidez.py --etapa=h1    volatilidad contra liquidez
    python hipotesis_liquidez.py --etapa=h2    volatilidad contra numero de perturbaciones

=========================================================================
PREDICCIONES DECLARADAS ANTES DE MIRAR NINGUNA CIFRA  (2026-08-23)
=========================================================================

H1 -- "la volatilidad es un reflejo de la presion del libro de liquidez".

  Se mide `D` = profundidad en el mejor nivel (`B + A`, en BTC) del bookTicker.
  El contraste es la regresion

        log sigma = a + b*log nu + c*log D

  y lo que decide es el SIGNO de `c`, controlando por el flujo.

  P1. **La microestructura estandar predice `c < 0`**: a flujo dado, un libro
      mas profundo absorbe la misma orden con menos desplazamiento, asi que mas
      profundidad va con MENOS volatilidad. Es la forma reducida de Kyle
      (lambda ~ 1/profundidad).
  P2. La intuicion del operador -- "sigma proporcional a medida que el libro se
      satura" -- implica `c > 0`. Las dos lecturas son incompatibles y el signo
      las separa. **Se declara que P1 es la esperada**; si sale `c > 0`, gana P2
      y hay que explicar el mecanismo.
  P3. `|c|` sera PEQUENA comparada con `b`. Solo hay nivel 1 (`bookTicker`), y
      la profundidad del mejor nivel es un proxy pobre de la liquidez total del
      libro. Un `|c|` grande seria sospechoso, no triunfal.

  [!] LIMITACION QUE NO SE PUEDE ARREGLAR CON ESTOS DATOS: `@bookTicker` da solo
      el TOPE del libro. "Saturacion del libro" en el sentido de ordenes en
      reposo a varios niveles NO es observable aqui. Lo que se contrasta es la
      cola del mejor precio, que en un mercado con el spread clavado en 1 tick
      -- y el de BTCUSDT lo esta, mediana y p90 = 0.1000 USD en los cuatro
      tramos -- es la variable de cola relevante, pero NO es el libro entero.

H2 -- "sistema perturbado": `nu` = numero de perturbaciones, `sigma` = respuesta.

  Si cada transaccion aporta una perturbacion de varianza independiente y
  constante, entonces sobre una ventana de duracion fija

        sigma^2  =  k * nu        =>        log sigma = const + 0.5*log nu

  P4. **La pendiente de `log sigma` contra `log nu` sera ~0.5** si la hipotesis
      es literal. Una pendiente por encima de 0.5 significa que las
      perturbaciones no son independientes (se agrupan y se refuerzan); por
      debajo, que se cancelan entre si.
  P5. **La varianza POR TRANSACCION sera mas estable que la varianza por
      segundo.** Es la prediccion que decide si se puede modelar UNA sola
      variable en vez de dos: si `sigma^2/nu` es casi constante, `nu` contiene a
      `sigma` y sobra una.
  P6. Si P5 se cumple, el residuo de `log sigma` tras quitar `log nu` **no
      tendra estructura diurna**: no quedara nada que la hora explique.

  [!] ESTO YA TIENE NOMBRE, Y HAY QUE DECIRLO ANTES DE MEDIRLO. La subordinacion
      del precio al reloj de transacciones es Clark (1973), y la Hipotesis de
      Mezcla de Distribuciones es Tauchen & Pitts (1983). Este proyecto YA
      tropezo con la version de 1994 (Jones, Kaul & Lipson) al medir `q̄`, y la
      registro como literatura, no como hallazgo. **Si H2 sale confirmada, lo
      que se ha hecho es reproducir un hecho estilizado de los mejor
      replicados de la microestructura** -- lo cual valida la tuberia y da un
      modelo con el que trabajar, pero no es un descubrimiento.
"""

from __future__ import annotations

import argparse
import datetime as DT
import glob
import os

import numpy as np
import pyarrow.parquet as pq

import curvas_estacional as C
import estacionalidad as E
import horizonte as H

log, titulo = H.log, H.titulo
CACHE_LIBRO = "telemetria/libro_binado.npz"


# ===========================================================================
# Serie de libro binada
# ===========================================================================

def construir_libro(bin_s: float = E.BIN_S) -> dict:
    """Profundidad y spread medios por casilla, en streaming sobre el libro."""
    ventanas = []
    for k in C._tramos_cache():
        m = np.load(C.CACHE % k)
        t = m["t"]
        ventanas.append((float(t[0]), float(t[-1])))
        del m
    acc = {}
    for f, ra, rb, _ in C._indice("libro_"):
        if not any(rb >= a - 1 and ra <= b + 1 for a, b in ventanas):
            continue
        try:
            tb = pq.read_table(f, columns=["t", "b", "B", "a", "A"])
        except Exception:
            continue
        bt = tb["t"].to_numpy().astype(np.float64)
        bb = tb["b"].to_numpy().astype(np.float64)
        qb = tb["B"].to_numpy().astype(np.float64)
        ba = tb["a"].to_numpy().astype(np.float64)
        qa = tb["A"].to_numpy().astype(np.float64)
        del tb
        ok = (np.isfinite(bb) & (bb > 0) & np.isfinite(ba) & (ba > 0)
              & np.isfinite(qb) & (qb > 0) & np.isfinite(qa) & (qa > 0))
        bt, bb, qb, ba, qa = bt[ok], bb[ok], qb[ok], ba[ok], qa[ok]
        if bt.size == 0:
            continue
        dentro = np.zeros(bt.size, bool)
        for a, b in ventanas:
            dentro |= (bt >= a) & (bt <= b)
        if not dentro.any():
            continue
        bt, bb, qb, ba, qa = bt[dentro], bb[dentro], qb[dentro], ba[dentro], qa[dentro]
        # [!] VECTORIZADO A PROPOSITO. La primera version recorria las filas con
        # un `for` de Python: son 223 M de actualizaciones de libro y habria
        # tardado horas. `bincount` sobre la celda hace lo mismo en un barrido.
        cel = np.floor(bt / bin_s).astype(np.int64)
        prof = qb + qa
        des = (qb - qa) / prof
        spr = (ba - bb) * 1e4 / (0.5 * (ba + bb))          # spread en pb
        base = int(cel.min())
        idx = cel - base
        m = int(idx.max()) + 1
        n_ = np.bincount(idx, minlength=m)
        sp = np.bincount(idx, weights=prof, minlength=m)
        sd = np.bincount(idx, weights=des, minlength=m)
        ss = np.bincount(idx, weights=spr, minlength=m)
        for j in np.flatnonzero(n_):
            r = acc.setdefault(base + int(j), [0, 0.0, 0.0, 0.0])
            r[0] += int(n_[j])
            r[1] += float(sp[j])
            r[2] += float(sd[j])
            r[3] += float(ss[j])
    cels = np.array(sorted(acc), dtype=np.int64)
    n = np.array([acc[int(c)][0] for c in cels], float)
    d = {"t": (cels + 0.5) * bin_s,
         "n_upd": n,
         "prof": np.array([acc[int(c)][1] for c in cels]) / n,
         "desq": np.array([acc[int(c)][2] for c in cels]) / n,
         "spread_pb": np.array([acc[int(c)][3] for c in cels]) / n}
    return d


def etapa_libro(args) -> int:
    titulo("SERIE DE LIBRO BINADA (nivel 1)")
    d = construir_libro()
    np.savez_compressed(CACHE_LIBRO, **d)
    log("  casillas: %d   cache: %s" % (d["t"].size, CACHE_LIBRO))
    log("  profundidad B+A [BTC] : p10 %.3f  MED %.3f  p90 %.3f"
        % tuple(np.percentile(d["prof"], [10, 50, 90])))
    log("  spread [pb]           : p10 %.4f  MED %.4f  p90 %.4f"
        % tuple(np.percentile(d["spread_pb"], [10, 50, 90])))
    log("  actualizaciones/casilla: MED %.0f" % np.median(d["n_upd"]))
    return 0


def unir() -> dict:
    """Casillas donde hay a la vez precio/flujo y libro."""
    a = E.serie_binada()
    b = dict(np.load(CACHE_LIBRO))
    ca = np.round(a["t"] / E.BIN_S).astype(np.int64)
    cb = np.round(b["t"] / E.BIN_S).astype(np.int64)
    com, ia, ib = np.intersect1d(ca, cb, return_indices=True)
    d = {k: v[ia] for k, v in a.items()}
    for k in ("prof", "desq", "spread_pb", "n_upd"):
        d[k] = b[k][ib]
    return d


# ===========================================================================
# Utilidades de regresion con nulo y fuera de muestra
# ===========================================================================

def _ols(X, y):
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    return c


def parcial(y, x, ctrl):
    """Correlacion parcial de `y` con `x` controlando por las columnas `ctrl`."""
    Z = np.column_stack([np.ones_like(y)] + list(ctrl))
    ry = y - Z @ _ols(Z, y)
    rx = x - Z @ _ols(Z, x)
    return float(np.corrcoef(ry, rx)[0, 1])


def barajar_en_dia(v, dias, rng):
    o = v.copy()
    for q in np.unique(dias):
        s = dias == q
        o[s] = rng.permutation(o[s])
    return o


# ===========================================================================
# H1
# ===========================================================================

def etapa_h1(args) -> int:
    titulo("H1 -- LA VOLATILIDAD CONTRA LA PRESION DEL LIBRO")
    d = unir()
    ls, ln = d["log_sigma"], d["log_nu"]
    lp = np.log(d["prof"])
    log("")
    log("  %d casillas de %.0f s con precio y libro" % (ls.size, E.BIN_S))
    log("")
    log("  --- perfil diurno de la PROFUNDIDAD (la pregunta '?cuando se satura?') ---")
    log("  hora UTC |  habil: prof[BTC]  spread[pb] |  finde: prof[BTC]  spread[pb]")
    for h in range(0, 24, 2):
        s = (d["hora"] >= h) & (d["hora"] < h + 2)
        for etq, mask in (("", s & (d["finde"] == 0)),):
            pa = np.median(d["prof"][mask]) if mask.sum() else np.nan
            sa = np.median(d["spread_pb"][mask]) if mask.sum() else np.nan
            mb = s & (d["finde"] == 1)
            pb_ = np.median(d["prof"][mb]) if mb.sum() else np.nan
            sb = np.median(d["spread_pb"][mb]) if mb.sum() else np.nan
            log("    %2d-%2d  |        %8.3f    %7.4f  |       %8.3f    %7.4f"
                % (h, h + 2, pa, sa, pb_, sb))
    imax = int(np.argmax([np.median(d["prof"][(d["hora"] >= h) & (d["hora"] < h + 1)])
                          if ((d["hora"] >= h) & (d["hora"] < h + 1)).sum() else -1
                          for h in range(24)]))
    imin = int(np.argmin([np.median(d["prof"][(d["hora"] >= h) & (d["hora"] < h + 1)])
                          if ((d["hora"] >= h) & (d["hora"] < h + 1)).sum() else 1e18
                          for h in range(24)]))
    log("    -> libro MAS profundo en UTC %d, MENOS profundo en UTC %d" % (imax, imin))

    log("")
    log("  --- el contraste: log sigma = a + b*log nu + c*log D ---")
    X = np.column_stack([np.ones_like(ls), ln, lp])
    c = _ols(X, ls)
    res = ls - X @ c
    r2 = 1 - np.var(res) / np.var(ls)
    log("    b (log nu) = %+.4f     c (log D) = %+.4f     R2 dentro de muestra = %.4f"
        % (c[1], c[2], r2))
    log("    corr parcial de log sigma con log D, controlando log nu = %+.4f"
        % parcial(ls, lp, [ln]))

    rng = np.random.default_rng(5)
    nul = [parcial(ls, barajar_en_dia(lp, d["dia"], rng), [ln]) for _ in range(200)]
    q = np.percentile(nul, [2.5, 97.5])
    log("    NULO (profundidad barajada dentro de cada dia): IC95 [%+.4f, %+.4f]" % tuple(q))
    log("")
    if c[2] < 0 and parcial(ls, lp, [ln]) < q[0]:
        log("    *** P1 CONFIRMADA: c < 0 y fuera del nulo. Mas libro -> menos")
        log("        volatilidad a flujo dado. Es la forma reducida de Kyle. ***")
    elif c[2] > 0 and parcial(ls, lp, [ln]) > q[1]:
        log("    *** P2 GANA: c > 0 y fuera del nulo. La intuicion del operador se")
        log("        sostiene y P1 (microestructura estandar) queda refutada aqui. ***")
    else:
        log("    *** NI P1 NI P2: la parcial cae dentro del nulo. La profundidad de")
        log("        nivel 1 no aporta sobre el flujo. ***")
    log("")
    log("  [!] Solo hay NIVEL 1. 'Saturacion del libro' a varios niveles NO es")
    log("      observable con `@bookTicker`; esto contrasta la COLA DEL MEJOR PRECIO.")
    return 0


# ===========================================================================
# H2
# ===========================================================================

def etapa_h2(args) -> int:
    titulo("H2 -- SISTEMA PERTURBADO: ?basta UNA variable?")
    d = E.serie_binada()
    ls, ln, dias = d["log_sigma"], d["log_nu"], d["dia"]
    log("")
    log("  %d casillas de %.0f s sobre %d dias" % (ls.size, E.BIN_S, np.unique(dias).size))

    log("")
    log("  --- P4: pendiente de log sigma contra log nu (0.5 = perturbaciones independientes) ---")
    X = np.column_stack([np.ones_like(ln), ln])
    c = _ols(X, ls)
    r = float(np.corrcoef(ln, ls)[0, 1])
    log("    TODO           pendiente = %+.4f    corr = %+.4f" % (c[1], r))
    for etq, y_, x_ in (("ENTRE dias", None, None), ("DENTRO del dia", None, None)):
        pass
    med_s = np.array([ls[dias == q].mean() for q in np.unique(dias)])
    med_n = np.array([ln[dias == q].mean() for q in np.unique(dias)])
    ce = _ols(np.column_stack([np.ones_like(med_n), med_n]), med_s)
    log("    ENTRE dias     pendiente = %+.4f    corr = %+.4f  (n = %d dias)"
        % (ce[1], float(np.corrcoef(med_n, med_s)[0, 1]), med_n.size))
    ws, wn = ls.copy(), ln.copy()
    for q in np.unique(dias):
        s = dias == q
        ws[s] -= ws[s].mean()
        wn[s] -= wn[s].mean()
    cw = _ols(np.column_stack([np.ones_like(wn), wn]), ws)
    log("    DENTRO del dia pendiente = %+.4f    corr = %+.4f"
        % (cw[1], float(np.corrcoef(wn, ws)[0, 1])))

    log("")
    log("  --- P5: ?es la varianza POR TRANSACCION mas estable que la varianza por segundo? ---")
    # sigma es por sub-paso de SUB_S s; nu en tx/s. Varianza por transaccion:
    #   sigma^2 / (nu * SUB_S)
    lv_seg = 2 * ls                              # log de la varianza por sub-paso
    lv_tx = 2 * ls - (ln + np.log(E.SUB_S))      # log de la varianza por transaccion
    for etq, v in (("varianza por segundo   ", lv_seg),
                   ("varianza por transaccion", lv_tx)):
        ent = float(np.var([v[dias == q].mean() for q in np.unique(dias)]))
        w = v.copy()
        for q in np.unique(dias):
            w[dias == q] -= w[dias == q].mean()
        log("    %s  var TOTAL = %.4f   entre dias = %.4f   dentro = %.4f"
            % (etq, float(np.var(v)), ent, float(np.var(w))))
    red = 1 - np.var(lv_tx) / np.var(lv_seg)
    log("    -> normalizar por nu reduce la varianza de log(varianza) en %+.1f %%" % (100 * red))
    if red > 0.5:
        log("    *** P5 CONFIRMADA: nu contiene a sigma. Sobra una variable. ***")
    elif red > 0.0:
        log("    *** P5 PARCIAL: nu explica parte, no basta con una. ***")
    else:
        log("    *** P5 REFUTADA: normalizar por nu EMPEORA la estabilidad. ***")

    log("")
    log("  --- P6: ?queda estructura diurna en el residuo de log sigma dado log nu? ---")
    Xn = np.column_stack([np.ones_like(ln), ln])
    resid = ls - Xn @ _ols(Xn, ls)
    dr = dict(d)
    dr["resid"] = resid
    r2_res = E.r2_oos_por_dias(E.diseno(d, E.K_ARM, True), resid, dias)
    r2_sig = E.r2_oos_por_dias(E.diseno(d, E.K_ARM, True), ls, dias)
    log("    R2 fuera de muestra del modelo diurno sobre log sigma  : %+.4f" % r2_sig)
    log("    R2 fuera de muestra del modelo diurno sobre el RESIDUO : %+.4f" % r2_res)

    log("")
    log("  --- lo que decide en la practica: ?predice nu la sigma de la casilla SIGUIENTE? ---")
    o = np.argsort(d["t"])
    t_, ls_, ln_, di_ = d["t"][o], ls[o], ln[o], dias[o]
    cons = np.isclose(np.diff(t_), E.BIN_S)
    i0 = np.flatnonzero(cons)
    y = ls_[i0 + 1]
    for etq, X_ in (("constante            ", np.ones((i0.size, 1))),
                    ("log nu(t)            ", np.column_stack([np.ones(i0.size), ln_[i0]])),
                    ("log sigma(t)         ", np.column_stack([np.ones(i0.size), ls_[i0]])),
                    ("log nu(t) + log sigma(t)",
                     np.column_stack([np.ones(i0.size), ln_[i0], ls_[i0]]))):
        log("    %-24s R2 fuera de muestra = %+.4f"
            % (etq, E.r2_oos_por_dias(X_, y, di_[i0])))
    rng = np.random.default_rng(11)
    Xn_ = np.column_stack([np.ones(i0.size), barajar_en_dia(ln_[i0], di_[i0], rng)])
    log("    %-24s R2 fuera de muestra = %+.4f  <- NULO"
        % ("log nu barajado en el dia", E.r2_oos_por_dias(Xn_, y, di_[i0])))
    return 0


# ===========================================================================
# Controles
# ===========================================================================

def _autotest() -> int:
    titulo("CONTROLES DE hipotesis_liquidez.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-56s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    rng = np.random.default_rng(2)
    n = 4000
    dias = np.repeat(np.arange(20), n // 20)
    ln = rng.normal(0, 1, n)
    lp = rng.normal(0, 1, n)

    # verdad: log sigma = 0.5 log nu - 0.3 log D + ruido
    ls = 0.5 * ln - 0.3 * lp + 0.1 * rng.normal(size=n)
    X = np.column_stack([np.ones(n), ln, lp])
    c = _ols(X, ls)
    chk(abs(c[1] - 0.5) < 0.02, "recupera la pendiente 0.5 de log nu", "%+.4f" % c[1])
    chk(abs(c[2] + 0.3) < 0.02, "recupera el coeficiente -0.3 de log D", "%+.4f" % c[2])

    pc = parcial(ls, lp, [ln])
    chk(pc < -0.5, "la parcial con log D sale negativa y grande", "%+.4f" % pc)

    # nulo: barajar D dentro del dia debe anular la parcial
    nul = [parcial(ls, barajar_en_dia(lp, dias, rng), [ln]) for _ in range(100)]
    chk(abs(np.mean(nul)) < 0.05, "el nulo de barajado anula la parcial",
        "media %+.4f" % np.mean(nul))
    chk(pc < np.percentile(nul, 2.5), "la parcial real cae fuera del nulo")

    # P5 sintetica: si sigma^2 = k*nu exacto, normalizar debe reducir la varianza
    ls2 = 0.5 * ln + 0.02 * rng.normal(size=n)
    lv_seg = 2 * ls2
    lv_tx = 2 * ls2 - ln
    red = 1 - np.var(lv_tx) / np.var(lv_seg)
    chk(red > 0.9, "con sigma^2 = k*nu exacto, normalizar reduce >90 % la varianza",
        "%.1f %%" % (100 * red))

    # y si son independientes, normalizar debe EMPEORAR
    ls3 = 0.02 * rng.normal(size=n)
    red2 = 1 - np.var(2 * ls3 - ln) / np.var(2 * ls3)
    chk(red2 < 0, "con sigma independiente de nu, normalizar EMPEORA", "%.1f %%" % (100 * red2))

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
    if a.etapa == "libro":
        return etapa_libro(a)
    if a.etapa == "h1":
        return etapa_h1(a)
    if a.etapa == "h2":
        return etapa_h2(a)
    log("usa --autotest, --etapa=libro, =h1 o =h2")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
