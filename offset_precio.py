# -*- coding: utf-8 -*-
"""
offset_precio.py -- el OFFSET `P_ref`: cuanto del nivel de precio lo carga el flujo.

    python offset_precio.py --autotest
    python offset_precio.py --etapa=offset [--fuente=estacional] [--zona=ny]

LA PREGUNTA DEL OPERADOR (2026-09-05), en sus terminos: si el flujo de mercado
mueve el precio, entonces hay una parte del precio que el flujo explica y un
resto -- el OFFSET -- que no. Si ese offset **no cambia mucho entre dias**, el
nivel del activo lo sostiene el volumen; si cambia tanto como el precio, el
flujo no sostiene nada.

⚠ LA FORMULA LITERAL NO CIERRA DIMENSIONALMENTE, y conviene decirlo antes de
nada. `P_ref = S - Q_neto` resta BTC (flujo neto) a USD/BTC (precio), y el paso
intermedio `P/Q_neto = P_ref + P_var/Q_neto = 1 + P_ref/Q_neto` no es una
identidad. Lo que SI cierra, y es exactamente la misma idea, sale de integrar la
premisa del propio operador (`dP/dt` proporcional al flujo neto):

    P(t) = P_ref(t) + lambda * CumQ(t)          CumQ = suma de Q_neto hasta t

con `lambda` en [USD/BTC por BTC] -- el coeficiente de impacto permanente, la
`lambda` de Kyle. De ahi

    P_ref(t) = P(t) - lambda * CumQ(t)          <- EL OFFSET

que ya es un precio, en USD/BTC, y se puede comparar con `P` dia a dia. La
pregunta "cuanto cambia P_ref entre dias" pasa a tener respuesta numerica.

⚠ `lambda` SE ESTIMA SOBRE INCREMENTOS, NUNCA SOBRE NIVELES. `P` y `CumQ` son
las dos series integradas: regresar una sobre otra es el caso de manual de
regresion espuria, y este proyecto ya se lo encontro (2026-08-08 f: "la
correlacion con el NIVEL de precio es espuria por construccion, dos series no
estacionarias"). Se ajusta `dP_k = lambda * Q_k` sobre la rejilla, que es
estacionaria, y `P_ref` se DERIVA de ahi. No es un ajuste sobre el nivel.

⚠ SIN INTERCEPTO, Y ES UNA DECISION. Un intercepto en la regresion de
incrementos es una deriva constante por casilla: se comeria exactamente la
tendencia que estamos intentando atribuir, y `P_ref` saldria plano por
construccion. Es el mismo modo de fallo que la v4.2 Sec.4 tuvo que corregir tres
veces (la deriva colandose por cuatro puertas). Se ajusta sin intercepto y se
reporta al lado cuanto habria valido, como diagnostico.

⚠ EL NULO ES OBLIGATORIO. `P - lambda*CumQ` es una diferencia de dos series
integradas: con `Q` ROTADO circularmente sigue saliendo un `P_ref`, y sigue
teniendo un recorrido. La pregunta no es "cuanto se mueve P_ref" sino "cuanto
menos se mueve que el que saldria por azar". Por eso todo se reporta contra
rotacion circular de `Q` -- rotacion y no barajado, porque `Q` tiene memoria
larga (`N_eff` de los signos = 46.5) y barajarla daria un suelo demasiado
estrecho (la Adenda C midio ese sesgo en 75x a 292x).
"""

from __future__ import annotations

import argparse
import io
import json
import os

import numpy as np

import flujo_omega as F
import horizonte as H

log, titulo = H.log, H.titulo
N_SORTEOS = 200


# ===========================================================================
# 1. lambda y el offset
# ===========================================================================

def ols_sin_intercepto(x, y):
    """`y = b*x`. Devuelve `b`, `R2`, y el intercepto que HABRIA salido."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 8 or not np.any(x != 0):
        return {"b": np.nan, "R2": np.nan, "a_diag": np.nan, "n": int(x.size)}
    b = float(np.dot(x, y) / np.dot(x, x))
    res = y - b * x
    # R2 contra el modelo nulo `y = 0`, que es el correcto sin intercepto:
    # la pregunta es que fraccion de la varianza del incremento explica el flujo.
    sst = float(np.dot(y, y))
    r2 = 1.0 - float(np.dot(res, res)) / sst if sst > 0 else np.nan
    vx = float(np.var(x))
    a = (float(np.mean(y) - b * np.mean(x))) if vx > 0 else np.nan
    return {"b": b, "R2": r2, "a_diag": a, "n": int(x.size)}


def rachas(valida):
    """Indices agrupados en tramos de casillas validas CONSECUTIVAS."""
    idx = np.flatnonzero(np.asarray(valida, bool))
    if idx.size < 3:
        return []
    corte = np.flatnonzero(np.diff(idx) != 1)
    return [g for g in np.split(idx, corte + 1) if g.size >= 3]


def incrementos(r, regresor="Q"):
    """(`dP`, `x`) sobre las rachas validas. `x` es `Q_neto` o `Q_neto*nu`."""
    dP, X, POS = [], [], []
    for g in rachas(r["valida"]):
        p = r["precio"][g]
        d = np.diff(p)
        q = r["q_neto"][g][1:]
        if regresor == "Qnu":
            q = q * r["nu"][g][1:]
        m = np.isfinite(d) & np.isfinite(q)
        dP.append(d[m]); X.append(q[m]); POS.append(g[1:][m])
    if not dP:
        return np.empty(0), np.empty(0), np.empty(0, np.int64)
    return (np.concatenate(dP), np.concatenate(X), np.concatenate(POS))


def offset(r, lam, regresor="Q"):
    """`P_ref = P - lambda*CumQ` por racha, cada una anclada en su propio inicio.

    ⚠ El anclaje por racha no es cosmetico: entre dos rachas hay un hueco de la
    captura, y arrastrar `CumQ` a traves de el sumaria un flujo que nunca se
    observo. Cada racha arranca con `CumQ = 0`, asi que `P_ref` de cada racha
    empieza en el precio con el que la racha empieza y lo que se lee es su
    DERIVA dentro de la racha, no su nivel absoluto entre rachas.
    """
    out = np.full(r["precio"].size, np.nan)
    for g in rachas(r["valida"]):
        q = r["q_neto"][g]
        if regresor == "Qnu":
            q = q * r["nu"][g]
        cum = np.concatenate([[0.0], np.cumsum(q[1:])])
        out[g] = r["precio"][g] - lam * cum
    return out


def _rota(x, k):
    return np.roll(x, k)


def suelo_offset(r, regresor="Q", n=N_SORTEOS, semilla=31):
    """Distribucion de `lambda`, `R2` y recorrido de `P_ref` bajo `Q` ROTADO."""
    rng = np.random.default_rng(semilla)
    dP, X, _ = incrementos(r, regresor)
    m = X.size
    if m < 50:
        return {}
    lam, r2 = np.empty(n), np.empty(n)
    for i in range(n):
        k = int(rng.integers(max(3, m // 20), m - max(3, m // 20)))
        f = ols_sin_intercepto(_rota(X, k), dP)
        lam[i], r2[i] = f["b"], f["R2"]
    return {"lam_q95": float(np.nanpercentile(np.abs(lam), 95)),
            "R2_q95": float(np.nanpercentile(r2, 95)),
            "R2_med": float(np.nanmedian(r2)),
            "lam_med": float(np.nanmedian(lam))}


# ===========================================================================
# 2. Por dia
# ===========================================================================

def offset_en(r, g, lam, regresor="Q"):
    """`P_ref` de UNA racha, anclada en su inicio."""
    q = r["q_neto"][g]
    if regresor == "Qnu":
        q = q * r["nu"][g]
    cum = np.concatenate([[0.0], np.cumsum(q[1:])])
    return r["precio"][g] - lam * cum


def tabla_por_dia(r, g, lam, regresor="Q", zona="ny"):
    """Una fila por dia DENTRO de una racha continua.

    ⚠ POR RACHA Y NO SOBRE TODA LA REJILLA, y esto fue un fallo propio que la
    corrida de verificacion destapo. `P_ref` se ancla al inicio de cada racha
    (arrastrar `CumQ` por un hueco sumaria un flujo que nunca se observo), asi
    que un dia que CONTIENE una frontera de racha lleva dos anclas distintas y
    su `P_ref` da un salto de nivel artificial. En la corrida sintetica ese dia
    marcaba 460 pb de recorrido contra 3 pb de los demas, y peor: el recorrido
    ENTRE dias acababa midiendo otra vez el precio, porque cada racha arranca su
    `P_ref` en el precio con el que arranca la racha. Dentro de una racha el
    ancla es UNA, y ahi la comparacion entre dias si significa algo.
    """
    P = r["precio"][g]
    PR = offset_en(r, g, lam, regresor)
    t = r["t"][g]
    cal = F.calendario_ny(t) if zona == "ny" else None
    dia = cal["dia"] if cal is not None else np.floor(t / 86400.0).astype(np.int64)
    finde = cal["finde"] if cal is not None else np.zeros(dia.size, bool)
    Q = r["q_neto"][g]
    filas = []
    for d in np.unique(dia):
        m = dia == d
        if m.sum() < 10:
            continue
        pm = np.median(P[m])
        filas.append({
            "dia": int(d),
            "fecha": (F._fecha_ny(d) if cal is not None else F._fecha(d)),
            "finde": int(bool(finde[m][0])), "n": int(m.sum()),
            "P_med": float(pm),
            "P_rango_pb": float(1e4 * (np.max(P[m]) - np.min(P[m])) / pm),
            "Pref_med": float(np.median(PR[m])),
            "Pref_rango_pb": float(1e4 * (np.max(PR[m]) - np.min(PR[m])) / pm),
            "q_neto": float(np.sum(Q[m])),
            "sd_dP": float(np.std(np.diff(P[m]))),
            "sd_dPref": float(np.std(np.diff(PR[m]))),
        })
    return filas


def recorridos(filas):
    """Recorrido del precio y del offset ENTRE dias, los dos en pb."""
    if len(filas) < 2:
        return np.nan, np.nan
    Pm = np.array([x["P_med"] for x in filas], float)
    pm = np.array([x["Pref_med"] for x in filas], float)
    ref = float(np.median(Pm))
    return (1e4 * (np.max(Pm) - np.min(Pm)) / ref,
            1e4 * (np.max(pm) - np.min(pm)) / ref)


def nulo_recorrido(r, g, lam, regresor, zona, n=60, semilla=97):
    """Recorrido del offset entre dias con `Q` ROTADO dentro de la racha."""
    rng = np.random.default_rng(semilla)
    out = []
    m = g.size
    if m < 400:
        return np.array([])
    r2 = {k: r[k] for k in ("t", "precio", "q_neto", "nu", "valida")}
    for _ in range(n):
        k = int(rng.integers(m // 20, m - m // 20))
        q = r["q_neto"].copy()
        q[g] = np.roll(r["q_neto"][g], k)
        r2["q_neto"] = q
        fn = tabla_por_dia(r2, g, lam, regresor, zona)
        if len(fn) >= 2:
            out.append(recorridos(fn)[1])
    return np.array(out, float)


# ===========================================================================
# 3. Etapa
# ===========================================================================

def etapa_offset(args) -> int:
    r = F.cargar_rejilla(args.fuente)
    titulo("EL OFFSET  P_ref = P - lambda*CumQ   -- fuente `%s`" % args.fuente)
    log("  casillas validas: %d de %d" % (r["valida"].sum(), r["valida"].size))
    log("")
    log("  --- lambda, ajustada sobre INCREMENTOS y sin intercepto ---")
    log("  %-12s %14s %10s %10s %12s %10s"
        % ("regresor", "lambda", "R2", "R2 nulo", "n", "interc."))
    ajustes = {}
    for reg, nom in (("Q", "Q_neto"), ("Qnu", "Q_neto*nu")):
        dP, X, _ = incrementos(r, reg)
        f = ols_sin_intercepto(X, dP)
        su = suelo_offset(r, reg, n=args.sorteos)
        ajustes[reg] = (f, su)
        log("  %-12s %14.6e %10.5f %10.5f %12d %10.3e"
            % (nom, f["b"], f["R2"], su.get("R2_q95", np.nan), f["n"], f["a_diag"]))
    log("")
    log("  [!] El `R2` es CONTEMPORANEO (convencion del proyecto): explica el")
    log("      movimiento del MISMO intervalo, no predice el siguiente. No es")
    log("      comparable con nada de la Adenda C ni del Sec.1 de la v4.1.")
    log("")
    mejor = max(ajustes, key=lambda k: (ajustes[k][0]["R2"]
                                        if np.isfinite(ajustes[k][0]["R2"]) else -9))
    f, su = ajustes[mejor]
    lam = f["b"]
    log("  regresor con mas R2: %s   ->  lambda = %.6e USD/BTC por %s"
        % (mejor, lam, "BTC" if mejor == "Q" else "BTC/s"))
    if np.isfinite(su.get("R2_q95", np.nan)) and f["R2"] <= su["R2_q95"]:
        log("  ⚠ El R2 medido NO supera su suelo de rotacion: el offset que sigue")
        log("    NO es legible, cualquier serie rotada da lo mismo.")

    titulo("VARIANZA DEL PRECIO ABSORBIDA POR EL FLUJO")
    gs = rachas(r["valida"])
    dPR = np.concatenate([np.diff(offset_en(r, g, lam, mejor)) for g in gs])
    dPP = np.concatenate([np.diff(r["precio"][g]) for g in gs])
    vr = float(np.var(dPR) / np.var(dPP)) if np.var(dPP) > 0 else np.nan
    log("  sd(dP)     = %.6f USD/BTC" % np.std(dPP))
    log("  sd(dP_ref) = %.6f USD/BTC" % np.std(dPR))
    log("  razon de varianzas Var(dP_ref)/Var(dP) = %.5f   (1 - R2 = %.5f)"
        % (vr, 1.0 - f["R2"]))
    log("")
    log("  Por debajo de 1 el flujo absorbe varianza del precio; en 1 no absorbe")
    log("  nada; POR ENCIMA de 1 restarle el flujo empeora, o sea que `lambda` no")
    log("  describe la relacion.")

    # --- por RACHA continua ------------------------------------------------
    dur = [(g, (r["t"][g[-1]] - r["t"][g[0]]) / 86400.0) for g in gs]
    dur = sorted([x for x in dur if x[1] >= args.min_dias], key=lambda z: -z[1])
    titulo("EL OFFSET DIA A DIA -- por RACHA CONTINUA, reloj de %s"
           % ("NUEVA YORK" if args.zona == "ny" else "UTC"))
    log("  rachas continuas de >= %.1f dias: %d de %d"
        % (args.min_dias, len(dur), len(gs)))
    if not dur:
        log("")
        log("  ⚠ Ninguna racha llega a %.1f dias. La pregunta 'cuanto cambia el"
            % args.min_dias)
        log("    offset entre dias' NO se puede responder: `P_ref` se ancla al")
        log("    inicio de cada racha, asi que comparar dias de rachas distintas")
        log("    compara dos anclas arbitrarias, no dos offsets.")
        return 2
    resumen = []
    for ki, (g, d) in enumerate(dur):
        filas = tabla_por_dia(r, g, lam, mejor, args.zona)
        if len(filas) < 2:
            continue
        log("")
        log("  --- racha %d: %.2f dias, %d casillas, %d dias con dato ---"
            % (ki, d, g.size, len(filas)))
        log("  %-16s %3s %6s %12s %12s %12s %12s %11s"
            % ("fecha", "fin", "n", "P mediana", "P rango pb", "Pref rango pb",
               "q_neto", "sd dPref/dP"))
        for x in filas:
            log("  %-16s %3s %6d %12.2f %12.1f %12.1f %+12.2f %11.4f"
                % (x["fecha"], "SI" if x["finde"] else "-", x["n"], x["P_med"],
                   x["P_rango_pb"], x["Pref_rango_pb"], x["q_neto"],
                   x["sd_dPref"] / max(x["sd_dP"], 1e-12)))
        rp, rr = recorridos(filas)
        nul = nulo_recorrido(r, g, lam, mejor, args.zona,
                             n=min(args.sorteos, 60))
        log("")
        log("  recorrido ENTRE dias:  precio %8.1f pb   offset %8.1f pb   razon %.4f"
            % (rp, rr, rr / max(rp, 1e-12)))
        if nul.size:
            log("  NULO (Q rotado dentro de la racha): offset p05 %.1f  MED %.1f"
                "  p95 %.1f pb" % tuple(np.percentile(nul, [5, 50, 95])))
            log("  razon medido/nulo mediano: %.4f" % (rr / max(np.median(nul), 1e-12)))
        resumen.append({"racha": ki, "dias": d, "n_dias": len(filas),
                        "rango_P_pb": rp, "rango_Pref_pb": rr,
                        "razon": rr / max(rp, 1e-12),
                        "nulo_med": float(np.median(nul)) if nul.size else np.nan,
                        "filas": filas})

    titulo("CUANTO CAMBIA EL OFFSET ENTRE DIAS -- la pregunta del operador")
    if not resumen:
        log("  ninguna racha con al menos dos dias.")
        return 2
    log("  %-8s %8s %8s %14s %14s %10s %14s"
        % ("racha", "dias", "n_dias", "P rango pb", "Pref rango pb", "razon",
           "nulo MED pb"))
    for x in resumen:
        log("  %-8d %8.2f %8d %14.1f %14.1f %10.4f %14.1f"
            % (x["racha"], x["dias"], x["n_dias"], x["rango_P_pb"],
               x["rango_Pref_pb"], x["razon"], x["nulo_med"]))
    log("")
    log("  LECTURA. `razon` cerca de 0: el nivel del activo lo carga el flujo")
    log("  acumulado, y el offset es casi una constante. Cerca de 1 o mayor: el")
    log("  flujo no sostiene el nivel. Y lo que decide no es `razon` sino la")
    log("  comparacion contra el NULO -- `P - lambda*CumQ` es una diferencia de")
    log("  dos series integradas y SIEMPRE tiene recorrido, tambien por azar.")

    # --- estabilidad de lambda dia a dia -----------------------------------
    titulo("ESTABILIDAD DE lambda -- si el coeficiente no se sostiene, la"
           " descomposicion tampoco")
    t = r["t"]
    cal = F.calendario_ny(t) if args.zona == "ny" else None
    dia_c = cal["dia"] if cal is not None else np.floor(t / 86400.0).astype(np.int64)
    log("")
    log("  %-16s %3s %14s %10s %10s" % ("fecha", "fin", "lambda", "R2", "n"))
    lams = []
    for d in np.unique(dia_c[r["valida"]]):
        sel = {k: r[k] for k in ("t", "precio", "q_neto", "nu", "valida")}
        sel["valida"] = r["valida"] & (dia_c == d)
        dPd, Xd, _ = incrementos(sel, mejor)
        if Xd.size < 100:
            continue
        fd = ols_sin_intercepto(Xd, dPd)
        fin = bool(cal["finde"][dia_c == d][0]) if cal is not None else False
        lams.append(fd["b"])
        log("  %-16s %3s %14.6e %10.5f %10d"
            % (F._fecha_ny(d) if cal is not None else F._fecha(d),
               "SI" if fin else "-", fd["b"], fd["R2"], fd["n"]))
    if len(lams) > 1:
        lams = np.array(lams, float)
        log("")
        log("  lambda entre dias: min %.4e  MED %.4e  max %.4e   recorrido %.2fx"
            % (np.min(lams), np.median(lams), np.max(lams),
               np.max(lams) / max(np.min(lams), 1e-30)))
        log("  [!] Si `lambda` no se sostiene entre dias, `P_ref` calculada con la")
        log("      `lambda` global no es un offset: es el residuo de un modelo que")
        log("      no aplica a todos los dias por igual.")

    json.dump({"fuente": args.fuente, "regresor": mejor, "lambda": lam,
               "R2_contemporaneo": f["R2"], "R2_nulo_q95": su.get("R2_q95"),
               "var_ratio": vr,
               "rachas": [{k: v for k, v in x.items() if k != "filas"}
                          for x in resumen]},
              io.open(F.ruta_salida("offset_%s.json" % args.fuente), "w",
                      encoding="ascii"))
    csv = F.ruta_salida("offset_dias_%s.csv" % args.fuente)
    cab = ["racha", "fecha", "dia", "finde", "n", "P_med", "P_rango_pb", "Pref_med",
           "Pref_rango_pb", "q_neto", "sd_dP", "sd_dPref"]
    with io.open(csv, "w", encoding="ascii") as fh:
        fh.write(",".join(cab) + "\n")
        for x in resumen:
            for y in x["filas"]:
                y = dict(y, racha=x["racha"])
                fh.write(",".join(("%s" % y[c]) if c == "fecha" else ("%.10g" % y[c])
                                  for c in cab) + "\n")
    log("")
    log("  tabla por dia en %s" % csv)
    return 0


# ===========================================================================
# 4. Controles con verdad conocida
# ===========================================================================

def _rejilla_sintetica(n, lam_true, ruido, semilla=3, deriva=0.0):
    """Rejilla con `P = P0 + lam_true*CumQ + ruido + deriva`. Verdad conocida."""
    rng = np.random.default_rng(semilla)
    q = rng.normal(0, 2.0, n)
    dP = lam_true * q + rng.normal(0, ruido, n) + deriva
    P = 63000.0 + np.cumsum(dP)
    return {"t": np.arange(n) * 10.0, "precio": P, "q_neto": q,
            "nu": np.full(n, 20.0), "valida": np.ones(n, bool)}


def _autotest() -> int:
    titulo("offset_precio.py -- CONTROLES")
    n_ok = n_tot = 0

    def chk(ok, msg, det=""):
        nonlocal n_ok, n_tot
        n_tot += 1
        n_ok += bool(ok)
        log("  [%s] %s%s" % ("OK " if ok else "FALLA", msg,
                             ("   %s" % det) if det else ""))

    # --- 1. OLS sin intercepto, exacto sobre un caso sin ruido -------------
    x = np.linspace(-5, 5, 200)
    f = ols_sin_intercepto(x, 3.7 * x)
    chk(abs(f["b"] - 3.7) < 1e-12 and abs(f["R2"] - 1.0) < 1e-12,
        "1 OLS sin intercepto exacto sobre y = 3.7x",
        "b = %.12f  R2 = %.12f" % (f["b"], f["R2"]))

    # --- 2. Recupera una lambda conocida y aplana el offset -----------------
    lam_v = 0.35
    r = _rejilla_sintetica(6000, lam_v, ruido=0.05)
    dP, X, _ = incrementos(r, "Q")
    f2 = ols_sin_intercepto(X, dP)
    pref = offset(r, f2["b"])
    sd_P = float(np.std(np.diff(r["precio"])))
    sd_R = float(np.std(np.diff(pref)))
    chk(abs(f2["b"] - lam_v) / lam_v < 0.02 and sd_R < 0.2 * sd_P,
        "2 recupera lambda conocida y aplana el offset",
        "lambda %.5f (verdad %.2f), sd(dP_ref)/sd(dP) = %.4f"
        % (f2["b"], lam_v, sd_R / sd_P))

    # --- 3. Sin acoplamiento, lambda ~ 0 y el offset NO se aplana -----------
    r3 = _rejilla_sintetica(6000, 0.0, ruido=0.30, semilla=11)
    dP3, X3, _ = incrementos(r3, "Q")
    f3 = ols_sin_intercepto(X3, dP3)
    pref3 = offset(r3, f3["b"])
    vr3 = float(np.var(np.diff(pref3)) / np.var(np.diff(r3["precio"])))
    chk(abs(f3["b"]) < 0.02 and abs(vr3 - 1.0) < 0.02,
        "3 sin acoplamiento: lambda ~ 0 y la razon de varianzas ~ 1",
        "lambda %.5f, Var(dP_ref)/Var(dP) = %.5f" % (f3["b"], vr3))

    # --- 4. La identidad Var(dP_ref)/Var(dP) = 1 - R2 ----------------------
    pref2 = offset(r, f2["b"])
    vr = float(np.var(np.diff(pref2)) / np.var(np.diff(r["precio"])))
    chk(abs(vr - (1.0 - f2["R2"])) < 5e-3,
        "4 la razon de varianzas es 1 - R2", "%.6f contra %.6f" % (vr, 1 - f2["R2"]))

    # --- 5. El offset RECONSTRUYE el precio exactamente ---------------------
    cum = np.concatenate([[0.0], np.cumsum(r["q_neto"][1:])])
    chk(np.max(np.abs(pref2 + f2["b"] * cum - r["precio"])) < 1e-9,
        "5 P_ref + lambda*CumQ reconstruye P exactamente",
        "error maximo %.2e" % np.max(np.abs(pref2 + f2["b"] * cum - r["precio"])))

    # --- 6. Una DERIVA PURA no se cuela en lambda --------------------------
    r6 = _rejilla_sintetica(6000, 0.0, ruido=0.05, semilla=23, deriva=0.5)
    dP6, X6, _ = incrementos(r6, "Q")
    f6 = ols_sin_intercepto(X6, dP6)
    chk(abs(f6["b"]) < 0.02,
        "6 una deriva pura no se cuela en lambda (sin intercepto no la absorbe)",
        "lambda %.5f con deriva de 0.5 USD/casilla, intercepto diag. %.4f"
        % (f6["b"], f6["a_diag"]))

    # --- 7. El nulo por rotacion destruye el acoplamiento -------------------
    su = suelo_offset(r, "Q", n=40)
    chk(su["R2_q95"] < 0.1 * f2["R2"] and abs(su["lam_med"]) < 0.1 * lam_v,
        "7 la rotacion destruye el acoplamiento",
        "R2 %.4f -> q95 nulo %.4f ; lambda %.4f -> nulo %.4f"
        % (f2["R2"], su["R2_q95"], f2["b"], su["lam_med"]))

    # --- 8. Las rachas cortan en los huecos --------------------------------
    r8 = _rejilla_sintetica(300, 0.3, 0.05, semilla=7)
    r8["valida"][100:120] = False
    g = rachas(r8["valida"])
    chk(len(g) == 2 and not any(np.any((x >= 100) & (x < 120)) for x in g),
        "8 las rachas cortan en el hueco y no lo cruzan",
        "%d rachas de tamanos %s" % (len(g), [x.size for x in g]))
    # y el offset se ancla en CADA racha, no arrastra CumQ por el hueco
    p8 = offset(r8, 0.3)
    chk(abs(p8[g[1][0]] - r8["precio"][g[1][0]]) < 1e-9,
        "8b cada racha ancla su CumQ en cero",
        "P_ref al inicio de la 2a racha == P")

    # --- 9. El recorrido del offset entre dias, con verdad conocida --------
    r9 = _rejilla_sintetica(8640 * 3, 0.35, ruido=0.02, semilla=13)
    r9["t"] = 1755000000.0 + np.arange(r9["precio"].size) * 10.0
    dP9, X9, _ = incrementos(r9, "Q")
    f9 = ols_sin_intercepto(X9, dP9)
    fil = tabla_por_dia(r9, rachas(r9["valida"])[0], f9["b"], "Q", "ny")
    if len(fil) >= 2:
        Pm = np.array([x["P_med"] for x in fil]); pm = np.array([x["Pref_med"] for x in fil])
        rp = 1e4 * (np.max(Pm) - np.min(Pm)) / np.median(Pm)
        rr = 1e4 * (np.max(pm) - np.min(pm)) / np.median(Pm)
        chk(rr < 0.1 * rp, "9 con acoplamiento real el offset recorre mucho menos"
            " que el precio", "precio %.1f pb contra offset %.1f pb (%d dias NY)"
            % (rp, rr, len(fil)))
    else:
        chk(False, "9 la tabla por dia produce al menos dos dias",
            "%d dias" % len(fil))

    # --- 10. Un hueco NO fabrica recorrido de offset -----------------------
    # Este es el control del fallo propio que la corrida de verificacion
    # destapo: con `P_ref` anclada por racha, un dia que contiene una frontera
    # de racha marcaba 460 pb de recorrido contra 3 pb de los demas. Trabajando
    # DENTRO de cada racha eso no puede pasar, y aqui se comprueba.
    r10 = _rejilla_sintetica(8640 * 4, 0.35, ruido=0.02, semilla=29)
    r10["t"] = 1755000000.0 + np.arange(r10["precio"].size) * 10.0
    r10["valida"][8640 * 2 - 300: 8640 * 2 + 300] = False   # hueco a mitad
    dP10, X10, _ = incrementos(r10, "Q")
    f10 = ols_sin_intercepto(X10, dP10)
    peor = 0.0
    n_r = 0
    for g10 in rachas(r10["valida"]):
        fl = tabla_por_dia(r10, g10, f10["b"], "Q", "ny")
        if fl:
            n_r += 1
            peor = max(peor, max(x["Pref_rango_pb"] for x in fl))
    chk(n_r == 2 and peor < 20.0,
        "10 un hueco no fabrica recorrido de offset (el fallo del 2026-09-05)",
        "%d rachas, peor recorrido de Pref dentro de un dia %.2f pb" % (n_r, peor))

    log("")
    log("  %d / %d" % (n_ok, n_tot))
    return 0 if n_ok == n_tot else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", choices=("offset",))
    ap.add_argument("--fuente", default="estacional", choices=("estacional", "v33"))
    ap.add_argument("--zona", default="ny", choices=("ny", "utc"))
    F._anadir_rutas(ap)
    ap.add_argument("--sorteos", type=int, default=N_SORTEOS)
    ap.add_argument("--min-dias", dest="min_dias", type=float, default=1.5,
                    help="racha continua minima para comparar dias entre si")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    F.configurar_rutas(a)
    if a.etapa == "offset":
        return etapa_offset(a)
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
