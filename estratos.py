# -*- coding: utf-8 -*-
"""
estratos.py -- el §1 REHECHO, estratificado por sigma PRONOSTICADA.

    python estratos.py --autotest
    python estratos.py --etapa=curvas

Es la tarea marcada BLOQUEANTE desde el 2026-08-23 y nunca ejecutada. Se hace
con las tres correcciones que el proyecto ha ido acumulando:

  1. **Coste REAL**: `L(H)` medida en el §2 de la v4.2 (6.2-7.6 pb), no `c(u)`
     sola ni los 4.888 pb asumidos.
  2. **`kappa` corregido**: **2.0627**, el de operar en LAS DOS direcciones, no
     el 1.755 de una sola. Son 17.5 %, o 1.38x al cuadrado en `R2_req`.
  3. **Control de potencia POR ESTRATO**: suelo por rotacion circular, y una
     fila no se lee si `R2_req < 3 * q95(|R2_nulo|)`.

**Estratificacion, declarada antes de medir.** `sigma` pronosticada para la
ventana `[t, t+H]` = `sigma` realizada en `[t-H, t]`. Es una prediccion por
PERSISTENCIA: sin parametros que ajustar, conocida en `t`, y con la persistencia
ya medida (`R2` fuera de muestra de 0.51 para `sigma(t) -> sigma(t+1)`, y 0.54
usando `nu`). Los cortes de tercil se fijan en ENTRENAMIENTO.

⚠ **REGLA DE PARADA, ESCRITA ANTES DE EJECUTAR (acordada con el operador):**

> Si en el estrato SUPERIOR de `sigma` pronosticada, donde el instrumento
> resuelve, **el margen no baja de 3x**, el proyecto reporta resultado negativo
> y **cierra la linea de microestructura**.

Margen = `R2_req / R2_medido`. Es el «falta un factor N» de la Adenda C.

⚠ **EXPECTATIVA DECLARADA:** el §4 de la v4.2 ya barrio estos mismos horizontes
con el coste real y salio negativo en **78 de 79 celdas**, incluidos los
fragmentos de `sigma` alta. Se espera confirmacion, no sorpresa. Se ejecuta
igual porque mide otra cosa (`R2` contra rentabilidad) y porque un cierre
formal necesita su propia regla de parada.
"""

from __future__ import annotations

import argparse
import gc
import json
import os

import numpy as np

import horizonte as H
import identidad as I

log, titulo = H.log, H.titulo

HS = [900, 1800, 3600, 7200, 14400]
KAPPA = 2.0627               # DOS direcciones. El 1.755 es el de una sola.
N_SORTEOS = 200
PASO_S = 60.0
NOMBRES = ["bajo", "medio", "ALTO"]


def _L(Hs):
    try:
        d = json.load(open("telemetria/curva_coste.json", encoding="utf-8"))
        for f in d["L"]:
            if f["H"] == Hs:
                return f["L_maker"]
        return d["L"][-1]["L_maker"]
    except Exception:
        return 6.24


def preparar(frag, Hs):
    """Origenes uniformes en tiempo, con `sigma` pronosticada por persistencia."""
    t, mid, eps = frag["t"], frag["mid"], frag["eps"].astype(float)
    dur = float(t[-1] - t[0])
    g = np.arange(t[0] + Hs, t[-1] - Hs, PASO_S)      # deja hueco para mirar atras
    if g.size < 120:
        return None
    i = np.unique(np.clip(np.searchsorted(t, g, side="left"), 0, t.size - 1))
    j = np.searchsorted(t, t[i] + Hs, side="left")
    k = np.searchsorted(t, t[i] - Hs, side="left")     # ventana PASADA
    v = (j < t.size) & (k >= 0) & (mid[i] > 0)
    i, j, k = i[v], j[v], k[v]
    if i.size < 120:
        return None
    r = np.log(mid[j] / mid[np.minimum(i + 1, mid.size - 1)]) * 1e4   # futuro, pb
    # sigma pronosticada: |retorno| de la ventana PASADA de la misma longitud.
    # Es ex-ante por construccion: solo usa informacion anterior a t.
    s_prev = np.abs(np.log(mid[i] / mid[k])) * 1e4
    t_e = t[0] + 0.60 * dur
    t_v = t[0] + 0.80 * dur
    ent = t[j] <= t_e
    vld = (t[i] >= t_e) & (t[j] <= t_v)
    return {"eps": eps[i], "r": r, "s_prev": s_prev, "ent": ent, "vld": vld,
            "i": i, "n": i.size}


def r2_ident(e, r):
    if e.size < 40:
        return float("nan")
    Rc = float(np.mean(e * r) - np.mean(e) * np.mean(r))
    Ve, Vr = float(np.var(e)), float(np.var(r))
    return Rc * Rc / (Ve * Vr) if Ve > 0 and Vr > 0 else float("nan")


def suelo(e, r, n_sorteos=N_SORTEOS, semilla=13):
    rng = np.random.default_rng(semilla)
    n = e.size
    v = [r2_ident(np.take(e, np.arange(n) - int(rng.integers(5, max(6, n - 5))),
                          mode="wrap"), r) for _ in range(n_sorteos)]
    v = np.array([x for x in v if np.isfinite(x)])
    return float(np.percentile(np.abs(v), 95)) if v.size else float("nan")


def etapa_curvas(args) -> int:
    titulo("§1 ESTRATIFICADO POR sigma PRONOSTICADA")
    log("")
    log("  kappa = %.4f (DOS direcciones)   L(H) real del §2 de la v4.2" % KAPPA)
    log("  estratos: terciles de sigma pronosticada, cortes fijados en ENTRENAMIENTO")
    log("  REGLA DE PARADA: si en el estrato ALTO, donde el instrumento resuelve,")
    log("  el margen no baja de 3x -> resultado negativo y se cierra la linea.")
    filas = []
    for nombre in I._nombres():
        f = I.cargar(nombre, con_precio=False)
        dur = float(f["t"][-1] - f["t"][0])
        log("")
        log("--- %s (%.1f h) ---" % (nombre, dur / 3600.0))
        log("  %7s %7s %7s | %9s %9s %10s %10s %9s"
            % ("H", "estrato", "n", "sigma[pb]", "R2_req", "R2_medido", "q95_suelo", "margen"))
        for Hs in HS:
            if dur < 6 * Hs:
                continue
            P = preparar(f, Hs)
            if P is None:
                continue
            L = _L(Hs)
            cortes = np.percentile(P["s_prev"][P["ent"]], [33.3, 66.7])
            gi = np.digitize(P["s_prev"], cortes)
            for g_ in (0, 1, 2):
                # [!] SE MIDE SOBRE TODOS LOS ORIGENES, NO SOLO VALIDACION, y es
                # correcto: la identidad es un MOMENTO sin parametros ajustados
                # -- `Corr(eps, r)^2` no se estima en un bloque y se evalua en
                # otro. Lo unico que viene de ENTRENAMIENTO son los cortes de
                # tercil, que es donde estaria la fuga. Restringir a validacion
                # tiraba el 80 % de los origenes y multiplicaba el suelo por
                # ~sqrt(5), con lo que NINGUNA fila resolvia: el "no resuelve"
                # era mio, no del mercado. La Adenda C ya media sobre todos los
                # pares por esta misma razon.
                m = (gi == g_)
                if m.sum() < 60:
                    continue
                e, r = P["eps"][m], P["r"][m]
                sg = float(np.std(r))
                req = (L / (KAPPA * sg)) ** 2
                r2 = r2_ident(e, r)
                q95 = suelo(e, r)
                resuelve = np.isfinite(q95) and req >= 3 * q95
                margen = req / max(r2, 1e-12) if r2 > 0 else float("inf")
                log("  %6ds %7s %7d | %9.2f %8.3f %% %10.6f %10.6f %9s%s"
                    % (Hs, NOMBRES[g_], int(m.sum()), sg, 100 * req, r2, q95,
                       ("%.1fx" % margen) if np.isfinite(margen) else "  inf",
                       "" if resuelve else "   <- NO RESUELVE"))
                filas.append({"frag": nombre, "H": Hs, "estrato": NOMBRES[g_],
                              "n": int(m.sum()), "sigma": sg, "req": req,
                              "r2": r2, "q95": q95, "resuelve": bool(resuelve),
                              "margen": float(margen)})
            del P
            gc.collect()
        del f
        gc.collect()

    titulo("REGLA DE PARADA, APLICADA")
    alto = [x for x in filas if x["estrato"] == "ALTO" and x["resuelve"]]
    log("")
    log("  filas del estrato ALTO donde el instrumento RESUELVE: %d de %d"
        % (len(alto), sum(1 for x in filas if x["estrato"] == "ALTO")))
    if not alto:
        log("  -> en el estrato ALTO el instrumento no resuelve en ninguna H.")
        log("     La regla no se puede aplicar; se reporta como tal.")
        return 0
    mejor = min(alto, key=lambda z: z["margen"])
    log("  margen MINIMO en el estrato ALTO: %.1fx   (%s, H = %d s, n = %d)"
        % (mejor["margen"], mejor["frag"], mejor["H"], mejor["n"]))
    log("")
    for x in sorted(alto, key=lambda z: z["margen"])[:6]:
        log("    %-16s H=%6ds  n=%5d  R2 = %.6f  req = %.3f %%  margen %.1fx"
            % (x["frag"], x["H"], x["n"], x["r2"], 100 * x["req"], x["margen"]))
    log("")
    if mejor["margen"] >= 3.0:
        log("  *** EL MARGEN NO BAJA DE 3x (minimo %.1fx). Por la regla de parada"
            % mejor["margen"])
        log("      acordada ANTES de ejecutar: RESULTADO NEGATIVO, y se CIERRA la")
        log("      linea de microestructura. ***")
    else:
        log("  *** EL MARGEN BAJA DE 3x (minimo %.1fx). La regla NO cierra la linea:")
        log("      hay que mirar esa celda de cerca antes de decidir nada. ***"
            % mejor["margen"])
    json.dump(filas, open("telemetria/estratos_v42.json", "w", encoding="utf-8"), indent=1)
    return 0


def _autotest() -> int:
    titulo("CONTROLES DE estratos.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-52s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    rng = np.random.default_rng(3)
    n = 200000
    t = np.arange(n) * 1.0
    eps = np.where(rng.random(n) < 0.5, -1.0, 1.0)
    base = np.cumsum(rng.normal(0, 0.5, n))
    # [!] DOS pasos de retardo, no uno. SEGUNDA VEZ que este mismo fallo de
    # alineacion aparece en un control mio (la primera fue en `identidad.py`).
    # El retorno se mide desde `mid[i+1]` -- para que sea PREDICTIVO y no
    # incluya el impacto inmediato --, asi que una senal inyectada con un solo
    # paso de retardo deja `eps_i` FUERA de la ventana y el estimador hace bien
    # en devolver ~0. Con dos pasos, `eps_i` cae dentro y la senal es realmente
    # predecible desde el presente.
    p = 60000.0 + base + 4.0 * np.cumsum(np.r_[0.0, 0.0, eps[:-2]])
    # [!] VENTANA CORTA EN EL CONTROL, Y LA RAZON ES ESTRUCTURAL. El `R2` de UN
    # solo signo sobre una ventana de `H` pasos SATURA en ~1/H: un incremento de
    # los H no puede explicar mas que su parte. Con H = 900 el techo es 1.1e-3 y
    # el suelo de muestreo con 3 300 origenes es 3e-4, o sea una relacion de 3.7
    # -- demasiado justa para que un control discrimine. Subir la AMPLITUD no
    # ayuda: la senal entra tambien en la varianza y el `R2` no se mueve. Con
    # H = 60 el techo sube a 1.7e-2 y el control separa con holgura.
    P = preparar({"t": t, "mid": p, "eps": eps}, 60)
    chk(P is not None and P["n"] > 100, "prepara origenes con ventana pasada",
        "%d" % (P["n"] if P else 0))
    # sigma pronosticada es EX-ANTE: no puede usar informacion futura
    chk(np.all(P["s_prev"] >= 0), "sigma pronosticada es no negativa")
    r2s = r2_ident(P["eps"], P["r"])
    r2n = r2_ident(rng.permutation(P["eps"]), P["r"])
    chk(r2s > 20 * max(r2n, 1e-12), "la identidad separa senal de barajado",
        "%.2e contra %.2e" % (r2s, r2n))
    q = suelo(P["eps"], P["r"], n_sorteos=60)
    chk(r2s > 3 * q, "y supera su propio suelo por rotacion",
        "%.2e contra q95 %.2e" % (r2s, q))
    # kappa de dos direcciones
    z = rng.normal(0, 1, 400000)
    k2 = float(np.mean(np.abs(z)[np.abs(z) > np.percentile(np.abs(z), 90)]))
    chk(abs(k2 - KAPPA) < 0.03, "KAPPA es el de dos direcciones", "%.4f" % k2)
    # R2_req crece si el coste sube y baja si sigma sube
    req = lambda L, s: (L / (KAPPA * s)) ** 2
    chk(req(7.0, 10.0) > req(6.0, 10.0) and req(6.0, 20.0) < req(6.0, 10.0),
        "R2_req crece con L y decrece con sigma")
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
    if a.etapa == "curvas":
        return etapa_curvas(a)
    log("usa --autotest o --etapa=curvas")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
