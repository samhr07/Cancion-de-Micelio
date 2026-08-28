# -*- coding: utf-8 -*-
"""
micro_v42.py -- v4.1 §6.2: reajuste de M0/M1/M2 sobre el MICRO-PRECIO.

    python micro_v42.py --autotest   controles con verdad conocida
    python micro_v42.py --etapa=ajuste

Habilitado por el §7.3 de la v4.2 (2026-08-28): el micro-precio de Stoikov da
**x7.2** en observaciones informativas frente al punto medio (238.4 contra 33.3
cambios por cada mil). El §6.2 de la v4.1 condicionaba este reajuste
exactamente a eso.

⚠ **EL CONJUNTO DE PRUEBA DE LA v3.2 NO SE ABRE.** Ya se abrio dos veces. Todo
lo de aqui se ajusta en ENTRENAMIENTO y se evalua en VALIDACION, y por eso las
cifras **no son comparables** con las que la v3.2 publico sobre prueba. La
columna del punto medio se recalcula en validacion **con el mismo codigo** para
que la comparacion entre observables sea interna y limpia.

⚠ **LO QUE ESTE REAJUSTE PUEDE Y NO PUEDE CAMBIAR.** El paso 3 de la v3.2 fallo
porque `q90(|mu|) = 11.30` contra `1.5*c(u) = 39.05` USD/BTC. Un observable con
mas resolucion puede subir `mu` -- corrige el sesgo de atenuacion por
discretizacion --, pero el peaje no se mueve. Y ahora el peaje esta MEDIDO, no
asumido: `c(u)` maker+maker = 4.00 pb LEIDA de la cuenta, y el lastre completo
`L(H)` es 6.2-7.6 pb (v4.2 §2). Se reporta contra los dos.
"""

from __future__ import annotations

import argparse
import os

import numpy as np

import coste as CO
import horizonte as H
import migracion_v32 as M
from experimento_v32 import mu_horizonte

log, titulo = H.log, H.titulo

RUTA = "telemetria/muestra_v32.npz"
EMBARGO = 2233          # el medido por la v3.2 (H*_ticks), no el piso heredado
DELTA = 0.0             # el §5.1 de la v3.2: delta = 0 gana en validacion


def micro_precio(bid, ask, qb, qa):
    s = qa + qb
    return np.where(s > 0, (ask * qb + bid * qa) / np.maximum(s, 1e-30),
                    0.5 * (bid + ask))


def observables(d) -> dict:
    """Los tres, con la MISMA alineacion que uso la v3.2.

    Verificado contra el `y_mid` guardado: `y = [0, diff(precio)]` reproduce el
    vector del cache con error exactamente 0.
    """
    bid, ask = d["bid"].astype(float), d["ask"].astype(float)
    qb, qa = d["qb"].astype(float), d["qa"].astype(float)
    mid = 0.5 * (bid + ask)
    mic = micro_precio(bid, ask, qb, qa)
    pr = d["precio"].astype(float)
    return {"punto medio": np.r_[0.0, np.diff(mid)],
            "MICRO-PRECIO": np.r_[0.0, np.diff(mic)],
            "precio transac.": np.r_[0.0, np.diff(pr)]}


def etapa_ajuste(args) -> int:
    titulo("v4.1 §6.2 -- M0/M1/M2 REAJUSTADOS SOBRE EL MICRO-PRECIO")
    if not os.path.exists(RUTA):
        log("  falta %s" % RUTA)
        return 2
    d = np.load(RUTA)
    eps, v, t = d["eps"].astype(float), d["v"].astype(float), d["t"].astype(float)
    n = eps.size
    K = 2 * EMBARGO
    part = M.particionar(n, EMBARGO)
    # [!] `particionar` devuelve INDICES, no mascaras. Mezclarlos con la mascara
    # de calentamiento reventaba por difusion de formas, y antes de reventar la
    # linea de log imprimia `ent.sum()` -- la SUMA DE LOS INDICES -- como si
    # fuera un recuento: 191 389 586 086 "observaciones". Un numero absurdo que
    # el formato no delataba.
    ent = np.zeros(n, bool); ent[part["entrenamiento"]] = True
    val = np.zeros(n, bool); val[part["validacion"]] = True
    log("")
    log("  n = %d   K = %d   embargo = %d   delta = %.1f" % (n, K, EMBARGO, DELTA))
    log("  entrenamiento %d   validacion %d   PRUEBA: NO SE ABRE" % (ent.sum(), val.sum()))
    cal = M.mascara_calentamiento(t, K)
    pmed = float(np.median(d["precio"]))
    c_u_pb = CO.c_u("maker_maker", 0.0)["comision_pb"]
    c_u = c_u_pb * pmed / 1e4                       # USD/BTC
    L_pb = 6.24                                     # v4.2 §2, H = 15 min
    L_u = L_pb * pmed / 1e4
    log("  c(u) maker+maker LEIDA = %.4f pb = %.4f USD/BTC" % (c_u_pb, c_u))
    log("  lastre completo L(15 min) = %.2f pb = %.4f USD/BTC  (v4.2 §2)" % (L_pb, L_u))

    obs = observables(d)
    resumen = {}
    for nom, y in obs.items():
        log("")
        log("--- %s ---   ceros: %.2f %%" % (nom, 100 * np.mean(y == 0)))
        m0 = M.ajustar_M0(y[ent])
        m1 = M.ajustar(y[ent], eps[ent], v[ent], K, DELTA, beta_libre=False)
        m2 = M.ajustar(y[ent], eps[ent], v[ent], K, DELTA, beta_libre=True,
                       con_suelo=True)
        vv = val & cal
        base = None
        log("  %-6s %14s %14s | %10s %10s %10s"
            % ("modelo", "LL/N valid.", "dLL/N vs M0", "G0", "beta", "f_inf"))
        for nm, mod in (("M0", m0), ("M1", m1), ("M2", m2)):
            ll = float(np.mean(M.ll_por_obs(mod, y[vv], eps[vv], v[vv])))
            if base is None:
                base = ll
            log("  %-6s %14.8f %14.3e | %10.6f %10.4f %10.4f"
                % (nm, ll, ll - base, mod["G0"], mod["beta"], mod.get("f_inf", 0.0)))
            resumen.setdefault(nom, {})[nm] = {"ll": ll, "d": ll - base,
                                               "mod": mod}
        # --- paso 3 sobre VALIDACION
        mu = mu_horizonte(m2, eps[vv], v[vv], EMBARGO)[K:]
        am = np.abs(mu)
        if am.size:
            q90 = float(np.percentile(am, 90))
            log("  |mu| de M2 a H = %d ticks:  q50 %.4f  q90 %.4f  q99 %.4f  max %.4f USD/BTC"
                % (EMBARGO, np.percentile(am, 50), q90, np.percentile(am, 99), am.max()))
            log("     contra 1.5*c(u) = %.4f  ->  %s" % (1.5 * c_u,
                                                         "PASA" if q90 >= 1.5 * c_u else "falla"))
            log("     contra 1.5*L    = %.4f  ->  %s" % (1.5 * L_u,
                                                         "PASA" if q90 >= 1.5 * L_u else "falla"))
            resumen[nom]["q90"] = q90
        # --- residuo blanco
        pred = M.predecir(m2, eps[vv], v[vv])
        res = (y[vv] - pred)[K:]
        if res.size > 200:
            r1 = float(np.corrcoef(res[:-1], res[1:])[0, 1])
            log("  residuo de M2: rho1 = %+.4f   %s"
                % (r1, "blanco" if abs(r1) < 0.05 else "*** ESTRUCTURADO ***"))

    log("")
    titulo("COMPARACION ENTRE OBSERVABLES (todo en VALIDACION)")
    log("  %-18s %10s %14s %14s %12s"
        % ("observable", "% ceros", "dLL/N M2-M0", "q90(|mu|)", "q90/(1.5c)"))
    for nom, y in obs.items():
        r = resumen[nom]
        q = r.get("q90", float("nan"))
        log("  %-18s %9.2f %% %14.3e %14.4f %12.3f"
            % (nom, 100 * np.mean(y == 0), r["M2"]["d"], q, q / (1.5 * c_u)))
    log("")
    log("  [!] Estas cifras son de VALIDACION. Las de la v3.2 eran de PRUEBA y no")
    log("      son comparables; la columna del punto medio se recalcula aqui con")
    log("      el mismo codigo para que la comparacion entre observables sea interna.")
    return 0


def _autotest() -> int:
    titulo("CONTROLES DE micro_v42.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-52s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    # la alineacion reproduce el y_mid del cache con error 0
    if os.path.exists(RUTA):
        d = np.load(RUTA)
        o = observables(d)
        err = float(np.max(np.abs(o["punto medio"] - d["y_mid"].astype(float))))
        chk(err == 0.0, "la alineacion reproduce `y_mid` del cache EXACTO",
            "error maximo %.3g" % err)
        chk(np.mean(o["MICRO-PRECIO"] == 0) < np.mean(o["punto medio"] == 0),
            "el micro-precio tiene menos ceros que el punto medio",
            "%.2f %% contra %.2f %%" % (100 * np.mean(o["MICRO-PRECIO"] == 0),
                                        100 * np.mean(o["punto medio"] == 0)))
    # micro con colas simetricas = punto medio
    b = np.array([100.0]); a = np.array([100.1])
    chk(abs(micro_precio(b, a, np.array([5.0]), np.array([5.0]))[0] - 100.05) < 1e-12,
        "colas simetricas -> micro = punto medio")
    # mu es CERO identicamente con f_inf = 1 (M1), que es el contenido de MkII
    mod = {"G0": 1.0, "K": 10, "tau0": 5.0, "beta": 0.0, "f_inf": 1.0,
           "delta": 0.0, "v_mediana": 1.0, "omega_G": 0.0, "phi": 0.0}
    mu = mu_horizonte(mod, np.ones(50), np.ones(50), 5)
    chk(float(np.max(np.abs(mu))) < 1e-12,
        "con impacto permanente (f_inf=1) mu es CERO identicamente",
        "max |mu| = %.3g" % float(np.max(np.abs(mu))))
    # y NO es cero con nucleo transitorio
    mod2 = dict(mod); mod2["f_inf"] = 0.0; mod2["beta"] = 0.5
    mu2 = mu_horizonte(mod2, np.ones(50), np.ones(50), 5)
    chk(float(np.max(np.abs(mu2))) > 1e-6, "y NO es cero con nucleo transitorio",
        "max |mu| = %.4f" % float(np.max(np.abs(mu2))))
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
    if a.etapa == "ajuste":
        return etapa_ajuste(a)
    log("usa --autotest o --etapa=ajuste")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
