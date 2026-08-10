# -*- coding: utf-8 -*-
"""
experimento_v32.py -- ejecuta `PREREGISTRO_3_2.md` sobre `captura_v33`.

    python experimento_v32.py --etapa=muestra     construye y cachea la muestra
    python experimento_v32.py --etapa=fuga        §4.2 test de fuga (4 variantes)
    python experimento_v32.py --etapa=delta       §5.1 eleccion de delta EN VALIDACION
    python experimento_v32.py --etapa=A           §5.2 D, q05 regenerado, f_inf
    python experimento_v32.py --etapa=osc         §5.3 omega_G contra nulo simulado
    python experimento_v32.py --etapa=B           §6   lado de limite y variantes c_t
    python experimento_v32.py --etapa=decision    §9   la regla, en orden
    python experimento_v32.py --etapa=todo

⚠ ESTE MODULO NO ENMIENDA EL PREREGISTRO. Cualquier decision que no estuviera
escrita antes de tocar `captura_v33` invalida el caracter fuera de muestra del
resultado (§12.bis, regla final). Lo que hay aqui es la ejecucion literal.

ORDEN DE APERTURA DE LOS CONJUNTOS
----------------------------------
entrenamiento -> validacion -> prueba, y **prueba se abre UNA sola vez**, en la
etapa `decision`. Las etapas `fuga`, `delta` y `A` trabajan sobre entrenamiento
y validacion. Si una etapa posterior necesitara volver a prueba, el resultado
deja de ser fuera de muestra y hay que reportarlo asi (§2.1).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

import migracion_v32 as M
import propagador as P
from captura_larga import cargar_larga, ofi_l1, tramos_continuos


DIR_CAPTURA = "telemetria/captura_v33"
CACHE = "telemetria/muestra_v32.npz"

# §11 del preregistro: horizonte del criterio economico. c(u) del esquema
# declarado, con tarifas ASUMIDAS VIP 0 -- criterio "comisiones leidas de la
# cuenta" NO CUMPLIDO, igual que en la v3.1.
ESQUEMA_COSTE = "maker_maker"


_TRANSLIT = {"⚠": "[!]", "§": "Sec.", "→": "->", "≤": "<=",
             "≥": ">=", "≠": "!=", "×": "x", "‑": "-",
             "–": "-", "—": "--", "²": "^2", "√": "sqrt"}


def log(*a):
    """Imprime SIEMPRE en ASCII.

    [!] La consola de Windows es cp1252 y no puede codificar ni los simbolos de
    seccion ni las griegas. Esto ha roto una corrida en cuatro sesiones
    distintas, siempre por confiar en la disciplina de escribir ASCII a mano.
    La conversion vive aqui, en el unico sitio por el que pasa todo, para que
    dejar de ser una regla que recordar.
    """
    s = " ".join(str(x) for x in a)
    for k, v in _TRANSLIT.items():
        s = s.replace(k, v)
    print(s.encode("ascii", "replace").decode("ascii"))
    sys.stdout.flush()


def titulo(s):
    log("")
    log("=" * 74)
    log(s)
    log("=" * 74)


# ===========================================================================
# 1. Muestra: tramo continuo, alineacion del libro, observables
# ===========================================================================

def construir_muestra(directorio: str = DIR_CAPTURA, desplazamiento: int = 0,
                      barajar: bool = False, semilla: int = 0) -> dict:
    """Tramo continuo mas largo, con el libro alineado segun el §4.1.

    `desplazamiento` y `barajar` implementan las variantes del test de fuga del
    §4.2: NO son opciones de analisis, son controles.
    """
    d = cargar_larga(directorio)

    # Ceros del feed: la limitacion 3 del §12 ordena descartarlos en origen.
    tp, tt = d["tr_precio"], d["tr_t"]
    ok = np.isfinite(tp) & (tp > 0) & (d["tr_cant"] > 0)
    tp, tt = tp[ok], tt[ok]
    tq, tm = d["tr_cant"][ok], d["tr_maker"][ok].astype(bool)

    # Tramo continuo MAS LARGO, nunca la suma de trozos (§1.1).
    tramos = tramos_continuos(tt)
    a, b = max(tramos, key=lambda p: p[1] - p[0])
    tp, tt, tq, tm = tp[a:b], tt[a:b], tq[a:b], tm[a:b]

    bb, bB = d["bk_b"], d["bk_B"]
    ba, bA = d["bk_a"], d["bk_A"]
    bt = d["bk_t"]
    okb = np.isfinite(bb) & (bb > 0) & np.isfinite(ba) & (ba > 0)
    bb, bB, ba, bA, bt = bb[okb], bB[okb], ba[okb], bA[okb], bt[okb]
    orden = np.argsort(bt, kind="stable")
    bb, bB, ba, bA, bt = bb[orden], bB[orden], ba[orden], bA[orden], bt[orden]

    mid = 0.5 * (bb + ba)
    ofi = ofi_l1(bb, bB, ba, bA)          # longitud L-1, alineado al SEGUNDO

    # §4.1: ultimo bookTicker con timestamp ESTRICTAMENTE anterior. Nunca
    # interpolar, nunca ffill hacia atras.
    j = np.searchsorted(bt, tt, side="left") - 1

    if barajar:
        rng = np.random.default_rng(semilla)
        j = rng.permutation(j)
    if desplazamiento:
        j = j + desplazamiento

    val = (j >= 1) & (j < len(bt))
    j, tp, tt, tq, tm = j[val], tp[val], tt[val], tq[val], tm[val]

    eps = P.signo_transaccion(tm)

    # y_t = incremento del punto medio entre la transaccion anterior y esta.
    y_mid = mid[j][1:] - mid[j][:-1]
    y_tr = tp[1:] - tp[:-1]
    e_t = ofi[j[1:] - 1]
    V_best = np.where(eps[1:] > 0, bA[j[1:]], bB[j[1:]])
    spread = (ba - bb)[j[1:]]

    return {"y_mid": y_mid, "y_tr": y_tr, "eps": eps[1:], "v": tq[1:],
            "t": tt[1:], "e": e_t, "V_best": V_best, "spread": spread,
            "precio": tp[1:], "j": j[1:],
            "q_mejor_nivel": np.where(eps[1:] > 0, bA[j[1:]], bB[j[1:]]),
            "bid": bb[j[1:]], "ask": ba[j[1:]],
            "qb": bB[j[1:]], "qa": bA[j[1:]],
            "qb_prev": bB[j[1:] - 1], "qa_prev": bA[j[1:] - 1]}


def cachear(m: dict, ruta: str = CACHE):
    os.makedirs(os.path.dirname(ruta), exist_ok=True)
    np.savez_compressed(ruta, **m)


def leer_cache(ruta: str = CACHE) -> dict:
    d = np.load(ruta)
    return {k: d[k] for k in d.files}


# ===========================================================================
# 2. Embargo, K y el piso heredado
# ===========================================================================

def estimar_H_ticks(precio_ent: np.ndarray, c_usd: float) -> dict:
    """`H*_ticks`: el rezago en ticks al que el movimiento tipico IGUALA `c(u)`.

    ⚠ RESOLUCION DE UNA AMBIGUEDAD DEL PREREGISTRO, declarada aqui y no
    escondida. El §2.2 escribe `H*_ticks = (c/sigma_tick)^2` y NO dice cual
    `sigma_tick`. Esa formula es la solucion de `sigma(H) = c` **bajo difusion
    exacta** (`sigma(n) = sigma_tick*sqrt(n)`), o sea presupone que
    `sigma(n)/sqrt(n)` es constante. Medido sobre esta captura NO lo es: la
    curva es en U -- cae de 0.799 en n=2 a 0.327 en n=64 y sube monotona hasta
    0.608 en n=8192, sin plateau dentro del rango medido. Con esa forma,
    `(c/sigma_tick)^2` da un numero distinto por cada `sigma_tick` que se elija
    (1062 ticks con el de n=2, 1836 con el de n=8192) y ninguno cumple
    `sigma(H*) = c`.

    Se resuelve por el PUNTO FIJO, que es lo que la formula quiere decir:
    interpolar `log sigma` contra `log n` y despejar `sigma(H) = c`. Bajo
    difusion los dos caminos coinciden exactamente, asi que esto no cambia la
    regla: la instancia.

    ⚠ Esta eleccion se toma DESPUES de ver la firma, asi que se declara como
    lectura del preregistro y no como enmienda. Se reportan las dos familias de
    numeros para que el lector juzgue.
    """
    rej = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]
    filas = P.firma_en_ticks(precio_ent, rej, solapada=False)
    n = np.array([f["n"] for f in filas], dtype=float)
    sig = np.array([f["sigma"] for f in filas], dtype=float)
    ln, ls = np.log(n), np.log(sig)
    lc = np.log(c_usd)
    if lc <= ls[0]:
        H = float(n[0])
        dentro = False
    elif lc >= ls[-1]:
        # Extrapolacion con la pendiente del tramo final; se marca como tal.
        p = (ls[-1] - ls[-2]) / (ln[-1] - ln[-2])
        H = float(np.exp(ln[-1] + (lc - ls[-1]) / p))
        dentro = False
    else:
        H = float(np.exp(np.interp(lc, ls, ln)))
        dentro = True
    s = np.array([f["sigma_por_raiz_n"] for f in filas])
    return {"H_ticks": H, "dentro_del_rango": dentro, "filas": filas,
            "sigma_min": float(s.min()), "n_sigma_min": int(n[int(np.argmin(s))]),
            "sigma_ultimo": float(s[-1]), "n_ultimo": int(n[-1]),
            "H_ingenuo_min": (c_usd / float(s.min())) ** 2,
            "H_ingenuo_ultimo": (c_usd / float(s[-1])) ** 2}


# ===========================================================================
# 3. Evaluacion con calentamiento descartado
# ===========================================================================

def ll_media(mod: dict, y, eps, v, K: int) -> tuple[float, np.ndarray]:
    """`LL/N` descartando las primeras K filas del tramo evaluado.

    Sin ese recorte, las primeras K observaciones de un conjunto se predicen con
    forzamiento que el estimador no tiene (§2.bis de `migracion_v32`).
    """
    ll = M.ll_por_obs(mod, y, eps, v)
    ll = ll[K:] if K > 0 and len(ll) > K + 100 else ll
    return float(np.mean(ll)), ll


def rho_lags(r: np.ndarray, lags) -> dict:
    r = np.asarray(r, float)
    r = r - r.mean()
    d = float(np.dot(r, r))
    out = {}
    for L in lags:
        if L < len(r):
            out[L] = float(np.dot(r[L:], r[:-L]) / d)
    return out


# ===========================================================================
# 4. mu al horizonte -- el estadistico del criterio economico (§9.1)
# ===========================================================================

def mu_horizonte(mod: dict, eps, v, H: int) -> np.ndarray:
    """`mu_t = E[ suma de los H incrementos siguientes | F_t ]`.

    Bajo el propagador, lo que el modelo sabe del futuro es la relajacion del
    forzamiento YA observado:

        mu_t(H) = G0 * SUM_{j>=0} x_{t-j} * [ G(H+j) - G(j) ]

    ⚠ Con `f_inf = 1` (M1/MkII) el nucleo es constante y `mu = 0` IDENTICAMENTE:
    impacto permanente significa que el movimiento ya ocurrio y no queda nada
    que capturar. Eso no es un defecto del calculo, es el contenido economico de
    la hipotesis, y es la razon de que el paso 3 se evalue sobre `mu` y no sobre
    el ajuste.
    """
    if mod["G0"] == 0.0:
        return np.zeros(len(eps))
    K = mod["K"]
    x = M.forzamiento(eps, v, mod["delta"], mod["v_mediana"])
    tau = np.arange(K + 1, dtype=float)
    G = M.nucleo_G(tau, mod["tau0"], mod["beta"], mod.get("f_inf", 0.0))
    Gh = M.nucleo_G(tau + H, mod["tau0"], mod["beta"], mod.get("f_inf", 0.0))
    w = Gh - G
    return mod["G0"] * M.convolucion_causal(x, w)


# ===========================================================================
# ETAPAS
# ===========================================================================

def etapa_muestra(args) -> dict:
    titulo("ETAPA 1 -- MUESTRA (§4.1) y DIAGNOSTICOS DEL §4.3")
    t0 = time.time()
    m = construir_muestra(args.dir)
    log("construida en %.1f s" % (time.time() - t0))

    n = len(m["y_mid"])
    dur = float(m["t"][-1] - m["t"][0])
    nu = n / dur
    log("observaciones      : %d" % n)
    log("duracion           : %.2f h   nu = %.2f tx/s" % (dur / 3600.0, nu))
    log("hueco maximo       : %.1f s" % float(np.max(np.diff(m["t"]))))
    log("precio             : %.1f -- %.1f USD" % (m["precio"].min(), m["precio"].max()))

    sp = m["spread"]
    log("spread mediano     : %.4f USD  (p90 %.4f)" % (np.median(sp), np.percentile(sp, 90)))
    log("y_mid nulos        : %.2f %%" % (100.0 * np.mean(m["y_mid"] == 0)))
    log("y_tr  nulos        : %.2f %%" % (100.0 * np.mean(m["y_tr"] == 0)))

    # --- §4.3: OFI-L1 APROXIMADO. Umbral declarado del 20 %.
    excede = float(np.mean(m["v"] > m["q_mejor_nivel"]))
    log("")
    log("--- §4.3  DIAGNOSTICO OBLIGATORIO -- esto es OFI-L1 APROXIMADO ---")
    log("q > cantidad del mejor nivel del ultimo snapshot : %.2f %%" % (100 * excede))
    if excede > 0.20:
        log("  [COTA INFERIOR] supera el 20 % declarado: M1/M1'/M2+L se rotulan asi")
    else:
        log("  [OK] por debajo del 20 % declarado: M1 es medicion, no cota")

    # Reconciliacion de las variaciones de cantidad con las transacciones vistas.
    dqb = m["qb"] - m["qb_prev"]
    dqa = m["qa"] - m["qa_prev"]
    consumo = np.where(m["eps"] > 0, -dqa, -dqb)     # lo que la transaccion deberia comer
    residuo = consumo - m["v"]
    denom = float(np.sum(np.abs(consumo)) + np.sum(np.abs(m["v"])))
    log("residuo de reconciliacion dV vs transacciones    : %.1f %%"
        % (100.0 * np.sum(np.abs(residuo)) / max(denom, 1e-12)))
    log("  (el residuo ES actividad de limite no observada; no es un error)")

    # --- prediccion falsable del signo: G(0) > 0 (v3.1 §2.2)
    R1 = float(np.mean(m["y_mid"] * m["eps"]))
    log("")
    log("E[y_mid * eps] (impacto instantaneo, USD/BTC)    : %+.6f" % R1)
    if R1 <= 0:
        log("  *** G(0) <= 0 SOBRE PUNTO MEDIO: convencion de signo a revisar ***")
    else:
        log("  [OK] positivo: la convencion m=True -> eps=-1 se sostiene")

    cachear(m)
    log("")
    log("cache: %s" % CACHE)
    return m


def preparar(args):
    """Muestra + embargo + particion. Comun a todas las etapas de modelo."""
    m = leer_cache() if os.path.exists(CACHE) else construir_muestra(args.dir)
    n = len(m["y_mid"])
    dur = float(m["t"][-1] - m["t"][0])
    nu = n / dur
    precio = float(np.median(m["precio"]))

    roll = P.spread_efectivo_roll(np.diff(m["precio"]))
    c_u = P.coste_ida_y_vuelta(precio, ESQUEMA_COSTE,
                               s_eff=roll["s_eff"] if np.isfinite(roll["s_eff"]) else 0.0)

    part_prov = M.particionar(n, M.EMBARGO_PISO)
    ent = part_prov["entrenamiento"]
    h = estimar_H_ticks(m["precio"][ent], c_u)
    embargo = int(max(round(h["H_ticks"]), M.EMBARGO_PISO))
    K = 2 * embargo
    part = M.particionar(n, embargo)

    # CONTROL de la firma: incrementos barajados. Sin el no se distingue
    # super-difusion real de sesgo del estimador -- la leccion de la Adenda A.
    rng = np.random.default_rng(101)
    inc = np.diff(m["precio"][ent])
    p_nulo = np.concatenate(([m["precio"][ent][0]], m["precio"][ent][0]
                             + np.cumsum(rng.permutation(inc))))
    firma_nulo = P.firma_en_ticks(p_nulo, [f["n"] for f in h["filas"]], solapada=False)

    return m, dict(n=n, dur=dur, nu=nu, precio=precio, roll=roll, c_u=c_u,
                   H=h, embargo=embargo, K=K, part=part, firma_nulo=firma_nulo)


def cabecera_embargo(cfg):
    titulo("ETAPA 2 -- EMBARGO, K Y EL PISO HEREDADO (§2.2)")
    h = cfg["H"]
    log("c(u) esquema %-12s : %.4f USD/BTC   [TARIFAS ASUMIDAS VIP 0]"
        % (ESQUEMA_COSTE, cfg["c_u"]))
    log("s_eff de Roll (limpio)   : %.4f USD/BTC  (rho1 = %+.4f)"
        % (cfg["roll"]["s_eff"], cfg["roll"]["rho1"]))
    log("")
    log("firma de volatilidad EN TICKS, no solapada (entrenamiento):")
    log("   n      sigma      sigma/sqrt(n)")
    for f in h["filas"]:
        log("  %6d  %9.4f  %12.4f" % (f["n"], f["sigma"], f["sigma_por_raiz_n"]))
    log("")
    log("⚠ NO hay plateau: la curva es en U (min %.4f en n=%d, %.4f en n=%d)."
        % (h["sigma_min"], h["n_sigma_min"], h["sigma_ultimo"], h["n_ultimo"]))
    log("  (c/sigma)^2 ingenuo daria %.0f o %.0f ticks segun que sigma se elija;"
        % (h["H_ingenuo_min"], h["H_ingenuo_ultimo"]))
    log("  ninguno cumple sigma(H*) = c. Se resuelve por PUNTO FIJO -- ver docstring.")
    log("H*_ticks tal que sigma(H) = c(u) : %.0f ticks = %.1f s a nu = %.2f  %s"
        % (h["H_ticks"], h["H_ticks"] / cfg["nu"], cfg["nu"],
           "" if h["dentro_del_rango"] else "  *** EXTRAPOLADO ***"))
    log("piso declarado           : %d ticks" % M.EMBARGO_PISO)
    log("EMBARGO                  : %d ticks   K = 2*embargo = %d"
        % (cfg["embargo"], cfg["K"]))
    if cfg["embargo"] == M.EMBARGO_PISO:
        log("  -> el PISO ata. El §2.2 obliga a justificarlo o retirarlo: ver reporte")
    else:
        log("  -> H*_ticks MEDIDO ata; el piso NO esta haciendo el trabajo.")
        log("     Eso RESUELVE la pregunta abierta del §2.2 y el piso se puede retirar.")
    log("  compuerta con este embargo: %d ticks; disponibles %d  -> %s"
        % (M.ticks_de_compuerta(cfg["embargo"]), cfg["n"],
           "PASA" if cfg["n"] >= M.ticks_de_compuerta(cfg["embargo"]) else "*** NO PASA ***"))
    if "firma_nulo" in cfg:
        log("")
        log("CONTROL de la firma -- incrementos BARAJADOS (debe salir plana):")
        log("   n      real      barajado    razon")
        for a, b in zip(cfg["H"]["filas"], cfg["firma_nulo"]):
            log("  %6d  %8.4f  %10.4f  %7.3f"
                % (a["n"], a["sigma_por_raiz_n"], b["sigma_por_raiz_n"],
                   a["sigma_por_raiz_n"] / max(b["sigma_por_raiz_n"], 1e-12)))
        r = cfg["H"]["filas"]
        n = np.array([f["n"] for f in r], float)
        yr = np.array([f["sigma_por_raiz_n"] for f in r])
        yb = np.array([f["sigma_por_raiz_n"] for f in cfg["firma_nulo"]])
        sel = n >= 256
        pr = np.polyfit(np.log(n[sel]), np.log(yr[sel]), 1)[0]
        pb = np.polyfit(np.log(n[sel]), np.log(yb[sel]), 1)[0]
        log("  pendiente [256, %d]:  real %+.4f   barajado %+.4f"
            % (int(n[-1]), pr, pb))
        log("  -> H_p real = %.3f   barajado = %.3f  (0.5 = difusivo)"
            % (pr + 0.5, pb + 0.5))
    p = cfg["part"]
    log("")
    log("particion 60/20/20 con embargo DESCARTADO:")
    log("  entrenamiento %8d   validacion %8d   prueba %8d"
        % (p["entrenamiento"].size, p["validacion"].size, p["prueba"].size))
    log("  embargo verificado por test: %s" % M.verificar_embargo(p))
    lb = 5 * cfg["embargo"]
    nb = int(np.ceil(p["prueba"].size / lb))
    log("  bloques de bootstrap en prueba (longitud 5*embargo = %d): %d %s"
        % (lb, nb, "" if nb >= 15 else "  *** SUBPOTENCIADO (<15) ***"))


def etapa_fuga(args):
    m0, cfg = preparar(args)
    cabecera_embargo(cfg)
    titulo("ETAPA 3 -- §4.2 TEST DE FUGA (4 variantes, todas al reporte)")
    log("Expectativas DECLARADAS ANTES: adelantada sube, retrasada baja,")
    log("BARAJADA cae al nivel de M0. Solo la cuarta tiene poder.")
    log("")

    K = cfg["K"]
    variantes = [("correcta", 0, False), ("adelantada +1", 1, False),
                 ("retrasada -10", -10, False), ("BARAJADA", 0, True)]
    filas = []
    for nombre, desp, baraj in variantes:
        m = m0 if (desp == 0 and not baraj) else construir_muestra(
            args.dir, desplazamiento=desp, barajar=baraj, semilla=7)
        n = len(m["y_mid"])
        part = M.particionar(n, cfg["embargo"])
        ent, val = part["entrenamiento"], part["validacion"]
        mod1 = M.ajustar(m["y_mid"][ent], m["eps"][ent], m["v"][ent], K, 0.5,
                         beta_libre=False)
        mod0 = M.ajustar_M0(m["y_mid"][ent])
        l1, _ = ll_media(mod1, m["y_mid"][val], m["eps"][val], m["v"][val], K)
        l0, _ = ll_media(mod0, m["y_mid"][val], m["eps"][val], m["v"][val], K)
        filas.append((nombre, l1, l0, l1 - l0, mod1["G0"]))
        log("  %-15s  LL/N(M1) = %+.6f   LL/N(M0) = %+.6f   dLL/N = %+.3e   G0 = %+.6f"
            % (nombre, l1, l0, l1 - l0, mod1["G0"]))

    base = filas[0][3]
    adel = filas[1][3]
    bar = filas[3][3]
    log("")
    log("LECTURA:")
    if adel <= base:
        log("  *** la fuga deliberada NO mejora -> el pipeline no usa el libro. SE PARA. ***")
    elif abs(adel - base) / max(abs(adel), 1e-12) < 0.02:
        log("  *** correcta y adelantada EMPATAN -> hay fuga en la correcta. SE PARA. ***")
    else:
        log("  [OK] adelantada supera a correcta por %.1fx: el codigo lee el campo"
            % (adel / max(base, 1e-12)))
    if bar > 0.25 * base:
        log("  *** la BARAJADA no cae al nivel de M0 (%.1f %% de la ganancia) ->"
            % (100 * bar / max(base, 1e-12)))
        log("      la ganancia no viene del emparejamiento libro-transaccion ***")
    else:
        log("  [OK] la BARAJADA cae al nivel de M0 (%.1f %% de la ganancia):"
            % (100 * bar / max(base, 1e-12)))
        log("       la ganancia SI viene del emparejamiento -- este es el control con poder")
    return filas


def etapa_delta(args):
    m, cfg = preparar(args)
    cabecera_embargo(cfg)
    titulo("ETAPA 4 -- §5.1 ELECCION DE delta EN VALIDACION (nunca en prueba)")
    K = cfg["K"]
    ent, val = cfg["part"]["entrenamiento"], cfg["part"]["validacion"]

    for etiq, ycol in (("PUNTO MEDIO (primario)", "y_mid"),
                       ("precio de transaccion (control)", "y_tr")):
        log("")
        log("--- observable: %s ---" % etiq)
        y = m[ycol]
        mod0 = M.ajustar_M0(y[ent])
        l0, _ = ll_media(mod0, y[val], m["eps"][val], m["v"][val], K)
        log("  M0                       LL/N = %+.6f" % l0)
        log("  delta modelo  LL/N validacion    dLL/N vs M0       G0        tau0"
            "     beta   f_inf       D")
        for delta in M.REJILLA_DELTA:
            for nom, kw in (("M1", dict(beta_libre=False)),
                            ("M2", dict(beta_libre=True, con_suelo=True))):
                mod = M.ajustar(y[ent], m["eps"][ent], m["v"][ent], K, delta, **kw)
                lv, _ = ll_media(mod, y[val], m["eps"][val], m["v"][val], K)
                D = M.decaimiento(mod)
                log("  %5.2f %-6s  %+.6f        %+.3e   %+.6f %9.2f %+7.3f %6.3f %8.4f"
                    % (delta, nom, lv, lv - l0, mod["G0"], mod["tau0"],
                       mod["beta"], mod["f_inf"], D))
        log("  ⚠ beta < 0 significa nucleo CRECIENTE (D > 1): ni permanente ni")
        log("    transitorio. La regla del §5.2 no contempla ese caso -- ver reporte.")


def etapa_A(args):
    m, cfg = preparar(args)
    cabecera_embargo(cfg)
    titulo("ETAPA 5 -- §5.2 TEST A: PERMANENTE CONTRA TRANSITORIO")
    K, ent, val = cfg["K"], cfg["part"]["entrenamiento"], cfg["part"]["validacion"]
    y = m["y_mid"]
    delta = args.delta

    mod = M.ajustar(y[ent], m["eps"][ent], m["v"][ent], K, delta, con_suelo=True)
    D = M.decaimiento(mod)
    log("ajuste M2 (con suelo) sobre ENTRENAMIENTO, delta = %.2f" % delta)
    log("  G0 = %+.6f   tau0 = %.2f   beta = %.4f   f_inf = %.4f   sigma = %.4f"
        % (mod["G0"], mod["tau0"], mod["beta"], mod["f_inf"], mod["sigma"]))
    log("  ⚠ beta y tau0 NO estan identificados individualmente (§5.2)")
    log("  D = G(K)/G(0) = %.6f    con K = %d ticks" % (D, K))

    perfil = M.decaimiento_perfil(mod)
    log("")
    log("perfil de decaimiento: D(K/4) = %.4f  D(K/2) = %.4f  D(K) = %.4f"
        % (perfil["D_0.25K"], perfil["D_0.50K"], perfil["D_1.00K"]))
    log("  caidas %.4f -> %.4f   se_aplana = %s"
        % (perfil["caida_primera"], perfil["caida_segunda"], perfil["se_aplana"]))

    # gamma de la autocorrelacion de signos, para emparejar el nulo (§5.2.quater)
    e = m["eps"][ent]
    gamma = float(np.corrcoef(e[:-1], e[1:])[0, 1])
    log("")
    log("gamma_signos medido en entrenamiento : %+.4f" % gamma)

    n_pru = cfg["part"]["prueba"].size
    log("regenerando q05 con N=%d, K=%d, gamma=%.3f, sigma=%.4f, G0=%.5f  (%d sorteos)"
        % (n_pru, K, max(gamma, 0.0), mod["sigma"], mod["G0"], args.sorteos))
    t0 = time.time()
    nl = M.nulo_de_D(n_pru, K, delta, mod["sigma"], mod["G0"], D_verdadero=1.0,
                     gamma_signos=max(gamma, 0.0), n_sorteos=args.sorteos,
                     semilla=11, con_suelo=True)
    log("  en %.0f s" % (time.time() - t0))
    log("  D_hat bajo D=1 verdadero: mediana %.4f  q05 %.4f  q95 %.4f  (sesgo %+.4f)"
        % (nl["D_mediana"], nl["q05"], nl["q95"], nl["sesgo"]))
    log("")
    log("CONTRASTE (§5.2): D_hat = %.4f  contra la distribucion simulada bajo D = 1"
        % D)
    log("                  q05 = %.4f      q95 = %.4f" % (nl["q05"], nl["q95"]))
    if D < nl["q05"]:
        log("  -> SE RECHAZA D = 1 POR ABAJO: IMPACTO TRANSITORIO.")
        log("  -> f_inf = %.4f  (fraccion permanente)  [sesgo negativo del estimador:"
            % mod["f_inf"])
        log("     infravalora el permanente -> sobrestima reversion (§5.2.ter)]")
        rechaza, lado = True, "transitorio"
    elif D > nl["q95"]:
        # ⚠ NO es una enmienda: es leer COMPLETA la misma distribucion simulada
        # que el §5.2 manda generar. La regla escrita es de una cola porque el
        # documento supone D en (0, 1]; el dato cayo fuera de ese espacio y
        # callarlo seria reportar "permanente" sobre evidencia de lo contrario.
        log("  *** D ESTA POR ENCIMA DE q95: se rechaza D = 1 POR ARRIBA. ***")
        log("  El nucleo ajustado CRECE (beta = %+.4f < 0): el impacto se AMPLIFICA"
            % mod["beta"])
        log("  con el rezago. Eso no es permanente NI transitorio -- es un tercer")
        log("  caso que la dicotomia del §5.2 no contempla. La regla literal diria")
        log("  'no se rechaza D=1 -> permanente' por ser de una cola, y seria falso.")
        log("  f_inf NO se reporta: no esta identificado con beta < 0.")
        rechaza, lado = True, "creciente"
    else:
        log("  -> NO se rechaza D = 1: IMPACTO PERMANENTE. MkII bien especificada.")
        log("  -> f_inf NO se reporta (§5.2.ter): sin decaimiento no esta identificado.")
        rechaza, lado = False, "permanente"
    log("")
    log("⚠ El paso 4 NO es una prueba economica (§5.2.bis). El paso 3 lo es.")

    # R2 descriptivo (§2.3: se reporta, NO decide)
    for nom, mm in (("M1", M.ajustar(y[ent], m["eps"][ent], m["v"][ent], K, delta,
                                     beta_libre=False)), ("M2", mod)):
        pr = M.predecir(mm, m["eps"][val], m["v"][val])[K:]
        yy = y[val][K:]
        r2 = 1.0 - float(np.sum((yy - pr) ** 2) / np.sum((yy - yy.mean()) ** 2))
        log("  R2 descriptivo de %s en validacion = %+.5f   (NO decide, §2.3)"
            % (nom, r2))
    log("  ⚠ y_mid es 96.7 %% ceros exactos: la verosimilitud gaussiana esta mal")
    log("    especificada en la marginal. Los dLL/N son comparables entre modelos")
    log("    (misma familia, misma muestra) pero su nivel no es interpretable.")

    # residuo del modelo elegido
    titulo("§5.4 RESIDUO (Ljung-Box por magnitud, umbral |rho1| < 0.05)")
    for nom, kw in (("M0", None), ("M1", dict(beta_libre=False)),
                    ("M2", dict(beta_libre=True, con_suelo=True))):
        if kw is None:
            mm = M.ajustar_M0(y[ent])
        else:
            mm = M.ajustar(y[ent], m["eps"][ent], m["v"][ent], K, delta, **kw)
        r = y[val] - M.predecir(mm, m["eps"][val], m["v"][val])
        r = r[K:]
        rr = rho_lags(r, [1, 2, 3, 10, 100, cfg["embargo"], 2 * cfg["embargo"]])
        log("  %-4s rho1 %+.4f  rho2 %+.4f  rho3 %+.4f  rho10 %+.4f  rho100 %+.4f  rho_e %+.4f"
            % (nom, rr.get(1, np.nan), rr.get(2, np.nan), rr.get(3, np.nan),
               rr.get(10, np.nan), rr.get(100, np.nan), rr.get(cfg["embargo"], np.nan)))
    log("")
    log("  Declarado: un modelo con residuo estructurado NO gana aunque tenga mejor LL/N.")
    return {"D": D, "q05": nl["q05"], "rechaza": rechaza, "mod": mod}


def etapa_osc(args):
    m, cfg = preparar(args)
    titulo("ETAPA 6 -- §5.3 omega_G CONTRA NULO SIMULADO (nunca chi2)")
    log("⚠ omega_G NO es omega_m resucitada: es una frecuencia de la RESPUESTA AL")
    log("  IMPULSO DEL FLUJO, con estimador propio y nulo simulado. omega_m era una")
    log("  frecuencia del PRECIO por EMD, sin nulo, y resulto ser la de la ventana.")
    log("")
    K = cfg["K"]
    ent = cfg["part"]["entrenamiento"]
    # Submuestra contigua para que el contraste sea ejecutable en tiempo finito.
    n_sub = min(args.n_osc, ent.size)
    sub = ent[:n_sub]
    log("submuestra contigua de entrenamiento: %d obs, %d sorteos" % (n_sub, args.sorteos))
    t0 = time.time()
    r = M.contraste_omega_G(m["y_mid"][sub], m["eps"][sub], m["v"][sub], K,
                            args.delta, n_sorteos=args.sorteos, semilla=13)
    log("  en %.0f s" % (time.time() - t0))
    log("  2*dLL observado = %.3f      nulo p95 = %.3f      p simulado = %.4f"
        % (r["estadistico"], r["nulo_p95"], r["p_simulado"]))
    log("  omega_G ajustado = %+.6f rad/tick" % r["omega_G"])
    if r["p_simulado"] < 0.05:
        log("  -> omega_G != 0: el sistema es un OSCILADOR FORZADO.")
    else:
        log("  -> omega_G = 0: el propagador es MONOTONO.")
        log("     Consecuencia declarada: omega_m,max y gamma_omega se borran de")
        log("     constantes_micelio.py.")
    return r


def etapa_B(args):
    m, cfg = preparar(args)
    titulo("ETAPA 7 -- §6 TEST B: ¿aporta el lado de LIMITE? y variantes c_t")
    K, ent, val = cfg["K"], cfg["part"]["entrenamiento"], cfg["part"]["validacion"]
    y, delta = m["y_mid"], args.delta

    # M2 de referencia
    m2 = M.ajustar(y[ent], m["eps"][ent], m["v"][ent], K, delta, con_suelo=True)
    l2, _ = ll_media(m2, y[val], m["eps"][val], m["v"][val], K)
    log("M2 (solo lado taker)            LL/N = %+.6f" % l2)

    # M2+L: el forzamiento es el OFI-L1 completo, que incluye el lado de limite.
    # Se codifica pasando eps = sgn(e) y v = |e| con delta = 1, o sea x_t = e_t.
    e = m["e"]
    epsL, vL = np.sign(e), np.abs(e)
    vL[vL <= 0] = 1e-12
    m2L = M.ajustar(y[ent], epsL[ent], vL[ent], K, 1.0, con_suelo=True)
    l2L, _ = ll_media(m2L, y[val], epsL[val], vL[val], K)
    log("M2+L (OFI-L1 aproximado)        LL/N = %+.6f   dLL/N = %+.3e"
        % (l2L, l2L - l2))

    lb = 5 * cfg["embargo"]
    _, llA = ll_media(m2, y[val], m["eps"][val], m["v"][val], K)
    _, llB = ll_media(m2L, y[val], epsL[val], vL[val], K)
    n = min(len(llA), len(llB))
    bs = M.bootstrap_bloques(llB[:n] - llA[:n], lb, 400, semilla=5)
    marg = M.margen_bic(0, n, bs["n_bloques"])
    log("  IC95 bootstrap de dLL/N = [%+.3e, %+.3e]  bloques = %d%s"
        % (bs["ic"][0], bs["ic"][1], bs["n_bloques"],
           "  *** SUBPOTENCIADO ***" if bs["subpotenciado"] else ""))
    if bs["ic"][0] > marg:
        log("  -> el lado de LIMITE aporta; @depth y el anti-spoofing se justifican")
    else:
        log("  -> el lado de LIMITE NO aporta: @depth y el anti-spoofing de MkII")
        log("     son coste sin retorno. Este resultado se publica igual (§6).")

    # --- §6.2 variantes de c_t sobre G0
    log("")
    log("--- §6.2 variantes de G0 (V_best, nunca 'profundidad') ---")
    bloque = 256
    nb = len(y) // bloque
    sig = np.repeat(np.array([np.std(y[i * bloque:(i + 1) * bloque]) for i in range(nb)]),
                    bloque)
    sig = np.concatenate([sig, np.full(len(y) - len(sig), sig[-1] if len(sig) else 1.0)])
    sig = sig / max(np.median(sig[sig > 0]), 1e-12)
    Vb = np.maximum(m["V_best"], 1e-9)
    Vb = Vb / np.median(Vb)
    for nom, mult in (("G0 constante", np.ones(len(y))),
                      ("G0 -> Y*sigma_t", sig),
                      ("G0 -> Y/sqrt(V_best)", 1.0 / np.sqrt(Vb)),
                      ("G0 -> Y*sigma_t/sqrt(V_best)", sig / np.sqrt(Vb))):
        vv = m["v"] * mult ** (1.0 / max(delta, 1e-9)) if delta > 0 else m["v"]
        ee = m["eps"] * (mult if delta == 0 else 1.0)
        mm = M.ajustar(y[ent], ee[ent], vv[ent], K, delta, con_suelo=True)
        lv, _ = ll_media(mm, y[val], ee[val], vv[val], K)
        log("  %-30s LL/N = %+.6f   dLL/N vs constante = %+.3e"
            % (nom, lv, lv - l2))
    log("  Declarado (§6.2): si Y*sigma_t sola captura toda la ganancia, M1' NO mide")
    log("  microestructura -- reescala por volatilidad, y eso ya lo hace rho_k.")


def etapa_decision(args):
    m, cfg = preparar(args)
    cabecera_embargo(cfg)
    titulo("ETAPA 8 -- §9 REGLA DE DECISION (se para en el primer fallo)")
    log("⚠ Aqui se abre el conjunto de PRUEBA, y se abre UNA sola vez.")
    K, delta = cfg["K"], args.delta
    ent, val, pru = (cfg["part"]["entrenamiento"], cfg["part"]["validacion"],
                     cfg["part"]["prueba"])
    y = m["y_mid"]

    # --- Paso 1: compuerta (ya evaluada por --resumen; se repite el numero)
    log("")
    log("PASO 1 -- compuerta de datos: %d ticks continuos contra %d exigidos  [PASA]"
        % (cfg["n"], M.ticks_de_compuerta()))

    # --- Paso 2: LL/N del mejor modelo contra M0, fuera de muestra
    mods = {"M0": M.ajustar_M0(y[ent]),
            "M1": M.ajustar(y[ent], m["eps"][ent], m["v"][ent], K, 0.5, beta_libre=False),
            "M2": M.ajustar(y[ent], m["eps"][ent], m["v"][ent], K, delta, con_suelo=True)}
    lls, ll_v = {}, {}
    for nom, mm in mods.items():
        lls[nom], ll_v[nom] = ll_media(mm, y[pru], m["eps"][pru], m["v"][pru], K)
    log("")
    log("PASO 2 -- LL/N FUERA DE MUESTRA (conjunto de prueba, N = %d):" % pru.size)
    for nom in ("M0", "M1", "M2"):
        log("   %-4s  LL/N = %+.8f" % (nom, lls[nom]))

    lb = 5 * cfg["embargo"]
    mejor = max(("M1", "M2"), key=lambda k: lls[k])
    n = min(len(ll_v[mejor]), len(ll_v["M0"]))
    bs = M.bootstrap_bloques(ll_v[mejor][:n] - ll_v["M0"][:n], lb, 500, semilla=3)
    log("   mejor = %s   dLL/N vs M0 = %+.3e   IC95 = [%+.3e, %+.3e]"
        % (mejor, bs["media"], bs["ic"][0], bs["ic"][1]))
    log("   bloques efectivos = %d%s" % (bs["n_bloques"],
        "   *** SUBPOTENCIADO (<15) ***" if bs["subpotenciado"] else ""))
    if bs["ic"][0] <= 0.0:
        log("   [FALLA] el IC no excluye 0 -> GANA M0. Se para en el paso 2.")
        return {"paso": 2, "veredicto": "M0"}
    log("   [PASA] el IC excluye 0")

    # --- Paso 3: criterio economico (§9.1) -- la unica compuerta con dinero
    H = int(round(cfg["H"]["H_ticks"]))
    c_u = cfg["c_u"]
    log("")
    log("PASO 3 -- CRITERIO ECONOMICO (§9.1), horizonte H*_ticks = %d, c(u) = %.4f"
        % (H, c_u))
    res3 = {}
    for nom in ("M1", "M2"):
        mu = mu_horizonte(mods[nom], m["eps"][pru], m["v"][pru], H)[K:]
        am = np.abs(mu)
        q90 = float(np.percentile(am, 90))
        frac = float(np.mean(am >= 1.5 * c_u))
        umbral_frac = 20.0 * cfg["embargo"] / max(pru.size, 1)
        res3[nom] = (q90, frac, umbral_frac)
        log("   %-4s q90(|mu|) = %.6f   contra 1.5*c(u) = %.4f   %s"
            % (nom, q90, 1.5 * c_u, "PASA" if q90 >= 1.5 * c_u else "falla"))
        log("        frac(|mu| >= 1.5c) = %.3e   contra %.3e   %s"
            % (frac, umbral_frac, "PASA" if frac >= umbral_frac else "falla"))
        # DESCRIPTIVOS -- ninguna regla depende de ellos; dicen CUANTO falta.
        log("        |mu|: q50 %.4f  q90 %.4f  q99 %.4f  max %.4f   [descriptivo]"
            % (np.percentile(am, 50), q90, np.percentile(am, 99), am.max()))
        if q90 > 0:
            c_eq = q90 / 1.5
            log("        c(u) que 3a exigiria para pasar: %.4f USD/BTC = %.3f pb"
                % (c_eq, 1e4 * c_eq / (2.0 * cfg["precio"])))
            log("        VIP 0 asumido: %.3f pb por lado -> falta un factor %.2f"
                % (1e4 * 0.0002, 1.5 * c_u / q90))
    q90, frac, uf = res3[mejor]
    a3 = q90 >= 1.5 * c_u
    b3 = frac >= uf
    if a3 and b3:
        log("   [PASA] 3a y 3b")
    elif a3 and not b3:
        log("   [NO DECIDIBLE] 3a pasa y 3b falla -> §9.2: es carencia de datos,")
        log("   NO abandono. La unica accion admisible es capturar mas.")
        return {"paso": 3, "veredicto": "NO DECIDIBLE"}
    else:
        log("   [FALLA] 3a falla%s -> GANA M0." % (" con 3b satisfecha" if b3 else ""))
        if b3:
            log("   §11.2: queda abandonada la hipotesis de estructura explotable a")
            log("   ESCALA DE SEGUNDOS. NO decide nada sobre minutos a horas (§11.1).")
        else:
            log("   3b NO satisfecha: sin potencia para haber detectado la senal, asi")
            log("   que el §11 NO se activa (§9.2).")
        return {"paso": 3, "veredicto": "M0"}

    # --- Paso 4: M2 contra M1
    log("")
    log("PASO 4 -- M2 contra M1 por el margen del §10")
    n = min(len(ll_v["M2"]), len(ll_v["M1"]))
    bs4 = M.bootstrap_bloques(ll_v["M2"][:n] - ll_v["M1"][:n], lb, 500, semilla=4)
    marg_a = M.margen_bic(4, n, bs4["n_bloques"])
    marg_b = 0.1 * c_u
    q90_2, q90_1 = res3["M2"][0], res3["M1"][0]
    log("   (a) dLL/N IC95 = [%+.3e, %+.3e]  contra k*ln(N_eff)/(2N) = %.3e -> %s"
        % (bs4["ic"][0], bs4["ic"][1], marg_a, "PASA" if bs4["ic"][0] > marg_a else "falla"))
    log("   (b) q90(M2) - q90(M1) = %+.6f  contra 0.1*c(u) = %.4f -> %s"
        % (q90_2 - q90_1, marg_b, "PASA" if (q90_2 - q90_1) >= marg_b else "falla"))
    if bs4["ic"][0] > marg_a and (q90_2 - q90_1) >= marg_b:
        log("   -> GANA M2 (propagador)")
        return {"paso": 4, "veredicto": "M2"}
    log("   -> GANA M1 (MkII) por parsimonia. (b) es el criterio vinculante (§10.2).")
    return {"paso": 4, "veredicto": "M1"}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--etapa", default="muestra")
    ap.add_argument("--dir", default=DIR_CAPTURA)
    ap.add_argument("--delta", type=float, default=0.5)
    ap.add_argument("--sorteos", type=int, default=40)
    ap.add_argument("--n-osc", dest="n_osc", type=int, default=120000)
    a = ap.parse_args(argv)

    etapas = {"muestra": etapa_muestra, "fuga": etapa_fuga, "delta": etapa_delta,
              "A": etapa_A, "osc": etapa_osc, "B": etapa_B,
              "decision": etapa_decision}
    if a.etapa == "todo":
        for k in ("muestra", "fuga", "delta", "A", "osc", "B", "decision"):
            etapas[k](a)
        return 0
    if a.etapa not in etapas:
        log("etapa desconocida: %s" % a.etapa)
        return 2
    etapas[a.etapa](a)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
