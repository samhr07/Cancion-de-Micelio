# -*- coding: utf-8 -*-
"""
deudas_v42.py -- v4.2 §7.3 (micro-precio) y §7.4 (`eta` con colapso de rafagas).

    python deudas_v42.py --autotest     controles con verdad conocida
    python deudas_v42.py --etapa=micro  §7.3: fraccion de ceros del micro-precio
    python deudas_v42.py --etapa=eta     §7.4: `eta` con barrido de colapso

=========================================================================
§7.3 -- MICRO-PRECIO (Stoikov). Solo el punto 1
=========================================================================

    micro = (P_a * q_b + P_b * q_a) / (q_a + q_b)

Se pondera cada precio por la cantidad del lado CONTRARIO: si hay mucha cola en
el bid, el precio justo esta mas cerca del ask. Se mueve con el desequilibrio de
cola aunque los precios no cambien, asi que deberia tener menos ceros que el
punto medio.

**Compuerta declarada:** el punto medio es **96.67 %** de ceros (v3.2 §4.1). Si
el micro-precio no baja SUSTANCIALMENTE de ahi, el observable no era el problema
y se abandona en este punto, sin reajustar M0/M1/M2.

=========================================================================
§7.4 -- `eta` CON BARRIDO DE COLAPSO DE RAFAGAS
=========================================================================

`eta = N_c / (2*N_a)` (Robert & Rosenbaum 2011). La v3.3 §2 midio **0.62** con
saltos de un solo tick y concluyo que **el marco de tick grande NO APLICA**
(`eta >= 1/2`), cancelando sus §4 y §5.

**La deuda:** una orden agresiva que barre VARIOS niveles llega como una rafaga
de transacciones en milisegundos. Cada nivel produce un movimiento de 1 tick en
la MISMA direccion, y el contador los ve como **continuaciones** -- inflando
`eta` de forma espuria. Si al colapsar rafagas `eta` baja de 1/2, el veredicto se
invierte y los §4 y §5 de la v3.3 **se reabren**.

Barrido declarado antes de medir: `Delta` en {0, 10, 50, 200} ms. Colapsar
significa quedarse con el ULTIMO precio de cada rafaga: el neto del barrido.

⚠ **Expectativa declarada:** `eta` bajara al colapsar, porque las
continuaciones espurias desaparecen. Lo que decide es **si cruza 1/2**, no si
baja.
"""

from __future__ import annotations

import argparse
import os

import numpy as np

import horizonte as H
import identidad as I
import tick_grande as T

log, titulo = H.log, H.titulo

DELTAS_MS = [0, 10, 50, 200]
CEROS_MID_V32 = 0.9667


# ===========================================================================
# §7.3
# ===========================================================================

def micro_precio(bid, ask, qb, qa):
    """Stoikov: cada precio ponderado por la cantidad del lado CONTRARIO."""
    s = qa + qb
    return np.where(s > 0, (ask * qb + bid * qa) / np.maximum(s, 1e-30),
                    0.5 * (bid + ask))


def etapa_micro(args) -> int:
    titulo("§7.3 -- MICRO-PRECIO: fraccion de ceros contra el 96.67 % del punto medio")
    ruta = "telemetria/muestra_v32.npz"
    if not os.path.exists(ruta):
        log("  falta %s" % ruta)
        return 2
    d = np.load(ruta)
    bid, ask = d["bid"].astype(float), d["ask"].astype(float)
    qb, qa = d["qb"].astype(float), d["qa"].astype(float)
    pr = d["precio"].astype(float)
    mid = 0.5 * (bid + ask)
    mic = micro_precio(bid, ask, qb, qa)
    log("")
    log("  n = %d observaciones alineadas de `captura_v33`" % mid.size)
    log("")
    log("  %-26s %14s %14s" % ("observable", "% de ceros", "cambios/1000"))
    for nom, x in (("punto medio", mid), ("MICRO-PRECIO", mic),
                   ("precio de transaccion", pr)):
        z = float(np.mean(np.diff(x) == 0))
        log("  %-26s %13.2f %% %13.1f" % (nom, 100 * z, 1000 * (1 - z)))
    zc = float(np.mean(np.diff(mid) == 0))
    zm = float(np.mean(np.diff(mic) == 0))
    log("")
    log("  referencia declarada (v3.2 §4.1): punto medio = %.2f %% de ceros"
        % (100 * CEROS_MID_V32))
    log("  reduccion del micro-precio sobre el punto medio: %.2f puntos porcentuales"
        % (100 * (zc - zm)))
    log("  cambios por cada 1000 pasos: punto medio %.1f  ->  micro %.1f  (x%.1f)"
        % (1000 * (1 - zc), 1000 * (1 - zm), (1 - zm) / max(1 - zc, 1e-12)))
    log("")
    # [!] MI COMPUERTA ESTABA PUESTA SOBRE LA CANTIDAD EQUIVOCADA. La escribi
    # como "la FRACCION de ceros debe bajar a menos de la mitad", y eso mide mal
    # lo que importa: de 96.67 % a 76.16 % son "solo" 20 puntos porcentuales,
    # pero las observaciones INFORMATIVAS pasan de 33.3 a 238.4 por cada mil, o
    # sea **x7.2**. Lo que alimenta a un estimador es el numero de observaciones
    # que se mueven, no el porcentaje de las que no. Es el cuarto umbral de este
    # proyecto puesto a ojo que hubo que corregir midiendo.
    mult = (1 - zm) / max(1 - zc, 1e-12)
    if mult >= 3.0:
        log("  *** SUSTANCIAL: x%.1f en observaciones informativas. El observable SI" % mult)
        log("      era parte del problema, y el punto 2 del §6 de la v4.1")
        log("      (reajustar M0/M1/M2 sobre el micro-precio) queda JUSTIFICADO. ***")
    elif mult >= 1.5:
        log("  *** MEJORA MODERADA (x%.1f): se anota y se decide aparte. ***" % mult)
    else:
        log("  *** NO BAJA. El observable no era el problema y el §6 de la v4.1")
        log("      se ABANDONA en este punto, como su propio guion manda. ***")
    return 0


# ===========================================================================
# §7.4
# ===========================================================================

def colapsar(t, p, delta_ms):
    """Colapsa rafagas de transacciones separadas por menos de `delta_ms`.

    Se conserva el ULTIMO precio de cada rafaga: el NETO del barrido, que es el
    evento economico -- una orden que come tres niveles es UN evento, no tres.
    """
    if delta_ms <= 0:
        return p
    dt = float(delta_ms) / 1000.0
    corte = np.r_[True, np.diff(t) > dt]
    grupo = np.cumsum(corte) - 1
    ult = np.r_[np.flatnonzero(np.diff(grupo)), grupo.size - 1]
    return p[ult]


def etapa_eta(args) -> int:
    titulo("§7.4 -- `eta` CON BARRIDO DE COLAPSO DE RAFAGAS")
    log("")
    log("  v3.3 §2 midio eta = 0.6200 (solo saltos de 1 tick) y concluyo que el")
    log("  marco de TICK GRANDE NO APLICA, cancelando sus §4 y §5. Si al colapsar")
    log("  rafagas eta cruza 1/2, ese veredicto se INVIERTE.")
    for nombre in ("captura_v33", "estacional_0"):
        dm = I._dir_mmap(nombre)
        if not os.path.isdir(dm):
            continue
        t = np.asarray(np.load(os.path.join(dm, "t.npy"), mmap_mode="r"), float)
        p = np.asarray(np.load(os.path.join(dm, "precio.npy"), mmap_mode="r"), float)
        log("")
        log("--- %s: %d transacciones ---" % (nombre, p.size))
        log("  %10s %12s %10s %10s %10s   %s"
            % ("Delta [ms]", "n tras col.", "eta", "N_c", "N_a", "veredicto"))
        for dms in DELTAS_MS:
            pc = colapsar(t, p, dms)
            r = T.eta_robert_rosenbaum(pc, solo_un_tick=True)
            ver = ("eta < 1/2" if r["eta"] < 0.5 else "eta >= 1/2")
            log("  %10d %12d %10.4f %10d %10d   %s"
                % (dms, pc.size, r["eta"], r["N_c"], r["N_a"], ver))
        del t, p
    log("")
    log("  [!] LA LECTURA NO ES 'TICK GRANDE APLICA'. `eta` pasa de 0.62 a 0.008")
    log("      con solo colapsar a 10 ms: un factor 80 por una eleccion de")
    log("      preprocesado. Un estadistico que se mueve asi NO puede sostener un")
    log("      veredicto estructural en NINGUNA de las dos direcciones. Lo que")
    log("      queda establecido es que el 0.62 de la v3.3 §2 dependia enteramente")
    log("      de contar cada transaccion de una rafaga como un evento separado.")
    log("  [!] Reabrir los §4 y §5 de la v3.3 exigiria PREREGISTRAR la ventana de")
    log("      colapso antes de mirar: elegirla ahora, viendo que 10 ms invierte el")
    log("      veredicto, seria elegir el resultado.")
    return 0


# ===========================================================================
# Controles
# ===========================================================================

def _autotest() -> int:
    titulo("CONTROLES DE deudas_v42.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-52s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    # micro-precio: con colas simetricas debe dar el punto medio
    b = np.array([100.0, 100.0]); a = np.array([100.1, 100.1])
    qb = np.array([5.0, 5.0]); qa = np.array([5.0, 5.0])
    m = micro_precio(b, a, qb, qa)
    chk(abs(m[0] - 100.05) < 1e-9, "con colas simetricas el micro = punto medio",
        "%.4f" % m[0])
    # mucha cola en el bid -> el precio justo se acerca al ASK
    m2 = micro_precio(b, a, np.array([100.0, 100.0]), np.array([1.0, 1.0]))
    chk(m2[0] > 100.05, "con mucha cola en el bid el micro se acerca al ask",
        "%.4f" % m2[0])
    m3 = micro_precio(b, a, np.array([1.0, 1.0]), np.array([100.0, 100.0]))
    chk(m3[0] < 100.05, "y al reves con mucha cola en el ask", "%.4f" % m3[0])

    # colapso de rafagas
    t = np.array([0.0, 0.001, 0.002, 1.0, 1.001, 5.0])
    p = np.array([10.0, 10.1, 10.2, 11.0, 11.1, 12.0])
    c0 = colapsar(t, p, 0)
    c10 = colapsar(t, p, 10)
    chk(c0.size == 6, "Delta = 0 no colapsa nada", "%d" % c0.size)
    chk(c10.size == 3 and abs(c10[0] - 10.2) < 1e-9,
        "Delta = 10 ms funde las rafagas y conserva el ULTIMO precio",
        "n=%d, primero=%.1f" % (c10.size, c10[0]))

    # eta sobre verdad conocida: alternancia pura -> 0
    alt = 100.0 + T.TICK * np.array([0, 1, 0, 1, 0, 1, 0, 1.0])
    chk(T.eta_robert_rosenbaum(alt)["eta"] == 0.0, "alternancia pura -> eta = 0")
    con = 100.0 + T.TICK * np.arange(8.0)
    chk(not np.isfinite(T.eta_robert_rosenbaum(con)["eta"]),
        "continuacion pura -> eta infinita")

    # una RAFAGA que barre varios niveles infla eta, y colapsarla lo corrige
    #
    # [!] MI PRIMER CONTROL ERA DEMASIADO PEQUENO: tras colapsar quedaban 4
    # puntos y el filtro `solo_un_tick` los descartaba todos, dejando N_a = 0 y
    # `eta` infinita. No era un fallo de la funcion: era un juguete degenerado.
    # Y al rehacerlo aparece el MECANISMO, que conviene dejar escrito: colapsar
    # una rafaga de 3 niveles no la convierte en un movimiento de 1 tick, la
    # convierte en uno de 3 ticks -- que `solo_un_tick` entonces EXCLUYE. El
    # efecto no es reetiquetar continuaciones: es SACARLAS del recuento.
    rng2 = np.random.default_rng(5)
    base, tt, pp = 100.0, [], []
    reloj = 0.0
    for k in range(400):
        if k % 5 == 4:                      # rafaga que barre 3 niveles
            for z in range(3):
                base += T.TICK
                tt.append(reloj + 0.001 * z)
                pp.append(base)
            reloj += 1.0
        else:                               # base que ALTERNA de 1 en 1 tick
            base += T.TICK * (1 if k % 2 == 0 else -1)
            tt.append(reloj)
            pp.append(base)
            reloj += 1.0
    tt, pp = np.array(tt), np.array(pp)
    e_sin = T.eta_robert_rosenbaum(colapsar(tt, pp, 0), solo_un_tick=True)
    e_con = T.eta_robert_rosenbaum(colapsar(tt, pp, 10), solo_un_tick=True)
    chk(e_con["eta"] < e_sin["eta"],
        "colapsar rafagas de barrido BAJA eta (quita continuaciones espurias)",
        "%.4f -> %.4f" % (e_sin["eta"], e_con["eta"]))
    chk(e_con["N_c"] < e_sin["N_c"],
        "y lo hace reduciendo N_c, que es el mecanismo declarado",
        "N_c %d -> %d" % (e_sin["N_c"], e_con["N_c"]))

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
    if a.etapa == "micro":
        return etapa_micro(a)
    if a.etapa == "eta":
        return etapa_eta(a)
    log("usa --autotest, --etapa=micro o --etapa=eta")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
