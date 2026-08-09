# -*- coding: utf-8 -*-
"""
cola.py -- que se puede saber de la ejecucion PASIVA sin poner una sola orden.

    python cola.py --dir=telemetria/captura_v32
    python cola.py --autotest

⚠ ESTO NO ES LA v3.2 Y NO LA TOCA. La decision M0/M1/M1'/M2 esta congelada por
`PREREGISTRO_3_2.md` hasta que la compuerta pase, y se hara sobre `captura_v33`.
Lo de aqui es otra pregunta -- si entrase pasivo, ¿me llenarian, cuando, y a que
precio se movio el mercado despues? -- y se mide sobre `captura_v32`, que quedo
partida por la suspension y ya estaba apartada como banco de pruebas. Ningun
resultado de este modulo entra en ninguna regla de decision de la v3.2.

LAS CUATRO MEDICIONES
---------------------
1. **Supervivencia de la cola.** Insercion hipotetica en `t` al final de la cola
   del bid: ¿cuando me llenan? Da `P(sin llenar tras tau | Q)` sin operar.
2. **Markout sintetico.** En cada llenado hipotetico, ¿que hace el mid a 1, 5 y
   50 s? Es la seleccion adversa que habria sufrido estando en esa cola.
3. **Flujos de la cola.** Cuanto de la desaparicion de cola es transaccion
   observada y cuanto no. ⚠ Con L1 **cancelacion y reposicion no se separan**,
   solo su neto, asi que se reporta lo identificable y se nombra lo que no.
4. **Horquillado de la posicion en cola.** Con L1 se ve bajar `B` sin saber si
   fue transaccion (visible en `@trade`) o cancelacion (invisible). Atribuyendo
   TODAS las cancelaciones a delante sale una cota optimista, y todas a detras
   una pesimista. El ancho del horquillado dice cuanta incertidumbre hay ANTES
   de gastar en `@depth`.

CONVENCION DE LADO
------------------
`m = True` -> el comprador es maker -> el agresor es el VENDEDOR -> la
transaccion consume el BID. Es la misma convencion del PREREGISTRO_3_1 §4, ya
verificada alli (799/799 rezagos con el signo correcto).
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

import captura_larga as cl


TICK = 0.10          # tickSize de BTCUSDT en Binance futuros (medido en la v1.3)
TOL = TICK / 2.0     # tolerancia para casar el precio de un trade con un nivel


def construir_eventos(d: dict) -> dict:
    """Funde libro y transacciones en una linea temporal unica y ordenada.

    Se ordena por tiempo con `kind='stable'` y con el libro ANTES que el trade a
    igualdad de marca: el estado del libro utilizable para una transaccion es el
    inmediatamente anterior, que es la misma regla del §3.1 de la v3.2.
    """
    bt, bb, bB = d["bk_t"], d["bk_b"], d["bk_B"]
    ba, bA = d["bk_a"], d["bk_A"]
    tt, tp, tq = d["tr_t"], d["tr_precio"], d["tr_cant"]
    tm = d["tr_maker"].astype(bool)

    n_b, n_t = len(bt), len(tt)
    t_all = np.concatenate([bt, tt])
    # 0 = actualizacion de libro, 1 = transaccion. El desempate favorece al libro.
    tipo = np.concatenate([np.zeros(n_b, np.int8), np.ones(n_t, np.int8)])
    orden = np.lexsort((tipo, t_all))
    return {"t": t_all[orden], "tipo": tipo[orden],
            "idx": np.concatenate([np.arange(n_b), np.arange(n_t)])[orden],
            "bk_t": bt, "bk_b": bb, "bk_B": bB, "bk_a": ba, "bk_A": bA,
            "tr_t": tt, "tr_p": tp, "tr_q": tq, "tr_maker": tm,
            "n_b": n_b, "n_t": n_t}


def _mid(ev: dict, i_b: int) -> float:
    return 0.5 * (ev["bk_b"][i_b] + ev["bk_a"][i_b])


def simular_insercion(ev: dict, k0: int, horizonte_s: float = 60.0) -> dict:
    """Una insercion hipotetica al FINAL de la cola del bid en el evento `k0`.

    Devuelve el desenlace y, cuando lo hay, el instante de llenado bajo las dos
    atribuciones del horquillado (§4).

    Desenlaces:
      LLENADO      la cola por delante se agoto y me tocaria
      NIVEL_BARRIDO el bid bajo de p0 sin haberme llenado -- el precio se fue
      SUPERADO     el bid subio por encima de p0: quedo detras del mejor nivel
      CENSURADO    se acabo el horizonte sin desenlace
    """
    t = ev["t"]; tipo = ev["tipo"]; idx = ev["idx"]
    # Estado del libro en k0: hace falta que k0 sea una actualizacion de libro.
    i_b = int(idx[k0])
    p0 = float(ev["bk_b"][i_b])
    Q0 = float(ev["bk_B"][i_b])
    if not np.isfinite(p0) or p0 <= 0 or Q0 <= 0:
        return {"desenlace": "INVALIDO"}

    t0 = float(t[k0])
    t_fin = t0 + horizonte_s
    consumido = 0.0      # volumen transado contra el bid a p0
    cancelado = 0.0      # bajadas de B no explicadas por transacciones
    B_prev = Q0
    p_prev = p0
    fill_opt = None      # todas las cancelaciones DELANTE (optimista)
    fill_pes = None      # todas las cancelaciones DETRAS (pesimista)

    k = k0 + 1
    n = len(t)
    while k < n and t[k] <= t_fin:
        if tipo[k] == 1:                                  # transaccion
            j = int(idx[k])
            if ev["tr_maker"][j] and abs(ev["tr_p"][j] - p0) <= TOL:
                consumido += float(ev["tr_q"][j])
        else:                                             # actualizacion de libro
            i = int(idx[k])
            p_now = float(ev["bk_b"][i])
            B_now = float(ev["bk_B"][i])
            if abs(p_now - p0) <= TOL:
                # ⚠ EL HORQUILLADO, planteado sin inventar cantidades.
                #
                # PESIMISTA: toda cancelacion ocurre DETRAS de mi, asi que solo
                # me adelanta lo transado -> cola por delante = Q0 - consumido.
                #
                # OPTIMISTA: toda cancelacion ocurre DELANTE. Entonces mi
                # posicion no puede ser peor que la cantidad TOTAL exhibida:
                # si solo quedan B_now BTC en el nivel, como mucho tengo eso por
                # delante -> cola por delante = min(B_now, Q0 - consumido).
                #
                # Asi el optimista sale de una cantidad OBSERVADA (`B_now`) y no
                # de una "cancelacion" estimada, que con L1 no es identificable.
                cancelado = max(0.0, (Q0 - consumido) - B_now)
                B_prev = B_now
            elif p_now < p0 - TOL:
                return _cierre("NIVEL_BARRIDO", t, k, t0, p0, Q0, consumido,
                               cancelado, fill_opt, fill_pes, ev, i)
            else:
                return _cierre("SUPERADO", t, k, t0, p0, Q0, consumido,
                               cancelado, fill_opt, fill_pes, ev, i)
            p_prev = p_now

        # ¿Me tocaria ya, bajo cada atribucion?
        if fill_pes is None and consumido >= Q0:
            fill_pes = float(t[k])
        if fill_opt is None and (consumido + cancelado) >= Q0:
            fill_opt = float(t[k])
        if fill_pes is not None:
            return _cierre("LLENADO", t, k, t0, p0, Q0, consumido, cancelado,
                           fill_opt, fill_pes, ev, None)
        k += 1

    return _cierre("CENSURADO", t, min(k, n - 1), t0, p0, Q0, consumido,
                   cancelado, fill_opt, fill_pes, ev, None)


def _cierre(desenlace, t, k, t0, p0, Q0, consumido, cancelado,
            fill_opt, fill_pes, ev, i_b):
    return {"desenlace": desenlace, "t0": t0, "t_fin": float(t[k]),
            "duracion": float(t[k] - t0), "p0": p0, "Q0": Q0,
            "consumido": consumido, "cancelado": cancelado,
            "t_fill_optimista": fill_opt, "t_fill_pesimista": fill_pes,
            "espera_optimista": (fill_opt - t0) if fill_opt else None,
            "espera_pesimista": (fill_pes - t0) if fill_pes else None,
            "k_fin": int(k)}


def mid_en(ev: dict, t_obj: float) -> float:
    """Mid vigente en `t_obj`: ultimo snapshot ESTRICTAMENTE anterior."""
    i = int(np.searchsorted(ev["bk_t"], t_obj, side="right") - 1)
    if i < 0:
        return float("nan")
    i = min(i, len(ev["bk_t"]) - 1)
    return 0.5 * (ev["bk_b"][i] + ev["bk_a"][i])


def markout(ev: dict, t_fill: float, p0: float, horizontes=(1.0, 5.0, 50.0)) -> dict:
    """Markout de una COMPRA pasiva llenada en `t_fill` al precio `p0`.

    Positivo = el precio subio despues de comprar (favorable). Negativo = te
    llenaron justo antes de que el mercado se fuera en contra, que es la firma
    de la seleccion adversa y lo que se paga por estar al frente de la cola.
    """
    out = {}
    for h in horizontes:
        m = mid_en(ev, t_fill + h)
        out["mk_%gs" % h] = float(m - p0) if np.isfinite(m) else float("nan")
    return out


def descomponer_agotamiento(ev: dict, n_max: int = 2000000) -> dict:
    """Flujos de la cola del bid, separando lo IDENTIFICABLE de lo que no lo es.

    ⚠ CORRECCION DE UNA VERSION ANTERIOR DE ESTA MISMA FUNCION. Atribuia a
    transaccion `min(caida_de_B, trades_pendientes)`, o sea capaba la atribucion
    por la bajada NETA observada entre snapshots. Como el nivel se repone
    constantemente, la bajada neta es mucho menor que el volumen transado y las
    transacciones se perdian: daba 54.3 BTC contra los 2058.1 realmente transados
    al mejor bid en `captura_v32`, **un factor 38**, y de ahi salia un "99.1 % por
    cancelacion" que era artefacto entero.

    ⚠ Y el fondo del problema no es de codigo: con L1 **cancelaciones y
    reposiciones NO son separables**, solo su neto. Sobre un tramo a precio
    constante,

        B_final = B_inicial - transado - cancelado + anadido

    de donde solo se despeja `(cancelado - anadido)`. Llamar "cancelacion" a esa
    diferencia es confundir un flujo neto con uno bruto. Aqui se reportan las
    dos cantidades que SI se identifican y se nombra la que no.
    """
    t, tipo, idx = ev["t"], ev["tipo"], ev["idx"]
    n = min(len(t), n_max)
    v_trans = 0.0            # transado al mejor bid: EXACTO, viene de @trade
    bajadas = 0.0            # suma de bajadas de B a precio constante (bruto)
    subidas = 0.0            # suma de subidas de B a precio constante (bruto)
    n_barridos = n_bajadas = 0
    B_prev = p_prev = None

    for k in range(n):
        if tipo[k] == 1:
            j = int(idx[k])
            if ev["tr_maker"][j] and p_prev is not None \
               and abs(ev["tr_p"][j] - p_prev) <= TOL:
                v_trans += float(ev["tr_q"][j])      # SIN capar
            continue
        i = int(idx[k])
        p_now = float(ev["bk_b"][i]); B_now = float(ev["bk_B"][i])
        if p_prev is not None and abs(p_now - p_prev) <= TOL:
            dif = B_now - B_prev
            if dif < 0:
                bajadas += -dif
                n_bajadas += 1
            elif dif > 0:
                subidas += dif
        elif p_prev is not None:
            n_barridos += 1
        B_prev, p_prev = B_now, p_now

    neto_retirado = bajadas - subidas       # (cancelado - anadido) + transado
    no_explicado = neto_retirado - v_trans  # (cancelado - anadido), NO "cancelado"
    return {"transado_al_bid": v_trans,
            "bajadas_brutas": bajadas, "subidas_brutas": subidas,
            "retirada_neta": neto_retirado,
            "no_explicado_por_trades": no_explicado,
            "frac_transaccion_sobre_bajadas": (v_trans / bajadas) if bajadas > 0 else float("nan"),
            "rotacion": (bajadas / max(v_trans, 1e-9)),
            "n_bajadas": n_bajadas, "n_cambios_de_precio": n_barridos,
            "eventos_recorridos": n}


def estudio(d: dict, n_inserciones: int = 1500, horizonte_s: float = 60.0,
            semilla: int = 0) -> dict:
    """Las cuatro mediciones sobre una captura."""
    ev = construir_eventos(d)
    rng = np.random.default_rng(semilla)
    # Solo se insertan en eventos de LIBRO, y con margen para el horizonte.
    candidatos = np.flatnonzero(ev["tipo"] == 0)
    t_lim = ev["t"][-1] - horizonte_s - 1.0
    candidatos = candidatos[ev["t"][candidatos] < t_lim]
    if candidatos.size == 0:
        return {}
    sel = rng.choice(candidatos, size=int(min(n_inserciones, candidatos.size)),
                     replace=False)
    sel.sort()

    res, marks = [], []
    for k0 in sel:
        r = simular_insercion(ev, int(k0), horizonte_s)
        if r.get("desenlace") == "INVALIDO":
            continue
        res.append(r)
        if r["desenlace"] == "LLENADO" and r["t_fill_pesimista"]:
            marks.append(markout(ev, r["t_fill_pesimista"], r["p0"]))
    return {"ev": ev, "resultados": res, "markouts": marks,
            "agotamiento": descomponer_agotamiento(ev)}


def _fmt_pct(x):
    return "n/a" if not np.isfinite(x) else "%.1f %%" % (100 * x)


def informe(est: dict, horizonte_s: float) -> None:
    res = est["resultados"]
    n = len(res)
    if n == 0:
        print("sin inserciones simuladas")
        return

    print("-" * 74)
    print("1. SUPERVIVENCIA DE LA COLA (%d inserciones, horizonte %.0f s)"
          % (n, horizonte_s))
    print("-" * 74)
    desen = {}
    for r in res:
        desen[r["desenlace"]] = desen.get(r["desenlace"], 0) + 1
    for k in ("LLENADO", "NIVEL_BARRIDO", "SUPERADO", "CENSURADO"):
        c = desen.get(k, 0)
        print("   %-14s %5d  (%s)" % (k, c, _fmt_pct(c / n)))

    llenos = [r for r in res if r["desenlace"] == "LLENADO"]
    if llenos:
        esp = np.array([r["espera_pesimista"] for r in llenos], float)
        print("   espera hasta el llenado [s]: p25=%.2f  MED=%.2f  p75=%.2f  max=%.2f"
              % (*np.percentile(esp, [25, 50, 75]), esp.max()))

    # P(sin llenar tras tau) estratificada por tamano de cola por delante.
    Q = np.array([r["Q0"] for r in res], float)
    cortes = np.percentile(Q, [33.3, 66.7])
    etiquetas = ("cola PEQUENA", "cola MEDIA", "cola GRANDE")
    print("")
    print("   P(SIN llenar tras tau | Q), por tercil de cola por delante:")
    print("   %-14s %8s %9s %9s %9s %9s"
          % ("estrato", "Q med", "tau=1s", "tau=5s", "tau=15s", "tau=60s"))
    for g, etq in enumerate(etiquetas):
        if g == 0:
            m = Q <= cortes[0]
        elif g == 1:
            m = (Q > cortes[0]) & (Q <= cortes[1])
        else:
            m = Q > cortes[1]
        sub = [r for r, mm in zip(res, m) if mm]
        if not sub:
            continue
        fila = []
        for tau in (1.0, 5.0, 15.0, 60.0):
            sin_llenar = sum(
                1 for r in sub
                if not (r["desenlace"] == "LLENADO"
                        and r["espera_pesimista"] is not None
                        and r["espera_pesimista"] <= tau))
            fila.append(sin_llenar / len(sub))
        print("   %-14s %8.3f %9s %9s %9s %9s"
              % (etq, float(np.median([r["Q0"] for r in sub])),
                 *[_fmt_pct(x) for x in fila]))

    print("")
    print("-" * 74)
    print("2. MARKOUT SINTETICO de los llenados (compra pasiva, USD/BTC)")
    print("-" * 74)
    mk = est["markouts"]
    if not mk:
        print("   sin llenados con markout evaluable")
    else:
        print("   %-8s %9s %9s %9s %9s" % ("h", "media", "MEDIANA", "p25", "p75"))
        for h in (1.0, 5.0, 50.0):
            v = np.array([m["mk_%gs" % h] for m in mk], float)
            v = v[np.isfinite(v)]
            if v.size < 3:
                continue
            print("   %-8s %+9.4f %+9.4f %+9.4f %+9.4f"
                  % ("%gs" % h, v.mean(), np.median(v), *np.percentile(v, [25, 75])))
        print("   AVISO: la mediana de +%.2f a horizonte corto NO es beneficio."
              % (TICK / 2))
        print("          Es la MEDIA HORQUILLA: compraste al bid y el mid esta")
        print("          medio tick por encima por definicion. Lo informativo es")
        print("          que la MEDIA sea muy negativa -- cola de seleccion")
        print("          adversa -- y que a 50 s la mediana ya lo sea tambien.")

    print("")
    print("-" * 74)
    print("3. DESCOMPOSICION DEL AGOTAMIENTO DE LA COLA DEL BID")
    print("-" * 74)
    ag = est["agotamiento"]
    print("   IDENTIFICABLE:")
    print("     transado al mejor bid (de @trade) : %10.1f BTC" % ag["transado_al_bid"])
    print("     bajadas brutas de B               : %10.1f BTC" % ag["bajadas_brutas"])
    print("     subidas brutas de B               : %10.1f BTC" % ag["subidas_brutas"])
    print("     retirada NETA                     : %10.1f BTC" % ag["retirada_neta"])
    print("   NO IDENTIFICABLE con L1:")
    print("     (cancelado - anadido)             : %10.1f BTC" % ag["no_explicado_por_trades"])
    print("     no se puede separar cancelacion de reposicion: solo su neto.")
    print("   lecturas:")
    print("     transado / bajadas brutas         : %s"
          % _fmt_pct(ag["frac_transaccion_sobre_bajadas"]))
    print("     rotacion de cola (bajadas/transado): %10.1f x" % ag["rotacion"])
    print("     bajadas de B a precio constante   : %d" % ag["n_bajadas"])
    print("     cambios de precio del bid         : %d" % ag["n_cambios_de_precio"])

    print("")
    print("-" * 74)
    print("4. HORQUILLADO DE LA POSICION EN COLA (cuanto cuesta NO tener @depth)")
    print("-" * 74)
    con_ambas = [r for r in res
                 if r["espera_optimista"] is not None
                 and r["espera_pesimista"] is not None]
    solo_opt = [r for r in res
                if r["espera_optimista"] is not None
                and r["espera_pesimista"] is None]
    print("   llenado bajo AMBAS atribuciones      : %d" % len(con_ambas))
    print("   llenado SOLO si las cancelaciones van")
    print("   por delante (optimista)              : %d" % len(solo_opt))
    if con_ambas:
        o = np.array([r["espera_optimista"] for r in con_ambas], float)
        p = np.array([r["espera_pesimista"] for r in con_ambas], float)
        anchura = p - o
        print("   espera optimista  [s]: MED %.3f" % np.median(o))
        print("   espera pesimista  [s]: MED %.3f" % np.median(p))
        print("   ANCHURA del horquillado [s]: MED %.3f  p90 %.3f  max %.3f"
              % (np.median(anchura), np.percentile(anchura, 90), anchura.max()))
    n_amb = len(con_ambas) + len(solo_opt)
    if n_amb:
        print("   fraccion de llenados que DEPENDEN de la atribucion: %s"
              % _fmt_pct(len(solo_opt) / n_amb))
        print("   -> esa fraccion es la incertidumbre que compraria `@depth`.")


# ===========================================================================
# Autotest -- verdad conocida
# ===========================================================================

def _autotest() -> int:
    fallos = 0

    def ok(nombre, cond, detalle=""):
        nonlocal fallos
        if not cond:
            fallos += 1
        print("  [%s] %s %s" % ("OK  " if cond else "FALLA", nombre, detalle))

    print("== 1. CONTROL POSITIVO: cola consumida solo por transacciones ==")
    # Bid a 100.0 con 3 BTC por delante. Llegan trades vendedores de 1 BTC cada
    # segundo: el llenado debe caer al tercer trade, o sea a los ~3 s.
    bk_t = np.arange(0, 11, 1.0)
    bk_b = np.full(11, 100.0)
    bk_B = np.array([3.0, 3.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    bk_a = np.full(11, 100.1)
    bk_A = np.full(11, 5.0)
    tr_t = np.array([1.5, 2.5, 3.5])
    tr_p = np.full(3, 100.0)
    tr_q = np.full(3, 1.0)
    tr_m = np.ones(3, np.uint8)          # m=True -> agresor vendedor -> pega al bid
    d = {"bk_t": bk_t, "bk_b": bk_b, "bk_B": bk_B, "bk_a": bk_a, "bk_A": bk_A,
         "tr_t": tr_t, "tr_precio": tr_p, "tr_cant": tr_q, "tr_maker": tr_m}
    ev = construir_eventos(d)
    k0 = int(np.flatnonzero((ev["tipo"] == 0) & (ev["t"] == 0.0))[0])
    r = simular_insercion(ev, k0, horizonte_s=10.0)
    print("     desenlace=%s  espera_pesimista=%s" % (r["desenlace"], r["espera_pesimista"]))
    ok("llena", r["desenlace"] == "LLENADO")
    ok("a los ~3.5 s (el tercer trade)", abs(r["espera_pesimista"] - 3.5) < 0.6)

    print("== 2. CONTROL NEGATIVO: la cola se CANCELA, no se transa ==")
    # Misma bajada de B pero SIN transacciones. Bajo la atribucion pesimista
    # (cancelaciones detras) NO deberia llenarse; bajo la optimista si.
    d2 = dict(d)
    d2["tr_t"] = np.array([]); d2["tr_precio"] = np.array([])
    d2["tr_cant"] = np.array([]); d2["tr_maker"] = np.array([], np.uint8)
    ev2 = construir_eventos(d2)
    k0 = int(np.flatnonzero((ev2["tipo"] == 0) & (ev2["t"] == 0.0))[0])
    r2 = simular_insercion(ev2, k0, horizonte_s=10.0)
    print("     desenlace=%s | optimista=%s | pesimista=%s"
          % (r2["desenlace"], r2["espera_optimista"], r2["espera_pesimista"]))
    ok("el pesimista NO llena", r2["espera_pesimista"] is None)
    ok("el optimista SI llena", r2["espera_optimista"] is not None)
    print("     (ese es el horquillado: mismo dato, dos respuestas)")

    print("== 3. El precio se va y el nivel queda barrido ==")
    bk_b3 = np.full(11, 100.0); bk_b3[4:] = 99.9
    d3 = dict(d); d3["bk_b"] = bk_b3
    d3["bk_B"] = np.full(11, 3.0)
    d3["tr_t"] = np.array([]); d3["tr_precio"] = np.array([])
    d3["tr_cant"] = np.array([]); d3["tr_maker"] = np.array([], np.uint8)
    ev3 = construir_eventos(d3)
    k0 = int(np.flatnonzero((ev3["tipo"] == 0) & (ev3["t"] == 0.0))[0])
    r3 = simular_insercion(ev3, k0, horizonte_s=10.0)
    print("     desenlace=%s a los %.1f s" % (r3["desenlace"], r3["duracion"]))
    ok("detecta el barrido", r3["desenlace"] == "NIVEL_BARRIDO")

    print("== 4. El markout mide el signo correcto ==")
    # Compra a 100.0 y el mid sube a 100.5: markout POSITIVO.
    bk_t4 = np.arange(0, 61, 1.0)
    d4 = {"bk_t": bk_t4, "bk_b": np.where(bk_t4 < 10, 100.0, 100.4),
          "bk_B": np.full(61, 3.0), "bk_a": np.where(bk_t4 < 10, 100.1, 100.6),
          "bk_A": np.full(61, 3.0), "tr_t": np.array([]),
          "tr_precio": np.array([]), "tr_cant": np.array([]),
          "tr_maker": np.array([], np.uint8)}
    ev4 = construir_eventos(d4)
    mk = markout(ev4, t_fill=1.0, p0=100.0)
    print("     markout 1s=%+.3f  5s=%+.3f  50s=%+.3f"
          % (mk["mk_1s"], mk["mk_5s"], mk["mk_50s"]))
    ok("positivo a 50 s tras subir el mid", mk["mk_50s"] > 0.3)
    ok("cerca de cero a 1 s (aun no se movio)", abs(mk["mk_1s"]) < 0.1)

    print("")
    print("RESULTADO: %d fallo(s)" % fallos)
    return fallos


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="telemetria/captura_v32")
    ap.add_argument("--inserciones", type=int, default=1500)
    ap.add_argument("--horizonte", type=float, default=60.0)
    ap.add_argument("--semilla", type=int, default=0)
    ap.add_argument("--autotest", action="store_true")
    a = ap.parse_args()
    if a.autotest:
        return 1 if _autotest() else 0

    d = cl.cargar_larga(a.dir)
    if "bk_B" not in d:
        print("La captura %s no tiene cantidades de libro: nada que medir aqui."
              % a.dir)
        return 1
    print("=" * 74)
    print("EJECUCION PASIVA -- que se puede saber sin poner una orden")
    print("=" * 74)
    print("captura: %s | %d trades | %d snapshots de libro"
          % (a.dir, len(d["tr_t"]), len(d["bk_t"])))
    print("AVISO: esto NO es la v3.2 y no toca su decision ni su conjunto de prueba.")
    print("")
    est = estudio(d, a.inserciones, a.horizonte, a.semilla)
    if not est:
        print("no hay eventos suficientes")
        return 1
    informe(est, a.horizonte)
    return 0


if __name__ == "__main__":
    sys.exit(main())
