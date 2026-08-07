"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: test_contaminacion_emd.py — Test §5.2 de ORDEN_TRABAJO_RELOJES_2_0

    python test_contaminacion_emd.py --capturar=420   # captura y analiza
    python test_contaminacion_emd.py                  # reanaliza la captura

QUE DECIDE ESTE TEST
--------------------
El §5.2 pide tamizar LA MISMA serie de mercado por dos vias y comparar `f_hz` y
`C`. Si difieren de forma material, queda demostrado que `ω_m` estaba
contaminado EN ORIGEN y **el resultado del A/B de la v1.3 se anula por una
segunda causa independiente** de la propagacion entre correcciones.

La pregunta de fondo: la sincronizacion con el mercado, ¿se arreglo de forma
GLOBAL o solo LOCAL? El §4 puso el FILTRO en reloj de transacciones. Si el EMD
—el otro consumidor de datos de mercado, y el que produce `ω_m` y `C`— seguia en
reloj de pared, la cadena estaria sincronizada a medias y `A_arm` recibiria una
frecuencia calculada sobre una escalera.

LAS TRES VIAS QUE SE COMPARAN
-----------------------------
Se construyen desde EL MISMO flujo de trades capturado, para que la unica
diferencia sea el muestreo:

  A. ESCALERA (v1.3)   muestreo cada 0.5 s de un precio que solo se refresca a
                       la cadencia de PAQUETE (~1 Hz). Cerca de la mitad de las
                       muestras son duplicados literales. Es lo que el bot hacia.
  B. PARED SIN HUECOS  muestreo cada 0.5 s de un precio refrescado a la cadencia
                       REAL de transacciones. Sin duplicados, pero uniforme en
                       SEGUNDOS, no en ticks.
  C. TICKS (v2.0)      una muestra cada K transacciones. Uniforme en TICKS y sin
                       duplicados a la vez.

B existe para separar las dos causas: si A y B difieren, el dano lo hacian los
DUPLICADOS; si B y C difieren, el dano lo hace el ESPACIO de muestreo. Sin B, un
resultado positivo no diria cual de las dos cosas hay que arreglar.
"""

from __future__ import annotations

import asyncio
import json
import os
import sys
import time

import numpy as np

import hht

ARCHIVO_CAPTURA = os.path.join("telemetria", "captura_trades.npz")
W_VENTANA = 384  # muestras por ventana de EMD, igual que P_W_MAX del bot
PERIODO_PARED = 0.5  # [s] el PERIODO_MUESTREO_EMD del bot
CADENCIA_PAQUETE_V13 = 1.0  # [s] el sondeo REST de la v1.3


# ==============================================================================
# CAPTURA
# ==============================================================================
async def capturar(segundos: float, symbol: str = "btcusdt"):
    """Guarda (precio, T_trade) de cada transaccion real. Sin agregar nada."""
    import aiohttp

    url = f"wss://fstream.binance.com/ws/{symbol}@trade"
    precios, tiempos = [], []
    t0 = time.time()
    async with aiohttp.ClientSession() as s:
        async with s.ws_connect(url, heartbeat=30.0) as ws:
            while time.time() - t0 < segundos:
                try:
                    m = await asyncio.wait_for(ws.receive(), 15.0)
                except asyncio.TimeoutError:
                    break
                if m.type is not aiohttp.WSMsgType.TEXT:
                    continue
                d = json.loads(m.data)
                if d.get("e") != "trade":
                    continue
                precios.append(float(d["p"]))
                tiempos.append(float(d["T"]) / 1000.0)
                if len(precios) % 2000 == 0:
                    print(
                        f"    ... {len(precios)} trades en {time.time()-t0:.0f} s",
                        flush=True,
                    )
    os.makedirs(os.path.dirname(ARCHIVO_CAPTURA), exist_ok=True)
    np.savez(ARCHIVO_CAPTURA, precios=np.array(precios), tiempos=np.array(tiempos))
    return np.array(precios), np.array(tiempos)


# ==============================================================================
# LAS TRES VIAS DE MUESTREO
# ==============================================================================
def muestrear_escalera(precios, tiempos, w, periodo_pared, cadencia_paquete):
    """Via A — la v1.3: reloj de pared sobre un precio que se refresca lento.

    Reproduce el defecto exacto del §1.1: se toma `P_spot` cada `periodo_pared`
    "estuviera o no cambiado", y `P_spot` solo se refrescaba cada
    `cadencia_paquete`. El resultado es uniforme en t y constante a trozos en S.
    """
    t_ini = tiempos[0]
    # El precio visible se congela entre paquetes.
    t_pub, precio_visible, visibles = t_ini, precios[0], []
    j = 0
    t = t_ini
    while len(visibles) < w and t <= tiempos[-1]:
        # Avanza hasta el ultimo trade con T <= t
        while j + 1 < len(tiempos) and tiempos[j + 1] <= t:
            j += 1
        if (t - t_pub) >= cadencia_paquete:
            precio_visible = precios[j]
            t_pub = t
        visibles.append(precio_visible)
        t += periodo_pared
    return np.array(visibles), periodo_pared


def muestrear_pared_sin_huecos(precios, tiempos, w, periodo_pared):
    """Via B — reloj de pared pero con el precio SIEMPRE fresco.

    Aisla el efecto de los duplicados del efecto del espacio de muestreo.
    """
    t_ini = tiempos[0]
    muestras, j, t = [], 0, t_ini
    while len(muestras) < w and t <= tiempos[-1]:
        while j + 1 < len(tiempos) and tiempos[j + 1] <= t:
            j += 1
        muestras.append(precios[j])
        t += periodo_pared
    return np.array(muestras), periodo_pared


def muestrear_ticks(precios, tiempos, w, k):
    """Via C — la v2.0: una muestra cada K transacciones.

    Malla uniforme en TICKS. `dt` equivalente en segundos se deriva del span real
    cubierto, que es lo que `hht.analizar_ventana` necesita para devolver Hz.
    """
    idx = np.arange(0, min(len(precios), w * k), k)[:w]
    muestras = precios[idx]
    if len(idx) < 2:
        return muestras, 0.0
    dt_equiv = float((tiempos[idx[-1]] - tiempos[idx[0]]) / (len(idx) - 1))
    return muestras, dt_equiv


# ==============================================================================
# ANALISIS
# ==============================================================================
def analizar(nombre, muestras, dt, span_referencia=None):
    if len(muestras) < 64 or dt <= 0.0:
        print(f"  {nombre:<22} INSUFICIENTE ({len(muestras)} muestras, dt={dt})")
        return None
    duplicados = float(np.mean(np.diff(muestras) == 0.0))
    r = hht.analizar_ventana(muestras, dt)
    span = len(muestras) * dt
    print(
        f"  {nombre:<22} n={len(muestras)} dt={dt:.3f}s span={span:6.1f}s "
        f"duplicados={duplicados:5.1%} | f_hz={r['f_hz']:.5f} "
        f"(periodo {1/r['f_hz']:7.1f}s) C={r['C']:.3f} imfs={r['n_imfs']} "
        f"valido={r['valido']}"
    )
    return {
        "nombre": nombre, "f_hz": r["f_hz"], "C": r["C"],
        "duplicados": duplicados, "span": span, "n_imfs": r["n_imfs"],
        "valido": r["valido"],
    }


def barrido_banda(precios, tiempos, k):
    """§2.1 de la v2.1: elegir `ciclos_min` y `muestras_por_ciclo_min` POR EVIDENCIA.

    El documento los propone en 3.0 y 5.0 "de la literatura de EMD, no una
    medicion de este sistema", y exige barrerlos contra esta captura.

    Criterio de eleccion, en este orden:
      1. DEBE rechazar las vias A y B (el defecto de tres sesiones) y aceptar C.
         Sin esto la guarda no sirve para lo que se creo.
      2. Entre las que cumplen 1, preferir la que ACEPTE MAS ventanas reales de
         la via C: una guarda que rechaza todo es trivialmente segura e inutil.
    """
    print()
    print("=" * 92)
    print("BARRIDO DE LA BANDA DE RESOLUBILIDAD (§2.1 de la v2.1)")
    print("=" * 92)

    # Las tres vias, con su periodo y su espaciado.
    m_a, dt_a = muestrear_escalera(
        precios, tiempos, W_VENTANA, PERIODO_PARED, CADENCIA_PAQUETE_V13
    )
    m_b, dt_b = muestrear_pared_sin_huecos(precios, tiempos, W_VENTANA, PERIODO_PARED)
    m_c, dt_c = muestrear_ticks(precios, tiempos, W_VENTANA, k)
    vias = []
    for nombre, m, dt in (("A", m_a, dt_a), ("B", m_b, dt_b), ("C", m_c, dt_c)):
        r = hht.analizar_ventana(m, dt)
        vias.append((nombre, r["periodo_s"], dt, len(m)))

    # Ventanas independientes de la via C, para medir cuantas sobreviven.
    ventanas = []
    paso = W_VENTANA * k // 2  # 50 % de solape entre ventanas consecutivas
    for ini in range(0, max(1, len(precios) - W_VENTANA * k), paso):
        idx = np.arange(ini, ini + W_VENTANA * k, k)[:W_VENTANA]
        if len(idx) < W_VENTANA or idx[-1] >= len(precios):
            continue
        dt_v = float((tiempos[idx[-1]] - tiempos[idx[0]]) / (len(idx) - 1))
        r = hht.analizar_ventana(precios[idx], dt_v)
        if r["valido"]:
            ventanas.append((r["periodo_s"], dt_v))

    print(f"  Ventanas independientes de la via C evaluadas: {len(ventanas)}")
    print()
    print(f"  {'ciclos_min':>10} {'m/ciclo_min':>12} | {'A':>6} {'B':>6} {'C':>6} | "
          f"{'ventanas C aceptadas':>22} | veredicto")
    print("  " + "-" * 88)

    mejores = []
    for ciclos_min in (2.0, 2.5, 3.0, 3.5, 4.0, 5.0):
        for mpc_min in (3.0, 4.0, 5.0, 6.0, 8.0):
            acepta = {}
            for nombre, periodo, dt, W in vias:
                acepta[nombre] = hht.omega_en_banda(periodo, dt, W, ciclos_min, mpc_min)
            n_ok = sum(
                1 for p, d in ventanas
                if hht.omega_en_banda(p, d, W_VENTANA, ciclos_min, mpc_min)
            )
            frac = n_ok / max(1, len(ventanas))
            cumple = (not acepta["A"]) and (not acepta["B"]) and acepta["C"]
            marca = "CUMPLE" if cumple else ""
            print(
                f"  {ciclos_min:10.1f} {mpc_min:12.1f} | "
                f"{'si' if acepta['A'] else 'NO':>6} {'si' if acepta['B'] else 'NO':>6} "
                f"{'si' if acepta['C'] else 'NO':>6} | {n_ok:6d}/{len(ventanas):<6d} "
                f"({frac:5.1%})      | {marca}"
            )
            if cumple:
                mejores.append((frac, ciclos_min, mpc_min))

    print()
    if not mejores:
        print("  >> NINGUNA combinacion rechaza A y B aceptando C.")
        print("     La guarda de banda no basta por si sola sobre esta captura.")
        return None

    fracs = {round(f, 6) for f, _, _ in mejores}
    print(f"  {len(mejores)} combinaciones cumplen el criterio 1 (rechazar A y B, aceptar C).")
    if len(fracs) == 1:
        # ⚠ HONESTIDAD DEL BARRIDO. Con pocas ventanas independientes el criterio
        # 2 no separa: todas empatan. Elegir "la mejor" de un empate seria elegir
        # por el orden del desempate y llamarlo evidencia — exactamente el tipo de
        # conclusion sin sustento que este proyecto persigue.
        print("  >> EL CRITERIO 2 NO DISCRIMINA: todas las que cumplen 1 empatan en")
        print(f"     {mejores[0][0]:.1%} de ventanas aceptadas. Con "
              f"{len(ventanas)} ventanas independientes no hay resolucion para")
        print("     separarlas. Se CONSERVAN los valores propuestos por el documento")
        print(f"     (ciclos_min={hht.CICLOS_MIN_VENTANA}, "
              f"muestras_por_ciclo_min={hht.MUESTRAS_POR_CICLO_MIN}), que estan")
        print("     dentro de la region valida, en vez de elegir un extremo del empate.")
        print("     Para discriminar hace falta una captura mucho mas larga.")
        elegido = (hht.CICLOS_MIN_VENTANA, hht.MUESTRAS_POR_CICLO_MIN)
        dentro = any(
            abs(c - elegido[0]) < 1e-9 and abs(m - elegido[1]) < 1e-9
            for _, c, m in mejores
        )
        if not dentro:
            print("     [!] Pero la propuesta del documento NO esta entre las que")
            print("         cumplen el criterio 1. Hay que revisarla.")
        return elegido

    mejores.sort(reverse=True)
    frac, ciclos_min, mpc_min = mejores[0]
    print(f"  >> ELECCION POR EVIDENCIA: ciclos_min={ciclos_min}, "
          f"muestras_por_ciclo_min={mpc_min}")
    print(f"     Rechaza A y B, acepta C, y deja pasar el {frac:.1%} de las ventanas")
    print(f"     reales de la via C, contra {sorted(fracs)[0]:.1%} de la peor que cumple.")
    return ciclos_min, mpc_min


def main(argv):
    segundos = 0.0
    hacer_barrido = "--barrido" in argv
    for a in argv[1:]:
        if a.startswith("--capturar="):
            segundos = float(a.split("=", 1)[1])

    if segundos > 0:
        print(f"Capturando {segundos:.0f} s de btcusdt@trade real...")
        precios, tiempos = asyncio.run(capturar(segundos))
    else:
        if not os.path.exists(ARCHIVO_CAPTURA):
            print(f"[ERROR] no hay captura en {ARCHIVO_CAPTURA}. Usa --capturar=420")
            return 1
        d = np.load(ARCHIVO_CAPTURA)
        precios, tiempos = d["precios"], d["tiempos"]

    dur = tiempos[-1] - tiempos[0]
    nu = len(precios) / dur
    k = max(1, int(round(nu * PERIODO_PARED)))
    print()
    print("=" * 92)
    print("TEST §5.2 - CONTAMINACION DEL EMD POR MUESTREO DE RELOJ DE PARED")
    print("=" * 92)
    print(
        f"Serie real: {len(precios)} transacciones en {dur:.0f} s "
        f"-> nu = {nu:.1f} tx/s   |   K = round(nu*{PERIODO_PARED}) = {k} ticks/muestra"
    )
    print(f"Ventana: W = {W_VENTANA} muestras en las tres vias.")
    print()

    a = analizar(
        "A escalera (v1.3)",
        *muestrear_escalera(precios, tiempos, W_VENTANA, PERIODO_PARED, CADENCIA_PAQUETE_V13),
    )
    b = analizar(
        "B pared sin huecos",
        *muestrear_pared_sin_huecos(precios, tiempos, W_VENTANA, PERIODO_PARED),
    )
    c = analizar("C ticks (v2.0)", *muestrear_ticks(precios, tiempos, W_VENTANA, k))

    print()
    print("-" * 92)
    print("VEREDICTO")
    print("-" * 92)
    if not (a and b and c):
        print("  Captura insuficiente para las tres vias. Repetir con mas segundos.")
        return 1

    def dif(x, y, campo):
        vx, vy = x[campo], y[campo]
        if vx <= 0 or vy <= 0:
            return float("inf")
        return abs(vx - vy) / max(vx, vy)

    d_ac_f, d_ac_C = dif(a, c, "f_hz"), abs(a["C"] - c["C"])
    d_bc_f, d_bc_C = dif(b, c, "f_hz"), abs(b["C"] - c["C"])
    d_ab_f = dif(a, b, "f_hz")

    print(f"  A vs C (v1.3 contra v2.0):  f_hz difiere {d_ac_f:6.1%}   C difiere {d_ac_C:.3f}")
    print(f"  A vs B (efecto DUPLICADOS): f_hz difiere {d_ab_f:6.1%}")
    print(f"  B vs C (efecto ESPACIO):    f_hz difiere {d_bc_f:6.1%}   C difiere {d_bc_C:.3f}")
    print()

    MATERIAL = 0.20  # 20 % en f_hz, o 0.15 en C
    contaminado = d_ac_f > MATERIAL or d_ac_C > 0.15
    if contaminado:
        print("  >> CONTAMINACION MATERIAL CONFIRMADA.")
        print("     `omega_m` y `C` de la v1.3 estaban calculados sobre una escalera, no")
        print("     sobre la serie de mercado. EL RESULTADO DEL A/B DE LA v1.3 QUEDA")
        print("     ANULADO POR UNA SEGUNDA CAUSA INDEPENDIENTE de la propagacion entre")
        print("     correcciones. Anotarlo en CLAUDE.md.")
        if d_ab_f > d_bc_f:
            print("     Causa dominante: los DUPLICADOS del muestreo por reloj de pared.")
        else:
            print("     Causa dominante: el ESPACIO de muestreo (segundos contra ticks).")
    else:
        print("  >> Sin contaminacion material en esta ventana.")
        print("     El muestreo por reloj de pared NO altero f_hz ni C de forma relevante")
        print("     sobre esta serie. El veredicto del A/B de la v1.3 sigue invalidado por")
        print("     la causa ya conocida (0.76 mediciones/s), pero NO por una segunda.")
    print()
    print("  Nota: un solo tramo de mercado no generaliza. Repetir en regimenes")
    print("  distintos antes de dar la contaminacion por descartada del todo.")

    if hacer_barrido:
        barrido_banda(precios, tiempos, k)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
