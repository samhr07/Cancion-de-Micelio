"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: captura_dual.py — Captura simultánea de `@trade` y `@bookTicker`

    python captura_dual.py --segundos=900

Sirve a tres criterios de aceptación de la v2.1 a la vez:
  §3.4  fracción de duplicados de `mid` contra precio de transacción, a igual K
  §3.3  fracción de cambios de `mid` SIN transacción intermedia
  §2.1  ventanas suficientes para que el barrido de la banda discrimine

⚠ POR QUE UNA SOLA CAPTURA Y NO DOS. Las dos series tienen que ser del MISMO
tramo de mercado o la comparación no significa nada: si `mid` se captura en un
régimen y `trade` en otro, cualquier diferencia de duplicados podría ser del
régimen y no del observable. Es el mismo principio que hizo que el §5.2 de la
v2.0 construyera sus tres vías desde el mismo flujo.

El reloj lo sigue marcando `@trade` (§3.2): `mid` es un ESTADO que se lee, no un
evento que se cuenta, así que no introduce un segundo reloj.
"""

from __future__ import annotations

import asyncio
import json
import os
import sys
import time

import numpy as np

ARCHIVO = os.path.join("telemetria", "captura_dual.npz")


async def capturar(segundos: float, symbol: str = "btcusdt"):
    import aiohttp

    url = (
        f"wss://fstream.binance.com/stream?streams="
        f"{symbol}@trade/{symbol}@bookTicker"
    )
    # Transacciones
    tr_precio, tr_t, tr_id = [], [], []
    # Libro: se guarda el mid Y el instante, para poder reconstruir el mid
    # vigente en cada transaccion.
    bk_mid, bk_t = [], []
    t0 = time.time()
    async with aiohttp.ClientSession() as s:
        async with s.ws_connect(url, heartbeat=30.0) as ws:
            while time.time() - t0 < segundos:
                try:
                    m = await asyncio.wait_for(ws.receive(), 20.0)
                except asyncio.TimeoutError:
                    break
                if m.type is not aiohttp.WSMsgType.TEXT:
                    continue
                env = json.loads(m.data)
                d = env.get("data", env)
                e = d.get("e")
                if e == "trade":
                    tr_precio.append(float(d["p"]))
                    tr_t.append(float(d["T"]) / 1000.0)
                    tr_id.append(int(d["t"]))
                elif e == "bookTicker":
                    bid, ask = float(d["b"]), float(d["a"])
                    if bid > 0.0 and ask > 0.0:
                        bk_mid.append(0.5 * (bid + ask))
                        bk_t.append(float(d.get("T", d.get("E", 0))) / 1000.0)
                if (len(tr_precio) + len(bk_mid)) % 20000 == 0 and tr_precio:
                    print(
                        f"    ... {len(tr_precio)} trades, {len(bk_mid)} book en "
                        f"{time.time()-t0:.0f} s",
                        flush=True,
                    )
    os.makedirs(os.path.dirname(ARCHIVO), exist_ok=True)
    np.savez(
        ARCHIVO,
        tr_precio=np.array(tr_precio),
        tr_t=np.array(tr_t),
        tr_id=np.array(tr_id, dtype=np.int64),
        bk_mid=np.array(bk_mid),
        bk_t=np.array(bk_t),
    )
    return ARCHIVO


def main(argv):
    segundos = 900.0
    for a in argv[1:]:
        if a.startswith("--segundos="):
            segundos = float(a.split("=", 1)[1])
    print(f"Capturando {segundos:.0f} s de @trade + @bookTicker...")
    ruta = asyncio.run(capturar(segundos))
    d = np.load(ruta)
    dur = d["tr_t"][-1] - d["tr_t"][0]
    print()
    print(f"Guardado en {ruta}")
    print(f"  transacciones : {len(d['tr_precio']):6d}  ({len(d['tr_precio'])/dur:6.1f} tx/s)")
    print(f"  book updates  : {len(d['bk_mid']):6d}  ({len(d['bk_mid'])/dur:6.1f} msg/s)")
    print(f"  duracion      : {dur:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
