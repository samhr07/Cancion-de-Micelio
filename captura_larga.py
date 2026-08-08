"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: captura_larga.py — Captura continua por bloques para el §4 de la v2.2

    python captura_larga.py --horas=48
    python captura_larga.py --horas=48 --dir=telemetria/captura48

⚠ VOLCADO POR BLOQUES, VERIFICADO CONTRA PARADA NO LIMPIA. El §4.4 lo exige
explícitamente: "el defecto operativo de la v2.0 —volcado solo al llenar, con
pérdida en parada no limpia— costó varias corridas enteras". Aquí cada bloque se
cierra a disco cada `PERIODO_BLOQUE` segundos, así que una parada abrupta pierde
como mucho ese último tramo.

Dimensionado (§4.4): a 102 tx/s, 48 h son ~17.6 M transacciones. Se guardan
float64 de precio y tiempo más int64 de id = 24 B/trade -> ~420 MB, más el libro.
Para acotarlo, del libro se guarda SOLO cuando el mid CAMBIA: medido en la v2.1,
únicamente el 1.1 % de los mensajes lo mueven, así que se conserva toda la
información de movimiento a 1/90 del volumen.

⚠ Este módulo NO toca `Micelio.py` ni el modelo. El §8 lo exige: "un experimento
que modifica lo que mide no mide nada".
"""

from __future__ import annotations

import asyncio
import json
import os
import sys
import time

import numpy as np

DIR_POR_DEFECTO = os.path.join("telemetria", "captura_larga")
PERIODO_BLOQUE = 300.0  # [s] cada 5 min se cierra un bloque a disco


class Bloque:
    """Acumulador en memoria de un bloque, con volcado atómico."""

    def __init__(self, directorio: str):
        self.dir = directorio
        os.makedirs(directorio, exist_ok=True)
        self.n = 0
        self.reiniciar()

    def reiniciar(self):
        self.tr_precio, self.tr_t, self.tr_id = [], [], []
        self.bk_mid, self.bk_t = [], []
        self.tr_cant = []
        # ⚠ `es_maker` se captura desde la v2.0 en el bot pero NO se persistia en
        # las capturas de analisis. Sin el no hay forzamiento medido y el §2 de
        # la v3.1 es inejecutable. Omision propia, corregida aqui.
        self.tr_maker = []

    def volcar(self):
        if not self.tr_precio:
            return None
        ruta = os.path.join(self.dir, f"bloque_{self.n:05d}.npz")
        tmp = ruta + ".tmp"
        # Escritura a temporal y rename: un corte a media escritura deja el
        # bloque anterior intacto en vez de un .npz truncado que rompe la carga.
        # ⚠ Se pasa un DESCRIPTOR DE ARCHIVO, no una ruta: `np.savez_compressed`
        # ANADE ".npz" al nombre si no lo lleva, asi que con la ruta ".tmp" se
        # escribia en "bloque.npz.tmp.npz" y el rename fallaba con WinError 2.
        # La captura corrio media hora sin guardar un solo bloque.
        with open(tmp, "wb") as fh:
            np.savez_compressed(
                fh,
                tr_precio=np.array(self.tr_precio, dtype=np.float64),
                tr_t=np.array(self.tr_t, dtype=np.float64),
                tr_id=np.array(self.tr_id, dtype=np.int64),
                tr_cant=np.array(self.tr_cant, dtype=np.float64),
                tr_maker=np.array(self.tr_maker, dtype=np.uint8),
                bk_mid=np.array(self.bk_mid, dtype=np.float64),
                bk_t=np.array(self.bk_t, dtype=np.float64),
            )
        os.replace(tmp, ruta)
        n_tr = len(self.tr_precio)
        self.n += 1
        self.reiniciar()
        return ruta, n_tr


async def capturar(horas: float, directorio: str, symbol: str = "btcusdt"):
    import aiohttp

    url = (
        f"wss://fstream.binance.com/stream?streams="
        f"{symbol}@trade/{symbol}@bookTicker"
    )
    bloque = Bloque(directorio)
    t0 = time.time()
    fin = t0 + horas * 3600.0
    t_ultimo_bloque = t0
    mid_prev = 0.0
    n_total = 0
    n_invalidos = 0
    reintentos = 0

    while time.time() < fin:
        try:
            async with aiohttp.ClientSession() as s:
                async with s.ws_connect(url, heartbeat=30.0) as ws:
                    reintentos = 0
                    while time.time() < fin:
                        try:
                            m = await asyncio.wait_for(ws.receive(), 20.0)
                        except asyncio.TimeoutError:
                            raise ConnectionError("20 s sin datos: socket estancado")
                        if m.type is not aiohttp.WSMsgType.TEXT:
                            if m.type in (
                                aiohttp.WSMsgType.CLOSED,
                                aiohttp.WSMsgType.ERROR,
                            ):
                                raise ConnectionError(f"WS {m.type.name}")
                            continue
                        env = json.loads(m.data)
                        d = env.get("data", env)
                        e = d.get("e")
                        if e == "trade":
                            # Filtro de validez del feed: ~0.2 % de los mensajes
                            # traen p=0 y q=0 con id valido (ver mercado.tick_valido).
                            if float(d["p"]) <= 0.0 or float(d["q"]) <= 0.0:
                                n_invalidos += 1
                                continue
                            bloque.tr_precio.append(float(d["p"]))
                            bloque.tr_t.append(float(d["T"]) / 1000.0)
                            bloque.tr_id.append(int(d["t"]))
                            bloque.tr_cant.append(float(d["q"]))
                            bloque.tr_maker.append(1 if d.get("m") else 0)
                            n_total += 1
                        elif e == "bookTicker":
                            b, a = float(d["b"]), float(d["a"])
                            if b > 0.0 and a > 0.0:
                                mid = 0.5 * (b + a)
                                # Solo se guarda cuando CAMBIA: el 98.9 % de los
                                # mensajes no mueven el mid (medido en la v2.1).
                                if mid != mid_prev:
                                    bloque.bk_mid.append(mid)
                                    bloque.bk_t.append(
                                        float(d.get("T", d.get("E", 0))) / 1000.0
                                    )
                                    mid_prev = mid

                        ahora = time.time()
                        if (ahora - t_ultimo_bloque) >= PERIODO_BLOQUE:
                            r = bloque.volcar()
                            t_ultimo_bloque = ahora
                            if r:
                                print(
                                    f"    [{(ahora-t0)/3600:6.2f} h] {r[0]} "
                                    f"({r[1]} trades) | total {n_total} "
                                    f"({n_total/(ahora-t0):.1f} tx/s)",
                                    flush=True,
                                )
        except asyncio.CancelledError:
            raise
        except Exception as err:
            # El watchdog aqui es mas simple que el de produccion a proposito:
            # esto es captura, no control. Lo unico que importa es no perder el
            # bloque acumulado y reconectar.
            bloque.volcar()
            t_ultimo_bloque = time.time()
            espera = min(60.0, 0.5 * (2**reintentos))
            print(f"    [!] {type(err).__name__}: {err}. Reintento en {espera:.0f} s",
                  flush=True)
            await asyncio.sleep(espera)
            reintentos += 1

    bloque.volcar()
    return n_total


def cargar_larga(directorio: str = DIR_POR_DEFECTO):
    """Concatena todos los bloques. Tolera un bloque final truncado."""
    import glob

    rutas = sorted(glob.glob(os.path.join(directorio, "bloque_*.npz")))
    if not rutas:
        raise FileNotFoundError(f"sin bloques en {directorio}")
    campos = ("tr_precio", "tr_t", "tr_id", "tr_cant", "tr_maker", "bk_mid", "bk_t")
    acum = {c: [] for c in campos}
    for r in rutas:
        try:
            d = np.load(r)
        except Exception:
            print(f"[AVISO] bloque ilegible, se omite: {r}")
            continue
        for c in campos:
            if c in d.files:
                acum[c].append(d[c])
    return {c: np.concatenate(v) for c, v in acum.items() if v}


def main(argv):
    horas = 48.0
    directorio = DIR_POR_DEFECTO
    for a in argv[1:]:
        if a.startswith("--horas="):
            horas = float(a.split("=", 1)[1])
        elif a.startswith("--dir="):
            directorio = a.split("=", 1)[1]
    print(f"Captura continua de {horas:.1f} h en {directorio}")
    print(f"Bloques cada {PERIODO_BLOQUE:.0f} s; una parada abrupta pierde como "
          f"mucho ese tramo.")
    n = asyncio.run(capturar(horas, directorio))
    print(f"\nTerminado: {n} transacciones.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
