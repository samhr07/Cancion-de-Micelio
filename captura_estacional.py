# -*- coding: utf-8 -*-
"""
captura_estacional.py -- captura de 14-21 dias para el §1 de la v4.0.

    python captura_estacional.py --dias=21
    python captura_estacional.py --cobertura        resumen sin capturar

⚠ ES LO UNICO DE LA v4.0 QUE EL CALENDARIO NO PERDONA. La estacionalidad exige
semanas de reloj de pared y no se acelera con computo ni con mejor metodo. Todo
lo demas del documento se ejecuta sobre datos que se van acumulando.

⚠ NO TOCA `captura_v33`. Es un proceso independiente con su propia conexion. La
decision de la v3.2 tiene PRIORIDAD y su compuerta es de continuidad, asi que si
las dos capturas se estorban **la que se para es esta**. Implementado: cada
volcado comprueba la frescura de `captura_v33` y esta captura **se detiene sola**
si aquella lleva mas de `UMBRAL_V33` sin datos (§8, fallo 1).

LA COMPUERTA AQUI ES DISTINTA
-----------------------------
La v3.2 exige un TRAMO CONTINUO porque su estimador necesita bloques contiguos.
La estacionalidad necesita otra cosa: **COBERTURA**. Un corte de 40 min no la
invalida; que falte sistematicamente la franja 03:00-05:00 UTC si.

    compuerta = min sobre las 168 casillas (hora UTC x dia de la semana) de
                los minutos observados en esa casilla  >=  30

⚠ Con 14 dias son 2 observaciones por casilla hora x dia y 14 por casilla de
hora. Eso basta para el perfil horario marginal y **NO** basta para la
interaccion hora x dia: el perfil semanal es EXPLORATORIO hasta >= 4 semanas.

TODO EN UTC
-----------
Persistencia y analisis. El horario de verano es un fallo silencioso clasico
(§8, fallo 5): los efectos se desplazan una hora a mitad de registro y no se ve.
"""

from __future__ import annotations

import argparse
import asyncio
import datetime as dt
import json
import os
import shutil
import sys
import time

import numpy as np

import captura_larga as cl


DIR_POR_DEFECTO = os.path.join("telemetria", "estacional")
PERIODO_VOLCADO = 600.0     # [s] 10 min; con rotacion DIARIA de fichero
UMBRAL_V33 = 1800.0         # [s] sin datos en captura_v33 -> esta se para
DIR_V33 = os.path.join("telemetria", "captura_v33")


def dia_utc(t: float) -> str:
    return dt.datetime.fromtimestamp(t, dt.timezone.utc).strftime("%Y%m%d")


class BufferDia:
    """Acumula en memoria y vuelca a Parquet, con rotacion DIARIA en UTC."""

    def __init__(self, directorio: str):
        self.dir = directorio
        os.makedirs(directorio, exist_ok=True)
        self.limpiar()

    def limpiar(self):
        self.tr = {k: [] for k in ("t", "precio", "cant", "id", "maker")}
        self.bk = {k: [] for k in ("t", "b", "B", "a", "A", "u")}

    def n(self) -> int:
        return len(self.tr["t"]) + len(self.bk["t"])

    def volcar(self):
        """Un fichero por dia UTC y tipo. Se ANEXA por particiones numeradas:
        reescribir un Parquet de un dia entero cada 10 min seria O(n^2)."""
        import pyarrow as pa
        import pyarrow.parquet as pq

        escritos = []
        for nombre, col in (("trades", self.tr), ("libro", self.bk)):
            if not col["t"]:
                continue
            d = dia_utc(col["t"][0])
            sub = os.path.join(self.dir, "%s_%s" % (nombre, d))
            os.makedirs(sub, exist_ok=True)
            k = len([f for f in os.listdir(sub) if f.endswith(".parquet")])
            tabla = pa.table({c: pa.array(v) for c, v in col.items()})
            ruta = os.path.join(sub, "parte_%05d.parquet" % k)
            tmp = ruta + ".tmp"
            # Escritura a temporal y rename, como en captura_larga: un corte a
            # media escritura deja intacto lo anterior.
            pq.write_table(tabla, tmp, compression="zstd")
            os.replace(tmp, ruta)
            escritos.append((ruta, len(col["t"]), os.path.getsize(ruta)))
        self.limpiar()
        return escritos


def v33_esta_viva(dir_v33: str = DIR_V33, umbral: float = UMBRAL_V33) -> tuple:
    """Frescura del ultimo dato de `captura_v33`. La v3.2 tiene prioridad."""
    try:
        import glob
        rutas = sorted(glob.glob(os.path.join(dir_v33, "bloque_*.npz")))
        if not rutas:
            return True, float("nan")     # aun no arranco; no se penaliza
        d = np.load(rutas[-1])
        edad = time.time() - float(d["tr_t"][-1])
        return edad <= umbral, edad
    except Exception:
        return True, float("nan")


def matriz_cobertura(directorio: str) -> dict:
    """Minutos observados en cada casilla (hora UTC x dia de la semana).

    Se cuenta un minuto como cubierto si hay al menos una transaccion en el.
    Es la compuerta del §1.2, y se reporta ENTERA -- no un porcentaje agregado.
    """
    import glob
    import pyarrow.parquet as pq

    minutos = set()
    for sub in sorted(glob.glob(os.path.join(directorio, "trades_*"))):
        for f in sorted(glob.glob(os.path.join(sub, "*.parquet"))):
            try:
                t = pq.read_table(f, columns=["t"]).column("t").to_numpy()
            except Exception:
                continue
            minutos.update(np.unique((t // 60).astype(np.int64)).tolist())
    if not minutos:
        return {}

    M = np.zeros((7, 24), dtype=np.int64)   # dia de la semana x hora
    for m in minutos:
        u = dt.datetime.fromtimestamp(m * 60, dt.timezone.utc)
        M[u.weekday(), u.hour] += 1
    return {"matriz": M, "n_minutos": len(minutos),
            "minimo_casilla": int(M.min()), "casillas_vacias": int((M == 0).sum())}


def informe_cobertura(directorio: str) -> int:
    c = matriz_cobertura(directorio)
    if not c:
        print("aun no hay datos en %s" % directorio)
        return 1
    M = c["matriz"]
    dias = ("lun", "mar", "mie", "jue", "vie", "sab", "dom")
    print("=" * 78)
    print("COBERTURA -- minutos observados por casilla (hora UTC x dia)")
    print("=" * 78)
    print("     " + " ".join("%3d" % h for h in range(24)))
    for i, nom in enumerate(dias):
        print("%-4s " % nom + " ".join("%3d" % v for v in M[i]))
    print("")
    print("minutos totales    : %d (%.2f dias equivalentes)"
          % (c["n_minutos"], c["n_minutos"] / 1440.0))
    print("minimo por casilla : %d   (la compuerta del Sec. 1.2 pide >= 30)"
          % c["minimo_casilla"])
    print("casillas vacias    : %d de 168" % c["casillas_vacias"])
    print("")
    # Perfil horario marginal: es lo unico no exploratorio con 14 dias.
    porh = M.sum(axis=0)
    print("perfil HORARIO marginal (suma sobre dias):")
    print("     " + " ".join("%4d" % h for h in range(24)))
    print("min  " + " ".join("%4d" % v for v in porh))
    print("")
    print("AVISO: con 14 dias hay 2 observaciones por casilla hora x dia. El")
    print("       perfil semanal es EXPLORATORIO hasta >= 4 semanas; el horario")
    print("       marginal si es utilizable.")
    ok = c["minimo_casilla"] >= 30
    print("")
    print("COMPUERTA DE COBERTURA: %s" % ("PASA" if ok else "NO PASA -- seguir capturando"))
    return 0 if ok else 2


async def capturar(dias: float, directorio: str, symbol: str = "btcusdt"):
    import aiohttp

    url = ("wss://fstream.binance.com/stream?streams="
           f"{symbol}@trade/{symbol}@bookTicker")
    buf = BufferDia(directorio)
    t0 = time.time()
    fin = t0 + dias * 86400.0
    t_ultimo = t0
    t_ultimo_resumen = t0
    libro_prev = (0.0, 0.0, 0.0, 0.0)
    n_tr = n_inval = n_desconexiones = 0
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
                            if m.type in (aiohttp.WSMsgType.CLOSED,
                                          aiohttp.WSMsgType.ERROR):
                                raise ConnectionError("WS %s" % m.type.name)
                            continue
                        env = json.loads(m.data)
                        d = env.get("data", env)
                        e = d.get("e")
                        if e == "trade":
                            p, q = float(d["p"]), float(d["q"])
                            if p <= 0.0 or q <= 0.0:
                                n_inval += 1
                                continue
                            buf.tr["t"].append(float(d["T"]) / 1000.0)
                            buf.tr["precio"].append(p)
                            buf.tr["cant"].append(q)
                            buf.tr["id"].append(int(d["t"]))
                            buf.tr["maker"].append(bool(d.get("m")))
                            n_tr += 1
                        elif e == "bookTicker":
                            b, a = float(d["b"]), float(d["a"])
                            qb, qa = float(d["B"]), float(d["A"])
                            if b > 0 and a > 0 and qb > 0 and qa > 0 \
                               and (b, qb, a, qa) != libro_prev:
                                buf.bk["t"].append(
                                    float(d.get("T", d.get("E", 0))) / 1000.0)
                                buf.bk["b"].append(b); buf.bk["B"].append(qb)
                                buf.bk["a"].append(a); buf.bk["A"].append(qa)
                                buf.bk["u"].append(int(d.get("u", 0)))
                                libro_prev = (b, qb, a, qa)

                        ahora = time.time()
                        if (ahora - t_ultimo) >= PERIODO_VOLCADO:
                            esc = buf.volcar()
                            t_ultimo = ahora
                            # PRIORIDAD DE LA v3.2 (§8, fallo 1)
                            viva, edad = v33_esta_viva()
                            if not viva:
                                print("    [!!] captura_v33 lleva %.0f s sin datos."
                                      " ESTA captura se detiene: la v3.2 tiene"
                                      " prioridad." % edad, flush=True)
                                return n_tr
                            if (ahora - t_ultimo_resumen) >= 86400.0:
                                _resumen_diario(directorio, t0, n_tr,
                                                n_desconexiones, n_inval)
                                t_ultimo_resumen = ahora
                            if esc:
                                tot = sum(x[2] for x in esc)
                                print("    [%5.2f d] %d trades | %.1f tx/s | "
                                      "%.2f MB volcados | v33 hace %.0f s"
                                      % ((ahora - t0) / 86400.0, n_tr,
                                         n_tr / (ahora - t0), tot / 1e6, edad),
                                      flush=True)
        except asyncio.CancelledError:
            raise
        except Exception as err:
            buf.volcar()
            t_ultimo = time.time()
            n_desconexiones += 1
            espera = min(60.0, 0.5 * (2 ** min(reintentos, 12)))
            print("    [!] %s: %s. Reintento en %.0f s"
                  % (type(err).__name__, err, espera), flush=True)
            await asyncio.sleep(espera)
            reintentos += 1

    buf.volcar()
    _resumen_diario(directorio, t0, n_tr, n_desconexiones, n_inval)
    return n_tr


def _resumen_diario(directorio, t0, n_tr, n_desc, n_inval):
    """§1.3: si el proceso muere, la perdida debe detectarse en horas."""
    libre = shutil.disk_usage(".")[2] / 1e9
    horas = (time.time() - t0) / 3600.0
    c = matriz_cobertura(directorio)
    print("=" * 70, flush=True)
    print("RESUMEN DIARIO  %s UTC"
          % dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M"), flush=True)
    print("  horas capturadas   : %.1f" % horas, flush=True)
    print("  transacciones      : %d (%.1f tx/s)"
          % (n_tr, n_tr / max(horas * 3600, 1)), flush=True)
    print("  desconexiones      : %d" % n_desc, flush=True)
    print("  ticks invalidos    : %d" % n_inval, flush=True)
    if c:
        print("  minimo por casilla : %d (compuerta >= 30)" % c["minimo_casilla"],
              flush=True)
        print("  casillas vacias    : %d de 168" % c["casillas_vacias"], flush=True)
    print("  disco libre        : %.1f GB" % libre, flush=True)
    print("=" * 70, flush=True)


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dias", type=float, default=21.0)
    ap.add_argument("--dir", default=DIR_POR_DEFECTO)
    ap.add_argument("--cobertura", action="store_true")
    a = ap.parse_args(argv[1:])

    if a.cobertura:
        return informe_cobertura(a.dir)

    libre = shutil.disk_usage(".")[2] / 1e9
    print("Captura de estacionalidad -- Sec. 1 de la v4.0")
    print("arranque UTC : %s"
          % dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S"))
    print("duracion     : %.1f dias | volcado cada %.0f s | rotacion DIARIA"
          % (a.dias, PERIODO_VOLCADO))
    print("destino      : %s (Parquet + zstd)" % a.dir)
    print("disco libre  : %.1f GB (la v4.0 estima 2-4 GB)" % libre)
    if libre < 10.0:
        print("ABORTA: menos de 10 GB libres.")
        return 1
    viva, edad = v33_esta_viva()
    print("captura_v33  : %s (ultimo dato hace %.0f s)"
          % ("viva" if viva else "SIN DATOS", edad))
    print("AVISO: la v3.2 tiene prioridad. Si captura_v33 lleva mas de %.0f s"
          % UMBRAL_V33)
    print("       sin datos, ESTA captura se detiene sola.")
    if cl.impedir_suspension():
        print("Suspension del sistema IMPEDIDA.")
    n = asyncio.run(capturar(a.dias, a.dir))
    print("\nTerminado: %d transacciones." % n)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
