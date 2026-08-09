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

⚠ CAMBIO v3.2 — SE PERSISTEN LAS CANTIDADES DEL LIBRO
------------------------------------------------------
Hasta aquí se guardaba **sólo el mid**, y sólo cuando cambiaba, justificándolo
con la medición de la v2.1 de que el 98.9 % de los mensajes no lo mueven. Para
reconstruir el mid eso es exacto; para el OFI está **exactamente del revés**: el
OFI de nivel 1 mide agotamiento y reposición de cola, así que los mensajes
informativos son precisamente los que cambian cantidades sin mover el precio.

Con `bk_mid` solo, `e_t` no es computable y el §1 entero de la v3.2 (M1, M1',
M2+L, el test de fuga del §3.2 y el diagnóstico del §3.3) es inejecutable. Es la
misma omisión que costó el §2 de la v3.1 con `tr_maker`, por segunda vez: el bot
lleva el dato, la captura no lo persistía.

Ahora se guardan los cuatro campos del mejor nivel —`b`, `B`, `a`, `A`— más el
`u` (updateId) para detectar huecos igual que `tr_id` los detecta en trades.

**Criterio de compresión, y esta vez es una identidad, no una heurística.** Con
la definición de Cont-Kukanov-Stoikov,

    e_n =  1{Pb_n >= Pb_{n-1}}*qb_n - 1{Pb_n <= Pb_{n-1}}*qb_{n-1}
         - 1{Pa_n <= Pa_{n-1}}*qa_n + 1{Pa_n >= Pa_{n-1}}*qa_{n-1}

si los cuatro campos son idénticos entre dos mensajes consecutivos, ambos
indicadores de cada lado valen 1 y los términos se cancelan: `e_n = 0`
exactamente. Guardar sólo los mensajes en que cambia **alguno** de los cuatro es
por tanto SIN PÉRDIDA para el OFI, no una aproximación.

Volumen esperado: ~674 msg/s de libro (medido en la v2.1) x 48 B/registro son
~32 MB/min. 8 h -> ~1.5 GB antes de comprimir. Se reporta el tamaño real en
disco cada bloque para que no sorprenda.

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
        self.tr_cant = []
        # ⚠ `es_maker` se captura desde la v2.0 en el bot pero NO se persistia en
        # las capturas de analisis. Sin el no hay forzamiento medido y el §2 de
        # la v3.1 es inejecutable. Omision propia, corregida aqui.
        self.tr_maker = []
        # Mejor nivel del libro: precio y CANTIDAD de ambos lados, mas el
        # updateId. `bk_mid` ya no se almacena -- es 0.5*(b+a) y duplicarlo son
        # 8 B por registro sobre un stream de ~674 msg/s. `cargar_larga` lo
        # sintetiza al leer, asi que nada aguas abajo lo nota.
        self.bk_b, self.bk_B = [], []
        self.bk_a, self.bk_A = [], []
        self.bk_u, self.bk_t = [], []

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
                bk_b=np.array(self.bk_b, dtype=np.float64),
                bk_B=np.array(self.bk_B, dtype=np.float64),
                bk_a=np.array(self.bk_a, dtype=np.float64),
                bk_A=np.array(self.bk_A, dtype=np.float64),
                bk_u=np.array(self.bk_u, dtype=np.int64),
                bk_t=np.array(self.bk_t, dtype=np.float64),
            )
        os.replace(tmp, ruta)
        n_tr = len(self.tr_precio)
        n_bk = len(self.bk_b)
        tam = os.path.getsize(ruta)
        self.n += 1
        self.reiniciar()
        return ruta, n_tr, n_bk, tam


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
    libro_prev = (0.0, 0.0, 0.0, 0.0)
    n_total = 0
    n_invalidos = 0
    n_bk_recibidos = 0
    n_bk_guardados = 0
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
                            qb, qa = float(d["B"]), float(d["A"])
                            n_bk_recibidos += 1
                            # Mismo filtro de validez que en trades: un lado a
                            # cero no es un estado del libro.
                            if b > 0.0 and a > 0.0 and qb > 0.0 and qa > 0.0:
                                # Se guarda cuando cambia ALGUNO de los cuatro.
                                # Si los cuatro son iguales, e_n = 0 exactamente
                                # (ver el encabezado): omitirlo es sin perdida.
                                if (b, qb, a, qa) != libro_prev:
                                    bloque.bk_b.append(b)
                                    bloque.bk_B.append(qb)
                                    bloque.bk_a.append(a)
                                    bloque.bk_A.append(qa)
                                    bloque.bk_u.append(int(d.get("u", 0)))
                                    bloque.bk_t.append(
                                        float(d.get("T", d.get("E", 0))) / 1000.0
                                    )
                                    libro_prev = (b, qb, a, qa)
                                    n_bk_guardados += 1

                        ahora = time.time()
                        if (ahora - t_ultimo_bloque) >= PERIODO_BLOQUE:
                            r = bloque.volcar()
                            t_ultimo_bloque = ahora
                            if r:
                                dt_h = (ahora - t0) / 3600.0
                                print(
                                    f"    [{dt_h:6.2f} h] {r[0]} "
                                    f"({r[1]} trades, {r[2]} libro, {r[3]/1e6:.1f} MB) "
                                    f"| total {n_total} ({n_total/(ahora-t0):.1f} tx/s) "
                                    f"| libro {n_bk_recibidos/(ahora-t0):.0f} msg/s, "
                                    f"guardado {n_bk_guardados/max(1,n_bk_recibidos):.1%}",
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
            # El exponente va acotado: con `reintentos` grande, 2**reintentos es
            # un entero enorme y `0.5 * (2**1024)` lanza OverflowError, o sea que
            # el propio watchdog mataria la captura tras un flapeo prolongado.
            espera = min(60.0, 0.5 * (2 ** min(reintentos, 12)))
            print(f"    [!] {type(err).__name__}: {err}. Reintento en {espera:.0f} s",
                  flush=True)
            await asyncio.sleep(espera)
            reintentos += 1

    bloque.volcar()
    return n_total


CAMPOS = ("tr_precio", "tr_t", "tr_id", "tr_cant", "tr_maker",
          "bk_b", "bk_B", "bk_a", "bk_A", "bk_u", "bk_t", "bk_mid")


def cargar_larga(directorio: str = DIR_POR_DEFECTO):
    """Concatena todos los bloques. Tolera un bloque final truncado.

    `bk_mid` se sintetiza a partir de `bk_b` y `bk_a` cuando no está en el
    fichero, para que las capturas nuevas —que ya no lo almacenan— se lean
    igual que las viejas.
    """
    import glob

    rutas = sorted(glob.glob(os.path.join(directorio, "bloque_*.npz")))
    if not rutas:
        raise FileNotFoundError(f"sin bloques en {directorio}")
    acum = {c: [] for c in CAMPOS}
    for r in rutas:
        try:
            d = np.load(r)
        except Exception:
            print(f"[AVISO] bloque ilegible, se omite: {r}")
            continue
        for c in CAMPOS:
            if c in d.files:
                acum[c].append(d[c])
            elif c == "bk_mid" and "bk_b" in d.files:
                acum[c].append(0.5 * (d["bk_b"] + d["bk_a"]))
    return {c: np.concatenate(v) for c, v in acum.items() if v}


def ofi_l1(bk_b, bk_B, bk_a, bk_A):
    """OFI de nivel 1 (Cont, Kukanov & Stoikov 2014) entre snapshots consecutivos.

        e_n =  1{Pb_n >= Pb_{n-1}}*qb_n - 1{Pb_n <= Pb_{n-1}}*qb_{n-1}
             - 1{Pa_n <= Pa_{n-1}}*qa_n + 1{Pa_n >= Pa_{n-1}}*qa_{n-1}

    Devuelve un vector de longitud len-1 alineado al SEGUNDO elemento de cada par.

    ⚠ Esto es **OFI-L1 aproximado desde `bookTicker`**, no el OFI real, y así hay
    que etiquetarlo en todo reporte (§3.3 de la v3.2): `@bookTicker` es un
    snapshot con estrangulamiento, así que las cancelaciones entre
    actualizaciones son invisibles. El OFI real exige el flujo diferencial de
    `@depth`, que hoy no se ingiere.
    """
    Pb, qb = np.asarray(bk_b, float), np.asarray(bk_B, float)
    Pa, qa = np.asarray(bk_a, float), np.asarray(bk_A, float)
    sube_b = (Pb[1:] >= Pb[:-1]).astype(float)
    baja_b = (Pb[1:] <= Pb[:-1]).astype(float)
    baja_a = (Pa[1:] <= Pa[:-1]).astype(float)
    sube_a = (Pa[1:] >= Pa[:-1]).astype(float)
    return (sube_b * qb[1:] - baja_b * qb[:-1]
            - baja_a * qa[1:] + sube_a * qa[:-1])


HUECO_MAX = 300.0        # [s] por encima de esto el tramo continuo se rompe
EMBARGO_PISO = 1950      # ticks; ver PREREGISTRO_3_2 Sec. 2.2
TICKS_COMPUERTA = 380 * EMBARGO_PISO   # 741 000; ver migracion_v32.ticks_de_compuerta


def tramos_continuos(t: np.ndarray, hueco_max: float = HUECO_MAX) -> list:
    """Parte la serie en tramos SIN huecos mayores que `hueco_max`.

    ⚠ La compuerta es sobre el tramo continuo MAS LARGO, no sobre el total.
    La captura larga de la v2.2 sufrió un corte de DNS de 36 442 s y quedó
    partida en dos; sumar los dos trozos habría dado por buena una muestra que
    no existe como serie temporal.
    """
    if t.size < 2:
        return [(0, t.size)]
    cortes = np.flatnonzero(np.diff(t) > hueco_max) + 1
    ini = np.concatenate(([0], cortes))
    fin = np.concatenate((cortes, [t.size]))
    return [(int(a), int(b)) for a, b in zip(ini, fin)]


def resumen(directorio: str) -> int:
    """Comprueba la compuerta de datos del Sec. 2.1 de la v3.2 y lo dice sin adornos."""
    d = cargar_larga(directorio)
    n_tr = len(d.get("tr_t", []))
    if n_tr == 0:
        print("sin transacciones en %s" % directorio)
        return 1
    t = d["tr_t"]
    dur = float(t[-1] - t[0])
    huecos = np.diff(t)
    hueco_max = float(huecos.max()) if huecos.size else 0.0

    tramos = tramos_continuos(t)
    largos = [(b - a) for a, b in tramos]
    k_mejor = int(np.argmax(largos))
    a_m, b_m = tramos[k_mejor]
    n_cont = b_m - a_m
    dur_cont = float(t[b_m - 1] - t[a_m])
    nu_cont = n_cont / dur_cont if dur_cont > 0 else float("nan")

    print("=" * 70)
    print("COMPUERTA DE DATOS -- Sec. 2.1 de la v3.2")
    print("=" * 70)
    print("directorio      : %s" % directorio)
    print("transacciones   : %d en %.2f h  ->  nu = %.2f tx/s"
          % (n_tr, dur / 3600.0, n_tr / dur if dur > 0 else float("nan")))
    print("hueco maximo    : %.1f s" % hueco_max)
    print("tramos continuos: %d (corte a %.0f s)" % (len(tramos), HUECO_MAX))
    print("TRAMO MAS LARGO : %d ticks en %.2f h (nu = %.2f tx/s)  <- lo que cuenta"
          % (n_cont, dur_cont / 3600.0, nu_cont))
    falta = TICKS_COMPUERTA - n_cont
    if falta > 0 and nu_cont > 0:
        print("faltan          : %d ticks = %.1f h mas a la nu actual"
              % (falta, falta / nu_cont / 3600.0))
    if len(tramos) > 1:
        print("AVISO: la serie esta PARTIDA. Solo cuenta el tramo mas largo; los")
        print("       demas no se suman porque no son la misma serie temporal.")

    tiene_maker = "tr_maker" in d
    tiene_cant = "tr_cant" in d
    tiene_libro = all(c in d for c in ("bk_b", "bk_B", "bk_a", "bk_A"))
    n_bk = len(d.get("bk_t", []))

    if tiene_libro:
        e = ofi_l1(d["bk_b"], d["bk_B"], d["bk_a"], d["bk_A"])
        print("libro           : %d snapshots (%.1f msg/s guardados)"
              % (n_bk, n_bk / dur if dur > 0 else float("nan")))
        print("OFI-L1          : %d valores | mediana |e| = %.4f | %.1f %% no nulos"
              % (e.size, float(np.median(np.abs(e))),
                 100.0 * float(np.mean(e != 0.0))))
        if "bk_u" in d and d["bk_u"].size > 1:
            du = np.diff(d["bk_u"])
            # ⚠ Un salto > 1 en `u` NO es un mensaje perdido: `u` es el
            # updateId del libro COMPLETO y avanza con cualquier cambio a
            # cualquier nivel, no solo en el mejor. Medido: el 97.7 % de los
            # pares consecutivos saltan, y eso es lo normal. Lo que si seria un
            # defecto es un `u` que retrocede (mensajes desordenados) o repetido
            # (mismo estado contado dos veces).
            print("updateId        : %d retrocesos, %d repetidos (ambos deben ser 0)"
                  % (int(np.count_nonzero(du < 0)), int(np.count_nonzero(du == 0))))
    else:
        print("libro           : %d entradas, SIN cantidades (captura antigua)" % n_bk)

    print("")
    reqs = [
        ("tramo continuo >= %d ticks" % TICKS_COMPUERTA, n_cont >= TICKS_COMPUERTA,
         "%d ticks (%.2f h)" % (n_cont, dur_cont / 3600.0)),
        ("bloques de bootstrap >= 15", (int(0.20 * n_cont) - EMBARGO_PISO)
         // (5 * EMBARGO_PISO) >= 15,
         "%d bloques" % max(0, (int(0.20 * n_cont) - EMBARGO_PISO) // (5 * EMBARGO_PISO))),
        ("tr_maker persistido", tiene_maker, "si" if tiene_maker else "NO"),
        ("q persistido", tiene_cant, "si" if tiene_cant else "NO"),
        ("bookTicker con cantidades (e_t computable)", tiene_libro,
         "si" if tiene_libro else "NO"),
    ]
    for nombre, ok, detalle in reqs:
        print("  [%s] %-42s %s" % ("OK  " if ok else "FALLA", nombre, detalle))
    todas = all(ok for _, ok, _ in reqs)
    print("")
    print("COMPUERTA: %s" % ("PASA -- se puede ejecutar el §4 en adelante" if todas
                             else "NO PASA -- capturar mas antes de decidir nada"))
    return 0 if todas else 2


def main(argv):
    horas = 48.0
    directorio = DIR_POR_DEFECTO
    solo_resumen = False
    for a in argv[1:]:
        if a.startswith("--horas="):
            horas = float(a.split("=", 1)[1])
        elif a.startswith("--dir="):
            directorio = a.split("=", 1)[1]
        elif a == "--resumen":
            solo_resumen = True
    if solo_resumen:
        return resumen(directorio)
    print(f"Captura continua de {horas:.1f} h en {directorio}")
    print(f"Bloques cada {PERIODO_BLOQUE:.0f} s; una parada abrupta pierde como "
          f"mucho ese tramo.")
    print("Libro: se guardan b, B, a, A y u en cada cambio de cualquiera de los "
          "cuatro (sin perdida para el OFI; ver encabezado).")
    n = asyncio.run(capturar(horas, directorio))
    print(f"\nTerminado: {n} transacciones.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
