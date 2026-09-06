# -*- coding: utf-8 -*-
"""
descarga_binance.py -- volcados historicos de Binance, en el formato del proyecto.

    python descarga_binance.py --listar [--dataset=bookTicker]
    python descarga_binance.py --bajar --desde=2024-03-01 --hasta=2024-03-30 --dir=RUTA
    python descarga_binance.py --autotest

QUE RESUELVE. Hasta ahora la unica fuente de dato de este proyecto era la captura
en vivo (`captura_estacional.py`), que va a ~6.5 GB por 25 dias y solo avanza en
tiempo real. Binance publica volcados historicos gratuitos, y **uno de ellos es
exactamente el stream que se estaba capturando**.

⚠ COBERTURA MEDIDA el 2026-09-06 (paginando el listado de S3, no de memoria):

    dataset      periodo          ultimo dia    estado
    bookTicker   diario            2023-05-16 -> 2024-03-30   *** DISCONTINUADO ***
    bookTicker   mensual           2023-05    -> 2024-04      *** DISCONTINUADO ***
    trades       diario            2019-09-08 -> ayer         vivo
    aggTrades    diario            2019-12-31 -> ayer         vivo
    bookDepth    diario            2023-01-01 -> ayer         vivo

**El `bookTicker` se corto en marzo de 2024.** O sea:

  - para el periodo 2023-05-16 .. 2024-03-30 hay **320 dias** de L1 con
    cantidades, que es 13x la captura propia y sirve entero para `phi` y `Omega`;
  - para HOY no hay `bookTicker`, asi que **la captura en vivo sigue siendo
    necesaria** si se quiere `tau_0` de ahora;
  - `bookDepth` SI esta vivo, pero es otra cosa: profundidad agregada en bandas
    de +-1/2/3/4/5 % cada ~10 s, no el nivel 1. Es mas PROFUNDO que `bookTicker`
    y menos fino. Ver la nota al final.

COLUMNAS, verificadas descargando (no de memoria):

    bookTicker: update_id, best_bid_price, best_bid_qty, best_ask_price,
                best_ask_qty, transaction_time, event_time
    trades:     id, price, qty, quote_qty, time, is_buyer_maker

que corresponden **una a una** con lo que escribe `captura_estacional.BufferDia`:
`t, b, B, a, A, u` y `t, precio, cant, id, maker`. `is_buyer_maker` es el campo
`m` del que sale `eps` en `propagador.signo_transaccion`.

⚠ SE DEDUPLICA EL LIBRO POR OMISION, y es una decision con consecuencia. La
captura propia solo escribe un snapshot cuando `(b, qb, a, qa)` CAMBIA; el
volcado trae todas las actualizaciones. Sin deduplicar, `tau_upd` y el piso de
parpadeo de `flujo_omega` no serian comparables entre las dos fuentes. Con
`--sin-dedup` se conserva todo.

⚠ NO MEZCLAR EPOCAS. `phi = Q_neto*P/tau_0` combina volumen, precio y libro **del
mismo instante**; juntar precio de 2024 con flujo de 2026 no es una serie, es un
artefacto. Cada dia descargado es autoconsistente (libro y transacciones del
mismo dia y del mismo simbolo), y `flujo_omega` ya corta en rachas continuas, asi
que el peligro solo aparece si alguien apunta `--datos` a un directorio que
mezcle descargas y captura propia. **No hacerlo.**
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import io
import os
import re
import urllib.parse
import urllib.request
import zipfile

import numpy as np

import horizonte as H

log, titulo = H.log, H.titulo

# El CDN `data.binance.vision` puede estar bloqueado por politica de red; el
# bucket de S3 sirve el mismo contenido y es el que se usa para listar.
HOST = "https://s3-ap-northeast-1.amazonaws.com/data.binance.vision"
MERCADO = "data/futures/um"
DATASETS = ("bookTicker", "trades", "aggTrades", "bookDepth", "klines")


def _url(dataset, simbolo, dia, periodo="daily"):
    d = dia.strftime("%Y-%m-%d" if periodo == "daily" else "%Y-%m")
    return "%s/%s/%s/%s/%s/%s-%s-%s.zip" % (
        HOST, MERCADO, periodo, dataset, simbolo, simbolo, dataset, d)


def listar(dataset, simbolo="BTCUSDT", periodo="daily"):
    """Todas las claves del prefijo, PAGINANDO.

    ⚠ S3 devuelve 1000 claves por pagina y hay que seguir `IsTruncated`. Sin
    paginar, la primera version dio "trades hasta 2021-01-19" cuando en realidad
    llega hasta ayer: estaba leyendo la ultima clave de la PRIMERA pagina.
    """
    pref = "%s/%s/%s/%s/" % (MERCADO, periodo, dataset, simbolo)
    claves, marcador = [], ""
    while True:
        q = {"prefix": pref}
        if marcador:
            q["marker"] = marcador
        with urllib.request.urlopen(HOST + "?" + urllib.parse.urlencode(q),
                                    timeout=120) as r:
            x = r.read().decode()
        k = re.findall(r"<Key>([^<]+)</Key>", x)
        claves += [y for y in k if not y.endswith(".CHECKSUM")]
        if "<IsTruncated>true</IsTruncated>" in x and k:
            marcador = k[-1]
        else:
            return claves


def _fecha_de(clave, dataset):
    m = re.search(re.escape(dataset) + r"-(\d{4}-\d{2}(?:-\d{2})?)\.zip$", clave)
    return m.group(1) if m else None


def cobertura(simbolo="BTCUSDT"):
    out = {}
    for ds in ("bookTicker", "trades", "aggTrades", "bookDepth"):
        for per in ("daily", "monthly"):
            try:
                k = listar(ds, simbolo, per)
            except Exception as e:
                out[(ds, per)] = ("ERROR", str(e), 0)
                continue
            f = [_fecha_de(x, ds) for x in k]
            f = sorted(x for x in f if x)
            out[(ds, per)] = (f[0], f[-1], len(f)) if f else (None, None, 0)
    return out


def bajar(url, verificar=True):
    """Descarga a memoria y comprueba el CHECKSUM publicado, si lo hay."""
    with urllib.request.urlopen(url, timeout=600) as r:
        datos = r.read()
    if verificar:
        try:
            with urllib.request.urlopen(url + ".CHECKSUM", timeout=60) as r:
                esperado = r.read().decode().split()[0].lower()
            real = hashlib.sha256(datos).hexdigest()
            if real != esperado:
                raise IOError("CHECKSUM no coincide para %s" % url)
        except urllib.error.HTTPError:
            pass          # no todos los ficheros publican checksum
    return datos


# ===========================================================================
# Conversion al formato de `captura_estacional`
# ===========================================================================

def _csv_a_columnas(datos, esperadas):
    """Del ZIP a un dict de arrays. Tolera que el volcado traiga o no cabecera.

    ⚠ Se lee con `pyarrow.csv`, no con `numpy.genfromtxt`. Un dia de
    `bookTicker` son 7.4 M de filas por 7 columnas y `genfromtxt` es de Python
    puro: tardaria minutos por dia, y aqui se procesan decenas. `pyarrow.csv` es
    multihilo y ya es dependencia del proyecto porque escribe el parquet.
    """
    with zipfile.ZipFile(io.BytesIO(datos)) as z:
        nombre = [n for n in z.namelist() if n.endswith(".csv")][0]
        crudo = z.read(nombre)
    cab = crudo[:400].split(b"\n", 1)[0].decode("utf-8", "replace").strip()
    tiene_cab = cab.split(",")[0] == esperadas[0]
    from pyarrow import csv as _csv
    lee = _csv.ReadOptions(
        column_names=None if tiene_cab else list(esperadas),
        autogenerate_column_names=False)
    tb = _csv.read_csv(io.BytesIO(crudo), read_options=lee)
    return {c: tb[c].to_numpy(zero_copy_only=False) for c in tb.column_names}


COLS_BT = ("update_id", "best_bid_price", "best_bid_qty", "best_ask_price",
           "best_ask_qty", "transaction_time", "event_time")
COLS_TR = ("id", "price", "qty", "quote_qty", "time", "is_buyer_maker")


def _booleano(x):
    if x.dtype.kind == "b":
        return x.astype(bool)
    return np.char.lower(x.astype(str)) == "true"


def convertir_libro(datos, dedup=True):
    c = _csv_a_columnas(datos, COLS_BT)
    t = c["transaction_time"].astype(np.float64) / 1000.0
    b = c["best_bid_price"].astype(np.float64)
    B = c["best_bid_qty"].astype(np.float64)
    a = c["best_ask_price"].astype(np.float64)
    A = c["best_ask_qty"].astype(np.float64)
    u = c["update_id"].astype(np.int64)
    ok = (b > 0) & (a > 0) & (B > 0) & (A > 0) & np.isfinite(t)
    t, b, B, a, A, u = t[ok], b[ok], B[ok], a[ok], A[ok], u[ok]
    o = np.argsort(t, kind="stable")
    t, b, B, a, A, u = t[o], b[o], B[o], a[o], A[o], u[o]
    if dedup and t.size > 1:
        # misma regla que `captura_estacional`: solo se conserva el snapshot
        # cuando la tupla (b, qb, a, qa) cambia respecto al anterior
        cam = np.r_[True, (b[1:] != b[:-1]) | (B[1:] != B[:-1])
                    | (a[1:] != a[:-1]) | (A[1:] != A[:-1])]
        t, b, B, a, A, u = t[cam], b[cam], B[cam], a[cam], A[cam], u[cam]
    return {"t": t, "b": b, "B": B, "a": a, "A": A, "u": u}


def convertir_trades(datos):
    c = _csv_a_columnas(datos, COLS_TR)
    t = c["time"].astype(np.float64) / 1000.0
    p = c["price"].astype(np.float64)
    q = c["qty"].astype(np.float64)
    i = c["id"].astype(np.int64)
    m = _booleano(c["is_buyer_maker"])
    ok = (p > 0) & (q > 0) & np.isfinite(t)
    t, p, q, i, m = t[ok], p[ok], q[ok], i[ok], m[ok]
    o = np.argsort(t, kind="stable")
    return {"t": t[o], "precio": p[o], "cant": q[o], "id": i[o], "maker": m[o]}


def escribir_parquet(destino, sub, cols):
    import pyarrow as pa
    import pyarrow.parquet as pq
    d = os.path.join(destino, sub)
    os.makedirs(d, exist_ok=True)
    ruta = os.path.join(d, "parte_00000.parquet")
    tmp = ruta + ".tmp"
    pq.write_table(pa.table({k: pa.array(v) for k, v in cols.items()}), tmp,
                   compression="zstd")
    os.replace(tmp, ruta)       # atomico, como la captura
    return ruta, os.path.getsize(ruta)


# ===========================================================================
# Etapas
# ===========================================================================

def etapa_listar(args) -> int:
    titulo("COBERTURA DE LOS VOLCADOS DE BINANCE -- %s, futuros USD-M" % args.simbolo)
    log("  (paginando el listado de S3; no es de memoria)")
    log("")
    log("  %-12s %-9s %8s   %-12s %-12s" % ("dataset", "periodo", "ficheros",
                                            "primero", "ultimo"))
    cob = cobertura(args.simbolo)
    for (ds, per), (a_, b_, n) in cob.items():
        log("  %-12s %-9s %8d   %-12s %-12s" % (ds, per, n, a_ or "-", b_ or "-"))
    log("")
    log("  [!] `bookTicker` es el UNICO que trae las cantidades de nivel 1")
    log("      (`best_bid_qty` / `best_ask_qty`), o sea el unico que sirve para")
    log("      `tau_0` y por tanto para `phi` y `Omega`. Y esta DISCONTINUADO.")
    log("      `bookDepth` sigue vivo pero es profundidad agregada en bandas de")
    log("      porcentaje cada ~10 s: mas profundo que L1 y menos fino.")
    return 0


def _dias(desde, hasta):
    a = dt.datetime.strptime(desde, "%Y-%m-%d").date()
    b = dt.datetime.strptime(hasta, "%Y-%m-%d").date()
    if b < a:
        raise SystemExit("--hasta es anterior a --desde")
    return [a + dt.timedelta(days=k) for k in range((b - a).days + 1)]


def etapa_bajar(args) -> int:
    titulo("DESCARGA -> formato de `captura_estacional`  (%s)" % args.simbolo)
    dias = _dias(args.desde, args.hasta)
    os.makedirs(args.dir, exist_ok=True)
    log("  destino: %s" % os.path.abspath(args.dir))
    log("  dias solicitados: %d  (%s .. %s)" % (len(dias), dias[0], dias[-1]))
    log("  deduplicacion del libro: %s" % ("NO" if args.sin_dedup else "SI"))
    log("")
    log("  %-12s %12s %12s %10s %10s" % ("dia", "libro filas", "trades filas",
                                         "libro MB", "trades MB"))
    ok = fallos = 0
    for d in dias:
        sufijo = d.strftime("%Y%m%d")
        if (os.path.exists(os.path.join(args.dir, "libro_" + sufijo,
                                        "parte_00000.parquet"))
                and not args.rehacer):
            log("  %-12s  ya estaba, se omite" % d)
            ok += 1
            continue
        try:
            lb = convertir_libro(bajar(_url("bookTicker", args.simbolo, d)),
                                 dedup=not args.sin_dedup)
            tr = convertir_trades(bajar(_url("trades", args.simbolo, d)))
        except Exception as e:
            log("  %-12s  FALLA: %s" % (d, str(e)[:70]))
            fallos += 1
            continue
        if lb["t"].size == 0 or tr["t"].size == 0:
            log("  %-12s  vacio, se omite" % d)
            fallos += 1
            continue
        _, n1 = escribir_parquet(args.dir, "libro_" + sufijo, lb)
        _, n2 = escribir_parquet(args.dir, "trades_" + sufijo, tr)
        log("  %-12s %12d %12d %10.1f %10.1f"
            % (d, lb["t"].size, tr["t"].size, n1 / 1e6, n2 / 1e6))
        ok += 1
    log("")
    log("  dias listos: %d   fallidos: %d" % (ok, fallos))
    log("")
    log("  ahora:  python flujo_omega.py --etapa=serie --datos=%s" % args.dir)
    return 0 if fallos == 0 else 2


def _autotest() -> int:
    titulo("descarga_binance.py -- CONTROLES")
    n_ok = n_tot = 0

    def chk(cond, msg, det=""):
        nonlocal n_ok, n_tot
        n_tot += 1
        n_ok += bool(cond)
        log("  [%s] %s%s" % ("OK " if cond else "FALLA", msg,
                             ("   %s" % det) if det else ""))

    def _zip(txt):
        b = io.BytesIO()
        with zipfile.ZipFile(b, "w") as z:
            z.writestr("x.csv", txt)
        return b.getvalue()

    # --- 1. libro: columnas, unidades y filtrado ---------------------------
    csv_bt = ("update_id,best_bid_price,best_bid_qty,best_ask_price,"
              "best_ask_qty,transaction_time,event_time\n"
              "1,100.0,2.0,100.1,3.0,1700000000000,1700000000005\n"
              "2,100.0,2.0,100.1,3.0,1700000001000,1700000001005\n"   # duplicado
              "3,100.0,4.0,100.1,3.0,1700000002000,1700000002005\n"
              "4,0.0,4.0,100.1,3.0,1700000003000,1700000003005\n")     # invalido
    lb = convertir_libro(_zip(csv_bt))
    chk(lb["t"].size == 2 and abs(lb["t"][0] - 1700000000.0) < 1e-9
        and lb["B"][1] == 4.0,
        "1 libro: ms -> s, filtra precio <= 0 y DEDUPLICA como la captura",
        "3 validas -> %d tras dedup" % lb["t"].size)
    chk(convertir_libro(_zip(csv_bt), dedup=False)["t"].size == 3,
        "1b --sin-dedup conserva las repetidas", "3")

    # --- 2. trades: el flag de maker y la convencion de eps ----------------
    csv_tr = ("id,price,qty,quote_qty,time,is_buyer_maker\n"
              "10,100.0,0.5,50.0,1700000000000,true\n"
              "11,100.0,0.25,25.0,1700000001000,false\n"
              "12,-1.0,0.25,25.0,1700000002000,false\n")
    tr = convertir_trades(_zip(csv_tr))
    import propagador as P
    eps = P.signo_transaccion(tr["maker"])
    chk(tr["t"].size == 2 and list(tr["maker"]) == [True, False]
        and list(eps) == [-1, 1],
        "2 trades: is_buyer_maker -> maker -> eps con la convencion del proyecto",
        "eps = %s (maker=True debe dar -1)" % list(eps))

    # --- 3. sin cabecera, el orden documentado ----------------------------
    sin_cab = "1,100.0,2.0,100.1,3.0,1700000000000,1700000000005\n"
    chk(convertir_libro(_zip(sin_cab))["t"].size == 1,
        "3 tolera un volcado SIN cabecera usando el orden documentado")

    # --- 4. orden temporal garantizado ------------------------------------
    desord = ("id,price,qty,quote_qty,time,is_buyer_maker\n"
              "2,100.0,1.0,100.0,1700000005000,false\n"
              "1,100.0,1.0,100.0,1700000001000,true\n")
    t4 = convertir_trades(_zip(desord))["t"]
    chk(np.all(np.diff(t4) >= 0), "4 se ordena por tiempo aunque el volcado no lo este")

    # --- 5. el parquet queda donde `curvas_estacional` lo busca -----------
    import tempfile, shutil
    import curvas_estacional as C
    d = tempfile.mkdtemp()
    guarda = C.DIR
    try:
        escribir_parquet(d, "libro_20240330", lb)
        escribir_parquet(d, "trades_20240330", tr)
        C.DIR = d
        for f in os.listdir(d):
            if f.startswith("_indice"):
                os.remove(os.path.join(d, f))
        i1, i2 = C._indice("libro_"), C._indice("trades_")
        chk(len(i1) == 1 and len(i2) == 1 and i1[0][3] == lb["t"].size,
            "5 el parquet cae donde `curvas_estacional._indice` lo encuentra",
            "libro %d parte(s), trades %d parte(s)" % (len(i1), len(i2)))
    finally:
        C.DIR = guarda
        shutil.rmtree(d, ignore_errors=True)

    # --- 6. la fecha se extrae de la clave de S3 --------------------------
    k = "data/futures/um/daily/bookTicker/BTCUSDT/BTCUSDT-bookTicker-2024-03-30.zip"
    chk(_fecha_de(k, "bookTicker") == "2024-03-30", "6 fecha extraida de la clave")

    # --- 7. la URL se arma como el volcado real ---------------------------
    u = _url("bookTicker", "BTCUSDT", dt.date(2024, 3, 30))
    chk(u.endswith("BTCUSDT/BTCUSDT-bookTicker-2024-03-30.zip"),
        "7 URL bien formada", u.split("/data/")[-1])

    log("")
    log("  %d / %d" % (n_ok, n_tot))
    return 0 if n_ok == n_tot else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--listar", action="store_true")
    ap.add_argument("--bajar", action="store_true")
    ap.add_argument("--simbolo", default="BTCUSDT")
    ap.add_argument("--desde", help="YYYY-MM-DD")
    ap.add_argument("--hasta", help="YYYY-MM-DD")
    ap.add_argument("--dir", default="telemetria/descarga")
    ap.add_argument("--sin-dedup", dest="sin_dedup", action="store_true")
    ap.add_argument("--rehacer", action="store_true")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    if a.listar:
        return etapa_listar(a)
    if a.bajar:
        if not (a.desde and a.hasta):
            raise SystemExit("--bajar necesita --desde y --hasta")
        return etapa_bajar(a)
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
