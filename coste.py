# -*- coding: utf-8 -*-
"""
coste.py -- `c(u)`, con PROCEDENCIA por componente. Fuente unica.

    python coste.py --autotest    controles
    python coste.py --leer        lee lo que es publico y lo cachea
    python coste.py --informe     la tabla de c(u) y que criterio falta

QUE ES `c(u)`, PORQUE HAY TRES `c` EN ESTE PROYECTO Y SE CONFUNDEN:

    c(u)     COSTE DE TRANSACCION de ida y vuelta por unidad, en pb.  <- ESTE
             La `u` esta porque depende de la ACCION: maker o taker, y tamano.
             Es el peaje contra el que se compara la senal (paso 3 de la v3.2,
             banda muerta |alpha| > c de la decision de diseno 2).
    sigma^2  coeficiente difusivo de la EDP de Loeper (Sec. 7.4 del PDF). NO es
             c(u). Es la misma sigma que se pronostica en `retroalimentacion.py`.
    c2_vol   `k*omega_m*nu` de la Sec. 1, unidades 1/Anos. Tampoco es c(u).

=========================================================================
LAS CUATRO COMPONENTES, Y CUAL SE PUEDE CERRAR
=========================================================================

    c(u, H) = comision(esquema) + cruce_de_spread(esquema)
              + seleccion_adversa(esquema) + financiacion(H)

| componente        | procedencia         | estado |
|-------------------|---------------------|--------|
| comision          | endpoint FIRMADO    | ⛔ ASUMIDA. Es el 81.8 % del total |
| cruce de spread   | la captura          | ✅ MEDIDO: 0.0146 pb (34.8 M ticks) |
| seleccion adversa | v4.0 Sec.5          | ✅ MEDIDA: 0.888 pb |
| financiacion      | endpoint PUBLICO    | ✅ LEIDA: 0.426 pb / 8 h |

⚠ **EL CRITERIO DEL Sec.8 SIGUE SIN CUMPLIRSE, Y NO SE PUEDE FINGIR.**
`/fapi/v1/commissionRate` es firmado. Sin credenciales de MAINNET con permiso de
LECTURA no hay forma de leer el escalon, y la comision es el **81.8 %** de `c(u)`
-- el cruce del spread medido es **274x menor** que la comision maker de ida y
vuelta. Cerrar las otras tres componentes NO cierra el criterio.

⚠ **Y NO SIRVEN LAS CREDENCIALES DE TESTNET.** El escalon de comisiones es una
propiedad de la CUENTA (nivel VIP, descuento BNB, referido), y la de Testnet no
es la de Mainnet. Este proyecto ya se quemo con esa diferencia: la v1.3 midio
que el `stepSize` de Testnet es 10x mas fino que el de Mainnet y dejo escrito
que calibrar contra el entorno equivocado produce un sistema que funciona en
pruebas y se degrada en produccion. Lo mismo aplica aqui.

⚠ **`c(u)` NO ES UNA CONSTANTE, y esa es la tercera pata.** Si la orden se pone
como maker no hay garantia de que se llene:

    c_efectivo = p * c_maker + (1 - p) * C_respaldo

con `p` = probabilidad de llenado y `C_respaldo` = lo que cuesta el plan B
(cruzar). `cola.py` ya tiene esa estructura; `p` y `C_respaldo` son el Sec.5 de
la v4.1, que sigue pendiente y que a su vez espera a que el Sec.1 diga a que
horizonte. **Todo el Sec.1 se midio contra `c_maker` puro, o sea suponiendo
llenado maker perfecto en las dos patas: es el mejor caso posible.**
"""

from __future__ import annotations

import argparse
import datetime as DT
import json
import os
import ssl
import urllib.request

import numpy as np

import horizonte as H

log, titulo = H.log, H.titulo
CACHE = "telemetria/coste_publico.json"
BASE = "https://fapi.binance.com"

# --- componentes ASUMIDAS -------------------------------------------------
# Binance USD-M futures, VIP 0. NO leidas de la cuenta.
COMISION_MAKER_ASUMIDA = 0.000200      # 0.0200 %
COMISION_TAKER_ASUMIDA = 0.000500      # 0.0500 %  <- CONSERVADORA, ver abajo
COMISIONES_LEIDAS = False              # <- la unica bandera que importa

# --- LEIDO DE TESTNET el 2026-08-23, con credenciales de la cuenta demo -----
# [!] NO CIERRA EL CRITERIO DEL Sec.8 y no se usa por omision. El escalon es
#     propiedad de la CUENTA de Mainnet; una cuenta de Testnet nace siempre en
#     `feeTier = 0` (verificado: `/fapi/v2/account` devolvio 0) y no sabe nada
#     del descuento BNB, del nivel VIP real ni de un referido.
#
# ⚠ PERO CONTRADICE AL REPO EN LA TAKER, y hay que decidirlo:
#     leido en Testnet : maker 0.000200   taker 0.000400
#     asumido en 4 modulos: maker 0.000200   taker 0.000500
#   La maker coincide exacta. La taker asumida es un 25 % MAS ALTA que la leida
#   -- probablemente una tarifa vieja: Binance bajo la taker de futuros de
#   0.0500 % a 0.0400 % en algun momento. Se conserva 0.0005 por omision porque
#   es la CONSERVADORA (mas coste = requisito mas duro = conclusion negativa mas
#   robusta), y se expone la leida para que la proxima sesion decida con el dato
#   de Mainnet delante, no con este.
COMISION_MAKER_TESTNET = 0.000200
COMISION_TAKER_TESTNET = 0.000400
TESTNET_FEE_TIER = 0

# --- LEIDO DE TESTNET: `leverageBracket`, que la v1.3 no pudo leer ----------
# La v1.3 dejo escrito que `mmr` no se puede leer sin credenciales y que
# `mercado.leer_mmr` devuelve el valor por defecto INFLADO por un factor de
# seguridad de 2x. El valor que asumia era `MMR_PRIMER_TRAMO_BTCUSDT = 0.004`
# para nocional <= 50 000, y la lectura lo confirma EXACTO. El factor 2x era
# conservadurismo puro, no ignorancia.
MMR_TRAMOS_TESTNET = [(50000, 0.0040), (250000, 0.0050), (3000000, 0.0100),
                      (20000000, 0.0250), (40000000, 0.0500)]

# --- componentes MEDIDAS --------------------------------------------------
# `curvas_estacional`, 34 812 523 ticks alineados de los cuatro tramos largos.
SPREAD_MEDIANO_PB = 0.0146
# v4.0 Sec.5, markout tras llenado maker.
SELECCION_ADVERSA_PB = 1e4 * 5.78 / 65076.0     # 0.888 pb

ESQUEMAS = ("maker_maker", "maker_taker", "taker_taker")


# ===========================================================================
# Lectura de lo publico
# ===========================================================================

def _get(url, timeout=20):
    req = urllib.request.Request(url, headers={"User-Agent": "micelio/1.0"})
    with urllib.request.urlopen(req, timeout=timeout,
                                context=ssl.create_default_context()) as r:
        return json.loads(r.read().decode())


def leer_financiacion(limite: int = 1000) -> dict:
    """Tasa de financiacion del perpetuo. Endpoint PUBLICO, sin firma."""
    fr = _get("%s/fapi/v1/fundingRate?symbol=BTCUSDT&limit=%d" % (BASE, limite))
    r = np.array([float(x["fundingRate"]) for x in fr])
    t = np.array([int(x["fundingTime"]) for x in fr]) / 1000.0
    reciente = t > (t.max() - 30 * 86400)
    return {"n": int(r.size),
            "desde": DT.datetime.utcfromtimestamp(float(t.min())).isoformat(),
            "hasta": DT.datetime.utcfromtimestamp(float(t.max())).isoformat(),
            "media_pb": float(1e4 * r.mean()),
            "abs_mediana_pb": float(1e4 * np.median(np.abs(r))),
            "p90_abs_pb": float(1e4 * np.percentile(np.abs(r), 90)),
            "abs_mediana_30d_pb": float(1e4 * np.median(np.abs(r[reciente]))),
            "leido_en": DT.datetime.utcnow().isoformat()}


def leer_comision_firmado(clave: str, secreto: str, symbol="BTCUSDT") -> dict:
    """`/fapi/v1/commissionRate`, FIRMADO. Requiere credenciales de MAINNET.

    No se invoca en ningun flujo automatico: hace falta que el operador pase las
    credenciales explicitamente. Una clave de solo LECTURA basta -- este endpoint
    no necesita permiso de trading, y pedir mas permisos de los necesarios es
    exactamente el riesgo que no hay que correr.
    """
    import hashlib
    import hmac
    import time
    q = "symbol=%s&timestamp=%d&recvWindow=5000" % (symbol, int(time.time() * 1000))
    firma = hmac.new(secreto.encode(), q.encode(), hashlib.sha256).hexdigest()
    req = urllib.request.Request("%s/fapi/v1/commissionRate?%s&signature=%s" % (BASE, q, firma),
                                 headers={"X-MBX-APIKEY": clave,
                                          "User-Agent": "micelio/1.0"})
    with urllib.request.urlopen(req, timeout=20,
                                context=ssl.create_default_context()) as r:
        d = json.loads(r.read().decode())
    return {"maker": float(d["makerCommissionRate"]),
            "taker": float(d["takerCommissionRate"]),
            "symbol": d.get("symbol", symbol),
            "leido_en": DT.datetime.utcnow().isoformat()}


def cache() -> dict:
    if os.path.exists(CACHE):
        try:
            return json.load(open(CACHE, encoding="utf-8"))
        except Exception:
            pass
    return {}


# ===========================================================================
# c(u)
# ===========================================================================

def financiacion_pb(H_s: float, por_8h_pb: float = None) -> float:
    """Financiacion ESPERADA de una tenencia de `H_s` segundos.

    [!] LINEAL DESDE CERO, y esto corrige un defecto de `horizonte.py`. Alli se
    escribe `max(H_s - 3600, 0)/(8*3600)`, o sea financiacion CERO por debajo de
    una hora. No tiene justificacion: la financiacion se cobra en marcas fijas
    cada 8 h, asi que una tenencia de `H` cruza una marca con probabilidad
    `H/(8h)` y su coste ESPERADO es lineal en `H` desde 0.
    Numericamente da igual -- a 300 s son 0.0044 pb contra 4 pb de comision,
    900x menos -- pero un cero puesto a ojo es un cero puesto a ojo.
    """
    if por_8h_pb is None:
        por_8h_pb = cache().get("financiacion", {}).get("abs_mediana_pb", 0.426)
    return por_8h_pb * max(H_s, 0.0) / (8 * 3600.0)


def c_u(esquema: str = "maker_maker", H_s: float = 0.0,
        maker: float = None, taker: float = None) -> dict:
    """`c(u)` en pb, desglosado, con su procedencia."""
    if esquema not in ESQUEMAS:
        raise ValueError("esquema desconocido: %s" % esquema)
    m = COMISION_MAKER_ASUMIDA if maker is None else maker
    t = COMISION_TAKER_ASUMIDA if taker is None else taker
    n_taker = {"maker_maker": 0, "maker_taker": 1, "taker_taker": 2}[esquema]
    comision = 1e4 * ((2 - n_taker) * m + n_taker * t)
    # cada pata taker cruza medio spread; la maker no cruza pero sufre markout
    cruce = n_taker * 0.5 * SPREAD_MEDIANO_PB
    adversa = (2 - n_taker) * 0.5 * SELECCION_ADVERSA_PB * 2 / 2.0
    fin = financiacion_pb(H_s)
    total = comision + cruce + adversa + fin
    return {"esquema": esquema, "H_s": H_s, "total_pb": total,
            "comision_pb": comision, "cruce_pb": cruce,
            "seleccion_adversa_pb": adversa, "financiacion_pb": fin,
            "comisiones_leidas": bool(COMISIONES_LEIDAS or maker is not None)}


def criterio_sec8() -> dict:
    """?Se cumple "escalon de comisiones LEIDO de la cuenta, no asumido"?"""
    c = c_u("maker_maker")
    return {"cumplido": bool(COMISIONES_LEIDAS),
            "fraccion_de_c_que_es_asumida": c["comision_pb"] / c["total_pb"],
            "bloqueante": "credenciales de MAINNET con permiso de LECTURA"}


# ===========================================================================
# Etapas
# ===========================================================================

def etapa_leer(args) -> int:
    titulo("LECTURA DE LO QUE ES PUBLICO")
    d = cache()
    try:
        f = leer_financiacion()
        d["financiacion"] = f
        log("  financiacion: %d periodos, %s -> %s" % (f["n"], f["desde"][:10], f["hasta"][:10]))
        log("    media          %+7.4f pb / 8 h   [!] es una TRANSFERENCIA con"
            " signo: un corto la COBRA" % f["media_pb"])
        log("    |tasa| mediana  %7.4f pb / 8 h" % f["abs_mediana_pb"])
        log("    |tasa| p90      %7.4f pb / 8 h" % f["p90_abs_pb"])
        log("    |tasa| 30 dias  %7.4f pb / 8 h" % f["abs_mediana_30d_pb"])
    except Exception as e:
        log("  financiacion FALLO: %r" % (e,))
    os.makedirs(os.path.dirname(CACHE), exist_ok=True)
    json.dump(d, open(CACHE, "w", encoding="utf-8"), indent=1)
    log("  cache: %s" % CACHE)
    log("")
    log("  [!] La comision NO se lee aqui: exige firma. Ver `leer_comision_firmado`.")
    return 0


def etapa_informe(args) -> int:
    titulo("c(u) -- DESGLOSE Y PROCEDENCIA")
    f = cache().get("financiacion")
    log("")
    if f:
        log("  financiacion LEIDA el %s: |tasa| mediana %.4f pb / 8 h  (%d periodos)"
            % (f["leido_en"][:10], f["abs_mediana_pb"], f["n"]))
        log("    constante ASUMIDA que sustituye: 0.9681 pb / 8 h  ->  factor %.2fx"
            % (0.9681 / f["abs_mediana_pb"]))
    else:
        log("  financiacion NO leida todavia: corre `python coste.py --leer`")
    log("")
    log("  c(u) en pb, por esquema y horizonte de tenencia:")
    log("    %-12s %8s | %9s %8s %9s %8s"
        % ("esquema", "H", "comision", "cruce", "adversa", "financ."))
    log("    " + "-" * 62)
    for esq in ESQUEMAS:
        for Hs in (300.0, 3600.0, 14400.0):
            c = c_u(esq, Hs)
            log("    %-12s %6.0f s | %9.4f %8.4f %9.4f %8.4f   TOTAL %8.4f"
                % (esq, Hs, c["comision_pb"], c["cruce_pb"],
                   c["seleccion_adversa_pb"], c["financiacion_pb"], c["total_pb"]))
    log("")
    log("  --- lectura de TESTNET (2026-08-23), que NO cierra el criterio ---")
    log("    maker: leida %.6f  contra asumida %.6f   -> %s"
        % (COMISION_MAKER_TESTNET, COMISION_MAKER_ASUMIDA,
           "COINCIDE" if abs(COMISION_MAKER_TESTNET - COMISION_MAKER_ASUMIDA) < 1e-9
           else "DIFIERE"))
    log("    taker: leida %.6f  contra asumida %.6f   -> %s"
        % (COMISION_TAKER_TESTNET, COMISION_TAKER_ASUMIDA,
           "COINCIDE" if abs(COMISION_TAKER_TESTNET - COMISION_TAKER_ASUMIDA) < 1e-9
           else "*** DIFIERE en %.0f %% ***"
                % (100 * (COMISION_TAKER_ASUMIDA / COMISION_TAKER_TESTNET - 1))))
    ct_a = c_u("taker_taker", 300.0)["total_pb"]
    ct_l = c_u("taker_taker", 300.0, maker=COMISION_MAKER_TESTNET,
               taker=COMISION_TAKER_TESTNET)["total_pb"]
    log("    c(u) taker+taker a 300 s: asumida %.4f pb  contra leida %.4f pb"
        " -> R2_req x %.2f" % (ct_a, ct_l, (ct_l / ct_a) ** 2))
    log("    `feeTier` de la cuenta de Testnet: %d (nace siempre en VIP 0)" % TESTNET_FEE_TIER)
    log("")
    log("  --- mmr, hueco abierto desde la v1.3 y ahora LEIDO (Testnet) ---")
    log("    tramo 1 (nocional <= 50 000): mmr = %.4f" % MMR_TRAMOS_TESTNET[0][1])
    log("    `mercado.MMR_PRIMER_TRAMO_BTCUSDT` asumia 0.0040  ->  COINCIDE EXACTO")
    log("    el factor de seguridad 2x de `mercado.leer_mmr` era conservadurismo,")
    log("    no ignorancia: el numero asumido estaba bien.")
    log("")
    cr = criterio_sec8()
    log("  --- CRITERIO DEL Sec.8: comisiones leidas de la cuenta ---")
    log("    cumplido: %s" % ("SI" if cr["cumplido"] else "*** NO ***"))
    log("    fraccion de c(u) que descansa en un numero ASUMIDO: %.1f %%"
        % (100 * cr["fraccion_de_c_que_es_asumida"]))
    log("    bloqueante: %s" % cr["bloqueante"])
    log("")
    base = c_u("maker_maker", 300.0)["total_pb"]
    log("  --- lo que esto le hace al R2 requerido del Sec.1 ---")
    log("      R2_req ~ c^2, y el Sec.1 se midio entero con maker+maker:")
    for esq in ESQUEMAS:
        t = c_u(esq, 300.0)["total_pb"]
        log("        %-12s c = %7.4f pb  ->  R2_req x %.2f   (a H = 300 s, tramo 1:"
            " %.2f %% -> %.2f %%)" % (esq, t, (t / base) ** 2, 2.04, 2.04 * (t / base) ** 2))
    log("")
    log("  [!] Y falta la TERCERA PATA: c(u) no es constante. Con llenado maker")
    log("      incierto, c_efectivo = p*c_maker + (1-p)*C_respaldo. Todo el Sec.1")
    log("      se midio contra c_maker puro, o sea el MEJOR CASO posible.")
    return 0


# ===========================================================================
# Controles
# ===========================================================================

def _autotest() -> int:
    titulo("CONTROLES DE coste.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-54s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    c = c_u("maker_maker", 0.0)
    chk(abs(c["comision_pb"] - 4.0) < 1e-9, "maker+maker da 4.00 pb de comision",
        "%.4f" % c["comision_pb"])
    ct = c_u("taker_taker", 0.0)
    chk(abs(ct["comision_pb"] - 10.0) < 1e-9, "taker+taker da 10.00 pb (0.05 % por lado)",
        "%.4f" % ct["comision_pb"])
    chk(ct["comision_pb"] > 2 * c["comision_pb"] - 1e-9,
        "taker cuesta mas del doble que maker")
    chk(abs(c["cruce_pb"]) < 1e-12, "maker+maker NO cruza el spread")
    chk(ct["cruce_pb"] > 0, "taker+taker si cruza", "%.4f pb" % ct["cruce_pb"])
    chk(ct["cruce_pb"] < 0.01 * ct["comision_pb"],
        "el cruce es despreciable frente a la comision",
        "%.1fx menor" % (ct["comision_pb"] / max(ct["cruce_pb"], 1e-12)))

    # financiacion lineal desde cero, no con escalon a 1 h
    f1 = financiacion_pb(300.0, 0.426)
    f2 = financiacion_pb(600.0, 0.426)
    chk(abs(f2 - 2 * f1) < 1e-12, "la financiacion es LINEAL en H desde cero",
        "%.6f y %.6f pb" % (f1, f2))
    chk(financiacion_pb(0.0, 0.426) == 0.0, "financiacion nula a H = 0")
    chk(abs(financiacion_pb(8 * 3600.0, 0.426) - 0.426) < 1e-12,
        "a 8 h la financiacion es la tasa entera")

    # la bandera no se puede falsear sin pasar tarifas
    chk(c_u("maker_maker")["comisiones_leidas"] is False,
        "sin credenciales, `comisiones_leidas` es False")
    chk(c_u("maker_maker", maker=0.0001, taker=0.0003)["comisiones_leidas"] is True,
        "pasando tarifas explicitas, se marca como leidas")
    cr = criterio_sec8()
    chk(cr["cumplido"] is False, "el criterio del Sec.8 se reporta INCUMPLIDO")
    chk(cr["fraccion_de_c_que_es_asumida"] > 0.75,
        "y la parte asumida domina c(u)", "%.1f %%" % (100 * cr["fraccion_de_c_que_es_asumida"]))

    log("")
    log("RESULTADO: %d fallo(s)" % fallos)
    return 1 if fallos else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--leer", action="store_true")
    ap.add_argument("--informe", action="store_true")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    if a.leer:
        return etapa_leer(a)
    if a.informe:
        return etapa_informe(a)
    return etapa_informe(a)


if __name__ == "__main__":
    raise SystemExit(main())
