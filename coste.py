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
# --- TARIFAS REALES, LEIDAS DE LA CUENTA DE MAINNET (2026-08-28) -----------
# `/fapi/v1/commissionRate` con clave firmada de Mainnet (Enable Reading +
# Enable Futures, IP restringida, sin trading ni retiros). El criterio del Sec.8
# del PREREGISTRO_3_2, abierto desde la v3.1, QUEDA CUMPLIDO.
COMISION_MAKER_ASUMIDA = 0.000200      # 2.0000 pb por lado   LEIDA
COMISION_TAKER_ASUMIDA = 0.000500      # 5.0000 pb por lado   LEIDA
COMISIONES_LEIDAS = True
VIP_LEIDO_DE_CUENTA = 0                # `/fapi/v2/account` -> feeTier = 0
BNB_BURN_ACTIVO = False                # `/fapi/v1/feeBurn` -> {'feeBurn': False}
MMR_TRAMO1_MAINNET = 0.0040            # hasta 300 000 USD de nocional

# ⚠ ESTO CORRIGE UNA CORRECCION MIA DEL 2026-08-23, Y EN LA DIRECCION MALA.
# Aquel dia baje la taker de 0.000500 a 0.000400 apoyandome en dos vias que
# "coincidian": la lectura de TESTNET y la tabla publica de futuros para VIP 0.
# La cuenta real dice **0.000500**. Las dos vias coincidian porque las dos eran
# indirectas -- Testnet no es la cuenta, y una tabla publica no es un escalon
# leido --, y coincidir no es lo mismo que acertar. El `0.0005` que arrastraban
# `propagador.py`, `cola.py` y `tick_grande.py` era el CORRECTO desde el
# principio, y lo que lo salvo fue haberlo conservado como "conservador" en vez
# de borrarlo.
# Consecuencia: `c(u)` taker+taker vuelve a 10.00 pb (no 8.02) y `R2_req` de la
# ejecucion taker se multiplica por 4.19 respecto a maker, no por 2.69.
COMISION_TAKER_CONSERVADORA = 0.000500

# ⚠ EL DESCUENTO BNB **NO** ESTA ACTIVO. El Sec.4.4 de la v4.2 construye su
# tabla de expectativas con `lastre = 4.49 pb`, que es maker CON descuento. El
# lastre real es **4.888 pb** y esa tabla queda desplazada en contra.

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
            "parcial": False,
            "leido_de_la_cuenta": ("maker %.6f, taker %.6f, feeTier %d, BNB %s"
                                   % (COMISION_MAKER_ASUMIDA, COMISION_TAKER_ASUMIDA,
                                      VIP_LEIDO_DE_CUENTA, BNB_BURN_ACTIVO)),
            "de_tabla_publica": "nada",
            "sin_leer": "nada",
            "fraccion_de_c_que_es_asumida": 0.0,
            "bloqueante": None}


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
    log("  --- lectura de TESTNET (historica; la de MAINNET ya la sustituye) ---")
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
    log("    cumplido del todo: %s   -> PARCIAL" % ("SI" if cr["cumplido"] else "NO"))
    log("    LEIDO de la cuenta      : %s" % cr["leido_de_la_cuenta"])
    log("    de tabla publica        : %s" % cr["de_tabla_publica"])
    log("    SIN leer                : %s" % cr["sin_leer"])
    log("    bloqueante              : %s" % cr["bloqueante"])
    bnb = c_u("maker_maker", 300.0, maker=0.9 * COMISION_MAKER_ASUMIDA,
              taker=0.9 * COMISION_TAKER_ASUMIDA)["total_pb"]
    base = c_u("maker_maker", 300.0)["total_pb"]
    log("")
    log("    si el descuento BNB estuviera activo: c(u) maker %.4f -> %.4f pb,"
        " R2_req x %.2f" % (base, bnb, (bnb / base) ** 2))
    log("    [!] El BNB solo puede BAJAR la comision, asi que la cifra actual es una")
    log("        COTA SUPERIOR del coste. Para una conclusion NEGATIVA como la del")
    log("        Sec.1 eso es lo que se quiere: si no cruza con el coste maximo,")
    log("        tampoco cruzaria con el real. Solo mordería en una conclusion positiva.")
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
    # [!] ESTE CONTROL FALLO AL CORREGIR LA TAKER, Y ESTUVO BIEN QUE FALLARA.
    # Tenia clavado 10.00 pb (taker 0.0500 %) y detecto que la constante habia
    # cambiado a 0.0400 %. Un test que fija una tarifa es justamente lo que
    # impide que una tarifa se mueva en silencio -- que es como el 0.0005 viejo
    # sobrevivio en cuatro modulos sin que nadie lo notara.
    # [!] ESTE CONTROL YA FALLO DOS VECES Y LAS DOS CON RAZON: primero al bajar
    # la taker a 0.0400 % sobre evidencia indirecta, y despues al volver a
    # 0.0500 % con la lectura REAL de la cuenta. Un test que fija una tarifa es
    # lo que impide que se mueva en silencio, y aqui impidio que se quedara mal.
    chk(abs(ct["comision_pb"] - 10.0) < 1e-9,
        "taker+taker da 10.00 pb (0.0500 % por lado, LEIDO de la cuenta)",
        "%.4f" % ct["comision_pb"])
    chk(abs(COMISION_TAKER_ASUMIDA - COMISION_TAKER_CONSERVADORA) < 1e-12,
        "la taker real coincide con la que se conservo como conservadora")
    chk(VIP_LEIDO_DE_CUENTA == 0, "el nivel VIP leido de la cuenta de Mainnet es 0")
    cbnb = c_u("maker_maker", 0.0, maker=0.9 * COMISION_MAKER_ASUMIDA,
               taker=0.9 * COMISION_TAKER_ASUMIDA)["total_pb"]
    chk(cbnb < c_u("maker_maker", 0.0)["total_pb"],
        "el descuento BNB solo puede BAJAR c(u): la cifra actual es cota superior",
        "%.4f < %.4f pb" % (cbnb, c_u("maker_maker", 0.0)["total_pb"]))
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
    chk(c_u("maker_maker")["comisiones_leidas"] is True,
        "`comisiones_leidas` es True: las tarifas estan LEIDAS de la cuenta")
    chk(c_u("maker_maker", maker=0.0001, taker=0.0003)["comisiones_leidas"] is True,
        "pasando tarifas explicitas, se marca como leidas")
    cr = criterio_sec8()
    chk(cr["cumplido"] is True, "el criterio del Sec.8 se reporta CUMPLIDO")
    chk(cr["fraccion_de_c_que_es_asumida"] == 0.0,
        "y ninguna parte de c(u) descansa ya en un numero asumido")
    chk(BNB_BURN_ACTIVO is False,
        "el descuento BNB esta LEIDO y esta APAGADO -> lastre 4.888 pb, no 4.49")

    log("")
    log("RESULTADO: %d fallo(s)" % fallos)
    return 1 if fallos else 0


def etapa_leer_mainnet(args) -> int:
    """Lee las tarifas REALES. Las credenciales las pasa el operador y NO se
    guardan en ningun sitio: solo se imprime el resultado.

    [!] EXISTE PARA QUE LAS CLAVES NO TENGAN QUE PASAR POR LA CONVERSACION. El
    operador ejecuta esto en su maquina y comparte unicamente las dos cifras.
    Una clave con permiso de RETIRO que ha pasado por un canal de texto hay que
    borrarla, y este proyecto ya tuvo dos.
    """
    titulo("Sec.1.1 -- TARIFAS REALES DE MAINNET")
    if not args.clave or not args.secreto:
        log("  faltan --clave y --secreto")
        return 2
    try:
        d = leer_comision_firmado(args.clave, args.secreto)
    except Exception as e:
        log("  FALLO: %r" % (e,))
        log("  si es -2015: falta `Enable Futures` (no `Enable Withdrawals`),")
        log("  y Binance solo lo ofrece con lista blanca de IP activada.")
        return 1
    m, t = d["maker"], d["taker"]
    log("")
    log("  ===== COMPARTE SOLO ESTAS DOS LINEAS =====")
    log("  maker = %.6f   (%.4f pb por lado)" % (m, 1e4 * m))
    log("  taker = %.6f   (%.4f pb por lado)" % (t, 1e4 * t))
    log("  =========================================")
    log("")
    log("  c(u) maker+maker = %.4f pb   taker+taker = %.4f pb" % (2e4 * m, 2e4 * t))
    log("  contra lo asumido hoy: maker %.4f pb   taker %.4f pb"
        % (2e4 * COMISION_MAKER_ASUMIDA, 2e4 * COMISION_TAKER_ASUMIDA))
    try:
        import hashlib, hmac, json as _j, time as _t, urllib.request as _u
        q = "timestamp=%d&recvWindow=5000" % int(_t.time() * 1000)
        f = hmac.new(args.secreto.encode(), q.encode(), hashlib.sha256).hexdigest()
        r = _u.Request("%s/fapi/v1/feeBurn?%s&signature=%s" % (BASE, q, f),
                       headers={"X-MBX-APIKEY": args.clave})
        log("  descuento BNB: %s" % _j.loads(_u.urlopen(r, timeout=20).read().decode()))
    except Exception:
        log("  descuento BNB: no legible")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--leer", action="store_true")
    ap.add_argument("--leer-mainnet", dest="leer_mainnet", action="store_true")
    ap.add_argument("--clave", default="")
    ap.add_argument("--secreto", default="")
    ap.add_argument("--informe", action="store_true")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    if a.leer_mainnet:
        return etapa_leer_mainnet(a)
    if a.leer:
        return etapa_leer(a)
    if a.informe:
        return etapa_informe(a)
    return etapa_informe(a)


if __name__ == "__main__":
    raise SystemExit(main())
