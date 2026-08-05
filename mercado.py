"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: mercado.py — Modo de operación, filtros del instrumento y feed público

Implementa la PRECONDICIÓN y la Sección A de ORDEN_TRABAJO_RIESGO_1_3.md.

RESPONSABILIDADES
-----------------
  1. El tri-estado `MODO ∈ {LECTURA, TESTNET, MAINNET}` que sustituye al booleano
     `IS_TESTNET` (Sec. B.6 del orden de trabajo).
  2. La lectura de `exchangeInfo` — stepSize, minQty, tickSize, minNotional — de
     cualquiera de los dos entornos. NINGUNO de esos valores puede volver a ser
     un literal en el código (Sec. A.2).
  3. El feed público de solo lectura de Mainnet, sin credenciales.

POR QUÉ EL TRI-ESTADO
---------------------
`IS_TESTNET` decidía a la vez TRES cosas que deben moverse por separado: el
generador de precios, el modelo de λ y el modelo de fills. Con un solo booleano
no se puede pedir "precios reales de Mainnet pero sin ejecutar", que es
exactamente lo que las Secciones D y E de la v1.3 necesitan para calibrarse.

    LECTURA  — precios REALES de Mainnet. Sin credenciales. Ejecución
               físicamente imposible: no hay clave con la que firmar y
               `assert_ejecucion_permitida` aborta antes de intentarlo.
    TESTNET  — precios y libro simulados de Binance Futures Testnet. λ del
               proceso OU (Sec. 8.1.1). Valida PLOMERÍA, no calidad de ejecución.
    MAINNET  — dinero real. Fuera del alcance de la v1.3 (Sec. G).

SEGURIDAD
---------
Este módulo NO lee, NO almacena y NO transmite credenciales. Todos los endpoints
que usa son públicos y sin firma. Esa es la propiedad que hace del Modo LECTURA
un modo seguro por construcción y no por disciplina.
"""

from __future__ import annotations

import asyncio
import json
import math
import time
import urllib.request
from dataclasses import dataclass
from enum import Enum

SYMBOL_DEFECTO = "BTCUSDT"


# ==============================================================================
# 1. MODO DE OPERACIÓN (Sec. B.6)
# ==============================================================================
class Modo(str, Enum):
    """Tri-estado que sustituye a `IS_TESTNET`."""

    LECTURA = "LECTURA"
    TESTNET = "TESTNET"
    MAINNET = "MAINNET"


# Códigos numéricos para poder viajar por el bloque de hot-reloading, que es un
# array de float64. El orden es el de riesgo creciente, a propósito.
CODIGO_MODO = {Modo.LECTURA: 0.0, Modo.TESTNET: 1.0, Modo.MAINNET: 2.0}
MODO_DESDE_CODIGO = {v: k for k, v in CODIGO_MODO.items()}


def modo_desde_codigo(codigo: float) -> Modo:
    """Inversa de CODIGO_MODO. Ante un código desconocido devuelve LECTURA.

    La degradación es deliberadamente hacia el modo MENOS capaz: un bloque de
    memoria compartida corrupto o a medio inicializar (todo ceros) debe dejar el
    sistema incapaz de ejecutar, nunca al revés.
    """
    return MODO_DESDE_CODIGO.get(round(float(codigo), 6), Modo.LECTURA)


def ejecucion_permitida(modo: Modo) -> bool:
    """¿Este modo puede firmar y enviar órdenes al exchange?"""
    return modo in (Modo.TESTNET, Modo.MAINNET)


class ErrorDeModo(RuntimeError):
    """Se intentó una acción prohibida en el modo vigente."""


def assert_ejecucion_permitida(modo: Modo, accion: str = "enviar orden") -> None:
    """Compuerta dura antes de cualquier ruta que toque la cuenta.

    ⚠ ES UNA ASERCIÓN, NO UN `if`. En Modo LECTURA la ejecución debe ser
    IMPOSIBLE, no improbable: si alguna vez este camino se alcanza, es un bug de
    control de flujo y hay que verlo reventar, no silenciarlo.
    """
    if not ejecucion_permitida(modo):
        raise ErrorDeModo(
            f"MODO={modo.value} prohibe la accion '{accion}'. El modo LECTURA no "
            f"tiene credenciales cargadas y la ejecucion es fisicamente imposible."
        )


# ==============================================================================
# 2. ENDPOINTS PÚBLICOS (sin firma, sin credenciales)
# ==============================================================================
# LECTURA y MAINNET comparten infraestructura: son el mismo mercado. Lo que los
# distingue no es de dónde se leen los precios, sino si existe una ruta de
# ejecución. Por eso `Modo.LECTURA` apunta a los endpoints de producción.
REST_BASE = {
    Modo.LECTURA: "https://fapi.binance.com",
    Modo.MAINNET: "https://fapi.binance.com",
    Modo.TESTNET: "https://testnet.binancefuture.com",
}

WS_BASE = {
    Modo.LECTURA: "wss://fstream.binance.com/ws",
    Modo.MAINNET: "wss://fstream.binance.com/ws",
    Modo.TESTNET: "wss://fstream.binancefuture.com/ws",
}


def _get_json(url: str, timeout: float = 15.0):
    """GET público. Sin cabeceras de autenticación, por construcción."""
    with urllib.request.urlopen(url, timeout=timeout) as respuesta:
        return json.loads(respuesta.read())


# ==============================================================================
# 3. FILTROS DEL INSTRUMENTO (Sec. A.2)
# ==============================================================================
@dataclass(frozen=True)
class FiltrosInstrumento:
    """Filtros comerciales vigentes de un símbolo. Medidos, nunca asumidos.

    CLAUDE.md documentaba que `Micelio.py` llamaba a
    `apply_filters(u_c, 1e-5, 1e-5, 10.0, P_spot)` — los cuatro valores
    inventados, y los cuatro equivocados por órdenes de magnitud. Un `stepSize`
    de 1e-5 contra el real de 1e-3 hace que el `floor` sea prácticamente la
    identidad, y entonces NADA en las pruebas revela que el controlador continuo
    se convierte en un interruptor al llegar al exchange de verdad.
    """

    symbol: str
    modo: Modo
    step_size: float  # [BTC]     LOT_SIZE.stepSize
    min_qty: float  # [BTC]     LOT_SIZE.minQty
    max_qty: float  # [BTC]     LOT_SIZE.maxQty
    tick_size: float  # [USD]     PRICE_FILTER.tickSize
    min_notional: float  # [USDT]    MIN_NOTIONAL.notional
    ts_lectura: float  # [s] reloj de pared de la lectura

    def lotes_por_nocional(self, nocional_usd: float, precio: float) -> float:
        """Cuántos lotes de `stepSize` caben en un tope de nocional dado."""
        if self.step_size <= 0.0 or precio <= 0.0:
            return 0.0
        return nocional_usd / (self.step_size * precio)

    def cantidad_minima_legal(self, precio: float) -> float:
        """La orden legal más pequeña: satisface a la vez minQty y minNotional.

        No es `min_qty` a secas. Con minNotional = 50 USDT y BTC a 40 000, un
        lote de 0.001 BTC vale 40 USD y es ILEGAL pese a cumplir minQty: hace
        falta redondear hacia arriba al siguiente múltiplo de stepSize que supere
        el nocional mínimo. Esta cantidad SUBE cuando el precio BAJA.
        """
        if precio <= 0.0 or self.step_size <= 0.0:
            return float("inf")
        por_nocional = self.min_notional / precio
        objetivo = max(self.min_qty, por_nocional)
        n_lotes = math.ceil(objetivo / self.step_size - 1e-12)
        return n_lotes * self.step_size

    def cuantizar(self, u_bruto: float, precio: float) -> float:
        """Cuantización comercial (Sec. 8.5). Paso 4 de la cadena de la Sec. B.4.

        ORDEN CORRECTO: se discretiza con `floor` y SOLO DESPUÉS se revalidan
        minQty y minNotional. `apply_filters` arrastraba el bug contrario —
        validaba minNotional ANTES del floor— así que podía devolver una orden
        que tras discretizar ya no cumplía el nocional mínimo y el exchange
        rechazaba.

        Devuelve 0.0 cuando la orden cae por debajo de la resolución del
        instrumento. Ese cero NO es "no operar": es "esta orden es ilegal", y el
        llamador debe contabilizarlo aparte (paso 6 de la Sec. B.4) porque una
        racha de ceros es la firma de un tope por orden mal dimensionado.
        """
        if u_bruto <= 0.0 or self.step_size <= 0.0 or precio <= 0.0:
            return 0.0
        u = math.floor(u_bruto / self.step_size) * self.step_size
        # El floor puede dejar residuos de coma flotante en el último dígito;
        # se redondea a la precisión del propio stepSize.
        decimales = max(0, int(round(-math.log10(self.step_size)))) if self.step_size < 1 else 0
        u = round(u, decimales)
        if u < self.min_qty:
            return 0.0
        if u * precio < self.min_notional:
            return 0.0
        if self.max_qty > 0.0 and u > self.max_qty:
            u = math.floor(self.max_qty / self.step_size) * self.step_size
        return u

    def resumen(self, precio: float) -> str:
        return (
            f"[{self.modo.value}] {self.symbol}: stepSize={self.step_size:g} BTC "
            f"minQty={self.min_qty:g} tickSize={self.tick_size:g} USD "
            f"minNotional={self.min_notional:g} USDT | a S={precio:.2f}: "
            f"1 lote = {self.step_size * precio:.2f} USD, orden legal minima = "
            f"{self.cantidad_minima_legal(precio):g} BTC "
            f"({self.cantidad_minima_legal(precio) * precio:.2f} USD)"
        )


def leer_filtros(
    modo: Modo, symbol: str = SYMBOL_DEFECTO, timeout: float = 15.0
) -> FiltrosInstrumento:
    """Lee `exchangeInfo` del entorno indicado. Endpoint público, sin firma.

    Sec. A.2, tarea 1: "No hardcodear." Se levanta una excepción si el símbolo o
    alguno de los tres filtros no aparece; degradar a valores por defecto sería
    reintroducir el literal por la puerta de atrás.
    """
    datos = _get_json(f"{REST_BASE[modo]}/fapi/v1/exchangeInfo", timeout=timeout)
    info = next((s for s in datos.get("symbols", []) if s["symbol"] == symbol), None)
    if info is None:
        raise ValueError(f"{symbol} no aparece en exchangeInfo de {modo.value}")

    por_tipo = {f["filterType"]: f for f in info.get("filters", [])}
    try:
        lote = por_tipo["LOT_SIZE"]
        precio_f = por_tipo["PRICE_FILTER"]
        nocional = por_tipo["MIN_NOTIONAL"]
    except KeyError as err:
        raise ValueError(
            f"exchangeInfo de {modo.value} no trae el filtro {err}; no se puede "
            f"operar sin conocerlo y no hay valor por defecto defendible."
        ) from err

    return FiltrosInstrumento(
        symbol=symbol,
        modo=modo,
        step_size=float(lote["stepSize"]),
        min_qty=float(lote["minQty"]),
        max_qty=float(lote["maxQty"]),
        tick_size=float(precio_f["tickSize"]),
        min_notional=float(nocional["notional"]),
        ts_lectura=time.time(),
    )


# ------------------------------------------------------------------------------
# MEDICIÓN DEL 2026-08-04 — la Sección A.1 razona sobre cifras desactualizadas
# ------------------------------------------------------------------------------
# Leído de `exchangeInfo` de ambos entornos, con BTC ≈ 63 767 USD:
#
#   | entorno | stepSize | minQty | tickSize | minNotional | 1 lote  |
#   |---------|----------|--------|----------|-------------|---------|
#   | MAINNET | 0.001    | 0.001  | 0.10     | 50 USDT     | 63.8 USD|
#   | TESTNET | 0.0001   | 0.0001 | 0.10     | 50 USDT     |  6.4 USD|
#
# Dos correcciones al documento, ambas materiales:
#
#  1. **minNotional es 50 USDT, no 100.** La primera fila de la tabla de A.1
#     ("ninguna orden es legal") NO se cumple hoy: con 1 lote = 63.8 USD > 50, la
#     orden mínima de Mainnet es legal. Pero el margen es estrecho — por debajo de
#     BTC = 50 000 USD un solo lote deja de alcanzar el nocional mínimo y la orden
#     legal más pequeña pasa a ser 0.002 BTC. Por eso `cantidad_minima_legal`
#     depende del precio y se reevalúa, en vez de fijarse una vez.
#
#  2. ⚠ **TESTNET ES 10× MÁS FINO QUE MAINNET.** Esto es una trampa, y de las
#     silenciosas: la guarda de resolución de la Sec. A.3 evaluada contra Testnet
#     da 476 lotes y pasa cómodamente, mientras que la misma configuración en
#     Mainnet da 47 y va mucho más justa. Calibrar el dimensionamiento contra
#     Testnet produciría un sistema que funciona en pruebas y se degrada a un
#     interruptor en producción — exactamente el fallo que la Sección A existe
#     para prevenir. Por eso `verificar_resolucion_control` se evalúa SIEMPRE
#     contra el stepSize de MAINNET, sea cual sea el modo vigente.
#
# TODO(revisión periódica): minNotional y stepSize los cambia Binance sin previo
# aviso. Por eso se refrescan en caliente (Sec. A.2, tarea 2) y no se congelan.


def leer_filtros_ambos_entornos(
    symbol: str = SYMBOL_DEFECTO, timeout: float = 15.0
) -> dict:
    """Sec. A.2, tarea 1: leer `exchangeInfo` de Testnet **y** de Mainnet.

    Se leen ambos aunque solo se opere en uno, porque la comparación es lo que
    delata la discrepancia de granularidad documentada arriba.
    """
    salida = {}
    for modo in (Modo.MAINNET, Modo.TESTNET):
        try:
            salida[modo] = leer_filtros(modo, symbol=symbol, timeout=timeout)
        except Exception as err:  # Red caída: se reporta, no se inventa un valor
            salida[modo] = err
    return salida


# ------------------------------------------------------------------------------
# Maintenance margin rate (Sec. B.3)
# ------------------------------------------------------------------------------
# `GET /fapi/v1/leverageBracket` es un endpoint FIRMADO (USER_DATA): exige clave
# y HMAC. En Modo LECTURA no hay credenciales por diseño, así que ese dato NO se
# puede obtener y el documento pide explícitamente "LEERLO ..., no asumir el valor
# por defecto".
#
# Resolución: `leer_mmr` acepta un lector inyectado. Cuando hay credenciales
# (TESTNET/MAINNET) el Motor de Red le pasa una función que consulta el endpoint
# firmado; cuando no las hay, devuelve el valor por defecto del primer tramo de
# BTCUSDT y DEJA CONSTANCIA. La guarda de la Sec. B.3 se vuelve conservadora ante
# la duda: un mmr subestimado agranda el colchón calculado y podría dejar pasar
# una configuración que en realidad no muerde antes de la liquidación, así que el
# valor por defecto se infla con un factor de seguridad.
MMR_PRIMER_TRAMO_BTCUSDT = 0.004  # Documentado por Binance para el tramo <= 50k USD
FACTOR_SEGURIDAD_MMR = 2.0  # Se usa 2·mmr cuando el valor es asumido, no leído


def leer_mmr(lector_firmado=None) -> tuple[float, bool]:
    """Devuelve (mmr, fue_leido). `fue_leido=False` significa asumido y ya inflado.

    El factor de seguridad va aquí y no en la guarda para que la guarda reciba un
    número honesto y no tenga que saber de dónde salió.
    """
    if lector_firmado is not None:
        try:
            return float(lector_firmado()), True
        except Exception:
            pass
    return MMR_PRIMER_TRAMO_BTCUSDT * FACTOR_SEGURIDAD_MMR, False


# ==============================================================================
# 4. FEED PÚBLICO DE SOLO LECTURA
# ==============================================================================
# El Modo LECTURA es la PRECONDICIÓN del documento entero: las Secciones D y E se
# calibran sobre datos reales de Mainnet, no sobre Testnet (libro simulado) ni
# sobre mocks (los pusimos nosotros, así que no prueban nada).
#
# Se usa `aggTrade`, no `bookTicker`: la Sec. 1.1 del PDF define ΣQ como volumen
# TRANSADO, y aggTrade lo trae en cada mensaje. Con bookTicker habría que
# inventarse el volumen, que es justo lo que el mock hacía.


@dataclass
class TickMercado:
    """Una transaccion real. `es_sintetico` marca las que no vienen del exchange.

    v2.0 §3.1: lleva IDENTIDAD (`trade_id`) y lado del taker (`es_maker`). Sin la
    identidad no se puede deduplicar (§3.3) ni detectar huecos (§3.4), y hasta la
    v1.3 el sistema sencillamente no la tenia.
    """

    precio: float  # [USD/BTC]
    cantidad: float  # [BTC] volumen de la transaccion
    ts_evento: float  # [s] reloj del EXCHANGE, no el local
    trade_id: int = 0  # 'a' de aggTrade o 't' de trade; 0 = sin identidad
    es_maker: bool = False  # 'm': el comprador es maker -> el taker vendio
    es_sintetico: bool = False


class EstancamientoFeed(RuntimeError):
    """El socket está conectado pero no entrega datos. Ver `FeedPublico`."""


# ------------------------------------------------------------------------------
# ⚠ EL STREAM DE TRANSACCIONES: POR QUE `@trade` Y NO `@aggTrade`
# ------------------------------------------------------------------------------
# La v1.3 concluyó que el WebSocket de futuros estaba filtrado por la red, porque
# `btcusdt@aggTrade` conectaba y callaba mientras el de spot funcionaba. **Esa
# conclusión era incorrecta**, y sondearlo en serio (2026-08-04) lo demuestra:
#
#   EN EL MISMO HOST fstream.binance.com, EN EL MISMO TIPO DE SOCKET:
#     btcusdt@bookTicker    ->  91-288 msgs/s   HABLA
#     btcusdt@trade         ->    26-32 msgs/s  HABLA
#     btcusdt@depth@100ms   ->      9.8 msgs/s  HABLA
#     btcusdt@aggTrade      ->        0         MUDO
#     btcusdt@markPrice     ->        0         MUDO
#     btcusdt@kline_1m      ->        0         MUDO
#     btcusdt@ticker        ->        0         MUDO
#
# Y con suscripción explícita por mensaje al mismo socket, el servidor
# **CONFIRMA** la suscripción — `{"result":["btcusdt@aggTrade"],"id":99}` — y aun
# así no manda un solo dato. Es decir: el socket es bidireccional, la ruta no está
# bloqueada, Binance nos oye y nos contesta. Descartado el filtrado de ISP, y con
# él la necesidad de montar un VPS.
#
# Los certificados son legítimos (DigiCert/GeoTrust de Binance) y no hay proxy en
# el entorno, así que tampoco es inspección TLS. La causa última del silencio de
# `@aggTrade` queda SIN EXPLICAR, y se deja anotada como tal en vez de inventarle
# una: lo que importa operativamente es que hay un stream equivalente que sí
# fluye.
#
# `@trade` da todo lo que `@aggTrade` daba, sin agregar: precio, cantidad, hora
# del exchange, identidad (`t`) y lado del taker (`m`). Medido sobre 45 s:
#     26.0 transacciones/s, tau_d p50 = 0.198 s (p90 0.274, max 0.522)
#     0 huecos en 1170 ids consecutivos
#     cobertura 0.99x -> tiempo real, no replay
# Contra las 0.76 "mediciones"/s de la v1.3, son **34x mas datos**.
#
# TODO(revisar): si algun dia `@aggTrade` vuelve a fluir, `@trade` sigue siendo
# preferible para el reloj de transacciones — agregar trades es justo lo que
# `aggTrade` hace y lo que §1.2 quiere deshacer.
STREAM_TRANSACCIONES = "trade"


# Segundos sin un solo mensaje que bastan para declarar el socket muerto.
# BTCUSDT perpetuo hace decenas de trades por segundo a cualquier hora, así que
# un silencio de este orden no es un mercado tranquilo: es un socket colgado.
TIMEOUT_ESTANCAMIENTO = 12.0


class FeedPublico:
    """Cliente de mercado público. Sin credenciales, sin ruta de ejecución.

    Se implementa con `aiohttp` (disponible en el entorno) en vez de `websockets`
    (no instalado).

    ⚠ FALLO REAL MEDIDO EL 2026-08-04, y de los peligrosos. Desde esta red,
    `wss://fstream.binance.com` (futuros) **completa el handshake y luego no
    entrega ni un solo mensaje**, mientras `wss://stream.binance.com:9443` (spot)
    funciona con normalidad. Es decir: la conexión se establece, no hay
    excepción, no hay cierre, no hay error — y el bot se queda ciego para
    siempre creyendo que el mercado está en calma.

    Ese modo de fallo derrota al watchdog de la Sec. 8.4.2 tal como estaba
    escrito, porque ese watchdog solo reacciona a EXCEPCIONES. Un socket sano que
    no habla no levanta ninguna. En la corrida de verificación se vio el efecto
    completo: `paquete_valido=False` en todos los ciclos, el NIS en `nan`, el
    burn-in eternamente en racha 0 y Tr(P) creciendo sin límite hasta 6e9.

    Defensa, en dos capas:
      1. DETECTOR DE ESTANCAMIENTO. Si pasan `TIMEOUT_ESTANCAMIENTO` segundos sin
         un mensaje, se levanta `EstancamientoFeed` — que el watchdog SÍ ve.
      2. DEGRADACIÓN A SONDEO REST tras varios estancamientos seguidos. El REST
         de futuros responde con normalidad desde la misma red, así que el bot
         sigue operando sobre datos reales de Mainnet en vez de quedarse ciego o
         caer al mock.
    """

    def __init__(self, modo: Modo, symbol: str = SYMBOL_DEFECTO, max_estancamientos: int = 2):
        self.modo = modo
        self.symbol = symbol.lower()
        self.symbol_mayus = symbol.upper()
        self.url_ws = f"{WS_BASE[modo]}/{self.symbol}@{STREAM_TRANSACCIONES}"
        self.url_precio = (
            f"{REST_BASE[modo]}/fapi/v1/ticker/price?symbol={self.symbol_mayus}"
        )
        self.url_trades = (
            f"{REST_BASE[modo]}/fapi/v1/aggTrades?symbol={self.symbol_mayus}&limit=500"
        )
        self.max_estancamientos = int(max_estancamientos)
        self.estancamientos = 0
        self.modo_degradado = False

    async def escuchar(self, al_recibir, cancelado=None) -> None:
        """Bucle de ingesta. `al_recibir(TickMercado)` se invoca por cada tick.

        No reimplementa el watchdog: la reconexión con backoff exponencial de la
        Sec. 8.4.2 vive en el Motor de Red, que es quien debe marcar el flag de
        dropout. Aquí una excepción se PROPAGA a propósito.
        """
        if self.modo_degradado:
            await self._escuchar_rest(al_recibir, cancelado)
            return
        try:
            import aiohttp
        except ImportError:
            self.modo_degradado = True
            await self._escuchar_rest(al_recibir, cancelado)
            return

        try:
            async with aiohttp.ClientSession() as sesion:
                async with sesion.ws_connect(self.url_ws, heartbeat=30.0) as ws:
                    while True:
                        if cancelado is not None and cancelado():
                            return
                        try:
                            mensaje = await asyncio.wait_for(
                                ws.receive(), TIMEOUT_ESTANCAMIENTO
                            )
                        except asyncio.TimeoutError:
                            self.estancamientos += 1
                            if self.estancamientos >= self.max_estancamientos:
                                self.modo_degradado = True
                            raise EstancamientoFeed(
                                f"{TIMEOUT_ESTANCAMIENTO:.0f} s sin datos sobre un "
                                f"socket conectado a {self.url_ws} "
                                f"(estancamiento {self.estancamientos}/"
                                f"{self.max_estancamientos})"
                            ) from None
                        if mensaje.type in (
                            aiohttp.WSMsgType.CLOSED,
                            aiohttp.WSMsgType.CLOSING,
                            aiohttp.WSMsgType.ERROR,
                        ):
                            raise ConnectionError(f"WS cerrado: {mensaje.type}")
                        if mensaje.type is not aiohttp.WSMsgType.TEXT:
                            continue
                        dato = json.loads(mensaje.data)
                        if dato.get("e") not in ("trade", "aggTrade"):
                            continue
                        self.estancamientos = 0
                        # 't' en `@trade`, 'a' en `@aggTrade`: mismo papel de
                        # identidad para deduplicar y detectar huecos (§3.3/§3.4).
                        tid = int(dato.get("t", dato.get("a", 0)))
                        al_recibir(
                            TickMercado(
                                precio=float(dato["p"]),
                                cantidad=float(dato["q"]),
                                # 'T' es el instante del trade en el reloj del
                                # exchange. Usar el local aquí falsearía τ_d de la
                                # Sec. 6.5, que es precisamente la latencia que se
                                # quiere medir.
                                ts_evento=float(dato["T"]) / 1000.0,
                                trade_id=tid,
                                es_maker=bool(dato.get("m", False)),
                            )
                        )
        finally:
            pass

    # --- Sondeo REST: elección de endpoint, MEDIDA no supuesta ----------------
    # Los tres candidatos, medidos el 2026-08-04 con sesión reutilizada, sobre
    # `tau_d` = (reloj local al recibir) − (timestamp que trae el propio dato),
    # que es exactamente el retardo del sensor de la Sec. 6.5:
    #
    #   | endpoint              | peso | tau_d mediana | volumen |
    #   |-----------------------|------|---------------|---------|
    #   | ticker/price          |   1  |    7.64 s     |   no    |
    #   | ticker/bookTicker     |   2  |    1.54 s     |   no    |
    #   | aggTrades             |  20  |    0.70 s     |   SI    |
    #
    # ⚠ `ticker/price` ERA LA ELECCIÓN OBVIA Y ES LA EQUIVOCADA. Su campo `time`
    # no es la hora del servidor sino la del último cambio de precio publicado
    # para el símbolo, y llega con ~7.6 s de rezago. Contra el τ_max = 2 s de la
    # Sec. 6.5 eso significa que **el filtro rechaza el 100 % de los paquetes**:
    # en la corrida de verificación el NIS quedó en `nan` de principio a fin, el
    # burn-in eternamente en racha 0 y Tr(P) creció hasta 6e9 sin una sola
    # corrección. El síntoma no apuntaba al endpoint por ninguna parte.
    #
    # Se usa `aggTrades`, que además trae el volumen: sin él ΣQ queda en cero y
    # con ello mueren Φ, Ψ y Ω, o sea todo el acoplamiento endógeno de la
    # Sec. 1.4 — la tesis del sistema.
    #
    # Presupuesto: peso 20 a 1 Hz = 1 200/min contra 2 400/min disponibles. El
    # 50 %, y deja sitio para las órdenes (peso 20 cada una, ~1/s).
    PERIODO_TRADES_REST = 1.0  # [s]

    async def _escuchar_rest(self, al_recibir, cancelado=None) -> None:
        """Degradación por sondeo de `aggTrades`. Precio Y volumen REALES.

        Se emite un tick por cada trade agregado nuevo, con SU propio precio y SU
        propio timestamp: no se promedia ni se resume. Así `n_ticks` sigue
        contando trades reales —que es lo que ν mide (Sec. 4.6)— y ΣQ acumula el
        volumen real.

        Se usa `aiohttp` y no `urllib` porque la sesión reutiliza la conexión TLS:
        medido, 0.27 s por llamada contra ~0.9 s abriendo un socket nuevo cada
        vez. Con el sondeo a 1 Hz, esa diferencia es la que decide si el retardo
        del sensor cabe o no bajo τ_max.
        """
        try:
            import aiohttp
        except ImportError:
            await self._escuchar_rest_urllib(al_recibir, cancelado)
            return

        ultimo_agg_id = 0
        async with aiohttp.ClientSession() as sesion:
            while True:
                if cancelado is not None and cancelado():
                    return
                async with sesion.get(self.url_trades) as respuesta:
                    trades = await respuesta.json()
                nuevos = [t for t in trades if int(t["a"]) > ultimo_agg_id]
                if trades:
                    ultimo_agg_id = max(int(t["a"]) for t in trades)
                for tr in nuevos:
                    al_recibir(
                        TickMercado(
                            precio=float(tr["p"]),
                            cantidad=float(tr["q"]),
                            ts_evento=float(tr["T"]) / 1000.0,
                            trade_id=int(tr["a"]),
                            es_maker=bool(tr.get("m", False)),
                            # No es sintético: precio y volumen son reales. Lo
                            # degradado es la resolución temporal, no la
                            # procedencia del dato.
                            es_sintetico=False,
                        )
                    )
                await asyncio.sleep(self.PERIODO_TRADES_REST)

    async def _escuchar_rest_urllib(self, al_recibir, cancelado=None) -> None:
        """Último recurso sin aiohttp. Mismo endpoint, sin reutilizar conexión."""
        lazo = asyncio.get_running_loop()
        ultimo_agg_id = 0
        while True:
            if cancelado is not None and cancelado():
                return
            trades = await lazo.run_in_executor(None, _get_json, self.url_trades)
            nuevos = [t for t in trades if int(t["a"]) > ultimo_agg_id]
            if trades:
                ultimo_agg_id = max(int(t["a"]) for t in trades)
            for tr in nuevos:
                al_recibir(
                    TickMercado(
                        precio=float(tr["p"]),
                        cantidad=float(tr["q"]),
                        ts_evento=float(tr["T"]) / 1000.0,
                        trade_id=int(tr["a"]),
                        es_maker=bool(tr.get("m", False)),
                        es_sintetico=True,  # el retardo extra sí es artefacto
                    )
                )
            await asyncio.sleep(self.PERIODO_TRADES_REST)


async def offset_reloj(modo: Modo, timeout: float = 10.0) -> float:
    """Guarda 7 de la Sec. B.5: |offset| contra `/fapi/v1/time`, en segundos.

    Se mide con el round-trip descontado (mitad de la latencia ida y vuelta), que
    es la corrección estándar de NTP. Sin ella, una conexión lenta se confundiría
    con un reloj desviado y la guarda dispararía por la razón equivocada.
    """
    lazo = asyncio.get_running_loop()
    t0 = time.time()
    dato = await lazo.run_in_executor(
        None, _get_json, f"{REST_BASE[modo]}/fapi/v1/time", timeout
    )
    t1 = time.time()
    servidor = float(dato["serverTime"]) / 1000.0
    return (t0 + t1) / 2.0 - servidor
