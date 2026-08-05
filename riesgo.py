"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: riesgo.py — Capa de riesgo de cuenta

Implementa la Sección B de ORDEN_TRABAJO_RIESGO_1_3.md.

UBICACIÓN Y CONTRATO (Sec. B.4)
-------------------------------
Esta capa vive en el MOTOR DE RED, como último filtro antes de firmar, y debe
poder detener el sistema SIN COOPERACIÓN del Hilo Rápido. Si viviera dentro del
NMPC no sería una capa de seguridad: sería parte de lo que debe vigilar.

Ninguna de sus constantes aparece en una ecuación del PDF. El requisito de diseño
es que siga funcionando aunque el modelo esté completamente equivocado — y por eso
no importa `Micelio.py` ni depende de nada del modelo salvo el precio y el
inventario.

SEMÁNTICA DEL HALT: CERRAR Y PARAR, NO CONGELAR (Sec. B.5)
----------------------------------------------------------
Congelar deja exposición abierta sin supervisión. Cerrar tiene además la virtud de
ejercitar la ruta de cierre, que de otro modo es código que nunca se prueba. Si el
cierre falla: reintentar con backoff y escalar la alerta; JAMÁS dar el halt por
completado sin posición plana confirmada por el exchange.
"""

from __future__ import annotations

import math
import time
from collections import deque
from enum import IntEnum

import constantes_micelio as CTE


# ==============================================================================
# B.5 — LAS SIETE GUARDAS, CADA UNA CON SU CAUSA DISTINGUIBLE
# ==============================================================================
class CausaHalt(IntEnum):
    """Cada guarda es un motivo de halt DISTINTO y debe registrarse como tal.

    La distinción no es cosmética: la compuerta de reaprovisionamiento de la
    Sec. C.3 solo recarga sobre la guarda 1. Las guardas 2 a 7 son fallos de
    sistema, y recargar sobre ellas es tapar el bug con dinero.
    """

    NINGUNA = 0
    DRAWDOWN = 1  # equity − equity_inicio <= −PERDIDA_MAX_EPISODIO
    EXPOSICION = 2  # |inventario|·S > nocional_max_posicion(S)
    TASA_ORDENES = 3  # > N_ORDENES_MAX_MIN en ventana de 60 s
    RECHAZOS = 4  # M_RECHAZOS_MAX consecutivos del exchange
    ESTADO_DESCONOCIDO = 5  # órdenes en vuelo sin confirmar > T_CONFIRM_MAX
    DESINCRONIA = 6  # |inv_local − inv_exchange| > TOL_INV
    RELOJ = 7  # |offset| > OFFSET_RELOJ_MAX contra /fapi/v1/time


DESCRIPCION_CAUSA = {
    CausaHalt.NINGUNA: "sin halt",
    CausaHalt.DRAWDOWN: "drawdown de equity (perdida del modelo, no fallo de sistema)",
    CausaHalt.EXPOSICION: "exposicion por encima de la abrazadera de posicion",
    CausaHalt.TASA_ORDENES: "tasa de ordenes anomala (probable bucle)",
    CausaHalt.RECHAZOS: "rechazos consecutivos del exchange",
    CausaHalt.ESTADO_DESCONOCIDO: "ordenes en vuelo sin confirmacion",
    CausaHalt.DESINCRONIA: "inventario local desincronizado del exchange",
    CausaHalt.RELOJ: "reloj local desviado del servidor",
}

# Solo la guarda 1 describe una pérdida del MODELO. El resto son fallos de
# sistema, y el documento es explícito: recargar sobre ellas tapa el bug.
CAUSAS_DE_SISTEMA = frozenset(
    {
        CausaHalt.EXPOSICION,
        CausaHalt.TASA_ORDENES,
        CausaHalt.RECHAZOS,
        CausaHalt.ESTADO_DESCONOCIDO,
        CausaHalt.DESINCRONIA,
        CausaHalt.RELOJ,
    }
)


class MotivoNoEnvio(IntEnum):
    """Por qué una orden candidata no llegó a enviarse. Se contabiliza aparte."""

    ENVIADA = 0
    SIN_SENAL = 1  # el NMPC no pidió nada
    BAJO_RESOLUCION = 2  # paso 6 de la Sec. B.4: el floor la dejó en cero
    CLAMP_POSICION = 3  # la abrazadera de exposición la anuló
    HALT = 4  # el sistema ya estaba detenido
    MODO_SIN_EJECUCION = 5  # Modo LECTURA


# ==============================================================================
# CAPA DE RIESGO
# ==============================================================================
class CapaRiesgo:
    """Abrazadera externa entre el Ring Buffer y la firma de la orden.

    Estado mínimo y explícito: cada guarda tiene su propio contador, y ninguno se
    deriva de otro. Un estado compartido entre guardas haría que un bug en una
    enmascarara a las demás.
    """

    def __init__(
        self,
        filtros,
        nocional_max_orden: float = CTE.NOCIONAL_MAX_ORDEN,
        perdida_max: float = CTE.PERDIDA_MAX_EPISODIO,
        n_ordenes_max_min: int = CTE.N_ORDENES_MAX_MIN,
        m_rechazos_max: int = CTE.M_RECHAZOS_MAX,
        t_confirm_max: float = CTE.T_CONFIRM_MAX,
        tol_inv: float = CTE.TOL_INV,
        offset_reloj_max: float = CTE.OFFSET_RELOJ_MAX,
    ):
        self.filtros = filtros
        self.nocional_max_orden = float(nocional_max_orden)
        self.perdida_max = float(perdida_max)
        self.n_ordenes_max_min = int(n_ordenes_max_min)
        self.m_rechazos_max = int(m_rechazos_max)
        self.t_confirm_max = float(t_confirm_max)
        self.tol_inv = float(tol_inv)
        self.offset_reloj_max = float(offset_reloj_max)
        self.reiniciar()

    # ------------------------------------------------------------------
    # Reset limpio (Sec. C.5): el estado que sobrevive a un halt es fuente
    # segura de confusión al analizar el episodio siguiente.
    # ------------------------------------------------------------------
    def reiniciar(self, equity_inicio: float | None = None) -> None:
        self.equity_inicio = equity_inicio
        self.equity_actual = equity_inicio
        self.ts_equity = 0.0
        self.inventario_local = 0.0
        self.ts_ordenes = deque()  # marcas de tiempo para la ventana de 60 s
        self.rechazos_consecutivos = 0
        self.en_vuelo = {}  # id_cliente -> ts de emisión
        self.offset_reloj = 0.0
        self.ts_offset = 0.0
        self.desvio_inventario = 0.0
        self.causa_halt = CausaHalt.NINGUNA
        self.detalle_halt = ""
        self.n_enviadas = 0
        self.n_rechazadas = 0
        self.n_bajo_resolucion = 0

    # ------------------------------------------------------------------
    # Alimentación de estado
    # ------------------------------------------------------------------
    def registrar_equity(self, equity: float, ts: float | None = None) -> None:
        """Equity del `ACCOUNT_UPDATE` del User Data Stream.

        ⚠ LA PÉRDIDA SE MIDE SOBRE EQUITY, NO SOBRE PnL REALIZADO (Sec. B.5). Si
        solo se cuenta lo realizado, una posición abierta sangra indefinidamente
        más allá del tope sin disparar nada: el drawdown no existe hasta que
        cierras, y para entonces ya te pasaste.

        ⚠ Y VIENE DEL EXCHANGE, NO DE UN CÁLCULO LOCAL. El cálculo local es
        precisamente lo que la guarda 6 (desincronía) existe para desconfiar; si
        el equity saliera de él, un bug de contabilidad desactivaría a la vez el
        kill switch y la guarda que debía detectarlo.
        """
        self.equity_actual = float(equity)
        self.ts_equity = ts if ts is not None else time.time()
        if self.equity_inicio is None:
            self.equity_inicio = self.equity_actual

    def registrar_envio(self, id_cliente: str, ts: float | None = None) -> None:
        t = ts if ts is not None else time.time()
        self.ts_ordenes.append(t)
        self.en_vuelo[id_cliente] = t
        self.n_enviadas += 1

    def registrar_confirmacion(self, id_cliente: str, inventario_nuevo: float) -> None:
        self.en_vuelo.pop(id_cliente, None)
        self.rechazos_consecutivos = 0
        self.inventario_local = float(inventario_nuevo)

    def registrar_rechazo(self, id_cliente: str) -> None:
        self.en_vuelo.pop(id_cliente, None)
        self.rechazos_consecutivos += 1
        self.n_rechazadas += 1

    def registrar_reconciliacion(self, inv_exchange: float) -> None:
        self.desvio_inventario = abs(self.inventario_local - float(inv_exchange))

    def registrar_offset_reloj(self, offset: float, ts: float | None = None) -> None:
        self.offset_reloj = float(offset)
        self.ts_offset = ts if ts is not None else time.time()

    # ------------------------------------------------------------------
    # B.5 — evaluación de las siete guardas
    # ------------------------------------------------------------------
    def evaluar(self, S: float, ahora: float | None = None) -> CausaHalt:
        """Evalúa las siete guardas. Devuelve la primera que dispare.

        El orden NO es arbitrario: se evalúan primero las que describen daño ya
        ocurrido (1, 2) y después las que describen un sistema en mal estado
        (3-7). Así, si el bot está a la vez sangrando y desincronizado, la causa
        registrada es la pérdida, que es la que el operador debe mirar primero.
        """
        t = ahora if ahora is not None else time.time()
        if self.causa_halt != CausaHalt.NINGUNA:
            return self.causa_halt

        # 1. Drawdown de equity.
        if self.equity_inicio is not None and self.equity_actual is not None:
            delta = self.equity_actual - self.equity_inicio
            if delta <= -self.perdida_max:
                return self._disparar(
                    CausaHalt.DRAWDOWN,
                    f"equity {self.equity_actual:.2f} contra inicio "
                    f"{self.equity_inicio:.2f} -> drawdown {delta:.2f} USD "
                    f"(cap {self.perdida_max:.0f})",
                )

        # 2. Exposición, tras reconciliación.
        tope = CTE.nocional_max_posicion(S)
        expuesto = abs(self.inventario_local) * S
        if S > 0.0 and expuesto > tope:
            return self._disparar(
                CausaHalt.EXPOSICION,
                f"|inv|*S = {expuesto:.0f} USD supera la abrazadera {tope:.0f} USD "
                f"(inv={self.inventario_local:.6f} BTC a S={S:.2f})",
            )

        # 3. Tasa de órdenes en ventana de 60 s.
        while self.ts_ordenes and (t - self.ts_ordenes[0]) > 60.0:
            self.ts_ordenes.popleft()
        if len(self.ts_ordenes) > self.n_ordenes_max_min:
            return self._disparar(
                CausaHalt.TASA_ORDENES,
                f"{len(self.ts_ordenes)} ordenes en 60 s (max "
                f"{self.n_ordenes_max_min}). Probable bucle, no un NMPC malo.",
            )

        # 4. Rechazos consecutivos.
        if self.rechazos_consecutivos >= self.m_rechazos_max:
            return self._disparar(
                CausaHalt.RECHAZOS,
                f"{self.rechazos_consecutivos} rechazos consecutivos del exchange "
                f"(max {self.m_rechazos_max})",
            )

        # 5. Estado desconocido: órdenes en vuelo sin confirmación.
        for id_cliente, ts_envio in self.en_vuelo.items():
            if (t - ts_envio) > self.t_confirm_max:
                return self._disparar(
                    CausaHalt.ESTADO_DESCONOCIDO,
                    f"orden {id_cliente} lleva {t - ts_envio:.1f} s sin confirmar "
                    f"(max {self.t_confirm_max:.1f} s). No se sabe si hay posicion.",
                )

        # 6. Desincronía de inventario.
        if self.desvio_inventario > self.tol_inv:
            return self._disparar(
                CausaHalt.DESINCRONIA,
                f"|inv_local - inv_exchange| = {self.desvio_inventario:.6f} BTC "
                f"supera la tolerancia {self.tol_inv:.6f} BTC",
            )

        # 7. Reloj.
        if abs(self.offset_reloj) > self.offset_reloj_max:
            return self._disparar(
                CausaHalt.RELOJ,
                f"offset de reloj {self.offset_reloj*1e3:.0f} ms supera "
                f"{self.offset_reloj_max*1e3:.0f} ms contra /fapi/v1/time",
            )

        return CausaHalt.NINGUNA

    def _disparar(self, causa: CausaHalt, detalle: str) -> CausaHalt:
        self.causa_halt = causa
        self.detalle_halt = detalle
        return causa

    @property
    def en_halt(self) -> bool:
        return self.causa_halt != CausaHalt.NINGUNA

    # ------------------------------------------------------------------
    # B.4 — cadena de operaciones por orden candidata
    # ------------------------------------------------------------------
    def preparar_orden(self, u_compra: float, u_venta: float, S: float, inventario: float):
        """Pasos 2 a 6 de la Sec. B.4. Devuelve (u_c, u_v, motivo).

        Orden estricto:
            1. leer u del Ring Buffer                     (lo hace el llamador)
            2. clamp de nocional  -> u <- min(u, NOCIONAL_MAX_ORDEN / S)
            3. clamp de posicion  -> u <- min(u, tope_posicion/S - |inv|)
            4. apply_filters (floor a stepSize, minQty, minNotional)
            5. RE-validar nocional y posicion DESPUES de cuantizar
            6. si u == 0 tras el floor -> no enviar, contar "bajo resolucion"
            7. firmar y enviar                            (lo hace el llamador)

        El paso 5 no es redundante. `apply_filters` ya arrastraba el bug de
        validar minNotional antes del floor; aquí la cuantización solo puede
        reducir, así que no puede violar un tope superior — pero la revalidación
        es lo que documenta que se pensó, y protege contra un cambio futuro que
        haga la cuantización redondear hacia arriba.
        """
        if self.en_halt:
            return 0.0, 0.0, MotivoNoEnvio.HALT
        if S <= 0.0:
            return 0.0, 0.0, MotivoNoEnvio.SIN_SENAL
        u_c = max(0.0, float(u_compra))
        u_v = max(0.0, float(u_venta))
        if u_c <= 0.0 and u_v <= 0.0:
            return 0.0, 0.0, MotivoNoEnvio.SIN_SENAL

        # --- Paso 2: clamp de nocional por orden ---
        tope_orden_btc = self.nocional_max_orden / S
        u_c = min(u_c, tope_orden_btc)
        u_v = min(u_v, tope_orden_btc)

        # --- Paso 3: clamp de posición ---
        # NOTA DE INTERPRETACION: la Sec. B.4 escribe el clamp como
        #     u <- min(u, nocional_max_posicion(S)/S - |inv|)
        # sin distinguir compra de venta. Tomado literalmente, al llegar al tope
        # de exposición ese término vale 0 y anula AMBAS componentes — incluida
        # la que REDUCE la posición. El sistema quedaría atrapado en el límite,
        # incapaz de deshacer, que es lo contrario de lo que una abrazadera de
        # riesgo debe hacer (y rompería la ruta de cierre del halt, que necesita
        # emitir exactamente esas órdenes).
        # Se aplica por tanto de forma DIRECCIONAL, sobre el inventario resultante:
        #     inv + u_c <= +tope        y        inv - u_v >= -tope
        # Con inventario largo esto reproduce la fórmula del documento para la
        # compra, y deja la venta libre, que es la intención evidente.
        tope_pos_btc = CTE.nocional_max_posicion(S) / S
        margen_compra = max(0.0, tope_pos_btc - inventario)
        margen_venta = max(0.0, tope_pos_btc + inventario)
        u_c = min(u_c, margen_compra)
        u_v = min(u_v, margen_venta)
        if u_c <= 0.0 and u_v <= 0.0:
            return 0.0, 0.0, MotivoNoEnvio.CLAMP_POSICION

        # --- Paso 4: cuantización comercial ---
        u_c_q = self.filtros.cuantizar(u_c, S)
        u_v_q = self.filtros.cuantizar(u_v, S)

        # --- Paso 5: re-validación DESPUÉS de cuantizar ---
        if u_c_q * S > self.nocional_max_orden + 1e-9 or u_c_q > margen_compra + 1e-12:
            u_c_q = 0.0
        if u_v_q * S > self.nocional_max_orden + 1e-9 or u_v_q > margen_venta + 1e-12:
            u_v_q = 0.0

        # --- Paso 6: bajo resolución ---
        if u_c_q <= 0.0 and u_v_q <= 0.0:
            self.n_bajo_resolucion += 1
            return 0.0, 0.0, MotivoNoEnvio.BAJO_RESOLUCION

        return u_c_q, u_v_q, MotivoNoEnvio.ENVIADA


# ==============================================================================
# RUTA DE CIERRE (Sec. B.5, semántica del halt)
# ==============================================================================
class ErrorDeCierre(RuntimeError):
    """El cierre no pudo confirmarse plano. NUNCA se degrada a éxito."""


class RutaDeCierre:
    """Cierra la posición y NO reporta éxito sin confirmación de posición plana.

    Es deliberadamente pesimista. El fallo que esta clase previene es el más caro
    del sistema: dar un halt por completado mientras queda exposición abierta y
    sin supervisión, porque a partir de ese momento nadie está mirando.

    `cerrar_fn(inv)` debe emitir la orden de cierre; `leer_posicion_fn()` debe
    devolver la posición SEGÚN EL EXCHANGE, no según la contabilidad local.
    """

    def __init__(
        self,
        cerrar_fn,
        leer_posicion_fn,
        alertar_fn=None,
        tolerancia: float = CTE.TOL_INV,
        max_intentos: int = 5,
        backoff_base: float = 0.5,
        dormir_fn=time.sleep,
    ):
        self.cerrar_fn = cerrar_fn
        self.leer_posicion_fn = leer_posicion_fn
        self.alertar_fn = alertar_fn
        self.tolerancia = float(tolerancia)
        self.max_intentos = int(max_intentos)
        self.backoff_base = float(backoff_base)
        self.dormir_fn = dormir_fn
        self.intentos = 0
        self.escalado = False

    def ejecutar(self) -> float:
        """Cierra y devuelve la posición final confirmada. Levanta si no queda plana."""
        self.intentos = 0
        self.escalado = False
        ultimo_error = None
        for intento in range(self.max_intentos):
            self.intentos = intento + 1
            try:
                posicion = float(self.leer_posicion_fn())
                if abs(posicion) <= self.tolerancia:
                    return posicion
                self.cerrar_fn(posicion)
                posicion = float(self.leer_posicion_fn())
                if abs(posicion) <= self.tolerancia:
                    return posicion
                ultimo_error = f"posicion residual {posicion:.6f} BTC"
            except Exception as err:
                ultimo_error = f"{type(err).__name__}: {err}"
            # Backoff exponencial, mismo criterio que el watchdog de la Sec. 8.4.2.
            if intento < self.max_intentos - 1:
                self.dormir_fn(self.backoff_base * (2**intento))

        self.escalado = True
        mensaje = (
            f"CIERRE FALLIDO tras {self.intentos} intentos ({ultimo_error}). "
            f"QUEDA EXPOSICION ABIERTA SIN SUPERVISION - intervencion manual."
        )
        if self.alertar_fn is not None:
            # La alerta va FUERA de la ruta crítica (Sec. C.5): su fallo jamás
            # debe tumbar el proceso, mismo principio que `log()`.
            try:
                self.alertar_fn(mensaje)
            except Exception:
                pass
        raise ErrorDeCierre(mensaje)


# ==============================================================================
# CUENTA DE PAPEL — solo para MODO LECTURA
# ==============================================================================
class CuentaPapel:
    """Contabilidad simulada sobre precios REALES. Exclusiva del Modo LECTURA.

    ⚠ ESTO NO ES EL `ACCOUNT_UPDATE` DEL USER DATA STREAM, y la distinción es la
    misma que la Sec. B.5 exige respetar: en TESTNET y MAINNET el equity DEBE
    venir del exchange, porque el cálculo local es precisamente lo que la guarda
    de desincronía existe para desconfiar. Aquí no hay exchange que consultar —no
    hay cuenta— así que el cálculo local es lo único que existe, y por eso esta
    clase está atada a LECTURA y no se usa en ningún otro modo.

    Qué justifica que exista: sin ella, en LECTURA el inventario sería siempre
    cero, y entonces las guardas 1, 2 y 6 nunca se ejercitarían. El criterio de
    aceptación más importante del documento —el test de disparo forzado del kill
    switch— necesita una cuenta que pueda perder dinero.

    Marca a mercado con el precio real: equity = efectivo + inventario·S.
    """

    def __init__(self, equity_inicial: float, comision: float = 4.0e-4):
        self.efectivo = float(equity_inicial)
        self.inventario = 0.0
        self.comision = float(comision)  # taker de futuros USDⓈ-M
        self.equity_inicial = float(equity_inicial)
        self.n_fills = 0

    def ejecutar(self, u_compra: float, u_venta: float, precio: float) -> float:
        """Aplica un fill al precio real. Devuelve el inventario resultante."""
        neto = float(u_compra) - float(u_venta)
        bruto = abs(float(u_compra)) + abs(float(u_venta))
        if bruto <= 0.0 or precio <= 0.0:
            return self.inventario
        self.efectivo -= neto * precio
        self.efectivo -= bruto * precio * self.comision
        self.inventario += neto
        self.n_fills += 1
        return self.inventario

    def equity(self, precio: float) -> float:
        return self.efectivo + self.inventario * float(precio)

    def perdida_ficticia(self, monto: float) -> None:
        """Inyecta una pérdida. SOLO para el test de disparo forzado (Sec. F).

        "Un freno que nunca se probó no es un freno." Este método es el mecanismo
        que permite ejercitar la cadena completa —deteccion, cierre, posicion
        plana confirmada, halt, alerta, volcado— sin esperar a que el mercado
        cause la pérdida por su cuenta.
        """
        self.efectivo -= abs(float(monto))


def formatear_estado_riesgo(capa: CapaRiesgo, S: float) -> str:
    """Línea de diagnóstico en ASCII para el arranque y el resumen de episodio."""
    equity = capa.equity_actual if capa.equity_actual is not None else float("nan")
    inicio = capa.equity_inicio if capa.equity_inicio is not None else float("nan")
    dd = equity - inicio if math.isfinite(equity) and math.isfinite(inicio) else float("nan")
    return (
        f"[RIESGO] equity={equity:.2f} (inicio {inicio:.2f}, dd {dd:+.2f} contra cap "
        f"{capa.perdida_max:.0f}) inv={capa.inventario_local:.6f} BTC "
        f"expuesto={abs(capa.inventario_local)*S:.0f}/{CTE.nocional_max_posicion(S):.0f} USD "
        f"env={capa.n_enviadas} rech={capa.n_rechazadas} "
        f"bajo_resolucion={capa.n_bajo_resolucion} "
        f"halt={capa.causa_halt.name}"
    )
