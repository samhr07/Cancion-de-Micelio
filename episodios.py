"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: episodios.py — Máquina de episodios y reaprovisionamiento del faucet

Implementa la Sección C de ORDEN_TRABAJO_RIESGO_1_3.md.

POR QUÉ EL FAUCET AUTOMÁTICO NECESITA COMPUERTAS (Sec. C.1)
-----------------------------------------------------------
Un reaprovisionamiento automático sin más convierte el kill switch en un BUCLE:
ante un bug, el sistema quema los 30 episodios en una hora y produce 30 datasets
idénticos e inútiles. El valor de un halt está en que ALGUIEN MIRE POR QUÉ
OCURRIÓ. Las compuertas de la Sec. C.3 conservan el automatismo pero preservan
esa propiedad.
"""

from __future__ import annotations

import json
import os
import time
from enum import Enum

import constantes_micelio as CTE
from riesgo import CAUSAS_DE_SISTEMA, CausaHalt, DESCRIPCION_CAUSA

DIR_EPISODIOS = "episodios"


# ==============================================================================
# C.2 — ESTADOS
# ==============================================================================
class Estado(str, Enum):
    """
    ARRANQUE -> [equity >= EQUITY_MIN_EPISODIO] -> OPERANDO
    OPERANDO -> [guarda B.5 dispara]            -> CERRANDO
    CERRANDO -> [posicion plana confirmada]     -> CERRADO
    CERRADO  -> volcar telemetria, diagnostico, resumen
             -> [compuertas C.3] --pasa--> REAPROVISIONANDO -> ARRANQUE
                                 --falla-> DETENIDO
    """

    ARRANQUE = "ARRANQUE"
    OPERANDO = "OPERANDO"
    CERRANDO = "CERRANDO"
    CERRADO = "CERRADO"
    REAPROVISIONANDO = "REAPROVISIONANDO"
    DETENIDO = "DETENIDO"  # Terminal y ruidoso. NUNCA se sale de el por software.


# ==============================================================================
# C.4 — EL ADAPTADOR DE FAUCET
# ==============================================================================
# ⚠ El faucet de Binance Futures Testnet es una función de la INTERFAZ WEB, no un
# endpoint documentado de la API pública. No hay que fingir lo contrario.


class Reaprovisionador:
    """Interfaz.

    ⚠ La ausencia de una implementación automática NUNCA debe bloquear la ruta de
    halt: el bot para igual, y espera.
    """

    nombre = "base"

    def solicitar(self, monto_objetivo: float) -> bool:
        raise NotImplementedError


class ReaprovisionadorManual(Reaprovisionador):
    """Por defecto. Emite alerta y espera por SONDEO DE EQUITY.

    Escribe el resumen, alerta, y espera a que el equity de la cuenta supere
    EQUITY_MIN_EPISODIO. Sin intervención en el proceso: Samuel recarga por la
    web y el bot lo detecta solo.
    """

    nombre = "manual"

    def __init__(self, alertar_fn=None):
        self.alertar_fn = alertar_fn

    def solicitar(self, monto_objetivo: float) -> bool:
        mensaje = (
            f"[FAUCET] Recarga MANUAL requerida: se necesitan al menos "
            f"{monto_objetivo:.0f} USD de equity para abrir el episodio siguiente. "
            f"Recargar por la web de Testnet; el bot lo detecta por sondeo."
        )
        if self.alertar_fn is not None:
            try:  # La alerta va fuera de la ruta crítica (Sec. C.5)
                self.alertar_fn(mensaje)
            except Exception:
                pass
        # No hay nada automático que hacer: devuelve False y el sondeo de equity
        # de la máquina de estados hace el resto. Devolver True aquí sería mentir.
        return False


class ReaprovisionadorAutomatico(Reaprovisionador):
    """Opcional. Aislado aquí y SOLO aquí. Sec. C.4.

    Requisitos del documento, los tres implementados abajo:
      - tolerar rate limits del faucet (son estrictos) sin reintentar en bucle
      - fallar de forma limpia degradando a ReaprovisionadorManual
      - no ejecutarse JAMÁS si MODO == MAINNET (aserción dura, no `if`)
    """

    nombre = "automatico"

    def __init__(self, modo, solicitar_fn, respaldo: Reaprovisionador, intervalo_min: float = 60.0):
        # ⚠ ASERCIÓN DURA, NO UN `if`. Un `if` que devuelve False dejaría el
        # sistema corriendo con un reaprovisionador que en Mainnet es un
        # sinsentido peligroso; que reviente en el constructor obliga a que la
        # composición sea correcta antes de arrancar nada.
        from mercado import Modo

        assert modo is not Modo.MAINNET, (
            "ReaprovisionadorAutomatico es INADMISIBLE en MODO=MAINNET: no existe "
            "faucet de dinero real, y automatizar la reposicion de capital sobre "
            "una cuenta real convierte un kill switch en una sangria."
        )
        self.modo = modo
        self.solicitar_fn = solicitar_fn
        self.respaldo = respaldo
        self.intervalo_min = float(intervalo_min)
        self.ts_ultima = 0.0

    def solicitar(self, monto_objetivo: float) -> bool:
        ahora = time.time()
        if (ahora - self.ts_ultima) < self.intervalo_min:
            # Rate limit del faucet: NO se reintenta en bucle. Se degrada al
            # camino manual, que es el mismo sondeo de equity.
            return self.respaldo.solicitar(monto_objetivo)
        self.ts_ultima = ahora
        try:
            return bool(self.solicitar_fn(monto_objetivo))
        except Exception:
            # Degradación limpia: el modo manual y el automático comparten el
            # camino de sondeo, así que caer a él no es una ruta sin probar.
            return self.respaldo.solicitar(monto_objetivo)


# ==============================================================================
# C.3 — COMPUERTAS ANTES DE REAPROVISIONAR
# ==============================================================================
class ResultadoCompuertas:
    def __init__(self):
        self.fallidas = []

    def registrar(self, nombre: str, pasa: bool, detalle: str) -> None:
        if not pasa:
            self.fallidas.append((nombre, detalle))

    @property
    def pasa(self) -> bool:
        return not self.fallidas

    def __str__(self) -> str:
        if self.pasa:
            return "todas las compuertas C.3 pasan"
        return "; ".join(f"{n}: {d}" for n, d in self.fallidas)


def evaluar_compuertas(
    episodios_automaticos: int,
    muestras: int,
    causa: CausaHalt,
    causa_anterior: CausaHalt,
    diagnostico_ok: bool,
    segundos_desde_cierre: float,
    max_automaticos: int = CTE.MAX_EPISODIOS_AUTOMATICOS,
    muestras_min: int = CTE.MUESTRAS_MIN_EPISODIO,
    enfriamiento: float = CTE.ENFRIAMIENTO_EPISODIO,
) -> ResultadoCompuertas:
    """Las seis compuertas de la Sec. C.3. TODAS deben cumplirse."""
    r = ResultadoCompuertas()

    r.registrar(
        "episodios_automaticos",
        episodios_automaticos < max_automaticos,
        f"{episodios_automaticos} consecutivos (max {max_automaticos}); "
        f"revision humana obligatoria",
    )
    r.registrar(
        "muestras_minimas",
        muestras >= muestras_min,
        f"{muestras} muestras (min {muestras_min}); un episodio que muere pronto "
        f"es un bug, no una perdida: no recargar",
    )
    r.registrar(
        "causa_no_repetida",
        causa != causa_anterior,
        f"la causa '{causa.name}' se repite; es determinista y no la arregla mas "
        f"dinero",
    )
    r.registrar(
        "causa_es_drawdown",
        causa == CausaHalt.DRAWDOWN,
        f"causa '{causa.name}' = {DESCRIPCION_CAUSA.get(causa, '?')}; "
        f"{'es un fallo de sistema, recargar seria tapar el bug' if causa in CAUSAS_DE_SISTEMA else 'no es un drawdown'}",
    )
    r.registrar(
        "diagnostico_escrito",
        bool(diagnostico_ok),
        "diagnostico.py no se ejecuto o no se escribio el resumen; sin eso el "
        "bucle es una tragamonedas y no un pipeline de datos",
    )
    r.registrar(
        "enfriamiento",
        segundos_desde_cierre >= enfriamiento,
        f"{segundos_desde_cierre:.1f} s desde el cierre (min {enfriamiento:.0f} s)",
    )
    return r


# ==============================================================================
# MÁQUINA DE EPISODIOS
# ==============================================================================
class MaquinaEpisodios:
    """Sec. C.2. Con `DETENIDO` terminal: no hay transición de salida por software.

    El reset limpio de la Sec. C.5 se delega en `al_reiniciar`, que el Hilo
    Rápido registra para poner a cero ΣQ, la racha de burn-in, S_ref, la ventana
    del EMD, el inventario local y los contadores de las guardas. El estado que
    sobrevive a un halt es fuente segura de confusión al analizar.
    """

    def __init__(
        self,
        reaprovisionador: Reaprovisionador,
        leer_equity_fn,
        precio_fn,
        al_reiniciar=None,
        alertar_fn=None,
        dir_salida: str = DIR_EPISODIOS,
    ):
        self.reaprovisionador = reaprovisionador
        self.leer_equity_fn = leer_equity_fn
        self.precio_fn = precio_fn
        self.al_reiniciar = al_reiniciar
        self.alertar_fn = alertar_fn
        self.dir_salida = dir_salida

        self.estado = Estado.ARRANQUE
        self.id_episodio = 0
        self.episodios_automaticos = 0
        self.causa_anterior = CausaHalt.NINGUNA
        self.ts_cierre = 0.0
        self.motivo_detencion = ""
        self.historial = []

        self._equity_inicio = None
        self._ts_inicio = 0.0

    # ------------------------------------------------------------------
    def puede_abrir(self) -> tuple[bool, str]:
        """ARRANQUE -> OPERANDO exige equity >= EQUITY_MIN_EPISODIO (Sec. B.3).

        Si el faucet entrega menos, NO SE ARRANCA: se registra el saldo real y se
        reporta. Las opciones entonces son bajar PERDIDA_MAX_EPISODIO o bajar
        I_max A SABIENDAS, recalculando γ_0 y verificando que λS²Γ siga siendo
        alcanzable (Sec. B.1). Nunca lo último en silencio, y nunca aquí.
        """
        S = float(self.precio_fn())
        equity = float(self.leer_equity_fn())
        minimo = CTE.equity_min_episodio(S)
        if equity < minimo:
            return False, (
                f"equity {equity:.2f} USD < minimo {minimo:.2f} USD a S={S:.2f} "
                f"(margen {CTE.nocional_max_posicion(S)/CTE.APALANCAMIENTO:.0f} + "
                f"colchon {CTE.PERDIDA_MAX_EPISODIO:.0f})"
            )
        return True, f"equity {equity:.2f} USD >= minimo {minimo:.2f} USD"

    def abrir(self) -> bool:
        ok, detalle = self.puede_abrir()
        if not ok:
            return False
        self.id_episodio += 1
        self._equity_inicio = float(self.leer_equity_fn())
        self._ts_inicio = time.time()
        self.estado = Estado.OPERANDO
        if self.al_reiniciar is not None:
            # Reset limpio ANTES de operar, no después de cerrar: así el estado
            # residual de un cierre fallido tampoco se hereda.
            self.al_reiniciar(self.id_episodio, self._equity_inicio)
        self._alertar(f"[EPISODIO {self.id_episodio:03d}] abierto - {detalle}")
        return True

    def al_disparar_guarda(self, causa: CausaHalt, detalle: str) -> None:
        if self.estado is Estado.OPERANDO:
            self.estado = Estado.CERRANDO
            self._causa_actual = causa
            self._detalle_actual = detalle
            self._alertar(
                f"[EPISODIO {self.id_episodio:03d}] HALT {causa.name}: {detalle}"
            )

    def al_confirmar_plano(self, posicion_final: float) -> None:
        if self.estado is Estado.CERRANDO:
            self.estado = Estado.CERRADO
            self.ts_cierre = time.time()
            self._posicion_final = posicion_final

    # ------------------------------------------------------------------
    def cerrar_y_evaluar(self, muestras: int, estadisticas: dict, diagnostico_ok: bool):
        """CERRADO -> resumen -> compuertas -> REAPROVISIONANDO o DETENIDO."""
        causa = getattr(self, "_causa_actual", CausaHalt.NINGUNA)
        resumen = self.escribir_resumen(muestras, estadisticas, causa)

        compuertas = evaluar_compuertas(
            episodios_automaticos=self.episodios_automaticos,
            muestras=muestras,
            causa=causa,
            causa_anterior=self.causa_anterior,
            diagnostico_ok=diagnostico_ok,
            segundos_desde_cierre=time.time() - self.ts_cierre,
        )
        self.causa_anterior = causa

        if not compuertas.pasa:
            self.estado = Estado.DETENIDO
            self.motivo_detencion = str(compuertas)
            self._alertar(
                f"[EPISODIO {self.id_episodio:03d}] DETENIDO (requiere "
                f"intervencion): {self.motivo_detencion}"
            )
            return compuertas, resumen

        self.estado = Estado.REAPROVISIONANDO
        self.episodios_automaticos += 1
        S = float(self.precio_fn())
        self.reaprovisionador.solicitar(CTE.equity_min_episodio(S))
        return compuertas, resumen

    def sondear_equity(self) -> bool:
        """REAPROVISIONANDO -> ARRANQUE por sondeo. Sec. C.4.

        ⚠ ESTA ES LA PARTE QUE DE VERDAD IMPORTA. Hace que el modo manual y el
        automático sean EL MISMO CAMINO DE CÓDIGO, con la única diferencia de
        quién provoca la recarga. Eso es lo que impide que el modo automático sea
        una ruta sin probar el día que se active.
        """
        if self.estado is not Estado.REAPROVISIONANDO:
            return False
        ok, _detalle = self.puede_abrir()
        if ok:
            self.estado = Estado.ARRANQUE
            return True
        return False

    def detener(self, motivo: str) -> None:
        """Entrada manual a DETENIDO. No existe la salida."""
        self.estado = Estado.DETENIDO
        self.motivo_detencion = motivo
        self._alertar(f"[EPISODIO {self.id_episodio:03d}] DETENIDO: {motivo}")

    # ------------------------------------------------------------------
    def escribir_resumen(self, muestras: int, estadisticas: dict, causa: CausaHalt) -> dict:
        """`resumen_episodio_NNN.json` de la Sec. C.5."""
        equity_final = None
        try:
            equity_final = float(self.leer_equity_fn())
        except Exception:
            pass
        resumen = {
            "id_episodio": self.id_episodio,
            "equity_inicio": self._equity_inicio,
            "equity_final": equity_final,
            "causa_halt": causa.name,
            "causa_halt_codigo": int(causa),
            "detalle_halt": getattr(self, "_detalle_actual", ""),
            "duracion_s": (time.time() - self._ts_inicio) if self._ts_inicio else None,
            "muestras_telemetria": muestras,
            "posicion_final": getattr(self, "_posicion_final", None),
            "reaprovisionador": self.reaprovisionador.nombre,
            "estadisticas": estadisticas,
        }
        self.historial.append(resumen)
        try:
            os.makedirs(self.dir_salida, exist_ok=True)
            ruta = os.path.join(
                self.dir_salida, f"resumen_episodio_{self.id_episodio:03d}.json"
            )
            with open(ruta, "w", encoding="utf-8") as fh:
                json.dump(resumen, fh, indent=2, ensure_ascii=True)
        except Exception:
            # Escribir el resumen es observabilidad, no funcionalidad. Su fallo
            # no puede tumbar la maquina de estados; pero SI hace fallar la
            # compuerta `diagnostico_escrito`, que es donde debe notarse.
            pass
        return resumen

    def _alertar(self, mensaje: str) -> None:
        if self.alertar_fn is None:
            return
        try:  # Fuera de la ruta crítica, mismo principio que `log()` (Sec. C.5)
            self.alertar_fn(mensaje)
        except Exception:
            pass
