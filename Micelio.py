"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Repositorio: samhr07/Cancion-de-Micelio
Módulo: main.py (Orquestador Principal)

Descripción:
Inicializa la topología de procesos aislados (Red, Hilo Rápido, Hilo Lento),
asigna la memoria compartida a nivel de SO y gestiona el apagado seguro.
"""

import multiprocessing as mp
from multiprocessing import shared_memory
import signal
import sys
import time
import numpy as np
from dataclasses import dataclass, field

# ==============================================================================
# 1. ESTRUCTURAS DE CONFIGURACIÓN (DATACLASSES)
# ==============================================================================


@dataclass
class KalmanConfig:
    """Parámetros de inicialización del EAKF."""

    r_S_base: float = 0.5  # Varianza mínima del tick size
    r_EMD: float = 1.0  # Varianza histórica de Hilbert
    q_base: list = field(
        default_factory=lambda: [0.01, 0.05, 0.1]
    )  # Varianzas térmicas base
    beta_jitter: float = 0.02  # Factor de sensibilidad al retardo de red


@dataclass
class NMPCConfig:
    """Matrices de ponderación cuadrática del optimizador CasADi/acados."""

    N_horizon: int = 20  # Nodos de predicción
    q_delta: float = 10.0  # Rigidez de seguimiento de cobertura
    q_base_inv: float = 2.0  # Costo asintótico por mantenimiento de inventario
    mu_inv: float = 5.0  # Hiperparámetro de sensibilidad al riesgo (Inventario)
    r_base: float = 0.05  # Costo estático por comisiones (Taker fee)
    kappa_risk: float = 8.0  # Aversión a la inestabilidad del Micelio


@dataclass
class MicelioConfig:
    """Parámetros para la EDP de Loeper y EMD."""

    eta_impact: float = 0.15  # Sensibilidad estructural al impacto
    theta_ou: float = 0.5  # Tasa de reversión a la media (Ornstein-Uhlenbeck)
    sigma_ou: float = 0.2  # Volatilidad del spread simulado
    w_min: int = 100  # Ventana mínima EMD
    w_max: int = 1000  # Ventana máxima EMD


@dataclass
class BotConfig:
    """Agrupación maestra de configuración."""

    symbol: str = "BTCUSDT"
    is_testnet: bool = True
    kalman: KalmanConfig = field(default_factory=KalmanConfig)
    nmpc: NMPCConfig = field(default_factory=NMPCConfig)
    micelio: MicelioConfig = field(default_factory=MicelioConfig)


# ==============================================================================
# 2. DEFINICIONES DE GESTIÓN DE MEMORIA COMPARTIDA (IPC LATENCIA CERO)
# ==============================================================================


def allocate_shared_memory(name: str, size: int) -> shared_memory.SharedMemory:
    """
    Crea un bloque de memoria compartida. Si ya existe (por un cierre abrupto previo),
    lo vincula y lo sobreescribe para evitar el error 'File exists'.
    """
    try:
        shm = shared_memory.SharedMemory(name=name, create=True, size=size)
        print(f"[OK] Memoria asignada: {name} ({size} bytes)")
    except FileExistsError:
        print(f"[WARN] Memoria {name} ya existía. Reasignando (limpiando zombis)...")
        shm = shared_memory.SharedMemory(name=name, create=False)
        shm.unlink()  # Purga inmediata
        shm = shared_memory.SharedMemory(name=name, create=True, size=size)
    return shm


def cleanup_shared_memory(shm_list: list):
    """Cierra y purga los bloques de RAM al detener el sistema."""
    for shm in shm_list:
        try:
            shm.close()
            shm.unlink()
            print(f"[OK] Memoria liberada: {shm.name}")
        except Exception as e:
            print(f"[ERROR] Fallo al liberar memoria {shm.name}: {e}")


# ==============================================================================
# 3. ESQUELETOS DE PROCESOS AISLADOS (EVASIÓN DEL GIL)
# ==============================================================================


def network_engine_process(
    config: BotConfig, shm_env_name: str, shm_actuator_name: str
):
    """
    PROCESO 1: I/O Bound.
    Responsable del WebSocket, Token Bucket, Cuantización y Watchdog.
    NO contiene cálculos matemáticos pesados.
    """
    print("[INIT] Motor de Red iniciado.")
    # TODO: Conectar a shared_memory mediante shm_env_name
    # TODO: Iniciar Event Loop de asyncio
    # TODO: Implementar WebSocket con Binance y Rate Limiter

    try:
        while True:
            time.sleep(1)  # Simulación de loop asíncrono
    except KeyboardInterrupt:
        print("[SHUTDOWN] Motor de Red detenido.")


def fast_thread_process(
    config: BotConfig, shm_env_name: str, shm_micelio_name: str, shm_actuator_name: str
):
    """
    PROCESO 2: CPU/GPU Bound.
    Núcleo táctico. Ejecuta EAKF, interrumpe a la GPU (Loeper) y resuelve NMPC (CasADi).
    """
    print("[INIT] Hilo Rápido iniciado (NMPC + EAKF).")
    # TODO: Conectar a todas las memorias compartidas
    # TODO: Configurar JAX/XLA para las matrices de Kalman
    # TODO: Inicializar solver acados/CasADi

    try:
        while True:
            time.sleep(0.01)  # Simulación de interrupción de microsegundos
    except KeyboardInterrupt:
        print("[SHUTDOWN] Hilo Rápido detenido.")


def slow_thread_process(config: BotConfig, shm_env_name: str, shm_micelio_name: str):
    """
    PROCESO 3: CPU Bound.
    Análisis microestructural. Ejecuta EMD, Hilbert y simulación OU para testnet.
    """
    print("[INIT] Hilo Lento iniciado (Micelio + HHT).")
    # TODO: Inicializar buffers para EMD

    try:
        while True:
            time.sleep(0.1)  # Simulación de cálculo denso
    except KeyboardInterrupt:
        print("[SHUTDOWN] Hilo Lento detenido.")


# ==============================================================================
# 4. ORQUESTADOR PRINCIPAL Y GRACEFUL SHUTDOWN
# ==============================================================================


def main():
    # 4.1. Configuración maestra
    config = BotConfig()

    # 4.2. Asignación de tamaños de memoria (Bytes)
    # (En producción, calcular con np.dtype.itemsize * dimensiones)
    BYTES_ENV = 128  # [P_spot, V, Dropout_Flag]
    BYTES_MICELIO = 256  # [Omega, w_m, lambda_sim, R_n]
    BYTES_ACTUATOR = 64  # [u_compra, u_venta]

    shm_env = None
    shm_micelio = None
    shm_actuator = None
    processes = []

    def signal_handler(sig, frame):
        """Atrapa Ctrl+C o comandos del SO para matar el bot suavemente."""
        print("\n[ALERTA] Señal de apagado recibida. Iniciando Graceful Shutdown...")

        # 1. Terminar procesos hijos
        for p in processes:
            if p.is_alive():
                p.terminate()
                p.join()

        # 2. Purgar memoria compartida
        shm_list = [
            shm for shm in [shm_env, shm_micelio, shm_actuator] if shm is not None
        ]
        cleanup_shared_memory(shm_list)

        print("[EXIT] Sistema apagado correctamente. Natura Recusat Silentium.")
        sys.exit(0)

    # Enganchar el handler
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    try:
        # 4.3. Reserva de Memoria Compartida
        print("=== INICIALIZANDO CANCIÓN DEL MICELIO ===")
        shm_env = allocate_shared_memory("shm_entorno", BYTES_ENV)
        shm_micelio = allocate_shared_memory("shm_estructural", BYTES_MICELIO)
        shm_actuator = allocate_shared_memory("shm_control", BYTES_ACTUATOR)

        # 4.4. Instanciación de Procesos
        p_net = mp.Process(
            target=network_engine_process,
            args=(config, shm_env.name, shm_actuator.name),
        )
        p_fast = mp.Process(
            target=fast_thread_process,
            args=(config, shm_env.name, shm_micelio.name, shm_actuator.name),
        )
        p_slow = mp.Process(
            target=slow_thread_process, args=(config, shm_env.name, shm_micelio.name)
        )

        processes = [p_net, p_fast, p_slow]

        # 4.5. Arranque concurrente
        for p in processes:
            p.start()

        # El proceso padre se queda dormido supervisando
        while True:
            time.sleep(1)

    except Exception as e:
        print(f"[ERROR FATAL] {e}")
        signal_handler(None, None)


# ==============================================================================
# PUNTO DE ENTRADA (MÉTODO SPAWN OBLIGATORIO)
# ==============================================================================
if __name__ == "__main__":
    # Blindaje contra colapso del contexto de CUDA al bifurcar procesos
    mp.set_start_method("spawn", force=True)
    main()
