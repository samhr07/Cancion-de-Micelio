"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Repositorio: samhr07/Cancion-de-Micelio
Módulo: main.py (Orquestador Principal)

Descripción:
Inicializa la topología de procesos aislados, asigna la memoria compartida a nivel de SO,
configura el vector de hiperparámetros dinámicos (Hot-Reloading) y gestiona el apagado seguro.

==============================================================================
VARIABLES DINÁMICAS DE MODIFICACIÓN (HOT-RELOADING)
Buscar con Ctrl+F usando los siguientes identificadores numéricos:

[Filtro de Kalman]
0: PARAM_R_S_BASE      -> Varianza mínima del tick size
1: PARAM_R_EMD         -> Varianza histórica de Hilbert
2: PARAM_BETA_JITTER   -> Factor de sensibilidad al retardo de red

[NMPC - Optimizador CasADi]
3: PARAM_Q_DELTA       -> Rigidez de seguimiento de cobertura (Delta)
4: PARAM_Q_BASE_INV    -> Costo asintótico por mantenimiento de inventario
5: PARAM_MU_INV        -> Hiperparámetro de sensibilidad al riesgo de exposición
6: PARAM_R_BASE        -> Costo estático base por comisiones (Taker fee)
7: PARAM_KAPPA_RISK    -> Aversión a la inestabilidad de volumen del Micelio

[Micelio - Simulación y Loeper]
8: PARAM_ETA_IMPACT    -> Sensibilidad estructural al impacto
9: PARAM_THETA_OU      -> Tasa de reversión a la media (Ornstein-Uhlenbeck)
10: PARAM_SIGMA_OU     -> Volatilidad intrínseca del spread simulado
==============================================================================
"""

import multiprocessing as mp
from multiprocessing import shared_memory
import signal
import sys
import time
import numpy as np

# ==============================================================================
# 1. CONFIGURACIÓN BASE Y MAPA DE PARÁMETROS DINÁMICOS
# ==============================================================================

# Variables estáticas del sistema (No se alteran en caliente)
SYMBOL = "BTCUSDT"
IS_TESTNET = True
N_HORIZON = 20
W_MIN_EMD = 100
W_MAX_EMD = 1000

# Tamaño del vector de hiperparámetros dinámicos
TOTAL_PARAMS = 11


def initialize_default_parameters() -> np.ndarray:
    """Genera el vector inicial con los hiperparámetros de control por defecto."""
    params = np.zeros(TOTAL_PARAMS, dtype=np.float64)

    # [Filtro de Kalman]
    params[0] = 0.5  # PARAM_R_S_BASE
    params[1] = 1.0  # PARAM_R_EMD
    params[2] = 0.02  # PARAM_BETA_JITTER

    # [NMPC]
    params[3] = 10.0  # PARAM_Q_DELTA
    params[4] = 2.0  # PARAM_Q_BASE_INV
    params[5] = 5.0  # PARAM_MU_INV
    params[6] = 0.05  # PARAM_R_BASE
    params[7] = 8.0  # PARAM_KAPPA_RISK

    # [Micelio]
    params[8] = 0.15  # PARAM_ETA_IMPACT
    params[9] = 0.5  # PARAM_THETA_OU
    params[10] = 0.2  # PARAM_SIGMA_OU

    return params


# ==============================================================================
# 2. DEFINICIONES DE GESTIÓN DE MEMORIA COMPARTIDA (IPC LATENCIA CERO)
# ==============================================================================


def allocate_shared_memory(name: str, size: int) -> shared_memory.SharedMemory:
    """
    Crea un bloque de memoria compartida. Si ya existe (cierre abrupto previo),
    lo vincula y lo sobreescribe para evitar errores.
    """
    try:
        shm = shared_memory.SharedMemory(name=name, create=True, size=size)
        print(f"[OK] Memoria asignada: {name} ({size} bytes)")
    except FileExistsError:
        print(f"[WARN] Memoria {name} ya existía. Reasignando y purgando zombis...")
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
    shm_env_name: str, shm_actuator_name: str, shm_params_name: str
):
    """
    PROCESO 1: I/O Bound. WebSocket, Token Bucket, Cuantización y Watchdog.
    """
    print("[INIT] Motor de Red iniciado.")
    # TODO: Lógica de reconexión infinita asíncrona y bypass de Hot-Reloading
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        pass


def fast_thread_process(
    shm_env_name: str,
    shm_micelio_name: str,
    shm_actuator_name: str,
    shm_params_name: str,
):
    """
    PROCESO 2: CPU/GPU Bound. EAKF, Loeper (GPU) y NMPC (CasADi).
    """
    print("[INIT] Hilo Rápido iniciado (NMPC + EAKF).")
    # Conexión al vector de parámetros en caliente
    existing_shm = shared_memory.SharedMemory(name=shm_params_name)
    params = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=existing_shm.buf)

    try:
        while True:
            # Ejemplo de lectura dinámica: q_delta_actual = params[3]
            time.sleep(0.01)
    except KeyboardInterrupt:
        pass
    finally:
        existing_shm.close()


def slow_thread_process(shm_env_name: str, shm_micelio_name: str, shm_params_name: str):
    """
    PROCESO 3: CPU Bound. EMD, Hilbert y simulación OU.
    """
    print("[INIT] Hilo Lento iniciado (Micelio + HHT).")
    existing_shm = shared_memory.SharedMemory(name=shm_params_name)
    params = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=existing_shm.buf)

    try:
        while True:
            time.sleep(0.1)
    except KeyboardInterrupt:
        pass
    finally:
        existing_shm.close()


# ==============================================================================
# 4. ORQUESTADOR PRINCIPAL Y GRACEFUL SHUTDOWN
# ==============================================================================


def main():
    # 4.1. Cálculo dimensional de memoria (Bytes)
    BYTES_ENV = 128  # [P_spot, V, Dropout_Flag]
    BYTES_MICELIO = 256  # [Omega, w_m, lambda_sim, R_n]
    BYTES_ACTUATOR = 64  # [u_compra, u_venta]
    BYTES_PARAMS = TOTAL_PARAMS * np.dtype(np.float64).itemsize  # Bloque Hot-Reloading

    shm_env = None
    shm_micelio = None
    shm_actuator = None
    shm_params = None
    processes = []

    def signal_handler(sig, frame):
        """Atrapa Ctrl+C o comandos del SO para matar el bot suavemente."""
        print("\n[ALERTA] Señal de apagado recibida. Iniciando Graceful Shutdown...")

        for p in processes:
            if p.is_alive():
                p.terminate()
                p.join()

        shm_list = [
            shm
            for shm in [shm_env, shm_micelio, shm_actuator, shm_params]
            if shm is not None
        ]
        cleanup_shared_memory(shm_list)

        print("[EXIT] Sistema apagado correctamente. Natura Recusat Silentium.")
        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    try:
        print("=== INICIALIZANDO CANCIÓN DEL MICELIO ===")
        # 4.2. Reserva de Memoria Compartida
        shm_env = allocate_shared_memory("shm_entorno", BYTES_ENV)
        shm_micelio = allocate_shared_memory("shm_estructural", BYTES_MICELIO)
        shm_actuator = allocate_shared_memory("shm_control", BYTES_ACTUATOR)
        shm_params = allocate_shared_memory("shm_parametros", BYTES_PARAMS)

        # 4.3. Inicialización del vector de hiperparámetros
        params_array = np.ndarray(
            (TOTAL_PARAMS,), dtype=np.float64, buffer=shm_params.buf
        )
        params_array[:] = initialize_default_parameters()[:]

        # 4.4. Instanciación y arranque de Procesos
        p_net = mp.Process(
            target=network_engine_process,
            args=(shm_env.name, shm_actuator.name, shm_params.name),
        )
        p_fast = mp.Process(
            target=fast_thread_process,
            args=(shm_env.name, shm_micelio.name, shm_actuator.name, shm_params.name),
        )
        p_slow = mp.Process(
            target=slow_thread_process,
            args=(shm_env.name, shm_micelio.name, shm_params.name),
        )

        processes = [p_net, p_fast, p_slow]

        for p in processes:
            p.start()

        # El proceso padre se queda dormido supervisando
        while True:
            time.sleep(1)

    except Exception as e:
        print(f"[ERROR FATAL] {e}")
        signal_handler(None, None)


# ==============================================================================
# PUNTO DE ENTRADA (MÉTODO SPAWN OBLIGATORIO PARA CUDA)
# ==============================================================================
if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    main()
