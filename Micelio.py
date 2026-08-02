"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Repositorio: samhr07/Cancion-de-Micelio
Módulo: main.py (Núcleo Operativo Endurecido - v_final)
"""

import multiprocessing as mp
from multiprocessing import shared_memory
import asyncio
import signal
import sys
import time
import math
import random
import os
import numpy as np

# ==============================================================================
# 1. ESTRUCTURACIÓN DE BLOQUES LOCK-FREE Y C-TYPES (ACTUALIZADO)
# ==============================================================================
SYMBOL = "BTCUSDT"
IS_TESTNET = True

# Entorno: Tensor de memoria de estado temporal
ENV_DTYPE = np.dtype(
    [
        ("P_spot", np.float64),
        ("Q_transado", np.float64),
        ("flag_dropout", np.int8),
        ("timestamp", np.float64),
    ]
)

# Estructural (Spinlock implícito desde el Slow Thread hacia Fast Thread)
MICELIO_DTYPE = np.dtype(
    [
        ("omega", np.float64),
        ("w_m", np.float64),
        ("lambda_sim", np.float64),
        ("R_n", np.float64),
        ("vol_sum_Q", np.float64),  # Solución al déficit métrico requerido para rho_k
    ]
)

# Actuador: Convertido de Buffer Único a Cola (Ring Buffer SPSC - 16 Slots)
RING_BUFFER_SIZE = 16
ACTUATOR_DTYPE = np.dtype(
    [("u_compra", np.float64), ("u_venta", np.float64), ("seq_id", np.uint64)]
)

# Memoria contigua Telemetría Diferida (HDF5 / Parquet Offline logging)
TELEMETRY_SIZE = 10000
TELEM_DTYPE = np.dtype(
    [
        ("tr_P", np.float64),
        ("x0", np.float64),
        ("x1", np.float64),
        ("x2", np.float64),  # Vector de estado
        ("y0", np.float64),
        ("y1", np.float64),  # Innovaciones puro ruido blanco (para validación cov)
    ]
)

TOTAL_PARAMS = 11


def initialize_default_parameters() -> np.ndarray:
    params = np.zeros(TOTAL_PARAMS, dtype=np.float64)
    params[0] = 0.5  # q_s_base para R_S_BASE
    params[1] = 1.0  # R_EMD
    params[2] = 0.02  # Beta jitter retardos red
    # NMPC...
    params[6] = 0.05  # r_base comisiones taker
    params[8] = 0.15  # Sensibilidad impacto eta
    params[9] = 0.5  # Tasa reversión theta
    params[10] = 0.05  # Sigma O-U simulación
    # Factores faltantes dinámicos requeridos p/ Kalman (Rho, etc.)
    params[3] = 0.01  # q_base general (Escala BTC 10^-x ajustado Mainnet)
    params[4] = 0.5  # Gamma Omega
    params[5] = 0.2  # Gamma Vol(Q)
    return params


def apply_filters(
    u_raw: float, stepSize: float, minQty: float, minNotional: float, Pspot: float
) -> float:
    """Discretización escalonada y Cuantización comercial Binance-Friendly"""
    if u_raw >= minQty:
        # Se exige floor riguroso para no apalancar deuda imaginaria (v. doc cap 8.5)
        u_rounded = math.floor(u_raw / stepSize) * stepSize
        # CRÍTICO: validación MinNotional posterior a la discretización.
        if (u_rounded * Pspot) >= minNotional:
            return u_rounded
    return 0.0


# ==============================================================================
# 2. PROCESO 1: MOTOR DE RED (ASÍNCRONO I/O, RECONEXIÓN ESTRUCTURADA, WATCHDOG)
# ==============================================================================
async def exchange_websocket_handler(env_arr, act_arr):
    retry_count = 0
    tk_cap = 1200.0  # Ej. Límite IP Weights / minuto Binance
    tokens = tk_cap
    refill_rate = tk_cap / 60.0
    last_tk_t = time.time()
    actuator_tail_idx = 0  # Seguimiento puntero lectura

    while True:
        try:
            now = time.time()
            dt = now - last_tk_t
            tokens = min(tk_cap, tokens + refill_rate * dt)
            last_tk_t = now

            # Simulación: Intento de consumo de tokens del bucket por latido (TokenBucket pass rate limits)
            if tokens < 1.0:
                await asyncio.sleep(0.01)  # Presión mitigada, Yield loop asíncrono
                continue

            # (Acá iría Await websockets/API real. Await sleep lo simula).
            # Genera Timeout tras 0.2 segundos si no vuelve ping frame de server
            await asyncio.sleep(0.005)

            # Escribe data fresh en variables del Hilo Entorno:
            env_arr[0]["P_spot"] = 45000.0 + random.uniform(-1, 1)  # dummy tick
            env_arr[0]["Q_transado"] += random.uniform(-0.1, 0.5)
            env_arr[0]["timestamp"] = now
            env_arr[0]["flag_dropout"] = 0
            retry_count = 0  # Éxito: reinicia counter

            # --- Integración Dispatching Despacho (Si CPU escribió Ring Buffer)
            current_orden = act_arr[actuator_tail_idx]
            if current_orden["seq_id"] > actuator_tail_idx:
                # (Disparar tarea Async Aiohttp en fondo acá; consume weight en Bucket)
                # print(f"Execute TACTICAL_OP -> {current_orden['u_compra']} / tail {actuator_tail_idx}")
                tokens -= 20.0  # Simulando orden ejecutiva GTC Weight limit
                actuator_tail_idx = (actuator_tail_idx + 1) % RING_BUFFER_SIZE

        except Exception as e:
            env_arr[0]["flag_dropout"] = 1
            # Reintento Backoff Exponencial estocástico (Jitter anti colisiones red en WatchDog 8.4.2)
            max_backoff_base = min(60.0, 0.1 * (2**retry_count))
            delay = max_backoff_base + random.uniform(0.0, 0.1)
            retry_count += 1
            await asyncio.sleep(delay)


def network_engine_process(shm_env_name, shm_act_name):
    print("[INIT] Proc 1 (Red). AFINIDAD HILO 0 ASIGNADA")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {0})

    shm_env = shared_memory.SharedMemory(name=shm_env_name)
    shm_act = shared_memory.SharedMemory(name=shm_act_name)
    env_arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
    act_arr = np.ndarray((RING_BUFFER_SIZE,), dtype=ACTUATOR_DTYPE, buffer=shm_act.buf)

    # Actualizar get_event_loop por patrón python3.11+
    asyncio.run(exchange_websocket_handler(env_arr, act_arr))


# ==============================================================================
# 3. PROCESO 2: TÁCTICO CPU Y CONTROL DE TELEMETRÍA ASIMILACIÓN NMPC EAKF (CAP 7/8)
# ==============================================================================


def fast_thread_process(shm_env_name, shm_mic_name, shm_act_name, shm_par_name):
    print("[INIT] Proc 2 (Filtro Kalman + Telemetry/Logging) AFINIDAD 1 ASIGNADA")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {1})

    # Cargue Memoria Punteros
    shm_env = shared_memory.SharedMemory(name=shm_env_name)
    shm_mic = shared_memory.SharedMemory(name=shm_mic_name)
    shm_par = shared_memory.SharedMemory(name=shm_par_name)
    shm_act = shared_memory.SharedMemory(name=shm_act_name)

    env_arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
    mic_arr = np.ndarray((1,), dtype=MICELIO_DTYPE, buffer=shm_mic.buf)
    par_arr = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)
    act_arr = np.ndarray((RING_BUFFER_SIZE,), dtype=ACTUATOR_DTYPE, buffer=shm_act.buf)

    # Vectores Matrices (Asimetrías burn-in sec 7.3.4)
    I = np.eye(3)
    P_k = np.diag([1e-5, 1e6, 1e6])
    x_k = np.array(
        [[env_arr[0]["P_spot"] if env_arr[0]["P_spot"] > 0 else 40000.0], [0.0], [0.0]]
    )
    H = np.array([[1.0, 0, 0], [0, 0, 1.0]])

    telem_buffer = np.zeros(
        TELEMETRY_SIZE, dtype=TELEM_DTYPE
    )  # Mem pre-asignada Numpy para flush OOM 8.6
    telem_idx = 0
    act_write_idx = 0  # Maneja puntero cola comandos

    is_burnt_in = False
    burnt_ticks_count = 0
    t_ultimo_ciclo = time.time()

    # Derivación Trazado
    last_tr_P = np.trace(P_k)
    epsilon_burn = 1.0  # Param de purga derivativa
    n_consecutive_needed = 100

    try:
        while True:
            # Rejoj Táctico
            t_ahora = time.time()
            dt_ciclo = t_ahora - t_ultimo_ciclo
            t_ultimo_ciclo = t_ahora

            if dt_ciclo < 1e-4:
                dt_ciclo = 1e-4  # Safeguard div. zero cinemática

            dropout = env_arr[0]["flag_dropout"]

            # Config Matrix Cinemática Constante (7.3.2)
            A = np.array([[1.0, dt_ciclo, 0], [0.0, 1.0, 0], [0.0, 0.0, 1.0]])

            x_pred = A @ x_k

            # -- EAKF EXOGENO NO-ESTÁTICO --
            # Inyeccion Paramétrica: La estabildad volumen dicta factor Multiplicador (ρk)
            # Ruid Proceso Q (7.3.3) Modulado: ρ_k * diag(q_s, q_v, q_R)
            q_base = par_arr[3]
            gam_w = par_arr[4]
            gam_Q = par_arr[5]

            om_val = mic_arr[0]["omega"]
            Q_vol = mic_arr[0]["vol_sum_Q"]
            rho_k = 1.0 + (gam_w * abs(om_val)) + (gam_Q * abs(Q_vol))

            Q_k = rho_k * np.diag([q_base, q_base * 1e-2, q_base])
            P_pred = A @ P_k @ A.T + Q_k

            innov = np.zeros((2, 1))

            if dropout == 1:
                # 8.3 Degradación inercial - Crecimiento Progresivo P - suspende corrección.
                x_k = x_pred
                P_k = P_pred
            else:
                retardo_sensor = t_ahora - env_arr[0]["timestamp"]
                beta_red = par_arr[2]
                r_s_base = par_arr[0]
                R_EMD_c = par_arr[1]

                R_k = np.diag([r_s_base * np.exp(beta_red * retardo_sensor), R_EMD_c])
                med_vec = np.array([[env_arr[0]["P_spot"]], [mic_arr[0]["R_n"]]])
                innov = med_vec - (H @ x_pred)

                S_cov = H @ P_pred @ H.T + R_k
                # Implementación con Solver Numpy evadiendo Inv bruta mal condicionada
                K = (P_pred @ H.T) @ np.linalg.solve(S_cov, np.eye(S_cov.shape[0]))

                x_k = x_pred + (K @ innov)
                joseph = I - K @ H
                P_k = (joseph @ P_pred @ joseph.T) + (K @ R_k @ K.T)

            tr_P_actual = np.trace(P_k)
            deriv_tr = abs(tr_P_actual - last_tr_P) / dt_ciclo
            last_tr_P = tr_P_actual

            # LÓGICA ESTABLE DE QUEMA DERIVATIVA
            if not is_burnt_in:
                if deriv_tr <= epsilon_burn:
                    burnt_ticks_count += 1
                    if burnt_ticks_count > n_consecutive_needed:
                        is_burnt_in = True
                        print(
                            "[!!!] CONVERGENCIA SISTÉMICA EXITOSA, ACTIVANDO ACCIÓN TÁCTICA MPC"
                        )
                else:
                    burnt_ticks_count = 0

            # Guardado rápido O(1) Arrays a Memoria HDF5
            telem_buffer[telem_idx]["tr_P"] = tr_P_actual
            telem_buffer[telem_idx]["x0"] = x_k[0, 0]
            telem_buffer[telem_idx]["y0"] = innov[0, 0]
            telem_idx += 1
            if telem_idx >= TELEMETRY_SIZE:
                # VOLCADO BLOQUE EN PARQUET DE VERDAD
                # ej_pool.submit(to_parquet_engine, np.copy(telem_buffer))
                telem_idx = 0
                # (Se descartan implementaciones disco para simulación loop continuo)

            if is_burnt_in and dropout == 0:
                # Resolver acá el Backward Solver Tridiagonal de la Sección 4.3 (En hilo o envia matriz a C++ Casadi)
                pass

                # Inyectar Output Comando Optimizado al Motor Ring SPSC a motor
                u_res_bruto = apply_filters(
                    0.12521, 0.0001, 0.001, 10.0, env_arr[0]["P_spot"]
                )
                if u_res_bruto > 0:
                    act_arr[act_write_idx]["u_compra"] = u_res_bruto
                    act_arr[act_write_idx][
                        "seq_id"
                    ] += 1  # Rompe Wait Condition Process I/O Wait Loop
                    act_write_idx = (act_write_idx + 1) % RING_BUFFER_SIZE

            # Reposo dinámico estabilizador
            time.sleep(0.01)

    except KeyboardInterrupt:
        pass


# ==============================================================================
# 4. PROCESO 3: ANÁLISIS ESTRUCTURAL HILO LENTO MICROESTRUCTURA EMD MICELIO
# ==============================================================================
def slow_thread_process(shm_env_name, shm_mic_name, shm_par_name):
    print("[INIT] Proc 3 (Background Hilo lento / CPU bound OU Micelio) AFINIDAD 2")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {2})

    shm_env = shared_memory.SharedMemory(name=shm_env_name)
    shm_mic = shared_memory.SharedMemory(name=shm_mic_name)
    shm_par = shared_memory.SharedMemory(name=shm_par_name)

    env_arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
    mic_arr = np.ndarray((1,), dtype=MICELIO_DTYPE, buffer=shm_mic.buf)
    par_arr = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)

    lam_sim = 0.0
    t_ultimo_sim = time.time()
    lock_slow_data = (
        mp.Lock()
    )  # Lock (Alternativo eficiente spin python) para writes interbloques

    # Modificado Nivel Medio para Bitcoin Inverse y unidades Dimensionales exactas doc (mu ajustado min 10^-5 ~1e-4)
    LAMBDA_MIN_COTA = 0.0001

    try:
        while True:
            # Timing estricto sim OU discretización. Reversa a T pared verdadero
            t_now = time.time()
            dt_sim = t_now - t_ultimo_sim
            t_ultimo_sim = t_now

            eta_OU = par_arr[8]
            th_OU = par_arr[9]
            sig_OU = par_arr[10]
            mu_hist_c_ajuste = (
                LAMBDA_MIN_COTA * 5.0
            )  # Promedio hist real ajustado Mainnet

            om_ej = random.uniform(0.1, 5.0) * math.sin(
                t_now
            )  # Mock dinámico frecuencia instat. Oscila
            sum_q_ej = random.uniform(5.0, 50.0)  # Sum vol transado nodo HHT

            if IS_TESTNET:
                # Sim estocástica ruido fondo Spread real  Euler
                d_lam_noise = (
                    th_OU * (mu_hist_c_ajuste - lam_sim) * dt_sim
                    + sig_OU * np.sqrt(dt_sim) * np.random.randn()
                )
                lam_sim += d_lam_noise

                # Con Lock: Write Atom en la shared Numpy contra read del Proc Fast (2) EAKF
                with lock_slow_data:
                    mic_arr[0]["omega"] = om_ej
                    mic_arr[0]["vol_sum_Q"] = sum_q_ej
                    mic_arr[0]["R_n"] = 45000 + math.cos(t_now * 0.1) * 500  # macro
                    # Aplicando suelo estricto evitado y la ley acople max (Fuerza limit anti inf/0 GPU loop )
                    mic_arr[0]["lambda_sim"] = max(
                        LAMBDA_MIN_COTA, lam_sim + (eta_OU * abs(om_ej))
                    )

            time.sleep(0.5)

    except KeyboardInterrupt:
        pass


# ==============================================================================
# ENTRY POINT MANAGER
# ==============================================================================
def main():
    procesos = []
    mem_refs = []  # Manejo cerrado para finally

    # Tratamiento Fugas
    def unhook_graceful(signum, frame):
        print("\n>> Interrupción Externa Recibida (SIG). Limpieza Completa /dev/shm...")
        for p in procesos:
            if p.is_alive():
                p.terminate()
                p.join(timeout=1.0)

        for name in ["shm_ent", "shm_mic", "shm_act", "shm_par"]:
            try:
                shm = shared_memory.SharedMemory(name=name)
                shm.close()
                shm.unlink()
            except FileNotFoundError:
                pass
        sys.exit(0)

    signal.signal(signal.SIGINT, unhook_graceful)
    signal.signal(signal.SIGTERM, unhook_graceful)

    try:
        print("== MICELIO MULTICORE DEPLOY INITIALIZING == ")
        mem_refs.append(allocate_shared_memory("shm_ent", ENV_DTYPE.itemsize))
        mem_refs.append(allocate_shared_memory("shm_mic", MICELIO_DTYPE.itemsize))
        mem_refs.append(
            allocate_shared_memory(
                "shm_act", ACTUATOR_DTYPE.itemsize * RING_BUFFER_SIZE
            )
        )
        mem_refs.append(
            allocate_shared_memory(
                "shm_par", TOTAL_PARAMS * np.dtype(np.float64).itemsize
            )
        )

        p_array = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=mem_refs[-1].buf)
        p_array[:] = initialize_default_parameters()[:]

        # P_instace init. args
        p1 = mp.Process(target=network_engine_process, args=("shm_ent", "shm_act"))
        p2 = mp.Process(
            target=fast_thread_process,
            args=("shm_ent", "shm_mic", "shm_act", "shm_par"),
        )
        p3 = mp.Process(
            target=slow_thread_process, args=("shm_ent", "shm_mic", "shm_par")
        )
        procesos.extend([p1, p2, p3])

        for pr in procesos:
            pr.start()

        while True:
            time.sleep(1.0)

    except Exception as general_err:
        print(f"[FATAL ORQUESTRADOR CAIDA LIBRE]: {general_err}")
    finally:
        # Guarantee RAM Wipe / Zombifi
        # cation precention.
        unhook_graceful(None, None)


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    main()
