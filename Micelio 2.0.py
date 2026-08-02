"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Repositorio: samhr07/Cancion-de-Micelio
Módulo: Micelio.py (Núcleo Operativo — Orquestador consolidado)

Documento de diseño: Trading_Bot.pdf (Samuel Hoyos R, Julio 2026).
Diagnóstico de referencia: CLAUDE.md.

------------------------------------------------------------------------------
CAMBIOS DE ESTA REVISIÓN (según el "Orden de trabajo recomendado" de CLAUDE.md)
------------------------------------------------------------------------------
 1. INVERSIÓN TEMPORAL DE LA EDP DE LOEPER.  El esquema de la Sec. 7.4.3 marcha
    hacia adelante en t desde t0; de la Ec. 4 la parábola es *backward* y marchar
    hacia adelante la vuelve antidifusiva (blowup ×5 por paso documentado en
    CLAUDE.md, con D = 0.999998 — es decir, NO era la Singularidad de Loeper).
    Ahora se define condición terminal en t0 + Δτ_pred y se integra hacia t0.
    Ver `resolver_malla_loeper`.
 2. TIEMPO DUAL (Sec. 4.5) explícito: el eje difusivo avanza en ticks
    (ΔT̄ = ν·Δt) y el descuento financiero en años. Se elimina la ambigüedad de
    unidades señalada en CLAUDE.md.
 3. Bloque de parámetros de hot-reloading ampliado, indexado por constantes con
    nombre y LEÍDO en todo el programa (antes se ignoraba por completo).
    Incluye κ (Sec. 6.1), que el PDF nunca declara, y μ_OU (antes hardcodeado).
 4. Se añade `allocate_shared_memory`, que `main()` invocaba sin estar definida.
 5. Seqlock real sobre el bloque del Micelio (Sec. 7.6.2). El `mp.Lock()` que
    se creaba dentro del Hilo Lento era local al proceso y no protegía nada.
 6. Ring Buffer SPSC con secuencia monótona global (antes se comparaba
    `seq_id` contra un índice de slot, lo cual nunca era correcto).
 7. Ω deja de ser ruido: se calcula por la Ec. de la Sec. 1.4 a partir de
    Φ, Ψ y ΔS. ρ_k pasa a usar ω_m (Sec. 7.3.3), no Ω.
 8. Telemetría completa (x_k y la innovación ỹ_k, Sec. 8.6.1) con volcado real
    a disco en hilo secundario de I/O (Sec. 8.6.2).
 9. NMPC conectado a la malla: el actuador ya no emite un 0.12521 hardcodeado.

Toda suposición que rellena un hueco del PDF va marcada con
`# NOTA DE INTERPRETACION:`.

DEPENDENCIAS: este archivo solo requiere NumPy. Los núcleos que el PDF asigna a
CUDA (@cuda.jit, Sec. 7.4) y a acados/CasADi (Sec. 7.5) están implementados como
referencia CPU con layout de arrays idéntico, para que el porte sea mecánico.
"""

import multiprocessing as mp
from multiprocessing import shared_memory
import asyncio
import signal
import sys
import threading
import time
import math
import random
import os
import numpy as np

# ==============================================================================
# 1. ESTRUCTURACIÓN DE BLOQUES LOCK-FREE Y C-TYPES
# ==============================================================================
SYMBOL = "BTCUSDT"
IS_TESTNET = True

SEGUNDOS_POR_ANIO = 31_557_600.0  # Año juliano; escala Ticks/s -> Ticks/Año (4.6)

# Entorno: tensor de memoria de estado temporal escrito por el Motor de Red.
ENV_DTYPE = np.dtype(
    [
        ("P_spot", np.float64),
        ("Q_transado", np.float64),  # Volumen acumulado bruto (USD)
        ("n_ticks", np.uint64),  # Contador monótono: alimenta ν en el Hilo Lento
        ("flag_dropout", np.int8),
        ("timestamp", np.float64),
        ("inv_confirmado", np.float64),  # Inventario real ejecutado (BTC), Sec. 4.5
        ("seq_fill", np.uint64),  # Publicación del fill (orden de escritura)
    ]
)

# Estructural: escrito por el Hilo Lento, leído por el Hilo Rápido.
# El campo `seq` implementa el seqlock exigido por la Sec. 7.6.2 (spinlocks en
# lugar de mutexes) para que las escrituras multi-campo no sean carreras de datos.
MICELIO_DTYPE = np.dtype(
    [
        ("seq", np.uint64),  # Par = consistente, impar = escritura en curso
        ("Omega", np.float64),  # Ω  estabilidad del modelo   [BTC/Ticks²]
        ("w_m", np.float64),  # ω_m frecuencia de mercado   [1/Ticks]
        ("nu", np.float64),  # ν   tasa de llegada de ticks [Ticks/Años]
        ("lambda_sim", np.float64),  # λ_sim falta de liquidez      [1/BTC]
        ("R_n", np.float64),  # Residuo macro de la EMD      [USD/BTC]
        ("S_ref", np.float64),  # Precio anclado al último nodo de fase (2.6)
        ("delta_S", np.float64),  # ΔS = S - S_ref               [USD/BTC]
        ("vol_sum_Q", np.float64),  # ΣQ del ciclo vigente (6.3)   [USD]
    ]
)

# Actuador: Ring Buffer SPSC lock-free (Sec. 7.6.2).
RING_BUFFER_SIZE = 16
ACTUATOR_DTYPE = np.dtype(
    [
        ("u_compra", np.float64),
        ("u_venta", np.float64),
        ("ts_emision", np.float64),
        ("seq_id", np.uint64),  # Secuencia global monótona; 0 = slot vacío
    ]
)

# Telemetría diferida (Sec. 8.6.1): estado, innovación e incertidumbre.
TELEMETRY_SIZE = 10000
TELEM_DTYPE = np.dtype(
    [
        ("t_wall", np.float64),
        ("tr_P", np.float64),
        ("x0", np.float64),  # S   (precio verdadero)
        ("x1", np.float64),  # v   (velocidad instantánea)
        ("x2", np.float64),  # R_n (residuo macro)
        ("y0", np.float64),  # ỹ_k componente precio
        ("y1", np.float64),  # ỹ_k componente residuo  -> Covariance Matching (3.5)
    ]
)

DIR_TELEMETRIA = "telemetria"

# ------------------------------------------------------------------------------
# 1.1 Bloque de parámetros de HOT-RELOADING
# ------------------------------------------------------------------------------
# CLAUDE.md: "Los 11 parámetros de hot-reloading documentados en el encabezado no
# se leen en ninguna parte del programa". Se amplía el bloque a todo lo que el
# modelo necesita y se indexa por constante con nombre (los índices desnudos
# `par_arr[3]` eran la causa de que nadie los conectara).
P_R_S_BASE = 0  # r_S,base      varianza mínima del tick size   (7.3.3)
P_R_EMD = 1  # r_EMD         varianza histórica de Hilbert   (7.3.3)
P_BETA_JITTER = 2  # β             sensibilidad a la latencia      (6.5)
P_Q_BASE = 3  # q_S,q_v,q_Rn  varianzas base del modelo térmico (7.3.3)
P_GAMMA_OMEGA = 4  # γ_ω           peso de ω_m en ρ_k              (7.3.3)
P_GAMMA_Q = 5  # γ_Q           peso de ΣQ  en ρ_k              (7.3.3)
P_R_BASE = 6  # R_base        comisión taker estática          (6.1)
P_MU_OU = 7  # μ_OU          media histórica de falta de liquidez (8.1.1)
P_ETA_IMPACTO = 8  # η             sensibilidad estructural         (8.1.2)
P_THETA_OU = 9  # θ_OU          tasa de reversión a la media     (8.1.1)
P_SIGMA_OU = 10  # σ_OU          volatilidad intrínseca del spread(8.1.1)
P_KAPPA = 11  # κ             escalado de aversión R(Ω)        (6.1)
P_MU_INV = 12  # μ             sensibilidad al riesgo de exposición (6.2)
P_Q_INV_BASE = 13  # q_base        costo asintótico de inventario   (6.2)
P_Q_DELTA = 14  # q_Δ           rigidez del seguimiento de cobertura (7.5.2)
P_P_DELTA = 15  # p_Δ           costo terminal de cobertura      (7.5.2)
P_P_INV = 16  # p_inv         costo terminal de inventario     (7.5.2)
P_OMEGA_CRIT = 17  # Ω_crit        umbral de transición de fase    (6.2)
P_K_VOL = 18  # k             c²_vol = k·ω_m·ν                (4.5)
P_R_USD = 19  # r_USD         tasa libre de riesgo USD        (4.4)
P_Q_BTC = 20  # q_BTC         tasa de conveniencia BTC        (4.4)
P_LAMBDA_MIN = 21  # λ_min         suelo de fricción               (8.1.2)
P_GAMMA_0 = 22  # γ₀            curvatura de la condición terminal
P_DS_MAX = 23  # ΔS_max        margen relativo del dominio      (7.4.1)
P_TAU_PRED = 24  # Δτ_pred       horizonte de predicción [AÑOS]  (4.5 / 7.4.1)
P_EPS_BURN = 25  # ε_burn        banda muerta de la traza        (7.1)
P_N_BURN = 26  # N             ticks continuos de convergencia (7.1)
P_U_MAX = 27  # U_max         riesgo máximo por operación      (6.2)
P_C_MAX = 28  # C_max         capital disponible en billetera  (6.2)
P_EPS_SING = 29  # ε             margen de seguridad de Loeper   (4.4.3)
P_TAU_MAX = 30  # τ_max         tolerancia máxima de retardo [s] (6.5)
TOTAL_PARAMS = 31

# Discretización de la malla de Loeper (Sec. 7.4.1).
N_S_MALLA = 61  # Impar: garantiza que S_k caiga exactamente en un nodo
N_T_MAX = 256  # Techo de pasos temporales por ciclo de control
N_NMPC = 8  # N pasos del horizonte de predicción (Sec. 5 / 6.1)
ITER_NMPC_MAX = 80  # Iteraciones del descenso proyectado
T_MAX_SOLVER = 0.020  # t_max: fallback determinista U_t = 0 (Sec. 4.5)

PERIODO_HILO_RAPIDO = 0.010  # Cadencia nominal del ciclo de control (s)
PERIODO_HILO_LENTO = 0.500  # Cadencia del análisis estructural (s)

# Parámetros del generador sintético de Testnet (solo mock; ver Sec. 8.1).
PRECIO_BASE_MOCK = 45000.0  # [USD/BTC]
AMPLITUD_CICLO_MOCK = 150.0  # [USD/BTC]
PERIODO_CICLO_MOCK = 40.0  # [s] por ciclo estructural completo


def initialize_default_parameters() -> np.ndarray:
    """Valores por defecto del bloque de hot-reloading.

    NOTA DE INTERPRETACION: el PDF no tabula valores numéricos para casi ningún
    hiperparámetro. Los de abajo están escogidos para ser DIMENSIONALMENTE
    coherentes con la Sec. 4.4.1 y para que λ·S²·Γ quede lejos de 1 a precios de
    BTC (CLAUDE.md: con μ_OU = 0.5 /BTC la singularidad se dispara de inmediato).
    Comprobación con los valores de abajo, S = 45 000 USD/BTC:
        λ·S²·Γ = 5e-4 · 2.03e9 · 2e-8 ≈ 0.020  ✓ (D ≈ 0.98, sin singularidad)
    Todos los marcados [CALIBRAR] deben venir de métricas reales de Mainnet.
    """
    params = np.zeros(TOTAL_PARAMS, dtype=np.float64)

    # --- EAKF (Sec. 7.3.3) ---
    params[P_R_S_BASE] = 0.5  # ~ (tick size)² del par BTCUSDT
    params[P_R_EMD] = 1.0
    params[P_BETA_JITTER] = 0.02
    params[P_Q_BASE] = 0.01
    # NOTA DE INTERPRETACION: la Sec. 7.3.3 escribe ρ_k = 1 + γ_ω|ω_m| + γ_Q|ΣQ|
    # pero NO declara las unidades de γ_ω ni de γ_Q (mismo hueco que κ en la 6.1).
    # Para que ρ_k sea adimensional con [ω_m] = 1/Ticks y [ΣQ] = USD (Sec. 1.1),
    # necesariamente [γ_ω] = Ticks y [γ_Q] = 1/USD. Se calibran para que cada
    # término aporte ~1.0 en régimen típico; con γ_Q = 0.2 y ΣQ ~ 1e5 USD daba
    # ρ_k ≈ 2e4, y el burn-in de la Sec. 7.1 nunca alcanzaba la banda muerta.
    params[P_GAMMA_OMEGA] = 500.0  # [Ticks]  -> ω_m = 2e-3 aporta 1.0
    params[P_GAMMA_Q] = 1.0e-5  # [1/USD]  -> ΣQ = 100 000 USD aporta 1.0
    # Sec. 8: desde Colombia la señal tarda ~300 ms hasta Tokio/Europa; τ_max debe
    # dejar margen sobre ese piso físico sin admitir datos tácticamente obsoletos.
    params[P_TAU_MAX] = 2.0  # [s]

    # --- Costos y aversión (Sec. 6.1 / 6.2) ---
    params[P_R_BASE] = 0.05
    params[P_KAPPA] = 1.0  # [Ticks⁴/BTC²]  <- κ NO está declarado en el PDF
    params[P_MU_INV] = 1.0
    params[P_Q_INV_BASE] = 0.05
    params[P_Q_DELTA] = 1.0
    params[P_P_DELTA] = 20.0  # p_Δ >> q_Δ exigido por la Sec. 7.5.2
    params[P_P_INV] = 5.0  # p_inv >> q_base
    params[P_OMEGA_CRIT] = 2.0  # [CALIBRAR]

    # --- Ornstein-Uhlenbeck de liquidez (Sec. 8.1.1) ---
    params[P_MU_OU] = 5.0e-4  # [1/BTC] [CALIBRAR con métricas de Mainnet]
    # NOTA DE INTERPRETACION: la Sec. 8.1.2 declara [η] = Ticks²/BTC² pero no su
    # magnitud, y en λ_sim = máx(λ_min, λ_ruido + η|Ω|) el término η|Ω| debe ser
    # COMPARABLE a λ_ruido ~ μ_OU, no dominarlo. Con η = 0.15 y |Ω| ~ O(1) salía
    # λ_sim ≈ 0.15 /BTC, que da λ·S²·Γ ≈ 6 y mete al solver en singularidad
    # permanente. Se escala η para que |Ω| ~ 1 aporte del orden de μ_OU.
    params[P_ETA_IMPACTO] = 5.0e-4  # [CALIBRAR contra profundidad real del libro]
    params[P_THETA_OU] = 0.5
    # NOTA DE INTERPRETACION: la desviación estacionaria de un OU es σ/√(2θ). Para
    # que λ_ruido fluctúe ALREDEDOR de μ_OU y no lo sepulte hay que exigir
    # σ/√(2θ) ≪ μ_OU. Con σ = 0.05 y θ = 0.5 daba σ_est = 0.05, cien veces μ_OU:
    # λ paseaba por ±0.05 /BTC y cruzaba a negativo constantemente (el suelo λ_min
    # lo tapaba, pero el sesgo positivo resultante disparaba la singularidad).
    # σ = 1e-4 deja σ_est = 1e-4, o sea ~20 % de μ_OU. [CALIBRAR]
    params[P_SIGMA_OU] = 1.0e-4
    params[P_LAMBDA_MIN] = 1.0e-4  # Suelo estricto: evita fricciones negativas

    # --- Minimax-Loeper (Sec. 4.4 / 4.5 / 7.4) ---
    params[P_K_VOL] = 4.0e-8  # k tal que c²_vol ~ 0.5 /Año  [CALIBRAR]
    params[P_R_USD] = 0.05  # 1/Años
    params[P_Q_BTC] = 0.00  # 1/Años
    params[P_GAMMA_0] = 2.0e-8  # [BTC³/USD²]  curvatura de la cobertura
    params[P_DS_MAX] = 0.02  # ±2 % alrededor de S_k
    params[P_TAU_PRED] = 30.0 / SEGUNDOS_POR_ANIO  # 30 s expresados en AÑOS
    params[P_EPS_SING] = 1.0e-3

    # --- Burn-in (Sec. 7.1) ---
    # [CALIBRAR — REQUERIDO ANTES DE OPERAR EN ENTORNO REAL]
    # La Sec. 7.1 define el criterio |Tr(P_k)-Tr(P_{k-1})|/Δt ≤ ε_burn sostenido N
    # ticks continuos, pero NO fija ε_burn ni N. El valor de abajo se midió contra
    # el generador sintético de Testnet: en régimen estacionario dTr/dt vive en
    # ~0.01-0.03 y el burn-in cierra en ~215 ciclos (≈2.4 s).
    #
    # ADVERTENCIA: ε_burn se está tratando aquí como una CONSTANTE ESTÁTICA, y muy
    # probablemente NO lo sea al entrar al mercado real. La derivada de la traza
    # depende de Δt (que en producción no es la cadencia nominal sino la latencia
    # real medida: desde Colombia ~300 ms, Sec. 8), de ρ_k —que el Micelio modula
    # tick a tick vía ω_m y ΣQ (Sec. 7.3.3)— y del régimen de volatilidad vigente.
    # Un umbral fijo calibrado en un régimen tranquilo bloqueará el burn-in en uno
    # agitado, y uno calibrado en régimen agitado dejará pasar un filtro que aún no
    # convergió. Lo correcto en Mainnet es hacerlo ADAPTATIVO: normalizarlo contra
    # la escala de Tr(P) (criterio relativo) o contra un percentil móvil de la
    # propia derivada, en lugar de compararlo contra un absoluto.
    params[P_EPS_BURN] = 1.0
    params[P_N_BURN] = 100.0

    # --- Límites operativos (Sec. 6.2) ---
    params[P_U_MAX] = 0.05  # BTC por operación
    params[P_C_MAX] = 0.50  # BTC de capital desplegable
    return params


# ==============================================================================
# 1.2 PRIMITIVAS DE CONCURRENCIA
# ==============================================================================
def allocate_shared_memory(nombre: str, n_bytes: int) -> shared_memory.SharedMemory:
    """Reserva (o recicla) un bloque de /dev/shm. Sec. 7.7.2.

    Esta función se invocaba en `main()` sin estar definida: el orquestador moría
    con NameError antes de arrancar ningún proceso.
    """
    try:
        shm = shared_memory.SharedMemory(name=nombre, create=True, size=n_bytes)
    except FileExistsError:
        # Residuo de una ejecución anterior mal terminada: se recicla y se limpia.
        previo = shared_memory.SharedMemory(name=nombre)
        previo.close()
        previo.unlink()
        shm = shared_memory.SharedMemory(name=nombre, create=True, size=n_bytes)
    shm.buf[:n_bytes] = b"\x00" * n_bytes
    return shm


def escribir_micelio_seqlock(mic_arr, campos: dict) -> None:
    """Escritura atómica por seqlock del bloque estructural (Sec. 7.6.2).

    Protocolo: seq impar mientras se escribe, par cuando el registro es coherente.
    NOTA DE INTERPRETACION: en x86-64 el modelo de memoria (TSO) no reordena
    store-store, que es la única garantía que este protocolo necesita del lado
    del escritor. Una barrera de liberación explícita exigiría una extensión en C.
    """
    mic_arr[0]["seq"] += 1  # -> impar: escritura en curso
    for clave, valor in campos.items():
        mic_arr[0][clave] = valor
    mic_arr[0]["seq"] += 1  # -> par: registro consistente


def leer_micelio_seqlock(mic_arr, reintentos: int = 8):
    """Lectura lock-free del bloque estructural. Devuelve (snapshot, es_valido)."""
    for _ in range(reintentos):
        s1 = int(mic_arr[0]["seq"])
        if s1 & 1:  # Escritor dentro de la sección crítica
            continue
        snap = mic_arr[0:1].copy()
        if int(mic_arr[0]["seq"]) == s1:
            return snap[0], True
    # Degradación: se devuelve la última lectura aunque pueda estar rasgada.
    return mic_arr[0:1].copy()[0], False


def apply_filters(
    u_raw: float, stepSize: float, minQty: float, minNotional: float, Pspot: float
) -> float:
    """Cuantización comercial del vector de control (Sec. 8.5).

    Orden correcto: se discretiza con floor y SOLO DESPUÉS se revalidan minQty y
    minNotional, porque el floor puede sacar la orden por debajo del nocional.
    """
    if u_raw <= 0.0 or stepSize <= 0.0 or Pspot <= 0.0:
        return 0.0
    u_rounded = math.floor(u_raw / stepSize) * stepSize
    if u_rounded < minQty:
        return 0.0
    if (u_rounded * Pspot) < minNotional:
        return 0.0
    return u_rounded


# ==============================================================================
# 2. EDP DE MINIMAX-LOEPER (Sec. 7.4) — MARCHA INVERTIDA
# ==============================================================================
def resolver_malla_loeper(S_k, tr_P, c2_vol, lam, S_ref, par):
    """Resuelve la EDP de Minimax-Loeper y devuelve las superficies M_Γ y M_Δ.

    ---------------------------------------------------------------------------
    CORRECCIÓN CENTRAL (CLAUDE.md, "Hallazgo crítico")
    ---------------------------------------------------------------------------
    La Sec. 7.4.3 escribe el paso discreto como

        U[i,j+1] = U[i,j] - dt·( ½·σ²·S_i²·Γ[i,j]/D[i,j] - costos )

    es decir, marchando HACIA ADELANTE en t desde t0. Pero de la Ec. 4,

        ∂U/∂t = -½·c²_vol·S²·Γ/D - (r_USD - q_BTC)·S·∂U/∂S + q_BTC·U

    es una parábola BACKWARD: se resuelve retrocediendo desde una condición
    terminal, como todo pricing de opciones. Marchar hacia adelante equivale a
    resolver la ecuación del calor con el tiempo invertido: es antidifusiva y por
    tanto mal planteada, y amplifica el ruido de escala de malla con factor ~5
    por paso (evidencia experimental en CLAUDE.md, con D = 0.999998 — o sea que
    aquel blowup NO era la Singularidad de Loeper).

    Aquí se integra en tiempo restante τ = t_terminal - t, que avanza HACIA
    ADELANTE mientras t retrocede:

        ∂U/∂τ = +½·c²_vol·S²·Γ/D + (r_USD - q_BTC)·S·∂U/∂S - q_BTC·U

    Con esta orientación el coeficiente difusivo es positivo, el problema es bien
    planteado y la condición CFL recupera su sentido (CLAUDE.md: "CFL presupone
    un problema bien planteado de entrada").

    ---------------------------------------------------------------------------
    TIEMPO DUAL (Sec. 4.5) — resuelve la ambigüedad de unidades de CLAUDE.md
    ---------------------------------------------------------------------------
    La Sec. 4.5 separa ν·∂U/∂T̄ + ∂U/∂t. Sobre la trayectoria del reloj del
    sistema, la derivada total es dU/dt = ∂U/∂t + ν·∂U/∂T̄, de modo que un ÚNICO
    eje temporal en AÑOS es suficiente si —y solo si— el término difusivo entra
    ya escalado a años. Eso es exactamente lo que hace c²_vol = k·ω_m·ν
    ([1/Ticks]·[Ticks/Años] = [1/Años]). Por eso aquí el eje t se mide en años y
    ν entra únicamente a través de c²_vol: no hay dos ejes que sincronizar.

    Devuelve
    --------
    dict con:
        M_Gamma : (N_S, n_t)  convexidad Γ(S,t)      [BTC³/USD²]
        M_Delta : (N_S, n_t)  cobertura ∂U/∂S        [BTC²/USD]
        S_malla : (N_S,)      eje de precios         [USD/BTC]
        D_min   : float       mínimo de 1 - λS²Γ sobre toda la malla
        singular: bool        True si D_min <= ε  (Sec. 4.4.3)
        tau_eff : float       horizonte efectivamente integrado [Años]
    """
    n_s = N_S_MALLA
    ds_max_rel = par[P_DS_MAX]
    gamma_0 = par[P_GAMMA_0]
    r_usd = par[P_R_USD]
    q_btc = par[P_Q_BTC]
    eps_sing = par[P_EPS_SING]
    tau_pred = par[P_TAU_PRED]

    # --- Eje de precio (Sec. 7.4.1) -----------------------------------------
    S_min = S_k * (1.0 - ds_max_rel)
    S_sup = S_k * (1.0 + ds_max_rel)
    ds = (S_sup - S_min) / (n_s - 1)
    S = S_min + ds * np.arange(n_s, dtype=np.float64)
    S2 = S * S
    idx_centro = (n_s - 1) // 2  # S[idx_centro] == S_k exactamente (n_s impar)

    # --- Varianza instantánea efectiva --------------------------------------
    # Sec. 8.3.1: "Este aumento de P se transmite directamente a la EDP de Loeper
    # a través del término de varianza instantánea σ²", pero el PDF no da la
    # fórmula del acople.
    # NOTA DE INTERPRETACION: se ensancha c²_vol con la incertidumbre relativa de
    # precio, p_S/S², que es adimensional ([USD²/BTC²] / [USD²/BTC²]). Durante un
    # dropout Tr(P) crece y el cono de riesgo se abre solo, como pide la Sec. 8.3.
    c2_eff = c2_vol * (1.0 + max(0.0, tr_P) / (S_k * S_k))

    # --- Condición terminal --------------------------------------------------
    # NOTA DE INTERPRETACION: la Sec. 7.4 nunca especifica U en el borde temporal
    # (CLAUDE.md, contradicción interna #3: "el solver de GPU es irreproducible").
    # Se asume el payoff cuadrático de cobertura U = ½·γ₀·(S - S_ref)², anclado no
    # al precio actual sino al NODO DE FASE S_ref de la Sec. 2.6. Anclarlo en S_k
    # daría Δ(S_k) ≡ 0 por simetría y el objetivo del NMPC sería trivialmente
    # nulo; anclarlo en S_ref hace que la cobertura siga el desplazamiento
    # estructural ΔS = S - S_ref, que es la variable que el modelo de la Canción
    # del Micelio define como significativa.
    #   Unidades: [γ₀] = BTC³/USD², luego [U] = BTC ✓ y [Γ] = BTC³/USD² ✓ (4.4.1).
    U = 0.5 * gamma_0 * (S - S_ref) ** 2

    # --- Condición CFL (Sec. 7.4.1), ahora sí legítima ----------------------
    # Cota a-priori de la difusividad usando la Γ terminal (γ₀), que es la máxima
    # curvatura del perfil: la integración backward solo la suaviza.
    D_ref = 1.0 - lam * (S_sup * S_sup) * gamma_0
    if D_ref <= eps_sing:
        # Singularidad de Loeper ya en la condición terminal: no hay malla que
        # calcular. Sec. 4.4.3 -> inacción total.
        return {
            "M_Gamma": np.zeros((n_s, 1)),
            "M_Delta": np.zeros((n_s, 1)),
            "S_malla": S,
            "idx_centro": idx_centro,
            "D_min": D_ref,
            "singular": True,
            "tau_eff": 0.0,
            "n_t": 1,
        }

    a_max = 0.5 * c2_eff * (S_sup * S_sup) / D_ref
    dt_cfl = 0.45 * ds * ds / a_max if a_max > 0.0 else tau_pred

    n_t = int(math.ceil(tau_pred / dt_cfl)) + 1 if dt_cfl > 0.0 else 2
    n_t = int(min(max(n_t, 2), N_T_MAX))
    dt = tau_pred / (n_t - 1)
    if dt > dt_cfl:
        # Techo de pasos alcanzado. Se prefiere integrar un horizonte MÁS CORTO
        # de forma estable antes que violar CFL: el NMPC recibe tau_eff y sabe
        # hasta dónde puede confiar en la malla.
        dt = dt_cfl
    tau_eff = dt * (n_t - 1)

    M_Gamma = np.zeros((n_s, n_t), dtype=np.float64)
    M_Delta = np.zeros((n_s, n_t), dtype=np.float64)
    D_min = 1.0
    singular = False

    Gamma = np.zeros(n_s, dtype=np.float64)
    Delta = np.zeros(n_s, dtype=np.float64)

    # --- Marcha BACKWARD: j = n_t-1 (terminal) -> j = 0 (t0) ------------------
    for j in range(n_t - 1, -1, -1):
        # Γ por diferencias centrales (Sec. 7.4.2).
        Gamma[1:-1] = (U[2:] - 2.0 * U[1:-1] + U[:-2]) / (ds * ds)
        # Condición de linealidad en los bordes: Γ = 0 (ya está a cero).
        Gamma[0] = 0.0
        Gamma[-1] = 0.0

        Delta[1:-1] = (U[2:] - U[:-2]) / (2.0 * ds)
        Delta[0] = (U[1] - U[0]) / ds
        Delta[-1] = (U[-1] - U[-2]) / ds

        M_Gamma[:, j] = Gamma
        M_Delta[:, j] = Delta

        # Factor de fricción D = 1 - λ·S²·Γ (Sec. 7.4.3 / 8.1.3).
        D = 1.0 - lam * S2 * Gamma
        D_min = min(D_min, float(D[1:-1].min()))
        if D_min <= eps_sing:
            # Sec. 4.4.3: el solver detecta la singularidad y no busca gradiente.
            singular = True
            M_Gamma[:, :j] = Gamma[:, None]
            M_Delta[:, :j] = Delta[:, None]
            break

        if j == 0:
            break

        # Paso explícito en τ (equivale a retroceder un dt en t).
        difusion = 0.5 * c2_eff * S2 * Gamma / D
        deriva = (r_usd - q_btc) * S * Delta
        descuento = q_btc * U
        U = U + dt * (difusion + deriva - descuento)

    return {
        "M_Gamma": M_Gamma,
        "M_Delta": M_Delta,
        "S_malla": S,
        "idx_centro": idx_centro,
        "D_min": D_min,
        "singular": singular,
        "tau_eff": tau_eff,
        "n_t": n_t,
    }


# ==============================================================================
# 3. CONTROLADOR NMPC (Sec. 6.1 / 6.2 / 7.5)
# ==============================================================================
def resolver_nmpc(I_0, malla, S_k, Omega, lam, par, t_inicio):
    """Optimización de trayectoria sobre N pasos. Devuelve (u_compra_0, u_venta_0).

    Formulación (reconciliación de la contradicción #2 de CLAUDE.md):
    la Sec. 7.5.1 define x^c_k ∈ R^(1×1) = [I_k] pero la 7.5.2 define
    Q(Ω) ∈ R^(2×2), y la 7.5.3 escribe (x^c_k)ᵀ Q(Ω) x^c_k — un escalar no se
    multiplica por una 2×2. El vector penalizado es [e_k ; I_k] ∈ R², con
    Q = diag(q_Δ, q_inv(Ω)) y e_k = I_k - Δ_k (relación algebraica ya dada en
    7.5.1). Eso además reconcilia la 6.1, que formula el costo sobre e_k.

        mín J = Σ_{k=0}^{N-1} [ q_Δ·e_k² + q_inv(Ω)·I_k² + R_eff,k·(u_c,k² + u_v,k²) ]
                + p_Δ·e_N² + p_inv·I_N²
        s.a.  I_{k+1} = I_k + u_c,k - u_v,k        (7.5.1)
              0 <= u_k <= ū(Ω, I)                  (6.2)

    NOTA DE INTERPRETACION: el PDF asigna esta resolución a acados/CasADi por SQP
    (Sec. 7.5.3). Aquí se implementa un descenso de gradiente proyectado en NumPy:
    la función de costo es una cuadrática estrictamente convexa con restricciones
    de caja, así que la proyección es exacta y el óptimo es único. Sirve como
    referencia numérica verificable y como fallback cuando acados no está
    disponible; el reemplazo por `solver.solve()` es directo.
    """
    q_delta = par[P_Q_DELTA]
    q_inv = par[P_Q_INV_BASE] + par[P_MU_INV] * Omega * Omega  # q_inv(Ω), Sec. 6.2
    p_delta = par[P_P_DELTA]
    p_inv = par[P_P_INV]
    R_base = par[P_R_BASE]
    kappa = par[P_KAPPA]
    eps_sing = par[P_EPS_SING]

    # --- Freno de emergencia por singularidad (Sec. 4.4.3) -------------------
    if malla["singular"] or malla["D_min"] <= eps_sing:
        return 0.0, 0.0, True

    M_Delta = malla["M_Delta"]
    M_Gamma = malla["M_Gamma"]
    ic = malla["idx_centro"]
    n_t = malla["M_Delta"].shape[1]

    # --- Δ_k y Γ_k como TVP interpolados de la malla (Sec. 6.4 / 7.5.1) ------
    # NOTA DE INTERPRETACION: la Sec. 6.4 pide interpolación por B-splines; aquí
    # se usa interpolación lineal sobre el eje temporal, que basta porque el
    # NMPC evalúa el nodo central del eje de precios (S_k) y no una coordenada
    # arbitraria. Además la Sec. 7.4.4 solo devuelve M_Γ, pero la 7.5.1 exige Δ_k
    # como TVP: la malla debe devolver AMBAS superficies (hueco del documento).
    j_nodos = np.linspace(0.0, n_t - 1, N_NMPC + 1)
    delta_grid = M_Delta[ic, :]
    gamma_grid = M_Gamma[ic, :]
    delta_raw = np.interp(j_nodos, np.arange(n_t), delta_grid)
    gamma_nodo = np.interp(j_nodos, np.arange(n_t), gamma_grid)

    # Reconciliación dimensional: [Δ] = ∂U/∂S = BTC²/USD (Sec. 4.4.1), pero
    # e_k = I_k - Δ_k compara contra un inventario en BTC (Sec. 7.5.1).
    # NOTA DE INTERPRETACION: la cantidad de cobertura denominada en BTC es
    # S·∂U/∂S ([USD/BTC]·[BTC²/USD] = BTC). Se usa esa conversión.
    Delta_obj = S_k * delta_raw  # [BTC]

    # --- R_eff,k = R_base + κΩ² + f_Loeper(S², λ, Γ_k)  (Sec. 7.5.2) ---------
    # NOTA DE INTERPRETACION: el PDF nombra f_Loeper pero no la define. Se toma
    # f_Loeper = R_base·(1/D_k - 1), es decir, el costo de ejecución se infla por
    # el mismo factor de amplificación 1/(1-λS²Γ) con que Loeper amplifica la
    # varianza. Es homogénea con R_base, vale 0 cuando no hay impacto (Γ→0) y
    # diverge exactamente en la singularidad, que es el comportamiento que la
    # Sec. 7.5.2 describe ("dispara el valor de R_eff,k en ese nodo específico,
    # forzando al optimizador a diferir la compra/venta a un nodo futuro").
    D_nodo = np.maximum(1.0 - lam * (S_k * S_k) * gamma_nodo[:N_NMPC], eps_sing)
    f_loeper = R_base * (1.0 / D_nodo - 1.0)
    R_eff = R_base + kappa * Omega * Omega + f_loeper  # (N_NMPC,)

    # --- Restricciones de caja dinámicas asimétricas (Sec. 6.2) --------------
    U_max = par[P_U_MAX]
    C_max = par[P_C_MAX]
    if abs(Omega) < par[P_OMEGA_CRIT]:
        ub_c = min(C_max, U_max)
        ub_v = min(C_max, U_max)
    else:
        # Transición de fase: vaciado direccional forzado.
        if I_0 > 0.0:
            ub_c, ub_v = 0.0, U_max
        elif I_0 < 0.0:
            ub_c, ub_v = U_max, 0.0
        else:
            ub_c, ub_v = 0.0, 0.0
    if ub_c <= 0.0 and ub_v <= 0.0:
        return 0.0, 0.0, False

    # --- Descenso de gradiente proyectado ------------------------------------
    u_c = np.zeros(N_NMPC, dtype=np.float64)
    u_v = np.zeros(N_NMPC, dtype=np.float64)

    # Cota de Lipschitz del gradiente respecto a d = u_c - u_v (conservadora).
    L = 2.0 * N_NMPC * (N_NMPC * (q_delta + q_inv) + p_delta + p_inv)
    L += 2.0 * float(R_eff.max())
    paso = 1.0 / L if L > 0.0 else 1e-3

    for _ in range(ITER_NMPC_MAX):
        if (time.perf_counter() - t_inicio) > T_MAX_SOLVER:
            # Fallback determinista (Sec. 4.5): mantener posición.
            return 0.0, 0.0, False

        d = u_c - u_v
        I = I_0 + np.concatenate(([0.0], np.cumsum(d)))  # I[0..N], I[0] = I_0
        e = I - Delta_obj

        # dJ/dI_k para k = 1..N (I_0 no es variable de decisión).
        dJ_dI = np.empty(N_NMPC + 1, dtype=np.float64)
        dJ_dI[:N_NMPC] = 2.0 * (q_delta * e[:N_NMPC] + q_inv * I[:N_NMPC])
        dJ_dI[N_NMPC] = 2.0 * (p_delta * e[N_NMPC] + p_inv * I[N_NMPC])

        # d_m influye en I_k para todo k > m  ->  suma de cola.
        cola = np.cumsum(dJ_dI[::-1])[::-1]  # cola[k] = Σ_{i>=k} dJ_dI[i]
        g_d = cola[1:]  # (N_NMPC,)

        g_c = g_d + 2.0 * R_eff * u_c
        g_v = -g_d + 2.0 * R_eff * u_v

        u_c = np.clip(u_c - paso * g_c, 0.0, ub_c)
        u_v = np.clip(u_v - paso * g_v, 0.0, ub_v)

    # Horizonte recedente (Sec. 6.4): solo se inyecta U_0, el resto se descarta.
    return float(u_c[0]), float(u_v[0]), False


# ==============================================================================
# 4. PROCESO 1: MOTOR DE RED (ASÍNCRONO I/O, WATCHDOG, TOKEN BUCKET)
# ==============================================================================
class TokenBucket:
    """Limitador de tasa por pesos acumulativos (Sec. 8.2.1).

    T_k = mín(C_max, T_{k-1} + ρ·Δt) - w_k
    """

    def __init__(self, capacidad: float = 1200.0, ventana_seg: float = 60.0):
        self.capacidad = capacidad
        self.rho = capacidad / ventana_seg
        self.tokens = capacidad
        self.t_ultimo = time.monotonic()

    def _recargar(self) -> None:
        ahora = time.monotonic()
        self.tokens = min(
            self.capacidad, self.tokens + self.rho * (ahora - self.t_ultimo)
        )
        self.t_ultimo = ahora

    async def adquirir(self, peso: float, alta_prioridad: bool = False) -> bool:
        """Sec. 8.2.2: baja prioridad se descarta, alta prioridad espera τ_espera."""
        self._recargar()
        if self.tokens >= peso:
            self.tokens -= peso
            return True
        if not alta_prioridad:
            return False  # Sondeo/sincronización: se descarta para no asfixiar la cola
        tau_espera = (peso - self.tokens) / self.rho
        await asyncio.sleep(tau_espera)
        self._recargar()
        self.tokens -= peso
        return True


async def exchange_websocket_handler(env_arr, act_arr, par_arr):
    """Watchdog con reconexión infinita y backoff exponencial (Sec. 8.4.2).

    El `except` vive DENTRO del `while True`: en la revisión anterior estaba fuera,
    de modo que el primer TimeoutError marcaba dropout y terminaba la corrutina,
    dejando al bot ciego para siempre.
    """
    bucket = TokenBucket(capacidad=1200.0, ventana_seg=60.0)
    retry_count = 0
    seq_esperada = 1  # Las secuencias emitidas arrancan en 1 (0 = slot vacío)
    n_ticks = 0
    inventario = 0.0
    seq_fill = 0
    t_origen = time.time()

    while True:
        try:
            # --- Ingesta de mercado (peso bajo, prioridad baja) --------------
            if not await bucket.adquirir(1.0, alta_prioridad=False):
                await asyncio.sleep(0.01)  # Presión mitigada, yield del event loop
                continue

            # (Aquí iría el `await ws.recv()` real con timeout; el sleep lo simula.)
            await asyncio.sleep(0.005)

            ahora = time.time()
            n_ticks += 1
            # PENDIENTE (mock de Testnet): serie sintética con un modo cíclico
            # dominante más microestructura. NO es ruido blanco a propósito: el
            # modelo de la Sec. 1 supone un mercado que revierte a nodos de precio,
            # y el detector de nodos de fase (Sec. 2.6) necesita una señal con
            # ciclo real. Con ruido blanco puro el nodo dispara en casi todos los
            # ticks, ΣQ se reinicia sin parar y el burn-in de la Sec. 7.1 nunca
            # alcanza la banda muerta. NO USAR EN MAINNET.
            fase = 2.0 * math.pi * (ahora - t_origen) / PERIODO_CICLO_MOCK
            env_arr[0]["P_spot"] = (
                PRECIO_BASE_MOCK
                + AMPLITUD_CICLO_MOCK * math.sin(fase)
                + random.gauss(0.0, 2.0)
            )
            env_arr[0]["Q_transado"] += abs(random.gauss(0.0, 1.0)) * 500.0
            env_arr[0]["n_ticks"] = n_ticks
            env_arr[0]["timestamp"] = ahora
            env_arr[0]["flag_dropout"] = 0
            retry_count = 0

            # --- Consumo del Ring Buffer SPSC (Sec. 7.6.2) -------------------
            # La secuencia es un contador global monótono, no un índice de slot:
            # antes se comparaba `seq_id > actuator_tail_idx`, lo que mezclaba dos
            # espacios de numeración distintos y nunca despachaba correctamente.
            # Detección de sobrepaso (overrun). El Hilo Rápido publica a ~100 Hz
            # pero el presupuesto de pesos de Binance (Sec. 8.2) solo sostiene del
            # orden de 1 orden/s, así que el productor LAPEA al consumidor y
            # sobrescribe slots no consumidos. Sin esta guarda el consumidor queda
            # esperando para siempre una secuencia que ya no existe y el lazo de
            # control se rompe en silencio (el inventario deja de actualizarse).
            slot = int(seq_esperada % RING_BUFFER_SIZE)
            seq_slot = int(act_arr[slot]["seq_id"])
            if seq_slot > seq_esperada:
                # Sec. 8.2.3: una orden vieja se purga "para evitar que el bot
                # ejecute operaciones basadas en un estado obsoleto del mercado", y
                # por el horizonte recedente (Sec. 6.4) solo el U_0 más fresco es
                # válido. Se resincroniza al ÚLTIMO publicado, no al más antiguo.
                seq_nueva = int(act_arr["seq_id"].max())
                descartadas = seq_nueva - seq_esperada
                seq_esperada = seq_nueva
                slot = int(seq_esperada % RING_BUFFER_SIZE)
                seq_slot = int(act_arr[slot]["seq_id"])
                if descartadas > 0:
                    print(
                        f"[RING] Sobrepaso SPSC: {descartadas} órdenes obsoletas "
                        f"purgadas; resincronizado en seq={seq_esperada}."
                    )

            if seq_slot == seq_esperada:
                u_c = float(act_arr[slot]["u_compra"])
                u_v = float(act_arr[slot]["u_venta"])
                seq_esperada += 1

                if (u_c > 0.0 or u_v > 0.0) and await bucket.adquirir(
                    20.0, alta_prioridad=True
                ):
                    # Sec. 4.5: órdenes IOC; el inventario se actualiza con el
                    # volumen REALMENTE ejecutado, no con el calculado.
                    fill_c = u_c * random.uniform(0.8, 1.0) if IS_TESTNET else u_c
                    fill_v = u_v * random.uniform(0.8, 1.0) if IS_TESTNET else u_v
                    inventario += fill_c - fill_v
                    seq_fill += 1
                    env_arr[0]["inv_confirmado"] = inventario
                    env_arr[0]["seq_fill"] = seq_fill

        except asyncio.CancelledError:
            raise
        except Exception:
            # Señalización de ceguera: el Hilo Rápido anulará K y dejará crecer P.
            env_arr[0]["flag_dropout"] = 1
            # τ_reintento,n = τ_base·2ⁿ + U(0, j_max)   (Sec. 8.2.3)
            base = min(60.0, 0.1 * (2**retry_count))
            await asyncio.sleep(base + random.uniform(0.0, 0.1))
            retry_count += 1


def network_engine_process(shm_env_name, shm_act_name, shm_par_name):
    print("[INIT] Proc 1 (Motor de Red) — afinidad núcleo 0")
    if hasattr(os, "sched_setaffinity"):  # No existe en Windows
        os.sched_setaffinity(0, {0})

    shm_env = shared_memory.SharedMemory(name=shm_env_name)
    shm_act = shared_memory.SharedMemory(name=shm_act_name)
    shm_par = shared_memory.SharedMemory(name=shm_par_name)
    env_arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
    act_arr = np.ndarray((RING_BUFFER_SIZE,), dtype=ACTUATOR_DTYPE, buffer=shm_act.buf)
    par_arr = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)

    try:
        asyncio.run(exchange_websocket_handler(env_arr, act_arr, par_arr))
    except KeyboardInterrupt:
        pass
    finally:
        shm_env.close()
        shm_act.close()
        shm_par.close()


# ==============================================================================
# 5. PROCESO 2: HILO RÁPIDO (EAKF + LOEPER + NMPC + TELEMETRÍA)
# ==============================================================================
def volcar_telemetria(bloque: np.ndarray, etiqueta: str) -> None:
    """Volcado a disco en hilo secundario de I/O (Sec. 8.6.2).

    NOTA DE INTERPRETACION: la Sec. 8.6.2 pide .parquet o .h5. pyarrow no está
    instalado en este entorno, así que se intenta parquet y se degrada a .npy
    (formato binario nativo, sin pérdida) dejando constancia en consola.
    """
    try:
        os.makedirs(DIR_TELEMETRIA, exist_ok=True)
        try:
            import pyarrow  # noqa: F401
            import pyarrow.parquet as pq

            tabla = pyarrow.table({n: bloque[n] for n in bloque.dtype.names})
            pq.write_table(
                tabla, os.path.join(DIR_TELEMETRIA, f"telem_{etiqueta}.parquet")
            )
        except ImportError:
            np.save(os.path.join(DIR_TELEMETRIA, f"telem_{etiqueta}.npy"), bloque)
    except Exception as err:  # El logging jamás debe tumbar el hilo de control
        print(f"[TELEMETRIA] volcado fallido ({etiqueta}): {err}")


def fast_thread_process(shm_env_name, shm_mic_name, shm_act_name, shm_par_name):
    print("[INIT] Proc 2 (EAKF + NMPC + Telemetría) — afinidad núcleo 1")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {1})

    shm_env = shared_memory.SharedMemory(name=shm_env_name)
    shm_mic = shared_memory.SharedMemory(name=shm_mic_name)
    shm_par = shared_memory.SharedMemory(name=shm_par_name)
    shm_act = shared_memory.SharedMemory(name=shm_act_name)

    env_arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
    mic_arr = np.ndarray((1,), dtype=MICELIO_DTYPE, buffer=shm_mic.buf)
    par_arr = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)
    act_arr = np.ndarray((RING_BUFFER_SIZE,), dtype=ACTUATOR_DTYPE, buffer=shm_act.buf)

    # --- Matrices del EAKF (Sec. 7.3) ---------------------------------------
    I3 = np.eye(3)
    H = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])  # Sec. 7.3.2
    P_k = np.diag([1e-5, 1e6, 1e6])  # Sec. 7.3.4: p_S≈0, p_v,p_Rn≈1e6
    S0 = float(env_arr[0]["P_spot"]) if env_arr[0]["P_spot"] > 0 else 40000.0
    x_k = np.array([[S0], [0.0], [0.0]])

    telem_buffer = np.zeros(TELEMETRY_SIZE, dtype=TELEM_DTYPE)
    telem_idx = 0
    n_volcado = 0

    is_burnt_in = False
    burnt_ticks = 0
    last_tr_P = float(np.trace(P_k))
    t_ultimo_ciclo = time.time()

    seq_emision = 1  # 0 queda reservado como "slot vacío"
    ultimo_seq_fill = 0
    inventario = 0.0
    n_ciclos = 0
    en_singularidad = False

    try:
        while True:
            n_ciclos += 1
            t_ahora = time.time()
            dt_ciclo = t_ahora - t_ultimo_ciclo
            t_ultimo_ciclo = t_ahora
            if dt_ciclo < 1e-4:
                dt_ciclo = 1e-4  # Guarda contra división por cero cinemática

            # --- Lectura lock-free del bloque estructural -------------------
            mic, _mic_ok = leer_micelio_seqlock(mic_arr)
            Omega = float(mic["Omega"])
            w_m = float(mic["w_m"])
            nu = float(mic["nu"])
            lam = float(mic["lambda_sim"])
            R_n_med = float(mic["R_n"])
            S_ref = float(mic["S_ref"])
            vol_sum_Q = float(mic["vol_sum_Q"])

            dropout = int(env_arr[0]["flag_dropout"])

            # ================= EAKF: PREDICCIÓN (Sec. 7.3.5) =================
            # Δt medido por ciclo: el 0.001 hardcodeado en A no coincidía con el
            # sleep de 0.01, y desde Colombia la señal tarda ~300 ms (Sec. 8), así
            # que ningún valor fijo es defendible.
            A = np.array([[1.0, dt_ciclo, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            x_pred = A @ x_k

            # Q_k = ρ_k·diag(q_S, q_v, q_Rn),  ρ_k = 1 + γ_ω|ω_m| + γ_Q|ΣQ|  (7.3.3)
            # ρ_k va con ω_m (frecuencia de mercado), NO con Ω (estabilidad).
            q_base = par_arr[P_Q_BASE]
            rho_k = (
                1.0
                + par_arr[P_GAMMA_OMEGA] * abs(w_m)
                + par_arr[P_GAMMA_Q] * abs(vol_sum_Q)
            )
            Q_k = rho_k * np.diag([q_base, q_base * 1e-2, q_base])
            P_pred = A @ P_k @ A.T + Q_k

            innov = np.zeros((2, 1))

            # ================= EAKF: CORRECCIÓN =============================
            # Sec. 6.5: "Si τ_d supera un umbral de tolerancia máximo τ_max, el
            # paquete se descarta por pérdida de relevancia táctica." Este umbral
            # nunca se había implementado; sin él, en el primer ciclo el timestamp
            # aún vale 0, τ_d ≈ 1.75e9 s y e^(β·τ_d) desborda a OverflowError.
            ts_paquete = float(env_arr[0]["timestamp"])
            if ts_paquete <= 0.0:
                retardo_sensor = float("inf")  # Aún no ha llegado ningún paquete
            else:
                # τ_d puede salir levemente NEGATIVO: t_ahora se muestrea al inicio
                # del ciclo y el Motor de Red publica un timestamp más reciente
                # mientras este hilo trabaja. Es granularidad de reloj, no un
                # paquete inválido; se satura a 0. Rechazarlo anulaba K, impedía
                # que P se contrajera y rompía la racha del burn-in (Sec. 7.1).
                retardo_sensor = max(0.0, t_ahora - ts_paquete)
            paquete_valido = retardo_sensor <= par_arr[P_TAU_MAX]  # Criterio Sec. 6.5

            if dropout == 1 or not paquete_valido:
                # Sec. 8.3.1: K = 0, la incertidumbre crece con Q_k (z_k = ∅).
                x_k = x_pred
                P_k = P_pred
            else:
                # R_eff^(1,1) = r_S,base · e^(β·τ_d)   (Sec. 6.5 / 7.3.3)
                R_k = np.diag(
                    [
                        par_arr[P_R_S_BASE]
                        * math.exp(par_arr[P_BETA_JITTER] * retardo_sensor),
                        par_arr[P_R_EMD],
                    ]
                )
                z_k = np.array([[float(env_arr[0]["P_spot"])], [R_n_med]])
                innov = z_k - (H @ x_pred)

                S_cov = H @ P_pred @ H.T + R_k
                # Se evita np.linalg.inv: solve es más estable si S_cov se degrada.
                K = (P_pred @ H.T) @ np.linalg.solve(S_cov, np.eye(S_cov.shape[0]))

                x_k = x_pred + (K @ innov)
                joseph = I3 - K @ H
                P_k = (joseph @ P_pred @ joseph.T) + (K @ R_k @ K.T)  # Forma de Joseph

            # ================= BURN-IN POR DERIVADA DE LA TRAZA (7.1) ========
            tr_P = float(np.trace(P_k))
            deriv_tr = abs(tr_P - last_tr_P) / dt_ciclo
            last_tr_P = tr_P

            if not is_burnt_in:
                if deriv_tr <= par_arr[P_EPS_BURN]:
                    burnt_ticks += 1
                    if burnt_ticks >= int(par_arr[P_N_BURN]):
                        is_burnt_in = True
                        print(
                            f"[BURN-IN] Traza en banda muerta tras {n_ciclos} ciclos. "
                            f"Tr(P)={tr_P:.4f}. Se libera el NMPC."
                        )
                else:
                    burnt_ticks = 0  # La ventana debe ser CONTINUA (Sec. 7.1)
                # Latido de diagnóstico: en shadow mode el bot no opera, y sin esto
                # no hay forma de distinguir "convergiendo" de "atascado".
                if n_ciclos % 200 == 0:
                    print(
                        f"[BURN-IN] ciclo={n_ciclos} Tr(P)={tr_P:.4f} "
                        f"dTr/dt={deriv_tr:.4f} (eps={par_arr[P_EPS_BURN]:.3f}) "
                        f"racha={burnt_ticks}/{int(par_arr[P_N_BURN])} "
                        f"dt={dt_ciclo*1e3:.1f}ms dropout={dropout} "
                        f"tau_d={retardo_sensor:.3f}s valido={paquete_valido}"
                    )

            # ================= TELEMETRÍA (Sec. 8.6.1) =======================
            reg = telem_buffer[telem_idx]
            reg["t_wall"] = t_ahora
            reg["tr_P"] = tr_P
            reg["x0"] = x_k[0, 0]
            reg["x1"] = x_k[1, 0]
            reg["x2"] = x_k[2, 0]
            reg["y0"] = innov[0, 0]  # Sin ỹ_k no hay Covariance Matching offline,
            reg["y1"] = innov[1, 0]  # que es el propósito declarado del logging.
            telem_idx += 1
            if telem_idx >= TELEMETRY_SIZE:
                # Se delega una COPIA a un hilo de I/O; el ciclo de control no se
                # detiene y no se pierden datos (antes el buffer se vaciaba sin
                # escribir nada a disco).
                threading.Thread(
                    target=volcar_telemetria,
                    args=(telem_buffer.copy(), f"{os.getpid()}_{n_volcado:05d}"),
                    daemon=True,
                ).start()
                n_volcado += 1
                telem_idx = 0

            # ================= SINCRONIZACIÓN DE INVENTARIO (Sec. 4.5) =======
            seq_fill = int(env_arr[0]["seq_fill"])
            if seq_fill != ultimo_seq_fill:
                inventario = float(env_arr[0]["inv_confirmado"])
                ultimo_seq_fill = seq_fill

            # ================= LOEPER (GPU) + NMPC (CPU) =====================
            if is_burnt_in and dropout == 0 and nu > 0.0 and S_ref > 0.0:
                t_solver = time.perf_counter()
                S_filtrado = float(x_k[0, 0])

                # c²_vol = k·ω_m·ν  ->  [1/Ticks]·[Ticks/Años] = [1/Años]  (Sec. 4.5)
                c2_vol = par_arr[P_K_VOL] * abs(w_m) * nu

                malla = resolver_malla_loeper(
                    S_k=S_filtrado,
                    tr_P=tr_P,
                    c2_vol=c2_vol,
                    lam=lam,
                    S_ref=S_ref,
                    par=par_arr,
                )

                u_c, u_v, freno = resolver_nmpc(
                    I_0=inventario,
                    malla=malla,
                    S_k=S_filtrado,
                    Omega=Omega,
                    lam=lam,
                    par=par_arr,
                    t_inicio=t_solver,
                )

                # Se registra la TRANSICIÓN de estado, no cada ciclo: en
                # singularidad sostenida el log inundaba la consola a 100 Hz.
                if freno != en_singularidad:
                    en_singularidad = freno
                    if freno:
                        print(
                            f"[LOEPER] Singularidad: D_min={malla['D_min']:.4f} <= eps "
                            f"(lambda={lam:.3e}/BTC, Omega={Omega:.3e}). "
                            f"Inacción total (Sec. 4.4.3)."
                        )
                    else:
                        print(
                            f"[LOEPER] Singularidad despejada: D_min={malla['D_min']:.4f}. "
                            f"Se reanuda la ejecución."
                        )

                # Cuantización comercial (Sec. 8.5) antes de publicar la orden.
                P_spot = float(env_arr[0]["P_spot"])
                u_c_q = apply_filters(u_c, 1e-5, 1e-5, 10.0, P_spot)
                u_v_q = apply_filters(u_v, 1e-5, 1e-5, 10.0, P_spot)

                if u_c_q > 0.0 or u_v_q > 0.0:
                    slot = int(seq_emision % RING_BUFFER_SIZE)
                    act_arr[slot]["u_compra"] = u_c_q
                    act_arr[slot]["u_venta"] = u_v_q
                    act_arr[slot]["ts_emision"] = t_ahora
                    act_arr[slot]["seq_id"] = seq_emision  # Publicación (último write)
                    seq_emision += 1

            time.sleep(PERIODO_HILO_RAPIDO)

    except KeyboardInterrupt:
        pass
    finally:
        if telem_idx > 0:
            volcar_telemetria(
                telem_buffer[:telem_idx].copy(), f"{os.getpid()}_{n_volcado:05d}_final"
            )
        shm_env.close()
        shm_mic.close()
        shm_par.close()
        shm_act.close()


# ==============================================================================
# 6. PROCESO 3: HILO LENTO (MICELIO, OU DE LIQUIDEZ, ESTRUCTURA)
# ==============================================================================
def slow_thread_process(shm_env_name, shm_mic_name, shm_par_name):
    print("[INIT] Proc 3 (Micelio / EMD / OU) — afinidad núcleo 2")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {2})

    shm_env = shared_memory.SharedMemory(name=shm_env_name)
    shm_mic = shared_memory.SharedMemory(name=shm_mic_name)
    shm_par = shared_memory.SharedMemory(name=shm_par_name)

    env_arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
    mic_arr = np.ndarray((1,), dtype=MICELIO_DTYPE, buffer=shm_mic.buf)
    par_arr = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)

    lam_ruido = float(par_arr[P_MU_OU])  # Arranca en la media: evita transitorio
    t_ultimo = time.time()

    # Estado del ciclo estructural vigente (Sec. 2.6 / 6.3)
    S_ref = 0.0
    Q_ref = 0.0  # Q_transado en el último nodo de fase
    ticks_ref = 0  # n_ticks en el último nodo de fase
    signo_prev = 0

    # Historia para las derivadas respecto a T̄ de la Sec. 1.4
    T_prev = 0.0
    sumQ_prev = 0.0
    w_m_prev = 0.0
    dS_prev = 0.0
    Omega = 0.0

    try:
        while True:
            t_now = time.time()
            # dt_sim medido contra el reloj de pared: con dt_sim = 0.1 fijo y un
            # sleep de 0.5 s el OU revertía a la media 5× más lento de lo calibrado.
            dt_sim = t_now - t_ultimo
            t_ultimo = t_now
            if dt_sim <= 0.0:
                dt_sim = PERIODO_HILO_LENTO

            S_actual = float(env_arr[0]["P_spot"])
            Q_bruto = float(env_arr[0]["Q_transado"])
            n_ticks = int(env_arr[0]["n_ticks"])
            if S_actual <= 0.0:
                time.sleep(PERIODO_HILO_LENTO)
                continue

            # --- ν: tasa de llegada de ticks (Sec. 4.6) ---------------------
            # NOTA DE INTERPRETACION: la Sec. 4.6 pide un buffer circular de
            # timestamps con filtro EMA. Aquí se diferencia el contador monótono
            # n_ticks del Motor de Red, que es equivalente y O(1).
            T_actual = float(n_ticks)
            nu = ((T_actual - T_prev) / dt_sim) * SEGUNDOS_POR_ANIO  # [Ticks/Años]

            # --- ω_m: frecuencia de mercado ---------------------------------
            # PENDIENTE: la cadena EMD -> Hilbert -> colapso espectral (Sec. 2.5)
            # está validada aparte (7-12 ms por llamada) pero no vive en este
            # archivo. Mientras tanto, en Testnet se inyecta un oscilador de
            # prueba. NO USAR EN MAINNET.
            w_m = 2.0e-3 * (1.0 + 0.5 * math.sin(t_now * 0.17))  # [1/Ticks]

            # --- Nodo de fase y ΔS (Sec. 2.6) -------------------------------
            # NOTA DE INTERPRETACION: el nodo real es el cruce por cero de la fase
            # desenrollada (θ ≡ 0 mód π). Sin la Transformada de Hilbert se
            # aproxima por el cambio de signo de S - S_ref.
            if S_ref <= 0.0:
                S_ref = S_actual
                Q_ref = Q_bruto
                ticks_ref = n_ticks
            delta_S = S_actual - S_ref
            signo = 1 if delta_S > 0 else (-1 if delta_S < 0 else 0)
            if signo != 0 and signo_prev != 0 and signo != signo_prev:
                S_ref = S_actual  # Nuevo nodo estructural
                Q_ref = Q_bruto  # ΣQ se reinicia a cero (Sec. 6.3)
                ticks_ref = n_ticks
                delta_S = 0.0
            if signo != 0:
                signo_prev = signo

            # ΣQ del ciclo vigente: integral por tramos reiniciada en cada nodo.
            sum_Q = Q_bruto - Q_ref

            # --- Ω: estabilidad del modelo (Sec. 1.4) -----------------------
            #   Ω = (1/ΔS)(ω_m·∂ΣQ/∂T̄ + ΣQ·∂ω_m/∂T̄) - (ΣQ·ω_m/ΔS²)·∂ΔS/∂T̄
            # Sustituye al `random.uniform(0.1,5.0)*sin(t)` que hacía de Ω ruido
            # puro y dejaba muerto todo el acoplamiento endógeno.
            dT = T_actual - T_prev
            if dT > 0.0 and abs(delta_S) > 1e-9:
                d_sumQ = (sum_Q - sumQ_prev) / dT
                d_w_m = (w_m - w_m_prev) / dT
                d_dS = (delta_S - dS_prev) / dT
                Omega = (w_m * d_sumQ + sum_Q * d_w_m) / delta_S - (
                    sum_Q * w_m / (delta_S * delta_S)
                ) * d_dS
                if not math.isfinite(Omega):
                    Omega = 0.0
            T_prev, sumQ_prev, w_m_prev, dS_prev = T_actual, sum_Q, w_m, delta_S

            # --- λ_sim: fricción sintética (Sec. 8.1.1 / 8.1.2) -------------
            if IS_TESTNET:
                theta = par_arr[P_THETA_OU]
                mu_OU = par_arr[P_MU_OU]  # Ya no está hardcodeado: viene del bloque
                sigma = par_arr[P_SIGMA_OU]
                # λ_ruido,k = λ_{k-1} + θ(μ - λ_{k-1})Δt + σ√Δt·N(0,1)
                lam_ruido += (
                    theta * (mu_OU - lam_ruido) * dt_sim
                    + sigma * math.sqrt(dt_sim) * np.random.randn()
                )
                # λ_sim = máx(λ_min, λ_ruido + η|Ω|) — el suelo evita fricciones
                # negativas computacionales; el OU con σ alto sí cruza a negativo.
                lam_sim = max(
                    par_arr[P_LAMBDA_MIN],
                    lam_ruido + par_arr[P_ETA_IMPACTO] * abs(Omega),
                )
            else:
                # Mainnet: λ = ω_m/Ψ, con Ψ = Φ/ΔS  (Sec. 4.5 / 4.6)
                Phi = sum_Q * w_m
                Psi = Phi / delta_S if abs(delta_S) > 1e-9 else 0.0
                lam_sim = (
                    max(par_arr[P_LAMBDA_MIN], w_m / Psi)
                    if abs(Psi) > 1e-12
                    else par_arr[P_LAMBDA_MIN]
                )

            # --- R_n: residuo macro de la EMD -------------------------------
            # PENDIENTE: proviene del residuo de la EMD (Sec. 2.1). Mock de Testnet.
            R_n = S_actual + math.cos(t_now * 0.1) * 50.0

            # --- Publicación atómica por seqlock (Sec. 7.6.2) ---------------
            escribir_micelio_seqlock(
                mic_arr,
                {
                    "Omega": Omega,
                    "w_m": w_m,
                    "nu": nu,
                    "lambda_sim": lam_sim,
                    "R_n": R_n,
                    "S_ref": S_ref,
                    "delta_S": delta_S,
                    "vol_sum_Q": sum_Q,
                },
            )

            time.sleep(PERIODO_HILO_LENTO)

    except KeyboardInterrupt:
        pass
    finally:
        shm_env.close()
        shm_mic.close()
        shm_par.close()


# ==============================================================================
# 7. ORQUESTADOR (Sec. 7.7)
# ==============================================================================
NOMBRES_SHM = ("shm_ent", "shm_mic", "shm_act", "shm_par")


def main():
    procesos = []
    mem_refs = []
    cerrado = threading.Event()

    def limpiar_recursos():
        """Apagado controlado (Sec. 7.7.4). Idempotente."""
        if cerrado.is_set():
            return
        cerrado.set()
        for p in procesos:
            if p.is_alive():
                p.terminate()
                p.join(timeout=1.0)
        for shm in mem_refs:
            try:
                shm.close()
                shm.unlink()
            except FileNotFoundError:
                pass
        # Barrido de seguridad por si quedó algún bloque de una corrida anterior.
        for nombre in NOMBRES_SHM:
            try:
                resto = shared_memory.SharedMemory(name=nombre)
                resto.close()
                resto.unlink()
            except FileNotFoundError:
                pass

    def unhook_graceful(signum, frame):
        print("\n>> Interrupción externa (SIG). Limpieza de /dev/shm...")
        limpiar_recursos()
        sys.exit(0)

    signal.signal(signal.SIGINT, unhook_graceful)
    signal.signal(signal.SIGTERM, unhook_graceful)

    try:
        print("== MICELIO MULTICORE DEPLOY INITIALIZING ==")
        shm_env = allocate_shared_memory("shm_ent", ENV_DTYPE.itemsize)
        shm_mic = allocate_shared_memory("shm_mic", MICELIO_DTYPE.itemsize)
        # El actuador reserva los 16 slots del Ring Buffer, no uno solo.
        shm_act = allocate_shared_memory(
            "shm_act", ACTUATOR_DTYPE.itemsize * RING_BUFFER_SIZE
        )
        shm_par = allocate_shared_memory(
            "shm_par", TOTAL_PARAMS * np.dtype(np.float64).itemsize
        )
        mem_refs.extend([shm_env, shm_mic, shm_act, shm_par])

        p_array = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)
        p_array[:] = initialize_default_parameters()

        p1 = mp.Process(
            target=network_engine_process, args=("shm_ent", "shm_act", "shm_par")
        )
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
            if any(not pr.is_alive() for pr in procesos):
                print("[SUPERVISOR] Un proceso hijo murió. Apagando el sistema.")
                break

    except KeyboardInterrupt:
        pass
    except Exception as err:
        print(f"[FATAL ORQUESTADOR]: {err}")
    finally:
        # try...finally garantiza el unlink incluso ante excepciones no-señal;
        # antes se dependía solo del signal handler y cualquier fallo filtraba
        # memoria compartida en /dev/shm.
        limpiar_recursos()


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)  # Obligatorio por el contexto CUDA (7.7.1)
    main()
