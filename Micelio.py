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

------------------------------------------------------------------------------
FASE DE CALIBRACIÓN (ORDEN_TRABAJO_CALIBRACION_1.1.md)
------------------------------------------------------------------------------
 A. SECCIÓN 0 — Las constantes de acoplamiento (γ_0, γ_ω, γ_Q, κ, μ) dejan de
    ser valores libres: se DERIVAN de límites estructurales mediante fórmulas en
    `constantes_micelio.py`. Este archivo ya no define ninguna de ellas; solo
    almacena los LÍMITES en el bloque de hot-reloading y llama a las funciones.
    Eso es lo que impide que Ω_crit, κ y μ se desincronicen (Fase 3).
 B. SECCIÓN 0.5 — Guardas dimensionales al arranque del Hilo Rápido.
 C. FASE 1.1 — NIS (ε_k = ỹᵀS⁻¹ỹ) calculado en línea y publicado en memoria
    compartida, con doble consumo: telemetría y ventana adaptativa del EMD.
 D. FASE 1.2 — El burn-in pasa a criterio NIS. Ver `DIVERGE DEL PDF (Sec. 7.1)`.
 E. FASE 1.1 — Ventana adaptativa W_k de la Sec. 2.2.1, nunca implementada.

Toda suposición que rellena un hueco del PDF va marcada con
`# NOTA DE INTERPRETACION:`, y toda decisión que contradice al PDF con
`# DIVERGE DEL PDF (Sec X.Y):`.

DEPENDENCIAS: NumPy y SciPy (`scipy.stats.chi2`, solo para las cotas del NIS, que
se evalúan UNA vez al arranque y nunca en el lazo de control). Los núcleos que el
PDF asigna a CUDA (@cuda.jit, Sec. 7.4) y a acados/CasADi (Sec. 7.5) están
implementados como referencia CPU con layout de arrays idéntico, para que el
porte sea mecánico.
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

# Las constantes de acoplamiento viven en un único sitio (Sec. 0 del orden de
# trabajo: "El resto del código importa desde ahí; nadie más las define").
import constantes_micelio as CTE

# Cadena EMD -> Hilbert (Sección 2 del PDF). Prioridad 1 de la v1.2.
import hht

# --- v1.3: riesgo de cuenta, episodios y modelo oscilatorio -------------------
import dinamica
import episodios
import mercado
import riesgo
from episodios import Estado
from mercado import Modo

# La consola de Windows usa cp1252 y no puede codificar letras griegas. Un print
# fallido lanzaría UnicodeEncodeError dentro del lazo de control y tumbaría el
# proceso, así que se degrada a reemplazo en vez de excepción.
try:
    sys.stdout.reconfigure(errors="replace")
    sys.stderr.reconfigure(errors="replace")
except (AttributeError, OSError):
    pass


def log(mensaje: str) -> None:
    """Emisión de diagnóstico que NUNCA puede tumbar el lazo de control.

    Mismo principio que el volcado de telemetría: el logging es observabilidad,
    no funcionalidad, y jamás debe propagar una excepción hacia el ciclo de
    control. Bajo `spawn` los hijos heredan el handle de stdout del padre, así que
    si el consumidor de ese pipe se cierra (un rotador de logs, un `head`, una
    terminal que muere) el siguiente print lanza BrokenPipeError y mata el
    proceso en silencio — el supervisor solo ve morir un hijo, sin causa.
    """
    try:
        print(mensaje, flush=True)
    except Exception:
        pass

# ==============================================================================
# 1. ESTRUCTURACIÓN DE BLOQUES LOCK-FREE Y C-TYPES
# ==============================================================================
SYMBOL = "BTCUSDT"

# ------------------------------------------------------------------------------
# MODO TRI-ESTADO (Sec. B.6 de ORDEN_TRABAJO_RIESGO_1_3)
# ------------------------------------------------------------------------------
# Sustituye al booleano `IS_TESTNET`, que decidía A LA VEZ tres cosas que deben
# moverse por separado: el generador de precios, el modelo de λ y el modelo de
# fills. Con un solo booleano no se puede pedir "precios reales de Mainnet pero
# sin ejecutar", que es exactamente lo que las Secciones D y E necesitan.
#
# Por defecto LECTURA: el modo que no puede tocar la cuenta. Elevar el modo tiene
# que ser un acto explícito del operador, nunca el valor por omisión.
MODO = Modo(os.environ.get("MICELIO_MODO", Modo.LECTURA.value).upper())

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
        # FASE 1.1: NIS ε_k = ỹᵀS⁻¹ỹ. Único escritor: el HILO RÁPIDO (no el Motor
        # de Red). Vive aquí y no en MICELIO_DTYPE porque el flujo es al revés que
        # el resto del bloque estructural: Rápido -> Lento. Al ser un escalar de
        # 8 bytes alineado con un solo escritor, no necesita seqlock.
        ("nis_eps", np.float64),
        # --- v1.3: estado de cuenta y de episodio, escrito por el Motor de Red ---
        # El equity viene del ACCOUNT_UPDATE del User Data Stream (Sec. B.5), NO de
        # un cálculo local: el cálculo local es justo lo que la guarda 6 existe
        # para desconfiar.
        ("equity", np.float64),  # [USD] equity de cuenta reportado por el exchange
        ("id_episodio", np.uint64),  # Sec. C.5: separa los ~30 episodios en el análisis
        ("causa_halt", np.int8),  # riesgo.CausaHalt vigente (0 = ninguna)
        # Muestras de telemetría del episodio vigente. Único escritor: el HILO
        # RÁPIDO, igual que `nis_eps`. Alimenta la compuerta de muestras mínimas
        # de la Sec. C.3, que vive en el Motor de Red: sin cruzar el dato, esa
        # compuerta no tendría forma de saber si el episodio dio para algo.
        ("muestras_episodio", np.uint64),
        # Latido del orquestador, para detectar una SEGUNDA INSTANCIA (ver
        # `verificar_instancia_unica`). Escrito por `main()` una vez por segundo.
        ("latido", np.float64),
        ("pid_orquestador", np.uint64),
        # --- v2.0 §3.4: huecos del feed. Modo de fallo NUEVO y distinto del
        # socket mudo: el feed funciona y aun asi falta informacion.
        ("n_huecos", np.uint64),
        ("trades_perdidos", np.uint64),
        ("ts_ultimo_hueco", np.float64),
        # --- v2.0 §6.3: reloj de volumen. MONOTONO y distinto de Q_transado/ΣQ,
        # que se reinicia en cada nodo de fase (Sec. 6.3 del PDF). Un reloj que
        # retrocede a cero varias veces por hora no es un reloj.
        ("Q_acumulado_total", np.float64),
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
        ("W_k", np.float64),  # Ventana adaptativa del EMD (2.2.1) [muestras]
        # --- v1.3, Sec. D.2: DOS variables de frecuencia para DOS usos ---
        # ω_m [1/Ticks] alimenta ρ_k (7.3.3) y c²_vol (4.5) y NO cambia.
        # ω_ang [rad/s] alimenta A_arm, y SOLO A_arm. Se publican por separado
        # precisamente para que nadie tenga que convertir en el punto de uso: las
        # dos trampas de la Sec. D (el 2π y el factor 125) son errores de
        # conversión en el consumidor, y ambas fallan en silencio.
        ("w_ang", np.float64),  # ω_ang = 2π·f_hz              [rad/s]
        # v2.0 §4.1: ω ANGULAR en rad/TICK, que es lo que consume A_arm bajo
        # Δn = 1. Se publica APARTE de w_ang [rad/s] y de w_m [ciclos/Tick]
        # por la regla del §2.4: cada reloj su nombre, y prohibido convertir en
        # el punto de uso — ahi es donde el 2π se ha colado ya dos veces.
        ("omega_ang_rad_tick", np.float64),
        ("C_espectral", np.float64),  # Concentración de la IMF dominante [0,1]
    ]
)

# ------------------------------------------------------------------------------
# Ring Buffer de DATOS DE MERCADO (§3.1 de ORDEN_TRABAJO_RELOJES_2_0)
# ------------------------------------------------------------------------------
# ⚠ EL CAMBIO ESTRUCTURAL DE LA v2.0. Hasta ahora `ENV_DTYPE` tenia una UNICA
# casilla escalar `P_spot`, asi que el bloque de entorno no podia transportar un
# lote: de cada sondeo de `aggTrades` —que devuelve TODAS las transacciones del
# intervalo— sobrevivia la ultima y el resto se perdia sin dejar rastro.
#
# MEDIDO EL 2026-08-04, sondeando a 1 Hz durante 61 s:
#     transacciones nuevas por lote: p50=3, p95=6, max=1000 (el limite del API)
#     tasa real de transacciones:    18.8 tx/s por REST, 26 tx/s por WebSocket
#     el bot registraba:             0.76 "mediciones"/s
#     -> FACTOR DE LOTE 24.9x, o sea que se DESCARTABA EL 96.0 % de lo recibido.
#
# Las "0.76 mediciones/s" del reporte de la v1.3 eran la tasa de PAQUETES, no la
# de informacion. El sistema no estaba escaso de datos: estaba tirandolos.
#
# Espejo del anillo de actuacion: SPSC lock-free, un productor (Motor de Red) y
# dos consumidores (Hilos Rapido y Lento).
# ⚠ DIMENSIONADO: LA TASA DE TRANSACCIONES VARIA POR UN FACTOR 20.
# Medido el mismo dia, sobre el mismo par:
#     18.8 tx/s  (REST, sondeo a 1 Hz, momento tranquilo)
#     26-32 tx/s (WebSocket @trade, momento tranquilo)
#    517.9 tx/s  (WebSocket @trade, media hora despues)
# Cualquier constante calibrada contra "la tasa tipica" es por tanto sospechosa
# de origen. Las que dependen de ν —el K de muestreo del EMD (§5.1) y el avance
# del reloj (§6.3)— se derivan de ν en ejecucion y no se fijan aqui.
# El anillo se dimensiona contra el PICO, no contra la mediana: 16384 entradas
# son ~31 s de colchon a 518 tx/s, y ~10 min en mercado tranquilo. El coste es
# 16384 x 41 bytes = 0.7 MB de memoria compartida, que es irrelevante frente a
# perder transacciones en la primera rafaga fuerte.
RING_MERCADO_SIZE = 16384
MERCADO_DTYPE = np.dtype(
    [
        # Identidad del trade. Es lo que permite DEDUPLICAR (§3.3) y DETECTAR
        # HUECOS (§3.4), y hasta la v2.0 el sistema sencillamente no la tenia.
        ("trade_id", np.uint64),
        ("precio", np.float64),  # [USD/BTC]
        ("cantidad", np.float64),  # [BTC]
        # Hora del trade EN EL EXCHANGE, no de recepcion. Usar el reloj local
        # falsearia tau_d de la Sec. 6.5, que es precisamente lo que mide.
        ("T_trade", np.float64),  # [s] epoch
        ("es_maker", np.uint8),  # lado del taker, para el desbalance de flujo
        ("seq", np.uint64),  # Secuencia global monotona; 0 = slot vacio
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
        ("y1", np.float64),  # ỹ_k componente residuo
        # FASE 1.1: se registra el ESCALAR ε_k, no la matriz S_k. Logear matrices
        # es caro e innecesario si el escalar ya se computa en línea.
        ("nis", np.float64),  # ε_k = ỹᵀS⁻¹ỹ  -> diagnóstico de consistencia
        ("dtr_dt", np.float64),  # Criterio legacy de la Sec. 7.1, para contraste
        # --- v1.3, Sec. C.5 y D.4.3 ---
        ("id_episodio", np.uint64),  # Separa los ~30 episodios en el análisis
        ("rama_A", np.int8),  # 0 = velocidad constante, 1 = armónico
        # --- v1.3, Sec. E.1: el EAKF SOMBRA ---
        # Corre sobre EL MISMO flujo de mediciones con la otra matriz A y no
        # alimenta el control. Es la única comparación honesta: dos corridas
        # distintas verían mercados distintos, y sobre un mercado real el
        # experimento no se puede repetir.
        ("y0_sombra", np.float64),
        ("y1_sombra", np.float64),
        ("nis_sombra", np.float64),
        # ⚠ IMPRESCINDIBLE DESDE QUE EL FILTRO CORRIGE POR PAQUETE NUEVO.
        # La telemetría sigue registrando una fila por CICLO DE CONTROL (~90 Hz),
        # porque tr_P, x_k y rama_A evolucionan en cada ciclo. Pero la innovación
        # solo existe cuando llegó un paquete nuevo, y en los demás ciclos vale
        # cero POR RELLENO, no por acierto del filtro.
        # Sin esta bandera, Ljung-Box y el NIS se calculan sobre una serie que es
        # ~99 % ceros: ρ₁ sale 0.0000 y el veredicto del A/B es basura con
        # apariencia de dato. Medido: 178 innovaciones reales entre 20 000 filas.
        ("hay_medicion", np.int8),
        # Entradas de la conmutación de rama (Sec. D.4). Sin ellas, `rama_A` dice
        # QUÉ rama estuvo activa pero no POR QUÉ, y un chatter no se distingue de
        # un régimen que cambia de verdad.
        ("C_espectral", np.float64),
        ("w_ang", np.float64),
    ]
)

DIR_TELEMETRIA = "telemetria"

# ⚠ VOLCADO TAMBIEN POR TIEMPO, no solo por bloque lleno.
# Desde la v2.0 hay una fila de telemetria por TRANSACCION, y la tasa de
# transacciones varia por un factor 20 con la actividad del mercado (medido:
# 26 -> 518 tx/s el mismo dia). Con volcado solo al llenar el bloque, el retardo
# hasta tener datos en disco pasa de ~20 s a ~6 min sin previo aviso, y una
# parada no limpia se lleva TODO lo acumulado desde el ultimo volcado. Eso hizo
# perder varias corridas de verificacion enteras.
PERIODO_VOLCADO_TELEMETRIA = 60.0  # [s]

# ------------------------------------------------------------------------------
# 1.1 Bloque de parámetros de HOT-RELOADING
# ------------------------------------------------------------------------------
# CLAUDE.md: "Los 11 parámetros de hot-reloading documentados en el encabezado no
# se leen en ninguna parte del programa". Se amplía el bloque a todo lo que el
# modelo necesita y se indexa por constante con nombre (los índices desnudos
# `par_arr[3]` eran la causa de que nadie los conectara).
# REORGANIZADO POR LA SECCIÓN 0 DEL ORDEN DE TRABAJO. El bloque ya NO contiene
# γ_0, γ_ω, γ_Q, κ ni μ: esas se derivan en `constantes_micelio`. Lo que se
# almacena son los LÍMITES ESTRUCTURALES de los que dependen, de modo que un
# hot-reload de un límite recalcula todas sus derivadas a la vez.

# --- Ruido del EAKF (Sec. 7.3.3) ---
P_R_S_BASE = 0  # r_S,base      varianza mínima del tick size    (7.3.3)
P_R_EMD = 1  # r_EMD         varianza histórica de Hilbert    (7.3.3)
P_Q_BASE = 2  # q_S,q_v,q_Rn  varianzas base del modelo térmico(7.3.3)

# --- Constantes de ruido CONGELADAS hasta tener datos de la cuenta demo ---
# Sec. 0.6 del orden de trabajo: NO tocar estos valores. Si una tarea parece
# exigirlo, dejar un `# TODO(calibración demo):` y seguir adelante.
P_BETA_JITTER = 3  # β             sensibilidad a la latencia       (6.5)
P_MU_OU = 4  # μ_OU          media histórica de iliquidez     (8.1.1)
P_ETA_IMPACTO = 5  # η             sensibilidad estructural         (8.1.2)
P_THETA_OU = 6  # θ_OU          tasa de reversión a la media     (8.1.1)
P_SIGMA_OU = 7  # σ_OU          volatilidad intrínseca del spread(8.1.1)
P_LAMBDA_MIN = 8  # λ_min         suelo de fricción                (8.1.2)

# --- Límites estructurales (Sec. 0.2) — de aquí se deriva todo lo demás ---
P_I_MAX = 9  # I_max         exposición máx. de inventario  [BTC]
P_K_USD = 10  # C_max         capital operativo total        [USD]
P_OMEGA_M_MAX = 11  # ω_m,max       mayor frecuencia modal HHT [1/Ticks]
P_DELTA_S_REF = 12  # ΔS_max        desplazamiento máx. vs nodo [USD/BTC]
P_OMEGA_CRIT = 13  # Ω_crit        umbral de crisis        [BTC/Ticks²]
P_SUMA_Q_MAX = 14  # max(ΣQ_T̄)     máximo histórico de volumen    [USD]
P_C_ESCALA = 15  # c             factor de escala de κ    [adimensional]
P_C_ESCALA_PRIMA = 16  # c'            factor de escala de μ    [adimensional]

# --- Adimensionales por declaración (Sec. 0.7) ---
# q_Δ, q_base y R_base son ADIMENSIONALES. Es lo que hace conmensurables a los
# tres términos de J: una vez e_k está en BTC, los tres quedan en BTC²·adim.
P_Q_DELTA = 17  # q_Δ           rigidez de seguimiento de cobertura (7.5.2)
P_Q_INV_BASE = 18  # q_base        costo asintótico de inventario      (6.2)
P_R_BASE = 19  # R_base        comisión taker estática             (6.1)

# --- Costos terminales (Sec. 7.5.2) ---
P_P_DELTA = 20  # p_Δ           costo terminal de cobertura
P_P_INV = 21  # p_inv         costo terminal de inventario

# --- Minimax-Loeper (Sec. 4.4 / 4.5 / 7.4) ---
P_K_VOL = 22  # k             c²_vol = k·ω_m·ν                 (4.5)
P_R_USD = 23  # r_USD         tasa libre de riesgo USD         (4.4)
P_Q_BTC = 24  # q_BTC         tasa de conveniencia BTC         (4.4)
P_MARGEN_MALLA = 25  # margen relativo ±% del dominio espacial   (7.4.1)
P_TAU_PRED = 26  # Δτ_pred       horizonte de predicción [AÑOS]   (4.5)
P_EPS_SING = 27  # ε             margen de seguridad de Loeper  (4.4.3)

# --- Operativos ---
P_U_MAX = 28  # U_max         riesgo máximo por operación      (6.2)
P_TAU_MAX = 29  # τ_max         tolerancia máx. de retardo [s]   (6.5)

# --- FASE 1: diagnóstico por NIS ---
P_NIS_VENTANA = 30  # n             muestras de la media móvil de ε_k
P_NIS_CONF = 31  # nivel de confianza de la banda χ² (dos colas)
P_GAMMA_NIS = 32  # γ             sensibilidad de α = e^(-γ·ε_k)  (2.2.1)
P_W_MIN = 33  # W_min         ventana mínima del EMD           (2.2.1)
P_W_MAX = 34  # W_max         ventana máxima del EMD           (2.2.1)

# --- Burn-in por derivada de la traza: criterio SECUNDARIO (Sec. 7.1) ---
# Sustituido por el criterio NIS (Fase 1.2). Se conserva para contrastarlos.
P_EPS_BURN = 35  # ε_burn        banda muerta de la traza         (7.1)
P_N_BURN = 36  # N             ticks continuos de convergencia  (7.1)

# --- Seccion C de la v1.2: varianza de proceso relativa al precio ---
P_SIGMA_REL = 37  # σ_rel         q_S = (σ_rel·S_k)²   [adimensional]

# --- v1.3 Sec. A.2: FILTROS REALES DEL INSTRUMENTO, leídos de exchangeInfo ---
# ⚠ Estos cuatro slots NO tienen valor por defecto defendible. Antes de la v1.3 el
# código llamaba `apply_filters(u_c, 1e-5, 1e-5, 10.0, P_spot)` con los cuatro
# valores inventados: un stepSize de 1e-5 contra el real de 1e-3 hace que el
# `floor` sea prácticamente la identidad, y entonces NADA en las pruebas revela
# que el controlador continuo se convierte en un interruptor al llegar al
# exchange de verdad.
# Se rellenan desde `mercado.leer_filtros` al arrancar y se refrescan en caliente
# (Sec. A.2, tarea 2): Binance cambia minNotional y stepSize sin previo aviso.
P_STEP_SIZE = 38  # stepSize del LOT_SIZE            [BTC]
P_MIN_QTY = 39  # minQty del LOT_SIZE              [BTC]
P_MIN_NOTIONAL = 40  # notional del MIN_NOTIONAL        [USDT]
P_TICK_SIZE = 41  # tickSize del PRICE_FILTER        [USD]
P_MODO = 42  # Código de mercado.Modo           [0=LECTURA,1=TESTNET,2=MAINNET]

# --- v1.3 Sec. B: capa de riesgo de cuenta ---
# NO son parámetros del modelo: no aparecen en ninguna ecuación del PDF. Viajan
# por el bloque de hot-reloading para poder endurecerlos en caliente sin reiniciar
# —endurecer un límite de riesgo nunca debe exigir parar el bot— pero se DERIVAN
# de `constantes_micelio`, que sigue siendo su única definición.
P_NOCIONAL_MAX_ORDEN = 43  # [USD] techo por orden individual        (A.3)
P_PERDIDA_MAX_EPISODIO = 44  # [USD] drawdown de equity por episodio   (B.2)
P_APALANCAMIENTO = 45  # [adim]                                  (B.3)
P_EPS_HOLGURA_POS = 46  # [adim] holgura de la abrazadera sobre I_max
P_MMR = 47  # maintenance margin rate                 (B.3)
P_N_ORDENES_MAX_MIN = 48  # guarda 3 de la Sec. B.5
P_M_RECHAZOS_MAX = 49  # guarda 4
P_T_CONFIRM_MAX = 50  # guarda 5  [s]
P_TOL_INV = 51  # guarda 6  [BTC]
P_OFFSET_RELOJ_MAX = 52  # guarda 7  [s]

# --- v1.3 Sec. D.4.2: conmutación de la rama de A por concentración espectral ---
P_C_ON = 53  # C_ON   -> entra al armónico
P_C_OFF = 54  # C_OFF  -> vuelve a velocidad constante (C_OFF < C_ON: histéresis)

# --- v2.0 §6.6: selector de reloj del sistema ---
# 0 = TICKS (transacciones, primario), 1 = VOLUMEN (cuantos de DELTA_Q*).
# Vive en el bloque de hot-reloading para poder alternarlo entre episodios sin
# reiniciar: es un experimento, no una configuracion, y la v2.0 NO decide cual
# gana — eso se decide con datos que aun no existen.
P_RELOJ = 55
P_DELTA_Q_ESTRELLA = 56  # [USD] cuanto de volumen, medido (§6.3)

TOTAL_PARAMS = 57

# Discretización de la malla de Loeper (Sec. 7.4.1).
N_S_MALLA = 61  # Impar: garantiza que S_k caiga exactamente en un nodo
N_T_MAX = 256  # Techo de pasos temporales por ciclo de control
N_NMPC = 8  # N pasos del horizonte de predicción (Sec. 5 / 6.1)
ITER_NMPC_MAX = 80  # Iteraciones del descenso proyectado
T_MAX_SOLVER = 0.020  # t_max: fallback determinista U_t = 0 (Sec. 4.5)

PERIODO_HILO_RAPIDO = 0.010  # Cadencia nominal del ciclo de control (s)
PERIODO_HILO_LENTO = 0.500  # Cadencia del análisis estructural (s)

# Período de muestreo del buffer que alimenta la EMD (Sec. 2.2). Se declara
# aparte de la cadencia del Hilo Lento para que ambas puedan divergir.
#
# DIMENSIONAMIENTO MEDIDO. La ventana del EMD está sujeta a DOS restricciones que
# tiran en direcciones opuestas:
#   (a) ciclos por ventana — con menos de ~2 períodos dentro, el tamizado no
#       aísla el modo y se lo traga el residuo (ω_m ≡ 0, ningún nodo detectado);
#   (b) muestras por ciclo — los splines cúbicos de la Sec. 2.1 necesitan
#       resolución para construir las envolventes.
# Se barrieron dt ∈ {0.25, 0.5, 1.0} s y W ∈ {128..512} sobre el ciclo de 40 s del
# generador, 12 semillas de ruido cada punto. Resultado: manda (b) con diferencia.
#     40 muestras/ciclo (dt=1.0 s)  -> error mediano 30-58 %
#     80 muestras/ciclo (dt=0.5 s)  -> error mediano  5-17 %   <- óptimo
#    160 muestras/ciclo (dt=0.25 s) -> error mediano  9-48 %   (pierde ciclos)
# Un intento previo de muestrear a 2 s "para abarcar más ciclos" degradó el error
# al 61 %: la intuición de más-ciclos-es-mejor era incorrecta.
# TODO(calibración demo): T_muestreo ≈ T_ciclo/80 contra el período real medido.
PERIODO_MUESTREO_EMD = 0.5  # [s]  -> 80 muestras por ciclo estructural de 40 s

# Parámetros del generador sintético de Testnet (solo mock; ver Sec. 8.1).
PRECIO_BASE_MOCK = 45000.0  # [USD/BTC]
AMPLITUD_CICLO_MOCK = 150.0  # [USD/BTC]
PERIODO_CICLO_MOCK = 40.0  # [s] por ciclo estructural completo


def initialize_default_parameters() -> np.ndarray:
    """Valores por defecto del bloque de hot-reloading.

    Tras la Sección 0 del orden de trabajo, este bloque contiene LÍMITES
    ESTRUCTURALES y constantes de ruido, NO constantes de acoplamiento. γ_0, γ_ω,
    γ_Q, κ y μ se derivan en `constantes_micelio` a partir de estos límites y no
    aparecen en ninguna parte como literal.
    """
    params = np.zeros(TOTAL_PARAMS, dtype=np.float64)

    # --- Ruido del EAKF (Sec. 7.3.3) ---
    params[P_R_S_BASE] = 0.5  # ~ (tick size)² del par BTCUSDT
    params[P_R_EMD] = 1.0
    params[P_Q_BASE] = 0.01

    # --- Constantes de ruido: CONGELADAS (Sec. 0.6 del orden de trabajo) ---
    # "Estas quedan deliberadamente sin revisar hasta tener datos reales de la
    # cuenta demo, que es de donde saldrá su calibración. Claude Code no debe
    # cambiar sus valores." Se conservan tal cual quedaron en la sesión anterior.
    # La guarda de estacionariedad del OU (σ/√(2θ) ≪ μ_OU) se evalúa al arranque
    # del Hilo Rápido y ADVIERTE sin corregir.
    params[P_BETA_JITTER] = 0.02
    params[P_MU_OU] = 5.0e-4  # [1/BTC]
    params[P_ETA_IMPACTO] = 5.0e-4
    params[P_THETA_OU] = 0.5
    params[P_SIGMA_OU] = 1.0e-4
    params[P_LAMBDA_MIN] = 1.0e-4

    # --- Límites estructurales (Sec. 0.2) — únicos valores libres que quedan ---
    # Se leen del módulo de constantes para que exista UNA sola definición.
    params[P_I_MAX] = CTE.I_MAX
    params[P_K_USD] = CTE.K_USD
    params[P_OMEGA_M_MAX] = CTE.OMEGA_M_MAX
    params[P_DELTA_S_REF] = CTE.DELTA_S_REF
    params[P_OMEGA_CRIT] = CTE.OMEGA_CRIT  # TODO(Fase 3): búsqueda de orden cero
    params[P_SUMA_Q_MAX] = CTE.SUMA_Q_MAX
    params[P_C_ESCALA] = CTE.C_ESCALA_KAPPA
    params[P_C_ESCALA_PRIMA] = CTE.C_ESCALA_MU

    # --- Adimensionales por declaración (Sec. 0.7) ---
    params[P_Q_DELTA] = 1.0
    params[P_Q_INV_BASE] = 0.05
    params[P_R_BASE] = 0.05

    # --- Costos terminales (Sec. 7.5.2: p_Δ ≫ q_Δ y p_inv ≫ q_base) ---
    params[P_P_DELTA] = 20.0
    params[P_P_INV] = 5.0

    # --- Minimax-Loeper (Sec. 4.4 / 4.5 / 7.4) ---
    params[P_K_VOL] = 4.0e-8  # k tal que c²_vol ~ 0.5 /Año  [CALIBRAR]
    params[P_R_USD] = 0.05  # 1/Años
    params[P_Q_BTC] = 0.00  # 1/Años
    params[P_MARGEN_MALLA] = CTE.MARGEN_MALLA_REL
    params[P_TAU_PRED] = 30.0 / SEGUNDOS_POR_ANIO  # 30 s expresados en AÑOS
    params[P_EPS_SING] = 1.0e-3

    # --- Operativos (Sec. 6.2) ---
    params[P_U_MAX] = 0.05  # BTC por operación
    # Sec. 8: desde Colombia la señal tarda ~300 ms hasta Tokio/Europa.
    params[P_TAU_MAX] = 2.0  # [s]

    # --- FASE 1: criterio de consistencia por NIS ---
    # El NIS es adimensional y auto-normalizado: su distribución de referencia
    # (χ² con m g.l.) no depende de Δt, ni de ρ_k, ni del régimen de volatilidad.
    params[P_NIS_VENTANA] = 100.0  # muestras de la media móvil
    params[P_NIS_CONF] = 0.95  # banda de dos colas al 95 %
    # γ de α = e^(-γ·ε_k) (Sec. 2.2.1). Con ε_k ~ m = 2 en régimen consistente,
    # γ = 0.35 deja α ≈ e^(-0.7) ≈ 0.50: la ventana se sitúa a mitad de rango
    # cuando el filtro está sano, con margen para expandirse y contraerse.
    params[P_GAMMA_NIS] = 0.35
    # W_min y W_max fijados por el barrido documentado en PERIODO_MUESTREO_EMD:
    # a 0.5 s de muestreo, W=192 abarca 2.4 ciclos estructurales (error mediano
    # 4.9 %) y W=384 abarca 4.8 (p90 del 20 %, el mejor de la tabla). Latencia
    # medida 6-10 ms por llamada, contra los 500 ms del Hilo Lento.
    params[P_W_MIN] = 192.0
    params[P_W_MAX] = 384.0

    # --- Burn-in por derivada de la traza: criterio SECUNDARIO ---
    # DIVERGE DEL PDF (Sec. 7.1): el criterio primario pasa a ser el NIS; ver la
    # justificación completa en `burnin_por_nis`. Estos dos parámetros se
    # conservan únicamente para registrar el criterio viejo en telemetría y poder
    # contrastar ambos durante un tiempo, como pide la Fase 1.2.
    params[P_EPS_BURN] = 1.0
    params[P_N_BURN] = 100.0

    # --- Seccion C de la v1.2: varianza de proceso relativa ---
    params[P_SIGMA_REL] = CTE.SIGMA_REL

    # --- v1.3 Sec. A.2: filtros del instrumento ---
    # Se dejan en CERO a propósito. No es un valor por defecto: es un centinela.
    # `cargar_filtros_en_parametros` los rellena leyendo `exchangeInfo`, y las
    # guardas de arranque abortan si siguen en cero — así, olvidarse de leerlos
    # falla ruidosamente en vez de operar con un literal inventado.
    params[P_STEP_SIZE] = 0.0
    params[P_MIN_QTY] = 0.0
    params[P_MIN_NOTIONAL] = 0.0
    params[P_TICK_SIZE] = 0.0
    params[P_MODO] = mercado.CODIGO_MODO[MODO]

    # --- v1.3 Sec. B: capa de riesgo de cuenta ---
    params[P_NOCIONAL_MAX_ORDEN] = CTE.NOCIONAL_MAX_ORDEN
    params[P_PERDIDA_MAX_EPISODIO] = CTE.PERDIDA_MAX_EPISODIO
    params[P_APALANCAMIENTO] = float(CTE.APALANCAMIENTO)
    params[P_EPS_HOLGURA_POS] = CTE.EPS_HOLGURA_POSICION
    params[P_MMR] = mercado.MMR_PRIMER_TRAMO_BTCUSDT * mercado.FACTOR_SEGURIDAD_MMR
    params[P_N_ORDENES_MAX_MIN] = float(CTE.N_ORDENES_MAX_MIN)
    params[P_M_RECHAZOS_MAX] = float(CTE.M_RECHAZOS_MAX)
    params[P_T_CONFIRM_MAX] = CTE.T_CONFIRM_MAX
    params[P_TOL_INV] = CTE.TOL_INV
    params[P_OFFSET_RELOJ_MAX] = CTE.OFFSET_RELOJ_MAX

    # --- v1.3 Sec. D.4.2: histéresis de la rama de A ---
    params[P_C_ON] = CTE.C_ON_ARMONICO
    params[P_C_OFF] = CTE.C_OFF_ARMONICO

    # --- v2.0 §6.6: reloj vigente ---
    params[P_RELOJ] = CTE.CODIGO_RELOJ[
        CTE.reloj_desde_codigo(
            CTE.CODIGO_RELOJ.get(os.environ.get("MICELIO_RELOJ", CTE.RELOJ_TICKS), 0.0)
        )
    ]
    params[P_DELTA_Q_ESTRELLA] = CTE.DELTA_Q_ESTRELLA
    return params


def cargar_filtros_en_parametros(par_arr, modo=None, verboso: bool = True):
    """Sec. A.2: lee `exchangeInfo` y publica los filtros en el hot-reloading.

    Se leen AMBOS entornos aunque solo se opere en uno. La comparación es lo que
    delata la discrepancia de granularidad medida el 2026-08-04: Testnet usa
    stepSize = 0.0001 y Mainnet 0.001, o sea que **Testnet es 10× más fino**.
    Calibrar el dimensionamiento contra Testnet produciría un sistema que pasa la
    guarda de resolución en pruebas (476 lotes) y se degrada a un interruptor en
    producción (47 lotes) — exactamente el fallo que la Sección A previene.

    Devuelve (filtros_del_modo, filtros_de_mainnet). El segundo es el que la
    guarda de resolución debe usar, sea cual sea el modo vigente.
    """
    modo = modo if modo is not None else MODO
    lecturas = mercado.leer_filtros_ambos_entornos(SYMBOL)
    for entorno, resultado in lecturas.items():
        if isinstance(resultado, Exception):
            log(f"[EXCHANGE] Fallo leyendo exchangeInfo de {entorno.value}: {resultado}")

    f_mainnet = lecturas.get(Modo.MAINNET)
    if isinstance(f_mainnet, Exception) or f_mainnet is None:
        raise RuntimeError(
            "No se pudo leer exchangeInfo de MAINNET. La guarda de resolucion de "
            "control (Sec. A.3) DEBE evaluarse contra el stepSize de Mainnet, y no "
            "hay valor por defecto defendible para el. Sin red, no se arranca."
        )

    # LECTURA comparte instrumento con MAINNET: es el mismo mercado.
    clave = Modo.MAINNET if modo in (Modo.LECTURA, Modo.MAINNET) else Modo.TESTNET
    filtros = lecturas.get(clave)
    if isinstance(filtros, Exception) or filtros is None:
        raise RuntimeError(f"No se pudo leer exchangeInfo de {clave.value}: {filtros}")

    par_arr[P_STEP_SIZE] = filtros.step_size
    par_arr[P_MIN_QTY] = filtros.min_qty
    par_arr[P_MIN_NOTIONAL] = filtros.min_notional
    par_arr[P_TICK_SIZE] = filtros.tick_size

    if verboso:
        for resultado in lecturas.values():
            if not isinstance(resultado, Exception):
                log("[EXCHANGE] " + resultado.resumen(CTE.PRECIO_REFERENCIA))
    return filtros, f_mainnet


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
        # Residuo de una ejecución anterior mal terminada.
        #
        # ⚠ DIFERENCIA DE PLATAFORMA. En POSIX `unlink()` borra la entrada de
        # /dev/shm y un `create=True` posterior funciona. En Windows NO existe
        # unlink real: el bloque vive mientras algún proceso lo tenga mapeado, así
        # que si quedó un huérfano vivo, destruir y recrear falla otra vez con
        # WinError 183. La ruta anterior hacía exactamente eso y reventaba el
        # arranque tras una caída sucia.
        #
        # Estrategia robusta en ambos: ADJUNTARSE al bloque existente y reutilizarlo
        # si tiene tamaño suficiente. Se pone a cero igual, así que el estado previo
        # no se hereda.
        previo = shared_memory.SharedMemory(name=nombre)
        if previo.size >= n_bytes:
            log(
                f"[SHM] Reciclando bloque huerfano '{nombre}' "
                f"({previo.size} bytes disponibles, se necesitan {n_bytes})."
            )
            shm = previo
        else:
            # Demasiado pequeño para reutilizarlo: no queda más que destruirlo.
            previo.close()
            previo.unlink()
            shm = shared_memory.SharedMemory(name=nombre, create=True, size=n_bytes)
    shm.buf[:n_bytes] = b"\x00" * n_bytes
    return shm


PERIODO_LATIDO = 1.0  # [s] cadencia con que el orquestador marca que sigue vivo
TOLERANCIA_LATIDO = 5.0  # [s] sin latido para dar por muerta a la otra instancia


class ErrorSegundaInstancia(RuntimeError):
    """Ya hay un Micelio vivo sobre estos mismos bloques de memoria compartida."""


def verificar_instancia_unica(nombre_shm: str = "shm_ent") -> None:
    """Impide que dos bots compartan silenciosamente la memoria. v1.3.

    ⚠ DEFECTO REAL, ENCONTRADO EN EJECUCIÓN EL 2026-08-04. La recuperación de
    bloques huérfanos que la v1.2 añadió para Windows (donde `unlink()` es un
    no-op) tiene un efecto de segundo orden que no se vio entonces: si el bloque
    existe porque **hay otro Micelio VIVO**, la segunda instancia se ADJUNTA a él
    en vez de fallar. Entonces dos Hilos Lentos escriben el mismo `shm_mic` bajo
    el mismo seqlock, dos Hilos Rápidos publican en el mismo Ring Buffer, y los
    tres procesos de cada bot leen una mezcla de ambos mundos.

    Se manifestó como CHATTER en la conmutación de la rama de A: 294 cambios en
    10 000 ciclos, permanencia mediana de 0.27 s. Con una sola instancia, sobre el
    mismo mercado y la misma configuración: **12 cambios en 20 000 ciclos**. El
    síntoma apuntaba a la histéresis de la Sec. D.4.2 —que estaba bien— y no a la
    memoria compartida por ningún lado.

    Detección por latido y no por PID: un PID puede reciclarse, y en Windows no
    hay forma barata y portable de preguntar si un PID concreto sigue siendo el
    mismo proceso. Un timestamp que se refresca cada segundo no tiene esa
    ambigüedad — o está fresco, o no lo está.
    """
    try:
        previo = shared_memory.SharedMemory(name=nombre_shm)
    except FileNotFoundError:
        return  # No hay bloque: camino limpio
    try:
        if previo.size < ENV_DTYPE.itemsize:
            return  # Bloque de otro esquema; `allocate_shared_memory` lo maneja
        arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=previo.buf)
        latido = float(arr[0]["latido"])
        pid = int(arr[0]["pid_orquestador"])
        edad = time.time() - latido
        if latido > 0.0 and edad < TOLERANCIA_LATIDO:
            raise ErrorSegundaInstancia(
                f"Ya hay un Micelio vivo (pid={pid}, ultimo latido hace "
                f"{edad:.1f} s) usando '{nombre_shm}'. Arrancar un segundo bot "
                f"haria que ambos escribieran la MISMA memoria compartida y los "
                f"datos de los dos serian basura silenciosa. Cierra el otro "
                f"proceso, o espera {TOLERANCIA_LATIDO:.0f} s si ya murio."
            )
    finally:
        previo.close()


class ConsumidorMercado:
    """Lector SPSC del anillo de mercado, con deteccion de sobrepaso. §3.

    Cada consumidor (Hilo Rapido, Hilo Lento) mantiene su propia secuencia, de
    modo que van a su ritmo sin bloquearse entre si — que es la razon de que el
    anillo tenga un solo productor y varios lectores independientes.

    ⚠ EL SOBREPASO NO SE PUEDE IGNORAR. El productor publica a ~26 tx/s en
    rafagas de hasta 203; si un consumidor se retrasa mas de RING_MERCADO_SIZE
    entradas, el productor lo lapea y esas transacciones se pierden. Eso es
    indistinguible de un hueco del exchange (§3.4) salvo porque este es CULPA
    NUESTRA, y por eso se cuenta aparte: confundirlos haria buscar el problema en
    Binance.
    """

    def __init__(self, mer_arr, nombre: str = ""):
        self.mer_arr = mer_arr
        self.nombre = nombre
        self.seq_esperada = 1  # 0 queda reservado como "slot vacio"
        self.perdidos_por_sobrepaso = 0
        self.n_sobrepasos = 0

    def leer_lote(self, maximo: int = RING_MERCADO_SIZE):
        """Devuelve la lista de trades nuevos, en orden. Puede venir vacia."""
        lote = []
        for _ in range(maximo):
            slot = int(self.seq_esperada % RING_MERCADO_SIZE)
            seq_slot = int(self.mer_arr[slot]["seq"])
            if seq_slot == self.seq_esperada:
                lote.append(self.mer_arr[slot].copy())
                self.seq_esperada += 1
                continue
            if seq_slot > self.seq_esperada:
                # Lapeado: se resincroniza al mas antiguo AUN VIVO, no al mas
                # nuevo. Aqui, a diferencia del anillo de actuacion, los datos
                # viejos SI valen: son mediciones de mercado, no ordenes que
                # caducan por horizonte recedente (Sec. 8.2.3).
                seq_max = int(self.mer_arr["seq"].max())
                seq_min_viva = max(1, seq_max - RING_MERCADO_SIZE + 1)
                self.perdidos_por_sobrepaso += seq_min_viva - self.seq_esperada
                self.n_sobrepasos += 1
                self.seq_esperada = seq_min_viva
                continue
            break  # seq_slot < esperada -> slot aun sin escribir: no hay mas
        return lote


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

    ⚠ v1.3: esta función ya NO es el punto de cuantización del sistema. La cadena
    de la Sec. B.4 pasa por `riesgo.CapaRiesgo.preparar_orden`, que aplica los
    clamps de nocional y posición ANTES del floor y los revalida después. Se
    conserva porque el Hilo Rápido cuantiza para decidir si vale la pena publicar
    en el Ring Buffer, pero los cuatro filtros llegan ahora del bloque de
    hot-reloading (leídos de `exchangeInfo`), no como literales.
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
# 1.3 CONSISTENCIA DEL FILTRO POR NIS (Fase 1.1 y 1.2 del orden de trabajo)
# ==============================================================================
# La Secuencia de Innovación Normalizada al Cuadrado es
#
#     ε_k = ỹ_kᵀ S_k⁻¹ ỹ_k
#
# y ESTO YA ESTÁ EN EL PDF: la Sec. 2.2.1 lo define con ese mismo símbolo y esa
# misma fórmula para modular la ventana adaptativa del EMD. Es el mismo escalar
# con dos usos — diagnóstico de consistencia del filtro Y driver de la ventana —
# que hasta ahora no se calculaba para ninguno de los dos.
DIM_MEDICION = 2  # m = rango del vector de medición z_k (Sec. 7.3.1)


def nis_escalar(innov: np.ndarray, S_cov: np.ndarray) -> float:
    """ε_k = ỹᵀ S⁻¹ ỹ. Se resuelve el sistema, no se invierte S."""
    try:
        return float(innov.T @ np.linalg.solve(S_cov, innov))
    except np.linalg.LinAlgError:
        return float("nan")  # S_cov singular: el ciclo lo trata como no informativo


def cotas_nis(n_ventana: int, confianza: float, m: int = DIM_MEDICION):
    """Banda de dos colas para la MEDIA MÓVIL de ε_k sobre `n_ventana` muestras.

    Si el filtro es consistente, ε_k ~ χ²_m, luego E[ε_k] = m = 2. La suma de n
    muestras independientes es χ²_{n·m}, así que la media móvil tiene por banda

        [ χ²_{n·m}(α/2) / n ,  χ²_{n·m}(1−α/2) / n ]

    ⚠ La independencia entre muestras es exactamente lo que testea Ljung-Box
    (Fase 1.3). Si la innovación está autocorrelacionada, esta banda es optimista
    y el burn-in cerrará antes de tiempo — otra razón por la que Ljung-Box es la
    compuerta de entrada a la Fase 2.

    Se evalúa UNA vez al arranque de cada proceso: `scipy` nunca entra al lazo de
    control. Por eso el import es perezoso.
    """
    from scipy.stats import chi2  # Import perezoso: fuera del lazo de control

    alfa = 1.0 - confianza
    gl = n_ventana * m
    return (
        float(chi2.ppf(alfa / 2.0, gl) / n_ventana),
        float(chi2.ppf(1.0 - alfa / 2.0, gl) / n_ventana),
    )


def ventana_adaptativa(eps_k: float, gamma: float, w_min: float, w_max: float) -> int:
    """W_k = ⌊W_min + α·(W_max − W_min)⌉  con  α = e^(−γ·ε_k).   Sec. 2.2.1.

    Nunca se había implementado: el EMD usaba ventana fija. Comportamiento que
    describe el PDF:
      - Mercado estable (ε_k → 0) ⟹ α → 1 ⟹ W_k → W_max: más resolución para
        detectar ciclos de baja frecuencia.
      - Salto abrupto (ε_k ≫ 0) ⟹ α → 0 ⟹ W_k → W_min: se descarta la memoria
        pasada para que el offset histórico no corrompa el nuevo régimen.
    """
    if not math.isfinite(eps_k) or eps_k < 0.0:
        eps_k = 0.0
    alfa = math.exp(-gamma * eps_k)
    return int(round(w_min + alfa * (w_max - w_min)))


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
    margen_rel = par[P_MARGEN_MALLA]
    r_usd = par[P_R_USD]
    q_btc = par[P_Q_BTC]
    eps_sing = par[P_EPS_SING]
    tau_pred = par[P_TAU_PRED]

    # γ_0 = I_max/(S·ΔS_max)  (Sec. 0.1 del orden de trabajo). DERIVADA, no leída:
    # depende de S, así que se reevalúa cada ciclo con el precio filtrado. Fijarla
    # como literal es precisamente el error que la fórmula elimina — la versión con
    # unidades de Delta daba λS²Γ ≈ 25 y singularidad permanente desde el arranque.
    gamma_0 = CTE.gamma_0(S_k, par[P_I_MAX], par[P_DELTA_S_REF])

    # --- Eje de precio (Sec. 7.4.1) -----------------------------------------
    S_min = S_k * (1.0 - margen_rel)
    S_sup = S_k * (1.0 + margen_rel)
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
    p_delta = par[P_P_DELTA]
    p_inv = par[P_P_INV]
    R_base = par[P_R_BASE]
    eps_sing = par[P_EPS_SING]

    # κ y μ DERIVADOS de Ω_crit (Sec. 0.1). Es crítico que se evalúen aquí y no se
    # lean de un slot propio: la Fase 3 calibra Ω_crit por búsqueda de orden cero,
    # y "los tres deben recalcularse juntos — Ω_crit no puede quedar como valor
    # suelto en un sitio y κ hardcodeado en otro". Al derivarlos, mover Ω_crit
    # mueve a la vez la pendiente de R(Ω) y la de q_inv(Ω), por construcción.
    omega_crit = par[P_OMEGA_CRIT]
    kappa = CTE.kappa(R_base, omega_crit, par[P_C_ESCALA])
    mu_inv = CTE.mu_inventario(par[P_Q_INV_BASE], omega_crit, par[P_C_ESCALA_PRIMA])
    q_inv = par[P_Q_INV_BASE] + mu_inv * Omega * Omega  # q_inv(Ω), Sec. 6.2

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
    # El cap de capital de la Sec. 6.2 está en BTC, luego es I_max (ver la nota
    # sobre la colisión de nombres C_max en `constantes_micelio`).
    I_max = par[P_I_MAX]
    if abs(Omega) < omega_crit:
        ub_c = min(I_max, U_max)
        ub_v = min(I_max, U_max)
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


async def ingesta_mercado(env_arr, modo, estado):
    """Alimenta el bloque de entorno. Modo LECTURA = feed público REAL de Mainnet.

    PRECONDICIÓN de la v1.3: las Secciones D y E se calibran sobre datos reales de
    Mainnet, no sobre Testnet (libro simulado) ni sobre mocks (los pusimos
    nosotros, así que no prueban nada sobre el mercado).

    En LECTURA y MAINNET se consume `btcusdt@aggTrade`, que trae precio Y volumen
    por trade — el volumen es lo que la Sec. 1.1 llama ΣQ y lo que el mock tenía
    que inventarse.
    """
    if modo is Modo.TESTNET:
        await _ingesta_mock(env_arr, estado)
        return

    # El feed se CONSERVA entre reconexiones del watchdog. Si se reconstruyera en
    # cada reintento, el contador de estancamientos volvería a cero y el bot
    # reintentaría eternamente un socket que ya se demostró mudo, sin llegar
    # nunca a degradar al sondeo REST.
    feed = estado.get("feed")
    if feed is None:
        feed = mercado.FeedPublico(modo, SYMBOL)
        estado["feed"] = feed

    await feed.escuchar(lambda tick: publicar_trade(env_arr, estado, tick))


def publicar_trade(env_arr, estado, tick) -> None:
    """Publica UNA transaccion en el anillo de mercado. §3.2, §3.3 y §3.4.

    Reemplaza al `env_arr[0]["P_spot"] = tick.precio` de la v1.3, que sobrescribia
    la casilla escalar y perdia el 96 % del lote.
    """
    mer_arr = estado["mer_arr"]

    # --- §3.3 DEDUPLICACION OBLIGATORIA ---------------------------------------
    # El sondeo REST devuelve ventanas SOLAPADAS: sin deduplicar, los mismos
    # trades entran al filtro varias veces. Eso es exactamente el bug de las 90
    # correcciones con otro disfraz, y seria la tercera aparicion de la misma
    # familia en tres sesiones.
    # El contador NO se reinicia entre reconexiones, a proposito: reiniciarlo
    # reintroduciria el duplicado justo despues de cada caida del feed, que es
    # cuando menos se mira.
    tid = int(tick.trade_id)
    if tid != 0 and tid <= estado["ultimo_trade_id"]:
        estado["n_duplicados"] += 1
        return

    # --- §3.4 DETECCION DE HUECOS ---------------------------------------------
    # Si el id salta mas de 1, hubo transacciones que NO se vieron. Es un modo de
    # fallo distinto del socket mudo: el feed funciona y aun asi falta
    # informacion, y hasta la v2.0 el sistema no podia detectarlo porque no tenia
    # identidad de trade.
    if tid != 0 and estado["ultimo_trade_id"] > 0:
        salto = tid - estado["ultimo_trade_id"]
        if salto > 1:
            estado["trades_perdidos"] += salto - 1
            estado["n_huecos"] += 1
            env_arr[0]["n_huecos"] = estado["n_huecos"]
            env_arr[0]["trades_perdidos"] = estado["trades_perdidos"]
            # Marca el instante para que el analisis offline pueda EXCLUIR el
            # tramo: un rho_1 calculado sobre una serie con huecos silenciosos es
            # basura con apariencia de dato, misma razon que `hay_medicion`.
            env_arr[0]["ts_ultimo_hueco"] = time.time()
    if tid != 0:
        estado["ultimo_trade_id"] = tid

    # --- Publicacion en el anillo (ultimo write: `seq`) -----------------------
    seq = estado["seq_mercado"]
    slot = int(seq % RING_MERCADO_SIZE)
    mer_arr[slot]["trade_id"] = tid
    mer_arr[slot]["precio"] = tick.precio
    mer_arr[slot]["cantidad"] = tick.cantidad
    mer_arr[slot]["T_trade"] = tick.ts_evento
    mer_arr[slot]["es_maker"] = 1 if getattr(tick, "es_maker", False) else 0
    mer_arr[slot]["seq"] = seq  # publicacion: debe ser la ULTIMA escritura
    estado["seq_mercado"] = seq + 1

    # --- §3.2 n_ticks incrementa POR TRANSACCION, no por paquete --------------
    # Con esto ν pasa a ser la tasa real de transacciones y se corrigen a la vez
    # todas las filas de la tabla del §1.2: omega_m dejaba de estar sobreestimada
    # por el factor de lote, rho_k dejaba de inflar Q por esa via, y el piso de
    # Nyquist deja de aplicarse sobre el tick equivocado.
    estado["n_ticks"] += 1
    env_arr[0]["n_ticks"] = estado["n_ticks"]

    # Se conserva P_spot como "ultimo precio conocido" para los consumidores de
    # reloj de pared (guardas de riesgo, clamps, NMPC), que no quieren un lote.
    env_arr[0]["P_spot"] = tick.precio
    # ΣQ en USD: la Sec. 1.1 lo define como volumen transado, y el feed lo
    # entrega en BTC, así que se convierte con el precio del propio trade.
    dq_usd = tick.cantidad * tick.precio
    env_arr[0]["Q_transado"] += dq_usd
    # §6.3: acumulador MONOTONO de volumen, DISTINTO de ΣQ_T̄. Este no se
    # reinicia en los nodos de fase, porque un reloj debe ser monotono y ΣQ
    # retrocede a cero varias veces por hora.
    env_arr[0]["Q_acumulado_total"] += dq_usd
    # Reloj del EXCHANGE, no el local: τ_d de la Sec. 6.5 es precisamente la
    # latencia entre ambos, y usar el local la haría idénticamente cero.
    env_arr[0]["timestamp"] = tick.ts_evento
    env_arr[0]["flag_dropout"] = 0
    estado["ultimo_precio"] = tick.precio


async def _ingesta_mock(env_arr, estado):
    """Generador sintético. SOLO para TESTNET mientras no haya credenciales.

    Serie con un modo cíclico dominante más microestructura. NO es ruido blanco a
    propósito: el modelo de la Sec. 1 supone un mercado que revierte a nodos de
    precio, y el detector de nodos de fase (Sec. 2.6) necesita una señal con ciclo
    real. Con ruido blanco puro el nodo dispara en casi todos los ticks, ΣQ se
    reinicia sin parar y el burn-in nunca alcanza la banda muerta.

    ⚠ NO USAR EN MAINNET, y tampoco para concluir nada sobre el modelo: la
    Sec. E.2 exige Modo LECTURA sobre Mainnet porque una comparación de modelos de
    precio contra una serie que fabricamos nosotros no significa nada.
    """
    t_origen = time.time()
    while True:
        await asyncio.sleep(0.005)
        ahora = time.time()
        estado["n_ticks"] += 1
        fase = 2.0 * math.pi * (ahora - t_origen) / PERIODO_CICLO_MOCK
        precio = (
            PRECIO_BASE_MOCK
            + AMPLITUD_CICLO_MOCK * math.sin(fase)
            + random.gauss(0.0, 2.0)
        )
        env_arr[0]["P_spot"] = precio
        env_arr[0]["Q_transado"] += abs(random.gauss(0.0, 1.0)) * 500.0
        env_arr[0]["n_ticks"] = estado["n_ticks"]
        env_arr[0]["timestamp"] = ahora
        env_arr[0]["flag_dropout"] = 0
        estado["ultimo_precio"] = precio


async def vigilancia_reloj(capa, modo, periodo: float = 30.0):
    """Guarda 7 de la Sec. B.5: |offset| contra `/fapi/v1/time`.

    Corre aparte del lazo de actuación porque su cadencia es tres órdenes de
    magnitud más lenta y porque un fallo de red aquí no debe frenar las órdenes:
    un sondeo fallido deja el último offset conocido, que es información vieja
    pero no falsa.
    """
    while True:
        try:
            offset = await mercado.offset_reloj(modo)
            capa.registrar_offset_reloj(offset)
        except Exception:
            pass
        await asyncio.sleep(periodo)


async def lazo_actuacion(env_arr, act_arr, par_arr, modo, capa, cuenta, estado):
    """Consumo del Ring Buffer y CADENA DE LA SEC. B.4 antes de firmar.

    Este es el único punto del sistema por el que una orden puede salir, y vive en
    el Motor de Red a propósito (Sec. B.4): la capa de riesgo debe poder detener el
    sistema SIN COOPERACIÓN del Hilo Rápido. Si viviera dentro del NMPC no sería
    una capa de seguridad, sería parte de lo que debe vigilar.
    """
    bucket = TokenBucket(capacidad=1200.0, ventana_seg=60.0)
    seq_esperada = 1  # Las secuencias emitidas arrancan en 1 (0 = slot vacío)
    seq_fill = 0
    n_orden = 0

    while True:
        await asyncio.sleep(0.005)
        S = float(env_arr[0]["P_spot"])
        if S <= 0.0:
            continue

        # La cuenta de papel se dimensiona con el PRIMER PRECIO REAL, no con un
        # ancla de referencia: EQUITY_MIN_EPISODIO va con S, y arrancar con un
        # precio inventado dejaría el episodio mal capitalizado desde el minuto
        # cero. Solo aplica a LECTURA (Sec. B.5).
        if cuenta is None and estado.get("modo_lectura"):
            cuenta = riesgo.CuentaPapel(CTE.equity_min_episodio(S) * 1.10)
            estado["cuenta"] = cuenta
            capa.registrar_equity(cuenta.equity_inicial)
            log(
                f"[CUENTA] Modo LECTURA: contabilidad de PAPEL sobre precios "
                f"reales, equity inicial {cuenta.equity_inicial:.0f} USD a "
                f"S={S:.2f}. NO es ACCOUNT_UPDATE."
            )

        # --- Equity y evaluación de las siete guardas (Sec. B.5) -------------
        if cuenta is not None:
            capa.registrar_equity(cuenta.equity(S))
            env_arr[0]["equity"] = capa.equity_actual
        causa = capa.evaluar(S)
        env_arr[0]["causa_halt"] = int(causa)
        if causa != riesgo.CausaHalt.NINGUNA:
            if not estado.get("halt_reportado"):
                estado["halt_reportado"] = True
                log(
                    f"[HALT] causa={causa.name} ({riesgo.DESCRIPCION_CAUSA[causa]}): "
                    f"{capa.detalle_halt}"
                )
                estado["al_disparar"](causa, capa.detalle_halt)
            # Semántica del halt (Sec. B.5): CERRAR Y PARAR, no congelar.
            # Congelar deja exposición abierta sin supervisión.
            await asyncio.sleep(0.1)
            continue

        # --- Consumo del Ring Buffer SPSC (Sec. 7.6.2) -----------------------
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
                log(
                    f"[RING] Sobrepaso SPSC: {descartadas} órdenes obsoletas "
                    f"purgadas; resincronizado en seq={seq_esperada}."
                )

        if seq_slot != seq_esperada:
            continue

        # --- Paso 1 de la Sec. B.4: leer u del Ring Buffer -------------------
        u_c_bruto = float(act_arr[slot]["u_compra"])
        u_v_bruto = float(act_arr[slot]["u_venta"])
        seq_esperada += 1

        # --- Pasos 2 a 6: clamps, cuantización y re-validación ---------------
        u_c, u_v, motivo = capa.preparar_orden(
            u_c_bruto, u_v_bruto, S, capa.inventario_local
        )
        if motivo is riesgo.MotivoNoEnvio.BAJO_RESOLUCION and (
            capa.n_bajo_resolucion % 100 == 1
        ):
            # Una racha de estos es la firma de un tope por orden mal
            # dimensionado (Sec. A.3), no de un mercado tranquilo.
            log(
                f"[RESOLUCION] {capa.n_bajo_resolucion} ordenes anuladas por el "
                f"floor: u_c={u_c_bruto:.6f} u_v={u_v_bruto:.6f} BTC contra "
                f"stepSize={capa.filtros.step_size:g} y minNotional="
                f"{capa.filtros.min_notional:g} USDT a S={S:.2f}."
            )
        if motivo is not riesgo.MotivoNoEnvio.ENVIADA:
            continue

        if not await bucket.adquirir(20.0, alta_prioridad=True):
            continue

        # --- Paso 7: firmar y enviar -----------------------------------------
        n_orden += 1
        id_cliente = f"mic-{estado['id_episodio']:03d}-{n_orden:06d}"
        capa.registrar_envio(id_cliente)

        if mercado.ejecucion_permitida(modo):
            # Aquí va la petición firmada. `assert_ejecucion_permitida` es
            # redundante con el `if` a propósito: es la compuerta dura que hace
            # IMPOSIBLE —no improbable— que una orden salga en Modo LECTURA, y
            # debe seguir en pie aunque alguien reordene este bloque.
            mercado.assert_ejecucion_permitida(modo, "enviar orden")
            # TODO(credenciales): POST /fapi/v1/order firmado con HMAC. Hasta
            # entonces TESTNET simula el fill igual que LECTURA, y por eso las
            # estadisticas de fill de Testnet NO son evidencia sobre calidad de
            # ejecucion (Sec. B.6): el libro es delgado y erratico, y Testnet
            # valida PLOMERIA -- firma, filtros, reconciliacion, kill switch.
            fill_c = u_c * random.uniform(0.8, 1.0)
            fill_v = u_v * random.uniform(0.8, 1.0)
        else:
            # Modo LECTURA: contabilidad de papel sobre precios REALES. No se
            # abre ningún socket de ejecución ni existe credencial con la que
            # firmar. Sec. 4.5: órdenes IOC, el inventario se actualiza con el
            # volumen realmente ejecutado.
            fill_c, fill_v = u_c, u_v

        inventario = cuenta.ejecutar(fill_c, fill_v, S) if cuenta is not None else 0.0
        capa.registrar_confirmacion(id_cliente, inventario)
        # Guarda 6: en LECTURA la "reconciliación" es trivialmente exacta porque
        # no hay dos contabilidades. Se invoca igualmente para que la ruta esté
        # ejercitada el día que haya un exchange al otro lado.
        capa.registrar_reconciliacion(inventario)
        seq_fill += 1
        env_arr[0]["inv_confirmado"] = inventario
        env_arr[0]["seq_fill"] = seq_fill


def correr_diagnostico_episodio(id_episodio: int, dir_salida: str) -> bool:
    """Compuerta 5 de la Sec. C.3: `diagnostico.py` ejecutado y resumen escrito.

    Es lo que convierte el bucle de episodios en un pipeline de datos y no en una
    tragamonedas. Devuelve False si no se pudo generar, y entonces la compuerta
    bloquea el reaprovisionamiento — que es exactamente lo que debe pasar: sin
    diagnóstico, recargar es tirar dinero a ciegas.
    """
    try:
        import diagnostico

        datos = diagnostico.cargar_telemetria(DIR_TELEMETRIA)
        if "id_episodio" in datos.dtype.names:
            mascara = datos["id_episodio"].astype(int) == id_episodio
            if mascara.any():
                datos = datos[mascara]
        if len(datos) < 100:
            return False
        os.makedirs(dir_salida, exist_ok=True)
        ruta = os.path.join(dir_salida, f"diagnostico_episodio_{id_episodio:03d}.txt")
        with open(ruta, "w", encoding="utf-8") as fh:
            fh.write(diagnostico.generar_reporte(datos))
        return True
    except Exception as err:
        log(f"[DIAGNOSTICO] Episodio {id_episodio}: no se pudo generar ({err}).")
        return False


async def lazo_episodios(env_arr, par_arr, modo, capa, estado, periodo: float = 2.0):
    """Máquina de episodios de la Sec. C.2, viviendo en el Motor de Red.

    Aquí y no en el orquestador porque la capa de riesgo, la cuenta y la ruta de
    cierre ya están en este proceso: mover la máquina fuera obligaría a
    sincronizar el cierre a través de memoria compartida justo en el momento en
    que menos se puede confiar en el estado.
    """
    maquina = estado["maquina"]
    while True:
        await asyncio.sleep(periodo)
        cuenta = estado.get("cuenta")
        if cuenta is None:
            continue
        S = float(env_arr[0]["P_spot"])
        if S <= 0.0:
            continue

        if maquina.estado is Estado.ARRANQUE:
            if maquina.abrir():
                estado["id_episodio"] = maquina.id_episodio
                env_arr[0]["id_episodio"] = maquina.id_episodio
                estado["halt_reportado"] = False
            continue

        if maquina.estado is Estado.CERRANDO:
            # Semántica del halt (Sec. B.5): cerrar y parar. La ruta de cierre
            # reintenta con backoff y jamás reporta éxito sin posición plana.
            ruta = riesgo.RutaDeCierre(
                cerrar_fn=lambda pos: cuenta.ejecutar(
                    max(0.0, -pos), max(0.0, pos), S
                ),
                leer_posicion_fn=lambda: cuenta.inventario,
                alertar_fn=log,
                dormir_fn=lambda s: None,
            )
            try:
                final = ruta.ejecutar()
                maquina.al_confirmar_plano(final)
                log(f"[EPISODIO] Cierre confirmado plano: {final:.6f} BTC.")
            except riesgo.ErrorDeCierre as err:
                # No se degrada a "cerrado": queda exposición abierta y hay que
                # decirlo. La máquina se detiene y exige intervención.
                maquina.detener(str(err))
            continue

        if maquina.estado is Estado.CERRADO:
            muestras = int(env_arr[0]["muestras_episodio"])
            diag_ok = correr_diagnostico_episodio(
                maquina.id_episodio, maquina.dir_salida
            )
            compuertas, _resumen = maquina.cerrar_y_evaluar(
                muestras, {"ordenes_enviadas": capa.n_enviadas}, diag_ok
            )
            log(f"[EPISODIO] Compuertas C.3: {compuertas}")
            continue

        if maquina.estado is Estado.REAPROVISIONANDO:
            maquina.sondear_equity()
            continue

        if maquina.estado is Estado.DETENIDO:
            # Terminal y ruidoso. Nunca se sale de aquí por software.
            await asyncio.sleep(30.0)
            log(f"[DETENIDO] {maquina.motivo_detencion}")


async def exchange_websocket_handler(env_arr, act_arr, par_arr, modo, filtros, cuenta, estado):
    """Watchdog con reconexión infinita y backoff exponencial (Sec. 8.4.2).

    El `except` vive DENTRO del `while True`: en una revisión anterior estaba
    fuera, de modo que el primer TimeoutError marcaba dropout y terminaba la
    corrutina, dejando al bot ciego para siempre.

    v1.3: la ingesta, la actuación y la vigilancia del reloj corren como tareas
    independientes. Que la ingesta caiga NO debe detener la evaluación de las
    guardas — al contrario, un feed muerto es justo cuando más falta hace que
    alguien esté mirando la exposición abierta.
    """
    capa = riesgo.CapaRiesgo(
        filtros,
        nocional_max_orden=float(par_arr[P_NOCIONAL_MAX_ORDEN]),
        perdida_max=float(par_arr[P_PERDIDA_MAX_EPISODIO]),
        n_ordenes_max_min=int(par_arr[P_N_ORDENES_MAX_MIN]),
        m_rechazos_max=int(par_arr[P_M_RECHAZOS_MAX]),
        t_confirm_max=float(par_arr[P_T_CONFIRM_MAX]),
        tol_inv=float(par_arr[P_TOL_INV]),
        offset_reloj_max=float(par_arr[P_OFFSET_RELOJ_MAX]),
    )
    if cuenta is not None:
        capa.registrar_equity(cuenta.equity_inicial)
    estado["capa"] = capa

    # --- Máquina de episodios (Sec. C) --------------------------------------
    def al_reiniciar(id_ep, equity_inicio):
        """Reset limpio al ABRIR episodio (Sec. C.5).

        Se hace al abrir y no al cerrar, a propósito: así el estado residual de
        un cierre fallido tampoco se hereda. Lo que este proceso puede reiniciar
        —capa de riesgo, cuenta— se reinicia aquí; ΣQ, la racha de burn-in, S_ref
        y la ventana del EMD viven en otros procesos y se reinician solos al
        detectar el cambio de `id_episodio` en memoria compartida.
        """
        capa.reiniciar(equity_inicio=equity_inicio)
        estado["muestras_episodio"] = 0
        env_arr[0]["causa_halt"] = 0
        log(f"[EPISODIO {id_ep:03d}] Reset limpio; equity inicial {equity_inicio:.2f} USD.")

    estado["maquina"] = episodios.MaquinaEpisodios(
        reaprovisionador=episodios.ReaprovisionadorManual(log),
        leer_equity_fn=lambda: (
            estado["cuenta"].equity(float(env_arr[0]["P_spot"]))
            if estado.get("cuenta") is not None
            else 0.0
        ),
        precio_fn=lambda: float(env_arr[0]["P_spot"]) or CTE.PRECIO_REFERENCIA,
        al_reiniciar=al_reiniciar,
        alertar_fn=log,
    )
    estado["al_disparar"] = estado["maquina"].al_disparar_guarda

    tareas = [
        asyncio.create_task(lazo_actuacion(env_arr, act_arr, par_arr, modo, capa, cuenta, estado)),
        asyncio.create_task(vigilancia_reloj(capa, modo)),
        asyncio.create_task(lazo_episodios(env_arr, par_arr, modo, capa, estado)),
    ]
    retry_count = 0
    try:
        while True:
            try:
                await ingesta_mercado(env_arr, modo, estado)
                retry_count = 0
            except asyncio.CancelledError:
                raise
            except Exception as err:
                # Señalización de ceguera: el Hilo Rápido anulará K y dejará
                # crecer P (Sec. 8.3.1).
                env_arr[0]["flag_dropout"] = 1
                # τ_reintento,n = τ_base·2ⁿ + U(0, j_max)   (Sec. 8.2.3)
                base = min(60.0, 0.1 * (2**retry_count))
                if retry_count % 8 == 0 or isinstance(err, mercado.EstancamientoFeed):
                    log(f"[WATCHDOG] Feed caido ({type(err).__name__}: {err}); "
                        f"reintento en {base:.1f} s.")
                    feed = estado.get("feed")
                    if feed is not None and feed.modo_degradado:
                        log(
                            "[WATCHDOG] Degradando a sondeo REST: el WebSocket de "
                            "futuros conecta pero no entrega datos desde esta red. "
                            "Precio y volumen siguen siendo REALES; lo que baja es "
                            "la resolucion temporal (4 Hz)."
                        )
                await asyncio.sleep(base + random.uniform(0.0, 0.1))
                retry_count += 1
    finally:
        for t in tareas:
            t.cancel()


def network_engine_process(shm_env_name, shm_act_name, shm_par_name, shm_mer_name):
    log("[INIT] Proc 1 (Motor de Red) — afinidad núcleo 0")
    if hasattr(os, "sched_setaffinity"):  # No existe en Windows
        os.sched_setaffinity(0, {0})

    shm_env = shared_memory.SharedMemory(name=shm_env_name)
    shm_act = shared_memory.SharedMemory(name=shm_act_name)
    shm_par = shared_memory.SharedMemory(name=shm_par_name)
    shm_mer = shared_memory.SharedMemory(name=shm_mer_name)
    env_arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
    mer_arr = np.ndarray((RING_MERCADO_SIZE,), dtype=MERCADO_DTYPE, buffer=shm_mer.buf)
    act_arr = np.ndarray((RING_BUFFER_SIZE,), dtype=ACTUATOR_DTYPE, buffer=shm_act.buf)
    par_arr = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)

    modo = mercado.modo_desde_codigo(par_arr[P_MODO])
    filtros = mercado.FiltrosInstrumento(
        symbol=SYMBOL,
        modo=modo,
        step_size=float(par_arr[P_STEP_SIZE]),
        min_qty=float(par_arr[P_MIN_QTY]),
        max_qty=0.0,
        tick_size=float(par_arr[P_TICK_SIZE]),
        min_notional=float(par_arr[P_MIN_NOTIONAL]),
        ts_lectura=time.time(),
    )
    if filtros.step_size <= 0.0 or filtros.min_notional <= 0.0:
        # Centinela de la Sec. A.2: los filtros no se leyeron de exchangeInfo.
        log(
            "[FATAL] Los filtros del instrumento siguen a cero: no se leyeron de "
            "exchangeInfo. No hay valor por defecto defendible (Sec. A.2)."
        )
        return

    # Cuenta de papel: solo en LECTURA, donde no hay cuenta que consultar. En
    # TESTNET/MAINNET el equity DEBE venir del ACCOUNT_UPDATE (Sec. B.5).
    # Se construye de forma PEREZOSA, en `lazo_actuacion`, con el primer precio
    # real: dimensionarla aquí exigiría un precio que todavía no existe.
    cuenta = None

    estado = {
        "n_ticks": 0,
        "ultimo_precio": 0.0,
        "id_episodio": int(env_arr[0]["id_episodio"]) or 1,
        "halt_reportado": False,
        "al_disparar": lambda causa, detalle: None,
        "cuenta": None,
        "modo_lectura": modo is Modo.LECTURA,
        "feed": None,
        # --- v2.0 §3: estado del productor del anillo de mercado ---
        "mer_arr": mer_arr,
        "seq_mercado": 1,  # 0 queda reservado como "slot vacio"
        # El contador de deduplicacion NO se reinicia entre reconexiones (§3.3):
        # reiniciarlo reintroduciria el duplicado justo despues de cada caida.
        "ultimo_trade_id": 0,
        "n_duplicados": 0,
        "n_huecos": 0,
        "trades_perdidos": 0,
    }

    try:
        asyncio.run(
            exchange_websocket_handler(env_arr, act_arr, par_arr, modo, filtros, cuenta, estado)
        )
    except KeyboardInterrupt:
        pass
    finally:
        shm_env.close()
        shm_act.close()
        shm_par.close()
        shm_mer.close()


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
        log(f"[TELEMETRIA] volcado fallido ({etiqueta}): {err}")


def fast_thread_process(shm_env_name, shm_mic_name, shm_act_name, shm_par_name, shm_mer_name):
    log("[INIT] Proc 2 (EAKF + NMPC + Telemetría) — afinidad núcleo 1")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {1})

    shm_env = shared_memory.SharedMemory(name=shm_env_name)
    shm_mic = shared_memory.SharedMemory(name=shm_mic_name)
    shm_par = shared_memory.SharedMemory(name=shm_par_name)
    shm_act = shared_memory.SharedMemory(name=shm_act_name)
    shm_mer = shared_memory.SharedMemory(name=shm_mer_name)

    env_arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
    mer_arr = np.ndarray((RING_MERCADO_SIZE,), dtype=MERCADO_DTYPE, buffer=shm_mer.buf)
    consumidor = ConsumidorMercado(mer_arr, "rapido")
    mic_arr = np.ndarray((1,), dtype=MICELIO_DTYPE, buffer=shm_mic.buf)
    par_arr = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)
    act_arr = np.ndarray((RING_BUFFER_SIZE,), dtype=ACTUATOR_DTYPE, buffer=shm_act.buf)

    # --- Matrices del EAKF (Sec. 7.3) ---------------------------------------
    I3 = np.eye(3)
    H = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])  # Sec. 7.3.2
    P_k = np.diag([1e-5, 1e6, 1e6])  # Sec. 7.3.4: p_S≈0, p_v,p_Rn≈1e6

    # v1.3: se ESPERA al primer tick real antes de inicializar el estado.
    # Antes se caía a un literal de 40 000 USD si el Motor de Red aún no había
    # publicado, y ese literal contaminaba TODO el arranque: γ_0 ∝ 1/S, la
    # abrazadera de posición ∝ S, el equity mínimo del episodio, la guarda de
    # resolución de control y el λS²Γ que se reporta. Con BTC a ~63 800, arrancar
    # creyendo 40 000 desplazaba el nocional máximo de posición un 36 %.
    t_espera = time.time()
    while float(env_arr[0]["P_spot"]) <= 0.0:
        if (time.time() - t_espera) > 30.0:
            log(
                "[FATAL] 30 s sin un solo precio del Motor de Red. Las guardas de "
                "arranque necesitan un precio REAL: evaluarlas contra un literal "
                "las vuelve decorativas."
            )
            return
        time.sleep(0.05)
    S0 = float(env_arr[0]["P_spot"])
    log(f"[INIT] Primer precio real recibido: S={S0:.2f} USD/BTC. Evaluando guardas.")
    x_k = np.array([[S0], [0.0], [0.0]])

    telem_buffer = np.zeros(TELEMETRY_SIZE, dtype=TELEM_DTYPE)
    telem_idx = 0
    t_ultimo_volcado = time.time()
    n_volcado = 0

    is_burnt_in = False
    burnt_ticks = 0
    last_tr_P = float(np.trace(P_k))
    t_ultimo_ciclo = time.time()

    seq_emision = 1  # 0 queda reservado como "slot vacío"
    ultimo_seq_fill = 0
    ultimo_id_episodio = -1  # Reset limpio entre episodios (Sec. C.5)
    # --- v2.0 §4: estado del filtro en reloj de transacciones ---
    # Con Δn = 1 la matriz de velocidad constante es literalmente [[1,1],[0,1]]:
    # un tick de avance, no un Δt. Se construye una vez.
    A_VEL_CONSTANTE_TICK = dinamica.matriz_A_velocidad_constante(1.0)
    A_cacheada = A_VEL_CONSTANTE_TICK
    omega_cacheada = -1.0  # fuerza la primera construccion de A_arm
    n_trades_procesados = 0
    tr_P = float(np.trace(P_k))
    deriv_tr = 0.0
    eps_nis = float('nan')
    innov = np.zeros((2, 1))
    innov_sombra = np.zeros((2, 1))
    ultima_seq_mic = -1  # frescura de R_n: actualizacion secuencial multi-tasa
    # Fila de H que observa SOLO el precio, para los ticks en que R_n es viejo.
    H_SOLO_PRECIO = np.array([[1.0, 0.0, 0.0]])
    nis_sombra = float('nan')
    hay_medicion = False
    paquete_valido = False
    muestras_episodio = 0  # Compuerta de muestras minimas de la Sec. C.3
    inventario = 0.0
    n_ciclos = 0
    en_singularidad = False

    # --- FASE 1.2: estado del criterio de burn-in por NIS --------------------
    n_ventana_nis = max(2, int(par_arr[P_NIS_VENTANA]))
    nis_lo, nis_hi = cotas_nis(n_ventana_nis, par_arr[P_NIS_CONF])
    buffer_nis = np.zeros(n_ventana_nis, dtype=np.float64)
    idx_nis = 0
    n_nis = 0  # Muestras acumuladas; la media no es válida hasta llenar la ventana
    burnt_ticks_traza = 0  # Criterio secundario, solo para contraste

    # --- v1.3 Sec. D: conmutador de la rama de A y filtro SOMBRA -------------
    conmutador = dinamica.ConmutadorRamaA(
        c_on=float(par_arr[P_C_ON]), c_off=float(par_arr[P_C_OFF])
    )
    # Sec. E.1: el sombra usa SIEMPRE la rama contraria a la que gobierna el
    # control en el arranque, de modo que la comparación existe desde el primer
    # ciclo. Comparte x0, P0, z_k, R_k y Q_k con el de producción: cualquier otra
    # diferencia contaminaría el A/B.
    sombra = dinamica.EAKFSombra(x_k, P_k, H, usa_armonico=True)
    coste_sombra_acum = 0.0
    n_coste_sombra = 0

    # --- SECCIÓN 0.5/0.6: guardas dimensionales al arranque ------------------
    # Dos asserts baratos que delatan un error dimensional en el primer ciclo en
    # vez de a los veinte minutos.
    S_guarda = float(x_k[0, 0])
    g0_guarda = CTE.gamma_0(S_guarda, par_arr[P_I_MAX], par_arr[P_DELTA_S_REF])
    CTE.verificar_cobertura_acotada(
        S_guarda, g0_guarda, par_arr[P_DELTA_S_REF], par_arr[P_I_MAX]
    )
    CTE.verificar_friccion_subcritica(S_guarda, par_arr[P_MU_OU], g0_guarda)
    log(
        CTE.resumen_constantes(
            S=S_guarda,
            r_base=par_arr[P_R_BASE],
            q_base_inv=par_arr[P_Q_INV_BASE],
            sigma_ou=par_arr[P_SIGMA_OU],
            theta_ou=par_arr[P_THETA_OU],
            mu_ou=par_arr[P_MU_OU],
        )
    )
    log(
        f"[NIS] Banda de la media movil (n={n_ventana_nis}, m={DIM_MEDICION}, "
        f"conf={par_arr[P_NIS_CONF]:.0%}): [{nis_lo:.4f}, {nis_hi:.4f}] "
        f"alrededor del valor teorico {DIM_MEDICION}."
    )

    # --- v1.3 Secs. A.3 y B.3: guardas de la capa de riesgo al arranque ------
    # Van aquí, en el Hilo Rápido, junto a las guardas dimensionales de la
    # Sec. 0.5, aunque la capa de riesgo viva en el Motor de Red: son guardas de
    # CONFIGURACIÓN, y su sitio es donde el sistema falla más ruidosamente y antes
    # de que se emita ninguna orden.
    step_mainnet = float(par_arr[P_STEP_SIZE])
    CTE.verificar_resolucion_control(
        S_guarda, step_mainnet, float(par_arr[P_NOCIONAL_MAX_ORDEN])
    )
    CTE.verificar_cap_antes_de_liquidacion(
        S_guarda,
        float(par_arr[P_APALANCAMIENTO]),
        float(par_arr[P_PERDIDA_MAX_EPISODIO]),
        float(par_arr[P_MMR]),
    )
    log(CTE.resumen_riesgo(S_guarda, step_mainnet, float(par_arr[P_MMR]), False))
    # Criterio de aceptación de la Sec. F: λS²Γ reportado al arranque. Si alguien
    # bajara I_max "para que quepa" en el presupuesto de la cuenta (Sec. B.1),
    # este número se desploma y el freno de Loeper queda desconectado en silencio.
    prod_loeper, alcanzable = CTE.loeper_alcanzable(S_guarda, float(par_arr[P_MU_OU]))
    log(
        f"[LOEPER] lambda*S^2*gamma_0 = {prod_loeper:.4f} contra un umbral de 1.0 "
        f"-> freno {'ALCANZABLE' if alcanzable else 'ESTRUCTURALMENTE INALCANZABLE'}"
    )
    if not alcanzable:
        log(
            "[!] El freno de singularidad de Loeper (Sec. 4.4.3) no puede "
            "dispararse con esta configuracion. Sospechoso principal: I_max "
            "recortado como si fuera un parametro de riesgo (Sec. B.1)."
        )

    try:
        while True:
            n_ciclos += 1
            t_ahora = time.time()
            dt_ciclo = t_ahora - t_ultimo_ciclo
            t_ultimo_ciclo = t_ahora
            if dt_ciclo < 1e-4:
                dt_ciclo = 1e-4  # Guarda contra división por cero cinemática

            # --- Reset limpio de episodio (Sec. C.5) ------------------------
            # El estado que sobrevive a un halt es fuente segura de confusión al
            # analizar el episodio siguiente. Cada proceso reinicia LO SUYO al
            # detectar el cambio de `id_episodio` en memoria compartida; no hace
            # falta señalización adicional porque el campo ya es monótono y lo
            # escribe un único productor.
            id_ep = int(env_arr[0]["id_episodio"])
            if id_ep != ultimo_id_episodio:
                ultimo_id_episodio = id_ep
                if id_ep > 0:
                    P_k = np.diag([1e-5, 1e6, 1e6])  # Sec. 7.3.4
                    x_k = np.array([[float(env_arr[0]["P_spot"])], [0.0], [0.0]])
                    is_burnt_in = False
                    burnt_ticks = 0
                    burnt_ticks_traza = 0
                    buffer_nis[:] = 0.0
                    idx_nis = 0
                    n_nis = 0
                    inventario = 0.0
                    ultimo_seq_fill = int(env_arr[0]["seq_fill"])
                    en_singularidad = False
                    conmutador.rama = dinamica.RAMA_VELOCIDAD_CONSTANTE
                    sombra = dinamica.EAKFSombra(x_k, P_k, H, usa_armonico=True)
                    last_tr_P = float(np.trace(P_k))
                    muestras_episodio = 0
                    env_arr[0]["muestras_episodio"] = 0
                    log(
                        f"[EPISODIO {id_ep:03d}] Hilo Rapido reiniciado: burn-in, "
                        f"NIS, P, inventario y rama de A a estado inicial."
                    )

            # --- Lectura lock-free del bloque estructural -------------------
            mic, _mic_ok = leer_micelio_seqlock(mic_arr)
            Omega = float(mic["Omega"])
            w_m = float(mic["w_m"])
            nu = float(mic["nu"])
            lam = float(mic["lambda_sim"])
            R_n_med = float(mic["R_n"])
            S_ref = float(mic["S_ref"])
            vol_sum_Q = float(mic["vol_sum_Q"])
            # Sec. D.2: ω_ang viene YA en rad/s desde el Hilo Lento. NO se deriva
            # aquí de w_m ni de f_hz — ese es exactamente el punto de uso donde
            # las dos trampas (el 2π y el factor 125) se cuelan sin dar síntoma.
            w_ang = float(mic["w_ang"])
            C_esp = float(mic["C_espectral"])

            dropout = int(env_arr[0]["flag_dropout"])

            # ================================================================
            # v2.0 §4 — EL FILTRO BAJO Δn = 1
            # ================================================================
            # UN TICK = UN PASO DE PREDICCION + UNA CORRECCION. Es la condicion
            # prioritaria de la v2.0 y sustituye al muestreo de `P_spot` por
            # ciclo de reloj.
            #
            # Lo que esto ELIMINA estructuralmente —no parchea: vuelve imposible
            # de expresar— es:
            #   - corregir varias veces con la misma medicion (v1.3),
            #   - descartar el resto del lote (v2.0 §1.2: era el 96 %),
            #   - meter al planificador dentro de Q (§1.3),
            #   - la inconsistencia de Δt arrastrada desde la v1.1: se disuelve,
            #     porque deja de haber Δt DENTRO del filtro.
            #
            # ⚠ x[1] CAMBIA DE UNIDADES: de USD/BTC por SEGUNDO a USD/BTC por
            # TRANSACCION. Auditoria de consumidores hecha: el unico es la
            # telemetria (`reg["x1"]`). Ninguna formula la consume, asi que el
            # cambio es seguro. La constante que la escala, `q_base*1e-2`, si
            # queda mal calibrada y pasa a ser [CALIBRAR] de la Fase 2.
            lote_trades = consumidor.leer_lote()

            # Rama de A: se decide UNA VEZ por ciclo, no por trade. `C` y `w_ang`
            # los publica el Hilo Lento a su cadencia, asi que dentro de un lote
            # no cambian y reevaluarlos por trade solo gastaria tiempo.
            rama = conmutador.actualizar(C_esp, w_ang)
            # §4.1: ω ANGULAR en rad/TICK. La conversion vive en el sitio unico
            # de `constantes_micelio` (§2.4) y NO se hace aqui: este es
            # exactamente el punto de uso donde la trampa del 2π se cuela.
            omega_ang_rad_tick = float(mic["omega_ang_rad_tick"])
            if rama == dinamica.RAMA_ARMONICO:
                # Δn = 1: `A_arm` es CONSTANTE entre cambios de regimen y se
                # cachea; solo se reconstruye cuando ω_m cambia (§4.1).
                if omega_ang_rad_tick != omega_cacheada:
                    A_cacheada = dinamica.matriz_A_armonica(omega_ang_rad_tick, 1.0)
                    omega_cacheada = omega_ang_rad_tick
                A = A_cacheada
            else:
                A = A_VEL_CONSTANTE_TICK

            # Q_k = ρ_k·diag(q_S, q_v, q_Rn),  ρ_k = 1 + γ_ω|ω_m| + γ_Q|ΣQ|  (7.3.3)
            # ρ_k va con ω_m (frecuencia de mercado), NO con Ω (estabilidad).
            # γ_ω y γ_Q DERIVADOS de los límites estructurales (Sec. 0.1), nunca
            # literales: γ_ω = 1/ω_m,max acota el primer término a ≤1, y
            # γ_Q = 1/max(ΣQ_max, C_max) acota el segundo — con max, no min, porque
            # ΣQ es volumen de mercado y tomar el mínimo maximizaría γ_Q.
            # SECCIÓN C de la v1.2: q_S pasa a ser una varianza RELATIVA al nivel
            # de precio, q_S = (σ_rel·S_k)², con las mismas unidades (USD/BTC)²
            # que ya tenía — el álgebra del filtro no cambia. Con q_S absoluto, un
            # movimiento de BTC de 45k a 90k lo dejaba mal escalado por 4×.
            q_base = par_arr[P_Q_BASE]
            rho_k = (
                1.0
                + CTE.gamma_omega(par_arr[P_OMEGA_M_MAX]) * abs(w_m)
                + CTE.gamma_Q(par_arr[P_SUMA_Q_MAX], par_arr[P_K_USD])
                * abs(vol_sum_Q)
            )

            # --- §4.2: Q POR TICK, no por paso de bucle ---------------------
            # ⚠ ESTO ERA UN DEFECTO REAL, no una mejora cosmetica. `Q` no
            # dependia de Δt, asi que Σ AⁱQAⁱᵀ crecia con el NUMERO DE PASOS y no
            # con el tiempo transcurrido: la tasa del bucle cambiaba el filtro.
            # El ruido de proceso inyectado era proporcional a cuantas veces
            # desperto el planificador, y `q_base` estaba calibrada en silencio
            # contra ~90 Hz — mover `PERIODO_HILO_RAPIDO` recalibraba el filtro
            # sin que nada lo dijera.
            #
            # `Q` por tick es varianza inyectada POR TRANSACCION, que es una
            # propiedad del mercado. `Q` por paso de bucle era varianza por
            # despertar del planificador, que es una propiedad del sistema
            # operativo. El reloj de ticks lo saca del presupuesto de ruido.
            #
            # NOTA DE INTERPRETACION: el §4.2 solo especifica `q_S_tick`. Se
            # aplica el mismo principio a las otras dos entradas —varianza por
            # tick = varianza por segundo / ν, porque los incrementos entre ticks
            # se suponen iid y la VARIANZA es lo aditivo— para que la diagonal
            # entera quede en el mismo reloj. Mezclar relojes dentro de una misma
            # matriz de covarianza seria peor que tenerla mal escalada entera.
            nu_ticks_por_s = CTE.nu_ticks_por_s_desde_anio(nu)
            sigma_rel_tick = CTE.sigma_rel_tick_desde_s(
                par_arr[P_SIGMA_REL], nu_ticks_por_s
            )
            q_S_tick = CTE.q_S_relativa(float(x_k[0, 0]), sigma_rel_tick)
            escala_tick = 1.0 / nu_ticks_por_s if nu_ticks_por_s > 0.0 else 1.0
            Q_k = rho_k * np.diag(
                [q_S_tick, q_base * 1e-2 * escala_tick, q_base * escala_tick]
            )

            # FORMA AFÍN (Sec. D.3): x_pred = x_ref + A·(x_k − x_ref), con
            # x_ref = [S_ref, 0, 0]ᵀ. Conserva x[0] = S ABSOLUTO, así que no hay
            # que auditar a los consumidores (q_S = (σ_rel·S)², γ_0(S), el centro
            # de la malla de Loeper, z₀ = P_spot, la telemetría, el NMPC).
            # Y el manejo de nodos de fase sale gratis: cuando S_ref salta, el
            # offset se aplica en la predicción siguiente sin discontinuidad en x
            # ni en P — se elimina la cascada innovación espuria -> pico de NIS ->
            # ventana del EMD a W_min -> racha de burn-in rota.
            # Con la rama de velocidad constante es algebraicamente idéntico a
            # P_pred = A·P·Aᵀ + Q y x_pred = A·x, porque [1,Δt;0,1] deja invariante
            # a [S_ref,0,0]ᵀ salvo por el propio offset.
            S_ref_pred = S_ref if S_ref > 0.0 else float(x_k[0, 0])
            x_pred, P_pred = dinamica.predecir_afin(A, x_k, P_k, Q_k, S_ref_pred)

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

            # --- v1.3: NOVEDAD DEL PAQUETE ---------------------------------
            # DIVERGE DEL PDF (Sec. 7.3): el filtro corrige UNA VEZ POR PAQUETE
            # NUEVO, no una vez por ciclo de control.
            #
            # El Hilo Rápido corre a ~90 Hz y el feed entrega a una cadencia
            # menor —2 Hz para R_n desde el Hilo Lento, y hasta 1 Hz para el
            # precio si el WebSocket degrada a sondeo REST. Corrigiendo cada
            # ciclo, la MISMA medición entra al filtro decenas de veces
            # seguidas: P se contrae como si hubiera decenas de observaciones
            # independientes, el filtro se declara mucho más seguro de lo que
            # está, y la innovación queda autocorrelacionada POR CONSTRUCCIÓN.
            #
            # Eso no es una hipótesis: es la explicación medida del ρ₁ = 0.87 de
            # `y1` que la v1.2 dejó como pregunta abierta para este documento
            # (retención de orden cero, Hilo Lento a 2 Hz contra Hilo Rápido a
            # 89 Hz). La Sec. E.3 lo resuelve excluyendo `y1` de la compuerta;
            # esto lo resuelve en el origen, y además impide que el canal `y0`
            # herede el mismo artefacto al degradar el feed — lo que habría
            # contaminado justo la serie sobre la que se decide A_arm.
            #
            # Sin paquete nuevo, z_k = ∅ y aplica la Sec. 8.3.1 tal cual: K = 0 y
            # la incertidumbre crece con Q_k. Es la semántica correcta, no un
            # apaño: entre dos observaciones el filtro efectivamente no sabe más.
            # ================= BUCLE POR TRANSACCION (§4) ====================
            # Cada trade del lote produce EXACTAMENTE un paso de prediccion y una
            # correccion. Si el lote viene vacio el filtro NO avanza, y eso es lo
            # correcto en reloj de transacciones: sin transacciones no hay
            # informacion nueva NI evolucion del estado. El tiempo del filtro es
            # el mercado, no el reloj de pared.
            for trade in lote_trades:
                # τ_d por TRANSACCION, con el reloj del exchange (Sec. 6.5). Antes
                # se leia un unico `timestamp` del bloque escalar y todo el lote
                # compartia el mismo retardo, que es falso: el primer trade de una
                # rafaga de 203 es bastante mas viejo que el ultimo.
                ts_trade = float(trade["T_trade"])
                retardo_sensor = (
                    max(0.0, t_ahora - ts_trade) if ts_trade > 0.0 else float("inf")
                )
                paquete_valido = retardo_sensor <= par_arr[P_TAU_MAX]
                hay_medicion = not (dropout == 1 or not paquete_valido)

                # --- Prediccion afin con Δn = 1 (Sec. D.3 de la v1.3) -------
                x_pred, P_pred = dinamica.predecir_afin(A, x_k, P_k, Q_k, S_ref_pred)

                # --- ACTUALIZACION SECUENCIAL MULTI-TASA -------------------
                # ⚠ EL Δn = 1 DEL §4 ARREGLA EL CANAL DE PRECIO PERO NO EL DE
                # R_n, y medirlo lo delata: con el filtro a ~94 ticks/s y el Hilo
                # Lento publicando `R_n` a 2 Hz, `y1` se repetia 11.6 veces y su
                # ρ₁ salia +0.82. Es EXACTAMENTE el bug de las 90 correcciones de
                # la v1.3, desplazado de canal — la tercera aparicion de la misma
                # familia, que es justo lo que el §3.3 advierte.
                #
                # El vector de medicion tiene DOS componentes con cadencias
                # distintas, asi que la correccion tiene que ser secuencial:
                #   R_n FRESCO  -> H completa (2 filas), m = 2
                #   R_n VIEJO   -> solo la fila del precio, m = 1
                # Asimilar un R_n retenido como si fuera nuevo contrae P por
                # informacion que no existe.
                #
                # La frescura se detecta con el contador del seqlock, que el Hilo
                # Lento incrementa en cada publicacion: no hace falta comparar
                # valores, que confundiria "no cambio" con "no llego".
                seq_mic = int(mic["seq"])
                R_n_fresco = seq_mic != ultima_seq_mic
                if R_n_fresco:
                    ultima_seq_mic = seq_mic

                if hay_medicion:
                    r_precio = par_arr[P_R_S_BASE] * math.exp(
                        par_arr[P_BETA_JITTER] * retardo_sensor
                    )
                    if R_n_fresco:
                        H_ef = H
                        R_k_ab = np.diag([r_precio, par_arr[P_R_EMD]])
                        # El precio de ESTE trade, no el ultimo del bloque.
                        z_k_ab = np.array([[float(trade["precio"])], [R_n_med]])
                    else:
                        H_ef = H_SOLO_PRECIO
                        R_k_ab = np.array([[r_precio]])
                        z_k_ab = np.array([[float(trade["precio"])]])
                else:
                    H_ef = H
                    R_k_ab = np.eye(2)
                    z_k_ab = np.zeros((2, 1))

                # ------ Sec. E.1: EAKF SOMBRA sobre el MISMO flujo ----------
                # Solo uno alimenta el control; el otro corre en sombra. Es la
                # unica comparacion honesta, porque dos corridas distintas verian
                # mercados distintos y sobre un mercado real el experimento no se
                # repite. El sombra usa la rama CONTRARIA a la vigente, asi que
                # el A/B tiene datos de ambas sin importar cual gobierne.
                sombra.usa_armonico = rama == dinamica.RAMA_VELOCIDAD_CONSTANTE
                t_sombra = time.perf_counter()
                innov_sombra, nis_sombra = sombra.paso(
                    z_k_ab, R_k_ab, Q_k, 1.0, omega_ang_rad_tick, S_ref_pred,
                    hay_medicion, H_ef,
                )
                coste_sombra_acum += time.perf_counter() - t_sombra
                n_coste_sombra += 1

                if not hay_medicion:
                    # Sec. 8.3.1: K = 0, la incertidumbre crece con Q_k (z_k = ∅).
                    x_k, P_k = x_pred, P_pred
                    innov = np.zeros((2, 1))
                    # Sin medicion no hay innovacion que evaluar: el NIS de este
                    # tick NO EXISTE (no es 0, que seria "perfectamente
                    # consistente"). Se deja fuera de la ventana movil para no
                    # sesgarla durante un dropout, que es justo cuando el filtro
                    # NO esta siendo validado.
                    eps_nis = float("nan")
                else:
                    # R_eff^(1,1) = r_S,base · e^(β·τ_d)   (Sec. 6.5 / 7.3.3)
                    innov = z_k_ab - (H_ef @ x_pred)
                    S_cov = H_ef @ P_pred @ H_ef.T + R_k_ab
                    # np.linalg.solve y no inv: mas estable si S_cov se degrada.
                    K = (P_pred @ H_ef.T) @ np.linalg.solve(
                        S_cov, np.eye(S_cov.shape[0])
                    )
                    x_k = x_pred + (K @ innov)
                    joseph = I3 - K @ H_ef
                    P_k = (joseph @ P_pred @ joseph.T) + (K @ R_k_ab @ K.T)  # Joseph

                    eps_nis = nis_escalar(innov, S_cov)
                    if math.isfinite(eps_nis):
                        buffer_nis[idx_nis] = eps_nis
                        idx_nis = (idx_nis + 1) % n_ventana_nis
                        n_nis += 1
                        # Para la ventana adaptativa del EMD (Sec. 2.2.1).
                        env_arr[0]["nis_eps"] = eps_nis

                tr_P = float(np.trace(P_k))
                deriv_tr = abs(tr_P - last_tr_P) / max(dt_ciclo, 1e-6)
                last_tr_P = tr_P

                # --- TELEMETRIA POR TICK (Sec. 8.6.1) -----------------------
                # Una fila por TRANSACCION, no por ciclo de control. Con esto la
                # tasa de observaciones del reporte pasa a ser la tasa real de
                # mercado, y `hay_medicion` deja de marcar el 99 % de relleno.
                reg = telem_buffer[telem_idx]
                reg["t_wall"] = t_ahora
                reg["tr_P"] = tr_P
                reg["x0"] = x_k[0, 0]
                reg["x1"] = x_k[1, 0]  # ⚠ ahora en USD/BTC por TRANSACCION
                reg["x2"] = x_k[2, 0]
                reg["y0"] = innov[0, 0] if innov.size else float("nan")
                # ⚠ `y1` va a NaN cuando R_n NO se asimilo en este tick. No es
                # cero: cero significaria "innovacion nula", que es lo contrario.
                # `diagnostico.ljung_box` descarta los no finitos, asi que la
                # serie de y1 queda con las ~2 observaciones/s REALES en vez de
                # con las 94 retenidas — que es lo que inflaba su rho_1 a +0.82.
                reg["y1"] = innov[1, 0] if innov.shape[0] > 1 else float("nan")
                reg["nis"] = eps_nis
                reg["dtr_dt"] = deriv_tr
                reg["id_episodio"] = int(env_arr[0]["id_episodio"])
                reg["rama_A"] = int(rama)
                # Mismo criterio que en el canal de control: cuando la
                # actualizacion fue de una sola fila, la componente de R_n no se
                # observo y va a NaN, no a cero.
                reg["y0_sombra"] = (
                    innov_sombra[0, 0] if innov_sombra.size else float("nan")
                )
                reg["y1_sombra"] = (
                    innov_sombra[1, 0] if innov_sombra.shape[0] > 1 else float("nan")
                )
                reg["nis_sombra"] = nis_sombra
                reg["hay_medicion"] = 1 if hay_medicion else 0
                reg["C_espectral"] = C_esp
                reg["w_ang"] = omega_ang_rad_tick
                telem_idx += 1
                muestras_episodio += 1
                env_arr[0]["muestras_episodio"] = muestras_episodio
                lleno = telem_idx >= TELEMETRY_SIZE
                vencido = (
                    telem_idx > 0
                    and (t_ahora - t_ultimo_volcado) >= PERIODO_VOLCADO_TELEMETRIA
                )
                if lleno or vencido:
                    # Se delega una COPIA a un hilo de I/O; el ciclo de control no
                    # se detiene. Al volcar por tiempo el bloque va PARCIAL, que
                    # es lo correcto: mejor un bloque corto en disco que un bloque
                    # largo que se pierde en el proximo corte.
                    threading.Thread(
                        target=volcar_telemetria,
                        args=(
                            telem_buffer[:telem_idx].copy(),
                            f"{os.getpid()}_{n_volcado:05d}",
                        ),
                        daemon=True,
                    ).start()
                    n_volcado += 1
                    telem_idx = 0
                    t_ultimo_volcado = t_ahora

            n_trades_procesados += len(lote_trades)

            # ================= BURN-IN: CRITERIO NIS (Fase 1.2) ==============
            # DIVERGE DEL PDF (Sec. 7.1): el PDF cierra el burn-in cuando la
            # derivada de la traza entra en una banda muerta estática ε_burn. Ese
            # criterio NO puede ser constante — CLAUDE.md ya lo dejaba abierto:
            # |ΔTr(P)|/Δt depende de Δt (en producción la latencia real, ~300 ms
            # desde Colombia, no la cadencia nominal), de ρ_k —que el Micelio
            # modula tick a tick— y del régimen de volatilidad vigente. Un umbral
            # fijo bloquea el burn-in en mercado agitado y lo deja pasar
            # prematuramente en mercado tranquilo.
            #
            # El NIS resuelve esto de raíz porque es ADIMENSIONAL y
            # AUTO-NORMALIZADO: su distribución de referencia (χ² con m g.l.) es
            # la misma en todo momento, sin importar Δt, ρ_k ni el régimen.
            #
            #   burn-in completo ⟺ media móvil de ε_k dentro de la banda χ²_m
            #                       durante N ticks continuos
            #
            # La derivada de la traza se conserva como criterio SECUNDARIO y se
            # registra en telemetría junto al NIS, para poder contrastarlos.
            media_nis = (
                float(buffer_nis.mean()) if n_nis >= n_ventana_nis else float("nan")
            )
            nis_consistente = (
                math.isfinite(media_nis) and nis_lo <= media_nis <= nis_hi
            )

            if not is_burnt_in:
                if nis_consistente:
                    burnt_ticks += 1
                    if burnt_ticks >= int(par_arr[P_N_BURN]):
                        is_burnt_in = True
                        log(
                            f"[BURN-IN] NIS consistente tras {n_ciclos} ciclos: "
                            f"media={media_nis:.4f} en [{nis_lo:.4f}, {nis_hi:.4f}] "
                            f"(teorico {DIM_MEDICION}). Tr(P)={tr_P:.4f}. "
                            f"Se libera el NMPC."
                        )
                else:
                    burnt_ticks = 0  # La ventana debe ser CONTINUA (Sec. 7.1)

                # Criterio legacy en paralelo: no gobierna nada, solo se mide.
                if deriv_tr <= par_arr[P_EPS_BURN]:
                    burnt_ticks_traza += 1
                else:
                    burnt_ticks_traza = 0

                # Latido de diagnóstico: en shadow mode el bot no opera, y sin esto
                # no hay forma de distinguir "convergiendo" de "atascado".
                if n_ciclos % 200 == 0:
                    log(
                        f"[BURN-IN] ciclo={n_ciclos} NIS_medio="
                        f"{media_nis:.4f} banda=[{nis_lo:.3f},{nis_hi:.3f}] "
                        f"racha={burnt_ticks}/{int(par_arr[P_N_BURN])} | "
                        f"(legacy dTr/dt={deriv_tr:.4f} racha={burnt_ticks_traza}) "
                        f"Tr(P)={tr_P:.4f} dt={dt_ciclo*1e3:.1f}ms "
                        f"dropout={dropout} valido={paquete_valido}"
                    )

            # v1.3 Sec. E.1: sobrecosto del filtro sombra. El criterio de
            # aceptación exige < 0.3 ms/ciclo contra los 2.2 ms medidos de
            # Loeper+NMPC; si se dispara, el A/B se está pagando con latencia de
            # control y hay que replantearlo.
            if n_ciclos % 1000 == 0 and n_coste_sombra > 0:
                coste_ms = 1e3 * coste_sombra_acum / n_coste_sombra
                if w_ang > 0.0:
                    ciclo = f"w_ang={w_ang:.5f} rad/s (periodo {2*math.pi/w_ang:.1f} s)"
                else:
                    ciclo = "sin ciclo dominante"
                log(
                    f"[A/B] rama_control={'armonico' if rama else 'vel_const'} "
                    f"C={C_esp:.3f} {ciclo} conmutaciones="
                    f"{conmutador.n_conmutaciones} | NIS_control={eps_nis:.3f} "
                    f"NIS_sombra={nis_sombra:.3f} | sobrecosto_sombra="
                    f"{coste_ms:.4f} ms/ciclo (max 0.3)"
                )

            # La TELEMETRIA ya no vive aqui: paso al bucle por transaccion (§4),
            # porque desde la v2.0 hay una fila por TICK y no una por ciclo de
            # control. Registrar por ciclo volveria a llenar la serie de ceros de
            # relleno, que es justo lo que `hay_medicion` existe para evitar.

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
                        log(
                            f"[LOEPER] Singularidad: D_min={malla['D_min']:.4f} <= eps "
                            f"(lambda={lam:.3e}/BTC, Omega={Omega:.3e}). "
                            f"Inacción total (Sec. 4.4.3)."
                        )
                    else:
                        log(
                            f"[LOEPER] Singularidad despejada: D_min={malla['D_min']:.4f}. "
                            f"Se reanuda la ejecución."
                        )

                # Cuantización comercial (Sec. 8.5) antes de publicar la orden.
                # v1.3 Sec. A.2: los cuatro filtros vienen del bloque de
                # hot-reloading, leídos de `exchangeInfo`. Antes eran
                # `apply_filters(u_c, 1e-5, 1e-5, 10.0, P_spot)` — los cuatro
                # valores inventados, y los cuatro equivocados por órdenes de
                # magnitud contra los reales (0.001 / 0.001 / 50 / 0.10).
                # Esta cuantización es solo un pre-filtro para no llenar el Ring
                # Buffer de órdenes que morirán en el floor; la cadena vinculante
                # es la de la Sec. B.4, en el Motor de Red.
                P_spot = float(env_arr[0]["P_spot"])
                u_c_q = apply_filters(
                    u_c,
                    par_arr[P_STEP_SIZE],
                    par_arr[P_MIN_QTY],
                    par_arr[P_MIN_NOTIONAL],
                    P_spot,
                )
                u_v_q = apply_filters(
                    u_v,
                    par_arr[P_STEP_SIZE],
                    par_arr[P_MIN_QTY],
                    par_arr[P_MIN_NOTIONAL],
                    P_spot,
                )

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
        shm_mer.close()


# ==============================================================================
# 6. PROCESO 3: HILO LENTO (MICELIO, OU DE LIQUIDEZ, ESTRUCTURA)
# ==============================================================================
def slow_thread_process(shm_env_name, shm_mic_name, shm_par_name, shm_mer_name):
    log("[INIT] Proc 3 (Micelio / EMD / OU) — afinidad núcleo 2")
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {2})

    shm_env = shared_memory.SharedMemory(name=shm_env_name)
    shm_mic = shared_memory.SharedMemory(name=shm_mic_name)
    shm_par = shared_memory.SharedMemory(name=shm_par_name)

    shm_mer = shared_memory.SharedMemory(name=shm_mer_name)
    env_arr = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
    mic_arr = np.ndarray((1,), dtype=MICELIO_DTYPE, buffer=shm_mic.buf)
    par_arr = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)
    mer_arr = np.ndarray((RING_MERCADO_SIZE,), dtype=MERCADO_DTYPE, buffer=shm_mer.buf)
    consumidor_lento = ConsumidorMercado(mer_arr, "lento")

    modo = mercado.modo_desde_codigo(par_arr[P_MODO])
    lam_ruido = float(par_arr[P_MU_OU])  # Arranca en la media: evita transitorio
    t_ultimo = time.time()

    # Estado del ciclo estructural vigente (Sec. 2.6 / 6.3)
    S_ref = 0.0
    Q_ref = 0.0  # Q_transado en el último nodo de fase
    t_ultimo_nodo = 0.0  # Para el período refractario del detector de nodos
    n_nodos = 0
    signo_imf_prev = 0  # Signo de la IMF dominante en t0, iteración anterior

    # Ventana deslizante de precios que alimenta la EMD (Sec. 2.2).
    buffer_precios = []
    t_ultima_muestra = 0.0
    ticks_desde_muestra = 0  # §5.1: muestreo del EMD cada K transacciones
    # --- v2.0 §6: reloj de volumen y factor phi' ---
    # `phi_prima_acum` es el acumulador MONOTONO del reloj vigente. Distinto de
    # SigmaQ (`sum_Q`), que se reinicia en cada nodo de fase: un reloj que
    # retrocede a cero varias veces por hora no es un reloj.
    reloj_vigente = CTE.reloj_desde_codigo(par_arr[P_RELOJ])
    phi_prima_acum = 0.0
    phi_prima_ultimo_nodo = 0.0
    Q_acum_prev = float(env_arr[0]["Q_acumulado_total"])
    n_ticks_prev_phi = int(env_arr[0]["n_ticks"])
    dt_muestreo = 0.0  # Intervalo real de muestreo, suavizado por EMA
    w_m = 0.0  # [1/Ticks] hasta que el HHT produzca la primera estimación
    R_n = 0.0  # [USD/BTC] residuo macro de la EMD
    f_hz_actual = 0.0  # Frecuencia dominante en Hz, para el refractario
    # --- v1.3 Sec. D: frecuencia angular y concentración espectral ---
    f_hz_ema = 0.0  # f suavizada por EMA (Sec. D.4.1), antes de derivar ω_ang
    w_ang = 0.0  # [rad/s] = 2π·f_hz_ema. Para el reloj de pared (Sec. D.2)
    omega_ang_rad_tick = 0.0  # [rad/Tick] = 2π·ω_m. Para A_arm bajo Δn=1 (§4.1)
    C_espectral = 0.0  # Fracción de energía de la IMF dominante, [0,1]

    # Historia para las derivadas respecto a T̄ de la Sec. 1.4
    T_prev = 0.0
    sumQ_prev = 0.0
    w_m_prev = 0.0
    dS_prev = 0.0
    Omega = 0.0
    ultimo_id_episodio = -1  # Reset limpio entre episodios (Sec. C.5)

    try:
        while True:
            t_now = time.time()
            # dt_sim medido contra el reloj de pared: con dt_sim = 0.1 fijo y un
            # sleep de 0.5 s el OU revertía a la media 5× más lento de lo calibrado.
            dt_sim = t_now - t_ultimo
            t_ultimo = t_now
            if dt_sim <= 0.0:
                dt_sim = PERIODO_HILO_LENTO

            # --- Reset limpio de episodio (Sec. C.5) ------------------------
            # ΣQ, S_ref y la ventana del EMD son exactamente el estado que el
            # documento nombra: arrastrarlo de un episodio al siguiente haría que
            # el nodo de fase y el volumen acumulado vinieran de un mercado que
            # ya no es el que se está midiendo.
            id_ep = int(env_arr[0]["id_episodio"])
            if id_ep != ultimo_id_episodio:
                ultimo_id_episodio = id_ep
                if id_ep > 0:
                    buffer_precios.clear()
                    t_ultima_muestra = 0.0
                    ticks_desde_muestra = 0
                    phi_prima_acum = 0.0
                    phi_prima_ultimo_nodo = 0.0
                    Q_acum_prev = float(env_arr[0]["Q_acumulado_total"])
                    n_ticks_prev_phi = int(env_arr[0]["n_ticks"])
                    S_ref = 0.0
                    Q_ref = float(env_arr[0]["Q_transado"])
                    t_ultimo_nodo = 0.0
                    n_nodos = 0
                    signo_imf_prev = 0
                    w_m = 0.0
                    f_hz_ema = 0.0
                    w_ang = 0.0
                    omega_ang_rad_tick = 0.0
                    C_espectral = 0.0
                    R_n = 0.0
                    log(
                        f"[EPISODIO {id_ep:03d}] Hilo Lento reiniciado: SigmaQ, "
                        f"S_ref, ventana del EMD y w_m a estado inicial."
                    )

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

            # --- Ventana adaptativa del EMD (Sec. 2.2.1) --------------------
            # FASE 1.1: segundo consumo del NIS. El escalar ε_k lo publica el Hilo
            # Rápido; aquí modula el tamaño de la ventana temporal del EMD, que
            # hasta ahora era fija. Estable (ε_k→0) => W_k→W_max; salto abrupto
            # (ε_k≫0) => W_k→W_min, descartando la memoria pasada.
            eps_nis = float(env_arr[0]["nis_eps"])
            W_k = ventana_adaptativa(
                eps_nis, par_arr[P_GAMMA_NIS], par_arr[P_W_MIN], par_arr[P_W_MAX]
            )

            # --- CADENA EMD -> HILBERT (Secs. 2.1 a 2.6) --------------------
            # PRIORIDAD 1 de la v1.2: sustituye los mocks de ω_m y R_n. Mientras
            # R_n fuera un coseno puro, el NIS y el Ljung-Box estaban midiendo el
            # generador de mocks y no el filtro, porque R_n entra directamente en
            # el vector de medición z_k (Sec. 7.3.1).
            # ================================================================
            # v2.0 §5 — EL EMD DEJA DE TAMIZAR UNA ESCALERA
            # ================================================================
            # ⚠ EL DEFECTO QUE ESTO CORRIGE. Hasta la v1.3 el buffer se llenaba
            # con RELOJ DE PARED: se tomaba `P_spot` cada 0.5 s "estuviera o no
            # cambiado". Con el feed a 0.76 Hz, cerca de la mitad de las muestras
            # eran DUPLICADOS LITERALES.
            #
            # El comentario que habia aqui se preocupaba de que el espaciado
            # IRREGULAR rompiera los splines de la EMD — y resolvio la
            # uniformidad de la malla TEMPORAL justo mientras la malla de VALORES
            # se volvia una escalera: uniforme en t, constante a trozos en S. Y la
            # transformada de Hilbert de un escalon tiene contenido en TODO el
            # espectro, que es la singularidad que el propio docstring de `hht.py`
            # advierte. O sea que ω_m y C podian estar contaminados EN ORIGEN.
            #
            # La solucion satisface a la vez las dos exigencias que estaban en
            # conflicto: se muestrea CADA K TRANSACCIONES, lo que da malla
            # uniforme (en ticks) Y ausencia de duplicados. En reloj de
            # transacciones el jitter del planificador sencillamente no existe,
            # asi que la preocupacion original desaparece en vez de resolverse.
            #
            # K conserva la cobertura de ventana actual: con W = 384 muestras y
            # periodos modales de ~795 ticks, K ≈ ν · PERIODO_MUESTREO_EMD es el
            # punto de partida. [CALIBRAR] contra la distribucion real de ν.
            nu_ticks_por_s_lento = CTE.nu_ticks_por_s_desde_anio(nu)

            # --- §6.3: avance del reloj phi' en este ciclo -------------------
            # dphi' = dQ es el INCREMENTO, nunca el acumulado SigmaQ. Se lee del
            # acumulador monotono `Q_acumulado_total`, que el Motor de Red
            # incrementa por transaccion y que NO se reinicia en los nodos.
            Q_acum_ahora = float(env_arr[0]["Q_acumulado_total"])
            dq_usd = max(0.0, Q_acum_ahora - Q_acum_prev)
            Q_acum_prev = Q_acum_ahora
            dn_ticks = max(0, n_ticks - n_ticks_prev_phi)
            n_ticks_prev_phi = n_ticks
            phi_prima_acum += CTE.avance_de_reloj(dn_ticks, dq_usd, reloj_vigente)

            k_muestreo = max(1, int(round(nu_ticks_por_s_lento * PERIODO_MUESTREO_EMD)))
            lote_lento = consumidor_lento.leer_lote()
            for trade in lote_lento:
                ticks_desde_muestra += 1
                if ticks_desde_muestra >= k_muestreo:
                    ticks_desde_muestra = 0
                    buffer_precios.append(float(trade["precio"]))
                    if len(buffer_precios) > int(par_arr[P_W_MAX]):
                        buffer_precios.pop(0)
            # `dt_muestreo` pasa a ser el espaciado de la malla EN SEGUNDOS que
            # corresponde a K ticks, que es lo que `hht.analizar_ventana` necesita
            # para devolver f en Hz. La malla es uniforme en TICKS; esta es su
            # traduccion al reloj de pared, y por eso se calcula aqui y no se
            # mide con un cronometro.
            if nu_ticks_por_s_lento > 0.0:
                dt_muestreo_nuevo = k_muestreo / nu_ticks_por_s_lento
                dt_muestreo = (
                    0.9 * dt_muestreo + 0.1 * dt_muestreo_nuevo
                    if dt_muestreo > 0
                    else dt_muestreo_nuevo
                )

            hht_ok = False
            if len(buffer_precios) >= max(32, min(W_k, int(par_arr[P_W_MAX]))):
                ventana = np.asarray(buffer_precios[-W_k:], dtype=np.float64)
                res_hht = hht.analizar_ventana(ventana, dt_muestreo)
                if res_hht["valido"]:
                    hht_ok = True
                    # Sec. 0.4: el HHT entrega f en Hz (tiempo físico) y el modelo
                    # exige ω_m en 1/Ticks. La conversión pasa por ν. Sin esto,
                    # γ_ω·ω_m no es adimensional en ejecución aunque lo sea en el
                    # papel. Ojo: ν se almacena en Ticks/AÑOS.
                    if nu > 0.0:
                        w_m_bruto = CTE.omega_m_desde_hz(res_hht["f_hz"], nu)
                        # NOTA DE INTERPRETACION: la Sec. 2.5 define ω_m como una
                        # cantidad instantánea, pero el estimador conserva un p90
                        # de ~20 % de error incluso con el promediado terminal del
                        # colapso espectral. Como ω_m alimenta ρ_k (Sec. 7.3.3) y
                        # c²_vol (Sec. 4.5), ese ruido se propagaría al filtro y a
                        # la malla de Loeper. Se suaviza con una EMA ligera: el
                        # ciclo estructural cambia en decenas de segundos, así que
                        # la constante de tiempo del filtro es muy inferior a la
                        # escala de la señal que sigue.
                        w_m = 0.7 * w_m + 0.3 * w_m_bruto if w_m > 0.0 else w_m_bruto
                    R_n = res_hht["R_n"]  # Residuo REAL de la EMD, no un coseno
                    f_hz_actual = res_hht["f_hz"]

                    # --- v1.3 Sec. D.4.1: EMA sobre f_hz ANTES de derivar ω_ang.
                    # `A_arm` asciende la frecuencia de modulador de Q a
                    # DETERMINANTE DE LA DINÁMICA DEL ESTADO: hoy un ω_m ruidoso
                    # solo ensancha la incertidumbre, con el armónico haría
                    # ruidosa la transición misma. La mediana sobre las 12
                    # muestras terminales ya vive en `hht.py`; esta EMA va encima.
                    f_hz_ema = (
                        0.7 * f_hz_ema + 0.3 * f_hz_actual
                        if f_hz_ema > 0.0
                        else f_hz_actual
                    )
                    # Sec. D.1 y D.2: ω_ang [rad/s] = 2π·f_hz, publicada APARTE de
                    # ω_m [1/Ticks]. Dos variables distintas para dos usos
                    # distintos, sin conversión en el punto de uso — que es donde
                    # el 2π y el factor 125 se cuelan sin dar síntoma.
                    w_ang = dinamica.omega_angular_desde_hz(f_hz_ema)
                    # §4.1: la variante en rad/TICK, para A_arm bajo Δn = 1. La
                    # conversion pasa por `omega_m_desde_hz` (ciclos/tick) y
                    # luego por el 2π, cada una en su funcion unica del §2.4.
                    omega_ang_rad_tick = CTE.omega_ang_rad_tick_desde_ciclos(
                        CTE.omega_m_desde_hz(f_hz_ema, nu) if nu > 0.0 else 0.0
                    )
                    C_espectral = float(res_hht.get("C", 0.0))

            # --- Nodo de fase y ΔS (Sec. 2.6) -------------------------------
            if S_ref <= 0.0:
                S_ref = S_actual
                Q_ref = Q_bruto

            nodo = False
            if hht_ok:
                # Nodo = cruce por cero de la IMF dominante (Sec. 2.6). Se detecta
                # por CAMBIO DE SIGNO EN t0 ENTRE ITERACIONES, no comparando las
                # dos últimas muestras de una misma ventana: el borde de la EMD se
                # recalcula en cada llamada y esa comparación intra-ventana
                # resultaba demasiado inestable para disparar de forma fiable.
                imf_t0 = res_hht["imf_dom_t0"]
                signo_imf = 1 if imf_t0 > 0 else (-1 if imf_t0 < 0 else 0)
                if signo_imf != 0 and signo_imf_prev != 0 and signo_imf != signo_imf_prev:
                    # ------------------------------------------------------------
                    # v2.0 §6.4 — REFRACTARIO EN φ', NO EN SEGUNDOS
                    # ------------------------------------------------------------
                    # El principio: "para ponderar no uses el transcurso del
                    # tiempo, usa el mercado que pasó". Este es el caso más claro
                    # de todos: no quieres dos nodos de fase separados por poco
                    # TIEMPO, quieres que estén separados por poco MERCADO.
                    #
                    # Con un refractario temporal, en un mercado muerto disparas
                    # nodos espurios —pasa el tiempo sin que pase nada— y en una
                    # ráfaga los suprimes justo cuando son reales. Medir la
                    # separación en transacciones (o en cuantos de volumen)
                    # elimina las dos fallas a la vez.
                    #
                    # El umbral se deriva del mismo criterio de antes —0.8 medios
                    # ciclos— pero expresado en el reloj vigente: medio ciclo son
                    # (0.5/f_hz) segundos, que a ν ticks/s son (0.5/f_hz)·ν ticks.
                    if f_hz_actual > 0.0 and nu_ticks_por_s_lento > 0.0:
                        medio_ciclo_s = 0.5 / f_hz_actual
                        umbral_refractario = 0.8 * CTE.avance_de_reloj(
                            dn_ticks=int(medio_ciclo_s * nu_ticks_por_s_lento),
                            dq_usd=medio_ciclo_s * nu_ticks_por_s_lento * CTE.DELTA_Q_ESTRELLA,
                            reloj=reloj_vigente,
                        )
                    else:
                        # Sin estimación de frecuencia todavía: se cae al criterio
                        # temporal de la v1.3, convertido al reloj vigente.
                        umbral_refractario = CTE.avance_de_reloj(
                            dn_ticks=int(5.0 * PERIODO_HILO_LENTO * max(nu_ticks_por_s_lento, 1.0)),
                            dq_usd=5.0 * PERIODO_HILO_LENTO * max(nu_ticks_por_s_lento, 1.0) * CTE.DELTA_Q_ESTRELLA,
                            reloj=reloj_vigente,
                        )
                    if (phi_prima_acum - phi_prima_ultimo_nodo) >= umbral_refractario:
                        nodo = True
                if signo_imf != 0:
                    signo_imf_prev = signo_imf

            if nodo:
                phi_prima_ultimo_nodo = phi_prima_acum
                S_ref = S_actual  # Sec. 2.6: S_ref = S(T̄_nodo)
                Q_ref = Q_bruto  # ΣQ se reinicia a cero (Sec. 6.3)
                t_ultimo_nodo = t_now
                n_nodos += 1

            delta_S = S_actual - S_ref

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
            # v1.3 Sec. B.6: en TESTNET λ sigue viniendo del proceso OU, y es
            # correcto — no hay impacto de mercado sobre un libro simulado.
            #
            # NOTA DE INTERPRETACION — CONSECUENCIA QUE HAY QUE DEJAR ESCRITA:
            # **Testnet NO puede calibrar μ_OU, θ_OU, σ_OU, η ni λ_min.** Esas
            # cinco salen del feed de solo lectura de Mainnet. Si una corrida de
            # Testnet pareciera "confirmarlas", estaría confirmando el propio
            # proceso OU que las genera, que es una tautología.
            # Con I_max intacto (Sec. B.1), λS²Γ se mantiene en rango realista y
            # el freno de singularidad SÍ se ejercita: una razón más para no
            # tocar I_max.
            if modo is Modo.TESTNET:
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
            # R_n ya no es un mock: sale del residuo real de la EMD unas líneas
            # más arriba (`res_hht["R_n"]`). Hasta que el buffer se llene lo
            # suficiente para el primer tamizado, se usa el precio actual como
            # mejor estimación de la tendencia — es lo que el residuo converge a
            # ser cuando no hay historia que descomponer.
            if not hht_ok and R_n <= 0.0:
                R_n = S_actual

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
                    "W_k": float(W_k),
                    "w_ang": w_ang,
                    "omega_ang_rad_tick": omega_ang_rad_tick,
                    "C_espectral": C_espectral,
                },
            )

            time.sleep(PERIODO_HILO_LENTO)

    except KeyboardInterrupt:
        pass
    finally:
        shm_env.close()
        shm_mic.close()
        shm_par.close()
        shm_mer.close()


# ==============================================================================
# 7. ORQUESTADOR (Sec. 7.7)
# ==============================================================================
NOMBRES_SHM = ("shm_ent", "shm_mic", "shm_act", "shm_par", "shm_mer")


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
        log("\n>> Interrupción externa (SIG). Limpieza de /dev/shm...")
        limpiar_recursos()
        sys.exit(0)

    signal.signal(signal.SIGINT, unhook_graceful)
    signal.signal(signal.SIGTERM, unhook_graceful)

    try:
        log("== MICELIO MULTICORE DEPLOY INITIALIZING ==")
        # ANTES de reservar nada: si hay otro Micelio vivo, esto aborta. Va aquí y
        # no dentro de `allocate_shared_memory` porque esa función tiene que poder
        # reciclar bloques huérfanos —su razón de existir en Windows— y la
        # distinción entre "huérfano" y "de alguien que sigue trabajando" es
        # justo lo que el latido resuelve.
        verificar_instancia_unica("shm_ent")
        shm_env = allocate_shared_memory("shm_ent", ENV_DTYPE.itemsize)
        shm_mic = allocate_shared_memory("shm_mic", MICELIO_DTYPE.itemsize)
        # El actuador reserva los 16 slots del Ring Buffer, no uno solo.
        shm_act = allocate_shared_memory(
            "shm_act", ACTUATOR_DTYPE.itemsize * RING_BUFFER_SIZE
        )
        # v2.0 §3.1: anillo de datos de mercado. Es el bloque que sustituye a la
        # casilla escalar P_spot y deja de tirar el 96 % de las transacciones.
        shm_mer = allocate_shared_memory(
            "shm_mer", MERCADO_DTYPE.itemsize * RING_MERCADO_SIZE
        )
        shm_par = allocate_shared_memory(
            "shm_par", TOTAL_PARAMS * np.dtype(np.float64).itemsize
        )
        mem_refs.extend([shm_env, shm_mic, shm_act, shm_par, shm_mer])

        p_array = np.ndarray((TOTAL_PARAMS,), dtype=np.float64, buffer=shm_par.buf)
        p_array[:] = initialize_default_parameters()

        # --- v1.3 Sec. A.2: filtros REALES antes de arrancar nada ------------
        # Se leen aquí, en el padre, y no en cada hijo: bajo `spawn` cada proceso
        # pagaría su propia latencia de red y —peor— podrían quedarse con lecturas
        # distintas si Binance cambiara un filtro entre medias. Una sola lectura,
        # publicada en el bloque compartido, es lo que garantiza que los tres
        # procesos operen contra el mismo instrumento.
        log(f"== MODO DE OPERACION: {MODO.value} ==")
        if MODO is Modo.LECTURA:
            log(
                "   Feed publico REAL de Mainnet. Sin credenciales cargadas: la "
                "ejecucion es fisicamente imposible (precondicion de la v1.3)."
            )
        filtros, filtros_mainnet = cargar_filtros_en_parametros(p_array, MODO)
        # ⚠ La guarda de resolucion se evalua contra el stepSize de MAINNET sea
        # cual sea el modo: Testnet es 10x mas fino y dejaria pasar una
        # configuracion que en produccion degrada el NMPC a un interruptor.
        p_array[P_STEP_SIZE] = max(filtros.step_size, filtros_mainnet.step_size)
        if filtros.step_size != filtros_mainnet.step_size:
            log(
                f"[EXCHANGE] stepSize de {MODO.value} ({filtros.step_size:g}) difiere "
                f"del de MAINNET ({filtros_mainnet.step_size:g}). Se publica el mas "
                f"GRUESO para que la guarda de resolucion (Sec. A.3) y la "
                f"cuantizacion sean conservadoras contra produccion."
            )

        p1 = mp.Process(
            target=network_engine_process,
            args=("shm_ent", "shm_act", "shm_par", "shm_mer"),
        )
        p2 = mp.Process(
            target=fast_thread_process,
            args=("shm_ent", "shm_mic", "shm_act", "shm_par", "shm_mer"),
        )
        p3 = mp.Process(
            target=slow_thread_process,
            args=("shm_ent", "shm_mic", "shm_par", "shm_mer"),
        )
        procesos.extend([p1, p2, p3])
        for pr in procesos:
            pr.start()

        nombres = {id(p1): "Motor de Red", id(p2): "Hilo Rapido", id(p3): "Hilo Lento"}
        env_latido = np.ndarray((1,), dtype=ENV_DTYPE, buffer=shm_env.buf)
        env_latido[0]["pid_orquestador"] = os.getpid()
        while True:
            # Latido: marca que esta instancia sigue viva, para que un segundo
            # arranque se niegue en vez de compartir la memoria en silencio.
            env_latido[0]["latido"] = time.time()
            time.sleep(PERIODO_LATIDO)
            muertos = [pr for pr in procesos if not pr.is_alive()]
            if muertos:
                for pr in muertos:
                    log(
                        f"[SUPERVISOR] Murio el proceso '{nombres.get(id(pr), pr.name)}' "
                        f"(pid={pr.pid}, exitcode={pr.exitcode}). Apagando el sistema."
                    )
                break

    except KeyboardInterrupt:
        pass
    except ErrorSegundaInstancia as err:
        # No se limpia la memoria compartida: es de la OTRA instancia, que sigue
        # trabajando. Liberarla aquí sería exactamente el daño que se evita.
        log(f"[ABORTADO] {err}")
        cerrado.set()
        return
    except Exception as err:
        log(f"[FATAL ORQUESTADOR]: {err}")
    finally:
        # try...finally garantiza el unlink incluso ante excepciones no-señal;
        # antes se dependía solo del signal handler y cualquier fallo filtraba
        # memoria compartida en /dev/shm.
        limpiar_recursos()


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)  # Obligatorio por el contexto CUDA (7.7.1)
    main()
