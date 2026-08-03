"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: constantes_micelio.py — Constantes de acoplamiento del modelo

Implementa la Sección 0 de ORDEN_TRABAJO_CALIBRACION_1.1.md.

REGLA DURA DEL ORDEN DE TRABAJO
-------------------------------
Las constantes de interconexión (γ_0, γ_ω, γ_Q, κ, μ) **no son valores libres**:
cada una se DERIVA de un límite físico o de capital del sistema. Se implementan
como FÓRMULAS evaluadas a partir de esos límites, nunca como literales numéricos.
Un número mágico en el código vuelve a abrir el agujero que esto cierra.

Nadie más define estas constantes. El resto del código las importa de aquí.

ORGANIZACIÓN
------------
  1. Límites estructurales (Sec. 0.2) — los únicos valores libres que quedan.
  2. Constantes derivadas (Sec. 0.1) — funciones puras de los límites.
  3. Conversión de unidades de ω_m (Sec. 0.4).
  4. Guardas de ejecución (Sec. 0.5).
  5. Adimensionales por declaración (Sec. 0.7).

Las funciones derivadas son PURAS: reciben los límites como argumentos explícitos,
de modo que pueden evaluarse tanto con los valores por defecto de este módulo como
con los que viven en el bloque de memoria compartida de hot-reloading. Eso es lo
que garantiza que κ, μ y Ω_crit no puedan desincronizarse (Fase 3 del orden de
trabajo: "los tres deben recalcularse juntos").
"""

import math

SEGUNDOS_POR_ANIO = 31_557_600.0  # Año juliano


# ==============================================================================
# 1. LÍMITES ESTRUCTURALES (Sec. 0.2)
# ==============================================================================
# Estos SÍ son valores libres, y son los únicos. Todo lo demás se deriva.

# I_max — exposición máxima absoluta de inventario.
#   Unidades: BTC
#   Es también el techo de las restricciones de caja de la Sec. 6.2 del PDF: la
#   cantidad máxima de criptoactivo que el sistema puede llegar a tener abierta.
I_MAX = 0.50  # [BTC]

# Precio de referencia para expresar límites de capital en USD.
#   Unidades: USD/BTC
#   No entra en ninguna fórmula del modelo; solo convierte I_MAX a K_USD.
PRECIO_REFERENCIA = 45_000.0  # [USD/BTC]

# K_USD — capital operativo total.
#   Unidades: USD
#   RENOMBRADO POR LA v1.2 (Sec. A.1). La v1.1 lo llamaba `C_max`, que colisionaba
#   con el `C_max` de la Sec. 6.2 del PDF — allí es la capacidad de billetera que
#   acota ū_compra/ū_venta, y esas están en BTC. El PDF llegó primero y es la
#   fuente, así que `C_max` queda RESERVADO para su significado del PDF y el
#   capital en USD pasa a llamarse K_USD.
#   Las dos magnitudes siguen siendo la misma cantidad en distintas unidades, así
#   que se deriva una de otra: el cap en BTC de las cajas es I_MAX, y K_USD es su
#   equivalente en USD.
K_USD = I_MAX * PRECIO_REFERENCIA  # [USD]

# ω_m,max — mayor frecuencia modal observada vía HHT.
#   Unidades: 1/Ticks
#   MEDIDO con la cadena `hht.py` ya integrada (25 semillas sobre el ciclo
#   estructural de 40 s del generador, con ν ≈ 20 ticks/s):
#       f observada: mediana 0.02511 Hz, p95 0.02695, max 0.02775  (real 0.025)
#       ω_m = f/ν  : mediana 1.26e-3, p95 1.35e-3, max 1.39e-3  [1/Ticks]
#   Se conserva 3.0e-3 —algo más del doble del máximo observado— a propósito:
#   ese máximo corresponde a UN solo régimen (el ciclo de 40 s del mock), y
#   γ_ω = 1/ω_m,max debe acotar el término también en regímenes más rápidos.
#   Con este valor, γ_ω·ω_m ≈ 0.42 en régimen típico, dejando margen sin saturar.
#   TODO(calibración demo): confirmar contra el rango real de regímenes de
#   Mainnet; si aparecen ciclos más rápidos que ~20 s, subir esta cota.
#   Piso teórico (Sec. 0.3): ω_m es la inversa de los ticks para completar un
#   ciclo, y el mínimo son 2 ticks (Nyquist en espacio de ticks), luego
#   ω_m,max ≤ 0.5 /Ticks  →  γ_ω ≥ 2 Ticks. Ver `verificar_nyquist_omega_m`.
OMEGA_M_MAX = 3.0e-3  # [1/Ticks]

# ΔS_ref — desplazamiento absoluto esperado respecto al nodo de fase (Sec. 2.6).
#   Unidades: USD/BTC
#   RENOMBRADO POR LA v1.2 (Sec. A.1), y este era el renombre peligroso. La v1.1
#   lo llamaba `ΔS_max`, mismo símbolo que la Sec. 7.4.1 del PDF usa para el
#   MARGEN PORCENTUAL que fija el ancho del dominio de la malla de Loeper. Son
#   cantidades de naturaleza distinta —una absoluta en USD/BTC, otra una fracción
#   adimensional— y confundirlas escala mal γ_0 por varios órdenes de magnitud.
#   `ΔS_max` queda RESERVADO para el margen de malla del PDF, que vive en
#   `MARGEN_MALLA_REL`.
#   Coherencia exigida: el dominio de la malla debe cubrir el desplazamiento
#   esperado, es decir  MARGEN_MALLA_REL · S ≥ ΔS_ref. La guarda
#   `verificar_dominio_malla` se mantiene precisamente para vigilar esa relación.
#   TODO(calibración demo): debe salir de la distribución real de |S − S_ref|.
DELTA_S_REF = 150.0  # [USD/BTC]

# Margen relativo del dominio espacial de la malla de Loeper (Sec. 7.4.1).
#   Unidades: adimensional (fracción de S_k)
MARGEN_MALLA_REL = 0.02  # ±2 %

# Ω_crit — umbral de crisis estructural.
#   Unidades: BTC/Ticks²
#   TODO(Fase 3): se calibra por búsqueda de orden cero sobre datos de Testnet.
#   NO es un valor cualquiera: κ y μ se derivan de él, así que moverlo mueve
#   simultáneamente la pendiente de R(Ω) y la de q_inv(Ω). Por eso está aquí y
#   no suelto en el código: los tres se recalculan juntos por construcción.
OMEGA_CRIT = 2.0  # [BTC/Ticks²]

# max(ΣQ_T̄) — máximo histórico del volumen neto acumulado por ciclo estructural.
#   Unidades: USD
#   TODO(calibración demo): máximo empírico de la integral por tramos de la
#   Sec. 6.3, reiniciada en cada nodo de fase.
SUMA_Q_MAX = 1.3e5  # [USD]

# c, c' — factores de escala adimensionales de κ y μ (Sec. 0.2).
#   En Ω = Ω_crit el costo de ejecución vale (1+c)× el fee base: con c = 100 son
#   101×, de modo que el solver difiere órdenes ANTES de que las restricciones de
#   caja se activen. Es la "aversión preventiva" de la Sec. 6.2 del PDF.
C_ESCALA_KAPPA = 100.0  # c   [adimensional]
C_ESCALA_MU = 100.0  # c'  [adimensional]

# σ_rel — volatilidad relativa de precio por paso del Hilo Rápido.
#   Unidades: adimensional (fracción de S)
#   Sección C de la v1.2. Parametriza la varianza de proceso del precio como
#   RELATIVA en vez de absoluta:  q_S = (σ_rel · S_k)².
#   Las unidades resultantes son (USD/BTC)², idénticas a las que q_S ya debía
#   tener, así que el álgebra del filtro no cambia en absoluto. Lo que cambia es
#   que q_S deja de depender del NIVEL de precio: con q_S fijo, un movimiento de
#   BTC de 45k a 90k lo deja mal escalado por un factor de 4; con σ_rel fijo se
#   corrige solo.
#   [CALIBRAR] debe anclarse a la volatilidad realizada del par por intervalo de
#   control. El valor de abajo reproduce el q_S = 0.01 anterior a S = 45 000
#   (sqrt(0.01)/45000 = 2.22e-6), es decir, no cambia el comportamiento actual:
#   solo lo vuelve invariante al nivel de precio.
SIGMA_REL = 2.22e-6  # [adimensional]


# ==============================================================================
# 2. CONSTANTES DERIVADAS (Sec. 0.1) — funciones puras de los límites
# ==============================================================================
def gamma_0(S: float, i_max: float = I_MAX, delta_s_ref: float = DELTA_S_REF) -> float:
    """γ_0 — curvatura del payoff terminal de cobertura.

        γ_0 = I_max / (S · ΔS_max)                       [BTC³/USD²]

    Sección del PDF: 7.4 (condición terminal U = ½γ_0(S − S_ref)²).

    Derivación (Sec. 0.3 del orden de trabajo): se exige que la cobertura
    denominada en BTC, S·∂U/∂S, sature I_max cuando el precio se aleja ΔS_max del
    nodo de fase:

        S · ∂U/∂S |_{S−S_ref = ΔS_max} = S · γ_0 · ΔS_max = I_max
        ⟹ γ_0 = I_max / (S · ΔS_max)

    Comprobación dimensional: [BTC] / ([USD/BTC]·[USD/BTC]) = BTC³/USD² ✓, que es
    exactamente [Γ] en la Sec. 4.4.1 — porque Γ = ∂²U/∂S² = γ_0 exactamente.

    ⚠ γ_0 ES UN GAMMA, NO UN DELTA. Una versión anterior de esta constante tenía
    unidades de Delta (BTC²/USD); daba λS²Γ ≈ 25, o sea singularidad de Loeper
    permanente desde el primer ciclo, y una cobertura implícita de 169 BTC contra
    un límite de 0.5 BTC.

    ⚠ DEPENDE DE S. No es un parámetro estático: hay que reevaluarla en cada ciclo
    de control con el precio filtrado S_k. Guardarla como literal reintroduce
    exactamente el error que esta fórmula elimina.
    """
    if S <= 0.0 or delta_s_ref <= 0.0:
        raise ValueError(f"gamma_0: S y ΔS_max deben ser positivos (S={S}, ΔS_max={delta_s_ref})")
    return i_max / (S * delta_s_ref)


def q_S_relativa(S: float, sigma_rel: float = SIGMA_REL) -> float:
    """q_S — varianza de proceso del precio, parametrizada de forma RELATIVA.

        q_S = (σ_rel · S)²                                [(USD/BTC)²]

    Sección C de la v1.2 del orden de trabajo. Entra en la diagonal de Q_k de la
    Sec. 7.3.3 del PDF, en la posición del estado S.

    ⚠ ESTO NO ES UNA ADIMENSIONALIZACIÓN DEL FILTRO, y la distinción importa. Se
    evaluó y se DESCARTÓ normalizar el EAKF completo: escalar Q y R sin escalar a
    la vez x, P, z y H no es un cambio de variables sino un cambio de modelo, y
    lleva el NIS a ~1e11 con Tr(P) colapsando a 1e-10 — el filtro declarando
    certeza casi perfecta mientras está equivocado. La transformación completa sí
    sería exacta, pero tampoco mejora el condicionamiento: cond(P₀) queda igual,
    porque el mal condicionamiento no viene de la escala sino de la asimetría de
    confianza deliberada de la Sec. 7.3.4 (p_S ≈ 0 contra p_v, p_Rn ≈ 1e6), que es
    información y no un defecto numérico.

    Aquí solo se reparametriza UNA entrada de Q con las mismas unidades que ya
    tenía. El álgebra del filtro es idéntica.
    """
    return (sigma_rel * S) ** 2


def gamma_omega(omega_m_max: float = OMEGA_M_MAX) -> float:
    """γ_ω — peso de ω_m en el escalar de expansión ρ_k.

        γ_ω = 1 / ω_m,max                                 [Ticks]

    Sección del PDF: 7.3.3 (ρ_k = 1 + γ_ω|ω_m| + γ_Q|ΣQ|).

    Con [ω_m] = 1/Ticks, γ_ω debe estar en Ticks para que el término sea
    adimensional, y normalizarlo por el máximo lo acota a ≤ 1.
    """
    if omega_m_max <= 0.0:
        raise ValueError(f"gamma_omega: ω_m,max debe ser positivo (recibido {omega_m_max})")
    return 1.0 / omega_m_max


def gamma_Q(suma_q_max: float = SUMA_Q_MAX, k_usd: float = K_USD) -> float:
    """γ_Q — peso de ΣQ en el escalar de expansión ρ_k.

        γ_Q = 1 / max( max(ΣQ_T̄), C_max )                 [1/USD]

    Sección del PDF: 7.3.3.

    ⚠ EL DENOMINADOR LLEVA max, NO min (Sec. 0.3 del orden de trabajo). ΣQ es
    volumen de MERCADO (Secs. 1.1 y 6.3), órdenes de magnitud mayor que C_max;
    tomar el mínimo del denominador *maximiza* γ_Q, que es lo contrario de acotar.
    Con max, el término queda ≤ 1 mientras ΣQ no supere su máximo histórico, y
    C_max actúa como piso del denominador durante el arranque en frío, cuando
    todavía no hay historial de volumen.
    """
    denominador = max(suma_q_max, k_usd)
    if denominador <= 0.0:
        raise ValueError("gamma_Q: el denominador max(ΣQ_max, C_max) debe ser positivo")
    return 1.0 / denominador


def kappa(
    r_base: float, omega_crit: float = OMEGA_CRIT, c: float = C_ESCALA_KAPPA
) -> float:
    """κ — factor de escalado de la aversión a la inestabilidad.

        κ = c · R_base / Ω_crit²                          [Ticks⁴/BTC²]

    Sección del PDF: 6.1 (R(Ω) = R_base + κ·Ω²).

    El PDF nunca declara κ — era el hueco original de CLAUDE.md.

    Anclaje a Ω_crit y no a un "Ω_max" (Sec. 0.3): Ω no tiene cota superior
    natural, porque por la Sec. 2.6 en el nodo de fase ΔS → 0 ⟹ Ψ → ∞ ⟹ Ω → ∞.
    Y eso es intencional: la Sec. 6.1 exige que R(Ω) diverja para forzar la
    inacción. El papel de κ no es normalizar a ≤ 1 sino fijar DÓNDE MUERDE la
    divergencia, y ese punto ya lo define Ω_crit:

        R(Ω_crit) = R_base + c·R_base = (1 + c)·R_base
    """
    if omega_crit <= 0.0:
        raise ValueError(f"kappa: Ω_crit debe ser positivo (recibido {omega_crit})")
    return c * r_base / (omega_crit * omega_crit)


def mu_inventario(
    q_base: float, omega_crit: float = OMEGA_CRIT, c_prima: float = C_ESCALA_MU
) -> float:
    """μ — sensibilidad al riesgo de exposición de inventario.

        μ = c' · q_base / Ω_crit²                         [Ticks⁴/BTC²]

    Sección del PDF: 6.2 (q_inv(Ω) = q_base + μ·Ω²).

    Mismo anclaje que κ: en Ω_crit la penalización por mantener inventario vale
    (1 + c')× su base, forzando el vaciado preventivo ANTES de que las
    restricciones de caja asimétricas se activen.

    `q_base` es el costo base de mantenimiento de inventario de la Sec. 6.2
    (NO la varianza base q_S del ruido de proceso de la Sec. 7.3.3, que comparte
    nombre y no tiene nada que ver).
    """
    if omega_crit <= 0.0:
        raise ValueError(f"mu_inventario: Ω_crit debe ser positivo (recibido {omega_crit})")
    return c_prima * q_base / (omega_crit * omega_crit)


# ==============================================================================
# 3. CONVERSIÓN DE UNIDADES DE ω_m (Sec. 0.4)
# ==============================================================================
def omega_m_desde_hz(f_hz: float, nu_ticks_por_anio: float) -> float:
    """Convierte la frecuencia del HHT de Hz (tiempo físico) a 1/Ticks.

        ω_m = f / ν        con ν expresada en Ticks/segundo

    Sección del PDF: 4.5 (definición de ν) y 2.5 (colapso espectral, que entrega f).

    ⚠ TRAMPA DE UNIDADES: la cadena EMD → Hilbert calcula f en Hz = 1/s, pero el
    modelo exige ω_m en 1/Ticks. Y ν se almacena en Ticks/AÑOS (que es la unidad
    que necesita c²_vol = k·ω_m·ν para quedar en 1/Años). Dividir f por ν tal cual
    da (1/s)·(Años/Ticks), que NO es 1/Ticks: hay que pasar ν a Ticks/segundo
    primero. Sin esta conversión, γ_ω·ω_m no es adimensional en ejecución aunque
    lo sea en el papel.
    """
    if nu_ticks_por_anio <= 0.0:
        return 0.0
    nu_por_segundo = nu_ticks_por_anio / SEGUNDOS_POR_ANIO  # [Ticks/s]
    if nu_por_segundo <= 0.0:
        return 0.0
    return f_hz / nu_por_segundo


# ==============================================================================
# 4. GUARDAS DE EJECUCIÓN (Sec. 0.5 y 0.6)
# ==============================================================================
class ErrorDimensional(AssertionError):
    """Un acoplamiento del modelo está mal escalado. Nombra la constante sospechosa."""


def verificar_cobertura_acotada(
    S: float, g0: float, delta_s_ref: float = DELTA_S_REF, i_max: float = I_MAX,
    tolerancia: float = 1e-6,
) -> None:
    """Guarda 1 (Sec. 0.5): la cobertura no puede exceder el límite de inventario.

        S · γ_0 · ΔS_max ≤ I_max

    Por construcción de `gamma_0` esto es una IGUALDAD, así que la guarda detecta
    que alguien reintrodujo un γ_0 literal o lo escaló mal.
    """
    cobertura = S * g0 * delta_s_ref
    if cobertura > i_max * (1.0 + tolerancia):
        # Mensajes en ASCII: la consola de Windows (cp1252) no puede codificar
        # letras griegas y un print fallido tumbaría el proceso.
        raise ErrorDimensional(
            f"gamma_0 MAL ESCALADO: la cobertura implicita es {cobertura:.4f} BTC "
            f"contra un limite I_max = {i_max:.4f} BTC "
            f"(S={S:.2f}, gamma_0={g0:.6e}, dS_max={delta_s_ref:.2f}). "
            f"Revisar si gamma_0 quedo con unidades de Delta (BTC^2/USD) en vez de "
            f"Gamma (BTC^3/USD^2) - ver constantes_micelio.gamma_0."
        )


def verificar_friccion_subcritica(S: float, lam: float, g0: float) -> None:
    """Guarda 2 (Sec. 0.5): al arranque no puede haber singularidad por construcción.

        λ · S² · γ_0 < 1

    Si esto falla en el primer ciclo, la singularidad de Loeper (Sec. 4.4.2) viene
    del escalado de las constantes, no de las condiciones de mercado.
    """
    producto = lam * S * S * g0
    if producto >= 1.0:
        raise ErrorDimensional(
            f"SINGULARIDAD DE LOEPER POR CONSTRUCCION: lambda*S^2*gamma_0 = "
            f"{producto:.4f} >= 1 en el primer ciclo (lambda={lam:.6e} /BTC, "
            f"S={S:.2f}, gamma_0={g0:.6e}). El denominador D = 1 - lambda*S^2*Gamma "
            f"nace no positivo. Sospechosos, en orden: gamma_0 "
            f"(constantes_micelio.gamma_0) y lambda (mu_OU / eta del bloque de ruido)."
        )


def verificar_nyquist_omega_m(omega_m_max: float = OMEGA_M_MAX) -> str | None:
    """Sanity check del piso teórico de γ_ω (Sec. 0.3).

    ω_m es la inversa de los ticks para completar un ciclo y el mínimo son 2 ticks
    (Nyquist en espacio de ticks), luego ω_m,max ≤ 0.5 /Ticks y γ_ω ≥ 2 Ticks.
    Devuelve un mensaje de advertencia, o None si está en rango.
    """
    if omega_m_max > 0.5:
        return (
            f"w_m,max = {omega_m_max:.4e} /Ticks viola el limite de Nyquist en espacio "
            f"de ticks (0.5 /Ticks = un ciclo cada 2 ticks). gamma_w = "
            f"{1.0/omega_m_max:.4f} Ticks queda por debajo de su piso teorico de "
            f"2 Ticks. El estimador HHT esta devolviendo algo fuera de rango."
        )
    return None


def verificar_estacionariedad_ou(
    sigma_ou: float, theta_ou: float, mu_ou: float, factor_holgura: float = 0.5
) -> str | None:
    """Guarda del proceso OU de liquidez (Sec. 0.6). ADVIERTE, NO CORRIGE.

    La desviación estacionaria de un Ornstein-Uhlenbeck es σ/√(2θ). Para que λ
    fluctúe ALREDEDOR de μ_OU y no lo sepulte hace falta σ/√(2θ) ≪ μ_OU.

    Las constantes de ruido (σ_OU, θ_OU, μ_OU, λ_min, η, β) están deliberadamente
    congeladas hasta tener datos de la cuenta demo (Sec. 0.6 del orden de trabajo),
    así que esta función NO ajusta nada: solo devuelve el mensaje de advertencia.
    """
    if theta_ou <= 0.0:
        return f"theta_OU = {theta_ou} no es positivo: el proceso OU no revierte a la media."
    sigma_estacionaria = sigma_ou / math.sqrt(2.0 * theta_ou)
    if mu_ou <= 0.0:
        return f"mu_OU = {mu_ou} no es positivo."
    if sigma_estacionaria > factor_holgura * mu_ou:
        return (
            f"sigma_OU/sqrt(2*theta_OU) = {sigma_estacionaria:.4e} no es << mu_OU = "
            f"{mu_ou:.4e} (razon {sigma_estacionaria/mu_ou:.1f}x). lambda no fluctua "
            f"alrededor de su media: pasea sobre ella y cruza a negativo, y el suelo "
            f"lambda_min sesga el resultado hacia arriba. NO SE CORRIGE AQUI - las "
            f"constantes de ruido esperan datos de la cuenta demo (Sec. 0.6)."
        )
    return None


# ==============================================================================
# 5. ADIMENSIONALES POR DECLARACIÓN (Sec. 0.7)
# ==============================================================================
# q_Δ, q_base (inventario) y R_base son ADIMENSIONALES. No es cosmético: es lo que
# hace que los tres términos de la función de costo J de la Sec. 6.1 sean
# conmensurables. Una vez e_k está en BTC (vía la conversión S·∂U/∂S de la
# Sec. 0.3), los tres sumandos de J quedan en BTC²·adimensional:
#
#     e_kᵀ Q e_k        →  [BTC²]·[adim]
#     U_kᵀ R(Ω) U_k     →  [BTC²]·[adim]     con R(Ω) = R_base + κΩ²,
#                                             [κ][Ω²] = (Ticks⁴/BTC²)(BTC²/Ticks⁴) ✓
#     q_inv(Ω)·I_k²     →  [adim]·[BTC²]     con q_inv = q_base + μΩ², igual ✓
#
# Sumarlos es legítimo sólo bajo esta declaración.
ADIMENSIONALES = ("q_Delta", "q_base_inventario", "R_base")


def resumen_constantes(
    S: float, r_base: float, q_base_inv: float, sigma_ou: float, theta_ou: float,
    mu_ou: float,
) -> str:
    """Vuelca las constantes derivadas y sus guardas. Para el arranque y los tests."""
    g0 = gamma_0(S)
    lineas = [
        "--- Constantes de acoplamiento derivadas (Sec. 0.1) ---",
        f"  gamma_0 = I_max/(S*dS_max)       = {g0:.6e}  [BTC^3/USD^2]  @ S={S:.2f}",
        f"  gamma_w = 1/w_m,max              = {gamma_omega():.6e}  [Ticks]",
        f"  gamma_Q = 1/max(SQ_max, C_max)   = {gamma_Q():.6e}  [1/USD]",
        f"  kappa   = c*R_base/W_crit^2      = {kappa(r_base):.6e}  [Ticks^4/BTC^2]",
        f"  mu      = c'*q_base/W_crit^2     = {mu_inventario(q_base_inv):.6e}  [Ticks^4/BTC^2]",
        "--- Limites estructurales (Sec. 0.2) ---",
        f"  I_max={I_MAX} BTC   C_max={K_USD:.0f} USD   dS_max={DELTA_S_REF} USD/BTC",
        f"  w_m,max={OMEGA_M_MAX:.4e} 1/Ticks   W_crit={OMEGA_CRIT}   SQ_max={SUMA_Q_MAX:.3e} USD",
        f"  cobertura implicita en dS_max = {S*g0*DELTA_S_REF:.4f} BTC (tope I_max={I_MAX})",
    ]
    for aviso in (
        verificar_nyquist_omega_m(),
        verificar_estacionariedad_ou(sigma_ou, theta_ou, mu_ou),
        verificar_dominio_malla(S),
    ):
        if aviso:
            lineas.append(f"  [!] {aviso}")
    return "\n".join(lineas)


def verificar_dominio_malla(
    S: float, margen_rel: float = MARGEN_MALLA_REL, delta_s_ref: float = DELTA_S_REF
) -> str | None:
    """El dominio de la malla debe cubrir el desplazamiento esperado del nodo de fase.

    Resuelve la colisión de nombres documentada arriba: el ΔS_max relativo de la
    Sec. 7.4.1 (ancho de malla) tiene que ser ≥ el ΔS_max absoluto de la Sec. 0.2
    (desplazamiento respecto al nodo), o el NMPC estará interpolando fuera de la
    región donde la cobertura es significativa.
    """
    ancho_absoluto = margen_rel * S
    if ancho_absoluto < delta_s_ref:
        return (
            f"El dominio de la malla (+/-{margen_rel:.1%} de S = +/-{ancho_absoluto:.1f} "
            f"USD) NO cubre el desplazamiento maximo esperado respecto al nodo de fase "
            f"(dS_max = {delta_s_ref:.1f} USD). El NMPC interpolaria fuera de la region "
            f"util. Subir MARGEN_MALLA_REL o bajar DELTA_S_REF."
        )
    return None
