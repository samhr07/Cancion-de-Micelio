"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: dinamica.py — Matriz de transición A del EAKF y protocolo A/B

Implementa las Secciones D y E de ORDEN_TRABAJO_RIESGO_1_3.md.

QUÉ RESUELVE
------------
La matriz `A` del EAKF (Sec. 7.3.5 del PDF) asume VELOCIDAD CONSTANTE sobre un
precio que el propio modelo describe como oscilatorio (Sec. 1). Medido en la v1.2
(§1.6), con un 10 % de error en ω — el caso realista, no el ideal:

    | Modelo                        | NIS  | ρ₁     |
    |-------------------------------|------|--------|
    | Velocidad constante (actual)  | 5.24 | +0.841 |
    | Sinusoide prescrita, ω −10 %  | 2.57 | +0.591 |
    | Oscilador armónico, ω −10 %   | 1.43 | −0.050 |

Estaba aparcado en la v1.2 "hasta datos de cuenta demo" porque validarlo contra
un mock sinusoidal que pusimos nosotros no probaba nada. Con el Modo LECTURA
disponible se levanta el aparcamiento.

⚠ ESTE MÓDULO EXISTE PORQUE LOS DOS ERRORES QUE ACECHAN AQUÍ FALLAN EN SILENCIO.
Ni la trampa del 2π (Sec. D.1) ni la de las unidades de Δt (Sec. D.2) producen
excepción, `nan` ni log: producen un modelo que "no aporta", y la conclusión
equivocada sería que la propuesta no sirve. Aislarlos en funciones puras es lo
que permite testearlos (Sec. F).
"""

from __future__ import annotations

import math

import numpy as np

# Códigos de la rama activa, registrados en telemetría (Sec. D.4.3) para poder
# condicionar el análisis del A/B sobre ella.
RAMA_VELOCIDAD_CONSTANTE = 0
RAMA_ARMONICO = 1

# Umbral de la rama de Taylor. Por debajo de |ωΔt| = 1e-6, sin(ωΔt)/ω pierde
# dígitos significativos por cancelación catastrófica y en ω = 0 exacto es 0/0.
UMBRAL_TAYLOR = 1.0e-6

# ⚠ El jacobiano necesita un umbral MAYOR que `A_arm`, y el valor está MEDIDO.
# `sin(ω)/ω` pierde poca precisión cerca de cero (numerador y denominador son
# ambos O(ω)), pero su derivada `(ω·cos ω − sin ω)/ω²` resta dos cantidades casi
# iguales de orden ω para dar un resultado de orden ω³: cancelación catastrófica.
# La rama de Taylor tiene el problema opuesto — su error de truncamiento crece
# con ω. El umbral óptimo es donde ambos se cruzan, y eso no se razona, se mide:
#
#     ω       error relativo de la fórmula exacta   error relativo de Taylor
#     1e-6            7.8e-05                              1.0e-13
#     1e-5            1.4e-06                              1.0e-11
#     1e-4            1.1e-08                              1.0e-09
#     3e-4            ~1e-08                               ~3e-08   <- cruce
#     1e-3            5.6e-11                              1.0e-07
#     1e-2            8.5e-13                              1.0e-05
#
# Se toma 3e-4: el peor de los dos errores queda en ~1e-8 relativo. Un intento
# previo de fijarlo en 1e-3 "porque ahí se cruzan" era razonamiento sin medición,
# y dejaba un salto del 0.2 % en el umbral.
UMBRAL_TAYLOR_JACOBIANO = 3.0e-4


# ==============================================================================
# D.1 y D.2 — LAS DOS TRAMPAS DE UNIDADES
# ==============================================================================
def omega_angular_desde_hz(f_hz: float) -> float:
    """ω_ang [rad/s] = 2π·f [Hz].  Sec. D.1 y D.2.

    ⚠ TRAMPA DEL 2π (Sec. D.1). `hht.frecuencia_instantanea` divide por 2π, así
    que devuelve `f` en Hz = ciclos/s. Y `constantes_micelio.omega_m_desde_hz`
    devuelve ω_m = f/ν en ciclos/tick. Ambas son frecuencias ORDINARIAS.
    Pero `A_arm` sale de s̈ = −ω²s, cuya solución es cos(ωt): ahí ω es ANGULAR.
    Usar ω_m o f directamente introduce un factor 2π ≈ 6.28 de error, o sea un
    628 %. Con solo −30 % de error el armónico ya se degrada a ρ₁ = +0.544, así
    que con 2π sería PEOR que velocidad constante y concluirías, equivocadamente,
    que la propuesta no sirve.

    ⚠ TRAMPA DE LAS UNIDADES DE Δt (Sec. D.2). ω_m está en 1/Ticks y el Δt del
    Hilo Rápido en segundos: ω_m·Δt NO es adimensional. Con ω_m ≈ 1.26e-3 y
    Δt = 0.01 s daría 1.26e-5 en vez del 1.57e-3 correcto — un factor 125 — y
    `A_arm` degeneraría a velocidad constante SIN AVISAR.

    Resolución (Sec. D.2): dos variables de frecuencia distintas para dos usos
    distintos, publicadas ambas, sin convertir en el punto de uso.

        ω_m   [1/Ticks]  ->  ρ_k (7.3.3) y c²_vol = k·ω_m·ν (4.5).  SIN CAMBIOS.
        ω_ang [rad/s]    ->  A_arm, y SOLO A_arm.

    Magnitud de control: ciclo de 40 s -> f = 0.025 Hz -> ω_ang = 0.157 rad/s ->
    ω_ang·Δt = 1.57e-3 con Δt = 0.01 s.
    """
    if not math.isfinite(f_hz) or f_hz <= 0.0:
        return 0.0
    return 2.0 * math.pi * f_hz


def periodo_implicito_ticks(A: np.ndarray) -> float:
    """Período en TICKS que la matriz A codifica de verdad.  §4.1 / §8.

    Igual que `periodo_implicito` pero con Δn = 1, que es el paso del filtro en
    reloj de transacciones. Sigue siendo la ÚNICA defensa contra la trampa del
    2π, que en espacio de ticks reaparece intacta: `omega_m_desde_hz` devuelve
    ciclos/tick (ordinaria) y `A_arm` necesita rad/tick (angular).

    Para pasar a segundos: `periodo_s = periodo_ticks / nu_ticks_por_s`.
    """
    return periodo_implicito(A, 1.0)


def periodo_implicito(A: np.ndarray, dt: float) -> float:
    """Período [s] que la matriz A codifica de verdad. Instrumento de test.

    Es la ÚNICA defensa contra D.1 y D.2, que no producen ningún síntoma
    observable. Se invierte la construcción de `matriz_A_armonica`:

        A[0,0] = cos(ωΔt)        A[0,1] = sin(ωΔt)/ω
        A[1,0] = −ω·sin(ωΔt)

    luego  sin²(ωΔt) = A[0,1]·(−A[1,0])  y  ωΔt = atan2(|sin|, cos).

    Devuelve inf si la matriz es de velocidad constante (ω = 0), que es el
    período de un ciclo que nunca vuelve.
    """
    cos_wdt = float(A[0, 0])
    sin2 = float(A[0, 1]) * (-float(A[1, 0]))
    sin_wdt = math.sqrt(max(0.0, sin2))
    wdt = math.atan2(sin_wdt, cos_wdt)
    if wdt <= 0.0 or dt <= 0.0:
        return float("inf")
    w = wdt / dt
    return 2.0 * math.pi / w


# ==============================================================================
# D.3 — MATRICES DE TRANSICIÓN
# ==============================================================================
def matriz_A_velocidad_constante(dt: float) -> np.ndarray:
    """A actual (Sec. 7.3.5): posición-velocidad con velocidad constante.

        x = [S, v, R_n]ᵀ,   S_{k+1} = S_k + v_k·Δt
    """
    return np.array(
        [[1.0, dt, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64
    )


def matriz_A_armonica(w_ang: float, dt: float) -> np.ndarray:
    """A del oscilador armónico. Sec. D.3.

            ⎡  cos(ωΔt)      sin(ωΔt)/ω    0 ⎤
    A_arm = ⎢ -ω·sin(ωΔt)    cos(ωΔt)      0 ⎥     ω ≡ ω_ang [rad/s], Δt [s]
            ⎣  0             0             1 ⎦

    Es la solución exacta de s̈ = −ω²s discretizada, no una aproximación: para el
    modo dominante, propaga el estado sin error de truncamiento.

    La tercera fila/columna deja `R_n` intacto, exactamente igual que en la matriz
    de velocidad constante.

    ⚠ RAMA DE TAYLOR OBLIGATORIA. sin(ωΔt)/ω es 0/0 cuando ω → 0. Sin ella, un
    régimen sin ciclo dominante produce `nan`, y un `nan` en A contamina P de
    forma IRREVERSIBLE: P_pred = A·P·Aᵀ + Q propaga el nan a las nueve entradas y
    ninguna medición posterior lo limpia. El filtro no se recupera nunca.
        |ωΔt| < 1e-6:   sin(ωΔt)/ω ≈ Δt·(1 − (ωΔt)²/6)
                        −ω·sin(ωΔt) ≈ −ω²Δt
    Con ω = 0 exacto esto reproduce [[1, Δt], [0, 1]] BIT A BIT, que es lo que
    exige el test de equivalencia de la Sec. F.
    """
    w = float(w_ang)
    if not math.isfinite(w) or not math.isfinite(dt):
        return matriz_A_velocidad_constante(dt if math.isfinite(dt) else 0.0)

    wdt = w * dt
    if abs(wdt) < UMBRAL_TAYLOR:
        cos_wdt = 1.0 - 0.5 * wdt * wdt
        sinc = dt * (1.0 - wdt * wdt / 6.0)  # sin(ωΔt)/ω
        menos_w_sin = -w * w * dt  # −ω·sin(ωΔt)
        if w == 0.0:
            # Reproducción exacta, sin residuos de coma flotante.
            cos_wdt, sinc, menos_w_sin = 1.0, dt, 0.0
    else:
        cos_wdt = math.cos(wdt)
        sinc = math.sin(wdt) / w
        menos_w_sin = -w * math.sin(wdt)

    return np.array(
        [
            [cos_wdt, sinc, 0.0],
            [menos_w_sin, cos_wdt, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )


def jacobiano_A_respecto_omega(w_ang_rad_tick: float) -> np.ndarray:
    """J = ∂A_arm/∂ω evaluada en Δn = 1.  §4.3 de la v2.1.

        ∂cos(ω)/∂ω      = −sin(ω)
        ∂(sin(ω)/ω)/∂ω  = (ω·cos(ω) − sin(ω))/ω²
        ∂(−ω·sin(ω))/∂ω = −sin(ω) − ω·cos(ω)

    ⚠ RAMA DE TAYLOR OBLIGATORIA, misma disciplina que `matriz_A_armonica` pero
    con umbral PROPIO Y MAYOR (ver `UMBRAL_TAYLOR_JACOBIANO`): el término central
    es 0/0 en ω = 0 y sufre cancelación catastrófica mucho antes de llegar ahí.
    Para |ω| < UMBRAL_TAYLOR_JACOBIANO (medido en 3e-4, ver la constante):
        −sin(ω)                ≈ −ω
        (ω·cos ω − sin ω)/ω²   ≈ −ω/3
        −sin ω − ω·cos ω       ≈ −2ω

    Los tres tienden a CERO cuando ω → 0, luego `Q_omega → 0`: **sin ciclo, no
    hay incertidumbre de ciclo.** Esa es la comprobación de coherencia de toda la
    construcción, y tiene test propio.
    """
    w = float(w_ang_rad_tick)
    if not math.isfinite(w):
        return np.zeros((3, 3), dtype=np.float64)
    if abs(w) < UMBRAL_TAYLOR_JACOBIANO:
        d_cos = -w
        d_sinc = -w / 3.0
        d_menos_w_sin = -2.0 * w
    else:
        d_cos = -math.sin(w)
        d_sinc = (w * math.cos(w) - math.sin(w)) / (w * w)
        d_menos_w_sin = -math.sin(w) - w * math.cos(w)
    return np.array(
        [
            [d_cos, d_sinc, 0.0],
            [d_menos_w_sin, d_cos, 0.0],
            [0.0, 0.0, 0.0],  # R_n no depende de ω
        ],
        dtype=np.float64,
    )


def Q_omega(w_ang_rad_tick: float, x: np.ndarray, sigma_omega: float,
            S_ref: float = 0.0) -> np.ndarray:
    """Ruido de proceso por INCERTIDUMBRE EN ω.  §4.3 de la v2.1.

        Q_ω = (J·x)(J·x)ᵀ · σ_ω²        con  J = ∂A_arm/∂ω

    `DIVERGE DEL PDF (Sec. 7.3.3)`. Hoy `ρ_k = 1 + γ_ω·|ω_m| + γ_Q·|ΣQ|` infla `Q`
    con la MAGNITUD de ω — dice "los ciclos rápidos son más inciertos". Lo que se
    quiere es que infle con la INCERTIDUMBRE — "los ciclos mal conocidos son más
    inciertos". Con ω promovida a determinante de la dinámica en `A_arm`, lo
    segundo es lo físicamente correcto, y además es propagación de incertidumbre
    paramétrica derivable **sin ninguna constante libre**.

    Se evalúa sobre la DESVIACIÓN `x − x_ref`, coherente con la forma afín de la
    predicción (Sec. D.3 de la v1.3): es esa desviación la que `A_arm` propaga.

    ⚠ NO SUSTITUYE a `ρ_k` todavía. El §4.3 pide implementar AMBAS variantes,
    registrarlas en paralelo y **no decidir aquí**: decide el dato.
    """
    if not math.isfinite(sigma_omega) or sigma_omega <= 0.0:
        return np.zeros((3, 3), dtype=np.float64)
    J = jacobiano_A_respecto_omega(w_ang_rad_tick)
    x_ref = np.array([[S_ref], [0.0], [0.0]], dtype=np.float64)
    v = J @ (np.asarray(x, dtype=np.float64).reshape(3, 1) - x_ref)
    return (v @ v.T) * (sigma_omega * sigma_omega)


def acumular_Q(A_un_paso: np.ndarray, Q_tick: np.ndarray, n_pasos: int) -> np.ndarray:
    """Ruido de proceso acumulado sobre `n_pasos` de Δn = 1.  §4.2 / §8.

        Q_N = Σ_{i=0}^{N-1}  Aⁱ · Q_tick · (Aⁱ)ᵀ

    ⚠ NO es `N · Q_tick`, y la diferencia importa. Solo coinciden si A = I: en
    cuanto A propaga (velocidad hacia posición, u oscilación), el ruido inyectado
    en el paso i se transporta i pasos más y su contribución a la covarianza
    final crece. Suponer `N·Q` sobreestima la certeza en la posición.

    Esta función existe sobre todo para el TEST DE INVARIANCIA A LA TASA del §8:
    N pasos de Δn=1 deben dar el mismo `x` y la misma `P` que un paso de Δn=N con
    esta `Q_N`. Ese test es el que demuestra que §1.3 quedó cerrado — que el ruido
    de proceso dejó de depender de cuántas veces despertó el planificador.
    """
    Q_acum = np.zeros_like(Q_tick)
    Ai = np.eye(A_un_paso.shape[0])
    for _ in range(int(n_pasos)):
        Q_acum = Q_acum + Ai @ Q_tick @ Ai.T
        Ai = A_un_paso @ Ai
    return Q_acum


def predecir_afin(A: np.ndarray, x: np.ndarray, P: np.ndarray, Q: np.ndarray, S_ref: float):
    """Paso de predicción del EAKF en FORMA AFÍN. Sec. D.3.

        x_ref  = [S_ref, 0, 0]ᵀ
        x_pred = x_ref + A·(x_k − x_ref)
        P_pred = A·P_k·Aᵀ + Q_k          (un offset no afecta a la covarianza)

    POR QUÉ AFÍN Y NO SOBRE s = S − S_ref
    -------------------------------------
    La v1.2 formulaba `A_arm` sobre la desviación s. Hacerlo literalmente cambia
    el significado de x[0] de precio absoluto a desviación, y obliga a auditar
    TODOS los consumidores: q_S = (σ_rel·S)², γ_0(S), el centro de la malla de
    Loeper, z₀ = P_spot, la telemetría, el NMPC. Es una superficie de bug grande
    y silenciosa. La forma afín es algebraicamente equivalente y conserva
    x[0] = S absoluto.

    BENEFICIO ADICIONAL: el manejo de nodos de fase sale gratis. La v1.2 exigía
    `s ← s + (S_ref_viejo − S_ref_nuevo)` con P sin tocar. En coordenadas
    absolutas eso es exactamente NO HACER NADA: S_ref salta, el offset se aplica
    en la predicción siguiente, y no hay discontinuidad ni en x ni en P. Se
    elimina toda una clase de bugs, incluida la cascada
    innovación espuria -> pico de NIS -> ventana del EMD a W_min -> racha de
    burn-in rota.

    La tercera componente de x_ref es 0, así que R_n pasa intacto por la fila
    [0, 0, 1] tanto en la rama armónica como en la de velocidad constante.
    """
    x_ref = np.array([[S_ref], [0.0], [0.0]], dtype=np.float64)
    x_pred = x_ref + A @ (x - x_ref)
    P_pred = A @ P @ A.T + Q
    return x_pred, P_pred


# ==============================================================================
# D.4 — MITIGACIONES DEL RIESGO QUE INTRODUCE A_arm
# ==============================================================================
class ConmutadorRamaA:
    """Elige la rama de A por concentración espectral, con histéresis. Sec. D.4.2.

    `A_arm` asciende ω_m de modulador de Q a DETERMINANTE DE LA DINÁMICA DEL
    ESTADO. Hoy un ω_m ruidoso solo ensancha la incertidumbre; con el armónico
    haría ruidosa la transición misma. Las tres mitigaciones del documento:

      1. EMA sobre f_hz en el Hilo Lento antes de derivar ω_ang (vive allí: la
         mediana sobre las últimas 12 muestras ya existe en `hht.py`, la EMA va
         encima).
      2. Esta conmutación: C ≥ C_ON -> armónico; C ≤ C_OFF -> velocidad
         constante. C_OFF < C_ON mata el chatter.
      3. `rama_A` en telemetría.

    ⚠ NO ESCALAR ω POR C. Sería la tentación obvia ("si el ciclo es difuso, usa
    menos oscilación") y es incorrecta: sesga la frecuencia a la baja INCLUSO
    cuando el ciclo es nítido, porque C nunca llega a 1 exactamente. La decisión
    es binaria con histéresis, y ω se usa tal cual o no se usa.
    """

    def __init__(self, c_on: float, c_off: float):
        if c_off >= c_on:
            raise ValueError(
                f"ConmutadorRamaA: C_OFF={c_off} debe ser < C_ON={c_on} o no hay "
                f"histeresis y la rama oscila con el ruido de C (chatter)."
            )
        self.c_on = float(c_on)
        self.c_off = float(c_off)
        self.rama = RAMA_VELOCIDAD_CONSTANTE  # Arranque conservador
        self.n_conmutaciones = 0

    def actualizar(self, concentracion: float, w_ang: float) -> int:
        """Devuelve la rama vigente. `w_ang <= 0` fuerza velocidad constante.

        La condición sobre w_ang no es redundante con la de C: durante el warm-up
        del EMD (~96 s) todavía no hay estimación de frecuencia, y C podría venir
        de una descomposición degenerada.
        """
        c = float(concentracion) if math.isfinite(concentracion) else 0.0
        if w_ang <= 0.0 or not math.isfinite(w_ang):
            nueva = RAMA_VELOCIDAD_CONSTANTE
        elif self.rama == RAMA_ARMONICO:
            nueva = RAMA_VELOCIDAD_CONSTANTE if c <= self.c_off else RAMA_ARMONICO
        else:
            nueva = RAMA_ARMONICO if c >= self.c_on else RAMA_VELOCIDAD_CONSTANTE
        if nueva != self.rama:
            self.n_conmutaciones += 1
            self.rama = nueva
        return self.rama

    def matriz(self, dt: float, w_ang: float) -> np.ndarray:
        if self.rama == RAMA_ARMONICO:
            return matriz_A_armonica(w_ang, dt)
        return matriz_A_velocidad_constante(dt)


# ==============================================================================
# E. PROTOCOLO A/B — DOS EAKF SOBRE EL MISMO FLUJO DE MEDICIONES
# ==============================================================================
class EAKFSombra:
    """EAKF de contraste que corre en paralelo sin alimentar el control. Sec. E.1.

    Es la ÚNICA comparación honesta: dos corridas distintas verían mercados
    distintos, y sobre un mercado real no se puede repetir el experimento. Son
    matrices 3×3, así que el costo es despreciable contra los 2.2 ms/ciclo
    medidos de Loeper+NMPC (criterio de aceptación: sobrecosto < 0.3 ms/ciclo).

    El filtro sombra replica el de producción salvo en la matriz A. Comparte
    exactamente el mismo z_k, el mismo R_k y el mismo Q_k: cualquier otra
    diferencia contaminaría la comparación.
    """

    def __init__(self, x0: np.ndarray, P0: np.ndarray, H: np.ndarray, usa_armonico: bool):
        self.x = np.array(x0, dtype=np.float64, copy=True)
        self.P = np.array(P0, dtype=np.float64, copy=True)
        self.H = H
        self.usa_armonico = bool(usa_armonico)
        self.I3 = np.eye(3)
        self.innov = np.zeros((2, 1))
        self.nis = float("nan")

    def paso(self, z, R_k, Q_k, dt, w_ang, S_ref, hay_medicion: bool, H_ef=None):
        """Un ciclo completo. Devuelve (innovacion, nis).

        `H_ef` permite la ACTUALIZACION SECUENCIAL MULTI-TASA: cuando `R_n` no
        es fresco, el filtro de produccion corrige solo con la fila del precio, y
        el sombra tiene que hacer EXACTAMENTE lo mismo. Si el sombra asimilara
        siempre las dos componentes, la comparacion del A/B dejaria de aislar la
        matriz `A` — que es su unico proposito.
        """
        H_ef = self.H if H_ef is None else H_ef
        A = (
            matriz_A_armonica(w_ang, dt)
            if self.usa_armonico
            else matriz_A_velocidad_constante(dt)
        )
        x_pred, P_pred = predecir_afin(A, self.x, self.P, Q_k, S_ref)

        if not hay_medicion:
            self.x, self.P = x_pred, P_pred
            self.innov = np.zeros((z.shape[0], 1))
            self.nis = float("nan")
            return self.innov, self.nis

        innov = z - (H_ef @ x_pred)
        S_cov = H_ef @ P_pred @ H_ef.T + R_k
        try:
            K = (P_pred @ H_ef.T) @ np.linalg.solve(S_cov, np.eye(S_cov.shape[0]))
            self.nis = float(innov.T @ np.linalg.solve(S_cov, innov))
        except np.linalg.LinAlgError:
            self.x, self.P = x_pred, P_pred
            self.innov, self.nis = innov, float("nan")
            return self.innov, self.nis

        self.x = x_pred + K @ innov
        joseph = self.I3 - K @ H_ef
        self.P = (joseph @ P_pred @ joseph.T) + (K @ R_k @ K.T)  # Forma de Joseph
        self.innov = innov
        return self.innov, self.nis
