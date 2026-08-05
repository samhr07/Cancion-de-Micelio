"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: hht.py — Transformada de Hilbert-Huang (Sección 2 del PDF)

Implementa la Prioridad 1 (Sección B) de ORDEN_TRABAJO_CALIBRACION_1.2.md:
la cadena EMD → Hilbert que hasta ahora vivía fuera de `Micelio.py` y cuya
ausencia bloqueaba las Fases 2 y 3.

POR QUÉ ES LA PRECONDICIÓN DE TODO
-----------------------------------
Mientras `R_n` fuera un coseno puro y el precio una sinusoide de 40 s, el NIS y
el Ljung-Box estaban midiendo el GENERADOR DE MOCKS, no el filtro. El NIS = 13.93
que reportó la fase anterior es real como medida, pero su causa era artificial:
`R_n` entra directamente en el vector de medición z_k (Sec. 7.3.1), así que un
`R_n` sintético contamina la innovación de la que dependen los tres tests.

CONTENIDO
---------
  2.1  Tamizado (sifting) y descomposición en IMFs.
  2.3  Mitigación del efecto de borde: extensión por reflexión + modelo AR.
  2.4  Señal analítica, desenrollado de fase y frecuencia instantánea.
  2.5  Colapso espectral: frecuencia ponderada por energía.
  2.6  Detección de nodos de fase.

DEPENDENCIAS: NumPy y SciPy (`signal.hilbert`, `interpolate.CubicSpline`).

Los mensajes que se imprimen van en ASCII (consola cp1252); los comentarios y
docstrings sí llevan acentos y símbolos.
"""

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.signal import hilbert

# Parámetros del tamizado. No son constantes de acoplamiento del modelo (no van
# en `constantes_micelio`): son tolerancias numéricas del algoritmo de la Sec. 2.1.
SD_UMBRAL = 0.2  # Criterio de parada de Huang: desviación estándar entre tamizados
MAX_ITER_TAMIZADO = 20  # Techo de iteraciones por IMF (presupuesto de latencia)
MAX_IMFS = 8  # Techo de modos extraídos
MIN_EXTREMOS = 3  # Debajo de esto el residuo se considera monótono
# Muestras del segmento terminal sobre el que se promedia la frecuencia
# instantánea. Ver la nota de robustez en `colapso_espectral`.
N_PROMEDIO_TERMINAL = 12


# ==============================================================================
# 2.3 MITIGACIÓN DEL EFECTO DE BORDE
# ==============================================================================
def _extender_por_reflexion(x: np.ndarray, m: int) -> tuple[np.ndarray, int]:
    """Extiende la serie por reflexión simétrica en ambos bordes.

    Sec. 2.3 del PDF: las envolventes construidas por splines cúbicos carecen de
    extremos locales futuros donde anclarse en el límite derecho del dominio (el
    instante presente t0). Si la Transformada de Hilbert se aplica sobre IMFs
    afectadas por esa divergencia, la fase instantánea y por tanto ω_m resultan
    erráticas justo en el último tick — que es el único que el control usa.

    Devuelve (serie_extendida, desplazamiento_del_origen).
    """
    n = len(x)
    m = int(min(m, n - 1))
    if m <= 0:
        return x.copy(), 0
    # Reflexión sobre el valor del extremo (no simple espejo del índice): mantiene
    # continuidad de valor y de primera derivada en la unión.
    izq = 2.0 * x[0] - x[m:0:-1]
    der = 2.0 * x[-1] - x[-2 : -m - 2 : -1]
    return np.concatenate([izq, x, der]), m


def _extrapolar_ar(x: np.ndarray, m: int, orden: int = 4) -> np.ndarray:
    """Predicción por modelo autoregresivo de bajo orden hacia el pseudo-futuro.

    Sec. 2.3: la extensión sintética combina reflexión simétrica con una predicción
    AR. La reflexión sola impone una simetría que puede no existir; el AR aporta la
    tendencia local. Se resuelve por mínimos cuadrados sobre la propia ventana.

    Si el ajuste sale inestable (raíces fuera del círculo unidad, o la predicción
    se dispara fuera del rango histórico), se devuelve `None` y el llamador cae a
    reflexión pura. Es preferible un borde conservador a uno que diverge.
    """
    n = len(x)
    if n < 4 * orden or m <= 0:
        return None
    # Matriz de diseño de un AR(orden) sobre la serie centrada.
    mu = float(x.mean())
    z = x - mu
    filas = n - orden
    A = np.empty((filas, orden), dtype=np.float64)
    for j in range(orden):
        A[:, j] = z[orden - 1 - j : n - 1 - j]
    b = z[orden:]
    try:
        coef, *_ = np.linalg.lstsq(A, b, rcond=None)
    except np.linalg.LinAlgError:
        return None
    # Estabilidad: raíces del polinomio característico dentro del círculo unidad.
    raices = np.roots(np.concatenate([[1.0], -coef]))
    if np.any(np.abs(raices) > 1.0):
        return None

    hist = list(z[-orden:])
    pred = []
    for _ in range(m):
        siguiente = float(np.dot(coef, hist[::-1][:orden]))
        pred.append(siguiente)
        hist.append(siguiente)
        hist.pop(0)
    pred = np.asarray(pred) + mu
    # Guarda anti-divergencia: la predicción no puede salirse del rango observado
    # más de media amplitud; si lo hace, el AR no es de fiar en esta ventana.
    amplitud = float(x.max() - x.min())
    if amplitud > 0 and (
        pred.max() > x.max() + 0.5 * amplitud or pred.min() < x.min() - 0.5 * amplitud
    ):
        return None
    return pred


def extender_senal(x: np.ndarray, m: int | None = None, usar_ar: bool = True):
    """Construye S_ext = [reflexión | señal | pseudo-futuro] de la Sec. 2.3.

    Devuelve (extendida, desplazamiento, longitud_original).
    """
    n = len(x)
    if m is None:
        # Un cuarto de la ventana basta para que el "aleteo" del spline ocurra en
        # el dominio sintético que después se descarta.
        m = max(4, n // 4)
    ext, desp = _extender_por_reflexion(x, m)
    if usar_ar and desp > 0:
        pred = _extrapolar_ar(x, desp)
        if pred is not None:
            ext[desp + n :] = pred  # El borde derecho (pseudo-futuro) pasa a ser AR
    return ext, desp, n


# ==============================================================================
# 2.1 TAMIZADO Y DESCOMPOSICIÓN EMPÍRICA DE MODOS
# ==============================================================================
def _indices_extremos(x: np.ndarray):
    """Índices de máximos y mínimos locales."""
    d = np.diff(x)
    # Se ignoran las mesetas: solo cambios estrictos de signo de la derivada.
    signo = np.sign(d)
    # Rellena ceros con el signo previo para no perder extremos en mesetas planas.
    for i in range(1, len(signo)):
        if signo[i] == 0:
            signo[i] = signo[i - 1]
    cambio = np.diff(signo)
    maximos = np.where(cambio < 0)[0] + 1
    minimos = np.where(cambio > 0)[0] + 1
    return maximos, minimos


def _envolventes(x: np.ndarray):
    """Envolventes superior e inferior por splines cúbicos (Sec. 2.1)."""
    n = len(x)
    maximos, minimos = _indices_extremos(x)
    if len(maximos) < 2 or len(minimos) < 2:
        return None, None

    # Se anclan los bordes al primer/último extremo para que el spline no dispare.
    ix_max = np.concatenate([[0], maximos, [n - 1]])
    iy_max = np.concatenate([[x[maximos[0]]], x[maximos], [x[maximos[-1]]]])
    ix_min = np.concatenate([[0], minimos, [n - 1]])
    iy_min = np.concatenate([[x[minimos[0]]], x[minimos], [x[minimos[-1]]]])

    ix_max, idx = np.unique(ix_max, return_index=True)
    iy_max = iy_max[idx]
    ix_min, idx = np.unique(ix_min, return_index=True)
    iy_min = iy_min[idx]
    if len(ix_max) < 4 or len(ix_min) < 4:
        return None, None

    t = np.arange(n)
    try:
        sup = CubicSpline(ix_max, iy_max)(t)
        inf = CubicSpline(ix_min, iy_min)(t)
    except ValueError:
        return None, None
    return sup, inf


def _tamizar(x: np.ndarray):
    """Extrae una IMF por tamizado iterativo (Sec. 2.1).

    Resta iterativamente la media m(T̄) de las envolventes hasta que la media sea
    ~cero, que es la condición de Función de Modo Intrínseco.
    """
    h = x.copy()
    for _ in range(MAX_ITER_TAMIZADO):
        sup, inf = _envolventes(h)
        if sup is None:
            return None
        media = 0.5 * (sup + inf)
        h_nuevo = h - media
        # Criterio de parada de Huang: desviación estándar normalizada.
        denom = np.sum(h * h)
        if denom <= 0.0:
            return h_nuevo
        sd = float(np.sum((h - h_nuevo) ** 2) / denom)
        h = h_nuevo
        if sd < SD_UMBRAL:
            break
    return h


def emd(x: np.ndarray, max_imfs: int = MAX_IMFS):
    """Descomposición Empírica de Modos. Devuelve (lista_de_IMFs, residuo).

    El residuo final monótono es R_n(T̄), la "marea" u offset dinámico de la
    analogía del barco de la Sec. 2.1.
    """
    x = np.asarray(x, dtype=np.float64)
    residuo = x.copy()
    imfs = []
    for _ in range(max_imfs):
        maximos, minimos = _indices_extremos(residuo)
        if len(maximos) + len(minimos) < MIN_EXTREMOS:
            break  # Residuo monótono: se acabó la descomposición
        imf = _tamizar(residuo)
        if imf is None:
            break
        imfs.append(imf)
        residuo = residuo - imf
    return imfs, residuo


# ==============================================================================
# 2.4 SEÑAL ANALÍTICA, FASE DESENROLLADA Y FRECUENCIA INSTANTÁNEA
# ==============================================================================
def fase_y_amplitud(imf: np.ndarray):
    """Z(T̄) = A(T̄)·e^{iθ(T̄)} vía Transformada de Hilbert (Sec. 2.4.1).

    Devuelve (amplitud_instantánea, fase_desenrollada).

    El desenrollado (Sec. 2.4.2) es obligatorio: el arcotangente devuelve la fase
    acotada en (−π, π], y al completarse el ciclo salta de π a −π. Derivar ese
    salto da una frecuencia negativa infinita — una singularidad matemática que
    inyectaría ruido en Ω y provocaría falsos bloqueos del NMPC.
    """
    z = hilbert(imf)
    amplitud = np.abs(z)
    fase = np.unwrap(np.angle(z))
    return amplitud, fase


def frecuencia_instantanea(fase: np.ndarray, dt: float) -> np.ndarray:
    """f(T̄) = (1/2π)·dθ_unwrap/dT̄   (Sec. 2.4.3).  Devuelve Hz si dt está en s."""
    return np.gradient(fase, dt) / (2.0 * np.pi)


# ==============================================================================
# 2.5 COLAPSO ESPECTRAL — frecuencia instantánea ponderada por energía
# ==============================================================================
def colapso_espectral(analiticas, n_prom: int = N_PROMEDIO_TERMINAL) -> float:
    """ω_m escalar consolidado, en las MISMAS unidades que 1/dt (Hz si dt en s).

        ω_m(T̄) = Σ_{i=2}^{n-1} A_i²(T̄)·f_i(T̄) / Σ_{i=2}^{n-1} A_i²(T̄)

    Sec. 2.5. Sumar las frecuencias de cada IMF directamente violaría el principio
    de superposición no lineal, así que se extrae el CENTROIDE ENERGÉTICO: la
    dinámica real del mercado está dictada por el ciclo que moviliza más capital,
    y eso se refleja en la amplitud.

    Se descartan sistemáticamente los extremos del espectro:
      - IMF1, cuyo comportamiento se asocia casi enteramente al ruido blanco de
        alta frecuencia del bid-ask spread.
      - El residuo R_n, que es la tendencia asintótica y carece de fase cíclica.

    `analiticas` es la lista [(A_i, θ_i, f_i), ...] YA TRUNCADA al dominio real,
    pero calculada sobre el vector extendido (ver `analizar_ventana`).

    ⚠ ROBUSTEZ: la Sec. 2.5 define ω_m(T̄) como una cantidad INSTANTÁNEA, y la
    lectura literal sería evaluarla en la última muestra. Medido, eso da un
    estimador inutilizable: f_i sale de derivar numéricamente la fase desenrollada
    y una sola muestra es dominada por el ruido de esa derivada — sobre la misma
    señal sintética, cambiar la semilla del ruido movía el error del 2 % al 45 %.
    Se toma la MEDIANA de f_i sobre las últimas `n_prom` muestras y la energía
    media sobre ese mismo segmento. La mediana, no la media, porque los saltos
    residuales del desenrollado son atípicos aislados, no ruido gaussiano.
    """
    if len(analiticas) < 2:
        return 0.0
    numerador = 0.0
    denominador = 0.0
    for amplitud, _fase, f in analiticas[1:]:  # [0] es IMF1: ruido de spread
        k = int(min(max(2, n_prom), len(f)))
        f_seg = f[-k:]
        f_seg = f_seg[np.isfinite(f_seg) & (f_seg > 0.0)]
        if f_seg.size == 0:
            continue
        f_i = float(np.median(f_seg))
        energia = float(np.mean(amplitud[-k:] ** 2))  # E_i = A_i²
        if not np.isfinite(f_i) or not np.isfinite(energia):
            continue
        numerador += energia * f_i
        denominador += energia
    if denominador <= 0.0:
        return 0.0
    return numerador / denominador


def concentracion_espectral(analiticas, n_prom: int = N_PROMEDIO_TERMINAL) -> float:
    """C ∈ [0,1]: fracción de la energía estructural que aporta la IMF dominante.

    Sec. D.4.2 de ORDEN_TRABAJO_RIESGO_1_3. Es el conmutador que decide si la
    matriz A del EAKF usa el oscilador armónico o velocidad constante:

        C ≥ C_ON   -> hay UN ciclo dominante nítido; A_arm(ω_ang) tiene sentido
        C ≤ C_OFF  -> la energía está repartida; ω_ang no describe la dinámica

    Se calcula sobre las mismas energías E_i = A_i² del colapso espectral de la
    Sec. 2.5 y con la misma exclusión de IMF1 (ruido de spread) y del residuo (sin
    fase cíclica). Reutilizar ese denominador no es economía: si C se midiera
    sobre un conjunto de modos distinto al que produce ω_m, estaría autorizando
    una frecuencia que no describe.

    ⚠ NO USAR C PARA ESCALAR ω. La tentación obvia —"si el ciclo es difuso, oscila
    menos"— sesga la frecuencia a la baja incluso cuando el ciclo es nítido,
    porque C nunca llega a 1. La decisión es binaria con histéresis.
    """
    if len(analiticas) < 2:
        return 0.0
    energias = []
    for amplitud, _fase, _f in analiticas[1:]:  # [0] es IMF1: ruido de spread
        k = int(min(max(2, n_prom), len(amplitud)))
        e = float(np.mean(amplitud[-k:] ** 2))
        if np.isfinite(e) and e > 0.0:
            energias.append(e)
    if not energias:
        return 0.0
    total = sum(energias)
    if total <= 0.0:
        return 0.0
    return max(energias) / total


def imf_dominante(analiticas, indice: int = -1) -> int:
    """Índice de la IMF estructural con más energía en `indice` (excluye IMF1)."""
    mejor, mejor_e = -1, -1.0
    for i, (amplitud, _fase, _f) in enumerate(analiticas):
        if i == 0:
            continue  # IMF1 = ruido de spread
        e = float(amplitud[indice] ** 2)
        if e > mejor_e:
            mejor, mejor_e = i, e
    return mejor


# ==============================================================================
# 2.6 NODOS DE FASE
# ==============================================================================
def detectar_nodo_fase(imf_dominante_serie: np.ndarray) -> bool:
    """¿El último instante cruzó un nodo de fase? (Sec. 2.6)

    NOTA DE INTERPRETACION — INCONSISTENCIA EN LA PROPIA SEC. 2.6.
    El texto define el nodo en prosa como el instante "en el cual la fase
    desenrollada de la dinámica dominante cruza el eje horizontal (amplitud
    cero)", y refuerza esa lectura más abajo: "El retorno del precio exactamente
    al nodo (S(T̄) → S_ref) fuerza a ΔS → 0". Pero la formaliza como

        θ_unwrap(T̄_nodo) ≡ 0   (mód π)

    y eso NO es amplitud cero. Con la convención estándar de la señal analítica
    Z = A·e^{iθ}, se tiene Re(Z) = IMF = A·cos θ, luego:
        θ ≡ 0   (mód π)  ⟺  IMF = ±A   -> la oscilación está en un EXTREMO
        θ ≡ π/2 (mód π)  ⟺  IMF = 0    -> la oscilación cruza CERO
    Es decir, la fórmula detecta los picos, no los cruces por cero: justo lo
    contrario de lo que pide la prosa y de lo que el modelo necesita.

    Se implementa lo que exige la FÍSICA del modelo, que es lo que la prosa
    describe: el nodo es el CRUCE POR CERO de la IMF dominante, porque ahí es
    donde el precio vuelve a su nivel de referencia y ΔS → 0. Además es mucho
    más estable numéricamente: el valor de la IMF en t0 está bien definido,
    mientras que el offset absoluto de la fase desenrollada cambia de una ventana
    deslizante a la siguiente.

    Corregir en el PDF: la condición formal debe ser θ ≡ π/2 (mód π), o bien
    reescribirse directamente como cruce por cero de la IMF dominante.
    """
    if len(imf_dominante_serie) < 2:
        return False
    a, b = float(imf_dominante_serie[-2]), float(imf_dominante_serie[-1])
    if not (np.isfinite(a) and np.isfinite(b)):
        return False
    return (a < 0.0 < b) or (b < 0.0 < a) or b == 0.0


# ==============================================================================
# API DE ALTO NIVEL — lo que consume el Hilo Lento
# ==============================================================================
def analizar_ventana(precios: np.ndarray, dt: float, usar_ar: bool = True) -> dict:
    """Cadena HHT completa sobre la ventana deslizante del Hilo Lento.

    Flujo algorítmico de la Sec. 2.3, en este orden exacto:
      1. Se extiende la serie S a S_ext.
      2. Se ejecuta EMD y Hilbert sobre el vector COMPLETO S_ext.
      3. Se truncan los tensores resultantes descartando los M puntos sintéticos,
         extrayendo las variables analíticas exactas para el instante verdadero t0.

    Esto fuerza a que el aleteo del spline ocurra en el dominio que se descarta, y
    entrega al Kalman y al NMPC métricas de fase limpias de distorsión de contorno.

    Devuelve un dict con:
        f_hz        frecuencia de mercado consolidada [Hz]  (convertir con
                    CTE.omega_m_desde_hz antes de usarla como ω_m en 1/Ticks)
        R_n         residuo macroeconómico en t0        [USD/BTC]
        nodo        True si t0 cruzó un nodo de fase (Sec. 2.6)
        fase        fase desenrollada de la IMF dominante, truncada
        n_imfs      número de IMFs extraídas
        valido      False si la ventana no dio para descomponer
    """
    x = np.asarray(precios, dtype=np.float64)
    n = len(x)
    vacio = {
        "f_hz": 0.0, "R_n": float(x[-1]) if n else 0.0, "nodo": False,
        "imf_dom_t0": 0.0, "fase": np.zeros(0), "n_imfs": 0, "valido": False,
        # C = 0 sin descomposición: sin evidencia de ciclo, la conmutación de la
        # Sec. D.4.2 debe quedarse en velocidad constante, que es la rama segura.
        "C": 0.0,
    }
    if n < 16:
        return vacio

    ext, desp, n_orig = extender_senal(x, usar_ar=usar_ar)
    imfs_ext, residuo_ext = emd(ext)
    if len(imfs_ext) < 2:
        return vacio

    ini, fin = desp, desp + n_orig

    # ORDEN CRÍTICO (Sec. 2.3, paso 2): la Transformada de Hilbert se ejecuta
    # sobre el vector COMPLETO S_ext, y solo DESPUÉS se truncan los tensores
    # resultantes. Aplicarla a las IMFs ya truncadas reintroduce un borde duro
    # exactamente en t0, que es justo el instante que el control usa: medido,
    # eso daba errores de frecuencia del 200-600 %.
    analiticas = []
    for imf in imfs_ext:
        amplitud, fase = fase_y_amplitud(imf)
        f = frecuencia_instantanea(fase, dt)
        analiticas.append((amplitud[ini:fin], fase[ini:fin], f[ini:fin]))

    imfs = [imf[ini:fin] for imf in imfs_ext]
    residuo = residuo_ext[ini:fin]

    f_hz = colapso_espectral(analiticas)

    # Nodo de fase sobre la IMF dominante (la de mayor energía, excluyendo IMF1).
    i_dom = imf_dominante(analiticas, indice=-1)
    if i_dom >= 0:
        fase_dom = analiticas[i_dom][1]
        nodo = detectar_nodo_fase(imfs[i_dom])
    else:
        fase_dom = np.zeros(0)
        nodo = False

    return {
        "f_hz": float(f_hz),
        "R_n": float(residuo[-1]),
        "nodo": bool(nodo),
        # Valor de la IMF dominante en t0. El consumidor debería preferir ESTE
        # campo al booleano `nodo` para detectar el cruce: comparar imf[-2] contra
        # imf[-1] dentro de una misma ventana depende del borde, que se recalcula
        # en cada llamada; seguir el SIGNO en t0 entre llamadas sucesivas usa un
        # único valor por ventana y es mucho más estable. Ver el Hilo Lento.
        "imf_dom_t0": float(imfs[i_dom][-1]) if i_dom >= 0 else 0.0,
        "fase": fase_dom,
        "n_imfs": len(analiticas),
        # Concentración espectral para la conmutación de la rama de A (Sec. D.4.2
        # de la v1.3). Se devuelve aquí y no se recalcula en el Hilo Lento porque
        # depende de las mismas `analiticas` que ya se descompusieron.
        "C": float(concentracion_espectral(analiticas)),
        "valido": True,
    }
