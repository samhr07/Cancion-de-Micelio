# -*- coding: utf-8 -*-
"""
ssa.py -- Analisis Espectral Singular (SSA) como sustituto de la cadena EMD -> Hilbert.

POR QUE SE CAMBIA DE HERRAMIENTA
--------------------------------
La v2.2 dejo demostrado que la EMD **no tiene hipotesis nula**: sobre un paseo
aleatorio puro devuelve un "periodo" de 118.1 s que se parte a 37.5 s al partir la
ventana en dos. El modo dominante lo fija la ventana, no la senal, asi que `omega_m`
era salida del algoritmo y no del mercado.

SSA no arregla eso por si solo -- tambien puede inventar componentes suaves sobre
ruido rojo. Lo que si tiene, y la EMD no, es un **test nulo estandar**: Monte Carlo
SSA (Allen & Smith 1996), que compara los autovalores observados contra los que
produce un AR(1) ajustado a los mismos datos. Por eso el cambio vale la pena: no
porque descomponga mejor, sino porque se puede FALSAR.

ESTRUCTURA DEL METODO
---------------------
1. Encaje: matriz de trayectoria X (L x K), K = N - L + 1, X[i,j] = x[i+j].
2. Descomposicion: SVD de X. Se obtiene por autodescomposicion de la matriz de Gram
   G = X X^T (L x L), que se construye en O(L*N) con sumas prefijas -- nunca se
   materializa X, que a L=1024 y N=8192 serian 58 MB.
3. Series elementales: la media diagonal de sigma_i * U_i V_i^T es exactamente
   conv(U_i, X^T U_i) / cuentas, porque la suma sobre la antidiagonal i+j=k de un
   producto exterior ES la convolucion. Sale exacta y barata.
4. Agrupamiento: se decide con la matriz de w-correlacion.

SOBRE LA CONDICION DE PARADA DEL BARRIDO DE L
---------------------------------------------
Buscar `ortogonalidad == 0` no termina nunca: el ruido del sensor (el mismo que
alimenta R en el EAKF) siempre filtra energia entre componentes, asi que el minimo
alcanzable es estrictamente positivo y desconocido de antemano. Aqui el barrido es
una **rejilla finita** -- no hay bucle -- y la eleccion es por **minimo local**
(`elegir_minimo_local`), con dos salvaguardas:
  - si el minimo cae en un extremo de la rejilla, se reporta que el criterio NO
    discrimina en vez de devolver el extremo como si fuera una eleccion;
  - se acompana siempre del mismo estadistico medido sobre sustitutos barajados,
    porque una metrica que baja igual sobre datos sin estructura no esta midiendo
    separacion, esta midiendo la rejilla.

Convenciones del proyecto: comentarios en espanol; todo texto que se IMPRIME va en
ASCII (la consola es cp1252); las suposiciones que rellenan huecos van marcadas con
`# NOTA DE INTERPRETACION:`.
"""

from __future__ import annotations

import numpy as np
from scipy import signal as sig_sp
from scipy import linalg as lin_sp


# ---------------------------------------------------------------------------
# Umbrales de decision.
#
# NOTA DE INTERPRETACION: ninguno de estos numeros sale del PDF ni de los ordenes
# de trabajo -- SSA entra al proyecto en esta sesion. Son los criterios habituales
# de la literatura de SSA (Golyandina & Zhigljavsky) y estan aqui, agrupados y con
# nombre, precisamente para que se puedan recalibrar sin buscarlos por el codigo.
# ---------------------------------------------------------------------------

# Un par oscilatorio genuino tiene autovalores casi iguales (seno y coseno de la
# misma frecuencia llevan la misma energia).
RAZON_LAMBDA_PAR = 0.80

# Las dos frecuencias dominantes del par deben caer dentro de un bin de Rayleigh
# (1/L ciclos por muestra), que es la resolucion real del autovector.
BINS_RAYLEIGH_PAR = 1.0

# Fraccion minima de la energia del par concentrada en su frecuencia comun.
CONCENTRACION_PAR = 0.60

# Cuadratura: el desfase mediano entre los dos autovectores debe ser ~90 grados.
# 0.35 rad = 20 grados de tolerancia.
TOLERANCIA_CUADRATURA = 0.35

# Clasificacion del color de ruido por la pendiente beta de S(f) ~ f^(-beta).
# La banda de beta negativa hace falta y no es academica: el residuo que deja el
# SSA tiene las bajas frecuencias extraidas, asi que su espectro CRECE con f y
# beta sale muy negativa. Sin esta banda se etiquetaba como "BLANCO", que es
# justo lo contrario de lo que mide.
BANDAS_COLOR = (
    (-np.inf, -0.5, "AZUL (crece con f)"),
    (-0.5, 0.5, "BLANCO"),
    (0.5, 1.5, "ROSA (1/f)"),
    (1.5, 2.5, "ROJO (browniano)"),
    (2.5, np.inf, "NEGRO"),
)

# Por debajo de este R2 la ley de potencias no describe el espectro y la etiqueta
# de color no significa nada.
R2_MINIMO_LEY_POTENCIAS = 0.70


# ===========================================================================
# 1. Descomposicion
# ===========================================================================

def matriz_gram(x: np.ndarray, L: int) -> np.ndarray:
    """G = X X^T con X la matriz de trayectoria (L x K), sin materializar X.

    G[i, i+m] = sum_{k=0}^{K-1} x[i+k] * x[i+m+k]. Para cada desfase m eso es una
    suma deslizante de longitud K sobre el producto x[n]*x[n+m], o sea una
    diferencia de sumas prefijas. Coste O(L*N) en vez de O(L^2*K).
    """
    x = np.asarray(x, dtype=np.float64)
    N = x.size
    K = N - L + 1
    if K < L:
        raise ValueError("L=%d es demasiado grande para N=%d (hace falta K >= L)" % (L, N))

    G = np.empty((L, L), dtype=np.float64)
    for m in range(L):
        producto = x[: N - m] * x[m:]
        prefijo = np.empty(producto.size + 1, dtype=np.float64)
        prefijo[0] = 0.0
        np.cumsum(producto, out=prefijo[1:])
        i = np.arange(L - m)
        val = prefijo[i + K] - prefijo[i]
        G[i, i + m] = val
        G[i + m, i] = val
    return G


def cuentas_diagonales(L: int, K: int) -> np.ndarray:
    """Numero de terminos de cada antidiagonal i+j=k. Es el peso w_k de la
    w-correlacion y a la vez el divisor de la media diagonal."""
    return np.convolve(np.ones(L), np.ones(K))


def descomponer(x: np.ndarray, L: int, d: int = 30, centrar: bool = True) -> dict:
    """Descomposicion SSA de `x` con longitud de ventana `L`, quedandose con las
    `d` componentes principales.

    Devuelve un diccionario con:
      lambdas      : espectro de autovalores COMPLETO (L valores, descendente)
      U            : autovectores (L x d) -- las EOF empiricas
      V            : componentes principales (K x d)
      elementales  : series elementales reconstruidas (d x N)
      cuentas      : pesos w_k
      media        : la media que se resto (0.0 si centrar=False)
      residuo      : x - suma de las d elementales
    """
    x = np.asarray(x, dtype=np.float64)
    N = x.size
    K = N - L + 1
    media = float(np.mean(x)) if centrar else 0.0
    xc = x - media

    G = matriz_gram(xc, L)
    # eigh devuelve ascendente; se invierte.
    lam, U = np.linalg.eigh(G)
    orden = np.argsort(lam)[::-1]
    lam = lam[orden]
    U = U[:, orden]
    lam = np.clip(lam, 0.0, None)

    d = int(min(d, L))
    Ud = U[:, :d]
    sigma = np.sqrt(lam[:d])

    # V_i = X^T U_i / sigma_i.  (X^T U_i)[j] = sum_k x[j+k] U_i[k]  -> correlacion.
    V = np.empty((K, d), dtype=np.float64)
    elementales = np.empty((d, N), dtype=np.float64)
    cuentas = cuentas_diagonales(L, K)

    for i in range(d):
        proy = sig_sp.correlate(xc, Ud[:, i], mode="valid", method="fft")
        V[:, i] = proy / sigma[i] if sigma[i] > 0 else 0.0
        # media diagonal de sigma_i U_i V_i^T = conv(U_i, X^T U_i) / cuentas
        elementales[i] = sig_sp.fftconvolve(Ud[:, i], proy)[:N] / cuentas

    residuo = xc - elementales.sum(axis=0)

    return {
        "L": L,
        "K": K,
        "N": N,
        "d": d,
        "lambdas": lam,
        "U": Ud,
        "V": V,
        "elementales": elementales,
        "cuentas": cuentas,
        "media": media,
        "residuo": residuo,
    }


# ===========================================================================
# 2. w-correlacion y metrica de ortogonalidad
# ===========================================================================

def matriz_wcorrelacion(elementales: np.ndarray, cuentas: np.ndarray) -> np.ndarray:
    """Matriz de w-correlacion entre series elementales.

    <x,y>_w = sum_k w_k x_k y_k, con w_k = numero de veces que el elemento k
    aparece en la matriz de trayectoria. Dos componentes w-ortogonales son
    separables; una w-correlacion alta significa que la descomposicion partio un
    unico modo fisico en dos trozos (o mezclo dos modos en uno).
    """
    w = cuentas
    prod = elementales * w  # (d x N)
    F = prod @ elementales.T  # <x_i, x_j>_w
    norma = np.sqrt(np.clip(np.diag(F), 1e-300, None))
    W = F / np.outer(norma, norma)
    return np.clip(W, -1.0, 1.0)


def metrica_ortogonalidad(W: np.ndarray, d_eval: int | None = None) -> dict:
    """Reduce la matriz de w-correlacion a escalares comparables entre L distintos.

    Se evalua sobre las primeras `d_eval` componentes: comparar matrices de
    tamano distinto no tendria sentido, y las componentes de cola son ruido cuya
    w-correlacion mutua no informa sobre la separacion de los modos.
    """
    d = W.shape[0] if d_eval is None else int(min(d_eval, W.shape[0]))
    sub = np.abs(W[:d, :d])
    fuera = sub[~np.eye(d, dtype=bool)]
    return {
        "media_abs": float(np.mean(fuera)),
        "max_abs": float(np.max(fuera)),
        "p90_abs": float(np.percentile(fuera, 90)),
        "d_eval": d,
    }


# ===========================================================================
# 3. Barrido de L y eleccion por minimo local
# ===========================================================================

def sustituto_barajado(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Baraja los INCREMENTOS y reintegra. Conserva exactamente la marginal de
    incrementos (curtosis, masa en cero) y destruye todo el orden temporal.
    Es el mismo nulo que la v2.2 uso contra la EMD."""
    dx = np.diff(x)
    return np.concatenate(([x[0]], x[0] + np.cumsum(rng.permutation(dx))))


def ajustar_ar1(x: np.ndarray) -> tuple[float, float, float]:
    """AR(1) sobre la serie: x_t = mu + a (x_{t-1} - mu) + eps. Devuelve (a, mu, s_eps)."""
    x = np.asarray(x, dtype=np.float64)
    mu = float(np.mean(x))
    xc = x - mu
    num = float(np.dot(xc[:-1], xc[1:]))
    den = float(np.dot(xc[:-1], xc[:-1]))
    a = num / den if den > 0 else 0.0
    a = float(np.clip(a, -0.999, 0.999))
    res = xc[1:] - a * xc[:-1]
    return a, mu, float(np.std(res, ddof=1))


def sustituto_ar1_incrementos(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Nulo ARIMA(1,1,0): AR(1) sobre los INCREMENTOS, reintegrado.

    Es el nulo que hace falta aqui y no el AR(1) sobre el nivel. El precio de
    transaccion tiene rho_1 = -0.216 en los retornos (rebote bid-ask, medido en
    la v3.0), y un AR(1) ajustado al NIVEL no puede reproducir eso: sale con
    a ~ 0.998 e innovaciones independientes. Contrastar contra el nulo
    equivocado convierte microestructura en "componente significativa".

    Es la misma leccion de la sesion 2026-08-08 (e): el control solo vale si su
    nulo reproduce las propiedades del dato que importan.
    """
    dx = np.diff(np.asarray(x, dtype=np.float64))
    a, mu, s = ajustar_ar1(dx)
    n = dx.size
    e = np.empty(n, dtype=np.float64)
    var_est = s * s / max(1.0 - a * a, 1e-12)
    e[0] = rng.normal(0.0, np.sqrt(var_est))
    ruido = rng.normal(0.0, s, size=n)
    for t in range(1, n):
        e[t] = a * e[t - 1] + ruido[t]
    return np.concatenate(([x[0]], x[0] + np.cumsum(e + mu)))


def sustituto_ar1(n: int, a: float, mu: float, s: float, rng: np.random.Generator) -> np.ndarray:
    """Realizacion de un AR(1) con los parametros dados, arrancado en estacionario."""
    y = np.empty(n, dtype=np.float64)
    var_est = s * s / max(1.0 - a * a, 1e-12)
    y[0] = rng.normal(0.0, np.sqrt(var_est))
    ruido = rng.normal(0.0, s, size=n)
    for t in range(1, n):
        y[t] = a * y[t - 1] + ruido[t]
    return y + mu


def barrer_L(x: np.ndarray, Ls, d: int = 30, d_eval: int = 10,
             n_sustitutos: int = 3, semilla: int = 0) -> dict:
    """Barre la rejilla `Ls` midiendo la ortogonalidad de la descomposicion.

    Devuelve, por cada L, la metrica sobre los datos y la misma metrica sobre
    `n_sustitutos` series barajadas. Sin esa segunda columna una curva
    descendente no se puede interpretar: podria ser separacion real o podria ser
    que la metrica baja con L pase lo que pase.
    """
    rng = np.random.default_rng(semilla)
    Ls = [int(v) for v in Ls]
    media, maximo, p90 = [], [], []
    media_sus, sd_sus = [], []

    sustitutos = [sustituto_barajado(x, rng) for _ in range(n_sustitutos)]

    for L in Ls:
        des = descomponer(x, L, d=d)
        W = matriz_wcorrelacion(des["elementales"], des["cuentas"])
        m = metrica_ortogonalidad(W, d_eval)
        media.append(m["media_abs"])
        maximo.append(m["max_abs"])
        p90.append(m["p90_abs"])

        vals = []
        for s in sustitutos:
            ds = descomponer(s, L, d=d)
            Ws = matriz_wcorrelacion(ds["elementales"], ds["cuentas"])
            vals.append(metrica_ortogonalidad(Ws, d_eval)["media_abs"])
        media_sus.append(float(np.mean(vals)))
        sd_sus.append(float(np.std(vals)) if len(vals) > 1 else 0.0)

    return {
        "Ls": np.array(Ls),
        "media_abs": np.array(media),
        "max_abs": np.array(maximo),
        "p90_abs": np.array(p90),
        "media_abs_sustituto": np.array(media_sus),
        "sd_sustituto": np.array(sd_sus),
        "d": d,
        "d_eval": d_eval,
        "n_sustitutos": n_sustitutos,
    }


def elegir_minimo_local(Ls, metrica, suavizado: int = 1) -> dict:
    """Elige L por MINIMO LOCAL de la metrica, no por cero absoluto.

    El cero absoluto es inalcanzable -- el ruido de medicion siempre filtra
    energia entre componentes -- asi que una condicion de parada `== 0` no
    termina nunca. Aqui la rejilla es finita, de modo que no hay bucle; lo que
    hay que decidir es cual de los valores medidos es una eleccion y cual es un
    artefacto del borde.

    Reglas:
      - solo cuentan minimos INTERIORES (m[i-1] > m[i] < m[i+1]);
      - si el minimo global cae en un extremo y no hay ninguno interior, se
        devuelve `hay_minimo_local = False` y `L = None`: el barrido no
        discrimina y decir lo contrario seria inventar precision.
    """
    Ls = np.asarray(Ls)
    m = np.asarray(metrica, dtype=np.float64)
    if suavizado > 1:
        nucleo = np.ones(suavizado) / suavizado
        m_s = np.convolve(m, nucleo, mode="same")
    else:
        m_s = m

    interiores = []
    for i in range(1, len(m_s) - 1):
        if m_s[i] <= m_s[i - 1] and m_s[i] <= m_s[i + 1] and (m_s[i] < m_s[i - 1] or m_s[i] < m_s[i + 1]):
            # Prominencia: cuanto sube la curva a ambos lados antes de volver a bajar.
            izq = np.max(m_s[:i]) - m_s[i]
            der = np.max(m_s[i + 1:]) - m_s[i]
            interiores.append((i, float(min(izq, der))))

    if not interiores:
        return {
            "hay_minimo_local": False,
            "L": None,
            "indice": None,
            "valor": None,
            "prominencia": None,
            "argmin_global": int(np.argmin(m_s)),
            "L_argmin_global": int(Ls[int(np.argmin(m_s))]),
            "en_borde": bool(np.argmin(m_s) in (0, len(m_s) - 1)),
            "motivo": "no hay minimo interior en la rejilla",
        }

    # De los minimos interiores se toma el mas profundo, desempatando por prominencia.
    interiores.sort(key=lambda par: (m_s[par[0]], -par[1]))
    i, prom = interiores[0]
    return {
        "hay_minimo_local": True,
        "L": int(Ls[i]),
        "indice": int(i),
        "valor": float(m[i]),
        "prominencia": prom,
        "argmin_global": int(np.argmin(m_s)),
        "L_argmin_global": int(Ls[int(np.argmin(m_s))]),
        "en_borde": bool(np.argmin(m_s) in (0, len(m_s) - 1)),
        "motivo": "minimo local interior",
    }


# ===========================================================================
# 4. Deteccion de pares oscilatorios
# ===========================================================================

def espectro_autovector(u: np.ndarray, relleno: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """Periodograma de un autovector, con relleno de ceros para afinar el argmax.
    El relleno NO anade resolucion (esa sigue siendo 1/L); solo interpola."""
    L = u.size
    n = int(relleno * L)
    esp = np.abs(np.fft.rfft(u - np.mean(u), n=n)) ** 2
    f = np.fft.rfftfreq(n, d=1.0)
    return f, esp


def frecuencia_dominante(u: np.ndarray) -> tuple[float, float]:
    """Frecuencia dominante (ciclos por muestra) y fraccion de energia en ella."""
    f, esp = espectro_autovector(u)
    if esp.size < 2:
        return 0.0, 0.0
    k = int(np.argmax(esp[1:]) + 1)  # se excluye la continua
    total = float(np.sum(esp[1:]))
    # Concentracion: energia en un entorno de +-1 bin de Rayleigh alrededor del pico.
    L = u.size
    ancho = max(1, int(round(esp.size / L)))
    lo, hi = max(1, k - ancho), min(esp.size, k + ancho + 1)
    conc = float(np.sum(esp[lo:hi]) / total) if total > 0 else 0.0
    return float(f[k]), conc


def desfase_mediano(u1: np.ndarray, u2: np.ndarray) -> float:
    """Desfase mediano entre dos autovectores via senal analitica. Un par
    oscilatorio genuino esta en cuadratura: |desfase| ~ pi/2."""
    z1 = sig_sp.hilbert(u1 - np.mean(u1))
    z2 = sig_sp.hilbert(u2 - np.mean(u2))
    dphi = np.angle(z2 * np.conj(z1))
    # Se recorta el 10 % de cada borde: la transformada de Hilbert tiene efecto
    # de borde y ahi el desfase no significa nada.
    m = len(dphi) // 10
    if len(dphi) - 2 * m > 8:
        dphi = dphi[m:-m]
    return float(np.median(np.abs(dphi)))


def amortiguamiento_ar2(serie: np.ndarray) -> dict:
    """Ajusta AR(2) a una componente reconstruida y traduce a ciclo estocastico
    de Harvey: rho = sqrt(-phi2), lambda = arccos(phi1 / (2 sqrt(-phi2))).

    Es la misma correspondencia del §4.1 de la v3.0. Alli fallo aplicada al
    precio ENTERO (phi2 = +0.216 > 0, sin parametrizacion valida). Aplicada a una
    componente SSA ya aislada puede tener solucion o puede no tenerla, y que no
    la tenga es informacion: significa que esa componente no es un ciclo.
    """
    y = np.asarray(serie, dtype=np.float64)
    y = y - np.mean(y)
    if y.size < 20:
        return {"valido": False, "motivo": "serie demasiado corta"}
    Z = np.column_stack([y[1:-1], y[:-2]])
    obj = y[2:]
    coef, *_ = np.linalg.lstsq(Z, obj, rcond=None)
    phi1, phi2 = float(coef[0]), float(coef[1])
    if phi2 >= 0:
        return {"valido": False, "phi1": phi1, "phi2": phi2,
                "motivo": "phi2 >= 0: raices reales, no es un ciclo"}
    rho = float(np.sqrt(-phi2))
    arg = phi1 / (2.0 * np.sqrt(-phi2))
    if abs(arg) > 1.0:
        return {"valido": False, "phi1": phi1, "phi2": phi2, "rho": rho,
                "motivo": "|phi1/(2 sqrt(-phi2))| > 1: raices reales"}
    lam = float(np.arccos(arg))
    return {
        "valido": True,
        "phi1": phi1,
        "phi2": phi2,
        "rho": rho,
        "lambda_rad": lam,
        "periodo_muestras": float(2.0 * np.pi / lam) if lam > 0 else np.inf,
    }


def detectar_pares(des: dict, dt_seg: float = 1.0, d_max: int | None = None) -> list:
    """Recorre las componentes por parejas adyacentes y marca las que forman un
    par oscilatorio. Devuelve una lista de diccionarios, ordenada por energia.

    Se prueban solo parejas ADYACENTES en el espectro porque un modo oscilatorio
    produce dos autovalores casi iguales y por tanto contiguos. Si un par
    genuino aparece separado, la energia esta mal ordenada y eso ya lo delata la
    w-correlacion.
    """
    lam = des["lambdas"]
    U = des["U"]
    elem = des["elementales"]
    d = des["d"] if d_max is None else int(min(d_max, des["d"]))
    total = float(np.sum(lam)) if np.sum(lam) > 0 else 1.0

    pares = []
    i = 0
    while i < d - 1:
        l1, l2 = lam[i], lam[i + 1]
        razon = float(l2 / l1) if l1 > 0 else 0.0
        f1, c1 = frecuencia_dominante(U[:, i])
        f2, c2 = frecuencia_dominante(U[:, i + 1])
        L = des["L"]
        cerca = abs(f1 - f2) <= BINS_RAYLEIGH_PAR / L
        conc = 0.5 * (c1 + c2)
        dphi = desfase_mediano(U[:, i], U[:, i + 1])
        cuadratura = abs(dphi - np.pi / 2.0) <= TOLERANCIA_CUADRATURA

        es_par = (razon >= RAZON_LAMBDA_PAR and cerca and
                  conc >= CONCENTRACION_PAR and cuadratura)

        if es_par:
            f = 0.5 * (f1 + f2)
            serie = elem[i] + elem[i + 1]
            ar2 = amortiguamiento_ar2(serie)
            pares.append({
                "indices": (i, i + 1),
                "lambda1": float(l1),
                "lambda2": float(l2),
                "razon_lambda": razon,
                "energia_frac": float((l1 + l2) / total),
                "f_ciclos_por_muestra": float(f),
                "periodo_muestras": float(1.0 / f) if f > 0 else np.inf,
                "periodo_seg": float(dt_seg / f) if f > 0 else np.inf,
                "concentracion": float(conc),
                "desfase_rad": float(dphi),
                "ar2": ar2,
            })
            i += 2
        else:
            i += 1

    pares.sort(key=lambda p: -p["energia_frac"])
    return pares


def escalera_de_ventana(des: dict, k_max: int = 12) -> dict:
    """¿Los autovectores son modos del MERCADO o armonicos de la VENTANA?

    Sobre un paseo aleatorio, la base que produce el SSA no depende de la serie:
    son sinusoides con periodos T_k = 2L/k, la escalera armonica de la propia
    ventana. Sobre una senal con ciclos, los primeros autovectores llevan el
    periodo de la senal y NO caen en esa escalera.

    Devuelve el error relativo mediano entre el periodo de cada autovector y el
    escalon 2L/k mas cercano. Cerca de cero = la descomposicion esta describiendo
    la ventana; lejos = hay algo mas.

    Es el equivalente en SSA del test que la v2.2 uso contra la EMD (el periodo
    escalaba con la ventana), pero con la ventaja de que aqui no hace falta
    partir la ventana en dos: la escalera se ve entera de una vez.
    """
    U = des["U"]
    L = des["L"]
    lam = des["lambdas"]
    d = int(min(k_max, U.shape[1]))
    escalones = np.array([2.0 * L / k for k in range(1, 4 * d + 2)])

    periodos, errores = [], []
    for i in range(d):
        f, _ = frecuencia_dominante(U[:, i])
        T = (1.0 / f) if f > 0 else np.inf
        periodos.append(T)
        if np.isfinite(T):
            errores.append(float(np.min(np.abs(escalones - T)) / T))
        else:
            errores.append(np.nan)

    errores = np.array(errores, dtype=np.float64)
    fin = errores[np.isfinite(errores)]
    total = float(np.sum(lam)) if np.sum(lam) > 0 else 1.0
    return {
        "L": L,
        "periodos": np.array(periodos),
        "T_sobre_L": np.array(periodos) / L,
        "error_a_escalera": errores,
        "error_mediano": float(np.median(fin)) if fin.size else float("nan"),
        "energia": lam[:d] / total,
        # Energia fuera del primer autovector: sobre un paseo aleatorio casi todo
        # se va al primero y el resto es la escalera con energia despreciable.
        "energia_fuera_del_primero": float(np.sum(lam[1:d]) / total),
    }


# ===========================================================================
# 5. Monte Carlo SSA -- la hipotesis nula que la EMD no tenia
# ===========================================================================

def monte_carlo_ssa(x: np.ndarray, L: int, d: int = 30, n_sustitutos: int = 100,
                    nulo: str = "ar1", semilla: int = 0, centrar: bool = True) -> dict:
    """Monte Carlo SSA (Allen & Smith 1996).

    Se ajusta un AR(1) a los datos, se generan `n_sustitutos` realizaciones y se
    proyecta la matriz de Gram de cada una sobre los AUTOVECTORES DE LOS DATOS:
        lambda_i^surr = U_i^T G_surr U_i
    Un autovalor observado por encima del percentil 95 de esa distribucion es una
    componente que el ruido rojo no explica.

    ADVERTENCIA (sesgo conocido del metodo): proyectar sobre las EOF de los
    propios datos favorece a los datos, porque esas direcciones se eligieron para
    maximizar su varianza. El test es por tanto ANTICONSERVADOR y su p-valor no
    es exacto. Se reporta ademas cuantas componentes superan el umbral frente a
    las ~5 % esperadas por azar, que es la lectura que la v2.2 aplico al
    multitaper y la unica defendible aqui.
    """
    x = np.asarray(x, dtype=np.float64)
    N = x.size
    des = descomponer(x, L, d=d, centrar=centrar)
    U = des["U"]
    lam_obs = des["lambdas"][:des["d"]]

    rng = np.random.default_rng(semilla)
    a, mu, s = ajustar_ar1(x)

    lam_nulo = np.empty((n_sustitutos, des["d"]), dtype=np.float64)
    for j in range(n_sustitutos):
        if nulo == "ar1":
            y = sustituto_ar1(N, a, mu, s, rng)
        elif nulo == "ar1_incrementos":
            y = sustituto_ar1_incrementos(x, rng)
        elif nulo == "barajado":
            y = sustituto_barajado(x, rng)
        else:
            raise ValueError("nulo desconocido: %s" % nulo)
        if centrar:
            y = y - np.mean(y)
        Gy = matriz_gram(y, L)
        # lambda_i = U_i^T G U_i, sin volver a diagonalizar.
        lam_nulo[j] = np.einsum("ij,jk,ki->i", U.T, Gy, U)

    p5 = np.percentile(lam_nulo, 5, axis=0)
    p95 = np.percentile(lam_nulo, 95, axis=0)
    sobre = lam_obs > p95
    # p-valor empirico por componente (fraccion de sustitutos que igualan o superan).
    pvals = np.array([float(np.mean(lam_nulo[:, i] >= lam_obs[i])) for i in range(des["d"])])

    return {
        "L": L,
        "d": des["d"],
        "nulo": nulo,
        "ar1_a": a,
        "ar1_sigma": s,
        "lambda_obs": lam_obs,
        "p5": p5,
        "p95": p95,
        "sobre_p95": sobre,
        "n_sobre": int(np.sum(sobre)),
        "frac_sobre": float(np.mean(sobre)),
        "esperado_por_azar": 0.05,
        "pvalor": pvals,
        "n_sustitutos": n_sustitutos,
    }


# ===========================================================================
# 6. Color de ruido
# ===========================================================================

def color_de_ruido(x: np.ndarray, dt_seg: float = 1.0, nperseg: int | None = None,
                   banda: tuple[float, float] = (0.02, 0.40)) -> dict:
    """Estima beta en S(f) ~ f^(-beta) por regresion en log-log y clasifica.

    `banda` esta en fraccion de Nyquist: se excluyen las frecuencias mas bajas
    (donde Welch tiene pocos grados de libertad y el sesgo de tendencia domina) y
    la region de Nyquist (donde el aliasing y el suelo de cuantizacion mandan).

    NOTA DE INTERPRETACION: las fronteras de BANDAS_COLOR son las convencionales
    (blanco 0, rosa 1, browniano 2). Si R2 < R2_MINIMO_LEY_POTENCIAS la etiqueta
    se marca como no valida: un espectro que no es una recta en log-log no tiene
    un "color", y ponerle uno seria una lectura inventada.
    """
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    if nperseg is None:
        nperseg = int(min(2048, max(64, n // 8)))
    f, S = sig_sp.welch(x - np.mean(x), fs=1.0 / dt_seg, nperseg=nperseg,
                        detrend="constant")
    fny = 0.5 / dt_seg
    mask = (f >= banda[0] * fny) & (f <= banda[1] * fny) & (S > 0)
    if np.count_nonzero(mask) < 8:
        return {"valido": False, "motivo": "muy pocas frecuencias en la banda"}

    lf = np.log(f[mask])
    lS = np.log(S[mask])
    A = np.column_stack([np.ones_like(lf), lf])
    coef, *_ = np.linalg.lstsq(A, lS, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((lS - pred) ** 2))
    ss_tot = float(np.sum((lS - np.mean(lS)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    beta = -float(coef[1])

    # Error estandar de la pendiente (OLS, homocedastico -- Welch aproxima esto
    # razonablemente en log-log porque el log de una chi2 tiene varianza constante).
    gl = max(1, len(lf) - 2)
    s2 = ss_res / gl
    cov = s2 * np.linalg.inv(A.T @ A)
    se_beta = float(np.sqrt(cov[1, 1]))

    etiqueta = "INDETERMINADO"
    for lo, hi, nombre in BANDAS_COLOR:
        if lo <= beta < hi:
            etiqueta = nombre
            break
    valido = r2 >= R2_MINIMO_LEY_POTENCIAS

    a_ar1, _, _ = ajustar_ar1(x)
    return {
        "valido": bool(valido),
        "beta": beta,
        "se_beta": se_beta,
        "r2": float(r2),
        "color": etiqueta if valido else "NO ES LEY DE POTENCIAS",
        "ar1_a": float(a_ar1),
        "f": f,
        "S": S,
        "banda_fny": banda,
        "n_frec_ajuste": int(np.count_nonzero(mask)),
    }


# ===========================================================================
# 7. De los pares a la matriz A del EAKF
# ===========================================================================

def bloque_ciclo(omega_rad_por_paso: float, rho: float = 1.0) -> np.ndarray:
    """Bloque 2x2 de rotacion amortiguada: un ciclo estocastico de Harvey.

        [ s ]      [ cos w  -sin w ] [ s ]
        [ s*] <- rho [ sin w   cos w ] [ s*]

    Con rho = 1 es un oscilador sin amortiguar; con rho < 1 el ciclo se apaga si
    no lo excita el ruido de proceso, que es el comportamiento que la Sec. 1 del
    PDF describe cualitativamente y que `A_arm` no podia representar.
    """
    w = float(omega_rad_por_paso)
    c, s = np.cos(w), np.sin(w)
    return float(rho) * np.array([[c, -s], [s, c]], dtype=np.float64)


def matriz_A_estructural(pares: list, con_tendencia: bool = True,
                         dt_paso: float = 1.0, usar_ar2: bool = True) -> dict:
    """Ensambla A en bloques a partir de los pares oscilatorios detectados.

    Estado: [nivel, (velocidad), s_1, s*_1, s_2, s*_2, ...]
    Observacion: H selecciona el nivel y la primera coordenada de cada ciclo.

    Esto es lo que sustituye a `A_arm`: en vez de asumir UNA frecuencia unica
    salida de un colapso espectral, se montan tantos bloques como ciclos haya
    sobrevivido al test nulo. Si no sobrevive ninguno, A queda reducida a la
    tendencia -- y eso es un resultado, no un fallo.
    """
    bloques = []
    nombres = []
    if con_tendencia:
        bloques.append(np.array([[1.0, dt_paso], [0.0, 1.0]], dtype=np.float64))
        nombres += ["nivel", "velocidad"]

    detalle = []
    for k, p in enumerate(pares):
        # omega en rad por paso de muestreo del SSA.
        w = 2.0 * np.pi * p["f_ciclos_por_muestra"]
        rho = 1.0
        fuente = "sin amortiguamiento (rho = 1)"
        if usar_ar2 and p["ar2"].get("valido"):
            rho = float(p["ar2"]["rho"])
            # Comprobacion cruzada: la omega del AR(2) contra la del autovector.
            w_ar2 = float(p["ar2"]["lambda_rad"])
            fuente = "AR(2) de la componente (w_eig=%.5f, w_ar2=%.5f)" % (w, w_ar2)
        bloques.append(bloque_ciclo(w, rho))
        nombres += ["ciclo%d" % k, "ciclo%d*" % k]
        detalle.append({"k": k, "omega_rad_por_paso": w, "rho": rho,
                        "periodo_muestras": p["periodo_muestras"],
                        "periodo_seg": p["periodo_seg"], "fuente_rho": fuente})

    if not bloques:
        A = np.eye(1)
        nombres = ["nivel"]
        H = np.array([[1.0]])
    else:
        A = lin_sp.block_diag(*bloques)
        H = np.zeros((1, A.shape[0]))
        idx = 0
        if con_tendencia:
            H[0, 0] = 1.0
            idx = 2
        for _ in pares:
            H[0, idx] = 1.0
            idx += 2

    return {"A": A, "H": H, "estado": nombres, "ciclos": detalle,
            "n_ciclos": len(pares), "con_tendencia": con_tendencia}


# ===========================================================================
# 8. Autotest -- controles POSITIVOS: dada una verdad conocida, ¿la recupera?
# ===========================================================================

def _autotest() -> int:
    """Controles positivos y negativos del propio estimador.

    Este proyecto lleva cinco sesiones ejecutando controles NEGATIVOS con
    disciplina (barajados, IAAFT, AR(1)). La mitad que faltaba son los
    positivos: con una verdad conocida, ¿el estimador la recupera? Eso es lo que
    hay aqui, y es la razon de que este fichero traiga su propio `--autotest`.
    """
    fallos = 0
    rng = np.random.default_rng(7)

    def ok(nombre, cond, detalle=""):
        nonlocal fallos
        estado = "OK  " if cond else "FALLA"
        if not cond:
            fallos += 1
        print("  [%s] %s %s" % (estado, nombre, detalle))

    print("== 1. La descomposicion es EXACTA (suma de elementales = serie) ==")
    x = np.cumsum(rng.normal(0, 1, 2000))
    des = descomponer(x, L=200, d=200)
    err = float(np.max(np.abs(des["elementales"].sum(axis=0) + des["media"] - x)))
    ok("reconstruccion exacta", err < 1e-8, "error_max = %.3e" % err)

    print("== 2. CONTROL POSITIVO: dos periodos conocidos (120 y 37 muestras) ==")
    n = 4000
    t = np.arange(n)
    verdad = [120.0, 37.0]
    y = (3.0 * np.sin(2 * np.pi * t / 120.0 + 0.4)
         + 1.5 * np.sin(2 * np.pi * t / 37.0)
         + 0.02 * t
         + rng.normal(0, 0.5, n))
    des2 = descomponer(y, L=400, d=20)
    pares = detectar_pares(des2, dt_seg=1.0)
    periodos = sorted([p["periodo_muestras"] for p in pares], reverse=True)
    print("     periodos recuperados: %s" % ["%.2f" % v for v in periodos[:4]])
    hallados = []
    for v in verdad:
        cerca = [q for q in periodos if abs(q - v) / v < 0.05]
        hallados.append(len(cerca) > 0)
        if cerca:
            print("     verdad %.1f -> %.2f (error %.2f %%)"
                  % (v, cerca[0], 100 * abs(cerca[0] - v) / v))
        else:
            print("     verdad %.1f -> NO RECUPERADO" % v)
    ok("recupera ambos periodos al 5 %", all(hallados))

    print("== 3. La w-correlacion detecta la mezcla ==")
    W2 = matriz_wcorrelacion(des2["elementales"], des2["cuentas"])
    m_bien = metrica_ortogonalidad(W2, 6)["media_abs"]
    des_mal = descomponer(y, L=25, d=20)   # L << periodo largo: no puede separar
    W_mal = matriz_wcorrelacion(des_mal["elementales"], des_mal["cuentas"])
    m_mal = metrica_ortogonalidad(W_mal, 6)["media_abs"]
    ok("L adecuada separa mejor que L corta", m_bien < m_mal,
       "media|wcorr| = %.4f (L=400) contra %.4f (L=25)" % (m_bien, m_mal))

    print("== 4. Color de ruido: blanco, rojo y AR(1) ==")
    blanco = rng.normal(0, 1, 20000)
    c_b = color_de_ruido(blanco)
    ok("ruido blanco -> beta ~ 0", abs(c_b["beta"]) < 0.25,
       "beta = %+.3f +- %.3f, R2 = %.3f, etiqueta = %s"
       % (c_b["beta"], c_b["se_beta"], c_b["r2"], c_b["color"]))
    rojo = np.cumsum(rng.normal(0, 1, 20000))
    c_r = color_de_ruido(rojo)
    ok("paseo aleatorio -> beta ~ 2", abs(c_r["beta"] - 2.0) < 0.35,
       "beta = %+.3f +- %.3f, R2 = %.3f, etiqueta = %s"
       % (c_r["beta"], c_r["se_beta"], c_r["r2"], c_r["color"]))

    print("== 5. CONTROL NEGATIVO: Monte Carlo SSA sobre AR(1) puro ==")
    a_ver = 0.6
    ruido_ar = sustituto_ar1(6000, a_ver, 0.0, 1.0, rng)
    mc = monte_carlo_ssa(ruido_ar, L=200, d=20, n_sustitutos=40, semilla=3)
    print("     a estimada = %.3f (verdad %.2f); componentes sobre p95: %d/%d = %.1f %%"
          % (mc["ar1_a"], a_ver, mc["n_sobre"], mc["d"], 100 * mc["frac_sobre"]))
    # Anticonservador por construccion, pero no puede marcar casi todo.
    ok("no marca mas de la mitad de las componentes", mc["frac_sobre"] <= 0.50)

    print("== 6. CONTROL POSITIVO: MC-SSA detecta una senal enterrada en AR(1) ==")
    n6 = 6000
    t6 = np.arange(n6)
    senal = 1.2 * np.sin(2 * np.pi * t6 / 80.0) + sustituto_ar1(n6, a_ver, 0.0, 1.0, rng)
    mc6 = monte_carlo_ssa(senal, L=200, d=20, n_sustitutos=40, semilla=4)
    print("     componentes sobre p95: %d/%d; pvalor de la primera = %.3f"
          % (mc6["n_sobre"], mc6["d"], mc6["pvalor"][0]))
    ok("marca la componente dominante", mc6["sobre_p95"][0])

    print("== 7. El minimo local NO devuelve el borde ni busca el cero ==")
    Ls = np.array([10, 20, 30, 40, 50, 60, 70])
    curva_mono = np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3])  # sin minimo interior
    r_mono = elegir_minimo_local(Ls, curva_mono)
    ok("curva monotona -> no discrimina", not r_mono["hay_minimo_local"],
       "motivo: %s" % r_mono["motivo"])
    curva_v = np.array([0.9, 0.7, 0.4, 0.2, 0.45, 0.6, 0.8])
    r_v = elegir_minimo_local(Ls, curva_v)
    ok("curva en V -> elige el fondo", r_v["hay_minimo_local"] and r_v["L"] == 40,
       "L elegida = %s, prominencia = %.3f" % (r_v["L"], r_v["prominencia"] or 0.0))

    print("== 8. CONTROL NEGATIVO: paseo aleatorio puro -- lo que hundio a la EMD ==")
    pa = np.cumsum(rng.normal(0, 1, 4000))
    for Lp in (400, 200):
        dp = descomponer(pa, L=Lp, d=20)
        pp = detectar_pares(dp, dt_seg=1.0)
        if pp:
            print("     L=%d -> %d par(es); periodo dominante %.1f muestras (%.2f*L)"
                  % (Lp, len(pp), pp[0]["periodo_muestras"], pp[0]["periodo_muestras"] / Lp))
        else:
            print("     L=%d -> ningun par oscilatorio" % Lp)
    mc8 = monte_carlo_ssa(pa, L=400, d=20, n_sustitutos=40, semilla=5)
    print("     MC-SSA contra AR(1): %d/%d componentes sobre p95 (%.1f %%)"
          % (mc8["n_sobre"], mc8["d"], 100 * mc8["frac_sobre"]))
    print("     (informativo: si el periodo escala con L, SSA hereda el defecto de la EMD)")

    print("== 9. La ESCALERA DE VENTANA separa senal de armonicos de la ventana ==")
    # Sobre un paseo aleatorio los autovectores son sinusoides de periodo 2L/k:
    # la base es de la VENTANA, no de la serie. Sobre una senal con ciclos, los
    # primeros autovectores llevan el periodo de la senal y rompen la escalera.
    # Medido sobre mercado real (v31b y captura_larga, nu de 6 a 109 tx/s) el
    # error a la escalera es 0.0000-0.0335 y el primer EOF se lleva 84-99 % de
    # la energia: indistinguible del paseo aleatorio.
    esc_pa = escalera_de_ventana(descomponer(pa, L=400, d=20))
    esc_se = escalera_de_ventana(des2)
    print("     paseo aleatorio : error a 2L/k = %.4f | E(EOF1) = %.1f %%"
          % (esc_pa["error_mediano"], 100 * esc_pa["energia"][0]))
    print("     senal con ciclos: error a 2L/k = %.4f | E(EOF1) = %.1f %%"
          % (esc_se["error_mediano"], 100 * esc_se["energia"][0]))
    ok("el paseo aleatorio CAE en la escalera", esc_pa["error_mediano"] < 0.06)
    ok("la senal con ciclos la ROMPE", esc_se["error_mediano"] > 0.10)

    print("== 10. La matriz A se ensambla y es estable ==")
    A_est = matriz_A_estructural(pares, con_tendencia=True)
    A = A_est["A"]
    radios = np.abs(np.linalg.eigvals(A))
    print("     estado: %s" % A_est["estado"])
    print("     radio espectral = %.6f" % float(np.max(radios)))
    ok("A tiene la dimension esperada", A.shape[0] == 2 + 2 * len(pares),
       "shape = %s con %d ciclos" % (A.shape, len(pares)))

    print("")
    print("RESULTADO: %d fallo(s)" % fallos)
    return fallos


if __name__ == "__main__":
    import sys
    if "--autotest" in sys.argv:
        raise SystemExit(1 if _autotest() else 0)
    print("uso: python ssa.py --autotest")
