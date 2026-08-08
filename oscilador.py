"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: oscilador.py — Primitivas del oscilador forzado (ORDEN_TRABAJO_OSCILADOR_3_0)

⚠ ESTE MODULO NO SE IMPORTA DESDE `Micelio.py`. El §8 de la v3.0 lo exige:
"Micelio.py sin cambios de modelo. Todo en código de análisis aparte, como en la
v2.2". Aquí viven las primitivas medibles —`k`, `m`, `γ`, `Q`— y su test con
hipótesis nula.

EL CAMBIO DE PLANTEAMIENTO
--------------------------
Hasta la v2.2 se estimaba `ω_m` como cantidad PRIMARIA, extrayendo el modo
dominante de una ventana. Ahí vive el artefacto: la v2.2 demostró que sobre un
paseo aleatorio la EMD devuelve un "ciclo" cuyo período escala con la ventana
(118.1 s con W=384, 37.5 s con W=192).

`ω₀ = √(k/m)` es DERIVADA. Este módulo mide `k`, `m` y `γ` por regresión sobre
todos los datos, sin ventana espectral, sin banda de resolubilidad, y —lo que
importa— **con un nulo que es exactamente el paseo aleatorio**.
"""

from __future__ import annotations

import math

import numpy as np


# ==============================================================================
# §2.3 — DERIVACION DIMENSIONAL (obligatoria antes de implementar nada)
# ==============================================================================
# PASO 1 — Unidades en tiempo de ticks, con la nomenclatura del §2.4 de la v2.0.
#
#   x       = S − S_ref            [USD/BTC]
#   τ                              [Ticks]  (contador, adimensional en magnitud)
#   ẋ  = dx/dτ                     [USD/BTC/Tick]
#   ẍ  = d²x/dτ²                   [USD/BTC/Tick²]
#
# PASO 2 — Unidades REALES de λ, EXTRAIDAS DEL PDF (Sec. 4.4.1, línea 503 de la
# transcripción): "las unidades obligatorias de la falta de liquidez son
# [λ] = 1/BTC. Representa el impacto en el mercado por cada unidad de
# criptoactivo transada." Se verifica con la propia condición del PDF:
#
#   [λ·S²·Γ] = (1/BTC)·(USD/BTC)²·(BTC³/USD²) = 1   ✓ adimensional
#
# PASO 3 — Los cuatro términos de `m ẍ + γ ẋ + k x = F` con unidades idénticas.
# Llamando [k] a las unidades de la constante recuperadora:
#
#   [k·x]  = [k]·USD/BTC
#   [γ·ẋ]  = [γ]·USD/BTC/Tick     ⟹  [γ] = [k]·Ticks
#   [m·ẍ]  = [m]·USD/BTC/Tick²    ⟹  [m] = [k]·Ticks²
#   [F]    = [k]·USD/BTC
#
# PASO 4 — Q adimensional:
#
#   [Q] = [√(k·m)]/[γ] = √([k]·[k]·Ticks²) / ([k]·Ticks) = 1   ✓
#
# ⚠ Y ESTO ES LO QUE HACE UTIL AL TEST: `Q` sale adimensional CUALQUIERA que sea
# [k]. Por eso el §3.1 puede identificar `Q` sin conocer la escala común
# `D = m + γ + k`, que el forzamiento fija y el AR(2) no puede separar.
#
# ⚠ SOBRE `m ∝ 1/λ` (§2.2). [1/λ] = BTC, pero [m] = [k]·Ticks². La proporción
# NO es directa: exige un factor de conversión de unidades [k]·Ticks²/BTC. La
# intuición física es buena —libro profundo (λ pequeña) cuesta más de mover, o
# sea más inercia— pero la identificación **necesita ese factor y no puede
# escribirse como igualdad**. Se deja anotado en vez de colarlo.
UNIDADES = {
    "x": "USD/BTC",
    "xdot": "USD/BTC/Tick",
    "xddot": "USD/BTC/Tick^2",
    "lambda": "1/BTC",  # Sec. 4.4.1 del PDF, verificado con λS^2Γ adimensional
    "k": "[k]",  # escala libre: el AR(2) solo identifica razones
    "gamma": "[k]*Tick",
    "m": "[k]*Tick^2",
    "F": "[k]*USD/BTC",
    "Q": "1",  # adimensional para cualquier [k]
}


def verificar_dimensiones() -> dict:
    """PASO 5 del §2.3: test automático que falla si algo del anterior se rompe.

    Las trampas del 2π y del factor 125 fueron errores dimensionales
    silenciosos; aquí hay cuatro identificaciones nuevas a la vez.
    """
    # Se representan las unidades como exponentes de (USD, BTC, Tick).
    def u(usd=0, btc=0, tick=0):
        return np.array([usd, btc, tick], dtype=int)

    x = u(1, -1, 0)  # USD/BTC
    xdot = x + u(0, 0, -1)
    xddot = x + u(0, 0, -2)
    k = u(0, 0, 0)  # escala libre; se toma adimensional como referencia
    gamma = k + u(0, 0, 1)
    m = k + u(0, 0, 2)

    t_k = k + x
    t_gamma = gamma + xdot
    t_m = m + xddot
    assert np.array_equal(t_k, t_gamma), (t_k, t_gamma)
    assert np.array_equal(t_k, t_m), (t_k, t_m)

    # Q = sqrt(k*m)/gamma
    q = (k + m) / 2.0 - gamma
    assert np.allclose(q, 0.0), q

    # λS²Γ adimensional (la condición que el PDF usa para fijar [λ])
    lam = u(0, -1, 0)  # 1/BTC
    S2 = 2 * u(1, -1, 0)
    Gamma = u(-2, 3, 0)  # BTC³/USD²
    assert np.allclose(lam + S2 + Gamma, 0.0), lam + S2 + Gamma

    return {
        "termino_comun": tuple(t_k),
        "Q_adimensional": True,
        "lambda_S2_Gamma_adimensional": True,
    }


# ==============================================================================
# §3 — EL OSCILADOR DISCRETIZADO ES UN AR(2)
# ==============================================================================
def primitivas_desde_phi(phi1: float, phi2: float) -> dict:
    """`m`, `γ`, `k` y `Q` a partir de los coeficientes AR(2). §3.1.

        x_t = φ₁ x_{t-1} + φ₂ x_{t-2} + ε_t

    y, salvo la escala común `D = m + γ + k` que el forzamiento fija:

        m ∝ −φ₂        γ ∝ φ₁ + 2φ₂        k ∝ 1 − φ₁ − φ₂

    `Q = √(k·m)/γ` es identificable SIN conocer `D`, que es lo que lo hace útil.
    Condición de oscilación: raíces complejas, `φ₁² + 4φ₂ < 0`, equivalente a
    `Q > ½`.
    """
    m = -phi2
    gamma = phi1 + 2.0 * phi2
    k = 1.0 - phi1 - phi2
    disc = phi1 * phi1 + 4.0 * phi2
    km = k * m
    if km > 0.0 and gamma != 0.0:
        Q = math.sqrt(km) / gamma
    else:
        Q = float("nan")
    return {
        "m": m, "gamma": gamma, "k": k, "Q": Q,
        "discriminante": disc,
        "raices_complejas": bool(disc < 0.0),
    }


def ajustar_ar2(x: np.ndarray, con_intercepto: bool = True):
    """AR(2) por minimos cuadrados. Devuelve (phi1, phi2, residuos, X, y)."""
    x = np.asarray(x, dtype=float)
    y = x[2:]
    X = np.column_stack([x[1:-1], x[:-2]])
    if con_intercepto:
        X = np.column_stack([X, np.ones(len(y))])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    return float(beta[0]), float(beta[1]), resid, X, y


def ajustar_ar2_huber(x: np.ndarray, c: float = 1.345, n_iter: int = 50):
    """AR(2) robusto por IRLS con función de peso de Huber. §3.3.

    Con curtosis de incrementos ~1180 (medida tras limpiar los ceros del feed),
    mínimos cuadrados es ineficiente y sus errores estándar son optimistas. Si
    MCO y Huber difieren de forma material, **prevalece el robusto**.
    """
    x = np.asarray(x, dtype=float)
    y = x[2:]
    X = np.column_stack([x[1:-1], x[:-2], np.ones(len(y))])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    for _ in range(n_iter):
        r = y - X @ beta
        s = 1.4826 * np.median(np.abs(r - np.median(r)))  # MAD -> sigma robusta
        if s <= 0:
            break
        u = r / s
        w = np.where(np.abs(u) <= c, 1.0, c / np.maximum(np.abs(u), 1e-12))
        W = np.sqrt(w)
        beta_new, *_ = np.linalg.lstsq(X * W[:, None], y * W, rcond=None)
        if np.max(np.abs(beta_new - beta)) < 1e-12:
            beta = beta_new
            break
        beta = beta_new
    return float(beta[0]), float(beta[1])


# ==============================================================================
# §3.3 — EL CONTRASTE H0: k = 0
# ==============================================================================
# ⚠ PUNTO ESTADISTICO QUE LA ORDEN DE TRABAJO NO MENCIONA Y QUE DECIDE LA
# VALIDEZ DEL RESULTADO.
#
# `k = 1 − φ₁ − φ₂`, luego `H₀: k = 0` es exactamente `φ₁ + φ₂ = 1`, o sea
# **una raíz unitaria**. Bajo esa hipótesis la serie NO es estacionaria y el
# estadístico t de `k` **no sigue una distribución normal ni t de Student**:
# sigue la distribución de Dickey-Fuller, que está desplazada hacia la izquierda
# y tiene colas mucho más anchas. Un p-valor normal daría significancia
# espuria de forma sistemática.
#
# Es el mismo tipo de error que el §4.3 de la propia orden advierte para el
# contraste de razón de verosimilitudes en la frontera. Aquí se evita del mismo
# modo: **se construye el nulo por simulación** en vez de asumir su forma.
#
# El nulo se genera con el sustituto BARAJADO de la v2.2 —permutar los
# incrementos observados— que conserva EXACTAMENTE la distribución marginal
# (curtosis, retícula de tickSize, masa en cero) y destruye solo el orden
# temporal. Es decir: el nulo es un paseo aleatorio con los MISMOS incrementos
# que los datos reales, que es la hipótesis que la v2.2 no pudo rechazar.
def nulo_por_barajado(x: np.ndarray, n_sim: int, rng) -> np.ndarray:
    """Distribución de `k̂` bajo H₀, por permutación de incrementos."""
    d = np.diff(np.asarray(x, dtype=float))
    ks = np.empty(n_sim)
    for i in range(n_sim):
        sim = np.concatenate([[x[0]], x[0] + np.cumsum(rng.permutation(d))])
        p1, p2, *_ = ajustar_ar2(sim)
        ks[i] = 1.0 - p1 - p2
    return ks


def bootstrap_bloques(x: np.ndarray, n_boot: int, rng, largo_bloque: int = None):
    """Errores estándar de φ₁, φ₂, k y Q por bootstrap POR BLOQUES. §3.3.

    Por bloques y no i.i.d. porque el bootstrap ingenuo destruiría la dependencia
    serial, que es justo lo que el AR(2) está estimando: daría errores estándar
    demasiado optimistas por la misma razón que corregir 90 veces con la misma
    medición contraía P de más.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    if largo_bloque is None:
        largo_bloque = max(50, int(n ** (1.0 / 3.0)))
    n_bloques = max(1, n // largo_bloque)
    out = {"phi1": [], "phi2": [], "k": [], "Q": []}
    for _ in range(n_boot):
        inicios = rng.integers(0, max(1, n - largo_bloque), size=n_bloques)
        sim = np.concatenate([x[i : i + largo_bloque] for i in inicios])
        p1, p2, *_ = ajustar_ar2(sim)
        pr = primitivas_desde_phi(p1, p2)
        out["phi1"].append(p1)
        out["phi2"].append(p2)
        out["k"].append(pr["k"])
        out["Q"].append(pr["Q"])
    return {kk: np.asarray(v, dtype=float) for kk, v in out.items()}


def analizar_serie(x: np.ndarray, rng, n_sim_nulo: int = 400, n_boot: int = 400,
                   etiqueta: str = ""):
    """Ajuste completo del §3 sobre una serie ya muestreada a Δn = 1 tick."""
    x = np.asarray(x, dtype=float)
    p1, p2, resid, _, _ = ajustar_ar2(x)
    pr = primitivas_desde_phi(p1, p2)

    h1, h2 = ajustar_ar2_huber(x)
    pr_h = primitivas_desde_phi(h1, h2)

    ks_nulo = nulo_por_barajado(x, n_sim_nulo, rng)
    # p-valor de una cola: ¿con qué frecuencia el NULO produce una k tan grande
    # como la observada? Una cola porque k < 0 no tiene sentido físico (sería
    # una fuerza que aleja del equilibrio, o sea explosiva).
    p_valor = float(np.mean(ks_nulo >= pr["k"]))

    boot = bootstrap_bloques(x, n_boot, rng)

    return {
        "etiqueta": etiqueta,
        "n": len(x),
        "phi1": p1, "phi2": p2,
        "mco": pr,
        "huber": pr_h,
        "phi1_huber": h1, "phi2_huber": h2,
        "p_valor_k": p_valor,
        "k_nulo_p50": float(np.median(ks_nulo)),
        "k_nulo_p95": float(np.percentile(ks_nulo, 95)),
        "se_k": float(np.std(boot["k"], ddof=1)),
        "se_phi1": float(np.std(boot["phi1"], ddof=1)),
        "se_phi2": float(np.std(boot["phi2"], ddof=1)),
        "Q_boot_p05": float(np.nanpercentile(boot["Q"], 5)),
        "Q_boot_p95": float(np.nanpercentile(boot["Q"], 95)),
        "resid_std": float(np.std(resid)),
    }


def guarda_masa(res: dict) -> str | None:
    """`m < 0` es masa negativa: la identificacion del oscilador NO se sostiene.

    ⚠ LA ORDEN DE TRABAJO PIDE UNA GUARDA PARA `γ < 0` PERO NO PARA `m < 0`, Y EL
    MODO DE FALLO REAL FUE ESTE. Con `m = −φ₂`, una masa negativa significa
    `φ₂ > 0`, y eso tiene una lectura mecanica exacta:

        si  x_t = x_{t-1} + r_t   con   r_t = a·r_{t-1} + ε   (retornos AR(1))
        entonces  x_t = (1+a)·x_{t-1} − a·x_{t-2} + ε
        o sea     φ₁ = 1+a ,  φ₂ = −a ,  y  φ₁+φ₂ = 1 EXACTAMENTE

    Es decir: un paseo aleatorio cuyos RETORNOS estan autocorrelacionados produce
    un AR(2) con `k = 0` por construccion y `m = a`. Con rebote bid-ask (`a < 0`,
    medido −0.164 en la v2.0) sale `m < 0`.

    Masa negativa no es "poca inercia": es la firma de que el segundo rezago
    viene de la microestructura de los retornos y no de una inercia del precio.
    """
    m = res["mco"]["m"]
    if m < 0.0:
        return (
            f"m = {m:.6f} < 0 en [{res['etiqueta']}]: MASA NEGATIVA. El segundo "
            f"rezago viene de la autocorrelacion de RETORNOS (rebote bid-ask), no "
            f"de inercia. La identificacion del oscilador no se sostiene y Q no "
            f"esta definida."
        )
    return None


def descomponer_rebote(x: np.ndarray) -> dict:
    """Contrasta el AR(2) medido contra la prediccion del rebote bid-ask.

    Si `φ₁ ≈ 1+a` y `φ₂ ≈ −a` con `a` la autocorrelacion de retornos, entonces
    TODO el AR(2) es paseo aleatorio + microestructura, y `k = 0` no es un
    resultado marginal sino una identidad.
    """
    x = np.asarray(x, dtype=float)
    r = np.diff(x)
    a = float(np.corrcoef(r[:-1], r[1:])[0, 1])
    p1, p2, *_ = ajustar_ar2(x)
    return {
        "a_retornos": a,
        "phi1": p1, "phi2": p2,
        "phi1_pred": 1.0 + a, "phi2_pred": -a,
        "err_phi1": abs(p1 - (1.0 + a)),
        "err_phi2": abs(p2 - (-a)),
        "suma_pred": (1.0 + a) + (-a),
    }


def harvey_desde_phi(phi1: float, phi2: float) -> dict:
    """Correspondencia del §4.1: `ρ = √(−φ₂)`, `λ = arccos(φ₁/(2√(−φ₂)))`.

    Se evalua como SEGUNDA VIA INDEPENDIENTE y es gratis: si `φ₂ > 0` entonces
    `−φ₂ < 0` y `ρ` no es real, o sea que **el ciclo de Harvey no tiene solucion
    valida sobre estos datos**. No hace falta ajustar nada por maxima
    verosimilitud para saberlo.
    """
    if phi2 >= 0.0:
        return {
            "rho": float("nan"), "lambda": float("nan"),
            "existe": False,
            "motivo": (
                f"phi2 = {phi2:+.6f} >= 0 -> rho = sqrt(-phi2) no es real. "
                f"El ciclo estocastico de Harvey no tiene parametrizacion valida."
            ),
        }
    rho = math.sqrt(-phi2)
    arg = phi1 / (2.0 * rho)
    if abs(arg) > 1.0:
        return {"rho": rho, "lambda": float("nan"), "existe": False,
                "motivo": f"|phi1/(2*rho)| = {abs(arg):.4f} > 1 -> arccos indefinido"}
    return {"rho": rho, "lambda": math.acos(arg), "existe": True, "motivo": ""}


def guarda_gamma(res: dict) -> str | None:
    """§3.3: `γ < 0` es inyección neta de energía. Es una GUARDA, no un resultado.

    Si sale, el sistema es inestable o la discretización está mal: hay que parar
    y revisar antes de interpretar nada.
    """
    g = res["mco"]["gamma"]
    if g < 0.0:
        return (
            f"gamma = {g:.6f} < 0 en [{res['etiqueta']}]: INYECCION NETA DE "
            f"ENERGIA. Parar y revisar la discretizacion antes de interpretar."
        )
    return None
