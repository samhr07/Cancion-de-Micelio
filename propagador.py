"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: propagador.py — Horizonte derivado y propagador (ORDEN_TRABAJO_PROPAGADOR_3_1)

⚠ NO SE IMPORTA DESDE `Micelio.py`. §9 de la v3.1: este documento MIDE; si el
resultado lo justifica, la integración es otra versión.

QUE CAMBIA RESPECTO A LA v3.0
-----------------------------
El AR(2) de la v3.0 metía el forzamiento `F` DENTRO del residuo. El flujo de
órdenes tiene memoria larga (ley de potencias, por fragmentación), así que era un
regresor omitido y autocorrelado: sesga y confunde los coeficientes. El veredicto
`k = 0` sigue en pie para la parte NO FORZADA, pero no era un test de la premisa
de la v3.0, que era un oscilador **forzado**.

Aquí `F` está medido: `es_maker` y `q` llevan capturándose desde la v2.0 sin
usarse. Y la respuesta al impulso `G(τ)` PUEDE ser oscilatoria aunque el AR(2)
libre no lo sea — un impacto que sobrepasa y revierte ES un oscilador
amortiguado, solo que forzado por algo observado.
"""

from __future__ import annotations

import math

import numpy as np


# ==============================================================================
# §1 — EL HORIZONTE SE DERIVA DE LOS COSTES
# ==============================================================================
# ⚠ TARIFAS ASUMIDAS, NO LEIDAS. El §8 exige "escalón de comisiones leído de la
# cuenta, no asumido", y el endpoint que lo da (`/fapi/v1/commissionRate`) es
# FIRMADO: el Modo LECTURA no tiene credenciales por diseño. Se usan las tarifas
# públicas VIP 0 de USDⓈ-M y **ese criterio de aceptación queda NO CUMPLIDO**,
# declarado como tal en el preregistro en vez de fingido.
# El orden de magnitud no cambia con el escalón; el número sí.
COMISION_MAKER_VIP0 = 0.0002  # 0.0200 %
COMISION_TAKER_VIP0 = 0.0005  # 0.0500 %
COMISIONES_LEIDAS_DE_LA_CUENTA = False


def coste_ida_y_vuelta(precio: float, esquema: str = "taker_taker",
                       s_eff: float = 0.0, impacto: float = 0.0) -> float:
    """`c(u)` — coste de ida y vuelta en USD/BTC.  §1.2.

    comisiones + spread efectivo cruzado (si taker) + impacto.

    ⚠ `c` DEPENDE DEL TAMAÑO `u` a través del impacto, luego `H* = H*(u)`:
    órdenes mayores necesitan horizontes mayores. Eso acopla `H*` con la variable
    de decisión del NMPC y con la capa de Loeper, **que ya existen**. No es un
    parámetro nuevo, es una relación entre dos que ya están.
    """
    tarifas = {
        "maker_maker": 2 * COMISION_MAKER_VIP0,
        "maker_taker": COMISION_MAKER_VIP0 + COMISION_TAKER_VIP0,
        "taker_taker": 2 * COMISION_TAKER_VIP0,
    }
    comision = tarifas[esquema] * precio
    # El spread efectivo solo se cruza en las patas taker.
    n_taker = {"maker_maker": 0, "maker_taker": 1, "taker_taker": 2}[esquema]
    return comision + n_taker * s_eff + impacto


def spread_efectivo_roll(retornos: np.ndarray) -> dict:
    """Spread efectivo por el modelo de Roll: `s_eff = 2·σ_r·√(−ρ₁)`.  §1.2.

    ⚠ `σ_r` DEBE medirse sobre datos LIMPIOS. Con los ceros del feed dentro, el
    |incremento| máximo pasa de 11.8 a 65 245 USD y σ_r queda inflada por un
    factor enorme (§0 de la v3.1).

    El modelo de Roll es exactamente el rebote bid-ask que la v3.0 identificó
    como explicación completa del AR(2): aquí se usa para ponerle precio.
    """
    r = np.asarray(retornos, dtype=float)
    r = r[np.isfinite(r)]
    if len(r) < 100:
        return {"s_eff": float("nan"), "rho1": float("nan"), "sigma_r": float("nan")}
    rho1 = float(np.corrcoef(r[:-1], r[1:])[0, 1])
    sigma_r = float(np.std(r, ddof=1))
    # Roll clásico: s = 2·√(−cov(r_t, r_{t-1})). Equivale a 2σ√(−ρ₁).
    s_eff = 2.0 * sigma_r * math.sqrt(-rho1) if rho1 < 0 else float("nan")
    return {"s_eff": s_eff, "rho1": rho1, "sigma_r": sigma_r}


def firma_de_volatilidad(precios: np.ndarray, tiempos: np.ndarray,
                         horizontes_s=(0.5, 1, 2, 5, 10, 30, 60, 120, 300, 600)):
    """σ(H) contra H — el GRAFICO DE FIRMA DE VOLATILIDAD del §8.

    Dice EMPIRICAMENTE a qué Δτ deja de dominar el ruido de microestructura, que
    es la misma pregunta que la elección de `K` y que hasta ahora se resolvía por
    conjetura. Si σ(H)/√H cae al crecer H, el ruido de microestructura domina en
    la parte corta; se aplana cuando deja de hacerlo.
    """
    p = np.asarray(precios, dtype=float)
    t = np.asarray(tiempos, dtype=float)
    filas = []
    for H in horizontes_s:
        # Muestreo por reloj de pared al horizonte H, con el precio vigente.
        bordes = np.arange(t[0], t[-1], H)
        if len(bordes) < 30:
            continue
        idx = np.clip(np.searchsorted(t, bordes, side="right") - 1, 0, len(p) - 1)
        s = p[idx]
        r = np.diff(s)
        if len(r) < 20:
            continue
        sigma = float(np.std(r, ddof=1))
        filas.append({
            "H": float(H),
            "n": len(r),
            "sigma": sigma,
            "sigma_por_raiz_H": sigma / math.sqrt(H),
        })
    return filas


def horizonte_derivado(filas_firma, c_usd_btc: float) -> dict:
    """`H*` resuelve `σ(H*) = c(u)`.  §1.2.

    Se interpola sobre la firma de volatilidad medida en vez de suponer
    `σ ∝ √H`: el §1.1 da la tabla bajo esa suposición, pero la firma real puede
    apartarse de ella justo en la escala corta, que es donde vive el problema.
    """
    if not filas_firma:
        return {"H_estrella": float("nan"), "metodo": "sin datos"}
    H = np.array([f["H"] for f in filas_firma])
    S = np.array([f["sigma"] for f in filas_firma])
    if S[-1] < c_usd_btc:
        return {
            "H_estrella": float("inf"),
            "metodo": "extrapolado",
            "aviso": (
                f"sigma no alcanza c = {c_usd_btc:.2f} USD/BTC dentro del rango "
                f"medido (max sigma = {S[-1]:.2f} a H = {H[-1]:.0f} s). H* esta "
                f"FUERA de los datos: el sistema deberia operar MENOS, no mas."
            ),
        }
    i = int(np.searchsorted(S, c_usd_btc))
    if i == 0:
        return {"H_estrella": float(H[0]), "metodo": "por debajo del rango"}
    # Interpolación log-log, que es la natural para una ley de potencias.
    lo, hi = i - 1, i
    w = (math.log(c_usd_btc) - math.log(S[lo])) / (math.log(S[hi]) - math.log(S[lo]))
    H_est = math.exp(math.log(H[lo]) + w * (math.log(H[hi]) - math.log(H[lo])))
    return {"H_estrella": float(H_est), "metodo": "interpolado log-log"}


# ==============================================================================
# §2 — EL PROPAGADOR CON FORZAMIENTO MEDIDO
# ==============================================================================
def signo_transaccion(es_maker) -> np.ndarray:
    """ε a partir del campo `m` de Binance.  §2.2.

    ⚠ CONVENCION, Y ES LA TERCERA VEZ QUE ESTE PROYECTO SE JUEGA UN RESULTADO EN
    UN SIGNO (tras el 2π y el factor 125). El campo `m` indica si el COMPRADOR
    fue el maker:

        m = True   -> comprador maker  -> el taker VENDIA   -> eps = -1
        m = False  -> comprador taker  -> el taker COMPRABA -> eps = +1

    La predicción falsable está en el preregistro: **`G(0) > 0`**. Si sale
    negativo, el signo está invertido y SE PARA; no se le da la vuelta y se sigue.
    """
    m = np.asarray(es_maker).astype(bool)
    return np.where(m, -1.0, 1.0)


def propagador_parametrico(tau, G0, tau0, beta):
    """`G(τ) = G₀·(1 + τ/τ₀)^(−β)`.  §2.3.

    Tres parámetros en vez de ~3 300 rezagos libres: a ν ≈ 102 tx/s y H* ≈ 32 s,
    cubrir el horizonte con rezagos libres sería inestimable. Es la dirección que
    pide el proyecto — menos constantes, y medidas en vez de elegidas.
    """
    return G0 * np.power(1.0 + np.asarray(tau, dtype=float) / tau0, -beta)


def ajustar_respuesta_impulso(delta_p, eps, vol, max_rezago=None, delta=0.0):
    """Respuesta media del precio a un impulso de signo, rezago a rezago.

    Estima `R(τ) = E[(p_{t+τ} − p_t) · ε_t] / E[|f(v_t)|]`, la respuesta
    acumulada, que es la forma empírica estándar del propagador y no exige
    invertir ninguna matriz grande.

    `f(v) = v^δ` con `δ` barrida (§2.3).
    """
    dp = np.asarray(delta_p, dtype=float)  # precio (nivel), no incrementos
    e = np.asarray(eps, dtype=float)
    v = np.asarray(vol, dtype=float)
    f = np.power(np.maximum(v, 1e-12), delta) if delta != 0.0 else np.ones_like(v)
    peso = e * f
    n = len(dp)
    if max_rezago is None:
        max_rezago = min(2000, n // 10)
    norm = float(np.mean(np.abs(f)))
    R = np.empty(max_rezago + 1)
    for tau in range(max_rezago + 1):
        # E[(p_{t+tau} - p_t) * eps_t]
        R[tau] = float(np.mean((dp[tau : n] - dp[0 : n - tau]) * peso[0 : n - tau]))
    return R / max(norm, 1e-12)


def ajustar_ley_potencias(tau, y, tau_min=1):
    """Ajuste log-log de `y ≈ A·τ^(−b)`. Devuelve (A, b, R²)."""
    tau = np.asarray(tau, dtype=float)
    y = np.asarray(y, dtype=float)
    sel = (tau >= tau_min) & (y > 0) & np.isfinite(y)
    if sel.sum() < 5:
        return float("nan"), float("nan"), float("nan")
    X = np.column_stack([np.log(tau[sel]), np.ones(sel.sum())])
    b, *_ = np.linalg.lstsq(X, np.log(y[sel]), rcond=None)
    pred = X @ b
    ss_res = float(np.sum((np.log(y[sel]) - pred) ** 2))
    ss_tot = float(np.sum((np.log(y[sel]) - np.mean(np.log(y[sel]))) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(math.exp(b[1])), float(-b[0]), r2


# ==============================================================================
# §3 SUSTITUIDO POR LA ADENDA A — LA PENDIENTE DE LA FIRMA DE VOLATILIDAD
# ==============================================================================
# ⚠ EL SIMULADOR SE ELIMINA. El A.5 lo marca como motivo de rechazo del
# entregable, y con razon: el §3 original contrastaba una condicion ASINTOTICA DE
# BAJA FRECUENCIA (`β = (1−γ)/2`, que sale de `S_p(f) ~ f^(2β+γ−3)` cuando
# `f → 0`) usando `ρ₁`, un estadistico de REZAGO 1 dominado por las frecuencias
# ALTAS. Estadistico de un regimen para contrastar una condicion del otro.
#
# El sintoma que lo delato: el "sesgo de tamano finito" CAMBIABA DE SIGNO
# (+0.113, +0.063, −0.042 contra la gamma medida). Un sesgo de tamano finito no
# hace eso. El cruce no seguia a la teoria, seguia a otra cosa.
#
# LO QUE LO REEMPLAZA, del mismo desarrollo espectral: el exponente de Hurst del
# precio es `H_p = 1 − β − γ/2`, y como `σ(H) ~ H^(H_p)`:
#
#       d log[σ(H)/√H] / d log H  =  H_p − ½  =  (1−γ)/2 − β
#
# **La pendiente de la firma de volatilidad ES la desviacion respecto de la
# condicion de difusividad.** Directamente, sobre datos reales, sin simulador y
# sin umbral calibrado.
#
#       pendiente = 0  ->  H_p = 0.50  ->  difusivo, SIN VENTAJA EXPLOTABLE
#       pendiente < 0  ->  H_p < 0.50  ->  subdifusivo -> REVERSION
#       pendiente > 0  ->  H_p > 0.50  ->  superdifusivo -> MOMENTUM
#
# Para la DIRECCION no hace falta estimar ni `β` ni `γ`. La pendiente sola basta.
# `β` y `γ` siguen haciendo falta para saber si el propagador es el MECANISMO,
# pero no para el veredicto. Dos niveles, y el que decide es el barato.
def firma_solapada(precios, tiempos, horizontes_s):
    """σ(H) por ventanas SOLAPADAS. A.3.1.

    Mas eficiente que las no solapadas; la dependencia que introduce el solape la
    absorbe el bootstrap de bloques moviles del A.3.3, que es donde debe tratarse.
    """
    p = np.asarray(precios, dtype=float)
    t = np.asarray(tiempos, dtype=float)
    filas = []
    for H in horizontes_s:
        # Para cada t_i, el precio a t_i + H. Ventanas solapadas: una por tick.
        j = np.searchsorted(t, t + H, side="right") - 1
        val = (j > np.arange(len(t))) & (j < len(t))
        if val.sum() < 30:
            continue
        r = p[j[val]] - p[val]
        n_indep = max(1, int((t[-1] - t[0]) / H))
        filas.append({
            "H": float(H),
            "n_solapadas": int(val.sum()),
            "n_indep": n_indep,
            "sigma": float(np.std(r, ddof=1)),
            "sigma_por_raiz_H": float(np.std(r, ddof=1) / math.sqrt(H)),
        })
    return filas


def determinar_H_lo(filas):
    """`H_lo` = donde `σ(H)/√H` DEJA DE CRECER. A.3.2, punto 1.

    ⚠ Se determina POR LA FIRMA MISMA, no por criterio, y se commitea ANTES del
    ajuste. El rango no puede elegirse mirando donde sale la pendiente que gusta:
    ese es exactamente el defecto que la v2.2 documento al medir `C = 0.783`
    sobre datos barajados.
    """
    if len(filas) < 3:
        return float("nan")
    H = np.array([f["H"] for f in filas])
    y = np.array([f["sigma_por_raiz_H"] for f in filas])
    i_max = int(np.argmax(y))
    return float(H[i_max])


def pendiente_firma(filas, H_lo, H_hi):
    """Regresion de `log[σ(H)/√H]` sobre `log H` en el rango dado."""
    H = np.array([f["H"] for f in filas])
    y = np.array([f["sigma_por_raiz_H"] for f in filas])
    sel = (H >= H_lo) & (H <= H_hi) & (y > 0)
    if sel.sum() < 3:
        return {"pendiente": float("nan"), "n": int(sel.sum())}
    X = np.column_stack([np.log(H[sel]), np.ones(sel.sum())])
    b, *_ = np.linalg.lstsq(X, np.log(y[sel]), rcond=None)
    return {
        "pendiente": float(b[0]),
        "H_p": float(b[0] + 0.5),
        "n": int(sel.sum()),
        "rango": (float(H_lo), float(H_hi)),
    }


def bootstrap_pendiente(precios, tiempos, horizontes_s, H_lo, H_hi, rng,
                        n_boot=300, largo_bloque=None):
    """IC de la pendiente por BOOTSTRAP DE BLOQUES MOVILES. A.3.3.

    ⚠ NO se usa la formula de Lo-MacKinlay ni ningun error estandar gaussiano.
    Con curtosis 1179.7 esa asintotica es inservible, y usarla seria la TERCERA
    aparicion del mismo patron en este proyecto (tras los umbrales chi^2 sobre el
    NIS y la razon de verosimilitudes de Harvey en la frontera).

    Longitud de bloque >= 5*H_hi para preservar la memoria larga del flujo.
    """
    p = np.asarray(precios, dtype=float)
    t = np.asarray(tiempos, dtype=float)
    dur = t[-1] - t[0]
    if largo_bloque is None:
        largo_bloque = 5.0 * H_hi
    n_bloques = max(1, int(dur / largo_bloque))
    pend = []
    for _ in range(n_boot):
        trozos_p, trozos_t, t_acum = [], [], 0.0
        for _ in range(n_bloques):
            ini = float(rng.uniform(t[0], max(t[0], t[-1] - largo_bloque)))
            sel = (t >= ini) & (t < ini + largo_bloque)
            if sel.sum() < 100:
                continue
            # Se reancla el tiempo para poder concatenar sin fabricar huecos.
            tt = t[sel] - t[sel][0] + t_acum
            trozos_t.append(tt)
            trozos_p.append(p[sel])
            t_acum = tt[-1] + 1e-3
        if len(trozos_p) < 2:
            continue
        pb = np.concatenate(trozos_p)
        tb = np.concatenate(trozos_t)
        f = firma_solapada(pb, tb, horizontes_s)
        r = pendiente_firma(f, H_lo, H_hi)
        if math.isfinite(r["pendiente"]):
            pend.append(r["pendiente"])
    pend = np.asarray(pend)
    if pend.size < 20:
        return {"n_boot": int(pend.size), "ic95": (float("nan"), float("nan"))}
    return {
        "n_boot": int(pend.size),
        "n_bloques": n_bloques,
        "largo_bloque_s": float(largo_bloque),
        "ic95": (float(np.percentile(pend, 2.5)), float(np.percentile(pend, 97.5))),
        # p-valor bootstrap de dos colas contra pendiente = 0
        "p_contra_cero": float(2.0 * min(np.mean(pend >= 0.0), np.mean(pend <= 0.0))),
        "mediana": float(np.median(pend)),
    }


def razon_de_varianzas(filas, H_ref):
    """`VR(H) = [σ(H)/√H]² / [σ(H_ref)/√H_ref]²`. Lectura directa del A.2."""
    d = {f["H"]: f["sigma_por_raiz_H"] for f in filas}
    if H_ref not in d:
        return {}
    base = d[H_ref] ** 2
    return {H: (v ** 2) / base for H, v in d.items() if base > 0}
