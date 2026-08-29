"""
Proyecto: Canción del Micelio — Cribado de activos (ORDEN_TRABAJO_CRIBADO_4_0)
Módulo: cribado_activos.py — Etapa 1, cribado sin captura.

⚠ NO SE IMPORTA DESDE `Micelio.py`. Igual que `propagador.py` y `oscilador.py`:
este módulo MIDE. Si el resultado lo justifica, la integración es otra versión.

QUE HACE
--------
Mide, por cada par USDⓈ-M, los seis números que deciden si vale la pena capturarlo,
y aplica la **compuerta 1** de `PREREGISTRO_CRIBADO_4_0.md` §4 con los umbrales
declarados ANTES de correr. Después contrasta la predicción falsable del §2:
`σ₁` y `R²` predictivo NO deben estar negativamente correlacionados entre activos.

POR QUE ESTE CRIBADO Y NO OTRO
------------------------------
De la identidad económica (§3.1 del preregistro), con la comisión de futuros
USDⓈ-M **idéntica en todos los pares**, la única variable libre entre pares es
`σ₁`. Todo lo demás que se mide aquí está para impedir que un `σ₁` alto engañe:

  - **Roll** — σ de velas de 1 min SOBREESTIMA por rebote bid-ask, y sobreestima
    MAS cuanto más ancha es la horquilla, que es justo la dirección que favorece a
    los pares malos. La compuerta usa el σ₁ CORREGIDO (§4, C1.1).
  - **fracción del movimiento en el 1 % de los minutos** y **curtosis** — un activo
    cuya volatilidad es toda saltos discretos no sirve: el decil superior no se
    captura con órdenes maker (§6.3).
  - **horquilla y granularidad** — un margen que existe pero no cabe entre el
    tick y el spread no es un margen.
  - **control barajado** — bajo barajado `H_p` debe salir 0.5. Si no sale, el
    estimador tiene sesgo propio y su `H_p` no se reporta como medida (§6.2).
    Es el control del A.3.4, que ya salvó dos análisis en el proyecto viejo.

Sin dependencias más allá de NumPy: la suite del proyecto corre en cualquier
entorno donde corra el bot, y eso incluye no tener SciPy.
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
import urllib.error
import urllib.request

import numpy as np


# ==============================================================================
# CONSTANTES DECLARADAS EN EL PREREGISTRO — NO SE TOCAN AQUI
# ==============================================================================
# Cualquier cambio a estos valores es un cambio al preregistro y necesita commit
# fechado que diga que cambio y por que (§ cabecera del PREREGISTRO_CRIBADO_4_0).

KAPPA = math.sqrt(2.0 / math.pi)          # E[y*signo(y_hat)] = kappa*rho*sigma_H
COMISION_MAKER_PB = 2.000                  # USDⓈ-M VIP 0, 0.0200 %  -> pb por lado
COMISION_TAKER_PB = 5.000                  # USDⓈ-M VIP 0, 0.0500 %  -> pb por lado
COMISIONES_LEIDAS_DE_LA_CUENTA = False     # §7.3: /fapi/v1/commissionRate es firmado

H_REFERENCIA_S = 3600.0                    # el horizonte al que se declara el umbral

# --- Compuerta 1 (§4 del preregistro). Umbrales declarados antes de medir. -----
SIGMA1_EQUILIBRIO_PB = 1.30                # sigma1 que hace R2_req = R2_medido en BTC
MARGEN_COMPUERTA = 2.0                     # "con margen de 2x"
C1_1_SIGMA1_MIN_PB = SIGMA1_EQUILIBRIO_PB * MARGEN_COMPUERTA   # = 2.60
C1_2_HORQUILLA_MAX_PB = 1.00
C1_3_GRANULARIDAD_MAX_USD = 5.00
C1_3_LOTES_MIN = 40                        # Sec. A.3 de la v1.3
C1_4_FRACCION_SALTO_MAX = 0.35
C1_5_CURTOSIS_MAX = 60.0
C1_6_COBERTURA_MIN = 0.95

# --- Predicción falsable (§2 del preregistro). ---------------------------------
SPEARMAN_MATA_HIPOTESIS = -0.30
SPEARMAN_ALFA = 0.05

# --- Control barajado (§6.2). --------------------------------------------------
TOLERANCIA_HP_BARAJADO = 0.05              # |H_p_barajado - 0.5| > 0.05 -> sesgo propio

# --- Nocional de orden, para la comprobación de lotes de C1.3. -----------------
NOCIONAL_MAX_ORDEN_USD = 3000.0            # el de la v1.3, Sec. B

# --- Centinelas: las cifras de BTC del §3.2 NO estan en este arbol (§7.2). ------
# Nacen en CERO a proposito. Olvidarse de fijarlas falla ruidosamente en vez de
# operar con un literal inventado, que es la regla de la Sec. A de la v1.3.
BTC_R2_DIRECCIONAL = 0.0                   # [DECLARAR] R2 fuera de muestra, banda 0.002-0.008
BTC_SIGMA1_ESTACIONAL_0_PB = 0.0           # [DECLARAR]
BTC_SIGMA1_ESTACIONAL_1_PB = 0.0           # [DECLARAR]

HORIZONTES_S = (60, 120, 180, 300, 600, 900, 1800, 3600)
H_PREDICCION_MIN = 60                      # 1 h, el horizonte de referencia del §3.3


# ==============================================================================
# §A — LA IDENTIDAD ECONOMICA
# ==============================================================================
def r2_requerido(sigma1_pb: float, h_p: float, horizonte_s: float = H_REFERENCIA_S,
                 c_lado_pb: float = COMISION_MAKER_PB) -> float:
    """`R²_req` — el R² fuera de muestra que hace nulo el margen.  §3.1.

        E[y_H * signo(y_hat)] = kappa * rho * sigma_H  =  2 * c_lado
        sigma_H = sigma_1 * H**H_p

    Es **lineal en `c_lado²` e inversa en `σ₁²`**: por eso, con la comisión igual
    en todos los pares USDⓈ-M, el cribado se ordena por `σ₁` y por nada más.
    """
    sigma_h = sigma1_pb * (horizonte_s ** h_p)
    if sigma_h <= 0.0:
        return float("inf")
    return (2.0 * c_lado_pb / (KAPPA * sigma_h)) ** 2


def comision_de_equilibrio_pb(sigma1_pb: float, h_p: float, r2_medido: float,
                              horizonte_s: float = H_REFERENCIA_S) -> float:
    """`c_lado` que anula el margen dado el `R²` medido.  Inversa de `r2_requerido`."""
    if r2_medido <= 0.0:
        return 0.0
    return KAPPA * math.sqrt(r2_medido) * sigma1_pb * (horizonte_s ** h_p) / 2.0


# ==============================================================================
# §B — LIMPIEZA (§6.1 del preregistro)
# ==============================================================================
def limpiar_klines(k: np.ndarray) -> dict:
    """Descarta velas vacias ANTES de medir varianza.  §6.1.

    ⚠ El analogo de los trades con `p = 0` de la v3.0. Y la leccion de la v3.0 es
    que los datos sucios ENMASCARAN: al limpiar, la curtosis SUBE (601.8 -> 1179.7)
    porque los ceros inflaban sigma^4 y aplastaban el cociente. Se limpia antes de
    cualquier estadistico de dispersion, no despues.

    `k` es (n, 7): [t_apertura_ms, apertura, max, min, cierre, volumen, n_trades,
    volumen_taker_compra].  (n, 8) en realidad; ver `KLINE_COLS`.
    """
    if k.size == 0:
        return {"k": k, "n_bruto": 0, "n_limpio": 0, "cobertura": 0.0}
    n_bruto = len(k)
    cierre = k[:, 4]
    volumen = k[:, 5]
    n_trades = k[:, 6]
    valido = (
        np.isfinite(cierre) & (cierre > 0.0)
        & np.isfinite(volumen) & (volumen > 0.0)
        & np.isfinite(n_trades) & (n_trades > 0.0)
    )
    kl = k[valido]
    # Cobertura: minutos observados contra minutos esperados en el tramo.
    if len(kl) >= 2:
        span_min = (kl[-1, 0] - kl[0, 0]) / 60000.0 + 1.0
        cobertura = len(kl) / span_min if span_min > 0 else 0.0
    else:
        cobertura = 0.0
    return {
        "k": kl,
        "n_bruto": n_bruto,
        "n_limpio": len(kl),
        "descartadas": n_bruto - len(kl),
        "cobertura": float(min(cobertura, 1.0)),
    }


# ==============================================================================
# §C — FIRMA DE VOLATILIDAD: sigma_1 y H_p
# ==============================================================================
def firma_sigma(log_precios: np.ndarray, paso_s: float = 60.0,
                horizontes_s=HORIZONTES_S) -> list:
    """`σ(H)` contra `H`, en pb, sobre retornos NO SOLAPADOS.

    No solapados a proposito: la Adenda A del proyecto viejo midio que la version
    solapada SOBREPONDERA los tramos de alta actividad —que son tambien los de alta
    volatilidad— y que solapada y no solapada **no estiman la misma cantidad**
    (discrepaban hasta en el signo de la pendiente). Aqui se necesita una sola
    cantidad bien definida, asi que se usa la no solapada.
    """
    p = np.asarray(log_precios, dtype=float)
    filas = []
    for H in horizontes_s:
        paso = int(round(H / paso_s))
        if paso < 1:
            continue
        s = p[::paso]
        r = np.diff(s)
        r = r[np.isfinite(r)]
        if len(r) < 30:
            continue
        filas.append({
            "H": float(H),
            "n": int(len(r)),
            "sigma_pb": float(np.std(r, ddof=1) * 1e4),
        })
    return filas


def ajustar_ley_sigma(filas: list) -> dict:
    """Ajusta `log σ = log σ₁ + H_p · log H` por MCO.  Devuelve σ₁ [pb·s^(−H_p)].

    Se reporta el `R²` del propio ajuste: si la ley de potencias no describe la
    firma, `σ₁` y `H_p` no significan nada por separado y hay que verlo.
    """
    filas = [f for f in filas if f["sigma_pb"] > 0 and np.isfinite(f["sigma_pb"])]
    if len(filas) < 3:
        return {"sigma1_pb": float("nan"), "H_p": float("nan"),
                "r2_ajuste": float("nan"), "n_puntos": len(filas)}
    x = np.log(np.array([f["H"] for f in filas]))
    y = np.log(np.array([f["sigma_pb"] for f in filas]))
    A = np.column_stack([np.ones_like(x), x])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))

    # ⚠ `sigma_1` ES UNA EXTRAPOLACION, y hay que reportarla como tal.
    # El ajuste vive en [60, 3600] s (rango declarado en el preregistro) y `sigma_1`
    # es su valor en H = 1 s: SESENTA VECES por debajo del punto mas corto ajustado.
    # Un error `d` en H_p se convierte en un factor 60^(-d) sobre sigma_1, o sea que
    # d = 0.012 ya son ~5 % en sigma_1. Medido sobre paseos sinteticos con verdad
    # conocida, la dispersion de sigma_1 entre semillas es de ~3.5 % con H_p bueno a
    # 0.008 (ver `test_v40_firma_recupera_sigma1_y_Hp_de_verdad_conocida`).
    # No se cambia el rango de ajuste —esta declarado— pero se publica el error, para
    # que un par cerca del umbral C1.1 se vea que esta cerca en vez de parecer nitido.
    gl = len(filas) - 2
    if gl > 0 and ss_res > 0:
        s2 = ss_res / gl
        cov = s2 * np.linalg.inv(A.T @ A)
        ee_intercepto = float(math.sqrt(max(cov[0, 0], 0.0)))
        ee_pendiente = float(math.sqrt(max(cov[1, 1], 0.0)))
    else:
        ee_intercepto = ee_pendiente = float("nan")
    return {
        "sigma1_pb": float(math.exp(coef[0])),
        "H_p": float(coef[1]),
        # sigma_1 = exp(b0)  =>  error relativo de sigma_1 ~ error absoluto de b0.
        # ⚠ COTA INFERIOR: los residuos de horizontes anidados salen de los MISMOS
        # datos y estan correlacionados, asi que el error de MCO subestima.
        "sigma1_ee_rel": ee_intercepto,
        "H_p_ee": ee_pendiente,
        "r2_ajuste": float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan"),
        "n_puntos": len(filas),
    }


# ==============================================================================
# §D — ROLL: quitar el rebote bid-ask de sigma
# ==============================================================================
def roll_desde_retornos(retornos: np.ndarray) -> dict:
    """Spread efectivo relativo por Roll: `s = 2·√(−cov₁)`.  Devuelve pb.

    Mismo estimador que `propagador.spread_efectivo_roll`, pero sobre retornos
    RELATIVOS (para poder comparar pares de precios muy distintos) y devolviendo
    tambien la contribucion de varianza que hay que restar.

    Roll: `cov(r_t, r_{t-1}) = −(s/2)²`, y el rebote inyecta en la varianza de
    CUALQUIER horizonte una cantidad fija `s²/2` (es un error de medida iid en los
    dos extremos del intervalo). Por eso muerde mas cuanto mas corto es H.
    """
    r = np.asarray(retornos, dtype=float)
    r = r[np.isfinite(r)]
    if len(r) < 100:
        return {"s_eff_pb": float("nan"), "rho1": float("nan"),
                "sigma_r_pb": float("nan"), "var_rebote_pb2": float("nan")}
    rc = r - r.mean()
    cov1 = float(np.dot(rc[1:], rc[:-1]) / (len(rc) - 1))
    var = float(np.dot(rc, rc) / (len(rc) - 1))
    rho1 = cov1 / var if var > 0 else float("nan")
    if cov1 < 0.0:
        s_eff = 2.0 * math.sqrt(-cov1)          # relativo
        var_rebote = s_eff * s_eff / 2.0        # relativo^2
    else:
        # rho1 >= 0: el modelo de Roll no aplica. NO se inventa un spread.
        s_eff = float("nan")
        var_rebote = 0.0
    return {
        "s_eff_pb": float(s_eff * 1e4),
        "rho1": float(rho1),
        "sigma_r_pb": float(math.sqrt(var) * 1e4),
        "var_rebote_pb2": float(var_rebote * 1e8),   # (rel^2) -> pb^2
    }


def corregir_firma_por_roll(filas: list, var_rebote_pb2: float) -> dict:
    """Resta `s²/2` de la varianza de CADA horizonte y reajusta la ley.

    Si a algun horizonte la varianza corregida sale <= 0, ese horizonte esta
    DOMINADO por el rebote y se descarta con aviso, en vez de producir un sigma
    imaginario o un cero silencioso.
    """
    if not np.isfinite(var_rebote_pb2) or var_rebote_pb2 <= 0.0:
        return {"filas": list(filas), "descartados": [], "aplicada": False}
    corregidas, descartados = [], []
    for f in filas:
        v = f["sigma_pb"] ** 2 - var_rebote_pb2
        if v <= 0.0:
            descartados.append(f["H"])
            continue
        g = dict(f)
        g["sigma_pb"] = math.sqrt(v)
        corregidas.append(g)
    return {"filas": corregidas, "descartados": descartados, "aplicada": True}


# ==============================================================================
# §E — R2 PREDICTIVO FUERA DE MUESTRA (proxy declarado, §6.4)
# ==============================================================================
def rasgos_y_objetivo(k: np.ndarray, ventana_min: int = H_PREDICCION_MIN) -> tuple:
    """Agrega por bloques NO SOLAPADOS de `ventana_min` minutos.

    Rasgos del bloque `i` (todos observables al cerrar el bloque):
      0. desequilibrio de flujo   (2*vol_taker_compra - vol) / vol      <- el eps del propagador
      1. retorno del bloque
      2. retorno del bloque anterior
      3. rango medio (max-min)/cierre
      4. log(n_trades)
      5. log(volumen)
    Objetivo: retorno del bloque `i+1`.

    ⚠ NO SOLAPADO a proposito. Con objetivo a 60 min sobre datos de 1 min, el
    solapamiento inflaria el R2 por muestras casi repetidas — la misma familia del
    defecto de "corregir 90 veces con la misma medicion" de la v1.3.
    """
    cierre, maxi, mini = k[:, 4], k[:, 2], k[:, 3]
    volumen, n_trades, vol_taker = k[:, 5], k[:, 6], k[:, 7]
    n_bloques = len(k) // ventana_min
    if n_bloques < 20:
        return np.zeros((0, 6)), np.zeros(0)
    corte = n_bloques * ventana_min
    rs = lambda a: a[:corte].reshape(n_bloques, ventana_min)
    c_bloque, v_bloque = rs(cierre), rs(volumen)
    vt_bloque, nt_bloque = rs(vol_taker), rs(n_trades)
    hi_bloque, lo_bloque = rs(maxi), rs(mini)

    cierre_fin = c_bloque[:, -1]
    # ⚠ El bloque 0 no tiene cierre anterior con el que medir su retorno. Usar su
    # propia primera vela daria un retorno CORTO (le falta el primer minuto) y lo
    # mezclaria con los demas como si fuera comparable. Se marca NaN y la mascara
    # de abajo lo descarta, junto con el bloque 1 (que lo lleva como rezago).
    # Cuesta 2 bloques de ~2160; inventarse la apertura cuesta un sesgo silencioso.
    apertura = np.concatenate([[np.nan], cierre_fin[:-1]])
    with np.errstate(divide="ignore", invalid="ignore"):
        ret = np.log(cierre_fin / apertura)

    vol_tot = v_bloque.sum(axis=1)
    vt_tot = vt_bloque.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        ofi = np.where(vol_tot > 0, (2.0 * vt_tot - vol_tot) / vol_tot, 0.0)
        rango = np.mean((hi_bloque - lo_bloque) / np.maximum(c_bloque, 1e-12), axis=1)
    nt_tot = nt_bloque.sum(axis=1)

    X = np.column_stack([
        ofi,
        ret,
        np.concatenate([[np.nan], ret[:-1]]),
        rango,
        np.log(np.maximum(nt_tot, 1.0)),
        np.log(np.maximum(vol_tot, 1e-12)),
    ])
    y = np.concatenate([ret[1:], [np.nan]])      # objetivo = retorno del bloque siguiente
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    return X[ok], y[ok]


def r2_predictivo_oos(X: np.ndarray, y: np.ndarray, fraccion_train: float = 0.70) -> dict:
    """`R²` fuera de muestra por particion TEMPORAL 70/30 sin solapamiento.  §6.4.

    El denominador usa la media de ENTRENAMIENTO, no la de test: usar la de test
    seria mirar el futuro para fijar el baseline, y es la trampa clasica del R2
    fuera de muestra. Con esa convencion `R²_oos` puede salir NEGATIVO, y debe
    poder: significa que el predictor es peor que la media conocida de antemano.
    """
    n = len(y)
    if n < 40:
        return {"r2_oos": float("nan"), "n_train": 0, "n_test": 0}
    corte = int(n * fraccion_train)
    Xtr, ytr = X[:corte], y[:corte]
    Xte, yte = X[corte:], y[corte:]
    if len(yte) < 10:
        return {"r2_oos": float("nan"), "n_train": len(ytr), "n_test": len(yte)}
    # Estandarizar con estadisticos de TRAIN (de nuevo: nada de test entra al ajuste).
    mu, sd = Xtr.mean(axis=0), Xtr.std(axis=0)
    sd = np.where(sd > 0, sd, 1.0)
    Atr = np.column_stack([np.ones(len(Xtr)), (Xtr - mu) / sd])
    Ate = np.column_stack([np.ones(len(Xte)), (Xte - mu) / sd])
    coef, *_ = np.linalg.lstsq(Atr, ytr, rcond=None)
    pred = Ate @ coef
    ss_res = float(np.sum((yte - pred) ** 2))
    ss_tot = float(np.sum((yte - ytr.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {"r2_oos": float(r2), "n_train": len(ytr), "n_test": len(yte)}


# ==============================================================================
# §F — SALTO CONTRA VOLATILIDAD (C1.4 y C1.5)
# ==============================================================================
def metricas_de_salto(retornos: np.ndarray, cola: float = 0.01) -> dict:
    """Cuanto del movimiento vive en el 1 % de los minutos, y curtosis.

    ⚠ La distincion que da sentido a la compuerta: volatilidad **difusa** se puede
    capturar con ordenes maker; volatilidad que es toda **salto** no, porque el
    salto se lleva por delante la orden pasiva en vez de llenarla al precio puesto.
    Dos activos con el mismo sigma_1 y distinta concentracion NO valen lo mismo.
    """
    r = np.asarray(retornos, dtype=float)
    r = r[np.isfinite(r)]
    if len(r) < 200:
        return {"fraccion_cola": float("nan"), "curtosis": float("nan"), "n": len(r)}
    a = np.abs(r)
    total = float(a.sum())
    k_cola = max(1, int(round(len(a) * cola)))
    top = np.sort(a)[-k_cola:]
    frac = float(top.sum() / total) if total > 0 else float("nan")
    rc = r - r.mean()
    var = float(np.mean(rc ** 2))
    curt = float(np.mean(rc ** 4) / var ** 2 - 3.0) if var > 0 else float("nan")
    return {"fraccion_cola": frac, "curtosis": curt, "n": int(len(r)),
            "n_cola": k_cola}


# ==============================================================================
# §G — ESTADISTICA SIN SCIPY
# ==============================================================================
def _betacf(a: float, b: float, x: float, itmax: int = 300, eps: float = 3e-12) -> float:
    """Fraccion continua de Lentz para la beta incompleta.  Numerical Recipes 6.4."""
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < 1e-300:
        d = 1e-300
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-300:
            d = 1e-300
        c = 1.0 + aa / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def betainc_reg(a: float, b: float, x: float) -> float:
    """Beta incompleta regularizada `I_x(a,b)`."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _betacf(a, b, x) / a
    return 1.0 - front * _betacf(b, a, 1.0 - x) / b


def t_sf(t: float, df: float) -> float:
    """Cola superior bilateral de una t de Student: `P(|T| > |t|)`."""
    if df <= 0 or not math.isfinite(t):
        return float("nan")
    return betainc_reg(df / 2.0, 0.5, df / (df + t * t))


def spearman(x, y) -> dict:
    """`ρ` de Spearman con p-valor bilateral por la aproximacion t.

    Es el estadistico declarado en el §2 del preregistro para la prediccion
    falsable. De ORDEN a proposito: el `R²` del cribado es un proxy (§6.4) y no
    se puede comparar en valor absoluto contra el de BTC, pero su ORDEN entre
    pares si es informativo.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 5:
        return {"rho": float("nan"), "p": float("nan"), "n": n}
    rx, ry = _rangos(x), _rangos(y)
    rxc, ryc = rx - rx.mean(), ry - ry.mean()
    den = math.sqrt(float(np.dot(rxc, rxc)) * float(np.dot(ryc, ryc)))
    if den <= 0:
        return {"rho": float("nan"), "p": float("nan"), "n": n}
    rho = float(np.dot(rxc, ryc) / den)
    if abs(rho) >= 1.0:
        return {"rho": rho, "p": 0.0, "n": n}
    t = rho * math.sqrt((n - 2) / (1.0 - rho * rho))
    return {"rho": rho, "p": float(t_sf(t, n - 2)), "n": n,
            "rho_minimo_detectable": float(_rho_min_detectable(n))}


def _rangos(a: np.ndarray) -> np.ndarray:
    """Rangos con promedio en los empates (convencion estandar de Spearman)."""
    orden = np.argsort(a, kind="mergesort")
    rangos = np.empty(len(a), dtype=float)
    rangos[orden] = np.arange(1, len(a) + 1, dtype=float)
    # Promediar empates.
    vals = a[orden]
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[j + 1] == vals[i]:
            j += 1
        if j > i:
            rangos[orden[i:j + 1]] = rangos[orden[i:j + 1]].mean()
        i = j + 1
    return rangos


def _rho_min_detectable(n: int, alfa: float = SPEARMAN_ALFA) -> float:
    """`|ρ|` minimo que alcanza significacion con `n` pares.

    Se reporta SIEMPRE junto a rho, por la misma razon que `_rho_minimo_detectable`
    en `diagnostico.py`: un "no rechaza" con n pequeno significa "no se detecto",
    NO "no lo hay", y el numero lo hace explicito.
    """
    if n <= 2:
        return float("nan")
    # t critico bilateral por biseccion sobre la cola.
    lo, hi = 0.0, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if t_sf(mid, n - 2) > alfa:
            lo = mid
        else:
            hi = mid
    tc = 0.5 * (lo + hi)
    return float(tc / math.sqrt(tc * tc + (n - 2)))


# ==============================================================================
# §H — COMPUERTA 1
# ==============================================================================
def evaluar_compuerta(m: dict) -> dict:
    """Aplica los seis criterios del §4 del preregistro. NO se rebaja ninguno.

    Devuelve el detalle completo Y la primera que falla, para que el reporte
    pueda decir POR QUE cayo un par y no solo que cayo.
    """
    lotes = (NOCIONAL_MAX_ORDEN_USD / m["granularidad_usd"]
             if m.get("granularidad_usd", 0) > 0 else 0.0)
    criterios = [
        ("C1.1 sigma1_roll >= %.2f pb" % C1_1_SIGMA1_MIN_PB,
         _ge(m.get("sigma1_roll_pb"), C1_1_SIGMA1_MIN_PB), m.get("sigma1_roll_pb")),
        ("C1.2 horquilla <= %.2f pb" % C1_2_HORQUILLA_MAX_PB,
         _le(m.get("horquilla_pb"), C1_2_HORQUILLA_MAX_PB), m.get("horquilla_pb")),
        ("C1.3 granularidad <= %.0f USD y >= %d lotes" % (C1_3_GRANULARIDAD_MAX_USD, C1_3_LOTES_MIN),
         _le(m.get("granularidad_usd"), C1_3_GRANULARIDAD_MAX_USD) and lotes >= C1_3_LOTES_MIN,
         m.get("granularidad_usd")),
        ("C1.4 fraccion_cola <= %.2f" % C1_4_FRACCION_SALTO_MAX,
         _le(m.get("fraccion_cola"), C1_4_FRACCION_SALTO_MAX), m.get("fraccion_cola")),
        ("C1.5 curtosis <= %.0f" % C1_5_CURTOSIS_MAX,
         _le(m.get("curtosis"), C1_5_CURTOSIS_MAX), m.get("curtosis")),
        ("C1.6 cobertura >= %.2f" % C1_6_COBERTURA_MIN,
         _ge(m.get("cobertura"), C1_6_COBERTURA_MIN), m.get("cobertura")),
    ]
    fallos = [(nombre, valor) for nombre, ok, valor in criterios if not ok]
    # ⚠ MARCA DE MARGINALIDAD. `sigma_1` es una extrapolacion (ver
    # `ajustar_ley_sigma`) y su error tipico medido sobre verdad conocida es de
    # varios por ciento. Un par cuyo sigma_1 cae a menos de 2 errores tipicos del
    # umbral NO esta decidido por el dato: esta decidido por el ruido del estimador.
    # No cambia el veredicto —el umbral es el declarado y no se toca— pero se marca,
    # que es lo contrario de silenciarlo.
    s1, ee = m.get("sigma1_roll_pb"), m.get("sigma1_ee_rel")
    marginal = False
    if s1 is not None and np.isfinite(s1) and ee is not None and np.isfinite(ee):
        marginal = abs(s1 - C1_1_SIGMA1_MIN_PB) < 2.0 * ee * s1
    return {
        "pasa": len(fallos) == 0,
        "marginal_C1_1": bool(marginal),
        "criterios": [{"criterio": n, "ok": ok, "valor": v} for n, ok, v in criterios],
        "primera_falla": fallos[0][0] if fallos else None,
        "n_fallos": len(fallos),
        "lotes": lotes,
    }


def _ge(v, u):
    return v is not None and np.isfinite(v) and v >= u


def _le(v, u):
    return v is not None and np.isfinite(v) and v <= u


# ==============================================================================
# §I — CRIBADO DE UN PAR (sin red: recibe los datos ya descargados)
# ==============================================================================
def cribar_par(simbolo: str, klines: np.ndarray, horquilla_pb: float,
               granularidad_usd: float, semilla: int = 0) -> dict:
    """Todas las metricas de un par. Funcion PURA: no toca la red.

    Separada de la descarga a proposito, para que la suite pueda ejercitarla contra
    datos sinteticos con verdad conocida — que es la unica forma de saber que el
    estimador mide lo que dice medir.
    """
    lim = limpiar_klines(klines)
    k = lim["k"]
    if len(k) < 2000:
        return {"simbolo": simbolo, "error": "datos insuficientes (%d velas)" % len(k),
                "cobertura": lim["cobertura"]}

    log_p = np.log(k[:, 4])
    r1 = np.diff(log_p)

    filas = firma_sigma(log_p, paso_s=60.0)
    crudo = ajustar_ley_sigma(filas)
    roll = roll_desde_retornos(r1)
    corr = corregir_firma_por_roll(filas, roll["var_rebote_pb2"])
    ajustado = ajustar_ley_sigma(corr["filas"]) if corr["filas"] else dict(crudo)

    # Control barajado (§6.2): bajo barajado los retornos son iid -> H_p debe dar 0.5.
    rng = np.random.default_rng(semilla)
    r_baraj = rng.permutation(r1)
    log_p_baraj = np.concatenate([[log_p[0]], log_p[0] + np.cumsum(r_baraj)])
    baraj = ajustar_ley_sigma(firma_sigma(log_p_baraj, paso_s=60.0))
    sesgo_hp = abs(baraj["H_p"] - 0.5) if np.isfinite(baraj["H_p"]) else float("nan")
    control_ok = np.isfinite(sesgo_hp) and sesgo_hp <= TOLERANCIA_HP_BARAJADO

    X, y = rasgos_y_objetivo(k)
    pred = r2_predictivo_oos(X, y)
    salto = metricas_de_salto(r1)

    m = {
        "simbolo": simbolo,
        "n_velas": len(k),
        "descartadas": lim["descartadas"],
        "cobertura": lim["cobertura"],
        "sigma1_crudo_pb": crudo["sigma1_pb"],
        "H_p_crudo": crudo["H_p"],
        "sigma1_roll_pb": ajustado["sigma1_pb"],
        "H_p": ajustado["H_p"],
        "r2_ajuste": ajustado["r2_ajuste"],
        "sigma1_ee_rel": ajustado.get("sigma1_ee_rel", float("nan")),
        "H_p_ee": ajustado.get("H_p_ee", float("nan")),
        # sigma al horizonte de referencia: NO es extrapolacion, cae dentro del
        # rango ajustado, asi que es el numero robusto de los dos. Se reporta al
        # lado de sigma_1 para poder ver cuando discrepan.
        "sigma_Href_pb": (ajustado["sigma1_pb"] * H_REFERENCIA_S ** ajustado["H_p"]
                          if np.isfinite(ajustado["sigma1_pb"]) and np.isfinite(ajustado["H_p"])
                          else float("nan")),
        "delta_roll_pct": (100.0 * (crudo["sigma1_pb"] - ajustado["sigma1_pb"])
                           / crudo["sigma1_pb"]
                           if np.isfinite(crudo["sigma1_pb"]) and crudo["sigma1_pb"] > 0
                           else float("nan")),
        "s_eff_roll_pb": roll["s_eff_pb"],
        "rho1_1min": roll["rho1"],
        "horizontes_descartados": corr["descartados"],
        "H_p_barajado": baraj["H_p"],
        "control_barajado_ok": bool(control_ok),
        "horquilla_pb": horquilla_pb,
        "granularidad_usd": granularidad_usd,
        "r2_pred_oos": pred["r2_oos"],
        "n_test": pred.get("n_test", 0),
        "fraccion_cola": salto["fraccion_cola"],
        "curtosis": salto["curtosis"],
    }
    m["r2_requerido"] = r2_requerido(m["sigma1_roll_pb"], m["H_p"]) \
        if np.isfinite(m["sigma1_roll_pb"]) and np.isfinite(m["H_p"]) else float("nan")
    m["compuerta"] = evaluar_compuerta(m)
    return m


# ==============================================================================
# §J — DESCARGA (la unica parte que toca la red)
# ==============================================================================
BASE = "https://fapi.binance.com"
DIR_CACHE = "cribado_cache"

# Layout de `klines` que consume todo lo de arriba. Los indices son los de la
# respuesta cruda de Binance, mapeados a 8 columnas.
KLINE_COLS = ("t_apertura_ms", "apertura", "maximo", "minimo", "cierre",
              "volumen", "n_trades", "volumen_taker_compra")
_IDX_BINANCE = (0, 1, 2, 3, 4, 5, 8, 9)

# Token bucket del peso de la API (Sec. 8.2 del PDF). Futuros: 2400/min.
# `klines` con limit > 1000 pesa 10; `exchangeInfo` 1; `bookTicker` sin simbolo 5.
PESO_POR_MINUTO = 2400
_gasto = []


def _gastar(peso: int) -> None:
    """Token bucket real, no un `pass`.  (La v1.3 ya aprendio esta leccion.)"""
    ahora = time.time()
    _gasto[:] = [(t, w) for t, w in _gasto if ahora - t < 60.0]
    usado = sum(w for _, w in _gasto)
    if usado + peso > PESO_POR_MINUTO * 0.85:      # 85 % del limite, margen deliberado
        espera = 60.0 - (ahora - _gasto[0][0]) + 0.5
        if espera > 0:
            log("  [peso] %d/%d usado; esperando %.1f s" % (usado, PESO_POR_MINUTO, espera))
            time.sleep(espera)
        _gasto.clear()
    _gasto.append((time.time(), peso))


def log(msg: str) -> None:
    """Todo el diagnostico pasa por aqui y NUNCA propaga.

    Leccion de la v1.2: bajo `spawn` un `print` fallido (BrokenPipeError) puede
    matar un proceso en silencio. Aqui no hay procesos hijos, pero la consola de
    Windows es cp1252 y un caracter no ASCII revienta igual.
    """
    try:
        print(msg, flush=True)
    except Exception:
        pass


def _get(ruta: str, peso: int, reintentos: int = 4):
    """GET con backoff exponencial (Sec. 8.4.2). Reintenta red, NO politica."""
    espera = 2.0
    for intento in range(reintentos):
        _gastar(peso)
        try:
            req = urllib.request.Request(BASE + ruta,
                                         headers={"User-Agent": "cancion-micelio/4.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (418, 429):               # rate limit: esperar de verdad
                log("  [%d] limitado por el exchange; esperando 60 s" % e.code)
                time.sleep(60.0)
                continue
            if e.code == 403:
                raise RuntimeError(
                    "403 del proxy/exchange en %s. Si corre en una sesion remota, el "
                    "egreso a Binance puede estar denegado por politica: correr en local."
                    % ruta) from e
            raise
        except Exception as e:
            if intento == reintentos - 1:
                raise
            log("  reintento %d en %s: %s" % (intento + 1, ruta, type(e).__name__))
            time.sleep(espera)
            espera *= 2
    raise RuntimeError("agotados los reintentos en " + ruta)


def descargar_klines(simbolo: str, dias: int = 90, usar_cache: bool = True) -> np.ndarray:
    """90 dias de velas de 1 minuto, paginadas, con cache en disco.

    El cache hace la corrida REANUDABLE: son ~87 peticiones por par y ~30 pares, y
    perder todo por un corte a mitad ya paso en este proyecto (la captura de 48 h de
    la v2.2 corrio media hora sin guardar un solo bloque).
    """
    os.makedirs(DIR_CACHE, exist_ok=True)
    ruta = os.path.join(DIR_CACHE, "%s_%dd.npz" % (simbolo, dias))
    if usar_cache and os.path.exists(ruta):
        with np.load(ruta) as z:
            return z["k"]

    fin_ms = int(time.time() * 1000)
    inicio_ms = fin_ms - dias * 86400 * 1000
    trozos, cursor = [], inicio_ms
    while cursor < fin_ms:
        datos = _get("/fapi/v1/klines?symbol=%s&interval=1m&startTime=%d&limit=1500"
                     % (simbolo, cursor), peso=10)
        if not datos:
            break
        arr = np.array([[float(f[i]) for i in _IDX_BINANCE] for f in datos], dtype=float)
        trozos.append(arr)
        nuevo = int(datos[-1][0]) + 60000
        if nuevo <= cursor:
            break
        cursor = nuevo
        if len(datos) < 1500:
            break
    if not trozos:
        return np.zeros((0, 8))
    k = np.vstack(trozos)
    _, unicos = np.unique(k[:, 0], return_index=True)     # deduplicar por t_apertura
    k = k[np.sort(unicos)]
    # Escritura atomica: descriptor abierto, no ruta. `savez_compressed` ANADE
    # ".npz" al nombre si se le pasa una ruta, y eso ya rompio una captura entera
    # en la v2.2 (`bloque.npz.tmp.npz` y `os.replace` con WinError 2).
    tmp = ruta + ".tmp"
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, k=k)
    os.replace(tmp, ruta)
    return k


def pares_por_volumen(n: int = 30, excluir_sufijos=("UP", "DOWN", "BULL", "BEAR")) -> list:
    """Los `n` pares USDⓈ-M perpetuos con mas volumen nocional en 24 h."""
    info = _get("/fapi/v1/exchangeInfo", peso=1)
    vivos = {}
    for s in info["symbols"]:
        if (s.get("status") == "TRADING" and s.get("quoteAsset") == "USDT"
                and s.get("contractType") == "PERPETUAL"
                and not any(s["symbol"].startswith(b + "UP") or s["symbol"].startswith(b + "DOWN")
                            for b in ("",))
                and not any(s["baseAsset"].endswith(x) for x in excluir_sufijos)):
            vivos[s["symbol"]] = s
    tick = _get("/fapi/v1/ticker/24hr", peso=40)
    filas = [(float(t["quoteVolume"]), t["symbol"]) for t in tick
             if t["symbol"] in vivos]
    filas.sort(reverse=True)
    return [s for _, s in filas[:n]], vivos


def filtros_del_par(spec: dict) -> dict:
    """`minQty`, `stepSize`, `tickSize`, `minNotional` desde `exchangeInfo`.

    ⚠ La v1.3 midio que Testnet es 10x MAS FINO que Mainnet y que calibrar contra
    Testnet produce un sistema que se degrada a interruptor en produccion. Aqui se
    lee siempre de Mainnet.
    """
    f = {x["filterType"]: x for x in spec.get("filters", [])}
    return {
        "min_qty": float(f.get("LOT_SIZE", {}).get("minQty", 0.0)),
        "step_size": float(f.get("LOT_SIZE", {}).get("stepSize", 0.0)),
        "tick_size": float(f.get("PRICE_FILTER", {}).get("tickSize", 0.0)),
        "min_notional": float(f.get("MIN_NOTIONAL", {}).get("notional", 0.0)),
    }


def muestrear_horquillas(simbolos: list, n_muestras: int = 60,
                         intervalo_s: float = 5.0) -> dict:
    """Horquilla relativa MEDIANA por par, muestreando `bookTicker` de todos a la vez.

    Una sola lectura no sirve: la horquilla se ensancha y se estrecha con la
    actividad, y lo que decide la compuerta C1.2 es la mediana, no una foto.
    Sin `symbol` la llamada devuelve TODOS los pares y pesa 5 — por eso se muestrea
    el conjunto entero en cada pasada en vez de par por par.
    """
    acum = {s: [] for s in simbolos}
    objetivo = set(simbolos)
    for i in range(n_muestras):
        try:
            datos = _get("/fapi/v1/ticker/bookTicker", peso=5)
        except Exception as e:
            log("  horquilla: pasada %d fallo (%s)" % (i, type(e).__name__))
            time.sleep(intervalo_s)
            continue
        for d in datos:
            s = d.get("symbol")
            if s not in objetivo:
                continue
            bid, ask = float(d.get("bidPrice", 0)), float(d.get("askPrice", 0))
            if bid > 0 and ask > bid:
                acum[s].append((ask - bid) / (0.5 * (ask + bid)) * 1e4)   # pb
        if i < n_muestras - 1:
            time.sleep(intervalo_s)
    return {s: (float(np.median(v)) if v else float("nan")) for s, v in acum.items()}


# ==============================================================================
# §K — REPORTE
# ==============================================================================
def _fmt(v, ancho=9, dec=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "n/d".rjust(ancho)
    return ("%*.*f" % (ancho, dec, v))


def reporte(resultados: list) -> str:
    """La TABLA COMPLETA, ordenada por sigma1. No los tres primeros.  §5."""
    out = []
    A = out.append
    buenos = [r for r in resultados if "error" not in r]
    malos = [r for r in resultados if "error" in r]
    buenos.sort(key=lambda r: (-r["sigma1_roll_pb"]
                               if np.isfinite(r["sigma1_roll_pb"]) else 0.0))

    A("=" * 126)
    A("CRIBADO DE ACTIVOS v4.0 - ETAPA 1 - TABLA COMPLETA (%d pares)" % len(buenos))
    A("Umbrales del PREREGISTRO_CRIBADO_4_0 Sec. 4, declarados antes de medir.")
    A("Comisiones %s (maker %.3f pb/lado)."
      % ("LEIDAS" if COMISIONES_LEIDAS_DE_LA_CUENTA else "ASUMIDAS VIP 0",
         COMISION_MAKER_PB))
    A("=" * 126)
    A("")
    cab = ("%-13s %9s %7s %9s %8s %7s %8s %9s %8s %8s %7s  %s"
           % ("par", "sig1_roll", "+-%", "sig1_crudo", "d_roll%", "H_p", "horq_pb",
              "R2_req", "R2_pred", "fr_cola", "curt", "compuerta"))
    A(cab)
    A("-" * 126)
    for r in buenos:
        c = r["compuerta"]
        veredicto = "PASA" if c["pasa"] else "cae: " + c["primera_falla"].split()[0]
        if c.get("marginal_C1_1"):
            veredicto += "  [MARGINAL en C1.1]"
        A("%-13s %s %s %s %s %s %s %s %s %s %s  %s"
          % (r["simbolo"],
             _fmt(r["sigma1_roll_pb"], 9, 4),
             _fmt(100.0 * r.get("sigma1_ee_rel", float("nan")), 7, 1),
             _fmt(r["sigma1_crudo_pb"], 9, 4),
             _fmt(r["delta_roll_pct"], 8, 1), _fmt(r["H_p"], 7, 3),
             _fmt(r["horquilla_pb"], 8, 3), _fmt(r["r2_requerido"], 9, 5),
             _fmt(r["r2_pred_oos"], 8, 5), _fmt(r["fraccion_cola"], 8, 3),
             _fmt(r["curtosis"], 7, 1), veredicto))
    A("-" * 126)
    A("'+-%' es el error tipico RELATIVO de sigma1, propagado del ajuste. Es COTA")
    A("INFERIOR: sigma1 se extrapola a H=1 s desde un ajuste en [60,3600] s y los")
    A("residuos de horizontes anidados estan correlacionados. Ver ajustar_ley_sigma.")
    if malos:
        A("")
        A("SIN MEDIR (%d):" % len(malos))
        for r in malos:
            A("  %-13s %s" % (r["simbolo"], r["error"]))

    # --- Control barajado (§6.2) ------------------------------------------------
    A("")
    A("CONTROL BARAJADO (Sec. 6.2): bajo barajado H_p debe dar 0.500 +- %.3f"
      % TOLERANCIA_HP_BARAJADO)
    sesgados = [r for r in buenos if not r["control_barajado_ok"]]
    if sesgados:
        A("  !! %d par(es) con SESGO PROPIO del estimador; su H_p NO es una medida:"
          % len(sesgados))
        for r in sesgados[:10]:
            A("    %-13s H_p_barajado = %s" % (r["simbolo"], _fmt(r["H_p_barajado"], 7, 4)))
    else:
        A("  OK: los %d pares medidos pasan el control." % len(buenos))

    # --- La prediccion falsable (§2) --------------------------------------------
    A("")
    A("=" * 126)
    A("PREDICCION FALSABLE (PREREGISTRO Sec. 2), declarada antes de medir:")
    A("  'sigma1 y R2 predictivo NO estaran negativamente correlacionados entre activos'")
    A("  Muere si  rho_Spearman <= %.2f  con  p < %.2f."
      % (SPEARMAN_MATA_HIPOTESIS, SPEARMAN_ALFA))
    sp = spearman([r["sigma1_roll_pb"] for r in buenos],
                  [r["r2_pred_oos"] for r in buenos])
    A("")
    A("  rho_Spearman = %s   p = %s   n = %d   |rho| minimo detectable = %s"
      % (_fmt(sp["rho"], 7, 4), _fmt(sp["p"], 7, 4), sp["n"],
         _fmt(sp.get("rho_minimo_detectable"), 6, 3)))
    if np.isfinite(sp["rho"]) and sp["rho"] <= SPEARMAN_MATA_HIPOTESIS and sp["p"] < SPEARMAN_ALFA:
        A("  VEREDICTO: LA HIPOTESIS MUERE EN LA ETAPA 1. No se captura nada.")
    elif not np.isfinite(sp["rho"]):
        A("  VEREDICTO: SIN DATOS SUFICIENTES para evaluar la prediccion.")
    else:
        A("  VEREDICTO: la hipotesis SOBREVIVE (no se confirma: sobrevive).")
        A("  !! Con n = %d, un 'no rechaza' significa NO SE DETECTO, no 'no lo hay'."
          % sp["n"])

    # --- Compuerta 1 -------------------------------------------------------------
    pasan = [r for r in buenos if r["compuerta"]["pasa"]]
    A("")
    A("COMPUERTA 1: %d de %d pares la superan." % (len(pasan), len(buenos)))
    if pasan:
        g = max(pasan, key=lambda r: r["sigma1_roll_pb"])
        A("  Pasan: " + ", ".join(r["simbolo"] for r in pasan))
        A("  GANADOR para la etapa 2 (mayor sigma1 corregido): %s" % g["simbolo"])
        A("  Regla de parada Sec. 5: NO se cambia de par a mitad de la etapa 2.")
    else:
        A("  Ningun par la supera. La etapa 2 NO se ejecuta.")
        A("  Sec. 8: no se rebaja ningun umbral declarado en el Sec. 4.")
    # Motivos de caida, agregados: dice DONDE esta el cuello de botella.
    motivos = {}
    for r in buenos:
        for c in r["compuerta"]["criterios"]:
            if not c["ok"]:
                motivos[c["criterio"]] = motivos.get(c["criterio"], 0) + 1
    if motivos:
        A("")
        A("  Motivos de caida (un par puede fallar varios):")
        for crit, n in sorted(motivos.items(), key=lambda kv: -kv[1]):
            A("    %3d x  %s" % (n, crit))
    A("=" * 126)
    return "\n".join(out)


# ==============================================================================
# §L — MAIN
# ==============================================================================
def main(argv) -> int:
    n_pares = 30
    dias = 90
    muestras_horquilla = 60
    sin_cache = "--sin-cache" in argv
    for a in argv:
        if a.startswith("--pares="):
            n_pares = int(a.split("=")[1])
        elif a.startswith("--dias="):
            dias = int(a.split("=")[1])
        elif a.startswith("--horquilla="):
            muestras_horquilla = int(a.split("=")[1])

    log("== CRIBADO v4.0 - ETAPA 1 ==")
    log("pares=%d  dias=%d  muestras de horquilla=%d" % (n_pares, dias, muestras_horquilla))
    log("Umbrales: C1.1 sigma1 >= %.2f pb | C1.2 horquilla <= %.2f pb | "
        "C1.4 cola <= %.2f | C1.5 curtosis <= %.0f"
        % (C1_1_SIGMA1_MIN_PB, C1_2_HORQUILLA_MAX_PB,
           C1_4_FRACCION_SALTO_MAX, C1_5_CURTOSIS_MAX))
    log("")

    try:
        simbolos, specs = pares_por_volumen(n_pares)
    except RuntimeError as e:
        log("FALLO DE RED: %s" % e)
        return 2
    log("Pares por volumen 24 h: %s" % ", ".join(simbolos))

    log("")
    log("Muestreando horquillas (%d pasadas)..." % muestras_horquilla)
    horquillas = muestrear_horquillas(simbolos, muestras_horquilla)

    resultados = []
    for i, s in enumerate(simbolos, 1):
        log("[%2d/%2d] %s" % (i, len(simbolos), s))
        try:
            k = descargar_klines(s, dias, usar_cache=not sin_cache)
        except Exception as e:
            resultados.append({"simbolo": s, "error": "descarga: %s" % type(e).__name__})
            continue
        f = filtros_del_par(specs[s])
        precio = float(k[-1, 4]) if len(k) else 0.0
        gran = f["min_qty"] * precio
        # `minNotional` puede atar mas que `minQty`: manda el mayor de los dos.
        gran = max(gran, f["min_notional"])
        resultados.append(cribar_par(s, k, horquillas.get(s, float("nan")), gran))

    txt = reporte(resultados)
    log("")
    log(txt)
    with open("cribado_etapa1.txt", "w", encoding="ascii", errors="replace") as fh:
        fh.write(txt + "\n")
    log("")
    log("Reporte escrito en cribado_etapa1.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
