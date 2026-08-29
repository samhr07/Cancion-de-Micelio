"""
Proyecto: Canción del Micelio — Identidad de covarianza multivariante (§C.3.5)
Módulo: identidad_multivariante.py — la tarea residual de BTC.

⚠ NO SE IMPORTA DESDE `Micelio.py`.

POR QUE ESTA TAREA VA PRIMERA
-----------------------------
Es la única cosa del proyecto de BTC que quedó **sin contrastar**, es barata, y su
resultado **fija el `R²` alcanzable que entra en el factor 3.4×** del cribado. Si la
multivariante rinde el doble que la univariante, el factor necesario baja a 2.4× y
todo el cribado se vuelve más fácil (§ del preregistro: el margen va con `√R²`, así
que duplicar `R²` divide el factor de `σ` requerido por `√2 = 1.414`, y
3.4 / 1.414 = 2.40).

LA IDENTIDAD
------------
Para el mejor predictor lineal de `y` sobre los rasgos `X`:

    R2_max = c' * Sigma^-1 * c / var(y)        c = Cov(X, y),  Sigma = Cov(X, X)

Es un **techo**: ningún predictor lineal sobre esos rasgos lo supera. Y es
**superior por construcción** a cualquier `R²` univariante, porque el univariante es
el mismo cociente restringido a una coordenada — está anidado.

⚠ EL TECHO ESTA SESGADO AL ALZA EN MUESTRA FINITA, y el sesgo es del orden de `p/n`.
Por eso la identidad **no se lee sola**: se lee contra un **suelo** medido, y ese
suelo es el que decide si el número es falsable.

POR QUE ROTACION CIRCULAR Y NO BARAJADO
---------------------------------------
El barajado destruye la autocorrelación de las series, y estos rasgos son series de
tiempo fuertemente autocorreladas. Un suelo por barajado saldría demasiado bajo y
haría falsable cualquier cosa. La **rotación circular** desplaza `y` respecto de `X`
conservando **exactamente** la marginal Y la estructura de autocorrelación de ambas,
y destruyendo solo el alineamiento entre ellas — que es justo la hipótesis a probar.
Es el mismo principio que el nulo por sustitutos de la v2.2, corregido para series
con memoria.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np


# La marca falsable/no-falsable (corrección del 2026-08-29).
PERCENTIL_SUELO = 95.0
N_ROTACIONES = 500


# ==============================================================================
# §1 — LA IDENTIDAD
# ==============================================================================
def identidad_covarianza(X: np.ndarray, y: np.ndarray, ridge: float = 0.0) -> dict:
    """`R²_max = c' Σ⁻¹ c / var(y)`, con el techo por fila y el diagnostico numerico.

    ⚠ `Σ⁻¹` NUNCA por `np.linalg.inv`. Es la nota de "menores" del diagnostico
    original del proyecto y aqui muerde de verdad: con rasgos casi colineales
    (volumen y n_trades, por ejemplo) `Σ` queda mal condicionada y la inversa
    explota el `R²` hacia arriba sin avisar. Se resuelve por sistema lineal y se
    publica el numero de condicion y el rango efectivo.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    X, y = X[ok], y[ok]
    n, p = X.shape
    if n <= p + 2:
        return {"r2_max": float("nan"), "n": n, "p": p,
                "aviso": "n <= p + 2: la identidad no esta definida"}

    # ⚠ Un rasgo CONSTANTE (varianza cero) hace singular a `Sigma` y su contribucion
    # al techo no significa nada: no puede explicar variacion porque no varia. Es un
    # caso degenerado distinto de la colinealidad genuina, y se resuelve quitandolo,
    # no avisando y siguiendo. Se publica cual se quito.
    var_col = X.var(axis=0)
    escala = float(np.max(var_col)) if len(var_col) else 0.0
    vivos = var_col > escala * 1e-12 if escala > 0 else np.zeros(p, dtype=bool)
    descartados = [int(i) for i in np.flatnonzero(~vivos)]
    if descartados:
        X = X[:, vivos]
        p = X.shape[1]
        if p == 0:
            return {"r2_max": float("nan"), "n": n, "p": 0,
                    "aviso": "todos los rasgos son constantes"}

    Xc = X - X.mean(axis=0)
    yc = y - y.mean()
    var_y = float(yc @ yc / (n - 1))
    if var_y <= 0.0:
        return {"r2_max": float("nan"), "n": n, "p": p, "aviso": "var(y) = 0"}

    Sigma = (Xc.T @ Xc) / (n - 1)
    c = (Xc.T @ yc) / (n - 1)
    if ridge > 0.0:
        Sigma = Sigma + ridge * np.trace(Sigma) / p * np.eye(p)

    evals = np.linalg.eigvalsh(Sigma)
    cond = float(evals[-1] / evals[0]) if evals[0] > 0 else float("inf")
    rango_efectivo = int(np.sum(evals > evals[-1] * 1e-10))

    try:
        beta = np.linalg.solve(Sigma, c)
    except np.linalg.LinAlgError:
        beta = np.linalg.lstsq(Sigma, c, rcond=None)[0]
    r2_max = float(c @ beta / var_y)

    # Techo POR FILA: el R2 univariante de cada rasgo por separado.
    diag = np.diag(Sigma)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2_fila = np.where(diag > 0, c ** 2 / (diag * var_y), np.nan)

    if descartados:
        completo = np.full(len(vivos), np.nan)
        completo[vivos] = r2_fila
        r2_fila = completo
    return {
        "r2_max": r2_max,
        "r2_por_fila": r2_fila,
        "rasgos_constantes": descartados,
        "r2_univariante_max": float(np.nanmax(r2_fila)) if p else float("nan"),
        "suma_univariantes": float(np.nansum(r2_fila)),
        "beta": beta,
        "n": int(n), "p": int(p),
        "cond_Sigma": cond,
        "rango_efectivo": rango_efectivo,
        "sesgo_esperado_aprox": float(p / n),   # sesgo al alza de un R2 en muestra
    }


# ==============================================================================
# §2 — EL SUELO POR ROTACION CIRCULAR
# ==============================================================================
def suelo_por_rotacion(X: np.ndarray, y: np.ndarray, n_rot: int = N_ROTACIONES,
                       semilla: int = 0, margen: float = 0.05) -> dict:
    """Distribucion de `R²_max` bajo desalineamiento, conservando la memoria.

    Se rota `y` respecto de `X` por desplazamientos aleatorios, evitando los
    desplazamientos casi nulos (`margen`) que dejarian las series todavia alineadas.
    Devuelve el percentil que actua de suelo.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    X, y = X[ok], y[ok]
    n = len(y)
    if n < 50:
        return {"suelo": float("nan"), "n_rot": 0}
    rng = np.random.default_rng(semilla)
    lo = max(1, int(n * margen))
    valores = []
    for _ in range(n_rot):
        k = int(rng.integers(lo, n - lo)) if n - 2 * lo > 1 else 1
        res = identidad_covarianza(X, np.roll(y, k))
        if np.isfinite(res.get("r2_max", np.nan)):
            valores.append(res["r2_max"])
    if not valores:
        return {"suelo": float("nan"), "n_rot": 0}
    v = np.array(valores)
    return {
        "suelo": float(np.percentile(v, PERCENTIL_SUELO)),
        "mediana": float(np.median(v)),
        "p05": float(np.percentile(v, 5.0)),
        "max": float(v.max()),
        "n_rot": len(v),
        "percentil": PERCENTIL_SUELO,
    }


def marcar_falsable(r2_max: float, suelo: dict) -> dict:
    """La marca falsable / NO falsable (correccion del 2026-08-29).

    ⚠ NO es un p-valor y no se debe leer como tal. Dice si el techo medido se
    distingue del techo que sale por puro tamano de muestra. Un `R²` por debajo del
    suelo **no es un `R²` pequeno: es un `R²` que no significa nada**, y presentarlo
    como una medida seria exactamente el error que esta marca existe para impedir.
    """
    s = suelo.get("suelo", float("nan"))
    if not (np.isfinite(r2_max) and np.isfinite(s)):
        return {"falsable": False, "motivo": "sin suelo medible", "exceso": float("nan")}
    if r2_max <= s:
        return {"falsable": False,
                "motivo": "R2_max (%.6f) <= suelo p%.0f (%.6f): indistinguible del "
                          "desalineamiento" % (r2_max, suelo["percentil"], s),
                "exceso": float(r2_max - s)}
    return {"falsable": True,
            "motivo": "R2_max (%.6f) supera el suelo p%.0f (%.6f) por %.1fx"
                      % (r2_max, suelo["percentil"], s, r2_max / s if s > 0 else float("inf")),
            "exceso": float(r2_max - s),
            "razon": float(r2_max / s) if s > 0 else float("inf")}


# ==============================================================================
# §3 — LA CONSECUENCIA ECONOMICA: que le hace al factor 3.4x
# ==============================================================================
def factor_sigma_requerido(r2_referencia: float, r2_nuevo: float,
                           factor_actual: float) -> float:
    """Como cambia el factor de `σ` necesario al mejorar el `R²`.

    El margen va con `ρ = √R²`, luego el factor de `σ` requerido va con `1/√R²`:

        factor_nuevo = factor_actual / sqrt(r2_nuevo / r2_referencia)

    Con `r2_nuevo = 2 * r2_referencia`:  3.4 / sqrt(2) = **2.40**, que es
    exactamente la cifra que la orden de trabajo anticipa.
    """
    if r2_referencia <= 0 or r2_nuevo <= 0:
        return float("nan")
    return factor_actual / math.sqrt(r2_nuevo / r2_referencia)


# ==============================================================================
# §4 — CARGA DE DATOS Y ADAPTADOR DE RASGOS
# ==============================================================================
def cargar_captura(ruta: str) -> dict:
    """Carga una captura `.npz` (formato de `captura_larga.py`) o un directorio de bloques."""
    if os.path.isdir(ruta):
        bloques = sorted(f for f in os.listdir(ruta) if f.endswith(".npz"))
        if not bloques:
            raise FileNotFoundError("no hay bloques .npz en " + ruta)
        acum = {}
        for b in bloques:
            with np.load(os.path.join(ruta, b)) as z:
                for k in z.files:
                    acum.setdefault(k, []).append(z[k])
        return {k: np.concatenate(v) for k, v in acum.items()}
    with np.load(ruta) as z:
        return {k: z[k] for k in z.files}


def rasgos_desde_captura(d: dict, ventana_tx: int = 200) -> tuple:
    """Construye `X` e `y` por bloques NO SOLAPADOS de `ventana_tx` transacciones.

    # NOTA DE INTERPRETACION: el §C.3.5 y el conjunto de rasgos «estacionales» que
    # midieron `R² = 0.615` en volatilidad NO estan en este arbol (ver el §7.2 del
    # PREREGISTRO_CRIBADO_4_0). Este adaptador reconstruye la familia de rasgos que
    # SI esta documentada en el proyecto —desequilibrio de signo del propagador
    # (v3.1 §2), volumen, tasa de transacciones y volatilidad realizada— para que el
    # modulo sea ejecutable de punta a punta. **No pretende ser el conjunto de
    # rasgos del §C.3.5**: al correr sobre los datos reales hay que sustituirlo por
    # el suyo, y la identidad del §1 no cambia al hacerlo.

    Objetivo por omision: la **volatilidad realizada del bloque siguiente**, que es
    el observable con `R² = 0.615` y no el direccional. Se devuelven ambos.
    """
    p = np.asarray(d["tr_precio"], dtype=float)
    t = np.asarray(d["tr_t"], dtype=float)
    q = np.asarray(d.get("tr_cant", np.ones_like(p)), dtype=float)
    m = np.asarray(d.get("tr_maker", np.zeros_like(p)), dtype=float)

    # ⚠ Los ceros del feed (v3.0): ~0.2 % de trades con p = 0. Fuera ANTES de medir
    # varianza; recordar que limpiar SUBE la curtosis, no la baja.
    val = np.isfinite(p) & (p > 0) & np.isfinite(t)
    p, t, q, m = p[val], t[val], q[val], m[val]
    if len(p) < 10 * ventana_tx:
        return np.zeros((0, 6)), np.zeros(0), np.zeros(0)

    eps = np.where(m > 0.5, -1.0, 1.0)          # v3.1 §2.2: m=True -> el taker vendia
    r = np.diff(np.log(p))
    nb = (len(p) - 1) // ventana_tx
    corte = nb * ventana_tx
    rs = lambda a: a[:corte].reshape(nb, ventana_tx)
    rb, qb, eb = rs(r), rs(q[1:]), rs(eps[1:])
    tb = rs(t[1:])

    ret = rb.sum(axis=1)
    vol_real = np.sqrt((rb ** 2).sum(axis=1))                    # volatilidad realizada
    dur = np.maximum(tb[:, -1] - tb[:, 0], 1e-9)
    X = np.column_stack([
        (eb * qb).sum(axis=1) / np.maximum(qb.sum(axis=1), 1e-12),   # desequilibrio firmado
        eb.mean(axis=1),                                             # desequilibrio en n
        np.log(np.maximum(qb.sum(axis=1), 1e-12)),                   # volumen
        np.log(ventana_tx / dur),                                    # nu, tasa de transacciones
        vol_real,                                                    # volatilidad del bloque
        ret,                                                         # retorno del bloque
    ])
    # Objetivos: bloque SIGUIENTE. NaN al final para que la mascara lo descarte.
    y_vol = np.concatenate([vol_real[1:], [np.nan]])
    y_dir = np.concatenate([ret[1:], [np.nan]])
    return X, y_vol, y_dir


# ==============================================================================
# §5 — REPORTE
# ==============================================================================
def analizar(X, y, etiqueta: str, semilla: int = 0, n_rot: int = N_ROTACIONES) -> dict:
    res = identidad_covarianza(X, y)
    suelo = suelo_por_rotacion(X, y, n_rot=n_rot, semilla=semilla)
    marca = marcar_falsable(res.get("r2_max", float("nan")), suelo)
    return {"etiqueta": etiqueta, "identidad": res, "suelo": suelo, "marca": marca}


def reporte(analisis: dict, nombres_rasgos=None) -> str:
    r, s, mk = analisis["identidad"], analisis["suelo"], analisis["marca"]
    out = ["=" * 92,
           "IDENTIDAD DE COVARIANZA MULTIVARIANTE (Sec. C.3.5) - %s" % analisis["etiqueta"],
           "=" * 92]
    A = out.append
    if not np.isfinite(r.get("r2_max", float("nan"))):
        A("SIN RESULTADO: %s" % r.get("aviso", "desconocido"))
        A("=" * 92)
        return "\n".join(out)

    A("n = %d observaciones, p = %d rasgos" % (r["n"], r["p"]))
    A("cond(Sigma) = %.3e, rango efectivo %d/%d" % (r["cond_Sigma"], r["rango_efectivo"], r["p"]))
    if r.get("rasgos_constantes"):
        A("  !! rasgos CONSTANTES descartados (varianza cero): %s"
          % ", ".join(str(i) for i in r["rasgos_constantes"]))
        A("     Un rasgo que no varia no puede explicar variacion; se quita, no se avisa.")
    if r["rango_efectivo"] < r["p"]:
        A("  !! Sigma es DEFICIENTE EN RANGO: hay rasgos colineales y el techo esta")
        A("     inflado. Quitar rasgos redundantes antes de leer el numero.")
    A("")
    A("TECHO POR FILA (R2 univariante de cada rasgo):")
    for i, v in enumerate(r["r2_por_fila"]):
        nom = nombres_rasgos[i] if nombres_rasgos and i < len(nombres_rasgos) else "rasgo_%d" % i
        A("  %-28s R2 = %.6f" % (nom, v))
    A("")
    A("  max univariante   = %.6f" % r["r2_univariante_max"])
    A("  suma univariantes = %.6f  (no es un techo: los rasgos estan correlacionados)"
      % r["suma_univariantes"])
    A("  R2_MAX MULTIVARIANTE = %.6f" % r["r2_max"])
    ganancia = (r["r2_max"] / r["r2_univariante_max"]
                if r["r2_univariante_max"] > 0 else float("nan"))
    A("  ganancia sobre el mejor univariante = %.3fx" % ganancia)
    A("  sesgo al alza esperado por muestra finita ~ p/n = %.6f" % r["sesgo_esperado_aprox"])
    A("")
    A("SUELO POR ROTACION CIRCULAR (%d rotaciones):" % s.get("n_rot", 0))
    A("  p05 = %.6f | mediana = %.6f | p%.0f = %.6f | max = %.6f"
      % (s.get("p05", float("nan")), s.get("mediana", float("nan")),
         s.get("percentil", 95), s.get("suelo", float("nan")), s.get("max", float("nan"))))
    A("")
    A("MARCA: %s" % ("FALSABLE" if mk["falsable"] else "NO FALSABLE"))
    A("  %s" % mk["motivo"])
    if not mk["falsable"]:
        A("  Un R2 por debajo del suelo NO es un R2 pequeno: es un R2 que no significa")
        A("  nada. No entra en el criterio economico ni en el factor de sigma.")
    A("=" * 92)
    return "\n".join(out)


NOMBRES_RASGOS = ("desequilibrio firmado por volumen", "desequilibrio en numero",
                  "log volumen", "log nu (tasa tx)", "volatilidad realizada",
                  "retorno del bloque")


def main(argv) -> int:
    ruta = None
    ventana = 200
    n_rot = N_ROTACIONES
    for a in argv:
        if a.startswith("--datos="):
            ruta = a.split("=", 1)[1]
        elif a.startswith("--ventana="):
            ventana = int(a.split("=")[1])
        elif a.startswith("--rotaciones="):
            n_rot = int(a.split("=")[1])
    if not ruta:
        print("uso: python identidad_multivariante.py --datos=<captura.npz | dir_bloques>")
        print("     [--ventana=200] [--rotaciones=500]")
        print("")
        print("Los datos de BTC NO estan en este repositorio (ver el Sec. 7.2 del")
        print("PREREGISTRO_CRIBADO_4_0). Apuntar a la captura local.")
        return 1
    try:
        d = cargar_captura(ruta)
    except Exception as e:
        print("no se pudo cargar %s: %s" % (ruta, e))
        return 2

    X, y_vol, y_dir = rasgos_desde_captura(d, ventana)
    if len(y_vol) == 0:
        print("captura demasiado corta para ventana=%d" % ventana)
        return 3

    print("Captura: %d transacciones -> %d bloques de %d tx" % (len(d["tr_precio"]), len(y_vol), ventana))
    print("")
    for etiqueta, y in (("VOLATILIDAD del bloque siguiente", y_vol),
                        ("DIRECCION del bloque siguiente", y_dir)):
        an = analizar(X, y, etiqueta, n_rot=n_rot)
        print(reporte(an, NOMBRES_RASGOS))
        print("")
        r2 = an["identidad"].get("r2_max", float("nan"))
        uni = an["identidad"].get("r2_univariante_max", float("nan"))
        if an["marca"]["falsable"] and np.isfinite(uni) and uni > 0:
            for base, f in (("estacional_0", 3.4), ("estacional_1", 1.2)):
                nuevo = factor_sigma_requerido(uni, r2, f)
                print("  CONSECUENCIA: factor de sigma sobre %s: %.2fx -> %.2fx"
                      % (base, f, nuevo))
            print("")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
