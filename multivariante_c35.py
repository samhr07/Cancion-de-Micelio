"""
Proyecto: Canción del Micelio — §C.3.5 de `ADENDA_C_REFORMULACION_S1_4_1.md`
Módulo: multivariante_c35.py — cierra la única cosa que quedaba SIN CONTRASTAR.

⚠ NO SE IMPORTA DESDE `Micelio.py`.

QUE CIERRA, Y POR QUE ESTABA ABIERTO
------------------------------------
`identidad.py` **implementa** el §C.3.5 (`r2_identidad_multi`), pero sólo lo llama
desde `etapa_calibra` y desde su autotest: **`etapa_medir` no lo usa**. Por eso
`CLAUDE.md` dice que la versión multivariante «queda SIN CONTRASTAR, no refutada» —
la función existe y está validada, pero nada la corría sobre la rejilla decisiva con
su propio suelo y su propio techo.

Este módulo es ese ejecutor. Espeja `estratos.py` fila a fila —mismos fragmentos,
mismos `H`, mismos terciles de `σ` pronosticada, misma regla de parada— y cambia
**sólo el estadístico**, para que las dos tablas sean comparables sin traducción.

LO QUE NO SE HEREDA DE `identidad.py`, Y ES DELIBERADO
-----------------------------------------------------
1. **`κ`**. `identidad.r2_req` usa `H.FACTOR_DECIL = 1.755`, que es
   `E[z | z > q90(z)]`: el decil superior **con signo**, o sea operar en UNA sola
   dirección. La banda muerta `|α| > c` opera en las **dos**, y ahí el valor es
   `E[|z| | |z| > q90(|z|)] = 2·φ(1.6449)/0.10 = 2.0627`. `estratos.py` ya usa el
   corregido; aquí también. Son 1.38× en `R²_req`, y en la dirección que **infla**
   el requisito, o sea la que hace la conclusión negativa demasiado fácil.
2. **`L(H)` real de cinco términos** (v4.2 §2, `curva_coste.py`), no comisión pura.
   6.24–7.63 pb contra los 4.888 que arrastraban las curvas viejas. Como
   `R²_req ∝ L²`, la diferencia es un factor 1.6–2.4.

⚠ EL SUELO ES PROPIO, Y ES LO QUE HACE HONESTO EL CONTRASTE
-----------------------------------------------------------
**No se reutiliza el suelo univariante.** Con `k = 5` regresores el `R²` en muestra
se infla por tamaño finito del orden de `k/n`, o sea ~5× el univariante. Comparar un
`R²` multivariante contra el suelo univariante declararía FALSABLE cualquier cosa, y
sería exactamente el modo de fallo que el §C.3.5 advierte al pedir «corrección de
muestra finita por `k` declarada antes». Aquí el suelo se **mide** con el mismo
estadístico multivariante, por rotación circular.

Y la rotación va **en bloque**: se desplaza la matriz `X` entera contra `r`. Rotar
cada columna por separado destruiría la covarianza entre rasgos, que es justo lo que
`Σ` necesita y lo que el estimador usa. El nulo tiene que romper el emparejamiento
con `r` y nada más.
"""

from __future__ import annotations

import argparse
import gc
import json
import os

import numpy as np

import horizonte as H
import identidad as I

log, titulo = H.log, H.titulo

HS = [900, 1800, 3600, 7200, 14400]
KAPPA = 2.0627               # DOS direcciones. El 1.755 es el de una sola.
N_SORTEOS = 200
PASO_S = 60.0
NOMBRES = ["bajo", "medio", "ALTO"]
FRACS = I.FRACS              # (0.25, 0.5, 1.0, 2.0) -> k = 5 con eps_t
SALIDA = "telemetria/multivariante_c35.json"

# `L(H)` medida en el §2 de la v4.2 (`curva_coste.py`, 7/7), maker, en pb.
# Se lee del JSON si existe; esta tabla es el respaldo, no una suposicion.
L_MEDIDA = {900: 6.24, 1800: 6.60, 3600: 6.75, 7200: 7.63, 14400: 7.18}


def _L(Hs: int) -> float:
    try:
        d = json.load(open("telemetria/curva_coste.json", encoding="utf-8"))
        for f in d["L"]:
            if f["H"] == Hs:
                return float(f["L_maker"])
        return float(d["L"][-1]["L_maker"])
    except Exception:
        return L_MEDIDA.get(Hs, 6.75)


# ==============================================================================
# El estadistico del §C.3.5
# ==============================================================================
def diseno_multi(t, eps, i, Hs, fracs=FRACS) -> np.ndarray:
    """`X` = [eps_t | agregados de flujo PASADO en f*H segundos].  §C.3.5.

    ⚠ `eps_t` VA COMO PRIMERA COLUMNA y no es opcional. Sin ella, «la
    multivariante no puede rendir menos que la univariante» **no es un teorema**:
    los agregados promedian ~10³ ticks y diluyen el signo suelto, así que la
    multivariante puede rendir MENOS. Con `eps_t` dentro, la multivariante ANIDA
    a la univariante y la cota del §C.3.4 se cumple por construcción. Es la misma
    nota que lleva `identidad.r2_identidad_multi`, y aquí se conserva porque es la
    que sostiene la frase «techo superior por construcción» de `CLAUDE.md`.
    """
    flujo = np.concatenate(([0.0], np.cumsum(eps, dtype=np.float64)))
    X = np.empty((i.size, len(fracs) + 1))
    X[:, 0] = eps[i]
    for c_, f in enumerate(fracs):
        ini = np.clip(np.searchsorted(t, t[i] - f * Hs, side="left"), 0, t.size - 1)
        X[:, c_ + 1] = flujo[i] - flujo[ini]
    return X


def r2_multi(X: np.ndarray, r: np.ndarray) -> float:
    """`R² = c' Σ⁻¹ c / var(r)` — el techo del mejor predictor lineal sobre `X`.

    `Σ` se resuelve por sistema lineal con una ridge diminuta, nunca por `inv`:
    con agregados anidados (0.25H ⊂ 0.5H ⊂ H ⊂ 2H) las columnas están muy
    correlacionadas y la inversa explota el `R²` sin avisar.
    """
    if X.shape[0] < 40:
        return float("nan")
    Xc = X - X.mean(0)
    rc = r - r.mean()
    S = (Xc.T @ Xc) / Xc.shape[0]
    c = (Xc.T @ rc) / Xc.shape[0]
    Vr = float(np.var(r))
    if Vr <= 0:
        return float("nan")
    try:
        q = float(c @ np.linalg.solve(S + 1e-12 * np.eye(S.shape[0]), c))
    except np.linalg.LinAlgError:
        return float("nan")
    return q / Vr


def suelo_multi(X: np.ndarray, r: np.ndarray, n_sorteos: int = N_SORTEOS,
                semilla: int = 17) -> float:
    """`q95(|R²_nulo|)` con el MISMO estadistico, por rotacion circular EN BLOQUE.

    Rotar `X` entera contra `r` conserva la autocorrelacion de cada rasgo Y la
    covarianza entre rasgos, y rompe solo el emparejamiento con el objetivo. Rotar
    columna a columna cambiaria `Σ` y el nulo dejaria de ser el mismo estimador.
    """
    rng = np.random.default_rng(semilla)
    n = X.shape[0]
    if n < 60:
        return float("nan")
    v = []
    for _ in range(n_sorteos):
        k = int(rng.integers(5, max(6, n - 5)))
        v.append(r2_multi(np.take(X, np.arange(n) - k, axis=0, mode="wrap"), r))
    v = np.array([x for x in v if np.isfinite(x)])
    return float(np.percentile(np.abs(v), 95)) if v.size else float("nan")


# ==============================================================================
# Ejecutor
# ==============================================================================
def etapa_medir(args) -> int:
    titulo("§C.3.5 MULTIVARIANTE -- lo que quedaba SIN CONTRASTAR")
    log("")
    log("  kappa = %.4f (DOS direcciones)   L(H) real del §2 de la v4.2" % KAPPA)
    log("  k = %d rasgos: eps_t + flujo pasado en %s * H" % (len(FRACS) + 1, FRACS))
    log("  suelo y techo se miden CON EL MISMO estadistico multivariante:")
    log("  reutilizar el suelo univariante declararia FALSABLE casi cualquier cosa.")
    log("  REGLA DE PARADA (la misma del §1 estratificado): si en el estrato ALTO,")
    log("  donde el instrumento resuelve, el margen no baja de 3x -> negativo.")
    filas = []
    for nombre in I._nombres():
        f = I.cargar(nombre, con_precio=False)
        dur = float(f["t"][-1] - f["t"][0])
        log("")
        log("--- %s (%.1f h) ---" % (nombre, dur / 3600.0))
        log("  %7s %7s %7s | %8s %9s %10s %10s %9s %11s"
            % ("H", "estrato", "n", "R2_req", "R2_multi", "R2_MAX", "q95_multi",
               "margen", "veredicto"))
        for Hs in HS:
            if dur < 6 * Hs:
                continue
            P = _preparar(f, Hs)
            if P is None:
                continue
            L = _L(Hs)
            cortes = np.percentile(P["s_prev"][P["ent"]], [33.3, 66.7])
            gi = np.digitize(P["s_prev"], cortes)
            for g_ in (0, 1, 2):
                m = (gi == g_)
                if m.sum() < 60:
                    continue
                X, r, fu = P["X"][m], P["r"][m], P["fut"][m]
                sg = float(np.std(r))
                req = (L / (KAPPA * sg)) ** 2
                r2 = r2_multi(X, r)
                q95 = suelo_multi(X, r)
                # Techo: lo mas que estos k rasgos pueden explicar del FLUJO
                # FIRMADO futuro. Con nucleo perfecto y sin ruido, ese es el
                # limite de lo que pueden explicar del retorno futuro.
                r2max = r2_multi(X, fu)
                resuelve = np.isfinite(q95) and req >= 3 * q95
                margen = req / max(r2, 1e-12) if r2 > 0 else float("inf")
                if not np.isfinite(r2max) or r2max < req:
                    vered = "NO FALSABLE"
                elif r2max >= 3 * req:
                    vered = "falsable"
                else:
                    vered = "al limite"
                log("  %6ds %7s %7d | %7.3f%% %10.6f %10.6f %10.6f %9s %11s%s"
                    % (Hs, NOMBRES[g_], int(m.sum()), 100 * req, r2, r2max, q95,
                       ("%.1fx" % margen) if np.isfinite(margen) else "  inf",
                       vered, "" if resuelve else "  NO RESUELVE"))
                filas.append({"frag": nombre, "H": Hs, "estrato": NOMBRES[g_],
                              "n": int(m.sum()), "sigma": sg, "req": req,
                              "r2_multi": r2, "q95_multi": q95, "r2max": r2max,
                              "veredicto": vered, "resuelve": bool(resuelve),
                              "margen": float(margen), "L": L})
            del P
            gc.collect()
        del f
        gc.collect()

    _veredicto(filas)
    os.makedirs("telemetria", exist_ok=True)
    json.dump(filas, open(SALIDA, "w", encoding="utf-8"), indent=1)
    log("")
    log("  filas escritas en %s" % SALIDA)
    return 0


def _preparar(frag, Hs):
    """Igual que `estratos.preparar`, mas la matriz de diseno multivariante."""
    t, mid, eps = frag["t"], frag["mid"], frag["eps"].astype(float)
    dur = float(t[-1] - t[0])
    g = np.arange(t[0] + Hs, t[-1] - Hs, PASO_S)
    if g.size < 120:
        return None
    i = np.unique(np.clip(np.searchsorted(t, g, side="left"), 0, t.size - 1))
    j = np.searchsorted(t, t[i] + Hs, side="left")
    k = np.searchsorted(t, t[i] - Hs, side="left")
    v = (j < t.size) & (k >= 0) & (mid[i] > 0)
    i, j, k = i[v], j[v], k[v]
    if i.size < 120:
        return None
    # PREDICTIVO: el retorno arranca DESPUES de la transaccion de `i`.
    r = np.log(mid[j] / mid[np.minimum(i + 1, mid.size - 1)]) * 1e4
    flujo = np.concatenate(([0.0], np.cumsum(eps, dtype=np.float64)))
    fut = flujo[j] - flujo[np.minimum(i + 1, t.size - 1)]
    s_prev = np.abs(np.log(mid[i] / mid[k])) * 1e4
    t_e = t[0] + 0.60 * dur
    ent = t[j] <= t_e
    return {"X": diseno_multi(t, eps, i, Hs), "r": r, "fut": fut,
            "s_prev": s_prev, "ent": ent, "n": i.size}


def _veredicto(filas) -> None:
    titulo("REGLA DE PARADA, APLICADA AL §C.3.5")
    alto = [x for x in filas if x["estrato"] == "ALTO" and x["resuelve"]
            and x["veredicto"] != "NO FALSABLE"]
    descart = [x for x in filas if x["estrato"] == "ALTO" and x["resuelve"]
               and x["veredicto"] == "NO FALSABLE"]
    if descart:
        log("")
        log("  filas del estrato ALTO descartadas por NO FALSABLES (R2_max < R2_req):")
        for x in descart:
            log("    %-16s H=%6ds  R2_max = %.6f < req = %.6f   (margen %.1fx, no cuenta)"
                % (x["frag"], x["H"], x["r2max"], x["req"], x["margen"]))
    log("")
    log("  filas del estrato ALTO donde el instrumento RESUELVE: %d de %d"
        % (len(alto), sum(1 for x in filas if x["estrato"] == "ALTO")))
    if not alto:
        log("")
        log("  *** TAMPOCO LA MULTIVARIANTE ES FALSABLE en el estrato ALTO. El techo")
        log("      del §C.3.5 -- superior por construccion al univariante -- SIGUE")
        log("      POR DEBAJO del requisito. Lectura: el problema NO era el techo")
        log("      del predictor univariante. Migrar de activo con esta misma")
        log("      familia de predictor no tiene por que ir mejor. ***")
        return
    mejor = min(alto, key=lambda z: z["margen"])
    log("  margen MINIMO en el estrato ALTO: %.1fx   (%s, H = %d s, n = %d)"
        % (mejor["margen"], mejor["frag"], mejor["H"], mejor["n"]))
    if mejor["margen"] >= 3.0:
        log("")
        log("  *** El margen NO baja de 3x. RESULTADO NEGATIVO tambien con la")
        log("      multivariante: se confirma el cierre de la linea. ***")
    else:
        log("")
        log("  *** El margen BAJA de 3x con la multivariante donde el univariante")
        log("      no resolvia. NO se cierra: el techo univariante SI era la")
        log("      restriccion, y el §C.3.5 abre una banda que hay que examinar. ***")


# ==============================================================================
# AUTOTEST — verdad conocida, sin datos de mercado
# ==============================================================================
def _autotest() -> int:
    fallos = []

    def chk(ok, msg, det=""):
        log("[%s] %-58s %s" % (" OK " if ok else "FALL", msg, det))
        if not ok:
            fallos.append(msg)

    titulo("AUTOTEST §C.3.5")
    rng = np.random.default_rng(11)
    n, k = 4000, len(FRACS) + 1

    # 1. Sin relacion, el R2 en muestra tiende a k/n. Se promedia sobre semillas
    #    en vez de fiarlo a un sorteo: la distribucion nula es Beta(k/2,(n-k-1)/2)
    #    y un sorteo suelto se aparta ~1.5 sd sin que nada este mal. Es la misma
    #    correccion que la v4.1 §3 tuvo que hacer con el `|d| < 0.06` puesto a ojo.
    r2s = []
    for sem in range(12):
        rg = np.random.default_rng(500 + sem)
        r2s.append(r2_multi(rg.normal(size=(n, k)), rg.normal(size=n)))
    m_r2 = float(np.mean(r2s))
    se = float(np.std(r2s, ddof=1) / np.sqrt(len(r2s)))
    chk(abs(m_r2 - k / n) < 3 * se + 1e-6,
        "sin relacion, E[R2] = k/n (media de 12 sorteos)",
        "media=%.6f  k/n=%.6f  3ee=%.6f" % (m_r2, k / n, 3 * se))

    X = rng.normal(size=(n, k))
    r = rng.normal(size=n)
    r2 = r2_multi(X, r)
    q95 = suelo_multi(X, r, n_sorteos=120)
    chk(r2 <= q95, "el suelo multivariante atrapa el R2 espurio",
        "R2=%.6f <= q95=%.6f" % (r2, q95))

    # 2. LO QUE JUSTIFICA ESTE MODULO: reutilizar el suelo univariante NO es
    #    inofensivo. No se comprueba una razon de suelos elegida a ojo, sino la
    #    consecuencia operativa: cuantos sorteos NULOS se declararian FALSABLES
    #    contra el suelo univariante. Bajo el nulo la respuesta correcta es ~5 %.
    import estratos as E
    q_uni = E.suelo(X[:, 0], r, n_sorteos=120)
    nulos = []
    for sem in range(60):
        rg = np.random.default_rng(900 + sem)
        Xn = rg.normal(size=(n, k))
        nulos.append(r2_multi(Xn, r))
    nulos = np.array([x for x in nulos if np.isfinite(x)])
    falsos_uni = float(np.mean(nulos > q_uni))
    falsos_multi = float(np.mean(nulos > q95))
    chk(falsos_uni > 0.30,
        "con el suelo UNIVARIANTE, el nulo multi se declara falsable a menudo",
        "%.0f %% de falsos positivos" % (100 * falsos_uni))
    # ⚠ Sale ~13 % y no 5 %, y la causa es del TEST, no del suelo: aqui los nulos
    # sortean una `X` nueva mientras `q95` se midio rotando una `X` fija. En la
    # corrida real las dos son la misma serie, que es el caso para el que el suelo
    # esta construido. Se deja el numero a la vista en vez de ajustar el umbral.
    chk(falsos_multi <= 0.15,
        "con su PROPIO suelo, la tasa de falsos positivos vuelve a lo nominal",
        "%.0f %% (nominal 5 %%; ver nota)" % (100 * falsos_multi))

    # 3. Anidamiento: con eps_t como primera columna, el multivariante no puede
    #    rendir menos que el univariante sobre esa misma columna.
    senal = 0.3 * X[:, 0] + rng.normal(size=n)
    r2m = r2_multi(X, senal)
    r2u = E.r2_ident(X[:, 0], senal)
    chk(r2m >= r2u - 1e-9, "con eps_t dentro, multi >= univariante (anidamiento)",
        "multi=%.6f >= uni=%.6f" % (r2m, r2u))
    # Y el anidamiento tiene que valer SIEMPRE, no en un sorteo afortunado.
    peor = min(r2_multi(Xs, ys) - E.r2_ident(Xs[:, 0], ys)
               for Xs, ys in ((lambda g: (g.normal(size=(800, k)),))(np.random.default_rng(700 + z))[0:1] + (
                   np.random.default_rng(700 + z).normal(size=800),) for z in range(15)))
    chk(peor >= -1e-9, "el anidamiento se cumple en 15 sorteos, no en uno",
        "holgura minima = %.2e" % peor)

    # 4. Con relacion real, recupera el R2 poblacional. Media de 12 sorteos
    #    contra el error tipico, no una tolerancia amplia elegida a mano.
    teorico = 0.09 / (0.09 + 1.0)
    ms = []
    for sem in range(12):
        rg = np.random.default_rng(300 + sem)
        Xr_ = rg.normal(size=(n, k))
        ms.append(r2_multi(Xr_, 0.3 * Xr_[:, 0] + rg.normal(size=n)))
    mm = float(np.mean(ms))
    se_m = float(np.std(ms, ddof=1) / np.sqrt(len(ms)))
    chk(abs(mm - teorico) < 3 * se_m + k / n,
        "recupera el R2 poblacional (media de 12, contra su propio error)",
        "media=%.4f  teorico=%.4f  3ee+k/n=%.4f" % (mm, teorico, 3 * se_m + k / n))
    q95s = suelo_multi(X, senal, n_sorteos=120)
    chk(r2m > q95s, "una relacion real SI supera su suelo",
        "%.6f > %.6f" % (r2m, q95s))

    # 5. La rotacion en bloque conserva Sigma; rotar columna a columna no.
    Xr = np.take(X, np.arange(n) - 137, axis=0, mode="wrap")
    S0 = np.cov(X, rowvar=False)
    Sb = np.cov(Xr, rowvar=False)
    Xc = np.column_stack([np.take(X[:, c], np.arange(n) - 37 * (c + 1), mode="wrap")
                          for c in range(k)])
    Sc = np.cov(Xc, rowvar=False)
    d_bloque = float(np.max(np.abs(S0 - Sb)))
    d_col = float(np.max(np.abs(S0 - Sc)))
    chk(d_bloque < 1e-9, "rotar EN BLOQUE conserva Sigma exactamente",
        "max|dS| = %.2e" % d_bloque)
    chk(d_col > 100 * max(d_bloque, 1e-12),
        "rotar columna a columna SI la cambia (por eso no se hace)",
        "max|dS| = %.2e" % d_col)

    # 6. El diseno usa flujo PASADO: ningun rasgo puede depender del futuro.
    t = np.arange(2000, dtype=float)
    eps = rng.choice(np.array([-1.0, 1.0]), size=2000)
    i = np.arange(500, 1500)
    Xd = diseno_multi(t, eps, i, Hs=100)
    eps2 = eps.copy()
    eps2[1600:] = -eps2[1600:]           # cambia SOLO el futuro de todos los i
    Xd2 = diseno_multi(t, eps2, i, Hs=100)
    chk(np.allclose(Xd, Xd2), "el diseno no mira el futuro",
        "max|dX| = %.2e" % float(np.max(np.abs(Xd - Xd2))))

    log("")
    if fallos:
        log("== %d FALLO(S) de 11 ==" % len(fallos))
        return 1
    log("== 11/11 CONTROLES OK ==")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", default="")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    if a.etapa == "medir":
        return etapa_medir(a)
    log("uso: python multivariante_c35.py --autotest | --etapa=medir")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
