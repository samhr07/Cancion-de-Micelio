# -*- coding: utf-8 -*-
"""
retroalimentacion.py -- como `sigma` y `nu` se combinan para dictar `sigma` futura.

    python retroalimentacion.py --autotest   controles con verdad conocida
    python retroalimentacion.py --etapa=var  el VAR: quien empuja a quien
    python retroalimentacion.py --etapa=etapas  la prediccion por etapas del operador

=========================================================================
PARTE A -- EL VAR.  x_t = [log nu_t, log sigma_t]
=========================================================================

    x_t = c + SUM_i A_i x_{t-i} + e_t

Lo que decide son los terminos CRUZADOS de `A_1`:

    A[sigma <- nu]   la actividad de ahora empuja la volatilidad de despues
    A[nu <- sigma]   la volatilidad de ahora empuja la actividad de despues

Si solo el primero es distinto de cero, `nu` es el motor y `sigma` la respuesta:
hay UNA variable de estado. Si los dos lo son, hay retroalimentacion genuina y
hacen falta las dos. La direccion se contrasta **fuera de muestra por dias
enteros**, no con un test F dentro de muestra: este proyecto ya sabe que un
estadistico dentro de muestra sobre series con memoria larga no significa lo que
parece.

=========================================================================
PARTE B -- LA PREDICCION POR ETAPAS
=========================================================================

Propuesta del operador: predecir el PRESENTE con datos del pasado, medir el
porcentaje de error, y comprobar si ese porcentaje cuadra con el de la
prediccion del FUTURO a partir del presente.

[!] PRECISION NECESARIA, PORQUE SI NO EL TEST ES VACUO. "Predecir sigma_t desde
    t-1" y "predecir sigma_{t+1} desde t" son LA MISMA OPERACION desplazada un
    paso. Bajo estacionariedad sus errores coinciden por construccion y
    comprobarlo no aporta nada. Lo que si tiene contenido -- y es lo que aqui se
    mide -- es la version util de la misma idea:

      (1) ?el error medido en la ventana RECIENTE predice el error de la ventana
          SIGUIENTE?  Si si, hay una barra de error VIVA y autocalibrada.
      (2) ?el error normalizado por su desviacion predicha tiene varianza 1?
          Es el NIS de la Sec.7.1, ya en el vocabulario del proyecto: mide si la
          incertidumbre declarada es honesta y no solo pequena.
      (3) ?predice `nu` el TAMANO del error, ademas del valor?  Si el error es
          mayor cuando hay mas actividad, `nu` condiciona a la vez el nivel y la
          confianza.

    Eso es lo que hace falta para el Sec.1: `R2_req ~ sigma^-2`, asi que una
    `sigma` futura con barra de error honesta convierte la curva requerida en una
    senal condicional conocida POR ADELANTADO.
"""

from __future__ import annotations

import argparse

import numpy as np

import estacionalidad as E
import horizonte as H

log, titulo = H.log, H.titulo
P_LAGS = 3          # 3 casillas de 5 min = 15 min de historia


# ===========================================================================
# Construccion de rezagos sobre tramos CONTIGUOS
# ===========================================================================

def matriz_rezagos(d: dict, p: int = P_LAGS, exog=False):
    """`Y` (objetivo en t) y `X` (rezagos 1..p), solo donde hay p+1 contiguas."""
    o = np.argsort(d["t"])
    t = d["t"][o]
    ln = d["log_nu"][o]
    ls = d["log_sigma"][o]
    dia = d["dia"][o]
    hora = d["hora"][o]
    finde = d["finde"][o]
    n = t.size
    ok = np.ones(n, bool)
    ok[:p] = False
    for k in range(1, p + 1):
        dif = np.full(n, np.inf)
        dif[k:] = t[k:] - t[:-k]
        ok &= np.isclose(dif, k * E.BIN_S)
    idx = np.flatnonzero(ok)
    col = [np.ones(idx.size)]
    for k in range(1, p + 1):
        col.append(ln[idx - k])
        col.append(ls[idx - k])
    if exog:
        col += E._armonicos(hora[idx], E.K_ARM)
        col.append(finde[idx])
    X = np.column_stack(col)
    return {"X": X, "y_nu": ln[idx], "y_sigma": ls[idx], "dia": dia[idx],
            "t": t[idx], "hora": hora[idx], "idx": idx, "p": p,
            "ln": ln, "ls": ls}


def pred_oos(X, y, dias, n_pliegues=4, semilla=7):
    """Predicciones fuera de muestra con pliegues de DIAS COMPLETOS."""
    pred = np.full(y.size, np.nan)
    base = np.full(y.size, np.nan)
    u = np.unique(dias)
    rng = np.random.default_rng(semilla)
    for p_ in np.array_split(rng.permutation(u), n_pliegues):
        te = np.isin(dias, p_)
        tr = ~te
        if tr.sum() < X.shape[1] + 5 or te.sum() < 3:
            continue
        c, *_ = np.linalg.lstsq(X[tr], y[tr], rcond=None)
        pred[te] = X[te] @ c
        base[te] = y[tr].mean()
    return pred, base


def r2_de(pred, base, y):
    m = np.isfinite(pred) & np.isfinite(base)
    if m.sum() < 5:
        return float("nan")
    return 1.0 - float(np.sum((y[m] - pred[m]) ** 2)) / float(np.sum((y[m] - base[m]) ** 2))


# ===========================================================================
# Parte A -- el VAR
# ===========================================================================

def etapa_var(args) -> int:
    titulo("PARTE A -- EL VAR: ?quien empuja a quien?")
    d = E.serie_binada()
    R = matriz_rezagos(d, P_LAGS)
    X, dias = R["X"], R["dia"]
    log("")
    log("  %d observaciones con %d rezagos contiguos (%.0f min de historia), %d dias"
        % (X.shape[0], P_LAGS, P_LAGS * E.BIN_S / 60, np.unique(dias).size))

    # --- coeficientes del rezago 1, ajustados sobre TODO (descriptivo)
    log("")
    log("  --- A_1, la matriz de transicion del primer rezago ---")
    A = np.zeros((2, 2))
    for j, (nom, y) in enumerate((("log nu   ", R["y_nu"]), ("log sigma", R["y_sigma"]))):
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        A[j, 0], A[j, 1] = c[1], c[2]          # col 1 = ln(t-1), col 2 = ls(t-1)
        log("    %s(t)  <-  %+.4f * log nu(t-1)  %+.4f * log sigma(t-1)"
            % (nom, c[1], c[2]))
    ev = np.linalg.eigvals(A)
    log("    autovalores de A_1: %s" % np.array2string(ev, precision=4))
    rho = float(np.max(np.abs(ev)))
    if 0 < rho < 1:
        log("    radio espectral %.4f  ->  vida media %.2f casillas = %.1f min"
            % (rho, np.log(0.5) / np.log(rho), np.log(0.5) / np.log(rho) * E.BIN_S / 60))

    # --- direccion, FUERA DE MUESTRA
    log("")
    log("  --- direccion, fuera de muestra por dias enteros ---")
    p = P_LAGS
    col_nu = [0] + [1 + 2 * k for k in range(p)]
    col_si = [0] + [2 + 2 * k for k in range(p)]
    todos = list(range(X.shape[1]))
    log("    %-42s %10s %10s" % ("modelo para el objetivo", "objetivo", "R2 fuera"))
    for nom_obj, y in (("log nu(t)", R["y_nu"]), ("log sigma(t)", R["y_sigma"])):
        for nom, cols in (("solo rezagos de nu", col_nu),
                          ("solo rezagos de sigma", col_si),
                          ("AMBOS (VAR completo)", todos)):
            pr, ba = pred_oos(X[:, cols], y, dias)
            log("    %-42s %10s %+10.4f" % (nom, nom_obj, r2_de(pr, ba, y)))
        log("")
    log("    LECTURA: si anadir los rezagos de nu mejora la prediccion de sigma pero")
    log("    no al reves, nu es el motor y sigma la respuesta -> UNA variable de estado.")

    # --- respuesta al impulso
    log("")
    log("  --- respuesta al impulso (choque unitario, propagado con A_1) ---")
    log("    pasos                    1      2      3      4      6      8     12")
    for nom, e0 in (("choque en nu    -> sigma", np.array([1.0, 0.0])),
                    ("choque en sigma -> nu   ", np.array([0.0, 1.0]))):
        fila, x = [], e0.copy()
        for k in range(1, 13):
            x = A @ x
            if k in (1, 2, 3, 4, 6, 8, 12):
                fila.append(x[1] if "-> sigma" in nom else x[0])
        log("    %-22s %s" % (nom, " ".join("%+6.3f" % v for v in fila)))
    return 0


# ===========================================================================
# Parte B -- la prediccion por etapas
# ===========================================================================

def etapa_etapas(args) -> int:
    titulo("PARTE B -- PREDICCION POR ETAPAS Y BARRA DE ERROR VIVA")
    d = E.serie_binada()
    R = matriz_rezagos(d, P_LAGS)
    X, y, dias = R["X"], R["y_sigma"], R["dia"]

    pred, base = pred_oos(X, y, dias)
    m = np.isfinite(pred)
    err = np.full(y.size, np.nan)
    err[m] = y[m] - pred[m]
    log("")
    log("  --- las dos etapas del operador, medidas ---")
    log("  ETAPA 1  predecir sigma(t) con datos hasta t-1  (nowcast)")
    log("  ETAPA 2  predecir sigma(t+1) con datos hasta t  (forecast, incluye nu(t))")
    rmse = float(np.sqrt(np.nanmean(err ** 2)))
    # error relativo en la escala original: sigma es lognormal, e^err - 1
    rel = float(np.nanmedian(np.abs(np.exp(err) - 1.0)))
    log("")
    log("    RMSE en log sigma           : %.4f" % rmse)
    log("    error relativo mediano      : %.2f %%" % (100 * rel))
    log("    R2 fuera de muestra         : %+.4f" % r2_de(pred, base, y))
    # etapa 2 = la misma operacion desplazada. Se mide y se compara.
    o = np.argsort(R["t"])
    t_, e_, di_ = R["t"][o], err[o], dias[o]
    cons = np.isclose(np.diff(t_), E.BIN_S)
    i0 = np.flatnonzero(cons)
    e1, e2 = e_[i0], e_[i0 + 1]
    v = np.isfinite(e1) & np.isfinite(e2)
    r1 = float(np.sqrt(np.mean(e1[v] ** 2)))
    r2_ = float(np.sqrt(np.mean(e2[v] ** 2)))
    log("")
    log("    RMSE etapa 1 (presente) = %.4f     RMSE etapa 2 (futuro) = %.4f" % (r1, r2_))
    log("    razon futuro/presente   = %.4f" % (r2_ / max(r1, 1e-12)))
    log("    correlacion de errores consecutivos = %+.4f" % float(np.corrcoef(e1[v], e2[v])[0, 1]))
    log("")
    log("    [!] Que la razon salga ~1 NO es un hallazgo: las dos etapas son la misma")
    log("        operacion desplazada un paso. Lo que decide es lo de abajo.")

    # --- (1) el error reciente predice el error siguiente
    log("")
    log("  --- (1) ?el error de la ventana RECIENTE predice el de la SIGUIENTE? ---")
    for W in (3, 6, 12, 24):
        rms_p, rms_f, dd = [], [], []
        for q in np.unique(di_):
            s = np.flatnonzero(di_ == q)
            ee = e_[s]
            for a in range(0, ee.size - 2 * W, W):
                p1, p2 = ee[a:a + W], ee[a + W:a + 2 * W]
                if np.isfinite(p1).sum() < W - 1 or np.isfinite(p2).sum() < W - 1:
                    continue
                rms_p.append(np.sqrt(np.nanmean(p1 ** 2)))
                rms_f.append(np.sqrt(np.nanmean(p2 ** 2)))
                dd.append(q)
        rms_p, rms_f, dd = np.array(rms_p), np.array(rms_f), np.array(dd)
        if rms_p.size < 20:
            continue
        cc = float(np.corrcoef(np.log(rms_p), np.log(rms_f))[0, 1])
        Xw = np.column_stack([np.ones(rms_p.size), np.log(rms_p)])
        pr, ba = pred_oos(Xw, np.log(rms_f), dd)
        log("    W = %2d casillas (%3.0f min)  n = %4d   corr = %+.4f   R2 fuera = %+.4f"
            % (W, W * E.BIN_S / 60, rms_p.size, cc, r2_de(pr, ba, np.log(rms_f))))

    # --- (2) NIS: ?es honesta la barra de error?
    log("")
    log("  --- (2) NIS: ?la incertidumbre declarada es HONESTA? ---")
    log("      z = error / sd predicha. Si la sd es honesta, var(z) ~ 1.")
    # sd predicha: se modela log|err| con los mismos rezagos, fuera de muestra
    # [!] CORRECCION DE JENSEN, Y SIN ELLA EL NIS SALE MAL POR UN FACTOR 2.27.
    # Se modela log|e|, no |e|. Para e ~ N(0, s):  log|e| = log s + log|z|  con
    # E[log|z|] = -(gamma + log 2)/2 = -0.63518. Asi que la desviacion es
    #     s = exp( E[log|e|] + 0.63518 )
    # y NO exp(E[log|e|])*sqrt(pi/2), que es la conversion valida para E|e| y no
    # para exp(E[log|e|]). El error era de exp(0.6352-0.2258) = 1.506x en `s`,
    # o sea 2.27x en var(z): lo bastante para declarar deshonesta una barra de
    # error que no lo era tanto.
    CORR_LOGNORMAL = 0.6351814227307392
    le = np.log(np.abs(err) + 1e-6)
    prl, bal = pred_oos(X, le, dias)
    sd_pred = np.exp(prl + CORR_LOGNORMAL)
    z = err / sd_pred
    zz = z[np.isfinite(z)]
    zc = (err / rmse)[np.isfinite(err)]
    log("      var(z) con sd CONSTANTE  = %.4f   (1.0 por construccion, no informa)"
        % float(np.var(zc)))
    log("      var(z) con sd PREDICHA   = %.4f   (1.0 = honesta)" % float(np.var(zz)))
    log("      curtosis: sd constante %.2f  ->  sd predicha %.2f"
        % (float(np.mean((zc - zc.mean()) ** 4) / np.var(zc) ** 2),
           float(np.mean((zz - zz.mean()) ** 4) / np.var(zz) ** 2)))
    log("      fraccion |z| > 3 : constante %.3f %%   predicha %.3f %%   (normal: 0.270 %%)"
        % (100 * float(np.mean(np.abs(zc) > 3)), 100 * float(np.mean(np.abs(zz) > 3))))

    # --- (3) nu predice el TAMANO del error
    log("")
    log("  --- (3) ?predice nu el TAMANO del error, ademas del valor? ---")
    ln_t1 = X[:, 1]                                   # log nu(t-1)
    for nom, Xe in (("constante        ", np.ones((le.size, 1))),
                    ("log nu(t-1)      ", np.column_stack([np.ones(le.size), ln_t1])),
                    ("rezagos completos", X)):
        pr, ba = pred_oos(Xe, le, dias)
        log("    %-18s -> R2 fuera de muestra sobre log|error| = %+.4f"
            % (nom, r2_de(pr, ba, le)))
    log("")
    log("    Si nu predice el tamano del error, condiciona a la vez el nivel de sigma")
    log("    y la confianza en el -- que es lo que el Sec.1 necesita, porque")
    log("    R2_req ~ sigma^-2 se vuelve una senal conocida POR ADELANTADO.")
    return 0


# ===========================================================================
# Controles
# ===========================================================================

def _autotest() -> int:
    titulo("CONTROLES DE retroalimentacion.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-56s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    rng = np.random.default_rng(4)
    n_dias, por_dia = 24, 200
    A_v = np.array([[0.60, 0.00],      # nu(t) <- 0.6 nu(t-1),  0 sigma(t-1)
                    [0.35, 0.40]])     # sigma(t) <- 0.35 nu(t-1) + 0.4 sigma(t-1)
    x = np.zeros((n_dias * por_dia, 2))
    for i in range(1, x.shape[0]):
        x[i] = A_v @ x[i - 1] + 0.3 * rng.normal(size=2)
    t = np.arange(x.shape[0]) * E.BIN_S
    dia = np.repeat(np.arange(n_dias), por_dia)
    # los dias no son contiguos entre si: se corta el enlace
    t = t + dia * 10 * E.BIN_S
    d = {"t": t, "log_nu": x[:, 0], "log_sigma": x[:, 1], "dia": dia,
         "hora": np.mod(t / 3600.0, 24.0), "finde": np.zeros(t.size)}
    R = matriz_rezagos(d, 1)
    X = R["X"]
    cn, *_ = np.linalg.lstsq(X, R["y_nu"], rcond=None)
    cs, *_ = np.linalg.lstsq(X, R["y_sigma"], rcond=None)
    log("     A recuperada: [[%+.3f %+.3f] [%+.3f %+.3f]]  verdad [[%.2f %.2f] [%.2f %.2f]]"
        % (cn[1], cn[2], cs[1], cs[2], A_v[0, 0], A_v[0, 1], A_v[1, 0], A_v[1, 1]))
    chk(abs(cn[1] - 0.60) < 0.05, "recupera A[nu<-nu] = 0.60", "%+.3f" % cn[1])
    chk(abs(cn[2] - 0.00) < 0.05, "recupera A[nu<-sigma] = 0 (sin retroalimentacion)",
        "%+.3f" % cn[2])
    chk(abs(cs[0 + 1] - 0.35) < 0.05, "recupera A[sigma<-nu] = 0.35", "%+.3f" % cs[1])
    chk(abs(cs[2] - 0.40) < 0.05, "recupera A[sigma<-sigma] = 0.40", "%+.3f" % cs[2])

    # la direccion se detecta fuera de muestra
    pr, ba = pred_oos(X[:, [0, 1]], R["y_sigma"], R["dia"])
    r_nu = r2_de(pr, ba, R["y_sigma"])
    pr, ba = pred_oos(X[:, [0, 2]], R["y_sigma"], R["dia"])
    r_si = r2_de(pr, ba, R["y_sigma"])
    chk(r_nu > 0.05, "nu solo ya predice sigma fuera de muestra", "%+.4f" % r_nu)
    pr, ba = pred_oos(X[:, [0, 2]], R["y_nu"], R["dia"])
    r_cruz = r2_de(pr, ba, R["y_nu"])
    chk(r_cruz < 0.05, "sigma NO predice nu (la direccion se detecta)", "%+.4f" % r_cruz)

    # --- barra de error viva: con varianza que cambia por bloques, el error
    #     reciente debe predecir el siguiente; con varianza constante, no.
    for etq, esperado, sd in (("varianza por bloques", True,
                               np.repeat(rng.uniform(0.2, 2.0, 200), 24)),
                              ("varianza constante  ", False,
                               np.ones(200 * 24))):
        e = sd * rng.normal(size=sd.size)
        W = 12
        a = np.array([np.sqrt(np.mean(e[i:i + W] ** 2))
                      for i in range(0, e.size - 2 * W, W)])
        b = np.array([np.sqrt(np.mean(e[i + W:i + 2 * W] ** 2))
                      for i in range(0, e.size - 2 * W, W)])
        cc = float(np.corrcoef(np.log(a), np.log(b))[0, 1])
        chk((cc > 0.3) == esperado,
            "error reciente predice el siguiente: %s" % etq, "corr %+.4f" % cc)

    log("")
    log("RESULTADO: %d fallo(s)" % fallos)
    return 1 if fallos else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", default="")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    if a.etapa == "var":
        return etapa_var(a)
    if a.etapa == "etapas":
        return etapa_etapas(a)
    log("usa --autotest, --etapa=var o --etapa=etapas")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
