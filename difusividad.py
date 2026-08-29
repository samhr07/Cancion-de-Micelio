# -*- coding: utf-8 -*-
"""
difusividad.py -- ¿sigue siendo difusivo el precio MAS ALLA de 12 minutos?

    python difusividad.py
    python difusividad.py --autotest

POR QUE IMPORTA, Y MUCHO
------------------------
La enmienda 9 del `PREREGISTRO_3_2.md` convirtio `H*` de "horizonte de operacion"
en "suelo por debajo del cual operar es imposible", y apoyo esa lectura en una
tabla de R2 requerido que usa `sigma(H) = sigma_1 * sqrt(H)` hasta 1 h, 4 h y
dias. Esa extrapolacion **descansa en que el precio sea difusivo a esas escalas**,
y la difusividad esta medida solo hasta 8 192 ticks (~12 min).

Si a 1 h el precio revierte, `sigma(H)` crece menos que `sqrt(H)`, el R2 exigido
sube y la banda que la enmienda declaro "linea principal" se estrecha o
desaparece. Si tiende (momentum), baja aun mas. **La enmienda es tan buena como
esta comprobacion**, y por eso se hace antes de construir nada encima.

QUE SE MIDE
-----------
La firma de volatilidad **en tiempo de ticks**, `sigma(n)/sqrt(n)` contra `n`.
Plana = difusivo. La sesion 2026-08-08 (e) establecio que en reloj de PARED el
estimador tiene sesgo propio y los muestreos solapado y no solapado discrepan en
signo; en reloj de ticks coinciden a tres decimales. Aqui se usan los dos y se
comprueba que siguen coincidiendo a `n` grande.

⚠ LO QUE LIMITA EL ALCANCE NO ES EL DATO TOTAL, SON LAS VENTANAS INDEPENDIENTES.
A `n` grande caben pocas ventanas no solapadas y la estimacion se vuelve ruido.
Se reporta `n_indep` en cada fila y se marca la frontera a partir de la cual el
numero cae por debajo de un minimo declarado.

⚠ RETRACTACION (2026-08-09, posterior al commit b33fb85). Este modulo reporto
"tres regimenes" -- reversion en `captura_larga`, difusion en `v31b`, momentum en
`v32` -- y **esa lectura queda retirada**. Dos defectos, y el segundo es fatal
para cualquier conclusion:

1. **El rango de ajuste estaba fijo en TICKS, asi que la banda en SEGUNDOS
   difiere por factor 7 entre capturas.** A `nu = 39` el ajuste de
   `captura_larga` arrancaba en 6.6 s, dentro de la microestructura donde el
   rebote bid-ask produce reversion (`R(1) = -0.0052`, ya medido en la v3.1).
   Con banda comun [44, 202] s la reversion **desaparece**: -0.2311 -> +0.0158.
   Quinta aparicion del mismo patron: mezclar relojes en un estimador.

2. **Y sobre banda comun la pendiente no es estimable.** Cambiando solo la
   densidad de la rejilla, con los mismos datos y la misma banda:

       captura_larga:  +0.016 (4 pts) -> +1.074 (11) -> +1.519 (23)
       captura_v31b:   +0.129        -> +0.072       -> +0.107
       captura_v32:    +0.207        -> +0.120       -> +0.169

   Una pendiente que se mueve de +0.02 a +1.52 con el numero de puntos no es una
   medicion. La banda comun cubre 4.6x en escala (log-rango 1.52) y la curva
   tiene curvatura: no hay brazo de palanca para una ley de potencias.

**Lo que sobrevive:** el alcance verificado sigue siendo ~12 min, el estimador no
tiene sesgo propio (controles barajados en ~0), y la difusividad a escalas largas
**sigue sin verificarse** -- que era la pregunta. Lo que NO sobrevive es la
afirmacion de que hay tres regimenes distintos.

**Lo que haria falta**, y sale de la captura estacional: bandas de al menos una
decada en escala, y ~100 estimaciones independientes a 10 min para tener
distribucion muestral empirica del estimador en vez de discutir sobre tres puntos.

⚠ NO SE TOCA `captura_v33`. Es el dato de la decision de la v3.2 y mirarlo aqui
seria contaminarlo. Se usan las capturas ya gastadas como banco de pruebas.
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np

import captura_larga as cl
import propagador as pr


# ⚠ El limite REAL son 30 ventanas no solapadas, no 8: `pr.firma_en_ticks`
# descarta toda fila con menos de 30 muestras, y en muestreo no solapado las
# muestras SON las ventanas. Este 8 solo acota la rejilla; el filtro que manda
# es aquel, y por eso el alcance se queda en ~12 min pese a tener 447 698 ticks
# continuos: lo que escasea no es el dato, son las ventanas independientes.
MIN_INDEP = 8          # acota la rejilla; el filtro efectivo son 30 ventanas
CAPTURAS = ("captura_larga", "captura_v31b", "captura_v32")


def tramo_mayor(nombre: str) -> dict:
    """Devuelve el tramo continuo mas largo de una captura."""
    d = cl.cargar_larga("telemetria/" + nombre)
    t = d["tr_t"]
    tramos = cl.tramos_continuos(t)
    L = [b - a for a, b in tramos]
    k = int(np.argmax(L))
    a, b = tramos[k]
    dur = float(t[b - 1] - t[a])
    return {"nombre": nombre, "precio": d["tr_precio"][a:b], "t": t[a:b],
            "n": b - a, "dur_s": dur, "nu": (b - a) / dur if dur > 0 else np.nan}


def malla(n_max: int) -> list:
    """Rejilla geometrica de `n` en ticks, hasta `n_max`."""
    ns, n = [], 2
    while n <= n_max:
        ns.append(n)
        n *= 2
    return ns


def firma(tr: dict, min_indep: int = MIN_INDEP) -> list:
    """Firma en ticks con las dos formas de muestreo y el numero de ventanas."""
    ns = malla(tr["n"] // min_indep)
    sol = {f["n"]: f for f in pr.firma_en_ticks(tr["precio"], ns, solapada=True)}
    nos = {f["n"]: f for f in pr.firma_en_ticks(tr["precio"], ns, solapada=False)}
    filas = []
    for n in ns:
        if n not in sol or n not in nos:
            continue
        ind = max(1, (tr["n"] - n) // n)
        filas.append({"n": n, "n_indep": ind,
                      "seg": n / tr["nu"],
                      "sol": sol[n]["sigma_por_raiz_n"],
                      "nosol": nos[n]["sigma_por_raiz_n"],
                      "sigma": sol[n]["sigma"],
                      "fiable": ind >= min_indep})
    return filas


def pendiente(filas: list, n_lo: int, n_hi: int, clave: str = "nosol") -> dict:
    """Pendiente de log[sigma(n)/sqrt(n)] contra log n. Cero = difusivo."""
    n = np.array([f["n"] for f in filas], float)
    y = np.array([f[clave] for f in filas], float)
    m = (n >= n_lo) & (n <= n_hi) & (y > 0)
    if m.sum() < 3:
        return {"pendiente": float("nan"), "puntos": int(m.sum())}
    X = np.column_stack([np.log(n[m]), np.ones(int(m.sum()))])
    b, res, *_ = np.linalg.lstsq(X, np.log(y[m]), rcond=None)
    pred = X @ b
    ss = float(np.sum((np.log(y[m]) - pred) ** 2))
    gl = max(1, int(m.sum()) - 2)
    cov = (ss / gl) * np.linalg.inv(X.T @ X)
    return {"pendiente": float(b[0]), "ee": float(np.sqrt(cov[0, 0])),
            "H_p": float(b[0] + 0.5), "puntos": int(m.sum())}


def barajar(precios: np.ndarray, semilla: int = 0) -> np.ndarray:
    """Control: incrementos barajados. Difusivo por construccion."""
    rng = np.random.default_rng(semilla)
    dp = np.diff(np.asarray(precios, float))
    return np.concatenate(([precios[0]], precios[0] + np.cumsum(rng.permutation(dp))))


def analizar(nombre: str, n_lo_frac: float = 0.125) -> dict:
    tr = tramo_mayor(nombre)
    f = firma(tr)
    if not f:
        return {}
    fiables = [x for x in f if x["fiable"]]
    n_max = max(x["n"] for x in fiables) if fiables else 0
    n_lo = max(256, int(n_max * n_lo_frac))

    tr_b = dict(tr)
    tr_b["precio"] = barajar(tr["precio"], 7)
    fb = firma(tr_b)

    return {"tramo": tr, "filas": f, "filas_barajado": fb,
            "n_max_fiable": n_max, "n_lo": n_lo,
            "pend_nosol": pendiente(f, n_lo, n_max, "nosol"),
            "pend_sol": pendiente(f, n_lo, n_max, "sol"),
            "pend_barajado": pendiente(fb, n_lo, n_max, "nosol"),
            "seg_max": n_max / tr["nu"] if tr["nu"] > 0 else float("nan")}


def informe(res: dict) -> None:
    tr = res["tramo"]
    print("-" * 78)
    print("%s | %d ticks continuos en %.2f h | nu = %.2f tx/s"
          % (tr["nombre"], tr["n"], tr["dur_s"] / 3600.0, tr["nu"]))
    print("-" * 78)
    print("%9s %9s %11s %11s %11s %8s"
          % ("n [ticks]", "= seg", "sigma", "sol/sqrt(n)", "nosol", "ventanas"))
    for f in res["filas"]:
        marca = "" if f["fiable"] else "  <- pocas"
        print("%9d %9.0f %11.4f %11.4f %11.4f %8d%s"
              % (f["n"], f["seg"], f["sigma"], f["sol"], f["nosol"],
                 f["n_indep"], marca))
    print("")
    print("rango de ajuste: n en [%d, %d] ticks = [%.0f, %.0f] s"
          % (res["n_lo"], res["n_max_fiable"],
             res["n_lo"] / tr["nu"], res["seg_max"]))
    for etq, k in (("no solapada", "pend_nosol"), ("solapada", "pend_sol"),
                   ("BARAJADO (control)", "pend_barajado")):
        p = res[k]
        print("   pendiente %-20s %+.4f +- %.4f   (H_p = %.3f)"
              % (etq, p["pendiente"], p.get("ee", float("nan")), p["H_p"]))
    d = abs(res["pend_nosol"]["pendiente"] - res["pend_sol"]["pendiente"])
    print("   |solapada - no solapada| = %.4f  (en reloj de PARED discrepaban"
          " en SIGNO)" % d)
    print("   ALCANCE VERIFICADO: hasta %.0f s = %.1f min"
          % (res["seg_max"], res["seg_max"] / 60.0))


def _autotest() -> int:
    fallos = 0

    def ok(nombre, cond, detalle=""):
        nonlocal fallos
        if not cond:
            fallos += 1
        print("  [%s] %s %s" % ("OK  " if cond else "FALLA", nombre, detalle))

    rng = np.random.default_rng(3)
    n = 200000

    print("== 1. CONTROL POSITIVO: paseo aleatorio -> pendiente 0 ==")
    p = np.cumsum(rng.normal(0, 0.2, n))
    tr = {"nombre": "sintetico", "precio": p, "t": np.arange(n) / 10.0,
          "n": n, "dur_s": n / 10.0, "nu": 10.0}
    f = firma(tr)
    nm = max(x["n"] for x in f if x["fiable"])
    pe = pendiente(f, 256, nm, "nosol")
    print("     pendiente = %+.4f +- %.4f (verdad 0)" % (pe["pendiente"], pe["ee"]))
    ok("difusivo -> pendiente ~ 0", abs(pe["pendiente"]) < 0.03)

    print("== 2. CONTROL POSITIVO: proceso con REVERSION -> pendiente < 0 ==")
    # AR(1) sobre el nivel: revierte, asi que sigma(n) satura y la firma cae.
    a = 0.9995
    y = np.empty(n); y[0] = 0.0
    e = rng.normal(0, 0.2, n)
    for i in range(1, n):
        y[i] = a * y[i - 1] + e[i]
    trr = dict(tr); trr["precio"] = y
    fr = firma(trr)
    nmr = max(x["n"] for x in fr if x["fiable"])
    per = pendiente(fr, 256, nmr, "nosol")
    print("     pendiente = %+.4f +- %.4f (deberia ser < 0)"
          % (per["pendiente"], per["ee"]))
    ok("reversion -> pendiente negativa", per["pendiente"] < -0.05)

    print("== 3. CONTROL POSITIVO: momentum de RANGO CORTO ==")
    # ⚠ Este control estaba mal disenado en su primera version: se ajustaba desde
    # n = 256 un proceso cuyo momentum vive en los primeros rezagos, y salia
    # plano. No era un fallo del estimador -- era buscar el efecto donde no
    # esta. Ahora se comprueban LAS DOS cosas, que es ademas la situacion real:
    # rebote/momentum a rezago corto y difusion mas alla.
    e2 = rng.normal(0, 0.2, n)
    r = np.empty(n); r[0] = e2[0]
    for i in range(1, n):
        r[i] = 0.25 * r[i - 1] + e2[i]
    trm = dict(tr); trm["precio"] = np.cumsum(r)
    fm = firma(trm)
    nmm = max(x["n"] for x in fm if x["fiable"])
    corto = pendiente(fm, 2, 64, "nosol")
    largo = pendiente(fm, 256, nmm, "nosol")
    print("     rango CORTO  [2, 64]   : %+.4f +- %.4f  (debe ser > 0)"
          % (corto["pendiente"], corto["ee"]))
    print("     rango LARGO [256, %d] : %+.4f +- %.4f  (debe ser ~ 0)"
          % (nmm, largo["pendiente"], largo["ee"]))
    ok("detecta el momentum donde vive", corto["pendiente"] > 0.02)
    ok("y no lo extiende a donde no esta", abs(largo["pendiente"]) < 0.06)

    print("== 4. Los dos muestreos coinciden en reloj de TICKS ==")
    ps = pendiente(f, 256, nm, "sol")
    print("     no solapada %+.4f | solapada %+.4f | dif %.4f"
          % (pe["pendiente"], ps["pendiente"],
             abs(pe["pendiente"] - ps["pendiente"])))
    ok("difieren menos de 0.02", abs(pe["pendiente"] - ps["pendiente"]) < 0.02)

    print("")
    print("RESULTADO: %d fallo(s)" % fallos)
    return fallos


def main(argv) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--capturas", default=",".join(CAPTURAS))
    a = ap.parse_args(argv[1:])
    if a.autotest:
        return 1 if _autotest() else 0

    print("=" * 78)
    print("DIFUSIVIDAD MAS ALLA DE 12 MINUTOS")
    print("=" * 78)
    print("La enmienda 9 del PREREGISTRO_3_2 extrapola sigma(H) = sigma_1*sqrt(H)")
    print("hasta 1 h, 4 h y dias. Esto comprueba hasta donde aguanta.")
    print("AVISO: captura_v33 NO se toca -- es el dato de la decision v3.2.")
    print("")
    resultados = {}
    for nom in a.capturas.split(","):
        try:
            r = analizar(nom.strip())
        except Exception as err:
            print("%s: no analizable (%s)" % (nom, err))
            continue
        if r:
            informe(r)
            print("")
            resultados[nom.strip()] = r

    print("=" * 78)
    print("RESUMEN")
    print("=" * 78)
    print("%-16s %10s %12s %14s %14s"
          % ("captura", "nu", "alcance", "pendiente", "barajado"))
    for nom, r in resultados.items():
        print("%-16s %10.2f %9.0f s %+14.4f %+14.4f"
              % (nom, r["tramo"]["nu"], r["seg_max"],
                 r["pend_nosol"]["pendiente"], r["pend_barajado"]["pendiente"]))
    if resultados:
        alcance = max(r["seg_max"] for r in resultados.values())
        print("")
        print("ALCANCE MAXIMO VERIFICADO: %.0f s = %.1f min" % (alcance, alcance / 60))
        print("Mas alla de eso la tabla de R2 de la enmienda 9 SIGUE SIENDO")
        print("EXTRAPOLACION, y asi debe citarse.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
