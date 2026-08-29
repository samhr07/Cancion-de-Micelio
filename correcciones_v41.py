# -*- coding: utf-8 -*-
"""v4.1 Sec.3 -- Correcciones al reporte de la v3.3.

    python correcciones_v41.py

Ninguna cambia un veredicto. Todas cambian un numero que otros documentos ya
estan citando. Las cuatro que pide el Sec.3:

    3.1  tercer estimador de `gamma` (GPH y Whittle local), contra la
         prediccion de sesgo congelada en PREDICCION_SESGO_GAMMA_4_1.md
    3.2  tabla de `N_eff` con la definicion COHERENTE `N/inflacion^2`
    3.3  `gamma` y `H` en la MISMA ventana de escala, y `beta` implicita
         como CURVA contra la ventana, no como punto
    3.4  reclasificacion de `D = 11.53` (es texto, va al reporte)

[!] Se corre sobre DOS series:
      A) `captura_v33` entrenamiento -- la que la v3.3 midio. Es la unica
         comparacion valida para «corregir» sus cifras.
      B) el tramo continuo de 23.33 h de `captura_estacional` -- replica
         INDEPENDIENTE, otro tramo y otro regimen de `nu`.
    La v3.3 dejo escrito que la heterogeneidad entre tramos era el problema de
    fondo; con una sola serie no se puede saber si una correccion transporta.

`Micelio.py` no se toca.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

import tick_grande as T
from tick_grande import log, titulo

CACHE = "telemetria/muestra_v32.npz"
DIR_EST = "telemetria/estacional"
CORTE = 300.0            # s, el mismo criterio de tramo continuo de la v3.2


# ---------------------------------------------------------------------------
# Carga
# ---------------------------------------------------------------------------

def serie_v33():
    """Entrenamiento de `captura_v33`: exactamente la muestra que la v3.3 uso."""
    d = np.load(CACHE)
    eps, precio, t = d["eps"], d["precio"], d["t"]
    i1 = int(eps.size * 0.60)
    ent = np.arange(0, i1)
    nu = eps.size / float(t[-1] - t[0])
    return {"nombre": "captura_v33 (entrenamiento)", "eps": eps[ent].astype(float),
            "precio": precio[ent], "nu": nu}


def serie_estacional():
    """Tramo continuo mas largo de `captura_estacional`.

    [!] Se usa UN SOLO tramo continuo. Concatenar tramos separados por horas
    fabricaria un salto de precio en la costura que la firma de volatilidad
    leeria como un incremento gigante, y la ACF de signos como una
    discontinuidad. Es el mismo error que costo la captura larga de la v2.2.
    """
    import pyarrow.parquet as pq
    trozos = []
    for sub in sorted(d for d in os.listdir(DIR_EST) if d.startswith("trades_")):
        for f in sorted(glob.glob(os.path.join(DIR_EST, sub, "*.parquet"))):
            tb = pq.read_table(f, columns=["t", "precio", "maker"])
            trozos.append((tb["t"].to_numpy().astype(float),
                           tb["precio"].to_numpy().astype(float),
                           tb["maker"].to_numpy().astype(bool)))
    t = np.concatenate([x[0] for x in trozos])
    p = np.concatenate([x[1] for x in trozos])
    mk = np.concatenate([x[2] for x in trozos])
    o = np.argsort(t, kind="stable")
    t, p, mk = t[o], p[o], mk[o]
    val = (p > 0)                     # el feed emite p = 0 (v3.0, hallazgo de datos)
    t, p, mk = t[val], p[val], mk[val]

    corte = np.where(np.diff(t) > CORTE)[0]
    bordes = np.concatenate(([0], corte + 1, [t.size]))
    dur = np.array([t[bordes[i + 1] - 1] - t[bordes[i]]
                    for i in range(len(bordes) - 1)])
    k = int(np.argmax(dur))
    sl = slice(bordes[k], bordes[k + 1])
    t, p, mk = t[sl], p[sl], mk[sl]
    # Convencion de la v3.2 Sec.4.2, verificada alli: maker=True -> eps = -1.
    eps = np.where(mk, -1.0, +1.0)
    return {"nombre": "captura_estacional (tramo continuo mas largo)",
            "eps": eps, "precio": p, "nu": p.size / float(t[-1] - t[0]),
            "horas": (t[-1] - t[0]) / 3600.0}


# ---------------------------------------------------------------------------
# Sec.3.1
# ---------------------------------------------------------------------------

def sec_3_1(s):
    titulo("Sec.3.1 -- gamma medida por TRES vias, contra la prediccion congelada")
    log("Prediccion en PREDICCION_SESGO_GAMMA_4_1.md, commit f1fd84f, escrita")
    log("antes de calcular ningun GPH ni ningun Whittle local.")
    log("")
    eps = s["eps"]

    a = T.gamma_de_signos(eps)                       # A: log-log de C(l)
    esc = T.ic_memoria_larga(eps)                    # B: escalado de Var
    g = T.gph(eps)                                   # C1: GPH
    w = T.whittle_local(eps)                         # C2: Whittle local

    log("  metodo                                      gamma      d       dominio")
    log("  ---------------------------------------------------------------------")
    log("  A  regresion log-log de C(l), [10, 2000]   %+.4f  %+.4f   tiempo"
        % (a["gamma"], T.d_desde_gamma(a["gamma"])))
    log("  B  escalado de Var(media de bloque)        %+.4f  %+.4f   tiempo"
        % (esc["gamma_eff"], T.d_desde_gamma(esc["gamma_eff"])))
    log("  C1 GPH            (m = N^0.5 = %5d)      %+.4f  %+.4f   frecuencia"
        % (g["m"], g["gamma"], g["d"]))
    log("  C2 Whittle local  (m = N^0.5 = %5d)      %+.4f  %+.4f   frecuencia"
        % (w["m"], w["gamma"], w["d"]))
    log("     errores estandar asintoticos: GPH +-%.4f   LW +-%.4f  (en gamma)"
        % (g["se_gamma"], w["se_gamma"]))
    log("     C(1) (lo que la v3.2 llamo gamma_hat) = %+.4f  -- NO es el exponente"
        % a["C1"])

    log("")
    log("  --- P1: direccion.  gamma_espectral < %.4f ? ---" % a["gamma"])
    p1 = (g["gamma"] < a["gamma"]) and (w["gamma"] < a["gamma"])
    log("      GPH %.4f %s %.4f    LW %.4f %s %.4f   ->  P1 %s"
        % (g["gamma"], "<" if g["gamma"] < a["gamma"] else ">=", a["gamma"],
           w["gamma"], "<" if w["gamma"] < a["gamma"] else ">=", a["gamma"],
           "SE CUMPLE" if p1 else "REFUTADA"))

    log("")
    log("  --- P2: magnitud.  gamma_espectral en [0.28, 0.46] ? ---")
    p2 = all(0.28 <= v <= 0.46 for v in (g["gamma"], w["gamma"]))
    log("      GPH %.4f   LW %.4f   ->  P2 %s"
        % (g["gamma"], w["gamma"], "SE CUMPLE" if p2 else "REFUTADA"))

    log("")
    log("  --- P3: dependencia de la banda.  gamma CRECE con m ? ---")
    bar = T.barrido_banda(eps)
    log("      alpha      m     gamma GPH  (+-se)    gamma LW  (+-se)")
    for f in bar:
        log("      %.1f   %7d   %+.4f (%.4f)   %+.4f (%.4f)"
            % (f["alpha"], f["m"], f["gamma_gph"], f["se_gph"],
               f["gamma_lw"], f["se_lw"]))
    gg = [f["gamma_gph"] for f in bar]
    gw = [f["gamma_lw"] for f in bar]
    mono = (all(np.diff(gg) > 0) and all(np.diff(gw) > 0))
    log("      recorrido: GPH %.4f -> %.4f (%.4f)   LW %.4f -> %.4f (%.4f)"
        % (gg[0], gg[-1], gg[-1] - gg[0], gw[0], gw[-1], gw[-1] - gw[0]))
    log("      P3 %s (monotona creciente en las dos)"
        % ("SE CUMPLE" if mono else "REFUTADA"))

    # -------------------------------------------------------------------
    # EL CONTROL QUE DECIDE COMO SE LEE EL FALLO DE P3.
    #
    # Que `gamma_hat` se mueva con la banda admite dos lecturas incompatibles:
    #   (i)  mis estimadores son sensibles a `m` y el numero no vale;
    #   (ii) los estimadores estan bien y la SERIE no tiene un regimen de
    #        escala unico, con lo que no existe «la» gamma que estimar.
    # Se separan corriendo EL MISMO BARRIDO sobre fGn de exponente CONOCIDO y
    # de la MISMA longitud. Si ahi el barrido sale plano, el fallo es del dato.
    # -------------------------------------------------------------------
    log("")
    log("  --- CONTROL DE P3: el mismo barrido sobre fGn de gamma CONOCIDA ---")
    H_v = 0.75                       # -> gamma verdadera = 2 - 2H = 0.50
    ctrl = T.barrido_banda(T._fgn(eps.size, H_v, np.random.default_rng(41)))
    log("      verdad: H = %.2f  ->  gamma = %.4f" % (H_v, 2 - 2 * H_v))
    log("      alpha      m     gamma GPH    gamma LW      error max")
    for f in ctrl:
        log("      %.1f   %7d   %+.4f     %+.4f     %.4f"
            % (f["alpha"], f["m"], f["gamma_gph"], f["gamma_lw"],
               max(abs(f["gamma_gph"] - 0.5), abs(f["gamma_lw"] - 0.5))))
    rec_ctrl = max(abs(np.diff([f["gamma_gph"] for f in ctrl])).max(),
                   abs(np.diff([f["gamma_lw"] for f in ctrl])).max())
    rec_real = max(abs(np.diff(gg)).max(), abs(np.diff(gw)).max())
    log("      salto maximo entre bandas:  fGn %.4f   REAL %.4f   razon %.1fx"
        % (rec_ctrl, rec_real, rec_real / max(rec_ctrl, 1e-9)))
    if rec_real > 3 * rec_ctrl:
        log("      ->  el barrido es PLANO sobre verdad conocida y se derrumba sobre")
        log("          el dato. El fallo de P3 es DE LA SERIE, no del estimador: no")
        log("          hay un regimen de escala unico, asi que no existe «la» gamma.")
    else:
        log("      ->  el control se mueve tanto como el dato: el barrido no")
        log("          discrimina y el fallo de P3 no es interpretable.")

    return {"A": a["gamma"], "B": esc["gamma_eff"], "gph": g["gamma"],
            "lw": w["gamma"], "C1": a["C1"], "barrido": bar, "ctrl": ctrl,
            "P1": p1, "P2": p2, "P3": mono,
            "el_dato_no_tiene_regimen": rec_real > 3 * rec_ctrl}


# ---------------------------------------------------------------------------
# Sec.3.2
# ---------------------------------------------------------------------------

def sec_3_2(s):
    titulo("Sec.3.2 -- N_eff con la definicion COHERENTE")
    log("`inflacion = sqrt(N/N_eff)` debe satisfacerse. `N^gamma_eff` descarta la")
    log("constante del escalado y no la satisface en ninguna fila.")
    log("")
    log("  serie                    N        infl   N_eff=N/infl^2   N^gamma (v3.3)")
    log("  --------------------------------------------------------------------------")
    out = {}
    for etiq, x in (("signos eps", s["eps"]),
                    ("incrementos de precio", np.diff(s["precio"]))):
        r = T.ic_memoria_larga(np.asarray(x, dtype=float))
        out[etiq] = r
        log("  %-22s %8d  %8.2f   %14.1f   %12.0f"
            % (etiq, r["n"], r["inflacion_vs_iid"], r["n_eff"],
               r["n_eff_potencia"]))
        if r["antipersistente"]:
            log("      [!] N_eff > N -> ANTIPERSISTENCIA, no «perdida». Es la firma")
            log("          del rebote bid-ask (rho1(retornos) = -0.216, v3.0), y la")
            log("          tabla de la v3.3 la escondia tras una columna de perdida.")
        else:
            log("      perdida efectiva %.1fx" % (r["n"] / max(r["n_eff"], 1.0)))
    log("  %-22s %8d  %8.2f   %14.1f   %12.0f    <- v3.3, ya publicada"
        % ("dLL de prueba", 199532, 3.15, 199532 / 3.15 ** 2, 8574))
    log("      [!] Esta fila NO se recalcula: sale por ARITMETICA sobre la inflacion")
    log("          3.15 que la v3.3 ya publico. El Sec.7 prohibe reabrir el conjunto")
    log("          de prueba y aqui no hace falta -- no se ajusta ni se elige nada.")
    log("          El paso 2 se recalculo con la inflacion, que es la columna")
    log("          correcta, asi que su veredicto no cambia.")
    log("")
    log("  [!] REGLA DEL PROYECTO, del propio Sec.3.2: `Var(media) ~ N^(-gamma)` es")
    log("      el resultado de Beran PARA LA MEDIA MUESTRAL. No se transporta a")
    log("      estimadores de exponentes, a pendientes de regresion ni a razones de")
    log("      verosimilitud. N_eff SE MIDE POR ESTADISTICO, midiendo el escalado de")
    log("      ese estadistico. No se importa de la serie de signos.")
    return out


# ---------------------------------------------------------------------------
# Sec.3.3
# ---------------------------------------------------------------------------

def sec_3_3(s):
    titulo("Sec.3.3 -- gamma y H en la MISMA ventana, y beta implicita como CURVA")
    log("[!] IDENTIDAD ESTRUCTURAL, que el Sec.3.1 obliga a escribir en el reporte:")
    log("")
    log("        beta_implicita = (2 - gamma)/2 - H = 0   <=>   gamma = 2 - 2H")
    log("")
    log("    que es la relacion entre exponente de ACF y exponente de Hurst del")
    log("    RUIDO GAUSSIANO FRACCIONARIO. `beta` implicita NO mide el nucleo: mide")
    log("    cuanto se desvia el par (signos, precio) de la relacion fGn. Si los dos")
    log("    comparten proceso de memoria larga, sale 0 MECANICAMENTE. Control 7 de")
    log("    `tick_grande.py` lo demuestra sobre fGn puro, donde no hay propagador.")
    log("")
    cur = T.curva_beta_implicita(s["eps"], s["precio"])
    log("  ventana         gamma  (R2)   pts    H     pts    beta     beta_dif  g=2-2H")
    log("  ------------------------------------------------------------------------------")
    for f in cur:
        log("  [%5d,%6d]  %+.4f (%.2f) %3d  %+.4f  %3d  %+.4f   %+.4f  %+.4f"
            % (f["lo"], f["hi"], f["gamma"], f["r2_gamma"], f["n_pts_gamma"],
               f["H"], f["n_pts_H"], f["beta"], f["beta_dif"],
               f["gamma_que_anula_beta"]))
    comun = [f for f in cur if (f["lo"], f["hi"]) == (256, 2000)][0]
    v33 = [f for f in cur if (f["lo"], f["hi"]) == (10, 2000)][0]
    v33H = [f for f in cur if (f["lo"], f["hi"]) == (256, 16384)][0]
    log("")
    log("  la v3.3 combino gamma de [10, 2000] con H de [256, 16384]:")
    log("     gamma = %+.4f (de la 1a)   H = %+.4f (de la 2a)   ->  beta = %+.4f"
        % (v33["gamma"], v33H["H"], T.beta_implicita(v33["gamma"], v33H["H"])))
    log("  en la ventana COMUN [256, 2000] las dos se ajustan juntas:")
    log("     gamma = %+.4f              H = %+.4f              ->  beta = %+.4f"
        % (comun["gamma"], comun["H"], comun["beta"]))

    betas = np.array([f["beta"] for f in cur])
    rec = float(betas.max() - betas.min())
    dist = float(np.mean([f["beta_dif"] for f in cur]))
    log("")
    log("  --- REGLA DE LECTURA DEL Sec.3.3, declarada antes ---")
    log("  recorrido de beta al variar la ventana : %.4f" % rec)
    log("  distancia entre las dos hipotesis      : %.4f  (0 contra beta_dif medio)"
        % dist)
    if rec > dist:
        log("  ->  beta implicita se mueve MAS que la distancia entre las hipotesis.")
        log("      EL Sec.1 DE LA v3.3 NO ESTA MIDIENDO NADA, y asi se reporta.")
    else:
        log("  ->  beta implicita se mueve MENOS que la distancia entre hipotesis;")
        log("      la ventana no la domina y el punto de la v3.3 es defendible.")
    return {"curva": cur, "recorrido": rec, "distancia": dist,
            "comun": comun, "domina_la_ventana": rec > dist}


# ---------------------------------------------------------------------------

def main() -> int:
    series = [serie_v33()]
    try:
        e = serie_estacional()
        series.append(e)
    except Exception as exc:
        log("[!] no se pudo cargar captura_estacional: %s: %s"
            % (type(exc).__name__, exc))

    for s in series:
        titulo("SERIE: %s" % s["nombre"])
        log("  N = %d ticks   nu = %.2f tx/s%s"
            % (s["eps"].size, s["nu"],
               "   duracion %.2f h" % s["horas"] if "horas" in s else ""))
        log("  precio %.1f - %.1f   mediana %.1f"
            % (s["precio"].min(), s["precio"].max(), np.median(s["precio"])))
        sec_3_1(s)
        sec_3_2(s)
        sec_3_3(s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
