# -*- coding: utf-8 -*-
"""
barrido_ssa.py -- TOMA DE DATOS para el analisis SSA.

Este script no decide nada sobre el modelo. Su unico trabajo es producir un log
completo y reutilizable: por cada ventana de muestreo guarda los datos crudos
(tiempo, precio, id de tick, cantidad y volumen neto firmado) y, por cada
longitud de ventana L del barrido, los autovalores, los autovectores y la matriz
de w-correlacion. Todo queda en disco para graficarlo despues sin volver a
calcular.

QUE SE MIDE POR VENTANA
-----------------------
  - ortogonalidad de la descomposicion (media y maximo de |w-correlacion| fuera
    de la diagonal), que es lo que dice cuan limpia es esa ventana;
  - la misma metrica sobre sustitutos barajados, como linea base;
  - la L elegida por MINIMO LOCAL -- nunca por cero absoluto, que es inalcanzable
    porque el ruido de medicion siempre filtra energia entre componentes;
  - Monte Carlo SSA contra un nulo AR(1): que componentes no explica el ruido rojo;
  - color de ruido (blanco / rosa / rojo / negro) del nivel, de los incrementos y
    del residuo que deja el SSA.

OBSERVABLES
-----------
  precio           nivel de precio de transaccion (reloj de ticks, Delta n = 1)
  vol              cantidad por transaccion (BTC)
  vol_neto         volumen firmado: eps * q, con eps = -1 si el comprador es maker
  vol_neto_acum    desequilibrio de flujo acumulado (suma de vol_neto)

El volumen neto necesita el campo `m` de Binance, que solo persisten las capturas
posteriores a la v3.1 (`tr_maker`). Si no esta, se avisa y se sigue sin el en vez
de fabricar un signo.

CONTROLES
---------
Con `--controles` se anaden al MISMO log tres ventanas sinteticas, para que en las
graficas aparezcan al lado de las reales:
  ctrl_positivo : dos periodos conocidos (120 y 37 ticks) sobre tendencia y ruido
  ctrl_paseo    : paseo aleatorio con la sigma de los datos -- lo que hundio a la EMD
  ctrl_barajado : incrementos reales barajados (conserva la marginal exacta)

Uso:
  python barrido_ssa.py --captura=telemetria/captura_v31b --controles
  python barrido_ssa.py --captura=telemetria/captura_larga --n=16384 --ventanas=6
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import time

import numpy as np

import ssa


RAIZ_LOG = "ssa_log"

# Rejilla de L por omision. Geometrica: la resolucion en frecuencia va como 1/L,
# asi que pasos multiplicativos cubren el rango con muestreo uniforme en log.
LS_POR_OMISION = [64, 96, 128, 192, 256, 384, 512, 768, 1024]


# ---------------------------------------------------------------------------
# Carga
# ---------------------------------------------------------------------------

def cargar_captura(ruta: str) -> dict:
    """Lee un .npz suelto o un directorio de bloques y los concatena en orden."""
    if os.path.isdir(ruta):
        ficheros = sorted(glob.glob(os.path.join(ruta, "*.npz")))
    else:
        ficheros = [ruta]
    if not ficheros:
        raise SystemExit("no hay ficheros en %s" % ruta)

    campos = {}
    tiene_maker = True
    for f in ficheros:
        d = np.load(f)
        if "tr_maker" not in d.files:
            tiene_maker = False
        for k in ("tr_t", "tr_precio", "tr_id", "tr_cant", "tr_maker"):
            if k in d.files:
                campos.setdefault(k, []).append(d[k])

    out = {k: np.concatenate(v) for k, v in campos.items()}
    if not tiene_maker:
        out.pop("tr_maker", None)

    # Orden temporal y filtro de los ceros del feed (~0.2 %): un precio de 0
    # USD/BTC no es una transaccion, y hasta la v3.0 se colaban sin filtrar.
    orden = np.argsort(out["tr_t"], kind="stable")
    for k in list(out):
        out[k] = out[k][orden]
    valido = (out["tr_precio"] > 0)
    if "tr_cant" in out:
        valido &= (out["tr_cant"] > 0)
    n_ceros = int(np.count_nonzero(~valido))
    for k in list(out):
        out[k] = out[k][valido]

    out["_n_ceros_filtrados"] = n_ceros
    out["_n_ficheros"] = len(ficheros)
    out["_tiene_maker"] = tiene_maker
    return out


def construir_observables(d: dict, ini: int, fin: int) -> dict:
    """Recorta una ventana y arma los observables. Todo en reloj de ticks."""
    t = d["tr_t"][ini:fin]
    p = d["tr_precio"][ini:fin]
    ident = d["tr_id"][ini:fin] if "tr_id" in d else np.arange(ini, fin)
    q = d["tr_cant"][ini:fin] if "tr_cant" in d else np.zeros(fin - ini)

    if "tr_maker" in d:
        # Convencion del proyecto (PREREGISTRO_3_1 §4): m = True -> el comprador
        # es maker -> el agresor es el vendedor -> eps = -1.
        maker = d["tr_maker"][ini:fin].astype(bool)
        eps = np.where(maker, -1.0, 1.0)
    else:
        eps = np.full(fin - ini, np.nan)

    vol_neto = eps * q
    vol_neto_acum = np.cumsum(np.nan_to_num(vol_neto, nan=0.0))

    # Retornos. Hacen falta como observable propio: sobre un NIVEL casi de raiz
    # unitaria el primer autovector se lleva ~98 % de la energia por
    # construccion (es la tendencia), y el resto de la descomposicion queda
    # reducida a los armonicos de la ventana. Si hay ciclo, en retornos tiene
    # sitio para aparecer.
    retorno = np.concatenate(([0.0], np.diff(p)))
    log_retorno = np.concatenate(([0.0], np.diff(np.log(p))))

    dur = float(t[-1] - t[0]) if len(t) > 1 else 0.0
    dt_seg = dur / max(len(t) - 1, 1)

    return {
        "t": t, "precio": p, "id": ident, "cant": q, "signo": eps,
        "vol_neto": vol_neto, "vol_neto_acum": vol_neto_acum,
        "vol": q, "retorno": retorno, "log_retorno": log_retorno,
        "dt_seg": dt_seg, "dur_seg": dur,
        "nu_tx_s": (len(t) / dur) if dur > 0 else float("nan"),
    }


def serie_observable(obs: dict, nombre: str) -> np.ndarray:
    if nombre not in ("precio", "vol", "vol_neto", "vol_neto_acum",
                      "retorno", "log_retorno"):
        raise SystemExit("observable desconocido: %s" % nombre)
    x = np.asarray(obs[nombre], dtype=np.float64)
    if not np.all(np.isfinite(x)):
        return np.array([])
    return x


# ---------------------------------------------------------------------------
# Analisis de una ventana
# ---------------------------------------------------------------------------

def analizar_serie(x: np.ndarray, Ls, d: int, d_eval: int, n_sus: int,
                   n_mc: int, dt_seg: float, semilla: int, nulos) -> dict:
    """Barrido de L + analisis completo en la L elegida. Devuelve todo lo que
    haya que guardar, sin imprimir ni decidir."""
    res = {}

    # --- barrido, guardando la descomposicion entera de cada L ---
    bar = ssa.barrer_L(x, Ls, d=d, d_eval=d_eval, n_sustitutos=n_sus, semilla=semilla)
    res["barrido"] = bar

    por_L = {}
    for L in Ls:
        des = ssa.descomponer(x, L, d=d)
        W = ssa.matriz_wcorrelacion(des["elementales"], des["cuentas"])
        por_L[L] = {
            "lambdas": des["lambdas"],       # espectro COMPLETO (L valores)
            "U": des["U"],                   # L x d
            "wcorr": W,                      # d x d
        }
    res["por_L"] = por_L

    # --- eleccion por minimo local ---
    eleccion = ssa.elegir_minimo_local(bar["Ls"], bar["media_abs"])
    L_est = eleccion["L"] if eleccion["hay_minimo_local"] else eleccion["L_argmin_global"]
    eleccion["L_usada"] = int(L_est)
    res["eleccion"] = eleccion

    # --- analisis completo en L_est ---
    des = ssa.descomponer(x, L_est, d=d)
    W = ssa.matriz_wcorrelacion(des["elementales"], des["cuentas"])
    res["elegido"] = {
        "L": int(L_est),
        "lambdas": des["lambdas"],
        "U": des["U"],
        "V": des["V"],
        "elementales": des["elementales"],
        "wcorr": W,
        "residuo": des["residuo"],
        "media": des["media"],
    }
    res["metrica_elegido"] = ssa.metrica_ortogonalidad(W, d_eval)
    res["pares"] = ssa.detectar_pares(des, dt_seg=dt_seg)
    res["escalera"] = ssa.escalera_de_ventana(des)
    # Un MC-SSA por cada nulo. Un solo nulo no basta: AR(1) sobre el nivel no
    # reproduce el rebote bid-ask y marca como significativa la microestructura.
    res["mc"] = {}
    for j, nulo in enumerate(nulos):
        res["mc"][nulo] = ssa.monte_carlo_ssa(x, L_est, d=d, n_sustitutos=n_mc,
                                              nulo=nulo, semilla=semilla + 1 + j)
    # El nulo principal es el mas exigente disponible: el que reproduce el rebote.
    res["nulo_principal"] = ("ar1_incrementos" if "ar1_incrementos" in res["mc"]
                             else list(res["mc"])[0])

    # --- color de ruido: nivel, incrementos y residuo del SSA ---
    dx = np.diff(x)
    res["color"] = {
        "nivel": ssa.color_de_ruido(x, dt_seg=dt_seg),
        "incrementos": ssa.color_de_ruido(dx, dt_seg=dt_seg),
        "residuo_ssa": ssa.color_de_ruido(des["residuo"], dt_seg=dt_seg),
    }
    return res


# ---------------------------------------------------------------------------
# Escritura del log
# ---------------------------------------------------------------------------

def _limpiar(v):
    """Convierte a tipos serializables en JSON."""
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, dict):
        return {k: _limpiar(w) for k, w in v.items()}
    if isinstance(v, (list, tuple)):
        return [_limpiar(w) for w in v]
    return v


def guardar_ventana(dir_v: str, obs: dict, resultados: dict) -> None:
    os.makedirs(dir_v, exist_ok=True)

    # Datos crudos: tiempo, precio, tick y volumen. Es lo que pidio el operador y
    # lo que permite rehacer cualquier analisis sin volver a la captura.
    np.savez_compressed(
        os.path.join(dir_v, "datos.npz"),
        t=obs["t"], precio=obs["precio"], id=obs["id"], cant=obs["cant"],
        signo=obs["signo"], vol_neto=obs["vol_neto"],
        vol_neto_acum=obs["vol_neto_acum"],
        dt_seg=np.array([obs["dt_seg"]]), nu_tx_s=np.array([obs["nu_tx_s"]]),
    )

    for nombre, r in resultados.items():
        bar = r["barrido"]
        paquete = {
            "Ls": bar["Ls"],
            "media_abs": bar["media_abs"],
            "max_abs": bar["max_abs"],
            "p90_abs": bar["p90_abs"],
            "media_abs_sustituto": bar["media_abs_sustituto"],
            "sd_sustituto": bar["sd_sustituto"],
        }
        for L, blo in r["por_L"].items():
            paquete["L%04d_lambdas" % L] = blo["lambdas"]
            paquete["L%04d_U" % L] = blo["U"]
            paquete["L%04d_wcorr" % L] = blo["wcorr"]
        np.savez_compressed(os.path.join(dir_v, "barrido_%s.npz" % nombre), **paquete)

        el = r["elegido"]
        col = r["color"]
        extra = {}
        for nulo, mc in r["mc"].items():
            extra["mc_%s_p5" % nulo] = mc["p5"]
            extra["mc_%s_p95" % nulo] = mc["p95"]
            extra["mc_%s_pvalor" % nulo] = mc["pvalor"]
            extra["mc_%s_sobre" % nulo] = mc["sobre_p95"]
        mc0 = r["mc"][r["nulo_principal"]]
        np.savez_compressed(
            os.path.join(dir_v, "elegido_%s.npz" % nombre),
            L=np.array([el["L"]]),
            lambdas=el["lambdas"], U=el["U"], V=el["V"],
            elementales=el["elementales"], wcorr=el["wcorr"],
            residuo=el["residuo"], media=np.array([el["media"]]),
            mc_lambda_obs=mc0["lambda_obs"], mc_p5=mc0["p5"], mc_p95=mc0["p95"],
            mc_pvalor=mc0["pvalor"], mc_sobre=mc0["sobre_p95"],
            **extra,
            col_nivel_f=col["nivel"].get("f", np.array([])),
            col_nivel_S=col["nivel"].get("S", np.array([])),
            col_incr_f=col["incrementos"].get("f", np.array([])),
            col_incr_S=col["incrementos"].get("S", np.array([])),
            col_res_f=col["residuo_ssa"].get("f", np.array([])),
            col_res_S=col["residuo_ssa"].get("S", np.array([])),
        )


def resumen_de_ventana(etiqueta: str, obs: dict, resultados: dict, ini: int, fin: int) -> dict:
    r = {
        "ventana": etiqueta,
        "tick_ini": int(ini), "tick_fin": int(fin), "n": int(fin - ini),
        "t_ini": float(obs["t"][0]), "t_fin": float(obs["t"][-1]),
        "dur_seg": float(obs["dur_seg"]), "dt_seg": float(obs["dt_seg"]),
        "nu_tx_s": float(obs["nu_tx_s"]),
        "precio_ini": float(obs["precio"][0]), "precio_fin": float(obs["precio"][-1]),
        # Hueco temporal maximo dentro de la ventana. La captura larga de la v3.0
        # tuvo un corte de DNS de 10 h; una ventana que lo cruce mezcla dos tramos
        # y su espectro no significa nada.
        "gap_max_seg": float(np.max(np.diff(obs["t"]))) if len(obs["t"]) > 1 else 0.0,
        "vol_total_btc": float(np.nansum(obs["cant"])),
        "vol_neto_btc": float(np.nansum(obs["vol_neto"])),
        "observables": {},
    }
    for nombre, res in resultados.items():
        bar = res["barrido"]
        i_el = int(np.where(bar["Ls"] == res["eleccion"]["L_usada"])[0][0])
        pares = res["pares"]
        L_el = int(res["eleccion"]["L_usada"])
        # T/L del par dominante. Si el "ciclo" lo fija la ventana y no la senal,
        # esta razon es estable entre ventanas de L muy distinta -- que es
        # exactamente el sintoma con el que la v2.2 desenmascaro a la EMD.
        T_dom = pares[0]["periodo_muestras"] if pares else float("nan")
        r["observables"][nombre] = {
            "L_elegida": L_el,
            "T_dominante_ticks": float(T_dom),
            "T_dominante_seg": float(pares[0]["periodo_seg"]) if pares else float("nan"),
            "T_sobre_L": float(T_dom / L_el) if pares else float("nan"),
            "energia_par_dominante": float(pares[0]["energia_frac"]) if pares else 0.0,
            "escalera_error_mediano": float(res["escalera"]["error_mediano"]),
            "escalera_T_sobre_L": _limpiar(res["escalera"]["T_sobre_L"][:8]),
            "energia_primer_eof": float(res["escalera"]["energia"][0]),
            "mc_por_nulo": {k: {"n_sobre": int(v["n_sobre"]),
                                "frac_sobre": float(v["frac_sobre"])}
                            for k, v in res["mc"].items()},
            "nulo_principal": res["nulo_principal"],
            "hay_minimo_local": bool(res["eleccion"]["hay_minimo_local"]),
            "motivo_eleccion": res["eleccion"]["motivo"],
            "argmin_en_borde": bool(res["eleccion"]["en_borde"]),
            "ortogonalidad_media": float(res["metrica_elegido"]["media_abs"]),
            "ortogonalidad_max": float(res["metrica_elegido"]["max_abs"]),
            "ortogonalidad_sustituto": float(bar["media_abs_sustituto"][i_el]),
            "razon_real_sustituto": float(bar["media_abs"][i_el] /
                                          max(bar["media_abs_sustituto"][i_el], 1e-12)),
            "n_pares": len(pares),
            "pares": _limpiar(pares[:5]),
            "mc_n_sobre_p95": int(res["mc"][res["nulo_principal"]]["n_sobre"]),
            "mc_frac_sobre": float(res["mc"][res["nulo_principal"]]["frac_sobre"]),
            "mc_ar1_a": float(res["mc"][res["nulo_principal"]]["ar1_a"]),
            "color_nivel": _limpiar({k: v for k, v in res["color"]["nivel"].items()
                                     if k not in ("f", "S")}),
            "color_incrementos": _limpiar({k: v for k, v in res["color"]["incrementos"].items()
                                           if k not in ("f", "S")}),
            "color_residuo_ssa": _limpiar({k: v for k, v in res["color"]["residuo_ssa"].items()
                                           if k not in ("f", "S")}),
        }
    return r


# ---------------------------------------------------------------------------
# Ventanas de control
# ---------------------------------------------------------------------------

def ventanas_de_control(obs_ref: dict, n: int, semilla: int) -> list:
    """Tres ventanas sinteticas con verdad conocida, en el mismo formato que las
    reales para que atraviesen exactamente el mismo codigo."""
    rng = np.random.default_rng(semilla)
    t = obs_ref["t"][:n].copy()
    dt = obs_ref["dt_seg"]
    p = obs_ref["precio"][:n]
    sigma = float(np.std(np.diff(p), ddof=1))
    ident = np.arange(n)
    vol_neto_real = np.nan_to_num(obs_ref["vol_neto"][:n], nan=0.0)
    q_real = np.nan_to_num(obs_ref["cant"][:n], nan=0.0)

    def envolver(precio, vol_neto, cant, etiqueta):
        # Cada control lleva su PROPIO flujo. Compartirlo haria que la fila de
        # vol_neto_acum saliera identica en los tres y no midiera nada.
        return etiqueta, {
            "t": t, "precio": precio, "id": ident,
            "cant": cant, "signo": np.sign(vol_neto),
            "vol_neto": vol_neto, "vol_neto_acum": np.cumsum(vol_neto), "vol": cant,
            "retorno": np.concatenate(([0.0], np.diff(precio))),
            "log_retorno": np.concatenate(([0.0], np.diff(np.log(precio)))),
            "dt_seg": dt, "dur_seg": float(t[-1] - t[0]),
            "nu_tx_s": obs_ref["nu_tx_s"],
        }

    k = np.arange(n)
    base = float(np.mean(p))

    # POSITIVO: dos periodos conocidos, 120 y 37 ticks, sobre tendencia y ruido.
    # El flujo lleva el MISMO ciclo de 120 ticks, para que vol_neto_acum tambien
    # tenga una verdad conocida que recuperar.
    positivo = (base
                + 40.0 * np.sin(2 * np.pi * k / 120.0 + 0.4)
                + 15.0 * np.sin(2 * np.pi * k / 37.0)
                + 0.002 * k
                + rng.normal(0, 4.0 * sigma, n))
    q_pos = np.abs(rng.normal(0.01, 0.005, n))
    vn_pos = q_pos * (0.8 * np.sin(2 * np.pi * k / 120.0) + rng.normal(0, 0.5, n))

    # NEGATIVO: paseo aleatorio con la sigma real, y flujo de signos iid.
    paseo = base + np.cumsum(rng.normal(0, sigma, n))
    q_pa = np.abs(rng.normal(0.01, 0.005, n))
    vn_pa = q_pa * rng.choice([-1.0, 1.0], size=n)

    # NEGATIVO: incrementos reales barajados y volumen neto real barajado. Ambos
    # conservan la marginal exacta y destruyen solo el orden temporal.
    barajado = ssa.sustituto_barajado(p, rng)
    vn_bar = rng.permutation(vol_neto_real)

    return [envolver(positivo, vn_pos, q_pos, "ctrl_positivo"),
            envolver(paseo, vn_pa, q_pa, "ctrl_paseo"),
            envolver(barajado, vn_bar, q_real, "ctrl_barajado")]


# ---------------------------------------------------------------------------
# Principal
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--captura", default="telemetria/captura_v31b")
    ap.add_argument("--etiqueta", default=None)
    ap.add_argument("--n", type=int, default=8192, help="ticks por ventana")
    ap.add_argument("--ventanas", type=int, default=8)
    ap.add_argument("--obs", default="precio,vol_neto_acum")
    ap.add_argument("--Ls", default=",".join(str(v) for v in LS_POR_OMISION))
    ap.add_argument("--d", type=int, default=30, help="componentes conservadas")
    ap.add_argument("--d-eval", type=int, default=10,
                    help="componentes sobre las que se evalua la ortogonalidad")
    ap.add_argument("--sustitutos", type=int, default=3)
    ap.add_argument("--mc", type=int, default=100, help="sustitutos de Monte Carlo SSA")
    ap.add_argument("--nulos", default="ar1_incrementos,ar1,barajado",
                    help="nulos de MC-SSA; el principal es ar1_incrementos si esta")
    ap.add_argument("--controles", action="store_true")
    ap.add_argument("--semilla", type=int, default=1)
    args = ap.parse_args()

    Ls = [int(v) for v in args.Ls.split(",") if v.strip()]
    observables = [v.strip() for v in args.obs.split(",") if v.strip()]
    nulos = [v.strip() for v in args.nulos.split(",") if v.strip()]
    etiqueta = args.etiqueta or os.path.basename(os.path.normpath(args.captura))
    dir_log = os.path.join(RAIZ_LOG, etiqueta)
    os.makedirs(dir_log, exist_ok=True)

    print("=" * 74)
    print("BARRIDO SSA -- toma de datos")
    print("=" * 74)
    d = cargar_captura(args.captura)
    N = len(d["tr_t"])
    dur = float(d["tr_t"][-1] - d["tr_t"][0])
    print("captura      : %s (%d bloque(s))" % (args.captura, d["_n_ficheros"]))
    print("transacciones: %d en %.2f h  ->  nu = %.2f tx/s"
          % (N, dur / 3600.0, N / dur if dur > 0 else float("nan")))
    print("ceros del feed filtrados: %d" % d["_n_ceros_filtrados"])
    if not d["_tiene_maker"]:
        print("AVISO: la captura NO trae `tr_maker`; el volumen neto no se puede firmar.")
        observables = [o for o in observables if "neto" not in o]
        if not observables:
            observables = ["precio"]
        print("       observables reducidos a: %s" % ", ".join(observables))
    print("observables  : %s" % ", ".join(observables))
    print("rejilla de L : %s" % Ls)
    print("ventanas     : %d de %d ticks (d=%d, d_eval=%d)"
          % (args.ventanas, args.n, args.d, args.d_eval))
    print("log          : %s" % dir_log)
    print("")

    if max(Ls) * 2 > args.n:
        raise SystemExit("L maxima (%d) necesita al menos 2*L = %d ticks por ventana"
                         % (max(Ls), 2 * max(Ls)))

    # Ventanas NO solapadas, repartidas por toda la captura.
    n_pos = N // args.n
    n_ventanas = int(min(args.ventanas, n_pos))
    if n_ventanas < 1:
        raise SystemExit("la captura no da ni para una ventana de %d ticks" % args.n)
    paso = n_pos // n_ventanas
    cortes = [(k * paso * args.n, k * paso * args.n + args.n) for k in range(n_ventanas)]

    trabajo = []
    for k, (ini, fin) in enumerate(cortes):
        trabajo.append(("v%02d" % k, construir_observables(d, ini, fin), ini, fin))

    if args.controles:
        ref = trabajo[0][1]
        for nombre, obs_c in ventanas_de_control(ref, args.n, args.semilla):
            trabajo.append((nombre, obs_c, -1, -1))
        print("controles anadidos: ctrl_positivo, ctrl_paseo, ctrl_barajado")
        print("")

    resumenes = []
    t0 = time.time()
    for etiq_v, obs, ini, fin in trabajo:
        print("-" * 74)
        print("VENTANA %s | n=%d | dur=%.1f s | nu=%.2f tx/s | dt=%.4f s/tick"
              % (etiq_v, len(obs["t"]), obs["dur_seg"], obs["nu_tx_s"], obs["dt_seg"]))
        print("   precio %.2f -> %.2f | volumen %.4f BTC | volumen neto %+.4f BTC"
              % (obs["precio"][0], obs["precio"][-1],
                 float(np.nansum(obs["cant"])), float(np.nansum(obs["vol_neto"]))))
        hueco = float(np.max(np.diff(obs["t"]))) if len(obs["t"]) > 1 else 0.0
        if hueco > 60.0:
            print("   AVISO: hueco temporal de %.0f s dentro de la ventana; "
                  "el espectro mezcla dos tramos" % hueco)

        resultados = {}
        for nombre in observables:
            x = serie_observable(obs, nombre)
            if x.size == 0:
                print("   %-14s: sin datos utilizables, se omite" % nombre)
                continue
            ta = time.time()
            r = analizar_serie(x, Ls, args.d, args.d_eval, args.sustitutos,
                               args.mc, obs["dt_seg"], args.semilla, nulos)
            resultados[nombre] = r

            bar = r["barrido"]
            el = r["eleccion"]
            i_el = int(np.where(bar["Ls"] == el["L_usada"])[0][0])
            m = r["metrica_elegido"]
            pares = r["pares"]
            mc = r["mc"]
            col = r["color"]

            print("   %-14s: L=%d (%s) | media|wcorr|=%.4f  max=%.4f"
                  % (nombre, el["L_usada"],
                     "minimo local" if el["hay_minimo_local"] else "SIN minimo local",
                     m["media_abs"], m["max_abs"]))
            print("   %-14s  sustituto barajado=%.4f  ->  razon real/nulo = %.3f"
                  % ("", bar["media_abs_sustituto"][i_el],
                     bar["media_abs"][i_el] / max(bar["media_abs_sustituto"][i_el], 1e-12)))
            if pares:
                p0 = pares[0]
                print("   %-14s  %d par(es); dominante T=%.1f ticks = %.1f s "
                      "(%.2f%% energia, T/L=%.3f)"
                      % ("", len(pares), p0["periodo_muestras"], p0["periodo_seg"],
                         100 * p0["energia_frac"],
                         p0["periodo_muestras"] / el["L_usada"]))
            else:
                print("   %-14s  ningun par oscilatorio" % "")
            for nulo, m_ in mc.items():
                print("   %-14s  MC-SSA vs %-16s: %2d/%d sobre p95 (%.1f%%, azar 5%%)"
                      % ("", nulo, m_["n_sobre"], m_["d"], 100 * m_["frac_sobre"]))
            cn, ci, cr = col["nivel"], col["incrementos"], col["residuo_ssa"]
            print("   %-14s  color: nivel beta=%+.2f (%s) | incr beta=%+.2f (%s) | resid beta=%+.2f (%s)"
                  % ("", cn.get("beta", float("nan")), cn.get("color", "?"),
                     ci.get("beta", float("nan")), ci.get("color", "?"),
                     cr.get("beta", float("nan")), cr.get("color", "?")))
            print("   %-14s  [%.1f s]" % ("", time.time() - ta))

        dir_v = os.path.join(dir_log, etiq_v)
        guardar_ventana(dir_v, obs, resultados)
        resumenes.append(resumen_de_ventana(etiq_v, obs, resultados, ini, fin))

    meta = {
        "captura": args.captura,
        "etiqueta": etiqueta,
        "n_transacciones": int(N),
        "duracion_seg": dur,
        "nu_tx_s": (N / dur) if dur > 0 else None,
        "ceros_filtrados": int(d["_n_ceros_filtrados"]),
        "tiene_maker": bool(d["_tiene_maker"]),
        "observables": observables,
        "Ls": Ls,
        "n_por_ventana": int(args.n),
        "d": int(args.d), "d_eval": int(args.d_eval),
        "sustitutos_barrido": int(args.sustitutos),
        "sustitutos_mc": int(args.mc),
        "nulos": nulos,
        "semilla": int(args.semilla),
        "controles": bool(args.controles),
        "ventanas": resumenes,
        "segundos_totales": time.time() - t0,
    }
    with open(os.path.join(dir_log, "resumen.json"), "w", encoding="utf-8") as fh:
        json.dump(_limpiar(meta), fh, indent=2, ensure_ascii=False)

    # ---- tabla final ----
    print("")
    print("=" * 74)
    print("ORTOGONALIDAD POR VENTANA (menor = separacion mas limpia)")
    print("=" * 74)
    for nombre in observables:
        filas = [(r["ventana"], r["observables"][nombre])
                 for r in resumenes if nombre in r["observables"]]
        if not filas:
            continue
        print("")
        print("observable: %s" % nombre)
        print("%-14s %5s %8s %6s %7s %7s %7s %7s %s"
              % ("ventana", "L*", "media", "razon", "T/L", "E_par%", "E_EOF1%",
                 "escal.", "MC>p95 por nulo"))
        for etiq_v, o in filas:
            mcs = " ".join("%s=%d" % (k[:6], v["n_sobre"])
                           for k, v in o["mc_por_nulo"].items())
            print("%-14s %5d %8.4f %6.3f %7.3f %7.2f %7.2f %7.4f %s"
                  % (etiq_v, o["L_elegida"], o["ortogonalidad_media"],
                     o["razon_real_sustituto"], o["T_sobre_L"],
                     100 * o["energia_par_dominante"],
                     100 * o["energia_primer_eof"],
                     o["escalera_error_mediano"], mcs))
        reales = [o for e, o in filas if not e.startswith("ctrl")]
        if reales:
            med = np.median([o["ortogonalidad_media"] for o in reales])
            tsl = np.array([o["T_sobre_L"] for o in reales], dtype=float)
            tsl = tsl[np.isfinite(tsl)]
            esc = np.median([o["escalera_error_mediano"] for o in reales])
            print("%-14s %5s %8.4f %6s %7.3f %7s %7s %7.4f  <- medianas reales"
                  % ("", "", med, "", np.median(tsl) if tsl.size else float("nan"),
                     "", "", esc))
        print("   escal. = error mediano al escalon 2L/k. Cerca de 0 = los "
              "autovectores son armonicos de la VENTANA, no del mercado.")

    print("")
    print("log completo en %s (%.1f s)" % (dir_log, meta["segundos_totales"]))
    print("graficar con: python graficar_ssa.py --log=%s" % dir_log)


if __name__ == "__main__":
    main()
