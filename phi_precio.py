# -*- coding: utf-8 -*-
"""
phi_precio.py -- ¿que relacion hay entre phi' (ticks por volumen inyectado) y el precio?

QUE PREGUNTA RESPONDE
---------------------
El barrido SSA dejo un resultado incomodo: los autovectores de las ventanas
reales caen EXACTAMENTE en la escalera armonica de la ventana (T = 2L/k), con
error mediano 0.0000 y con el mismo T/L en ventanas de L entre 96 y 768 y con nu
entre 6 y 109 tx/s. Un patron identico en todas las ventanas es la firma de un
instrumento, no la de un mercado.

Este script mide la otra mitad. Si phi' y el precio se asocian de una forma que
CAMBIA de ventana en ventana mientras el periodo no cambia, entonces las dos
cosas estan separadas: la escalera es del aparato y la asociacion es del
mercado. Si la asociacion tambien sale constante, o si desaparece bajo el nulo,
entonces no hay tal separacion.

DEFINICION DE phi'
------------------
phi' = tasa de ticks por volumen inyectado [ticks/BTC]. Es la derivada del reloj
de transacciones respecto del reloj de volumen, los dos relojes que la v2.0
introdujo. Sobre un bloque de m ticks:

    phi'_b = m / sum(q_j)        j en el bloque b

Es el inverso del tamano medio de transaccion en ese bloque: phi' alta = muchas
operaciones pequenas; phi' baja = pocas y grandes. Se usa su logaritmo porque es
un cociente de colas pesadas (el volumen por transaccion tiene mediana 0.0060 BTC
y maximo 30.1420, medido en la v2.0).

BLOQUES NO SOLAPADOS
--------------------
Todas las correlaciones se estiman sobre bloques DISJUNTOS de m ticks. Con
ventanas solapadas la correlacion sale inflada por construccion, y este proyecto
ya lo pago tres veces: el filtro corrigiendo 90 veces con la misma medicion
(v1.3), R_n de ventanas EMD solapadas (v2.0) y la firma de volatilidad solapada
contra la no solapada discrepando EN SIGNO (sesion 2026-08-08 e).

EL NULO
-------
Se barajan los VOLUMENES dentro de la ventana dejando los precios intactos. Eso
conserva exactamente la marginal de q y toda la serie de precios, y destruye
solo el emparejamiento entre ambos -- que es justo la hipotesis a probar. Sin
este nulo, una correlacion no nula no significa nada: phi' y la volatilidad
comparten el conteo de ticks y se correlacionarian por aritmetica.

Uso:
  python phi_precio.py --captura=telemetria/captura_v31b
  python phi_precio.py --captura=telemetria/captura_larga --m=128
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np

import barrido_ssa as bs


RAIZ_LOG = "ssa_log"


def _comb(n: int, k: int) -> float:
    """Coeficiente binomial exacto, para la prueba de signos."""
    from math import comb
    return float(comb(n, k))


def bloques(obs: dict, m: int) -> dict:
    """Corta la ventana en bloques disjuntos de m ticks y calcula, por bloque,
    las magnitudes que se van a correlacionar."""
    p = obs["precio"]
    q = np.nan_to_num(obs["cant"], nan=0.0)
    vn = np.nan_to_num(obs["vol_neto"], nan=0.0)
    t = obs["t"]
    n = (len(p) // m) * m
    if n < 4 * m:
        return {}

    P = p[:n].reshape(-1, m)
    Q = q[:n].reshape(-1, m)
    VN = vn[:n].reshape(-1, m)
    T = t[:n].reshape(-1, m)

    vol_bloque = Q.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi = m / vol_bloque                      # [ticks/BTC]
    log_phi = np.log(np.where(vol_bloque > 0, phi, np.nan))

    # Cambio de precio DENTRO del bloque, y volatilidad realizada del bloque.
    dP = np.diff(P, axis=1)
    R = P[:, -1] - P[:, 0]
    V = np.sqrt((dP ** 2).sum(axis=1))            # volatilidad realizada [USD/BTC]
    dur = T[:, -1] - T[:, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        nu = np.where(dur > 0, m / dur, np.nan)   # tasa de ticks [tx/s]
    flujo = VN.sum(axis=1)                        # volumen neto firmado [BTC]

    ok = np.isfinite(log_phi) & np.isfinite(V) & np.isfinite(nu)
    return {
        "m": m,
        "n_bloques": int(np.count_nonzero(ok)),
        "log_phi": log_phi[ok],
        "phi": phi[ok],
        "precio": P[:, -1][ok],
        "R": R[ok],
        "absR": np.abs(R[ok]),
        "V": V[ok],
        "nu": nu[ok],
        "flujo": flujo[ok],
        "vol": vol_bloque[ok],
    }


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 4 or b.size < 4:
        return float("nan")
    a = a - np.mean(a)
    b = b - np.mean(b)
    da, db = np.linalg.norm(a), np.linalg.norm(b)
    if da == 0 or db == 0:
        return float("nan")
    return float(np.dot(a, b) / (da * db))


def _rango(x: np.ndarray) -> np.ndarray:
    """Rangos para Spearman. Con colas pesadas Pearson mide los tres outliers."""
    orden = np.argsort(np.argsort(x))
    return orden.astype(np.float64)


def relaciones(b: dict) -> dict:
    """Las correlaciones que interesan, cada una con su lectura.

    NOTA DE INTERPRETACION: `rho_nivel` (log phi' contra el NIVEL de precio) se
    calcula porque el operador pregunto por la asociacion con el precio, pero es
    una correlacion entre dos series no estacionarias y por tanto ESPURIA por
    construccion (Yule 1926). Se reporta marcada como tal y no se usa para
    concluir nada.
    """
    lp = b["log_phi"]
    # ⚠ Hay bloques con volatilidad EXACTAMENTE cero (64 ticks sin que el precio
    # se mueva): 2.82 % en captura_v31b y 0.62 % en captura_larga. Metidos en un
    # logaritmo se van a -27.6, unas 30 desviaciones fuera, y dominan cualquier
    # correlacion de Pearson ellos solos. Se excluyen de las correlaciones sobre
    # log(V) -- y solo de esas. Spearman no los necesita: el rango de un cero
    # esta perfectamente definido y es el mas bajo, que es justo lo que son.
    vpos = b["V"] > 0
    lpv = lp[vpos]
    lV = np.log(b["V"][vpos])
    return {
        "rho_nivel_ESPURIA": _corr(lp, b["precio"]),
        "rho_volatilidad": _corr(lpv, lV),
        "n_bloques_V_cero": int(np.count_nonzero(~vpos)),
        "rho_absR": _corr(lp, b["absR"]),
        "rho_retorno": _corr(lp, b["R"]),
        "rho_flujo": _corr(lp, b["flujo"]),
        "rho_nu": _corr(lp, np.log(b["nu"])),
        # Spearman de la relacion principal, por las colas pesadas del volumen.
        "spearman_volatilidad": _corr(_rango(lp), _rango(b["V"])),
        # Predictivo: phi' del bloque contra el retorno del bloque SIGUIENTE.
        "rho_predictivo_R": _corr(lp[:-1], b["R"][1:]),
        "rho_predictivo_absR": _corr(lp[:-1], b["absR"][1:]),
        "rho_predictivo_V": _corr(_rango(lp[:-1]), _rango(b["V"][1:])),
    }


def terciles(b: dict) -> dict:
    """Volatilidad realizada mediana en el tercil bajo y alto de phi'.

    Mas legible que una correlacion: dice CUANTO cambia la volatilidad entre
    mercado de operaciones grandes y mercado de operaciones pequenas, en las
    mismas unidades del precio.
    """
    lp = b["log_phi"]
    V = b["V"]
    if lp.size < 12:
        return {}
    q1, q2 = np.percentile(lp, [33.333, 66.667])
    bajo = V[lp <= q1]
    alto = V[lp >= q2]
    if bajo.size < 4 or alto.size < 4:
        return {}
    v_bajo = float(np.median(bajo))
    v_alto = float(np.median(alto))
    return {
        "V_phi_bajo": v_bajo,
        "V_phi_alto": v_alto,
        "razon_alto_bajo": float(v_alto / v_bajo) if v_bajo > 0 else float("nan"),
        "phi_bajo_med": float(np.median(np.exp(lp[lp <= q1]))),
        "phi_alto_med": float(np.median(np.exp(lp[lp >= q2]))),
    }


def nulo_por_desplazamiento(obs: dict, m: int, n_sorteos: int, rng) -> dict:
    """Desplaza CIRCULARMENTE la serie de volumenes contra la de precios.

    Es un nulo mas exigente que barajar: conserva intacta toda la estructura
    temporal del volumen -- incluido el agrupamiento de operaciones grandes, que
    el barajado destruye -- y rompe unicamente la alineacion entre volumen y
    precio, que es la hipotesis a probar. Barajar hace el nulo demasiado
    homogeneo y por tanto demasiado facil de superar.
    """
    claves = None
    acum = {}
    q0 = np.nan_to_num(obs["cant"], nan=0.0)
    vn0 = np.nan_to_num(obs["vol_neto"], nan=0.0)
    n = len(q0)
    # Se evitan los desplazamientos minusculos: un corrimiento de 3 ticks deja
    # los bloques practicamente alineados y el "nulo" seguiria conteniendo la
    # asociacion que se quiere destruir.
    minimo = max(4 * m, n // 20)
    for _ in range(n_sorteos):
        s = int(rng.integers(minimo, n - minimo))
        obs_n = dict(obs)
        obs_n["cant"] = np.roll(q0, s)
        obs_n["vol_neto"] = np.roll(vn0, s)
        bb = bloques(obs_n, m)
        if not bb:
            continue
        r = relaciones(bb)
        if claves is None:
            claves = list(r)
            acum = {k: [] for k in claves}
        for k in claves:
            acum[k].append(r[k])
    return {k: np.array(v) for k, v in acum.items()}


def nulo_por_barajado(obs: dict, m: int, n_sorteos: int, rng) -> dict:
    """Baraja los VOLUMENES dejando los precios intactos.

    Conserva la marginal exacta de q y la serie de precios entera; destruye solo
    el emparejamiento. Es el nulo correcto para "¿phi' se asocia al precio?".
    """
    claves = None
    acum = {}
    q0 = np.nan_to_num(obs["cant"], nan=0.0)
    vn0 = np.nan_to_num(obs["vol_neto"], nan=0.0)
    for _ in range(n_sorteos):
        perm = rng.permutation(len(q0))
        obs_n = dict(obs)
        obs_n["cant"] = q0[perm]
        obs_n["vol_neto"] = vn0[perm]
        bb = bloques(obs_n, m)
        if not bb:
            continue
        r = relaciones(bb)
        if claves is None:
            claves = list(r)
            acum = {k: [] for k in claves}
        for k in claves:
            acum[k].append(r[k])
    return {k: np.array(v) for k, v in acum.items()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--captura", default="telemetria/captura_v31b")
    ap.add_argument("--etiqueta", default=None)
    ap.add_argument("--n", type=int, default=8192, help="ticks por ventana")
    ap.add_argument("--ventanas", type=int, default=8)
    ap.add_argument("--m", type=int, default=64, help="ticks por bloque disjunto")
    ap.add_argument("--sorteos", type=int, default=300)
    ap.add_argument("--semilla", type=int, default=1)
    args = ap.parse_args()

    rng = np.random.default_rng(args.semilla)
    etiqueta = args.etiqueta or os.path.basename(os.path.normpath(args.captura))
    dir_log = os.path.join(RAIZ_LOG, "phi_%s" % etiqueta)
    os.makedirs(dir_log, exist_ok=True)

    print("=" * 78)
    print("phi' = ticks por volumen inyectado  CONTRA  precio")
    print("=" * 78)
    d = bs.cargar_captura(args.captura)
    N = len(d["tr_t"])
    dur = float(d["tr_t"][-1] - d["tr_t"][0])
    print("captura: %s | %d tx en %.2f h | nu = %.2f tx/s"
          % (args.captura, N, dur / 3600.0, N / dur))
    if not d["_tiene_maker"]:
        print("AVISO: sin `tr_maker`; rho_flujo no es interpretable.")
    print("ventanas: %d de %d ticks | bloques DISJUNTOS de %d ticks -> %d por ventana"
          % (args.ventanas, args.n, args.m, args.n // args.m))
    print("nulo: %d barajados del volumen dejando los precios intactos" % args.sorteos)
    print("")

    n_pos = N // args.n
    n_v = int(min(args.ventanas, n_pos))
    paso = max(1, n_pos // n_v)
    cortes = [(k * paso * args.n, k * paso * args.n + args.n) for k in range(n_v)]

    filas = []
    for k, (ini, fin) in enumerate(cortes):
        obs = bs.construir_observables(d, ini, fin)
        b = bloques(obs, args.m)
        if not b:
            continue
        r = relaciones(b)
        ter = terciles(b)
        nulos = {
            "desplazado": nulo_por_desplazamiento(obs, args.m, args.sorteos, rng),
            "barajado": nulo_por_barajado(obs, args.m, args.sorteos, rng),
        }
        pvals, z = {}, {}
        for nombre_n, nul in nulos.items():
            pvals[nombre_n], z[nombre_n] = {}, {}
            for clave, val in r.items():
                muestras = nul.get(clave, np.array([]))
                muestras = muestras[np.isfinite(muestras)]
                if muestras.size < 10 or not np.isfinite(val):
                    pvals[nombre_n][clave] = float("nan")
                    z[nombre_n][clave] = float("nan")
                    continue
                # p bilateral empirico y z contra la dispersion del nulo.
                pvals[nombre_n][clave] = float(np.mean(np.abs(muestras) >= abs(val)))
                s = float(np.std(muestras))
                z[nombre_n][clave] = (float((val - float(np.mean(muestras))) / s)
                                      if s > 0 else float("nan"))

        filas.append({
            "ventana": "v%02d" % k,
            "tick_ini": int(ini), "tick_fin": int(fin),
            "n_bloques": b["n_bloques"],
            "nu_tx_s": float(obs["nu_tx_s"]),
            "phi_mediana": float(np.median(b["phi"])),
            "phi_p10": float(np.percentile(b["phi"], 10)),
            "phi_p90": float(np.percentile(b["phi"], 90)),
            "rho": r, "p": pvals, "z": z, "terciles": ter,
        })

        print("-" * 78)
        print("v%02d | nu=%.1f tx/s | %d bloques | phi' med=%.1f tk/BTC "
              "(p10=%.1f p90=%.1f)"
              % (k, obs["nu_tx_s"], b["n_bloques"], np.median(b["phi"]),
                 np.percentile(b["phi"], 10), np.percentile(b["phi"], 90)))
        if ter:
            print("   terciles de phi': V mediana %.3f (phi'=%.1f) -> %.3f (phi'=%.1f)"
                  "  razon %.3f"
                  % (ter["V_phi_bajo"], ter["phi_bajo_med"],
                     ter["V_phi_alto"], ter["phi_alto_med"], ter["razon_alto_bajo"]))
        for clave in ("spearman_volatilidad", "rho_volatilidad", "rho_absR",
                      "rho_retorno", "rho_predictivo_V", "rho_nivel_ESPURIA"):
            print("   %-22s rho=%+.4f | desplazado p=%.4f z=%+6.2f | barajado p=%.4f"
                  % (clave, r[clave], pvals["desplazado"][clave],
                     z["desplazado"][clave], pvals["barajado"][clave]))

    # ---- lo que decide: ¿varia entre ventanas o no? ----
    print("")
    print("=" * 78)
    print("DISPERSION ENTRE VENTANAS -- la pregunta del operador")
    print("=" * 78)
    print("Si el periodo (T/L) es constante y esto NO lo es, son cosas distintas.")
    print("")
    print("%-24s %9s %9s %9s %9s %9s %8s %9s"
          % ("magnitud", "min", "mediana", "max", "rango", "sd", "signo", "p_signos"))
    claves = list(filas[0]["rho"]) if filas else []
    for clave in claves:
        v = np.array([f["rho"][clave] for f in filas], dtype=float)
        v = v[np.isfinite(v)]
        if v.size < 2:
            continue
        # Prueba de signos: bajo "no hay asociacion" cada ventana da signo + o -
        # con probabilidad 1/2. Es independiente de la magnitud, asi que
        # sobrevive aunque cada rho suelta sea pequena.
        n_pos = int(np.count_nonzero(v > 0))
        k_ext = max(n_pos, v.size - n_pos)
        p_signos = 2.0 * float(sum(_comb(v.size, i) for i in range(k_ext, v.size + 1))) / (2 ** v.size)
        p_signos = min(1.0, p_signos)
        print("%-24s %+9.4f %+9.4f %+9.4f %9.4f %9.4f %5d/%-2d %9.4f"
              % (clave, v.min(), np.median(v), v.max(), v.max() - v.min(), np.std(v),
                 n_pos, v.size, p_signos))

    razones = np.array([f["terciles"].get("razon_alto_bajo", np.nan)
                        for f in filas if f["terciles"]], dtype=float)
    razones = razones[np.isfinite(razones)]
    if razones.size:
        n_pos = int(np.count_nonzero(razones > 1.0))
        print("")
        print("terciles: V(phi' alto)/V(phi' bajo) = %.3f mediana "
              "(min %.3f, max %.3f); mayor que 1 en %d/%d ventanas"
              % (np.median(razones), razones.min(), razones.max(), n_pos, razones.size))

    phis = np.array([f["phi_mediana"] for f in filas])
    nus = np.array([f["nu_tx_s"] for f in filas])
    print("")
    print("%-24s %+9.4f %+9.4f %+9.4f %9.4f %9.4f"
          % ("phi' mediana [tk/BTC]", phis.min(), np.median(phis), phis.max(),
             phis.max() - phis.min(), np.std(phis)))
    print("%-24s %+9.4f %+9.4f %+9.4f %9.4f %9.4f"
          % ("nu [tx/s]", nus.min(), np.median(nus), nus.max(),
             nus.max() - nus.min(), np.std(nus)))

    # Comparacion con la escalera SSA, si hay log de la misma captura.
    ruta_ssa = os.path.join(RAIZ_LOG, etiqueta, "resumen.json")
    if os.path.exists(ruta_ssa):
        m_ssa = json.load(open(ruta_ssa, encoding="utf-8"))
        tsl = [v["observables"]["precio"]["T_sobre_L"]
               for v in m_ssa["ventanas"]
               if not v["ventana"].startswith("ctrl") and "precio" in v["observables"]]
        tsl = np.array([x for x in tsl if np.isfinite(x)])
        if tsl.size > 1:
            print("%-24s %+9.4f %+9.4f %+9.4f %9.4f %9.4f   <- del barrido SSA"
                  % ("T/L (periodo SSA)", tsl.min(), np.median(tsl), tsl.max(),
                     tsl.max() - tsl.min(), np.std(tsl)))

    # ---- figura ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        ncol = min(4, len(filas))
        nfil = int(np.ceil(len(filas) / ncol))
        fig, axs = plt.subplots(nfil + 1, ncol, figsize=(3.6 * ncol, 3.2 * (nfil + 1)),
                                squeeze=False)
        for k, (ini, fin) in enumerate(cortes[:len(filas)]):
            obs = bs.construir_observables(d, ini, fin)
            b = bloques(obs, args.m)
            ax = axs[k // ncol][k % ncol]
            vpos = b["V"] > 0
            ax.loglog(b["phi"][vpos], b["V"][vpos], ".", ms=3, alpha=0.55,
                      color="#1f77b4")
            n0 = int(np.count_nonzero(~vpos))
            if n0:
                # Los bloques de volatilidad cero existen y hay que verlos, pero
                # no pueden ir en un eje logaritmico: van al borde inferior.
                lo = float(np.min(b["V"][vpos])) * 0.5
                ax.loglog(b["phi"][~vpos], np.full(n0, lo), "x", ms=4,
                          color="#d62728", alpha=0.7, label="V=0 (%d)" % n0)
            t_ = filas[k]["terciles"]
            if t_:
                ax.plot([t_["phi_bajo_med"], t_["phi_alto_med"]],
                        [t_["V_phi_bajo"], t_["V_phi_alto"]], "o-", color="#d62728",
                        ms=8, lw=2, label="terciles x%.2f" % t_["razon_alto_bajo"])
                ax.legend(fontsize=7)
            ax.set_title("%s | rho_s=%+.3f (p=%.3f)"
                         % (filas[k]["ventana"],
                            filas[k]["rho"]["spearman_volatilidad"],
                            filas[k]["p"]["desplazado"]["spearman_volatilidad"]),
                         fontsize=9)
            ax.set_xlabel("phi' [ticks/BTC]")
            ax.set_ylabel("volatilidad realizada [USD/BTC]")
            ax.grid(alpha=0.3, which="both")
        for k in range(len(filas), nfil * ncol):
            axs[k // ncol][k % ncol].axis("off")

        # Fila final: la comparacion que decide.
        for j in range(ncol):
            axs[nfil][j].axis("off")
        axr = fig.add_subplot(nfil + 1, 1, nfil + 1)
        etiqs = [f["ventana"] for f in filas]
        xr = np.arange(len(etiqs))
        rho_s = [f["rho"]["spearman_volatilidad"] for f in filas]
        axr.bar(xr, rho_s, 0.6, color="#1f77b4", label="Spearman(log phi', V)")
        axr.axhline(0, color="k", lw=1.0)
        ruta_ssa2 = os.path.join(RAIZ_LOG, etiqueta, "resumen.json")
        if os.path.exists(ruta_ssa2):
            m2 = json.load(open(ruta_ssa2, encoding="utf-8"))
            tsl2 = [v["observables"]["precio"]["T_sobre_L"] for v in m2["ventanas"]
                    if not v["ventana"].startswith("ctrl") and "precio" in v["observables"]]
            tsl2 = [x for x in tsl2 if np.isfinite(x)]
            if len(tsl2) == len(etiqs):
                axr.plot(xr, tsl2, "s--", color="#d62728", ms=7,
                         label="T/L del SSA (constante = del instrumento)")
        axr.set_xticks(xr)
        axr.set_xticklabels(etiqs, fontsize=8)
        axr.set_ylabel("valor")
        axr.set_title("Lo que VARIA entre ventanas (phi'-volatilidad) contra lo que NO (T/L)",
                      fontsize=10)
        axr.legend(fontsize=8)
        axr.grid(alpha=0.3, axis="y")

        fig.tight_layout()
        ruta_fig = os.path.join(dir_log, "phi_precio.png")
        fig.savefig(ruta_fig, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print("figura en %s" % ruta_fig)
    except Exception as exc:  # una figura fallida no debe tirar el analisis
        print("no se pudo dibujar: %s" % exc)

    with open(os.path.join(dir_log, "phi_precio.json"), "w", encoding="utf-8") as fh:
        json.dump(bs._limpiar({
            "captura": args.captura, "m": args.m, "n": args.n,
            "sorteos": args.sorteos, "ventanas": filas,
        }), fh, indent=2, ensure_ascii=False)
    print("")
    print("log en %s" % os.path.join(dir_log, "phi_precio.json"))


if __name__ == "__main__":
    main()
