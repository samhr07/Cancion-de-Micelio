# -*- coding: utf-8 -*-
"""
graficar_ssa.py -- lee el log de `barrido_ssa.py` y dibuja. No calcula nada.

La separacion es deliberada: la toma de datos puede tardar minutos y no debe
repetirse cada vez que se quiera mirar el mismo resultado desde otro angulo.
Todo lo que se dibuja aqui salio ya del disco.

Figuras que produce, en <log>/figuras/:

  datos_<ventana>.png          precio, volumen y volumen neto acumulado crudos
  barrido_L_<obs>.png          ortogonalidad contra L, con la linea base barajada
  ortogonalidad_<obs>.png      barras por ventana: cual es limpia y cual no
  espectro_<obs>.png           autovalores (scree) con la banda p95 de MC-SSA
  autovectores_<obs>_<v>.png   los autovectores con su periodo -- para ver patrones
  wcorr_<obs>_<v>.png          matriz de w-correlacion
  pares_<obs>_<v>.png          diagramas de fase U_i contra U_{i+1}
  color_<obs>.png              densidad espectral en log-log con la pendiente beta

Uso:
  python graficar_ssa.py --log=ssa_log/captura_v31b
  python graficar_ssa.py --log=ssa_log/captura_v31b --ventana=v00 --obs=precio
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec


COLOR_REAL = "#1f77b4"
COLOR_NULO = "#d62728"
COLOR_CTRL = "#2ca02c"


def _es_control(nombre: str) -> bool:
    return nombre.startswith("ctrl")


def _color_de(nombre: str) -> str:
    return COLOR_CTRL if _es_control(nombre) else COLOR_REAL


def cargar(dir_log: str) -> dict:
    with open(os.path.join(dir_log, "resumen.json"), "r", encoding="utf-8") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------

def fig_datos(dir_log: str, dir_fig: str, v: dict) -> None:
    """Los datos crudos de la ventana: tiempo, precio, ticks y volumen."""
    ruta = os.path.join(dir_log, v["ventana"], "datos.npz")
    if not os.path.exists(ruta):
        return
    d = np.load(ruta)
    t = d["t"] - d["t"][0]
    n = np.arange(len(t))

    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(4, 1, height_ratios=[3, 2, 2, 1.4], hspace=0.28)

    ax = fig.add_subplot(gs[0])
    ax.plot(n, d["precio"], lw=0.7, color=_color_de(v["ventana"]))
    ax.set_ylabel("precio [USD/BTC]")
    ax.set_title("%s | %d ticks | %.1f s | nu = %.2f tx/s"
                 % (v["ventana"], v["n"], v["dur_seg"], v["nu_tx_s"]))
    ax.grid(alpha=0.3)

    ax2 = fig.add_subplot(gs[1], sharex=ax)
    ax2.plot(n, d["vol_neto_acum"], lw=0.8, color="#9467bd")
    ax2.set_ylabel("vol. neto acum. [BTC]")
    ax2.axhline(0, color="k", lw=0.5)
    ax2.grid(alpha=0.3)

    ax3 = fig.add_subplot(gs[2], sharex=ax)
    ax3.plot(n, d["cant"], lw=0.4, color="#8c564b")
    ax3.set_ylabel("cantidad [BTC]")
    ax3.set_yscale("log")
    ax3.grid(alpha=0.3)

    # El desfase entre los dos relojes: cuanto tiempo de pared cuesta cada tick.
    ax4 = fig.add_subplot(gs[3], sharex=ax)
    ax4.plot(n, t, lw=0.8, color="#7f7f7f")
    ax4.set_ylabel("t [s]")
    ax4.set_xlabel("tick (reloj de transacciones)")
    ax4.grid(alpha=0.3)

    fig.savefig(os.path.join(dir_fig, "datos_%s.png" % v["ventana"]), dpi=110,
                bbox_inches="tight")
    plt.close(fig)


def fig_barrido(dir_log: str, dir_fig: str, meta: dict, obs: str) -> None:
    """Ortogonalidad contra L. Es la figura que justifica la eleccion de L."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    hay = False
    for v in meta["ventanas"]:
        if obs not in v["observables"]:
            continue
        ruta = os.path.join(dir_log, v["ventana"], "barrido_%s.npz" % obs)
        if not os.path.exists(ruta):
            continue
        b = np.load(ruta)
        hay = True
        est = "--" if _es_control(v["ventana"]) else "-"
        lc = _color_de(v["ventana"])
        ax.plot(b["Ls"], b["media_abs"], est, marker="o", ms=3, lw=1.1,
                color=lc, alpha=0.75, label=v["ventana"])
        ax.plot(b["Ls"], b["media_abs_sustituto"], ":", lw=0.9,
                color=COLOR_NULO, alpha=0.35)
        Lel = v["observables"][obs]["L_elegida"]
        i = int(np.where(b["Ls"] == Lel)[0][0])
        marca = "*" if v["observables"][obs]["hay_minimo_local"] else "x"
        ax.plot([Lel], [b["media_abs"][i]], marca, ms=13, color="k", zorder=5)
        # razon real/nulo: por debajo de 1 hay separacion que el nulo no produce
        ax2.plot(b["Ls"], b["media_abs"] / np.maximum(b["media_abs_sustituto"], 1e-12),
                 est, marker="o", ms=3, lw=1.1, color=lc, alpha=0.75, label=v["ventana"])

    if not hay:
        plt.close(fig)
        return

    ax.set_xscale("log")
    ax.set_xlabel("L [ticks]")
    ax.set_ylabel("media |w-correlacion| fuera de la diagonal")
    ax.set_title("Ortogonalidad contra L -- %s\n(* = minimo local, x = sin minimo local; "
                 "punteado rojo = sustituto barajado)" % obs, fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, ncol=2)

    ax2.axhline(1.0, color="k", lw=1.0)
    ax2.set_xscale("log")
    ax2.set_xlabel("L [ticks]")
    ax2.set_ylabel("razon real / sustituto barajado")
    ax2.set_title("Por debajo de 1 = separacion que el nulo NO produce", fontsize=10)
    ax2.grid(alpha=0.3, which="both")

    fig.savefig(os.path.join(dir_fig, "barrido_L_%s.png" % obs), dpi=110,
                bbox_inches="tight")
    plt.close(fig)


def fig_ortogonalidad(dir_fig: str, meta: dict, obs: str) -> None:
    nombres, med, nulo = [], [], []
    for v in meta["ventanas"]:
        if obs not in v["observables"]:
            continue
        o = v["observables"][obs]
        nombres.append(v["ventana"])
        med.append(o["ortogonalidad_media"])
        nulo.append(o["ortogonalidad_sustituto"])
    if not nombres:
        return

    x = np.arange(len(nombres))
    fig, ax = plt.subplots(figsize=(max(7, 0.85 * len(nombres) + 3), 4.6))
    ax.bar(x - 0.2, med, 0.4, color=[_color_de(n) for n in nombres], label="real")
    ax.bar(x + 0.2, nulo, 0.4, color=COLOR_NULO, alpha=0.55, label="barajado")
    ax.set_xticks(x)
    ax.set_xticklabels(nombres, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("media |w-correlacion|")
    ax.set_title("Limpieza de cada ventana -- %s (menor = mejor separacion)" % obs,
                 fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    fig.savefig(os.path.join(dir_fig, "ortogonalidad_%s.png" % obs), dpi=110,
                bbox_inches="tight")
    plt.close(fig)


def fig_espectro(dir_log: str, dir_fig: str, meta: dict, obs: str) -> None:
    """Scree de autovalores con la banda del nulo AR(1). Un autovalor por encima
    de p95 es una componente que el ruido rojo no explica."""
    ventanas = [v for v in meta["ventanas"] if obs in v["observables"]]
    if not ventanas:
        return
    ncol = min(3, len(ventanas))
    nfil = int(np.ceil(len(ventanas) / ncol))
    fig, axs = plt.subplots(nfil, ncol, figsize=(4.6 * ncol, 3.4 * nfil), squeeze=False)

    for k, v in enumerate(ventanas):
        ax = axs[k // ncol][k % ncol]
        ruta = os.path.join(dir_log, v["ventana"], "elegido_%s.npz" % obs)
        if not os.path.exists(ruta):
            ax.axis("off")
            continue
        e = np.load(ruta)
        lam = e["mc_lambda_obs"]
        i = np.arange(1, len(lam) + 1)
        o = v["observables"][obs]
        principal = o.get("nulo_principal", "ar1")
        estilos = {"ar1_incrementos": "-", "ar1": "--", "barajado": ":"}
        for nulo in o.get("mc_por_nulo", {"ar1": None}):
            clave = "mc_%s_p95" % nulo
            if clave not in e.files:
                continue
            if nulo == principal:
                ax.fill_between(i, e["mc_%s_p5" % nulo], e[clave],
                                color=COLOR_NULO, alpha=0.20)
            ax.plot(i, e[clave], estilos.get(nulo, "-."), lw=1.0, color=COLOR_NULO,
                    alpha=0.85, label="p95 %s" % nulo)
        ax.plot(i, lam, "o-", ms=3, lw=1.0, color=_color_de(v["ventana"]),
                label="observado")
        sob = e["mc_sobre"].astype(bool)
        if np.any(sob):
            ax.plot(i[sob], lam[sob], "o", ms=7, mfc="none", mec="k", mew=1.3)
        ax.set_yscale("log")
        cuentas = " ".join("%s=%d" % (k[:6], w["n_sobre"])
                           for k, w in o.get("mc_por_nulo", {}).items())
        ax.set_title("%s | L=%d | sobre p95: %s"
                     % (v["ventana"], o["L_elegida"], cuentas), fontsize=8)
        ax.set_xlabel("componente")
        ax.set_ylabel("autovalor")
        ax.grid(alpha=0.3, which="both")
        if k == 0:
            ax.legend(fontsize=7)

    for k in range(len(ventanas), nfil * ncol):
        axs[k // ncol][k % ncol].axis("off")
    fig.suptitle("Monte Carlo SSA -- %s   (banda sombreada = nulo principal; "
                 "ar1_incrementos reproduce el rebote bid-ask y ar1 sobre el nivel NO)\n"
                 "el test es ANTICONSERVADOR: proyecta sobre las EOF de los propios datos"
                 % obs, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(dir_fig, "espectro_%s.png" % obs), dpi=110,
                bbox_inches="tight")
    plt.close(fig)


def fig_autovectores(dir_log: str, dir_fig: str, v: dict, obs: str, n_mostrar: int = 12) -> None:
    """Los autovectores con su periodo dominante. Es la figura para buscar
    patrones: una EOF de un ciclo genuino es una sinusoide limpia; una de ruido
    parece una funcion de Legendre o una onda que se deshace."""
    ruta = os.path.join(dir_log, v["ventana"], "elegido_%s.npz" % obs)
    if not os.path.exists(ruta):
        return
    e = np.load(ruta)
    U = e["U"]
    lam = e["lambdas"]
    L = U.shape[0]
    sob = e["mc_sobre"].astype(bool)
    total = float(np.sum(lam))
    d = min(n_mostrar, U.shape[1])

    ncol = 4
    nfil = int(np.ceil(d / ncol))
    fig, axs = plt.subplots(nfil, ncol, figsize=(3.2 * ncol, 2.0 * nfil), squeeze=False)
    for i in range(d):
        ax = axs[i // ncol][i % ncol]
        ax.plot(U[:, i], lw=0.9, color="k" if not sob[i] else COLOR_CTRL)
        # Periodo dominante del propio autovector.
        esp = np.abs(np.fft.rfft(U[:, i] - np.mean(U[:, i]), n=8 * L)) ** 2
        f = np.fft.rfftfreq(8 * L, d=1.0)
        kk = int(np.argmax(esp[1:]) + 1)
        T = 1.0 / f[kk] if f[kk] > 0 else np.inf
        marca = " *" if sob[i] else ""
        ax.set_title("EOF %d%s | %.1f%% | T=%.0f tk"
                     % (i + 1, marca, 100 * lam[i] / total, T), fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(alpha=0.2)
    for i in range(d, nfil * ncol):
        axs[i // ncol][i % ncol].axis("off")
    fig.suptitle("Autovectores (EOF) -- %s / %s / L=%d   "
                 "(verde y * = supera el p95 del nulo AR(1))"
                 % (v["ventana"], obs, int(e["L"][0])), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(dir_fig, "autovectores_%s_%s.png" % (obs, v["ventana"])),
                dpi=110, bbox_inches="tight")
    plt.close(fig)


def fig_wcorr(dir_log: str, dir_fig: str, v: dict, obs: str) -> None:
    ruta = os.path.join(dir_log, v["ventana"], "elegido_%s.npz" % obs)
    if not os.path.exists(ruta):
        return
    e = np.load(ruta)
    W = np.abs(e["wcorr"])
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    im = ax.imshow(W, cmap="gray_r", vmin=0, vmax=1, origin="upper")
    ax.set_title("|w-correlacion| -- %s / %s / L=%d\nmedia fuera diag = %.4f"
                 % (v["ventana"], obs, int(e["L"][0]),
                    v["observables"][obs]["ortogonalidad_media"]), fontsize=10)
    ax.set_xlabel("componente")
    ax.set_ylabel("componente")
    fig.colorbar(im, ax=ax, shrink=0.85)
    fig.savefig(os.path.join(dir_fig, "wcorr_%s_%s.png" % (obs, v["ventana"])),
                dpi=110, bbox_inches="tight")
    plt.close(fig)


def fig_pares(dir_log: str, dir_fig: str, v: dict, obs: str, n_pares: int = 6) -> None:
    """Diagrama de fase U_i contra U_{i+1}. Un par oscilatorio genuino dibuja un
    circulo (cuadratura); si sale una recta o una nube, no hay ciclo."""
    ruta = os.path.join(dir_log, v["ventana"], "elegido_%s.npz" % obs)
    if not os.path.exists(ruta):
        return
    e = np.load(ruta)
    U = e["U"]
    d = min(2 * n_pares, U.shape[1])
    detect = {tuple(p["indices"]): p for p in v["observables"][obs].get("pares", [])}

    ncol = 3
    n = d // 2
    nfil = int(np.ceil(n / ncol))
    fig, axs = plt.subplots(nfil, ncol, figsize=(3.0 * ncol, 3.0 * nfil), squeeze=False)
    for k in range(n):
        i = 2 * k
        ax = axs[k // ncol][k % ncol]
        es_par = (i, i + 1) in detect
        ax.plot(U[:, i], U[:, i + 1], lw=0.7,
                color=COLOR_CTRL if es_par else "0.4")
        ax.set_aspect("equal", "datalim")
        ax.set_xticks([])
        ax.set_yticks([])
        if es_par:
            p = detect[(i, i + 1)]
            ax.set_title("EOF %d-%d  PAR\nT=%.0f tk = %.1f s"
                         % (i + 1, i + 2, p["periodo_muestras"], p["periodo_seg"]),
                         fontsize=8)
        else:
            ax.set_title("EOF %d-%d" % (i + 1, i + 2), fontsize=8)
        ax.grid(alpha=0.2)
    for k in range(n, nfil * ncol):
        axs[k // ncol][k % ncol].axis("off")
    fig.suptitle("Diagramas de fase -- %s / %s (circulo = cuadratura = ciclo)"
                 % (v["ventana"], obs), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(dir_fig, "pares_%s_%s.png" % (obs, v["ventana"])),
                dpi=110, bbox_inches="tight")
    plt.close(fig)


def fig_color(dir_log: str, dir_fig: str, meta: dict, obs: str) -> None:
    """Densidad espectral en log-log de nivel, incrementos y residuo del SSA."""
    ventanas = [v for v in meta["ventanas"] if obs in v["observables"]]
    if not ventanas:
        return
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
    titulos = [("col_nivel", "nivel", "color_nivel"),
               ("col_incr", "incrementos", "color_incrementos"),
               ("col_res", "residuo del SSA", "color_residuo_ssa")]

    for j, (pref, titulo, clave) in enumerate(titulos):
        ax = axs[j]
        for v in ventanas:
            ruta = os.path.join(dir_log, v["ventana"], "elegido_%s.npz" % obs)
            if not os.path.exists(ruta):
                continue
            e = np.load(ruta)
            f, S = e["%s_f" % pref], e["%s_S" % pref]
            m = (f > 0) & (S > 0)
            if np.count_nonzero(m) < 4:
                continue
            c = v["observables"][obs][clave]
            est = "--" if _es_control(v["ventana"]) else "-"
            ax.loglog(f[m], S[m], est, lw=0.9, alpha=0.75, color=_color_de(v["ventana"]),
                      label="%s b=%+.2f %s" % (v["ventana"], c.get("beta", np.nan),
                                               "" if c.get("valido") else "(no ley pot.)"))
        ax.set_title("Densidad espectral -- %s" % titulo, fontsize=10)
        ax.set_xlabel("f [Hz]")
        ax.set_ylabel("S(f)")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=6, ncol=1)
    fig.suptitle("Color de ruido -- %s  (S(f) ~ f^-beta: 0 blanco, 1 rosa, 2 rojo)" % obs,
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(dir_fig, "color_%s.png" % obs), dpi=110, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", required=True)
    ap.add_argument("--obs", default=None, help="uno solo; por omision todos")
    ap.add_argument("--ventana", default=None, help="una sola; por omision todas")
    args = ap.parse_args()

    dir_log = args.log
    meta = cargar(dir_log)
    dir_fig = os.path.join(dir_log, "figuras")
    os.makedirs(dir_fig, exist_ok=True)

    observables = [args.obs] if args.obs else meta["observables"]
    ventanas = [v for v in meta["ventanas"]
                if args.ventana is None or v["ventana"] == args.ventana]

    print("log        : %s" % dir_log)
    print("captura    : %s (%d tx, nu = %.2f tx/s)"
          % (meta["captura"], meta["n_transacciones"], meta["nu_tx_s"] or float("nan")))
    print("ventanas   : %d | observables: %s" % (len(ventanas), ", ".join(observables)))
    print("")

    for v in ventanas:
        fig_datos(dir_log, dir_fig, v)
    print("datos crudos          -> datos_*.png")

    for obs in observables:
        fig_barrido(dir_log, dir_fig, meta, obs)
        fig_ortogonalidad(dir_fig, meta, obs)
        fig_espectro(dir_log, dir_fig, meta, obs)
        fig_color(dir_log, dir_fig, meta, obs)
        for v in ventanas:
            if obs not in v["observables"]:
                continue
            fig_autovectores(dir_log, dir_fig, v, obs)
            fig_wcorr(dir_log, dir_fig, v, obs)
            fig_pares(dir_log, dir_fig, v, obs)
        print("observable %-14s -> barrido_L, ortogonalidad, espectro, color, "
              "autovectores, wcorr, pares" % obs)

    n = len([f for f in os.listdir(dir_fig) if f.endswith(".png")])
    print("")
    print("%d figuras en %s" % (n, dir_fig))


if __name__ == "__main__":
    main()
