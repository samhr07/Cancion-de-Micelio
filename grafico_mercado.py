# -*- coding: utf-8 -*-
"""
grafico_mercado.py -- grafico de velas al estilo de una plataforma de trading,
sobre EXACTAMENTE los mismos datos con los que se hizo el estudio SSA.

No calcula nada nuevo ni decide nada: lee las capturas de `telemetria/` con el
mismo cargador que `barrido_ssa.py` (mismo filtro de los ceros del feed, mismo
orden temporal) y las dibuja. Sirve para mirar con los ojos el tramo de mercado
sobre el que se midieron la escalera de ventana y la relacion phi'-volatilidad.

CUATRO PANELES
--------------
  1. Velas OHLC + medias moviles.
  2. Volumen por vela, coloreado segun la vela suba o baje.
  3. Volumen neto firmado (delta) y su acumulado. Necesita `tr_maker`, que solo
     persisten las capturas posteriores a la v3.1.
  4. phi' = ticks por BTC inyectado, que es la variable que si resulto tener
     relacion medible con el mercado (sesion 2026-08-08 f).

LOS DOS RELOJES
---------------
`--reloj=tiempo` agrega por segundos de reloj de pared, que es lo que hace
cualquier broker. `--reloj=ticks` agrega cada N transacciones, que es el reloj en
el que trabaja el sistema desde la v2.0. Las dos vistas del mismo tramo se ven
distintas, y esa diferencia es justamente el asunto de media docena de sesiones:
con nu variando por un factor 20, una vela de 60 s puede llevar 300 o 6000
transacciones.

Uso:
  python grafico_mercado.py                                  # captura_v31b, velas de 60 s
  python grafico_mercado.py --intervalo=300                  # velas de 5 min
  python grafico_mercado.py --reloj=ticks --intervalo=256    # velas de 256 ticks
  python grafico_mercado.py --captura=telemetria/captura_larga
  python grafico_mercado.py --desde=0 --hasta=8192           # una ventana del estudio
  python grafico_mercado.py --guardar=grafico.png            # sin ventana interactiva
"""

from __future__ import annotations

import argparse
import datetime as dt
import os

import numpy as np

import barrido_ssa as bs


VERDE = "#26a69a"
ROJO = "#ef5350"
FONDO = "#131722"
REJILLA = "#2a2e39"
TEXTO = "#d1d4dc"


def construir_velas(t, p, q, eps, reloj: str, intervalo: float) -> dict:
    """Agrega las transacciones en velas OHLCV.

    Con reloj de tiempo se agrupa por intervalos de pared; con reloj de ticks,
    cada `intervalo` transacciones. En ambos casos las velas se dibujan
    equiespaciadas, como en cualquier plataforma: los huecos de mercado no se
    reservan espacio, se anotan en el eje.
    """
    if reloj == "tiempo":
        clave = np.floor((t - t[0]) / intervalo).astype(np.int64)
    elif reloj == "ticks":
        clave = (np.arange(len(t)) // int(intervalo)).astype(np.int64)
    else:
        raise SystemExit("reloj desconocido: %s" % reloj)

    # Fronteras de grupo. `clave` es monotona no decreciente en los dos casos.
    cortes = np.flatnonzero(np.diff(clave)) + 1
    ini = np.concatenate(([0], cortes))
    fin = np.concatenate((cortes, [len(t)]))

    n = len(ini)
    O = np.empty(n); H = np.empty(n); L = np.empty(n); C = np.empty(n)
    V = np.empty(n); NETO = np.empty(n); T0 = np.empty(n); NT = np.empty(n, dtype=np.int64)
    for k in range(n):
        a, b = ini[k], fin[k]
        tramo = p[a:b]
        O[k] = tramo[0]
        H[k] = tramo.max()
        L[k] = tramo.min()
        C[k] = tramo[-1]
        V[k] = q[a:b].sum()
        NETO[k] = np.nansum(eps[a:b] * q[a:b])
        T0[k] = t[a]
        NT[k] = b - a

    with np.errstate(divide="ignore", invalid="ignore"):
        phi = np.where(V > 0, NT / V, np.nan)   # ticks por BTC

    return {"O": O, "H": H, "L": L, "C": C, "V": V, "neto": NETO,
            "t": T0, "n_ticks": NT, "phi": phi, "n": n}


def media_movil(x: np.ndarray, w: int) -> np.ndarray:
    if w <= 1 or w > x.size:
        return np.full_like(x, np.nan)
    nucleo = np.ones(w) / w
    m = np.convolve(x, nucleo, mode="valid")
    return np.concatenate((np.full(w - 1, np.nan), m))


def dibujar(v: dict, titulo: str, tiene_maker: bool, ruta: str | None,
            reloj: str, intervalo: float) -> None:
    import matplotlib
    if ruta:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib import gridspec

    n = v["n"]
    x = np.arange(n)
    sube = v["C"] >= v["O"]
    colores = np.where(sube, VERDE, ROJO)

    filas = 4 if tiene_maker else 3
    alturas = [3.2, 1.0, 1.0, 1.0][:filas]
    fig = plt.figure(figsize=(15, 9), facecolor=FONDO)
    gs = gridspec.GridSpec(filas, 1, height_ratios=alturas, hspace=0.08)

    def estilo(ax):
        ax.set_facecolor(FONDO)
        ax.grid(color=REJILLA, lw=0.6, alpha=0.8)
        ax.tick_params(colors=TEXTO, labelsize=8)
        for lado in ax.spines.values():
            lado.set_color(REJILLA)
        ax.yaxis.label.set_color(TEXTO)

    # --- panel 1: velas ---
    ax = fig.add_subplot(gs[0])
    estilo(ax)
    # Mecha: una linea vertical de minimo a maximo.
    ax.vlines(x, v["L"], v["H"], color=colores, lw=0.8, zorder=2)
    # Cuerpo: rectangulo de apertura a cierre. Las velas planas (O == C) se
    # dibujan como una raya, que es lo que hace cualquier plataforma.
    ancho = 0.68
    for i in range(n):
        alto = abs(v["C"][i] - v["O"][i])
        base = min(v["O"][i], v["C"][i])
        if alto == 0:
            ax.hlines(base, i - ancho / 2, i + ancho / 2, color=colores[i], lw=1.2, zorder=3)
        else:
            ax.add_patch(Rectangle((i - ancho / 2, base), ancho, alto,
                                   facecolor=colores[i], edgecolor=colores[i],
                                   lw=0.5, zorder=3))
    for w, col in ((7, "#f0b90b"), (25, "#2962ff")):
        if n > w:
            ax.plot(x, media_movil(v["C"], w), color=col, lw=1.1,
                    label="MM%d" % w, zorder=4)
    ax.set_ylabel("precio [USD/BTC]")
    ax.set_title(titulo, color=TEXTO, fontsize=11)
    ax.legend(fontsize=8, facecolor=FONDO, edgecolor=REJILLA, labelcolor=TEXTO)
    ax.set_xlim(-1, n)

    # --- panel 2: volumen ---
    ax2 = fig.add_subplot(gs[1], sharex=ax)
    estilo(ax2)
    ax2.bar(x, v["V"], width=ancho, color=colores, alpha=0.85)
    ax2.set_ylabel("volumen [BTC]")

    # --- panel 3: volumen neto (delta) ---
    idx = 2
    if tiene_maker:
        ax3 = fig.add_subplot(gs[idx], sharex=ax)
        estilo(ax3)
        ax3.bar(x, v["neto"], width=ancho,
                color=np.where(v["neto"] >= 0, VERDE, ROJO), alpha=0.85)
        ax3.axhline(0, color=TEXTO, lw=0.6)
        ax3.set_ylabel("vol. neto [BTC]")
        axc = ax3.twinx()
        axc.plot(x, np.cumsum(np.nan_to_num(v["neto"])), color="#f0b90b", lw=1.0)
        axc.tick_params(colors=TEXTO, labelsize=7)
        axc.set_ylabel("acumulado", color=TEXTO, fontsize=8)
        for lado in axc.spines.values():
            lado.set_color(REJILLA)
        idx += 1

    # --- panel final: phi' ---
    ax4 = fig.add_subplot(gs[idx], sharex=ax)
    estilo(ax4)
    ax4.plot(x, v["phi"], color="#b39ddb", lw=0.9)
    ax4.set_yscale("log")
    ax4.set_ylabel("phi' [ticks/BTC]")

    # Eje temporal: etiquetas con la hora real, pero velas equiespaciadas.
    paso = max(1, n // 12)
    pos = x[::paso]
    etiq = [dt.datetime.fromtimestamp(v["t"][i]).strftime("%H:%M:%S") for i in pos]
    ax4.set_xticks(pos)
    ax4.set_xticklabels(etiq, rotation=30, ha="right")
    ax4.set_xlabel("hora local  (reloj de %s, %s por vela)"
                   % (reloj, ("%.0f s" % intervalo) if reloj == "tiempo"
                      else ("%d ticks" % intervalo)), color=TEXTO)
    for a in fig.axes[:-1]:
        a.tick_params(labelbottom=False)

    if ruta:
        fig.savefig(ruta, dpi=110, bbox_inches="tight", facecolor=FONDO)
        print("guardado en %s" % ruta)
    else:
        plt.show()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--captura", default="telemetria/captura_v31b")
    ap.add_argument("--reloj", default="tiempo", choices=["tiempo", "ticks"])
    ap.add_argument("--intervalo", type=float, default=60.0,
                    help="segundos por vela (reloj=tiempo) o ticks por vela (reloj=ticks)")
    ap.add_argument("--desde", type=int, default=0, help="primer tick a dibujar")
    ap.add_argument("--hasta", type=int, default=0, help="ultimo tick; 0 = hasta el final")
    ap.add_argument("--guardar", default=None, help="ruta PNG; sin esto abre ventana")
    args = ap.parse_args()

    d = bs.cargar_captura(args.captura)
    t = d["tr_t"]; p = d["tr_precio"]
    q = d["tr_cant"] if "tr_cant" in d else np.zeros_like(p)
    tiene_maker = "tr_maker" in d
    if tiene_maker:
        # Convencion del proyecto: m = True -> el comprador es maker -> el
        # agresor es el vendedor -> eps = -1.
        eps = np.where(d["tr_maker"].astype(bool), -1.0, 1.0)
    else:
        eps = np.zeros_like(p)

    a = max(0, args.desde)
    b = len(t) if args.hasta <= 0 else min(len(t), args.hasta)
    if b - a < 10:
        raise SystemExit("el tramo pedido tiene menos de 10 transacciones")
    t, p, q, eps = t[a:b], p[a:b], q[a:b], eps[a:b]

    v = construir_velas(t, p, q, eps, args.reloj, args.intervalo)
    dur = float(t[-1] - t[0])
    nu = len(t) / dur if dur > 0 else float("nan")

    print("captura      : %s" % args.captura)
    print("tramo        : ticks %d a %d (%d transacciones)" % (a, b, len(t)))
    print("duracion     : %.2f h  ->  nu = %.2f tx/s" % (dur / 3600.0, nu))
    print("velas        : %d (%s)"
          % (v["n"], ("%.0f s" % args.intervalo) if args.reloj == "tiempo"
             else ("%d ticks" % args.intervalo)))
    print("ticks por vela: p10=%d  mediana=%d  p90=%d  max=%d"
          % (np.percentile(v["n_ticks"], 10), np.median(v["n_ticks"]),
             np.percentile(v["n_ticks"], 90), v["n_ticks"].max()))
    print("precio       : %.2f -> %.2f  (min %.2f, max %.2f)"
          % (v["O"][0], v["C"][-1], v["L"].min(), v["H"].max()))
    print("volumen      : %.4f BTC total | neto %+.4f BTC%s"
          % (v["V"].sum(), np.nansum(v["neto"]),
             "" if tiene_maker else "  (sin tr_maker: el neto no es fiable)"))
    print("phi'         : mediana %.1f ticks/BTC (p10 %.1f, p90 %.1f)"
          % (np.nanmedian(v["phi"]), np.nanpercentile(v["phi"], 10),
             np.nanpercentile(v["phi"], 90)))
    if not tiene_maker:
        print("AVISO: esta captura no persiste `tr_maker`; se omite el panel de volumen neto.")

    titulo = ("%s  |  %d tx en %.2f h  |  nu = %.2f tx/s  |  %d velas de %s"
              % (os.path.basename(os.path.normpath(args.captura)), len(t),
                 dur / 3600.0, nu, v["n"],
                 ("%.0f s" % args.intervalo) if args.reloj == "tiempo"
                 else ("%d ticks" % args.intervalo)))
    dibujar(v, titulo, tiene_maker, args.guardar, args.reloj, args.intervalo)


if __name__ == "__main__":
    main()
