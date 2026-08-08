"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: experimento_v30.py — §3 de ORDEN_TRABAJO_OSCILADOR_3_0

    python experimento_v30.py

EL TEST QUE DECIDE. Mide `k`, `m`, `γ` y `Q` por AR(2) sobre todos los datos, sin
ventana espectral y sin banda de resolubilidad, y contrasta `H₀: k = 0` contra un
nulo que es **exactamente un paseo aleatorio con los incrementos observados** —la
hipótesis que la v2.2 no pudo rechazar.

⚠ NO TOCA `Micelio.py`. §8 de la v3.0.
"""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy import stats

import oscilador as OSC


def cargar_capturas():
    """Series limpias, con el filtro de validez del feed aplicado.

    ⚠ EL FILTRADO NO ES COSMETICO. El feed emite ~0.1-0.24 % de mensajes con
    p = 0 y q = 0 (ver `mercado.tick_valido`). Sin quitarlos, el |incremento|
    máximo pasa de 11.8 USD a 65 245 USD y cualquier estimación de varianza
    —incluido el AR(2)— queda dominada por ellos.
    """
    series = {}
    try:
        import captura_larga as CL

        d = CL.cargar_larga()
        p, t = d["tr_precio"], d["tr_t"]
        ok = p > 0
        series["captura_larga"] = (p[ok], t[ok], int((~ok).sum()))
    except Exception as err:
        print(f"[AVISO] captura_larga no disponible: {err}")
    ruta = os.path.join("telemetria", "captura_dual.npz")
    if os.path.exists(ruta):
        d = np.load(ruta)
        p, t = d["tr_precio"], d["tr_t"]
        ok = p > 0
        series["captura_dual"] = (p[ok], t[ok], int((~ok).sum()))
    return series


def tramos_continuos(t, hueco_max=60.0):
    """Índices de tramos sin huecos: un corte de red parte la serie en dos.

    Concatenar a través de un hueco de 10 h fabricaría un salto de precio que el
    AR(2) leería como un incremento gigante. La captura tuvo un corte de DNS de
    36 442 s, así que esto no es hipotético.
    """
    corte = np.flatnonzero(np.diff(t) > hueco_max) + 1
    ini = np.concatenate([[0], corte])
    fin = np.concatenate([corte, [len(t)]])
    return [(a, b) for a, b in zip(ini, fin) if b - a > 1000]


def informe(res, nulo_ref=None):
    m = res["mco"]
    h = res["huber"]
    print(f"  [{res['etiqueta']}]  n = {res['n']}")
    print(f"     phi1 = {res['phi1']:+.6f} (+-{res['se_phi1']:.6f})   "
          f"phi2 = {res['phi2']:+.6f} (+-{res['se_phi2']:.6f})")
    print(f"     MCO   : k = {m['k']:+.3e} (+-{res['se_k']:.1e})  "
          f"gamma = {m['gamma']:+.6f}  m = {m['m']:+.6f}  "
          f"Q = {m['Q']:.4f}  raices {'complejas' if m['raices_complejas'] else 'REALES'}")
    print(f"     Huber : k = {h['k']:+.3e}  gamma = {h['gamma']:+.6f}  "
          f"Q = {h['Q']:.4f}  raices {'complejas' if h['raices_complejas'] else 'REALES'}")
    print(f"     Q por bootstrap: [{res['Q_boot_p05']:.4f}, {res['Q_boot_p95']:.4f}] al 90 %")
    print(f"     H0: k = 0  ->  p = {res['p_valor_k']:.4f}   "
          f"(nulo: k p50 = {res['k_nulo_p50']:.3e}, p95 = {res['k_nulo_p95']:.3e})")
    aviso = OSC.guarda_gamma(res)
    if aviso:
        print(f"     [!] {aviso}")
    # Lectura
    if res["p_valor_k"] > 0.05:
        print("     >> NO se rechaza H0: k = 0. Sin fuerza recuperadora, SIN OSCILADOR.")
    elif m["Q"] > 0.5 and m["raices_complejas"]:
        print(f"     >> k > 0 con Q = {m['Q']:.3f} > 1/2 y raices complejas: OSCILADOR "
              f"SUBAMORTIGUADO.")
    else:
        print(f"     >> k > 0 pero Q = {m['Q']:.3f} <= 1/2: reversion SIN oscilacion "
              f"(sobreamortiguado).")
    print()


def por_bloques(x, t, rng, n_bloques_max=40, minutos=5.0):
    """§3.3: estabilidad temporal. ¿`k` y `Q` varían tanto como variaba `ω_m`?"""
    print("  ESTABILIDAD POR BLOQUES (§3.3)")
    bordes = np.arange(t[0], t[-1], minutos * 60.0)
    ks, Qs, gs, n_ok = [], [], [], 0
    for i in range(len(bordes) - 1):
        sel = (t >= bordes[i]) & (t < bordes[i + 1])
        if sel.sum() < 2000:
            continue
        p1, p2, *_ = OSC.ajustar_ar2(x[sel])
        pr = OSC.primitivas_desde_phi(p1, p2)
        ks.append(pr["k"])
        Qs.append(pr["Q"])
        gs.append(pr["gamma"])
        n_ok += 1
        if n_ok >= n_bloques_max:
            break
    if n_ok < 3:
        print(f"     solo {n_ok} bloques utiles: insuficiente")
        return
    ks, Qs, gs = np.array(ks), np.array(Qs), np.array(gs)
    Qf = Qs[np.isfinite(Qs)]
    print(f"     {n_ok} bloques de {minutos:.0f} min")
    print(f"     k     : p10={np.percentile(ks,10):+.3e} MED={np.median(ks):+.3e} "
          f"p90={np.percentile(ks,90):+.3e}   dispersion "
          f"{np.std(ks)/max(1e-30,abs(np.mean(ks))):.1%}")
    if Qf.size:
        print(f"     Q     : p10={np.percentile(Qf,10):.4f} MED={np.median(Qf):.4f} "
              f"p90={np.percentile(Qf,90):.4f}  (finitos {Qf.size}/{n_ok})")
    print(f"     gamma : negativa en {int((gs<0).sum())}/{n_ok} bloques")
    print(f"     bloques con raices complejas: "
          f"{int(np.sum(np.isfinite(Qs)))}/{n_ok}")
    print()


def main():
    print("=" * 92)
    print("§2.3 — VERIFICACION DIMENSIONAL")
    print("=" * 92)
    d = OSC.verificar_dimensiones()
    print(f"  los cuatro terminos comparten unidades (USD,BTC,Tick) = {d['termino_comun']}")
    print(f"  Q adimensional: {d['Q_adimensional']}   "
          f"lambda*S^2*Gamma adimensional: {d['lambda_S2_Gamma_adimensional']}")
    print(f"  [lambda] = {OSC.UNIDADES['lambda']}  (Sec. 4.4.1 del PDF, EXTRAIDA no asumida)")
    print("  [!] m  prop. a  1/lambda NO es una igualdad: [1/lambda]=BTC contra [m]=[k]*Tick^2.")
    print("    La proporcion exige un factor de conversion; se anota, no se cuela.")
    print()

    series = cargar_capturas()
    if not series:
        print("[ERROR] sin capturas")
        return 1
    rng = np.random.default_rng(3030)

    for nombre, (p, t, n_ceros) in series.items():
        dur = t[-1] - t[0]
        print("=" * 92)
        print(f"§3 — AR(2) SOBRE {nombre}")
        print("=" * 92)
        print(f"  {len(p)} trades validos en {dur/3600:.2f} h "
              f"(nu = {len(p)/dur:.1f} tx/s); {n_ceros} descartados por p=0 o q=0")
        tr = tramos_continuos(t)
        print(f"  tramos continuos (sin huecos > 60 s): {len(tr)} "
              f"-> {[b-a for a,b in tr]}")
        a, b = max(tr, key=lambda z: z[1] - z[0])
        x, tt = p[a:b], t[a:b]
        print(f"  se analiza el tramo mas largo: {len(x)} ticks, "
              f"{(tt[-1]-tt[0])/3600:.2f} h")
        print()
        res = OSC.analizar_serie(x, rng, etiqueta=f"{nombre}: precio de transaccion")
        informe(res)
        por_bloques(x, tt, rng)
    return 0


if __name__ == "__main__":
    sys.exit(main())
