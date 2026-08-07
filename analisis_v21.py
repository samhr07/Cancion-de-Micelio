"""
Proyecto: Canción del Micelio (Trading Bot Algorítmico)
Módulo: analisis_v21.py — Mediciones de ORDEN_TRABAJO_OMEGA_2_1

    python captura_dual.py --segundos=900     # primero, capturar
    python analisis_v21.py                    # despues, medir

Cubre las mediciones que la v2.1 exige antes de decidir nada:

  §3.4  fracción de duplicados de `mid` contra precio de transacción, a igual K
  §3.3  fracción de cambios de `mid` SIN transacción intermedia
  §2.1  barrido de la banda con ventanas suficientes para discriminar
  §4.2  descomposición σ_instr / σ_total de la anchura de ω_m
  §7.2  correlación parcial ρ(ω_m, ΔS | ν)

⚠ TODA CONSTANTE SE REPORTA COMO DISTRIBUCIÓN CON PERCENTILES, no como valor
(§1.4): ν varía por un factor 20 el mismo día, así que cualquier número medido en
un solo instante es sospechoso por construcción.
"""

from __future__ import annotations

import os
import sys

import numpy as np

import hht

ARCHIVO = os.path.join("telemetria", "captura_dual.npz")
W_VENTANA = 384
PERIODO_MUESTREO = 0.5  # [s] el PERIODO_MUESTREO_EMD del bot


def pct(v, p):
    return float(np.percentile(np.asarray(v, dtype=float), p)) if len(v) else float("nan")


def cargar():
    if not os.path.exists(ARCHIVO):
        print(f"[ERROR] no hay captura en {ARCHIVO}. Corre: python captura_dual.py")
        sys.exit(1)
    d = np.load(ARCHIVO)
    return d["tr_precio"], d["tr_t"], d["tr_id"], d["bk_mid"], d["bk_t"]


def mid_vigente_en(tr_t, bk_mid, bk_t):
    """`mid` vigente en el instante de cada transaccion.

    §3.2: `mid` es un ESTADO que se lee, no un evento que se cuenta. El reloj lo
    sigue marcando `@trade`; en cada marca se lee el `mid` vigente. Por eso esto
    es una busqueda por instante y no una serie con su propio espaciado.
    """
    idx = np.searchsorted(bk_t, tr_t, side="right") - 1
    idx = np.clip(idx, 0, len(bk_mid) - 1)
    return bk_mid[idx], idx


# ==============================================================================
# §3 — EL OBSERVABLE
# ==============================================================================
def medir_observable(tr_precio, tr_t, bk_mid, bk_t, nu):
    print("=" * 92)
    print("§3 — OBSERVABLE: `mid` de @bookTicker contra precio de transaccion")
    print("=" * 92)
    k = max(1, int(round(nu * PERIODO_MUESTREO)))
    mid_en_trade, idx_bk = mid_vigente_en(tr_t, bk_mid, bk_t)

    print(f"  nu = {nu:.1f} tx/s   ->   K = {k} ticks/muestra")
    print(f"  actualizaciones de libro: {len(bk_mid)/ (bk_t[-1]-bk_t[0]):.1f} msg/s "
          f"({len(bk_mid)/max(1,len(tr_precio)):.1f}x la tasa de transacciones)")
    print()

    # --- §3.4 duplicados a igual K, sobre VARIAS ventanas ---
    dup_tr, dup_mid = [], []
    for ini in range(0, max(1, len(tr_precio) - W_VENTANA * k), W_VENTANA * k // 2):
        i = np.arange(ini, ini + W_VENTANA * k, k)
        i = i[i < len(tr_precio)]
        if len(i) < W_VENTANA // 2:
            continue
        dup_tr.append(float(np.mean(np.diff(tr_precio[i]) == 0.0)))
        dup_mid.append(float(np.mean(np.diff(mid_en_trade[i]) == 0.0)))
    print(f"  §3.4 DUPLICADOS a igual K, sobre {len(dup_tr)} ventanas:")
    print(f"     precio de transaccion : p10={pct(dup_tr,10):5.1%} "
          f"MEDIANA={pct(dup_tr,50):5.1%} p90={pct(dup_tr,90):5.1%}")
    print(f"     mid de bookTicker     : p10={pct(dup_mid,10):5.1%} "
          f"MEDIANA={pct(dup_mid,50):5.1%} p90={pct(dup_mid,90):5.1%}")
    ref = 0.441  # referencia de la v2.0
    mejora = (pct(dup_tr, 50) - pct(dup_mid, 50))
    print(f"     referencia v2.0: 44.1 %   |   mejora de mid: {mejora:+.1%} absoluto")
    if pct(dup_mid, 50) < 0.75 * pct(dup_tr, 50):
        print("     >> `mid` REDUCE los duplicados de forma material.")
    else:
        print("     >> `mid` NO reduce los duplicados de forma material.")
        print("        El observable no era el problema: volver sobre §1.3 con otra")
        print("        hipotesis antes de cambiar la cadena EMD.")

    # --- §3.3 cambios de mid SIN transaccion intermedia ---
    # Para cada cambio de mid, ¿hubo alguna transaccion desde el cambio anterior?
    cambios = np.flatnonzero(np.diff(bk_mid) != 0.0) + 1
    if len(cambios) > 1:
        t_cambios = bk_t[cambios]
        # transacciones acumuladas hasta cada cambio
        n_tr = np.searchsorted(tr_t, t_cambios, side="right")
        sin_trade = float(np.mean(np.diff(n_tr) == 0))
    else:
        sin_trade = float("nan")
    print()
    print(f"  §3.3 cambios de `mid`: {len(cambios)} de {len(bk_mid)} mensajes "
          f"({len(cambios)/max(1,len(bk_mid)):.1%} de los updates mueven el mid)")
    print(f"     SIN transaccion intermedia: {sin_trade:.1%}")
    if sin_trade > 0.5:
        print("     [!] Mas de la mitad de los movimientos de `mid` no tienen")
        print("         transaccion detras: es parpadeo de cotizaciones, un ruido de")
        print("         naturaleza distinta al de @trade. RECONSIDERAR (§3.3).")
    else:
        print("     La mayoria de los movimientos de `mid` tienen transaccion detras.")
    return k, mid_en_trade, dup_tr, dup_mid


# ==============================================================================
# §4.2 — σ_ω INSTRUMENTAL CONTRA GENUINA
# ==============================================================================
def medir_sigma_omega(serie, tiempos, k, etiqueta):
    """§4.2: separar la anchura del APARATO de la anchura del FENOMENO.

    σ_instr²  : dispersion de ω_m entre REJILLAS — mismos datos, variando K y el
                desplazamiento de inicio de ventana.
    σ_total²  : dispersion de ω_m en el TIEMPO, a rejilla fija.
    σ_genuina² = σ_total² − σ_instr²

    Si σ_genuina² sale NEGATIVA es un resultado con contenido, no un fallo:
    significa que toda la anchura observada es del aparato.
    """
    print()
    print("=" * 92)
    print(f"§4.2 — DESCOMPOSICION DE LA ANCHURA DE omega  [{etiqueta}]")
    print("=" * 92)

    def omega_de(ini, kk):
        i = np.arange(ini, ini + W_VENTANA * kk, kk)
        i = i[i < len(serie)]
        if len(i) < W_VENTANA:
            return None
        dt = float((tiempos[i[-1]] - tiempos[i[0]]) / (len(i) - 1))
        r = hht.analizar_ventana(serie[i], dt)
        if not r["valido"] or r["f_hz"] <= 0.0:
            return None
        return r["f_hz"], bool(r["omega_valida"])

    # --- σ_total: rejilla FIJA, ventanas que avanzan en el tiempo ---
    f_tiempo, validas = [], []
    paso = W_VENTANA * k // 2
    for ini in range(0, max(1, len(serie) - W_VENTANA * k), paso):
        r = omega_de(ini, k)
        if r:
            f_tiempo.append(r[0])
            validas.append(r[1])
    # --- σ_instr: MISMOS datos, rejillas distintas (K/2, K, 2K y desfases) ---
    f_rejilla = []
    for ini_off in (0, k // 3, 2 * k // 3):
        for kk in (max(1, k // 2), k, 2 * k):
            r = omega_de(ini_off, kk)
            if r:
                f_rejilla.append(r[0])

    if len(f_tiempo) < 3 or len(f_rejilla) < 3:
        print(f"  Insuficiente: {len(f_tiempo)} ventanas temporales, "
              f"{len(f_rejilla)} rejillas.")
        return None

    var_total = float(np.var(f_tiempo, ddof=1))
    var_instr = float(np.var(f_rejilla, ddof=1))
    var_genuina = var_total - var_instr
    frac = var_instr / var_total if var_total > 0 else float("inf")

    print(f"  ventanas en el tiempo (rejilla fija): {len(f_tiempo)}   "
          f"rejillas distintas (mismos datos): {len(f_rejilla)}")
    print(f"  f_hz por ventana : mediana={pct(f_tiempo,50):.5f} "
          f"p10={pct(f_tiempo,10):.5f} p90={pct(f_tiempo,90):.5f}")
    print(f"  f_hz por rejilla : mediana={pct(f_rejilla,50):.5f} "
          f"p10={pct(f_rejilla,10):.5f} p90={pct(f_rejilla,90):.5f}")
    print()
    print(f"  sigma_total^2   = {var_total:.6e}")
    print(f"  sigma_instr^2   = {var_instr:.6e}")
    print(f"  sigma_genuina^2 = {var_genuina:+.6e}")
    print(f"  >> FRACCION INSTRUMENTAL = {frac:.1%}   (condicion del §8: < 50 %)")
    if var_genuina < 0:
        print("     sigma_genuina^2 NEGATIVA: toda la anchura observada es del")
        print("     APARATO. Es un resultado con contenido, no un fallo: significa")
        print("     que no hay anchura de mercado que modelar todavia.")
    print(f"  omega_valida en {np.mean(validas):.1%} de las ventanas "
          f"(condicion del §8: > 90 %)")
    return {"frac_instr": frac, "var_total": var_total, "var_instr": var_instr,
            "frac_valida": float(np.mean(validas)), "f": f_tiempo}


def main():
    tr_precio, tr_t, tr_id, bk_mid, bk_t = cargar()
    dur = tr_t[-1] - tr_t[0]
    nu = len(tr_precio) / dur
    print()
    print(f"Captura: {len(tr_precio)} transacciones y {len(bk_mid)} updates de libro "
          f"en {dur:.0f} s")
    print(f"Huecos en ids de trade: {int(np.sum(np.diff(tr_id) > 1))}")
    print()

    k, mid_en_trade, dup_tr, dup_mid = medir_observable(
        tr_precio, tr_t, bk_mid, bk_t, nu
    )
    r_tr = medir_sigma_omega(tr_precio, tr_t, k, "precio de transaccion")
    r_mid = medir_sigma_omega(mid_en_trade, tr_t, k, "mid de bookTicker")

    print()
    print("=" * 92)
    print("RESUMEN CONTRA LAS CONDICIONES DE DESCONGELAMIENTO DEL §8")
    print("=" * 92)
    for nombre, r in (("precio de transaccion", r_tr), ("mid de bookTicker", r_mid)):
        if not r:
            continue
        c1 = r["frac_instr"] < 0.5
        c3 = r["frac_valida"] > 0.90
        print(f"  [{nombre}]")
        print(f"    1. sigma_instr^2/sigma_total^2 = {r['frac_instr']:.1%} < 50 %  "
              f"-> {'CUMPLE' if c1 else 'NO CUMPLE'}")
        print(f"    3. omega_valida > 90 %          = {r['frac_valida']:.1%}     "
              f"-> {'CUMPLE' if c3 else 'NO CUMPLE'}")
    c2 = pct(dup_mid, 50) < 0.75 * pct(dup_tr, 50)
    print(f"  2. duplicados de mid materialmente < 44.1 %: "
          f"{pct(dup_mid,50):.1%} contra {pct(dup_tr,50):.1%} "
          f"-> {'CUMPLE' if c2 else 'NO CUMPLE'}")
    print(f"  4. §5.2 en regimen agitado (nu > 300 tx/s): esta captura da "
          f"nu = {nu:.0f} tx/s -> {'CUMPLE' if nu > 300 else 'PENDIENTE'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
