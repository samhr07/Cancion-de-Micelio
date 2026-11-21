#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verificar_horizonte_arfima_y_rachas.py

Auditoria cuantitativa reproducible para contrastar el §1 del dictamen de Argon
(1er_dictamen_enmienda_via_b_lilyth_Argon.md) frente a:
  1. La varianza exacta de truncado de innovaciones en ARFIMA(0, d, 0) mediante
     la autocovarianza cerrada de Hosking (1981) y la recurrencia de Durbin-Levinson.
  2. La longitud de racha util bajo corrupcion q en Z1b comparando corte duro
     en el primer NaN (r_max = 0, CORTE_SUFIJO) vs Puente Historico de Ito
     (r_max = 2, alfa_hueco = 0.01).
  3. El sesgo de nivel determinista por aplicar (1 - L)^d sin pre-centrar mu_serie.

Uso (Windows / Linux, solo biblioteca estandar Python 3):
  python verificar_horizonte_arfima_y_rachas.py
"""

from __future__ import annotations

import math


def pi_coefs(delta: float, N: int) -> list[float]:
    pi = [0.0] * (N + 1)
    pi[0] = 1.0
    for k in range(1, N + 1):
        pi[k] = pi[k - 1] * (k - 1.0 - delta) / k
    return pi


def hosking_acov(d: float, M: int) -> list[float]:
    """Autocovarianza relativa gamma_z(h)/sigma_eps^2 exacta de ARFIMA(0, d, 0)."""
    g0 = math.gamma(1.0 - 2.0 * d) / (math.gamma(1.0 - d) ** 2)
    g = [0.0] * (M + 1)
    g[0] = g0
    for h in range(1, M + 1):
        g[h] = g[h - 1] * (h - 1.0 + d) / (h - d)
    return g


def exact_pi_mse(d: float, M: int) -> float:
    """
    MSE relativo exacto Var(e_hat^(M)_t - eps_t) / sigma_eps^2 para e_hat^(M) = sum_{k=0}^M pi_k(d) z_{t-k}:
      Cov(e_hat^(M), eps_t) = pi_0 * psi_0 * sigma_eps^2 = sigma_eps^2
      => Var(e_hat^(M) - eps_t)/sigma_eps^2 = Var(e_hat^(M))/sigma_eps^2 - 1.
    Sin truncamiento asintotico.
    """
    pi = pi_coefs(d, M)
    g = hosking_acov(d, M)
    var = g[0] * sum(x * x for x in pi)
    for h in range(1, M + 1):
        s = 0.0
        for j in range(0, M - h + 1):
            s += pi[j] * pi[j + h]
        var += 2.0 * g[h] * s
    return var - 1.0


def levinson_durbin_mse(d: float, M: int) -> float:
    """Comprobacion numerica directa de Durbin-Levinson sobre gamma_z(0..M)."""
    g = hosking_acov(d, M)
    v = g[0]
    a: list[float] = []
    for k in range(1, M + 1):
        num = g[k]
        for j in range(1, k):
            num -= a[j - 1] * g[k - j]
        kappa = num / v
        a_new = [0.0] * k
        a_new[k - 1] = kappa
        for j in range(1, k):
            a_new[j - 1] = a[j - 1] - kappa * a[k - j - 1]
        a = a_new
        v *= (1.0 - kappa * kappa)
    return v - 1.0


def m_argon(delta: float, eps2: float, cap: int = 50000) -> tuple[int, bool]:
    """Criterio 1 - sum_{k=0}^M pi_k(delta)^2 / E(delta) <= eps2 de §1.2 de Argon."""
    E = math.gamma(1.0 + 2.0 * delta) / (math.gamma(1.0 + delta) ** 2)
    s = 1.0
    p = 1.0
    for k in range(1, cap + 1):
        if 1.0 - s / E <= eps2:
            return k - 1, False
        p = p * (k - 1.0 - delta) / k
        s += p * p
    return cap, True


def m_hosking_opt(d: float, eps2: float, cap: int = 50000) -> int:
    """Menor M tal que la varianza de innovacion optima finita v_M/sigma_eps^2 - 1 <= eps2."""
    v = math.gamma(1.0 - 2.0 * d) / (math.gamma(1.0 - d) ** 2)
    for k in range(1, cap + 1):
        if v - 1.0 <= eps2:
            return k - 1
        phi_kk = d / (k - d)
        v *= (1.0 - phi_kk * phi_kk)
    return cap


def m_pi_exact(d: float, eps2: float, max_m: int = 1500) -> int:
    """Menor M <= max_m tal que exact_pi_mse(d, M) <= eps2."""
    lo, hi = 1, max_m
    while lo < hi:
        mid = (lo + hi) // 2
        if exact_pi_mse(d, mid) <= eps2:
            hi = mid
        else:
            lo = mid + 1
    return lo if exact_pi_mse(d, lo) <= eps2 else -1


def main() -> None:
    z = 1.959963984540054  # norm_ppf(0.975)
    n = 4096
    alfa = 0.05
    eps2 = (z * z) / n

    print("================================================================================")
    print(" AUDITORIA DE HORIZONTE ARFIMA Y RACHAS BAJO NaN (RESPUESTA AL DICTAMEN ARGON)")
    print("================================================================================")
    print(f"Parametros: n = {n}, alfa = {alfa:.2f}, z_{{1-alfa/2}} = {z:.12f}, eps_Q^2 = {eps2:.8f}\n")

    print("--- SECCION 1: M_Argon(+d) en innovaciones() vs MSE exacto en ARFIMA(0, d, 0) ---")
    # Cross-check Levinson-Durbin vs Hosking closed form at M=64
    max_ld_diff = 0.0
    for d in (-0.35, -0.25, 0.25, 0.35):
        v_hosk = math.gamma(1.0 - 2.0 * d) / (math.gamma(1.0 - d) ** 2)
        for k in range(1, 65):
            v_hosk *= 1.0 - (d / (k - d)) ** 2
        mse_hosk = v_hosk - 1.0
        mse_ld = levinson_durbin_mse(d, 64)
        max_ld_diff = max(max_ld_diff, abs(mse_hosk - mse_ld))
    print(f"Verificacion cruzada Levinson-Durbin vs Hosking (M=64): max |dif| = {max_ld_diff:.3e}\n")

    hdr = (
        f"{'d':>6} | {'M_Arg(+d)':>10} | {'MSE(M_Arg)':>11} | {'x_epsQ2':>8} | "
        f"{'M_Arg(-d)':>10} | {'M_pi_exact':>10} | {'M_Hosk_opt':>10} | "
        f"{'MSE_pi(64)':>10} | {'MSE_pi(256)':>11}"
    )
    print(hdr)
    print("-" * len(hdr))
    for d in (-0.40, -0.35, -0.25, -0.10, +0.10, +0.25, +0.35, +0.40):
        ma_pos, cap_pos = m_argon(+d, eps2)
        ma_neg, cap_neg = m_argon(-d, eps2)
        m_pi = m_pi_exact(d, eps2, max_m=1500)
        m_hk = m_hosking_opt(d, eps2)
        mse64 = exact_pi_mse(d, 64)
        mse256 = exact_pi_mse(d, 256)
        if not cap_pos and ma_pos <= 1500:
            mse_arg = exact_pi_mse(d, ma_pos)
            mse_arg_s = f"{mse_arg:11.6f}"
            rat_s = f"{mse_arg / eps2:7.2f}x"
        else:
            mse_arg_s = f"{'n/a(>cap)':>11}"
            rat_s = f"{'n/a':>8}"
        s_pos = f">={ma_pos}" if cap_pos else f"{ma_pos}"
        s_neg = f">={ma_neg}" if cap_neg else f"{ma_neg}"
        print(
            f"{d:+6.2f} | {s_pos:>10} | {mse_arg_s} | {rat_s} | "
            f"{s_neg:>10} | {m_pi:10d} | {m_hk:10d} | "
            f"{mse64:10.6f} | {mse256:11.6f}"
        )

    print("\n--- SECCION 2: Racha util en Z1b: break en 1er NaN (r_max=0) vs Puente Ito (r_max=2, alfa_hueco=0.01) ---")
    alfa_h = 0.01
    hdr2 = (
        f"{'q':>6} | {'E[racha|r=0] iid':>16} | {'E[racha|Ito] iid':>16} | {'Ganancia':>8} | "
        f"{'P(r>=85) r=0 / Ito':>19} | {'P(r>=233) r=0 / Ito':>20} | {'P(r>=418) r=0 / Ito':>20} | "
        f"{'Rafaga r=2: r=0 vs Ito':>22}"
    )
    print(hdr2)
    print("-" * len(hdr2))
    for q in (0.005, 0.010, 0.020, 0.050):
        e_r0 = (1.0 - q) / q
        # En Bernoulli i.i.d.(q), un corte no puenteable ocurre si:
        #   - arranca una racha >= 3 de NaNs: prob (1-q)*q^3, o
        #   - un micro-hueco valido de longitud 1 o 2 dispara falso veto |Z_hueco| > z_max (alfa_h = 0.01)
        p_cut_iid = (1.0 - q) * (q ** 3) + alfa_h * (1.0 - q) * (q + q * q)
        e_ito = 1.0 / p_cut_iid
        gan = e_ito / e_r0

        def surv(M: int) -> str:
            p0 = (1.0 - q) ** M
            pi_ito = (1.0 - p_cut_iid) ** M
            return f"{p0:6.3f} / {pi_ito:6.3f}"

        # Modelo de rafagas aisladas r=2 de Nivel 1 / Fase 4 (tasa de inicio p_g = q / 2):
        p_burst_r0 = q / 2.0
        p_burst_ito = (q / 2.0) * alfa_h
        e_burst_0 = 1.0 / p_burst_r0
        e_burst_ito = 1.0 / p_burst_ito
        print(
            f"{q*100:5.1f}% | {e_r0:16.1f} | {e_ito:16.1f} | {gan:7.1f}x | "
            f"{surv(85):>19} | {surv(233):>20} | {surv(418):>20} | "
            f"{e_burst_0:6.0f} vs {e_burst_ito:7.0f}"
        )

    print("\n--- SECCION 3: Confirmacion cuantitativa de §1.1 (sesgo de nivel sin pre-centrar mu_serie) ---")
    mu_z = 50.0
    N_cal = 4096
    print(f"Escenario: mu_z = {mu_z:.1f}, N_cal = {N_cal}")
    for d in (-0.35, +0.25, +0.35):
        pi_dm1_cal = pi_coefs(d - 1.0, N_cal)  # sum_{k=0}^t pi_k(d) = pi_t(d - 1)
        mean_bias_cal = sum(pi_dm1_cal[:N_cal]) / N_cal
        for n_hist in (64, 256):
            # En el paso j=0 del perfil tras n_hist barras previas:
            bias_j0 = pi_dm1_cal[n_hist]
            delta_e0 = mu_z * (bias_j0 - mean_bias_cal)
            print(
                f"  d={d:+.2f}, n_hist={n_hist:3d}: pi_{n_hist}(d-1)={bias_j0:+.6f}, "
                f"prom_Ncal={mean_bias_cal:+.6f}, sesgo residual Delta_e0 = {delta_e0:+.6f}"
            )


if __name__ == "__main__":
    main()
