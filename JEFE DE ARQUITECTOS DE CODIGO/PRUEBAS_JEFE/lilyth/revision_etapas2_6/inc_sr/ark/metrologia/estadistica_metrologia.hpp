// estadistica_metrologia.hpp -- utilidades numericas del bloque ARK-METRO-01 (fuera del ciclo caliente salvo lo marcado noexcept).
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04.
#pragma once

#include "ark/metrologia/tipos_metrologia.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace ark::metrologia::est {

inline constexpr double PI = 3.14159265358979323846;

// Phi(x), CDF normal estandar.
inline double phi_cdf(double x) noexcept { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// Inversa de Phi: Acklam (2003) + una iteracion de Halley con erfc (error ~1e-15).
inline double norm_ppf(double p) noexcept {
    if (!(p > 0.0 && p < 1.0)) return p == 0.0 ? -INFINITY : (p == 1.0 ? INFINITY : NaN);
    static constexpr double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                                   1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static constexpr double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                                   6.680131188771972e+01, -1.328068155288572e+01};
    static constexpr double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                                   -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static constexpr double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
    const double pl = 0.02425;
    double x;
    if (p < pl) {
        const double q = std::sqrt(-2.0 * std::log(p));
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    } else if (p <= 1.0 - pl) {
        const double q = p - 0.5, r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
    } else {
        const double q = std::sqrt(-2.0 * std::log(1.0 - p));
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
    }
    const double e = phi_cdf(x) - p;
    const double u = e * std::sqrt(2.0 * PI) * std::exp(0.5 * x * x);
    return x - u / (1.0 + 0.5 * x * u);
}

// psi(x) y psi1(x) (digamma, trigamma) para x > 0: recurrencia hasta x >= 8 y serie asintotica.
inline double digamma(double x) noexcept {
    double r = 0.0;
    while (x < 8.0) { r -= 1.0 / x; x += 1.0; }
    const double f = 1.0 / (x * x);
    return r + std::log(x) - 0.5 / x - f * (1.0 / 12.0 - f * (1.0 / 120.0 - f * (1.0 / 252.0 - f * (1.0 / 240.0 - f / 132.0))));
}
inline double trigamma(double x) noexcept {
    double r = 0.0;
    while (x < 16.0) { r += 1.0 / (x * x); x += 1.0; }
    const double f = 1.0 / (x * x);
    return r + 1.0 / x + f / 2.0 + (1.0 / (x * x * x)) * (1.0 / 6.0 - f * (1.0 / 30.0 - f * (1.0 / 42.0 - f / 30.0)));
}

// Cuantil lineal (numpy.quantile method="linear"). Copia la entrada.
inline double cuantil(std::vector<double> v, double q) {
    if (v.empty()) return NaN;
    std::sort(v.begin(), v.end());
    const double pos = q * static_cast<double>(v.size() - 1);
    const std::size_t i = static_cast<std::size_t>(std::floor(pos));
    const double fr = pos - static_cast<double>(i);
    return i + 1 < v.size() ? v[i] + fr * (v[i + 1] - v[i]) : v.back();
}

// Momentos agrupados de varios segmentos contiguos: media, sd (n-1), ACF de rezago k agrupada (pares dentro de cada segmento).
struct Momentos { double mu = NaN, sd = NaN; std::size_t n = 0; };
inline Momentos momentos(const std::vector<std::vector<double>>& segs) {
    Momentos m{};
    double s = 0.0;
    for (const auto& g : segs) { for (double v : g) s += v; m.n += g.size(); }
    if (m.n < 2) return m;
    m.mu = s / static_cast<double>(m.n);
    double ss = 0.0;
    for (const auto& g : segs) for (double v : g) ss += (v - m.mu) * (v - m.mu);
    m.sd = std::sqrt(ss / static_cast<double>(m.n - 1));
    return m;
}
inline double acf_agrupada(const std::vector<std::vector<double>>& segs, double mu, double sd, std::size_t k, std::size_t& pares) {
    double num = 0.0;
    pares = 0;
    std::size_t n = 0;
    for (const auto& g : segs) {
        n += g.size();
        for (std::size_t i = k; i < g.size(); ++i) { num += (g[i] - mu) * (g[i - k] - mu); ++pares; }
    }
    if (pares == 0 || !(sd > 0.0)) return NaN;
    return num / (static_cast<double>(n) * sd * sd);   // estimador sesgado clasico (normaliza por n): estable en rezagos largos
}

// tau_int = 1 + 2 sum_{k=1}^{K-1} rho_k, K = primer rezago con |rho_k| <= z_{1-alfa/2}/sqrt(n) (dentro de la banda del nulo).
struct TauInt { double tau = NaN; double acf1 = NaN; std::size_t K = 0; };
inline TauInt tau_integrado(const std::vector<std::vector<double>>& segs, double mu, double sd, double alfa, std::size_t k_max) {
    TauInt t{};
    std::size_t n = 0;
    for (const auto& g : segs) n += g.size();
    if (n < 4) return t;
    const double banda = norm_ppf(1.0 - 0.5 * alfa) / std::sqrt(static_cast<double>(n));
    double suma = 0.0;
    for (std::size_t k = 1; k <= k_max; ++k) {
        std::size_t pares = 0;
        const double r = acf_agrupada(segs, mu, sd, k, pares);
        if (k == 1) t.acf1 = r;
        if (!std::isfinite(r) || std::abs(r) <= banda) { t.K = k; break; }
        suma += r;
        t.K = k + 1;
    }
    t.tau = 1.0 + 2.0 * suma;
    return t;
}

}  // namespace ark::metrologia::est
