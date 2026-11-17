// superestadistica_inversa.hpp -- inversion de Fredholm p = K f (ARK-METRO-01). Diagnostico fuera del ciclo caliente.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04. Prerregistro: PRUEBAS_JEFE/metrologia/EXPECTATIVAS.md §1.
//
// Modelo: p(y) = sum_j f_j N(y; 0, 1/beta_j), f_j >= 0, sum f_j = 1 (masas sobre una rejilla de M <= 32 betas log-espaciados).
// Datos: histograma EQUIPROBABLE con K = min(K_MAX, ceil(2 n^(1/3))) casillas (Rice) y extremos abiertos; nucleo integrado por
// casilla: K_ij = Phi(e_{i+1} sqrt(beta_j)) - Phi(e_i sqrt(beta_j)). Peso W_i = 1/sigma_i, sigma_i^2 = p^_i (1 - p^_i) / n.
// Penalizacion sobre las primeras diferencias de f (decision del Comandante: la misma norma que reconstruyo la serie):
//   L1_TV:       lambda ||D f||_1        (con f >= 0 y sum f = 1, ||f||_1 es constante: el L1 con sentido es sobre D f)
//   L2_TIKHONOV: lambda/2 ||D f||_2^2
// Solver: Chambolle-Pock (2011) con A = [W K ; D], G = indicadora del simplex (proyeccion exacta). lambda por discrepancia:
// chi^2 = ||W (K f - p^)||^2 = K (declarado, Morozov en conteos). Exporta ademas los valores singulares de W K (identificabilidad).
#pragma once

#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace ark::metrologia {

template <std::size_t K_MAX = 64, std::size_t M_MAX = 32>
struct ResultadoInversion {
    NormaRegularizador norma = NormaRegularizador::SIN_DATO;
    std::size_t K = 0, M = 0, n = 0;
    std::array<double, M_MAX> beta{}, f{};
    double lambda = NaN, chi2 = NaN, chi2_objetivo = NaN;
    bool discrepancia_alcanzable = false;
    double media_log_beta = NaN, var_log_beta = NaN, beta_media = NaN, q_tsallis = NaN;
    std::array<double, M_MAX> valores_singulares{};
    std::size_t iteraciones = 0;
};

namespace detalle_inv {

// proyeccion euclidea sobre el simplex {f >= 0, sum f = 1} (Held, Wolfe y Crowder; Duchi et al. 2008).
template <std::size_t M_MAX>
inline void proyectar_simplex(std::array<double, M_MAX>& v, std::size_t M) noexcept {
    std::array<double, M_MAX> u{};
    for (std::size_t i = 0; i < M; ++i) u[i] = v[i];
    std::sort(u.begin(), u.begin() + static_cast<std::ptrdiff_t>(M), [](double a, double b) { return a > b; });
    double cs = 0.0, theta = 0.0;
    for (std::size_t i = 0; i < M; ++i) {
        cs += u[i];
        const double t = (cs - 1.0) / static_cast<double>(i + 1);
        if (u[i] - t > 0.0) theta = t;
    }
    for (std::size_t i = 0; i < M; ++i) v[i] = std::max(0.0, v[i] - theta);
}

// valores propios de una simetrica M x M por Jacobi ciclico
template <std::size_t M_MAX>
inline void jacobi(std::array<std::array<double, M_MAX>, M_MAX> a, std::size_t M, std::array<double, M_MAX>& ev) noexcept {
    for (int barrido = 0; barrido < 100; ++barrido) {
        double off = 0.0;
        for (std::size_t p = 0; p < M; ++p) for (std::size_t q = p + 1; q < M; ++q) off += a[p][q] * a[p][q];
        if (off < 1e-30) break;
        for (std::size_t p = 0; p < M; ++p)
            for (std::size_t q = p + 1; q < M; ++q) {
                if (std::abs(a[p][q]) < 1e-300) continue;
                const double th = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
                const double t = (th >= 0 ? 1.0 : -1.0) / (std::abs(th) + std::sqrt(th * th + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (std::size_t k = 0; k < M; ++k) {
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < M; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
    }
    for (std::size_t i = 0; i < M; ++i) ev[i] = a[i][i];
}

}  // namespace detalle_inv

template <std::size_t K_MAX = 64, std::size_t M_MAX = 32>
class InversionSuperestadistica {
    std::size_t K_ = 0, M_ = 0, n_ = 0;
    std::array<double, M_MAX> beta_{};
    std::array<std::array<double, M_MAX>, K_MAX> WK_{};   // W K
    std::array<double, K_MAX> b_{};                       // W p^
    double normaA_ = 0.0;
    double c_ = 1.0;   // ||W K||: el bloque de datos se resuelve dividido por c_ (equilibra los dos bloques de A); lambda' = lambda / c_^2

    void aplicar_A(const std::array<double, M_MAX>& f, std::array<double, K_MAX>& u, std::array<double, M_MAX>& v) const noexcept {
        for (std::size_t i = 0; i < K_; ++i) { double s = 0.0; for (std::size_t j = 0; j < M_; ++j) s += WK_[i][j] * f[j]; u[i] = s; }
        for (std::size_t j = 0; j + 1 < M_; ++j) v[j] = f[j + 1] - f[j];
    }
    void aplicar_At(const std::array<double, K_MAX>& u, const std::array<double, M_MAX>& v, std::array<double, M_MAX>& f) const noexcept {
        for (std::size_t j = 0; j < M_; ++j) { double s = 0.0; for (std::size_t i = 0; i < K_; ++i) s += WK_[i][j] * u[i]; f[j] = s; }
        for (std::size_t j = 0; j + 1 < M_; ++j) { f[j + 1] += v[j]; f[j] -= v[j]; }
    }

public:
    double chi2(const std::array<double, M_MAX>& f) const noexcept {
        double c = 0.0;
        for (std::size_t i = 0; i < K_; ++i) { double s = -b_[i]; for (std::size_t j = 0; j < M_; ++j) s += WK_[i][j] * f[j]; c += s * s; }
        return c * c_ * c_;
    }

    // y: datos ya centrados (y - mu0). Rejilla de beta: M_MAX puntos log-espaciados en [beta_lo, beta_hi].
    bool preparar(std::vector<double> y, double beta_lo, double beta_hi) {
        y.erase(std::remove_if(y.begin(), y.end(), [](double v) { return !std::isfinite(v); }), y.end());
        n_ = y.size();
        if (n_ < 16 || !(beta_lo > 0.0) || !(beta_hi > beta_lo)) return false;
        std::sort(y.begin(), y.end());
        K_ = std::min<std::size_t>(K_MAX, static_cast<std::size_t>(std::ceil(2.0 * std::cbrt(static_cast<double>(n_)))));
        M_ = M_MAX;
        std::array<double, K_MAX + 1> e{};
        e[0] = -INFINITY; e[K_] = INFINITY;
        for (std::size_t i = 1; i < K_; ++i) e[i] = y[i * n_ / K_];
        std::array<double, K_MAX> p{};
        for (std::size_t i = 0; i < K_; ++i) {
            const auto lo = i == 0 ? y.begin() : std::lower_bound(y.begin(), y.end(), e[i]);
            const auto hi = i + 1 == K_ ? y.end() : std::lower_bound(y.begin(), y.end(), e[i + 1]);
            p[i] = static_cast<double>(hi - lo) / static_cast<double>(n_);
        }
        for (std::size_t j = 0; j < M_; ++j)
            beta_[j] = beta_lo * std::pow(beta_hi / beta_lo, static_cast<double>(j) / static_cast<double>(M_ - 1));
        for (std::size_t i = 0; i < K_; ++i) {
            const double sp = std::max(p[i] * (1.0 - p[i]), 1.0 / static_cast<double>(n_)) / static_cast<double>(n_);
            const double w = 1.0 / std::sqrt(sp);
            b_[i] = w * p[i];
            for (std::size_t j = 0; j < M_; ++j) {
                const double sb = std::sqrt(beta_[j]);
                WK_[i][j] = w * (est::phi_cdf(e[i + 1] * sb) - est::phi_cdf(e[i] * sb));
            }
        }
        // c_ = ||W K|| (potencia sobre (WK)^T (WK)); se reescala el bloque de datos para que los dos bloques de A sean comparables
        {
            std::array<double, M_MAX> f{}, g{};
            f.fill(1.0);
            double nr = 0.0;
            for (int it = 0; it < 300; ++it) {
                for (std::size_t j = 0; j < M_; ++j) g[j] = 0.0;
                for (std::size_t i = 0; i < K_; ++i) {
                    double ui = 0.0;
                    for (std::size_t j = 0; j < M_; ++j) ui += WK_[i][j] * f[j];
                    for (std::size_t j = 0; j < M_; ++j) g[j] += WK_[i][j] * ui;
                }
                double q = 0.0;
                for (std::size_t j = 0; j < M_; ++j) q += g[j] * g[j];
                nr = std::sqrt(q);
                for (std::size_t j = 0; j < M_; ++j) f[j] = g[j] / nr;
            }
            c_ = std::sqrt(nr);
            for (std::size_t i = 0; i < K_; ++i) { b_[i] /= c_; for (std::size_t j = 0; j < M_; ++j) WK_[i][j] /= c_; }
        }
        // ||A|| por potencia sobre A^T A
        std::array<double, M_MAX> f{}, g{};
        std::array<double, K_MAX> u{};
        std::array<double, M_MAX> v{};
        f.fill(1.0);
        double nrm = 0.0;
        for (int it = 0; it < 200; ++it) {
            aplicar_A(f, u, v); aplicar_At(u, v, g);
            double s = 0.0;
            for (std::size_t j = 0; j < M_; ++j) s += g[j] * g[j];
            s = std::sqrt(s);
            nrm = s;
            for (std::size_t j = 0; j < M_; ++j) f[j] = g[j] / s;
        }
        normaA_ = std::sqrt(nrm);
        return true;
    }

    // resolver para un lambda dado (Chambolle-Pock)
    std::size_t resolver(NormaRegularizador norma, double lambda, std::array<double, M_MAX>& f, std::size_t iter_max = 200000) const noexcept {
        const double tau = 0.99 / normaA_, sig = 0.99 / normaA_;
        lambda /= c_ * c_;   // misma solucion que el problema sin reescalar con este lambda
        std::array<double, K_MAX> y1{}, u{};
        std::array<double, M_MAX> y2{}, v{}, fb{}, at{}, fn{};
        for (std::size_t j = 0; j < M_; ++j) f[j] = 1.0 / static_cast<double>(M_);
        fb = f;
        std::size_t it = 0;
        for (; it < iter_max; ++it) {
            aplicar_A(fb, u, v);
            for (std::size_t i = 0; i < K_; ++i) y1[i] = (y1[i] + sig * u[i] - sig * b_[i]) / (1.0 + sig);
            for (std::size_t j = 0; j + 1 < M_; ++j) {
                const double w = y2[j] + sig * v[j];
                y2[j] = norma == NormaRegularizador::L1_TV ? std::clamp(w, -lambda, lambda) : w / (1.0 + sig / lambda);
            }
            aplicar_At(y1, y2, at);
            for (std::size_t j = 0; j < M_; ++j) fn[j] = f[j] - tau * at[j];
            detalle_inv::proyectar_simplex<M_MAX>(fn, M_);
            double dif = 0.0;
            for (std::size_t j = 0; j < M_; ++j) { dif = std::max(dif, std::abs(fn[j] - f[j])); fb[j] = 2.0 * fn[j] - f[j]; f[j] = fn[j]; }
            if (dif < 1e-12 && it > 50) break;
        }
        return it;
    }

    ResultadoInversion<K_MAX, M_MAX> invertir(NormaRegularizador norma) const {
        ResultadoInversion<K_MAX, M_MAX> r{};
        r.norma = norma; r.K = K_; r.M = M_; r.n = n_;
        r.beta = beta_;
        r.chi2_objetivo = static_cast<double>(K_);
        // biseccion en log lambda: chi2 crece con lambda
        double lo = std::log(1e-9), hi = std::log(1e9);
        std::array<double, M_MAX> f{};
        resolver(norma, std::exp(lo), f);
        if (chi2(f) > r.chi2_objetivo) {   // ni sin penalizacion se alcanza el objetivo: se reporta el minimo
            r.lambda = std::exp(lo);
            r.discrepancia_alcanzable = false;
        } else {
            for (int it = 0; it < 40; ++it) {
                const double mid = 0.5 * (lo + hi);
                resolver(norma, std::exp(mid), f);
                (chi2(f) > r.chi2_objetivo ? hi : lo) = mid;
            }
            r.lambda = std::exp(lo);
            r.discrepancia_alcanzable = true;
        }
        r.iteraciones = resolver(norma, r.lambda, f);
        r.f = f;
        r.chi2 = chi2(f);
        double m1 = 0, m2 = 0, b1 = 0, b2 = 0;
        for (std::size_t j = 0; j < M_; ++j) {
            const double lb = std::log(beta_[j]);
            m1 += f[j] * lb; m2 += f[j] * lb * lb; b1 += f[j] * beta_[j]; b2 += f[j] * beta_[j] * beta_[j];
        }
        r.media_log_beta = m1;
        r.var_log_beta = m2 - m1 * m1;
        r.beta_media = b1;
        r.q_tsallis = b2 / (b1 * b1);
        // valores singulares de W K
        std::array<std::array<double, M_MAX>, M_MAX> AtA{};
        for (std::size_t p = 0; p < M_; ++p)
            for (std::size_t q = 0; q < M_; ++q) { double s = 0.0; for (std::size_t i = 0; i < K_; ++i) s += WK_[i][p] * WK_[i][q]; AtA[p][q] = s; }
        std::array<double, M_MAX> ev{};
        detalle_inv::jacobi<M_MAX>(AtA, M_, ev);
        std::sort(ev.begin(), ev.begin() + static_cast<std::ptrdiff_t>(M_), [](double a, double b) { return a > b; });
        for (std::size_t j = 0; j < M_; ++j) r.valores_singulares[j] = c_ * std::sqrt(std::max(0.0, ev[j]));
        return r;
    }
};

}  // namespace ark::metrologia
