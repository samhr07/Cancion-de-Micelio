// regularizador_dual_fractal.hpp -- resolvente dual L2/L1 del Bloque Fractal Universal, zero-heap (ARK-FRAC-01).
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-03.
//
// Los tres solvers son COPIA LITERAL de ark_cxx/include/ark/coordinacion/regularizador_dual.hpp (Fisico de Cabecera):
//   resolver_tikhonov_tridiagonal  (Thomas TDMA O(N), Tikhonov de primeras diferencias),
//   calibrar_lambda_morozov        (discrepancia de Morozov: residuo cuadratico medio = s0^2),
//   resolver_total_variation_condat (Condat 2013, O(N), exacto).
// La prueba C2 los compara bit a bit contra la cabecera del Fisico y C3 certifica Condat por las condiciones KKT
// (acta §1: 20 000 problemas, igual a la transcripcion del C de Condat a 7e-16, certificado 3e-12).
//
// Lo que NO se copia, por decision del Comandante del 2026-11-03 (acta §2):
//   - la compuerta Q de colas de los precios (calcular_ratio_q): la norma la elige la innovacion del propio observable;
//   - lambda_TV = b/ln 2 del Fisico y lambda_TV = s0 sqrt(2 ln N) del mandato: la TV deja plana una ventana si y solo si
//     max_k |sum_{i<=k}(y_i - media)| <= lambda (certificado dual), y ese maximo crece como s0*sqrt(N) y mas con ACF1 0.94;
//     con Donoho-Johnstone el 99-100 % de las ventanas en reposo se parte en 7-18 escalones. Aqui lambda_TV es el cuantil
//     1 - alfa de ese maximo bajo el nulo calibrado del observable (lo calcula el observador);
//   - el "filtro de piso nulo" del Fisico (recorte a [0, mu0] bajo mu0 + 0.5 s0): no tiene sentido para H.
//
// Causalidad: el resolvente ve solo la ventana [t - N + 1, t] y se extrae x[N-1], el borde presente. No reasigna fechas.
#pragma once

#include "ark/fractal/tipos_fractal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace ark::fractal {

// ---------- copia literal del Fisico (regularizador_dual.hpp, ARK-COORD-01-CXX) ----------

inline void resolver_tikhonov_tridiagonal(
    const double* y,
    std::size_t n,
    double lam,
    double* x_out,
    double* c_prime,
    double* d_prime
) noexcept {
    if (n <= 1 || lam <= 1e-12) {
        for (std::size_t i = 0; i < n; ++i) x_out[i] = y[i];
        return;
    }

    const double sub = -lam;
    const double d0 = 1.0 + lam;

    c_prime[0] = sub / d0;
    d_prime[0] = y[0] / d0;

    for (std::size_t i = 1; i < n - 1; ++i) {
        const double di = 1.0 + 2.0 * lam;
        const double denom = di - sub * c_prime[i - 1];
        c_prime[i] = sub / denom;
        d_prime[i] = (y[i] - sub * d_prime[i - 1]) / denom;
    }

    const double d_last = 1.0 + lam;
    const double denom_last = d_last - sub * c_prime[n - 2];
    d_prime[n - 1] = (y[n - 1] - sub * d_prime[n - 2]) / denom_last;

    x_out[n - 1] = d_prime[n - 1];
    for (std::size_t i = n - 2; i < n; --i) {
        x_out[i] = d_prime[i] - c_prime[i] * x_out[i + 1];
    }
}

inline double calibrar_lambda_morozov(
    const double* y,
    std::size_t n,
    double s0,
    double* x_out,
    double* c_prime,
    double* d_prime,
    double* work_cand
) noexcept {
    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) sum += y[i];
    const double mean = sum / static_cast<double>(n);

    double var_y = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double diff = y[i] - mean;
        var_y += diff * diff;
    }
    var_y /= static_cast<double>(n);
    const double target_res = s0 * s0;

    if (var_y <= target_res) {
        for (std::size_t i = 0; i < n; ++i) x_out[i] = mean;
        return 1e4;
    }

    double lo = 1e-4;
    double hi = 1e4;

    for (int iter = 0; iter < 22; ++iter) {
        const double mid = std::exp(0.5 * (std::log(lo) + std::log(hi)));
        resolver_tikhonov_tridiagonal(y, n, mid, work_cand, c_prime, d_prime);

        double res = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double diff = work_cand[i] - y[i];
            res += diff * diff;
        }
        res /= static_cast<double>(n);

        if (res < target_res) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    const double lam_opt = std::exp(0.5 * (std::log(lo) + std::log(hi)));
    resolver_tikhonov_tridiagonal(y, n, lam_opt, x_out, c_prime, d_prime);
    return lam_opt;
}

inline void resolver_total_variation_condat(
    const double* y,
    std::size_t n,
    double lam,
    double* x_out
) noexcept {
    if (n <= 1 || lam <= 1e-12) {
        for (std::size_t i = 0; i < n; ++i) x_out[i] = y[i];
        return;
    }

    std::size_t k = 0, k0 = 0;
    double vmin = y[0] - lam;
    double vmax = y[0] + lam;
    double umin = lam;
    double umax = -lam;
    std::size_t kplus = 0, kminus = 0;

    while (true) {
        while (k == n - 1) {
            if (umin < 0.0) {
                for (std::size_t i = k0; i <= kminus; ++i) x_out[i] = vmin;
                k = k0 = kminus + 1;
                if (k >= n) return;
                vmin = y[k];
                vmax = y[k] + 2.0 * lam;
                umin = lam;
                umax = -lam;
                kplus = kminus = k;
            } else if (umax > 0.0) {
                for (std::size_t i = k0; i <= kplus; ++i) x_out[i] = vmax;
                k = k0 = kplus + 1;
                if (k >= n) return;
                vmin = y[k] - 2.0 * lam;
                vmax = y[k];
                umin = lam;
                umax = -lam;
                kplus = kminus = k;
            } else {
                const double val = vmin + umin / static_cast<double>(k - k0 + 1);
                for (std::size_t i = k0; i < n; ++i) x_out[i] = val;
                return;
            }
        }

        if (y[k + 1] + umin < vmin - lam) {
            for (std::size_t i = k0; i <= kminus; ++i) x_out[i] = vmin;
            k = k0 = kminus + 1;
            vmin = y[k];
            vmax = y[k] + 2.0 * lam;
            umin = lam;
            umax = -lam;
            kplus = kminus = k;
        } else if (y[k + 1] + umax > vmax + lam) {
            for (std::size_t i = k0; i <= kplus; ++i) x_out[i] = vmax;
            k = k0 = kplus + 1;
            vmin = y[k] - 2.0 * lam;
            vmax = y[k];
            umin = lam;
            umax = -lam;
            kplus = kminus = k;
        } else {
            ++k;
            umin += y[k] - vmin;
            umax += y[k] - vmax;
            if (umin >= lam) {
                vmin += (umin - lam) / static_cast<double>(k - k0 + 1);
                umin = lam;
                kminus = k;
            }
            if (umax <= -lam) {
                vmax += (umax + lam) / static_cast<double>(k - k0 + 1);
                umax = -lam;
                kplus = k;
            }
        }
    }
}

// ---------- lo propio del Bloque Fractal ----------

// Certificado dual de la TV: la solucion es plana (un solo segmento) si y solo si max_{k < n-1} |S_k| <= lambda,
// con S_k = sum_{i<=k} (y_i - media). Es la cantidad cuyo cuantil de reposo ancla lambda_TV.
[[nodiscard]] inline double maximo_puente(const double* y, std::size_t n) noexcept {
    if (n < 2) return 0.0;
    double suma = 0.0;
    for (std::size_t i = 0; i < n; ++i) suma += y[i];
    const double media = suma / static_cast<double>(n);
    double s = 0.0, mx = 0.0;
    for (std::size_t i = 0; i + 1 < n; ++i) {
        s += y[i] - media;
        mx = std::max(mx, std::abs(s));
    }
    return mx;
}

struct SalidaRegularizador {
    double borde = NaN;      // x[N-1]
    double lambda = NaN;
};

// Resolvente dual zero-heap sobre la ventana y[0..n-1] (n <= MAX_CAP).
template <std::size_t MAX_CAP = 120>
class RegularizadorDualFractal {
    static_assert(MAX_CAP >= 8, "MAX_CAP >= 8");
    std::array<double, MAX_CAP> x_{}, c_{}, d_{}, w_{};

public:
    // L2: Tikhonov con lambda por Morozov sobre el piso de ruido s0 del observable.
    SalidaRegularizador l2(const double* y, std::size_t n, double s0) noexcept {
        n = std::min(n, MAX_CAP);
        if (n == 0) return {};
        const double lam = calibrar_lambda_morozov(y, n, s0, x_.data(), c_.data(), d_.data(), w_.data());
        return {x_[n - 1], lam};
    }
    // L1: TV exacta de Condat con lambda_TV calibrado (cuantil del puente en reposo).
    SalidaRegularizador l1(const double* y, std::size_t n, double lambda_tv) noexcept {
        n = std::min(n, MAX_CAP);
        if (n == 0) return {};
        resolver_total_variation_condat(y, n, lambda_tv, x_.data());
        return {x_[n - 1], lambda_tv};
    }
    [[nodiscard]] const std::array<double, MAX_CAP>& solucion() const noexcept { return x_; }
};

}  // namespace ark::fractal
