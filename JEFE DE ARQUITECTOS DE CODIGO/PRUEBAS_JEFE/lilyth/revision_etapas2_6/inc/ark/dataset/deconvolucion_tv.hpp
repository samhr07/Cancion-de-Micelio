// deconvolucion_tv.hpp -- calca C++ de ARQUITECTO DE DATOS/deconvolucion_tv/motor_deconvolucion_tv.py (resolver_deconvolucion_tv_admm),
// la rama TV de motor_ventanas.desacoplar_rampa_convolucional_tv. Mismos pasos y mismas reglas; el producto de matrices y la
// factorizacion de Cholesky son los de un solver directo, no los de LAPACK/BLAS de numpy-scipy: aqui NO se exige identidad bit a
// bit sino tolerancia declarada (el contraste lo mide). Jefe de Arquitectos, 2026-11-03.
#pragma once

#include "ark/dataset/errores_py.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace ark::dataset::ventanas {

ARK_ERROR_PY(ErrorDeconvolucionTV, ErrorPy)
ARK_ERROR_PY(ErrorDimensionInsuficienteTV, ErrorDeconvolucionTV)
ARK_ERROR_PY(ErrorVarianzaDegeneradaTV, ErrorDeconvolucionTV)
ARK_ERROR_PY(ErrorParametroInvalidoTV, ErrorDeconvolucionTV)

enum class KernelDeconvolucion : uint8_t { RECTANGULAR_BOXCAR, EXPONENCIAL, TRIANGULAR_BARTLETT };

struct ConfiguracionTV {
    double lambda_tv = 0.40;
    double rho_admm = 0.50;
    int64_t max_iter = 100;
    double tol_convergencia = 1e-5;
    KernelDeconvolucion kernel_tipo = KernelDeconvolucion::RECTANGULAR_BOXCAR;
    int64_t ancho_kernel_W = 30;
    double parametro_decae_exp = 0.1;
    double umbral_salto_significativo = 0.10;
    int64_t min_muestras = 40;
    bool ajuste_estado_estacionario = true;
};

struct ResultadoTV {
    std::vector<double> serie_deconvolucionada;
    int64_t iteraciones = 0;
};

namespace detalle_tv {

// suma por pares de numpy (la que usa np.sum en las filas de A y del kernel)
inline double suma_pares(const double* a, size_t n) {
    if (n < 8) { double r = 0.0; for (size_t k = 0; k < n; ++k) r += a[k]; return r; }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; ++j) r[j] = a[j];
        size_t k = 8;
        for (; k < n - (n % 8); k += 8) for (int j = 0; j < 8; ++j) r[j] += a[k + j];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; k < n; ++k) res += a[k];
        return res;
    }
    size_t n2 = n / 2; n2 -= n2 % 8;
    return suma_pares(a, n2) + suma_pares(a + n2, n - n2);
}

inline std::vector<double> generar_kernel_causal(int64_t W, KernelDeconvolucion tipo, double param_exp) {
    if (W < 2) throw ErrorParametroInvalidoTV("ancho_kernel_W=" + std::to_string(W) + " debe ser >= 2.");
    std::vector<double> k(static_cast<size_t>(W));
    if (tipo == KernelDeconvolucion::RECTANGULAR_BOXCAR) {
        for (auto& x : k) x = 1.0 / static_cast<double>(W);
    } else if (tipo == KernelDeconvolucion::EXPONENCIAL) {
        for (int64_t i = 0; i < W; ++i) k[static_cast<size_t>(i)] = std::exp(-param_exp * static_cast<double>(i));
        const double s = 0.0 + suma_pares(k.data(), k.size());
        for (auto& x : k) x /= s;
    } else {
        for (int64_t i = 0; i < W; ++i) k[static_cast<size_t>(i)] = static_cast<double>(W) - static_cast<double>(i);
        const double s = 0.0 + suma_pares(k.data(), k.size());
        for (auto& x : k) x /= s;
    }
    return k;
}

}  // namespace detalle_tv

// min_x 0.5 ||A x - y||^2 + lambda ||D x||_1 por ADMM (Boyd et al. 2011), con las filas interiores A[W-1:, :].
inline ResultadoTV resolver_deconvolucion_tv_admm(const std::vector<double>& y, const ConfiguracionTV& cfg) {
    const int64_t N = static_cast<int64_t>(y.size());
    const int64_t W = cfg.ancho_kernel_W;
    if (cfg.lambda_tv <= 0.0) throw ErrorParametroInvalidoTV("lambda_tv debe ser estrictamente positivo (> 0).");
    if (cfg.rho_admm <= 0.0) throw ErrorParametroInvalidoTV("rho_admm debe ser estrictamente positivo (> 0).");
    if (cfg.max_iter < 1) throw ErrorParametroInvalidoTV("max_iter debe ser >= 1.");
    if (N < cfg.min_muestras || N <= W)
        throw ErrorDimensionInsuficienteTV("Longitud N=" + std::to_string(N) + " insuficiente para ventana W=" + std::to_string(W) +
                                           " y min_muestras=" + std::to_string(cfg.min_muestras) + ".");
    for (double v : y) if (!std::isfinite(v)) throw ErrorVarianzaDegeneradaTV("La serie observada contiene valores no finitos (NaN o Inf).");
    double mn = y[0], mx = y[0];
    for (double v : y) { if (v < mn) mn = v; if (v > mx) mx = v; }
    if (std::fabs(mn - mx) <= 1e-12 + 1e-5 * std::fabs(mx))   // np.isclose(min, max, atol=1e-12)
        throw ErrorVarianzaDegeneradaTV("Varianza nula en la serie observada (valor constante).");

    const size_t n = static_cast<size_t>(N);
    // A (N x N) Toeplitz triangular inferior y el ajuste de estado estacionario de sus primeras filas
    const auto ker = detalle_tv::generar_kernel_causal(W, cfg.kernel_tipo, cfg.parametro_decae_exp);
    std::vector<double> A(n * n, 0.0);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j <= i; ++j) {
            const size_t lag = i - j;
            if (lag < ker.size()) A[i * n + j] = ker[lag];
        }
    if (cfg.ajuste_estado_estacionario) {
        for (size_t i = 0; i < std::min<size_t>(n, static_cast<size_t>(W - 1)); ++i) {
            const double s = 0.0 + detalle_tv::suma_pares(&A[i * n], n);
            A[i * n] += 1.0 - s;
        }
    }
    const size_t f0 = static_cast<size_t>(W - 1), nf = n - f0;   // filas A_val
    const double rho = cfg.rho_admm;
    // M = A_val^T A_val + rho D^T D ;  Aty = A_val^T y_val
    std::vector<double> M(n * n, 0.0), Aty(n, 0.0);
    for (size_t r = 0; r < nf; ++r) {
        const double* fila = &A[(f0 + r) * n];
        const double yr = y[f0 + r];
        for (size_t a = 0; a < n; ++a) {
            if (fila[a] == 0.0) continue;
            Aty[a] += fila[a] * yr;
            for (size_t b = 0; b < n; ++b) M[a * n + b] += fila[a] * fila[b];
        }
    }
    for (size_t i = 0; i + 1 < n; ++i) {   // D^T D: tridiagonal (1, 2, ..., 2, 1) y -1 fuera de la diagonal
        M[i * n + i] += rho;
        M[(i + 1) * n + (i + 1)] += rho;
        M[i * n + (i + 1)] -= rho;
        M[(i + 1) * n + i] -= rho;
    }
    // Cholesky inferior
    std::vector<double> L(n * n, 0.0);
    for (size_t j = 0; j < n; ++j) {
        double s = M[j * n + j];
        for (size_t k = 0; k < j; ++k) s -= L[j * n + k] * L[j * n + k];
        if (!(s > 0.0)) throw ErrorParametroInvalidoTV("Matriz M no invertible en ADMM.");
        const double d = std::sqrt(s);
        L[j * n + j] = d;
        for (size_t i = j + 1; i < n; ++i) {
            double t = M[i * n + j];
            for (size_t k = 0; k < j; ++k) t -= L[i * n + k] * L[j * n + k];
            L[i * n + j] = t / d;
        }
    }
    auto resolver = [&](const std::vector<double>& rhs) {
        std::vector<double> z(n), x(n);
        for (size_t i = 0; i < n; ++i) {
            double s = rhs[i];
            for (size_t k = 0; k < i; ++k) s -= L[i * n + k] * z[k];
            z[i] = s / L[i * n + i];
        }
        for (size_t ii = n; ii-- > 0;) {
            double s = z[ii];
            for (size_t k = ii + 1; k < n; ++k) s -= L[k * n + ii] * x[k];
            x[ii] = s / L[ii * n + ii];
        }
        return x;
    };
    auto Dt = [&](const std::vector<double>& v) {   // D^T v, v de largo n-1
        std::vector<double> o(n, 0.0);
        for (size_t i = 0; i + 1 < n; ++i) { o[i] -= v[i]; o[i + 1] += v[i]; }
        return o;
    };
    std::vector<double> z(n - 1, 0.0), u(n - 1, 0.0), x = y;
    const double umbral = cfg.lambda_tv / rho;
    int64_t it_ej = 0;
    for (int64_t it = 0; it < cfg.max_iter; ++it) {
        ++it_ej;
        std::vector<double> zu(n - 1);
        for (size_t i = 0; i + 1 < n; ++i) zu[i] = z[i] - u[i];
        const auto dtz = Dt(zu);
        std::vector<double> rhs(n);
        for (size_t i = 0; i < n; ++i) rhs[i] = Aty[i] + rho * dtz[i];
        const auto x_new = resolver(rhs);
        std::vector<double> dx(n - 1), z_new(n - 1), rp(n - 1);
        for (size_t i = 0; i + 1 < n; ++i) {
            dx[i] = x_new[i + 1] - x_new[i];
            const double v = dx[i] + u[i];
            const double sg = v > 0.0 ? 1.0 : (v < 0.0 ? -1.0 : 0.0);
            z_new[i] = sg * std::max(std::fabs(v) - umbral, 0.0);
            rp[i] = dx[i] - z_new[i];
        }
        std::vector<double> dz(n - 1);
        for (size_t i = 0; i + 1 < n; ++i) dz[i] = z_new[i] - z[i];
        const auto sd = Dt(dz);
        for (size_t i = 0; i + 1 < n; ++i) u[i] += rp[i];
        z = z_new;
        double rpm = 0.0, rdm = 0.0, dxm = 0.0;
        for (double v : rp) rpm = std::max(rpm, std::fabs(v));
        for (double v : sd) rdm = std::max(rdm, std::fabs(rho * v));
        for (size_t i = 0; i < n; ++i) dxm = std::max(dxm, std::fabs(x_new[i] - x[i]));
        x = x_new;
        if (std::max({rpm, rdm, dxm}) < cfg.tol_convergencia && it >= 5) break;
    }
    return {std::move(x), it_ej};
}

}  // namespace ark::dataset::ventanas
