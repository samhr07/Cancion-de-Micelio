// modelo_lineal_nulo.hpp -- Operador de Lilyth, ETAPA 2: modelo lineal del nulo de calma (AR(p) por Levinson-Durbin con orden de
// Quenouille y ARFIMA(p, d, 0) opcional con d de Whittle exacto), filtro inverso A(L)(1-L)^d para extraer las innovaciones de un
// perfil real, y simulacion generativa con o sin inyeccion de innovaciones. C++20.
//
// Autoria: Argon (Auditoria Externa), sobre el plano e interfaces de Code (Jefe de Arquitectos de Codigo, Opus 5.5), 2026-11-06.
// Especificacion: PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md (E1.19, E1.22, §8.2, §9.2) y PLANO_ETAPA2_MODELO_LINEAL_NULO.md.
//
// Reglas cumplidas:
//   1. El modelo se ajusta sobre la serie que Lilyth mide en pasos de la barra base, con su propio K por est::tau_integrado
//      (tope k_max = n / 4, igual al primer termino de canal_nivel.hpp).
//   2. Momentos y ACF mediante est::momentos y est::acf_agrupada (sesgada, normaliza por n: matriz de Toeplitz definida positiva).
//   3. Levinson-Durbin sobre rho_1..rho_K y seleccion del orden p = mayor k <= K con |kappa_k| > z_{1-alfa/2} / sqrt(n) (Quenouille).
//   4. Residuos e_t = (x_t - mu) - sum_{i=1}^p phi_i (x_{t-i} - mu) solo para t >= p dentro de cada segmento (nunca cruza un corte).
//   5. ARFIMA(p, d, 0): exige d_externa finita con |d| < 1/2; diferencia cada segmento con elw::coeficientes(d, ...) (tipo II) antes
//      de estimar el AR(p), y en la simulacion integra con elw::coeficientes(-d, ...).
//   6. Cero reservas de memoria dinamica en `innovaciones` y `simular` (ciclo de evaluacion/simulacion).
#pragma once

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"
#include "ark/metrologia/whittle_exacto.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace ark::lilyth {

enum class InnovacionesNulo : std::uint8_t { GAUSS = 0, RESIDUOS_REMUESTREADOS = 1 };

struct ConfiguracionModeloLineal {
    ModoNuloLineal modo = ModoNuloLineal::AR_LEVINSON_DURBIN;
    double alfa = NaN;                                    // el de la familia; obligatorio
    InnovacionesNulo innovaciones = InnovacionesNulo::RESIDUOS_REMUESTREADOS;
    double d_externa = NaN;                               // ARFIMA: ResultadoELW.d de la serie modelada
};

struct ModeloLinealNulo {
    bool valido = false;
    const char* motivo = "sin ajustar";
    ModoNuloLineal modo_usado = ModoNuloLineal::AR_LEVINSON_DURBIN;
    InnovacionesNulo tipo_innovaciones = InnovacionesNulo::RESIDUOS_REMUESTREADOS;
    bool cayo_a_generativo = false;                       // E1.19: se pidio FASE_ALEATORIZADA_FFT
    std::size_t n = 0, K = 0, p = 0;
    double mu = NaN, sd = NaN, sigma2_innov = NaN, d = 0.0;
    std::vector<double> phi;                              // phi_1..phi_p
    std::vector<double> kappa;                            // kappa_1..kappa_K (todos, para la ficha)
    std::vector<double> residuos;                         // innovaciones de la calma (para RESIDUOS_REMUESTREADOS)
    double media_residuos = 0.0;                          // media de `residuos` para centrar en el remuestreo
    double radio_max = NaN;                               // |raiz| mayor de la recurrencia (inversa de las raices de A(z))
    std::size_t burn_in = 0;
};

namespace detalle_lineal {

inline constexpr std::size_t BUF_STACK_MAX = 16384;
inline constexpr std::size_t P_STACK_MAX = 1024;

// Modulo de la mayor raiz del polinomio caracteristico de la recurrencia z^p - phi_1 z^{p-1} - ... - phi_p = 0
// (equivalente al inverso de la menor raiz de A(L) = 1 - phi_1 L - ... - phi_p L^p).
// Para p >= 3 usa iteracion simultanea de Aberth-Ehrlich (con suma de reciprocos en lugar de producto, libre de overflow/underflow
// incluso para ordenes p > 100).
inline double calcular_radio_max(const std::vector<double>& phi) {
    const std::size_t p = phi.size();
    if (p == 0) return 0.0;
    if (p == 1) return std::abs(phi[0]);
    if (p == 2) {
        const double disc = phi[0] * phi[0] + 4.0 * phi[1];
        if (disc >= 0.0) return 0.5 * (std::abs(phi[0]) + std::sqrt(disc));
        return std::sqrt(std::max(0.0, -phi[1]));
    }
    const double r0 = (std::abs(phi.back()) > 1e-12)
                          ? std::pow(std::abs(phi.back()), 1.0 / static_cast<double>(p))
                          : 0.85;
    std::vector<std::complex<double>> z(p);
    for (std::size_t k = 0; k < p; ++k) {
        const double ang = 2.0 * ark::metrologia::est::PI * static_cast<double>(k) / static_cast<double>(p) + 0.13579;
        z[k] = std::polar(r0, ang);
    }
    for (int it = 0; it < 80; ++it) {
        double max_w = 0.0;
        for (std::size_t k = 0; k < p; ++k) {
            const std::complex<double> zk = z[k];
            std::complex<double> ratio{0.0, 0.0};
            if (std::abs(zk) >= 0.5) {
                const std::complex<double> u = 1.0 / zk;
                std::complex<double> A{1.0, 0.0};
                std::complex<double> B{static_cast<double>(p), 0.0};
                std::complex<double> up = u;
                for (std::size_t i = 1; i <= p; ++i) {
                    A -= phi[i - 1] * up;
                    B -= static_cast<double>(p - i) * phi[i - 1] * up;
                    up *= u;
                }
                if (std::abs(B) > 1e-30) ratio = zk * (A / B);
            } else {
                std::complex<double> P{1.0, 0.0};
                std::complex<double> dP{0.0, 0.0};
                for (std::size_t i = 0; i < p; ++i) {
                    dP = dP * zk + P;
                    P = P * zk - phi[i];
                }
                if (std::abs(dP) > 1e-30) ratio = P / dP;
            }
            std::complex<double> suma_inv{0.0, 0.0};
            for (std::size_t j = 0; j < p; ++j) {
                if (j == k) continue;
                const std::complex<double> dif = zk - z[j];
                if (std::abs(dif) > 1e-15) suma_inv += 1.0 / dif;
            }
            const std::complex<double> denom = 1.0 - ratio * suma_inv;
            const std::complex<double> w = (std::abs(denom) > 1e-15) ? (ratio / denom) : ratio;
            z[k] -= w;
            max_w = std::max(max_w, std::abs(w));
        }
        if (max_w <= 1e-13) break;
    }
    double rmax = 0.0;
    for (const auto& r : z) {
        if (std::isfinite(r.real()) && std::isfinite(r.imag())) {
            rmax = std::max(rmax, std::abs(r));
        }
    }
    return rmax;
}

}  // namespace detalle_lineal

// Tamano de buffer externo opcional para `simular` cuando se usa ARFIMA en paralelo (burn_in + n para pi y burn_in + n para u).
[[nodiscard]] inline std::size_t tamano_trabajo_simular(const ModeloLinealNulo& m, std::size_t n) noexcept {
    if (m.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE && m.d != 0.0) {
        return 2 * (m.burn_in + n);
    }
    return m.p;
}

// Fuera del ciclo caliente: ajusta el modelo lineal sobre los segmentos contiguos de calma.
inline ModeloLinealNulo ajustar(const std::vector<std::vector<double>>& segs, const ConfiguracionModeloLineal& cfg) {
    ModeloLinealNulo m{};
    m.modo_usado = cfg.modo;
    m.tipo_innovaciones = cfg.innovaciones;

    if (!(std::isfinite(cfg.alfa) && cfg.alfa > 0.0 && cfg.alfa < 1.0)) {
        m.motivo = "alfa obligatorio en (0, 1)";
        return m;
    }

    if (cfg.modo == ModoNuloLineal::FASE_ALEATORIZADA_FFT) {
        // E1.19: FASE_ALEATORIZADA_FFT solo sirve para el regimen espontaneo (a); para el modelo generativo cae a AR.
        m.cayo_a_generativo = true;
        m.modo_usado = ModoNuloLineal::AR_LEVINSON_DURBIN;
    } else if (cfg.modo == ModoNuloLineal::ARFIMA_WHITTLE) {
        if (!std::isfinite(cfg.d_externa)) {
            m.motivo = "ARFIMA sin d_externa finita";
            return m;
        }
        if (!(std::abs(cfg.d_externa) < 0.5)) {
            m.motivo = "ARFIMA con |d_externa| >= 0.5 (no estacionario o no invertible)";
            return m;
        }
        m.d = cfg.d_externa;
    }

    // Extraer tramos contiguos finitos sin mezclar a traves de cortes.
    std::vector<std::vector<double>> limpios;
    std::size_t max_seg = 0;
    for (const auto& g : segs) {
        std::vector<double> cur;
        for (double v : g) {
            if (std::isfinite(v)) {
                cur.push_back(v);
            } else if (!cur.empty()) {
                max_seg = std::max(max_seg, cur.size());
                limpios.push_back(std::move(cur));
                cur.clear();
            }
        }
        if (!cur.empty()) {
            max_seg = std::max(max_seg, cur.size());
            limpios.push_back(std::move(cur));
        }
    }
    if (limpios.empty()) {
        m.motivo = "sin segmentos validos";
        return m;
    }

    // Si es ARFIMA(p, d, 0), diferenciar cada segmento con elw::coeficientes(d, ...) truncado al inicio del segmento (tipo II).
    std::vector<std::vector<double>> segs_modelo;
    if (m.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE && m.d != 0.0) {
        std::vector<double> pi(max_seg, 0.0);
        ark::metrologia::elw::coeficientes(m.d, max_seg, pi.data());
        segs_modelo.reserve(limpios.size());
        for (const auto& g : limpios) {
            std::vector<double> gd(g.size(), 0.0);
            for (std::size_t t = 0; t < g.size(); ++t) {
                double s = 0.0;
                for (std::size_t k = 0; k <= t; ++k) {
                    s += pi[k] * g[t - k];
                }
                gd[t] = s;
            }
            segs_modelo.push_back(std::move(gd));
        }
    } else {
        segs_modelo = limpios;
    }

    const auto mo = ark::metrologia::est::momentos(segs_modelo);
    m.n = mo.n;
    m.mu = mo.mu;
    m.sd = mo.sd;
    if (m.n < 8 || !(m.sd > 0.0) || !std::isfinite(m.sd)) {
        m.motivo = "historia insuficiente o varianza nula";
        return m;
    }

    const std::size_t k_max = m.n / 4;
    if (k_max < 1) {
        m.motivo = "historia corta (n / 4 < 1)";
        return m;
    }
    const auto ti = ark::metrologia::est::tau_integrado(segs_modelo, m.mu, m.sd, cfg.alfa, k_max);
    if (!(ti.tau > 0.0) || !std::isfinite(ti.tau)) {
        m.motivo = "tau_int <= 0 o no finito";
        return m;
    }
    if (ti.K > k_max) {
        m.K = k_max;
        m.radio_max = 1.0;
        m.motivo = "raiz unitaria: la ACF no entra en la banda del nulo en n/4 rezagos";
        return m;
    }
    m.K = ti.K;

    std::vector<double> rho(m.K + 1, 0.0);
    rho[0] = 1.0;
    for (std::size_t k = 1; k <= m.K; ++k) {
        std::size_t pares = 0;
        rho[k] = ark::metrologia::est::acf_agrupada(segs_modelo, m.mu, m.sd, k, pares);
        if (!std::isfinite(rho[k])) {
            m.motivo = "ACF no finita en Levinson-Durbin";
            return m;
        }
    }

    // Levinson-Durbin sobre rho_1..rho_K y orden p por Quenouille (el MAYOR k <= K con |kappa_k| > z_{1-alfa/2} / sqrt(n)).
    const double umbral_q = ark::metrologia::est::norm_ppf(1.0 - 0.5 * cfg.alfa) / std::sqrt(static_cast<double>(m.n));
    double v = 1.0;
    double v_en_p = 1.0;
    std::vector<double> phi_prev, phi_k;
    m.kappa.assign(m.K, NaN);
    m.p = 0;
    m.phi.clear();

    for (std::size_t k = 1; k <= m.K; ++k) {
        double num = 0.0;
        for (std::size_t j = 1; j < k; ++j) {
            num += phi_prev[j - 1] * rho[k - j];
        }
        if (!(v > 0.0)) {
            m.radio_max = 1.0;
            m.motivo = "varianza de innovacion no positiva en Levinson-Durbin";
            return m;
        }
        const double kappa_k = (rho[k] - num) / v;
        if (!std::isfinite(kappa_k) || std::abs(kappa_k) >= 1.0) {
            m.radio_max = 1.0;
            m.motivo = "raiz unitaria o estacionariedad rota (|kappa_k| >= 1)";
            return m;
        }
        m.kappa[k - 1] = kappa_k;
        phi_k.assign(k, 0.0);
        phi_k[k - 1] = kappa_k;
        for (std::size_t j = 1; j < k; ++j) {
            phi_k[j - 1] = phi_prev[j - 1] - kappa_k * phi_prev[k - j - 1];
        }
        v *= (1.0 - kappa_k * kappa_k);
        phi_prev = phi_k;

        if (std::abs(kappa_k) > umbral_q) {
            m.p = k;
            m.phi = phi_k;
            v_en_p = v;
        }
    }

    m.sigma2_innov = m.sd * m.sd * v_en_p;
    if (!(m.sigma2_innov > 0.0) || !std::isfinite(m.sigma2_innov)) {
        m.motivo = "sigma2_innov no positiva";
        return m;
    }

    m.radio_max = detalle_lineal::calcular_radio_max(m.phi);
    if (!std::isfinite(m.radio_max) || m.radio_max >= 1.0) {
        m.motivo = "raiz unitaria o inestable (radio_max >= 1)";
        return m;
    }

    if (m.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE) {
        m.burn_in = max_seg;
    } else if (m.p > 0 && m.radio_max > 0.0) {
        const double eps = std::numeric_limits<double>::epsilon();
        m.burn_in = static_cast<std::size_t>(std::ceil(std::log(eps) / std::log(m.radio_max)));
    } else {
        m.burn_in = 0;
    }

    // Residuos de la calma: e_t = (x_t - mu) - sum_{i=1}^p phi_i (x_{t-i} - mu), solo con t >= p dentro de cada segmento.
    m.residuos.clear();
    double suma_res = 0.0;
    for (const auto& g : segs_modelo) {
        if (g.size() <= m.p) continue;
        for (std::size_t t = m.p; t < g.size(); ++t) {
            double e = g[t] - m.mu;
            for (std::size_t i = 1; i <= m.p; ++i) {
                e -= m.phi[i - 1] * (g[t - i] - m.mu);
            }
            m.residuos.push_back(e);
            suma_res += e;
        }
    }
    if (m.residuos.empty()) {
        m.motivo = "ningun segmento supera el orden p para extraer residuos";
        return m;
    }
    m.media_residuos = suma_res / static_cast<double>(m.residuos.size());
    m.valido = true;
    m.motivo = m.p == 0 ? "ajustado (p = 0, sin memoria AR)" : "ajustado";
    return m;
}

// Filtro inverso A(L)(1 - L)^d sobre el perfil, con condiciones iniciales tomadas de la historia previa. Sin reservas.
// e_out con l entradas. historia: la serie justo antes del perfil (al menos p valores; en ARFIMA, toda la que haya).
inline void innovaciones(const ModeloLinealNulo& m, const double* historia, std::size_t n_hist,
                         const double* perfil, std::size_t l, double* e_out) noexcept {
    if (!e_out || l == 0) return;
    if (!m.valido || !perfil) {
        for (std::size_t j = 0; j < l; ++j) e_out[j] = NaN;
        return;
    }
    if (!historia) n_hist = 0;

    auto z_en = [&](std::size_t idx) noexcept -> double {
        return (idx < n_hist) ? historia[idx] : perfil[idx - n_hist];
    };

    if (m.modo_usado != ModoNuloLineal::ARFIMA_WHITTLE || m.d == 0.0) {
        for (std::size_t j = 0; j < l; ++j) {
            const std::size_t t = n_hist + j;
            double e = z_en(t) - m.mu;
            for (std::size_t i = 1; i <= m.p; ++i) {
                if (t >= i) {
                    e -= m.phi[i - 1] * (z_en(t - i) - m.mu);
                }
            }
            e_out[j] = e;
        }
        return;
    }

    // Caso ARFIMA(p, d, 0): primero u_t = ((1 - L)^d z)_t con elw::coeficientes(d, ...) y luego A(L) sobre (u_t - mu).
    const std::size_t n_tot = n_hist + l;
    std::array<double, detalle_lineal::BUF_STACK_MAX> pi_stack{};
    const std::size_t n_pi = std::min(n_tot, detalle_lineal::BUF_STACK_MAX);
    ark::metrologia::elw::coeficientes(m.d, n_pi, pi_stack.data());

    auto u_en = [&](std::size_t t) noexcept -> double {
        double s = 0.0;
        const std::size_t lim = std::min(t + 1, n_pi);
        for (std::size_t k = 0; k < lim; ++k) {
            s += pi_stack[k] * z_en(t - k);
        }
        if (t + 1 > n_pi) {
            double pik = pi_stack[n_pi - 1];
            for (std::size_t k = n_pi; k <= t; ++k) {
                pik *= (static_cast<double>(k) - 1.0 - m.d) / static_cast<double>(k);
                s += pik * z_en(t - k);
            }
        }
        return s;
    };

    for (std::size_t j = 0; j < l; ++j) {
        const std::size_t t = n_hist + j;
        double e = u_en(t) - m.mu;
        for (std::size_t i = 1; i <= m.p; ++i) {
            if (t >= i) {
                e -= m.phi[i - 1] * (u_en(t - i) - m.mu);
            }
        }
        e_out[j] = e;
    }
}

// Una realizacion de n pasos. Si e_iny no es nulo, en t_iny..t_iny+l-1 se usan esas innovaciones en lugar de las sorteadas;
// despues el filtro sigue solo. x_out con n entradas; el burn-in va antes y no se escribe. Sin reservas en el ciclo.
template <class Rng>
inline void simular(const ModeloLinealNulo& m, Rng& rng, const double* e_iny, std::size_t t_iny, std::size_t l,
                    double* x_out, std::size_t n, double* trabajo = nullptr, std::size_t n_trabajo = 0) noexcept {
    if (!x_out || n == 0) return;
    if (!m.valido) {
        for (std::size_t t = 0; t < n; ++t) x_out[t] = NaN;
        return;
    }

    const double sd_innov = std::sqrt(std::max(0.0, m.sigma2_innov));
    std::normal_distribution<double> dist_gauss(0.0, sd_innov);
    const std::size_t n_res = m.residuos.size();
    std::uniform_int_distribution<std::size_t> dist_idx(0, n_res > 0 ? n_res - 1 : 0);

    auto sortear_e = [&]() noexcept -> double {
        if (m.tipo_innovaciones == InnovacionesNulo::RESIDUOS_REMUESTREADOS && n_res > 0) {
            return m.residuos[dist_idx(rng)] - m.media_residuos;
        }
        return dist_gauss(rng);
    };

    // Si e_iny cubre desde t_iny == 0 toda la salida n (control L4 de ida y vuelta desde el origen), no se antepone burn-in
    // aleatorio para que x_out[0..n-1] sea funcion determinista pura de e_iny[0..n-1]. En cualquier otro caso corre m.burn_in.
    const std::size_t burn = (e_iny != nullptr && t_iny == 0 && l >= n) ? 0 : m.burn_in;
    const std::size_t p = m.p;

    if (m.modo_usado != ModoNuloLineal::ARFIMA_WHITTLE || m.d == 0.0) {
        // Caso AR(p): anillo para mantener los ultimos p valores centrados durante el burn-in.
        std::array<double, detalle_lineal::P_STACK_MAX> ring_stack{};
        double* ring = (trabajo && n_trabajo >= p) ? trabajo : ring_stack.data();
        const std::size_t cap_ring = (trabajo && n_trabajo >= p) ? p : std::min(p, detalle_lineal::P_STACK_MAX);
        for (std::size_t i = 0; i < cap_ring; ++i) ring[i] = 0.0;

        for (std::size_t s = 0; s < burn; ++s) {
            const double e = sortear_e();
            double u_c = e;
            for (std::size_t i = 1; i <= cap_ring; ++i) {
                if (s >= i) {
                    u_c += m.phi[i - 1] * ring[(s - i) % cap_ring];
                }
            }
            if (cap_ring > 0) ring[s % cap_ring] = u_c;
        }

        for (std::size_t t = 0; t < n; ++t) {
            const bool inyectar = (e_iny != nullptr) && (t >= t_iny) && (t < t_iny + l);
            const double e = inyectar ? e_iny[t - t_iny] : sortear_e();
            double u_c = e;
            for (std::size_t i = 1; i <= p; ++i) {
                if (t >= i) {
                    u_c += m.phi[i - 1] * (x_out[t - i] - m.mu);
                } else if (burn + t >= i && cap_ring > 0 && i <= cap_ring) {
                    const std::size_t s_prev = burn + t - i;
                    u_c += m.phi[i - 1] * ring[s_prev % cap_ring];
                }
            }
            x_out[t] = m.mu + u_c;
        }
        return;
    }

    // Caso ARFIMA(p, d, 0): primero se genera la serie estacionaria de memoria corta u_s (con su AR(p)) y luego se integra
    // con (1 - L)^(-d) usando elw::coeficientes(-m.d, ...).
    const std::size_t tot = burn + n;
    std::array<double, detalle_lineal::BUF_STACK_MAX> pi_stack{};
    std::array<double, detalle_lineal::BUF_STACK_MAX> u_stack{};

    double* pi_ptr = nullptr;
    double* u_ptr = nullptr;
    std::size_t cap_buf = 0;
    if (trabajo && n_trabajo >= 2 * tot) {
        pi_ptr = trabajo;
        u_ptr = trabajo + tot;
        cap_buf = tot;
    } else {
        pi_ptr = pi_stack.data();
        u_ptr = u_stack.data();
        cap_buf = std::min(tot, detalle_lineal::BUF_STACK_MAX);
    }

    ark::metrologia::elw::coeficientes(-m.d, cap_buf, pi_ptr);

    for (std::size_t s = 0; s < cap_buf; ++s) {
        const bool en_salida = (s >= burn);
        const std::size_t t = en_salida ? (s - burn) : 0;
        const bool inyectar = en_salida && (e_iny != nullptr) && (t >= t_iny) && (t < t_iny + l);
        const double e = inyectar ? e_iny[t - t_iny] : sortear_e();
        double u_c = e;
        for (std::size_t i = 1; i <= p; ++i) {
            if (s >= i) {
                u_c += m.phi[i - 1] * (u_ptr[s - i] - m.mu);
            }
        }
        u_ptr[s] = m.mu + u_c;
        if (en_salida && t < n) {
            double suma = 0.0;
            for (std::size_t k = 0; k <= s; ++k) {
                suma += pi_ptr[k] * u_ptr[s - k];
            }
            x_out[t] = suma;
        }
    }
}

}  // namespace ark::lilyth