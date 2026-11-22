// ============================================================================
// puente_historico_ito.hpp -- Reconstrucción Histórica por Puente Estocástico
// de Itô / Hoffman-Ribak (AR / ARFIMA / Browniano) con Doble Llave Derivada
// desde alfa, Curva de Potencia Mínima Detectable y Combinación Múltiple
// de Rubin (Fase 4 / Piloto Nivel 1). C++20. Cero reservas dinámicas.
// ============================================================================
#pragma once

#include "ark/dataset/tipos_dataset.hpp"
#include "ark/lilyth/modelo_lineal_nulo.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"
#include "ark/metrologia/whittle_exacto.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <vector>

namespace ark::lilyth::piloto {

inline constexpr std::size_t MAX_RACHA_PUENTE_ITO = 16;
inline constexpr std::size_t MAX_SEMILLAS_PUENTE_ITO = 32;
inline constexpr std::size_t MAX_HIST_CAUSAL_PUENTE = 512;

struct ConfiguracionPuenteHistoricoIto {
    ark::metrologia::Integracion integracion = ark::metrologia::Integracion::NIVEL;
    std::size_t r_max = 2;                 // Llave temporal: tope operativo derivado en NIVEL (r <= 2; en Z1b puede fijarse explícitamente)
    double alfa_hueco = 0.01;              // Nivel bilateral bajo H0 gaussiano (z_max = Phi^{-1}(1 - alfa_hueco/2))
    double z_max_hueco = ark::metrologia::NaN; // Si es NaN (por omisión), se deriva de alfa_hueco; cero constantes a ojo
    std::size_t n_semillas_rng = 16;       // Número de realizaciones RNG para medir dispersión inicio-fin y Rubin
    bool usar_recta_lineal_pura = false;   // true = recta determinista sin ruido; false = Puente de Itô/Schur
    std::uint64_t semilla_base = 20261118ULL;
};

struct ResultadoHuecoIto {
    bool reconstruido = false;
    bool excedio_r_max = false;
    bool choque_en_hueco = false;
    std::size_t r_hueco = 0;               // Barras faltantes r
    std::size_t delta_t = 0;               // r + 1
    double x_ini = ark::metrologia::NaN;
    double x_fin = ark::metrologia::NaN;
    double delta_x_obs = 0.0;              // X_fin - X_ini
    double delta_x_esp = 0.0;              // E[X_fin - X_ini | pasado]
    double sigma_dt = 0.0;                 // Desviación estándar condicional de (X_fin - X_ini) en dt
    double z_hueco = 0.0;                  // (delta_x_obs - delta_x_esp) / sigma_dt
    double z_max_efectivo = 0.0;           // Umbral derivado de alfa_hueco (o explícito si se proveyó)
    double delta_min_detectable_sigma_e = 0.0; // Salto mínimo (en sigma_e) detectable con potencia 0.80 para r_hueco
    std::uint8_t flags_emitidos = 0;       // FLAG_PUENTE_ITO, FLAG_CHOQUE_EN_HUECO o FLAG_PAUSA_UMBRAL

    // Trayectoria condicional media (suavizador exacto de 2 bordes) para k = 0..r-1
    std::array<double, MAX_RACHA_PUENTE_ITO> x_media_cond{};
    std::array<double, MAX_RACHA_PUENTE_ITO> y_est_media_cond{};

    // M realizaciones estocásticas del puente de Itô/ARFIMA (en niveles o flujo crudo y en representación estacionaria)
    std::size_t n_semillas_efectivas = 0;
    std::array<std::array<double, MAX_RACHA_PUENTE_ITO>, MAX_SEMILLAS_PUENTE_ITO> x_rng{};
    std::array<std::array<double, MAX_RACHA_PUENTE_ITO>, MAX_SEMILLAS_PUENTE_ITO> y_est_rng{};

    // Auditoría de inicio-fin del RNG (punto 2 del Comandante)
    double max_error_anclaje_rng = 0.0;    // max_m |sum_{s=1}^{r+1} dY^(m)_s - (X_fin - X_ini)|
    double sd_rng_media_trayectoria = 0.0; // Dispersión media entre las M semillas en el interior del hueco
};

// Resultado combinado de M imputaciones mediante reglas de Rubin (1987) y grados de libertad de Barnard-Rubin.
struct ResumenRubinIto {
    bool valido = false;
    std::size_t m_imputaciones = 0;
    double s_barra = ark::metrologia::NaN;            // Media de S^(m)
    double w_barra = ark::metrologia::NaN;            // Varianza media intra-imputación W_bar = mean(ee_m^2)
    double b_entre = ark::metrologia::NaN;            // Varianza entre imputaciones B_M
    double t_total = ark::metrologia::NaN;            // T_M = W_bar + (1 + 1/M) * B_M
    double ee_rubin = ark::metrologia::NaN;           // sqrt(T_M)
    double nu_barnard_rubin = ark::metrologia::NaN;   // Grados de libertad efectivos de Student-t
    double gamma_info_perdida = ark::metrologia::NaN; // Fracción de información perdida (1 + 1/M) * B_M / T_M
    double estadistico_t = ark::metrologia::NaN;      // s_barra / ee_rubin
    double p_bilateral = ark::metrologia::NaN;        // p bilateral exacto bajo t_{nu_M}
};

// Deriva el umbral |Z_hueco| desde alfa_hueco (Phi^{-1}(1 - alfa_hueco / 2)), sin números mágicos.
[[nodiscard]] inline double derivar_z_max_hueco(const ConfiguracionPuenteHistoricoIto& cfg) noexcept {
    if (std::isfinite(cfg.z_max_hueco) && cfg.z_max_hueco > 0.0) {
        return cfg.z_max_hueco;
    }
    const double a = (std::isfinite(cfg.alfa_hueco) && cfg.alfa_hueco > 0.0 && cfg.alfa_hueco < 1.0)
                         ? cfg.alfa_hueco
                         : 0.01;
    return ark::metrologia::est::norm_ppf(1.0 - 0.5 * a);
}

namespace detalle_ito {

[[nodiscard]] inline std::uint64_t splitmix64_ito(std::uint64_t x) noexcept {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Fracción continua de Lentz modificada para la función beta incompleta regularizada I_x(a, b).
[[nodiscard]] inline double betacf_ito(double a, double b, double x) noexcept {
    constexpr int MAXIT = 200;
    constexpr double EPS = 3.0e-14;
    constexpr double FPMIN = 1.0e-300;
    const double qab = a + b;
    const double qap = a + 1.0;
    const double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::abs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= MAXIT; ++m) {
        const double dm = static_cast<double>(m);
        const double m2 = 2.0 * dm;
        double aa = dm * (b - dm) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + dm) * (qab + dm) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::abs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (std::abs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (std::abs(del - 1.0) <= EPS) break;
    }
    return h;
}

[[nodiscard]] inline double betai_ito(double a, double b, double x) noexcept {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double bt = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                               a * std::log(x) + b * std::log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0)) {
        return bt * betacf_ito(a, b, x) / a;
    }
    return 1.0 - bt * betacf_ito(b, a, 1.0 - x) / b;
}

[[nodiscard]] inline double p_bilateral_student_t(double t_abs, double nu) noexcept {
    if (!std::isfinite(t_abs)) return ark::metrologia::NaN;
    if (t_abs <= 0.0) return 1.0;
    if (!std::isfinite(nu) || nu > 1.0e8) {
        return std::erfc(t_abs / std::sqrt(2.0));
    }
    const double x = nu / (nu + t_abs * t_abs);
    return std::clamp(betai_ito(0.5 * nu, 0.5, x), 0.0, 1.0);
}

// Calcula la respuesta al impulso de Wold psi[0..n_irf-1] del ModeloLinealNulo:
//   A(L) (1 - L)^d (z_t - mu) = eps_t  =>  z_t - mu = sum_{k>=0} psi_k eps_{t-k}, con psi_0 = 1.
inline void calcular_irf_wold(
    const ModeloLinealNulo& m,
    std::size_t n_irf,
    double* psi_out) noexcept {
    if (n_irf == 0) return;
    const std::size_t lim = std::min(n_irf, MAX_RACHA_PUENTE_ITO + 4);
    std::array<double, MAX_RACHA_PUENTE_ITO + 4> u_irf{};
    u_irf[0] = 1.0;
    if (m.valido && m.p > 0) {
        for (std::size_t k = 1; k < lim; ++k) {
            double acc = 0.0;
            for (std::size_t i = 1; i <= m.p && i <= k; ++i) {
                acc += m.phi[i - 1] * u_irf[k - i];
            }
            u_irf[k] = acc;
        }
    }
    if (!m.valido || m.modo_usado != ModoNuloLineal::ARFIMA_WHITTLE || m.d == 0.0) {
        for (std::size_t k = 0; k < lim; ++k) psi_out[k] = u_irf[k];
        return;
    }
    std::array<double, MAX_RACHA_PUENTE_ITO + 4> pi_inv{};
    ark::metrologia::elw::coeficientes(-m.d, lim, pi_inv.data());
    for (std::size_t k = 0; k < lim; ++k) {
        double s = 0.0;
        for (std::size_t j = 0; j <= k; ++j) {
            s += pi_inv[j] * u_irf[k - j];
        }
        psi_out[k] = s;
    }
}

// Predicción causal a 1 paso en representación estacionaria z_t usando historia previa buf[0..n_buf-1]
[[nodiscard]] inline double predecir_estacionario_un_paso(
    const ModeloLinealNulo& m,
    const double* buf,
    std::size_t n_buf) noexcept {
    if (!m.valido || n_buf == 0) {
        if (m.valido && m.enmienda_nivel2 && m.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE &&
            m.d != 0.0 && std::isfinite(m.mu_serie)) {
            return m.mu_serie + (std::isfinite(m.mu) ? m.mu : 0.0);
        }
        return (m.valido && std::isfinite(m.mu)) ? m.mu : 0.0;
    }
    if (m.modo_usado != ModoNuloLineal::ARFIMA_WHITTLE || m.d == 0.0) {
        double pred = m.mu;
        for (std::size_t i = 1; i <= m.p && i <= n_buf; ++i) {
            pred += m.phi[i - 1] * (buf[n_buf - i] - m.mu);
        }
        return pred;
    }
    std::array<double, MAX_HIST_CAUSAL_PUENTE> pi_loc{};
    std::array<double, MAX_HIST_CAUSAL_PUENTE> u_loc{};
    const std::size_t n_ef = std::min(n_buf + 1, pi_loc.size());
    const std::size_t offset = (n_buf + 1 > n_ef) ? (n_buf + 1 - n_ef) : 0;
    const double* sub = buf + offset;
    const std::size_t n_sub = n_buf - offset;

    const bool centrar_n2 = m.enmienda_nivel2 && std::isfinite(m.mu_serie);
    const double mu_pre = centrar_n2 ? m.mu_serie : 0.0;

    ark::metrologia::elw::coeficientes(m.d, n_ef, pi_loc.data());
    for (std::size_t t = 0; t < n_sub; ++t) {
        double s = 0.0;
        if (centrar_n2) {
            for (std::size_t k = 0; k <= t && k < n_ef; ++k) s += pi_loc[k] * (sub[t - k] - mu_pre);
        } else {
            for (std::size_t k = 0; k <= t && k < n_ef; ++k) s += pi_loc[k] * sub[t - k];
        }
        u_loc[t] = s;
    }
    double u_pred = m.mu;
    for (std::size_t i = 1; i <= m.p && i <= n_sub; ++i) {
        u_pred += m.phi[i - 1] * (u_loc[n_sub - i] - m.mu);
    }
    double cola = 0.0;
    if (centrar_n2) {
        for (std::size_t k = 1; k <= n_sub && k < n_ef; ++k) {
            cola += pi_loc[k] * (sub[n_sub - k] - mu_pre);
        }
        return mu_pre + u_pred - cola;
    }
    for (std::size_t k = 1; k <= n_sub && k < n_ef; ++k) {
        cola += pi_loc[k] * sub[n_sub - k];
    }
    return u_pred - cola;
}

} // namespace detalle_ito

// Curva analítica de potencia mínima detectable:
// Devuelve cuántos sigma_e necesita un choque oculto elemental para alcanzar `potencia_objetivo` en un hueco de longitud r_hueco.
// Si `escalas_fase` tiene longitud >= r_hueco + 1, pondera cada paso s = 0..r_hueco por su factor diurno causal kappa_s = sqrt(v_fase / v_media).
[[nodiscard]] inline double delta_minima_detectable_ito(
    const ConfiguracionPuenteHistoricoIto& cfg,
    const ModeloLinealNulo& m,
    std::size_t r_hueco,
    double potencia_objetivo = 0.80,
    std::span<const double> escalas_fase = {}) noexcept {
    if (r_hueco == 0) return 0.0;
    const std::size_t dt = r_hueco + 1;
    const std::size_t dt_cap = std::min(dt, MAX_RACHA_PUENTE_ITO + 2);
    std::array<double, MAX_RACHA_PUENTE_ITO + 4> psi{};
    detalle_ito::calcular_irf_wold(m, dt_cap, psi.data());

    std::array<double, MAX_RACHA_PUENTE_ITO + 4> kappa{};
    kappa.fill(1.0);
    if (escalas_fase.size() >= dt_cap) {
        for (std::size_t i = 0; i < dt_cap; ++i) {
            const double ki = escalas_fase[i];
            if (std::isfinite(ki) && ki > 1e-6) kappa[i] = ki;
        }
    }

    double var_rel = 0.0;
    if (cfg.integracion == ark::metrologia::Integracion::NIVEL) {
        if (dt <= dt_cap) {
            for (std::size_t j = 0; j < dt; ++j) {
                double s_psi = 0.0;
                for (std::size_t i = j; i < dt; ++i) s_psi += kappa[i] * psi[i - j];
                var_rel += s_psi * s_psi;
            }
        } else {
            double s_psi = 0.0;
            double k_med = 0.0;
            for (std::size_t k = 0; k < dt_cap; ++k) { s_psi += psi[k]; k_med += kappa[k]; }
            k_med /= static_cast<double>(dt_cap);
            const double h_ef = (m.valido && m.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE)
                                    ? std::clamp(0.5 + m.d, 0.1, 0.95) : 0.5;
            var_rel = (k_med * s_psi) * (k_med * s_psi) * std::pow(static_cast<double>(dt), 2.0 * h_ef);
        }
    } else {
        const double k_last = kappa[dt_cap - 1];
        for (std::size_t j = 0; j < dt_cap; ++j) var_rel += (k_last * psi[j]) * (k_last * psi[j]);
    }
    const double z_max_ef = derivar_z_max_hueco(cfg);
    const double p_obj = std::clamp(potencia_objetivo, 0.5, 0.9999);
    const double z_beta = ark::metrologia::est::norm_ppf(p_obj);
    return (z_max_ef + z_beta) * std::sqrt(std::max(1e-15, var_rel));
}

// Evalúa la doble llave (r <= r_max && |Z_hueco| <= z_max_ef) entre (x_ini, x_fin) separados por r barras faltantes
// (delta_t = r + 1) y, si pasa, construye el Puente Condicionado de Itô / Hoffman-Ribak en M semillas RNG.
// Cuando `escalas_fase` tiene longitud >= r_hueco + 1 (PurgaDiurna activa), aplica el escalado diurno causal exacto
// kappa_s = sqrt(v_fase(s) / v_media) = 1 / EspejoPurgaNulo::factor(fase_s, dia_s) tanto a Z_hueco como al puente.
[[nodiscard]] inline ResultadoHuecoIto evaluar_y_reconstruir_hueco_ito(
    const ConfiguracionPuenteHistoricoIto& cfg,
    const ModeloLinealNulo& m,
    std::span<const double> hist_est_previa,
    double x_ini,
    double x_fin,
    std::size_t r_hueco,
    std::uint64_t id_hueco = 0,
    std::span<const double> escalas_fase = {}) noexcept {

    ResultadoHuecoIto out{};
    out.r_hueco = r_hueco;
    out.delta_t = r_hueco + 1;
    out.x_ini = x_ini;
    out.x_fin = x_fin;
    out.delta_x_obs = x_fin - x_ini;
    out.z_max_efectivo = derivar_z_max_hueco(cfg);

    if (r_hueco == 0 || !std::isfinite(x_ini) || !std::isfinite(x_fin)) {
        out.flags_emitidos = ark::dataset::FLAG_CORRUPCION_CAIDA;
        return out;
    }

    const std::size_t dt = r_hueco + 1;
    if (r_hueco > cfg.r_max || r_hueco > MAX_RACHA_PUENTE_ITO) {
        out.excedio_r_max = true;
        out.flags_emitidos = ark::dataset::FLAG_PAUSA_UMBRAL | ark::dataset::FLAG_CORRUPCION_CAIDA;
    }

    const double sigma2_e = (m.valido && m.sigma2_innov > 1e-15) ? m.sigma2_innov : 1.0;
    const double sd_e = std::sqrt(sigma2_e);

    // 1. Respuesta al impulso de Wold psi[0..dt-1] y factores diurnos causales kappa[0..dt-1]
    std::array<double, MAX_RACHA_PUENTE_ITO + 4> psi{};
    const std::size_t dt_cap = std::min(dt, MAX_RACHA_PUENTE_ITO + 2);
    detalle_ito::calcular_irf_wold(m, dt_cap, psi.data());

    std::array<double, MAX_RACHA_PUENTE_ITO + 4> kappa{};
    kappa.fill(1.0);
    const bool usar_escalado_diurno = (escalas_fase.size() >= dt_cap);
    if (usar_escalado_diurno) {
        for (std::size_t i = 0; i < dt_cap; ++i) {
            const double ki = escalas_fase[i];
            if (std::isfinite(ki) && ki > 1e-6) kappa[i] = ki;
        }
    }

    // 2. Pronóstico causal incondicional mu_z[0..dt_cap-1] en coordenadas estacionarias (purgadas si aplica)
    std::array<double, MAX_HIST_CAUSAL_PUENTE> buf_hist{};
    const std::size_t cap_prev = MAX_HIST_CAUSAL_PUENTE - dt_cap;
    const std::size_t n_copy = std::min(hist_est_previa.size(), cap_prev);
    const std::size_t off_prev = hist_est_previa.size() - n_copy;
    for (std::size_t i = 0; i < n_copy; ++i) buf_hist[i] = hist_est_previa[off_prev + i];

    std::array<double, MAX_RACHA_PUENTE_ITO + 4> mu_z{};
    for (std::size_t s = 0; s < dt_cap; ++s) {
        const double zp = detalle_ito::predecir_estacionario_un_paso(m, buf_hist.data(), n_copy + s);
        mu_z[s] = zp;
        buf_hist[n_copy + s] = zp;
    }

    // 3. Pesos de condicionamiento w[0..dt_cap-1] y varianza del borde v según NIVEL vs FLUJO (exactos con kappa_s)
    std::array<double, MAX_RACHA_PUENTE_ITO + 4> w{};
    double v_borde = 0.0;

    if (cfg.integracion == ark::metrologia::Integracion::NIVEL) {
        // En NIVEL: dX_s = kappa_s * y_s, con restricción sum_{s=0}^{dt-1} kappa_s * y_s = x_fin - x_ini.
        if (dt <= dt_cap) {
            std::array<double, MAX_RACHA_PUENTE_ITO + 4> c_col{};
            for (std::size_t j = 0; j < dt; ++j) {
                double s_psi = 0.0;
                for (std::size_t i = j; i < dt; ++i) s_psi += kappa[i] * psi[i - j];
                c_col[j] = s_psi;
                v_borde += sigma2_e * s_psi * s_psi;
            }
            for (std::size_t i = 0; i < dt; ++i) {
                double wi_purg = 0.0;
                for (std::size_t j = 0; j <= i; ++j) {
                    wi_purg += sigma2_e * psi[i - j] * c_col[j];
                }
                w[i] = kappa[i] * wi_purg;
            }
            double sum_mu = 0.0;
            for (std::size_t s = 0; s < dt; ++s) sum_mu += kappa[s] * mu_z[s];
            out.delta_x_esp = sum_mu;
        } else {
            double s_psi = 0.0;
            double k_med = 0.0;
            for (std::size_t k = 0; k < dt_cap; ++k) { s_psi += psi[k]; k_med += kappa[k]; }
            k_med /= static_cast<double>(dt_cap);
            const double h_ef = (m.valido && m.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE)
                                    ? std::clamp(0.5 + m.d, 0.1, 0.95) : 0.5;
            v_borde = sigma2_e * (k_med * s_psi) * (k_med * s_psi) * std::pow(static_cast<double>(dt), 2.0 * h_ef);
            out.delta_x_esp = static_cast<double>(dt) * k_med * (m.valido ? m.mu : 0.0);
        }
        out.sigma_dt = std::sqrt(std::max(1e-15, v_borde));
        out.z_hueco = (out.delta_x_obs - out.delta_x_esp) / out.sigma_dt;
    } else {
        // En FLUJO: x_s = kappa_s * y_s, y el borde observado en s = dt - 1 es x_{dt-1} = x_fin.
        const std::size_t idx_last = dt_cap - 1;
        const double k_last = kappa[idx_last];
        for (std::size_t j = 0; j <= idx_last; ++j) {
            const double term = k_last * psi[j];
            v_borde += sigma2_e * term * term;
        }
        for (std::size_t i = 0; i <= idx_last; ++i) {
            double wi = 0.0;
            for (std::size_t j = 0; j <= i; ++j) {
                wi += sigma2_e * psi[i - j] * psi[idx_last - j];
            }
            w[i] = kappa[i] * k_last * wi;
        }
        out.delta_x_esp = k_last * mu_z[idx_last] - x_ini;
        out.sigma_dt = std::sqrt(std::max(1e-15, v_borde));
        out.z_hueco = (x_fin - k_last * mu_z[idx_last]) / out.sigma_dt;
    }

    const double z_beta_80 = ark::metrologia::est::norm_ppf(0.80);
    out.delta_min_detectable_sigma_e = (out.z_max_efectivo + z_beta_80) * (out.sigma_dt / std::max(1e-15, sd_e));

    if (std::abs(out.z_hueco) > out.z_max_efectivo) {
        out.choque_en_hueco = true;
        out.flags_emitidos |= ark::dataset::FLAG_CHOQUE_EN_HUECO | ark::dataset::FLAG_CORRUPCION_CAIDA;
    }

    // Doble llave: si excedió r_max o hubo choque en el hueco, NO se reconstruye
    if (out.excedio_r_max || out.choque_en_hueco) {
        return out;
    }

    out.reconstruido = true;
    out.flags_emitidos = ark::dataset::FLAG_PUENTE_ITO;

    // 4. Trayectoria media condicional exacta (Kalman / Schur smoother entre ambos bordes)
    if (cfg.usar_recta_lineal_pura) {
        const double paso_lin = (x_fin - x_ini) / static_cast<double>(dt);
        for (std::size_t k = 0; k < r_hueco; ++k) {
            const double xk = x_ini + static_cast<double>(k + 1) * paso_lin;
            out.x_media_cond[k] = xk;
            out.y_est_media_cond[k] = (cfg.integracion == ark::metrologia::Integracion::NIVEL) ? paso_lin : xk;
        }
    } else if (cfg.integracion == ark::metrologia::Integracion::NIVEL) {
        const double innov_borde = (out.delta_x_obs - out.delta_x_esp) / std::max(1e-15, v_borde);
        double niv_acum = x_ini;
        for (std::size_t k = 0; k < r_hueco; ++k) {
            const double dy_k = kappa[k] * mu_z[k] + innov_borde * w[k];
            niv_acum += dy_k;
            out.y_est_media_cond[k] = dy_k;
            out.x_media_cond[k] = niv_acum;
        }
    } else {
        const double k_last = kappa[dt - 1];
        const double innov_borde = (x_fin - k_last * mu_z[dt - 1]) / std::max(1e-15, v_borde);
        for (std::size_t k = 0; k < r_hueco; ++k) {
            const double xk = kappa[k] * mu_z[k] + innov_borde * w[k];
            out.x_media_cond[k] = xk;
            out.y_est_media_cond[k] = xk;
        }
    }

    // 5. M realizaciones estocásticas por proyección condicional de Hoffman-Ribak (Itô / ARFIMA bridge)
    const std::size_t n_sem = std::min(cfg.n_semillas_rng, MAX_SEMILLAS_PUENTE_ITO);
    out.n_semillas_efectivas = n_sem;
    double max_err_anc = 0.0;

    for (std::size_t m_idx = 0; m_idx < n_sem; ++m_idx) {
        if (cfg.usar_recta_lineal_pura) {
            for (std::size_t k = 0; k < r_hueco; ++k) {
                out.x_rng[m_idx][k] = out.x_media_cond[k];
                out.y_est_rng[m_idx][k] = out.y_est_media_cond[k];
            }
            continue;
        }

        std::mt19937_64 rng(detalle_ito::splitmix64_ito(cfg.semilla_base ^ (id_hueco * 0x9E3779B97F4A7C15ULL) ^ (m_idx + 1)));
        std::normal_distribution<double> nd(0.0, sd_e);

        std::array<double, MAX_RACHA_PUENTE_ITO + 4> eps_sim{};
        std::array<double, MAX_RACHA_PUENTE_ITO + 4> z_tilde{};
        for (std::size_t s = 0; s < dt; ++s) eps_sim[s] = nd(rng);
        for (std::size_t i = 0; i < dt; ++i) {
            double acc = mu_z[i];
            for (std::size_t j = 0; j <= i; ++j) {
                acc += psi[i - j] * eps_sim[j];
            }
            z_tilde[i] = acc;
        }

        if (cfg.integracion == ark::metrologia::Integracion::NIVEL) {
            double sum_tilde = 0.0;
            for (std::size_t s = 0; s < dt; ++s) sum_tilde += kappa[s] * z_tilde[s];
            const double corr_factor = (out.delta_x_obs - sum_tilde) / std::max(1e-15, v_borde);

            double niv = x_ini;
            double sum_dy_check = 0.0;
            for (std::size_t s = 0; s < dt; ++s) {
                const double dy_s = kappa[s] * z_tilde[s] + corr_factor * w[s];
                sum_dy_check += dy_s;
                if (s < r_hueco) {
                    niv += dy_s;
                    out.y_est_rng[m_idx][s] = dy_s;
                    out.x_rng[m_idx][s] = niv;
                }
            }
            const double err_anc = std::abs(sum_dy_check - out.delta_x_obs);
            if (err_anc > max_err_anc) max_err_anc = err_anc;
        } else {
            const double x_tilde_last = kappa[dt - 1] * z_tilde[dt - 1];
            const double corr_factor = (x_fin - x_tilde_last) / std::max(1e-15, v_borde);
            for (std::size_t k = 0; k < r_hueco; ++k) {
                const double xk = kappa[k] * z_tilde[k] + corr_factor * w[k];
                out.x_rng[m_idx][k] = xk;
                out.y_est_rng[m_idx][k] = xk;
            }
            const double x_end_check = x_tilde_last + corr_factor * w[dt - 1];
            const double err_anc = std::abs(x_end_check - x_fin);
            if (err_anc > max_err_anc) max_err_anc = err_anc;
        }
    }
    out.max_error_anclaje_rng = max_err_anc;

    // 6. Dispersión inter-semilla RNG en el interior del hueco
    if (n_sem > 1 && r_hueco > 0) {
        double suma_sd = 0.0;
        for (std::size_t k = 0; k < r_hueco; ++k) {
            double media_k = 0.0;
            for (std::size_t m_idx = 0; m_idx < n_sem; ++m_idx) media_k += out.y_est_rng[m_idx][k];
            media_k /= static_cast<double>(n_sem);
            double var_k = 0.0;
            for (std::size_t m_idx = 0; m_idx < n_sem; ++m_idx) {
                const double d = out.y_est_rng[m_idx][k] - media_k;
                var_k += d * d;
            }
            suma_sd += std::sqrt(var_k / static_cast<double>(n_sem - 1));
        }
        out.sd_rng_media_trayectoria = suma_sd / static_cast<double>(r_hueco);
    }

    return out;
}

// Combina M estimaciones (S^(m), ee_S^(m)) mediante las reglas de imputación múltiple de Rubin (1987).
// Cero reservas dinámicas.
[[nodiscard]] inline ResumenRubinIto combinar_rubin(
    std::span<const double> s_imputaciones,
    std::span<const double> ee_imputaciones) noexcept {
    ResumenRubinIto out{};
    const std::size_t M = std::min(s_imputaciones.size(), ee_imputaciones.size());
    if (M == 0) return out;

    double sum_s = 0.0;
    double sum_w = 0.0;
    std::size_t validas = 0;
    for (std::size_t m = 0; m < M; ++m) {
        if (!std::isfinite(s_imputaciones[m]) || !std::isfinite(ee_imputaciones[m]) || ee_imputaciones[m] < 0.0) {
            continue;
        }
        sum_s += s_imputaciones[m];
        sum_w += ee_imputaciones[m] * ee_imputaciones[m];
        ++validas;
    }
    if (validas == 0) return out;

    out.m_imputaciones = validas;
    const double dm = static_cast<double>(validas);
    out.s_barra = sum_s / dm;
    out.w_barra = sum_w / dm;

    if (validas == 1) {
        out.b_entre = 0.0;
        out.t_total = out.w_barra;
        out.ee_rubin = std::sqrt(std::max(0.0, out.t_total));
        out.nu_barnard_rubin = std::numeric_limits<double>::infinity();
        out.gamma_info_perdida = 0.0;
        out.estadistico_t = (out.ee_rubin > 0.0) ? (out.s_barra / out.ee_rubin) : 0.0;
        out.p_bilateral = detalle_ito::p_bilateral_student_t(std::abs(out.estadistico_t), out.nu_barnard_rubin);
        out.valido = std::isfinite(out.p_bilateral);
        return out;
    }

    double ss_entre = 0.0;
    for (std::size_t m = 0; m < M; ++m) {
        if (!std::isfinite(s_imputaciones[m]) || !std::isfinite(ee_imputaciones[m]) || ee_imputaciones[m] < 0.0) {
            continue;
        }
        const double d = s_imputaciones[m] - out.s_barra;
        ss_entre += d * d;
    }
    out.b_entre = ss_entre / static_cast<double>(validas - 1);
    const double factor_entre = (1.0 + 1.0 / dm) * out.b_entre;
    out.t_total = out.w_barra + factor_entre;
    out.ee_rubin = std::sqrt(std::max(0.0, out.t_total));

    if (out.t_total <= 0.0) {
        out.gamma_info_perdida = 0.0;
        out.nu_barnard_rubin = std::numeric_limits<double>::infinity();
        out.estadistico_t = 0.0;
        out.p_bilateral = 1.0;
        out.valido = true;
        return out;
    }

    out.gamma_info_perdida = factor_entre / out.t_total;
    if (factor_entre <= 1e-18) {
        out.nu_barnard_rubin = std::numeric_limits<double>::infinity();
    } else {
        const double ratio = 1.0 + out.w_barra / factor_entre;
        out.nu_barnard_rubin = static_cast<double>(validas - 1) * ratio * ratio;
    }
    out.estadistico_t = out.s_barra / out.ee_rubin;
    out.p_bilateral = detalle_ito::p_bilateral_student_t(std::abs(out.estadistico_t), out.nu_barnard_rubin);
    out.valido = std::isfinite(out.p_bilateral);
    return out;
}

// Selector seguro para la entrada de Z0 (ObservadorSaltos) en NIVEL vs FLUJO:
// - Si el hueco no fue reconstruido, devuelve NaN (corta limpio en Z0).
// - En NIVEL, inyecta la trayectoria imputada x_rng[m][k] preservando continuidad sin teletransporte (+Delta).
// - En FLUJO, mantiene Z0 causal usando la media estacionaria calibrada mu_calma (nunca 0.0 literal),
//   evitando desplazar t_cambio hacia el hueco y evitando escalones falsos cuando mu_calma != 0.
[[nodiscard]] inline double valor_z0_imputado(
    const ResultadoHuecoIto& res,
    ark::metrologia::Integracion integracion,
    std::size_t m_idx,
    std::size_t k_en_hueco,
    double mu_calma = 0.0) noexcept {
    if (!res.reconstruido || m_idx >= res.n_semillas_efectivas || k_en_hueco >= res.r_hueco) {
        return ark::metrologia::NaN;
    }
    if (integracion == ark::metrologia::Integracion::NIVEL) {
        return res.x_rng[m_idx][k_en_hueco];
    }
    const double reposo_flujo = std::isfinite(mu_calma) ? mu_calma : 0.0;
    return reposo_flujo;
}

} // namespace ark::lilyth::piloto
