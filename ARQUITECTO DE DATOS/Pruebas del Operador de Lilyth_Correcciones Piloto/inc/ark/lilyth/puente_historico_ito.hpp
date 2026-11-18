// ============================================================================
// puente_historico_ito.hpp -- Reconstrucción Histórica por Puente Estocástico
// de Itô / Hoffman-Ribak (AR / ARFIMA / Browniano) con Doble Llave y Etiquetas
// Fase 2 (Piloto Nivel 1, C++20). Cero reservas dinámicas en el núcleo de puente.
// ============================================================================
#pragma once

#include "ark/dataset/tipos_dataset.hpp"
#include "ark/lilyth/modelo_lineal_nulo.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"
#include "ark/metrologia/whittle_exacto.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

namespace ark::lilyth::piloto {

inline constexpr std::size_t MAX_RACHA_PUENTE_ITO = 16;
inline constexpr std::size_t MAX_SEMILLAS_PUENTE_ITO = 32;
inline constexpr std::size_t MAX_HIST_CAUSAL_PUENTE = 512;

struct ConfiguracionPuenteHistoricoIto {
    ark::metrologia::Integracion integracion = ark::metrologia::Integracion::NIVEL;
    std::size_t r_max = 3;                 // Llave temporal: solo reconstruye si 1 <= r <= r_max
    double z_max_hueco = 3.5;              // Llave de amplitud: solo reconstruye si |Z_hueco| <= z_max_hueco
    std::size_t n_semillas_rng = 16;       // Número de realizaciones RNG para medir dispersión inicio-fin
    bool usar_recta_lineal_pura = false;   // true = Idea 1 (recta determinista sin ruido); false = Puente de Itô/Schur
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

namespace detalle_ito {

[[nodiscard]] inline std::uint64_t splitmix64_ito(std::uint64_t x) noexcept {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
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
    if (!m.valido || n_buf == 0) return (m.valido && std::isfinite(m.mu)) ? m.mu : 0.0;
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

    ark::metrologia::elw::coeficientes(m.d, n_ef, pi_loc.data());
    for (std::size_t t = 0; t < n_sub; ++t) {
        double s = 0.0;
        for (std::size_t k = 0; k <= t && k < n_ef; ++k) s += pi_loc[k] * sub[t - k];
        u_loc[t] = s;
    }
    double u_pred = m.mu;
    for (std::size_t i = 1; i <= m.p && i <= n_sub; ++i) {
        u_pred += m.phi[i - 1] * (u_loc[n_sub - i] - m.mu);
    }
    double cola = 0.0;
    for (std::size_t k = 1; k <= n_sub && k < n_ef; ++k) {
        cola += pi_loc[k] * sub[n_sub - k];
    }
    return u_pred - cola;
}

} // namespace detalle_ito

// Evalúa la doble llave (r <= r_max && |Z_hueco| <= z_max) entre (x_ini, x_fin) separados por r barras faltantes
// (delta_t = r + 1) y, si pasa, construye el Puente Condicionado de Itô / Hoffman-Ribak en M semillas RNG.
[[nodiscard]] inline ResultadoHuecoIto evaluar_y_reconstruir_hueco_ito(
    const ConfiguracionPuenteHistoricoIto& cfg,
    const ModeloLinealNulo& m,
    std::span<const double> hist_est_previa,
    double x_ini,
    double x_fin,
    std::size_t r_hueco,
    std::uint64_t id_hueco = 0) noexcept {

    ResultadoHuecoIto out{};
    out.r_hueco = r_hueco;
    out.delta_t = r_hueco + 1;
    out.x_ini = x_ini;
    out.x_fin = x_fin;
    out.delta_x_obs = x_fin - x_ini;

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

    // 1. Respuesta al impulso de Wold psi[0..dt-1]
    std::array<double, MAX_RACHA_PUENTE_ITO + 4> psi{};
    const std::size_t dt_cap = std::min(dt, MAX_RACHA_PUENTE_ITO + 2);
    detalle_ito::calcular_irf_wold(m, dt_cap, psi.data());

    // 2. Pronóstico causal incondicional mu_z[0..dt_cap-1] desde hist_est_previa
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

    // 3. Pesos de condicionamiento w[0..dt_cap-1] y varianza del borde v según NIVEL vs FLUJO
    std::array<double, MAX_RACHA_PUENTE_ITO + 4> w{};
    double v_borde = 0.0;

    if (cfg.integracion == ark::metrologia::Integracion::NIVEL) {
        // En NIVEL: z_s = dX_s, y la restricción de borde derecho es sum_{s=0}^{dt-1} z_s = x_fin - x_ini.
        // Si r_hueco excede MAX_RACHA_PUENTE_ITO, escalamos la varianza asintótica para evaluar Z_hueco.
        if (dt <= dt_cap) {
            std::array<double, MAX_RACHA_PUENTE_ITO + 4> c_col{};
            for (std::size_t j = 0; j < dt; ++j) {
                double s_psi = 0.0;
                for (std::size_t i = j; i < dt; ++i) s_psi += psi[i - j];
                c_col[j] = s_psi;
                v_borde += sigma2_e * s_psi * s_psi;
            }
            for (std::size_t i = 0; i < dt; ++i) {
                double wi = 0.0;
                for (std::size_t j = 0; j <= i; ++j) {
                    wi += sigma2_e * psi[i - j] * c_col[j];
                }
                w[i] = wi;
            }
            double sum_mu = 0.0;
            for (std::size_t s = 0; s < dt; ++s) sum_mu += mu_z[s];
            out.delta_x_esp = sum_mu;
        } else {
            double s_psi = 0.0;
            for (std::size_t k = 0; k < dt_cap; ++k) s_psi += psi[k];
            const double h_ef = (m.valido && m.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE) ? std::clamp(0.5 + m.d, 0.1, 0.95) : 0.5;
            v_borde = sigma2_e * s_psi * s_psi * std::pow(static_cast<double>(dt), 2.0 * h_ef);
            out.delta_x_esp = static_cast<double>(dt) * (m.valido ? m.mu : 0.0);
        }
        out.sigma_dt = std::sqrt(std::max(1e-15, v_borde));
        out.z_hueco = (out.delta_x_obs - out.delta_x_esp) / out.sigma_dt;
    } else {
        // En FLUJO: z_s = x_s es estacionario, y el borde observado en s = dt - 1 es z_{dt-1} = x_fin.
        const std::size_t idx_last = dt_cap - 1;
        for (std::size_t j = 0; j <= idx_last; ++j) {
            v_borde += sigma2_e * psi[j] * psi[j];
        }
        for (std::size_t i = 0; i <= idx_last; ++i) {
            double wi = 0.0;
            for (std::size_t j = 0; j <= i; ++j) {
                wi += sigma2_e * psi[i - j] * psi[idx_last - j];
            }
            w[i] = wi;
        }
        out.delta_x_esp = mu_z[idx_last] - x_ini;
        out.sigma_dt = std::sqrt(std::max(1e-15, v_borde));
        out.z_hueco = (x_fin - mu_z[idx_last]) / out.sigma_dt;
    }

    if (std::abs(out.z_hueco) > cfg.z_max_hueco) {
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
            const double dy_k = mu_z[k] + innov_borde * w[k];
            niv_acum += dy_k;
            out.y_est_media_cond[k] = dy_k;
            out.x_media_cond[k] = niv_acum;
        }
    } else {
        const double innov_borde = (x_fin - mu_z[dt - 1]) / std::max(1e-15, v_borde);
        for (std::size_t k = 0; k < r_hueco; ++k) {
            const double xk = mu_z[k] + innov_borde * w[k];
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
            for (std::size_t s = 0; s < dt; ++s) sum_tilde += z_tilde[s];
            const double corr_factor = (out.delta_x_obs - sum_tilde) / std::max(1e-15, v_borde);

            double niv = x_ini;
            double sum_dy_check = 0.0;
            for (std::size_t s = 0; s < dt; ++s) {
                const double dy_s = z_tilde[s] + corr_factor * w[s];
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
            const double corr_factor = (x_fin - z_tilde[dt - 1]) / std::max(1e-15, v_borde);
            for (std::size_t k = 0; k < r_hueco; ++k) {
                const double xk = z_tilde[k] + corr_factor * w[k];
                out.x_rng[m_idx][k] = xk;
                out.y_est_rng[m_idx][k] = xk;
            }
            const double x_end_check = z_tilde[dt - 1] + corr_factor * w[dt - 1];
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

} // namespace ark::lilyth::piloto
