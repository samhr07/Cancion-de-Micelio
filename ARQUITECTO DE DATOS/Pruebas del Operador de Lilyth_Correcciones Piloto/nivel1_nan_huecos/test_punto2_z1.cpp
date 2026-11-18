// ============================================================================
// test_punto2_z1.cpp
// Corrección del Punto 2 (Auditoría Opus): Evaluación honesta de Zona Z1
// (hist_prev_est antes de t_cambio) usando EXCLUSIVAMENTE ark::lilyth::innovaciones()
// y ModeloLinealNulo::sigma2_innov reales sobre:
//   - Mundo M1 (FLUJO AR(1), phi=0.80) calibrado con AR_LEVINSON_DURBIN
//   - Mundo M2 (FLUJO ARFIMA(1, d=0.25, 0), phi1=0.40) calibrado con ARFIMA_WHITTLE
// Sin fórmulas sustitutas ni sabotajes artificiales, comparando ESTRICTO (<=0.15)
// y DULCE (<0.50) según EXPECTATIVAS_COTA_NIVEL1.md (sello 54aa457e...).
// ============================================================================

#include "banco_sintetico_comun.hpp"

#include <ark/lilyth/observador_lilyth.hpp>
#include <ark/lilyth/modelo_lineal_nulo.hpp>
#include <ark/metrologia/whittle_exacto.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;
using namespace ark::lilyth::piloto;

constexpr std::size_t T_MAX_P2 = 64;
using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, T_MAX_P2>;

enum class PoliticaZ1Honesta {
    COMPRESION_BASE = 0, // Lo que hace copiar_perfil_estacionario hoy: salta el inválido (continue)
    CORTE_SUFIJO = 1,    // Conserva solo el sufijo contiguo limpio inmediatamente anterior a t_cambio
    PUENTE_CAUSAL = 2    // Predicción causal a 1 paso (e_t = 0) usando el ModeloLinealNulo real (AR o ARFIMA)
};

// Predicción causal a 1 paso exacta (inversa algebraica de ark::lilyth::innovaciones con e_t = 0)
void aplicar_puente_causal_modelo_real(
    const ModeloLinealNulo& m,
    const std::vector<double>& prev0,
    const std::vector<std::uint8_t>& corrupto,
    std::vector<double>& prev_proc) {

    prev_proc.clear();
    prev_proc.reserve(prev0.size());

    if (m.modo_usado != ModoNuloLineal::ARFIMA_WHITTLE || m.d == 0.0) {
        for (std::size_t t = 0; t < prev0.size(); ++t) {
            if (!corrupto[t]) {
                prev_proc.push_back(prev0[t]);
            } else {
                double pred = m.mu;
                for (std::size_t i = 1; i <= m.p && i <= prev_proc.size(); ++i) {
                    pred += m.phi[i - 1] * (prev_proc[prev_proc.size() - i] - m.mu);
                }
                prev_proc.push_back(pred);
            }
        }
        return;
    }

    // Caso ARFIMA_WHITTLE real: u_t = sum_{k=0}^t pi_k z_{t-k}, con pi_0 = 1
    const std::size_t n = prev0.size();
    std::vector<double> pi(n, 0.0);
    elw::coeficientes(m.d, n, pi.data());
    std::vector<double> u_hist(n, 0.0);

    for (std::size_t t = 0; t < n; ++t) {
        double cola_frac = 0.0;
        for (std::size_t k = 1; k <= t; ++k) {
            cola_frac += pi[k] * prev_proc[t - k];
        }
        if (!corrupto[t]) {
            const double z_t = prev0[t];
            prev_proc.push_back(z_t);
            u_hist[t] = z_t + cola_frac;
        } else {
            double u_pred = m.mu;
            for (std::size_t i = 1; i <= m.p && i <= t; ++i) {
                u_pred += m.phi[i - 1] * (u_hist[t - i] - m.mu);
            }
            const double z_pred = u_pred - cola_frac;
            prev_proc.push_back(z_pred);
            u_hist[t] = u_pred;
        }
    }
}

struct MetricaCeldaZ1 {
    double rmse_rel = 0.0;
    double rmse_paso0_rel = 0.0;
    double max_abs_rel = 0.0;
    bool pasa_estricto = false; // <= 0.15
    bool pasa_dulce = false;    // <  0.50
};

MetricaCeldaZ1 medir_celda_z1(
    const ModeloLinealNulo& mod,
    const std::vector<std::vector<double>>& hist_prev_limpias,
    const std::vector<std::vector<double>>& perfiles_limpios,
    std::size_t delta_t_antes_cambio, // 1 = hueco termina en t_cambio - 1; 2 = termina en t_cambio - 2, etc.
    std::size_t racha_r,
    PoliticaZ1Honesta pol) {

    const double sd_innov = std::sqrt(std::max(1e-15, mod.sigma2_innov));
    double sum_err2 = 0.0;
    double sum_err0_2 = 0.0;
    double max_err = 0.0;
    std::size_t n_tot = 0;
    std::size_t n_ev_usados = 0;

    std::vector<double> e0;
    std::vector<double> ec;
    std::vector<double> prev_proc;

    for (std::size_t ev = 0; ev < perfiles_limpios.size(); ++ev) {
        const auto& prev0 = hist_prev_limpias[ev];
        const auto& perf0 = perfiles_limpios[ev];
        if (perf0.empty() || prev0.size() < delta_t_antes_cambio + racha_r + 4) continue;

        // 1. Innovaciones limpias con ark::lilyth::innovaciones() real
        e0.assign(perf0.size(), 0.0);
        ark::lilyth::innovaciones(mod, prev0.data(), prev0.size(), perf0.data(), perf0.size(), e0.data());

        // 2. Máscara de racha corrupta de longitud racha_r que termina en prev0.size() - delta_t_antes_cambio
        std::vector<std::uint8_t> corrupto(prev0.size(), 0u);
        const std::size_t idx_fin = prev0.size() - delta_t_antes_cambio; // inclusive
        const std::size_t idx_ini = idx_fin + 1 - racha_r;
        for (std::size_t k = idx_ini; k <= idx_fin; ++k) {
            corrupto[k] = 1u;
        }

        // 3. Construir prev_proc según la política honesta (sin ningún -prev0[i-1] falso)
        prev_proc.clear();
        if (pol == PoliticaZ1Honesta::COMPRESION_BASE) {
            for (std::size_t i = 0; i < prev0.size(); ++i) {
                if (!corrupto[i]) prev_proc.push_back(prev0[i]);
            }
        } else if (pol == PoliticaZ1Honesta::CORTE_SUFIJO) {
            std::size_t ini = 0;
            for (std::size_t i = 0; i < prev0.size(); ++i) {
                if (corrupto[i]) ini = i + 1;
            }
            for (std::size_t i = ini; i < prev0.size(); ++i) {
                prev_proc.push_back(prev0[i]);
            }
        } else {
            aplicar_puente_causal_modelo_real(mod, prev0, corrupto, prev_proc);
        }

        // 4. Innovaciones contaminadas con ark::lilyth::innovaciones() real
        ec.assign(perf0.size(), 0.0);
        ark::lilyth::innovaciones(
            mod,
            prev_proc.empty() ? nullptr : prev_proc.data(),
            prev_proc.size(),
            perf0.data(),
            perf0.size(),
            ec.data());

        for (std::size_t k = 0; k < perf0.size(); ++k) {
            const double d = (ec[k] - e0[k]) / sd_innov;
            sum_err2 += d * d;
            max_err = std::max(max_err, std::abs(d));
            ++n_tot;
            if (k == 0) sum_err0_2 += d * d;
        }
        ++n_ev_usados;
    }

    MetricaCeldaZ1 out{};
    out.rmse_rel = (n_tot > 0) ? std::sqrt(sum_err2 / static_cast<double>(n_tot)) : 0.0;
    out.rmse_paso0_rel = (n_ev_usados > 0) ? std::sqrt(sum_err0_2 / static_cast<double>(n_ev_usados)) : 0.0;
    out.max_abs_rel = max_err;
    out.pasa_estricto = (out.rmse_rel <= 0.15);
    out.pasa_dulce = (out.rmse_rel < 0.50);
    return out;
}

void evaluar_mundo_z1(
    const char* titulo_mundo,
    const SerieSinteticaCoordinada& mundo,
    ModoNuloLineal modo_lin) {

    ConfiguracionObservador cfg_d{};
    cfg_d.integracion = Integracion::FLUJO;
    cfg_d.agregacion = AgregacionEscala::SUMA;
    cfg_d.arl0_pasos_nativos = 450.0;
    cfg_d.alfa = ALFA_PROYECTO;
    cfg_d.medir_tau50 = false;
    cfg_d.escalas_max = 3;

    ObsDelta obs_d{};
    obs_d.calibrar(cfg_d, mundo.x_cal);

    ConfiguracionLilyth cfg_l{};
    cfg_l.alfa = ALFA_PROYECTO;
    cfg_l.b_max_nulo = 16;
    cfg_l.hilos_nulo = 4;
    cfg_l.modo_nulo_lineal = modo_lin;
    if (modo_lin == ModoNuloLineal::ARFIMA_WHITTLE) {
        const std::size_t m_freq = std::max<std::size_t>(32, static_cast<std::size_t>(std::pow(mundo.x_cal.size(), 0.65)));
        const auto est_d = elw::estimar_calibracion(mundo.x_cal.data(), mundo.x_cal.size(), m_freq, -0.45, 0.45);
        cfg_l.d_whittle_cruda = est_d.valido ? est_d.d : 0.25;
    }

    ObsLil lil{};
    const bool cal_ok = lil.calibrar(obs_d, cfg_l, cfg_d, mundo.x_cal);
    if (!cal_ok) {
        std::printf("ERROR: Fallo calibracion en %s\n", titulo_mundo);
        return;
    }

    obs_d.reiniciar_flujo();
    for (std::size_t t = 0; t < mundo.x_op_limpia.size(); ++t) {
        const auto sd = obs_d.paso(mundo.x_op_limpia[t], NaN, true, 0, 0);
        (void)lil.paso(sd, mundo.x_op_limpia[t], true, NaN, 0, 0);
    }

    std::vector<std::vector<double>> hist_prev;
    std::vector<std::vector<double>> perfiles;
    for (const auto& r : lil.cerradas_historico()) {
        if (r.n_perfil_cruda > 0 && r.n_hist_prev_cruda >= 32) {
            hist_prev.emplace_back(r.hist_prev_cruda_est.begin(), r.hist_prev_cruda_est.begin() + r.n_hist_prev_cruda);
            perfiles.emplace_back(r.perfil_cruda_est.begin(), r.perfil_cruda_est.begin() + r.n_perfil_cruda);
        }
    }

    const auto& mod = lil.ficha_nulo_cruda().modelo_lineal;
    const double sd_innov = std::sqrt(std::max(1e-15, mod.sigma2_innov));

    std::printf("\n====================================================================================================\n");
    std::printf(" %s\n", titulo_mundo);
    std::printf(" Modelo ajustado: modo=%d | d=%.4f | p=%zu | mu=%.4f | sd_serie=%.4f | sd_innov=%.4f | Eventos=%zu\n",
                static_cast<int>(mod.modo_usado), mod.d, mod.p, mod.mu, mod.sd, sd_innov, perfiles.size());
    std::printf("====================================================================================================\n");
    std::printf(" dt_antes | racha r | Politica Z1     | RMSE_total/sd | RMSE_paso0/sd | Max|de|/sd | ESTRICTO(<=0.15) | DULCE(<0.50)\n");
    std::printf("----------+---------+-----------------+---------------+---------------+------------+------------------+-------------\n");

    const std::array<std::size_t, 5> distancias = {1, 2, 3, 5, 10};
    const std::array<std::size_t, 4> rachas = {1, 2, 3, 5};

    for (std::size_t dt : distancias) {
        for (std::size_t r : rachas) {
            for (int p_idx = 0; p_idx < 3; ++p_idx) {
                const auto pol = static_cast<PoliticaZ1Honesta>(p_idx);
                const char* nom = (pol == PoliticaZ1Honesta::COMPRESION_BASE) ? "COMPRESION_BASE"
                                : (pol == PoliticaZ1Honesta::CORTE_SUFIJO)    ? "CORTE_SUFIJO"
                                                                              : "PUENTE_CAUSAL";
                const auto m = medir_celda_z1(mod, hist_prev, perfiles, dt, r, pol);
                std::printf(" %8zu | %7zu | %-15s | %13.4f | %13.4f | %10.4f | %-16s | %-11s\n",
                            dt, r, nom, m.rmse_rel, m.rmse_paso0_rel, m.max_abs_rel,
                            m.pasa_estricto ? "PASA" : "FALLA",
                            m.pasa_dulce    ? "PASA" : "FALLA");
            }
        }
        std::printf("----------+---------+-----------------+---------------+---------------+------------+------------------+-------------\n");
    }
}

} // namespace

int main() {
    ParametrosBancoComun params{};
    const auto m1 = generar_mundo_m1_ar1_flujo(0.80, params, 101);
    const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, params, 202);

    evaluar_mundo_z1("MUNDO M1: FLUJO AR(1) phi=0.80 calibrado con ModoNuloLineal::AR_LEVINSON_DURBIN",
                     m1, ModoNuloLineal::AR_LEVINSON_DURBIN);

    evaluar_mundo_z1("MUNDO M2: FLUJO ARFIMA(1, d=0.25, 0) phi1=0.40 calibrado con ModoNuloLineal::ARFIMA_WHITTLE",
                     m2, ModoNuloLineal::ARFIMA_WHITTLE);

    return 0;
}
