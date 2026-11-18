// ============================================================================
// test_fase2_pausa_y_puente_ito.cpp
// Verificación de la Fase 2:
//   - Máquina de estados en AdaptadorInputCausalLilyth (FLUJO_ACTIVO,
//     ESPERA_ACOTADA, EN_PAUSA_DURA, doble llave y evitar_teletransporte_nivel)
//   - Equivalencia exacta de la media condicional con la recta del Comandante
//     (100 -> 110, 120, 130, 140 -> 150) bajo paseo Browniano y conservación
//     exacta inicio-fin en las M realizaciones estocásticas de Itô
//   - Doble llave (r <= r_max && |Z_hueco| <= z_max) bloqueando choques ocultos
//     y cortes largos
//   - Puente condicionado de 2 bordes en memoria larga ARFIMA (Z1b)
// ============================================================================

#include "banco_sintetico_comun.hpp"
#include "protocolo_inputs_lilyth.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {

using namespace ark::dataset;
using namespace ark::lilyth;
using namespace ark::lilyth::piloto;
using namespace ark::metrologia;

int g_fallas = 0;

void verificar(const char* codigo, const char* desc, bool cond) {
    std::printf("[%s] %s: %s\n", cond ? "PASA" : "FALLA", codigo, desc);
    if (!cond) ++g_fallas;
}

} // namespace

int main() {
    std::printf("====================================================================================================\n");
    std::printf(" FASE 2: VERIFICACION DE MAQUINA DE PAUSA, DOBLE LLAVE Y PUENTE HISTORICO DE ITO\n");
    std::printf("====================================================================================================\n");

    // ------------------------------------------------------------------------
    // CF2.1: Máquina de estados en vivo (ESPERA_ACOTADA vs EN_PAUSA_DURA + evitar teletransporte en NIVEL)
    // ------------------------------------------------------------------------
    {
        AdaptadorInputCausalLilyth ad{};
        ad.integracion = Integracion::NIVEL;
        ad.racha_max_corrupcion = 3;
        ad.z_max_hueco = 3.5;
        ad.sigma_base = 1.0;
        ad.tau_int_base = 1.0;
        ad.h_hurst_base = 0.5;
        ad.evitar_teletransporte_nivel = true;

        const auto p0 = ad.preparar(100.0, ClasePasoLilyth::VALIDO);
        // Micro-corte de 2 barras (r = 2 <= 3): debe estar en ESPERA_ACOTADA y NO emitir valido=true con 100.0
        const auto p1 = ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
        const auto p2 = ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
        // Reconexión suave a 101.5 (dt = 3, sigma_dt = sqrt(3) = 1.732, Z = 1.5 / 1.732 = 0.866 <= 3.5)
        const auto p3 = ad.preparar(101.5, ClasePasoLilyth::VALIDO);

        const bool micro_ok = p0.valido_para_operador &&
                              (p1.estado == EstadoContinuidadLilyth::ESPERA_ACOTADA) && !p1.valido_para_operador &&
                              (p2.estado == EstadoContinuidadLilyth::ESPERA_ACOTADA) && !p2.valido_para_operador &&
                              p3.reconexion_tras_hueco && (p3.delta_t_hueco == 3) &&
                              !p3.choque_en_hueco_detectado && !p3.pausa_umbral_activada &&
                              p3.cortar_segmento_al_reenganchar;

        // Micro-corte de 2 barras con salto oculto gigante (101.5 -> 151.5, +50 en dt=3 => Z = 28.87 > 3.5)
        (void)ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
        (void)ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
        const auto p_chq = ad.preparar(151.5, ClasePasoLilyth::VALIDO);
        const bool chq_ok = p_chq.reconexion_tras_hueco && p_chq.choque_en_hueco_detectado &&
                            ((p_chq.flags_emitidos & FLAG_CHOQUE_EN_HUECO) != 0u);

        // Corte largo de 4 barras (r = 4 > r_max = 3): en la 4ta barra entra en EN_PAUSA_DURA
        (void)ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
        (void)ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
        const auto p_esp3 = ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
        const auto p_pau4 = ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
        const auto p_rec = ad.preparar(152.0, ClasePasoLilyth::VALIDO);

        const bool pausa_ok = (p_esp3.estado == EstadoContinuidadLilyth::ESPERA_ACOTADA) &&
                              (p_pau4.estado == EstadoContinuidadLilyth::EN_PAUSA_DURA) &&
                              p_pau4.pausa_umbral_activada &&
                              ((p_pau4.flags_emitidos & FLAG_PAUSA_UMBRAL) != 0u) &&
                              p_rec.reconexion_tras_hueco && p_rec.pausa_umbral_activada;

        verificar("CF2.1", "Maquina de estados en vivo (ESPERA_ACOTADA <= r_max, EN_PAUSA_DURA > r_max, deteccion Z_hueco y sin teletransporte en NIVEL)",
                  micro_ok && chq_ok && pausa_ok);
    }

    // ------------------------------------------------------------------------
    // CF2.2: Ejemplo exacto del Comandante (100 -> 4 barras perdidas -> 150 => 110, 120, 130, 140 + ruido de Itô)
    // ------------------------------------------------------------------------
    {
        ModeloLinealNulo m_bm{};
        m_bm.valido = true;
        m_bm.modo_usado = ModoNuloLineal::AR_LEVINSON_DURBIN;
        m_bm.p = 0;
        m_bm.d = 0.0;
        m_bm.mu = 0.0;
        m_bm.sigma2_innov = 100.0; // sd_innov = 10.0 para que +50 en dt=5 sea Z = 50 / (10*sqrt(5)) = 2.236 <= 3.5

        ConfiguracionPuenteHistoricoIto cfg{};
        cfg.integracion = Integracion::NIVEL;
        cfg.r_max = 4;
        cfg.z_max_hueco = 3.5;
        cfg.n_semillas_rng = 16;

        const auto res = evaluar_y_reconstruir_hueco_ito(cfg, m_bm, {}, 100.0, 150.0, 4, 1);

        const bool media_es_recta = res.reconstruido &&
                                    (std::abs(res.x_media_cond[0] - 110.0) < 1e-12) &&
                                    (std::abs(res.x_media_cond[1] - 120.0) < 1e-12) &&
                                    (std::abs(res.x_media_cond[2] - 130.0) < 1e-12) &&
                                    (std::abs(res.x_media_cond[3] - 140.0) < 1e-12);

        const bool anclaje_y_ruido_ok = (res.max_error_anclaje_rng < 1e-12) &&
                                        (res.sd_rng_media_trayectoria > 1.0) &&
                                        (res.flags_emitidos == FLAG_PUENTE_ITO);

        std::printf("  [Detalle CF2.2] 100 -> (r=4) -> 150: media_cond = [%.2f, %.2f, %.2f, %.2f] | max_err_anclaje_RNG = %.2e | sd_RNG = %.4f\n",
                    res.x_media_cond[0], res.x_media_cond[1], res.x_media_cond[2], res.x_media_cond[3],
                    res.max_error_anclaje_rng, res.sd_rng_media_trayectoria);

        verificar("CF2.2", "En 100 -> (r=4) -> 150, la esperanza del puente es [110, 120, 130, 140] y las 16 trayectorias de Ito conservan inicio-fin exacto (<1e-12)",
                  media_es_recta && anclaje_y_ruido_ok);
    }

    // ------------------------------------------------------------------------
    // CF2.3: Doble llave en evaluar_y_reconstruir_hueco_ito (bloqueo por |Z_hueco| > z_max o r > r_max)
    // ------------------------------------------------------------------------
    {
        ModeloLinealNulo m_calma{};
        m_calma.valido = true;
        m_calma.modo_usado = ModoNuloLineal::AR_LEVINSON_DURBIN;
        m_calma.p = 0;
        m_calma.d = 0.0;
        m_calma.mu = 0.0;
        m_calma.sigma2_innov = 1.0; // sd_innov = 1.0

        ConfiguracionPuenteHistoricoIto cfg{};
        cfg.integracion = Integracion::NIVEL;
        cfg.r_max = 3;
        cfg.z_max_hueco = 3.5;

        // Caso A: r = 2 <= 3, pero salto 100 -> 150 (+50 con sigma_dt = sqrt(3) => Z = 28.87 > 3.5) -> BLOQUEADO por llave de amplitud
        const auto res_salto = evaluar_y_reconstruir_hueco_ito(cfg, m_calma, {}, 100.0, 150.0, 2, 10);
        // Caso B: r = 5 > 3, movimiento pequeño 100 -> 101 -> BLOQUEADO por llave temporal
        const auto res_largo = evaluar_y_reconstruir_hueco_ito(cfg, m_calma, {}, 100.0, 101.0, 5, 11);

        const bool doble_llave_ok = !res_salto.reconstruido && res_salto.choque_en_hueco &&
                                    ((res_salto.flags_emitidos & FLAG_CHOQUE_EN_HUECO) != 0u) &&
                                    !res_largo.reconstruido && res_largo.excedio_r_max &&
                                    ((res_largo.flags_emitidos & FLAG_PAUSA_UMBRAL) != 0u);

        verificar("CF2.3", "Doble llave bloquea reconstruccion cuando |Z_hueco| > z_max (FLAG_CHOQUE_EN_HUECO) o cuando r > r_max (FLAG_PAUSA_UMBRAL)",
                  doble_llave_ok);
    }

    // ------------------------------------------------------------------------
    // CF2.4: Puente condicionado de 2 bordes en memoria larga ARFIMA (M2, Z1b: dt=2, r=3)
    // ------------------------------------------------------------------------
    {
        ParametrosBancoComun pb{};
        const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, pb, 202);
        const std::size_t m_freq = std::max<std::size_t>(32, static_cast<std::size_t>(std::pow(m2.x_cal.size(), 0.65)));
        const auto est_d = elw::estimar_calibracion(m2.x_cal.data(), m2.x_cal.size(), m_freq, -0.45, 0.45);

        using ObsDelta = ObservadorSaltos<24, 64, 120>;
        using ObsLil = ObservadorLilyth<ObsDelta, 8, 64>;
        ConfiguracionObservador cfg_d{};
        cfg_d.integracion = Integracion::FLUJO;
        cfg_d.agregacion = AgregacionEscala::SUMA;
        cfg_d.arl0_pasos_nativos = 450.0;
        cfg_d.alfa = ALFA_PROYECTO;
        cfg_d.medir_tau50 = false;
        cfg_d.escalas_max = 3;
        ObsDelta obs_d{};
        obs_d.calibrar(cfg_d, m2.x_cal);

        ConfiguracionLilyth cfg_l{};
        cfg_l.alfa = ALFA_PROYECTO;
        cfg_l.representacion = ModoRepresentacionLilyth::CRUDA;
        cfg_l.b_max_nulo = 16;
        cfg_l.hilos_nulo = 4;
        cfg_l.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
        cfg_l.d_whittle_cruda = est_d.valido ? est_d.d : 0.25;
        ObsLil lil{};
        (void)lil.calibrar(obs_d, cfg_l, cfg_d, m2.x_cal);

        const auto& mod = lil.ficha_nulo_cruda().modelo_lineal;
        const double sd_e = std::sqrt(std::max(1e-12, mod.sigma2_innov));

        ConfiguracionPuenteHistoricoIto cfg_f{};
        cfg_f.integracion = Integracion::FLUJO;
        cfg_f.r_max = 3;
        cfg_f.z_max_hueco = 4.5;
        cfg_f.n_semillas_rng = 16;

        const std::size_t w_prev = 256;
        const std::size_t l_perf = 16;
        const std::size_t dt_antes = 2;
        const std::size_t r_hueco = 3;

        double sum_err2_sufijo = 0.0, sum_err2_unilat = 0.0, sum_err2_ito_media = 0.0, sum_err2_ito_rng = 0.0;
        std::size_t n_tot = 0;
        double max_anc_arfima = 0.0;

        for (std::size_t idx_ev = 0; idx_ev < m2.tiempos_choque.size(); ++idx_ev) {
            const std::size_t tc = m2.tiempos_choque[idx_ev];
            if (tc < w_prev + 10 || tc + l_perf >= m2.x_op_limpia.size()) continue;

            std::vector<double> prev_limpia(w_prev);
            for (std::size_t i = 0; i < w_prev; ++i) prev_limpia[i] = m2.x_op_limpia[tc - w_prev + i];
            const double* perf = &m2.x_op_limpia[tc];

            std::vector<double> e_ref(l_perf, 0.0);
            innovaciones(mod, prev_limpia.data(), w_prev, perf, l_perf, e_ref.data());

            const std::size_t idx_fin = w_prev - dt_antes; // hueco en [idx_ini .. idx_fin]
            const std::size_t idx_ini = idx_fin + 1 - r_hueco;

            // 1) CORTE_SUFIJO
            std::vector<double> prev_suf(prev_limpia.begin() + static_cast<std::ptrdiff_t>(idx_fin + 1), prev_limpia.end());
            std::vector<double> e_suf(l_perf, 0.0);
            innovaciones(mod, prev_suf.data(), prev_suf.size(), perf, l_perf, e_suf.data());

            // 2) PUENTE_CAUSAL unilateral (1 borde)
            std::vector<double> prev_uni = prev_limpia;
            for (std::size_t k = 0; k < r_hueco; ++k) {
                prev_uni[idx_ini + k] = detalle_ito::predecir_estacionario_un_paso(mod, prev_uni.data(), idx_ini + k);
            }
            std::vector<double> e_uni(l_perf, 0.0);
            innovaciones(mod, prev_uni.data(), prev_uni.size(), perf, l_perf, e_uni.data());

            // 3) PUENTE_ITO condicionado a 2 bordes (usando idx_ini - 1 como x_ini e idx_fin + 1 como x_fin, ambos < tc)
            const std::span<const double> hist_antes_hueco(prev_limpia.data(), idx_ini);
            const double x_ini = prev_limpia[idx_ini - 1];
            const double x_fin = prev_limpia[idx_fin + 1];
            const auto p_ito = evaluar_y_reconstruir_hueco_ito(cfg_f, mod, hist_antes_hueco, x_ini, x_fin, r_hueco, idx_ev + 1);
            if (p_ito.max_error_anclaje_rng > max_anc_arfima) max_anc_arfima = p_ito.max_error_anclaje_rng;

            std::vector<double> prev_ito_med = prev_limpia;
            for (std::size_t k = 0; k < r_hueco; ++k) prev_ito_med[idx_ini + k] = p_ito.x_media_cond[k];
            std::vector<double> e_ito_med(l_perf, 0.0);
            innovaciones(mod, prev_ito_med.data(), prev_ito_med.size(), perf, l_perf, e_ito_med.data());

            // Promedio de innovaciones sobre las M=16 semillas RNG del puente de Itô
            std::vector<double> e_ito_rng_prom(l_perf, 0.0);
            for (std::size_t m_i = 0; m_i < p_ito.n_semillas_efectivas; ++m_i) {
                std::vector<double> prev_m = prev_limpia;
                for (std::size_t k = 0; k < r_hueco; ++k) prev_m[idx_ini + k] = p_ito.x_rng[m_i][k];
                std::vector<double> e_m(l_perf, 0.0);
                innovaciones(mod, prev_m.data(), prev_m.size(), perf, l_perf, e_m.data());
                for (std::size_t j = 0; j < l_perf; ++j) e_ito_rng_prom[j] += e_m[j] / static_cast<double>(p_ito.n_semillas_efectivas);
            }

            for (std::size_t j = 0; j < l_perf; ++j) {
                const double ds = (e_suf[j] - e_ref[j]) / sd_e;
                const double du = (e_uni[j] - e_ref[j]) / sd_e;
                const double dm = (e_ito_med[j] - e_ref[j]) / sd_e;
                const double dr = (e_ito_rng_prom[j] - e_ref[j]) / sd_e;
                sum_err2_sufijo += ds * ds;
                sum_err2_unilat += du * du;
                sum_err2_ito_media += dm * dm;
                sum_err2_ito_rng += dr * dr;
                ++n_tot;
            }
        }

        const double rmse_suf = std::sqrt(sum_err2_sufijo / static_cast<double>(n_tot));
        const double rmse_uni = std::sqrt(sum_err2_unilat / static_cast<double>(n_tot));
        const double rmse_ito_med = std::sqrt(sum_err2_ito_media / static_cast<double>(n_tot));
        const double rmse_ito_rng = std::sqrt(sum_err2_ito_rng / static_cast<double>(n_tot));

        std::printf("  [Detalle CF2.4 ARFIMA Z1b (dt=2, r=3)] RMSE/sd: CORTE_SUFIJO=%.4f | UNILATERAL=%.4f | ITO_2BORDES_MEDIA=%.4f | ITO_2BORDES_RNG16=%.4f\n",
                    rmse_suf, rmse_uni, rmse_ito_med, rmse_ito_rng);

        const bool arfima_ok = (rmse_ito_med <= 0.15) && (rmse_ito_med < rmse_suf) &&
                               (rmse_ito_med <= rmse_uni) && (max_anc_arfima < 1e-12);

        verificar("CF2.4", "En ARFIMA Z1b (dt=2, r=3), el puente de 2 bordes supera a CORTE_SUFIJO y a PUENTE_CAUSAL unilateral (RMSE <= 0.15, ESTRICTO)",
                  arfima_ok);
    }

    std::printf("----------------------------------------------------------------------------------------------------\n");
    std::printf(" FALLAS FASE 2: %d\n", g_fallas);
    std::printf("====================================================================================================\n");
    return g_fallas == 0 ? 0 : 1;
}
