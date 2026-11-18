// ============================================================================
// test_fase3_integracion_operadores.cpp
// Verificación de la Fase 3:
//   - CF3.1: Tubería completa BarraCausal -> EntradaPasoCxx ->
//            ObservadorSaltosFilas -> ObservadorLilyth::paso_filas() en M3
//   - CF3.2: Veto automático de FLAG_PUENTE_ITO en Z2 ([t_cambio, t_a])
//            bloqueando la vía (b) aunque x sea finito
//   - CF3.3: Veto automático de FLAG_PUENTE_ITO en Z3 (t > t_a) impidiendo
//            que un valor sorteado por RNG dispare un cruce falso eta <= 0.5
//   - CF3.4: PoliticaZ1Piloto::PUENTE_ITO_HISTORICO integrada en
//            ObservadorLilyth::copiar_perfil_estacionario() sobre M2 (ARFIMA)
// ============================================================================

#include "ark/metrologia/observador_saltos_filas.hpp"
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
using namespace ark::ventanas;

constexpr std::size_t T_MAX_F3 = 64;
using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsDeltaFilas = ObservadorSaltosFilas<24, 64, 120, 1440, 8, 256>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, T_MAX_F3>;

int g_fallas = 0;

void verificar(const char* codigo, const char* desc, bool cond) {
    std::printf("[%s] %s: %s\n", cond ? "PASA" : "FALLA", codigo, desc);
    if (!cond) ++g_fallas;
}

} // namespace

int main() {
    std::printf("====================================================================================================\n");
    std::printf(" FASE 3: VERIFICACION DE INTEGRACION EN OBSERVADORSALTOSFILAS Y OBSERVADORLILYTH\n");
    std::printf("====================================================================================================\n");

    // ------------------------------------------------------------------------
    // CF3.1: Tubería completa BarraCausal -> EntradaPasoCxx -> ObservadorSaltosFilas -> paso_filas() en M3
    // ------------------------------------------------------------------------
    {
        ParametrosBancoComun pb{};
        pb.n_op = 6000;
        const auto m3 = generar_mundo_m3_nivel_purga(0.65, pb, 303);

        // Construir filas de calibración desde m3.x_cal
        std::vector<EntradaPasoCxx> filas_cal(m3.x_cal.size());
        for (std::size_t t = 0; t < m3.x_cal.size(); ++t) {
            BarraCausal bc{};
            bc.cierre_real_utc_ms = (t + 1) * 60000ULL;
            bc.estampa_disponible_ms = bc.cierre_real_utc_ms;
            bc.close = m3.x_cal[t];
            bc.volume = 10.0;
            bc.valor_nominal = bc.close * 10.0;
            bc.es_sesion_con_dato = true;
            bc.es_imputado = false;
            bc.hay_dato = true;
            bc.estado_supervivencia = EstadoSupervivencia::COTIZANDO;
            bc.tipo_hueco = TipoHueco::NINGUNO;
            bc.flags_control = 0;
            filas_cal[t] = a_entrada_paso_cxx(bc, bc.close);
        }

        // Construir serie operativa M3 con q = 2% de MERCADO_QUIETO fuera de los pulsos
        std::vector<std::uint8_t> en_pulso(m3.x_op_limpia.size(), 0u);
        for (std::size_t t0 : m3.tiempos_choque) {
            for (std::size_t k = 0; k < 6 && (t0 + k) < en_pulso.size(); ++k) en_pulso[t0 + k] = 1u;
        }
        std::mt19937_64 rng(detalle_nulo::splitmix64(20261118ULL));
        std::uniform_real_distribution<double> u01(0.0, 1.0);

        std::vector<EntradaPasoCxx> filas_op(m3.x_op_limpia.size());
        double niv = 100.0;
        for (std::size_t t = 0; t < m3.x_op_limpia.size(); ++t) {
            const bool quieta = (t >= 1 && !en_pulso[t] && u01(rng) < 0.02);
            const double dx = (t == 0) ? (m3.x_op_limpia[0] - 100.0) : (m3.x_op_limpia[t] - m3.x_op_limpia[t - 1]);
            niv += quieta ? 0.0 : dx;

            BarraCausal bo{};
            bo.cierre_real_utc_ms = (m3.x_cal.size() + t + 1) * 60000ULL;
            bo.estampa_disponible_ms = bo.cierre_real_utc_ms;
            bo.close = niv;
            bo.volume = quieta ? 0.0 : 10.0;
            bo.valor_nominal = bo.close * bo.volume;
            bo.es_sesion_con_dato = !quieta;
            bo.es_imputado = quieta;
            bo.hay_dato = !quieta;
            bo.estado_supervivencia = EstadoSupervivencia::COTIZANDO;
            bo.tipo_hueco = TipoHueco::NINGUNO;
            bo.flags_control = quieta ? FLAG_FILA_INSERTADA : 0u;
            filas_op[t] = a_entrada_paso_cxx(bo, bo.close);
        }

        ObsDeltaFilas::Fuente fuente{};
        fuente.base.reloj = TipoRelojVentana::PARED;
        fuente.base.congelar_en_huecos = false;
        fuente.base.metodo_imputacion = MetodoImputacion::FORWARD_FILL;
        fuente.base.umbral_imputacion_maximo = 1.0;
        fuente.delta_nu = 1.0;
        fuente.agregacion = AgregacionBarra::ULTIMO;

        auto correr_modo = [&](bool activar_protocolo) {
            ConfiguracionObservador cd{};
            cd.integracion = Integracion::NIVEL;
            cd.agregacion = AgregacionEscala::SUMA;
            cd.arl0_pasos_nativos = 450.0;
            cd.alfa = ALFA_PROYECTO;
            cd.purga = true;
            cd.bins_fase = 24;
            cd.dias_minimos_purga = 20;
            cd.medir_tau50 = false;
            cd.escalas_max = 3;
            cd.usar_protocolo_huecos = activar_protocolo;

            ObsDeltaFilas odf{};
            odf.calibrar_filas(cd, fuente, filas_cal, m3.fase_cal, m3.dia_cal);

            ConfiguracionLilyth cl{};
            cl.alfa = ALFA_PROYECTO;
            cl.representacion = ModoRepresentacionLilyth::PURGADA;
            cl.b_max_nulo = 16;
            cl.hilos_nulo = 4;
            cl.usar_protocolo_filas_piloto = activar_protocolo;

            ObsLil lil{};
            (void)lil.calibrar(odf.observador(), cl, cd, m3.x_cal, {}, m3.fase_cal, m3.dia_cal);

            for (std::size_t t = 0; t < filas_op.size(); ++t) {
                const auto sf = odf.paso_fila(filas_op[t], m3.fase_op[t], m3.dia_op[t]);
                (void)lil.paso_filas(sf, false, m3.fase_op[t], m3.dia_op[t]);
            }
            std::size_t k_med = 0, k_cens_inv = 0;
            for (const auto& r : lil.cerradas_historico()) {
                if (r.purgada.desenlace == DesenlaceRanura::RESTITUIDA || r.purgada.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) ++k_med;
                if (r.purgada.desenlace == DesenlaceRanura::CENSURADA_DATO_INVALIDO) ++k_cens_inv;
            }
            return std::make_pair(k_med, k_cens_inv);
        };

        const auto [k_med_base, k_cens_base] = correr_modo(false);
        const auto [k_med_prot, k_cens_prot] = correr_modo(true);

        std::printf("  [Detalle CF3.1 paso_filas M3 q=2%%] BASE: medibles=%zu, cens_inv=%zu | PROTOCOLO: medibles=%zu, cens_inv=%zu\n",
                    k_med_base, k_cens_base, k_med_prot, k_cens_prot);

        verificar("CF3.1", "En paso_filas() sobre M3 (q=2% mercado quieto), PROTOCOLO conserva las ranuras (cens_inv=0) frente al colapso de BASE",
                  (k_med_prot > 0) && (k_cens_prot == 0) && (k_med_base == 0));
    }

    // ------------------------------------------------------------------------
    // CF3.2 y CF3.3: Veto automático de FLAG_PUENTE_ITO en Z2 y en Z3 dentro de ObservadorLilyth
    // ------------------------------------------------------------------------
    {
        ParametrosBancoComun pb{};
        pb.n_op = 6000;
        const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, pb, 202);
        const std::size_t m_freq = std::max<std::size_t>(32, static_cast<std::size_t>(std::pow(m2.x_cal.size(), 0.65)));
        const auto est_d = elw::estimar_calibracion(m2.x_cal.data(), m2.x_cal.size(), m_freq, -0.45, 0.45);

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.agregacion = AgregacionEscala::SUMA;
        cd.arl0_pasos_nativos = 450.0;
        cd.alfa = ALFA_PROYECTO;
        cd.medir_tau50 = false;
        cd.escalas_max = 3;

        ObsDelta od{};
        od.calibrar(cd, m2.x_cal);
        od.reiniciar_flujo();
        std::vector<ObsDelta::Salida> sd_m2(m2.x_op_limpia.size());
        for (std::size_t t = 0; t < m2.x_op_limpia.size(); ++t) {
            sd_m2[t] = od.paso(m2.x_op_limpia[t], NaN, true, 0, 0);
        }

        ConfiguracionLilyth cl{};
        cl.alfa = ALFA_PROYECTO;
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 16;
        cl.hilos_nulo = 4;
        cl.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
        cl.d_whittle_cruda = est_d.valido ? est_d.d : 0.25;
        cl.forzar_via_b = true;
        cl.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cl.r_max_z2_piloto = 0;
        cl.politica_z3_piloto = ConfiguracionLilyth::PoliticaZ3Piloto::RELOJ_FISICO_TOPE_R;
        cl.r_max_z3_piloto = 2;

        // Corrida limpia de referencia para localizar eventos con t_cambio < t_apertura y tau_rest > 2
        ObsLil lil_ref{};
        (void)lil_ref.calibrar(od, cl, cd, m2.x_cal);
        for (std::size_t t = 0; t < m2.x_op_limpia.size(); ++t) {
            (void)lil_ref.paso(sd_m2[t], m2.x_op_limpia[t], true, NaN, 0, 0);
        }
        lil_ref.completar_nulos();

        std::int64_t tc_target_z2 = -1, ta_target_z2 = -1;
        std::int64_t ta_target_z3 = -1;
        double tau_rest_ref_z3 = NaN;
        double y_ancla_ref_z3 = 0.0;

        for (const auto& r : lil_ref.cerradas_historico()) {
            if (tc_target_z2 < 0 && r.t_cambio >= 0 && r.t_cambio < r.t_apertura && r.cruda.nulo_onsager.evaluado) {
                tc_target_z2 = r.t_cambio;
                ta_target_z2 = r.t_apertura;
            } else if (ta_target_z3 < 0 && r.cruda.desenlace == DesenlaceRanura::RESTITUIDA &&
                       std::isfinite(r.cruda.tau_rest_consolidado) && r.cruda.tau_rest_consolidado >= 4.0) {
                ta_target_z3 = r.t_apertura;
                tau_rest_ref_z3 = r.cruda.tau_rest_consolidado;
                y_ancla_ref_z3 = r.cruda.y_ancla;
            }
        }

        // Correr inyectando:
        //  (1) En tc_target_z2: un valor finito pero con FLAG_PUENTE_ITO -> debe bloquear vía (b) en esa ranura
        //  (2) En ta_target_z3 + 1 (tau = 1): un valor falsamente restituido (y_ancla_ref_z3) con valido=true pero con FLAG_PUENTE_ITO
        //      -> NO debe disparar cruce en tau=1, sino conservar tau_rest_ref_z3 exacto!
        ObsLil lil_test{};
        (void)lil_test.calibrar(od, cl, cd, m2.x_cal);
        for (std::size_t t = 0; t < m2.x_op_limpia.size(); ++t) {
            double x_in = m2.x_op_limpia[t];
            std::uint8_t fl = 0;
            if (static_cast<std::int64_t>(t) == tc_target_z2) {
                fl = FLAG_PUENTE_ITO;
            }
            if (static_cast<std::int64_t>(t) == ta_target_z3 + 1) {
                x_in = y_ancla_ref_z3; // eta = 0.0 <= 0.5 si no hubiera veto de FLAG_PUENTE_ITO!
                fl = FLAG_PUENTE_ITO;
            }
            (void)lil_test.paso(sd_m2[t], x_in, true, NaN, 0, 0, NaN, {}, fl);
        }
        lil_test.completar_nulos();

        bool veto_z2_ok = false;
        bool veto_z3_ok = false;
        for (const auto& r : lil_test.cerradas_historico()) {
            if (r.t_apertura == ta_target_z2) {
                veto_z2_ok = r.perfil_corrupto_z2_cruda && r.cruda.nulo_onsager.sin_nulo_b && !r.cruda.nulo_onsager.evaluado;
            }
            if (r.t_apertura == ta_target_z3) {
                veto_z3_ok = std::isnan(r.cruda.eta[0]) &&
                             std::isfinite(r.cruda.tau_rest_consolidado) &&
                             (std::abs(r.cruda.tau_rest_consolidado - tau_rest_ref_z3) < 1e-12);
            }
        }

        verificar("CF3.2", "Veto automatico de FLAG_PUENTE_ITO en Z2 bloquea via (b) aunque la barra tenga valor finito",
                  veto_z2_ok);
        verificar("CF3.3", "Veto automatico de FLAG_PUENTE_ITO en Z3 impide cruce falso eta <= 0.5 y preserva tau_rest real",
                  veto_z3_ok);

        // --------------------------------------------------------------------
        // CF3.4: PoliticaZ1Piloto::PUENTE_ITO_HISTORICO integrada en ObservadorLilyth sobre M2 (Z1b: dt=2, r=3)
        // --------------------------------------------------------------------
        std::vector<std::uint8_t> val_z1b(m2.x_op_limpia.size(), 1u);
        for (const auto& r : lil_ref.cerradas_historico()) {
            const std::int64_t idx_fin = r.t_cambio - 2;
            const std::int64_t idx_ini = idx_fin + 1 - 3;
            for (std::int64_t t = idx_ini; t <= idx_fin; ++t) {
                if (t >= 0 && static_cast<std::size_t>(t) < val_z1b.size()) val_z1b[static_cast<std::size_t>(t)] = 0u;
            }
        }

        ConfiguracionLilyth cl_ito = cl;
        cl_ito.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO;
        cl_ito.r_max_z1_piloto = 3;
        cl_ito.z_max_hueco_piloto = 4.5;

        ObsLil lil_ito{};
        (void)lil_ito.calibrar(od, cl_ito, cd, m2.x_cal);
        for (std::size_t t = 0; t < m2.x_op_limpia.size(); ++t) {
            const bool ok = (val_z1b[t] != 0u);
            (void)lil_ito.paso(sd_m2[t], ok ? m2.x_op_limpia[t] : NaN, ok, NaN, 0, 0);
        }
        lil_ito.completar_nulos();

        const auto& mod = lil_ref.ficha_nulo_cruda().modelo_lineal;
        const double sd_e = std::sqrt(std::max(1e-15, mod.sigma2_innov));
        double sum_err2 = 0.0;
        std::size_t n_e = 0;
        for (std::size_t i = 0; i < lil_ref.cerradas_historico().size() && i < lil_ito.cerradas_historico().size(); ++i) {
            const auto& r0 = lil_ref.cerradas_historico()[i];
            const auto& ri = lil_ito.cerradas_historico()[i];
            if (r0.t_apertura != ri.t_apertura || r0.n_perfil_cruda == 0) continue;
            std::vector<double> e0(r0.n_perfil_cruda, 0.0), ei(ri.n_perfil_cruda, 0.0);
            innovaciones(mod, r0.hist_prev_cruda_est.data(), r0.n_hist_prev_cruda, r0.perfil_cruda_est.data(), r0.n_perfil_cruda, e0.data());
            innovaciones(mod, ri.hist_prev_cruda_est.data(), ri.n_hist_prev_cruda, ri.perfil_cruda_est.data(), ri.n_perfil_cruda, ei.data());
            for (std::size_t k = 0; k < e0.size(); ++k) {
                const double d = (ei[k] - e0[k]) / sd_e;
                sum_err2 += d * d;
                ++n_e;
            }
        }
        const double rmse_ito_op = (n_e > 0) ? std::sqrt(sum_err2 / static_cast<double>(n_e)) : 99.0;
        const double s0 = lil_ref.ficha_nulo_cruda().contraste_onsager.estadistico_s;
        const double ee0 = lil_ref.ficha_nulo_cruda().contraste_onsager.ee_nulo_s;
        const double si = lil_ito.ficha_nulo_cruda().contraste_onsager.estadistico_s;
        const double ds_ee = std::abs(si - s0) / std::max(1e-12, ee0);

        std::printf("  [Detalle CF3.4 PUENTE_ITO_HISTORICO en ObservadorLilyth (M2 Z1b dt=2, r=3)] RMSE(e)/sd = %.4f | |dS|/ee_S = %.4f\n",
                    rmse_ito_op, ds_ee);

        verificar("CF3.4", "PoliticaZ1Piloto::PUENTE_ITO_HISTORICO en ObservadorLilyth pasa ESTRICTO en M2 Z1b (RMSE <= 0.15 y |dS|/ee_S <= 0.25)",
                  (rmse_ito_op <= 0.15) && (ds_ee <= 0.25));
    }

    std::printf("----------------------------------------------------------------------------------------------------\n");
    std::printf(" FALLAS FASE 3: %d\n", g_fallas);
    std::printf("====================================================================================================\n");
    return g_fallas == 0 ? 0 : 1;
}
