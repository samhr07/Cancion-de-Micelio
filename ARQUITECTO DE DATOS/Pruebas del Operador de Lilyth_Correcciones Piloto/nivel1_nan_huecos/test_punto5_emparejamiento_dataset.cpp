// ============================================================================
// test_punto5_emparejamiento_dataset.cpp
// Corrección del Punto 5 (Auditoría Opus):
//   1. Conexión directa con ark::dataset::BarraCausal y ObservadorSaltosFilas
//      para evaluar ILIQUIDEZ_TRANSITORIA (Mercado Quieto) en M1, M2 y M3.
//   2. Emparejamiento estricto por clave física (t_apertura, escala_origen, signo)
//      frente al emparejamiento defectuoso por índice de vector (idx).
//   3. Evaluación bajo los regímenes sellados ESTRICTO y DULCE
//      (EXPECTATIVAS_COTA_NIVEL1.md, sello 54aa457e...).
// ============================================================================

#include "protocolo_inputs_lilyth.hpp"
#include "banco_sintetico_comun.hpp"

#include <ark/dataset/tipos_dataset.hpp>
#include <ark/metrologia/observador_saltos_filas.hpp>
#include <ark/metrologia/whittle_exacto.hpp>
#include <ark/lilyth/observador_lilyth.hpp>
#include <ark/lilyth/nulo_lilyth.hpp>
#include <ark/lilyth/modelo_lineal_nulo.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <tuple>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;
using namespace ark::lilyth::piloto;

constexpr std::size_t T_MAX_P5 = 64;
using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, T_MAX_P5>;
using ClaveEvento = std::tuple<std::int64_t, int, int>; // (t_apertura, escala_origen, signo)

struct EventoP5 {
    std::int64_t t_apertura = -1;
    int escala = 0;
    int signo = 0;
    DesenlaceRanura desenlace = DesenlaceRanura::INACTIVA;
    double tau_rest = NaN;
    std::vector<double> e_iny;
};

struct CorridaP5 {
    std::map<ClaveEvento, EventoP5> por_clave;
    std::vector<EventoP5> por_indice;
    std::size_t k_medibles = 0;       // RESTITUIDA || CENSURADA_HORIZONTE
    std::size_t n_cens_invalido = 0;
    std::size_t n_no_ancla = 0;
    double s_onsager = NaN;
    double ee_s_onsager = NaN;
    double p_onsager = NaN;
    double sd_innov = 1.0;
};

// Ejecuta el pipeline sobre un flujo de BarraCausal de ark/dataset/
CorridaP5 correr_flujo_barras_dataset(
    const SerieSinteticaCoordinada& mundo,
    Integracion integ,
    bool usar_purga,
    ModoNuloLineal modo_lin,
    const std::vector<ark::dataset::BarraCausal>& barras_op,
    const std::vector<double>& valores_op,
    bool usar_protocolo_nuevo) {

    ConfiguracionObservador cfg_d{};
    cfg_d.integracion = integ;
    cfg_d.agregacion = AgregacionEscala::SUMA;
    cfg_d.arl0_pasos_nativos = 450.0;
    cfg_d.alfa = ALFA_PROYECTO;
    cfg_d.purga = usar_purga;
    cfg_d.bins_fase = usar_purga ? 24 : 1;
    cfg_d.dias_minimos_purga = 20;
    cfg_d.medir_tau50 = false;
    cfg_d.escalas_max = 3;

    ObsDelta obs_d{};
    obs_d.calibrar(cfg_d, mundo.x_cal, {}, {},
                   usar_purga ? std::span<const std::size_t>(mundo.fase_cal) : std::span<const std::size_t>{},
                   usar_purga ? std::span<const std::int64_t>(mundo.dia_cal) : std::span<const std::int64_t>{});

    ConfiguracionLilyth cfg_l{};
    cfg_l.alfa = ALFA_PROYECTO;
    cfg_l.b_max_nulo = 16;
    cfg_l.hilos_nulo = 4;
    cfg_l.semilla_nulo = 20261117ULL;
    cfg_l.modo_nulo_lineal = modo_lin;
    cfg_l.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cfg_l.b_mc_resumen = 1999;
    if (modo_lin == ModoNuloLineal::ARFIMA_WHITTLE) {
        const std::size_t m_freq = std::max<std::size_t>(32, static_cast<std::size_t>(std::pow(mundo.x_cal.size(), 0.65)));
        const auto est_d = elw::estimar_calibracion(mundo.x_cal.data(), mundo.x_cal.size(), m_freq, -0.45, 0.45);
        cfg_l.d_whittle_cruda = est_d.valido ? est_d.d : 0.25;
    }

    ObsLil lil{};
    lil.calibrar(obs_d, cfg_l, cfg_d, mundo.x_cal, {},
                 usar_purga ? std::span<const std::size_t>(mundo.fase_cal) : std::span<const std::size_t>{},
                 usar_purga ? std::span<const std::int64_t>(mundo.dia_cal) : std::span<const std::int64_t>{});

    AdaptadorInputCausalLilyth adapt{};
    adapt.integracion = integ;
    adapt.admitir_mercado_quieto = usar_protocolo_nuevo;
    adapt.racha_max_corrupcion = 0;

    obs_d.reiniciar_flujo();
    for (std::size_t t = 0; t < barras_op.size(); ++t) {
        const ClasePasoLilyth clase = clasificar_barra_dataset(barras_op[t], valores_op[t]);
        double x_in = valores_op[t];
        bool val_in = true;
        if (usar_protocolo_nuevo) {
            const auto prep = adapt.preparar(valores_op[t], clase, 0.0);
            x_in = prep.x_para_operador;
            val_in = prep.valido_para_operador;
        } else {
            // Comportamiento actual en observador_saltos_filas.hpp (valido = !b.imputada)
            val_in = !barras_op[t].es_imputado && barras_op[t].hay_dato && std::isfinite(valores_op[t]);
            x_in = val_in ? valores_op[t] : NaN;
        }
        const std::size_t fa = usar_purga ? mundo.fase_op[t] : 0;
        const std::int64_t di = usar_purga ? mundo.dia_op[t] : 0;
        const auto sd = obs_d.paso(x_in, NaN, val_in, fa, di);
        (void)lil.paso(sd, x_in, val_in, NaN, fa, di);
    }

    lil.completar_nulos();

    CorridaP5 out{};
    const auto& fn = usar_purga ? lil.ficha_nulo_purgada() : lil.ficha_nulo_cruda();
    const auto& mod = fn.modelo_lineal;
    out.sd_innov = std::sqrt(std::max(1e-15, mod.sigma2_innov));
    out.s_onsager = fn.contraste_onsager.estadistico_s;
    out.ee_s_onsager = fn.contraste_onsager.ee_nulo_s;
    out.p_onsager = fn.contraste_onsager.p_bilateral;

    for (const auto& r : lil.cerradas_historico()) {
        const auto& rama = usar_purga ? r.purgada : r.cruda;
        if (rama.desenlace == DesenlaceRanura::RESTITUIDA || rama.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) {
            ++out.k_medibles;
        } else if (rama.desenlace == DesenlaceRanura::CENSURADA_DATO_INVALIDO) {
            ++out.n_cens_invalido;
        } else if (rama.desenlace == DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE) {
            ++out.n_no_ancla;
        }

        EventoP5 ev{};
        ev.t_apertura = r.t_apertura;
        ev.escala = r.escala_origen;
        ev.signo = r.signo;
        ev.desenlace = rama.desenlace;
        ev.tau_rest = rama.tau_rest_consolidado;

        const std::size_t n_perf = usar_purga ? r.n_perfil_purgada : r.n_perfil_cruda;
        const std::size_t n_prev = usar_purga ? r.n_hist_prev_purgada : r.n_hist_prev_cruda;
        const auto& perf = usar_purga ? r.perfil_purgada_est : r.perfil_cruda_est;
        const auto& prev = usar_purga ? r.hist_prev_purgada_est : r.hist_prev_cruda_est;

        if (n_perf > 0 && mod.valido) {
            ev.e_iny.assign(n_perf, 0.0);
            ark::lilyth::innovaciones(
                mod,
                n_prev > 0 ? prev.data() : nullptr,
                n_prev,
                perf.data(),
                n_perf,
                ev.e_iny.data());
        }
        out.por_indice.push_back(ev);
        out.por_clave[{r.t_apertura, r.escala_origen, r.signo}] = std::move(ev);
    }
    return out;
}

struct MetricasEmparejamientoP5 {
    std::size_t k_base = 0;
    std::size_t k_emparejados_clave = 0;
    double delta_censura = 0.0;
    double rmse_por_clave = 0.0;     // Correcto: emparejado por (t_apertura, escala, signo)
    double rmse_por_indice = 0.0;    // Defectuoso (1a corrida): emparejado por posición i en el vector
    double sesgo_tau_clave = 0.0;
    double ds_sobre_ee = 0.0;
    bool pasa_estricto = false;
    bool pasa_dulce = false;
};

MetricasEmparejamientoP5 evaluar_p5(const CorridaP5& base, const CorridaP5& eval) {
    MetricasEmparejamientoP5 m{};
    m.k_base = base.k_medibles;

    // 1. Emparejamiento CORRECTO por clave (t_apertura, escala_origen, signo)
    double sum_e2_clave = 0.0;
    std::size_t n_e_clave = 0;
    double sum_dtau = 0.0;
    std::size_t n_tau = 0;

    for (const auto& [cl, ev0] : base.por_clave) {
        if (ev0.desenlace != DesenlaceRanura::RESTITUIDA && ev0.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
        auto it = eval.por_clave.find(cl);
        if (it == eval.por_clave.end()) continue;
        const auto& evc = it->second;
        if (evc.desenlace == DesenlaceRanura::RESTITUIDA || evc.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) {
            ++m.k_emparejados_clave;
        }
        const std::size_t len = std::min(ev0.e_iny.size(), evc.e_iny.size());
        for (std::size_t k = 0; k < len; ++k) {
            const double d = (evc.e_iny[k] - ev0.e_iny[k]) / base.sd_innov;
            sum_e2_clave += d * d;
            ++n_e_clave;
        }
        if (std::isfinite(ev0.tau_rest) && std::isfinite(evc.tau_rest)) {
            sum_dtau += (evc.tau_rest - ev0.tau_rest);
            ++n_tau;
        }
    }

    // 2. Emparejamiento DEFECTUOSO por índice de vector (para demostrar el error del punto 5 de Opus)
    double sum_e2_idx = 0.0;
    std::size_t n_e_idx = 0;
    const std::size_t n_vec = std::min(base.por_indice.size(), eval.por_indice.size());
    for (std::size_t i = 0; i < n_vec; ++i) {
        const auto& ev0 = base.por_indice[i];
        const auto& evc = eval.por_indice[i];
        const std::size_t len = std::min(ev0.e_iny.size(), evc.e_iny.size());
        for (std::size_t k = 0; k < len; ++k) {
            const double d = (evc.e_iny[k] - ev0.e_iny[k]) / base.sd_innov;
            sum_e2_idx += d * d;
            ++n_e_idx;
        }
    }

    m.delta_censura = (m.k_base > 0)
                          ? (1.0 - static_cast<double>(m.k_emparejados_clave) / static_cast<double>(m.k_base))
                          : 0.0;
    m.rmse_por_clave = (n_e_clave > 0) ? std::sqrt(sum_e2_clave / static_cast<double>(n_e_clave)) : 0.0;
    m.rmse_por_indice = (n_e_idx > 0) ? std::sqrt(sum_e2_idx / static_cast<double>(n_e_idx)) : 0.0;
    m.sesgo_tau_clave = (n_tau > 0) ? (sum_dtau / static_cast<double>(n_tau)) : 0.0;

    const double ee = (std::isfinite(base.ee_s_onsager) && base.ee_s_onsager > 1e-12) ? base.ee_s_onsager : 0.05;
    m.ds_sobre_ee = (std::isfinite(eval.s_onsager) && std::isfinite(base.s_onsager))
                        ? (std::abs(eval.s_onsager - base.s_onsager) / ee)
                        : 99.0;

    m.pasa_estricto = (m.delta_censura <= 0.10) && (m.rmse_por_clave <= 0.15) &&
                      (std::abs(m.sesgo_tau_clave) <= 0.25) && (m.ds_sobre_ee <= 0.25);
    m.pasa_dulce    = (m.delta_censura <  0.50) && (m.rmse_por_clave <  0.50) &&
                      (std::abs(m.sesgo_tau_clave) <= 1.50) && (m.ds_sobre_ee <= 1.00);
    return m;
}

} // namespace

int main() {
    ParametrosBancoComun params{};
    const auto m1 = generar_mundo_m1_ar1_flujo(0.80, params, 101);
    const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, params, 202);
    const auto m3 = generar_mundo_m3_nivel_purga(0.65, params, 303);

    const std::array<double, 6> tasas_q = {0.0, 0.005, 0.01, 0.02, 0.05, 0.10};

    std::printf("========================================================================================================================\n");
    std::printf(" PUNTO 5: CONEXION CON ark::dataset::BarraCausal Y EMPAREJAMIENTO POR CLAVE (t_apertura, escala, signo) VS POR INDICE\n");
    std::printf("========================================================================================================================\n");

    // ------------------------------------------------------------------------
    // PARTE 5A: Neutralidad física de las barras imputadas por mercado quieto
    // La serie física ya tiene barras quietas (dx=0 en NIVEL, x=0 en FLUJO);
    // la única diferencia entre ACTUAL y PROTOCOLO es la etiqueta b.es_imputado
    // que ObservadorSaltosFilas hoy convierte en valido = !b.es_imputado.
    // ------------------------------------------------------------------------
    std::printf("\n[PARTE 5A] NEUTRALIDAD FISICA SOBRE LA MISMA TRAYECTORIA CON BARRAS DE MERCADO QUIETO (b.es_imputado = true):\n");
    std::printf("           Compara ACTUAL (valido = !b.es_imputado) vs PROTOCOLO (MERCADO_QUIETO -> valido = true)\n");
    std::printf("           contra la verdad física de esa misma serie, emparejando por (t_apertura, escala, signo).\n\n");

    auto correr_tabla_5a = [&](const char* titulo, const SerieSinteticaCoordinada& m, Integracion integ, bool purga, ModoNuloLineal mlin) {
        std::printf("--- %s ---\n", titulo);
        std::printf("  q(%%) | Modo      | K_base | K_emp | CensInv | NoAncla | dCens | RMSE_clave | RMSE_idx(err) | |dS|/ee | ESTRICTO | DULCE\n");
        std::printf("-------+-----------+--------+-------+---------+---------+-------+------------+---------------+---------+----------+-------\n");

        // Precalcular mascara de pulsos para no tocar los choques
        std::vector<std::uint8_t> en_pulso(m.x_op_limpia.size(), 0u);
        for (std::size_t t0 : m.tiempos_choque) {
            for (std::size_t k = 0; k < 6 && (t0 + k) < en_pulso.size(); ++k) en_pulso[t0 + k] = 1u;
        }

        for (double q : tasas_q) {
            std::mt19937_64 rng(detalle_nulo::splitmix64(20261117ULL ^ static_cast<std::uint64_t>(q * 100000.0)));
            std::uniform_real_distribution<double> u01(0.0, 1.0);

            std::vector<double> x_fisica = m.x_op_limpia;
            std::vector<ark::dataset::BarraCausal> barras_limpias(x_fisica.size());
            std::vector<ark::dataset::BarraCausal> barras_etiquetadas(x_fisica.size());

            double acum_nivel = 100.0;
            for (std::size_t t = 0; t < x_fisica.size(); ++t) {
                const bool quieta = (q > 0.0 && t >= 1 && !en_pulso[t] && u01(rng) < q);
                if (integ == Integracion::NIVEL) {
                    const double dx = (t == 0) ? (m.x_op_limpia[0] - 100.0) : (m.x_op_limpia[t] - m.x_op_limpia[t - 1]);
                    acum_nivel += quieta ? 0.0 : dx;
                    x_fisica[t] = acum_nivel;
                } else {
                    if (quieta) x_fisica[t] = 0.0;
                }

                barras_limpias[t].hay_dato = true;
                barras_limpias[t].es_imputado = false;
                barras_limpias[t].flags_control = 0u;
                barras_limpias[t].tipo_hueco = ark::dataset::TipoHueco::NINGUNO;

                barras_etiquetadas[t] = barras_limpias[t];
                if (quieta) {
                    barras_etiquetadas[t].hay_dato = false;
                    barras_etiquetadas[t].es_imputado = true; // Etiqueta de rejilla en mercado quieto
                }
            }

            const auto verdad_q = correr_flujo_barras_dataset(m, integ, purga, mlin, barras_limpias, x_fisica, true);
            const auto r_act    = correr_flujo_barras_dataset(m, integ, purga, mlin, barras_etiquetadas, x_fisica, false);
            const auto r_pro    = correr_flujo_barras_dataset(m, integ, purga, mlin, barras_etiquetadas, x_fisica, true);

            const auto m_act = evaluar_p5(verdad_q, r_act);
            const auto m_pro = evaluar_p5(verdad_q, r_pro);

            std::printf(" %5.1f | ACTUAL    | %6zu | %5zu | %7zu | %7zu | %5.2f | %10.4f | %13.4f | %7.4f | %-8s | %-5s\n",
                        q * 100.0, m_act.k_base, m_act.k_emparejados_clave, r_act.n_cens_invalido, r_act.n_no_ancla,
                        m_act.delta_censura, m_act.rmse_por_clave, m_act.rmse_por_indice, m_act.ds_sobre_ee,
                        m_act.pasa_estricto ? "PASA" : "FALLA", m_act.pasa_dulce ? "PASA" : "FALLA");
            std::printf(" %5.1f | PROTOCOLO | %6zu | %5zu | %7zu | %7zu | %5.2f | %10.4f | %13.4f | %7.4f | %-8s | %-5s\n",
                        q * 100.0, m_pro.k_base, m_pro.k_emparejados_clave, r_pro.n_cens_invalido, r_pro.n_no_ancla,
                        m_pro.delta_censura, m_pro.rmse_por_clave, m_pro.rmse_por_indice, m_pro.ds_sobre_ee,
                        m_pro.pasa_estricto ? "PASA" : "FALLA", m_pro.pasa_dulce ? "PASA" : "FALLA");
        }
        std::printf("\n");
    };

    correr_tabla_5a("Mundo M1 (FLUJO AR(1) phi=0.80)", m1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN);
    correr_tabla_5a("Mundo M2 (FLUJO ARFIMA(1, d=0.25, 0))", m2, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE);
    correr_tabla_5a("Mundo M3 (NIVEL AR(1) + PurgaDiurna 24 fases)", m3, Integracion::NIVEL, true, ModoNuloLineal::AR_LEVINSON_DURBIN);

    // ------------------------------------------------------------------------
    // PARTE 5B: Sensibilidad contra la serie activa sin quietud (q = 0)
    // Muestra cómo el emparejamiento por clave (t_apertura, escala, signo)
    // corrige el RMSE espurio del emparejamiento por índice cuando se compara
    // la serie con q% de quietud contra la serie activa q = 0%.
    // ------------------------------------------------------------------------
    std::printf("[PARTE 5B] COMPARACION CONTRA LA SERIE ACTIVA BASE (q = 0%%) EN M2 (ARFIMA_WHITTLE):\n");
    std::printf("           Demostracion directa de por que el emparejamiento por indice (RMSE_idx) inflaba el error a ~2.4-5.5\n");
    std::printf("           cuando cambiaba el conjunto de ranuras, frente al emparejamiento real por (t_apertura, escala, signo).\n\n");
    std::printf("  q(%%) | Modo      | K_q0 | K_emp_q0 | dCens | RMSE_clave(real) | RMSE_idx(espurio) | Sesgo(tau) | |dS|/ee | ESTRICTO | DULCE\n");
    std::printf("-------+-----------+------+----------+-------+------------------+-------------------+------------+---------+----------+-------\n");

    std::vector<ark::dataset::BarraCausal> barras_q0(m2.x_op_limpia.size());
    for (auto& b : barras_q0) { b.hay_dato = true; b.es_imputado = false; }
    const auto base_q0_m2 = correr_flujo_barras_dataset(m2, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE, barras_q0, m2.x_op_limpia, true);

    std::vector<std::uint8_t> en_pulso_m2(m2.x_op_limpia.size(), 0u);
    for (std::size_t t0 : m2.tiempos_choque) {
        for (std::size_t k = 0; k < 6 && (t0 + k) < en_pulso_m2.size(); ++k) en_pulso_m2[t0 + k] = 1u;
    }

    for (double q : tasas_q) {
        std::mt19937_64 rng(detalle_nulo::splitmix64(20261117ULL ^ static_cast<std::uint64_t>(q * 100000.0)));
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        std::vector<double> x_mod = m2.x_op_limpia;
        std::vector<ark::dataset::BarraCausal> barras_mod(x_mod.size());
        for (std::size_t t = 0; t < x_mod.size(); ++t) {
            const bool quieta = (q > 0.0 && t >= 1 && !en_pulso_m2[t] && u01(rng) < q);
            barras_mod[t].hay_dato = !quieta;
            barras_mod[t].es_imputado = quieta;
            if (quieta) x_mod[t] = 0.0;
        }
        const auto r_act = correr_flujo_barras_dataset(m2, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE, barras_mod, x_mod, false);
        const auto r_pro = correr_flujo_barras_dataset(m2, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE, barras_mod, x_mod, true);

        const auto m_act = evaluar_p5(base_q0_m2, r_act);
        const auto m_pro = evaluar_p5(base_q0_m2, r_pro);

        std::printf(" %5.1f | ACTUAL    | %4zu | %8zu | %5.2f | %16.4f | %17.4f | %10.4f | %7.4f | %-8s | %-5s\n",
                    q * 100.0, m_act.k_base, m_act.k_emparejados_clave, m_act.delta_censura,
                    m_act.rmse_por_clave, m_act.rmse_por_indice, m_act.sesgo_tau_clave, m_act.ds_sobre_ee,
                    m_act.pasa_estricto ? "PASA" : "FALLA", m_act.pasa_dulce ? "PASA" : "FALLA");
        std::printf(" %5.1f | PROTOCOLO | %4zu | %8zu | %5.2f | %16.4f | %17.4f | %10.4f | %7.4f | %-8s | %-5s\n",
                    q * 100.0, m_pro.k_base, m_pro.k_emparejados_clave, m_pro.delta_censura,
                    m_pro.rmse_por_clave, m_pro.rmse_por_indice, m_pro.sesgo_tau_clave, m_pro.ds_sobre_ee,
                    m_pro.pasa_estricto ? "PASA" : "FALLA", m_pro.pasa_dulce ? "PASA" : "FALLA");
    }

    return 0;
}
