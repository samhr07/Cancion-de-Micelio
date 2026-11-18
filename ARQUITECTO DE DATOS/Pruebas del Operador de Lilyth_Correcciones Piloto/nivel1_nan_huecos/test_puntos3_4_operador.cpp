// ============================================================================
// test_puntos3_4_operador.cpp
// Corrección de los Puntos 3, 4 y 5 (Auditoría Opus):
// Evaluación de Z1, Z2 y Z3 ejecutando EXCLUSIVAMENTE el operador real:
//   - ObservadorLilyth::paso() -> copiar_perfil_estacionario() -> avanzar()
//   - ObservadorLilyth::completar_nulos() con forzar_via_b = true
//   - detalle_nulo::ejecutar_via_b_evento()
//   - Emparejamiento estricto de eventos por clave (t_apertura, escala_origen, signo)
//   - Comparación simultánea contra régimen ESTRICTO y régimen DULCE
//     (EXPECTATIVAS_COTA_NIVEL1.md, sello 54aa457e...).
// ============================================================================

#include "banco_sintetico_comun.hpp"

#include <ark/lilyth/observador_lilyth.hpp>
#include <ark/lilyth/nulo_lilyth.hpp>
#include <ark/lilyth/modelo_lineal_nulo.hpp>
#include <ark/metrologia/whittle_exacto.hpp>

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

constexpr std::size_t T_MAX_OP = 64;
using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, T_MAX_OP>;
using ClaveEvento = std::tuple<std::int64_t, int, int>; // (t_apertura, escala_origen, signo)

struct EventoCapturado {
    std::int64_t t_apertura = -1;
    std::int64_t t_cambio = -1;
    int escala_origen = 0;
    int signo = 0;
    DesenlaceRanura desenlace = DesenlaceRanura::INACTIVA;
    double tau_rest = NaN;
    double min_eta_2 = NaN;
    double p_det_ini = NaN;
    double u_pit_ons = NaN;
    double delta_tau_ons = NaN;
    std::vector<double> e_iny;
};

struct ResumenEjecucionReal {
    std::map<ClaveEvento, EventoCapturado> eventos;
    std::size_t n_evaluados = 0;
    double s_onsager = NaN;
    double ee_s_onsager = NaN;
    double p_onsager = NaN;
    double p_holm_onsager = NaN;
    double sd_innov = 1.0;
};

// Corre ObservadorLilyth real sobre una serie de salidas delta precalculadas y una serie (x_op, val_op)
ResumenEjecucionReal correr_lilyth_real(
    const ObsDelta& obs_d_cal,
    const ConfiguracionObservador& cfg_d,
    const ConfiguracionLilyth& cfg_l,
    const std::vector<double>& x_cal,
    const std::vector<ObsDelta::Salida>& salidas_delta,
    const std::vector<double>& x_op,
    const std::vector<std::uint8_t>& val_op) {

    ObsLil lil{};
    lil.calibrar(obs_d_cal, cfg_l, cfg_d, x_cal);

    for (std::size_t t = 0; t < x_op.size(); ++t) {
        const bool ok = (val_op[t] != 0u) && std::isfinite(x_op[t]);
        (void)lil.paso(salidas_delta[t], ok ? x_op[t] : NaN, ok, NaN, 0, 0);
    }

    lil.completar_nulos();

    ResumenEjecucionReal res{};
    const auto& fn = lil.ficha_nulo_cruda();
    const auto& fa = lil.resumen_armonico_cruda();
    const auto& mod = fn.modelo_lineal;
    res.sd_innov = std::sqrt(std::max(1e-15, mod.sigma2_innov));
    res.s_onsager = fn.contraste_onsager.estadistico_s;
    res.ee_s_onsager = fn.contraste_onsager.ee_nulo_s;
    res.p_onsager = fn.contraste_onsager.p_bilateral;
    res.p_holm_onsager = fa.p_holm[1];

    for (const auto& r : lil.cerradas_historico()) {
        EventoCapturado ev{};
        ev.t_apertura = r.t_apertura;
        ev.t_cambio = r.t_cambio;
        ev.escala_origen = r.escala_origen;
        ev.signo = r.signo;
        ev.desenlace = r.cruda.desenlace;
        ev.tau_rest = r.cruda.tau_rest_consolidado;
        ev.min_eta_2 = r.cruda.min_eta_desde_2;
        ev.p_det_ini = r.cruda.nulo_onsager.p_det_inicial;
        ev.u_pit_ons = r.cruda.nulo_onsager.u_pit;
        ev.delta_tau_ons = r.cruda.nulo_onsager.delta_tau;

        if (r.n_perfil_cruda > 0 && mod.valido) {
            ev.e_iny.assign(r.n_perfil_cruda, 0.0);
            ark::lilyth::innovaciones(
                mod,
                r.n_hist_prev_cruda > 0 ? r.hist_prev_cruda_est.data() : nullptr,
                r.n_hist_prev_cruda,
                r.perfil_cruda_est.data(),
                r.n_perfil_cruda,
                ev.e_iny.data());
        }
        if (r.cruda.desenlace == DesenlaceRanura::RESTITUIDA || r.cruda.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) {
            ++res.n_evaluados;
        }
        res.eventos[{r.t_apertura, r.escala_origen, r.signo}] = std::move(ev);
    }
    return res;
}

struct ComparacionContraLimpia {
    std::size_t k_limpios = 0;
    std::size_t k_emparejados_eval = 0;
    double delta_censura = 0.0;
    double rmse_e_iny_rel = 0.0;
    double p_det_limpia = 0.0;
    double p_det_cont = 0.0;
    double caida_p_det = 0.0;
    double sesgo_tau_rest = 0.0;
    double max_abs_dtau_rest = 0.0;
    double dif_min_eta2 = 0.0;
    double ds_sobre_ee = 0.0;
    bool pasa_estricto = false;
    bool pasa_dulce = false;
};

ComparacionContraLimpia comparar_ejecuciones(
    const ResumenEjecucionReal& ref,
    const ResumenEjecucionReal& cont,
    bool evaluar_e_iny_y_pdet,
    bool evaluar_tau_rest) {

    ComparacionContraLimpia out{};
    out.k_limpios = ref.n_evaluados;

    double sum_err2 = 0.0;
    std::size_t n_e = 0;
    double sum_pdet0 = 0.0, sum_pdetc = 0.0;
    std::size_t n_pdet0 = 0, n_pdetc = 0;
    double sum_dtau = 0.0, max_dtau = 0.0;
    std::size_t n_tau = 0;
    double sum_dmin = 0.0;
    std::size_t n_min = 0;

    for (const auto& [clave, ev0] : ref.eventos) {
        if (ev0.desenlace != DesenlaceRanura::RESTITUIDA && ev0.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
        if (std::isfinite(ev0.p_det_ini)) { sum_pdet0 += ev0.p_det_ini; ++n_pdet0; }

        auto it = cont.eventos.find(clave);
        if (it == cont.eventos.end()) continue;
        const auto& evc = it->second;

        if (evc.desenlace == DesenlaceRanura::RESTITUIDA || evc.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) {
            ++out.k_emparejados_eval;
        }
        if (std::isfinite(evc.p_det_ini)) { sum_pdetc += evc.p_det_ini; ++n_pdetc; }

        if (!ev0.e_iny.empty() && !evc.e_iny.empty()) {
            const std::size_t m = std::max(ev0.e_iny.size(), evc.e_iny.size());
            for (std::size_t k = 0; k < m; ++k) {
                // Alinear por el final (t_apertura) para no penalizar por recorte izquierdo si el choque final coincide
                const double v0 = (k < ev0.e_iny.size()) ? ev0.e_iny[ev0.e_iny.size() - 1 - k] : 0.0;
                const double vc = (k < evc.e_iny.size()) ? evc.e_iny[evc.e_iny.size() - 1 - k] : 0.0;
                const double d = (vc - v0) / ref.sd_innov;
                sum_err2 += d * d;
                ++n_e;
            }
        }

        if (std::isfinite(ev0.tau_rest) && std::isfinite(evc.tau_rest)) {
            const double dt = evc.tau_rest - ev0.tau_rest;
            sum_dtau += dt;
            max_dtau = std::max(max_dtau, std::abs(dt));
            ++n_tau;
        }
        if (std::isfinite(ev0.min_eta_2) && std::isfinite(evc.min_eta_2)) {
            sum_dmin += std::abs(evc.min_eta_2 - ev0.min_eta_2);
            ++n_min;
        }
    }

    out.delta_censura = (out.k_limpios > 0)
                            ? (1.0 - static_cast<double>(out.k_emparejados_eval) / static_cast<double>(out.k_limpios))
                            : 0.0;
    out.rmse_e_iny_rel = (n_e > 0) ? std::sqrt(sum_err2 / static_cast<double>(n_e)) : 0.0;
    out.p_det_limpia = (n_pdet0 > 0) ? (sum_pdet0 / static_cast<double>(n_pdet0)) : 0.0;
    out.p_det_cont = (n_pdetc > 0) ? (sum_pdetc / static_cast<double>(n_pdetc)) : 0.0;
    out.caida_p_det = std::max(0.0, out.p_det_limpia - out.p_det_cont);
    out.sesgo_tau_rest = (n_tau > 0) ? (sum_dtau / static_cast<double>(n_tau)) : NaN;
    out.max_abs_dtau_rest = max_dtau;
    out.dif_min_eta2 = (n_min > 0) ? (sum_dmin / static_cast<double>(n_min)) : 0.0;

    const double ee = (std::isfinite(ref.ee_s_onsager) && ref.ee_s_onsager > 1e-12) ? ref.ee_s_onsager : 0.05;
    out.ds_sobre_ee = (std::isfinite(cont.s_onsager) && std::isfinite(ref.s_onsager))
                          ? (std::abs(cont.s_onsager - ref.s_onsager) / ee)
                          : 99.0;

    // Verificar umbrales sellados en EXPECTATIVAS_COTA_NIVEL1.md (Sección 2)
    bool est = (out.delta_censura <= 0.10) && (out.ds_sobre_ee <= 0.25);
    bool dul = (out.delta_censura <  0.50) && (out.ds_sobre_ee <= 1.00);
    if (evaluar_e_iny_y_pdet) {
        est = est && (out.rmse_e_iny_rel <= 0.15) && (out.caida_p_det <= 0.05);
        dul = dul && (out.rmse_e_iny_rel <  0.50) && (out.caida_p_det <= 0.25);
    }
    if (evaluar_tau_rest) {
        est = est && std::isfinite(out.sesgo_tau_rest) && (std::abs(out.sesgo_tau_rest) <= 0.25);
        dul = dul && std::isfinite(out.sesgo_tau_rest) && (std::abs(out.sesgo_tau_rest) <= 1.50);
    }
    out.pasa_estricto = est;
    out.pasa_dulce = dul;
    return out;
}

} // namespace

int main() {
    ParametrosBancoComun params{};
    params.n_op = 6000; // 13 choques bien separados para que vía (b) corra rápido y con precisión
    const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, params, 202);

    ConfiguracionObservador cfg_d{};
    cfg_d.integracion = Integracion::FLUJO;
    cfg_d.agregacion = AgregacionEscala::SUMA;
    cfg_d.arl0_pasos_nativos = 450.0;
    cfg_d.alfa = ALFA_PROYECTO;
    cfg_d.medir_tau50 = false;
    cfg_d.escalas_max = 3;

    ObsDelta obs_d{};
    obs_d.calibrar(cfg_d, m2.x_cal);

    // Precalcular salidas del detector delta sobre la serie limpia para que todas las políticas se evalúen
    // exactamente sobre las mismas aperturas físicas (t_apertura, escala_origen, signo)
    obs_d.reiniciar_flujo();
    std::vector<ObsDelta::Salida> salidas_delta(m2.x_op_limpia.size());
    for (std::size_t t = 0; t < m2.x_op_limpia.size(); ++t) {
        salidas_delta[t] = obs_d.paso(m2.x_op_limpia[t], NaN, true, 0, 0);
    }

    const std::size_t m_freq = std::max<std::size_t>(32, static_cast<std::size_t>(std::pow(m2.x_cal.size(), 0.65)));
    const auto est_d = elw::estimar_calibracion(m2.x_cal.data(), m2.x_cal.size(), m_freq, -0.45, 0.45);

    ConfiguracionLilyth cfg_base{};
    cfg_base.alfa = ALFA_PROYECTO;
    cfg_base.b_max_nulo = 16;
    cfg_base.hilos_nulo = 4;
    cfg_base.semilla_nulo = 20261117ULL;
    cfg_base.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
    cfg_base.d_whittle_cruda = est_d.valido ? est_d.d : 0.25;
    cfg_base.forzar_via_b = true; // Fuerza ejecutar_via_b_evento real en cada ranura cerrada
    cfg_base.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cfg_base.b_mc_resumen = 1999;

    const std::vector<std::uint8_t> val_limpia(m2.x_op_limpia.size(), 1u);
    const auto ref_limpia = correr_lilyth_real(obs_d, cfg_d, cfg_base, m2.x_cal, salidas_delta, m2.x_op_limpia, val_limpia);

    std::printf("====================================================================================================================\n");
    std::printf(" PUNTOS 3, 4 Y 5: EVALUACION CON EL OPERADOR REAL (ObservadorLilyth + ejecutar_via_b_evento + avanzar)\n");
    std::printf(" Mundo M2 (ARFIMA_WHITTLE d=%.4f) | Eventos evaluados base K=%zu | S0=%.4f | ee_S=%.4f | p_det0=%.4f\n",
                cfg_base.d_whittle_cruda, ref_limpia.n_evaluados, ref_limpia.s_onsager, ref_limpia.ee_s_onsager,
                comparar_ejecuciones(ref_limpia, ref_limpia, true, true).p_det_limpia);
    std::printf("====================================================================================================================\n");

    // ========================================================================
    // 1. ZONA Z1 (hist_prev_est) A TRAVÉS DE ejecutar_via_b_evento REAL
    // ========================================================================
    std::printf("\n--- [E2.Z1] HIST_PREV_EST (t < t_cambio) via ejecutar_via_b_evento() real ---\n");
    std::printf(" dt | r | Politica Z1     | RMSE(e)/sd | p_det_b | dp_det | |dS|/ee_S | dCens | ESTRICTO | DULCE\n");
    std::printf("----+---+-----------------+------------+---------+--------+-----------+-------+----------+-------\n");

    for (std::size_t dt : {std::size_t{1}, std::size_t{2}, std::size_t{5}}) {
        for (std::size_t r : {std::size_t{1}, std::size_t{2}, std::size_t{3}}) {
            // Construir máscara con r inválidos en [t_cambio - dt - r + 1, t_cambio - dt] de cada evento
            // (sin pisar el ancla t <= t_fin_ancla)
            std::vector<std::uint8_t> val_z1(m2.x_op_limpia.size(), 1u);
            for (const auto& [cl, ev] : ref_limpia.eventos) {
                const std::int64_t idx_fin = ev.t_cambio - static_cast<std::int64_t>(dt);
                const std::int64_t idx_ini = idx_fin + 1 - static_cast<std::int64_t>(r);
                for (std::int64_t t = idx_ini; t <= idx_fin; ++t) {
                    if (t >= 0 && static_cast<std::size_t>(t) < val_z1.size()) val_z1[static_cast<std::size_t>(t)] = 0u;
                }
            }

            for (int p_i = 0; p_i < 3; ++p_i) {
                ConfiguracionLilyth cfg_z1 = cfg_base;
                cfg_z1.politica_z1_piloto = static_cast<ConfiguracionLilyth::PoliticaZ1Piloto>(p_i);
                cfg_z1.r_max_z1_piloto = 5;
                const char* nom = (p_i == 0) ? "COMPRESION_BASE" : (p_i == 1) ? "CORTE_SUFIJO" : "PUENTE_CAUSAL";
                const auto corr = correr_lilyth_real(obs_d, cfg_d, cfg_z1, m2.x_cal, salidas_delta, m2.x_op_limpia, val_z1);
                const auto cmp = comparar_ejecuciones(ref_limpia, corr, true, false);
                std::printf(" %2zu | %1zu | %-15s | %10.4f | %7.4f | %6.4f | %9.4f | %5.2f | %-8s | %-5s\n",
                            dt, r, nom, cmp.rmse_e_iny_rel, cmp.p_det_cont, cmp.caida_p_det,
                            cmp.ds_sobre_ee, cmp.delta_censura,
                            cmp.pasa_estricto ? "PASA" : "FALLA",
                            cmp.pasa_dulce    ? "PASA" : "FALLA");
            }
        }
        std::printf("----+---+-----------------+------------+---------+--------+-----------+-------+----------+-------\n");
    }

    // ========================================================================
    // 2. ZONA Z2 (perfil_est [t_cambio, t_a]) A TRAVÉS DE ejecutar_via_b_evento REAL
    // ========================================================================
    std::printf("\n--- [E2.Z2] PERFIL_EST ([t_cambio, t_apertura]) via ejecutar_via_b_evento() real ---\n");
    std::printf("  r | Politica Z2          | RMSE(e)/sd | p_det_b | dp_det | |dS|/ee_S | dCens | ESTRICTO | DULCE\n");
    std::printf("----+----------------------+------------+---------+--------+-----------+-------+----------+-------\n");

    for (std::size_t r : {std::size_t{0}, std::size_t{1}, std::size_t{2}, std::size_t{3}}) {
        std::vector<std::uint8_t> val_z2(m2.x_op_limpia.size(), 1u);
        if (r > 0) {
            for (const auto& [cl, ev] : ref_limpia.eventos) {
                for (std::size_t k = 0; k < r; ++k) {
                    const std::int64_t t = ev.t_cambio + static_cast<std::int64_t>(k);
                    // No invalidar t_apertura mismo para medir el efecto puro sobre perfil_est en via (b)
                    if (t >= 0 && t < ev.t_apertura && static_cast<std::size_t>(t) < val_z2.size()) {
                        val_z2[static_cast<std::size_t>(t)] = 0u;
                    }
                }
            }
        }

        for (int p_i = 0; p_i < 2; ++p_i) {
            ConfiguracionLilyth cfg_z2 = cfg_base;
            cfg_z2.politica_z2_piloto = static_cast<ConfiguracionLilyth::PoliticaZ2Piloto>(p_i);
            cfg_z2.r_max_z2_piloto = 5;
            const char* nom = (p_i == 0) ? "COMPRESION_BASE" : "PUENTE_CAUSAL";
            const auto corr = correr_lilyth_real(obs_d, cfg_d, cfg_z2, m2.x_cal, salidas_delta, m2.x_op_limpia, val_z2);
            const auto cmp = comparar_ejecuciones(ref_limpia, corr, true, false);
            std::printf(" %2zu | %-20s | %10.4f | %7.4f | %6.4f | %9.4f | %5.2f | %-8s | %-5s\n",
                        r, nom, cmp.rmse_e_iny_rel, cmp.p_det_cont, cmp.caida_p_det,
                        cmp.ds_sobre_ee, cmp.delta_censura,
                        cmp.pasa_estricto ? "PASA" : "FALLA",
                        cmp.pasa_dulce    ? "PASA" : "FALLA");
        }
    }

    // ========================================================================
    // 3. ZONA Z3 (ObservadorLilyth::avanzar real en t > t_apertura)
    //    Z3a: Hueco antes del paso de restitución (en tau = 1..r)
    //    Z3b: Hueco justo encima del paso exacto de restitución (en tau = tau_rest..tau_rest+r-1)
    // ========================================================================
    for (int subzona = 0; subzona < 2; ++subzona) {
        std::printf("\n--- [E2.%s] AVANZAR() REAL (%s) ---\n",
                    subzona == 0 ? "Z3a" : "Z3b",
                    subzona == 0 ? "Hueco en tau=1..r antes del cruce de 0.5"
                                 : "Hueco encima del paso exacto de restitucion tau=tau_rest..tau_rest+r-1");
        std::printf("  r | Politica Z3              | dCens | Sesgo(tau_r) | Max|dtau| | d(min_eta) | |dS|/ee_S | ESTRICTO | DULCE\n");
        std::printf("----+--------------------------+-------+--------------+-----------+------------+-----------+----------+-------\n");

        for (std::size_t r : {std::size_t{0}, std::size_t{1}, std::size_t{2}, std::size_t{3}, std::size_t{5}}) {
            std::vector<std::uint8_t> val_z3(m2.x_op_limpia.size(), 1u);
            if (r > 0) {
                for (const auto& [cl, ev] : ref_limpia.eventos) {
                    std::int64_t tau_ini = 1;
                    if (subzona == 1) {
                        if (!std::isfinite(ev.tau_rest)) continue;
                        tau_ini = static_cast<std::int64_t>(ev.tau_rest);
                    } else {
                        // Solo si el evento restituye después de r+1 para que Z3a sea verdaderamente antes del cruce
                        if (std::isfinite(ev.tau_rest) && ev.tau_rest <= static_cast<double>(r + 1)) continue;
                    }
                    for (std::size_t k = 0; k < r; ++k) {
                        const std::int64_t t = ev.t_apertura + tau_ini + static_cast<std::int64_t>(k);
                        if (t >= 0 && static_cast<std::size_t>(t) < val_z3.size()) {
                            val_z3[static_cast<std::size_t>(t)] = 0u;
                        }
                    }
                }
            }

            for (int p_i = 0; p_i < 3; ++p_i) {
                ConfiguracionLilyth cfg_z3 = cfg_base;
                cfg_z3.politica_z3_piloto = static_cast<ConfiguracionLilyth::PoliticaZ3Piloto>(p_i);
                cfg_z3.r_max_z3_piloto = 64; // Sin tope artificial para medir la distorsión física real en cada r
                const char* nom = (p_i == 0) ? "CORTA_AL_PRIMERO"
                                : (p_i == 1) ? "MUTANTE_COMPRIME_TAU"
                                             : "RELOJ_FISICO_SIN_COMPR";
                const auto corr = correr_lilyth_real(obs_d, cfg_d, cfg_z3, m2.x_cal, salidas_delta, m2.x_op_limpia, val_z3);
                const auto cmp = comparar_ejecuciones(ref_limpia, corr, false, true);
                std::printf(" %2zu | %-24s | %5.2f | %12.4f | %9.2f | %10.4f | %9.4f | %-8s | %-5s\n",
                            r, nom, cmp.delta_censura, cmp.sesgo_tau_rest, cmp.max_abs_dtau_rest,
                            cmp.dif_min_eta2, cmp.ds_sobre_ee,
                            cmp.pasa_estricto ? "PASA" : "FALLA",
                            cmp.pasa_dulce    ? "PASA" : "FALLA");
            }
        }
    }

    return 0;
}
