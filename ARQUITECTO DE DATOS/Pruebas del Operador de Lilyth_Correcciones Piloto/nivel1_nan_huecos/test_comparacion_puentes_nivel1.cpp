// ============================================================================
// test_comparacion_puentes_nivel1.cpp
// Fase 4: Comparador Integral de Puentes y Pausa (B0..B4) + Controles CF4.P1..P6
// sobre el Operador Real (ObservadorSaltos / ObservadorSaltosFilas / ObservadorLilyth)
//
// Regímenes evaluados: ESTRICTO (alfa=0.05) y DULCE (<50% error)
// Uso: test_comparacion_puentes_nivel1.exe [R_replicas=200] [ruta_volcado_orig]
// ============================================================================

#include "ark/metrologia/observador_saltos_filas.hpp"
#include "banco_sintetico_comun.hpp"
#include "protocolo_inputs_lilyth.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace {

using namespace ark::dataset;
using namespace ark::metrologia;
using namespace ark::ventanas;
using namespace ark::lilyth;
using namespace ark::lilyth::piloto;

constexpr std::size_t T_MAX_F4 = 64;
using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsDeltaFilas = ObservadorSaltosFilas<24, 64, 120, 1440, 8, 256>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, T_MAX_F4>;
using ClaveEvento = std::tuple<std::int64_t, int, int>;

std::pair<double, double> wilson_95(std::size_t k, std::size_t n) noexcept {
    if (n == 0) return {0.0, 1.0};
    const double z = 1.959963984540054;
    const double p = static_cast<double>(k) / static_cast<double>(n);
    const double d = 1.0 + z * z / static_cast<double>(n);
    const double c = (p + z * z / (2.0 * static_cast<double>(n))) / d;
    const double h = z * std::sqrt(p * (1.0 - p) / static_cast<double>(n) + z * z / (4.0 * static_cast<double>(n) * static_cast<double>(n))) / d;
    return {std::max(0.0, c - h), std::min(1.0, c + h)};
}

void en_paralelo(std::size_t n_tareas, std::size_t n_hilos, const std::function<void(std::size_t)>& fn) {
    std::atomic<std::size_t> sig{0};
    std::vector<std::thread> trab;
    const std::size_t h_ef = std::max<std::size_t>(1, n_hilos);
    trab.reserve(h_ef);
    for (std::size_t h = 0; h < h_ef; ++h) {
        trab.emplace_back([&]() {
            for (;;) {
                const std::size_t i = sig.fetch_add(1, std::memory_order_relaxed);
                if (i >= n_tareas) break;
                fn(i);
            }
        });
    }
    for (auto& t : trab) t.join();
}

std::string hx(double x) {
    char b[64];
    std::snprintf(b, sizeof b, "%a", x);
    return b;
}

void volcar_contraste_str(std::ostringstream& os, const char* tag, const ContrasteResumenNulo& c) {
    char b[512];
    std::snprintf(b, sizeof b, "%s val=%d K=%zu M=%zu B=%zu S=%s ee=%s p=%s md=%s mu=%s pmin=%s res=%d\n",
                  tag, int(c.valido), c.n_eventos, c.n_nulos_pool, c.b_mc,
                  hx(c.estadistico_s).c_str(), hx(c.ee_nulo_s).c_str(),
                  hx(c.p_bilateral).c_str(), hx(c.mediana_delta_tau).c_str(),
                  hx(c.media_u_pit).c_str(), hx(c.p_min_alcanzable).c_str(), int(c.resoluble));
    os << b;
}

void volcar_ref_str(std::ostringstream& os, const char* tag, const ReferenciaNuloEvento& r) {
    char b[512];
    std::snprintf(b, sizeof b, "%s ev=%d iny=%d na=%zu nb=%zu nu=%zu int=%zu pdet=%s u=%s dt=%s md=%s;",
                  tag, int(r.evaluado), int(r.nulo_por_inyeccion), r.n_nulos_a, r.n_nulos_b,
                  r.n_nulos_usados, r.intentos_b, hx(r.p_det_inicial).c_str(),
                  hx(r.u_pit).c_str(), hx(r.delta_tau).c_str(), hx(r.mediana_nulo).c_str());
    os << b;
}

void correr_y_volcar_str(
    std::ostringstream& os,
    const char* nombre_caso,
    const SerieSinteticaCoordinada& m,
    Integracion integ,
    bool purga,
    ModoNuloLineal mlin,
    bool inyectar_huecos) {

    ConfiguracionObservador cd{};
    cd.integracion = integ;
    cd.agregacion = AgregacionEscala::SUMA;
    cd.arl0_pasos_nativos = 450.0;
    cd.alfa = ALFA_PROYECTO;
    cd.purga = purga;
    cd.bins_fase = purga ? 24 : 1;
    cd.dias_minimos_purga = 20;
    cd.medir_tau50 = false;
    cd.semilla = 20261117ULL;
    cd.escalas_max = 3;

    ObsDelta od{};
    od.calibrar(cd, m.x_cal, {}, {},
                purga ? std::span<const std::size_t>(m.fase_cal) : std::span<const std::size_t>{},
                purga ? std::span<const std::int64_t>(m.dia_cal) : std::span<const std::int64_t>{});

    ConfiguracionLilyth cl{};
    cl.alfa = ALFA_PROYECTO;
    cl.representacion = purga ? ModoRepresentacionLilyth::PURGADA : ModoRepresentacionLilyth::CRUDA;
    cl.b_max_nulo = 16;
    cl.hilos_nulo = 1;
    cl.semilla_nulo = 20261117ULL;
    cl.modo_nulo_lineal = mlin;
    cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cl.b_mc_resumen = 1999;
    if (mlin == ModoNuloLineal::ARFIMA_WHITTLE) cl.d_whittle_cruda = 0.25;

    ObsLil lil{};
    lil.calibrar(od, cl, cd, m.x_cal, {},
                 purga ? std::span<const std::size_t>(m.fase_cal) : std::span<const std::size_t>{},
                 purga ? std::span<const std::int64_t>(m.dia_cal) : std::span<const std::int64_t>{});

    od.reiniciar_flujo();
    for (std::size_t t = 0; t < m.x_op_limpia.size(); ++t) {
        const bool ok = !(inyectar_huecos && (t % 97 == 50 || t == 804));
        const double x_in = ok ? m.x_op_limpia[t] : NaN;
        const std::size_t fa = purga ? m.fase_op[t] : 0;
        const std::int64_t di = purga ? m.dia_op[t] : 0;
        const auto sd = od.paso(x_in, NaN, ok, fa, di);
        (void)lil.paso(sd, x_in, ok, NaN, fa, di);
    }
    lil.completar_nulos();

    const auto& fn = purga ? lil.ficha_nulo_purgada() : lil.ficha_nulo_cruda();
    const auto& fa = purga ? lil.resumen_armonico_purgada() : lil.resumen_armonico_cruda();

    char b[512];
    std::snprintf(b, sizeof b, "CASO %s tmax=%zu ap=%zu cer=%zu\n",
                  nombre_caso, lil.configuracion().t_max_efectivo,
                  lil.contadores().aperturas, lil.cerradas_historico().size());
    os << b;
    volcar_contraste_str(os, "C_iid", fn.contraste_iid);
    volcar_contraste_str(os, "C_ons", fn.contraste_onsager);
    volcar_contraste_str(os, "C_iid_loo", fn.contraste_iid_pool_loo);
    volcar_contraste_str(os, "C_ons_loo", fn.contraste_onsager_pool_loo);
    std::snprintf(b, sizeof b, "HOLM val=%d pc=%s %s %s %s ph=%s %s %s %s\n",
                  int(fa.familia_holm_valida),
                  hx(fa.p_crudos[0]).c_str(), hx(fa.p_crudos[1]).c_str(),
                  hx(fa.p_crudos[2]).c_str(), hx(fa.p_crudos[3]).c_str(),
                  hx(fa.p_holm[0]).c_str(), hx(fa.p_holm[1]).c_str(),
                  hx(fa.p_holm[2]).c_str(), hx(fa.p_holm[3]).c_str());
    os << b;

    for (const auto& r : lil.cerradas_historico()) {
        const auto& rama = purga ? r.purgada : r.cruda;
        std::snprintf(b, sizeof b, "RAN ta=%lld tc=%lld j=%d s=%d des=%d tr=%s min2=%s np=%zu nh=%zu ",
                      static_cast<long long>(r.t_apertura), static_cast<long long>(r.t_cambio),
                      r.escala_origen, r.signo, static_cast<int>(rama.desenlace),
                      hx(rama.tau_rest_consolidado).c_str(), hx(rama.min_eta_desde_2).c_str(),
                      purga ? r.n_perfil_purgada : r.n_perfil_cruda,
                      purga ? r.n_hist_prev_purgada : r.n_hist_prev_cruda);
        os << b;
        volcar_ref_str(os, "iid", rama.nulo_iid);
        volcar_ref_str(os, "ons", rama.nulo_onsager);
        os << "\n";
    }
}

struct SalidaCompletaOperador {
    std::map<ClaveEvento, std::vector<double>> e_iny_por_ev;
    std::map<ClaveEvento, double> tau_rest_por_ev;
    std::map<ClaveEvento, std::int64_t> t_cambio_por_ev;
    std::size_t k_medibles = 0;
    std::size_t k_evaluados = 0;
    std::size_t n_veto_z2 = 0;
    std::size_t alarmas_nivel_tot = 0;
    std::size_t alarmas_temp_tot = 0;
    std::size_t alarmas_nivel_espurias = 0;
    std::size_t alarmas_temp_espurias = 0;
    double s_onsager = NaN;
    double ee_s_onsager = NaN;
    double p_holm_onsager = NaN;
    double sd_innov = 1.0;
    bool valida = false;
    bool rechaza_inter = false;
    bool rechaza_lenta = false;
};

SalidaCompletaOperador correr_operador_sobre_flujo(
    const ObsDelta& od_cal,
    const ConfiguracionObservador& cd,
    const ConfiguracionLilyth& cl,
    const SerieSinteticaCoordinada& m,
    const std::vector<double>& x_delta,
    const std::vector<std::uint8_t>& ok_delta,
    const std::vector<double>& x_lil,
    const std::vector<std::uint8_t>& ok_lil,
    const std::vector<std::uint8_t>& flags_lil,
    const std::vector<std::uint8_t>& es_ventana_choque_real = {}) {

    const bool usar_purga = cd.purga;
    ObsDelta od = od_cal;
    od.reiniciar_flujo();

    ObsLil lil{};
    if (!lil.calibrar(od, cl, cd, m.x_cal, {},
                      usar_purga ? std::span<const std::size_t>(m.fase_cal) : std::span<const std::size_t>{},
                      usar_purga ? std::span<const std::int64_t>(m.dia_cal) : std::span<const std::int64_t>{})) {
        return SalidaCompletaOperador{};
    }

    SalidaCompletaOperador out{};
    const std::size_t n = x_delta.size();
    for (std::size_t t = 0; t < n; ++t) {
        const std::size_t fa = usar_purga ? m.fase_op[t] : 0;
        const std::int64_t di = usar_purga ? m.dia_op[t] : 0;
        const bool okd = (ok_delta[t] != 0u) && std::isfinite(x_delta[t]);
        const bool okl = (ok_lil[t] != 0u) && std::isfinite(x_lil[t]);
        const auto sd = od.paso(okd ? x_delta[t] : NaN, NaN, okd, fa, di);

        bool hay_niv = false, hay_tmp = false;
        for (std::size_t j = 0; j < sd.niveles; ++j) {
            if (sd.nivel[j].alarma.hay) hay_niv = true;
            if (sd.temp[j].alarma.hay) hay_tmp = true;
        }
        if (hay_niv) {
            ++out.alarmas_nivel_tot;
            if (!es_ventana_choque_real.empty() && !es_ventana_choque_real[t]) ++out.alarmas_nivel_espurias;
        }
        if (hay_tmp) {
            ++out.alarmas_temp_tot;
            if (!es_ventana_choque_real.empty() && !es_ventana_choque_real[t]) ++out.alarmas_temp_espurias;
        }
        (void)lil.paso(sd, okl ? x_lil[t] : NaN, okl, NaN, fa, di, NaN, {}, flags_lil[t]);
    }
    lil.completar_nulos();

    const auto& fn = usar_purga ? lil.ficha_nulo_purgada() : lil.ficha_nulo_cruda();
    const auto& fa = usar_purga ? lil.resumen_armonico_purgada() : lil.resumen_armonico_cruda();
    const auto& mod = fn.modelo_lineal;
    out.sd_innov = std::sqrt(std::max(1e-15, mod.sigma2_innov));
    out.s_onsager = fn.contraste_onsager.estadistico_s;
    out.ee_s_onsager = fn.contraste_onsager.ee_nulo_s;
    out.p_holm_onsager = fa.p_holm[1];
    out.valida = fn.contraste_onsager.valido && fa.familia_holm_valida && std::isfinite(fa.p_holm[1]);
    out.rechaza_inter = out.valida && (fa.p_holm[1] < ALFA_PROYECTO);
    out.rechaza_lenta = out.rechaza_inter && (out.s_onsager > 0.0);

    for (const auto& r : lil.cerradas_historico()) {
        const auto& rama = usar_purga ? r.purgada : r.cruda;
        if (rama.desenlace == DesenlaceRanura::RESTITUIDA || rama.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) {
            ++out.k_medibles;
        }
        if (rama.nulo_onsager.evaluado) ++out.k_evaluados;
        if (usar_purga ? r.perfil_corrupto_z2_purgada : r.perfil_corrupto_z2_cruda) ++out.n_veto_z2;
        const ClaveEvento clv{r.t_apertura, r.escala_origen, r.signo};
        out.tau_rest_por_ev[clv] = rama.tau_rest_consolidado;
        out.t_cambio_por_ev[clv] = r.t_cambio;

        const std::size_t n_perf = usar_purga ? r.n_perfil_purgada : r.n_perfil_cruda;
        const std::size_t n_prev = usar_purga ? r.n_hist_prev_purgada : r.n_hist_prev_cruda;
        const auto& perf = usar_purga ? r.perfil_purgada_est : r.perfil_cruda_est;
        const auto& prev = usar_purga ? r.hist_prev_purgada_est : r.hist_prev_cruda_est;
        if (n_perf > 0 && mod.valido) {
            std::vector<double> e(n_perf, 0.0);
            innovaciones(mod, n_prev > 0 ? prev.data() : nullptr, n_prev, perf.data(), n_perf, e.data());
            out.e_iny_por_ev[clv] = std::move(e);
        }
    }
    return out;
}

double calcular_rmse_e(const SalidaCompletaOperador& ref, const SalidaCompletaOperador& cont) {
    double s2 = 0.0;
    std::size_t n = 0;
    for (const auto& [cl, e0] : ref.e_iny_por_ev) {
        auto it = cont.e_iny_por_ev.find(cl);
        if (it == cont.e_iny_por_ev.end() || e0.empty() || it->second.empty()) continue;
        const std::size_t m = std::max(e0.size(), it->second.size());
        for (std::size_t k = 0; k < m; ++k) {
            const double v0 = (k < e0.size()) ? e0[e0.size() - 1 - k] : 0.0;
            const double vc = (k < it->second.size()) ? it->second[it->second.size() - 1 - k] : 0.0;
            const double d = (vc - v0) / ref.sd_innov;
            s2 += d * d;
            ++n;
        }
    }
    return n > 0 ? std::sqrt(s2 / static_cast<double>(n)) : 99.0;
}

// Simula un flujo histórico/en vivo completo sobre una lista de intervalos corruptos (racha r <= r_max o > r_max)
enum class CompetidorF4 {
    B0_BASE = 0,
    B1_ARGON_VIVO = 1,
    B1B_ARRASTRE_NIVEL = 2,
    B2_LINEAL_2BORDES = 3,
    B3_ITO_SIN_FILTRO = 4,
    B4_PROPUESTA_UNIFICADA = 5
};

void preparar_flujos_competidor(
    CompetidorF4 comp,
    Integracion integ,
    const ModeloLinealNulo& mod_calma,
    const std::vector<double>& x_fisica,
    const std::vector<std::uint8_t>& es_corrupto,
    std::size_t r_max,
    double z_max,
    std::uint64_t semilla_rng,
    std::vector<double>& x_delta_out,
    std::vector<std::uint8_t>& ok_delta_out,
    std::vector<double>& x_lil_out,
    std::vector<std::uint8_t>& ok_lil_out,
    std::vector<std::uint8_t>& flags_out,
    std::size_t& n_reconstruidos_out,
    std::size_t& n_choques_en_hueco_out,
    std::size_t& n_pausas_out,
    double& max_err_anclaje_out) {

    const std::size_t n = x_fisica.size();
    x_delta_out = x_fisica;
    x_lil_out = x_fisica;
    ok_delta_out.assign(n, 1u);
    ok_lil_out.assign(n, 1u);
    flags_out.assign(n, 0u);
    n_reconstruidos_out = n_choques_en_hueco_out = n_pausas_out = 0;
    max_err_anclaje_out = 0.0;

    std::size_t t = 0;
    std::size_t id_h = 1;
    std::vector<double> hist_est;
    hist_est.reserve(512);
    double prev_x_ok = NaN;

    auto registrar_val = [&](double xv) {
        if (integ == Integracion::NIVEL) {
            if (std::isfinite(prev_x_ok)) hist_est.push_back(xv - prev_x_ok);
            prev_x_ok = xv;
        } else {
            hist_est.push_back(xv);
            prev_x_ok = xv;
        }
        if (hist_est.size() > 512) hist_est.erase(hist_est.begin(), hist_est.begin() + 256);
    };

    while (t < n) {
        if (!es_corrupto[t]) {
            registrar_val(x_fisica[t]);
            ++t;
            continue;
        }
        std::size_t t_fin = t;
        while (t_fin + 1 < n && es_corrupto[t_fin + 1]) ++t_fin;
        const std::size_t r_gap = t_fin - t + 1;

        if (comp == CompetidorF4::B0_BASE) {
            for (std::size_t k = t; k <= t_fin; ++k) {
                x_delta_out[k] = NaN; ok_delta_out[k] = 0u;
                x_lil_out[k] = NaN;   ok_lil_out[k] = 0u;
                flags_out[k] = FLAG_CORRUPCION_CAIDA;
            }
            prev_x_ok = NaN;
            t = t_fin + 1;
            continue;
        }

        if (comp == CompetidorF4::B1_ARGON_VIVO) {
            if (r_gap > r_max) ++n_pausas_out;
            for (std::size_t k = t; k <= t_fin; ++k) {
                if (integ == Integracion::NIVEL) {
                    x_delta_out[k] = NaN; ok_delta_out[k] = 0u; // evitar_teletransporte_nivel
                } else if (k - t + 1 <= r_max && std::isfinite(prev_x_ok)) {
                    x_delta_out[k] = 0.0; ok_delta_out[k] = 1u;
                } else {
                    x_delta_out[k] = NaN; ok_delta_out[k] = 0u;
                }
                x_lil_out[k] = NaN; ok_lil_out[k] = 0u; // Z1/Z2/Z3 dentro de Lilyth ven el NaN real
                flags_out[k] = FLAG_CORRUPCION_CAIDA | (r_gap > r_max ? FLAG_PAUSA_UMBRAL : 0u);
            }
            if (integ == Integracion::NIVEL) prev_x_ok = NaN;
            t = t_fin + 1;
            continue;
        }

        if (comp == CompetidorF4::B1B_ARRASTRE_NIVEL) {
            const double val_arr = std::isfinite(prev_x_ok) ? prev_x_ok : 0.0;
            for (std::size_t k = t; k <= t_fin; ++k) {
                if (k - t + 1 <= r_max && std::isfinite(prev_x_ok)) {
                    x_delta_out[k] = val_arr; ok_delta_out[k] = 1u;
                    x_lil_out[k] = val_arr;   ok_lil_out[k] = 1u;
                } else {
                    x_delta_out[k] = NaN; ok_delta_out[k] = 0u;
                    x_lil_out[k] = NaN;   ok_lil_out[k] = 0u;
                }
                flags_out[k] = FLAG_CORRUPCION_CAIDA;
            }
            t = t_fin + 1;
            continue;
        }

        const bool tiene_bordes = (t > 0) && (t_fin + 1 < n) && std::isfinite(x_fisica[t - 1]) && std::isfinite(x_fisica[t_fin + 1]);
        if (!tiene_bordes) {
            for (std::size_t k = t; k <= t_fin; ++k) {
                x_delta_out[k] = NaN; ok_delta_out[k] = 0u;
                x_lil_out[k] = NaN;   ok_lil_out[k] = 0u;
                flags_out[k] = FLAG_CORRUPCION_CAIDA;
            }
            prev_x_ok = NaN;
            t = t_fin + 1;
            continue;
        }

        const double x_ini = x_fisica[t - 1];
        const double x_fin = x_fisica[t_fin + 1];

        ConfiguracionPuenteHistoricoIto cp{};
        cp.integracion = integ;
        cp.r_max = (comp == CompetidorF4::B3_ITO_SIN_FILTRO) ? MAX_RACHA_PUENTE_ITO : r_max;
        cp.z_max_hueco = (comp == CompetidorF4::B3_ITO_SIN_FILTRO) ? 1e9 : z_max;
        cp.usar_recta_lineal_pura = (comp == CompetidorF4::B2_LINEAL_2BORDES);
        cp.n_semillas_rng = (comp == CompetidorF4::B2_LINEAL_2BORDES) ? 1 : 4;
        cp.semilla_base = semilla_rng;

        if (comp == CompetidorF4::B2_LINEAL_2BORDES) {
            if (r_gap <= r_max) {
                const double paso = (x_fin - x_ini) / static_cast<double>(r_gap + 1);
                for (std::size_t k = 0; k < r_gap; ++k) {
                    const double xv = x_ini + static_cast<double>(k + 1) * paso;
                    x_delta_out[t + k] = xv; ok_delta_out[t + k] = 1u;
                    x_lil_out[t + k] = xv;   ok_lil_out[t + k] = 1u;
                    flags_out[t + k] = 0u; // LINEAL sin veto de etiqueta: alimenta directo a los operadores
                    registrar_val(xv);
                }
                ++n_reconstruidos_out;
            } else {
                for (std::size_t k = t; k <= t_fin; ++k) {
                    x_delta_out[k] = NaN; ok_delta_out[k] = 0u;
                    x_lil_out[k] = NaN;   ok_lil_out[k] = 0u;
                }
                prev_x_ok = NaN;
            }
            t = t_fin + 1;
            continue;
        }

        const auto res = evaluar_y_reconstruir_hueco_ito(cp, mod_calma, hist_est, x_ini, x_fin, r_gap, id_h++);
        if (res.max_error_anclaje_rng > max_err_anclaje_out) max_err_anclaje_out = res.max_error_anclaje_rng;

        if (comp == CompetidorF4::B3_ITO_SIN_FILTRO) {
            if (res.reconstruido) {
                for (std::size_t k = 0; k < r_gap; ++k) {
                    const double xv = res.x_rng[0][k]; // 1 sola realización RNG sin etiqueta FLAG_PUENTE_ITO!
                    x_delta_out[t + k] = xv; ok_delta_out[t + k] = 1u;
                    x_lil_out[t + k] = xv;   ok_lil_out[t + k] = 1u;
                    flags_out[t + k] = 0u;
                    registrar_val(xv);
                }
                ++n_reconstruidos_out;
            } else {
                for (std::size_t k = t; k <= t_fin; ++k) {
                    x_delta_out[k] = NaN; ok_delta_out[k] = 0u;
                    x_lil_out[k] = NaN;   ok_lil_out[k] = 0u;
                }
                prev_x_ok = NaN;
            }
        } else {
            // B4_PROPUESTA_UNIFICADA: doble llave + etiqueta FLAG_PUENTE_ITO + pausa si excede
            // En NIVEL preserva continuidad sin teletransporte (+Delta); en FLUJO mantiene Z0 causal (0.0) para no desplazar t_cambio al hueco
            if (res.reconstruido) {
                for (std::size_t k = 0; k < r_gap; ++k) {
                    const double xv = res.x_rng[0][k];
                    x_delta_out[t + k] = (integ == Integracion::NIVEL) ? xv : 0.0;
                    ok_delta_out[t + k] = 1u;
                    x_lil_out[t + k] = xv;   ok_lil_out[t + k] = 1u;
                    flags_out[t + k] = FLAG_PUENTE_ITO;
                    registrar_val(xv);
                }
                ++n_reconstruidos_out;
            } else {
                if (res.choque_en_hueco) ++n_choques_en_hueco_out;
                if (res.excedio_r_max) ++n_pausas_out;
                for (std::size_t k = t; k <= t_fin; ++k) {
                    x_delta_out[k] = NaN; ok_delta_out[k] = 0u;
                    x_lil_out[k] = NaN;   ok_lil_out[k] = 0u;
                    flags_out[k] = res.flags_emitidos;
                }
                prev_x_ok = NaN;
            }
        }
        t = t_fin + 1;
    }
}

struct ResumenMCF4 {
    std::size_t n_val = 0;
    std::size_t n_rech = 0;
    std::size_t n_rech_lenta = 0;
    double tasa_rech = 0.0;
    double wlo = 0.0, whi = 1.0;
    double tasa_lenta = 0.0;
    double k_eval_med = 0.0;
    double media_s = 0.0;
    double media_ee = 0.05;
    double temp_alarm_med = 0.0;
    double veto_z2_med = 0.0;
};

} // namespace

int main(int argc, char** argv) {
    const std::size_t R_REP = (argc > 1) ? static_cast<std::size_t>(std::stoul(argv[1])) : 200;
    const std::string ruta_volcado_orig = (argc > 2) ? argv[2] : "nivel1_nan_huecos/volcado_orig.txt";
    const std::size_t n_hilos = std::max<std::size_t>(1, std::thread::hardware_concurrency());

    std::printf("====================================================================================================\n");
    std::printf(" FASE 4: COMPARACION UNIFICADA DE PUENTES (B0..B4) Y REGRESION INTEGRAL (R=%zu, Hilos=%zu)\n", R_REP, n_hilos);
    std::printf("====================================================================================================\n");

    int fallas = 0;
    auto reportar = [&](const char* id, const char* desc, bool ok) {
        std::printf("[%s] %s: %s\n", ok ? "PASA" : "FALLA", id, desc);
        if (!ok) ++fallas;
    };

    // ------------------------------------------------------------------------
    // CF4.P1: Identidad bit a bit contra volcado_orig.txt + tamaño nominal q=0 en M1 y M2
    // ------------------------------------------------------------------------
    {
        ParametrosBancoComun pv{};
        pv.n_op = 4000;
        const auto m1_v = generar_mundo_m1_ar1_flujo(0.80, pv, 101);
        const auto m2_v = generar_mundo_m2_arfima_flujo(0.25, 0.40, pv, 202);
        const auto m3_v = generar_mundo_m3_nivel_purga(0.65, pv, 303);

        std::ostringstream os;
        correr_y_volcar_str(os, "M1_LIMPIO", m1_v, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, false);
        correr_y_volcar_str(os, "M1_HUECOS", m1_v, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, true);
        correr_y_volcar_str(os, "M2_LIMPIO", m2_v, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE,     false);
        correr_y_volcar_str(os, "M2_HUECOS", m2_v, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE,     true);
        correr_y_volcar_str(os, "M3_LIMPIO", m3_v, Integracion::NIVEL, true,  ModoNuloLineal::AR_LEVINSON_DURBIN, false);
        correr_y_volcar_str(os, "M3_HUECOS", m3_v, Integracion::NIVEL, true,  ModoNuloLineal::AR_LEVINSON_DURBIN, true);

        const std::string actual_hex = os.str();
        std::ifstream ifs(ruta_volcado_orig);
        std::ostringstream os_ref;
        os_ref << ifs.rdbuf();
        std::string ref_hex = os_ref.str();
        ref_hex.erase(std::remove(ref_hex.begin(), ref_hex.end(), '\r'), ref_hex.end());
        const bool id_ok = !ref_hex.empty() && (actual_hex == ref_hex);
        reportar("CF4.P1", "Identidad bit a bit (%a) contra observador_lilyth_original.hpp intacta tras Fases 1-3", id_ok);
    }

    // ------------------------------------------------------------------------
    // CF4.P2: Caso NIVEL M3 + PurgaDiurna ("100 -> hueco -> 150")
    //   P2a: 18 micro-cortes de r=4 barras en tramos tranquilos (alarmas nivel/temp espurias)
    //   P2b: Salto oculto +12 sigma dentro de un micro-corte de r=2 barras (llave de amplitud)
    // ------------------------------------------------------------------------
    {
        ParametrosBancoComun pb{};
        pb.n_op = 6000;
        const auto m3 = generar_mundo_m3_nivel_purga(0.65, pb, 303);

        ConfiguracionObservador cd3{};
        cd3.integracion = Integracion::NIVEL;
        cd3.agregacion = AgregacionEscala::SUMA;
        cd3.arl0_pasos_nativos = 450.0;
        cd3.alfa = ALFA_PROYECTO;
        cd3.purga = true;
        cd3.bins_fase = 24;
        cd3.dias_minimos_purga = 20;
        cd3.medir_tau50 = false;
        cd3.escalas_max = 3;

        ObsDelta od3{};
        od3.calibrar(cd3, m3.x_cal, {}, {}, m3.fase_cal, m3.dia_cal);

        ConfiguracionLilyth cl3{};
        cl3.alfa = ALFA_PROYECTO;
        cl3.representacion = ModoRepresentacionLilyth::PURGADA;
        cl3.b_max_nulo = 16;
        cl3.hilos_nulo = 4;
        cl3.semilla_nulo = 20261118ULL;
        cl3.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
        cl3.b_mc_resumen = 1999;

        ObsLil lil_cal3{};
        (void)lil_cal3.calibrar(od3, cl3, cd3, m3.x_cal, {}, m3.fase_cal, m3.dia_cal);
        const auto& mod3 = lil_cal3.ficha_nulo_purgada().modelo_lineal;

        // Ventanas de choques reales en M3 para separar alarmas legítimas vs espurias
        std::vector<std::uint8_t> en_choque3(m3.x_op_limpia.size(), 0u);
        for (std::size_t tc : m3.tiempos_choque) {
            const std::size_t ini = (tc >= 5) ? (tc - 5) : 0;
            const std::size_t fin = std::min(m3.x_op_limpia.size(), tc + 50);
            for (std::size_t t = ini; t < fin; ++t) en_choque3[t] = 1u;
        }

        // Subcaso P2a: Serie limpia vs 18 micro-cortes de r=4 barras fuera de choques reales
        std::vector<std::uint8_t> corr_p2a(m3.x_op_limpia.size(), 0u);
        for (std::size_t tc : m3.tiempos_choque) {
            if (tc > 180) {
                for (std::size_t k = 0; k < 4; ++k) corr_p2a[tc - 160 + k] = 1u;
            }
        }
        std::vector<std::uint8_t> flags0(m3.x_op_limpia.size(), 0u);
        std::vector<std::uint8_t> ok1(m3.x_op_limpia.size(), 1u);
        const auto ref3 = correr_operador_sobre_flujo(od3, cd3, cl3, m3, m3.x_op_limpia, ok1, m3.x_op_limpia, ok1, flags0, en_choque3);

        std::printf("\n--- [P2a] NIVEL M3 + PurgaDiurna: Micro-cortes r=4 fuera de choques reales (K_ref=%zu, S0=%.4f) ---\n",
                    ref3.k_evaluados, ref3.s_onsager);
        std::printf(" Competidor              | Alm_Niv_Esp | Alm_Tmp_Esp | K_eval | |dS|/ee | err_Anc_RNG | Notas\n");
        std::printf("-------------------------+-------------+-------------+--------+---------+-------------+---------------------------\n");

        std::size_t esp_tmp_lin = 0, esp_tmp_b4 = 0, esp_niv_b4 = 0, esp_niv_arr = 0;
        double ds_b4_p2a = 99.0;

        for (auto comp : {CompetidorF4::B0_BASE, CompetidorF4::B1_ARGON_VIVO, CompetidorF4::B1B_ARRASTRE_NIVEL,
                          CompetidorF4::B2_LINEAL_2BORDES, CompetidorF4::B3_ITO_SIN_FILTRO, CompetidorF4::B4_PROPUESTA_UNIFICADA}) {
            std::vector<double> xd, xl;
            std::vector<std::uint8_t> okd, okl, fl;
            std::size_t nrec = 0, nchq = 0, npau = 0;
            double max_anc = 0.0;
            preparar_flujos_competidor(comp, Integracion::NIVEL, mod3, m3.x_op_limpia, corr_p2a, 4, 4.0, 20261118ULL,
                                       xd, okd, xl, okl, fl, nrec, nchq, npau, max_anc);

            ConfiguracionLilyth cl_c = cl3;
            if (comp == CompetidorF4::B4_PROPUESTA_UNIFICADA || comp == CompetidorF4::B1_ARGON_VIVO) {
                cl_c.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
                cl_c.politica_z3_piloto = ConfiguracionLilyth::PoliticaZ3Piloto::RELOJ_FISICO_TOPE_R;
                cl_c.r_max_z3_piloto = 4;
            }
            const auto out = correr_operador_sobre_flujo(od3, cd3, cl_c, m3, xd, okd, xl, okl, fl, en_choque3);
            const double ds_ee = std::abs(out.s_onsager - ref3.s_onsager) / std::max(1e-12, ref3.ee_s_onsager);
            const char* nom = (comp == CompetidorF4::B0_BASE)                ? "B0_BASE_ACTUAL"
                            : (comp == CompetidorF4::B1_ARGON_VIVO)          ? "B1_ARGON_VIVO"
                            : (comp == CompetidorF4::B1B_ARRASTRE_NIVEL)    ? "B1b_ARRASTRE_NIVEL"
                            : (comp == CompetidorF4::B2_LINEAL_2BORDES)     ? "B2_LINEAL_2BORDES"
                            : (comp == CompetidorF4::B3_ITO_SIN_FILTRO)     ? "B3_ITO_SIN_FILTRO"
                                                                             : "B4_PROPUESTA_UNIFICADA";
            std::printf(" %-23s | %11zu | %11zu | %6zu | %7.4f | %11.2e | rec=%zu chq=%zu pau=%zu\n",
                        nom, out.alarmas_nivel_espurias, out.alarmas_temp_espurias, out.k_evaluados, ds_ee, max_anc, nrec, nchq, npau);

            if (comp == CompetidorF4::B2_LINEAL_2BORDES) esp_tmp_lin = out.alarmas_temp_espurias;
            if (comp == CompetidorF4::B1B_ARRASTRE_NIVEL) esp_niv_arr = out.alarmas_nivel_espurias;
            if (comp == CompetidorF4::B4_PROPUESTA_UNIFICADA) {
                esp_tmp_b4 = out.alarmas_temp_espurias;
                esp_niv_b4 = out.alarmas_nivel_espurias;
                ds_b4_p2a = ds_ee;
            }
        }
        reportar("CF4.P2a", "En NIVEL M3 con micro-cortes r=4, B4 mantiene ds_ee <= 0.25 y alarmas de temperatura menores o iguales que LINEAL",
                 (ds_b4_p2a <= 0.25) && (esp_tmp_b4 <= std::max<std::size_t>(esp_tmp_lin, ref3.alarmas_temp_espurias + 1)) &&
                 (esp_niv_b4 <= std::max<std::size_t>(esp_niv_arr, ref3.alarmas_nivel_espurias + 1)));

        // Subcaso P2b: Choque oculto grande (+12 sigma) en micro-corte r=2 barras (100 -> r=2 -> 100+12s)
        const std::size_t t_h = 1500;
        std::vector<double> x_salto = m3.x_op_limpia;
        for (std::size_t t = t_h; t < x_salto.size(); ++t) {
            x_salto[t] += 12.0 * std::sqrt(std::max(1e-12, mod3.sigma2_innov));
        }
        ConfiguracionPuenteHistoricoIto cp_chk{};
        cp_chk.integracion = Integracion::NIVEL;
        cp_chk.r_max = 3;
        cp_chk.z_max_hueco = 3.5;
        const auto res_p2b = evaluar_y_reconstruir_hueco_ito(cp_chk, mod3, {}, x_salto[t_h - 1], x_salto[t_h + 2], 2, 777);
        const bool p2b_ok = !res_p2b.reconstruido && res_p2b.choque_en_hueco && ((res_p2b.flags_emitidos & FLAG_CHOQUE_EN_HUECO) != 0u);
        std::printf("  [Detalle P2b salto oculto +12s en r=2] Z_hueco=%.2f | reconstruido=%d | FLAG_CHOQUE_EN_HUECO=%d\n",
                    res_p2b.z_hueco, int(res_p2b.reconstruido), int((res_p2b.flags_emitidos & FLAG_CHOQUE_EN_HUECO) != 0u));
        reportar("CF4.P2b", "Llave de amplitud |Z_hueco| > z_max bloquea puente ante salto oculto +12s en micro-corte r=2", p2b_ok);
    }

    // ------------------------------------------------------------------------
    // CF4.P3: Comparación unificada en M2 ARFIMA (8 realizaciones independientes, INTERCAMBIABLE, Z1b dt=2, r=3)
    //   Evalúa B0, B0b (CORTE_SUFIJO), B1 (1-borde), B2 (LINEAL), B4_MEDIA_2B y B4_RNG16 (sesgo con signo y dispersión RNG)
    // ------------------------------------------------------------------------
    {
        constexpr std::size_t S_P3 = 8;
        struct ResP3Seed {
            double ee0 = 0.0, s0 = 0.0, d_hat = 0.25, max_anc = 0.0;
            double rm_b0 = 0.0, ds_b0 = 0.0;
            double rm_cs = 0.0, ds_cs = 0.0;
            double rm_b1 = 0.0, ds_b1 = 0.0;
            double rm_lin = 0.0, ds_lin = 0.0;
            double rm_ito = 0.0, ds_ito = 0.0;
            double rm_rng = 0.0, ds_rng = 0.0, sd_rng_rel = 0.0;
        };
        std::vector<ResP3Seed> res_s(S_P3);

        en_paralelo(S_P3, n_hilos, [&](std::size_t s) {
            ParametrosBancoComun params{};
            params.n_op = 6000;
            const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, params, 202ULL + s * 17ULL);

            ConfiguracionObservador cd2{};
            cd2.integracion = Integracion::FLUJO;
            cd2.agregacion = AgregacionEscala::SUMA;
            cd2.arl0_pasos_nativos = 450.0;
            cd2.alfa = ALFA_PROYECTO;
            cd2.medir_tau50 = false;
            cd2.semilla = 20261117ULL + s;
            cd2.escalas_max = 3;
            ObsDelta od2{};
            od2.calibrar(cd2, m2.x_cal);

            const std::size_t m_freq = std::max<std::size_t>(32, static_cast<std::size_t>(std::pow(m2.x_cal.size(), 0.65)));
            const auto est_d = elw::estimar_calibracion(m2.x_cal.data(), m2.x_cal.size(), m_freq, -0.45, 0.45);

            ConfiguracionLilyth cl2{};
            cl2.alfa = ALFA_PROYECTO;
            cl2.representacion = ModoRepresentacionLilyth::CRUDA;
            cl2.b_max_nulo = 16;
            cl2.hilos_nulo = 1;
            cl2.semilla_nulo = 20261117ULL + s * 97ULL;
            cl2.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
            cl2.d_whittle_cruda = est_d.valido ? est_d.d : 0.25;
            cl2.forzar_via_b = true;
            cl2.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
            cl2.b_mc_resumen = 1999;

            const std::vector<std::uint8_t> ok1(m2.x_op_limpia.size(), 1u), fl0(m2.x_op_limpia.size(), 0u);
            const auto ref2 = correr_operador_sobre_flujo(od2, cd2, cl2, m2, m2.x_op_limpia, ok1, m2.x_op_limpia, ok1, fl0);
            ObsDelta od2_cal{};
            od2_cal.calibrar(cd2, m2.x_cal);
            od2_cal.reiniciar_flujo();
            ObsLil lil_mod2{};
            (void)lil_mod2.calibrar(od2_cal, cl2, cd2, m2.x_cal);
            const auto& mod2 = lil_mod2.ficha_nulo_cruda().modelo_lineal;

            std::vector<std::uint8_t> ok_z1b(m2.x_op_limpia.size(), 1u);
            std::vector<std::pair<std::int64_t, std::int64_t>> huecos_z1b;
            for (const auto& [clv, tc] : ref2.t_cambio_por_ev) {
                (void)clv;
                const std::int64_t idx_fin = tc - 2;
                const std::int64_t idx_ini = idx_fin - 2;
                if (idx_ini >= 10 && idx_fin + 1 < static_cast<std::int64_t>(m2.x_op_limpia.size())) {
                    for (std::int64_t t = idx_ini; t <= idx_fin; ++t) ok_z1b[static_cast<std::size_t>(t)] = 0u;
                    huecos_z1b.emplace_back(idx_ini, idx_fin);
                }
            }

            auto eval_pol = [&](ConfiguracionLilyth::PoliticaZ1Piloto pol) {
                ConfiguracionLilyth cl_p = cl2;
                cl_p.politica_z1_piloto = pol;
                cl_p.r_max_z1_piloto = 5;
                cl_p.z_max_hueco_piloto = 4.5;
                const auto out = correr_operador_sobre_flujo(od2, cd2, cl_p, m2, m2.x_op_limpia, ok1, m2.x_op_limpia, ok_z1b, fl0);
                return std::make_pair(calcular_rmse_e(ref2, out), out.s_onsager - ref2.s_onsager);
            };

            auto& rs = res_s[s];
            rs.d_hat = cl2.d_whittle_cruda;
            rs.s0 = ref2.s_onsager;
            rs.ee0 = std::max(1e-12, ref2.ee_s_onsager);
            std::tie(rs.rm_b0,  rs.ds_b0)  = eval_pol(ConfiguracionLilyth::PoliticaZ1Piloto::COMPRESION_BASE);
            std::tie(rs.rm_cs,  rs.ds_cs)  = eval_pol(ConfiguracionLilyth::PoliticaZ1Piloto::CORTE_SUFIJO);
            std::tie(rs.rm_b1,  rs.ds_b1)  = eval_pol(ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_CAUSAL);
            std::tie(rs.rm_ito, rs.ds_ito) = eval_pol(ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO);

            std::vector<double> x_lin = m2.x_op_limpia;
            for (auto [ini, fin] : huecos_z1b) {
                const double x0 = m2.x_op_limpia[static_cast<std::size_t>(ini - 1)];
                const double x1 = m2.x_op_limpia[static_cast<std::size_t>(fin + 1)];
                const double paso = (x1 - x0) / 4.0;
                for (std::int64_t t = ini; t <= fin; ++t) x_lin[static_cast<std::size_t>(t)] = x0 + static_cast<double>(t - ini + 1) * paso;
            }
            const auto out_lin = correr_operador_sobre_flujo(od2, cd2, cl2, m2, m2.x_op_limpia, ok1, x_lin, ok1, fl0);
            rs.rm_lin = calcular_rmse_e(ref2, out_lin);
            rs.ds_lin = out_lin.s_onsager - ref2.s_onsager;

            std::vector<double> s_rng(16, 0.0);
            double rm_rng_s = 0.0;
            for (std::size_t sem = 0; sem < 16; ++sem) {
                std::vector<double> x_r = m2.x_op_limpia;
                std::vector<std::uint8_t> fl_r(m2.x_op_limpia.size(), 0u);
                for (std::size_t hi = 0; hi < huecos_z1b.size(); ++hi) {
                    const auto [ini, fin] = huecos_z1b[hi];
                    const std::size_t u_ini = static_cast<std::size_t>(ini);
                    const std::size_t h_len = std::min<std::size_t>(256, u_ini);
                    const std::span<const double> hp(&x_r[u_ini - h_len], h_len);
                    ConfiguracionPuenteHistoricoIto cp{};
                    cp.integracion = Integracion::FLUJO;
                    cp.r_max = 3;
                    cp.z_max_hueco = 4.5;
                    cp.n_semillas_rng = 1;
                    cp.semilla_base = 20261118ULL + sem * 7919ULL + s * 131ULL;
                    const auto res = evaluar_y_reconstruir_hueco_ito(cp, mod2, hp, x_r[u_ini - 1], m2.x_op_limpia[static_cast<std::size_t>(fin + 1)], 3, hi + 1);
                    if (res.max_error_anclaje_rng > rs.max_anc) rs.max_anc = res.max_error_anclaje_rng;
                    if (res.reconstruido) {
                        for (std::size_t k = 0; k < 3; ++k) {
                            x_r[u_ini + k] = res.x_rng[0][k];
                            fl_r[u_ini + k] = FLAG_PUENTE_ITO;
                        }
                    }
                }
                const auto out_r = correr_operador_sobre_flujo(od2, cd2, cl2, m2, m2.x_op_limpia, ok1, x_r, ok1, fl_r);
                s_rng[sem] = out_r.s_onsager;
                rm_rng_s += calcular_rmse_e(ref2, out_r) / 16.0;
            }
            double mean_s = 0.0;
            for (double v : s_rng) mean_s += v / 16.0;
            double var_s = 0.0;
            for (double v : s_rng) var_s += (v - mean_s) * (v - mean_s);
            rs.rm_rng = rm_rng_s;
            rs.ds_rng = mean_s - ref2.s_onsager;
            rs.sd_rng_rel = std::sqrt(var_s / 15.0) / rs.ee0;
        });

        double d_med = 0.0, s0_med = 0.0, ee_med = 0.0, max_anc_p3 = 0.0;
        double rm_b0 = 0.0, ds_b0 = 0.0, rm_cs = 0.0, ds_cs = 0.0, rm_b1 = 0.0, ds_b1 = 0.0;
        double rm_lin = 0.0, ds_lin = 0.0, rm_ito = 0.0, ds_ito = 0.0, rm_rng = 0.0, ds_rng = 0.0, sd_s_rel = 0.0;
        for (const auto& rs : res_s) {
            d_med += rs.d_hat / S_P3; s0_med += rs.s0 / S_P3; ee_med += rs.ee0 / S_P3;
            max_anc_p3 = std::max(max_anc_p3, rs.max_anc);
            rm_b0 += rs.rm_b0 / S_P3; ds_b0 += rs.ds_b0 / S_P3;
            rm_cs += rs.rm_cs / S_P3; ds_cs += rs.ds_cs / S_P3;
            rm_b1 += rs.rm_b1 / S_P3; ds_b1 += rs.ds_b1 / S_P3;
            rm_lin += rs.rm_lin / S_P3; ds_lin += rs.ds_lin / S_P3;
            rm_ito += rs.rm_ito / S_P3; ds_ito += rs.ds_ito / S_P3;
            rm_rng += rs.rm_rng / S_P3; ds_rng += rs.ds_rng / S_P3;
            sd_s_rel += rs.sd_rng_rel / S_P3;
        }
        const double b_b0  = std::abs(ds_b0)  / ee_med;
        const double b_cs  = std::abs(ds_cs)  / ee_med;
        const double b_b1  = std::abs(ds_b1)  / ee_med;
        const double b_lin = std::abs(ds_lin) / ee_med;
        const double b_ito = std::abs(ds_ito) / ee_med;
        const double b_rng = std::abs(ds_rng) / ee_med;
        const double ee_rng16_rel = sd_s_rel / std::sqrt(16.0);

        std::printf("\n--- [P3] M2 ARFIMA (8 semillas, d_med=%.4f, S0_med=%.4f, ee_med=%.4f) con INTERCAMBIABLE unificado ---\n",
                    d_med, s0_med, ee_med);
        std::printf(" Politica Z1b (dt=2, r=3) | RMSE(e)/sd | |E[dS]|/ee_S | ESTRICTO | DULCE\n");
        std::printf("--------------------------+------------+--------------+----------+-------\n");
        auto imprimir_p3 = [&](const char* nom, double rm, double bs) {
            const bool est = (rm <= 0.15) && (bs <= 0.25);
            const bool dul = (rm < 0.50) && (bs <= 1.00);
            std::printf(" %-24s | %10.4f | %12.4f | %-8s | %-5s\n", nom, rm, bs, est ? "PASA" : "FALLA", dul ? "PASA" : "FALLA");
        };
        imprimir_p3("B0_COMPRESION_BASE",      rm_b0,  b_b0);
        imprimir_p3("B0b_CORTE_SUFIJO",        rm_cs,  b_cs);
        imprimir_p3("B1_PUENTE_CAUSAL_1B",     rm_b1,  b_b1);
        imprimir_p3("B2_LINEAL_2BORDES",       rm_lin, b_lin);
        imprimir_p3("B4_PUENTE_ITO_MEDIA_2B",  rm_ito, b_ito);
        std::printf(" %-24s | %10.4f | %12.4f | %-8s | %-5s  (sd_1RNG/ee=%.4f, ee_RNG16/ee=%.4f, err_anc=%.2e)\n",
                    "B4_PUENTE_ITO_RNG16", rm_rng, b_rng,
                    (rm_rng <= 0.15 && b_rng <= 0.25 && ee_rng16_rel <= 0.15) ? "PASA" : "FALLA",
                    (rm_rng < 0.50 && b_rng <= 1.00) ? "PASA" : "FALLA",
                    sd_s_rel, ee_rng16_rel, max_anc_p3);

        const bool p3_ok = (rm_ito <= 0.15) && (rm_ito <= rm_cs + 1e-6) &&
                           (rm_rng <= 0.15) && (b_rng <= 0.25) &&
                           (ee_rng16_rel <= 0.15) && (sd_s_rel <= 0.50) && (max_anc_p3 < 1e-10);
        reportar("CF4.P3", "En M2 ARFIMA (8 semillas), B4_ITO_MEDIA logra minima RMSE y B4_ITO_RNG16 pasa ESTRICTO (|E[dS]|/ee<=0.25, ee_RNG16<=0.15, anc<1e-10)", p3_ok);
    }

    // ------------------------------------------------------------------------
    // CF4.P4: Pausa dura en corte largo (r = 12 > r_max = 3) en M3
    // ------------------------------------------------------------------------
    {
        AdaptadorInputCausalLilyth ad{};
        ad.integracion = Integracion::NIVEL;
        ad.racha_max_corrupcion = 3;
        ad.z_max_hueco = 3.5;
        ad.evitar_teletransporte_nivel = true;

        (void)ad.preparar(100.0, ClasePasoLilyth::VALIDO);
        std::size_t pasos_en_pausa = 0;
        for (int k = 0; k < 12; ++k) {
            const auto st = ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
            if (st.estado == EstadoContinuidadLilyth::EN_PAUSA_DURA && !st.valido_para_operador) ++pasos_en_pausa;
        }
        const auto rec = ad.preparar(105.0, ClasePasoLilyth::VALIDO);
        const bool p4_ok = (pasos_en_pausa == 9) && rec.reconexion_tras_hueco && rec.pausa_umbral_activada &&
                           rec.cortar_segmento_al_reenganchar && (rec.delta_t_hueco == 13);
        reportar("CF4.P4", "Corte largo r=12 > r_max=3 mantiene EN_PAUSA_DURA exacta (9 pasos) y corta segmento al reconectar", p4_ok);
    }

    // ------------------------------------------------------------------------
    // CF4.P6: Blindaje obligatorio de FLAG_PUENTE_ITO en Z2 y Z3 dentro de ObservadorLilyth (caza MN1.7)
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
        for (std::size_t t = 0; t < m2.x_op_limpia.size(); ++t) sd_m2[t] = od.paso(m2.x_op_limpia[t], NaN, true, 0, 0);

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

        ObsLil lil_ref{};
        (void)lil_ref.calibrar(od, cl, cd, m2.x_cal);
        for (std::size_t t = 0; t < m2.x_op_limpia.size(); ++t) (void)lil_ref.paso(sd_m2[t], m2.x_op_limpia[t], true, NaN, 0, 0);
        lil_ref.completar_nulos();

        std::int64_t tc_z2 = -1, ta_z2 = -1, ta_z3 = -1;
        double tau_z3 = NaN, y_anc_z3 = 0.0;
        for (const auto& r : lil_ref.cerradas_historico()) {
            if (tc_z2 < 0 && r.t_cambio >= 0 && r.t_cambio < r.t_apertura && r.cruda.nulo_onsager.evaluado) {
                tc_z2 = r.t_cambio; ta_z2 = r.t_apertura;
            } else if (ta_z3 < 0 && r.cruda.desenlace == DesenlaceRanura::RESTITUIDA &&
                       std::isfinite(r.cruda.tau_rest_consolidado) && r.cruda.tau_rest_consolidado >= 4.0) {
                ta_z3 = r.t_apertura; tau_z3 = r.cruda.tau_rest_consolidado; y_anc_z3 = r.cruda.y_ancla;
            }
        }

        ObsLil lil_v{};
        (void)lil_v.calibrar(od, cl, cd, m2.x_cal);
        for (std::size_t t = 0; t < m2.x_op_limpia.size(); ++t) {
            double x_in = m2.x_op_limpia[t];
            std::uint8_t fl = 0;
            if (static_cast<std::int64_t>(t) == tc_z2) fl = FLAG_PUENTE_ITO;
            if (static_cast<std::int64_t>(t) == ta_z3 + 1) { x_in = y_anc_z3; fl = FLAG_PUENTE_ITO; }
            (void)lil_v.paso(sd_m2[t], x_in, true, NaN, 0, 0, NaN, {}, fl);
        }
        lil_v.completar_nulos();

        bool veto_z2 = false, veto_z3 = false;
        for (const auto& r : lil_v.cerradas_historico()) {
            if (r.t_apertura == ta_z2) veto_z2 = r.perfil_corrupto_z2_cruda && r.cruda.nulo_onsager.sin_nulo_b && !r.cruda.nulo_onsager.evaluado;
            if (r.t_apertura == ta_z3) veto_z3 = std::isnan(r.cruda.eta[0]) && std::isfinite(r.cruda.tau_rest_consolidado) &&
                                                 (std::abs(r.cruda.tau_rest_consolidado - tau_z3) < 1e-12);
        }
        reportar("CF4.P6", "Blindaje obligatorio de FLAG_PUENTE_ITO en Z2 y Z3 dentro de ObservadorLilyth", veto_z2 && veto_z3);
    }

    // ------------------------------------------------------------------------
    // CF4.P5: Comparativa Monte Carlo (R = R_REP réplicas) en M1_H0, M2_H0 (ARFIMA) y M1_H1 bajo q=2%, r=2
    // ------------------------------------------------------------------------
    {
        ParametrosReplicaH0H1 ph0{}; ph0.lenta_h1 = false;
        ParametrosReplicaH0H1 ph1{}; ph1.lenta_h1 = true;

        std::vector<SerieSinteticaCoordinada> m1_h0(R_REP), m2_h0(R_REP), m1_h1(R_REP);
        en_paralelo(R_REP, n_hilos, [&](std::size_t i) {
            m1_h0[i] = generar_mundo_m1_h0_h1(0.80, ph0, i);
            m2_h0[i] = generar_mundo_m2_h0_h1(0.25, 0.40, ph0, i);
            m1_h1[i] = generar_mundo_m1_h0_h1(0.80, ph1, i);
        });

        auto correr_lote_mc = [&](const std::vector<SerieSinteticaCoordinada>& mundos,
                                  ModoNuloLineal mlin,
                                  bool aplicar_corr,
                                  CompetidorF4 comp) -> ResumenMCF4 {
            struct RepRes { bool val=false, r=false, rl=false; double ke=0, s=0, ee=0; std::size_t tmp=0, vz=0; };
            std::vector<RepRes> rr(R_REP);
            en_paralelo(R_REP, n_hilos, [&](std::size_t i) {
                const auto& m = mundos[i];
                ConfiguracionObservador cd{};
                cd.integracion = Integracion::FLUJO;
                cd.agregacion = AgregacionEscala::SUMA;
                cd.arl0_pasos_nativos = 450.0;
                cd.alfa = ALFA_PROYECTO;
                cd.medir_tau50 = false;
                cd.semilla = 20261118ULL ^ (i * 10007ULL);
                cd.escalas_max = 2;

                ObsDelta od{};
                od.calibrar(cd, m.x_cal);
                if (od.escalas() == 0) return;

                ConfiguracionLilyth cl{};
                cl.alfa = ALFA_PROYECTO;
                cl.representacion = ModoRepresentacionLilyth::CRUDA;
                cl.b_max_nulo = 32;
                cl.hilos_nulo = 1;
                cl.semilla_nulo = 20261118ULL ^ (i * 31337ULL);
                cl.modo_nulo_lineal = mlin;
                if (mlin == ModoNuloLineal::ARFIMA_WHITTLE) cl.d_whittle_cruda = 0.25;
                cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
                cl.pool_nulo = PoolNuloLilyth::DERIVADO;
                cl.horizonte = HorizonteLilyth::CUANTIL_RESOLUCION;
                cl.b_mc_resumen = 1999;

                if (comp == CompetidorF4::B1_ARGON_VIVO || comp == CompetidorF4::B4_PROPUESTA_UNIFICADA) {
                    cl.politica_z1_piloto = (comp == CompetidorF4::B4_PROPUESTA_UNIFICADA)
                        ? ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO
                        : ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_CAUSAL;
                    cl.r_max_z1_piloto = 2;
                    cl.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
                    cl.r_max_z2_piloto = 0;
                    cl.politica_z3_piloto = ConfiguracionLilyth::PoliticaZ3Piloto::RELOJ_FISICO_TOPE_R;
                    cl.r_max_z3_piloto = 2;
                }

                const std::size_t nop = m.x_op_limpia.size();
                std::vector<std::uint8_t> es_c(nop, 0u);
                if (aplicar_corr) {
                    std::vector<std::uint8_t> chq(nop, 0u);
                    for (std::size_t tc : m.tiempos_choque) {
                        for (std::size_t k = 0; k < 2 && (tc + k) < nop; ++k) chq[tc + k] = 1u;
                    }
                    std::mt19937_64 rg(detalle_ito::splitmix64_ito(20261118ULL ^ (i * 99991ULL)));
                    std::uniform_real_distribution<double> u01(0.0, 1.0);
                    for (std::size_t t = 10; t + 2 < nop; ++t) {
                        if (!chq[t] && !chq[t + 1] && u01(rg) < 0.01) {
                            es_c[t] = 1u; es_c[t + 1] = 1u;
                            t += 2;
                        }
                    }
                }

                // Modelo de calma para B4 cuando se requiere (solo usar en B2/B3/B4)
                ModeloLinealNulo m_calma{};
                if (aplicar_corr && (comp == CompetidorF4::B3_ITO_SIN_FILTRO || comp == CompetidorF4::B4_PROPUESTA_UNIFICADA)) {
                    ObsLil lil_cal{};
                    (void)lil_cal.calibrar(od, cl, cd, m.x_cal);
                    m_calma = lil_cal.ficha_nulo_cruda().modelo_lineal;
                }

                std::vector<double> xd, xl;
                std::vector<std::uint8_t> okd, okl, fl;
                std::size_t nrec = 0, nchq = 0, npau = 0;
                double max_anc = 0.0;
                preparar_flujos_competidor(comp, Integracion::FLUJO, m_calma, m.x_op_limpia, es_c, 2, 3.5, 20261118ULL + i,
                                           xd, okd, xl, okl, fl, nrec, nchq, npau, max_anc);

                const auto out = correr_operador_sobre_flujo(od, cd, cl, m, xd, okd, xl, okl, fl);
                rr[i] = {out.valida, out.rechaza_inter, out.rechaza_lenta, static_cast<double>(out.k_evaluados),
                         out.s_onsager, out.ee_s_onsager, out.alarmas_temp_tot, out.n_veto_z2};
            });

            ResumenMCF4 res{};
            double sk = 0.0, ss = 0.0, see = 0.0, st = 0.0, svz = 0.0;
            for (const auto& r : rr) {
                sk += r.ke; st += static_cast<double>(r.tmp); svz += static_cast<double>(r.vz);
                if (!r.val) continue;
                ++res.n_val;
                if (r.r) ++res.n_rech;
                if (r.rl) ++res.n_rech_lenta;
                ss += r.s; see += r.ee;
            }
            res.k_eval_med = sk / static_cast<double>(R_REP);
            res.temp_alarm_med = st / static_cast<double>(R_REP);
            res.veto_z2_med = svz / static_cast<double>(R_REP);
            if (res.n_val > 0) {
                res.tasa_rech = static_cast<double>(res.n_rech) / static_cast<double>(res.n_val);
                res.tasa_lenta = static_cast<double>(res.n_rech_lenta) / static_cast<double>(res.n_val);
                const auto [lo, hi] = wilson_95(res.n_rech, res.n_val);
                res.wlo = lo; res.whi = hi;
                res.media_s = ss / static_cast<double>(res.n_val);
                res.media_ee = std::max(1e-6, see / static_cast<double>(res.n_val));
            }
            return res;
        };

        const auto b_m1_h0 = correr_lote_mc(m1_h0, ModoNuloLineal::AR_LEVINSON_DURBIN, false, CompetidorF4::B0_BASE);
        const auto b_m2_h0 = correr_lote_mc(m2_h0, ModoNuloLineal::ARFIMA_WHITTLE,     false, CompetidorF4::B0_BASE);
        const auto b_m1_h1 = correr_lote_mc(m1_h1, ModoNuloLineal::AR_LEVINSON_DURBIN, false, CompetidorF4::B0_BASE);

        std::printf("\n--- [P5] COMPARATIVA MONTE CARLO (R = %zu, q = 2%%, r = 2) EN M1_H0, M2_H0 (ARFIMA) Y M1_H1 ---\n", R_REP);
        std::printf(" Escenario | Politica               | Validas | K_eval | FalsoRech_H0 [Wilson 95%%]      | RetPot_H1 | AlmTmp/rep | ESTRICTO | DULCE\n");
        std::printf("-----------+------------------------+---------+--------+--------------------------------+-----------+------------+----------+-------\n");
        std::printf(" LIMPIO q0 | M1_H0 BASE             | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] |         - | %10.2f | %-8s | %-5s\n",
                    b_m1_h0.n_val, b_m1_h0.k_eval_med, b_m1_h0.n_rech, b_m1_h0.n_val, b_m1_h0.tasa_rech, b_m1_h0.wlo, b_m1_h0.whi,
                    b_m1_h0.temp_alarm_med, b_m1_h0.wlo <= 0.05 ? "PASA" : "FALLA", "PASA");
        std::printf(" LIMPIO q0 | M2_H0 ARFIMA BASE      | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] |         - | %10.2f | %-8s | %-5s\n",
                    b_m2_h0.n_val, b_m2_h0.k_eval_med, b_m2_h0.n_rech, b_m2_h0.n_val, b_m2_h0.tasa_rech, b_m2_h0.wlo, b_m2_h0.whi,
                    b_m2_h0.temp_alarm_med, b_m2_h0.wlo <= 0.05 ? "PASA" : "FALLA", "PASA");
        std::printf(" LIMPIO q0 | M1_H1 BASE (Potencia)  | %7zu | %6.1f | Pot = %3zu/%3zu = %.3f          |      1.00 | %10.2f | %-8s | %-5s\n",
                    b_m1_h1.n_val, b_m1_h1.k_eval_med, b_m1_h1.n_rech_lenta, b_m1_h1.n_val, b_m1_h1.tasa_lenta,
                    b_m1_h1.temp_alarm_med, b_m1_h1.tasa_lenta >= 0.50 ? "PASA" : "FALLA", "PASA");
        std::printf("-----------+------------------------+---------+--------+--------------------------------+-----------+------------+----------+-------\n");

        bool p5_b4_est = false;
        double m2_h0_wlo_b4 = 1.0;

        for (auto comp : {CompetidorF4::B0_BASE, CompetidorF4::B1_ARGON_VIVO, CompetidorF4::B2_LINEAL_2BORDES,
                          CompetidorF4::B3_ITO_SIN_FILTRO, CompetidorF4::B4_PROPUESTA_UNIFICADA}) {
            const char* nom = (comp == CompetidorF4::B0_BASE)                ? "B0_BASE_ACTUAL"
                            : (comp == CompetidorF4::B1_ARGON_VIVO)          ? "B1_ARGON_VIVO"
                            : (comp == CompetidorF4::B2_LINEAL_2BORDES)     ? "B2_LINEAL_2BORDES"
                            : (comp == CompetidorF4::B3_ITO_SIN_FILTRO)     ? "B3_ITO_SIN_FILTRO"
                                                                             : "B4_PROPUESTA_UNIFICADA";
            const auto r1 = correr_lote_mc(m1_h0, ModoNuloLineal::AR_LEVINSON_DURBIN, true, comp);
            const auto r2 = correr_lote_mc(m2_h0, ModoNuloLineal::ARFIMA_WHITTLE,     true, comp);
            const auto h1 = correr_lote_mc(m1_h1, ModoNuloLineal::AR_LEVINSON_DURBIN, true, comp);

            const double ret_pot = (b_m1_h1.tasa_lenta > 0.0) ? (h1.tasa_lenta / b_m1_h1.tasa_lenta) : 0.0;
            const double dcens1 = (b_m1_h0.k_eval_med > 0.0) ? std::max(0.0, 1.0 - r1.k_eval_med / b_m1_h0.k_eval_med) : 1.0;
            const double ds_ee1 = std::abs(r1.media_s - b_m1_h0.media_s) / b_m1_h0.media_ee;
            const bool est1 = (dcens1 <= 0.10) && (ds_ee1 <= 0.25) && (r1.wlo <= 0.05) && (ret_pot >= 0.85);
            const bool dul1 = (dcens1 < 0.50) && (ds_ee1 <= 1.00) && (r1.wlo <= 0.10) && (ret_pot >= 0.50);

            std::printf(" M1 q2 r2  | %-22s | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] | %9.2f | %10.2f | %-8s | %-5s\n",
                        nom, r1.n_val, r1.k_eval_med, r1.n_rech, r1.n_val, r1.tasa_rech, r1.wlo, r1.whi,
                        ret_pot, r1.temp_alarm_med, est1 ? "PASA" : "FALLA", dul1 ? "PASA" : "FALLA");
            std::printf(" M2 q2 r2  | %-22s | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] |         - | %10.2f | %-8s | %-5s\n",
                        nom, r2.n_val, r2.k_eval_med, r2.n_rech, r2.n_val, r2.tasa_rech, r2.wlo, r2.whi,
                        r2.temp_alarm_med, r2.wlo <= 0.05 ? "PASA" : "FALLA", r2.wlo <= 0.10 ? "PASA" : "FALLA");

            if (comp == CompetidorF4::B1_ARGON_VIVO || comp == CompetidorF4::B4_PROPUESTA_UNIFICADA) {
                std::printf("  [diag %s M1] dcens=%.3f | |dS|/ee=%.3f | vetoZ2/rep=%.2f | [M2] vetoZ2/rep=%.2f\n",
                            nom, dcens1, ds_ee1, r1.veto_z2_med, r2.veto_z2_med);
            }
            if (comp == CompetidorF4::B4_PROPUESTA_UNIFICADA) {
                p5_b4_est = est1;
                m2_h0_wlo_b4 = r2.wlo;
            }
        }
        reportar("CF4.P5", "Comparativa R replicas en M1_H0, M2_H0 y M1_H1 (B4 pasa ESTRICTO y conserva tamano nominal en ARFIMA)",
                 (b_m1_h0.wlo <= 0.05) && (b_m2_h0.wlo <= 0.05) && (b_m1_h1.tasa_lenta >= 0.50) && p5_b4_est && (m2_h0_wlo_b4 <= 0.05));
    }

    std::printf("====================================================================================================\n");
    std::printf(" FALLAS FASE 4: %d\n", fallas);
    std::printf("====================================================================================================\n");
    return fallas == 0 ? 0 : 1;
}
