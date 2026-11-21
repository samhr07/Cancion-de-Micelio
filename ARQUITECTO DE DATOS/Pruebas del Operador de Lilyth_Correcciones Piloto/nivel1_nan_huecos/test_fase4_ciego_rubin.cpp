// ============================================================================
// test_fase4_ciego_rubin.cpp -- Batería Ciega Sellada de Fase 4 (Fusión Limpia)
// Contrato sellado en: nivel1_nan_huecos/EXPECTATIVAS_FASE4.md
// Sello SHA-256 en:    nivel1_nan_huecos/SELLO_EXPECTATIVAS_FASE4.txt
//
// Controles:
//   CF4.P1  : Identidad bit a bit (%a) real contra volcado_orig.txt
//   CF4.P2a : NIVEL M3 (semilla ciega) micro-cortes r=4, B0..B4_RUBIN (M=16)
//   CF4.P2b : Matriz 4x5 choque oculto en NIVEL M3 con historia causal (K=20)
//   CF4.P2c : Varianza intra-hueco real (B2 vs B4) + s^2 de temperatura + choques +10s (B3 vs B4)
//   CF4.P3  : M2 ARFIMA (16 semillas ciegas, forzar_via_b=true, Z1b dt=2 r=3, Rubin M=16 por serie)
//   CF4.P4  : Máquina de pausa dura EN_PAUSA_DURA en r=12 > r_max=3
//   CF4.P5  : Bloque ciego Monte Carlo R=200 (M1_H0, M2_H0, M3_H0, M1_H1), B0..B4_RUBIN (M=5 + Holm)
//   CF4.P6  : Blindaje obligatorio de FLAG_PUENTE_ITO en Z2 && Z3 dentro de ObservadorLilyth
//   CF4.P7  : Deriva no nula en FLUJO (mu = +5.0*sigma_e) vía valor_z0_imputado (caza MN1.10)
//
// Uso: test_fase4_ciego_rubin.exe [R_replicas=200] [ruta_volcado_orig]
// ============================================================================

#include "ark/holm.hpp"
#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/lilyth/puente_historico_ito.hpp"
#include "ark/metrologia/observador_saltos_filas.hpp"
#include "banco_sintetico_comun.hpp"
#include "protocolo_inputs_lilyth.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
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

inline constexpr std::size_t T_MAX_F4 = 64;
inline constexpr std::size_t M_RUBIN_FINO = 16;
inline constexpr std::size_t M_RUBIN_MC = 5;
inline constexpr double ALFA_HUECO_F4 = 0.01;
inline constexpr std::uint64_t SEMILLA_CIEGA_BASE = 20261119ULL;

using ObsDelta = ObservadorSaltos<24, 64, 120>;
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
    for (auto& th : trab) th.join();
}

double mediana_vector(std::vector<double> v) {
    if (v.empty()) return NaN;
    std::sort(v.begin(), v.end());
    const std::size_t m = v.size();
    return (m % 2 == 1) ? v[m / 2] : 0.5 * (v[m / 2 - 1] + v[m / 2]);
}

double varianza_muestral(const double* x, std::size_t n) noexcept {
    if (n < 2) return 0.0;
    double mean = 0.0;
    for (std::size_t i = 0; i < n; ++i) mean += x[i];
    mean /= static_cast<double>(n);
    double ss = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double d = x[i] - mean;
        ss += d * d;
    }
    return ss / static_cast<double>(n - 1);
}

// ----------------------------------------------------------------------------
// Utilidades de volcado bit a bit para CF4.P1 (idénticas a volcado_identidad_nivel1.cpp)
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
// Ejecutor del operador real (ObservadorSaltos + ObservadorLilyth)
// ----------------------------------------------------------------------------
struct SalidaCompletaOperador {
    std::map<ClaveEvento, std::vector<double>> e_iny_por_ev;
    std::map<ClaveEvento, std::int64_t> t_cambio_por_ev;
    std::size_t k_medibles = 0;
    std::size_t k_evaluados = 0;
    std::size_t n_veto_z2 = 0;
    double l_z2_suma = 0.0;
    std::size_t l_z2_cuenta = 0;
    std::size_t alarmas_nivel_tot = 0;
    std::size_t alarmas_temp_tot = 0;
    std::size_t alarmas_nivel_espurias = 0;
    std::size_t alarmas_temp_espurias = 0;
    double s_onsager = NaN;
    double ee_s_onsager = NaN;
    std::array<double, 4> p_crudos{NaN, NaN, NaN, NaN};
    std::array<double, 4> p_holm{NaN, NaN, NaN, NaN};
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
    const std::vector<std::uint8_t>& es_ventana_choque_real = {},
    std::vector<std::pair<std::size_t, double>>* s2_temp_out = nullptr) {

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

        if (s2_temp_out != nullptr && sd.niveles > 0 && sd.temp[0].emitida && std::isfinite(sd.temp[0].log_beta)) {
            s2_temp_out->emplace_back(t, std::exp(-sd.temp[0].log_beta));
        }

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
    out.p_crudos = fa.p_crudos;
    out.p_holm = fa.p_holm;
    out.valida = fn.contraste_onsager.valido && fa.familia_holm_valida && std::isfinite(fa.p_holm[1]);
    out.rechaza_inter = out.valida && (fa.p_holm[1] < ALFA_PROYECTO);
    out.rechaza_lenta = out.rechaza_inter && (out.s_onsager > 0.0);

    for (const auto& r : lil.cerradas_historico()) {
        const auto& rama = usar_purga ? r.purgada : r.cruda;
        if (rama.desenlace == DesenlaceRanura::RESTITUIDA || rama.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) {
            ++out.k_medibles;
        }
        if (rama.nulo_onsager.evaluado) {
            ++out.k_evaluados;
            if (r.t_apertura >= r.t_cambio) {
                out.l_z2_suma += static_cast<double>(r.t_apertura - r.t_cambio + 1);
                ++out.l_z2_cuenta;
            }
        }
        if (usar_purga ? r.perfil_corrupto_z2_purgada : r.perfil_corrupto_z2_cruda) ++out.n_veto_z2;

        const ClaveEvento clv{r.t_apertura, r.escala_origen, r.signo};
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
    double alfa_hueco,
    std::uint64_t semilla_rng,
    std::size_t m_idx_rubin,
    std::size_t n_semillas_total,
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

    const double mu_calma = (mod_calma.valido && std::isfinite(mod_calma.mu)) ? mod_calma.mu : 0.0;

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
                    x_delta_out[k] = NaN; ok_delta_out[k] = 0u;
                } else if (k - t + 1 <= r_max && std::isfinite(prev_x_ok)) {
                    x_delta_out[k] = mu_calma; ok_delta_out[k] = 1u;
                } else {
                    x_delta_out[k] = NaN; ok_delta_out[k] = 0u;
                }
                x_lil_out[k] = NaN; ok_lil_out[k] = 0u;
                flags_out[k] = FLAG_CORRUPCION_CAIDA | (r_gap > r_max ? FLAG_PAUSA_UMBRAL : 0u);
            }
            if (integ == Integracion::NIVEL) prev_x_ok = NaN;
            t = t_fin + 1;
            continue;
        }

        if (comp == CompetidorF4::B1B_ARRASTRE_NIVEL) {
            const double val_arr = std::isfinite(prev_x_ok) ? prev_x_ok : mu_calma;
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

        if (comp == CompetidorF4::B2_LINEAL_2BORDES) {
            if (r_gap <= r_max) {
                const double paso = (x_fin - x_ini) / static_cast<double>(r_gap + 1);
                for (std::size_t k = 0; k < r_gap; ++k) {
                    const double xv = x_ini + static_cast<double>(k + 1) * paso;
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
            t = t_fin + 1;
            continue;
        }

        ConfiguracionPuenteHistoricoIto cp{};
        cp.integracion = integ;
        cp.r_max = (comp == CompetidorF4::B3_ITO_SIN_FILTRO) ? MAX_RACHA_PUENTE_ITO : r_max;
        cp.alfa_hueco = alfa_hueco;
        if (comp == CompetidorF4::B3_ITO_SIN_FILTRO) cp.z_max_hueco = 1e9;
        cp.n_semillas_rng = std::max<std::size_t>(1, n_semillas_total);
        cp.semilla_base = semilla_rng;

        const auto res = evaluar_y_reconstruir_hueco_ito(cp, mod_calma, hist_est, x_ini, x_fin, r_gap, id_h++);
        if (res.max_error_anclaje_rng > max_err_anclaje_out) max_err_anclaje_out = res.max_error_anclaje_rng;

        const std::size_t m_ef = (res.n_semillas_efectivas > 0) ? std::min(m_idx_rubin, res.n_semillas_efectivas - 1) : 0;

        if (comp == CompetidorF4::B3_ITO_SIN_FILTRO) {
            if (res.reconstruido) {
                for (std::size_t k = 0; k < r_gap; ++k) {
                    const double xv = res.x_rng[m_ef][k];
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
            if (res.reconstruido) {
                for (std::size_t k = 0; k < r_gap; ++k) {
                    const double xv = res.x_rng[m_ef][k];
                    const double xd_val = valor_z0_imputado(res, integ, m_ef, k, mu_calma);
                    x_delta_out[t + k] = xd_val;
                    ok_delta_out[t + k] = std::isfinite(xd_val) ? 1u : 0u;
                    x_lil_out[t + k] = xv;
                    ok_lil_out[t + k] = 1u;
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
    double l_z2_med = 1.0;
    double media_s = 0.0;
    double media_ee = 0.05;
    double temp_alarm_med = 0.0;
    double veto_z2_med = 0.0;
    double media_gamma_rubin = 0.0;
};

} // namespace

int main(int argc, char** argv) {
    const std::size_t R_REP = (argc > 1) ? static_cast<std::size_t>(std::max(10, std::atoi(argv[1]))) : 200;
    const std::string ruta_volcado_orig = (argc > 2) ? argv[2] : "nivel1_nan_huecos/volcado_orig.txt";
    const std::size_t n_hilos = std::max<std::size_t>(1, std::thread::hardware_concurrency());

    int fallas = 0;
    auto reportar = [&](const char* id, const char* desc, bool ok) {
        if (!ok) ++fallas;
        std::printf("[%s] %s: %s\n", ok ? "PASA" : "FALLA", id, desc);
    };

    const double z_max_sel = derivar_z_max_hueco(ConfiguracionPuenteHistoricoIto{});
    std::printf("====================================================================================================\n");
    std::printf(" FASE 4 SELLADA (EXPECTATIVAS_FASE4.md): alfa_hueco=%.2f (z_max=%.4f), M_fino=%zu, M_mc=%zu, R=%zu, Hilos=%zu\n",
                ALFA_HUECO_F4, z_max_sel, M_RUBIN_FINO, M_RUBIN_MC, R_REP, n_hilos);
    std::printf("====================================================================================================\n");

    // ------------------------------------------------------------------------
    // CF4.P1: Identidad bit a bit REAL (%a) contra volcado_orig.txt (observador_lilyth_original.hpp)
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
        reportar("CF4.P1", "Identidad bit a bit REAL (%a) contra volcado_orig.txt (observador_lilyth_original.hpp) intacta", id_ok);
    }

    // ------------------------------------------------------------------------
    // CF4.P2a, CF4.P2b y CF4.P2c: NIVEL M3 + PurgaDiurna sobre semilla ciega 20261119ULL
    // ------------------------------------------------------------------------
    {
        ParametrosBancoComun pb{};
        pb.n_op = 6000;
        const auto m3 = generar_mundo_m3_nivel_purga(0.65, pb, SEMILLA_CIEGA_BASE);

        ConfiguracionObservador cd3{};
        cd3.integracion = Integracion::NIVEL;
        cd3.agregacion = AgregacionEscala::SUMA;
        cd3.arl0_pasos_nativos = 450.0;
        cd3.alfa = ALFA_PROYECTO;
        cd3.purga = true;
        cd3.bins_fase = 24;
        cd3.dias_minimos_purga = 20;
        cd3.medir_tau50 = false;
        cd3.semilla = SEMILLA_CIEGA_BASE;
        cd3.escalas_max = 3;

        ObsDelta od3{};
        od3.calibrar(cd3, m3.x_cal, {}, {}, m3.fase_cal, m3.dia_cal);

        ConfiguracionLilyth cl3{};
        cl3.alfa = ALFA_PROYECTO;
        cl3.representacion = ModoRepresentacionLilyth::PURGADA;
        cl3.b_max_nulo = 16;
        cl3.hilos_nulo = 4;
        cl3.semilla_nulo = SEMILLA_CIEGA_BASE;
        cl3.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
        cl3.b_mc_resumen = 1999;

        ObsLil lil_cal3{};
        (void)lil_cal3.calibrar(od3, cl3, cd3, m3.x_cal, {}, m3.fase_cal, m3.dia_cal);
        const ModeloLinealNulo mod3 = lil_cal3.ficha_nulo_purgada().modelo_lineal;
        const double sd_e3 = std::sqrt(std::max(1e-12, mod3.sigma2_innov));

        const std::size_t nop = m3.x_op_limpia.size();
        std::vector<std::uint8_t> en_choque3(nop, 0u);
        for (std::size_t tc : m3.tiempos_choque) {
            const std::size_t ini = (tc >= 5) ? (tc - 5) : 0;
            const std::size_t fin = std::min(nop, tc + 55);
            for (std::size_t t = ini; t < fin; ++t) en_choque3[t] = 1u;
        }
        const std::vector<std::uint8_t> ok1(nop, 1u), fl0(nop, 0u);
        const auto ref3 = correr_operador_sobre_flujo(od3, cd3, cl3, m3, m3.x_op_limpia, ok1, m3.x_op_limpia, ok1, fl0, en_choque3);

        // --- CF4.P2a: micro-cortes limpios r=4 fuera de choques reales; tabla completa B0..B4_RUBIN (M=16) ---
        std::vector<std::uint8_t> corr_p2a(nop, 0u);
        for (std::size_t tc : m3.tiempos_choque) {
            if (tc > 180) {
                bool libre = true;
                for (std::size_t k = 0; k < 4; ++k) if (en_choque3[tc - 160 + k]) libre = false;
                if (libre) for (std::size_t k = 0; k < 4; ++k) corr_p2a[tc - 160 + k] = 1u;
            }
        }

        std::printf("\n--- [P2a] NIVEL M3 + PurgaDiurna (semilla ciega %llu, micro-cortes r=4, K_ref=%zu, S0=%.4f) ---\n",
                    static_cast<unsigned long long>(SEMILLA_CIEGA_BASE), ref3.k_evaluados, ref3.s_onsager);
        std::printf(" Competidor              | Alm_Niv_Esp | Alm_Tmp_Esp | K_eval | |dS|/ee | err_Anc_RNG | Notas\n");
        std::printf("-------------------------+-------------+-------------+--------+---------+-------------+-----------------------------\n");

        std::size_t esp_tmp_b2 = 0;
        for (auto comp : {CompetidorF4::B0_BASE, CompetidorF4::B1_ARGON_VIVO, CompetidorF4::B1B_ARRASTRE_NIVEL,
                          CompetidorF4::B2_LINEAL_2BORDES, CompetidorF4::B3_ITO_SIN_FILTRO}) {
            std::vector<double> xd, xl;
            std::vector<std::uint8_t> okd, okl, fl;
            std::size_t nrec = 0, nchq = 0, npau = 0;
            double ma = 0.0;
            preparar_flujos_competidor(comp, Integracion::NIVEL, mod3, m3.x_op_limpia, corr_p2a, 4, ALFA_HUECO_F4,
                                       SEMILLA_CIEGA_BASE, 0, M_RUBIN_FINO, xd, okd, xl, okl, fl, nrec, nchq, npau, ma);
            ConfiguracionLilyth cl_c = cl3;
            if (comp == CompetidorF4::B1_ARGON_VIVO) {
                cl_c.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
                cl_c.politica_z3_piloto = ConfiguracionLilyth::PoliticaZ3Piloto::RELOJ_FISICO_TOPE_R;
                cl_c.r_max_z3_piloto = 4;
            }
            const auto out = correr_operador_sobre_flujo(od3, cd3, cl_c, m3, xd, okd, xl, okl, fl, en_choque3);
            const double ds_ee = std::abs(out.s_onsager - ref3.s_onsager) / std::max(1e-12, ref3.ee_s_onsager);
            const char* nom = (comp == CompetidorF4::B0_BASE)             ? "B0_BASE_ACTUAL"
                            : (comp == CompetidorF4::B1_ARGON_VIVO)       ? "B1_ARGON_VIVO"
                            : (comp == CompetidorF4::B1B_ARRASTRE_NIVEL)  ? "B1b_ARRASTRE_NIVEL"
                            : (comp == CompetidorF4::B2_LINEAL_2BORDES)   ? "B2_LINEAL_2BORDES"
                                                                          : "B3_ITO_SIN_FILTRO";
            std::printf(" %-23s | %11zu | %11zu | %6zu | %7.4f | %11.2e | rec=%zu chq=%zu pau=%zu\n",
                        nom, out.alarmas_nivel_espurias, out.alarmas_temp_espurias, out.k_evaluados, ds_ee, ma, nrec, nchq, npau);
            if (comp == CompetidorF4::B2_LINEAL_2BORDES) esp_tmp_b2 = out.alarmas_temp_espurias;
        }

        // B4_RUBIN16: M=16 imputaciones reales combinadas mediante reglas de Rubin
        std::array<double, M_RUBIN_FINO> s_m_p2a{}, ee_m_p2a{};
        double max_anc_p2a = 0.0;
        std::size_t max_niv_b4 = 0, max_tmp_b4 = 0, ke_b4_sum = 0, rec_b4 = 0, chq_b4 = 0, pau_b4 = 0;
        for (std::size_t m_idx = 0; m_idx < M_RUBIN_FINO; ++m_idx) {
            std::vector<double> xd, xl;
            std::vector<std::uint8_t> okd, okl, fl;
            std::size_t nrec = 0, nchq = 0, npau = 0;
            double ma = 0.0;
            preparar_flujos_competidor(CompetidorF4::B4_PROPUESTA_UNIFICADA, Integracion::NIVEL, mod3, m3.x_op_limpia,
                                       corr_p2a, 4, ALFA_HUECO_F4, SEMILLA_CIEGA_BASE, m_idx, M_RUBIN_FINO,
                                       xd, okd, xl, okl, fl, nrec, nchq, npau, ma);
            max_anc_p2a = std::max(max_anc_p2a, ma);
            if (m_idx == 0) { rec_b4 = nrec; chq_b4 = nchq; pau_b4 = npau; }
            ConfiguracionLilyth cl_c = cl3;
            cl_c.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
            cl_c.politica_z3_piloto = ConfiguracionLilyth::PoliticaZ3Piloto::RELOJ_FISICO_TOPE_R;
            cl_c.r_max_z3_piloto = 4;
            const auto out = correr_operador_sobre_flujo(od3, cd3, cl_c, m3, xd, okd, xl, okl, fl, en_choque3);
            s_m_p2a[m_idx] = out.s_onsager;
            ee_m_p2a[m_idx] = out.ee_s_onsager;
            max_niv_b4 = std::max(max_niv_b4, out.alarmas_nivel_espurias);
            max_tmp_b4 = std::max(max_tmp_b4, out.alarmas_temp_espurias);
            ke_b4_sum += out.k_evaluados;
        }
        const auto rub_p2a = combinar_rubin(s_m_p2a, ee_m_p2a);
        const double ds_rub_p2a = rub_p2a.valido ? std::abs(rub_p2a.s_barra - ref3.s_onsager) / std::max(1e-12, rub_p2a.ee_rubin) : 99.0;
        std::printf(" %-23s | %11zu | %11zu | %6.1f | %7.4f | %11.2e | rec=%zu chq=%zu pau=%zu gamma=%.4f\n",
                    "B4_RUBIN16", max_niv_b4, max_tmp_b4, static_cast<double>(ke_b4_sum) / M_RUBIN_FINO,
                    ds_rub_p2a, max_anc_p2a, rec_b4, chq_b4, pau_b4, rub_p2a.gamma_info_perdida);
        reportar("CF4.P2a", "NIVEL M3 (semilla ciega) micro-cortes r=4 bajo Rubin M=16 (|dS|/ee_Rubin<=0.25, AlmTmp<=B2, anc<=1e-10)",
                 rub_p2a.valido && (ds_rub_p2a <= 0.25) && (max_tmp_b4 <= esp_tmp_b2) && (max_anc_p2a <= 1e-10));

        // --- CF4.P2b: Matriz 4x5 choque oculto en NIVEL M3 con historia causal real (K=20 posiciones tranquilas) ---
        std::vector<std::size_t> pos_calmas;
        for (std::size_t t = 200; t + 20 < nop && pos_calmas.size() < 20; t += 240) {
            bool libre = true;
            for (std::size_t k = 0; k < 10; ++k) if (en_choque3[t + k]) libre = false;
            if (libre) pos_calmas.push_back(t);
        }
        const std::size_t K_p2b = pos_calmas.size();
        bool p2b_ok = (K_p2b == 20);
        std::printf("\n--- [P2b] Matriz choque oculto NIVEL M3 (alfa_hueco=%.2f, z_max=%.4f, hist=64, K=%zu posiciones) ---\n",
                    ALFA_HUECO_F4, z_max_sel, K_p2b);
        std::printf(" r | Delta_min(80%%) |    +4s (bloq)   |    +6s (bloq)   |    +8s (bloq)   |   +10s (bloq)   |   +12s (bloq)\n");
        std::printf("---+----------------+-----------------+-----------------+-----------------+-----------------+-----------------\n");
        for (std::size_t r = 1; r <= 4; ++r) {
            ConfiguracionPuenteHistoricoIto cp_m{};
            cp_m.integracion = Integracion::NIVEL;
            cp_m.r_max = 4;
            cp_m.alfa_hueco = ALFA_HUECO_F4;
            cp_m.n_semillas_rng = 1;
            const double dmin = delta_minima_detectable_ito(cp_m, mod3, r, 0.80);
            std::printf(" %zu | %13.2fs |", r, dmin);
            for (double amp : {4.0, 6.0, 8.0, 10.0, 12.0}) {
                double z_sum = 0.0;
                std::size_t n_bloq = 0;
                for (std::size_t pi = 0; pi < K_p2b; ++pi) {
                    const std::size_t th = pos_calmas[pi];
                    std::vector<double> hr(64, 0.0);
                    for (std::size_t j = 0; j < 64; ++j) hr[j] = m3.x_op_limpia[th - 64 + j] - m3.x_op_limpia[th - 65 + j];
                    const auto res = evaluar_y_reconstruir_hueco_ito(
                        cp_m, mod3, hr, m3.x_op_limpia[th - 1], m3.x_op_limpia[th + r] + amp * sd_e3, r, 700 + pi);
                    z_sum += res.z_hueco;
                    if (!res.reconstruido && res.choque_en_hueco && ((res.flags_emitidos & FLAG_CHOQUE_EN_HUECO) != 0u)) {
                        ++n_bloq;
                    }
                }
                const double z_med = (K_p2b > 0) ? (z_sum / static_cast<double>(K_p2b)) : 0.0;
                const double pct = (K_p2b > 0) ? (100.0 * static_cast<double>(n_bloq) / static_cast<double>(K_p2b)) : 0.0;
                std::printf(" %5.2f (%5.1f%%) |", z_med, pct);
                if (amp >= dmin && (n_bloq * 5 < K_p2b * 4)) p2b_ok = false;
                if ((r == 1 || r == 2) && amp >= 12.0 && (n_bloq != K_p2b)) p2b_ok = false;
            }
            std::printf("\n");
        }
        reportar("CF4.P2b", "Matriz 4x5 choque oculto (K=20, hist=64): bloqueo >=80% en Delta>=Delta_min y 100% en +12s r=1..2", p2b_ok);

        // --- CF4.P2c: Varianza intra-hueco real (B2 vs B4) + telemetría continua s^2 temp + choque oculto +10s (B3 vs B4) ---
        std::vector<std::size_t> pos_p2c;
        for (std::size_t t = 220; t + 20 < nop && pos_p2c.size() < 24; t += 200) {
            bool libre = true;
            for (std::size_t k = 0; k < 8; ++k) if (en_choque3[t + k]) libre = false;
            if (libre) pos_p2c.push_back(t);
        }
        constexpr std::size_t r_p2c = 4;
        double var_phys_sum = 0.0, var_b2_sum = 0.0, var_b4_sum = 0.0;
        std::size_t n_rec_p2c = 0;
        std::vector<std::uint8_t> corr_p2c_limp(nop, 0u);
        for (std::size_t idx_g = 0; idx_g < pos_p2c.size(); ++idx_g) {
            const std::size_t th = pos_p2c[idx_g];
            for (std::size_t k = 0; k < r_p2c; ++k) corr_p2c_limp[th + k] = 1u;
            std::vector<double> hr(64, 0.0);
            for (std::size_t j = 0; j < 64; ++j) hr[j] = m3.x_op_limpia[th - 64 + j] - m3.x_op_limpia[th - 65 + j];

            ConfiguracionPuenteHistoricoIto cpd{};
            cpd.integracion = Integracion::NIVEL;
            cpd.r_max = 4;
            cpd.alfa_hueco = ALFA_HUECO_F4;
            cpd.n_semillas_rng = M_RUBIN_FINO;
            cpd.semilla_base = SEMILLA_CIEGA_BASE + 500ULL + idx_g;

            const double x_ini = m3.x_op_limpia[th - 1];
            const double x_fin = m3.x_op_limpia[th + r_p2c];
            const auto rb4 = evaluar_y_reconstruir_hueco_ito(cpd, mod3, hr, x_ini, x_fin, r_p2c, idx_g + 1);
            if (!rb4.reconstruido) continue;

            std::array<double, r_p2c + 1> d_phys{}, d_b2{};
            const double paso_b2 = (x_fin - x_ini) / static_cast<double>(r_p2c + 1);
            for (std::size_t s = 0; s <= r_p2c; ++s) {
                d_phys[s] = m3.x_op_limpia[th + s] - m3.x_op_limpia[th + s - 1];
                d_b2[s] = paso_b2;
            }
            var_phys_sum += varianza_muestral(d_phys.data(), r_p2c + 1);
            var_b2_sum += varianza_muestral(d_b2.data(), r_p2c + 1);

            double var_b4_g = 0.0;
            for (std::size_t m_idx = 0; m_idx < rb4.n_semillas_efectivas; ++m_idx) {
                std::array<double, r_p2c + 1> d_b4{};
                for (std::size_t s = 0; s < r_p2c; ++s) d_b4[s] = rb4.y_est_rng[m_idx][s];
                d_b4[r_p2c] = x_fin - rb4.x_rng[m_idx][r_p2c - 1];
                var_b4_g += varianza_muestral(d_b4.data(), r_p2c + 1);
            }
            var_b4_sum += var_b4_g / static_cast<double>(rb4.n_semillas_efectivas);
            ++n_rec_p2c;
        }
        const double var_ratio_b2 = (var_phys_sum > 1e-15) ? (var_b2_sum / var_phys_sum) : 99.0;
        const double var_ratio_b4 = (var_phys_sum > 1e-15) ? (var_b4_sum / var_phys_sum) : 0.0;

        // Telemetría continua informativa del canal de temperatura (sd.temp[0].s2) sobre ventanas tocadas por huecos r=4
        auto medir_s2_huecos = [&](CompetidorF4 comp) {
            std::vector<double> xd, xl;
            std::vector<std::uint8_t> okd, okl, fl;
            std::size_t nr = 0, nc = 0, np = 0; double ma = 0.0;
            if (comp == CompetidorF4::B0_BASE) {
                xd = m3.x_op_limpia; xl = m3.x_op_limpia; okd = ok1; okl = ok1; fl = fl0;
            } else {
                preparar_flujos_competidor(comp, Integracion::NIVEL, mod3, m3.x_op_limpia, corr_p2c_limp, 4,
                                           ALFA_HUECO_F4, SEMILLA_CIEGA_BASE, 0, M_RUBIN_FINO, xd, okd, xl, okl, fl, nr, nc, np, ma);
            }
            std::vector<std::pair<std::size_t, double>> s2_v;
            (void)correr_operador_sobre_flujo(od3, cd3, cl3, m3, xd, okd, xl, okl, fl, en_choque3, &s2_v);
            double sum_s2 = 0.0; std::size_t cnt_s2 = 0;
            for (auto [t_em, s2_val] : s2_v) {
                bool toca = false;
                for (std::size_t th : pos_p2c) {
                    if (t_em >= th && t_em <= th + 32) { toca = true; break; }
                }
                if (toca && std::isfinite(s2_val)) { sum_s2 += s2_val; ++cnt_s2; }
            }
            return (cnt_s2 > 0) ? (sum_s2 / static_cast<double>(cnt_s2)) : NaN;
        };
        const double s2_ref = medir_s2_huecos(CompetidorF4::B0_BASE);
        const double s2_b1b = medir_s2_huecos(CompetidorF4::B1B_ARRASTRE_NIVEL);
        const double s2_b2  = medir_s2_huecos(CompetidorF4::B2_LINEAL_2BORDES);
        const double s2_b4  = medir_s2_huecos(CompetidorF4::B4_PROPUESTA_UNIFICADA);

        // Parte (ii) de P2c: 10 micro-cortes r=2 con choques ocultos +10 sigma_e en NIVEL M3 (B3 vs B4 en ObservadorLilyth)
        std::vector<double> x_m3_ocultos = m3.x_op_limpia;
        std::vector<std::uint8_t> corr_ocultos(nop, 0u);
        std::size_t n_iny_ocultos = 0;
        for (std::size_t i = 0; i < pos_calmas.size() && n_iny_ocultos < 10; i += 2) {
            const std::size_t th = pos_calmas[i];
            corr_ocultos[th] = 1u;
            corr_ocultos[th + 1] = 1u;
            for (std::size_t t = th; t < nop; ++t) x_m3_ocultos[t] += 10.0 * sd_e3;
            ++n_iny_ocultos;
        }
        std::vector<double> xd_b3, xl_b3, xd_b4, xl_b4;
        std::vector<std::uint8_t> okd_b3, okl_b3, fl_b3, okd_b4, okl_b4, fl_b4;
        std::size_t rec_o3 = 0, chq_o3 = 0, pau_o3 = 0, rec_o4 = 0, chq_o4 = 0, pau_o4 = 0;
        double ma3 = 0.0, ma4 = 0.0;
        preparar_flujos_competidor(CompetidorF4::B3_ITO_SIN_FILTRO, Integracion::NIVEL, mod3, x_m3_ocultos, corr_ocultos, 2,
                                   ALFA_HUECO_F4, SEMILLA_CIEGA_BASE, 0, M_RUBIN_FINO, xd_b3, okd_b3, xl_b3, okl_b3, fl_b3, rec_o3, chq_o3, pau_o3, ma3);
        preparar_flujos_competidor(CompetidorF4::B4_PROPUESTA_UNIFICADA, Integracion::NIVEL, mod3, x_m3_ocultos, corr_ocultos, 2,
                                   ALFA_HUECO_F4, SEMILLA_CIEGA_BASE, 0, M_RUBIN_FINO, xd_b4, okd_b4, xl_b4, okl_b4, fl_b4, rec_o4, chq_o4, pau_o4, ma4);
        ConfiguracionLilyth cl_o4 = cl3;
        cl_o4.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cl_o4.politica_z3_piloto = ConfiguracionLilyth::PoliticaZ3Piloto::RELOJ_FISICO_TOPE_R;
        cl_o4.r_max_z3_piloto = 2;
        const auto out_o3 = correr_operador_sobre_flujo(od3, cd3, cl3, m3, xd_b3, okd_b3, xl_b3, okl_b3, fl_b3, en_choque3);
        const auto out_o4 = correr_operador_sobre_flujo(od3, cd3, cl_o4, m3, xd_b4, okd_b4, xl_b4, okl_b4, fl_b4, en_choque3);

        std::printf("\n--- [P2c] Discriminante Calmas Inventadas (K=%zu r=4) y Choques Ocultos +10s (K=%zu r=2) ---\n", n_rec_p2c, n_iny_ocultos);
        std::printf("  Varianza intra-hueco / fisica : B2_LINEAL = %.4f | B4_ITO_RNG = %.4f (rango [0.60, 1.40])\n", var_ratio_b2, var_ratio_b4);
        std::printf("  Telemetria continua s^2(temp) : LIMPIO=%.4f | B1b=%.4f (%.2fx) | B2=%.4f (%.2fx) | B4=%.4f (%.2fx)\n",
                    s2_ref, s2_b1b, s2_b1b / std::max(1e-12, s2_ref), s2_b2, s2_b2 / std::max(1e-12, s2_ref),
                    s2_b4, s2_b4 / std::max(1e-12, s2_ref));
        std::printf("  Choques ocultos +10s (r=2)    : B3 bloquea %zu/%zu (K_eval=%zu) | B4 bloquea %zu/%zu (K_eval=%zu)\n",
                    chq_o3, n_iny_ocultos, out_o3.k_evaluados, chq_o4, n_iny_ocultos, out_o4.k_evaluados);

        const bool p2c_ok = (n_rec_p2c >= 16) && (var_ratio_b2 < 0.10) && (var_ratio_b4 >= 0.60 && var_ratio_b4 <= 1.40) &&
                            (chq_o3 == 0) && (chq_o4 >= 8);
        reportar("CF4.P2c", "Discriminacion efectiva: B2 colapsa varianza intra-hueco (<0.10), B4 la preserva [0.60,1.40] y bloquea >=8/10 choques +10s vs 0/10 de B3", p2c_ok);
    }

    // ------------------------------------------------------------------------
    // CF4.P3: M2 ARFIMA (16 semillas ciegas, forzar_via_b=true, Z1b dt=2 r=3, Rubin M=16 POR SERIE)
    // ------------------------------------------------------------------------
    {
        constexpr std::size_t S_P3 = 16;
        struct ResP3Seed {
            bool ok = false;
            double ee0 = 0.0, s0 = 0.0, d_hat = 0.25, max_anc = 0.0;
            double rm_b0 = 0.0, ds_b0 = 0.0;
            double rm_cs = 0.0, ds_cs = 0.0;
            double rm_b1 = 0.0, ds_b1 = 0.0;
            double rm_lin = 0.0, ds_lin = 0.0;
            double rm_ito = 0.0, ds_ito = 0.0;
            double rm_rub = 0.0, ds_rub = 0.0, ee_rub = 0.0, gamma_rub = 0.0;
        };
        std::vector<ResP3Seed> res_s(S_P3);

        en_paralelo(S_P3, n_hilos, [&](std::size_t s) {
            const std::uint64_t sem_s = SEMILLA_CIEGA_BASE + 31337ULL * s;
            ParametrosBancoComun params{};
            params.n_op = 6000;
            const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, params, sem_s);

            ConfiguracionObservador cd2{};
            cd2.integracion = Integracion::FLUJO;
            cd2.agregacion = AgregacionEscala::SUMA;
            cd2.arl0_pasos_nativos = 450.0;
            cd2.alfa = ALFA_PROYECTO;
            cd2.medir_tau50 = false;
            cd2.semilla = sem_s;
            cd2.escalas_max = 3;
            ObsDelta od2{};
            od2.calibrar(cd2, m2.x_cal);
            if (od2.escalas() == 0) return;

            const std::size_t m_freq = std::max<std::size_t>(32, static_cast<std::size_t>(std::pow(m2.x_cal.size(), 0.65)));
            const auto est_d = elw::estimar_calibracion(m2.x_cal.data(), m2.x_cal.size(), m_freq, -0.45, 0.45);

            ConfiguracionLilyth cl2{};
            cl2.alfa = ALFA_PROYECTO;
            cl2.representacion = ModoRepresentacionLilyth::CRUDA;
            cl2.b_max_nulo = 16;
            cl2.hilos_nulo = 1;
            cl2.semilla_nulo = sem_s ^ 0x51A7E123ULL;
            cl2.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
            cl2.d_whittle_cruda = est_d.valido ? est_d.d : 0.25;
            cl2.forzar_via_b = true;
            cl2.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
            cl2.b_mc_resumen = 1999;

            const std::size_t nop = m2.x_op_limpia.size();
            const std::vector<std::uint8_t> ok1(nop, 1u), fl0(nop, 0u);
            const auto ref2 = correr_operador_sobre_flujo(od2, cd2, cl2, m2, m2.x_op_limpia, ok1, m2.x_op_limpia, ok1, fl0);
            if (!ref2.valida) return;

            ObsDelta od2_cal = od2;
            od2_cal.reiniciar_flujo();
            ObsLil lil_mod2{};
            (void)lil_mod2.calibrar(od2_cal, cl2, cd2, m2.x_cal);
            const auto& mod2 = lil_mod2.ficha_nulo_cruda().modelo_lineal;
            if (!mod2.valido) return;

            std::vector<std::uint8_t> ok_z1b(nop, 1u);
            std::vector<std::pair<std::int64_t, std::int64_t>> huecos_z1b;
            for (const auto& [clv, tc] : ref2.t_cambio_por_ev) {
                (void)clv;
                const std::int64_t idx_fin = tc - 2;
                const std::int64_t idx_ini = idx_fin - 2;
                if (idx_ini >= 10 && idx_fin + 1 < static_cast<std::int64_t>(nop)) {
                    for (std::int64_t t = idx_ini; t <= idx_fin; ++t) ok_z1b[static_cast<std::size_t>(t)] = 0u;
                    huecos_z1b.emplace_back(idx_ini, idx_fin);
                }
            }

            auto eval_pol = [&](ConfiguracionLilyth::PoliticaZ1Piloto pol) {
                ConfiguracionLilyth cl_p = cl2;
                cl_p.politica_z1_piloto = pol;
                cl_p.r_max_z1_piloto = 5;
                cl_p.alfa_hueco_piloto = ALFA_HUECO_F4;
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

            // B4_RUBIN16: M=16 imputaciones reales por hueco, combinadas con Rubin DENTRO de esta serie s
            std::array<double, M_RUBIN_FINO> s_m{}, ee_m{};
            double rm_rng_sum = 0.0;
            for (std::size_t m_idx = 0; m_idx < M_RUBIN_FINO; ++m_idx) {
                std::vector<double> x_r = m2.x_op_limpia;
                std::vector<std::uint8_t> fl_r(nop, 0u);
                for (std::size_t hi = 0; hi < huecos_z1b.size(); ++hi) {
                    const auto [ini, fin] = huecos_z1b[hi];
                    const std::size_t u_ini = static_cast<std::size_t>(ini);
                    const std::size_t h_len = std::min<std::size_t>(256, u_ini);
                    const std::span<const double> hp(&x_r[u_ini - h_len], h_len);
                    ConfiguracionPuenteHistoricoIto cp{};
                    cp.integracion = Integracion::FLUJO;
                    cp.r_max = 3;
                    cp.alfa_hueco = ALFA_HUECO_F4;
                    cp.n_semillas_rng = M_RUBIN_FINO;
                    cp.semilla_base = sem_s ^ 0x9E3779B9ULL;
                    const auto res = evaluar_y_reconstruir_hueco_ito(cp, mod2, hp, x_r[u_ini - 1], m2.x_op_limpia[static_cast<std::size_t>(fin + 1)], 3, hi + 1);
                    if (res.max_error_anclaje_rng > rs.max_anc) rs.max_anc = res.max_error_anclaje_rng;
                    if (res.reconstruido) {
                        for (std::size_t k = 0; k < 3; ++k) {
                            x_r[u_ini + k] = res.x_rng[m_idx][k];
                            fl_r[u_ini + k] = FLAG_PUENTE_ITO;
                        }
                    }
                }
                const auto out_r = correr_operador_sobre_flujo(od2, cd2, cl2, m2, m2.x_op_limpia, ok1, x_r, ok1, fl_r);
                s_m[m_idx] = out_r.s_onsager;
                ee_m[m_idx] = out_r.ee_s_onsager;
                rm_rng_sum += calcular_rmse_e(ref2, out_r);
            }
            const auto rub_s = combinar_rubin(s_m, ee_m);
            if (!rub_s.valido) return;
            rs.rm_rub = rm_rng_sum / static_cast<double>(M_RUBIN_FINO);
            rs.ds_rub = rub_s.s_barra - ref2.s_onsager;
            rs.ee_rub = rub_s.ee_rubin;
            rs.gamma_rub = rub_s.gamma_info_perdida;
            rs.ok = true;
        });

        std::size_t n_ok_p3 = 0;
        double d_med = 0.0, s0_med = 0.0, ee0_med = 0.0, ee_rub_med = 0.0, gamma_med = 0.0, max_anc_p3 = 0.0;
        double rm_b0 = 0.0, ds_b0 = 0.0, rm_cs = 0.0, ds_cs = 0.0, rm_b1 = 0.0, ds_b1 = 0.0;
        double rm_lin = 0.0, ds_lin = 0.0, rm_ito = 0.0, ds_ito = 0.0, rm_rub = 0.0, ds_rub = 0.0;
        for (const auto& rs : res_s) {
            if (!rs.ok) continue;
            ++n_ok_p3;
            d_med += rs.d_hat; s0_med += rs.s0; ee0_med += rs.ee0; ee_rub_med += rs.ee_rub; gamma_med += rs.gamma_rub;
            max_anc_p3 = std::max(max_anc_p3, rs.max_anc);
            rm_b0 += rs.rm_b0; ds_b0 += rs.ds_b0;
            rm_cs += rs.rm_cs; ds_cs += rs.ds_cs;
            rm_b1 += rs.rm_b1; ds_b1 += rs.ds_b1;
            rm_lin += rs.rm_lin; ds_lin += rs.ds_lin;
            rm_ito += rs.rm_ito; ds_ito += rs.ds_ito;
            rm_rub += rs.rm_rub; ds_rub += rs.ds_rub;
        }
        const double inv_n = (n_ok_p3 > 0) ? (1.0 / static_cast<double>(n_ok_p3)) : 0.0;
        d_med *= inv_n; s0_med *= inv_n; ee0_med = std::max(1e-12, ee0_med * inv_n);
        ee_rub_med = std::max(1e-12, ee_rub_med * inv_n); gamma_med *= inv_n;
        rm_b0 *= inv_n; ds_b0 *= inv_n; rm_cs *= inv_n; ds_cs *= inv_n; rm_b1 *= inv_n; ds_b1 *= inv_n;
        rm_lin *= inv_n; ds_lin *= inv_n; rm_ito *= inv_n; ds_ito *= inv_n; rm_rub *= inv_n; ds_rub *= inv_n;

        const double b_b0  = std::abs(ds_b0)  / ee0_med;
        const double b_cs  = std::abs(ds_cs)  / ee0_med;
        const double b_b1  = std::abs(ds_b1)  / ee0_med;
        const double b_lin = std::abs(ds_lin) / ee0_med;
        const double b_ito = std::abs(ds_ito) / ee0_med;
        const double b_rub = std::abs(ds_rub) / ee_rub_med;

        std::printf("\n--- [P3] M2 ARFIMA (%zu/%zu semillas ciegas, d_med=%.4f, S0_med=%.4f, ee0=%.4f, ee_Rubin=%.4f) ---\n",
                    n_ok_p3, S_P3, d_med, s0_med, ee0_med, ee_rub_med);
        std::printf(" Politica Z1b (dt=2, r=3) | RMSE(e)/sd | |E[dS]|/ee   | ESTRICTO | DULCE\n");
        std::printf("--------------------------+------------+--------------+----------+-------\n");
        auto imp_p3 = [&](const char* nom, double rm, double bs) {
            std::printf(" %-24s | %10.4f | %12.4f | %-8s | %-5s\n", nom, rm, bs,
                        (rm <= 0.15 && bs <= 0.25) ? "PASA" : "FALLA",
                        (rm < 0.50 && bs <= 1.00) ? "PASA" : "FALLA");
        };
        imp_p3("B0_COMPRESION_BASE",      rm_b0,  b_b0);
        imp_p3("B0b_CORTE_SUFIJO",        rm_cs,  b_cs);
        imp_p3("B1_PUENTE_CAUSAL_1B",     rm_b1,  b_b1);
        imp_p3("B2_LINEAL_2BORDES",       rm_lin, b_lin);
        imp_p3("B4_PUENTE_ITO_MEDIA_2B",  rm_ito, b_ito);
        std::printf(" %-24s | %10.4f | %12.4f | %-8s | %-5s  (gamma_med=%.4f, err_anc=%.2e)\n",
                    "B4_PUENTE_ITO_RUBIN16", rm_rub, b_rub,
                    (rm_rub <= 0.15 && b_rub <= 0.25 && gamma_med > 1e-6 && gamma_med <= 0.25) ? "PASA" : "FALLA",
                    (rm_rub < 0.50 && b_rub <= 1.00) ? "PASA" : "FALLA",
                    gamma_med, max_anc_p3);

        const bool p3_ok = (n_ok_p3 >= 12) && (rm_rub <= 0.15) && (rm_rub <= rm_cs + 1e-12) &&
                           (b_rub <= 0.25) && (gamma_med > 1e-6 && gamma_med <= 0.25) && (max_anc_p3 <= 1e-10);
        reportar("CF4.P3", "M2 ARFIMA (16 semillas ciegas, Rubin M=16 por serie): B4_RUBIN16 pasa ESTRICTO y gamma en (1e-6, 0.25]", p3_ok);
    }

    // ------------------------------------------------------------------------
    // CF4.P4: Pausa dura exacta en corte largo r = 12 > r_max = 3
    // ------------------------------------------------------------------------
    {
        AdaptadorInputCausalLilyth ad{};
        ad.integracion = Integracion::NIVEL;
        ad.racha_max_corrupcion = 3;
        ad.z_max_hueco = z_max_sel;
        ad.evitar_teletransporte_nivel = true;

        (void)ad.preparar(100.0, ClasePasoLilyth::VALIDO);
        std::size_t pasos_en_pausa = 0;
        for (int k = 0; k < 12; ++k) {
            const auto st = ad.preparar(NaN, ClasePasoLilyth::CORRUPCION_O_CAIDA);
            if (st.estado == EstadoContinuidadLilyth::EN_PAUSA_DURA &&
                !st.valido_para_operador &&
                ((st.flags_emitidos & FLAG_PAUSA_UMBRAL) != 0u)) {
                ++pasos_en_pausa;
            }
        }
        const auto rec = ad.preparar(105.0, ClasePasoLilyth::VALIDO);
        const bool p4_ok = (pasos_en_pausa == 9) && rec.reconexion_tras_hueco && rec.pausa_umbral_activada &&
                           rec.cortar_segmento_al_reenganchar && (rec.delta_t_hueco == 13);
        reportar("CF4.P4", "Corte largo r=12 > r_max=3 mantiene EN_PAUSA_DURA exacta (9 pasos) y corta segmento al reconectar", p4_ok);
    }

    // ------------------------------------------------------------------------
    // CF4.P6: Blindaje obligatorio de FLAG_PUENTE_ITO en Z2 && Z3 dentro de ObservadorLilyth (caza MN1.7)
    // ------------------------------------------------------------------------
    {
        ParametrosBancoComun pb{};
        pb.n_op = 6000;
        const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, pb, SEMILLA_CIEGA_BASE + 7777ULL);
        const std::size_t m_freq = std::max<std::size_t>(32, static_cast<std::size_t>(std::pow(m2.x_cal.size(), 0.65)));
        const auto est_d = elw::estimar_calibracion(m2.x_cal.data(), m2.x_cal.size(), m_freq, -0.45, 0.45);

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.agregacion = AgregacionEscala::SUMA;
        cd.arl0_pasos_nativos = 450.0;
        cd.alfa = ALFA_PROYECTO;
        cd.medir_tau50 = false;
        cd.semilla = SEMILLA_CIEGA_BASE + 7777ULL;
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
        cl.semilla_nulo = SEMILLA_CIEGA_BASE + 7777ULL;
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
        reportar("CF4.P6", "Blindaje obligatorio de FLAG_PUENTE_ITO en Z2 && Z3 dentro de ObservadorLilyth", veto_z2 && veto_z3);
    }

    // ------------------------------------------------------------------------
    // CF4.P7: Deriva no nula en FLUJO (mu = +5.0*sigma_e) vía valor_z0_imputado (caza MN1.10)
    // ------------------------------------------------------------------------
    {
        const double phi_d = 0.65;
        const double sd_e_d = std::sqrt(1.0 - phi_d * phi_d);
        const double mu_d = 5.0 * sd_e_d;
        const std::size_t n_cal_d = 3200, n_op_d = 3200;
        std::mt19937_64 rng_d(detalle_ito::splitmix64_ito(SEMILLA_CIEGA_BASE ^ 0xF00DCAFEULL));
        std::normal_distribution<double> nd_d(0.0, sd_e_d);

        std::vector<double> x_cal_d(n_cal_d), x_op_d(n_op_d);
        double u = 0.0;
        for (std::size_t s = 0; s < 200; ++s) u = phi_d * u + nd_d(rng_d);
        for (std::size_t t = 0; t < n_cal_d; ++t) { u = phi_d * u + nd_d(rng_d); x_cal_d[t] = mu_d + u; }
        for (std::size_t t = 0; t < n_op_d; ++t)  { u = phi_d * u + nd_d(rng_d); x_op_d[t] = mu_d + u; }

        ConfiguracionObservador cd_d{};
        cd_d.integracion = Integracion::FLUJO;
        cd_d.agregacion = AgregacionEscala::SUMA;
        cd_d.arl0_pasos_nativos = 450.0;
        cd_d.alfa = ALFA_PROYECTO;
        cd_d.medir_tau50 = false;
        cd_d.semilla = SEMILLA_CIEGA_BASE ^ 0xF00DCAFEULL;
        cd_d.escalas_max = 2;
        ObsDelta od_d{};
        od_d.calibrar(cd_d, x_cal_d);

        ConfiguracionLilyth cl_d{};
        cl_d.alfa = ALFA_PROYECTO;
        cl_d.representacion = ModoRepresentacionLilyth::CRUDA;
        cl_d.b_max_nulo = 16;
        cl_d.hilos_nulo = 1;
        ObsLil lil_d{};
        (void)lil_d.calibrar(od_d, cl_d, cd_d, x_cal_d);
        const ModeloLinealNulo mod_d = lil_d.ficha_nulo_cruda().modelo_lineal;

        std::vector<std::uint8_t> es_c_d(n_op_d, 0u), cerca_hueco(n_op_d, 0u);
        for (std::size_t k = 0; k < 12; ++k) {
            const std::size_t th = 400 + k * 200;
            es_c_d[th] = 1u;
            es_c_d[th + 1] = 1u;
            for (std::size_t w = th; w < std::min(n_op_d, th + 10); ++w) cerca_hueco[w] = 1u;
        }
        std::vector<double> xd_d, xl_d;
        std::vector<std::uint8_t> okd_d, okl_d, fl_d;
        std::size_t nrec_d = 0, nchq_d = 0, npau_d = 0; double ma_d = 0.0;
        preparar_flujos_competidor(CompetidorF4::B4_PROPUESTA_UNIFICADA, Integracion::FLUJO, mod_d, x_op_d, es_c_d, 2,
                                   ALFA_HUECO_F4, SEMILLA_CIEGA_BASE ^ 0xF00DCAFEULL, 0, 4, xd_d, okd_d, xl_d, okl_d, fl_d,
                                   nrec_d, nchq_d, npau_d, ma_d);

        od_d.reiniciar_flujo();
        std::size_t alarmas_neg_hueco = 0;
        for (std::size_t t = 0; t < n_op_d; ++t) {
            const bool ok = (okd_d[t] != 0u) && std::isfinite(xd_d[t]);
            const auto sd = od_d.paso(ok ? xd_d[t] : NaN, NaN, ok, 0, 0);
            if (cerca_hueco[t]) {
                for (std::size_t j = 0; j < sd.niveles; ++j) {
                    if (sd.nivel[j].alarma.hay && sd.nivel[j].alarma.signo < 0) ++alarmas_neg_hueco;
                }
            }
        }
        std::printf("  [Detalle P7 deriva FLUJO mu=%.3f (cal=%.3f)] huecos_rec=%zu | alarmas_neg_en_hueco=%zu\n",
                    mu_d, mod_d.mu, nrec_d, alarmas_neg_hueco);
        reportar("CF4.P7", "Deriva no nula en FLUJO (mu=+5s) via valor_z0_imputado: 0 alarmas negativas espurias en 12 huecos r=2",
                 (nrec_d >= 10) && (alarmas_neg_hueco == 0));
    }

    // ------------------------------------------------------------------------
    // CF4.P5: Bloque ciego Monte Carlo R=R_REP (M1_H0, M2_H0, M3_H0, M1_H1) con Rubin M=5 + ark::holm()
    // ------------------------------------------------------------------------
    {
        ParametrosReplicaH0H1 ph0{};
        ph0.lenta_h1 = false;
        ph0.semilla_base = SEMILLA_CIEGA_BASE;
        ParametrosReplicaH0H1 ph1{};
        ph1.lenta_h1 = true;
        ph1.semilla_base = SEMILLA_CIEGA_BASE;

        std::vector<SerieSinteticaCoordinada> m1_h0(R_REP), m2_h0(R_REP), m3_h0(R_REP), m1_h1(R_REP);
        en_paralelo(R_REP, n_hilos, [&](std::size_t i) {
            m1_h0[i] = generar_mundo_m1_h0_h1(0.80, ph0, i);
            m2_h0[i] = generar_mundo_m2_h0_h1(0.25, 0.40, ph0, i);
            m3_h0[i] = generar_mundo_m3_h0_h1(0.65, ph0, i);
            m1_h1[i] = generar_mundo_m1_h0_h1(0.80, ph1, i);
        });

        auto correr_lote_mc = [&](const std::vector<SerieSinteticaCoordinada>& mundos,
                                  Integracion integ,
                                  bool purga,
                                  ModoNuloLineal mlin,
                                  bool aplicar_corr,
                                  CompetidorF4 comp) -> ResumenMCF4 {
            struct RepRes {
                bool val = false, r = false, rl = false;
                double ke = 0.0, lz2_sum = 0.0;
                std::size_t lz2_cnt = 0;
                double s = 0.0, ee = 0.0, gam = 0.0;
                std::size_t tmp = 0, vz = 0;
            };
            std::vector<RepRes> rr(R_REP);
            en_paralelo(R_REP, n_hilos, [&](std::size_t i) {
                const auto& m = mundos[i];
                ConfiguracionObservador cd{};
                cd.integracion = integ;
                cd.agregacion = AgregacionEscala::SUMA;
                cd.arl0_pasos_nativos = 450.0;
                cd.alfa = ALFA_PROYECTO;
                cd.purga = purga;
                cd.bins_fase = purga ? 24 : 1;
                cd.dias_minimos_purga = 20;
                cd.medir_tau50 = false;
                cd.semilla = SEMILLA_CIEGA_BASE ^ (i * 100003ULL);
                cd.escalas_max = 2;

                ObsDelta od{};
                od.calibrar(cd, m.x_cal, {}, {},
                            purga ? std::span<const std::size_t>(m.fase_cal) : std::span<const std::size_t>{},
                            purga ? std::span<const std::int64_t>(m.dia_cal) : std::span<const std::int64_t>{});
                if (od.escalas() == 0) return;

                ConfiguracionLilyth cl{};
                cl.alfa = ALFA_PROYECTO;
                cl.representacion = purga ? ModoRepresentacionLilyth::PURGADA : ModoRepresentacionLilyth::CRUDA;
                cl.b_max_nulo = 32;
                cl.hilos_nulo = 1;
                cl.semilla_nulo = SEMILLA_CIEGA_BASE ^ (i * 31337ULL);
                cl.modo_nulo_lineal = mlin;
                if (mlin == ModoNuloLineal::ARFIMA_WHITTLE) cl.d_whittle_cruda = 0.25;
                cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
                cl.pool_nulo = PoolNuloLilyth::DERIVADO;
                cl.horizonte = HorizonteLilyth::CUANTIL_RESOLUCION;
                cl.b_mc_resumen = 1999;

                if (aplicar_corr && (comp == CompetidorF4::B1_ARGON_VIVO || comp == CompetidorF4::B4_PROPUESTA_UNIFICADA)) {
                    cl.politica_z1_piloto = (comp == CompetidorF4::B4_PROPUESTA_UNIFICADA)
                        ? ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO
                        : ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_CAUSAL;
                    cl.r_max_z1_piloto = 2;
                    cl.alfa_hueco_piloto = ALFA_HUECO_F4;
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
                    std::mt19937_64 rg(detalle_ito::splitmix64_ito(SEMILLA_CIEGA_BASE ^ (i * 99991ULL)));
                    std::uniform_real_distribution<double> u01(0.0, 1.0);
                    for (std::size_t t = 10; t + 2 < nop; ++t) {
                        if (!chq[t] && !chq[t + 1] && u01(rg) < 0.01) {
                            es_c[t] = 1u; es_c[t + 1] = 1u;
                            t += 2;
                        }
                    }
                }

                ModeloLinealNulo m_calma{};
                if (aplicar_corr && (comp == CompetidorF4::B1_ARGON_VIVO ||
                                     comp == CompetidorF4::B3_ITO_SIN_FILTRO ||
                                     comp == CompetidorF4::B4_PROPUESTA_UNIFICADA)) {
                    ObsLil lil_cal{};
                    (void)lil_cal.calibrar(od, cl, cd, m.x_cal, {},
                                           purga ? std::span<const std::size_t>(m.fase_cal) : std::span<const std::size_t>{},
                                           purga ? std::span<const std::int64_t>(m.dia_cal) : std::span<const std::int64_t>{});
                    m_calma = purga ? lil_cal.ficha_nulo_purgada().modelo_lineal : lil_cal.ficha_nulo_cruda().modelo_lineal;
                }

                const bool usar_rubin = aplicar_corr && (comp == CompetidorF4::B4_PROPUESTA_UNIFICADA);
                if (!usar_rubin) {
                    std::vector<double> xd, xl;
                    std::vector<std::uint8_t> okd, okl, fl;
                    std::size_t nrec = 0, nchq = 0, npau = 0; double ma = 0.0;
                    preparar_flujos_competidor(comp, integ, m_calma, m.x_op_limpia, es_c, 2, ALFA_HUECO_F4,
                                               SEMILLA_CIEGA_BASE + i, 0, M_RUBIN_MC, xd, okd, xl, okl, fl, nrec, nchq, npau, ma);
                    const auto out = correr_operador_sobre_flujo(od, cd, cl, m, xd, okd, xl, okl, fl);
                    rr[i] = {out.valida, out.rechaza_inter, out.rechaza_lenta, static_cast<double>(out.k_evaluados),
                             out.l_z2_suma, out.l_z2_cuenta, out.s_onsager, out.ee_s_onsager, 0.0,
                             out.alarmas_temp_tot, out.n_veto_z2};
                    return;
                }

                // B4_RUBIN (M = M_RUBIN_MC imputaciones combinadas por Rubin + ark::holm())
                std::array<double, M_RUBIN_MC> sm_v{}, eem_v{};
                std::vector<double> pc0_v, pc2_v, pc3_v;
                double ke_acc = 0.0;
                std::size_t tmp_acc = 0, vz_acc = 0, n_ok_m = 0;
                for (std::size_t m_idx = 0; m_idx < M_RUBIN_MC; ++m_idx) {
                    std::vector<double> xd, xl;
                    std::vector<std::uint8_t> okd, okl, fl;
                    std::size_t nrec = 0, nchq = 0, npau = 0; double ma = 0.0;
                    preparar_flujos_competidor(CompetidorF4::B4_PROPUESTA_UNIFICADA, integ, m_calma, m.x_op_limpia, es_c, 2,
                                               ALFA_HUECO_F4, SEMILLA_CIEGA_BASE + i, m_idx, M_RUBIN_MC,
                                               xd, okd, xl, okl, fl, nrec, nchq, npau, ma);
                    const auto out = correr_operador_sobre_flujo(od, cd, cl, m, xd, okd, xl, okl, fl);
                    if (out.valida) {
                        sm_v[m_idx] = out.s_onsager;
                        eem_v[m_idx] = out.ee_s_onsager;
                        pc0_v.push_back(out.p_crudos[0]);
                        pc2_v.push_back(out.p_crudos[2]);
                        pc3_v.push_back(out.p_crudos[3]);
                        ke_acc += static_cast<double>(out.k_evaluados);
                        tmp_acc += out.alarmas_temp_tot;
                        vz_acc += out.n_veto_z2;
                        ++n_ok_m;
                    } else {
                        sm_v[m_idx] = NaN;
                        eem_v[m_idx] = NaN;
                    }
                }
                if (n_ok_m == 0) return;
                const auto rub = combinar_rubin(sm_v, eem_v);
                if (!rub.valido) return;
                const std::vector<double> pc_comb{mediana_vector(pc0_v), rub.p_bilateral, mediana_vector(pc2_v), mediana_vector(pc3_v)};
                const auto ph_comb = ark::holm(pc_comb);
                const bool holm_ok = (ph_comb.size() == 4) && std::isfinite(ph_comb[1]);
                const bool rech = holm_ok && (ph_comb[1] < ALFA_PROYECTO);
                const bool rech_l = rech && (rub.s_barra > 0.0);
                rr[i] = {holm_ok, rech, rech_l, ke_acc / static_cast<double>(n_ok_m), 0.0, 0,
                         rub.s_barra, rub.ee_rubin, rub.gamma_info_perdida,
                         tmp_acc / n_ok_m, vz_acc / n_ok_m};
            });

            ResumenMCF4 res{};
            double sk = 0.0, ss = 0.0, see = 0.0, st = 0.0, svz = 0.0, sgam = 0.0, slz = 0.0;
            std::size_t clz = 0;
            for (const auto& r : rr) {
                sk += r.ke; st += static_cast<double>(r.tmp); svz += static_cast<double>(r.vz);
                slz += r.lz2_sum; clz += r.lz2_cnt;
                if (!r.val) continue;
                ++res.n_val;
                if (r.r) ++res.n_rech;
                if (r.rl) ++res.n_rech_lenta;
                ss += r.s; see += r.ee; sgam += r.gam;
            }
            res.k_eval_med = sk / static_cast<double>(R_REP);
            res.temp_alarm_med = st / static_cast<double>(R_REP);
            res.veto_z2_med = svz / static_cast<double>(R_REP);
            res.l_z2_med = (clz > 0) ? (slz / static_cast<double>(clz)) : 1.0;
            if (res.n_val > 0) {
                res.tasa_rech = static_cast<double>(res.n_rech) / static_cast<double>(res.n_val);
                res.tasa_lenta = static_cast<double>(res.n_rech_lenta) / static_cast<double>(res.n_val);
                const auto [lo, hi] = wilson_95(res.n_rech, res.n_val);
                res.wlo = lo; res.whi = hi;
                res.media_s = ss / static_cast<double>(res.n_val);
                res.media_ee = std::max(1e-6, see / static_cast<double>(res.n_val));
                res.media_gamma_rubin = sgam / static_cast<double>(res.n_val);
            }
            return res;
        };

        const auto b_m1 = correr_lote_mc(m1_h0, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, false, CompetidorF4::B0_BASE);
        const auto b_m2 = correr_lote_mc(m2_h0, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE,     false, CompetidorF4::B0_BASE);
        const auto b_m3 = correr_lote_mc(m3_h0, Integracion::NIVEL, true,  ModoNuloLineal::AR_LEVINSON_DURBIN, false, CompetidorF4::B0_BASE);
        const auto b_h1 = correr_lote_mc(m1_h1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, false, CompetidorF4::B0_BASE);

        const double p_teor_m1 = 1.0 - std::pow(1.0 - 0.02 / 2.0, b_m1.l_z2_med + 2.0 - 1.0);
        const double cota_cens_m1 = p_teor_m1 + 1.96 * std::sqrt(p_teor_m1 * (1.0 - p_teor_m1) / std::max(1.0, static_cast<double>(R_REP) * b_m1.k_eval_med));

        std::printf("\n--- [P5] BLOQUE CIEGO MONTE CARLO (R=%zu, q=2%%, r=2, M_rubin=%zu + Holm) | l_bar_Z2(M1)=%.2f -> cota_cens=%.3f ---\n",
                    R_REP, M_RUBIN_MC, b_m1.l_z2_med, cota_cens_m1);
        std::printf(" Escenario | Politica               | Validas | K_eval | FalsoRech_H0 [Wilson 95%%]      | RetPot_H1 | |dS|/ee | AlmTmp | ESTRICTO\n");
        std::printf("-----------+------------------------+---------+--------+--------------------------------+-----------+---------+--------+---------\n");
        std::printf(" LIMPIO q0 | M1_H0 BASE             | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] |         - |   0.000 | %6.2f | %-8s\n",
                    b_m1.n_val, b_m1.k_eval_med, b_m1.n_rech, b_m1.n_val, b_m1.tasa_rech, b_m1.wlo, b_m1.whi, b_m1.temp_alarm_med, b_m1.wlo <= 0.05 ? "PASA" : "FALLA");
        std::printf(" LIMPIO q0 | M2_H0 ARFIMA BASE      | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] |         - |   0.000 | %6.2f | %-8s\n",
                    b_m2.n_val, b_m2.k_eval_med, b_m2.n_rech, b_m2.n_val, b_m2.tasa_rech, b_m2.wlo, b_m2.whi, b_m2.temp_alarm_med, b_m2.wlo <= 0.05 ? "PASA" : "FALLA");
        std::printf(" LIMPIO q0 | M3_H0 NIVEL BASE       | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] |         - |   0.000 | %6.2f | %-8s\n",
                    b_m3.n_val, b_m3.k_eval_med, b_m3.n_rech, b_m3.n_val, b_m3.tasa_rech, b_m3.wlo, b_m3.whi, b_m3.temp_alarm_med, b_m3.wlo <= 0.05 ? "PASA" : "FALLA");
        std::printf(" LIMPIO q0 | M1_H1 BASE (Potencia)  | %7zu | %6.1f | Pot = %3zu/%3zu = %.3f          |      1.00 |       - | %6.2f | %-8s\n",
                    b_h1.n_val, b_h1.k_eval_med, b_h1.n_rech_lenta, b_h1.n_val, b_h1.tasa_lenta, b_h1.temp_alarm_med, b_h1.tasa_lenta >= 0.50 ? "PASA" : "FALLA");
        std::printf("-----------+------------------------+---------+--------+--------------------------------+-----------+---------+--------+---------\n");

        bool p5_b4_ok = false;
        for (auto comp : {CompetidorF4::B0_BASE, CompetidorF4::B1_ARGON_VIVO, CompetidorF4::B2_LINEAL_2BORDES,
                          CompetidorF4::B3_ITO_SIN_FILTRO, CompetidorF4::B4_PROPUESTA_UNIFICADA}) {
            const char* nom = (comp == CompetidorF4::B0_BASE)                ? "B0_BASE_ACTUAL"
                            : (comp == CompetidorF4::B1_ARGON_VIVO)          ? "B1_ARGON_VIVO"
                            : (comp == CompetidorF4::B2_LINEAL_2BORDES)     ? "B2_LINEAL_2BORDES"
                            : (comp == CompetidorF4::B3_ITO_SIN_FILTRO)     ? "B3_ITO_SIN_FILTRO"
                                                                             : "B4_RUBIN5_HOLM";
            const auto r1 = correr_lote_mc(m1_h0, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, true, comp);
            const auto r2 = correr_lote_mc(m2_h0, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE,     true, comp);
            const auto r3 = correr_lote_mc(m3_h0, Integracion::NIVEL, true,  ModoNuloLineal::AR_LEVINSON_DURBIN, true, comp);
            const auto h1 = correr_lote_mc(m1_h1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, true, comp);

            const double ret_pot = (b_h1.tasa_lenta > 0.0) ? (h1.tasa_lenta / b_h1.tasa_lenta) : 0.0;
            const double dcens1 = (b_m1.k_eval_med > 0.0) ? std::max(0.0, 1.0 - r1.k_eval_med / b_m1.k_eval_med) : 1.0;
            const double ds1 = std::abs(r1.media_s - b_m1.media_s) / std::max(1e-6, r1.media_ee);
            const double ds2 = std::abs(r2.media_s - b_m2.media_s) / std::max(1e-6, r2.media_ee);
            const double ds3 = std::abs(r3.media_s - b_m3.media_s) / std::max(1e-6, r3.media_ee);

            const bool ok1 = (r1.wlo <= 0.05) && (ret_pot >= 0.85) && (ds1 <= 0.25) && (dcens1 <= cota_cens_m1);
            const bool ok2 = (r2.wlo <= 0.05) && (ds2 <= 0.25);
            const bool ok3 = (r3.wlo <= 0.05) && (ds3 <= 0.25);

            std::printf(" M1 q2 r2  | %-22s | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] | %9.2f | %7.3f | %6.2f | %-8s\n",
                        nom, r1.n_val, r1.k_eval_med, r1.n_rech, r1.n_val, r1.tasa_rech, r1.wlo, r1.whi, ret_pot, ds1, r1.temp_alarm_med, ok1 ? "PASA" : "FALLA");
            std::printf(" M2 q2 r2  | %-22s | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] |         - | %7.3f | %6.2f | %-8s\n",
                        nom, r2.n_val, r2.k_eval_med, r2.n_rech, r2.n_val, r2.tasa_rech, r2.wlo, r2.whi, ds2, r2.temp_alarm_med, ok2 ? "PASA" : "FALLA");
            std::printf(" M3 q2 r2  | %-22s | %7zu | %6.1f | %3zu/%3zu = %.3f [%.3f, %.3f] |         - | %7.3f | %6.2f | %-8s\n",
                        nom, r3.n_val, r3.k_eval_med, r3.n_rech, r3.n_val, r3.tasa_rech, r3.wlo, r3.whi, ds3, r3.temp_alarm_med, ok3 ? "PASA" : "FALLA");

            if (comp == CompetidorF4::B1_ARGON_VIVO || comp == CompetidorF4::B4_PROPUESTA_UNIFICADA) {
                std::printf("  [diag %s] M1: dcens=%.3f (cota=%.3f) vz=%.2f gam=%.4f | M2: vz=%.2f gam=%.4f | M3: vz=%.2f gam=%.4f\n",
                            nom, dcens1, cota_cens_m1, r1.veto_z2_med, r1.media_gamma_rubin,
                            r2.veto_z2_med, r2.media_gamma_rubin, r3.veto_z2_med, r3.media_gamma_rubin);
            }
            if (comp == CompetidorF4::B4_PROPUESTA_UNIFICADA) {
                p5_b4_ok = ok1 && ok2 && ok3;
            }
        }
        reportar("CF4.P5", "Bloque ciego R=200 bajo Rubin M=5 + ark::holm() en M1_H0, M2_H0, M3_H0 y M1_H1",
                 (b_m1.wlo <= 0.05) && (b_m2.wlo <= 0.05) && (b_m3.wlo <= 0.05) && (b_h1.tasa_lenta >= 0.50) && p5_b4_ok);
    }

    std::printf("====================================================================================================\n");
    std::printf(" FALLAS FASE 4: %d\n", fallas);
    std::printf("====================================================================================================\n");
    return fallas == 0 ? 0 : 1;
}
