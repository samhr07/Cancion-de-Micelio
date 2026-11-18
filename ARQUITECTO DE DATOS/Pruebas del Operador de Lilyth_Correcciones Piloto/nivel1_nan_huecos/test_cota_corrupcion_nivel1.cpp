// ============================================================================
// test_cota_corrupcion_nivel1.cpp
// Capa 3 (Batería Unificada Nivel 1 corregida según la 2a revisión de Opus):
//   - CN1.1: Compara el volcado hexadecimal (%a) en M1, M2 y M3 (limpios y con huecos)
//     contra volcado_orig.txt generado por observador_lilyth_original.hpp (identidad bit a bit).
//   - CN1.2: Tamaño nominal bajo H0 en q = 0 con parámetros sellados (a=8.0, CUANTIL_RESOLUCION).
//   - CN1.3: Supervivencia en M3 bajo MERCADO_QUIETO (q=2%) y corte en FRONTERA_ESTRUCTURAL.
//   - CN1.4: Escenario sellado (25% de eventos con Z2 corrupto r=3 y 50% con corte en Z3):
//     verifica exclusión limpia de ambos y caza tanto MN1.6 como MN1.3.
//   - CN1.5: Derivación de cotas con parada monótona estricta (break al 1er fallo) y
//     verificación exclusiva de invariantes estructurales a priori (sin números post-hoc).
// ============================================================================

#include "protocolo_inputs_lilyth.hpp"
#include "banco_sintetico_comun.hpp"

#include <ark/dataset/tipos_dataset.hpp>
#include <ark/lilyth/observador_lilyth.hpp>
#include <ark/lilyth/nulo_lilyth.hpp>
#include <ark/lilyth/modelo_lineal_nulo.hpp>
#include <ark/metrologia/whittle_exacto.hpp>
#include <ark/holm.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;
using namespace ark::lilyth::piloto;

constexpr std::size_t T_MAX_U = 64;
using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, T_MAX_U>;
using ClaveEvento = std::tuple<std::int64_t, int, int>;

inline std::pair<double, double> wilson_95(std::size_t k, std::size_t n) noexcept {
    if (n == 0) return {0.0, 1.0};
    const double z = 1.959963984540054;
    const double p = static_cast<double>(k) / static_cast<double>(n);
    const double d = 1.0 + z * z / static_cast<double>(n);
    const double c = (p + z * z / (2.0 * static_cast<double>(n))) / d;
    const double h = z * std::sqrt(p * (1.0 - p) / static_cast<double>(n) + z * z / (4.0 * static_cast<double>(n) * static_cast<double>(n))) / d;
    return {std::max(0.0, c - h), std::min(1.0, c + h)};
}

inline void en_paralelo(std::size_t n_tareas, std::size_t n_hilos, const std::function<void(std::size_t)>& fn) {
    std::atomic<std::size_t> sig{0};
    std::vector<std::thread> trabajadores;
    const std::size_t h_ef = std::max<std::size_t>(1, n_hilos);
    trabajadores.reserve(h_ef);
    for (std::size_t h = 0; h < h_ef; ++h) {
        trabajadores.emplace_back([&]() {
            for (;;) {
                const std::size_t i = sig.fetch_add(1, std::memory_order_relaxed);
                if (i >= n_tareas) break;
                fn(i);
            }
        });
    }
    for (auto& t : trabajadores) t.join();
}

static std::string hx(double x) {
    char b[64];
    std::snprintf(b, sizeof b, "%a", x);
    return b;
}

static void volcar_contraste_str(std::ostringstream& os, const char* tag, const ContrasteResumenNulo& c) {
    char b[512];
    std::snprintf(b, sizeof b, "%s val=%d K=%zu M=%zu B=%zu S=%s ee=%s p=%s md=%s mu=%s pmin=%s res=%d\n",
                  tag, int(c.valido), c.n_eventos, c.n_nulos_pool, c.b_mc,
                  hx(c.estadistico_s).c_str(), hx(c.ee_nulo_s).c_str(),
                  hx(c.p_bilateral).c_str(), hx(c.mediana_delta_tau).c_str(),
                  hx(c.media_u_pit).c_str(), hx(c.p_min_alcanzable).c_str(), int(c.resoluble));
    os << b;
}

static void volcar_ref_str(std::ostringstream& os, const char* tag, const ReferenciaNuloEvento& r) {
    char b[512];
    std::snprintf(b, sizeof b, "%s ev=%d iny=%d na=%zu nb=%zu nu=%zu int=%zu pdet=%s u=%s dt=%s md=%s;",
                  tag, int(r.evaluado), int(r.nulo_por_inyeccion), r.n_nulos_a, r.n_nulos_b,
                  r.n_nulos_usados, r.intentos_b, hx(r.p_det_inicial).c_str(),
                  hx(r.u_pit).c_str(), hx(r.delta_tau).c_str(), hx(r.mediana_nulo).c_str());
    os << b;
}

static void correr_y_volcar_str(
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

struct EventoMedido {
    std::int64_t t_cambio = -1;
    DesenlaceRanura desenlace = DesenlaceRanura::INACTIVA;
    bool evaluado_ons = false;
    double tau_rest = NaN;
    double p_det_ini = NaN;
    std::vector<double> e_iny;
};

struct ResumenCorrida {
    std::map<ClaveEvento, EventoMedido> eventos;
    std::size_t k_medibles = 0;
    std::size_t k_evaluados_ons = 0;
    std::size_t n_cens_invalido = 0;
    std::size_t n_no_ancla = 0;
    double s_onsager = NaN;
    double ee_s_onsager = NaN;
    double p_holm_onsager = NaN;
    double sd_innov = 1.0;
};

ResumenCorrida ejecutar_corrida_operador(
    const ObsDelta& obs_d_cal,
    const ConfiguracionObservador& cfg_d,
    const ConfiguracionLilyth& cfg_l,
    const SerieSinteticaCoordinada& mundo,
    const std::vector<ObsDelta::Salida>& salidas_delta,
    const std::vector<double>& x_op,
    const std::vector<std::uint8_t>& val_op) {

    const bool usar_purga = cfg_d.purga;
    ObsLil lil{};
    lil.calibrar(obs_d_cal, cfg_l, cfg_d, mundo.x_cal, {},
                 usar_purga ? std::span<const std::size_t>(mundo.fase_cal) : std::span<const std::size_t>{},
                 usar_purga ? std::span<const std::int64_t>(mundo.dia_cal) : std::span<const std::int64_t>{});

    for (std::size_t t = 0; t < x_op.size(); ++t) {
        const bool ok = (val_op[t] != 0u) && std::isfinite(x_op[t]);
        const std::size_t fa = usar_purga ? mundo.fase_op[t] : 0;
        const std::int64_t di = usar_purga ? mundo.dia_op[t] : 0;
        (void)lil.paso(salidas_delta[t], ok ? x_op[t] : NaN, ok, NaN, fa, di);
    }

    lil.completar_nulos();

    ResumenCorrida res{};
    const auto& fn = usar_purga ? lil.ficha_nulo_purgada() : lil.ficha_nulo_cruda();
    const auto& fa = usar_purga ? lil.resumen_armonico_purgada() : lil.resumen_armonico_cruda();
    const auto& mod = fn.modelo_lineal;
    res.sd_innov = std::sqrt(std::max(1e-15, mod.sigma2_innov));
    res.s_onsager = fn.contraste_onsager.estadistico_s;
    res.ee_s_onsager = fn.contraste_onsager.ee_nulo_s;
    res.p_holm_onsager = fa.p_holm[1];

    for (const auto& r : lil.cerradas_historico()) {
        const auto& rama = usar_purga ? r.purgada : r.cruda;
        if (rama.desenlace == DesenlaceRanura::RESTITUIDA || rama.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) {
            ++res.k_medibles;
        } else if (rama.desenlace == DesenlaceRanura::CENSURADA_DATO_INVALIDO) {
            ++res.n_cens_invalido;
        } else if (rama.desenlace == DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE) {
            ++res.n_no_ancla;
        }
        if (rama.nulo_onsager.evaluado) ++res.k_evaluados_ons;

        EventoMedido ev{};
        ev.t_cambio = r.t_cambio;
        ev.desenlace = rama.desenlace;
        ev.evaluado_ons = rama.nulo_onsager.evaluado;
        ev.tau_rest = rama.tau_rest_consolidado;
        ev.p_det_ini = rama.nulo_onsager.p_det_inicial;

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
        res.eventos[{r.t_apertura, r.escala_origen, r.signo}] = std::move(ev);
    }
    return res;
}

struct ComparacionCelda {
    double delta_censura = 0.0;
    double rmse_e_iny_rel = 0.0;
    double caida_p_det = 0.0;
    double sesgo_tau_rest = 0.0;
    double ds_sobre_ee = 0.0;
    bool pasa_estricto = false;
    bool pasa_dulce = false;
};

ComparacionCelda comparar_contra_base(
    const ResumenCorrida& ref,
    const ResumenCorrida& cont,
    bool evaluar_e_iny_y_pdet,
    bool evaluar_tau_rest) {

    ComparacionCelda out{};
    std::size_t k_emp = 0;
    double sum_err2 = 0.0;
    std::size_t n_e = 0;
    double sum_pdet0 = 0.0, sum_pdetc = 0.0;
    std::size_t n_p0 = 0, n_pc = 0;
    double sum_dtau = 0.0;
    std::size_t n_tau = 0;

    for (const auto& [cl, ev0] : ref.eventos) {
        if (ev0.desenlace != DesenlaceRanura::RESTITUIDA && ev0.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
        if (std::isfinite(ev0.p_det_ini)) { sum_pdet0 += ev0.p_det_ini; ++n_p0; }

        auto it = cont.eventos.find(cl);
        if (it == cont.eventos.end()) continue;
        const auto& evc = it->second;

        if (evc.desenlace == DesenlaceRanura::RESTITUIDA || evc.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) {
            ++k_emp;
        }
        if (std::isfinite(evc.p_det_ini)) { sum_pdetc += evc.p_det_ini; ++n_pc; }

        if (!ev0.e_iny.empty() && !evc.e_iny.empty()) {
            const std::size_t m = std::max(ev0.e_iny.size(), evc.e_iny.size());
            for (std::size_t k = 0; k < m; ++k) {
                const double v0 = (k < ev0.e_iny.size()) ? ev0.e_iny[ev0.e_iny.size() - 1 - k] : 0.0;
                const double vc = (k < evc.e_iny.size()) ? evc.e_iny[evc.e_iny.size() - 1 - k] : 0.0;
                const double d = (vc - v0) / ref.sd_innov;
                sum_err2 += d * d;
                ++n_e;
            }
        }
        if (std::isfinite(ev0.tau_rest) && std::isfinite(evc.tau_rest)) {
            sum_dtau += (evc.tau_rest - ev0.tau_rest);
            ++n_tau;
        }
    }

    out.delta_censura = (ref.k_medibles > 0)
                            ? (1.0 - static_cast<double>(k_emp) / static_cast<double>(ref.k_medibles))
                            : 0.0;
    out.rmse_e_iny_rel = (n_e > 0) ? std::sqrt(sum_err2 / static_cast<double>(n_e)) : 0.0;
    const double p0 = (n_p0 > 0) ? (sum_pdet0 / static_cast<double>(n_p0)) : 0.0;
    const double pc = (n_pc > 0) ? (sum_pdetc / static_cast<double>(n_pc)) : 0.0;
    out.caida_p_det = std::max(0.0, p0 - pc);
    out.sesgo_tau_rest = (n_tau > 0) ? (sum_dtau / static_cast<double>(n_tau)) : NaN;

    const double ee = (std::isfinite(ref.ee_s_onsager) && ref.ee_s_onsager > 1e-12) ? ref.ee_s_onsager : 0.05;
    out.ds_sobre_ee = (std::isfinite(cont.s_onsager) && std::isfinite(ref.s_onsager))
                          ? (std::abs(cont.s_onsager - ref.s_onsager) / ee)
                          : 99.0;

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

int main(int argc, char** argv) {
    const std::size_t n_replicas_h0 = (argc > 1) ? static_cast<std::size_t>(std::stoul(argv[1])) : 80;
    const std::string ruta_volcado_orig = (argc > 2) ? argv[2] : "nivel1_nan_huecos/volcado_orig.txt";
    const std::size_t n_hilos = std::max<std::size_t>(1, std::thread::hardware_concurrency());

    std::printf("====================================================================================================\n");
    std::printf(" BATERIA UNIFICADA NIVEL 1 (CAPA 3) - EXPECTATIVAS SELLADAS (54aa457e...) + ANEXO D1-D5\n");
    std::printf(" Replicas H0 = %zu | Hilos = %zu | Regimenes: ESTRICTO (alfa=0.05) y DULCE (<50%% error)\n",
                n_replicas_h0, n_hilos);
    std::printf("====================================================================================================\n");

    int fallas_condiciones = 0;
    auto reportar_condicion = [&](const char* id, const char* desc, bool ok) {
        std::printf("[%s] condicion %s: %s\n", ok ? "PASA" : "FALLA", id, desc);
        if (!ok) ++fallas_condiciones;
    };

    // ------------------------------------------------------------------------
    // CN1.1 REAL: Volcado hexadecimal bit a bit contra observador_lilyth_original.hpp
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
        // Normalizar CRLF a LF para comparar contenido exacto
        ref_hex.erase(std::remove(ref_hex.begin(), ref_hex.end(), '\r'), ref_hex.end());

        const bool cn1_1_ok = !ref_hex.empty() && (actual_hex == ref_hex);
        reportar_condicion("CN1.1", "Volcado hexadecimal (%a) en M1, M2 y M3 (limpios y con huecos) identico bit a bit a observador_lilyth_original.hpp", cn1_1_ok);
    }

    // ------------------------------------------------------------------------
    // Preparar mundos base M2 (ARFIMA) y M3 (NIVEL + PurgaDiurna)
    // ------------------------------------------------------------------------
    ParametrosBancoComun params{};
    params.n_op = 6000;
    const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, params, 202);
    const auto m3 = generar_mundo_m3_nivel_purga(0.65, params, 303);

    ConfiguracionObservador cfg_d_m2{};
    cfg_d_m2.integracion = Integracion::FLUJO;
    cfg_d_m2.agregacion = AgregacionEscala::SUMA;
    cfg_d_m2.arl0_pasos_nativos = 450.0;
    cfg_d_m2.alfa = ALFA_PROYECTO;
    cfg_d_m2.medir_tau50 = false;
    cfg_d_m2.escalas_max = 3;

    ObsDelta obs_d_m2{};
    obs_d_m2.calibrar(cfg_d_m2, m2.x_cal);
    obs_d_m2.reiniciar_flujo();
    std::vector<ObsDelta::Salida> salidas_d_m2(m2.x_op_limpia.size());
    for (std::size_t t = 0; t < m2.x_op_limpia.size(); ++t) {
        salidas_d_m2[t] = obs_d_m2.paso(m2.x_op_limpia[t], NaN, true, 0, 0);
    }

    const std::size_t m_freq = std::max<std::size_t>(32, static_cast<std::size_t>(std::pow(m2.x_cal.size(), 0.65)));
    const auto est_d = elw::estimar_calibracion(m2.x_cal.data(), m2.x_cal.size(), m_freq, -0.45, 0.45);

    ConfiguracionLilyth cfg_l_m2{};
    cfg_l_m2.alfa = ALFA_PROYECTO;
    cfg_l_m2.representacion = ModoRepresentacionLilyth::CRUDA;
    cfg_l_m2.b_max_nulo = 16;
    cfg_l_m2.hilos_nulo = 4;
    cfg_l_m2.semilla_nulo = 20261117ULL;
    cfg_l_m2.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
    cfg_l_m2.d_whittle_cruda = est_d.valido ? est_d.d : 0.25;
    cfg_l_m2.forzar_via_b = true;
    cfg_l_m2.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cfg_l_m2.b_mc_resumen = 1999;

    const std::vector<std::uint8_t> val_limpia_m2(m2.x_op_limpia.size(), 1u);
    const auto ref_m2 = ejecutar_corrida_operador(obs_d_m2, cfg_d_m2, cfg_l_m2, m2, salidas_d_m2, m2.x_op_limpia, val_limpia_m2);

    // ------------------------------------------------------------------------
    // CN1.3 y Derivación MONÓTONA de q_max para MERCADO_QUIETO en M3
    // ------------------------------------------------------------------------
    double q_max_quieto_estricto = 0.0;
    double q_max_quieto_dulce = 0.0;
    {
        ConfiguracionObservador cfg_d_m3{};
        cfg_d_m3.integracion = Integracion::NIVEL;
        cfg_d_m3.agregacion = AgregacionEscala::SUMA;
        cfg_d_m3.arl0_pasos_nativos = 450.0;
        cfg_d_m3.alfa = ALFA_PROYECTO;
        cfg_d_m3.purga = true;
        cfg_d_m3.bins_fase = 24;
        cfg_d_m3.dias_minimos_purga = 20;
        cfg_d_m3.medir_tau50 = false;
        cfg_d_m3.escalas_max = 3;

        ObsDelta obs_d_m3{};
        obs_d_m3.calibrar(cfg_d_m3, m3.x_cal, {}, {}, m3.fase_cal, m3.dia_cal);

        ConfiguracionLilyth cfg_l_m3{};
        cfg_l_m3.alfa = ALFA_PROYECTO;
        cfg_l_m3.representacion = ModoRepresentacionLilyth::PURGADA;
        cfg_l_m3.b_max_nulo = 16;
        cfg_l_m3.hilos_nulo = 4;
        cfg_l_m3.semilla_nulo = 20261117ULL;
        cfg_l_m3.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
        cfg_l_m3.b_mc_resumen = 1999;

        std::vector<std::uint8_t> en_pulso_m3(m3.x_op_limpia.size(), 0u);
        for (std::size_t t0 : m3.tiempos_choque) {
            for (std::size_t k = 0; k < 6 && (t0 + k) < en_pulso_m3.size(); ++k) en_pulso_m3[t0 + k] = 1u;
        }

        bool cn1_3_q2_ok = false;
        bool vivo_est = true, vivo_dul = true;
        for (double q : {0.005, 0.01, 0.02, 0.05, 0.10}) {
            std::mt19937_64 rng(detalle_nulo::splitmix64(20261117ULL ^ static_cast<std::uint64_t>(q * 100000.0)));
            std::uniform_real_distribution<double> u01(0.0, 1.0);

            std::vector<double> x_fis(m3.x_op_limpia.size(), 100.0);
            std::vector<ark::dataset::BarraCausal> barras(m3.x_op_limpia.size());
            double niv = 100.0;
            for (std::size_t t = 0; t < x_fis.size(); ++t) {
                const bool quieta = (t >= 1 && !en_pulso_m3[t] && u01(rng) < q);
                const double dx = (t == 0) ? (m3.x_op_limpia[0] - 100.0) : (m3.x_op_limpia[t] - m3.x_op_limpia[t - 1]);
                niv += quieta ? 0.0 : dx;
                x_fis[t] = niv;
                barras[t].hay_dato = !quieta;
                barras[t].es_imputado = quieta;
            }

            obs_d_m3.reiniciar_flujo();
            std::vector<ObsDelta::Salida> sd_fis(x_fis.size());
            for (std::size_t t = 0; t < x_fis.size(); ++t) {
                sd_fis[t] = obs_d_m3.paso(x_fis[t], NaN, true, m3.fase_op[t], m3.dia_op[t]);
            }
            const std::vector<std::uint8_t> val_todos(x_fis.size(), 1u);
            const auto verdad = ejecutar_corrida_operador(obs_d_m3, cfg_d_m3, cfg_l_m3, m3, sd_fis, x_fis, val_todos);

            std::vector<std::uint8_t> val_act(x_fis.size(), 1u);
            for (std::size_t t = 0; t < x_fis.size(); ++t) val_act[t] = barras[t].es_imputado ? 0u : 1u;
            const auto r_act = ejecutar_corrida_operador(obs_d_m3, cfg_d_m3, cfg_l_m3, m3, sd_fis, x_fis, val_act);

            AdaptadorInputCausalLilyth adapt{};
            adapt.integracion = Integracion::NIVEL;
            adapt.admitir_mercado_quieto = true;
            std::vector<double> x_pro(x_fis.size());
            std::vector<std::uint8_t> val_pro(x_fis.size());
            for (std::size_t t = 0; t < x_fis.size(); ++t) {
                const auto cl = clasificar_barra_dataset(barras[t], x_fis[t]);
                const auto pr = adapt.preparar(x_fis[t], cl, 0.0);
                x_pro[t] = pr.x_para_operador;
                val_pro[t] = pr.valido_para_operador ? 1u : 0u;
            }
            const auto r_pro = ejecutar_corrida_operador(obs_d_m3, cfg_d_m3, cfg_l_m3, m3, sd_fis, x_pro, val_pro);
            const auto cmp_pro = comparar_contra_base(verdad, r_pro, true, true);
            const auto cmp_act = comparar_contra_base(verdad, r_act, true, true);

            if (vivo_est && cmp_pro.pasa_estricto) q_max_quieto_estricto = q; else vivo_est = false;
            if (vivo_dul && cmp_pro.pasa_dulce)    q_max_quieto_dulce = q;    else vivo_dul = false;

            if (std::abs(q - 0.02) < 1e-9) {
                cn1_3_q2_ok = (r_pro.n_cens_invalido == 0) && (cmp_pro.delta_censura == 0.0) && (cmp_act.delta_censura > 0.50);
            }
        }

        ark::dataset::BarraCausal b_finde{};
        b_finde.hay_dato = false;
        b_finde.es_imputado = false;
        b_finde.tipo_hueco = ark::dataset::TipoHueco::FIN_DE_SEMANA;
        AdaptadorInputCausalLilyth ad_f{};
        ad_f.integracion = Integracion::NIVEL;
        (void)ad_f.preparar(100.0, ClasePasoLilyth::VALIDO);
        const auto pf = ad_f.preparar(NaN, clasificar_barra_dataset(b_finde, NaN));
        const bool front_ok = (!pf.valido_para_operador) && (pf.clase == ClasePasoLilyth::FRONTERA_ESTRUCTURAL);

        reportar_condicion("CN1.3", "En M3 (q=2% mercado quieto), PROTOCOLO tiene 0 cortes (dCens=0) y ACTUAL pierde >50%; FRONTERA_ESTRUCTURAL corta limpio",
                           cn1_3_q2_ok && front_ok);
    }

    // ------------------------------------------------------------------------
    // Derivación MONÓTONA (break al 1er fallo) de r_max en Z1, Z2 y Z3
    // ------------------------------------------------------------------------
    std::size_t r_max_z1_borde_estricto = 0, r_max_z1_borde_dulce = 0;
    std::size_t r_max_z1_int_estricto = 0,   r_max_z1_int_dulce = 0;
    std::size_t r_max_z2_estricto = 0,       r_max_z2_dulce = 0;
    std::size_t r_max_z3a_estricto = 0,      r_max_z3a_dulce = 0;
    std::size_t r_max_z3b_estricto = 0,      r_max_z3b_dulce = 0;
    double abs_sesgo_z3a_reloj_r2 = 99.0;
    double abs_sesgo_z3a_comprime_r2 = 0.0;

    for (std::size_t dt : {std::size_t{1}, std::size_t{2}}) {
        bool vivo_e = true, vivo_d = true;
        for (std::size_t r : {std::size_t{1}, std::size_t{2}, std::size_t{3}}) {
            std::vector<std::uint8_t> val_z1 = val_limpia_m2;
            for (const auto& [cl, ev] : ref_m2.eventos) {
                const std::int64_t idx_fin = ev.t_cambio - static_cast<std::int64_t>(dt);
                const std::int64_t idx_ini = idx_fin + 1 - static_cast<std::int64_t>(r);
                for (std::int64_t t = idx_ini; t <= idx_fin; ++t) {
                    if (t >= 0 && static_cast<std::size_t>(t) < val_z1.size()) val_z1[static_cast<std::size_t>(t)] = 0u;
                }
            }
            ConfiguracionLilyth cfg_z1 = cfg_l_m2;
            cfg_z1.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_CAUSAL;
            cfg_z1.r_max_z1_piloto = 5;
            const auto corr = ejecutar_corrida_operador(obs_d_m2, cfg_d_m2, cfg_z1, m2, salidas_d_m2, m2.x_op_limpia, val_z1);
            const auto cmp = comparar_contra_base(ref_m2, corr, true, false);
            if (dt == 1) {
                if (vivo_e && cmp.pasa_estricto) r_max_z1_borde_estricto = r; else vivo_e = false;
                if (vivo_d && cmp.pasa_dulce)    r_max_z1_borde_dulce = r;    else vivo_d = false;
            } else {
                if (vivo_e && cmp.pasa_estricto) r_max_z1_int_estricto = r;   else vivo_e = false;
                if (vivo_d && cmp.pasa_dulce)    r_max_z1_int_dulce = r;      else vivo_d = false;
            }
        }
    }

    {
        bool vivo_e = true, vivo_d = true;
        for (std::size_t r : {std::size_t{1}, std::size_t{2}, std::size_t{3}}) {
            std::vector<std::uint8_t> val_z2 = val_limpia_m2;
            for (const auto& [cl, ev] : ref_m2.eventos) {
                const std::int64_t t_ap = std::get<0>(cl);
                for (std::size_t k = 0; k < r; ++k) {
                    const std::int64_t t = ev.t_cambio + static_cast<std::int64_t>(k);
                    if (t >= 0 && t < t_ap && static_cast<std::size_t>(t) < val_z2.size()) val_z2[static_cast<std::size_t>(t)] = 0u;
                }
            }
            ConfiguracionLilyth cfg_z2 = cfg_l_m2;
            cfg_z2.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::PUENTE_CAUSAL;
            cfg_z2.r_max_z2_piloto = 5;
            const auto corr = ejecutar_corrida_operador(obs_d_m2, cfg_d_m2, cfg_z2, m2, salidas_d_m2, m2.x_op_limpia, val_z2);
            const auto cmp = comparar_contra_base(ref_m2, corr, true, false);
            if (vivo_e && cmp.pasa_estricto) r_max_z2_estricto = r; else vivo_e = false;
            if (vivo_d && cmp.pasa_dulce)    r_max_z2_dulce = r;    else vivo_d = false;
        }
    }

    for (int sub = 0; sub < 2; ++sub) {
        bool vivo_e = true, vivo_d = true;
        for (std::size_t r : {std::size_t{1}, std::size_t{2}, std::size_t{3}}) {
            std::vector<std::uint8_t> val_z3 = val_limpia_m2;
            for (const auto& [cl, ev] : ref_m2.eventos) {
                const std::int64_t t_ap = std::get<0>(cl);
                std::int64_t tau_ini = 1;
                if (sub == 1) {
                    if (!std::isfinite(ev.tau_rest)) continue;
                    tau_ini = static_cast<std::int64_t>(ev.tau_rest);
                } else {
                    if (std::isfinite(ev.tau_rest) && ev.tau_rest <= static_cast<double>(r + 1)) continue;
                }
                for (std::size_t k = 0; k < r; ++k) {
                    const std::int64_t t = t_ap + tau_ini + static_cast<std::int64_t>(k);
                    if (t >= 0 && static_cast<std::size_t>(t) < val_z3.size()) val_z3[static_cast<std::size_t>(t)] = 0u;
                }
            }
            ConfiguracionLilyth cfg_z3 = cfg_l_m2;
            cfg_z3.politica_z3_piloto = ConfiguracionLilyth::PoliticaZ3Piloto::RELOJ_FISICO_TOPE_R;
            cfg_z3.r_max_z3_piloto = 64;
            const auto corr = ejecutar_corrida_operador(obs_d_m2, cfg_d_m2, cfg_z3, m2, salidas_d_m2, m2.x_op_limpia, val_z3);
            const auto cmp = comparar_contra_base(ref_m2, corr, false, true);

            if (sub == 0) {
                if (vivo_e && cmp.pasa_estricto) r_max_z3a_estricto = r; else vivo_e = false;
                if (vivo_d && cmp.pasa_dulce)    r_max_z3a_dulce = r;    else vivo_d = false;
                if (r == 2) {
                    abs_sesgo_z3a_reloj_r2 = std::isfinite(cmp.sesgo_tau_rest) ? std::abs(cmp.sesgo_tau_rest) : 99.0;
                    ConfiguracionLilyth cfg_comp = cfg_l_m2;
                    cfg_comp.politica_z3_piloto = ConfiguracionLilyth::PoliticaZ3Piloto::MUTANTE_COMPRIME_TAU;
                    const auto corr_comp = ejecutar_corrida_operador(obs_d_m2, cfg_d_m2, cfg_comp, m2, salidas_d_m2, m2.x_op_limpia, val_z3);
                    const auto cmp_comp = comparar_contra_base(ref_m2, corr_comp, false, true);
                    abs_sesgo_z3a_comprime_r2 = std::isfinite(cmp_comp.sesgo_tau_rest) ? std::abs(cmp_comp.sesgo_tau_rest) : 0.0;
                }
            } else {
                if (vivo_e && cmp.pasa_estricto) r_max_z3b_estricto = r; else vivo_e = false;
                if (vivo_d && cmp.pasa_dulce)    r_max_z3b_dulce = r;    else vivo_d = false;
            }
        }
    }

    // ------------------------------------------------------------------------
    // CN1.2 y CN1.4 sobre réplicas H0 (M1_H0 con CUANTIL_RESOLUCION)
    // En CN1.4 se corrompe el 25% de eventos en Z2 (r=3) y el 50% en Z3 (r=1)
    // ------------------------------------------------------------------------
    ParametrosReplicaH0H1 p_h0{};
    p_h0.lenta_h1 = false;
    std::atomic<std::size_t> n_val_q0{0}, n_rech_q0{0};
    std::atomic<std::size_t> n_val_z2{0}, n_rech_z2{0};
    std::atomic<std::size_t> n_eval_via_b_z2_corrupta{0}, n_eval_z3_censurada{0};

    en_paralelo(n_replicas_h0, n_hilos, [&](std::size_t i) {
        const auto m_rep = generar_mundo_m1_h0_h1(0.80, p_h0, i);
        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.agregacion = AgregacionEscala::SUMA;
        cd.arl0_pasos_nativos = 450.0;
        cd.alfa = ALFA_PROYECTO;
        cd.medir_tau50 = false;
        cd.semilla = 20261117ULL + 10007ULL * i;
        cd.escalas_max = 2;

        ObsDelta od{};
        od.calibrar(cd, m_rep.x_cal);
        if (od.escalas() == 0) return;

        od.reiniciar_flujo();
        std::vector<ObsDelta::Salida> sd_rep(m_rep.x_op_limpia.size());
        for (std::size_t t = 0; t < m_rep.x_op_limpia.size(); ++t) {
            sd_rep[t] = od.paso(m_rep.x_op_limpia[t], NaN, true, 0, 0);
        }

        ConfiguracionLilyth cl{};
        cl.alfa = ALFA_PROYECTO;
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 32;
        cl.hilos_nulo = 1;
        cl.semilla_nulo = 20261117ULL + 31337ULL * i;
        cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
        cl.pool_nulo = PoolNuloLilyth::DERIVADO;
        cl.horizonte = HorizonteLilyth::CUANTIL_RESOLUCION;
        cl.b_mc_resumen = 1999;

        const std::vector<std::uint8_t> val0(m_rep.x_op_limpia.size(), 1u);
        const auto r0 = ejecutar_corrida_operador(od, cd, cl, m_rep, sd_rep, m_rep.x_op_limpia, val0);
        if (std::isfinite(r0.p_holm_onsager)) {
            ++n_val_q0;
            if (r0.p_holm_onsager < ALFA_PROYECTO) ++n_rech_q0;
        }

        // CN1.4 (Escenario sellado + D3): 25% de eventos con Z2 corrupto (r=3) y 50% con corte en Z3 (t_ap + 1)
        std::vector<std::uint8_t> val_z2 = val0;
        std::vector<ClaveEvento> claves_z2_corruptas;
        std::vector<ClaveEvento> claves_z3_cortadas;
        std::size_t idx_ev = 0;
        for (const auto& [k_ev, ev] : r0.eventos) {
            const std::int64_t t_ap = std::get<0>(k_ev);
            if ((idx_ev % 4 == 0) && ev.t_cambio >= 0 && ev.t_cambio < t_ap) {
                for (std::int64_t k = 0; k < 3 && (ev.t_cambio + k) < t_ap; ++k) {
                    val_z2[static_cast<std::size_t>(ev.t_cambio + k)] = 0u;
                }
                claves_z2_corruptas.push_back(k_ev);
            } else if ((idx_ev % 2 == 1) && (t_ap + 1 < static_cast<std::int64_t>(val_z2.size()))) {
                // Cortar en Z3 en t_ap + 1 (antes de restituir si tau_rest > 1)
                if (!std::isfinite(ev.tau_rest) || ev.tau_rest > 1.0) {
                    val_z2[static_cast<std::size_t>(t_ap + 1)] = 0u;
                    claves_z3_cortadas.push_back(k_ev);
                }
            }
            ++idx_ev;
        }

        ConfiguracionLilyth cl_z2 = cl;
        cl_z2.forzar_via_b = true;
        cl_z2.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cl_z2.r_max_z2_piloto = 0;
        const auto rz2 = ejecutar_corrida_operador(od, cd, cl_z2, m_rep, sd_rep, m_rep.x_op_limpia, val_z2);
        if (std::isfinite(rz2.p_holm_onsager)) {
            ++n_val_z2;
            if (rz2.p_holm_onsager < ALFA_PROYECTO) ++n_rech_z2;
        }
        for (const auto& kz2 : claves_z2_corruptas) {
            auto it = rz2.eventos.find(kz2);
            if (it != rz2.eventos.end() && it->second.evaluado_ons) ++n_eval_via_b_z2_corrupta;
        }
        for (const auto& kz3 : claves_z3_cortadas) {
            auto it = rz2.eventos.find(kz3);
            if (it != rz2.eventos.end() && it->second.desenlace == DesenlaceRanura::CENSURADA_DATO_INVALIDO && it->second.evaluado_ons) {
                ++n_eval_z3_censurada;
            }
        }
    });

    const auto [w_lo_q0, w_hi_q0] = wilson_95(n_rech_q0.load(), n_val_q0.load());
    const auto [w_lo_z2, w_hi_z2] = wilson_95(n_rech_z2.load(), n_val_z2.load());

    reportar_condicion("CN1.2", ("Tamano nominal H0 en q=0: " + std::to_string(n_rech_q0.load()) + "/" +
                                 std::to_string(n_val_q0.load()) + " [Wilson lo=" + std::to_string(w_lo_q0) + " <= 0.05]").c_str(),
                       (n_val_q0.load() > 0 && w_lo_q0 <= ALFA_PROYECTO));

    reportar_condicion("CN1.4", ("Exclusion limpia de 25% Z2 corrupto (eval_z2=" + std::to_string(n_eval_via_b_z2_corrupta.load()) +
                                 ") y 50% Z3 censurado (eval_z3=" + std::to_string(n_eval_z3_censurada.load()) +
                                 ") con Wilson lo=" + std::to_string(w_lo_z2) + " <= 0.05").c_str(),
                       (n_eval_via_b_z2_corrupta.load() == 0 && n_eval_z3_censurada.load() == 0 &&
                        n_val_z2.load() > 0 && w_lo_z2 <= ALFA_PROYECTO));

    // CN1.5: Solo invariantes estructurales a priori (monotonía ESTRICTO <= DULCE y reloj físico con menor sesgo que comprimir tau)
    const bool monotonia_regimenes = (q_max_quieto_estricto <= q_max_quieto_dulce) &&
                                     (r_max_z1_borde_estricto <= r_max_z1_borde_dulce) &&
                                     (r_max_z1_int_estricto <= r_max_z1_int_dulce) &&
                                     (r_max_z2_estricto <= r_max_z2_dulce) &&
                                     (r_max_z3a_estricto <= r_max_z3a_dulce) &&
                                     (r_max_z3b_estricto <= r_max_z3b_dulce);
    const bool reloj_mejor_que_comprimir = (abs_sesgo_z3a_reloj_r2 < abs_sesgo_z3a_comprime_r2);
    const bool cn1_5_ok = monotonia_regimenes && reloj_mejor_que_comprimir;

    reportar_condicion("CN1.5", "Invariantes estructurales a priori (cota_ESTRICTO <= cota_DULCE en toda zona y |sesgo_reloj| < |sesgo_comprime| en Z3a)", cn1_5_ok);

    std::printf("\n====================================================================================================\n");
    std::printf(" COTAS MONOTONAS DERIVADAS DESDE LA MALLA MEDIDA (r in {1,2,3}, q in {0.5%%..10%%}):\n");
    std::printf("  1. MERCADO_QUIETO (M3 trayectoria):       q_max_ESTRICTO = %.1f%% | q_max_DULCE = %.1f%%\n",
                q_max_quieto_estricto * 100.0, q_max_quieto_dulce * 100.0);
    std::printf("  2. ZONA Z1a (Borde inmediato t_cambio-1): r_max_ESTRICTO = %zu     | r_max_DULCE = %zu (malla r<=3)\n",
                r_max_z1_borde_estricto, r_max_z1_borde_dulce);
    std::printf("  3. ZONA Z1b (Interior t <= t_cambio-2):   r_max_ESTRICTO = %zu     | r_max_DULCE = %zu (malla r<=3)\n",
                r_max_z1_int_estricto, r_max_z1_int_dulce);
    std::printf("  4. ZONA Z2  (Perfil choque [t_cambio,t_a]):r_max_ESTRICTO = %zu     | r_max_DULCE = %zu (malla r<=3)\n",
                r_max_z2_estricto, r_max_z2_dulce);
    std::printf("  5. ZONA Z3a (Avanzar fuera de tau_rest):  r_max_ESTRICTO = %zu     | r_max_DULCE = %zu (malla r<=3)\n",
                r_max_z3a_estricto, r_max_z3a_dulce);
    std::printf("  6. ZONA Z3b (Avanzar encima de tau_rest): r_max_ESTRICTO = %zu     | r_max_DULCE = %zu (malla r<=3)\n",
                r_max_z3b_estricto, r_max_z3b_dulce);
    std::printf("====================================================================================================\n");
    std::printf(" FALLAS DE CONDICIONES: %d\n", fallas_condiciones);
    std::printf("====================================================================================================\n");

    return fallas_condiciones == 0 ? 0 : 1;
}
