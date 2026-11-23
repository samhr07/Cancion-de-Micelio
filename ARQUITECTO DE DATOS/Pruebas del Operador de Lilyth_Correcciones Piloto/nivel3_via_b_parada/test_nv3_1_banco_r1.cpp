// ============================================================================
// test_nv3_1_banco_r1.cpp -- Control sellado NV3.1 (T3.5):
//   Banco R1: AR(1) phi=0.80 estacionario, 40 semillas ciegas,
//   forzar_via_b = true, usar_enmienda_nivel2 = true, usar_enmienda_nivel3 = true.
// Semilla ciega exclusiva: SEMILLA_CIEGA_NIVEL3 = 20261124ULL.
// ============================================================================

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <random>
#include <string_view>
#include <thread>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;

using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, 64>;

constexpr std::uint64_t SEMILLA_CIEGA_NIVEL3 = 20261124ULL;
constexpr std::size_t N_SEM = 40;
constexpr std::size_t N_CAL = 4096;
constexpr std::size_t N_OP = 6144;
constexpr double PHI_R1 = 0.80;
constexpr double SD_EPS_R1 = 1.0;

std::uint64_t sem_cal(std::size_t s) noexcept {
    return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL3 ^ (0x9E3779B97F4A7C15ULL * (s + 1ULL)));
}

std::uint64_t sem_op(std::size_t s) noexcept {
    return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL3 ^ (0xBF58476D1CE4E5B9ULL * (s + 101ULL)));
}

std::uint64_t sem_nulo(std::size_t s) noexcept {
    return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL3 ^ (0x94D049BB133111EBULL * (s + 211ULL)));
}

std::vector<double> generar_ar1_estacionario(
    std::size_t n, double phi, double sd_eps, std::uint64_t semilla) {
    std::vector<double> x(n, 0.0);
    if (n == 0) return x;
    std::mt19937_64 rng(semilla);
    std::normal_distribution<double> n01(0.0, 1.0);
    const double std_est = sd_eps / std::sqrt(std::max(1e-12, 1.0 - phi * phi));
    double prev = std_est * n01(rng);
    x[0] = prev;
    for (std::size_t t = 1; t < n; ++t) {
        prev = phi * prev + sd_eps * n01(rng);
        x[t] = prev;
    }
    return x;
}

struct ResultadoSemilla {
    std::size_t semilla_idx = 0;
    bool evaluada = false;
    double estadistico_s = 0.0;
    double p_bilateral = 1.0;
    bool rechazada = false;
    std::size_t n_eventos_contraste = 0;
    std::size_t m_min = 0;
    std::size_t t_max_ef = 0;
    std::size_t k_med = 0;
    std::size_t k_pdet_pos = 0;
    std::size_t k_eval = 0;
    std::size_t n_nulos_b = 0;
    std::size_t intentos_b = 0;
    std::size_t paradas_tempranas_b = 0;
    std::size_t solapes_b = 0;
    std::size_t con_ajena_previa_b = 0;
    std::size_t con_reapertura_posterior_b = 0;
    std::size_t escala_coincide_b = 0;
    bool calibracion_ok = false;
};

ResultadoSemilla ejecutar_semilla(std::size_t s) {
    ResultadoSemilla res{};
    res.semilla_idx = s;

    const auto x_cal = generar_ar1_estacionario(N_CAL, PHI_R1, SD_EPS_R1, sem_cal(s));
    const auto x_op = generar_ar1_estacionario(N_OP, PHI_R1, SD_EPS_R1, sem_op(s));

    ConfiguracionObservador cd{};
    cd.integracion = Integracion::FLUJO;
    cd.agregacion = AgregacionEscala::SUMA;
    cd.arl0_pasos_nativos = 450.0;
    cd.alfa = ALFA_PROYECTO;
    cd.purga = false;
    cd.bins_fase = 1;
    cd.medir_tau50 = false;
    cd.escalas_max = 3;
    cd.semilla = sem_cal(s);

    auto od = std::make_unique<ObsDelta>();
    od->calibrar(cd, x_cal);

    ConfiguracionLilyth cl{};
    cl.alfa = ALFA_PROYECTO;
    cl.representacion = ModoRepresentacionLilyth::CRUDA;
    cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
    cl.modo_armonico = ModoArmonico::DESACTIVADO;
    cl.b_max_nulo = 16;
    cl.hilos_nulo = 1;
    cl.b_mc_resumen = 1999;
    cl.semilla_nulo = sem_nulo(s);
    cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cl.forzar_via_b = true;
    cl.usar_enmienda_nivel2 = true;
    cl.usar_enmienda_nivel3 = true;

    auto lil = std::make_unique<ObsLil>();
    const bool cal_lil = lil->calibrar(*od, cl, cd, x_cal);
    res.calibracion_ok = cal_lil;
    if (!cal_lil) return res;

    od->reiniciar_flujo();
    lil->reiniciar_flujo_conservando_purga();

    for (std::size_t t = 0; t < N_OP; ++t) {
        const auto sd = od->paso(x_op[t], NaN, true, 0, 0);
        (void)lil->paso(sd, x_op[t], true, NaN, 0, 0);
    }
    lil->completar_nulos();

    const auto& fn = lil->ficha_nulo_cruda();
    res.m_min = fn.m_min;
    res.t_max_ef = lil->configuracion().t_max_efectivo;
    res.n_eventos_contraste = fn.contraste_onsager.n_eventos;

    if (fn.contraste_onsager.valido &&
        std::isfinite(fn.contraste_onsager.estadistico_s) &&
        std::isfinite(fn.contraste_onsager.p_bilateral)) {
        res.evaluada = true;
        res.estadistico_s = fn.contraste_onsager.estadistico_s;
        res.p_bilateral = fn.contraste_onsager.p_bilateral;
        res.rechazada = (res.p_bilateral < 0.05);
    }

    for (const auto& r : lil->cerradas_historico()) {
        const auto d = r.cruda.desenlace;
        if (d != DesenlaceRanura::RESTITUIDA && d != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
        ++res.k_med;
        const auto& ref = r.cruda.nulo_onsager;
        if (std::isfinite(ref.p_det_inicial) && ref.p_det_inicial > 0.0) ++res.k_pdet_pos;
        if (ref.evaluado && ref.n_nulos_b >= fn.m_min) ++res.k_eval;
        res.n_nulos_b += ref.n_nulos_b;
        res.intentos_b += ref.intentos_b;
        res.paradas_tempranas_b += ref.paradas_tempranas_b;
        res.solapes_b += ref.solapes_b;
        res.con_ajena_previa_b += ref.con_ajena_previa_b;
        res.con_reapertura_posterior_b += ref.con_reapertura_posterior_b;
        res.escala_coincide_b += ref.escala_coincide_b;
    }

    return res;
}

std::size_t parsear_hilos(int argc, char** argv) noexcept {
    const unsigned hw = std::thread::hardware_concurrency();
    std::size_t hilos = hw > 0 ? static_cast<std::size_t>(hw) : 4u;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        constexpr std::string_view pref = "--hilos=";
        if (arg.starts_with(pref)) {
            const std::string_view val = arg.substr(pref.size());
            std::size_t parsed = 0;
            bool ok = !val.empty();
            for (char c : val) {
                if (c < '0' || c > '9') {
                    ok = false;
                    break;
                }
                parsed = parsed * 10u + static_cast<std::size_t>(c - '0');
            }
            if (ok && parsed > 0) hilos = parsed;
        }
    }
    return std::clamp<std::size_t>(hilos, 1u, N_SEM);
}

}  // namespace

int main(int argc, char** argv) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t hilos = parsear_hilos(argc, argv);

    std::printf("================================================================================\n");
    std::printf("NV3.1 (T3.5) -- Control Sellado Banco R1 (AR(1) phi=%.2f, forzar_via_b=true)\n", PHI_R1);
    std::printf("  SEMILLA_CIEGA_NIVEL3 = %lluULL | N_SEM = %zu | N_CAL = %zu | N_OP = %zu | hilos = %zu\n",
                static_cast<unsigned long long>(SEMILLA_CIEGA_NIVEL3), N_SEM, N_CAL, N_OP, hilos);
    std::printf("================================================================================\n");

    std::vector<ResultadoSemilla> resultados(N_SEM);
    std::atomic<std::size_t> siguiente{0};

    std::vector<std::thread> trabajadores;
    trabajadores.reserve(hilos);
    for (std::size_t w = 0; w < hilos; ++w) {
        trabajadores.emplace_back([&]() {
            while (true) {
                const std::size_t s = siguiente.fetch_add(1u, std::memory_order_relaxed);
                if (s >= N_SEM) break;
                resultados[s] = ejecutar_semilla(s);
            }
        });
    }
    for (auto& th : trabajadores) th.join();

    const auto t1 = std::chrono::steady_clock::now();
    const double segs = std::chrono::duration<double>(t1 - t0).count();

    std::printf("\n%3s %3s %4s %5s %4s %8s %7s %4s | %4s %4s %4s | %5s %5s %4s %4s %4s %4s %5s\n",
                "sem", "cal", "eval", "t_ef", "mmin", "S", "p_bilat", "rech",
                "kmed", "kpdt", "kevl",
                "nulB", "intB", "parB", "solB", "ajPr", "reap", "escC");
    std::printf("---------------------------------------------------------------------------------------------------------\n");

    std::size_t n_cal_ok = 0;
    std::size_t n_sem_eval = 0;
    std::size_t n_rech = 0;
    std::size_t k_med_total = 0;
    std::size_t k_pdet_pos_total = 0;
    std::size_t k_eval_total = 0;
    std::size_t n_nulos_b_total = 0;
    std::size_t intentos_b_total = 0;
    std::size_t paradas_tempranas_b_total = 0;
    std::size_t solapes_b_total = 0;
    std::size_t con_ajena_previa_b_total = 0;
    std::size_t con_reapertura_posterior_b_total = 0;
    std::size_t escala_coincide_b_total = 0;
    std::size_t n_ev_contraste_total = 0;

    std::vector<double> vals_s;
    vals_s.reserve(N_SEM);

    for (const auto& r : resultados) {
        if (r.calibracion_ok) ++n_cal_ok;
        if (r.evaluada) {
            ++n_sem_eval;
            vals_s.push_back(r.estadistico_s);
            if (r.rechazada) ++n_rech;
        }
        k_med_total += r.k_med;
        k_pdet_pos_total += r.k_pdet_pos;
        k_eval_total += r.k_eval;
        n_nulos_b_total += r.n_nulos_b;
        intentos_b_total += r.intentos_b;
        paradas_tempranas_b_total += r.paradas_tempranas_b;
        solapes_b_total += r.solapes_b;
        con_ajena_previa_b_total += r.con_ajena_previa_b;
        con_reapertura_posterior_b_total += r.con_reapertura_posterior_b;
        escala_coincide_b_total += r.escala_coincide_b;
        n_ev_contraste_total += r.n_eventos_contraste;

        std::printf("%3zu %3d %4d %5zu %4zu %+8.4f %7.4f %4d | %4zu %4zu %4zu | %5zu %5zu %4zu %4zu %4zu %4zu %5zu\n",
                    r.semilla_idx,
                    r.calibracion_ok ? 1 : 0,
                    r.evaluada ? 1 : 0,
                    r.t_max_ef,
                    r.m_min,
                    r.estadistico_s,
                    r.p_bilateral,
                    r.rechazada ? 1 : 0,
                    r.k_med,
                    r.k_pdet_pos,
                    r.k_eval,
                    r.n_nulos_b,
                    r.intentos_b,
                    r.paradas_tempranas_b,
                    r.solapes_b,
                    r.con_ajena_previa_b,
                    r.con_reapertura_posterior_b,
                    r.escala_coincide_b);
    }

    double S_medio = 0.0;
    for (double v : vals_s) S_medio += v;
    if (!vals_s.empty()) S_medio /= static_cast<double>(vals_s.size());

    double sum_sq = 0.0;
    for (double v : vals_s) {
        const double d = v - S_medio;
        sum_sq += d * d;
    }
    const double var_S = (vals_s.size() > 1) ? (sum_sq / static_cast<double>(vals_s.size() - 1)) : 0.0;
    const double ee_S = (vals_s.size() > 1) ? std::sqrt(var_S / static_cast<double>(vals_s.size())) : 0.0;
    const double z_S = (ee_S > 0.0) ? (S_medio / ee_S) : 0.0;

    const double tasa_rech = (n_sem_eval > 0)
        ? (static_cast<double>(n_rech) / static_cast<double>(n_sem_eval))
        : 1.0;
    const double p_det_global = (intentos_b_total > 0)
        ? (static_cast<double>(n_nulos_b_total) / static_cast<double>(intentos_b_total))
        : 0.0;
    const double surv_s3 = (k_med_total > 0)
        ? (static_cast<double>(k_pdet_pos_total) / static_cast<double>(k_med_total))
        : 0.0;
    const double cob_s3 = (k_med_total > 0)
        ? (static_cast<double>(k_eval_total) / static_cast<double>(k_med_total))
        : 0.0;
    const double tasa_escala_coincide = (n_nulos_b_total > 0)
        ? (static_cast<double>(escala_coincide_b_total) / static_cast<double>(n_nulos_b_total))
        : 0.0;

    std::printf("\n================================================================================\n");
    std::printf("RESUMEN AGREGADO NV3.1 (tiempo=%.2f s)\n", segs);
    std::printf("  calibradas                  = %zu / %zu\n", n_cal_ok, N_SEM);
    std::printf("  n_sem_eval                  = %zu / %zu  (ev_contraste_total=%zu)\n",
                n_sem_eval, N_SEM, n_ev_contraste_total);
    std::printf("  S_medio                     = %+.6f | ee_S = %.6f | z_S = %+.6f\n", S_medio, ee_S, z_S);
    std::printf("  tasa_rech (p_bilat < 0.05)  = %.6f  (%zu / %zu)\n", tasa_rech, n_rech, n_sem_eval);
    std::printf("  p_det_global                = %.6f  (%zu / %zu intentos_b)\n",
                p_det_global, n_nulos_b_total, intentos_b_total);
    std::printf("  paradas_tempranas_b_total   = %zu  (%.4f%% de intentos_b)\n",
                paradas_tempranas_b_total,
                intentos_b_total > 0
                    ? 100.0 * static_cast<double>(paradas_tempranas_b_total) / static_cast<double>(intentos_b_total)
                    : 0.0);
    std::printf("  surv_s3 (p_det_inicial > 0) = %.6f  (%zu / %zu ranuras)\n",
                surv_s3, k_pdet_pos_total, k_med_total);
    std::printf("  cob_s3 (n_nulos_b >= m_min) = %.6f  (%zu / %zu ranuras)\n",
                cob_s3, k_eval_total, k_med_total);
    std::printf("  solapes_b_total             = %zu\n", solapes_b_total);
    std::printf("  con_ajena_previa_b_total    = %zu\n", con_ajena_previa_b_total);
    std::printf("  con_reapertura_post_b_total = %zu\n", con_reapertura_posterior_b_total);
    std::printf("  Telemetria Nivel 4:\n");
    std::printf("    escala_coincide_b_total   = %zu / %zu  (tasa=%.6f)\n",
                escala_coincide_b_total, n_nulos_b_total, tasa_escala_coincide);
    std::printf("--------------------------------------------------------------------------------\n");

    const bool ok_sem = (n_sem_eval >= 38);
    const bool ok_s = (n_sem_eval >= 38 && ee_S > 0.0 && std::isfinite(z_S) && std::abs(z_S) <= 2.0);
    const bool ok_rec = (n_sem_eval >= 38 && tasa_rech >= 0.0 && tasa_rech <= 0.10);
    const bool ok_pdt = (p_det_global >= 0.40 && paradas_tempranas_b_total > 0);
    const bool ok_cob = (surv_s3 >= 0.95 && cob_s3 >= 0.55);
    const bool ok_sim = (solapes_b_total == 0 && con_ajena_previa_b_total > 0);

    std::size_t subfallas = 0;
    auto reportar_sub = [&](const char* tag, bool ok, const char* criterio) {
        if (!ok) ++subfallas;
        std::printf("  [%s] %-6s : %s\n", ok ? "PASA" : "FALLA", tag, criterio);
    };
    reportar_sub("ok_sem", ok_sem, "n_sem_eval >= 38");
    reportar_sub("ok_s",   ok_s,   "|z_S| <= 2.0");
    reportar_sub("ok_rec", ok_rec, "tasa_rech in [0.0, 0.10]");
    reportar_sub("ok_pdt", ok_pdt, "p_det_global >= 0.40 && paradas_tempranas_b_total > 0");
    reportar_sub("ok_cob", ok_cob, "surv_s3 >= 0.95 && cob_s3 >= 0.55");
    reportar_sub("ok_sim", ok_sim, "solapes_b_total == 0 && con_ajena_previa_b_total > 0");

    const bool ok_nv3_1 = (n_cal_ok == N_SEM) && ok_sem && ok_s && ok_rec && ok_pdt && ok_cob && ok_sim;
    const std::size_t fallas = ok_nv3_1 ? 0u : 1u;
    std::printf("\n[%s] NV3.1 (subfallas=%zu)\n", ok_nv3_1 ? "PASA" : "FALLA", subfallas);
    std::printf("FALLAS NV3.1: %zu\n", fallas);
    std::fflush(stdout);
    return fallas == 0 ? 0 : 1;
}
