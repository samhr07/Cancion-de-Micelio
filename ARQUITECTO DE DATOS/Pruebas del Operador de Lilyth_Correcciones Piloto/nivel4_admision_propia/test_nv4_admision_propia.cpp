// ============================================================================
// test_nv4_admision_propia.cpp -- Bateria Confirmatoria Sellada de Nivel 4
//   Decision formal: "No-Op Demostrado" (P0_NIVEL3_BASE mantenida intacta).
//   Contrato pre-registrado: nivel4_admision_propia/EXPECTATIVAS_NIVEL4.md
//   Controles falsables:
//     NV4.1 : E2E AR(1) phi=0.80 (24 semillas, 3 escalas, completar_nulos)
//     NV4.2 : E2E ARFIMA(1, d=+0.35, 0) (16 semillas, 3 escalas, completar_nulos)
//     NV4.3 : Contraste pareado ciego P0..P3 + verificacion Regla 6
//   Opciones CLI:
//     --hilos=N        : numero de hilos externos para paralelismo por semilla/evento
//     --solo-operador  : omite NV4.3 y ejecuta unicamente NV4.1 y NV4.2 (uso en mutaciones)
//   Semilla ciega exclusiva e inedita: SEMILLA_CIEGA_NIVEL4 = 20261125ULL.
// ============================================================================

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <random>
#include <span>
#include <string_view>
#include <thread>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;

constexpr std::uint64_t SEMILLA_CIEGA_NIVEL4 = 20261125ULL;

struct OpcionesCli {
    std::size_t hilos = 1;
    bool solo_operador = false;
};

OpcionesCli parsear_cli(int argc, char** argv) noexcept {
    OpcionesCli opt{};
    const unsigned hw = std::thread::hardware_concurrency();
    opt.hilos = hw > 0 ? static_cast<std::size_t>(hw) : 4u;
    for (int i = 1; i < argc; ++i) {
        const std::string_view arg(argv[i]);
        constexpr std::string_view pref = "--hilos=";
        if (arg == "--solo-operador") {
            opt.solo_operador = true;
        } else if (arg.starts_with(pref)) {
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
            if (ok && parsed > 0u) opt.hilos = parsed;
        }
    }
    opt.hilos = std::clamp<std::size_t>(opt.hilos, 1u, 64u);
    return opt;
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

std::vector<double> generar_arfima1d0(
    std::size_t n, double phi1, double d, double mu, double sigma_e, std::uint64_t semilla, std::size_t burn = 1024) {
    const std::size_t tot = burn + n;
    std::vector<double> pi_inv(tot, 0.0);
    elw::coeficientes(-d, tot, pi_inv.data());
    std::vector<double> u(tot, 0.0);
    std::vector<double> out(n, 0.0);
    std::mt19937_64 rng(semilla);
    std::normal_distribution<double> nd(0.0, sigma_e);
    double u_prev = 0.0;
    for (std::size_t s = 0; s < tot; ++s) {
        const double u_t = phi1 * u_prev + nd(rng);
        u[s] = u_t;
        u_prev = u_t;
        if (s >= burn) {
            double acc = 0.0;
            for (std::size_t k = 0; k <= s; ++k) acc += pi_inv[k] * u[s - k];
            out[s - burn] = mu + acc;
        }
    }
    return out;
}

struct ResultadoSemillaN4 {
    std::size_t semilla_idx = 0;
    bool calibracion_ok = false;
    bool evaluada = false;
    std::size_t escalas = 0;
    std::size_t m_min = 0;
    std::size_t t_max_ef = 0;
    double estadistico_s = 0.0;
    double p_bilateral = 1.0;
    bool rechazada = false;
    std::size_t n_eventos_contraste = 0;
    std::array<std::size_t, 3> k_med_j{};
    std::size_t k_med = 0;
    std::size_t k_pdet_pos = 0;
    std::size_t k_eval = 0;
    std::size_t n_nulos_b = 0;
    std::size_t intentos_b = 0;
    std::size_t paradas_tempranas_b = 0;
    std::size_t solapes_b = 0;
    std::size_t escala_coincide_b = 0;
    std::size_t n_nulos_b_j12 = 0;
    std::size_t intentos_b_j12 = 0;
    std::size_t escala_coincide_b_j12 = 0;
};

template <class ObsLil>
void acumular_ranuras_cruda(const ObsLil& lil, ResultadoSemillaN4& res) {
    const auto& fn = lil.ficha_nulo_cruda();
    res.m_min = fn.m_min;
    res.t_max_ef = lil.configuracion().t_max_efectivo;
    res.n_eventos_contraste = fn.contraste_onsager.n_eventos;

    if (fn.contraste_onsager.valido &&
        std::isfinite(fn.contraste_onsager.estadistico_s) &&
        std::isfinite(fn.contraste_onsager.p_bilateral)) {
        res.evaluada = true;
        res.estadistico_s = fn.contraste_onsager.estadistico_s;
        res.p_bilateral = fn.contraste_onsager.p_bilateral;
        res.rechazada = (res.p_bilateral < 0.05);
    }

    for (const auto& r : lil.cerradas_historico()) {
        const auto d = r.cruda.desenlace;
        if (d != DesenlaceRanura::RESTITUIDA && d != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
        ++res.k_med;
        if (r.escala_origen >= 0 && r.escala_origen < 3) {
            ++res.k_med_j[static_cast<std::size_t>(r.escala_origen)];
        }
        const auto& ref = r.cruda.nulo_onsager;
        if (std::isfinite(ref.p_det_inicial) && ref.p_det_inicial > 0.0) ++res.k_pdet_pos;
        if (ref.evaluado && ref.n_nulos_b >= fn.m_min) ++res.k_eval;
        res.n_nulos_b += ref.n_nulos_b;
        res.intentos_b += ref.intentos_b;
        res.paradas_tempranas_b += ref.paradas_tempranas_b;
        res.solapes_b += ref.solapes_b;
        res.escala_coincide_b += ref.escala_coincide_b;
        if (r.escala_origen == 1 || r.escala_origen == 2) {
            res.n_nulos_b_j12 += ref.n_nulos_b;
            res.intentos_b_j12 += ref.intentos_b;
            res.escala_coincide_b_j12 += ref.escala_coincide_b;
        }
    }
}

struct ResumenBancoN4 {
    std::size_t n_sem = 0;
    std::size_t n_cal_ok = 0;
    std::size_t min_esc = 0;
    std::size_t n_sem_eval = 0;
    std::size_t n_rech = 0;
    double S_medio = 0.0;
    double ee_S = 0.0;
    double z_S = 0.0;
    std::array<std::size_t, 3> k_med_j{};
    std::size_t k_med = 0;
    std::size_t k_pdet_pos = 0;
    std::size_t k_eval = 0;
    std::size_t n_nulos_b = 0;
    std::size_t intentos_b = 0;
    std::size_t paradas_tempranas_b = 0;
    std::size_t solapes_b = 0;
    std::size_t escala_coincide_b = 0;
    std::size_t n_nulos_b_j12 = 0;
    std::size_t intentos_b_j12 = 0;
    std::size_t escala_coincide_b_j12 = 0;
    std::size_t escala_cruzada_j12 = 0;
    double p_det = 0.0;
    double surv_s3 = 0.0;
    double cob_s3 = 0.0;
    double p_det_j12 = 0.0;
};

ResumenBancoN4 consolidar_resultados(const std::vector<ResultadoSemillaN4>& lista) {
    ResumenBancoN4 out{};
    out.n_sem = lista.size();
    out.min_esc = lista.empty() ? 0u : 99u;
    std::vector<double> vals_s;
    vals_s.reserve(lista.size());
    for (const auto& r : lista) {
        if (r.calibracion_ok) ++out.n_cal_ok;
        out.min_esc = std::min(out.min_esc, r.escalas);
        if (r.evaluada) {
            ++out.n_sem_eval;
            vals_s.push_back(r.estadistico_s);
            if (r.rechazada) ++out.n_rech;
        }
        for (std::size_t j = 0; j < 3; ++j) out.k_med_j[j] += r.k_med_j[j];
        out.k_med += r.k_med;
        out.k_pdet_pos += r.k_pdet_pos;
        out.k_eval += r.k_eval;
        out.n_nulos_b += r.n_nulos_b;
        out.intentos_b += r.intentos_b;
        out.paradas_tempranas_b += r.paradas_tempranas_b;
        out.solapes_b += r.solapes_b;
        out.escala_coincide_b += r.escala_coincide_b;
        out.n_nulos_b_j12 += r.n_nulos_b_j12;
        out.intentos_b_j12 += r.intentos_b_j12;
        out.escala_coincide_b_j12 += r.escala_coincide_b_j12;
    }
    for (double v : vals_s) out.S_medio += v;
    if (!vals_s.empty()) out.S_medio /= static_cast<double>(vals_s.size());
    double sum_sq = 0.0;
    for (double v : vals_s) {
        const double d = v - out.S_medio;
        sum_sq += d * d;
    }
    const double var_S = (vals_s.size() > 1u) ? (sum_sq / static_cast<double>(vals_s.size() - 1u)) : 0.0;
    out.ee_S = (vals_s.size() > 1u) ? std::sqrt(var_S / static_cast<double>(vals_s.size())) : 0.0;
    out.z_S = (out.ee_S > 0.0) ? (out.S_medio / out.ee_S) : 0.0;

    out.p_det = out.intentos_b > 0u ? static_cast<double>(out.n_nulos_b) / static_cast<double>(out.intentos_b) : 0.0;
    out.surv_s3 = out.k_med > 0u ? static_cast<double>(out.k_pdet_pos) / static_cast<double>(out.k_med) : 0.0;
    out.cob_s3 = out.k_med > 0u ? static_cast<double>(out.k_eval) / static_cast<double>(out.k_med) : 0.0;
    out.p_det_j12 = out.intentos_b_j12 > 0u ? static_cast<double>(out.n_nulos_b_j12) / static_cast<double>(out.intentos_b_j12) : 0.0;
    out.escala_cruzada_j12 = (out.n_nulos_b_j12 >= out.escala_coincide_b_j12)
                                 ? (out.n_nulos_b_j12 - out.escala_coincide_b_j12)
                                 : 0u;
    return out;
}

// ============================================================================
// NV4.1: AR(1) phi=0.80, 24 semillas, E2E completar_nulos()
// ============================================================================
std::size_t ejecutar_nv4_1(std::size_t hilos_max) {
    constexpr std::size_t N_SEM = 24;
    constexpr std::size_t N_CAL = 12000;
    constexpr std::size_t N_OP = 12000;
    constexpr double PHI = 0.80;

    auto sem_cal = [](std::size_t s) noexcept {
        return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL4 ^ (0x9E3779B97F4A7C15ULL * (s + 1ULL)));
    };
    auto sem_op = [](std::size_t s) noexcept {
        return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL4 ^ (0xBF58476D1CE4E5B9ULL * (s + 101ULL)));
    };
    auto sem_nulo = [](std::size_t s) noexcept {
        return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL4 ^ (0x94D049BB133111EBULL * (s + 211ULL)));
    };

    const auto t0 = std::chrono::steady_clock::now();
    std::vector<ResultadoSemillaN4> res(N_SEM);
    std::atomic<std::size_t> sig{0};
    const std::size_t nw = std::clamp<std::size_t>(hilos_max, 1u, N_SEM);

    std::vector<std::thread> ths;
    ths.reserve(nw);
    for (std::size_t w = 0; w < nw; ++w) {
        ths.emplace_back([&]() {
            using ObsDelta = ObservadorSaltos<24, 64, 120>;
            using ObsLil = ObservadorLilyth<ObsDelta, 8, 64>;
            while (true) {
                const std::size_t s = sig.fetch_add(1u, std::memory_order_relaxed);
                if (s >= N_SEM) break;
                ResultadoSemillaN4 r{};
                r.semilla_idx = s;

                const auto x_cal = generar_ar1_estacionario(N_CAL, PHI, 1.0, sem_cal(s));
                const auto x_op = generar_ar1_estacionario(N_OP, PHI, 1.0, sem_op(s));

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
                r.escalas = od->escalas();

                ConfiguracionLilyth cl{};
                cl.alfa = ALFA_PROYECTO;
                cl.representacion = ModoRepresentacionLilyth::CRUDA;
                cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
                cl.modo_armonico = ModoArmonico::DESACTIVADO;
                cl.b_max_nulo = 8;
                cl.hilos_nulo = 1;
                cl.b_mc_resumen = 1999;
                cl.semilla_nulo = sem_nulo(s);
                cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
                cl.forzar_via_b = true;
                cl.usar_enmienda_nivel2 = true;
                cl.usar_enmienda_nivel3 = true;

                auto lil = std::make_unique<ObsLil>();
                r.calibracion_ok = lil->calibrar(*od, cl, cd, x_cal);
                if (r.calibracion_ok) {
                    od->reiniciar_flujo();
                    lil->reiniciar_flujo_conservando_purga();
                    for (std::size_t t = 0; t < N_OP; ++t) {
                        const auto sd = od->paso(x_op[t], NaN, true, 0, 0);
                        (void)lil->paso(sd, x_op[t], true, NaN, 0, 0);
                    }
                    lil->completar_nulos();
                    acumular_ranuras_cruda(*lil, r);
                }
                res[s] = r;
            }
        });
    }
    for (auto& t : ths) t.join();
    const double segs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    const ResumenBancoN4 b = consolidar_resultados(res);
    const bool ok_esc = (b.n_cal_ok == N_SEM && b.min_esc == 3u &&
                         b.k_med_j[0] > 0u && b.k_med_j[1] > 0u && b.k_med_j[2] > 0u &&
                         b.n_sem_eval >= 23u);
    const bool ok_cal = (b.n_sem_eval >= 23u && b.ee_S > 0.0 && std::isfinite(b.z_S) &&
                         std::abs(b.z_S) <= 2.0 && b.n_rech <= 3u);
    const bool ok_cov = (b.p_det >= 0.38 && b.surv_s3 >= 0.90 && b.cob_s3 >= 0.60 && b.solapes_b == 0u);
    const bool ok_j12 = (b.p_det_j12 >= 0.42 && b.escala_cruzada_j12 > 0u);

    std::printf("--------------------------------------------------------------------------------\n");
    std::printf("NV4.1 -- E2E AR(1) phi=0.80, 24 semillas, 3 escalas (tiempo=%.2f s, hilos=%zu)\n", segs, nw);
    std::printf("  cal=%zu/%zu min_esc=%zu | eval=%zu/%zu k_med=[%zu,%zu,%zu] (tot=%zu)\n",
                b.n_cal_ok, N_SEM, b.min_esc, b.n_sem_eval, N_SEM,
                b.k_med_j[0], b.k_med_j[1], b.k_med_j[2], b.k_med);
    std::printf("  S_medio=%+.6f ee_S=%.6f z_S=%+.6f | rech=%zu/%zu\n",
                b.S_medio, b.ee_S, b.z_S, b.n_rech, b.n_sem_eval);
    std::printf("  p_det=%.6f (%zu/%zu) | surv_s3=%.6f (%zu/%zu) | cob_s3=%.6f (%zu/%zu) | parB=%zu solapes=%zu\n",
                b.p_det, b.n_nulos_b, b.intentos_b, b.surv_s3, b.k_pdet_pos, b.k_med,
                b.cob_s3, b.k_eval, b.k_med, b.paradas_tempranas_b, b.solapes_b);
    std::printf("  j12: p_det_j12=%.6f (%zu/%zu) | escala_coincide_j12=%zu | escala_cruzada_j12=%zu\n",
                b.p_det_j12, b.n_nulos_b_j12, b.intentos_b_j12, b.escala_coincide_b_j12, b.escala_cruzada_j12);
    std::printf("  [%s] ok_esc : cal==24/24 && min_esc==3 && k_med_j[0..2]>0 && n_sem_eval>=23\n", ok_esc ? "PASA" : "FALLA");
    std::printf("  [%s] ok_cal : |z_S|<=2.0 && rech<=3/24\n", ok_cal ? "PASA" : "FALLA");
    std::printf("  [%s] ok_cov : p_det>=0.38 && surv_s3>=0.90 && cob_s3>=0.60 && solapes==0\n", ok_cov ? "PASA" : "FALLA");
    std::printf("  [%s] ok_j12 : p_det_j12>=0.42 && escala_cruzada_j12>0\n", ok_j12 ? "PASA" : "FALLA");

    const bool ok = ok_esc && ok_cal && ok_cov && ok_j12;
    std::printf("FALLAS NV4.1: %zu\n", ok ? std::size_t{0} : std::size_t{1});
    std::fflush(stdout);
    return ok ? 0u : 1u;
}

// ============================================================================
// NV4.2: ARFIMA(1, d=+0.35, 0) phi1=0.40, mu=50.0, 16 semillas, E2E completar_nulos()
// ============================================================================
std::size_t ejecutar_nv4_2(std::size_t hilos_max) {
    constexpr std::size_t N_SEM = 16;
    constexpr std::size_t N_CAL = 12000;
    constexpr std::size_t N_OP = 8000;
    constexpr double PHI1 = 0.40;
    constexpr double D_ARF = 0.35;
    constexpr double MU_SERIE = 50.0;

    auto sem_cal = [](std::size_t s) noexcept {
        return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL4 ^ (0x517CC1B727220A95ULL * (s + 1ULL)));
    };
    auto sem_op = [](std::size_t s) noexcept {
        return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL4 ^ (0xD6E8FEB86659FD93ULL * (s + 307ULL)));
    };
    auto sem_nulo = [](std::size_t s) noexcept {
        return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL4 ^ (0xA0761D6478BD642FULL * (s + 409ULL)));
    };

    const auto t0 = std::chrono::steady_clock::now();
    std::vector<ResultadoSemillaN4> res(N_SEM);
    std::atomic<std::size_t> sig{0};
    const std::size_t nw = std::clamp<std::size_t>(hilos_max, 1u, N_SEM);

    std::vector<std::thread> ths;
    ths.reserve(nw);
    for (std::size_t w = 0; w < nw; ++w) {
        ths.emplace_back([&]() {
            using ObsDelta = ObservadorSaltos<24, 64, 120>;
            using ObsLil128 = ObservadorLilyth<ObsDelta, 8, 128>;
            while (true) {
                const std::size_t s = sig.fetch_add(1u, std::memory_order_relaxed);
                if (s >= N_SEM) break;
                ResultadoSemillaN4 r{};
                r.semilla_idx = s;

                const auto x_cal = generar_arfima1d0(N_CAL, PHI1, D_ARF, MU_SERIE, 1.0, sem_cal(s));
                const auto x_op = generar_arfima1d0(N_OP, PHI1, D_ARF, MU_SERIE, 1.0, sem_op(s));

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
                r.escalas = od->escalas();

                ConfiguracionLilyth cl{};
                cl.alfa = ALFA_PROYECTO;
                cl.representacion = ModoRepresentacionLilyth::CRUDA;
                cl.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
                cl.d_whittle_cruda = D_ARF;
                cl.modo_armonico = ModoArmonico::DESACTIVADO;
                cl.b_max_nulo = 4;
                cl.hilos_nulo = 1;
                cl.b_mc_resumen = 1999;
                cl.semilla_nulo = sem_nulo(s);
                cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
                cl.forzar_via_b = true;
                cl.usar_enmienda_nivel2 = true;
                cl.usar_enmienda_nivel3 = true;

                auto lil = std::make_unique<ObsLil128>();
                r.calibracion_ok = lil->calibrar(*od, cl, cd, x_cal);
                if (r.calibracion_ok) {
                    od->reiniciar_flujo();
                    lil->reiniciar_flujo_conservando_purga();
                    for (std::size_t t = 0; t < N_OP; ++t) {
                        const auto sd = od->paso(x_op[t], NaN, true, 0, 0);
                        (void)lil->paso(sd, x_op[t], true, NaN, 0, 0);
                    }
                    lil->completar_nulos();
                    acumular_ranuras_cruda(*lil, r);
                }
                res[s] = r;
            }
        });
    }
    for (auto& t : ths) t.join();
    const double segs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    const ResumenBancoN4 b = consolidar_resultados(res);
    const bool ok_esc = (b.n_cal_ok == N_SEM && b.min_esc == 3u &&
                         b.k_med_j[0] > 0u && b.k_med_j[1] > 0u && b.k_med_j[2] > 0u &&
                         b.n_sem_eval >= 14u);
    const bool ok_cal = (b.n_sem_eval >= 14u && b.ee_S > 0.0 && std::isfinite(b.z_S) &&
                         std::abs(b.z_S) <= 2.0 && b.n_rech <= 2u);
    const bool ok_cov = (b.p_det >= 0.16 && b.cob_s3 >= 0.25 && b.solapes_b == 0u);
    const bool ok_j12 = (b.p_det_j12 >= 0.15 && b.escala_cruzada_j12 > 0u);

    std::printf("--------------------------------------------------------------------------------\n");
    std::printf("NV4.2 -- E2E ARFIMA(1, d=+0.35, 0), 16 semillas, 3 escalas (tiempo=%.2f s, hilos=%zu)\n", segs, nw);
    std::printf("  cal=%zu/%zu min_esc=%zu | eval=%zu/%zu k_med=[%zu,%zu,%zu] (tot=%zu)\n",
                b.n_cal_ok, N_SEM, b.min_esc, b.n_sem_eval, N_SEM,
                b.k_med_j[0], b.k_med_j[1], b.k_med_j[2], b.k_med);
    std::printf("  S_medio=%+.6f ee_S=%.6f z_S=%+.6f | rech=%zu/%zu\n",
                b.S_medio, b.ee_S, b.z_S, b.n_rech, b.n_sem_eval);
    std::printf("  p_det=%.6f (%zu/%zu) | surv_s3=%.6f (%zu/%zu) | cob_s3=%.6f (%zu/%zu) | parB=%zu solapes=%zu\n",
                b.p_det, b.n_nulos_b, b.intentos_b, b.surv_s3, b.k_pdet_pos, b.k_med,
                b.cob_s3, b.k_eval, b.k_med, b.paradas_tempranas_b, b.solapes_b);
    std::printf("  j12: p_det_j12=%.6f (%zu/%zu) | escala_coincide_j12=%zu | escala_cruzada_j12=%zu\n",
                b.p_det_j12, b.n_nulos_b_j12, b.intentos_b_j12, b.escala_coincide_b_j12, b.escala_cruzada_j12);
    std::printf("  [%s] ok_esc : cal==16/16 && min_esc==3 && k_med_j[0..2]>0 && n_sem_eval>=14\n", ok_esc ? "PASA" : "FALLA");
    std::printf("  [%s] ok_cal : |z_S|<=2.0 && rech<=2/16\n", ok_cal ? "PASA" : "FALLA");
    std::printf("  [%s] ok_cov : p_det>=0.16 && cob_s3>=0.25 && solapes==0\n", ok_cov ? "PASA" : "FALLA");
    std::printf("  [%s] ok_j12 : p_det_j12>=0.15 && escala_cruzada_j12>0\n", ok_j12 ? "PASA" : "FALLA");

    const bool ok = ok_esc && ok_cal && ok_cov && ok_j12;
    std::printf("FALLAS NV4.2: %zu\n", ok ? std::size_t{0} : std::size_t{1});
    std::fflush(stdout);
    return ok ? 0u : 1u;
}

// ============================================================================
// NV4.3: Contraste pareado ciego P0..P3 + Regla 6
// ============================================================================
struct AcumPolitica {
    std::size_t adm_global = 0;
    std::size_t adm_j2 = 0;
    double suma_gdif = 0.0;
};

struct ResEventoNV43 {
    std::size_t intentos = 0;
    std::size_t intentos_j2 = 0;
    AcumPolitica p0{};
    AcumPolitica p1a{};
    AcumPolitica p2{};
    AcumPolitica p3{};
};

std::size_t ejecutar_nv4_3(std::size_t hilos_max) {
    using ObsDelta = ObservadorSaltos<24, 64, 120>;
    using ObsLil = ObservadorLilyth<ObsDelta, 8, 64>;

    constexpr std::uint64_t SEM_43 = SEMILLA_CIEGA_NIVEL4 + 9001ULL;
    constexpr std::size_t N_CAL = 12000;
    constexpr std::size_t N_OP = 120000;
    constexpr std::size_t MAX_EV_ESCALA = 16;
    constexpr std::size_t N_MC = 150;
    constexpr double PHI = 0.80;

    const auto t0 = std::chrono::steady_clock::now();
    const auto x_cal = generar_ar1_estacionario(N_CAL, PHI, 1.0, detalle_nulo::splitmix64(SEM_43 + 1ULL));
    const auto x_op = generar_ar1_estacionario(N_OP, PHI, 1.0, detalle_nulo::splitmix64(SEM_43 + 2ULL));

    ConfiguracionObservador cd{};
    cd.integracion = Integracion::FLUJO;
    cd.agregacion = AgregacionEscala::SUMA;
    cd.arl0_pasos_nativos = 450.0;
    cd.alfa = ALFA_PROYECTO;
    cd.purga = false;
    cd.bins_fase = 1;
    cd.medir_tau50 = false;
    cd.escalas_max = 3;
    cd.semilla = detalle_nulo::splitmix64(SEM_43 + 3ULL);

    auto od = std::make_unique<ObsDelta>();
    od->calibrar(cd, x_cal);

    ConfiguracionLilyth cl{};
    cl.alfa = ALFA_PROYECTO;
    cl.representacion = ModoRepresentacionLilyth::CRUDA;
    cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
    cl.modo_armonico = ModoArmonico::DESACTIVADO;
    cl.b_max_nulo = 4;
    cl.hilos_nulo = 1;
    cl.b_mc_resumen = 1999;
    cl.semilla_nulo = detalle_nulo::splitmix64(SEM_43 + 4ULL);
    cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cl.forzar_via_b = true;
    cl.usar_enmienda_nivel2 = true;
    cl.usar_enmienda_nivel3 = true;

    auto lil = std::make_unique<ObsLil>();
    const bool cal_ok = lil->calibrar(*od, cl, cd, x_cal);
    const std::size_t escalas_det = od->escalas();
    const std::size_t t_max_ef = lil->configuracion().t_max_efectivo;
    const auto& mod_lin = lil->ficha_nulo_cruda().modelo_lineal;

    od->reiniciar_flujo();
    lil->reiniciar_flujo_conservando_purga();
    for (std::size_t t = 0; t < N_OP; ++t) {
        const auto sd = od->paso(x_op[t], NaN, true, 0, 0);
        (void)lil->paso(sd, x_op[t], true, NaN, 0, 0);
    }

    struct EventoSel {
        ObsLil::Ranura ranura{};
        std::vector<double> hist_prev{};
        std::vector<double> perfil{};
    };
    std::vector<EventoSel> seleccionados;
    std::array<std::size_t, 3> cuenta_j{};
    for (const auto& r : lil->cerradas_historico()) {
        const auto d = r.cruda.desenlace;
        if (d != DesenlaceRanura::RESTITUIDA && d != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
        if (r.perfil_corrupto_z2_cruda || r.escala_origen < 0 || r.escala_origen >= 3) continue;
        if (!std::isfinite(r.g_disparo) || r.n_hist_prev_cruda < mod_lin.p) continue;
        const auto pf = lil->perfil_evento_n2(r, false);
        if (pf.empty()) continue;
        const std::size_t j = static_cast<std::size_t>(r.escala_origen);
        if (cuenta_j[j] >= MAX_EV_ESCALA) continue;
        ++cuenta_j[j];
        EventoSel ev{};
        ev.ranura = r;
        ev.hist_prev.assign(r.hist_prev_cruda_est.begin(),
                            r.hist_prev_cruda_est.begin() + static_cast<std::ptrdiff_t>(r.n_hist_prev_cruda));
        ev.perfil.assign(pf.begin(), pf.end());
        seleccionados.push_back(std::move(ev));
    }

    const std::size_t n_ev = seleccionados.size();
    std::vector<ResEventoNV43> res_ev(n_ev);
    std::atomic<std::size_t> sig{0};
    const std::size_t nw = n_ev > 0 ? std::clamp<std::size_t>(hilos_max, 1u, n_ev) : 1u;

    std::vector<std::thread> ths;
    ths.reserve(nw);
    for (std::size_t w = 0; w < nw; ++w) {
        ths.emplace_back([&]() {
            auto lil_w = std::make_unique<ObsLil>(*lil);
            while (true) {
                const std::size_t k = sig.fetch_add(1u, std::memory_order_relaxed);
                if (k >= n_ev) break;
                const auto& ev = seleccionados[k];
                const auto& r = ev.ranura;
                const auto dims = dimensiones_simulacion_via_b(
                    r.w_ancla, r.guarda_ancla, ev.perfil.size(), t_max_ef, Integracion::FLUJO);
                ResEventoNV43 out_k{};
                out_k.intentos = N_MC;
                if (r.escala_origen == 2) out_k.intentos_j2 = N_MC;

                for (std::size_t mc = 0; mc < N_MC; ++mc) {
                    const std::uint64_t sem_mc = detalle_nulo::splitmix64(
                        SEM_43 ^ (0x9E3779B97F4A7C15ULL * (k + 1ULL)) ^ (0xBF58476D1CE4E5B9ULL * (mc + 1ULL)));
                    std::mt19937_64 rng_t(sem_mc);
                    const std::int64_t t_iny = dims.t_iny_min + static_cast<std::int64_t>(rng_t() % dims.w_sorteo);
                    const std::array<std::int64_t, 1> ts{t_iny};

                    const auto r_mc = lil_w->evaluar_inyecciones_realizacion(
                        false, ev.hist_prev, ev.perfil, r.signo, r.w_ancla, r.guarda_ancla,
                        ts, dims.n_sim, detalle_nulo::splitmix64(sem_mc + 0x517CC1B727220A95ULL),
                        false, true, r.escala_origen);

                    if (r_mc.eventos_admitidos.size() == 1u) {
                        const auto& e_n = r_mc.eventos_admitidos[0];
                        if (!std::isfinite(e_n.g_disparo)) continue;
                        const double gdif = e_n.g_disparo - r.g_disparo;
                        const bool mismo_j = (e_n.escala_origen == r.escala_origen);
                        const bool g_excede = (e_n.g_disparo >= r.g_disparo);

                        auto anotar = [&](AcumPolitica& ac, bool adm) {
                            if (!adm) return;
                            ++ac.adm_global;
                            if (r.escala_origen == 2) ++ac.adm_j2;
                            ac.suma_gdif += gdif;
                        };
                        anotar(out_k.p0,  true);
                        anotar(out_k.p1a, mismo_j);
                        anotar(out_k.p2,  g_excede);
                        anotar(out_k.p3,  mismo_j && g_excede);
                    }
                }
                res_ev[k] = out_k;
            }
        });
    }
    for (auto& t : ths) t.join();
    const double segs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::size_t tot_int = 0, tot_int_j2 = 0;
    AcumPolitica tot_p0{}, tot_p1a{}, tot_p2{}, tot_p3{};
    auto sumar = [](AcumPolitica& dst, const AcumPolitica& src) noexcept {
        dst.adm_global += src.adm_global;
        dst.adm_j2 += src.adm_j2;
        dst.suma_gdif += src.suma_gdif;
    };
    for (const auto& rk : res_ev) {
        tot_int += rk.intentos;
        tot_int_j2 += rk.intentos_j2;
        sumar(tot_p0, rk.p0);
        sumar(tot_p1a, rk.p1a);
        sumar(tot_p2, rk.p2);
        sumar(tot_p3, rk.p3);
    }

    auto calc_pdet = [&](const AcumPolitica& ac) {
        return tot_int > 0 ? static_cast<double>(ac.adm_global) / static_cast<double>(tot_int) : 0.0;
    };
    auto calc_pdet_j2 = [&](const AcumPolitica& ac) {
        return tot_int_j2 > 0 ? static_cast<double>(ac.adm_j2) / static_cast<double>(tot_int_j2) : 0.0;
    };
    auto calc_egdif = [&](const AcumPolitica& ac) {
        return ac.adm_global > 0 ? (ac.suma_gdif / static_cast<double>(ac.adm_global)) : NaN;
    };

    const double pdet_p0 = calc_pdet(tot_p0), pdet_p1a = calc_pdet(tot_p1a), pdet_p2 = calc_pdet(tot_p2), pdet_p3 = calc_pdet(tot_p3);
    const double j2_p0 = calc_pdet_j2(tot_p0), j2_p1a = calc_pdet_j2(tot_p1a), j2_p2 = calc_pdet_j2(tot_p2), j2_p3 = calc_pdet_j2(tot_p3);
    const double eg_p0 = calc_egdif(tot_p0), eg_p1a = calc_egdif(tot_p1a), eg_p2 = calc_egdif(tot_p2), eg_p3 = calc_egdif(tot_p3);

    const ConfiguracionLilyth cfg_def{};
    const bool ok_prep = (cal_ok && escalas_det == 3u && cuenta_j[0] > 0u && cuenta_j[1] > 0u && cuenta_j[2] > 0u && tot_int > 0u && tot_int_j2 > 0u);
    const bool ok_p2 = ok_prep && (pdet_p0 > 0.0) && (pdet_p2 <= 0.70 * pdet_p0) &&
                       std::isfinite(eg_p0) && std::isfinite(eg_p2) && (eg_p2 > eg_p0 + 0.05);
    const bool ok_p1 = ok_prep && (j2_p0 > 0.0) && (j2_p1a <= 0.75 * j2_p0) && (j2_p3 <= 0.50 * j2_p0);
    const bool ok_r6 = (!cfg_def.usar_enmienda_nivel2 && !cfg_def.usar_enmienda_nivel3);

    std::printf("--------------------------------------------------------------------------------\n");
    std::printf("NV4.3 -- Contraste pareado P0..P3 + Regla 6 (tiempo=%.2f s, hilos=%zu)\n", segs, nw);
    std::printf("  cal=%d J=%zu eventos=[%zu,%zu,%zu] (tot=%zu, N_MC=%zu -> int=%zu, int_j2=%zu)\n",
                cal_ok ? 1 : 0, escalas_det, cuenta_j[0], cuenta_j[1], cuenta_j[2], n_ev, N_MC, tot_int, tot_int_j2);
    std::printf("  P0  : p_det=%.6f (%zu/%zu) | p_det_j2=%.6f (%zu/%zu) | Eg_dif=%+.6f\n",
                pdet_p0, tot_p0.adm_global, tot_int, j2_p0, tot_p0.adm_j2, tot_int_j2, eg_p0);
    std::printf("  P1a : p_det=%.6f (%zu/%zu) | p_det_j2=%.6f (%zu/%zu) | Eg_dif=%+.6f\n",
                pdet_p1a, tot_p1a.adm_global, tot_int, j2_p1a, tot_p1a.adm_j2, tot_int_j2, eg_p1a);
    std::printf("  P2  : p_det=%.6f (%zu/%zu) | p_det_j2=%.6f (%zu/%zu) | Eg_dif=%+.6f\n",
                pdet_p2, tot_p2.adm_global, tot_int, j2_p2, tot_p2.adm_j2, tot_int_j2, eg_p2);
    std::printf("  P3  : p_det=%.6f (%zu/%zu) | p_det_j2=%.6f (%zu/%zu) | Eg_dif=%+.6f\n",
                pdet_p3, tot_p3.adm_global, tot_int, j2_p3, tot_p3.adm_j2, tot_int_j2, eg_p3);
    std::printf("  Regla 6 defaults: usar_enmienda_nivel2=%d, usar_enmienda_nivel3=%d\n",
                cfg_def.usar_enmienda_nivel2 ? 1 : 0, cfg_def.usar_enmienda_nivel3 ? 1 : 0);
    std::printf("  [%s] ok_p2 : p_det(P2)<=0.70*p_det(P0) && Eg_dif(P2)>Eg_dif(P0)+0.05\n", ok_p2 ? "PASA" : "FALLA");
    std::printf("  [%s] ok_p1 : p_det_j2(P1a)<=0.75*p_det_j2(P0) && p_det_j2(P3)<=0.50*p_det_j2(P0)\n", ok_p1 ? "PASA" : "FALLA");
    std::printf("  [%s] ok_r6 : !usar_enmienda_nivel2 && !usar_enmienda_nivel3 por omision\n", ok_r6 ? "PASA" : "FALLA");

    const bool ok = ok_p2 && ok_p1 && ok_r6;
    std::printf("FALLAS NV4.3: %zu\n", ok ? std::size_t{0} : std::size_t{1});
    std::fflush(stdout);
    return ok ? 0u : 1u;
}

}  // namespace

int main(int argc, char** argv) {
    const OpcionesCli opt = parsear_cli(argc, argv);

    std::printf("================================================================================\n");
    std::printf(" BATERIA SELLADA NIVEL 4 -- NO-OP DEMOSTRADO (P0_NIVEL3_BASE)\n");
    std::printf("  SEMILLA_CIEGA_NIVEL4 = %lluULL | hilos = %zu | modo = %s\n",
                static_cast<unsigned long long>(SEMILLA_CIEGA_NIVEL4),
                opt.hilos,
                opt.solo_operador ? "solo-operador (NV4.1 + NV4.2)" : "completo (NV4.1..NV4.3)");
    std::printf("================================================================================\n");

    std::size_t fallas = 0;
    fallas += ejecutar_nv4_1(opt.hilos);
    fallas += ejecutar_nv4_2(opt.hilos);
    if (!opt.solo_operador) {
        fallas += ejecutar_nv4_3(opt.hilos);
    } else {
        std::printf("--------------------------------------------------------------------------------\n");
        std::printf("NV4.3 -- OMITIDO (--solo-operador)\n");
    }

    std::printf("================================================================================\n");
    std::printf("[%s] BATERIA NIVEL 4 COMPLETA\n", (fallas == 0u) ? "PASA" : "FALLA");
    std::printf("FALLAS NV4: %zu\n", fallas);
    std::fflush(stdout);
    return (fallas == 0u) ? 0 : 1;
}
