// ============================================================================
// test_nv3_2_banco_r2.cpp -- Control sellado NV3.2 (T3.6):
//   Banco R2: AR(1) phi=0.88 + Hawkes 2 signos bajo H0 exacto,
//   24 semillas ciegas (R2a s=0..11, R2b s=12..23), forzar_via_b = true,
//   usar_enmienda_nivel2 = true, usar_enmienda_nivel3 = true.
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
#include <memory>
#include <random>
#include <span>
#include <string_view>
#include <thread>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;

using ObsH = ObservadorSaltos<24, 4, 120>;
constexpr std::size_t TMAXC = 128;
using LilH = ObservadorLilyth<ObsH, 8, TMAXC, 2048>;

constexpr std::uint64_t SEMILLA_CIEGA_NIVEL3 = 20261124ULL;
constexpr std::size_t N_SEM = 24;
constexpr std::size_t N_CAL = 8192;
constexpr std::size_t N_TEST = 6000;
constexpr double PHI_R2 = 0.88;
constexpr double TASA_MADRE_R2 = 1.0 / 300.0;
constexpr std::size_t D0_R2 = 6;
constexpr double TAU_RET_R2 = 4.0;

std::uint64_t sem_cal(std::size_t s) noexcept {
    return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL3 ^ (0x517CC1B727220A95ULL * (s + 1ULL)));
}

std::uint64_t sem_op(std::size_t s) noexcept {
    return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL3 ^ (0xD6E8FEB86659FD93ULL * (s + 307ULL)));
}

std::uint64_t sem_nulo(std::size_t s) noexcept {
    return detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL3 ^ (0xA0761D6478BD642FULL * (s + 409ULL)));
}

std::vector<double> generar_ar1_limpio_r2(std::size_t n, double phi, double sd_e, std::uint64_t semilla) {
    std::vector<double> x(n, 0.0);
    if (n == 0) return x;
    std::mt19937_64 rng(semilla);
    std::normal_distribution<double> n01(0.0, 1.0);
    std::normal_distribution<double> nd_e(0.0, sd_e);
    double prev = n01(rng);
    x[0] = prev;
    for (std::size_t t = 1; t < n; ++t) {
        prev = phi * prev + nd_e(rng);
        x[t] = prev;
    }
    return x;
}

struct GanchoHijasHawkes {
    double tasa_madre = TASA_MADRE_R2;
    double n_mismo = 0.0;
    double n_contrario = 0.0;
    std::size_t d0 = D0_R2;
    double tau_ret = TAU_RET_R2;
    double amp = 0.0;

    void preparar(std::mt19937_64& rng, std::vector<double>& eps_nulo, std::size_t burn, std::size_t n_sim) const {
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        for (std::size_t t = 0; t < n_sim; ++t) {
            if (u01(rng) < tasa_madre) {
                const double sgn = (u01(rng) < 0.5) ? +1.0 : -1.0;
                eps_nulo[burn + t] += sgn * amp;
            }
        }
    }

    template <class SalidaLil, class Lil>
    void en_paso(std::size_t t, const SalidaLil& o, const Lil& lil_local, std::mt19937_64& rng,
                 std::vector<double>& eps, std::vector<double>& eps_nulo, std::size_t burn, std::size_t n_sim) const {
        if (!o.abrio_ranura || o.indice_ranura_abierta >= lil_local.ranuras().size()) return;
        const auto& ran = lil_local.ranuras()[o.indice_ranura_abierta];
        if (!ran.ocupada || ran.t_apertura != o.t) return;

        std::geometric_distribution<int> geo(1.0 - std::exp(-1.0 / tau_ret));
        std::poisson_distribution<int> pm(std::max(n_mismo, 1e-300));
        std::poisson_distribution<int> pc(std::max(n_contrario, 1e-300));
        const int km = (n_mismo > 0.0) ? pm(rng) : 0;
        const int kc = (n_contrario > 0.0) ? pc(rng) : 0;
        for (int k = 0; k < km + kc; ++k) {
            const std::size_t th = t + d0 + static_cast<std::size_t>(geo(rng));
            if (th < n_sim) {
                const int sig_h = (k < km) ? ran.signo : -ran.signo;
                const double delta_e = static_cast<double>(sig_h) * amp;
                eps[burn + th] += delta_e;
                eps_nulo[burn + th] += delta_e;
            }
        }
    }
};

struct ResultadoSemilla {
    std::size_t semilla_idx = 0;
    bool es_r2b = false;
    bool calibracion_ok = false;
    bool evaluada = false;
    double estadistico_s = 0.0;
    double p_bilateral = 1.0;
    bool rechazada = false;
    std::size_t n_eventos_contraste = 0;
    std::size_t m_min = 0;
    std::size_t t_max_ef = 0;
    std::size_t cortadas_real = 0;
    std::size_t k_med = 0;
    std::size_t k_pdet_pos = 0;
    std::size_t k_eval = 0;
    std::size_t n_nulos_b = 0;
    std::size_t intentos_b = 0;
    std::size_t cortadas_b = 0;
    std::size_t paradas_tempranas_b = 0;
    std::size_t solapes_b = 0;
    std::size_t con_ajena_previa_b = 0;
    std::size_t con_reapertura_posterior_b = 0;
    std::size_t escala_coincide_b = 0;
};

ResultadoSemilla ejecutar_semilla(std::size_t s) {
    ResultadoSemilla res{};
    res.semilla_idx = s;
    res.es_r2b = (s >= 12u);

    const double n_mismo = res.es_r2b ? 0.45 : 0.30;
    const double n_contrario = res.es_r2b ? 0.25 : 0.15;
    const double sd_e = std::sqrt(1.0 - PHI_R2 * PHI_R2);
    const double amp = 9.5 * sd_e;

    // 1. Calibracion limpia
    const auto x_cal = generar_ar1_limpio_r2(N_CAL, PHI_R2, sd_e, sem_cal(s));

    ConfiguracionObservador cd{};
    cd.integracion = Integracion::FLUJO;
    cd.agregacion = AgregacionEscala::SUMA;
    cd.arl0_pasos_nativos = 800.0;
    cd.escalas_max = 1;
    cd.alfa = ALFA_PROYECTO;
    cd.purga = false;
    cd.bins_fase = 1;
    cd.medir_tau50 = false;
    cd.semilla = sem_cal(s);

    auto od = std::make_unique<ObsH>();
    od->calibrar(cd, x_cal);

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

    auto lil = std::make_unique<LilH>();
    const bool cal_ok = lil->calibrar(*od, cl, cd, x_cal);
    res.calibracion_ok = cal_ok;
    if (!cal_ok) return res;

    auto od_cal = std::make_unique<ObsH>(*od);
    auto lil_plantilla = std::make_unique<LilH>(*lil);

    // 2. Tramo operativo AR(1) + Hawkes 2 signos bajo H0 exacto
    od->reiniciar_flujo();
    lil->reiniciar_flujo_conservando_purga();

    {
        std::mt19937_64 rng_op(sem_op(s));
        std::normal_distribution<double> nd_e(0.0, sd_e);
        std::normal_distribution<double> n01(0.0, 1.0);
        std::uniform_real_distribution<double> u01(0.0, 1.0);
        std::geometric_distribution<int> geo(1.0 - std::exp(-1.0 / TAU_RET_R2));
        std::poisson_distribution<int> pm(std::max(n_mismo, 1e-300));
        std::poisson_distribution<int> pc(std::max(n_contrario, 1e-300));

        std::vector<double> eps_op(N_TEST, 0.0);
        for (std::size_t t = 0; t < N_TEST; ++t) {
            eps_op[t] = nd_e(rng_op);
            if (u01(rng_op) < TASA_MADRE_R2) {
                const double sgn = (u01(rng_op) < 0.5) ? +1.0 : -1.0;
                eps_op[t] += sgn * amp;
            }
        }

        double prev = n01(rng_op);
        for (std::size_t t = 0; t < N_TEST; ++t) {
            prev = PHI_R2 * prev + eps_op[t];
            const double xt = prev;
            const auto sd = od->paso(xt, NaN, true, 0, 0);
            const auto o = lil->paso(sd, xt, true, NaN, 0, 0);
            if (o.abrio_ranura && o.indice_ranura_abierta < lil->ranuras().size()) {
                const auto& ran = lil->ranuras()[o.indice_ranura_abierta];
                if (ran.ocupada && ran.t_apertura == o.t) {
                    const int km = (n_mismo > 0.0) ? pm(rng_op) : 0;
                    const int kc = (n_contrario > 0.0) ? pc(rng_op) : 0;
                    for (int k = 0; k < km + kc; ++k) {
                        const std::size_t th = t + D0_R2 + static_cast<std::size_t>(geo(rng_op));
                        if (th < N_TEST) {
                            const int sig_h = (k < km) ? ran.signo : -ran.signo;
                            eps_op[th] += static_cast<double>(sig_h) * amp;
                        }
                    }
                }
            }
        }
    }

    // 3. Evaluacion Via (b) con GanchoHijasHawkes + contraste INTERCAMBIABLE
    const auto& fn = lil->ficha_nulo_cruda();
    const std::size_t t_max_ef = lil->configuracion().t_max_efectivo;
    res.m_min = fn.m_min;
    res.t_max_ef = t_max_ef;

    const GanchoHijasHawkes gancho{TASA_MADRE_R2, n_mismo, n_contrario, D0_R2, TAU_RET_R2, amp};
    const auto& cerradas = lil->cerradas_historico();
    std::vector<std::vector<double>> estratos;
    std::size_t solapes_seed = 0;

    for (std::size_t idx = 0; idx < cerradas.size(); ++idx) {
        const auto& r = cerradas[idx];
        const auto d = r.cruda.desenlace;
        if (d == DesenlaceRanura::CORTADA_POR_REPLICA) {
            ++res.cortadas_real;
            continue;
        }
        if (d != DesenlaceRanura::RESTITUIDA && d != DesenlaceRanura::CENSURADA_HORIZONTE) {
            continue;
        }
        ++res.k_med;

        if (r.perfil_corrupto_z2_cruda) continue;
        if (!fn.modelo_lineal.valido || r.n_hist_prev_disp_cruda < fn.modelo_lineal.p) continue;
        const std::span<const double> pf = lil->perfil_evento_n2(r, false);
        const std::span<const double> einy_n2 = lil->e_iny_evento_n2(r, false);
        if (pf.empty()) continue;
        if (einy_n2.size() != pf.size() && r.n_hist_prev_disp_cruda > MAX_PERFIL_EVENTO_LILYTH) continue;

        const std::span<const double> hp(r.hist_prev_cruda_est.data(), r.n_hist_prev_cruda);
        ReferenciaNuloEvento ref{};
        std::vector<EventoNuloLilyth<TMAXC>> ev_b;
        detalle_nulo::ejecutar_via_b_evento<TMAXC, RasgosDelta<ObsH>::b_max>(
            *od_cal, *lil_plantilla, Integracion::FLUJO, false, false,
            cl.nulo_sin_memoria, fn.segmentos_calma, fn.modelo_lineal,
            lil_plantilla->contexto_purga_nulo(), NaN, hp, pf,
            r.signo, r.w_ancla, r.guarda_ancla, fn.m_min, t_max_ef,
            cl.semilla_nulo ^ (0x9E3779B97F4A7C15ULL * (idx + 1ULL)),
            ref, ev_b, solapes_seed, einy_n2, true, r.escala_origen, gancho);

        if (std::isfinite(ref.p_det_inicial) && ref.p_det_inicial > 0.0) ++res.k_pdet_pos;
        res.n_nulos_b += ref.n_nulos_b;
        res.intentos_b += ref.intentos_b;
        res.cortadas_b += ref.cortadas_b;
        res.paradas_tempranas_b += ref.paradas_tempranas_b;
        res.solapes_b += ref.solapes_b;
        res.con_ajena_previa_b += ref.con_ajena_previa_b;
        res.con_reapertura_posterior_b += ref.con_reapertura_posterior_b;
        res.escala_coincide_b += ref.escala_coincide_b;

        if (ev_b.size() >= fn.m_min) {
            ++res.k_eval;
            std::vector<double> estrato;
            estrato.reserve(1u + ev_b.size());
            estrato.push_back(detalle_nulo::codificar_tau(r.cruda.tau_rest_consolidado, t_max_ef));
            for (const auto& ev : ev_b) {
                estrato.push_back(detalle_nulo::codificar_tau(ev.tau_rest_consolidado, t_max_ef));
            }
            estratos.push_back(std::move(estrato));
        }
    }

    const auto ctr = detalle_nulo::evaluar_contraste_intercambiable(
        estratos, cl.b_mc_resumen, cl.semilla_nulo ^ 0x9876543210FEDCBAULL, cl.alfa / 4.0);
    res.n_eventos_contraste = ctr.n_eventos;
    if (ctr.valido && std::isfinite(ctr.estadistico_s) && std::isfinite(ctr.p_bilateral)) {
        res.evaluada = true;
        res.estadistico_s = ctr.estadistico_s;
        res.p_bilateral = ctr.p_bilateral;
        res.rechazada = (res.p_bilateral < 0.05);
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

struct AcumGrupo {
    std::size_t n_sem = 0;
    std::size_t n_cal_ok = 0;
    std::size_t n_eval = 0;
    std::size_t n_s_neg = 0;
    std::size_t n_rech = 0;
    std::size_t ev_ctr = 0;
    std::size_t cort_real = 0;
    std::size_t k_med = 0;
    std::size_t k_pdet_pos = 0;
    std::size_t k_eval = 0;
    std::size_t n_nulos_b = 0;
    std::size_t intentos_b = 0;
    std::size_t cortadas_b = 0;
    std::size_t paradas_b = 0;
    std::size_t solapes_b = 0;
    std::size_t ajena_prev_b = 0;
    std::size_t reap_post_b = 0;
    std::size_t escala_coinc_b = 0;
    std::vector<double> s_vals;

    void agregar(const ResultadoSemilla& r) {
        ++n_sem;
        if (r.calibracion_ok) ++n_cal_ok;
        if (r.evaluada) {
            ++n_eval;
            s_vals.push_back(r.estadistico_s);
            if (r.estadistico_s < 0.0) ++n_s_neg;
            if (r.rechazada) ++n_rech;
        }
        ev_ctr += r.n_eventos_contraste;
        cort_real += r.cortadas_real;
        k_med += r.k_med;
        k_pdet_pos += r.k_pdet_pos;
        k_eval += r.k_eval;
        n_nulos_b += r.n_nulos_b;
        intentos_b += r.intentos_b;
        cortadas_b += r.cortadas_b;
        paradas_b += r.paradas_tempranas_b;
        solapes_b += r.solapes_b;
        ajena_prev_b += r.con_ajena_previa_b;
        reap_post_b += r.con_reapertura_posterior_b;
        escala_coinc_b += r.escala_coincide_b;
    }

    double s_medio() const noexcept {
        if (s_vals.empty()) return 0.0;
        double sum = 0.0;
        for (double v : s_vals) sum += v;
        return sum / static_cast<double>(s_vals.size());
    }

    double ee_s() const noexcept {
        if (s_vals.size() < 2) return 0.0;
        const double m = s_medio();
        double ss = 0.0;
        for (double v : s_vals) {
            const double d = v - m;
            ss += d * d;
        }
        return std::sqrt(ss / static_cast<double>(s_vals.size() * (s_vals.size() - 1u)));
    }
};

}  // namespace

int main(int argc, char** argv) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t hilos = parsear_hilos(argc, argv);

    std::printf("================================================================================\n");
    std::printf("NV3.2 (T3.6) -- Control Sellado Banco R2 (AR(1) phi=%.2f + Hawkes 2 signos)\n", PHI_R2);
    std::printf("  SEMILLA_CIEGA_NIVEL3 = %lluULL | N_SEM = %zu | N_CAL = %zu | N_TEST = %zu | hilos = %zu\n",
                static_cast<unsigned long long>(SEMILLA_CIEGA_NIVEL3), N_SEM, N_CAL, N_TEST, hilos);
    std::printf("  R2a (s=0..11): n_mismo=0.30, n_contr=0.15 | R2b (s=12..23): n_mismo=0.45, n_contr=0.25\n");
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

    std::printf("\n%3s %3s %3s %4s %4s %4s %8s %7s %4s | %4s %4s %4s %4s | %5s %5s %5s %4s %4s %4s %4s %5s\n",
                "sem", "grp", "cal", "eval", "t_ef", "mmin", "S", "p_bilat", "rech",
                "cRea", "kmed", "kpdt", "kevl",
                "nulB", "intB", "cortB", "parB", "solB", "ajPr", "reap", "escC");
    std::printf("--------------------------------------------------------------------------------------------------------------------\n");

    AcumGrupo g_r2a{}, g_r2b{}, g_tot{};
    for (const auto& r : resultados) {
        if (r.es_r2b) g_r2b.agregar(r);
        else g_r2a.agregar(r);
        g_tot.agregar(r);

        std::printf("%3zu %3s %3d %4d %4zu %4zu %+8.4f %7.4f %4d | %4zu %4zu %4zu %4zu | %5zu %5zu %5zu %4zu %4zu %4zu %4zu %5zu\n",
                    r.semilla_idx,
                    r.es_r2b ? "R2b" : "R2a",
                    r.calibracion_ok ? 1 : 0,
                    r.evaluada ? 1 : 0,
                    r.t_max_ef,
                    r.m_min,
                    r.estadistico_s,
                    r.p_bilateral,
                    r.rechazada ? 1 : 0,
                    r.cortadas_real,
                    r.k_med,
                    r.k_pdet_pos,
                    r.k_eval,
                    r.n_nulos_b,
                    r.intentos_b,
                    r.cortadas_b,
                    r.paradas_tempranas_b,
                    r.solapes_b,
                    r.con_ajena_previa_b,
                    r.con_reapertura_posterior_b,
                    r.escala_coincide_b);
    }

    auto imprimir_sub = [](const char* nom, const AcumGrupo& g) {
        const double sm = g.s_medio();
        const double ee = g.ee_s();
        const double zs = (ee > 0.0) ? (sm / ee) : 0.0;
        const double pdet = (g.intentos_b > 0) ? static_cast<double>(g.n_nulos_b) / static_cast<double>(g.intentos_b) : 0.0;
        std::printf("  %-3s : eval=%zu/%zu S<0=%zu/%zu rech=%zu/%zu | S_medio=%+.5f ee=%.5f z=%+.3f | "
                    "kmed=%zu kevl=%zu cRea=%zu | p_det=%.4f (%zu/%zu) cortB=%zu parB=%zu solB=%zu ajPr=%zu reap=%zu\n",
                    nom, g.n_eval, g.n_sem, g.n_s_neg, g.n_sem, g.n_rech, g.n_eval,
                    sm, ee, zs, g.k_med, g.k_eval, g.cort_real,
                    pdet, g.n_nulos_b, g.intentos_b, g.cortadas_b, g.paradas_b,
                    g.solapes_b, g.ajena_prev_b, g.reap_post_b);
    };

    const double S_medio = g_tot.s_medio();
    const double ee_S = g_tot.ee_s();
    const double z_S = (ee_S > 0.0) ? (S_medio / ee_S) : 0.0;
    const double tasa_rech = (g_tot.n_eval > 0)
        ? (static_cast<double>(g_tot.n_rech) / static_cast<double>(g_tot.n_eval))
        : 1.0;
    const double p_det_global = (g_tot.intentos_b > 0)
        ? (static_cast<double>(g_tot.n_nulos_b) / static_cast<double>(g_tot.intentos_b))
        : 0.0;
    const double surv_s3 = (g_tot.k_med > 0)
        ? (static_cast<double>(g_tot.k_pdet_pos) / static_cast<double>(g_tot.k_med))
        : 0.0;
    const double cob_s3 = (g_tot.k_med > 0)
        ? (static_cast<double>(g_tot.k_eval) / static_cast<double>(g_tot.k_med))
        : 0.0;
    const double tasa_escala_coincide = (g_tot.n_nulos_b > 0)
        ? (static_cast<double>(g_tot.escala_coinc_b) / static_cast<double>(g_tot.n_nulos_b))
        : 0.0;

    std::printf("\n================================================================================\n");
    std::printf("SUBGRUPOS NV3.2\n");
    imprimir_sub("R2a", g_r2a);
    imprimir_sub("R2b", g_r2b);
    std::printf("\nRESUMEN AGREGADO NV3.2 (tiempo=%.2f s)\n", segs);
    std::printf("  calibradas                  = %zu / %zu\n", g_tot.n_cal_ok, N_SEM);
    std::printf("  n_sem_eval                  = %zu / %zu  (ev_contraste_total=%zu)\n",
                g_tot.n_eval, N_SEM, g_tot.ev_ctr);
    std::printf("  semillas con S < 0          = %zu / %zu  (R2a=%zu/12, R2b=%zu/12)\n",
                g_tot.n_s_neg, g_tot.n_eval, g_r2a.n_s_neg, g_r2b.n_s_neg);
    std::printf("  S_medio                     = %+.6f | ee_S = %.6f | z_S = %+.6f\n", S_medio, ee_S, z_S);
    std::printf("  falsos rechazos (p < 0.05)  = %zu / %zu  (tasa_rech=%.6f)\n",
                g_tot.n_rech, g_tot.n_eval, tasa_rech);
    std::printf("  p_det_global                = %.6f  (%zu / %zu intentos_b)\n",
                p_det_global, g_tot.n_nulos_b, g_tot.intentos_b);
    std::printf("  surv_s3 (p_det_inicial > 0) = %.6f  (%zu / %zu ranuras)\n",
                surv_s3, g_tot.k_pdet_pos, g_tot.k_med);
    std::printf("  cob_s3 (n_nulos_b >= m_min) = %.6f  (%zu / %zu ranuras)\n",
                cob_s3, g_tot.k_eval, g_tot.k_med);
    std::printf("  cortadas_real_total         = %zu\n", g_tot.cort_real);
    std::printf("  cortadas_b_total            = %zu\n", g_tot.cortadas_b);
    std::printf("  paradas_tempranas_b_total   = %zu\n", g_tot.paradas_b);
    std::printf("  solapes_b_total             = %zu\n", g_tot.solapes_b);
    std::printf("  con_ajena_previa_b_total    = %zu\n", g_tot.ajena_prev_b);
    std::printf("  con_reapertura_post_b_total = %zu\n", g_tot.reap_post_b);
    std::printf("  Telemetria Nivel 4:\n");
    std::printf("    escala_coincide_b_total   = %zu / %zu  (tasa=%.6f)\n",
                g_tot.escala_coinc_b, g_tot.n_nulos_b, tasa_escala_coincide);
    std::printf("--------------------------------------------------------------------------------\n");

    const bool ok_sem = (g_tot.n_cal_ok == N_SEM && g_tot.n_eval == N_SEM);
    const bool ok_sneg = (g_tot.n_eval == N_SEM && g_tot.n_s_neg >= 7u && g_tot.n_s_neg <= 17u);
    const bool ok_rec = (g_tot.n_eval == N_SEM && g_tot.n_rech <= 3u);
    const bool ok_act = (g_tot.cortadas_b > 0u &&
                         g_tot.reap_post_b > 0u &&
                         g_tot.solapes_b == 0u &&
                         g_tot.ajena_prev_b > 0u);

    std::size_t subfallas = 0;
    auto reportar_sub = [&](const char* tag, bool ok, const char* criterio) {
        if (!ok) ++subfallas;
        std::printf("  [%s] %-7s : %s\n", ok ? "PASA" : "FALLA", tag, criterio);
    };
    reportar_sub("ok_sem",  ok_sem,  "n_cal_ok == 24 && n_sem_eval == 24");
    reportar_sub("ok_sneg", ok_sneg, "n_s_neg in [7, 17] sobre 24 semillas");
    reportar_sub("ok_rec",  ok_rec,  "falsos rechazos (p_bilat < 0.05) <= 3 / 24");
    reportar_sub("ok_act",  ok_act,  "cortadas_b_total > 0 && reap_post_b > 0 && solapes_b == 0 && ajena_prev_b > 0");

    const bool ok_nv3_2 = ok_sem && ok_sneg && ok_rec && ok_act;
    const std::size_t fallas = ok_nv3_2 ? 0u : 1u;
    std::printf("\n[%s] NV3.2 (subfallas=%zu)\n", ok_nv3_2 ? "PASA" : "FALLA", subfallas);
    std::printf("FALLAS NV3.2: %zu\n", fallas);
    std::fflush(stdout);
    return fallas == 0 ? 0 : 1;
}
