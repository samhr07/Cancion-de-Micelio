// ============================================================================
// test_nivel2_horizonte.cpp -- Bateria Sellada Completa de Nivel 2 (T2.8)
// Contrato: EXPECTATIVAS_NIVEL2.md (bloque ciego SEMILLA_CIEGA_NIVEL2 = 20261123ULL)
// Controles: NV2.1, NV2.2 (a, b, c), NV2.3 (a y b 1..6), NV2.4 (K10, l>256) y NV2.5.
// Cero cambios post-sello a umbrales, formulas o semillas.
// ============================================================================

#include <atomic>
#include <cstdlib>
#include <new>

namespace {
std::atomic<bool> g_contar_new{false};
std::atomic<std::size_t> g_news{0};
}  // namespace

void* operator new(std::size_t sz) {
    if (g_contar_new.load(std::memory_order_relaxed)) {
        g_news.fetch_add(1, std::memory_order_relaxed);
    }
    void* p = std::malloc(sz == 0 ? 1 : sz);
    if (!p) throw std::bad_alloc();
    return p;
}

void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <random>
#include <span>
#include <string_view>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;

constexpr std::uint64_t SEMILLA_CIEGA_NIVEL2 = 20261123ULL;

std::vector<double> simular_arp(
    std::size_t n, std::span<const double> phi, double mu, double sigma_e, std::uint64_t sem) {
    const std::size_t burn = 512;
    const std::size_t tot = burn + n;
    std::vector<double> u(tot, 0.0);
    std::vector<double> out(n, 0.0);
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> nd(0.0, sigma_e);
    for (std::size_t s = 0; s < tot; ++s) {
        double val = nd(rng);
        for (std::size_t i = 1; i <= phi.size() && i <= s; ++i) {
            val += phi[i - 1] * u[s - i];
        }
        u[s] = val;
        if (s >= burn) out[s - burn] = mu + val;
    }
    return out;
}

struct SerieArfima {
    std::vector<double> x;
    std::vector<double> eps;
};

SerieArfima simular_arfima1d0_con_eps(
    std::size_t n, double phi1, double d, double mu, double sigma_e, std::uint64_t sem, std::size_t burn = 2048) {
    const std::size_t tot = burn + n;
    std::vector<double> pi_inv(tot, 0.0);
    elw::coeficientes(-d, tot, pi_inv.data());
    std::vector<double> u(tot, 0.0);
    SerieArfima out{};
    out.x.assign(n, 0.0);
    out.eps.assign(n, 0.0);
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> nd(0.0, sigma_e);
    double u_prev = 0.0;
    for (std::size_t s = 0; s < tot; ++s) {
        const double e_t = nd(rng);
        const double u_t = phi1 * u_prev + e_t;
        u[s] = u_t;
        u_prev = u_t;
        if (s >= burn) {
            double acc = 0.0;
            for (std::size_t k = 0; k <= s; ++k) acc += pi_inv[k] * u[s - k];
            out.x[s - burn] = mu + acc;
            out.eps[s - burn] = e_t;
        }
    }
    return out;
}

template <class ObsLT>
typename ObsLT::Ranura tomar_ranura_evento(
    const ObsLT& lil,
    const typename ObsLT::SalidaPaso& out,
    bool& enc) {
    enc = false;
    if (out.abrio_ranura) {
        enc = true;
        return lil.ranuras()[out.indice_ranura_abierta];
    }
    if (out.ranuras_cerradas > 0) {
        enc = true;
        return lil.cerrada(out.ranuras_cerradas - 1);
    }
    return typename ObsLT::Ranura{};
}

}  // namespace

int main() {
    std::size_t fallas = 0;
    auto verificar = [&](const char* codigo, bool ok, const char* detalle) {
        std::printf("[%s] %s: %s\n", ok ? "PASA" : "FALLA", codigo, detalle);
        std::fflush(stdout);
        if (!ok) ++fallas;
    };

    std::printf("================================================================================\n");
    std::printf(" BATERIA SELLADA NIVEL 2 (EXPECTATIVAS_NIVEL2.md | SEMILLA = %lluULL)\n",
                static_cast<unsigned long long>(SEMILLA_CIEGA_NIVEL2));
    std::printf("================================================================================\n\n");

    // =========================================================================
    // NV2.1: Exactitud Markov en AR(p) (d = 0, p in {1, 4, 8})
    // =========================================================================
    {
        const std::vector<std::vector<double>> phis = {
            {0.70},
            {0.25, 0.05, 0.05, 0.15},
            {0.18, 0.12, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15},
        };
        const std::size_t ps[3] = {1, 4, 8};
        bool ok_nv2_1 = true;
        double worst_dif = 0.0;
        for (std::size_t idx = 0; idx < 3; ++idx) {
            const std::size_t p_esp = ps[idx];
            const auto z = simular_arp(4096, phis[idx], 12.5, 1.0, SEMILLA_CIEGA_NIVEL2 + 101ULL * p_esp);
            ConfiguracionModeloLineal cfg{};
            cfg.modo = ModoNuloLineal::AR_LEVINSON_DURBIN;
            cfg.alfa = 0.05;
            cfg.usar_enmienda_nivel2 = true;
            const auto m = ajustar({z}, cfg);
            std::vector<double> e_p(32, 0.0), e_2048(32, 0.0);
            innovaciones(m, z.data() + (2500 - m.n_hist_req), m.n_hist_req, z.data() + 2500, 32, e_p.data());
            innovaciones(m, z.data() + (2500 - 2048), 2048, z.data() + 2500, 32, e_2048.data());
            double max_d = 0.0;
            for (std::size_t j = 0; j < 32; ++j) max_d = std::max(max_d, std::abs(e_p[j] - e_2048[j]));
            worst_dif = std::max(worst_dif, max_d);
            const bool ok_i = m.valido && m.enmienda_nivel2 && m.p == p_esp && m.n_hist_req == p_esp &&
                              m.m_inv_arfima == 0 && m.m_sim_arfima == 0 && !m.burn_in_topado && max_d <= 1e-15;
            if (!ok_i) ok_nv2_1 = false;
            std::printf("  NV2.1 AR(%zu): p_fit=%zu n_req=%zu m_inv=%zu m_sim=%zu top=%d max|dif|=%.3e\n",
                        p_esp, m.p, m.n_hist_req, m.m_inv_arfima, m.m_sim_arfima, m.burn_in_topado ? 1 : 0, max_d);
        }
        char buf[192];
        std::snprintf(buf, sizeof(buf), "AR(1/4/8) Markov exacto con n_hist_req == p, worst_max|dif| = %.3e <= 1e-15", worst_dif);
        verificar("NV2.1", ok_nv2_1, buf);
    }

    // =========================================================================
    // NV2.2: Hosking analitico + pre-centrado mu_serie + ida/vuelta + precision cola K=512
    // =========================================================================
    {
        const double z_q = est::norm_ppf(0.975);
        const double eps_q = z_q / std::sqrt(4096.0);
        const double eps_q2 = eps_q * eps_q;
        const double d_tabla[6] = {-0.40, -0.35, -0.25, +0.25, +0.35, +0.40};
        const std::size_t m_tabla[6] = {338, 211, 82, 85, 233, 418};
        bool ok_a = true;
        for (std::size_t i = 0; i < 6; ++i) {
            const std::size_t m_calc = detalle_lineal::horizonte_inv_arfima_hosking(d_tabla[i], eps_q2, 4096);
            if (m_calc != m_tabla[i]) ok_a = false;
        }

        const double ds[3] = {-0.35, +0.25, +0.35};
        const std::size_t m_esps[3] = {211, 85, 233};
        const bool top_esps[3] = {false, true, true};
        bool ok_b = true;
        bool ok_c = true;
        double max_dif_tras_all = 0.0;
        double max_rt_all = 0.0;
        double max_rmse_ratio = 0.0;

        for (std::size_t idx = 0; idx < 3; ++idx) {
            const double d = ds[idx];
            const auto s50 = simular_arfima1d0_con_eps(4096, 0.40, d, 50.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 1009ULL * idx, 2048);
            const std::vector<double>& z50 = s50.x;
            std::vector<double> z1050 = z50;
            for (double& v : z1050) v += 1000.0;

            ConfiguracionModeloLineal cfg{};
            cfg.modo = ModoNuloLineal::ARFIMA_WHITTLE;
            cfg.alfa = 0.05;
            cfg.d_externa = d;
            cfg.usar_enmienda_nivel2 = true;
            const auto m50 = ajustar({z50}, cfg);
            const auto m1050 = ajustar({z1050}, cfg);

            if (!m50.valido || !m1050.valido || m50.m_inv_arfima != m_esps[idx] ||
                m50.n_hist_req != m50.p + m50.m_inv_arfima || m50.burn_in_topado != top_esps[idx]) {
                ok_a = false;
            }

            // (b) Traslacion + ida y vuelta
            std::vector<double> e50(64, 0.0), e1050(64, 0.0);
            const std::size_t nh50 = std::min<std::size_t>(m50.n_hist_req, 2500u);
            const std::size_t nh1050 = std::min<std::size_t>(m1050.n_hist_req, 2500u);
            innovaciones(m50, z50.data() + (2500 - nh50), nh50, z50.data() + 2500, 64, e50.data());
            innovaciones(m1050, z1050.data() + (2500 - nh1050), nh1050, z1050.data() + 2500, 64, e1050.data());
            double dif_tras = 0.0;
            for (std::size_t j = 0; j < 64; ++j) dif_tras = std::max(dif_tras, std::abs(e1050[j] - e50[j]));
            max_dif_tras_all = std::max(max_dif_tras_all, dif_tras);

            auto rt = [&](const ModeloLinealNulo& mod, const std::vector<double>& z) {
                constexpr std::size_t N_RT = 512;
                std::vector<double> e_rt(N_RT, 0.0), x_rec(N_RT, 0.0);
                std::vector<double> trab(tamano_trabajo_simular(mod, N_RT), 0.0);
                innovaciones(mod, nullptr, 0, z.data(), N_RT, e_rt.data());
                std::mt19937_64 rng_rt(SEMILLA_CIEGA_NIVEL2 + 99ULL + idx);
                simular(mod, rng_rt, e_rt.data(), 0, N_RT, x_rec.data(), N_RT, trab.data(), trab.size());
                double mx = 0.0;
                for (std::size_t t = 0; t < N_RT; ++t) mx = std::max(mx, std::abs(x_rec[t] - z[t]));
                return mx;
            };
            const double rt50 = rt(m50, z50);
            const double rt1050 = rt(m1050, z1050);
            max_rt_all = std::max({max_rt_all, rt50, rt1050});
            if (dif_tras > 1e-11 || rt50 > 1e-12 || rt1050 > 1e-12) ok_b = false;

            // (c) Precision de cola sobre K = 512 realizaciones estacionarias independientes (oracle params + n_hist_exact = m50.n_hist_req)
            ModeloLinealNulo m_eval{};
            m_eval.valido = true;
            m_eval.modo_usado = ModoNuloLineal::ARFIMA_WHITTLE;
            m_eval.enmienda_nivel2 = true;
            m_eval.p = 1;
            m_eval.phi = {0.40};
            m_eval.d = d;
            m_eval.mu_serie = 50.0;
            m_eval.mu = 0.0;
            m_eval.sigma2_innov = 1.0;

            constexpr std::size_t K_MC = 512;
            constexpr std::size_t BURN_MC = 4096;
            const std::size_t n_hist_exact = std::min<std::size_t>(m50.n_hist_req, 2048u);
            const std::size_t n_hist_m8 = m_eval.p + 8u;
            const std::size_t n_need = std::max(n_hist_exact, n_hist_m8) + 1u;

            std::vector<double> psi(BURN_MC + 1, 0.0);
            elw::coeficientes(-d, BURN_MC + 1, psi.data());
            std::mt19937_64 rng_mc(SEMILLA_CIEGA_NIVEL2 + 1009ULL * idx + 77777ULL);
            std::normal_distribution<double> n01(0.0, 1.0);
            std::vector<double> u_mc(BURN_MC + 1, 0.0);
            std::vector<double> z_tail(n_need, 0.0);

            double sse_ex = 0.0;
            double sse_m8 = 0.0;
            for (std::size_t k_rep = 0; k_rep < K_MC; ++k_rep) {
                double u_prev = 0.0;
                double eps_last = 0.0;
                for (std::size_t s = 0; s <= BURN_MC; ++s) {
                    eps_last = n01(rng_mc);
                    const double u_t = 0.40 * u_prev + eps_last;
                    u_mc[s] = u_t;
                    u_prev = u_t;
                }
                for (std::size_t j = 0; j < n_need; ++j) {
                    const std::size_t t_idx = (BURN_MC + 1u - n_need) + j;
                    double acc = 0.0;
                    for (std::size_t k = 0; k <= t_idx; ++k) acc += psi[k] * u_mc[t_idx - k];
                    z_tail[j] = 50.0 + acc;
                }
                double e_ex = NaN, e_8 = NaN;
                const double z_obj = z_tail.back();
                innovaciones(m_eval, &z_tail[n_need - 1u - n_hist_exact], n_hist_exact, &z_obj, 1, &e_ex);
                innovaciones(m_eval, &z_tail[n_need - 1u - n_hist_m8], n_hist_m8, &z_obj, 1, &e_8);
                const double de_ex = e_ex - eps_last;
                const double de_8 = e_8 - eps_last;
                sse_ex += de_ex * de_ex;
                sse_m8 += de_8 * de_8;
            }
            const double rmse_ex = std::sqrt(sse_ex / static_cast<double>(K_MC));
            const double rmse_m8 = std::sqrt(sse_m8 / static_cast<double>(K_MC));
            max_rmse_ratio = std::max(max_rmse_ratio, rmse_ex / eps_q);
            if (!(rmse_ex <= 1.15 * eps_q && rmse_ex < rmse_m8)) ok_c = false;
            std::printf("  NV2.2(d=%+.2f): m_inv=%zu top=%d dif_tras=%.3e rt=%.3e | RMSE_exact=%.5f (%.2f*epsQ) < RMSE_M8=%.5f\n",
                        d, m50.m_inv_arfima, m50.burn_in_topado ? 1 : 0, dif_tras, std::max(rt50, rt1050),
                        rmse_ex, rmse_ex / eps_q, rmse_m8);
        }
        char buf[220];
        std::snprintf(buf, sizeof(buf),
                      "Analitico=%d Trasl<=%.2e RT<=%.2e Cola<=%.2f*epsQ (lim=1.15)",
                      ok_a ? 1 : 0, max_dif_tras_all, max_rt_all, max_rmse_ratio);
        verificar("NV2.2", ok_a && ok_b && ok_c, buf);
    }

    // =========================================================================
    // NV2.3a: Supervivencia Z1b con q=2.0% en rafagas r=2 (ARFIMA d=+0.35, n_hist_req=234, 64 ventanas)
    // =========================================================================
    {
        using ObsS = ObservadorSaltos<24, 64, 60>;
        using ObsL = ObservadorLilyth<ObsS, 8, 64>;
        const auto s_cal = simular_arfima1d0_con_eps(4096, 0.40, 0.35, 25.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 2300ULL, 2048);
        const auto mom_cal = est::momentos(std::vector<std::vector<double>>{s_cal.x});

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 1.0e5;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<ObsS>();
        od->calibrar(cd, s_cal.x);

        ConfiguracionLilyth cl_ito{};
        cl_ito.representacion = ModoRepresentacionLilyth::CRUDA;
        cl_ito.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
        cl_ito.d_whittle_cruda = 0.35;
        cl_ito.b_max_nulo = 2;
        cl_ito.usar_enmienda_nivel2 = true;
        cl_ito.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO;
        cl_ito.r_max_z1_piloto = 2;
        cl_ito.alfa_hueco_piloto = 0.01;
        cl_ito.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cl_ito.r_max_z2_piloto = 0;

        ConfiguracionLilyth cl_corte = cl_ito;
        cl_corte.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::CORTE_SUFIJO;
        cl_corte.r_max_z1_piloto = 0;

        auto lil_ito = std::make_unique<ObsL>();
        auto lil_limp = std::make_unique<ObsL>();
        auto lil_corte = std::make_unique<ObsL>();
        const bool cal_ok = lil_ito->calibrar(*od, cl_ito, cd, s_cal.x) &&
                            lil_limp->calibrar(*od, cl_ito, cd, s_cal.x) &&
                            lil_corte->calibrar(*od, cl_corte, cd, s_cal.x);
        const auto& mod_cal = lil_ito->ficha_nulo_cruda().modelo_lineal;
        const std::size_t n_req_esp = mod_cal.p + 233u;

        constexpr std::size_t N_SERIES = 16;
        constexpr std::size_t N_LEN = 1800;
        const std::array<std::size_t, 4> t_evals{450, 850, 1250, 1650};
        const double salto = 8.0 * mom_cal.sd;

        std::size_t n_ventanas = 0;
        std::size_t surv_ito = 0;
        std::size_t surv_corte = 0;
        std::size_t total_puentes_z1b = 0;
        double sse_ito_vs_eps = 0.0;
        double sse_ito_vs_limpio = 0.0;
        std::size_t n_rmse = 0;

        for (std::size_t s_idx = 0; s_idx < N_SERIES; ++s_idx) {
            const std::uint64_t sem_s = SEMILLA_CIEGA_NIVEL2 + 7919ULL * s_idx;
            const auto sop = simular_arfima1d0_con_eps(N_LEN, 0.40, 0.35, 25.0, 1.0, sem_s, 2048);

            // Generar mascara q=2% en rafagas r=2 (probabilidad 0.01 de iniciar pareja [t, t+1])
            std::vector<std::uint8_t> val(N_LEN, 1u);
            auto protegido = [&](std::size_t t) {
                for (std::size_t te : t_evals) {
                    if (t + 1u >= te && t <= te + 2u) return true;
                }
                return false;
            };
            std::mt19937_64 rng_mask(sem_s ^ 0xA5A5A5A5A5A5A5A5ULL);
            std::uniform_real_distribution<double> u01(0.0, 1.0);
            for (std::size_t t = 5; t + 2 < N_LEN;) {
                if (!protegido(t) && !protegido(t + 1) && u01(rng_mask) < 0.01) {
                    val[t] = 0u;
                    val[t + 1] = 0u;
                    t += 3;
                } else {
                    ++t;
                }
            }

            lil_ito->reiniciar_flujo_conservando_purga();
            lil_limp->reiniciar_flujo_conservando_purga();
            lil_corte->reiniciar_flujo_conservando_purga();

            std::size_t k_ev = 0;
            bool esperando_cierre_espejo = false;
            EventoAlarma ult_alarma{};
            for (std::size_t t = 0; t < N_LEN; ++t) {
                typename ObsS::Salida sal{};
                sal.t = static_cast<std::int64_t>(t);
                sal.niveles = 1;
                double x_in = sop.x[t];
                const bool es_ev = (k_ev < t_evals.size() && t == t_evals[k_ev]);
                if (es_ev) {
                    x_in += salto;
                    sal.nivel[0].alarma.hay = true;
                    sal.nivel[0].alarma.signo = 1;
                    sal.nivel[0].alarma.t = static_cast<std::int64_t>(t);
                    sal.nivel[0].alarma.t_cambio = static_cast<std::int64_t>(t);
                    sal.nivel[0].alarma.l = 1;
                    sal.nivel[0].alarma.g = 25.0;
                    sal.nivel[0].alarma.amplitud = salto;
                    sal.nivel[0].alarma.se_amplitud = 1.0;
                    ult_alarma = sal.nivel[0].alarma;
                    esperando_cierre_espejo = true;
                    ++k_ev;
                } else if (esperando_cierre_espejo && static_cast<std::int64_t>(t) >= ult_alarma.t + 5) {
                    sal.n_eventos = 1;
                    sal.eventos[0].canal = 0;
                    sal.eventos[0].escala_origen = 0;
                    sal.eventos[0].alarma = ult_alarma;
                    esperando_cierre_espejo = false;
                }
                const auto o_ito = lil_ito->paso(sal, x_in, val[t] != 0u, NaN, 0, 0);
                const auto o_lmp = lil_limp->paso(sal, x_in, true, NaN, 0, 0);
                const auto o_crt = lil_corte->paso(sal, x_in, val[t] != 0u, NaN, 0, 0);

                if (es_ev) {
                    ++n_ventanas;
                    bool e1 = false, e2 = false, e3 = false;
                    const auto r_i = tomar_ranura_evento(*lil_ito, o_ito, e1);
                    const auto r_l = tomar_ranura_evento(*lil_limp, o_lmp, e2);
                    const auto r_c = tomar_ranura_evento(*lil_corte, o_crt, e3);
                    if (e1 && r_i.n_hist_prev_disp_cruda == mod_cal.n_hist_req) ++surv_ito;
                    if (e3 && r_c.n_hist_prev_disp_cruda == mod_cal.n_hist_req) ++surv_corte;
                    total_puentes_z1b += r_i.n_huecos_puenteados_z1b_cruda;
                    const auto ei_i = lil_ito->e_iny_evento_n2(r_i, false);
                    const auto ei_l = lil_limp->e_iny_evento_n2(r_l, false);
                    if (e1 && e2 && !ei_i.empty() && !ei_l.empty()) {
                        const double e_hat_ito = ei_i[0] - salto;
                        const double e_hat_lmp = ei_l[0] - salto;
                        const double d_eps = e_hat_ito - sop.eps[t];
                        const double d_lmp = e_hat_ito - e_hat_lmp;
                        sse_ito_vs_eps += d_eps * d_eps;
                        sse_ito_vs_limpio += d_lmp * d_lmp;
                        ++n_rmse;
                    }
                }
            }
        }
        const double tasa_ito = static_cast<double>(surv_ito) / static_cast<double>(n_ventanas);
        const double tasa_crt = static_cast<double>(surv_corte) / static_cast<double>(n_ventanas);
        const double rmse_eps = std::sqrt(sse_ito_vs_eps / static_cast<double>(std::max<std::size_t>(1, n_rmse)));
        const double rmse_lmp = std::sqrt(sse_ito_vs_limpio / static_cast<double>(std::max<std::size_t>(1, n_rmse)));
        const bool ok_nv2_3a = cal_ok && mod_cal.n_hist_req == n_req_esp && n_req_esp == 234 &&
                               n_ventanas == 64 && tasa_ito >= 0.85 && total_puentes_z1b > 0 &&
                               rmse_eps <= 0.10 && tasa_crt < 0.35;
        std::printf("  NV2.3a (q=2%% r=2, n_req=%zu, %zu vent): surv_Ito=%.1f%% (%zu/%zu) puentes=%zu RMSE(vs eps)=%.4f (vs limpio=%.4f) | surv_Corte=%.1f%%\n",
                    mod_cal.n_hist_req, n_ventanas, 100.0 * tasa_ito, surv_ito, n_ventanas, total_puentes_z1b,
                    rmse_eps, rmse_lmp, 100.0 * tasa_crt);
        char buf[200];
        std::snprintf(buf, sizeof(buf),
                      "Ito=%.1f%% (>=85%%, puentes=%zu), RMSE=%.4f (<=0.10), Corte=%.1f%% (<35%%)",
                      100.0 * tasa_ito, total_puentes_z1b, rmse_eps, 100.0 * tasa_crt);
        verificar("NV2.3a", ok_nv2_3a, buf);
    }

    // =========================================================================
    // NV2.3b: Cortes duros (r=4, choque +12sigma, hist_seg_purgada_, Z1a, Z2, n_hist < p)
    // =========================================================================
    {
        using ObsS = ObservadorSaltos<24, 64, 60>;
        using ObsL = ObservadorLilyth<ObsS, 8, 64>;
        constexpr std::uint64_t SEM_3B = SEMILLA_CIEGA_NIVEL2 + 5003ULL;
        const auto s_cal = simular_arfima1d0_con_eps(1500, 0.40, 0.35, 25.0, 1.0, SEM_3B + 1ULL, 1024);
        const auto mom_c = est::momentos(std::vector<std::vector<double>>{s_cal.x});

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 1.0e6;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<ObsS>();
        od->calibrar(cd, s_cal.x);

        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
        cl.d_whittle_cruda = 0.35;
        cl.b_max_nulo = 2;
        cl.usar_enmienda_nivel2 = true;
        cl.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO;
        cl.r_max_z1_piloto = 2;
        cl.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cl.r_max_z2_piloto = 0;

        auto lil = std::make_unique<ObsL>();
        const bool cal_3b = lil->calibrar(*od, cl, cd, s_cal.x);
        const double sd_e = std::sqrt(lil->ficha_nulo_cruda().modelo_lineal.sigma2_innov);

        // (1) r = 4 > r_max = 2 -> exact d_cut
        bool ok_b1 = false;
        {
            auto sop = simular_arfima1d0_con_eps(560, 0.40, 0.35, 25.0, 1.0, SEM_3B + 10ULL, 1024);
            for (std::size_t t = 500; t < sop.x.size(); ++t) sop.x[t] += 10.0 * mom_c.sd;
            std::vector<std::uint8_t> val(sop.x.size(), 1u);
            for (std::size_t k = 460; k <= 463; ++k) val[k] = 0u;
            od->reiniciar_flujo();
            lil->reiniciar_flujo_conservando_purga();
            for (std::size_t t = 0; t < sop.x.size(); ++t) {
                auto s = od->paso(sop.x[t], NaN, val[t] != 0u, 0, 0);
                if (t < 495) for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
                const auto o = lil->paso(s, sop.x[t], val[t] != 0u, NaN, 0, 0);
                if (o.abrio_ranura) {
                    bool enc = false;
                    const auto r = tomar_ranura_evento(*lil, o, enc);
                    const std::int64_t t_ini_perf = r.t_apertura - static_cast<std::int64_t>(r.len_perfil_cruda_n2) + 1;
                    const std::size_t d_cut = static_cast<std::size_t>(t_ini_perf - 464);
                    ok_b1 = (r.n_hist_prev_disp_cruda == d_cut && r.hist_prev_truncada_cruda);
                    break;
                }
            }
        }

        // (2) Choque oculto +12*sigma en r=2 [620..621]
        bool ok_b2 = false;
        {
            auto sop = simular_arfima1d0_con_eps(760, 0.40, 0.35, 25.0, 1.0, SEM_3B + 20ULL, 1024);
            for (std::size_t t = 700; t < sop.x.size(); ++t) sop.x[t] += 10.0 * mom_c.sd;
            std::vector<std::uint8_t> val(sop.x.size(), 1u);
            val[620] = 0u;
            val[621] = 0u;
            sop.x[622] += 12.0 * sd_e;
            od->reiniciar_flujo();
            lil->reiniciar_flujo_conservando_purga();
            for (std::size_t t = 0; t < sop.x.size(); ++t) {
                auto s = od->paso(sop.x[t], NaN, val[t] != 0u, 0, 0);
                if (t < 695) for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
                const auto o = lil->paso(s, sop.x[t], val[t] != 0u, NaN, 0, 0);
                if (o.abrio_ranura) {
                    bool enc = false;
                    const auto r = tomar_ranura_evento(*lil, o, enc);
                    const std::int64_t t_ini_perf = r.t_apertura - static_cast<std::int64_t>(r.len_perfil_cruda_n2) + 1;
                    const std::size_t d_shock = static_cast<std::size_t>(t_ini_perf - 622);
                    ok_b2 = (r.n_hist_prev_disp_cruda == d_shock && r.n_huecos_puenteados_z1b_cruda == 0 &&
                             r.hist_prev_truncada_cruda);
                    break;
                }
            }
        }

        // (3) Cambio de hist_seg_purgada_ en PURGADA vs puente r=2 en CRUDA
        bool ok_b3 = false;
        {
            constexpr std::size_t BINS = 4;
            constexpr std::size_t N_NIV = 1500;
            const auto dy_cal = simular_arfima1d0_con_eps(N_NIV, 0.40, 0.35, 0.0, 1.0, SEM_3B + 30ULL, 1024);
            std::vector<double> x_cn(N_NIV, 100.0);
            std::vector<std::size_t> f_cn(N_NIV, 0);
            std::vector<std::int64_t> d_cn(N_NIV, 0);
            for (std::size_t t = 1; t < N_NIV; ++t) {
                f_cn[t] = t % BINS;
                d_cn[t] = static_cast<std::int64_t>(t / BINS);
                x_cn[t] = x_cn[t - 1] + (1.0 + 0.25 * static_cast<double>(f_cn[t])) * dy_cal.x[t];
            }
            ConfiguracionObservador cd_niv{};
            cd_niv.integracion = Integracion::NIVEL;
            cd_niv.purga = true;
            cd_niv.bins_fase = BINS;
            cd_niv.dias_minimos_purga = 5;
            cd_niv.arl0_pasos_nativos = 1.0e6;
            cd_niv.medir_tau50 = false;
            cd_niv.escalas_max = 1;
            auto od_n = std::make_unique<ObsS>();
            od_n->calibrar(cd_niv, x_cn, {}, {}, f_cn, d_cn);
            ConfiguracionLilyth cl_n = cl;
            cl_n.representacion = ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA;
            cl_n.d_whittle_purgada = 0.35;
            cl_n.d_whittle_cruda = 0.35;
            auto lil_n = std::make_unique<ObsL>();
            if (lil_n->calibrar(*od_n, cl_n, cd_niv, x_cn, {}, f_cn, d_cn)) {
                const auto dy_op = simular_arfima1d0_con_eps(600, 0.40, 0.35, 0.0, 1.0, SEM_3B + 31ULL, 1024);
                std::vector<double> x_op(600, 100.0);
                for (std::size_t t = 1; t < 600; ++t) {
                    const double amp = 1.0 + 0.25 * static_cast<double>(t % BINS);
                    x_op[t] = x_op[t - 1] + amp * (dy_op.x[t] + ((t >= 520) ? 10.0 : 0.0));
                }
                od_n->reiniciar_flujo();
                lil_n->reiniciar_flujo_conservando_purga();
                for (std::size_t t = 0; t < 600; ++t) {
                    const std::size_t fa = t % BINS;
                    const std::int64_t di = static_cast<std::int64_t>(t / BINS);
                    const bool ok = (t != 480);
                    auto s = od_n->paso(x_op[t], NaN, ok, fa, di);
                    if (t < 515) for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
                    const auto o = lil_n->paso(s, x_op[t], ok, NaN, fa, di);
                    if (o.abrio_ranura) {
                        bool enc = false;
                        const auto r = tomar_ranura_evento(*lil_n, o, enc);
                        const std::int64_t tip = r.t_apertura - static_cast<std::int64_t>(r.len_perfil_purgada_n2) + 1;
                        const std::size_t d_seg = static_cast<std::size_t>(tip - 482);
                        const std::size_t req_c = lil_n->ficha_nulo_cruda().modelo_lineal.n_hist_req;
                        ok_b3 = (r.n_hist_prev_disp_purgada == d_seg && r.hist_prev_truncada_purgada &&
                                 r.n_hist_prev_disp_cruda == req_c && r.n_huecos_puenteados_z1b_cruda == 1 &&
                                 !r.hist_prev_truncada_cruda);
                        break;
                    }
                }
            }
        }

        // (4) Z1a seguro (t_cambio - 1): nunca usa t_cambio como borde derecho del puente
        bool ok_b4 = false;
        {
            const auto sop = simular_arfima1d0_con_eps(520, 0.40, 0.35, 25.0, 1.0, SEM_3B + 40ULL, 1024);
            const double dA = 8.0 * mom_c.sd;
            const double dB = 22.0 * mom_c.sd;
            auto correr_z1a = [&](double dv, double& htail, double& e0) {
                od->reiniciar_flujo();
                lil->reiniciar_flujo_conservando_purga();
                for (std::size_t t = 0; t < sop.x.size(); ++t) {
                    const bool ok = (t != 499);
                    const double xt = sop.x[t] + ((t >= 500) ? dv : 0.0);
                    auto s = od->paso(xt, NaN, ok, 0, 0);
                    for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
                    if (t == 500) {
                        s.nivel[0].alarma.hay = true;
                        s.nivel[0].alarma.signo = 1;
                        s.nivel[0].alarma.t = 500;
                        s.nivel[0].alarma.t_cambio = 500;
                        s.nivel[0].alarma.l = 1;
                        s.nivel[0].alarma.g = 25.0;
                        s.nivel[0].alarma.amplitud = dv;
                        s.nivel[0].alarma.se_amplitud = 1.0;
                    }
                    const auto o = lil->paso(s, xt, ok, NaN, 0, 0);
                    if (o.abrio_ranura) {
                        bool enc = false;
                        const auto rr = tomar_ranura_evento(*lil, o, enc);
                        htail = rr.n_hist_prev_cruda > 0 ? rr.hist_prev_cruda_est[rr.n_hist_prev_cruda - 1] : NaN;
                        const auto ei = lil->e_iny_evento_n2(rr, false);
                        e0 = ei.empty() ? NaN : ei[0];
                        return true;
                    }
                }
                return false;
            };
            double hA = NaN, hB = NaN, eA = NaN, eB = NaN;
            if (correr_z1a(dA, hA, eA) && correr_z1a(dB, hB, eB)) {
                ok_b4 = (std::abs(hA - hB) == 0.0 && std::abs((eB - eA) - (dB - dA)) <= 1e-12);
            }
        }

        // (5) Veto Z2 (forzar_via_b true vs false con via (a) disponible)
        bool ok_b5 = false;
        {
            const std::array<double, 1> phi1{0.40};
            const auto xc = simular_arp(6000, phi1, 0.0, 1.0, SEM_3B + 50ULL);
            ConfiguracionObservador cdz{};
            cdz.integracion = Integracion::FLUJO;
            cdz.arl0_pasos_nativos = 80.0;
            cdz.medir_tau50 = false;
            cdz.escalas_max = 1;
            auto odz = std::make_unique<ObsS>();
            odz->calibrar(cdz, xc);

            auto correr_z2 = [&](bool forzar_b) {
                ConfiguracionLilyth clz{};
                clz.representacion = ModoRepresentacionLilyth::CRUDA;
                clz.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
                clz.b_max_nulo = 16;
                clz.usar_enmienda_nivel2 = true;
                clz.forzar_via_b = forzar_b;
                clz.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO;
                clz.r_max_z1_piloto = 2;
                clz.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
                clz.r_max_z2_piloto = 0;
                auto lilz = std::make_unique<ObsL>();
                (void)lilz->calibrar(*odz, clz, cdz, xc);
                auto xop = simular_arp(620, phi1, 0.0, 1.0, SEM_3B + 51ULL);
                for (std::size_t t = 500; t <= 505; ++t) xop[t] += 5.0;
                odz->reiniciar_flujo();
                lilz->reiniciar_flujo_conservando_purga();
                for (std::size_t t = 0; t < xop.size(); ++t) {
                    auto s = odz->paso(xop[t], NaN, true, 0, 0);
                    for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
                    if (t == 505) {
                        s.nivel[0].alarma.hay = true;
                        s.nivel[0].alarma.signo = 1;
                        s.nivel[0].alarma.t = 505;
                        s.nivel[0].alarma.t_cambio = 500;
                        s.nivel[0].alarma.l = 6;
                        s.nivel[0].alarma.g = 2.0;
                        s.nivel[0].alarma.amplitud = 5.0;
                        s.nivel[0].alarma.se_amplitud = 1.0;
                    }
                    const std::uint8_t fl = (t == 502) ? ark::dataset::FLAG_PUENTE_ITO : 0u;
                    (void)lilz->paso(s, xop[t], true, NaN, 0, 0, NaN, {}, fl);
                }
                lilz->completar_nulos();
                return lilz->cerradas_historico().empty() ? ObsL::Ranura{} : lilz->cerradas_historico().front();
            };
            const auto rb = correr_z2(true);
            const auto ra = correr_z2(false);
            ok_b5 = (rb.perfil_corrupto_z2_cruda && rb.cruda.nulo_onsager.sin_nulo_b && !rb.cruda.nulo_onsager.evaluado &&
                     ra.perfil_corrupto_z2_cruda && ra.cruda.desenlace != DesenlaceRanura::NO_MEDIDA_PERFIL_CORRUPTO &&
                     ra.cruda.nulo_onsager.evaluado && !ra.cruda.nulo_onsager.nulo_por_inyeccion &&
                     !ra.cruda.nulo_onsager.sin_nulo_b);
        }

        // (6) Veto n_hist_prev_disp < mod.p en AR(4)
        bool ok_b6 = false;
        {
            const std::array<double, 4> phi4{0.25, 0.05, 0.05, 0.15};
            const auto x4 = simular_arp(4096, phi4, 0.0, 1.0, SEM_3B + 60ULL);
            const auto m4 = est::momentos(std::vector<std::vector<double>>{x4});
            od->calibrar(cd, x4);
            ConfiguracionLilyth cl4{};
            cl4.representacion = ModoRepresentacionLilyth::CRUDA;
            cl4.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
            cl4.b_max_nulo = 2;
            cl4.usar_enmienda_nivel2 = true;
            cl4.forzar_via_b = true;
            cl4.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::CORTE_SUFIJO;
            auto lil4 = std::make_unique<ObsL>();
            (void)lil4->calibrar(*od, cl4, cd, x4);
            auto xop4 = simular_arp(620, phi4, 0.0, 1.0, SEM_3B + 61ULL);
            for (std::size_t t = 500; t <= 503; ++t) xop4[t] += 8.0 * m4.sd;
            od->reiniciar_flujo();
            lil4->reiniciar_flujo_conservando_purga();
            for (std::size_t t = 0; t < xop4.size(); ++t) {
                const bool ok = (t != 497);
                auto s = od->paso(xop4[t], NaN, ok, 0, 0);
                for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
                if (t == 500) {
                    s.nivel[0].alarma.hay = true;
                    s.nivel[0].alarma.signo = 1;
                    s.nivel[0].alarma.t = 500;
                    s.nivel[0].alarma.t_cambio = 500;
                    s.nivel[0].alarma.l = 1;
                    s.nivel[0].alarma.g = 20.0;
                    s.nivel[0].alarma.amplitud = 8.0 * m4.sd;
                    s.nivel[0].alarma.se_amplitud = 1.0;
                }
                (void)lil4->paso(s, xop4[t], ok, NaN, 0, 0);
            }
            lil4->completar_nulos();
            const std::size_t cer4 = lil4->cerradas_historico().size();
            const auto ref4 = cer4 > 0 ? lil4->cerradas_historico().front().cruda.nulo_onsager : ReferenciaNuloEvento{};
            std::printf("  NV2.3b(6) detalle: p_fit=%zu cer=%zu n_disp=%zu sin_b=%d eval=%d trunc=%d motivo='%s'\n",
                        lil4->ficha_nulo_cruda().modelo_lineal.p, cer4, ref4.n_hist_prev_disp,
                        ref4.sin_nulo_b ? 1 : 0, ref4.evaluado ? 1 : 0, ref4.hist_prev_truncada_en_corte ? 1 : 0,
                        ref4.motivo);
            if (cer4 > 0) {
                const std::size_t p_fit4 = lil4->ficha_nulo_cruda().modelo_lineal.p;
                ok_b6 = (p_fit4 >= 4 && ref4.n_hist_prev_disp == 2 && ref4.n_hist_prev_disp < p_fit4 &&
                         ref4.sin_nulo_b && !ref4.evaluado && ref4.hist_prev_truncada_en_corte &&
                         std::string_view(ref4.motivo) == "via (b): historia contigua menor que p");
            }
        }

        std::printf("  NV2.3b subcontroles: (1)r=4:%d (2)choque:%d (3)seg_purgada:%d (4)Z1a:%d (5)Z2:%d (6)<p:%d\n",
                    ok_b1 ? 1 : 0, ok_b2 ? 1 : 0, ok_b3 ? 1 : 0, ok_b4 ? 1 : 0, ok_b5 ? 1 : 0, ok_b6 ? 1 : 0);
        verificar("NV2.3b", cal_3b && ok_b1 && ok_b2 && ok_b3 && ok_b4 && ok_b5 && ok_b6,
                  "6/6 cortes duros y vetos Z1b/Z1a/Z2/<p exactos");
    }

    // =========================================================================
    // NV2.4: K10 Cero-Heap en paso() y Memoria Adaptable con l=320 > 256
    // =========================================================================
    {
        using ObsS3 = ObservadorSaltos<24, 120, 120>;
        using ObsL3 = ObservadorLilyth<ObsS3, 8, 64>;
        const std::array<double, 1> phi1{0.40};
        const auto x_j3 = simular_arp(10000, phi1, 10.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 9001ULL);
        const auto s_arf3 = simular_arfima1d0_con_eps(1400, 0.40, 0.35, 10.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 9002ULL, 1024);

        ConfiguracionObservador cd3{};
        cd3.integracion = Integracion::FLUJO;
        cd3.arl0_pasos_nativos = 450.0;
        cd3.medir_tau50 = false;
        cd3.escalas_max = 3;
        auto od3 = std::make_unique<ObsS3>();
        od3->calibrar(cd3, x_j3);

        ConfiguracionLilyth cl3{};
        cl3.representacion = ModoRepresentacionLilyth::CRUDA;
        cl3.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
        cl3.d_whittle_cruda = 0.35;
        cl3.b_max_nulo = 2;
        cl3.usar_enmienda_nivel2 = true;
        cl3.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO;
        cl3.r_max_z1_piloto = 2;
        cl3.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cl3.r_max_z2_piloto = 0;

        auto lil3 = std::make_unique<ObsL3>();
        const bool cal3 = lil3->calibrar(*od3, cl3, cd3, s_arf3.x);
        const auto sq3 = simular_arfima1d0_con_eps(1300, 0.40, 0.35, 10.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 9003ULL, 1024);

        od3->reiniciar_flujo();
        lil3->reiniciar_flujo_conservando_purga();
        ObsL3::Ranura r3{};
        bool abrio3 = false;

        g_news.store(0, std::memory_order_relaxed);
        g_contar_new.store(true, std::memory_order_release);
        for (std::size_t t = 0; t < sq3.x.size(); ++t) {
            auto s = od3->paso(sq3.x[t], NaN, true, 0, 0);
            for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            const double x_t = sq3.x[t] + ((t >= 800 && t <= 1119) ? 4.0 : 0.0);
            if (t == 1119) {
                s.nivel[2].alarma.hay = true;
                s.nivel[2].alarma.signo = 1;
                s.nivel[2].alarma.t = 1119;
                s.nivel[2].alarma.t_cambio = 800;
                s.nivel[2].alarma.l = 80;
                s.nivel[2].alarma.g = 20.0;
                s.nivel[2].alarma.amplitud = 4.0;
                s.nivel[2].alarma.se_amplitud = 0.5;
            }
            const auto o = lil3->paso(s, x_t, true, NaN, 0, 0);
            if (o.abrio_ranura) {
                r3 = lil3->ranuras()[o.indice_ranura_abierta];
                abrio3 = true;
            }
        }
        g_contar_new.store(false, std::memory_order_release);
        const std::size_t news_paso = g_news.load(std::memory_order_relaxed);

        const auto pf3 = lil3->perfil_evento_n2(r3, false);
        const auto ei3 = lil3->e_iny_evento_n2(r3, false);
        const auto& mod3 = lil3->ficha_nulo_cruda().modelo_lineal;
        std::vector<double> e_ref(320, 0.0);
        if (r3.n_hist_prev_disp_cruda > 0 && pf3.size() == 320) {
            const std::size_t nh_ref = std::min<std::size_t>(r3.n_hist_prev_disp_cruda, 800u);
            innovaciones(mod3, &sq3.x[800 - nh_ref], nh_ref, pf3.data(), 320, e_ref.data());
        }
        double max_dif_ei = 0.0;
        for (std::size_t i = 0; i < ei3.size() && i < e_ref.size(); ++i) {
            max_dif_ei = std::max(max_dif_ei, std::abs(ei3[i] - e_ref[i]));
        }
        ObsL3 lil_copia(*lil3);
        const bool ok_nv2_4 = cal3 && od3->escalas() == 3 && abrio3 && lil3->cap_perfil_adaptable() == 480 &&
                              news_paso == 0 && pf3.size() == 320 && ei3.size() == 320 &&
                              max_dif_ei <= 1e-13 && lil3->contadores().arena_perfiles_desbordada == 0 &&
                              sizeof(ObsL3::Ranura) < 16384 &&
                              lil3->arena_perfiles_usados() > 0 &&
                              lil_copia.arena_perfiles_capacidad() == 0 && lil_copia.arena_perfiles_usados() == 0;
        char buf[200];
        std::snprintf(buf, sizeof(buf),
                      "news=%zu cap=%zu len_pf=%zu len_ei=%zu dif_einy=%.2e sizeof(Ranura)=%zu copia_cap=%zu",
                      news_paso, lil3->cap_perfil_adaptable(), pf3.size(), ei3.size(), max_dif_ei,
                      sizeof(ObsL3::Ranura), lil_copia.arena_perfiles_capacidad());
        verificar("NV2.4", ok_nv2_4, buf);
    }

    // =========================================================================
    // NV2.5 (en C++): Opt-in Regla 6 por omision == false + 0 arena reservada
    // =========================================================================
    {
        ConfiguracionModeloLineal cm_def{};
        ConfiguracionLilyth cl_def{};
        cl_def.t_max_efectivo = 32;
        using ObsS = ObservadorSaltos<24, 64, 60>;
        using ObsL = ObservadorLilyth<ObsS, 8, 64>;
        const std::array<double, 1> phi1{0.40};
        const auto x_def = simular_arp(2000, phi1, 0.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 9500ULL);
        ConfiguracionObservador cd_def{};
        cd_def.integracion = Integracion::FLUJO;
        cd_def.arl0_pasos_nativos = 450.0;
        cd_def.medir_tau50 = false;
        cd_def.escalas_max = 1;
        auto od_def = std::make_unique<ObsS>();
        od_def->calibrar(cd_def, x_def);
        auto lil_def = std::make_unique<ObsL>();
        const bool cfg_ok = lil_def->configurar(cl_def, cd_def, *od_def, x_def);
        const bool ok_nv2_5 = !cm_def.usar_enmienda_nivel2 && !cl_def.usar_enmienda_nivel2 && !cl_def.usar_enmienda_nivel3 &&
                              cfg_ok && lil_def->arena_perfiles_capacidad() == 0;
        verificar("NV2.5", ok_nv2_5, "usar_enmienda_nivel2=false por omision y arena_cap=0 (verificacion bit a bit CN1.1 externa)");
    }

    std::printf("\n--------------------------------------------------------------------------------\n");
    std::printf("FALLAS NIVEL 2: %zu\n", fallas);
    return fallas == 0 ? 0 : 1;
}
