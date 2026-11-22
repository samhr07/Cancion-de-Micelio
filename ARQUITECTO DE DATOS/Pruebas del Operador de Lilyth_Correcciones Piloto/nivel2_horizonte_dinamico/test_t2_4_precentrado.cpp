// ============================================================================
// test_t2_4_precentrado.cpp -- Micro-test de verificacion de T2.4:
//   1. AR(p) (p in {1, 4, 8}) con usar_enmienda_nivel2 = true: n_hist_req == p,
//      m_inv_arfima == 0, m_sim_arfima == 0, y exactitud bit a bit en innovaciones().
//   2. ARFIMA(1, d, 0) (d in {-0.35, +0.25, +0.35}) con pre-centrado mu_serie:
//      horizontes Hosking {211, 85, 233}, burn_in_topado, invariancia ante traslacion
//      mu_z = 50 -> 1050 (<= 1e-11), e ida y vuelta innovaciones -> simular (<= 1e-12).
// ============================================================================

#include "ark/lilyth/modelo_lineal_nulo.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

namespace {

using namespace ark::lilyth;

std::vector<double> simular_arp(std::size_t n, const std::vector<double>& phi, double mu, std::uint64_t sem) {
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> n01(0.0, 1.0);
    const std::size_t burn = 512;
    const std::size_t tot = burn + n;
    std::vector<double> u(tot, 0.0);
    for (std::size_t t = 0; t < tot; ++t) {
        double acc = n01(rng);
        for (std::size_t i = 1; i <= phi.size() && i <= t; ++i) {
            acc += phi[i - 1] * u[t - i];
        }
        u[t] = acc;
    }
    std::vector<double> out(n, 0.0);
    for (std::size_t t = 0; t < n; ++t) out[t] = mu + u[burn + t];
    return out;
}

std::vector<double> simular_arfima1d0(std::size_t n, double phi1, double d, double mu, std::uint64_t sem) {
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> n01(0.0, 1.0);
    const std::size_t burn = 2048;
    const std::size_t tot = burn + n;
    std::vector<double> u(tot, 0.0);
    for (std::size_t t = 0; t < tot; ++t) {
        u[t] = n01(rng) + (t >= 1 ? phi1 * u[t - 1] : 0.0);
    }
    std::vector<double> psi(tot, 0.0);
    ark::metrologia::elw::coeficientes(-d, tot, psi.data());
    std::vector<double> out(n, 0.0);
    for (std::size_t t = 0; t < n; ++t) {
        const std::size_t idx = burn + t;
        double s = 0.0;
        for (std::size_t k = 0; k <= idx; ++k) s += psi[k] * u[idx - k];
        out[t] = mu + s;
    }
    return out;
}

}  // namespace

int main() {
    std::size_t fallas = 0;
    auto check = [&](bool cond, const char* msg) {
        if (!cond) {
            ++fallas;
            std::printf("  [FALLA] %s\n", msg);
        }
    };

    constexpr std::uint64_t SEM_BASE = 20261123ULL;

    std::printf("================================================================================\n");
    std::printf(" TEST T2.4: PRE-CENTRADO mu_serie E INVERSION EXACTA EN modelo_lineal_nulo.hpp\n");
    std::printf("================================================================================\n\n");

    // 1. AR(p), p in {1, 4, 8}
    std::printf("--- Check 1: Exactitud Markov AR(p) con usar_enmienda_nivel2 = true ---\n");
    const std::vector<std::vector<double>> phis_ar = {
        {0.70},
        {0.25, 0.05, 0.05, 0.15},
        {0.18, 0.12, 0.09, 0.08, 0.07, 0.06, 0.05, 0.15},
    };
    const std::size_t p_esps[3] = {1, 4, 8};

    for (std::size_t idx = 0; idx < 3; ++idx) {
        const std::size_t p_esp = p_esps[idx];
        const auto z = simular_arp(4096, phis_ar[idx], 12.5, SEM_BASE + 101ULL * p_esp);
        ConfiguracionModeloLineal cfg{};
        cfg.modo = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cfg.alfa = 0.05;
        cfg.usar_enmienda_nivel2 = true;

        const auto m = ajustar({z}, cfg);
        std::vector<double> e_p(32, 0.0), e_2048(32, 0.0);
        const std::size_t t0 = 2500;
        innovaciones(m, z.data() + (t0 - m.n_hist_req), m.n_hist_req, z.data() + t0, 32, e_p.data());
        innovaciones(m, z.data() + (t0 - 2048), 2048, z.data() + t0, 32, e_2048.data());
        double max_dif = 0.0;
        for (std::size_t j = 0; j < 32; ++j) max_dif = std::max(max_dif, std::abs(e_p[j] - e_2048[j]));

        std::printf("  AR(%zu): valido=%d, p_fit=%zu, n_hist_req=%zu, m_inv=%zu, m_sim=%zu, top=%d, max|dif|=%.3e\n",
                    p_esp, int(m.valido), m.p, m.n_hist_req, m.m_inv_arfima, m.m_sim_arfima, int(m.burn_in_topado), max_dif);
        check(m.valido && m.enmienda_nivel2, "AR ajustado valido con enmienda Nivel 2");
        check(m.p == p_esp, "Quenouille recupera p exacto");
        check(m.n_hist_req == m.p, "n_hist_req == p");
        check(m.m_inv_arfima == 0 && m.m_sim_arfima == 0 && !m.burn_in_topado, "ARFIMA campos en 0 para AR(p)");
        check(max_dif <= 1e-15, "Exactitud bit a bit innovaciones(p) vs innovaciones(2048)");
    }

    // 2. ARFIMA(1, d, 0) con mu_serie y traslacion +1000
    std::printf("\n--- Check 2: ARFIMA(1, d, 0) mu_serie + Hosking + Ida y Vuelta ---\n");
    const double ds[3] = {-0.35, +0.25, +0.35};
    const std::size_t m_inv_esps[3] = {211, 85, 233};
    const bool top_esps[3] = {false, true, true};

    for (std::size_t idx = 0; idx < 3; ++idx) {
        const double d = ds[idx];
        const auto z50 = simular_arfima1d0(4096, 0.40, d, 50.0, SEM_BASE + 1009ULL * idx);
        std::vector<double> z1050 = z50;
        for (double& v : z1050) v += 1000.0;

        ConfiguracionModeloLineal cfg{};
        cfg.modo = ModoNuloLineal::ARFIMA_WHITTLE;
        cfg.alfa = 0.05;
        cfg.d_externa = d;
        cfg.usar_enmienda_nivel2 = true;

        const auto m50 = ajustar({z50}, cfg);
        const auto m1050 = ajustar({z1050}, cfg);

        const std::size_t l_perf = 64;
        const std::size_t t0 = 2500;
        std::vector<double> e50(l_perf, 0.0), e1050(l_perf, 0.0);
        innovaciones(m50, z50.data() + (t0 - m50.n_hist_req), m50.n_hist_req, z50.data() + t0, l_perf, e50.data());
        innovaciones(m1050, z1050.data() + (t0 - m1050.n_hist_req), m1050.n_hist_req, z1050.data() + t0, l_perf, e1050.data());
        double max_dif_tras = 0.0;
        for (std::size_t j = 0; j < l_perf; ++j) max_dif_tras = std::max(max_dif_tras, std::abs(e1050[j] - e50[j]));

        auto check_rt = [&](const ModeloLinealNulo& m, const std::vector<double>& z) {
            constexpr std::size_t N_RT = 512;
            std::vector<double> e_rt(N_RT, 0.0), x_rec(N_RT, 0.0);
            std::vector<double> trab(tamano_trabajo_simular(m, N_RT), 0.0);
            innovaciones(m, nullptr, 0, z.data(), N_RT, e_rt.data());
            std::mt19937_64 rng(999ULL);
            simular(m, rng, e_rt.data(), 0, N_RT, x_rec.data(), N_RT, trab.data(), trab.size());
            double mx = 0.0;
            for (std::size_t t = 0; t < N_RT; ++t) mx = std::max(mx, std::abs(x_rec[t] - z[t]));
            return mx;
        };
        const double err_rt50 = check_rt(m50, z50);
        const double err_rt1050 = check_rt(m1050, z1050);

        std::printf("  d=%+.2f | p=%zu m_inv=%zu (esp=%zu) top=%d (esp=%d) | dif_tras(mu+1000)=%.3e | rt50=%.3e rt1050=%.3e\n",
                    d, m50.p, m50.m_inv_arfima, m_inv_esps[idx], int(m50.burn_in_topado), int(top_esps[idx]),
                    max_dif_tras, err_rt50, err_rt1050);

        check(m50.valido && m1050.valido, "ARFIMA ajustado valido");
        check(m50.m_inv_arfima == m_inv_esps[idx] && m50.n_hist_req == m50.p + m50.m_inv_arfima, "Horizonte Hosking exacto");
        check(m50.burn_in_topado == top_esps[idx], "Bandera burn_in_topado coincide");
        check(std::abs((m1050.mu_serie - m50.mu_serie) - 1000.0) <= 1e-11, "mu_serie traslada exactamente +1000");
        check(max_dif_tras <= 1e-11, "Invariancia traslacional mu -> mu+1000 en innovaciones()");
        check(err_rt50 <= 1e-12 && err_rt1050 <= 1e-12, "Ida y vuelta innovaciones -> simular <= 1e-12");
    }

    std::printf("\nFALLAS T2.4: %zu\n", fallas);
    return fallas == 0 ? 0 : 1;
}
