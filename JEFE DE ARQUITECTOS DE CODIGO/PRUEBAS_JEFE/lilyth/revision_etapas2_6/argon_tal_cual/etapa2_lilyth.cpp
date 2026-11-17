// etapa2_lilyth.cpp -- Controles de la Etapa 2 del Operador de Lilyth: `modelo_lineal_nulo.hpp` (Argon, 2026-11-06).
// Construido sobre la plantilla de `etapa1_lilyth.cpp` de Code y `PLANO_ETAPA2_MODELO_LINEAL_NULO.md`.
//
// Controles:
//   L1   Levinson = Yule-Walker resuelto directamente (eliminacion de Gauss con pivoteo sobre la matriz de Toeplitz), a 1e-12
//   L2   AR(2) conocido (phi 0.6, -0.3): phi estimados dentro de 2 ee; p = 2; burn_in derivado exacto y estacionariedad en t = 0
//   L2b  AR(3) con kappa_2 = 0 (phi = 0.5, -0.2, 0.4; kappa = 0.5, 0.0, 0.4): Quenouille elige el ULTIMO kappa significativo (p = 3)
//   L3   ruido iid: p = 0, burn_in = 0
//   L4   ida y vuelta: `simular` con innovaciones fijas -> `innovaciones` sobre la serie simulada -> las mismas innovaciones a 1e-12
//   L5   el filtro inverso reconstruye las innovaciones de un AR(1) cuyo perfil empieza a mitad de serie (con su historia)
//   L6   ARFIMA(0, d, 0) con d = 0.3: `simular` -> diferenciar con d -> ACF compatible con ruido blanco
//   L7   ARFIMA con |d| >= 1/2, sin d, o con raiz unitaria: valido = false con su motivo
//   L8   FASE_ALEATORIZADA_FFT -> cayo_a_generativo = true
//   L9   residuos que no cruzan los cortes entre segmentos (= referencia independiente tramo a tramo)
//   L10  cero reservas de memoria dinamica en `innovaciones` y `simular`
//   P8   contra el `ObservadorSaltos` real y `ObservadorLilyth`: la inyeccion en una serie fija restituye en tau = 1 en FLUJO;
//        la del modelo generativo se propaga y restituye con la memoria del AR (> 1 paso)

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <random>
#include <vector>

#include "ark/lilyth/modelo_lineal_nulo.hpp"

static std::size_t g_reservas = 0;
static bool g_contar = false;
void* operator new(std::size_t n) { if (g_contar) ++g_reservas; if (void* p = std::malloc(n ? n : 1)) return p; throw std::bad_alloc{}; }
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

using namespace ark::metrologia;
using namespace ark::lilyth;

namespace {

int fallos = 0;
void chk(bool ok, const char* id, const char* txt) {
    std::printf("[%s] %-4s %s\n", ok ? "OK   " : "FALLA", id, txt);
    if (!ok) ++fallos;
}

// Eliminacion de Gauss con pivoteo parcial para resolver R * phi = r (matriz de Toeplitz simetrica R_{i,j} = rho_{|i-j|}).
std::vector<double> gauss_toeplitz_referencia(const std::vector<double>& rho, std::size_t m) {
    if (m == 0) return {};
    std::vector<std::vector<double>> A(m, std::vector<double>(m + 1, 0.0));
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            const std::size_t d = (i >= j) ? (i - j) : (j - i);
            A[i][j] = rho[d];
        }
        A[i][m] = rho[i + 1];
    }
    for (std::size_t col = 0; col < m; ++col) {
        std::size_t piv = col;
        for (std::size_t fila = col + 1; fila < m; ++fila) {
            if (std::abs(A[fila][col]) > std::abs(A[piv][col])) piv = fila;
        }
        std::swap(A[col], A[piv]);
        const double diag = A[col][col];
        for (std::size_t j = col; j <= m; ++j) A[col][j] /= diag;
        for (std::size_t fila = 0; fila < m; ++fila) {
            if (fila == col) continue;
            const double f = A[fila][col];
            for (std::size_t j = col; j <= m; ++j) A[fila][j] -= f * A[col][j];
        }
    }
    std::vector<double> sol(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) sol[i] = A[i][m];
    return sol;
}

}  // namespace

int main() {
    // -------------------------------------------------------------------------
    // L1: Levinson = Yule-Walker resuelto por eliminacion de Gauss sobre la Toeplitz, a 1e-12
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(20261106);
        std::normal_distribution<double> nd(0.0, 1.0);
        const std::size_t N = 8192;
        std::vector<double> s(N, 0.0);
        const std::array<double, 4> phi_v{0.55, -0.25, 0.20, -0.15};
        for (std::size_t t = 4; t < N; ++t) {
            s[t] = 2.5 + phi_v[0] * (s[t - 1] - 2.5) + phi_v[1] * (s[t - 2] - 2.5) +
                   phi_v[2] * (s[t - 3] - 2.5) + phi_v[3] * (s[t - 4] - 2.5) + nd(rng);
        }
        std::vector<std::vector<double>> segs{s};
        ConfiguracionModeloLineal cfg{};
        cfg.alfa = ALFA_PROYECTO;
        const auto m = ajustar(segs, cfg);

        const auto mo = est::momentos(segs);
        std::vector<double> rho(m.K + 1, 1.0);
        for (std::size_t k = 1; k <= m.K; ++k) {
            std::size_t pares = 0;
            rho[k] = est::acf_agrupada(segs, mo.mu, mo.sd, k, pares);
        }

        bool ok = m.valido && m.p >= 4;
        const auto phi_gauss_p = gauss_toeplitz_referencia(rho, m.p);
        double max_err_phi = 0.0, max_err_kappa = 0.0;
        for (std::size_t i = 0; i < m.p; ++i) {
            max_err_phi = std::max(max_err_phi, std::abs(m.phi[i] - phi_gauss_p[i]));
        }
        double v_ref = 1.0;
        for (std::size_t k = 1; k <= m.K; ++k) {
            const auto sol_k = gauss_toeplitz_referencia(rho, k);
            const double kk_ref = sol_k.back();
            max_err_kappa = std::max(max_err_kappa, std::abs(m.kappa[k - 1] - kk_ref));
            if (k <= m.p) v_ref *= (1.0 - kk_ref * kk_ref);
        }
        const double err_sigma2 = std::abs(m.sigma2_innov - mo.sd * mo.sd * v_ref);
        ok = ok && (max_err_phi <= 1e-12) && (max_err_kappa <= 1e-12) && (err_sigma2 <= 1e-12);
        std::printf("      L1: K=%zu p=%zu max_err_phi=%.2e max_err_kappa=%.2e err_sigma2=%.2e\n",
                    m.K, m.p, max_err_phi, max_err_kappa, err_sigma2);
        chk(ok, "L1", "Levinson-Durbin = eliminacion de Gauss directa sobre la Toeplitz (phi, kappa y sigma2_innov a <= 1e-12)");
    }

    // -------------------------------------------------------------------------
    // L2 y L2b: AR(2) conocido (phi 0.6, -0.3) -> phi dentro de 2 ee, p = 2, burn_in efectivo en t = 0;
    // y AR(3) con kappa_2 = 0 (phi = 0.5, -0.2, 0.4) -> Quenouille elige el ultimo significativo (p = 3)
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(1002);
        std::normal_distribution<double> nd(0.0, 1.0);
        const std::size_t N = 16384;
        std::vector<double> x(N, 0.0);
        for (std::size_t t = 2; t < N; ++t) {
            x[t] = 0.6 * x[t - 1] - 0.3 * x[t - 2] + nd(rng);
        }
        ConfiguracionModeloLineal cfg{};
        cfg.alfa = ALFA_PROYECTO;
        const auto m = ajustar({x}, cfg);
        const double ee = std::sqrt((1.0 - 0.3 * 0.3) / static_cast<double>(N));
        const std::size_t burn_esp = static_cast<std::size_t>(
            std::ceil(std::log(std::numeric_limits<double>::epsilon()) / std::log(m.radio_max)));

        // Verificamos que el burn-in realmente corre en `simular`: en un AR(1) con phi = 0.92, sin burn-in la varianza en t = 0
        // seria sigma2_innov ~ 0.15 sd^2; con burn-in es ~ 1.0 sd^2.
        std::vector<double> x_lento(N, 0.0);
        for (std::size_t t = 1; t < N; ++t) x_lento[t] = 0.92 * x_lento[t - 1] + nd(rng);
        const auto m_lento = ajustar({x_lento}, cfg);
        double s2_t0 = 0.0;
        const std::size_t R_burn = 256;
        std::array<double, 4> buf_t0{};
        for (std::size_t r = 0; r < R_burn; ++r) {
            simular(m_lento, rng, nullptr, 0, 0, buf_t0.data(), buf_t0.size());
            const double d0 = buf_t0[0] - m_lento.mu;
            s2_t0 += d0 * d0;
        }
        s2_t0 /= static_cast<double>(R_burn);
        const double razon_var_t0 = s2_t0 / (m_lento.sd * m_lento.sd);

        const bool ok = m.valido && m.p == 2 &&
                        std::abs(m.phi[0] - 0.6) <= 2.0 * ee &&
                        std::abs(m.phi[1] - (-0.3)) <= 2.0 * ee &&
                        m.burn_in == burn_esp && m.burn_in > 0 &&
                        m_lento.valido && m_lento.burn_in > 0 && razon_var_t0 > 0.6;
        std::printf("      L2: p=%zu phi1=%.4f (0.6, 2ee=%.4f) phi2=%.4f (-0.3) radio=%.4f burn_in=%zu var_t0/sd2=%.3f\n",
                    m.p, m.p > 0 ? m.phi[0] : NaN, 2.0 * ee, m.p > 1 ? m.phi[1] : NaN, m.radio_max, m.burn_in, razon_var_t0);
        chk(ok, "L2", "AR(2) conocido (phi 0.6, -0.3): phi dentro de 2 ee, p = 2 y burn_in calienta la varianza en t = 0");

        // L2b: kappa_1 = 0.5, kappa_2 = 0.0 (no significativo), kappa_3 = 0.4 -> phi = (0.5, -0.2, 0.4).
        std::vector<double> x3(N, 0.0);
        for (std::size_t t = 3; t < N; ++t) {
            x3[t] = 0.5 * x3[t - 1] - 0.2 * x3[t - 2] + 0.4 * x3[t - 3] + nd(rng);
        }
        const auto m3 = ajustar({x3}, cfg);
        const double umbral_q = est::norm_ppf(1.0 - 0.5 * cfg.alfa) / std::sqrt(static_cast<double>(N));
        const bool ok3 = m3.valido && m3.K >= 3 &&
                         std::abs(m3.kappa[1]) <= umbral_q &&
                         std::abs(m3.kappa[2]) > umbral_q &&
                         m3.p >= 3 && std::abs(m3.phi[2] - 0.4) < 0.05;
        std::printf("      L2b: K=%zu p=%zu kappa1=%.3f kappa2=%.3f (<= %.3f) kappa3=%.3f\n",
                    m3.K, m3.p, m3.kappa[0], m3.kappa[1], umbral_q, m3.kappa[2]);
        chk(ok3, "L2b", "Quenouille elige el ULTIMO kappa significativo (p >= 3 aunque kappa_2 no sea significativo)");
    }

    // -------------------------------------------------------------------------
    // L3: ruido iid -> p = 0
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(2026);
        std::normal_distribution<double> nd(5.0, 2.0);
        std::vector<double> x(8192);
        for (double& v : x) v = nd(rng);
        ConfiguracionModeloLineal cfg{};
        cfg.alfa = ALFA_PROYECTO;
        const auto m = ajustar({x}, cfg);
        chk(m.valido && m.p == 0 && m.phi.empty() && m.burn_in == 0 && m.radio_max == 0.0,
            "L3", "ruido iid: p = 0, phi vacio, radio_max = 0 y burn_in = 0");
    }

    // -------------------------------------------------------------------------
    // L4: ida y vuelta: `simular` con innovaciones fijas -> `innovaciones` sobre la serie simulada -> las mismas a 1e-12
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(3004);
        std::normal_distribution<double> nd(0.0, 1.0);
        std::vector<double> cal(4096, 0.0);
        for (std::size_t t = 2; t < cal.size(); ++t) cal[t] = 1.5 + 0.65 * (cal[t - 1] - 1.5) - 0.25 * (cal[t - 2] - 1.5) + nd(rng);
        ConfiguracionModeloLineal cfg{};
        cfg.alfa = ALFA_PROYECTO;
        const auto m = ajustar({cal}, cfg);

        const std::size_t n_sim = 256;
        std::vector<double> e_fijas(n_sim), x_sim(n_sim), e_rec(n_sim);
        for (double& e : e_fijas) e = nd(rng);

        simular(m, rng, e_fijas.data(), 0, n_sim, x_sim.data(), n_sim);
        innovaciones(m, nullptr, 0, x_sim.data(), n_sim, e_rec.data());

        double max_dif = 0.0;
        for (std::size_t i = 0; i < n_sim; ++i) max_dif = std::max(max_dif, std::abs(e_fijas[i] - e_rec[i]));

        ConfiguracionModeloLineal cfg_arf = cfg;
        cfg_arf.modo = ModoNuloLineal::ARFIMA_WHITTLE;
        cfg_arf.d_externa = 0.25;
        const auto marf = ajustar({cal}, cfg_arf);
        std::vector<double> x_arf(n_sim), e_rec_arf(n_sim);
        simular(marf, rng, e_fijas.data(), 0, n_sim, x_arf.data(), n_sim);
        innovaciones(marf, nullptr, 0, x_arf.data(), n_sim, e_rec_arf.data());
        for (std::size_t i = 0; i < n_sim; ++i) max_dif = std::max(max_dif, std::abs(e_fijas[i] - e_rec_arf[i]));

        std::printf("      L4: max_dif (AR y ARFIMA) = %.2e\n", max_dif);
        chk(max_dif <= 1e-12, "L4", "ida y vuelta: simular con innovaciones fijas -> innovaciones = las mismas a <= 1e-12");
    }

    // -------------------------------------------------------------------------
    // L5: el filtro inverso reconstruye las innovaciones de un AR(1) cuyo perfil empieza a mitad de serie (con su historia)
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(5005);
        std::normal_distribution<double> nd(0.0, 1.0);
        const std::size_t N = 4096;
        std::vector<double> cal(N, 0.0);
        for (std::size_t t = 1; t < N; ++t) cal[t] = -3.0 + 0.8 * (cal[t - 1] - (-3.0)) + nd(rng);
        ConfiguracionModeloLineal cfg{};
        cfg.alfa = ALFA_PROYECTO;
        const auto m = ajustar({cal}, cfg);

        const std::size_t n_sim = 400, t_iny = 200, l_perf = 40;
        std::vector<double> e_perf(l_perf), x_sim(n_sim), e_rec(l_perf);
        for (std::size_t j = 0; j < l_perf; ++j) e_perf[j] = 3.0 + 0.2 * static_cast<double>(j);

        simular(m, rng, e_perf.data(), t_iny, l_perf, x_sim.data(), n_sim);
        innovaciones(m, x_sim.data(), t_iny, x_sim.data() + t_iny, l_perf, e_rec.data());

        double max_dif = 0.0;
        for (std::size_t j = 0; j < l_perf; ++j) max_dif = std::max(max_dif, std::abs(e_perf[j] - e_rec[j]));
        std::printf("      L5: p=%zu max_dif con historia previa = %.2e\n", m.p, max_dif);
        chk(m.valido && m.p >= 1 && max_dif <= 1e-12,
            "L5", "filtro inverso reconstruye las innovaciones de un AR(1) a mitad de serie usando su historia previa (<= 1e-12)");
    }

    // -------------------------------------------------------------------------
    // L6: ARFIMA(0, d, 0) con d = 0.3: simular -> diferenciar con d -> ACF compatible con ruido blanco
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(6008);
        std::normal_distribution<double> nd(0.0, 1.0);
        const std::size_t N = 4096;
        const double d_verdad = 0.3;
        std::vector<double> u_cal(N), x_cal(N, 0.0), pi_neg(N), pi_pos(N);
        for (double& v : u_cal) v = nd(rng);
        elw::coeficientes(-d_verdad, N, pi_neg.data());
        elw::coeficientes(d_verdad, N, pi_pos.data());
        for (std::size_t t = 0; t < N; ++t) {
            double s = 0.0;
            for (std::size_t k = 0; k <= t; ++k) s += pi_neg[k] * u_cal[t - k];
            x_cal[t] = s;
        }

        ConfiguracionModeloLineal cfg{};
        cfg.modo = ModoNuloLineal::ARFIMA_WHITTLE;
        cfg.alfa = ALFA_PROYECTO;
        cfg.d_externa = d_verdad;
        cfg.innovaciones = InnovacionesNulo::GAUSS;
        const auto m = ajustar({x_cal}, cfg);

        std::vector<double> x_sim(N, 0.0), u_dif(N, 0.0);
        simular(m, rng, nullptr, 0, 0, x_sim.data(), N);
        for (std::size_t t = 0; t < N; ++t) {
            double s = 0.0;
            for (std::size_t k = 0; k <= t; ++k) s += pi_pos[k] * x_sim[t - k];
            u_dif[t] = s;
        }
        const auto mo_x = est::momentos({x_sim});
        std::size_t pares = 0;
        const double acf1_x = est::acf_agrupada({x_sim}, mo_x.mu, mo_x.sd, 1, pares);
        const auto mo_u = est::momentos({u_dif});
        const double acf1_u = est::acf_agrupada({u_dif}, mo_u.mu, mo_u.sd, 1, pares);
        const double cota_blanca = 3.0 / std::sqrt(static_cast<double>(N));
        const bool ok = m.valido && m.p == 0 && m.burn_in == N && acf1_x > 0.25 && std::abs(acf1_u) <= cota_blanca;
        std::printf("      L6: p=%zu burn_in=%zu acf1(x_sim)=%.4f acf1(Delta^d x_sim)=%.4f (cota %.4f)\n",
                    m.p, m.burn_in, acf1_x, acf1_u, cota_blanca);
        chk(ok, "L6", "ARFIMA(0, 0.3, 0): simular genera memoria larga (acf1 > 0.25) y al diferenciar con d = 0.3 vuelve a ruido blanco");
    }

    // -------------------------------------------------------------------------
    // L7: ARFIMA con |d| >= 1/2, sin d, o con raiz unitaria: valido = false con su motivo
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(7007);
        std::normal_distribution<double> nd(0.0, 1.0);
        std::vector<double> x_cal(2048), paseo(2048, 0.0);
        for (std::size_t t = 0; t < x_cal.size(); ++t) {
            x_cal[t] = nd(rng);
            paseo[t] = (t > 0 ? paseo[t - 1] : 0.0) + x_cal[t];
        }
        ConfiguracionModeloLineal cfg{};
        cfg.alfa = ALFA_PROYECTO;
        cfg.modo = ModoNuloLineal::ARFIMA_WHITTLE;
        cfg.d_externa = NaN;
        const auto m_sin_d = ajustar({x_cal}, cfg);

        cfg.d_externa = 0.5;
        const auto m_d_alto = ajustar({x_cal}, cfg);

        cfg.d_externa = -0.6;
        const auto m_d_bajo = ajustar({x_cal}, cfg);

        cfg.modo = ModoNuloLineal::AR_LEVINSON_DURBIN;
        const auto m_rw = ajustar({paseo}, cfg);

        const bool ok = !m_sin_d.valido && !m_d_alto.valido && !m_d_bajo.valido && !m_rw.valido &&
                        std::strcmp(m_sin_d.motivo, "sin ajustar") != 0 &&
                        std::strcmp(m_d_alto.motivo, "sin ajustar") != 0 &&
                        std::strcmp(m_rw.motivo, "sin ajustar") != 0;
        std::printf("      L7 motivos: sin_d='%s' | d=0.5='%s' | rw='%s'\n", m_sin_d.motivo, m_d_alto.motivo, m_rw.motivo);
        chk(ok, "L7", "ARFIMA sin d o con |d| >= 0.5, y paseo aleatorio (raiz unitaria) dan valido = false con su motivo");
    }

    // -------------------------------------------------------------------------
    // L8: FASE_ALEATORIZADA_FFT -> cayo_a_generativo = true
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(8008);
        std::normal_distribution<double> nd(0.0, 1.0);
        std::vector<double> x(1024);
        for (double& v : x) v = nd(rng);
        ConfiguracionModeloLineal cfg{};
        cfg.alfa = ALFA_PROYECTO;
        cfg.modo = ModoNuloLineal::FASE_ALEATORIZADA_FFT;
        const auto m = ajustar({x}, cfg);
        chk(m.valido && m.cayo_a_generativo && m.modo_usado == ModoNuloLineal::AR_LEVINSON_DURBIN,
            "L8", "FASE_ALEATORIZADA_FFT marca cayo_a_generativo = true y usa AR_LEVINSON_DURBIN (E1.19)");
    }

    // -------------------------------------------------------------------------
    // L9: residuos que no cruzan los cortes entre segmentos
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(9009);
        std::normal_distribution<double> nd(0.0, 1.0);
        std::vector<double> seg1(1000, 0.0), seg2(1000, 0.0);
        for (std::size_t t = 1; t < seg1.size(); ++t) seg1[t] = 0.7 * seg1[t - 1] + nd(rng);
        for (std::size_t t = 1; t < seg2.size(); ++t) seg2[t] = 0.7 * seg2[t - 1] + nd(rng);
        // Valores extremos opuestos justo en los bordes del corte para cazar cualquier residuo que cruce el corte:
        seg1.back() = 25.0;
        seg2.front() = -25.0;

        std::vector<std::vector<double>> segs{seg1, seg2};
        ConfiguracionModeloLineal cfg{};
        cfg.alfa = ALFA_PROYECTO;
        const auto m = ajustar(segs, cfg);

        std::vector<double> ref_res;
        for (const auto& g : segs) {
            for (std::size_t t = m.p; t < g.size(); ++t) {
                double e = g[t] - m.mu;
                for (std::size_t i = 1; i <= m.p; ++i) e -= m.phi[i - 1] * (g[t - i] - m.mu);
                ref_res.push_back(e);
            }
        }
        bool ok = m.valido && m.p >= 1 && m.residuos.size() == (seg1.size() - m.p) + (seg2.size() - m.p) &&
                  m.residuos.size() == ref_res.size();
        for (std::size_t i = 0; ok && i < ref_res.size(); ++i) {
            if (std::abs(m.residuos[i] - ref_res[i]) > 1e-12) ok = false;
        }
        chk(ok, "L9", "los residuos se calculan con t >= p dentro de cada segmento y nunca cruzan el corte entre segmentos");
    }

    // -------------------------------------------------------------------------
    // L10: cero reservas en `innovaciones` y `simular`
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(1010);
        std::normal_distribution<double> nd(0.0, 1.0);
        std::vector<double> cal(2048, 0.0);
        for (std::size_t t = 1; t < cal.size(); ++t) cal[t] = 0.7 * cal[t - 1] + nd(rng);
        ConfiguracionModeloLineal cfg{};
        cfg.alfa = ALFA_PROYECTO;
        const auto m_ar = ajustar({cal}, cfg);

        cfg.modo = ModoNuloLineal::ARFIMA_WHITTLE;
        cfg.d_externa = 0.25;
        const auto m_arf = ajustar({cal}, cfg);

        std::array<double, 512> x_out{}, e_out{};
        g_reservas = 0;
        g_contar = true;
        simular(m_ar, rng, nullptr, 0, 0, x_out.data(), x_out.size());
        innovaciones(m_ar, x_out.data(), 100, x_out.data() + 100, 50, e_out.data());
        simular(m_arf, rng, nullptr, 0, 0, x_out.data(), x_out.size());
        innovaciones(m_arf, x_out.data(), 100, x_out.data() + 100, 50, e_out.data());
        g_contar = false;
        chk(g_reservas == 0, "L10", "cero reservas de memoria dinamica en innovaciones y simular (AR y ARFIMA)");
    }

    // -------------------------------------------------------------------------
    // P8: contra el ObservadorSaltos real y ObservadorLilyth: la inyeccion en una serie fija restituye en tau = 1
    // en FLUJO; la del modelo generativo se propaga y restituye con la memoria del AR (> 1 paso)
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 512, 4096>;
        std::mt19937_64 rng(8888);
        std::normal_distribution<double> nd(0.0, 1.0);

        const std::size_t N_CAL = 1 << 14;
        const double phi_real = 0.85;
        std::vector<double> xc(N_CAL, 0.0);
        for (std::size_t t = 1; t < N_CAL; ++t) xc[t] = phi_real * xc[t - 1] + nd(rng);

        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 5000.0;
        auto obs_base = std::make_unique<Obs>();
        obs_base->calibrar(cd, xc);

        ConfiguracionModeloLineal cfg_m{};
        cfg_m.alfa = ALFA_PROYECTO;
        cfg_m.innovaciones = InnovacionesNulo::RESIDUOS_REMUESTREADOS;
        const auto mod = ajustar({xc}, cfg_m);

        const std::size_t N_PRUEBA = 1200;
        std::vector<double> x_real(N_PRUEBA, 0.0);
        for (std::size_t t = 1; t < N_PRUEBA; ++t) {
            const double shock = (t == 500) ? 14.0 : 0.0;
            x_real[t] = phi_real * x_real[t - 1] + shock + 0.25 * nd(rng);
        }

        auto obs_real = std::make_unique<Obs>(*obs_base);
        obs_real->reiniciar_flujo();
        std::size_t t_a = 0, l_perf = 0, t_cam = 0;
        for (std::size_t t = 0; t < N_PRUEBA; ++t) {
            const auto s = obs_real->paso(x_real[t]);
            if (s.nivel[0].alarma.hay && t >= 500) {
                t_a = t;
                l_perf = s.nivel[0].alarma.l;
                t_cam = static_cast<std::size_t>(s.nivel[0].alarma.t_cambio);
                break;
            }
        }
        std::vector<double> e_perfil(l_perf, 0.0);
        innovaciones(mod, x_real.data(), t_cam, x_real.data() + t_cam, l_perf, e_perfil.data());

        const std::size_t t_iny = 500;
        std::vector<double> x_fija(N_PRUEBA, 0.0);
        simular(mod, rng, nullptr, 0, 0, x_fija.data(), N_PRUEBA);
        for (std::size_t j = 0; j < l_perf; ++j) {
            x_fija[t_iny + j] = x_real[t_cam + j];
        }

        auto medir_tau_primer_evento = [&](const std::vector<double>& serie) -> double {
            auto d = std::make_unique<Obs>(*obs_base);
            d->reiniciar_flujo();
            auto l = std::make_unique<Lil>();
            ConfiguracionLilyth cl{};
            cl.representacion = ModoRepresentacionLilyth::CRUDA;
            cl.t_max_efectivo = 100;
            l->configurar(cl, cd, *d, xc);
            for (std::size_t t = 0; t < serie.size(); ++t) {
                const auto sd = d->paso(serie[t]);
                const auto sl = l->paso(sd, serie[t], true);
                for (std::size_t i = 0; i < sl.ranuras_cerradas; ++i) {
                    const auto& r = l->cerrada(i);
                    if (r.t_apertura >= static_cast<std::int64_t>(t_iny)) {
                        return r.cruda.tau_rest_consolidado;
                    }
                }
            }
            return NaN;
        };

        const double tau_fija = medir_tau_primer_evento(x_fija);

        std::vector<double> x_gen(N_PRUEBA, 0.0);
        simular(mod, rng, e_perfil.data(), t_iny, l_perf, x_gen.data(), N_PRUEBA);
        const double tau_gen = medir_tau_primer_evento(x_gen);

        std::printf("      P8: t_a=%zu l*=%zu tau_serie_fija=%.1f (esperado 1.0), tau_modelo_generativo=%.1f (esperado >= 3.0)\n",
                    t_a, l_perf, tau_fija, tau_gen);
        chk(tau_fija == 1.0 && std::isfinite(tau_gen) && tau_gen >= 3.0,
            "P8", "contra ObservadorSaltos real: inyeccion en serie fija restituye en tau = 1; el modelo generativo se propaga (tau >= 3)");
    }

    std::printf("\n%s: %d fallas\n", fallos ? "FALLA" : "PASA", fallos);
    return fallos ? 1 : 0;
}