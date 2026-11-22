// ============================================================================
// test_hosking_unitario.cpp -- Verificacion unitaria de T2.3 (funciones puras
// de horizonte causal de Hosking y Gauss en modelo_lineal_nulo.hpp).
// ============================================================================

#include "ark/lilyth/modelo_lineal_nulo.hpp"

#include <cmath>
#include <cstddef>
#include <cstdio>

namespace {

using namespace ark::lilyth::detalle_lineal;
using ark::metrologia::est::norm_ppf;

}  // namespace

int main() {
    std::size_t fallas = 0;
    auto check = [&](bool cond, const char* msg) {
        if (!cond) {
            ++fallas;
            std::printf("  [FALLA] %s\n", msg);
        }
    };

    const double z = norm_ppf(1.0 - 0.5 * 0.05);
    const double eps2 = (z * z) / 4096.0;

    std::printf("================================================================================\n");
    std::printf(" TEST UNITARIO T2.3: HORIZONTE DE HOSKING (mse_pi, M_pi_exact, M_sim_gauss)\n");
    std::printf("================================================================================\n");
    std::printf("n=4096, alfa=0.05, z=%.12f, eps_Q^2=%.8f\n\n", z, eps2);

    struct Caso {
        double d;
        std::size_t m_inv_esp;
        std::size_t m_sim_esp;
        bool sim_top_esp;
        double mse64_esp;
        double mse256_esp;
    };
    const Caso casos[8] = {
        {-0.40, 338,    8, false, 0.004442, 0.001218},
        {-0.35, 211,    8, false, 0.002925, 0.000776},
        {-0.25,  82,    9, false, 0.001192, 0.000305},
        {-0.10,  11,    5, false, 0.000160, 0.000040},
        {+0.10,  11,   28, false, 0.000161, 0.000040},
        {+0.25,  85, 4096,  true, 0.001236, 0.000310},
        {+0.35, 233, 4096,  true, 0.003399, 0.000853},
        {+0.40, 418, 4096,  true, 0.006094, 0.001529},
    };

    for (const auto& c : casos) {
        const std::size_t m_inv = horizonte_inv_arfima_hosking(c.d, eps2, 4096);
        bool top = false;
        const std::size_t m_sim = horizonte_sim_gauss(-c.d, eps2, 4096, top);
        const double mse64 = mse_pi_exacto_hosking(c.d, 64);
        const double mse256 = mse_pi_exacto_hosking(c.d, 256);
        const double mse_m = mse_pi_exacto_hosking(c.d, m_inv);
        const double mse_mm1 = m_inv > 1 ? mse_pi_exacto_hosking(c.d, m_inv - 1) : 1.0;

        // Verificar monotonicidad estricta decreciente y fuerza bruta en M = 1..450
        bool mono_ok = true;
        std::size_t m_bruto = 0;
        double prev = mse_pi_exacto_hosking(c.d, 1);
        if (prev <= eps2) m_bruto = 1;
        for (std::size_t M = 2; M <= 450; ++M) {
            const double cur = mse_pi_exacto_hosking(c.d, M);
            if (!(cur <= prev + 1e-15)) mono_ok = false;
            if (m_bruto == 0 && cur <= eps2) m_bruto = M;
            prev = cur;
        }

        std::printf("  d=%+5.2f | M_inv=%3zu (esp=%3zu) | M_sim(-d)=%4zu top=%d (esp=%4zu/%d) | MSE(64)=%.6f MSE(256)=%.6f | mono=%d\n",
                    c.d, m_inv, c.m_inv_esp, m_sim, int(top), c.m_sim_esp, int(c.sim_top_esp),
                    mse64, mse256, int(mono_ok));

        check(m_inv == c.m_inv_esp, "M_inv coincide exactamente con tabla de verificacion");
        check(m_inv == m_bruto, "Galope+binaria coincide exactamente con busqueda por fuerza bruta");
        check(mse_m <= eps2 && mse_mm1 > eps2, "M_inv es el minimo estricto bajo eps_Q^2");
        check(m_sim == c.m_sim_esp && top == c.sim_top_esp, "M_sim(-d) y bandera topado coinciden");
        check(std::abs(mse64 - c.mse64_esp) <= 1e-6, "MSE(64) coincide con Python");
        check(std::abs(mse256 - c.mse256_esp) <= 1e-6, "MSE(256) coincide con Python");
        check(mono_ok, "MSE_pi(d, M) monotona no creciente en M=1..450");
    }

    // Verificacion sin tope en d = +0.25 -> M_sim(-0.25) = 18891
    bool top25 = true;
    const std::size_t m_sim_25k = horizonte_sim_gauss(-0.25, eps2, 25000, top25);
    std::printf("\n  d=+0.25 sin tope (cap=25000): M_sim(-0.25) = %zu, topado = %d (esp = 18891/0)\n",
                m_sim_25k, int(top25));
    check(m_sim_25k == 18891 && !top25, "M_sim(-0.25) == 18891 exacto");

    // Casos borde
    bool top_edge = true;
    check(horizonte_inv_arfima_hosking(0.0, eps2, 4096) == 0, "d == 0 devuelve 0");
    check(horizonte_inv_arfima_hosking(0.5, eps2, 4096) == 0, "d == 0.5 devuelve 0");
    check(horizonte_inv_arfima_hosking(-0.5, eps2, 4096) == 0, "d == -0.5 devuelve 0");
    check(horizonte_inv_arfima_hosking(0.25, eps2, 0) == 0, "m_tope == 0 devuelve 0");
    check(horizonte_inv_arfima_hosking(0.25, 0.0, 4096) == 0, "eps_var == 0 devuelve 0");
    check(mse_pi_exacto_hosking(0.0, 10) == 0.0, "mse(0, 10) == 0.0");
    check(std::isnan(mse_pi_exacto_hosking(0.5, 10)), "mse(0.5, 10) == NaN");
    check(horizonte_sim_gauss(0.0, eps2, 4096, top_edge) == 0 && !top_edge, "sim_gauss(0.0) == 0");

    std::printf("\nFALLAS T2.3: %zu\n", fallas);
    return fallas == 0 ? 0 : 1;
}
