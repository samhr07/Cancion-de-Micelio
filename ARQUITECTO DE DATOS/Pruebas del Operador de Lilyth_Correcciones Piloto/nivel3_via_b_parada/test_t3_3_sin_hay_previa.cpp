// ============================================================================
// test_t3_3_sin_hay_previa.cpp -- Pruebas unitarias aisladas de T3.3:
//   C1. Inyeccion simple con aperturas espontaneas en ancla (ARL0=150):
//       usar_enmienda_nivel3=true -> solapes==0 y rescata det_n3 > det_base;
//       usar_enmienda_nivel3=false -> solapes == con_ajena_previa (Regla 6).
//   C2. Inyecciones multiples i != j bajo ARL0=150:
//       en borde t_max_ef -> solapes==2; separadas -> solapes==0 pese a hay_previa>0.
//   C3. E2E completar_nulos() via (b) con forzar_via_b=true (n3=true vs false):
//       solapes_b_n3 == 0, con_ajena_previa_b > 0 en ReferenciaNuloEvento vs solapes_b_base > 0.
// Semilla ciega exclusiva: SEMILLA_CIEGA_NIVEL3 = 20261124ULL.
// ============================================================================

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;

constexpr std::uint64_t SEMILLA_CIEGA_NIVEL3 = 20261124ULL;

std::vector<double> generar_ar1_unidad(std::size_t n, double phi, double mu, std::uint64_t sem) {
    const std::size_t burn = 512;
    const double sd_e = std::sqrt(std::max(1e-12, 1.0 - phi * phi));
    std::mt19937_64 rng(detalle_nulo::splitmix64(sem));
    std::normal_distribution<double> nd(0.0, sd_e);
    std::vector<double> out(n, 0.0);
    double u = 0.0;
    for (std::size_t t = 0; t < burn + n; ++t) {
        u = phi * u + nd(rng);
        if (t >= burn) out[t - burn] = mu + u;
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

    std::printf("================================================================================\n");
    std::printf(" TEST T3.3: RETIRADA DE hay_previa EN INYECCION SIMPLE Y SOLAPE SOLO ENTRE i!=j\n");
    std::printf("================================================================================\n\n");

    // -------------------------------------------------------------------------
    // C1: Inyeccion simple con aperturas espontaneas en ancla (n3=true vs false)
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 64>;

        const auto x_cal = generar_ar1_unidad(4096, 0.80, 0.0, SEMILLA_CIEGA_NIVEL3 + 3001ULL);
        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 150.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<Obs>();
        od->calibrar(cd, x_cal);

        ConfiguracionLilyth cl_n3{};
        cl_n3.representacion = ModoRepresentacionLilyth::CRUDA;
        cl_n3.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl_n3.b_max_nulo = 4;
        cl_n3.hilos_nulo = 1;
        cl_n3.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 3002ULL;
        cl_n3.usar_enmienda_nivel2 = true;
        cl_n3.usar_enmienda_nivel3 = true;

        ConfiguracionLilyth cl_base = cl_n3;
        cl_base.usar_enmienda_nivel3 = false;

        auto lil_n3 = std::make_unique<Lil>();
        auto lil_base = std::make_unique<Lil>();
        check(lil_n3->calibrar(*od, cl_n3, cd, x_cal), "C1 calibrar lil_n3");
        check(lil_base->calibrar(*od, cl_base, cd, x_cal), "C1 calibrar lil_base");

        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od->ficha_nivel(0).nulo.tau_int)));
        const std::vector<double> hist_prev(lil_n3->ficha_nulo_cruda().modelo_lineal.p + 8u, 0.0);
        const std::vector<double> perfil{4.0, 6.0, 7.5, 8.5};
        const std::array<std::int64_t, 1> t_arr{250};
        const std::size_t n_sim = 600;

        std::size_t prev_base = 0, sol_base = 0, det_base = 0;
        std::size_t prev_n3 = 0, sol_n3 = 0, det_n3 = 0;
        for (std::size_t k = 0; k < 80; ++k) {
            const std::uint64_t sem = SEMILLA_CIEGA_NIVEL3 + 3100ULL + k;
            const auto rb = lil_base->evaluar_inyecciones_realizacion(
                false, hist_prev, perfil, +1, w_ancla, guarda, t_arr, n_sim, sem, false, true);
            const auto r3 = lil_n3->evaluar_inyecciones_realizacion(
                false, hist_prev, perfil, +1, w_ancla, guarda, t_arr, n_sim, sem, false, true);
            prev_base += rb.con_ajena_previa;
            sol_base += rb.solapes;
            det_base += rb.detectadas;
            prev_n3 += r3.con_ajena_previa;
            sol_n3 += r3.solapes;
            det_n3 += r3.detectadas;
        }
        std::printf("  [C1 80 iny] n3=0(prev=%zu sol=%zu det=%zu) | n3=1(prev=%zu sol=%zu det=%zu)\n",
                    prev_base, sol_base, det_base, prev_n3, sol_n3, det_n3);
        check(prev_base > 0 && sol_base == prev_base, "C1 Regla 6: n3=0 mantiene solapes == con_ajena_previa");
        check(prev_n3 > 0 && sol_n3 == 0 && det_n3 > det_base,
              "C1 n3=1 retira veto hay_previa en inyeccion simple y rescata detecciones");
    }

    // -------------------------------------------------------------------------
    // C2: Inyecciones multiples i != j bajo ARL0=150
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 128>;

        const auto x_cal = generar_ar1_unidad(4096, 0.85, 0.0, SEMILLA_CIEGA_NIVEL3 + 3011ULL);
        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 150.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<Obs>();
        od->calibrar(cd, x_cal);

        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl.b_max_nulo = 8;
        cl.hilos_nulo = 1;
        cl.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 3012ULL;
        cl.usar_enmienda_nivel2 = true;
        cl.usar_enmienda_nivel3 = true;

        auto lil = std::make_unique<Lil>();
        check(lil->calibrar(*od, cl, cd, x_cal), "C2 calibrar");

        const std::size_t t_max_ef = lil->configuracion().t_max_efectivo;
        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od->ficha_nivel(0).nulo.tau_int)));
        const std::size_t l = 4;
        const std::size_t dist_borde = w_ancla + guarda + l + t_max_ef;
        const std::size_t dist_entre = w_ancla + guarda + l + t_max_ef + 16;
        check(t_max_ef + 16 < 128, "C2 t_max_ef + 16 < 128");

        const std::vector<double> hist_prev(lil->ficha_nulo_cruda().modelo_lineal.p + 8u, 0.0);
        const std::vector<double> perfil{4.0, 6.0, 7.5, 8.5};
        const std::array<std::int64_t, 2> t_borde{250, static_cast<std::int64_t>(250 + dist_borde)};
        const std::array<std::int64_t, 2> t_entre{250, static_cast<std::int64_t>(250 + dist_entre)};
        const std::size_t n_sim = 1200;

        const auto rb = lil->evaluar_inyecciones_realizacion(
            false, hist_prev, perfil, +1, w_ancla, guarda, t_borde, n_sim, SEMILLA_CIEGA_NIVEL3 + 3199ULL, false, true);
        std::size_t sol_tot = 0, prev_tot = 0, adm_tot = 0;
        for (std::size_t s = 0; s < 20; ++s) {
            const auto rs = lil->evaluar_inyecciones_realizacion(
                false, hist_prev, perfil, +1, w_ancla, guarda, t_entre, n_sim, SEMILLA_CIEGA_NIVEL3 + 3200ULL + s, false, true);
            sol_tot += rs.solapes;
            prev_tot += rs.con_ajena_previa;
            adm_tot += rs.eventos_admitidos.size();
        }
        std::printf("  [C2] t_max_ef=%zu | borde(sol=%zu adm=%zu) | 20 pares separados(sol=%zu prev=%zu adm=%zu)\n",
                    t_max_ef, rb.solapes, rb.eventos_admitidos.size(), sol_tot, prev_tot, adm_tot);
        check(rb.solapes == 2 && rb.eventos_admitidos.empty(),
              "C2 solape geometrico entre inyecciones i!=j sigue vetando en borde");
        check(sol_tot == 0 && prev_tot > 0 && adm_tot > 0,
              "C2 pares separados tienen solapes==0 pese a hay_previa>0");
    }

    // -------------------------------------------------------------------------
    // C3: E2E completar_nulos() via (b) con forzar_via_b=true (n3=true vs false)
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 64>;

        const auto x_cal = generar_ar1_unidad(4096, 0.80, 0.0, SEMILLA_CIEGA_NIVEL3 + 3021ULL);
        auto x_op = generar_ar1_unidad(4096, 0.80, 0.0, SEMILLA_CIEGA_NIVEL3 + 3022ULL);
        const auto mom = est::momentos(std::vector<std::vector<double>>{x_cal});
        for (std::size_t t = 1200; t <= 1202; ++t) x_op[t] += 9.0 * mom.sd;
        for (std::size_t t = 2400; t <= 2402; ++t) x_op[t] += 9.0 * mom.sd;

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 200.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<Obs>();
        od->calibrar(cd, x_cal);

        auto correr_e2e = [&](bool usar_n3, std::size_t& n_eval, std::size_t& sol_b, std::size_t& prev_b, std::size_t& post_b) {
            ConfiguracionLilyth cl{};
            cl.representacion = ModoRepresentacionLilyth::CRUDA;
            cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
            cl.b_max_nulo = 4;
            cl.hilos_nulo = 1;
            cl.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 3023ULL;
            cl.forzar_via_b = true;
            cl.usar_enmienda_nivel2 = true;
            cl.usar_enmienda_nivel3 = usar_n3;

            auto lil = std::make_unique<Lil>();
            check(lil->calibrar(*od, cl, cd, x_cal), "C3 calibrar");
            auto od_run = std::make_unique<Obs>(*od);
            od_run->reiniciar_flujo();
            lil->reiniciar_flujo_conservando_purga();
            for (double v : x_op) {
                const auto sd = od_run->paso(v, NaN, true, 0, 0);
                (void)lil->paso(sd, v, true, NaN, 0, 0);
            }
            lil->completar_nulos();
            n_eval = 0;
            sol_b = 0;
            prev_b = 0;
            post_b = 0;
            for (const auto& r : lil->cerradas_historico()) {
                const auto d = r.cruda.desenlace;
                if (d != DesenlaceRanura::RESTITUIDA && d != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
                const auto& ref = r.cruda.nulo_onsager;
                if (ref.evaluado) ++n_eval;
                sol_b += ref.solapes_b;
                prev_b += ref.con_ajena_previa_b;
                post_b += ref.con_reapertura_posterior_b;
            }
        };

        std::size_t ev3 = 0, sol3 = 0, prev3 = 0, post3 = 0;
        std::size_t ev0 = 0, sol0 = 0, prev0 = 0, post0 = 0;
        correr_e2e(true,  ev3, sol3, prev3, post3);
        correr_e2e(false, ev0, sol0, prev0, post0);

        std::printf("  [C3 E2E] n3=1(eval=%zu sol_b=%zu prev_b=%zu post_b=%zu) | n3=0(eval=%zu sol_b=%zu prev_b=%zu post_b=%zu)\n",
                    ev3, sol3, prev3, post3, ev0, sol0, prev0, post0);
        check(ev3 > 0 && sol3 == 0 && prev3 > 0, "C3 n3=1 evalua eventos con solapes_b==0 y telemetria prev_b>0");
        check(sol0 > 0, "C3 Regla 6: n3=0 mantiene solapes_b>0 bajo ARL0=200");
    }

    std::printf("\nFALLAS T3.3: %zu\n", fallas);
    return fallas == 0 ? 0 : 1;
}
