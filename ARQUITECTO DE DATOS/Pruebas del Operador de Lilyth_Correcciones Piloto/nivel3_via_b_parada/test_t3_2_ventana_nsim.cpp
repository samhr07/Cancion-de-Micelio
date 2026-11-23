// ============================================================================
// test_t3_2_ventana_nsim.cpp -- Pruebas unitarias aisladas de T3.2:
//   C1. dimensiones_simulacion_via_b() en FLUJO, NIVEL (+1) y w_ancla=0.
//   C2. S2_N3: ventana con t_max_ef real (< T_MAX=128) vs T_MAX fijo bajo
//       usar_enmienda_nivel3 true/false (y passthrough en evaluar_inyecciones_realizacion).
//   C3. Alineacion ARFIMA mu_serie + mu en ejecutar_inyecciones_en_realizacion:
//       invariancia exacta ante traslacion mu_z = 50 -> 1050.
//   C4. E2E ejecutar_via_b_evento() con n_sim derivado en FLUJO y NIVEL.
// Semilla ciega exclusiva: SEMILLA_CIEGA_NIVEL3 = 20261124ULL.
// ============================================================================

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <random>
#include <span>
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

std::vector<double> generar_arfima1(
    std::size_t n, double phi1, double d, double mu_serie, double sigma_e, std::uint64_t sem) {
    const std::size_t burn = 1024;
    const std::size_t tot = burn + n;
    std::vector<double> pi_inv(tot, 0.0);
    elw::coeficientes(-d, tot, pi_inv.data());
    std::vector<double> u(tot, 0.0);
    std::vector<double> x(n, 0.0);
    std::mt19937_64 rng(detalle_nulo::splitmix64(sem));
    std::normal_distribution<double> nd(0.0, sigma_e);
    double u_prev = 0.0;
    for (std::size_t s = 0; s < tot; ++s) {
        const double u_cur = phi1 * u_prev + nd(rng);
        u[s] = u_cur;
        u_prev = u_cur;
        if (s >= burn) {
            double suma = 0.0;
            for (std::size_t k = 0; k <= s; ++k) suma += pi_inv[k] * u[s - k];
            x[s - burn] = mu_serie + suma;
        }
    }
    return x;
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
    std::printf(" TEST T3.2: VENTANA DERIVADA t_max_ef, n_sim SIN LITERALES Y ARFIMA mu_serie\n");
    std::printf("================================================================================\n\n");

    // -------------------------------------------------------------------------
    // C1: dimensiones_simulacion_via_b() en FLUJO, NIVEL (+1) y w_ancla=0
    // -------------------------------------------------------------------------
    {
        const auto d_fl = dimensiones_simulacion_via_b(120, 14, 5, 22, Integracion::FLUJO);
        const auto d_nv = dimensiones_simulacion_via_b(120, 14, 5, 22, Integracion::NIVEL);
        const auto d_0  = dimensiones_simulacion_via_b(0, 8, 3, 10, Integracion::NIVEL);

        const bool ok_fl = d_fl.w_pre == 134 && d_fl.w_sorteo == 120 && d_fl.w_post == 27 &&
                           d_fl.n_sim == 281 && d_fl.t_iny_min == 134 && d_fl.t_iny_max == 253 &&
                           static_cast<std::size_t>(d_fl.t_iny_max) + 5 + 22 == d_fl.n_sim - 1;
        const bool ok_nv = d_nv.w_pre == 135 && d_nv.w_sorteo == 120 && d_nv.w_post == 27 &&
                           d_nv.n_sim == 282 && d_nv.t_iny_min == 135 && d_nv.t_iny_max == 254 &&
                           static_cast<std::size_t>(d_nv.t_iny_max) + 5 + 22 == d_nv.n_sim - 1;
        const bool ok_0  = d_0.w_pre == 9 && d_0.w_sorteo == 1 && d_0.w_post == 13 &&
                           d_0.n_sim == 23 && d_0.t_iny_min == 9 && d_0.t_iny_max == 9 &&
                           static_cast<std::size_t>(d_0.t_iny_max) + 3 + 10 == d_0.n_sim - 1;

        std::printf("  [C1] FLUJO(n_sim=%zu, t_iny=[%lld,%lld]) NIVEL(w_pre=%zu, n_sim=%zu) w0(n_sim=%zu)\n",
                    d_fl.n_sim, static_cast<long long>(d_fl.t_iny_min), static_cast<long long>(d_fl.t_iny_max),
                    d_nv.w_pre, d_nv.n_sim, d_0.n_sim);
        check(ok_fl && ok_nv && ok_0, "C1 dimensiones_simulacion_via_b exactas y t_iny_max + l + t_max_ef == n_sim - 1");
    }

    // -------------------------------------------------------------------------
    // C2: S2_N3 con t_max_ef real (< T_MAX=128) bajo ARL0=1e6 y evaluar_inyecciones_realizacion
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 128>;
        const auto x_cal = generar_ar1_unidad(4096, 0.85, 0.0, SEMILLA_CIEGA_NIVEL3 + 101ULL);

        // Derivar horizonte fisico t_max_ef < 128
        ConfiguracionObservador cd_ref{};
        cd_ref.integracion = Integracion::FLUJO;
        cd_ref.arl0_pasos_nativos = 450.0;
        cd_ref.medir_tau50 = false;
        cd_ref.escalas_max = 1;
        auto od_ref = std::make_unique<Obs>();
        od_ref->calibrar(cd_ref, x_cal);

        ConfiguracionLilyth cl_ref{};
        cl_ref.representacion = ModoRepresentacionLilyth::CRUDA;
        cl_ref.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl_ref.b_max_nulo = 8;
        cl_ref.hilos_nulo = 1;
        cl_ref.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 102ULL;
        cl_ref.usar_enmienda_nivel2 = true;
        cl_ref.usar_enmienda_nivel3 = true;

        auto lil_ref = std::make_unique<Lil>();
        check(lil_ref->calibrar(*od_ref, cl_ref, cd_ref, x_cal), "C2 calibrar lil_ref");

        ConfiguracionLilyth cl_ref_base = cl_ref;
        cl_ref_base.usar_enmienda_nivel3 = false;
        auto lil_ref_base = std::make_unique<Lil>();
        check(lil_ref_base->calibrar(*od_ref, cl_ref_base, cd_ref, x_cal), "C2 calibrar lil_ref_base");

        const std::size_t t_max_ef = lil_ref->configuracion().t_max_efectivo;
        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od_ref->ficha_nivel(0).nulo.tau_int)));
        const std::size_t l = 4;
        const std::size_t dist_borde = w_ancla + guarda + l + t_max_ef;
        const std::size_t dist_entre = w_ancla + guarda + l + t_max_ef + 16;
        const std::size_t dist_tmax  = w_ancla + guarda + l + 128;

        check(t_max_ef + 16 < 128 && dist_entre < dist_tmax, "C2 horizonte derivado t_max_ef + 16 < 128");

        // C2a: detector ARL0=1e6 (sin falsas alarmas espontaneas) + usar_enmienda_nivel3 true vs false
        ConfiguracionObservador cd_hi = cd_ref;
        cd_hi.arl0_pasos_nativos = 1.0e6;
        auto od_1e6 = std::make_unique<Obs>();
        od_1e6->calibrar(cd_hi, x_cal);

        ConfiguracionLilyth cl_1e6 = cl_ref;
        cl_1e6.t_max_efectivo = t_max_ef;
        auto lil_1e6 = std::make_unique<Lil>();
        check(lil_1e6->configurar(cl_1e6, cd_hi, *od_1e6, x_cal), "C2a configurar lil_1e6");

        const auto& fn = lil_ref->ficha_nulo_cruda();
        const std::vector<double> hist_prev(fn.modelo_lineal.p + 8u, 0.0);
        const std::vector<double> perfil{4.0, 6.0, 7.5, 8.5};
        std::vector<double> e_iny(l, 0.0);
        innovaciones(fn.modelo_lineal, hist_prev.data(), hist_prev.size(), perfil.data(), l, e_iny.data());

        const std::array<std::int64_t, 2> t_borde{250, static_cast<std::int64_t>(250 + dist_borde)};
        const std::array<std::int64_t, 2> t_entre{250, static_cast<std::int64_t>(250 + dist_entre)};
        const std::size_t n_sim = 1200;

        auto correr_c2a = [&](std::span<const std::int64_t> ts, bool usar_n3, std::uint64_t sem) {
            std::mt19937_64 rng(detalle_nulo::splitmix64(sem));
            return detalle_nulo::ejecutar_inyecciones_en_realizacion<128, RasgosDelta<Obs>::b_max>(
                *od_1e6, *lil_1e6, Integracion::FLUJO, false, false, cl_1e6.nulo_sin_memoria,
                fn.segmentos_calma, fn.modelo_lineal, lil_1e6->contexto_purga_nulo(), NaN,
                perfil, e_iny, +1, w_ancla, guarda, t_max_ef, ts, n_sim, rng, true, usar_n3);
        };

        const auto ra_borde = correr_c2a(t_borde, true,  SEMILLA_CIEGA_NIVEL3 + 201ULL);
        const auto ra_ent_3 = correr_c2a(t_entre, true,  SEMILLA_CIEGA_NIVEL3 + 202ULL);
        const auto ra_ent_0 = correr_c2a(t_entre, false, SEMILLA_CIEGA_NIVEL3 + 202ULL);

        std::printf("  [C2a] t_max_ef=%zu (T_MAX=128) | dist_borde=%zu(sol=%zu adm=%zu) | dist_entre=%zu n3=1(sol=%zu adm=%zu) n3=0(sol=%zu adm=%zu)\n",
                    t_max_ef, dist_borde, ra_borde.solapes, ra_borde.eventos_admitidos.size(),
                    dist_entre, ra_ent_3.solapes, ra_ent_3.eventos_admitidos.size(),
                    ra_ent_0.solapes, ra_ent_0.eventos_admitidos.size());
        check(ra_borde.solapes == 2 && ra_borde.eventos_admitidos.empty(),
              "C2a borde t_max_ef con n3=true descarta ambos eventos");
        check(ra_ent_3.solapes == 0 && ra_ent_3.eventos_admitidos.size() == 2,
              "C2a intermedio (< T_MAX=128) con n3=true admite ambos eventos");
        check(ra_ent_0.solapes == 2 && ra_ent_0.eventos_admitidos.empty(),
              "C2a Regla 6: mismo intermedio con n3=false sigue descartando con T_MAX=128");

        // C2b: passthrough completo via ObservadorLilyth::evaluar_inyecciones_realizacion
        bool hallada = false;
        std::uint64_t sem_limpia = 0;
        ResultadoInyeccionesRealizacion<128> rb_3{}, rb_0{};
        for (std::uint64_t s = 0; s < 32; ++s) {
            const std::uint64_t sem = SEMILLA_CIEGA_NIVEL3 + 2200ULL + s;
            const auto probe = lil_ref->evaluar_inyecciones_realizacion(
                false, hist_prev, perfil, +1, w_ancla, guarda, t_entre, n_sim, sem, false, true);
            if (probe.con_ajena_previa == 0 && probe.solapes == 0 && probe.eventos_admitidos.size() == 2) {
                hallada = true;
                sem_limpia = sem;
                rb_3 = probe;
                rb_0 = lil_ref_base->evaluar_inyecciones_realizacion(
                    false, hist_prev, perfil, +1, w_ancla, guarda, t_entre, n_sim, sem, false, true);
                break;
            }
        }
        std::printf("  [C2b] passthrough ObservadorLilyth (sem=%llu): hallada=%d | n3=1(sol=%zu adm=%zu) vs n3=0(sol=%zu adm=%zu)\n",
                    static_cast<unsigned long long>(sem_limpia), hallada ? 1 : 0,
                    rb_3.solapes, rb_3.eventos_admitidos.size(), rb_0.solapes, rb_0.eventos_admitidos.size());
        check(hallada && rb_3.solapes == 0 && rb_3.eventos_admitidos.size() == 2 &&
              rb_0.solapes == 2 && rb_0.eventos_admitidos.empty(),
              "C2b passthrough en ObservadorLilyth::evaluar_inyecciones_realizacion");
    }

    // -------------------------------------------------------------------------
    // C3: Alineacion ARFIMA mu_serie + mu en ejecutar_inyecciones_en_realizacion
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 64>;

        const auto z50 = generar_arfima1(4096, 0.40, 0.35, 50.0, 1.0, SEMILLA_CIEGA_NIVEL3 + 301ULL);
        std::vector<double> z1050 = z50;
        for (double& v : z1050) v += 1000.0;
        const auto mom50 = est::momentos(std::vector<std::vector<double>>{z50});

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 1.0e6;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od50 = std::make_unique<Obs>();
        auto od1050 = std::make_unique<Obs>();
        od50->calibrar(cd, z50);
        od1050->calibrar(cd, z1050);

        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
        cl.d_whittle_cruda = 0.35;
        cl.b_max_nulo = 2;
        cl.hilos_nulo = 1;
        cl.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 302ULL;
        cl.usar_enmienda_nivel2 = true;
        cl.usar_enmienda_nivel3 = true;

        auto lil50 = std::make_unique<Lil>();
        auto lil1050 = std::make_unique<Lil>();
        check(lil50->calibrar(*od50, cl, cd, z50), "C3 calibrar lil50");
        check(lil1050->calibrar(*od1050, cl, cd, z1050), "C3 calibrar lil1050");

        const auto& m50 = lil50->ficha_nulo_cruda().modelo_lineal;
        const auto& m1050 = lil1050->ficha_nulo_cruda().modelo_lineal;
        const std::size_t nh = std::min<std::size_t>(m50.n_hist_req, 1500);
        const std::vector<double> hp50(z50.begin() + static_cast<std::ptrdiff_t>(2000 - nh), z50.begin() + 2000);
        const std::vector<double> hp1050(z1050.begin() + static_cast<std::ptrdiff_t>(2000 - nh), z1050.begin() + 2000);

        constexpr std::size_t L = 4;
        std::vector<double> pf50(L, 0.0), pf1050(L, 0.0);
        for (std::size_t k = 0; k < L; ++k) {
            const double bump = (6.0 + 2.0 * static_cast<double>(k)) * mom50.sd;
            pf50[k] = z50[2000 + k] + bump;
            pf1050[k] = pf50[k] + 1000.0;
        }

        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od50->ficha_nivel(0).nulo.tau_int)));
        const std::array<std::int64_t, 1> t_arr{450};
        const std::size_t n_sim = 900;

        const auto r50 = lil50->evaluar_inyecciones_realizacion(
            false, hp50, pf50, +1, w_ancla, guarda, t_arr, n_sim, SEMILLA_CIEGA_NIVEL3 + 303ULL, false, true);
        const auto r1050 = lil1050->evaluar_inyecciones_realizacion(
            false, hp1050, pf1050, +1, w_ancla, guarda, t_arr, n_sim, SEMILLA_CIEGA_NIVEL3 + 303ULL, false, true);

        bool ok_c3 = (m50.valido && m1050.valido && r50.detectadas == 1 && r1050.detectadas == 1 &&
                      r50.eventos_admitidos.size() == 1 && r1050.eventos_admitidos.size() == 1);
        double dif_anc = NaN, dif_dy = NaN, dif_tau = NaN;
        if (ok_c3) {
            const auto& ev50 = r50.eventos_admitidos.front();
            const auto& ev1050 = r1050.eventos_admitidos.front();
            dif_anc = std::abs(ev1050.y_ancla - (ev50.y_ancla + 1000.0));
            dif_dy = std::abs(ev1050.dy - ev50.dy);
            const bool ambos_nan = std::isnan(ev50.tau_rest_consolidado) && std::isnan(ev1050.tau_rest_consolidado);
            dif_tau = ambos_nan ? 0.0 : std::abs(ev1050.tau_rest_consolidado - ev50.tau_rest_consolidado);
            ok_c3 = ok_c3 && (dif_anc <= 1e-8) && (dif_dy <= 1e-8) &&
                    (ambos_nan ? (ev50.tau_fin == ev1050.tau_fin) : (dif_tau == 0.0));
        }
        std::printf("  [C3] ARFIMA d=+0.35 mu_serie(%.3f -> %.3f): det=(%zu,%zu) dif_anc=%.3e dif_dy=%.3e dif_tau=%.3e\n",
                    m50.mu_serie, m1050.mu_serie, r50.detectadas, r1050.detectadas, dif_anc, dif_dy, dif_tau);
        check(ok_c3, "C3 traslacion mu_z=50->1050 en ejecutar_inyecciones_en_realizacion invariante <= 1e-8");
    }

    // -------------------------------------------------------------------------
    // C4: End-to-end ejecutar_via_b_evento() con n_sim derivado en FLUJO y NIVEL
    // -------------------------------------------------------------------------
    {
        auto probar_modo = [&](Integracion integ, const char* nom, std::uint64_t sem_base) {
            using Obs = ObservadorSaltos<24, 64, 120>;
            using Lil = ObservadorLilyth<Obs, 8, 64>;

            const auto inc_cal = generar_ar1_unidad(4096, 0.80, 0.0, sem_base);
            const auto inc_op  = generar_ar1_unidad(2048, 0.80, 0.0, sem_base + 1ULL);
            const auto mom = est::momentos(std::vector<std::vector<double>>{inc_cal});

            std::vector<double> x_cal(4096, 0.0), x_op(2048, 0.0);
            if (integ == Integracion::FLUJO) {
                x_cal = inc_cal;
                x_op = inc_op;
                for (std::size_t t = 600; t <= 602; ++t) x_op[t] += 9.0 * mom.sd;
                for (std::size_t t = 1300; t <= 1302; ++t) x_op[t] += 9.0 * mom.sd;
            } else {
                double lv = 100.0;
                for (std::size_t t = 0; t < 4096; ++t) { lv += inc_cal[t]; x_cal[t] = lv; }
                lv = 100.0;
                for (std::size_t t = 0; t < 2048; ++t) {
                    double dx = inc_op[t];
                    if (t == 600 || t == 601) dx += 10.0 * mom.sd;
                    if (t == 1300 || t == 1301) dx += 10.0 * mom.sd;
                    lv += dx;
                    x_op[t] = lv;
                }
            }

            ConfiguracionObservador cd{};
            cd.integracion = integ;
            cd.arl0_pasos_nativos = 450.0;
            cd.medir_tau50 = false;
            cd.escalas_max = 1;
            auto od = std::make_unique<Obs>();
            od->calibrar(cd, x_cal);

            ConfiguracionLilyth cl{};
            cl.representacion = ModoRepresentacionLilyth::CRUDA;
            cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
            cl.b_max_nulo = 4;
            cl.hilos_nulo = 1;
            cl.semilla_nulo = sem_base + 2ULL;
            cl.forzar_via_b = true;
            cl.usar_enmienda_nivel2 = true;
            cl.usar_enmienda_nivel3 = true;

            auto lil = std::make_unique<Lil>();
            check(lil->calibrar(*od, cl, cd, x_cal), "C4 calibrar");

            od->reiniciar_flujo();
            lil->reiniciar_flujo_conservando_purga();
            for (double v : x_op) {
                const auto sd = od->paso(v, NaN, true, 0, 0);
                (void)lil->paso(sd, v, true, NaN, 0, 0);
            }
            lil->completar_nulos();

            const std::size_t m_min = lil->ficha_nulo_cruda().m_min;
            bool ok_ev = false;
            ReferenciaNuloEvento ref_ok{};
            for (const auto& r : lil->cerradas_historico()) {
                const auto d = r.cruda.desenlace;
                if (d != DesenlaceRanura::RESTITUIDA && d != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
                const auto& ref = r.cruda.nulo_onsager;
                if (ref.evaluado && ref.nulo_por_inyeccion && ref.n_nulos_b >= m_min && std::isfinite(ref.u_pit)) {
                    ok_ev = true;
                    ref_ok = ref;
                    break;
                }
            }
            std::printf("  [C4 %s] t_max_ef=%zu cerradas=%zu eval=%d iny=%d n_b=%zu (m_min=%zu) intentos=%zu u_pit=%.4f\n",
                        nom, lil->configuracion().t_max_efectivo, lil->cerradas_historico().size(),
                        ref_ok.evaluado ? 1 : 0, ref_ok.nulo_por_inyeccion ? 1 : 0,
                        ref_ok.n_nulos_b, m_min, ref_ok.intentos_b, ref_ok.u_pit);
            check(ok_ev, nom);
        };

        probar_modo(Integracion::FLUJO, "FLUJO", SEMILLA_CIEGA_NIVEL3 + 401ULL);
        probar_modo(Integracion::NIVEL, "NIVEL", SEMILLA_CIEGA_NIVEL3 + 411ULL);
    }

    std::printf("\nFALLAS T3.2: %zu\n", fallas);
    return fallas == 0 ? 0 : 1;
}
