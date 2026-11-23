// ============================================================================
// test_t3_4_parada_v4.cpp -- Pruebas unitarias aisladas de T3.4:
//   C1. Corte efectivo de cola V4 en AR(1) multi-barra l=20: paradas_tempranas>0,
//       tau_med_v4 < tau_med_base, escala_coincide>0 (base paradas==0).
//   C2. ARFIMA(1, d=+0.35, 0) con V4 activado: invariancia mu_z=50->1050 <= 1e-10.
//   C3. PasoCrudoNulo bit a bit vs llevar_a_crudo + coherencia V4 en NIVEL+Purga.
//   C4. GanchoHijas en linea en ejecutar_inyecciones_en_realizacion: replica
//       contraria en t+2 corta la ranura (cortadas==1, detectadas==0).
//   C5. E2E completar_nulos() con V4 activada (AR(1) escalas=3, forzar_via_b=true).
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

struct GanchoReplicaContraria {
    bool disparado = false;
    double amp_sig = -40.0;
    std::size_t t_min = 250;

    void preparar(std::mt19937_64&, std::vector<double>&, std::size_t, std::size_t) const noexcept {}

    template <class SalidaLil, class Lil>
    void en_paso(std::size_t t, const SalidaLil& o, const Lil&, std::mt19937_64&,
                 std::vector<double>& eps, std::vector<double>& eps_nulo, std::size_t burn, std::size_t n_sim) noexcept {
        if (!disparado && t >= t_min && o.abrio_ranura && (t + 6u) < n_sim) {
            disparado = true;
            eps[burn + t + 5u] += amp_sig;
            eps_nulo[burn + t + 5u] += amp_sig;
        }
    }
};

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
    std::printf(" TEST T3.4: PARADA TEMPRANA V4_OPCION_B_PARADA_TEMP, PasoCrudoNulo Y GanchoHijas\n");
    std::printf("================================================================================\n\n");

    // -------------------------------------------------------------------------
    // C1: Corte efectivo de cola V4 en AR(1) multi-barra l=20
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 64>;

        const auto x_cal = generar_ar1_unidad(4096, 0.80, 0.0, SEMILLA_CIEGA_NIVEL3 + 4001ULL);
        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 450.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 3;
        auto od = std::make_unique<Obs>();
        od->calibrar(cd, x_cal);

        ConfiguracionLilyth cl_v4{};
        cl_v4.representacion = ModoRepresentacionLilyth::CRUDA;
        cl_v4.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl_v4.b_max_nulo = 4;
        cl_v4.hilos_nulo = 1;
        cl_v4.t_max_efectivo = 16;
        cl_v4.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 4002ULL;
        cl_v4.usar_enmienda_nivel2 = true;
        cl_v4.usar_enmienda_nivel3 = true;

        ConfiguracionLilyth cl_base = cl_v4;
        cl_base.usar_enmienda_nivel3 = false;

        auto lil_v4 = std::make_unique<Lil>();
        auto lil_base = std::make_unique<Lil>();
        check(lil_v4->calibrar(*od, cl_v4, cd, x_cal), "C1 calibrar lil_v4");
        check(lil_base->calibrar(*od, cl_base, cd, x_cal), "C1 calibrar lil_base");

        const std::size_t t_max_ef = lil_v4->configuracion().t_max_efectivo;
        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od->ficha_nivel(0).nulo.tau_int)));
        const std::vector<double> hist_prev(lil_v4->ficha_nulo_cruda().modelo_lineal.p + 8u, 0.0);
        std::vector<double> perfil(20, 12.0);
        for (std::size_t k = 0; k < 6; ++k) perfil[k] = 4.0 + 0.5 * static_cast<double>(k);

        const std::array<std::int64_t, 1> t_arr{250};
        const std::size_t n_sim = 600;
        std::size_t paradas_v4 = 0, paradas_base = 0, esc_v4 = 0;
        std::vector<double> taus_v4, taus_base;

        for (std::size_t k = 0; k < 40; ++k) {
            const std::uint64_t sem = SEMILLA_CIEGA_NIVEL3 + 4100ULL + k;
            const auto r_v4 = lil_v4->evaluar_inyecciones_realizacion(
                false, hist_prev, perfil, +1, w_ancla, guarda, t_arr, n_sim, sem, false, true, 0);
            const auto r_bs = lil_base->evaluar_inyecciones_realizacion(
                false, hist_prev, perfil, +1, w_ancla, guarda, t_arr, n_sim, sem, false, true, 0);
            paradas_v4 += r_v4.paradas_tempranas;
            paradas_base += r_bs.paradas_tempranas;
            esc_v4 += r_v4.escala_coincide;
            for (const auto& ev : r_v4.eventos_admitidos) taus_v4.push_back(detalle_nulo::codificar_tau(ev.tau_rest_consolidado, t_max_ef));
            for (const auto& ev : r_bs.eventos_admitidos) taus_base.push_back(detalle_nulo::codificar_tau(ev.tau_rest_consolidado, t_max_ef));
        }

        const double med_v4 = taus_v4.empty() ? NaN : est::cuantil(taus_v4, 0.5);
        const double med_base = taus_base.empty() ? NaN : est::cuantil(taus_base, 0.5);
        std::printf("  [C1 l=20] t_max_ef=%zu | V4(paradas=%zu adm=%zu med=%.2f esc=%zu) | Base(paradas=%zu adm=%zu med=%.2f)\n",
                    t_max_ef, paradas_v4, taus_v4.size(), med_v4, esc_v4, paradas_base, taus_base.size(), med_base);
        check(paradas_v4 > 0 && paradas_base == 0, "C1 paradas_tempranas activa solo con n3=true");
        check(std::isfinite(med_v4) && std::isfinite(med_base) && med_v4 < med_base,
              "C1 V4 reduce tau_med respecto a la inyeccion completa");
        check(esc_v4 > 0, "C1 telemetria escala_coincide > 0");
    }

    // -------------------------------------------------------------------------
    // C2: ARFIMA(1, d=+0.35, 0) con V4 activado: invariancia mu_z=50->1050 <= 1e-10
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 64>;

        const auto z50 = generar_arfima1(4096, 0.40, 0.35, 50.0, 1.0, SEMILLA_CIEGA_NIVEL3 + 4201ULL);
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
        cl.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 4202ULL;
        cl.usar_enmienda_nivel2 = true;
        cl.usar_enmienda_nivel3 = true;

        auto lil50 = std::make_unique<Lil>();
        auto lil1050 = std::make_unique<Lil>();
        check(lil50->calibrar(*od50, cl, cd, z50), "C2 calibrar lil50");
        check(lil1050->calibrar(*od1050, cl, cd, z1050), "C2 calibrar lil1050");

        const std::size_t nh = std::min<std::size_t>(lil50->ficha_nulo_cruda().modelo_lineal.n_hist_req, 1500);
        const std::vector<double> hp50(z50.begin() + static_cast<std::ptrdiff_t>(2000 - nh), z50.begin() + 2000);
        const std::vector<double> hp1050(z1050.begin() + static_cast<std::ptrdiff_t>(2000 - nh), z1050.begin() + 2000);

        constexpr std::size_t L = 12;
        std::vector<double> pf50(L, 0.0), pf1050(L, 0.0);
        for (std::size_t k = 0; k < L; ++k) {
            const double bump = (k < 4 ? (6.0 + 3.0 * static_cast<double>(k)) : 15.0) * mom50.sd;
            pf50[k] = z50[2000 + k] + bump;
            pf1050[k] = pf50[k] + 1000.0;
        }

        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od50->ficha_nivel(0).nulo.tau_int)));
        const std::array<std::int64_t, 1> t_arr{450};
        const std::size_t n_sim = 900;

        const auto r50 = lil50->evaluar_inyecciones_realizacion(
            false, hp50, pf50, +1, w_ancla, guarda, t_arr, n_sim, SEMILLA_CIEGA_NIVEL3 + 4203ULL, false, true, 0);
        const auto r1050 = lil1050->evaluar_inyecciones_realizacion(
            false, hp1050, pf1050, +1, w_ancla, guarda, t_arr, n_sim, SEMILLA_CIEGA_NIVEL3 + 4203ULL, false, true, 0);

        bool ok_c2 = (r50.paradas_tempranas == 1 && r1050.paradas_tempranas == 1 &&
                      r50.detectadas == 1 && r1050.detectadas == 1 &&
                      r50.eventos_admitidos.size() == 1 && r1050.eventos_admitidos.size() == 1);
        double dif_anc = NaN, dif_dy = NaN, max_dif_eta = NaN;
        if (ok_c2) {
            const auto& e50 = r50.eventos_admitidos.front();
            const auto& e1050 = r1050.eventos_admitidos.front();
            dif_anc = std::abs(e1050.y_ancla - (e50.y_ancla + 1000.0));
            dif_dy = std::abs(e1050.dy - e50.dy);
            max_dif_eta = 0.0;
            ok_c2 = ok_c2 && (e50.pasos_observados == e1050.pasos_observados) && (e50.tau_fin == e1050.tau_fin);
            for (std::size_t tau = 0; tau < e50.pasos_observados && tau < e50.eta.size(); ++tau) {
                const double d = std::abs(e1050.eta[tau] - e50.eta[tau]);
                if (d > max_dif_eta) max_dif_eta = d;
            }
            const bool ambos_nan = std::isnan(e50.tau_rest_consolidado) && std::isnan(e1050.tau_rest_consolidado);
            ok_c2 = ok_c2 && (dif_anc <= 1e-10) && (dif_dy <= 1e-10) && (max_dif_eta <= 1e-10) &&
                    (ambos_nan || e50.tau_rest_consolidado == e1050.tau_rest_consolidado);
        }
        std::printf("  [C2 ARFIMA+V4] paradas=(%zu,%zu) det=(%zu,%zu) dif_anc=%.3e dif_dy=%.3e max|dif_eta|=%.3e\n",
                    r50.paradas_tempranas, r1050.paradas_tempranas, r50.detectadas, r1050.detectadas,
                    dif_anc, dif_dy, max_dif_eta);
        check(ok_c2, "C2 invariancia mu_z=50->1050 <= 1e-10 en ARFIMA + V4");
    }

    // -------------------------------------------------------------------------
    // C3: PasoCrudoNulo bit a bit vs llevar_a_crudo + coherencia V4 en NIVEL+Purga
    // -------------------------------------------------------------------------
    {
        constexpr std::size_t BINS = 4;
        using Obs = ObservadorSaltos<24, 64, 120, 8>;
        using Lil = ObservadorLilyth<Obs, 8, 64>;

        const auto inc_cal = generar_ar1_unidad(4096, 0.70, 0.0, SEMILLA_CIEGA_NIVEL3 + 4301ULL);
        const std::array<double, BINS> escala_f{0.8, 1.3, 1.0, 1.6};
        std::vector<double> x_cal(4096, 0.0);
        std::vector<std::size_t> f_cal(4096, 0);
        std::vector<std::int64_t> d_cal(4096, 0);
        double lvl = 100.0;
        for (std::size_t t = 0; t < 4096; ++t) {
            f_cal[t] = t % BINS;
            d_cal[t] = static_cast<std::int64_t>(t / BINS);
            if (t > 0) lvl += inc_cal[t] * escala_f[f_cal[t]];
            x_cal[t] = lvl;
        }

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::NIVEL;
        cd.purga = true;
        cd.bins_fase = BINS;
        cd.dias_minimos_purga = 5;
        cd.arl0_pasos_nativos = 450.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<Obs>();
        od->calibrar(cd, x_cal, {}, {}, f_cal, d_cal);

        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::PURGADA;
        cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl.b_max_nulo = 2;
        cl.hilos_nulo = 1;
        cl.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 4302ULL;
        cl.usar_enmienda_nivel2 = true;
        cl.usar_enmienda_nivel3 = true;

        auto lil = std::make_unique<Lil>();
        check(lil->calibrar(*od, cl, cd, x_cal, {}, f_cal, d_cal), "C3 calibrar NIVEL+Purga");

        // 3a. Equivalencia bit a bit PasoCrudoNulo vs llevar_a_crudo
        const auto est = generar_ar1_unidad(512, 0.70, 0.0, SEMILLA_CIEGA_NIVEL3 + 4303ULL);
        std::vector<double> x_lac;
        std::vector<std::size_t> f_lac;
        std::vector<std::int64_t> d_lac;
        std::size_t sf_lac = 0;
        detalle_nulo::llevar_a_crudo<8>(
            Integracion::NIVEL, true, lil->contexto_purga_nulo(), est, 100.0, 3, 0, x_lac, f_lac, d_lac, sf_lac);

        detalle_nulo::PasoCrudoNulo<8> pc(Integracion::NIVEL, true, lil->contexto_purga_nulo(), 100.0, 3, 0);
        std::size_t sf_pc = 0;
        double max_dif_pc = 0.0;
        bool fd_ok = true;
        for (std::size_t t = 0; t < est.size(); ++t) {
            std::size_t ft = 0;
            std::int64_t dt = 0;
            const double xp = pc.paso(t, est[t], ft, dt, sf_pc);
            if (ft != f_lac[t] || dt != d_lac[t]) fd_ok = false;
            const double d = std::abs(xp - x_lac[t]);
            if (!std::isfinite(d) || d > max_dif_pc) max_dif_pc = d;
        }
        check(fd_ok && sf_pc == sf_lac && max_dif_pc == 0.0, "C3a PasoCrudoNulo identico bit a bit a llevar_a_crudo");

        // 3b. Inyeccion multi-barra NIVEL+Purga con parada temprana V4
        const std::vector<double> hp(lil->ficha_nulo_purgada().modelo_lineal.p + 8u, 0.0);
        const std::vector<double> pf_lvl{5.0, 7.0, 9.0, 11.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0};
        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od->ficha_nivel(0).nulo.tau_int)));
        const std::array<std::int64_t, 1> t_arr{450};
        const auto r_np = lil->evaluar_inyecciones_realizacion(
            true, hp, pf_lvl, +1, w_ancla, guarda, t_arr, 900, SEMILLA_CIEGA_NIVEL3 + 4304ULL, false, true, 0);

        bool fin_np = (r_np.paradas_tempranas > 0 && r_np.detectadas == 1 &&
                       r_np.sin_factor_purga == 0 && r_np.eventos_admitidos.size() == 1);
        if (fin_np) {
            const auto& ev = r_np.eventos_admitidos.front();
            fin_np = std::isfinite(ev.y_ancla) && std::isfinite(ev.dy) && ev.pasos_observados > 0;
            for (std::size_t tau = 0; tau < ev.pasos_observados && fin_np; ++tau) {
                fin_np = std::isfinite(ev.eta[tau]);
            }
        }
        std::printf("  [C3 NIVEL+Purga] max_dif_pc=%.3e sf=(%zu,%zu) | paradas=%zu det=%zu sf_sim=%zu\n",
                    max_dif_pc, sf_lac, sf_pc, r_np.paradas_tempranas, r_np.detectadas, r_np.sin_factor_purga);
        check(fin_np, "C3b V4 en NIVEL+Purga coherente y finito con sin_factor_purga==0");
    }

    // -------------------------------------------------------------------------
    // C4: GanchoHijas en linea en ejecutar_inyecciones_en_realizacion
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 4, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 128>;

        const auto x_cal = generar_ar1_unidad(4096, 0.92, 0.0, SEMILLA_CIEGA_NIVEL3 + 4401ULL);
        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 1.0e6;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<Obs>();
        od->calibrar(cd, x_cal);

        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl.b_max_nulo = 2;
        cl.hilos_nulo = 1;
        cl.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 4402ULL;
        cl.usar_enmienda_nivel2 = true;
        cl.usar_enmienda_nivel3 = true;

        auto lil = std::make_unique<Lil>();
        check(lil->calibrar(*od, cl, cd, x_cal), "C4 calibrar");

        const auto& fn = lil->ficha_nulo_cruda();
        const std::size_t t_max_ef = lil->configuracion().t_max_efectivo;
        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od->ficha_nivel(0).nulo.tau_int)));
        const std::vector<double> hp(fn.modelo_lineal.p + 8u, 0.0);
        const std::vector<double> pf{5.0, 7.0, 9.0, 11.0};
        std::vector<double> ei(pf.size(), 0.0);
        innovaciones(fn.modelo_lineal, hp.data(), hp.size(), pf.data(), pf.size(), ei.data());
        const std::array<std::int64_t, 1> ts{250};
        const std::size_t n_sim = 900;

        std::mt19937_64 rng_lim(detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL3 + 4403ULL));
        const auto r_lim = detalle_nulo::ejecutar_inyecciones_en_realizacion<128, RasgosDelta<Obs>::b_max>(
            *od, *lil, Integracion::FLUJO, false, false, cl.nulo_sin_memoria,
            fn.segmentos_calma, fn.modelo_lineal, lil->contexto_purga_nulo(), NaN,
            pf, ei, +1, w_ancla, guarda, t_max_ef, ts, n_sim, rng_lim, true, true, 0);

        std::mt19937_64 rng_h(detalle_nulo::splitmix64(SEMILLA_CIEGA_NIVEL3 + 4403ULL));
        GanchoReplicaContraria hook{false, -40.0, 250};
        const auto r_h = detalle_nulo::ejecutar_inyecciones_en_realizacion<128, RasgosDelta<Obs>::b_max>(
            *od, *lil, Integracion::FLUJO, false, false, cl.nulo_sin_memoria,
            fn.segmentos_calma, fn.modelo_lineal, lil->contexto_purga_nulo(), NaN,
            pf, ei, +1, w_ancla, guarda, t_max_ef, ts, n_sim, rng_h, true, true, 0, hook);

        std::printf("  [C4 GanchoHijas] limpio(det=%zu cort=%zu) | replica_contraria(det=%zu cort=%zu)\n",
                    r_lim.detectadas, r_lim.cortadas, r_h.detectadas, r_h.cortadas);
        check(r_lim.detectadas == 1 && r_lim.cortadas == 0, "C4 sin hijas admite el evento");
        check(r_h.cortadas == 1 && r_h.detectadas == 0, "C4 replica contraria en linea corta la ranura");
    }

    // -------------------------------------------------------------------------
    // C5: E2E completar_nulos() con V4 en AR(1) phi=0.80, escalas=3, forzar_via_b
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 64>;

        const auto x_cal = generar_ar1_unidad(4096, 0.80, 0.0, SEMILLA_CIEGA_NIVEL3 + 4501ULL);
        auto x_op = generar_ar1_unidad(4096, 0.80, 0.0, SEMILLA_CIEGA_NIVEL3 + 4502ULL);
        const auto mom = est::momentos(std::vector<std::vector<double>>{x_cal});
        for (std::size_t t = 1000; t <= 1006; ++t) x_op[t] += (4.0 + 0.5 * static_cast<double>(t - 1000)) * mom.sd;
        for (std::size_t t = 2500; t <= 2506; ++t) x_op[t] -= (4.0 + 0.5 * static_cast<double>(t - 2500)) * mom.sd;

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 450.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 3;
        auto od = std::make_unique<Obs>();
        od->calibrar(cd, x_cal);

        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl.b_max_nulo = 8;
        cl.hilos_nulo = 1;
        cl.semilla_nulo = SEMILLA_CIEGA_NIVEL3 + 4503ULL;
        cl.forzar_via_b = true;
        cl.usar_enmienda_nivel2 = true;
        cl.usar_enmienda_nivel3 = true;

        auto lil = std::make_unique<Lil>();
        check(lil->calibrar(*od, cl, cd, x_cal), "C5 calibrar");
        od->reiniciar_flujo();
        lil->reiniciar_flujo_conservando_purga();
        for (double v : x_op) {
            const auto sd = od->paso(v, NaN, true, 0, 0);
            (void)lil->paso(sd, v, true, NaN, 0, 0);
        }
        lil->completar_nulos();

        std::size_t n_eval = 0, paradas_b = 0, esc_b = 0, sol_b = 0;
        for (const auto& r : lil->cerradas_historico()) {
            const auto d = r.cruda.desenlace;
            if (d != DesenlaceRanura::RESTITUIDA && d != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
            const auto& ref = r.cruda.nulo_onsager;
            if (ref.evaluado) ++n_eval;
            paradas_b += ref.paradas_tempranas_b;
            esc_b += ref.escala_coincide_b;
            sol_b += ref.solapes_b;
        }
        std::printf("  [C5 E2E] n_eval=%zu paradas_b=%zu esc_coincide_b=%zu sol_b=%zu\n",
                    n_eval, paradas_b, esc_b, sol_b);
        check(n_eval > 0 && paradas_b > 0 && esc_b > 0 && sol_b == 0,
              "C5 E2E completar_nulos() ejerce V4 y telemetria con sol_b==0");
    }

    std::printf("\nFALLAS T3.4: %zu\n", fallas);
    return fallas == 0 ? 0 : 1;
}
