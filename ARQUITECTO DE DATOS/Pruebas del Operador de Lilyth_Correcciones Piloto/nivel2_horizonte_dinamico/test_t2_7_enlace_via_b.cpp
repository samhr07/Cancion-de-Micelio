// ============================================================================
// test_t2_7_enlace_via_b.cpp -- Pruebas unitarias aisladas de T2.7:
//   1. tamano_trabajo_simular() + construir_estacionaria_nula() con n=18000 (> BUF_STACK_MAX).
//   2. Enlace E2E completar_nulos() via (b) en ARFIMA(1,+0.35,0) (mu_serie=50.0, Z1b r=2).
//   3. Veto de via (b) cuando n_hist_prev_disp < mod.p (AR(4), n_hist_prev_disp=2 < 4).
//   4. Veto Z2 exclusivo a via (b) sin bloquear via (a) (forzar_via_b true vs false).
//   5. Via (b) consumiendo e_iny_evento_n2() precalculado con perfil largo l=280 (>256).
// Semilla ciega: SEMILLA_CIEGA_NIVEL2 = 20261123ULL.
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
#include <string_view>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;

constexpr std::uint64_t SEMILLA_CIEGA_NIVEL2 = 20261123ULL;

std::vector<double> generar_arfima1(
    std::size_t n, double phi1, double d, double mu_serie, double sigma_e, std::uint64_t sem) {
    const std::size_t burn = 1024;
    const std::size_t tot = burn + n;
    std::vector<double> pi_inv(tot, 0.0);
    elw::coeficientes(-d, tot, pi_inv.data());
    std::vector<double> u(tot, 0.0);
    std::vector<double> x(n, 0.0);
    std::mt19937_64 rng(sem);
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

std::vector<double> generar_ar_orden(
    std::size_t n, std::span<const double> phi, double mu_serie, double sigma_e, std::uint64_t sem) {
    const std::size_t burn = 512;
    const std::size_t tot = burn + n;
    std::vector<double> u(tot, 0.0);
    std::vector<double> x(n, 0.0);
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> nd(0.0, sigma_e);
    for (std::size_t s = 0; s < tot; ++s) {
        double val = nd(rng);
        for (std::size_t i = 1; i <= phi.size() && i <= s; ++i) {
            val += phi[i - 1] * u[s - i];
        }
        u[s] = val;
        if (s >= burn) x[s - burn] = mu_serie + val;
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
    std::printf(" TEST T2.7: ENLACE EN evaluar_una_referencia(), VETO <p Y BUFFER 2*tot\n");
    std::printf("================================================================================\n\n");

    // -------------------------------------------------------------------------
    // Check 1: tamano_trabajo_simular() + construir_estacionaria_nula() n = 18000 (> BUF_STACK_MAX)
    // -------------------------------------------------------------------------
    {
        const auto x_cal = generar_arfima1(1024, 0.40, 0.35, 50.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 17ULL);
        ConfiguracionModeloLineal cm_arf{};
        cm_arf.modo = ModoNuloLineal::ARFIMA_WHITTLE;
        cm_arf.d_externa = 0.35;
        cm_arf.alfa = 0.05;
        cm_arf.usar_enmienda_nivel2 = true;
        const auto m_arf = ajustar({x_cal}, cm_arf);

        ConfiguracionModeloLineal cm_ar{};
        cm_ar.modo = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cm_ar.alfa = 0.05;
        cm_ar.usar_enmienda_nivel2 = true;
        const auto m_ar = ajustar({x_cal}, cm_ar);

        constexpr std::size_t N_GRANDE = 18000;
        const std::size_t work_arf = tamano_trabajo_simular(m_arf, N_GRANDE);
        const std::size_t work_ar = tamano_trabajo_simular(m_ar, N_GRANDE);

        std::mt19937_64 rng1(SEMILLA_CIEGA_NIVEL2 + 18ULL);
        std::vector<double> est;
        detalle_nulo::construir_estacionaria_nula(true, LILYTH_NULO_PERMUTACION, {x_cal}, m_arf, N_GRANDE, rng1, est);

        bool todo_finito = (est.size() == N_GRANDE);
        for (double v : est) {
            if (!std::isfinite(v)) { todo_finito = false; break; }
        }
        double sum_cola = 0.0;
        const std::size_t ini_cola = detalle_lineal::BUF_STACK_MAX;
        const std::size_t n_cola = N_GRANDE - ini_cola;
        for (std::size_t t = ini_cola; t < N_GRANDE; ++t) sum_cola += est[t];
        const double mean_cola = sum_cola / static_cast<double>(n_cola);
        double var_cola = 0.0;
        for (std::size_t t = ini_cola; t < N_GRANDE; ++t) {
            const double d = est[t] - mean_cola;
            var_cola += d * d;
        }
        const double sd_cola = std::sqrt(var_cola / static_cast<double>(n_cola));
        std::printf("  [C1] tam_trab_arf=%zu (esp=%zu) tam_trab_ar=%zu | cola[16384..17999]: mu=%.3f (mu_serie=%.3f), sd=%.3f\n",
                    work_arf, 2u * (m_arf.burn_in + N_GRANDE), work_ar, mean_cola, m_arf.mu_serie, sd_cola);
        check(m_arf.valido && work_arf == 2u * (m_arf.burn_in + N_GRANDE), "C1 tamano_trabajo_simular ARFIMA == 2*(burn+n)");
        check(m_ar.valido && work_ar == m_ar.p, "C1 tamano_trabajo_simular AR == p");
        check(todo_finito && sd_cola > 0.5 && std::abs(mean_cola - m_arf.mu_serie) < 5.0,
              "C1 cola > BUF_STACK_MAX=16384 finita, no truncada y centrada en mu_serie");
    }

    // -------------------------------------------------------------------------
    // Check 2: E2E completar_nulos() via (b) en ARFIMA(1,+0.35,0) con mu_serie=50 y puente Z1b r=2
    // -------------------------------------------------------------------------
    {
        using ObsS = ObservadorSaltos<24, 64, 60>;
        using ObsL = ObservadorLilyth<ObsS, 8, 64>;
        const auto x_cal = generar_arfima1(1500, 0.40, 0.35, 50.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 21ULL);
        const auto mom = est::momentos(std::vector<std::vector<double>>{x_cal});

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 1.0e5;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<ObsS>();
        od->calibrar(cd, x_cal);
        check(od->escalas() >= 1, "C2 calibrar ObsS");

        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
        cl.d_whittle_cruda = 0.35;
        cl.b_max_nulo = 2;
        cl.usar_enmienda_nivel2 = true;
        cl.forzar_via_b = true;
        cl.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO;
        cl.r_max_z1_piloto = 2;
        cl.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cl.r_max_z2_piloto = 0;

        auto lil = std::make_unique<ObsL>();
        check(lil->calibrar(*od, cl, cd, x_cal), "C2 calibrar ObsL ARFIMA d=+0.35");

        auto seq = generar_arfima1(620, 0.40, 0.35, 50.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 22ULL);
        for (std::size_t t = 500; t <= 503; ++t) seq[t] += 8.0 * mom.sd; // Pulso corto -> restituye y cierra
        std::vector<std::uint8_t> val(seq.size(), 1u);
        val[460] = 0u;
        val[461] = 0u; // micro-hueco Z1b r=2

        od->reiniciar_flujo();
        lil->reiniciar_flujo_conservando_purga();
        for (std::size_t t = 0; t < seq.size(); ++t) {
            auto s = od->paso(seq[t], NaN, val[t] != 0u, 0, 0);
            if (t < 495 || t > 510) {
                for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            }
            (void)lil->paso(s, seq[t], val[t] != 0u, NaN, 0, 0);
        }
        lil->completar_nulos();
        const auto& mod = lil->ficha_nulo_cruda().modelo_lineal;
        const auto& hist = lil->cerradas_historico();
        bool hay_ev = false;
        ReferenciaNuloEvento ref_ons{};
        for (const auto& rc : hist) {
            if (rc.t_cambio >= 495 &&
                (rc.cruda.desenlace == DesenlaceRanura::RESTITUIDA || rc.cruda.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE)) {
                hay_ev = true;
                ref_ons = rc.cruda.nulo_onsager;
                break;
            }
        }
        std::printf("  [C2] ARFIMA E2E (Z1b r=2): hay_ev=%d evaluado=%d iny=%d sin_b=%d n_b=%zu hist_disp=%zu (req=%zu) puentes=%zu trunc=%d\n",
                    hay_ev ? 1 : 0, ref_ons.evaluado ? 1 : 0, ref_ons.nulo_por_inyeccion ? 1 : 0,
                    ref_ons.sin_nulo_b ? 1 : 0, ref_ons.n_nulos_b,
                    ref_ons.n_hist_prev_disp, mod.n_hist_req,
                    ref_ons.n_huecos_puenteados_z1b, ref_ons.hist_prev_truncada_en_corte ? 1 : 0);
        check(hay_ev && ref_ons.nulo_por_inyeccion && ref_ons.evaluado && !ref_ons.sin_nulo_b,
              "C2 via (b) evalua con exito el evento ARFIMA con mu_serie=50");
        check(ref_ons.n_nulos_b >= lil->ficha_nulo_cruda().m_min, "C2 n_nulos_b >= m_min");
        check(ref_ons.n_hist_prev_disp == mod.n_hist_req && ref_ons.n_huecos_puenteados_z1b == 1 &&
              !ref_ons.hist_prev_truncada_en_corte,
              "C2 metadatos ReferenciaNuloEvento exactos");
    }

    // -------------------------------------------------------------------------
    // Check 3: Veto n_hist_prev_disp < mod.p (AR(4), p_fit == 4, historia disponible == 2)
    // -------------------------------------------------------------------------
    {
        using ObsS = ObservadorSaltos<24, 64, 60>;
        using ObsL = ObservadorLilyth<ObsS, 8, 64>;
        const std::array<double, 4> phi4{0.25, 0.05, 0.05, 0.15};
        const auto x_cal = generar_ar_orden(4096, phi4, 0.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 31ULL);
        const auto mom = est::momentos(std::vector<std::vector<double>>{x_cal});

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 1.0e5;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<ObsS>();
        od->calibrar(cd, x_cal);

        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl.b_max_nulo = 2;
        cl.usar_enmienda_nivel2 = true;
        cl.forzar_via_b = true;
        cl.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::CORTE_SUFIJO;
        cl.r_max_z1_piloto = 0;
        cl.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cl.r_max_z2_piloto = 0;

        auto lil = std::make_unique<ObsL>();
        check(lil->calibrar(*od, cl, cd, x_cal), "C3 calibrar AR(4)");
        const auto& mod4 = lil->ficha_nulo_cruda().modelo_lineal;
        check(mod4.valido && mod4.p == 4, "C3 modelo AR(4) con p == 4");

        auto seq = generar_ar_orden(620, phi4, 0.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 32ULL);
        for (std::size_t t = 500; t <= 503; ++t) seq[t] += 8.0 * mom.sd;
        od->reiniciar_flujo();
        lil->reiniciar_flujo_conservando_purga();
        for (std::size_t t = 0; t < seq.size(); ++t) {
            const bool ok = (t != 497); // Deja exactamente 2 barras [498, 499] antes de t_cambio=500
            auto s = od->paso(seq[t], NaN, ok, 0, 0);
            for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            if (t == 500) {
                s.nivel[0].alarma.hay = true;
                s.nivel[0].alarma.signo = 1;
                s.nivel[0].alarma.t = 500;
                s.nivel[0].alarma.t_cambio = 500;
                s.nivel[0].alarma.l = 1;
                s.nivel[0].alarma.g = 20.0;
                s.nivel[0].alarma.amplitud = 8.0 * mom.sd;
                s.nivel[0].alarma.se_amplitud = 1.0;
            }
            (void)lil->paso(s, seq[t], ok, NaN, 0, 0);
        }
        lil->completar_nulos();
        check(!lil->cerradas_historico().empty(), "C3 evento cerrado");
        const auto ref3 = lil->cerradas_historico().empty() ? ReferenciaNuloEvento{}
                                                            : lil->cerradas_historico().front().cruda.nulo_onsager;
        std::printf("  [C3] Veto hist<p (p=%zu): n_disp=%zu evaluado=%d sin_b=%d trunc=%d motivo='%s'\n",
                    mod4.p, ref3.n_hist_prev_disp, ref3.evaluado ? 1 : 0,
                    ref3.sin_nulo_b ? 1 : 0, ref3.hist_prev_truncada_en_corte ? 1 : 0, ref3.motivo);
        check(ref3.n_hist_prev_disp == 2 && ref3.sin_nulo_b && !ref3.evaluado && ref3.hist_prev_truncada_en_corte &&
              std::string_view(ref3.motivo) == "via (b): historia contigua menor que p",
              "C3 veta exactamente con motivo 'via (b): historia contigua menor que p'");
    }

    // -------------------------------------------------------------------------
    // Check 4: Veto Z2 exclusivo a via (b) sin bloquear via (a)
    // -------------------------------------------------------------------------
    {
        using ObsS = ObservadorSaltos<24, 64, 60>;
        using ObsL = ObservadorLilyth<ObsS, 8, 64>;
        const std::array<double, 1> phi1{0.40};
        const auto x_cal = generar_ar_orden(6000, phi1, 0.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 41ULL);

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 80.0; // Excedencia abundante (>= m_min) en via (a) con b_max_nulo = 16
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<ObsS>();
        od->calibrar(cd, x_cal);

        auto correr_z2_contraste = [&](bool forzar_b) {
            ConfiguracionLilyth cl{};
            cl.representacion = ModoRepresentacionLilyth::CRUDA;
            cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
            cl.b_max_nulo = 16;
            cl.usar_enmienda_nivel2 = true;
            cl.forzar_via_b = forzar_b;
            cl.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO;
            cl.r_max_z1_piloto = 2;
            cl.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
            cl.r_max_z2_piloto = 0;

            auto lil_z = std::make_unique<ObsL>();
            check(lil_z->calibrar(*od, cl, cd, x_cal), "C4 calibrar");

            auto seq = generar_ar_orden(620, phi1, 0.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 42ULL);
            for (std::size_t t = 500; t <= 505; ++t) seq[t] += 5.0; // Salto modesto para que esp.size() >= m_min
            od->reiniciar_flujo();
            lil_z->reiniciar_flujo_conservando_purga();
            for (std::size_t t = 0; t < seq.size(); ++t) {
                auto s = od->paso(seq[t], NaN, true, 0, 0);
                for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
                if (t == 505) {
                    s.nivel[0].alarma.hay = true;
                    s.nivel[0].alarma.signo = 1;
                    s.nivel[0].alarma.t = 505;
                    s.nivel[0].alarma.t_cambio = 500;
                    s.nivel[0].alarma.l = 6;
                    s.nivel[0].alarma.g = 2.0; // g moderado para asegurar excedencia espontanea en (a)
                    s.nivel[0].alarma.amplitud = 5.0;
                    s.nivel[0].alarma.se_amplitud = 1.0;
                }
                const std::uint8_t fl = (t == 502) ? ark::dataset::FLAG_PUENTE_ITO : 0u;
                (void)lil_z->paso(s, seq[t], true, NaN, 0, 0, NaN, {}, fl);
            }
            lil_z->completar_nulos();
            check(!lil_z->cerradas_historico().empty(), "C4 evento cerrado");
            return lil_z->cerradas_historico().empty() ? ObsL::Ranura{} : lil_z->cerradas_historico().front();
        };

        const auto r_b = correr_z2_contraste(true);
        const auto r_a = correr_z2_contraste(false);
        std::printf("  [C4] Z2 FLAG_PUENTE_ITO: forzar_b=true(sin_b=%d eval=%d motivo='%s') | forzar_b=false(eval=%d iny=%d sin_b=%d n_a=%zu)\n",
                    r_b.cruda.nulo_onsager.sin_nulo_b ? 1 : 0, r_b.cruda.nulo_onsager.evaluado ? 1 : 0, r_b.cruda.nulo_onsager.motivo,
                    r_a.cruda.nulo_onsager.evaluado ? 1 : 0, r_a.cruda.nulo_onsager.nulo_por_inyeccion ? 1 : 0,
                    r_a.cruda.nulo_onsager.sin_nulo_b ? 1 : 0, r_a.cruda.nulo_onsager.n_nulos_a);
        check(r_b.perfil_corrupto_z2_cruda && r_b.cruda.nulo_onsager.sin_nulo_b && !r_b.cruda.nulo_onsager.evaluado &&
              std::string_view(r_b.cruda.nulo_onsager.motivo) == "via (b): perfil Z2 corrupto",
              "C4 Z2 corrupto veta via (b) bajo forzar_via_b");
        check(r_a.perfil_corrupto_z2_cruda && r_a.cruda.desenlace != DesenlaceRanura::NO_MEDIDA_PERFIL_CORRUPTO &&
              r_a.cruda.nulo_onsager.evaluado && !r_a.cruda.nulo_onsager.nulo_por_inyeccion &&
              !r_a.cruda.nulo_onsager.sin_nulo_b && r_a.cruda.nulo_onsager.n_nulos_a >= 8,
              "C4 Z2 corrupto NO bloquea via (a) cuando hay excedencia suficiente");
    }

    // -------------------------------------------------------------------------
    // Check 5: Via (b) con perfil largo l=280 (>256) desde arena N2
    // -------------------------------------------------------------------------
    {
        using ObsS3 = ObservadorSaltos<24, 120, 120>;
        using ObsL3 = ObservadorLilyth<ObsS3, 8, 64>;
        const std::array<double, 1> phi1{0.40};
        const auto x_cal = generar_ar_orden(10000, phi1, 10.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 51ULL);

        ConfiguracionObservador cd3{};
        cd3.integracion = Integracion::FLUJO;
        cd3.arl0_pasos_nativos = 450.0;
        cd3.medir_tau50 = false;
        cd3.escalas_max = 3;
        auto od3 = std::make_unique<ObsS3>();
        od3->calibrar(cd3, x_cal);
        check(od3->escalas() == 3, "C5 calibrar J=3");

        ConfiguracionLilyth cl3{};
        cl3.representacion = ModoRepresentacionLilyth::CRUDA;
        cl3.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl3.b_max_nulo = 2;
        cl3.usar_enmienda_nivel2 = true;
        cl3.forzar_via_b = true;
        cl3.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO;
        cl3.r_max_z1_piloto = 2;
        cl3.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cl3.r_max_z2_piloto = 0;

        auto lil3 = std::make_unique<ObsL3>();
        check(lil3->calibrar(*od3, cl3, cd3, x_cal), "C5 calibrar ObsL3");

        auto sq = generar_ar_orden(1200, phi1, 10.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 52ULL);
        for (std::size_t t = 800; t <= 1079; ++t) sq[t] += 3.0; // l = 280 > 256, luego vuelve a calma y restituye

        od3->reiniciar_flujo();
        lil3->reiniciar_flujo_conservando_purga();
        for (std::size_t t = 0; t < sq.size(); ++t) {
            auto s = od3->paso(sq[t], NaN, true, 0, 0);
            for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            if (t == 1079) {
                s.nivel[2].alarma.hay = true;
                s.nivel[2].alarma.signo = 1;
                s.nivel[2].alarma.t = 1079;
                s.nivel[2].alarma.t_cambio = 800;
                s.nivel[2].alarma.l = 70;
                s.nivel[2].alarma.g = 20.0;
                s.nivel[2].alarma.amplitud = 3.0;
                s.nivel[2].alarma.se_amplitud = 0.5;
            }
            (void)lil3->paso(s, sq[t], true, NaN, 0, 0);
        }
        check(!lil3->cerradas_historico().empty(), "C5 evento largo cerrado");
        const std::size_t len_pf = lil3->cerradas_historico().empty()
            ? 0
            : lil3->perfil_evento_n2(lil3->cerradas_historico().front(), false).size();
        const std::size_t len_ei = lil3->cerradas_historico().empty()
            ? 0
            : lil3->e_iny_evento_n2(lil3->cerradas_historico().front(), false).size();

        lil3->completar_nulos();
        const auto ref5 = lil3->cerradas_historico().empty()
            ? ReferenciaNuloEvento{}
            : lil3->cerradas_historico().front().cruda.nulo_onsager;
        const auto& mod5 = lil3->ficha_nulo_cruda().modelo_lineal;
        std::printf("  [C5] Perfil l=%zu (>256) via (b): len_ei=%zu iny=%d eval=%d intentos=%zu n_b=%zu hist_disp=%zu (req=%zu) motivo='%s'\n",
                    len_pf, len_ei, ref5.nulo_por_inyeccion ? 1 : 0, ref5.evaluado ? 1 : 0,
                    ref5.intentos_b, ref5.n_nulos_b, ref5.n_hist_prev_disp, mod5.n_hist_req, ref5.motivo);
        check(len_pf == 280 && len_ei == 280, "C5 arena N2 almacena l=280 completo para perfil y e_iny");
        check(ref5.nulo_por_inyeccion && ref5.evaluado && !ref5.sin_nulo_b &&
              ref5.intentos_b >= lil3->ficha_nulo_cruda().m_min &&
              ref5.n_nulos_b >= lil3->ficha_nulo_cruda().m_min &&
              ref5.n_hist_prev_disp == mod5.n_hist_req,
              "C5 completar_nulos() ejecuta via (b) sobre perfil largo l=280 sin truncamiento");
    }

    std::printf("\nFALLAS T2.7: %zu\n", fallas);
    return fallas == 0 ? 0 : 1;
}
