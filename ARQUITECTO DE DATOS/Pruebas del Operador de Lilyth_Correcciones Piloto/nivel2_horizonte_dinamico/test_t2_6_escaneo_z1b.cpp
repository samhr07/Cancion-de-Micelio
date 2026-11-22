// ============================================================================
// test_t2_6_escaneo_z1b.cpp -- Pruebas unitarias aisladas de T2.6:
//   1. Pre-centrado de mu_serie en puente_historico_ito.hpp (un paso e Ito bridge).
//   2. Protocolo hibrido Z1b (puente Ito r=2 en ARFIMA d=+0.35), e_iny exacto y
//      aislamiento de copia ArenaNoCopiableN2 en lil_local.
//   3. Cortes duros en Z1b: r=4 > r_max=2 y choque oculto +12*sigma en r=2.
//   4. Corte duro exacto por cambio de hist_seg_purgada_ en PURGADA (NIVEL+purga)
//      mientras CRUDA puentea r=2 en el mismo salto.
//   5. Seguridad estricta Z1a (t_ini_perf - 1): jamas lee t_cambio.
//   6. Veto estricto Z2 (r=1, FLAG_PUENTE_ITO y seg_purgada) exclusivo a via (b).
//   7. Perfil largo l=320 > 256 con K10 cero-heap (operator new == 0 en paso()).
// Semilla ciega: SEMILLA_CIEGA_NIVEL2 = 20261123ULL.
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

constexpr std::uint64_t SEMILLA_CIEGA_NIVEL2 = 20261123ULL;

// Genera realizacion estacionaria ARFIMA(1, d, 0) centrada en mu_serie con burn-in largo
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

template <class ObsLT>
typename ObsLT::Ranura buscar_ranura_evento(
    const ObsLT& lil,
    const typename ObsLT::SalidaPaso& out,
    bool& encontrada) {
    encontrada = false;
    if (out.abrio_ranura) {
        encontrada = true;
        return lil.ranuras()[out.indice_ranura_abierta];
    }
    if (out.ranuras_cerradas > 0) {
        encontrada = true;
        return lil.cerrada(out.ranuras_cerradas - 1);
    }
    return typename ObsLT::Ranura{};
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
    std::printf(" TEST T2.6: ESCANEO DINAMICO Z1b, Z1a SEGURO, VETO Z2 Y PRECALCULO e_iny\n");
    std::printf("================================================================================\n\n");

    // -------------------------------------------------------------------------
    // Test 1: Pre-centrado mu_serie en predecir_estacionario_un_paso & puente Ito
    // -------------------------------------------------------------------------
    {
        const auto z0 = generar_arfima1(1024, 0.40, 0.35, 50.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 11ULL);
        std::vector<double> z1(z0.size());
        for (std::size_t i = 0; i < z0.size(); ++i) z1[i] = z0[i] + 1000.0;

        ConfiguracionModeloLineal cm{};
        cm.modo = ModoNuloLineal::ARFIMA_WHITTLE;
        cm.d_externa = 0.35;
        cm.alfa = 0.05;
        cm.usar_enmienda_nivel2 = true;

        const auto m0 = ajustar({z0}, cm);
        const auto m1 = ajustar({z1}, cm);

        std::vector<double> h0(z0.begin() + 600, z0.begin() + 900);
        std::vector<double> h1(z1.begin() + 600, z1.begin() + 900);
        const double pred0 = piloto::detalle_ito::predecir_estacionario_un_paso(m0, h0.data(), h0.size());
        const double pred1 = piloto::detalle_ito::predecir_estacionario_un_paso(m1, h1.data(), h1.size());

        piloto::ConfiguracionPuenteHistoricoIto ci{};
        ci.integracion = Integracion::FLUJO;
        ci.r_max = 2;
        ci.alfa_hueco = 0.01;
        ci.n_semillas_rng = 0;
        const double xf0 = h0.back() + 0.5;
        const double xf1 = h1.back() + 0.5;
        const auto r0 = piloto::evaluar_y_reconstruir_hueco_ito(ci, m0, h0, h0.back(), xf0, 2);
        const auto r1 = piloto::evaluar_y_reconstruir_hueco_ito(ci, m1, h1, h1.back(), xf1, 2);
        const double dif_pred = std::abs((pred1 - 1000.0) - pred0);
        const double dif_z = std::abs(r1.z_hueco - r0.z_hueco);
        double max_dif_x = 0.0;
        for (std::size_t k = 0; k < 2; ++k) {
            max_dif_x = std::max(max_dif_x, std::abs((r1.x_media_cond[k] - 1000.0) - r0.x_media_cond[k]));
        }
        std::printf("  [T1] Pre-centrado puente Ito (mu=50 -> 1050): dif_pred=%.3e, dif_Z=%.3e, max_dif_x=%.3e\n",
                    dif_pred, dif_z, max_dif_x);
        check(m0.valido && m1.valido && r0.reconstruido && r1.reconstruido, "T1 modelos y puente validos");
        check(dif_pred <= 1e-11 && dif_z <= 1e-11 && max_dif_x <= 1e-11, "T1 invariancia traslacion <= 1e-11");
    }

    // Configuracion base compartida para ObservadorSaltos<24, 64, 60> en ARFIMA(1,+0.35,0)
    using ObsS = ObservadorSaltos<24, 64, 60>;
    using ObsL = ObservadorLilyth<ObsS, 8, 64>;

    const auto x_cal = generar_arfima1(1500, 0.40, 0.35, 25.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 23ULL);
    ConfiguracionObservador cd{};
    cd.integracion = Integracion::FLUJO;
    cd.arl0_pasos_nativos = 1.0e6;
    cd.medir_tau50 = false;
    cd.escalas_max = 1;
    auto od = std::make_unique<ObsS>();
    od->calibrar(cd, x_cal);
    check(od->escalas() >= 1, "Calibracion base ObsS");
    const auto cal_mom = est::momentos(std::vector<std::vector<double>>{x_cal});

    ConfiguracionLilyth cl_base{};
    cl_base.representacion = ModoRepresentacionLilyth::CRUDA;
    cl_base.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
    cl_base.d_whittle_cruda = 0.35;
    cl_base.b_max_nulo = 2;
    cl_base.usar_enmienda_nivel2 = true;
    cl_base.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_ITO_HISTORICO;
    cl_base.r_max_z1_piloto = 2;
    cl_base.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
    cl_base.r_max_z2_piloto = 0;

    auto lil = std::make_unique<ObsL>();
    check(lil->calibrar(*od, cl_base, cd, x_cal), "Calibracion base ObsL ARFIMA d=0.35");
    const auto& mod_c = lil->ficha_nulo_cruda().modelo_lineal;
    const double z_q_cal = est::norm_ppf(1.0 - 0.5 * cl_base.alfa);
    const std::size_t m_esp_cal = detalle_lineal::horizonte_inv_arfima_hosking(
        mod_c.d, (z_q_cal * z_q_cal) / static_cast<double>(mod_c.n), 4096);
    std::printf("  Base ARFIMA calibrada: n=%zu, p=%zu, m_inv=%zu (esp=%zu), n_hist_req=%zu, mu_serie=%.3f, arena_cap=%zu\n",
                mod_c.n, mod_c.p, mod_c.m_inv_arfima, m_esp_cal, mod_c.n_hist_req, mod_c.mu_serie, lil->arena_perfiles_capacidad());
    check(mod_c.valido && mod_c.m_inv_arfima == m_esp_cal && mod_c.n_hist_req == mod_c.p + m_esp_cal,
          "Modelo ARFIMA n_hist_req == p + M_pi_exact(n)");

    // Helper para correr una serie operativa con calma [0..499] + salto fuerte en t>=500
    auto generar_flujo_con_salto = [&](std::uint64_t sem_op, double salto_sigma, std::size_t n_tot = 560) {
        auto seq = generar_arfima1(n_tot, 0.40, 0.35, 25.0, 1.0, sem_op);
        for (std::size_t t = 500; t < n_tot; ++t) {
            seq[t] += salto_sigma * cal_mom.sd;
        }
        return seq;
    };

    // -------------------------------------------------------------------------
    // Test 2: Puente Z1b (r=2) mantiene n_hist_prev_disp == n_hist_req, e_iny exacto y lil_local = 0 copy
    // -------------------------------------------------------------------------
    {
        auto seq = generar_flujo_con_salto(SEMILLA_CIEGA_NIVEL2 + 101ULL, 10.0);
        std::vector<std::uint8_t> val(seq.size(), 1u);
        // Micro-hueco r=2 en [460, 461] dentro de Z1b (distancia ~40 < n_hist_req=87)
        val[460] = 0u;
        val[461] = 0u;

        od->reiniciar_flujo();
        lil->reiniciar_flujo_conservando_purga();
        ObsL::Ranura r_ev{};
        bool abrio = false;
        for (std::size_t t = 0; t < seq.size(); ++t) {
            auto s = od->paso(seq[t], NaN, val[t] != 0u, 0, 0);
            if (t < 495) {
                for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            }
            const auto o = lil->paso(s, seq[t], val[t] != 0u, NaN, 0, 0);
            if (!abrio && o.abrio_ranura) {
                r_ev = buscar_ranura_evento(*lil, o, abrio);
                break;
            }
        }
        const auto pf = lil->perfil_evento_n2(r_ev, false);
        const auto ei = lil->e_iny_evento_n2(r_ev, false);
        std::printf("  [T2] Puente Z1b r=2: t_c=%lld t_a=%lld n_disp=%zu (req=%zu) puentes=%zu trunc=%d len_pf=%zu len_ei=%zu\n",
                    static_cast<long long>(r_ev.t_cambio), static_cast<long long>(r_ev.t_apertura),
                    r_ev.n_hist_prev_disp_cruda, mod_c.n_hist_req,
                    r_ev.n_huecos_puenteados_z1b_cruda, r_ev.hist_prev_truncada_cruda ? 1 : 0,
                    pf.size(), ei.size());
        check(abrio && r_ev.t_cambio >= 480, "T2 evento real abierto cerca del salto");
        check(!r_ev.perfil_corrupto_z2_cruda, "T2 Z2 limpio");
        check(r_ev.n_hist_prev_disp_cruda == mod_c.n_hist_req, "T2 n_hist_prev_disp == n_hist_req tras puente Z1b r=2");
        check(r_ev.n_huecos_puenteados_z1b_cruda == 1, "T2 n_huecos_puenteados_z1b == 1");
        check(!r_ev.hist_prev_truncada_cruda, "T2 hist_prev_truncada == false");
        check(pf.size() > 0 && ei.size() == pf.size() &&
              std::isfinite(ei[0]) && std::isfinite(ei.back()),
              "T2 perfil y e_iny finitos almacenados en arena");

        // Verificar aislamiento de copia en lil_local cuando arena tiene datos escritos (> 0)
        ObsL lil_copia(*lil);
        std::printf("  [T2] Copia ArenaNoCopiableN2: main usados=%zu cap=%zu | copia usados=%zu cap=%zu\n",
                    lil->arena_perfiles_usados(), lil->arena_perfiles_capacidad(),
                    lil_copia.arena_perfiles_usados(), lil_copia.arena_perfiles_capacidad());
        check(lil->arena_perfiles_usados() > 0 && lil_copia.arena_perfiles_usados() == 0 &&
              lil_copia.arena_perfiles_capacidad() == 0,
              "T2 copia en lil_local copia 0 bytes de arena");
    }

    // -------------------------------------------------------------------------
    // Test 3: Cortes duros en Z1b (3a: r=4 > r_max=2; 3b: choque +12*sigma en r=2)
    // -------------------------------------------------------------------------
    {
        auto seq = generar_flujo_con_salto(SEMILLA_CIEGA_NIVEL2 + 202ULL, 10.0);
        std::vector<std::uint8_t> val(seq.size(), 1u);
        // Hueco largo r=4 en [460, 463] dentro de Z1b
        for (std::size_t k = 460; k <= 463; ++k) val[k] = 0u;

        od->reiniciar_flujo();
        lil->reiniciar_flujo_conservando_purga();
        ObsL::Ranura r_ev{};
        bool abrio = false;
        for (std::size_t t = 0; t < seq.size(); ++t) {
            auto s = od->paso(seq[t], NaN, val[t] != 0u, 0, 0);
            if (t < 495) {
                for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            }
            const auto o = lil->paso(s, seq[t], val[t] != 0u, NaN, 0, 0);
            if (!abrio && o.abrio_ranura) {
                r_ev = buscar_ranura_evento(*lil, o, abrio);
                break;
            }
        }
        const std::size_t l_ev = r_ev.len_perfil_cruda_n2;
        const std::int64_t t_ini_perf = r_ev.t_apertura - static_cast<std::int64_t>(l_ev) + 1;
        const std::size_t d_cut = abrio ? static_cast<std::size_t>(t_ini_perf - 464) : 0;
        std::printf("  [T3a] Corte duro r=4 [460..463]: t_ini_perf=%lld n_disp=%zu d_cut=%zu trunc=%d\n",
                    static_cast<long long>(t_ini_perf), r_ev.n_hist_prev_disp_cruda, d_cut,
                    r_ev.hist_prev_truncada_cruda ? 1 : 0);
        check(abrio && r_ev.hist_prev_truncada_cruda && r_ev.n_hist_prev_disp_cruda == d_cut,
              "T3a corte duro exacto en borde derecho de r=4");
    }

    {
        auto seq = generar_arfima1(760, 0.40, 0.35, 25.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 303ULL);
        for (std::size_t t = 700; t < seq.size(); ++t) seq[t] += 10.0 * cal_mom.sd;
        std::vector<std::uint8_t> val(seq.size(), 1u);
        // Micro-hueco r=2 [620..621] y choque puntual +12*sigma_e justo en t=622; ventana GLR L=64 cierra antes de t=700
        val[620] = 0u;
        val[621] = 0u;
        const double sd_e = std::sqrt(mod_c.sigma2_innov);
        seq[622] += 12.0 * sd_e;

        od->reiniciar_flujo();
        lil->reiniciar_flujo_conservando_purga();
        ObsL::Ranura r_ev{};
        bool abrio = false;
        for (std::size_t t = 0; t < seq.size(); ++t) {
            auto s = od->paso(seq[t], NaN, val[t] != 0u, 0, 0);
            if (t < 695) {
                for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            }
            const auto o = lil->paso(s, seq[t], val[t] != 0u, NaN, 0, 0);
            if (!abrio && o.abrio_ranura) {
                r_ev = buscar_ranura_evento(*lil, o, abrio);
                break;
            }
        }
        const std::size_t l_ev = r_ev.len_perfil_cruda_n2;
        const std::int64_t t_ini_perf = r_ev.t_apertura - static_cast<std::int64_t>(l_ev) + 1;
        const std::size_t d_shock = abrio ? static_cast<std::size_t>(t_ini_perf - 622) : 0;
        std::printf("  [T3b] Corte por choque +12*sigma en r=2 [620..621]: t_ini_perf=%lld n_disp=%zu d_shock=%zu puentes=%zu trunc=%d\n",
                    static_cast<long long>(t_ini_perf), r_ev.n_hist_prev_disp_cruda, d_shock,
                    r_ev.n_huecos_puenteados_z1b_cruda, r_ev.hist_prev_truncada_cruda ? 1 : 0);
        check(abrio && r_ev.n_hist_prev_disp_cruda == d_shock && r_ev.n_huecos_puenteados_z1b_cruda == 0 &&
              r_ev.hist_prev_truncada_cruda,
              "T3b break exacto en borde del choque |Z_hueco| > z_max");
    }

    // -------------------------------------------------------------------------
    // Test 4: Corte duro por hist_seg_purgada_ en PURGADA mientras CRUDA puentea r=2 (NIVEL + purga)
    // -------------------------------------------------------------------------
    {
        constexpr std::size_t BINS = 4;
        constexpr std::size_t N_CAL_NIV = 1500;
        const auto dy_cal = generar_arfima1(N_CAL_NIV, 0.40, 0.35, 0.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 404ULL);
        std::vector<double> x_cal_niv(N_CAL_NIV, 100.0);
        std::vector<std::size_t> f_cal(N_CAL_NIV, 0);
        std::vector<std::int64_t> d_cal(N_CAL_NIV, 0);
        for (std::size_t t = 1; t < N_CAL_NIV; ++t) {
            f_cal[t] = t % BINS;
            d_cal[t] = static_cast<std::int64_t>(t / BINS);
            const double amp_f = 1.0 + 0.25 * static_cast<double>(f_cal[t]);
            x_cal_niv[t] = x_cal_niv[t - 1] + amp_f * dy_cal[t];
        }

        ConfiguracionObservador cd_niv{};
        cd_niv.integracion = Integracion::NIVEL;
        cd_niv.purga = true;
        cd_niv.bins_fase = BINS;
        cd_niv.dias_minimos_purga = 5;
        cd_niv.arl0_pasos_nativos = 1.0e6;
        cd_niv.medir_tau50 = false;
        cd_niv.escalas_max = 1;
        auto od_niv = std::make_unique<ObsS>();
        od_niv->calibrar(cd_niv, x_cal_niv, {}, {}, f_cal, d_cal);
        check(od_niv->escalas() >= 1, "T4 calibrar ObsS NIVEL+purga");

        ConfiguracionLilyth cl_niv = cl_base;
        cl_niv.representacion = ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA;
        cl_niv.d_whittle_purgada = 0.35;
        cl_niv.d_whittle_cruda = 0.35;
        auto lil_niv = std::make_unique<ObsL>();
        check(lil_niv->calibrar(*od_niv, cl_niv, cd_niv, x_cal_niv, {}, f_cal, d_cal), "T4 calibrar ObsL NIVEL+purga");

        const auto dy_op = generar_arfima1(600, 0.40, 0.35, 0.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 405ULL);
        std::vector<double> x_op(600, 100.0);
        std::vector<std::uint8_t> val_op(600, 1u);
        for (std::size_t t = 1; t < 600; ++t) {
            const std::size_t fa = t % BINS;
            const double amp_f = 1.0 + 0.25 * static_cast<double>(fa);
            const double step_sig = (t >= 520) ? 10.0 : 0.0;
            x_op[t] = x_op[t - 1] + amp_f * (dy_op[t] + step_sig);
        }
        // Un NaN en x_op[480] crea r=2 en diferencias de CRUDA ([480..481]) y cambia hist_seg_purgada_ a S+1 en t>=480
        val_op[480] = 0u;

        od_niv->reiniciar_flujo();
        lil_niv->reiniciar_flujo_conservando_purga();
        ObsL::Ranura r_ev{};
        bool abrio = false;
        for (std::size_t t = 0; t < 600; ++t) {
            const std::size_t fa = t % BINS;
            const std::int64_t di = static_cast<std::int64_t>(t / BINS);
            auto s = od_niv->paso(x_op[t], NaN, val_op[t] != 0u, fa, di);
            if (t < 515) {
                for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            }
            const auto o = lil_niv->paso(s, x_op[t], val_op[t] != 0u, NaN, fa, di);
            if (!abrio && o.abrio_ranura) {
                r_ev = buscar_ranura_evento(*lil_niv, o, abrio);
                break;
            }
        }
        const std::size_t l_p = r_ev.len_perfil_purgada_n2;
        const std::int64_t t_ini_perf = r_ev.t_apertura - static_cast<std::int64_t>(l_p) + 1;
        const std::size_t d_seg = abrio ? static_cast<std::size_t>(t_ini_perf - 482) : 0;
        const std::size_t req_c = lil_niv->ficha_nulo_cruda().modelo_lineal.n_hist_req;
        std::printf("  [T4] NIVEL+purga NaN x[480]: t_ini_perf=%lld PURGADA(n=%zu d_seg=%zu trunc=%d) | CRUDA(n=%zu req=%zu puentes=%zu)\n",
                    static_cast<long long>(t_ini_perf),
                    r_ev.n_hist_prev_disp_purgada, d_seg, r_ev.hist_prev_truncada_purgada ? 1 : 0,
                    r_ev.n_hist_prev_disp_cruda, req_c, r_ev.n_huecos_puenteados_z1b_cruda);
        check(abrio && l_p > 0, "T4 evento NIVEL+purga abierto con perfil");
        check(r_ev.n_hist_prev_disp_purgada == d_seg && r_ev.hist_prev_truncada_purgada,
              "T4 PURGADA corta duro exacto por hist_seg_purgada_");
        check(r_ev.n_hist_prev_disp_cruda == req_c && r_ev.n_huecos_puenteados_z1b_cruda == 1 &&
              !r_ev.hist_prev_truncada_cruda,
              "T4 CRUDA puentea r=2 limpiamente hasta n_hist_req");

        // Sub-chequeo de T6 en NIVEL+purga: cambio de hist_seg_purgada_ dentro de Z2 ([520..525] con NaN en 522)
        od_niv->reiniciar_flujo();
        lil_niv->reiniciar_flujo_conservando_purga();
        bool z2_seg_corrupto = false;
        bool z2_seg_no_cierra = false;
        for (std::size_t t = 0; t < 560; ++t) {
            const std::size_t fa = t % BINS;
            const std::int64_t di = static_cast<std::int64_t>(t / BINS);
            const bool ok_z2 = (t != 522);
            auto s = od_niv->paso(x_op[t], NaN, ok_z2, fa, di);
            for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            if (t == 525) {
                s.nivel[0].alarma.hay = true;
                s.nivel[0].alarma.signo = 1;
                s.nivel[0].alarma.t = 525;
                s.nivel[0].alarma.t_cambio = 520;
                s.nivel[0].alarma.l = 6;
                s.nivel[0].alarma.g = 20.0;
                s.nivel[0].alarma.amplitud = 10.0;
                s.nivel[0].alarma.se_amplitud = 1.0;
            }
            const auto o = lil_niv->paso(s, x_op[t], ok_z2, NaN, fa, di);
            if (o.abrio_ranura) {
                bool enc = false;
                const auto rr = buscar_ranura_evento(*lil_niv, o, enc);
                z2_seg_corrupto = rr.perfil_corrupto_z2_purgada &&
                                  lil_niv->perfil_evento_n2(rr, true).empty() &&
                                  lil_niv->e_iny_evento_n2(rr, true).empty();
                z2_seg_no_cierra = rr.purgada.desenlace != DesenlaceRanura::NO_MEDIDA_PERFIL_CORRUPTO;
                break;
            }
        }
        std::printf("  [T6-seg] Z2 cambio hist_seg_purgada_: corrupto=%d sigue_abierta=%d\n",
                    z2_seg_corrupto ? 1 : 0, z2_seg_no_cierra ? 1 : 0);
        check(z2_seg_corrupto && z2_seg_no_cierra, "T6 cambio hist_seg_purgada_ en Z2 veta solo via (b) sin cerrar ranura");
    }

    // -------------------------------------------------------------------------
    // Test 5: Seguridad estricta en Z1a (t_cambio - 1): cero filtracion de t_cambio
    // -------------------------------------------------------------------------
    {
        const auto base = generar_arfima1(520, 0.40, 0.35, 25.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 505ULL);
        auto seqA = base;
        auto seqB = base;
        const double saltoA = 8.0 * cal_mom.sd;
        const double saltoB = 22.0 * cal_mom.sd;
        for (std::size_t t = 500; t < 520; ++t) {
            seqA[t] += saltoA;
            seqB[t] += saltoB;
        }
        auto correr_z1a = [&](const std::vector<double>& sq, double& hist_tail, double& e0) -> std::pair<std::int64_t, std::int64_t> {
            od->reiniciar_flujo();
            lil->reiniciar_flujo_conservando_purga();
            for (std::size_t t = 0; t < sq.size(); ++t) {
                const bool ok = (t != 499); // Hueco justo en Z1a (t_cambio - 1)
                auto s = od->paso(sq[t], NaN, ok, 0, 0);
                for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
                if (t == 500) {
                    s.nivel[0].alarma.hay = true;
                    s.nivel[0].alarma.signo = 1;
                    s.nivel[0].alarma.t = 500;
                    s.nivel[0].alarma.t_cambio = 500;
                    s.nivel[0].alarma.l = 1;
                    s.nivel[0].alarma.g = 25.0;
                    s.nivel[0].alarma.amplitud = sq[500] - 25.0;
                    s.nivel[0].alarma.se_amplitud = 1.0;
                }
                const auto o = lil->paso(s, sq[t], ok, NaN, 0, 0);
                if (o.abrio_ranura) {
                    bool enc = false;
                    const auto rr = buscar_ranura_evento(*lil, o, enc);
                    hist_tail = rr.n_hist_prev_cruda > 0 ? rr.hist_prev_cruda_est[rr.n_hist_prev_cruda - 1] : NaN;
                    const auto ei = lil->e_iny_evento_n2(rr, false);
                    e0 = ei.empty() ? NaN : ei[0];
                    return {rr.t_cambio, rr.t_apertura};
                }
            }
            return {-1, -1};
        };
        double htailA = NaN, htailB = NaN, e0A = NaN, e0B = NaN;
        const auto ta = correr_z1a(seqA, htailA, e0A);
        const auto tb = correr_z1a(seqB, htailB, e0B);
        const double dif_hist_z1a = std::abs(htailA - htailB);
        const double dif_e0 = std::abs((e0B - e0A) - (saltoB - saltoA));
        std::printf("  [T5] Z1a seguro (t_c=500, hueco 499): ta=(%lld,%lld) tb=(%lld,%lld) dif_hist=%.3e dif_de0=%.3e\n",
                    static_cast<long long>(ta.first), static_cast<long long>(ta.second),
                    static_cast<long long>(tb.first), static_cast<long long>(tb.second),
                    dif_hist_z1a, dif_e0);
        check(ta.first == 500 && ta.second == 500 && tb.first == 500 && tb.second == 500, "T5 Z1a ventanas iguales");
        check(dif_hist_z1a == 0.0, "T5 Z1a jamas filtra t_cambio (hist[499] bit a bit identico)");
        check(dif_e0 <= 1e-12, "T5 e_iny_B[0] - e_iny_A[0] == xB - xA exactamente");
    }

    // -------------------------------------------------------------------------
    // Test 6: Veto estricto Z2 (r=1 y FLAG_PUENTE_ITO) exclusivo a via (b)
    // -------------------------------------------------------------------------
    {
        auto correr_z2 = [&](bool con_nan_z2, bool con_flag_ito_z2) {
            const auto sq = generar_arfima1(530, 0.40, 0.35, 25.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 606ULL);
            od->reiniciar_flujo();
            lil->reiniciar_flujo_conservando_purga();
            for (std::size_t t = 0; t < sq.size(); ++t) {
                const bool ok = !(con_nan_z2 && t == 502);
                const std::uint8_t fl = (con_flag_ito_z2 && t == 502) ? ark::dataset::FLAG_PUENTE_ITO : 0u;
                auto s = od->paso(sq[t], NaN, ok, 0, 0);
                for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
                double x_t = sq[t];
                if (t >= 500 && t <= 505) x_t += 6.0 * cal_mom.sd;
                if (t == 505) {
                    s.nivel[0].alarma.hay = true;
                    s.nivel[0].alarma.signo = 1;
                    s.nivel[0].alarma.t = 505;
                    s.nivel[0].alarma.t_cambio = 500;
                    s.nivel[0].alarma.l = 6;
                    s.nivel[0].alarma.g = 20.0;
                    s.nivel[0].alarma.amplitud = 6.0 * cal_mom.sd;
                    s.nivel[0].alarma.se_amplitud = 1.0;
                }
                const auto o = lil->paso(s, x_t, ok, NaN, 0, 0, NaN, {}, fl);
                if (o.abrio_ranura) {
                    bool enc = false;
                    const auto rr = buscar_ranura_evento(*lil, o, enc);
                    const bool sin_perf = lil->perfil_evento_n2(rr, false).empty() && lil->e_iny_evento_n2(rr, false).empty();
                    const bool abierta_o_medible = rr.cruda.desenlace != DesenlaceRanura::NO_MEDIDA_PERFIL_CORRUPTO;
                    return std::pair<bool, bool>{rr.perfil_corrupto_z2_cruda && sin_perf, abierta_o_medible};
                }
            }
            return std::pair<bool, bool>{false, false};
        };
        const auto res_nan = correr_z2(true, false);
        const auto res_ito = correr_z2(false, true);
        std::printf("  [T6] Z2 veto r=1 NaN: corrupto=%d sigue_abierta=%d | FLAG_PUENTE_ITO: corrupto=%d sigue_abierta=%d\n",
                    res_nan.first ? 1 : 0, res_nan.second ? 1 : 0,
                    res_ito.first ? 1 : 0, res_ito.second ? 1 : 0);
        check(res_nan.first && res_nan.second, "T6 r=1 en Z2 veta solo via (b) sin cerrar ranura");
        check(res_ito.first && res_ito.second, "T6 FLAG_PUENTE_ITO en Z2 veta solo via (b) sin cerrar ranura");
    }

    // -------------------------------------------------------------------------
    // Test 7: Perfil largo l=320 > 256 con ObservadorSaltos<24,120,120>, J=3 (cap=480) y K10
    // -------------------------------------------------------------------------
    {
        using ObsS3 = ObservadorSaltos<24, 120, 120>;
        using ObsL3 = ObservadorLilyth<ObsS3, 8, 64>;
        const auto x_cal3 = generar_arfima1(1400, 0.40, 0.35, 10.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 707ULL);
        std::vector<double> x_cal_j3(10000, 0.0);
        {
            std::mt19937_64 r3(SEMILLA_CIEGA_NIVEL2 + 708ULL);
            std::normal_distribution<double> n01(0.0, 1.0);
            double cur = 0.0;
            for (std::size_t i = 0; i < x_cal_j3.size(); ++i) { cur = 0.40 * cur + n01(r3); x_cal_j3[i] = 10.0 + cur; }
        }
        ConfiguracionObservador cd3{};
        cd3.integracion = Integracion::FLUJO;
        cd3.arl0_pasos_nativos = 450.0;
        cd3.medir_tau50 = false;
        cd3.escalas_max = 3;
        auto od3 = std::make_unique<ObsS3>();
        od3->calibrar(cd3, x_cal_j3);
        check(od3->escalas() == 3, "T7 calibrate J=3");

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
        check(lil3->calibrar(*od3, cl3, cd3, x_cal3), "T7 calibrar ObsL l>256");

        const auto sq3 = generar_arfima1(1200, 0.40, 0.35, 10.0, 1.0, SEMILLA_CIEGA_NIVEL2 + 709ULL);
        od3->reiniciar_flujo();
        lil3->reiniciar_flujo_conservando_purga();

        ObsL3::Ranura r3{};
        bool abrio3 = false;
        g_news.store(0, std::memory_order_relaxed);
        g_contar_new.store(true, std::memory_order_release);
        for (std::size_t t = 0; t < sq3.size(); ++t) {
            auto s = od3->paso(sq3[t], NaN, true, 0, 0);
            for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            double x_t = sq3[t] + (t >= 800 ? 4.0 : 0.0);
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
        const std::size_t n_req3 = r3.n_hist_prev_disp_cruda;
        std::vector<double> e_ref(320, 0.0);
        if (n_req3 > 0 && pf3.size() == 320) {
            std::vector<double> h_ref(sq3.begin() + (800 - static_cast<std::ptrdiff_t>(n_req3)), sq3.begin() + 800);
            innovaciones(mod3, h_ref.data(), h_ref.size(), pf3.data(), pf3.size(), e_ref.data());
        }
        double max_dif_ei = 0.0;
        for (std::size_t i = 0; i < ei3.size() && i < e_ref.size(); ++i) {
            max_dif_ei = std::max(max_dif_ei, std::abs(ei3[i] - e_ref[i]));
        }
        std::printf("  [T7] Perfil l=320 (>256): cap_adapt=%zu len_pf=%zu len_ei=%zu news=%zu desborde=%zu max_dif_einy=%.3e\n",
                    lil3->cap_perfil_adaptable(), pf3.size(), ei3.size(), news_paso,
                    lil3->contadores().arena_perfiles_desbordada, max_dif_ei);
        check(abrio3 && lil3->cap_perfil_adaptable() == 480, "T7 cap_perfil_adaptable == 480");
        check(pf3.size() == 320 && ei3.size() == 320, "T7 perfil y e_iny completos l=320 sin truncar en 256");
        check(news_paso == 0, "T7 K10 cero reservas dinamicas en paso()");
        check(lil3->contadores().arena_perfiles_desbordada == 0, "T7 sin desborde de arena");
        check(max_dif_ei <= 1e-13, "T7 e_iny precalculado coincide con innovaciones() directa <= 1e-13");
        check(sizeof(ObsL3::Ranura) < 16384, "T7 sizeof(Ranura) < 16384");
    }

    std::printf("\nFALLAS T2.6: %zu\n", fallas);
    return fallas == 0 ? 0 : 1;
}
