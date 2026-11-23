// ============================================================================
// sonda_admision_propia_via_b.cpp -- Fase 4.0: Sonda Empirica Multiescala
//   Nivel 4 (solo lectura sobre el operador actual, sin modificar cabeceras):
//   Evalua P0 (Nivel 3 base), P1a (escala en admision), P1b (escala en V4 + adm),
//   P2 (g_disparo >= g_real) y P3 (escala + g >= g_real) en trayectorias acopladas
//   con 3 escalas reales activas (n_cal = 12000, od->escalas() == 3, j in {0,1,2}):
//     Banco 1 (M1_AR1_3ESC)   : AR(1) phi=0.80, 24 semillas, S3 + INTERCAMBIABLE
//     Banco 2 (M2_ARFIMA_3ESC): ARFIMA(1,d=+0.35,0) phi1=0.40 mu=50, 16 sem, S3 + INTERCAMBIABLE
//     Banco 3 (M3_DIAG_GRUESAS): Diagnostico alta precision N_mc=200 por escala j=0,1,2
// Semilla exclusiva de sonda: SEMILLA_SONDA_NIVEL4 = 2026112500ULL.
// ============================================================================

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <algorithm>
#include <array>
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

using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsLil64 = ObservadorLilyth<ObsDelta, 8, 64>;
using ObsLil128 = ObservadorLilyth<ObsDelta, 8, 128>;

constexpr std::uint64_t SEMILLA_SONDA_NIVEL4 = 2026112500ULL;

enum PolIdx : std::size_t {
    P0_NIVEL3_BASE = 0,
    P1A_ESCALA_ADM = 1,
    P1B_ESCALA_V4_ADM = 2,
    P2_EXCEDENCIA_G = 3,
    P3_ESCALA_Y_G = 4,
    N_POL = 5
};

constexpr std::array<const char*, N_POL> NOMBRE_POL{
    "P0_NIVEL3_BASE",
    "P1a_ESCALA_ADM",
    "P1b_ESCALA_V4_ADM",
    "P2_EXCEDENCIA_G",
    "P3_ESCALA_Y_G"
};

struct SlotCerradoSim {
    std::int64_t t_ap = -1;
    int signo = 0;
    int escala = -1;
    double g = NaN;
    double y_ancla = NaN;
    std::uint8_t des = 0;
    double tau_cod = NaN;
};

struct SalidaTrayectoria {
    bool v4_disparo = false;
    bool vio_otra_escala_sig = false;
    double x_pre_sig = 0.0;
    std::vector<SlotCerradoSim> slots;
};

std::vector<double> generar_ar1_estacionario(
    std::size_t n, double phi, double sd_eps, std::uint64_t sem) {
    std::vector<double> x(n, 0.0);
    if (n == 0) return x;
    std::mt19937_64 rng(detalle_nulo::splitmix64(sem));
    std::normal_distribution<double> n01(0.0, 1.0);
    const double std_est = sd_eps / std::sqrt(std::max(1e-12, 1.0 - phi * phi));
    double prev = std_est * n01(rng);
    x[0] = prev;
    for (std::size_t t = 1; t < n; ++t) {
        prev = phi * prev + sd_eps * n01(rng);
        x[t] = prev;
    }
    return x;
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

template <class ObsLilT>
SalidaTrayectoria simular_trayectoria(
    const ObsDelta& od_cal,
    const ObsLilT& lil_plantilla,
    const ModeloLinealNulo& m_lin,
    std::span<const double> e_iny,
    std::span<const double> eps_nulo_in,
    std::span<const double> pi_inv,
    bool es_arfima,
    double mu_ref,
    std::size_t burn,
    std::size_t n_sim,
    std::int64_t t0,
    int signo_real,
    int escala_real,
    std::size_t t_max_ef,
    bool v4_exige_misma_escala) {

    SalidaTrayectoria st{};
    const std::size_t l = e_iny.size();
    const std::size_t total = burn + n_sim;
    std::vector<double> eps(eps_nulo_in.begin(), eps_nulo_in.end());
    const std::size_t u_t0 = static_cast<std::size_t>(t0);
    for (std::size_t k = 0; k < l && (u_t0 + k) < n_sim; ++k) {
        eps[burn + u_t0 + k] = e_iny[k];
    }

    std::vector<double> u(total, 0.0);
    for (std::size_t idx = 0; idx < burn; ++idx) {
        double acc = eps[idx];
        const std::size_t imax = std::min(m_lin.p, idx);
        for (std::size_t i = 1; i <= imax; ++i) acc += m_lin.phi[i - 1] * u[idx - i];
        u[idx] = acc;
    }

    const bool usar_mu_serie = es_arfima && std::isfinite(m_lin.mu_serie);
    const std::int64_t t1 = t0 + static_cast<std::int64_t>(l);
    const std::int64_t t_fin = t1 - 1;

    auto d_local = std::make_unique<ObsDelta>(od_cal);
    d_local->reiniciar_flujo();
    auto lil_local = std::make_unique<ObsLilT>(lil_plantilla);
    lil_local->reiniciar_flujo_conservando_purga();

    std::vector<double> x_sim(n_sim, 0.0);
    bool parada = false;

    for (std::size_t t = 0; t < n_sim; ++t) {
        const std::size_t idx = burn + t;
        double acc = eps[idx];
        const std::size_t imax = std::min(m_lin.p, idx);
        for (std::size_t i = 1; i <= imax; ++i) acc += m_lin.phi[i - 1] * u[idx - i];
        u[idx] = acc;

        double xt = 0.0;
        if (es_arfima) {
            double s_arf = 0.0;
            if (usar_mu_serie) {
                for (std::size_t k = 0; k <= idx; ++k) s_arf += pi_inv[k] * (m_lin.mu + u[idx - k]);
                xt = m_lin.mu_serie + s_arf;
            } else {
                for (std::size_t k = 0; k <= idx; ++k) s_arf += pi_inv[k] * u[idx - k];
                xt = m_lin.mu + s_arf;
            }
        } else {
            xt = m_lin.mu + u[idx];
        }
        x_sim[t] = xt;

        const auto sd = d_local->paso(xt, NaN, true, 0, 0);
        const auto o = lil_local->paso(sd, xt, true, NaN, 0, 0);

        for (std::size_t c = 0; c < o.ranuras_cerradas; ++c) {
            const auto& rc = lil_local->cerrada(c);
            SlotCerradoSim sl{};
            sl.t_ap = rc.t_apertura;
            sl.signo = rc.signo;
            sl.escala = rc.escala_origen;
            sl.g = rc.g_disparo;
            sl.y_ancla = rc.cruda.y_ancla;
            sl.des = static_cast<std::uint8_t>(rc.cruda.desenlace);
            if (sl.des == 2 || sl.des == 4) {
                sl.tau_cod = detalle_nulo::codificar_tau(rc.cruda.tau_rest_consolidado, t_max_ef);
            }
            st.slots.push_back(sl);
        }

        if (o.abrio_ranura && o.indice_ranura_abierta < lil_local->ranuras().size()) {
            const auto& ran_ab = lil_local->ranuras()[o.indice_ranura_abierta];
            if (ran_ab.ocupada && ran_ab.t_apertura == o.t && ran_ab.signo == signo_real) {
                if (o.t >= t0 && o.t <= t1 && ran_ab.escala_origen != escala_real) {
                    st.vio_otra_escala_sig = true;
                }
                const bool esc_ok = (!v4_exige_misma_escala) || (ran_ab.escala_origen == escala_real);
                if (!parada && esc_ok && o.t >= t0 && o.t < t_fin) {
                    parada = true;
                    st.v4_disparo = true;
                    const std::size_t tr_ini = static_cast<std::size_t>(o.t) + 1u;
                    const std::size_t tr_fin = std::min<std::size_t>(u_t0 + l, n_sim);
                    for (std::size_t tr = tr_ini; tr < tr_fin; ++tr) {
                        eps[burn + tr] = eps_nulo_in[burn + tr];
                    }
                }
            }
        }
    }

    for (const auto& r : lil_local->ranuras()) {
        if (!r.ocupada) continue;
        SlotCerradoSim sl{};
        sl.t_ap = r.t_apertura;
        sl.signo = r.signo;
        sl.escala = r.escala_origen;
        sl.g = r.g_disparo;
        sl.y_ancla = r.cruda.y_ancla;
        sl.des = static_cast<std::uint8_t>(r.cruda.desenlace);
        st.slots.push_back(sl);
    }

    if (t0 >= 1 && static_cast<std::size_t>(t0 - 1) < n_sim) {
        st.x_pre_sig = (x_sim[static_cast<std::size_t>(t0 - 1)] - mu_ref) * static_cast<double>(signo_real);
    }
    return st;
}

struct SeleccionSlot {
    bool hay_propia = false;
    bool admitida = false;
    int escala = -1;
    double g = NaN;
    double y_ancla = NaN;
    double tau_cod = NaN;
};

SeleccionSlot evaluar_politica_slot(
    const SalidaTrayectoria& st,
    std::int64_t t0,
    std::size_t l,
    int signo_real,
    int escala_real,
    double g_real,
    PolIdx pol) noexcept {

    SeleccionSlot out{};
    const std::int64_t t1 = t0 + static_cast<std::int64_t>(l);
    const bool buscar_con_escala = (pol == P1B_ESCALA_V4_ADM);

    const SlotCerradoSim* propia = nullptr;
    for (const auto& ci : st.slots) {
        if (ci.t_ap >= t0 && ci.t_ap <= t1 && ci.signo == signo_real) {
            if (buscar_con_escala && ci.escala != escala_real) continue;
            if (propia == nullptr || ci.t_ap < propia->t_ap) propia = &ci;
        }
    }
    if (propia == nullptr) return out;

    out.hay_propia = true;
    out.escala = propia->escala;
    out.g = propia->g;
    out.y_ancla = propia->y_ancla;
    out.tau_cod = propia->tau_cod;

    const bool extraible = (propia->des == 2u || propia->des == 4u) && std::isfinite(propia->tau_cod);
    if (!extraible) return out;

    switch (pol) {
        case P0_NIVEL3_BASE:
        case P1B_ESCALA_V4_ADM:
            out.admitida = true;
            break;
        case P1A_ESCALA_ADM:
            out.admitida = (propia->escala == escala_real);
            break;
        case P2_EXCEDENCIA_G:
            out.admitida = (std::isfinite(propia->g) && propia->g >= g_real);
            break;
        case P3_ESCALA_Y_G:
            out.admitida = (propia->escala == escala_real && std::isfinite(propia->g) && propia->g >= g_real);
            break;
        default:
            break;
    }
    return out;
}

struct AcumEscalaPol {
    std::size_t intentos_tot = 0;
    std::size_t k_med = 0;
    std::size_t k_pdet_pos = 0;
    std::size_t k_eval = 0;
    std::size_t n_adm = 0;
    std::size_t n_esc_match = 0;
    double sum_x_pre_sig = 0.0;
    double sum_g_diff = 0.0;
    double sum_u_ev = 0.0;
    std::size_t n_u_ev = 0;

    void sumar(const AcumEscalaPol& o) noexcept {
        intentos_tot += o.intentos_tot;
        k_med += o.k_med;
        k_pdet_pos += o.k_pdet_pos;
        k_eval += o.k_eval;
        n_adm += o.n_adm;
        n_esc_match += o.n_esc_match;
        sum_x_pre_sig += o.sum_x_pre_sig;
        sum_g_diff += o.sum_g_diff;
        sum_u_ev += o.sum_u_ev;
        n_u_ev += o.n_u_ev;
    }
};

struct ResumenSemillaPol {
    bool evaluada = false;
    double s_val = 0.0;
    double p_bil = 1.0;
    bool rech = false;
    std::size_t n_ev_ctr = 0;
    std::array<AcumEscalaPol, 3> por_esc{};
};

struct ResultadoSemilla {
    bool cal_ok = false;
    std::size_t escalas = 0;
    std::size_t t_max_ef = 0;
    std::size_t burn_in = 0;
    std::array<std::size_t, 3> n_reales_esc{};
    // Matriz 3x3 bajo P0 sobre intentos donde hay_propia (incluye extraible des in {2,4})
    std::array<std::array<std::size_t, 3>, 3> mat_esc_adm{};
    std::array<std::array<std::size_t, 3>, 3> mat_esc_propia{};
    double sum_x_incond = 0.0;
    std::size_t n_incond = 0;
    std::array<ResumenSemillaPol, N_POL> pols{};
};

template <class ObsLilT>
ResultadoSemilla correr_semilla_banco(
    bool es_arfima,
    double phi1,
    double d_arf,
    double mu_serie,
    std::size_t n_cal,
    std::size_t n_op,
    std::size_t b_max_nulo,
    std::uint64_t sc,
    std::uint64_t so,
    std::uint64_t sn) {

    ResultadoSemilla out{};
    const auto x_cal = es_arfima
        ? generar_arfima1(n_cal, phi1, d_arf, mu_serie, 1.0, sc)
        : generar_ar1_estacionario(n_cal, phi1, 1.0, sc);
    const auto x_op = es_arfima
        ? generar_arfima1(n_op, phi1, d_arf, mu_serie, 1.0, so)
        : generar_ar1_estacionario(n_op, phi1, 1.0, so);

    ConfiguracionObservador cd{};
    cd.integracion = Integracion::FLUJO;
    cd.agregacion = AgregacionEscala::SUMA;
    cd.arl0_pasos_nativos = 450.0;
    cd.alfa = ALFA_PROYECTO;
    cd.purga = false;
    cd.bins_fase = 1;
    cd.medir_tau50 = false;
    cd.escalas_max = 3;
    cd.semilla = sc;

    auto od = std::make_unique<ObsDelta>();
    od->calibrar(cd, x_cal);
    out.escalas = od->escalas();

    ConfiguracionLilyth cl{};
    cl.alfa = ALFA_PROYECTO;
    cl.representacion = ModoRepresentacionLilyth::CRUDA;
    cl.modo_nulo_lineal = es_arfima ? ModoNuloLineal::ARFIMA_WHITTLE : ModoNuloLineal::AR_LEVINSON_DURBIN;
    if (es_arfima) cl.d_whittle_cruda = d_arf;
    cl.modo_armonico = ModoArmonico::DESACTIVADO;
    cl.b_max_nulo = b_max_nulo;
    cl.hilos_nulo = 1;
    cl.b_mc_resumen = 1999;
    cl.semilla_nulo = sn;
    cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cl.forzar_via_b = true;
    cl.usar_enmienda_nivel2 = true;
    cl.usar_enmienda_nivel3 = true;

    auto lil = std::make_unique<ObsLilT>();
    if (!lil->calibrar(*od, cl, cd, x_cal)) return out;
    out.cal_ok = true;

    auto od_cal = std::make_unique<ObsDelta>(*od);
    auto lil_plantilla = std::make_unique<ObsLilT>(*lil);

    od->reiniciar_flujo();
    lil->reiniciar_flujo_conservando_purga();
    for (std::size_t t = 0; t < n_op; ++t) {
        const auto sd = od->paso(x_op[t], NaN, true, 0, 0);
        (void)lil->paso(sd, x_op[t], true, NaN, 0, 0);
    }

    const auto& fn = lil->ficha_nulo_cruda();
    const auto& m_lin = fn.modelo_lineal;
    const std::size_t t_max_ef = lil->configuracion().t_max_efectivo;
    out.t_max_ef = t_max_ef;
    out.burn_in = m_lin.burn_in;
    if (!m_lin.valido) return out;

    const double mu_ref = (es_arfima && std::isfinite(m_lin.mu_serie)) ? m_lin.mu_serie : m_lin.mu;
    const std::size_t m_min = fn.m_min;
    const auto& cerradas = lil->cerradas_historico();

    std::array<std::vector<std::vector<double>>, N_POL> estratos_pol;
    std::vector<double> pi_inv_cache;

    for (std::size_t idx = 0; idx < cerradas.size(); ++idx) {
        const auto& r = cerradas[idx];
        const auto des = r.cruda.desenlace;
        if (des != DesenlaceRanura::RESTITUIDA && des != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
        if (r.perfil_corrupto_z2_cruda || r.n_hist_prev_disp_cruda < m_lin.p) continue;
        if (r.escala_origen < 0 || r.escala_origen >= 3) continue;

        const std::span<const double> pf = lil->perfil_evento_n2(r, false);
        const std::span<const double> einy = lil->e_iny_evento_n2(r, false);
        if (pf.empty() || einy.size() != pf.size()) continue;
        bool einy_ok = true;
        for (double v : einy) if (!std::isfinite(v)) { einy_ok = false; break; }
        if (!einy_ok) continue;

        const std::size_t je = static_cast<std::size_t>(r.escala_origen);
        ++out.n_reales_esc[je];

        const double tau_real_cod = detalle_nulo::codificar_tau(r.cruda.tau_rest_consolidado, t_max_ef);
        const auto dim = dimensiones_simulacion_via_b(r.w_ancla, r.guarda_ancla, pf.size(), t_max_ef, Integracion::FLUJO);
        const std::size_t burn = m_lin.burn_in;
        const std::size_t n_sim = dim.n_sim;
        const std::size_t total = burn + n_sim;

        if (es_arfima && pi_inv_cache.size() < total) {
            pi_inv_cache.assign(total, 0.0);
            elw::coeficientes(-m_lin.d, total, pi_inv_cache.data());
        }

        const std::uint64_t sem_ev = sn ^ (0x9E3779B97F4A7C15ULL * (idx + 1ULL));
        std::mt19937_64 rng(detalle_nulo::splitmix64(sem_ev ^ 0x6C62272E07BB0142ULL));
        std::uniform_int_distribution<std::int64_t> dist_t(dim.t_iny_min, dim.t_iny_max);

        struct EstadoS3 {
            bool done = false;
            bool p_det_fijado = false;
            std::size_t tope = MAX_INTENTOS_VIA_B;
            double p_det_ini = 0.0;
            std::size_t intentos = 0;
            std::vector<double> taus;
            std::size_t n_esc_match = 0;
            double sum_x = 0.0;
            double sum_gd = 0.0;
        };
        std::array<EstadoS3, N_POL> s3{};
        std::vector<double> eps_nulo(total, 0.0);

        for (std::size_t intento = 1; intento <= MAX_INTENTOS_VIA_B; ++intento) {
            bool todos_listos = true;
            for (std::size_t p = 0; p < N_POL; ++p) {
                if (!s3[p].done) { todos_listos = false; break; }
            }
            if (todos_listos) break;

            const std::int64_t t0 = dist_t(rng);
            if (!m_lin.residuos.empty()) {
                std::uniform_int_distribution<std::size_t> dist_res(0, m_lin.residuos.size() - 1u);
                for (std::size_t k = 0; k < total; ++k) {
                    eps_nulo[k] = m_lin.residuos[dist_res(rng)] - m_lin.media_residuos;
                }
            } else {
                std::normal_distribution<double> dist_n(0.0, std::sqrt(m_lin.sigma2_innov));
                for (std::size_t k = 0; k < total; ++k) eps_nulo[k] = dist_n(rng);
            }

            const SalidaTrayectoria st_a = simular_trayectoria(
                *od_cal, *lil_plantilla, m_lin, einy, eps_nulo, pi_inv_cache, es_arfima,
                mu_ref, burn, n_sim, t0, r.signo, r.escala_origen, t_max_ef, false);

            out.sum_x_incond += st_a.x_pre_sig;
            ++out.n_incond;

            SalidaTrayectoria st_b{};
            const SalidaTrayectoria* ptr_b = &st_a;
            if (!s3[P1B_ESCALA_V4_ADM].done && st_a.vio_otra_escala_sig) {
                st_b = simular_trayectoria(
                    *od_cal, *lil_plantilla, m_lin, einy, eps_nulo, pi_inv_cache, es_arfima,
                    mu_ref, burn, n_sim, t0, r.signo, r.escala_origen, t_max_ef, true);
                ptr_b = &st_b;
            }

            // Registrar matriz de confusión de escala bajo P0 para todas las realizaciones activas
            {
                const auto sel_p0 = evaluar_politica_slot(st_a, t0, pf.size(), r.signo, r.escala_origen, r.g_disparo, P0_NIVEL3_BASE);
                if (sel_p0.hay_propia && sel_p0.escala >= 0 && sel_p0.escala < 3) {
                    ++out.mat_esc_propia[je][static_cast<std::size_t>(sel_p0.escala)];
                    if (sel_p0.admitida) {
                        ++out.mat_esc_adm[je][static_cast<std::size_t>(sel_p0.escala)];
                    }
                }
            }

            for (std::size_t p = 0; p < N_POL; ++p) {
                auto& ep = s3[p];
                if (ep.done) continue;
                ep.intentos = intento;
                const SalidaTrayectoria& st = (p == P1B_ESCALA_V4_ADM) ? *ptr_b : st_a;
                const auto sel = evaluar_politica_slot(
                    st, t0, pf.size(), r.signo, r.escala_origen, r.g_disparo, static_cast<PolIdx>(p));
                if (sel.admitida) {
                    ep.taus.push_back(sel.tau_cod);
                    if (sel.escala == r.escala_origen) ++ep.n_esc_match;
                    ep.sum_x += st.x_pre_sig;
                    if (std::isfinite(sel.g)) ep.sum_gd += (sel.g - r.g_disparo);
                }

                if (intento >= m_min && !ep.p_det_fijado) {
                    ep.p_det_fijado = true;
                    ep.p_det_ini = static_cast<double>(ep.taus.size()) / static_cast<double>(intento);
                    if (ep.taus.empty()) {
                        ep.done = true;
                    } else {
                        const double tc = std::ceil(static_cast<double>(m_min) / ep.p_det_ini);
                        ep.tope = (tc > static_cast<double>(MAX_INTENTOS_VIA_B))
                            ? MAX_INTENTOS_VIA_B
                            : static_cast<std::size_t>(tc);
                    }
                }
                if (ep.taus.size() >= m_min) {
                    if (!ep.p_det_fijado) {
                        ep.p_det_fijado = true;
                        ep.p_det_ini = 1.0;
                    }
                    ep.done = true;
                } else if (ep.p_det_fijado && intento >= ep.tope) {
                    ep.done = true;
                }
            }
        }

        for (std::size_t p = 0; p < N_POL; ++p) {
            const auto& ep = s3[p];
            auto& ac = out.pols[p].por_esc[je];
            ++ac.k_med;
            ac.intentos_tot += ep.intentos;
            ac.n_adm += ep.taus.size();
            ac.n_esc_match += ep.n_esc_match;
            ac.sum_x_pre_sig += ep.sum_x;
            ac.sum_g_diff += ep.sum_gd;
            if (ep.p_det_ini > 0.0) ++ac.k_pdet_pos;
            if (ep.taus.size() >= m_min) {
                ++ac.k_eval;
                const double u_ev = detalle_nulo::pit_punto_medio(tau_real_cod, ep.taus);
                ac.sum_u_ev += u_ev;
                ++ac.n_u_ev;

                std::vector<double> estrato;
                estrato.reserve(1u + ep.taus.size());
                estrato.push_back(tau_real_cod);
                estrato.insert(estrato.end(), ep.taus.begin(), ep.taus.end());
                estratos_pol[p].push_back(std::move(estrato));
            }
        }
    }

    for (std::size_t p = 0; p < N_POL; ++p) {
        const auto ctr = detalle_nulo::evaluar_contraste_intercambiable(
            estratos_pol[p], cl.b_mc_resumen, cl.semilla_nulo ^ 0x9876543210FEDCBAULL, cl.alfa / 4.0);
        out.pols[p].n_ev_ctr = ctr.n_eventos;
        if (ctr.valido && std::isfinite(ctr.estadistico_s) && std::isfinite(ctr.p_bilateral)) {
            out.pols[p].evaluada = true;
            out.pols[p].s_val = ctr.estadistico_s;
            out.pols[p].p_bil = ctr.p_bilateral;
            out.pols[p].rech = (ctr.p_bilateral < 0.05);
        }
    }

    return out;
}

void imprimir_matriz_confusion(
    const char* titulo,
    const std::array<std::array<std::size_t, 3>, 3>& mat) {
    std::printf("  %s (filas = escala_real 0,1,2 | cols = escala_nulo 0,1,2):\n", titulo);
    for (std::size_t jr = 0; jr < 3; ++jr) {
        const std::size_t tot = mat[jr][0] + mat[jr][1] + mat[jr][2];
        const double p_match = (tot > 0) ? (100.0 * static_cast<double>(mat[jr][jr]) / static_cast<double>(tot)) : 0.0;
        std::printf("    j_real=%zu : [%5zu, %5zu, %5zu]  (tot=%5zu, match=%6.2f%%)\n",
                    jr, mat[jr][0], mat[jr][1], mat[jr][2], tot, p_match);
    }
}

void imprimir_resumen_banco(
    const char* titulo,
    const std::vector<ResultadoSemilla>& ress) {

    std::size_t n_cal_ok = 0;
    std::size_t min_esc = 999;
    std::size_t min_tef = 999, max_tef = 0;
    std::array<std::size_t, 3> reales_esc{};
    std::array<std::array<std::size_t, 3>, 3> mat_propia{}, mat_adm{};
    double sum_x_incond = 0.0;
    std::size_t n_incond = 0;

    for (const auto& r : ress) {
        if (r.cal_ok) {
            ++n_cal_ok;
            min_esc = std::min(min_esc, r.escalas);
            min_tef = std::min(min_tef, r.t_max_ef);
            max_tef = std::max(max_tef, r.t_max_ef);
        }
        for (std::size_t j = 0; j < 3; ++j) {
            reales_esc[j] += r.n_reales_esc[j];
            for (std::size_t k = 0; k < 3; ++k) {
                mat_propia[j][k] += r.mat_esc_propia[j][k];
                mat_adm[j][k] += r.mat_esc_adm[j][k];
            }
        }
        sum_x_incond += r.sum_x_incond;
        n_incond += r.n_incond;
    }

    std::printf("\n========================================================================================================\n");
    std::printf("%s | semillas=%zu cal_ok=%zu min_escalas=%zu t_max_ef=[%zu,%zu]\n",
                titulo, ress.size(), n_cal_ok, min_esc, min_tef, max_tef);
    std::printf("  Eventos reales medibles por escala_real: j=0: %zu | j=1: %zu | j=2: %zu (tot=%zu)\n",
                reales_esc[0], reales_esc[1], reales_esc[2], reales_esc[0] + reales_esc[1] + reales_esc[2]);
    std::printf("  Incondicional E[(x_{t0-1}-mu)*signo] = %+.5f (n=%zu intentos P0)\n",
                n_incond > 0 ? sum_x_incond / static_cast<double>(n_incond) : 0.0, n_incond);
    imprimir_matriz_confusion("Matriz escala (propia antes de filtrar)", mat_propia);
    imprimir_matriz_confusion("Matriz escala (admitida P0)            ", mat_adm);

    std::printf("\n  %-18s | %4s | %8s %7s %7s | %4s %4s | %6s %6s %6s | %7s %7s %6s\n",
                "Politica", "eval", "S_medio", "ee_S", "z_S", "S<0", "rech",
                "p_det", "survS3", "cobS3", "Expre_s", "Eg_dif", "esc==");
    std::printf("  ------------------------------------------------------------------------------------------------------\n");

    for (std::size_t p = 0; p < N_POL; ++p) {
        std::vector<double> sv;
        std::size_t n_sneg = 0, n_rech = 0;
        AcumEscalaPol tot{};
        for (const auto& r : ress) {
            const auto& rp = r.pols[p];
            if (rp.evaluada) {
                sv.push_back(rp.s_val);
                if (rp.s_val < 0.0) ++n_sneg;
                if (rp.rech) ++n_rech;
            }
            for (std::size_t j = 0; j < 3; ++j) tot.sumar(rp.por_esc[j]);
        }
        double sm = 0.0, ee = 0.0;
        for (double v : sv) sm += v;
        if (!sv.empty()) sm /= static_cast<double>(sv.size());
        if (sv.size() >= 2) {
            double ss = 0.0;
            for (double v : sv) ss += (v - sm) * (v - sm);
            ee = std::sqrt(ss / static_cast<double>(sv.size() * (sv.size() - 1u)));
        }
        const double zs = (ee > 0.0) ? (sm / ee) : 0.0;
        const double pdet = tot.intentos_tot > 0 ? static_cast<double>(tot.n_adm) / static_cast<double>(tot.intentos_tot) : 0.0;
        const double surv = tot.k_med > 0 ? static_cast<double>(tot.k_pdet_pos) / static_cast<double>(tot.k_med) : 0.0;
        const double cob = tot.k_med > 0 ? static_cast<double>(tot.k_eval) / static_cast<double>(tot.k_med) : 0.0;
        const double xpre = tot.n_adm > 0 ? tot.sum_x_pre_sig / static_cast<double>(tot.n_adm) : 0.0;
        const double gdif = tot.n_adm > 0 ? tot.sum_g_diff / static_cast<double>(tot.n_adm) : 0.0;
        const double emat = tot.n_adm > 0 ? 100.0 * static_cast<double>(tot.n_esc_match) / static_cast<double>(tot.n_adm) : 0.0;

        std::printf("  %-18s | %2zu/%-2zu| %+8.5f %7.5f %+7.3f | %2zu/%-2zu %2zu  | %6.4f %6.4f %6.4f | %+7.4f %+7.3f %5.1f%%\n",
                    NOMBRE_POL[p], sv.size(), ress.size(), sm, ee, zs, n_sneg, sv.size(), n_rech,
                    pdet, surv, cob, xpre, gdif, emat);
    }

    std::printf("\n  Desglose por escala_real (j = 0, 1, 2):\n");
    for (std::size_t j = 0; j < 3; ++j) {
        std::printf("   [j_real=%zu]\n", j);
        for (std::size_t p = 0; p < N_POL; ++p) {
            AcumEscalaPol ej{};
            for (const auto& r : ress) ej.sumar(r.pols[p].por_esc[j]);
            const double pdet = ej.intentos_tot > 0 ? static_cast<double>(ej.n_adm) / static_cast<double>(ej.intentos_tot) : 0.0;
            const double surv = ej.k_med > 0 ? static_cast<double>(ej.k_pdet_pos) / static_cast<double>(ej.k_med) : 0.0;
            const double cob = ej.k_med > 0 ? static_cast<double>(ej.k_eval) / static_cast<double>(ej.k_med) : 0.0;
            const double s_ev = ej.n_u_ev > 0 ? (ej.sum_u_ev / static_cast<double>(ej.n_u_ev) - 0.5) : NaN;
            const double xpre = ej.n_adm > 0 ? ej.sum_x_pre_sig / static_cast<double>(ej.n_adm) : NaN;
            const double gdif = ej.n_adm > 0 ? ej.sum_g_diff / static_cast<double>(ej.n_adm) : NaN;
            const double emat = ej.n_adm > 0 ? 100.0 * static_cast<double>(ej.n_esc_match) / static_cast<double>(ej.n_adm) : 0.0;
            std::printf("    %-18s : kmed=%4zu kevl=%4zu | pdet=%.4f surv=%.4f cob=%.4f | S_ev=%+7.4f Ex=%+7.4f Eg=%+7.3f esc==%5.1f%%\n",
                        NOMBRE_POL[p], ej.k_med, ej.k_eval, pdet, surv, cob, s_ev, xpre, gdif, emat);
        }
    }
}

void ejecutar_diag_gruesas_preciso(std::size_t hilos) {
    constexpr std::size_t N_CAL = 12000;
    constexpr std::size_t N_OP = 120000;
    constexpr double PHI = 0.80;
    constexpr std::size_t MAX_EV_ESC = 20;
    constexpr std::size_t N_MC = 200;

    const std::uint64_t sc = detalle_nulo::splitmix64(SEMILLA_SONDA_NIVEL4 + 90001ULL);
    const std::uint64_t so = detalle_nulo::splitmix64(SEMILLA_SONDA_NIVEL4 + 90002ULL);
    const std::uint64_t sn = detalle_nulo::splitmix64(SEMILLA_SONDA_NIVEL4 + 90003ULL);

    const auto x_cal = generar_ar1_estacionario(N_CAL, PHI, 1.0, sc);
    const auto x_op = generar_ar1_estacionario(N_OP, PHI, 1.0, so);

    ConfiguracionObservador cd{};
    cd.integracion = Integracion::FLUJO;
    cd.agregacion = AgregacionEscala::SUMA;
    cd.arl0_pasos_nativos = 450.0;
    cd.alfa = ALFA_PROYECTO;
    cd.purga = false;
    cd.bins_fase = 1;
    cd.medir_tau50 = false;
    cd.escalas_max = 3;
    cd.semilla = sc;

    auto od = std::make_unique<ObsDelta>();
    od->calibrar(cd, x_cal);

    ConfiguracionLilyth cl{};
    cl.alfa = ALFA_PROYECTO;
    cl.representacion = ModoRepresentacionLilyth::CRUDA;
    cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
    cl.modo_armonico = ModoArmonico::DESACTIVADO;
    cl.b_max_nulo = 8;
    cl.hilos_nulo = 1;
    cl.semilla_nulo = sn;
    cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cl.forzar_via_b = true;
    cl.usar_enmienda_nivel2 = true;
    cl.usar_enmienda_nivel3 = true;

    auto lil = std::make_unique<ObsLil64>();
    const bool cal_ok = lil->calibrar(*od, cl, cd, x_cal);
    auto od_cal = std::make_unique<ObsDelta>(*od);
    auto lil_plantilla = std::make_unique<ObsLil64>(*lil);

    od->reiniciar_flujo();
    lil->reiniciar_flujo_conservando_purga();
    for (std::size_t t = 0; t < N_OP; ++t) {
        const auto sd = od->paso(x_op[t], NaN, true, 0, 0);
        (void)lil->paso(sd, x_op[t], true, NaN, 0, 0);
    }

    const auto& fn = lil->ficha_nulo_cruda();
    const auto& m_lin = fn.modelo_lineal;
    const std::size_t t_max_ef = lil->configuracion().t_max_efectivo;

    struct EvSel {
        std::size_t idx = 0;
        int escala = 0;
        int signo = 0;
        double g = 0.0;
        double y_ancla = 0.0;
        double tau_real_cod = 0.0;
        std::size_t w_ancla = 0;
        std::size_t guarda = 0;
        std::vector<double> einy;
    };

    std::array<std::size_t, 3> total_op_esc{};
    std::vector<EvSel> sels;
    std::array<std::size_t, 3> cuenta_sel{};
    const auto& cerradas = lil->cerradas_historico();
    for (std::size_t idx = 0; idx < cerradas.size(); ++idx) {
        const auto& r = cerradas[idx];
        const auto des = r.cruda.desenlace;
        if (des != DesenlaceRanura::RESTITUIDA && des != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
        if (r.perfil_corrupto_z2_cruda || r.n_hist_prev_disp_cruda < m_lin.p) continue;
        if (r.escala_origen < 0 || r.escala_origen >= 3) continue;
        const auto pf = lil->perfil_evento_n2(r, false);
        const auto ei = lil->e_iny_evento_n2(r, false);
        if (pf.empty() || ei.size() != pf.size()) continue;
        const std::size_t j = static_cast<std::size_t>(r.escala_origen);
        ++total_op_esc[j];
        if (cuenta_sel[j] < MAX_EV_ESC) {
            ++cuenta_sel[j];
            EvSel es{};
            es.idx = idx;
            es.escala = r.escala_origen;
            es.signo = r.signo;
            es.g = r.g_disparo;
            es.y_ancla = r.cruda.y_ancla;
            es.tau_real_cod = detalle_nulo::codificar_tau(r.cruda.tau_rest_consolidado, t_max_ef);
            es.w_ancla = r.w_ancla;
            es.guarda = r.guarda_ancla;
            es.einy.assign(ei.begin(), ei.end());
            sels.push_back(std::move(es));
        }
    }

    struct ResEv {
        int esc = 0;
        double y_ancla_real = 0.0;
        std::array<std::size_t, 3> mat_adm{};
        double sum_x_incond = 0.0;
        std::array<std::size_t, N_POL> n_adm{};
        std::array<std::size_t, N_POL> n_esc_match{};
        std::array<double, N_POL> u_ev{};
        std::array<bool, N_POL> tiene_u{};
        std::array<double, N_POL> sum_x{};
        std::array<double, N_POL> sum_g_dif{};
        std::array<double, N_POL> sum_ya{};
        std::array<double, N_POL> sum_ya2{};
    };

    std::vector<ResEv> res_evs(sels.size());
    std::atomic<std::size_t> next_ev{0};
    const std::size_t hw = std::clamp<std::size_t>(hilos, 1u, std::max<std::size_t>(1u, sels.size()));
    std::vector<std::thread> ths;
    ths.reserve(hw);
    for (std::size_t w = 0; w < hw; ++w) {
        ths.emplace_back([&]() {
            while (true) {
                const std::size_t k = next_ev.fetch_add(1u, std::memory_order_relaxed);
                if (k >= sels.size()) break;
                const auto& ev = sels[k];
                auto& re = res_evs[k];
                re.esc = ev.escala;
                re.y_ancla_real = ev.y_ancla;

                const auto dim = dimensiones_simulacion_via_b(ev.w_ancla, ev.guarda, ev.einy.size(), t_max_ef, Integracion::FLUJO);
                const std::size_t burn = m_lin.burn_in;
                const std::size_t n_sim = dim.n_sim;
                const std::size_t total = burn + n_sim;
                std::vector<double> pi_vacio;
                std::vector<double> eps_nulo(total, 0.0);

                std::mt19937_64 rng(detalle_nulo::splitmix64(sn ^ (0x9E3779B97F4A7C15ULL * (ev.idx + 1ULL))));
                std::uniform_int_distribution<std::int64_t> dist_t(dim.t_iny_min, dim.t_iny_max);
                std::array<std::vector<double>, N_POL> taus_p;

                for (std::size_t mc = 0; mc < N_MC; ++mc) {
                    const std::int64_t t0 = dist_t(rng);
                    if (!m_lin.residuos.empty()) {
                        std::uniform_int_distribution<std::size_t> dr(0, m_lin.residuos.size() - 1u);
                        for (std::size_t q = 0; q < total; ++q) eps_nulo[q] = m_lin.residuos[dr(rng)] - m_lin.media_residuos;
                    } else {
                        std::normal_distribution<double> dn(0.0, std::sqrt(m_lin.sigma2_innov));
                        for (std::size_t q = 0; q < total; ++q) eps_nulo[q] = dn(rng);
                    }

                    const auto st_a = simular_trayectoria(
                        *od_cal, *lil_plantilla, m_lin, ev.einy, eps_nulo, pi_vacio, false,
                        m_lin.mu, burn, n_sim, t0, ev.signo, ev.escala, t_max_ef, false);
                    re.sum_x_incond += st_a.x_pre_sig;

                    SalidaTrayectoria st_b{};
                    const SalidaTrayectoria* ptr_b = &st_a;
                    if (st_a.vio_otra_escala_sig) {
                        st_b = simular_trayectoria(
                            *od_cal, *lil_plantilla, m_lin, ev.einy, eps_nulo, pi_vacio, false,
                            m_lin.mu, burn, n_sim, t0, ev.signo, ev.escala, t_max_ef, true);
                        ptr_b = &st_b;
                    }

                    for (std::size_t p = 0; p < N_POL; ++p) {
                        const SalidaTrayectoria& st = (p == P1B_ESCALA_V4_ADM) ? *ptr_b : st_a;
                        const auto sel = evaluar_politica_slot(
                            st, t0, ev.einy.size(), ev.signo, ev.escala, ev.g, static_cast<PolIdx>(p));
                        if (sel.admitida) {
                            ++re.n_adm[p];
                            if (p == P0_NIVEL3_BASE && sel.escala >= 0 && sel.escala < 3) {
                                ++re.mat_adm[static_cast<std::size_t>(sel.escala)];
                            }
                            if (sel.escala == ev.escala) ++re.n_esc_match[p];
                            re.sum_x[p] += st.x_pre_sig;
                            if (std::isfinite(sel.g)) re.sum_g_dif[p] += (sel.g - ev.g);
                            if (std::isfinite(sel.y_ancla)) {
                                re.sum_ya[p] += sel.y_ancla;
                                re.sum_ya2[p] += sel.y_ancla * sel.y_ancla;
                            }
                            taus_p[p].push_back(sel.tau_cod);
                        }
                    }
                }

                for (std::size_t p = 0; p < N_POL; ++p) {
                    if (!taus_p[p].empty()) {
                        re.tiene_u[p] = true;
                        re.u_ev[p] = detalle_nulo::pit_punto_medio(ev.tau_real_cod, taus_p[p]);
                    }
                }
            }
        });
    }
    for (auto& th : ths) th.join();

    std::printf("\n========================================================================================================\n");
    std::printf("BANCO 3: M3_DIAG_GRUESAS (AR(1) phi=0.80, n_cal=%zu, n_op=%zu, cal_ok=%d, escalas=%zu, t_max_ef=%zu, N_mc=%zu)\n",
                N_CAL, N_OP, cal_ok ? 1 : 0, od->escalas(), t_max_ef, N_MC);
    std::printf("  Eventos reales medibles totales: j=0: %zu (sel=%zu) | j=1: %zu (sel=%zu) | j=2: %zu (sel=%zu)\n",
                total_op_esc[0], cuenta_sel[0], total_op_esc[1], cuenta_sel[1], total_op_esc[2], cuenta_sel[2]);

    for (std::size_t j = 0; j < 3; ++j) {
        std::array<std::size_t, 3> mat_j{};
        std::vector<double> ya_reales;
        double sx_inc = 0.0;
        std::size_t n_ev_j = 0;
        for (const auto& re : res_evs) {
            if (static_cast<std::size_t>(re.esc) != j) continue;
            ++n_ev_j;
            for (std::size_t q = 0; q < 3; ++q) mat_j[q] += re.mat_adm[q];
            if (std::isfinite(re.y_ancla_real)) ya_reales.push_back(re.y_ancla_real);
            sx_inc += re.sum_x_incond;
        }
        double sd_ya_real = NaN;
        if (ya_reales.size() >= 2) {
            double m = 0.0;
            for (double v : ya_reales) m += v;
            m /= static_cast<double>(ya_reales.size());
            double ss = 0.0;
            for (double v : ya_reales) ss += (v - m) * (v - m);
            sd_ya_real = std::sqrt(ss / static_cast<double>(ya_reales.size() - 1u));
        }
        const std::size_t mtot = mat_j[0] + mat_j[1] + mat_j[2];
        const double ex_inc = (n_ev_j > 0) ? (sx_inc / static_cast<double>(n_ev_j * N_MC)) : 0.0;
        std::printf("\n  --- j_real=%zu (eventos=%zu, intents=%zu, E_incond[x*s]=%+.4f, SD(y_anc_real)=%.4f) ---\n",
                    j, n_ev_j, n_ev_j * N_MC, ex_inc, sd_ya_real);
        std::printf("  P0 P(escala_nulo | j_real=%zu): j'=0:%5zu (%.1f%%) | j'=1:%5zu (%.1f%%) | j'=2:%5zu (%.1f%%)\n",
                    j,
                    mat_j[0], mtot > 0 ? 100.0 * static_cast<double>(mat_j[0]) / static_cast<double>(mtot) : 0.0,
                    mat_j[1], mtot > 0 ? 100.0 * static_cast<double>(mat_j[1]) / static_cast<double>(mtot) : 0.0,
                    mat_j[2], mtot > 0 ? 100.0 * static_cast<double>(mat_j[2]) / static_cast<double>(mtot) : 0.0);

        for (std::size_t p = 0; p < N_POL; ++p) {
            std::size_t na = 0, nm = 0, nu = 0;
            double su = 0.0, sx = 0.0, sgd = 0.0, sya = 0.0, sya2 = 0.0;
            for (const auto& re : res_evs) {
                if (static_cast<std::size_t>(re.esc) != j) continue;
                na += re.n_adm[p];
                nm += re.n_esc_match[p];
                sx += re.sum_x[p];
                sgd += re.sum_g_dif[p];
                sya += re.sum_ya[p];
                sya2 += re.sum_ya2[p];
                if (re.tiene_u[p]) { su += re.u_ev[p]; ++nu; }
            }
            const std::size_t n_tot = n_ev_j * N_MC;
            const double pdet = n_tot > 0 ? static_cast<double>(na) / static_cast<double>(n_tot) : 0.0;
            const double s_med = nu > 0 ? (su / static_cast<double>(nu) - 0.5) : NaN;
            const double ex = na > 0 ? sx / static_cast<double>(na) : NaN;
            const double eg = na > 0 ? sgd / static_cast<double>(na) : NaN;
            const double sd_ya_nul = (na >= 2) ? std::sqrt(std::max(0.0, (sya2 - sya * sya / static_cast<double>(na)) / static_cast<double>(na - 1u))) : NaN;
            const double emat = na > 0 ? 100.0 * static_cast<double>(nm) / static_cast<double>(na) : 0.0;
            std::printf("    %-18s : pdet=%.4f (%4zu/%4zu) | S=%+7.4f (nu=%2zu) | Ex=%+7.4f | Eg=%+7.3f | SD_anc=%.4f | esc==%5.1f%%\n",
                        NOMBRE_POL[p], pdet, na, n_tot, s_med, nu, ex, eg, sd_ya_nul, emat);
        }
    }
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
                if (c < '0' || c > '9') { ok = false; break; }
                parsed = parsed * 10u + static_cast<std::size_t>(c - '0');
            }
            if (ok && parsed > 0) hilos = parsed;
        }
    }
    return std::max<std::size_t>(1u, hilos);
}

}  // namespace

int main(int argc, char** argv) {
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t hilos = parsear_hilos(argc, argv);

    std::printf("========================================================================================================\n");
    std::printf("FASE 4.0 -- SONDA EMPIRICA MULTIESCALA DE NIVEL 4 (escala_origen y g_disparo >= g_real en Via (b))\n");
    std::printf("  SEMILLA_SONDA_NIVEL4 = %lluULL | hilos = %zu\n",
                static_cast<unsigned long long>(SEMILLA_SONDA_NIVEL4), hilos);
    std::printf("========================================================================================================\n");
    std::fflush(stdout);

    // Banco 1: M1_AR1_3ESC (AR(1) phi=0.80, 24 semillas, n_cal=12000, n_op=12000, b_max=8)
    {
        constexpr std::size_t N_SEM = 24;
        std::vector<ResultadoSemilla> res(N_SEM);
        std::atomic<std::size_t> sig{0};
        const std::size_t nh = std::min(hilos, N_SEM);
        std::vector<std::thread> ths;
        for (std::size_t w = 0; w < nh; ++w) {
            ths.emplace_back([&]() {
                while (true) {
                    const std::size_t s = sig.fetch_add(1u, std::memory_order_relaxed);
                    if (s >= N_SEM) break;
                    const std::uint64_t sc = detalle_nulo::splitmix64(SEMILLA_SONDA_NIVEL4 ^ (0x9E3779B97F4A7C15ULL * (s + 1ULL)));
                    const std::uint64_t so = detalle_nulo::splitmix64(SEMILLA_SONDA_NIVEL4 ^ (0xBF58476D1CE4E5B9ULL * (s + 101ULL)));
                    const std::uint64_t sn = detalle_nulo::splitmix64(SEMILLA_SONDA_NIVEL4 ^ (0x94D049BB133111EBULL * (s + 211ULL)));
                    res[s] = correr_semilla_banco<ObsLil64>(false, 0.80, 0.0, 0.0, 12000, 12000, 8, sc, so, sn);
                }
            });
        }
        for (auto& th : ths) th.join();
        imprimir_resumen_banco("BANCO 1: M1_AR1_3ESC (AR(1) phi=0.80, 24 sem, n_cal=12000, n_op=12000, b_max=8)", res);
        std::fflush(stdout);
    }

    // Banco 2: M2_ARFIMA_3ESC (ARFIMA(1,d=+0.35,0) phi1=0.40 mu=50.0, 16 sem, n_cal=12000, n_op=8000, b_max=4)
    {
        constexpr std::size_t N_SEM = 16;
        std::vector<ResultadoSemilla> res(N_SEM);
        std::atomic<std::size_t> sig{0};
        const std::size_t nh = std::min(hilos, N_SEM);
        std::vector<std::thread> ths;
        for (std::size_t w = 0; w < nh; ++w) {
            ths.emplace_back([&]() {
                while (true) {
                    const std::size_t s = sig.fetch_add(1u, std::memory_order_relaxed);
                    if (s >= N_SEM) break;
                    const std::uint64_t sc = detalle_nulo::splitmix64(SEMILLA_SONDA_NIVEL4 ^ (0x517CC1B727220A95ULL * (s + 1001ULL)));
                    const std::uint64_t so = detalle_nulo::splitmix64(SEMILLA_SONDA_NIVEL4 ^ (0xD6E8FEB86659FD93ULL * (s + 1307ULL)));
                    const std::uint64_t sn = detalle_nulo::splitmix64(SEMILLA_SONDA_NIVEL4 ^ (0xA0761D6478BD642FULL * (s + 1409ULL)));
                    res[s] = correr_semilla_banco<ObsLil128>(true, 0.40, 0.35, 50.0, 12000, 8000, 4, sc, so, sn);
                }
            });
        }
        for (auto& th : ths) th.join();
        imprimir_resumen_banco("BANCO 2: M2_ARFIMA_3ESC (ARFIMA d=+0.35 phi1=0.40 mu=50, 16 sem, n_cal=12000, n_op=8000, b_max=4)", res);
        std::fflush(stdout);
    }

    // Banco 3: M3_DIAG_GRUESAS (N_mc = 200 por evento en j = 0, 1, 2)
    ejecutar_diag_gruesas_preciso(hilos);

    const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("\n========================================================================================================\n");
    std::printf("SONDA NIVEL 4 COMPLETADA (tiempo=%.2f s)\n", dt);
    std::printf("========================================================================================================\n");
    std::fflush(stdout);
    return 0;
}
