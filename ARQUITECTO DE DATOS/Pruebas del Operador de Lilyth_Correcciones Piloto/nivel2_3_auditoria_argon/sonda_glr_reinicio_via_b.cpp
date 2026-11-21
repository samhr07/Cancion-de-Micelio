// ============================================================================
// sonda_glr_reinicio_via_b.cpp
//
// Sonda reproducible para responder al Dictamen de Argon (§2.1 y §3.1.C):
//   1. Mide el efecto real de DetectorGLR::reiniciar() (canal_nivel.hpp:54-68)
//      cuando Argon §3.1.C suprime aperturas de ObservadorLilyth en [t0, t0+l-2]
//      sin detener ni preservar a ObservadorSaltos (d_local).
//   2. Cuantifica el sesgo de condicionamiento en x_{t0-1} y el impacto de retirar
//      hay_previa (nulo_lilyth.hpp:1079-1087).
//   3. Contrasta 5 politicas de via (b) en banco multi-semilla R1 (AR(1), phi=0.80):
//        V0_BASE_ACTUAL            : [t0, t0+l], hay_previa=true (codigo actual)
//        V1_BASE_SIN_HAY_PREVIA    : [t0, t0+l], hay_previa=false
//        V2_ARGON_3_1_C_LITERAL    : suprime lil en [t0, t0+l-2], d_local corre intacto, exige t_fin
//        V3_OPCION_A_EXACTA_TFIN   : (d_local, lil_local) intactos, acepta solo t_ap == t0+l-1 (sin hay_previa)
//        V4_OPCION_B_PARADA_TEMP   : (d_local, lil_local) intactos, corta e_iny tras apertura temprana (sin hay_previa)
//
// Compilacion (Windows):
//   compilar.bat nivel2_3_auditoria_argon\sonda_glr_reinicio_via_b.cpp sonda_glr_reinicio_via_b.exe
// Compilacion (Linux / POSIX):
//   g++ -std=c++20 -O2 -pthread -Iinc nivel2_3_auditoria_argon/sonda_glr_reinicio_via_b.cpp -o sonda_glr_reinicio_via_b
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
#include <thread>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;

inline constexpr std::size_t J_CAP = 24;
inline constexpr std::size_t L_GLR = 64;
inline constexpr std::size_t N_REG_S = 120;
inline constexpr std::size_t T_MAX_CAP = 64;

using ObsDelta = ObservadorSaltos<J_CAP, L_GLR, N_REG_S>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, T_MAX_CAP>;
using RanuraLil = typename ObsLil::Ranura;

std::vector<double> generar_ar1_estacionaria(std::size_t n, double phi, double sd_eps, std::mt19937_64& rng) {
    std::vector<double> x(n, 0.0);
    if (n == 0) return x;
    std::normal_distribution<double> n01(0.0, 1.0);
    const double sd_est = sd_eps / std::sqrt(std::max(1e-12, 1.0 - phi * phi));
    double prev = sd_est * n01(rng);
    x[0] = prev;
    for (std::size_t t = 1; t < n; ++t) {
        prev = phi * prev + sd_eps * n01(rng);
        x[t] = prev;
    }
    return x;
}

struct EventoCapturado {
    std::int64_t t_cambio = -1;
    std::int64_t t_apertura = -1;
    std::size_t l = 0;
    int escala = 0;
    int signo = 0;
    std::size_t w_ancla = 0;
    std::size_t guarda = 0;
    double tau_real_cod = NaN;
    std::vector<double> e_iny;
    std::vector<double> perfil;
};

struct CalibracionSonda {
    ConfiguracionObservador cd{};
    ConfiguracionLilyth cl{};
    std::unique_ptr<ObsDelta> od_base;
    std::unique_ptr<ObsLil> lil_base;
    ModeloLinealNulo m_lin{};
    std::size_t t_max_ef = T_MAX_CAP;
};

CalibracionSonda calibrar_sonda_ar1(std::uint64_t sem, std::size_t n_cal, double phi) {
    CalibracionSonda cs{};
    cs.cd.integracion = Integracion::FLUJO;
    cs.cd.agregacion = AgregacionEscala::SUMA;
    cs.cd.arl0_pasos_nativos = 450.0;
    cs.cd.alfa = ALFA_PROYECTO;
    cs.cd.purga = false;
    cs.cd.bins_fase = 1;
    cs.cd.medir_tau50 = false;
    cs.cd.semilla = sem;
    cs.cd.escalas_max = 3;

    cs.cl.alfa = ALFA_PROYECTO;
    cs.cl.representacion = ModoRepresentacionLilyth::CRUDA;
    cs.cl.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
    cs.cl.modo_armonico = ModoArmonico:: DESACTIVADO;
    cs.cl.b_max_nulo = 16;
    cs.cl.hilos_nulo = 1;
    cs.cl.semilla_nulo = sem ^ 0x9E3779B97F4A7C15ULL;
    cs.cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;

    std::mt19937_64 rng_cal(detalle_nulo::splitmix64(sem ^ 0x517CC1B727220A95ULL));
    const auto x_cal = generar_ar1_estacionaria(n_cal, phi, 1.0, rng_cal);

    cs.od_base = std::make_unique<ObsDelta>();
    cs.od_base->calibrar(cs.cd, x_cal);

    cs.lil_base = std::make_unique<ObsLil>();
    cs.lil_base->calibrar(*cs.od_base, cs.cl, cs.cd, x_cal);

    cs.m_lin = cs.lil_base->ficha_nulo_cruda().modelo_lineal;
    cs.t_max_ef = cs.lil_base->configuracion().t_max_efectivo;
    return cs;
}

std::vector<EventoCapturado> extraer_eventos_reales(
    const CalibracionSonda& cs, std::uint64_t sem_op, std::size_t n_op, double phi) {
    std::mt19937_64 rng_op(detalle_nulo::splitmix64(sem_op));
    const auto x_op = generar_ar1_estacionaria(n_op, phi, 1.0, rng_op);

    auto od = std::make_unique<ObsDelta>(*cs.od_base);
    od->reiniciar_flujo();
    auto lil = std::make_unique<ObsLil>(*cs.lil_base);
    lil->reiniciar_flujo_conservando_purga();

    for (std::size_t t = 0; t < n_op; ++t) {
        const auto sd = od->paso(x_op[t], NaN, true, 0, 0);
        (void)lil->paso(sd, x_op[t], true, NaN, 0, 0);
    }

    std::vector<EventoCapturado> out;
    for (const auto& r : lil->cerradas_historico()) {
        const auto d = r.cruda.desenlace;
        if (d != DesenlaceRanura::RESTITUIDA && d != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
        if (r.n_perfil_cruda == 0 || !cs.m_lin.valido) continue;
        EventoCapturado ev{};
        ev.t_cambio = r.t_cambio;
        ev.t_apertura = r.t_apertura;
        ev.l = r.n_perfil_cruda;
        ev.escala = r.escala_origen;
        ev.signo = r.signo;
        ev.w_ancla = r.w_ancla;
        ev.guarda = r.guarda_ancla;
        ev.tau_real_cod = detalle_nulo::codificar_tau(r.cruda.tau_rest_consolidado, cs.t_max_ef);
        ev.perfil.assign(r.perfil_cruda_est.begin(), r.perfil_cruda_est.begin() + static_cast<std::ptrdiff_t>(r.n_perfil_cruda));
        ev.e_iny.assign(ev.l, 0.0);
        innovaciones(cs.m_lin,
                     r.n_hist_prev_cruda > 0 ? r.hist_prev_cruda_est.data() : nullptr,
                     r.n_hist_prev_cruda,
                     ev.perfil.data(), ev.l, ev.e_iny.data());
        bool fin = true;
        for (double v : ev.e_iny) fin = fin && std::isfinite(v);
        if (fin) out.push_back(std::move(ev));
    }
    return out;
}

struct CierreSim {
    std::int64_t t_ap = -1;
    int signo = 0;
    int escala = -1;
    std::uint8_t des = 0;
    double tau_cod = NaN;
};

struct CorreccionPaso {
    std::vector<CierreSim> cierres;
    std::vector<std::int64_t> aperturas_todas;
    bool glr_prem_cualq = false;
    bool glr_prem_misma = false;
    bool glr_en_tfin_misma = false;
    bool glr_en_tfin_cualq = false;
    // Trayectoria Argon 3.1.C literal: d_local corre intacto, lil_local ignora aperturas < t_fin,
    // y evalua apertura en t_fin si d_local alarma en escala_real con signo_real
    bool arg_admitida = false;
    double arg_tau_cod = NaN;
};

CorreccionPaso ejecutar_sim_unificada(
    const CalibracionSonda& cs,
    const EventoCapturado& ev,
    const std::vector<double>& eps_tot,
    std::size_t burn,
    std::size_t n_sim,
    std::int64_t t0,
    bool modo_parada_temprana,
    std::vector<double>& x_sim_out) {

    CorreccionPaso cp{};
    const std::size_t total = burn + n_sim;
    const std::size_t l = ev.l;
    const std::int64_t t_fin = t0 + static_cast<std::int64_t>(l) - 1;

    x_sim_out.assign(n_sim, 0.0);
    std::vector<double> u(total, 0.0);

    // Precalentamiento hasta t0 - 1
    const std::size_t idx_t0 = burn + static_cast<std::size_t>(t0);
    for (std::size_t idx = 0; idx < total; ++idx) {
        if (!modo_parada_temprana && idx >= idx_t0 && idx < idx_t0 + l) {
            double acc = ev.e_iny[idx - idx_t0];
            const std::size_t imax = std::min(cs.m_lin.p, idx);
            for (std::size_t i = 1; i <= imax; ++i) acc += cs.m_lin.phi[i - 1] * u[idx - i];
            u[idx] = acc;
        } else {
            double acc = eps_tot[idx];
            const std::size_t imax = std::min(cs.m_lin.p, idx);
            for (std::size_t i = 1; i <= imax; ++i) acc += cs.m_lin.phi[i - 1] * u[idx - i];
            u[idx] = acc;
        }
    }

    auto od = std::make_unique<ObsDelta>(*cs.od_base);
    od->reiniciar_flujo();
    auto lil = std::make_unique<ObsLil>(*cs.lil_base);
    lil->reiniciar_flujo_conservando_purga();

    // Tambien mantenemos un espejo paralelo de lil para medir Argon 3.1.C literal sobre el MISMO flujo x_sim
    std::unique_ptr<ObsLil> lil_arg;
    if (!modo_parada_temprana) {
        lil_arg = std::make_unique<ObsLil>(*cs.lil_base);
        lil_arg->reiniciar_flujo_conservando_purga();
    }

    bool parada_activada = false;
    std::vector<ObsDelta::Salida> salidas_d;
    if (!modo_parada_temprana) salidas_d.resize(n_sim);

    for (std::size_t t = 0; t < n_sim; ++t) {
        const std::size_t idx = burn + t;
        if (modo_parada_temprana && !parada_activada && static_cast<std::int64_t>(t) >= t0 && static_cast<std::int64_t>(t) <= t_fin) {
            double acc = ev.e_iny[t - static_cast<std::size_t>(t0)];
            const std::size_t imax = std::min(cs.m_lin.p, idx);
            for (std::size_t i = 1; i <= imax; ++i) acc += cs.m_lin.phi[i - 1] * u[idx - i];
            u[idx] = acc;
        } else if (modo_parada_temprana) {
            double acc = eps_tot[idx];
            const std::size_t imax = std::min(cs.m_lin.p, idx);
            for (std::size_t i = 1; i <= imax; ++i) acc += cs.m_lin.phi[i - 1] * u[idx - i];
            u[idx] = acc;
        }
        const double xt = cs.m_lin.mu + u[idx];
        x_sim_out[t] = xt;

        const auto sd = od->paso(xt, NaN, true, 0, 0);
        if (!modo_parada_temprana) {
            salidas_d[t] = sd;
            const std::int64_t st = static_cast<std::int64_t>(t);
            if (st >= t0 && st < t_fin) {
                for (std::size_t j = 0; j < sd.niveles; ++j) {
                    if (sd.nivel[j].alarma.hay) {
                        cp.glr_prem_cualq = true;
                        if (static_cast<int>(j) == ev.escala && sd.nivel[j].alarma.signo == ev.signo) {
                            cp.glr_prem_misma = true;
                        }
                    }
                }
            } else if (st == t_fin) {
                for (std::size_t j = 0; j < sd.niveles; ++j) {
                    if (sd.nivel[j].alarma.hay && sd.nivel[j].alarma.signo == ev.signo) {
                        cp.glr_en_tfin_cualq = true;
                        if (static_cast<int>(j) == ev.escala) cp.glr_en_tfin_misma = true;
                    }
                }
            }
        }

        const auto o = lil->paso(sd, xt, true, NaN, 0, 0);
        if (o.abrio_ranura) {
            cp.aperturas_todas.push_back(o.t);
            if (modo_parada_temprana && static_cast<std::int64_t>(t) >= t0 && static_cast<std::int64_t>(t) <= t_fin) {
                // Comprobar si la nueva ranura tiene el signo_real
                for (const auto& r : lil->ranuras()) {
                    if (r.ocupada && r.t_apertura == o.t && r.signo == ev.signo) {
                        parada_activada = true;
                        break;
                    }
                }
            }
        }
        for (std::size_t c = 0; c < o.ranuras_cerradas; ++c) {
            const auto& rc = lil->cerrada(c);
            CierreSim ci{};
            ci.t_ap = rc.t_apertura;
            ci.signo = rc.signo;
            ci.escala = rc.escala_origen;
            ci.des = static_cast<std::uint8_t>(rc.cruda.desenlace);
            if (ci.des == 2 || ci.des == 4) {
                ci.tau_cod = detalle_nulo::codificar_tau(rc.cruda.tau_rest_consolidado, cs.t_max_ef);
            }
            cp.cierres.push_back(ci);
        }

        // Simulacion exacta de Argon 3.1.C literal en lil_arg:
        // Durante t < t_fin, alimentamos a lil_arg con salidas vacias (sin abrir ranuras prematuras),
        // manteniendo identicos el reloj y la historia del anillo; en t == t_fin dejamos pasar solo la
        // alarma de escala_real con signo_real, y desde t > t_fin corre libremente.
        if (!modo_parada_temprana) {
            ObsDelta::Salida sd_arg = sd;
            const std::int64_t st = static_cast<std::int64_t>(t);
            if (st >= t0 && st < t_fin) {
                for (std::size_t j = 0; j < sd_arg.niveles; ++j) {
                    sd_arg.nivel[j].alarma = EventoAlarma{};
                }
            } else if (st == t_fin) {
                // Argon §3.1.C l.218: "si el operador delta presenta alarma en escala_real con signo == signo_real ... abre"
                bool tiene_escala_real = false;
                EventoAlarma al_real{};
                if (ev.escala >= 0 && static_cast<std::size_t>(ev.escala) < sd.niveles &&
                    sd.nivel[static_cast<std::size_t>(ev.escala)].alarma.hay &&
                    sd.nivel[static_cast<std::size_t>(ev.escala)].alarma.signo == ev.signo) {
                    tiene_escala_real = true;
                    al_real = sd.nivel[static_cast<std::size_t>(ev.escala)].alarma;
                }
                for (std::size_t j = 0; j < sd_arg.niveles; ++j) sd_arg.nivel[j].alarma = EventoAlarma{};
                if (tiene_escala_real) {
                    // Restaurar t_cambio = t0 para respetar ancla [t0 - guarda - w_ancla, t0 - guarda - 1] (§2.1)
                    al_real.t_cambio = t0;
                    sd_arg.nivel[static_cast<std::size_t>(ev.escala)].alarma = al_real;
                }
            }
            const auto o_arg = lil_arg->paso(sd_arg, xt, true, NaN, 0, 0);
            for (std::size_t c = 0; c < o_arg.ranuras_cerradas; ++c) {
                const auto& rc = lil_arg->cerrada(c);
                if (rc.t_apertura == t_fin && rc.signo == ev.signo) {
                    const auto d = static_cast<std::uint8_t>(rc.cruda.desenlace);
                    if (d == 2 || d == 4) {
                        cp.arg_admitida = true;
                        cp.arg_tau_cod = detalle_nulo::codificar_tau(rc.cruda.tau_rest_consolidado, cs.t_max_ef);
                    }
                }
            }
        }
    }
    return cp;
}

struct MetricasPol {
    std::size_t intentos = 0;
    std::size_t admitidas = 0;
    std::size_t bloqueadas_hay_previa = 0;
    double sum_tau = 0.0;
    double sum_x_pre_sig = 0.0;
    double sum_u_pit = 0.0;
    std::size_t n_u = 0;
};

void registrar_diagnostico_grupo(
    const CalibracionSonda& cs,
    const std::vector<EventoCapturado>& eventos,
    std::size_t max_ev,
    std::size_t n_mc,
    std::uint64_t sem_base,
    const char* titulo) {

    std::vector<const EventoCapturado*> sel;
    for (const auto& ev : eventos) {
        if (sel.size() >= max_ev) break;
        sel.push_back(&ev);
    }
    if (sel.empty()) {
        std::printf("\n[%s] Sin eventos disponibles.\n", titulo);
        return;
    }

    MetricasPol m_base{}, m_nohp{}, m_arg{}, m_opA{}, m_opB{};
    std::size_t tot_real = 0;
    std::size_t glr_prem_cualq = 0, glr_prem_misma = 0;
    std::size_t glr_prem_pierde_tfin_misma = 0, glr_prem_pierde_tfin_cualq = 0;
    std::size_t glr_no_prem_tfin_misma = 0, tot_no_prem = 0;
    double sum_x_all = 0.0, sum_x_prem = 0.0;

    std::mt19937_64 rng_noise(detalle_nulo::splitmix64(sem_base));
    std::uniform_int_distribution<std::size_t> dist_res(0, cs.m_lin.residuos.size() - 1);

    for (std::size_t ei = 0; ei < sel.size(); ++ei) {
        const auto& ev = *sel[ei];
        const std::size_t burn = cs.m_lin.burn_in;
        // Margen estacionario real antes del ancla para que alarmas espontaneas en [win_ini, t0) tengan ancla >= 0
        const std::size_t w_base = 2 * (ev.w_ancla + ev.guarda) + L_GLR;
        const std::size_t w_sorteo = 64;
        const std::size_t w_post = ev.l + cs.t_max_ef + 8;
        const std::size_t n_sim = w_base + w_sorteo + w_post;
        const std::size_t total = burn + n_sim;
        std::uniform_int_distribution<std::int64_t> dist_t(
            static_cast<std::int64_t>(w_base),
            static_cast<std::int64_t>(w_base + w_sorteo - 1));

        std::vector<double> taus_base, taus_nohp, taus_arg, taus_opA, taus_opB;
        std::vector<double> eps_tot(total, 0.0), x_sim;

        for (std::size_t r = 0; r < n_mc; ++r) {
            const std::int64_t t0 = dist_t(rng_noise);
            const std::int64_t t1 = t0 + static_cast<std::int64_t>(ev.l);
            const std::int64_t t_fin = t1 - 1;
            const std::int64_t win_ini = t0 - static_cast<std::int64_t>(ev.w_ancla + ev.guarda);
            for (std::size_t k = 0; k < total; ++k) {
                eps_tot[k] = cs.m_lin.residuos[dist_res(rng_noise)] - cs.m_lin.media_residuos;
            }

            const auto cp1 = ejecutar_sim_unificada(cs, ev, eps_tot, burn, n_sim, t0, false, x_sim);
            const double x_pre_sig = x_sim[static_cast<std::size_t>(t0 - 1)] * static_cast<double>(ev.signo);
            ++tot_real;
            sum_x_all += x_pre_sig;

            if (cp1.glr_prem_cualq) {
                ++glr_prem_cualq;
                sum_x_prem += x_pre_sig;
                if (!cp1.glr_en_tfin_cualq) ++glr_prem_pierde_tfin_cualq;
            } else {
                ++tot_no_prem;
                if (cp1.glr_en_tfin_misma) ++glr_no_prem_tfin_misma;
            }
            if (cp1.glr_prem_misma) {
                ++glr_prem_misma;
                if (!cp1.glr_en_tfin_misma) ++glr_prem_pierde_tfin_misma;
            }

            bool hay_previa = false;
            for (std::int64_t ta : cp1.aperturas_todas) {
                if (ta >= win_ini && ta < t0) { hay_previa = true; break; }
            }

            const CierreSim* propia = nullptr;
            for (const auto& ci : cp1.cierres) {
                if (ci.t_ap >= t0 && ci.t_ap <= t1 && ci.signo == ev.signo) {
                    if (propia == nullptr || ci.t_ap < propia->t_ap) propia = &ci;
                }
            }

            ++m_base.intentos;
            ++m_nohp.intentos;
            ++m_arg.intentos;
            ++m_opA.intentos;

            if (propia != nullptr && (propia->des == 2 || propia->des == 4)) {
                if (hay_previa) {
                    ++m_base.bloqueadas_hay_previa;
                } else {
                    ++m_base.admitidas;
                    m_base.sum_tau += propia->tau_cod;
                    m_base.sum_x_pre_sig += x_pre_sig;
                    taus_base.push_back(propia->tau_cod);
                }
                ++m_nohp.admitidas;
                m_nohp.sum_tau += propia->tau_cod;
                m_nohp.sum_x_pre_sig += x_pre_sig;
                taus_nohp.push_back(propia->tau_cod);

                if (propia->t_ap == t_fin) {
                    ++m_opA.admitidas;
                    m_opA.sum_tau += propia->tau_cod;
                    m_opA.sum_x_pre_sig += x_pre_sig;
                    taus_opA.push_back(propia->tau_cod);
                }
            }
            if (cp1.arg_admitida) {
                ++m_arg.admitidas;
                m_arg.sum_tau += cp1.arg_tau_cod;
                m_arg.sum_x_pre_sig += x_pre_sig;
                taus_arg.push_back(cp1.arg_tau_cod);
            }

            // Opcion B: parada temprana online sobre el MISMO eps_tot
            const auto cp2 = ejecutar_sim_unificada(cs, ev, eps_tot, burn, n_sim, t0, true, x_sim);
            ++m_opB.intentos;
            const CierreSim* propia_b = nullptr;
            for (const auto& ci : cp2.cierres) {
                if (ci.t_ap >= t0 && ci.t_ap <= t_fin && ci.signo == ev.signo) {
                    if (propia_b == nullptr || ci.t_ap < propia_b->t_ap) propia_b = &ci;
                }
            }
            if (propia_b != nullptr && (propia_b->des == 2 || propia_b->des == 4)) {
                ++m_opB.admitidas;
                m_opB.sum_tau += propia_b->tau_cod;
                m_opB.sum_x_pre_sig += x_pre_sig;
                taus_opB.push_back(propia_b->tau_cod);
            }
        }

        auto acc_u = [&](MetricasPol& m, const std::vector<double>& taus) {
            if (taus.size() >= 8) {
                m.sum_u_pit += detalle_nulo::pit_punto_medio(ev.tau_real_cod, taus);
                ++m.n_u;
            }
        };
        acc_u(m_base, taus_base);
        acc_u(m_nohp, taus_nohp);
        acc_u(m_arg, taus_arg);
        acc_u(m_opA, taus_opA);
        acc_u(m_opB, taus_opB);
    }

    double l_med = 0.0;
    for (const auto* e : sel) l_med += static_cast<double>(e->l);
    l_med /= static_cast<double>(sel.size());

    std::printf("\n=== %s (K_ev = %zu, l_medio = %.2f, N_mc/ev = %zu, tot = %zu) ===\n",
                titulo, sel.size(), l_med, n_mc, tot_real);
    std::printf("  Disparos prematuros d_local en [t0, t0+l-2] (cualq escala): %zu/%zu (%.2f%%)\n",
                glr_prem_cualq, tot_real, 100.0 * static_cast<double>(glr_prem_cualq) / static_cast<double>(tot_real));
    std::printf("  Disparos prematuros d_local en escala_real+signo_real     : %zu/%zu (%.2f%%)\n",
                glr_prem_misma, tot_real, 100.0 * static_cast<double>(glr_prem_misma) / static_cast<double>(tot_real));
    if (glr_prem_misma > 0) {
        std::printf("  Entre prematuros (misma esc+sig), NO re-dispara en t_fin por reiniciar(): %zu/%zu (%.2f%%)\n",
                    glr_prem_pierde_tfin_misma, glr_prem_misma,
                    100.0 * static_cast<double>(glr_prem_pierde_tfin_misma) / static_cast<double>(glr_prem_misma));
    }
    if (tot_no_prem > 0) {
        std::printf("  Entre NO prematuros, si dispara en t_fin (misma esc+sig): %zu/%zu (%.2f%%)\n",
                    glr_no_prem_tfin_misma, tot_no_prem,
                    100.0 * static_cast<double>(glr_no_prem_tfin_misma) / static_cast<double>(tot_no_prem));
    }
    std::printf("  E[x_{t0-1} * signo] incondicional = %+.4f | en prematuros = %+.4f\n",
                sum_x_all / static_cast<double>(tot_real),
                glr_prem_cualq > 0 ? sum_x_prem / static_cast<double>(glr_prem_cualq) : NaN);

    auto imprimir_pol = [&](const char* nom, const MetricasPol& m) {
        const double pdet = m.intentos > 0 ? static_cast<double>(m.admitidas) / static_cast<double>(m.intentos) : 0.0;
        const double tmed = m.admitidas > 0 ? m.sum_tau / static_cast<double>(m.admitidas) : NaN;
        const double xmed = m.admitidas > 0 ? m.sum_x_pre_sig / static_cast<double>(m.admitidas) : NaN;
        const double smed = m.n_u > 0 ? (m.sum_u_pit / static_cast<double>(m.n_u) - 0.5) : NaN;
        std::printf("  %-24s | p_det=%.4f (%4zu/%4zu) | bloq_hp=%3zu | tau_med=%.3f | E[x_{t0-1}s]=%+.4f | S=%.4f (n=%zu)\n",
                    nom, pdet, m.admitidas, m.intentos, m.bloqueadas_hay_previa, tmed, xmed, smed, m.n_u);
    };
    imprimir_pol("V0_BASE_ACTUAL", m_base);
    imprimir_pol("V1_BASE_SIN_HAY_PREVIA", m_nohp);
    imprimir_pol("V2_ARGON_3_1_C_LITERAL", m_arg);
    imprimir_pol("V3_OPCION_A_EXACTA_TFIN", m_opA);
    imprimir_pol("V4_OPCION_B_PARADA_TEMP", m_opB);
}

struct AcumSemilla {
    double s[5]{};
    bool val[5]{};
    std::size_t k_ev[5]{};
    double p_det_sum[5]{};
    std::size_t n_pdet[5]{};
    std::size_t hp_bloq = 0;
    std::size_t tot_int = 0;
};

AcumSemilla correr_semilla_r1(std::size_t idx_sem) {
    AcumSemilla out{};
    const std::uint64_t sem_cal = 2026112100ULL + static_cast<std::uint64_t>(idx_sem) * 101ULL;
    const std::uint64_t sem_op = 2026120000ULL + static_cast<std::uint64_t>(idx_sem) * 313ULL;
    const auto cs = calibrar_sonda_ar1(sem_cal, 4096, 0.80);
    const auto eventos = extraer_eventos_reales(cs, sem_op, 4096, 0.80);
    if (eventos.empty()) return out;

    const std::size_t m_min = ark::metrologia::familia::m_min(ALFA_PROYECTO);  // 8
    const std::size_t max_int = 180;
    std::mt19937_64 rng(detalle_nulo::splitmix64(sem_op ^ 0xA5A5A5A5A5A5A5A5ULL));
    std::uniform_int_distribution<std::size_t> dist_res(0, cs.m_lin.residuos.size() - 1);

    std::vector<double> u_list[5];
    std::vector<double> x_sim;

    for (const auto& ev : eventos) {
        const std::size_t burn = cs.m_lin.burn_in;
        const std::size_t w_base = 2 * (ev.w_ancla + ev.guarda) + L_GLR;
        const std::size_t w_sorteo = 64;
        const std::size_t w_post = ev.l + cs.t_max_ef + 8;
        const std::size_t n_sim = w_base + w_sorteo + w_post;
        const std::size_t total = burn + n_sim;
        std::uniform_int_distribution<std::int64_t> dist_t(
            static_cast<std::int64_t>(w_base),
            static_cast<std::int64_t>(w_base + w_sorteo - 1));

        std::vector<double> taus[5];
        std::size_t int_pol[5]{};
        std::vector<double> eps_tot(total, 0.0);

        for (std::size_t it = 1; it <= max_int; ++it) {
            bool todos_llenos = true;
            for (int p = 0; p < 5; ++p) {
                if (taus[p].size() < m_min) { todos_llenos = false; break; }
            }
            if (todos_llenos) break;

            const std::int64_t t0 = dist_t(rng);
            const std::int64_t t1 = t0 + static_cast<std::int64_t>(ev.l);
            const std::int64_t t_fin = t1 - 1;
            const std::int64_t win_ini = t0 - static_cast<std::int64_t>(ev.w_ancla + ev.guarda);
            for (std::size_t k = 0; k < total; ++k) {
                eps_tot[k] = cs.m_lin.residuos[dist_res(rng)] - cs.m_lin.media_residuos;
            }

            if (taus[0].size() < m_min || taus[1].size() < m_min || taus[2].size() < m_min || taus[3].size() < m_min) {
                for (int p = 0; p < 4; ++p) if (taus[p].size() < m_min) ++int_pol[p];
                const auto cp1 = ejecutar_sim_unificada(cs, ev, eps_tot, burn, n_sim, t0, false, x_sim);
                bool hay_previa = false;
                for (std::int64_t ta : cp1.aperturas_todas) {
                    if (ta >= win_ini && ta < t0) { hay_previa = true; break; }
                }
                const CierreSim* propia = nullptr;
                for (const auto& ci : cp1.cierres) {
                    if (ci.t_ap >= t0 && ci.t_ap <= t1 && ci.signo == ev.signo) {
                        if (propia == nullptr || ci.t_ap < propia->t_ap) propia = &ci;
                    }
                }
                ++out.tot_int;
                if (propia != nullptr && (propia->des == 2 || propia->des == 4)) {
                    if (hay_previa) {
                        ++out.hp_bloq;
                    } else if (taus[0].size() < m_min) {
                        taus[0].push_back(propia->tau_cod);
                    }
                    if (taus[1].size() < m_min) taus[1].push_back(propia->tau_cod);
                    if (propia->t_ap == t_fin && taus[3].size() < m_min) taus[3].push_back(propia->tau_cod);
                }
                if (cp1.arg_admitida && taus[2].size() < m_min) {
                    taus[2].push_back(cp1.arg_tau_cod);
                }
            }

            if (taus[4].size() < m_min) {
                ++int_pol[4];
                const auto cp2 = ejecutar_sim_unificada(cs, ev, eps_tot, burn, n_sim, t0, true, x_sim);
                const CierreSim* propia_b = nullptr;
                for (const auto& ci : cp2.cierres) {
                    if (ci.t_ap >= t0 && ci.t_ap <= t_fin && ci.signo == ev.signo) {
                        if (propia_b == nullptr || ci.t_ap < propia_b->t_ap) propia_b = &ci;
                    }
                }
                if (propia_b != nullptr && (propia_b->des == 2 || propia_b->des == 4)) {
                    taus[4].push_back(propia_b->tau_cod);
                }
            }
        }

        for (int p = 0; p < 5; ++p) {
            if (int_pol[p] > 0) {
                out.p_det_sum[p] += static_cast<double>(taus[p].size()) / static_cast<double>(int_pol[p]);
                ++out.n_pdet[p];
            }
            if (taus[p].size() >= m_min) {
                u_list[p].push_back(detalle_nulo::pit_punto_medio(ev.tau_real_cod, taus[p]));
            }
        }
    }

    for (int p = 0; p < 5; ++p) {
        out.k_ev[p] = u_list[p].size();
        if (!u_list[p].empty()) {
            double su = 0.0;
            for (double u : u_list[p]) su += u;
            out.s[p] = su / static_cast<double>(u_list[p].size()) - 0.5;
            out.val[p] = true;
        }
    }
    return out;
}

}  // namespace

int main() {
    const auto t_ini_clock = std::chrono::steady_clock::now();
    std::printf("================================================================================\n");
    std::printf(" SONDA EMPIRICA: DetectorGLR::reiniciar() Y hay_previa EN VIA (b) (AR(1) phi=0.8)\n");
    std::printf("================================================================================\n");

    const auto cs1 = calibrar_sonda_ar1(20261121ULL, 4096, 0.80);
    std::printf("Calibracion 4096 barras: t_max_ef=%zu, p=%zu, phi_1=%.4f, sd_innov=%.4f, radio_max=%.4f\n",
                cs1.t_max_ef, cs1.m_lin.p, cs1.m_lin.p > 0 ? cs1.m_lin.phi[0] : 0.0,
                std::sqrt(cs1.m_lin.sigma2_innov), cs1.m_lin.radio_max);

    const auto eventos1 = extraer_eventos_reales(cs1, 20261199ULL, 16384, 0.80);
    std::vector<EventoCapturado> ev_l1, ev_l2p;
    std::size_t max_l = 0;
    for (const auto& ev : eventos1) {
        if (ev.l > max_l) max_l = ev.l;
        if (ev.l == 1) ev_l1.push_back(ev);
        else ev_l2p.push_back(ev);
    }
    std::printf("Serie operativa 16384 barras: %zu eventos cerrados medibles (l==1: %zu, l>=2: %zu, l_max=%zu)\n",
                eventos1.size(), ev_l1.size(), ev_l2p.size(), max_l);

    registrar_diagnostico_grupo(cs1, ev_l2p, 12, 500, 10001ULL, "DIAGNOSTICO 1A: PERFILES MULTI-BARRA (l >= 2)");
    registrar_diagnostico_grupo(cs1, ev_l1, 12, 500, 20002ULL, "DIAGNOSTICO 1B: CONTROL PERFILES 1-BARRA (l == 1)");

    // Diagnostico 2: 16 semillas independientes R1 en paralelo
    constexpr std::size_t N_SEM = 16;
    std::vector<AcumSemilla> acs(N_SEM);
    std::atomic<std::size_t> sig{0};
    const std::size_t n_hilos = std::min<std::size_t>(N_SEM, std::max<std::size_t>(1, std::thread::hardware_concurrency()));
    std::vector<std::thread> trab;
    for (std::size_t h = 0; h < n_hilos; ++h) {
        trab.emplace_back([&]() {
            for (;;) {
                const std::size_t i = sig.fetch_add(1, std::memory_order_relaxed);
                if (i >= N_SEM) break;
                acs[i] = correr_semilla_r1(i);
            }
        });
    }
    for (auto& th : trab) th.join();

    const char* nombres[5] = {
        "V0_BASE_ACTUAL",
        "V1_BASE_SIN_HAY_PREVIA",
        "V2_ARGON_3_1_C_LITERAL",
        "V3_OPCION_A_EXACTA_TFIN",
        "V4_OPCION_B_PARADA_TEMP"
    };
    std::printf("\n=== DIAGNOSTICO 2: BANCO MULTI-SEMILLA R1 (AR(1) phi=0.80, %zu semillas x 4096 barras, m_min=%zu) ===\n",
                N_SEM, ark::metrologia::familia::m_min(ALFA_PROYECTO));
    for (int p = 0; p < 5; ++p) {
        std::vector<double> ss;
        std::size_t k_tot = 0, s_neg = 0;
        double pd_sum = 0.0;
        std::size_t pd_n = 0;
        for (const auto& a : acs) {
            if (a.val[p]) {
                ss.push_back(a.s[p]);
                if (a.s[p] < 0.0) ++s_neg;
            }
            k_tot += a.k_ev[p];
            pd_sum += a.p_det_sum[p];
            pd_n += a.n_pdet[p];
        }
        double mean = 0.0, ee = 0.0;
        for (double v : ss) mean += v;
        if (!ss.empty()) mean /= static_cast<double>(ss.size());
        if (ss.size() >= 2) {
            double vv = 0.0;
            for (double v : ss) vv += (v - mean) * (v - mean);
            ee = std::sqrt(vv / static_cast<double>(ss.size() * (ss.size() - 1)));
        }
        const double z = ee > 0.0 ? mean / ee : 0.0;
        std::printf("  %-24s | n_sem=%2zu K_tot=%3zu | S_medio=%+.5f ee=%.5f z=%+.2f | S<0: %2zu/%2zu | p_det=%.4f\n",
                    nombres[p], ss.size(), k_tot, mean, ee, z, s_neg, ss.size(),
                    pd_n > 0 ? pd_sum / static_cast<double>(pd_n) : 0.0);
    }
    std::size_t hp_tot = 0, int_tot = 0;
    for (const auto& a : acs) { hp_tot += a.hp_bloq; int_tot += a.tot_int; }
    std::printf("  Descartes exclusivos por hay_previa en V0: %zu / %zu intentos efectivos (%.2f%%)\n",
                hp_tot, int_tot, int_tot > 0 ? 100.0 * static_cast<double>(hp_tot) / static_cast<double>(int_tot) : 0.0);

    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_ini_clock).count();
    std::printf("\nTiempo total: %.2f s\n", elapsed);
    return 0;
}
