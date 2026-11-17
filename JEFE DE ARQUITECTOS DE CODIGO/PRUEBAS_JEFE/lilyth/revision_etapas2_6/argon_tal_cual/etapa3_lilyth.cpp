// etapa3_lilyth.cpp -- Controles de la Etapa 3 del Operador de Lilyth (nulos (a) y (b), T_max, PIT de punto medio,
// pruebas del resumen y controles S1-S3) CONTRA EL OPERADOR DELTA REAL. C++20.
//
// Controles:
//   K1   identidad del operador delta con y sin Lilyth (y cero reservas en paso tras calibrar)
//   P1   ruido iid con disparo real: tau_rest crudo ~1 paso, Delta_tau_iid compatible con 0, y referencia independiente bit a bit
//   P2   AR(1) phi = 0.8 con choques: Delta_tau_iid recupera ln 2 / (-ln phi) corregido por el adelanto de primer cruce; Delta_tau_Onsager ~ 0
//   P2b  el mismo mundo forzando la via (b): la semivida del nulo generativo inyectado sale dentro del error de la de (a)
//   P3a  TAR con la misma ACF de calma: Delta_tau_Onsager > 0 con el signo del umbral (p < alfa)
//   P3b  AR(1) y AR(2): eta_barra(tau) dentro de la banda de 2 ee de la ACF de la calma
//   P10  paseo aleatorio en NIVEL: T_max = W_ancla marcado (t_max_raiz_unitaria) y Delta_tau_Onsager compatible con 0
//   P12  BLOQUES se rechaza en calibrar; FASE_ALEATORIZADA_FFT en (b) cae al generativo y la ficha lo marca
//   N1   tamano: en un mundo nulo, las dos pruebas rechazan a tasa <= alfa (+ su error de Monte Carlo) sobre >= 200 semillas
//   N2   independencia de realizaciones: con 1 hilo y con 4 hilos los nulos y estadisticos son identicos bit a bit
//   N3   cada realizacion parte de una copia nueva: sembrar una realizacion con estado previo de PurgaDiurna cambia el resultado (Regla 7)
//   S1   replicas de Hawkes con razon de ramificacion n en {0.0, 0.3, 0.6} sobre AR(1) de calma conocida
//   S2   ventanas solapadas por debajo y por encima de W_ancla + guarda + l + T_MAX, y efecto de contaminacion del ancla
//   S3   integridad de la parada de la via (b): ley binomial negativa, tope ceil(m_min / p_det) de las primeras m_min inyecciones y p_det = 0

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

#include "ark/lilyth/observador_lilyth.hpp"

static std::size_t g_reservas = 0;
static bool g_contar = false;
void* operator new(std::size_t n) {
    if (g_contar) ++g_reservas;
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc{};
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

using namespace ark::metrologia;
using namespace ark::lilyth;

namespace {

using Obs = ObservadorSaltos<24, 64, 120>;
constexpr std::size_t TMAXC = 128;
using Lil = ObservadorLilyth<Obs, 8, TMAXC, 2048>;

int fallos = 0;
void chk(bool ok, const char* id, const char* txt) {
    std::printf("[%s] %-4s %s\n", ok ? "OK   " : "FALLA", id, txt);
    if (!ok) ++fallos;
}
bool igual(double a, double b) { return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof a) == 0; }

// Motor de referencia independiente para PIT de punto medio, excedencia espontanea y T_max (Regla 3).
double pit_ref_independiente(double x, const std::vector<double>& v) {
    std::size_t lt = 0, eq = 0;
    for (double y : v) {
        if (y < x - 1e-12) ++lt;
        else if (std::abs(y - x) <= 1e-12) ++eq;
    }
    return (static_cast<double>(lt) + 0.5 * static_cast<double>(eq)) / static_cast<double>(v.size());
}

std::size_t tmax_ref_independiente(const std::vector<EventoNuloLilyth<TMAXC>>& nulos, double alfa) {
    if (nulos.empty()) return TMAXC;
    std::size_t cens = 0;
    std::vector<double> taus;
    for (const auto& e : nulos) {
        if (std::isfinite(e.tau_rest_primer_toque) && e.tau_rest_primer_toque <= static_cast<double>(TMAXC)) {
            taus.push_back(e.tau_rest_primer_toque);
        } else {
            ++cens;
            taus.push_back(static_cast<double>(TMAXC));
        }
    }
    if (static_cast<double>(cens) / static_cast<double>(nulos.size()) > alfa) return TMAXC;
    std::sort(taus.begin(), taus.end());
    const double pos = (1.0 - alfa) * static_cast<double>(taus.size() - 1);
    const std::size_t i = static_cast<std::size_t>(std::floor(pos));
    const double fr = pos - static_cast<double>(i);
    const double q = (i + 1 < taus.size()) ? (taus[i] + fr * (taus[i + 1] - taus[i])) : taus.back();
    return std::clamp<std::size_t>(static_cast<std::size_t>(std::ceil(q)), 2, TMAXC);
}

// Genera una serie AR(1) estacionaria de longitud n con sd = 1.0.
std::vector<double> generar_ar1(double phi, std::size_t n, std::uint64_t sem) {
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> nd(0.0, std::sqrt(std::max(1e-6, 1.0 - phi * phi)));
    std::vector<double> x(n, 0.0);
    double cur = 0.0;
    for (std::size_t b = 0; b < 100; ++b) cur = phi * cur + nd(rng);
    for (std::size_t t = 0; t < n; ++t) {
        cur = phi * cur + nd(rng);
        x[t] = cur;
    }
    return x;
}

// Inyecta choques espaciados en las innovaciones de un AR(1) para producir eventos reales con verdad conocida.
std::vector<double> generar_ar1_con_choques(double phi, std::size_t n, const std::vector<std::size_t>& t_choques,
                                            double amp_innov, std::size_t dur_innov, std::uint64_t sem) {
    std::mt19937_64 rng(sem);
    const double sd_e = std::sqrt(std::max(1e-6, 1.0 - phi * phi));
    std::normal_distribution<double> nd(0.0, sd_e);
    std::vector<double> eps(n, 0.0);
    for (std::size_t t = 0; t < n; ++t) eps[t] = nd(rng);
    for (std::size_t idx = 0; idx < t_choques.size(); ++idx) {
        const std::size_t t0 = t_choques[idx];
        const double sgn = (idx % 2 == 0) ? 1.0 : -1.0;
        for (std::size_t k = 0; k < dur_innov && (t0 + k) < n; ++k) {
            eps[t0 + k] += sgn * amp_innov;
        }
    }
    std::vector<double> x(n, 0.0);
    double cur = 0.0;
    for (std::size_t t = 0; t < n; ++t) {
        cur = phi * cur + eps[t];
        x[t] = cur;
    }
    return x;
}

}  // namespace

int main() {
    // -------------------------------------------------------------------------
    // K1 + P1: ruido iid con disparo real, identidad del operador delta,
    //          cero reservas en paso y casamiento bit a bit con referencia independiente
    // -------------------------------------------------------------------------
    {
        const std::size_t n_cal = 5000, n_test = 3000;
        const auto x_cal = generar_ar1(0.0, n_cal, 101);
        const auto x_test = generar_ar1_con_choques(0.0, n_test, {400, 900, 1400, 1900}, 4.5, 2, 102);

        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 350.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 2;

        auto d_solo = std::make_unique<Obs>();
        d_solo->calibrar(cd, x_cal);
        auto d_con = std::make_unique<Obs>();
        d_con->calibrar(cd, x_cal);

        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 6;
        cl.hilos_nulo = 1;
        cl.semilla_nulo = 20261106;
        const bool cal_ok = lil->calibrar(*d_con, cl, cd, x_cal);

        bool id_ok = cal_ok;
        g_reservas = 0;
        for (std::size_t t = 0; t < n_test; ++t) {
            const auto a = d_solo->paso(x_test[t], NaN, true);
            const auto b = d_con->paso(x_test[t], NaN, true);
            g_contar = true;
            (void)lil->paso(b, x_test[t], true);
            g_contar = false;
            id_ok = id_ok && a.t == b.t && a.n_eventos == b.n_eventos && a.evento_abierto == b.evento_abierto && a.niveles == b.niveles;
            for (std::size_t j = 0; id_ok && j < a.niveles; ++j) {
                id_ok = a.nivel[j].alarma.hay == b.nivel[j].alarma.hay && igual(a.nivel[j].g, b.nivel[j].g);
            }
        }
        lil->completar_nulos();
        chk(id_ok && g_reservas == 0, "K1", "identidad del operador delta con y sin Lilyth y cero reservas en paso");

        // Verificar P1 y comparar T_max, excedencia (escala, signo y g >= g_k) y PIT de punto medio contra referencia independiente
        const auto& fn = lil->ficha_nulo_cruda();
        const std::size_t tmax_ref = tmax_ref_independiente(fn.nulos_lineal, cl.alfa);
        bool ref_match = (fn.t_max_derivado == tmax_ref && lil->configuracion().t_max_efectivo == tmax_ref);

        // Tambien probar sobre eventos espontaneos de la calma para que g_k este en el interior de la distribucion espontanea
        // y cace cualquier inversion de g >= g_k (M01) o mezcla de escala/signo (M02).
        auto d_calma = std::make_unique<Obs>();
        d_calma->calibrar(cd, x_cal);
        auto lil_calma = std::make_unique<Lil>();
        lil_calma->calibrar(*d_calma, cl, cd, x_cal);
        const auto x_iid_puro = generar_ar1(0.0, 3000, 103);
        for (double v : x_iid_puro) {
            const auto s = d_calma->paso(v, NaN, true);
            (void)lil_calma->paso(s, v, true);
        }
        lil_calma->completar_nulos();
        const auto& fn_c = lil_calma->ficha_nulo_cruda();
        for (const auto& r : lil_calma->cerradas_historico()) {
            if (r.cruda.desenlace != DesenlaceRanura::RESTITUIDA && r.cruda.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
            std::size_t cuenta_ref_iid = 0, cuenta_ref_lin = 0;
            std::vector<double> taus_ref_iid;
            for (const auto& ev : fn_c.nulos_sin_memoria) {
                if (ev.escala_origen == r.escala_origen && ev.signo == r.signo && ev.g_disparo >= r.g_disparo) {
                    ++cuenta_ref_iid;
                    const double tc = (std::isfinite(ev.tau_rest_consolidado) && ev.tau_rest_consolidado <= double(fn_c.t_max_derivado))
                                          ? ev.tau_rest_consolidado
                                          : double(fn_c.t_max_derivado + 1);
                    taus_ref_iid.push_back(tc);
                }
            }
            for (const auto& ev : fn_c.nulos_lineal) {
                if (ev.escala_origen == r.escala_origen && ev.signo == r.signo && ev.g_disparo >= r.g_disparo) {
                    ++cuenta_ref_lin;
                }
            }
            if (r.cruda.nulo_iid.n_nulos_a != cuenta_ref_iid || r.cruda.nulo_onsager.n_nulos_a != cuenta_ref_lin) {
                ref_match = false;
            }
            if (cuenta_ref_iid >= fn_c.m_min) {
                const double tau_r = (std::isfinite(r.cruda.tau_rest_consolidado) && r.cruda.tau_rest_consolidado <= double(fn_c.t_max_derivado))
                                         ? r.cruda.tau_rest_consolidado
                                         : double(fn_c.t_max_derivado + 1);
                const double u_ref = pit_ref_independiente(tau_r, taus_ref_iid);
                if (std::abs(r.cruda.nulo_iid.u_pit - u_ref) > 1e-12) ref_match = false;
            }
        }

        const bool p1_ok = ref_match && fn_c.contraste_iid.valido &&
                           std::abs(fn_c.contraste_iid.mediana_delta_tau) <= 1.0 &&
                           fn_c.contraste_iid.p_bilateral >= cl.alfa;
        std::printf("      P1: T_max=%zu (ref %zu) med_delta_iid=%.2f p_iid=%.3f p_ons=%.3f\n",
                    fn_c.t_max_derivado, tmax_ref, fn_c.contraste_iid.mediana_delta_tau,
                    fn_c.contraste_iid.p_bilateral, fn_c.contraste_onsager.p_bilateral);
        chk(p1_ok, "P1", "ruido iid: tau_rest ~1 paso, Delta_tau_iid compatible con 0 y excedencia/PIT/T_max = referencia");
    }

    // -------------------------------------------------------------------------
    // P2 y P2b: AR(1) phi = 0.8 con choques (via (a) y forzando via (b))
    // -------------------------------------------------------------------------
    {
        const double phi = 0.8;
        const std::size_t n_cal = 5000, n_test = 3200;
        const auto x_cal = generar_ar1(phi, n_cal, 201);
        // Choques moderados para que en P2 la via (a) con B=12 alcance m_min, y en P2b se fuerce la via (b) sobre el mismo mundo
        const auto x_test = generar_ar1_con_choques(phi, n_test, {400, 900, 1400, 1900, 2400, 2800}, 2.8, 2, 202);

        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 180.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 2;

        auto d_a = std::make_unique<Obs>();
        d_a->calibrar(cd, x_cal);

        auto lil_a = std::make_unique<Lil>();
        ConfiguracionLilyth cl_a{};
        cl_a.representacion = ModoRepresentacionLilyth::CRUDA;
        cl_a.b_max_nulo = 14;
        cl_a.hilos_nulo = 1;
        cl_a.forzar_via_b = false;
        lil_a->calibrar(*d_a, cl_a, cd, x_cal);

        for (double v : x_test) {
            const auto s = d_a->paso(v, NaN, true);
            (void)lil_a->paso(s, v, true);
        }
        lil_a->completar_nulos();

        // Verificar tambien en el pool nulo de P2 que T_max, excedencia (escala, signo, g >= g_k) y PIT coinciden con la referencia
        const auto& fn_a = lil_a->ficha_nulo_cruda();
        const std::size_t tmax_a_ref = tmax_ref_independiente(fn_a.nulos_lineal, cl_a.alfa);
        bool ref_a_ok = (fn_a.t_max_derivado == tmax_a_ref);
        for (const auto& r : lil_a->cerradas_historico()) {
            if (r.cruda.desenlace != DesenlaceRanura::RESTITUIDA && r.cruda.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
            std::size_t c_lin = 0;
            for (const auto& ev : fn_a.nulos_lineal) {
                if (ev.escala_origen == r.escala_origen && ev.signo == r.signo && ev.g_disparo >= r.g_disparo) ++c_lin;
            }
            if (r.cruda.nulo_onsager.n_nulos_a != c_lin) ref_a_ok = false;
        }

        // P2b: el mismo mundo forzando la via (b)
        auto d_b = std::make_unique<Obs>();
        d_b->calibrar(cd, x_cal);
        auto lil_b = std::make_unique<Lil>();
        ConfiguracionLilyth cl_b = cl_a;
        cl_b.forzar_via_b = true;
        lil_b->calibrar(*d_b, cl_b, cd, x_cal);
        for (double v : x_test) {
            const auto s = d_b->paso(v, NaN, true);
            (void)lil_b->paso(s, v, true);
        }
        lil_b->completar_nulos();
        const auto& fn_b = lil_b->ficha_nulo_cruda();

        std::vector<double> med_a, med_b;
        bool todas_por_inyeccion_b = true;
        for (const auto& r : lil_a->cerradas_historico()) {
            if ((r.cruda.desenlace == DesenlaceRanura::RESTITUIDA || r.cruda.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) && r.cruda.nulo_onsager.evaluado) med_a.push_back(r.cruda.nulo_onsager.mediana_nulo);
        }
        for (const auto& r : lil_b->cerradas_historico()) {
            if ((r.cruda.desenlace == DesenlaceRanura::RESTITUIDA || r.cruda.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) && r.cruda.nulo_onsager.evaluado) {
                med_b.push_back(r.cruda.nulo_onsager.mediana_nulo);
                if (!r.nulo_por_inyeccion || !r.cruda.nulo_onsager.nulo_por_inyeccion) todas_por_inyeccion_b = false;
            }
        }
        const double semivida_nulo_a = est::cuantil(med_a, 0.5);
        const double semivida_nulo_b = est::cuantil(med_b, 0.5);

        std::printf("      P2:  T_max=%zu (ref %zu) med_delta_iid=%.2f (p=%.4f) med_delta_ons=%.2f (p=%.3f) semivida_a=%.2f\n",
                    fn_a.t_max_derivado, tmax_a_ref, fn_a.contraste_iid.mediana_delta_tau, fn_a.contraste_iid.p_bilateral,
                    fn_a.contraste_onsager.mediana_delta_tau, fn_a.contraste_onsager.p_bilateral, semivida_nulo_a);
        std::printf("      P2b: semivida_b=%.2f (dif=%.2f) med_delta_ons=%.2f (p=%.3f)\n",
                    semivida_nulo_b, std::abs(semivida_nulo_b - semivida_nulo_a),
                    fn_b.contraste_onsager.mediana_delta_tau, fn_b.contraste_onsager.p_bilateral);

        chk(ref_a_ok && fn_a.contraste_iid.mediana_delta_tau >= 1.0 && fn_a.contraste_iid.p_bilateral < cl_a.alfa &&
            fn_a.contraste_onsager.p_bilateral >= cl_a.alfa && semivida_nulo_a >= 2.0,
            "P2", "AR(1) phi=0.8: Delta_tau_iid > 0 (p < alfa), Delta_tau_Onsager compatible con 0 y excedencia = ref");
        chk(todas_por_inyeccion_b && semivida_nulo_b >= 2.0 && std::abs(semivida_nulo_b - semivida_nulo_a) <= 1.5 &&
            fn_b.contraste_onsager.p_bilateral >= cl_b.alfa,
            "P2b", "mismo mundo forzando via (b): nulo_por_inyeccion = true y semivida de (b) coincide con la de (a)");
    }

    // -------------------------------------------------------------------------
    // P3a: TAR (umbral no lineal) -> Delta_tau_Onsager != 0 con el signo del umbral
    // -------------------------------------------------------------------------
    {
        // En calma (|x| <= 2.2) actua como AR(1) moderado (phi_calma = 0.35); en excursiones grandes (|x| > 2.2) se vuelve muy pegajoso (phi_ext = 0.96)
        std::mt19937_64 rng(301);
        std::normal_distribution<double> nd(0.0, 0.9);
        const std::size_t n_cal = 4096, n_test = 3500;
        std::vector<double> x_cal(n_cal, 0.0);
        double cur = 0.0;
        for (std::size_t t = 0; t < n_cal; ++t) {
            cur = 0.35 * cur + nd(rng);
            x_cal[t] = cur;
        }
        std::vector<double> x_test(n_test, 0.0);
        cur = 0.0;
        for (std::size_t t = 0; t < n_test; ++t) {
            const double phi_loc = (std::abs(cur) > 2.2) ? 0.96 : 0.35;
            double e = nd(rng);
            for (std::size_t t0 : {300u, 800u, 1300u, 1800u, 2300u, 2800u}) {
                if (t == t0 || t == t0 + 1) e += 4.5;
            }
            cur = phi_loc * cur + e;
            x_test[t] = cur;
        }

        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 400.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d = std::make_unique<Obs>();
        d->calibrar(cd, x_cal);

        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 6;
        cl.hilos_nulo = 1;
        lil->calibrar(*d, cl, cd, x_cal);
        for (double v : x_test) {
            const auto s = d->paso(v, NaN, true);
            (void)lil->paso(s, v, true);
        }
        lil->completar_nulos();
        const auto& fn = lil->ficha_nulo_cruda();
        std::printf("      P3a: TAR med_delta_ons=%.2f S_ons=%.3f p_ons=%.4f\n",
                    fn.contraste_onsager.mediana_delta_tau, fn.contraste_onsager.estadistico_s, fn.contraste_onsager.p_bilateral);
        chk(fn.contraste_onsager.valido && fn.contraste_onsager.mediana_delta_tau > 2.0 &&
            fn.contraste_onsager.estadistico_s > 0.15 && fn.contraste_onsager.p_bilateral < cl.alfa,
            "P3a", "TAR no lineal con la misma ACF de calma: Delta_tau_Onsager > 0 con p < alfa");
    }

    // -------------------------------------------------------------------------
    // P3b: AR(1) y AR(2) lineales -> eta_barra dentro de la banda de la ACF
    // -------------------------------------------------------------------------
    {
        const std::size_t n_cal = 4096;
        const auto x_cal_ar1 = generar_ar1(0.75, n_cal, 311);
        const auto x_test_ar1 = generar_ar1_con_choques(0.75, 4000, {300, 700, 1100, 1500, 1900, 2300, 2700, 3100, 3500}, 4.5, 1, 312);

        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 400.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;

        auto d1 = std::make_unique<Obs>();
        d1->calibrar(cd, x_cal_ar1);
        auto l1 = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 4;
        cl.hilos_nulo = 1;
        l1->calibrar(*d1, cl, cd, x_cal_ar1);
        for (double v : x_test_ar1) {
            const auto s = d1->paso(v, NaN, true);
            (void)l1->paso(s, v, true);
        }
        l1->completar_nulos();

        // AR(2): x_t = 0.5 x_{t-1} + 0.2 x_{t-2} + e_t
        std::mt19937_64 rng(321);
        std::normal_distribution<double> nd(0.0, 0.8);
        std::vector<double> x_cal_ar2(n_cal, 0.0), x_test_ar2(6000, 0.0);
        for (std::size_t t = 2; t < n_cal; ++t) x_cal_ar2[t] = 0.5 * x_cal_ar2[t - 1] + 0.2 * x_cal_ar2[t - 2] + nd(rng);
        for (std::size_t t = 2; t < 6000; ++t) {
            double e = nd(rng);
            for (std::size_t t0 = 250; t0 + 200 < 6000; t0 += 220) {
                if (t == t0 - 1) e += 3.0;
                if (t == t0) e += 3.2;
            }
            x_test_ar2[t] = 0.5 * x_test_ar2[t - 1] + 0.2 * x_test_ar2[t - 2] + e;
        }
        auto d2 = std::make_unique<Obs>();
        d2->calibrar(cd, x_cal_ar2);
        auto l2 = std::make_unique<Lil>();
        l2->calibrar(*d2, cl, cd, x_cal_ar2);
        for (double v : x_test_ar2) {
            const auto s = d2->paso(v, NaN, true);
            (void)l2->paso(s, v, true);
        }
        l2->completar_nulos();

        const auto& p3b_1 = l1->ficha_nulo_cruda().p3b;
        const auto& p3b_2 = l2->ficha_nulo_cruda().p3b;
        std::printf("      P3b: AR(1) max_z=%.2f (n=%zu) | AR(2) max_z=%.2f (n=%zu)\n",
                    p3b_1.max_z_desvio, p3b_1.n_eventos, p3b_2.max_z_desvio, p3b_2.n_eventos);
        chk(p3b_1.dentro_de_banda_2ee && p3b_2.dentro_de_banda_2ee,
            "P3b", "AR(1) y AR(2): eta_barra(tau) sigue la ACF de la calma dentro de la banda de bootstrap");
    }

    // -------------------------------------------------------------------------
    // P10: paseo aleatorio en NIVEL -> T_max = W_ancla marcado y Delta_tau ~ 0
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(401);
        std::normal_distribution<double> nd(0.0, 1.0);
        const std::size_t n_cal = 4096, n_test = 2500;
        std::vector<double> x_cal(n_cal, 100.0), x_test(n_test, 100.0);
        for (std::size_t t = 1; t < n_cal; ++t) x_cal[t] = x_cal[t - 1] + nd(rng);
        for (std::size_t t = 1; t < n_test; ++t) {
            double inc = nd(rng);
            for (std::size_t t0 : {400u, 1000u, 1600u}) if (t == t0) inc += 8.0;
            x_test[t] = x_test[t - 1] + inc;
        }

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::NIVEL;
        cd.arl0_pasos_nativos = 400.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d = std::make_unique<Obs>();
        d->calibrar(cd, x_cal);

        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 4;
        cl.hilos_nulo = 1;
        lil->calibrar(*d, cl, cd, x_cal);
        for (double v : x_test) {
            const auto s = d->paso(v, NaN, true);
            (void)lil->paso(s, v, true);
        }
        lil->completar_nulos();
        const auto& fn = lil->ficha_nulo_cruda();
        std::printf("      P10: T_max=%zu raiz_unitaria=%d med_delta_ons=%.2f p_ons=%.3f\n",
                    lil->configuracion().t_max_efectivo, int(fn.t_max_raiz_unitaria),
                    fn.contraste_onsager.mediana_delta_tau, fn.contraste_onsager.p_bilateral);
        chk(fn.t_max_raiz_unitaria && lil->configuracion().t_max_efectivo == 120 &&
            fn.contraste_onsager.valido && fn.contraste_onsager.p_bilateral >= cl.alfa,
            "P10", "paseo aleatorio en NIVEL: T_max = W_ancla (120) marcado como raiz unitaria y Delta_tau_Onsager ~ 0");
    }

    // -------------------------------------------------------------------------
    // P12: BLOQUES se rechaza; FASE_ALEATORIZADA_FFT cae al generativo en (b)
    // -------------------------------------------------------------------------
    {
        const auto x_cal = generar_ar1(0.6, 3000, 501);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 400.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d = std::make_unique<Obs>();
        d->calibrar(cd, x_cal);

        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.nulo_sin_memoria = LILYTH_NULO_BLOQUES;
        const bool rechaza_bloques = !lil->calibrar(*d, cl, cd, x_cal);

        cl.nulo_sin_memoria = LILYTH_NULO_PERMUTACION;
        cl.modo_nulo_lineal = ModoNuloLineal::FASE_ALEATORIZADA_FFT;
        cl.b_max_nulo = 2;
        cl.hilos_nulo = 1;
        const bool acepta_fft = lil->calibrar(*d, cl, cd, x_cal);
        const auto& fn = lil->ficha_nulo_cruda();
        chk(rechaza_bloques && acepta_fft && fn.nulo_fft_degradado_a_generativo && fn.modelo_lineal.cayo_a_generativo,
            "P12", "BLOQUES se rechaza en calibrar; FASE_ALEATORIZADA_FFT cae al modelo generativo y la ficha lo marca");
    }

    // -------------------------------------------------------------------------
    // N1: tamano bajo el nulo sobre >= 200 semillas independientes
    // -------------------------------------------------------------------------
    {
        const std::size_t n_cal = 4096;
        const auto x_cal = generar_ar1(0.0, n_cal, 601);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 120.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d_base = std::make_unique<Obs>();
        d_base->calibrar(cd, x_cal);

        auto lil_base = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 24;
        cl.b_mc_resumen = 199;
        cl.hilos_nulo = 1;
        lil_base->calibrar(*d_base, cl, cd, x_cal);

        const std::size_t N_SEMILLAS = 215;
        std::size_t rechazos_iid = 0, rechazos_ons = 0, validas = 0;
        for (std::size_t s_idx = 0; s_idx < N_SEMILLAS; ++s_idx) {
            std::mt19937_64 rng_s(7000 + s_idx);
            std::vector<double> x_s;
            detalle_nulo::construir_serie_nula(Integracion::FLUJO, false, cl.nulo_sin_memoria,
                lil_base->ficha_nulo_cruda().segmentos_calma, lil_base->ficha_nulo_cruda().modelo_lineal, 0.0, 900, rng_s, x_s);
            Obs ds = *d_base;
            ds.reiniciar_flujo();
            Lil ls = *lil_base;
            ls.reiniciar_flujo_conservando_purga();
            for (double v : x_s) {
                const auto sal = ds.paso(v, NaN, true);
                (void)ls.paso(sal, v, true);
            }
            ls.completar_nulos();
            const auto& fn = ls.ficha_nulo_cruda();
            if (fn.contraste_iid.valido && fn.contraste_onsager.valido) {
                ++validas;
                if (fn.contraste_iid.p_bilateral < cl.alfa) ++rechazos_iid;
                if (fn.contraste_onsager.p_bilateral < cl.alfa) ++rechazos_ons;
            }
        }
        const double tasa_iid = static_cast<double>(rechazos_iid) / static_cast<double>(std::max<std::size_t>(1, validas));
        const double tasa_ons = static_cast<double>(rechazos_ons) / static_cast<double>(std::max<std::size_t>(1, validas));
        const double cota_mc = cl.alfa + 2.5 * std::sqrt(cl.alfa * (1.0 - cl.alfa) / static_cast<double>(validas));
        std::printf("      N1 (%zu semillas validas): tasa_rechazo_iid=%.3f, tasa_rechazo_ons=%.3f (cota %.3f)\n",
                    validas, tasa_iid, tasa_ons, cota_mc);
        chk(validas >= 200 && tasa_iid <= cota_mc && tasa_ons <= cota_mc,
            "N1", "tamano bajo el nulo (>= 200 semillas): las dos pruebas rechazan a tasa <= alfa + su error de Monte Carlo");
    }

    // -------------------------------------------------------------------------
    // N2: independencia de realizaciones (1 hilo vs 4 hilos identicos bit a bit)
    // -------------------------------------------------------------------------
    {
        const auto x_cal = generar_ar1(0.6, 3000, 701);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 250.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d = std::make_unique<Obs>();
        d->calibrar(cd, x_cal);

        auto l1 = std::make_unique<Lil>();
        auto l4 = std::make_unique<Lil>();
        ConfiguracionLilyth c1{}, c4{};
        c1.representacion = c4.representacion = ModoRepresentacionLilyth::CRUDA;
        c1.b_max_nulo = c4.b_max_nulo = 8;
        c1.semilla_nulo = c4.semilla_nulo = 99999;
        c1.hilos_nulo = 1;
        c4.hilos_nulo = 4;
        l1->calibrar(*d, c1, cd, x_cal);
        l4->calibrar(*d, c4, cd, x_cal);

        const auto& f1 = l1->ficha_nulo_cruda();
        const auto& f4 = l4->ficha_nulo_cruda();
        bool id_hilos = (f1.t_max_derivado == f4.t_max_derivado) &&
                        (f1.nulos_sin_memoria.size() == f4.nulos_sin_memoria.size()) &&
                        (f1.nulos_lineal.size() == f4.nulos_lineal.size()) &&
                        !f1.nulos_lineal.empty();
        for (std::size_t i = 0; id_hilos && i < f1.nulos_lineal.size(); ++i) {
            id_hilos = igual(f1.nulos_lineal[i].g_disparo, f4.nulos_lineal[i].g_disparo) &&
                       igual(f1.nulos_lineal[i].dy, f4.nulos_lineal[i].dy) &&
                       igual(f1.nulos_lineal[i].tau_rest_consolidado, f4.nulos_lineal[i].tau_rest_consolidado);
        }
        chk(id_hilos, "N2", "independencia de realizaciones: 1 hilo y 4 hilos producen poblaciones nulas identicas bit a bit");
    }

    // -------------------------------------------------------------------------
    // N3: Regla 7 (cada realizacion parte de una copia nueva del operador delta)
    // -------------------------------------------------------------------------
    {
        std::mt19937_64 rng(801);
        std::normal_distribution<double> nd(0.0, 1.0);
        const std::size_t dia_pasos = 100, n_cal = 4600;
        std::vector<double> xc(n_cal);
        std::vector<std::uint8_t> vc(n_cal, 1);
        std::vector<std::size_t> fc(n_cal);
        std::vector<std::int64_t> dc(n_cal);
        for (std::size_t t = 0; t < n_cal; ++t) {
            fc[t] = (t % dia_pasos) / 25;
            dc[t] = static_cast<std::int64_t>(t / dia_pasos);
            xc[t] = (fc[t] == 0 ? 4.0 : 1.0) * nd(rng);
        }

        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 250.0;
        cd.purga = true;
        cd.bins_fase = 4;
        cd.dias_minimos_purga = 20;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;

        auto d = std::make_unique<Obs>();
        d->calibrar(cd, xc, {}, vc, fc, dc);

        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 6;
        cl.hilos_nulo = 1;
        cl.semilla_nulo = 80808;
        lil->calibrar(*d, cl, cd, xc, vc, fc, dc);

        const auto& fn = lil->ficha_nulo_cruda();
        ConfiguracionLilyth cl_init = cl;
        cl_init.t_max_efectivo = TMAXC;
        Lil lil_plantilla;
        lil_plantilla.configurar(cl_init, cd, *d, xc, vc, fc, dc);

        std::vector<double> g_copia_nueva, g_reutilizada;
        Obs d_reuso = *d;
        for (std::size_t b = 0; b < fn.B; ++b) {
            std::uint64_t sem_b = detalle_nulo::splitmix64((cl.semilla_nulo ^ 0x1111111111111111ULL) ^ (0x9E3779B97F4A7C15ULL * (b + 1)));
            std::mt19937_64 rng_1(sem_b), rng_2(sem_b);
            std::vector<double> x1, x2;
            detalle_nulo::construir_serie_nula(cd.integracion, false, cl.nulo_sin_memoria, fn.segmentos_calma, fn.modelo_lineal, 0.0, n_cal, rng_1, x1);
            detalle_nulo::construir_serie_nula(cd.integracion, false, cl.nulo_sin_memoria, fn.segmentos_calma, fn.modelo_lineal, 0.0, n_cal, rng_2, x2);

            Obs d_nueva = *d;
            d_nueva.reiniciar_flujo();
            d_reuso.reiniciar_flujo();

            Lil l_nueva = lil_plantilla;
            l_nueva.reiniciar_flujo_conservando_purga();
            Lil l_reuso = lil_plantilla;
            l_reuso.reiniciar_flujo_conservando_purga();

            for (std::size_t t = 0; t < n_cal; ++t) {
                const auto s1 = d_nueva.paso(x1[t], NaN, true, fc[t], dc[t]);
                const auto o1 = l_nueva.paso(s1, x1[t], true, NaN, fc[t], dc[t]);
                for (std::size_t i = 0; i < o1.ranuras_cerradas; ++i) {
                    EventoNuloLilyth<TMAXC> ev{};
                    if (detalle_nulo::extraer_evento_nulo<TMAXC>(l_nueva.cerrada(i), false, ev)) g_copia_nueva.push_back(ev.g_disparo);
                }
                const auto s2 = d_reuso.paso(x2[t], NaN, true, fc[t], dc[t]);
                const auto o2 = l_reuso.paso(s2, x2[t], true, NaN, fc[t], dc[t]);
                for (std::size_t i = 0; i < o2.ranuras_cerradas; ++i) {
                    EventoNuloLilyth<TMAXC> ev{};
                    if (detalle_nulo::extraer_evento_nulo<TMAXC>(l_reuso.cerrada(i), false, ev)) g_reutilizada.push_back(ev.g_disparo);
                }
            }
        }
        bool casa_con_copia_nueva = (fn.nulos_sin_memoria.size() == g_copia_nueva.size()) && !g_copia_nueva.empty();
        for (std::size_t i = 0; casa_con_copia_nueva && i < g_copia_nueva.size(); ++i) {
            casa_con_copia_nueva = igual(fn.nulos_sin_memoria[i].g_disparo, g_copia_nueva[i]);
        }
        bool difiere_de_reuso = (g_copia_nueva.size() != g_reutilizada.size());
        for (std::size_t i = 0; !difiere_de_reuso && i < g_copia_nueva.size(); ++i) {
            if (!igual(g_copia_nueva[i], g_reutilizada[i])) difiere_de_reuso = true;
        }
        std::printf("      N3: B=%zu n_nueva=%zu (modulo=%zu) n_reuso=%zu casa=%d difiere=%d\n",
                    fn.B, g_copia_nueva.size(), fn.nulos_sin_memoria.size(), g_reutilizada.size(), int(casa_con_copia_nueva), int(difiere_de_reuso));
        chk(casa_con_copia_nueva && difiere_de_reuso, "N3",
            "Regla 7: cada realizacion parte de una copia nueva del operador delta (= ref bit a bit, distinto de reusar PurgaDiurna)");
    }

    // -------------------------------------------------------------------------
    // S1: replicas de Hawkes en n = {0.0, 0.3, 0.6} sobre AR(1) de calma conocida
    // -------------------------------------------------------------------------
    {
        using ObsH = ObservadorSaltos<24, 4, 120>;
        using LilH = ObservadorLilyth<ObsH, 8, TMAXC, 2048>;

        const double phi = 0.88;
        const std::size_t n_cal = 4096, n_test = 5000;
        const auto x_cal = generar_ar1(phi, n_cal, 901);

        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 800.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d_base = std::make_unique<ObsH>();
        d_base->calibrar(cd, x_cal);

        auto lil_base = std::make_unique<LilH>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 6;
        cl.hilos_nulo = 1;
        cl.forzar_via_b = true;
        lil_base->calibrar(*d_base, cl, cd, x_cal);

        const std::array<double, 3> razones_n{0.0, 0.3, 0.6};
        std::array<double, 3> cif_rep_fin{};
        std::array<double, 3> p_ons_n{};
        std::array<double, 3> p_det_medio_n{};
        bool s1_ok = true;

        for (std::size_t idx_n = 0; idx_n < 3; ++idx_n) {
            const double n_hawkes = razones_n[idx_n];
            std::mt19937_64 rng(910 + idx_n);
            const double sd_e = std::sqrt(1.0 - phi * phi);
            std::normal_distribution<double> nd(0.0, sd_e);
            std::uniform_real_distribution<double> unif(0.0, 1.0);

            std::vector<double> eps(n_test, 0.0);
            for (double& e : eps) e = nd(rng);

            // Madres Poisson espaciadas; cada madre genera con razon de ramificacion n_hawkes hijas en dt = 5..6 pasos
            // (despues de que el espejo de L=4 barras cerro en t_a + 4, pero antes de que el AR(1) phi=0.88 restituya en tau ~ 7..10)
            for (std::size_t t_madre = 300; t_madre + 120 < n_test; t_madre += 300) {
                const int sgn_m = ((t_madre / 300) % 2 == 0) ? +1 : -1;
                eps[t_madre] += sgn_m * 9.5 * sd_e;
                eps[t_madre + 1] += sgn_m * 4.5 * sd_e;
                std::size_t t_actual = t_madre + 1;
                int sgn_actual = sgn_m;
                std::size_t n_hijas = 0;
                while (unif(rng) < n_hawkes && n_hijas < 3) {
                    ++n_hijas;
                    const std::size_t dt_hija = 5;
                    t_actual += dt_hija;
                    if (t_actual + 50 >= n_test) break;
                    // La primera hija corta con signo contrario (replica competitiva); la siguiente continua con mismo signo
                    sgn_actual = (n_hijas == 1) ? -sgn_actual : sgn_actual;
                    eps[t_actual] += sgn_actual * 28.0 * sd_e;
                }
            }

            std::vector<double> x_test(n_test, 0.0);
            double cur = 0.0;
            for (std::size_t t = 0; t < n_test; ++t) {
                cur = phi * cur + eps[t];
                x_test[t] = cur;
            }

            ObsH d_run = *d_base;
            d_run.reiniciar_flujo();
            LilH l_run = *lil_base;
            l_run.reiniciar_flujo_conservando_purga();

            struct Apertura { std::int64_t t, t_cambio; int signo, j; std::size_t l; };
            std::vector<Apertura> cand_ref;
            bool esp_r = false;
            for (std::size_t t = 0; t < n_test; ++t) {
                const auto s = d_run.paso(x_test[t], NaN, true);
                if (!esp_r) {
                    for (std::size_t j = 0; j < s.niveles && j < d_run.escalas(); ++j) {
                        if (!s.nivel[j].alarma.hay) continue;
                        const auto& a = s.nivel[j].alarma;
                        cand_ref.push_back({a.t, a.t_cambio, a.signo, static_cast<int>(j), a.l});
                        esp_r = true;
                        break;
                    }
                }
                for (std::size_t i = 0; i < s.n_eventos; ++i) {
                    if (s.eventos[i].canal == 0) esp_r = false;
                }
                (void)l_run.paso(s, x_test[t], true);
            }
            l_run.completar_nulos();

            // Reconstruir continuaciones y replicas en motor independiente
            std::size_t cont_ref = 0, rep_ref = 0;
            bool viva_ref = false, rest_ref = false;
            int sgn_ref = 0;
            std::int64_t ta_ref = -1;
            double ya_ref = 0.0, dy_ref = 0.0;
            std::size_t ic = 0;
            const std::size_t tm = l_run.configuracion().t_max_efectivo;
            for (std::size_t t = 0; t < n_test; ++t) {
                const bool hay_cand_hoy = (ic < cand_ref.size() && cand_ref[ic].t == static_cast<std::int64_t>(t));
                const bool cand_opuesto = hay_cand_hoy && (cand_ref[ic].signo != sgn_ref);
                if (viva_ref && static_cast<std::int64_t>(t) > ta_ref) {
                    const std::size_t tau = t - static_cast<std::size_t>(ta_ref);
                    if (!cand_opuesto && (x_test[t] - ya_ref) / dy_ref <= 0.5) rest_ref = true;
                    if (tau >= tm) viva_ref = false;
                }
                if (ic < cand_ref.size() && cand_ref[ic].t == static_cast<std::int64_t>(t)) {
                    const auto ap = cand_ref[ic++];
                    bool es_c = false;
                    if (viva_ref) {
                        if (ap.signo == sgn_ref && !rest_ref) { es_c = true; ++cont_ref; }
                        else { ++rep_ref; viva_ref = false; }
                    }
                    if (!es_c) {
                        const auto& f = d_run.ficha_nivel(static_cast<std::size_t>(ap.j));
                        const std::int64_t m = std::max<std::int64_t>(1, f.m);
                        const std::int64_t guarda = static_cast<std::int64_t>(std::ceil(std::max(1.0, f.nulo.tau_int))) * m;
                        const std::int64_t tf = ap.t_cambio - guarda - 1;
                        const std::int64_t ti = tf - 120 * m + 1;
                        if (ti >= 0 && tf >= ti) {
                            std::vector<double> w(x_test.begin() + ti, x_test.begin() + tf + 1);
                            ya_ref = est::cuantil(w, 0.5);
                            dy_ref = x_test[t] - ya_ref;
                            if (dy_ref * ap.signo > 0.0) {
                                viva_ref = true; rest_ref = false; sgn_ref = ap.signo; ta_ref = ap.t;
                            }
                        }
                    }
                }
            }

            const auto& k = l_run.contadores();
            const bool e121_ok = (k.continuaciones == cont_ref && k.replicas == rep_ref);

            const auto& aj = l_run.aalen_johansen_cruda().total;
            const double frac_rep_emp = aj.n_eventos_registrados > 0
                ? static_cast<double>(aj.n_cortados_replica) / static_cast<double>(aj.n_eventos_registrados)
                : 0.0;
            cif_rep_fin[idx_n] = aj.cif_replica[tm - 1];
            const bool cif_ok = std::abs(cif_rep_fin[idx_n] - frac_rep_emp) <= 0.05;

            const auto& fn = l_run.ficha_nulo_cruda();
            p_ons_n[idx_n] = fn.contraste_onsager.p_bilateral;
            const bool pons_ok = fn.contraste_onsager.valido && p_ons_n[idx_n] >= cl.alfa;

            double suma_pdet = 0.0;
            std::size_t n_ev_b = 0, n_completos_b = 0;
            for (const auto& r : l_run.cerradas_historico()) {
                if ((r.cruda.desenlace == DesenlaceRanura::RESTITUIDA || r.cruda.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) &&
                    r.cruda.nulo_onsager.evaluado && r.amplitud_glr != 0.0) {
                    if (r.cruda.nulo_onsager.n_nulos_b >= fn.m_min) {
                        ++n_completos_b;
                        if ((r.t_apertura % 300) <= 5) {
                            suma_pdet += r.cruda.nulo_onsager.p_det_inicial;
                            ++n_ev_b;
                        }
                    }
                }
            }
            p_det_medio_n[idx_n] = n_ev_b > 0 ? suma_pdet / static_cast<double>(n_ev_b) : 0.0;
            std::printf("         [S1 n=%.1f] cont=%zu/%zu rep=%zu/%zu cif=%.3f/%.3f p_ons=%.3f p_det=%.2f comp=%zu/%zu\n",
                        n_hawkes, k.continuaciones, cont_ref, k.replicas, rep_ref, cif_rep_fin[idx_n], frac_rep_emp,
                        p_ons_n[idx_n], p_det_medio_n[idx_n], n_completos_b, n_ev_b);
            if (!e121_ok || !cif_ok || !pons_ok || n_completos_b == 0) s1_ok = false;
        }

        const bool mono_cif = cif_rep_fin[2] > cif_rep_fin[0] && cif_rep_fin[1] > cif_rep_fin[0];
        const bool pdet_estable = std::abs(p_det_medio_n[2] - p_det_medio_n[0]) <= 0.15 && p_det_medio_n[0] >= 0.65;
        std::printf("      S1: CIF_rep(n=0,0.3,0.6)=(%.3f, %.3f, %.3f) p_ons=(%.3f, %.3f, %.3f) p_det=(%.2f, %.2f, %.2f)\n",
                    cif_rep_fin[0], cif_rep_fin[1], cif_rep_fin[2], p_ons_n[0], p_ons_n[1], p_ons_n[2],
                    p_det_medio_n[0], p_det_medio_n[1], p_det_medio_n[2]);
        chk(s1_ok && mono_cif && pdet_estable, "S1",
            "Hawkes n in {0, 0.3, 0.6}: E1.21 = verdad, CIF_replica dentro de banda, Delta_tau_Onsager ~ 0 y p_det estable");
    }

    // -------------------------------------------------------------------------
    // S2: ventanas solapadas por debajo y por encima de W_ancla + guarda + l + T_MAX
    //     y efecto fisico de contaminacion del ancla
    // -------------------------------------------------------------------------
    {
        const double phi = 0.85;
        const auto x_cal = generar_ar1(phi, 4096, 1001);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 2000.0;   // alto para evitar falsas alarmas espontaneas en la prueba de separacion exacta
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d = std::make_unique<Obs>();
        d->calibrar(cd, x_cal);

        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 2;
        cl.hilos_nulo = 1;
        lil->calibrar(*d, cl, cd, x_cal);

        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, d->ficha_nivel(0).nulo.tau_int)));
        const std::size_t l = 4;
        const std::size_t dist_min = w_ancla + guarda + l + TMAXC;   // separacion exacta del borde de solape

        const std::array<double, 8> hist_prev{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const std::array<double, 4> perfil{4.0, 6.0, 7.5, 8.5};
        const std::size_t n_sim = 1200;

        // 1. Por debajo de la ventana (separacion = dist_min - 20): las 2 inyecciones se solapan -> solapes == 2, admitidas == 0
        const std::array<std::int64_t, 2> t_sub{250, static_cast<std::int64_t>(250 + dist_min - 20)};
        const auto r_sub = lil->evaluar_inyecciones_realizacion(false, hist_prev, perfil, +1, w_ancla, guarda, t_sub, n_sim, 1010);

        // 2. En el borde exacto (separacion == dist_min): ini_2 == fin_1 -> solapes == 2, admitidas == 0
        const std::array<std::int64_t, 2> t_borde{250, static_cast<std::int64_t>(250 + dist_min)};
        const auto r_borde = lil->evaluar_inyecciones_realizacion(false, hist_prev, perfil, +1, w_ancla, guarda, t_borde, n_sim, 1011);

        // 3. Por encima de la ventana (separacion = dist_min + 20): 0 solapes y ambas entran al nulo
        const std::array<std::int64_t, 2> t_sup{250, static_cast<std::int64_t>(250 + dist_min + 20)};
        const auto r_sup = lil->evaluar_inyecciones_realizacion(false, hist_prev, perfil, +1, w_ancla, guarda, t_sup, n_sim, 1012);

        // 4. Mundo donde la contaminacion SI cambia el resultado: si no hubiera guarda, una inyeccion previa a distancia 50 pasos
        //    deja su cola en la ventana del ancla Y* de la segunda inyeccion y sesga tau_rest.
        std::vector<double> taus_con_guarda, taus_contaminados;
        for (std::size_t k = 0; k < 20; ++k) {
            const std::array<std::int64_t, 1> t_limpio{300};
            const auto rl = lil->evaluar_inyecciones_realizacion(false, hist_prev, perfil, +1, w_ancla, guarda, t_limpio, 800, 1100 + k);
            for (const auto& ev : rl.eventos_admitidos) {
                if (std::isfinite(ev.tau_rest_consolidado)) taus_con_guarda.push_back(ev.tau_rest_consolidado);
            }
            // Si una segunda inyeccion cae apenas 60 pasos despues de una primera que ya restituyo a mitad pero sigue elevada en la ventana del ancla:
            std::vector<double> perfil_contaminante(50, 5.0);
            std::vector<double> e_cont(50, 0.0);
            innovaciones(lil->ficha_nulo_cruda().modelo_lineal, hist_prev.data(), hist_prev.size(), perfil_contaminante.data(), 50, e_cont.data());
            // Medir el sesgo de Y* cuando el ancla contiene 50 barras a +5.0:
            if (!rl.eventos_admitidos.empty()) {
                const auto& ev0 = rl.eventos_admitidos.front();
                const double y_ancla_sesgada = ev0.y_ancla + 3.5;
                const double dy_sesgado = (ev0.y_ancla + ev0.dy) - y_ancla_sesgada;
                double tau_sesg = double(TMAXC);
                for (std::size_t tau = 1; tau <= ev0.pasos_observados; ++tau) {
                    const double y_t = ev0.y_ancla + ev0.eta[tau - 1] * ev0.dy;
                    if ((y_t - y_ancla_sesgada) / dy_sesgado <= 0.5) {
                        tau_sesg = static_cast<double>(tau);
                        break;
                    }
                }
                taus_contaminados.push_back(tau_sesg);
            }
        }
        const double med_limpia = est::cuantil(taus_con_guarda, 0.5);
        const double med_contam = est::cuantil(taus_contaminados, 0.5);

        std::printf("      S2: dist_min=%zu | sub(sol=%zu, adm=%zu) borde(sol=%zu, adm=%zu) sup(sol=%zu, adm=%zu) | med_limpia=%.1f vs contam=%.1f\n",
                    dist_min, r_sub.solapes, r_sub.eventos_admitidos.size(), r_borde.solapes, r_borde.eventos_admitidos.size(),
                    r_sup.solapes, r_sup.eventos_admitidos.size(), med_limpia, med_contam);

        chk(r_sub.solapes == 2 && r_sub.eventos_admitidos.empty() &&
            r_borde.solapes == 2 && r_borde.eventos_admitidos.empty() &&
            r_sup.solapes == 0 && r_sup.eventos_admitidos.size() == 2 &&
            (med_limpia - med_contam) >= 1.5,
            "S2", "ventanas solapadas: subcritica y borde se descartan (solapes=2, adm=0); supercritica entra (solapes=0, adm=2)");
    }

    // -------------------------------------------------------------------------
    // S3: integridad de la parada de la via (b)
    //     - p_det intermedio cerca del umbral: sigue Binomial Negativa y respeta ceil(m_min / p_det) de las primeras m_min
    //     - p_det = 0: para en m_min intentos, marca sin_nulo_b = true, sin bucle infinito
    // -------------------------------------------------------------------------
    {
        const auto x_cal = generar_ar1(0.5, 4096, 1201);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 1000.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d = std::make_unique<Obs>();
        d->calibrar(cd, x_cal);

        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 2;
        cl.hilos_nulo = 1;
        cl.forzar_via_b = true;
        lil->calibrar(*d, cl, cd, x_cal);

        const auto& fn = lil->ficha_nulo_cruda();
        const std::size_t mm = fn.m_min;   // 8 con alfa = 0.05
        const std::array<double, 4> hist_prev{0.0, 0.0, 0.0, 0.0};

        // 1. Caso p_det = 0 (perfil de amplitud 0.02 sd, muy por debajo del umbral del GLR)
        const std::array<double, 3> perfil_cero{0.01, 0.02, 0.01};
        ReferenciaNuloEvento ref_cero{};
        std::vector<EventoNuloLilyth<TMAXC>> ev_cero;
        std::size_t sol_cero = 0;
        detalle_nulo::ejecutar_via_b_evento<TMAXC>(
            *d, *lil, Integracion::FLUJO, false, fn.modelo_lineal, hist_prev, perfil_cero, +1,
            120, 4, mm, lil->configuracion().t_max_efectivo, 55555, ref_cero, ev_cero, sol_cero);

        const bool pdet_cero_ok = ref_cero.sin_nulo_b && ref_cero.intentos_b == mm &&
                                  ref_cero.n_nulos_b == 0 && ref_cero.p_det_inicial == 0.0 && ev_cero.empty();

        // 2. Caso p_det intermedio (~0.5, perfil cerca del umbral del detector) sobre R = 60 semillas
        const std::array<double, 2> perfil_umbral{3.70, 3.95};
        std::size_t corridas_truncadas_en_tope = 0;
        bool respeta_tope_inicial = true;
        double suma_intentos = 0.0, suma_pdet = 0.0;
        const std::size_t R_S3 = 40;

        for (std::size_t r = 0; r < R_S3; ++r) {
            ReferenciaNuloEvento ref_r{};
            std::vector<EventoNuloLilyth<TMAXC>> ev_r;
            std::size_t sol_r = 0;
            detalle_nulo::ejecutar_via_b_evento<TMAXC>(
                *d, *lil, Integracion::FLUJO, false, fn.modelo_lineal, hist_prev, perfil_umbral, +1,
                120, 4, mm, lil->configuracion().t_max_efectivo, 80000 + r, ref_r, ev_r, sol_r);

            suma_intentos += static_cast<double>(ref_r.intentos_b);
            suma_pdet += ref_r.p_det_inicial;

            // Verificar que p_det_inicial es un multiplo exacto de 1 / mm (calculado en exactamente mm inyecciones)
            const double k0_est = ref_r.p_det_inicial * static_cast<double>(mm);
            if (std::abs(k0_est - std::round(k0_est)) > 1e-9) {
                respeta_tope_inicial = false;
            }
            if (ref_r.p_det_inicial > 0.0) {
                const std::size_t tope_esp = std::min<std::size_t>(
                    MAX_INTENTOS_VIA_B, static_cast<std::size_t>(std::ceil(static_cast<double>(mm) / ref_r.p_det_inicial)));
                if (ref_r.tope_intentos_b != tope_esp || ref_r.intentos_b > tope_esp) {
                    respeta_tope_inicial = false;
                }
                if (ref_r.n_nulos_b < mm && ref_r.intentos_b == tope_esp) {
                    ++corridas_truncadas_en_tope;
                }
            } else {
                if (ref_r.intentos_b != mm || !ref_r.sin_nulo_b) respeta_tope_inicial = false;
            }
        }
        const double p_det_prom = suma_pdet / static_cast<double>(R_S3);
        const double int_prom = suma_intentos / static_cast<double>(R_S3);
        const double media_nb_teor = static_cast<double>(mm) / std::max(0.1, p_det_prom);

        std::printf("      S3: p_det=0 -> intentos=%zu sin_nulo_b=%d | umbral -> p_det_prom=%.2f int_prom=%.1f (NB=%.1f) truncadas=%zu/%zu\n",
                    ref_cero.intentos_b, int(ref_cero.sin_nulo_b), p_det_prom, int_prom, media_nb_teor, corridas_truncadas_en_tope, R_S3);

        chk(pdet_cero_ok && respeta_tope_inicial && p_det_prom > 0.20 && p_det_prom < 0.85 &&
            corridas_truncadas_en_tope >= 3 && std::abs(int_prom - media_nb_teor) <= 5.0,
            "S3", "integridad de la parada: p_det=0 para en m_min con sin_nulo_b; p_det intermedio respeta ceil(m_min/p_det_inicial)");
    }

    std::printf("\n%s: %d fallas\n", fallos ? "FALLA" : "PASA", fallos);
    return fallos ? 1 : 0;
}