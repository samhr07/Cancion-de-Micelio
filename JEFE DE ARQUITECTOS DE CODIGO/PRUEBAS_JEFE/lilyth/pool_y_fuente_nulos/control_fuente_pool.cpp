// control_fuente_pool.cpp -- Ultima fase del Operador de Lilyth (Code, Opus 5.5, 2026-11-16): control con verdad conocida de la FUENTE de
// los nulos (CALMA, VIENTO, COMPLETA) y del POOL DERIVADO. Expectativas selladas en EXPECTATIVAS_POOL_Y_FUENTE.md (sha256 9ece06cd...).
//
// Secciones:
//   pool    : EC6a, EC6b, EC6c (mundo lineal AR(1) phi 0.82, K = 1, M_MIN contra DERIVADO; 400 replicas por caso)
//   fuente  : EC8, pasada 1 del bufer (EC9), EC1-EC4 (mundos M0, M1, M2 de umbral; 40 replicas por mundo)
// Uso: control_fuente_pool.exe hilos [pool|fuente|todo] [replicas_pool] [replicas_fuente]
#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <mutex>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <vector>

namespace al = ark::lilyth;
namespace me = ark::metrologia;
using ObsReal = me::ObservadorSaltos<24, 64, 120>;
using Lil = al::ObservadorLilyth<ObsReal, 8, 64>;
constexpr double ALFA = 0.05;

static int g_fallas = 0;
static std::mutex g_mtx;
static void reporta(bool ok, const char* nombre, const std::string& detalle, bool condicion) {
    std::lock_guard<std::mutex> g(g_mtx);
    std::printf("[%s] %s %s: %s\n", condicion ? (ok ? "PASA" : "FALLA") : (ok ? "SE CUMPLE" : "NO SE CUMPLE"), condicion ? "condicion" : "expectativa",
                nombre, detalle.c_str());
    if (condicion && !ok) ++g_fallas;
}
static std::string f4(double x) { char b[64]; std::snprintf(b, sizeof b, "%.4f", x); return b; }

// Intervalo de Wilson al 95 % para k de n.
static std::pair<double, double> wilson(std::size_t k, std::size_t n) {
    if (n == 0) return {0.0, 1.0};
    const double z = 1.959963984540054, p = double(k) / double(n), d = 1.0 + z * z / n;
    const double c = (p + z * z / (2.0 * n)) / d, h = z * std::sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / d;
    return {std::max(0.0, c - h), std::min(1.0, c + h)};
}
static std::string tasa(std::size_t k, std::size_t n) {
    const auto w = wilson(k, n);
    return std::to_string(k) + "/" + std::to_string(n) + " = " + f4(n ? double(k) / n : 0.0) + " [" + f4(w.first) + ", " + f4(w.second) + "]";
}

static void en_paralelo(std::size_t n, std::size_t hilos, const std::function<void(std::size_t)>& f) {
    std::atomic<std::size_t> sig{0};
    std::vector<std::thread> th;
    for (std::size_t h = 0; h < std::max<std::size_t>(1, hilos); ++h)
        th.emplace_back([&] { for (std::size_t i; (i = sig++) < n;) f(i); });
    for (auto& t : th) t.join();
}

static me::ConfiguracionObservador cfg_delta(double arl0, std::uint64_t sem) {
    me::ConfiguracionObservador c{};
    c.integracion = me::Integracion::FLUJO;
    c.escalas_max = 2;
    c.medir_tau50 = false;
    c.arl0_pasos_nativos = arl0;
    c.semilla = sem;
    return c;
}
static al::ConfiguracionLilyth cfg_lil(al::PoolNuloLilyth pool, std::uint64_t sem) {
    al::ConfiguracionLilyth c{};
    c.representacion = al::ModoRepresentacionLilyth::CRUDA;
    c.nulo_sin_memoria = al::LILYTH_NULO_PERMUTACION;
    c.modo_nulo_lineal = al::ModoNuloLineal::AR_LEVINSON_DURBIN;
    c.modo_armonico = al::ModoArmonico::PRONY_CURVA_MEDIA;
    c.hilos_nulo = 1;
    c.b_mc_resumen = 1999;
    c.semilla_nulo = sem;
    c.contraste_resumen = al::ContrasteResumenLilyth::INTERCAMBIABLE;
    c.pool_nulo = pool;
    return c;
}

// Lectura por serie: "rechaza lenta" en la referencia ref (0 sin memoria, 1 Onsager) = p_holm[ref] < alfa y S > 0 (como E7b).
struct Lectura {
    bool valida = false;
    std::array<bool, 2> rech{}, lenta{};
    std::array<double, 2> p_holm{NaN_(), NaN_()}, S{NaN_(), NaN_()};
    static double NaN_() { return std::numeric_limits<double>::quiet_NaN(); }
};
static Lectura leer(const Lil& lil) {
    Lectura L;
    const auto& fn = lil.ficha_nulo_cruda();
    const auto& h = lil.resumen_armonico_cruda();
    L.valida = h.familia_holm_valida;
    const std::array<const al::ContrasteResumenNulo*, 2> c{&fn.contraste_iid, &fn.contraste_onsager};
    for (int r = 0; r < 2; ++r) {
        L.p_holm[r] = h.p_holm[r];
        L.S[r] = c[r]->estadistico_s;
        L.rech[r] = L.valida && std::isfinite(h.p_holm[r]) && h.p_holm[r] < ALFA;
        L.lenta[r] = L.rech[r] && c[r]->estadistico_s > 0.0;
    }
    return L;
}

// ===================================================================================================================== POOL (EC6)
// Mundo lineal: AR(1) phi 0.82; la tormenta es una innovacion grande (amplitud 12) en t = 200 de un flujo de 400; con `lenta` lleva
// ademas z = 4 * 0.97^tau. Delta y Lilyth calibrados en calma (ARL0 2000 para que K = 1 sea lo comun).
static std::vector<double> ar1(std::size_t n, double phi, std::mt19937_64& g) {
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<double> x(n, 0.0);
    for (std::size_t t = 1; t < n; ++t) x[t] = phi * x[t - 1] + nd(g);
    return x;
}
// DECLARADO TRAS SELLAR (Code): p minimo exacto por enumeracion completa, independiente de la cabecera. Estrato k = {real} U nulos (codigos
// de tau); u de cada miembro = PIT de punto medio contra el resto; S = media de (u - 0.5) eligiendo un miembro por estrato; p_min =
// max(#configuraciones con |S| >= max|S| / total, 1 / (B + 1)). Solo para total de configuraciones <= 2e6.
static double pmin_enumerado(const std::vector<std::vector<double>>& est, std::size_t B, bool& posible) {
    posible = false;
    double total = 1.0;
    for (const auto& e : est) total *= double(e.size());
    if (est.empty() || total > 2e6) return std::numeric_limits<double>::quiet_NaN();
    std::vector<std::vector<double>> u(est.size());
    for (std::size_t k = 0; k < est.size(); ++k)
        for (std::size_t i = 0; i < est[k].size(); ++i) {
            double menor = 0.0, igual = 0.0;
            for (std::size_t j = 0; j < est[k].size(); ++j) if (j != i) { menor += est[k][j] < est[k][i]; igual += est[k][j] == est[k][i]; }
            u[k].push_back((menor + 0.5 * igual) / double(est[k].size() - 1));
        }
    std::vector<double> sumas{0.0};
    for (const auto& uk : u) {
        std::vector<double> nv;
        nv.reserve(sumas.size() * uk.size());
        for (double a : sumas) for (double b : uk) nv.push_back(a + (b - 0.5));
        sumas.swap(nv);
    }
    double mx = 0.0;
    for (double v : sumas) mx = std::max(mx, std::abs(v) / est.size());
    std::size_t c = 0;
    for (double v : sumas) c += std::abs(v) / est.size() >= mx - 1e-12;
    posible = true;
    return std::max(double(c) / double(sumas.size()), 1.0 / double(B + 1));
}

struct ResPool {
    bool ok = false;
    std::size_t K = 0, obj = 0, eventos_eval = 0;
    bool resoluble = false, coherente = true, piso_ok = true, llega = false;
    bool incoherencia_por_empate = false;   // DECLARADO tras sellar: la discrepancia se explica por un empate en el extremo de un estrato
    bool necesaria_ok = true;               // (retirada: mi atajo con el producto tampoco era necesario; ver EC6a'')
    bool pmin_ok = true, pmin_posible = false;
    bool mecanica_ok = true;               // EC6m (anadido antes de los mutantes): objetivo, K, via (b) intentada y completada hasta n*   // EC6a'': p_min de la cabecera = enumeracion independiente; resoluble = (p_min <= alfa/4)
    double pmin = NaN_();
    Lectura L{};
    static double NaN_() { return std::numeric_limits<double>::quiet_NaN(); }
};
// La propagacion de la innovacion: x ya es AR(1); sumar 12 en x[200] sin propagar seria un pico de un paso. Se propaga aparte.
static void propagar(std::vector<double>& x, std::size_t t0, double a, double phi) {
    double v = a;
    for (std::size_t t = t0; t < x.size(); ++t) { x[t] += v; v *= phi; }
}

static ResPool una_pool_bien(std::uint64_t sem, bool lenta, al::PoolNuloLilyth pool) {
    ResPool R;
    std::mt19937_64 g(sem);
    const auto x_cal = ar1(3000, 0.82, g);
    const auto cd = cfg_delta(2000.0, sem ^ 0x1234ULL);
    ObsReal delta;
    delta.calibrar(cd, x_cal);
    if (delta.escalas() == 0) return R;
    Lil lil;
    if (!lil.calibrar(delta, cfg_lil(pool, sem ^ 0xBEEFULL), cd, x_cal) || !lil.calibrada_nulo()) return R;
    auto x = ar1(400, 0.82, g);
    propagar(x, 200, 12.0, 0.82);
    if (lenta) propagar(x, 200, 4.0, 0.97);
    delta.reiniciar_flujo();
    for (std::size_t t = 0; t < x.size(); ++t) (void)lil.paso(delta.paso(x[t]), x[t], true);
    lil.completar_nulos();
    const auto& fn = lil.ficha_nulo_cruda();
    R.ok = true;
    R.K = fn.k_candidatos_pool;
    R.obj = fn.objetivo_pool;
    R.L = leer(lil);
    const auto& c = fn.contraste_onsager;
    R.resoluble = c.resoluble;
    R.pmin = c.p_min_alcanzable;
    double prod = 1.0;
    std::size_t n_eval = 0;
    bool llega = true;
    for (const auto& rn : lil.cerradas_historico()) {
        const auto& ref = rn.cruda.nulo_onsager;
        if (!ref.evaluado) continue;
        ++n_eval;
        prod *= double(ref.n_nulos_usados + 1);
        if (ref.n_nulos_usados < fn.m_min) R.piso_ok = false;
        if (R.obj > 0 && ref.n_nulos_usados < R.obj) llega = false;
    }
    R.eventos_eval = n_eval;
    R.llega = n_eval > 0 && llega;
    if (pool == al::PoolNuloLilyth::DERIVADO) {   // EC6m
        std::size_t kc = 0;
        for (const auto& rn : lil.cerradas_historico())
            kc += rn.cruda.desenlace == al::DesenlaceRanura::RESTITUIDA || rn.cruda.desenlace == al::DesenlaceRanura::CENSURADA_HORIZONTE;
        std::size_t obj = 0;
        if (kc > 0) { obj = 1; while (std::pow(double(obj + 1), double(kc)) < 8.0 / ALFA) ++obj; }   // menor n con (n + 1)^K >= 160
        if (!fn.pool_derivado || fn.k_candidatos_pool != kc || fn.objetivo_pool != obj) R.mecanica_ok = false;
        const std::size_t umb = std::max(fn.m_min, obj);
        for (const auto& rn : lil.cerradas_historico()) {
            const auto& ref = rn.cruda.nulo_onsager;
            if (!ref.evaluado) continue;
            if (ref.n_nulos_usados < umb && ref.intentos_b == 0) R.mecanica_ok = false;   // bajo el objetivo sin intentar la via (b)
            if (ref.nulo_por_inyeccion && std::string(ref.motivo).find("completado") != std::string::npos && ref.n_nulos_usados < umb) R.mecanica_ok = false;
        }
    } else if (fn.pool_derivado || fn.objetivo_pool != 0) {
        R.mecanica_ok = false;
    }
    if (n_eval > 0 && c.valido) {
        const bool cond = prod >= 8.0 / ALFA - 1e-9;
        R.coherente = (cond == c.resoluble);
        R.necesaria_ok = !c.resoluble || cond;
        {   // EC6a'': estratos desde los nulos emparejados publicos y el tau del real
            const auto& emp = lil.nulos_emparejados(false, false);
            const std::size_t tm = lil.configuracion().t_max_efectivo;
            const auto& cer = lil.cerradas_historico();
            std::vector<std::vector<double>> est;
            for (std::size_t k = 0; k < cer.size() && k < emp.size(); ++k) {
                if (!cer[k].cruda.nulo_onsager.evaluado) continue;
                std::vector<double> e{al::detalle_nulo::codificar_tau(cer[k].cruda.tau_rest_consolidado, tm)};
                for (const auto& ev : emp[k]) e.push_back(al::detalle_nulo::codificar_tau(ev.tau_rest_consolidado, tm));
                est.push_back(std::move(e));
            }
            bool pos = false;
            const double pe = pmin_enumerado(est, lil.configuracion().b_mc_resumen, pos);
            R.pmin_posible = pos;
            if (pos) R.pmin_ok = std::abs(pe - c.p_min_alcanzable) <= 1e-12 && (c.resoluble == (c.p_min_alcanzable <= ALFA / 4.0));
        }
        if (!R.coherente) {   // explicacion declarada: el real empata con nulos en su extremo (censura en T_max + 1 o tau iguales)
            const auto& emp = lil.nulos_emparejados(false, false);
            const std::size_t tm = lil.configuracion().t_max_efectivo;
            const auto& cer = lil.cerradas_historico();
            for (std::size_t k = 0; k < cer.size() && k < emp.size(); ++k) {
                if (!cer[k].cruda.nulo_onsager.evaluado) continue;
                const double real = al::detalle_nulo::codificar_tau(cer[k].cruda.tau_rest_consolidado, tm);
                std::size_t empates = 0;
                for (const auto& e : emp[k]) empates += al::detalle_nulo::codificar_tau(e.tau_rest_consolidado, tm) == real;
                double mx = real, mn = real;
                for (const auto& e : emp[k]) { const double v = al::detalle_nulo::codificar_tau(e.tau_rest_consolidado, tm); mx = std::max(mx, v); mn = std::min(mn, v); }
                if (empates > 0 && (real == mx || real == mn)) R.incoherencia_por_empate = true;
            }
        }
    }
    return R;
}

static void seccion_pool(std::size_t hilos, std::size_t reps) {
    std::printf("== POOL DERIVADO (EC6): mundo lineal AR(1) phi 0.82, tormenta 12 en t = 200, %zu replicas por caso ==\n", reps);
    struct Caso { const char* nm; bool lenta; al::PoolNuloLilyth pool; };
    const Caso casos[] = {{"tamano M_MIN", false, al::PoolNuloLilyth::M_MIN}, {"tamano DERIVADO", false, al::PoolNuloLilyth::DERIVADO},
                          {"potencia M_MIN", true, al::PoolNuloLilyth::M_MIN}, {"potencia DERIVADO", true, al::PoolNuloLilyth::DERIVADO}};
    std::vector<std::vector<ResPool>> res(4, std::vector<ResPool>(reps));
    for (int c = 0; c < 4; ++c)
        en_paralelo(reps, hilos, [&](std::size_t i) { res[c][i] = una_pool_bien(9100000ULL + i, casos[c].lenta, casos[c].pool); });
    for (int c = 0; c < 4; ++c) {
        std::size_t ok = 0, val = 0, rech = 0, lenta = 0, k1 = 0, resol = 0, llega = 0, incoh = 0, piso = 0, pmin_mal = 0, incoh_emp = 0, nec_mal = 0, pm_mal = 0, pm_n = 0, mec_mal = 0;
        std::array<std::size_t, 6> hk{};
        for (const auto& r : res[c]) {
            if (!r.ok) continue;
            ++ok;
            ++hk[std::min<std::size_t>(r.K, 5)];
            k1 += r.K == 1;
            if (!r.L.valida) continue;
            ++val;
            rech += r.L.rech[1];
            lenta += r.L.lenta[1];
            resol += r.resoluble;
            llega += r.llega;
            incoh += !r.coherente;
            incoh_emp += !r.coherente && r.incoherencia_por_empate;
            nec_mal += !r.necesaria_ok;
            pm_n += r.pmin_posible;
            mec_mal += !r.mecanica_ok;
            pm_mal += r.pmin_posible && !r.pmin_ok;
            piso += !r.piso_ok;
            if (r.resoluble && !(r.pmin <= ALFA / 4.0 + 1e-15)) ++pmin_mal;
        }
        std::printf("  %-18s ok %zu, familia valida %zu; K: 0:%zu 1:%zu 2:%zu 3:%zu 4:%zu 5+:%zu; resolubles %zu; llegan a n* %zu; rechazo Onsager %s, por lenta %zu\n",
                    casos[c].nm, ok, val, hk[0], hk[1], hk[2], hk[3], hk[4], hk[5], resol, llega, tasa(rech, val).c_str(), lenta);
        reporta(mec_mal == 0 && al::ConfiguracionLilyth{}.pool_nulo == al::PoolNuloLilyth::M_MIN, "EC6m (mecanica, anadida antes de los mutantes)",
                std::string(casos[c].nm) + ": series con objetivo, K o via (b) incoherentes con el diseno A2: " + std::to_string(mec_mal) +
                    "; omision M_MIN: " + (al::ConfiguracionLilyth{}.pool_nulo == al::PoolNuloLilyth::M_MIN ? "si" : "NO"), true);
        if (casos[c].pool == al::PoolNuloLilyth::DERIVADO) {
            reporta(incoh == 0 && piso == 0 && pmin_mal == 0, casos[c].lenta ? "EC6a (potencia)" : "EC6a (tamano)",
                    "series con resoluble distinto de prod(n_k+1) >= 160: " + std::to_string(incoh) + "; con un evento evaluado bajo m_min: " +
                        std::to_string(piso) + "; resolubles con p_min > alfa/4: " + std::to_string(pmin_mal), true);
            reporta(incoh == incoh_emp && nec_mal == 0, casos[c].lenta ? "EC6a' DECLARADO TRAS SELLAR, RETIRADO (potencia)" : "EC6a' DECLARADO TRAS SELLAR, RETIRADO (tamano)",
                    "de las " + std::to_string(incoh) + " discrepancias, explicadas por un empate del real en su extremo: " + std::to_string(incoh_emp) +
                        "; resolubles con prod < 160: " + std::to_string(nec_mal) + " (el atajo tambien era falso: con extremos asimetricos p_min = 1/prod)", false);
            reporta(pm_mal == 0 && pm_n > 0, casos[c].lenta ? "EC6a'' DECLARADO TRAS SELLAR (potencia)" : "EC6a'' DECLARADO TRAS SELLAR (tamano)",
                    "p_min de la cabecera = enumeracion completa independiente (<= 1e-12) y resoluble = (p_min <= alfa/4): fallan " + std::to_string(pm_mal) + " de " +
                        std::to_string(pm_n) + " series enumerables", true);
        }
    }
    {
        std::size_t val = 0, rech = 0;
        for (const auto& r : res[1]) if (r.ok && r.L.valida) { ++val; rech += r.L.rech[1]; }
        reporta(wilson(rech, val).first <= ALFA, "EC6b", "DERIVADO, tamano en la referencia Onsager: " + tasa(rech, val) + " (Wilson inferior <= 0.05)", true);
    }
    {
        std::size_t vd = 0, ld = 0, vm = 0, lm = 0;
        for (const auto& r : res[3]) if (r.ok && r.L.valida) { ++vd; ld += r.L.lenta[1]; }
        for (const auto& r : res[2]) if (r.ok && r.L.valida) { ++vm; lm += r.L.lenta[1]; }
        reporta(vd > 0 && 2 * ld >= vd && lm == 0, "EC6c", "potencia por lenta, DERIVADO " + tasa(ld, vd) + " (>= 50 %), M_MIN " + tasa(lm, vm) + " (0 %)", false);
    }
}

// ===================================================================================================================== FUENTE
// Proceso de umbral: y_t = phi(|y_{t-1}|) y_{t-1} + e_t, phi = 0.5 si |y| < 4, 0.92 si no. Tormenta en t0: e_{t0} += 12 * signo. Mecanismo
// lento: x = y + z, z = 4 * 0.97^tau desde t0 (no entra en el umbral).
struct Tormenta { std::size_t t0; int signo; bool z; };
static std::vector<double> umbral(std::size_t n, const std::vector<Tormenta>& tor, std::mt19937_64& g) {
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<double> y(n, 0.0), x(n, 0.0);
    std::vector<double> shock(n, 0.0);
    for (const auto& t : tor) if (t.t0 < n) shock[t.t0] += 12.0 * t.signo;
    for (std::size_t t = 1; t < n; ++t) {
        const double phi = std::abs(y[t - 1]) < 4.0 ? 0.5 : 0.92;
        y[t] = phi * y[t - 1] + nd(g) + shock[t];
    }
    for (std::size_t t = 0; t < n; ++t) x[t] = y[t];
    for (const auto& t : tor)
        if (t.z)
            for (std::size_t s = t.t0; s < n && s < t.t0 + 400; ++s) x[s] += 4.0 * t.signo * std::pow(0.97, double(s - t.t0));
    return x;
}
// Catalogo: tormenta cada `cada` desde `ini`, signo alterno; mundo 0: ninguna con z; 1: todas; 2: una de cada cuatro.
static std::vector<Tormenta> catalogo(std::size_t n, std::size_t ini, std::size_t cada, int mundo) {
    std::vector<Tormenta> v;
    std::size_t k = 0;
    for (std::size_t t = ini; t + 60 < n; t += cada, ++k) v.push_back({t, (k % 2 == 0) ? 1 : -1, mundo == 1 || (mundo == 2 && k % 4 == 1)});
    return v;
}

constexpr std::size_t N_CAL_DELTA = 3000, N_HIST = 6000, CADA_HIST = 300, N_OP = 3200, CADA_OP = 200;
constexpr double ARL0_FUENTE = 250.0;

// Binomial unilateral P(X >= k), X ~ Bin(n, p).
static double binom_sup(std::size_t k, std::size_t n, double p) {
    if (k == 0) return 1.0;
    double s = 0.0;
    for (std::size_t i = k; i <= n; ++i)
        s += std::exp(std::lgamma(n + 1.0) - std::lgamma(i + 1.0) - std::lgamma(n - i + 1.0) + i * std::log(p) + (n - i) * std::log1p(-p));
    return std::min(1.0, s);
}

struct Bufer { std::size_t pre = 0, post = 0; bool post_censurado = false; std::size_t n_tormentas = 0, n_medidas = 0; double tasa_calma = 0.0; };

// Pasada 1 (nota del bufer): Lilyth de la Etapa 1 (configurar, horizonte T_MAX, sin nulos) sobre la historia de calibracion como flujo.
static Bufer pasada1(int mundo, std::size_t reps, std::size_t hilos) {
    constexpr std::size_t TM = 64;
    std::mutex m;
    std::vector<double> taus;
    std::vector<std::size_t> cuenta_lag(TM + 1, 0);
    std::size_t n_tor = 0, aper_calma = 0, pasos_calma = 0;
    en_paralelo(reps, hilos, [&](std::size_t r) {
        std::mt19937_64 g(7700000ULL + 1000ULL * mundo + r);
        const auto calma = umbral(N_CAL_DELTA, {}, g);
        const auto cd = cfg_delta(ARL0_FUENTE, 5500ULL + r);
        ObsReal delta;
        delta.calibrar(cd, calma);
        if (delta.escalas() == 0) return;
        const auto tor = catalogo(N_HIST, 150, CADA_HIST, mundo);
        const auto hist = umbral(N_HIST, tor, g);
        auto cl = cfg_lil(al::PoolNuloLilyth::M_MIN, 1ULL);
        cl.t_max_efectivo = TM;
        Lil lil;
        if (!lil.configurar(cl, cd, delta, calma)) return;
        delta.reiniciar_flujo();
        for (std::size_t t = 0; t < hist.size(); ++t) (void)lil.paso(delta.paso(hist[t]), hist[t], true);
        // tasa de apertura en calma: la misma Lilyth sobre la serie de calma
        Lil lc;
        (void)lc.configurar(cl, cd, delta, calma);
        delta.reiniciar_flujo();
        std::mt19937_64 g2(8800000ULL + r);
        const auto calma2 = umbral(N_HIST, {}, g2);
        for (std::size_t t = 0; t < calma2.size(); ++t) (void)lc.paso(delta.paso(calma2[t]), calma2[t], true);
        std::vector<double> tl;
        std::vector<std::size_t> cl_lag(TM + 1, 0);
        for (const auto& tt : tor) {
            bool medido = false;
            for (const auto& rn : lil.cerradas_historico()) {
                if (rn.t_apertura < 0) continue;
                const auto ta = std::size_t(rn.t_apertura);
                if (!medido && ta >= tt.t0 && ta <= tt.t0 + 5) {
                    const auto d = rn.cruda.desenlace;
                    if (d == al::DesenlaceRanura::RESTITUIDA) { tl.push_back(rn.cruda.tau_rest_consolidado); medido = true; }
                    else if (d == al::DesenlaceRanura::CENSURADA_HORIZONTE) { tl.push_back(double(TM) + 1.0); medido = true; }
                }
                if (ta < tt.t0 && tt.t0 - ta <= TM) cl_lag[tt.t0 - ta] = 1;   // apertura en el rezago l (una por tormenta y rezago)
            }
            std::lock_guard<std::mutex> gg(m);
            for (std::size_t l = 1; l <= TM; ++l) cuenta_lag[l] += cl_lag[l], cl_lag[l] = 0;
        }
        std::lock_guard<std::mutex> gg(m);
        n_tor += tor.size();
        taus.insert(taus.end(), tl.begin(), tl.end());
        aper_calma += lc.contadores().aperturas;
        pasos_calma += calma2.size();
    });
    Bufer B;
    B.n_tormentas = n_tor;
    B.n_medidas = taus.size();
    B.tasa_calma = pasos_calma ? double(aper_calma) / double(pasos_calma) : 0.0;
    if (!taus.empty()) {
        std::sort(taus.begin(), taus.end());
        const std::size_t i = std::min(taus.size() - 1, std::size_t(std::ceil((1.0 - ALFA) * taus.size())) - 1);
        const double q = taus[i];
        B.post_censurado = q > double(TM);
        B.post = B.post_censurado ? TM : std::size_t(std::ceil(q));
    }
    // bufer previo: binomial unilateral por rezago contra 1/ARL0, Holm sobre los rezagos
    std::vector<std::pair<double, std::size_t>> ps;
    for (std::size_t l = 1; l <= TM; ++l) ps.push_back({binom_sup(cuenta_lag[l], n_tor, 1.0 / ARL0_FUENTE), l});
    std::sort(ps.begin(), ps.end());
    for (std::size_t i = 0; i < ps.size(); ++i) {
        if (ps[i].first * double(ps.size() - i) >= ALFA) break;
        B.pre = std::max(B.pre, ps[i].second);
    }
    return B;
}

struct ResFuente {
    bool ok = false;
    std::array<Lectura, 3> L{};                       // CALMA, VIENTO, COMPLETA
    std::array<std::vector<double>, 3> u_sin_z{}, u_con_z{};
    std::array<std::size_t, 3> n_fuente{}, tmax{};
    bool ec8_n = true, ec8_ident = true, ec8_disjunto = true, ec8_seg = true;
};

static ResFuente una_fuente(int mundo, std::size_t r, const Bufer& B, bool verificar_ident) {
    ResFuente R;
    std::mt19937_64 g(6600000ULL + 1000ULL * mundo + r);
    const auto calma = umbral(N_CAL_DELTA, {}, g);
    const auto cd = cfg_delta(ARL0_FUENTE, 4400ULL + r);
    ObsReal delta;
    delta.calibrar(cd, calma);
    if (delta.escalas() == 0) return R;
    const auto tor_h = catalogo(N_HIST, 150, CADA_HIST, mundo);
    const auto hist = umbral(N_HIST, tor_h, g);
    const auto tor_o = catalogo(N_OP, 200, CADA_OP, mundo);
    const auto op = umbral(N_OP, tor_o, g);
    // mascaras: VIENTO = [t0 - pre, t0 + post]; CALMA = su complemento; COMPLETA = vacia
    std::vector<std::uint8_t> viento(N_HIST, 0), calm(N_HIST, 1);
    for (const auto& t : tor_h)
        for (std::size_t s = (t.t0 >= B.pre ? t.t0 - B.pre : 0); s <= t.t0 + B.post && s < N_HIST; ++s) viento[s] = 1, calm[s] = 0;
    const std::array<std::span<const std::uint8_t>, 3> mask{std::span<const std::uint8_t>(calm), std::span<const std::uint8_t>(viento), {}};
    // D2: la historia de calibracion y el flujo de operacion son series distintas (disjuntas por construccion).
    R.ec8_disjunto = (hist.data() != op.data());
    const std::uint64_t sem_l = 3300000ULL + 1000ULL * mundo + r;
    for (int f = 0; f < 3; ++f) {
        Lil lil;
        if (!lil.calibrar(delta, cfg_lil(al::PoolNuloLilyth::M_MIN, sem_l), cd, std::span<const double>(hist), {}, {}, {}, {}, mask[f]) ||
            !lil.calibrada_nulo())
            return R;
        delta.reiniciar_flujo();
        for (std::size_t t = 0; t < op.size(); ++t) (void)lil.paso(delta.paso(op[t]), op[t], true);
        lil.completar_nulos();
        R.L[f] = leer(lil);
        const auto& fn = lil.ficha_nulo_cruda();
        R.n_fuente[f] = fn.n_fuente;
        R.tmax[f] = lil.configuracion().t_max_efectivo;
        std::size_t esperado = 0;
        for (std::size_t i = 0; i < N_HIST; ++i) esperado += (mask[f].empty() || mask[f][i]) && std::isfinite(hist[i]);
        if (fn.n_fuente != esperado || fn.fuente_externa != !mask[f].empty()) R.ec8_n = false;
        {   // segmentos esperados: tramos maximales de 1 en (mascara y valido); anadido antes de los mutantes
            std::size_t seg = 0;
            bool dentro = false;
            for (std::size_t i = 0; i < N_HIST; ++i) {
                const bool en = (mask[f].empty() || mask[f][i]) && std::isfinite(hist[i]);
                seg += en && !dentro;
                dentro = en;
            }
            if (fn.segmentos_calma.size() != seg) R.ec8_seg = false;
        }
        for (const auto& tt : tor_o)
            for (const auto& rn : lil.cerradas_historico()) {
                if (rn.t_apertura < 0) continue;
                const auto ta = std::size_t(rn.t_apertura);
                if (ta >= tt.t0 && ta <= tt.t0 + 5 && rn.cruda.nulo_onsager.evaluado) {
                    (tt.z ? R.u_con_z : R.u_sin_z)[f].push_back(rn.cruda.nulo_onsager.u_pit);
                    break;
                }
            }
        if (verificar_ident && f == 2) {   // EC8: mascara toda en 1 = mascara vacia, bit a bit en lo que se lee
            std::vector<std::uint8_t> unos(N_HIST, 1);
            Lil l1;
            if (!l1.calibrar(delta, cfg_lil(al::PoolNuloLilyth::M_MIN, sem_l), cd, std::span<const double>(hist), {}, {}, {}, {},
                             std::span<const std::uint8_t>(unos))) { R.ec8_ident = false; }
            else {
                delta.reiniciar_flujo();
                for (std::size_t t = 0; t < op.size(); ++t) (void)l1.paso(delta.paso(op[t]), op[t], true);
                l1.completar_nulos();
                const auto a = lil.resumen_armonico_cruda(), b = l1.resumen_armonico_cruda();
                const auto& fa = lil.ficha_nulo_cruda();
                const auto& fb = l1.ficha_nulo_cruda();
                auto igual = [](double x, double y) { return (std::isnan(x) && std::isnan(y)) || std::memcmp(&x, &y, sizeof x) == 0; };
                for (int k = 0; k < 4; ++k) R.ec8_ident = R.ec8_ident && igual(a.p_crudos[k], b.p_crudos[k]) && igual(a.p_holm[k], b.p_holm[k]);
                R.ec8_ident = R.ec8_ident && igual(fa.contraste_onsager.estadistico_s, fb.contraste_onsager.estadistico_s) &&
                              igual(fa.contraste_iid.estadistico_s, fb.contraste_iid.estadistico_s) && fa.B == fb.B &&
                              fa.t_max_derivado == fb.t_max_derivado && fa.n_fuente == fb.n_fuente;
            }
        }
    }
    R.ok = true;
    return R;
}

static double media(const std::vector<double>& v) { return v.empty() ? std::numeric_limits<double>::quiet_NaN() : std::accumulate(v.begin(), v.end(), 0.0) / v.size(); }

static void seccion_fuente(std::size_t hilos, std::size_t reps) {
    const char* NM[3] = {"M0 solo agrupamiento", "M1 mecanismo universal", "M2 mecanismo en 1 de 4"};
    const char* FU[3] = {"CALMA", "VIENTO", "COMPLETA"};
    std::array<Bufer, 3> buf{};
    std::array<std::array<std::size_t, 3>, 3> lenta_ons{}, val_ons{}, lenta_iid{};
    std::array<std::array<std::vector<double>, 3>, 3> usz{}, ucz{};
    bool ec8n = true, ec8i = true, ec8d = true, ec8s = true;
    for (int w = 0; w < 3; ++w) {
        buf[w] = pasada1(w, 5, hilos);
        std::printf("== %s: pasada 1 del bufer: %zu tormentas, %zu medidas; posterior %zu%s, previo %zu; tasa de apertura en calma %s por paso (1/ARL0 = %s) ==\n",
                    NM[w], buf[w].n_tormentas, buf[w].n_medidas, buf[w].post, buf[w].post_censurado ? " (CENSURADO en T_MAX)" : "", buf[w].pre,
                    f4(buf[w].tasa_calma).c_str(), f4(1.0 / ARL0_FUENTE).c_str());
        std::vector<ResFuente> res(reps);
        en_paralelo(reps, hilos, [&](std::size_t r) { res[r] = una_fuente(w, r, buf[w], r < 3); });
        std::array<std::size_t, 3> tm_sum{}, nf{};
        std::size_t ok = 0;
        for (const auto& R : res) {
            if (!R.ok) continue;
            ++ok;
            ec8n = ec8n && R.ec8_n;
            ec8i = ec8i && R.ec8_ident;
            ec8d = ec8d && R.ec8_disjunto;
            ec8s = ec8s && R.ec8_seg;
            for (int f = 0; f < 3; ++f) {
                if (R.L[f].valida) { ++val_ons[w][f]; lenta_ons[w][f] += R.L[f].lenta[1]; lenta_iid[w][f] += R.L[f].lenta[0]; }
                usz[w][f].insert(usz[w][f].end(), R.u_sin_z[f].begin(), R.u_sin_z[f].end());
                ucz[w][f].insert(ucz[w][f].end(), R.u_con_z[f].begin(), R.u_con_z[f].end());
                tm_sum[f] += R.tmax[f];
                nf[f] += R.n_fuente[f];
            }
        }
        for (int f = 0; f < 3; ++f)
            std::printf("  %-9s T_max medio %s, n_fuente medio %s; Onsager rechaza por lenta %s; sin memoria por lenta %zu; u_pit Onsager medio: sin z %s (n %zu), con z %s (n %zu)\n",
                        FU[f], f4(ok ? double(tm_sum[f]) / ok : 0).c_str(), f4(ok ? double(nf[f]) / ok : 0).c_str(), tasa(lenta_ons[w][f], val_ons[w][f]).c_str(),
                        lenta_iid[w][f], f4(media(usz[w][f])).c_str(), usz[w][f].size(), f4(media(ucz[w][f])).c_str(), ucz[w][f].size());
    }
    reporta(ec8n && ec8i && ec8d && ec8s, "EC8", std::string("n_fuente = muestras validas en la mascara: ") + (ec8n ? "si" : "NO") + "; mascara de unos = vacia bit a bit: " +
                                               (ec8i ? "si" : "NO") + "; historia y operacion disjuntas: " + (ec8d ? "si" : "NO") +
                                               "; segmentos = tramos de la mascara (anadido antes de los mutantes): " + (ec8s ? "si" : "NO"), true);
    auto wl = [&](int w, int f) { return wilson(lenta_ons[w][f], val_ons[w][f]); };
    auto fr = [&](int w, int f) { return val_ons[w][f] ? double(lenta_ons[w][f]) / val_ons[w][f] : 0.0; };
    reporta(wl(0, 1).first <= ALFA, "EC1", "M0 VIENTO " + tasa(lenta_ons[0][1], val_ons[0][1]) + " (Wilson inferior <= 0.05)", false);
    reporta(wl(0, 0).first > ALFA, "EC2", "M0 CALMA " + tasa(lenta_ons[0][0], val_ons[0][0]) + " (Wilson inferior > 0.05)", false);
    reporta(fr(0, 2) <= fr(0, 0) && fr(0, 2) >= fr(0, 1), "EC1b", "M0 COMPLETA " + f4(fr(0, 2)) + " entre VIENTO " + f4(fr(0, 1)) + " y CALMA " + f4(fr(0, 0)), false);
    reporta(fr(1, 0) >= 0.8 && fr(1, 2) >= 0.8 && wl(1, 1).first <= ALFA, "EC3",
            "M1 CALMA " + f4(fr(1, 0)) + ", COMPLETA " + f4(fr(1, 2)) + " (>= 0.80); VIENTO " + tasa(lenta_ons[1][1], val_ons[1][1]) + " (Wilson inferior <= 0.05)", false);
    const double v_sz = media(usz[2][1]), v_cz = media(ucz[2][1]), c_sz = media(usz[2][0]);
    reporta(v_sz >= 0.4 && v_sz <= 0.6 && v_cz > v_sz && c_sz > 0.6, "EC4",
            "M2 u_pit Onsager medio: VIENTO sin z " + f4(v_sz) + " (en [0.4, 0.6]), con z " + f4(v_cz) + " (mayor); CALMA sin z " + f4(c_sz) + " (> 0.6)", false);
    reporta(buf[0].pre == 0 && buf[1].pre == 0 && buf[2].pre == 0 && buf[1].post > buf[0].post, "EC9",
            "bufer previo M0/M1/M2 = " + std::to_string(buf[0].pre) + "/" + std::to_string(buf[1].pre) + "/" + std::to_string(buf[2].pre) + " (0); posterior M1 " +
                std::to_string(buf[1].post) + " > M0 " + std::to_string(buf[0].post), false);
}

#ifndef CONTROL_SIN_MAIN
int main(int argc, char** argv) {
    const std::size_t hilos = argc > 1 ? std::stoul(argv[1]) : 8;
    const std::string que = argc > 2 ? argv[2] : "todo";
    const std::size_t rp = argc > 3 ? std::stoul(argv[3]) : 400;
    const std::size_t rf = argc > 4 ? std::stoul(argv[4]) : 40;
    if (que == "pool" || que == "todo") seccion_pool(hilos, rp);
    if (que == "fuente" || que == "todo") seccion_fuente(hilos, rf);
    std::printf("FALLAS DE CONDICIONES: %d\n", g_fallas);
    return g_fallas ? 1 : 0;
}
#endif
