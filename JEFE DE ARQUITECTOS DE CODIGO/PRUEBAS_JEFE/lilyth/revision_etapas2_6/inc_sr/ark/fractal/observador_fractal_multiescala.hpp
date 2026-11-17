// observador_fractal_multiescala.hpp -- Capa 2 y piramide del Bloque Fractal Universal (ARK-FRAC-01), C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-03.
//
// ObservadorFractalUniversal<T, W, S, N_REG>: ingiere CUALQUIER senal causal x_t (precio, OFI, volumen, profundidad L3, spread,
// volatilidad; T aritmetico), calibra su propio nulo y emite en cada paso del sensor el estado del borde presente:
//   nu_t = (H^_t - mu0(X)) / s0(X),  CUSUM bilateral (k, h_crit),  norma L2 o L1 elegida por la propia innovacion,
//   H regularizado = x[N-1] del resolvente (Morozov L2 o Condat L1),  tau50(X) medido.
//
// Decisiones del Comandante (2026-11-03, acta JEFE DE ARQUITECTOS DE CODIGO/AUDITORIA_CAUSAL/BLOQUE_FRACTAL_UNIVERSAL_CXX.md):
//   lambda_TV(X) = cuantil 1 - alfa de max|S_k| (puente de sumas parciales de H^) en ventanas de N_REG pasos del nulo de X;
//                  Donoho-Johnstone s0 sqrt(2 ln N) se exporta solo como telemetria (parte el reposo en escalones);
//   h_crit       = el que da el ARL0 objetivo sobre la innovacion del nulo de X (Siegmund iid queda como referencia: con
//                  ACF1 ~ 0.94 da ARL0 realizados 10-60 veces menores que el nominal);
//   k            = Delta/2, Delta = ganancia del sensor para el dH objetivo (medida por inyeccion) / s0;
//   lambda_L2    = Morozov del Fisico sobre s0.
//
// Regla de la norma, por omision (decision del Comandante, 2026-11-03): CERTIFICADO_PUENTE, L1 si la ventana de H rechaza la
// meseta al nivel alfa de su propio nulo (max|S_k| > lambda_TV), L2 en otro caso. La literal de la doctrina (SORPRESA_CUSUM: L1
// mientras S+ > 0 o S- > 0) queda disponible; con ACF1 ~ 0.94 da L1 el 96 % del reposo (acta §5).
//
// Ciclo caliente (paso): cero reservas de memoria (prueba C5). La calibracion es fuera de linea y si reserva (std::vector).
//
// DECLARACION METROLOGICA (EL ABISMO - COTA I & II):
// Observador de estado metrologico del medio. Queda terminantemente prohibido
// usar las salidas de H en escalas < 1h como gatillos reactivos directos taker.
// Su funcion es la telemetria del espacio de fases, la modulacion de matrices de covarianza
// y la calibracion de Kelly en horizontes macro (> 4 horas).
#pragma once

#include "ark/fractal/regularizador_dual_fractal.hpp"
#include "ark/fractal/sensor_dfa_causal.hpp"
#include "ark/fractal/tipos_fractal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace ark::fractal {

// Regla que elige la norma del borde presente (las dos se gobiernan por la propia senal, sin variables externas).
//  SORPRESA_CUSUM (literal de la doctrina §3.C): L1 mientras S+ > 0 o S- > 0 (|nu| > k sostenido o frente de choque).
//  CERTIFICADO_PUENTE: L1 si la ventana de H rechaza la meseta al nivel alfa, max|S_k| > lambda_TV (el mismo certificado dual
//                      que ancla lambda_TV); en reposo vale L1 en una fraccion alfa de los pasos por construccion.
enum class ReglaNorma : std::uint8_t { SORPRESA_CUSUM, CERTIFICADO_PUENTE };

struct ConfiguracionFractal {
    Integracion integracion = Integracion::NIVEL;
    ReglaNorma regla_norma = ReglaNorma::CERTIFICADO_PUENTE;   // decision del Comandante, 2026-11-03 (acta §5)
    NuloCalibracion nulo = NuloCalibracion::HISTORICO;
    Agregacion agregacion = Agregacion::ULTIMO;
    double delta_h_objetivo = 0.25;        // salto de H que el detector debe ver (EL ABISMO §1)
    double arl0_objetivo_barras = 0.0;     // OBLIGATORIO: falsas alarmas esperadas una vez cada ARL0 barras de la escala
    double alfa_lambda = 0.05;             // nivel nominal del proyecto: fraccion de ventanas en reposo que la TV puede partir
    std::uint64_t semilla = 20261103;
    std::size_t replicas_nulo_signos = 8;  // mundos de signos para el nulo SIGNOS y para la ganancia
    std::size_t replicas_salto = 64;       // saltos inyectados para tau50
    std::size_t horizonte_salto_barras = 0;  // 0 = 20*W barras tras el salto
    std::size_t alarmas_minimas_arl = 30;  // ARL(h) se lee con >= 30 alarmas; por encima se extrapola log-lineal (declarado)
    ConfiguracionDFA dfa{};
};

namespace detalle_obs {

inline double cuantil_lineal(std::vector<double> v, double q) {   // numpy.quantile(method="linear")
    if (v.empty()) return NaN;
    std::sort(v.begin(), v.end());
    const double pos = q * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    return v[lo] + (pos - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

// FFT radix 2 (Cooley-Tukey iterativo), n potencia de 2.
inline void fft(std::vector<std::complex<double>>& a) {
    const std::size_t n = a.size();
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    const double pi = 3.14159265358979323846;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double ang = -2.0 * pi / static_cast<double>(len);
        const std::complex<double> wl(std::cos(ang), std::sin(ang));
        for (std::size_t i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (std::size_t j = 0; j < len / 2; ++j) {
                const auto u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wl;
            }
        }
    }
}

// Ruido gaussiano fraccionario de varianza 1 (Davies-Harte), longitud n (se genera la potencia de 2 siguiente y se trunca).
inline std::vector<double> fgn(std::size_t n, double H, std::mt19937_64& rng) {
    std::normal_distribution<double> N01(0.0, 1.0);
    std::vector<double> out(n);
    if (std::abs(H - 0.5) < 1e-12) { for (auto& v : out) v = N01(rng); return out; }
    std::size_t m = 1;
    while (m < n) m <<= 1;
    const std::size_t L = 2 * m;
    std::vector<double> g(m + 1);
    for (std::size_t k = 0; k <= m; ++k) {
        const double kk = static_cast<double>(k);
        g[k] = 0.5 * (std::pow(kk + 1.0, 2 * H) - 2.0 * std::pow(kk, 2 * H) + std::pow(std::abs(kk - 1.0), 2 * H));
    }
    std::vector<std::complex<double>> c(L);
    for (std::size_t k = 0; k <= m; ++k) c[k] = g[k];
    for (std::size_t k = m + 1; k < L; ++k) c[k] = g[L - k];
    fft(c);
    std::vector<std::complex<double>> w(L);
    for (std::size_t k = 0; k < L; ++k) {
        const double lam = std::max(0.0, c[k].real());
        w[k] = std::sqrt(lam / (2.0 * static_cast<double>(L))) * std::complex<double>(N01(rng), N01(rng));
    }
    fft(w);
    for (std::size_t i = 0; i < n; ++i) out[i] = w[i].real() * std::sqrt(2.0);
    return out;
}

// CUSUM bilateral con reinicio de los dos lados tras cada alarma. Devuelve el numero de alarmas.
inline std::size_t alarmas_cusum(const std::vector<double>& nu, double k, double h) {
    double sp = 0.0, sn = 0.0;
    std::size_t al = 0;
    for (double v : nu) {
        sp = std::max(0.0, sp + v - k);
        sn = std::max(0.0, sn - v - k);
        if (sp >= h || sn >= h) { ++al; sp = sn = 0.0; }
    }
    return al;
}

struct ResultadoH { double h = NaN; bool extrapolado = false; std::size_t alarmas = 0; };

// h tal que ARL(h) = arl0 (pasos) sobre los segmentos nulos `segs` (cada uno empieza en reposo).
inline ResultadoH h_por_arl(const std::vector<std::vector<double>>& segs, double k, double arl0, std::size_t min_alarmas) {
    std::size_t total = 0;
    for (const auto& s : segs) total += s.size();
    auto arl = [&](double h, std::size_t& al) {
        al = 0;
        for (const auto& s : segs) al += alarmas_cusum(s, k, h);
        return al ? static_cast<double>(total) / static_cast<double>(al) : std::numeric_limits<double>::infinity();
    };
    std::vector<double> hs, la;
    std::vector<std::size_t> als;
    for (int j = 1; j <= 4000; ++j) {        // rejilla h = 0.05 j
        const double h = 0.05 * j;
        std::size_t al = 0;
        const double a = arl(h, al);
        if (al < min_alarmas) break;
        hs.push_back(h); la.push_back(std::log(a)); als.push_back(al);
        if (a >= arl0) break;
    }
    ResultadoH r{};
    if (hs.empty()) return r;
    const double objetivo = std::log(arl0);
    if (la.back() >= objetivo) {
        std::size_t i = la.size() - 1;
        if (i == 0) { r.h = hs[0]; r.alarmas = als[0]; return r; }
        const double t = (objetivo - la[i - 1]) / (la[i] - la[i - 1]);
        r.h = hs[i - 1] + t * (hs[i] - hs[i - 1]);
        std::size_t al = 0; arl(r.h, al); r.alarmas = al;
        return r;
    }
    // extrapolacion log-lineal (declarada en la calibracion) con los puntos de la mitad superior del rango observable de h:
    // los ultimos puntos de la rejilla cubren 0.05 de h y suelen tener el mismo recuento de alarmas (pendiente nula)
    const double h_medio = 0.5 * hs.back();
    std::size_t i0 = 0;
    while (i0 + 1 < hs.size() && hs[i0] < h_medio) ++i0;
    const std::size_t m = hs.size() - i0;
    if (m < 2) return r;
    double mx = 0, my = 0;
    for (std::size_t i = i0; i < hs.size(); ++i) { mx += hs[i]; my += la[i]; }
    mx /= static_cast<double>(m); my /= static_cast<double>(m);
    double sxx = 0, sxy = 0;
    for (std::size_t i = i0; i < hs.size(); ++i) { sxx += (hs[i] - mx) * (hs[i] - mx); sxy += (hs[i] - mx) * (la[i] - my); }
    if (sxx <= 0 || sxy <= 0) return r;
    const double b = sxy / sxx;
    r.h = hs.back() + (objetivo - la.back()) / b;
    r.extrapolado = true;
    std::size_t al = 0; arl(r.h, al); r.alarmas = al;
    return r;
}

// Referencia iid de operador_decision.hpp (calcular_h_siegmund), transcrita. NO se usa para decidir.
inline double h_siegmund(double arl0, double k) {
    if (arl0 <= 1.0 || k <= 1e-6) return 5.0;
    const double target = 1.0 + 2.0 * k * k * arl0;
    double th = std::log(target);
    for (int it = 0; it < 15; ++it) {
        const double e = std::exp(th), f = e - th - target, d = f / (e - 1.0);
        th -= d;
        if (std::abs(d) < 1e-12) break;
    }
    return std::max(0.1, th / (2.0 * k) - 1.166);
}

}  // namespace detalle_obs

template <typename T, std::size_t W = 50, std::size_t S = 2, std::size_t N_REG = 120>
class ObservadorFractalUniversal {
    static_assert(std::is_arithmetic_v<T>, "el observable es un escalar aritmetico (precio, OFI, volumen, profundidad, spread...)");
    static_assert(N_REG >= 8, "la ventana del regularizador necesita al menos 8 pasos");

    ConfiguracionFractal cfg_{};
    CalibracionFractal cal_{};
    SensorDFACausal<W, S> sensor_{};
    RegularizadorDualFractal<N_REG> reg_{};
    std::array<double, N_REG> anillo_h_{};
    std::array<double, N_REG> ventana_h_{};
    std::size_t n_h_ = 0, cab_h_ = 0;
    double sp_ = 0.0, sn_ = 0.0;
    std::int64_t barras_ = 0, pasos_desde_alarma_ = -1;

    // serie del sensor (H) sobre una serie del observable, en pasos validos contiguos (un NaN corta el segmento)
    std::vector<std::vector<double>> sensor_sobre(const std::vector<double>& x) const {
        SensorDFACausal<W, S> s(cfg_.integracion, cfg_.dfa);
        std::vector<std::vector<double>> segs(1);
        for (double v : x)
            if (auto r = s.paso(v)) {
                if (r->valido) segs.back().push_back(r->h);
                else if (!segs.back().empty()) segs.emplace_back();
            }
        if (segs.back().empty()) segs.pop_back();
        return segs;
    }

    // incrementos del observable (lo que el mundo de signos aleatoriza) y reconstruccion de la serie
    std::vector<double> incrementos(const std::vector<double>& x) const {
        std::vector<double> z;
        if (cfg_.integracion == Integracion::NIVEL) {
            for (std::size_t i = 1; i < x.size(); ++i) z.push_back(x[i] - x[i - 1]);
        } else {
            double m = 0.0;
            for (double v : x) m += v;
            m /= static_cast<double>(x.size());
            for (double v : x) z.push_back(v - m);
        }
        return z;
    }
    std::vector<double> reconstruir(const std::vector<double>& z, double ancla) const {
        std::vector<double> x(z.size());
        if (cfg_.integracion == Integracion::NIVEL) {
            double a = ancla;
            for (std::size_t i = 0; i < z.size(); ++i) { a += z[i]; x[i] = a; }
        } else {
            for (std::size_t i = 0; i < z.size(); ++i) x[i] = ancla + z[i];
        }
        return x;
    }
    // mundo de signos: magnitudes |z| de X en su orden (heterocedasticidad real), signo iid (H = 0.5) o de un fGn(H1)
    std::vector<double> mundo(const std::vector<double>& az, std::size_t ini, std::size_t n, double H, std::mt19937_64& rng) const {
        std::vector<double> g = detalle_obs::fgn(n, H, rng), z(n);
        for (std::size_t i = 0; i < n; ++i) z[i] = (g[i] >= 0.0 ? 1.0 : -1.0) * az[(ini + i) % az.size()];
        return z;
    }

    static void momentos(const std::vector<std::vector<double>>& segs, double& mu, double& s0, double& acf1, std::size_t& n) {
        double sm = 0.0; n = 0;
        for (const auto& s : segs) for (double v : s) { sm += v; ++n; }
        mu = n ? sm / static_cast<double>(n) : NaN;
        double sv = 0.0, num = 0.0, den = 0.0;
        for (const auto& s : segs) for (double v : s) sv += (v - mu) * (v - mu);
        s0 = n ? std::sqrt(sv / static_cast<double>(n)) : NaN;
        for (const auto& s : segs) for (std::size_t i = 1; i < s.size(); ++i) num += (s[i - 1] - mu) * (s[i] - mu);
        den = sv;
        acf1 = den > 0 ? num / den : NaN;
    }

public:
    ObservadorFractalUniversal() = default;
    explicit ObservadorFractalUniversal(const ConfiguracionFractal& c) { configurar(c); }

    void configurar(const ConfiguracionFractal& c) {
        cfg_ = c;
        sensor_.configurar(c.integracion, c.dfa);
        cal_ = CalibracionFractal{};
        reiniciar_flujo();
    }
    void reiniciar_flujo() noexcept {
        sensor_.reiniciar();
        n_h_ = cab_h_ = 0;
        sp_ = sn_ = 0.0;
        barras_ = 0;
        pasos_desde_alarma_ = -1;
    }
    [[nodiscard]] const CalibracionFractal& calibracion() const noexcept { return cal_; }
    [[nodiscard]] const ConfiguracionFractal& configuracion() const noexcept { return cfg_; }
    // Para las pruebas: inyecta una calibracion conocida.
    void fijar_calibracion(const CalibracionFractal& c) noexcept { cal_ = c; }

    // Calibra el nulo del observable con su historia causal (todo <= el instante de calibracion). Fuera del ciclo caliente.
    const CalibracionFractal& calibrar(std::span<const T> historia) {
        if (!(cfg_.arl0_objetivo_barras > 0.0)) throw std::invalid_argument("ConfiguracionFractal.arl0_objetivo_barras es obligatorio");
        CalibracionFractal c{};
        c.nulo = cfg_.nulo;
        std::vector<double> x(historia.size());
        for (std::size_t i = 0; i < x.size(); ++i) x[i] = static_cast<double>(historia[i]);
        if (x.size() < W + N_REG * S) return cal_ = c;   // historia insuficiente: lista = false
        std::mt19937_64 rng(cfg_.semilla);

        const std::vector<double> z = incrementos(x);
        std::vector<double> az(z.size());
        for (std::size_t i = 0; i < z.size(); ++i) az[i] = std::abs(z[i]);
        const double ancla = cfg_.integracion == Integracion::NIVEL ? x.front() : (x.empty() ? 0.0 : [&] {
            double m = 0.0; for (double v : x) m += v; return m / static_cast<double>(x.size()); }());

        // mundo de signos de la historia: nulo SIGNOS y referencia H0 de la ganancia y de tau50
        std::vector<std::vector<double>> segs_signos, segs_h1;
        for (std::size_t r = 0; r < cfg_.replicas_nulo_signos; ++r) {
            for (auto& s : sensor_sobre(reconstruir(mundo(az, 0, az.size(), 0.5, rng), ancla))) segs_signos.push_back(std::move(s));
            for (auto& s : sensor_sobre(reconstruir(mundo(az, 0, az.size(), 0.5 + cfg_.delta_h_objetivo, rng), ancla))) segs_h1.push_back(std::move(s));
        }
        double mu_w, s0_w, acf_w, mu_1, s_1, acf_1;
        std::size_t n_w, n_1;
        momentos(segs_signos, mu_w, s0_w, acf_w, n_w);
        momentos(segs_h1, mu_1, s_1, acf_1, n_1);
        c.ganancia_sensor = mu_1 - mu_w;

        // nulo vivo
        std::vector<std::vector<double>> segs = cfg_.nulo == NuloCalibracion::SIGNOS ? segs_signos : sensor_sobre(x);
        momentos(segs, c.mu0, c.s0, c.acf1, c.pasos_calibracion);
        if (!(c.s0 > 0.0) || !(c.ganancia_sensor > 0.0)) return cal_ = c;

        // lambda_TV: cuantil 1 - alfa del maximo del puente sobre ventanas de N_REG pasos del nulo vivo
        std::vector<double> mx_disj, mx_sol;
        for (const auto& s : segs)
            for (std::size_t a = 0; a + N_REG <= s.size(); ++a) {
                const double m = maximo_puente(s.data() + a, N_REG);
                mx_sol.push_back(m);
                if (a % N_REG == 0) mx_disj.push_back(m);
            }
        const bool bastan = static_cast<double>(mx_disj.size()) >= std::ceil(1.0 / cfg_.alfa_lambda);
        c.ventanas_lambda = bastan ? mx_disj.size() : mx_sol.size();
        c.lambda_tv = detalle_obs::cuantil_lineal(bastan ? mx_disj : mx_sol, 1.0 - cfg_.alfa_lambda);
        c.lambda_dj = c.s0 * std::sqrt(2.0 * std::log(static_cast<double>(N_REG)));

        // k = Delta/2 y h por ARL0 sobre la innovacion del nulo vivo
        c.delta = c.ganancia_sensor / c.s0;
        c.k = 0.5 * c.delta;
        c.arl0_objetivo_pasos = cfg_.arl0_objetivo_barras / static_cast<double>(S);
        std::vector<std::vector<double>> nus;
        for (const auto& s : segs) { std::vector<double> v(s.size()); for (std::size_t i = 0; i < s.size(); ++i) v[i] = (s[i] - c.mu0) / c.s0; nus.push_back(std::move(v)); }
        const auto rh = detalle_obs::h_por_arl(nus, c.k, c.arl0_objetivo_pasos, cfg_.alarmas_minimas_arl);
        c.h_crit = rh.h; c.h_extrapolado = rh.extrapolado; c.alarmas_nulo_en_h = rh.alarmas;
        c.h_siegmund = detalle_obs::h_siegmund(c.arl0_objetivo_pasos, c.k);
        if (!std::isfinite(c.h_crit)) return cal_ = c;

        // tau50: saltos 0.5 -> 0.5 + dH inyectados en el mundo de signos de X; detector del mismo mundo (mu_w, s0_w, k_w, h_w)
        const double k_w = 0.5 * c.ganancia_sensor / s0_w;
        std::vector<std::vector<double>> nus_w;
        for (const auto& s : segs_signos) { std::vector<double> v(s.size()); for (std::size_t i = 0; i < s.size(); ++i) v[i] = (s[i] - mu_w) / s0_w; nus_w.push_back(std::move(v)); }
        const double h_w = cfg_.nulo == NuloCalibracion::SIGNOS ? c.h_crit
                                                                 : detalle_obs::h_por_arl(nus_w, k_w, c.arl0_objetivo_pasos, cfg_.alarmas_minimas_arl).h;
        const std::size_t t_pre = W + N_REG * S;                 // reposo previo (barras)
        const std::size_t t_post = cfg_.horizonte_salto_barras ? cfg_.horizonte_salto_barras : 20 * W;
        std::vector<double> taus;
        std::uniform_int_distribution<std::size_t> U(0, az.size() - 1);
        for (std::size_t r = 0; r < cfg_.replicas_salto && std::isfinite(h_w); ++r) {
            const std::size_t ini = U(rng);
            std::vector<double> zz = mundo(az, ini, t_pre, 0.5, rng);
            const std::vector<double> post = mundo(az, ini + t_pre, t_post, 0.5 + cfg_.delta_h_objetivo, rng);
            zz.insert(zz.end(), post.begin(), post.end());
            const std::vector<double> xs = reconstruir(zz, ancla);
            SensorDFACausal<W, S> s(cfg_.integracion, cfg_.dfa);
            double sp = 0.0, sn = 0.0, tau = std::numeric_limits<double>::infinity();
            std::size_t barra = 0;
            for (double v : xs) {
                ++barra;
                auto res = s.paso(v);
                if (!res || !res->valido) continue;
                const double nu = (res->h - mu_w) / s0_w;
                sp = std::max(0.0, sp + nu - k_w);
                sn = std::max(0.0, sn - nu - k_w);
                const bool al_p = sp >= h_w, al_n = sn >= h_w;
                if (al_p || al_n) { sp = sn = 0.0; }
                if (al_p && barra > t_pre) { tau = static_cast<double>(barra - t_pre); break; }   // primera alarma S+ tras el salto
            }
            taus.push_back(tau);
        }
        if (!taus.empty()) {
            std::size_t nd = 0;
            for (double t : taus) nd += !std::isfinite(t);
            c.fraccion_sin_detectar = static_cast<double>(nd) / static_cast<double>(taus.size());
            std::vector<double> ord = taus;
            std::sort(ord.begin(), ord.end());
            const std::size_t m = ord.size();
            c.tau50_barras = (m % 2) ? ord[m / 2] : 0.5 * (ord[m / 2 - 1] + ord[m / 2]);
            c.tau50_pasos = c.tau50_barras / static_cast<double>(S);
        }
        c.lista = true;
        return cal_ = c;
    }

    // Ciclo caliente: una muestra del observable; emite el estado cuando el sensor da paso. Cero reservas.
    std::optional<EstadoFractal> paso(T x_t) noexcept {
        ++barras_;
        auto r = sensor_.paso(static_cast<double>(x_t));
        if (!r) return std::nullopt;
        EstadoFractal e{};
        e.barra = barras_ - 1;
        e.sensor = *r;
        e.tau50_barras = cal_.tau50_barras;
        if (pasos_desde_alarma_ >= 0) ++pasos_desde_alarma_;
        if (!r->valido) {   // ventana degenerada: se corta la ventana del regularizador (sin imputar) y el CUSUM no avanza
            n_h_ = cab_h_ = 0;
            e.norma = cal_.lista ? Norma::CALENTANDO : Norma::SIN_CALIBRAR;
            e.cusum_pos = sp_; e.cusum_neg = sn_;
            e.pasos_desde_alarma = pasos_desde_alarma_;
            return e;
        }
        anillo_h_[(cab_h_ + n_h_) % N_REG] = r->h;
        if (n_h_ < N_REG) ++n_h_; else cab_h_ = (cab_h_ + 1) % N_REG;
        if (!cal_.lista) { e.norma = Norma::SIN_CALIBRAR; return e; }

        e.nu = (r->h - cal_.mu0) / cal_.s0;
        sp_ = std::max(0.0, sp_ + e.nu - cal_.k);
        sn_ = std::max(0.0, sn_ - e.nu - cal_.k);
        e.alarma_pos = sp_ >= cal_.h_crit;
        e.alarma_neg = sn_ >= cal_.h_crit;
        bool sorpresa = sp_ > 0.0 || sn_ > 0.0;    // se lee ANTES del reinicio: la alarma es sorpresa
        e.cusum_pos = sp_; e.cusum_neg = sn_;
        if (e.alarma_pos || e.alarma_neg) { sp_ = sn_ = 0.0; pasos_desde_alarma_ = 0; }
        e.pasos_desde_alarma = pasos_desde_alarma_;

        if (n_h_ < N_REG) { e.norma = Norma::CALENTANDO; return e; }
        for (std::size_t i = 0; i < N_REG; ++i) ventana_h_[i] = anillo_h_[(cab_h_ + i) % N_REG];
        if (cfg_.regla_norma == ReglaNorma::CERTIFICADO_PUENTE) sorpresa = maximo_puente(ventana_h_.data(), N_REG) > cal_.lambda_tv;
        const SalidaRegularizador s = sorpresa ? reg_.l1(ventana_h_.data(), N_REG, cal_.lambda_tv) : reg_.l2(ventana_h_.data(), N_REG, cal_.s0);
        e.norma = sorpresa ? Norma::L1_TV : Norma::L2_TIKHONOV;
        e.h_regularizado = s.borde;
        e.lambda_usado = s.lambda;
        e.n_ventana_regularizador = N_REG;
        return e;
    }
};

// Agregador causal de una escala: cierra la barra de la escala cuando llega la primera observacion del bucket siguiente
// (o cuando el reloj pasa su final, con cerrar_hasta). Nunca emite un bucket abierto; un bucket vacio no se emite (sin imputar).
class AgregadorEscala {
    std::int64_t seg_ = 60;
    std::int64_t bucket_ = std::numeric_limits<std::int64_t>::min();     // bucket abierto (si n_ > 0)
    std::int64_t cerrado_ = std::numeric_limits<std::int64_t>::min();    // ultimo bucket ya emitido
    double acc_ = 0.0;
    std::size_t n_ = 0;
    Agregacion ag_ = Agregacion::ULTIMO;

public:
    AgregadorEscala() = default;
    AgregadorEscala(std::int64_t segundos, Agregacion ag) : seg_(segundos), ag_(ag) {}
    static std::int64_t piso_div(std::int64_t a, std::int64_t b) noexcept { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }
    // t: segundos de epoca UTC de la observacion (su instante de disponibilidad). Devuelve el valor de la barra que se cierra.
    std::optional<double> paso(std::int64_t t, double x) noexcept {
        const std::int64_t b = piso_div(t, seg_);
        if (b <= cerrado_ || (n_ > 0 && b < bucket_)) return std::nullopt;   // atrasada: el flujo debe venir ordenado
        std::optional<double> cerrada;
        if (n_ > 0 && b > bucket_) { cerrada = valor(); cerrado_ = bucket_; n_ = 0; }
        if (n_ == 0) { bucket_ = b; acc_ = 0.0; }
        if (ag_ == Agregacion::ULTIMO) acc_ = x; else acc_ += x;
        ++n_;
        return cerrada;
    }
    std::optional<double> cerrar_hasta(std::int64_t t_reloj) noexcept {
        if (n_ == 0 || (bucket_ + 1) * seg_ > t_reloj) return std::nullopt;
        const double v = valor();
        cerrado_ = bucket_;
        n_ = 0;
        return v;
    }
    [[nodiscard]] double valor() const noexcept { return ag_ == Agregacion::MEDIA ? acc_ / static_cast<double>(n_) : acc_; }
};

// Piramide 1m..D1 de un mismo observable: seis observadores, cada uno con su calibracion y su tau50.
template <typename T, std::size_t W = 50, std::size_t S = 2, std::size_t N_REG = 120>
class PiramideFractal {
    std::array<ObservadorFractalUniversal<double, W, S, N_REG>, NUM_ESCALAS> obs_{};
    std::array<AgregadorEscala, NUM_ESCALAS> ag_{};

public:
    using Salida = std::array<std::optional<EstadoFractal>, NUM_ESCALAS>;
    explicit PiramideFractal(const std::array<ConfiguracionFractal, NUM_ESCALAS>& cfg) {
        for (std::size_t e = 0; e < NUM_ESCALAS; ++e) {
            obs_[e].configurar(cfg[e]);
            ag_[e] = AgregadorEscala(SEGUNDOS_ESCALA[e], cfg[e].agregacion);
        }
    }
    [[nodiscard]] const ObservadorFractalUniversal<double, W, S, N_REG>& escala(Escala e) const { return obs_[static_cast<std::size_t>(e)]; }

    // Calibra cada escala con la historia agregada a su resolucion (solo buckets cerrados de la historia).
    void calibrar(std::span<const std::int64_t> t, std::span<const T> x) {
        for (std::size_t e = 0; e < NUM_ESCALAS; ++e) {
            AgregadorEscala a(SEGUNDOS_ESCALA[e], obs_[e].configuracion().agregacion);
            std::vector<double> serie;
            for (std::size_t i = 0; i < x.size(); ++i)
                if (auto v = a.paso(t[i], static_cast<double>(x[i]))) serie.push_back(*v);
            obs_[e].calibrar(std::span<const double>(serie));
            obs_[e].reiniciar_flujo();
        }
    }
    // Ciclo caliente: una observacion de la barra base. Emite el estado de cada escala cuya barra se cierra ahora.
    Salida paso(std::int64_t t, T x) noexcept {
        Salida s{};
        for (std::size_t e = 0; e < NUM_ESCALAS; ++e)
            if (auto v = ag_[e].paso(t, static_cast<double>(x))) s[e] = obs_[e].paso(*v);
        return s;
    }
};

}  // namespace ark::fractal
