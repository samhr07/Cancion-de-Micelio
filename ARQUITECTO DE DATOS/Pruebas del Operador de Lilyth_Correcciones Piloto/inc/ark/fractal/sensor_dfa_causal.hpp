// sensor_dfa_causal.hpp -- Capa 1 del Bloque Fractal Universal: DFA-1 solapado causal, zero-heap (ARK-FRAC-01).
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-03.
//
// Oraculo: ARQUITECTO DE DATOS/dfa_core (calcular_dfa con ConfiguracionDFA() por omision = DFA_OVERLAPPING_1, s_min 6,
// s_max = max(s_min + 4, W // 4), 12 escalas, residuo cero si hay >= 4 divisores, paso 1), que es el sensor de EL ABISMO y del
// manual de viabilidad (W1 = 50, paso 2). Identidad exigida: |alpha_C++ - alpha_Py| <= 1e-10 (Python hace la proyeccion con
// inv(X'X) y np.mean por pares; aqui el ajuste lineal es cerrado, asi que el redondeo difiere en el ultimo bit).
//
// DIVERGE DEL ORACULO (declarado, configurable):
//   la ventana "constante" en Python es np.isclose(min, max, atol=1e-12), que lleva rtol = 1e-5 relativo a max: una ventana
//   de log-precios (~11) con recorrido < 1.1e-4 o de precios en 60 000 con recorrido < 0.6 se declara degenerada. Aqui, por
//   omision, solo se exige recorrido > 1e-12; `degenerado_relativo_python = true` reproduce el oraculo (lo usa la prueba C1).
//
// Causalidad: cada salida usa solo las ultimas W muestras recibidas (la ventana termina en el borde presente). El primer
// paso sale con W muestras y despues cada S, como MotorDFATTF.flujo_causal_streaming.
//
// DECLARACION METROLOGICA (EL ABISMO - COTA I & II):
// Observador de estado metrologico del medio. Queda terminantemente prohibido
// usar las salidas de H en escalas < 1h como gatillos reactivos directos taker.
#pragma once

#include "ark/fractal/tipos_fractal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>

namespace ark::fractal {

struct ConfiguracionDFA {
    int s_min = 6;
    int s_max = 0;                       // 0 = max(s_min + 4, W / 4), como el oraculo
    int num_escalas = 12;
    bool forzar_residuo_cero = true;
    int min_muestras = 32;
    bool degenerado_relativo_python = false;
};

inline constexpr std::size_t MAX_ESCALAS_DFA = 64;

// generar_lista_escalas de estimadores_dfa.py. Devuelve el numero de escalas (0 si el rango es invalido).
inline std::size_t generar_escalas(int W, int s_min, int s_max, int num, bool residuo_cero, std::array<int, MAX_ESCALAS_DFA>& out) {
    if (s_min < 4 || s_max > W / 2 || s_min >= s_max || num < 2 || num > static_cast<int>(MAX_ESCALAS_DFA)) return 0;
    std::size_t n = 0;
    if (residuo_cero) {
        std::array<int, MAX_ESCALAS_DFA * 4> divs{};
        std::size_t nd = 0;
        for (int s = s_min; s <= s_max && nd < divs.size(); ++s)
            if (W % s == 0) divs[nd++] = s;
        if (nd >= 4) {
            if (nd <= static_cast<std::size_t>(num)) {
                for (std::size_t i = 0; i < nd && i < MAX_ESCALAS_DFA; ++i) out[n++] = divs[i];
                return n;
            }
            // muestreo cuasi-logaritmico sobre los divisores: argmin |log d - t|, t = linspace(log d0, log d_last, num)
            const double a = std::log(static_cast<double>(divs[0])), b = std::log(static_cast<double>(divs[nd - 1]));
            const double paso = (b - a) / static_cast<double>(num - 1);
            std::array<int, MAX_ESCALAS_DFA> sel{};
            for (int i = 0; i < num; ++i) {
                const double t = (i == num - 1) ? b : a + paso * static_cast<double>(i);
                std::size_t mejor = 0;
                double dm = std::abs(std::log(static_cast<double>(divs[0])) - t);
                for (std::size_t j = 1; j < nd; ++j) {
                    const double d = std::abs(std::log(static_cast<double>(divs[j])) - t);
                    if (d < dm) { dm = d; mejor = j; }
                }
                sel[static_cast<std::size_t>(i)] = divs[mejor];
            }
            std::sort(sel.begin(), sel.begin() + num);
            n = static_cast<std::size_t>(std::unique(sel.begin(), sel.begin() + num) - sel.begin());
            if (n >= 4) { std::copy(sel.begin(), sel.begin() + static_cast<std::ptrdiff_t>(n), out.begin()); return n; }
            // colapso por duplicidad: indices np.round(linspace(0, nd - 1, num))
            const double pi = static_cast<double>(nd - 1) / static_cast<double>(num - 1);
            for (int i = 0; i < num; ++i) {
                const double v = (i == num - 1) ? static_cast<double>(nd - 1) : pi * static_cast<double>(i);
                sel[static_cast<std::size_t>(i)] = divs[static_cast<std::size_t>(std::nearbyint(v))];
            }
            std::sort(sel.begin(), sel.begin() + num);
            n = static_cast<std::size_t>(std::unique(sel.begin(), sel.begin() + num) - sel.begin());
            std::copy(sel.begin(), sel.begin() + static_cast<std::ptrdiff_t>(n), out.begin());
            return n;
        }
    }
    // np.unique(np.round(np.logspace(log10 s_min, log10 s_max, num)).astype(int))
    const double a = std::log10(static_cast<double>(s_min)), b = std::log10(static_cast<double>(s_max));
    const double paso = (b - a) / static_cast<double>(num - 1);
    std::array<int, MAX_ESCALAS_DFA> sel{};
    for (int i = 0; i < num; ++i) {
        const double e = (i == num - 1) ? b : a + paso * static_cast<double>(i);
        sel[static_cast<std::size_t>(i)] = static_cast<int>(std::nearbyint(std::pow(10.0, e)));
    }
    std::sort(sel.begin(), sel.begin() + num);
    n = static_cast<std::size_t>(std::unique(sel.begin(), sel.begin() + num) - sel.begin());
    if (n < 4) return 0;
    std::copy(sel.begin(), sel.begin() + static_cast<std::ptrdiff_t>(n), out.begin());
    return n;
}

// DFA-1 solapado sobre una ventana x[0..W-1] (calcular_dfa_vectorizado con solapamiento, paso 1, orden 1).
// `y` es memoria de trabajo de al menos W dobles. No reserva memoria.
inline ResultadoSensor dfa1_solapado(const double* x, std::size_t W, const std::array<int, MAX_ESCALAS_DFA>& esc, std::size_t n_esc,
                                     Integracion d, const ConfiguracionDFA& cfg, double* y) noexcept {
    ResultadoSensor r{};
    if (n_esc < 2 || W < static_cast<std::size_t>(cfg.min_muestras)) return r;
    double mn = x[0], mx = x[0], suma = 0.0;
    for (std::size_t i = 0; i < W; ++i) {
        if (!std::isfinite(x[i])) return r;
        mn = std::min(mn, x[i]); mx = std::max(mx, x[i]);
        suma += x[i];
    }
    const double tol = 1e-12 + (cfg.degenerado_relativo_python ? 1e-5 * std::abs(mx) : 0.0);
    if (std::abs(mn - mx) <= tol) return r;
    const double media = suma / static_cast<double>(W);
    double acc = 0.0;
    for (std::size_t i = 0; i < W; ++i) { acc += x[i] - media; y[i] = acc; }   // perfil integrado

    std::array<double, MAX_ESCALAS_DFA> ls{}, lf{};
    for (std::size_t e = 0; e < n_esc; ++e) {
        const std::size_t s = static_cast<std::size_t>(esc[e]);
        const double tb = 0.5 * static_cast<double>(s - 1);
        const double stt = static_cast<double>(s) * (static_cast<double>(s) * static_cast<double>(s) - 1.0) / 12.0;
        double ss = 0.0;
        for (std::size_t i = 0; i + s <= W; ++i) {           // cajas deslizantes paso 1
            double ym = 0.0;
            for (std::size_t j = 0; j < s; ++j) ym += y[i + j];
            ym /= static_cast<double>(s);
            double sty = 0.0;
            for (std::size_t j = 0; j < s; ++j) sty += (static_cast<double>(j) - tb) * (y[i + j] - ym);
            const double b = sty / stt;
            for (std::size_t j = 0; j < s; ++j) {
                const double res = y[i + j] - ym - b * (static_cast<double>(j) - tb);
                ss += res * res;
            }
        }
        const double f2 = ss / static_cast<double>((W - s + 1) * s);
        ls[e] = std::log(static_cast<double>(s));
        lf[e] = std::log(std::sqrt(std::max(1e-16, f2)));
    }
    double msx = 0.0, mfy = 0.0;
    for (std::size_t e = 0; e < n_esc; ++e) { msx += ls[e]; mfy += lf[e]; }
    msx /= static_cast<double>(n_esc); mfy /= static_cast<double>(n_esc);
    double sxx = 0.0, sxy = 0.0, syy = 0.0;
    for (std::size_t e = 0; e < n_esc; ++e) {
        sxx += (ls[e] - msx) * (ls[e] - msx);
        sxy += (ls[e] - msx) * (lf[e] - mfy);
        syy += (lf[e] - mfy) * (lf[e] - mfy);
    }
    if (sxx <= 0.0) return r;
    const double alpha = sxy / sxx;
    const double c = mfy - alpha * msx;
    double ssr = 0.0;
    for (std::size_t e = 0; e < n_esc; ++e) { const double q = lf[e] - (c + alpha * ls[e]); ssr += q * q; }
    const double df = std::max(1.0, static_cast<double>(n_esc) - 2.0);
    r.alpha = alpha;
    r.h = alpha - static_cast<double>(static_cast<int>(d));
    r.r2 = syy > 0.0 ? 1.0 - ssr / syy : 0.0;
    r.se = std::sqrt((ssr / df) / sxx);
    r.valido = true;
    return r;
}

// Sensor rodante: anillo de W muestras en la pila; evalua con W muestras y despues cada S.
template <std::size_t W = 50, std::size_t S = 2>
class SensorDFACausal {
    static_assert(W >= 20, "MotorDFATTF exige lookback_W1 >= 20");
    static_assert(S >= 1, "stride_S >= 1");

    std::array<double, W> anillo_{};
    std::array<double, W> ventana_{};
    std::array<double, W> perfil_{};
    std::array<int, MAX_ESCALAS_DFA> escalas_{};
    std::size_t n_escalas_ = 0;
    std::size_t recibidas_ = 0;
    Integracion d_ = Integracion::NIVEL;
    ConfiguracionDFA cfg_{};

public:
    SensorDFACausal() { configurar(Integracion::NIVEL, ConfiguracionDFA{}); }
    explicit SensorDFACausal(Integracion d, ConfiguracionDFA cfg = {}) { configurar(d, cfg); }

    void configurar(Integracion d, ConfiguracionDFA cfg) {
        d_ = d;
        cfg_ = cfg;
        cfg_.min_muestras = std::min(cfg_.min_muestras, static_cast<int>(W));   // MotorDFATTF: min(min_muestras, W1)
        const int smax = cfg_.s_max > 0 ? cfg_.s_max : std::max(cfg_.s_min + 4, static_cast<int>(W) / 4);
        n_escalas_ = generar_escalas(static_cast<int>(W), cfg_.s_min, smax, cfg_.num_escalas, cfg_.forzar_residuo_cero, escalas_);
        reiniciar();
    }
    void reiniciar() noexcept { recibidas_ = 0; }

    [[nodiscard]] std::size_t n_escalas() const noexcept { return n_escalas_; }
    [[nodiscard]] const std::array<int, MAX_ESCALAS_DFA>& escalas() const noexcept { return escalas_; }
    [[nodiscard]] std::size_t recibidas() const noexcept { return recibidas_; }

    // Empuja una muestra; devuelve el resultado si en esta muestra toca evaluar.
    std::optional<ResultadoSensor> paso(double x) noexcept {
        anillo_[recibidas_ % W] = x;
        ++recibidas_;
        if (recibidas_ < W || (recibidas_ - W) % S != 0) return std::nullopt;
        const std::size_t ini = recibidas_ % W;            // la mas antigua de las ultimas W
        for (std::size_t i = 0; i < W; ++i) ventana_[i] = anillo_[(ini + i) % W];
        return dfa1_solapado(ventana_.data(), W, escalas_, n_escalas_, d_, cfg_, perfil_.data());
    }
};

}  // namespace ark::fractal
