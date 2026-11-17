// whittle_exacto.hpp -- sensor de memoria por Whittle local EXACTO con media desconocida (ARK-METRO-01), C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04, por orden del Comandante. Prerregistro:
// JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/metrologia/EXPECTATIVAS_WHITTLE.md (escrito antes de este codigo).
//
// Fuente (SECUNDARIA, verificada): Blevins (2025), arXiv:2511.15689v4, secs. 2.1, 2.4, 2.5, que reexpone a Shimotsu y Phillips
// (2005, Ann. Statist. 33) y Shimotsu (2010, Econometric Theory 26(2)):
//   Delta^d X_t = sum_{k=0}^{t-1} pi_k(d) X_{t-k},  pi_0 = 1,  pi_k = pi_{k-1} (k - 1 - d) / k          (tipo II)
//   R(d) = log( (1/m) sum_{j=1}^m I_{Delta^d (X - mu~(d))}(lambda_j) ) - (2d/m) sum_{j=1}^m log lambda_j,  lambda_j = 2 pi j / n
//   mu~(d) = w(d) Xbar + (1 - w(d)) X_1;  w = 1 (d <= 1/2), (1/2)[1 + cos(4 pi d)] (1/2 < d < 3/4), 0 (d >= 3/4)
//   valido para d0 en (-1/2, 2); ancho de busqueda <= 9/2; sqrt(m)(d^ - d0) -> N(0, 1/4).
// DIVERGE DE SHIMOTSU (2010) (declarado en el prerregistro): minimizacion GLOBAL del objetivo modificado (rejilla de paso se(m) y
// refinamiento aureo en la celda del minimo) en vez de los dos pasos con taper, cuyo taper la fuente no especifica.
//
// Ciclo caliente (SensorWhittleExacto::paso): cero reservas. La derivacion de la banda (calibracion) reserva.
#pragma once

#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace ark::metrologia {

// ------------------------------------------------------------------------------------------------------------------------
namespace elw {

inline constexpr double PI = 3.14159265358979323846;
inline constexpr double ANCHO_MAX = 4.5;   // Shimotsu-Phillips: ancho del intervalo de busqueda <= 9/2

// Peso de Shimotsu (2010) para la media desconocida.
inline double peso_media(double d) noexcept {
    if (d <= 0.5) return 1.0;
    if (d >= 0.75) return 0.0;
    return 0.5 * (1.0 + std::cos(4.0 * PI * d));
}

// pi_k(d), k = 0..n-1 (coeficientes de (1-L)^d). `pi` con al menos n entradas.
inline void coeficientes(double d, std::size_t n, double* pi) noexcept {
    if (n == 0) return;
    pi[0] = 1.0;
    for (std::size_t k = 1; k < n; ++k) pi[k] = pi[k - 1] * (static_cast<double>(k) - 1.0 - d) / static_cast<double>(k);
}

// FFT radix 2 in situ (n potencia de 2). sentido = -1 directa (e^{-i}), +1 inversa sin normalizar.
inline void fft(std::complex<double>* a, std::size_t n, int sentido) noexcept {
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double ang = 2.0 * PI / static_cast<double>(len) * static_cast<double>(sentido);
        const std::complex<double> wl(std::cos(ang), std::sin(ang));
        for (std::size_t i = 0; i < n; i += len) {
            std::complex<double> w(1.0, 0.0);
            for (std::size_t j = 0; j < len / 2; ++j) {
                if ((j & 63) == 0 && j) w = std::complex<double>(std::cos(ang * static_cast<double>(j)), std::sin(ang * static_cast<double>(j)));
                const std::complex<double> u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wl;
            }
        }
    }
}

inline std::size_t potencia2_mayor_igual(std::size_t n) noexcept { std::size_t p = 1; while (p < n) p <<= 1; return p; }

// Media de log lambda_j, j = 1..m, con n puntos.
inline double media_log_lambda(std::size_t m, std::size_t n) noexcept {
    double s = 0.0;
    for (std::size_t j = 1; j <= m; ++j) s += std::log(2.0 * PI * static_cast<double>(j) / static_cast<double>(n));
    return s / static_cast<double>(m);
}

// Error estandar asintotico de d^ con m frecuencias.
inline double se_asintotico(std::size_t m) noexcept { return m ? 1.0 / (2.0 * std::sqrt(static_cast<double>(m))) : NaN; }

// Objetivo R(d) por el camino del FLUJO: convolucion directa y DFT con la tabla de raices (raiz[k] = e^{i 2 pi k / n}).
// x: ventana (n); trabajo: pi (n), u (n). No reserva.
inline double objetivo_directo(const double* x, std::size_t n, std::size_t m, double d, double xbar, double mlog,
                               const std::complex<double>* raiz, double* pi, double* u) noexcept {
    coeficientes(d, n, pi);
    const double w = peso_media(d), mu = w * xbar + (1.0 - w) * x[0];
    for (std::size_t i = 0; i < n; ++i) {               // u_t = sum_{k=0}^{t-1} pi_k (x_{t-k} - mu), t = i + 1
        double s = 0.0;
        for (std::size_t k = 0; k <= i; ++k) s += pi[k] * (x[i - k] - mu);
        u[i] = s;
    }
    double suma = 0.0;
    for (std::size_t j = 1; j <= m; ++j) {
        double re = 0.0, im = 0.0;
        std::size_t idx = j % n;                          // t = 1: indice j * 1
        for (std::size_t i = 0; i < n; ++i) {
            re += u[i] * raiz[idx].real();
            im += u[i] * raiz[idx].imag();
            idx += j; if (idx >= n) idx -= n;
        }
        suma += re * re + im * im;
    }
    const double I = suma / (2.0 * PI * static_cast<double>(n) * static_cast<double>(m));
    return std::log(std::max(I, 1e-300)) - 2.0 * d * mlog;
}

struct Minimo { double d = NaN; double R = NaN; bool en_borde = false; bool valido = false; };

// Minimizacion global sobre [d_min, d_max]: rejilla de paso se(m) y refinamiento aureo en la celda del mejor punto.
template <typename F>
inline Minimo minimizar(F&& R, double d_min, double d_max, std::size_t m) noexcept {
    Minimo r{};
    const double h = se_asintotico(m);
    if (!(h > 0.0) || !(d_max > d_min)) return r;
    const std::size_t G = static_cast<std::size_t>(std::ceil((d_max - d_min) / h)) + 1;
    const double paso = (d_max - d_min) / static_cast<double>(G - 1);
    std::size_t b = 0;
    double Rb = std::numeric_limits<double>::infinity();
    for (std::size_t g = 0; g < G; ++g) {
        const double v = R(d_min + paso * static_cast<double>(g));
        if (std::isfinite(v) && v < Rb) { Rb = v; b = g; }
    }
    if (!std::isfinite(Rb)) return r;
    double a = d_min + paso * static_cast<double>(b == 0 ? 0 : b - 1);
    double c = d_min + paso * static_cast<double>(b + 1 >= G ? G - 1 : b + 1);
    const double gr = 0.5 * (std::sqrt(5.0) - 1.0);
    double x1 = c - gr * (c - a), x2 = a + gr * (c - a);
    double f1 = R(x1), f2 = R(x2);
    const double tol = 1e-7;
    while (c - a > tol) {
        if (f1 <= f2) { c = x2; x2 = x1; f2 = f1; x1 = c - gr * (c - a); f1 = R(x1); }
        else { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (c - a); f2 = R(x2); }
    }
    r.d = 0.5 * (a + c);
    r.R = R(r.d);
    if (Rb < r.R) { r.d = d_min + paso * static_cast<double>(b); r.R = Rb; }   // el aureo nunca empeora la rejilla
    r.en_borde = (r.d - d_min) < paso * 0.5 || (d_max - r.d) < paso * 0.5;
    r.valido = std::isfinite(r.R) && !r.en_borde;
    return r;
}

// Objetivo para TODOS los m = 1..M a la vez (camino de la calibracion, FFT): devuelve en S[j-1] = sum_{l<=j} I_l(d) (sin la
// constante 2 pi n) para j = 1..M. x0: FFT (longitud P >= 2n) de la ventana relleno con ceros, ya calculada. Reserva en `trabajo`.
inline void sumas_periodograma(const double* x, std::size_t n, const std::vector<std::complex<double>>& X_fft, std::size_t P,
                               double d, std::size_t M, std::vector<double>& S, std::vector<std::complex<double>>& trabajo,
                               std::vector<std::complex<double>>& trabajo_n, std::vector<double>& pi, double xbar) {
    pi.resize(n);
    coeficientes(d, n, pi.data());
    trabajo.assign(P, {0.0, 0.0});
    for (std::size_t k = 0; k < n; ++k) trabajo[k] = pi[k];
    fft(trabajo.data(), P, -1);
    for (std::size_t k = 0; k < P; ++k) trabajo[k] *= X_fft[k];
    fft(trabajo.data(), P, +1);                          // conv(pi, x)[i] * P
    const double w = peso_media(d), mu = w * xbar + (1.0 - w) * x[0];
    trabajo_n.assign(n, {0.0, 0.0});
    double cum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        cum += pi[i];                                    // sum_{k<=i} pi_k
        trabajo_n[i] = trabajo[i].real() / static_cast<double>(P) - mu * cum;
    }
    fft(trabajo_n.data(), n, -1);                        // |.| igual al convenio t = 1..n con e^{+i}
    S.resize(M);
    double s = 0.0;
    for (std::size_t j = 1; j <= M; ++j) { s += std::norm(trabajo_n[j]); S[j - 1] = s; }
}

}  // namespace elw

// ------------------------------------------------------------------------------------------------------------------------
struct ResultadoELW {
    double d = NaN;        // d^ de la serie tal cual
    double h = NaN;        // d^ + 1/2 - orden
    double se = NaN;       // 1 / (2 sqrt m)
    double R = NaN;        // objetivo en el minimo
    bool en_borde = false; // minimo pegado a d_min o d_max: no es una medicion
    bool valido = false;
};

struct ConfiguracionSensorELW {
    std::size_t W = 0;           // ventana (<= W_MAX)
    std::size_t m = 0;           // frecuencias (1 <= m <= (W - 1) / 2)
    std::size_t S = 1;           // paso de salida
    Integracion orden = Integracion::FLUJO;
    double d_min = -0.5, d_max = 2.0;
};

// Sensor rodante causal: evalua con W muestras y despues cada S; la ventana termina en el borde presente.
template <std::size_t W_MAX = 1024>
class SensorWhittleExacto {
    ConfiguracionSensorELW cfg_{};
    std::array<double, W_MAX> anillo_{}, ventana_{}, pi_{}, u_{};
    std::array<std::complex<double>, W_MAX> raiz_{};
    double mlog_ = NaN;
    std::size_t recibidas_ = 0;

public:
    SensorWhittleExacto() = default;
    explicit SensorWhittleExacto(const ConfiguracionSensorELW& c) { configurar(c); }

    void configurar(const ConfiguracionSensorELW& c) {
        if (c.W < 4 || c.W > W_MAX) throw std::invalid_argument("SensorWhittleExacto: W fuera de [4, W_MAX]");
        if (c.m < 1 || c.m > (c.W - 1) / 2) throw std::invalid_argument("SensorWhittleExacto: m fuera de [1, (W-1)/2]");
        if (c.S < 1) throw std::invalid_argument("SensorWhittleExacto: S >= 1");
        if (!(c.d_max > c.d_min) || c.d_max - c.d_min > elw::ANCHO_MAX) throw std::invalid_argument("SensorWhittleExacto: ancho de d > 9/2");
        cfg_ = c;
        for (std::size_t k = 0; k < c.W; ++k) {
            const double a = 2.0 * elw::PI * static_cast<double>(k) / static_cast<double>(c.W);
            raiz_[k] = {std::cos(a), std::sin(a)};
        }
        mlog_ = elw::media_log_lambda(c.m, c.W);
        reiniciar();
    }
    void reiniciar() noexcept { recibidas_ = 0; }
    [[nodiscard]] const ConfiguracionSensorELW& configuracion() const noexcept { return cfg_; }

    // Estima sobre una ventana dada (n = W). No reserva.
    ResultadoELW estimar(const double* x) noexcept {
        ResultadoELW r{};
        const std::size_t n = cfg_.W;
        double mn = x[0], mx = x[0], s = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(x[i])) return r;
            mn = std::min(mn, x[i]); mx = std::max(mx, x[i]); s += x[i];
        }
        if (!(mx - mn > 1e-12 * (1.0 + std::abs(mx)))) return r;   // ventana constante: no hay periodograma
        const double xbar = s / static_cast<double>(n);
        auto R = [&](double d) { return elw::objetivo_directo(x, n, cfg_.m, d, xbar, mlog_, raiz_.data(), pi_.data(), u_.data()); };
        const auto mn_ = elw::minimizar(R, cfg_.d_min, cfg_.d_max, cfg_.m);
        r.d = mn_.d;
        r.R = mn_.R;
        r.en_borde = mn_.en_borde;
        r.valido = mn_.valido;
        r.se = elw::se_asintotico(cfg_.m);
        r.h = r.d + 0.5 - static_cast<double>(static_cast<int>(cfg_.orden));
        return r;
    }

    // Empuja una muestra (NaN = invalida); devuelve resultado si toca evaluar.
    std::optional<ResultadoELW> paso(double x) noexcept {
        anillo_[recibidas_ % cfg_.W] = x;
        ++recibidas_;
        if (recibidas_ < cfg_.W || (recibidas_ - cfg_.W) % cfg_.S != 0) return std::nullopt;
        const std::size_t ini = recibidas_ % cfg_.W;
        for (std::size_t i = 0; i < cfg_.W; ++i) ventana_[i] = anillo_[(ini + i) % cfg_.W];
        return estimar(ventana_.data());
    }
};

// ------------------------------------------------------------------------------------------------------------------------
// Derivacion de la banda (calibracion; reserva). Todo derivado de: ruido_objetivo (s*), ventana_interna (w), W_MAX, alfa.
struct ConfiguracionBandaELW {
    double ruido_objetivo = NaN;     // OBLIGATORIO: sd de H^ por ventana exigida (suelo de ruido)
    std::size_t ventana_interna = 1; // w de la estadistica rodante que construye el observable (1: no es rodante)
    double alfa = ALFA_PROYECTO;
    double d_min = -0.5, d_max = 2.0;
};

struct FichaBandaELW {
    bool viable = false;
    const char* motivo = "sin calibrar";
    double ruido_objetivo = NaN;
    std::size_t m_suelo = 0;            // ceil(1 / (4 s*^2))
    std::size_t W_ref = 0;              // ventana de referencia del perfil (mayor potencia de 2 <= W_MAX)
    std::size_t ventanas_perfil = 0;    // ventanas disjuntas de W_ref
    std::size_t m_estrella = 0;         // techo de localidad en W_ref
    double lambda_estrella = NaN;       // 2 pi m* / W_ref
    double lambda_interna = NaN;        // 2 pi / w (techo de la ventana interna)
    std::size_t W = 0, m = 0;           // elegidos
    double se = NaN;                    // 1 / (2 sqrt m)
    double frecuencia_corte = NaN;      // 2 pi m / W de la banda elegida
    // perfil de la deriva pareada (para la ficha): m, D(m), se(D), cota superior
    std::vector<std::size_t> perfil_m;
    std::vector<double> perfil_D, perfil_se, perfil_cota;
    double d_ref_medio = NaN;           // media de d^(m_suelo; W_ref)
    // enmienda 1: sesgo medido EN cada W candidata (sub-ventanas pareadas contra d^(m_suelo; W_ref) de su ventana)
    std::vector<std::size_t> perfil_W;
    std::vector<double> perfil_W_D, perfil_W_se, perfil_W_cota;
    double sesgo_W = NaN, se_sesgo_W = NaN;   // de la W elegida (ya incluye el de la referencia)
    // enmienda 2: sesgo de la propia referencia por la ley (m/n)^2 del sesgo de segundo orden
    double sesgo_referencia = NaN, se_sesgo_referencia = NaN;
    // enmienda 3 (decision del Comandante): el suelo de ruido se DERIVA DEL NULO. m se itera hasta que la sd medida de d^(m; W)
    // sobre ventanas disjuntas del nulo de ruido quede <= s*; m_suelo queda como el punto de partida asintotico.
    std::size_t m_asintotico = 0;          // ceil(1 / (4 s*^2))
    double ruido_nulo = NaN;               // sd medida de d^ en la (W, m) elegida sobre el nulo de ruido
    std::size_t ventanas_ruido = 0;
    std::vector<std::size_t> historial_m, historial_W;
    std::vector<double> historial_ruido;
};

namespace elw {

// d^(m) de una ventana con el MISMO camino que el sensor (convolucion directa, rejilla + aureo). Calibracion: reserva.
inline Minimo estimar_calibracion(const double* x, std::size_t n, std::size_t m, double d_min, double d_max) {
    Minimo r{};
    double mn = x[0], mx = x[0], s = 0.0;
    for (std::size_t i = 0; i < n; ++i) { if (!std::isfinite(x[i])) return r; mn = std::min(mn, x[i]); mx = std::max(mx, x[i]); s += x[i]; }
    if (!(mx - mn > 1e-12 * (1.0 + std::abs(mx)))) return r;
    std::vector<std::complex<double>> raiz(n);
    for (std::size_t k = 0; k < n; ++k) {   // igual que SensorWhittleExacto::configurar
        const double a = 2.0 * PI * static_cast<double>(k) / static_cast<double>(n);
        raiz[k] = {std::cos(a), std::sin(a)};
    }
    std::vector<double> pi(n), u(n);
    const double xbar = s / static_cast<double>(n), mlog = media_log_lambda(m, n);
    return minimizar([&](double d) { return objetivo_directo(x, n, m, d, xbar, mlog, raiz.data(), pi.data(), u.data()); }, d_min, d_max, m);
}

// d^(m) para todos los m = m_lo..M de una ventana (FFT): rejilla fina de paso se(M)/2 y parabola en los tres puntos del minimo.
inline void perfil_ventana(const double* x, std::size_t n, std::size_t m_lo, std::size_t M, double d_min, double d_max,
                           std::vector<double>& d_out) {
    const std::size_t P = potencia2_mayor_igual(2 * n);
    std::vector<std::complex<double>> X(P, {0.0, 0.0}), tr, trn;
    double xbar = 0.0;
    for (std::size_t i = 0; i < n; ++i) { X[i] = x[i]; xbar += x[i]; }
    xbar /= static_cast<double>(n);
    fft(X.data(), P, -1);
    const double h = 0.5 * se_asintotico(M);
    const std::size_t G = static_cast<std::size_t>(std::ceil((d_max - d_min) / h)) + 1;
    const double paso = (d_max - d_min) / static_cast<double>(G - 1);
    std::vector<double> mlog(M + 1, 0.0);
    {
        double s = 0.0;
        for (std::size_t j = 1; j <= M; ++j) { s += std::log(2.0 * PI * static_cast<double>(j) / static_cast<double>(n)); mlog[j] = s / static_cast<double>(j); }
    }
    std::vector<double> Rg(G * (M + 1), NaN), S, pi;
    for (std::size_t g = 0; g < G; ++g) {
        const double d = d_min + paso * static_cast<double>(g);
        sumas_periodograma(x, n, X, P, d, M, S, tr, trn, pi, xbar);
        for (std::size_t m = m_lo; m <= M; ++m) Rg[g * (M + 1) + m] = std::log(std::max(S[m - 1] / static_cast<double>(m), 1e-300)) - 2.0 * d * mlog[m];
    }
    d_out.assign(M + 1, NaN);
    for (std::size_t m = m_lo; m <= M; ++m) {
        std::size_t b = 0;
        double Rb = std::numeric_limits<double>::infinity();
        for (std::size_t g = 0; g < G; ++g) { const double v = Rg[g * (M + 1) + m]; if (std::isfinite(v) && v < Rb) { Rb = v; b = g; } }
        if (!std::isfinite(Rb) || b == 0 || b + 1 >= G) continue;   // en el borde: no es medicion
        const double y0 = Rg[(b - 1) * (M + 1) + m], y1 = Rb, y2 = Rg[(b + 1) * (M + 1) + m];
        const double den = y0 - 2.0 * y1 + y2;
        const double off = den > 0.0 ? 0.5 * (y0 - y2) / den : 0.0;
        d_out[m] = d_min + paso * (static_cast<double>(b) + std::clamp(off, -1.0, 1.0));
    }
}

}  // namespace elw

// segs: tramos contiguos validos del observable tal cual (calibracion). W_MAX: tope de la ventana.
// segs_nulo: los mismos tramos en el NULO DE RUIDO (enmienda 3, obligatorio; lo construye el llamador con el nulo que elija el
// usuario para su variable). Sin el (nullptr) el suelo queda asintotico y la ficha lo dice.
inline FichaBandaELW derivar_banda(const std::vector<std::vector<double>>& segs, std::size_t W_MAX, const ConfiguracionBandaELW& c,
                                   const std::vector<std::vector<double>>* segs_nulo = nullptr) {
    FichaBandaELW f{};
    f.ruido_objetivo = c.ruido_objetivo;
    if (!(c.ruido_objetivo > 0.0)) throw std::invalid_argument("ConfiguracionBandaELW.ruido_objetivo es obligatorio");
    if (c.ventana_interna < 1) throw std::invalid_argument("ventana_interna >= 1");
    const std::size_t m0 = static_cast<std::size_t>(std::ceil(1.0 / (4.0 * c.ruido_objetivo * c.ruido_objetivo)));
    f.m_asintotico = m0;
    f.m_suelo = m0;
    f.lambda_interna = 2.0 * elw::PI / static_cast<double>(c.ventana_interna);
    std::size_t Wr = 1;
    while (Wr * 2 <= W_MAX) Wr *= 2;
    f.W_ref = Wr;
    const std::size_t M = (Wr - 1) / 2;
    if (M <= m0) { f.motivo = "el suelo de ruido no cabe ni en W_MAX (m_suelo >= Nyquist)"; return f; }
    // ventanas disjuntas de W_ref, en orden, agrupadas por segmento (para el tau_int de las diferencias); el perfil cubre m >= m0
    struct VentanaRef { const double* x; std::size_t grupo; double d_ref; };
    std::vector<VentanaRef> refs;
    std::vector<std::vector<std::vector<double>>> d_por_seg;   // [grupo][ventana] -> d^(m), m = 0..M
    for (const auto& g : segs) {
        std::vector<std::vector<double>> dv;
        for (std::size_t a = 0; a + Wr <= g.size(); a += Wr) {
            const double* x = g.data() + a;
            double mn = x[0], mx = x[0];
            for (std::size_t i = 0; i < Wr; ++i) { mn = std::min(mn, x[i]); mx = std::max(mx, x[i]); }
            if (!(mx - mn > 1e-12 * (1.0 + std::abs(mx)))) continue;
            std::vector<double> d;
            elw::perfil_ventana(x, Wr, m0, M, c.d_min, c.d_max, d);
            dv.push_back(std::move(d));
            refs.push_back({x, d_por_seg.size(), NaN});
        }
        if (!dv.empty()) d_por_seg.push_back(std::move(dv));
    }
    f.ventanas_perfil = refs.size();
    const std::size_t nmin = static_cast<std::size_t>(std::ceil(1.0 / c.alfa));
    if (refs.size() < nmin) { f.motivo = "historia insuficiente (< 1/alfa ventanas disjuntas de W_ref)"; return f; }
    const double z = est::norm_ppf(1.0 - c.alfa);
    auto resumen = [&](std::vector<std::vector<double>>& diffs, double& D, double& se, double& cota) -> bool {
        diffs.erase(std::remove_if(diffs.begin(), diffs.end(), [](const std::vector<double>& v) { return v.empty(); }), diffs.end());
        const auto mo = est::momentos(diffs);
        if (mo.n < nmin || !(mo.sd >= 0.0)) return false;
        double tau = 1.0;
        if (mo.sd > 0.0) {
            const auto ti = est::tau_integrado(diffs, mo.mu, mo.sd, c.alfa, std::max<std::size_t>(1, mo.n / 4));
            if (std::isfinite(ti.tau)) tau = std::max(1.0, ti.tau);
        }
        D = mo.mu;
        se = mo.sd * std::sqrt(tau / static_cast<double>(mo.n));
        cota = std::abs(D) + z * se;
        return true;
    };

    // Banda para un m dado (enmiendas 1 y 2): referencia d^(m; W_ref), perfil de localidad relativo a m, sesgo de la referencia y
    // busqueda de la menor W con sesgo total acotado por s*. Devuelve false (con motivo) si no hay W.
    auto banda_para = [&](std::size_t m, FichaBandaELW& g) -> bool {
        g.m_suelo = m;
        g.perfil_m.clear(); g.perfil_D.clear(); g.perfil_se.clear(); g.perfil_cota.clear();
        g.perfil_W.clear(); g.perfil_W_D.clear(); g.perfil_W_se.clear(); g.perfil_W_cota.clear();
        double s = 0.0;
        std::size_t k = 0;
        for (auto& rw : refs) {
            const auto r = elw::estimar_calibracion(rw.x, Wr, m, c.d_min, c.d_max);   // referencia: el camino del sensor
            rw.d_ref = r.valido ? r.d : NaN;
            if (r.valido) { s += r.d; ++k; }
        }
        g.d_ref_medio = k ? s / static_cast<double>(k) : NaN;
        // (diagnostico) perfil de localidad en W_ref: deriva pareada d^(mm) - d^(m)
        std::size_t m_estrella = m, m_ult = 0;
        double D_ult = NaN, se_ult = NaN;
        for (std::size_t mm = m + 1; mm <= M; ++mm) {
            std::vector<std::vector<double>> diffs;
            for (const auto& dv : d_por_seg) {
                std::vector<double> q;
                for (const auto& d : dv) if (std::isfinite(d[mm]) && std::isfinite(d[m])) q.push_back(d[mm] - d[m]);
                diffs.push_back(std::move(q));
            }
            double D, se, cota;
            if (!resumen(diffs, D, se, cota)) break;
            m_ult = mm; D_ult = D; se_ult = se;
            const bool registrar = mm <= m + 8 || (mm & (mm - 1)) == 0 || mm == M || cota > c.ruido_objetivo;
            if (registrar) { g.perfil_m.push_back(mm); g.perfil_D.push_back(D); g.perfil_se.push_back(se); g.perfil_cota.push_back(cota); }
            if (cota > c.ruido_objetivo) break;
            m_estrella = mm;
        }
        g.m_estrella = m_estrella;
        g.lambda_estrella = 2.0 * elw::PI * static_cast<double>(m_estrella) / static_cast<double>(Wr);
        // (enmienda 2) sesgo de la referencia por la ley (m/n)^2: D(mm) = c (mm^2 - m^2) -> b_ref = c m^2
        if (!(m_ult > m && std::isfinite(D_ult))) { g.motivo = "sin perfil en W_ref: el sesgo de la referencia no se puede acotar"; return false; }
        {
            const double x = static_cast<double>(m_ult) * static_cast<double>(m_ult) - static_cast<double>(m) * static_cast<double>(m);
            const double m2 = static_cast<double>(m) * static_cast<double>(m);
            g.sesgo_referencia = D_ult * m2 / x;
            g.se_sesgo_referencia = se_ult * m2 / x;
        }
        // (enmienda 1) sesgo EN la W candidata: sub-ventanas de W pareadas contra d^(m; W_ref) de su ventana
        auto evaluar = [&](std::size_t W, double& D, double& se) -> bool {
            const std::size_t q = Wr / W;
            std::vector<std::vector<double>> diffs(d_por_seg.size());
            for (const auto& rw : refs) {
                if (!std::isfinite(rw.d_ref)) continue;
                double sv = 0.0;
                std::size_t kv = 0;
                for (std::size_t b = 0; b < q; ++b) {
                    const auto r = elw::estimar_calibracion(rw.x + b * W, W, m, c.d_min, c.d_max);
                    if (r.valido) { sv += r.d; ++kv; }
                }
                if (kv == q) diffs[rw.grupo].push_back(sv / static_cast<double>(kv) - rw.d_ref);
            }
            double cota = NaN;
            const bool ok = resumen(diffs, D, se, cota);
            if (ok) {
                D += g.sesgo_referencia;
                se = std::sqrt(se * se + g.se_sesgo_referencia * g.se_sesgo_referencia);
                cota = std::abs(D) + z * se;
            }
            g.perfil_W.push_back(W); g.perfil_W_D.push_back(ok ? D : NaN); g.perfil_W_se.push_back(ok ? se : NaN); g.perfil_W_cota.push_back(ok ? cota : NaN);
            return ok && cota <= c.ruido_objetivo;
        };
        const std::size_t W0 = std::max(2 * m + 1, c.ventana_interna * (m + 1));
        if (W0 > Wr) { g.motivo = "suelo de ruido y techo de la ventana interna exigen W > W_ref"; return false; }
        std::size_t fallo = 0, W = W0;
        double D = NaN, se = NaN;
        bool paso = false;
        while (true) {
            double Dc, sc;
            if (evaluar(W, Dc, sc)) { paso = true; D = Dc; se = sc; break; }
            fallo = W;
            if (W == Wr) break;
            W = std::min(Wr, std::max(W + 1, static_cast<std::size_t>(std::ceil(static_cast<double>(W) * std::pow(2.0, 0.25)))));
        }
        if (!paso) { g.motivo = "ninguna W <= W_ref tiene el sesgo dentro del suelo de ruido"; return false; }
        if (fallo) {
            std::size_t lo = fallo, hi = W;
            while (hi - lo > 1) {
                const std::size_t mid = lo + (hi - lo) / 2;
                double Dc, sc;
                if (evaluar(mid, Dc, sc)) { hi = mid; D = Dc; se = sc; } else lo = mid;
            }
            W = hi;
        }
        g.W = W;
        g.m = m;
        g.se = elw::se_asintotico(m);
        g.sesgo_W = D;
        g.se_sesgo_W = se;
        g.frecuencia_corte = 2.0 * elw::PI * static_cast<double>(m) / static_cast<double>(W);
        return true;
    };

    // (enmienda 3) suelo derivado del nulo: m <- ceil(m (sd/s*)^2) hasta que la sd medida en el nulo quede <= s*
    std::size_t m = m0;
    for (;;) {
        if (!banda_para(m, f)) return f;
        f.historial_m.push_back(m);
        f.historial_W.push_back(f.W);
        if (!segs_nulo) {
            f.historial_ruido.push_back(NaN);
            f.viable = true;
            f.motivo = "viable (suelo ASINTOTICO: sin nulo de ruido)";
            return f;
        }
        std::vector<double> e;
        for (const auto& g : *segs_nulo)
            for (std::size_t a = 0; a + f.W <= g.size(); a += f.W) {
                const auto r = elw::estimar_calibracion(g.data() + a, f.W, m, c.d_min, c.d_max);
                if (r.valido) e.push_back(r.d);
            }
        f.ventanas_ruido = e.size();
        if (e.size() < nmin) { f.historial_ruido.push_back(NaN); f.motivo = "el nulo de ruido no da 1/alfa ventanas disjuntas de W"; return f; }
        double mu = 0.0;
        for (double v : e) mu += v;
        mu /= static_cast<double>(e.size());
        double ss = 0.0;
        for (double v : e) ss += (v - mu) * (v - mu);
        const double sd = std::sqrt(ss / static_cast<double>(e.size() - 1));
        f.ruido_nulo = sd;
        f.historial_ruido.push_back(sd);
        if (sd <= c.ruido_objetivo) { f.viable = true; f.motivo = "viable"; return f; }
        const std::size_t m_nuevo = std::max(m + 1, static_cast<std::size_t>(std::ceil(static_cast<double>(m) * (sd / c.ruido_objetivo) * (sd / c.ruido_objetivo))));
        if (m_nuevo >= M) { f.motivo = "el ruido medido en el nulo exige m >= Nyquist de W_ref"; return f; }
        m = m_nuevo;
    }
}

}  // namespace ark::metrologia
