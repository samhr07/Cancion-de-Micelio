// lazo_shannon.hpp -- LAZO CAUSAL RECURSIVO DE SHANNON Y FEEDBACK DE VOLATILIDAD (ARK-WIN-01-CXX)
// Reactor ARK / Canción del Micelio -- ARQUITECTO DE DATOS & COMANDANTE
// Corregido por el Jefe de Arquitectos de Código (2026-11-01):
//   - Transaccional: calcular() no toca el estado y confirmar() lo avanza; el motor solo confirma cuando el paso entero es
//     admisible (B4/E8).
//   - Arranque en frío (E2): no se emite d_nu hasta que la historia tenga k_min = min(max(2 s_min, 8), W_innov) retornos reales.
//   - La historia crece de a uno y se recorta a W_innov (historia_retornos.append + pop(0) del Python).
// Paso 3 de CORRECCION_SHANNON_LAPLACE_CXX.md (2026-11-02, órdenes del Comandante y de Antigravity; especificación 1-8):
//   - Mandato P <= 0: un precio no positivo o no finito es una celda enmascarada que no toca el estado.
//   - Escala gaussiana con PRIOR CONJUGADO (inversa-chi2 escalada, media desconocida), en lugar del suelo max(1e-4, tick/P), que
//     en las 26 series diarias no actuaba nunca:
//         sigma^2 = max((nu0 s0^2 + n sh^2) / (nu0 + n - 1), Sheppard),  s0^2 = max(M2 / N, Sheppard)   (nu0 = inf: s0^2)
//     nu0 por trigamma: Var(log b_j) = psi1((B-1)/2) + psi1(nu0/2), b_j las varianzas de bloques NO solapados de B = 8 retornos;
//     el exceso se declara resuelto con z_(1-alfa) * SE (SE = V sqrt(2/(m-1))); si no, nu0 = inf. Sheppard = (tick/P_(t-1))^2/6.
//   - Escala de Laplace con prior conjugado (gamma inversa sobre b), centro en la mediana de la ventana:
//         b = max((alfa0 b0 + n bh) / (alfa0 + n), sqrt(Sheppard / 2)),  b0 = media |x - mediana| de la historia larga
//     alfa0: Var(log b_j) = VLAP + psi1(alfa0), VLAP tabulada (bloques de 8, mediana estimada).
//   - Régimen: Laplace si Q = sd / (MAD / 0.67449) de la historia supera el cuantil 1 - alfa (alfa = 0.05) de su nulo gaussiano
//     para ese N (tabla de Monte Carlo hasta 512, normal asintótica después); exacta mientras la historia quepa en CAP, por
//     bloques de 8 cuando no. Reloj: gauss z^2 / (2 ln 2); Laplace |z_b| / (2 ln 2) (opción (a): la misma tasa media).
//   - E15: W_innov = floor(clamp(2.5 tau_1/2(phi), lo, W_max) + 0.5), phi = (rho10 / rho1)^(1/9) de |r| con historia
//     expansiva y guarda de resolución (Bartlett con el modelo geométrico); sin resolver, W_max.
//   - W_max derivado por activo: el mayor n en [8, W] cuya latencia esperada (trayectoria esperada del perfil de exceso de
//     varianza tras los eventos de la cola 1 - alfa/2) sea menor que su vida útil (mediana de la masa del exceso). Sin
//     1/alfa = 20 eventos completos, 8; con nu0 = inf (la volatilidad no se mueve de forma detectable), W.
//   El oráculo independiente es PRUEBAS_JEFE/escala_prior/oraculo_lazo.py; las tablas, lazo_tablas.hpp (generadas).
//   DIVERGE del Python canónico de Antigravity (motor_ventanas.py) en todo el reloj de Shannon hasta que él lo sincronice.
#pragma once

#include "ark/ventanas/contratos_ventanas.hpp"
#include "ark/ventanas/lazo_tablas.hpp"
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <array>
#include <algorithm>
#include <limits>

namespace ark::ventanas {

namespace detalle_lazo {

inline constexpr double LN2 = 0.693147180559945309417;
inline constexpr double Z1A = 1.6448536269514722;   // Phi^-1(1 - 0.05)
inline constexpr double C_MAD = 0.6744897501960817;
inline constexpr std::size_t B = tablas_lazo::B;
inline constexpr std::size_t H = 64;                  // 2 W_max del paso 2
inline constexpr std::size_t EV_MIN_OBS = 40;         // 1 / (alfa / 2)
inline constexpr std::int64_t EV_MIN_COMPLETOS = 20;  // 1 / alfa
inline constexpr double INF = std::numeric_limits<double>::infinity();

// Suma compensada de Neumaier (la misma en el oráculo).
struct Neumaier {
    double s{0.0}, c{0.0};
    void add(double x) noexcept {
        const double t = s + x;
        if (std::fabs(s) >= std::fabs(x)) c += (s - t) + x;
        else c += (x - t) + s;
        s = t;
    }
    [[nodiscard]] double val() const noexcept { return s + c; }
};

// Trigamma: recurrencia hasta x >= 10 y serie asintótica hasta x^-13.
[[nodiscard]] inline double psi1(double x) noexcept {
    double s = 0.0;
    while (x < 10.0) {
        s += 1.0 / (x * x);
        x += 1.0;
    }
    const double x2 = 1.0 / (x * x);
    return s + 1.0 / x + x2 / 2.0 +
           (x2 / x) * (1.0 / 6.0 + x2 * (-1.0 / 30.0 + x2 * (1.0 / 42.0 + x2 * (-1.0 / 30.0 + x2 * (5.0 / 66.0 + x2 * (-691.0 / 2730.0))))));
}

[[nodiscard]] inline double psi1_inv(double y) noexcept {
    double lo = -30.0, hi = 30.0;
    for (int i = 0; i < 100; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (psi1(std::exp(mid)) > y) lo = mid;
        else hi = mid;
    }
    return std::exp(0.5 * (lo + hi));
}

// Cuantil 1 - alfa del nulo gaussiano de Q para N retornos.
[[nodiscard]] inline double q95_nulo(double N) noexcept {
    const auto& qn = tablas_lazo::Q_N;
    const auto& q95 = tablas_lazo::Q_95;
    if (N > qn.back()) return std::exp(tablas_lazo::Q_MU_INF + Z1A * std::sqrt(tablas_lazo::Q_V_ASINTOTICA / N));
    const std::size_t j = static_cast<std::size_t>(std::lower_bound(qn.begin(), qn.end(), N) - qn.begin());
    if (qn[j] == N) return q95[j];
    const std::size_t a = j - 1, b = j;
    const double w = (std::log(N) - std::log(qn[a])) / (std::log(qn[b]) - std::log(qn[a]));
    return std::exp(std::log(q95[a]) + w * (std::log(q95[b]) - std::log(q95[a])));
}

// Mediana de v[0..k) ya ordenado.
[[nodiscard]] inline double mediana_ordenada(const double* v, std::size_t k) noexcept {
    return (k % 2) ? v[k / 2] : 0.5 * (v[k / 2 - 1] + v[k / 2]);
}

// Mediana de v[0..k) por selección (desordena v), el mismo valor que numpy.median.
[[nodiscard]] inline double mediana_seleccion(double* v, std::size_t k) noexcept {
    std::nth_element(v, v + k / 2, v + k);
    const double a = v[k / 2];
    if (k % 2) return a;
    const double b = *std::max_element(v, v + k / 2);
    return 0.5 * (b + a);
}

}  // namespace detalle_lazo

// =======================================================================
// LazoFeedbackShannonCxx: reloj de Shannon causal con prior conjugado, régimen gauss/Laplace, phi y W_max derivado.
// Cero reservas en el ciclo caliente: todo vive en std::array (unos 4 CAP doubles por activo).
// =======================================================================

template <size_t W, size_t CAP = 4096>
class LazoFeedbackShannonCxx {
    static_assert(W >= 8, "lookback_W debe ser al menos 8 muestras (el piso de W_innov).");
    static_assert(CAP >= 64, "la historia larga necesita al menos 64 retornos.");
    using N_ = detalle_lazo::Neumaier;

public:
    struct Propuesta {
        double d_nu{0.0};
        double ret_k{0.0};
        double sigma{0.0};
        bool primera{false};           // la primera observación válida solo fija val_prev
        bool precio_invalido{false};   // P <= 0 o no finito: celda enmascarada que no toca el estado
        bool laplace{false};           // régimen con el que se midió la sorpresa
    };

    // Mandato P <= 0 (2026-11-02).
    [[nodiscard]] static bool precio_valido(double p) noexcept { return std::isfinite(p) && p > 0.0; }

private:
    static constexpr double LN2 = detalle_lazo::LN2;
    static constexpr double INF = detalle_lazo::INF;
    static constexpr size_t H = detalle_lazo::H;
    static constexpr size_t B = detalle_lazo::B;

    ConfiguracionVentanaCausalCxx<W> config_{};
    double tick_{std::numeric_limits<double>::quiet_NaN()};
    size_t ventana_fija_{0};
    bool modo_feedback_{true};

    // Ventana (los W retornos reales más recientes; la historia usa los len_historia_ últimos)
    std::array<double, W> ring_retornos_{};
    size_t head_historia_{0};
    size_t len_historia_{0};
    size_t llenos_{0};

    double val_prev_{0.0};
    bool inicializado_{false};

    // Welford de toda la historia real
    std::int64_t n_reales_{0};
    double media_{0.0}, M2_{0.0};

    // Historia larga: los últimos CAP retornos (anillo), ordenados, y sus |r| ordenados; y una zona de trabajo
    std::array<double, CAP> hist_{};
    size_t hist_head_{0}, hist_n_{0};
    std::array<double, CAP> ord_{};
    std::array<double, CAP> ord_abs_{};
    std::array<double, CAP> trabajo_{};

    // Eventos y perfil del exceso
    std::array<std::int64_t, H> ev_ini_{};
    std::array<double, H> ev_s0_{};
    size_t ev_n_{0};
    std::array<N_, H + 1> sum_rat_{};
    std::array<std::int64_t, H + 1> cnt_{};
    N_ sum_d0_{};
    std::int64_t cnt_d0_{0};
    std::int64_t cnt_w_max_{-1};
    std::int64_t k_{0};

    // Bloques de B
    std::array<double, B> bloque_{};
    size_t nb_{0};
    N_ SG_{}, SG2_{}, SL_{}, SL2_{}, SQ_{};
    std::int64_t mG_{0}, mL_{0}, mQ_{0};
    double nu0_{INF}, alfa0_{INF};

    // Autocorrelación de |r| (rezagos 1 y 10)
    std::int64_t uN_{0};
    N_ S1_{}, S2_{}, P1_{}, P10_{};
    std::array<double, 10> primeros_{};
    size_t n_primeros_{0};
    std::array<double, 10> ultimos_{};   // anillo de los 10 últimos u
    size_t ult_head_{0}, n_ultimos_{0};
    bool phi_ok_{false};
    double phi_{std::numeric_limits<double>::quiet_NaN()};
    double w_phi_{std::numeric_limits<double>::quiet_NaN()};

    bool laplace_{false};
    size_t w_max_{8};
    size_t tau_v_{0};

    // Telemetría
    size_t w_innov_actual_{0};
    double sigma_actual_{0.0};
    double nu_acumulado_{0.0};

    [[nodiscard]] size_t w_inicial() const noexcept {
        return ventana_fija_ > 0 ? ventana_fija_ : (modo_feedback_ ? std::min<size_t>(8, W) : std::min<size_t>(30, W));
    }
    [[nodiscard]] double ventana_en(size_t i, size_t n) const noexcept {   // i-ésimo (0 = el más viejo) de los n últimos
        return ring_retornos_[(head_historia_ + W - n + i) % W];
    }
    [[nodiscard]] double ultimo_u(size_t atras) const noexcept {   // atras = 1: el último; 10: el de hace 10
        return ultimos_[(ult_head_ + 10 - atras) % 10];
    }
    [[nodiscard]] double b0_historia() const noexcept {
        const double m = detalle_lazo::mediana_ordenada(ord_.data(), hist_n_);
        N_ a;
        for (size_t i = 0; i < hist_n_; ++i) a.add(std::fabs(ord_[i] - m));
        return a.val() / static_cast<double>(hist_n_);
    }

public:
    constexpr LazoFeedbackShannonCxx() noexcept = default;

    explicit LazoFeedbackShannonCxx(const ConfiguracionVentanaCausalCxx<W>& config) noexcept { configurar(config); }

    // tick NaN = el de la configuración
    void configurar(const ConfiguracionVentanaCausalCxx<W>& config, double tick = std::numeric_limits<double>::quiet_NaN()) noexcept {
        config_ = config;
        tick_ = std::isnan(tick) ? config.tick_size : tick;
        ventana_fija_ = config.ventana_innovacion_shannon;
        modo_feedback_ = config.lazo_feedback_shannon && (config.ventana_innovacion_shannon == 0);
        reset();
    }

    [[nodiscard]] double tick() const noexcept { return tick_; }

    // d_nu del paso (val_k, mask_k) con el estado anterior (F_(t-1)), sin modificar nada.
    [[nodiscard]] Propuesta calcular(double val_k, bool mask_k) const noexcept {
        Propuesta p{};
        if (!precio_valido(val_k)) {
            p.precio_invalido = true;
            p.sigma = sigma_actual_;
            return p;
        }
        if (!inicializado_) {
            p.primera = true;
            return p;
        }
        const double ret_raw = std::log(val_k) - std::log(val_prev_);
        p.ret_k = (config_.congelar_en_huecos && mask_k) ? 0.0 : ret_raw;
        p.sigma = sigma_actual_;
        p.laplace = laplace_;
        if (mask_k) return p;
        const double qs = tick_ / val_prev_;
        const double sh = qs * qs / 6.0;
        const size_t n = len_historia_;
        const size_t k_min = std::min(config_.k_min_shannon(), w_innov_actual_);
        if (std::abs(p.ret_k) < 1e-12 || n < k_min) return p;
        const double r = p.ret_k, nn = static_cast<double>(n);
        if (!laplace_) {
            N_ s;
            for (size_t i = 0; i < n; ++i) s.add(ventana_en(i, n));
            const double mu = s.val() / nn;
            N_ q;
            for (size_t i = 0; i < n; ++i) {
                const double d = ventana_en(i, n) - mu;
                q.add(d * d);
            }
            const double v = q.val() / nn;
            const double s02 = std::max(M2_ / static_cast<double>(n_reales_), sh);
            const double s2 = std::isinf(nu0_) ? s02 : (nu0_ * s02 + nn * v) / (nu0_ + nn - 1.0);
            const double sg = std::sqrt(std::max(s2, sh));
            const double z = (r - mu) / sg;
            p.d_nu = (z * z) / (2.0 * LN2);
            p.sigma = sg;
        } else {
            std::array<double, W> x{};
            for (size_t i = 0; i < n; ++i) x[i] = ventana_en(i, n);
            std::sort(x.begin(), x.begin() + static_cast<std::ptrdiff_t>(n));
            const double m = detalle_lazo::mediana_ordenada(x.data(), n);
            N_ s;
            for (size_t i = 0; i < n; ++i) s.add(std::fabs(ventana_en(i, n) - m));
            const double bh = s.val() / nn;
            const double bfl = std::sqrt(sh / 2.0);
            const double b0 = std::max(b0_historia(), bfl);
            double b = std::isinf(alfa0_) ? b0 : (alfa0_ * b0 + nn * bh) / (alfa0_ + nn);
            b = std::max(b, bfl);
            p.d_nu = std::fabs((r - m) / b) / (2.0 * LN2);
            p.sigma = std::sqrt(2.0) * b;
        }
        return p;
    }

    // Avanza el estado con el paso ya calculado (mismos val_k, mask_k que calcular()).
    void confirmar(double val_k, bool mask_k, const Propuesta& p) noexcept {
        if (p.precio_invalido) return;   // d_nu = 0 y ningún estado se mueve
        if (p.primera) {
            val_prev_ = val_k;
            inicializado_ = true;
            return;
        }
        nu_acumulado_ += p.d_nu;
        if (mask_k) {   // Axioma 4e: las filas imputadas no entran en ninguna historia
            val_prev_ = val_k;
            return;
        }
        sigma_actual_ = p.sigma;
        const double r = p.ret_k;
        const double qs = tick_ / val_prev_;
        const double sh = qs * qs / 6.0;
        const double u = std::fabs(r);
        const double s0_pre = n_reales_ > 0 ? M2_ / static_cast<double>(n_reales_) : 0.0;
        // a. evento
        bool evento = false;
        if (hist_n_ >= detalle_lazo::EV_MIN_OBS && s0_pre > 0.0) {
            const double kq = static_cast<double>(hist_n_);
            const size_t iq = std::min(hist_n_ - 1, static_cast<size_t>(std::ceil(0.975 * kq)) - 1);
            evento = u >= ord_abs_[iq];
        }
        // b. perfil de los eventos activos (en orden de alta)
        size_t quedan = 0;
        for (size_t e = 0; e < ev_n_; ++e) {
            const size_t h = static_cast<size_t>(k_ - ev_ini_[e]);
            sum_rat_[h].add(r * r / ev_s0_[e]);
            ++cnt_[h];
            if (h != H) {
                ev_ini_[quedan] = ev_ini_[e];
                ev_s0_[quedan] = ev_s0_[e];
                ++quedan;
            }
        }
        ev_n_ = quedan;
        // c. alta
        if (evento) {
            ev_ini_[ev_n_] = k_;
            ev_s0_[ev_n_] = s0_pre;
            ++ev_n_;
            sum_d0_.add(r * r / s0_pre);
            ++cnt_d0_;
        }
        // d. Welford
        ++n_reales_;
        const double d = r - media_;
        media_ += d / static_cast<double>(n_reales_);
        M2_ += d * (r - media_);
        // e. historia larga
        insertar_historia(r);
        // f. bloques
        bloque_[nb_++] = r;
        bool bloque_nuevo = false;
        if (nb_ == B) {
            bloque_nuevo = true;
            cerrar_bloque(sh);
            nb_ = 0;
        }
        // g. autocorrelación de |r|
        actualizar_phi(u);
        // h. régimen
        actualizar_regimen();
        // i. W_max
        if (modo_feedback_ && cnt_[H] >= detalle_lazo::EV_MIN_COMPLETOS && (bloque_nuevo || cnt_[H] != cnt_w_max_)) {
            cnt_w_max_ = cnt_[H];
            recalcular_w_max();
        }
        // j. W_innov
        if (ventana_fija_ > 0) {
            w_innov_actual_ = ventana_fija_;
        } else if (!modo_feedback_) {
            w_innov_actual_ = std::min<size_t>(30, W);
        } else {
            const size_t lo = std::min(std::max(config_.s_min * 2, static_cast<size_t>(8)), w_max_);
            w_innov_actual_ = phi_ok_ ? static_cast<size_t>(std::floor(std::min(std::max(w_phi_, static_cast<double>(lo)),
                                                                                 static_cast<double>(w_max_)) + 0.5))
                                      : w_max_;
        }
        // k. ventana
        ring_retornos_[head_historia_] = r;
        head_historia_ = (head_historia_ + 1) % W;
        llenos_ = std::min(llenos_ + 1, W);
        len_historia_ = std::min(len_historia_ + 1, w_innov_actual_);
        ++k_;
        val_prev_ = val_k;
    }

    // Comodidad: calcular y confirmar de una vez.
    [[nodiscard]] double paso(double val_k, bool mask_k) noexcept {
        const Propuesta p = calcular(val_k, mask_k);
        confirmar(val_k, mask_k, p);
        return p.d_nu;
    }

    // Telemetría causal
    [[nodiscard]] constexpr size_t w_innov_actual() const noexcept { return w_innov_actual_; }
    [[nodiscard]] double tau_half_actual() const noexcept { return phi_ok_ ? w_phi_ / 2.5 : std::numeric_limits<double>::quiet_NaN(); }
    [[nodiscard]] constexpr double sigma_actual() const noexcept { return sigma_actual_; }
    [[nodiscard]] constexpr double rho_actual() const noexcept { return phi_; }
    [[nodiscard]] constexpr double nu_acumulado() const noexcept { return nu_acumulado_; }
    [[nodiscard]] constexpr size_t conteo_reales() const noexcept { return static_cast<size_t>(n_reales_); }
    [[nodiscard]] constexpr size_t longitud_historia() const noexcept { return len_historia_; }
    [[nodiscard]] constexpr bool inicializado() const noexcept { return inicializado_; }
    [[nodiscard]] constexpr bool regimen_laplace() const noexcept { return laplace_; }
    [[nodiscard]] constexpr size_t w_max() const noexcept { return w_max_; }
    [[nodiscard]] constexpr double nu0() const noexcept { return nu0_; }
    [[nodiscard]] constexpr double alfa0() const noexcept { return alfa0_; }
    [[nodiscard]] constexpr bool phi_resuelto() const noexcept { return phi_ok_; }
    [[nodiscard]] constexpr size_t vida_util() const noexcept { return tau_v_; }

    // Reinicio limpio de estado sin desasignaciones
    void reset() noexcept {
        ring_retornos_.fill(0.0);
        head_historia_ = len_historia_ = llenos_ = 0;
        val_prev_ = 0.0;
        inicializado_ = false;
        n_reales_ = 0;
        media_ = M2_ = 0.0;
        hist_head_ = hist_n_ = 0;
        ev_n_ = 0;
        sum_rat_.fill(N_{});
        cnt_.fill(0);
        sum_d0_ = N_{};
        cnt_d0_ = 0;
        cnt_w_max_ = -1;
        k_ = 0;
        nb_ = 0;
        SG_ = SG2_ = SL_ = SL2_ = SQ_ = N_{};
        mG_ = mL_ = mQ_ = 0;
        nu0_ = alfa0_ = INF;
        uN_ = 0;
        S1_ = S2_ = P1_ = P10_ = N_{};
        n_primeros_ = ult_head_ = n_ultimos_ = 0;
        phi_ok_ = false;
        phi_ = w_phi_ = std::numeric_limits<double>::quiet_NaN();
        laplace_ = false;
        w_max_ = ventana_fija_ > 0 ? ventana_fija_ : std::min<size_t>(8, W);
        tau_v_ = 0;
        w_innov_actual_ = w_inicial();
        sigma_actual_ = 0.0;
        nu_acumulado_ = 0.0;
    }

private:
    void insertar_historia(double r) noexcept {
        auto insertar = [](std::array<double, CAP>& a, size_t n, double x) {   // insort_right
            const size_t pos = static_cast<size_t>(std::upper_bound(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(n), x) - a.begin());
            for (size_t i = n; i > pos; --i) a[i] = a[i - 1];
            a[pos] = x;
        };
        auto borrar = [](std::array<double, CAP>& a, size_t n, double x) {   // borra la primera aparición (bisect_left)
            const size_t pos = static_cast<size_t>(std::lower_bound(a.begin(), a.begin() + static_cast<std::ptrdiff_t>(n), x) - a.begin());
            for (size_t i = pos; i + 1 < n; ++i) a[i] = a[i + 1];
        };
        if (hist_n_ == CAP) {
            const double viejo = hist_[hist_head_];
            borrar(ord_, hist_n_, viejo);
            borrar(ord_abs_, hist_n_, std::fabs(viejo));
            --hist_n_;
        }
        insertar(ord_, hist_n_, r);
        insertar(ord_abs_, hist_n_, std::fabs(r));
        hist_[hist_head_] = r;
        hist_head_ = (hist_head_ + 1) % CAP;
        ++hist_n_;
    }

    void cerrar_bloque(double sh) noexcept {
        const double bd = static_cast<double>(B);
        N_ s;
        for (size_t i = 0; i < B; ++i) s.add(bloque_[i]);
        const double bm = s.val() / bd;
        N_ q;
        for (size_t i = 0; i < B; ++i) q.add((bloque_[i] - bm) * (bloque_[i] - bm));
        const double bv = q.val() / bd;
        const double lg = std::log(std::max(bv, sh));
        SG_.add(lg);
        SG2_.add(lg * lg);
        ++mG_;
        std::array<double, B> xs = bloque_;
        std::sort(xs.begin(), xs.end());
        const double med = detalle_lazo::mediana_ordenada(xs.data(), B);
        N_ a;
        for (size_t i = 0; i < B; ++i) a.add(std::fabs(bloque_[i] - med));
        const double bb = a.val() / bd;
        const double lb = std::log(std::max(bb, std::sqrt(sh / 2.0)));
        SL_.add(lb);
        SL2_.add(lb * lb);
        ++mL_;
        std::array<double, B> dv{};
        for (size_t i = 0; i < B; ++i) dv[i] = std::fabs(bloque_[i] - med);
        std::sort(dv.begin(), dv.end());
        const double mad = detalle_lazo::mediana_ordenada(dv.data(), B) / detalle_lazo::C_MAD;
        if (mad > 0.0) {
            SQ_.add(std::log(std::sqrt(bv) / mad));
            ++mQ_;
        }
        if (mG_ >= 3) {
            const double m = static_cast<double>(mG_);
            const double sg = SG_.val();
            const double V = (SG2_.val() - sg * sg / m) / (m - 1.0);
            const double E = V - detalle_lazo::psi1((bd - 1.0) / 2.0);
            nu0_ = E > detalle_lazo::Z1A * V * std::sqrt(2.0 / (m - 1.0)) ? 2.0 * detalle_lazo::psi1_inv(E) : INF;
            const double ml = static_cast<double>(mL_);
            const double sl = SL_.val();
            const double VL = (SL2_.val() - sl * sl / ml) / (ml - 1.0);
            const double EL = VL - tablas_lazo::VLAP_MEDIANA;
            alfa0_ = EL > detalle_lazo::Z1A * VL * std::sqrt(2.0 / (ml - 1.0)) ? detalle_lazo::psi1_inv(EL) : INF;
        }
    }

    void actualizar_phi(double u) noexcept {
        if (uN_ >= 1) P1_.add(u * ultimo_u(1));
        if (uN_ >= 10) P10_.add(u * ultimo_u(10));
        ++uN_;
        S1_.add(u);
        S2_.add(u * u);
        if (n_primeros_ < 10) primeros_[n_primeros_++] = u;
        ultimos_[ult_head_] = u;
        ult_head_ = (ult_head_ + 1) % 10;
        n_ultimos_ = std::min<size_t>(n_ultimos_ + 1, 10);
        phi_ok_ = false;
        if (uN_ < 11) return;
        const double Nn = static_cast<double>(uN_);
        const double s1 = S1_.val();
        const double ub = s1 / Nn;
        const double c0 = S2_.val() / Nn - ub * ub;
        auto ck = [&](size_t k, const N_& P) {
            N_ a, b;
            for (size_t i = 0; i < k; ++i) a.add(primeros_[i]);
            for (size_t i = 0; i < k; ++i) b.add(ultimo_u(k - i));   // los k últimos en orden temporal
            const double A = s1 - a.val(), Bk = s1 - b.val();
            return (P.val() - ub * (A + Bk) + (Nn - static_cast<double>(k)) * ub * ub) / Nn;
        };
        const double c1 = ck(1, P1_), c10 = ck(10, P10_);
        if (!(c0 > 0.0)) return;
        const double r1 = c1 / c0, r10 = c10 / c0;
        if (!(r1 > 0.0 && r10 > 0.0 && r10 < r1)) return;
        const double phi = std::pow(r10 / r1, 1.0 / 9.0);
        double acc = 1.0;
        for (int j = 1; j < 10; ++j) {
            const double rj = r1 * std::pow(phi, static_cast<double>(j - 1));
            acc += 2.0 * rj * rj;
        }
        if (r10 > detalle_lazo::Z1A * std::sqrt(acc / Nn)) {
            phi_ok_ = true;
            phi_ = phi;
            w_phi_ = 2.5 * (-LN2 / std::log(phi));
        }
    }

    void actualizar_regimen() noexcept {
        if (static_cast<size_t>(n_reales_) <= CAP) {
            if (n_reales_ < 8) return;
            const double med = detalle_lazo::mediana_ordenada(ord_.data(), hist_n_);
            for (size_t i = 0; i < hist_n_; ++i) trabajo_[i] = std::fabs(ord_[i] - med);
            const double mad = detalle_lazo::mediana_seleccion(trabajo_.data(), hist_n_) / detalle_lazo::C_MAD;
            if (mad > 0.0) {
                const double Q = std::sqrt(M2_ / static_cast<double>(n_reales_)) / mad;
                laplace_ = Q >= detalle_lazo::q95_nulo(static_cast<double>(n_reales_));
            }
        } else if (mQ_ >= 2) {
            const double m = static_cast<double>(mQ_);
            const double z = (SQ_.val() / m - tablas_lazo::QBLOQ_MEDIA) / (tablas_lazo::QBLOQ_SD / std::sqrt(m));
            laplace_ = z >= detalle_lazo::Z1A;
        }
    }

    void recalcular_w_max() noexcept {
        if (std::isinf(nu0_)) {   // la volatilidad no se mueve de forma detectable: no hay régimen que perseguir
            w_max_ = W;
            return;
        }
        std::array<double, H> ex{}, pos{};
        N_ tot_n;
        double exmax = -INF;
        for (size_t h = 1; h <= H; ++h) {
            ex[h - 1] = sum_rat_[h].val() / static_cast<double>(cnt_[h]) - 1.0;
            pos[h - 1] = std::max(ex[h - 1], 0.0);
            tot_n.add(pos[h - 1]);
            exmax = std::max(exmax, ex[h - 1]);
        }
        const double tot = tot_n.val();
        if (tot <= 0.0 || exmax <= 0.0) {
            w_max_ = W;
            return;
        }
        N_ acu;
        size_t tau_v = H;
        for (size_t h = 0; h < H; ++h) {
            acu.add(pos[h]);
            if (acu.val() >= 0.5 * tot) {
                tau_v = h + 1;
                break;
            }
        }
        const double d0 = sum_d0_.val() / static_cast<double>(cnt_d0_);
        const double umbral = 1.0 + exmax / 2.0;
        size_t mejor = 8;
        for (size_t n = 8; n <= W; ++n) {
            const double nd = static_cast<double>(n);
            size_t lat = H + 1;
            for (size_t h = 1; h <= H; ++h) {
                N_ acc;
                for (std::int64_t idx = static_cast<std::int64_t>(h) - static_cast<std::int64_t>(n); idx < static_cast<std::int64_t>(h); ++idx)
                    acc.add(idx < 0 ? 1.0 : (idx == 0 ? d0 : std::max(1.0 + ex[static_cast<size_t>(idx - 1)], 0.0)));
                const double vw = acc.val() / nd;
                const double est = (nu0_ + nd * vw) / (nu0_ + nd - 1.0);
                if (est >= umbral) {
                    lat = h;
                    break;
                }
            }
            if (lat < tau_v) mejor = n;
        }
        w_max_ = mejor;
        tau_v_ = tau_v;
    }
};

} // namespace ark::ventanas
