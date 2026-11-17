// canal_nivel.hpp -- canal de nivel de una escala (ARK-METRO-01): GLR de Lorden para un cambio de media de tamano desconocido,
// cota de Lorden, tau50 medido, y el resolvente dual del Bloque Fractal (certificado del puente + Condat L1 / Tikhonov L2).
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04. Prerregistro: PRUEBAS_JEFE/metrologia/EXPECTATIVAS.md §1.
//
// GLR (Lorden 1971; Siegmund y Venkatraman 1995) con varianza de largo plazo s0^2 tau_int:
//     G_t = max_{1 <= l <= L} (sum de los ultimos l z)^2 / (2 l tau_int),   z = (y - mu0)/s0,
// reinicio tras cada alarma y tras cada barra no valida (sin imputar). h por ARL0 sobre el nulo de la propia escala.
// Lo REUTILIZADO del Bloque Fractal, sin tocarlo: maximo_puente, RegularizadorDualFractal (Condat exacto, TDMA + Morozov) y la
// regla consolidada de la norma (certificado del puente).
//
// Ciclo caliente (paso): cero reservas. La calibracion es fuera de linea y reserva (std::vector).
#pragma once

#include "ark/fractal/regularizador_dual_fractal.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace ark::metrologia {

template <std::size_t L>
class DetectorGLR {
    static_assert(L >= 2, "L >= 2");
    std::array<double, L> z_{};
    std::size_t cab_ = 0, n_ = 0;
    double mu_ = NaN, s_ = NaN, tau_ = NaN, h_ = NaN;

public:
    struct Paso { bool alarma = false; double g = NaN; std::size_t l = 0; double suma = NaN; double t_pendiente = NaN; };
    void configurar(double mu, double s, double tau, double h) noexcept { mu_ = mu; s_ = s; tau_ = tau; h_ = h; reiniciar(); }
    void reiniciar() noexcept { n_ = 0; cab_ = 0; }
    [[nodiscard]] std::size_t lleno() const noexcept { return n_; }

    Paso paso(double y) noexcept {
        Paso p{};
        z_[(cab_ + n_) % L] = (y - mu_) / s_;
        if (n_ < L) ++n_; else cab_ = (cab_ + 1) % L;
        double S = 0.0, gmax = -1.0, smax = 0.0;
        std::size_t lmax = 0;
        for (std::size_t l = 1; l <= n_; ++l) {
            S += z_[(cab_ + n_ - l) % L];
            const double g = S * S / (2.0 * static_cast<double>(l) * tau_);
            if (g > gmax) { gmax = g; lmax = l; smax = S; }
        }
        p.g = gmax; p.l = lmax; p.suma = smax;
        if (gmax >= h_) {
            p.alarma = true;
            // forma del cambio: pendiente MCO de los l* ultimos z contra el tiempo, con la varianza de largo plazo tau
            if (lmax >= 3) {
                const double kb = 0.5 * (static_cast<double>(lmax) + 1.0), zb = smax / static_cast<double>(lmax);
                double sxx = 0.0, sxy = 0.0;
                for (std::size_t k = 1; k <= lmax; ++k) {   // k = 1 el mas antiguo de los l*
                    const double zk = z_[(cab_ + n_ - lmax + k - 1) % L];
                    sxx += (k - kb) * (k - kb);
                    sxy += (k - kb) * (zk - zb);
                }
                p.t_pendiente = (sxy / sxx) / std::sqrt(tau_ / sxx);
            }
            reiniciar();
        }
        return p;
    }
};

namespace detalle_nivel {

// Alarmas del GLR (con reinicio) sobre segmentos contiguos; cada segmento empieza con la memoria vacia.
template <std::size_t L>
inline std::size_t alarmas(const std::vector<std::vector<double>>& segs, double mu, double s, double tau, double h, std::size_t& total) {
    DetectorGLR<L> d;
    std::size_t al = 0;
    total = 0;
    for (const auto& g : segs) {
        d.configurar(mu, s, tau, h);
        for (double v : g) al += d.paso(v).alarma;
        total += g.size();
    }
    return al;
}

struct ResultadoH { double h = NaN; bool extrapolado = false; std::size_t alarmas = 0; };

// h con ARL(h) = arl0 sobre el nulo. ARL(h) se lee con >= min_al alarmas; mas alla, extrapolacion log-lineal (declarada) con 8
// puntos en [h_b/2, h_b], h_b el mayor h legible. `eval(h, al)` devuelve ARL(h) (pasos / alarmas) y escribe las alarmas en `al`:
// el detector disjunto y el deslizante (detector_deslizante.hpp) usan la MISMA busqueda con su propia construccion del nulo.
template <class Eval>
inline ResultadoH buscar_h_con(Eval&& eval, double arl0, std::size_t min_al) {
    ResultadoH r{};
    double h_prev = 0.0, h = 1.0;
    std::size_t al = 0;
    double a = eval(h, al);
    // con memoria fuerte (tau_int > L) el estadistico es chico: se baja h hasta tener alarmas legibles (antes se rendia en h = 1)
    while (al < min_al) {
        h *= 0.5;
        if (h < 1e-8) return r;   // ni asi hay alarmas legibles: historia insuficiente
        a = eval(h, al);
    }
    while (a < arl0) {
        h_prev = h;
        h *= 1.5;
        a = eval(h, al);
        if (al < min_al) break;
        if (h > 1e4) return r;
    }
    if (al >= min_al && a >= arl0) {   // biseccion directa en log ARL
        double lo = h_prev, hi = h;
        for (int it = 0; it < 25; ++it) {
            const double mid = 0.5 * (lo + hi);
            std::size_t a2 = 0;
            const double am = eval(mid, a2);
            (am >= arl0 ? hi : lo) = mid;
        }
        r.h = hi;
        eval(r.h, r.alarmas);
        return r;
    }
    // frontera de legibilidad
    double lo = h_prev, hi = h;
    for (int it = 0; it < 25; ++it) {
        const double mid = 0.5 * (lo + hi);
        std::size_t a2 = 0;
        eval(mid, a2);
        (a2 >= min_al ? lo : hi) = mid;
    }
    const double hb = lo;
    double mx = 0, my = 0, sxx = 0, sxy = 0;
    std::array<double, 8> xs{}, ys{};
    for (int i = 0; i < 8; ++i) {
        xs[i] = hb * (0.5 + 0.5 * i / 7.0);
        std::size_t a2 = 0;
        ys[i] = std::log(eval(xs[i], a2));
        mx += xs[i]; my += ys[i];
    }
    mx /= 8; my /= 8;
    for (int i = 0; i < 8; ++i) { sxx += (xs[i] - mx) * (xs[i] - mx); sxy += (xs[i] - mx) * (ys[i] - my); }
    if (!(sxx > 0.0) || !(sxy > 0.0)) return r;
    const double b = sxy / sxx;
    r.h = hb + (std::log(arl0) - (my + b * (hb - mx))) / b;
    r.extrapolado = true;
    eval(r.h, r.alarmas);
    return r;
}

template <std::size_t L>
inline ResultadoH buscar_h(const std::vector<std::vector<double>>& segs, double mu, double s, double tau, double arl0, std::size_t min_al) {
    std::size_t total = 0;
    auto eval = [&](double h, std::size_t& al) {
        al = alarmas<L>(segs, mu, s, tau, h, total);
        return al ? static_cast<double>(total) / static_cast<double>(al) : std::numeric_limits<double>::infinity();
    };
    return buscar_h_con(eval, arl0, min_al);
}

}  // namespace detalle_nivel

struct ConfiguracionNivel {
    double arl0_pasos = 0.0;          // OBLIGATORIO: ARL0 de esta escala, en pasos de la escala
    double alfa = ALFA_PROYECTO;
    std::size_t alarmas_minimas = 30; // ARL(h) se lee con >= 30 alarmas (declarado, el del Bloque Fractal)
    bool medir_tau50 = true;
    std::size_t replicas_tau50 = 64;  // declarado (replicas_salto del Bloque Fractal)
    std::uint64_t semilla = 20261104;
};

struct SalidaNivel {
    bool emitida = false;
    double g = NaN;
    EventoAlarma alarma{};
    NormaRegularizador norma = NormaRegularizador::SIN_DATO;
    double borde = NaN;           // x[N-1] del resolvente
    double derivada = NaN;        // x[N-1] - x[N-2]
    int ultimo_salto = -1;        // en L1: indice (0..N-1) del ultimo cambio de nivel de la solucion; -1 si no hay
};

template <std::size_t L = 120, std::size_t N_REG = 120>
class CanalNivel {
    static_assert(N_REG >= 8, "N_REG >= 8");
    FichaEscala f_{};
    DetectorGLR<L> glr_{};
    ark::fractal::RegularizadorDualFractal<N_REG> reg_{};
    std::array<double, N_REG> anillo_{}, ventana_{};
    std::size_t cab_ = 0, n_ = 0;

public:
    [[nodiscard]] const FichaEscala& ficha() const noexcept { return f_; }
    void fijar_ficha(const FichaEscala& f) noexcept { f_ = f; reiniciar_flujo(); }
    void reiniciar_flujo() noexcept {
        glr_.configurar(f_.nulo.mu0, f_.nulo.s0, f_.nulo.tau_int, f_.det.h);
        n_ = cab_ = 0;
    }

    // Calibra con segmentos contiguos de y (la historia de esta escala). Fuera del ciclo caliente.
    const FichaEscala& calibrar(int j, std::int64_t m, const std::vector<std::vector<double>>& segs, const ConfiguracionNivel& cfg) {
        FichaEscala f{};
        f.j = j; f.m = m;
        f.det.L = L; f.det.alfa = cfg.alfa; f.det.arl0_pasos = cfg.arl0_pasos;
        const auto mo = est::momentos(segs);
        f.nulo.pasos = mo.n; f.nulo.mu0 = mo.mu; f.nulo.s0 = mo.sd;
        // ventanas disjuntas para lambda_TV
        std::vector<double> mx;
        for (const auto& g : segs)
            for (std::size_t a = 0; a + N_REG <= g.size(); a += N_REG) mx.push_back(ark::fractal::maximo_puente(g.data() + a, N_REG));
        f.det.ventanas_lambda = mx.size();
        const std::size_t minimas = static_cast<std::size_t>(std::ceil(1.0 / cfg.alfa));
        if (mx.size() < minimas || !(mo.sd > 0.0)) { f.motivo = "historia corta: menos de ceil(1/alfa) ventanas disjuntas de N_REG"; return f_ = f; }
        f.det.lambda_tv = est::cuantil(mx, 1.0 - cfg.alfa);
        const auto ti = est::tau_integrado(segs, mo.mu, mo.sd, cfg.alfa, std::min<std::size_t>(mo.n / 4, 10 * L));
        f.nulo.tau_int = ti.tau; f.nulo.acf1 = ti.acf1; f.nulo.rezago_corte = ti.K;
        if (!(ti.tau > 0.0)) { f.motivo = "tau_int <= 0"; return f_ = f; }
        if (!(cfg.arl0_pasos > 1.0)) { f.motivo = "ARL0 de la escala <= 1 paso"; return f_ = f; }
        const auto rh = detalle_nivel::buscar_h<L>(segs, mo.mu, mo.sd, ti.tau, cfg.arl0_pasos, cfg.alarmas_minimas);
        f.det.h = rh.h; f.det.h_extrapolado = rh.extrapolado; f.det.alarmas_nulo = rh.alarmas;
        if (!std::isfinite(rh.h)) { f.motivo = "h no finito: el nulo no da alarmas legibles"; return f_ = f; }
        for (std::size_t i = 0; i < DELTA_REPORTE.size(); ++i) {
            const double d = DELTA_REPORTE[i];
            f.det.t_lorden[i] = std::log(cfg.arl0_pasos) * 2.0 * ti.tau / (d * d);
            f.det.tau50[i] = NaN;
            f.det.fraccion_sin_detectar[i] = NaN;
        }
        f.detectable = true;
        f.motivo = "calibrada";
        if (cfg.medir_tau50) medir_tau50(f, segs, cfg);
        f_ = f;
        reiniciar_flujo();
        return f_;
    }

    // Ciclo caliente: una barra de la escala. t_nativo = indice nativo de su ultima barra nativa.
    SalidaNivel paso(double y, bool valido, std::int64_t t_nativo) noexcept {
        SalidaNivel s{};
        s.emitida = true;
        if (!f_.detectable) return s;
        if (!valido || !std::isfinite(y)) { glr_.reiniciar(); n_ = cab_ = 0; return s; }
        const auto p = glr_.paso(y);
        s.g = p.g;
        if (p.alarma) {
            s.alarma.hay = true;
            s.alarma.signo = p.suma > 0.0 ? 1 : -1;
            s.alarma.t = t_nativo;
            s.alarma.l = p.l;
            s.alarma.t_cambio = t_nativo - static_cast<std::int64_t>(p.l) * f_.m + 1;
            s.alarma.amplitud = p.suma / static_cast<double>(p.l) * f_.nulo.s0;
            s.alarma.se_amplitud = f_.nulo.s0 * std::sqrt(f_.nulo.tau_int / static_cast<double>(p.l));
            s.alarma.g = p.g;
            s.alarma.t_pendiente = p.t_pendiente;
        }
        anillo_[(cab_ + n_) % N_REG] = y;
        if (n_ < N_REG) ++n_; else cab_ = (cab_ + 1) % N_REG;
        if (n_ < N_REG) return s;
        for (std::size_t i = 0; i < N_REG; ++i) ventana_[i] = anillo_[(cab_ + i) % N_REG];
        const bool l1 = ark::fractal::maximo_puente(ventana_.data(), N_REG) > f_.det.lambda_tv;   // regla consolidada
        const auto r = l1 ? reg_.l1(ventana_.data(), N_REG, f_.det.lambda_tv) : reg_.l2(ventana_.data(), N_REG, f_.nulo.s0);
        const auto& x = reg_.solucion();
        s.norma = l1 ? NormaRegularizador::L1_TV : NormaRegularizador::L2_TIKHONOV;
        s.borde = r.borde;
        s.derivada = x[N_REG - 1] - x[N_REG - 2];
        if (l1) for (std::size_t i = N_REG - 1; i >= 1; --i) if (x[i] != x[i - 1]) { s.ultimo_salto = static_cast<int>(i); break; }
        return s;
    }

private:
    // tau50: escalon delta*s0 inyectado en la propia historia nula (conserva su memoria y su marginal); retardo = primera alarma
    // de signo + en o tras t0, en pasos de la escala. Horizonte: 20 t* o el fin del segmento.
    void medir_tau50(FichaEscala& f, const std::vector<std::vector<double>>& segs, const ConfiguracionNivel& cfg) const {
        std::vector<std::size_t> largos;
        for (std::size_t i = 0; i < segs.size(); ++i) if (segs[i].size() >= 2 * L) largos.push_back(i);
        if (largos.empty()) return;
        std::mt19937_64 rng(cfg.semilla ^ static_cast<std::uint64_t>(f.j));
        for (std::size_t k = 0; k < DELTA_REPORTE.size(); ++k) {
            const double sh = DELTA_REPORTE[k] * f.nulo.s0;
            std::vector<double> taus;
            for (std::size_t r = 0; r < cfg.replicas_tau50; ++r) {
                const auto& g = segs[largos[std::uniform_int_distribution<std::size_t>(0, largos.size() - 1)(rng)]];
                const std::size_t t0 = std::uniform_int_distribution<std::size_t>(L, g.size() - L)(rng);
                const std::size_t fin = std::min(g.size(), t0 + static_cast<std::size_t>(std::ceil(20.0 * f.det.t_lorden[k])) + 1);
                DetectorGLR<L> d;
                d.configurar(f.nulo.mu0, f.nulo.s0, f.nulo.tau_int, f.det.h);
                double tau = std::numeric_limits<double>::infinity();
                for (std::size_t i = t0 - L; i < fin; ++i) {
                    const auto p = d.paso(g[i] + (i >= t0 ? sh : 0.0));
                    if (p.alarma && i >= t0 && p.suma > 0.0) { tau = static_cast<double>(i - t0 + 1); break; }
                }
                taus.push_back(tau);
            }
            std::size_t nd = 0;
            for (double t : taus) nd += !std::isfinite(t);
            f.det.fraccion_sin_detectar[k] = static_cast<double>(nd) / static_cast<double>(taus.size());
            std::sort(taus.begin(), taus.end());
            const std::size_t m = taus.size();
            f.det.tau50[k] = (m % 2) ? taus[m / 2] : 0.5 * (taus[m / 2 - 1] + taus[m / 2]);
        }
    }
};

}  // namespace ark::metrologia
