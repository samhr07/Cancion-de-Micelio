// detector_deslizante.hpp -- modo DESLIZANTE del operador delta (ARK-METRO-01): la deteccion no espera a que cierre la barra de la
// escala ni la ventana de temperatura. Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-05, por orden del Comandante ("haga el
// detector deslizante, pero no descarte el que va por ventanas independientes"). Prerregistro: PRUEBAS_JEFE/deslizante/EXPECTATIVAS_DESLIZANTE.md.
//
// Que se reutiliza (no se rehace):
//   - Acumulador / combinar / valor de la escalera disjunta: la escalera deslizante es el MISMO arbol binario en cada paso nativo,
//     b_j(t) = combinar(b_{j-1}(t - 2^{j-1}), b_{j-1}(t)); en los instantes alineados coincide bit a bit con la disjunta.
//   - el nulo (mu0, s0, tau_int) ya estimado por los canales disjuntos: elementos separados D son los disjuntos.
//   - la busqueda de h por ARL0 (detalle_nivel::buscar_h_con), con la construccion deslizante del nulo.
//   - RingBufferEspejoCxx de ark/ventanas: la ventana contigua de las ultimas T* barras del canal de temperatura.
//
// Detector espaciado: G_t = max_{1<=l<=L} (sum_{i<l} z(t - i D))^2 / (2 l tau), z = (y - mu0)/s0, evaluado en CADA entrada. Con
// D = 1 es exactamente DetectorGLR (mismas operaciones). Tras una alarma espera D - 1 entradas: el siguiente elemento se forma solo con
// datos posteriores a la alarma (lo que hace el reinicio del disjunto). Una entrada no valida reinicia la memoria (sin imputar).
//
// Ciclo caliente: cero reservas (los anillos se reservan en la calibracion).
#pragma once

#include "ark/metrologia/canal_nivel.hpp"
#include "ark/metrologia/canal_temperatura.hpp"
#include "ark/metrologia/escalera_diadica.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"
#include "ark/ventanas/contratos_ventanas.hpp"   // RingBufferEspejoCxx

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

namespace ark::metrologia {

// ---------------------------------------------------------------------------------------------------------------------------------
// Escalera deslizante: las J escalas en cada paso nativo. Anillo del nivel k (0..J-2): las ultimas 2^k barras del nivel k.
template <std::size_t J_MAX>
class EscaleraDeslizante {
    static_assert(J_MAX >= 1 && J_MAX <= 40, "J_MAX entre 1 y 40");
    AgregacionEscala ag_ = AgregacionEscala::SUMA;
    std::size_t J_ = 0;
    std::vector<Acumulador> anillo_;
    std::array<std::size_t, J_MAX> off_{};
    std::uint64_t t_ = 0;

public:
    struct Salida {
        std::size_t niveles = 0;
        std::array<double, J_MAX> valor{};
        std::array<bool, J_MAX> valido{};
    };
    // Fuera del ciclo caliente: reserva sum_{k < J-1} 2^k acumuladores.
    void configurar(AgregacionEscala ag, std::size_t J) {
        ag_ = ag;
        J_ = std::min(J, J_MAX);
        std::size_t tot = 0;
        for (std::size_t k = 0; k + 1 < J_; ++k) { off_[k] = tot; tot += std::size_t{1} << k; }
        anillo_.assign(tot, Acumulador{});
        t_ = 0;
    }
    void reiniciar() noexcept { t_ = 0; }
    [[nodiscard]] std::size_t escalas() const noexcept { return J_; }

    Salida paso(double x, double den, bool valido) noexcept {
        Salida s{};
        s.niveles = J_;
        if (J_ == 0) return s;
        Acumulador cur = acumulador_nativo(ag_, x, den, valido);
        s.valor[0] = valor(ag_, cur);
        s.valido[0] = cur.valido && std::isfinite(s.valor[0]);
        for (std::size_t j = 1; j < J_; ++j) {
            const std::size_t mk = std::size_t{1} << (j - 1);
            Acumulador& slot = anillo_[off_[j - 1] + static_cast<std::size_t>(t_ % mk)];
            const Acumulador viejo = slot;   // b_{j-1}(t - 2^{j-1})
            slot = cur;                      // b_{j-1}(t)
            if (t_ + 1 >= (std::uint64_t{1} << j)) {   // b_j(t) usa las barras nativas t - 2^j + 1 .. t
                cur = combinar(ag_, viejo, cur);
                s.valor[j] = valor(ag_, cur);
                s.valido[j] = cur.valido && std::isfinite(s.valor[j]);
            } else {
                cur = Acumulador{};
                cur.valido = false;
                s.valor[j] = NaN;
                s.valido[j] = false;
            }
        }
        ++t_;
        return s;
    }
};

// ---------------------------------------------------------------------------------------------------------------------------------
template <std::size_t L>
class DetectorGLREspaciado {
    static_assert(L >= 2, "L >= 2");
    std::vector<double> z_;
    std::size_t D_ = 1, cap_ = 0, c_ = 0, espera_ = 0;
    std::uint64_t t_ = 0;
    double mu_ = NaN, s_ = NaN, tau_ = NaN, h_ = NaN;

    [[nodiscard]] double atras(std::size_t q) const noexcept { return z_[static_cast<std::size_t>((t_ - 1 - static_cast<std::uint64_t>(q) * D_) % cap_)]; }

public:
    using Paso = typename DetectorGLR<L>::Paso;
    // Fuera del ciclo caliente: reserva L D entradas.
    void configurar(double mu, double s, double tau, double h, std::size_t D) {
        mu_ = mu; s_ = s; tau_ = tau; h_ = h;
        D_ = D ? D : 1;
        cap_ = L * D_;
        z_.assign(cap_, 0.0);
        t_ = 0;
        reiniciar();
    }
    void reiniciar() noexcept { c_ = 0; espera_ = 0; }
    [[nodiscard]] std::size_t espaciado() const noexcept { return D_; }

    Paso paso(double y, bool valido = true) noexcept {
        Paso p{};
        // una entrada no valida reinicia todo: el siguiente elemento valido ya se forma con datos posteriores (sin imputar)
        if (!valido || !std::isfinite(y)) { c_ = 0; espera_ = 0; return p; }
        if (espera_ > 0) { --espera_; return p; }
        z_[static_cast<std::size_t>(t_ % cap_)] = (y - mu_) / s_;
        ++t_;
        ++c_;
        const std::size_t n = std::min(L, (c_ - 1) / D_ + 1);
        double S = 0.0, gmax = -1.0, smax = 0.0;
        std::size_t lmax = 0;
        for (std::size_t l = 1; l <= n; ++l) {
            S += atras(l - 1);
            const double g = S * S / (2.0 * static_cast<double>(l) * tau_);
            if (g > gmax) { gmax = g; lmax = l; smax = S; }
        }
        p.g = gmax; p.l = lmax; p.suma = smax;
        if (gmax >= h_) {
            p.alarma = true;
            if (lmax >= 3) {
                const double kb = 0.5 * (static_cast<double>(lmax) + 1.0), zb = smax / static_cast<double>(lmax);
                double sxx = 0.0, sxy = 0.0;
                for (std::size_t k = 1; k <= lmax; ++k) {   // k = 1 el mas antiguo de los l*
                    const double zk = atras(lmax - k);
                    sxx += (k - kb) * (k - kb);
                    sxy += (k - kb) * (zk - zb);
                }
                p.t_pendiente = (sxy / sxx) / std::sqrt(tau_ / sxx);
            }
            if (reiniciar_en_alarma_) {
                c_ = 0;
                espera_ = D_ - 1;
            }
        }
        return p;
    }
    // Sin reinicio automatico (temperatura por acumulacion): `alarma` = G >= h en cada entrada; el canal cierra con cerrar().
    void reiniciar_en_alarma(bool r) noexcept { reiniciar_en_alarma_ = r; }
    // Cierre explicito: como tras una alarma (el siguiente elemento, solo con datos posteriores).
    void cerrar() noexcept { c_ = 0; espera_ = D_ - 1; }

private:
    bool reiniciar_en_alarma_ = true;
};

// Plazo de confirmacion de la temperatura deslizante, en elementos del GLR (cada uno abarca T* barras de la escala): el tiempo de Lorden
// con el umbral CALIBRADO de confirmacion, t* = 2 h_c tau / delta^2, mas uno (enmienda 1, punto 2 del Comandante). Lo usan el canal y el
// observador (una sola funcion). Sin h_c finito o sin amplitud: L elementos (la memoria del GLR, declarado).
inline std::int64_t plazo_elementos_temperatura(double h_c, double tau, double amplitud, double s0, std::size_t L) noexcept {
    const double d = std::abs(amplitud) / s0;
    const double ts = 2.0 * h_c * tau / (d * d);
    if (!std::isfinite(ts) || !(ts >= 0.0)) return static_cast<std::int64_t>(L);
    return static_cast<std::int64_t>(std::ceil(std::min(ts, 1e15))) + 1;
}

namespace detalle_desl {

// ARL del detector espaciado sobre una serie (NaN = no valida): pasos validos / alarmas.
template <std::size_t L>
inline double arl(const std::vector<double>& serie, double mu, double s, double tau, double h, std::size_t D, std::size_t& al) {
    DetectorGLREspaciado<L> d;
    d.configurar(mu, s, tau, h, D);
    std::size_t total = 0;
    al = 0;
    for (double v : serie) {
        const bool ok = std::isfinite(v);
        total += ok;
        al += d.paso(v, ok).alarma;
    }
    return al ? static_cast<double>(total) / static_cast<double>(al) : std::numeric_limits<double>::infinity();
}

struct Tramo { std::size_t ini = 0, n = 0; };
inline std::vector<Tramo> tramos(const std::vector<double>& serie) {
    std::vector<Tramo> out;
    for (std::size_t i = 0; i < serie.size();) {
        if (!std::isfinite(serie[i])) { ++i; continue; }
        std::size_t k = i;
        while (k < serie.size() && std::isfinite(serie[k])) ++k;
        out.push_back({i, k - i});
        i = k;
    }
    return out;
}

// La ventana de T barras (de la mas antigua a la presente) con el MISMO Welford del canal de temperatura disjunto.
inline double s2_welford(const double* y, std::size_t T) noexcept {
    double media = 0.0, m2 = 0.0;
    for (std::size_t i = 0; i < T; ++i) {
        const double d = y[i] - media;
        media += d / static_cast<double>(i + 1);
        m2 += d * (y[i] - media);
    }
    return m2 / static_cast<double>(T - 1);
}

}  // namespace detalle_desl

// ---------------------------------------------------------------------------------------------------------------------------------
// Canal de nivel deslizante de la escala j: GLR espaciado D = 2^j sobre la barra deslizante, evaluado en cada paso nativo.
template <std::size_t L = 120>
class CanalNivelDeslizante {
    FichaEscala f_{};
    DetectorGLREspaciado<L> glr_{};

public:
    [[nodiscard]] const FichaEscala& ficha() const noexcept { return f_; }
    void reiniciar_flujo() noexcept {
        if (f_.detectable) glr_.configurar(f_.nulo.mu0, f_.nulo.s0, f_.nulo.tau_int, f_.det.h, static_cast<std::size_t>(f_.m));
    }

    // base: la ficha del canal disjunto de la escala (nulo y lambda_TV). serie: la barra deslizante de la escala en cada paso nativo de
    // la historia (NaN = no valida). ag: la agregacion (para la inyeccion del tau50). Fuera del ciclo caliente.
    const FichaEscala& calibrar(const FichaEscala& base, const std::vector<double>& serie, AgregacionEscala ag, const ConfiguracionNivel& cfg) {
        FichaEscala f = base;
        if (!base.detectable) { f_ = f; return f_; }
        const std::size_t D = static_cast<std::size_t>(base.m);
        const double arl0_nativo = cfg.arl0_pasos * static_cast<double>(base.m);   // cada canal: ARL0 2J pasos nativos
        auto eval = [&](double h, std::size_t& al) { return detalle_desl::arl<L>(serie, f.nulo.mu0, f.nulo.s0, f.nulo.tau_int, h, D, al); };
        const auto rh = detalle_nivel::buscar_h_con(eval, arl0_nativo, cfg.alarmas_minimas);
        f.det.h = rh.h; f.det.h_extrapolado = rh.extrapolado; f.det.alarmas_nulo = rh.alarmas;
        for (std::size_t i = 0; i < DELTA_REPORTE.size(); ++i) { f.det.tau50[i] = NaN; f.det.fraccion_sin_detectar[i] = NaN; }
        if (!std::isfinite(rh.h)) { f.detectable = false; f.motivo = "h deslizante no finito: el nulo no da alarmas legibles"; f_ = f; return f_; }
        f.motivo = "calibrada (deslizante)";
        if (cfg.medir_tau50) medir_tau50(f, serie, ag, cfg);
        f_ = f;
        reiniciar_flujo();
        return f_;
    }

    // Ciclo caliente: la barra deslizante de la escala en el paso nativo t.
    SalidaNivel paso(double y, bool valido, std::int64_t t_nativo) noexcept {
        SalidaNivel s{};
        s.emitida = true;
        if (!f_.detectable) return s;
        const auto p = glr_.paso(y, valido);
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
        return s;
    }

private:
    // tau50 deslizante: escalon delta*s0 de la barra, inyectado en una FASE al azar de la historia nula. La barra deslizante lo recibe
    // en rampa: fraccion min(1, (i - t0 + 1)/2^j) (exacta para SUMA y MEDIA, linealizada para las demas -declarado-; ULTIMO entero).
    // Retardo en pasos de la escala = pasos nativos / 2^j. Mismo procedimiento que el disjunto (con 2^j = 1 coincide bit a bit).
    void medir_tau50(FichaEscala& f, const std::vector<double>& serie, AgregacionEscala ag, const ConfiguracionNivel& cfg) const {
        const std::size_t m = static_cast<std::size_t>(f.m), W0 = L * m;
        const auto tr = detalle_desl::tramos(serie);
        std::vector<std::size_t> largos;
        for (std::size_t i = 0; i < tr.size(); ++i) if (tr[i].n >= 2 * W0) largos.push_back(i);
        if (largos.empty()) return;
        std::mt19937_64 rng(cfg.semilla ^ static_cast<std::uint64_t>(f.j));
        const bool entero = m == 1 || ag == AgregacionEscala::ULTIMO;
        for (std::size_t k = 0; k < DELTA_REPORTE.size(); ++k) {
            const double sh = DELTA_REPORTE[k] * f.nulo.s0;
            std::vector<double> taus;
            for (std::size_t r = 0; r < cfg.replicas_tau50; ++r) {
                const auto& g = tr[largos[std::uniform_int_distribution<std::size_t>(0, largos.size() - 1)(rng)]];
                const std::size_t t0 = std::uniform_int_distribution<std::size_t>(W0, g.n - W0)(rng);
                const std::size_t fin = std::min(g.n, t0 + static_cast<std::size_t>(std::ceil(20.0 * f.det.t_lorden[k])) * m + 1);
                DetectorGLREspaciado<L> d;
                d.configurar(f.nulo.mu0, f.nulo.s0, f.nulo.tau_int, f.det.h, m);
                double tau = std::numeric_limits<double>::infinity();
                for (std::size_t i = t0 - W0; i < fin; ++i) {
                    const double fr = entero ? 1.0 : std::min(1.0, static_cast<double>(i + 1 - t0) / static_cast<double>(m));
                    const auto p = d.paso(serie[g.ini + i] + (i >= t0 ? sh * fr : 0.0), true);
                    if (p.alarma && i >= t0 && p.suma > 0.0) { tau = static_cast<double>(i - t0 + 1) / static_cast<double>(m); break; }
                }
                taus.push_back(tau);
            }
            std::size_t nd = 0;
            for (double t : taus) nd += !std::isfinite(t);
            f.det.fraccion_sin_detectar[k] = static_cast<double>(nd) / static_cast<double>(taus.size());
            std::sort(taus.begin(), taus.end());
            const std::size_t n = taus.size();
            f.det.tau50[k] = (n % 2) ? taus[n / 2] : 0.5 * (taus[n / 2 - 1] + taus[n / 2]);
        }
    }
};

// ---------------------------------------------------------------------------------------------------------------------------------
// Canal de temperatura deslizante de la escala j: la ventana de T* barras avanza una barra; GLR espaciado D = T* sobre log beta^.
template <std::size_t L = 120>
class CanalTemperaturaDeslizante {
    static constexpr std::size_t T_MAX = REJILLA_T.back();
    FichaEscala f_{};
    DetectorGLREspaciado<L> glr_{};
    ark::ventanas::RingBufferEspejoCxx<double, T_MAX> ventana_{};
    std::size_t T_ = 0, racha_ = 0;
    // evento abierto por acumulacion (enmienda 1): signo, barras desde la apertura y su limite (el plazo)
    bool abierto_ = false;
    int signo_ = 0;
    std::int64_t barras_ = 0, limite_ = 0;

    void llenar(EventoAlarma& a, const typename DetectorGLREspaciado<L>::Paso& p, std::int64_t t_nativo) const noexcept {
        a.hay = true;
        a.signo = p.suma > 0.0 ? 1 : -1;     // + : beta sube (se enfria), - : beta baja (se calienta)
        a.t = t_nativo;
        a.l = p.l;
        a.t_cambio = t_nativo - static_cast<std::int64_t>(p.l) * f_.m + 1;
        a.amplitud = p.suma / static_cast<double>(p.l) * f_.nulo.s0;
        a.se_amplitud = f_.nulo.s0 * std::sqrt(f_.nulo.tau_int / static_cast<double>(p.l));
        a.g = p.g;
        a.t_pendiente = p.t_pendiente;
    }

public:
    [[nodiscard]] const FichaEscala& ficha_glr() const noexcept { return f_; }
    [[nodiscard]] std::size_t T() const noexcept { return T_; }
    [[nodiscard]] bool evento_abierto() const noexcept { return abierto_; }
    void reiniciar_flujo() noexcept {
        racha_ = 0;
        abierto_ = false;
        if (f_.detectable) {
            glr_.configurar(f_.nulo.mu0, f_.nulo.s0, f_.nulo.tau_int, f_.det.h, T_);
            glr_.reiniciar_en_alarma(false);
        }
    }

    // base: el canal disjunto ya calibrado (T* y el nulo de log beta^ en ventanas disjuntas). segs: la historia de barras de la
    // escala. arl0_barras: ARL0 del canal en barras de la escala. Fuera del ciclo caliente.
    const FichaEscala& calibrar(const CanalTemperatura<L>& base, const std::vector<std::vector<double>>& segs, double arl0_barras,
                                const ConfiguracionTemperatura& cfg) {
        f_ = base.ficha_glr();
        T_ = base.T();
        racha_ = 0;
        if (T_ == 0 || !f_.detectable) return f_;
        std::vector<double> lb;
        for (const auto& g : segs) {
            for (std::size_t i = 0; i + T_ <= g.size(); ++i) {
                const double s2 = detalle_desl::s2_welford(g.data() + i, T_);
                lb.push_back(s2 > 0.0 ? -std::log(s2) : NaN);
            }
            lb.push_back(NaN);   // los segmentos no se tocan
        }
        auto eval = [&](double h, std::size_t& al) { return detalle_desl::arl<L>(lb, f_.nulo.mu0, f_.nulo.s0, f_.nulo.tau_int, h, T_, al); };
        const auto rh = detalle_nivel::buscar_h_con(eval, arl0_barras, cfg.alarmas_minimas);
        f_.det.arl0_pasos = arl0_barras;
        f_.det.h = rh.h; f_.det.h_extrapolado = rh.extrapolado; f_.det.alarmas_nulo = rh.alarmas;
        for (std::size_t i = 0; i < DELTA_REPORTE.size(); ++i)
            f_.det.t_lorden[i] = std::log(arl0_barras) * 2.0 * f_.nulo.tau_int / (DELTA_REPORTE[i] * DELTA_REPORTE[i]);
        if (!std::isfinite(rh.h)) { f_.detectable = false; f_.motivo = "h deslizante no finito: el nulo no da alarmas legibles"; return f_; }
        // umbral de confirmacion por acumulacion: ARL(h_c) = ARL0 / alfa sobre el mismo nulo y la misma construccion (con su reinicio);
        // en tasa, una de cada 1/alfa alarmas falsas llega a confirmarse. Sin h_c finito no se confirma nada (TRANSITORIO, declarado).
        const auto rc = detalle_nivel::buscar_h_con(eval, arl0_barras / cfg.alfa, cfg.alarmas_minimas);
        f_.det.h_confirmacion = rc.h;
        f_.motivo = "calibrada (deslizante, confirmacion por acumulacion)";
        reiniciar_flujo();
        return f_;
    }

    // Ciclo caliente: una barra de la escala (la disjunta). emitida = hay ventana completa de T* barras validas.
    // Confirmacion por acumulacion (enmienda 1): la apertura (G >= h, sin evento abierto) NO reinicia el GLR; desde la barra siguiente,
    // G >= h_c con el mismo signo emite la alarma de confirmacion (s.confirmacion = true) y cierra; si vence el plazo, cierra sin
    // emitir (el observador lo da por TRANSITORIO en el mismo paso). Cerrar = reinicio y espera de T* - 1 barras.
    SalidaTemperatura paso(double y, bool valido, std::int64_t t_nativo) noexcept {
        SalidaTemperatura s{};
        if (T_ == 0) return s;
        if (!valido || !std::isfinite(y)) { racha_ = 0; abierto_ = false; glr_.reiniciar(); return s; }
        ventana_.push(y);
        if (racha_ < T_) ++racha_;
        if (racha_ < T_) return s;
        const double s2 = detalle_desl::s2_welford(ventana_.ventana_actual() + (T_MAX - T_), T_);
        s.emitida = true;
        if (!(s2 > 0.0)) { abierto_ = false; glr_.reiniciar(); return s; }   // ventana constante: no se estima (declarado)
        s.log_beta = -std::log(s2);
        if (!f_.detectable) return s;
        const auto p = glr_.paso(s.log_beta, true);
        s.g = p.g;
        if (!abierto_) {
            if (p.alarma) {
                llenar(s.alarma, p, t_nativo);
                abierto_ = true;
                signo_ = s.alarma.signo;
                barras_ = 0;
                limite_ = static_cast<std::int64_t>(T_) *
                          plazo_elementos_temperatura(f_.det.h_confirmacion, f_.nulo.tau_int, s.alarma.amplitud, f_.nulo.s0, L);
            }
            return s;
        }
        ++barras_;
        if (p.g >= f_.det.h_confirmacion && (p.suma > 0.0 ? 1 : -1) == signo_) {
            llenar(s.alarma, p, t_nativo);
            s.confirmacion = true;
            abierto_ = false;
            glr_.cerrar();
        } else if (barras_ >= limite_) {
            abierto_ = false;
            glr_.cerrar();
        }
        return s;
    }
};

}  // namespace ark::metrologia
