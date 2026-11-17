// escalera_diadica.hpp -- escalera causal de escalas 2^j y purga diurna causal (ARK-METRO-01). Zero-heap.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04.
//
// La escalera trabaja sobre la REPRESENTACION ESTACIONARIA del observable (un nivel ya llega convertido en incrementos nativos y se
// agrega por SUMA). El nivel j emite una barra de 2^j barras nativas cuando se cierra su par (arbol binario): nunca emite una barra
// abierta y una barra solo es valida si sus dos hijas lo son (sin imputar).
#pragma once

#include "ark/metrologia/tipos_metrologia.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace ark::metrologia {

// Acumulador de una barra. Su significado depende de la agregacion (ver valor()).
struct Acumulador {
    double a = 0.0, b = 0.0;
    std::int64_t n = 0;
    bool valido = true;
};

inline Acumulador acumulador_nativo(AgregacionEscala ag, double x, double den, bool valido) noexcept {
    Acumulador c{};
    c.n = 1;
    c.valido = valido && std::isfinite(x) && (ag != AgregacionEscala::RAZON_DE_SUMAS || std::isfinite(den));
    switch (ag) {
        case AgregacionEscala::RAIZ_SUMA_CUADRADOS: c.a = x * x; break;
        case AgregacionEscala::RAZON_DE_SUMAS: c.a = x; c.b = den; break;   // x = numerador, den = denominador
        default: c.a = x; break;
    }
    return c;
}

inline Acumulador combinar(AgregacionEscala ag, const Acumulador& p, const Acumulador& q) noexcept {
    Acumulador c{};
    c.n = p.n + q.n;
    c.valido = p.valido && q.valido;
    if (ag == AgregacionEscala::ULTIMO) { c.a = q.a; return c; }
    c.a = p.a + q.a;
    c.b = p.b + q.b;
    return c;
}

inline double valor(AgregacionEscala ag, const Acumulador& c) noexcept {
    switch (ag) {
        case AgregacionEscala::MEDIA: return c.a / static_cast<double>(c.n);
        case AgregacionEscala::RAZON_DE_SUMAS: return c.b != 0.0 ? c.a / c.b : NaN;
        case AgregacionEscala::LOG_DE_SUMA: return c.a > 0.0 ? std::log(c.a) : NaN;
        case AgregacionEscala::RAIZ_SUMA_CUADRADOS: return std::sqrt(c.a);
        default: return c.a;   // ULTIMO, SUMA
    }
}

template <std::size_t J_MAX>
class EscaleraDiadica {
    static_assert(J_MAX >= 1 && J_MAX <= 40, "J_MAX entre 1 y 40");
    AgregacionEscala ag_ = AgregacionEscala::SUMA;
    std::array<Acumulador, J_MAX> pend_{};    // pend_[j]: barra del nivel j-1 que espera su pareja (j >= 1)
    std::array<bool, J_MAX> hay_{};

public:
    struct Salida {
        std::size_t niveles = 0;                    // emiten los niveles 0..niveles-1 en este paso
        std::array<double, J_MAX> valor{};
        std::array<bool, J_MAX> valido{};
    };
    EscaleraDiadica() = default;
    explicit EscaleraDiadica(AgregacionEscala ag) : ag_(ag) {}
    void reiniciar() noexcept { hay_.fill(false); }

    Salida paso(double x, double den, bool valido) noexcept {
        Salida s{};
        Acumulador cur = acumulador_nativo(ag_, x, den, valido);
        s.valor[0] = valor(ag_, cur);
        s.valido[0] = cur.valido && std::isfinite(s.valor[0]);
        s.niveles = 1;
        for (std::size_t j = 1; j < J_MAX; ++j) {
            if (!hay_[j]) { pend_[j] = cur; hay_[j] = true; break; }
            cur = combinar(ag_, pend_[j], cur);
            hay_[j] = false;
            s.valor[j] = valor(ag_, cur);
            s.valido[j] = cur.valido && std::isfinite(s.valor[j]);
            s.niveles = j + 1;
        }
        return s;
    }
};

// Purga exogena: perfil diurno CAUSAL de la varianza. y <- y sqrt(v_media / v_fase), con v_fase la media de y^2 de esa fase en los
// dias ANTERIORES (lo del dia en curso se acumula aparte y se incorpora al cambiar de dia). Una fase se activa cuando tiene dato de
// dias_minimos dias distintos (ceil(1/alfa) = 20, declarado); antes su salida es NO valida. La fase la codifica quien llama: p. ej.
// fase = tipo_dia * horas + hora separa el perfil laborable del de fin de semana (cada tipo madura con sus propios dias).
template <std::size_t B_MAX>
class PurgaDiurna {
    std::size_t bins_ = 1;
    std::size_t dias_minimos_ = 20;
    std::array<double, B_MAX> s_{}, sh_{};      // sumas de y^2: comprometidas (dias previos) / del dia en curso
    std::array<std::int64_t, B_MAX> n_{}, nh_{}, dias_b_{};
    std::int64_t dia_ = -1;
    std::size_t dias_ = 0;
    double vm_ = NaN;                            // media entre fases de v_fase (se actualiza al comprometer un dia)

    void comprometer() noexcept {
        double vm = 0.0;
        std::size_t nb = 0;
        for (std::size_t b = 0; b < bins_; ++b) {
            if (nh_[b] > 0) ++dias_b_[b];
            s_[b] += sh_[b]; n_[b] += nh_[b]; sh_[b] = 0.0; nh_[b] = 0;
            if (n_[b] > 0) { vm += s_[b] / static_cast<double>(n_[b]); ++nb; }
        }
        vm_ = nb ? vm / static_cast<double>(nb) : NaN;
    }

public:
    PurgaDiurna() = default;
    PurgaDiurna(std::size_t bins, std::size_t dias_minimos) : bins_(bins < B_MAX ? (bins ? bins : 1) : B_MAX), dias_minimos_(dias_minimos) {}
    [[nodiscard]] bool activa(std::size_t b) const noexcept { return dias_b_[b % bins_] >= static_cast<std::int64_t>(dias_minimos_); }

    // y: valor de la fase `fase` (0..bins-1) del dia `dia` (entero monotono). Devuelve el valor purgado; `valido` = false si aun no.
    double paso(double y, std::size_t fase, std::int64_t dia, bool& valido) noexcept {
        if (dia != dia_) { if (dia_ >= 0) { comprometer(); ++dias_; } dia_ = dia; }
        const std::size_t b = fase % bins_;
        valido = false;
        double out = NaN;
        if (activa(b) && n_[b] > 0) {
            const double vb = s_[b] / static_cast<double>(n_[b]);
            if (vb > 0.0 && std::isfinite(y)) { out = y * std::sqrt(vm_ / vb); valido = true; }
        }
        if (std::isfinite(y)) { sh_[b] += y * y; ++nh_[b]; }
        return out;
    }
};

}  // namespace ark::metrologia
