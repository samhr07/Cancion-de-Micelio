// flujo_barras.hpp -- fuente de FILAS comun a los operadores y flujo de barras intrinsecas una a una (ARK-WIN-01-CXX), C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04, por orden del Comandante ("conectelos"): el operador fractal y el
// ObservadorSaltos sobre el observable crudo reciben filas por la misma aduana y el mismo reloj intrinseco que el sensor de memoria.
//
// FuenteFilasCxx declara, una sola vez para cualquier operador: la configuracion del motor de pared que hace de aduana (reloj PARED,
// NU_TRADES, NU_VOLUMEN, NU_NOCIONAL o NU_INFO_SHANNON, tick, presupuesto, sesiones), el ancho de barra delta_nu (NaN = nu medio
// por fila con sesion y sin imputar de la historia de calibracion), la agregacion de la barra y la FASE de la purga diurna (desde el
// calendario del activo; inactiva = la entrega el llamador).
//
// FlujoBarrasIntrinsecasCxx: el motor intrinseco con W = 1 y S = 1 (cada barra que se cierra se entrega sola). Lo usan los
// operadores con ventana propia (DFA del fractal, escalera del ObservadorSaltos); el sensor de memoria usa el motor con su W.
// Una barra con alguna fila imputada se entrega marcada (`imputada`): el operador decide que hacer (los operadores de la base no
// estiman sobre ella; ver los README).
#pragma once

#include "ark/dataset/fase_calendario.hpp"
#include "ark/ventanas/motor_reloj_intrinseco.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <vector>

namespace ark::ventanas {

template <size_t W_MAX>
struct FuenteFilasCxx {
    ConfiguracionVentanaCausalCxx<W_MAX> base{};                     // aduana y reloj
    double delta_nu = std::numeric_limits<double>::quiet_NaN();      // ancho de barra en nu; NaN = delta_nu_por_fila de la historia
    AgregacionBarra agregacion = AgregacionBarra::ULTIMO;            // ULTIMO niveles, SUMA flujos, MEDIA_NU intensivos
    ark::dataset::FaseCalendarioCxx fase{};                          // activa: fase y dia desde el calendario y la hora local
};

// Una barra intrinseca ya cerrada.
struct BarraIntrinsecaCxx {
    double valor = std::numeric_limits<double>::quiet_NaN();
    bool imputada = false;          // contiene alguna fila con sesion imputada
    double nu = 0.0;                // nu de sus filas
    uint32_t filas = 0;             // filas con sesion
    uint32_t fronteras = 1;         // fronteras k delta_nu que cruzo (> 1: fronteras saltadas por una fila grande)
    uint64_t t_cierre_ms = 0;       // marca de la ultima fila
    uint64_t t_disponible_ms = 0;   // maxima disponibilidad de sus filas (B5)
    size_t fase = 0;                // de la purga (calendario si la fuente lo activa; si no, 0)
    int64_t dia = 0;
};

// delta_nu efectivo de una fuente: el declarado o el nu medio por fila de la historia.
template <size_t W_MAX>
[[nodiscard]] inline double resolver_delta_nu(const FuenteFilasCxx<W_MAX>& f, std::span<const EntradaPasoCxx> filas) noexcept {
    if (std::isfinite(f.delta_nu)) return f.delta_nu;
    return delta_nu_por_fila<W_MAX>(f.base, filas.data(), filas.size());
}

template <size_t W_MAX>
class FlujoBarrasIntrinsecasCxx {
    FuenteFilasCxx<W_MAX> fuente_{};
    MotorVentanasIntrinsecoCxx<W_MAX> motor_{};
    VentanaIntrinsecaCxx<W_MAX> v_{};
    ErrorMotorIntrinseco error_{ErrorMotorIntrinseco::BASE_INVALIDA};
    ResultadoPaso ultimo_{ResultadoPaso::SIN_EMISION};

public:
    // Configura con un delta_nu ya resuelto (resolver_delta_nu). La W de la aduana no limita: el flujo emite cada barra (W = 1).
    ErrorMotorIntrinseco configurar(const FuenteFilasCxx<W_MAX>& f, double delta_nu) noexcept {
        fuente_ = f;
        ConfiguracionMotorIntrinsecoCxx<W_MAX> c{};
        c.base = f.base;
        c.base.min_muestras_requeridas = 1;   // la ventana propia la lleva el operador; aqui se entrega barra a barra
        c.W = 1;
        c.S = 1;
        c.delta_nu = delta_nu;
        c.agregacion = f.agregacion;
        error_ = motor_.configurar(c);
        return error_;
    }
    [[nodiscard]] ErrorMotorIntrinseco error_configuracion() const noexcept { return error_; }
    [[nodiscard]] const FuenteFilasCxx<W_MAX>& fuente() const noexcept { return fuente_; }
    [[nodiscard]] const MotorVentanasIntrinsecoCxx<W_MAX>& motor() const noexcept { return motor_; }
    [[nodiscard]] ResultadoPaso ultimo_resultado() const noexcept { return ultimo_; }   // error de la aduana si la fila se rechazo
    void reset() noexcept { motor_.reset(); ultimo_ = ResultadoPaso::SIN_EMISION; }

    // Una fila. Devuelve la barra si en esta fila se cerro una. Cero reservas.
    std::optional<BarraIntrinsecaCxx> paso(const EntradaPasoCxx& e) noexcept {
        if (error_ != ErrorMotorIntrinseco::NINGUNO) { ultimo_ = ResultadoPaso::ERROR_CONFIGURACION_INVALIDA; return std::nullopt; }
        const PasoIntrinsecoCxx p = motor_.paso(e, v_);
        ultimo_ = p.resultado;
        if (p.resultado != ResultadoPaso::VENTANA_EMITIDA) return std::nullopt;
        BarraIntrinsecaCxx b{};
        b.valor = v_.valores[0];
        b.imputada = v_.mascara[0];
        b.nu = v_.nu_barra[0];
        b.filas = v_.filas_barra[0];
        b.fronteras = v_.fronteras_barra[0];
        b.t_cierre_ms = v_.t_cierre_ms;
        b.t_disponible_ms = v_.t_actual;
        if (fuente_.fase.activa) {
            const auto fd = ark::dataset::fase_de(fuente_.fase, static_cast<int64_t>(b.t_cierre_ms));
            b.fase = fd.fase;
            b.dia = fd.dia;
        }
        return b;
    }
};

// Barras de una historia de filas (calibracion: reserva). Mismo motor y misma aduana que el ciclo caliente.
template <size_t W_MAX>
[[nodiscard]] inline std::vector<BarraIntrinsecaCxx> barras_de_historia(const FuenteFilasCxx<W_MAX>& f, double delta_nu,
                                                                         std::span<const EntradaPasoCxx> filas) {
    std::vector<BarraIntrinsecaCxx> out;
    auto fl = std::make_unique<FlujoBarrasIntrinsecasCxx<W_MAX>>();
    if (fl->configurar(f, delta_nu) != ErrorMotorIntrinseco::NINGUNO) return out;
    for (const auto& e : filas)
        if (auto b = fl->paso(e)) out.push_back(*b);
    return out;
}

}  // namespace ark::ventanas
