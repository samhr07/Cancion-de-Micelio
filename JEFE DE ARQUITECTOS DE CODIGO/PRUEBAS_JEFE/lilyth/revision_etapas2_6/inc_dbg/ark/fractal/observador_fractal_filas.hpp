// observador_fractal_filas.hpp -- el operador fractal (ARK-FRAC-01) con entrada por FILAS en reloj intrinseco, C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04, por orden del Comandante ("conectelos").
//
// filas -> aduana del motor de pared auditado -> barra intrinseca (FlujoBarrasIntrinsecasCxx, delta_nu de la fuente) ->
// ObservadorFractalUniversal (sensor DFA de W barras, nulo propio, CUSUM, dualidad L2/L1, tau50). El observador no se toca: este
// adaptador solo cambia de donde vienen sus muestras.
//
// Imputadas (declarado):
//   ciclo caliente  una barra con alguna fila imputada entra como NaN: las W ventanas del DFA que la contienen salen no validas y
//                   el regularizador corta su ventana (la misma regla que el sensor de memoria: no se estima sobre datos inventados);
//   calibracion     la historia entra con el valor que la aduana dio a esas barras (arrastre del ultimo valido en ULTIMO; 0 de
//                   aportacion en SUMA): los mundos del nulo reconstruyen el nivel por suma de incrementos y un NaN lo anularia
//                   desde ese punto. Se reporta cuantas barras imputadas tenia la historia.
// Ciclo caliente: cero reservas. La calibracion reserva.
#pragma once

#include "ark/fractal/observador_fractal_multiescala.hpp"
#include "ark/ventanas/flujo_barras.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

namespace ark::fractal {

template <std::size_t W = 50, std::size_t S = 2, std::size_t N_REG = 120, std::size_t W_FLUJO = 256>
class ObservadorFractalFilas {
public:
    using Fuente = ark::ventanas::FuenteFilasCxx<W_FLUJO>;
    using Observador = ObservadorFractalUniversal<double, W, S, N_REG>;
    struct Salida {
        bool barra_cerrada = false;
        ark::ventanas::BarraIntrinsecaCxx barra{};
        std::optional<EstadoFractal> estado{};   // cuando el sensor evalua en esta barra
    };

private:
    Observador obs_{};
    ark::ventanas::FlujoBarrasIntrinsecasCxx<W_FLUJO> flujo_{};
    double delta_nu_ = NaN;
    std::size_t barras_calibracion_ = 0, imputadas_calibracion_ = 0;
    bool lista_ = false;

public:
    [[nodiscard]] bool lista() const noexcept { return lista_; }
    [[nodiscard]] double delta_nu() const noexcept { return delta_nu_; }
    [[nodiscard]] const Observador& observador() const noexcept { return obs_; }
    [[nodiscard]] const CalibracionFractal& calibracion() const noexcept { return obs_.calibracion(); }
    [[nodiscard]] const ark::ventanas::FlujoBarrasIntrinsecasCxx<W_FLUJO>& flujo() const noexcept { return flujo_; }
    [[nodiscard]] std::size_t barras_calibracion() const noexcept { return barras_calibracion_; }
    [[nodiscard]] std::size_t imputadas_calibracion() const noexcept { return imputadas_calibracion_; }

    // Serie de barras de una historia (la misma que ve la calibracion).
    [[nodiscard]] std::vector<ark::ventanas::BarraIntrinsecaCxx> barras(const Fuente& f, double delta_nu,
                                                                        std::span<const ark::ventanas::EntradaPasoCxx> filas) const {
        return ark::ventanas::barras_de_historia<W_FLUJO>(f, delta_nu, filas);
    }

    const CalibracionFractal& calibrar_filas(const ConfiguracionFractal& cfg, const Fuente& f, std::span<const ark::ventanas::EntradaPasoCxx> filas) {
        lista_ = false;
        delta_nu_ = ark::ventanas::resolver_delta_nu<W_FLUJO>(f, filas);
        if (!(delta_nu_ > 0.0)) throw std::invalid_argument("delta_nu: no se pudo derivar de la historia (sin filas con sesion)");
        const auto b = barras(f, delta_nu_, filas);
        std::vector<double> x(b.size());
        imputadas_calibracion_ = 0;
        for (std::size_t i = 0; i < b.size(); ++i) { x[i] = b[i].valor; imputadas_calibracion_ += b[i].imputada; }
        barras_calibracion_ = b.size();
        obs_.configurar(cfg);
        obs_.calibrar(std::span<const double>(x));
        obs_.reiniciar_flujo();
        if (flujo_.configurar(f, delta_nu_) != ark::ventanas::ErrorMotorIntrinseco::NINGUNO) return obs_.calibracion();
        lista_ = obs_.calibracion().lista;
        return obs_.calibracion();
    }

    // Ciclo caliente: una fila. Cero reservas.
    Salida paso_fila(const ark::ventanas::EntradaPasoCxx& e) noexcept {
        Salida s{};
        auto b = flujo_.paso(e);
        if (!b) return s;
        s.barra_cerrada = true;
        s.barra = *b;
        s.estado = obs_.paso(b->imputada ? NaN : b->valor);
        return s;
    }
};

}  // namespace ark::fractal
