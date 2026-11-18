// observador_saltos_filas.hpp -- el bloque de saltos y temperatura (ARK-METRO-01) sobre el observable CRUDO con entrada por FILAS en
// reloj intrinseco, C++20. Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04, por orden del Comandante ("conectelos").
//
// filas -> aduana del motor de pared auditado -> barra intrinseca (FlujoBarrasIntrinsecasCxx) -> ObservadorSaltos::paso(valor de la
// barra, den, valido, fase, dia). El observador no se toca (su escalera 2^j agrega barras intrinsecas, no filas de pared).
//   valido = la barra no contiene filas imputadas (la escalera ya trata una barra no valida como hueco: un nivel solo es valido si
//            sus dos hijas lo son);
//   den    = nu de la barra (declarado): con agregacion RAZON_DE_SUMAS la escalera da sum(valor)/sum(nu), el observable por unidad de
//            reloj intrinseco (por ejemplo OFI por unidad de volumen). Con las otras agregaciones den no se usa;
//   fase y dia = los del calendario del activo en el cierre de la barra (fuente.fase activa; bins_fase se fija a 2 bins_dia), o los
//            que entregue el llamador para la fila que cierra la barra si la fuente no los activa.
// Ciclo caliente: cero reservas. La calibracion reserva.
#pragma once

#include "ark/metrologia/observador_saltos.hpp"
#include "ark/ventanas/flujo_barras.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

namespace ark::metrologia {

template <std::size_t J_MAX = 24, std::size_t L = 120, std::size_t N_REG = 120, std::size_t B_MAX = 1440, std::size_t EV_MAX = 8,
          std::size_t W_FLUJO = 256>
class ObservadorSaltosFilas {
public:
    using Fuente = ark::ventanas::FuenteFilasCxx<W_FLUJO>;
    using Delta = ObservadorSaltos<J_MAX, L, N_REG, B_MAX, EV_MAX>;
    struct Salida {
        bool barra_cerrada = false;
        ark::ventanas::BarraIntrinsecaCxx barra{};
        typename Delta::Salida delta{};   // solo si barra_cerrada
    };

private:
    Delta obs_{};
    ark::ventanas::FlujoBarrasIntrinsecasCxx<W_FLUJO> flujo_{};
    ConfiguracionObservador cfg_{};
    double delta_nu_ = NaN;
    std::size_t barras_calibracion_ = 0;
    bool lista_ = false;

public:
    [[nodiscard]] bool lista() const noexcept { return lista_; }
    [[nodiscard]] double delta_nu() const noexcept { return delta_nu_; }
    [[nodiscard]] const Delta& observador() const noexcept { return obs_; }
    [[nodiscard]] const ConfiguracionObservador& configuracion() const noexcept { return cfg_; }
    [[nodiscard]] std::size_t barras_calibracion() const noexcept { return barras_calibracion_; }
    // Otro ARL0 de la familia sobre la misma historia (ObservadorSaltos::recalibrar_arl0); lo usa el observador conjunto.
    void recalibrar_arl0(double arl0) { cfg_.arl0_pasos_nativos = arl0; obs_.recalibrar_arl0(arl0); }

    // fase y dia por FILA (solo si la fuente no activa el calendario): se toman los de la fila que cierra cada barra.
    void calibrar_filas(const ConfiguracionObservador& cfg, const Fuente& f, std::span<const ark::ventanas::EntradaPasoCxx> filas,
                        std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        lista_ = false;
        delta_nu_ = ark::ventanas::resolver_delta_nu<W_FLUJO>(f, filas);
        if (!(delta_nu_ > 0.0)) throw std::invalid_argument("delta_nu: no se pudo derivar de la historia (sin filas con sesion)");
        cfg_ = cfg;
        if (f.fase.activa) cfg_.bins_fase = f.fase.bins_fase();
        std::vector<double> x, den;
        std::vector<std::uint8_t> ok;
        std::vector<std::size_t> fs;
        std::vector<std::int64_t> ds;
        {
            auto fl = std::make_unique<ark::ventanas::FlujoBarrasIntrinsecasCxx<W_FLUJO>>();
            if (fl->configurar(f, delta_nu_) != ark::ventanas::ErrorMotorIntrinseco::NINGUNO) return;
            for (std::size_t i = 0; i < filas.size(); ++i) {
                auto b = fl->paso(filas[i]);
                if (!b) continue;
                x.push_back(b->valor);
                den.push_back(b->nu);
                const bool ok_b = cfg.usar_protocolo_huecos
                    ? ((!b->imputada || b->solo_mercado_quieto) && ((b->flags_hueco & (2u | 64u | 128u)) == 0u))
                    : (!b->imputada);
                ok.push_back(ok_b ? 1 : 0);
                fs.push_back(f.fase.activa ? b->fase : (fase.empty() ? 0 : fase[i]));
                ds.push_back(f.fase.activa ? b->dia : (dia.empty() ? 0 : dia[i]));
            }
        }
        barras_calibracion_ = x.size();
        obs_.calibrar(cfg_, x, den, ok, fs, ds);
        if (flujo_.configurar(f, delta_nu_) != ark::ventanas::ErrorMotorIntrinseco::NINGUNO) return;
        lista_ = obs_.escalas() > 0;
    }

    // Ciclo caliente: una fila (fase y dia solo cuentan si la fuente no activa el calendario). Cero reservas.
    Salida paso_fila(const ark::ventanas::EntradaPasoCxx& e, std::size_t fase = 0, std::int64_t dia = 0) noexcept {
        Salida s{};
        auto b = flujo_.paso(e);
        if (!b) return s;
        s.barra_cerrada = true;
        s.barra = *b;
        const bool cal = flujo_.fuente().fase.activa;
        const bool ok_b = cfg_.usar_protocolo_huecos
            ? ((!b->imputada || b->solo_mercado_quieto) && ((b->flags_hueco & (2u | 64u | 128u)) == 0u))
            : (!b->imputada);
        s.delta = obs_.paso(b->valor, b->nu, ok_b, cal ? b->fase : fase, cal ? b->dia : dia);
        return s;
    }
};

}  // namespace ark::metrologia
