// observador_memoria_dfa.hpp -- lente de memoria con el sensor DFA-1 solapado (ARK-FRAC-01) sobre cualquier observable, con el
// operador delta (DISJUNTO o DESLIZANTE) sobre su H^. Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-05, por orden del Comandante:
// "puedes usar Whittle, pero no descartes el DFA: falta calibrarlo respecto a otra variable que no sean los retornos".
//
// Es la misma tuberia que ObservadorMemoria (ELW), con el sensor del Bloque Fractal en lugar del ELW, para que las dos lentes sean
// intercambiables en el observador conjunto (observador_delta_fractal.hpp):
//   observable crudo -> SensorDFACausal<W, S> (ark/fractal, sin tocar) -> H^ cada S muestras -> ObservadorSaltos (FLUJO, MEDIA).
// Reutiliza: el sensor y dfa1_solapado del fractal, los sustitutos del crudo, los tramos y los estadisticos de ObservadorMemoria
// (memoria_detalle), el ObservadorSaltos, la fuente comun de filas y el motor de reloj intrinseco (la ventana del motor = W del DFA).
//
// NULOS = ENTRADAS: nulo_sensor (con que H^ se calibra el operador; NINGUNO = el H^ real; PERMUTACION / SIGNOS / BLOQUES b = W =
// H^ de un sustituto del CRUDO) y nulos_contraste (H^ real contra los nulos pedidos). Los umbrales del operador se autocalibran sobre
// ese H^ (mu0, s0, tau_int y h por ARL0): no hay numero puesto a mano para el observable que sea.
// Diferencias con el ELW, declaradas: W y S son los del DFA (plantilla; 50 y 2 los del Bloque Fractal), no se derivan de una banda;
// sin suelo de ruido derivado (el ruido de H^ lo mide el propio operador: s0 de su nulo).
// Ciclo caliente (paso, paso_fila): cero reservas.
#pragma once

#include "ark/fractal/sensor_dfa_causal.hpp"
#include "ark/metrologia/observador_memoria.hpp"
#include "ark/metrologia/observador_saltos.hpp"
#include "ark/ventanas/flujo_barras.hpp"
#include "ark/ventanas/motor_reloj_intrinseco.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

namespace ark::metrologia {

struct ConfiguracionMemoriaDFA {
    Integracion orden = Integracion::FLUJO;          // del observable crudo (FLUJO: H = alfa del DFA; NIVEL: H = alfa - 1)
    ark::fractal::ConfiguracionDFA dfa{};            // escalas del DFA (los valores del Bloque Fractal por omision)
    std::uint32_t nulo_sensor = NULO_NINGUNO;        // uno solo: con que H^ se calibra el operador
    std::uint32_t nulos_contraste = NULO_NINGUNO;    // mascara: contra que nulos se compara el H^ real
    std::uint64_t semilla = 20261105;
    ConfiguracionObservador delta{};                 // arl0_pasos_nativos OBLIGATORIO (pasos de H^); integracion y agregacion se fuerzan
};

template <std::size_t W = 50, std::size_t S = 2, std::size_t J_MAX = 24, std::size_t L = 120, std::size_t N_REG = 120,
          std::size_t B_MAX = 1440, std::size_t EV_MAX = 8>
class ObservadorMemoriaDFA {
public:
    using Delta = ObservadorSaltos<J_MAX, L, N_REG, B_MAX, EV_MAX>;
    using Configuracion = ConfiguracionMemoriaDFA;
    static constexpr const char* NOMBRE_SENSOR = "DFA-1 solapado";
    struct Salida {
        bool emitio = false;                     // el sensor evaluo en esta muestra
        ark::fractal::ResultadoSensor sensor{};
        typename Delta::Salida delta{};          // solo si emitio
    };
    using FuenteFilas = ark::ventanas::FuenteFilasCxx<W>;

private:
    ConfiguracionMemoriaDFA cfg_{};
    ark::fractal::ConfiguracionDFA cdfa_{};      // la que usa el sensor (min_muestras ajustado a W, como MotorDFATTF)
    ark::fractal::SensorDFACausal<W, S> sen_{};
    std::array<double, W> perfil_{};
    Delta delta_{};
    FichaFlujoH real_{}, calibracion_{};
    std::vector<FichaFlujoH> contraste_;
    bool lista_ = false;
    // camino por filas
    FuenteFilas fuente_{};
    double delta_nu_ = NaN;
    bool por_filas_ = false;
    ark::ventanas::MotorVentanasIntrinsecoCxx<W> motor_{};
    ark::ventanas::VentanaIntrinsecaCxx<W> vent_{};

    [[nodiscard]] ark::fractal::Integracion orden_fractal() const noexcept {
        return cfg_.orden == Integracion::NIVEL ? ark::fractal::Integracion::NIVEL : ark::fractal::Integracion::FLUJO;
    }

public:
    [[nodiscard]] bool lista() const noexcept { return lista_; }
    [[nodiscard]] const Delta& delta() const noexcept { return delta_; }
    [[nodiscard]] const FichaFlujoH& flujo_real() const noexcept { return real_; }
    [[nodiscard]] const FichaFlujoH& flujo_calibracion() const noexcept { return calibracion_; }
    [[nodiscard]] const std::vector<FichaFlujoH>& contraste() const noexcept { return contraste_; }
    [[nodiscard]] std::size_t ventana() const noexcept { return W; }
    [[nodiscard]] std::size_t paso_sensor() const noexcept { return S; }
    [[nodiscard]] double delta_nu() const noexcept { return delta_nu_; }

    // Serie de H^ del sensor sobre una serie cruda (reinicia el sensor). Misma forma que la del ELW.
    SerieMemoria serie(std::span<const double> x, std::span<const std::uint8_t> ok, std::span<const std::size_t> fase,
                       std::span<const std::int64_t> dia) {
        SerieMemoria s;
        sen_.reiniciar();
        for (std::size_t i = 0; i < x.size(); ++i) {
            const bool v = (ok.empty() || ok[i]) && std::isfinite(x[i]);
            const auto r = sen_.paso(v ? x[i] : NaN);
            if (!r) continue;
            s.h.push_back(r->valido ? r->h : NaN);
            s.ok.push_back(r->valido ? 1 : 0);
            s.fase.push_back(fase.empty() ? 0 : fase[i]);
            s.dia.push_back(dia.empty() ? 0 : dia[i]);
        }
        sen_.reiniciar();
        return s;
    }

    void calibrar(const ConfiguracionMemoriaDFA& cfg, std::span<const double> x, std::span<const std::uint8_t> ok = {},
                  std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        if (!(cfg.delta.arl0_pasos_nativos > 1.0)) throw std::invalid_argument("ConfiguracionMemoriaDFA.delta.arl0_pasos_nativos es obligatorio");
        if (cfg.nulo_sensor != NULO_NINGUNO && cfg.nulo_sensor != NULO_PERMUTACION && cfg.nulo_sensor != NULO_BLOQUES && cfg.nulo_sensor != NULO_SIGNOS)
            throw std::invalid_argument("nulo_sensor: uno solo");
        cfg_ = cfg;
        lista_ = false;
        contraste_.clear();
        sen_.configurar(orden_fractal(), cfg.dfa);
        cdfa_ = cfg.dfa;
        cdfa_.min_muestras = std::min(cdfa_.min_muestras, static_cast<int>(W));   // lo mismo que hace el sensor al configurarse
        std::mt19937_64 rng(cfg.semilla);
        const auto sr = serie(x, ok, fase, dia);
        real_ = memoria_detalle::estadisticos(sr, NULO_NINGUNO, cfg.delta.alfa);
        SerieMemoria sc;
        if (cfg.nulo_sensor == NULO_NINGUNO) sc = sr;
        else {
            const auto y = memoria_detalle::sustituto_crudo(x, ok, cfg.orden, cfg.nulo_sensor, W, rng);
            sc = serie(y, ok, fase, dia);
        }
        calibracion_ = memoria_detalle::estadisticos(sc, cfg.nulo_sensor, cfg.delta.alfa);
        for (std::uint32_t nb : NULOS_TODOS) {
            if (!(cfg.nulos_contraste & nb)) continue;
            const auto y = memoria_detalle::sustituto_crudo(x, ok, cfg.orden, nb, W, rng);
            auto f = memoria_detalle::estadisticos(serie(y, ok, fase, dia), nb, cfg.delta.alfa);
            if (f.validos > 1 && real_.validos > 1) {
                f.dif_real_menos_nulo = real_.mu - f.mu;
                const double vr = real_.sd * real_.sd * (std::isfinite(real_.tau_int) ? real_.tau_int : 1.0) / static_cast<double>(real_.validos);
                const double vn = f.sd * f.sd * (std::isfinite(f.tau_int) ? f.tau_int : 1.0) / static_cast<double>(f.validos);
                f.se_dif = std::sqrt(vr + vn);
            }
            contraste_.push_back(f);
        }
        ConfiguracionObservador cd = cfg.delta;
        cd.integracion = Integracion::FLUJO;
        cd.agregacion = AgregacionEscala::MEDIA;
        delta_.calibrar(cd, sc.h, {}, sc.ok, sc.fase, sc.dia);
        sen_.reiniciar();
        delta_.reiniciar_flujo();
        lista_ = delta_.escalas() > 0;
    }

    void reiniciar_flujo() noexcept { sen_.reiniciar(); delta_.reiniciar_flujo(); }
    void recalibrar_arl0(double arl0_pasos_h) {
        if (!lista_) return;
        cfg_.delta.arl0_pasos_nativos = arl0_pasos_h;
        delta_.recalibrar_arl0(arl0_pasos_h);
        sen_.reiniciar();
    }

    Salida paso(double x, bool valido = true, std::size_t fase = 0, std::int64_t dia = 0) noexcept {
        Salida s{};
        if (!lista_) return s;
        const auto r = sen_.paso(valido && std::isfinite(x) ? x : NaN);
        if (!r) return s;
        s.emitio = true;
        s.sensor = *r;
        s.delta = delta_.paso(r->valido ? r->h : NaN, NaN, r->valido, fase, dia);
        return s;
    }

    // ---------------- entrada por filas (la fuente comun y el motor de reloj intrinseco, como el ELW) ----------------
    void calibrar_filas(const ConfiguracionMemoriaDFA& cfg, const FuenteFilas& fuente, std::span<const ark::ventanas::EntradaPasoCxx> filas,
                        std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        fuente_ = fuente;
        por_filas_ = false;
        delta_nu_ = ark::ventanas::resolver_delta_nu<W>(fuente, filas);
        if (!(delta_nu_ > 0.0)) throw std::invalid_argument("delta_nu: no se pudo derivar de la historia (sin filas con sesion)");
        std::vector<double> x;
        std::vector<std::uint8_t> ok;
        std::vector<std::size_t> fs;
        std::vector<std::int64_t> ds;
        {   // barras de la historia (barra a barra con el mismo motor y la misma aduana); fase y dia de la fila que cierra la barra
            auto fl = std::make_unique<ark::ventanas::FlujoBarrasIntrinsecasCxx<W>>();
            if (fl->configurar(fuente, delta_nu_) != ark::ventanas::ErrorMotorIntrinseco::NINGUNO) return;
            for (std::size_t i = 0; i < filas.size(); ++i) {
                auto b = fl->paso(filas[i]);
                if (!b) continue;
                x.push_back(b->valor);
                ok.push_back(b->imputada ? 0 : 1);
                fs.push_back(fuente.fase.activa ? b->fase : (fase.empty() ? 0 : fase[i]));
                ds.push_back(fuente.fase.activa ? b->dia : (dia.empty() ? 0 : dia[i]));
            }
        }
        ConfiguracionMemoriaDFA c = cfg;
        if (fuente.fase.activa) c.delta.bins_fase = fuente.fase.bins_fase();
        calibrar(c, x, ok, fs, ds);
        if (!lista_) return;
        ark::ventanas::ConfiguracionMotorIntrinsecoCxx<W> cm{};
        cm.base = fuente.base;
        cm.W = W;
        cm.S = S;
        cm.delta_nu = delta_nu_;
        cm.agregacion = fuente.agregacion;
        if (motor_.configurar(cm) != ark::ventanas::ErrorMotorIntrinseco::NINGUNO) { lista_ = false; return; }
        por_filas_ = true;
    }

    // Ciclo caliente por filas: aduana + barra + ventana del motor -> DFA (si la ventana no tiene barras imputadas) -> operador delta.
    Salida paso_fila(const ark::ventanas::EntradaPasoCxx& e, std::size_t fase = 0, std::int64_t dia = 0) noexcept {
        Salida s{};
        if (!lista_ || !por_filas_) return s;
        const auto p = motor_.paso(e, vent_);
        if (p.resultado != ark::ventanas::ResultadoPaso::VENTANA_EMITIDA) return s;
        s.emitio = true;
        ark::fractal::ResultadoSensor r{};
        if (!vent_.hay_imputadas)
            r = ark::fractal::dfa1_solapado(vent_.valores.data(), W, sen_.escalas(), sen_.n_escalas(), orden_fractal(), cdfa_, perfil_.data());
        s.sensor = r;
        if (fuente_.fase.activa) {
            const auto fd = ark::dataset::fase_de(fuente_.fase, static_cast<std::int64_t>(vent_.t_cierre_ms));
            fase = fd.fase;
            dia = fd.dia;
        }
        s.delta = delta_.paso(r.valido ? r.h : NaN, NaN, r.valido, fase, dia);
        return s;
    }
};

}  // namespace ark::metrologia
