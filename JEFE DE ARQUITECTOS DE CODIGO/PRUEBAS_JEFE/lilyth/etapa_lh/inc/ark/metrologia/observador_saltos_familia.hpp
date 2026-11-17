// observador_saltos_familia.hpp -- el operador delta sobre un crudo (ObservadorSaltos) con presupuesto de FAMILIA, C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-05, por orden del Comandante: la calibracion de familia de la coordinacion llevada a
// todos los operadores, sin eliminar la arquitectura antigua. Prerregistro:
// JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/calibracion_familia/EXPECTATIVAS_FAMILIA_GENERICA.md (etapa 1).
//
// Presupuestos (Configuracion.presupuesto):
//   BONFERRONI_CANALES  el ObservadorSaltos de siempre con arl0 = ARL0 (el reparto por canal dentro del operador, extrapolado canal a
//                       canal): identico bit a bit a usar el ObservadorSaltos solo.
//   FAMILIA_EPISODIOS   la familia son los 2J canales del operador (una subfamilia, ARL0_sub = ARL0). La falsa alarma es el EPISODIO
//                       (abre cuando el operador abre un evento, de nivel o de temperatura; cierra cuando no queda ninguno abierto). La
//                       rareza comun A por canal (barras) se mide sobre B realizaciones del nulo juntadas como tramos y se biseca
//                       (metrologia/calibrador_familia.hpp).
// Nulo del crudo (decision del Comandante): la realizacion 1 es la HISTORIA tal cual; las 2..B son sustitutos por BLOQUES de la
//   historia (la tupla x, den, valido, fase, dia se permuta por bloques enteros; con NIVEL se permutan los incrementos y se reintegra,
//   y un incremento no valido deja el paso no valido). Bloque derivado b = N_REG 2^(J-1), con J las escalas que calibra la historia
//   sola: la convencion del propio operador para su nulo BLOQUES en la escala mas gruesa; la regla de J garantiza >= 1/alfa bloques.
//   Declarado: las realizaciones no son independientes (la misma historia reordenada). El operador interno conserva las J escalas de
//   la historia sola (escalas_max): el nulo juntado habilitaria escalas mas gruesas, sin estructura en los sustitutos por bloques.
// Ciclo caliente (paso): cero reservas. La calibracion reserva.
#pragma once

#include "ark/metrologia/calibrador_familia.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

namespace ark::metrologia {

template <std::size_t J_MAX = 24, std::size_t L = 120, std::size_t N_REG = 120, std::size_t B_MAX = 1440, std::size_t EV_MAX = 8>
class ObservadorSaltosFamilia {
public:
    using Delta = ObservadorSaltos<J_MAX, L, N_REG, B_MAX, EV_MAX>;
    struct Configuracion {
        double arl0_familia_pasos = 0.0;          // OBLIGATORIO: una falsa alarma (episodio) del operador cada ARL0 pasos nativos
        ConfiguracionObservador delta{};          // el operador (su arl0 lo fija la calibracion)
        PresupuestoFamilia presupuesto = PresupuestoFamilia::BONFERRONI_CANALES;
        ParametrosFamilia familia{};              // solo FAMILIA_EPISODIOS
        std::uint64_t semilla_nulo = 20261105;    // de los sustitutos por BLOQUES
        std::size_t bloque_nulo = 0;              // 0: N_REG 2^(J-1), J de la historia sola
    };
    struct Salida {
        typename Delta::Salida delta{};
        bool episodio_abierto = false, episodio_nuevo = false;
    };

private:
    Configuracion cfg_{};
    std::unique_ptr<Delta> obs_ = std::make_unique<Delta>();
    FichaSubfamilia fsub_{};
    double arl0_pasos_ = NaN;                     // el arl0 del operador interno (pasos nativos)
    std::size_t bloque_ = 0, j_historia_ = 0;
    std::vector<double> nx_, nden_;
    std::vector<std::uint8_t> nval_;
    std::vector<std::size_t> nfase_;
    std::vector<std::int64_t> ndia_;
    std::vector<std::size_t> ini_, fin_;
    bool lista_ = false, epi_ = false;

    // Una realizacion r >= 1: la historia permutada por bloques (con NIVEL, sus incrementos), anadida al nulo juntado.
    void anadir_bloques(std::span<const double> x, std::span<const double> den, std::span<const std::uint8_t> val,
                        std::span<const std::size_t> fase, std::span<const std::int64_t> dia, std::mt19937_64& rng) {
        const std::size_t n = x.size();
        const bool nivel = cfg_.delta.integracion == Integracion::NIVEL;
        auto v_ok = [&](std::size_t i) { return (val.empty() || val[i]) && std::isfinite(x[i]); };
        if (!nivel) {
            const auto idx = familia::permutacion_bloques(n, bloque_, rng);
            for (std::size_t k = 0; k < n; ++k) {
                const std::size_t i = idx[k];
                nx_.push_back(x[i]);
                if (!den.empty()) nden_.push_back(den[i]);
                nval_.push_back(val.empty() ? 1 : val[i]);
                nfase_.push_back(fase.empty() ? 0 : fase[i]);
                ndia_.push_back(dia.empty() ? 0 : dia[i]);
            }
            return;
        }
        // NIVEL: incrementos d_i = x_i - x_{i-1} (i >= 1), permutados por bloques y reintegrados desde x_0
        const auto idx = familia::permutacion_bloques(n > 0 ? n - 1 : 0, bloque_, rng);
        double acc = n ? x[0] : 0.0;
        nx_.push_back(acc);
        if (!den.empty()) nden_.push_back(den[0]);
        nval_.push_back(n && v_ok(0) ? 1 : 0);
        nfase_.push_back(fase.empty() ? 0 : fase[0]);
        ndia_.push_back(dia.empty() ? 0 : dia[0]);
        for (std::size_t k = 1; k < n; ++k) {
            const std::size_t s = 1 + idx[k - 1];
            const bool ok = v_ok(s) && v_ok(s - 1);
            if (ok) acc += x[s] - x[s - 1];
            nx_.push_back(acc);
            if (!den.empty()) nden_.push_back(den[s]);
            nval_.push_back(ok ? 1 : 0);
            nfase_.push_back(fase.empty() ? 0 : fase[s]);
            ndia_.push_back(dia.empty() ? 0 : dia[s]);
        }
    }

    void anadir_union() {
        nx_.push_back(NaN);
        if (!nden_.empty()) nden_.push_back(NaN);
        nval_.push_back(0);
        nfase_.push_back(0);
        ndia_.push_back(0);
    }

    // Corre el operador interno (con sus umbrales actuales) sobre el nulo juntado, cada realizacion desde reiniciar_flujo.
    void correr_nulo(Delta& d, std::vector<std::uint8_t>& act) const {
        act.assign(nx_.size(), 0);
        for (std::size_t r = 0; r < ini_.size(); ++r) {
            d.reiniciar_flujo();
            for (std::size_t p = ini_[r]; p < fin_[r]; ++p) {
                const auto s = d.paso(nx_[p], nden_.empty() ? NaN : nden_[p], nval_[p] != 0, nfase_[p], ndia_[p]);
                act[p] = (s.evento_abierto || s.n_eventos > 0) ? 1 : 0;
            }
        }
        d.reiniciar_flujo();
    }

public:
    [[nodiscard]] bool lista() const noexcept { return lista_; }
    [[nodiscard]] const Delta& delta() const noexcept { return *obs_; }
    [[nodiscard]] const Configuracion& configuracion() const noexcept { return cfg_; }
    [[nodiscard]] const FichaSubfamilia& ficha_subfamilia() const noexcept { return fsub_; }
    [[nodiscard]] double arl0_pasos() const noexcept { return arl0_pasos_; }   // el del operador interno
    [[nodiscard]] std::size_t bloque() const noexcept { return bloque_; }
    [[nodiscard]] std::size_t escalas_historia() const noexcept { return j_historia_; }
    // el nulo juntado (solo FAMILIA_EPISODIOS) y sus tramos [ini, fin)
    [[nodiscard]] const std::vector<double>& nulo_x() const noexcept { return nx_; }
    [[nodiscard]] const std::vector<double>& nulo_den() const noexcept { return nden_; }
    [[nodiscard]] const std::vector<std::uint8_t>& nulo_valido() const noexcept { return nval_; }
    [[nodiscard]] const std::vector<std::size_t>& nulo_fase() const noexcept { return nfase_; }
    [[nodiscard]] const std::vector<std::int64_t>& nulo_dia() const noexcept { return ndia_; }
    [[nodiscard]] const std::vector<std::size_t>& nulo_inicios() const noexcept { return ini_; }
    [[nodiscard]] const std::vector<std::size_t>& nulo_fines() const noexcept { return fin_; }

    // Calibra con la historia nativa (todo <= el instante de calibracion). den, valido, fase y dia pueden ir vacios.
    void calibrar(const Configuracion& cfg, std::span<const double> x, std::span<const double> den = {}, std::span<const std::uint8_t> valido = {},
                  std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        if (!(cfg.arl0_familia_pasos > 1.0)) throw std::invalid_argument("Configuracion.arl0_familia_pasos es obligatorio");
        const std::size_t n = x.size();
        if ((!den.empty() && den.size() != n) || (!valido.empty() && valido.size() != n) || (!fase.empty() && fase.size() != n) ||
            (!dia.empty() && dia.size() != n))
            throw std::invalid_argument("ObservadorSaltosFamilia: longitudes");
        cfg_ = cfg;
        lista_ = false;
        fsub_ = FichaSubfamilia{};
        nx_.clear(); nden_.clear(); nval_.clear(); nfase_.clear(); ndia_.clear(); ini_.clear(); fin_.clear();
        bloque_ = 0;
        j_historia_ = 0;
        auto c0 = cfg.delta;
        c0.arl0_pasos_nativos = cfg.arl0_familia_pasos;
        if (cfg.presupuesto == PresupuestoFamilia::BONFERRONI_CANALES) {
            obs_->calibrar(c0, x, den, valido, fase, dia);
            arl0_pasos_ = cfg.arl0_familia_pasos;
            lista_ = obs_->escalas() > 0;
            epi_ = false;
            return;
        }
        const auto& pf = cfg.familia;
        if (pf.realizaciones_max < 1) throw std::invalid_argument("ParametrosFamilia.realizaciones_max >= 1");
        if (!(pf.alfa_certificacion > 0.0 && pf.alfa_certificacion < 1.0)) throw std::invalid_argument("ParametrosFamilia.alfa_certificacion en (0, 1)");
        // J de la historia sola (fija el bloque)
        {
            auto h = std::make_unique<Delta>();
            h->calibrar(c0, x, den, valido, fase, dia);
            j_historia_ = h->escalas();
        }
        if (j_historia_ == 0) { arl0_pasos_ = NaN; return; }
        bloque_ = cfg.bloque_nulo ? cfg.bloque_nulo : (N_REG << (j_historia_ - 1));
        const std::size_t mm = familia::m_min(pf.alfa_certificacion);
        const double arl0_sub = cfg.arl0_familia_pasos;   // una sola subfamilia
        double b_pedido = 0.0;
        const std::size_t B = familia::realizaciones(mm, arl0_sub, n, pf.realizaciones_max, b_pedido);
        // nulo juntado: historia + B-1 sustitutos por BLOQUES, con una union no valida entre realizaciones
        std::mt19937_64 rng(cfg.semilla_nulo);
        for (std::size_t r = 0; r < B; ++r) {
            if (r) anadir_union();
            ini_.push_back(nx_.size());
            if (r == 0) {
                nx_.insert(nx_.end(), x.begin(), x.end());
                if (!den.empty()) nden_.insert(nden_.end(), den.begin(), den.end());
                for (std::size_t i = 0; i < n; ++i) {
                    nval_.push_back(valido.empty() ? 1 : valido[i]);
                    nfase_.push_back(fase.empty() ? 0 : fase[i]);
                    ndia_.push_back(dia.empty() ? 0 : dia[i]);
                }
            } else {
                anadir_bloques(x, den, valido, fase, dia, rng);
            }
            fin_.push_back(nx_.size());
        }
        // las escalas de la historia sola: las mas gruesas que habilitaria el nulo juntado no tienen estructura en los sustitutos
        auto cj = c0;
        cj.escalas_max = cj.escalas_max ? std::min(cj.escalas_max, j_historia_) : j_historia_;
        obs_->calibrar(cj, nx_, nden_, nval_, nfase_, ndia_);   // umbrales provisionales: la busqueda los recalibra
        const std::size_t J = obs_->escalas();
        if (J == 0) { arl0_pasos_ = NaN; return; }
        auto& f = fsub_;
        familia::preparar_ficha(f, B, mm, b_pedido, pf.realizaciones_max, n, n, arl0_sub);
        f.presente = true;
        f.componentes = 1;
        f.canales = 2 * J;
        f.a_bonferroni = arl0_sub * static_cast<double>(f.canales);
        const double suelo = 2.0 * 2.0 * static_cast<double>(J);   // a >= 2 (la convencion de la coordinacion, declarada)
        std::vector<std::vector<std::uint8_t>> act(1);
        const std::vector<std::uint8_t> miembro{1};
        familia::buscar_a_estrella(f, suelo, pf, [&](double A) {
            obs_->recalibrar_arl0(A / (2.0 * static_cast<double>(J)));
            correr_nulo(*obs_, act[0]);
            return episodios_coalescencia(act, miembro, ini_, fin_);
        });
        arl0_pasos_ = f.a_estrella / (2.0 * static_cast<double>(J));
        obs_->recalibrar_arl0(arl0_pasos_);
        lista_ = true;
        epi_ = false;
    }

    void reiniciar_flujo() noexcept { obs_->reiniciar_flujo(); epi_ = false; }

    // Ciclo caliente: una observacion nativa. Cero reservas.
    Salida paso(double x, double den = NaN, bool valido = true, std::size_t fase = 0, std::int64_t dia = 0) noexcept {
        Salida s{};
        s.delta = obs_->paso(x, den, valido, fase, dia);
        const bool u = s.delta.evento_abierto || s.delta.n_eventos > 0;
        s.episodio_nuevo = u && !epi_;
        epi_ = u;
        s.episodio_abierto = epi_;
        return s;
    }
};

}  // namespace ark::metrologia
