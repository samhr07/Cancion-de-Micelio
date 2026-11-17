// observador_saltos.hpp -- orquestador del bloque de saltos y temperatura (ARK-METRO-01), C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04. Prerregistro: PRUEBAS_JEFE/metrologia/EXPECTATIVAS.md.
//
// theta(t) nativo -> representacion estacionaria y (NIVEL: incrementos) -> purga diurna causal (opcional) -> escalera 2^j ->
// por escala: canal de nivel (GLR + Lorden + resolvente dual del fractal) y canal de temperatura (superestadistica directa + GLR
// sobre log beta^) -> clasificador entre escalas con ficha formal:
//   DERIVA      (enmiendas 3 y 4) un evento que persiste hasta la escala mas gruesa sigue abierto hasta juntar L barras de su escala
//               de origen (la memoria del GLR); si la pendiente de z en esas barras es significativa (alfa, dos colas) en el
//               sentido del cambio -> DERIVA, si no SALTO;
//   confirmacion EN CASCADA (enmienda 1 del prerregistro): se espera la escala j0+1; si confirma (mismo signo, antes del plazo) se
//   espera la siguiente, y asi hasta la mas gruesa calibrada;
//   SALTO       persistio hasta la escala mas gruesa (con j0 la mas gruesa: confirmado = false);
//   TRANSITORIO una confirmacion no llego a tiempo; escala_confirmacion = la escala de persistencia (-1 si ninguna);
//   CAMBIO_TEMPERATURA (enmienda 2) alarma de temperatura confirmada EN EL TIEMPO (el GLR reiniciado vuelve a alarmar, mismo
//               signo y escala, antes de su plazo); precedencia sobre el evento de nivel abierto. Sin confirmar: TRANSITORIO, canal 1.
// Plazo = t_alarma + m_{j0+1} (ceil(kappa t*_{j0+1}(delta^)) + 1), delta^ la amplitud llevada a j0+1, kappa = tau50/t* medido.
// ARL0 en pasos nativos, repartido entre los 2J canales (Bonferroni de tasas): ARL0_j = ARL0 2J / m_j (pasos de la escala j).
// Modo (2026-11-05): DISJUNTO (el de siempre) o DESLIZANTE (detector_deslizante.hpp): el nivel se evalua en cada paso nativo en todas
// las escalas y la temperatura en cada barra de su escala, con h recalibrado sobre la misma construccion y el mismo reparto del ARL0.
// En DESLIZANTE la superestadistica, el regularizador (borde, derivada, norma) y la ventana de forma siguen sobre barras disjuntas.
//
// Ciclo caliente (paso): cero reservas. La calibracion reserva.
#pragma once

#include "ark/metrologia/canal_nivel.hpp"
#include "ark/metrologia/canal_temperatura.hpp"
#include "ark/metrologia/detector_deslizante.hpp"
#include "ark/metrologia/escalera_diadica.hpp"
#include "ark/metrologia/sustitutos.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace ark::metrologia {

struct ConfiguracionObservador {
    Integracion integracion = Integracion::FLUJO;
    AgregacionEscala agregacion = AgregacionEscala::SUMA;   // con NIVEL se fuerza SUMA sobre los incrementos
    double arl0_pasos_nativos = 0.0;   // OBLIGATORIO: una falsa alarma de la familia cada ARL0 pasos nativos
    double alfa = ALFA_PROYECTO;
    bool purga = false;                // purga diurna causal (no se aplica a RAZON_DE_SUMAS)
    std::size_t bins_fase = 1;         // fases por dia de la purga
    std::size_t dias_minimos_purga = 20;   // ceil(1/alfa)
    bool medir_tau50 = true;
    std::size_t replicas_tau50 = 64;
    std::uint64_t semilla = 20261104;
    // NULOS = ENTRADAS del usuario (decision del Comandante, 2026-11-04): ningun nulo va fijo dentro del operador.
    std::uint32_t nulos_separacion = NULO_NINGUNO;   // contra que se contrasta la separacion de escalas (vacio: no contrastada)
    std::size_t replicas_nulo = 199;                 // sustitutos por nulo
    std::uint32_t nulo_nivel = NULO_NINGUNO;         // nulo de calibracion del canal de nivel: NINGUNO = la historia tal cual;
                                                     // PERMUTACION o SIGNOS = la historia sustituida (BLOQUES con b = N_REG).
                                                     // NINGUNO: "evento" = mas raro que una alarma cada ARL0 en la propia historia
                                                     // (las alarmas de la calibracion son el presupuesto). PERMUTACION y SIGNOS
                                                     // SOLO para series sin memoria: con memoria el ARL0 no se cumple (x28-30 en un
                                                     // AR(1) phi 0.88). Entrada del usuario segun la serie (2026-11-11, README).
    ModoDeteccion modo = ModoDeteccion::DISJUNTO;    // DESLIZANTE: detector_deslizante.hpp (los dos quedan)
    std::size_t escalas_max = 0;                     // 0: todas las que la historia soporte; > 0: tope (la familia del crudo lo fija en
                                                     // las escalas de la historia sola: sus sustitutos por bloques no sostienen mas)
};

template <std::size_t J_MAX = 24, std::size_t L = 120, std::size_t N_REG = 120, std::size_t B_MAX = 1440, std::size_t EV_MAX = 8>
class ObservadorSaltos {
public:
    struct Salida {
        std::int64_t t = -1;
        std::size_t niveles = 0;                         // niveles emitidos en este paso (0..niveles-1, solo los calibrados)
        std::array<SalidaNivel, J_MAX> nivel{};
        std::array<SalidaTemperatura, J_MAX> temp{};
        std::array<FichaEvento, EV_MAX> eventos{};
        std::size_t n_eventos = 0;
        bool evento_abierto = false;
    };

private:
    ConfiguracionObservador cfg_{};
    EscaleraDiadica<J_MAX> esc_{};
    PurgaDiurna<B_MAX> purga_{};
    std::array<CanalNivel<L, N_REG>, J_MAX> niv_{};
    std::array<CanalTemperatura<L>, J_MAX> tmp_{};
    EscaleraDeslizante<J_MAX> escd_{};                         // solo en DESLIZANTE
    std::array<CanalNivelDeslizante<L>, J_MAX> nivd_{};
    std::array<CanalTemperaturaDeslizante<L>, J_MAX> tmpd_{};
    std::array<double, J_MAX> arl0_{};
    std::size_t J_ = 0;
    double x_prev_ = NaN;
    bool hay_prev_ = false;
    std::int64_t t_ = -1;
    bool abierto_ = false, tabierto_ = false;
    FichaEvento ev_{}, tev_{};
    double fk_ = 0, fy_ = 0, fky_ = 0, fkk_ = 0;   // acumuladores de la forma (escala de origen, desde la alarma)
    std::size_t fn_ = 0;
    bool persistente_ = false;                     // la cascada llego a la escala mas gruesa: se espera la ventana de forma

    [[nodiscard]] bool deslizante() const noexcept { return cfg_.modo == ModoDeteccion::DESLIZANTE; }

    // cierre de un evento que persistio: SALTO o DERIVA segun la pendiente acumulada
    void cerrar_persistente(Salida& s) noexcept {
        ev_.confirmado = ev_.escala_confirmacion >= 0;   // false si la alarma nacio en la escala mas gruesa
        ev_.t_cierre = t_;
        ev_.barras_forma = fn_;
        ev_.estado = EstadoSaltos::SALTO;
        if (fn_ >= 3) {
            const double n = static_cast<double>(fn_), kb = fk_ / n, yb = fy_ / n;
            const double sxx = fkk_ - n * kb * kb, sxy = fky_ - n * kb * yb;
            const double tau = ficha_nivel(static_cast<std::size_t>(ev_.escala_origen)).nulo.tau_int;
            if (sxx > 0.0) {
                ev_.t_pendiente_cierre = (sxy / sxx) / std::sqrt(tau / sxx);
                const double zc = est::norm_ppf(1.0 - 0.5 * cfg_.alfa);
                if (std::abs(ev_.t_pendiente_cierre) > zc && (ev_.t_pendiente_cierre > 0.0) == (ev_.alarma.signo > 0)) ev_.estado = EstadoSaltos::DERIVA;
            }
        }
        emitir(s, ev_);
        abierto_ = false;
    }
    std::vector<std::vector<std::vector<double>>> hist_;   // historia de calibracion por escala (solo calibracion)

    AgregacionEscala ag_efectiva() const noexcept {
        return cfg_.integracion == Integracion::NIVEL ? AgregacionEscala::SUMA : cfg_.agregacion;
    }

    // representacion estacionaria nativa (con estado: incremento y purga)
    bool representar(double x, double den, bool valido, std::size_t fase, std::int64_t dia, double& y, double& dy) noexcept {
        dy = den;
        if (cfg_.integracion == Integracion::NIVEL) {
            const bool ok = valido && std::isfinite(x);
            y = (ok && hay_prev_) ? x - x_prev_ : NaN;
            const bool v = ok && hay_prev_;
            hay_prev_ = ok;
            if (ok) x_prev_ = x;
            valido = v;
        } else {
            y = x;
            valido = valido && std::isfinite(x);
        }
        if (cfg_.purga && ag_efectiva() != AgregacionEscala::RAZON_DE_SUMAS) {
            bool vp = false;
            const double yp = purga_.paso(valido ? y : NaN, fase, dia, vp);
            y = yp;
            valido = valido && vp;
        }
        return valido;
    }

    void emitir(Salida& s, const FichaEvento& e) const noexcept {
        if (s.n_eventos < EV_MAX) s.eventos[s.n_eventos++] = e;
    }

    // plazo de confirmacion en la escala j0+1 para una amplitud por barra `amplitud` de la escala j0, contado desde t0
    void fijar_plazo(FichaEvento& e, std::size_t j0, double amplitud, std::int64_t t0) const noexcept {
        const std::size_t j1 = j0 + 1;
        e.escala_esperada = static_cast<int>(j1);
        const auto& f1 = ficha_nivel(j1);
        const double factor = ag_efectiva() == AgregacionEscala::SUMA ? 2.0 : 1.0;
        const double delta = std::abs(amplitud) * factor / f1.nulo.s0;
        const double ts = std::log(arl0_[j1]) * 2.0 * f1.nulo.tau_int / (delta * delta);
        std::size_t best = DELTA_REPORTE.size();
        double dist = INFINITY;
        for (std::size_t i = 0; i < DELTA_REPORTE.size(); ++i) {
            if (!std::isfinite(f1.det.tau50[i]) || !(f1.det.t_lorden[i] > 0.0)) continue;
            const double d = std::abs(std::log(DELTA_REPORTE[i] / delta));
            if (d < dist) { dist = d; best = i; }
        }
        e.t_lorden_confirmacion = ts;
        if (best < DELTA_REPORTE.size()) {
            e.kappa = std::max(1.0, f1.det.tau50[best] / f1.det.t_lorden[best]);
            e.plazo = t0 + f1.m * (static_cast<std::int64_t>(std::ceil(e.kappa * ts)) + 1);
        } else {
            e.kappa = NaN;
            e.plazo = t0 + f1.m * static_cast<std::int64_t>(L);   // sin tau50 medido: la memoria del GLR (declarado)
        }
    }

    // Calibracion del modo DESLIZANTE, tras la disjunta (reusa su nulo). yn, dyn, vn: la representacion nativa de la historia.
    void calibrar_deslizante(const std::vector<double>& yn, const std::vector<double>& dyn, const std::vector<std::uint8_t>& vn) {
        escd_.configurar(ag_efectiva(), J_);
        // tramos nativos validos (los sustitutos actuan dentro de cada tramo, como en el disjunto)
        std::vector<std::vector<double>> tramos_idx;
        double mu_y = 0.0;
        std::size_t ny = 0;
        for (std::size_t i = 0; i < yn.size(); ++i) {
            if (!vn[i]) continue;
            if (i == 0 || !vn[i - 1]) tramos_idx.emplace_back();
            tramos_idx.back().push_back(static_cast<double>(i));
            mu_y += yn[i]; ++ny;
        }
        if (ny) mu_y /= static_cast<double>(ny);
        std::vector<double> serie;
        for (std::size_t j = 0; j < J_; ++j) {
            const std::int64_t m = std::int64_t{1} << j;
            // historia nativa del nulo: la propia (NINGUNO) o su sustituto POR TRAMOS NATIVOS (DIVERGE del disjunto, declarado: el
            // disjunto sustituye las barras de la escala con b = N_REG; aqui los pasos nativos con b = N_REG 2^j, y numerador y
            // denominador viajan juntos). Despues, la barra deslizante con el mismo codigo que el ciclo caliente.
            std::vector<double> ys, ds;
            std::vector<std::uint8_t> vs;
            if (cfg_.nulo_nivel == NULO_NINGUNO) {
                ys = yn; ds = dyn; vs = vn;
            } else {
                std::mt19937_64 rng(cfg_.semilla ^ (0xA24BAED4963EE407ULL * (j + 1)));
                const bool signos = cfg_.nulo_nivel == NULO_SIGNOS;
                const auto orden = signos ? tramos_idx : sustituto(tramos_idx, cfg_.nulo_nivel, N_REG * static_cast<std::size_t>(m), 0.0, rng);
                std::bernoulli_distribution B(0.5);
                for (const auto& g : orden) {
                    for (double di : g) {
                        const std::size_t i = static_cast<std::size_t>(di);
                        ys.push_back(signos ? mu_y + (B(rng) ? 1.0 : -1.0) * (yn[i] - mu_y) : yn[i]);
                        ds.push_back(dyn[i]);
                        vs.push_back(1);
                    }
                    ys.push_back(NaN); ds.push_back(NaN); vs.push_back(0);   // los tramos no se tocan
                }
            }
            EscaleraDeslizante<J_MAX> ed;
            ed.configurar(ag_efectiva(), j + 1);
            serie.assign(ys.size(), NaN);
            for (std::size_t i = 0; i < ys.size(); ++i) {
                const auto d = ed.paso(ys[i], ds[i], vs[i] != 0);
                serie[i] = d.valido[j] ? d.valor[j] : NaN;
            }
            ConfiguracionNivel cn{};
            cn.arl0_pasos = arl0_[j]; cn.alfa = cfg_.alfa; cn.medir_tau50 = cfg_.medir_tau50;
            cn.replicas_tau50 = cfg_.replicas_tau50; cn.semilla = cfg_.semilla;
            nivd_[j].calibrar(niv_[j].ficha(), serie, ag_efectiva(), cn);
            ConfiguracionTemperatura ct{};
            ct.arl0_pasos = arl0_[j]; ct.alfa = cfg_.alfa; ct.semilla = cfg_.semilla;
            ct.nulos_separacion = cfg_.nulos_separacion; ct.replicas_nulo = cfg_.replicas_nulo;
            tmpd_[j].calibrar(tmp_[j], hist_[j], arl0_[j], ct);
        }
    }

public:
    [[nodiscard]] std::size_t escalas() const noexcept { return J_; }
    [[nodiscard]] const CanalNivel<L, N_REG>& nivel(std::size_t j) const noexcept { return niv_[j]; }
    [[nodiscard]] const CanalTemperatura<L>& temperatura(std::size_t j) const noexcept { return tmp_[j]; }
    [[nodiscard]] const std::vector<std::vector<double>>& historia(std::size_t j) const { return hist_.at(j); }
    [[nodiscard]] double arl0_escala(std::size_t j) const noexcept { return arl0_[j]; }
    [[nodiscard]] ModoDeteccion modo() const noexcept { return cfg_.modo; }
    [[nodiscard]] const CanalNivelDeslizante<L>& nivel_deslizante(std::size_t j) const noexcept { return nivd_[j]; }
    [[nodiscard]] const CanalTemperaturaDeslizante<L>& temperatura_deslizante(std::size_t j) const noexcept { return tmpd_[j]; }
    // La ficha con que el modo en curso detecta (la que entra en los plazos y en las fichas de los eventos).
    [[nodiscard]] const FichaEscala& ficha_nivel(std::size_t j) const noexcept { return deslizante() ? nivd_[j].ficha() : niv_[j].ficha(); }
    [[nodiscard]] const FichaEscala& ficha_temperatura(std::size_t j) const noexcept { return deslizante() ? tmpd_[j].ficha_glr() : tmp_[j].ficha_glr(); }

    // Calibra con la historia nativa (todo <= el instante de calibracion). den, valido, fase y dia pueden ir vacios.
    void calibrar(const ConfiguracionObservador& cfg, std::span<const double> x, std::span<const double> den = {},
                  std::span<const std::uint8_t> valido = {}, std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        if (!(cfg.arl0_pasos_nativos > 1.0)) throw std::invalid_argument("ConfiguracionObservador.arl0_pasos_nativos es obligatorio");
        cfg_ = cfg;
        esc_ = EscaleraDiadica<J_MAX>(ag_efectiva());
        purga_ = PurgaDiurna<B_MAX>(cfg.bins_fase, cfg.dias_minimos_purga);
        hay_prev_ = false;
        t_ = -1;
        hist_.assign(J_MAX, {});
        std::array<bool, J_MAX> corriendo{};
        auto& yn = yn_;   // representacion nativa (solo DESLIZANTE: la construccion deslizante del nulo; se guarda para recalibrar)
        auto& dyn = dyn_;
        auto& vn = vn_;
        yn.clear(); dyn.clear(); vn.clear();
        const bool desl = cfg.modo == ModoDeteccion::DESLIZANTE;
        if (desl) { yn.reserve(x.size()); dyn.reserve(x.size()); vn.reserve(x.size()); }
        for (std::size_t i = 0; i < x.size(); ++i) {
            double y, dy;
            const bool v = representar(x[i], den.empty() ? NaN : den[i], valido.empty() ? true : valido[i] != 0,
                                       fase.empty() ? 0 : fase[i], dia.empty() ? 0 : dia[i], y, dy);
            if (desl) { yn.push_back(y); dyn.push_back(dy); vn.push_back(v); }
            const auto s = esc_.paso(y, dy, v);
            for (std::size_t j = 0; j < s.niveles; ++j) {
                if (!s.valido[j]) { corriendo[j] = false; continue; }
                if (!corriendo[j]) { hist_[j].emplace_back(); corriendo[j] = true; }
                hist_[j].back().push_back(s.valor[j]);
            }
        }
        const std::size_t minimo = static_cast<std::size_t>(std::ceil(1.0 / cfg.alfa)) * N_REG;
        J_ = 0;
        for (std::size_t j = 0; j < J_MAX; ++j) {
            std::size_t n = 0;
            for (const auto& g : hist_[j]) n += g.size();
            if (n < minimo) break;
            J_ = j + 1;
        }
        if (cfg.escalas_max > 0 && J_ > cfg.escalas_max) J_ = cfg.escalas_max;
        calibrar_canales();
        reiniciar_flujo();
    }

    // Recalibra los umbrales de todos los canales con otro ARL0 de la familia, sobre la MISMA historia ya preparada (la escalera, la
    // purga y el nulo no cambian). Igual, bit a bit, a calibrar de nuevo con ese ARL0. Lo usa el observador conjunto para repartir el
    // ARL0 entre lentes por numero de canales. Fuera del ciclo caliente.
    void recalibrar_arl0(double arl0) {
        if (!(arl0 > 1.0)) throw std::invalid_argument("recalibrar_arl0: ARL0 > 1");
        cfg_.arl0_pasos_nativos = arl0;
        calibrar_canales();
        reiniciar_flujo();
    }

private:
    std::vector<double> yn_, dyn_;
    std::vector<std::uint8_t> vn_;

    // Umbrales, nulos y fichas de cada canal desde la historia preparada (hist_, J_ y, en DESLIZANTE, yn_/dyn_/vn_).
    void calibrar_canales() {
        const auto& cfg = cfg_;
        const bool desl = cfg.modo == ModoDeteccion::DESLIZANTE;
        for (std::size_t j = 0; j < J_; ++j) {
            const std::int64_t m = std::int64_t{1} << j;
            arl0_[j] = cfg.arl0_pasos_nativos * 2.0 * static_cast<double>(J_) / static_cast<double>(m);
            ConfiguracionNivel cn{};
            cn.arl0_pasos = arl0_[j]; cn.alfa = cfg.alfa; cn.medir_tau50 = cfg.medir_tau50;
            cn.replicas_tau50 = cfg.replicas_tau50; cn.semilla = cfg.semilla;
            if (cfg.nulo_nivel == NULO_NINGUNO) {
                niv_[j].calibrar(static_cast<int>(j), m, hist_[j], cn);
            } else {
                std::mt19937_64 rng(cfg.semilla ^ (0xA24BAED4963EE407ULL * (j + 1)));
                const auto mo = est::momentos(hist_[j]);
                niv_[j].calibrar(static_cast<int>(j), m, sustituto(hist_[j], cfg.nulo_nivel, N_REG, mo.mu, rng), cn);
            }
            ConfiguracionTemperatura ct{};
            ct.arl0_pasos = arl0_[j]; ct.alfa = cfg.alfa; ct.semilla = cfg.semilla;
            ct.nulos_separacion = cfg.nulos_separacion; ct.replicas_nulo = cfg.replicas_nulo;
            tmp_[j].calibrar(static_cast<int>(j), m, hist_[j], ct);
        }
        if (desl) calibrar_deslizante(yn_, dyn_, vn_);
    }

public:

    // Vacia la memoria de los detectores y de la escalera (la purga conserva su perfil aprendido). El indice nativo vuelve a 0.
    void reiniciar_flujo() noexcept {
        esc_.reiniciar();
        for (std::size_t j = 0; j < J_; ++j) { niv_[j].reiniciar_flujo(); tmp_[j].reiniciar_flujo(); }
        if (deslizante()) {
            escd_.reiniciar();
            for (std::size_t j = 0; j < J_; ++j) { nivd_[j].reiniciar_flujo(); tmpd_[j].reiniciar_flujo(); }
        }
        hay_prev_ = false;
        abierto_ = false;
        tabierto_ = false;
        persistente_ = false;
        t_ = -1;
    }

    // Ciclo caliente: una observacion nativa. Cero reservas.
    Salida paso(double x, double den = NaN, bool valido = true, std::size_t fase = 0, std::int64_t dia = 0) noexcept {
        Salida s{};
        s.t = ++t_;
        double y, dy;
        const bool v = representar(x, den, valido, fase, dia, y, dy);
        const auto e = esc_.paso(y, dy, v);
        if (!deslizante()) {
            s.niveles = std::min(e.niveles, J_);
            for (std::size_t j = 0; j < s.niveles; ++j) {
                s.nivel[j] = niv_[j].paso(e.valor[j], e.valido[j], t_);
                s.temp[j] = tmp_[j].paso(e.valor[j], e.valido[j], t_);
            }
        } else {
            // nivel: todas las escalas en cada paso nativo; regularizador y temperatura: en la barra disjunta de la escala (la que cierra)
            const auto d = escd_.paso(y, dy, v);
            s.niveles = J_;
            const std::size_t ne = std::min(e.niveles, J_);
            for (std::size_t j = 0; j < J_; ++j) {
                s.nivel[j] = nivd_[j].paso(d.valor[j], d.valido[j], t_);
                if (j < ne) {
                    const auto r = niv_[j].paso(e.valor[j], e.valido[j], t_);   // su GLR disjunto corre pero no se usa (declarado)
                    s.nivel[j].norma = r.norma; s.nivel[j].borde = r.borde; s.nivel[j].derivada = r.derivada; s.nivel[j].ultimo_salto = r.ultimo_salto;
                    s.temp[j] = tmpd_[j].paso(e.valor[j], e.valido[j], t_);
                }
            }
        }
        if (abierto_) {
            const std::size_t j0 = static_cast<std::size_t>(ev_.escala_origen);
            if (j0 < e.niveles && e.valido[j0]) {   // la ventana de forma va sobre barras DISJUNTAS de la escala de origen
                const auto& f0 = ficha_nivel(j0);
                const double z = (e.valor[j0] - f0.nulo.mu0) / f0.nulo.s0, k = static_cast<double>(++fn_);
                fk_ += k; fy_ += z; fky_ += k * z; fkk_ += k * k;
            }
        }
        // 1. alarmas de nivel, de fina a gruesa (una confirmacion en j0+1 puede llegar en el mismo paso)
        for (std::size_t j = 0; j < s.niveles; ++j) {
            const auto& a = s.nivel[j].alarma;
            if (!a.hay) continue;
            if (abierto_) {
                if (static_cast<int>(j) == ev_.escala_esperada && a.signo == ev_.alarma.signo) {
                    ev_.escala_confirmacion = static_cast<int>(j);
                    if (j + 1 >= J_ || !ficha_nivel(j + 1).detectable) {
                        persistente_ = true;           // se espera la ventana de forma (L barras de la escala de origen)
                        ev_.escala_esperada = -1;
                    } else {
                        fijar_plazo(ev_, j, a.amplitud, t_);   // cascada: se espera la siguiente escala
                    }
                }
                continue;
            }
            FichaEvento n{};
            n.canal = 0;
            n.escala_origen = static_cast<int>(j);
            n.alarma = a;
            n.ficha_origen = ficha_nivel(j);
            ev_ = n;
            abierto_ = true;
            if (j + 1 >= J_ || !ficha_nivel(j + 1).detectable) {
                persistente_ = true;   // la escala mas gruesa no tiene a quien pedir confirmacion: solo la ventana de forma (enmienda 5)
                ev_.escala_esperada = -1;
            } else {
                fijar_plazo(ev_, j, a.amplitud, a.t);
                persistente_ = false;
            }
            fk_ = fy_ = fky_ = fkk_ = 0.0;
            fn_ = 0;
        }
        // 2. temperatura: confirmacion EN EL TIEMPO (enmienda 2). Una alarma abre el evento; se confirma si el GLR, reiniciado,
        //    vuelve a alarmar en la misma escala y con el mismo signo antes del plazo. Confirmado: CAMBIO_TEMPERATURA, con
        //    precedencia sobre el evento de nivel abierto. Sin confirmar: TRANSITORIO del canal de temperatura.
        //    DESLIZANTE (enmienda 1 del deslizante): el canal acumula y emite aparte la alarma de confirmacion (s.temp[j].confirmacion);
        //    solo esa confirma, y una apertura con un evento abierto se ignora (como en el disjunto).
        for (std::size_t j = 0; j < s.niveles; ++j) {
            const auto& a = s.temp[j].alarma;
            if (!a.hay) continue;
            const bool conf = deslizante() && s.temp[j].confirmacion;
            if (tabierto_ || conf) {
                if (tabierto_ && (!deslizante() || conf) && static_cast<int>(j) == tev_.escala_origen && a.signo == tev_.alarma.signo) {
                    tev_.estado = EstadoSaltos::CAMBIO_TEMPERATURA;
                    tev_.confirmado = true;
                    tev_.escala_confirmacion = static_cast<int>(j);
                    tev_.t_cierre = t_;
                    if (abierto_) { ev_.estado = EstadoSaltos::CAMBIO_TEMPERATURA; ev_.t_cierre = t_; emitir(s, ev_); abierto_ = false; }
                    emitir(s, tev_);
                    tabierto_ = false;
                }
                continue;
            }
            const auto& fg = ficha_temperatura(j);
            FichaEvento n{};
            n.canal = 1;
            n.escala_origen = static_cast<int>(j);
            n.escala_esperada = static_cast<int>(j);
            n.alarma = a;
            n.ficha_origen = fg;
            const double delta = std::abs(a.amplitud) / fg.nulo.s0;
            n.kappa = 1.0;   // sin tau50 medido para el canal de temperatura (declarado)
            if (!deslizante()) {
                const double ts = std::log(fg.det.arl0_pasos) * 2.0 * fg.nulo.tau_int / (delta * delta);
                n.t_lorden_confirmacion = ts;
                n.plazo = t_ + fg.m * (static_cast<std::int64_t>(std::ceil(ts)) + 1);
            } else {   // plazo con el umbral CALIBRADO de confirmacion (la misma funcion que usa el canal)
                n.t_lorden_confirmacion = 2.0 * fg.det.h_confirmacion * fg.nulo.tau_int / (delta * delta);
                n.plazo = t_ + fg.m * plazo_elementos_temperatura(fg.det.h_confirmacion, fg.nulo.tau_int, a.amplitud, fg.nulo.s0, L);
            }
            tev_ = n;
            tabierto_ = true;
        }
        // 3. plazos vencidos; ventana de forma completa
        if (abierto_ && persistente_ && fn_ >= L) cerrar_persistente(s);
        if (abierto_ && !persistente_ && t_ >= ev_.plazo) {
            ev_.estado = EstadoSaltos::TRANSITORIO;
            ev_.t_cierre = t_;
            emitir(s, ev_);
            abierto_ = false;
        }
        if (tabierto_ && t_ >= tev_.plazo) {
            tev_.estado = EstadoSaltos::TRANSITORIO;
            tev_.t_cierre = t_;
            emitir(s, tev_);
            tabierto_ = false;
        }
        s.evento_abierto = abierto_ || tabierto_;
        return s;
    }
};

}  // namespace ark::metrologia
