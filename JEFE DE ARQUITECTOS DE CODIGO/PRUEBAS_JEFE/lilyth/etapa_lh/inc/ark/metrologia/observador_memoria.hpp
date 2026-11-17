// observador_memoria.hpp -- operador delta y temperatura sobre la memoria H(t) de cualquier observable (ARK-METRO-01), C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04, por orden del Comandante. Prerregistro:
// JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/metrologia/EXPECTATIVAS_WHITTLE.md.
//
// observable crudo -> SensorWhittleExacto (W, m derivados por derivar_banda: suelo de ruido, techo de la ventana interna, techo de
// localidad) -> H^(t) cada S muestras -> ObservadorSaltos (canal de nivel = saltos de H; canal de temperatura = regimenes de la
// dispersion de H), con FLUJO y agregacion MEDIA (H es intensivo).
//
// NULOS = ENTRADAS (decision del Comandante, 2026-11-04):
//   nulo_sensor     historia con que se calibra el operador sobre H^. NINGUNO = H^ de la historia real. PERMUTACION / SIGNOS /
//                   BLOQUES (b = W) = H^ de un sustituto del OBSERVABLE CRUDO (no de H^: permutar H^ destruiria su propia
//                   autocorrelacion de solape). Con NIVEL el sustituto actua sobre los incrementos y se reintegra.
//   nulos_contraste H^ en cada nulo pedido contra el H^ real (media, sd, ACF1, tau_int y la diferencia con su se).
//   delta.nulos_separacion / delta.nulo_nivel siguen siendo los del operador (sobre la serie H^ que recibe).
// Ciclo caliente (paso): cero reservas. La calibracion reserva.
#pragma once

#include "ark/metrologia/observador_saltos.hpp"
#include "ark/metrologia/sustitutos.hpp"
#include "ark/metrologia/whittle_exacto.hpp"
#include "ark/ventanas/flujo_barras.hpp"
#include "ark/ventanas/motor_reloj_intrinseco.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <stdexcept>
#include <vector>

namespace ark::metrologia {

struct ConfiguracionMemoria {
    Integracion orden = Integracion::FLUJO;    // del observable crudo: H = d^ + 1/2 - orden
    ConfiguracionBandaELW banda{};             // ruido_objetivo OBLIGATORIO; ventana_interna del observable
    std::size_t paso_sensor = 1;               // S: resolucion temporal de H^ (entrada del usuario)
    std::uint32_t nulo_ruido = NULO_NINGUNO;   // OBLIGATORIO (enmienda 3): PERMUTACION o SIGNOS; de el se deriva el suelo de ruido
    std::uint32_t nulo_sensor = NULO_NINGUNO;  // uno solo: con que H^ se calibra el operador
    std::uint32_t nulos_contraste = NULO_NINGUNO;   // mascara: contra que nulos se compara el H^ real
    std::uint64_t semilla = 20261104;
    ConfiguracionObservador delta{};           // arl0_pasos_nativos OBLIGATORIO (en pasos de H^); integracion y agregacion se fuerzan
};

// Estadisticos de un flujo de H^ (real o de un nulo).
struct FichaFlujoH {
    std::uint32_t nulo = NULO_NINGUNO;
    std::size_t pasos = 0, validos = 0, en_borde = 0;
    double mu = NaN, sd = NaN, acf1 = NaN, tau_int = NaN;
    double dif_real_menos_nulo = NaN, se_dif = NaN;   // solo en los nulos de contraste
};

// Serie de H^ que produce un sensor sobre una serie cruda (para calibrar y contrastar).
struct SerieMemoria { std::vector<double> h; std::vector<std::uint8_t> ok; std::vector<std::size_t> fase; std::vector<std::int64_t> dia; std::size_t en_borde = 0; };

namespace memoria_detalle {

// Tramos contiguos validos (indices de inicio y fin exclusivos).
inline std::vector<std::pair<std::size_t, std::size_t>> tramos(std::span<const double> x, std::span<const std::uint8_t> ok) {
    std::vector<std::pair<std::size_t, std::size_t>> t;
    std::size_t i = 0;
    while (i < x.size()) {
        while (i < x.size() && !((ok.empty() || ok[i]) && std::isfinite(x[i]))) ++i;
        const std::size_t a = i;
        while (i < x.size() && (ok.empty() || ok[i]) && std::isfinite(x[i])) ++i;
        if (i > a) t.emplace_back(a, i);
    }
    return t;
}

// Sustituto del observable crudo (copia de x). NIVEL: sobre los incrementos de cada tramo y reintegrado desde su primer valor.
inline std::vector<double> sustituto_crudo(std::span<const double> x, std::span<const std::uint8_t> ok, Integracion orden, std::uint32_t nulo,
                                          std::size_t bloque, std::mt19937_64& rng) {
    std::vector<double> y(x.begin(), x.end());
    if (nulo == NULO_NINGUNO) return y;
    const auto tr = tramos(x, ok);
    std::vector<std::vector<double>> segs;
    double s = 0.0;
    std::size_t n = 0;
    for (const auto& [a, b] : tr) {
        std::vector<double> g;
        if (orden == Integracion::NIVEL) for (std::size_t i = a + 1; i < b; ++i) g.push_back(x[i] - x[i - 1]);
        else for (std::size_t i = a; i < b; ++i) g.push_back(x[i]);
        for (double v : g) { s += v; ++n; }
        segs.push_back(std::move(g));
    }
    const double mu = n ? s / static_cast<double>(n) : 0.0;
    const auto sus = sustituto(segs, nulo, bloque, mu, rng);
    for (std::size_t q = 0; q < tr.size(); ++q) {
        const auto [a, b] = tr[q];
        if (orden == Integracion::NIVEL) {
            double acc = x[a];
            for (std::size_t i = a + 1; i < b; ++i) { acc += sus[q][i - a - 1]; y[i] = acc; }
        } else {
            for (std::size_t i = a; i < b; ++i) y[i] = sus[q][i - a];
        }
    }
    return y;
}

inline FichaFlujoH estadisticos(const SerieMemoria& s, std::uint32_t nulo, double alfa) {
    FichaFlujoH f{};
    f.nulo = nulo;
    f.pasos = s.h.size();
    f.en_borde = s.en_borde;
    std::vector<std::vector<double>> segs;
    bool corriendo = false;
    for (std::size_t i = 0; i < s.h.size(); ++i) {
        if (!s.ok[i]) { corriendo = false; continue; }
        if (!corriendo) { segs.emplace_back(); corriendo = true; }
        segs.back().push_back(s.h[i]);
        ++f.validos;
    }
    const auto mo = est::momentos(segs);
    f.mu = mo.mu;
    f.sd = mo.sd;
    if (mo.n >= 4 && mo.sd > 0.0) {
        const auto ti = est::tau_integrado(segs, mo.mu, mo.sd, alfa, std::min<std::size_t>(mo.n / 4, 20000));
        f.acf1 = ti.acf1;
        f.tau_int = std::isfinite(ti.tau) ? std::max(1.0, ti.tau) : NaN;
    }
    return f;
}

}  // namespace memoria_detalle

template <std::size_t W_MAX = 1024, std::size_t J_MAX = 24, std::size_t L = 120, std::size_t N_REG = 120, std::size_t B_MAX = 1440,
          std::size_t EV_MAX = 8>
class ObservadorMemoria {
public:
    using Delta = ObservadorSaltos<J_MAX, L, N_REG, B_MAX, EV_MAX>;
    struct Salida {
        bool emitio = false;                 // el sensor evaluo en esta muestra
        ResultadoELW sensor{};
        typename Delta::Salida delta{};      // solo si emitio
    };

private:
    ConfiguracionMemoria cfg_{};
    FichaBandaELW banda_{};
    SensorWhittleExacto<W_MAX> sen_{};
    Delta delta_{};
    FichaFlujoH real_{}, calibracion_{};
    std::vector<FichaFlujoH> contraste_;
    bool lista_ = false;

public:
    [[nodiscard]] bool lista() const noexcept { return lista_; }
    [[nodiscard]] const FichaBandaELW& banda() const noexcept { return banda_; }
    [[nodiscard]] const FichaFlujoH& flujo_real() const noexcept { return real_; }
    [[nodiscard]] const FichaFlujoH& flujo_calibracion() const noexcept { return calibracion_; }   // el que calibra el operador
    [[nodiscard]] const std::vector<FichaFlujoH>& contraste() const noexcept { return contraste_; }
    [[nodiscard]] const Delta& delta() const noexcept { return delta_; }
    [[nodiscard]] const ConfiguracionSensorELW& sensor() const noexcept { return sen_.configuracion(); }
    // Interfaz comun de las lentes de memoria (2026-11-05, observador_delta_fractal.hpp): configuracion, ventana W en barras, paso S,
    // y recalibrar el ARL0 del operador delta sobre H^ sin repetir el sensor (el observador conjunto reparte el ARL0 entre lentes).
    using Configuracion = ConfiguracionMemoria;
    static constexpr const char* NOMBRE_SENSOR = "Whittle local exacto";
    [[nodiscard]] std::size_t ventana() const noexcept { return banda_.W; }
    [[nodiscard]] std::size_t paso_sensor() const noexcept { return cfg_.paso_sensor; }
    // Vacia el sensor y el operador delta (camino de serie).
    void reiniciar_flujo() noexcept { sen_.reiniciar(); delta_.reiniciar_flujo(); }
    void recalibrar_arl0(double arl0_pasos_h) {
        if (!lista_) return;
        cfg_.delta.arl0_pasos_nativos = arl0_pasos_h;
        delta_.recalibrar_arl0(arl0_pasos_h);
        sen_.reiniciar();
    }

    // Serie de H^ de un sensor ya configurado sobre una serie cruda (reinicia el sensor).
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
            s.en_borde += r->en_borde ? 1 : 0;
            s.fase.push_back(fase.empty() ? 0 : fase[i]);
            s.dia.push_back(dia.empty() ? 0 : dia[i]);
        }
        sen_.reiniciar();
        return s;
    }

    void calibrar(const ConfiguracionMemoria& cfg, std::span<const double> x, std::span<const std::uint8_t> ok = {},
                  std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        if (!(cfg.delta.arl0_pasos_nativos > 1.0)) throw std::invalid_argument("ConfiguracionMemoria.delta.arl0_pasos_nativos es obligatorio");
        if (cfg.nulo_sensor != NULO_NINGUNO && cfg.nulo_sensor != NULO_PERMUTACION && cfg.nulo_sensor != NULO_BLOQUES && cfg.nulo_sensor != NULO_SIGNOS)
            throw std::invalid_argument("nulo_sensor: uno solo");
        if (cfg.nulo_ruido != NULO_PERMUTACION && cfg.nulo_ruido != NULO_SIGNOS)
            throw std::invalid_argument("nulo_ruido es obligatorio: PERMUTACION o SIGNOS (BLOQUES conserva la memoria dentro de la ventana y no mide ruido)");
        cfg_ = cfg;
        lista_ = false;
        contraste_.clear();
        std::mt19937_64 rng(cfg.semilla);
        // 1. banda sobre la historia cruda (tramos validos, tal cual); el suelo de ruido sale del nulo de ruido (mismos tramos)
        std::vector<std::vector<double>> segs, segs_nulo;
        const auto tr = memoria_detalle::tramos(x, ok);
        const auto yr = memoria_detalle::sustituto_crudo(x, ok, cfg.orden, cfg.nulo_ruido, 1, rng);
        for (const auto& [a, b] : tr) {
            segs.emplace_back(x.begin() + static_cast<std::ptrdiff_t>(a), x.begin() + static_cast<std::ptrdiff_t>(b));
            segs_nulo.emplace_back(yr.begin() + static_cast<std::ptrdiff_t>(a), yr.begin() + static_cast<std::ptrdiff_t>(b));
        }
        banda_ = derivar_banda(segs, W_MAX, cfg.banda, &segs_nulo);
        if (!banda_.viable) return;
        ConfiguracionSensorELW cs{};
        cs.W = banda_.W;
        cs.m = banda_.m;
        cs.S = cfg.paso_sensor;
        cs.orden = cfg.orden;
        cs.d_min = cfg.banda.d_min;
        cs.d_max = cfg.banda.d_max;
        sen_.configurar(cs);
        // 2. H^ real y de los nulos (sustitutos del crudo)
        // (el mismo generador sigue: nulo_sensor y nulos_contraste)
        const auto sr = serie(x, ok, fase, dia);
        real_ = memoria_detalle::estadisticos(sr, NULO_NINGUNO, cfg.banda.alfa);
        SerieMemoria sc;
        if (cfg.nulo_sensor == NULO_NINGUNO) sc = sr;
        else {
            const auto y = memoria_detalle::sustituto_crudo(x, ok, cfg.orden, cfg.nulo_sensor, banda_.W, rng);
            sc = serie(y, ok, fase, dia);
        }
        calibracion_ = memoria_detalle::estadisticos(sc, cfg.nulo_sensor, cfg.banda.alfa);
        for (std::uint32_t nb : NULOS_TODOS) {
            if (!(cfg.nulos_contraste & nb)) continue;
            const auto y = memoria_detalle::sustituto_crudo(x, ok, cfg.orden, nb, banda_.W, rng);
            auto f = memoria_detalle::estadisticos(serie(y, ok, fase, dia), nb, cfg.banda.alfa);
            if (f.validos > 1 && real_.validos > 1) {
                f.dif_real_menos_nulo = real_.mu - f.mu;
                const double vr = real_.sd * real_.sd * (std::isfinite(real_.tau_int) ? real_.tau_int : 1.0) / static_cast<double>(real_.validos);
                const double vn = f.sd * f.sd * (std::isfinite(f.tau_int) ? f.tau_int : 1.0) / static_cast<double>(f.validos);
                f.se_dif = std::sqrt(vr + vn);
            }
            contraste_.push_back(f);
        }
        // 3. operador delta y temperatura sobre H^ (FLUJO, MEDIA)
        ConfiguracionObservador cd = cfg.delta;
        cd.integracion = Integracion::FLUJO;
        cd.agregacion = AgregacionEscala::MEDIA;
        delta_.calibrar(cd, sc.h, {}, sc.ok, sc.fase, sc.dia);
        sen_.reiniciar();
        delta_.reiniciar_flujo();
        lista_ = true;
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

    // ------------------------------------------------------------------------------------------------------------------------
    // Entrada por el motor de ventanas de la base auditada (2026-11-04, orden del Comandante): las filas pasan por la aduana de
    // ark::ventanas, se agregan en barras de RELOJ INTRINSECO y la ventana del sensor es la del motor, con la W que el operador
    // deriva de su banda (entrada del motor) y su paso S. El anillo propio del sensor no se usa en este camino.
    // Validez para el estimador: la ventana no contiene ninguna barra con filas imputadas (el estimador no ve datos inventados). Los
    // huecos de calendario no invalidan: en reloj intrinseco son ausencia de nu, no un defecto (declarado).
    // La fuente es la comun de la base (ark/ventanas/flujo_barras.hpp): aduana y reloj, delta_nu (NaN = nu medio por fila),
    // agregacion y fase de la purga desde el calendario (si fase.activa, fase y dia salen del cierre de cada barra en hora local
    // del mercado y bins_fase del operador delta se fija a 2 bins_dia; si no, los entrega el llamador por fila).
    using FuenteFilas = ark::ventanas::FuenteFilasCxx<W_MAX>;
    struct SerieBarras { std::vector<double> x; std::vector<std::uint8_t> ok; std::vector<std::size_t> fase; std::vector<std::int64_t> dia; };

private:
    FuenteFilas fuente_{};
    double delta_nu_ = NaN;
    bool por_filas_ = false;
    ark::ventanas::MotorVentanasIntrinsecoCxx<W_MAX> motor_{};
    ark::ventanas::VentanaIntrinsecaCxx<W_MAX> vent_{};

public:
    [[nodiscard]] double delta_nu() const noexcept { return delta_nu_; }
    [[nodiscard]] const ark::ventanas::MotorVentanasIntrinsecoCxx<W_MAX>& motor() const noexcept { return motor_; }

    // Serie de barras intrinsecas de una historia de filas (calibracion; reserva). fase y dia: los de la fila que cierra cada barra.
    SerieBarras serie_barras(const FuenteFilas& fuente, double delta_nu, std::span<const ark::ventanas::EntradaPasoCxx> filas,
                             std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        SerieBarras s;
        auto fl = std::make_unique<ark::ventanas::FlujoBarrasIntrinsecasCxx<W_MAX>>();   // el mismo motor y la misma aduana, barra a barra
        if (fl->configurar(fuente, delta_nu) != ark::ventanas::ErrorMotorIntrinseco::NINGUNO) return s;
        for (std::size_t i = 0; i < filas.size(); ++i) {
            auto b = fl->paso(filas[i]);
            if (!b) continue;
            s.x.push_back(b->valor);
            s.ok.push_back(b->imputada ? 0 : 1);
            s.fase.push_back(fuente.fase.activa ? b->fase : (fase.empty() ? 0 : fase[i]));
            s.dia.push_back(fuente.fase.activa ? b->dia : (dia.empty() ? 0 : dia[i]));
        }
        return s;
    }

    // Calibra con una historia de FILAS: barras intrinsecas -> banda (W, m) -> el motor queda con la W del operador y su paso S.
    void calibrar_filas(const ConfiguracionMemoria& cfg, const FuenteFilas& fuente, std::span<const ark::ventanas::EntradaPasoCxx> filas,
                        std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        fuente_ = fuente;
        por_filas_ = false;
        delta_nu_ = ark::ventanas::resolver_delta_nu<W_MAX>(fuente, filas);
        if (!(delta_nu_ > 0.0)) throw std::invalid_argument("delta_nu: no se pudo derivar de la historia (sin filas con sesion)");
        const auto b = serie_barras(fuente, delta_nu_, filas, fase, dia);
        ConfiguracionMemoria c = cfg;
        if (fuente.fase.activa) c.delta.bins_fase = fuente.fase.bins_fase();
        calibrar(c, b.x, b.ok, b.fase, b.dia);
        if (!lista_) return;
        ark::ventanas::ConfiguracionMotorIntrinsecoCxx<W_MAX> cm{};
        cm.base = fuente.base;
        cm.W = banda_.W;            // la W la entrega el operador: la de su banda derivada
        cm.S = cfg.paso_sensor;
        cm.delta_nu = delta_nu_;
        cm.agregacion = fuente.agregacion;
        if (motor_.configurar(cm) != ark::ventanas::ErrorMotorIntrinseco::NINGUNO) { lista_ = false; return; }
        por_filas_ = true;
    }

    // Ciclo caliente por filas: aduana + barra intrinseca + ventana del motor -> ELW -> operador delta. Cero reservas.
    Salida paso_fila(const ark::ventanas::EntradaPasoCxx& e, std::size_t fase = 0, std::int64_t dia = 0) noexcept {
        Salida s{};
        if (!lista_ || !por_filas_) return s;
        const auto p = motor_.paso(e, vent_);
        if (p.resultado != ark::ventanas::ResultadoPaso::VENTANA_EMITIDA) return s;
        s.emitio = true;
        ResultadoELW r{};
        if (!vent_.hay_imputadas) r = sen_.estimar(vent_.valores.data());
        else r.se = elw::se_asintotico(banda_.m);
        s.sensor = r;
        if (fuente_.fase.activa) {   // fase de la purga: la del cierre de la ultima barra de la ventana, en hora local del mercado
            const auto fd = ark::dataset::fase_de(fuente_.fase, static_cast<std::int64_t>(vent_.t_cierre_ms));
            fase = fd.fase;
            dia = fd.dia;
        }
        s.delta = delta_.paso(r.valido ? r.h : NaN, NaN, r.valido, fase, dia);
        return s;
    }
};

}  // namespace ark::metrologia
