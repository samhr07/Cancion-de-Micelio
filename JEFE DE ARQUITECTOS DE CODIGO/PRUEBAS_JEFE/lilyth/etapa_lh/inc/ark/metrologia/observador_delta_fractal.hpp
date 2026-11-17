// observador_delta_fractal.hpp -- un observable, dos lentes: el operador delta sobre el CRUDO (saltos de nivel y temperatura) y el
// operador delta sobre su MEMORIA H^ (lente ELW o DFA), con un solo presupuesto de falsas alarmas y una etiqueta conjunta. Jefe de
// Arquitectos de Codigo (Opus 5.5), 2026-11-05, especificacion aprobada por el Comandante ("reutilizando las ventanas deslizantes e
// independientes, nulos autocalibrados y el pipeline ya hecho"). Prerregistro: PRUEBAS_JEFE/delta_fractal/EXPECTATIVAS_DELTA_FRACTAL.md.
//
// Lentes (las dos reutilizadas, sin tocar su logica):
//   crudo    ObservadorSaltos (serie) / ObservadorSaltosFilas (filas), modo DISJUNTO o DESLIZANTE segun su configuracion.
//   memoria  Memoria = ObservadorMemoria (Whittle local exacto) u ObservadorMemoriaDFA (DFA-1 solapado), con su propio operador delta
//            sobre H^ (tambien DISJUNTO o DESLIZANTE). Nulos de cada lente = entradas, autocalibrados sobre su historia.
// Presupuesto (sin numeros a mano): una falsa alarma de la familia cada ARL0 barras, repartida por igual entre TODOS los canales de las
//   dos lentes (2 J_crudo + 2 J_memoria), cada uno con la misma tasa por barra. La lente de memoria cuenta en pasos de H^ (S barras):
//   ARL0_crudo = ARL0 (J_c + J_m) / J_c barras;  ARL0_memoria = ARL0 (J_c + J_m) / (J_m S) pasos de H^.
//   Se calibra cada lente, se leen sus J y se recalibran los umbrales con su parte (recalibrar_arl0: misma historia, otro ARL0).
// Etiqueta conjunta de cada evento de la lente de memoria: MEMORIA_CON_VOLATILIDAD si el canal de temperatura de la ESCALA NATIVA del
//   crudo (enmienda 1: la dispersion de las escalas gruesas cambia tambien con la memoria) abrio una alarma
//   dentro de las barras que sostienen la evidencia del evento de memoria ([ventana del primer H^ posterior al cambio estimado, ahora]);
//   MEMORIA_PURA si no. Razon (medido en real, MEDIDO_DELTA_FRACTAL.md): H^ se mueve con la volatilidad; asi un mismo fenomeno no cuenta
//   como dos testigos. Se recuerdan las ultimas K_TEMP aperturas (declarado).
// Ciclo caliente: cero reservas. La calibracion reserva.
#pragma once

#include "ark/metrologia/observador_memoria.hpp"
#include "ark/metrologia/observador_memoria_dfa.hpp"
#include "ark/metrologia/observador_saltos.hpp"
#include "ark/metrologia/observador_saltos_filas.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace ark::metrologia {

enum class EtiquetaMemoria : std::uint8_t { MEMORIA_PURA = 0, MEMORIA_CON_VOLATILIDAD = 1 };
inline constexpr std::array<const char*, 2> NOMBRE_ETIQUETA{"MEMORIA_PURA", "MEMORIA_CON_VOLATILIDAD"};

struct EventoConjunto {
    FichaEvento memoria{};                 // el evento de la lente de memoria (canal 0: cambio de H; canal 1: regimen de dispersion de H)
    EtiquetaMemoria etiqueta = EtiquetaMemoria::MEMORIA_PURA;
    std::int64_t t_barra = -1;             // barra en que se clasifico
    std::int64_t inicio_barra = -1;        // primera barra de la evidencia del evento de memoria
    std::int64_t barra_temperatura = -1;   // la apertura de temperatura del crudo que lo acompana (-1: ninguna)
};

struct RepartoARL0 {
    std::size_t canales_crudo = 0, canales_memoria = 0;
    double arl0_familia_barras = NaN;
    double arl0_crudo_barras = NaN;        // familia del crudo
    double arl0_memoria_pasos_h = NaN;     // familia de la memoria, en pasos de H^
    double arl0_canal_barras = NaN;        // ARL0 de cada canal (los dos lentes), en barras
};

template <class Destino, std::size_t A>
[[nodiscard]] inline Destino fuente_como(const ark::ventanas::FuenteFilasCxx<A>& f) noexcept {
    // FuenteFilasCxx<W> solo usa W en un static_assert: los campos son los mismos (declarado y comprobado aqui).
    static_assert(sizeof(ark::ventanas::FuenteFilasCxx<A>) == sizeof(Destino));
    static_assert(std::is_trivially_copyable_v<ark::ventanas::FuenteFilasCxx<A>> && std::is_trivially_copyable_v<Destino>);
    Destino g{};
    std::memcpy(static_cast<void*>(&g), static_cast<const void*>(&f), sizeof g);
    return g;
}

template <class Memoria, std::size_t J_MAX = 24, std::size_t L = 120, std::size_t N_REG = 120, std::size_t B_MAX = 1440,
          std::size_t EV_MAX = 8, std::size_t W_FLUJO = 256, std::size_t K_TEMP = 64>
class ObservadorDeltaFractal {
public:
    using Crudo = ObservadorSaltos<J_MAX, L, N_REG, B_MAX, EV_MAX>;
    using CrudoFilas = ObservadorSaltosFilas<J_MAX, L, N_REG, B_MAX, EV_MAX, W_FLUJO>;
    using ConfigMemoria = typename Memoria::Configuracion;
    struct Configuracion {
        double arl0_familia_barras = 0.0;   // OBLIGATORIO: una falsa alarma de las DOS lentes juntas cada ARL0 barras
        ConfiguracionObservador crudo{};    // su arl0 lo fija el reparto
        ConfigMemoria memoria{};            // su delta.arl0 lo fija el reparto
    };
    struct Salida {
        std::int64_t t = -1;                // barra (indice desde el reinicio)
        bool barra = false;                 // camino por filas: se cerro una barra
        typename Crudo::Salida crudo{};
        typename Memoria::Salida memoria{};
        std::array<EventoConjunto, EV_MAX> conjuntos{};
        std::size_t n_conjuntos = 0;
    };

private:
    Configuracion cfg_{};
    Crudo crudo_{};
    CrudoFilas crudof_{};
    Memoria mem_{};
    RepartoARL0 rep_{};
    bool por_filas_ = false, lista_ = false;
    std::int64_t t_ = -1;
    std::array<std::int64_t, K_TEMP> temp_{};   // barras de las ultimas aperturas de temperatura del crudo
    std::size_t nt_ = 0;

    [[nodiscard]] const Crudo& crudo_activo() const noexcept { return por_filas_ ? crudof_.observador() : crudo_; }

    void repartir() {
        const std::size_t jc = crudo_activo().escalas();
        const std::size_t jm = mem_.lista() ? mem_.delta().escalas() : 0;
        rep_ = RepartoARL0{};
        rep_.canales_crudo = 2 * jc;
        rep_.canales_memoria = 2 * jm;
        rep_.arl0_familia_barras = cfg_.arl0_familia_barras;
        const double tot = static_cast<double>(jc + jm);
        rep_.arl0_canal_barras = cfg_.arl0_familia_barras * 2.0 * tot;
        if (jc > 0) {
            rep_.arl0_crudo_barras = cfg_.arl0_familia_barras * tot / static_cast<double>(jc);
            if (por_filas_) crudof_.recalibrar_arl0(rep_.arl0_crudo_barras); else crudo_.recalibrar_arl0(rep_.arl0_crudo_barras);
        }
        if (jm > 0) {
            rep_.arl0_memoria_pasos_h = cfg_.arl0_familia_barras * tot / (static_cast<double>(jm) * static_cast<double>(mem_.paso_sensor()));
            mem_.recalibrar_arl0(rep_.arl0_memoria_pasos_h);
        }
        lista_ = jc + jm > 0;
        reiniciar_flujo();
    }

    // Registro de aperturas de temperatura del crudo y etiqueta de los eventos de memoria de este paso.
    void conjugar(Salida& s) noexcept {
        // solo la escala NATIVA (enmienda 1): la dispersion de las escalas gruesas tambien cambia con la memoria (la varianza de una suma
        // de 2^j barras va como 2^{2Hj}), asi que su temperatura no separa volatilidad de memoria
        if (s.crudo.niveles > 0 && s.crudo.temp[0].alarma.hay && !s.crudo.temp[0].confirmacion) temp_[nt_++ % K_TEMP] = t_;
        if (!s.memoria.emitio) return;
        const std::int64_t W = static_cast<std::int64_t>(mem_.ventana()), S = static_cast<std::int64_t>(mem_.paso_sensor());
        for (std::size_t k = 0; k < s.memoria.delta.n_eventos && s.n_conjuntos < EV_MAX; ++k) {
            EventoConjunto c{};
            c.memoria = s.memoria.delta.eventos[k];
            c.t_barra = t_;
            // el H^ de paso q es la ventana de barras [q S, q S + W - 1]: la evidencia empieza en la ventana del paso t_cambio
            const std::int64_t q = c.memoria.alarma.t_cambio > 0 ? c.memoria.alarma.t_cambio : 0;
            c.inicio_barra = q * S;
            (void)W;
            const std::size_t n = nt_ < K_TEMP ? nt_ : K_TEMP;
            for (std::size_t i = 0; i < n; ++i) {
                const std::int64_t b = temp_[(nt_ - 1 - i) % K_TEMP];   // de la mas reciente a la mas antigua
                if (b >= c.inicio_barra && b <= t_) { c.etiqueta = EtiquetaMemoria::MEMORIA_CON_VOLATILIDAD; c.barra_temperatura = b; break; }
            }
            s.conjuntos[s.n_conjuntos++] = c;
        }
    }

public:
    [[nodiscard]] bool lista() const noexcept { return lista_; }
    [[nodiscard]] const RepartoARL0& reparto() const noexcept { return rep_; }
    [[nodiscard]] const Crudo& crudo() const noexcept { return crudo_activo(); }
    [[nodiscard]] const Memoria& memoria() const noexcept { return mem_; }

    // Camino de serie (barras ya hechas). Fuera del ciclo caliente.
    void calibrar(const Configuracion& cfg, std::span<const double> x, std::span<const double> den = {}, std::span<const std::uint8_t> ok = {},
                  std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        if (!(cfg.arl0_familia_barras > 1.0)) throw std::invalid_argument("Configuracion.arl0_familia_barras es obligatorio");
        cfg_ = cfg;
        por_filas_ = false;
        auto cc = cfg.crudo;
        cc.arl0_pasos_nativos = cfg.arl0_familia_barras;   // provisional: el reparto lo recalibra
        crudo_.calibrar(cc, x, den, ok, fase, dia);
        auto cm = cfg.memoria;
        cm.delta.arl0_pasos_nativos = cfg.arl0_familia_barras;
        mem_.calibrar(cm, x, ok, fase, dia);
        repartir();
    }

    // Camino por filas: las dos lentes sobre la MISMA fuente (la misma aduana, el mismo reloj, el mismo delta_nu: las mismas barras).
    void calibrar_filas(const Configuracion& cfg, const ark::ventanas::FuenteFilasCxx<W_FLUJO>& fuente,
                        std::span<const ark::ventanas::EntradaPasoCxx> filas, std::span<const std::size_t> fase = {},
                        std::span<const std::int64_t> dia = {}) {
        if (!(cfg.arl0_familia_barras > 1.0)) throw std::invalid_argument("Configuracion.arl0_familia_barras es obligatorio");
        cfg_ = cfg;
        por_filas_ = true;
        auto cc = cfg.crudo;
        cc.arl0_pasos_nativos = cfg.arl0_familia_barras;
        crudof_.calibrar_filas(cc, fuente, filas, fase, dia);
        auto cm = cfg.memoria;
        cm.delta.arl0_pasos_nativos = cfg.arl0_familia_barras;
        const auto fm = fuente_como<typename Memoria::FuenteFilas>(fuente);
        mem_.calibrar_filas(cm, fm, filas, fase, dia);
        repartir();
    }

    // Vacia la memoria de las dos lentes (camino de serie; el de filas se reinicia volviendo a calibrar o copiando el calibrado).
    void reiniciar_flujo() noexcept {
        t_ = -1;
        nt_ = 0;
        if (!por_filas_) { crudo_.reiniciar_flujo(); mem_.reiniciar_flujo(); }
    }

    // Ciclo caliente, camino de serie: una barra. Cero reservas.
    Salida paso(double x, double den = NaN, bool valido = true, std::size_t fase = 0, std::int64_t dia = 0) noexcept {
        Salida s{};
        s.t = ++t_;
        s.barra = true;
        s.crudo = crudo_.paso(x, den, valido, fase, dia);
        s.memoria = mem_.paso(x, valido, fase, dia);
        conjugar(s);
        return s;
    }

    // Ciclo caliente, camino por filas: una fila. Cero reservas.
    Salida paso_fila(const ark::ventanas::EntradaPasoCxx& e, std::size_t fase = 0, std::int64_t dia = 0) noexcept {
        Salida s{};
        const auto c = crudof_.paso_fila(e, fase, dia);
        s.memoria = mem_.paso_fila(e, fase, dia);
        if (!c.barra_cerrada) { s.t = t_; return s; }
        s.t = ++t_;
        s.barra = true;
        s.crudo = c.delta;
        conjugar(s);
        return s;
    }
};

}  // namespace ark::metrologia
