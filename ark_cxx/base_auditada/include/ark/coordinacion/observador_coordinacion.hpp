// observador_coordinacion.hpp -- operador de coordinacion de un PAR: el sensor de coordinacion (theta^ cada S barras) vigilado por
// el operador delta, componente a componente, mas el operador delta sobre el CRUDO de cada activo para la etiqueta, con un solo
// presupuesto de falsas alarmas. Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-05, diseno acordado con el Comandante.
// Prerregistro: JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/coordinacion_sensor/EXPECTATIVAS_COORDINACION.md.
//
// Sin CUSUM, sin h de Siegmund, sin compuerta Q, sin regularizador ni constantes de EL ABISMO (el modulo heredado, NO AUDITADO, no se
// toca: lo usa solo el ARK 1.0). El detector es el ObservadorSaltos de la base (DISJUNTO o DESLIZANTE, segun su configuracion).
//
// Nulos (autocalibrados sobre la historia del par): cada componente calibra su operador sobre la serie de theta^ de un sustituto
//   DIRECCION y NORMALIZADA  SIGNOS de centro cero en los retornos base de X y de Y (independientes): conserva RV_j y su co-movimiento
//   MAGNITUD                 ROTACION de Y n/2 contra X: conserva la memoria de cada serie, rompe el emparejamiento
// Presupuesto: una falsa alarma de la familia cada ARL0 barras, repartida por igual entre TODOS los canales (2 J del crudo de X, 2 J
//   del crudo de Y y 2 J_c de cada componente; los de los componentes cuentan en pasos de S):
//   ARL0_o = ARL0 tot / (J_o S_o), tot = J_x + J_y + suma J_c; cada canal ARL0 2 tot barras.
// Etiqueta de cada evento de un componente: COORDINACION_CON_VOLATILIDAD si la temperatura de la ESCALA NATIVA del crudo de X o de Y
//   (la regla de la enmienda 1 del observador conjunto) abrio una alarma en las barras que sostienen la evidencia ([ventana del primer
//   theta^ posterior al cambio estimado, ahora]); COORDINACION_PURA si no. Se recuerdan las ultimas K_TEMP aperturas (declarado).
// Guarda de pares: la clase del par (EJECUTABLE, REZAGO, NO_VERIFICADO) va sellada en cada evento; fuera de EJECUTABLE no se lee como
//   transmision.
// Camino de serie (barras alineadas en reloj de pared). Ciclo caliente (paso): cero reservas. La calibracion reserva.
//
// Dos presupuestos (Configuracion.presupuesto):
//   BONFERRONI_CANALES (omision)  el reparto de arriba, sin cambios.
//   FAMILIA_EPISODIOS (fase 2 de la calibracion de h; reglas de Argon aprobadas por el Comandante, 2026-11-05; prerregistro
//     PRUEBAS_JEFE/calibracion_h/EXPECTATIVAS_FASE2.md). Los crudos solo etiquetan: cada uno con ARL0 = ARL0_familia, fuera del
//     presupuesto. Los componentes se agrupan por su nulo: S1 = DIRECCION y NORMALIZADA (SIGNOS), S2 = MAGNITUD (ROTACION), cada una
//     con ARL0_sub = 2 ARL0_familia en EPISODIOS (coalescencia: abre cuando algun componente de la subfamilia abre, cierra cuando todos
//     cerraron). Igual rareza A por canal dentro de la subfamilia (a_c = A / (2 J_c S) pasos). A* se MIDE: se corre la maquina de
//     estados de cada componente sobre B realizaciones del sustituto juntadas como tramos y se biseca log A. B = min(B_max,
//     ceil(m_min ARL0_sub / N)), m_min = ceil((z_{1-alfa/2}/ln 2)^2) (intervalo de factor 2). Si B queda topado, el objetivo es
//     arl0_max_certificable = B N / m_min y la ficha lo marca. No exporta rareza (fase 1).
#pragma once

#include "ark/coordinacion/sensor_coordinacion.hpp"
#include "ark/metrologia/calibrador_familia.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/observador_saltos.hpp"
#include "ark/ventanas/contratos_ventanas.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <mutex>
#include <span>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace ark::coordinacion {

enum class EtiquetaCoordinacion : std::uint8_t { COORDINACION_PURA = 0, COORDINACION_CON_VOLATILIDAD = 1 };
inline constexpr std::array<const char*, 2> NOMBRE_ETIQUETA_COORD{"COORDINACION_PURA", "COORDINACION_CON_VOLATILIDAD"};

struct EventoCoordinacion {
    std::size_t componente = 0;
    ark::metrologia::FichaEvento evento{};
    EtiquetaCoordinacion etiqueta = EtiquetaCoordinacion::COORDINACION_PURA;
    std::int64_t t_barra = -1;              // barra en que se clasifico
    std::int64_t inicio_barra = -1;         // primera barra de la evidencia
    std::int64_t barra_temperatura = -1;    // la apertura de temperatura que lo acompana (-1: ninguna)
    int activo_temperatura = -1;            // 0 X, 1 Y
    ark::ventanas::ClasePar clase_par = ark::ventanas::ClasePar::NO_VERIFICADO;
};

struct RepartoCoordinacion {
    std::size_t canales_crudo_x = 0, canales_crudo_y = 0, canales_coordinacion = 0;
    double arl0_familia_barras = ark::metrologia::NaN;
    double arl0_canal_barras = ark::metrologia::NaN;
    double arl0_crudo_x_barras = ark::metrologia::NaN, arl0_crudo_y_barras = ark::metrologia::NaN;
    std::vector<double> arl0_componente_pasos;   // NaN si el componente no se pudo calibrar
};

struct FichaComponente {
    std::size_t pasos = 0, validos_real = 0, validos_nulo = 0;
    double mu_real = ark::metrologia::NaN, sd_real = ark::metrologia::NaN, mu_nulo = ark::metrologia::NaN, sd_nulo = ark::metrologia::NaN;
};

// Lo generico vive en metrologia/calibrador_familia.hpp (lo comparten todos los operadores delta).
using PresupuestoCoordinacion = ark::metrologia::PresupuestoFamilia;
using ark::metrologia::FichaSubfamilia;
using ark::metrologia::episodios_coalescencia;

// Subfamilias por nulo (modo FAMILIA_EPISODIOS).
inline constexpr std::size_t SUB_SIGNOS = 0, SUB_ROTACION = 1;
inline constexpr std::array<const char*, 2> NOMBRE_SUBFAMILIA{"S1_SIGNOS", "S2_ROTACION"};

template <std::size_t W, std::size_t J_MAX = 24, std::size_t L = 120, std::size_t N_REG = 120, std::size_t B_MAX = 1440,
          std::size_t EV_MAX = 8, std::size_t EVC_MAX = 64, std::size_t K_TEMP = 64>
class ObservadorCoordinacion {
public:
    using Sensor = SensorCoordinacion<W>;
    using Delta = ark::metrologia::ObservadorSaltos<J_MAX, L, N_REG, B_MAX, EV_MAX>;
    struct Configuracion {
        double arl0_familia_barras = 0.0;                       // OBLIGATORIO
        ConfiguracionSensorCoordinacion sensor{};                // paso_sensor y k_vecinos OBLIGATORIOS
        ark::metrologia::ConfiguracionObservador crudo{};        // el de X y el de Y (su arl0 lo fija el reparto)
        ark::metrologia::ConfiguracionObservador delta{};        // el de cada componente (FLUJO y MEDIA se fuerzan; arl0 del reparto)
        ark::ventanas::ClasePar clase_par = ark::ventanas::ClasePar::NO_VERIFICADO;
        std::uint64_t semilla_nulos = 20261105;
        PresupuestoCoordinacion presupuesto = PresupuestoCoordinacion::BONFERRONI_CANALES;
        // solo FAMILIA_EPISODIOS:
        std::size_t realizaciones_max = 10;     // B_max (presupuesto de computo declarado)
        double alfa_certificacion = 0.05;       // alfa del proyecto, fija m_min
        std::size_t iteraciones_max = 12;       // tope de la biseccion
        std::size_t expansiones_max = 10;       // tope del acotamiento (factor 4 cada una)
        std::size_t hilos = 0;                  // 0: hardware_concurrency
        // RAMA (2026-11-14, paso 2): POR_CANAL parte cada subfamilia de nulo (SIGNOS, ROTACION) en sus dos canales (nivel, temperatura),
        // cada uno con presupuesto 2 ARL0_sub (Bonferroni sobre 2) y su rareza A*, calibradas en alterna. CONJUNTO = lo de siempre.
        ark::metrologia::DefinicionEpisodio definicion_episodio = ark::metrologia::DefinicionEpisodio::CONJUNTO;
        std::size_t rondas_alternas_max = 6;    // POR_CANAL: tope de rondas de la alterna (declarado)
    };
    struct Salida {
        std::int64_t t = -1;
        typename Delta::Salida crudo_x{}, crudo_y{};
        bool emitio = false;
        typename Sensor::Resultado theta{};
        std::array<EventoCoordinacion, EVC_MAX> eventos{};
        std::size_t n_eventos = 0;
        std::size_t eventos_perdidos = 0;                         // eventos que no cupieron en EVC_MAX
        // episodios de coalescencia por subfamilia (en los dos presupuestos; solo cambian en los pasos en que el sensor emite)
        std::array<bool, 2> episodio_abierto{};
        std::array<bool, 2> episodio_nuevo{};
        // por canal (las dos ramas), por subfamilia de nulo: su union es el episodio de arriba
        std::array<bool, 2> episodio_nivel_abierto{}, episodio_nivel_nuevo{};
        std::array<bool, 2> episodio_temperatura_abierto{}, episodio_temperatura_nuevo{};
    };

private:
    Configuracion cfg_{};
    Delta cx_{}, cy_{};
    Sensor sen_{};
    std::vector<std::unique_ptr<Delta>> comp_;
    std::vector<std::uint8_t> comp_ok_;                           // el componente quedo calibrado
    std::vector<FichaComponente> fichas_;
    RepartoCoordinacion rep_{};
    // series nulas guardadas (para auditar la identidad componente = operador solo)
    std::vector<std::vector<double>> nulo_theta_;
    std::vector<std::uint8_t> nulo_ok_flip_, nulo_ok_rot_;
    std::vector<std::size_t> nulo_fase_;
    std::vector<std::int64_t> nulo_dia_;
    bool lista_ = false;
    std::int64_t t_ = -1;
    std::array<std::int64_t, K_TEMP> temp_t_{};
    std::array<int, K_TEMP> temp_a_{};
    std::size_t nt_ = 0;
    // subfamilias y episodios
    std::vector<std::uint8_t> sub_;                               // subfamilia de cada componente
    std::array<FichaSubfamilia, 2> fsub_{};                       // CONJUNTO: la subfamilia; POR_CANAL: su canal de NIVEL
    std::array<FichaSubfamilia, 2> fsubt_{};                      // POR_CANAL: el canal de TEMPERATURA de cada subfamilia
    std::array<std::size_t, 2> rondas_{};
    std::array<bool, 2> convergio_{};
    std::vector<double> arl0t_comp_;                              // POR_CANAL: el arl0 de temperatura de cada componente (pasos)
    std::vector<std::size_t> real_ini_, real_fin_;               // tramos [ini, fin) del nulo juntado, en pasos de theta^
    std::array<bool, 2> epi_{}, epn_{}, ept_{};

    [[nodiscard]] bool es_magnitud(std::size_t c) const noexcept { return sen_.descriptor(c).observable == OBS_MAGNITUD; }

    [[nodiscard]] ark::metrologia::ParametrosFamilia parametros_familia() const noexcept {
        ark::metrologia::ParametrosFamilia p{};
        p.realizaciones_max = cfg_.realizaciones_max;
        p.alfa_certificacion = cfg_.alfa_certificacion;
        p.iteraciones_max = cfg_.iteraciones_max;
        p.expansiones_max = cfg_.expansiones_max;
        p.hilos = cfg_.hilos;
        return p;
    }

    template <class F>
    void en_paralelo(std::size_t n, F&& f) const { ark::metrologia::familia::en_paralelo(n, cfg_.hilos, std::forward<F>(f)); }

    // Corre la maquina de estados del componente c (con sus umbrales actuales) sobre el nulo juntado, cada realizacion desde
    // reiniciar_flujo, igual que el ciclo caliente (theta^ si el paso es valido, NaN si no). act[p] = evento abierto o emitido.
    void correr_nulo(std::size_t c, std::vector<std::uint8_t>& act, ark::metrologia::CanalEpisodio canal = ark::metrologia::CanalEpisodio::CONJUNTO) {
        auto& d = *comp_[c];
        const auto& th = nulo_theta_[c];
        const auto& ok = nulo_ok(c);
        act.assign(th.size(), 0);
        for (std::size_t r = 0; r < real_ini_.size(); ++r) {
            d.reiniciar_flujo();
            for (std::size_t p = real_ini_[r]; p < real_fin_[r]; ++p) {
                const bool v = ok[p] != 0;
                const auto s = d.paso(v ? th[p] : ark::metrologia::NaN, ark::metrologia::NaN, v, nulo_fase_[p], nulo_dia_[p]);
                act[p] = ark::metrologia::actividad_canal(s, canal) ? 1 : 0;
            }
        }
        d.reiniciar_flujo();
    }

    // POR_CANAL: episodios del canal `canal` de la subfamilia q, con la rareza A en ese canal y `otro` fija en el otro.
    std::size_t evaluar_canal(std::size_t q, ark::metrologia::CanalEpisodio canal, double A, double otro, std::vector<std::vector<std::uint8_t>>& act) {
        const double S = static_cast<double>(sen_.paso_sensor());
        std::vector<std::size_t> cs;
        for (std::size_t c = 0; c < comp_.size(); ++c) if (comp_ok_[c] && sub_[c] == q) cs.push_back(c);
        const bool nivel = canal == ark::metrologia::CanalEpisodio::NIVEL;
        en_paralelo(cs.size(), [&](std::size_t i) {
            const std::size_t c = cs[i];
            const double k = 2.0 * static_cast<double>(comp_[c]->escalas()) * S;
            comp_[c]->recalibrar_arl0((nivel ? A : otro) / k, (nivel ? otro : A) / k);
            correr_nulo(c, act[c], canal);
        });
        return contar_episodios(q, act);
    }
    void fijar_subfamilia_por_canal(std::size_t q, double An, double At) {
        const double S = static_cast<double>(sen_.paso_sensor());
        std::vector<std::size_t> cs;
        for (std::size_t c = 0; c < comp_.size(); ++c) if (comp_ok_[c] && sub_[c] == q) cs.push_back(c);
        en_paralelo(cs.size(), [&](std::size_t i) {
            const std::size_t c = cs[i];
            const double k = 2.0 * static_cast<double>(comp_[c]->escalas()) * S;
            comp_[c]->recalibrar_arl0(An / k, At / k);
            rep_.arl0_componente_pasos[c] = An / k;
            arl0t_comp_[c] = At / k;
        });
    }

    [[nodiscard]] std::size_t contar_episodios(std::size_t q, const std::vector<std::vector<std::uint8_t>>& act) const {
        std::vector<std::uint8_t> m(comp_.size(), 0);
        for (std::size_t c = 0; c < comp_.size(); ++c) m[c] = (comp_ok_[c] && sub_[c] == q) ? 1 : 0;
        return episodios_coalescencia(act, m, real_ini_, real_fin_);
    }

    // Episodios de la subfamilia q con rareza A por canal (barras): recalibra y corre sus componentes en paralelo.
    std::size_t evaluar(std::size_t q, double A, std::vector<std::vector<std::uint8_t>>& act) {
        const double S = static_cast<double>(sen_.paso_sensor());
        std::vector<std::size_t> cs;
        for (std::size_t c = 0; c < comp_.size(); ++c) if (comp_ok_[c] && sub_[c] == q) cs.push_back(c);
        en_paralelo(cs.size(), [&](std::size_t i) {
            const std::size_t c = cs[i];
            comp_[c]->recalibrar_arl0(A / (2.0 * static_cast<double>(comp_[c]->escalas()) * S));
            correr_nulo(c, act[c]);
        });
        return contar_episodios(q, act);
    }

    void fijar_subfamilia(std::size_t q, double A) {
        const double S = static_cast<double>(sen_.paso_sensor());
        std::vector<std::size_t> cs;
        for (std::size_t c = 0; c < comp_.size(); ++c) if (comp_ok_[c] && sub_[c] == q) cs.push_back(c);
        en_paralelo(cs.size(), [&](std::size_t i) {
            const std::size_t c = cs[i];
            const double a = A / (2.0 * static_cast<double>(comp_[c]->escalas()) * S);
            comp_[c]->recalibrar_arl0(a);
            rep_.arl0_componente_pasos[c] = a;
        });
    }

    // Busca A* de la subfamilia q: acotamiento por factor 4 desde el A de Bonferroni y biseccion en log A.
    void buscar_a_estrella(std::size_t q, std::vector<std::vector<std::uint8_t>>& act) {
        auto& f = fsub_[q];
        const double S = static_cast<double>(sen_.paso_sensor());
        std::size_t jmax = 0;
        for (std::size_t c = 0; c < comp_.size(); ++c) if (comp_ok_[c] && sub_[c] == q) jmax = std::max(jmax, comp_[c]->escalas());
        const double suelo = 2.0 * 2.0 * static_cast<double>(jmax) * S;   // a_c >= 2 en el componente de mas escalas (declarado)
        ark::metrologia::familia::buscar_a_estrella(f, suelo, parametros_familia(), [&](double A) { return evaluar(q, A, act); });
    }

    void calibrar_familia(std::span<const double> x, std::span<const double> y, std::span<const std::uint8_t> okx,
                          std::span<const std::uint8_t> oky, std::span<const std::size_t> fase, std::span<const std::int64_t> dia) {
        const auto& cfg = cfg_;
        if (cfg.realizaciones_max < 1) throw std::invalid_argument("Configuracion.realizaciones_max >= 1");
        if (!(cfg.alfa_certificacion > 0.0 && cfg.alfa_certificacion < 1.0)) throw std::invalid_argument("Configuracion.alfa_certificacion en (0, 1)");
        const std::size_t n = x.size();
        auto cc = cfg.crudo;
        cc.arl0_pasos_nativos = cfg.arl0_familia_barras;   // los crudos solo etiquetan
        cx_.calibrar(cc, x, {}, okx, fase, dia);
        cy_.calibrar(cc, y, {}, oky, fase, dia);
        sen_.configurar(cfg.sensor);
        std::vector<std::uint8_t> ok(n, 1), oy(n, 1), ox(n, 1);
        for (std::size_t i = 0; i < n; ++i) {
            ox[i] = (okx.empty() || okx[i]) && std::isfinite(x[i]) ? 1 : 0;
            oy[i] = (oky.empty() || oky[i]) && std::isfinite(y[i]) ? 1 : 0;
            ok[i] = ox[i] && oy[i] ? 1 : 0;
        }
        const std::size_t m_min = ark::metrologia::familia::m_min(cfg.alfa_certificacion);
        const bool por_canal = cfg.definicion_episodio == ark::metrologia::DefinicionEpisodio::POR_CANAL;
        const double arl0_sub = 2.0 * cfg.arl0_familia_barras;
        const double arl0_canal = 2.0 * arl0_sub;   // POR_CANAL: Bonferroni sobre los 2 canales de cada subfamilia
        double b_pedido = 0.0;
        const std::size_t B = ark::metrologia::familia::realizaciones(m_min, por_canal ? arl0_canal : arl0_sub, n, cfg.realizaciones_max, b_pedido);
        // realizaciones: SIGNOS (x y luego y, del mismo generador; la primera es la de BONFERRONI) y ROTACION Delta_b = floor(b n / (B+1))
        GeneradorU01 g(cfg.semilla_nulos);
        std::vector<std::vector<double>> xs(B), ys(B), yr(B);
        std::vector<std::vector<std::uint8_t>> oyr(B, std::vector<std::uint8_t>(n, 0));
        for (std::size_t b = 0; b < B; ++b) { xs[b] = nulo_signos_cero(x, g); ys[b] = nulo_signos_cero(y, g); }
        for (std::size_t b = 0; b < B; ++b) {
            const std::size_t d = ((b + 1) * n) / (B + 1);
            yr[b] = nulo_rotacion(y, d);
            const auto rm = nulo_rotacion(std::span<const std::uint8_t>(oy), d);
            for (std::size_t i = 0; i < n; ++i) oyr[b][i] = (ox[i] && rm[i]) ? 1 : 0;
        }
        // series del sensor: 0 real, 1..B SIGNOS, B+1..2B ROTACION (cada una con su sensor, en paralelo)
        std::vector<typename Sensor::Serie> ser(2 * B + 1);
        en_paralelo(2 * B + 1, [&](std::size_t i) {
            auto s = std::make_unique<Sensor>();
            s->configurar(cfg.sensor);
            if (i == 0) ser[i] = s->serie(x, y, ok);
            else if (i <= B) ser[i] = s->serie(xs[i - 1], ys[i - 1], ok);
            else ser[i] = s->serie(x, yr[i - B - 1], oyr[i - B - 1]);
        });
        const std::size_t Lr = ser[1].ok.size();
        for (std::size_t i = 1; i < ser.size(); ++i)
            if (ser[i].ok.size() != Lr) throw std::logic_error("observador de coordinacion: realizaciones de distinta longitud");
        const std::size_t P = B * Lr + (B - 1);
        real_ini_.assign(B, 0);
        real_fin_.assign(B, 0);
        nulo_ok_flip_.assign(P, 0);
        nulo_ok_rot_.assign(P, 0);
        nulo_fase_.assign(P, 0);
        nulo_dia_.assign(P, 0);
        for (std::size_t b = 0; b < B; ++b) {
            real_ini_[b] = b * (Lr + 1);
            real_fin_[b] = real_ini_[b] + Lr;
            for (std::size_t k = 0; k < Lr; ++k) {
                const std::size_t p = real_ini_[b] + k;
                nulo_ok_flip_[p] = ser[1 + b].ok[k];
                nulo_ok_rot_[p] = ser[1 + B + b].ok[k];
                const auto br = static_cast<std::size_t>(ser[1 + b].barra[k]);
                nulo_fase_[p] = fase.empty() ? 0 : fase[br];
                nulo_dia_[p] = dia.empty() ? 0 : dia[br];
            }
        }
        const std::size_t C = sen_.componentes();
        sub_.assign(C, 0);
        nulo_theta_.assign(C, {});
        comp_.clear();
        comp_ok_.assign(C, 0);
        fichas_.assign(C, FichaComponente{});
        for (std::size_t c = 0; c < C; ++c) {
            const bool mag = es_magnitud(c);
            sub_[c] = static_cast<std::uint8_t>(mag ? SUB_ROTACION : SUB_SIGNOS);
            auto& th = nulo_theta_[c];
            th.assign(P, ark::metrologia::NaN);
            for (std::size_t b = 0; b < B; ++b)
                for (std::size_t k = 0; k < Lr; ++k) th[real_ini_[b] + k] = (mag ? ser[1 + B + b] : ser[1 + b]).theta[c][k];
            fichas_[c] = momentos(ser[0].theta[c], ser[0].ok, true, fichas_[c]);
            fichas_[c] = momentos(th, mag ? nulo_ok_rot_ : nulo_ok_flip_, false, fichas_[c]);
            comp_.push_back(std::make_unique<Delta>());
        }
        auto cd = cfg.delta;
        cd.integracion = ark::metrologia::Integracion::FLUJO;
        cd.agregacion = ark::metrologia::AgregacionEscala::MEDIA;
        cd.arl0_pasos_nativos = cfg.arl0_familia_barras;   // provisional: la busqueda de A* lo recalibra
        en_paralelo(C, [&](std::size_t c) {
            comp_[c]->calibrar(cd, nulo_theta_[c], {}, nulo_ok(c), nulo_fase_, nulo_dia_);
            comp_ok_[c] = comp_[c]->escalas() > 0 ? 1 : 0;
        });
        // reparto: crudos fuera del presupuesto
        rep_ = RepartoCoordinacion{};
        rep_.canales_crudo_x = 2 * cx_.escalas();
        rep_.canales_crudo_y = 2 * cy_.escalas();
        rep_.arl0_familia_barras = cfg.arl0_familia_barras;
        rep_.arl0_crudo_x_barras = cx_.escalas() ? cfg.arl0_familia_barras : ark::metrologia::NaN;
        rep_.arl0_crudo_y_barras = cy_.escalas() ? cfg.arl0_familia_barras : ark::metrologia::NaN;
        rep_.arl0_componente_pasos.assign(C, ark::metrologia::NaN);
        arl0t_comp_.assign(C, ark::metrologia::NaN);
        rondas_ = {};
        convergio_ = {};
        fsubt_ = {};
        std::size_t jc = 0;
        for (std::size_t c = 0; c < C; ++c) if (comp_ok_[c]) jc += comp_[c]->escalas();
        rep_.canales_coordinacion = 2 * jc;
        // A* por subfamilia
        std::vector<std::vector<std::uint8_t>> act(C);
        if (por_canal) {
            const double S = static_cast<double>(sen_.paso_sensor());
            for (std::size_t q = 0; q < 2; ++q) {
                FichaSubfamilia f0{};
                std::size_t J = 0, jmax = 0;
                for (std::size_t c = 0; c < C; ++c)
                    if (comp_ok_[c] && sub_[c] == q) { ++f0.componentes; J += comp_[c]->escalas(); jmax = std::max(jmax, comp_[c]->escalas()); }
                ark::metrologia::familia::preparar_ficha(f0, B, m_min, b_pedido, cfg.realizaciones_max, n, Lr, arl0_canal);
                if (q == SUB_ROTACION) f0.separacion_rotacion_barras = static_cast<double>(n) / static_cast<double>(B + 1);
                f0.canales = J;
                fsub_[q] = fsubt_[q] = f0;
                if (f0.componentes == 0) continue;
                f0.presente = true;
                f0.a_bonferroni = arl0_canal * static_cast<double>(J);
                const double suelo = 2.0 * 2.0 * static_cast<double>(jmax) * S;
                double an = f0.a_bonferroni, at = f0.a_bonferroni, an_prev = ark::metrologia::NaN, at_prev = ark::metrologia::NaN;
                std::size_t r = 0;
                for (r = 1; r <= cfg.rondas_alternas_max; ++r) {
                    fsub_[q] = f0;
                    ark::metrologia::familia::buscar_a_estrella(fsub_[q], suelo, parametros_familia(), [&](double A) {
                        return evaluar_canal(q, ark::metrologia::CanalEpisodio::NIVEL, A, at, act); });
                    an = fsub_[q].a_estrella;
                    fsubt_[q] = f0;
                    ark::metrologia::familia::buscar_a_estrella(fsubt_[q], suelo, parametros_familia(), [&](double A) {
                        return evaluar_canal(q, ark::metrologia::CanalEpisodio::TEMPERATURA, A, an, act); });
                    at = fsubt_[q].a_estrella;
                    if (an == an_prev && at == at_prev) { convergio_[q] = true; break; }
                    an_prev = an;
                    at_prev = at;
                }
                rondas_[q] = std::min(r, cfg.rondas_alternas_max);
                fijar_subfamilia_por_canal(q, an, at);
            }
            lista_ = (cx_.escalas() + cy_.escalas() + jc) > 0;
            reiniciar_flujo();
            return;
        }
        for (std::size_t q = 0; q < 2; ++q) {
            auto& f = fsub_[q];
            f = FichaSubfamilia{};
            for (std::size_t c = 0; c < C; ++c)
                if (comp_ok_[c] && sub_[c] == q) { ++f.componentes; f.canales += 2 * comp_[c]->escalas(); }
            ark::metrologia::familia::preparar_ficha(f, B, m_min, b_pedido, cfg.realizaciones_max, n, Lr, arl0_sub);
            if (q == SUB_ROTACION) f.separacion_rotacion_barras = static_cast<double>(n) / static_cast<double>(B + 1);
            if (f.componentes == 0) continue;
            f.presente = true;
            f.a_bonferroni = arl0_sub * static_cast<double>(f.canales);
            buscar_a_estrella(q, act);
            fijar_subfamilia(q, f.a_estrella);
        }
        lista_ = (cx_.escalas() + cy_.escalas() + jc) > 0;
        reiniciar_flujo();
    }

    void repartir() {
        const std::size_t jx = cx_.escalas(), jy = cy_.escalas();
        std::size_t jc = 0;
        for (std::size_t c = 0; c < comp_.size(); ++c) if (comp_ok_[c]) jc += comp_[c]->escalas();
        rep_ = RepartoCoordinacion{};
        rep_.canales_crudo_x = 2 * jx;
        rep_.canales_crudo_y = 2 * jy;
        rep_.canales_coordinacion = 2 * jc;
        rep_.arl0_familia_barras = cfg_.arl0_familia_barras;
        const double tot = static_cast<double>(jx + jy + jc);
        rep_.arl0_canal_barras = cfg_.arl0_familia_barras * 2.0 * tot;
        if (jx > 0) { rep_.arl0_crudo_x_barras = cfg_.arl0_familia_barras * tot / static_cast<double>(jx); cx_.recalibrar_arl0(rep_.arl0_crudo_x_barras); }
        if (jy > 0) { rep_.arl0_crudo_y_barras = cfg_.arl0_familia_barras * tot / static_cast<double>(jy); cy_.recalibrar_arl0(rep_.arl0_crudo_y_barras); }
        rep_.arl0_componente_pasos.assign(comp_.size(), ark::metrologia::NaN);
        const double S = static_cast<double>(sen_.paso_sensor());
        for (std::size_t c = 0; c < comp_.size(); ++c) {
            if (!comp_ok_[c]) continue;
            const double a = cfg_.arl0_familia_barras * tot / (static_cast<double>(comp_[c]->escalas()) * S);
            rep_.arl0_componente_pasos[c] = a;
            comp_[c]->recalibrar_arl0(a);
        }
        lista_ = tot > 0.0;
        reiniciar_flujo();
    }

    static FichaComponente momentos(const std::vector<double>& v, const std::vector<std::uint8_t>& ok, bool real, FichaComponente f) {
        double s = 0.0, q = 0.0;
        std::size_t n = 0;
        for (std::size_t i = 0; i < v.size(); ++i) if (ok[i] && std::isfinite(v[i])) { s += v[i]; ++n; }
        const double mu = n ? s / static_cast<double>(n) : ark::metrologia::NaN;
        for (std::size_t i = 0; i < v.size(); ++i) if (ok[i] && std::isfinite(v[i])) q += (v[i] - mu) * (v[i] - mu);
        const double sd = n > 1 ? std::sqrt(q / static_cast<double>(n - 1)) : ark::metrologia::NaN;
        f.pasos = v.size();
        if (real) { f.validos_real = n; f.mu_real = mu; f.sd_real = sd; } else { f.validos_nulo = n; f.mu_nulo = mu; f.sd_nulo = sd; }
        return f;
    }

    void conjugar(Salida& s, std::size_t c, const typename Delta::Salida& d) noexcept {
        const std::int64_t S = static_cast<std::int64_t>(sen_.paso_sensor());
        for (std::size_t k = 0; k < d.n_eventos; ++k) {
            if (s.n_eventos >= EVC_MAX) { ++s.eventos_perdidos; continue; }
            EventoCoordinacion e{};
            e.componente = c;
            e.evento = d.eventos[k];
            e.t_barra = t_;
            e.clase_par = cfg_.clase_par;
            // el theta^ de paso q es la ventana de barras [q S, q S + W - 1]
            const std::int64_t q = e.evento.alarma.t_cambio > 0 ? e.evento.alarma.t_cambio : 0;
            e.inicio_barra = q * S;
            const std::size_t n = nt_ < K_TEMP ? nt_ : K_TEMP;
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t p = (nt_ - 1 - i) % K_TEMP;   // de la mas reciente a la mas antigua
                if (temp_t_[p] >= e.inicio_barra && temp_t_[p] <= t_) {
                    e.etiqueta = EtiquetaCoordinacion::COORDINACION_CON_VOLATILIDAD;
                    e.barra_temperatura = temp_t_[p];
                    e.activo_temperatura = temp_a_[p];
                    break;
                }
            }
            s.eventos[s.n_eventos++] = e;
        }
    }

public:
    [[nodiscard]] bool lista() const noexcept { return lista_; }
    [[nodiscard]] const Sensor& sensor() const noexcept { return sen_; }
    [[nodiscard]] const Delta& crudo_x() const noexcept { return cx_; }
    [[nodiscard]] const Delta& crudo_y() const noexcept { return cy_; }
    [[nodiscard]] std::size_t componentes() const noexcept { return comp_.size(); }
    [[nodiscard]] bool componente_calibrado(std::size_t c) const noexcept { return comp_ok_[c] != 0; }
    [[nodiscard]] const Delta& componente(std::size_t c) const noexcept { return *comp_[c]; }
    [[nodiscard]] const FichaComponente& ficha(std::size_t c) const noexcept { return fichas_[c]; }
    [[nodiscard]] const RepartoCoordinacion& reparto() const noexcept { return rep_; }
    [[nodiscard]] const Configuracion& configuracion() const noexcept { return cfg_; }
    // la serie nula con que se calibro el componente c (theta^, ok, fase y dia por paso)
    [[nodiscard]] const std::vector<double>& nulo_theta(std::size_t c) const noexcept { return nulo_theta_[c]; }
    [[nodiscard]] const std::vector<std::uint8_t>& nulo_ok(std::size_t c) const noexcept { return es_magnitud(c) ? nulo_ok_rot_ : nulo_ok_flip_; }
    [[nodiscard]] const std::vector<std::size_t>& nulo_fase() const noexcept { return nulo_fase_; }
    [[nodiscard]] const std::vector<std::int64_t>& nulo_dia() const noexcept { return nulo_dia_; }
    // modo FAMILIA_EPISODIOS: fichas de cada subfamilia y tramos [ini, fin) del nulo juntado (vacios en BONFERRONI_CANALES)
    [[nodiscard]] const FichaSubfamilia& ficha_subfamilia(std::size_t q) const noexcept { return fsub_[q]; }
    [[nodiscard]] const FichaSubfamilia& ficha_subfamilia_temperatura(std::size_t q) const noexcept { return fsubt_[q]; }
    [[nodiscard]] std::size_t rondas_alternas(std::size_t q) const noexcept { return rondas_[q]; }
    [[nodiscard]] bool alterna_convergio(std::size_t q) const noexcept { return convergio_[q]; }
    [[nodiscard]] double arl0_componente_temperatura(std::size_t c) const noexcept { return c < arl0t_comp_.size() ? arl0t_comp_[c] : ark::metrologia::NaN; }
    [[nodiscard]] std::size_t subfamilia(std::size_t c) const noexcept { return sub_[c]; }
    [[nodiscard]] const std::vector<std::size_t>& nulo_inicios() const noexcept { return real_ini_; }
    [[nodiscard]] const std::vector<std::size_t>& nulo_fines() const noexcept { return real_fin_; }

    // x, y: retornos base del par, alineados barra a barra. Fuera del ciclo caliente.
    void calibrar(const Configuracion& cfg, std::span<const double> x, std::span<const double> y, std::span<const std::uint8_t> okx = {},
                  std::span<const std::uint8_t> oky = {}, std::span<const std::size_t> fase = {}, std::span<const std::int64_t> dia = {}) {
        if (!(cfg.arl0_familia_barras > 1.0)) throw std::invalid_argument("Configuracion.arl0_familia_barras es obligatorio");
        if (x.size() != y.size()) throw std::invalid_argument("observador de coordinacion: x e y de distinta longitud");
        const std::size_t n = x.size();
        cfg_ = cfg;
        lista_ = false;
        fsub_ = {};
        real_ini_.clear();
        real_fin_.clear();
        if (cfg.presupuesto == PresupuestoCoordinacion::FAMILIA_EPISODIOS) { calibrar_familia(x, y, okx, oky, fase, dia); return; }
        auto cc = cfg.crudo;
        cc.arl0_pasos_nativos = cfg.arl0_familia_barras;   // provisional: el reparto lo recalibra
        cx_.calibrar(cc, x, {}, okx, fase, dia);
        cy_.calibrar(cc, y, {}, oky, fase, dia);
        sen_.configurar(cfg.sensor);
        std::vector<std::uint8_t> ok(n, 1);
        for (std::size_t i = 0; i < n; ++i)
            ok[i] = (okx.empty() || okx[i]) && (oky.empty() || oky[i]) && std::isfinite(x[i]) && std::isfinite(y[i]) ? 1 : 0;
        // nulos
        GeneradorU01 g(cfg.semilla_nulos);
        const auto xs = nulo_signos_cero(x, g);
        const auto ys = nulo_signos_cero(y, g);
        const auto yr = nulo_rotacion(y);
        std::vector<std::uint8_t> oyr(n, 1);
        {
            std::vector<std::uint8_t> oy(n, 1);
            for (std::size_t i = 0; i < n; ++i) oy[i] = (oky.empty() || oky[i]) && std::isfinite(y[i]) ? 1 : 0;
            const auto r = nulo_rotacion(std::span<const std::uint8_t>(oy));
            for (std::size_t i = 0; i < n; ++i) oyr[i] = ((okx.empty() || okx[i]) && std::isfinite(x[i]) && r[i]) ? 1 : 0;
        }
        const auto sr = sen_.serie(x, y, ok);
        const auto sf = sen_.serie(xs, ys, ok);
        const auto so = sen_.serie(x, yr, oyr);
        const std::size_t C = sen_.componentes();
        nulo_fase_.assign(sf.barra.size(), 0);
        nulo_dia_.assign(sf.barra.size(), 0);
        for (std::size_t p = 0; p < sf.barra.size(); ++p) {
            const auto b = static_cast<std::size_t>(sf.barra[p]);
            nulo_fase_[p] = fase.empty() ? 0 : fase[b];
            nulo_dia_[p] = dia.empty() ? 0 : dia[b];
        }
        nulo_ok_flip_ = sf.ok;
        nulo_ok_rot_ = so.ok;
        nulo_theta_.assign(C, {});
        comp_.clear();
        comp_ok_.assign(C, 0);
        fichas_.assign(C, FichaComponente{});
        auto cd = cfg.delta;
        cd.integracion = ark::metrologia::Integracion::FLUJO;
        cd.agregacion = ark::metrologia::AgregacionEscala::MEDIA;
        cd.arl0_pasos_nativos = cfg.arl0_familia_barras;
        for (std::size_t c = 0; c < C; ++c) {
            const bool mag = es_magnitud(c);
            nulo_theta_[c] = mag ? so.theta[c] : sf.theta[c];
            fichas_[c] = momentos(sr.theta[c], sr.ok, true, fichas_[c]);
            fichas_[c] = momentos(nulo_theta_[c], mag ? so.ok : sf.ok, false, fichas_[c]);
            comp_.push_back(std::make_unique<Delta>());
            comp_[c]->calibrar(cd, nulo_theta_[c], {}, mag ? so.ok : sf.ok, nulo_fase_, nulo_dia_);
            comp_ok_[c] = comp_[c]->escalas() > 0 ? 1 : 0;
        }
        sub_.assign(C, 0);
        for (std::size_t c = 0; c < C; ++c) sub_[c] = static_cast<std::uint8_t>(es_magnitud(c) ? SUB_ROTACION : SUB_SIGNOS);
        repartir();
    }

    // Solo BONFERRONI_CANALES: en FAMILIA_EPISODIOS el A* depende de B, que depende del ARL0: se calibra de nuevo.
    void recalibrar_arl0(double arl0_familia_barras) {
        if (cfg_.presupuesto == PresupuestoCoordinacion::FAMILIA_EPISODIOS)
            throw std::logic_error("recalibrar_arl0: en FAMILIA_EPISODIOS se llama a calibrar de nuevo");
        if (!(arl0_familia_barras > 1.0)) throw std::invalid_argument("recalibrar_arl0: ARL0 > 1");
        cfg_.arl0_familia_barras = arl0_familia_barras;
        repartir();
    }

    void reiniciar_flujo() noexcept {
        t_ = -1;
        nt_ = 0;
        cx_.reiniciar_flujo();
        cy_.reiniciar_flujo();
        sen_.reiniciar();
        for (std::size_t c = 0; c < comp_.size(); ++c) if (comp_ok_[c]) comp_[c]->reiniciar_flujo();
        epi_ = {};
        epn_ = {};
        ept_ = {};
    }

    // Ciclo caliente: una barra del par. Cero reservas.
    Salida paso(double x, double y, bool valido_x = true, bool valido_y = true, std::size_t fase = 0, std::int64_t dia = 0) noexcept {
        Salida s{};
        s.t = ++t_;
        if (!lista_) return s;
        s.crudo_x = cx_.paso(x, ark::metrologia::NaN, valido_x, fase, dia);
        s.crudo_y = cy_.paso(y, ark::metrologia::NaN, valido_y, fase, dia);
        // solo la escala NATIVA (enmienda 1 del observador conjunto)
        if (s.crudo_x.niveles > 0 && s.crudo_x.temp[0].alarma.hay && !s.crudo_x.temp[0].confirmacion) { temp_t_[nt_ % K_TEMP] = t_; temp_a_[nt_ % K_TEMP] = 0; ++nt_; }
        if (s.crudo_y.niveles > 0 && s.crudo_y.temp[0].alarma.hay && !s.crudo_y.temp[0].confirmacion) { temp_t_[nt_ % K_TEMP] = t_; temp_a_[nt_ % K_TEMP] = 1; ++nt_; }
        const auto r = sen_.paso(x, y, valido_x && valido_y);
        s.episodio_abierto = epi_;
        s.episodio_nivel_abierto = epn_;
        s.episodio_temperatura_abierto = ept_;
        if (!r) return s;
        s.emitio = true;
        s.theta = *r;
        std::array<bool, 2> u{}, un{}, ut{};
        for (std::size_t c = 0; c < comp_.size(); ++c) {
            if (!comp_ok_[c]) continue;
            const auto d = comp_[c]->paso(r->valido ? r->theta[c] : ark::metrologia::NaN, ark::metrologia::NaN, r->valido, fase, dia);
            if (d.evento_abierto || d.n_eventos > 0) u[sub_[c]] = true;
            if (ark::metrologia::actividad_canal(d, ark::metrologia::CanalEpisodio::NIVEL)) un[sub_[c]] = true;
            if (ark::metrologia::actividad_canal(d, ark::metrologia::CanalEpisodio::TEMPERATURA)) ut[sub_[c]] = true;
            if (d.n_eventos) conjugar(s, c, d);
        }
        for (std::size_t q = 0; q < 2; ++q) {
            s.episodio_nuevo[q] = u[q] && !epi_[q]; epi_[q] = u[q];
            s.episodio_nivel_nuevo[q] = un[q] && !epn_[q]; epn_[q] = un[q];
            s.episodio_temperatura_nuevo[q] = ut[q] && !ept_[q]; ept_[q] = ut[q];
        }
        s.episodio_abierto = epi_;
        s.episodio_nivel_abierto = epn_;
        s.episodio_temperatura_abierto = ept_;
        return s;
    }
};

}  // namespace ark::coordinacion
