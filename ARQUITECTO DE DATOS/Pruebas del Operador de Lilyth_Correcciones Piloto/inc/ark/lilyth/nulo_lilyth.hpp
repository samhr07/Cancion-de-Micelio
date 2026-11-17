// nulo_lilyth.hpp -- Etapa 3 del Operador de Lilyth: poblaciones nulas (sin memoria y lineal de Onsager) por via (a) espontanea
// y via (b) de inyeccion generativa conectadas al ObservadorSaltos real, derivacion de T_max (E1.9), PIT de punto medio,
// contrastes del resumen (Delta_tau_iid y Delta_tau_Onsager) y verificacion de P3b. C++20.
//
// Autoria: Argon (Auditoria Externa), 2026-11-06, bajo el plano de Code (Opus 5.5) y las decisiones del Comandante (2026-11-06).
// Especificacion: PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md (E1.1-E1.22, seccion 8.2) y PLANO_ETAPAS_3_A_8.md (seccion 3).
// Prerregistro sellado de S1-S3: README_ETAPA3_LILYTH.md (SHA-256: b0c2e81d08e1fb2e01f5b6a1ca7a0ee66f4ad76078694bcb6ec360001a3a2e78).
#pragma once

#include "ark/metrologia/observador_saltos.hpp"
#include "ark/metrologia/whittle_exacto.hpp"
#include "ark/lilyth/modelo_lineal_nulo.hpp"
#include "ark/metrologia/calibrador_familia.hpp"
#include "ark/metrologia/escalera_diadica.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/sustitutos.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <random>
#include <span>
#include <thread>
#include <vector>

namespace ark::lilyth {

using ark::metrologia::NaN;

// Tope de memoria declarado para los intentos de inyeccion de la via (b) por evento (Regla 5).
inline constexpr std::size_t MAX_INTENTOS_VIA_B = 256;
// Tope de memoria para el perfil estacionario guardado de un evento al abrir ranura.
inline constexpr std::size_t MAX_PERFIL_EVENTO_LILYTH = 256;

// Ventana de exclusion de una inyeccion en la via (b) (seccion 3.3):
//   [t_iny - (w_ancla + guarda), t_iny + l + t_max_ventana]
struct VentanaInyeccion {
    std::int64_t t_iny = 0;
    std::size_t l = 0;
    std::int64_t ini = 0;
    std::int64_t fin = 0;
};

[[nodiscard]] constexpr VentanaInyeccion construir_ventana_inyeccion(
    std::int64_t t_iny, std::size_t l, std::size_t w_ancla, std::size_t guarda, std::size_t t_max_ventana) noexcept {
    const std::int64_t pre = static_cast<std::int64_t>(w_ancla + guarda);
    const std::int64_t post = static_cast<std::int64_t>(l + t_max_ventana);
    return VentanaInyeccion{t_iny, l, t_iny - pre, t_iny + post};
}

[[nodiscard]] constexpr bool ventanas_se_solapan(const VentanaInyeccion& a, const VentanaInyeccion& b) noexcept {
    return a.ini <= b.fin && b.ini <= a.fin;
}

// Evento nulo compacto capturado en una realizacion de (a) o (b), con la misma ficha que un evento real (seccion 3.1).
template <std::size_t T_MAX>
struct EventoNuloLilyth {
    int signo = 0;
    int escala_origen = 0;
    double g_disparo = NaN;
    double dy = NaN;
    double y_ancla = NaN;
    std::size_t w_ancla = 0;
    std::size_t guarda_ancla = 0;
    std::uint8_t desenlace = 0;             // static_cast<uint8_t>(DesenlaceRanura)
    std::size_t pasos_observados = 0;
    std::size_t tau_fin = 0;
    double tau_rest_primer_toque = NaN;
    double tau_rest_consolidado = NaN;
    double min_eta_desde_2 = NaN;
    std::array<double, T_MAX> eta{};
};

// Referencia nula de un evento real k frente a una poblacion nula (sin memoria o lineal).
struct ReferenciaNuloEvento {
    bool evaluado = false;
    bool nulo_por_inyeccion = false;        // true si uso la via (b)
    bool sin_nulo_b = false;                // true si p_det == 0 en las primeras m_min inyecciones o modelo invalido
    bool tope_memoria_b = false;            // true si el tope de intentos fue acotado por MAX_INTENTOS_VIA_B
    std::size_t n_nulos_a = 0;              // eventos espontaneos en (a) con misma escala, mismo signo y g >= g_k
    std::size_t n_nulos_b = 0;              // eventos detectados en (b)
    std::size_t n_nulos_usados = 0;         // tamano del conjunto nulo usado para PIT y mediana
    std::size_t intentos_b = 0;             // inyecciones intentadas en (b)
    std::size_t tope_intentos_b = 0;        // ceil(m_min / p_det_inicial) fijado en las primeras m_min inyecciones
    std::size_t no_detectados_b = 0;        // inyecciones donde Lilyth no abrio ranura en [t_iny, t_iny + l] con el signo de k
    std::size_t solapes_b = 0;              // inyecciones descartadas por solape de ventana
    std::size_t cortadas_b = 0;             // Code: detectadas y cortadas por una replica
    std::size_t no_medidas_b = 0;           // Code: detectadas sin medida
    double p_det_inicial = NaN;             // fraccion detectada en las primeras m_min inyecciones
    double tasa_no_detectados = NaN;        // no_detectados_b / intentos_b
    double tau_real = NaN;                  // tau_rest del evento real (NaN si censurado)
    bool real_censurado = false;
    double mediana_nulo = NaN;              // mediana del nulo (censurados en T_max codificados en T_max + 1)
    double q25_nulo = NaN, q75_nulo = NaN, iqr_nulo = NaN;   // recorrido intercuartilico al lado
    double min_nulo = NaN, max_nulo = NaN;                   // recorrido completo al lado
    double delta_tau = NaN;                 // tau_real_efectivo - mediana_nulo
    double u_pit = NaN;                     // PIT de punto medio en [0, 1]
    const char* motivo = "";
};

// Code (2026-11-15, orden del Comandante): rama del contraste del resumen, elegida por ConfiguracionLilyth::contraste_resumen.
//   POOL_LOO       : la de Argon (por omision, bit a bit): el PIT del real contra sus n nulos, comparado con un pool de PIT de cada nulo
//                    contra n - 1 vecinos. Con pools pequenos (m_min = 8) y tau empatadas es ANTICONSERVADORA: un real por encima de todo
//                    su nulo cae fuera del soporte del pool (PRUEBAS_JEFE/red/MEDIDO_E7B.md, sonda_e7b_censura).
//   INTERCAMBIABLE : el real es un miembro mas de su estrato {real} U {sus nulos emparejados}; bajo H0 es uno uniforme de ellos
//                    (PRUEBAS_JEFE/lilyth/bifurcacion_intercambiable/EXPECTATIVAS_BIFURCACION.md, sello e9cfa2c9...).
enum class ContrasteResumenLilyth : std::uint8_t { POOL_LOO = 0, INTERCAMBIABLE = 1 };

// Code (2026-11-16, orden del Comandante): tamano del pool de nulos emparejados por evento. M_MIN (por omision) = el de hoy, m_min(alfa).
// DERIVADO: con K eventos candidatos en la rama, objetivo n* = ceil((2 * 4 / alfa)^(1/K)) - 1 por evento (la condicion de resolucion del
// primer escalon del Holm de 4 de Lilyth, bilateral: prod (n_k + 1) >= 8 / alfa), sin ningun numero nuevo. Solo agrega nulos.
enum class PoolNuloLilyth : std::uint8_t { M_MIN = 0, DERIVADO = 1 };

// Code (2026-11-16, orden del Comandante: enmienda de E1.9). Horizonte T_max del nulo lineal:
//   CUANTIL_ALFA (por omision, E1.9 de Argon): cuantil 1 - alfa de tau_rest del nulo lineal; guarda de censura en alfa.
//   CUANTIL_RESOLUCION: cuantil 1 - alfa/8 y guarda en alfa/8. Por que: con E1.9 queda ~alfa de los nulos censurados en T_max + 1 y un
//   evento real censurado empata con todos ellos, asi que su p minimo es del orden de esa fraccion (c en una cola, 2c en dos) por grande que
//   sea el pool. Resolver el primer escalon del Holm de 4 de Lilyth (alfa/4) en la prueba bilateral pide c <= alfa/(2 * 4) = alfa/8.
//   2 = bilateral y 4 = la familia de Lilyth: los mismos del pool derivado; ningun numero nuevo.
enum class HorizonteLilyth : std::uint8_t { CUANTIL_ALFA = 0, CUANTIL_RESOLUCION = 1 };

// Nivel con que se deriva T_max (cuantil 1 - nivel y guarda de censura en nivel).
[[nodiscard]] inline double nivel_horizonte(HorizonteLilyth h, double alfa) noexcept {
    return h == HorizonteLilyth::CUANTIL_RESOLUCION ? alfa / (2.0 * 4.0) : alfa;
}

// Objetivo n* de PoolNuloLilyth::DERIVADO para K eventos candidatos (0 si K == 0 o alfa invalido).
[[nodiscard]] inline std::size_t objetivo_pool_derivado(std::size_t K, double alfa) noexcept {
    if (K == 0 || !(alfa > 0.0 && alfa < 1.0)) return 0;
    const double base = 2.0 * 4.0 / alfa;
    const double r = std::ceil(std::pow(base, 1.0 / static_cast<double>(K)) - 1e-12);
    return r > 1.0 ? static_cast<std::size_t>(r) - 1 : 0;
}

// Contraste global del resumen para una familia nula (Delta_tau_iid o Delta_tau_Onsager, seccion 3.5).
struct ContrasteResumenNulo {
    bool valido = false;
    std::size_t n_eventos = 0;
    std::size_t n_nulos_pool = 0;
    std::size_t b_mc = 0;
    double estadistico_s = NaN;             // media de (u_k - 0.5) sobre los eventos reales
    double ee_nulo_s = NaN;                 // desviacion tipica de S bajo el Monte Carlo nulo
    double p_bilateral = NaN;               // p bilateral por Monte Carlo sacando conjuntos del mismo tamano del nulo
    double mediana_delta_tau = NaN;         // mediana de los excesos por evento (tau_rest_k - mediana_nulo_k)
    double media_u_pit = NaN;               // media de u_k
    // Code (2026-11-15): solo los llena la rama INTERCAMBIABLE.
    bool intercambiable = false;
    std::size_t pool_min_estrato = 0;       // el menor numero de nulos de un estrato
    double p_min_alcanzable = NaN;          // p exacto de la configuracion mas extrema, con el piso 1 / (B + 1) del Monte Carlo
    bool resoluble = false;                 // p_min_alcanzable <= umbral (alfa / 4: primer escalon del Holm de 4 de E1.16)
};

// Verificacion P3b (seccion 3.5): eta_barra(tau) contra la ACF de la calma con banda de 2 ee por bootstrap de eventos.
template <std::size_t T_MAX>
struct CurvaOnsagerP3b {
    std::size_t tau_comparados = 0;
    std::size_t n_eventos = 0;
    std::array<double, T_MAX> eta_barra{};
    std::array<double, T_MAX> acf_calma{};
    std::array<double, T_MAX> ee_bootstrap{};
    double max_z_desvio = NaN;              // max_{tau} |eta_barra(tau) - acf_calma(tau)| / max(ee, 1e-6)
    bool dentro_de_banda_2ee = false;
};

// Ficha de calibracion de la Etapa 3 para una representacion (PURGADA o CRUDA).
template <std::size_t T_MAX>
struct FichaNuloRepresentacion {
    bool activa = false;
    bool calibrada = false;
    std::size_t m_min = 0;
    std::size_t B = 0;
    double b_pedido = NaN;
    bool b_topado = false;
    std::size_t t_max_derivado = 0;
    bool t_max_censurado = false;
    bool t_max_raiz_unitaria = false;
    bool nulo_fft_degradado_a_generativo = false;
    double fraccion_censurada_lineal = NaN;
    std::size_t sin_factor_purga = 0;          // Code (revision D7): barras nulas de la rama purgada con fase inmadura
    ModeloLinealNulo modelo_lineal{};
    std::vector<std::vector<double>> segmentos_calma{};
    std::vector<EventoNuloLilyth<T_MAX>> nulos_sin_memoria{};
    std::vector<EventoNuloLilyth<T_MAX>> nulos_lineal{};
    ContrasteResumenNulo contraste_iid{};
    ContrasteResumenNulo contraste_onsager{};
    CurvaOnsagerP3b<T_MAX> p3b{};
    // Code (2026-11-15): con ContrasteResumenLilyth::INTERCAMBIABLE, la rama vieja al lado (A/B en la misma corrida); vacios con POOL_LOO.
    ContrasteResumenNulo contraste_iid_pool_loo{};
    ContrasteResumenNulo contraste_onsager_pool_loo{};
    // Code (2026-11-16): fuente de los segmentos del nulo (mascara `fuente_nulo` de `calibrar`) y pool derivado.
    bool fuente_externa = false;               // true si `calibrar` recibio una mascara de fuente
    std::size_t n_fuente = 0;                  // muestras de la representacion que entraron a los segmentos del nulo
    bool pool_derivado = false;                // PoolNuloLilyth::DERIVADO
    std::size_t k_candidatos_pool = 0;         // K con que se fijo el objetivo (eventos RESTITUIDA o CENSURADA_HORIZONTE)
    std::size_t objetivo_pool = 0;             // n* = objetivo_pool_derivado(K, alfa); 0 con M_MIN
    double nivel_t_max = NaN;                  // Code (2026-11-16): nivel con que se derivo T_max (alfa con E1.9; alfa/8 con CUANTIL_RESOLUCION)
};

// Resultado de un lote de inyecciones en una misma realizacion (usado por la via (b) y por el control S2).
template <std::size_t T_MAX>
struct ResultadoInyeccionesRealizacion {
    std::size_t candidatas = 0;
    std::size_t solapes = 0;
    std::size_t detectadas = 0;
    std::size_t no_detectadas = 0;
    std::size_t cortadas = 0;              // Code: detectadas cuya ranura corto una replica (fuera, como los reales cortados)
    std::size_t no_medidas = 0;            // Code: detectadas sin medida (NO_MEDIDA_* o abiertas al final)
    std::size_t solapes_sin_guarda = 0;    // Code: solo con aplicar_guarda = false (control S2)
    std::size_t sin_factor_purga = 0;      // Code: barras de la rama purgada con fase inmadura (quedan sin purgar)
    std::size_t con_ajena_previa = 0;      // Code: inyecciones con una apertura en la historia del ancla (contaminan: solape)
    std::size_t con_reapertura_posterior = 0;   // Code: inyecciones con aperturas despues del perfil (no contaminan la medida)
    std::vector<EventoNuloLilyth<T_MAX>> eventos_admitidos{};
    std::vector<std::size_t> indice_inyeccion{};   // Code: de que inyeccion (indice en tiempos_iny) es cada evento admitido
};

namespace detalle_nulo {

[[nodiscard]] inline std::uint64_t splitmix64(std::uint64_t x) noexcept {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Code (revision D2): destino de una inyeccion detectada (su ranura propia). Admitida si su rama tiene desenlace medible (RESTITUIDA o
// CENSURADA_HORIZONTE); cortada si la corto una replica (3, fuera, como los reales cortados); si no, sin medida.
enum class DestinoInyeccion : std::uint8_t { ADMITIDA = 0, CORTADA = 1, NO_MEDIDA = 2 };
[[nodiscard]] constexpr DestinoInyeccion destino_inyeccion(bool extraible, std::uint8_t desenlace) noexcept {
    if (extraible) return DestinoInyeccion::ADMITIDA;
    if (desenlace == 3u) return DestinoInyeccion::CORTADA;
    return DestinoInyeccion::NO_MEDIDA;
}

// Code (revision D7, 2026-11-07): espejo INVERTIBLE de ark::metrologia::PurgaDiurna. Misma aritmetica, en el mismo orden, mas
// `factor(fase, dia)`: el sqrt(vm / vb) que `paso` aplicaria a la proxima barra, sin mutar el estado. Sirve para construir una serie
// nula en el espacio purgado y llevarla al espacio crudo (y = p / factor) de modo que la purga del operador delta la devuelva a p.
// Control C-D7a: paso() = PurgaDiurna::paso() bit a bit sobre la historia de calibracion y su continuacion.
template <std::size_t B_MAX>
class EspejoPurgaNulo {
    std::size_t bins_ = 1;
    std::size_t dias_minimos_ = 20;
    std::array<double, B_MAX> s_{}, sh_{};
    std::array<std::int64_t, B_MAX> n_{}, nh_{}, dias_b_{};
    std::int64_t dia_ = -1;
    std::size_t dias_ = 0;
    double vm_ = NaN;

    // estado comprometido tras un cambio de dia, sin tocar *this (para factor())
    void comprometido(std::array<double, B_MAX>& s, std::array<std::int64_t, B_MAX>& n, std::array<std::int64_t, B_MAX>& db,
                      double& vm) const noexcept {
        double acc = 0.0;
        std::size_t nb = 0;
        for (std::size_t b = 0; b < bins_; ++b) {
            if (nh_[b] > 0) ++db[b];
            s[b] += sh_[b]; n[b] += nh_[b];
            if (n[b] > 0) { acc += s[b] / static_cast<double>(n[b]); ++nb; }
        }
        vm = nb ? acc / static_cast<double>(nb) : NaN;
    }

public:
    EspejoPurgaNulo() = default;
    EspejoPurgaNulo(std::size_t bins, std::size_t dias_minimos)
        : bins_(bins < B_MAX ? (bins ? bins : 1) : B_MAX), dias_minimos_(dias_minimos) {}

    double paso(double y, std::size_t fase, std::int64_t dia, bool& valido) noexcept {
        if (dia != dia_) {
            if (dia_ >= 0) {
                double vm = NaN;
                comprometido(s_, n_, dias_b_, vm);
                for (std::size_t b = 0; b < bins_; ++b) { sh_[b] = 0.0; nh_[b] = 0; }
                vm_ = vm;
                ++dias_;
            }
            dia_ = dia;
        }
        const std::size_t b = fase % bins_;
        valido = false;
        double out = NaN;
        if (dias_b_[b] >= static_cast<std::int64_t>(dias_minimos_) && n_[b] > 0) {
            const double vb = s_[b] / static_cast<double>(n_[b]);
            if (vb > 0.0 && std::isfinite(y)) { out = y * std::sqrt(vm_ / vb); valido = true; }
        }
        if (std::isfinite(y)) { sh_[b] += y * y; ++nh_[b]; }
        return out;
    }

    // Factor que aplicaria `paso` a la proxima barra (fase, dia); NaN si la fase aun no esta madura.
    [[nodiscard]] double factor(std::size_t fase, std::int64_t dia) const noexcept {
        auto s = s_; auto n = n_; auto db = dias_b_;
        double vm = vm_;
        if (dia != dia_ && dia_ >= 0) comprometido(s, n, db, vm);
        const std::size_t b = fase % bins_;
        if (!(db[b] >= static_cast<std::int64_t>(dias_minimos_) && n[b] > 0)) return NaN;
        const double vb = s[b] / static_cast<double>(n[b]);
        return vb > 0.0 ? std::sqrt(vm / vb) : NaN;
    }
    [[nodiscard]] std::int64_t ultimo_dia() const noexcept { return dia_; }
};

// Code (revision D7): contexto con el que las realizaciones nulas alimentan al operador delta. Sin purga: fase 0 y dia 0 (como el
// valor por omision de ObservadorSaltos::paso). Con purga: la secuencia de fases de la calibracion, con los dias desplazados para que
// sigan a la calibracion (monotonos), y, en la rama PURGADA, la serie nula se construye en el espacio purgado y se lleva al crudo con el
// espejo invertible.
template <std::size_t B_MAX>
struct ContextoPurgaNulo {
    bool purga_activa = false;
    EspejoPurgaNulo<B_MAX> espejo_calibrado{};       // estado tras la historia de calibracion
    std::vector<std::size_t> fases{};                 // fases de la calibracion (vacio sin purga)
    std::vector<std::int64_t> dias{};                 // dias de la calibracion
    std::int64_t desplazamiento_dias = 0;             // ultimo dia de calibracion - primero + 1

    [[nodiscard]] std::size_t fase(std::size_t i) const noexcept { return fases.empty() ? 0 : fases[i % fases.size()]; }
    [[nodiscard]] std::int64_t dia(std::size_t i, std::size_t vuelta) const noexcept {
        if (dias.empty()) return 0;
        return dias[i % dias.size()] + desplazamiento_dias * static_cast<std::int64_t>(vuelta + 1 + i / dias.size());
    }
};

// Lleva una serie nula del espacio de la representacion al espacio crudo que recibe el operador delta.
//   est: valores estacionarios (FLUJO: y; NIVEL: incrementos), en el espacio purgado si `en_purgado`, crudo si no.
//   Escribe x_out (FLUJO: y crudo; NIVEL: nivel crudo reintegrado desde x0), fase_out y dia_out (lo que se pasa a paso()).
//   inicio: indice de la secuencia de fases de calibracion donde empieza la realizacion (para no empezar siempre en la fase 0).
template <std::size_t B_MAX>
inline void llevar_a_crudo(ark::metrologia::Integracion integracion, bool en_purgado, const ContextoPurgaNulo<B_MAX>& ctx,
                           const std::vector<double>& est, double x0, std::size_t inicio, std::size_t vuelta,
                           std::vector<double>& x_out, std::vector<std::size_t>& fase_out, std::vector<std::int64_t>& dia_out,
                           std::size_t& sin_factor) {
    const std::size_t n = est.size();
    x_out.assign(n, 0.0); fase_out.assign(n, 0); dia_out.assign(n, 0);
    sin_factor = 0;
    EspejoPurgaNulo<B_MAX> esp = ctx.espejo_calibrado;
    double nivel = std::isfinite(x0) ? x0 : 0.0;
    for (std::size_t t = 0; t < n; ++t) {
        const std::size_t f = ctx.purga_activa ? ctx.fase(inicio + t) : 0;
        const std::int64_t d = ctx.purga_activa ? ctx.dia(inicio + t, vuelta) : 0;
        fase_out[t] = f; dia_out[t] = d;
        // En NIVEL la primera barra no tiene incremento (el operador delta y Lilyth la ven como NaN): no hay y que llevar.
        const bool hay_y = !(integracion == ark::metrologia::Integracion::NIVEL && t == 0);
        double y = hay_y ? est[t] : NaN;
        if (hay_y && ctx.purga_activa && en_purgado) {
            const double fa = esp.factor(f, d);
            if (std::isfinite(fa) && fa > 0.0) y = est[t] / fa;
            else ++sin_factor;   // fase inmadura: el valor queda sin purgar (se cuenta)
        }
        if (ctx.purga_activa) { bool v = false; (void)esp.paso(y, f, d, v); }
        if (integracion == ark::metrologia::Integracion::NIVEL) {
            if (hay_y) nivel += y;
            x_out[t] = nivel;
        } else {
            x_out[t] = y;
        }
    }
}

// Codifica tau_rest con los censurados en T_max como el mayor valor (T_max + 1, seccion 3.5).
[[nodiscard]] inline double codificar_tau(double tau_rest, std::size_t t_max_ef) noexcept {
    if (std::isfinite(tau_rest) && tau_rest >= 1.0 && tau_rest <= static_cast<double>(t_max_ef)) {
        return tau_rest;
    }
    return static_cast<double>(t_max_ef + 1);
}

// PIT de punto medio (seccion 3.5): (#{v < x} + 0.5 * #{v == x}) / m.
[[nodiscard]] inline double pit_punto_medio(double x, std::span<const double> muestra_nula) noexcept {
    if (muestra_nula.empty() || !std::isfinite(x)) return NaN;
    std::size_t menores = 0;
    std::size_t iguales = 0;
    for (double v : muestra_nula) {
        if (v < x - 1e-12) {
            ++menores;
        } else if (std::abs(v - x) <= 1e-12) {
            ++iguales;
        }
    }
    const double peso_empate = 0.5;
    const double num = static_cast<double>(menores) + peso_empate * static_cast<double>(iguales);
    return num / static_cast<double>(muestra_nula.size());
}

// Deriva T_max segun E1.9 y la seccion 3.4 del plano.
template <std::size_t T_MAX>
inline std::size_t derivar_t_max(
    ark::metrologia::Integracion integracion,
    std::size_t w_ancla_j0,
    const ModeloLinealNulo& m_lin,
    std::span<const EventoNuloLilyth<T_MAX>> nulos_lin,
    double alfa,
    bool& out_censurado,
    bool& out_raiz_unitaria,
    double& out_frac_cens) {
    out_censurado = false;
    out_raiz_unitaria = false;
    out_frac_cens = 0.0;

    // Caso raiz unitaria: radio_max >= 1 en la Etapa 2, modelo no estacionario, o NIVEL con incrementos sin memoria.
    const bool nivel_paseo_iid = (integracion == ark::metrologia::Integracion::NIVEL && m_lin.valido && m_lin.p == 0 && m_lin.d == 0.0);
    if (!m_lin.valido || m_lin.radio_max >= 1.0 || nivel_paseo_iid) {
        out_raiz_unitaria = true;
        return std::clamp<std::size_t>(w_ancla_j0, 2, T_MAX);
    }

    if (nulos_lin.empty()) {
        out_censurado = true;
        out_frac_cens = 1.0;
        return T_MAX;
    }

    std::size_t cens = 0;
    std::vector<double> taus;
    taus.reserve(nulos_lin.size());
    for (const auto& ev : nulos_lin) {
        if (std::isfinite(ev.tau_rest_primer_toque) && ev.tau_rest_primer_toque <= static_cast<double>(T_MAX)) {
            taus.push_back(ev.tau_rest_primer_toque);
        } else {
            ++cens;
            taus.push_back(static_cast<double>(T_MAX));
        }
    }
    out_frac_cens = static_cast<double>(cens) / static_cast<double>(nulos_lin.size());
    if (out_frac_cens > alfa) {
        out_censurado = true;
        return T_MAX;
    }

    const double q = ark::metrologia::est::cuantil(taus, 1.0 - alfa);
    if (!std::isfinite(q)) {
        out_censurado = true;
        return T_MAX;
    }
    return std::clamp<std::size_t>(static_cast<std::size_t>(std::ceil(q)), 2, T_MAX);
}

// Empareja los eventos nulos espontaneos de la via (a) para un evento (escala_origen, signo, g_disparo), omitiendo opcionalmente `idx_omitir`.
template <std::size_t T_MAX>
inline std::vector<const EventoNuloLilyth<T_MAX>*> filtrar_espontaneos_excedencia(
    int escala_real,
    int signo_real,
    double g_real,
    std::span<const EventoNuloLilyth<T_MAX>> pool,
    std::size_t idx_omitir = std::numeric_limits<std::size_t>::max()) {
    std::vector<const EventoNuloLilyth<T_MAX>*> out;
    for (std::size_t i = 0; i < pool.size(); ++i) {
        if (i == idx_omitir) continue;
        const auto& ev = pool[i];
        const bool misma_escala_y_signo = (ev.escala_origen == escala_real && ev.signo == signo_real);
        const bool cumple_excedencia = (ev.g_disparo >= g_real);
        if (misma_escala_y_signo && cumple_excedencia) {
            out.push_back(&ev);
        }
    }
    return out;
}

template <std::size_t T_MAX>
inline std::vector<double> construir_muestra_nula_emparejada(
    int escala_real,
    int signo_real,
    double g_real,
    std::span<const EventoNuloLilyth<T_MAX>> pool,
    std::size_t m_min,
    std::size_t t_max_ef,
    std::size_t idx_omitir,
    std::size_t& out_n_a) {
    auto esp = filtrar_espontaneos_excedencia<T_MAX>(escala_real, signo_real, g_real, pool, idx_omitir);
    out_n_a = esp.size();
    std::vector<double> taus;
    taus.reserve( std::max(esp.size(), m_min) );
    for (const auto* p : esp) {
        taus.push_back(codificar_tau(p->tau_rest_consolidado, t_max_ef));
    }
    // Code (revision D1b): la entrega completaba aqui la muestra con eventos de g < g_k, despues con eventos de cualquier escala y
    // signo, y si no habia ninguno con valores 1.0. Eso rompe la excedencia y fabrica el nulo. Ahora devuelve SOLO la excedencia; si no
    // llega a m_min, quien llama usa la via (b) o declara el evento sin nulo. (m_min se conserva en la firma para el llamador.)
    (void)m_min;
    return taus;
}

// Code (revision D1b): PIT de punto medio de cada evento nulo contra su excedencia SIN el mismo (misma escala, mismo signo, g >= g_i),
// solo si esa excedencia tiene al menos m_min eventos (la misma regla que los eventos reales en la via (a)). Es el pool del Monte Carlo
// de la prueba del resumen.
template <std::size_t T_MAX>
inline void pool_pit_loo(std::span<const EventoNuloLilyth<T_MAX>> pool, std::size_t m_min, std::size_t t_max_ef, std::vector<double>& u_pool) {
    for (std::size_t i = 0; i < pool.size(); ++i) {
        const auto& ev = pool[i];
        std::size_t n_a = 0;
        const auto taus_i = construir_muestra_nula_emparejada<T_MAX>(ev.escala_origen, ev.signo, ev.g_disparo, pool, m_min, t_max_ef, i, n_a);
        if (taus_i.size() < m_min) continue;
        u_pool.push_back(pit_punto_medio(codificar_tau(ev.tau_rest_consolidado, t_max_ef), taus_i));
    }
}

// Rellena las estadisticas de ReferenciaNuloEvento a partir de los valores codificados del nulo.
inline void llenar_estadisticas_referencia(
    ReferenciaNuloEvento& ref,
    double tau_real_raw,
    std::size_t t_max_ef,
    const std::vector<double>& taus_nulo_codificados) {
    if (taus_nulo_codificados.empty()) {
        ref.evaluado = false;
        return;
    }
    ref.evaluado = true;
    ref.n_nulos_usados = taus_nulo_codificados.size();
    ref.tau_real = tau_real_raw;
    ref.real_censurado = !(std::isfinite(tau_real_raw) && tau_real_raw <= static_cast<double>(t_max_ef));
    const double tau_real_cod = codificar_tau(tau_real_raw, t_max_ef);

    ref.mediana_nulo = ark::metrologia::est::cuantil(taus_nulo_codificados, 0.5);
    ref.q25_nulo = ark::metrologia::est::cuantil(taus_nulo_codificados, 0.25);
    ref.q75_nulo = ark::metrologia::est::cuantil(taus_nulo_codificados, 0.75);
    ref.iqr_nulo = ref.q75_nulo - ref.q25_nulo;
    const auto [it_min, it_max] = std::minmax_element(taus_nulo_codificados.begin(), taus_nulo_codificados.end());
    ref.min_nulo = *it_min;
    ref.max_nulo = *it_max;
    ref.delta_tau = tau_real_cod - ref.mediana_nulo;
    ref.u_pit = pit_punto_medio(tau_real_cod, taus_nulo_codificados);
}

// Prueba del resumen (seccion 3.5): S = media de (u_k - 0.5), distribucion nula por Monte Carlo sacando conjuntos del mismo tamano K
// de la poblacion nula (cada uno con su PIT de punto medio contra el resto del nulo).
inline ContrasteResumenNulo evaluar_contraste_resumen(
    std::span<const double> u_reales,
    std::span<const double> delta_taus_reales,
    std::span<const double> u_pool_nulo,
    std::size_t b_mc,
    std::uint64_t semilla) {
    ContrasteResumenNulo res{};
    if (u_reales.empty() || u_pool_nulo.size() < 2) return res;

    const std::size_t K = u_reales.size();
    const std::size_t M = u_pool_nulo.size();
    res.n_eventos = K;
    res.n_nulos_pool = M;
    res.b_mc = std::max<std::size_t>(99, b_mc);

    double suma_u = 0.0;
    for (double u : u_reales) suma_u += u;
    res.media_u_pit = suma_u / static_cast<double>(K);
    res.estadistico_s = res.media_u_pit - 0.5;
    res.mediana_delta_tau = ark::metrologia::est::cuantil(
        std::vector<double>(delta_taus_reales.begin(), delta_taus_reales.end()), 0.5);

    std::mt19937_64 rng(splitmix64(semilla ^ 0xD1B54A32D192ED03ULL));
    std::uniform_int_distribution<std::size_t> dist(0, M - 1);
    std::size_t extremos = 0;
    double suma_s2 = 0.0;
    const double abs_s = std::abs(res.estadistico_s);

    for (std::size_t b = 0; b < res.b_mc; ++b) {
        double sb = 0.0;
        for (std::size_t k = 0; k < K; ++k) {
            sb += (u_pool_nulo[dist(rng)] - 0.5);
        }
        sb /= static_cast<double>(K);
        suma_s2 += sb * sb;
        if (std::abs(sb) >= abs_s - 1e-15) ++extremos;
    }

    res.ee_nulo_s = std::sqrt(suma_s2 / static_cast<double>(res.b_mc));
    res.p_bilateral = static_cast<double>(1 + extremos) / static_cast<double>(1 + res.b_mc);
    res.valido = true;
    return res;
}

// Code (2026-11-15, orden del Comandante): prueba del resumen INTERCAMBIABLE (ContrasteResumenLilyth::INTERCAMBIABLE).
// estratos[k][0] = valor del evento real k; estratos[k][1..] = sus nulos emparejados, en el mismo codigo que su PIT (tau codificado en
// T_max + 1 si censurado, o el min de eta del sobrepaso). Bajo H0 el real es un miembro uniforme de su estrato: u de cada miembro = PIT de
// punto medio contra el RESTO del estrato (el u del real es el u_pit de su referencia); S = media(u_real - 0.5); el Monte Carlo sortea un
// miembro uniforme por estrato (estratos independientes: supuesto declarado). Un estrato con un miembro no finito o sin nulos no entra.
// p_min_alcanzable es exacto: solo la configuracion con todos los maximos (o todos los minimos) alcanza max |S|.
inline ContrasteResumenNulo evaluar_contraste_intercambiable(
    const std::vector<std::vector<double>>& estratos,
    std::size_t b_mc,
    std::uint64_t semilla,
    double umbral_resolucion) {
    ContrasteResumenNulo res{};
    res.intercambiable = true;
    std::vector<std::vector<double>> u;
    std::vector<double> dt, resto;
    std::size_t total_nulos = 0;
    std::size_t pool_min = std::numeric_limits<std::size_t>::max();
    for (const auto& e : estratos) {
        if (e.size() < 2) continue;
        bool finitos = true;
        for (double v : e) finitos = finitos && std::isfinite(v);
        if (!finitos) continue;
        std::vector<double> uk(e.size());
        for (std::size_t m = 0; m < e.size(); ++m) {
            resto.clear();
            for (std::size_t q = 0; q < e.size(); ++q) if (q != m) resto.push_back(e[q]);
            uk[m] = pit_punto_medio(e[m], resto);
        }
        dt.push_back(e[0] - ark::metrologia::est::cuantil(std::vector<double>(e.begin() + 1, e.end()), 0.5));
        total_nulos += e.size() - 1;
        pool_min = std::min(pool_min, e.size() - 1);
        u.push_back(std::move(uk));
    }
    if (u.empty()) return res;

    const std::size_t K = u.size();
    res.n_eventos = K;
    res.n_nulos_pool = total_nulos;
    res.pool_min_estrato = pool_min;
    res.b_mc = std::max<std::size_t>(99, b_mc);

    double suma_u = 0.0;
    for (const auto& uk : u) suma_u += uk[0];
    res.media_u_pit = suma_u / static_cast<double>(K);
    res.estadistico_s = res.media_u_pit - 0.5;
    res.mediana_delta_tau = ark::metrologia::est::cuantil(dt, 0.5);

    std::mt19937_64 rng(splitmix64(semilla ^ 0xD1B54A32D192ED03ULL));
    std::vector<std::uniform_int_distribution<std::size_t>> dist;
    dist.reserve(K);
    for (const auto& uk : u) dist.emplace_back(0, uk.size() - 1);
    std::size_t extremos = 0;
    double suma_s2 = 0.0;
    const double abs_s = std::abs(res.estadistico_s);
    for (std::size_t b = 0; b < res.b_mc; ++b) {
        double sb = 0.0;
        for (std::size_t k = 0; k < K; ++k) sb += u[k][dist[k](rng)] - 0.5;
        sb /= static_cast<double>(K);
        suma_s2 += sb * sb;
        if (std::abs(sb) >= abs_s - 1e-15) ++extremos;
    }
    res.ee_nulo_s = std::sqrt(suma_s2 / static_cast<double>(res.b_mc));
    res.p_bilateral = static_cast<double>(1 + extremos) / static_cast<double>(1 + res.b_mc);

    // p exacto de la configuracion mas extrema: P(|S_b| >= max |S|) = [cola superior] prod(#max / n) + [cola inferior] prod(#min / n).
    double a = 0.0, c = 0.0, pa = 1.0, pc = 1.0;
    for (const auto& uk : u) {
        const auto [it_mn, it_mx] = std::minmax_element(uk.begin(), uk.end());
        std::size_t n_mx = 0, n_mn = 0;
        for (double v : uk) { n_mx += v >= *it_mx - 1e-12; n_mn += v <= *it_mn + 1e-12; }
        a += *it_mx - 0.5;
        c += 0.5 - *it_mn;
        pa *= static_cast<double>(n_mx) / static_cast<double>(uk.size());
        pc *= static_cast<double>(n_mn) / static_cast<double>(uk.size());
    }
    a /= static_cast<double>(K);
    c /= static_cast<double>(K);
    const double s_max = std::max(a, c);
    double p_ext = 1.0;
    if (s_max > 1e-15) {
        p_ext = 0.0;
        if (a >= s_max - 1e-15) p_ext += pa;
        if (c >= s_max - 1e-15) p_ext += pc;
    }
    res.p_min_alcanzable = std::max(p_ext, 1.0 / static_cast<double>(1 + res.b_mc));
    res.resoluble = res.p_min_alcanzable <= umbral_resolucion;
    res.valido = true;
    return res;
}

// Evalua P3b: eta_barra(tau) del real contra la ACF de la calma (rho_tau en FLUJO; en NIVEL, respuesta acumulada del AR de incrementos),
// con banda de 2 ee por bootstrap de eventos.
template <std::size_t T_MAX>
inline CurvaOnsagerP3b<T_MAX> evaluar_curva_p3b(
    ark::metrologia::Integracion integracion,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    std::span<const std::array<double, T_MAX>> etas_eventos,
    std::span<const std::size_t> pasos_eventos,
    std::size_t tau_hasta,
    std::uint64_t semilla,
    std::size_t b_boot,
    double alfa) {
    CurvaOnsagerP3b<T_MAX> out{};
    out.eta_barra.fill(NaN);
    out.acf_calma.fill(NaN);
    out.ee_bootstrap.fill(NaN);
    if (etas_eventos.empty() || segs_calma.empty()) return out;

    const std::size_t K = etas_eventos.size();
    out.n_eventos = K;
    const std::size_t H = std::min<std::size_t>(tau_hasta, T_MAX);
    out.tau_comparados = H;

    if (integracion == ark::metrologia::Integracion::FLUJO) {
        for (std::size_t tau = 1; tau <= H; ++tau) {
            std::size_t pares = 0;
            out.acf_calma[tau - 1] = ark::metrologia::est::acf_agrupada(segs_calma, m_lin.mu, m_lin.sd, tau, pares);
        }
    } else {
        std::vector<double> psi(H + 1, 0.0);
        psi[0] = 1.0;
        double acum = 1.0;
        for (std::size_t tau = 1; tau <= H; ++tau) {
            double p_tau = 0.0;
            for (std::size_t i = 1; i <= m_lin.p && i <= tau; ++i) {
                p_tau += m_lin.phi[i - 1] * psi[tau - i];
            }
            psi[tau] = p_tau;
            acum += p_tau;
            out.acf_calma[tau - 1] = acum;
        }
    }

    std::mt19937_64 rng(splitmix64(semilla ^ 0xA5A5A5A55A5A5A5AULL));
    // Code (revision D6): B del bootstrap = el presupuesto de Monte Carlo de la configuracion (antes 200 escrito).
    const std::size_t B_boot = std::max<std::size_t>(2, b_boot);
    double max_z = 0.0;
    for (std::size_t tau = 1; tau <= H; ++tau) {
        std::vector<double> vals;
        vals.reserve(K);
        for (std::size_t e = 0; e < K; ++e) {
            if (tau <= pasos_eventos[e] && std::isfinite(etas_eventos[e][tau - 1])) {
                vals.push_back(etas_eventos[e][tau - 1]);
            }
        }
        if (vals.empty()) continue;
        double s = 0.0;
        for (double v : vals) s += v;
        const double media = s / static_cast<double>(vals.size());
        out.eta_barra[tau - 1] = media;

        double ee = 0.0;
        if (vals.size() >= 2) {
            std::uniform_int_distribution<std::size_t> dist(0, vals.size() - 1);
            std::vector<double> medias_b(B_boot, 0.0);
            double sb = 0.0;
            for (std::size_t b = 0; b < B_boot; ++b) {
                double m_b = 0.0;
                for (std::size_t j = 0; j < vals.size(); ++j) m_b += vals[dist(rng)];
                m_b /= static_cast<double>(vals.size());
                medias_b[b] = m_b;
                sb += m_b;
            }
            const double mb_prom = sb / static_cast<double>(B_boot);
            double ss = 0.0;
            for (double v : medias_b) ss += (v - mb_prom) * (v - mb_prom);
            ee = std::sqrt(ss / static_cast<double>(B_boot - 1));
        }
        // Code (revision D6): sin suelo escrito (antes max(ee, 0.25 / sqrt(n))). Con ee = 0 y la media sobre la ACF, z = 0.
        out.ee_bootstrap[tau - 1] = ee;

        if (std::isfinite(out.acf_calma[tau - 1])) {
            const double dif = std::abs(media - out.acf_calma[tau - 1]);
            const double z = ee > 0.0 ? dif / ee : (dif == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
            if (z > max_z) max_z = z;
        }
    }
    out.max_z_desvio = max_z;
    // Code (revision D6): la banda sale de alfa con Bonferroni sobre los H rezagos comparados (antes 2.5 escrito).
    const double z_banda = ark::metrologia::est::norm_ppf(1.0 - alfa / (2.0 * static_cast<double>(std::max<std::size_t>(H, 1))));
    out.dentro_de_banda_2ee = (max_z <= z_banda);
    return out;
}

// Convierte una RanuraLilyth cerrada en un EventoNuloLilyth compacto para la rama indicada.
template <std::size_t T_MAX, class Ranura>
inline bool extraer_evento_nulo(const Ranura& r, bool rama_purgada, EventoNuloLilyth<T_MAX>& out) noexcept {
    const auto& rama = rama_purgada ? r.purgada : r.cruda;
    const auto d = static_cast<std::uint8_t>(rama.desenlace);
    // Solo eventos con desenlace de restitucion resuelto en el horizonte: RESTITUIDA (2) o CENSURADA_HORIZONTE (4)
    if (d != 2 && d != 4) return false;
    out.signo = r.signo;
    out.escala_origen = r.escala_origen;
    out.g_disparo = r.g_disparo;
    out.dy = rama.dy;
    out.y_ancla = rama.y_ancla;
    out.w_ancla = r.w_ancla;
    out.guarda_ancla = r.guarda_ancla;
    out.desenlace = d;
    out.pasos_observados = rama.pasos_observados;
    out.tau_fin = rama.tau_fin;
    out.tau_rest_primer_toque = rama.tau_rest_primer_toque;
    out.tau_rest_consolidado = rama.tau_rest_consolidado;
    out.min_eta_desde_2 = rama.min_eta_desde_2;
    out.eta = rama.eta;
    return true;
}

// Genera la parte ESTACIONARIA (FLUJO: y; NIVEL: incrementos) de longitud n de una realizacion nula, en el espacio de la
// representacion: sin memoria por sustituto de los segmentos de calma, o lineal por ModeloLinealNulo::simular. Code (revision): la
// reintegracion y el paso al crudo los hace `llevar_a_crudo` (antes se reintegraba aqui sin purga ni fases).
inline void construir_estacionaria_nula(
    bool es_lineal,
    std::uint32_t codigo_sustituto,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    std::size_t n,
    std::mt19937_64& rng,
    std::vector<double>& est) {
    est.assign(n, 0.0);
    if (n == 0 || segs_calma.empty()) return;
    if (es_lineal && m_lin.valido) {
        std::vector<double> buf(n + m_lin.burn_in + m_lin.p + 16, 0.0);
        simular(m_lin, rng, nullptr, 0, 0, est.data(), n, buf.data(), buf.size());
        return;
    }
    // Sin memoria: un sustituto (PERMUTACION o SIGNOS) por vuelta sobre los segmentos; si la realizacion es mas larga que la historia
    // se saca un sustituto NUEVO en cada vuelta (antes se repetia el mismo).
    const auto mo = ark::metrologia::est::momentos(segs_calma);
    std::size_t pos = 0;
    while (pos < n) {
        const auto segs_sust = ark::metrologia::sustituto(segs_calma, codigo_sustituto, 1, mo.mu, rng);
        const std::size_t antes = pos;
        for (const auto& g : segs_sust) {
            for (double v : g) {
                if (pos >= n) break;
                est[pos++] = v;
            }
            if (pos >= n) break;
        }
        if (pos == antes) break;
    }
}

// Corre B realizaciones espontaneas de la via (a) en paralelo (seccion 3.2 y reglas 7 y N2/N3). Code (revision D1): sin relleno
// fabricado; la serie entra al operador delta en el espacio crudo y con las fases de la calibracion (D7).
template <std::size_t T_MAX, std::size_t B_PURGA, class Obs, class Lil, class ConfigLil>
inline std::vector<EventoNuloLilyth<T_MAX>> correr_poblacion_via_a(
    const Obs& delta_original,
    const ark::metrologia::ConfiguracionObservador& cfg_delta,
    const ConfigLil& cfg_lil,
    const Lil& lil_plantilla,
    bool es_lineal,
    bool rama_purgada,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    const ContextoPurgaNulo<B_PURGA>& ctx,
    double x0_nivel,
    std::size_t n_cal,
    std::size_t B,
    std::size_t hilos,
    std::uint64_t semilla_base,
    std::size_t& sin_factor_total) {
    sin_factor_total = 0;
    if (B == 0 || n_cal == 0 || segs_calma.empty()) return {};

    std::vector<std::vector<EventoNuloLilyth<T_MAX>>> por_realizacion(B);
    std::vector<std::size_t> sin_factor_b(B, 0);

    std::mt19937_64 rng_compartido(semilla_base);
    std::mutex mtx_compartido;
    Obs delta_compartido = delta_original;

    ark::metrologia::familia::en_paralelo(B, hilos, [&](std::size_t b) {
        const bool semilla_independiente_por_indice = true;
        std::uint64_t sem_b = splitmix64(semilla_base ^ (0x9E3779B97F4A7C15ULL * (b + 1)));
        if (!semilla_independiente_por_indice) {
            std::lock_guard<std::mutex> lk(mtx_compartido);
            sem_b = rng_compartido() ^ static_cast<std::uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id()));
        }
        std::mt19937_64 rng_b(sem_b);

        std::vector<double> est, x_sim;
        std::vector<std::size_t> fases;
        std::vector<std::int64_t> dias;
        construir_estacionaria_nula(es_lineal, cfg_lil.nulo_sin_memoria, segs_calma, m_lin, n_cal, rng_b, est);
        llevar_a_crudo<B_PURGA>(cfg_delta.integracion, rama_purgada, ctx, est, x0_nivel, 0, b, x_sim, fases, dias, sin_factor_b[b]);

        // Regla 7: cada realizacion parte de una copia NUEVA del original (la PurgaDiurna acumula estado).
        const bool copia_nueva_por_realizacion = true;
        Obs d_local = delta_original;
        if (!copia_nueva_por_realizacion) {
            std::lock_guard<std::mutex> lk(mtx_compartido);
            d_local = delta_compartido;
        }
        d_local.reiniciar_flujo();

        Lil lil_local = lil_plantilla;
        lil_local.reiniciar_flujo_conservando_purga();

        for (std::size_t t = 0; t < n_cal; ++t) {
            const auto s = d_local.paso(x_sim[t], NaN, true, fases[t], dias[t]);
            const auto o = lil_local.paso(s, x_sim[t], true, NaN, fases[t], dias[t]);
            for (std::size_t i = 0; i < o.ranuras_cerradas; ++i) {
                EventoNuloLilyth<T_MAX> ev{};
                if (extraer_evento_nulo<T_MAX>(lil_local.cerrada(i), rama_purgada, ev)) {
                    por_realizacion[b].push_back(ev);
                }
            }
        }
        if (!copia_nueva_por_realizacion) {
            std::lock_guard<std::mutex> lk(mtx_compartido);
            delta_compartido = d_local;
        }
    });

    std::vector<EventoNuloLilyth<T_MAX>> total;
    for (std::size_t b = 0; b < B; ++b) {
        total.insert(total.end(), por_realizacion[b].begin(), por_realizacion[b].end());
        sin_factor_total += sin_factor_b[b];
    }
    return total;
}

// Ejecuta una realizacion nula con una o varias inyecciones candidatas en los tiempos `tiempos_iny`, verificando la guarda de
// solapamiento de ventanas (seccion 3.3 y control S2) y exigiendo que el operador delta dispare.
// Code (revision 2026-11-07):
//   - `sin_memoria` (opcion (b) aprobada por el Comandante): la base es un sustituto sin memoria de la calma y el perfil REAL (en el
//     espacio de la representacion) se escribe tal cual en t_iny..t_iny+l-1. Que la serie no responda es la hipotesis del nulo sin
//     memoria (en el nulo lineal ese mismo gesto es el mutante P8). En el nulo lineal se inyectan las innovaciones del perfil, como antes.
//   - D2: una apertura posterior al perfil mientras la ranura de la propia inyeccion sigue viva es una replica de ESE choque, no un
//     evento ajeno: no cuenta como solape. Ajenas son las aperturas dentro de la ventana fuera de [t_iny, t_iny + l] que no son replicas
//     de la ranura propia.
//   - Una inyeccion detectada cuya ranura acaba cortada por una replica no es "no detectada": se cuenta en `cortadas` y queda fuera,
//     igual que los eventos reales cortados (evaluar_nulos_representacion solo usa RESTITUIDA y CENSURADA_HORIZONTE en los dos lados).
//   - D7: la serie entra en el espacio crudo con las fases de la calibracion (llevar_a_crudo); NIVEL parte de x0, no de 100.
//   - `aplicar_guarda` = false solo para el control S2 (mide lo que pasaria sin la guarda con el MISMO instrumento).
template <std::size_t T_MAX, std::size_t B_PURGA, class Obs, class Lil>
inline ResultadoInyeccionesRealizacion<T_MAX> ejecutar_inyecciones_en_realizacion(
    const Obs& delta_original,
    const Lil& lil_plantilla,
    ark::metrologia::Integracion integracion,
    bool rama_purgada,
    bool sin_memoria,
    std::uint32_t codigo_sustituto,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    const ContextoPurgaNulo<B_PURGA>& ctx,
    double x0_nivel,
    std::span<const double> perfil_est,
    std::span<const double> e_iny,
    int signo_real,
    std::size_t w_ancla,
    std::size_t guarda,
    std::size_t t_max_ef,
    std::span<const std::int64_t> tiempos_iny,
    std::size_t n_sim,
    std::mt19937_64& rng,
    bool aplicar_guarda = true) {
    ResultadoInyeccionesRealizacion<T_MAX> res{};
    res.candidatas = tiempos_iny.size();
    const std::size_t l = perfil_est.size();
    if (l == 0 || tiempos_iny.empty() || n_sim == 0) return res;
    if (!sin_memoria && (!m_lin.valido || e_iny.size() != l)) return res;
    if (sin_memoria && segs_calma.empty()) return res;

    // 1. Ventanas de exclusion [t_iny - (W_ancla + guarda), t_iny + l + T_MAX] y solapes entre inyecciones.
    std::vector<VentanaInyeccion> ventanas(tiempos_iny.size());
    std::vector<bool> solapada(tiempos_iny.size(), false);
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        ventanas[i] = construir_ventana_inyeccion(tiempos_iny[i], l, w_ancla, guarda, T_MAX);
    }
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        for (std::size_t j = i + 1; j < tiempos_iny.size(); ++j) {
            if (ventanas_se_solapan(ventanas[i], ventanas[j])) { solapada[i] = true; solapada[j] = true; }
        }
    }

    // 2. Parte estacionaria con la inyeccion.
    std::vector<double> x_est(n_sim, 0.0);
    if (sin_memoria) {
        construir_estacionaria_nula(false, codigo_sustituto, segs_calma, m_lin, n_sim, rng, x_est);
        for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
            if (tiempos_iny[i] < 0) continue;
            const std::size_t t0 = static_cast<std::size_t>(tiempos_iny[i]);
            for (std::size_t k = 0; k < l && (t0 + k) < n_sim; ++k) x_est[t0 + k] = perfil_est[k];
        }
    } else {
        const std::size_t burn = m_lin.burn_in;
        const std::size_t total = burn + n_sim;
        std::vector<double> eps(total, 0.0);
        if (!m_lin.residuos.empty()) {
            std::uniform_int_distribution<std::size_t> dist(0, m_lin.residuos.size() - 1);
            for (std::size_t t = 0; t < total; ++t) eps[t] = m_lin.residuos[dist(rng)] - m_lin.media_residuos;
        } else {
            std::normal_distribution<double> dist(0.0, std::sqrt(m_lin.sigma2_innov));
            for (std::size_t t = 0; t < total; ++t) eps[t] = dist(rng);
        }
        const bool inyectar_en_innovaciones = true;
        if (inyectar_en_innovaciones) {
            for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
                if (tiempos_iny[i] < 0) continue;
                const std::size_t t0 = static_cast<std::size_t>(tiempos_iny[i]);
                for (std::size_t k = 0; k < l && (t0 + k) < n_sim; ++k) eps[burn + t0 + k] = e_iny[k];
            }
        }
        std::vector<double> u(total, 0.0);
        for (std::size_t t = 0; t < total; ++t) {
            double acc = eps[t];
            const std::size_t i_max = std::min(m_lin.p, t);
            for (std::size_t i = 1; i <= i_max; ++i) acc += m_lin.phi[i - 1] * u[t - i];
            u[t] = acc;
        }
        if (m_lin.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE && m_lin.d != 0.0) {
            std::vector<double> pi_inv(total, 0.0);
            ark::metrologia::elw::coeficientes(-m_lin.d, total, pi_inv.data());
            for (std::size_t t = 0; t < n_sim; ++t) {
                const std::size_t idx = burn + t;
                double s = 0.0;
                for (std::size_t k = 0; k <= idx; ++k) s += pi_inv[k] * u[idx - k];
                x_est[t] = m_lin.mu + s;
            }
        } else {
            for (std::size_t t = 0; t < n_sim; ++t) x_est[t] = m_lin.mu + u[burn + t];
        }
        if (!inyectar_en_innovaciones) {
            // Mutante M08 (P8): superpone el perfil sobre la serie fija sin propagarlo por el filtro lineal
            for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
                if (tiempos_iny[i] < 0) continue;
                const std::size_t t0 = static_cast<std::size_t>(tiempos_iny[i]);
                for (std::size_t k = 0; k < l && (t0 + k) < n_sim; ++k) x_est[t0 + k] = perfil_est[k];
            }
        }
    }

    // 3. Al espacio crudo con las fases de la calibracion (D7).
    std::vector<double> x_sim;
    std::vector<std::size_t> fases;
    std::vector<std::int64_t> dias;
    std::size_t sin_factor = 0;
    std::size_t inicio = 0;
    if (ctx.purga_activa && !ctx.fases.empty()) {
        std::uniform_int_distribution<std::size_t> d_ini(0, ctx.fases.size() - 1);
        inicio = d_ini(rng);
    }
    llevar_a_crudo<B_PURGA>(integracion, rama_purgada, ctx, x_est, x0_nivel, inicio, 0, x_sim, fases, dias, sin_factor);
    res.sin_factor_purga = sin_factor;

    // 4. Copia limpia del operador delta y de ObservadorLilyth sobre la realizacion; se registran todas las ranuras.
    Obs d_local = delta_original;
    d_local.reiniciar_flujo();
    Lil lil_local = lil_plantilla;
    lil_local.reiniciar_flujo_conservando_purga();

    struct Cierre {
        std::int64_t t_apertura = -1;
        std::int64_t t_cierre = -1;
        int signo = 0;
        bool extraible = false;
        std::uint8_t desenlace = 0;
        EventoNuloLilyth<T_MAX> ev{};
    };
    std::vector<Cierre> cierres;
    std::vector<std::int64_t> aperturas_todas;
    for (std::size_t t = 0; t < n_sim; ++t) {
        const auto s = d_local.paso(x_sim[t], NaN, true, fases[t], dias[t]);
        const auto o = lil_local.paso(s, x_sim[t], true, NaN, fases[t], dias[t]);
        if (o.abrio_ranura) aperturas_todas.push_back(o.t);
        for (std::size_t c = 0; c < o.ranuras_cerradas; ++c) {
            const auto& rc = lil_local.cerrada(c);
            Cierre ci{};
            ci.t_apertura = rc.t_apertura;
            ci.t_cierre = o.t;
            ci.signo = rc.signo;
            const auto& rama = rama_purgada ? rc.purgada : rc.cruda;
            ci.desenlace = static_cast<std::uint8_t>(rama.desenlace);
            ci.extraible = extraer_evento_nulo<T_MAX>(rc, rama_purgada, ci.ev);
            cierres.push_back(ci);
        }
    }
    // Ranuras aun abiertas al final: siguen vivas hasta n_sim.
    for (const auto& r : lil_local.ranuras()) {
        if (!r.ocupada) continue;
        Cierre ci{};
        ci.t_apertura = r.t_apertura;
        ci.t_cierre = static_cast<std::int64_t>(n_sim);
        ci.signo = r.signo;
        cierres.push_back(ci);
    }

    // 5. Ranura propia de cada inyeccion: la primera abierta en [t_iny, t_iny + l] con el signo del evento real.
    // D2: solo una apertura ANTERIOR a la inyeccion (en [ini, t_iny), la historia del ancla) contamina la medida. Una apertura posterior
    // al perfil o bien llega con la ranura propia viva (es una replica que la corta: se cuenta en `cortadas`, o una continuacion, que no
    // abre) o bien llega con la ranura propia ya cerrada (su medida termino y no la toca; suele ser el detector realarmando sobre el
    // mismo choque inyectado). El solape entre inyecciones lo decide la geometria de sus ventanas (arriba).
    const bool solo_contamina_lo_previo = true;   // el mutante C05 de Code lo pone en false (la regla de la entrega)
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        const std::int64_t t0 = tiempos_iny[i];
        const std::int64_t t1 = t0 + static_cast<std::int64_t>(l);
        const Cierre* propia = nullptr;
        for (const auto& ci : cierres) {
            if (ci.t_apertura >= t0 && ci.t_apertura <= t1 && ci.signo == signo_real) {
                if (propia == nullptr || ci.t_apertura < propia->t_apertura) propia = &ci;
            }
        }
        // Eventos ajenos dentro de la ventana.
        const VentanaInyeccion& vi = ventanas[i];
        bool hay_previa = false, hay_posterior = false;
        for (std::int64_t ta : aperturas_todas) {
            if (ta >= t0 && ta <= t1) continue;
            if (ta < vi.ini || ta > vi.fin) continue;
            if (ta < t0) hay_previa = true; else hay_posterior = true;
        }
        res.con_ajena_previa += hay_previa;
        res.con_reapertura_posterior += hay_posterior;
        if (hay_previa || (hay_posterior && !solo_contamina_lo_previo)) solapada[i] = true;
        if (solapada[i] && aplicar_guarda) {
            ++res.solapes;
            continue;
        }
        if (solapada[i]) ++res.solapes_sin_guarda;
        const bool exige_disparo_detector = true;
        if (propia != nullptr) {
            const DestinoInyeccion dest = destino_inyeccion(propia->extraible, propia->desenlace);
            if (dest == DestinoInyeccion::ADMITIDA && propia->extraible) {
                ++res.detectadas;
                res.eventos_admitidos.push_back(propia->ev);
                res.indice_inyeccion.push_back(i);
            } else if (dest == DestinoInyeccion::CORTADA) {
                ++res.cortadas;          // detectada y cortada por una replica: fuera, igual que los reales cortados
            } else {
                ++res.no_medidas;        // detectada pero sin medida (ancla, signo, horizonte abierto al final)
            }
        } else if (!exige_disparo_detector) {
            // Mutante M03 (Argon): mide la restitucion desde t_iny + l - 1 sin exigir que el detector dispare
            const std::size_t ta = static_cast<std::size_t>(std::max<std::int64_t>(1, t0 + static_cast<std::int64_t>(l) - 1));
            if (ta + t_max_ef < n_sim && ta >= w_ancla + guarda) {
                std::vector<double> w(x_sim.begin() + static_cast<std::ptrdiff_t>(ta - w_ancla - guarda),
                                      x_sim.begin() + static_cast<std::ptrdiff_t>(ta - guarda));
                const double y_anc = ark::metrologia::est::cuantil(w, 0.5);
                const double dy = x_sim[ta] - y_anc;
                if (std::isfinite(dy) && dy * static_cast<double>(signo_real) > 0.0) {
                    EventoNuloLilyth<T_MAX> ev_forzado{};
                    ev_forzado.signo = signo_real;
                    ev_forzado.y_ancla = y_anc;
                    ev_forzado.dy = dy;
                    ev_forzado.eta.fill(NaN);
                    for (std::size_t tau = 1; tau <= t_max_ef && (ta + tau) < n_sim; ++tau) {
                        const double e = (x_sim[ta + tau] - y_anc) / dy;
                        ev_forzado.eta[tau - 1] = e;
                        ev_forzado.pasos_observados = tau;
                        if (e <= 0.5 && !std::isfinite(ev_forzado.tau_rest_primer_toque)) {
                            ev_forzado.tau_rest_primer_toque = static_cast<double>(tau);
                            ev_forzado.tau_rest_consolidado = static_cast<double>(tau);
                        }
                    }
                    ev_forzado.desenlace = std::isfinite(ev_forzado.tau_rest_consolidado) ? 2u : 4u;
                    ev_forzado.tau_fin = std::isfinite(ev_forzado.tau_rest_consolidado)
                                             ? static_cast<std::size_t>(ev_forzado.tau_rest_consolidado)
                                             : t_max_ef;
                    ++res.detectadas;
                    res.eventos_admitidos.push_back(ev_forzado);
                    res.indice_inyeccion.push_back(i);
                } else {
                    ++res.no_detectadas;
                }
            } else {
                ++res.no_detectadas;
            }
        } else {
            ++res.no_detectadas;
        }
    }
    return res;
}

// Ejecuta el bucle secuencial de la via (b) para un evento real k (seccion 3.3 y control S3).
template <std::size_t T_MAX, std::size_t B_PURGA, class Obs, class Lil>
inline void ejecutar_via_b_evento(
    const Obs& delta_original,
    const Lil& lil_plantilla,
    ark::metrologia::Integracion integracion,
    bool rama_purgada,
    bool sin_memoria,
    std::uint32_t codigo_sustituto,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    const ContextoPurgaNulo<B_PURGA>& ctx,
    double x0_nivel,
    std::span<const double> hist_prev_est,
    std::span<const double> perfil_est,
    int signo_real,
    std::size_t w_ancla,
    std::size_t guarda,
    std::size_t m_min,
    std::size_t t_max_ef,
    std::uint64_t semilla_evento,
    ReferenciaNuloEvento& ref,
    std::vector<EventoNuloLilyth<T_MAX>>& nulos_b_out,
    std::size_t& cont_solapes_global) {
    nulos_b_out.clear();
    ref.nulo_por_inyeccion = true;
    const std::size_t l = perfil_est.size();
    if ((!sin_memoria && !m_lin.valido) || (sin_memoria && segs_calma.empty()) || l == 0 || m_min == 0) {
        ref.sin_nulo_b = true;
        ref.motivo = sin_memoria ? "via (b) sin memoria: sin calma o perfil vacio" : "via (b): modelo lineal invalido o perfil vacio";
        return;
    }

    // Lineal: innovaciones del perfil (filtro inverso). Sin memoria: el perfil se escribe tal cual, no hacen falta.
    std::vector<double> e_iny(sin_memoria ? 0 : l, 0.0);
    if (!sin_memoria) innovaciones(m_lin, hist_prev_est.data(), hist_prev_est.size(), perfil_est.data(), l, e_iny.data());
    if (!sin_memoria && !std::isfinite(e_iny[0])) {
        ref.sin_nulo_b = true;
        ref.motivo = "via (b): innovaciones no finitas";
        return;
    }

    const std::size_t margen_pre = w_ancla + guarda + 8;
    const std::size_t margen_post = l + t_max_ef + 8;
    const std::size_t n_sim = margen_pre + margen_post + 16;
    std::mt19937_64 rng(splitmix64(semilla_evento ^ 0x6C62272E07BB0142ULL));
    std::uniform_int_distribution<std::int64_t> dist_t(
        static_cast<std::int64_t>(w_ancla + guarda + 4),
        static_cast<std::int64_t>(n_sim - margen_post));

    const bool fijar_p_det_en_m_min = true;
    const bool respetar_tope_parada = true;
    std::size_t tope_intentos = MAX_INTENTOS_VIA_B;
    bool p_det_fijado = false;

    for (std::size_t intento = 1; intento <= MAX_INTENTOS_VIA_B; ++intento) {
        const std::int64_t t_iny = dist_t(rng);
        const std::array<std::int64_t, 1> t_arr{t_iny};
        auto r_lote = ejecutar_inyecciones_en_realizacion<T_MAX, B_PURGA>(
            delta_original, lil_plantilla, integracion, rama_purgada, sin_memoria, codigo_sustituto, segs_calma, m_lin, ctx, x0_nivel,
            perfil_est, e_iny, signo_real, w_ancla, guarda, t_max_ef, t_arr, n_sim, rng);

        ++ref.intentos_b;
        ref.solapes_b += r_lote.solapes;
        ref.cortadas_b += r_lote.cortadas;
        ref.no_medidas_b += r_lote.no_medidas;
        cont_solapes_global += r_lote.solapes;
        ref.no_detectados_b += r_lote.no_detectadas;
        for (const auto& ev : r_lote.eventos_admitidos) {
            nulos_b_out.push_back(ev);
            ++ref.n_nulos_b;
        }

        if (intento >= m_min && (!p_det_fijado || !fijar_p_det_en_m_min)) {
            p_det_fijado = true;
            ref.p_det_inicial = static_cast<double>(ref.n_nulos_b) / static_cast<double>(intento);
            if (ref.n_nulos_b == 0) {
                ref.sin_nulo_b = true;
                ref.tope_intentos_b = intento;
                ref.motivo = "via (b): p_det = 0 en las primeras m_min inyecciones";
                if (respetar_tope_parada) break;
            } else {
                const double tope_calc = std::ceil(static_cast<double>(m_min) / ref.p_det_inicial);
                if (tope_calc > static_cast<double>(MAX_INTENTOS_VIA_B)) {
                    ref.tope_memoria_b = true;
                    tope_intentos = MAX_INTENTOS_VIA_B;
                } else {
                    tope_intentos = static_cast<std::size_t>(tope_calc);
                }
                ref.tope_intentos_b = tope_intentos;
            }
        }

        if (ref.n_nulos_b >= m_min) {
            if (!p_det_fijado) {
                ref.p_det_inicial = 1.0;
                ref.tope_intentos_b = m_min;
            }
            ref.motivo = "via (b): completado m_min nulos por inyeccion";
            break;
        }
        if (p_det_fijado && respetar_tope_parada && intento >= tope_intentos) {
            ref.motivo = "via (b): tope de intentos ceil(m_min / p_det) alcanzado";
            break;
        }
    }

    if (ref.intentos_b > 0) {
        ref.tasa_no_detectados = static_cast<double>(ref.no_detectados_b) / static_cast<double>(ref.intentos_b);
    }
}

}  // namespace detalle_nulo

}  // namespace ark::lilyth