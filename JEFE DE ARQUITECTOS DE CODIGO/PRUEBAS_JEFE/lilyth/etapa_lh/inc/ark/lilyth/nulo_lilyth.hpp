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
#include "ark/lilyth/replicas_lilyth.hpp"
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
    // E1.24 (nulo L+H): biblioteca de marcas de las hijas, por tipo (0: signo +, 1: signo -). Cada entrada son las innovaciones
    // (filtro inverso A(L) del modelo lineal de ESTA representacion) del perfil de una alarma de la calibracion reproducida.
    std::array<std::vector<std::vector<double>>, 2> biblioteca_innov{};
    std::vector<double> respuesta_impulso{};   // E1.24: respuesta de la serie de la representacion a una innovacion unitaria
    bool lh_utilizable = false;                // replicas DECIDIBLE y biblioteca >= m_min en los dos tipos
    const char* motivo_lh = "";
    ContrasteResumenNulo contraste_lh{};
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
    std::size_t hijas = 0;                 // E1.24: hijas inyectadas (nulo L+H)
    std::size_t cumulos_desbordados = 0;   // E1.24: inyecciones cuyo cumulo paso el tope de memoria (sin medida)
};

namespace detalle_nulo {

// E1.24: el generador de hijas del nulo L+H. Lo que necesita la inyeccion: el Hawkes ajustado y la biblioteca de marcas por tipo.
struct GeneradorHijasNulo {
    const AjusteHawkesReplicas* hawkes = nullptr;
    const std::array<std::vector<std::vector<double>>, 2>* biblioteca = nullptr;
    const std::vector<double>* respuesta = nullptr;   // respuesta al impulso del modelo lineal (h[0] = 1)
};

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

// Paso a paso de `llevar_a_crudo` (E1.24): el simulador en linea del nulo L+H lo necesita porque las hijas se suman a la serie
// despues de que el detector vea la alarma que las dispara. `llevar_a_crudo` es el mismo bucle sobre este paso: misma aritmetica,
// mismo orden (identidad por construccion, control R8).
template <std::size_t B_MAX>
class PasoCrudo {
    ark::metrologia::Integracion integ_ = ark::metrologia::Integracion::FLUJO;
    bool en_purgado_ = false;
    const ContextoPurgaNulo<B_MAX>* ctx_ = nullptr;
    EspejoPurgaNulo<B_MAX> esp_{};
    double nivel_ = 0.0;
    std::size_t inicio_ = 0, vuelta_ = 0;

public:
    PasoCrudo(ark::metrologia::Integracion integracion, bool en_purgado, const ContextoPurgaNulo<B_MAX>& ctx, double x0, std::size_t inicio,
              std::size_t vuelta)
        : integ_(integracion), en_purgado_(en_purgado), ctx_(&ctx), esp_(ctx.espejo_calibrado), nivel_(std::isfinite(x0) ? x0 : 0.0),
          inicio_(inicio), vuelta_(vuelta) {}

    // Valor crudo de la barra t a partir del valor estacionario est_t; devuelve la fase y el dia que se pasan a paso().
    double paso(std::size_t t, double est_t, std::size_t& fase, std::int64_t& dia, bool& sin_factor) noexcept {
        const std::size_t f = ctx_->purga_activa ? ctx_->fase(inicio_ + t) : 0;
        const std::int64_t d = ctx_->purga_activa ? ctx_->dia(inicio_ + t, vuelta_) : 0;
        fase = f; dia = d; sin_factor = false;
        // En NIVEL la primera barra no tiene incremento (el operador delta y Lilyth la ven como NaN): no hay y que llevar.
        const bool hay_y = !(integ_ == ark::metrologia::Integracion::NIVEL && t == 0);
        double y = hay_y ? est_t : NaN;
        if (hay_y && ctx_->purga_activa && en_purgado_) {
            const double fa = esp_.factor(f, d);
            if (std::isfinite(fa) && fa > 0.0) y = est_t / fa;
            else sin_factor = true;   // fase inmadura: el valor queda sin purgar (se cuenta)
        }
        if (ctx_->purga_activa) { bool v = false; (void)esp_.paso(y, f, d, v); }
        if (integ_ == ark::metrologia::Integracion::NIVEL) {
            if (hay_y) nivel_ += y;
            return nivel_;
        }
        return y;
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
    PasoCrudo<B_MAX> pc(integracion, en_purgado, ctx, x0, inicio, vuelta);
    for (std::size_t t = 0; t < n; ++t) {
        bool sf = false;
        x_out[t] = pc.paso(t, est[t], fase_out[t], dia_out[t], sf);
        sin_factor += sf;
    }
}

// Respuesta de la serie de la representacion a una innovacion unitaria en t = 0 (E1.24): g del AR(p) por su recurrencia y, en
// ARFIMA, convolucionada con los pesos de (1 - L)^(-d), como `simular`. h[0] = 1.
inline std::vector<double> respuesta_impulso(const ModeloLinealNulo& m, std::size_t n) {
    std::vector<double> g(n, 0.0);
    if (!m.valido || n == 0) return g;
    g[0] = 1.0;
    for (std::size_t t = 1; t < n; ++t) {
        double acc = 0.0;
        for (std::size_t i = 1; i <= m.p && i <= t; ++i) acc += m.phi[i - 1] * g[t - i];
        g[t] = acc;
    }
    if (m.modo_usado != ModoNuloLineal::ARFIMA_WHITTLE || m.d == 0.0) return g;
    std::vector<double> pi(n, 0.0);
    ark::metrologia::elw::coeficientes(-m.d, n, pi.data());
    std::vector<double> h(n, 0.0);
    for (std::size_t t = 0; t < n; ++t) {
        double s = 0.0;
        for (std::size_t k = 0; k <= t; ++k) s += pi[k] * g[t - k];
        h[t] = s;
    }
    return h;
}

// Corre una copia del operador delta y de Lilyth sobre x_est (espacio de la representacion) llevado al crudo PASO A PASO. Con `gen`,
// cada alarma del espejo dispara hijas de primera generacion (E1.24): su marca (innovaciones de una alarma de calibracion de su tipo)
// entra a la serie desde t + retardo por la respuesta al impulso. `visitar(t, x_crudo, salida_lilyth, lilyth)` ve cada paso.
// Devuelve false si las hijas pasan el tope (la realizacion no es el nulo: desborde).
template <std::size_t B_PURGA, class Obs, class Lil, class Visitante>
inline bool correr_en_linea(Obs& d, Lil& lil, ark::metrologia::Integracion integracion, bool rama_purgada, const ContextoPurgaNulo<B_PURGA>& ctx,
                            double x0, std::size_t inicio, std::size_t vuelta, std::vector<double>& x_est, const GeneradorHijasNulo* gen,
                            std::mt19937_64& rng, std::size_t tope_hijas, std::size_t& hijas, std::size_t& sin_factor, Visitante&& visitar) {
    const std::size_t n = x_est.size();
    PasoCrudo<B_PURGA> pc(integracion, rama_purgada, ctx, x0, inicio, vuelta);
    std::vector<std::pair<std::int64_t, std::uint8_t>> hv;
    const bool con_hijas = gen != nullptr && gen->hawkes != nullptr && gen->biblioteca != nullptr && gen->respuesta != nullptr;
    bool desborde = false;
    for (std::size_t t = 0; t < n; ++t) {
        std::size_t f = 0;
        std::int64_t dd = 0;
        bool sf = false;
        const double x = pc.paso(t, x_est[t], f, dd, sf);
        sin_factor += sf;
        const auto s = d.paso(x, NaN, true, f, dd);
        const auto o = lil.paso(s, x, true, NaN, f, dd);
        visitar(t, x, o, lil);
        if (!con_hijas || !o.alarma || desborde) continue;
        detalle_replicas::hijas_directas(*gen->hawkes, static_cast<std::uint8_t>(o.signo_alarma > 0 ? 0 : 1), rng, hv);
        for (const auto& [ret, tipo] : hv) {
            const auto& bib = (*gen->biblioteca)[tipo & 1u];
            if (bib.empty()) continue;
            std::uniform_int_distribution<std::size_t> u(0, bib.size() - 1);
            const auto& e_h = bib[u(rng)];
            const std::size_t ini = t + static_cast<std::size_t>(ret);
            if (ini >= n) continue;
            if (++hijas > tope_hijas) { desborde = true; break; }
            const auto& h = *gen->respuesta;
            for (std::size_t t2 = ini; t2 < n; ++t2) {
                const std::size_t dt = t2 - ini;
                if (dt >= h.size()) break;
                double acc = 0.0;
                for (std::size_t k = 0; k < e_h.size() && k <= dt; ++k) acc += e_h[k] * h[dt - k];
                x_est[t2] += acc;
            }
        }
    }
    return !desborde;
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
        // Code (2026-11-11): en ARFIMA `simular` necesita 2 (burn + n) de trabajo; con menos cae a su pila de 16 384 y deja ceros
        // despues (la via (a) con una calibracion larga daba un nulo truncado). Control D-ARFIMA de la bateria L+H.
        std::vector<double> buf(2 * (n + m_lin.burn_in) + m_lin.p + 16, 0.0);
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
    bool aplicar_guarda = true,
    const GeneradorHijasNulo* hijas = nullptr) {
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
    std::vector<bool> desbordado(tiempos_iny.size(), false);
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
    if (hijas == nullptr || sin_memoria) {
        llevar_a_crudo<B_PURGA>(integracion, rama_purgada, ctx, x_est, x0_nivel, inicio, 0, x_sim, fases, dias, sin_factor);
        res.sin_factor_purga = sin_factor;
    }

    // 4. Copia limpia del operador delta y de ObservadorLilyth sobre la realizacion; se registran todas las ranuras.
    Obs d_local = delta_original;
    d_local.reiniciar_flujo();
    Lil lil_local = lil_plantilla;
    lil_local.reiniciar_flujo_conservando_purga();
    const bool hijas_en_linea = hijas != nullptr && !sin_memoria;

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
    auto registrar = [&](std::size_t, double, const auto& o, const auto& lil_v) {
        if (o.abrio_ranura) aperturas_todas.push_back(o.t);
        for (std::size_t c = 0; c < o.ranuras_cerradas; ++c) {
            const auto& rc = lil_v.cerrada(c);
            Cierre ci{};
            ci.t_apertura = rc.t_apertura;
            ci.t_cierre = o.t;
            ci.signo = rc.signo;
            const auto& rama = rama_purgada ? rc.purgada : rc.cruda;
            ci.desenlace = static_cast<std::uint8_t>(rama.desenlace);
            ci.extraible = extraer_evento_nulo<T_MAX>(rc, rama_purgada, ci.ev);
            cierres.push_back(ci);
        }
    };
    if (hijas_en_linea) {
        // E1.24 (opcion (b) del Comandante): las hijas las disparan las alarmas de ESTA realizacion, en linea.
        std::size_t sf2 = 0, nh = 0;
        if (!correr_en_linea<B_PURGA>(d_local, lil_local, integracion, rama_purgada, ctx, x0_nivel, inicio, 0, x_est, hijas, rng,
                                      REPLICAS_HIJAS_MAX, nh, sf2, registrar)) {
            std::fill(desbordado.begin(), desbordado.end(), true);
            ++res.cumulos_desbordados;
        }
        res.hijas += nh;
        res.sin_factor_purga = sf2;
    } else {
        for (std::size_t t = 0; t < n_sim; ++t) {
            const auto s = d_local.paso(x_sim[t], NaN, true, fases[t], dias[t]);
            const auto o = lil_local.paso(s, x_sim[t], true, NaN, fases[t], dias[t]);
            registrar(t, x_sim[t], o, lil_local);
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
        if (desbordado[i]) { ++res.no_medidas; continue; }   // E1.24: un cumulo truncado no es el nulo
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
        } else if (!exige_disparo_detector && !x_sim.empty()) {
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
    std::size_t& cont_solapes_global,
    const GeneradorHijasNulo* hijas = nullptr) {
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
            perfil_est, e_iny, signo_real, w_ancla, guarda, t_max_ef, t_arr, n_sim, rng, true, hijas);

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


// E1.24 (opcion (b) del Comandante, 2026-11-11): INFERENCIA INDIRECTA del proceso de replicas.
// El Hawkes de las ALARMAS (el auxiliar) mezcla las replicas con las realarmas del propio detector sobre un mismo choque (sonda: 1.76
// alarmas por choque aislado, el 89 % en 10 pasos). El nulo ya corre el detector, asi que inyectar ese Hawkes contaria las realarmas
// dos veces. El modelo ESTRUCTURAL es el nulo de Onsager (calma = residuos remuestreados y el modelo lineal) mas hijas que disparan las
// alarmas de la propia simulacion, en linea; se busca su nucleo a[m][s][j] tal que el Hawkes auxiliar de las alarmas simuladas
// (calibracion simulada, mismo detector, mismas fases) iguale al de las alarmas reales: punto fijo a <- a + (a_real - media a_sim),
// con R = m_min simulaciones por iteracion y numeros aleatorios comunes. Con a = 0 el modelo es Onsager (control R8).
// Se igualan las CUATRO MASAS N (lo que decide cuantas replicas hay); la FORMA del nucleo de cada par (madre, hija) es la del auxiliar
// real, normalizada. Igualar componente a componente con el recorte en 0 no converge: donde la simulacion pone la masa en otra tau que
// la real, los componentes negativos se recortan y la masa total SUBE hasta explotar (sonda_ii.log, 2026-11-11).
//   H.4: el nulo del contraste es el mismo simulador con a = 0 (Onsager con su detector): sus realarmas entran en el nulo.
//   Cota de rho: bootstrap parametrico bajo el estructural ajustado; el intervalo basico del auxiliar se traslada al estructural.
inline constexpr std::size_t REPLICAS_ITER_INDIRECTA_MAX = 64;   // tope de iteraciones del punto fijo (contado)

template <std::size_t T_MAX, std::size_t B_PURGA, class Obs, class Lil>
inline FichaReplicas ajustar_replicas_indirecto(
    const Obs& delta, const Lil& lil_plantilla, ark::metrologia::Integracion integracion, bool rama_purgada, const ModeloLinealNulo& m_lin,
    const std::array<std::vector<std::vector<double>>, 2>& biblioteca, const std::vector<double>& h, const ContextoPurgaNulo<B_PURGA>& ctx,
    double x0, std::size_t n_cal, std::span<const AlarmaReplica> alarmas_reales, const std::vector<std::uint32_t>& bin_de_paso,
    const ConfiguracionReplicas& cfg) {
    namespace dr = detalle_replicas;
    FichaReplicas f{};
    f.n_pasos = n_cal;
    f.m_min = ark::metrologia::familia::m_min(cfg.alfa);
    const bool base_por_fase = true;               // mutante M2: base constante
    f.bins = base_por_fase ? std::max<std::size_t>(1, cfg.bins) : 1;
    f.n_alarmas = alarmas_reales.size();
    for (const auto& a : alarmas_reales) ++f.n_tipo[a.tipo & 1u];
    if (!m_lin.valido || n_cal == 0 || bin_de_paso.size() != n_cal) {
        f.estado = EstadoReplicas::NO_DECIDIBLE; f.motivo = "sin modelo lineal o sin historia"; return f;
    }
    if (f.n_alarmas < f.m_min) {
        f.estado = EstadoReplicas::NO_DECIDIBLE; f.motivo = "alarmas insuficientes en la calibracion: el ajuste no puede superar su nulo"; return f;
    }
    if (biblioteca[0].size() < f.m_min || biblioteca[1].size() < f.m_min) {
        f.estado = EstadoReplicas::NO_DECIDIBLE; f.motivo = "biblioteca de marcas por debajo de m_min en algun signo"; return f;
    }
    std::vector<std::uint32_t> bins_paso(n_cal);
    for (std::size_t i = 0; i < n_cal; ++i) {
        if (bin_de_paso[i] >= std::max<std::size_t>(1, cfg.bins)) { f.estado = EstadoReplicas::NO_DECIDIBLE; f.motivo = "fase fuera de rango"; return f; }
        bins_paso[i] = base_por_fase ? bin_de_paso[i] : 0u;
    }
    std::array<double, REPLICAS_J_MAX> tau{};
    const std::size_t J = dr::escalera_tau(cfg.t_max, tau);
    std::vector<double> pasos_por_bin(f.bins, 0.0);
    for (std::size_t i = 0; i < n_cal; ++i) pasos_por_bin[bins_paso[i]] += 1.0;
    const bool tipos_separados = true;              // mutante M3: un solo tipo que funde los dos signos

    struct Aux { AjusteHawkesReplicas hk{}, po{}; double lr = NaN; bool ok = false; };
    auto auxiliar = [&](std::vector<AlarmaReplica>& al) {
        Aux a{};
        for (auto& x : al) x.bin = bins_paso[static_cast<std::size_t>(x.t)];
        const auto d = dr::preparar_datos(al, n_cal, pasos_por_bin, J, tau, !tipos_separados);
        a.hk = dr::ajustar_hawkes(d, tau);
        a.po = dr::ajustar_poisson(d, tau);
        a.lr = 2.0 * std::max(0.0, a.hk.ll - a.po.ll);
        a.ok = a.hk.valido && a.po.valido;
        return a;
    };
    std::vector<AlarmaReplica> reales(alarmas_reales.begin(), alarmas_reales.end());
    const Aux aux_real = auxiliar(reales);
    f.auxiliar = aux_real.hk;
    f.poisson = aux_real.po;
    f.lr = aux_real.lr;
    if (!aux_real.ok) { f.estado = EstadoReplicas::NO_DECIDIBLE; f.motivo = "ajuste auxiliar no finito"; return f; }

    const std::size_t tope_hijas = 64 * (f.n_alarmas + 1);
    const std::size_t burn = m_lin.burn_in;
    // Una calibracion simulada bajo el estructural `est` (nullptr = Onsager): devuelve las alarmas; false si las hijas desbordan.
    auto simular_cal = [&](const AjusteHawkesReplicas* est, std::uint64_t semilla, std::vector<AlarmaReplica>& al) {
        al.clear();
        std::mt19937_64 rng(semilla);
        std::vector<double> x_est(n_cal, 0.0);
        std::vector<double> buf(2 * (burn + n_cal) + m_lin.p + 16, 0.0);
        simular(m_lin, rng, nullptr, 0, 0, x_est.data(), n_cal, buf.data(), buf.size());
        auto dl = std::make_unique<Obs>(delta);
        dl->reiniciar_flujo();
        auto ll = std::make_unique<Lil>(lil_plantilla);
        ll->reiniciar_flujo_conservando_purga();
        const GeneradorHijasNulo gen{est, &biblioteca, &h};
        std::size_t nh = 0, sf = 0;
        auto visitar = [&](std::size_t t, double, const auto& o, const auto&) {
            if (!o.alarma) return;
            AlarmaReplica a{};
            a.t = static_cast<std::int64_t>(t);
            a.tipo = o.signo_alarma > 0 ? 0 : 1;
            a.escala = static_cast<std::uint8_t>(std::max(0, o.escala_alarma));
            al.push_back(a);
        };
        return correr_en_linea<B_PURGA>(*dl, *ll, integracion, rama_purgada, ctx, x0, 0, 0, x_est, est != nullptr ? &gen : nullptr, rng,
                                        tope_hijas, nh, sf, visitar);
    };

    // H.4: el contraste contra Onsager con su detector (a = 0).
    f.b_boot = static_cast<std::size_t>(std::ceil(static_cast<double>(f.m_min) / cfg.alfa)) - 1;
    const std::size_t B = f.b_boot;
    std::vector<double> lr_b(B, 0.0);
    ark::metrologia::familia::en_paralelo(B, cfg.hilos, [&](std::size_t b) {
        std::vector<AlarmaReplica> al;
        (void)simular_cal(nullptr, dr::semilla(cfg.semilla, 2 * b + 1), al);
        lr_b[b] = auxiliar(al).lr;
    });
    std::size_t ge = 0;
    for (double v : lr_b) ge += (v >= f.lr);
    f.p_lr = static_cast<double>(1 + ge) / static_cast<double>(B + 1);
    const bool exigir_su_nulo = true;               // mutante M9: L+H usado sin pasar su nulo
    if (exigir_su_nulo && !(f.p_lr <= cfg.alfa)) {
        f.estado = EstadoReplicas::NO_DECIDIBLE;
        f.motivo = "las alarmas no se agrupan mas que en el nulo de Onsager con su detector";
        return f;
    }

    // Error de muestreo del AUXILIAR REAL: bootstrap parametrico del propio Hawkes auxiliar (B historias simuladas con sus
    // parametros, sin detector), sd de cada masa N. Es la ponderacion clasica de la inferencia indirecta (la varianza del auxiliar
    // estimada sobre los datos). Hace falta porque en el borde 0 la sd bajo el estructural es casi 0: una cruzada real de 0.008-0.013
    // (3 de 10 semillas sin replicas cruzadas, sonda_ii_semillas.log) quedaba fuera de tolerancia y el punto fijo la perseguia sin
    // fin (0.009 -> 0.136). Un contraste por masa contra Onsager tampoco sirve: Onsager da casi siempre 0 exacto y cualquier
    // cruzada positiva salia p = 1/(B+1).
    std::array<std::array<double, 2>, 2> sd_real{};
    {
        std::vector<std::vector<std::int64_t>> pasos_de_bin(f.bins);
        for (std::size_t i = 0; i < n_cal; ++i) pasos_de_bin[bins_paso[i]].push_back(static_cast<std::int64_t>(i));
        std::vector<std::array<std::array<double, 2>, 2>> nb(B);
        std::vector<std::uint8_t> ok(B, 0);
        const bool sd_del_auxiliar_real = true;     // mutante M11: sin el error del auxiliar real (vuelve la deriva de la cruzada)
        if (sd_del_auxiliar_real && aux_real.hk.rho < 1.0) {
            ark::metrologia::familia::en_paralelo(B, cfg.hilos, [&](std::size_t b) {
                std::mt19937_64 rng(dr::semilla(cfg.semilla, 5000011 + b));
                std::vector<AlarmaReplica> sim;
                if (!dr::simular_historia(aux_real.hk, true, n_cal, pasos_de_bin, tope_hijas, rng, sim, bins_paso)) return;
                const Aux a = auxiliar(sim);
                if (!a.ok) return;
                for (std::size_t m = 0; m < 2; ++m) for (std::size_t s = 0; s < 2; ++s) nb[b][m][s] = a.hk.N[m][s];
                ok[b] = 1;
            });
        }
        for (std::size_t m = 0; m < 2; ++m)
            for (std::size_t s = 0; s < 2; ++s) {
                double sm = 0.0, s2 = 0.0, k = 0.0;
                for (std::size_t b = 0; b < B; ++b) if (ok[b]) { sm += nb[b][m][s]; s2 += nb[b][m][s] * nb[b][m][s]; k += 1.0; }
                sd_real[m][s] = k > 1.0 ? std::sqrt(std::max(0.0, s2 / k - (sm / k) * (sm / k)) * k / (k - 1.0)) : 0.0;
                f.sd_masa_real[m][s] = sd_real[m][s];
            }
    }

    // Punto fijo de la inferencia indirecta sobre las masas (numeros aleatorios comunes: las mismas R semillas en cada iteracion).
    const std::size_t R = f.m_min;
    AjusteHawkesReplicas est = aux_real.hk;
    // forma del nucleo de cada par: w[m][s][j] = (a_aux / (1 - q_j)) / N_aux, la contribucion de cada tau a la masa
    std::array<std::array<std::array<double, REPLICAS_J_MAX>, 2>, 2> forma{};
    std::array<std::array<double, 2>, 2> masa{};
    for (std::size_t m = 0; m < 2; ++m)
        for (std::size_t s = 0; s < 2; ++s) {
            masa[m][s] = aux_real.hk.N[m][s];
            for (std::size_t j = 0; j < J; ++j) {
                const double c = aux_real.hk.a[m][s][j] / (1.0 - std::exp(-1.0 / tau[j]));
                forma[m][s][j] = aux_real.hk.N[m][s] > 0.0 ? c / aux_real.hk.N[m][s] : (j == 0 ? 1.0 : 0.0);
            }
        }
    auto fijar_desde_masa = [&]() {
        for (std::size_t m = 0; m < 2; ++m)
            for (std::size_t s = 0; s < 2; ++s)
                for (std::size_t j = 0; j < J; ++j) est.a[m][s][j] = masa[m][s] * forma[m][s][j] * (1.0 - std::exp(-1.0 / tau[j]));
        dr::calcular_masas(est);
    };
    fijar_desde_masa();
    std::vector<Aux> sims(R);
    std::vector<std::uint8_t> expl(R, 0);
    // Ligadura b(theta) para rho, la mayor masa de mismo signo y la mayor de signo contrario: valor del estructural y media del
    // auxiliar simulado, en la primera iteracion no explosiva y en la ultima (numeros aleatorios comunes). Da la pendiente secante.
    struct PuntoLigadura { double est[3] = {NaN, NaN, NaN}; double sim[3] = {NaN, NaN, NaN}; std::array<std::vector<double>, 3> v{}; bool ok = false; };
    PuntoLigadura lig_1{}, lig_f{};
    auto punto_ligadura = [&]() {
        PuntoLigadura q{};
        q.est[0] = est.rho;
        q.est[1] = std::max(est.N[0][0], est.N[1][1]);
        q.est[2] = std::max(est.N[0][1], est.N[1][0]);
        for (std::size_t g = 0; g < 3; ++g) { q.sim[g] = 0.0; q.v[g].reserve(R); }
        for (const auto& x : sims) {
            const double val[3] = {x.hk.rho, std::max(x.hk.N[0][0], x.hk.N[1][1]), std::max(x.hk.N[0][1], x.hk.N[1][0])};
            for (std::size_t g = 0; g < 3; ++g) { q.sim[g] += val[g] / double(R); q.v[g].push_back(val[g]); }
        }
        q.ok = true;
        return q;
    };
    for (std::size_t it = 1; it <= REPLICAS_ITER_INDIRECTA_MAX; ++it) {
        dr::calcular_masas(est);
        ark::metrologia::familia::en_paralelo(R, cfg.hilos, [&](std::size_t r) {
            std::vector<AlarmaReplica> al;
            expl[r] = simular_cal(&est, dr::semilla(cfg.semilla, 1000003 + r), al) ? 0 : 1;
            sims[r] = auxiliar(al);
        });
        f.iter_indirecta = it;
        std::size_t ne = 0;
        for (auto e : expl) ne += e;
        f.sims_explosivas_ii += ne;
        if (ne > 0) {   // el estructural es demasiado fuerte para esta historia: se amortigua y se sigue
            for (auto& x : masa) for (auto& y : x) y *= 0.5;
            fijar_desde_masa();
            continue;
        }
        // convergencia sobre las masas N del auxiliar: |N_real - media N_sim| <= 2 sd de UNA realizacion, con la sd tomada como la
        // mayor de sus dos estimaciones: bajo el estructural (la de las R simulaciones) y la del auxiliar real (bootstrap de arriba).
        // La calibracion real es una sola historia: basta con que sea plausible. Con el ee de la media (sd / sqrt(R)) una masa pegada
        // al borde 0, donde el auxiliar casi no responde al estructural, seguia sumando sin fin (sonda_ii2.log).
        bool conv = true;
        double se_max = 0.0;
        for (std::size_t m = 0; m < 2; ++m)
            for (std::size_t s = 0; s < 2; ++s) {
                double sm = 0.0, s2 = 0.0;
                for (const auto& x : sims) { sm += x.hk.N[m][s]; s2 += x.hk.N[m][s] * x.hk.N[m][s]; }
                const double med = sm / double(R);
                const double sd = std::sqrt(std::max(0.0, s2 / double(R) - med * med) * double(R) / double(R - 1));
                const double tol = 2.0 * std::max(sd, sd_real[m][s]);
                se_max = std::max(se_max, std::max(sd, sd_real[m][s]));
                if (std::abs(aux_real.hk.N[m][s] - med) > tol) conv = false;
            }
        f.se_max_ii = se_max;
#ifdef ARK_LILYTH_DEPURAR_II
        {
            double nm[2][2] = {{0, 0}, {0, 0}};
            for (const auto& x : sims) for (std::size_t m = 0; m < 2; ++m) for (std::size_t s = 0; s < 2; ++s) nm[m][s] += x.hk.N[m][s] / double(R);
            std::printf("      [ii %2zu] real N %.3f %.3f %.3f %.3f | sim N %.3f %.3f %.3f %.3f | est N %.3f %.3f %.3f %.3f | se %.3f\n", it,
                        aux_real.hk.N[0][0], aux_real.hk.N[0][1], aux_real.hk.N[1][0], aux_real.hk.N[1][1], nm[0][0], nm[0][1], nm[1][0], nm[1][1],
                        est.N[0][0], est.N[0][1], est.N[1][0], est.N[1][1], se_max);
        }
#endif
        if (!lig_1.ok) lig_1 = punto_ligadura();
        lig_f = punto_ligadura();
        if (conv) { f.convergio = true; break; }
        // Solo se mueve la masa que esta fuera de su tolerancia. Moverlas todas mientras otra no convergia arrastraba la cruzada, que
        // ya era plausible, de 0.009 a 0.136 (sonda_ii3.log): el auxiliar casi no responde a ella y cada paso le sumaba ruido.
        for (std::size_t m = 0; m < 2; ++m)
            for (std::size_t s = 0; s < 2; ++s) {
                double sm = 0.0, s2 = 0.0;
                for (const auto& x : sims) { sm += x.hk.N[m][s]; s2 += x.hk.N[m][s] * x.hk.N[m][s]; }
                const double med = sm / double(R);
                const double sd = std::sqrt(std::max(0.0, s2 / double(R) - med * med) * double(R) / double(R - 1));
                if (std::abs(aux_real.hk.N[m][s] - med) <= 2.0 * std::max(sd, sd_real[m][s])) continue;
                masa[m][s] = std::max(0.0, masa[m][s] + (aux_real.hk.N[m][s] - med));
            }
        fijar_desde_masa();
    }
    dr::calcular_masas(est);
    f.hawkes = est;
    f.rho = est.rho;
    f.n_mismo_signo = std::max(est.N[0][0], est.N[1][1]);
    f.n_signo_contrario = std::max(est.N[0][1], est.N[1][0]);
    if (!f.convergio) {
        f.estado = EstadoReplicas::NO_DECIDIBLE; f.motivo = "la inferencia indirecta no converge"; return f;
    }

    // Pendiente de la ligadura: secante entre el estructural final y el mismo con las masas a la mitad (la misma amortiguacion que tras
    // un desborde), con los mismos numeros aleatorios. Siempre el mismo segundo punto, para que las tres pendientes salgan del mismo
    // par; la primera iteracion del punto fijo queda de respaldo si la mitad desborda.
    if (lig_f.ok) {
        const auto masa_guardada = masa;
        for (auto& x : masa) for (auto& y : x) y *= 0.5;
        fijar_desde_masa();
        bool desbordo = false;
        ark::metrologia::familia::en_paralelo(R, cfg.hilos, [&](std::size_t r) {
            std::vector<AlarmaReplica> al;
            expl[r] = simular_cal(&est, dr::semilla(cfg.semilla, 1000003 + r), al) ? 0 : 1;
            sims[r] = auxiliar(al);
        });
        for (auto e : expl) desbordo = desbordo || e != 0;
        if (!desbordo) lig_1 = punto_ligadura();
        masa = masa_guardada;
        fijar_desde_masa();
    }
    // pendiente secante de cada estadistico (rho, mayor N de mismo signo, mayor N contrario): d media_sim / d estructural
    double pend[3] = {NaN, NaN, NaN};
    const bool invertir_ligadura = true;            // mutante M12: pendiente 1 (traslado del intervalo del auxiliar tal cual)
    // La pendiente solo cuenta si la respuesta del auxiliar esta RESUELTA: la diferencia pareada (mismos numeros aleatorios) entre los
    // dos puntos supera 2 errores estandar. Una pendiente positiva pero no resuelta dividia y daba masas de 1e7 y 1e29
    // (etapa_lh_rapido2.log, R2 n = 0.6).
    for (std::size_t g = 0; g < 3; ++g) {
        const double de = lig_f.est[g] - lig_1.est[g];
        double dm = 0.0, d2 = 0.0;
        const std::size_t k = std::min(lig_f.v[g].size(), lig_1.v[g].size());
        for (std::size_t r = 0; r < k; ++r) { const double dd = lig_f.v[g][r] - lig_1.v[g][r]; dm += dd; d2 += dd * dd; }
        bool resuelta = false;
        if (k > 1) {
            dm /= double(k);
            const double sd = std::sqrt(std::max(0.0, d2 / double(k) - dm * dm) * double(k) / double(k - 1));
            resuelta = de > 0.0 && dm > 2.0 * sd / std::sqrt(double(k));
        }
        pend[g] = invertir_ligadura ? (resuelta ? (lig_f.sim[g] - lig_1.sim[g]) / de : NaN) : 1.0;
    }
    f.pendiente_ligadura = pend[0];

    // Cotas: bootstrap bajo el estructural; el intervalo basico del auxiliar se traslada al estructural invirtiendo la ligadura
    // lineal: theta_sup = theta + (b_real - q_alfa(b*)) / pendiente. El punto se corrige con el mismo paso de Newton hacia la
    // mediana: theta_c = theta + (b_real - mediana(b*)) / pendiente. Sin pendiente (el auxiliar no responde, o es <= 0) rho no es
    // decidible; una masa sin pendiente queda en [0, inf): no esta identificada por este detector.
    std::vector<double> rho_b(B, NaN), nd_b(B, NaN), nc_b(B, NaN);
    std::vector<std::uint8_t> ex_b(B, 0);
    if (est.rho < 1.0) {
        ark::metrologia::familia::en_paralelo(B, cfg.hilos, [&](std::size_t b) {
            std::vector<AlarmaReplica> al;
            if (!simular_cal(&est, dr::semilla(cfg.semilla, 2 * b + 2), al)) {
                ex_b[b] = 1; rho_b[b] = nd_b[b] = nc_b[b] = std::numeric_limits<double>::infinity(); return;
            }
            const Aux a = auxiliar(al);
            rho_b[b] = a.hk.rho;
            nd_b[b] = std::max(a.hk.N[0][0], a.hk.N[1][1]);
            nc_b[b] = std::max(a.hk.N[0][1], a.hk.N[1][0]);
        });
        for (auto e : ex_b) f.boot_explosivos += e;
        auto cuantil = [](std::vector<double> v, double q) {
            std::sort(v.begin(), v.end());
            const double pos = q * static_cast<double>(v.size() - 1);
            const std::size_t i = static_cast<std::size_t>(std::floor(pos));
            if (i + 1 >= v.size()) return v.back();
            if (!std::isfinite(v[i + 1])) return v[i + 1];
            return v[i] + (pos - static_cast<double>(i)) * (v[i + 1] - v[i]);
        };
        const double real[3] = {aux_real.hk.rho, std::max(aux_real.hk.N[0][0], aux_real.hk.N[1][1]),
                                std::max(aux_real.hk.N[0][1], aux_real.hk.N[1][0])};
        const std::vector<double>* bs[3] = {&rho_b, &nd_b, &nc_b};
        double punto[3] = {f.rho, f.n_mismo_signo, f.n_signo_contrario};
        double inf_[3], sup_[3];
        const double a_lo[3] = {cfg.alfa, cfg.alfa / 2.0, cfg.alfa / 2.0};
        for (std::size_t g = 0; g < 3; ++g) {
            const double s = pend[g];
            if (!(std::isfinite(s) && s > 0.0)) {
                inf_[g] = 0.0; sup_[g] = std::numeric_limits<double>::infinity();
                continue;
            }
            const double c = punto[g];
            punto[g] = c + (real[g] - cuantil(*bs[g], 0.5)) / s;
            sup_[g] = c + (real[g] - cuantil(*bs[g], a_lo[g])) / s;
            inf_[g] = c + (real[g] - cuantil(*bs[g], 1.0 - a_lo[g])) / s;
        }
        if (!(std::isfinite(pend[0]) && pend[0] > 0.0)) {
            f.estado = EstadoReplicas::NO_DECIDIBLE;
            f.motivo = "la ligadura de rho no responde: el auxiliar no identifica la ramificacion";
            return f;
        }
        f.rho = punto[0]; f.rho_sup = sup_[0]; f.rho_inf = inf_[0];
        f.n_mismo_signo = punto[1]; f.n_mismo_signo_sup = sup_[1]; f.n_mismo_signo_inf = inf_[1];
        f.n_signo_contrario = punto[2]; f.n_contrario_sup = sup_[2]; f.n_contrario_inf = inf_[2];
    } else {
        f.rho_sup = std::numeric_limits<double>::infinity();
    }
    const bool rechazar_ramificacion_critica = true; // mutante M7: rho >= 1 aceptado
    if (rechazar_ramificacion_critica && (!(f.rho < 1.0) || !(f.rho_sup < 1.0))) {
        f.estado = EstadoReplicas::NO_IDENTIFICADO;
        f.motivo = "ramificacion critica: rho(N) o su cota superior >= 1";
    } else {
        f.estado = EstadoReplicas::DECIDIBLE;
        f.motivo = "replicas por encima de Onsager con su detector; estructural por inferencia indirecta";
    }
    return f;
}

}  // namespace detalle_nulo

}  // namespace ark::lilyth