// contratos_ventanas.hpp -- CONTRATOS Y ESTRUCTURAS ZERO-HEAP DE VENTANAS CAUSALES (ARK-WIN-01-CXX)
// Reactor ARK / Canción del Micelio -- ARQUITECTO DE DATOS & COMANDANTE
// Corregido por el Jefe de Arquitectos de Código (2026-11-01), acta AUDITORIA_CAUSAL/AUDITORIA_VENTANAS_CXX.md:
//   E1  suelo de sigma dinámico: max(sigma_floor_minimo, tick/P_t), con el tick del instrumento (sin literal 0.01) y P_t el
//       precio de cada paso; el reloj de Shannon exige el tick declarado.
//   E10 la configuración se valida como el __post_init__ del Python (W >= min_muestras_requeridas, stride >= 1, ...).
//   E3  las celdas sin sesión no cuentan en el presupuesto de imputación (ni arriba ni abajo).
//   B5  la ventana lleva el sello de disponibilidad (t_actual) y, aparte, el de cierre (t_cierre_ms).
//   E7  un activo cerrado sin historia vale NaN, no 0.0.
//   E13 comentarios ajustados al código.
#pragma once

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <array>
#include <span>
#include <string_view>
#include <algorithm>
#include <cstring>
#include <limits>

namespace ark::ventanas {

// =======================================================================
// 1. Enumeraciones Ontológicas Compactas (Alineadas en uint8_t)
// =======================================================================

enum class TipoRelojVentana : uint8_t {
    PARED = 0,             // Reloj de calendario (1.0 por barra con sesión)
    NU_INFO_SHANNON = 1,   // Sorpresa gaussiana en bits: d_nu = z^2 / (2 ln 2), z = (r - mu) / sigma con la historia real
    NU_VOLUMEN = 2,        // Reloj acumulado por volumen de intercambio
    NU_TRADES = 3,         // Reloj acumulado por número de transacciones
    NU_NOCIONAL = 4        // Reloj acumulado por volumen nocional (USD)
};

enum class MetodoImputacion : uint8_t {
    FORWARD_FILL = 0,      // Propagación causal estricta del último valor conocido
    CONGELAR_RELOJ = 1,    // Forward-fill congelando rígidamente el reloj intrínseco (d_nu = 0)
    PROHIBIR_INTERPOLACION = 2 // Rechazo estricto: aborta ante cualquier NaN, hueco o imputación
};

enum class PoliticaSincronizacion : uint8_t {
    INTERSECCION_ESTRICTA = 0, // Solo filas en las que los M activos tienen sesión (las demás se descartan sin tocar el estado)
    UNION_CON_MASCARA = 1      // Rejilla de unión causal con arrastre y market_open_mask
};

enum class ClasePar : uint8_t {
    EJECUTABLE = 0,        // La información de u está completa (+ latencia) antes de que empiece el retorno de v
    REZAGO = 1,            // La información de u llega después de que empiece el retorno de v: usarla es mirar el futuro
    NO_VERIFICADO = 2,     // Sin horas de cierre declaradas, o intervalos que se solapan
    SIN_SESION = 3         // El objetivo v no tiene sesión en la fila: no hay retorno que predecir
};

// =======================================================================
// 2. Funciones Utilitarias constexpr de Nombres
// =======================================================================

[[nodiscard]] constexpr std::string_view nombre_reloj(TipoRelojVentana r) noexcept {
    switch (r) {
        case TipoRelojVentana::PARED: return "PARED";
        case TipoRelojVentana::NU_INFO_SHANNON: return "NU_INFO_SHANNON";
        case TipoRelojVentana::NU_VOLUMEN: return "NU_VOLUMEN";
        case TipoRelojVentana::NU_TRADES: return "NU_TRADES";
        case TipoRelojVentana::NU_NOCIONAL: return "NU_NOCIONAL";
    }
    return "DESCONOCIDO";
}

[[nodiscard]] constexpr std::string_view nombre_clase_par(ClasePar c) noexcept {
    switch (c) {
        case ClasePar::EJECUTABLE: return "EJECUTABLE";
        case ClasePar::REZAGO: return "REZAGO";
        case ClasePar::NO_VERIFICADO: return "NO_VERIFICADO";
        case ClasePar::SIN_SESION: return "SIN_SESION";
    }
    return "DESCONOCIDO";
}

// Motivo por el que una configuración no es admisible (E10). NINGUNO = válida.
enum class ErrorConfiguracion : uint8_t {
    NINGUNO = 0,
    W_MENOR_QUE_MINIMO,          // W < min_muestras_requeridas (Python: lookback_W < min_muestras_requeridas)
    STRIDE_NULO,                 // stride_S < 1
    UMBRAL_FUERA_DE_RANGO,       // umbral_imputacion_maximo fuera de [0, 1]
    S_MIN_INVALIDO,              // s_min < 2
    RATIO_ESCALAS_INVALIDO,      // max_ratio_escalas <= 1
    TICK_INVALIDO,               // tick declarado <= 0
    TICK_NO_DECLARADO,           // reloj de Shannon sin tick: el suelo de sigma no puede ser dinámico (E1)
    SUELO_SIGMA_INVALIDO,        // sigma_floor_minimo <= 0
    VENTANA_INNOVACION_INVALIDA, // ventana_innovacion_shannon fuera de {0} U [2, W]
    DELTA_T_INVALIDO,            // delta_t_barra_segundos < 0
    JITTER_INVALIDO,             // tolerancia < 0, o >= delta_t nominal
    LATENCIA_INVALIDA            // latencia_guarda_s < 0 o no finita
};

[[nodiscard]] constexpr std::string_view nombre_error_configuracion(ErrorConfiguracion e) noexcept {
    switch (e) {
        case ErrorConfiguracion::NINGUNO: return "NINGUNO";
        case ErrorConfiguracion::W_MENOR_QUE_MINIMO: return "W_MENOR_QUE_MINIMO";
        case ErrorConfiguracion::STRIDE_NULO: return "STRIDE_NULO";
        case ErrorConfiguracion::UMBRAL_FUERA_DE_RANGO: return "UMBRAL_FUERA_DE_RANGO";
        case ErrorConfiguracion::S_MIN_INVALIDO: return "S_MIN_INVALIDO";
        case ErrorConfiguracion::RATIO_ESCALAS_INVALIDO: return "RATIO_ESCALAS_INVALIDO";
        case ErrorConfiguracion::TICK_INVALIDO: return "TICK_INVALIDO";
        case ErrorConfiguracion::TICK_NO_DECLARADO: return "TICK_NO_DECLARADO";
        case ErrorConfiguracion::SUELO_SIGMA_INVALIDO: return "SUELO_SIGMA_INVALIDO";
        case ErrorConfiguracion::VENTANA_INNOVACION_INVALIDA: return "VENTANA_INNOVACION_INVALIDA";
        case ErrorConfiguracion::DELTA_T_INVALIDO: return "DELTA_T_INVALIDO";
        case ErrorConfiguracion::JITTER_INVALIDO: return "JITTER_INVALIDO";
        case ErrorConfiguracion::LATENCIA_INVALIDA: return "LATENCIA_INVALIDA";
    }
    return "DESCONOCIDO";
}

// Suelo dinámico de sigma (Axioma 4, E1): sigma_floor(t) = max(sigma_floor_minimo, tick / P_t). Es el del Python canónico
// (ConfiguracionVentanaCausal.obtener_sigma_floor) y el de la bitácora del Arquitecto: un retorno menor que un tick
// relativo no es información, y el tick relativo cambia con el precio de cada paso. Sin tick declarado (NaN), o con un precio
// no positivo, queda el suelo mínimo.
[[nodiscard]] inline double sigma_floor_dinamico(double sigma_floor_minimo, double tick, double precio_actual) noexcept {
    if (precio_actual > 0.0 && tick > 0.0 && std::isfinite(tick)) {
        return std::max(sigma_floor_minimo, tick / precio_actual);
    }
    return sigma_floor_minimo;
}

// =======================================================================
// 3. Configuración Causal de Ventana (Plantilla de Tamaño Fijo W)
// =======================================================================

template <size_t W>
struct ConfiguracionVentanaCausalCxx {
    static_assert(W >= 4, "lookback_W debe ser al menos 4 barras (y al menos min_muestras_requeridas en tiempo de ejecución).");

    size_t min_muestras_requeridas{32};         // Python: W_min = 2^(J+2) = 32; el motor rechaza W menor (E10)
    size_t stride_S{1};
    size_t s_min{4};
    double max_ratio_escalas{4.0};
    bool forzar_divisores_exactos{true};
    bool bidireccional_kantelhardt{true};

    TipoRelojVentana reloj{TipoRelojVentana::PARED};
    MetodoImputacion metodo_imputacion{MetodoImputacion::FORWARD_FILL};

    double umbral_imputacion_maximo{0.20};      // Máximo 20 % de imputación sobre las celdas CON SESIÓN (Decisión #2, E3)
    bool abortar_si_excede_imputacion{false};
    bool congelar_en_huecos{true};

    bool estampar_al_cierre{false};            // Ingesta R1 es la autoridad; False por defecto
    int64_t delta_t_barra_segundos{0};         // 0 = sin cadencia nominal (los huecos solo llegan declarados por el calendario)
    double tolerancia_jitter_segundos{0.0};    // Python: 0.0 (E10); debe ser < delta_t nominal

    // Reloj de Shannon y suelo dinámico de sigma (E1)
    double sigma_floor_minimo{1e-4};
    double tick_size{std::numeric_limits<double>::quiet_NaN()};  // tick del instrumento; NaN = no declarado
    size_t ventana_innovacion_shannon{0};      // 0 = adaptativa por feedback tau_1/2
    bool lazo_feedback_shannon{true};

    // Guarda de pares: latencia (segundos) que se exige entre el final de la información de u y el inicio del retorno de v.
    // 0.0 = diaria (ventanas_informacion.LATENCIA_DIARIA_S); en barras intradía, la medida (0.208 s).
    double latencia_guarda_s{0.0};

    PoliticaSincronizacion politica_sincronizacion{PoliticaSincronizacion::INTERSECCION_ESTRICTA};

    // Suelo del paso (mono-activo: el tick de la configuración).
    [[nodiscard]] double obtener_sigma_floor(double precio_actual) const noexcept {
        return sigma_floor_dinamico(sigma_floor_minimo, tick_size, precio_actual);
    }

    // Retornos reales mínimos en la historia antes de emitir d_nu (E2): las dos primeras sorpresas no se miden contra el suelo.
    // Se deriva del soporte mínimo de la ventana de innovación, max(2 s_min, 8).
    [[nodiscard]] constexpr size_t k_min_shannon() const noexcept {
        return std::max<size_t>(2 * s_min, 8);
    }

    // Validación (E10). exigir_tick: el reloj de Shannon sin tick declarado no tiene suelo dinámico (E1).
    [[nodiscard]] ErrorConfiguracion validar(bool exigir_tick = true) const noexcept {
        if (W < min_muestras_requeridas) return ErrorConfiguracion::W_MENOR_QUE_MINIMO;
        if (stride_S < 1) return ErrorConfiguracion::STRIDE_NULO;
        if (!(umbral_imputacion_maximo >= 0.0 && umbral_imputacion_maximo <= 1.0)) return ErrorConfiguracion::UMBRAL_FUERA_DE_RANGO;
        if (s_min < 2) return ErrorConfiguracion::S_MIN_INVALIDO;
        if (!(max_ratio_escalas > 1.0)) return ErrorConfiguracion::RATIO_ESCALAS_INVALIDO;
        if (!std::isnan(tick_size) && !(tick_size > 0.0 && std::isfinite(tick_size))) return ErrorConfiguracion::TICK_INVALIDO;
        if (exigir_tick && reloj == TipoRelojVentana::NU_INFO_SHANNON && std::isnan(tick_size)) {
            return ErrorConfiguracion::TICK_NO_DECLARADO;
        }
        if (!(sigma_floor_minimo > 0.0)) return ErrorConfiguracion::SUELO_SIGMA_INVALIDO;
        if (ventana_innovacion_shannon != 0 && (ventana_innovacion_shannon < 2 || ventana_innovacion_shannon > W)) {
            return ErrorConfiguracion::VENTANA_INNOVACION_INVALIDA;
        }
        if (delta_t_barra_segundos < 0) return ErrorConfiguracion::DELTA_T_INVALIDO;
        if (!(tolerancia_jitter_segundos >= 0.0)) return ErrorConfiguracion::JITTER_INVALIDO;
        if (delta_t_barra_segundos > 0 && tolerancia_jitter_segundos >= static_cast<double>(delta_t_barra_segundos)) {
            return ErrorConfiguracion::JITTER_INVALIDO;
        }
        if (!(latencia_guarda_s >= 0.0 && std::isfinite(latencia_guarda_s))) return ErrorConfiguracion::LATENCIA_INVALIDA;
        return ErrorConfiguracion::NINGUNO;
    }
};

// =======================================================================
// 4. Ventana Causal Certificada Mono-Activo (Zero-Heap W)
// =======================================================================

template <size_t W>
struct VentanaCausalCertificadaCxx {
    static_assert(W >= 4, "W debe ser al menos 4 muestras.");

    // Sello de la ventana = estampa de DISPONIBILIDAD de su última fila (B5): el instante a partir del cual todos sus datos
    // existen. Si el que llama no la declara, es el sello de cierre y disponibilidad_declarada = false.
    uint64_t t_actual{0};
    uint64_t t_cierre_ms{0};                      // Sello de cierre de la última fila (el ts que entró, + offset si se estampa)
    bool disponibilidad_declarada{false};
    std::array<double, W> valores{};              // Valores observables [t - W + 1, t]
    std::array<bool, W> mascara_huecos{};         // true = imputado/rellenado, false = dato real
    std::array<bool, W> con_sesion{};             // false = fila de mercado cerrado (no cuenta en el presupuesto, E3)
    double reloj_acumulado{0.0};                  // Reloj acumulado nu(t)
    double delta_reloj_ventana{0.0};              // Suma de d_nu en las W barras de la ventana

    // Certificados Medidos Reales (No literales)
    size_t num_divisores_exactos{0};
    std::array<size_t, 32> divisores_exactos{};   // Catálogo exacto de sub-escalas s
    size_t residuo_temporal_maximo{0};            // max(W % s) medido analíticamente
    double porcentaje_descarte_causal{0.0};       // 0.0% con Kantelhardt; medido sin él
    double lookahead_score{0.0};                  // Lo registra el que audita con TestigoCausalPipelineCxx

    bool es_valida{true};
    double ratio_imputacion{0.0};                 // (imputadas con sesión + omitidas) / (celdas con sesión + omitidas); NaN si 0/0
    uint32_t barras_omitidas{0};                  // Barras perdidas dentro de la ventana (calendario o cadencia nominal)

    TipoRelojVentana tipo_reloj{TipoRelojVentana::PARED};
    uint8_t venue_id{0};                          // Identificador numérico del exchange
    char nivel[8]{"L1"};
    char simbolo[16]{"SYNTH"};

    // Métodos en línea O(1)
    [[nodiscard]] constexpr double valor_presente() const noexcept {
        return valores[W - 1];
    }

    [[nodiscard]] constexpr bool tiene_discontinuidad() const noexcept {
        return barras_omitidas > 0;
    }

    // Normalización Z-score causal estrictamente local O(W)
    void normalizar_zscore_causal(double* out_z) const noexcept {
        if (!out_z) return;
        double sum = 0.0;
        for (size_t i = 0; i < W; ++i) sum += valores[i];
        double mean = sum / static_cast<double>(W);

        double var_sum = 0.0;
        for (size_t i = 0; i < W; ++i) {
            double d = valores[i] - mean;
            var_sum += d * d;
        }
        double std_dev = (W > 1) ? std::sqrt(var_sum / static_cast<double>(W - 1)) : 0.0;

        if (std_dev < 1e-12) {
            for (size_t i = 0; i < W; ++i) out_z[i] = 0.0;
        } else {
            for (size_t i = 0; i < W; ++i) out_z[i] = (valores[i] - mean) / std_dev;
        }
    }

    // Eliminación de tendencia lineal (Proyector aniquilador M_X grado 1)
    void eliminar_tendencia_lineal(double* out_detrend) const noexcept {
        if (!out_detrend) return;
        const double w = static_cast<double>(W);
        const double sum_t = 0.5 * static_cast<double>((W - 1) * W);
        const double sum_t2 = static_cast<double>((W - 1) * W * (2 * W - 1)) / 6.0;
        double sum_y = 0.0, sum_ty = 0.0;
        for (size_t t = 0; t < W; ++t) {
            sum_y += valores[t];
            sum_ty += static_cast<double>(t) * valores[t];
        }
        const double denom = w * sum_t2 - sum_t * sum_t;
        const double b = (denom != 0.0) ? (w * sum_ty - sum_t * sum_y) / denom : 0.0;
        const double a = (sum_y - b * sum_t) / w;
        for (size_t t = 0; t < W; ++t) {
            out_detrend[t] = valores[t] - (a + b * static_cast<double>(t));
        }
    }

    // Eliminación de tendencia cuadrática (Proyector aniquilador M_X grado 2 vía polinomios ortogonales de Gram)
    void eliminar_tendencia_cuadratica(double* out_detrend) const noexcept {
        if (!out_detrend) return;
        const double w = static_cast<double>(W);
        const double t_bar = 0.5 * static_cast<double>(W - 1);
        const double s_p1_sq = w * (w * w - 1.0) / 12.0;
        const double s_p2_sq = w * (w * w - 1.0) * (w * w - 4.0) / 180.0;
        const double p2_shift = (w * w - 1.0) / 12.0;

        double sum_y = 0.0;
        double sum_y_p1 = 0.0;
        double sum_y_p2 = 0.0;

        for (size_t t = 0; t < W; ++t) {
            double dt = static_cast<double>(t) - t_bar;
            double p1 = dt;
            double p2 = dt * dt - p2_shift;
            double y = valores[t];
            sum_y += y;
            sum_y_p1 += y * p1;
            sum_y_p2 += y * p2;
        }

        double c0 = sum_y / w;
        double c1 = sum_y_p1 / s_p1_sq;
        double c2 = sum_y_p2 / s_p2_sq;

        for (size_t t = 0; t < W; ++t) {
            double dt = static_cast<double>(t) - t_bar;
            double p1 = dt;
            double p2 = dt * dt - p2_shift;
            out_detrend[t] = valores[t] - (c0 + c1 * p1 + c2 * p2);
        }
    }
};

// =======================================================================
// 5. Ventana Causal Certificada Multiactivo (Zero-Heap W x M)
// =======================================================================

template <size_t W, size_t M>
struct VentanaCausalCertificadaMultiactivoCxx {
    static_assert(W >= 4, "W debe ser al menos 4 muestras.");
    static_assert(M >= 1 && M <= 16, "M debe estar entre 1 y 16 activos.");

    uint64_t t_actual{0};                                // Disponibilidad de la ventana: la mayor de su última fila (B5)
    uint64_t t_cierre_ms{0};                             // Sello de la rejilla (cierre) de la última fila
    bool disponibilidad_declarada{false};
    std::array<std::array<double, M>, W> valores{};      // Matriz (W, M): [fila_tiempo][col_activo]
    std::array<std::array<bool, M>, W> mascara_huecos{}; // Matriz (W, M) de imputaciones
    std::array<std::array<bool, M>, W> market_open_mask{};// Matriz (W, M) de mercado abierto

    // Tiempos propios nu independientes por activo (Axioma 4 / Propuesta 5: 3b)
    std::array<double, M> delta_reloj_por_activo{};
    std::array<double, M> reloj_acum_por_activo{};
    double delta_reloj_ventana{0.0};                     // Media agregada de la ventana
    double reloj_acumulado{0.0};                         // Media acumulada

    // Certificados medidos reales
    size_t num_divisores_exactos{0};
    std::array<size_t, 32> divisores_exactos{};
    size_t residuo_temporal_maximo{0};
    double porcentaje_descarte_causal{0.0};
    double lookahead_score{0.0};

    bool es_valida{true};
    double ratio_imputacion_global{0.0};                 // sobre las celdas con sesión (E3); NaN si la ventana no tiene ninguna
    size_t celdas_con_sesion{0};
    size_t celdas_imputadas_con_sesion{0};
    uint32_t barras_omitidas{0};

    // Guarda de sincronización por pares dirigidos (u -> v) (Axioma 3 / 3a), con su margen en segundos (NaN si no se pudo)
    std::array<std::array<ClasePar, M>, M> clasificacion_pares{};
    std::array<std::array<double, M>, M> margen_pares_s{};
    char simbolos[M][16]{};
    TipoRelojVentana tipo_reloj{TipoRelojVentana::PARED};

    // Propiedades de guarda causal
    [[nodiscard]] bool tiene_rezago_causal() const noexcept {
        for (size_t u = 0; u < M; ++u) {
            for (size_t v = 0; v < M; ++v) {
                if (u != v && clasificacion_pares[u][v] == ClasePar::REZAGO) {
                    return true;
                }
            }
        }
        return false;
    }

    [[nodiscard]] size_t num_pares_rezago() const noexcept {
        size_t count = 0;
        for (size_t u = 0; u < M; ++u) {
            for (size_t v = 0; v < M; ++v) {
                if (u != v && clasificacion_pares[u][v] == ClasePar::REZAGO) {
                    ++count;
                }
            }
        }
        return count;
    }

    [[nodiscard]] constexpr double valor_presente(size_t activo_idx) const noexcept {
        return (activo_idx < M) ? valores[W - 1][activo_idx] : std::numeric_limits<double>::quiet_NaN();
    }

    // Diagnóstico de asimetría temporal bajo UNION_CON_MASCARA (Propuesta 5: 3c)
    [[nodiscard]] bool tiene_asimetria_union() const noexcept {
        for (size_t t = 0; t < W; ++t) {
            for (size_t m = 0; m < M; ++m) {
                if (!market_open_mask[t][m]) return true;
            }
        }
        return false;
    }

    [[nodiscard]] size_t barras_inactivas_activo(size_t m) const noexcept {
        if (m >= M) return 0;
        size_t count = 0;
        for (size_t t = 0; t < W; ++t) {
            if (!market_open_mask[t][m]) ++count;
        }
        return count;
    }
};

// =======================================================================
// 6. Ring-Buffer Espejo Zero-Heap O(1) de Acceso Contiguo (2W)
// =======================================================================

template <typename T, size_t W>
class RingBufferEspejoCxx {
    static_assert(W >= 4, "W debe ser al menos 4.");

private:
    std::array<T, 2 * W> ring_{};
    size_t head_{0};
    size_t total_muestras_{0};

public:
    constexpr RingBufferEspejoCxx() noexcept = default;

    // Inserción en tiempo constante O(1) con escritura doble espejo
    void push(const T& val) noexcept {
        size_t pos = head_ % W;
        ring_[pos] = val;
        ring_[pos + W] = val;
        head_ = pos + 1;
        ++total_muestras_;
    }

    // Puntero directo al bloque contiguo de las W muestras más recientes (la más vieja en [0], la presente en [W-1])
    [[nodiscard]] const T* ventana_actual() const noexcept {
        return &ring_[head_ % W];
    }

    [[nodiscard]] constexpr size_t capacidad() const noexcept { return W; }
    [[nodiscard]] constexpr size_t total() const noexcept { return total_muestras_; }
    [[nodiscard]] constexpr bool llena() const noexcept { return total_muestras_ >= W; }

    void clear() noexcept {
        head_ = 0;
        total_muestras_ = 0;
        ring_.fill(T{});
    }
};

// =======================================================================
// 7. Cálculo constexpr de Divisores Exactos y Residuos Multiescala
// =======================================================================

template <size_t W>
[[nodiscard]] constexpr size_t calcular_divisores_exactos_cxx(
    size_t s_min,
    double max_ratio,
    std::array<size_t, 32>& out_divisores
) noexcept {
    size_t count = 0;
    size_t s_max = (max_ratio > 0.0) ? static_cast<size_t>(static_cast<double>(W) / max_ratio) : W;
    if (s_max < s_min) s_max = s_min;

    for (size_t s = s_min; s <= s_max && count < out_divisores.size(); ++s) {
        if (W % s == 0) {
            out_divisores[count++] = s;
        }
    }
    return count;
}

template <size_t W>
[[nodiscard]] constexpr size_t calcular_residuo_maximo_cxx(
    size_t s_min,
    double max_ratio,
    bool forzar_divisores_exactos
) noexcept {
    if (forzar_divisores_exactos) {
        return 0; // Por definición analítica, W % s == 0 para todo divisor exacto
    }
    size_t s_max = (max_ratio > 0.0) ? static_cast<size_t>(static_cast<double>(W) / max_ratio) : W;
    if (s_max < s_min) s_max = s_min;

    size_t r_max = 0;
    for (size_t s = s_min; s <= s_max; ++s) {
        size_t r = W % s;
        if (r > r_max) r_max = r;
    }
    return r_max;
}

} // namespace ark::ventanas
