// motor_ventanas.hpp -- MOTOR STREAMING MONO-ACTIVO ZERO-HEAP Y ADUANA AXIOMA 2 (ARK-WIN-01-CXX)
// Reactor ARK / Canción del Micelio -- ARQUITECTO DE DATOS & COMANDANTE
// Corregido por el Jefe de Arquitectos de Código (2026-11-01), acta AUDITORIA_CAUSAL/AUDITORIA_VENTANAS_CXX.md:
//   B4/E8 paso transaccional: primero se valida y se calcula todo; el estado solo cambia si el paso es admisible. Tras
//         cualquier error el motor queda exactamente como antes de la llamada.
//   B5    la ventana se sella con la estampa de DISPONIBILIDAD (R1) de sus filas; el cierre va aparte.
//   E3    las filas sin sesión (mercado cerrado) no cuentan como imputación: ni en el numerador ni en el denominador.
//   E4    el jitter no inventa huecos: omitidas = floor((dt + tol) / nominal) - 1.
//   E5    los huecos no dependen del reloj: los declara el calendario (sesiones_omitidas) o, si no, la cadencia nominal.
//   E7    una fila cerrada sin historia vale NaN.
//   E10   configuración validada: un motor con configuración inválida rechaza todos los pasos.
#pragma once

#include "ark/ventanas/contratos_ventanas.hpp"
#include "ark/ventanas/lazo_shannon.hpp"
#include "ark/ventanas/desacoplamiento_artefactos.hpp"
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <array>
#include <span>
#include <algorithm>
#include <cstring>
#include <limits>

namespace ark::ventanas {

// =======================================================================
// 1. Estado de Retorno del Motor en Cada Paso de Streaming
// =======================================================================

enum class ResultadoPaso : uint8_t {
    SIN_EMISION = 0,               // Muestra absorbida correctamente; no correspondía emitir ventana (stride S)
    VENTANA_EMITIDA = 1,           // Ventana causal W completada y emitida exitosamente
    ERROR_INVERSION_TEMPORAL = 2,  // Violación Axioma 2: t_k < t_{k-1} (paquete desordenado)
    ERROR_TIMESTAMP_DUPLICADO = 3, // Violación Axioma 2: t_k == t_{k-1} (colisión de eventos simultáneos)
    ERROR_COLISION_SUB_BARRA = 4,  // Violación Axioma 2: delta_t < delta_t_nominal - jitter
    ERROR_NAN_INICIAL = 5,         // Violación Axioma 1: NaN en k=0 sin dato previo (prohíbe bfill)
    ERROR_IMPUTACION_PROHIBIDA = 6,// Violación Axioma 1: Dato no finito bajo política PROHIBIR_INTERPOLACION
    ERROR_PRESUPUESTO_EXCEDIDO = 7,// Violación Axioma 2: Imputación + gaps > umbral bajo abortar_si_excede=true
    ERROR_CONFIGURACION_INVALIDA = 8, // E10: la configuración no pasa validar()
    ERROR_DISPONIBILIDAD_INVALIDA = 9,// B5: disponibilidad declarada anterior al cierre, o que retrocede
    FILA_DESCARTADA_INTERSECCION = 10 // INTERSECCION_ESTRICTA: la fila no tiene sesión en todos los activos (no es error)
};

[[nodiscard]] constexpr bool es_error(ResultadoPaso r) noexcept {
    return r != ResultadoPaso::SIN_EMISION && r != ResultadoPaso::VENTANA_EMITIDA && r != ResultadoPaso::FILA_DESCARTADA_INTERSECCION;
}

[[nodiscard]] constexpr std::string_view nombre_resultado_paso(ResultadoPaso r) noexcept {
    switch (r) {
        case ResultadoPaso::SIN_EMISION: return "SIN_EMISION";
        case ResultadoPaso::VENTANA_EMITIDA: return "VENTANA_EMITIDA";
        case ResultadoPaso::ERROR_INVERSION_TEMPORAL: return "ERROR_INVERSION_TEMPORAL";
        case ResultadoPaso::ERROR_TIMESTAMP_DUPLICADO: return "ERROR_TIMESTAMP_DUPLICADO";
        case ResultadoPaso::ERROR_COLISION_SUB_BARRA: return "ERROR_COLISION_SUB_BARRA";
        case ResultadoPaso::ERROR_NAN_INICIAL: return "ERROR_NAN_INICIAL";
        case ResultadoPaso::ERROR_IMPUTACION_PROHIBIDA: return "ERROR_IMPUTACION_PROHIBIDA";
        case ResultadoPaso::ERROR_PRESUPUESTO_EXCEDIDO: return "ERROR_PRESUPUESTO_EXCEDIDO";
        case ResultadoPaso::ERROR_CONFIGURACION_INVALIDA: return "ERROR_CONFIGURACION_INVALIDA";
        case ResultadoPaso::ERROR_DISPONIBILIDAD_INVALIDA: return "ERROR_DISPONIBILIDAD_INVALIDA";
        case ResultadoPaso::FILA_DESCARTADA_INTERSECCION: return "FILA_DESCARTADA_INTERSECCION";
    }
    return "DESCONOCIDO";
}

// Entrada de un paso mono-activo. Los campos opcionales con su valor por omisión reproducen la API anterior.
struct EntradaPasoCxx {
    uint64_t ts_ms{0};                 // sello de cierre de la barra (el de la rejilla)
    double valor{std::numeric_limits<double>::quiet_NaN()};
    bool mascara{false};               // la fuente ya la marca como imputada
    bool con_sesion{true};             // false = mercado cerrado en esta fila (no cuenta en el presupuesto, d_nu = 0)
    uint64_t disponible_ms{0};         // estampa de disponibilidad R1 (BarraCausal::estampa_disponible_ms); 0 = no declarada
    int32_t sesiones_omitidas{-1};     // sesiones del calendario sin fila entre la anterior y esta; -1 = deducir de la cadencia
    double volumen{0.0};
    double trades{0.0};
    double nocional{0.0};
};

// Ultima fila ADMITIDA por la aduana (lectura; 2026-11-04, Jefe de Arquitectos): la usa el motor de reloj intrinseco como compuerta
// fila a fila sin duplicar la aduana. No cambia ninguna salida del motor de pared.
struct FilaAdmitidaCxx {
    uint64_t ts_ms{0};               // sello de cierre (con el offset si se estampa al cierre)
    uint64_t disponible_ms{0};       // disponibilidad efectiva (la declarada, o el cierre)
    bool disponibilidad_declarada{false};
    double valor{std::numeric_limits<double>::quiet_NaN()};   // tras la imputacion causal (forward-fill o E7)
    bool mascara{false};             // imputada o sin sesion
    bool con_sesion{true};
    double d_nu{0.0};                // incremento del reloj intrinseco de esta fila
    uint32_t omitidas{0};            // barras omitidas entre la fila anterior y esta
};

namespace detalle {
// Barras omitidas entre dos sellos con cadencia nominal y tolerancia de jitter (E4). dt >= nominal - tol (lo comprueba el que llama).
[[nodiscard]] constexpr uint32_t omitidas_por_cadencia(int64_t dt_ms, int64_t dt_nom_ms, int64_t tol_ms) noexcept {
    const int64_t q = (dt_ms + tol_ms) / dt_nom_ms;
    return q > 1 ? static_cast<uint32_t>(q - 1) : 0u;
}
}  // namespace detalle

// =======================================================================
// 2. Motor de Ventanas Causales Mono-Activo CXX (Zero-Heap W)
// =======================================================================

template <size_t W>
class MotorVentanasCausalesCxx {
    static_assert(W >= 4, "lookback_W debe ser al menos 4 muestras.");

private:
    ConfiguracionVentanaCausalCxx<W> config_{};
    ErrorConfiguracion error_config_{ErrorConfiguracion::NINGUNO};

    // Búferes circulares espejo de memoria contigua 2W (Axioma 5: cero heap en streaming)
    RingBufferEspejoCxx<double, W> ring_datos_{};
    RingBufferEspejoCxx<bool, W> ring_mask_{};
    RingBufferEspejoCxx<bool, W> ring_sesion_{};
    RingBufferEspejoCxx<double, W> ring_dnu_{};
    RingBufferEspejoCxx<uint64_t, W> ring_disp_{};
    RingBufferEspejoCxx<bool, W> ring_disp_decl_{};
    RingBufferEspejoCxx<uint32_t, W> ring_gaps_{};

    // Lazo causal recursivo de Shannon (Axioma 4)
    LazoFeedbackShannonCxx<W> lazo_shannon_;

    // Catálogo estático multiescala de divisores y residuos medidos (Axioma 6)
    size_t num_divisores_{0};
    std::array<size_t, 32> divisores_{};
    size_t residuo_max_{0};
    double pct_descarte_{0.0};

    double ultimo_lookahead_score_{0.0};

    // Estado interno del cursor temporal
    size_t cursor_global_{0};
    size_t siguiente_emision_{W - 1};
    uint64_t ts_prev_{0};
    uint64_t disp_prev_{0};
    double ultimo_val_valido_{0.0};
    bool tiene_valor_previo_{false};
    double nu_acum_actual_{0.0};

    // Métricas de auditoría acumuladas
    uint32_t total_ventanas_emitidas_{0};
    uint32_t total_ventanas_invalidadas_{0};
    uint32_t total_gaps_detectados_{0};
    FilaAdmitidaCxx ultima_{};

public:
    explicit MotorVentanasCausalesCxx(const ConfiguracionVentanaCausalCxx<W>& config = {}) noexcept
        : config_(config),
          error_config_(config.validar()),
          lazo_shannon_(config_) {
        num_divisores_ = calcular_divisores_exactos_cxx<W>(config_.s_min, config_.max_ratio_escalas, divisores_);
        residuo_max_ = calcular_residuo_maximo_cxx<W>(config_.s_min, config_.max_ratio_escalas, config_.forzar_divisores_exactos);
        if (config_.forzar_divisores_exactos || config_.bidireccional_kantelhardt || residuo_max_ == 0) {
            pct_descarte_ = 0.0;
        } else {
            pct_descarte_ = (static_cast<double>(residuo_max_) / static_cast<double>(W)) * 100.0;
        }
    }

    [[nodiscard]] ErrorConfiguracion error_configuracion() const noexcept { return error_config_; }

    void registrar_auditoria_testigo(double score) noexcept { ultimo_lookahead_score_ = score; }

    // Paso completo. Transaccional: si devuelve un error, el estado no cambió.
    [[nodiscard]] ResultadoPaso paso(const EntradaPasoCxx& e, VentanaCausalCertificadaCxx<W>& out_ventana) noexcept {
        if (error_config_ != ErrorConfiguracion::NINGUNO) return ResultadoPaso::ERROR_CONFIGURACION_INVALIDA;

        // 1. Estampado causal al cierre si está configurado (Propuesta 3: False por defecto)
        const int64_t offset_cierre_ms = (config_.estampar_al_cierre && config_.delta_t_barra_segundos > 0)
            ? (config_.delta_t_barra_segundos * 1000) : 0;
        const uint64_t ts_k = e.ts_ms + static_cast<uint64_t>(offset_cierre_ms);

        // 2. Aduana de Axioma 2 (monotonía, resolución, huecos) -- sin tocar el estado
        uint32_t gap_barras_k = 0;
        if (cursor_global_ > 0) {
            if (ts_k < ts_prev_) return ResultadoPaso::ERROR_INVERSION_TEMPORAL;
            if (ts_k == ts_prev_) return ResultadoPaso::ERROR_TIMESTAMP_DUPLICADO;
            if (config_.delta_t_barra_segundos > 0) {
                const int64_t dt_ms = static_cast<int64_t>(ts_k - ts_prev_);
                const int64_t dt_nom_ms = config_.delta_t_barra_segundos * 1000;
                const int64_t tol_ms = static_cast<int64_t>(std::llround(config_.tolerancia_jitter_segundos * 1000.0));
                if (dt_ms < (dt_nom_ms - tol_ms)) return ResultadoPaso::ERROR_COLISION_SUB_BARRA;
                if (e.sesiones_omitidas < 0) gap_barras_k = detalle::omitidas_por_cadencia(dt_ms, dt_nom_ms, tol_ms);
            }
            if (e.sesiones_omitidas > 0) gap_barras_k = static_cast<uint32_t>(e.sesiones_omitidas);
        }

        // 3. Disponibilidad (B5): nunca antes del cierre ni hacia atrás
        const bool disp_decl = e.disponible_ms != 0;
        const uint64_t disp_k = disp_decl ? e.disponible_ms : ts_k;
        if (disp_decl && disp_k < ts_k) return ResultadoPaso::ERROR_DISPONIBILIDAD_INVALIDA;
        if (cursor_global_ > 0 && disp_k < disp_prev_) return ResultadoPaso::ERROR_DISPONIBILIDAD_INVALIDA;

        // 4. Aduana de Axioma 1: imputación causal pura sin bfill retrospectivo
        double val_k = e.valor;
        bool mask_k = e.mascara;
        bool actualiza_valido = false;
        if (!e.con_sesion) {
            val_k = tiene_valor_previo_ ? ultimo_val_valido_ : std::numeric_limits<double>::quiet_NaN();   // E7
            mask_k = true;
        } else if (!std::isfinite(e.valor)) {
            if (config_.metodo_imputacion == MetodoImputacion::PROHIBIR_INTERPOLACION) return ResultadoPaso::ERROR_IMPUTACION_PROHIBIDA;
            if (!tiene_valor_previo_) return ResultadoPaso::ERROR_NAN_INICIAL;   // Prohibido mirar al futuro (bfill)
            val_k = ultimo_val_valido_;
            mask_k = true;
        } else {
            if (mask_k && config_.metodo_imputacion == MetodoImputacion::PROHIBIR_INTERPOLACION) {
                return ResultadoPaso::ERROR_IMPUTACION_PROHIBIDA;
            }
            actualiza_valido = true;
        }

        // 5. Reloj intrínseco (Axioma 4), calculado sin avanzar el lazo
        double d_nu = 0.0;
        typename LazoFeedbackShannonCxx<W>::Propuesta prop{};
        const bool usa_lazo = config_.reloj == TipoRelojVentana::NU_INFO_SHANNON && e.con_sesion;
        if (!e.con_sesion) {
            d_nu = 0.0;   // congelamiento rígido en sesiones cerradas
        } else {
            switch (config_.reloj) {
                case TipoRelojVentana::PARED:
                    d_nu = (config_.congelar_en_huecos && mask_k) ? 0.0 : 1.0;
                    break;
                case TipoRelojVentana::NU_INFO_SHANNON:
                    prop = lazo_shannon_.calcular(val_k, mask_k);
                    d_nu = prop.d_nu;
                    break;
                case TipoRelojVentana::NU_VOLUMEN:
                    d_nu = (config_.congelar_en_huecos && mask_k) ? 0.0 : std::max(e.volumen, 0.0);
                    break;
                case TipoRelojVentana::NU_TRADES:
                    d_nu = (config_.congelar_en_huecos && mask_k) ? 0.0 : std::max(e.trades, 0.0);
                    break;
                case TipoRelojVentana::NU_NOCIONAL: {
                    const double noc = (e.nocional > 0.0) ? e.nocional : (val_k * e.volumen);
                    d_nu = (config_.congelar_en_huecos && mask_k) ? 0.0 : std::max(noc, 0.0);
                    break;
                }
            }
        }

        // 6. Presupuesto de la ventana que se emitiría (E3), con la fila nueva y sin empujarla todavía.
        //    ventana_actual() de los búferes viejos: [0] es la fila que sale; [1..W-1] siguen.
        const bool emite = cursor_global_ >= siguiente_emision_;
        size_t imputadas = 0, con_sesion = 0;
        uint32_t omitidas = 0;
        double ratio = 0.0;
        bool es_valida = true;
        if (emite) {
            const bool* mp = ring_mask_.ventana_actual();
            const bool* sp = ring_sesion_.ventana_actual();
            const uint32_t* gp = ring_gaps_.ventana_actual();
            for (size_t i = 1; i < W; ++i) {
                con_sesion += sp[i] ? 1u : 0u;
                imputadas += (sp[i] && mp[i]) ? 1u : 0u;
            }
            for (size_t i = 2; i < W; ++i) omitidas += gp[i];   // el hueco de la fila más vieja es anterior a la ventana
            con_sesion += e.con_sesion ? 1u : 0u;
            imputadas += (e.con_sesion && mask_k) ? 1u : 0u;
            omitidas += gap_barras_k;
            const double den = static_cast<double>(con_sesion + omitidas);
            ratio = den > 0.0 ? static_cast<double>(imputadas + omitidas) / den : std::numeric_limits<double>::quiet_NaN();
            if (!(ratio <= config_.umbral_imputacion_maximo)) {   // NaN (ventana sin ninguna sesión) tampoco es válida
                if (config_.abortar_si_excede_imputacion) return ResultadoPaso::ERROR_PRESUPUESTO_EXCEDIDO;
                es_valida = false;
            }
        }

        // 7. Confirmación: a partir de aquí el paso es admisible y el estado avanza
        ts_prev_ = ts_k;
        disp_prev_ = disp_k;
        if (actualiza_valido) {
            ultimo_val_valido_ = val_k;
            tiene_valor_previo_ = true;
        }
        if (usa_lazo) lazo_shannon_.confirmar(val_k, mask_k, prop);
        total_gaps_detectados_ += gap_barras_k;
        nu_acum_actual_ += d_nu;
        ring_datos_.push(val_k);
        ring_mask_.push(mask_k);
        ring_sesion_.push(e.con_sesion);
        ring_dnu_.push(d_nu);
        ring_disp_.push(disp_k);
        ring_disp_decl_.push(disp_decl);
        ring_gaps_.push(gap_barras_k);
        ultima_ = FilaAdmitidaCxx{ts_k, disp_k, disp_decl, val_k, mask_k, e.con_sesion, d_nu, gap_barras_k};

        if (!emite) {
            ++cursor_global_;
            return ResultadoPaso::SIN_EMISION;
        }

        // 8. Emisión de la ventana certificada
        if (!es_valida) ++total_ventanas_invalidadas_;
        const uint64_t* dp = ring_disp_.ventana_actual();
        const bool* ddp = ring_disp_decl_.ventana_actual();
        uint64_t t_disp = 0;
        bool todas_decl = true;
        for (size_t i = 0; i < W; ++i) {
            t_disp = std::max(t_disp, dp[i]);
            todas_decl = todas_decl && ddp[i];
        }
        out_ventana.t_actual = t_disp;
        out_ventana.t_cierre_ms = ts_k;
        out_ventana.disponibilidad_declarada = todas_decl;
        std::memcpy(out_ventana.valores.data(), ring_datos_.ventana_actual(), W * sizeof(double));
        std::memcpy(out_ventana.mascara_huecos.data(), ring_mask_.ventana_actual(), W * sizeof(bool));
        std::memcpy(out_ventana.con_sesion.data(), ring_sesion_.ventana_actual(), W * sizeof(bool));

        const double* dnu_ptr = ring_dnu_.ventana_actual();
        double delta_reloj = 0.0;
        for (size_t i = 0; i < W; ++i) delta_reloj += dnu_ptr[i];
        out_ventana.delta_reloj_ventana = delta_reloj;
        out_ventana.reloj_acumulado = nu_acum_actual_;

        out_ventana.num_divisores_exactos = num_divisores_;
        out_ventana.divisores_exactos = divisores_;
        out_ventana.residuo_temporal_maximo = residuo_max_;
        out_ventana.porcentaje_descarte_causal = pct_descarte_;
        out_ventana.lookahead_score = ultimo_lookahead_score_;
        out_ventana.es_valida = es_valida;
        out_ventana.ratio_imputacion = ratio;
        out_ventana.barras_omitidas = omitidas;
        out_ventana.tipo_reloj = config_.reloj;

        siguiente_emision_ += config_.stride_S;
        ++total_ventanas_emitidas_;
        ++cursor_global_;
        return ResultadoPaso::VENTANA_EMITIDA;
    }

    // API anterior (posicional): fila con sesión, disponibilidad no declarada, huecos por cadencia.
    [[nodiscard]] ResultadoPaso paso(uint64_t ts_ms, double val_raw, bool mask_raw, VentanaCausalCertificadaCxx<W>& out_ventana,
                                     double volumen = 0.0, double trades = 0.0, double nocional = 0.0) noexcept {
        EntradaPasoCxx e{};
        e.ts_ms = ts_ms;
        e.valor = val_raw;
        e.mascara = mask_raw;
        e.volumen = volumen;
        e.trades = trades;
        e.nocional = nocional;
        return paso(e, out_ventana);
    }

    // Procesamiento causal de arrays planos: la fila k solo ve la entrada k.
    template <typename Callback>
    size_t procesar_lote(const EntradaPasoCxx* entradas, size_t N, Callback&& on_ventana) noexcept {
        if (!entradas) return 0;
        VentanaCausalCertificadaCxx<W> win{};
        size_t emitidas = 0;
        for (size_t k = 0; k < N; ++k) {
            const ResultadoPaso res = paso(entradas[k], win);
            if (res == ResultadoPaso::VENTANA_EMITIDA) {
                on_ventana(win);
                ++emitidas;
            } else if (es_error(res) && config_.abortar_si_excede_imputacion) {
                break;
            }
        }
        return emitidas;
    }

    template <typename Callback>
    size_t procesar_lote(const uint64_t* timestamps_ms, const double* precios, const bool* mascaras, size_t N,
                         Callback&& on_ventana) noexcept {
        if (!timestamps_ms || !precios) return 0;
        VentanaCausalCertificadaCxx<W> win{};
        size_t emitidas = 0;
        for (size_t k = 0; k < N; ++k) {
            const bool m = (mascaras) ? mascaras[k] : false;
            const ResultadoPaso res = paso(timestamps_ms[k], precios[k], m, win);
            if (res == ResultadoPaso::VENTANA_EMITIDA) {
                on_ventana(win);
                ++emitidas;
            } else if (es_error(res) && config_.abortar_si_excede_imputacion) {
                break;
            }
        }
        return emitidas;
    }

    // Telemetría y Contadores de Auditoría
    [[nodiscard]] constexpr size_t cursor_global() const noexcept { return cursor_global_; }
    [[nodiscard]] constexpr uint32_t ventanas_emitidas() const noexcept { return total_ventanas_emitidas_; }
    [[nodiscard]] constexpr uint32_t ventanas_invalidadas() const noexcept { return total_ventanas_invalidadas_; }
    [[nodiscard]] constexpr uint32_t gaps_detectados() const noexcept { return total_gaps_detectados_; }
    [[nodiscard]] constexpr double nu_acumulado() const noexcept { return nu_acum_actual_; }
    [[nodiscard]] const FilaAdmitidaCxx& ultima_fila() const noexcept { return ultima_; }   // valida solo tras un paso admitido
    [[nodiscard]] const LazoFeedbackShannonCxx<W>& lazo_shannon() const noexcept { return lazo_shannon_; }
    [[nodiscard]] const ConfiguracionVentanaCausalCxx<W>& configuracion() const noexcept { return config_; }

    void reset() noexcept {
        ring_datos_.clear();
        ring_mask_.clear();
        ring_sesion_.clear();
        ring_dnu_.clear();
        ring_disp_.clear();
        ring_disp_decl_.clear();
        ring_gaps_.clear();
        lazo_shannon_.reset();
        ultimo_lookahead_score_ = 0.0;
        cursor_global_ = 0;
        siguiente_emision_ = W - 1;
        ts_prev_ = 0;
        disp_prev_ = 0;
        ultimo_val_valido_ = 0.0;
        tiene_valor_previo_ = false;
        nu_acum_actual_ = 0.0;
        total_ventanas_emitidas_ = 0;
        total_ventanas_invalidadas_ = 0;
        total_gaps_detectados_ = 0;
        ultima_ = FilaAdmitidaCxx{};
    }
};

} // namespace ark::ventanas
