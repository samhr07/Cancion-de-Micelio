// motor_multiactivo.hpp -- MOTOR STREAMING MULTIACTIVO ZERO-HEAP Y GUARDA DE PARES (ARK-WIN-01-CXX)
// Reactor ARK / Canción del Micelio -- ARQUITECTO DE DATOS & COMANDANTE
// Corregido por el Jefe de Arquitectos de Código (2026-11-01), acta AUDITORIA_CAUSAL/AUDITORIA_VENTANAS_CXX.md:
//   B3  guarda de pares con las horas de cierre reales: la arista u -> v es EJECUTABLE solo si la información de u (el
//       intervalo de cierre de la barra cuyo valor lleva la fila) termina, más la latencia, antes de que empiece el retorno
//       de v (el cierre de su barra real anterior). Es clasificar_ventanas de ventanas_informacion / ark::dataset. Sin horas
//       declaradas: NO_VERIFICADO. Se quita la regla «fuente cerrada => REZAGO»: el último cierre de una fuente cerrada ya es
//       pasado. Antes, con los dos abiertos, todo salía EJECUTABLE (BTC -> SPX el mismo día es REZAGO: 246 de 246 días de 2025).
//   B4  fila transaccional: se validan las M entradas y se calcula todo antes de tocar el estado; tras un error no cambia nada.
//   B5  la ventana se sella con la mayor estampa de disponibilidad de sus celdas.
//   E3  presupuesto sobre las celdas con sesión: una celda cerrada no es una imputación.
//   E4/E5 huecos por calendario (sesiones_omitidas_rejilla) o por cadencia con tolerancia, con cualquier reloj.
//   E6  INTERSECCION_ESTRICTA implementada: la fila que no tiene sesión en todos los activos se descarta sin tocar el estado.
//   E7  activo cerrado sin historia: NaN.
//   E10 configuración validada, con el tick de cada activo cuando el reloj es de Shannon (E1).
#pragma once

#include "ark/ventanas/contratos_ventanas.hpp"
#include "ark/ventanas/lazo_shannon.hpp"
#include "ark/ventanas/motor_ventanas.hpp"
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <array>
#include <span>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <limits>
#include <string_view>

namespace ark::ventanas {

// Núcleo de la guarda (idéntico a ventanas_informacion.clasificar_ventanas y a ark::dataset::clasificar_ventanas):
// información de u en [info_lo, info_hi], retorno de v que empieza en [ini_lo, ini_hi], latencia L (s), en ms UTC.
//   info_hi + L <  ini_lo -> EJECUTABLE;  info_lo + L >= ini_hi -> REZAGO;  si no, NO_VERIFICADO.
// Se compara en nanosegundos enteros, como pd.Timedelta, para la misma frontera. Margen = ini_lo - (info_hi + L), en s.
struct GuardaVentanasCxx {
    ClasePar clase{ClasePar::NO_VERIFICADO};
    double margen_s{std::numeric_limits<double>::quiet_NaN()};
};

[[nodiscard]] inline GuardaVentanasCxx clasificar_ventanas_cxx(int64_t info_lo_ms, int64_t info_hi_ms, int64_t ini_lo_ms,
                                                               int64_t ini_hi_ms, double latencia_s) noexcept {
    GuardaVentanasCxx r;
    if (!(latencia_s >= 0.0)) return r;
    const int64_t lat_ns = static_cast<int64_t>(std::llround(latencia_s * 1e9));
    const int64_t ilo = info_lo_ms * 1000000, ihi = info_hi_ms * 1000000, slo = ini_lo_ms * 1000000, shi = ini_hi_ms * 1000000;
    if (ihi + lat_ns < slo) r.clase = ClasePar::EJECUTABLE;
    else if (ilo + lat_ns >= shi) r.clase = ClasePar::REZAGO;
    else r.clase = ClasePar::NO_VERIFICADO;
    r.margen_s = static_cast<double>(slo - (ihi + lat_ns)) / 1e9;
    return r;
}

// Entrada de un activo en una fila de la rejilla.
struct EntradaActivoCxx {
    double valor{std::numeric_limits<double>::quiet_NaN()};
    bool mascara{false};
    bool con_sesion{true};         // market_open: false = cerrado en esta fila
    uint64_t disponible_ms{0};     // estampa de disponibilidad R1 de la fila (0 = no declarada)
    int64_t cierre_lo_ms{0};       // intervalo del cierre de la barra de esta fila (BarraCausal::cierre_lo/hi_utc_ms);
    int64_t cierre_hi_ms{0};       // 0 = no declarado (la guarda queda NO_VERIFICADO)
    double volumen{0.0};
    double trades{0.0};
    double nocional{0.0};
};

// =======================================================================
// Motor de Ventanas Causales Multiactivo CXX (Zero-Heap W x M)
// =======================================================================

template <size_t W, size_t M>
class MotorVentanasMultiactivoCxx {
    static_assert(W >= 4, "W debe ser al menos 4 muestras.");
    static_assert(M >= 1 && M <= 16, "M debe estar entre 1 y 16 activos para stack zero-heap.");

private:
    ConfiguracionVentanaCausalCxx<W> config_{};
    ErrorConfiguracion error_config_{ErrorConfiguracion::NINGUNO};
    std::array<double, M> ticks_{};

    // Búferes circulares espejo zero-heap por cada activo (Axioma 5: cero heap en streaming)
    std::array<RingBufferEspejoCxx<double, W>, M> ring_datos_{};
    std::array<RingBufferEspejoCxx<bool, W>, M> ring_mask_{};
    std::array<RingBufferEspejoCxx<bool, W>, M> ring_market_open_{};
    std::array<RingBufferEspejoCxx<double, W>, M> ring_dnu_{};
    std::array<RingBufferEspejoCxx<uint64_t, W>, M> ring_disp_{};
    std::array<RingBufferEspejoCxx<bool, W>, M> ring_disp_decl_{};

    // Eje temporal unificado de la rejilla
    RingBufferEspejoCxx<uint64_t, W> ring_ts_{};
    RingBufferEspejoCxx<uint32_t, W> ring_gaps_{};

    std::array<LazoFeedbackShannonCxx<W>, M> lazo_shannon_{};

    std::array<double, M> nu_acum_por_activo_{};
    std::array<double, M> ultimo_val_valido_{};
    std::array<bool, M> tiene_valor_previo_{};
    std::array<uint64_t, M> disp_prev_{};

    // Guarda: intervalo de cierre de la última barra real de cada activo (0/0 = no declarado)
    std::array<int64_t, M> ult_cierre_lo_{};
    std::array<int64_t, M> ult_cierre_hi_{};

    // Restricciones externas (p. ej. verificadas con Helmholtz-Hodge): solo pueden ENDURECER la guarda, nunca dar EJECUTABLE.
    std::array<std::array<ClasePar, M>, M> restriccion_par_{};
    std::array<std::array<bool, M>, M> hay_restriccion_{};

    char simbolos_[M][16]{};

    size_t num_divisores_{0};
    std::array<size_t, 32> divisores_{};
    size_t residuo_max_{0};
    double pct_descarte_{0.0};

    double ultimo_lookahead_score_{0.0};

    size_t cursor_global_{0};
    size_t siguiente_emision_{W - 1};
    uint64_t ts_prev_{0};

    uint32_t total_ventanas_emitidas_{0};
    uint32_t total_ventanas_invalidadas_{0};
    uint32_t total_gaps_detectados_{0};
    uint32_t total_advertencias_rezago_{0};
    uint32_t total_asimetrias_union_{0};
    uint32_t total_filas_descartadas_{0};

    void revalidar() noexcept {
        error_config_ = config_.validar(false);
        if (error_config_ != ErrorConfiguracion::NINGUNO) return;
        for (size_t m = 0; m < M; ++m) {
            const double t = ticks_[m];
            if (!std::isnan(t) && !(t > 0.0 && std::isfinite(t))) { error_config_ = ErrorConfiguracion::TICK_INVALIDO; return; }
            if (config_.reloj == TipoRelojVentana::NU_INFO_SHANNON && std::isnan(t)) {
                error_config_ = ErrorConfiguracion::TICK_NO_DECLARADO;
                return;
            }
        }
    }

public:
    explicit MotorVentanasMultiactivoCxx(const ConfiguracionVentanaCausalCxx<W>& config = {}) noexcept
        : config_(config) {
        for (size_t m = 0; m < M; ++m) {
            ticks_[m] = config_.tick_size;
            lazo_shannon_[m].configurar(config_, ticks_[m]);
            std::snprintf(simbolos_[m], sizeof(simbolos_[m]), "ACT_%zu", m);
        }
        revalidar();
        num_divisores_ = calcular_divisores_exactos_cxx<W>(config_.s_min, config_.max_ratio_escalas, divisores_);
        residuo_max_ = calcular_residuo_maximo_cxx<W>(config_.s_min, config_.max_ratio_escalas, config_.forzar_divisores_exactos);
        if (config_.forzar_divisores_exactos || config_.bidireccional_kantelhardt || residuo_max_ == 0) {
            pct_descarte_ = 0.0;
        } else {
            pct_descarte_ = (static_cast<double>(residuo_max_) / static_cast<double>(W)) * 100.0;
        }
    }

    [[nodiscard]] ErrorConfiguracion error_configuracion() const noexcept { return error_config_; }

    void establecer_simbolo(size_t m, const char* sym) noexcept {
        if (m < M && sym) {
            size_t len = 0;
            while (sym[len] != '\0' && len < sizeof(simbolos_[m]) - 1) {
                simbolos_[m][len] = sym[len];
                ++len;
            }
            simbolos_[m][len] = '\0';
        }
    }

    // Tick del activo m (para su suelo dinámico de sigma). Debe llamarse antes del primer paso.
    void establecer_tick(size_t m, double tick) noexcept {
        if (m >= M) return;
        ticks_[m] = tick;
        lazo_shannon_[m].configurar(config_, tick);
        revalidar();
    }

    // Restricción externa sobre la arista u -> v: REZAGO o NO_VERIFICADO endurecen lo que diga la guarda; EJECUTABLE y
    // SIN_SESION se ignoran (una afirmación externa no puede convertir en ejecutable lo que el reloj dice que no lo es).
    void establecer_clasificacion_par(size_t u, size_t v, ClasePar clase) noexcept {
        if (u < M && v < M) {
            const bool endurece = clase == ClasePar::REZAGO || clase == ClasePar::NO_VERIFICADO;
            hay_restriccion_[u][v] = endurece;
            restriccion_par_[u][v] = clase;
        }
    }

    void registrar_auditoria_testigo(double score) noexcept { ultimo_lookahead_score_ = score; }

    // Paso de una fila. Transaccional: si devuelve error (o descarta la fila), el estado no cambió.
    // sesiones_omitidas_rejilla: filas de la rejilla que faltan entre la anterior y esta según el calendario (-1 = cadencia).
    [[nodiscard]] ResultadoPaso paso(uint64_t ts_ms, const EntradaActivoCxx* entradas,
                                     VentanaCausalCertificadaMultiactivoCxx<W, M>& out_ventana,
                                     int32_t sesiones_omitidas_rejilla = -1) noexcept {
        if (error_config_ != ErrorConfiguracion::NINGUNO) return ResultadoPaso::ERROR_CONFIGURACION_INVALIDA;
        if (!entradas) return ResultadoPaso::ERROR_NAN_INICIAL;

        // 0. Intersección estricta: solo las filas con sesión en todos los activos son de la rejilla (E6)
        if (config_.politica_sincronizacion == PoliticaSincronizacion::INTERSECCION_ESTRICTA) {
            for (size_t m = 0; m < M; ++m) {
                if (!entradas[m].con_sesion) {
                    ++total_filas_descartadas_;
                    return ResultadoPaso::FILA_DESCARTADA_INTERSECCION;
                }
            }
        }

        const int64_t offset_cierre_ms = (config_.estampar_al_cierre && config_.delta_t_barra_segundos > 0)
            ? (config_.delta_t_barra_segundos * 1000) : 0;
        const uint64_t ts_k = ts_ms + static_cast<uint64_t>(offset_cierre_ms);

        // 1. Aduana de Axioma 2 en la rejilla maestra
        uint32_t gap_barras_k = 0;
        if (cursor_global_ > 0) {
            if (ts_k < ts_prev_) return ResultadoPaso::ERROR_INVERSION_TEMPORAL;
            if (ts_k == ts_prev_) return ResultadoPaso::ERROR_TIMESTAMP_DUPLICADO;
            if (config_.delta_t_barra_segundos > 0) {
                const int64_t dt_ms = static_cast<int64_t>(ts_k - ts_prev_);
                const int64_t dt_nom_ms = config_.delta_t_barra_segundos * 1000;
                const int64_t tol_ms = static_cast<int64_t>(std::llround(config_.tolerancia_jitter_segundos * 1000.0));
                if (dt_ms < (dt_nom_ms - tol_ms)) return ResultadoPaso::ERROR_COLISION_SUB_BARRA;
                if (sesiones_omitidas_rejilla < 0) gap_barras_k = detalle::omitidas_por_cadencia(dt_ms, dt_nom_ms, tol_ms);
            }
            if (sesiones_omitidas_rejilla > 0) gap_barras_k = static_cast<uint32_t>(sesiones_omitidas_rejilla);
        }

        // 2. Validación y cálculo de las M celdas, sin tocar nada
        std::array<double, M> val{};
        std::array<bool, M> mask{};
        std::array<bool, M> real{};
        std::array<uint64_t, M> disp{};
        std::array<bool, M> disp_decl{};
        std::array<double, M> dnu{};
        std::array<typename LazoFeedbackShannonCxx<W>::Propuesta, M> prop{};
        for (size_t m = 0; m < M; ++m) {
            const EntradaActivoCxx& e = entradas[m];
            disp_decl[m] = e.disponible_ms != 0;
            disp[m] = disp_decl[m] ? e.disponible_ms : ts_k;
            if (disp_decl[m] && disp[m] < ts_k) return ResultadoPaso::ERROR_DISPONIBILIDAD_INVALIDA;
            if (cursor_global_ > 0 && disp[m] < disp_prev_[m]) return ResultadoPaso::ERROR_DISPONIBILIDAD_INVALIDA;

            val[m] = e.valor;
            mask[m] = e.mascara;
            real[m] = false;
            if (!e.con_sesion) {
                // Cerrado: arrastre del último cierre (ffill causal) con máscara; sin historia, NaN (E7)
                val[m] = tiene_valor_previo_[m] ? ultimo_val_valido_[m] : std::numeric_limits<double>::quiet_NaN();
                mask[m] = true;
            } else if (!std::isfinite(e.valor)) {
                if (config_.metodo_imputacion == MetodoImputacion::PROHIBIR_INTERPOLACION ||
                    config_.politica_sincronizacion == PoliticaSincronizacion::INTERSECCION_ESTRICTA) {
                    return ResultadoPaso::ERROR_IMPUTACION_PROHIBIDA;   // la intersección prohíbe NaN, como el Python
                }
                if (!tiene_valor_previo_[m]) return ResultadoPaso::ERROR_NAN_INICIAL;   // prohibido el bfill
                val[m] = ultimo_val_valido_[m];
                mask[m] = true;
            } else {
                if (mask[m] && config_.metodo_imputacion == MetodoImputacion::PROHIBIR_INTERPOLACION) {
                    return ResultadoPaso::ERROR_IMPUTACION_PROHIBIDA;
                }
                real[m] = true;
            }

            // Tiempo propio d_nu (Axioma 4 / 3b): congelado si cerrado o (con congelar_en_huecos) imputado
            if (!e.con_sesion) {
                dnu[m] = 0.0;
            } else if (config_.reloj == TipoRelojVentana::NU_INFO_SHANNON) {
                // Como el Python: el lazo ve toda celda con sesión, con su máscara (una imputada da d_nu = 0 y no entra
                // en la historia de retornos).
                prop[m] = lazo_shannon_[m].calcular(val[m], mask[m]);
                dnu[m] = prop[m].d_nu;
            } else if (config_.congelar_en_huecos && mask[m]) {
                dnu[m] = 0.0;
            } else {
                switch (config_.reloj) {
                    case TipoRelojVentana::PARED: dnu[m] = 1.0; break;
                    case TipoRelojVentana::NU_VOLUMEN: dnu[m] = std::max(e.volumen, 0.0); break;
                    case TipoRelojVentana::NU_TRADES: dnu[m] = std::max(e.trades, 0.0); break;
                    case TipoRelojVentana::NU_NOCIONAL: {
                        const double noc = (e.nocional > 0.0) ? e.nocional : (val[m] * e.volumen);
                        dnu[m] = std::max(noc, 0.0);
                        break;
                    }
                    default: dnu[m] = 0.0; break;
                }
            }
        }
        // 3. Presupuesto de la ventana que se emitiría (E3), sin empujar la fila
        const bool emite = cursor_global_ >= siguiente_emision_;
        size_t con_sesion = 0, imputadas = 0;
        uint32_t omitidas = 0;
        double ratio = 0.0;
        bool es_valida = true;
        if (emite) {
            for (size_t m = 0; m < M; ++m) {
                const bool* mp = ring_mask_[m].ventana_actual();
                const bool* op = ring_market_open_[m].ventana_actual();
                for (size_t i = 1; i < W; ++i) {
                    con_sesion += op[i] ? 1u : 0u;
                    imputadas += (op[i] && mp[i]) ? 1u : 0u;
                }
                con_sesion += entradas[m].con_sesion ? 1u : 0u;
                imputadas += (entradas[m].con_sesion && mask[m]) ? 1u : 0u;
            }
            const uint32_t* gp = ring_gaps_.ventana_actual();
            for (size_t i = 2; i < W; ++i) omitidas += gp[i];
            omitidas += gap_barras_k;
            const double den = static_cast<double>(con_sesion) + static_cast<double>(omitidas) * static_cast<double>(M);
            const double num = static_cast<double>(imputadas) + static_cast<double>(omitidas) * static_cast<double>(M);
            ratio = den > 0.0 ? num / den : std::numeric_limits<double>::quiet_NaN();
            if (!(ratio <= config_.umbral_imputacion_maximo)) {
                if (config_.abortar_si_excede_imputacion) return ResultadoPaso::ERROR_PRESUPUESTO_EXCEDIDO;
                es_valida = false;
            }
        }

        // 4. Guarda de pares sobre la fila nueva (antes de actualizar los cierres)
        std::array<int64_t, M> info_lo{}, info_hi{};
        for (size_t m = 0; m < M; ++m) {
            const bool cierre_decl = entradas[m].cierre_hi_ms != 0;
            if (real[m] && cierre_decl) {
                info_lo[m] = entradas[m].cierre_lo_ms;
                info_hi[m] = entradas[m].cierre_hi_ms;
            } else {
                info_lo[m] = ult_cierre_lo_[m];
                info_hi[m] = ult_cierre_hi_[m];
            }
        }

        // 5. Confirmación
        ts_prev_ = ts_k;
        total_gaps_detectados_ += gap_barras_k;
        for (size_t m = 0; m < M; ++m) {
            disp_prev_[m] = disp[m];
            if (real[m]) {
                ultimo_val_valido_[m] = val[m];
                tiene_valor_previo_[m] = true;
            }
            if (entradas[m].con_sesion && config_.reloj == TipoRelojVentana::NU_INFO_SHANNON) {
                lazo_shannon_[m].confirmar(val[m], mask[m], prop[m]);
            }
            nu_acum_por_activo_[m] += dnu[m];
            ring_datos_[m].push(val[m]);
            ring_mask_[m].push(mask[m]);
            ring_market_open_[m].push(entradas[m].con_sesion);
            ring_dnu_[m].push(dnu[m]);
            ring_disp_[m].push(disp[m]);
            ring_disp_decl_[m].push(disp_decl[m]);
        }
        ring_ts_.push(ts_k);
        ring_gaps_.push(gap_barras_k);

        if (!emite) {
            // Los cierres se actualizan igual: la guarda de la próxima fila los necesita
            for (size_t m = 0; m < M; ++m) {
                if (real[m] && entradas[m].cierre_hi_ms != 0) {
                    ult_cierre_lo_[m] = entradas[m].cierre_lo_ms;
                    ult_cierre_hi_[m] = entradas[m].cierre_hi_ms;
                }
            }
            ++cursor_global_;
            return ResultadoPaso::SIN_EMISION;
        }

        // 6. Emisión
        if (!es_valida) ++total_ventanas_invalidadas_;
        out_ventana.t_cierre_ms = ts_k;
        out_ventana.ratio_imputacion_global = ratio;
        out_ventana.celdas_con_sesion = con_sesion;
        out_ventana.celdas_imputadas_con_sesion = imputadas;
        out_ventana.es_valida = es_valida;
        out_ventana.barras_omitidas = omitidas;
        out_ventana.tipo_reloj = config_.reloj;
        out_ventana.lookahead_score = ultimo_lookahead_score_;
        out_ventana.num_divisores_exactos = num_divisores_;
        out_ventana.divisores_exactos = divisores_;
        out_ventana.residuo_temporal_maximo = residuo_max_;
        out_ventana.porcentaje_descarte_causal = pct_descarte_;
        for (size_t m = 0; m < M; ++m) std::memcpy(out_ventana.simbolos[m], simbolos_[m], sizeof(simbolos_[m]));

        uint64_t t_disp = 0;
        bool todas_decl = true;
        double sum_delta_reloj_total = 0.0;
        double sum_reloj_acum_total = 0.0;
        for (size_t m = 0; m < M; ++m) {
            const double* datos_m = ring_datos_[m].ventana_actual();
            const bool* mask_m = ring_mask_[m].ventana_actual();
            const bool* open_m = ring_market_open_[m].ventana_actual();
            const double* dnu_m = ring_dnu_[m].ventana_actual();
            const uint64_t* disp_m = ring_disp_[m].ventana_actual();
            const bool* decl_m = ring_disp_decl_[m].ventana_actual();
            double delta_nu_m = 0.0;
            for (size_t t = 0; t < W; ++t) {
                out_ventana.valores[t][m] = datos_m[t];
                out_ventana.mascara_huecos[t][m] = mask_m[t];
                out_ventana.market_open_mask[t][m] = open_m[t];
                delta_nu_m += dnu_m[t];
                t_disp = std::max(t_disp, disp_m[t]);
                todas_decl = todas_decl && decl_m[t];
            }
            out_ventana.delta_reloj_por_activo[m] = delta_nu_m;
            out_ventana.reloj_acum_por_activo[m] = nu_acum_por_activo_[m];
            sum_delta_reloj_total += delta_nu_m;
            sum_reloj_acum_total += nu_acum_por_activo_[m];
        }
        out_ventana.t_actual = t_disp;
        out_ventana.disponibilidad_declarada = todas_decl;
        out_ventana.delta_reloj_ventana = sum_delta_reloj_total / static_cast<double>(M);
        out_ventana.reloj_acumulado = sum_reloj_acum_total / static_cast<double>(M);

        // Guarda de sincronización causal de pares dirigidos (Axioma 3 / 3a), sobre la última fila
        for (size_t u = 0; u < M; ++u) {
            for (size_t v = 0; v < M; ++v) {
                if (u == v) {
                    out_ventana.clasificacion_pares[u][v] = ClasePar::EJECUTABLE;
                    out_ventana.margen_pares_s[u][v] = std::numeric_limits<double>::quiet_NaN();
                    continue;
                }
                GuardaVentanasCxx g{};
                if (!entradas[v].con_sesion) {
                    g.clase = ClasePar::SIN_SESION;
                } else if (info_hi[u] != 0 && ult_cierre_hi_[v] != 0) {
                    // el retorno de v en esta fila empieza en el cierre de su barra real anterior (ult_cierre_* de v)
                    g = clasificar_ventanas_cxx(info_lo[u], info_hi[u], ult_cierre_lo_[v], ult_cierre_hi_[v],
                                                config_.latencia_guarda_s);
                }
                if (hay_restriccion_[u][v] && g.clase == ClasePar::EJECUTABLE) g.clase = restriccion_par_[u][v];
                out_ventana.clasificacion_pares[u][v] = g.clase;
                out_ventana.margen_pares_s[u][v] = g.margen_s;
            }
        }
        for (size_t m = 0; m < M; ++m) {
            if (real[m] && entradas[m].cierre_hi_ms != 0) {
                ult_cierre_lo_[m] = entradas[m].cierre_lo_ms;
                ult_cierre_hi_[m] = entradas[m].cierre_hi_ms;
            }
        }

        if (out_ventana.tiene_rezago_causal()) ++total_advertencias_rezago_;
        if (out_ventana.tiene_asimetria_union()) ++total_asimetrias_union_;

        siguiente_emision_ += config_.stride_S;
        ++total_ventanas_emitidas_;
        ++cursor_global_;
        return ResultadoPaso::VENTANA_EMITIDA;
    }

    // API anterior: precios, máscaras y mercado abierto por activo; sin horas de cierre (la guarda da NO_VERIFICADO) ni
    // disponibilidad declarada.
    [[nodiscard]] ResultadoPaso paso(uint64_t ts_ms, const double* precios_raw, const bool* mascara_raw,
                                     const bool* market_open_raw, VentanaCausalCertificadaMultiactivoCxx<W, M>& out_ventana,
                                     const double* volumenes = nullptr, const double* trades = nullptr,
                                     const double* nocionales = nullptr) noexcept {
        if (!precios_raw || !market_open_raw) return ResultadoPaso::ERROR_NAN_INICIAL;
        std::array<EntradaActivoCxx, M> e{};
        for (size_t m = 0; m < M; ++m) {
            e[m].valor = precios_raw[m];
            e[m].mascara = mascara_raw ? mascara_raw[m] : false;
            e[m].con_sesion = market_open_raw[m];
            e[m].volumen = volumenes ? volumenes[m] : 0.0;
            e[m].trades = trades ? trades[m] : 0.0;
            e[m].nocional = nocionales ? nocionales[m] : 0.0;
        }
        return paso(ts_ms, e.data(), out_ventana);
    }

    [[nodiscard]] ResultadoPaso paso(uint64_t ts_ms, const std::array<double, M>& precios_raw, const std::array<bool, M>& mascara_raw,
                                     const std::array<bool, M>& market_open_raw,
                                     VentanaCausalCertificadaMultiactivoCxx<W, M>& out_ventana) noexcept {
        return paso(ts_ms, precios_raw.data(), mascara_raw.data(), market_open_raw.data(), out_ventana, nullptr, nullptr, nullptr);
    }

    // Telemetría y Contadores
    [[nodiscard]] constexpr size_t cursor_global() const noexcept { return cursor_global_; }
    [[nodiscard]] constexpr uint32_t ventanas_emitidas() const noexcept { return total_ventanas_emitidas_; }
    [[nodiscard]] constexpr uint32_t ventanas_invalidadas() const noexcept { return total_ventanas_invalidadas_; }
    [[nodiscard]] constexpr uint32_t gaps_detectados() const noexcept { return total_gaps_detectados_; }
    [[nodiscard]] constexpr uint32_t advertencias_rezago() const noexcept { return total_advertencias_rezago_; }
    [[nodiscard]] constexpr uint32_t asimetrias_union() const noexcept { return total_asimetrias_union_; }
    [[nodiscard]] constexpr uint32_t filas_descartadas() const noexcept { return total_filas_descartadas_; }
    [[nodiscard]] const std::array<double, M>& reloj_acum_por_activo() const noexcept { return nu_acum_por_activo_; }
    [[nodiscard]] const LazoFeedbackShannonCxx<W>& lazo_shannon(size_t m) const noexcept { return lazo_shannon_[m]; }

    void reset() noexcept {
        for (size_t m = 0; m < M; ++m) {
            ring_datos_[m].clear();
            ring_mask_[m].clear();
            ring_market_open_[m].clear();
            ring_dnu_[m].clear();
            ring_disp_[m].clear();
            ring_disp_decl_[m].clear();
            lazo_shannon_[m].reset();
            nu_acum_por_activo_[m] = 0.0;
            ultimo_val_valido_[m] = 0.0;
            tiene_valor_previo_[m] = false;
            disp_prev_[m] = 0;
            ult_cierre_lo_[m] = 0;
            ult_cierre_hi_[m] = 0;
        }
        ring_ts_.clear();
        ring_gaps_.clear();
        cursor_global_ = 0;
        siguiente_emision_ = W - 1;
        ts_prev_ = 0;
        ultimo_lookahead_score_ = 0.0;
        total_ventanas_emitidas_ = 0;
        total_ventanas_invalidadas_ = 0;
        total_gaps_detectados_ = 0;
        total_advertencias_rezago_ = 0;
        total_asimetrias_union_ = 0;
        total_filas_descartadas_ = 0;
    }
};

} // namespace ark::ventanas
