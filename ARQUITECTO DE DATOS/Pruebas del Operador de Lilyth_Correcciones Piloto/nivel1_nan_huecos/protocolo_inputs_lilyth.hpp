// protocolo_inputs_lilyth.hpp -- Protocolo unico de inputs para ObservadorLilyth apoyado en ark/dataset/
// Piloto Nivel 1 (ARQUITECTO DE DATOS / Auditoria Externa), C++20. Cero reservas dinamicas en ciclo caliente.
#pragma once

#include "ark/dataset/control_calidad.hpp"
#include "ark/dataset/tipos_dataset.hpp"
#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/lilyth/puente_historico_ito.hpp"
#include "ark/ventanas/flujo_barras.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace ark::lilyth::piloto {

// Ontologia operativa en el paso caliente de Lilyth, mapeada 1-a-1 desde ark::dataset::TipoHueco
enum class ClasePasoLilyth : std::uint8_t {
    VALIDO = 0,                 // Dato fisico presente y finito en sesion activa
    MERCADO_QUIETO = 1,         // ILIQUIDEZ_TRANSITORIA / fila imputada por ausencia de trades en sesion abierta
    FRONTERA_ESTRUCTURAL = 2,   // ESTRUCTURAL_CALENDARIO (fin de semana, feriado, pausa CME), INICIO_NIVEL o PURGA_INMADURA
    CORRUPCION_O_CAIDA = 3      // CORRUPCION_VALOR (NaN/Inf/<=0) o CAIDA_RED_SOCKET / CORTE_FEED / DESCONEXION
};

// Clasifica un paso de BarraCausal (ark/dataset/tipos_dataset.hpp) sin reservas dinamicas
[[nodiscard]] constexpr ClasePasoLilyth clasificar_barra_dataset(const ark::dataset::BarraCausal& b, double valor_extraido) noexcept {
    using TH = ark::dataset::TipoHueco;
    if (b.tipo_hueco == TH::CIERRE_PROGRAMADO || b.tipo_hueco == TH::FIN_DE_SEMANA || b.tipo_hueco == TH::FUERA_DE_CICLO) {
        return ClasePasoLilyth::FRONTERA_ESTRUCTURAL;
    }
    if ((b.flags_control & (ark::dataset::FLAG_FILA_VACIA | ark::dataset::FLAG_PAUSA_UMBRAL | ark::dataset::FLAG_CORRUPCION_CAIDA)) != 0u) {
        return ClasePasoLilyth::CORRUPCION_O_CAIDA;
    }
    if (b.tipo_hueco == TH::CORTE_FEED || b.tipo_hueco == TH::DESCONEXION) {
        return ClasePasoLilyth::CORRUPCION_O_CAIDA;
    }
    if (b.es_imputado) {
        if (!std::isfinite(valor_extraido)) {
            return ClasePasoLilyth::CORRUPCION_O_CAIDA;
        }
        // Fila insertada por rejilla en sesion activa con precio arrastrado: mercado quieto (ILIQUIDEZ_TRANSITORIA)
        return ClasePasoLilyth::MERCADO_QUIETO;
    }
    if (!b.hay_dato || !std::isfinite(valor_extraido)) {
        return ClasePasoLilyth::CORRUPCION_O_CAIDA;
    }
    return ClasePasoLilyth::VALIDO;
}

// Sintetiza los bits de flags_control / ontologia de hueco de una BarraCausal para propagarlos por EntradaPasoCxx -> BarraIntrinsecaCxx
[[nodiscard]] constexpr std::uint8_t etiquetar_flags_barra_causal(const ark::dataset::BarraCausal& b, double valor_extraido) noexcept {
    std::uint8_t f = b.flags_control;
    const ClasePasoLilyth cl = clasificar_barra_dataset(b, valor_extraido);
    if (cl == ClasePasoLilyth::MERCADO_QUIETO) {
        f |= ark::dataset::FLAG_MERCADO_QUIETO;
    } else if (cl == ClasePasoLilyth::CORRUPCION_O_CAIDA) {
        f |= ark::dataset::FLAG_CORRUPCION_CAIDA;
    }
    return f;
}

// Convierte una BarraCausal a EntradaPasoCxx preservando sin perdida la ontologia de huecos (Fase 1)
[[nodiscard]] inline ark::ventanas::EntradaPasoCxx a_entrada_paso_cxx(
    const ark::dataset::BarraCausal& b, double valor_extraido) noexcept {
    using TH = ark::dataset::TipoHueco;
    ark::ventanas::EntradaPasoCxx e{};
    e.ts_ms = b.cierre_real_utc_ms;
    e.valor = valor_extraido;
    const bool fuera_sesion = (b.tipo_hueco == TH::CIERRE_PROGRAMADO ||
                               b.tipo_hueco == TH::FIN_DE_SEMANA ||
                               b.tipo_hueco == TH::FUERA_DE_CICLO);
    e.con_sesion = !fuera_sesion;
    e.mascara = b.es_imputado || !b.hay_dato || !std::isfinite(valor_extraido);
    e.disponible_ms = b.estampa_disponible_ms;
    e.volumen = std::isfinite(b.volume) ? b.volume : 0.0;
    e.nocional = std::isfinite(b.valor_nominal) ? b.valor_nominal : 0.0;
    e.flags_hueco = etiquetar_flags_barra_causal(b, valor_extraido);
    return e;
}

// Clasifica una BarraIntrinsecaCxx ya cerrada leyendo sus etiquetas propagadas desde ark/ventanas (Fase 1)
[[nodiscard]] inline ClasePasoLilyth clasificar_barra_intrinseca(const ark::ventanas::BarraIntrinsecaCxx& b) noexcept {
    if (!std::isfinite(b.valor)) {
        return ClasePasoLilyth::CORRUPCION_O_CAIDA;
    }
    if ((b.flags_hueco & (ark::dataset::FLAG_FILA_VACIA | ark::dataset::FLAG_PAUSA_UMBRAL | ark::dataset::FLAG_CORRUPCION_CAIDA)) != 0u) {
        return ClasePasoLilyth::CORRUPCION_O_CAIDA;
    }
    if (b.imputada) {
        return b.solo_mercado_quieto ? ClasePasoLilyth::MERCADO_QUIETO : ClasePasoLilyth::CORRUPCION_O_CAIDA;
    }
    return ClasePasoLilyth::VALIDO;
}

[[nodiscard]] constexpr bool barra_tiene_puente_ito(const ark::ventanas::BarraIntrinsecaCxx& b) noexcept {
    return (b.flags_hueco & ark::dataset::FLAG_PUENTE_ITO) != 0u;
}

[[nodiscard]] constexpr bool barra_tiene_choque_en_hueco(const ark::ventanas::BarraIntrinsecaCxx& b) noexcept {
    return (b.flags_hueco & ark::dataset::FLAG_CHOQUE_EN_HUECO) != 0u;
}

[[nodiscard]] constexpr bool barra_tiene_pausa_umbral(const ark::ventanas::BarraIntrinsecaCxx& b) noexcept {
    return (b.flags_hueco & ark::dataset::FLAG_PAUSA_UMBRAL) != 0u;
}

// Clasifica desde el TipoHueco ontologico de control_calidad.hpp (namespace ark::dataset::ingesta)
[[nodiscard]] inline ClasePasoLilyth clasificar_hueco_ontologico(
    ark::dataset::ingesta::TipoHueco t, double valor, bool hay_hueco) noexcept {
    using THO = ark::dataset::ingesta::TipoHueco;
    if (!hay_hueco) {
        return std::isfinite(valor) ? ClasePasoLilyth::VALIDO : ClasePasoLilyth::CORRUPCION_O_CAIDA;
    }
    switch (t) {
        case THO::ILIQUIDEZ_TRANSITORIA:
            return std::isfinite(valor) ? ClasePasoLilyth::MERCADO_QUIETO : ClasePasoLilyth::CORRUPCION_O_CAIDA;
        case THO::ESTRUCTURAL_CALENDARIO:
        case THO::CENSURA_IZQUIERDA:
            return ClasePasoLilyth::FRONTERA_ESTRUCTURAL;
        case THO::CAIDA_RED_SOCKET:
        case THO::CORRUPCION_VALOR:
        default:
            return ClasePasoLilyth::CORRUPCION_O_CAIDA;
    }
}

// Estado de la maquina de pausa / espera acotada en el ciclo caliente (Fase 2)
enum class EstadoContinuidadLilyth : std::uint8_t {
    FLUJO_ACTIVO = 0,     // Operacion continua con datos validos o calma legitima
    ESPERA_ACOTADA = 1,   // Dentro del umbral (1 <= r <= racha_max_corrupcion): espera si vuelve antes de pausar
    EN_PAUSA_DURA = 2     // Excedio racha_max_corrupcion o frontera estructural: en pausa dura hasta que regrese dato
};

// Estado causal para preprocesar un flujo hacia ObservadorSaltos + ObservadorLilyth bajo politica controlada de huecos
struct AdaptadorInputCausalLilyth {
    ark::metrologia::Integracion integracion = ark::metrologia::Integracion::FLUJO;
    bool admitir_mercado_quieto = true;        // true: ILIQUIDEZ_TRANSITORIA no se marca como invalido
    std::size_t racha_max_corrupcion = 0;      // 0 = corta al 1er corrupto; > 0 = umbral temporal r_max
    double z_max_hueco = 3.5;                  // Llave de amplitud inicio-fin |Z_hueco| <= z_max_hueco
    double sigma_base = 1.0;                   // Escala base s0 / sd_innov calibrada en calma
    double tau_int_base = 1.0;                 // Tiempo integrado de autocorrelacion en calma
    double h_hurst_base = 0.5;                 // Exponente de escalado temporal H = 0.5 + d
    bool evitar_teletransporte_nivel = false;  // true (Fase 2): en NIVEL no arrastra ultimo_x_valido en CORRUPCION_O_CAIDA

    double ultimo_x_valido = ark::metrologia::NaN;
    bool tiene_ultimo = false;
    std::size_t racha_corrupta_actual = 0;
    EstadoContinuidadLilyth estado_actual = EstadoContinuidadLilyth::FLUJO_ACTIVO;

    // Contadores de auditoria del protocolo de entrada
    std::size_t n_validos = 0;
    std::size_t n_mercado_quieto = 0;
    std::size_t n_frontera_estructural = 0;
    std::size_t n_corrupcion_puenteada = 0;
    std::size_t n_corrupcion_cortada = 0;
    std::size_t n_pausas_umbral = 0;
    std::size_t n_choques_en_hueco = 0;
    std::size_t n_reconexiones_limpias = 0;

    struct PasoPreparado {
        double x_para_operador = ark::metrologia::NaN;
        bool valido_para_operador = false;
        ClasePasoLilyth clase = ClasePasoLilyth::VALIDO;
        bool puenteado = false;
        EstadoContinuidadLilyth estado = EstadoContinuidadLilyth::FLUJO_ACTIVO;
        std::uint8_t flags_emitidos = 0;
        bool reconexion_tras_hueco = false;
        bool choque_en_hueco_detectado = false;
        bool pausa_umbral_activada = false;
        bool cortar_segmento_al_reenganchar = false;
        std::size_t delta_t_hueco = 0;
        double delta_x_hueco = 0.0;
        double z_hueco = 0.0;
    };

    [[nodiscard]] double sigma_esperada_intervalo(std::size_t dt) const noexcept {
        const double s0 = (std::isfinite(sigma_base) && sigma_base > 1e-12) ? sigma_base : 1.0;
        const double ti = (std::isfinite(tau_int_base) && tau_int_base >= 1.0) ? tau_int_base : 1.0;
        if (integracion == ark::metrologia::Integracion::NIVEL) {
            const double h = std::clamp(h_hurst_base, 0.1, 0.95);
            return s0 * std::sqrt(ti) * std::pow(static_cast<double>( std::max<std::size_t>(1, dt) ), h);
        }
        return s0 * std::sqrt(ti);
    }

    PasoPreparado preparar(double x_crudo, ClasePasoLilyth clase, double valor_reposo_flujo = 0.0) noexcept {
        PasoPreparado out{};
        out.clase = clase;

        if (clase == ClasePasoLilyth::VALIDO && std::isfinite(x_crudo)) {
            if (racha_corrupta_actual > 0 && tiene_ultimo) {
                out.reconexion_tras_hueco = true;
                out.delta_t_hueco = racha_corrupta_actual + 1;
                out.delta_x_hueco = x_crudo - ultimo_x_valido;
                const double sig_dt = sigma_esperada_intervalo(out.delta_t_hueco);
                out.z_hueco = out.delta_x_hueco / std::max(1e-15, sig_dt);
                const bool excedio_r = (racha_corrupta_actual > racha_max_corrupcion);
                const bool excedio_z = (std::abs(out.z_hueco) > z_max_hueco);
                if (excedio_z) {
                    out.choque_en_hueco_detectado = true;
                    out.flags_emitidos |= ark::dataset::FLAG_CHOQUE_EN_HUECO;
                    ++n_choques_en_hueco;
                }
                if (excedio_r) {
                    out.pausa_umbral_activada = true;
                    out.flags_emitidos |= ark::dataset::FLAG_PAUSA_UMBRAL;
                }
                if (excedio_r || excedio_z || (evitar_teletransporte_nivel && integracion == ark::metrologia::Integracion::NIVEL)) {
                    out.cortar_segmento_al_reenganchar = true;
                } else {
                    ++n_reconexiones_limpias;
                }
            }
            ultimo_x_valido = x_crudo;
            tiene_ultimo = true;
            racha_corrupta_actual = 0;
            estado_actual = EstadoContinuidadLilyth::FLUJO_ACTIVO;
            out.estado = estado_actual;
            ++n_validos;
            out.x_para_operador = x_crudo;
            out.valido_para_operador = true;
            return out;
        }

        if (clase == ClasePasoLilyth::MERCADO_QUIETO) {
            ++n_mercado_quieto;
            racha_corrupta_actual = 0;
            estado_actual = EstadoContinuidadLilyth::FLUJO_ACTIVO;
            out.estado = estado_actual;
            out.flags_emitidos = ark::dataset::FLAG_MERCADO_QUIETO;
            if (!admitir_mercado_quieto) {
                out.x_para_operador = ark::metrologia::NaN;
                out.valido_para_operador = false;
                return out;
            }
            if (integracion == ark::metrologia::Integracion::NIVEL) {
                if (tiene_ultimo) {
                    out.x_para_operador = ultimo_x_valido; // dY = 0 exacto en mercado quieto
                    out.valido_para_operador = true;
                } else if (std::isfinite(x_crudo)) {
                    ultimo_x_valido = x_crudo;
                    tiene_ultimo = true;
                    out.x_para_operador = x_crudo;
                    out.valido_para_operador = true;
                }
            } else {
                out.x_para_operador = std::isfinite(x_crudo) ? x_crudo : valor_reposo_flujo;
                out.valido_para_operador = true;
            }
            return out;
        }

        if (clase == ClasePasoLilyth::FRONTERA_ESTRUCTURAL) {
            ++n_frontera_estructural;
            racha_corrupta_actual = 0;
            out.x_para_operador = ark::metrologia::NaN;
            out.valido_para_operador = false;
            return out;
        }

        // CORRUPCION_O_CAIDA
        ++racha_corrupta_actual;
        out.flags_emitidos = ark::dataset::FLAG_CORRUPCION_CAIDA;
        if (racha_corrupta_actual <= racha_max_corrupcion && tiene_ultimo) {
            estado_actual = EstadoContinuidadLilyth::ESPERA_ACOTADA;
            out.estado = estado_actual;
            // En NIVEL con evitar_teletransporte_nivel == true, NO se arrastra ultimo_x_valido durante CORRUPCION_O_CAIDA
            // para evitar fabricar dY = 0 durante el corte y concentrar todo el salto (0, 0, +Delta) al reconectar
            if (evitar_teletransporte_nivel && integracion == ark::metrologia::Integracion::NIVEL) {
                ++n_corrupcion_cortada;
                out.x_para_operador = ark::metrologia::NaN;
                out.valido_para_operador = false;
                return out;
            }
            ++n_corrupcion_puenteada;
            out.puenteado = true;
            out.x_para_operador = (integracion == ark::metrologia::Integracion::NIVEL) ? ultimo_x_valido : valor_reposo_flujo;
            out.valido_para_operador = true;
            return out;
        }
        if (estado_actual != EstadoContinuidadLilyth::EN_PAUSA_DURA) {
            ++n_pausas_umbral;
        }
        estado_actual = EstadoContinuidadLilyth::EN_PAUSA_DURA;
        out.estado = estado_actual;
        out.pausa_umbral_activada = true;
        out.flags_emitidos |= ark::dataset::FLAG_PAUSA_UMBRAL;
        ++n_corrupcion_cortada;
        out.x_para_operador = ark::metrologia::NaN;
        out.valido_para_operador = false;
        return out;
    }
};

} // namespace ark::lilyth::piloto
