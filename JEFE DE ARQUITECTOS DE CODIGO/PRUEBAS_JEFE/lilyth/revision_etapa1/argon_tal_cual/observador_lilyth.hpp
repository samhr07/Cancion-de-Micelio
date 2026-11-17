// observador_lilyth.hpp -- Operador de Lilyth (restitucion tras un evento), Etapa 1: tipos, trayectoria eta(tau),
// ranuras con ancla congelada (cruda, purgada y contraste de marea) y riesgos competitivos de Aalen-Johansen, C++20.
// Auditoria Externa (Argon) y Jefe de Arquitectos de Codigo (Code, Opus 5.5), 2026-11-05.
// Especificacion y prerregistro: README_OPERADOR_LILYTH.md (Secciones 1-9, Enmienda 1 E1.1-E1.16).
//
// Reglas de la base auditada cumplidas en esta cabecera:
//   1. No modifica ObservadorSaltos ni canal_nivel: envuelve o consume su salida permaneciendo identico bit a bit.
//   2. Ciclo caliente (paso): cero reservas de memoria dinamica (todo en std::array fijos en la instancia).
//   3. Causal por construccion: el ancla Y* se congela en t_cambio - guarda - 1 sobre ventana W_ancla = N_REG * 2^j0.
//   4. Soportadas las tres representaciones (PURGADA, CRUDA, AMBAS_CONTRASTE_MAREA), reloj doble y Aalen-Johansen.
#pragma once

#include "tipos_metrologia.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace ark::lilyth {

using ark::metrologia::ALFA_PROYECTO;
using ark::metrologia::EstadoSaltos;
using ark::metrologia::FichaEvento;
using ark::metrologia::Integracion;
using ark::metrologia::NaN;

// Codigos de nulos compatibles con ark/metrologia/sustitutos.hpp (sin acoplar cabecera si se compila aislado).
inline constexpr std::uint32_t LILYTH_NULO_NINGUNO     = 0u;
inline constexpr std::uint32_t LILYTH_NULO_PERMUTACION = 1u << 0;
inline constexpr std::uint32_t LILYTH_NULO_BLOQUES     = 1u << 1;
inline constexpr std::uint32_t LILYTH_NULO_SIGNOS      = 1u << 2;

// E1.14: Representacion sobre la que mide Lilyth (entrada del pipeline, decision del Comandante).
enum class ModoRepresentacionLilyth : std::uint8_t {
    PURGADA = 0,
    CRUDA = 1,
    AMBAS_CONTRASTE_MAREA = 2
};

// E1.13 y 8.2(b): Modelo lineal para el nulo de Onsager (entrada del pipeline).
enum class ModoNuloLineal : std::uint8_t {
    AR_LEVINSON_DURBIN = 0,
    ARFIMA_WHITTLE = 1,
    FASE_ALEATORIZADA_FFT = 2
};

// E1.8 y E1.13: Politica del detector de armonicos (entrada del pipeline).
enum class ModoArmonico : std::uint8_t {
    DESACTIVADO = 0,
    PRONY_CURVA_MEDIA = 1,
    PRONY_FASE_REESCALADA = 2,
    SOBREPASO_POR_EVENTO = 3,
    TODOS = 4
};

// E1.10 y E1.14: Estado y desenlace de cada rama de una ranura en el analisis de supervivencia competitiva.
enum class DesenlaceRanura : std::uint8_t {
    INACTIVA = 0,
    ABIERTA = 1,
    RESTITUIDA = 2,                  // cruzo eta(tau) <= 0.5 (consolidado sin calor) antes de replica y de T_max
    CORTADA_POR_REPLICA = 3,         // una nueva alarma del canal de nivel corto la trayectoria antes de restituir
    CENSURADA_HORIZONTE = 4,         // llego a t_max_efectivo sin restituir ni sufrir replica
    CENSURADA_DATO_INVALIDO = 5,     // hueco / barra invalida o incremento purgado invalido en NIVEL (E1.14)
    NO_MEDIDA_PURGA_INMADURA = 6,    // la purga diurna aun no cumplio dias_minimos (20) en t_cambio^- (E1.14)
    NO_MEDIDA_ANCLA_INSUFICIENTE = 7 // historia previa insuficiente para formar el ancla con su guarda (E1.3)
};

// Configuracion de entrada del Operador de Lilyth (E1.13, E1.14, E1.15).
struct ConfiguracionLilyth {
    Integracion integracion = Integracion::FLUJO;
    ModoRepresentacionLilyth representacion = ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA;
    std::uint32_t nulo_sin_memoria = LILYTH_NULO_PERMUTACION;   // E1.2, E1.13: solo PERMUTACION o SIGNOS
    ModoNuloLineal modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
    ModoArmonico modo_armonico = ModoArmonico::PRONY_CURVA_MEDIA;
    double alfa = ALFA_PROYECTO;
    std::size_t t_max_efectivo = 128;   // E1.9: cuantil 1-alfa del nulo lineal o W_ancla si hay raiz unitaria
    bool cortar_abiertas_con_replica = true;   // D8, E1.10: una alarma nueva corta las ranuras abiertas no restituidas
};

// Validacion de configuracion conforme a E1.2 y E1.13 (rechaza NULO_NINGUNO y NULO_BLOQUES en nulo_sin_memoria).
[[nodiscard]] constexpr bool configuracion_lilyth_valida(const ConfiguracionLilyth& cfg, std::size_t t_max_cap) noexcept {
    const bool nulo_ok = (cfg.nulo_sin_memoria == LILYTH_NULO_PERMUTACION) ||
                         (cfg.nulo_sin_memoria == LILYTH_NULO_SIGNOS);
    const bool alfa_ok = (cfg.alfa > 0.0) && (cfg.alfa < 1.0);
    const bool tmax_ok = (cfg.t_max_efectivo >= 2) && (cfg.t_max_efectivo <= t_max_cap);
    return nulo_ok && alfa_ok && tmax_ok;
}

// Trayectoria y metricas de una rama (PURGADA o CRUDA) dentro de una ranura de Lilyth.
template <std::size_t T_MAX>
struct RamaTrayectoriaLilyth {
    static_assert(T_MAX >= 2, "T_MAX >= 2");
    DesenlaceRanura desenlace = DesenlaceRanura::INACTIVA;
    double y_ancla = NaN;                     // Y* congelado en t_cambio - guarda - 1 (E1.3)
    double y_disparo = NaN;                   // Y(t_a) en el instante del disparo
    double dy = NaN;                          // dY = Y(t_a) - Y* (E1.4)
    double nu_punto_ancla = NaN;              // velocidad basal del reloj de actividad en W_ancla (E1.11)
    double nu_en_disparo = NaN;               // reloj de actividad acumulado en t_a
    std::size_t muestras_ancla = 0;           // muestras validas usadas en la mediana de Y*
    std::size_t pasos_observados = 0;         // ultimo tau registrado (1 .. T_MAX)
    std::size_t tau_fin = 0;                  // paso tau en que se cerro la rama (restitucion, replica o censura)

    // E1.4 y E1.6: tiempos de restitucion en reloj de pared (pasos desde t_a) y en reloj de actividad normalizado.
    double tau_rest_primer_toque = NaN;       // primer tau >= 1 con eta(tau) <= 0.5 (incluso con calor activo)
    double tau_rest_consolidado = NaN;        // primer tau >= 1 con eta(tau) <= 0.5 y sin calor activo
    double delta_nu_primer_toque = NaN;       // nu(t_a + tau_1) - nu(t_a)
    double delta_nu_consolidado = NaN;        // nu(t_a + tau_cons) - nu(t_a)
    double tau_act_primer_toque = NaN;        // delta_nu_primer_toque / nu_punto_ancla (barras equivalentes de calma)
    double tau_act_consolidado = NaN;         // delta_nu_consolidado / nu_punto_ancla

    // E1.8 y D10: seguimiento de sobrepaso individual (tau >= 2).
    double min_eta_desde_2 = NaN;             // min_{tau >= 2} eta(tau)
    std::size_t tau_min_eta = 0;              // tau donde se alcanzo min_eta_desde_2
    bool armonico_sospechoso = false;         // D10 / seccion 7.1: etiqueta (nunca veto) si min_eta_desde_2 < 0

    std::array<double, T_MAX> eta{};          // eta(tau) para tau = 1 .. pasos_observados (indice tau - 1)

    void limpiar() noexcept {
        desenlace = DesenlaceRanura::INACTIVA;
        y_ancla = y_disparo = dy = nu_punto_ancla = nu_en_disparo = NaN;
        muestras_ancla = pasos_observados = tau_fin = 0;
        tau_rest_primer_toque = tau_rest_consolidado = NaN;
        delta_nu_primer_toque = delta_nu_consolidado = NaN;
        tau_act_primer_toque = tau_act_consolidado = NaN;
        min_eta_desde_2 = NaN;
        tau_min_eta = 0;
        armonico_sospechoso = false;
        eta.fill(NaN);
    }
};

// Ranura completa de un evento seguida por Lilyth (E1.3, E1.4, E1.12, E1.14).
template <std::size_t T_MAX, std::size_t N_CAPA1 = 4>
struct RanuraLilyth {
    bool ocupada = false;
    std::uint64_t id_evento = 0;
    std::int64_t t_apertura = -1;             // paso nativo t_a en que se abrio la ranura (Flanco 7.2)
    std::int64_t t_cambio = -1;               // punto de cambio estimado por el GLR
    std::int64_t t_fin_ancla = -1;            // ultimo paso incluido en W_ancla (t_cambio - guarda - 1)
    std::size_t w_ancla = 0;                  // N_REG * 2^j0 barras nativas (E1.3)
    std::size_t guarda_ancla = 0;             // ceil(tau_int,j0) * 2^j0 barras nativas (E1.3)
    int signo = 0;                            // +1 o -1 del choque de nivel
    int escala_origen = 0;                    // j0
    double amplitud_glr = NaN;                // guardada aparte (E1.4)
    double se_amplitud_glr = NaN;
    double g_disparo = NaN;                   // estadistico GLR en el disparo (para estratificar excedencias g >= g_k)
    bool ancla_saturada_en_l = false;         // Flanco 7.1: true si alarma.l == L
    int signo_replica_corte = 0;              // si fue cortada por replica, signo de la replica que la corto

    // E1.12 (F3, P6): fotografia congelada por valor de las medianas de la Capa 1 (SensorEkaLilyth) en t_cambio^-.
    std::array<double, N_CAPA1> medianas_capa1_congeladas{};

    FichaEvento ficha_evento{};               // copia integra de la ficha del operador delta (D4)
    RamaTrayectoriaLilyth<T_MAX> purgada{};
    RamaTrayectoriaLilyth<T_MAX> cruda{};

    // E1.14: contraste directo de marea diurna en el evento cuando ambas ramas han restituido.
    double delta_tau_marea_primer_toque = NaN;   // tau_cruda - tau_purgada (primer toque)
    double delta_tau_marea_consolidado = NaN;    // tau_cruda - tau_purgada (consolidado)

    void limpiar() noexcept {
        ocupada = false;
        id_evento = 0;
        t_apertura = t_cambio = t_fin_ancla = -1;
        w_ancla = guarda_ancla = 0;
        signo = escala_origen = signo_replica_corte = 0;
        amplitud_glr = se_amplitud_glr = g_disparo = NaN;
        ancla_saturada_en_l = false;
        medianas_capa1_congeladas.fill(NaN);
        ficha_evento = FichaEvento{};
        purgada.limpiar();
        cruda.limpiar();
        delta_tau_marea_primer_toque = delta_tau_marea_consolidado = NaN;
    }

    [[nodiscard]] bool alguna_rama_abierta() const noexcept {
        return (purgada.desenlace == DesenlaceRanura::ABIERTA) ||
               (cruda.desenlace == DesenlaceRanura::ABIERTA);
    }
};

// =================================================================================================
// E1.10: Estimador no parametrico de Aalen-Johansen (riesgos competitivos en tiempo discreto)
// =================================================================================================
// Desenlaces competitivos en tau = 1 .. T_MAX:
//   - Evento 1 (restitucion consolidada en tau): d_rest[tau]
//   - Evento 2 (corte por replica antes de restituir en tau): d_replica[tau]
//   - Censura (horizonte T_max o dato invalido en tau_fin): sale del conjunto de riesgo tras tau_fin.
template <std::size_t T_MAX>
struct CurvaAalenJohansen {
    std::size_t n_eventos_registrados = 0;
    std::size_t n_restituidos = 0;
    std::size_t n_cortados_replica = 0;
    std::size_t n_censurados = 0;
    std::size_t t_max_efectivo = T_MAX;

    // Contadores crudos por paso tau = 1 .. T_MAX (almacenados en indice tau - 1):
    std::array<std::size_t, T_MAX> d_rest{};          // restituyen exactamente en tau
    std::array<std::size_t, T_MAX> d_replica{};       // cortados por replica exactamente en tau
    std::array<std::size_t, T_MAX> censurados_en{};   // censurados tras completar el paso tau (salen en tau + 1)
    std::size_t censurados_en_cero = 0;               // censurados en tau = 0 antes de dar un paso valido

    // Curvas estimadas para tau = 1 .. t_max_efectivo (indice tau - 1):
    std::array<std::size_t, T_MAX> en_riesgo{};       // n(tau)
    std::array<double, T_MAX> supervivientes{};       // S(tau)
    std::array<double, T_MAX> cif_rest{};             // CIF_rest(tau)
    std::array<double, T_MAX> cif_replica{};          // CIF_replica(tau)
    double rmst = NaN;                                // RMST(T_max) = sum_{tau=1}^{T_max} (1 - CIF_rest(tau))
    double tau_mediano_cif_rest = NaN;                // primer tau con CIF_rest(tau) >= 0.5 (NaN si no alcanza 0.5)

    void limpiar(std::size_t tmax_ef = T_MAX) noexcept {
        n_eventos_registrados = n_restituidos = n_cortados_replica = n_censurados = 0;
        t_max_efectivo = std::min(tmax_ef, T_MAX);
        d_rest.fill(0);
        d_replica.fill(0);
        censurados_en.fill(0);
        censurados_en_cero = 0;
        en_riesgo.fill(0);
        supervivientes.fill(NaN);
        cif_rest.fill(NaN);
        cif_replica.fill(NaN);
        rmst = tau_mediano_cif_rest = NaN;
    }

    // Registra el desenlace de una rama cerrada y recalcula la tabla en O(T_MAX) sin reservas.
    void registrar_rama(DesenlaceRanura desenlace, std::size_t tau_fin) noexcept {
        if (desenlace == DesenlaceRanura::INACTIVA ||
            desenlace == DesenlaceRanura::ABIERTA ||
            desenlace == DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA ||
            desenlace == DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE) {
            return;
        }
        const std::size_t tf = std::min(tau_fin, t_max_efectivo);
        ++n_eventos_registrados;
        if (desenlace == DesenlaceRanura::RESTITUIDA && tf >= 1) {
            ++d_rest[tf - 1];
            ++n_restituidos;
        } else if (desenlace == DesenlaceRanura::CORTADA_POR_REPLICA && tf >= 1) {
            ++d_replica[tf - 1];
            ++n_cortados_replica;
        } else {
            ++n_censurados;
            if (tf == 0) ++censurados_en_cero;
            else ++censurados_en[tf - 1];
        }
        recalcular();
    }

    void recalcular() noexcept {
        if (n_eventos_registrados == 0 || t_max_efectivo == 0) {
            rmst = tau_mediano_cif_rest = NaN;
            return;
        }
        std::size_t n_act = n_eventos_registrados >= censurados_en_cero
                                ? (n_eventos_registrados - censurados_en_cero)
                                : 0;
        double s_prev = 1.0, c_rest = 0.0, c_rep = 0.0, suma_rmst = 0.0;
        tau_mediano_cif_rest = NaN;

        for (std::size_t idx = 0; idx < t_max_efectivo; ++idx) {
            en_riesgo[idx] = n_act;
            const std::size_t dr = d_rest[idx];
            const std::size_t dp = d_replica[idx];
            if (n_act > 0) {
                const double denom = static_cast<double>(n_act);
                c_rest += s_prev * (static_cast<double>(dr) / denom);
                c_rep  += s_prev * (static_cast<double>(dp) / denom);
                const std::size_t salidas_ev = std::min(n_act, dr + dp);
                s_prev *= static_cast<double>(n_act - salidas_ev) / denom;
            }
            supervivientes[idx] = s_prev;
            cif_rest[idx] = c_rest;
            cif_replica[idx] = c_rep;
            suma_rmst += (1.0 - c_rest);
            if (!std::isfinite(tau_mediano_cif_rest) && c_rest >= 0.5) {
                tau_mediano_cif_rest = static_cast<double>(idx + 1);
            }
            const std::size_t bajas = dr + dp + censurados_en[idx];
            n_act = (n_act > bajas) ? (n_act - bajas) : 0;
        }
        rmst = suma_rmst;
    }
};

// Desglose por signo del choque (Flanco 7.4 y E1.10).
template <std::size_t T_MAX>
struct ResumenAalenJohansenSigno {
    CurvaAalenJohansen<T_MAX> total{};
    CurvaAalenJohansen<T_MAX> positivos{};   // signo == +1
    CurvaAalenJohansen<T_MAX> negativos{};   // signo == -1

    void limpiar(std::size_t tmax_ef = T_MAX) noexcept {
        total.limpiar(tmax_ef);
        positivos.limpiar(tmax_ef);
        negativos.limpiar(tmax_ef);
    }

    void registrar(int signo, DesenlaceRanura desenlace, std::size_t tau_fin) noexcept {
        total.registrar_rama(desenlace, tau_fin);
        if (signo > 0) positivos.registrar_rama(desenlace, tau_fin);
        else if (signo < 0) negativos.registrar_rama(desenlace, tau_fin);
    }

    [[nodiscard]] double asimetria_rmst() const noexcept {
        return negativos.rmst - positivos.rmst;   // Flanco 7.4: RMST^- - RMST^+
    }
};

// =================================================================================================
// Orquestador de Etapa 1: ObservadorLilyth<Obs, RANURAS, T_MAX, HIST_MAX, N_CAPA1>
// =================================================================================================
// Envuelve o acompana a `Obs` (p. ej. ObservadorSaltos) sin tocar una sola linea de `Obs`.
// El anillo causal de historia guarda los ultimos HIST_MAX pasos de:
//   - y_cruda(t): nivel crudo (o flujo crudo),
//   - y_purgada(t): nivel purgado reintegrado (o flujo purgado) y su bandera `purga_valida`,
//   - d_nu(t): incremento del reloj de actividad (NU_TRADES, NU_VOLUMEN o NU_NOCIONAL),
//   - medianas_capa1(t): estado del SensorEkaLilyth para congelar por copia en t_cambio^- (E1.12).
template <
    class Obs,
    std::size_t RANURAS = 8,
    std::size_t T_MAX = 512,
    std::size_t HIST_MAX = 4096,
    std::size_t N_CAPA1 = 4>
class ObservadorLilyth {
    static_assert(RANURAS >= 1, "RANURAS >= 1");
    static_assert(T_MAX >= 2, "T_MAX >= 2");
    static_assert(HIST_MAX >= 64, "HIST_MAX >= 64");

public:
    struct SalidaPaso {
        bool emitida = false;
        std::int64_t t_nativo = -1;
        bool abrio_ranura = false;
        std::size_t indice_ranura_abierta = 0;
        std::size_t ranuras_cerradas_en_paso = 0;
        std::array< RanuraLilyth<T_MAX, N_CAPA1>, RANURAS > cerradas{};
    };

private:
    Obs obs_{};
    ConfiguracionLilyth cfg_{};
    bool cfg_valida_ = false;

    // Anillo causal de historia pre-evento (cero reservas en ciclo caliente).
    std::array<double, HIST_MAX> hist_cruda_{};
    std::array<double, HIST_MAX> hist_purgada_{};
    std::array<double, HIST_MAX> hist_dnu_{};
    std::array<std::uint8_t, HIST_MAX> hist_val_cruda_{};
    std::array<std::uint8_t, HIST_MAX> hist_val_purgada_{};
    std::array<std::array<double, N_CAPA1>, HIST_MAX> hist_capa1_{};
    std::array<double, HIST_MAX> scratch_mediana_{};   // buffer fijo para nth_element (cero heap)

    std::size_t cab_hist_ = 0;
    std::size_t n_hist_ = 0;
    std::int64_t t_actual_ = -1;

    // Estado de reintegracion causal para Integracion::NIVEL (E1.14).
    bool nivel_crudo_iniciado_ = false;
    bool nivel_purgado_iniciado_ = false;
    bool nivel_purgado_roto_ = false;
    double acum_nivel_crudo_ = 0.0;
    double acum_nivel_purgado_ = 0.0;
    double y_cruda_anterior_ = NaN;
    double nu_acumulado_ = 0.0;

    // Estado termico modulado por canal_temperatura (D3, E1.6, E1.7).
    bool calor_activo_ = false;

    // Banco de ranuras activas y contadores de telemetria (D8).
    std::array<RanuraLilyth<T_MAX, N_CAPA1>, RANURAS> ranuras_{};
    std::uint64_t seq_eventos_ = 0;
    std::size_t desbordes_ranura_ = 0;
    std::size_t eventos_nivel_recibidos_ = 0;
    std::size_t eventos_temp_recibidos_ = 0;

    // Tablas acumuladas de Aalen-Johansen para la rama PURGADA y la rama CRUDA.
    ResumenAalenJohansenSigno<T_MAX> aj_purgada_{};
    ResumenAalenJohansenSigno<T_MAX> aj_cruda_{};

    // Curva media acumulada eta_barra(tau) (E1.5 P3b y E1.8).
    std::array<double, T_MAX> suma_eta_purgada_{};
    std::array<std::size_t, T_MAX> cuenta_eta_purgada_{};
    std::array<double, T_MAX> suma_eta_cruda_{};
    std::array<std::size_t, T_MAX> cuenta_eta_cruda_{};
    std::size_t excluidos_prony_reescalada_purgada_ = 0;   // E1.13: eventos con tau_rest == NaN excluidos
    std::size_t excluidos_prony_reescalada_cruda_ = 0;

public:
    ObservadorLilyth() noexcept { reiniciar_todo(); }

    [[nodiscard]] Obs& operador_delta() noexcept { return obs_; }
    [[nodiscard]] const Obs& operador_delta() const noexcept { return obs_; }
    [[nodiscard]] const ConfiguracionLilyth& configuracion() const noexcept { return cfg_; }
    [[nodiscard]] bool calibracion_valida() const noexcept { return cfg_valida_; }
    [[nodiscard]] bool calor_activo() const noexcept { return calor_activo_; }
    [[nodiscard]] std::size_t desbordes_ranura() const noexcept { return desbordes_ranura_; }
    [[nodiscard]] std::size_t eventos_nivel_recibidos() const noexcept { return eventos_nivel_recibidos_; }
    [[nodiscard]] std::size_t eventos_temperatura_recibidos() const noexcept { return eventos_temp_recibidos_; }
    [[nodiscard]] const std::array<RanuraLilyth<T_MAX, N_CAPA1>, RANURAS>& ranuras() const noexcept { return ranuras_; }
    [[nodiscard]] const ResumenAalenJohansenSigno<T_MAX>& aalen_johansen_purgada() const noexcept { return aj_purgada_; }
    [[nodiscard]] const ResumenAalenJohansenSigno<T_MAX>& aalen_johansen_cruda() const noexcept { return aj_cruda_; }
    [[nodiscard]] std::size_t excluidos_prony_reescalada_purgada() const noexcept { return excluidos_prony_reescalada_purgada_; }
    [[nodiscard]] std::size_t excluidos_prony_reescalada_cruda() const noexcept { return excluidos_prony_reescalada_cruda_; }

    [[nodiscard]] double eta_barra_purgada(std::size_t tau) const noexcept {
        if (tau < 1 || tau > T_MAX || cuenta_eta_purgada_[tau - 1] == 0) return NaN;
        return suma_eta_purgada_[tau - 1] / static_cast<double>(cuenta_eta_purgada_[tau - 1]);
    }

    [[nodiscard]] double eta_barra_cruda(std::size_t tau) const noexcept {
        if (tau < 1 || tau > T_MAX || cuenta_eta_cruda_[tau - 1] == 0) return NaN;
        return suma_eta_cruda_[tau - 1] / static_cast<double>(cuenta_eta_cruda_[tau - 1]);
    }

    bool configurar(const ConfiguracionLilyth& cfg) noexcept {
        cfg_ = cfg;
        cfg_valida_ = configuracion_lilyth_valida(cfg_, T_MAX);
        reiniciar_flujo();
        return cfg_valida_;
    }

    void fijar_calor_activo(bool activo) noexcept { calor_activo_ = activo; }

    void reiniciar_flujo() noexcept {
        cab_hist_ = 0;
        n_hist_ = 0;
        t_actual_ = -1;
        nivel_crudo_iniciado_ = false;
        nivel_purgado_iniciado_ = false;
        nivel_purgado_roto_ = false;
        acum_nivel_crudo_ = 0.0;
        acum_nivel_purgado_ = 0.0;
        y_cruda_anterior_ = NaN;
        nu_acumulado_ = 0.0;
        calor_activo_ = false;
        for (auto& r : ranuras_) r.limpiar();
        const std::size_t tmax_ef = cfg_valida_ ? cfg_.t_max_efectivo : T_MAX;
        aj_purgada_.limpiar(tmax_ef);
        aj_cruda_.limpiar(tmax_ef);
        suma_eta_purgada_.fill(0.0);
        cuenta_eta_purgada_.fill(0);
        suma_eta_cruda_.fill(0.0);
        cuenta_eta_cruda_.fill(0);
        excluidos_prony_reescalada_purgada_ = 0;
        excluidos_prony_reescalada_cruda_ = 0;
    }

    void reiniciar_todo() noexcept {
        cfg_ = ConfiguracionLilyth{};
        cfg_valida_ = configuracion_lilyth_valida(cfg_, T_MAX);
        seq_eventos_ = 0;
        desbordes_ranura_ = 0;
        eventos_nivel_recibidos_ = 0;
        eventos_temp_recibidos_ = 0;
        reiniciar_flujo();
    }

    // Ciclo caliente: ejecuta el paso en el operador delta `obs_` y alimenta Lilyth sin tocar `obs_`.
    // Devuelve por referencia la salida intacta de `obs_` (identidad bit a bit garantizada).
    template <class SalidaObs>
    SalidaPaso paso_con_salida_delta(
        double y_cruda_entrada,
        bool valido_crudo,
        double factor_purga_diurna,   // sqrt(v_media / v_fase) si la purga ya maduro (>= 20 dias); NaN si inmadura
        bool purga_madura,
        double d_nu_actividad,
        const SalidaObs& salida_delta,
        SalidaPaso& out_lilyth,
        const std::array<double, N_CAPA1>& medianas_capa1 = {}) noexcept
    {
        out_lilyth = SalidaPaso{};
        if (!cfg_valida_) return out_lilyth;

        ++t_actual_;
        out_lilyth.emitida = true;
        out_lilyth.t_nativo = t_actual_;

        const double dnu_limpio = (valido_crudo && std::isfinite(d_nu_actividad) && d_nu_actividad >= 0.0)
                                      ? d_nu_actividad
                                      : 0.0;
        nu_acumulado_ += dnu_limpio;

        // 1. Construccion causal de los valores Y_cruda(t) e Y_purgada(t) segun FLUJO o NIVEL (E1.14).
        double y_cruda_nivel = NaN;
        double y_purgada_nivel = NaN;
        bool ok_cruda = false;
        bool ok_purgada = false;

        if (cfg_.integracion == Integracion::FLUJO) {
            if (valido_crudo && std::isfinite(y_cruda_entrada)) {
                y_cruda_nivel = y_cruda_entrada;
                ok_cruda = true;
                if (purga_madura && std::isfinite(factor_purga_diurna) && factor_purga_diurna > 0.0) {
                    // Correccion de Code a E1.14: PurgaDiurna en escalera_diadica.hpp hace y * sqrt(v_media / v_fase).
                    y_purgada_nivel = y_cruda_entrada * factor_purga_diurna;
                    ok_purgada = std::isfinite(y_purgada_nivel);
                }
            }
        } else {
            // Integracion::NIVEL: la entrada es el nivel crudo Y(t); su incremento nativo es dY(t) = Y(t) - Y(t-1).
            if (valido_crudo && std::isfinite(y_cruda_entrada)) {
                y_cruda_nivel = y_cruda_entrada;
                ok_cruda = true;
                if (!nivel_crudo_iniciado_) {
                    nivel_crudo_iniciado_ = true;
                    acum_nivel_crudo_ = y_cruda_entrada;
                    y_cruda_anterior_ = y_cruda_entrada;
                    if (purga_madura && std::isfinite(factor_purga_diurna) && factor_purga_diurna > 0.0) {
                        nivel_purgado_iniciado_ = true;
                        nivel_purgado_roto_ = false;
                        acum_nivel_purgado_ = y_cruda_entrada;
                        y_purgada_nivel = acum_nivel_purgado_;
                        ok_purgada = true;
                    }
                } else {
                    const double dy_nat = y_cruda_entrada - y_cruda_anterior_;
                    y_cruda_anterior_ = y_cruda_entrada;
                    acum_nivel_crudo_ = y_cruda_entrada;
                    if (purga_madura && std::isfinite(factor_purga_diurna) && factor_purga_diurna > 0.0) {
                        if (!nivel_purgado_iniciado_) {
                            nivel_purgado_iniciado_ = true;
                            nivel_purgado_roto_ = false;
                            acum_nivel_purgado_ = y_cruda_entrada;
                        } else if (!nivel_purgado_roto_) {
                            acum_nivel_purgado_ += dy_nat * factor_purga_diurna;
                        }
                        if (!nivel_purgado_roto_ && std::isfinite(acum_nivel_purgado_)) {
                            y_purgada_nivel = acum_nivel_purgado_;
                            ok_purgada = true;
                        } else {
                            nivel_purgado_roto_ = true;
                            ok_purgada = false;
                        }
                    } else if (nivel_purgado_iniciado_) {
                        // Correccion de Code a E1.14: en NIVEL, un incremento purgado invalido rompe la reintegracion
                        // y censura las ranuras abiertas por fin de dato (nunca se rellena).
                        nivel_purgado_roto_ = true;
                        ok_purgada = false;
                    }
                }
            } else {
                y_cruda_anterior_ = NaN;
                nivel_crudo_iniciado_ = false;
                if (nivel_purgado_iniciado_) nivel_purgado_roto_ = true;
            }
        }

        // 2. Guardar en el anillo causal de historia ANTES de procesar el paso actual en las ranuras.
        guardar_en_historia(y_cruda_nivel, ok_cruda, y_purgada_nivel, ok_purgada, dnu_limpio, medianas_capa1);

        // 3. Extraer eventos del operador delta en este paso (nivel vs temperatura, E1.15).
        FichaEvento ev_nivel{};
        bool hay_ev_nivel = false;
        extraer_eventos_de_salida_delta(salida_delta, ev_nivel, hay_ev_nivel);

        // 4. Avanzar un paso (tau = t_actual_ - r.t_apertura) en todas las ranuras que ya estaban abiertas.
        for (std::size_t i = 0; i < RANURAS; ++i) {
            auto& r = ranuras_[i];
            if (!r.ocupada) continue;
            if (t_actual_ <= r.t_apertura) continue;
            const std::size_t tau = static_cast<std::size_t>(t_actual_ - r.t_apertura);

            avanzar_rama(r.purgada, tau, y_purgada_nivel, ok_purgada, true);
            avanzar_rama(r.cruda, tau, y_cruda_nivel, ok_cruda, false);

            // D8 y E1.10: si llega un nuevo evento de nivel (replica), corta las ramas abiertas que aun no han
            // consolidado su restitucion.
            if (hay_ev_nivel && cfg_.cortar_abiertas_con_replica) {
                r.signo_replica_corte = ev_nivel.alarma.signo;
                if (r.purgada.desenlace == DesenlaceRanura::ABIERTA &&
                    !std::isfinite(r.purgada.tau_rest_consolidado)) {
                    r.purgada.desenlace = DesenlaceRanura::CORTADA_POR_REPLICA;
                    r.purgada.tau_fin = tau;
                }
                if (r.cruda.desenlace == DesenlaceRanura::ABIERTA &&
                    !std::isfinite(r.cruda.tau_rest_consolidado)) {
                    r.cruda.desenlace = DesenlaceRanura::CORTADA_POR_REPLICA;
                    r.cruda.tau_fin = tau;
                }
            }

            // Si ambas ramas terminaron su seguimiento, cerramos la ranura y registramos en Aalen-Johansen.
            if (!r.alguna_rama_abierta()) {
                cerrar_y_archivar_ranura(r, out_lilyth);
            }
        }

        // 5. Si hubo un disparo de nivel en este paso, abrir una nueva ranura con t_apertura = t_actual_ (Flanco 7.2).
        if (hay_ev_nivel) {
            abrir_ranura_desde_evento(ev_nivel, y_cruda_nivel, ok_cruda, y_purgada_nivel, ok_purgada, purga_madura, out_lilyth);
        }

        return out_lilyth;
    }

    // Sobrecarga directa cuando `Obs` expone `paso(y, valido)` devolviendo una salida con `ficha_evento` o `eventos`.
    SalidaPaso paso(
        double y_cruda_entrada,
        bool valido_crudo,
        double factor_purga_diurna = 1.0,
        bool purga_madura = true,
        double d_nu_actividad = 1.0,
        const std::array<double, N_CAPA1>& medianas_capa1 = {}) noexcept
    {
        const auto salida_delta = obs_.paso(y_cruda_entrada, valido_crudo);
        SalidaPaso out{};
        paso_con_salida_delta(
            y_cruda_entrada, valido_crudo, factor_purga_diurna, purga_madura, d_nu_actividad, salida_delta, out, medianas_capa1);
        return out;
    }

private:
    void guardar_en_historia(
        double y_cruda,
        bool ok_cruda,
        double y_purgada,
        bool ok_purgada,
        double dnu,
        const std::array<double, N_CAPA1>& capa1) noexcept
    {
        const std::size_t pos = (cab_hist_ + n_hist_) % HIST_MAX;
        hist_cruda_[pos] = y_cruda;
        hist_purgada_[pos] = y_purgada;
        hist_dnu_[pos] = dnu;
        hist_val_cruda_[pos] = ok_cruda ? 1u : 0u;
        hist_val_purgada_[pos] = ok_purgada ? 1u : 0u;
        hist_capa1_[pos] = capa1;
        if (n_hist_ < HIST_MAX) ++n_hist_;
        else cab_hist_ = (cab_hist_ + 1) % HIST_MAX;
    }

    // Convierte un paso nativo t_buscado <= t_actual_ en su indice dentro del anillo si aun vive en HIST_MAX.
    [[nodiscard]] bool indice_en_anillo(std::int64_t t_buscado, std::size_t& idx_out) const noexcept {
        if (t_buscado < 0 || t_buscado > t_actual_ || n_hist_ == 0) return false;
        const std::int64_t atras = t_actual_ - t_buscado;
        if (static_cast<std::size_t>(atras) >= n_hist_) return false;
        const std::size_t offset = n_hist_ - 1 - static_cast<std::size_t>(atras);
        idx_out = (cab_hist_ + offset) % HIST_MAX;
        return true;
    }

    // Calcula la mediana causal en [t_ini, t_fin] sin reservar memoria dinamica.
    [[nodiscard]] double mediana_rango_anillo(
        std::int64_t t_ini,
        std::int64_t t_fin,
        bool usar_purgada,
        std::size_t& n_validos_out,
        double& nu_punto_out) noexcept
    {
        n_validos_out = 0;
        nu_punto_out = NaN;
        if (t_fin < t_ini) return NaN;
        double suma_nu = 0.0;
        std::size_t pasos_reloj = 0;
        for (std::int64_t t = t_ini; t <= t_fin; ++t) {
            std::size_t idx = 0;
            if (!indice_en_anillo(t, idx)) continue;
            suma_nu += hist_dnu_[idx];
            ++pasos_reloj;
            const bool val = usar_purgada ? (hist_val_purgada_[idx] != 0u) : (hist_val_cruda_[idx] != 0u);
            const double v = usar_purgada ? hist_purgada_[idx] : hist_cruda_[idx];
            if (val && std::isfinite(v) && n_validos_out < HIST_MAX) {
                scratch_mediana_[n_validos_out++] = v;
            }
        }
        if (pasos_reloj > 0) nu_punto_out = suma_nu / static_cast<double>(pasos_reloj);
        if (n_validos_out == 0) return NaN;
        const std::size_t mid = n_validos_out / 2;
        std::nth_element(scratch_mediana_.begin(), scratch_mediana_.begin() + mid, scratch_mediana_.begin() + n_validos_out);
        const double m_hi = scratch_mediana_[mid];
        if ((n_validos_out % 2u) == 1u) return m_hi;
        double m_lo = scratch_mediana_[0];
        for (std::size_t i = 1; i < mid; ++i) {
            if (scratch_mediana_[i] > m_lo) m_lo = scratch_mediana_[i];
        }
        return 0.5 * (m_lo + m_hi);
    }

    // Extrae eventos de nivel (y modula calor con temperatura) desde una FichaEvento o una salida del operador delta.
    void procesar_ficha_evento(const FichaEvento& fev, FichaEvento& ev_nivel_out, bool& hay_ev_nivel_out) noexcept {
        if (!fev.alarma.hay) return;
        if (fev.canal == 1 || fev.estado == EstadoSaltos::CAMBIO_TEMPERATURA) {
            ++eventos_temp_recibidos_;
            // En canal_temperatura.hpp: signo < 0 significa que beta baja (se calienta); signo > 0 que beta sube (se enfria).
            if (fev.alarma.signo < 0) calor_activo_ = true;
            else if (fev.alarma.signo > 0) calor_activo_ = false;
            return;
        }
        if (fev.canal == 0 &&
            (fev.estado == EstadoSaltos::SALTO ||
             fev.estado == EstadoSaltos::TRANSITORIO ||
             fev.estado == EstadoSaltos::DERIVA)) {
            ++eventos_nivel_recibidos_;
            if (!hay_ev_nivel_out || fev.alarma.g > ev_nivel_out.alarma.g) {
                ev_nivel_out = fev;
                hay_ev_nivel_out = true;
            }
        }
    }

    template <class SalidaObs>
    void extraer_eventos_de_salida_delta(const SalidaObs& s, FichaEvento& ev_nivel_out, bool& hay_ev_nivel_out) noexcept {
        hay_ev_nivel_out = false;
        if constexpr (requires { s.hay_evento; s.evento; }) {
            if (s.hay_evento) procesar_ficha_evento(s.evento, ev_nivel_out, hay_ev_nivel_out);
        } else if constexpr (requires { s.n_eventos; s.eventos; }) {
            for (std::size_t i = 0; i < s.n_eventos; ++i) {
                procesar_ficha_evento(s.eventos[i], ev_nivel_out, hay_ev_nivel_out);
            }
        } else if constexpr (requires { s.alarma; s.estado; s.canal; }) {
            procesar_ficha_evento(s, ev_nivel_out, hay_ev_nivel_out);
        }
    }

    void configurar_rama_en_apertura(
        RamaTrayectoriaLilyth<T_MAX>& rama,
        bool usar_purgada,
        bool habilitada_por_modo,
        bool purga_madura,
        std::int64_t t_ini_ancla,
        std::int64_t t_fin_ancla,
        double y_actual,
        bool y_actual_ok,
        int signo_alarma) noexcept
    {
        rama.limpiar();
        if (!habilitada_por_modo) {
            rama.desenlace = DesenlaceRanura::INACTIVA;
            return;
        }
        if (usar_purgada && (!purga_madura || nivel_purgado_roto_)) {
            rama.desenlace = DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA;
            return;
        }
        if (!y_actual_ok || !std::isfinite(y_actual)) {
            rama.desenlace = DesenlaceRanura::CENSURADA_DATO_INVALIDO;
            return;
        }
        std::size_t n_val = 0;
        double nu_punto = NaN;
        const double y_star = mediana_rango_anillo(t_ini_ancla, t_fin_ancla, usar_purgada, n_val, nu_punto);
        if (n_val == 0 || !std::isfinite(y_star)) {
            rama.desenlace = DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE;
            return;
        }
        const double dy = y_actual - y_star;
        // E1.4: dY = Y(t_a) - Y* debe tener el mismo signo que la alarma y ser no nulo.
        if (!( std::isfinite(dy) && (dy * static_cast<double>(signo_alarma) > 0.0) )) {
            rama.desenlace = DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE;
            return;
        }
        rama.desenlace = DesenlaceRanura::ABIERTA;
        rama.y_ancla = y_star;
        rama.y_disparo = y_actual;
        rama.dy = dy;
        rama.nu_punto_ancla = nu_punto;
        rama.nu_en_disparo = nu_acumulado_;
        rama.muestras_ancla = n_val;
    }

    void abrir_ranura_desde_evento(
        const FichaEvento& ev,
        double y_cruda_actual,
        bool ok_cruda,
        double y_purgada_actual,
        bool ok_purgada,
        bool purga_madura,
        SalidaPaso& out_lilyth) noexcept
    {
        std::size_t slot = RANURAS;
        for (std::size_t i = 0; i < RANURAS; ++i) {
            if (!ranuras_[i].ocupada) { slot = i; break; }
        }
        if (slot == RANURAS) {
            ++desbordes_ranura_;   // D8: desborde contado y declarado, cero reservas.
            return;
        }

        auto& r = ranuras_[slot];
        r.limpiar();
        r.ocupada = true;
        r.id_evento = ++seq_eventos_;
        r.t_apertura = t_actual_;
        r.t_cambio = (ev.alarma.t_cambio >= 0 && ev.alarma.t_cambio <= t_actual_) ? ev.alarma.t_cambio : t_actual_;
        r.signo = ev.alarma.signo;
        r.escala_origen = ev.escala_origen >= 0 ? ev.escala_origen : 0;
        r.amplitud_glr = ev.alarma.amplitud;
        r.se_amplitud_glr = ev.alarma.se_amplitud;
        r.g_disparo = ev.alarma.g;
        r.ficha_evento = ev;

        const std::size_t m_j0 = (ev.ficha_origen.m >= 1)
                                     ? static_cast<std::size_t>(ev.ficha_origen.m)
                                     : (std::size_t{1} << static_cast<std::size_t>(r.escala_origen));
        const std::size_t n_reg = (ev.ficha_origen.det.L >= 2) ? ev.ficha_origen.det.L : 120u;
        const double tau_j0 = ( std::isfinite(ev.ficha_origen.nulo.tau_int) && ev.ficha_origen.nulo.tau_int >= 1.0 )
                                  ? ev.ficha_origen.nulo.tau_int
                                  : 1.0;

        r.ancla_saturada_en_l = (ev.alarma.l >= n_reg);
        r.w_ancla = n_reg * m_j0;
        r.guarda_ancla = static_cast<std::size_t>(std::ceil(tau_j0)) * m_j0;
        r.t_fin_ancla = r.t_cambio - static_cast<std::int64_t>(r.guarda_ancla) - 1;
        const std::int64_t t_ini_ancla = r.t_fin_ancla - static_cast<std::int64_t>(r.w_ancla) + 1;

        // E1.12: congelar por valor las medianas de la Capa 1 vigentes en t_fin_ancla (o en t_cambio^-).
        std::size_t idx_c1 = 0;
        if (indice_en_anillo(std::max<std::int64_t>(0, r.t_fin_ancla), idx_c1)) {
            r.medianas_capa1_congeladas = hist_capa1_[idx_c1];
        }

        const bool medir_purgada = (cfg_.representacion == ModoRepresentacionLilyth::PURGADA) ||
                                   (cfg_.representacion == ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA);
        const bool medir_cruda   = (cfg_.representacion == ModoRepresentacionLilyth::CRUDA) ||
                                   (cfg_.representacion == ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA);

        configurar_rama_en_apertura(
            r.purgada, true, medir_purgada, purga_madura, t_ini_ancla, r.t_fin_ancla, y_purgada_actual, ok_purgada, r.signo);
        configurar_rama_en_apertura(
            r.cruda, false, medir_cruda, true, t_ini_ancla, r.t_fin_ancla, y_cruda_actual, ok_cruda, r.signo);

        if (!r.alguna_rama_abierta()) {
            cerrar_y_archivar_ranura(r, out_lilyth);
            return;
        }
        out_lilyth.abrio_ranura = true;
        out_lilyth.indice_ranura_abierta = slot;
    }

    void avanzar_rama(
        RamaTrayectoriaLilyth<T_MAX>& rama,
        std::size_t tau,
        double y_val,
        bool y_ok,
        bool es_purgada) noexcept
    {
        if (rama.desenlace != DesenlaceRanura::ABIERTA) return;
        if (!y_ok || !std::isfinite(y_val)) {
            rama.desenlace = DesenlaceRanura::CENSURADA_DATO_INVALIDO;
            rama.tau_fin = rama.pasos_observados;
            return;
        }
        if (tau < 1 || tau > cfg_.t_max_efectivo || tau > T_MAX) {
            rama.desenlace = std::isfinite(rama.tau_rest_consolidado)
                                 ? DesenlaceRanura::RESTITUIDA
                                 : DesenlaceRanura::CENSURADA_HORIZONTE;
            rama.tau_fin = std::isfinite(rama.tau_rest_consolidado)
                               ? static_cast<std::size_t>(rama.tau_rest_consolidado)
                               : rama.pasos_observados;
            return;
        }

        const double eta_t = (y_val - rama.y_ancla) / rama.dy;
        rama.eta[tau - 1] = eta_t;
        rama.pasos_observados = tau;

        // Acumular en la curva media eta_barra(tau) sin detenerse en 0.5 (E1.5 P3b y E1.8).
        if (es_purgada) {
            suma_eta_purgada_[tau - 1] += eta_t;
            ++cuenta_eta_purgada_[tau - 1];
        } else {
            suma_eta_cruda_[tau - 1] += eta_t;
            ++cuenta_eta_cruda_[tau - 1];
        }

        // E1.8: sobrepaso evaluado desde tau >= 2.
        if (tau >= 2) {
            if (!std::isfinite(rama.min_eta_desde_2) || eta_t < rama.min_eta_desde_2) {
                rama.min_eta_desde_2 = eta_t;
                rama.tau_min_eta = tau;
            }
            if (eta_t < 0.0) rama.armonico_sospechoso = true;   // D10 / 7.1: etiqueta, nunca veto.
        }

        // E1.4 y E1.6: cruce del 50 % (primer toque y consolidado sin calor).
        if (eta_t <= 0.5) {
            const double dnu = nu_acumulado_ - rama.nu_en_disparo;
            if (!std::isfinite(rama.tau_rest_primer_toque)) {
                rama.tau_rest_primer_toque = static_cast<double>(tau);
                rama.delta_nu_primer_toque = dnu;
                if (std::isfinite(rama.nu_punto_ancla) && rama.nu_punto_ancla > 0.0) {
                    rama.tau_act_primer_toque = dnu / rama.nu_punto_ancla;
                }
            }
            if (!calor_activo_ && !std::isfinite(rama.tau_rest_consolidado)) {
                rama.tau_rest_consolidado = static_cast<double>(tau);
                rama.delta_nu_consolidado = dnu;
                if (std::isfinite(rama.nu_punto_ancla) && rama.nu_punto_ancla > 0.0) {
                    rama.tau_act_consolidado = dnu / rama.nu_punto_ancla;
                }
            }
        }

        // Al llegar al final del horizonte t_max_efectivo, la rama se cierra.
        if (tau == cfg_.t_max_efectivo) {
            if (std::isfinite(rama.tau_rest_consolidado)) {
                rama.desenlace = DesenlaceRanura::RESTITUIDA;
                rama.tau_fin = static_cast<std::size_t>(rama.tau_rest_consolidado);
            } else {
                rama.desenlace = DesenlaceRanura::CENSURADA_HORIZONTE;
                rama.tau_fin = tau;
            }
        }
    }

    void cerrar_y_archivar_ranura(RanuraLilyth<T_MAX, N_CAPA1>& r, SalidaPaso& out_lilyth) noexcept {
        // E1.14: calcular diferencial de marea cuando ambas ramas midieron.
        if (std::isfinite(r.cruda.tau_rest_primer_toque) && std::isfinite(r.purgada.tau_rest_primer_toque)) {
            r.delta_tau_marea_primer_toque = r.cruda.tau_rest_primer_toque - r.purgada.tau_rest_primer_toque;
        }
        if (std::isfinite(r.cruda.tau_rest_consolidado) && std::isfinite(r.purgada.tau_rest_consolidado)) {
            r.delta_tau_marea_consolidado = r.cruda.tau_rest_consolidado - r.purgada.tau_rest_consolidado;
        }

        // E1.10: actualizar tablas competitivas de Aalen-Johansen y contadores de exclusion de Prony reescalado (E1.13).
        if (r.purgada.desenlace != DesenlaceRanura::INACTIVA &&
            r.purgada.desenlace != DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA &&
            r.purgada.desenlace != DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE) {
            aj_purgada_.registrar(r.signo, r.purgada.desenlace, r.purgada.tau_fin);
            if (!std::isfinite(r.purgada.tau_rest_consolidado)) ++excluidos_prony_reescalada_purgada_;
        }
        if (r.cruda.desenlace != DesenlaceRanura::INACTIVA &&
            r.cruda.desenlace != DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA &&
            r.cruda.desenlace != DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE) {
            aj_cruda_.registrar(r.signo, r.cruda.desenlace, r.cruda.tau_fin);
            if (!std::isfinite(r.cruda.tau_rest_consolidado)) ++excluidos_prony_reescalada_cruda_;
        }

        if (out_lilyth.ranuras_cerradas_en_paso < RANURAS) {
            out_lilyth.cerradas[out_lilyth.ranuras_cerradas_en_paso++] = r;
        }
        r.limpiar();
    }
};

}  // namespace ark::lilyth