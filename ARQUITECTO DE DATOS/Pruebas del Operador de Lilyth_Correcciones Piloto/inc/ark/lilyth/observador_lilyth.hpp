// observador_lilyth.hpp -- Operador de Lilyth (restitucion tras un evento del operador delta), ETAPAS 1, 2, 3.0 y 3:
// tipos, trayectoria eta(tau), ranuras con ancla congelada (purgada y cruda), regla de continuacion vs replica (E1.21),
// riesgos competitivos de Aalen-Johansen y calibracion/evaluacion de nulos (a) espontaneo y (b) generativo por inyeccion. C++20.
//
// Autoria: Argon (Auditoria Externa) y Code (Jefe de Arquitectos de Codigo, Opus 5.5), 2026-11-05 y 2026-11-06.
// Especificacion: JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md (secciones 1-9, E1.1-E1.22).
#pragma once

#include "ark/dataset/tipos_dataset.hpp"
#include "ark/lilyth/modelo_lineal_nulo.hpp"
#include "ark/lilyth/armonico_lilyth.hpp"
#include "ark/lilyth/fin_calor_lilyth.hpp"
#include "ark/lilyth/puente_historico_ito.hpp"
#include "ark/lilyth/reloj_marea_lilyth.hpp"
#include "ark/lilyth/nulo_lilyth.hpp"
#include "ark/metrologia/escalera_diadica.hpp"
#include "ark/metrologia/observador_saltos.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

namespace ark::lilyth {

using ark::metrologia::ALFA_PROYECTO;
using ark::metrologia::EventoAlarma;
using ark::metrologia::FichaEscala;
using ark::metrologia::FichaEvento;
using ark::metrologia::Integracion;
using ark::metrologia::NaN;

// Code: los codigos de nulo son los de la base (sustitutos.hpp), no una copia.
inline constexpr std::uint32_t LILYTH_NULO_NINGUNO = ark::metrologia::NULO_NINGUNO;
inline constexpr std::uint32_t LILYTH_NULO_PERMUTACION = ark::metrologia::NULO_PERMUTACION;
inline constexpr std::uint32_t LILYTH_NULO_BLOQUES = ark::metrologia::NULO_BLOQUES;
inline constexpr std::uint32_t LILYTH_NULO_SIGNOS = ark::metrologia::NULO_SIGNOS;

// Rasgos del operador delta: N_REG (ventana del ancla, E1.3), B_MAX (espejo de la purga) y J_MAX. Code: sale del tipo.
template <class Obs>
struct RasgosDelta;
template <std::size_t J_MAX, std::size_t L, std::size_t N_REG, std::size_t B_MAX, std::size_t EV_MAX>
struct RasgosDelta<ark::metrologia::ObservadorSaltos<J_MAX, L, N_REG, B_MAX, EV_MAX>> {
    static constexpr std::size_t j_max = J_MAX, l_glr = L, n_reg = N_REG, b_max = B_MAX;
};

// E1.14: representacion sobre la que mide Lilyth (entrada).
enum class ModoRepresentacionLilyth : std::uint8_t { PURGADA = 0, CRUDA = 1, AMBAS_CONTRASTE_MAREA = 2 };
#ifndef ARK_LILYTH_MODO_NULO_LINEAL_DEF
#define ARK_LILYTH_MODO_NULO_LINEAL_DEF
enum class ModoNuloLineal : std::uint8_t { AR_LEVINSON_DURBIN = 0, ARFIMA_WHITTLE = 1, FASE_ALEATORIZADA_FFT = 2 };
#endif
// E1.8, E1.13: politica del detector armonico (entrada; el contraste es de la Etapa 4).
enum class ModoArmonico : std::uint8_t { DESACTIVADO = 0, PRONY_CURVA_MEDIA = 1, PRONY_FASE_REESCALADA = 2, SOBREPASO_POR_EVENTO = 3, TODOS = 4 };

// E1.10, E1.14: desenlace de cada rama. Code: tres motivos nuevos de NO_MEDIDA para que cada perdida se cuente por su causa.
enum class DesenlaceRanura : std::uint8_t {
    INACTIVA = 0,
    ABIERTA = 1,
    RESTITUIDA = 2,                       // cruzo eta <= 0.5 sin calor (consolidado)
    CORTADA_POR_REPLICA = 3,              // una nueva apertura del operador delta llego antes de restituir
    CENSURADA_HORIZONTE = 4,              // llego a T_max sin restituir
    CENSURADA_DATO_INVALIDO = 5,          // barra invalida, o en NIVEL un incremento purgado invalido (E1.17)
    NO_MEDIDA_PURGA_INMADURA = 6,         // la purga no estaba madura (dias_minimos por fase) en t_a, o el operador delta no purga
    NO_MEDIDA_ANCLA_INSUFICIENTE = 7,     // no hay ancla: historia anterior al inicio, sin muestras validas o con un corte de segmento
    NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA = 8, // Code: W_ancla + guarda no cabe en HIST_MAX (antes se truncaba en silencio)
    NO_MEDIDA_SIGNO_INCONSISTENTE = 9,    // Code: dY = Y(t_a) - Y* no tiene el signo de la alarma (antes se contaba como ancla insuficiente)
    NO_MEDIDA_PERFIL_CORRUPTO = 10        // Piloto Nivel 1: perfil [t_cambio, t_a] con racha invalida > r_max_z2_piloto
};
inline constexpr std::array<const char*, 11> NOMBRE_DESENLACE{
    "INACTIVA", "ABIERTA", "RESTITUIDA", "CORTADA_POR_REPLICA", "CENSURADA_HORIZONTE", "CENSURADA_DATO_INVALIDO",
    "NO_MEDIDA_PURGA_INMADURA", "NO_MEDIDA_ANCLA_INSUFICIENTE", "NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA", "NO_MEDIDA_SIGNO_INCONSISTENTE",
    "NO_MEDIDA_PERFIL_CORRUPTO"};

[[nodiscard]] constexpr bool es_no_medida(DesenlaceRanura d) noexcept {
    return d == DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA || d == DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE ||
           d == DesenlaceRanura::NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA || d == DesenlaceRanura::NO_MEDIDA_SIGNO_INCONSISTENTE ||
           d == DesenlaceRanura::NO_MEDIDA_PERFIL_CORRUPTO;
}

// Configuracion de Lilyth (E1.13, E1.14, E1.15). Code: la integracion y la purga NO estan aqui: salen de la ConfiguracionObservador
// del operador delta (una sola fuente; con dos podian no coincidir).
struct ConfiguracionLilyth {
    ModoRepresentacionLilyth representacion = ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA;
    std::uint32_t nulo_sin_memoria = LILYTH_NULO_PERMUTACION;   // E1.2, E1.18: solo PERMUTACION o SIGNOS
    ModoNuloLineal modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
    ModoArmonico modo_armonico = ModoArmonico::PRONY_CURVA_MEDIA;
    double alfa = ALFA_PROYECTO;
    // E1.9: lo deriva `calibrar` en la Etapa 3 (cuantil 1 - alfa del nulo lineal, o W_ancla con raiz unitaria).
    std::size_t t_max_efectivo = 0;
    bool cortar_abiertas_con_replica = true;   // D8, E1.10

    // Parametros de la Etapa 3 (calibracion de nulos y Monte Carlo del resumen):
    std::size_t b_max_nulo = 16;               // tope B_max para familia::realizaciones en via (a)
    std::size_t hilos_nulo = 0;                // 0 = hardware_concurrency
    // Code (revision D6): presupuesto de Monte Carlo (conjuntos del resumen y bootstrap de P3b y de la marea), entrada declarada.
    // Por omision 1999, el valor del prerregistro de Argon (el codigo entregado usaba 499). Tiene que resolver el primer escalon de Holm
    // de la familia de 4 de E1.16: 1 / (B + 1) <= alfa / 4.
    std::size_t b_mc_resumen = 1999;
    std::uint64_t semilla_nulo = 20261106;
    double d_whittle_purgada = NaN;            // d externa del Whittle exacto cuando modo_nulo_lineal == ARFIMA_WHITTLE
    double d_whittle_cruda = NaN;
    bool forzar_via_b = false;                 // para el control P2b: fuerza la via (b) aunque (a) alcance m_min
    ark::ventanas::TipoRelojVentana reloj_actividad = ark::ventanas::TipoRelojVentana::PARED; // Etapa 6 (E1.11, D7): nunca NU_INFO_SHANNON
    // Code (2026-11-15, orden del Comandante): rama del contraste del resumen (Delta tau y sobrepaso por evento). Por omision la de Argon,
    // identica bit a bit; INTERCAMBIABLE la corrige con pools pequenos y deja la vieja al lado en la ficha (ver nulo_lilyth.hpp).
    ContrasteResumenLilyth contraste_resumen = ContrasteResumenLilyth::POOL_LOO;
    // Code (2026-11-16, orden del Comandante): tamano del pool de nulos emparejados (ver nulo_lilyth.hpp). Por omision el de hoy.
    PoolNuloLilyth pool_nulo = PoolNuloLilyth::M_MIN;
    // Code (2026-11-16, orden del Comandante): horizonte T_max (enmienda de E1.9; ver nulo_lilyth.hpp). Por omision el de Argon.
    HorizonteLilyth horizonte = HorizonteLilyth::CUANTIL_ALFA;

    // Piloto Nivel 1 (Huecos y cotas por zona Z1, Z2, Z3; por omision 0 = comportamiento base bit a bit):
    enum class PoliticaZ1Piloto : std::uint8_t { COMPRESION_BASE = 0, CORTE_SUFIJO = 1, PUENTE_CAUSAL = 2, PUENTE_ITO_HISTORICO = 3 };
    enum class PoliticaZ2Piloto : std::uint8_t { COMPRESION_BASE = 0, PUENTE_CAUSAL = 1, RECHAZAR_SI_CORRUPTO = 2 };
    enum class PoliticaZ3Piloto : std::uint8_t { CORTA_AL_PRIMERO = 0, MUTANTE_COMPRIME_TAU = 1, RELOJ_FISICO_TOPE_R = 2 };

    PoliticaZ1Piloto politica_z1_piloto = PoliticaZ1Piloto::COMPRESION_BASE;
    std::size_t r_max_z1_piloto = 0;
    PoliticaZ2Piloto politica_z2_piloto = PoliticaZ2Piloto::COMPRESION_BASE;
    std::size_t r_max_z2_piloto = 0;
    PoliticaZ3Piloto politica_z3_piloto = PoliticaZ3Piloto::CORTA_AL_PRIMERO;
    std::size_t r_max_z3_piloto = 0;
    double z_max_hueco_piloto = 3.5;
    bool usar_protocolo_filas_piloto = false;
};

[[nodiscard]] constexpr bool configuracion_lilyth_nulo_valida(const ConfiguracionLilyth& cfg) noexcept {
    const bool nulo_ok = cfg.nulo_sin_memoria == LILYTH_NULO_PERMUTACION || cfg.nulo_sin_memoria == LILYTH_NULO_SIGNOS;
    const bool alfa_ok = cfg.alfa > 0.0 && cfg.alfa < 1.0;
    const bool reloj_ok = detalle_reloj_marea::reloj_actividad_valido(cfg.reloj_actividad);
    const bool bmc_ok = alfa_ok && static_cast<double>(cfg.b_mc_resumen + 1) * (cfg.alfa / 4.0) >= 1.0;
    const bool contraste_ok = cfg.contraste_resumen == ContrasteResumenLilyth::POOL_LOO || cfg.contraste_resumen == ContrasteResumenLilyth::INTERCAMBIABLE;
    const bool pool_ok = cfg.pool_nulo == PoolNuloLilyth::M_MIN || cfg.pool_nulo == PoolNuloLilyth::DERIVADO;
    const bool hor_ok = cfg.horizonte == HorizonteLilyth::CUANTIL_ALFA || cfg.horizonte == HorizonteLilyth::CUANTIL_RESOLUCION;
    return nulo_ok && alfa_ok && reloj_ok && bmc_ok && contraste_ok && pool_ok && hor_ok;
}

[[nodiscard]] constexpr bool configuracion_lilyth_valida(const ConfiguracionLilyth& cfg, std::size_t t_max_cap) noexcept {
    const bool tmax_ok = cfg.t_max_efectivo >= 2 && cfg.t_max_efectivo <= t_max_cap;
    return configuracion_lilyth_nulo_valida(cfg) && tmax_ok;
}

// ------------------------------------------------------------------------------------------------- rama y ranura (Argon)

template <std::size_t T_MAX>
struct RamaTrayectoriaLilyth {
    static_assert(T_MAX >= 2, "T_MAX >= 2");
    DesenlaceRanura desenlace = DesenlaceRanura::INACTIVA;
    double y_ancla = NaN;                 // Y* (E1.3)
    double y_disparo = NaN;               // Y(t_a)
    double dy = NaN;                      // dY = Y(t_a) - Y* (E1.4)
    double nu_punto_ancla = NaN;          // nu por barra en W_ancla (E1.11); NaN sin reloj de actividad
    double nu_en_disparo = NaN;           // nu acumulado en t_a
    double pared_en_disparo = NaN;        // Code: marca de pared en t_a (ms en barras intrinsecas; NaN en el camino de serie)
    std::size_t muestras_ancla = 0;
    std::size_t pasos_observados = 0;
    std::size_t tau_fin = 0;

    double tau_rest_primer_toque = NaN;   // primer tau >= 1 con eta <= 0.5 (con o sin calor)
    double tau_rest_consolidado = NaN;    // primer tau >= 1 con eta <= 0.5 sin calor
    double delta_nu_primer_toque = NaN, delta_nu_consolidado = NaN;
    double tau_act_primer_toque = NaN, tau_act_consolidado = NaN;   // delta_nu / nu_punto_ancla (barras de calma)
    double pared_primer_toque = NaN, pared_consolidado = NaN;       // Code: tiempo de pared transcurrido (P4)

    double min_eta_desde_2 = NaN;
    std::size_t tau_min_eta = 0;
    bool armonico_sospechoso = false;
    std::size_t racha_invalida_z3 = 0;
    std::size_t pasos_validos_z3 = 0;

    std::array<double, T_MAX> eta{};      // eta(tau), indice tau - 1

    // Etapa 3: referencias contra el nulo sin memoria (iid) y el nulo lineal (Onsager) llenadas por completar_nulos().
    ReferenciaNuloEvento nulo_iid{};
    ReferenciaNuloEvento nulo_onsager{};

    void limpiar() noexcept {
        desenlace = DesenlaceRanura::INACTIVA;
        y_ancla = y_disparo = dy = nu_punto_ancla = nu_en_disparo = pared_en_disparo = NaN;
        muestras_ancla = pasos_observados = tau_fin = 0;
        tau_rest_primer_toque = tau_rest_consolidado = NaN;
        delta_nu_primer_toque = delta_nu_consolidado = NaN;
        tau_act_primer_toque = tau_act_consolidado = NaN;
        pared_primer_toque = pared_consolidado = NaN;
        min_eta_desde_2 = NaN;
        tau_min_eta = 0;
        armonico_sospechoso = false;
        racha_invalida_z3 = 0;
        pasos_validos_z3 = 0;
        eta.fill(NaN);
        nulo_iid = ReferenciaNuloEvento{};
        nulo_onsager = ReferenciaNuloEvento{};
    }
};

template <std::size_t T_MAX, std::size_t N_CAPA1 = 4>
struct RanuraLilyth {
    bool ocupada = false;
    std::uint64_t id_evento = 0;
    std::int64_t t_apertura = -1;          // t_a = alarma.t (paso nativo del operador delta)
    std::int64_t t_cambio = -1;
    std::int64_t t_fin_ancla = -1;         // t_cambio - guarda - 1
    std::size_t w_ancla = 0;               // N_REG * 2^j0 (E1.3)
    std::size_t guarda_ancla = 0;          // ceil(tau_int,j0) * 2^j0 (E1.3)
    int signo = 0;
    int escala_origen = 0;
    double amplitud_glr = NaN, se_amplitud_glr = NaN, g_disparo = NaN;   // E1.4: guardadas aparte
    bool ancla_saturada_en_l = false;      // la memoria elegida por el GLR llego a su tope L
    int signo_replica_corte = 0;           // signo de la apertura que corto la ranura (0 si no la corto una replica)
    bool nulo_por_inyeccion = false;       // Etapa 3 (seccion 3.3): true si el nulo lineal del evento se completo por via (b)
    std::array<double, N_CAPA1> medianas_capa1_congeladas{};   // E1.12: copia por valor en t_fin_ancla
    FichaEvento ficha_apertura{};          // Code: la ficha tal como el operador delta la arma al ABRIR (alarma, escala, ficha de la escala)
    RamaTrayectoriaLilyth<T_MAX> purgada{};
    RamaTrayectoriaLilyth<T_MAX> cruda{};
    double delta_tau_marea_primer_toque = NaN;   // E1.14: tau_cruda - tau_purgada
    double delta_tau_marea_consolidado = NaN;

    // Etapa 3: perfil estacionario [t_cambio, t_a] y su historia previa guardados en pila al abrir (para la via (b) fuera del ciclo caliente).
    std::size_t n_hist_prev_cruda = 0, n_perfil_cruda = 0;
    std::size_t n_hist_prev_purgada = 0, n_perfil_purgada = 0;
    bool perfil_corrupto_z2_cruda = false, perfil_corrupto_z2_purgada = false;
    std::array<double, MAX_PERFIL_EVENTO_LILYTH> hist_prev_cruda_est{}, perfil_cruda_est{};
    std::array<double, MAX_PERFIL_EVENTO_LILYTH> hist_prev_purgada_est{}, perfil_purgada_est{};

    void limpiar() noexcept {
        ocupada = false;
        id_evento = 0;
        t_apertura = t_cambio = t_fin_ancla = -1;
        w_ancla = guarda_ancla = 0;
        signo = escala_origen = signo_replica_corte = 0;
        amplitud_glr = se_amplitud_glr = g_disparo = NaN;
        ancla_saturada_en_l = false;
        nulo_por_inyeccion = false;
        medianas_capa1_congeladas.fill(NaN);
        ficha_apertura = FichaEvento{};
        purgada.limpiar();
        cruda.limpiar();
        delta_tau_marea_primer_toque = delta_tau_marea_consolidado = NaN;
        n_hist_prev_cruda = n_perfil_cruda = n_hist_prev_purgada = n_perfil_purgada = 0;
        perfil_corrupto_z2_cruda = perfil_corrupto_z2_purgada = false;
    }
    [[nodiscard]] bool alguna_rama_abierta() const noexcept {
        return purgada.desenlace == DesenlaceRanura::ABIERTA || cruda.desenlace == DesenlaceRanura::ABIERTA;
    }
    [[nodiscard]] bool ha_restituido() const noexcept {
        return (purgada.desenlace == DesenlaceRanura::ABIERTA && std::isfinite(purgada.tau_rest_consolidado)) ||
               (cruda.desenlace == DesenlaceRanura::ABIERTA && std::isfinite(cruda.tau_rest_consolidado));
    }
};

// ------------------------------------------------------------------------------------------------- Aalen-Johansen (Argon, sin cambios)
template <std::size_t T_MAX>
struct CurvaAalenJohansen {
    std::size_t n_eventos_registrados = 0, n_restituidos = 0, n_cortados_replica = 0, n_censurados = 0;
    std::size_t t_max_efectivo = T_MAX;
    std::array<std::size_t, T_MAX> d_rest{}, d_replica{}, censurados_en{};
    std::size_t censurados_en_cero = 0;
    std::array<std::size_t, T_MAX> en_riesgo{};
    std::array<double, T_MAX> supervivientes{}, cif_rest{}, cif_replica{};
    double rmst = NaN;
    double tau_mediano_cif_rest = NaN;

    void limpiar(std::size_t tmax_ef = T_MAX) noexcept {
        n_eventos_registrados = n_restituidos = n_cortados_replica = n_censurados = 0;
        t_max_efectivo = std::min(tmax_ef, T_MAX);
        d_rest.fill(0); d_replica.fill(0); censurados_en.fill(0);
        censurados_en_cero = 0;
        en_riesgo.fill(0);
        supervivientes.fill(NaN); cif_rest.fill(NaN); cif_replica.fill(NaN);
        rmst = tau_mediano_cif_rest = NaN;
    }
    void registrar_rama(DesenlaceRanura desenlace, std::size_t tau_fin) noexcept {
        if (desenlace == DesenlaceRanura::INACTIVA || desenlace == DesenlaceRanura::ABIERTA || es_no_medida(desenlace)) return;
        const std::size_t tf = std::min(tau_fin, t_max_efectivo);
        ++n_eventos_registrados;
        if (desenlace == DesenlaceRanura::RESTITUIDA && tf >= 1) { ++d_rest[tf - 1]; ++n_restituidos; }
        else if (desenlace == DesenlaceRanura::CORTADA_POR_REPLICA && tf >= 1) { ++d_replica[tf - 1]; ++n_cortados_replica; }
        else { ++n_censurados; if (tf == 0) ++censurados_en_cero; else ++censurados_en[tf - 1]; }
        recalcular();
    }
    void recalcular() noexcept {
        if (n_eventos_registrados == 0 || t_max_efectivo == 0) { rmst = tau_mediano_cif_rest = NaN; return; }
        std::size_t n_act = n_eventos_registrados >= censurados_en_cero ? n_eventos_registrados - censurados_en_cero : 0;
        double s_prev = 1.0, c_rest = 0.0, c_rep = 0.0, suma_rmst = 0.0;
        tau_mediano_cif_rest = NaN;
        for (std::size_t idx = 0; idx < t_max_efectivo; ++idx) {
            en_riesgo[idx] = n_act;
            const std::size_t dr = d_rest[idx], dp = d_replica[idx];
            if (n_act > 0) {
                const double denom = static_cast<double>(n_act);
                c_rest += s_prev * (static_cast<double>(dr) / denom);
                c_rep += s_prev * (static_cast<double>(dp) / denom);
                s_prev *= static_cast<double>(n_act - std::min(n_act, dr + dp)) / denom;
            }
            supervivientes[idx] = s_prev; cif_rest[idx] = c_rest; cif_replica[idx] = c_rep;
            suma_rmst += 1.0 - c_rest;
            if (!std::isfinite(tau_mediano_cif_rest) && c_rest >= 0.5) tau_mediano_cif_rest = static_cast<double>(idx + 1);
            const std::size_t bajas = dr + dp + censurados_en[idx];
            n_act = n_act > bajas ? n_act - bajas : 0;
        }
        rmst = suma_rmst;
    }
};

template <std::size_t T_MAX>
struct ResumenAalenJohansenSigno {
    CurvaAalenJohansen<T_MAX> total{}, positivos{}, negativos{};
    void limpiar(std::size_t tmax_ef = T_MAX) noexcept { total.limpiar(tmax_ef); positivos.limpiar(tmax_ef); negativos.limpiar(tmax_ef); }
    void registrar(int signo, DesenlaceRanura d, std::size_t tau_fin) noexcept {
        total.registrar_rama(d, tau_fin);
        if (signo > 0) positivos.registrar_rama(d, tau_fin);
        else if (signo < 0) negativos.registrar_rama(d, tau_fin);
    }
    [[nodiscard]] double asimetria_rmst() const noexcept { return negativos.rmst - positivos.rmst; }
};

// ------------------------------------------------------------------------------------------------- orquestador

template <class Obs, std::size_t RANURAS = 8, std::size_t T_MAX = 512, std::size_t HIST_MAX = 4096, std::size_t N_CAPA1 = 4>
class ObservadorLilyth {
    static_assert(RANURAS >= 1, "RANURAS >= 1");
    static_assert(T_MAX >= 2, "T_MAX >= 2");
    static_assert(HIST_MAX >= 64, "HIST_MAX >= 64");
    using Rasgos = RasgosDelta<Obs>;
    static constexpr std::size_t J_MAX = Rasgos::j_max;

public:
    using Ranura = RanuraLilyth<T_MAX, N_CAPA1>;
    using SalidaDelta = typename Obs::Salida;

    struct SalidaPaso {
        bool emitida = false;
        std::int64_t t = -1;
        bool abrio_ranura = false;
        std::size_t indice_ranura_abierta = 0;
        bool replica = false;
        bool continuacion = false;
        std::size_t ranuras_cerradas = 0;
    };

    struct Contadores {
        std::size_t pasos = 0, pasos_desalineados = 0, aperturas = 0, replicas = 0, continuaciones = 0, desbordes = 0;
        std::size_t cierres_espejo = 0;
        std::size_t espejo_inconsistente = 0;
        std::size_t solapes = 0;              // Etapa 3 (seccion 3.3): inyecciones descartadas por solape de ventana
        std::size_t historico_desbordado = 0; // Code: ranuras cerradas que no cupieron en el historico (tope de memoria)
        std::size_t no_medida[10]{};
    };

private:
    ConfiguracionLilyth cfg_{};
    ark::metrologia::ConfiguracionObservador cfg_delta_{};
    bool lista_ = false;
    Integracion integracion_ = Integracion::FLUJO;
    bool purga_activa_ = false;
    std::size_t escalas_ = 0;
    std::array<FichaEscala, J_MAX> fichas_{};
    ark::metrologia::PurgaDiurna<Rasgos::b_max> purga_{};
    ark::metrologia::PurgaDiurna<Rasgos::b_max> purga_calibrada_{};

    // historia (anillo) indexada por el t del operador delta
    std::array<double, HIST_MAX> hist_cruda_{}, hist_purgada_{}, hist_dnu_{};
    std::array<double, HIST_MAX> hist_est_cruda_{}, hist_est_purgada_{};
    std::array<std::uint8_t, HIST_MAX> hist_val_cruda_{}, hist_val_purgada_{};
    std::array<std::uint8_t, HIST_MAX> hist_val_est_cruda_{}, hist_val_est_purgada_{};
    std::array<std::uint8_t, HIST_MAX> hist_flags_hueco_{};
    std::array<std::uint32_t, HIST_MAX> hist_seg_purgada_{};
    std::array<std::array<double, N_CAPA1>, HIST_MAX> hist_capa1_{};
    std::array<double, HIST_MAX> scratch_{};
    std::size_t n_hist_ = 0;
    std::int64_t t_ = -1;

    // estado de representacion
    double x_prev_ = NaN;
    bool hay_prev_ = false;
    std::uint8_t flags_prev_ = 0;
    double nivel_purgado_ = NaN;
    bool purgado_vivo_ = false;
    std::uint32_t seg_purgado_ = 0;
    double nu_acum_ = 0.0;
    bool hay_reloj_ = false;
    double pared_ = NaN;

    bool calor_activo_ = false;
    GestorFinCalorLilyth<J_MAX, Rasgos::l_glr> gestor_calor_{};

    // espejo del evento de nivel del operador delta
    bool espejo_abierto_ = false;
    EventoAlarma espejo_alarma_{};
    int espejo_escala_ = -1;

    std::array<Ranura, RANURAS> ranuras_{};
    std::array<Ranura, RANURAS> archivo_{};
    std::uint64_t seq_ = 0;
    Contadores cont_{};
    ResumenAalenJohansenSigno<T_MAX> aj_purgada_{}, aj_cruda_{};
    std::array<double, T_MAX> suma_eta_purgada_{}, suma_eta_cruda_{};
    std::array<std::size_t, T_MAX> cuenta_eta_purgada_{}, cuenta_eta_cruda_{};
    std::size_t excluidos_reescalada_purgada_ = 0, excluidos_reescalada_cruda_ = 0;

    // Etapa 3: estado de calibracion de nulos y registro pre-reservado de ranuras cerradas para completar_nulos()
    bool calibrada_nulo_ = false;
    Obs delta_calibrado_{};
    FichaNuloRepresentacion<T_MAX> nulo_purgada_{}, nulo_cruda_{};
    detalle_nulo::ContextoPurgaNulo<Rasgos::b_max> ctx_purga_{};   // Code (revision D7)
    // Code (revision D1c): nulos EMPAREJADOS de cada evento real (indice del historico), los mismos de la Etapa 3; el armonico los usa
    // en vez de un pool filtrado por SNR y "enriquecido" con eventos fabricados.
    std::vector<std::vector<EventoNuloLilyth<T_MAX>>> emp_iid_cruda_{}, emp_ons_cruda_{}, emp_iid_purgada_{}, emp_ons_purgada_{};
    double x0_nivel_ = NaN;                                          // Code: primer x valido de la calibracion (NIVEL)
    std::size_t n_cal_ = 0;
    std::size_t espejo_purga_discrepancias_ = 0;                     // Code: control C-D7a (debe ser 0)
    ResumenArmonicoYHolm<T_MAX> armonico_purgada_{}, armonico_cruda_{};
    ResumenMareaLilyth resumen_marea_{};
    std::vector<Ranura> cerradas_historico_{};

    [[nodiscard]] bool medir_purgada() const noexcept {
        return cfg_.representacion != ModoRepresentacionLilyth::CRUDA;
    }
    [[nodiscard]] bool medir_cruda() const noexcept {
        return cfg_.representacion != ModoRepresentacionLilyth::PURGADA;
    }

public:
    ObservadorLilyth() = default;

    // Etapa 1: configurar con t_max_efectivo ya fijado.
    bool configurar(const ConfiguracionLilyth& cfg, const ark::metrologia::ConfiguracionObservador& cfg_delta, const Obs& delta,
                    std::span<const double> x_cal, std::span<const std::uint8_t> valido_cal = {},
                    std::span<const std::size_t> fase_cal = {}, std::span<const std::int64_t> dia_cal = {}) {
        cfg_ = cfg;
        cfg_delta_ = cfg_delta;
        lista_ = configuracion_lilyth_valida(cfg_, T_MAX);
        integracion_ = cfg_delta.integracion;
        const auto ag = integracion_ == Integracion::NIVEL ? ark::metrologia::AgregacionEscala::SUMA : cfg_delta.agregacion;
        purga_activa_ = cfg_delta.purga && ag != ark::metrologia::AgregacionEscala::RAZON_DE_SUMAS;
        purga_ = ark::metrologia::PurgaDiurna<Rasgos::b_max>(cfg_delta.bins_fase, cfg_delta.dias_minimos_purga);
        if (purga_activa_) {
            bool hp = false;
            double xp = NaN;
            for (std::size_t i = 0; i < x_cal.size(); ++i) {
                const bool ok = (valido_cal.empty() || valido_cal[i] != 0) && std::isfinite(x_cal[i]);
                double y = NaN;
                bool v = ok;
                if (integracion_ == Integracion::NIVEL) {
                    v = ok && hp;
                    y = v ? x_cal[i] - xp : NaN;
                    hp = ok;
                    if (ok) xp = x_cal[i];
                } else {
                    y = ok ? x_cal[i] : NaN;
                }
                bool vp = false;
                (void)purga_.paso(v ? y : NaN, fase_cal.empty() ? 0 : fase_cal[i], dia_cal.empty() ? 0 : dia_cal[i], vp);
            }
        }
        purga_calibrada_ = purga_;
        escalas_ = std::min(delta.escalas(), J_MAX);
        for (std::size_t j = 0; j < escalas_; ++j) fichas_[j] = delta.ficha_nivel(j);
        gestor_calor_.calibrar(delta, cfg_.t_max_efectivo, cfg_.alfa, ark::metrologia::familia::m_min(cfg_.alfa));   // Code: m_min(alfa), no 8
        cerradas_historico_.reserve(1024);
        reiniciar_flujo();
        return lista_;
    }

    // Etapa 3 (seccion 3.6): `calibrar` llama al `configurar` de la Etapa 1 (con horizonte inicial T_MAX), ajusta los modelos lineales
    // de la Etapa 2, corre la via (a), fija T_max (E1.9) en `cfg_.t_max_efectivo` y guarda las poblaciones nulas.
    bool calibrar(const Obs& delta, const ConfiguracionLilyth& cfg, const ark::metrologia::ConfiguracionObservador& cfg_delta,
                  std::span<const double> x_cal, std::span<const std::uint8_t> valido_cal = {},
                  std::span<const std::size_t> fase_cal = {}, std::span<const std::int64_t> dia_cal = {},
                  std::span<const double> nu_cal = {}, std::span<const std::uint8_t> fuente_nulo = {}) {
        (void)nu_cal;
        calibrada_nulo_ = false;
        espejo_purga_discrepancias_ = 0;
        nulo_purgada_ = FichaNuloRepresentacion<T_MAX>{};
        nulo_cruda_ = FichaNuloRepresentacion<T_MAX>{};
        armonico_purgada_ = ResumenArmonicoYHolm<T_MAX>{};
        armonico_cruda_ = ResumenArmonicoYHolm<T_MAX>{};
        resumen_marea_ = ResumenMareaLilyth{};
        if (!configuracion_lilyth_nulo_valida(cfg) || x_cal.size() < 16 || delta.escalas() == 0 ||
            (!fuente_nulo.empty() && fuente_nulo.size() != x_cal.size())) {
            lista_ = false;
            return false;
        }

        // 1. Configurar inicialmente con t_max_efectivo = T_MAX para medir las realizaciones de la via (a) hasta el tope T_MAX.
        ConfiguracionLilyth cfg_init = cfg;
        cfg_init.t_max_efectivo = T_MAX;
        if (!configurar(cfg_init, cfg_delta, delta, x_cal, valido_cal, fase_cal, dia_cal)) {
            return false;
        }
        delta_calibrado_ = delta;

        // 2. Extraer los segmentos estacionarios de calma para CRUDA y PURGADA sobre la historia de calibracion.
        std::vector<std::vector<double>> segs_cruda, segs_purg;
        {
            ark::metrologia::PurgaDiurna<Rasgos::b_max> pu(cfg_delta.bins_fase, cfg_delta.dias_minimos_purga);
            // Code (revision D7): espejo invertible de la misma purga, alimentado con lo mismo; control C-D7a (identidad bit a bit).
            detalle_nulo::EspejoPurgaNulo<Rasgos::b_max> esp(cfg_delta.bins_fase, cfg_delta.dias_minimos_purga);
            ctx_purga_ = detalle_nulo::ContextoPurgaNulo<Rasgos::b_max>{};
            ctx_purga_.purga_activa = purga_activa_;
            x0_nivel_ = NaN;
            bool hp = false;
            double xp = NaN;
            std::vector<double> cur_c, cur_p;
            // Code (2026-11-16): con mascara de fuente, un 0 corta el segmento del nulo (la purga y lo demas no cambian).
            std::size_t nf_c = 0, nf_p = 0;
            for (std::size_t i = 0; i < x_cal.size(); ++i) {
                const bool en_fuente = fuente_nulo.empty() || fuente_nulo[i] != 0;
                const bool ok = (valido_cal.empty() || valido_cal[i] != 0) && std::isfinite(x_cal[i]);
                double y = NaN;
                bool v = ok;
                if (integracion_ == Integracion::NIVEL) {
                    v = ok && hp;
                    y = v ? x_cal[i] - xp : NaN;
                    hp = ok;
                    if (ok) xp = x_cal[i];
                } else {
                    y = ok ? x_cal[i] : NaN;
                }
                if (v && std::isfinite(y) && en_fuente) {
                    cur_c.push_back(y);
                    ++nf_c;
                } else if (!cur_c.empty()) {
                    segs_cruda.push_back(std::move(cur_c));
                    cur_c.clear();
                }

                if (ok && !std::isfinite(x0_nivel_)) x0_nivel_ = x_cal[i];
                if (purga_activa_) {
                    bool vp = false;
                    bool ve = false;
                    const std::size_t fa_i = fase_cal.empty() ? 0 : fase_cal[i];
                    const std::int64_t di_i = dia_cal.empty() ? 0 : dia_cal[i];
                    const double yp = pu.paso(v ? y : NaN, fa_i, di_i, vp);
                    const double ye = esp.paso(v ? y : NaN, fa_i, di_i, ve);
                    if (ve != vp || !((std::isnan(yp) && std::isnan(ye)) || yp == ye)) ++espejo_purga_discrepancias_;
                    if (v && vp && std::isfinite(yp) && en_fuente) {
                        cur_p.push_back(yp);
                        ++nf_p;
                    } else if (!cur_p.empty()) {
                        segs_purg.push_back(std::move(cur_p));
                        cur_p.clear();
                    }
                }
            }
            if (!cur_c.empty()) segs_cruda.push_back(std::move(cur_c));
            if (!cur_p.empty()) segs_purg.push_back(std::move(cur_p));
            if (purga_activa_) {
                ctx_purga_.espejo_calibrado = esp;
                ctx_purga_.fases.assign(fase_cal.begin(), fase_cal.end());
                ctx_purga_.dias.assign(dia_cal.begin(), dia_cal.end());
                if (ctx_purga_.fases.empty()) ctx_purga_.fases.assign(x_cal.size(), 0);
                if (ctx_purga_.dias.empty()) ctx_purga_.dias.assign(x_cal.size(), 0);
                ctx_purga_.desplazamiento_dias = ctx_purga_.dias.back() - ctx_purga_.dias.front() + 1;
            }
            n_cal_ = x_cal.size();
            nulo_cruda_.fuente_externa = nulo_purgada_.fuente_externa = !fuente_nulo.empty();
            nulo_cruda_.n_fuente = nf_c;
            nulo_purgada_.n_fuente = nf_p;
        }

        const std::size_t mm = ark::metrologia::familia::m_min(cfg_.alfa);
        double b_ped = 0.0;
        // Code (revision D6): B = familia::realizaciones con el ARL0 de la familia del operador delta (una falsa alarma cada ARL0 pasos
        // nativos), sin el factor 2 ni el 1000 de respaldo de la entrega (el ARL0 es obligatorio en ConfiguracionObservador).
        if (!(cfg_delta.arl0_pasos_nativos > 1.0)) { lista_ = false; return false; }
        const std::size_t B = ark::metrologia::familia::realizaciones(
            mm, cfg_delta.arl0_pasos_nativos, x_cal.size(), std::max<std::size_t>(1, cfg_.b_max_nulo), b_ped);
        const bool b_top = (b_ped > static_cast<double>(std::max<std::size_t>(1, cfg_.b_max_nulo)));
        const std::size_t w_ancla_j0 = Rasgos::n_reg * static_cast<std::size_t>( std::max<std::int64_t>(1, fichas_[0].m) );

        auto preparar_ficha_rep = [&](FichaNuloRepresentacion<T_MAX>& fr, bool es_purgada,
                                      std::vector<std::vector<double>>&& segs, double d_ext, std::uint64_t salt) {
            fr.activa = true;
            fr.nivel_t_max = nivel_horizonte(cfg_.horizonte, cfg_.alfa);
            fr.m_min = mm;
            fr.B = B;
            fr.b_pedido = b_ped;
            fr.b_topado = b_top;
            fr.segmentos_calma = std::move(segs);
            if (fr.segmentos_calma.empty()) {
                fr.t_max_derivado = T_MAX;
                fr.t_max_censurado = true;
                return true;
            }
            ConfiguracionModeloLineal cfg_m{};
            cfg_m.alfa = cfg_.alfa;
            cfg_m.modo = cfg_.modo_nulo_lineal;
            cfg_m.d_externa = d_ext;
            fr.modelo_lineal = ajustar(fr.segmentos_calma, cfg_m);
            fr.nulo_fft_degradado_a_generativo = fr.modelo_lineal.cayo_a_generativo;
            if (cfg_.modo_nulo_lineal == ModoNuloLineal::ARFIMA_WHITTLE && !fr.modelo_lineal.valido) {
                return false;
            }

            // Correr via (a) para nulo sin memoria y nulo lineal
            std::size_t sf = 0;
            fr.nulos_sin_memoria = detalle_nulo::correr_poblacion_via_a<T_MAX, Rasgos::b_max>(
                delta_calibrado_, cfg_delta_, cfg_, *this, false, es_purgada, fr.segmentos_calma, fr.modelo_lineal,
                ctx_purga_, x0_nivel_, n_cal_, B, cfg_.hilos_nulo, cfg_.semilla_nulo ^ salt, sf);
            fr.sin_factor_purga += sf;

            if (fr.modelo_lineal.valido) {
                fr.nulos_lineal = detalle_nulo::correr_poblacion_via_a<T_MAX, Rasgos::b_max>(
                    delta_calibrado_, cfg_delta_, cfg_, *this, true, es_purgada, fr.segmentos_calma, fr.modelo_lineal,
                    ctx_purga_, x0_nivel_, n_cal_, B, cfg_.hilos_nulo, cfg_.semilla_nulo ^ (salt + 0x5555555555555555ULL), sf);
                fr.sin_factor_purga += sf;
            }

            fr.t_max_derivado = detalle_nulo::derivar_t_max<T_MAX>(
                integracion_, w_ancla_j0, fr.modelo_lineal, fr.nulos_lineal, fr.nivel_t_max,
                fr.t_max_censurado, fr.t_max_raiz_unitaria, fr.fraccion_censurada_lineal);
            fr.calibrada = true;
            return true;
        };

        if (medir_cruda()) {
            if (!preparar_ficha_rep(nulo_cruda_, false, std::move(segs_cruda), cfg_.d_whittle_cruda, 0x1111111111111111ULL)) {
                lista_ = false;
                return false;
            }
        }
        if (medir_purgada() && purga_activa_) {
            if (!preparar_ficha_rep(nulo_purgada_, true, std::move(segs_purg), cfg_.d_whittle_purgada, 0x2222222222222222ULL)) {
                lista_ = false;
                return false;
            }
        }

        // 3. Fijar cfg_.t_max_efectivo (seccion 3.4: deja de ser entrada y lo escribe calibrar).
        std::size_t tmax_final = 0;
        if (nulo_cruda_.calibrada && nulo_purgada_.calibrada) {
            tmax_final = std::max(nulo_cruda_.t_max_derivado, nulo_purgada_.t_max_derivado);
        } else if (nulo_cruda_.calibrada) {
            tmax_final = nulo_cruda_.t_max_derivado;
        } else if (nulo_purgada_.calibrada) {
            tmax_final = nulo_purgada_.t_max_derivado;
        } else {
            tmax_final = T_MAX;
        }
        cfg_.t_max_efectivo = std::clamp<std::size_t>(std::max(tmax_final, cfg.t_max_efectivo), 2, T_MAX);
        // Code (revision D6): las alarmas legibles minimas del fin del calor son m_min(alfa), no un 8 escrito.
        gestor_calor_.calibrar(delta_calibrado_, cfg_.t_max_efectivo, cfg_.alfa, ark::metrologia::familia::m_min(cfg_.alfa));
        lista_ = configuracion_lilyth_valida(cfg_, T_MAX);
        calibrada_nulo_ = lista_;
        reiniciar_flujo_conservando_purga();
        return lista_;
    }

    // Code (revision R13 / D6): se retira la sobrecarga que DEDUCIA la ConfiguracionObservador (NIVEL por la longitud de la historia,
    // ARL0 = arl0_escala / (2 J), 1000 de respaldo). La configuracion del operador delta tiene una sola fuente: se pasa explicita.

    void fijar_calor_activo(bool activo) noexcept { calor_activo_ = activo; }

    void reiniciar_flujo() noexcept {
        n_hist_ = 0; t_ = -1;
        x_prev_ = NaN; hay_prev_ = false; flags_prev_ = 0;
        nivel_purgado_ = NaN; purgado_vivo_ = false; seg_purgado_ = 0;
        nu_acum_ = 0.0; hay_reloj_ = false; pared_ = NaN;
        calor_activo_ = false;
        gestor_calor_.reiniciar_flujo();
        espejo_abierto_ = false; espejo_alarma_ = EventoAlarma{}; espejo_escala_ = -1;
        for (auto& r : ranuras_) r.limpiar();
        const std::size_t tm = lista_ ? cfg_.t_max_efectivo : T_MAX;
        aj_purgada_.limpiar(tm); aj_cruda_.limpiar(tm);
        suma_eta_purgada_.fill(0.0); suma_eta_cruda_.fill(0.0);
        cuenta_eta_purgada_.fill(0); cuenta_eta_cruda_.fill(0);
        excluidos_reescalada_purgada_ = excluidos_reescalada_cruda_ = 0;
        cerradas_historico_.clear();
        resumen_marea_ = ResumenMareaLilyth{};
        cont_ = Contadores{};
    }

    void reiniciar_flujo_conservando_purga() {
        purga_ = purga_calibrada_;
        reiniciar_flujo();
        cerradas_historico_.reserve(1024);
    }

    [[nodiscard]] bool lista() const noexcept { return lista_; }
    [[nodiscard]] bool calibrada_nulo() const noexcept { return calibrada_nulo_; }
    [[nodiscard]] const ConfiguracionLilyth& configuracion() const noexcept { return cfg_; }
    [[nodiscard]] const Contadores& contadores() const noexcept { return cont_; }
    [[nodiscard]] bool calor_activo() const noexcept { return calor_activo_ || gestor_calor_.calor_activo(); }
    [[nodiscard]] const GestorFinCalorLilyth<J_MAX, Rasgos::l_glr>& gestor_calor() const noexcept { return gestor_calor_; }
    [[nodiscard]] GestorFinCalorLilyth<J_MAX, Rasgos::l_glr>& gestor_calor_mutable() noexcept { return gestor_calor_; }
    [[nodiscard]] const std::array<Ranura, RANURAS>& ranuras() const noexcept { return ranuras_; }
    [[nodiscard]] const Ranura& cerrada(std::size_t i) const noexcept { return archivo_[i < RANURAS ? i : 0]; }
    [[nodiscard]] const std::vector<Ranura>& cerradas_historico() const noexcept { return cerradas_historico_; }
    [[nodiscard]] const FichaNuloRepresentacion<T_MAX>& ficha_nulo_purgada() const noexcept { return nulo_purgada_; }
    [[nodiscard]] const FichaNuloRepresentacion<T_MAX>& ficha_nulo_cruda() const noexcept { return nulo_cruda_; }
    [[nodiscard]] std::size_t espejo_purga_discrepancias() const noexcept { return espejo_purga_discrepancias_; }
    // Code (revision D1c): nulos emparejados de cada ranura del historico (indice), los de la Etapa 3; los usa el armonico.
    [[nodiscard]] const std::vector<std::vector<EventoNuloLilyth<T_MAX>>>& nulos_emparejados(bool purgada, bool sin_memoria) const noexcept {
        return purgada ? (sin_memoria ? emp_iid_purgada_ : emp_ons_purgada_) : (sin_memoria ? emp_iid_cruda_ : emp_ons_cruda_);
    }
    [[nodiscard]] const detalle_nulo::ContextoPurgaNulo<Rasgos::b_max>& contexto_purga_nulo() const noexcept { return ctx_purga_; }
    [[nodiscard]] const ResumenArmonicoYHolm<T_MAX>& resumen_armonico_purgada() const noexcept { return armonico_purgada_; }
    [[nodiscard]] const ResumenArmonicoYHolm<T_MAX>& resumen_armonico_cruda() const noexcept { return armonico_cruda_; }
    [[nodiscard]] const ResumenMareaLilyth& resumen_marea() const noexcept { return resumen_marea_; }
    [[nodiscard]] const ResumenAalenJohansenSigno<T_MAX>& aalen_johansen_purgada() const noexcept { return aj_purgada_; }
    [[nodiscard]] const ResumenAalenJohansenSigno<T_MAX>& aalen_johansen_cruda() const noexcept { return aj_cruda_; }
    [[nodiscard]] std::size_t excluidos_reescalada_purgada() const noexcept { return excluidos_reescalada_purgada_; }
    [[nodiscard]] std::size_t excluidos_reescalada_cruda() const noexcept { return excluidos_reescalada_cruda_; }
    [[nodiscard]] double eta_barra_purgada(std::size_t tau) const noexcept {
        return (tau < 1 || tau > T_MAX || !cuenta_eta_purgada_[tau - 1]) ? NaN : suma_eta_purgada_[tau - 1] / static_cast<double>(cuenta_eta_purgada_[tau - 1]);
    }
    [[nodiscard]] double eta_barra_cruda(std::size_t tau) const noexcept {
        return (tau < 1 || tau > T_MAX || !cuenta_eta_cruda_[tau - 1]) ? NaN : suma_eta_cruda_[tau - 1] / static_cast<double>(cuenta_eta_cruda_[tau - 1]);
    }

    // Ciclo caliente. Cero reservas de memoria.
    SalidaPaso paso(const SalidaDelta& s, double x, bool valido, double d_nu = NaN, std::size_t fase = 0, std::int64_t dia = 0,
                    double pared = NaN, const std::array<double, N_CAPA1>& capa1 = {}, std::uint8_t flags_hueco = 0) noexcept {
        SalidaPaso out{};
        if (!lista_) return out;
        ++cont_.pasos;
        if (s.t != t_ + 1) {
            ++cont_.pasos_desalineados;
            return out;
        }
        t_ = s.t;
        out.emitida = true;
        out.t = t_;
        if (std::isfinite(d_nu) && d_nu >= 0.0) { nu_acum_ += d_nu; hay_reloj_ = true; }
        if (std::isfinite(pared)) pared_ = pared;

        // 1. representaciones (espejo exacto de ObservadorSaltos::representar)
        const bool ok_x = valido && std::isfinite(x);
        double y_cruda = ok_x ? x : NaN;
        bool ok_cruda = ok_x;
        double y_rep = NaN;
        bool ok_rep = false;
        std::uint8_t flags_est = flags_hueco;
        if (integracion_ == Integracion::NIVEL) {
            ok_rep = ok_x && hay_prev_;
            y_rep = ok_rep ? x - x_prev_ : NaN;
            if (ok_rep) flags_est |= flags_prev_;
            hay_prev_ = ok_x;
            if (ok_x) { x_prev_ = x; flags_prev_ = flags_hueco; }
            else { flags_prev_ = 0; }
        } else {
            ok_rep = ok_x;
            y_rep = ok_rep ? x : NaN;
        }
        double y_purg = NaN;
        bool ok_purg = false;
        double yp_est = NaN;
        bool ok_yp_est = false;
        if (purga_activa_) {
            bool vp = false;
            const double yp = purga_.paso(ok_rep ? y_rep : NaN, fase, dia, vp);
            ok_yp_est = ok_rep && vp && std::isfinite(yp);
            yp_est = ok_yp_est ? yp : NaN;
            if (integracion_ == Integracion::NIVEL) {
                if (ok_rep && vp) {
                    if (!purgado_vivo_) { nivel_purgado_ = x - y_rep; purgado_vivo_ = true; }
                    nivel_purgado_ += yp;
                    y_purg = nivel_purgado_;
                    ok_purg = true;
                } else {
                    if (purgado_vivo_) ++seg_purgado_;
                    purgado_vivo_ = false;
                }
            } else {
                ok_purg = ok_rep && vp;
                y_purg = ok_purg ? yp : NaN;
            }
        }
        guardar(y_cruda, ok_cruda, y_purg, ok_purg, y_rep, ok_rep, yp_est, ok_yp_est, std::isfinite(d_nu) ? d_nu : 0.0, capa1, flags_est);

        // 2. apertura: la regla de ObservadorSaltos::paso (primera alarma de nivel, de fina a gruesa, con el evento de nivel cerrado)
        bool abre = false;
        EventoAlarma a_abre{};
        int j_abre = -1;
        if (!espejo_abierto_) {
            for (std::size_t j = 0; j < s.niveles && j < escalas_; ++j) {
                if (!s.nivel[j].alarma.hay) continue;
                abre = true; a_abre = s.nivel[j].alarma; j_abre = static_cast<int>(j);
                break;
            }
        }
        if (abre) { espejo_abierto_ = true; espejo_alarma_ = a_abre; espejo_escala_ = j_abre; }

        // 2b. Etapa 5 (E1.7): actualizar estado de calor multiescala desde las salidas del canal de temperatura
        gestor_calor_.paso(std::span<const ark::metrologia::SalidaTemperatura>(s.temp.data(), s.temp.size()), s.niveles, t_);

        // 3. avanzar las ranuras vivas en tau = t_ - r.t_apertura (si en este mismo paso abre una replica de signo opuesto,
        //    el cruce hacia el signo opuesto es la replica que corta la ranura, no una restitucion consolidada previa)
        // Fase 3: si la barra viene etiquetada con FLAG_PUENTE_ITO (reconstruida por RNG), Z3 NO la toma como cruce real medido
        const bool veto_z3_ito = ((flags_hueco & ark::dataset::FLAG_PUENTE_ITO) != 0u);
        const bool guarda_calor_prev = calor_activo_;
        const bool calor_detectado = gestor_calor_.calor_activo();
        for (auto& r : ranuras_) {
            if (!r.ocupada || t_ <= r.t_apertura) continue;
            const std::size_t tau = static_cast<std::size_t>(t_ - r.t_apertura);
            calor_activo_ = guarda_calor_prev || calor_detectado || (abre && a_abre.signo != r.signo);
            avanzar(r.purgada, tau, y_purg, !veto_z3_ito && ok_purg && hist_seg_actual_ok(r), true, r.w_ancla);
            avanzar(r.cruda, tau, y_cruda, !veto_z3_ito && ok_cruda, false, r.w_ancla);
            calor_activo_ = guarda_calor_prev;
            if (!r.alguna_rama_abierta()) cerrar(r, out);
        }

        // 4. apertura en Lilyth (E1.21): si hay alguna ranura viva, comparar contra la mas reciente (mayor id_evento)
        if (abre) {
            const Ranura* ref_viva = nullptr;
            for (const auto& r : ranuras_) {
                if (!r.ocupada) continue;
                if (ref_viva == nullptr || r.id_evento > ref_viva->id_evento) ref_viva = &r;
            }
            if (ref_viva != nullptr) {
                const bool mismo_signo = (a_abre.signo == ref_viva->signo);
                const bool sin_restituir = !ref_viva->ha_restituido();
                if (mismo_signo && sin_restituir) {
                    out.continuacion = true;
                    ++cont_.continuaciones;
                } else {
                    out.replica = true;
                    ++cont_.replicas;
                    if (cfg_.cortar_abiertas_con_replica) {
                        for (auto& r : ranuras_) {
                            if (!r.ocupada || t_ <= r.t_apertura) continue;
                            const std::size_t tau = static_cast<std::size_t>(t_ - r.t_apertura);
                            cortar_por_replica(r, tau, a_abre.signo);
                            if (!r.alguna_rama_abierta()) cerrar(r, out);
                        }
                    }
                    abrir(a_abre, j_abre, y_cruda, ok_cruda, y_purg, ok_purg, out);
                }
            } else {
                abrir(a_abre, j_abre, y_cruda, ok_cruda, y_purg, ok_purg, out);
            }
        }

        // 5. cierre del espejo: el operador delta emite el evento de nivel (canal 0) al clasificarlo.
        for (std::size_t i = 0; i < s.n_eventos; ++i) {
            const FichaEvento& e = s.eventos[i];
            if (e.canal != 0) continue;
            ++cont_.cierres_espejo;
            const bool casa = espejo_abierto_ && e.escala_origen == espejo_escala_ && e.alarma.t == espejo_alarma_.t &&
                              e.alarma.t_cambio == espejo_alarma_.t_cambio && e.alarma.signo == espejo_alarma_.signo &&
                              e.alarma.l == espejo_alarma_.l;
            if (!casa) ++cont_.espejo_inconsistente;
            espejo_abierto_ = false;
        }
        return out;
    }

    template <class SalidaFilas>
    SalidaPaso paso_filas(const SalidaFilas& sf, bool fase_de_barra, std::size_t fase = 0, std::int64_t dia = 0,
                          const std::array<double, N_CAPA1>& capa1 = {}) noexcept {
        if (!sf.barra_cerrada) return SalidaPaso{};
        const auto& b = sf.barra;
        const bool ok_barra = cfg_.usar_protocolo_filas_piloto
            ? ((!b.imputada || b.solo_mercado_quieto) && ((b.flags_hueco & (2u | 64u | 128u)) == 0u))
            : (!b.imputada);
        return paso(sf.delta, b.valor, ok_barra, b.nu, fase_de_barra ? b.fase : fase, fase_de_barra ? b.dia : dia,
                    static_cast<double>(b.t_cierre_ms), capa1, b.flags_hueco);
    }

    // Etapa 3 (seccion 3.6): evalua fuera del ciclo caliente todas las ranuras cerradas contra los nulos (a) y (b),
    // y calcula las pruebas del resumen (p_Delta_tau_iid, p_Delta_tau_Onsager) y la curva de Onsager P3b.
    void completar_nulos() {
        if (!calibrada_nulo_) return;
        if (nulo_cruda_.calibrada) {
            evaluar_nulos_representacion(nulo_cruda_, false);
            evaluar_armonico_representacion(nulo_cruda_, armonico_cruda_, false);
        }
        if (nulo_purgada_.calibrada) {
            evaluar_nulos_representacion(nulo_purgada_, true);
            evaluar_armonico_representacion(nulo_purgada_, armonico_purgada_, true);
        }
        (void)evaluar_marea();
    }

    // Etapa 6 (E1.14 [DECIDIDO]): evalua el resumen descriptivo de marea (Delta_tau_marea y Delta_RMST_marea con IC bootstrap de eventos,
    // fuera de la familia de 4 p de Holm).
    const ResumenMareaLilyth& evaluar_marea() {
        const std::size_t tm = lista_ ? cfg_.t_max_efectivo : T_MAX;
        resumen_marea_ = detalle_reloj_marea::evaluar_resumen_marea<Ranura>(
            cerradas_historico_, aj_cruda_.total.rmst, aj_purgada_.total.rmst,
            tm, cfg_.alfa, cfg_.b_mc_resumen, cfg_.semilla_nulo);
        return resumen_marea_;
    }

    // Expone la ejecucion de un lote de inyecciones en una misma realizacion para la verificacion directa de ventanas solapadas (S2).
    ResultadoInyeccionesRealizacion<T_MAX> evaluar_inyecciones_realizacion(
        bool rama_purgada,
        std::span<const double> hist_prev_est,
        std::span<const double> perfil_est,
        int signo_real,
        std::size_t w_ancla,
        std::size_t guarda,
        std::span<const std::int64_t> tiempos_iny,
        std::size_t n_sim,
        std::uint64_t semilla,
        bool sin_memoria = false,
        bool aplicar_guarda = true) {
        const auto& fn = rama_purgada ? nulo_purgada_ : nulo_cruda_;
        std::vector<double> e_iny(sin_memoria ? 0 : perfil_est.size(), 0.0);
        if (!sin_memoria) innovaciones(fn.modelo_lineal, hist_prev_est.data(), hist_prev_est.size(),
                                       perfil_est.data(), perfil_est.size(), e_iny.data());
        std::mt19937_64 rng(detalle_nulo::splitmix64(semilla));
        auto res = detalle_nulo::ejecutar_inyecciones_en_realizacion<T_MAX, Rasgos::b_max>(
            delta_calibrado_, *this, integracion_, rama_purgada, sin_memoria, cfg_.nulo_sin_memoria, fn.segmentos_calma,
            fn.modelo_lineal, ctx_purga_, x0_nivel_, perfil_est, e_iny,
            signo_real, w_ancla, guarda, cfg_.t_max_efectivo, tiempos_iny, n_sim, rng, aplicar_guarda);
        cont_.solapes += res.solapes;
        return res;
    }

private:
    // Code (revision D1, 2026-11-07, opcion (b) aprobada por el Comandante): para cada evento real y cada nulo,
    //   - via (a) si la excedencia (misma escala, mismo signo, g >= g_k) tiene al menos m_min eventos;
    //   - si no, via (b): en el nulo lineal, inyeccion de las innovaciones del perfil en el modelo; en el nulo SIN MEMORIA, el perfil real
    //     escrito en un sustituto sin memoria de la calma, con el mismo detector;
    //   - si (b) tampoco llega a m_min, el evento queda SIN nulo (no se rellena) y la ficha dice por que.
    // El pool del Monte Carlo usa la misma regla: el PIT de un evento nulo se calcula contra su excedencia sin el mismo, y solo si esta
    // tiene m_min eventos.
    void evaluar_una_referencia(FichaNuloRepresentacion<T_MAX>& fn, bool es_purgada, bool sin_memoria, Ranura& r, std::size_t idx,
                                ReferenciaNuloEvento& ref, std::vector<double>& u_pool_extra,
                                std::vector<EventoNuloLilyth<T_MAX>>& emparejados, std::size_t umbral) {
        // Code (2026-11-16): `umbral` = m_min con PoolNuloLilyth::M_MIN (hoy, bit a bit) o max(m_min, n*) con DERIVADO. Elige la via (a)
        // y es el objetivo de la via (b); el piso para evaluar el evento sigue siendo m_min.
        emparejados.clear();
        auto& rama = es_purgada ? r.purgada : r.cruda;
        const auto& pool = sin_memoria ? fn.nulos_sin_memoria : fn.nulos_lineal;
        auto esp = detalle_nulo::filtrar_espontaneos_excedencia<T_MAX>(r.escala_origen, r.signo, r.g_disparo, pool);
        ref.n_nulos_a = esp.size();
        std::vector<double> taus;
        auto usar_via_a = [&](const char* motivo) {
            ref.nulo_por_inyeccion = false;
            ref.motivo = motivo;
            for (const auto* p : esp) {
                taus.push_back(detalle_nulo::codificar_tau(p->tau_rest_consolidado, cfg_.t_max_efectivo));
                emparejados.push_back(*p);
            }
        };
        if (!cfg_.forzar_via_b && esp.size() >= umbral) {
            usar_via_a(sin_memoria ? "via (a) sin memoria: excedencia suficiente" : "via (a): excedencia espontanea suficiente");
        } else {
            const bool z2_corrupto = es_purgada ? r.perfil_corrupto_z2_purgada : r.perfil_corrupto_z2_cruda;
            if (z2_corrupto) {
                ref.sin_nulo_b = true;
                ref.evaluado = false;
                ref.motivo = "via (b): perfil Z2 corrupto";
                return;
            }
            const std::span<const double> hp = es_purgada
                ? std::span<const double>(r.hist_prev_purgada_est.data(), r.n_hist_prev_purgada)
                : std::span<const double>(r.hist_prev_cruda_est.data(), r.n_hist_prev_cruda);
            const std::span<const double> pf = es_purgada
                ? std::span<const double>(r.perfil_purgada_est.data(), r.n_perfil_purgada)
                : std::span<const double>(r.perfil_cruda_est.data(), r.n_perfil_cruda);
            std::vector<EventoNuloLilyth<T_MAX>> ev_b;
            const std::uint64_t sal = sin_memoria ? 0x51ED270B27C3A6E5ULL : 0x9E3779B97F4A7C15ULL;
            detalle_nulo::ejecutar_via_b_evento<T_MAX, Rasgos::b_max>(
                delta_calibrado_, *this, integracion_, es_purgada, sin_memoria, cfg_.nulo_sin_memoria, fn.segmentos_calma,
                fn.modelo_lineal, ctx_purga_, x0_nivel_, hp, pf, r.signo, r.w_ancla, r.guarda_ancla, umbral, cfg_.t_max_efectivo,
                cfg_.semilla_nulo ^ (sal * (idx + 1)), ref, ev_b, cont_.solapes);
            // Code (2026-11-16): con DERIVADO, si la via (b) no llega al objetivo y la excedencia de la via (a) ya llegaba a m_min y no es
            // menor, se queda la via (a): DERIVADO solo agrega nulos, nunca deja sin evaluar un evento que M_MIN evalua.
            if (!cfg_.forzar_via_b && ev_b.size() < umbral && esp.size() >= fn.m_min && esp.size() >= ev_b.size()) {
                usar_via_a(sin_memoria ? "via (a) sin memoria: la via (b) no llego al objetivo derivado"
                                       : "via (a): la via (b) no llego al objetivo derivado");
                ref.sin_nulo_b = false;
                detalle_nulo::llenar_estadisticas_referencia(ref, rama.tau_rest_consolidado, cfg_.t_max_efectivo, taus);
                return;
            }
            r.nulo_por_inyeccion = true;
            for (const auto& ev : ev_b) taus.push_back(detalle_nulo::codificar_tau(ev.tau_rest_consolidado, cfg_.t_max_efectivo));
            if (taus.size() < fn.m_min) {
                ref.sin_nulo_b = true;
                ref.evaluado = false;
                ref.motivo = sin_memoria ? "sin nulo: (a) y (b) sin memoria no llegan a m_min" : "sin nulo: (a) y (b) no llegan a m_min";
                return;
            }
            for (std::size_t i = 0; i < taus.size(); ++i) {
                std::vector<double> resto;
                for (std::size_t j = 0; j < taus.size(); ++j) if (j != i) resto.push_back(taus[j]);
                u_pool_extra.push_back(detalle_nulo::pit_punto_medio(taus[i], resto));
            }
            emparejados = std::move(ev_b);
        }
        detalle_nulo::llenar_estadisticas_referencia(ref, rama.tau_rest_consolidado, cfg_.t_max_efectivo, taus);
    }

    void evaluar_nulos_representacion(FichaNuloRepresentacion<T_MAX>& fn, bool es_purgada) {
        std::vector<double> u_iid, dt_iid, u_ons, dt_ons, u_pool_iid, u_pool_ons;
        std::vector<std::array<double, T_MAX>> etas_reales;
        std::vector<std::size_t> pasos_reales;
        detalle_nulo::pool_pit_loo<T_MAX>(fn.nulos_sin_memoria, fn.m_min, cfg_.t_max_efectivo, u_pool_iid);
        detalle_nulo::pool_pit_loo<T_MAX>(fn.nulos_lineal, fn.m_min, cfg_.t_max_efectivo, u_pool_ons);
        auto& emp_iid = es_purgada ? emp_iid_purgada_ : emp_iid_cruda_;
        auto& emp_ons = es_purgada ? emp_ons_purgada_ : emp_ons_cruda_;
        emp_iid.assign(cerradas_historico_.size(), {});
        emp_ons.assign(cerradas_historico_.size(), {});
        const bool inter = cfg_.contraste_resumen == ContrasteResumenLilyth::INTERCAMBIABLE;
        std::vector<std::vector<double>> est_iid, est_ons;

        // Code (2026-11-16): PoolNuloLilyth::DERIVADO. K = eventos candidatos de la rama (RESTITUIDA o CENSURADA_HORIZONTE), fijado antes
        // de evaluar (sin bucle); umbral por evento = max(m_min, n*). Con M_MIN el umbral es m_min (hoy).
        std::size_t umbral = fn.m_min;
        fn.pool_derivado = cfg_.pool_nulo == PoolNuloLilyth::DERIVADO;
        fn.k_candidatos_pool = 0;
        fn.objetivo_pool = 0;
        if (fn.pool_derivado) {
            for (const auto& rc : cerradas_historico_) {
                const auto& ra = es_purgada ? rc.purgada : rc.cruda;
                if (ra.desenlace == DesenlaceRanura::RESTITUIDA || ra.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) ++fn.k_candidatos_pool;
            }
            fn.objetivo_pool = objetivo_pool_derivado(fn.k_candidatos_pool, cfg_.alfa);
            umbral = std::max(fn.m_min, fn.objetivo_pool);
        }

        for (std::size_t idx = 0; idx < cerradas_historico_.size(); ++idx) {
            Ranura& r = cerradas_historico_[idx];
            auto& rama = es_purgada ? r.purgada : r.cruda;
            if (rama.desenlace != DesenlaceRanura::RESTITUIDA && rama.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
            etas_reales.push_back(rama.eta);
            pasos_reales.push_back(rama.pasos_observados);

            ReferenciaNuloEvento ref_iid{};
            evaluar_una_referencia(fn, es_purgada, true, r, idx, ref_iid, u_pool_iid, emp_iid[idx], umbral);
            rama.nulo_iid = ref_iid;
            if (ref_iid.evaluado) { u_iid.push_back(ref_iid.u_pit); dt_iid.push_back(ref_iid.delta_tau); }

            ReferenciaNuloEvento ref_ons{};
            if (fn.modelo_lineal.valido) evaluar_una_referencia(fn, es_purgada, false, r, idx, ref_ons, u_pool_ons, emp_ons[idx], umbral);
            else ref_ons.motivo = "sin nulo lineal: modelo invalido";
            rama.nulo_onsager = ref_ons;
            if (ref_ons.evaluado) { u_ons.push_back(ref_ons.u_pit); dt_ons.push_back(ref_ons.delta_tau); }

            // Code (2026-11-15): estratos de la rama INTERCAMBIABLE = {tau del real} U {tau de sus nulos emparejados}, en el codigo del PIT.
            if (inter) {
                auto estrato = [&](const std::vector<EventoNuloLilyth<T_MAX>>& emp) {
                    std::vector<double> e{detalle_nulo::codificar_tau(rama.tau_rest_consolidado, cfg_.t_max_efectivo)};
                    for (const auto& ev : emp) e.push_back(detalle_nulo::codificar_tau(ev.tau_rest_consolidado, cfg_.t_max_efectivo));
                    return e;
                };
                if (ref_iid.evaluado) est_iid.push_back(estrato(emp_iid[idx]));
                if (ref_ons.evaluado) est_ons.push_back(estrato(emp_ons[idx]));
            }
        }

        fn.contraste_iid = detalle_nulo::evaluar_contraste_resumen(
            u_iid, dt_iid, u_pool_iid, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ 0xABCDEF0123456789ULL);
        fn.contraste_onsager = detalle_nulo::evaluar_contraste_resumen(
            u_ons, dt_ons, u_pool_ons, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ 0x9876543210FEDCBAULL);
        if (inter) {
            fn.contraste_iid_pool_loo = fn.contraste_iid;
            fn.contraste_onsager_pool_loo = fn.contraste_onsager;
            fn.contraste_iid = detalle_nulo::evaluar_contraste_intercambiable(
                est_iid, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ 0xABCDEF0123456789ULL, cfg_.alfa / 4.0);
            fn.contraste_onsager = detalle_nulo::evaluar_contraste_intercambiable(
                est_ons, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ 0x9876543210FEDCBAULL, cfg_.alfa / 4.0);
        }
        fn.p3b = detalle_nulo::evaluar_curva_p3b<T_MAX>(
            integracion_, fn.segmentos_calma, fn.modelo_lineal, etas_reales, pasos_reales,
            // Code (revision D6): hasta el K de la calma (los rezagos con ACF fuera de su banda), no hasta 10 escrito.
            std::min<std::size_t>(std::max<std::size_t>(1, fn.modelo_lineal.K), cfg_.t_max_efectivo), cfg_.semilla_nulo, cfg_.b_mc_resumen, cfg_.alfa);
    }

    void evaluar_armonico_representacion(const FichaNuloRepresentacion<T_MAX>& fn, ResumenArmonicoYHolm<T_MAX>& arm, bool es_purgada) {
        arm = ResumenArmonicoYHolm<T_MAX>{};
        const std::size_t H = std::min<std::size_t>(cfg_.t_max_efectivo, T_MAX);
        if (H < 3) return;

        // 1. Curva media eta_barra(1..H) y Prony de orden 2 (PRONY_CURVA_MEDIA, E1.8)
        std::array<double, T_MAX> eta_barra{};
        eta_barra.fill(NaN);
        for (std::size_t tau = 1; tau <= H; ++tau) {
            eta_barra[tau - 1] = es_purgada ? eta_barra_purgada(tau) : eta_barra_cruda(tau);
        }
        arm.prony_curva_media = detalle_armonico::ajustar_prony(std::span<const double>(eta_barra.data(), H), H, 1.0);

        // 2. Recoger todas las ranuras medidas para la fase reescalada (E1.20), sobrepaso por evento y etiqueta ARMONICO_SOSPECHOSO
        std::vector<std::array<double, T_MAX>> etas_medidas;
        std::vector<std::size_t> pasos_medidos;
        std::vector<double> taus_rest_medidos;
        std::vector<std::size_t> idx_medidos;

        for (std::size_t idx = 0; idx < cerradas_historico_.size(); ++idx) {
            const auto& r = cerradas_historico_[idx];
            const auto& rama = es_purgada ? r.purgada : r.cruda;
            if (rama.desenlace == DesenlaceRanura::INACTIVA || rama.desenlace == DesenlaceRanura::ABIERTA || es_no_medida(rama.desenlace)) {
                continue;
            }
            etas_medidas.push_back(rama.eta);
            pasos_medidos.push_back(rama.pasos_observados);
            taus_rest_medidos.push_back(rama.tau_rest_consolidado);
            idx_medidos.push_back(idx);
        }

        arm.fase_reescalada = detalle_armonico::construir_fase_reescalada<T_MAX>(
            etas_medidas, pasos_medidos, taus_rest_medidos, H);

        // 3. Code (revision D1c): las curvas nulas son los nulos EMPAREJADOS de cada evento real en la Etapa 3 (excedencia en (a), o su
        //    via (b)), no un pool filtrado por |dY| >= 0.5 mediana (numero escrito) y completado con eventos fabricados sin detector. El
        //    min de eta_barra real se calcula sobre los MISMOS eventos que tienen nulo (RESTITUIDA o CENSURADA_HORIZONTE), para que las dos
        //    curvas promedien lo mismo.
        const auto& emp_iid = es_purgada ? emp_iid_purgada_ : emp_iid_cruda_;
        const auto& emp_ons = es_purgada ? emp_ons_purgada_ : emp_ons_cruda_;
        auto sets_de = [&](const std::vector<std::vector<EventoNuloLilyth<T_MAX>>>& emp,
                           std::vector<std::size_t>& idx_con_nulo) {
            std::vector<const std::vector<EventoNuloLilyth<T_MAX>>*> sets;
            idx_con_nulo.clear();
            for (std::size_t idx = 0; idx < emp.size() && idx < cerradas_historico_.size(); ++idx) {
                if (emp[idx].empty()) continue;
                sets.push_back(&emp[idx]);
                idx_con_nulo.push_back(idx);
            }
            return sets;
        };
        std::vector<std::size_t> idx_iid, idx_ons;
        const auto sets_iid = sets_de(emp_iid, idx_iid);
        const auto sets_ons = sets_de(emp_ons, idx_ons);
        auto min_eta_barra_de = [&](const std::vector<std::size_t>& idxs) {
            std::vector<double> suma(H, 0.0);
            std::vector<std::size_t> cuenta(H, 0);
            for (std::size_t idx : idxs) {
                const auto& rama = es_purgada ? cerradas_historico_[idx].purgada : cerradas_historico_[idx].cruda;
                const std::size_t pk = std::min<std::size_t>(rama.pasos_observados, H);
                for (std::size_t tau = 2; tau <= pk; ++tau) {
                    if (std::isfinite(rama.eta[tau - 1])) { suma[tau - 1] += rama.eta[tau - 1]; ++cuenta[tau - 1]; }
                }
            }
            double mn = NaN;
            for (std::size_t tau = 2; tau <= H; ++tau) {
                if (cuenta[tau - 1] == 0) continue;
                const double m = suma[tau - 1] / static_cast<double>(cuenta[tau - 1]);
                if (!std::isfinite(mn) || m < mn) mn = m;
            }
            return mn;
        };

        // 4. Contraste Monte Carlo del sobrepaso de la curva media contra ambos nulos (E1.8), emparejado evento a evento.
        arm.sobrepaso.min_eta_barra_real = arm.prony_curva_media.min_eta_desde_2;
        arm.sobrepaso.min_eta_barra_real_emparejado_iid = min_eta_barra_de(idx_iid);
        arm.sobrepaso.min_eta_barra_real_emparejado_onsager = min_eta_barra_de(idx_ons);
        detalle_armonico::contraste_min_eta_barra_emparejado<T_MAX>(
            arm.sobrepaso.min_eta_barra_real_emparejado_iid, sets_iid, H, cfg_.b_mc_resumen,
            cfg_.semilla_nulo ^ 0xA5A5A5A55A5A5A5AULL,
            arm.sobrepaso.mediana_min_eta_barra_iid, arm.sobrepaso.p_sobrepaso_curva_iid);
        detalle_armonico::contraste_min_eta_barra_emparejado<T_MAX>(
            arm.sobrepaso.min_eta_barra_real_emparejado_onsager, sets_ons, H, cfg_.b_mc_resumen,
            cfg_.semilla_nulo ^ 0x5A5A5A5AA5A5A5A5ULL,
            arm.sobrepaso.mediana_min_eta_barra_onsager, arm.sobrepaso.p_sobrepaso_curva_onsager);

        // 5. SOBREPASO_POR_EVENTO: s_k = min de eta_k en [tau_rest_k, T_max], PIT de punto medio contra SU conjunto emparejado; el pool
        //    del Monte Carlo es el PIT de cada nulo contra el resto de su conjunto.
        auto sob_de = [&](const EventoNuloLilyth<T_MAX>& ev) {
            return detalle_armonico::min_eta_post_restitucion<T_MAX>(ev.eta, ev.pasos_observados, ev.tau_rest_consolidado, H);
        };
        auto contraste_sob = [&](const std::vector<const std::vector<EventoNuloLilyth<T_MAX>>*>& sets,
                                 const std::vector<std::size_t>& idxs, std::uint64_t sal,
                                 std::vector<std::vector<double>>& estratos) {   // Code (2026-11-15): {s_k} U sn, para INTERCAMBIABLE
            estratos.clear();
            std::vector<double> u_k, d_k, u_pool;
            for (std::size_t q = 0; q < sets.size(); ++q) {
                std::vector<double> sn;
                for (const auto& ev : *sets[q]) { const double v = sob_de(ev); if (std::isfinite(v)) sn.push_back(v); }
                if (sn.size() < 2) continue;
                for (std::size_t i = 0; i < sn.size(); ++i) {
                    std::vector<double> resto;
                    for (std::size_t j = 0; j < sn.size(); ++j) if (j != i) resto.push_back(sn[j]);
                    u_pool.push_back(detalle_nulo::pit_punto_medio(sn[i], resto));
                }
                const auto& rama = es_purgada ? cerradas_historico_[idxs[q]].purgada : cerradas_historico_[idxs[q]].cruda;
                const double s_k = detalle_armonico::min_eta_post_restitucion<T_MAX>(rama.eta, rama.pasos_observados, rama.tau_rest_consolidado, H);
                if (!std::isfinite(s_k)) continue;
                u_k.push_back(detalle_nulo::pit_punto_medio(s_k, sn));
                d_k.push_back(s_k - ark::metrologia::est::cuantil(sn, 0.5));
                if (cfg_.contraste_resumen == ContrasteResumenLilyth::INTERCAMBIABLE) {
                    std::vector<double> e{s_k};
                    e.insert(e.end(), sn.begin(), sn.end());
                    estratos.push_back(std::move(e));
                }
            }
            return detalle_nulo::evaluar_contraste_resumen(u_k, d_k, u_pool, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ sal);
        };
        std::vector<std::vector<double>> est_sob_iid, est_sob_ons;
        arm.sobrepaso.sobrepaso_por_evento_iid = contraste_sob(sets_iid, idx_iid, 0x3333333333333333ULL, est_sob_iid);
        arm.sobrepaso.sobrepaso_por_evento_onsager = contraste_sob(sets_ons, idx_ons, 0x7777777777777777ULL, est_sob_ons);
        if (cfg_.contraste_resumen == ContrasteResumenLilyth::INTERCAMBIABLE) {
            arm.sobrepaso.sobrepaso_por_evento_iid_pool_loo = arm.sobrepaso.sobrepaso_por_evento_iid;
            arm.sobrepaso.sobrepaso_por_evento_onsager_pool_loo = arm.sobrepaso.sobrepaso_por_evento_onsager;
            arm.sobrepaso.sobrepaso_por_evento_iid = detalle_nulo::evaluar_contraste_intercambiable(
                est_sob_iid, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ 0x3333333333333333ULL, cfg_.alfa / 4.0);
            arm.sobrepaso.sobrepaso_por_evento_onsager = detalle_nulo::evaluar_contraste_intercambiable(
                est_sob_ons, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ 0x7777777777777777ULL, cfg_.alfa / 4.0);
        }

        // 6. Etiqueta ARMONICO_SOSPECHOSO por evento (Prony con >= 4 puntos y raices complejas Y min_eta_desde_2 < Q_alfa(min_eta_nulo)),
        //    con el cuantil sobre los nulos emparejados (cada nulo pesa por el evento real al que esta emparejado).
        auto juntar = [&](const std::vector<const std::vector<EventoNuloLilyth<T_MAX>>*>& sets) {
            std::vector<EventoNuloLilyth<T_MAX>> out;
            for (const auto* s : sets) out.insert(out.end(), s->begin(), s->end());
            return out;
        };
        const std::vector<EventoNuloLilyth<T_MAX>> pool_iid = juntar(sets_iid);
        const std::vector<EventoNuloLilyth<T_MAX>> pool_ons = juntar(sets_ons);
        auto extraer_mins_desde_2 = [&](const std::vector<EventoNuloLilyth<T_MAX>>& pool) {
            std::vector<double> mins;
            mins.reserve(pool.size());
            for (const auto& ev : pool) {
                const double mn = detalle_armonico::min_eta_en_rango<T_MAX>(ev.eta, ev.pasos_observados, 2, H);
                if (std::isfinite(mn)) mins.push_back(mn);
            }
            return mins;
        };
        auto cuantil_orden_inferior = [&](std::vector<double> v, double alfa_q) -> double {
            if (v.empty()) return NaN;
            std::sort(v.begin(), v.end());
            const std::size_t idx_q = std::min(v.size() - 1, static_cast<std::size_t>(std::floor(alfa_q * static_cast<double>(v.size()))));
            return v[idx_q];
        };
        const std::vector<double> mins_iid = extraer_mins_desde_2(pool_iid);
        const std::vector<double> mins_ons = extraer_mins_desde_2(pool_ons);
        arm.sobrepaso.cuantil_alfa_min_eta_iid = cuantil_orden_inferior(mins_iid, cfg_.alfa);
        arm.sobrepaso.cuantil_alfa_min_eta_onsager = cuantil_orden_inferior(mins_ons, cfg_.alfa);

        auto calcular_tasa_etiqueta_pool = [&](const std::vector<EventoNuloLilyth<T_MAX>>& pool, double q_alfa) {
            std::size_t n_eval = 0, n_etiq = 0;
            for (const auto& ev : pool) {
                const std::size_t pk = std::min({ev.pasos_observados, H, T_MAX});
                if (pk < 4) continue;
                const double mn = detalle_armonico::min_eta_en_rango<T_MAX>(ev.eta, pk, 2, H);
                if (!std::isfinite(mn)) continue;
                ++n_eval;
                if (detalle_armonico::evaluar_etiqueta_armonico_evento<T_MAX>(ev.eta, pk, mn, H, q_alfa)) {
                    ++n_etiq;
                }
            }
            return n_eval > 0 ? static_cast<double>(n_etiq) / static_cast<double>(n_eval) : 0.0;
        };
        arm.sobrepaso.tasa_etiqueta_nulo_iid = calcular_tasa_etiqueta_pool(pool_iid, arm.sobrepaso.cuantil_alfa_min_eta_iid);
        arm.sobrepaso.tasa_etiqueta_nulo_onsager = calcular_tasa_etiqueta_pool(pool_ons, arm.sobrepaso.cuantil_alfa_min_eta_onsager);

        const double q_alfa_real = std::isfinite(arm.sobrepaso.cuantil_alfa_min_eta_onsager)
                                       ? arm.sobrepaso.cuantil_alfa_min_eta_onsager
                                       : arm.sobrepaso.cuantil_alfa_min_eta_iid;
        std::size_t n_reales_eval = 0, n_reales_etiq = 0;
        for (std::size_t idx : idx_medidos) {
            auto& r = cerradas_historico_[idx];
            auto& rama = es_purgada ? r.purgada : r.cruda;
            const double mn = detalle_armonico::min_eta_en_rango<T_MAX>(rama.eta, rama.pasos_observados, 2, H);
            rama.min_eta_desde_2 = mn;
            const bool etiq = detalle_armonico::evaluar_etiqueta_armonico_evento<T_MAX>(
                rama.eta, rama.pasos_observados, mn, H, q_alfa_real);
            rama.armonico_sospechoso = etiq;
            if (std::min({rama.pasos_observados, H, T_MAX}) >= 4 && std::isfinite(mn)) {
                ++n_reales_eval;
                if (etiq) ++n_reales_etiq;
            }
            for (auto& arch : archivo_) {
                if (arch.id_evento == r.id_evento && arch.id_evento != 0) {
                    if (es_purgada) arch.purgada.armonico_sospechoso = etiq;
                    else arch.cruda.armonico_sospechoso = etiq;
                }
            }
        }
        arm.sobrepaso.eventos_etiquetados_reales = n_reales_etiq;
        arm.sobrepaso.tasa_etiqueta_reales = n_reales_eval > 0
                                                 ? static_cast<double>(n_reales_etiq) / static_cast<double>(n_reales_eval)
                                                 : 0.0;

        // 7. Familia de 4 hipotesis contrastada por ark::holm (E1.16)
        const double p_dt_iid = fn.contraste_iid.p_bilateral;
        const double p_dt_ons = fn.contraste_onsager.p_bilateral;
        const double p_sob_iid = (cfg_.modo_armonico == ModoArmonico::SOBREPASO_POR_EVENTO)
                                     ? arm.sobrepaso.sobrepaso_por_evento_iid.p_bilateral
                                     : arm.sobrepaso.p_sobrepaso_curva_iid;
        const double p_sob_ons = (cfg_.modo_armonico == ModoArmonico::SOBREPASO_POR_EVENTO)
                                     ? arm.sobrepaso.sobrepaso_por_evento_onsager.p_bilateral
                                     : arm.sobrepaso.p_sobrepaso_curva_onsager;
        detalle_armonico::aplicar_familia_holm_cuatro<T_MAX>(p_dt_iid, p_dt_ons, p_sob_iid, p_sob_ons, cfg_.alfa, arm);
        arm.sobrepaso.valido = std::isfinite(p_sob_iid) && std::isfinite(p_sob_ons);
        arm.valido = arm.prony_curva_media.valido && arm.familia_holm_valida;
    }


    [[nodiscard]] bool hist_seg_actual_ok(const Ranura& r) const noexcept {
        return r.purgada.desenlace != DesenlaceRanura::ABIERTA || segmento_de_rama(r) == seg_purgado_;
    }
    [[nodiscard]] std::uint32_t segmento_de_rama(const Ranura& r) const noexcept {
        std::size_t idx = 0;
        return indice(r.t_apertura, idx) ? hist_seg_purgada_[idx] : ~0u;
    }

    void guardar(double yc, bool okc, double yp, bool okp, double y_est_c, bool ok_est_c, double y_est_p, bool ok_est_p,
                 double dnu, const std::array<double, N_CAPA1>& c1, std::uint8_t flags_est = 0) noexcept {
        const std::size_t pos = static_cast<std::size_t>(t_) % HIST_MAX;
        hist_cruda_[pos] = yc; hist_val_cruda_[pos] = okc ? 1u : 0u;
        hist_purgada_[pos] = yp; hist_val_purgada_[pos] = okp ? 1u : 0u;
        hist_est_cruda_[pos] = y_est_c; hist_val_est_cruda_[pos] = ok_est_c ? 1u : 0u;
        hist_est_purgada_[pos] = y_est_p; hist_val_est_purgada_[pos] = ok_est_p ? 1u : 0u;
        hist_flags_hueco_[pos] = flags_est;
        hist_seg_purgada_[pos] = seg_purgado_;
        hist_dnu_[pos] = dnu;
        hist_capa1_[pos] = c1;
        if (n_hist_ < HIST_MAX) ++n_hist_;
    }
    [[nodiscard]] bool indice(std::int64_t t, std::size_t& idx) const noexcept {
        if (t < 0 || t > t_ || static_cast<std::size_t>(t_ - t) >= n_hist_) return false;
        idx = static_cast<std::size_t>(t) % HIST_MAX;
        return true;
    }

    DesenlaceRanura mediana(std::int64_t t_ini, std::int64_t t_fin, bool purgada, std::uint32_t seg, double& med, std::size_t& n,
                            double& nu_punto) noexcept {
        med = NaN; n = 0; nu_punto = NaN;
        if (t_fin < t_ini || t_ini < 0) return DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE;
        std::size_t i0 = 0;
        if (!indice(t_ini, i0)) return DesenlaceRanura::NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA;
        double snu = 0.0;
        for (std::int64_t t = t_ini; t <= t_fin; ++t) {
            std::size_t idx = 0;
            (void)indice(t, idx);
            snu += hist_dnu_[idx];
            if (purgada && hist_seg_purgada_[idx] != seg) return DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE;
            const bool v = purgada ? hist_val_purgada_[idx] != 0u : hist_val_cruda_[idx] != 0u;
            const double y = purgada ? hist_purgada_[idx] : hist_cruda_[idx];
            if (v && std::isfinite(y)) scratch_[n++] = y;
        }
        if (hay_reloj_) nu_punto = snu / static_cast<double>(t_fin - t_ini + 1);
        if (n == 0) return DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE;
        const std::size_t mid = n / 2;
        std::nth_element(scratch_.begin(), scratch_.begin() + static_cast<std::ptrdiff_t>(mid), scratch_.begin() + static_cast<std::ptrdiff_t>(n));
        const double hi = scratch_[mid];
        if (n % 2u) { med = hi; return DesenlaceRanura::ABIERTA; }
        const double lo = *std::max_element(scratch_.begin(), scratch_.begin() + static_cast<std::ptrdiff_t>(mid));
        med = 0.5 * (lo + hi);
        return DesenlaceRanura::ABIERTA;
    }

    void preparar_rama(RamaTrayectoriaLilyth<T_MAX>& rama, bool purgada, bool habilitada, std::int64_t t_ini, std::int64_t t_fin,
                       double y_a, bool ok_a, int signo) noexcept {
        rama.limpiar();
        if (!habilitada) return;
        if (purgada && (!purga_activa_ || !ok_a)) {
            rama.desenlace = DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA;
            ++cont_.no_medida[static_cast<int>(rama.desenlace)];
            return;
        }
        if (!ok_a || !std::isfinite(y_a)) { rama.desenlace = DesenlaceRanura::CENSURADA_DATO_INVALIDO; return; }
        double med = NaN, nu_p = NaN;
        std::size_t n = 0;
        const DesenlaceRanura d = mediana(t_ini, t_fin, purgada, seg_purgado_, med, n, nu_p);
        if (d != DesenlaceRanura::ABIERTA) { rama.desenlace = d; ++cont_.no_medida[static_cast<int>(d)]; return; }
        const double dy = y_a - med;
        if (!(std::isfinite(dy) && dy * static_cast<double>(signo) > 0.0)) {
            rama.desenlace = DesenlaceRanura::NO_MEDIDA_SIGNO_INCONSISTENTE;
            ++cont_.no_medida[static_cast<int>(rama.desenlace)];
            return;
        }
        rama.desenlace = DesenlaceRanura::ABIERTA;
        rama.y_ancla = med; rama.y_disparo = y_a; rama.dy = dy;
        rama.nu_punto_ancla = nu_p; rama.nu_en_disparo = hay_reloj_ ? nu_acum_ : NaN;
        rama.pared_en_disparo = pared_;
        rama.muestras_ancla = n;
    }

    // Prediccion causal a 1 paso sin reservas dinamicas (pila) usando el ModeloLinealNulo calibrado
    double predecir_un_paso_causal(const ModeloLinealNulo& m, const double* buf, std::size_t n_buf) const noexcept {
        if (!m.valido || n_buf == 0) return m.valido && std::isfinite(m.mu) ? m.mu : 0.0;
        if (m.modo_usado != ModoNuloLineal::ARFIMA_WHITTLE || m.d == 0.0) {
            double pred = m.mu;
            for (std::size_t i = 1; i <= m.p && i <= n_buf; ++i) {
                pred += m.phi[i - 1] * (buf[n_buf - i] - m.mu);
            }
            return pred;
        }
        std::array<double, MAX_PERFIL_EVENTO_LILYTH * 2> pi_loc{};
        std::array<double, MAX_PERFIL_EVENTO_LILYTH * 2> u_loc{};
        const std::size_t n_ef = std::min(n_buf + 1, pi_loc.size());
        ark::metrologia::elw::coeficientes(m.d, n_ef, pi_loc.data());
        for (std::size_t t = 0; t < n_buf && t < u_loc.size(); ++t) {
            double s = 0.0;
            for (std::size_t k = 0; k <= t && k < n_ef; ++k) s += pi_loc[k] * buf[t - k];
            u_loc[t] = s;
        }
        double u_pred = m.mu;
        for (std::size_t i = 1; i <= m.p && i <= n_buf; ++i) {
            u_pred += m.phi[i - 1] * (u_loc[n_buf - i] - m.mu);
        }
        double cola = 0.0;
        for (std::size_t k = 1; k <= n_buf && k < n_ef; ++k) {
            cola += pi_loc[k] * buf[n_buf - k];
        }
        return u_pred - cola;
    }

    void copiar_perfil_estacionario(Ranura& r) noexcept {
        const std::int64_t tc = std::max<std::int64_t>(0, std::min(r.t_cambio, r.t_apertura));
        const std::size_t l_raw = static_cast<std::size_t>(r.t_apertura - tc + 1);
        const std::size_t l = std::min<std::size_t>(l_raw, MAX_PERFIL_EVENTO_LILYTH);
        const std::int64_t t_ini_perf = r.t_apertura - static_cast<std::int64_t>(l) + 1;
        const std::int64_t t_prev_ini = std::max<std::int64_t>(0, t_ini_perf - static_cast<std::int64_t>(MAX_PERFIL_EVENTO_LILYTH));

        // Camino por omision (COMPRESION_BASE en Z1 y Z2): identico linea por linea al original
        if (cfg_.politica_z1_piloto == ConfiguracionLilyth::PoliticaZ1Piloto::COMPRESION_BASE &&
            cfg_.politica_z2_piloto == ConfiguracionLilyth::PoliticaZ2Piloto::COMPRESION_BASE) {
            r.n_perfil_cruda = r.n_perfil_purgada = 0;
            for (std::int64_t t = t_ini_perf; t <= r.t_apertura; ++t) {
                std::size_t idx = 0;
                if (!indice(t, idx)) continue;
                if (hist_val_est_cruda_[idx] != 0u && std::isfinite(hist_est_cruda_[idx])) {
                    r.perfil_cruda_est[r.n_perfil_cruda++] = hist_est_cruda_[idx];
                }
                if (hist_val_est_purgada_[idx] != 0u && std::isfinite(hist_est_purgada_[idx])) {
                    r.perfil_purgada_est[r.n_perfil_purgada++] = hist_est_purgada_[idx];
                }
            }

            r.n_hist_prev_cruda = r.n_hist_prev_purgada = 0;
            for (std::int64_t t = t_prev_ini; t < t_ini_perf; ++t) {
                std::size_t idx = 0;
                if (!indice(t, idx)) continue;
                if (hist_val_est_cruda_[idx] != 0u && std::isfinite(hist_est_cruda_[idx]) && r.n_hist_prev_cruda < MAX_PERFIL_EVENTO_LILYTH) {
                    r.hist_prev_cruda_est[r.n_hist_prev_cruda++] = hist_est_cruda_[idx];
                }
                if (hist_val_est_purgada_[idx] != 0u && std::isfinite(hist_est_purgada_[idx]) && r.n_hist_prev_purgada < MAX_PERFIL_EVENTO_LILYTH) {
                    r.hist_prev_purgada_est[r.n_hist_prev_purgada++] = hist_est_purgada_[idx];
                }
            }
            return;
        }

        // Politica instrumentada para una representacion (cruda o purgada): primero Z1 (hist_prev), luego Z2 (perfil)
        auto copiar_rep = [&](const std::array<double, HIST_MAX>& h_est,
                              const std::array<std::uint8_t, HIST_MAX>& h_val,
                              const ModeloLinealNulo& mod,
                              std::array<double, MAX_PERFIL_EVENTO_LILYTH>& out_prev,
                              std::size_t& n_prev,
                              std::array<double, MAX_PERFIL_EVENTO_LILYTH>& out_perf,
                              std::size_t& n_perf,
                              bool& flag_corrupto_z2) noexcept {
            n_prev = 0;
            n_perf = 0;
            flag_corrupto_z2 = false;

            // 1. Z1: hist_prev_est en [t_prev_ini, t_ini_perf)
            if (cfg_.politica_z1_piloto == ConfiguracionLilyth::PoliticaZ1Piloto::COMPRESION_BASE) {
                for (std::int64_t t = t_prev_ini; t < t_ini_perf; ++t) {
                    std::size_t idx = 0;
                    if (!indice(t, idx)) continue;
                    if (h_val[idx] != 0u && std::isfinite(h_est[idx]) && n_prev < MAX_PERFIL_EVENTO_LILYTH) {
                        out_prev[n_prev++] = h_est[idx];
                    }
                }
            } else if (cfg_.politica_z1_piloto == ConfiguracionLilyth::PoliticaZ1Piloto::CORTE_SUFIJO) {
                std::int64_t t_corte = t_prev_ini;
                for (std::int64_t t = t_prev_ini; t < t_ini_perf; ++t) {
                    std::size_t idx = 0;
                    if (!indice(t, idx) || h_val[idx] == 0u || !std::isfinite(h_est[idx])) {
                        t_corte = t + 1;
                    }
                }
                for (std::int64_t t = t_corte; t < t_ini_perf; ++t) {
                    std::size_t idx = 0;
                    if (indice(t, idx) && h_val[idx] != 0u && std::isfinite(h_est[idx]) && n_prev < MAX_PERFIL_EVENTO_LILYTH) {
                        out_prev[n_prev++] = h_est[idx];
                    }
                }
            } else if (cfg_.politica_z1_piloto == ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_CAUSAL) {
                // PUENTE_CAUSAL en Z1 hasta r_max_z1_piloto; si excede r_max_z1_piloto, corta al sufijo posterior
                std::size_t racha_inv = 0;
                for (std::int64_t t = t_prev_ini; t < t_ini_perf; ++t) {
                    std::size_t idx = 0;
                    const bool ok_t = indice(t, idx) && h_val[idx] != 0u && std::isfinite(h_est[idx]);
                    if (ok_t) {
                        racha_inv = 0;
                        if (n_prev < MAX_PERFIL_EVENTO_LILYTH) out_prev[n_prev++] = h_est[idx];
                    } else {
                        ++racha_inv;
                        if (racha_inv <= cfg_.r_max_z1_piloto && n_prev > 0 && n_prev < MAX_PERFIL_EVENTO_LILYTH) {
                            const double zp = predecir_un_paso_causal(mod, out_prev.data(), n_prev);
                            out_prev[n_prev++] = zp;
                        } else if (racha_inv > cfg_.r_max_z1_piloto) {
                            n_prev = 0;
                        }
                    }
                }
            } else {
                // PUENTE_ITO_HISTORICO en Z1: usa puente condicionado de 2 bordes con doble llave (r <= r_max && |Z_hueco| <= z_max)
                // cuando ambos bordes del hueco están dentro de Z1 (< t_ini_perf), y predicción causal en el borde t_ini_perf - 1
                std::int64_t t = t_prev_ini;
                while (t < t_ini_perf) {
                    std::size_t idx = 0;
                    const bool ok_t = indice(t, idx) && h_val[idx] != 0u && std::isfinite(h_est[idx]);
                    if (ok_t) {
                        if (n_prev < MAX_PERFIL_EVENTO_LILYTH) out_prev[n_prev++] = h_est[idx];
                        ++t;
                        continue;
                    }
                    // Medir longitud de la racha inválida dentro de Z1 [t, t_ini_perf)
                    std::int64_t t_fin_gap = t;
                    while (t_fin_gap + 1 < t_ini_perf) {
                        std::size_t idx_sig = 0;
                        if (indice(t_fin_gap + 1, idx_sig) && h_val[idx_sig] != 0u && std::isfinite(h_est[idx_sig])) break;
                        ++t_fin_gap;
                    }
                    const std::size_t r_gap = static_cast<std::size_t>(t_fin_gap - t + 1);
                    if (n_prev == 0 || r_gap > cfg_.r_max_z1_piloto) {
                        n_prev = 0;
                        t = t_fin_gap + 1;
                        continue;
                    }
                    if (t_fin_gap + 1 < t_ini_perf) {
                        // Hueco interior Z1b: ambos bordes pertenecen a la historia previa en calma (< t_ini_perf)
                        std::size_t idx_der = 0;
                        (void)indice(t_fin_gap + 1, idx_der);
                        const double x_ini = out_prev[n_prev - 1];
                        const double x_fin = h_est[idx_der];
                        piloto::ConfiguracionPuenteHistoricoIto cfg_ito{};
                        cfg_ito.integracion = Integracion::FLUJO; // h_est ya está en representación estacionaria
                        cfg_ito.r_max = cfg_.r_max_z1_piloto;
                        cfg_ito.z_max_hueco = cfg_.z_max_hueco_piloto;
                        cfg_ito.n_semillas_rng = 0; // En copiar_perfil_estacionario usamos la esperanza condicional exacta
                        const auto res_ito = piloto::evaluar_y_reconstruir_hueco_ito(
                            cfg_ito, mod, std::span<const double>(out_prev.data(), n_prev), x_ini, x_fin, r_gap);
                        if (res_ito.reconstruido) {
                            for (std::size_t k = 0; k < r_gap && n_prev < MAX_PERFIL_EVENTO_LILYTH; ++k) {
                                out_prev[n_prev++] = res_ito.x_media_cond[k];
                            }
                        } else {
                            n_prev = 0; // Falló la llave de amplitud |Z_hueco| <= z_max: corte limpio al sufijo
                        }
                    } else {
                        // Hueco en el borde inmediato Z1a (t_fin_gap == t_ini_perf - 1): no usa t_ini_perf para no filtrar el choque Z2
                        for (std::size_t k = 0; k < r_gap && n_prev < MAX_PERFIL_EVENTO_LILYTH; ++k) {
                            const double zp = predecir_un_paso_causal(mod, out_prev.data(), n_prev);
                            out_prev[n_prev++] = zp;
                        }
                    }
                    t = t_fin_gap + 1;
                }
            }

            // 2. Z2: perfil_est en [t_ini_perf, t_apertura]
            // Fase 3: cualquier barra con FLAG_PUENTE_ITO o FLAG_CHOQUE_EN_HUECO en Z2 también cuenta como perfil corrupto
            std::size_t racha_inv_z2 = 0;
            std::size_t max_racha_z2 = 0;
            for (std::int64_t t = t_ini_perf; t <= r.t_apertura; ++t) {
                std::size_t idx = 0;
                const bool tiene_idx = indice(t, idx);
                const bool veto_ito_z2 = tiene_idx && ((hist_flags_hueco_[idx] & (ark::dataset::FLAG_PUENTE_ITO | ark::dataset::FLAG_CHOQUE_EN_HUECO | ark::dataset::FLAG_PAUSA_UMBRAL | ark::dataset::FLAG_CORRUPCION_CAIDA)) != 0u);
                const bool ok_t = tiene_idx && !veto_ito_z2 && h_val[idx] != 0u && std::isfinite(h_est[idx]);
                if (ok_t) {
                    racha_inv_z2 = 0;
                    if (n_perf < MAX_PERFIL_EVENTO_LILYTH) out_perf[n_perf++] = h_est[idx];
                } else {
                    ++racha_inv_z2;
                    if (racha_inv_z2 > max_racha_z2) max_racha_z2 = racha_inv_z2;
                    if (cfg_.politica_z2_piloto == ConfiguracionLilyth::PoliticaZ2Piloto::PUENTE_CAUSAL &&
                        racha_inv_z2 <= cfg_.r_max_z2_piloto && n_perf < MAX_PERFIL_EVENTO_LILYTH) {
                        std::array<double, MAX_PERFIL_EVENTO_LILYTH * 2> concat_loc{};
                        std::size_t n_c = 0;
                        for (std::size_t k = 0; k < n_prev && n_c < concat_loc.size(); ++k) concat_loc[n_c++] = out_prev[k];
                        for (std::size_t k = 0; k < n_perf && n_c < concat_loc.size(); ++k) concat_loc[n_c++] = out_perf[k];
                        out_perf[n_perf++] = predecir_un_paso_causal(mod, concat_loc.data(), n_c);
                    }
                }
            }
            if (cfg_.politica_z2_piloto == ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO &&
                max_racha_z2 > cfg_.r_max_z2_piloto) {
                flag_corrupto_z2 = true;
            }
        };

        copiar_rep(hist_est_cruda_, hist_val_est_cruda_, nulo_cruda_.modelo_lineal,
                   r.hist_prev_cruda_est, r.n_hist_prev_cruda,
                   r.perfil_cruda_est, r.n_perfil_cruda, r.perfil_corrupto_z2_cruda);
        copiar_rep(hist_est_purgada_, hist_val_est_purgada_, nulo_purgada_.modelo_lineal,
                   r.hist_prev_purgada_est, r.n_hist_prev_purgada,
                   r.perfil_purgada_est, r.n_perfil_purgada, r.perfil_corrupto_z2_purgada);
    }

    void abrir(const EventoAlarma& a, int j, double yc, bool okc, double yp, bool okp, SalidaPaso& out) noexcept {
        ++cont_.aperturas;
        std::size_t slot = RANURAS;
        for (std::size_t i = 0; i < RANURAS; ++i) if (!ranuras_[i].ocupada) { slot = i; break; }
        if (slot == RANURAS) { ++cont_.desbordes; return; }
        Ranura& r = ranuras_[slot];
        r.limpiar();
        r.ocupada = true;
        r.id_evento = ++seq_;
        r.t_apertura = t_;
        r.t_cambio = a.t_cambio;
        r.signo = a.signo;
        r.escala_origen = j;
        r.amplitud_glr = a.amplitud; r.se_amplitud_glr = a.se_amplitud; r.g_disparo = a.g;
        const FichaEscala& f = fichas_[static_cast<std::size_t>(j)];
        r.ficha_apertura.canal = 0;
        r.ficha_apertura.escala_origen = j;
        r.ficha_apertura.alarma = a;
        r.ficha_apertura.ficha_origen = f;
        r.ancla_saturada_en_l = a.l >= f.det.L;
        const std::size_t m = static_cast<std::size_t>(f.m >= 1 ? f.m : 1);
        r.w_ancla = Rasgos::n_reg * m;
        const double ti = std::isfinite(f.nulo.tau_int) && f.nulo.tau_int >= 1.0 ? f.nulo.tau_int : 1.0;
        r.guarda_ancla = static_cast<std::size_t>(std::ceil(ti)) * m;
        r.t_fin_ancla = r.t_cambio - static_cast<std::int64_t>(r.guarda_ancla) - 1;
        const std::int64_t t_ini = r.t_fin_ancla - static_cast<std::int64_t>(r.w_ancla) + 1;
        std::size_t ic = 0;
        if (indice(r.t_fin_ancla, ic)) r.medianas_capa1_congeladas = hist_capa1_[ic];
        copiar_perfil_estacionario(r);
        preparar_rama(r.purgada, true, medir_purgada(), t_ini, r.t_fin_ancla, yp, okp, r.signo);
        preparar_rama(r.cruda, false, medir_cruda(), t_ini, r.t_fin_ancla, yc, okc, r.signo);
        if (!r.alguna_rama_abierta()) { cerrar(r, out); return; }
        out.abrio_ranura = true;
        out.indice_ranura_abierta = slot;
    }

    [[nodiscard]] double nu_punto_ventana_reciente(std::size_t w_anc) const noexcept {
        if (!hay_reloj_ || w_anc == 0 || t_ < 0) return NaN;
        const std::size_t n_disp = std::min<std::size_t>(w_anc, n_hist_);
        if (n_disp == 0) return NaN;
        double snu = 0.0;
        for (std::size_t k = 0; k < n_disp; ++k) {
            const std::size_t pos = static_cast<std::size_t>(t_ - static_cast<std::int64_t>(k)) % HIST_MAX;
            snu += hist_dnu_[pos];
        }
        return snu / static_cast<double>(n_disp);
    }

    void avanzar(RamaTrayectoriaLilyth<T_MAX>& rama, std::size_t tau, double y, bool ok, bool purgada, std::size_t w_anc = 120) noexcept {
        if (rama.desenlace != DesenlaceRanura::ABIERTA) return;
        if (!ok || !std::isfinite(y)) {
            if (cfg_.politica_z3_piloto == ConfiguracionLilyth::PoliticaZ3Piloto::CORTA_AL_PRIMERO) {
                rama.desenlace = std::isfinite(rama.tau_rest_consolidado) ? DesenlaceRanura::RESTITUIDA : DesenlaceRanura::CENSURADA_DATO_INVALIDO;
                rama.tau_fin = std::isfinite(rama.tau_rest_consolidado) ? static_cast<std::size_t>(rama.tau_rest_consolidado) : rama.pasos_observados;
                return;
            }
            if (cfg_.politica_z3_piloto == ConfiguracionLilyth::PoliticaZ3Piloto::MUTANTE_COMPRIME_TAU) {
                // Salta el paso invalido comprimiendo el indice de tau en los pasos siguientes; cierra si se alcanza el horizonte fisico
                if (tau >= cfg_.t_max_efectivo) {
                    const bool rest = std::isfinite(rama.tau_rest_consolidado);
                    rama.desenlace = rest ? DesenlaceRanura::RESTITUIDA : DesenlaceRanura::CENSURADA_HORIZONTE;
                    rama.tau_fin = rest ? static_cast<std::size_t>(rama.tau_rest_consolidado) : rama.pasos_observados;
                }
                return;
            }
            // RELOJ_FISICO_TOPE_R: mantiene el reloj fisico tau = t - t_a sin comprimir
            ++rama.racha_invalida_z3;
            if (rama.racha_invalida_z3 > cfg_.r_max_z3_piloto) {
                rama.desenlace = std::isfinite(rama.tau_rest_consolidado) ? DesenlaceRanura::RESTITUIDA : DesenlaceRanura::CENSURADA_DATO_INVALIDO;
                rama.tau_fin = std::isfinite(rama.tau_rest_consolidado) ? static_cast<std::size_t>(rama.tau_rest_consolidado) : rama.pasos_observados;
                return;
            }
            rama.eta[tau - 1] = NaN;
            rama.pasos_observados = tau;
            if (tau >= cfg_.t_max_efectivo) {
                const bool rest = std::isfinite(rama.tau_rest_consolidado);
                rama.desenlace = rest ? DesenlaceRanura::RESTITUIDA : DesenlaceRanura::CENSURADA_HORIZONTE;
                rama.tau_fin = rest ? static_cast<std::size_t>(rama.tau_rest_consolidado) : tau;
            }
            return;
        }
        rama.racha_invalida_z3 = 0;
        ++rama.pasos_validos_z3;
        const std::size_t tau_ef = (cfg_.politica_z3_piloto == ConfiguracionLilyth::PoliticaZ3Piloto::MUTANTE_COMPRIME_TAU)
                                       ? rama.pasos_validos_z3
                                       : tau;
        const double eta = (y - rama.y_ancla) / rama.dy;
        rama.eta[tau_ef - 1] = eta;
        rama.pasos_observados = tau_ef;
        if (purgada) { suma_eta_purgada_[tau_ef - 1] += eta; ++cuenta_eta_purgada_[tau_ef - 1]; }
        else { suma_eta_cruda_[tau_ef - 1] += eta; ++cuenta_eta_cruda_[tau_ef - 1]; }
        if (tau_ef >= 2) {
            if (!std::isfinite(rama.min_eta_desde_2) || eta < rama.min_eta_desde_2) { rama.min_eta_desde_2 = eta; rama.tau_min_eta = tau_ef; }
            if (eta < 0.0) rama.armonico_sospechoso = true;
        }
        if (eta <= 0.5) {
            const double dnu = hay_reloj_ ? nu_acum_ - rama.nu_en_disparo : NaN;
            const double dp = pared_ - rama.pared_en_disparo;
            const double nu_p_movil = nu_punto_ventana_reciente(w_anc);
            const double nu_p_ef = detalle_reloj_marea::seleccionar_nu_punto_ancla(rama.nu_punto_ancla, nu_p_movil);
            if (!std::isfinite(rama.tau_rest_primer_toque)) {
                rama.tau_rest_primer_toque = static_cast<double>(tau_ef);
                rama.delta_nu_primer_toque = dnu;
                if (nu_p_ef > 0.0) rama.tau_act_primer_toque = detalle_reloj_marea::calcular_tau_actividad(dnu, nu_p_ef, tau_ef);
                rama.pared_primer_toque = dp;
            }
            if (!calor_activo_ && !std::isfinite(rama.tau_rest_consolidado)) {
                rama.tau_rest_consolidado = static_cast<double>(tau_ef);
                rama.delta_nu_consolidado = dnu;
                if (nu_p_ef > 0.0) rama.tau_act_consolidado = detalle_reloj_marea::calcular_tau_actividad(dnu, nu_p_ef, tau_ef);
                rama.pared_consolidado = dp;
            }
        }
        if (tau >= cfg_.t_max_efectivo) {
            const bool rest = std::isfinite(rama.tau_rest_consolidado);
            rama.desenlace = rest ? DesenlaceRanura::RESTITUIDA : DesenlaceRanura::CENSURADA_HORIZONTE;
            rama.tau_fin = rest ? static_cast<std::size_t>(rama.tau_rest_consolidado) : tau_ef;
        }
    }

    void cortar_por_replica(Ranura& r, std::size_t tau, int signo_replica) noexcept {
        bool corto = false;
        for (auto* rama : {&r.purgada, &r.cruda}) {
            if (rama->desenlace != DesenlaceRanura::ABIERTA) continue;
            corto = true;
            if (std::isfinite(rama->tau_rest_consolidado)) {
                rama->desenlace = DesenlaceRanura::RESTITUIDA;
                rama->tau_fin = static_cast<std::size_t>(rama->tau_rest_consolidado);
            } else {
                rama->desenlace = DesenlaceRanura::CORTADA_POR_REPLICA;
                rama->tau_fin = tau;
            }
        }
        if (corto) r.signo_replica_corte = signo_replica;
    }

    void cerrar(Ranura& r, SalidaPaso& out) noexcept {
        r.delta_tau_marea_primer_toque = detalle_reloj_marea::diferencia_tau_marea(
            r.cruda.tau_rest_primer_toque, r.purgada.tau_rest_primer_toque);
        r.delta_tau_marea_consolidado = detalle_reloj_marea::diferencia_tau_marea(
            r.cruda.tau_rest_consolidado, r.purgada.tau_rest_consolidado);
        auto registrar = [&](ResumenAalenJohansenSigno<T_MAX>& aj, const RamaTrayectoriaLilyth<T_MAX>& rama, std::size_t& excl) {
            if (rama.desenlace == DesenlaceRanura::INACTIVA || es_no_medida(rama.desenlace)) return;
            aj.registrar(r.signo, rama.desenlace, rama.tau_fin);
            if (!std::isfinite(rama.tau_rest_consolidado)) ++excl;
        };
        registrar(aj_purgada_, r.purgada, excluidos_reescalada_purgada_);
        registrar(aj_cruda_, r.cruda, excluidos_reescalada_cruda_);
        if (out.ranuras_cerradas < RANURAS) archivo_[out.ranuras_cerradas++] = r;
        if (cerradas_historico_.size() < cerradas_historico_.capacity()) cerradas_historico_.push_back(r);
        else ++cont_.historico_desbordado;   // Code: tope de memoria del historico (antes se perdian en silencio)
        r.limpiar();
    }
};

}  // namespace ark::lilyth