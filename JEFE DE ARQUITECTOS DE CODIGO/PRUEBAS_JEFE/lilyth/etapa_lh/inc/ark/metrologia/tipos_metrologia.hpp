// tipos_metrologia.hpp -- contratos del bloque de saltos y temperatura (ARK-METRO-01), C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04, por orden del Comandante (conversacion del 2026-11-04). Prerregistro de
// los controles: JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/metrologia/EXPECTATIVAS.md.
//
// El bloque ingiere cualquier observable theta(t) a su resolucion nativa y responde: donde cambio, cuanto, de que forma (nivel o
// temperatura) y en que escalas. Es telemetria: no decide operar.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace ark::metrologia {

inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
inline constexpr double ALFA_PROYECTO = 0.05;   // nivel nominal del proyecto (EXPECTATIVAS §1)

// Orden de integracion del observable. NIVEL: precios, nu acumulado... -> se trabaja sobre sus incrementos nativos (SUMA).
enum class Integracion : std::uint8_t { FLUJO = 0, NIVEL = 1 };

// Como se agrega a escalas mayores la representacion estacionaria.
enum class AgregacionEscala : std::uint8_t {
    ULTIMO,               // intensivos muestreados (un nivel ya convertido no la usa)
    SUMA,                 // retornos, incrementos de un nivel, delta nu, volumen
    MEDIA,                // intensivos
    RAZON_DE_SUMAS,       // ofi_ratio, taker_buy_ratio, ticket: sum(num)/sum(den)
    LOG_DE_SUMA,          // volumen log: la entrada es el volumen BASE, se emite log(sum)
    RAIZ_SUMA_CUADRADOS   // volatilidad: la entrada es sigma, se emite sqrt(sum sigma^2)
};

enum class EstadoSaltos : std::uint8_t {
    SIN_CALIBRAR, NO_DETECTABLE, CALMA, SALTO, TRANSITORIO, DERIVA, CAMBIO_TEMPERATURA
};
inline constexpr std::array<const char*, 7> NOMBRE_ESTADO{
    "SIN_CALIBRAR", "NO_DETECTABLE", "CALMA", "SALTO", "TRANSITORIO", "DERIVA", "CAMBIO_TEMPERATURA"};

enum class NormaRegularizador : std::uint8_t { SIN_DATO, L2_TIKHONOV, L1_TV };

// Modo de deteccion del operador delta (orden del Comandante, 2026-11-05: los dos quedan).
//   DISJUNTO   cada escala alarma al cerrar su barra 2^j; la temperatura, al cerrar su ventana de T* barras (ventanas independientes).
//   DESLIZANTE el nivel se evalua en cada paso nativo en todas las escalas y la temperatura en cada barra de su escala
//              (detector_deslizante.hpp). La superestadistica, el regularizador y la ventana de forma siguen sobre barras disjuntas.
enum class ModoDeteccion : std::uint8_t { DISJUNTO = 0, DESLIZANTE = 1 };

// Rejilla de REPORTE de la curva de Lorden y del tau50 (delta en unidades de s0). No decide nada.
inline constexpr std::array<double, 5> DELTA_REPORTE{0.25, 0.5, 1.0, 2.0, 4.0};

// ---------------- fichas formales (obligatorias en cada estado) ----------------

struct FichaNulo {
    std::size_t pasos = 0;          // pasos validos del nulo de la escala
    double mu0 = NaN, s0 = NaN, acf1 = NaN;
    double tau_int = NaN;           // 1 + 2 sum_{k<K} rho_k
    std::size_t rezago_corte = 0;   // K: primer rezago con |rho_k| dentro de la banda del nulo
};

struct FichaDetector {
    std::size_t L = 0;                       // memoria maxima del GLR (pasos)
    double arl0_pasos = NaN;                 // ARL0 objetivo de ESTA escala (ya repartido entre escalas)
    double h = NaN;
    double h_confirmacion = NaN;             // solo la temperatura DESLIZANTE: umbral de confirmacion por acumulacion (ARL0 / alfa)
    bool h_extrapolado = false;
    std::size_t alarmas_nulo = 0;            // alarmas del nulo en h
    double lambda_tv = NaN;                  // cuantil 1 - alfa del maximo del puente en ventanas disjuntas del nulo
    std::size_t ventanas_lambda = 0;
    double alfa = ALFA_PROYECTO;
    std::array<double, DELTA_REPORTE.size()> t_lorden{};          // pasos de la escala
    std::array<double, DELTA_REPORTE.size()> tau50{};             // pasos de la escala, medido
    std::array<double, DELTA_REPORTE.size()> fraccion_sin_detectar{};
};

struct FichaEscala {
    int j = 0;                    // nivel de la escalera
    std::int64_t m = 1;           // barras nativas por barra de la escala (2^j)
    bool detectable = false;
    const char* motivo = "sin calibrar";
    FichaNulo nulo{};
    FichaDetector det{};
};

// Ficha de superestadistica de una escala (diagnostico de calibracion).
struct FichaSuperestadistica {
    int j = 0;
    std::size_t T = 0;                  // ventana de equilibrio local T* (0 si no hay)
    bool t1_equilibrio_local = false;   // curtosis local = gaussiana al nivel alfa en T*
    double exceso_curtosis = NaN;       // en T* (o en la menor T si T1 falla)
    double exceso_curtosis_umbral = NaN;
    bool t2_separacion = false;         // memoria de log beta entre ventanas (TELEMETRIA: no entra en el veredicto)
    double rho_log_beta = NaN, rho_umbral = NaN;
    double tau_int_log_beta = NaN;      // en ventanas
    double t_beta_pasos = NaN;          // T * tau_int
    std::size_t ventanas = 0;
    double media_log_beta = NaN, var_log_beta = NaN;   // corregidas por el error de la ventana
    double beta_media = NaN, q_tsallis = NaN;          // <beta>, <beta^2>/<beta>^2 corregidos
    // Separacion de escalas (decision del Comandante, 2026-11-04): con separacion, Var[log beta] de ventanas de L_k = T* 2^k barras
    // se mantiene (la mezcla de ventanas no cambia f mientras L_k << T_beta); sin ella cae. Estadistico: r_k = V_k / V_0, contra los
    // nulos ELEGIDOS por el usuario (cada uno con B sustitutos, p de una cola); Holm sobre los nulos elegidos en k = 1.
    static constexpr std::size_t K_SEP = 16;
    std::size_t k_sep = 0;                              // escalones k = 0..k_sep-1 medidos (>= ceil(1/alfa) ventanas cada uno)
    std::array<double, K_SEP> v_k{};                     // Var[log beta] corregida por la ventana, en L_k
    std::uint32_t nulos = 0;                            // mascara de nulos aplicados
    std::size_t replicas_nulo = 0;
    std::array<std::array<double, K_SEP>, 3> p_k{};      // por nulo (PERMUTACION, BLOQUES, SIGNOS) y k: p de s_k = V_{k+1}/V_k (enmienda 7:
                                                         // BLOQUES con b = L_k = T* 2^k; "beta persiste mas alla de L_k")
    std::array<std::array<double, K_SEP>, 3> r_nulo_mediana{};   // mediana de s_k bajo el nulo
    std::array<double, 3> p_holm{};                      // p de k = 0 ajustada por Holm entre los nulos elegidos
    std::array<std::size_t, 3> escala_separacion{};      // L_{k+1} del ultimo k rechazado contiguo desde k = 0 (0: no persiste mas alla de T*)
    bool separacion_contrastada = false;                // false si el usuario no eligio ningun nulo de separacion
    bool separacion = false;                            // todos los nulos elegidos rechazados (Holm, alfa) en k = 0
    bool es_superestadistica = false;                   // T1 y separacion (T2 queda como telemetria: la memoria de beta en t -> t+1
                                                        // la mide el operador de coordinacion)
};

struct EventoAlarma {
    bool hay = false;
    int signo = 0;
    std::int64_t t = -1;            // paso nativo de la alarma (ultima barra nativa incluida)
    std::int64_t t_cambio = -1;     // punto de cambio estimado (paso nativo)
    std::size_t l = 0;              // memoria elegida (pasos de la escala)
    double amplitud = NaN;          // cambio medio por barra de la escala (unidades de y); sesgado al alza (declarado)
    double se_amplitud = NaN;
    double g = NaN;                 // estadistico GLR en la alarma
    double t_pendiente = NaN;       // t de la pendiente MCO de los l* ultimos z (varianza de largo plazo): forma del cambio
};

// Evento clasificado entre escalas, con su ficha completa.
struct FichaEvento {
    EstadoSaltos estado = EstadoSaltos::CALMA;
    int canal = 0;                  // 0 nivel, 1 temperatura
    int escala_origen = -1;
    int escala_confirmacion = -1;   // escala de PERSISTENCIA: la mas gruesa que confirmo en cascada (-1: ninguna)
    int escala_esperada = -1;       // escala cuya confirmacion se espera (evento abierto)
    bool confirmado = false;        // SALTO: persistio hasta la escala mas gruesa calibrada
    EventoAlarma alarma{};          // la alarma que abrio el evento
    std::int64_t plazo = -1;        // paso nativo limite de confirmacion
    std::int64_t t_cierre = -1;     // paso nativo en que se clasifico
    double kappa = NaN;             // tau50/t* usado en el plazo
    double t_lorden_confirmacion = NaN;
    double t_pendiente_cierre = NaN;   // t de la pendiente de z (escala de origen) desde la alarma hasta el cierre: forma del cambio
    std::size_t barras_forma = 0;      // barras de la escala de origen usadas para esa pendiente
    FichaEscala ficha_origen{};
};

}  // namespace ark::metrologia
