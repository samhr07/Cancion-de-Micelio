// tipos_fractal.hpp -- contratos del Bloque Fractal Universal (ARK-FRAC-01), C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-03, por orden del Comandante: mensajes/MANDATO_CONSOLIDACION_BLOQUE_FRACTAL_CXX.md,
// mensajes/DOCTRINA_COMANDANTE_BLOQUE_FRACTAL_UNIVERSAL.md y mensajes/MANDATO_INTEGRACION_REGULARIZADOR_DUAL_FISICO_CXX.md, con las
// decisiones del Comandante del 2026-11-03 sobre las tres premisas medidas (acta
// JEFE DE ARQUITECTOS DE CODIGO/AUDITORIA_CAUSAL/BLOQUE_FRACTAL_UNIVERSAL_CXX.md).
//
// DECLARACION METROLOGICA (EL ABISMO - COTA I & II):
// Observador de estado metrologico del medio. Queda terminantemente prohibido
// usar las salidas de H en escalas < 1h como gatillos reactivos directos taker.
// Su funcion es la telemetria del espacio de fases, la modulacion de matrices de covarianza
// y la calibracion de Kelly en horizontes macro (> 4 horas).
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace ark::fractal {

inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

// Escalas canonicas de la piramide (segundos por barra; D1 = dia UTC).
enum class Escala : std::uint8_t { M1 = 0, M5, M15, H1, H4, D1 };
inline constexpr std::size_t NUM_ESCALAS = 6;
inline constexpr std::array<std::int64_t, NUM_ESCALAS> SEGUNDOS_ESCALA{60, 300, 900, 3600, 14400, 86400};
inline constexpr std::array<const char*, NUM_ESCALAS> NOMBRE_ESCALA{"1m", "5m", "15m", "1h", "4h", "D1"};

// Orden de integracion del observable (la doctrina: cualquier senal causal).
//  NIVEL (d = 1): precios, profundidad acumulada, etc. El DFA se aplica al nivel (como el sensor de EL ABISMO: calcular_dfa sobre
//                 precios) y H = alpha - 1.
//  FLUJO (d = 0): OFI, signo de transacciones, log-volumen, log-spread, log-volatilidad. El DFA se aplica a la serie y H = alpha.
enum class Integracion : std::uint8_t { FLUJO = 0, NIVEL = 1 };

// Como se agrega el observable al pasar de la barra base a una escala mayor (piramide).
//  ULTIMO: niveles (cierre de la barra). SUMA: flujos aditivos (OFI, volumen). MEDIA: intensivos (spread, volatilidad).
enum class Agregacion : std::uint8_t { ULTIMO, SUMA, MEDIA };

// De donde sale el nulo del observable.
//  HISTORICO: linea base congelada (K = 0) sobre la historia causal del propio observable: mu0, s0, lambda_TV y h_crit miden
//             el reposo de X tal como es (vale para observables cuyo reposo no es H = 0.5: volumen, profundidad).
//  SIGNOS:    nulo real de EL ABISMO: los incrementos del observable con signo aleatorio (H = 0.5 con la marginal y la
//             heterocedasticidad de X). Es el nulo de los precios.
enum class NuloCalibracion : std::uint8_t { HISTORICO, SIGNOS };

// Norma que gobierna el borde presente, elegida por la propia innovacion (doctrina §3.C).
enum class Norma : std::uint8_t { SIN_CALIBRAR, CALENTANDO, L2_TIKHONOV, L1_TV };
inline constexpr std::array<const char*, 4> NOMBRE_NORMA{"SIN_CALIBRAR", "CALENTANDO", "L2_TIKHONOV", "L1_TV"};

// Salida del sensor DFA de una ventana (Capa 1).
struct ResultadoSensor {
    double alpha = NaN;   // pendiente log F(s) - log s
    double h = NaN;       // alpha - d
    double r2 = NaN;      // R^2 del ajuste log-log
    double se = NaN;      // error estandar de alpha (OLS)
    bool valido = false;  // false si la ventana trae no finitos o varianza nula (ErrorDFA en Python)
};

// Calibracion del observable (todo derivado; ningun umbral a mano salvo los declarados en ConfiguracionFractal).
struct CalibracionFractal {
    bool lista = false;
    NuloCalibracion nulo = NuloCalibracion::HISTORICO;
    std::size_t pasos_calibracion = 0;   // pasos validos del sensor en el nulo de calibracion
    double mu0 = NaN;                    // media del sensor en reposo
    double s0 = NaN;                     // dispersion de reposo del sensor
    double acf1 = NaN;                   // autocorrelacion de rezago 1 del sensor en reposo
    double lambda_tv = NaN;              // q_{1-alfa} de max|S_k| (puente de sumas parciales) en ventanas disjuntas de reposo
    double lambda_dj = NaN;              // s0*sqrt(2 ln N) (Donoho-Johnstone): SOLO telemetria (parte el reposo, acta §2)
    std::size_t ventanas_lambda = 0;     // ventanas disjuntas usadas para el cuantil
    double ganancia_sensor = NaN;        // E[H^ | 0.5 + dH] - E[H^ | 0.5] medido por inyeccion (mundo de signos de X)
    double delta = NaN;                  // ganancia_sensor / s0 (unidades de nu)
    double k = NaN;                      // holgura CUSUM = delta / 2
    double h_crit = NaN;                 // umbral CUSUM bilateral para el ARL0 objetivo, sobre el nulo elegido
    double arl0_objetivo_pasos = NaN;
    bool h_extrapolado = false;          // true si el ARL0 objetivo excede lo observable en el nulo (extrapolacion log-lineal)
    std::size_t alarmas_nulo_en_h = 0;   // alarmas del nulo de calibracion en el h elegido
    double h_siegmund = NaN;             // referencia iid (operador_decision.hpp); NO se usa: con ACF1 ~ 0.94 subestima el ARL0
    double tau50_pasos = NaN;            // retardo mediano de la alarma tras un salto dH (mundo de signos con la marginal de X)
    double tau50_barras = NaN;
    double fraccion_sin_detectar = NaN;  // saltos no detectados en el horizonte simulado
};

// Estado del observador en el borde presente (Capa 2), emitido en cada paso del sensor.
struct EstadoFractal {
    std::int64_t barra = -1;           // indice (en barras de su escala) de la ultima barra incluida
    ResultadoSensor sensor{};          // Capa 1 en el borde presente
    double nu = NaN;                   // (H^ - mu0) / s0
    double cusum_pos = 0.0, cusum_neg = 0.0;
    bool alarma_pos = false, alarma_neg = false;
    Norma norma = Norma::SIN_CALIBRAR;
    double h_regularizado = NaN;       // x[N-1] del resolvente de la norma elegida (borde presente, sin mirar el futuro)
    double lambda_usado = NaN;         // Morozov (L2) o lambda_tv (L1)
    std::size_t n_ventana_regularizador = 0;
    std::int64_t pasos_desde_alarma = -1;   // -1: sin alarma aun
    double tau50_barras = NaN;         // retardo de transporte medido de esta escala, G(s) = exp(-s tau50)
};

// Referencia de EL ABISMO (§1, matriz de viabilidad, sensor sobre PRECIOS, dH = 0.25): retardo del detector CUSUM en barras,
// mundo sintetico / nulo real (ARL0 = 1 mes; D1 = 1 ano, punto medio de los rangos). No es el retardo de Condat (acta §1).
// Solo para comparar: el observador exporta su tau50 MEDIDO para su propio observable.
struct ReferenciaTau50 { const char* escala; double sintetico_barras; double real_barras; };
inline constexpr std::array<ReferenciaTau50, 5> TAU50_REFERENCIA_PRECIOS{{
    {"1m", 207.0, 464.0}, {"1h", 79.0, 155.0}, {"4h", 49.0, 91.0}, {"D1 SPX/GC/CL/JPY", 60.0, 93.5}, {"D1 VIX/BTC", 59.5, 104.0}}};

}  // namespace ark::fractal
