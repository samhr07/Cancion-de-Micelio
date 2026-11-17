// extractor_cxx.hpp -- OBSERVABLES SOBRE LA REJILLA DIARIA (ARK-DATA-01-CXX)
// Reactor ARK / Cancion del Micelio -- ARQUITECTO DE DATOS (Flash) y correccion del Jefe 2026-11-01
//
// Se calculan sobre la REJILLA (construir_rejilla_diaria), donde el cierre de una fila sin sesion va arrastrado:
//   RETORNOS_LOG     r_t = ln C_t - ln C_(t-1), r_0 = 0 (extractor_observables.RETORNOS_LOG). Una fila sin sesion da 0 y la
//                    sesion siguiente cobra el movimiento ENTERO desde el ultimo cierre con dato. Antes el retorno siguiente a
//                    una fila vacia se ponia a 0: en DXY se perdian 2 884 retornos, el 21.5 % de la suma de |r|.
//   LOG_PRECIO, VOLUMEN_BASE, VOLUMEN_LOG, VOLUMEN_NOCIONAL: identicos a Python (NaN donde Python da NaN).
//   RETORNO_APERTURA_A_CIERRE, VOLATILIDAD_PARKINSON, VOLATILIDAD_GARMAN_KLASS, RETORNO_GAP_NOCTURNO: SOLO en filas con sesion
//                    (0 en las demas). DIVERGE de Python a proposito: extractor_observables los calcula tambien sobre la fila
//                    arrastrada y REPITE el valor de la ultima sesion en cada feriado y fin de semana (un dato que no existio).
//                    En las filas con sesion los numeros son identicos a Python.
// El nocional es close * volume, como en Python. OJO: en BTC/ETH de Yahoo el volumen ya viene en USD, asi que ahi el nocional
// no es un nocional (USD^2): usar VOLUMEN_BASE.
#pragma once

#include "tipos_dataset.hpp"

#include <cmath>
#include <cstdint>
#include <cstddef>
#include <algorithm>

namespace ark::dataset {

enum class TipoObservableCxx : uint8_t {
    PRECIO_CLOSE = 0,
    RETORNOS_LOG = 1,
    RETORNO_APERTURA_A_CIERRE = 2,
    RETORNO_GAP_NOCTURNO = 3,
    LOG_PRECIO = 4,
    VOLUMEN_BASE = 5,
    VOLUMEN_LOG = 6,
    VOLUMEN_NOCIONAL = 7,
    VOLATILIDAD_PARKINSON = 8,
    VOLATILIDAD_GARMAN_KLASS = 9
};

class ExtractorObservablesCxx {
    static constexpr double LN2 = 0.693147180559945309417;

public:
    // Mandato P <= 0 (2026-11-02): un precio no positivo o no finito no tiene logaritmo. La rejilla ya lo deja en NaN y arrastra
    // el ultimo cierre valido (blindaje e. de construir_rejilla_diaria, igual que control_calidad.py en Python); aqui se blinda
    // ademas a quien llame al extractor sobre barras que no pasaron por la rejilla: NaN, nunca log(max(p, 1e-12)). Con 1e-12 el
    // WTI de -37.63 (2020-04-20) fabricaba retornos de -30.5 y +29.9.
    // DIVERGE de extractor_observables (np.maximum(p, 1e-12)) solo con p <= 0, p infinito o 0 < p < 1e-12, que el blindaje del
    // pipeline no deja llegar: con los precios de la rejilla el resultado es bit a bit el mismo (contraste_rejilla.py).
    [[nodiscard]] static bool precio_valido(double p) noexcept { return std::isfinite(p) && p > 0.0; }
    [[nodiscard]] static double lp(double p) noexcept { return precio_valido(p) ? std::log(p) : NAN_DATASET; }
    // ln(a / b) con la misma forma que numpy (el cociente antes del logaritmo, para la identidad bit a bit).
    [[nodiscard]] static double lr(double a, double b) noexcept {
        return (precio_valido(a) && precio_valido(b)) ? std::log(a / b) : NAN_DATASET;
    }

    static void extraer_retornos_log(const BarraCausal* b, size_t n, double* out) noexcept {
        if (!b || !out || n == 0) return;
        out[0] = 0.0;
        for (size_t i = 1; i < n; ++i) out[i] = lp(b[i].close) - lp(b[i - 1].close);
    }

    static void extraer_retorno_apertura_a_cierre(const BarraCausal* b, size_t n, double* out) noexcept {
        if (!b || !out) return;
        for (size_t i = 0; i < n; ++i)
            out[i] = b[i].es_sesion_con_dato ? lr(b[i].close, b[i].open) : 0.0;
    }

    // ln O_t - ln C_(t-1) en filas con sesion: el cierre anterior es el arrastrado, o sea el de la ultima sesion con dato.
    static void extraer_retorno_gap(const BarraCausal* b, size_t n, double* out) noexcept {
        if (!b || !out || n == 0) return;
        out[0] = 0.0;
        for (size_t i = 1; i < n; ++i)
            out[i] = b[i].es_sesion_con_dato ? lr(b[i].open, b[i - 1].close) : 0.0;
    }

    static void extraer_log_precio(const BarraCausal* b, size_t n, double* out) noexcept {
        if (!b || !out) return;
        for (size_t i = 0; i < n; ++i) out[i] = lp(b[i].close);
    }

    static void extraer_volumen_base(const BarraCausal* b, size_t n, double* out) noexcept {
        if (!b || !out) return;
        for (size_t i = 0; i < n; ++i) out[i] = b[i].volume;
    }

    static void extraer_volumen_log(const BarraCausal* b, size_t n, double* out) noexcept {
        if (!b || !out) return;
        for (size_t i = 0; i < n; ++i) out[i] = std::log(std::max(b[i].volume, 0.0) + 1.0);
    }

    static void extraer_volumen_nocional(const BarraCausal* b, size_t n, double* out) noexcept {
        if (!b || !out) return;
        for (size_t i = 0; i < n; ++i) out[i] = b[i].volume * b[i].close;
    }

    // Parkinson: sigma = sqrt(ln(H/L)^2 / (4 ln 2)), solo en filas con sesion.
    static void extraer_volatilidad_parkinson(const BarraCausal* b, size_t n, double* out) noexcept {
        if (!b || !out) return;
        const double denom = 4.0 * LN2;
        for (size_t i = 0; i < n; ++i) {
            if (!b[i].es_sesion_con_dato) { out[i] = 0.0; continue; }
            const double lhl = lr(b[i].high, b[i].low);
            out[i] = std::isnan(lhl) ? NAN_DATASET : std::sqrt(std::max((lhl * lhl) / denom, 0.0));
        }
    }

    // Garman-Klass: sigma^2 = 0.5 ln(H/L)^2 - (2 ln 2 - 1) ln(C/O)^2, solo en filas con sesion.
    static void extraer_volatilidad_garman_klass(const BarraCausal* b, size_t n, double* out) noexcept {
        if (!b || !out) return;
        const double k = 2.0 * LN2 - 1.0;
        for (size_t i = 0; i < n; ++i) {
            if (!b[i].es_sesion_con_dato) { out[i] = 0.0; continue; }
            const double lhl = lr(b[i].high, b[i].low);
            const double lco = lr(b[i].close, b[i].open);
            out[i] = (std::isnan(lhl) || std::isnan(lco)) ? NAN_DATASET : std::sqrt(std::max(0.5 * (lhl * lhl) - k * (lco * lco), 0.0));
        }
    }

    static void extraer(const BarraCausal* b, size_t n, TipoObservableCxx tipo, double* out) noexcept {
        switch (tipo) {
            case TipoObservableCxx::PRECIO_CLOSE:
                if (b && out) for (size_t i = 0; i < n; ++i) out[i] = b[i].close;
                break;
            case TipoObservableCxx::RETORNOS_LOG: extraer_retornos_log(b, n, out); break;
            case TipoObservableCxx::RETORNO_APERTURA_A_CIERRE: extraer_retorno_apertura_a_cierre(b, n, out); break;
            case TipoObservableCxx::RETORNO_GAP_NOCTURNO: extraer_retorno_gap(b, n, out); break;
            case TipoObservableCxx::LOG_PRECIO: extraer_log_precio(b, n, out); break;
            case TipoObservableCxx::VOLUMEN_BASE: extraer_volumen_base(b, n, out); break;
            case TipoObservableCxx::VOLUMEN_LOG: extraer_volumen_log(b, n, out); break;
            case TipoObservableCxx::VOLUMEN_NOCIONAL: extraer_volumen_nocional(b, n, out); break;
            case TipoObservableCxx::VOLATILIDAD_PARKINSON: extraer_volatilidad_parkinson(b, n, out); break;
            case TipoObservableCxx::VOLATILIDAD_GARMAN_KLASS: extraer_volatilidad_garman_klass(b, n, out); break;
            default: break;
        }
    }
};

}  // namespace ark::dataset
