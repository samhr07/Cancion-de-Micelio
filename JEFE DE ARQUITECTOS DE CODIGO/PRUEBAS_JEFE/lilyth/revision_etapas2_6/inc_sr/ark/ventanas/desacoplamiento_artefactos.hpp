// desacoplamiento_artefactos.hpp -- ALGEBRA MULTIESCALA, KANTELHARDT Y PROYECTORES ORTOGONALES (ARK-WIN-01-CXX)
// Reactor ARK / Canción del Micelio -- ARQUITECTO DE DATOS & COMANDANTE
// Corregido por el Jefe de Arquitectos de Código (2026-11-01), E9: un orden de tendencia no soportado ya no degrada en
// silencio a otro. El proyector y la fluctuación F(s) soportan los órdenes 0, 1 y 2; cualquier otro devuelve false / NaN.
// La fluctuación de orden 2 (DFA-2) quita la parábola de cada segmento con los polinomios de Gram de su longitud.
#pragma once

#include "ark/ventanas/contratos_ventanas.hpp"
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <array>
#include <span>
#include <algorithm>
#include <functional>
#include <limits>

namespace ark::ventanas {

// =======================================================================
// 1. Proyectores Ortogonales Aniquiladores de Tendencias (Axioma 6)
// Implementación cerrada O(W) vía polinomios discretos ortogonales de Gram
// =======================================================================

template <size_t W>
class ProyectorOrtogonalCxx {
    static_assert(W >= 3, "lookback_W debe ser al menos 3 para proyeccion ortogonal.");

private:
    static constexpr double W_d = static_cast<double>(W);
    static constexpr double t_bar = 0.5 * static_cast<double>(W - 1);
    static constexpr double p1_norm_sq = W_d * (W_d * W_d - 1.0) / 12.0;
    static constexpr double p2_shift = (W_d * W_d - 1.0) / 12.0;
    static constexpr double p2_norm_sq = W_d * (W_d * W_d - 1.0) * (W_d * W_d - 4.0) / 180.0;

public:
    // Aplica proyección aniquiladora M_X y = y - P y en tiempo O(W) sin asignaciones heap
    // orden = 0: remoción de media (P0); 1: lineal (P0 + P1); 2: cuadrática (P0 + P1 + P2).
    // Devuelve false (y no escribe) si el orden no está soportado o faltan punteros.
    static bool aplicar_detrending(
        const double* in_valores,
        double* out_detrend,
        int orden = 1
    ) noexcept {
        if (!in_valores || !out_detrend || orden < 0 || orden > 2) return false;

        double sum_y = 0.0;
        double sum_y_p1 = 0.0;
        double sum_y_p2 = 0.0;

        for (size_t t = 0; t < W; ++t) {
            const double y = in_valores[t];
            sum_y += y;
            if (orden >= 1) {
                const double dt = static_cast<double>(t) - t_bar;
                sum_y_p1 += y * dt;
                if (orden >= 2) {
                    const double p2 = dt * dt - p2_shift;
                    sum_y_p2 += y * p2;
                }
            }
        }

        const double c0 = sum_y / W_d;
        const double c1 = (orden >= 1 && p1_norm_sq > 0.0) ? sum_y_p1 / p1_norm_sq : 0.0;
        const double c2 = (orden >= 2 && p2_norm_sq > 0.0) ? sum_y_p2 / p2_norm_sq : 0.0;

        for (size_t t = 0; t < W; ++t) {
            double p_fit = c0;
            if (orden >= 1) {
                const double dt = static_cast<double>(t) - t_bar;
                p_fit += c1 * dt;
                if (orden >= 2) {
                    const double p2 = dt * dt - p2_shift;
                    p_fit += c2 * p2;
                }
            }
            out_detrend[t] = in_valores[t] - p_fit;
        }
        return true;
    }

    // Construye la matriz analítica completa M_X (W x W) en O(W^2)
    // Garantiza simetría e idempotencia exacta bit a bit: M_X^2 == M_X y M_X^T == M_X
    static bool construir_matriz_mx(
        std::array<std::array<double, W>, W>& out_mx,
        int orden = 1
    ) noexcept {
        if (orden < 0 || orden > 2) return false;
        for (size_t i = 0; i < W; ++i) {
            const double di = static_cast<double>(i) - t_bar;
            const double p2_i = di * di - p2_shift;

            for (size_t j = 0; j < W; ++j) {
                const double dj = static_cast<double>(j) - t_bar;
                const double p2_j = dj * dj - p2_shift;

                double p_ij = 1.0 / W_d;
                if (orden >= 1 && p1_norm_sq > 0.0) {
                    p_ij += (di * dj) / p1_norm_sq;
                }
                if (orden >= 2 && p2_norm_sq > 0.0) {
                    p_ij += (p2_i * p2_j) / p2_norm_sq;
                }

                const double delta_ij = (i == j) ? 1.0 : 0.0;
                out_mx[i][j] = delta_ij - p_ij;
            }
        }
        return true;
    }
};

// =======================================================================
// 2. Segmentación Causal de Kantelhardt Bidireccional (Axioma 6)
// Doble escaneo anclado al borde presente t con garantía de descarte = 0.0%
// =======================================================================

template <size_t W>
class SegmentadorKantelhardtCxx {
    static_assert(W >= 4, "W debe ser al menos 4 muestras.");

public:
    // Resultado estructurado zero-heap de la segmentación
    struct ResultadoSegmentacion {
        size_t s{0};
        size_t n_segmentos{0};
        size_t residuo{0};
        double porcentaje_descarte{0.0};
        bool bidireccional{true};
    };

    // Computa las propiedades de la segmentación sin asignar memoria
    [[nodiscard]] static constexpr ResultadoSegmentacion inspeccionar(
        size_t s,
        bool bidireccional = true
    ) noexcept {
        if (s < 2 || s > W) {
            return {s, 0, 0, 100.0, bidireccional};
        }
        const size_t n_s = W / s;
        const size_t res = W % s;

        if (res == 0 || !bidireccional) {
            const double desc = (res == 0) ? 0.0 : (static_cast<double>(res) / static_cast<double>(W) * 100.0);
            return {s, n_s, res, desc, bidireccional};
        }

        // Doble escaneo Kantelhardt: 2 * n_s segmentos y 0.0% de descarte causal en el presente
        return {s, 2 * n_s, res, 0.0, bidireccional};
    }

    // Itera causalmente sobre todos los segmentos hacia adelante y hacia atrás
    // invocando el callback fn(const double* seg, size_t len, bool es_reverso, size_t k)
    template <typename Callback>
    static void iterar_segmentos(
        const double* valores,
        size_t s,
        bool bidireccional,
        Callback&& fn
    ) noexcept {
        if (!valores || s < 2 || s > W) return;
        const size_t n_s = W / s;
        const size_t res = W % s;

        // 1. Escaneo hacia adelante [0 : n_s * s]
        for (size_t k = 0; k < n_s; ++k) {
            const double* ptr_fwd = &valores[k * s];
            fn(ptr_fwd, s, false, k);
        }

        // 2. Doble escaneo Kantelhardt (si hay residuo y está activado)
        // Anclado rígidamente al borde presente W - 1 retrocediendo en el tiempo
        if (bidireccional && res > 0) {
            for (size_t k = 0; k < n_s; ++k) {
                // El segmento k=0 cubre [W - s, W], culminando exactamente en valores[W - 1]
                const double* ptr_rev = &valores[W - (k + 1) * s];
                fn(ptr_rev, s, true, k);
            }
        }
    }

    // Calcula directamente la fluctuación cuadrática F(s) de DFA sobre la ventana, quitando a cada segmento su polinomio de
    // orden 0, 1 o 2 (mínimos cuadrados, sin heap). Otro orden, o un orden 2 sobre segmentos de menos de 3 puntos: NaN.
    template <size_t MAX_S = W>
    [[nodiscard]] static double computar_fluctuacion_f(
        const double* valores,
        size_t s,
        int orden_detrend = 1,
        bool bidireccional = true
    ) noexcept {
        if (!valores || s < 2 || s > W) return 0.0;
        if (orden_detrend < 0 || orden_detrend > 2 || (orden_detrend == 2 && s < 3)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        double sum_f2_total = 0.0;
        size_t conteo_segmentos = 0;

        iterar_segmentos(valores, s, bidireccional, [&](const double* seg_ptr, size_t len, bool, size_t) {
            const double len_d = static_cast<double>(len);
            if (orden_detrend == 1) {
                const double sum_t = 0.5 * static_cast<double>((len - 1) * len);
                const double sum_t2 = static_cast<double>((len - 1) * len * (2 * len - 1)) / 6.0;
                double sum_y = 0.0, sum_ty = 0.0;
                for (size_t t = 0; t < len; ++t) {
                    sum_y += seg_ptr[t];
                    sum_ty += static_cast<double>(t) * seg_ptr[t];
                }
                const double denom = len_d * sum_t2 - sum_t * sum_t;
                const double b = (denom != 0.0) ? (len_d * sum_ty - sum_t * sum_y) / denom : 0.0;
                const double a = (sum_y - b * sum_t) / len_d;

                double var_seg = 0.0;
                for (size_t t = 0; t < len; ++t) {
                    const double res = seg_ptr[t] - (a + b * static_cast<double>(t));
                    var_seg += res * res;
                }
                sum_f2_total += (var_seg / len_d);
            } else if (orden_detrend == 2) {
                // Polinomios de Gram de longitud len: p1 = t - t_bar, p2 = p1^2 - (len^2 - 1)/12
                const double tb = 0.5 * (len_d - 1.0);
                const double sh = (len_d * len_d - 1.0) / 12.0;
                const double n1 = len_d * (len_d * len_d - 1.0) / 12.0;
                const double n2 = len_d * (len_d * len_d - 1.0) * (len_d * len_d - 4.0) / 180.0;
                double s0 = 0.0, s1 = 0.0, s2 = 0.0;
                for (size_t t = 0; t < len; ++t) {
                    const double d = static_cast<double>(t) - tb;
                    s0 += seg_ptr[t];
                    s1 += seg_ptr[t] * d;
                    s2 += seg_ptr[t] * (d * d - sh);
                }
                const double c0 = s0 / len_d, c1 = s1 / n1, c2 = s2 / n2;
                double var_seg = 0.0;
                for (size_t t = 0; t < len; ++t) {
                    const double d = static_cast<double>(t) - tb;
                    const double res = seg_ptr[t] - (c0 + c1 * d + c2 * (d * d - sh));
                    var_seg += res * res;
                }
                sum_f2_total += (var_seg / len_d);
            } else {
                // Grado 0 (media)
                double sum_y = 0.0;
                for (size_t t = 0; t < len; ++t) sum_y += seg_ptr[t];
                const double mean = sum_y / len_d;

                double var_seg = 0.0;
                for (size_t t = 0; t < len; ++t) {
                    const double res = seg_ptr[t] - mean;
                    var_seg += res * res;
                }
                sum_f2_total += (var_seg / len_d);
            }
            ++conteo_segmentos;
        });

        if (conteo_segmentos == 0) return 0.0;
        const double f2_media = sum_f2_total / static_cast<double>(conteo_segmentos);
        return std::sqrt(std::max(0.0, f2_media));
    }
};

// =======================================================================
// 3. Alineación Espectral de Welch Causal Anclada al Presente (Axioma 6)
// =======================================================================

struct ResultadoAlineacionWelch {
    size_t indice_inicio{0};      // Índice donde arranca la partición útil anclada a t
    size_t longitud_util{0};      // Longitud total cubierta por los bloques enteros
    size_t n_bloques{0};          // Conteo de bloques completos
    size_t residuo_pasado{0};     // Muestras recortadas del pasado lejano [0 : residuo_pasado]
    bool es_valido{true};
};

[[nodiscard]] constexpr ResultadoAlineacionWelch alinear_welch_al_presente_cxx(
    size_t N,
    size_t nperseg,
    size_t noverlap
) noexcept {
    if (nperseg > N || noverlap >= nperseg) {
        return {0, 0, 0, 0, false};
    }

    const size_t step = nperseg - noverlap;
    const size_t n_bloques = (N - noverlap) / step;
    if (n_bloques < 1) {
        return {0, 0, 0, 0, false};
    }

    const size_t longitud_util = (n_bloques - 1) * step + nperseg;
    const size_t residuo_pasado = N - longitud_util;

    // La partición útil culmina exactamente en el índice N - 1 (instante presente)
    return {residuo_pasado, longitud_util, n_bloques, residuo_pasado, true};
}

} // namespace ark::ventanas
