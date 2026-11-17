// sincronizador_multiactivo.hpp -- RED MULTIACTIVO K_M, DESCOMPOSICION DE HELMHOLTZ-HODGE E INFORMACION MUTUA (ARK-DATA-01-CXX)
// Reactor ARK / Cancion del Micelio -- ARQUITECTO DE DATOS (Flash) y correccion del Jefe 2026-11-01
//
// Hodge en el grafo completo K_M (verificado): con B la incidencia arista x nodo, L = B^T B = M I - 1 1^T y su pseudoinversa
// L^+ = (I - 1 1^T / M) / M. La divergencia d = B^T f suma cero, asi que phi = L^+ d = d / M es EXACTO, sin SVD. El gradiente es
// F_grad(i, j) = phi_i - phi_j y el resto F - F_grad es ortogonal a el (<F_grad, F_rot> = 0 salvo redondeo). En K_M con todos
// sus triangulos no hay parte armonica: el resto es rotacional puro.
//
// Lo que corrige frente a la version anterior:
//   - la guarda de cada arista en LOS DOS sentidos y con los intervalos de cierre de ventanas_informacion (antes un solo
//     sentido y un instante NYSE para todo TradFi);
//   - la guarda NO entra en el flujo: las aristas en REZAGO o NO_VERIFICADO siguen dentro de Hodge y se marcan. Con este
//     flujo (r_u(t-1) r_v(t) - r_v(t-1) r_u(t)) un desfase de reloj NO es un gradiente exacto (I168), asi que "red desacoplada
//     de la asincronia horaria" no se puede afirmar; el potencial no se lee como "lider" sin un nulo;
//   - Miller-Madow verdadero: I_MM = I_plugin - (K_xy - K_x - K_y + 1) / (2 N ln 2) con las celdas OCUPADAS (antes suponia
//     ocupadas todas las K_x K_y), con cajas por RANGOS (regla del proyecto: sin cajas de igual ancho, que con colas pesadas
//     meten casi todo en una caja), sin recortar en 0 (el recorte sesga hacia arriba bajo la nula: era > 0 en el 52 % de los
//     sorteos independientes) y con un nulo por rotacion circular para poder leerla;
//   - divisores de la ventana desde s = 4 hasta W/4, como filtros_de_ventana.calcular_divisores_exactos (antes desde 6), y sin
//     los "certificados" literales descarte_causal_pct = 0.0 y lookahead_score = 0.0 (la causalidad de una ventana la certifica
//     el testigo de truncado y perturbacion del futuro, no un campo escrito a mano).
#pragma once

#include "tipos_dataset.hpp"
#include "guardas_causales.hpp"
#include "calendario.hpp"
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <array>
#include <span>
#include <algorithm>
#include <cstring>
#include <limits>
#include <numeric>
#include <utility>

namespace ark::dataset {

// =======================================================================
// 1. Arista de la red K_M
// =======================================================================
struct AristaRedMultiactivo {
    uint8_t u_idx;
    uint8_t v_idx;
    char simbolo_u[16];
    char simbolo_v[16];
    ClasificacionCausalPredictiva clasificacion;      // u(t) -> v(t+1), cierre a cierre
    ClasificacionCausalPredictiva clasificacion_vu;   // v(t) -> u(t+1)
    int64_t margen_ms;                                // del sentido u -> v: inicio_lo(v) - info_hi(u)
    int64_t margen_vu_ms;
    double flujo_raw;
    double flujo_ema;
    double flujo_grad;
    double flujo_rot;
};

// =======================================================================
// 2. Helmholtz-Hodge en K_M (sin memoria dinamica, O(E))
// =======================================================================
template <size_t M>
struct EstadoTopologicoHodge {
    static_assert(M >= 2, "La red multiactivo requiere al menos 2 nodos.");
    static constexpr size_t NUM_ARISTAS = M * (M - 1) / 2;

    std::array<double, M> potenciales_nodales{};
    std::array<double, NUM_ARISTAS> flujo_gradiente{};
    std::array<double, NUM_ARISTAS> flujo_rotacional{};
    double energia_gradiente{0.0};
    double energia_rotacional{0.0};
    double ratio_rotacional{0.0};          // E_rot / (E_grad + E_rot)
    double ortogonalidad_error{0.0};       // <F_grad, F_rot>
    // Fraccion que tendria la parte rotacional si el flujo no tuviera estructura: dim(rot) / dim(total) = (E - (M - 1)) / E.
    static constexpr double RATIO_SIN_ESTRUCTURA = static_cast<double>(NUM_ARISTAS - (M - 1)) / static_cast<double>(NUM_ARISTAS);

    void descomponer(const double* f_in, const std::pair<uint8_t, uint8_t>* edges) noexcept {
        potenciales_nodales.fill(0.0);
        flujo_gradiente.fill(0.0);
        flujo_rotacional.fill(0.0);
        energia_gradiente = energia_rotacional = ratio_rotacional = ortogonalidad_error = 0.0;
        if (!f_in || !edges) return;
        for (size_t e = 0; e < NUM_ARISTAS; ++e) {
            potenciales_nodales[edges[e].first] += f_in[e];
            potenciales_nodales[edges[e].second] -= f_in[e];
        }
        const double inv_m = 1.0 / static_cast<double>(M);
        for (size_t i = 0; i < M; ++i) potenciales_nodales[i] *= inv_m;
        double dot = 0.0;
        for (size_t e = 0; e < NUM_ARISTAS; ++e) {
            const double g = potenciales_nodales[edges[e].first] - potenciales_nodales[edges[e].second];
            const double r = f_in[e] - g;
            flujo_gradiente[e] = g;
            flujo_rotacional[e] = r;
            energia_gradiente += g * g;
            energia_rotacional += r * r;
            dot += g * r;
        }
        ortogonalidad_error = dot;
        const double t = energia_gradiente + energia_rotacional;
        if (t > 1e-300) ratio_rotacional = energia_rotacional / t;
    }
};

// =======================================================================
// 3. Informacion mutua discreta por rangos con Miller-Madow
// =======================================================================
class CorrectorMuestraFinita {
public:
    // Sesgo de Miller-Madow de la informacion mutua (bits) con las celdas OCUPADAS: (K_xy - K_x - K_y + 1) / (2 N ln 2).
    [[nodiscard]] static constexpr double correccion_miller_madow(size_t k_xy, size_t k_x, size_t k_y, size_t n) noexcept {
        if (n == 0) return 0.0;
        constexpr double ln2 = 0.693147180559945309417;
        const double num = static_cast<double>(k_xy) - static_cast<double>(k_x) - static_cast<double>(k_y) + 1.0;
        return num / (2.0 * static_cast<double>(n) * ln2);
    }

    // Caja por rangos con empates promediados: caja = floor(rango_medio * nb / n), rango_medio en [0, n-1]. NaN -> caja nb
    // (se excluye la pareja). `orden` es memoria del llamador (n indices); `cajas` recibe n cajas.
    static void cajas_por_rango(const double* x, size_t n, size_t nb, size_t* orden, uint8_t* cajas) noexcept {
        size_t m = 0;
        for (size_t i = 0; i < n; ++i)
            if (!std::isnan(x[i])) orden[m++] = i; else cajas[i] = static_cast<uint8_t>(nb);
        std::sort(orden, orden + m, [x](size_t a, size_t b) { return x[a] < x[b] || (x[a] == x[b] && a < b); });
        size_t i = 0;
        while (i < m) {
            size_t j = i + 1;
            while (j < m && x[orden[j]] == x[orden[i]]) ++j;
            const double rango = 0.5 * static_cast<double>(i + j - 1);
            size_t c = static_cast<size_t>(std::floor(rango * static_cast<double>(nb) / static_cast<double>(m)));
            if (c >= nb) c = nb - 1;
            for (size_t k = i; k < j; ++k) cajas[orden[k]] = static_cast<uint8_t>(c);
            i = j;
        }
    }

    // I_MM en bits sobre parejas sin NaN; puede ser negativa (por debajo del suelo de muestra finita). n_bins en [2, 16].
    // memoria: 2 n size_t + 2 n uint8_t que pone el llamador (cero heap aqui).
    [[nodiscard]] static double informacion_mutua_mm(const double* x, const double* y, size_t n, size_t n_bins, size_t* mem_orden,
                                                     uint8_t* mem_cajas) noexcept {
        if (!x || !y || n < 10 || n_bins < 2 || !mem_orden || !mem_cajas) return std::numeric_limits<double>::quiet_NaN();
        n_bins = std::min<size_t>(n_bins, 16);
        uint8_t* bx = mem_cajas;
        uint8_t* by = mem_cajas + n;
        cajas_por_rango(x, n, n_bins, mem_orden, bx);
        cajas_por_rango(y, n, n_bins, mem_orden + n, by);
        return mi_desde_cajas(bx, by, n, n_bins, 0);
    }

    // Atajo con memoria en la pila para n <= 4096 (la firma anterior). NaN si n es mayor.
    [[nodiscard]] static double informacion_mutua_mm(const double* x, const double* y, size_t n, size_t n_bins = 5) noexcept {
        if (n > 4096) return std::numeric_limits<double>::quiet_NaN();
        size_t orden[2 * 4096];
        uint8_t cajas[2 * 4096];
        return informacion_mutua_mm(x, y, n, n_bins, orden, cajas);
    }

    // I_MM con y rotada circularmente `k` posiciones (nulo por rotacion: conserva la memoria de cada serie).
    [[nodiscard]] static double mi_desde_cajas(const uint8_t* bx, const uint8_t* by, size_t n, size_t nb, size_t k) noexcept {
        uint32_t cxy[16][16]{}, cx[16]{}, cy[16]{};
        size_t m = 0;
        for (size_t i = 0; i < n; ++i) {
            const uint8_t a = bx[i], b = by[(i + k) % n];
            if (a >= nb || b >= nb) continue;
            ++cxy[a][b]; ++cx[a]; ++cy[b]; ++m;
        }
        if (m < 10) return std::numeric_limits<double>::quiet_NaN();
        constexpr double inv_ln2 = 1.44269504088896340736;
        const double inv_m = 1.0 / static_cast<double>(m);
        double mi = 0.0;
        size_t kxy = 0, kx = 0, ky = 0;
        for (size_t a = 0; a < nb; ++a) {
            if (cx[a]) ++kx;
            if (cy[a]) ++ky;
            for (size_t b = 0; b < nb; ++b) {
                if (!cxy[a][b]) continue;
                ++kxy;
                const double pxy = cxy[a][b] * inv_m;
                mi += pxy * std::log(pxy / ((cx[a] * inv_m) * (cy[b] * inv_m))) * inv_ln2;
            }
        }
        return mi - correccion_miller_madow(kxy, kx, ky, m);
    }

    // p del nulo por rotacion: B desplazamientos repartidos en [n/10, 9n/10]; p = (1 + #{nulo >= obs}) / (B + 1).
    [[nodiscard]] static double p_rotacion(const uint8_t* bx, const uint8_t* by, size_t n, size_t nb, double obs, size_t B = 199) noexcept {
        if (std::isnan(obs) || n < 20 || B == 0) return std::numeric_limits<double>::quiet_NaN();
        const size_t lo = n / 10, hi = n - n / 10;
        size_t ge = 0;
        for (size_t b = 0; b < B; ++b) {
            const size_t k = lo + (hi - lo) * b / B;
            if (mi_desde_cajas(bx, by, n, nb, k) >= obs) ++ge;
        }
        return static_cast<double>(1 + ge) / static_cast<double>(B + 1);
    }
};

// =======================================================================
// 4. Ventana para filtros_de_ventana (ARK-WIN-01): valores, mascara y divisores exactos
// =======================================================================
template <size_t W>
struct VentanaCausalCertificadaCxx {
    static_assert(W >= 12, "La ventana debe tener al menos 12 muestras.");

    uint64_t estampa_t_ms{0};                  // estampa de disponibilidad de la ultima barra de la ventana
    std::array<double, W> valores{};
    std::array<bool, W> es_imputado{};
    size_t num_divisores_exactos{0};
    std::array<size_t, 64> divisores_exactos{};
    bool es_valida{true};

    // Divisores s de W en [4, max(4, W / 4)] (filtros_de_ventana.calcular_divisores_exactos con s_min = 4, max_ratio = 4).
    void inicializar_catalogo_divisores() noexcept {
        num_divisores_exactos = 0;
        const size_t s_max = std::max<size_t>(4, W / 4);
        for (size_t s = 4; s <= s_max && num_divisores_exactos < divisores_exactos.size(); ++s)
            if (W % s == 0) divisores_exactos[num_divisores_exactos++] = s;
    }

    // Residuos de la recta de minimos cuadrados sobre la ventana. Es un AGREGADO de la ventana: el residuo de una muestra
    // interior depende de toda la ventana (tambien de las muestras posteriores a ella), asi que solo el de la ultima es una
    // cantidad del instante t.
    void eliminar_tendencia_lineal(double* out_detrend) const noexcept {
        if (!out_detrend) return;
        const double w = static_cast<double>(W);
        const double sum_t = 0.5 * static_cast<double>((W - 1) * W);
        const double sum_t2 = static_cast<double>((W - 1) * W * (2 * W - 1)) / 6.0;
        double sum_y = 0.0, sum_ty = 0.0;
        for (size_t t = 0; t < W; ++t) { sum_y += valores[t]; sum_ty += static_cast<double>(t) * valores[t]; }
        const double denom = w * sum_t2 - sum_t * sum_t;
        const double b = (denom != 0.0) ? (w * sum_ty - sum_t * sum_y) / denom : 0.0;
        const double a = (sum_y - b * sum_t) / w;
        for (size_t t = 0; t < W; ++t) out_detrend[t] = valores[t] - (a + b * static_cast<double>(t));
    }
};

// =======================================================================
// 5. Sincronizador multiactivo (K_M, sin memoria dinamica)
// =======================================================================
template <size_t M>
class SincronizadorMultiactivoCxx {
public:
    static constexpr size_t NUM_ARISTAS = M * (M - 1) / 2;

private:
    std::array<std::pair<uint8_t, uint8_t>, NUM_ARISTAS> aristas_{};
    std::array<AristaRedMultiactivo, NUM_ARISTAS> telemetria_aristas_{};
    EstadoTopologicoHodge<M> estado_hodge_{};
    char nombres_nodos_[M][16]{};
    Venue venues_nodos_[M]{};

public:
    SincronizadorMultiactivoCxx() noexcept {
        size_t idx = 0;
        for (uint8_t i = 0; i < static_cast<uint8_t>(M); ++i)
            for (uint8_t j = i + 1; j < static_cast<uint8_t>(M); ++j) {
                aristas_[idx] = {i, j};
                auto& a = telemetria_aristas_[idx];
                a = AristaRedMultiactivo{};
                a.u_idx = i;
                a.v_idx = j;
                a.clasificacion = a.clasificacion_vu = ClasificacionCausalPredictiva::NO_VERIFICADO;
                ++idx;
            }
    }

    void configurar_nodo(size_t idx, const char* nombre, Venue v) noexcept {
        if (idx >= M || !nombre) return;
        size_t i = 0;
        while (i < 15 && nombre[i] != '\0') { nombres_nodos_[idx][i] = nombre[i]; ++i; }
        nombres_nodos_[idx][i] = '\0';
        venues_nodos_[idx] = v;
        for (size_t e = 0; e < NUM_ARISTAS; ++e) {
            if (aristas_[e].first == idx) std::memcpy(telemetria_aristas_[e].simbolo_u, nombres_nodos_[idx], 16);
            if (aristas_[e].second == idx) std::memcpy(telemetria_aristas_[e].simbolo_v, nombres_nodos_[idx], 16);
        }
    }

    // Guarda de cada arista en los dos sentidos, cierre a cierre, con las barras de la MISMA fila de la rejilla: X(t) -> Y(t+1)
    // exige que la informacion de X este completa antes de que empiece el retorno de Y, que empieza en el cierre de Y en t.
    void auditar_asincronia_sesion(const BarraCausal* barras_sesion_por_nodo, double latencia_s = LATENCIA_DIARIA_S) noexcept {
        if (!barras_sesion_por_nodo) return;
        for (size_t e = 0; e < NUM_ARISTAS; ++e) {
            const BarraCausal& u = barras_sesion_por_nodo[aristas_[e].first];
            const BarraCausal& v = barras_sesion_por_nodo[aristas_[e].second];
            auto& a = telemetria_aristas_[e];
            const bool ok = u.cierre_hi_utc_ms != 0 && v.cierre_hi_utc_ms != 0;
            const ResultadoGuarda uv = clasificar_ventanas(u.cierre_lo_utc_ms, u.cierre_hi_utc_ms, v.cierre_lo_utc_ms, v.cierre_hi_utc_ms, latencia_s);
            const ResultadoGuarda vu = clasificar_ventanas(v.cierre_lo_utc_ms, v.cierre_hi_utc_ms, u.cierre_lo_utc_ms, u.cierre_hi_utc_ms, latencia_s);
            a.clasificacion = ok ? uv.clase : ClasificacionCausalPredictiva::NO_VERIFICADO;
            a.clasificacion_vu = ok ? vu.clase : ClasificacionCausalPredictiva::NO_VERIFICADO;
            a.margen_ms = static_cast<int64_t>(std::llround(uv.margen_s * 1000.0));
            a.margen_vu_ms = static_cast<int64_t>(std::llround(vu.margen_s * 1000.0));
        }
    }

    // Flujo de adelanto por arista F_uv(t) = r_u(t-1) r_v(t) - r_v(t-1) r_u(t), integrado con una EMA causal, y Hodge.
    // r_mat: T filas x M columnas. Filas con un NaN en u o v no entran en la EMA de esa arista.
    void procesar_flujos_y_hodge(const double* r_mat, size_t T, double alpha_ema = 0.05) noexcept {
        if (!r_mat || T < 3) return;
        std::array<double, NUM_ARISTAS> flujos_ema{};
        for (size_t e = 0; e < NUM_ARISTAS; ++e) {
            const uint8_t u = aristas_[e].first, v = aristas_[e].second;
            double ema = 0.0, ult = 0.0;
            bool primero = true;
            for (size_t t = 1; t < T; ++t) {
                const double f = r_mat[(t - 1) * M + u] * r_mat[t * M + v] - r_mat[(t - 1) * M + v] * r_mat[t * M + u];
                if (std::isnan(f)) continue;
                ema = primero ? f : (1.0 - alpha_ema) * ema + alpha_ema * f;
                primero = false;
                ult = f;
            }
            flujos_ema[e] = ema;
            telemetria_aristas_[e].flujo_raw = ult;
            telemetria_aristas_[e].flujo_ema = ema;
        }
        estado_hodge_.descomponer(flujos_ema.data(), aristas_.data());
        for (size_t e = 0; e < NUM_ARISTAS; ++e) {
            telemetria_aristas_[e].flujo_grad = estado_hodge_.flujo_gradiente[e];
            telemetria_aristas_[e].flujo_rot = estado_hodge_.flujo_rotacional[e];
        }
    }

    [[nodiscard]] const EstadoTopologicoHodge<M>& estado_hodge() const noexcept { return estado_hodge_; }
    [[nodiscard]] const AristaRedMultiactivo& arista(size_t e) const noexcept { return telemetria_aristas_[e]; }
    [[nodiscard]] size_t num_aristas() const noexcept { return NUM_ARISTAS; }
};

}  // namespace ark::dataset
