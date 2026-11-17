// replicas_lilyth.hpp -- proceso de replicas del nulo L+H de Lilyth (E1.24). C++20.
// Code (Jefe de Arquitectos de Codigo, Opus 5.5), 2026-11-11, por la directriz MENSAJE_ARGON_NULO_CON_REPLICAS.md (H.1-H.4 del
// Comandante, 2026-11-08).
//
// Que hace:
//   - ajusta un Hawkes de DOS TIPOS (alarma de signo + y de signo -) sobre la secuencia de alarmas del operador delta congelado
//     reproducido sobre la historia de calibracion (H.2, H.3: la escala es una marca, no un proceso aparte);
//   - tasa base por bin de fase (no constante: Filimonov y Sornette 2015);
//   - contrasta el ajuste contra su nulo (Poisson por fase, razon de verosimilitud con bootstrap parametrico) y lo declara
//     NO DECIDIBLE si no lo supera (H.4), y NO IDENTIFICADO si la ramificacion o su cota superior llegan a 1;
//   - genera los cumulos de hijas que el nulo L+H inyecta tras cada evento.
//
// Modelo, en tiempo DISCRETO (pasos de la barra base), que es como llegan las alarmas:
//   recuento de alarmas de tipo s en el paso t ~ Poisson(lambda_s(t)),
//   lambda_s(t) = mu_{s, b(t)} + sum_{k : t_k < t} sum_j a_{s_k, s, j} q_j^(t - t_k - 1),   q_j = exp(-1 / tau_j).
//   Masa de ramificacion madre s' -> hija s:  N_{s', s} = sum_j a_{s', s, j} / (1 - q_j).   Decide el radio espectral rho(N).
//   La log-verosimilitud es concava en (mu, a) con las tau fijas; se maximiza por EM (monotona).
//
// ENMIENDA DECLARADA (Code, 2026-11-11) sobre la directriz: las tau_j no son las de las escalas CALIBRADAS del operador delta (con una
// sola escala calibrada habria una sola tau de 1 paso), sino la escalera diadica del operador delta extendida hasta el horizonte en que
// Lilyth mide: tau_k = 2^k, k = 0..K, con 2^K >= T_max. Ningun numero escrito.

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <utility>
#include <vector>

#include "ark/metrologia/calibrador_familia.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

namespace ark::lilyth {

using ark::metrologia::NaN;

// Topes de memoria declarados (Regla 5): escalera de hasta 2^(J_MAX - 1) pasos, iteraciones de EM, eventos por simulacion.
inline constexpr std::size_t REPLICAS_J_MAX = 16;
inline constexpr std::size_t REPLICAS_ITER_MAX = 20000;
inline constexpr std::size_t REPLICAS_HIJAS_MAX = 4096;   // hijas por cumulo en el nulo (desborde contado)

enum class EstadoReplicas : std::uint8_t {
    SIN_AJUSTAR = 0,
    DECIDIBLE = 1,
    NO_DECIDIBLE = 2,        // H.4: no supera su nulo, o no hay alarmas suficientes para intentarlo
    NO_IDENTIFICADO = 3      // rho(N) >= 1 o su cota superior >= 1
};

// Una alarma del operador delta (la apertura del espejo de Lilyth): paso, tipo (0 = signo +, 1 = signo -), escala (marca) y bin de fase.
struct AlarmaReplica {
    std::int64_t t = 0;
    std::uint8_t tipo = 0;
    std::uint8_t escala = 0;
    std::uint32_t bin = 0;
};

// Parametros de un Hawkes de dos tipos con la escalera fija. mu[s * bins + b]; a[madre][hija][j].
struct AjusteHawkesReplicas {
    bool valido = false;
    std::size_t J = 0;
    std::array<double, REPLICAS_J_MAX> tau{};
    std::size_t bins = 1;
    std::vector<double> mu{};
    std::array<std::array<std::array<double, REPLICAS_J_MAX>, 2>, 2> a{};
    std::array<std::array<double, 2>, 2> N{};      // masa de ramificacion madre -> hija
    double rho = NaN;                               // radio espectral de N
    double ll = NaN;                                // log-verosimilitud (sin el termino constante de log n!)
    std::size_t iteraciones = 0;
    bool tope_iteraciones = false;
};

// Ficha del ajuste de replicas de una calibracion.
struct FichaReplicas {
    EstadoReplicas estado = EstadoReplicas::SIN_AJUSTAR;
    const char* motivo = "";
    std::size_t n_alarmas = 0, n_tipo[2]{0, 0};
    std::size_t n_pasos = 0;
    std::size_t bins = 1;
    std::size_t m_min = 0;
    AjusteHawkesReplicas hawkes{};
    AjusteHawkesReplicas poisson{};                 // nulo del ajuste: misma base por fase, a = 0
    double lr = NaN;                                // 2 (LL_Hawkes - LL_Poisson)
    double p_lr = NaN;                              // (1 + #{LR* >= LR}) / (B + 1) bajo el Poisson por fase ajustado
    std::size_t b_boot = 0;                         // B = ceil(m_min / alfa) - 1: m_min sorteos esperados mas alla del corte
    double rho = NaN, rho_sup = NaN;                // rho(N) y su cota superior basica 2 rho - q_alfa(rho*) del bootstrap parametrico
    std::size_t boot_explosivos = 0;                // sorteos del bootstrap de rho que desbordaron el tope (cuentan como rho = inf)
    // Excitacion separada (H.2): mismo signo (diagonal de N) y signo contrario (fuera de la diagonal), con su intervalo bootstrap.
    double n_mismo_signo = NaN, n_signo_contrario = NaN;
    double n_mismo_signo_inf = NaN, n_mismo_signo_sup = NaN;
    double n_contrario_inf = NaN, n_contrario_sup = NaN;
    double rho_inf = NaN;                            // cota inferior basica 2 rho - q_{1-alfa}(rho*) (intervalo de los controles)
    // Inferencia indirecta (opcion (b), ajustar_replicas_indirecto en nulo_lilyth.hpp): `hawkes` es el ESTRUCTURAL (las hijas que
    // disparan las alarmas de la simulacion) y `auxiliar` el Hawkes de las alarmas reales, al que la simulacion se iguala.
    AjusteHawkesReplicas auxiliar{};
    std::size_t iter_indirecta = 0;
    bool convergio = false;
    double se_max_ii = NaN;                          // mayor sd de UNA realizacion de las masas en la ultima iteracion
    std::size_t sims_explosivas_ii = 0;              // simulaciones del punto fijo que desbordaron (se amortigua el nucleo)
    double sd_masa_real[2][2] = {{NaN, NaN}, {NaN, NaN}}; // sd de cada masa del auxiliar real (bootstrap parametrico de si mismo)
    double pendiente_ligadura = NaN;                 // d rho_auxiliar / d rho_estructural (secante, numeros aleatorios comunes)
};

namespace detalle_replicas {

[[nodiscard]] inline std::uint64_t semilla(std::uint64_t base, std::uint64_t i) noexcept {
    std::uint64_t x = base ^ (0x9E3779B97F4A7C15ULL * (i + 1));
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Escalera diadica hasta el horizonte: tau_k = 2^k, k = 0..K, 2^K >= t_max (enmienda declarada arriba).
[[nodiscard]] inline std::size_t escalera_tau(std::size_t t_max, std::array<double, REPLICAS_J_MAX>& tau) noexcept {
    std::size_t K = 0;
    while ((std::size_t{1} << K) < std::max<std::size_t>(t_max, 1) && K + 1 < REPLICAS_J_MAX) ++K;
    const bool escalera_derivada_del_horizonte = true;   // mutante M8: escalera escrita a mano
    if (!escalera_derivada_del_horizonte) {
        tau.fill(0.0);
        tau[0] = 1.0; tau[1] = 10.0; tau[2] = 100.0;
        return 3;
    }
    for (std::size_t k = 0; k <= K; ++k) tau[k] = static_cast<double>(std::size_t{1} << k);
    return K + 1;
}

[[nodiscard]] inline double radio_espectral_2x2(const std::array<std::array<double, 2>, 2>& N) noexcept {
    const double tr = N[0][0] + N[1][1];
    const double det = N[0][0] * N[1][1] - N[0][1] * N[1][0];
    const double disc = std::max(0.0, tr * tr - 4.0 * det);
    return 0.5 * (tr + std::sqrt(disc));
}

inline void calcular_masas(AjusteHawkesReplicas& h) noexcept {
    for (std::size_t m = 0; m < 2; ++m) {
        for (std::size_t s = 0; s < 2; ++s) {
            double acc = 0.0;
            for (std::size_t j = 0; j < h.J; ++j) acc += h.a[m][s][j] / (1.0 - std::exp(-1.0 / h.tau[j]));
            h.N[m][s] = acc;
        }
    }
    h.rho = radio_espectral_2x2(h.N);
}

// Sufijos fijos del ajuste (no dependen de los parametros): R[i][m][j] = sum_{k tipo m, t_k < t_i} q_j^(t_i - t_k - 1) y la masa
// del compensador C[m][j] = sum_{k tipo m} (1 - q_j^(n - 1 - t_k)) / (1 - q_j); T[b] = pasos del bin b; recuentos por (tipo, bin).
struct DatosAjuste {
    std::size_t n = 0;                 // eventos
    std::size_t J = 0, bins = 1;
    std::vector<std::uint8_t> tipo{};
    std::vector<std::uint32_t> bin{};
    std::vector<double> R{};           // [i * 2J + m * J + j]
    std::array<std::array<double, REPLICAS_J_MAX>, 2> C{};
    std::vector<double> T{};           // pasos por bin
    std::vector<double> cuenta{};      // [s * bins + b]
};

// `alarmas` ordenadas por t. `pasos_por_bin` tiene `bins` entradas. n_pasos = longitud de la historia.
inline DatosAjuste preparar_datos(std::span<const AlarmaReplica> alarmas, std::size_t n_pasos, std::span<const double> pasos_por_bin,
                                  std::size_t J, const std::array<double, REPLICAS_J_MAX>& tau, bool fundir_tipos) {
    DatosAjuste d{};
    d.n = alarmas.size();
    d.J = J;
    d.bins = pasos_por_bin.size();
    d.T.assign(pasos_por_bin.begin(), pasos_por_bin.end());
    d.cuenta.assign(2 * d.bins, 0.0);
    d.tipo.resize(d.n);
    d.bin.resize(d.n);
    d.R.assign(d.n * 2 * J, 0.0);
    std::array<double, REPLICAS_J_MAX> q{};
    for (std::size_t j = 0; j < J; ++j) q[j] = std::exp(-1.0 / tau[j]);
    // V[m][j] = sum_{k tipo m, t_k < t_cur} q^(t_cur - t_k - 1), con t_cur el paso del evento en curso; se avanza por grupos de t.
    std::array<std::array<double, REPLICAS_J_MAX>, 2> V{};
    std::array<std::array<double, REPLICAS_J_MAX>, 2> pend{};   // eventos del paso anterior aun no incorporados
    std::int64_t t_cur = alarmas.empty() ? 0 : alarmas[0].t;
    std::size_t i = 0;
    while (i < d.n) {
        const std::int64_t t = alarmas[i].t;
        if (t > t_cur) {
            const double dt = static_cast<double>(t - t_cur);
            for (std::size_t m = 0; m < 2; ++m) {
                for (std::size_t j = 0; j < J; ++j) {
                    // los ya incorporados decaen dt pasos; los del paso t_cur entran con q^(t - t_cur - 1)
                    V[m][j] = V[m][j] * std::pow(q[j], dt) + pend[m][j] * std::pow(q[j], dt - 1.0);
                    pend[m][j] = 0.0;
                }
            }
            t_cur = t;
        }
        // todos los eventos del paso t ven V (solo lo estrictamente anterior)
        std::size_t k = i;
        while (k < d.n && alarmas[k].t == t) {
            const std::uint8_t s = fundir_tipos ? 0 : static_cast<std::uint8_t>(alarmas[k].tipo & 1u);
            d.tipo[k] = s;
            d.bin[k] = alarmas[k].bin;
            for (std::size_t m = 0; m < 2; ++m)
                for (std::size_t j = 0; j < J; ++j) d.R[k * 2 * J + m * J + j] = V[m][j];
            d.cuenta[s * d.bins + d.bin[k]] += 1.0;
            ++k;
        }
        for (std::size_t e = i; e < k; ++e) {
            for (std::size_t j = 0; j < J; ++j) pend[d.tipo[e]][j] += 1.0;
            const double D = static_cast<double>(static_cast<std::int64_t>(n_pasos) - 1 - alarmas[e].t);
            for (std::size_t j = 0; j < J; ++j) d.C[d.tipo[e]][j] += (1.0 - std::pow(q[j], std::max(0.0, D))) / (1.0 - q[j]);
        }
        i = k;
    }
    return d;
}

[[nodiscard]] inline double log_verosimilitud(const DatosAjuste& d, const AjusteHawkesReplicas& h) noexcept {
    double ll = 0.0;
    for (std::size_t i = 0; i < d.n; ++i) {
        const std::size_t s = d.tipo[i];
        double lam = h.mu[s * d.bins + d.bin[i]];
        for (std::size_t m = 0; m < 2; ++m)
            for (std::size_t j = 0; j < d.J; ++j) lam += h.a[m][s][j] * d.R[i * 2 * d.J + m * d.J + j];
        if (!(lam > 0.0)) return -std::numeric_limits<double>::infinity();
        ll += std::log(lam);
    }
    for (std::size_t s = 0; s < 2; ++s)
        for (std::size_t b = 0; b < d.bins; ++b) ll -= h.mu[s * d.bins + b] * d.T[b];
    for (std::size_t m = 0; m < 2; ++m)
        for (std::size_t s = 0; s < 2; ++s)
            for (std::size_t j = 0; j < d.J; ++j) ll -= h.a[m][s][j] * d.C[m][j];
    return ll;
}

// Poisson por fase: mu = recuento / pasos, a = 0 (el nulo del ajuste).
inline AjusteHawkesReplicas ajustar_poisson(const DatosAjuste& d, const std::array<double, REPLICAS_J_MAX>& tau) {
    AjusteHawkesReplicas h{};
    h.J = d.J; h.tau = tau; h.bins = d.bins;
    h.mu.assign(2 * d.bins, 0.0);
    for (std::size_t s = 0; s < 2; ++s)
        for (std::size_t b = 0; b < d.bins; ++b) h.mu[s * d.bins + b] = d.T[b] > 0.0 ? d.cuenta[s * d.bins + b] / d.T[b] : 0.0;
    calcular_masas(h);
    h.ll = log_verosimilitud(d, h);
    h.valido = std::isfinite(h.ll);
    return h;
}

// EM sobre (mu, a) con las tau fijas. Arranca de la mitad de la base de Poisson y una excitacion que explica la otra mitad.
inline AjusteHawkesReplicas ajustar_hawkes(const DatosAjuste& d, const std::array<double, REPLICAS_J_MAX>& tau) {
    AjusteHawkesReplicas h{};
    h.J = d.J; h.tau = tau; h.bins = d.bins;
    h.mu.assign(2 * d.bins, 0.0);
    for (std::size_t s = 0; s < 2; ++s)
        for (std::size_t b = 0; b < d.bins; ++b) h.mu[s * d.bins + b] = d.T[b] > 0.0 ? 0.5 * d.cuenta[s * d.bins + b] / d.T[b] : 0.0;
    std::array<double, 2> ns{0.0, 0.0};
    for (std::size_t i = 0; i < d.n; ++i) ns[d.tipo[i]] += 1.0;
    for (std::size_t m = 0; m < 2; ++m)
        for (std::size_t s = 0; s < 2; ++s)
            for (std::size_t j = 0; j < d.J; ++j)
                h.a[m][s][j] = (d.C[m][j] > 0.0 && ns[s] > 0.0) ? 0.5 * ns[s] / (2.0 * static_cast<double>(d.J) * d.C[m][j]) : 0.0;
    if (d.n == 0) { h.ll = 0.0; calcular_masas(h); h.valido = true; return h; }

    std::vector<double> num_mu(2 * d.bins);
    std::array<std::array<std::array<double, REPLICAS_J_MAX>, 2>, 2> num_a{};
    double ll_prev = log_verosimilitud(d, h);
    // tolerancia numerica relativa (doble precision), no un parametro del modelo
    const double tol = 1e-12;
    for (std::size_t it = 1; it <= REPLICAS_ITER_MAX; ++it) {
        std::fill(num_mu.begin(), num_mu.end(), 0.0);
        for (auto& x : num_a) for (auto& y : x) y.fill(0.0);
        for (std::size_t i = 0; i < d.n; ++i) {
            const std::size_t s = d.tipo[i];
            const double mu_i = h.mu[s * d.bins + d.bin[i]];
            double lam = mu_i;
            for (std::size_t m = 0; m < 2; ++m)
                for (std::size_t j = 0; j < d.J; ++j) lam += h.a[m][s][j] * d.R[i * 2 * d.J + m * d.J + j];
            if (!(lam > 0.0)) continue;
            num_mu[s * d.bins + d.bin[i]] += mu_i / lam;
            for (std::size_t m = 0; m < 2; ++m)
                for (std::size_t j = 0; j < d.J; ++j) num_a[m][s][j] += h.a[m][s][j] * d.R[i * 2 * d.J + m * d.J + j] / lam;
        }
        for (std::size_t s = 0; s < 2; ++s)
            for (std::size_t b = 0; b < d.bins; ++b) h.mu[s * d.bins + b] = d.T[b] > 0.0 ? num_mu[s * d.bins + b] / d.T[b] : 0.0;
        for (std::size_t m = 0; m < 2; ++m)
            for (std::size_t s = 0; s < 2; ++s)
                for (std::size_t j = 0; j < d.J; ++j) h.a[m][s][j] = d.C[m][j] > 0.0 ? num_a[m][s][j] / d.C[m][j] : 0.0;
        const double ll = log_verosimilitud(d, h);
        h.iteraciones = it;
        if (std::abs(ll - ll_prev) <= tol * (1.0 + std::abs(ll))) { ll_prev = ll; break; }
        ll_prev = ll;
        if (it == REPLICAS_ITER_MAX) h.tope_iteraciones = true;
    }
    h.ll = ll_prev;
    calcular_masas(h);
    h.valido = std::isfinite(h.ll);
    return h;
}

// Muestra Poisson(m) que no consume el generador cuando m == 0 (identidad bit a bit con a = 0, control R8).
template <class Rng>
[[nodiscard]] inline std::size_t poisson(double m, Rng& rng) {
    if (!(m > 0.0)) return 0;
    std::poisson_distribution<std::size_t> pd(m);
    return pd(rng);
}

// Retardo de una hija respecto de su madre: 1 + Geom(1 - q), en pasos (el nucleo discreto q^(d - 1)).
template <class Rng>
[[nodiscard]] inline std::int64_t retardo(double q, Rng& rng) {
    std::geometric_distribution<std::int64_t> gd(1.0 - q);
    return 1 + gd(rng);
}

// Hijas de PRIMERA generacion de una alarma de tipo `tipo_madre` (E1.24, opcion (b) del Comandante): (retardo >= 1, tipo). Las hijas
// de las hijas no se generan aqui: en el simulador en linea las dispara el detector cuando ve a cada hija (la cascada pasa por el
// detector, igual que en los datos). Con a = 0 no consume el generador.
template <class Rng>
inline void hijas_directas(const AjusteHawkesReplicas& h, std::uint8_t tipo_madre, Rng& rng,
                           std::vector<std::pair<std::int64_t, std::uint8_t>>& out) {
    out.clear();
    const bool hijas_con_su_tipo = true;          // mutante M4: las hijas heredan el tipo de su madre
    for (std::size_t s = 0; s < 2; ++s) {
        for (std::size_t j = 0; j < h.J; ++j) {
            const double q = std::exp(-1.0 / h.tau[j]);
            const std::size_t k = poisson(h.a[tipo_madre & 1u][s][j] / (1.0 - q), rng);
            for (std::size_t c = 0; c < k; ++c) {
                out.emplace_back(retardo(q, rng), hijas_con_su_tipo ? static_cast<std::uint8_t>(s) : static_cast<std::uint8_t>(tipo_madre & 1u));
            }
        }
    }
}

// Cumulo de descendientes de una madre (t_madre, tipo_madre) hasta el horizonte, en el orden en que se generan (BFS).
// Devuelve false si el cumulo desborda el tope de memoria. Hijas en [t_madre + 1, horizonte).
template <class Rng>
inline bool generar_cumulo(const AjusteHawkesReplicas& h, std::int64_t t_madre, std::uint8_t tipo_madre, std::int64_t horizonte,
                           Rng& rng, std::vector<std::pair<std::int64_t, std::uint8_t>>& hijas) {
    hijas.clear();
    std::vector<std::pair<std::int64_t, std::uint8_t>> cola{{t_madre, tipo_madre}};
    std::size_t cab = 0;
    const bool con_descendencia = true;           // mutante M5: solo la primera generacion
    const bool hijas_con_su_tipo = true;          // mutante M4: las hijas heredan el tipo de su madre
    while (cab < cola.size()) {
        const auto [tp, mp] = cola[cab++];
        for (std::size_t s = 0; s < 2; ++s) {
            for (std::size_t j = 0; j < h.J; ++j) {
                const double q = std::exp(-1.0 / h.tau[j]);
                const std::size_t k = poisson(h.a[mp][s][j] / (1.0 - q), rng);
                for (std::size_t c = 0; c < k; ++c) {
                    const std::int64_t th = tp + retardo(q, rng);
                    if (th >= horizonte) continue;
                    const std::uint8_t tipo_h = hijas_con_su_tipo ? static_cast<std::uint8_t>(s) : mp;
                    hijas.emplace_back(th, tipo_h);
                    if (hijas.size() > REPLICAS_HIJAS_MAX) return false;
                    if (con_descendencia) cola.emplace_back(th, tipo_h);
                }
            }
        }
    }
    return true;
}

// Simula la historia completa (base por fase + cumulos) sobre los mismos pasos y fases. `pasos_de_bin[b]` lista los pasos del bin b.
// Devuelve false si desborda el tope (cuenta como explosivo).
template <class Rng>
inline bool simular_historia(const AjusteHawkesReplicas& h, bool con_cumulos, std::size_t n_pasos,
                             const std::vector<std::vector<std::int64_t>>& pasos_de_bin, std::size_t tope_eventos, Rng& rng,
                             std::vector<AlarmaReplica>& out, const std::vector<std::uint32_t>& bin_de_paso) {
    out.clear();
    std::vector<std::pair<std::int64_t, std::uint8_t>> madres;
    for (std::size_t s = 0; s < 2; ++s) {
        for (std::size_t b = 0; b < h.bins && b < pasos_de_bin.size(); ++b) {
            if (pasos_de_bin[b].empty()) continue;
            const std::size_t k = poisson(h.mu[s * h.bins + b] * static_cast<double>(pasos_de_bin[b].size()), rng);
            std::uniform_int_distribution<std::size_t> u(0, pasos_de_bin[b].size() - 1);
            for (std::size_t c = 0; c < k; ++c) madres.emplace_back(pasos_de_bin[b][u(rng)], static_cast<std::uint8_t>(s));
        }
    }
    std::vector<std::pair<std::int64_t, std::uint8_t>> todos = madres;
    if (con_cumulos) {
        std::vector<std::pair<std::int64_t, std::uint8_t>> hijas;
        for (const auto& [t, s] : madres) {
            if (!generar_cumulo(h, t, s, static_cast<std::int64_t>(n_pasos), rng, hijas)) return false;
            todos.insert(todos.end(), hijas.begin(), hijas.end());
            if (todos.size() > tope_eventos) return false;
        }
    }
    std::sort(todos.begin(), todos.end());
    out.reserve(todos.size());
    for (const auto& [t, s] : todos) {
        AlarmaReplica a{};
        a.t = t; a.tipo = s;
        a.bin = bin_de_paso.empty() ? 0u : bin_de_paso[static_cast<std::size_t>(t)];
        out.push_back(a);
    }
    return true;
}

}  // namespace detalle_replicas

// Configuracion del ajuste: alfa del proyecto, horizonte T_max, bins de fase, semilla, hilos.
struct ConfiguracionReplicas {
    double alfa = 0.05;
    std::size_t t_max = 2;
    std::size_t bins = 1;
    std::uint64_t semilla = 20261111;
    std::size_t hilos = 0;
};

// Ajusta el proceso de replicas sobre las alarmas de calibracion y lo contrasta contra su nulo (H.4).
//   alarmas: ordenadas por t, en [0, bin_de_paso.size()); bin_de_paso: el bin de fase de cada paso de la historia de calibracion.
inline FichaReplicas ajustar_replicas(std::span<const AlarmaReplica> alarmas, const std::vector<std::uint32_t>& bin_de_paso,
                                      const ConfiguracionReplicas& cfg) {
    namespace dr = detalle_replicas;
    FichaReplicas f{};
    f.n_pasos = bin_de_paso.size();
    const bool base_por_fase = true;               // mutante M2: base constante
    f.bins = base_por_fase ? std::max<std::size_t>(1, cfg.bins) : 1;
    f.m_min = ark::metrologia::familia::m_min(cfg.alfa);
    f.n_alarmas = alarmas.size();
    for (const auto& a : alarmas) ++f.n_tipo[a.tipo & 1u];
    std::vector<std::uint32_t> bins_paso(bin_de_paso.size());
    for (std::size_t t = 0; t < bin_de_paso.size(); ++t) {
        if (bin_de_paso[t] >= std::max<std::size_t>(1, cfg.bins)) { f.estado = EstadoReplicas::NO_DECIDIBLE; f.motivo = "fase fuera de rango"; return f; }
        bins_paso[t] = base_por_fase ? bin_de_paso[t] : 0u;
    }
    std::vector<AlarmaReplica> al(alarmas.begin(), alarmas.end());
    for (auto& a : al) {
        if (a.t < 0 || static_cast<std::size_t>(a.t) >= f.n_pasos) { f.estado = EstadoReplicas::NO_DECIDIBLE; f.motivo = "alarma fuera de la historia"; return f; }
        a.bin = bins_paso[static_cast<std::size_t>(a.t)];
    }
    if (f.n_pasos == 0 || f.n_alarmas < f.m_min) {
        f.estado = EstadoReplicas::NO_DECIDIBLE;
        f.motivo = "alarmas insuficientes en la calibracion: el ajuste no puede superar su nulo";
        return f;
    }
    std::array<double, REPLICAS_J_MAX> tau{};
    const std::size_t J = dr::escalera_tau(cfg.t_max, tau);
    std::vector<double> pasos_por_bin(f.bins, 0.0);
    std::vector<std::vector<std::int64_t>> pasos_de_bin(f.bins);
    for (std::size_t t = 0; t < f.n_pasos; ++t) { pasos_por_bin[bins_paso[t]] += 1.0; pasos_de_bin[bins_paso[t]].push_back(static_cast<std::int64_t>(t)); }

    const bool tipos_separados = true;              // mutante M3: un solo tipo que funde los dos signos
    auto ajustar_par = [&](std::span<const AlarmaReplica> ev, AjusteHawkesReplicas& hk, AjusteHawkesReplicas& po) {
        const auto d = dr::preparar_datos(ev, f.n_pasos, pasos_por_bin, J, tau, !tipos_separados);
        hk = dr::ajustar_hawkes(d, tau);
        po = dr::ajustar_poisson(d, tau);
        return 2.0 * std::max(0.0, hk.ll - po.ll);
    };
    f.lr = ajustar_par(al, f.hawkes, f.poisson);
    f.rho = f.hawkes.rho;
    f.n_mismo_signo = std::max(f.hawkes.N[0][0], f.hawkes.N[1][1]);
    f.n_signo_contrario = std::max(f.hawkes.N[0][1], f.hawkes.N[1][0]);
    if (!f.hawkes.valido || !f.poisson.valido) { f.estado = EstadoReplicas::NO_DECIDIBLE; f.motivo = "ajuste no finito"; return f; }

    // B = ceil(m_min / alfa) - 1: en el corte de nivel alfa quedan, en media, m_min sorteos por encima.
    f.b_boot = static_cast<std::size_t>(std::ceil(static_cast<double>(f.m_min) / cfg.alfa)) - 1;
    const std::size_t B = f.b_boot;
    const std::size_t tope_eventos = 64 * (f.n_alarmas + 1);

    // 1. Nulo del ajuste: Poisson por fase ajustado -> LR*.
    std::vector<double> lr_b(B, 0.0);
    ark::metrologia::familia::en_paralelo(B, cfg.hilos, [&](std::size_t b) {
        std::mt19937_64 rng(dr::semilla(cfg.semilla, 2 * b + 1));
        std::vector<AlarmaReplica> sim;
        (void)dr::simular_historia(f.poisson, false, f.n_pasos, pasos_de_bin, tope_eventos, rng, sim, bins_paso);
        AjusteHawkesReplicas hk, po;
        lr_b[b] = ajustar_par(sim, hk, po);
    });
    std::size_t ge = 0;
    for (double v : lr_b) ge += (v >= f.lr);
    f.p_lr = static_cast<double>(1 + ge) / static_cast<double>(B + 1);

    // 2. Cota de rho: bootstrap parametrico bajo el Hawkes ajustado (percentil). Los que desbordan cuentan como rho = inf.
    std::vector<double> rho_b(B, NaN), nd_b(B, NaN), nc_b(B, NaN);
    std::vector<std::uint8_t> expl(B, 0);
    if (f.hawkes.rho < 1.0) {
        ark::metrologia::familia::en_paralelo(B, cfg.hilos, [&](std::size_t b) {
            std::mt19937_64 rng(dr::semilla(cfg.semilla, 2 * b + 2));
            std::vector<AlarmaReplica> sim;
            if (!dr::simular_historia(f.hawkes, true, f.n_pasos, pasos_de_bin, tope_eventos, rng, sim, bins_paso)) {
                expl[b] = 1; rho_b[b] = std::numeric_limits<double>::infinity();
                nd_b[b] = nc_b[b] = std::numeric_limits<double>::infinity();
                return;
            }
            AjusteHawkesReplicas hk, po;
            (void)ajustar_par(sim, hk, po);
            rho_b[b] = hk.rho;
            nd_b[b] = std::max(hk.N[0][0], hk.N[1][1]);
            nc_b[b] = std::max(hk.N[0][1], hk.N[1][0]);
        });
        for (auto e : expl) f.boot_explosivos += e;
        auto cuantil = [](std::vector<double> v, double p) {
            std::sort(v.begin(), v.end());
            const double pos = p * static_cast<double>(v.size() - 1);
            const std::size_t i = static_cast<std::size_t>(std::floor(pos));
            if (i + 1 >= v.size()) return v.back();
            if (!std::isfinite(v[i + 1])) return v[i + 1];
            return v[i] + (pos - static_cast<double>(i)) * (v[i + 1] - v[i]);
        };
        // Cotas BASICAS (2 theta - q): con una historia finita un proceso critico parece subcritico (el horizonte trunca los cumulos y
        // rho se sesga hacia abajo); el percentil hereda ese sesgo y la cota basica lo corrige (control R5). Con un estimador no
        // negativo y sesgado hacia arriba en el borde, el percentil tampoco puede contener el 0 y el basico si (control R2bu).
        f.rho_sup = 2.0 * f.rho - cuantil(rho_b, cfg.alfa);
        f.rho_inf = 2.0 * f.rho - cuantil(rho_b, 1.0 - cfg.alfa);
        f.n_mismo_signo_inf = 2.0 * f.n_mismo_signo - cuantil(nd_b, 1.0 - cfg.alfa / 2.0);
        f.n_mismo_signo_sup = 2.0 * f.n_mismo_signo - cuantil(nd_b, cfg.alfa / 2.0);
        f.n_contrario_inf = 2.0 * f.n_signo_contrario - cuantil(nc_b, 1.0 - cfg.alfa / 2.0);
        f.n_contrario_sup = 2.0 * f.n_signo_contrario - cuantil(nc_b, cfg.alfa / 2.0);
    } else {
        f.rho_sup = std::numeric_limits<double>::infinity();
    }

    const bool exigir_su_nulo = true;               // mutante M9: L+H usado sin pasar su nulo
    const bool rechazar_ramificacion_critica = true; // mutante M7: rho >= 1 aceptado
    if (exigir_su_nulo && !(f.p_lr <= cfg.alfa)) {
        f.estado = EstadoReplicas::NO_DECIDIBLE;
        f.motivo = "el ajuste de replicas no supera su nulo (Poisson por fase)";
    } else if (rechazar_ramificacion_critica && (!(f.rho < 1.0) || !(f.rho_sup < 1.0))) {
        f.estado = EstadoReplicas::NO_IDENTIFICADO;
        f.motivo = "ramificacion critica: rho(N) o su cota superior >= 1";
    } else {
        f.estado = EstadoReplicas::DECIDIBLE;
        f.motivo = "replicas ajustadas y por encima de su nulo";
    }
    return f;
}

}  // namespace ark::lilyth
