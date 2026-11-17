// armonico_lilyth.hpp -- Etapa 4 del Operador de Lilyth: detector armonico por Prony sobre la curva media (PRONY_CURVA_MEDIA, E1.8),
// fase reescalada (PRONY_FASE_REESCALADA, E1.20), sobrepaso por evento (SOBREPASO_POR_EVENTO), contrastes Monte Carlo contra los dos
// nulos (iid y lineal de Onsager), etiqueta ARMONICO_SOSPECHOSO por evento con su tasa en el nulo al lado, y ajuste de la familia de
// cuatro hipotesis por ark::holm (E1.16). C++20.
//
// Autoria: Argon (Auditoria Externa), 2026-11-06, bajo el plano de Code (Opus 5.5) y las decisiones del Comandante.
// Especificacion: PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md (E1.8, E1.13, E1.16, E1.20) y PLANO_ETAPAS_3_A_8.md (seccion 4).
#pragma once

#include "ark/holm.hpp"
#include "ark/lilyth/nulo_lilyth.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <vector>

namespace ark::lilyth {

using ark::metrologia::NaN;

// Resultado del ajuste de Prony de orden 2: eta(tau) = a1 * eta(tau - 1) + a2 * eta(tau - 2) para tau en [2, T_max], con eta(0) = 1.
struct ResultadoProny {
    bool valido = false;
    std::size_t puntos_usados = 0;          // numero de ecuaciones validas en [2, T_max]
    double a1 = NaN;
    double a2 = NaN;
    double discriminante = NaN;             // a1^2 + 4 * a2 (< 0 => raices complejas conjugadas)
    bool raices_complejas = false;          // oscilador subamortiguado (nunca se activa con rebote AR(1) phi < 0, cuya raiz es real negativa)
    bool raiz_real_negativa = false;        // rebote de periodo 2 barras (discriminante >= 0 y raiz dominante negativa)
    double modulo_raiz = NaN;               // sqrt(-a2) en raices complejas, o max(|r1|, |r2|) en raices reales
    double frecuencia_rad = NaN;            // omega = arccos(a1 / (2 * sqrt(-a2))) en (0, pi)
    double periodo = NaN;                   // 2 * pi / omega, SIN recortar (E1.8)
    double min_eta_desde_2 = NaN;           // min_{tau in [2, T_max]} eta(tau)
    std::size_t tau_min_eta = 0;
    double rmse_recurrencia = NaN;
};

// Resultado de la curva en fase reescalada tau_tilde = tau / tau_rest_k (E1.20).
template <std::size_t T_MAX>
struct ResultadoFaseReescalada {
    bool valido = false;
    std::size_t eventos_incluidos = 0;      // eventos con tau_rest finito <= T_max
    std::size_t excluidos_sin_restituir = 0;// E1.20: eventos con tau_rest = inf (no restituidos), contados y declarados
    double mediana_tau_rest = NaN;
    double paso_rejilla = NaN;              // 1 / mediana(tau_rest)
    std::size_t puntos_rejilla = 0;
    std::array<double, T_MAX> tau_tilde{};
    std::array<double, T_MAX> eta_reescalada{};
    std::array<std::size_t, T_MAX> cuenta_rejilla{};
    ResultadoProny prony{};
};

// Contraste de sobrepaso (curva media y por evento) contra los dos nulos (iid y lineal de Onsager).
struct ContrasteSobrepasoLilyth {
    bool valido = false;
    double min_eta_barra_real = NaN;
    double min_eta_barra_real_emparejado_iid = NaN;       // Code (revision D1c): sobre los eventos con nulo sin memoria
    double min_eta_barra_real_emparejado_onsager = NaN;   // Code (revision D1c): sobre los eventos con nulo lineal
    double mediana_min_eta_barra_iid = NaN;
    double mediana_min_eta_barra_onsager = NaN;
    double p_sobrepaso_curva_iid = NaN;     // p 3 de E1.16 (contra nulo sin memoria: hay oscilacion?)
    double p_sobrepaso_curva_onsager = NaN; // p 4 de E1.16 (contra nulo lineal: hay mas que en calma?)
    ContrasteResumenNulo sobrepaso_por_evento_iid{};
    ContrasteResumenNulo sobrepaso_por_evento_onsager{};
    double cuantil_alfa_min_eta_iid = NaN;
    double cuantil_alfa_min_eta_onsager = NaN;
    std::size_t eventos_etiquetados_reales = 0;
    double tasa_etiqueta_reales = NaN;
    double tasa_etiqueta_nulo_iid = NaN;    // tasa de ARMONICO_SOSPECHOSO en el nulo sin memoria (al lado)
    double tasa_etiqueta_nulo_onsager = NaN;// tasa de ARMONICO_SOSPECHOSO en el nulo lineal (al lado)
};

// Resumen completo de una representacion con la familia de 4 hipotesis contrastada por ark::holm (E1.16).
template <std::size_t T_MAX>
struct ResumenArmonicoYHolm {
    bool valido = false;
    ResultadoProny prony_curva_media{};
    ResultadoFaseReescalada<T_MAX> fase_reescalada{};
    ContrasteSobrepasoLilyth sobrepaso{};
    // Los 4 p crudos de E1.16:
    //   [0] = p_Delta_tau_iid
    //   [1] = p_Delta_tau_Onsager
    //   [2] = p_sobrepaso_iid
    //   [3] = p_sobrepaso_Onsager
    std::array<double, 4> p_crudos{NaN, NaN, NaN, NaN};
    std::array<double, 4> p_holm{NaN, NaN, NaN, NaN};
    std::array<bool, 4> rechaza_holm{false, false, false, false};
    bool familia_holm_valida = false;
};

namespace detalle_armonico {

// Ajusta Prony de orden 2 sobre eta(1..H) con eta(0) = eta0 (por omision 1.0):
//   eta(tau) = a1 * eta(tau - 1) + a2 * eta(tau - 2)  para tau en [2, H].
// Nota: `eta_1_a_H` esta indexado en 0..H-1 correspondiendo a tau = 1..H.
inline ResultadoProny ajustar_prony(std::span<const double> eta_1_a_H, std::size_t t_max_ef, double eta0 = 1.0) noexcept {
    ResultadoProny res{};
    const std::size_t H = std::min(eta_1_a_H.size(), t_max_ef);
    if (H < 3) return res;

    // Estadistico de sobrepaso: min de eta(tau) para tau en [2, H]
    for (std::size_t tau = 2; tau <= H; ++tau) {
        const double v = eta_1_a_H[tau - 1];
        if (!std::isfinite(v)) continue;
        if (!std::isfinite(res.min_eta_desde_2) || v < res.min_eta_desde_2) {
            res.min_eta_desde_2 = v;
            res.tau_min_eta = tau;
        }
    }

    auto leer_eta = [&](std::int64_t tau) -> double {
        if (tau == 0) return eta0;
        if (tau < 0 || static_cast<std::size_t>(tau) > H) return 0.0;
        return eta_1_a_H[static_cast<std::size_t>(tau - 1)];
    };

    const std::size_t tau_ini = 2;
    double g11 = 0.0, g12 = 0.0, g22 = 0.0, b1 = 0.0, b2 = 0.0;
    std::size_t m = 0;
    for (std::size_t tau = tau_ini; tau <= H; ++tau) {
        const double y = leer_eta(static_cast<std::int64_t>(tau));
        const double x1 = leer_eta(static_cast<std::int64_t>(tau) - 1);
        const double x2 = leer_eta(static_cast<std::int64_t>(tau) - 2);
        if (!std::isfinite(y) || !std::isfinite(x1) || !std::isfinite(x2)) continue;
        g11 += x1 * x1;
        g12 += x1 * x2;
        g22 += x2 * x2;
        b1 += x1 * y;
        b2 += x2 * y;
        ++m;
    }
    res.puntos_usados = m;
    if (m < 2) return res;

    const double det = g11 * g22 - g12 * g12;
    const double traza = g11 + g22;
    if (!(traza > 0.0) || std::abs(det) <= 1e-24 * traza * traza) {
        return res;
    }

    const double a1 = (b1 * g22 - b2 * g12) / det;
    const double a2 = (g11 * b2 - g12 * b1) / det;
    if (!std::isfinite(a1) || !std::isfinite(a2)) return res;

    res.a1 = a1;
    res.a2 = a2;
    const double disc = a1 * a1 + 4.0 * a2;
    res.discriminante = disc;

    double ss_err = 0.0;
    for (std::size_t tau = tau_ini; tau <= H; ++tau) {
        const double y = leer_eta(static_cast<std::int64_t>(tau));
        const double x1 = leer_eta(static_cast<std::int64_t>(tau) - 1);
        const double x2 = leer_eta(static_cast<std::int64_t>(tau) - 2);
        if (!std::isfinite(y) || !std::isfinite(x1) || !std::isfinite(x2)) continue;
        const double err = y - (a1 * x1 + a2 * x2);
        ss_err += err * err;
    }
    res.rmse_recurrencia = std::sqrt(ss_err / static_cast<double>(m));

    const bool recortar_periodo_corto = false;
    if (disc < 0.0 && a2 < 0.0) {
        res.raices_complejas = true;
        const double r = std::sqrt(-a2);
        res.modulo_raiz = r;
        const double cos_w = std::clamp(a1 / (2.0 * r), -1.0, 1.0);
        const double omega = std::acos(cos_w);
        res.frecuencia_rad = omega;
        if (omega > 0.0 && std::isfinite(omega)) {
            res.periodo = (2.0 * ark::metrologia::est::PI) / omega;
            if (recortar_periodo_corto && res.periodo <= 4.0) {
                res.raices_complejas = false;
                res.periodo = NaN;
            }
        }
    } else if (disc >= 0.0) {
        res.raices_complejas = false;
        const double s_disc = std::sqrt(std::max(0.0, disc));
        const double r1 = 0.5 * (a1 + s_disc);
        const double r2 = 0.5 * (a1 - s_disc);
        res.modulo_raiz = std::max(std::abs(r1), std::abs(r2));
        const double r_dom = (std::abs(r1) >= std::abs(r2)) ? r1 : r2;
        res.raiz_real_negativa = (r_dom < -1e-6);
    }

    res.valido = true;
    return res;
}

// Construye la curva media en fase reescalada tau_tilde = tau / tau_rest_k (E1.20) y ajusta Prony sobre ella.
template <std::size_t T_MAX>
inline ResultadoFaseReescalada<T_MAX> construir_fase_reescalada(
    std::span<const std::array<double, T_MAX>> etas_eventos,
    std::span<const std::size_t> pasos_eventos,
    std::span<const double> taus_rest_eventos,
    std::size_t t_max_ef) {
    ResultadoFaseReescalada<T_MAX> out{};
    out.tau_tilde.fill(NaN);
    out.eta_reescalada.fill(NaN);
    out.cuenta_rejilla.fill(0);

    const std::size_t K = etas_eventos.size();
    if (K == 0 || pasos_eventos.size() != K || taus_rest_eventos.size() != K) return out;

    const bool excluir_no_restituidos = true;
    std::vector<std::size_t> idx_inc;
    std::vector<double> taus_inc;
    for (std::size_t k = 0; k < K; ++k) {
        const double tr = taus_rest_eventos[k];
        const bool rest = std::isfinite(tr) && tr >= 1.0 && tr <= static_cast<double>(t_max_ef);
        if (rest) {
            idx_inc.push_back(k);
            taus_inc.push_back(tr);
        } else {
            ++out.excluidos_sin_restituir;
            if (!excluir_no_restituidos) {
                idx_inc.push_back(k);
                taus_inc.push_back(static_cast<double>(t_max_ef));
            }
        }
    }

    out.eventos_incluidos = idx_inc.size();
    if (idx_inc.empty()) return out;

    out.mediana_tau_rest = ark::metrologia::est::cuantil(taus_inc, 0.5);
    if (!(out.mediana_tau_rest > 0.0) || !std::isfinite(out.mediana_tau_rest)) return out;
    out.paso_rejilla = 1.0 / out.mediana_tau_rest;

    const std::size_t M_grid = std::min<std::size_t>(t_max_ef, T_MAX);
    std::size_t n_validos_grid = 0;

    for (std::size_t m = 1; m <= M_grid; ++m) {
        const double tt = static_cast<double>(m) * out.paso_rejilla;
        out.tau_tilde[m - 1] = tt;
        double suma = 0.0;
        std::size_t cuenta = 0;

        for (std::size_t ii = 0; ii < idx_inc.size(); ++ii) {
            const std::size_t k = idx_inc[ii];
            const double tr = taus_inc[ii];
            const double t_fis = tt * tr;
            const std::size_t pk = std::min<std::size_t>(pasos_eventos[k], T_MAX);
            // Sin extrapolar mas alla de los pasos observados pk del evento k
            if (t_fis > static_cast<double>(pk) + 1e-12) continue;

            const std::size_t t0 = static_cast<std::size_t>(std::floor(t_fis));
            const double alfa = t_fis - static_cast<double>(t0);
            const double y0 = (t0 == 0) ? 1.0 : etas_eventos[k][t0 - 1];
            if (!std::isfinite(y0)) continue;

            double val = y0;
            if (alfa > 1e-12) {
                if (t0 + 1 > pk) continue;
                const double y1 = etas_eventos[k][t0];
                if (!std::isfinite(y1)) continue;
                val = (1.0 - alfa) * y0 + alfa * y1;
            }
            suma += val;
            ++cuenta;
        }
        out.cuenta_rejilla[m - 1] = cuenta;
        if (cuenta > 0) {
            out.eta_reescalada[m - 1] = suma / static_cast<double>(cuenta);
            n_validos_grid = m;
        }
    }

    out.puntos_rejilla = n_validos_grid;
    if (n_validos_grid >= 3) {
        out.prony = ajustar_prony(std::span<const double>(out.eta_reescalada.data(), n_validos_grid), n_validos_grid, 1.0);
        out.valido = out.prony.valido;
    }
    return out;
}

// Calcula el minimo de eta_k(tau) en [tau_ini, min(pasos, t_max_ef)] para un evento.
template <std::size_t T_MAX>
[[nodiscard]] inline double min_eta_en_rango(
    const std::array<double, T_MAX>& eta,
    std::size_t pasos,
    std::size_t tau_ini,
    std::size_t t_max_ef) noexcept {
    const std::size_t H = std::min({pasos, t_max_ef, T_MAX});
    double mn = NaN;
    for (std::size_t tau = std::max<std::size_t>(1, tau_ini); tau <= H; ++tau) {
        const double v = eta[tau - 1];
        if (!std::isfinite(v)) continue;
        if (!std::isfinite(mn) || v < mn) mn = v;
    }
    return mn;
}

// Calcula el minimo de eta_k(tau) en [tau_rest_k, T_max] (o en [2, T_max] si no restituyo) para SOBREPASO_POR_EVENTO.
template <std::size_t T_MAX>
[[nodiscard]] inline double min_eta_post_restitucion(
    const std::array<double, T_MAX>& eta,
    std::size_t pasos,
    double tau_rest,
    std::size_t t_max_ef) noexcept {
    const std::size_t tau_ini = (std::isfinite(tau_rest) && tau_rest >= 1.0)
                                    ? static_cast<std::size_t>(tau_rest)
                                    : 2u;
    return min_eta_en_rango<T_MAX>(eta, pasos, tau_ini, t_max_ef);
}

// Contraste por Monte Carlo del min de eta_barra(tau) en [2, t_max_ef] contra una poblacion de curvas nulas.
template <std::size_t T_MAX>
inline void contraste_min_eta_barra_mc(
    double min_real,
    std::size_t n_eventos_reales,
    std::span<const EventoNuloLilyth<T_MAX>> pool_nulo,
    std::size_t t_max_ef,
    std::size_t b_mc,
    std::uint64_t semilla,
    double& out_mediana_nulo,
    double& out_p_unilateral) {
    out_mediana_nulo = NaN;
    out_p_unilateral = NaN;
    if (!std::isfinite(min_real) || n_eventos_reales == 0 || pool_nulo.empty()) return;

    const std::size_t M = pool_nulo.size();
    const std::size_t H = std::min<std::size_t>(t_max_ef, T_MAX);
    const std::size_t B = std::max<std::size_t>(99, b_mc);
    std::mt19937_64 rng(detalle_nulo::splitmix64(semilla));
    std::uniform_int_distribution<std::size_t> dist(0, M - 1);

    std::vector<double> mins_mc;
    mins_mc.reserve(B);
    std::size_t menores_o_iguales = 0;

    std::vector<double> suma(H, 0.0);
    std::vector<std::size_t> cuenta(H, 0);

    for (std::size_t b = 0; b < B; ++b) {
        std::fill(suma.begin(), suma.end(), 0.0);
        std::fill(cuenta.begin(), cuenta.end(), 0);
        for (std::size_t k = 0; k < n_eventos_reales; ++k) {
            const auto& ev = pool_nulo[dist(rng)];
            const std::size_t pk = std::min<std::size_t>(ev.pasos_observados, H);
            for (std::size_t tau = 2; tau <= pk; ++tau) {
                const double v = ev.eta[tau - 1];
                if (std::isfinite(v)) {
                    suma[tau - 1] += v;
                    ++cuenta[tau - 1];
                }
            }
        }
        double mn_b = NaN;
        for (std::size_t tau = 2; tau <= H; ++tau) {
            if (cuenta[tau - 1] > 0) {
                const double media = suma[tau - 1] / static_cast<double>(cuenta[tau - 1]);
                if (!std::isfinite(mn_b) || media < mn_b) mn_b = media;
            }
        }
        if (std::isfinite(mn_b)) {
            mins_mc.push_back(mn_b);
            if (mn_b <= min_real + 1e-15) ++menores_o_iguales;
        }
    }

    if (mins_mc.empty()) return;
    out_mediana_nulo = ark::metrologia::est::cuantil(mins_mc, 0.5);
    out_p_unilateral = static_cast<double>(1 + menores_o_iguales) / static_cast<double>(1 + mins_mc.size());
}

// Code (revision D1c): contraste por Monte Carlo del min de eta_barra(tau) en [2, t_max_ef], EMPAREJADO: en cada sorteo se toma, para
// cada evento real con nulo, un evento al azar de SU conjunto emparejado (Etapa 3) y se promedian las curvas. p unilateral (sobrepaso =
// min bajo): (1 + #{min_b <= min_real}) / (1 + B).
template <std::size_t T_MAX>
inline void contraste_min_eta_barra_emparejado(
    double min_real,
    const std::vector<const std::vector<EventoNuloLilyth<T_MAX>>*>& sets,
    std::size_t t_max_ef,
    std::size_t b_mc,
    std::uint64_t semilla,
    double& out_mediana_nulo,
    double& out_p_unilateral) {
    out_mediana_nulo = NaN;
    out_p_unilateral = NaN;
    if (!std::isfinite(min_real) || sets.empty()) return;
    const std::size_t H = std::min<std::size_t>(t_max_ef, T_MAX);
    const std::size_t B = std::max<std::size_t>(1, b_mc);
    std::mt19937_64 rng(detalle_nulo::splitmix64(semilla));
    std::vector<double> mins_mc;
    mins_mc.reserve(B);
    std::size_t menores_o_iguales = 0;
    std::vector<double> suma(H, 0.0);
    std::vector<std::size_t> cuenta(H, 0);
    for (std::size_t b = 0; b < B; ++b) {
        std::fill(suma.begin(), suma.end(), 0.0);
        std::fill(cuenta.begin(), cuenta.end(), 0);
        for (const auto* s : sets) {
            if (s == nullptr || s->empty()) continue;
            std::uniform_int_distribution<std::size_t> d(0, s->size() - 1);
            const auto& ev = (*s)[d(rng)];
            const std::size_t pk = std::min<std::size_t>(ev.pasos_observados, H);
            for (std::size_t tau = 2; tau <= pk; ++tau) {
                const double v = ev.eta[tau - 1];
                if (std::isfinite(v)) { suma[tau - 1] += v; ++cuenta[tau - 1]; }
            }
        }
        double mn_b = NaN;
        for (std::size_t tau = 2; tau <= H; ++tau) {
            if (cuenta[tau - 1] == 0) continue;
            const double m = suma[tau - 1] / static_cast<double>(cuenta[tau - 1]);
            if (!std::isfinite(mn_b) || m < mn_b) mn_b = m;
        }
        if (std::isfinite(mn_b)) {
            mins_mc.push_back(mn_b);
            if (mn_b <= min_real) ++menores_o_iguales;
        }
    }
    if (mins_mc.empty()) return;
    out_mediana_nulo = ark::metrologia::est::cuantil(mins_mc, 0.5);
    out_p_unilateral = static_cast<double>(1 + menores_o_iguales) / static_cast<double>(1 + mins_mc.size());
}

// Evalua si una trayectoria individual recibe la etiqueta ARMONICO_SOSPECHOSO:
//   al menos 4 puntos observados, raices complejas en su propio Prony y min_eta_desde_2 < cuantil_alfa_nulo.
template <std::size_t T_MAX>
[[nodiscard]] inline bool evaluar_etiqueta_armonico_evento(
    const std::array<double, T_MAX>& eta,
    std::size_t pasos_observados,
    double min_eta_desde_2,
    std::size_t t_max_ef,
    double cuantil_alfa_nulo) noexcept {
    const std::size_t H = std::min({pasos_observados, t_max_ef, T_MAX});
    if (H < 4 || !std::isfinite(min_eta_desde_2) || !std::isfinite(cuantil_alfa_nulo)) return false;
    const auto pr = ajustar_prony(std::span<const double>(eta.data(), H), H, 1.0);
    return pr.valido && pr.raices_complejas && (min_eta_desde_2 < cuantil_alfa_nulo);
}

// Construye y evalua la familia de 4 valores p con ark::holm (E1.16).
// Un p no valido hace fallar la familia entera (falla dura de holm.hpp), nunca se omite ni se sanea.
template <std::size_t T_MAX>
inline void aplicar_familia_holm_cuatro(
    double p_dt_iid,
    double p_dt_ons,
    double p_sob_iid,
    double p_sob_ons,
    double alfa,
    ResumenArmonicoYHolm<T_MAX>& res) {
    res.p_crudos = {p_dt_iid, p_dt_ons, p_sob_iid, p_sob_ons};
    res.p_holm = {NaN, NaN, NaN, NaN};
    res.rechaza_holm = {false, false, false, false};
    res.familia_holm_valida = false;

    const std::size_t n_familia_holm = 4;
    const bool sanear_nan_antes_de_holm = false;

    std::vector<double> vec_p;
    vec_p.reserve(n_familia_holm);
    for (std::size_t i = 0; i < n_familia_holm; ++i) {
        double p_i = res.p_crudos[i];
        if (sanear_nan_antes_de_holm && !ark::p_valido(p_i)) {
            p_i = 1.0;
        }
        vec_p.push_back(p_i);
    }

    const std::vector<double> h = ark::holm(vec_p);
    bool todos_validos = !h.empty();
    for (std::size_t i = 0; i < h.size() && i < 4; ++i) {
        res.p_holm[i] = h[i];
        if (!ark::p_valido(h[i])) {
            todos_validos = false;
        } else {
            res.rechaza_holm[i] = (h[i] < alfa);
        }
    }
    res.familia_holm_valida = todos_validos;
}

}  // namespace detalle_armonico

}  // namespace ark::lilyth