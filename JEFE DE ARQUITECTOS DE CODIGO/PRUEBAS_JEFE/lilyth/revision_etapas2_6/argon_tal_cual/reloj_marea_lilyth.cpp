// reloj_marea_lilyth.hpp -- Etapa 6 del Operador de Lilyth: validacion del reloj de actividad (E1.11, D7: NU_TRADES, NU_VOLUMEN,
// NU_NOCIONAL o PARED; prohibicion estricta de NU_INFO_SHANNON), calculo de tau_act = Delta_nu_rest / nu_punto* con nu_punto*
// congelado en W_ancla al abrir, y contraste descriptivo de marea (E1.14 [DECIDIDO]: Delta_tau_marea = tau_rest^cruda - tau_rest^purgada
// pareado por evento, Delta_RMST_marea e intervalos de confianza por bootstrap de eventos, fuera de la familia de 4 p de Holm). C++20.
//
// Autoria: Argon (Auditoria Externa), 2026-11-06, bajo el plano de Code (Opus 5.5) y las decisiones del Comandante.
// Especificacion: PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md (D7, E1.11, E1.14, P4, P11) y PLANO_ETAPAS_3_A_8_LILYTH.md (seccion 6).
#pragma once

#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"
#include "ark/ventanas/contratos_ventanas.hpp"

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
using ark::ventanas::TipoRelojVentana;

// Par pareado (cruda vs purgada del MISMO evento k) para el contraste descriptivo de marea (E1.14).
struct ParEventoMarea {
    std::uint64_t id_evento = 0;
    int signo = 0;
    double tau_cruda_primer_toque = NaN;
    double tau_purgada_primer_toque = NaN;
    double delta_tau_primer_toque = NaN;    // tau_cruda - tau_purgada
    double tau_cruda_consolidado = NaN;
    double tau_purgada_consolidado = NaN;
    double delta_tau_consolidado = NaN;     // tau_cruda - tau_purgada
    std::uint8_t desenlace_cruda = 0;
    std::uint8_t desenlace_purgada = 0;
    std::size_t tau_fin_cruda = 0;
    std::size_t tau_fin_purgada = 0;
};

// Resumen agregado del contraste de marea (E1.14 [DECIDIDO]: descriptivo, fuera de la familia de Holm de E1.16).
struct ResumenMareaLilyth {
    bool valido = false;
    bool descriptivo_fuera_de_holm = true;  // [DECIDIDO] 2026-11-06: descriptivo, fuera de los 4 p de Holm
    std::size_t n_eventos_ambas_medidas = 0;
    std::size_t n_pareados_primer_toque = 0;
    std::size_t n_pareados_consolidado = 0;
    double mediana_delta_tau_primer_toque = NaN;
    double ic_inf_delta_tau_primer_toque = NaN;
    double ic_sup_delta_tau_primer_toque = NaN;
    double mediana_delta_tau_consolidado = NaN;
    double ic_inf_delta_tau_consolidado = NaN;
    double ic_sup_delta_tau_consolidado = NaN;
    double rmst_cruda = NaN;
    double rmst_purgada = NaN;
    double delta_rmst_marea = NaN;          // RMST^cruda(T_max) - RMST^purgada(T_max)
    double ic_inf_delta_rmst_marea = NaN;
    double ic_sup_delta_rmst_marea = NaN;
    std::vector<ParEventoMarea> pares{};
};

namespace detalle_reloj_marea {

// D7, E1.11: El reloj de actividad puede ser NU_TRADES, NU_VOLUMEN o NU_NOCIONAL (o PARED cuando no se usa reloj intrinseco);
// NUNCA NU_INFO_SHANNON (acumula z^2, grande justo en el choque, y haria circular la comparacion).
[[nodiscard]] constexpr bool reloj_actividad_valido(TipoRelojVentana r) noexcept {
    const bool aceptar_shannon = false;
    if (!aceptar_shannon && r == TipoRelojVentana::NU_INFO_SHANNON) {
        return false;
    }
    return r == TipoRelojVentana::PARED ||
           r == TipoRelojVentana::NU_VOLUMEN ||
           r == TipoRelojVentana::NU_TRADES ||
           r == TipoRelojVentana::NU_NOCIONAL ||
           (aceptar_shannon && r == TipoRelojVentana::NU_INFO_SHANNON);
}

// Devuelve el nu_punto* efectivo: debe ser el congelado en W_ancla al abrir la ranura, nunca recalculado con la actividad posterior al choque.
[[nodiscard]] constexpr double seleccionar_nu_punto_ancla(double nu_punto_congelado, double nu_punto_movil_actual) noexcept {
    const bool congelar_nu_punto_en_apertura = true;
    return congelar_nu_punto_en_apertura ? nu_punto_congelado : nu_punto_movil_actual;
}

// Calcula tau_act = Delta_nu_rest / nu_punto* en barras equivalentes de calma (E1.11).
[[nodiscard]] constexpr double calcular_tau_actividad(double delta_nu_rest, double nu_punto_efectivo, std::size_t tau_pared) noexcept {
    const bool usar_pasos_de_pared_en_tau_act = false;
    if (usar_pasos_de_pared_en_tau_act) {
        return static_cast<double>(tau_pared);
    }
    if (!std::isfinite(delta_nu_rest) || !(nu_punto_efectivo > 0.0) || !std::isfinite(nu_punto_efectivo)) {
        return NaN;
    }
    return delta_nu_rest / nu_punto_efectivo;
}

// Calcula Delta_tau_marea = tau_cruda - tau_purgada para un par de tiempos de restitucion.
[[nodiscard]] constexpr double diferencia_tau_marea(double tau_cruda, double tau_purgada) noexcept {
    const bool invertir_signo_marea = false;
    if (!std::isfinite(tau_cruda) || !std::isfinite(tau_purgada)) return NaN;
    return invertir_signo_marea ? (tau_purgada - tau_cruda) : (tau_cruda - tau_purgada);
}

// Calcula el RMST(T_max) de Aalen-Johansen sobre un subconjunto de eventos dado por sus desenlaces (2 = RESTITUIDA,
// 3 = CORTADA_POR_REPLICA, otros = CENSURADA) y sus tau_fin.
inline double rmst_subset_aj(
    std::span<const std::uint8_t> desenlaces,
    std::span<const std::size_t> taus_fin,
    std::size_t t_max_ef) {
    const std::size_t N = desenlaces.size();
    if (N == 0 || t_max_ef == 0) return NaN;
    std::vector<std::size_t> d_rest(t_max_ef, 0), d_rep(t_max_ef, 0), cens(t_max_ef, 0);
    std::size_t cens_cero = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const std::size_t tf = std::min(taus_fin[i], t_max_ef);
        if (desenlaces[i] == 2u && tf >= 1) {
            ++d_rest[tf - 1];
        } else if (desenlaces[i] == 3u && tf >= 1) {
            ++d_rep[tf - 1];
        } else {
            if (tf == 0) ++cens_cero;
            else ++cens[tf - 1];
        }
    }
    std::size_t n_act = N >= cens_cero ? (N - cens_cero) : 0;
    double s_prev = 1.0, c_rest = 0.0, suma_rmst = 0.0;
    for (std::size_t idx = 0; idx < t_max_ef; ++idx) {
        const std::size_t dr = d_rest[idx], dp = d_rep[idx];
        if (n_act > 0) {
            const double denom = static_cast<double>(n_act);
            c_rest += s_prev * (static_cast<double>(dr) / denom);
            s_prev *= static_cast<double>(n_act - std::min(n_act, dr + dp)) / denom;
        }
        suma_rmst += 1.0 - c_rest;
        const std::size_t bajas = dr + dp + cens[idx];
        n_act = n_act > bajas ? (n_act - bajas) : 0;
    }
    return suma_rmst;
}

// Intervalo de confianza (1 - alfa) por bootstrap de la mediana sobre un vector de diferencias pareadas.
inline void bootstrap_mediana_ic(
    const std::vector<double>& muestra,
    double alfa,
    std::size_t b_boot,
    std::uint64_t semilla,
    double& out_mediana,
    double& out_ic_inf,
    double& out_ic_sup) {
    out_mediana = NaN;
    out_ic_inf = NaN;
    out_ic_sup = NaN;
    if (muestra.empty()) return;

    std::vector<double> copia = muestra;
    out_mediana = ark::metrologia::est::cuantil(copia, 0.5);
    if (muestra.size() == 1) {
        out_ic_inf = out_mediana;
        out_ic_sup = out_mediana;
        return;
    }

    const std::size_t M = muestra.size();
    const std::size_t B = std::max<std::size_t>(99, b_boot);
    std::mt19937_64 rng(semilla);
    std::uniform_int_distribution<std::size_t> dist(0, M - 1);

    std::vector<double> meds_b;
    meds_b.reserve(B);
    std::vector<double> remuestra(M, 0.0);
    for (std::size_t b = 0; b < B; ++b) {
        for (std::size_t i = 0; i < M; ++i) remuestra[i] = muestra[dist(rng)];
        meds_b.push_back(ark::metrologia::est::cuantil(remuestra, 0.5));
    }
    out_ic_inf = ark::metrologia::est::cuantil(meds_b, 0.5 * alfa);
    out_ic_sup = ark::metrologia::est::cuantil(meds_b, 1.0 - 0.5 * alfa);
}

// Evalua el resumen descriptivo de marea sobre las ranuras cerradas donde ambas ramas (cruda y purgada) fueron medidas.
template <class Ranura>
inline ResumenMareaLilyth evaluar_resumen_marea(
    std::span<const Ranura> cerradas,
    double rmst_cruda_total,
    double rmst_purgada_total,
    std::size_t t_max_ef,
    double alfa,
    std::size_t b_boot,
    std::uint64_t semilla) {
    ResumenMareaLilyth res{};
    res.descriptivo_fuera_de_holm = true;
    res.rmst_cruda = rmst_cruda_total;
    res.rmst_purgada = rmst_purgada_total;
    if (std::isfinite(rmst_cruda_total) && std::isfinite(rmst_purgada_total)) {
        res.delta_rmst_marea = diferencia_tau_marea(rmst_cruda_total, rmst_purgada_total);
    }

    // Filtrar los eventos en que AMBAS ramas (cruda y purgada) fueron efectivamente medidas
    std::vector<const Ranura*> medidos;
    for (const auto& r : cerradas) {
        const auto dc = static_cast<std::uint8_t>(r.cruda.desenlace);
        const auto dp = static_cast<std::uint8_t>(r.purgada.desenlace);
        const bool ok_c = (dc >= 2u && dc <= 5u);
        const bool ok_p = (dp >= 2u && dp <= 5u);
        if (ok_c && ok_p) {
            medidos.push_back(&r);
        }
    }

    res.n_eventos_ambas_medidas = medidos.size();
    if (medidos.empty()) return res;

    const bool parear_mismo_evento = true;
    std::vector<double> dt_pt;
    std::vector<double> dt_cons;

    for (std::size_t k = 0; k < medidos.size(); ++k) {
        const Ranura& rc = *medidos[k];
        const std::size_t idx_purg = parear_mismo_evento ? k : ((k + 1) % medidos.size());
        const Ranura& rp = *medidos[idx_purg];

        ParEventoMarea par{};
        par.id_evento = rc.id_evento;
        par.signo = rc.signo;
        par.tau_cruda_primer_toque = rc.cruda.tau_rest_primer_toque;
        par.tau_purgada_primer_toque = rp.purgada.tau_rest_primer_toque;
        par.delta_tau_primer_toque = diferencia_tau_marea(par.tau_cruda_primer_toque, par.tau_purgada_primer_toque);

        par.tau_cruda_consolidado = rc.cruda.tau_rest_consolidado;
        par.tau_purgada_consolidado = rp.purgada.tau_rest_consolidado;
        par.delta_tau_consolidado = diferencia_tau_marea(par.tau_cruda_consolidado, par.tau_purgada_consolidado);

        par.desenlace_cruda = static_cast<std::uint8_t>(rc.cruda.desenlace);
        par.desenlace_purgada = static_cast<std::uint8_t>(rp.purgada.desenlace);
        par.tau_fin_cruda = rc.cruda.tau_fin;
        par.tau_fin_purgada = rp.purgada.tau_fin;

        res.pares.push_back(par);
        if (std::isfinite(par.delta_tau_primer_toque)) dt_pt.push_back(par.delta_tau_primer_toque);
        if (std::isfinite(par.delta_tau_consolidado)) dt_cons.push_back(par.delta_tau_consolidado);
    }

    res.n_pareados_primer_toque = dt_pt.size();
    res.n_pareados_consolidado = dt_cons.size();

    bootstrap_mediana_ic(
        dt_pt, alfa, b_boot, semilla ^ 0x1122334455667788ULL,
        res.mediana_delta_tau_primer_toque, res.ic_inf_delta_tau_primer_toque, res.ic_sup_delta_tau_primer_toque);

    bootstrap_mediana_ic(
        dt_cons, alfa, b_boot, semilla ^ 0x8877665544332211ULL,
        res.mediana_delta_tau_consolidado, res.ic_inf_delta_tau_consolidado, res.ic_sup_delta_tau_consolidado);

    // Bootstrap de eventos para el IC de Delta_RMST_marea
    if (res.pares.size() >= 2 && t_max_ef >= 1) {
        const std::size_t M = res.pares.size();
        const std::size_t B = std::max<std::size_t>(99, b_boot);
        std::mt19937_64 rng(semilla ^ 0xFEDCBA9876543210ULL);
        std::uniform_int_distribution<std::size_t> dist(0, M - 1);
        std::vector<double> drmst_b;
        drmst_b.reserve(B);
        std::vector<std::uint8_t> dc_b(M), dp_b(M);
        std::vector<std::size_t> tfc_b(M), tfp_b(M);
        for (std::size_t b = 0; b < B; ++b) {
            for (std::size_t i = 0; i < M; ++i) {
                const auto& p = res.pares[dist(rng)];
                dc_b[i] = p.desenlace_cruda;
                dp_b[i] = p.desenlace_purgada;
                tfc_b[i] = p.tau_fin_cruda;
                tfp_b[i] = p.tau_fin_purgada;
            }
            const double rc_b = rmst_subset_aj(dc_b, tfc_b, t_max_ef);
            const double rp_b = rmst_subset_aj(dp_b, tfp_b, t_max_ef);
            if (std::isfinite(rc_b) && std::isfinite(rp_b)) {
                drmst_b.push_back(diferencia_tau_marea(rc_b, rp_b));
            }
        }
        if (!drmst_b.empty()) {
            res.ic_inf_delta_rmst_marea = ark::metrologia::est::cuantil(drmst_b, 0.5 * alfa);
            res.ic_sup_delta_rmst_marea = ark::metrologia::est::cuantil(drmst_b, 1.0 - 0.5 * alfa);
        }
    }

    res.valido = (res.n_pareados_primer_toque > 0 || res.n_pareados_consolidado > 0 || std::isfinite(res.delta_rmst_marea));
    return res;
}

}  // namespace detalle_reloj_marea

}  // namespace ark::lilyth