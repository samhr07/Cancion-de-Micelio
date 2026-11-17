// fin_calor_lilyth.hpp -- Etapa 5 del Operador de Lilyth: inicio del calor desde SalidaTemperatura (signo -1: baja beta = sube
// temperatura en canal_temperatura.hpp) y deteccion secuencial unilateral del fin del calor hacia mu_0,beta calibrada sobre el nulo de
// log beta^ desplazado a mu_calor (E1.7), con mu_calor = media de las l* ventanas posteriores a la alarma [DECIDIDO],
// ARL_0,caliente = T_max / alfa via detalle_nivel::buscar_h_con (rechazando h extrapolado), y regla multiescala (el calor termina
// cuando terminan todas las escalas calientes). Cero reservas de memoria en el ciclo caliente. C++20.
//
// Autoria: Argon (Auditoria Externa), 2026-11-06, bajo el plano de Code (Opus 5.5) y las decisiones del Comandante.
// Especificacion: PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md (E1.6, E1.7, P9) y PLANO_ETAPAS_3_A_8_LILYTH.md (seccion 5).
#pragma once

#include "ark/metrologia/canal_nivel.hpp"
#include "ark/metrologia/canal_temperatura.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
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

using ark::metrologia::NaN;

// Detector GLR secuencial unilateral sobre ventanas de log beta^:
//   Bajo el nulo caliente (media mu_calor < mu_0,beta, desviacion s0, correlacion integrada tau_int),
//   acumula z_w = (log_beta_w - mu_calor) / s0 y evalua el estadistico unilateral de subida hacia mu_0,beta (S_l > 0):
//     G_w^+ = max_{1 <= l <= n} [ max(0, S_l)^2 / (2 * l * tau_int) ].
//   Cero reservas de memoria en paso().
template <std::size_t L = 120>
class DetectorGLRUnilateralCalor {
    static_assert(L >= 2, "L >= 2");
    std::array<double, L> z_{};
    std::size_t cab_ = 0;
    std::size_t n_ = 0;
    double mu_base_ = NaN;
    double s0_ = NaN;
    double tau_int_ = NaN;
    double h_ = NaN;
    std::size_t l_star_base_ = 0; // 0 = mu_base fijo; > 0 = media muestral de l* ventanas (factor de varianza 1 + l / l*)
    bool invertir_hipotesis_f4_ = false;

public:
    struct Paso {
        bool alarma = false;
        double g = NaN;
        std::size_t l = 0;
        double suma = NaN;
    };

    void configurar(double mu_base, double s0, double tau_int, double h, bool invertir_f4 = false, std::size_t l_star_base = 0) noexcept {
        mu_base_ = mu_base;
        s0_ = s0;
        tau_int_ = tau_int;
        h_ = h;
        l_star_base_ = l_star_base;
        invertir_hipotesis_f4_ = invertir_f4;
        reiniciar();
    }

    void reiniciar() noexcept {
        cab_ = 0;
        n_ = 0;
    }

    [[nodiscard]] std::size_t lleno() const noexcept { return n_; }
    [[nodiscard]] double mu_base() const noexcept { return mu_base_; }
    [[nodiscard]] double h() const noexcept { return h_; }

    Paso paso(double log_beta) noexcept {
        Paso p{};
        if (!std::isfinite(log_beta) || !(s0_ > 0.0) || !(tau_int_ > 0.0) || !std::isfinite(mu_base_) || !std::isfinite(h_)) {
            reiniciar();
            return p;
        }
        const double z = (log_beta - mu_base_) / s0_;
        z_[(cab_ + n_) % L] = z;
        if (n_ < L) {
            ++n_;
        } else {
            cab_ = (cab_ + 1) % L;
        }

        if (invertir_hipotesis_f4_) {
            // Defecto F4 (hipotesis nula invertida): evalua contra mu_0,beta y declara fin del calor en cuanto una ventana no cae
            // por debajo de mu_0,beta - s0 (toca la banda de calma por puro ruido estando aun caliente).
            p.g = z;
            p.l = 1;
            p.suma = z;
            if (z >= -1.0) {
                p.alarma = true;
                reiniciar();
            }
            return p;
        }

        double S = 0.0;
        double gmax = 0.0;
        double smax = 0.0;
        std::size_t lmax = 1;
        for (std::size_t l = 1; l <= n_; ++l) {
            S += z_[(cab_ + n_ - l) % L];
            const double factor_mu = (l_star_base_ > 0) ? (1.0 + static_cast<double>(l) / static_cast<double>(l_star_base_)) : 1.0;
            const double g = (S > 0.0) ? (S * S) / (2.0 * static_cast<double>(l) * tau_int_ * factor_mu) : 0.0;
            if (g > gmax) {
                gmax = g;
                lmax = l;
                smax = S;
            }
        }
        p.g = gmax;
        p.l = lmax;
        p.suma = smax;
        if (gmax >= h_) {
            p.alarma = true;
            reiniciar();
        }
        return p;
    }
};

// Ficha de calibracion del detector de fin del calor para una escala j (E1.7).
struct FichaFinCalorEscala {
    int j = 0;
    std::int64_t m_escala = 1;
    std::size_t T_ventana = 0;
    bool temp_detectable = false;
    double mu0_log_beta = NaN;
    double s0_log_beta = NaN;
    double tau_int_log_beta = NaN;
    std::size_t ventanas_nulo = 0;
    double arl0_caliente = NaN;     // T_max / alfa (aproximacion geometrica, declarada)
    double h = NaN;                 // umbral unilateral calibrado con detalle_nivel::buscar_h_con
    bool h_extrapolado = false;     // true si buscar_h_con extrapolo por falta de alarmas legibles
    std::size_t alarmas_nulo = 0;
    bool calibrado = false;         // false si h_extrapolado == true (prohibido usar h extrapolado, seccion 5)
    const char* motivo = "";
};

// Estado en tiempo de ejecucion del episodio de calor de una escala j (cero reservas en ciclo caliente).
template <std::size_t L = 120>
struct EstadoCalorEscala {
    bool caliente = false;
    std::int64_t t_inicio_calor = -1;
    std::int64_t t_ultimo_fin_calor = -1;
    std::size_t l_star = 0;                 // memoria l* de la propia alarma de temperatura [DECIDIDO]
    std::size_t ventanas_mu_calor = 0;      // ventanas posteriores a la alarma ya acumuladas (de 0 a l_star)
    double suma_mu_calor = 0.0;
    bool mu_calor_listo = false;            // true cuando se completaron las l* ventanas posteriores a la alarma
    double mu_calor = NaN;                  // media de log_beta en las l* ventanas posteriores a la alarma
    double amplitud_alarma = NaN;           // guardada en la ficha de la alarma (no se usa para mu_calor salvo mutante M02)
    std::size_t ventanas_tras_mu_calor = 0; // ventanas evaluadas por el GLR unilateral tras cerrar las l* ventanas
    double ultima_g_fin = NaN;
    DetectorGLRUnilateralCalor<L> glr_fin{};

    void limpiar() noexcept {
        caliente = false;
        t_inicio_calor = -1;
        t_ultimo_fin_calor = -1;
        l_star = 0;
        ventanas_mu_calor = 0;
        suma_mu_calor = 0.0;
        mu_calor_listo = false;
        mu_calor = NaN;
        amplitud_alarma = NaN;
        ventanas_tras_mu_calor = 0;
        ultima_g_fin = NaN;
        glr_fin.reiniciar();
    }
};

struct ContadoresFinCalor {
    std::size_t inicios_calor = 0;
    std::size_t fines_calor = 0;
    std::size_t ventanas_calientes = 0;
    std::size_t alarmas_ignoradas_no_calibrado = 0;
};

namespace detalle_fin_calor {

// Cuenta las alarmas del detector unilateral sobre segmentos de log_beta desplazados a mu_calor_ref.
template <std::size_t L>
inline std::size_t alarmas_unilaterales(
    const std::vector<std::vector<double>>& segs_lb_desplazados,
    double mu_calor_ref,
    double s0,
    double tau_int,
    double h,
    std::size_t& total_ventanas) noexcept {
    DetectorGLRUnilateralCalor<L> d;
    std::size_t al = 0;
    total_ventanas = 0;
    for (const auto& g : segs_lb_desplazados) {
        d.configurar(mu_calor_ref, s0, tau_int, h, false);
        for (double v : g) {
            if (!std::isfinite(v)) {
                d.reiniciar();
                continue;
            }
            al += d.paso(v).alarma ? 1u : 0u;
            ++total_ventanas;
        }
    }
    return al;
}

// Construye los segmentos de log_beta^ en ventanas disjuntas de longitud T sobre la historia de la escala j,
// exactamente como CanalTemperatura::calibrar en canal_temperatura.hpp.
inline std::vector<std::vector<double>> extraer_segmentos_log_beta(
    const std::vector<std::vector<double>>& segs_escala,
    std::size_t T) {
    std::vector<std::vector<double>> lb;
    if (T < 2) return lb;
    const auto w = ark::metrologia::detalle_temp::ventanas(segs_escala, T);
    for (const auto& g : w) {
        std::vector<double> v;
        for (const auto& x : g) {
            if (!(x.s2 > 0.0)) continue;
            v.push_back(-std::log(x.s2));
        }
        if (!v.empty()) lb.push_back(std::move(v));
    }
    return lb;
}

// Calibra la FichaFinCalorEscala de la escala j usando detalle_nivel::buscar_h_con con ARL_0,caliente = T_max / alfa.
template <std::size_t L = 120>
inline FichaFinCalorEscala calibrar_escala_fin_calor(
    int j,
    std::int64_t m_escala,
    std::size_t T_ventana,
    const ark::metrologia::FichaEscala& ficha_temp,
    const std::vector<std::vector<double>>& segs_escala,
    std::size_t t_max_efectivo,
    double alfa,
    std::size_t alarmas_minimas = 8) {
    FichaFinCalorEscala f{};
    f.j = j;
    f.m_escala = m_escala;
    f.T_ventana = T_ventana;
    f.temp_detectable = ficha_temp.detectable;
    f.mu0_log_beta = ficha_temp.nulo.mu0;
    f.s0_log_beta = ficha_temp.nulo.s0;
    f.tau_int_log_beta = ficha_temp.nulo.tau_int;

    const bool dividir_arl_por_alfa = true;
    f.arl0_caliente = dividir_arl_por_alfa
                          ? (static_cast<double>(t_max_efectivo) / alfa)
                          : static_cast<double>(t_max_efectivo);

    if (T_ventana < 2 || !ficha_temp.detectable || !(f.s0_log_beta > 0.0) || !(f.tau_int_log_beta > 0.0) ||
        !(f.arl0_caliente > 1.0)) {
        f.motivo = "canal de temperatura no detectable o parametros invalidos";
        return f;
    }

    const auto lb = extraer_segmentos_log_beta(segs_escala, T_ventana);
    if (lb.empty()) {
        f.motivo = "sin ventanas de log_beta en la historia de la escala";
        return f;
    }

    // Desplazar el nulo de log_beta a un nivel caliente de referencia mu_calor_ref = mu0 - 2 * s0 (invariante en z).
    const double mu_calor_ref = f.mu0_log_beta - 2.0 * f.s0_log_beta;
    const double shift = mu_calor_ref - f.mu0_log_beta;
    std::vector<std::vector<double>> lb_desp = lb;
    std::size_t n_tot = 0;
    for (auto& g : lb_desp) {
        for (double& v : g) {
            v += shift;
            ++n_tot;
        }
    }
    f.ventanas_nulo = n_tot;

    std::size_t total_eval = 0;
    auto eval = [&](double h_cand, std::size_t& al) -> double {
        al = alarmas_unilaterales<L>(lb_desp, mu_calor_ref, f.s0_log_beta, f.tau_int_log_beta, h_cand, total_eval);
        return al ? (static_cast<double>(total_eval) / static_cast<double>(al))
                  : std::numeric_limits<double>::infinity();
    };

    const auto rh = ark::metrologia::detalle_nivel::buscar_h_con(eval, f.arl0_caliente, alarmas_minimas);
    f.h = rh.h;
    f.h_extrapolado = rh.extrapolado;
    f.alarmas_nulo = rh.alarmas;

    const bool permitir_h_extrapolado = false;
    if (!std::isfinite(rh.h)) {
        f.calibrado = false;
        f.motivo = "h no finito: historia insuficiente para fin del calor";
    } else if (rh.extrapolado && !permitir_h_extrapolado) {
        f.calibrado = false;
        f.motivo = "h extrapolado: fin del calor declarado no calibrado (no se usa)";
    } else {
        f.calibrado = true;
        f.motivo = "calibrado";
    }
    return f;
}

}  // namespace detalle_fin_calor

// Gestor multiescala de inicio y fin del calor para ObservadorLilyth.
// Declaracion multiescala (seccion 5): cuando hay varias escalas de temperatura calientes a la vez,
// el calor global permanece activo hasta que TODAS las escalas calientes hayan declarado su fin del calor.
template <std::size_t J_MAX = 24, std::size_t L = 120>
class GestorFinCalorLilyth {
    std::size_t escalas_ = 0;
    std::array<FichaFinCalorEscala, J_MAX> fichas_{};
    std::array<EstadoCalorEscala<L>, J_MAX> estados_{};
    ContadoresFinCalor cont_{};

public:
    GestorFinCalorLilyth() = default;

    template <class Obs>
    void calibrar(const Obs& delta, std::size_t t_max_efectivo, double alfa, std::size_t alarmas_minimas = 8) {
        escalas_ = std::min(delta.escalas(), J_MAX);
        for (std::size_t j = 0; j < escalas_; ++j) {
            const auto& ft = delta.ficha_temperatura(j);
            const std::size_t Tv = delta.temperatura(j).T();
            const std::int64_t m_esc = std::int64_t{1} << j;
            const std::size_t paso_ven = std::max<std::size_t>(1, static_cast<std::size_t>(m_esc) * std::max<std::size_t>(1, Tv));
            const std::size_t t_max_ven = std::max<std::size_t>(2, (t_max_efectivo + paso_ven - 1) / paso_ven);
            fichas_[j] = detalle_fin_calor::calibrar_escala_fin_calor<L>(
                static_cast<int>(j), m_esc, Tv, ft, delta.historia(j), t_max_ven, alfa, alarmas_minimas);
        }
        reiniciar_flujo();
    }

    void calibrar_desde_segmentos(
        std::size_t j,
        std::int64_t m_escala,
        std::size_t T_ventana,
        const ark::metrologia::FichaEscala& ficha_temp,
        const std::vector<std::vector<double>>& segs_escala,
        std::size_t t_max_efectivo,
        double alfa,
        std::size_t alarmas_minimas = 8) {
        if (j >= J_MAX) return;
        if (j + 1 > escalas_) escalas_ = j + 1;
        fichas_[j] = detalle_fin_calor::calibrar_escala_fin_calor<L>(
            static_cast<int>(j), m_escala, T_ventana, ficha_temp, segs_escala, t_max_efectivo, alfa, alarmas_minimas);
        estados_[j].limpiar();
    }

    void reiniciar_flujo() noexcept {
        for (auto& e : estados_) e.limpiar();
        cont_ = ContadoresFinCalor{};
    }

    [[nodiscard]] std::size_t escalas() const noexcept { return escalas_; }
    [[nodiscard]] const FichaFinCalorEscala& ficha(std::size_t j) const noexcept { return fichas_[j < J_MAX ? j : 0]; }
    [[nodiscard]] const EstadoCalorEscala<L>& estado(std::size_t j) const noexcept { return estados_[j < J_MAX ? j : 0]; }
    [[nodiscard]] const ContadoresFinCalor& contadores() const noexcept { return cont_; }

    [[nodiscard]] bool calor_activo() const noexcept {
        for (std::size_t j = 0; j < escalas_; ++j) {
            if (estados_[j].caliente) return true;
        }
        return false;
    }

    [[nodiscard]] std::size_t escalas_calientes() const noexcept {
        std::size_t c = 0;
        for (std::size_t j = 0; j < escalas_; ++j) {
            if (estados_[j].caliente) ++c;
        }
        return c;
    }

    // Ciclo caliente: procesa las salidas de temperatura `temp[0..niveles-1]` del paso nativo `t_nativo`. Cero reservas.
    void paso(std::span<const ark::metrologia::SalidaTemperatura> temp, std::size_t niveles, std::int64_t t_nativo) noexcept {
        const std::size_t n_esc = std::min({niveles, escalas_, temp.size()});
        const bool hipotesis_nula_invertida_f4 = false;
        const bool mu_calor_desde_amplitud = false;
        const bool terminar_con_primera_escala = false;

        bool alguna_termino_en_este_paso = false;

        for (std::size_t j = 0; j < n_esc; ++j) {
            const auto& st = temp[j];
            const auto& f = fichas_[j];
            auto& est = estados_[j];

            // 1. Inicio del calor: en canal_temperatura.hpp, log_beta = -log(s^2), por lo que un aumento de temperatura (s^2 sube,
            //    beta baja) produce una alarma con signo == -1 (p.suma < 0). Si la escala no tiene el fin del calor calibrado
            //    (p. ej. porque buscar_h_con extrapolo), no se usa ese detector no calibrado y se cuenta.
            if (!est.caliente && st.alarma.hay && st.alarma.signo < 0) {
                if (!f.calibrado) {
                    ++cont_.alarmas_ignoradas_no_calibrado;
                    continue;
                }
                est.caliente = true;
                est.t_inicio_calor = t_nativo;
                est.l_star = std::max<std::size_t>(1, st.alarma.l);
                est.ventanas_mu_calor = 0;
                est.suma_mu_calor = 0.0;
                est.mu_calor_listo = false;
                est.mu_calor = NaN;
                est.amplitud_alarma = st.alarma.amplitud;
                est.ventanas_tras_mu_calor = 0;
                est.ultima_g_fin = NaN;
                est.glr_fin.reiniciar();
                ++cont_.inicios_calor;
                continue;
            }

            // 2. Evolucion del episodio caliente en las ventanas posteriores a la alarma (t_nativo > t_inicio_calor)
            if (!est.caliente || !st.emitida || t_nativo <= est.t_inicio_calor) {
                continue;
            }
            if (!std::isfinite(st.log_beta)) {
                est.glr_fin.reiniciar();
                continue;
            }
            ++cont_.ventanas_calientes;

            // Mientras no se han cerrado las l* ventanas posteriores a la alarma, se acumula mu_calor y el fin del calor
            // no puede declararse (E1.7 [DECIDIDO]).
            if (!est.mu_calor_listo) {
                est.suma_mu_calor += st.log_beta;
                ++est.ventanas_mu_calor;
                if (est.ventanas_mu_calor >= est.l_star) {
                    est.mu_calor_listo = true;
                    if (mu_calor_desde_amplitud) {
                        est.mu_calor = f.mu0_log_beta + est.amplitud_alarma;
                    } else {
                        est.mu_calor = est.suma_mu_calor / static_cast<double>(est.ventanas_mu_calor);
                    }
                    const double mu_base_glr = hipotesis_nula_invertida_f4 ? f.mu0_log_beta : est.mu_calor;
                    est.glr_fin.configurar(mu_base_glr, f.s0_log_beta, f.tau_int_log_beta, f.h, hipotesis_nula_invertida_f4, est.l_star);
                }
                continue;
            }

            // Una vez cerradas las l* ventanas posteriores a la alarma, el detector unilateral evalua el regreso hacia mu_0,beta.
            ++est.ventanas_tras_mu_calor;
            const auto p = est.glr_fin.paso(st.log_beta);
            est.ultima_g_fin = p.g;
            if (p.alarma) {
                est.caliente = false;
                est.t_ultimo_fin_calor = t_nativo;
                ++cont_.fines_calor;
                alguna_termino_en_este_paso = true;
            }
        }

        if (terminar_con_primera_escala && alguna_termino_en_este_paso) {
            for (std::size_t j = 0; j < escalas_; ++j) {
                if (estados_[j].caliente) {
                    estados_[j].caliente = false;
                    estados_[j].t_ultimo_fin_calor = t_nativo;
                }
            }
        }
    }
};

}  // namespace ark::lilyth