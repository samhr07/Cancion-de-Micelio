// canal_temperatura.hpp -- canal de temperatura de una escala (ARK-METRO-01): superestadistica directa y deteccion de cambios de beta.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04. Prerregistro: PRUEBAS_JEFE/metrologia/EXPECTATIVAS.md §1.
//
// beta^ = 1/s^2 en ventanas disjuntas de T barras de la escala. Con equilibrio local gaussiano, (T-1) s^2 beta ~ chi^2_{T-1}, y de
// ahi las correcciones exactas por el error de la ventana (k = T - 1):
//   E[log beta]   = media(log beta^) + psi(k/2) - log(k/2),      Var[log beta] = Var(log beta^) - psi1(k/2),
//   <beta>        = media(beta^) (k-2)/k,                         <beta^2>      = media(beta^2) (k-2)(k-4)/k^2,   q = <beta^2>/<beta>^2.
// T* (Beck, Cohen y Swinney 2005, version corregida): el mayor T de la rejilla, contiguo desde el menor, cuya curtosis local media no
// excede la de una gaussiana del MISMO T (simulada) al nivel alfa (una cola).
// Veredicto (decision del Comandante, 2026-11-04): superestadistica = T1 (equilibrio local) y SEPARACION DE ESCALAS contra los nulos
// que ELIGE el usuario (Var[log beta] se mantiene al pasar de ventanas de T* a 2T*, Holm entre nulos). T2 (memoria de log beta entre
// ventanas adyacentes) queda como telemetria: esa memoria la mide el operador de coordinacion.
#pragma once

#include "ark/holm.hpp"   // la unica definicion de Holm del ARK (PRUEBAS_JEFE/holm_cxx, ARK-STAT-01)
#include "ark/metrologia/canal_nivel.hpp"
#include "ark/metrologia/sustitutos.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace ark::metrologia {

inline constexpr std::array<std::size_t, 11> REJILLA_T{8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256};   // candidatos (declarada)

namespace detalle_temp {

struct Ventana { double s2 = NaN, kurt_ex = NaN; };

inline Ventana estadisticos(const double* y, std::size_t T) noexcept {
    double m = 0.0;
    for (std::size_t i = 0; i < T; ++i) m += y[i];
    m /= static_cast<double>(T);
    double m2 = 0.0, m4 = 0.0;
    for (std::size_t i = 0; i < T; ++i) { const double d = (y[i] - m) * (y[i] - m); m2 += d; m4 += d * d; }
    Ventana v{};
    v.s2 = m2 / static_cast<double>(T - 1);
    const double mm2 = m2 / static_cast<double>(T);
    v.kurt_ex = mm2 > 0.0 ? (m4 / static_cast<double>(T)) / (mm2 * mm2) - 3.0 : NaN;
    return v;
}

// ventanas disjuntas consecutivas dentro de cada segmento
inline std::vector<std::vector<Ventana>> ventanas(const std::vector<std::vector<double>>& segs, std::size_t T) {
    std::vector<std::vector<Ventana>> out;
    for (const auto& g : segs) {
        std::vector<Ventana> w;
        for (std::size_t a = 0; a + T <= g.size(); a += T) w.push_back(estadisticos(g.data() + a, T));
        if (!w.empty()) out.push_back(std::move(w));
    }
    return out;
}

// V_k: Var[log beta^] corregida por la ventana (psi1) en ventanas disjuntas de T 2^k barras, k = 0..k_max-1, mientras cada
// escalon tenga >= minimas ventanas (si minimas = 0, se calculan exactamente k_max escalones).
inline std::vector<double> curva_separacion(const std::vector<std::vector<double>>& segs, std::size_t T, std::size_t minimas, std::size_t k_max) {
    std::vector<double> v;
    for (std::size_t k = 0; k < k_max; ++k) {
        const std::size_t Tk = T << k;
        double s1 = 0.0, s2 = 0.0;
        std::size_t n = 0;
        for (const auto& g : segs)
            for (std::size_t a = 0; a + Tk <= g.size(); a += Tk) {
                const double s2v = estadisticos(g.data() + a, Tk).s2;
                if (!(s2v > 0.0)) continue;
                const double lb = -std::log(s2v);
                s1 += lb; s2 += lb * lb; ++n;
            }
        if (n < std::max<std::size_t>(minimas, 2)) break;
        const double m = s1 / static_cast<double>(n);
        v.push_back((s2 - static_cast<double>(n) * m * m) / static_cast<double>(n - 1) - est::trigamma(0.5 * static_cast<double>(Tk - 1)));
    }
    return v;
}

}  // namespace detalle_temp

struct ConfiguracionTemperatura {
    double arl0_pasos = 0.0;          // ARL0 de la escala en pasos de la escala (se pasa a ventanas dividiendo por T)
    double alfa = ALFA_PROYECTO;
    std::size_t alarmas_minimas = 30;
    std::uint64_t semilla = 20261104;
    std::uint32_t nulos_separacion = NULO_NINGUNO;   // ENTRADA del usuario: nulos contra los que se contrasta la separacion
    std::size_t replicas_nulo = 199;                 // sustitutos por nulo (declarado: p minimo 1/200)
};

struct SalidaTemperatura {
    bool emitida = false;             // se cerro una ventana
    double log_beta = NaN;
    double g = NaN;
    EventoAlarma alarma{};
    bool confirmacion = false;        // solo DESLIZANTE: la alarma es la confirmacion por acumulacion del evento abierto
};

template <std::size_t L = 120>
class CanalTemperatura {
    FichaSuperestadistica sup_{};
    FichaEscala glr_ficha_{};          // ficha del GLR sobre log beta^ (m = T * m_escala)
    DetectorGLR<L> glr_{};
    std::size_t T_ = 0;
    std::int64_t m_ = 1;
    double media_ = 0.0, m2_ = 0.0;   // Welford de la ventana abierta
    std::size_t n_ = 0;

public:
    [[nodiscard]] const FichaSuperestadistica& superestadistica() const noexcept { return sup_; }
    [[nodiscard]] const FichaEscala& ficha_glr() const noexcept { return glr_ficha_; }
    [[nodiscard]] std::size_t T() const noexcept { return T_; }
    void reiniciar_flujo() noexcept {
        glr_.configurar(glr_ficha_.nulo.mu0, glr_ficha_.nulo.s0, glr_ficha_.nulo.tau_int, glr_ficha_.det.h);
        n_ = 0; media_ = m2_ = 0.0;
    }

    // Calibracion fuera de linea. segs: historia de la escala (contigua por segmentos).
    void calibrar(int j, std::int64_t m, const std::vector<std::vector<double>>& segs, const ConfiguracionTemperatura& cfg) {
        FichaSuperestadistica s{};
        s.j = j;
        m_ = m;
        const std::size_t minimas = static_cast<std::size_t>(std::ceil(1.0 / cfg.alfa));
        const double z1 = est::norm_ppf(1.0 - cfg.alfa);
        std::mt19937_64 rng(cfg.semilla ^ (0x9E3779B97F4A7C15ULL * static_cast<std::uint64_t>(j + 1)));
        std::normal_distribution<double> N01(0.0, 1.0);
        std::size_t T_star = 0, T_menor = 0;
        double ex_star = NaN, thr_star = NaN, ex_menor = NaN, thr_menor = NaN;
        std::vector<double> buf;
        for (std::size_t T : REJILLA_T) {
            const auto w = detalle_temp::ventanas(segs, T);
            std::size_t nw = 0;
            double ex = 0.0;
            for (const auto& g : w) for (const auto& v : g) if (std::isfinite(v.kurt_ex)) { ex += v.kurt_ex; ++nw; }
            if (nw < minimas) break;
            ex /= static_cast<double>(nw);
            const std::size_t nref = std::max<std::size_t>(4 * nw, 4000);
            buf.resize(T);
            double mr = 0.0, vr = 0.0;
            for (std::size_t r = 0; r < nref; ++r) {
                for (auto& b : buf) b = N01(rng);
                const double k = detalle_temp::estadisticos(buf.data(), T).kurt_ex;
                const double d = k - mr;
                mr += d / static_cast<double>(r + 1);
                vr += d * (k - mr);
            }
            vr /= static_cast<double>(nref - 1);
            const double thr = mr + z1 * std::sqrt(vr / static_cast<double>(nw) + vr / static_cast<double>(nref));
            if (T_menor == 0) { T_menor = T; ex_menor = ex; thr_menor = thr; }
            if (ex <= thr) { T_star = T; ex_star = ex; thr_star = thr; } else break;
        }
        s.t1_equilibrio_local = T_star > 0;
        T_ = T_star ? T_star : T_menor;
        s.T = T_star;
        s.exceso_curtosis = T_star ? ex_star : ex_menor;
        s.exceso_curtosis_umbral = T_star ? thr_star : thr_menor;
        if (T_ == 0) { sup_ = s; glr_ficha_ = FichaEscala{}; glr_ficha_.motivo = "historia corta para el canal de temperatura"; return; }

        // serie de log beta^ por segmento
        const auto w = detalle_temp::ventanas(segs, T_);
        std::vector<std::vector<double>> lb;
        double sb = 0.0, sb2 = 0.0;
        std::size_t nb = 0;
        for (const auto& g : w) {
            std::vector<double> v;
            for (const auto& x : g) {
                if (!(x.s2 > 0.0)) continue;
                const double b = 1.0 / x.s2;
                v.push_back(std::log(b));
                sb += b; sb2 += b * b; ++nb;
            }
            if (!v.empty()) lb.push_back(std::move(v));
        }
        const auto mo = est::momentos(lb);
        s.ventanas = mo.n;
        const double k = static_cast<double>(T_ - 1);
        if (mo.n >= 2 && mo.sd > 0.0) {
            s.media_log_beta = mo.mu + est::digamma(0.5 * k) - std::log(0.5 * k);
            s.var_log_beta = mo.sd * mo.sd - est::trigamma(0.5 * k);
            const double bm = sb / static_cast<double>(nb), b2m = sb2 / static_cast<double>(nb);
            s.beta_media = T_ > 3 ? bm * (k - 2.0) / k : NaN;
            const double b2 = T_ > 5 ? b2m * (k - 2.0) * (k - 4.0) / (k * k) : NaN;
            s.q_tsallis = b2 / (s.beta_media * s.beta_media);
            std::size_t pares = 0;
            s.rho_log_beta = est::acf_agrupada(lb, mo.mu, mo.sd, 1, pares);
            s.rho_umbral = pares ? z1 / std::sqrt(static_cast<double>(pares)) : NaN;
            s.t2_separacion = s.rho_log_beta > s.rho_umbral;
            const auto ti = est::tau_integrado(lb, mo.mu, mo.sd, cfg.alfa, mo.n / 4);
            s.tau_int_log_beta = ti.tau;
            s.t_beta_pasos = ti.tau * static_cast<double>(T_);
        }
        // separacion de escalas contra los nulos elegidos
        s.nulos = cfg.nulos_separacion;
        s.replicas_nulo = cfg.replicas_nulo;
        const auto v = detalle_temp::curva_separacion(segs, T_, minimas, FichaSuperestadistica::K_SEP);
        s.k_sep = v.size();
        for (std::size_t kk = 0; kk < v.size(); ++kk) s.v_k[kk] = v[kk];
        for (auto& f : s.p_k) f.fill(NaN);
        for (auto& f : s.r_nulo_mediana) f.fill(NaN);
        s.p_holm.fill(NaN);
        if (s.t1_equilibrio_local && v.size() >= 2 && v[0] > 0.0 && cfg.nulos_separacion != NULO_NINGUNO) {
            s.separacion_contrastada = true;
            double mu_y = 0.0;
            std::size_t ny = 0;
            for (const auto& g : segs) for (double x : g) { mu_y += x; ++ny; }
            mu_y /= static_cast<double>(ny);
            std::vector<double> p1;
            std::vector<std::size_t> idx;
            // Enmienda 7: cada escalon k con SU nulo. s_k = V_{k+1} / V_k pregunta si beta persiste mas alla de L_k = T* 2^k; el
            // nulo BLOQUES usa b = L_k (conserva lo que hay dentro de L_k, destruye el orden mas alla). PERMUTACION y SIGNOS no tienen
            // longitud: un mismo juego de sustitutos sirve para todos los k.
            const std::size_t nk = v.size() - 1;   // contrastes k = 0..nk-1
            for (std::size_t i = 0; i < 3; ++i) {
                const std::uint32_t nulo = NULOS_TODOS[i];
                if (!(cfg.nulos_separacion & nulo)) continue;
                std::vector<std::vector<double>> sn(nk);
                auto anotar = [&](const std::vector<double>& vn, std::size_t kk) {
                    sn[kk].push_back(kk + 1 < vn.size() && vn[kk] > 0.0 ? vn[kk + 1] / vn[kk] : -INFINITY);   // sin variacion de beta: no supera
                };
                if (nulo == NULO_BLOQUES) {
                    for (std::size_t kk = 0; kk < nk; ++kk)
                        for (std::size_t b = 0; b < cfg.replicas_nulo; ++b)
                            anotar(detalle_temp::curva_separacion(sustituto(segs, nulo, T_ << kk, mu_y, rng), T_, 0, kk + 2), kk);
                } else {
                    for (std::size_t b = 0; b < cfg.replicas_nulo; ++b) {
                        const auto vn = detalle_temp::curva_separacion(sustituto(segs, nulo, T_, mu_y, rng), T_, 0, v.size());
                        for (std::size_t kk = 0; kk < nk; ++kk) anotar(vn, kk);
                    }
                }
                std::size_t esc = 0;
                bool contiguo = true;
                for (std::size_t kk = 0; kk < nk; ++kk) {
                    const double r = v[kk] > 0.0 ? v[kk + 1] / v[kk] : -INFINITY;
                    std::size_t ge = 0;
                    for (double x : sn[kk]) ge += x >= r;
                    s.p_k[i][kk] = static_cast<double>(1 + ge) / static_cast<double>(1 + sn[kk].size());
                    s.r_nulo_mediana[i][kk] = est::cuantil(sn[kk], 0.5);
                    if (contiguo && s.p_k[i][kk] <= cfg.alfa) esc = T_ << (kk + 1); else contiguo = false;
                }
                s.escala_separacion[i] = esc;   // 0: beta no persiste mas alla de T*
                p1.push_back(s.p_k[i][0]);
                idx.push_back(i);
            }
            const auto ph = ark::holm(p1);
            s.separacion = !ph.empty();
            for (std::size_t q = 0; q < ph.size(); ++q) { s.p_holm[idx[q]] = ph[q]; s.separacion = s.separacion && ph[q] <= cfg.alfa; }
        }
        s.es_superestadistica = s.t1_equilibrio_local && s.separacion;
        sup_ = s;

        // GLR sobre log beta^ con su propio nulo
        CanalNivel<L, 8> aux;   // reutiliza la calibracion del canal de nivel (nulo, tau_int, h); sin tau50
        ConfiguracionNivel cn{};
        cn.arl0_pasos = cfg.arl0_pasos / static_cast<double>(T_);
        cn.alfa = cfg.alfa;
        cn.alarmas_minimas = cfg.alarmas_minimas;
        cn.medir_tau50 = false;
        glr_ficha_ = aux.calibrar(j, m * static_cast<std::int64_t>(T_), lb, cn);
        reiniciar_flujo();
    }

    // Ciclo caliente: una barra de la escala. t_nativo = indice nativo de su ultima barra nativa.
    SalidaTemperatura paso(double y, bool valido, std::int64_t t_nativo) noexcept {
        SalidaTemperatura s{};
        if (T_ == 0) return s;
        if (!valido || !std::isfinite(y)) { n_ = 0; media_ = m2_ = 0.0; glr_.reiniciar(); return s; }
        ++n_;
        const double d = y - media_;
        media_ += d / static_cast<double>(n_);
        m2_ += d * (y - media_);
        if (n_ < T_) return s;
        const double s2 = m2_ / static_cast<double>(T_ - 1);
        n_ = 0; media_ = m2_ = 0.0;
        s.emitida = true;
        if (!(s2 > 0.0)) return s;
        s.log_beta = -std::log(s2);
        if (!glr_ficha_.detectable) return s;
        const auto p = glr_.paso(s.log_beta);
        s.g = p.g;
        if (p.alarma) {
            s.alarma.hay = true;
            s.alarma.signo = p.suma > 0.0 ? 1 : -1;     // + : beta sube (se enfria), - : beta baja (se calienta)
            s.alarma.t = t_nativo;
            s.alarma.l = p.l;
            s.alarma.t_cambio = t_nativo - static_cast<std::int64_t>(p.l) * glr_ficha_.m + 1;
            s.alarma.amplitud = p.suma / static_cast<double>(p.l) * glr_ficha_.nulo.s0;
            s.alarma.se_amplitud = glr_ficha_.nulo.s0 * std::sqrt(glr_ficha_.nulo.tau_int / static_cast<double>(p.l));
            s.alarma.g = p.g;
            s.alarma.t_pendiente = p.t_pendiente;
        }
        return s;
    }
};

}  // namespace ark::metrologia
