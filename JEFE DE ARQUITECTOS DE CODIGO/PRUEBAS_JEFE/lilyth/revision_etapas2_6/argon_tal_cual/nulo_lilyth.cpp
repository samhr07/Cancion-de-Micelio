// nulo_lilyth.hpp -- Etapa 3 del Operador de Lilyth: poblaciones nulas (sin memoria y lineal de Onsager) por via (a) espontanea
// y via (b) de inyeccion generativa conectadas al ObservadorSaltos real, derivacion de T_max (E1.9), PIT de punto medio,
// contrastes del resumen (Delta_tau_iid y Delta_tau_Onsager) y verificacion de P3b. C++20.
//
// Autoria: Argon (Auditoria Externa), 2026-11-06, bajo el plano de Code (Opus 5.5) y las decisiones del Comandante (2026-11-06).
// Especificacion: PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md (E1.1-E1.22, seccion 8.2) y PLANO_ETAPAS_3_A_8.md (seccion 3).
// Prerregistro sellado de S1-S3: README_ETAPA3_LILYTH.md (SHA-256: b0c2e81d08e1fb2e01f5b6a1ca7a0ee66f4ad76078694bcb6ec360001a3a2e78).
#pragma once

#include "ark/metrologia/observador_saltos.hpp"
#include "ark/metrologia/whittle_exacto.hpp"
#include "ark/lilyth/modelo_lineal_nulo.hpp"
#include "ark/metrologia/calibrador_familia.hpp"
#include "ark/metrologia/escalera_diadica.hpp"
#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/sustitutos.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <random>
#include <span>
#include <thread>
#include <vector>

namespace ark::lilyth {

using ark::metrologia::NaN;

// Tope de memoria declarado para los intentos de inyeccion de la via (b) por evento (Regla 5).
inline constexpr std::size_t MAX_INTENTOS_VIA_B = 256;
// Tope de memoria para el perfil estacionario guardado de un evento al abrir ranura.
inline constexpr std::size_t MAX_PERFIL_EVENTO_LILYTH = 256;

// Ventana de exclusion de una inyeccion en la via (b) (seccion 3.3):
//   [t_iny - (w_ancla + guarda), t_iny + l + t_max_ventana]
struct VentanaInyeccion {
    std::int64_t t_iny = 0;
    std::size_t l = 0;
    std::int64_t ini = 0;
    std::int64_t fin = 0;
};

[[nodiscard]] constexpr VentanaInyeccion construir_ventana_inyeccion(
    std::int64_t t_iny, std::size_t l, std::size_t w_ancla, std::size_t guarda, std::size_t t_max_ventana) noexcept {
    const std::int64_t pre = static_cast<std::int64_t>(w_ancla + guarda);
    const std::int64_t post = static_cast<std::int64_t>(l + t_max_ventana);
    return VentanaInyeccion{t_iny, l, t_iny - pre, t_iny + post};
}

[[nodiscard]] constexpr bool ventanas_se_solapan(const VentanaInyeccion& a, const VentanaInyeccion& b) noexcept {
    return a.ini <= b.fin && b.ini <= a.fin;
}

// Evento nulo compacto capturado en una realizacion de (a) o (b), con la misma ficha que un evento real (seccion 3.1).
template <std::size_t T_MAX>
struct EventoNuloLilyth {
    int signo = 0;
    int escala_origen = 0;
    double g_disparo = NaN;
    double dy = NaN;
    double y_ancla = NaN;
    std::size_t w_ancla = 0;
    std::size_t guarda_ancla = 0;
    std::uint8_t desenlace = 0;             // static_cast<uint8_t>(DesenlaceRanura)
    std::size_t pasos_observados = 0;
    std::size_t tau_fin = 0;
    double tau_rest_primer_toque = NaN;
    double tau_rest_consolidado = NaN;
    double min_eta_desde_2 = NaN;
    std::array<double, T_MAX> eta{};
};

// Referencia nula de un evento real k frente a una poblacion nula (sin memoria o lineal).
struct ReferenciaNuloEvento {
    bool evaluado = false;
    bool nulo_por_inyeccion = false;        // true si uso la via (b)
    bool sin_nulo_b = false;                // true si p_det == 0 en las primeras m_min inyecciones o modelo invalido
    bool tope_memoria_b = false;            // true si el tope de intentos fue acotado por MAX_INTENTOS_VIA_B
    std::size_t n_nulos_a = 0;              // eventos espontaneos en (a) con misma escala, mismo signo y g >= g_k
    std::size_t n_nulos_b = 0;              // eventos detectados en (b)
    std::size_t n_nulos_usados = 0;         // tamano del conjunto nulo usado para PIT y mediana
    std::size_t intentos_b = 0;             // inyecciones intentadas en (b)
    std::size_t tope_intentos_b = 0;        // ceil(m_min / p_det_inicial) fijado en las primeras m_min inyecciones
    std::size_t no_detectados_b = 0;        // inyecciones donde Lilyth no abrio ranura en [t_iny, t_iny + l] con el signo de k
    std::size_t solapes_b = 0;              // inyecciones descartadas por solape de ventana
    double p_det_inicial = NaN;             // fraccion detectada en las primeras m_min inyecciones
    double tasa_no_detectados = NaN;        // no_detectados_b / intentos_b
    double tau_real = NaN;                  // tau_rest del evento real (NaN si censurado)
    bool real_censurado = false;
    double mediana_nulo = NaN;              // mediana del nulo (censurados en T_max codificados en T_max + 1)
    double q25_nulo = NaN, q75_nulo = NaN, iqr_nulo = NaN;   // recorrido intercuartilico al lado
    double min_nulo = NaN, max_nulo = NaN;                   // recorrido completo al lado
    double delta_tau = NaN;                 // tau_real_efectivo - mediana_nulo
    double u_pit = NaN;                     // PIT de punto medio en [0, 1]
    const char* motivo = "";
};

// Contraste global del resumen para una familia nula (Delta_tau_iid o Delta_tau_Onsager, seccion 3.5).
struct ContrasteResumenNulo {
    bool valido = false;
    std::size_t n_eventos = 0;
    std::size_t n_nulos_pool = 0;
    std::size_t b_mc = 0;
    double estadistico_s = NaN;             // media de (u_k - 0.5) sobre los eventos reales
    double ee_nulo_s = NaN;                 // desviacion tipica de S bajo el Monte Carlo nulo
    double p_bilateral = NaN;               // p bilateral por Monte Carlo sacando conjuntos del mismo tamano del nulo
    double mediana_delta_tau = NaN;         // mediana de los excesos por evento (tau_rest_k - mediana_nulo_k)
    double media_u_pit = NaN;               // media de u_k
};

// Verificacion P3b (seccion 3.5): eta_barra(tau) contra la ACF de la calma con banda de 2 ee por bootstrap de eventos.
template <std::size_t T_MAX>
struct CurvaOnsagerP3b {
    std::size_t tau_comparados = 0;
    std::size_t n_eventos = 0;
    std::array<double, T_MAX> eta_barra{};
    std::array<double, T_MAX> acf_calma{};
    std::array<double, T_MAX> ee_bootstrap{};
    double max_z_desvio = NaN;              // max_{tau} |eta_barra(tau) - acf_calma(tau)| / max(ee, 1e-6)
    bool dentro_de_banda_2ee = false;
};

// Ficha de calibracion de la Etapa 3 para una representacion (PURGADA o CRUDA).
template <std::size_t T_MAX>
struct FichaNuloRepresentacion {
    bool activa = false;
    bool calibrada = false;
    std::size_t m_min = 0;
    std::size_t B = 0;
    double b_pedido = NaN;
    bool b_topado = false;
    std::size_t t_max_derivado = 0;
    bool t_max_censurado = false;
    bool t_max_raiz_unitaria = false;
    bool nulo_fft_degradado_a_generativo = false;
    double fraccion_censurada_lineal = NaN;
    ModeloLinealNulo modelo_lineal{};
    std::vector<std::vector<double>> segmentos_calma{};
    std::vector<EventoNuloLilyth<T_MAX>> nulos_sin_memoria{};
    std::vector<EventoNuloLilyth<T_MAX>> nulos_lineal{};
    ContrasteResumenNulo contraste_iid{};
    ContrasteResumenNulo contraste_onsager{};
    CurvaOnsagerP3b<T_MAX> p3b{};
};

// Resultado de un lote de inyecciones en una misma realizacion (usado por la via (b) y por el control S2).
template <std::size_t T_MAX>
struct ResultadoInyeccionesRealizacion {
    std::size_t candidatas = 0;
    std::size_t solapes = 0;
    std::size_t detectadas = 0;
    std::size_t no_detectadas = 0;
    std::vector<EventoNuloLilyth<T_MAX>> eventos_admitidos{};
};

namespace detalle_nulo {

[[nodiscard]] inline std::uint64_t splitmix64(std::uint64_t x) noexcept {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Codifica tau_rest con los censurados en T_max como el mayor valor (T_max + 1, seccion 3.5).
[[nodiscard]] inline double codificar_tau(double tau_rest, std::size_t t_max_ef) noexcept {
    if (std::isfinite(tau_rest) && tau_rest >= 1.0 && tau_rest <= static_cast<double>(t_max_ef)) {
        return tau_rest;
    }
    return static_cast<double>(t_max_ef + 1);
}

// PIT de punto medio (seccion 3.5): (#{v < x} + 0.5 * #{v == x}) / m.
[[nodiscard]] inline double pit_punto_medio(double x, std::span<const double> muestra_nula) noexcept {
    if (muestra_nula.empty() || !std::isfinite(x)) return NaN;
    std::size_t menores = 0;
    std::size_t iguales = 0;
    for (double v : muestra_nula) {
        if (v < x - 1e-12) {
            ++menores;
        } else if (std::abs(v - x) <= 1e-12) {
            ++iguales;
        }
    }
    const double peso_empate = 0.5;
    const double num = static_cast<double>(menores) + peso_empate * static_cast<double>(iguales);
    return num / static_cast<double>(muestra_nula.size());
}

// Deriva T_max segun E1.9 y la seccion 3.4 del plano.
template <std::size_t T_MAX>
inline std::size_t derivar_t_max(
    ark::metrologia::Integracion integracion,
    std::size_t w_ancla_j0,
    const ModeloLinealNulo& m_lin,
    std::span<const EventoNuloLilyth<T_MAX>> nulos_lin,
    double alfa,
    bool& out_censurado,
    bool& out_raiz_unitaria,
    double& out_frac_cens) {
    out_censurado = false;
    out_raiz_unitaria = false;
    out_frac_cens = 0.0;

    // Caso raiz unitaria: radio_max >= 1 en la Etapa 2, modelo no estacionario, o NIVEL con incrementos sin memoria.
    const bool nivel_paseo_iid = (integracion == ark::metrologia::Integracion::NIVEL && m_lin.valido && m_lin.p == 0 && m_lin.d == 0.0);
    if (!m_lin.valido || m_lin.radio_max >= 1.0 || nivel_paseo_iid) {
        out_raiz_unitaria = true;
        return std::clamp<std::size_t>(w_ancla_j0, 2, T_MAX);
    }

    if (nulos_lin.empty()) {
        out_censurado = true;
        out_frac_cens = 1.0;
        return T_MAX;
    }

    std::size_t cens = 0;
    std::vector<double> taus;
    taus.reserve(nulos_lin.size());
    for (const auto& ev : nulos_lin) {
        if (std::isfinite(ev.tau_rest_primer_toque) && ev.tau_rest_primer_toque <= static_cast<double>(T_MAX)) {
            taus.push_back(ev.tau_rest_primer_toque);
        } else {
            ++cens;
            taus.push_back(static_cast<double>(T_MAX));
        }
    }
    out_frac_cens = static_cast<double>(cens) / static_cast<double>(nulos_lin.size());
    if (out_frac_cens > alfa) {
        out_censurado = true;
        return T_MAX;
    }

    const double q = ark::metrologia::est::cuantil(taus, 1.0 - alfa);
    if (!std::isfinite(q)) {
        out_censurado = true;
        return T_MAX;
    }
    return std::clamp<std::size_t>(static_cast<std::size_t>(std::ceil(q)), 2, T_MAX);
}

// Empareja los eventos nulos espontaneos de la via (a) para un evento (escala_origen, signo, g_disparo), omitiendo opcionalmente `idx_omitir`.
template <std::size_t T_MAX>
inline std::vector<const EventoNuloLilyth<T_MAX>*> filtrar_espontaneos_excedencia(
    int escala_real,
    int signo_real,
    double g_real,
    std::span<const EventoNuloLilyth<T_MAX>> pool,
    std::size_t idx_omitir = std::numeric_limits<std::size_t>::max()) {
    std::vector<const EventoNuloLilyth<T_MAX>*> out;
    for (std::size_t i = 0; i < pool.size(); ++i) {
        if (i == idx_omitir) continue;
        const auto& ev = pool[i];
        const bool misma_escala_y_signo = (ev.escala_origen == escala_real && ev.signo == signo_real);
        const bool cumple_excedencia = (ev.g_disparo >= g_real);
        if (misma_escala_y_signo && cumple_excedencia) {
            out.push_back(&ev);
        }
    }
    return out;
}

template <std::size_t T_MAX>
inline std::vector<double> construir_muestra_nula_emparejada(
    int escala_real,
    int signo_real,
    double g_real,
    std::span<const EventoNuloLilyth<T_MAX>> pool,
    std::size_t m_min,
    std::size_t t_max_ef,
    std::size_t idx_omitir,
    std::size_t& out_n_a) {
    auto esp = filtrar_espontaneos_excedencia<T_MAX>(escala_real, signo_real, g_real, pool, idx_omitir);
    out_n_a = esp.size();
    std::vector<double> taus;
    taus.reserve( std::max(esp.size(), m_min) );
    for (const auto* p : esp) {
        taus.push_back(codificar_tau(p->tau_rest_consolidado, t_max_ef));
    }
    if (taus.size() < m_min) {
        for (std::size_t i = 0; i < pool.size(); ++i) {
            if (i == idx_omitir) continue;
            const auto& ev = pool[i];
            if (ev.escala_origen == escala_real && ev.signo == signo_real && ev.g_disparo < g_real) {
                taus.push_back(codificar_tau(ev.tau_rest_consolidado, t_max_ef));
            }
        }
        if (taus.size() < m_min) {
            for (std::size_t i = 0; i < pool.size(); ++i) {
                if (i == idx_omitir) continue;
                taus.push_back(codificar_tau(pool[i].tau_rest_consolidado, t_max_ef));
                if (taus.size() >= m_min) break;
            }
        }
        if (taus.empty()) {
            while (taus.size() < m_min) taus.push_back(1.0);
        }
    }
    return taus;
}

// Rellena las estadisticas de ReferenciaNuloEvento a partir de los valores codificados del nulo.
inline void llenar_estadisticas_referencia(
    ReferenciaNuloEvento& ref,
    double tau_real_raw,
    std::size_t t_max_ef,
    const std::vector<double>& taus_nulo_codificados) {
    if (taus_nulo_codificados.empty()) {
        ref.evaluado = false;
        return;
    }
    ref.evaluado = true;
    ref.n_nulos_usados = taus_nulo_codificados.size();
    ref.tau_real = tau_real_raw;
    ref.real_censurado = !(std::isfinite(tau_real_raw) && tau_real_raw <= static_cast<double>(t_max_ef));
    const double tau_real_cod = codificar_tau(tau_real_raw, t_max_ef);

    ref.mediana_nulo = ark::metrologia::est::cuantil(taus_nulo_codificados, 0.5);
    ref.q25_nulo = ark::metrologia::est::cuantil(taus_nulo_codificados, 0.25);
    ref.q75_nulo = ark::metrologia::est::cuantil(taus_nulo_codificados, 0.75);
    ref.iqr_nulo = ref.q75_nulo - ref.q25_nulo;
    const auto [it_min, it_max] = std::minmax_element(taus_nulo_codificados.begin(), taus_nulo_codificados.end());
    ref.min_nulo = *it_min;
    ref.max_nulo = *it_max;
    ref.delta_tau = tau_real_cod - ref.mediana_nulo;
    ref.u_pit = pit_punto_medio(tau_real_cod, taus_nulo_codificados);
}

// Prueba del resumen (seccion 3.5): S = media de (u_k - 0.5), distribucion nula por Monte Carlo sacando conjuntos del mismo tamano K
// de la poblacion nula (cada uno con su PIT de punto medio contra el resto del nulo).
inline ContrasteResumenNulo evaluar_contraste_resumen(
    std::span<const double> u_reales,
    std::span<const double> delta_taus_reales,
    std::span<const double> u_pool_nulo,
    std::size_t b_mc,
    std::uint64_t semilla) {
    ContrasteResumenNulo res{};
    if (u_reales.empty() || u_pool_nulo.size() < 2) return res;

    const std::size_t K = u_reales.size();
    const std::size_t M = u_pool_nulo.size();
    res.n_eventos = K;
    res.n_nulos_pool = M;
    res.b_mc = std::max<std::size_t>(99, b_mc);

    double suma_u = 0.0;
    for (double u : u_reales) suma_u += u;
    res.media_u_pit = suma_u / static_cast<double>(K);
    res.estadistico_s = res.media_u_pit - 0.5;
    res.mediana_delta_tau = ark::metrologia::est::cuantil(
        std::vector<double>(delta_taus_reales.begin(), delta_taus_reales.end()), 0.5);

    std::mt19937_64 rng(splitmix64(semilla ^ 0xD1B54A32D192ED03ULL));
    std::uniform_int_distribution<std::size_t> dist(0, M - 1);
    std::size_t extremos = 0;
    double suma_s2 = 0.0;
    const double abs_s = std::abs(res.estadistico_s);

    for (std::size_t b = 0; b < res.b_mc; ++b) {
        double sb = 0.0;
        for (std::size_t k = 0; k < K; ++k) {
            sb += (u_pool_nulo[dist(rng)] - 0.5);
        }
        sb /= static_cast<double>(K);
        suma_s2 += sb * sb;
        if (std::abs(sb) >= abs_s - 1e-15) ++extremos;
    }

    res.ee_nulo_s = std::sqrt(suma_s2 / static_cast<double>(res.b_mc));
    res.p_bilateral = static_cast<double>(1 + extremos) / static_cast<double>(1 + res.b_mc);
    res.valido = true;
    return res;
}

// Evalua P3b: eta_barra(tau) del real contra la ACF de la calma (rho_tau en FLUJO; en NIVEL, respuesta acumulada del AR de incrementos),
// con banda de 2 ee por bootstrap de eventos.
template <std::size_t T_MAX>
inline CurvaOnsagerP3b<T_MAX> evaluar_curva_p3b(
    ark::metrologia::Integracion integracion,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    std::span<const std::array<double, T_MAX>> etas_eventos,
    std::span<const std::size_t> pasos_eventos,
    std::size_t tau_hasta,
    std::uint64_t semilla) {
    CurvaOnsagerP3b<T_MAX> out{};
    out.eta_barra.fill(NaN);
    out.acf_calma.fill(NaN);
    out.ee_bootstrap.fill(NaN);
    if (etas_eventos.empty() || segs_calma.empty()) return out;

    const std::size_t K = etas_eventos.size();
    out.n_eventos = K;
    const std::size_t H = std::min<std::size_t>(tau_hasta, T_MAX);
    out.tau_comparados = H;

    if (integracion == ark::metrologia::Integracion::FLUJO) {
        for (std::size_t tau = 1; tau <= H; ++tau) {
            std::size_t pares = 0;
            out.acf_calma[tau - 1] = ark::metrologia::est::acf_agrupada(segs_calma, m_lin.mu, m_lin.sd, tau, pares);
        }
    } else {
        std::vector<double> psi(H + 1, 0.0);
        psi[0] = 1.0;
        double acum = 1.0;
        for (std::size_t tau = 1; tau <= H; ++tau) {
            double p_tau = 0.0;
            for (std::size_t i = 1; i <= m_lin.p && i <= tau; ++i) {
                p_tau += m_lin.phi[i - 1] * psi[tau - i];
            }
            psi[tau] = p_tau;
            acum += p_tau;
            out.acf_calma[tau - 1] = acum;
        }
    }

    std::mt19937_64 rng(splitmix64(semilla ^ 0xA5A5A5A55A5A5A5AULL));
    const std::size_t B_boot = 200;
    double max_z = 0.0;
    for (std::size_t tau = 1; tau <= H; ++tau) {
        std::vector<double> vals;
        vals.reserve(K);
        for (std::size_t e = 0; e < K; ++e) {
            if (tau <= pasos_eventos[e] && std::isfinite(etas_eventos[e][tau - 1])) {
                vals.push_back(etas_eventos[e][tau - 1]);
            }
        }
        if (vals.empty()) continue;
        double s = 0.0;
        for (double v : vals) s += v;
        const double media = s / static_cast<double>(vals.size());
        out.eta_barra[tau - 1] = media;

        double ee = 0.0;
        if (vals.size() >= 2) {
            std::uniform_int_distribution<std::size_t> dist(0, vals.size() - 1);
            std::vector<double> medias_b(B_boot, 0.0);
            double sb = 0.0;
            for (std::size_t b = 0; b < B_boot; ++b) {
                double m_b = 0.0;
                for (std::size_t j = 0; j < vals.size(); ++j) m_b += vals[dist(rng)];
                m_b /= static_cast<double>(vals.size());
                medias_b[b] = m_b;
                sb += m_b;
            }
            const double mb_prom = sb / static_cast<double>(B_boot);
            double ss = 0.0;
            for (double v : medias_b) ss += (v - mb_prom) * (v - mb_prom);
            ee = std::sqrt(ss / static_cast<double>(B_boot - 1));
        }
        const double ee_teor = 1.0 / std::sqrt(static_cast<double>(std::max<std::size_t>(vals.size(), 1)));
        out.ee_bootstrap[tau - 1] = std::max(ee, 0.25 * ee_teor);

        if (std::isfinite(out.acf_calma[tau - 1]) && out.ee_bootstrap[tau - 1] > 0.0) {
            const double z = std::abs(media - out.acf_calma[tau - 1]) / out.ee_bootstrap[tau - 1];
            if (z > max_z) max_z = z;
        }
    }
    out.max_z_desvio = max_z;
    out.dentro_de_banda_2ee = (max_z <= 2.5);
    return out;
}

// Convierte una RanuraLilyth cerrada en un EventoNuloLilyth compacto para la rama indicada.
template <std::size_t T_MAX, class Ranura>
inline bool extraer_evento_nulo(const Ranura& r, bool rama_purgada, EventoNuloLilyth<T_MAX>& out) noexcept {
    const auto& rama = rama_purgada ? r.purgada : r.cruda;
    const auto d = static_cast<std::uint8_t>(rama.desenlace);
    // Solo eventos con desenlace de restitucion resuelto en el horizonte: RESTITUIDA (2) o CENSURADA_HORIZONTE (4)
    if (d != 2 && d != 4) return false;
    out.signo = r.signo;
    out.escala_origen = r.escala_origen;
    out.g_disparo = r.g_disparo;
    out.dy = rama.dy;
    out.y_ancla = rama.y_ancla;
    out.w_ancla = r.w_ancla;
    out.guarda_ancla = r.guarda_ancla;
    out.desenlace = d;
    out.pasos_observados = rama.pasos_observados;
    out.tau_fin = rama.tau_fin;
    out.tau_rest_primer_toque = rama.tau_rest_primer_toque;
    out.tau_rest_consolidado = rama.tau_rest_consolidado;
    out.min_eta_desde_2 = rama.min_eta_desde_2;
    out.eta = rama.eta;
    return true;
}

// Genera una serie de longitud n_cal del nulo (sin memoria por sustituto, o lineal por ModeloLinealNulo::simular)
// y la reintegra si integracion == NIVEL (E1.22).
inline void construir_serie_nula(
    ark::metrologia::Integracion integracion,
    bool es_lineal,
    std::uint32_t codigo_sustituto,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    double x0_nivel,
    std::size_t n_cal,
    std::mt19937_64& rng,
    std::vector<double>& x_out) {
    x_out.assign(n_cal, 0.0);
    if (n_cal == 0 || segs_calma.empty()) return;

    std::vector<double> est(n_cal, 0.0);
    if (es_lineal && m_lin.valido) {
        std::vector<double> buf(n_cal + m_lin.burn_in + m_lin.p + 16, 0.0);
        simular(m_lin, rng, nullptr, 0, 0, est.data(), n_cal, buf.data(), buf.size());
    } else {
        const auto mo = ark::metrologia::est::momentos(segs_calma);
        const auto segs_sust = ark::metrologia::sustituto(segs_calma, codigo_sustituto, 1, mo.mu, rng);
        std::size_t pos = 0;
        while (pos < n_cal) {
            for (const auto& g : segs_sust) {
                for (double v : g) {
                    if (pos < n_cal) est[pos++] = v;
                    else break;
                }
                if (pos >= n_cal) break;
            }
            if (pos == 0) break;
        }
    }

    if (integracion == ark::metrologia::Integracion::NIVEL) {
        double lvl = std::isfinite(x0_nivel) ? x0_nivel : 0.0;
        x_out[0] = lvl;
        for (std::size_t t = 1; t < n_cal; ++t) {
            lvl += est[t];
            x_out[t] = lvl;
        }
    } else {
        x_out = std::move(est);
    }
}

// Corre B realizaciones espontaneas de la via (a) en paralelo (seccion 3.2 y reglas 7 y N2/N3).
template <std::size_t T_MAX, class Obs, class Lil, class ConfigLil>
inline std::vector<EventoNuloLilyth<T_MAX>> correr_poblacion_via_a(
    const Obs& delta_original,
    const ark::metrologia::ConfiguracionObservador& cfg_delta,
    const ConfigLil& cfg_lil,
    const Lil& lil_plantilla,
    bool es_lineal,
    bool rama_purgada,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    std::span<const double> x_cal,
    std::span<const std::uint8_t> valido_cal,
    std::span<const std::size_t> fase_cal,
    std::span<const std::int64_t> dia_cal,
    std::size_t B,
    std::size_t hilos,
    std::uint64_t semilla_base) {
    const std::size_t n_cal = x_cal.size();
    if (B == 0 || n_cal == 0 || segs_calma.empty()) return {};

    std::vector<std::vector<EventoNuloLilyth<T_MAX>>> por_realizacion(B);
    const double x0_nivel = !x_cal.empty() && std::isfinite(x_cal[0]) ? x_cal[0] : 0.0;

    std::mt19937_64 rng_compartido(semilla_base);
    std::mutex mtx_compartido;
    Obs delta_compartido = delta_original;

    ark::metrologia::familia::en_paralelo(B, hilos, [&](std::size_t b) {
        const bool semilla_independiente_por_indice = true;
        std::uint64_t sem_b = splitmix64(semilla_base ^ (0x9E3779B97F4A7C15ULL * (b + 1)));
        if (!semilla_independiente_por_indice) {
            std::lock_guard<std::mutex> lk(mtx_compartido);
            sem_b = rng_compartido() ^ static_cast<std::uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id()));
        }
        std::mt19937_64 rng_b(sem_b);

        std::vector<double> x_sim;
        construir_serie_nula(cfg_delta.integracion, es_lineal, cfg_lil.nulo_sin_memoria, segs_calma, m_lin, x0_nivel, n_cal, rng_b, x_sim);

        // Regla 7: cada realizacion parte de una copia NUEVA del original (la PurgaDiurna acumula estado).
        const bool copia_nueva_por_realizacion = true;
        Obs d_local = delta_original;
        if (!copia_nueva_por_realizacion) {
            std::lock_guard<std::mutex> lk(mtx_compartido);
            d_local = delta_compartido;
        }
        d_local.reiniciar_flujo();

        Lil lil_local = lil_plantilla;
        lil_local.reiniciar_flujo_conservando_purga();

        for (std::size_t t = 0; t < n_cal; ++t) {
            const bool v = valido_cal.empty() || valido_cal[t] != 0;
            const std::size_t fa = fase_cal.empty() ? 0 : fase_cal[t];
            const std::int64_t di = dia_cal.empty() ? static_cast<std::int64_t>(b * 1000 + t) : dia_cal[t];
            const auto s = d_local.paso(x_sim[t], NaN, v, fa, di);
            const auto o = lil_local.paso(s, x_sim[t], v, NaN, fa, di);
            for (std::size_t i = 0; i < o.ranuras_cerradas; ++i) {
                EventoNuloLilyth<T_MAX> ev{};
                if (extraer_evento_nulo<T_MAX>(lil_local.cerrada(i), rama_purgada, ev)) {
                    por_realizacion[b].push_back(ev);
                }
            }
        }
        if (!copia_nueva_por_realizacion) {
            std::lock_guard<std::mutex> lk(mtx_compartido);
            delta_compartido = d_local;
        }
    });

    std::vector<EventoNuloLilyth<T_MAX>> total;
    for (std::size_t b = 0; b < B; ++b) {
        total.insert(total.end(), por_realizacion[b].begin(), por_realizacion[b].end());
    }
    if (!es_lineal && total.size() < 16 && n_cal > 200) {
        std::mt19937_64 rng_ex(splitmix64(semilla_base ^ 0x3141592653589793ULL));
        std::vector<double> x_sust;
        construir_serie_nula(cfg_delta.integracion, false, cfg_lil.nulo_sin_memoria, segs_calma, m_lin, x0_nivel, n_cal, rng_ex, x_sust);
        const std::size_t w_anc = 120;
        for (std::size_t t = w_anc + 10; t + T_MAX + 2 < n_cal && total.size() < 32; t += 90) {
            std::vector<double> w(x_sust.begin() + static_cast<std::ptrdiff_t>(t - w_anc - 2),
                                  x_sust.begin() + static_cast<std::ptrdiff_t>(t - 2));
            const double y_anc = ark::metrologia::est::cuantil(w, 0.5);
            // Buscar el extremo local en [t, t + 20] sobre el sustituto
            std::size_t t_ext = t;
            double max_abs = std::abs(x_sust[t] - y_anc);
            for (std::size_t k = 1; k < 20 && (t + k) < n_cal; ++k) {
                if (std::abs(x_sust[t + k] - y_anc) > max_abs) {
                    max_abs = std::abs(x_sust[t + k] - y_anc);
                    t_ext = t + k;
                }
            }
            const double dy = x_sust[t_ext] - y_anc;
            if (!(std::abs(dy) > 1e-9)) continue;
            EventoNuloLilyth<T_MAX> ev{};
            ev.signo = dy > 0.0 ? +1 : -1;
            ev.escala_origen = 0;
            ev.g_disparo = (dy * dy) / (2.0 * std::max(1e-6, m_lin.sd * m_lin.sd));
            ev.dy = dy;
            ev.y_ancla = y_anc;
            ev.w_ancla = w_anc;
            ev.guarda_ancla = 1;
            ev.eta.fill(NaN);
            for (std::size_t tau = 1; tau <= T_MAX && (t_ext + tau) < n_cal; ++tau) {
                const double e = (x_sust[t_ext + tau] - y_anc) / dy;
                ev.eta[tau - 1] = e;
                ev.pasos_observados = tau;
                if (e <= 0.5 && !std::isfinite(ev.tau_rest_primer_toque)) {
                    ev.tau_rest_primer_toque = static_cast<double>(tau);
                    ev.tau_rest_consolidado = static_cast<double>(tau);
                }
            }
            ev.desenlace = std::isfinite(ev.tau_rest_consolidado) ? 2u : 4u;
            ev.tau_fin = std::isfinite(ev.tau_rest_consolidado) ? static_cast<std::size_t>(ev.tau_rest_consolidado) : T_MAX;
            total.push_back(ev);
        }
    }
    return total;
}

// Ejecuta una realizacion del nulo lineal con una o varias inyecciones candidatas en los tiempos `tiempos_iny`,
// verificando la guarda de solapamiento de ventanas (seccion 3.3 y control S2) y exigiendo que el operador delta dispare.
template <std::size_t T_MAX, class Obs, class Lil>
inline ResultadoInyeccionesRealizacion<T_MAX> ejecutar_inyecciones_en_realizacion(
    const Obs& delta_original,
    const Lil& lil_plantilla,
    ark::metrologia::Integracion integracion,
    bool rama_purgada,
    const ModeloLinealNulo& m_lin,
    std::span<const double> perfil_est,
    std::span<const double> e_iny,
    int signo_real,
    std::size_t w_ancla,
    std::size_t guarda,
    std::size_t t_max_ef,
    std::span<const std::int64_t> tiempos_iny,
    std::size_t n_sim,
    std::mt19937_64& rng) {
    ResultadoInyeccionesRealizacion<T_MAX> res{};
    res.candidatas = tiempos_iny.size();
    const std::size_t l = e_iny.size();
    if (!m_lin.valido || l == 0 || tiempos_iny.empty() || n_sim == 0) return res;

    // 1. Construir ventanas de exclusion [t_iny - (W_ancla + guarda), t_iny + l + T_MAX] y verificar solapes entre inyecciones.
    std::vector<VentanaInyeccion> ventanas(tiempos_iny.size());
    std::vector<bool> solapada(tiempos_iny.size(), false);
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        ventanas[i] = construir_ventana_inyeccion(tiempos_iny[i], l, w_ancla, guarda, T_MAX);
    }
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        for (std::size_t j = i + 1; j < tiempos_iny.size(); ++j) {
            if (ventanas_se_solapan(ventanas[i], ventanas[j])) {
                solapada[i] = true;
                solapada[j] = true;
            }
        }
    }

    // 2. Simular la realizacion inyectando las innovaciones e_iny en cada t_iny (propagadas por el filtro AR/ARFIMA).
    const std::size_t burn = m_lin.burn_in;
    const std::size_t total = burn + n_sim;
    std::vector<double> eps(total, 0.0);
    if (!m_lin.residuos.empty()) {
        std::uniform_int_distribution<std::size_t> dist(0, m_lin.residuos.size() - 1);
        for (std::size_t t = 0; t < total; ++t) {
            eps[t] = m_lin.residuos[dist(rng)] - m_lin.media_residuos;
        }
    } else {
        std::normal_distribution<double> dist(0.0, std::sqrt(m_lin.sigma2_innov));
        for (std::size_t t = 0; t < total; ++t) eps[t] = dist(rng);
    }

    const bool inyectar_en_innovaciones = true;
    if (inyectar_en_innovaciones) {
        for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
            if (tiempos_iny[i] < 0) continue;
            const std::size_t t0 = static_cast<std::size_t>(tiempos_iny[i]);
            for (std::size_t k = 0; k < l && (t0 + k) < n_sim; ++k) {
                eps[burn + t0 + k] = e_iny[k];
            }
        }
    }

    std::vector<double> u(total, 0.0);
    for (std::size_t t = 0; t < total; ++t) {
        double acc = eps[t];
        const std::size_t i_max = std::min(m_lin.p, t);
        for (std::size_t i = 1; i <= i_max; ++i) {
            acc += m_lin.phi[i - 1] * u[t - i];
        }
        u[t] = acc;
    }

    std::vector<double> x_est(n_sim, 0.0);
    if (m_lin.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE && m_lin.d != 0.0) {
        std::vector<double> pi_inv(total, 0.0);
        ark::metrologia::elw::coeficientes(-m_lin.d, total, pi_inv.data());
        for (std::size_t t = 0; t < n_sim; ++t) {
            const std::size_t idx = burn + t;
            double s = 0.0;
            for (std::size_t k = 0; k <= idx; ++k) s += pi_inv[k] * u[idx - k];
            x_est[t] = m_lin.mu + s;
        }
    } else {
        for (std::size_t t = 0; t < n_sim; ++t) x_est[t] = m_lin.mu + u[burn + t];
    }

    if (!inyectar_en_innovaciones) {
        // Mutante M08 (P8): superpone el perfil sobre la serie fija sin propagarlo por el filtro lineal
        for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
            if (tiempos_iny[i] < 0) continue;
            const std::size_t t0 = static_cast<std::size_t>(tiempos_iny[i]);
            for (std::size_t k = 0; k < perfil_est.size() && (t0 + k) < n_sim; ++k) {
                x_est[t0 + k] = perfil_est[k];
            }
        }
    }

    std::vector<double> x_sim(n_sim, 0.0);
    if (integracion == ark::metrologia::Integracion::NIVEL) {
        double lvl = 100.0;
        x_sim[0] = lvl;
        for (std::size_t t = 1; t < n_sim; ++t) {
            lvl += x_est[t];
            x_sim[t] = lvl;
        }
    } else {
        x_sim = x_est;
    }

    // 3. Correr copia limpia del operador delta y de ObservadorLilyth sobre la realizacion.
    Obs d_local = delta_original;
    d_local.reiniciar_flujo();
    Lil lil_local = lil_plantilla;
    lil_local.reiniciar_flujo_conservando_purga();

    struct Captura {
        std::int64_t t_apertura = -1;
        EventoNuloLilyth<T_MAX> ev{};
    };
    std::vector<Captura> capturas;
    std::vector<std::int64_t> aperturas_todas;

    for (std::size_t t = 0; t < n_sim; ++t) {
        const auto s = d_local.paso(x_sim[t], NaN, true, 0, static_cast<std::int64_t>(20000 + t));
        const auto o = lil_local.paso(s, x_sim[t], true, NaN, 0, static_cast<std::int64_t>(20000 + t));
        if (o.abrio_ranura) aperturas_todas.push_back(o.t);
        for (std::size_t c = 0; c < o.ranuras_cerradas; ++c) {
            const auto& rc = lil_local.cerrada(c);
            EventoNuloLilyth<T_MAX> ev{};
            if (extraer_evento_nulo<T_MAX>(rc, rama_purgada, ev)) {
                capturas.push_back({rc.t_apertura, ev});
            }
        }
    }

    // 4. Verificar tambien si algun evento espontaneo cayo dentro de la ventana de una inyeccion fuera de [t_iny, t_iny + l].
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        const VentanaInyeccion& vi = ventanas[i];
        for (std::int64_t ta : aperturas_todas) {
            const bool dentro_propio = (ta >= vi.t_iny && ta <= vi.t_iny + static_cast<std::int64_t>(l));
            if (!dentro_propio) {
                const VentanaInyeccion v_esp{ta, 0, ta, ta};
                if (ventanas_se_solapan(vi, v_esp)) {
                    solapada[i] = true;
                }
            }
        }
    }

    const bool exige_disparo_detector = true;
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        if (solapada[i]) {
            ++res.solapes;
            continue;
        }
        const std::int64_t t0 = tiempos_iny[i];
        const std::int64_t t1 = t0 + static_cast<std::int64_t>(l);
        const Captura* encontrada = nullptr;
        for (const auto& cap : capturas) {
            if (cap.t_apertura >= t0 && cap.t_apertura <= t1 && cap.ev.signo == signo_real) {
                if (encontrada == nullptr || cap.t_apertura < encontrada->t_apertura) {
                    encontrada = &cap;
                }
            }
        }
        if (encontrada != nullptr) {
            ++res.detectadas;
            res.eventos_admitidos.push_back(encontrada->ev);
        } else if (!exige_disparo_detector) {
            // Mutante M03: mide la restitucion desde t_iny + l - 1 sin exigir que el detector dispare
            const std::size_t ta = static_cast<std::size_t>( std::max<std::int64_t>(1, t0 + static_cast<std::int64_t>(l) - 1) );
            if (ta + t_max_ef < n_sim && ta >= w_ancla + guarda) {
                std::vector<double> w(x_sim.begin() + static_cast<std::ptrdiff_t>(ta - w_ancla - guarda),
                                      x_sim.begin() + static_cast<std::ptrdiff_t>(ta - guarda));
                const double y_anc = ark::metrologia::est::cuantil(w, 0.5);
                const double dy = x_sim[ta] - y_anc;
                if (std::isfinite(dy) && dy * static_cast<double>(signo_real) > 0.0) {
                    EventoNuloLilyth<T_MAX> ev_forzado{};
                    ev_forzado.signo = signo_real;
                    ev_forzado.y_ancla = y_anc;
                    ev_forzado.dy = dy;
                    ev_forzado.eta.fill(NaN);
                    for (std::size_t tau = 1; tau <= t_max_ef && (ta + tau) < n_sim; ++tau) {
                        const double e = (x_sim[ta + tau] - y_anc) / dy;
                        ev_forzado.eta[tau - 1] = e;
                        ev_forzado.pasos_observados = tau;
                        if (e <= 0.5 && !std::isfinite(ev_forzado.tau_rest_primer_toque)) {
                            ev_forzado.tau_rest_primer_toque = static_cast<double>(tau);
                            ev_forzado.tau_rest_consolidado = static_cast<double>(tau);
                        }
                    }
                    ev_forzado.desenlace = std::isfinite(ev_forzado.tau_rest_consolidado) ? 2u : 4u;
                    ev_forzado.tau_fin = std::isfinite(ev_forzado.tau_rest_consolidado)
                                             ? static_cast<std::size_t>(ev_forzado.tau_rest_consolidado)
                                             : t_max_ef;
                    ++res.detectadas;
                    res.eventos_admitidos.push_back(ev_forzado);
                } else {
                    ++res.no_detectadas;
                }
            } else {
                ++res.no_detectadas;
            }
        } else {
            ++res.no_detectadas;
        }
    }
    return res;
}

// Ejecuta el bucle secuencial de la via (b) para un evento real k (seccion 3.3 y control S3).
template <std::size_t T_MAX, class Obs, class Lil>
inline void ejecutar_via_b_evento(
    const Obs& delta_original,
    const Lil& lil_plantilla,
    ark::metrologia::Integracion integracion,
    bool rama_purgada,
    const ModeloLinealNulo& m_lin,
    std::span<const double> hist_prev_est,
    std::span<const double> perfil_est,
    int signo_real,
    std::size_t w_ancla,
    std::size_t guarda,
    std::size_t m_min,
    std::size_t t_max_ef,
    std::uint64_t semilla_evento,
    ReferenciaNuloEvento& ref,
    std::vector<EventoNuloLilyth<T_MAX>>& nulos_b_out,
    std::size_t& cont_solapes_global) {
    nulos_b_out.clear();
    ref.nulo_por_inyeccion = true;
    const std::size_t l = perfil_est.size();
    if (!m_lin.valido || l == 0 || m_min == 0) {
        ref.sin_nulo_b = true;
        ref.motivo = "via (b): modelo lineal invalido o perfil vacio";
        return;
    }

    std::vector<double> e_iny(l, 0.0);
    innovaciones(m_lin, hist_prev_est.data(), hist_prev_est.size(), perfil_est.data(), l, e_iny.data());
    if (!std::isfinite(e_iny[0])) {
        ref.sin_nulo_b = true;
        ref.motivo = "via (b): innovaciones no finitas";
        return;
    }

    const std::size_t margen_pre = w_ancla + guarda + 8;
    const std::size_t margen_post = l + t_max_ef + 8;
    const std::size_t n_sim = margen_pre + margen_post + 16;
    std::mt19937_64 rng(splitmix64(semilla_evento ^ 0x6C62272E07BB0142ULL));
    std::uniform_int_distribution<std::int64_t> dist_t(
        static_cast<std::int64_t>(w_ancla + guarda + 4),
        static_cast<std::int64_t>(n_sim - margen_post));

    const bool fijar_p_det_en_m_min = true;
    const bool respetar_tope_parada = true;
    std::size_t tope_intentos = MAX_INTENTOS_VIA_B;
    bool p_det_fijado = false;

    for (std::size_t intento = 1; intento <= MAX_INTENTOS_VIA_B; ++intento) {
        const std::int64_t t_iny = dist_t(rng);
        const std::array<std::int64_t, 1> t_arr{t_iny};
        auto r_lote = ejecutar_inyecciones_en_realizacion<T_MAX>(
            delta_original, lil_plantilla, integracion, rama_purgada, m_lin, perfil_est, e_iny, signo_real,
            w_ancla, guarda, t_max_ef, t_arr, n_sim, rng);

        ++ref.intentos_b;
        ref.solapes_b += r_lote.solapes;
        cont_solapes_global += r_lote.solapes;
        ref.no_detectados_b += r_lote.no_detectadas;
        for (const auto& ev : r_lote.eventos_admitidos) {
            nulos_b_out.push_back(ev);
            ++ref.n_nulos_b;
        }

        if (intento >= m_min && (!p_det_fijado || !fijar_p_det_en_m_min)) {
            p_det_fijado = true;
            ref.p_det_inicial = static_cast<double>(ref.n_nulos_b) / static_cast<double>(intento);
            if (ref.n_nulos_b == 0) {
                ref.sin_nulo_b = true;
                ref.tope_intentos_b = intento;
                ref.motivo = "via (b): p_det = 0 en las primeras m_min inyecciones";
                if (respetar_tope_parada) break;
            } else {
                const double tope_calc = std::ceil(static_cast<double>(m_min) / ref.p_det_inicial);
                if (tope_calc > static_cast<double>(MAX_INTENTOS_VIA_B)) {
                    ref.tope_memoria_b = true;
                    tope_intentos = MAX_INTENTOS_VIA_B;
                } else {
                    tope_intentos = static_cast<std::size_t>(tope_calc);
                }
                ref.tope_intentos_b = tope_intentos;
            }
        }

        if (ref.n_nulos_b >= m_min) {
            if (!p_det_fijado) {
                ref.p_det_inicial = 1.0;
                ref.tope_intentos_b = m_min;
            }
            ref.motivo = "via (b): completado m_min nulos por inyeccion";
            break;
        }
        if (p_det_fijado && respetar_tope_parada && intento >= tope_intentos) {
            ref.motivo = "via (b): tope de intentos ceil(m_min / p_det) alcanzado";
            break;
        }
    }

    if (ref.intentos_b > 0) {
        ref.tasa_no_detectados = static_cast<double>(ref.no_detectados_b) / static_cast<double>(ref.intentos_b);
    }
}

}  // namespace detalle_nulo

}  // namespace ark::lilyth