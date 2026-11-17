// etapa5_lilyth.cpp -- Bateria de verificacion de la Etapa 5 del Operador de Lilyth (`fin_calor_lilyth.hpp` + `observador_lilyth.hpp`):
//   1. K1: Identidad bit a bit del ObservadorSaltos real con y sin Lilyth y cero reservas de memoria en `paso`.
//   2. Referencia independiente de DetectorGLRUnilateralCalor y de la calibracion de h con ARL_0,caliente = T_max / alfa.
//   3. Inicio del calor (signo == -1 en canal_temperatura.hpp), mu_calor de las l* ventanas posteriores a la alarma (E1.7 [DECIDIDO])
//      y rechazo estricto de h extrapolado (f.calibrado == false cuando buscar_h_con extrapola).
//   4. P9: Con el calor apagado en un t conocido, la tasa de fin declarado en falso antes de t en T_max es <= alfa, y el retraso
//      posterior a t se mide y se reporta.
//   5. K4b rehecho con calor detectado automaticamente por ObservadorSaltos + GestorFinCalorLilyth (sin fijar_calor_activo):
//      separa `tau_rest_primer_toque` (durante el calor) de `tau_rest_consolidado` (tras extinguirse el calor).
//   6. Regla multiescala: con varias escalas de temperatura calientes a la vez, el calor termina solo cuando terminan todas.
//
// Compilacion: g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wno-mismatched-new-delete -Werror -pthread -Iinclude etapa5_lilyth.cpp -o etapa5_lilyth

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/lilyth/fin_calor_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <random>
#include <span>
#include <vector>

namespace {

// ----------------------------------------------------------------------------- Contador global de operator new (Regla 6 / K10)
std::atomic<bool> g_contar_reservas{false};
std::atomic<std::size_t> g_reservas_ciclo_caliente{0};

}  // namespace

void* operator new(std::size_t sz) {
    if (g_contar_reservas.load(std::memory_order_relaxed)) {
        g_reservas_ciclo_caliente.fetch_add(1, std::memory_order_relaxed);
    }
    if (sz == 0) sz = 1;
    void* p = std::malloc(sz);
    if (!p) throw std::bad_alloc();
    return p;
}

void operator delete(void* p) noexcept {
    std::free(p);
}

void operator delete(void* p, std::size_t) noexcept {
    std::free(p);
}

namespace {

using ObsReal = ark::metrologia::ObservadorSaltos<24, 64, 120>;
using LilythReal = ark::lilyth::ObservadorLilyth<ObsReal, 8, 64>;
using ark::metrologia::NaN;

int g_fallas = 0;

void verificar(bool cond, const char* msg) {
    if (!cond) {
        std::cerr << "  [FALLA] " << msg << "\n";
        ++g_fallas;
    }
}

// Genera una serie superestadistica con equilibrio local gaussiano en bloques de `t_bloque` barras
// (Beck, Cohen y Swinney 2005): dentro de cada bloque de `t_bloque` barras la serie es gaussiana con sigma_b = exp(0.42 * xi_b),
// de modo que CanalTemperatura::calibrar identifica T* = t_bloque y kalibra el nulo de log_beta^ sin extrapolacion.
std::vector<double> generar_superestadistica_calma(std::size_t n, std::size_t t_bloque, std::uint64_t sem) {
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<double> x(n, 0.0);
    double sigma_b = 1.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (i % t_bloque == 0) sigma_b = std::exp(0.42 * nd(rng));
        x[i] = sigma_b * nd(rng);
    }
    return x;
}

// =============================================================================
// MOTOR DE REFERENCIA INDEPENDIENTE (Regla 3: sin compartir codigo con la cabecera)
// =============================================================================

struct PasoUnilateralRef {
    bool alarma = false;
    double g = 0.0;
};

template <std::size_t L = 64>
struct GLRUnilateralRef {
    std::vector<double> z_hist;
    double mu = 0.0, s0 = 1.0, tau = 1.0, h = 1.0;

    void configurar(double m, double s, double t, double thr) {
        mu = m;
        s0 = s;
        tau = t;
        h = thr;
        z_hist.clear();
    }

    PasoUnilateralRef paso(double lb) {
        PasoUnilateralRef out{};
        if (!std::isfinite(lb)) {
            z_hist.clear();
            return out;
        }
        z_hist.push_back((lb - mu) / s0);
        if (z_hist.size() > L) z_hist.erase(z_hist.begin());

        double suma = 0.0, gmax = 0.0;
        for (std::size_t l = 1; l <= z_hist.size(); ++l) {
            suma += z_hist[z_hist.size() - l];
            if (suma > 0.0) {
                const double g = (suma * suma) / (2.0 * static_cast<double>(l) * tau);
                if (g > gmax) gmax = g;
            }
        }
        out.g = gmax;
        if (gmax >= h) {
            out.alarma = true;
            z_hist.clear();
        }
        return out;
    }
};

// =============================================================================
// 1. K1: IDENTIDAD DEL OPERADOR DELTA Y CERO RESERVAS EN PASO (Reglas 2 y 6)
// =============================================================================
void prueba_k1_identidad_y_cero_reservas() {
    std::cout << "[1/6] K1: Identidad bit a bit del ObservadorSaltos real con/sin Lilyth y 0 reservas en paso...\n";
    const std::size_t N_CAL = 16000;
    const std::size_t N_FLUJO = 1200;
    std::mt19937_64 rng(2026110651ULL);
    std::normal_distribution<double> nd(0.0, 1.0);

    const auto x_cal = generar_superestadistica_calma(N_CAL, 8, 2026110651ULL);

    ark::metrologia::ConfiguracionObservador cfg_d{};
    cfg_d.integracion = ark::metrologia::Integracion::FLUJO;
    cfg_d.escalas_max = 1;
    cfg_d.medir_tau50 = false;
    cfg_d.arl0_pasos_nativos = 250.0;
    cfg_d.semilla = 5001;

    auto delta_solo = std::make_unique<ObsReal>();
    auto delta_con_lilyth = std::make_unique<ObsReal>();
    delta_solo->calibrar(cfg_d, x_cal);
    delta_con_lilyth->calibrar(cfg_d, x_cal);
    verificar(delta_solo->escalas() > 0 && delta_con_lilyth->escalas() > 0, "K1: ambos ObservadorSaltos deben calibrar");

    ark::lilyth::ConfiguracionLilyth cfg_l{};
    cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::CRUDA;
    cfg_l.b_max_nulo = 2;
    cfg_l.hilos_nulo = 2;

    auto lil = std::make_unique<LilythReal>();
    verificar(lil->calibrar(*delta_con_lilyth, cfg_l, cfg_d, x_cal), "K1: Lilyth debe calibrar");

    std::vector<double> x_flujo(N_FLUJO, 0.0);
    for (std::size_t t = 0; t < N_FLUJO; ++t) {
        const double sd = (t >= 200 && t < 360) ? 4.5 : 1.0;
        x_flujo[t] = sd * nd(rng);
        if (t == 280) x_flujo[t] += 9.0;
    }

    g_reservas_ciclo_caliente.store(0);
    g_contar_reservas.store(true);
    for (std::size_t t = 0; t < N_FLUJO; ++t) {
        const auto s1 = delta_solo->paso(x_flujo[t]);
        const auto s2 = delta_con_lilyth->paso(x_flujo[t]);
        (void)lil->paso(s2, x_flujo[t], true);
        bool eq = (s1.t == s2.t && s1.n_eventos == s2.n_eventos && s1.evento_abierto == s2.evento_abierto && s1.niveles == s2.niveles);
        for (std::size_t j = 0; eq && j < s1.niveles; ++j) {
            eq = (s1.nivel[j].alarma.hay == s2.nivel[j].alarma.hay &&
                  s1.temp[j].alarma.hay == s2.temp[j].alarma.hay &&
                  ((std::isnan(s1.nivel[j].g) && std::isnan(s2.nivel[j].g)) || s1.nivel[j].g == s2.nivel[j].g));
        }
        if (!eq) {
            verificar(false, "K1: Salida de ObservadorSaltos difiere bit a bit");
            break;
        }
    }
    g_contar_reservas.store(false);
    verificar(g_reservas_ciclo_caliente.load() == 0, "K10: Cero llamadas a operator new en el ciclo caliente paso()");
}

// =============================================================================
// 2. REFERENCIA INDEPENDIENTE, SIGNO DE TEMPERATURA (-1), mu_calor (l* VENTANAS) Y RECHAZO DE h EXTRAPOLADO
// =============================================================================
void prueba_referencia_mu_calor_y_no_extrapolacion() {
    std::cout << "[2/6] Referencia independiente, signo de inicio (-1), mu_calor en l* ventanas y rechazo de h extrapolado...\n";

    // 2a. Casamiento bit a bit del DetectorGLRUnilateralCalor contra el motor de referencia independiente
    {
        ark::lilyth::DetectorGLRUnilateralCalor<64> d_prod;
        GLRUnilateralRef<64> d_ref;
        const double mu_c = -2.4, s0 = 0.55, tau = 1.25, h = 4.2;
        d_prod.configurar(mu_c, s0, tau, h, false);
        d_ref.configurar(mu_c, s0, tau, h);

        std::mt19937_64 rng(2026110652ULL);
        std::normal_distribution<double> nd(mu_c, s0);
        for (std::size_t w = 0; w < 300; ++w) {
            const double val = (w >= 120 && w < 140) ? (mu_c + 1.8 * s0) : nd(rng);
            const auto p1 = d_prod.paso(val);
            const auto p2 = d_ref.paso(val);
            if (p1.alarma != p2.alarma || std::abs(p1.g - p2.g) > 1e-14) {
                verificar(false, "Referencia GLR unilateral: discrepancia contra motor independiente");
                break;
            }
        }
    }

    // 2b. Calibracion con historia suficiente (sin extrapolar) vs historia corta (extrapolada -> calibrado == false)
    {
        std::vector<std::vector<double>> segs_larga(1, generar_superestadistica_calma(20000, 8, 20261106));

        ark::metrologia::CanalTemperatura<64> ct_larga;
        ark::metrologia::ConfiguracionTemperatura cfg_t{};
        cfg_t.arl0_pasos = 2000.0;
        cfg_t.alfa = 0.05;
        ct_larga.calibrar(0, 1, segs_larga, cfg_t);

        const auto f_ok = ark::lilyth::detalle_fin_calor::calibrar_escala_fin_calor<64>(
            0, 1, ct_larga.T(), ct_larga.ficha_glr(), segs_larga, 5, 0.05, 8);
        verificar(std::abs(f_ok.arl0_caliente - (5.0 / 0.05)) <= 1e-12,
                  "Calibracion h: arl0_caliente == T_max / alfa == 100.0");
        verificar(std::isfinite(f_ok.h) && !f_ok.h_extrapolado && f_ok.calibrado,
                  "Calibracion h (historia suficiente): h finito, sin extrapolar y calibrado == true");

        // Historia corta: 2400 barras -> 300 ventanas de T* = 8 < 8 * 100 = 800 -> buscar_h_con extrapola -> calibrado == false
        std::vector<std::vector<double>> segs_corta(1, std::vector<double>(segs_larga[0].begin(), segs_larga[0].begin() + 2400));

        const auto f_ext = ark::lilyth::detalle_fin_calor::calibrar_escala_fin_calor<64>(
            0, 1, ct_larga.T(), ct_larga.ficha_glr(), segs_corta, 5, 0.05, 8);
        verificar(f_ext.h_extrapolado, "Calibracion h (historia corta): buscar_h_con debe marcar h_extrapolado == true");
        verificar(!f_ext.calibrado, "Calibracion h (historia corta): prohibido usar h extrapolado (calibrado debe ser false)");

        // Verificar que una escala con h extrapolado NO activa ni usa su detector no calibrado
        ark::lilyth::GestorFinCalorLilyth<24, 64> gestor_ext;
        gestor_ext.calibrar_desde_segmentos(0, 1, ct_larga.T(), ct_larga.ficha_glr(), segs_corta, 5, 0.05, 8);
        std::array<ark::metrologia::SalidaTemperatura, 1> st_al{};
        st_al[0].alarma.hay = true;
        st_al[0].alarma.signo = -1;
        st_al[0].alarma.l = 3;
        gestor_ext.paso(st_al, 1, 100);
        verificar(!gestor_ext.calor_activo() && gestor_ext.contadores().alarmas_ignoradas_no_calibrado == 1,
                  "Rechazo de h extrapolado: el fin del calor no calibrado no se usa");

        // 2c. Verificar signo == -1 (inicio de calor) vs signo == +1 (enfriamiento no abre calor), y mu_calor en las l* ventanas
        ark::lilyth::GestorFinCalorLilyth<24, 64> gestor_ok;
        gestor_ok.calibrar_desde_segmentos(0, 1, ct_larga.T(), ct_larga.ficha_glr(), segs_larga, 5, 0.05, 8);

        // Alarma con signo == +1 (beta sube = se enfria) NO debe iniciar episodio de calor
        st_al[0].alarma.hay = true;
        st_al[0].alarma.signo = +1;
        st_al[0].alarma.l = 3;
        gestor_ok.paso(st_al, 1, 10);
        verificar(!gestor_ok.calor_activo(), "Inicio del calor: signo +1 (enfriamiento) no debe activar calor");

        // Alarma con signo == -1 (beta baja = sube temperatura) con l* = 3 y amplitud_alarma = -4.50 (distinta de las ventanas posteriores)
        st_al[0].alarma.hay = true;
        st_al[0].alarma.signo = -1;
        st_al[0].alarma.l = 3;
        st_al[0].alarma.amplitud = -4.50;
        gestor_ok.paso(st_al, 1, 20);
        verificar(gestor_ok.calor_activo() && !gestor_ok.estado(0).mu_calor_listo,
                  "Inicio del calor: signo -1 activa calor y espera l* = 3 ventanas posteriores");

        // Enviar las 3 ventanas posteriores a la alarma: log_beta = {-2.10, -2.30, -2.20} -> media exacta = -2.20
        const std::array<double, 3> lb_post{-2.10, -2.30, -2.20};
        for (std::size_t k = 0; k < 3; ++k) {
            std::array<ark::metrologia::SalidaTemperatura, 1> st_w{};
            st_w[0].emitida = true;
            st_w[0].log_beta = lb_post[k];
            gestor_ok.paso(st_w, 1, static_cast<std::int64_t>(21 + k));
            if (k < 2) {
                verificar(gestor_ok.calor_activo() && !gestor_ok.estado(0).mu_calor_listo,
                          "mu_calor: mientras no hay l* ventanas cerradas, mu_calor_listo es false y el calor sigue activo");
            }
        }
        verificar(gestor_ok.calor_activo() && gestor_ok.estado(0).mu_calor_listo,
                  "mu_calor: al cerrar la l*-esima ventana posterior, mu_calor_listo pasa a true y el calor sigue activo");
        verificar(std::abs(gestor_ok.estado(0).mu_calor - (-2.20)) <= 1e-12,
                  "mu_calor [DECIDIDO]: debe ser exactamente la media (-2.20) de las l* ventanas posteriores, nunca de la amplitud de la alarma");
    }
}

// =============================================================================
// 3. P9: CALOR APAGADO EN UN t CONOCIDO -> TASA DE FIN EN FALSO ANTES DE t <= ALFA EN T_max Y RETRASO MEDIDO
// =============================================================================
void prueba_p9_tasa_falso_fin_y_retraso() {
    std::cout << "[3/6] P9: Calor apagado en t conocido -> tasa de fin en falso antes de t <= alfa en T_max y medicion del retraso...\n";
    const double alfa = 0.05;
    const std::size_t t_max_ventanas = 8; // horizonte T_max en ventanas del nulo caliente

    std::mt19937_64 rng(2026110654ULL);
    std::normal_distribution<double> nd(0.0, 1.0);

    // Calibrar sobre historia superestadistica de calma con 24000 barras (3000 ventanas de T* = 8 >= 8 * 160 = 1280)
    std::vector<std::vector<double>> segs_calma(1, generar_superestadistica_calma(24000, 8, 20261106));

    ark::metrologia::CanalTemperatura<64> ct;
    ark::metrologia::ConfiguracionTemperatura cfg_t{};
    cfg_t.arl0_pasos = 2000.0;
    cfg_t.alfa = alfa;
    ct.calibrar(0, 1, segs_calma, cfg_t);

    ark::lilyth::GestorFinCalorLilyth<24, 64> gestor;
    gestor.calibrar_desde_segmentos(0, 1, ct.T(), ct.ficha_glr(), segs_calma, t_max_ventanas, alfa, 8);
    const auto& f = gestor.ficha(0);
    verificar(f.calibrado && !f.h_extrapolado, "P9: escala 0 debe estar calibrada sin extrapolar");

    const std::size_t R = 400;
    const std::size_t l_star = 4;
    const std::size_t T_ven = ct.T();
    const double sigma_calor = 4.0;

    std::size_t falsos_fines_antes_de_t = 0;
    std::vector<double> retrasos_ventanas;
    retrasos_ventanas.reserve(R);

    for (std::size_t r = 0; r < R; ++r) {
        gestor.reiniciar_flujo();
        std::int64_t paso_t = 100;

        // Disparo de inicio de calor en t = 100
        std::array<ark::metrologia::SalidaTemperatura, 1> st_ini{};
        st_ini[0].alarma.hay = true;
        st_ini[0].alarma.signo = -1;
        st_ini[0].alarma.l = l_star;
        st_ini[0].alarma.amplitud = -2.8;
        gestor.paso(st_ini, 1, paso_t);

        bool falso_fin = false;
        // Fase caliente: l_star ventanas para mu_calor + t_max_ventanas ventanas antes de t_apagado
        // Cada ventana tiene la misma marginal del nulo de log_beta desplazada a sigma_calor
        for (std::size_t w = 0; w < l_star + t_max_ventanas; ++w) {
            const double sigma_b = sigma_calor * std::exp(0.42 * nd(rng));
            std::vector<double> buf(T_ven, 0.0);
            for (double& x : buf) x = sigma_b * nd(rng);
            const double s2 = ark::metrologia::detalle_temp::estadisticos(buf.data(), T_ven).s2;
            std::array<ark::metrologia::SalidaTemperatura, 1> st{};
            st[0].emitida = true;
            st[0].log_beta = -std::log(s2);
            paso_t += static_cast<std::int64_t>(T_ven);
            gestor.paso(st, 1, paso_t);
            if (!gestor.calor_activo()) {
                falso_fin = true;
                break;
            }
        }
        if (falso_fin) {
            ++falsos_fines_antes_de_t;
            continue;
        }

        // En este instante exacto (t_apagado) se apaga el calor: sigma vuelve al regimen de calma (mu_0,beta)
        for (std::size_t w_post = 1; w_post <= 24; ++w_post) {
            const double sigma_b = std::exp(0.42 * nd(rng));
            std::vector<double> buf(T_ven, 0.0);
            for (double& x : buf) x = sigma_b * nd(rng);
            const double s2 = ark::metrologia::detalle_temp::estadisticos(buf.data(), T_ven).s2;
            std::array<ark::metrologia::SalidaTemperatura, 1> st{};
            st[0].emitida = true;
            st[0].log_beta = -std::log(s2);
            paso_t += static_cast<std::int64_t>(T_ven);
            gestor.paso(st, 1, paso_t);
            if (!gestor.calor_activo()) {
                retrasos_ventanas.push_back(static_cast<double>(w_post));
                break;
            }
        }
    }

    const double tasa_falso_fin = static_cast<double>(falsos_fines_antes_de_t) / static_cast<double>(R);
    const double med_retraso = ark::metrologia::est::cuantil(retrasos_ventanas, 0.5);
    const double iqr_retraso = ark::metrologia::est::cuantil(retrasos_ventanas, 0.75) -
                               ark::metrologia::est::cuantil(retrasos_ventanas, 0.25);

    std::cout << "  [P9] tasa_falso_fin_en_T_max = " << std::fixed << std::setprecision(4) << tasa_falso_fin
              << " (alfa = " << alfa << "), retraso_mediano = " << std::setprecision(1) << med_retraso
              << " ventanas (IQR = " << iqr_retraso << " ventanas), detectados_post_apagado = "
              << retrasos_ventanas.size() << "/" << (R - falsos_fines_antes_de_t) << "\n";

    verificar(tasa_falso_fin <= alfa + 0.027,
              "P9: la tasa de fin del calor declarado en falso antes de t_apagado en T_max debe ser <= alfa (+ error MC)");
    verificar(retrasos_ventanas.size() >= (R - falsos_fines_antes_de_t) * 9 / 10 && med_retraso <= 4.0,
              "P9: al apagarse el calor en t_apagado, el detector declara fin del calor con retraso mediano <= 4 ventanas");
}

// =============================================================================
// 4. K4b REHECHO CON CALOR DETECTADO POR OBSERVADORSALTOS REAL (SIN fijar_calor_activo)
// =============================================================================
void prueba_k4b_calor_detectado_real() {
    std::cout << "[4/6] K4b rehecho con calor detectado automaticamente en ObservadorSaltos + ObservadorLilyth...\n";
    const std::size_t N_CAL = 20000;
    std::mt19937_64 rng(2026110655ULL);
    std::normal_distribution<double> nd(0.0, 1.0);

    const auto x_cal = generar_superestadistica_calma(N_CAL, 8, 20261106);

    ark::metrologia::ConfiguracionObservador cfg_d{};
    cfg_d.integracion = ark::metrologia::Integracion::FLUJO;
    cfg_d.escalas_max = 1;
    cfg_d.medir_tau50 = false;
    cfg_d.arl0_pasos_nativos = 250.0;
    cfg_d.semilla = 5005;

    auto delta = std::make_unique<ObsReal>();
    delta->calibrar(cfg_d, x_cal);
    verificar(delta->escalas() == 1, "K4b: ObservadorSaltos debe calibrar escala 0");

    ark::lilyth::ConfiguracionLilyth cfg_l{};
    cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::CRUDA;
    cfg_l.t_max_efectivo = 24; // ARL_0,caliente = 20 / 0.05 = 400 ventanas <= 4500 / 8 -> calibra sin extrapolar y cubre tau = 15

    auto lil = std::make_unique<LilythReal>();
    verificar(lil->configurar(cfg_l, cfg_d, *delta, x_cal), "K4b: ObservadorLilyth debe configurar con fin_calor calibrado");
    verificar(lil->gestor_calor().ficha(0).calibrado && !lil->gestor_calor().ficha(0).h_extrapolado,
              "K4b: fin_calor de escala 0 debe estar calibrado sin extrapolar");

    delta->reiniciar_flujo();
    lil->reiniciar_flujo();

    bool vio_calor_activo = false;
    bool vio_fin_calor = false;
    bool abrio_en_calor = false;

    for (std::size_t t = 0; t < 300; ++t) {
        double x = 0.0;
        if (t < 160) {
            x = 0.20 * nd(rng);
        } else if (t < 200) {
            // 5 ventanas de T* = 8 barras ([160..199]) con varianza alta (s^2 ~ 9.0) y suma exactamente 0 en cada par:
            // dispara s.temp[0].alarma (signo == -1) y completa las l* ventanas calientes sin disparar el canal de nivel
            x = (t % 2 == 0) ? 2.8 : -2.8;
        } else if (t == 200) {
            // Disparo de nivel positivo mientras el calor detectado ya tiene mu_calor_listo == true y sigue activo
            x = 14.0;
        } else if (t == 201) {
            // tau = 1: cruza eta <= 0.5 ESTANDO CALIENTE -> registra tau_rest_primer_toque = 1, pero NO consolidado
            x = 0.0;
        } else if (t >= 202 && t < 208) {
            // Mantiene el nivel por encima de 0.5 * dY y la ventana [200..207] caliente hasta cerrar t = 207
            x = 10.0 + ((t % 2 == 0) ? 3.5 : -3.5);
        } else {
            // Ventana [208..215]: el calor se apaga (s^2 ~ 0.0001 -> log_beta = +9.2 >> mu_calor) y el nivel vuelve a ~0.0 (eta <= 0.5).
            // Al cerrar la ventana en t = 215 (tau = 15 <= t_max_efectivo = 24), fin_calor declara fin del calor y consolida en tau = 15!
            x = (t % 2 == 0) ? 0.01 : -0.01;
        }

        const auto s = delta->paso(x);
        const auto o = lil->paso(s, x, true);
        if (lil->calor_activo()) vio_calor_activo = true;
        if (o.abrio_ranura && lil->calor_activo()) abrio_en_calor = true;
        if (vio_calor_activo && t >= 202 && !lil->calor_activo()) vio_fin_calor = true;
    }

    verificar(vio_calor_activo, "K4b: el canal de temperatura activo automaticamente calor_activo() sin usar fijar_calor_activo");
    verificar(abrio_en_calor, "K4b: la ranura de nivel abrio mientras el calor detectado estaba activo");
    verificar(vio_fin_calor, "K4b: el detector unilateral declaro automaticamente el fin del calor al enfriarse la serie");

    bool encontro_ranura_k4b = false;
    for (const auto& r : lil->cerradas_historico()) {
        if (std::isfinite(r.cruda.tau_rest_primer_toque) && std::isfinite(r.cruda.tau_rest_consolidado) &&
            r.cruda.tau_rest_consolidado > r.cruda.tau_rest_primer_toque) {
            encontro_ranura_k4b = true;
            std::cout << "  [K4b] Ranura id=" << r.id_evento << " t_a=" << r.t_apertura
                      << ": tau_rest_primer_toque = " << r.cruda.tau_rest_primer_toque
                      << ", tau_rest_consolidado = " << r.cruda.tau_rest_consolidado << "\n";
            break;
        }
    }
    verificar(encontro_ranura_k4b,
              "K4b: con calor detectado automaticamente, el cruce durante el calor se registra en primer_toque y no en consolidado");
}

// =============================================================================
// 5. MULTIESCALA: EL CALOR TERMINA SOLO CUANDO TERMINAN TODAS LAS ESCALAS CALIENTES
// =============================================================================
void prueba_multiescala_terminan_todas() {
    std::cout << "[5/6] Regla multiescala: con varias escalas calientes a la vez, el calor termina cuando terminan todas...\n";
    std::vector<std::vector<double>> segs0(1, generar_superestadistica_calma(20000, 8, 20261106));
    std::vector<std::vector<double>> segs1(1, generar_superestadistica_calma(20000, 8, 20261107));

    ark::metrologia::CanalTemperatura<64> ct0, ct1;
    ark::metrologia::ConfiguracionTemperatura cfg_t{};
    cfg_t.arl0_pasos = 2000.0;
    cfg_t.alfa = 0.05;
    ct0.calibrar(0, 1, segs0, cfg_t);
    ct1.calibrar(1, 2, segs1, cfg_t);

    ark::lilyth::GestorFinCalorLilyth<24, 64> gestor;
    gestor.calibrar_desde_segmentos(0, 1, ct0.T(), ct0.ficha_glr(), segs0, 5, 0.05, 8);
    gestor.calibrar_desde_segmentos(1, 2, ct1.T(), ct1.ficha_glr(), segs1, 5, 0.05, 8);
    verificar(gestor.ficha(0).calibrado && gestor.ficha(1).calibrado,
              "Multiescala: ambas escalas (0 y 1) deben estar calibradas");

    // 1. Disparo de calor simultaneo en las escalas 0 y 1 (l* = 2 en ambas)
    std::array<ark::metrologia::SalidaTemperatura, 2> st{};
    st[0].alarma.hay = true; st[0].alarma.signo = -1; st[0].alarma.l = 2; st[0].alarma.amplitud = -2.5;
    st[1].alarma.hay = true; st[1].alarma.signo = -1; st[1].alarma.l = 2; st[1].alarma.amplitud = -2.5;
    gestor.paso(st, 2, 10);
    verificar(gestor.calor_activo() && gestor.escalas_calientes() == 2,
              "Multiescala: deben quedar 2 escalas calientes activas");

    // 2. Completar las l* = 2 ventanas posteriores a la alarma en ambas escalas con log_beta = -2.5
    for (int k = 1; k <= 2; ++k) {
        std::array<ark::metrologia::SalidaTemperatura, 2> sw{};
        sw[0].emitida = true; sw[0].log_beta = -2.5;
        sw[1].emitida = true; sw[1].log_beta = -2.5;
        gestor.paso(sw, 2, 10 + k);
    }
    verificar(gestor.estado(0).mu_calor_listo && gestor.estado(1).mu_calor_listo,
              "Multiescala: mu_calor_listo en ambas escalas");

    // 3. Enfriar SOLO la escala 0 (log_beta = +1.5 >> mu_0,beta), manteniendo caliente la escala 1 (log_beta = -2.5)
    for (int k = 1; k <= 5; ++k) {
        std::array<ark::metrologia::SalidaTemperatura, 2> sw{};
        sw[0].emitida = true; sw[0].log_beta = +1.5;
        sw[1].emitida = true; sw[1].log_beta = -2.5;
        gestor.paso(sw, 2, 20 + k);
    }
    verificar(!gestor.estado(0).caliente && gestor.estado(1).caliente,
              "Multiescala: la escala 0 ya enfrio pero la escala 1 sigue caliente");
    verificar(gestor.calor_activo() && gestor.escalas_calientes() == 1,
              "Multiescala: el calor global NO debe terminar cuando termina solo la primera escala");

    // 4. Ahora enfriar tambien la escala 1 (log_beta = +1.5) -> el calor global termina
    for (int k = 1; k <= 5; ++k) {
        std::array<ark::metrologia::SalidaTemperatura, 2> sw{};
        sw[0].emitida = false;
        sw[1].emitida = true; sw[1].log_beta = +1.5;
        gestor.paso(sw, 2, 30 + k);
    }
    verificar(!gestor.estado(1).caliente && !gestor.calor_activo() && gestor.escalas_calientes() == 0,
              "Multiescala: el calor global termina cuando terminan todas las escalas calientes");
}

// =============================================================================
// 6. TESTIGO DETERMINISTA DE MUTANTES DE LA ETAPA 5
// =============================================================================
void prueba_testigo_mutantes_etapa5() {
    std::cout << "[6/6] Control testigo determinista de los mutantes de la Etapa 5...\n";
    std::vector<std::vector<double>> segs(1, generar_superestadistica_calma(20000, 8, 20261106));

    ark::metrologia::CanalTemperatura<64> ct;
    ark::metrologia::ConfiguracionTemperatura cfg_t{};
    cfg_t.arl0_pasos = 2000.0;
    cfg_t.alfa = 0.05;
    ct.calibrar(0, 1, segs, cfg_t);

    ark::lilyth::GestorFinCalorLilyth<24, 64> gestor;
    gestor.calibrar_desde_segmentos(0, 1, ct.T(), ct.ficha_glr(), segs, 6, 0.05, 8);
    verificar(std::abs(gestor.ficha(0).arl0_caliente - 120.0) <= 1e-12,
              "Testigo M03: arl0_caliente == T_max / alfa == 120.0");

    std::array<ark::metrologia::SalidaTemperatura, 1> st_al{};
    st_al[0].alarma.hay = true;
    st_al[0].alarma.signo = -1;
    st_al[0].alarma.l = 2;
    st_al[0].alarma.amplitud = -5.0;

    // Enviar mu_calor = -0.45 en las l* = 2 ventanas posteriores, y luego una ventana en -0.40:
    // Respecto de mu_calor = -0.45, z = +0.05 / s0 (G^+ ~ 0.001 << h, NO dispara),
    // pero bajo F4 invertida (M01, z >= -1.0 respecto de mu0 = +0.16), z = (-0.40 - 0.16)/0.98 = -0.57 >= -1.0 (SI dispara en falso!).
    gestor.paso(st_al, 1, 10);
    for (int k = 1; k <= 2; ++k) {
        std::array<ark::metrologia::SalidaTemperatura, 1> sw{};
        sw[0].emitida = true;
        sw[0].log_beta = -0.45;
        gestor.paso(sw, 1, 10 + k);
    }
    {
        std::array<ark::metrologia::SalidaTemperatura, 1> sw{};
        sw[0].emitida = true;
        sw[0].log_beta = -0.40;
        gestor.paso(sw, 1, 15);
        verificar(gestor.calor_activo(),
                  "Testigo M01 (F4): una fluctuacion leve dentro del calor (-0.45 -> -0.40) no puede declarar fin del calor");
    }
}

}  // namespace

int main(int argc, char** argv) {
    const bool modo_rapido = (argc > 1 && std::strcmp(argv[1], "--rapido") == 0);
    std::cout << "====================================================================\n";
    std::cout << " BATERIA ETAPA 5 -- OPERADOR DE LILYTH (FIN DEL CALOR, C++20)\n";
    std::cout << "====================================================================\n";

    prueba_referencia_mu_calor_y_no_extrapolacion();
    prueba_multiescala_terminan_todas();
    prueba_testigo_mutantes_etapa5();

    if (!modo_rapido) {
        prueba_k1_identidad_y_cero_reservas();
        prueba_p9_tasa_falso_fin_y_retraso();
        prueba_k4b_calor_detectado_real();
    }

    std::cout << "====================================================================\n";
    if (g_fallas == 0) {
        std::cout << " RESULTADO: TODOS LOS CONTROLES DE LA ETAPA 5 PASARON (0 fallas)\n";
        std::cout << "====================================================================\n";
        return 0;
    }
    std::cerr << " RESULTADO: " << g_fallas << " FALLA(S) DETECTADA(S)\n";
    std::cerr << "====================================================================\n";
    return 1;
}