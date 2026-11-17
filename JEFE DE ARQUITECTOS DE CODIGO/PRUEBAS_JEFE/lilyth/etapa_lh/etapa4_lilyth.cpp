// etapa4_lilyth.cpp -- Bateria de verificacion de la Etapa 4 del Operador de Lilyth (`armonico_lilyth.hpp` + `observador_lilyth.hpp`):
//   1. K1: Identidad bit a bit del ObservadorSaltos real con y sin Lilyth y cero reservas de memoria en `paso`.
//   2. Prony exacto sobre una recurrencia de orden 2 a 1e-12 (desde tau = 2 con eta(0) = 1) y motor de referencia independiente.
//   3. Periodo de un coseno amortiguado conocido (T = 8.0 y T = 3.2 <= 4 barras, verificando que NO se recorta el periodo, E1.8).
//   4. P5: AR(2) subamortiguado en el pipeline real -> raices complejas en eta_barra y rechazo de sobrepaso contra ambos nulos (iid y lineal) tras Holm.
//   5. P5b: AR(1) con phi < 0 (rebote bid-ask) -> raiz real negativa (discriminante >= 0), nunca oscilador.
//   6. Tasa de etiqueta ARMONICO_SOSPECHOSO en un AR(1) <= alfa (con su tasa en el nulo al lado, sin vetar ningun evento).
//   7. Fase reescalada (PRONY_FASE_REESCALADA, E1.20): exclusion y conteo exacto de no restituidos e interpolacion lineal sin extrapolar.
//   8. Familia de 4 hipotesis en ark::holm (E1.16) contra referencia independiente y falla dura ante cualquier p = NaN o fuera de [0, 1].
//
// Compilacion: g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wno-mismatched-new-delete -Werror -pthread -Iinclude etapa4_lilyth.cpp -o etapa4_lilyth

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/lilyth/armonico_lilyth.hpp"
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

// =============================================================================
// MOTORES DE REFERENCIA INDEPENDIENTES (Regla 3: sin compartir codigo con la cabecera)
// =============================================================================

struct PronyRef {
    bool valido = false;
    double a1 = NaN;
    double a2 = NaN;
    double disc = NaN;
    bool complejas = false;
    bool real_neg = false;
    double periodo = NaN;
    double min_eta = NaN;
    double rmse = NaN;
};

PronyRef prony_independiente_ref(const std::vector<double>& eta_1_a_H, double eta0 = 1.0) {
    PronyRef r{};
    const std::size_t H = eta_1_a_H.size();
    if (H < 3) return r;

    for (std::size_t tau = 2; tau <= H; ++tau) {
        const double v = eta_1_a_H[tau - 1];
        if (std::isfinite(v) && (!std::isfinite(r.min_eta) || v < r.min_eta)) {
            r.min_eta = v;
        }
    }

    auto at_tau = [&](std::size_t tau) -> double {
        return (tau == 0) ? eta0 : eta_1_a_H[tau - 1];
    };

    double s11 = 0.0, s12 = 0.0, s22 = 0.0, sy1 = 0.0, sy2 = 0.0;
    std::size_t n = 0;
    for (std::size_t tau = 2; tau <= H; ++tau) {
        const double y = at_tau(tau);
        const double u1 = at_tau(tau - 1);
        const double u2 = at_tau(tau - 2);
        if (!std::isfinite(y) || !std::isfinite(u1) || !std::isfinite(u2)) continue;
        s11 += u1 * u1;
        s12 += u1 * u2;
        s22 += u2 * u2;
        sy1 += u1 * y;
        sy2 += u2 * y;
        ++n;
    }
    if (n < 2) return r;
    const double det = s11 * s22 - s12 * s12;
    if (std::abs(det) <= 1e-24 * (s11 + s22) * (s11 + s22)) return r;

    r.a1 = (sy1 * s22 - sy2 * s12) / det;
    r.a2 = (s11 * sy2 - s12 * sy1) / det;
    r.disc = r.a1 * r.a1 + 4.0 * r.a2;

    double sse = 0.0;
    for (std::size_t tau = 2; tau <= H; ++tau) {
        const double e = at_tau(tau) - (r.a1 * at_tau(tau - 1) + r.a2 * at_tau(tau - 2));
        sse += e * e;
    }
    r.rmse = std::sqrt(sse / static_cast<double>(n));

    if (r.disc < 0.0 && r.a2 < 0.0) {
        r.complejas = true;
        const double mod_r = std::sqrt(-r.a2);
        const double cw = std::clamp(r.a1 / (2.0 * mod_r), -1.0, 1.0);
        const double w = std::acos(cw);
        if (w > 0.0) r.periodo = (2.0 * 3.14159265358979323846) / w;
    } else if (r.disc >= 0.0) {
        const double sd = std::sqrt(r.disc);
        const double lam1 = 0.5 * (r.a1 + sd);
        const double lam2 = 0.5 * (r.a1 - sd);
        const double dom = (std::abs(lam1) >= std::abs(lam2)) ? lam1 : lam2;
        r.real_neg = (dom < -1e-6);
    }
    r.valido = true;
    return r;
}

std::array<double, 4> holm_cuatro_independiente_ref(const std::array<double, 4>& p) {
    for (double x : p) {
        if (!std::isfinite(x) || x < 0.0 || x > 1.0) {
            return {NaN, NaN, NaN, NaN};
        }
    }
    std::array<std::size_t, 4> idx{0, 1, 2, 3};
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return p[a] < p[b]; });
    std::array<double, 4> out{};
    double mx = 0.0;
    for (std::size_t r = 0; r < 4; ++r) {
        const double v = std::min(1.0, static_cast<double>(4 - r) * p[idx[r]]);
        mx = std::max(mx, v);
        out[idx[r]] = mx;
    }
    return out;
}

// =============================================================================
// 1. K1: IDENTIDAD DEL OPERADOR DELTA Y CERO RESERVAS EN PASO (Reglas 2 y 6)
// =============================================================================
void prueba_k1_identidad_y_cero_reservas() {
    std::cout << "[1/8] K1: Identidad bit a bit del ObservadorSaltos real con/sin Lilyth y 0 reservas en paso...\n";
    const std::size_t N_CAL = 3000;
    const std::size_t N_FLUJO = 900;
    std::mt19937_64 rng(2026110601ULL);
    std::normal_distribution<double> nd(0.0, 1.0);

    std::vector<double> x_cal(N_CAL, 0.0);
    for (std::size_t t = 1; t < N_CAL; ++t) {
        x_cal[t] = 0.60 * x_cal[t - 1] + nd(rng);
    }

    ark::metrologia::ConfiguracionObservador cfg_d{};
    cfg_d.integracion = ark::metrologia::Integracion::FLUJO;
    cfg_d.escalas_max = 2;
    cfg_d.medir_tau50 = false;
    cfg_d.arl0_pasos_nativos = 250.0;
    cfg_d.semilla = 7777;

    ObsReal delta_solo;
    ObsReal delta_con_lilyth;
    delta_solo.calibrar(cfg_d, x_cal);
    verificar(delta_solo.escalas() > 0, "K1: delta_solo debe calibrar");
    delta_con_lilyth.calibrar(cfg_d, x_cal);
    verificar(delta_con_lilyth.escalas() > 0, "K1: delta_con_lilyth debe calibrar");

    ark::lilyth::ConfiguracionLilyth cfg_l{};
    cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::CRUDA;
    cfg_l.nulo_sin_memoria = ark::lilyth::LILYTH_NULO_PERMUTACION;
    cfg_l.modo_nulo_lineal = ark::lilyth::ModoNuloLineal::AR_LEVINSON_DURBIN;
    cfg_l.modo_armonico = ark::lilyth::ModoArmonico::PRONY_CURVA_MEDIA;
    cfg_l.b_max_nulo = 4;
    cfg_l.hilos_nulo = 2;

    LilythReal lil;
    verificar(lil.calibrar(delta_con_lilyth, cfg_l, cfg_d, x_cal), "K1: Lilyth debe calibrar");

    std::vector<double> x_flujo(N_FLUJO, 0.0);
    for (std::size_t t = 1; t < N_FLUJO; ++t) {
        x_flujo[t] = 0.60 * x_flujo[t - 1] + nd(rng);
        if (t == 220 || t == 540) x_flujo[t] += 7.5;
    }

    g_reservas_ciclo_caliente.store(0);
    g_contar_reservas.store(true);
    for (std::size_t t = 0; t < N_FLUJO; ++t) {
        const auto s1 = delta_solo.paso(x_flujo[t]);
        const auto s2 = delta_con_lilyth.paso(x_flujo[t]);
        (void)lil.paso(s2, x_flujo[t], true);
        bool eq = (s1.t == s2.t && s1.n_eventos == s2.n_eventos && s1.evento_abierto == s2.evento_abierto && s1.niveles == s2.niveles);
        for (std::size_t j = 0; eq && j < s1.niveles; ++j) {
            eq = (s1.nivel[j].alarma.hay == s2.nivel[j].alarma.hay &&
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
// 2. PRONY SOBRE RECURRENCIA EXACTA A 1e-12 Y PERIODO DE COSENO AMORTIGUADO (SIN RECORTE)
// =============================================================================
void prueba_prony_recurrencia_exacta_y_periodo() {
    std::cout << "[2/8] Prony exacto sobre recurrencia de orden 2 (<= 1e-12) y periodo sin recorte (E1.8)...\n";
    const std::size_t H = 24;

    // Caso A: Coseno amortiguado exacto eta(tau) = r^tau cos(omega * tau), T = 8.0 barras (omega = pi / 4)
    {
        const double r = 0.88;
        const double T_verd = 8.0;
        const double omega = (2.0 * ark::metrologia::est::PI) / T_verd;
        const double a1_verd = 2.0 * r * std::cos(omega);
        const double a2_verd = -(r * r);

        std::vector<double> eta(H, 0.0);
        for (std::size_t tau = 1; tau <= H; ++tau) {
            eta[tau - 1] = std::pow(r, static_cast<double>(tau)) * std::cos(omega * static_cast<double>(tau));
        }

        const auto pr = ark::lilyth::detalle_armonico::ajustar_prony(eta, H, 1.0);
        const auto ref = prony_independiente_ref(eta, 1.0);

        verificar(pr.valido && ref.valido, "Prony exacto: debe ser valido");
        verificar(pr.raices_complejas, "Prony exacto: debe detectar raices complejas (disc < 0)");
        verificar(!pr.raiz_real_negativa, "Prony exacto: no debe marcar raiz real negativa");
        verificar(std::abs(pr.a1 - a1_verd) <= 1e-12, "Prony exacto: |a1 - a1_verd| <= 1e-12 (desde tau = 2 con eta(0) = 1)");
        verificar(std::abs(pr.a2 - a2_verd) <= 1e-12, "Prony exacto: |a2 - a2_verd| <= 1e-12");
        verificar(pr.rmse_recurrencia <= 1e-12, "Prony exacto: rmse_recurrencia <= 1e-12");
        verificar(std::abs(pr.periodo - T_verd) <= 1e-12, "Prony exacto: periodo recupera T = 8.0 a <= 1e-12");
        verificar(std::abs(pr.a1 - ref.a1) <= 1e-15 && std::abs(pr.a2 - ref.a2) <= 1e-15 &&
                  std::abs(pr.periodo - ref.periodo) <= 1e-14,
                  "Prony exacto: coincide con motor de referencia independiente");
    }

    // Caso B: Oscilador subamortiguado de periodo corto T = 3.2 barras (<= 4 barras, E1.8: NO recortar el periodo)
    {
        const double r = 0.85;
        const double T_corto = 3.2;
        const double omega = (2.0 * ark::metrologia::est::PI) / T_corto;
        const double a1_verd = 2.0 * r * std::cos(omega);
        const double a2_verd = -(r * r);

        std::vector<double> eta(H, 0.0);
        for (std::size_t tau = 1; tau <= H; ++tau) {
            eta[tau - 1] = std::pow(r, static_cast<double>(tau)) * std::cos(omega * static_cast<double>(tau));
        }

        const auto pr = ark::lilyth::detalle_armonico::ajustar_prony(eta, H, 1.0);
        verificar(pr.valido && pr.raices_complejas, "Prony T=3.2: debe preservar raices_complejas sin recortar periodos <= 4 barras");
        verificar(std::isfinite(pr.periodo) && std::abs(pr.periodo - T_corto) <= 1e-12,
                  "Prony T=3.2: recupera periodo corto T = 3.2 a <= 1e-12 sin recorte");
        verificar(std::abs(pr.a1 - a1_verd) <= 1e-12 && std::abs(pr.a2 - a2_verd) <= 1e-12,
                  "Prony T=3.2: recupera (a1, a2) a <= 1e-12");
    }
}

// =============================================================================
// 3. P5: AR(2) SUBAMORTIGUADO -> RAICES COMPLEJAS Y RECHAZO CONTRA AMBOS NULOS
// =============================================================================
void prueba_p5_ar2_subamortiguado() {
    std::cout << "[3/8] P5: Respuesta AR(2) subamortiguada en ObservadorSaltos real -> raices complejas y rechazo contra ambos nulos...\n";
    // Mundo sintetico con verdad conocida (D10 / P5):
    // Calma AR(1) suave (phi_calma = 0.45, sigma = 1.0); tras cada choque (t_shock), una accion de control institucional
    // absorbe la cascada con dinamica AR(2) subamortiguada:
    //   r = 0.88, T = 7.0 barras (omega = 2*pi/7) => a1 = 2*r*cos(omega) = 1.0973, a2 = -r^2 = -0.7744 (disc = -1.8935 < 0).
    const std::size_t N_CAL = 3000;
    const std::size_t N_FLUJO = 3200;
    std::mt19937_64 rng(2026110603ULL);
    std::normal_distribution<double> nd(0.0, 1.0);

    std::vector<double> x_cal(N_CAL, 0.0);
    for (std::size_t t = 1; t < N_CAL; ++t) {
        x_cal[t] = 0.82 * x_cal[t - 1] + nd(rng);
    }

    ark::metrologia::ConfiguracionObservador cfg_d{};
    cfg_d.integracion = ark::metrologia::Integracion::FLUJO;
    cfg_d.escalas_max = 2;
    cfg_d.medir_tau50 = false;
    cfg_d.arl0_pasos_nativos = 250.0;
    cfg_d.semilla = 202604;

    ObsReal delta;
    delta.calibrar(cfg_d, x_cal);
    verificar(delta.escalas() > 0, "P5: ObservadorSaltos debe calibrar");

    ark::lilyth::ConfiguracionLilyth cfg_l{};
    cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::CRUDA;
    cfg_l.nulo_sin_memoria = ark::lilyth::LILYTH_NULO_PERMUTACION;
    cfg_l.modo_nulo_lineal = ark::lilyth::ModoNuloLineal::AR_LEVINSON_DURBIN;
    cfg_l.modo_armonico = ark::lilyth::ModoArmonico::PRONY_CURVA_MEDIA;
    cfg_l.b_max_nulo = 8;
    cfg_l.hilos_nulo = 2;
    cfg_l.b_mc_resumen = 499;
    cfg_l.semilla_nulo = 2026110603ULL;

    LilythReal lil;
    verificar(lil.calibrar(delta, cfg_l, cfg_d, x_cal), "P5: Lilyth debe calibrar");

    const double r_osc = 0.88;
    const double T_osc = 7.0;
    const double omega_osc = (2.0 * ark::metrologia::est::PI) / T_osc;
    const double A_shock = 11.0;

    std::vector<double> x_flujo(N_FLUJO, 0.0);
    for (std::size_t t = 1; t < N_FLUJO; ++t) {
        x_flujo[t] = 0.82 * x_flujo[t - 1] + 0.25 * nd(rng);
    }
    for (std::size_t t0 = 220; t0 + 48 < N_FLUJO; t0 += 200) {
        const double signo = ((t0 / 200) % 2 == 0) ? 1.0 : -1.0;
        for (std::size_t tau = 0; tau <= 36; ++tau) {
            const double env = std::pow(r_osc, static_cast<double>(tau)) * std::cos(omega_osc * static_cast<double>(tau));
            x_flujo[t0 + tau] += signo * A_shock * env;
        }
    }

    delta.reiniciar_flujo();
    for (std::size_t t = 0; t < N_FLUJO; ++t) {
        const auto s = delta.paso(x_flujo[t]);
        (void)lil.paso(s, x_flujo[t], true);
    }

    lil.completar_nulos();
    const auto& arm = lil.resumen_armonico_cruda();

    verificar(arm.valido, "P5: ResumenArmonicoYHolm debe ser valido");
    verificar(arm.prony_curva_media.valido, "P5: Prony sobre la curva media debe ser valido");
    verificar(arm.prony_curva_media.discriminante < 0.0 && arm.prony_curva_media.raices_complejas,
              "P5: AR(2) subamortiguado debe producir raices complejas en eta_barra (a1^2 + 4*a2 < 0)");
    verificar(std::isfinite(arm.prony_curva_media.periodo) &&
              std::abs(arm.prony_curva_media.periodo - T_osc) < 1.5,
              "P5: Prony sobre eta_barra debe recuperar el periodo cercano a T = 7.0 barras");
    verificar(arm.prony_curva_media.min_eta_desde_2 < -0.20,
              "P5: min_{tau >= 2} eta_barra(tau) debe ser claramente negativo (sobrepaso < -0.20)");
    verificar(arm.sobrepaso.p_sobrepaso_curva_iid < 0.01,
              "P5: contraste de sobrepaso contra el nulo iid debe rechazar (p < 0.01)");
    verificar(arm.sobrepaso.p_sobrepaso_curva_onsager < 0.01,
              "P5: contraste de sobrepaso contra el nulo lineal de calma debe rechazar (p < 0.01)");
    verificar(arm.familia_holm_valida && arm.rechaza_holm[2] && arm.rechaza_holm[3],
              "P5: ambos contrastes de sobrepaso (iid y Onsager) deben mantenerse rechazados tras ark::holm");
    verificar(arm.sobrepaso.sobrepaso_por_evento_iid.valido &&
              arm.sobrepaso.sobrepaso_por_evento_onsager.valido &&
              arm.sobrepaso.sobrepaso_por_evento_onsager.p_bilateral < 0.05,
              "P5: SOBREPASO_POR_EVENTO tambien detecta el sobrepaso contra el nulo");
    verificar(arm.sobrepaso.tasa_etiqueta_reales > 0.50,
              "P5: la mayoria de los eventos reales con respuesta subamortiguada reciben ARMONICO_SOSPECHOSO");

    // Code (revision D1c): el contraste de sobrepaso usa los nulos EMPAREJADOS de cada evento real (los de la Etapa 3) y el min de
    // eta_barra real sobre esos mismos eventos. Se recalcula aqui desde los datos publicos y tiene que dar lo mismo bit a bit.
    {
        const std::size_t H = std::min<std::size_t>(lil.configuracion().t_max_efectivo, 64);
        const auto& emp = lil.nulos_emparejados(false, false);
        std::vector<const std::vector<ark::lilyth::EventoNuloLilyth<64>>*> sets;
        std::vector<double> suma(H, 0.0);
        std::vector<std::size_t> cuenta(H, 0);
        for (std::size_t idx = 0; idx < emp.size(); ++idx) {
            if (emp[idx].empty()) continue;
            sets.push_back(&emp[idx]);
            const auto& rama = lil.cerradas_historico()[idx].cruda;
            const std::size_t pk = std::min<std::size_t>(rama.pasos_observados, H);
            for (std::size_t tau = 2; tau <= pk; ++tau)
                if (std::isfinite(rama.eta[tau - 1])) { suma[tau - 1] += rama.eta[tau - 1]; ++cuenta[tau - 1]; }
        }
        double mn = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t tau = 2; tau <= H; ++tau) {
            if (cuenta[tau - 1] == 0) continue;
            const double m = suma[tau - 1] / static_cast<double>(cuenta[tau - 1]);
            if (!std::isfinite(mn) || m < mn) mn = m;
        }
        double med = 0.0, pp = 0.0;
        ark::lilyth::detalle_armonico::contraste_min_eta_barra_emparejado<64>(
            mn, sets, H, cfg_l.b_mc_resumen, cfg_l.semilla_nulo ^ 0x5A5A5A5AA5A5A5A5ULL, med, pp);
        std::cout << "  [P5 emparejado] eventos con nulo lineal = " << sets.size() << ", min eta_barra = " << mn
                  << ", p = " << pp << " (resumen: " << arm.sobrepaso.p_sobrepaso_curva_onsager << ")\n";
        verificar(!sets.empty() && mn == arm.sobrepaso.min_eta_barra_real_emparejado_onsager && pp == arm.sobrepaso.p_sobrepaso_curva_onsager &&
                  med == arm.sobrepaso.mediana_min_eta_barra_onsager,
                  "P5 emparejado: el contraste de sobrepaso Onsager = recalculo con los nulos emparejados de la Etapa 3 (bit a bit)");
    }
}

// =============================================================================
// 4. P5b: AR(1) CON phi < 0 (REBOTE BID-ASK) -> RAIZ REAL NEGATIVA, NUNCA OSCILADOR
// =============================================================================
void prueba_p5b_ar1_rebote_negativo() {
    std::cout << "[4/8] P5b: AR(1) con phi < 0 (rebote bid-ask) -> raiz real negativa, nunca oscilador...\n";

    // 4a. Verificacion analitica directa sobre una curva con rebote real negativo (lambda_1 = -0.65, lambda_2 = 0.10 => disc > 0)
    {
        const std::size_t H = 20;
        const double l1 = -0.65;
        const double l2 = 0.10;
        std::vector<double> eta_rebote(H, 0.0);
        for (std::size_t tau = 1; tau <= H; ++tau) {
            eta_rebote[tau - 1] = 0.85 * std::pow(l1, static_cast<double>(tau)) + 0.15 * std::pow(l2, static_cast<double>(tau));
        }
        const auto pr = ark::lilyth::detalle_armonico::ajustar_prony(eta_rebote, H, 1.0);
        verificar(pr.valido, "P5b analitico: Prony debe ser valido");
        verificar(pr.discriminante >= 0.0, "P5b analitico: discriminante debe ser >= 0");
        verificar(!pr.raices_complejas, "P5b analitico: rebote phi < 0 NUNCA debe clasificarse como oscilador de raices complejas");
        verificar(pr.raiz_real_negativa, "P5b analitico: debe identificar raiz_real_negativa == true");
    }

    // 4b. Verificacion integrada en el pipeline con ObservadorSaltos sobre un AR(1) con phi = -0.65
    {
        const std::size_t N_CAL = 3000;
        const std::size_t N_FLUJO = 2600;
        const double phi_reb = -0.65;
        std::mt19937_64 rng(2026110604ULL);
        std::normal_distribution<double> nd(0.0, 1.0);

        std::vector<double> x_cal(N_CAL, 0.0);
        for (std::size_t t = 1; t < N_CAL; ++t) {
            x_cal[t] = phi_reb * x_cal[t - 1] + nd(rng);
        }

        ark::metrologia::ConfiguracionObservador cfg_d{};
        cfg_d.integracion = ark::metrologia::Integracion::FLUJO;
        cfg_d.escalas_max = 2;
    cfg_d.medir_tau50 = false;
        cfg_d.arl0_pasos_nativos = 250.0;
        cfg_d.semilla = 202605;

        ObsReal delta;
        delta.calibrar(cfg_d, x_cal);
    verificar(delta.escalas() > 0, "P5b pipeline: ObservadorSaltos debe calibrar");

        ark::lilyth::ConfiguracionLilyth cfg_l{};
        cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::CRUDA;
        cfg_l.nulo_sin_memoria = ark::lilyth::LILYTH_NULO_PERMUTACION;
        cfg_l.modo_nulo_lineal = ark::lilyth::ModoNuloLineal::AR_LEVINSON_DURBIN;
        cfg_l.b_max_nulo = 6;
        cfg_l.hilos_nulo = 2;
        cfg_l.t_max_efectivo = 16;

        LilythReal lil;
        verificar(lil.calibrar(delta, cfg_l, cfg_d, x_cal), "P5b pipeline: Lilyth debe calibrar");

        std::vector<double> x_flujo(N_FLUJO, 0.0);
        for (std::size_t t = 1; t < N_FLUJO; ++t) {
            x_flujo[t] = phi_reb * x_flujo[t - 1] + 0.25 * nd(rng);
        }
        for (std::size_t t0 = 200; t0 + 30 < N_FLUJO; t0 += 180) {
            const double sgn = ((t0 / 180) % 2 == 0) ? 1.0 : -1.0;
            for (std::size_t tau = 0; tau <= 20; ++tau) {
                x_flujo[t0 + tau] += sgn * 10.0 * std::pow(phi_reb, static_cast<double>(tau));
            }
        }

        delta.reiniciar_flujo();
        for (std::size_t t = 0; t < N_FLUJO; ++t) {
            const auto s = delta.paso(x_flujo[t]);
            (void)lil.paso(s, x_flujo[t], true);
        }
        lil.completar_nulos();
        const auto& arm = lil.resumen_armonico_cruda();
        verificar(arm.prony_curva_media.valido, "P5b pipeline: Prony de curva media valido");
        verificar(!arm.prony_curva_media.raices_complejas && arm.prony_curva_media.discriminante >= 0.0,
                  "P5b pipeline: rebote AR(1) phi < 0 tiene discriminante >= 0 y NUNCA da raices_complejas");
        verificar(arm.prony_curva_media.raiz_real_negativa,
                  "P5b pipeline: marca raiz_real_negativa == true");
    }
}

// =============================================================================
// 5. TASA DE ARMONICO_SOSPECHOSO EN UN AR(1) <= ALFA (SIN VETAR EVENTOS)
// =============================================================================
void prueba_tasa_etiqueta_ar1() {
    std::cout << "[5/8] Tasa de etiqueta ARMONICO_SOSPECHOSO en un AR(1) <= alfa (con su tasa en el nulo al lado, nunca veta)...\n";
    const std::size_t N_CAL = 3000;
    const std::size_t N_FLUJO = 3600;
    const double phi = 0.72;
    std::mt19937_64 rng(2026110605ULL);
    std::normal_distribution<double> nd(0.0, 1.0);

    std::vector<double> x_cal(N_CAL, 0.0);
    for (std::size_t t = 1; t < N_CAL; ++t) {
        x_cal[t] = phi * x_cal[t - 1] + nd(rng);
    }

    ark::metrologia::ConfiguracionObservador cfg_d{};
    cfg_d.integracion = ark::metrologia::Integracion::FLUJO;
    cfg_d.escalas_max = 2;
    cfg_d.medir_tau50 = false;
    cfg_d.arl0_pasos_nativos = 220.0;
    cfg_d.semilla = 202606;

    ObsReal delta;
    delta.calibrar(cfg_d, x_cal);
    verificar(delta.escalas() > 0, "Etiqueta AR(1): ObservadorSaltos debe calibrar");

    ark::lilyth::ConfiguracionLilyth cfg_l{};
    cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::CRUDA;
    cfg_l.nulo_sin_memoria = ark::lilyth::LILYTH_NULO_PERMUTACION;
    cfg_l.modo_nulo_lineal = ark::lilyth::ModoNuloLineal::AR_LEVINSON_DURBIN;
    cfg_l.alfa = 0.05;
    cfg_l.b_max_nulo = 8;
    cfg_l.hilos_nulo = 2;
    cfg_l.semilla_nulo = 2026110605ULL;

    LilythReal lil;
    verificar(lil.calibrar(delta, cfg_l, cfg_d, x_cal), "Etiqueta AR(1): Lilyth debe calibrar");

    std::vector<double> x_flujo(N_FLUJO, 0.0);
    for (std::size_t t = 1; t < N_FLUJO; ++t) {
        double eps = nd(rng);
        if (t >= 200 && t + 40 < N_FLUJO && (t % 160 == 0)) {
            eps += ((t / 160) % 2 == 0 ? 8.5 : -8.5);
        }
        x_flujo[t] = phi * x_flujo[t - 1] + eps;
    }

    delta.reiniciar_flujo();
    for (std::size_t t = 0; t < N_FLUJO; ++t) {
        const auto s = delta.paso(x_flujo[t]);
        (void)lil.paso(s, x_flujo[t], true);
    }

    lil.completar_nulos();
    const auto& arm = lil.resumen_armonico_cruda();

    verificar(std::isfinite(arm.sobrepaso.tasa_etiqueta_nulo_iid) &&
              arm.sobrepaso.tasa_etiqueta_nulo_iid <= cfg_l.alfa + 1e-12,
              "Etiqueta AR(1): tasa de ARMONICO_SOSPECHOSO en el nulo iid <= alfa");
    verificar(std::isfinite(arm.sobrepaso.tasa_etiqueta_nulo_onsager) &&
              arm.sobrepaso.tasa_etiqueta_nulo_onsager <= cfg_l.alfa + 1e-12,
              "Etiqueta AR(1): tasa de ARMONICO_SOSPECHOSO en el nulo lineal <= alfa");
    verificar(std::isfinite(arm.sobrepaso.tasa_etiqueta_reales) &&
              arm.sobrepaso.tasa_etiqueta_reales <= cfg_l.alfa + 1e-12,
              "Etiqueta AR(1): tasa de ARMONICO_SOSPECHOSO en eventos reales AR(1) <= alfa");
    // Nunca veta: todos los eventos cerrados permanecen registrados en Aalen-Johansen
    verificar(lil.aalen_johansen_cruda().total.n_eventos_registrados >= 10,
              "Etiqueta AR(1): ARMONICO_SOSPECHOSO nunca veta eventos");
}

// =============================================================================
// 6. PRONY_FASE_REESCALADA (E1.20): CONTEO Y EXCLUSION DE NO RESTITUIDOS
// =============================================================================
void prueba_fase_reescalada_e1_20() {
    std::cout << "[6/8] PRONY_FASE_REESCALADA (E1.20): exclusion y conteo de eventos no restituidos (tau_rest = inf)...\n";
    constexpr std::size_t T_MAX = 32;
    const std::size_t H = 24;

    // 6 eventos restituidos (con tau_rest en {4, 4, 6, 6, 8, 8} => mediana = 6, paso_rejilla = 1/6)
    // y 3 eventos NO restituidos (tau_rest = NaN, eta constante = 0.95).
    // Si un mutante incluye los no restituidos en la curva reescalada, contamina la curva y cambia eventos_incluidos.
    std::vector<std::array<double, T_MAX>> etas(9);
    std::vector<std::size_t> pasos(9, H);
    std::vector<double> taus_rest{4.0, 4.0, 6.0, 6.0, 8.0, 8.0, NaN, NaN, NaN};

    const double r_base = 0.86;
    const double omega_base = (2.0 * ark::metrologia::est::PI) / 8.0;
    for (std::size_t k = 0; k < 6; ++k) {
        etas[k].fill(NaN);
        const double escala_t = taus_rest[k] / 6.0;
        for (std::size_t tau = 1; tau <= H; ++tau) {
            const double t_eq = static_cast<double>(tau) / escala_t;
            etas[k][tau - 1] = std::pow(r_base, t_eq) * std::cos(omega_base * t_eq);
        }
    }
    for (std::size_t k = 6; k < 9; ++k) {
        etas[k].fill(0.95);
    }

    const auto fr = ark::lilyth::detalle_armonico::construir_fase_reescalada<T_MAX>(etas, pasos, taus_rest, H);
    verificar(fr.valido, "Fase reescalada: debe ser valida");
    verificar(fr.excluidos_sin_restituir == 3, "Fase reescalada (E1.20): debe contar exactamente 3 excluidos_sin_restituir");
    verificar(fr.eventos_incluidos == 6, "Fase reescalada (E1.20): debe incluir exactamente los 6 eventos restituidos");
    verificar(std::abs(fr.mediana_tau_rest - 6.0) <= 1e-12, "Fase reescalada: mediana(tau_rest) == 6.0");
    verificar(std::abs(fr.paso_rejilla - (1.0 / 6.0)) <= 1e-12, "Fase reescalada: paso_rejilla == 1 / mediana(tau_rest)");
    // En m = 6 (tau_tilde = 1.0), cada evento restituido k se evalua exactamente en t_fis = tau_rest_k,
    // donde t_eq = 6.0, por lo que eta_reescalada[5] == r_base^6 * cos(omega_base * 6)
    const double esperado_en_tau_tilde_1 = std::pow(r_base, 6.0) * std::cos(omega_base * 6.0);
    verificar(std::abs(fr.eta_reescalada[5] - esperado_en_tau_tilde_1) <= 1e-12,
              "Fase reescalada: en tau_tilde = 1.0 recupera exactamente r^6 cos(6*omega) sin contaminacion de no restituidos");
    verificar(fr.prony.valido && fr.prony.raices_complejas,
              "Fase reescalada: Prony sobre la curva reescalada detecta las raices complejas");
}

// =============================================================================
// 7. HOLM SOBRE LOS 6 P (E1.16 + E1.24, H.1) Y FALLA DURA ANTE NaN; L+H NO DECIDIBLE ENTRA CON p = 1 (H.4, control R9)
// =============================================================================
template <std::size_t M>
std::array<double, M> holm_m_independiente_ref(const std::array<double, M>& p) {
    std::array<double, M> nan_out{};
    nan_out.fill(NaN);
    for (double x : p) {
        if (!std::isfinite(x) || x < 0.0 || x > 1.0) return nan_out;
    }
    std::array<std::size_t, M> idx{};
    for (std::size_t k = 0; k < M; ++k) idx[k] = k;
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return p[a] < p[b]; });
    std::array<double, M> out{};
    double mx = 0.0;
    for (std::size_t r = 0; r < M; ++r) {
        const double v = std::min(1.0, static_cast<double>(M - r) * p[idx[r]]);
        mx = std::max(mx, v);
        out[idx[r]] = mx;
    }
    return out;
}

void prueba_holm_cuatro_y_falla_dura() {
    std::cout << "[7/8] Familia de 6 hipotesis en ark::holm (E1.16 + E1.24) y falla dura ante cualquier p = NaN...\n";
    constexpr std::size_t T_MAX = 32;

    // 7a. 6 valores p validos con L+H decidible -> coincide con la referencia de M = 6
    {
        ark::lilyth::ResumenArmonicoYHolm<T_MAX> res{};
        const std::array<double, 6> p_in{0.04, 0.20, 0.006, 0.015, 0.30, 0.007};
        ark::lilyth::detalle_armonico::aplicar_familia_holm_seis<T_MAX>(
            p_in[0], p_in[1], p_in[2], p_in[3], true, p_in[4], p_in[5], 0.05, res);
        const auto ref = holm_m_independiente_ref<6>(p_in);
        verificar(res.familia_holm_valida && res.n_familia == 6 && !res.lh_no_decidible, "Holm 6p: familia valida de 6 con L+H decidible");
        for (std::size_t i = 0; i < 6; ++i) {
            verificar(std::isfinite(res.p_holm[i]) && std::abs(res.p_holm[i] - ref[i]) <= 1e-15,
                      "Holm 6p: p_holm[i] debe coincidir con la referencia de 6 hipotesis");
        }
        // M = 6: 6 * 0.006 = 0.036 (rechaza), max(0.036, 5 * 0.007 = 0.035) = 0.036 (rechaza), 4 * 0.015 = 0.06 (no rechaza).
        verificar(res.rechaza_holm[2] && res.rechaza_holm[5] && !res.rechaza_holm[3] && !res.rechaza_holm[0] && !res.rechaza_holm[1] &&
                      !res.rechaza_holm[4],
                  "Holm 6p: decisiones de rechazo exactas para M = 6");
    }

    // 7a'. R9: L+H NO DECIDIBLE -> sus dos p entran como 1, la familia sigue siendo de 6 y las otras cuatro decisiones son las de
    //      una referencia escrita aparte sobre [p0, p1, p2, p3, 1, 1] (aunque se le pasen p de L+H que rechazarian).
    {
        ark::lilyth::ResumenArmonicoYHolm<T_MAX> res{};
        const std::array<double, 4> p4{0.04, 0.20, 0.006, 0.007};
        ark::lilyth::detalle_armonico::aplicar_familia_holm_seis<T_MAX>(p4[0], p4[1], p4[2], p4[3], false, 0.0001, 0.0001, 0.05, res);
        const auto ref = holm_m_independiente_ref<6>({p4[0], p4[1], p4[2], p4[3], 1.0, 1.0});
        bool ok = res.familia_holm_valida && res.n_familia == 6 && res.lh_no_decidible && res.p_crudos[4] == 1.0 && res.p_crudos[5] == 1.0;
        for (std::size_t i = 0; i < 6; ++i) ok = ok && std::abs(res.p_holm[i] - ref[i]) <= 1e-15 && res.rechaza_holm[i] == (ref[i] < 0.05);
        // con M = 6: 6 * 0.006 = 0.036 y max(0.036, 5 * 0.007 = 0.035) rechazan; con M = 4 tambien, pero 3 * 0.04 = 0.12 no: lo que
        // distingue M = 6 de M = 4 aqui es el p_holm, no la decision; se exige igualdad de p_holm.
        verificar(ok, "R9: L+H NO DECIDIBLE entra con p = 1, familia de 6, p_holm = referencia sobre [p0..p3, 1, 1]");
    }

    // 7b. Un p = NaN en cualquiera de las 6 posiciones (L+H decidible) -> falla dura de toda la familia
    for (std::size_t pos_nan = 0; pos_nan < 6; ++pos_nan) {
        std::array<double, 6> p_in{0.005, 0.010, 0.002, 0.008, 0.003, 0.004};
        p_in[pos_nan] = NaN;
        ark::lilyth::ResumenArmonicoYHolm<T_MAX> res{};
        ark::lilyth::detalle_armonico::aplicar_familia_holm_seis<T_MAX>(
            p_in[0], p_in[1], p_in[2], p_in[3], true, p_in[4], p_in[5], 0.05, res);

        verificar(!res.familia_holm_valida, "Holm NaN: familia_holm_valida debe ser false si un p es NaN");
        for (std::size_t i = 0; i < 6; ++i) {
            verificar(std::isnan(res.p_holm[i]), "Holm NaN: todos los p_holm[i] deben ser NaN (falla dura)");
            verificar(!res.rechaza_holm[i], "Holm NaN: ningun rechaza_holm[i] puede ser true ante falla dura");
        }
    }
}

// =============================================================================
// 8. CONTROL DE DISCRIMINACION DE MUTANTES (TESTIGO DETERMINISTA)
// =============================================================================
void prueba_testigo_mutantes_etapa4() {
    std::cout << "[8/8] Control testigo determinista de los 6 mutantes de la Etapa 4...\n";
    // Verifica en menos de 10 ms las 6 invariantes exactas que cazan M01..M06:
    // M01 (signo discriminante), M02 (recortar periodo <= 4), M03 (Prony desde tau = 1),
    // M04 (incluir no restituidos en fase reescalada), M05 (Holm sobre 3 p), M06 (omitir/sanear NaN en Holm).
    const std::size_t H = 16;
    const double r = 0.85;
    const double T = 3.5; // <= 4.0 barras (caza M02)
    const double w = (2.0 * ark::metrologia::est::PI) / T;
    std::vector<double> eta(H, 0.0);
    for (std::size_t tau = 1; tau <= H; ++tau) {
        eta[tau - 1] = std::pow(r, static_cast<double>(tau)) * std::cos(w * static_cast<double>(tau));
    }
    const auto pr = ark::lilyth::detalle_armonico::ajustar_prony(eta, H, 1.0);
    verificar(pr.valido && pr.raices_complejas && std::isfinite(pr.periodo) && std::abs(pr.periodo - T) <= 1e-12 &&
              pr.rmse_recurrencia <= 1e-12,
              "Testigo M01/M02/M03: Prony exacto con T = 3.5 barras");

    // Testigo M07: trayectoria con raiz real negativa (disc >= 0) y min_eta < cuantil_alfa NO debe etiquetarse ARMONICO_SOSPECHOSO
    std::array<double, 16> eta_real_neg{};
    for (std::size_t tau = 1; tau <= 16; ++tau) {
        eta_real_neg[tau - 1] = std::pow(-0.75, static_cast<double>(tau));
    }
    const double mn_rn = ark::lilyth::detalle_armonico::min_eta_en_rango<16>(eta_real_neg, 16, 2, 16);
    verificar(!ark::lilyth::detalle_armonico::evaluar_etiqueta_armonico_evento<16>(eta_real_neg, 16, mn_rn, 16, 0.0),
              "Testigo M07: evaluar_etiqueta_armonico_evento exige raices_complejas == true");

    // Testigo M08: min_eta_desde_2 evalua [2, H], nunca tau = 1
    std::vector<double> eta_tau1_bajo(8, 0.20);
    eta_tau1_bajo[0] = -0.95; // tau = 1
    eta_tau1_bajo[1] = -0.15; // tau = 2
    eta_tau1_bajo[2] = 0.05;
    const auto pr_m8 = ark::lilyth::detalle_armonico::ajustar_prony(eta_tau1_bajo, 8, 1.0);
    verificar(std::abs(pr_m8.min_eta_desde_2 - (-0.15)) <= 1e-12,
              "Testigo M08: min_eta_desde_2 excluye tau = 1 y toma el minimo en [2, T_max]");
}

}  // namespace

int main(int argc, char** argv) {
    const bool modo_rapido = (argc > 1 && std::strcmp(argv[1], "--rapido") == 0);
    std::cout << "====================================================================\n";
    std::cout << " BATERIA ETAPA 4 -- OPERADOR DE LILYTH (ARMONICO Y HOLM, C++20)\n";
    std::cout << "====================================================================\n";

    prueba_prony_recurrencia_exacta_y_periodo();
    prueba_fase_reescalada_e1_20();
    prueba_holm_cuatro_y_falla_dura();
    prueba_testigo_mutantes_etapa4();

    if (!modo_rapido) {
        prueba_k1_identidad_y_cero_reservas();
        prueba_p5_ar2_subamortiguado();
        prueba_p5b_ar1_rebote_negativo();
        prueba_tasa_etiqueta_ar1();
    } else {
        // En modo --rapido (para corrida agil de mutantes), ejecutamos tambien P5b analitico
        prueba_p5b_ar1_rebote_negativo();
    }

    std::cout << "====================================================================\n";
    if (g_fallas == 0) {
        std::cout << " RESULTADO: TODOS LOS CONTROLES DE LA ETAPA 4 PASARON (0 fallas)\n";
        std::cout << "====================================================================\n";
        return 0;
    }
    std::cerr << " RESULTADO: " << g_fallas << " FALLA(S) DETECTADA(S)\n";
    std::cerr << "====================================================================\n";
    return 1;
}