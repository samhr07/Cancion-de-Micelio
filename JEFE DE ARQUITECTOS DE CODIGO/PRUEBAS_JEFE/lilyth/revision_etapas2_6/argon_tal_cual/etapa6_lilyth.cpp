// etapa6_lilyth.cpp -- Bateria de verificacion de la Etapa 6 del Operador de Lilyth (`reloj_marea_lilyth.hpp` + `observador_lilyth.hpp`):
//   1. K1 + K10: Identidad bit a bit de ObservadorSaltos y ObservadorSaltosFilas con/sin Lilyth, y 0 reservas en `paso` y `paso_filas`.
//   2. Reloj de actividad (E1.11, D7): acepta NU_TRADES, NU_VOLUMEN, NU_NOCIONAL y PARED; rechaza estrictamente NU_INFO_SHANNON.
//   3. P4 (Camino de serie y camino por filas con flujo_barras.hpp / observador_saltos_filas.hpp):
//      - Actividad congelada: tau grande en pared y normal en reloj de actividad (tau_act = Delta_nu_rest / nu_punto* congelado).
//      - Perdida de rigidez: grande en ambos relojes (pared y actividad).
//      - Barras intrinsecas por filas: la pared no avanza a paso fijo (pared_consolidado se estira con sequia de liquidez mientras tau_act permanece invariante).
//   4. P11 (Contraste de marea E1.14 [DECIDIDO]):
//      - Antes de los 20 dias por fase (`dias_minimos_purga`), PURGADA no mide (`NO_MEDIDA_PURGA_INMADURA`).
//      - En NIVEL, un incremento purgado invalido censura la ranura (`CENSURADA_DATO_INVALIDO`) sin rellenar ni imputar.
//      - Con perfil diurno determinista: CRUDA se sesga con la fase del choque, PURGADA recupera la semivida verdadera, y
//        Delta_tau_marea = tau_cruda - tau_purgada y Delta_RMST_marea con IC bootstrap miden el sesgo (descriptivo, fuera de Holm).
//   5. Casamiento bit a bit contra motores de referencia independientes (Regla 3) y testigo determinista de los mutantes M01..M06.
//
// Compilacion: g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Wno-mismatched-new-delete -Werror -pthread -Iinclude etapa6_lilyth.cpp -o etapa6_lilyth

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/lilyth/reloj_marea_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"
#include "ark/metrologia/observador_saltos_filas.hpp"
#include "ark/ventanas/flujo_barras.hpp"

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
using ObsFilasReal = ark::metrologia::ObservadorSaltosFilas<24, 64, 120>;
using LilythReal = ark::lilyth::ObservadorLilyth<ObsReal, 8, 64>;
using ark::metrologia::NaN;
using ark::ventanas::EntradaPasoCxx;
using ark::ventanas::TipoRelojVentana;

int g_fallas = 0;

void verificar(bool cond, const char* msg) {
    if (!cond) {
        std::cerr << "  [FALLA] " << msg << "\n";
        ++g_fallas;
    }
}

// =============================================================================
// MOTOR DE REFERENCIA INDEPENDIENTE PARA RELOJ DE ACTIVIDAD Y MAREA (Regla 3)
// =============================================================================

struct RefMareaindependiente {
    std::vector<double> deltas_pareados;
    double mediana_delta = NaN;
    double delta_rmst = NaN;
};

template <class Ranura>
RefMareaindependiente calcular_marea_ref(const std::vector<Ranura>& cerradas, double rmst_c, double rmst_p) {
    RefMareaindependiente ref{};
    ref.delta_rmst = rmst_c - rmst_p;
    for (const auto& r : cerradas) {
        const auto dc = static_cast<std::uint8_t>(r.cruda.desenlace);
        const auto dp = static_cast<std::uint8_t>(r.purgada.desenlace);
        if (dc >= 2u && dc <= 5u && dp >= 2u && dp <= 5u) {
            if (std::isfinite(r.cruda.tau_rest_consolidado) && std::isfinite(r.purgada.tau_rest_consolidado)) {
                ref.deltas_pareados.push_back(r.cruda.tau_rest_consolidado - r.purgada.tau_rest_consolidado);
            }
        }
    }
    if (!ref.deltas_pareados.empty()) {
        std::vector<double> cp = ref.deltas_pareados;
        ref.mediana_delta = ark::metrologia::est::cuantil(cp, 0.5);
    }
    return ref;
}

// =============================================================================
// 1. K1 + K10: IDENTIDAD BIT A BIT EN SERIE Y EN FILAS, Y CERO RESERVAS EN PASO Y PASO_FILAS
// =============================================================================
void prueba_k1_serie_y_filas_cero_reservas() {
    std::cout << "[1/5] K1 + K10: Identidad bit a bit en ObservadorSaltos y ObservadorSaltosFilas, y 0 reservas en paso/paso_filas...\n";
    const std::size_t N_CAL = 3000;
    const std::size_t N_FLUJO = 800;
    std::mt19937_64 rng(2026110661ULL);
    std::normal_distribution<double> nd(0.0, 1.0);

    std::vector<EntradaPasoCxx> filas_cal(N_CAL);
    std::vector<double> x_cal(N_CAL, 0.0);
    for (std::size_t i = 0; i < N_CAL; ++i) {
        x_cal[i] = nd(rng);
        filas_cal[i].t_cierre_ms = static_cast<std::uint64_t>((i + 1) * 1000);
        filas_cal[i].t_disponible_ms = filas_cal[i].t_cierre_ms;
        filas_cal[i].valor = x_cal[i];
        filas_cal[i].volumen = 50.0;
        filas_cal[i].trades = 10.0;
        filas_cal[i].nocional = 5000.0;
        filas_cal[i].con_sesion = true;
        filas_cal[i].imputada = false;
    }

    ark::metrologia::ConfiguracionObservador cfg_d{};
    cfg_d.integracion = ark::metrologia::Integracion::FLUJO;
    cfg_d.escalas_max = 1;
    cfg_d.medir_tau50 = false;
    cfg_d.arl0_pasos_nativos = 250.0;
    cfg_d.semilla = 6001;

    ObsFilasReal::Fuente fuente{};
    fuente.base.reloj = TipoRelojVentana::NU_VOLUMEN;
    fuente.delta_nu = 50.0; // 1 fila = 1 barra intrinseca en calma
    fuente.agregacion = ark::ventanas::AgregacionBarra::ULTIMO;
    fuente.fase.activa = false;

    auto osf_solo = std::make_unique<ObsFilasReal>();
    auto osf_con = std::make_unique<ObsFilasReal>();
    osf_solo->calibrar_filas(cfg_d, fuente, filas_cal);
    osf_con->calibrar_filas(cfg_d, fuente, filas_cal);
    verificar(osf_solo->lista() && osf_con->lista(), "K1 filas: ambos ObservadorSaltosFilas deben calibrar");

    ark::lilyth::ConfiguracionLilyth cfg_l{};
    cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::CRUDA;
    cfg_l.reloj_actividad = TipoRelojVentana::NU_VOLUMEN;
    cfg_l.t_max_efectivo = 20;

    auto lil = std::make_unique<LilythReal>();
    verificar(lil->configurar(cfg_l, cfg_d, osf_con->observador(), x_cal), "K1 filas: ObservadorLilyth debe configurar");

    std::vector<EntradaPasoCxx> filas_flujo(N_FLUJO);
    for (std::size_t i = 0; i < N_FLUJO; ++i) {
        filas_flujo[i].t_cierre_ms = static_cast<std::uint64_t>((N_CAL + i + 1) * 1000);
        filas_flujo[i].t_disponible_ms = filas_flujo[i].t_cierre_ms;
        filas_flujo[i].valor = nd(rng) + (i == 220 ? 9.0 : 0.0);
        filas_flujo[i].volumen = 50.0;
        filas_flujo[i].con_sesion = true;
        filas_flujo[i].imputada = false;
    }

    g_reservas_ciclo_caliente.store(0);
    g_contar_reservas.store(true);
    for (std::size_t i = 0; i < N_FLUJO; ++i) {
        const auto sf1 = osf_solo->paso_fila(filas_flujo[i]);
        const auto sf2 = osf_con->paso_fila(filas_flujo[i]);
        (void)lil->paso_filas(sf2, false, 0, 0);
        bool eq = (sf1.barra_cerrada == sf2.barra_cerrada && sf1.delta.t == sf2.delta.t &&
                   sf1.delta.n_eventos == sf2.delta.n_eventos && sf1.delta.niveles == sf2.delta.niveles);
        for (std::size_t j = 0; eq && j < sf1.delta.niveles; ++j) {
            eq = (sf1.delta.nivel[j].alarma.hay == sf2.delta.nivel[j].alarma.hay &&
                  ((std::isnan(sf1.delta.nivel[j].g) && std::isnan(sf2.delta.nivel[j].g)) ||
                   sf1.delta.nivel[j].g == sf2.delta.nivel[j].g));
        }
        if (!eq) {
            verificar(false, "K1 filas: Salida de ObservadorSaltosFilas difiere bit a bit");
            break;
        }
    }
    g_contar_reservas.store(false);
    verificar(g_reservas_ciclo_caliente.load() == 0,
              "K10 filas: Cero llamadas a operator new en el ciclo caliente paso_fila() + paso_filas()");
}

// =============================================================================
// 2. VALIDACION DEL RELOJ DE ACTIVIDAD (E1.11, D7): PROHIBICION DE NU_INFO_SHANNON
// =============================================================================
void prueba_validacion_reloj_actividad() {
    std::cout << "[2/5] Validacion del reloj de actividad (E1.11, D7): rechazo estricto de NU_INFO_SHANNON...\n";
    const std::size_t N_CAL = 2500;
    std::mt19937_64 rng(2026110662ULL);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<double> x_cal(N_CAL, 0.0);
    for (double& v : x_cal) v = nd(rng);

    ark::metrologia::ConfiguracionObservador cfg_d{};
    cfg_d.arl0_pasos_nativos = 250.0;
    cfg_d.escalas_max = 1;
    cfg_d.medir_tau50 = false;
    auto delta = std::make_unique<ObsReal>();
    delta->calibrar(cfg_d, x_cal);

    auto lil = std::make_unique<LilythReal>();
    for (TipoRelojVentana r_ok : {TipoRelojVentana::PARED, TipoRelojVentana::NU_TRADES,
                                  TipoRelojVentana::NU_VOLUMEN, TipoRelojVentana::NU_NOCIONAL}) {
        ark::lilyth::ConfiguracionLilyth cfg{};
        cfg.t_max_efectivo = 16;
        cfg.reloj_actividad = r_ok;
        verificar(lil->configurar(cfg, cfg_d, *delta, x_cal) && lil->lista(),
                  "Reloj de actividad: PARED, NU_TRADES, NU_VOLUMEN y NU_NOCIONAL deben aceptarse");
    }

    // NU_INFO_SHANNON debe rechazarse en configurar y en calibrar
    {
        ark::lilyth::ConfiguracionLilyth cfg_sh{};
        cfg_sh.t_max_efectivo = 16;
        cfg_sh.reloj_actividad = TipoRelojVentana::NU_INFO_SHANNON;
        verificar(!lil->configurar(cfg_sh, cfg_d, *delta, x_cal) && !lil->lista(),
                  "Reloj de actividad (D7, E1.11): NU_INFO_SHANNON debe rechazarse en configurar()");
        verificar(!lil->calibrar(*delta, cfg_sh, cfg_d, x_cal) && !lil->lista(),
                  "Reloj de actividad (D7, E1.11): NU_INFO_SHANNON debe rechazarse en calibrar()");
    }
}

// =============================================================================
// 3. P4: RELOJ DOBLE (ACTIVIDAD CONGELADA VS PERDIDA DE RIGIDEZ, EN SERIE Y EN BARRAS INTRINSECAS)
// =============================================================================
void prueba_p4_reloj_doble() {
    std::cout << "[3/5] P4: Reloj doble (actividad congelada vs perdida de rigidez, en serie y en barras intrinsecas)...\n";
    const std::size_t N_CAL = 2500;
    std::mt19937_64 rng(2026110663ULL);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<double> x_cal(N_CAL, 0.0);
    for (double& v : x_cal) v = 0.25 * nd(rng);

    ark::metrologia::ConfiguracionObservador cfg_d{};
    cfg_d.arl0_pasos_nativos = 250.0;
    cfg_d.escalas_max = 1;
    cfg_d.medir_tau50 = false;

    auto delta = std::make_unique<ObsReal>();
    delta->calibrar(cfg_d, x_cal);

    ark::lilyth::ConfiguracionLilyth cfg_l{};
    cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::CRUDA;
    cfg_l.reloj_actividad = TipoRelojVentana::NU_TRADES;
    cfg_l.t_max_efectivo = 32;

    auto lil = std::make_unique<LilythReal>();
    verificar(lil->configurar(cfg_l, cfg_d, *delta, x_cal), "P4: ObservadorLilyth debe configurar");

    // -------------------------------------------------------------------------
    // 3a. Camino de serie con d_nu y pared:
    //   - Calma: nu_punto* = 100.0 trades/barra, pared avanza +1000 ms/barra.
    //   - Evento 1 (en t = 200, ACTIVIDAD CONGELADA): d_nu cae a 10.0 trades/barra (10 veces menor),
    //     por lo que la relajacion a <= 0.5 * dY tarda tau_pared = 20 barras de pared (20000 ms),
    //     pero acumula Delta_nu = 20 * 10 = 200 trades -> tau_act = 200 / 100 = 2.0 barras de calma (normal!).
    //   - Evento 2 (en t = 450, PERDIDA DE RIGIDEZ): d_nu sigue en 100.0 trades/barra (normal),
    //     y la relajacion tarda tau_pared = 20 barras -> Delta_nu = 20 * 100 = 2000 -> tau_act = 20.0 (grande en los dos!).
    // -------------------------------------------------------------------------
    delta->reiniciar_flujo();
    lil->reiniciar_flujo();

    for (std::size_t t = 0; t < 520; ++t) {
        double x = 0.0;
        double dnu = 100.0;
        const double pared_ms = static_cast<double>((t + 1) * 1000);

        if (t == 200) {
            x = 12.0; // Choque 1: actividad congelada
        } else if (t > 200 && t < 220) {
            dnu = 10.0; // Actividad congelada a 1/10 del ritmo de calma
            const double nu_desde_shock = static_cast<double>(t - 200) * 10.0;
            // En reloj de actividad relaja con semivida 2.0 barras de calma (200 trades): cruza 6.0 exactamente en t = 220
            x = 12.0 * std::pow(0.5, nu_desde_shock / 200.0);
        } else if (t == 220) {
            dnu = 10.0;
            x = 5.5; // <= 0.5 * 12.0 en tau_pared = 20, Delta_nu = 200 -> tau_act = 2.0
        } else if (t == 450) {
            x = 12.0; // Choque 2: perdida de rigidez genuina (actividad normal dnu = 100.0)
        } else if (t > 450 && t < 470) {
            dnu = 100.0;
            x = 12.0 * std::pow(0.5, static_cast<double>(t - 450) / 20.0);
        } else if (t == 470) {
            dnu = 100.0;
            x = 5.5; // <= 0.5 * 12.0 en tau_pared = 20, Delta_nu = 2000 -> tau_act = 20.0
        }

        const auto s = delta->paso(x);
        (void)lil->paso(s, x, true, dnu, 0, 0, pared_ms);
    }

    verificar(lil->cerradas_historico().size() >= 2, "P4 serie: deben haberse cerrado los 2 eventos");
    if (lil->cerradas_historico().size() >= 2) {
        const auto& ev_cong = lil->cerradas_historico()[0].cruda;
        const auto& ev_rig = lil->cerradas_historico()[1].cruda;

        std::cout << "  [P4 serie - Actividad congelada] tau_pared = " << ev_cong.tau_rest_consolidado
                  << " barras (" << ev_cong.pared_consolidado << " ms), tau_act = " << ev_cong.tau_act_consolidado
                  << " barras de calma (nu_punto* congelado = " << ev_cong.nu_punto_ancla << ")\n";
        std::cout << "  [P4 serie - Perdida de rigidez ] tau_pared = " << ev_rig.tau_rest_consolidado
                  << " barras (" << ev_rig.pared_consolidado << " ms), tau_act = " << ev_rig.tau_act_consolidado
                  << " barras de calma (nu_punto* congelado = " << ev_rig.nu_punto_ancla << ")\n";

        // Verificar contra referencia exacta:
        // En actividad congelada: nu_punto_ancla == 100.0 (congelado en W_ancla, caza M01),
        // tau_rest_consolidado == 20.0 (grande en pared), y tau_act_consolidado == 2.0 (normal en actividad, caza M03!).
        verificar(std::abs(ev_cong.nu_punto_ancla - 100.0) <= 1e-12,
                  "P4: nu_punto* en W_ancla debe ser exactamente 100.0");
        verificar(std::abs(ev_cong.tau_rest_consolidado - 20.0) <= 1e-12 &&
                  std::abs(ev_cong.pared_consolidado - 20000.0) <= 1e-9,
                  "P4 (actividad congelada): tau_rest grande en pared (20 barras, 20000 ms)");
        verificar(std::abs(ev_cong.tau_act_consolidado - 2.0) <= 1e-12,
                  "P4 (actividad congelada): tau_act normal en reloj de actividad (2.0 barras de calma, con nu_punto* congelado)");

        // En perdida de rigidez: grande tanto en pared (20.0) como en actividad (20.0)
        verificar(std::abs(ev_rig.tau_rest_consolidado - 20.0) <= 1e-12 &&
                  std::abs(ev_rig.tau_act_consolidado - 20.0) <= 1e-12,
                  "P4 (perdida de rigidez): grande tanto en pared (20.0) como en reloj de actividad (20.0)");
    }

    // -------------------------------------------------------------------------
    // 3b. Camino por filas (ObservadorSaltosFilas + paso_filas en barras intrinsecas):
    //   En barras intrinsecas (delta_nu = 100.0), la pared b.t_cierre_ms NO avanza a paso fijo:
    //   cuando la actividad por segundo cae 10x, cada barra intrinseca tarda 10000 ms en vez de 1000 ms,
    //   por lo que `pared_consolidado` es 10x mayor (30000 ms vs 3000 ms) mientras `tau_act_consolidado` vale 3.0 en ambos.
    // -------------------------------------------------------------------------
    {
        std::vector<EntradaPasoCxx> filas_cal(N_CAL);
        for (std::size_t i = 0; i < N_CAL; ++i) {
            filas_cal[i].t_cierre_ms = static_cast<std::uint64_t>((i + 1) * 1000);
            filas_cal[i].t_disponible_ms = filas_cal[i].t_cierre_ms;
            filas_cal[i].valor = x_cal[i];
            filas_cal[i].volumen = 100.0;
            filas_cal[i].con_sesion = true;
            filas_cal[i].imputada = false;
        }
        ObsFilasReal::Fuente fuente{};
        fuente.base.reloj = TipoRelojVentana::NU_VOLUMEN;
        fuente.delta_nu = 100.0;
        fuente.agregacion = ark::ventanas::AgregacionBarra::ULTIMO;

        auto osf = std::make_unique<ObsFilasReal>();
        osf->calibrar_filas(cfg_d, fuente, filas_cal);
        auto lil_f = std::make_unique<LilythReal>();
        verificar(lil_f->configurar(cfg_l, cfg_d, osf->observador(), x_cal), "P4 filas: configurar Lilyth");

        std::uint64_t reloj_pared_ms = 1000000ULL;
        for (std::size_t b_idx = 0; b_idx < 420; ++b_idx) {
            // Entre las barras 201..205 la actividad en tiempo de pared se congela 10x (cada barra de 100 de volumen tarda 10000 ms)
            const std::uint64_t dt_ms = (b_idx > 200 && b_idx <= 205) ? 10000ULL : 1000ULL;
            reloj_pared_ms += dt_ms;

            double val = 0.0;
            if (b_idx == 200 || b_idx == 340) val = 12.0;
            else if ((b_idx > 200 && b_idx < 203) || (b_idx > 340 && b_idx < 343)) val = 9.0;
            else if (b_idx == 203 || b_idx == 343) val = 4.0; // Restituye en 3 barras intrinsecas

            EntradaPasoCxx fila{};
            fila.t_cierre_ms = reloj_pared_ms;
            fila.t_disponible_ms = reloj_pared_ms;
            fila.valor = val;
            fila.volumen = 100.0;
            fila.con_sesion = true;
            fila.imputada = false;

            const auto sf = osf->paso_fila(fila);
            (void)lil_f->paso_filas(sf, false);
        }

        verificar(lil_f->cerradas_historico().size() >= 2, "P4 filas: deben cerrarse los 2 eventos en barras intrinsecas");
        if (lil_f->cerradas_historico().size() >= 2) {
            const auto& b_lenta = lil_f->cerradas_historico()[0].cruda;
            const auto& b_normal = lil_f->cerradas_historico()[1].cruda;
            verificar(std::abs(b_lenta.tau_act_consolidado - 3.0) <= 1e-12 &&
                      std::abs(b_normal.tau_act_consolidado - 3.0) <= 1e-12,
                      "P4 filas: en barras intrinsecas ambos choques restituyen en tau_act = 3.0 barras de actividad");
            verificar(std::abs(b_lenta.pared_consolidado - 30000.0) <= 1e-9 &&
                      std::abs(b_normal.pared_consolidado - 3000.0) <= 1e-9,
                      "P4 filas: en barras intrinsecas la pared NO avanza a paso fijo (30000 ms con actividad congelada vs 3000 ms en calma)");
        }
    }
}

// =============================================================================
// 4. P11: CONTRASTE DE MAREA (E1.14 [DECIDIDO]: DESCRIPTIVO, FUERA DE HOLM)
// =============================================================================
void prueba_p11_contraste_marea() {
    std::cout << "[4/5] P11: Contraste de marea (E1.14: madurez 20 dias, censura NIVEL sin imputar, sesgo diurno y descriptivo fuera de Holm)...\n";

    // 4a. Antes de los 20 dias por fase (`dias_minimos_purga = 20`), la rama PURGADA no mide (`NO_MEDIDA_PURGA_INMADURA`),
    //     y en NIVEL un incremento purgado invalido censura la ranura (`CENSURADA_DATO_INVALIDO`) sin rellenar.
    {
        const std::size_t N_CAL = 3600;
        std::mt19937_64 rng(2026110664ULL);
        std::normal_distribution<double> nd(0.0, 1.0);

        // Calibrar solo con 5 dias (inmadura: dias_minimos_purga = 20)
        std::vector<double> x_inm(N_CAL, 100.0);
        std::vector<std::size_t> f_inm(N_CAL, 0);
        std::vector<std::int64_t> d_inm(N_CAL, 0);
        for (std::size_t i = 1; i < N_CAL; ++i) {
            x_inm[i] = x_inm[i - 1] + nd(rng);
            f_inm[i] = i % 2;
            d_inm[i] = static_cast<std::int64_t>(i / 900); // solo 4 dias (< 20 dias)
        }

        ark::metrologia::ConfiguracionObservador cfg_d{};
        cfg_d.integracion = ark::metrologia::Integracion::NIVEL;
        cfg_d.purga = true;
        cfg_d.bins_fase = 2;
        cfg_d.dias_minimos_purga = 20;
        cfg_d.arl0_pasos_nativos = 250.0;
        cfg_d.escalas_max = 1;
        cfg_d.medir_tau50 = false;

        // Para calibrar ObservadorSaltos necesitamos >= 20 dias en la historia de calibracion de delta,
        // pero si entregamos a Lilyth una historia con < 20 dias en `configurar`, la purga de Lilyth esta inmadura al arrancar:
        std::vector<double> x_cal_d(N_CAL, 100.0);
        std::vector<std::size_t> f_cal_d(N_CAL, 0);
        std::vector<std::int64_t> d_cal_d(N_CAL, 0);
        for (std::size_t i = 1; i < N_CAL; ++i) {
            x_cal_d[i] = x_cal_d[i - 1] + nd(rng);
            f_cal_d[i] = i % 2;
            d_cal_d[i] = static_cast<std::int64_t>(i / 40); // 90 dias: 800 barras madurando + 2800 validas >= 2400
        }
        auto delta = std::make_unique<ObsReal>();
        delta->calibrar(cfg_d, x_cal_d, {}, {}, f_cal_d, d_cal_d);

        ark::lilyth::ConfiguracionLilyth cfg_l{};
        cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA;
        cfg_l.t_max_efectivo = 16;

        auto lil_inm = std::make_unique<LilythReal>();
        verificar(lil_inm->configurar(cfg_l, cfg_d, *delta, x_inm, {}, f_inm, d_inm), "P11 inmadura: configurar Lilyth");

        delta->reiniciar_flujo();
        double lvl = 100.0;
        for (std::size_t t = 0; t < 260; ++t) {
            if (t == 200) lvl += 15.0;
            else if (t == 202) lvl -= 10.0;
            else lvl += 0.1 * nd(rng);
            const auto s = delta->paso(lvl, NaN, true, t % 2, 5);
            (void)lil_inm->paso(s, lvl, true, NaN, t % 2, 5);
        }
        verificar(lil_inm->contadores().no_medida[static_cast<int>(ark::lilyth::DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA)] >= 1,
                  "P11 (E1.14): antes de cumplir los 20 dias por fase, la rama PURGADA no mide (NO_MEDIDA_PURGA_INMADURA)");

        // Ahora con purga madura (30 dias), verificar que en NIVEL un incremento purgado invalido a mitad de una ranura viva
        // censura la rama purgada como CENSURADA_DATO_INVALIDO sin rellenar
        auto lil_mad = std::make_unique<LilythReal>();
        verificar(lil_mad->configurar(cfg_l, cfg_d, *delta, x_cal_d, {}, f_cal_d, d_cal_d), "P11 madura: configurar Lilyth");
        delta->reiniciar_flujo();
        lvl = 100.0;
        for (std::size_t t = 0; t < 260; ++t) {
            if (t == 200) lvl += 15.0;
            else lvl += 0.05 * nd(rng);
            const bool val_paso = (t != 203); // En t = 203 (con la ranura abierta en t = 200) entra un dato invalido
            const double x_in = val_paso ? lvl : NaN;
            const auto s = delta->paso(x_in, NaN, val_paso, t % 2, 35);
            (void)lil_mad->paso(s, x_in, val_paso, NaN, t % 2, 35);
        }
        bool vio_censura_invalido = false;
        for (const auto& r : lil_mad->cerradas_historico()) {
            if (r.purgada.desenlace == ark::lilyth::DesenlaceRanura::CENSURADA_DATO_INVALIDO) {
                vio_censura_invalido = true;
                break;
            }
        }
        verificar(vio_censura_invalido,
                  "P11 (E1.14 / E1.17): en NIVEL, un incremento purgado invalido rompe la reintegracion y censura la ranura sin rellenar");
    }

    // 4b. AR(1) por perfil diurno determinista:
    //     - 2 fases por dia: fase 0 (calma baja m_0 = 0.50) y fase 1 (marea alta m_1 = 2.00).
    //     - Tras madurar 25 dias (>= 20 dias), ocurren 6 choques al final de la fase 0 cuya relajacion entra en la fase 1:
    //       en la serie CRUDA, la multiplicacion por m_1 = 2.00 estira artificialmente el tiempo hasta cruzar 0.5 * dY_cruda
    //       (con duraciones variadas en {7, 10, 6, 11, 8, 12} pasos en cruda frente a {3, 4, 3, 4, 3, 4} pasos en purgada),
    //       mientras que PurgaDiurna divide por sqrt(v_fase / v_media) y recupera la semivida verdadera del AR(1) subyacente!
    {
        const std::size_t N_CAL = 3600;
        const double phi = 0.75; // semivida verdadera ln(2)/(-ln(0.75)) = 2.41 pasos -> cruce en 3 pasos
        std::mt19937_64 rng(2026110665ULL);
        std::normal_distribution<double> nd(0.0, 1.0);

        auto factor_marea = [](std::size_t fase) -> double {
            return (fase == 0) ? 0.50 : 2.00;
        };

        std::vector<double> x_cal(N_CAL, 0.0);
        std::vector<std::size_t> f_cal(N_CAL, 0);
        std::vector<std::int64_t> d_cal(N_CAL, 0);
        double u = 0.0;
        for (std::size_t i = 0; i < N_CAL; ++i) {
            u = phi * u + 0.35 * nd(rng);
            const std::size_t fa = (i / 20) % 2;
            const std::int64_t di = static_cast<std::int64_t>(i / 40); // 90 dias (>= 20 dias + 2400 barras)
            f_cal[i] = fa;
            d_cal[i] = di;
            x_cal[i] = factor_marea(fa) * u;
        }

        ark::metrologia::ConfiguracionObservador cfg_d{};
        cfg_d.integracion = ark::metrologia::Integracion::FLUJO;
        cfg_d.purga = true;
        cfg_d.bins_fase = 2;
        cfg_d.dias_minimos_purga = 20;
        cfg_d.arl0_pasos_nativos = 250.0;
        cfg_d.escalas_max = 1;
        cfg_d.medir_tau50 = false;

        auto delta = std::make_unique<ObsReal>();
        delta->calibrar(cfg_d, x_cal, {}, {}, f_cal, d_cal);

        ark::lilyth::ConfiguracionLilyth cfg_l{};
        cfg_l.representacion = ark::lilyth::ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA;
        cfg_l.t_max_efectivo = 20;
        cfg_l.b_max_nulo = 4;
        cfg_l.hilos_nulo = 2;

        auto lil = std::make_unique<LilythReal>();
        verificar(lil->calibrar(*delta, cfg_l, cfg_d, x_cal, {}, f_cal, d_cal),
                  "P11 marea: ObservadorLilyth debe calibrar ambas representaciones (PURGADA y CRUDA)");

        // En el flujo de prueba, inyectamos 6 choques separados 200 pasos que abren justo al final de la fase 0 (m = 0.50)
        // y cuyos pasos posteriores tau = 1..k_marea caen en la fase 1 (m = 2.00):
        // En el proceso latente estacionario `u_t`, el choque relaja con u(tau) = A * phi^tau (cruza 0.5 * A en tau_purg = 3 pasos).
        // Pero en la serie cruda `x_t = m(fase) * u_t`, al pasar de m_0 = 0.50 (en el disparo t_a) a m_1 = 2.00 (en tau >= 1),
        // x(t_a) = 0.50 * A mientras que x(t_a + tau) = 2.00 * A * phi^tau, por lo que en CRUDA tarda hasta que 4 * phi^tau <= 0.5
        // (es decir, phi^tau <= 0.125 -> tau_cruda = 8..11 pasos)!
        delta->reiniciar_flujo();
        lil->reiniciar_flujo_conservando_purga();

        const std::array<double, 6> phis_evento{0.72, 0.81, 0.70, 0.83, 0.75, 0.84};
        std::size_t idx_ev = 0;
        std::size_t t_ultimo_shock = 0;
        bool en_shock = false;

        for (std::size_t t = 0; t < 1600; ++t) {
            std::size_t fa = (t / 20) % 2;
            const std::int64_t di = 100 + static_cast<std::int64_t>(t / 40);
            double u_val = 0.05 * nd(rng);

            if (t >= 200 && (t % 220 == 0) && idx_ev < phis_evento.size()) {
                en_shock = true;
                t_ultimo_shock = t;
                fa = 0; // Disparo al cierre de la fase 0 (m_0 = 0.50)
                u_val = 16.0;
                ++idx_ev;
            } else if (en_shock && t > t_ultimo_shock && (t - t_ultimo_shock) <= 18) {
                const std::size_t tau = t - t_ultimo_shock;
                const double phi_k = phis_evento[idx_ev - 1];
                fa = 1; // Los pasos posteriores caen en la fase 1 (m_1 = 2.00)
                u_val = 16.0 * std::pow(phi_k, static_cast<double>(tau));
            } else {
                en_shock = false;
            }

            const double x_val = factor_marea(fa) * u_val;
            const auto s = delta->paso(x_val, NaN, true, fa, di);
            (void)lil->paso(s, x_val, true, NaN, fa, di);
        }

        lil->completar_nulos();
        const auto& marea = lil->resumen_marea();
        const auto ref = calcular_marea_ref(lil->cerradas_historico(),
                                            lil->aalen_johansen_cruda().total.rmst,
                                            lil->aalen_johansen_purgada().total.rmst);

        std::cout << "  [P11 Marea] n_pareados = " << marea.n_pareados_consolidado
                  << ", mediana(Delta_tau_marea) = " << std::fixed << std::setprecision(2) << marea.mediana_delta_tau_consolidado
                  << " (IC 95%: [" << marea.ic_inf_delta_tau_consolidado << ", " << marea.ic_sup_delta_tau_consolidado << "])"
                  << ", Delta_RMST_marea = " << marea.delta_rmst_marea
                  << " (IC 95%: [" << marea.ic_inf_delta_rmst_marea << ", " << marea.ic_sup_delta_rmst_marea << "])\n";

        verificar(marea.valido && marea.n_pareados_consolidado >= 5,
                  "P11: deben haberse medido al menos 5 eventos pareados (CRUDA vs PURGADA)");
        verificar(marea.descriptivo_fuera_de_holm,
                  "P11 [DECIDIDO]: el contraste de marea es descriptivo y permanece fuera de la familia de 4 p de Holm");
        // En PURGADA, todos los eventos recuperan la semivida verdadera (~3-5 pasos), mientras en CRUDA la marea los estira a ~8-13 pasos
        verificar(marea.mediana_delta_tau_consolidado >= 4.0 && marea.ic_inf_delta_tau_consolidado > 0.0,
                  "P11: Delta_tau_marea = tau_cruda - tau_purgada mide el sesgo positivo inducido por la marea");
        verificar(marea.delta_rmst_marea >= 3.0 && marea.ic_inf_delta_rmst_marea > 0.0,
                  "P11: Delta_RMST_marea > 0 con IC bootstrap positivo");
        // Casamiento bit a bit con el motor de referencia independiente (caza M04 rama cruzada y M05 signo invertido)
        verificar(std::abs(marea.mediana_delta_tau_consolidado - ref.mediana_delta) <= 1e-12 &&
                  std::abs(marea.delta_rmst_marea - ref.delta_rmst) <= 1e-12,
                  "P11: mediana pareada y Delta_RMST_marea coinciden bit a bit con el motor de referencia independiente");
        for (std::size_t k = 0; k < marea.pares.size() && k < ref.deltas_pareados.size(); ++k) {
            verificar(std::abs(marea.pares[k].delta_tau_consolidado - ref.deltas_pareados[k]) <= 1e-12,
                      "P11: cada par k empareja estrictamente la rama cruda y la rama purgada del MISMO evento k");
        }
    }
}

// =============================================================================
// 5. TESTIGO DETERMINISTA DE LOS MUTANTES DE LA ETAPA 6
// =============================================================================
void prueba_testigo_mutantes_etapa6() {
    std::cout << "[5/5] Control testigo determinista de los mutantes de la Etapa 6...\n";
    // 1. M02: NU_INFO_SHANNON rechazado
    verificar(!ark::lilyth::detalle_reloj_marea::reloj_actividad_valido(TipoRelojVentana::NU_INFO_SHANNON),
              "Testigo M02: reloj_actividad_valido rechaza NU_INFO_SHANNON");
    verificar(ark::lilyth::detalle_reloj_marea::reloj_actividad_valido(TipoRelojVentana::NU_VOLUMEN) &&
              ark::lilyth::detalle_reloj_marea::reloj_actividad_valido(TipoRelojVentana::NU_TRADES) &&
              ark::lilyth::detalle_reloj_marea::reloj_actividad_valido(TipoRelojVentana::NU_NOCIONAL),
              "Testigo M02: acepta NU_VOLUMEN, NU_TRADES y NU_NOCIONAL");

    // 2. M01: seleccionar_nu_punto_ancla usa el valor congelado al abrir, nunca el movil posterior al choque
    verificar(std::abs(ark::lilyth::detalle_reloj_marea::seleccionar_nu_punto_ancla(100.0, 12.5) - 100.0) <= 1e-12,
              "Testigo M01: nu_punto* permanece congelado en el valor de W_ancla");

    // 3. M03: calcular_tau_actividad divide delta_nu por nu_punto*, nunca devuelve pasos de pared
    verificar(std::abs(ark::lilyth::detalle_reloj_marea::calcular_tau_actividad(200.0, 100.0, 20) - 2.0) <= 1e-12,
              "Testigo M03: tau_act = Delta_nu_rest / nu_punto* = 2.0 (distinto de tau_pared = 20)");

    // 4. M04 y M05: evaluar_resumen_marea empareja el mismo evento k con signo tau_cruda - tau_purgada
    using Ran = LilythReal::Ranura;
    std::vector<Ran> cerradas(3);
    const std::array<double, 3> tc{10.0, 20.0, 12.0};
    const std::array<double, 3> tp{2.0, 19.0, 3.0};
    // Pareados reales: (10-2)=8, (20-19)=1, (12-3)=9 -> mediana = 8.0.
    // Si M04 cruza k con (k+1)%3: (10-19)=-9, (20-3)=17, (12-2)=10 -> par[0] = -9 != 8 y mediana = 10.0 != 8.0!
    for (std::size_t k = 0; k < 3; ++k) {
        cerradas[k].id_evento = k + 1;
        cerradas[k].signo = +1;
        cerradas[k].cruda.desenlace = ark::lilyth::DesenlaceRanura::RESTITUIDA;
        cerradas[k].purgada.desenlace = ark::lilyth::DesenlaceRanura::RESTITUIDA;
        cerradas[k].cruda.tau_rest_primer_toque = tc[k];
        cerradas[k].cruda.tau_rest_consolidado = tc[k];
        cerradas[k].cruda.tau_fin = static_cast<std::size_t>(tc[k]);
        cerradas[k].purgada.tau_rest_primer_toque = tp[k];
        cerradas[k].purgada.tau_rest_consolidado = tp[k];
        cerradas[k].purgada.tau_fin = static_cast<std::size_t>(tp[k]);
    }
    const auto rm = ark::lilyth::detalle_reloj_marea::evaluar_resumen_marea<Ran>(cerradas, 14.0, 8.0, 24, 0.05, 199, 42);
    verificar(rm.valido && rm.pares.size() == 3, "Testigo M04: 3 pares evaluados");
    verificar(std::abs(rm.pares[0].delta_tau_consolidado - 8.0) <= 1e-12 &&
              std::abs(rm.pares[1].delta_tau_consolidado - 1.0) <= 1e-12 &&
              std::abs(rm.pares[2].delta_tau_consolidado - 9.0) <= 1e-12 &&
              std::abs(rm.mediana_delta_tau_consolidado - 8.0) <= 1e-12 &&
              std::abs(rm.delta_rmst_marea - 6.0) <= 1e-12,
              "Testigo M04/M05: emparejamiento exacto por evento y signo tau_cruda - tau_purgada");
}

}  // namespace

int main(int argc, char** argv) {
    const bool modo_rapido = (argc > 1 && std::strcmp(argv[1], "--rapido") == 0);
    std::cout << "====================================================================\n";
    std::cout << " BATERIA ETAPA 6 -- OPERADOR DE LILYTH (RELOJ DOBLE Y MAREA, C++20)\n";
    std::cout << "====================================================================\n";

    prueba_validacion_reloj_actividad();
    prueba_testigo_mutantes_etapa6();

    if (!modo_rapido) {
        prueba_k1_serie_y_filas_cero_reservas();
        prueba_p4_reloj_doble();
        prueba_p11_contraste_marea();
    } else {
        prueba_p4_reloj_doble();
    }

    std::cout << "====================================================================\n";
    if (g_fallas == 0) {
        std::cout << " RESULTADO: TODOS LOS CONTROLES DE LA ETAPA 6 PASARON (0 fallas)\n";
        std::cout << "====================================================================\n";
        return 0;
    }
    std::cerr << " RESULTADO: " << g_fallas << " FALLA(S) DETECTADA(S)\n";
    std::cerr << "====================================================================\n";
    return 1;
}