// prueba_lilyth_etapa1.cpp -- Bateria de verificacion y mutantes de la Etapa 1 del Operador de Lilyth.
// Cubre los 6 criterios exigidos por Code:
//   1. Identidad: el operador delta con Lilyth y sin ella da las mismas salidas, bit a bit.
//   2. Causalidad: truncado y perturbacion del futuro + mutante de 1 paso al futuro (que DEBE fallar).
//   3. eta(tau) y los dos tau (primer_toque, consolidado) + contraste de marea y censura en NIVEL sobre casos a mano.
//   4. Aalen-Johansen contra una referencia independiente escrita aparte, con replicas que cortan trayectorias.
//   5. Desborde de ranuras contado y CERO reservas de memoria dinamica en paso (operator new interceptado).
//   6. Rechazo de configuraciones invalidas (NULO_NINGUNO y NULO_BLOQUES en nulo_sin_memoria, E1.2 y E1.13).

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <new>
#include <vector>

#include "observador_lilyth.hpp"

namespace {

std::size_t g_reservas_ciclo_caliente = 0;
bool g_contar_reservas = false;

}  // namespace

void* operator new(std::size_t sz) {
    if (g_contar_reservas) ++g_reservas_ciclo_caliente;
    if (void* p = std::malloc(sz ? sz : 1)) return p;
    throw std::bad_alloc{};
}

void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

namespace {

using namespace ark::lilyth;
using namespace ark::metrologia;

// Operador delta sintetico / determinista para probar la envoltura bit a bit y los disparos exactos.
struct OperadorDeltaPrueba {
    struct Salida {
        bool emitida = false;
        double valor_interno = NaN;
        double estadistico_g = NaN;
        bool hay_evento = false;
        FichaEvento evento{};
    };

    std::int64_t paso_idx = -1;
    double acum = 0.0;
    std::int64_t t_disparo_1 = -1;
    std::int64_t t_disparo_2 = -1;
    double dy_disparo = 4.0;

    Salida paso(double y, bool valido) noexcept {
        ++paso_idx;
        Salida s{};
        s.emitida = valido;
        if (!valido || !std::isfinite(y)) return s;
        acum = 0.85 * acum + 0.15 * y;
        s.valor_interno = acum;
        s.estadistico_g = y * y;
        if (paso_idx == t_disparo_1 || paso_idx == t_disparo_2) {
            s.hay_evento = true;
            s.evento.estado = EstadoSaltos::SALTO;
            s.evento.canal = 0;
            s.evento.escala_origen = 0;
            s.evento.alarma.hay = true;
            s.evento.alarma.signo = (y >= 0.0) ? 1 : -1;
            s.evento.alarma.t = paso_idx;
            s.evento.alarma.l = 4;
            s.evento.alarma.t_cambio = paso_idx - 3;
            s.evento.alarma.amplitud = dy_disparo;
            s.evento.alarma.se_amplitud = 0.5;
            s.evento.alarma.g = 25.0;
            s.evento.ficha_origen.j = 0;
            s.evento.ficha_origen.m = 1;
            s.evento.ficha_origen.nulo.tau_int = 2.0;
            s.evento.ficha_origen.det.L = 16;   // W_ancla = 16, guarda = 2
        }
        return s;
    }
};

// Referencia independiente de Aalen-Johansen escrita desde cero por recorrido explicito de sujetos.
struct EventoRefAJ {
    int tipo = 0;   // 1 = restituido, 2 = replica, 0 = censurado
    std::size_t t = 0;
};

struct SalidaRefAJ {
    std::vector<std::size_t> n_riesgo;
    std::vector<double> S;
    std::vector<double> cif_rest;
    std::vector<double> cif_rep;
    double rmst = 0.0;
};

SalidaRefAJ aalen_johansen_referencia(const std::vector<EventoRefAJ>& evs, std::size_t tmax) {
    SalidaRefAJ out{};
    out.n_riesgo.assign(tmax, 0);
    out.S.assign(tmax, 1.0);
    out.cif_rest.assign(tmax, 0.0);
    out.cif_rep.assign(tmax, 0.0);
    double s_prev = 1.0, cr = 0.0, cp = 0.0, rmst = 0.0;
    for (std::size_t tau = 1; tau <= tmax; ++tau) {
        std::size_t n_at_risk = 0, d1 = 0, d2 = 0;
        for (const auto& e : evs) {
            const std::size_t tf = std::min(e.t, tmax);
            if (tf >= tau) ++n_at_risk;
            if (tf == tau && e.tipo == 1) ++d1;
            if (tf == tau && e.tipo == 2) ++d2;
        }
        out.n_riesgo[tau - 1] = n_at_risk;
        if (n_at_risk > 0) {
            cr += s_prev * (static_cast<double>(d1) / static_cast<double>(n_at_risk));
            cp += s_prev * (static_cast<double>(d2) / static_cast<double>(n_at_risk));
            s_prev *= static_cast<double>(n_at_risk - d1 - d2) / static_cast<double>(n_at_risk);
        }
        out.S[tau - 1] = s_prev;
        out.cif_rest[tau - 1] = cr;
        out.cif_rep[tau - 1] = cp;
        rmst += (1.0 - cr);
    }
    out.rmst = rmst;
    return out;
}

bool iguales_bit_a_bit(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return true;
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

}  // namespace

int main() {
    int fallos = 0;
    auto check = [&](bool cond, const char* msg) {
        if (!cond) {
            std::cerr << "[FALLO] " << msg << "\n";
            ++fallos;
        } else {
            std::cout << "[OK]    " << msg << "\n";
        }
    };

    // -------------------------------------------------------------------------
    // Control 0: Validacion de configuracion (E1.2 y E1.13)
    // -------------------------------------------------------------------------
    {
        ObservadorLilyth<OperadorDeltaPrueba, 4, 64> obs;
        ConfiguracionLilyth cfg{};
        cfg.nulo_sin_memoria = LILYTH_NULO_NINGUNO;
        check(!obs.configurar(cfg), "C0a: rechaza NULO_NINGUNO en nulo_sin_memoria (E1.2)");
        cfg.nulo_sin_memoria = LILYTH_NULO_BLOQUES;
        check(!obs.configurar(cfg), "C0b: rechaza NULO_BLOQUES en nulo_sin_memoria (E1.13)");
        cfg.nulo_sin_memoria = LILYTH_NULO_SIGNOS;
        cfg.t_max_efectivo = 16;
        check(obs.configurar(cfg), "C0c: acepta NULO_SIGNOS con t_max_efectivo valido");
    }

    // -------------------------------------------------------------------------
    // Control 1: Identidad bit a bit del operador delta solo vs con Lilyth
    // y CERO reservas de memoria dinamica en el ciclo caliente
    // -------------------------------------------------------------------------
    {
        OperadorDeltaPrueba solo{};
        solo.t_disparo_1 = 30;
        solo.t_disparo_2 = 55;

        ObservadorLilyth<OperadorDeltaPrueba, 4, 32> con_lilyth{};
        ConfiguracionLilyth cfg{};
        cfg.nulo_sin_memoria = LILYTH_NULO_PERMUTACION;
        cfg.t_max_efectivo = 12;
        con_lilyth.configurar(cfg);
        con_lilyth.operador_delta().t_disparo_1 = 30;
        con_lilyth.operador_delta().t_disparo_2 = 55;

        bool id_ok = true;
        g_reservas_ciclo_caliente = 0;
        g_contar_reservas = true;

        for (int t = 0; t < 100; ++t) {
            const double y = (t >= 27 && t <= 32) ? (10.0 - 1.5 * (t - 30)) : ((t % 2 == 0) ? 0.25 : -0.25);
            const auto s_solo = solo.paso(y, true);
            const auto s_conv = con_lilyth.operador_delta().paso(y, true);
            ObservadorLilyth<OperadorDeltaPrueba, 4, 32>::SalidaPaso s_lil{};
            con_lilyth.paso_con_salida_delta(y, true, 1.0, true, 2.0, s_conv, s_lil);

            if (s_solo.emitida != s_conv.emitida ||
                !iguales_bit_a_bit(s_solo.valor_interno, s_conv.valor_interno) ||
                !iguales_bit_a_bit(s_solo.estadistico_g, s_conv.estadistico_g) ||
                s_solo.hay_evento != s_conv.hay_evento) {
                id_ok = false;
            }
        }
        g_contar_reservas = false;
        check(id_ok, "C1a: identidad bit a bit del operador delta con y sin Lilyth (100 pasos)");
        check(g_reservas_ciclo_caliente == 0, "C1b: cero reservas de memoria dinamica en el ciclo caliente (operator new == 0)");
    }

    // -------------------------------------------------------------------------
    // Control 2: Caso a mano de eta(tau), tau_rest_primer_toque, tau_rest_consolidado
    // (con calor activo), reloj de actividad y contraste de marea (E1.3, E1.4, E1.6, E1.11, E1.14)
    // -------------------------------------------------------------------------
    {
        ObservadorLilyth<OperadorDeltaPrueba, 4, 16> lil{};
        ConfiguracionLilyth cfg{};
        cfg.integracion = Integracion::FLUJO;
        cfg.representacion = ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA;
        cfg.nulo_sin_memoria = LILYTH_NULO_PERMUTACION;
        cfg.t_max_efectivo = 8;
        lil.configurar(cfg);
        lil.operador_delta().t_disparo_1 = 25;

        // t = 0 .. 24: base constante Y = 10.0 en cruda; factor de purga = 0.5 -> Y_purgada = 5.0.
        // Velocidad del reloj de actividad en el ancla: d_nu = 4.0 por barra -> nu_punto* = 4.0.
        for (int t = 0; t < 25; ++t) {
            lil.paso(10.0, true, 0.5, true, 4.0);
        }

        // Activamos calor antes del disparo en t = 25:
        lil.fijar_calor_activo(true);

        // t = 25 (disparo t_a): Y_cruda(25) = 18.0 (dY_cruda = +8.0), factor_purga = 0.5 -> Y_purgada(25) = 9.0 (dY_purgada = +4.0).
        auto s25 = lil.paso(18.0, true, 0.5, true, 4.0);
        check(s25.abrio_ranura, "C2a: abre ranura en t_a = 25");

        // tau = 1 (t = 26): Y_cruda = 16.0 (eta_cruda = 6/8 = 0.75); factor_purga = 0.5 -> Y_purgada = 8.0 (eta_purgada = 3/4 = 0.75).
        lil.paso(16.0, true, 0.5, true, 8.0);

        // tau = 2 (t = 27): con calor aun activo, Y_cruda = 14.0 (eta_cruda = 4/8 = 0.50 -> PRIMER TOQUE en tau = 2, pero NO consolidado
        // porque hay calor!). En purgada cambiamos factor_purga = 0.6 -> Y_purgada = 8.4 (eta_purgada = 3.4/4 = 0.85 > 0.5, la marea
        // hizo parecer que la cruda restituia antes que la purgada).
        lil.paso(14.0, true, 0.6, true, 8.0);

        // tau = 3 (t = 28): se apaga el calor; Y_cruda = 13.0 (eta_cruda = 3/8 = 0.375 <= 0.5 -> CONSOLIDADO en tau = 3);
        // factor_purga = 0.5 -> Y_purgada = 6.5 (eta_purgada = 1.5/4 = 0.375 <= 0.5 -> primer toque y consolidado en tau = 3).
        lil.fijar_calor_activo(false);
        lil.paso(13.0, true, 0.5, true, 8.0);

        // tau = 4 (t = 29): sobrepaso al otro lado del ancla: Y_cruda = 8.0 (eta_cruda = -2/8 = -0.25 < 0 -> armonico_sospechoso = true).
        lil.paso(8.0, true, 0.5, true, 4.0);

        // Corremos tau = 5..8 para que cierre en t_max_efectivo = 8:
        ObservadorLilyth<OperadorDeltaPrueba, 4, 16>::SalidaPaso s_cierre{};
        for (int tau = 5; tau <= 8; ++tau) {
            s_cierre = lil.paso(10.0, true, 0.5, true, 4.0);
        }

        check(s_cierre.ranuras_cerradas_en_paso == 1, "C2b: la ranura cierra al cumplir t_max_efectivo = 8");
        const auto& rc = s_cierre.cerradas[0];
        check( std::abs(rc.cruda.y_ancla - 10.0) < 1e-12 && std::abs(rc.cruda.dy - 8.0) < 1e-12,
               "C2c: ancla Y*_cruda = 10.0 y dY_cruda = 8.0 exactos" );
        check( std::abs(rc.purgada.y_ancla - 5.0) < 1e-12 && std::abs(rc.purgada.dy - 4.0) < 1e-12,
               "C2d: ancla Y*_purgada = 5.0 y dY_purgada = 4.0 exactos" );
        check( std::abs(rc.cruda.tau_rest_primer_toque - 2.0) < 1e-12 &&
               std::abs(rc.cruda.tau_rest_consolidado - 3.0) < 1e-12,
               "C2e: cruda separa tau_rest_primer_toque (2, con calor) de tau_rest_consolidado (3, sin calor)" );
        // En tau = 1 y 2 el reloj de actividad sumo 8 + 8 = 16 unidades; con nu_punto* = 4, tau_act_primer_toque = 16 / 4 = 4.0 barras de calma:
        check( std::abs(rc.cruda.tau_act_primer_toque - 4.0) < 1e-12 &&
               std::abs(rc.cruda.tau_act_consolidado - 6.0) < 1e-12,
               "C2f: reloj de actividad normalizado tau_act = Delta_nu / nu_punto* exacto (4.0 y 6.0)" );
        check( std::abs(rc.delta_tau_marea_primer_toque - (-1.0)) < 1e-12,
               "C2g: diferencial de marea Delta_tau_marea = tau_cruda(2) - tau_purgada(3) = -1.0 exacto" );
        check( rc.cruda.armonico_sospechoso && std::abs(rc.cruda.min_eta_desde_2 - (-0.25)) < 1e-12,
               "C2h: etiqueta armonico_sospechoso activada en tau = 4 con min_eta = -0.25 sin vetar el evento" );
    }

    // -------------------------------------------------------------------------
    // Control 3: Integracion::NIVEL y censura cuando se rompe la purga (E1.14)
    // -------------------------------------------------------------------------
    {
        ObservadorLilyth<OperadorDeltaPrueba, 4, 16> lil_niv{};
        ConfiguracionLilyth cfg{};
        cfg.integracion = Integracion::NIVEL;
        cfg.representacion = ModoRepresentacionLilyth::AMBAS_CONTRASTE_MAREA;
        cfg.nulo_sin_memoria = LILYTH_NULO_SIGNOS;
        cfg.t_max_efectivo = 6;
        lil_niv.configurar(cfg);
        lil_niv.operador_delta().t_disparo_1 = 25;

        for (int t = 0; t < 25; ++t) lil_niv.paso(100.0, true, 1.0, true, 1.0);
        lil_niv.paso(110.0, true, 1.0, true, 1.0);   // disparo en t = 25
        lil_niv.paso(108.0, true, 1.0, true, 1.0);   // tau = 1: valido
        // En tau = 2 entra un incremento purgado invalido (purga_madura = false), mientras que el crudo sigue valido y restituye:
        lil_niv.paso(104.0, true, NaN, false, 1.0);  // tau = 2: cruda llega a 104 <= 105 (restituye), purgada se rompe!
        ObservadorLilyth<OperadorDeltaPrueba, 4, 16>::SalidaPaso sc{};
        for (int tau = 3; tau <= 6; ++tau) sc = lil_niv.paso(102.0, true, 1.0, true, 1.0);

        check(sc.ranuras_cerradas_en_paso == 1, "C3a: ranura NIVEL cierra correctamente");
        check(sc.cerradas[0].purgada.desenlace == DesenlaceRanura::CENSURADA_DATO_INVALIDO &&
              sc.cerradas[0].purgada.tau_fin == 1,
              "C3b: en NIVEL un incremento purgado invalido rompe la reintegracion y censura la rama purgada en tau = 1 sin rellenar (E1.14)");
        check(sc.cerradas[0].cruda.desenlace == DesenlaceRanura::RESTITUIDA &&
              std::abs(sc.cerradas[0].cruda.tau_rest_consolidado - 2.0) < 1e-12,
              "C3c: la rama cruda del mismo evento restituye normalmente en tau = 2");
    }

    // -------------------------------------------------------------------------
    // Control 4: Aalen-Johansen contra referencia escrita aparte + replicas que cortan trayectorias
    // -------------------------------------------------------------------------
    {
        constexpr std::size_t TMAX = 10;
        ResumenAalenJohansenSigno<TMAX> aj{};
        aj.limpiar(TMAX);

        // 8 eventos con mezcla de restituciones, cortes por replica y censuras, en signos + y -:
        struct Caso { int signo; DesenlaceRanura des; std::size_t tau; int tipo_ref; };
        const std::array<Caso, 8> casos{{
            {+1, DesenlaceRanura::RESTITUIDA, 2, 1},
            {+1, DesenlaceRanura::CORTADA_POR_REPLICA, 2, 2},
            {-1, DesenlaceRanura::RESTITUIDA, 4, 1},
            {-1, DesenlaceRanura::CENSURADA_DATO_INVALIDO, 3, 0},
            {+1, DesenlaceRanura::CORTADA_POR_REPLICA, 5, 2},
            {-1, DesenlaceRanura::RESTITUIDA, 5, 1},
            {+1, DesenlaceRanura::CENSURADA_HORIZONTE, 10, 0},
            {-1, DesenlaceRanura::RESTITUIDA, 7, 1}
        }};

        std::vector<EventoRefAJ> ref_tot, ref_pos, ref_neg;
        for (const auto& c : casos) {
            aj.registrar(c.signo, c.des, c.tau);
            ref_tot.push_back({c.tipo_ref, c.tau});
            if (c.signo > 0) ref_pos.push_back({c.tipo_ref, c.tau});
            else ref_neg.push_back({c.tipo_ref, c.tau});
        }

        const auto rt = aalen_johansen_referencia(ref_tot, TMAX);
        const auto rp = aalen_johansen_referencia(ref_pos, TMAX);
        const auto rn = aalen_johansen_referencia(ref_neg, TMAX);

        bool aj_ok = std::abs(aj.total.rmst - rt.rmst) < 1e-12 &&
                     std::abs(aj.positivos.rmst - rp.rmst) < 1e-12 &&
                     std::abs(aj.negativos.rmst - rn.rmst) < 1e-12;
        for (std::size_t i = 0; i < TMAX; ++i) {
            if (aj.total.en_riesgo[i] != rt.n_riesgo[i] ||
                std::abs(aj.total.supervivientes[i] - rt.S[i]) > 1e-12 ||
                std::abs(aj.total.cif_rest[i] - rt.cif_rest[i]) > 1e-12 ||
                std::abs(aj.total.cif_replica[i] - rt.cif_rep[i]) > 1e-12 ||
                std::abs((aj.total.supervivientes[i] + aj.total.cif_rest[i] + aj.total.cif_replica[i]) - 1.0) > 1e-12) {
                aj_ok = false;
            }
        }
        check(aj_ok, "C4a: Aalen-Johansen (total, + y -, S, CIF_rest, CIF_replica, RMST e identidad S+CIF1+CIF2=1) exacto contra referencia");

        // Prueba en vivo de una replica que corta una ranura abierta antes de restituir:
        ObservadorLilyth<OperadorDeltaPrueba, 4, 16> lil_rep{};
        ConfiguracionLilyth cfg{};
        cfg.nulo_sin_memoria = LILYTH_NULO_PERMUTACION;
        cfg.t_max_efectivo = 10;
        cfg.cortar_abiertas_con_replica = true;
        lil_rep.configurar(cfg);
        lil_rep.operador_delta().t_disparo_1 = 25;
        lil_rep.operador_delta().t_disparo_2 = 28;   // replica 3 pasos despues!

        for (int t = 0; t < 25; ++t) lil_rep.paso(0.0, true, 1.0, true, 1.0);
        lil_rep.paso(10.0, true, 1.0, true, 1.0);    // t = 25: abre ranura 0 (Y* = 0, dY = 10)
        lil_rep.paso(9.0, true, 1.0, true, 1.0);     // t = 26 (tau = 1): eta = 0.9 (sigue lejos de 0.5)
        lil_rep.paso(8.5, true, 1.0, true, 1.0);     // t = 27 (tau = 2): eta = 0.85
        auto s_rep = lil_rep.paso(15.0, true, 1.0, true, 1.0);   // t = 28 (tau = 3): llega REPLICA y corta la ranura 0!

        check(s_rep.ranuras_cerradas_en_paso == 1 &&
              s_rep.cerradas[0].cruda.desenlace == DesenlaceRanura::CORTADA_POR_REPLICA &&
              s_rep.cerradas[0].cruda.tau_fin == 3 &&
              lil_rep.aalen_johansen_cruda().total.n_cortados_replica == 1,
              "C4b: replica en t = 28 corta la ranura abierta en tau = 3 y suma 1 en d_replica[3] de Aalen-Johansen");
    }

    // -------------------------------------------------------------------------
    // Control 5: Desborde de ranuras contado (D8)
    // -------------------------------------------------------------------------
    {
        ObservadorLilyth<OperadorDeltaPrueba, 1, 16> lil_1{};   // RANURAS = 1
        ConfiguracionLilyth cfg{};
        cfg.nulo_sin_memoria = LILYTH_NULO_PERMUTACION;
        cfg.t_max_efectivo = 10;
        cfg.cortar_abiertas_con_replica = false;   // sin cortar la abierta para forzar que el banco de 1 ranura este lleno
        lil_1.configurar(cfg);
        lil_1.operador_delta().t_disparo_1 = 25;
        lil_1.operador_delta().t_disparo_2 = 27;

        for (int t = 0; t < 25; ++t) lil_1.paso(0.0, true);
        lil_1.paso(10.0, true);   // t = 25: ocupa la unica ranura
        lil_1.paso(9.0, true);    // t = 26
        auto s27 = lil_1.paso(12.0, true);   // t = 27: intenta abrir segunda ranura con RANURAS = 1 -> desborde!
        check(!s27.abrio_ranura && lil_1.desbordes_ranura() == 1, "C5: desborde de ranuras contado sin pisar ranura activa");
    }

    // -------------------------------------------------------------------------
    // Control 6: Causalidad (truncado / perturbacion del futuro) y Mutante de 1 paso al futuro
    // -------------------------------------------------------------------------
    {
        auto correr_serie = [](const std::vector<double>& serie, bool mutante_mira_futuro) {
            ObservadorLilyth<OperadorDeltaPrueba, 4, 16> lil{};
            ConfiguracionLilyth cfg{};
            cfg.nulo_sin_memoria = LILYTH_NULO_PERMUTACION;
            cfg.t_max_efectivo = 8;
            lil.configurar(cfg);
            lil.operador_delta().t_disparo_1 = 25;
            std::vector<double> huella;
            for (std::size_t t = 0; t < serie.size(); ++t) {
                double y_in = serie[t];
                if (mutante_mira_futuro && t + 1 < serie.size()) {
                    y_in += 1e-6 * serie[t + 1];   // mutante que espia 1 paso al futuro
                }
                const auto s = lil.paso(y_in, true, 1.0, true, 1.0);
                huella.push_back(static_cast<double>(s.abrio_ranura));
                for (const auto& r : lil.ranuras()) {
                    huella.push_back(r.ocupada ? r.cruda.y_ancla : -1.0);
                    huella.push_back(r.ocupada && r.cruda.pasos_observados > 0
                                         ? r.cruda.eta[r.cruda.pasos_observados - 1]
                                         : -2.0);
                }
            }
            return huella;
        };

        std::vector<double> s_base(40, 1.0), s_pert(40, 1.0);
        for (std::size_t t = 25; t < 40; ++t) {
            s_base[t] = 5.0 - 0.5 * static_cast<double>(t - 25);
            s_pert[t] = s_base[t];
        }
        // Perturbamos el futuro desde t = 30 en adelante:
        for (std::size_t t = 30; t < 40; ++t) s_pert[t] += 50.0;

        const auto h_base = correr_serie(s_base, false);
        const auto h_pert = correr_serie(s_pert, false);
        // En t = 0 .. 29 (primeros 30 * (1 + 4*2) = 270 doubles de la huella), deben coincidir bit a bit:
        const std::size_t n_pasado = 30 * (1 + 4 * 2);
        bool pasado_identico = true;
        for (std::size_t i = 0; i < n_pasado; ++i) {
            if (!iguales_bit_a_bit(h_base[i], h_pert[i])) pasado_identico = false;
        }
        check(pasado_identico, "C6a: causalidad estricta bit a bit (perturbar t >= 30 no altera nada en t <= 29)");

        // Ahora corremos el mutante de 1 paso al futuro y comprobamos que el arnes de causalidad LO CAZA en t = 29:
        const auto hm_base = correr_serie(s_base, true);
        const auto hm_pert = correr_serie(s_pert, true);
        bool mutante_cazado = false;
        for (std::size_t i = 0; i < n_pasado; ++i) {
            if (!iguales_bit_a_bit(hm_base[i], hm_pert[i])) mutante_cazado = true;
        }
        check(mutante_cazado, "C6b: el control de causalidad caza al mutante que mira 1 paso al futuro en t = 29");
    }

    if (fallos == 0) {
        std::cout << "\n=== ETAPA 1 LILYTH: TODOS LOS CONTROLES Y MUTANTES PASAN (0 FALLOS) ===\n";
        return 0;
    }
    std::cerr << "\n=== FALLOS DETECTADOS: " << fallos << " ===\n";
    return 1;
}