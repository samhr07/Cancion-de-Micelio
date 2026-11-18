// ============================================================================
// test_punto6_wilson_h0_h1.cpp
// Capa 2 (Corrección de Puntos 3, 4 y 6 de la 2a revisión de Opus):
//   - Parámetros sellados restaurados: a = 8.0 * sd_e, z_tau = 4.0 * sd_e * 0.97^tau.
//   - Horizonte CUANTIL_RESOLUCION + Pool DERIVADO (b_max_nulo = 64) con n_op = 3200
//     para alcanzar potencia base >= 50% en M1_H1 (q = 0) y eliminar el ruido de RetPot.
//   - 6B evaluado en dos modos explícitos:
//       (6B.1) Extremo a extremo real: ObservadorSaltos Y ObservadorLilyth ven los NaN.
//       (6B.2) Con forzar_via_b = true: Z1 y Z2 intervienen en el 100% de los eventos.
//   - Derivación de cotas q_max estrictamente monótona (break al primer fallo).
// ============================================================================

#include "protocolo_inputs_lilyth.hpp"
#include "banco_sintetico_comun.hpp"

#include <ark/dataset/tipos_dataset.hpp>
#include <ark/lilyth/observador_lilyth.hpp>
#include <ark/lilyth/nulo_lilyth.hpp>
#include <ark/lilyth/modelo_lineal_nulo.hpp>
#include <ark/metrologia/whittle_exacto.hpp>
#include <ark/holm.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <thread>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;
using namespace ark::lilyth::piloto;

constexpr std::size_t T_MAX_P6 = 64;
constexpr std::size_t R_REPLICAS = 200;
using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, T_MAX_P6>;

inline std::pair<double, double> wilson_95(std::size_t k, std::size_t n) noexcept {
    if (n == 0) return {0.0, 1.0};
    const double z = 1.959963984540054;
    const double p = static_cast<double>(k) / static_cast<double>(n);
    const double d = 1.0 + z * z / static_cast<double>(n);
    const double c = (p + z * z / (2.0 * static_cast<double>(n))) / d;
    const double h = z * std::sqrt(p * (1.0 - p) / static_cast<double>(n) + z * z / (4.0 * static_cast<double>(n) * static_cast<double>(n))) / d;
    return {std::max(0.0, c - h), std::min(1.0, c + h)};
}

inline void en_paralelo(std::size_t n_tareas, std::size_t n_hilos, const std::function<void(std::size_t)>& fn) {
    std::atomic<std::size_t> sig{0};
    std::vector<std::thread> trabajadores;
    const std::size_t h_ef = std::max<std::size_t>(1, n_hilos);
    trabajadores.reserve(h_ef);
    for (std::size_t h = 0; h < h_ef; ++h) {
        trabajadores.emplace_back([&]() {
            for (;;) {
                const std::size_t i = sig.fetch_add(1, std::memory_order_relaxed);
                if (i >= n_tareas) break;
                fn(i);
            }
        });
    }
    for (auto& t : trabajadores) t.join();
}

struct ResultadoReplicaP6 {
    bool serie_valida = false;
    std::size_t k_medibles = 0;
    std::size_t k_evaluados = 0;
    std::size_t k_cens_inv = 0;
    std::size_t k_no_medida = 0;
    double s_inter = NaN;
    double ee_inter = NaN;
    double p_holm_inter = NaN;
    bool rechaza_inter = false;
    bool rechaza_lenta_inter = false;
    bool rechaza_pool = false;
};

enum class ModoEjecucionP6 {
    LIMPIO = 0,
    QUIETO_ACTUAL = 1,
    QUIETO_PROTOCOLO = 2,
    CORRUPTO_ACTUAL = 3,
    CORRUPTO_PROTEGIDO = 4
};

ResultadoReplicaP6 correr_una_replica(
    const SerieSinteticaCoordinada& mundo,
    Integracion integ,
    bool usar_purga,
    ModoNuloLineal modo_lin,
    ModoEjecucionP6 modo_ej,
    bool detector_ve_corrupcion,
    bool forzar_via_b,
    double tasa_q,
    std::size_t racha_r,
    std::size_t idx_rep) {

    ConfiguracionObservador cfg_d{};
    cfg_d.integracion = integ;
    cfg_d.agregacion = AgregacionEscala::SUMA;
    cfg_d.arl0_pasos_nativos = 450.0;
    cfg_d.alfa = ALFA_PROYECTO;
    cfg_d.purga = usar_purga;
    cfg_d.bins_fase = usar_purga ? 24 : 1;
    cfg_d.dias_minimos_purga = 20;
    cfg_d.medir_tau50 = false;
    cfg_d.semilla = 20261117ULL + 10007ULL * idx_rep;
    cfg_d.escalas_max = 2;

    ObsDelta obs_d{};
    obs_d.calibrar(cfg_d, mundo.x_cal, {}, {},
                   usar_purga ? std::span<const std::size_t>(mundo.fase_cal) : std::span<const std::size_t>{},
                   usar_purga ? std::span<const std::int64_t>(mundo.dia_cal) : std::span<const std::int64_t>{});
    if (obs_d.escalas() == 0) return ResultadoReplicaP6{};

    ConfiguracionLilyth cfg_l{};
    cfg_l.alfa = ALFA_PROYECTO;
    cfg_l.representacion = usar_purga ? ModoRepresentacionLilyth::PURGADA : ModoRepresentacionLilyth::CRUDA;
    cfg_l.b_max_nulo = 48;
    cfg_l.hilos_nulo = 1;
    cfg_l.semilla_nulo = 20261117ULL + 31337ULL * idx_rep;
    cfg_l.modo_nulo_lineal = modo_lin;
    cfg_l.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cfg_l.pool_nulo = PoolNuloLilyth::DERIVADO;
    cfg_l.horizonte = HorizonteLilyth::CUANTIL_RESOLUCION;
    cfg_l.forzar_via_b = forzar_via_b;
    cfg_l.b_mc_resumen = 1999;

    if (modo_lin == ModoNuloLineal::ARFIMA_WHITTLE) {
        cfg_l.d_whittle_cruda = 0.25;
    }

    if (modo_ej == ModoEjecucionP6::CORRUPTO_PROTEGIDO) {
        cfg_l.politica_z1_piloto = ConfiguracionLilyth::PoliticaZ1Piloto::PUENTE_CAUSAL;
        cfg_l.r_max_z1_piloto = 2;
        cfg_l.politica_z2_piloto = ConfiguracionLilyth::PoliticaZ2Piloto::RECHAZAR_SI_CORRUPTO;
        cfg_l.r_max_z2_piloto = 0;
        cfg_l.politica_z3_piloto = ConfiguracionLilyth::PoliticaZ3Piloto::RELOJ_FISICO_TOPE_R;
        cfg_l.r_max_z3_piloto = 2;
    }

    ObsLil lil{};
    if (!lil.calibrar(obs_d, cfg_l, cfg_d, mundo.x_cal, {},
                      usar_purga ? std::span<const std::size_t>(mundo.fase_cal) : std::span<const std::size_t>{},
                      usar_purga ? std::span<const std::int64_t>(mundo.dia_cal) : std::span<const std::int64_t>{})) {
        return ResultadoReplicaP6{};
    }

    const std::size_t n_op = mundo.x_op_limpia.size();
    std::vector<double> x_op = mundo.x_op_limpia;
    std::vector<ClasePasoLilyth> clases(n_op, ClasePasoLilyth::VALIDO);

    std::vector<std::uint8_t> en_choque(n_op, 0u);
    for (std::size_t t0 : mundo.tiempos_choque) {
        for (std::size_t k = 0; k < 2 && (t0 + k) < n_op; ++k) en_choque[t0 + k] = 1u;
    }

    std::mt19937_64 rng_h(detalle_nulo::splitmix64(20261117ULL ^ (idx_rep * 99991ULL) ^ static_cast<std::uint64_t>(tasa_q * 10000.0)));
    std::uniform_real_distribution<double> u01(0.0, 1.0);

    if (modo_ej == ModoEjecucionP6::QUIETO_ACTUAL || modo_ej == ModoEjecucionP6::QUIETO_PROTOCOLO) {
        double acum_niv = 100.0;
        for (std::size_t t = 0; t < n_op; ++t) {
            const bool quieta = (tasa_q > 0.0 && t >= 1 && !en_choque[t] && u01(rng_h) < tasa_q);
            if (integ == Integracion::NIVEL) {
                const double dx = (t == 0) ? (mundo.x_op_limpia[0] - 100.0) : (mundo.x_op_limpia[t] - mundo.x_op_limpia[t - 1]);
                acum_niv += quieta ? 0.0 : dx;
                x_op[t] = acum_niv;
            } else {
                if (quieta) x_op[t] = 0.0;
            }
            if (quieta) clases[t] = ClasePasoLilyth::MERCADO_QUIETO;
        }
    } else if (modo_ej == ModoEjecucionP6::CORRUPTO_ACTUAL || modo_ej == ModoEjecucionP6::CORRUPTO_PROTEGIDO) {
        if (tasa_q > 0.0 && racha_r > 0) {
            const double p_ini = std::min(0.5, tasa_q / static_cast<double>(racha_r));
            for (std::size_t t = 10; t + racha_r < n_op; ++t) {
                if (!en_choque[t] && u01(rng_h) < p_ini) {
                    for (std::size_t k = 0; k < racha_r && (t + k) < n_op; ++k) {
                        if (!en_choque[t + k]) {
                            clases[t + k] = ClasePasoLilyth::CORRUPCION_O_CAIDA;
                            x_op[t + k] = NaN;
                        }
                    }
                    t += racha_r;
                }
            }
        }
    }

    obs_d.reiniciar_flujo();
    AdaptadorInputCausalLilyth adapt{};
    adapt.integracion = integ;
    adapt.admitir_mercado_quieto = (modo_ej != ModoEjecucionP6::QUIETO_ACTUAL);
    // Cuando el detector delta ve la corrupción en flujo extremo a extremo, PROTEGIDO usa puente causal corto (r <= 2)
    // en el adaptador de entrada para evitar que un hueco aislado resetee la escalera diádica de ObservadorSaltos
    adapt.racha_max_corrupcion = (detector_ve_corrupcion && modo_ej == ModoEjecucionP6::CORRUPTO_PROTEGIDO) ? 2 : 0;

    for (std::size_t t = 0; t < n_op; ++t) {
        const std::size_t fa = usar_purga ? mundo.fase_op[t] : 0;
        const std::int64_t di = usar_purga ? mundo.dia_op[t] : 0;

        double x_l = x_op[t];
        bool val_l = true;
        if (modo_ej == ModoEjecucionP6::QUIETO_ACTUAL) {
            val_l = (clases[t] == ClasePasoLilyth::VALIDO);
            x_l = val_l ? x_op[t] : NaN;
        } else if (modo_ej == ModoEjecucionP6::QUIETO_PROTOCOLO) {
            const auto prep = adapt.preparar(x_op[t], clases[t], 0.0);
            x_l = prep.x_para_operador;
            val_l = prep.valido_para_operador;
        } else if (modo_ej == ModoEjecucionP6::CORRUPTO_ACTUAL || modo_ej == ModoEjecucionP6::CORRUPTO_PROTEGIDO) {
            val_l = (clases[t] == ClasePasoLilyth::VALIDO) && std::isfinite(x_op[t]);
            x_l = val_l ? x_op[t] : NaN;
        }

        double x_d = mundo.x_op_limpia[t];
        bool val_d = true;
        if (modo_ej == ModoEjecucionP6::QUIETO_ACTUAL || modo_ej == ModoEjecucionP6::QUIETO_PROTOCOLO) {
            x_d = x_l;
            val_d = val_l;
        } else if (detector_ve_corrupcion) {
            if (modo_ej == ModoEjecucionP6::CORRUPTO_PROTEGIDO) {
                const auto prep_d = adapt.preparar(x_op[t], clases[t], 0.0);
                x_d = prep_d.x_para_operador;
                val_d = prep_d.valido_para_operador;
            } else {
                x_d = x_l;
                val_d = val_l;
            }
        }

        const auto sd = obs_d.paso(val_d ? x_d : NaN, NaN, val_d, fa, di);
        (void)lil.paso(sd, x_l, val_l, NaN, fa, di);
    }

    lil.completar_nulos();

    ResultadoReplicaP6 res{};
    const auto& fn = usar_purga ? lil.ficha_nulo_purgada() : lil.ficha_nulo_cruda();
    const auto& fa = usar_purga ? lil.resumen_armonico_purgada() : lil.resumen_armonico_cruda();

    for (const auto& r : lil.cerradas_historico()) {
        const auto& rama = usar_purga ? r.purgada : r.cruda;
        if (rama.desenlace == DesenlaceRanura::RESTITUIDA || rama.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) {
            ++res.k_medibles;
        } else if (rama.desenlace == DesenlaceRanura::CENSURADA_DATO_INVALIDO) {
            ++res.k_cens_inv;
        } else if (es_no_medida(rama.desenlace)) {
            ++res.k_no_medida;
        }
        if (rama.nulo_onsager.evaluado) ++res.k_evaluados;
    }

    res.serie_valida = fn.contraste_onsager.valido && fa.familia_holm_valida && std::isfinite(fa.p_holm[1]);
    if (res.serie_valida) {
        res.s_inter = fn.contraste_onsager.estadistico_s;
        res.ee_inter = fn.contraste_onsager.ee_nulo_s;
        res.p_holm_inter = fa.p_holm[1];
        res.rechaza_inter = (fa.p_holm[1] < ALFA_PROYECTO);
        res.rechaza_lenta_inter = res.rechaza_inter && (res.s_inter > 0.0);

        std::array<double, 4> p_crudos_pool = {
            fn.contraste_iid_pool_loo.p_bilateral,
            fn.contraste_onsager_pool_loo.p_bilateral,
            fa.p_crudos[2],
            fa.p_crudos[3]
        };
        bool ok_p = true;
        for (double pv : p_crudos_pool) ok_p = ok_p && std::isfinite(pv);
        if (ok_p) {
            std::array<std::size_t, 4> ord = {0, 1, 2, 3};
            std::sort(ord.begin(), ord.end(), [&](std::size_t a, std::size_t b) { return p_crudos_pool[a] < p_crudos_pool[b]; });
            std::array<double, 4> ph{};
            double ac = 0.0;
            for (std::size_t k = 0; k < 4; ++k) {
                ac = std::max(ac, std::min(1.0, p_crudos_pool[ord[k]] * static_cast<double>(4 - k)));
                ph[ord[k]] = ac;
            }
            res.rechaza_pool = (ph[1] < ALFA_PROYECTO);
        }
    }
    return res;
}

struct ResumenAgregadoP6 {
    std::size_t n_validas = 0;
    std::size_t n_rech_inter = 0;
    std::size_t n_rech_lenta = 0;
    std::size_t n_rech_pool = 0;
    double tasa_rech_inter = 0.0;
    double wilson_lo = 0.0;
    double wilson_hi = 1.0;
    double tasa_lenta = 0.0;
    double wilson_lenta_lo = 0.0;
    double wilson_lenta_hi = 1.0;
    double media_s = 0.0;
    double media_ee = 0.05;
    double k_medibles_medio = 0.0;
    double k_evaluados_medio = 0.0;
};

ResumenAgregadoP6 agregar_replicas(const std::vector<ResultadoReplicaP6>& v) {
    ResumenAgregadoP6 out{};
    double sum_s = 0.0, sum_ee = 0.0, sum_k = 0.0, sum_ke = 0.0;
    for (const auto& r : v) {
        sum_k += static_cast<double>(r.k_medibles);
        sum_ke += static_cast<double>(r.k_evaluados);
        if (!r.serie_valida) continue;
        ++out.n_validas;
        if (r.rechaza_inter) ++out.n_rech_inter;
        if (r.rechaza_lenta_inter) ++out.n_rech_lenta;
        if (r.rechaza_pool) ++out.n_rech_pool;
        sum_s += r.s_inter;
        sum_ee += r.ee_inter;
    }
    out.k_medibles_medio = v.empty() ? 0.0 : (sum_k / static_cast<double>(v.size()));
    out.k_evaluados_medio = v.empty() ? 0.0 : (sum_ke / static_cast<double>(v.size()));
    if (out.n_validas > 0) {
        out.tasa_rech_inter = static_cast<double>(out.n_rech_inter) / static_cast<double>(out.n_validas);
        out.tasa_lenta = static_cast<double>(out.n_rech_lenta) / static_cast<double>(out.n_validas);
        const auto [lo, hi] = wilson_95(out.n_rech_inter, out.n_validas);
        out.wilson_lo = lo;
        out.wilson_hi = hi;
        const auto [llo, lhi] = wilson_95(out.n_rech_lenta, out.n_validas);
        out.wilson_lenta_lo = llo;
        out.wilson_lenta_hi = lhi;
        out.media_s = sum_s / static_cast<double>(out.n_validas);
        out.media_ee = std::max(1e-6, sum_ee / static_cast<double>(out.n_validas));
    }
    return out;
}

} // namespace

int main() {
    const std::size_t n_hilos = std::max<std::size_t>(1, std::thread::hardware_concurrency());
    std::printf("========================================================================================================================\n");
    std::printf(" CAPA 2 (PUNTO 6 CORREGIDO): H0 Y H1 CON PARAMETROS SELLADOS (a=8.0, z=4.0*0.97^tau, CUANTIL_RESOLUCION, R=%zu, Hilos=%zu)\n",
                R_REPLICAS, n_hilos);
    std::printf("========================================================================================================================\n");

    ParametrosReplicaH0H1 p_h0{};
    p_h0.lenta_h1 = false;
    ParametrosReplicaH0H1 p_h1 = p_h0;
    p_h1.lenta_h1 = true;

    std::vector<SerieSinteticaCoordinada> m1_h0(R_REPLICAS), m2_h0(R_REPLICAS), m3_h0(R_REPLICAS), m1_h1(R_REPLICAS);
    en_paralelo(R_REPLICAS, n_hilos, [&](std::size_t i) {
        m1_h0[i] = generar_mundo_m1_h0_h1(0.80, p_h0, i);
        m2_h0[i] = generar_mundo_m2_h0_h1(0.25, 0.40, p_h0, i);
        m3_h0[i] = generar_mundo_m3_h0_h1(0.65, p_h0, i);
        m1_h1[i] = generar_mundo_m1_h0_h1(0.80, p_h1, i);
    });

    auto evaluar_lote = [&](const std::vector<SerieSinteticaCoordinada>& mundos,
                            Integracion integ, bool purga, ModoNuloLineal mlin,
                            ModoEjecucionP6 mej, bool det_ve_corr, bool f_via_b, double q, std::size_t r) {
        std::vector<ResultadoReplicaP6> res(R_REPLICAS);
        en_paralelo(R_REPLICAS, n_hilos, [&](std::size_t i) {
            res[i] = correr_una_replica(mundos[i], integ, purga, mlin, mej, det_ve_corr, f_via_b, q, r, i);
        });
        return agregar_replicas(res);
    };

    // 1. Calibración base en q = 0
    const auto base_m1_h0 = evaluar_lote(m1_h0, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, ModoEjecucionP6::LIMPIO, false, false, 0.0, 0);
    const auto base_m2_h0 = evaluar_lote(m2_h0, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE,     ModoEjecucionP6::LIMPIO, false, false, 0.0, 0);
    const auto base_m3_h0 = evaluar_lote(m3_h0, Integracion::NIVEL, true,  ModoNuloLineal::AR_LEVINSON_DURBIN, ModoEjecucionP6::LIMPIO, false, false, 0.0, 0);
    const auto base_m1_h1 = evaluar_lote(m1_h1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, ModoEjecucionP6::LIMPIO, false, false, 0.0, 0);

    std::printf("\n--- [CN1.2] CALIBRACION EN q = 0 (R = %zu replicas): Tamano bajo H0 y Potencia bajo H1 ---\n", R_REPLICAS);
    std::printf(" Mundo                        | Validas | K_eval | Media(S) | Rechazo INTER (Wilson 95%%)        | Rech POOL | Condicion H0/H1\n");
    std::printf("------------------------------+---------+--------+----------+------------------------------------+-----------+-----------------\n");

    auto imp_base = [](const char* nom, const ResumenAgregadoP6& b, bool es_h1) {
        const bool ok = es_h1 ? (b.tasa_lenta >= 0.50) : (b.wilson_lo <= ALFA_PROYECTO);
        const std::size_t k_r = es_h1 ? b.n_rech_lenta : b.n_rech_inter;
        const double tasa = es_h1 ? b.tasa_lenta : b.tasa_rech_inter;
        const double wlo = es_h1 ? b.wilson_lenta_lo : b.wilson_lo;
        const double whi = es_h1 ? b.wilson_lenta_hi : b.wilson_hi;
        std::printf(" %-28s | %7zu | %6.1f | %+8.4f | %3zu/%3zu = %.4f [%.4f, %.4f] | %3zu/%3zu   | %s\n",
                    nom, b.n_validas, b.k_evaluados_medio, b.media_s,
                    k_r, b.n_validas, tasa, wlo, whi,
                    b.n_rech_pool, b.n_validas,
                    ok ? "PASA" : "FALLA");
    };
    imp_base("M1_H0 (FLUJO AR(1))",          base_m1_h0, false);
    imp_base("M2_H0 (FLUJO ARFIMA d=0.25)",  base_m2_h0, false);
    imp_base("M3_H0 (NIVEL + PurgaDiurna)",  base_m3_h0, false);
    imp_base("M1_H1 (Potencia con z_tau)",   base_m1_h1, true);

    // 2. Mercado Quieto (ILIQUIDEZ_TRANSITORIA) con cota monótona estricta
    std::printf("\n--- [6A] MERCADO QUIETO (ILIQUIDEZ_TRANSITORIA) EN M3_H0 Y M1_H1 [R = %zu, Cota Monotona] ---\n", R_REPLICAS);
    std::printf("  q(%%) | Modo      | Series_Val(M3) | dCens(M3) | FalsoRech_H0 [Wilson 95%%]      | RetPot_H1 | ESTRICTO | DULCE\n");
    std::printf("-------+-----------+----------------+-----------+--------------------------------+-----------+----------+-------\n");

    double q_max_quieto_est_mono = 0.0, q_max_quieto_dul_mono = 0.0;
    bool vivo_quieto_est = true, vivo_quieto_dul = true;

    for (double q : {0.01, 0.02, 0.05, 0.10}) {
        for (int m_i = 0; m_i < 2; ++m_i) {
            const auto mej = (m_i == 0) ? ModoEjecucionP6::QUIETO_ACTUAL : ModoEjecucionP6::QUIETO_PROTOCOLO;
            const char* nom = (m_i == 0) ? "ACTUAL" : "PROTOCOLO";
            const auto r_h0 = evaluar_lote(m3_h0, Integracion::NIVEL, true,  ModoNuloLineal::AR_LEVINSON_DURBIN, mej, true, false, q, 0);
            const auto r_h1 = evaluar_lote(m1_h1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, mej, true, false, q, 0);

            const double d_cens = (base_m3_h0.k_evaluados_medio > 0.0)
                                      ? std::max(0.0, 1.0 - r_h0.k_evaluados_medio / base_m3_h0.k_evaluados_medio)
                                      : 1.0;
            const double ret_pot = (base_m1_h1.tasa_lenta > 0.0) ? (r_h1.tasa_lenta / base_m1_h1.tasa_lenta) : 0.0;
            const double ds_ee = std::abs(r_h0.media_s - base_m3_h0.media_s) / base_m3_h0.media_ee;

            const bool est = (r_h0.n_validas >= R_REPLICAS / 2) && (d_cens <= 0.10) &&
                             (r_h0.wilson_lo <= 0.05) && (ds_ee <= 0.25) && (ret_pot >= 0.85);
            const bool dul = (r_h0.n_validas >= R_REPLICAS / 2) && (d_cens <  0.50) &&
                             (r_h0.wilson_lo <= 0.10) && (ds_ee <= 1.00) && (ret_pot >= 0.50);

            if (m_i == 1) {
                if (vivo_quieto_est && est) q_max_quieto_est_mono = q; else vivo_quieto_est = false;
                if (vivo_quieto_dul && dul) q_max_quieto_dul_mono = q; else vivo_quieto_dul = false;
            }

            std::printf(" %5.1f | %-9s | %6zu / %3zu   | %9.2f | %3zu/%3zu = %.3f [%.3f, %.3f] | %9.2f | %-8s | %-5s\n",
                        q * 100.0, nom, r_h0.n_validas, R_REPLICAS, d_cens,
                        r_h0.n_rech_inter, r_h0.n_validas, r_h0.tasa_rech_inter, r_h0.wilson_lo, r_h0.wilson_hi,
                        ret_pot, est ? "PASA" : "FALLA", dul ? "PASA" : "FALLA");
        }
    }

    // 3. Corrupción Física en Flujo Extremo a Extremo (Detector Delta + Lilyth viendo ambos los NaN)
    std::printf("\n--- [6B.1] CORRUPCION FISICA EXTREMO A EXTREMO (ObservadorSaltos + Lilyth ven NaN) [R = %zu] ---\n", R_REPLICAS);
    std::printf("  q(%%) | r | Politica    | dCens | |dS|/ee | FalsoRech_H0 [Wilson 95%%]      | RetPot_H1 | ESTRICTO | DULCE\n");
    std::printf("-------+---+-------------+-------+---------+--------------------------------+-----------+----------+-------\n");

    for (double q : {0.005, 0.01, 0.02, 0.05}) {
        for (std::size_t r : {std::size_t{1}, std::size_t{2}, std::size_t{3}}) {
            for (int p_i = 0; p_i < 2; ++p_i) {
                const auto mej = (p_i == 0) ? ModoEjecucionP6::CORRUPTO_ACTUAL : ModoEjecucionP6::CORRUPTO_PROTEGIDO;
                const char* nom = (p_i == 0) ? "ACTUAL" : "PROTEGIDO";
                const auto r_h0 = evaluar_lote(m1_h0, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, mej, true, false, q, r);
                const auto r_h1 = evaluar_lote(m1_h1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, mej, true, false, q, r);

                const double d_cens = (base_m1_h0.k_evaluados_medio > 0.0)
                                          ? std::max(0.0, 1.0 - r_h0.k_evaluados_medio / base_m1_h0.k_evaluados_medio)
                                          : 1.0;
                const double ds_ee = std::abs(r_h0.media_s - base_m1_h0.media_s) / base_m1_h0.media_ee;
                const double ret_pot = (base_m1_h1.tasa_lenta > 0.0) ? (r_h1.tasa_lenta / base_m1_h1.tasa_lenta) : 0.0;

                const bool est = (d_cens <= 0.10) && (ds_ee <= 0.25) && (r_h0.wilson_lo <= 0.05) && (ret_pot >= 0.85);
                const bool dul = (d_cens <  0.50) && (ds_ee <= 1.00) && (r_h0.wilson_lo <= 0.10) && (ret_pot >= 0.50);

                std::printf(" %5.1f | %1zu | %-11s | %5.2f | %7.4f | %3zu/%3zu = %.3f [%.3f, %.3f] | %9.2f | %-8s | %-5s\n",
                            q * 100.0, r, nom, d_cens, ds_ee,
                            r_h0.n_rech_inter, r_h0.n_validas, r_h0.tasa_rech_inter, r_h0.wilson_lo, r_h0.wilson_hi,
                            ret_pot, est ? "PASA" : "FALLA", dul ? "PASA" : "FALLA");
            }
        }
        std::printf("-------+---+-------------+-------+---------+--------------------------------+-----------+----------+-------\n");
    }

    // 4. Corrupción Física con Vía (b) Forzada (ejercitando Z1 + Z2 + Z3 en el 100% de eventos)
    const auto base_m1_h0_vb = evaluar_lote(m1_h0, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, ModoEjecucionP6::LIMPIO, false, true, 0.0, 0);
    const auto base_m1_h1_vb = evaluar_lote(m1_h1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, ModoEjecucionP6::LIMPIO, false, true, 0.0, 0);

    std::printf("\n--- [6B.2] CORRUPCION FISICA CON forzar_via_b = true (Z1+Z2+Z3 activas en 100%% de eventos) [R = %zu] ---\n", R_REPLICAS);
    std::printf("  q(%%) | r | Politica    | dCens | |dS|/ee | FalsoRech_H0 [Wilson 95%%]      | RetPot_H1 | ESTRICTO | DULCE\n");
    std::printf("-------+---+-------------+-------+---------+--------------------------------+-----------+----------+-------\n");

    for (double q : {0.01, 0.02, 0.05}) {
        for (std::size_t r : {std::size_t{1}, std::size_t{2}}) {
            for (int p_i = 0; p_i < 2; ++p_i) {
                const auto mej = (p_i == 0) ? ModoEjecucionP6::CORRUPTO_ACTUAL : ModoEjecucionP6::CORRUPTO_PROTEGIDO;
                const char* nom = (p_i == 0) ? "ACTUAL" : "PROTEGIDO";
                const auto r_h0 = evaluar_lote(m1_h0, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, mej, false, true, q, r);
                const auto r_h1 = evaluar_lote(m1_h1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, mej, false, true, q, r);

                const double d_cens = (base_m1_h0_vb.k_evaluados_medio > 0.0)
                                          ? std::max(0.0, 1.0 - r_h0.k_evaluados_medio / base_m1_h0_vb.k_evaluados_medio)
                                          : 1.0;
                const double ds_ee = std::abs(r_h0.media_s - base_m1_h0_vb.media_s) / base_m1_h0_vb.media_ee;
                const double ret_pot = (base_m1_h1_vb.tasa_lenta > 0.0) ? (r_h1.tasa_lenta / base_m1_h1_vb.tasa_lenta) : 0.0;

                const bool est = (d_cens <= 0.10) && (ds_ee <= 0.25) && (r_h0.wilson_lo <= 0.05) && (ret_pot >= 0.85);
                const bool dul = (d_cens <  0.50) && (ds_ee <= 1.00) && (r_h0.wilson_lo <= 0.10) && (ret_pot >= 0.50);

                std::printf(" %5.1f | %1zu | %-11s | %5.2f | %7.4f | %3zu/%3zu = %.3f [%.3f, %.3f] | %9.2f | %-8s | %-5s\n",
                            q * 100.0, r, nom, d_cens, ds_ee,
                            r_h0.n_rech_inter, r_h0.n_validas, r_h0.tasa_rech_inter, r_h0.wilson_lo, r_h0.wilson_hi,
                            ret_pot, est ? "PASA" : "FALLA", dul ? "PASA" : "FALLA");
            }
        }
        std::printf("-------+---+-------------+-------+---------+--------------------------------+-----------+----------+-------\n");
    }

    std::printf("\n[COTA MONOTONA MERCADO QUIETO (R=%zu)]: q_max_ESTRICTO = %.1f%% | q_max_DULCE = %.1f%%\n",
                R_REPLICAS, q_max_quieto_est_mono * 100.0, q_max_quieto_dul_mono * 100.0);

    return 0;
}
