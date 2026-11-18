// ============================================================================
// volcado_identidad_nivel1.cpp
// Capa 1 (CN1.1 real): Volcado hexadecimal (%a) de ObservadorLilyth en M1, M2 y M3
// (tanto en flujo limpio q=0 como con huecos de mercado quieto y corrupción),
// para comparar byte a byte la cabecera original intacta (observador_lilyth_original.hpp)
// contra la cabecera piloto instrumentada (observador_lilyth.hpp) con valores por omisión.
// ============================================================================

#ifdef USAR_CABECERA_ORIGINAL
#include <ark/lilyth/observador_lilyth_original.hpp>
#else
#include <ark/lilyth/observador_lilyth.hpp>
#endif

#include "banco_sintetico_comun.hpp"
#include <ark/metrologia/whittle_exacto.hpp>

#include <cstdio>
#include <string>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;
using namespace ark::lilyth::piloto;

constexpr std::size_t T_MAX_V = 64;
using ObsDelta = ObservadorSaltos<24, 64, 120>;
using ObsLil = ObservadorLilyth<ObsDelta, 8, T_MAX_V>;

static std::string hx(double x) {
    char b[64];
    std::snprintf(b, sizeof b, "%a", x);
    return b;
}

static void volcar_contraste(const char* tag, const ContrasteResumenNulo& c) {
    std::printf("%s val=%d K=%zu M=%zu B=%zu S=%s ee=%s p=%s md=%s mu=%s pmin=%s res=%d\n",
                tag, int(c.valido), c.n_eventos, c.n_nulos_pool, c.b_mc,
                hx(c.estadistico_s).c_str(), hx(c.ee_nulo_s).c_str(),
                hx(c.p_bilateral).c_str(), hx(c.mediana_delta_tau).c_str(),
                hx(c.media_u_pit).c_str(), hx(c.p_min_alcanzable).c_str(), int(c.resoluble));
}

static void volcar_ref(const char* tag, const ReferenciaNuloEvento& r) {
    std::printf("%s ev=%d iny=%d na=%zu nb=%zu nu=%zu int=%zu pdet=%s u=%s dt=%s md=%s;",
                tag, int(r.evaluado), int(r.nulo_por_inyeccion), r.n_nulos_a, r.n_nulos_b,
                r.n_nulos_usados, r.intentos_b, hx(r.p_det_inicial).c_str(),
                hx(r.u_pit).c_str(), hx(r.delta_tau).c_str(), hx(r.mediana_nulo).c_str());
}

static void correr_y_volcar(
    const char* nombre_caso,
    const SerieSinteticaCoordinada& m,
    Integracion integ,
    bool purga,
    ModoNuloLineal mlin,
    bool inyectar_huecos) {

    ConfiguracionObservador cd{};
    cd.integracion = integ;
    cd.agregacion = AgregacionEscala::SUMA;
    cd.arl0_pasos_nativos = 450.0;
    cd.alfa = ALFA_PROYECTO;
    cd.purga = purga;
    cd.bins_fase = purga ? 24 : 1;
    cd.dias_minimos_purga = 20;
    cd.medir_tau50 = false;
    cd.semilla = 20261117ULL;
    cd.escalas_max = 3;

    ObsDelta od{};
    od.calibrar(cd, m.x_cal, {}, {},
                purga ? std::span<const std::size_t>(m.fase_cal) : std::span<const std::size_t>{},
                purga ? std::span<const std::int64_t>(m.dia_cal) : std::span<const std::int64_t>{});

    ConfiguracionLilyth cl{};
    cl.alfa = ALFA_PROYECTO;
    cl.representacion = purga ? ModoRepresentacionLilyth::PURGADA : ModoRepresentacionLilyth::CRUDA;
    cl.b_max_nulo = 16;
    cl.hilos_nulo = 1; // 1 hilo para determinismo exacto bit a bit
    cl.semilla_nulo = 20261117ULL;
    cl.modo_nulo_lineal = mlin;
    cl.contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;
    cl.b_mc_resumen = 1999;
    if (mlin == ModoNuloLineal::ARFIMA_WHITTLE) {
        cl.d_whittle_cruda = 0.25;
    }

    ObsLil lil{};
    lil.calibrar(od, cl, cd, m.x_cal, {},
                 purga ? std::span<const std::size_t>(m.fase_cal) : std::span<const std::size_t>{},
                 purga ? std::span<const std::int64_t>(m.dia_cal) : std::span<const std::int64_t>{});

    od.reiniciar_flujo();
    for (std::size_t t = 0; t < m.x_op_limpia.size(); ++t) {
        const bool ok = !(inyectar_huecos && (t % 97 == 50 || t == 804));
        const double x_in = ok ? m.x_op_limpia[t] : NaN;
        const std::size_t fa = purga ? m.fase_op[t] : 0;
        const std::int64_t di = purga ? m.dia_op[t] : 0;
        const auto sd = od.paso(x_in, NaN, ok, fa, di);
        (void)lil.paso(sd, x_in, ok, NaN, fa, di);
    }

    lil.completar_nulos();

    const auto& fn = purga ? lil.ficha_nulo_purgada() : lil.ficha_nulo_cruda();
    const auto& fa = purga ? lil.resumen_armonico_purgada() : lil.resumen_armonico_cruda();

    std::printf("CASO %s tmax=%zu ap=%zu cer=%zu\n",
                nombre_caso, lil.configuracion().t_max_efectivo,
                lil.contadores().aperturas, lil.cerradas_historico().size());
    volcar_contraste("C_iid", fn.contraste_iid);
    volcar_contraste("C_ons", fn.contraste_onsager);
    volcar_contraste("C_iid_loo", fn.contraste_iid_pool_loo);
    volcar_contraste("C_ons_loo", fn.contraste_onsager_pool_loo);
    std::printf("HOLM val=%d pc=%s %s %s %s ph=%s %s %s %s\n",
                int(fa.familia_holm_valida),
                hx(fa.p_crudos[0]).c_str(), hx(fa.p_crudos[1]).c_str(),
                hx(fa.p_crudos[2]).c_str(), hx(fa.p_crudos[3]).c_str(),
                hx(fa.p_holm[0]).c_str(), hx(fa.p_holm[1]).c_str(),
                hx(fa.p_holm[2]).c_str(), hx(fa.p_holm[3]).c_str());

    for (const auto& r : lil.cerradas_historico()) {
        const auto& rama = purga ? r.purgada : r.cruda;
        std::printf("RAN ta=%lld tc=%lld j=%d s=%d des=%d tr=%s min2=%s np=%zu nh=%zu ",
                    static_cast<long long>(r.t_apertura), static_cast<long long>(r.t_cambio),
                    r.escala_origen, r.signo, static_cast<int>(rama.desenlace),
                    hx(rama.tau_rest_consolidado).c_str(), hx(rama.min_eta_desde_2).c_str(),
                    purga ? r.n_perfil_purgada : r.n_perfil_cruda,
                    purga ? r.n_hist_prev_purgada : r.n_hist_prev_cruda);
        volcar_ref("iid", rama.nulo_iid);
        volcar_ref("ons", rama.nulo_onsager);
        std::printf("\n");
    }
}

} // namespace

int main() {
    ParametrosBancoComun p{};
    p.n_op = 4000;
    const auto m1 = generar_mundo_m1_ar1_flujo(0.80, p, 101);
    const auto m2 = generar_mundo_m2_arfima_flujo(0.25, 0.40, p, 202);
    const auto m3 = generar_mundo_m3_nivel_purga(0.65, p, 303);

    correr_y_volcar("M1_LIMPIO", m1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, false);
    correr_y_volcar("M1_HUECOS", m1, Integracion::FLUJO, false, ModoNuloLineal::AR_LEVINSON_DURBIN, true);
    correr_y_volcar("M2_LIMPIO", m2, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE,     false);
    correr_y_volcar("M2_HUECOS", m2, Integracion::FLUJO, false, ModoNuloLineal::ARFIMA_WHITTLE,     true);
    correr_y_volcar("M3_LIMPIO", m3, Integracion::NIVEL, true,  ModoNuloLineal::AR_LEVINSON_DURBIN, false);
    correr_y_volcar("M3_HUECOS", m3, Integracion::NIVEL, true,  ModoNuloLineal::AR_LEVINSON_DURBIN, true);
    return 0;
}
