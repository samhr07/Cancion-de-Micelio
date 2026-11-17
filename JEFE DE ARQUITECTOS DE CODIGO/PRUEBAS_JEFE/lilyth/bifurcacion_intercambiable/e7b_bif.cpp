// e7b_bif.cpp -- XB.2 y XB.3 de la bifurcacion del contraste del resumen de Lilyth (Code, 2026-11-15). Expectativas selladas en
// EXPECTATIVAS_BIFURCACION.md (sha256 e9cfa2c9...). El acople de E7b (PRUEBAS_JEFE/red/e7b_red.cpp) en MUNDOS NUEVOS (IID 8101+, HETERO 8201+,
// RAFAGAS 8301+, 10 cada uno; CAMBIO 8401+, 5, solo lectura) con Lilyth en INTERCAMBIABLE; la rama vieja sale de los *_pool_loo de la misma
// corrida (A/B), con su Holm de 4 rehecho por detalle_armonico::aplicar_familia_holm_cuatro. La recomposicion se lee con AMBAS referencias:
// por referencia (sin memoria: p_holm[0]; Onsager: p_holm[1]) y por rama, Holm de red sobre las series validas.
// Uso: e7b_bif.exe [hilos] [mundos_reposo] [mundos_rafaga] [mundos_cambio]
#define E7B_SIN_MAIN
#include "e7b_red.cpp"

#include <map>
#include <mutex>

namespace da = ark::lilyth::detalle_armonico;
namespace dn = ark::lilyth::detalle_nulo;

struct Extra {
    // [rama][referencia]: rama 0 = nueva (INTERCAMBIABLE), 1 = vieja (POOL_LOO); referencia 0 = sin memoria, 1 = Onsager
    std::array<std::array<double, 2>, 2> p_holm{}, S{};
    std::array<bool, 2> valida{false, false};
    std::array<bool, 2> resoluble{false, false};
    std::array<std::size_t, 2> pool_min{0, 0}, eventos{0, 0};
    double seg_inter = 0.0;
};
static std::mutex MTX;
static std::map<const Serie*, Extra> EXTRA;

static void gancho(Serie& S, const Lil& lil) {
    Extra e;
    const auto& fn = lil.ficha_nulo_cruda();
    const auto& h = lil.resumen_armonico_cruda();
    e.valida[0] = h.familia_holm_valida;
    e.p_holm[0] = {h.p_holm[0], h.p_holm[1]};
    e.S[0] = {fn.contraste_iid.estadistico_s, fn.contraste_onsager.estadistico_s};
    e.resoluble = {fn.contraste_iid.resoluble, fn.contraste_onsager.resoluble};
    e.pool_min = {fn.contraste_iid.pool_min_estrato, fn.contraste_onsager.pool_min_estrato};
    e.eventos = {fn.contraste_iid.n_eventos, fn.contraste_onsager.n_eventos};
    al::ResumenArmonicoYHolm<1024> v{};
    da::aplicar_familia_holm_cuatro<1024>(fn.contraste_iid_pool_loo.p_bilateral, fn.contraste_onsager_pool_loo.p_bilateral, h.p_crudos[2], h.p_crudos[3],
                                          lil.configuracion().alfa, v);
    e.valida[1] = v.familia_holm_valida;
    e.p_holm[1] = {v.p_holm[0], v.p_holm[1]};
    e.S[1] = {fn.contraste_iid_pool_loo.estadistico_s, fn.contraste_onsager_pool_loo.estadistico_s};
    // XB.3: coste de la rama nueva = rehacer sus dos contrastes del Delta tau sobre los mismos estratos (medido aparte)
    const auto t0 = std::chrono::steady_clock::now();
    const std::size_t tm = lil.configuracion().t_max_efectivo;
    const auto& cer = lil.cerradas_historico();
    for (int sm = 0; sm < 2; ++sm) {
        const auto& emp = lil.nulos_emparejados(false, sm == 1);
        std::vector<std::vector<double>> est;
        for (std::size_t idx = 0; idx < cer.size() && idx < emp.size(); ++idx) {
            const auto& ref = sm == 1 ? cer[idx].cruda.nulo_iid : cer[idx].cruda.nulo_onsager;
            if (!ref.evaluado) continue;
            std::vector<double> x{dn::codificar_tau(cer[idx].cruda.tau_rest_consolidado, tm)};
            for (const auto& ev : emp[idx]) x.push_back(dn::codificar_tau(ev.tau_rest_consolidado, tm));
            est.push_back(std::move(x));
        }
        volatile double sink = dn::evaluar_contraste_intercambiable(est, lil.configuracion().b_mc_resumen, 1, 0.0125).p_bilateral;
        (void)sink;
    }
    e.seg_inter = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::lock_guard<std::mutex> g(MTX);
    EXTRA[&S] = e;
}
static void config_inter(al::ConfiguracionLilyth& c) { c.contraste_resumen = al::ContrasteResumenLilyth::INTERCAMBIABLE; }

struct Lectura {
    std::size_t mundos = 0;
    std::array<std::array<std::size_t, 2>, 2> red{}, red_lenta{}, series{}, series_rech{};
    std::array<std::size_t, 2> no_resolubles{}, con_evento{};
    std::vector<double> S_ons_nueva, pool_ons;
    double seg_inter = 0.0, seg_series = 0.0;
};

static void un_mundo(const char* nm, int clase, std::uint64_t sem, Lectura& L) {
    std::size_t n_raf = 0;
    std::vector<std::array<double, N>> r;
    if (clase == 0) r = retornos(T_BARRAS(), sem, Tipo::IID);
    else if (clase == 1) r = retornos(T_BARRAS(), sem, Tipo::HETERO);
    else if (clase == 2) r = retornos_rafagas(sem, n_raf);
    else r = retornos_cambio(T_BARRAS(), B0() + CADA, sem, Tipo::IID, Tipo::REDUNDANTE);
    const auto C = correr(a_mundo(r), series_por_omision(), HIST_E7B + PRUEBA_E7B);
    const auto M = medir_mundo(C);
    ++L.mundos;
    std::string d;
    const auto sel = series_por_omision();
    std::array<std::array<bool, 2>, 2> rech{}, lenta{};
    for (int b = 0; b < 2; ++b)
        for (int ref = 0; ref < 2; ++ref) {
            std::vector<double> p, s;
            for (const auto& S : M.s) {
                const auto it = EXTRA.find(&S);
                if (it == EXTRA.end() || !S.lista) continue;
                const Extra& e = it->second;
                if (!e.valida[b] || !std::isfinite(e.p_holm[b][ref])) continue;
                p.push_back(e.p_holm[b][ref]);
                s.push_back(e.S[b][ref]);
                ++L.series[b][ref];
                L.series_rech[b][ref] += e.p_holm[b][ref] < me::ALFA_PROYECTO;
            }
            if (p.empty()) continue;
            const auto adj = ark::holm(p);
            for (std::size_t i = 0; i < adj.size(); ++i)
                if (adj[i] < me::ALFA_PROYECTO) { rech[b][ref] = true; lenta[b][ref] = lenta[b][ref] || s[i] > 0.0; }
            L.red[b][ref] += rech[b][ref];
            L.red_lenta[b][ref] += lenta[b][ref];
        }
    for (std::size_t k = 0; k < M.s.size(); ++k) {
        const auto it = EXTRA.find(&M.s[k]);
        if (it == EXTRA.end()) continue;
        const Extra& e = it->second;
        L.seg_inter += e.seg_inter;
        L.seg_series += M.s[k].seg_cal;
        for (int ref = 0; ref < 2; ++ref) if (e.eventos[ref] > 0) { ++L.con_evento[ref]; L.no_resolubles[ref] += !e.resoluble[ref]; }
        if (std::isfinite(e.S[0][1])) { L.S_ons_nueva.push_back(e.S[0][1]); L.pool_ons.push_back(double(e.pool_min[1])); }
        d += nombre(sel[k]) + "{pH ons " + cien(e.p_holm[0][1]) + "/" + cien(e.p_holm[1][1]) + " iid " + cien(e.p_holm[0][0]) + "/" + cien(e.p_holm[1][0]) + " So " + cien(e.S[0][1]) +
             (e.valida[0] ? "" : " NV") + "} ";
    }
    auto si = [&](int b, int ref) { return std::string(rech[b][ref] ? (lenta[b][ref] ? "RECHAZA lenta" : "RECHAZA") : "no"); };
    dato(std::string(nm) + " semilla " + std::to_string(sem) + ": nueva Onsager " + si(0, 1) + ", sin memoria " + si(0, 0) + "; vieja Onsager " + si(1, 1) + ", sin memoria " + si(1, 0) +
         "; " + d);
    std::lock_guard<std::mutex> g(MTX);
    EXTRA.clear();
}

static void leer(const char* nm, const Lectura& L) {
    for (int ref = 1; ref >= 0; --ref)
        lectura((std::string("XB.2 ") + nm + (ref ? " Onsager" : " sin memoria")).c_str(),
                "Holm de red: nueva " + ts(tasa(L.red[0][ref], L.mundos)) + " (lenta " + std::to_string(L.red_lenta[0][ref]) + "), vieja " + ts(tasa(L.red[1][ref], L.mundos)) + " (lenta " +
                    std::to_string(L.red_lenta[1][ref]) + "); series que rechazan su Holm de 4: nueva " + ts(tasa(L.series_rech[0][ref], L.series[0][ref])) + ", vieja " +
                    ts(tasa(L.series_rech[1][ref], L.series[1][ref])) + "; series con eventos no resolubles (p_min > alfa/4) " +
                    ts(tasa(L.no_resolubles[ref], L.con_evento[ref])));
    lectura((std::string("XB.2 ") + nm + " S").c_str(), "S Onsager de la nueva, mediana " + cien(mediana(L.S_ons_nueva)) + "; pool minimo por estrato, mediana " + cien(mediana(L.pool_ons)));
}

int main(int argc, char** argv) {
    HILOS = argc > 1 ? std::stoul(argv[1]) : std::max(1u, std::thread::hardware_concurrency());
    const std::size_t mr = argc > 2 ? std::stoul(argv[2]) : 10, mf = argc > 3 ? std::stoul(argv[3]) : 10, mc = argc > 4 ? std::stoul(argv[4]) : 5;
    GANCHO_CONFIG_LILYTH = config_inter;
    GANCHO_POSTERIOR = gancho;
    std::printf("Bifurcacion del contraste del resumen de Lilyth, XB.2 sobre E7b en mundos nuevos (reposo %zu + %zu, rafagas %zu, cambio %zu; hilos %zu)\n", mr, mr, mf, mc, HILOS);
    const auto t0 = std::chrono::steady_clock::now();
    Lectura I, H, R, Cb;
    for (std::size_t q = 0; q < mr; ++q) un_mundo("IID", 0, 8101 + q, I);
    leer("IID", I);
    for (std::size_t q = 0; q < mr; ++q) un_mundo("HETERO", 1, 8201 + q, H);
    leer("HETERO", H);
    for (std::size_t q = 0; q < mf; ++q) un_mundo("RAFAGAS", 2, 8301 + q, R);
    leer("RAFAGAS", R);
    for (std::size_t q = 0; q < mc; ++q) un_mundo("CAMBIO", 3, 8401 + q, Cb);
    if (mc) leer("CAMBIO", Cb);

    const double alfa = me::ALFA_PROYECTO;
    const auto rep_nueva = tasa(I.red[0][1] + H.red[0][1], I.mundos + H.mundos);
    reporta(!(rep_nueva.lo > alfa), "XB.2 reposo Onsager", "rama nueva, Holm de red contra Onsager en IID + HETERO: " + ts(rep_nueva) + " (Wilson inferior <= alfa)");
    reporta(R.red[0][1] >= std::min<std::size_t>(8, R.mundos) && R.red_lenta[0][1] == R.red[0][1], "XB.2 rafagas Onsager",
            "rama nueva: " + std::to_string(R.red[0][1]) + " de " + std::to_string(R.mundos) + " mundos rechazan, " + std::to_string(R.red_lenta[0][1]) + " por lenta (>= 8, todos por lenta)");
    const auto rep_iid = I.red[0][0] + H.red[0][0];
    std::printf("[%s] XB.2e sin memoria en reposo: rama nueva rechaza en %zu de %zu (se esperaba >= 5: la memoria del solape)\n", rep_iid >= 5 ? "SE CUMPLE" : "NO SE CUMPLE", rep_iid,
                I.mundos + H.mundos);
    const auto rep_vieja = I.red[1][1] + H.red[1][1];
    std::printf("[%s] XB.2e' rama vieja contra Onsager en reposo: %zu de %zu, la nueva %zu (se esperaba vieja >= nueva)\n", rep_vieja >= I.red[0][1] + H.red[0][1] ? "SE CUMPLE" : "NO SE CUMPLE",
                rep_vieja, I.mundos + H.mundos, std::size_t(I.red[0][1] + H.red[0][1]));
    const double frac = (I.seg_inter + H.seg_inter + R.seg_inter) / std::max(1e-9, I.seg_series + H.seg_series + R.seg_series);
    lectura("XB.3", "coste de los dos contrastes intercambiables del Delta tau (rehechos aparte) sobre el tiempo de la serie: " + cien(100.0 * frac) + " %");
    std::printf("[DATO] %.0f s\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    std::printf("%s: %d fallas\n", fallas ? "FALLA" : "OK", fallas);
    return fallas ? 1 : 0;
}
