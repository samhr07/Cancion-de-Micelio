// volcado_ab.cpp -- XB.0 (ii), XB.0b y XB.0c de la bifurcacion del contraste del resumen de Lilyth (Code, 2026-11-15).
// Corre el acople de E7b (PRUEBAS_JEFE/red/e7b_red.cpp) en 3 mundos (IID 7101, RAFAGAS 7301, CAMBIO 7401) y vuelca en hexadecimal, por serie,
// todo lo de Lilyth que la bifurcacion pudo tocar. Modos:
//   previa : compilado contra base_previa/ (la base de antes)                         -> volcado A
//   nueva  : compilado contra la base nueva con /DBASE_NUEVA, POOL_LOO por omision     -> volcado A   (XB.0 ii: identico a previa)
//   inter  : base nueva, INTERCAMBIABLE; vuelca los *_pool_loo en el lugar de los contrastes y aparta las lineas de Holm (dependen de la
//            rama) con el prefijo "H"                                                -> identico a nueva salvo lineas H (XB.0b), mas XB.0c
// Uso: volcado_ab.exe hilos modo > archivo
#define E7B_SIN_MAIN
#include "e7b_red.cpp"

#include <mutex>

static std::string MODO;
static std::mutex MTX;
static std::vector<std::string> BLOQUES;
static std::atomic<std::size_t> U_MAL{0}, N_MAL{0}, P_MAL{0}, ESTRATOS{0};

static std::string hx(double x) { char b[64]; std::snprintf(b, sizeof b, "%a", x); return b; }
static std::string contraste(const char* nm, const al::ContrasteResumenNulo& c) {
    return std::string(nm) + " v" + std::to_string(int(c.valido)) + " K" + std::to_string(c.n_eventos) + " M" + std::to_string(c.n_nulos_pool) + " B" + std::to_string(c.b_mc) + " S" +
           hx(c.estadistico_s) + " ee" + hx(c.ee_nulo_s) + " p" + hx(c.p_bilateral) + " md" + hx(c.mediana_delta_tau) + " mu" + hx(c.media_u_pit) + "\n";
}
static std::string referencia(const char* nm, const al::ReferenciaNuloEvento& r) {
    return std::string(nm) + " e" + std::to_string(int(r.evaluado)) + " n" + std::to_string(r.n_nulos_usados) + " u" + hx(r.u_pit) + " dt" + hx(r.delta_tau) + " md" + hx(r.mediana_nulo) +
           " mx" + hx(r.max_nulo) + " c" + std::to_string(int(r.real_censurado)) + ";";
}

static void volcar(Serie&, const Lil& lil) {
    const auto& fn = lil.ficha_nulo_cruda();
    const auto& h = lil.resumen_armonico_cruda();
    const auto& sb = h.sobrepaso;
    const bool inter = MODO == "inter";
    std::string s = "SERIE tmax " + std::to_string(lil.configuracion().t_max_efectivo) + " aperturas " + std::to_string(lil.contadores().aperturas) + " cerradas " +
                    std::to_string(lil.cerradas_historico().size()) + "\n";
#ifdef BASE_NUEVA
    s += contraste("C iid", inter ? fn.contraste_iid_pool_loo : fn.contraste_iid);
    s += contraste("C ons", inter ? fn.contraste_onsager_pool_loo : fn.contraste_onsager);
    s += contraste("SE iid", inter ? sb.sobrepaso_por_evento_iid_pool_loo : sb.sobrepaso_por_evento_iid);
    s += contraste("SE ons", inter ? sb.sobrepaso_por_evento_onsager_pool_loo : sb.sobrepaso_por_evento_onsager);
#else
    s += contraste("C iid", fn.contraste_iid);
    s += contraste("C ons", fn.contraste_onsager);
    s += contraste("SE iid", sb.sobrepaso_por_evento_iid);
    s += contraste("SE ons", sb.sobrepaso_por_evento_onsager);
#endif
    s += "SC p" + hx(sb.p_sobrepaso_curva_iid) + " " + hx(sb.p_sobrepaso_curva_onsager) + " m" + hx(sb.mediana_min_eta_barra_iid) + " " + hx(sb.mediana_min_eta_barra_onsager) + " r" +
         hx(sb.min_eta_barra_real) + "\n";
    s += "P3 z" + hx(fn.p3b.max_z_desvio) + " d" + std::to_string(int(fn.p3b.dentro_de_banda_2ee)) + " B" + std::to_string(fn.B) + " bp" + hx(fn.b_pedido) + " T" +
         std::to_string(fn.t_max_derivado) + "\n";
    s += "H pc";
    for (double p : h.p_crudos) s += " " + hx(p);
    s += " ph";
    for (double p : h.p_holm) s += " " + hx(p);
    s += " val" + std::to_string(int(h.familia_holm_valida)) + "\n";
    for (const auto& rn : lil.cerradas_historico()) {
        const auto& c = rn.cruda;
        s += "R " + std::to_string(rn.t_apertura) + " j" + std::to_string(rn.escala_origen) + " s" + std::to_string(rn.signo) + " " + al::NOMBRE_DESENLACE[static_cast<int>(c.desenlace)] + " tr" +
             hx(c.tau_rest_consolidado) + " " + referencia("iid", c.nulo_iid) + " " + referencia("ons", c.nulo_onsager) + "\n";
    }
#ifdef BASE_NUEVA
    if (inter) {   // XB.0c
        namespace dn = ark::lilyth::detalle_nulo;
        const std::size_t tm = lil.configuracion().t_max_efectivo;
        const auto& cer = lil.cerradas_historico();
        for (int sm = 0; sm < 2; ++sm) {
            const auto& emp = lil.nulos_emparejados(false, sm == 1);
            for (std::size_t idx = 0; idx < cer.size() && idx < emp.size(); ++idx) {
                const auto& ref = sm == 1 ? cer[idx].cruda.nulo_iid : cer[idx].cruda.nulo_onsager;
                if (!ref.evaluado) continue;
                ++ESTRATOS;
                N_MAL += emp[idx].size() != ref.n_nulos_usados;
                std::vector<double> resto;
                for (const auto& e : emp[idx]) resto.push_back(dn::codificar_tau(e.tau_rest_consolidado, tm));
                U_MAL += std::abs(dn::pit_punto_medio(dn::codificar_tau(cer[idx].cruda.tau_rest_consolidado, tm), resto) - ref.u_pit) > 1e-12;
            }
        }
        // los p 1 y 2 del Holm de 4 son los p_bilateral de los contrastes nuevos, y los contrastes vienen de la rama nueva
        P_MAL += !igual_bits(h.p_crudos[0], fn.contraste_iid.p_bilateral) || !igual_bits(h.p_crudos[1], fn.contraste_onsager.p_bilateral) ||
                 !fn.contraste_iid.intercambiable || !fn.contraste_onsager.intercambiable;
        s += "HI nuevo iid S" + hx(fn.contraste_iid.estadistico_s) + " p" + hx(fn.contraste_iid.p_bilateral) + " pmin" + hx(fn.contraste_iid.p_min_alcanzable) + " res" +
             std::to_string(int(fn.contraste_iid.resoluble)) + " ons S" + hx(fn.contraste_onsager.estadistico_s) + " p" + hx(fn.contraste_onsager.p_bilateral) + " pmin" +
             hx(fn.contraste_onsager.p_min_alcanzable) + " res" + std::to_string(int(fn.contraste_onsager.resoluble)) + "\n";
    }
#endif
    std::lock_guard<std::mutex> g(MTX);
    BLOQUES.push_back(std::move(s));
}

#ifdef BASE_NUEVA
static void config_inter(al::ConfiguracionLilyth& c) { c.contraste_resumen = al::ContrasteResumenLilyth::INTERCAMBIABLE; }
#endif

int main(int argc, char** argv) {
    HILOS = argc > 1 ? std::stoul(argv[1]) : 8;
    MODO = argc > 2 ? argv[2] : "nueva";
#ifdef BASE_NUEVA
    if (MODO == "previa") { std::fprintf(stderr, "el modo previa se compila sin BASE_NUEVA\n"); return 2; }
    if (MODO == "inter") GANCHO_CONFIG_LILYTH = config_inter;
#else
    if (MODO != "previa") { std::fprintf(stderr, "sin BASE_NUEVA solo el modo previa\n"); return 2; }
#endif
    GANCHO_POSTERIOR = volcar;
    struct Mundo { const char* nm; int clase; std::uint64_t sem; };
    const std::string solo = argc > 3 ? argv[3] : "";   // opcional: mundos separados por coma (los mutantes corren solo IID)
    for (const Mundo& w : {Mundo{"IID", 0, 7101}, Mundo{"RAFAGAS", 2, 7301}, Mundo{"CAMBIO", 3, 7401}}) {
        if (!solo.empty() && ("," + solo + ",").find("," + std::string(w.nm) + ",") == std::string::npos) continue;
        std::size_t n_raf = 0;
        std::vector<std::array<double, N>> r;
        if (w.clase == 0) r = retornos(T_BARRAS(), w.sem, Tipo::IID);
        else if (w.clase == 2) r = retornos_rafagas(w.sem, n_raf);
        else r = retornos_cambio(T_BARRAS(), B0() + CADA, w.sem, Tipo::IID, Tipo::REDUNDANTE);
        const auto C = correr(a_mundo(r), series_por_omision(), HIST_E7B + PRUEBA_E7B);
        BLOQUES.clear();
        (void)medir_mundo(C);
        std::sort(BLOQUES.begin(), BLOQUES.end());
        std::printf("MUNDO %s\n", w.nm);
        for (const auto& b : BLOQUES) std::fputs(b.c_str(), stdout);
        std::fflush(stdout);
    }
    if (MODO == "inter")
        std::fprintf(stderr, "XB.0c: estratos %zu; u del real distinto del u_pit %zu; estrato distinto de n_nulos_usados %zu; series con p del Holm de 4 distinto del contraste nuevo %zu\n",
                     ESTRATOS.load(), U_MAL.load(), N_MAL.load(), P_MAL.load());
    return 0;
}
