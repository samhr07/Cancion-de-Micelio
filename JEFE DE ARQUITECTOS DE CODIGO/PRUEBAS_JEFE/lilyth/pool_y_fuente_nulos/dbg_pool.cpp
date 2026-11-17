// dbg_pool.cpp -- diagnostico (Code): una replica del mundo lineal de EC6 con el detalle por evento.
#define CONTROL_SIN_MAIN
#include "control_fuente_pool.cpp"
int main(int argc, char** argv) {
    const bool lenta = argc > 1 && std::string(argv[1]) == "lenta";
    for (std::uint64_t i = 0; i < 6; ++i) {
        const std::uint64_t sem = 9100000ULL + i;
        std::mt19937_64 g(sem);
        const auto x_cal = ar1(3000, 0.82, g);
        const auto cd = cfg_delta(2000.0, sem ^ 0x1234ULL);
        ObsReal delta; delta.calibrar(cd, x_cal);
        Lil lil;
        if (!lil.calibrar(delta, cfg_lil(al::PoolNuloLilyth::DERIVADO, sem ^ 0xBEEFULL), cd, x_cal)) { std::printf("no calibra\n"); continue; }
        auto x = ar1(400, 0.82, g);
        propagar(x, 200, 12.0, 0.82);
        if (lenta) propagar(x, 200, 4.0, 0.97);
        delta.reiniciar_flujo();
        for (std::size_t t = 0; t < x.size(); ++t) (void)lil.paso(delta.paso(x[t]), x[t], true);
        lil.completar_nulos();
        const auto& fn = lil.ficha_nulo_cruda();
        const auto& c = fn.contraste_onsager;
        const auto& h = lil.resumen_armonico_cruda();
        std::printf("rep %llu Tmax %zu K %zu obj %zu | ons S %.4f p %.4f pmin %.5f res %d pool_min %zu | ph %.4f %.4f %.4f %.4f val %d\n", (unsigned long long)i,
                    lil.configuracion().t_max_efectivo, fn.k_candidatos_pool, fn.objetivo_pool, c.estadistico_s, c.p_bilateral, c.p_min_alcanzable, int(c.resoluble), c.pool_min_estrato,
                    h.p_holm[0], h.p_holm[1], h.p_holm[2], h.p_holm[3], int(h.familia_holm_valida));
        const auto& emp = lil.nulos_emparejados(false, false);
        for (std::size_t k = 0; k < lil.cerradas_historico().size(); ++k) {
            const auto& rn = lil.cerradas_historico()[k];
            const auto& r = rn.cruda.nulo_onsager;
            std::vector<double> t;
            if (k < emp.size()) for (const auto& e : emp[k]) t.push_back(al::detalle_nulo::codificar_tau(e.tau_rest_consolidado, lil.configuracion().t_max_efectivo));
            std::sort(t.begin(), t.end());
            std::size_t mx = 0; for (double v : t) mx += (v == t.back());
            std::printf("   ev t_a %lld %s tau %.1f | eval %d n %zu (a %zu, b %zu, intentos %zu) u %.3f med %.1f max %.1f #max %zu | %s\n", (long long)rn.t_apertura,
                        al::NOMBRE_DESENLACE[int(rn.cruda.desenlace)], rn.cruda.tau_rest_consolidado, int(r.evaluado), r.n_nulos_usados, r.n_nulos_a, r.n_nulos_b, r.intentos_b, r.u_pit,
                        r.mediana_nulo, t.empty() ? -1.0 : t.back(), mx, r.motivo);
        }
    }
}
