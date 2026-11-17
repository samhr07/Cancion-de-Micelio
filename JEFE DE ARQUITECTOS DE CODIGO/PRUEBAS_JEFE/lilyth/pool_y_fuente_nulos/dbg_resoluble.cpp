// dbg_resoluble.cpp -- diagnostico (Code): series con resoluble verdadero y prod (n_k + 1) < 160.
#define CONTROL_SIN_MAIN
#include "control_fuente_pool.cpp"
int main() {
    int mostrados = 0;
    for (std::uint64_t i = 0; i < 400 && mostrados < 3; ++i) {
        const std::uint64_t sem = 9100000ULL + i;
        std::mt19937_64 g(sem);
        const auto x_cal = ar1(3000, 0.82, g);
        const auto cd = cfg_delta(2000.0, sem ^ 0x1234ULL);
        ObsReal delta; delta.calibrar(cd, x_cal);
        Lil lil;
        if (!lil.calibrar(delta, cfg_lil(al::PoolNuloLilyth::DERIVADO, sem ^ 0xBEEFULL), cd, x_cal)) continue;
        auto x = ar1(400, 0.82, g);
        propagar(x, 200, 12.0, 0.82);
        delta.reiniciar_flujo();
        for (std::size_t t = 0; t < x.size(); ++t) (void)lil.paso(delta.paso(x[t]), x[t], true);
        lil.completar_nulos();
        const auto& fn = lil.ficha_nulo_cruda();
        const auto& c = fn.contraste_onsager;
        double prod = 1.0; std::size_t ne = 0;
        for (const auto& rn : lil.cerradas_historico()) if (rn.cruda.nulo_onsager.evaluado) { prod *= rn.cruda.nulo_onsager.n_nulos_usados + 1.0; ++ne; }
        if (!(c.resoluble && prod < 160.0)) continue;
        ++mostrados;
        std::printf("rep %llu: valido %d K %zu (eventos con onsager evaluado %zu) n_nulos_pool %zu pool_min %zu S %.4f p %.4f pmin %.6f res %d prod %.0f inter %d\n",
                    (unsigned long long)i, int(c.valido), c.n_eventos, ne, c.n_nulos_pool, c.pool_min_estrato, c.estadistico_s, c.p_bilateral, c.p_min_alcanzable,
                    int(c.resoluble), prod, int(c.intercambiable));
        for (const auto& rn : lil.cerradas_historico())
            std::printf("   t_a %lld %s eval %d n %zu u %.3f | %s\n", (long long)rn.t_apertura, al::NOMBRE_DESENLACE[int(rn.cruda.desenlace)], int(rn.cruda.nulo_onsager.evaluado),
                        rn.cruda.nulo_onsager.n_nulos_usados, rn.cruda.nulo_onsager.u_pit, rn.cruda.nulo_onsager.motivo);
    }
}
