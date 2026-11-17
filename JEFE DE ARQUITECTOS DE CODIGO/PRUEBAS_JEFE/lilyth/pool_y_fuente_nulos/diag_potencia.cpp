// diag_potencia.cpp -- DIAGNOSTICO POSTERIOR declarado (Code, 2026-11-16): por que la potencia de EC6c (un evento lento, DERIVADO) es baja.
// Mide, en las 400 replicas de potencia: T_max, si el real sale censurado, fraccion de nulos censurados en su estrato (empates en T_max + 1),
// y cuantas series tienen p_min <= alfa/4 contra cuantas rechazan.
#define CONTROL_SIN_MAIN
#include "control_fuente_pool.cpp"
int main() {
    std::vector<double> tm, frac_cens_null;
    std::size_t reales = 0, cens = 0, cens_y_empata = 0, resol = 0, resol_cens = 0, rech = 0, u_alto = 0;
    std::mutex m;
    en_paralelo(400, 12, [&](std::size_t i) {
        const std::uint64_t sem = 9100000ULL + i;
        std::mt19937_64 g(sem);
        const auto x_cal = ar1(3000, 0.82, g);
        const auto cd = cfg_delta(2000.0, sem ^ 0x1234ULL);
        ObsReal delta; delta.calibrar(cd, x_cal);
        Lil lil;
        if (!lil.calibrar(delta, cfg_lil(al::PoolNuloLilyth::DERIVADO, sem ^ 0xBEEFULL), cd, x_cal)) return;
        auto x = ar1(400, 0.82, g);
        propagar(x, 200, 12.0, 0.82);
        propagar(x, 200, 4.0, 0.97);
        delta.reiniciar_flujo();
        for (std::size_t t = 0; t < x.size(); ++t) (void)lil.paso(delta.paso(x[t]), x[t], true);
        lil.completar_nulos();
        const std::size_t T = lil.configuracion().t_max_efectivo;
        const auto& emp = lil.nulos_emparejados(false, false);
        const auto& cer = lil.cerradas_historico();
        const auto& c = lil.ficha_nulo_cruda().contraste_onsager;
        const auto L = leer(lil);
        std::lock_guard<std::mutex> gg(m);
        tm.push_back(double(T));
        resol += c.resoluble;
        rech += L.lenta[1];
        for (std::size_t k = 0; k < cer.size() && k < emp.size(); ++k) {
            const auto& rn = cer[k];
            if (!rn.cruda.nulo_onsager.evaluado || rn.t_apertura < 199 || rn.t_apertura > 205) continue;
            ++reales;
            const bool rc = rn.cruda.desenlace == al::DesenlaceRanura::CENSURADA_HORIZONTE;
            cens += rc;
            u_alto += rn.cruda.nulo_onsager.u_pit > 0.9;
            std::size_t nc = 0;
            for (const auto& e : emp[k]) nc += !std::isfinite(e.tau_rest_consolidado) || e.tau_rest_consolidado > double(T);
            frac_cens_null.push_back(emp[k].empty() ? 0.0 : double(nc) / emp[k].size());
            cens_y_empata += rc && nc > 0;
            resol_cens += rc && c.resoluble;
        }
    });
    std::sort(tm.begin(), tm.end());
    std::sort(frac_cens_null.begin(), frac_cens_null.end());
    auto q = [](const std::vector<double>& v, double p) { return v.empty() ? NAN : v[std::min(v.size() - 1, std::size_t(p * v.size()))]; };
    std::printf("series %zu; T_max p10/p50/p90 = %.0f/%.0f/%.0f\n", tm.size(), q(tm, .1), q(tm, .5), q(tm, .9));
    std::printf("tormentas lentas evaluadas %zu: censuradas en T_max %zu (%.3f); u_pit > 0.9: %zu; censuradas que empatan con nulos censurados %zu\n",
                reales, cens, reales ? double(cens) / reales : 0.0, u_alto, cens_y_empata);
    std::printf("fraccion de nulos censurados en el estrato del real p10/p50/p90 = %.4f/%.4f/%.4f (alfa/8 = %.4f)\n", q(frac_cens_null, .1), q(frac_cens_null, .5),
                q(frac_cens_null, .9), ALFA / 8.0);
    std::printf("series resolubles %zu (de ellas con el real censurado %zu); rechazan por lenta %zu\n", resol, resol_cens, rech);
}
