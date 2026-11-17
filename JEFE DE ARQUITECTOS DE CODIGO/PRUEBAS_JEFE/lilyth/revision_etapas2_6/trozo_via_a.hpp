// Genera la parte ESTACIONARIA (FLUJO: y; NIVEL: incrementos) de longitud n de una realizacion nula, en el espacio de la
// representacion: sin memoria por sustituto de los segmentos de calma, o lineal por ModeloLinealNulo::simular. Code (revision): la
// reintegracion y el paso al crudo los hace `llevar_a_crudo` (antes se reintegraba aqui sin purga ni fases).
inline void construir_estacionaria_nula(
    bool es_lineal,
    std::uint32_t codigo_sustituto,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    std::size_t n,
    std::mt19937_64& rng,
    std::vector<double>& est) {
    est.assign(n, 0.0);
    if (n == 0 || segs_calma.empty()) return;
    if (es_lineal && m_lin.valido) {
        std::vector<double> buf(n + m_lin.burn_in + m_lin.p + 16, 0.0);
        simular(m_lin, rng, nullptr, 0, 0, est.data(), n, buf.data(), buf.size());
        return;
    }
    // Sin memoria: un sustituto (PERMUTACION o SIGNOS) por vuelta sobre los segmentos; si la realizacion es mas larga que la historia
    // se saca un sustituto NUEVO en cada vuelta (antes se repetia el mismo).
    const auto mo = ark::metrologia::est::momentos(segs_calma);
    std::size_t pos = 0;
    while (pos < n) {
        const auto segs_sust = ark::metrologia::sustituto(segs_calma, codigo_sustituto, 1, mo.mu, rng);
        const std::size_t antes = pos;
        for (const auto& g : segs_sust) {
            for (double v : g) {
                if (pos >= n) break;
                est[pos++] = v;
            }
            if (pos >= n) break;
        }
        if (pos == antes) break;
    }
}

// Corre B realizaciones espontaneas de la via (a) en paralelo (seccion 3.2 y reglas 7 y N2/N3). Code (revision D1): sin relleno
// fabricado; la serie entra al operador delta en el espacio crudo y con las fases de la calibracion (D7).
template <std::size_t T_MAX, std::size_t B_PURGA, class Obs, class Lil, class ConfigLil>
inline std::vector<EventoNuloLilyth<T_MAX>> correr_poblacion_via_a(
    const Obs& delta_original,
    const ark::metrologia::ConfiguracionObservador& cfg_delta,
    const ConfigLil& cfg_lil,
    const Lil& lil_plantilla,
    bool es_lineal,
    bool rama_purgada,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    const ContextoPurgaNulo<B_PURGA>& ctx,
    double x0_nivel,
    std::size_t n_cal,
    std::size_t B,
    std::size_t hilos,
    std::uint64_t semilla_base,
    std::size_t& sin_factor_total) {
    sin_factor_total = 0;
    if (B == 0 || n_cal == 0 || segs_calma.empty()) return {};

    std::vector<std::vector<EventoNuloLilyth<T_MAX>>> por_realizacion(B);
    std::vector<std::size_t> sin_factor_b(B, 0);

    std::mt19937_64 rng_compartido(semilla_base);
    std::mutex mtx_compartido;
    Obs delta_compartido = delta_original;

    ark::metrologia::familia::en_paralelo(B, hilos, [&](std::size_t b) {
        const bool semilla_independiente_por_indice = true;
        std::uint64_t sem_b = splitmix64(semilla_base ^ (0x9E3779B97F4A7C15ULL * (b + 1)));
        if (!semilla_independiente_por_indice) {
            std::lock_guard<std::mutex> lk(mtx_compartido);
            sem_b = rng_compartido() ^ static_cast<std::uint64_t>(std::hash<std::thread::id>{}(std::this_thread::get_id()));
        }
        std::mt19937_64 rng_b(sem_b);

        std::vector<double> est, x_sim;
        std::vector<std::size_t> fases;
        std::vector<std::int64_t> dias;
        construir_estacionaria_nula(es_lineal, cfg_lil.nulo_sin_memoria, segs_calma, m_lin, n_cal, rng_b, est);
        llevar_a_crudo<B_PURGA>(cfg_delta.integracion, rama_purgada, ctx, est, x0_nivel, 0, b, x_sim, fases, dias, sin_factor_b[b]);

        // Regla 7: cada realizacion parte de una copia NUEVA del original (la PurgaDiurna acumula estado).
        const bool copia_nueva_por_realizacion = true;
        Obs d_local = delta_original;
        if (!copia_nueva_por_realizacion) {
            std::lock_guard<std::mutex> lk(mtx_compartido);
            d_local = delta_compartido;
        }
        d_local.reiniciar_flujo();

        Lil lil_local = lil_plantilla;
        lil_local.reiniciar_flujo_conservando_purga();

        for (std::size_t t = 0; t < n_cal; ++t) {
            const auto s = d_local.paso(x_sim[t], NaN, true, fases[t], dias[t]);
            const auto o = lil_local.paso(s, x_sim[t], true, NaN, fases[t], dias[t]);
            for (std::size_t i = 0; i < o.ranuras_cerradas; ++i) {
                EventoNuloLilyth<T_MAX> ev{};
                if (extraer_evento_nulo<T_MAX>(lil_local.cerrada(i), rama_purgada, ev)) {
                    por_realizacion[b].push_back(ev);
                }
            }
        }
        if (!copia_nueva_por_realizacion) {
            std::lock_guard<std::mutex> lk(mtx_compartido);
            delta_compartido = d_local;
        }
    });

    std::vector<EventoNuloLilyth<T_MAX>> total;
    for (std::size_t b = 0; b < B; ++b) {
        total.insert(total.end(), por_realizacion[b].begin(), por_realizacion[b].end());
        sin_factor_total += sin_factor_b[b];
    }
    return total;
}
