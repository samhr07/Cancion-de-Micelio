// Ejecuta una realizacion nula con una o varias inyecciones candidatas en los tiempos `tiempos_iny`, verificando la guarda de
// solapamiento de ventanas (seccion 3.3 y control S2) y exigiendo que el operador delta dispare.
// Code (revision 2026-11-07):
//   - `sin_memoria` (opcion (b) aprobada por el Comandante): la base es un sustituto sin memoria de la calma y el perfil REAL (en el
//     espacio de la representacion) se escribe tal cual en t_iny..t_iny+l-1. Que la serie no responda es la hipotesis del nulo sin
//     memoria (en el nulo lineal ese mismo gesto es el mutante P8). En el nulo lineal se inyectan las innovaciones del perfil, como antes.
//   - D2: una apertura posterior al perfil mientras la ranura de la propia inyeccion sigue viva es una replica de ESE choque, no un
//     evento ajeno: no cuenta como solape. Ajenas son las aperturas dentro de la ventana fuera de [t_iny, t_iny + l] que no son replicas
//     de la ranura propia.
//   - Una inyeccion detectada cuya ranura acaba cortada por una replica no es "no detectada": se cuenta en `cortadas` y queda fuera,
//     igual que los eventos reales cortados (evaluar_nulos_representacion solo usa RESTITUIDA y CENSURADA_HORIZONTE en los dos lados).
//   - D7: la serie entra en el espacio crudo con las fases de la calibracion (llevar_a_crudo); NIVEL parte de x0, no de 100.
//   - `aplicar_guarda` = false solo para el control S2 (mide lo que pasaria sin la guarda con el MISMO instrumento).
template <std::size_t T_MAX, std::size_t B_PURGA, class Obs, class Lil>
inline ResultadoInyeccionesRealizacion<T_MAX> ejecutar_inyecciones_en_realizacion(
    const Obs& delta_original,
    const Lil& lil_plantilla,
    ark::metrologia::Integracion integracion,
    bool rama_purgada,
    bool sin_memoria,
    std::uint32_t codigo_sustituto,
    const std::vector<std::vector<double>>& segs_calma,
    const ModeloLinealNulo& m_lin,
    const ContextoPurgaNulo<B_PURGA>& ctx,
    double x0_nivel,
    std::span<const double> perfil_est,
    std::span<const double> e_iny,
    int signo_real,
    std::size_t w_ancla,
    std::size_t guarda,
    std::size_t t_max_ef,
    std::span<const std::int64_t> tiempos_iny,
    std::size_t n_sim,
    std::mt19937_64& rng,
    bool aplicar_guarda = true) {
    ResultadoInyeccionesRealizacion<T_MAX> res{};
    res.candidatas = tiempos_iny.size();
    const std::size_t l = perfil_est.size();
    if (l == 0 || tiempos_iny.empty() || n_sim == 0) return res;
    if (!sin_memoria && (!m_lin.valido || e_iny.size() != l)) return res;
    if (sin_memoria && segs_calma.empty()) return res;

    // 1. Ventanas de exclusion [t_iny - (W_ancla + guarda), t_iny + l + T_MAX] y solapes entre inyecciones.
    std::vector<VentanaInyeccion> ventanas(tiempos_iny.size());
    std::vector<bool> solapada(tiempos_iny.size(), false);
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        ventanas[i] = construir_ventana_inyeccion(tiempos_iny[i], l, w_ancla, guarda, T_MAX);
    }
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        for (std::size_t j = i + 1; j < tiempos_iny.size(); ++j) {
            if (ventanas_se_solapan(ventanas[i], ventanas[j])) { solapada[i] = true; solapada[j] = true; }
        }
    }

    // 2. Parte estacionaria con la inyeccion.
    std::vector<double> x_est(n_sim, 0.0);
    if (sin_memoria) {
        construir_estacionaria_nula(false, codigo_sustituto, segs_calma, m_lin, n_sim, rng, x_est);
        for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
            if (tiempos_iny[i] < 0) continue;
            const std::size_t t0 = static_cast<std::size_t>(tiempos_iny[i]);
            for (std::size_t k = 0; k < l && (t0 + k) < n_sim; ++k) x_est[t0 + k] = perfil_est[k];
        }
    } else {
        const std::size_t burn = m_lin.burn_in;
        const std::size_t total = burn + n_sim;
        std::vector<double> eps(total, 0.0);
        if (!m_lin.residuos.empty()) {
            std::uniform_int_distribution<std::size_t> dist(0, m_lin.residuos.size() - 1);
            for (std::size_t t = 0; t < total; ++t) eps[t] = m_lin.residuos[dist(rng)] - m_lin.media_residuos;
        } else {
            std::normal_distribution<double> dist(0.0, std::sqrt(m_lin.sigma2_innov));
            for (std::size_t t = 0; t < total; ++t) eps[t] = dist(rng);
        }
        const bool inyectar_en_innovaciones = true;
        if (inyectar_en_innovaciones) {
            for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
                if (tiempos_iny[i] < 0) continue;
                const std::size_t t0 = static_cast<std::size_t>(tiempos_iny[i]);
                for (std::size_t k = 0; k < l && (t0 + k) < n_sim; ++k) eps[burn + t0 + k] = e_iny[k];
            }
        }
        std::vector<double> u(total, 0.0);
        for (std::size_t t = 0; t < total; ++t) {
            double acc = eps[t];
            const std::size_t i_max = std::min(m_lin.p, t);
            for (std::size_t i = 1; i <= i_max; ++i) acc += m_lin.phi[i - 1] * u[t - i];
            u[t] = acc;
        }
        if (m_lin.modo_usado == ModoNuloLineal::ARFIMA_WHITTLE && m_lin.d != 0.0) {
            std::vector<double> pi_inv(total, 0.0);
            ark::metrologia::elw::coeficientes(-m_lin.d, total, pi_inv.data());
            for (std::size_t t = 0; t < n_sim; ++t) {
                const std::size_t idx = burn + t;
                double s = 0.0;
                for (std::size_t k = 0; k <= idx; ++k) s += pi_inv[k] * u[idx - k];
                x_est[t] = m_lin.mu + s;
            }
        } else {
            for (std::size_t t = 0; t < n_sim; ++t) x_est[t] = m_lin.mu + u[burn + t];
        }
        if (!inyectar_en_innovaciones) {
            // Mutante M08 (P8): superpone el perfil sobre la serie fija sin propagarlo por el filtro lineal
            for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
                if (tiempos_iny[i] < 0) continue;
                const std::size_t t0 = static_cast<std::size_t>(tiempos_iny[i]);
                for (std::size_t k = 0; k < l && (t0 + k) < n_sim; ++k) x_est[t0 + k] = perfil_est[k];
            }
        }
    }

    // 3. Al espacio crudo con las fases de la calibracion (D7).
    std::vector<double> x_sim;
    std::vector<std::size_t> fases;
    std::vector<std::int64_t> dias;
    std::size_t sin_factor = 0;
    std::size_t inicio = 0;
    if (ctx.purga_activa && !ctx.fases.empty()) {
        std::uniform_int_distribution<std::size_t> d_ini(0, ctx.fases.size() - 1);
        inicio = d_ini(rng);
    }
    llevar_a_crudo<B_PURGA>(integracion, rama_purgada, ctx, x_est, x0_nivel, inicio, 0, x_sim, fases, dias, sin_factor);
    res.sin_factor_purga = sin_factor;

    // 4. Copia limpia del operador delta y de ObservadorLilyth sobre la realizacion; se registran todas las ranuras.
    Obs d_local = delta_original;
    d_local.reiniciar_flujo();
    Lil lil_local = lil_plantilla;
    lil_local.reiniciar_flujo_conservando_purga();

    struct Cierre {
        std::int64_t t_apertura = -1;
        std::int64_t t_cierre = -1;
        int signo = 0;
        bool extraible = false;
        std::uint8_t desenlace = 0;
        EventoNuloLilyth<T_MAX> ev{};
    };
    std::vector<Cierre> cierres;
    std::vector<std::int64_t> aperturas_todas;
    for (std::size_t t = 0; t < n_sim; ++t) {
        const auto s = d_local.paso(x_sim[t], NaN, true, fases[t], dias[t]);
        const auto o = lil_local.paso(s, x_sim[t], true, NaN, fases[t], dias[t]);
        if (o.abrio_ranura) aperturas_todas.push_back(o.t);
        for (std::size_t c = 0; c < o.ranuras_cerradas; ++c) {
            const auto& rc = lil_local.cerrada(c);
            Cierre ci{};
            ci.t_apertura = rc.t_apertura;
            ci.t_cierre = o.t;
            ci.signo = rc.signo;
            const auto& rama = rama_purgada ? rc.purgada : rc.cruda;
            ci.desenlace = static_cast<std::uint8_t>(rama.desenlace);
            ci.extraible = extraer_evento_nulo<T_MAX>(rc, rama_purgada, ci.ev);
            cierres.push_back(ci);
        }
    }
    // Ranuras aun abiertas al final: siguen vivas hasta n_sim.
    for (const auto& r : lil_local.ranuras()) {
        if (!r.ocupada) continue;
        Cierre ci{};
        ci.t_apertura = r.t_apertura;
        ci.t_cierre = static_cast<std::int64_t>(n_sim);
        ci.signo = r.signo;
        cierres.push_back(ci);
    }

    // 5. Ranura propia de cada inyeccion: la primera abierta en [t_iny, t_iny + l] con el signo del evento real.
    const bool replicas_propias_son_de_la_inyeccion = true;   // D2 (el mutante M13 de Code lo pone en false)
    for (std::size_t i = 0; i < tiempos_iny.size(); ++i) {
        const std::int64_t t0 = tiempos_iny[i];
        const std::int64_t t1 = t0 + static_cast<std::int64_t>(l);
        const Cierre* propia = nullptr;
        for (const auto& ci : cierres) {
            if (ci.t_apertura >= t0 && ci.t_apertura <= t1 && ci.signo == signo_real) {
                if (propia == nullptr || ci.t_apertura < propia->t_apertura) propia = &ci;
            }
        }
        // Eventos ajenos dentro de la ventana.
        const VentanaInyeccion& vi = ventanas[i];
        for (std::int64_t ta : aperturas_todas) {
            if (ta >= t0 && ta <= t1) continue;
            if (ta < vi.ini || ta > vi.fin) continue;
            const bool replica_propia = propia != nullptr && ta > t1 && ta <= propia->t_cierre;
            if (replica_propia && replicas_propias_son_de_la_inyeccion) continue;
            solapada[i] = true;
        }
        if (solapada[i] && aplicar_guarda) {
            ++res.solapes;
            continue;
        }
        if (solapada[i]) ++res.solapes_sin_guarda;
        const bool exige_disparo_detector = true;
        if (propia != nullptr) {
            if (propia->extraible) {
                ++res.detectadas;
                res.eventos_admitidos.push_back(propia->ev);
            } else if (propia->desenlace == static_cast<std::uint8_t>(3)) {
                ++res.cortadas;          // detectada y cortada por una replica: fuera, igual que los reales cortados
            } else {
                ++res.no_medidas;        // detectada pero sin medida (ancla, signo, horizonte abierto al final)
            }
        } else if (!exige_disparo_detector) {
            // Mutante M03 (Argon): mide la restitucion desde t_iny + l - 1 sin exigir que el detector dispare
            const std::size_t ta = static_cast<std::size_t>(std::max<std::int64_t>(1, t0 + static_cast<std::int64_t>(l) - 1));
            if (ta + t_max_ef < n_sim && ta >= w_ancla + guarda) {
                std::vector<double> w(x_sim.begin() + static_cast<std::ptrdiff_t>(ta - w_ancla - guarda),
                                      x_sim.begin() + static_cast<std::ptrdiff_t>(ta - guarda));
                const double y_anc = ark::metrologia::est::cuantil(w, 0.5);
                const double dy = x_sim[ta] - y_anc;
                if (std::isfinite(dy) && dy * static_cast<double>(signo_real) > 0.0) {
                    EventoNuloLilyth<T_MAX> ev_forzado{};
                    ev_forzado.signo = signo_real;
                    ev_forzado.y_ancla = y_anc;
                    ev_forzado.dy = dy;
                    ev_forzado.eta.fill(NaN);
                    for (std::size_t tau = 1; tau <= t_max_ef && (ta + tau) < n_sim; ++tau) {
                        const double e = (x_sim[ta + tau] - y_anc) / dy;
                        ev_forzado.eta[tau - 1] = e;
                        ev_forzado.pasos_observados = tau;
                        if (e <= 0.5 && !std::isfinite(ev_forzado.tau_rest_primer_toque)) {
                            ev_forzado.tau_rest_primer_toque = static_cast<double>(tau);
                            ev_forzado.tau_rest_consolidado = static_cast<double>(tau);
                        }
                    }
                    ev_forzado.desenlace = std::isfinite(ev_forzado.tau_rest_consolidado) ? 2u : 4u;
                    ev_forzado.tau_fin = std::isfinite(ev_forzado.tau_rest_consolidado)
                                             ? static_cast<std::size_t>(ev_forzado.tau_rest_consolidado)
                                             : t_max_ef;
                    ++res.detectadas;
                    res.eventos_admitidos.push_back(ev_forzado);
                } else {
                    ++res.no_detectadas;
                }
            } else {
                ++res.no_detectadas;
            }
        } else {
            ++res.no_detectadas;
        }
    }
    return res;
}
