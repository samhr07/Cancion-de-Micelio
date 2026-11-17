    // S2: ventanas solapadas por debajo y por encima de W_ancla + guarda + l + T_MAX
    //     y efecto fisico de contaminacion del ancla.
    // Code (revision D2/D5, 2026-11-07):
    //   - el mundo de Argon (ARL0 = 2000 sobre realizaciones de 1200 pasos) da alarmas espontaneas con probabilidad ~ 0.45: una de ellas
    //     (t = 470) caia en la ventana del ancla de la segunda inyeccion "supercritica" y la guarda la descartaba CON RAZON. El mundo pasa
    //     a ARL0 = 1e6 (sin espontaneas en la practica) para que la separacion sea lo unico que decida;
    //   - la parte fisica corre el MISMO instrumento con y sin la guarda (aplicar_guarda = false), en lugar de mover el ancla a mano.
    // -------------------------------------------------------------------------
    {
        const double phi = 0.85;
        const auto x_cal = generar_ar1(phi, 4096, 1001);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 1.0e6;   // sin falsas alarmas espontaneas en las realizaciones de la prueba
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d = std::make_unique<Obs>();
        d->calibrar(cd, x_cal);

        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 2;
        cl.hilos_nulo = 1;
        lil->calibrar(*d, cl, cd, x_cal);

        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, d->ficha_nivel(0).nulo.tau_int)));
        const std::size_t l = 4;
        const std::size_t dist_min = w_ancla + guarda + l + TMAXC;   // separacion exacta del borde de solape

        const std::array<double, 8> hist_prev{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const std::array<double, 4> perfil{4.0, 6.0, 7.5, 8.5};
        const std::size_t n_sim = 1200;

        // Geometria exacta de las ventanas (sin instrumento).
        const auto v1 = detalle_nulo::construir_ventana_inyeccion(250, l, w_ancla, guarda, TMAXC);
        const bool geo_ok =
            detalle_nulo::ventanas_se_solapan(v1, detalle_nulo::construir_ventana_inyeccion(250 + static_cast<std::int64_t>(dist_min) - 20, l, w_ancla, guarda, TMAXC)) &&
            detalle_nulo::ventanas_se_solapan(v1, detalle_nulo::construir_ventana_inyeccion(250 + static_cast<std::int64_t>(dist_min), l, w_ancla, guarda, TMAXC)) &&
            !detalle_nulo::ventanas_se_solapan(v1, detalle_nulo::construir_ventana_inyeccion(250 + static_cast<std::int64_t>(dist_min) + 20, l, w_ancla, guarda, TMAXC));

        const std::array<std::int64_t, 2> t_sub{250, static_cast<std::int64_t>(250 + dist_min - 20)};
        const auto r_sub = lil->evaluar_inyecciones_realizacion(false, hist_prev, perfil, +1, w_ancla, guarda, t_sub, n_sim, 1010);
        const std::array<std::int64_t, 2> t_borde{250, static_cast<std::int64_t>(250 + dist_min)};
        const auto r_borde = lil->evaluar_inyecciones_realizacion(false, hist_prev, perfil, +1, w_ancla, guarda, t_borde, n_sim, 1011);
        const std::array<std::int64_t, 2> t_sup{250, static_cast<std::int64_t>(250 + dist_min + 20)};
        const auto r_sup = lil->evaluar_inyecciones_realizacion(false, hist_prev, perfil, +1, w_ancla, guarda, t_sup, n_sim, 1012);

        // Contaminacion fisica: una meseta larga (30 barras a +5 sd) y una segunda inyeccion igual 60 barras despues, cuando la primera
        // ya restituyo: la ventana del ancla de la segunda contiene la meseta. Con la guarda la segunda se descarta; sin ella entra con su
        // ancla sesgada. Se compara su tau con el de inyecciones aisladas del mismo perfil.
        std::vector<double> meseta(30, 5.0);
        std::vector<double> hp_meseta(8, 0.0);
        std::vector<double> tau_limpio, tau_contaminado;
        std::size_t descartadas_con_guarda = 0, pares = 0;
        for (std::size_t k = 0; k < 40; ++k) {
            const std::array<std::int64_t, 1> t_solo{300};
            const auto rl = lil->evaluar_inyecciones_realizacion(false, hp_meseta, meseta, +1, w_ancla, guarda, t_solo, 900, 1100 + k);
            for (const auto& ev : rl.eventos_admitidos) if (std::isfinite(ev.tau_rest_consolidado)) tau_limpio.push_back(ev.tau_rest_consolidado);

            const std::array<std::int64_t, 2> t_par{300, 360};
            const auto rg = lil->evaluar_inyecciones_realizacion(false, hp_meseta, meseta, +1, w_ancla, guarda, t_par, 900, 2100 + k, false, true);
            const auto rs = lil->evaluar_inyecciones_realizacion(false, hp_meseta, meseta, +1, w_ancla, guarda, t_par, 900, 2100 + k, false, false);
            ++pares;
            descartadas_con_guarda += rg.solapes;
            for (std::size_t q = 0; q < rs.eventos_admitidos.size(); ++q) {
                if (rs.indice_inyeccion[q] == 1 && std::isfinite(rs.eventos_admitidos[q].tau_rest_consolidado))
                    tau_contaminado.push_back(rs.eventos_admitidos[q].tau_rest_consolidado);
            }
        }
        const double med_limpia = tau_limpio.empty() ? NaN : est::cuantil(tau_limpio, 0.5);
        const double med_contam = tau_contaminado.empty() ? NaN : est::cuantil(tau_contaminado, 0.5);

        std::printf("      S2: dist_min=%zu geo=%d | sub(sol=%zu, adm=%zu) borde(sol=%zu, adm=%zu) sup(sol=%zu, adm=%zu)\n",
                    dist_min, int(geo_ok), r_sub.solapes, r_sub.eventos_admitidos.size(), r_borde.solapes, r_borde.eventos_admitidos.size(),
                    r_sup.solapes, r_sup.eventos_admitidos.size());
        std::printf("      S2 contaminacion: con guarda %zu de %zu pares descartados (2 por par); sin guarda segunda inyeccion tau_med=%.1f (n=%zu) "
                    "contra aislada %.1f (n=%zu)\n",
                    descartadas_con_guarda, pares, med_contam, tau_contaminado.size(), med_limpia, tau_limpio.size());

        chk(geo_ok && r_sub.solapes == 2 && r_sub.eventos_admitidos.empty() &&
            r_borde.solapes == 2 && r_borde.eventos_admitidos.empty() &&
            r_sup.solapes == 0 && r_sup.eventos_admitidos.size() == 2 &&
            descartadas_con_guarda == 2 * pares && !tau_contaminado.empty() && !tau_limpio.empty() &&
            med_contam != med_limpia,
            "S2", "ventanas solapadas: subcritica y borde se descartan, supercritica entra; sin guarda la contaminacion cambia tau");
    }
