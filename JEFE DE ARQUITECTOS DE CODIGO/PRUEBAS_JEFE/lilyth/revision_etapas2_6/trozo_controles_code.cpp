    // =========================================================================================================================
    // Controles de Code (revision 2026-11-07)
    // =========================================================================================================================
    // D1: AR(1) phi = 0.8 con choques (el mundo de P2). Ningun evento se evalua con menos de m_min nulos; cada evaluado por (a) tiene
    //     su excedencia exacta (escala, signo, g >= g_k); los que no llegan en (a) usan la via (b) y lo marcan. Y la via (b) SIN MEMORIA
    //     se ejerce: el nulo sin memoria de ese mundo no llega a m_min en (a) (el detector, calibrado con memoria, casi no alarma sobre
    //     el sustituto permutado).
    {
        const double phi = 0.8;
        const auto x_cal = generar_ar1(phi, 1 << 13, 4401);
        std::vector<std::size_t> tch;
        for (std::size_t t = 400; t + 200 < 6000; t += 300) tch.push_back(t);
        const auto x_test = generar_ar1_con_choques(phi, 6000, tch, 3.0, 4, 4402);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 2000.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d = std::make_unique<Obs>();
        d->calibrar(cd, x_cal);
        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.b_max_nulo = 8;
        cl.hilos_nulo = 4;
        cl.b_mc_resumen = 999;
        lil->calibrar(*d, cl, cd, x_cal);
        for (double v : x_test) { const auto s = d->paso(v, NaN, true); (void)lil->paso(s, v, true); }
        lil->completar_nulos();
        const auto& fn = lil->ficha_nulo_cruda();
        std::size_t ev_iid_a = 0, ev_iid_b = 0, ev_iid_sin = 0, ev_ons_a = 0, ev_ons_b = 0, ev_ons_sin = 0, bajo_mmin = 0, excedencia_mal = 0;
        for (const auto& r : lil->cerradas_historico()) {
            if (r.cruda.desenlace != DesenlaceRanura::RESTITUIDA && r.cruda.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
            for (int q = 0; q < 2; ++q) {
                const auto& ref = q == 0 ? r.cruda.nulo_iid : r.cruda.nulo_onsager;
                const auto& pool = q == 0 ? fn.nulos_sin_memoria : fn.nulos_lineal;
                if (!ref.evaluado) { ++(q == 0 ? ev_iid_sin : ev_ons_sin); continue; }
                if (ref.n_nulos_usados < fn.m_min) ++bajo_mmin;
                if (ref.nulo_por_inyeccion) { ++(q == 0 ? ev_iid_b : ev_ons_b); continue; }
                ++(q == 0 ? ev_iid_a : ev_ons_a);
                std::size_t n_ref = 0;
                for (const auto& e : pool) if (e.escala_origen == r.escala_origen && e.signo == r.signo && e.g_disparo >= r.g_disparo) ++n_ref;
                if (n_ref != ref.n_nulos_usados) ++excedencia_mal;
            }
        }
        std::printf("      D1: iid (a=%zu b=%zu sin=%zu) Onsager (a=%zu b=%zu sin=%zu) | nulos en (a): iid=%zu lineal=%zu | bajo m_min=%zu "
                    "excedencia mal=%zu | p_iid=%.4f p_ons=%.4f\n",
                    ev_iid_a, ev_iid_b, ev_iid_sin, ev_ons_a, ev_ons_b, ev_ons_sin, fn.nulos_sin_memoria.size(), fn.nulos_lineal.size(),
                    bajo_mmin, excedencia_mal, fn.contraste_iid.p_bilateral, fn.contraste_onsager.p_bilateral);
        chk(bajo_mmin == 0 && excedencia_mal == 0 && ev_iid_b > 0 && (ev_iid_a + ev_iid_b) > 0 && (ev_ons_a + ev_ons_b) > 0 &&
            fn.contraste_iid.valido && fn.contraste_iid.p_bilateral < cl.alfa,
            "D1", "sin relleno: ningun nulo con < m_min, excedencia exacta en (a), via (b) sin memoria ejercida y Delta_tau_iid > 0 se ve");
    }

    // D1b: la via (b) sin memoria escribe el perfil real en un sustituto permutado y exige el disparo. Un perfil fuerte se detecta y
    //      restituye en tau = 1 (la serie sin memoria no responde); un perfil invisible da p_det = 0 y queda sin nulo (no se fabrica).
    {
        const auto x_cal = generar_ar1(0.6, 4096, 4501);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 1.0e5;
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
        const auto& fn = lil->ficha_nulo_cruda();
        const std::array<double, 4> hp{0.0, 0.0, 0.0, 0.0};
        const std::array<double, 4> fuerte{6.0, 8.0, 9.0, 9.5};
        const std::array<double, 3> invisible{0.01, 0.02, 0.01};
        ReferenciaNuloEvento rf{}, ri{};
        std::vector<EventoNuloLilyth<TMAXC>> ef, ei;
        std::size_t sol = 0;
        detalle_nulo::ejecutar_via_b_evento<TMAXC>(*d, *lil, Integracion::FLUJO, false, true, cl.nulo_sin_memoria, fn.segmentos_calma,
            fn.modelo_lineal, lil->contexto_purga_nulo(), 0.0, hp, fuerte, +1, 120, 2, fn.m_min, lil->configuracion().t_max_efectivo,
            4502, rf, ef, sol);
        detalle_nulo::ejecutar_via_b_evento<TMAXC>(*d, *lil, Integracion::FLUJO, false, true, cl.nulo_sin_memoria, fn.segmentos_calma,
            fn.modelo_lineal, lil->contexto_purga_nulo(), 0.0, hp, invisible, +1, 120, 2, fn.m_min, lil->configuracion().t_max_efectivo,
            4503, ri, ei, sol);
        std::size_t tau1 = 0;
        for (const auto& e : ef) if (e.tau_rest_consolidado == 1.0) ++tau1;
        std::printf("      D1b: fuerte n_b=%zu intentos=%zu tau=1 en %zu | invisible n_b=%zu intentos=%zu sin_nulo_b=%d\n",
                    rf.n_nulos_b, rf.intentos_b, tau1, ri.n_nulos_b, ri.intentos_b, int(ri.sin_nulo_b));
        chk(rf.n_nulos_b >= fn.m_min && tau1 * 2 >= ef.size() && ri.sin_nulo_b && ri.n_nulos_b == 0 && ri.intentos_b == fn.m_min,
            "D1b", "via (b) sin memoria: perfil fuerte se detecta y la serie no responde; perfil invisible queda sin nulo (sin fabricar)");
    }

    // D7: con purga activa, (a) el espejo invertible da lo mismo que PurgaDiurna bit a bit en la calibracion y en una continuacion,
    //     (b) llevar una serie del espacio purgado al crudo y volver a purgarla la devuelve (<= 1e-12 relativo), y (c) Lilyth calibra los
    //     nulos de la rama PURGADA sin discrepancias del espejo.
    {
        std::mt19937_64 rng(4601);
        std::normal_distribution<double> nd(0.0, 1.0);
        const std::size_t dia_pasos = 100, n_cal = 4600;
        std::vector<double> xc(n_cal);
        std::vector<std::uint8_t> vc(n_cal, 1);
        std::vector<std::size_t> fc(n_cal);
        std::vector<std::int64_t> dc(n_cal);
        double ar = 0.0;
        for (std::size_t t = 0; t < n_cal; ++t) {
            fc[t] = (t % dia_pasos) / 25;
            dc[t] = static_cast<std::int64_t>(t / dia_pasos);
            ar = 0.7 * ar + nd(rng);
            xc[t] = (fc[t] == 0 ? 3.0 : 1.0) * ar;
        }
        // (a) espejo contra PurgaDiurna
        ark::metrologia::PurgaDiurna<1440> pd(4, 20);
        detalle_nulo::EspejoPurgaNulo<1440> ep(4, 20);
        std::size_t dif_a = 0, n_a = 0;
        for (std::size_t t = 0; t < n_cal + 2000; ++t) {
            const std::size_t f = t < n_cal ? fc[t] : (t % dia_pasos) / 25;
            const std::int64_t dd = t < n_cal ? dc[t] : static_cast<std::int64_t>(t / dia_pasos);
            const double y = t < n_cal ? xc[t] : nd(rng);
            const double fa = ep.factor(f, dd);
            bool v1 = false, v2 = false;
            const double a = pd.paso(y, f, dd, v1);
            const double b = ep.paso(y, f, dd, v2);
            ++n_a;
            if (v1 != v2 || !igual(a, b)) ++dif_a;
            if (v1 && !(std::abs(y * fa - a) <= 1e-12 * std::abs(a))) ++dif_a;   // factor() = lo que paso() aplica
        }
        // (c) Lilyth con PURGADA
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 2400.0;
        cd.purga = true;
        cd.bins_fase = 4;
        cd.dias_minimos_purga = 20;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto d = std::make_unique<Obs>();
        d->calibrar(cd, xc, {}, vc, fc, dc);
        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::PURGADA;
        cl.b_max_nulo = 4;
        cl.hilos_nulo = 1;
        const bool cal = lil->calibrar(*d, cl, cd, xc, vc, fc, dc);
        // (b) ida y vuelta con el contexto de Lilyth
        const auto& ctx = lil->contexto_purga_nulo();
        const auto& fn = lil->ficha_nulo_purgada();
        std::mt19937_64 rng_b(4602);
        std::vector<double> est, x;
        std::vector<std::size_t> ff;
        std::vector<std::int64_t> ddd;
        std::size_t sf = 0;
        detalle_nulo::construir_estacionaria_nula(false, cl.nulo_sin_memoria, fn.segmentos_calma, fn.modelo_lineal, 3000, rng_b, est);
        detalle_nulo::llevar_a_crudo(Integracion::FLUJO, true, ctx, est, 0.0, 0, 0, x, ff, ddd, sf);
        auto esp = ctx.espejo_calibrado;
        double max_rel = 0.0;
        std::size_t validos = 0;
        for (std::size_t t = 0; t < x.size(); ++t) {
            bool v = false;
            const double p = esp.paso(x[t], ff[t], ddd[t], v);
            if (!v) continue;
            ++validos;
            max_rel = std::max(max_rel, std::abs(p - est[t]) / std::max(1e-300, std::abs(est[t])));
        }
        std::printf("      D7: espejo vs PurgaDiurna %zu difs en %zu | ida y vuelta max_rel=%.2e (%zu validos, sin_factor=%zu) | "
                    "Lilyth PURGADA calibra=%d discrepancias=%zu nulos(a): iid=%zu lin=%zu\n",
                    dif_a, n_a, max_rel, validos, sf, int(cal), lil->espejo_purga_discrepancias(),
                    fn.nulos_sin_memoria.size(), fn.nulos_lineal.size());
        chk(dif_a == 0 && validos == x.size() && max_rel <= 1e-12 && cal && lil->espejo_purga_discrepancias() == 0 && fn.calibrada,
            "D7", "purga en los nulos: espejo = PurgaDiurna bit a bit, ida y vuelta exacta y la rama PURGADA calibra sus nulos");
    }

