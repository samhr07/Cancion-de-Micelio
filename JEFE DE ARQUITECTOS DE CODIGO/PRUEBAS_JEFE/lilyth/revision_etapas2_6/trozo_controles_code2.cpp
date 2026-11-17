    // D1c: pool del Monte Carlo = PIT de punto medio de cada nulo contra su excedencia SIN el mismo, solo con m_min (referencia
    //      escrita aparte sobre un pool sintetico con escalas, signos, g y tau conocidos).
    {
        std::mt19937_64 rng(4701);
        std::uniform_int_distribution<int> di(0, 1);
        std::uniform_real_distribution<double> dg(1.0, 10.0);
        std::uniform_int_distribution<int> dt(1, 12);
        std::vector<EventoNuloLilyth<TMAXC>> pool(120);
        for (auto& e : pool) {
            e.escala_origen = di(rng);
            e.signo = di(rng) ? +1 : -1;
            e.g_disparo = dg(rng);
            e.tau_rest_consolidado = static_cast<double>(dt(rng));
        }
        const std::size_t mm = 8, tm = 10;   // tau > 10 se codifica como 11 (censurado)
        std::vector<double> u_mod;
        detalle_nulo::pool_pit_loo<TMAXC>(pool, mm, tm, u_mod);
        std::vector<double> u_ref;
        for (std::size_t i = 0; i < pool.size(); ++i) {
            std::vector<double> v;
            for (std::size_t j = 0; j < pool.size(); ++j) {
                if (j == i) continue;
                if (pool[j].escala_origen == pool[i].escala_origen && pool[j].signo == pool[i].signo && pool[j].g_disparo >= pool[i].g_disparo)
                    v.push_back(std::min(pool[j].tau_rest_consolidado, static_cast<double>(tm + 1)));
            }
            if (v.size() < mm) continue;
            const double x = std::min(pool[i].tau_rest_consolidado, static_cast<double>(tm + 1));
            std::size_t lt = 0, eq = 0;
            for (double y : v) { if (y < x) ++lt; else if (y == x) ++eq; }
            u_ref.push_back((static_cast<double>(lt) + 0.5 * static_cast<double>(eq)) / static_cast<double>(v.size()));
        }
        bool igual_todo = u_mod.size() == u_ref.size() && !u_ref.empty();
        for (std::size_t i = 0; igual_todo && i < u_ref.size(); ++i) igual_todo = std::abs(u_mod[i] - u_ref[i]) <= 1e-15;
        std::printf("      D1c: pool %zu nulos -> %zu PIT (referencia %zu)\n", pool.size(), u_mod.size(), u_ref.size());
        chk(igual_todo, "D1c", "pool del Monte Carlo = PIT contra la excedencia sin el mismo, solo con m_min (= referencia)");
    }

