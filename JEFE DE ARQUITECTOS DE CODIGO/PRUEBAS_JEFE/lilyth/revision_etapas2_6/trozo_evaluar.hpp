    // Code (revision D1, 2026-11-07, opcion (b) aprobada por el Comandante): para cada evento real y cada nulo,
    //   - via (a) si la excedencia (misma escala, mismo signo, g >= g_k) tiene al menos m_min eventos;
    //   - si no, via (b): en el nulo lineal, inyeccion de las innovaciones del perfil en el modelo; en el nulo SIN MEMORIA, el perfil real
    //     escrito en un sustituto sin memoria de la calma, con el mismo detector;
    //   - si (b) tampoco llega a m_min, el evento queda SIN nulo (no se rellena) y la ficha dice por que.
    // El pool del Monte Carlo usa la misma regla: el PIT de un evento nulo se calcula contra su excedencia sin el mismo, y solo si esta
    // tiene m_min eventos.
    void evaluar_una_referencia(FichaNuloRepresentacion<T_MAX>& fn, bool es_purgada, bool sin_memoria, Ranura& r, std::size_t idx,
                                ReferenciaNuloEvento& ref, std::vector<double>& u_pool_extra) {
        auto& rama = es_purgada ? r.purgada : r.cruda;
        const auto& pool = sin_memoria ? fn.nulos_sin_memoria : fn.nulos_lineal;
        auto esp = detalle_nulo::filtrar_espontaneos_excedencia<T_MAX>(r.escala_origen, r.signo, r.g_disparo, pool);
        ref.n_nulos_a = esp.size();
        std::vector<double> taus;
        if (!cfg_.forzar_via_b && esp.size() >= fn.m_min) {
            ref.nulo_por_inyeccion = false;
            ref.motivo = sin_memoria ? "via (a) sin memoria: excedencia suficiente" : "via (a): excedencia espontanea suficiente";
            for (const auto* p : esp) taus.push_back(detalle_nulo::codificar_tau(p->tau_rest_consolidado, cfg_.t_max_efectivo));
        } else {
            const std::span<const double> hp = es_purgada
                ? std::span<const double>(r.hist_prev_purgada_est.data(), r.n_hist_prev_purgada)
                : std::span<const double>(r.hist_prev_cruda_est.data(), r.n_hist_prev_cruda);
            const std::span<const double> pf = es_purgada
                ? std::span<const double>(r.perfil_purgada_est.data(), r.n_perfil_purgada)
                : std::span<const double>(r.perfil_cruda_est.data(), r.n_perfil_cruda);
            std::vector<EventoNuloLilyth<T_MAX>> ev_b;
            const std::uint64_t sal = sin_memoria ? 0x51ED270B27C3A6E5ULL : 0x9E3779B97F4A7C15ULL;
            detalle_nulo::ejecutar_via_b_evento<T_MAX, Rasgos::b_max>(
                delta_calibrado_, *this, integracion_, es_purgada, sin_memoria, cfg_.nulo_sin_memoria, fn.segmentos_calma,
                fn.modelo_lineal, ctx_purga_, x0_nivel_, hp, pf, r.signo, r.w_ancla, r.guarda_ancla, fn.m_min, cfg_.t_max_efectivo,
                cfg_.semilla_nulo ^ (sal * (idx + 1)), ref, ev_b, cont_.solapes);
            r.nulo_por_inyeccion = true;
            for (const auto& ev : ev_b) taus.push_back(detalle_nulo::codificar_tau(ev.tau_rest_consolidado, cfg_.t_max_efectivo));
            if (taus.size() < fn.m_min) {
                ref.sin_nulo_b = true;
                ref.evaluado = false;
                ref.motivo = sin_memoria ? "sin nulo: (a) y (b) sin memoria no llegan a m_min" : "sin nulo: (a) y (b) no llegan a m_min";
                return;
            }
            for (std::size_t i = 0; i < taus.size(); ++i) {
                std::vector<double> resto;
                for (std::size_t j = 0; j < taus.size(); ++j) if (j != i) resto.push_back(taus[j]);
                u_pool_extra.push_back(detalle_nulo::pit_punto_medio(taus[i], resto));
            }
        }
        detalle_nulo::llenar_estadisticas_referencia(ref, rama.tau_rest_consolidado, cfg_.t_max_efectivo, taus);
    }

    void pool_loo(const std::vector<EventoNuloLilyth<T_MAX>>& pool, std::size_t m_min, std::vector<double>& u_pool) const {
        for (std::size_t i = 0; i < pool.size(); ++i) {
            const auto& ev = pool[i];
            std::size_t n_a = 0;
            const auto taus_i = detalle_nulo::construir_muestra_nula_emparejada<T_MAX>(
                ev.escala_origen, ev.signo, ev.g_disparo, pool, m_min, cfg_.t_max_efectivo, i, n_a);
            if (taus_i.size() < m_min) continue;
            const double tc = detalle_nulo::codificar_tau(ev.tau_rest_consolidado, cfg_.t_max_efectivo);
            u_pool.push_back(detalle_nulo::pit_punto_medio(tc, taus_i));
        }
    }

    void evaluar_nulos_representacion(FichaNuloRepresentacion<T_MAX>& fn, bool es_purgada) {
        std::vector<double> u_iid, dt_iid, u_ons, dt_ons, u_pool_iid, u_pool_ons;
        std::vector<std::array<double, T_MAX>> etas_reales;
        std::vector<std::size_t> pasos_reales;
        pool_loo(fn.nulos_sin_memoria, fn.m_min, u_pool_iid);
        pool_loo(fn.nulos_lineal, fn.m_min, u_pool_ons);

        for (std::size_t idx = 0; idx < cerradas_historico_.size(); ++idx) {
            Ranura& r = cerradas_historico_[idx];
            auto& rama = es_purgada ? r.purgada : r.cruda;
            if (rama.desenlace != DesenlaceRanura::RESTITUIDA && rama.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE) continue;
            etas_reales.push_back(rama.eta);
            pasos_reales.push_back(rama.pasos_observados);

            ReferenciaNuloEvento ref_iid{};
            evaluar_una_referencia(fn, es_purgada, true, r, idx, ref_iid, u_pool_iid);
            rama.nulo_iid = ref_iid;
            if (ref_iid.evaluado) { u_iid.push_back(ref_iid.u_pit); dt_iid.push_back(ref_iid.delta_tau); }

            ReferenciaNuloEvento ref_ons{};
            if (fn.modelo_lineal.valido) evaluar_una_referencia(fn, es_purgada, false, r, idx, ref_ons, u_pool_ons);
            else ref_ons.motivo = "sin nulo lineal: modelo invalido";
            rama.nulo_onsager = ref_ons;
            if (ref_ons.evaluado) { u_ons.push_back(ref_ons.u_pit); dt_ons.push_back(ref_ons.delta_tau); }
        }

        fn.contraste_iid = detalle_nulo::evaluar_contraste_resumen(
            u_iid, dt_iid, u_pool_iid, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ 0xABCDEF0123456789ULL);
        fn.contraste_onsager = detalle_nulo::evaluar_contraste_resumen(
            u_ons, dt_ons, u_pool_ons, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ 0x9876543210FEDCBAULL);
        fn.p3b = detalle_nulo::evaluar_curva_p3b<T_MAX>(
            integracion_, fn.segmentos_calma, fn.modelo_lineal, etas_reales, pasos_reales,
            // Code (revision D6): hasta el K de la calma (los rezagos con ACF fuera de su banda), no hasta 10 escrito.
            std::min<std::size_t>(std::max<std::size_t>(1, fn.modelo_lineal.K), cfg_.t_max_efectivo), cfg_.semilla_nulo, cfg_.b_mc_resumen, cfg_.alfa);
    }

