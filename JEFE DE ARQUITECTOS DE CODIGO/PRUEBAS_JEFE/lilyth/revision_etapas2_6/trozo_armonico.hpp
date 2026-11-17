        // 3. Code (revision D1c): las curvas nulas son los nulos EMPAREJADOS de cada evento real en la Etapa 3 (excedencia en (a), o su
        //    via (b)), no un pool filtrado por |dY| >= 0.5 mediana (numero escrito) y completado con eventos fabricados sin detector. El
        //    min de eta_barra real se calcula sobre los MISMOS eventos que tienen nulo (RESTITUIDA o CENSURADA_HORIZONTE), para que las dos
        //    curvas promedien lo mismo.
        const auto& emp_iid = es_purgada ? emp_iid_purgada_ : emp_iid_cruda_;
        const auto& emp_ons = es_purgada ? emp_ons_purgada_ : emp_ons_cruda_;
        auto sets_de = [&](const std::vector<std::vector<EventoNuloLilyth<T_MAX>>>& emp,
                           std::vector<std::size_t>& idx_con_nulo) {
            std::vector<const std::vector<EventoNuloLilyth<T_MAX>>*> sets;
            idx_con_nulo.clear();
            for (std::size_t idx = 0; idx < emp.size() && idx < cerradas_historico_.size(); ++idx) {
                if (emp[idx].empty()) continue;
                sets.push_back(&emp[idx]);
                idx_con_nulo.push_back(idx);
            }
            return sets;
        };
        std::vector<std::size_t> idx_iid, idx_ons;
        const auto sets_iid = sets_de(emp_iid, idx_iid);
        const auto sets_ons = sets_de(emp_ons, idx_ons);
        auto min_eta_barra_de = [&](const std::vector<std::size_t>& idxs) {
            std::vector<double> suma(H, 0.0);
            std::vector<std::size_t> cuenta(H, 0);
            for (std::size_t idx : idxs) {
                const auto& rama = es_purgada ? cerradas_historico_[idx].purgada : cerradas_historico_[idx].cruda;
                const std::size_t pk = std::min<std::size_t>(rama.pasos_observados, H);
                for (std::size_t tau = 2; tau <= pk; ++tau) {
                    if (std::isfinite(rama.eta[tau - 1])) { suma[tau - 1] += rama.eta[tau - 1]; ++cuenta[tau - 1]; }
                }
            }
            double mn = NaN;
            for (std::size_t tau = 2; tau <= H; ++tau) {
                if (cuenta[tau - 1] == 0) continue;
                const double m = suma[tau - 1] / static_cast<double>(cuenta[tau - 1]);
                if (!std::isfinite(mn) || m < mn) mn = m;
            }
            return mn;
        };

        // 4. Contraste Monte Carlo del sobrepaso de la curva media contra ambos nulos (E1.8), emparejado evento a evento.
        arm.sobrepaso.min_eta_barra_real = arm.prony_curva_media.min_eta_desde_2;
        arm.sobrepaso.min_eta_barra_real_emparejado_iid = min_eta_barra_de(idx_iid);
        arm.sobrepaso.min_eta_barra_real_emparejado_onsager = min_eta_barra_de(idx_ons);
        detalle_armonico::contraste_min_eta_barra_emparejado<T_MAX>(
            arm.sobrepaso.min_eta_barra_real_emparejado_iid, sets_iid, H, cfg_.b_mc_resumen,
            cfg_.semilla_nulo ^ 0xA5A5A5A55A5A5A5AULL,
            arm.sobrepaso.mediana_min_eta_barra_iid, arm.sobrepaso.p_sobrepaso_curva_iid);
        detalle_armonico::contraste_min_eta_barra_emparejado<T_MAX>(
            arm.sobrepaso.min_eta_barra_real_emparejado_onsager, sets_ons, H, cfg_.b_mc_resumen,
            cfg_.semilla_nulo ^ 0x5A5A5A5AA5A5A5A5ULL,
            arm.sobrepaso.mediana_min_eta_barra_onsager, arm.sobrepaso.p_sobrepaso_curva_onsager);

        // 5. SOBREPASO_POR_EVENTO: s_k = min de eta_k en [tau_rest_k, T_max], PIT de punto medio contra SU conjunto emparejado; el pool
        //    del Monte Carlo es el PIT de cada nulo contra el resto de su conjunto.
        auto sob_de = [&](const EventoNuloLilyth<T_MAX>& ev) {
            return detalle_armonico::min_eta_post_restitucion<T_MAX>(ev.eta, ev.pasos_observados, ev.tau_rest_consolidado, H);
        };
        auto contraste_sob = [&](const std::vector<const std::vector<EventoNuloLilyth<T_MAX>>*>& sets,
                                 const std::vector<std::size_t>& idxs, std::uint64_t sal) {
            std::vector<double> u_k, d_k, u_pool;
            for (std::size_t q = 0; q < sets.size(); ++q) {
                std::vector<double> sn;
                for (const auto& ev : *sets[q]) { const double v = sob_de(ev); if (std::isfinite(v)) sn.push_back(v); }
                if (sn.size() < 2) continue;
                for (std::size_t i = 0; i < sn.size(); ++i) {
                    std::vector<double> resto;
                    for (std::size_t j = 0; j < sn.size(); ++j) if (j != i) resto.push_back(sn[j]);
                    u_pool.push_back(detalle_nulo::pit_punto_medio(sn[i], resto));
                }
                const auto& rama = es_purgada ? cerradas_historico_[idxs[q]].purgada : cerradas_historico_[idxs[q]].cruda;
                const double s_k = detalle_armonico::min_eta_post_restitucion<T_MAX>(rama.eta, rama.pasos_observados, rama.tau_rest_consolidado, H);
                if (!std::isfinite(s_k)) continue;
                u_k.push_back(detalle_nulo::pit_punto_medio(s_k, sn));
                d_k.push_back(s_k - ark::metrologia::est::cuantil(sn, 0.5));
            }
            return detalle_nulo::evaluar_contraste_resumen(u_k, d_k, u_pool, cfg_.b_mc_resumen, cfg_.semilla_nulo ^ sal);
        };
        arm.sobrepaso.sobrepaso_por_evento_iid = contraste_sob(sets_iid, idx_iid, 0x3333333333333333ULL);
        arm.sobrepaso.sobrepaso_por_evento_onsager = contraste_sob(sets_ons, idx_ons, 0x7777777777777777ULL);

        // 6. Etiqueta ARMONICO_SOSPECHOSO por evento (Prony con >= 4 puntos y raices complejas Y min_eta_desde_2 < Q_alfa(min_eta_nulo)),
        //    con el cuantil sobre los nulos emparejados (cada nulo pesa por el evento real al que esta emparejado).
        auto juntar = [&](const std::vector<const std::vector<EventoNuloLilyth<T_MAX>>*>& sets) {
            std::vector<EventoNuloLilyth<T_MAX>> out;
            for (const auto* s : sets) out.insert(out.end(), s->begin(), s->end());
            return out;
        };
        const std::vector<EventoNuloLilyth<T_MAX>> pool_iid = juntar(sets_iid);
        const std::vector<EventoNuloLilyth<T_MAX>> pool_ons = juntar(sets_ons);
        auto extraer_mins_desde_2 = [&](const std::vector<EventoNuloLilyth<T_MAX>>& pool) {
            std::vector<double> mins;
            mins.reserve(pool.size());
            for (const auto& ev : pool) {
                const double mn = detalle_armonico::min_eta_en_rango<T_MAX>(ev.eta, ev.pasos_observados, 2, H);
                if (std::isfinite(mn)) mins.push_back(mn);
            }
            return mins;
        };
        auto cuantil_orden_inferior = [&](std::vector<double> v, double alfa_q) -> double {
            if (v.empty()) return NaN;
            std::sort(v.begin(), v.end());
            const std::size_t idx_q = std::min(v.size() - 1, static_cast<std::size_t>(std::floor(alfa_q * static_cast<double>(v.size()))));
            return v[idx_q];
        };
        const std::vector<double> mins_iid = extraer_mins_desde_2(pool_iid);
        const std::vector<double> mins_ons = extraer_mins_desde_2(pool_ons);
        arm.sobrepaso.cuantil_alfa_min_eta_iid = cuantil_orden_inferior(mins_iid, cfg_.alfa);
        arm.sobrepaso.cuantil_alfa_min_eta_onsager = cuantil_orden_inferior(mins_ons, cfg_.alfa);

        auto calcular_tasa_etiqueta_pool = [&](const std::vector<EventoNuloLilyth<T_MAX>>& pool, double q_alfa) {
            std::size_t n_eval = 0, n_etiq = 0;
            for (const auto& ev : pool) {
                const std::size_t pk = std::min({ev.pasos_observados, H, T_MAX});
                if (pk < 4) continue;
                const double mn = detalle_armonico::min_eta_en_rango<T_MAX>(ev.eta, pk, 2, H);
                if (!std::isfinite(mn)) continue;
                ++n_eval;
                if (detalle_armonico::evaluar_etiqueta_armonico_evento<T_MAX>(ev.eta, pk, mn, H, q_alfa)) {
                    ++n_etiq;
                }
            }
            return n_eval > 0 ? static_cast<double>(n_etiq) / static_cast<double>(n_eval) : 0.0;
        };
        arm.sobrepaso.tasa_etiqueta_nulo_iid = calcular_tasa_etiqueta_pool(pool_iid, arm.sobrepaso.cuantil_alfa_min_eta_iid);
        arm.sobrepaso.tasa_etiqueta_nulo_onsager = calcular_tasa_etiqueta_pool(pool_ons, arm.sobrepaso.cuantil_alfa_min_eta_onsager);

