// adaptadores_fuentes.hpp -- PORTE TAL CUAL DE LOS ADAPTADORES PYTHON DE ARQUITECTO DE DATOS (ARK-DATA-02-CXX)
// Reactor ARK / Cancion del Micelio -- Jefe de Arquitectos de Codigo, 2026-11-02
//
// Calca, funcion por funcion y en el mismo orden de operaciones, de:
//   recopilacion_binance/adaptador_binance.py     cargar_binance_um_1min, cargar_binance_parquet, cargar_binance_diario_corte_ny,
//                                                 cargar_datos_binance
//   recopilacion_bitfinex/adaptador_bitfinex.py   normalizar_simbolo_bitfinex, cargar_bitfinex_1s, ordenar_eventos_l3_causal,
//                                                 stream_ordenes_l3_chunks, ReconstructorLibroL3, cargar_datos_bitfinex
//   recopilacion_bloomberg/adaptador_bloomberg.py cargar_bloomberg_1min, cargar_bloomberg_diario, cargar_datos_bloomberg
//   recopilacion_unificada/adaptador_yahoo.py     cargar_datos_yahoo_disco (con _migrar_a_cierre_crudo)
// Las reglas, umbrales y estampas son las del Python (decisiones ya tomadas; aqui no se cambia ninguna). La identidad celda a
// celda la comprueba tests/contraste_adaptadores.py contra el Python real.
// Soporte sin equivalente en el adaptador Python (que recibe los registros ya leidos): leer_jsonl_gz_l3, para alimentar el
// flujo L3 desde las capturas *.jsonl.gz (nunca abre un .tmp).
#pragma once
#include "ark/dataset/formato_parquet.hpp"
#include "ark/dataset/tabla_ingesta.hpp"
#include <functional>
#include <numeric>
#include <optional>
#include <tuple>
#include <unordered_map>

namespace ark::dataset::ingesta {

namespace fmt = ark::dataset::formatos;

// ===================================================================================================================== BINANCE
inline std::string dir_cache_cuadrivector() { return unir(raiz_datos(), "cachas/cuadrivector_oos"); }
inline std::string dir_cache_vrp_dvol() { return unir(raiz_datos(), "cachas/vrp_dvol"); }
inline std::string dir_data_parquet() { return unir(raiz_datos(), "data"); }

inline std::string limpiar_simbolo(const std::string& s) {
    return mayusculas(reemplazar(reemplazar(reemplazar(s, "/", ""), "-", ""), "_", ""));
}

inline std::vector<std::string> encontrar_ficheros_npz_simbolo(const std::string& simbolo) {
    const std::string clean = limpiar_simbolo(simbolo);
    std::vector<std::string> f;
    const std::string p_vrp = unir(dir_cache_vrp_dvol(), "k1m_ext_" + clean + ".npz");
    if (es_fichero(p_vrp)) f.push_back(p_vrp);
    for (const auto& x : glob_simple(dir_cache_cuadrivector(), "k1m_oos_" + clean + "_", ".npz")) f.push_back(x);
    return f;
}

inline Tabla cargar_binance_um_1min(const std::string& simbolo, const std::string& fecha_inicio = "", const std::string& fecha_fin = "",
                                    bool estampa_cierre = true, bool incluir_flujo = true) {
    const auto ficheros = encontrar_ficheros_npz_simbolo(simbolo);
    if (ficheros.empty())
        throw SimboloNoEncontradoError("INCOMPATIBILIDAD DE CONTRATO (Venue=BINANCE_UM, Símbolo='" + simbolo + "'): No se encontraron cachés .npz en "
                                       + dir_cache_cuadrivector() + " ni en " + dir_cache_vrp_dvol() + ".");
    std::vector<int64_t> t_total;
    std::vector<std::array<double, 10>> datos_total;
    std::vector<int64_t> t_falt;
    for (const auto& f : ficheros) {
        const auto d = fmt::cargar_npz(f);
        auto it_t = d.find("t"), it_d = d.find("datos");
        if (it_t == d.end() || it_d == d.end()) throw ValueError("KeyError: npz sin 't' o 'datos' en " + f);
        std::vector<int64_t> t_sec = it_t->second.a_int64();
        for (auto& x : t_sec) x *= 60;
        const std::vector<double> dat = it_d->second.a_double();
        const size_t n = t_sec.size(), nc = it_d->second.forma.size() > 1 ? it_d->second.forma[1] : 1;
        if (nc < 10 || dat.size() != n * nc) throw ValueError("datos de Binance con forma inesperada en " + f);
        for (size_t k = 1; k < n; ++k)
            if (t_sec[k] - t_sec[k - 1] <= 0)
                throw ErrorMonotoniaTiempo("VIOLACIÓN DE MONOTONÍA TEMPORAL EN TRAMO BINANCE (" + f
                                           + "): Se detectaron timestamps invertidos o duplicados dentro del archivo .npz.");
        auto it_f = d.find("faltan");
        if (it_f != d.end() && it_f->second.elementos() > 0) {
            for (int64_t x : it_f->second.a_int64())
                if (x >= 0 && x < static_cast<int64_t>(n)) t_falt.push_back(t_sec[static_cast<size_t>(x)]);
        }
        for (size_t k = 0; k < n; ++k) {
            std::array<double, 10> fila{};
            for (size_t c = 0; c < 10; ++c) fila[c] = dat[k * nc + c];
            t_total.push_back(t_sec[k]);
            datos_total.push_back(fila);
        }
    }
    // np.unique(return_index=True): orden por valor, primera aparicion (argsort estable).
    std::vector<size_t> idx(t_total.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return t_total[a] < t_total[b]; });
    std::vector<int64_t> t;
    std::vector<std::array<double, 10>> D;
    for (size_t k = 0; k < idx.size(); ++k) {
        if (k > 0 && t_total[idx[k]] == t_total[idx[k - 1]]) continue;
        t.push_back(t_total[idx[k]]);
        D.push_back(datos_total[idx[k]]);
    }
    std::sort(t_falt.begin(), t_falt.end());
    t_falt.erase(std::unique(t_falt.begin(), t_falt.end()), t_falt.end());
    std::vector<int64_t> imp(t.size());
    for (size_t k = 0; k < t.size(); ++k) imp[k] = std::binary_search(t_falt.begin(), t_falt.end(), t[k]) ? 1 : 0;
    if (estampa_cierre) for (auto& x : t) x += 60;

    Tabla df;
    std::vector<int64_t> ts(t.size());
    for (size_t k = 0; k < t.size(); ++k) ts[k] = t[k] * NS_S;
    df.poner_ts("fecha_utc", ts);
    const char* nombres[10] = {"open", "high", "low", "close", "volume", "close_time_ms", "quote_volume", "trades_count", "taker_buy_base", "taker_buy_quote"};
    for (size_t c = 0; c < 10; ++c) {
        std::vector<double> v(D.size());
        for (size_t k = 0; k < D.size(); ++k) v[k] = D[k][c];
        df.poner_f64(nombres[c], v);
    }
    df.poner_f64("precio", df.col("close").f);
    df.poner_i64("es_imputado", imp);
    if (incluir_flujo) {
        const auto& v = df.col("volume").f;
        const auto& tbv = df.col("taker_buy_base").f;
        const auto& qv = df.col("quote_volume").f;
        const auto& n = df.col("trades_count").f;
        std::vector<double> ofi(v.size()), tick(v.size()), tbr(v.size());
        for (size_t k = 0; k < v.size(); ++k) {
            const double dv = v[k] > 0 ? v[k] : 1.0;
            ofi[k] = v[k] > 0 ? (2.0 * tbv[k] - v[k]) / dv : 0.0;
            const double dn = n[k] > 0 ? n[k] : 1.0;
            tick[k] = n[k] > 0 ? qv[k] / dn : 0.0;
            tbr[k] = v[k] > 0 ? tbv[k] / dv : 0.5;
        }
        df.poner_f64("ofi_ratio", ofi);
        df.poner_f64("ticket_promedio_usd", tick);
        df.poner_f64("taker_buy_ratio", tbr);
    }
    filtrar_rango(df, "fecha_utc", fecha_inicio, fecha_fin);
    return df;
}

inline int64_t segundos_resolucion(const std::string& r) {
    static const std::map<std::string, int64_t> M = {{"1m", 60}, {"5m", 300}, {"15m", 900}, {"30m", 1800}, {"1h", 3600}, {"4h", 14400}, {"1d", 86400}};
    auto it = M.find(r);
    return it == M.end() ? 900 : it->second;
}

// Columna de parquet -> columna de la tabla (DOUBLE/FLOAT -> F64; TIMESTAMP -> TS en ns; enteros -> I64).
inline Columna columna_desde_parquet(const fmt::ColumnaParquet& c) {
    Columna o;
    o.nombre = c.nombre;
    if (c.tipo == fmt::TipoFisico::DOUBLE || c.tipo == fmt::TipoFisico::FLOAT) { o.tipo = TipoColumna::F64; o.f = c.dobles; return o; }
    if (!c.unidad_timestamp.empty()) {
        const int64_t m = c.unidad_timestamp == "ns" ? 1 : (c.unidad_timestamp == "us" ? 1000 : 1000000);
        o.tipo = TipoColumna::TS;
        o.i.resize(c.enteros.size());
        for (size_t k = 0; k < c.enteros.size(); ++k) o.i[k] = c.nulo[k] ? NAT : c.enteros[k] * m;
        return o;
    }
    o.tipo = TipoColumna::I64;
    o.i = c.enteros;
    return o;
}

// pd.to_datetime(col, utc=True): TS se queda; numerica se lee como ns (unidad por omision de pandas).
inline std::vector<int64_t> a_datetime_utc(const Columna& c) {
    if (c.tipo == TipoColumna::TS || c.tipo == TipoColumna::I64) return c.i;
    if (c.tipo == TipoColumna::F64) { std::vector<int64_t> v; for (double x : c.f) v.push_back(std::isfinite(x) ? static_cast<int64_t>(x) : NAT); return v; }
    std::vector<int64_t> v;
    for (const auto& s : c.s) v.push_back(ns_desde_iso(s));
    return v;
}

inline std::vector<int64_t> imputado_por_precios(const Tabla& df, const std::vector<std::string>& candidatos) {
    std::vector<int64_t> imp(df.filas(), 0);
    for (const auto& n : candidatos) {
        if (!df.tiene(n)) continue;
        const auto& v = df.col(n).f;
        for (size_t k = 0; k < v.size(); ++k) if (std::isnan(v[k]) || std::isinf(v[k]) || v[k] <= 0.0) imp[k] = 1;
    }
    return imp;
}

inline Tabla cargar_binance_parquet(const std::string& simbolo, const std::string& resolucion = "15m", const std::string& fecha_inicio = "",
                                    const std::string& fecha_fin = "", bool estampa_cierre = true) {
    const std::string clean = limpiar_simbolo(simbolo);
    std::string nombre;
    if (clean.size() >= 4 && clean.compare(clean.size() - 4, 4, "USDT") == 0) nombre = clean.substr(0, clean.size() - 4) + "_USDT_" + resolucion + ".parquet";
    else nombre = clean + "_" + resolucion + ".parquet";
    const std::string ruta = unir(dir_data_parquet(), nombre);
    if (!es_fichero(ruta))
        throw SimboloNoEncontradoError("INCOMPATIBILIDAD DE CONTRATO (Venue=BINANCE_UM, Símbolo='" + simbolo + "', Resolución='" + resolucion
                                       + "'): No se encontró el fichero parquet en " + ruta + ".");
    const fmt::TablaParquet tp = fmt::leer_parquet(ruta);
    // Indice de pandas: metadatos 'pandas' -> index_columns (nombres de columnas que forman el indice; un dict = RangeIndex).
    std::vector<std::string> indices;
    for (const auto& [k, v] : tp.metadatos) {
        if (k != "pandas") continue;
        const JSON j = parsear_json(v);
        if (const JSON* ic = j.get("index_columns")) for (const auto& e : ic->a) if (e.t == JSON::T::STR) indices.push_back(e.s);
    }
    Tabla df;
    std::optional<Columna> indice;
    for (const auto& c : tp.columnas) {
        const bool es_indice = std::find(indices.begin(), indices.end(), c.nombre) != indices.end();
        if (es_indice && !indice) { indice = columna_desde_parquet(c); continue; }
        if (es_indice) continue;
        df.columnas.push_back(columna_desde_parquet(c));
    }
    if (indice) {
        const std::string nl = minusculas(indice->nombre);
        const bool fecha = indice->tipo == TipoColumna::TS || nl == "timestamp" || nl == "time" || nl == "datetime" || nl == "date" || nl == "open_time";
        if (fecha) df.poner_ts("open_time", a_datetime_utc(*indice));
    }
    for (auto& c : df.columnas) c.nombre = minusculas(c.nombre);
    if (df.tiene("close_time") && estampa_cierre) {
        df.poner_ts("fecha_utc", a_datetime_utc(df.col("close_time")));
    } else {
        std::string col_fecha;
        for (const char* c : {"open_time", "timestamp", "time", "datetime", "date"}) if (df.tiene(c)) { col_fecha = c; break; }
        if (col_fecha.empty()) throw ValueError("No se encontró columna o índice temporal en parquet " + ruta + ".");
        std::vector<int64_t> ap = a_datetime_utc(df.col(col_fecha));
        if (estampa_cierre) { const int64_t d = segundos_resolucion(resolucion) * NS_S; for (auto& x : ap) if (x != NAT) x += d; }
        df.poner_ts("fecha_utc", ap);
    }
    if (resolucion == "1d" && estampa_cierre) {
        std::vector<int64_t> ses = df.col("fecha_utc").i, fu(ses.size());
        for (size_t k = 0; k < ses.size(); ++k) { ses[k] = normalizar_dia(ses[k] - 1000000LL); fu[k] = ses[k] + NS_DIA; }
        df.poner_ts("fecha_sesion", ses);
        df.poner_ts("fecha_utc", fu);
    }
    df.poner_f64("precio", df.col("close").f);
    df.poner_i64("es_imputado", imputado_por_precios(df, {"open", "high", "low", "close", "precio"}));
    filtrar_rango(df, "fecha_utc", fecha_inicio, fecha_fin);
    auto [mono, inc] = verificar_monotonicidad(df.col("fecha_utc").i);
    if (!mono)
        throw ErrorMonotoniaTiempo("VIOLACIÓN DE MONOTONÍA TEMPORAL EN PARQUET BINANCE (" + simbolo + "): "
                                   + (inc.empty() ? std::string("Índice temporal no monótono") : inc.front().descripcion) + ".");
    return df;
}

inline constexpr int HORA_CORTE_NY = 16;

inline Tabla cargar_binance_diario_corte_ny(const std::string& simbolo, const std::string& fecha_inicio = "", const std::string& fecha_fin = "") {
    std::string ini_min;
    if (!fecha_inicio.empty()) ini_min = texto_timestamp(ns_desde_iso(fecha_inicio) - 3 * NS_DIA);
    Tabla m = cargar_binance_um_1min(simbolo, ini_min, fecha_fin, true, false);
    if (m.vacia()) return m;
    const auto& fu = m.col("fecha_utc").i;
    const int64_t H16 = static_cast<int64_t>(HORA_CORTE_NY) * 3600 * NS_S;
    // clave = (fin_local - 16 h - 1 us).tz_localize(None).normalize() + 1 dia (resta absoluta, luego hora de pared de NY)
    std::vector<int64_t> clave(fu.size());
    for (size_t k = 0; k < fu.size(); ++k) {
        const int64_t x = fu[k] - H16 - 1000;
        const int64_t pared = x + offset_ny_desde_utc_s(piso_div(x, NS_S)) * NS_S;
        clave[k] = normalizar_dia(pared) + NS_DIA;
    }
    // groupby(clave, sort=True): first/last saltan NaN, max/min/sum ignoran NaN, mean del flag.
    std::vector<size_t> ord(fu.size());
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return clave[a] < clave[b]; });
    const auto& o = m.col("open").f; const auto& h = m.col("high").f; const auto& l = m.col("low").f; const auto& c = m.col("close").f;
    const auto& v = m.col("volume").f; const auto& qv = m.col("quote_volume").f; const auto& n = m.col("trades_count").f;
    const auto& tb = m.col("taker_buy_base").f; const auto& tq = m.col("taker_buy_quote").f; const auto& imp = m.col("es_imputado").i;
    std::vector<int64_t> dia, ses, disp, impd;
    std::vector<double> O, Hh, L, C, V, QV, N, TB, TQ, FR;
    size_t a = 0;
    while (a < ord.size()) {
        size_t b = a;
        while (b < ord.size() && clave[ord[b]] == clave[ord[a]]) ++b;
        double fo = NaN, fc = NaN, mx = NaN, mn = NaN;
        SumaKahanPandas sv, sq, sn, stb, stq, simp;   // groupby sum / mean de pandas: suma compensada
        int64_t ultimo = NAT;
        for (size_t k = a; k < b; ++k) {
            const size_t r = ord[k];
            if (std::isnan(fo) && !std::isnan(o[r])) fo = o[r];
            if (!std::isnan(c[r])) fc = c[r];
            if (!std::isnan(h[r]) && (std::isnan(mx) || h[r] > mx)) mx = h[r];
            if (!std::isnan(l[r]) && (std::isnan(mn) || l[r] < mn)) mn = l[r];
            sv.anadir(v[r]); sq.anadir(qv[r]); sn.anadir(n[r]); stb.anadir(tb[r]); stq.anadir(tq[r]);
            simp.anadir(static_cast<double>(imp[r]));
            if (ultimo == NAT || fu[r] > ultimo) ultimo = fu[r];
        }
        const int64_t key = clave[ord[a]];
        // corte = tz_localize(dia, NY) + 16 h (suma absoluta)
        const int64_t pared_s = piso_div(key, NS_S);
        const int64_t corte = (pared_s - offset_s(Zona::NUEVA_YORK, pared_s)) * NS_S + H16;
        if (ultimo == corte) {
            dia.push_back(key);
            O.push_back(fo); Hh.push_back(mx); L.push_back(mn); C.push_back(fc);
            V.push_back(sv.sumx); QV.push_back(sq.sumx); N.push_back(sn.sumx); TB.push_back(stb.sumx); TQ.push_back(stq.sumx);
            const double fr = simp.sumx / static_cast<double>(simp.nobs);
            FR.push_back(fr);
            impd.push_back(fr >= 1.0 ? 1 : 0);
        }
        a = b;
    }
    Tabla df;
    for (int64_t d : dia) { ses.push_back(d); disp.push_back(d + NS_DIA); }
    df.poner_ts("fecha_utc", disp);
    df.poner_ts("fecha_sesion", ses);
    df.poner_f64("open", O); df.poner_f64("high", Hh); df.poner_f64("low", L); df.poner_f64("close", C);
    df.poner_f64("volume", V); df.poner_f64("quote_volume", QV); df.poner_f64("trades_count", N);
    df.poner_f64("taker_buy_base", TB); df.poner_f64("taker_buy_quote", TQ);
    df.poner_f64("precio", C);
    df.poner_f64("frac_minutos_imputados", FR);
    df.poner_i64("es_imputado", impd);
    df.poner_str("corte_diario", std::vector<std::string>(dia.size(), "NY_16"));
    filtrar_rango(df, "fecha_utc", fecha_inicio, fecha_fin);
    return df;
}

inline Tabla cargar_datos_binance(const std::string& simbolo, const std::string& fecha_inicio, const std::string& fecha_fin,
                                  const std::string& resolucion = "1m", bool estampa_cierre = true) {
    if (resolucion == "1m" || resolucion == "intradia" || resolucion == "minuto")
        return cargar_binance_um_1min(simbolo, fecha_inicio, fecha_fin, estampa_cierre, true);
    if (resolucion == "1d_ny16") return cargar_binance_diario_corte_ny(simbolo, fecha_inicio, fecha_fin);
    if (resolucion == "15m" || resolucion == "1h" || resolucion == "4h" || resolucion == "1d")
        return cargar_binance_parquet(simbolo, resolucion, fecha_inicio, fecha_fin, estampa_cierre);
    throw ValueError("Resolución de Binance '" + resolucion + "' no soportada (use '1m', '15m', '1h', '4h', '1d').");
}

// ==================================================================================================================== BITFINEX
inline std::string dir_cache_bitfinex_1s() { return unir(raiz_datos(), "cachas/bitfinex_1s"); }

inline std::string normalizar_simbolo_bitfinex(const std::string& simbolo) {
    std::string s = recortar(simbolo);
    if (s.empty() || s[0] != 't') s = "t" + s;
    if (s.size() >= 4 && s.compare(s.size() - 4, 4, "USDT") == 0) s = s.substr(0, s.size() - 4) + "USD";
    return s;
}

inline Tabla cargar_bitfinex_1s(const std::string& simbolo, const std::string& ventana = "W1", const std::string& fecha_inicio = "",
                                const std::string& fecha_fin = "") {
    const std::string sym = normalizar_simbolo_bitfinex(simbolo);
    std::string ruta = unir(dir_cache_bitfinex_1s(), ventana + "_" + sym + "_1s.npz");
    if (!es_fichero(ruta)) {
        const std::string alt = unir(dir_cache_bitfinex_1s(), "W1P_" + sym + "_1s.npz");
        if (es_fichero(alt)) ruta = alt;
        else throw SimboloNoEncontradoError("INCOMPATIBILIDAD DE CONTRATO (Venue=BITFINEX, Símbolo='" + sym
                                            + "'): No se encontró la caché de 1s de Bitfinex en " + dir_cache_bitfinex_1s() + ".");
    }
    const auto d = fmt::cargar_npz(ruta);
    auto arr = [&](const char* k) -> const fmt::ArregloNpy& {
        auto it = d.find(k);
        if (it == d.end()) throw ValueError(std::string("KeyError: ") + k);
        return it->second;
    };
    const std::vector<int64_t> ts_ms = arr("ts_1s").a_int64();
    const std::vector<double> bb = arr("bb_1s").a_double(), ba = arr("ba_1s").a_double();
    const std::vector<uint8_t> cob = arr("v_1s").a_bool();
    const std::vector<double> lfin = arr("lfin_1s").a_double();
    for (size_t k = 1; k < ts_ms.size(); ++k)
        if (ts_ms[k] - ts_ms[k - 1] <= 0)
            throw ErrorMonotoniaTiempo("VIOLACIÓN DE MONOTONÍA TEMPORAL EN CACHÉ BITFINEX 1S (" + sym
                                       + "): Se detectaron timestamps invertidos o duplicados en la fuente cruda.");
    const size_t n = ts_ms.size();
    std::vector<double> mid(n), spread(n);
    std::vector<int64_t> imp(n), ts(n), cobv(n);
    for (size_t k = 0; k < n; ++k) {
        mid[k] = (bb[k] + ba[k]) / 2.0;
        spread[k] = mid[k] > 0 ? ((ba[k] - bb[k]) / mid[k]) * 10000.0 : 0.0;
        const bool corrupto = !std::isfinite(bb[k]) || !std::isfinite(ba[k]) || !std::isfinite(mid[k]) || !std::isfinite(lfin[k])
                              || bb[k] <= 0.0 || ba[k] <= 0.0 || mid[k] <= 0.0 || ba[k] < bb[k] || !cob[k];
        imp[k] = corrupto ? 1 : 0;
        ts[k] = ts_ms[k] * 1000000LL;
        cobv[k] = cob[k];
    }
    Tabla df;
    df.poner_ts("fecha_utc", ts);
    df.poner_f64("best_bid", bb);
    df.poner_f64("best_ask", ba);
    df.poner_f64("mid_price", mid);
    df.poner_f64("precio", mid);
    df.poner_bool("cobertura_valida", cobv);
    df.poner_f64("spread_pb", spread);
    df.poner_i64("es_imputado", imp);
    filtrar_rango(df, "fecha_utc", fecha_inicio, fecha_fin);
    return df;
}

// Registro L3 (un dict de Python): se conserva el objeto JSON entero y se exponen los campos que usa el adaptador.
struct RegistroL3 {
    JSON obj;
    bool tiene_s() const { const JSON* v = obj.get("s"); return v && v->t != JSON::T::NULO; }
    int64_t s() const { return obj.get("s")->i; }
    int64_t cx() const { const JSON* v = obj.get("cx"); return v ? v->i : 0; }   // r.get("cx", 0)
    int64_t cen() const { const JSON* v = obj.get("cen"); return v ? v->i : 0; }
    std::tuple<int64_t, int64_t> clave() const { return {cx(), s()}; }
};

inline std::pair<std::vector<RegistroL3>, int64_t> ordenar_eventos_l3_causal(const std::vector<RegistroL3>& registros) {
    std::vector<RegistroL3> con;
    int64_t sin = 0;
    for (const auto& r : registros) { if (r.tiene_s()) con.push_back(r); else ++sin; }
    std::stable_sort(con.begin(), con.end(), [](const RegistroL3& a, const RegistroL3& b) { return a.clave() < b.clave(); });
    return {con, sin};
}

struct ContadorL3 {
    int64_t sin_seq = 0, censuradas = 0;
};

// stream_ordenes_l3_chunks: bloques de chunk_size registros, ordenados por (cx, s) dentro del bloque, con la guarda de orden
// global entre bloques y el veto a censuradas. 'emitir' recibe cada bloque limpio no vacio (el yield del generador).
template <class Iterable>
inline void stream_ordenes_l3_chunks(const Iterable& flujo, const std::function<void(std::vector<RegistroL3>&)>& emitir,
                                     size_t chunk_size = 5000, bool descartar_censuradas = true, ContadorL3* contador = nullptr) {
    if (chunk_size < 1) throw ValueError("chunk_size=" + std::to_string(chunk_size) + " debe ser >= 1.");
    ContadorL3 local;
    ContadorL3& cuenta = contador ? *contador : local;
    std::optional<std::tuple<int64_t, int64_t>> ultima;
    std::vector<RegistroL3> buf;
    auto procesar = [&](std::vector<RegistroL3>& b) {
        auto [ord, sin] = ordenar_eventos_l3_causal(b);
        cuenta.sin_seq += sin;
        if (!ord.empty()) {
            const auto primera = ord.front().clave();
            if (ultima && primera <= *ultima)
                throw ErrorMonotoniaTiempo("ORDEN L3 NO GLOBAL: el bloque empieza en (cx, s) = (" + std::to_string(std::get<0>(primera)) + ", "
                                           + std::to_string(std::get<1>(primera)) + "), no posterior al último emitido ("
                                           + std::to_string(std::get<0>(*ultima)) + ", " + std::to_string(std::get<1>(*ultima))
                                           + "). Un evento llegó fuera de orden por encima del tamaño de bloque.");
            ultima = ord.back().clave();
        }
        if (descartar_censuradas) {
            std::vector<RegistroL3> val;
            for (auto& r : ord) if (r.cen() != 1) val.push_back(std::move(r));
            cuenta.censuradas += static_cast<int64_t>(ord.size() - val.size());
            return val;
        }
        return ord;
    };
    for (const auto& r : flujo) {
        buf.push_back(r);
        if (buf.size() >= chunk_size) {
            std::vector<RegistroL3> limpio = procesar(buf);
            buf.clear();
            if (!limpio.empty()) emitir(limpio);
        }
    }
    if (!buf.empty()) {
        std::vector<RegistroL3> limpio = procesar(buf);
        if (!limpio.empty()) emitir(limpio);
    }
}

class ReconstructorLibroL3 {
public:
    std::unordered_map<int64_t, std::tuple<double, double, double>> bids, asks;   // order_id -> (precio, cantidad, t)

    void procesar_evento(const RegistroL3& r) {
        auto num = [&](const char* k) { const JSON* v = r.obj.get(k); return v ? v->n : 0.0; };
        const JSON* kk = r.obj.get("k");
        const std::string k = (kk && kk->t == JSON::T::STR) ? kk->s : std::string();
        const JSON* oo = r.obj.get("o");
        const double p = num("p"), c = num("c"), t = num("t");
        if (k == "snap") {
            bids.clear();
            asks.clear();
            if (oo && oo->t == JSON::T::ARR) {
                for (const auto& orden : oo->a) {
                    const int64_t oid = orden.a.at(0).entero ? orden.a.at(0).i : static_cast<int64_t>(orden.a.at(0).n);
                    const double px = orden.a.at(1).n, cant = orden.a.at(2).n;
                    if (cant > 0) bids[oid] = {px, cant, t};
                    else if (cant < 0) asks[oid] = {px, cant, t};
                }
            }
        } else if (k == "a" || k == "m") {
            const int64_t oid = oo ? oo->i : 0;
            if (c > 0) { bids[oid] = {p, c, t}; asks.erase(oid); }
            else if (c < 0) { asks[oid] = {p, c, t}; bids.erase(oid); }
        } else if (k == "b") {
            const int64_t oid = oo ? oo->i : 0;
            bids.erase(oid);
            asks.erase(oid);
        }
    }
    std::pair<double, double> mejor_bid_ask() const {
        double mb = NaN, ma = NaN;
        for (const auto& [id, v] : bids) if (std::isnan(mb) || std::get<0>(v) > mb) mb = std::get<0>(v);
        for (const auto& [id, v] : asks) if (std::isnan(ma) || std::get<0>(v) < ma) ma = std::get<0>(v);
        return {mb, ma};
    }
};

// Soporte: registros de una captura *.jsonl.gz (una linea = un dict). Se niega a abrir un .tmp.
inline std::vector<RegistroL3> leer_jsonl_gz_l3(const std::string& ruta) {
    if (ruta.size() >= 4 && ruta.compare(ruta.size() - 4, 4, ".tmp") == 0) throw ValueError("prohibido abrir un .tmp: " + ruta);
    const std::vector<uint8_t> crudo = fmt::desgzip(fmt::leer_fichero(ruta));
    std::vector<RegistroL3> v;
    size_t a = 0;
    for (size_t b = 0; b <= crudo.size(); ++b) {
        if (b == crudo.size() || crudo[b] == '\n') {
            if (b > a) {
                const std::string lin(reinterpret_cast<const char*>(&crudo[a]), b - a);
                if (!recortar(lin).empty()) v.push_back(RegistroL3{parsear_json(lin)});
            }
            a = b + 1;
        }
    }
    return v;
}

enum class NivelDato : uint8_t { L1_OHLCV, L1_ENRIQUECIDO, L2_DEPTH, L3_MBO, TOPS };

inline Tabla cargar_datos_bitfinex(const std::string& simbolo, const std::string& fecha_inicio, const std::string& fecha_fin,
                                   NivelDato nivel = NivelDato::L1_OHLCV, const std::string& ventana = "W1") {
    if (nivel == NivelDato::L1_OHLCV || nivel == NivelDato::L2_DEPTH) return cargar_bitfinex_1s(simbolo, ventana, fecha_inicio, fecha_fin);
    if (nivel == NivelDato::L3_MBO)
        throw IncompatibilidadContratoError("INCOMPATIBILIDAD DE MODO (Venue=BITFINEX, Nivel=L3_MBO): Para L3_MBO está prohibida la carga tabular en RAM "
                                            "(Axioma 5); utilice directamente el generador 'stream_ordenes_l3_chunks' sobre la fuente de eventos.");
    throw IncompatibilidadContratoError("Nivel de datos no soportado en Bitfinex.");
}

// =================================================================================================================== BLOOMBERG
inline std::string dir_cache_bbg_1min() { return unir(raiz_datos(), "cachas/bbg_1min"); }
inline std::string dir_cache_macro_d1() { return unir(raiz_datos(), "cachas/macro_d1"); }
inline std::string ruta_macro_btc_diario() { return unir(raiz_datos(), "Bloomberg y modelo Macro/bloomberg_macro_btc_diario.csv"); }

inline std::string limpiar_ticker_bbg(const std::string& t) {
    return recortar(reemplazar(reemplazar(reemplazar(t, " Index", ""), " Comdty", ""), " Curncy", ""));
}

inline Tabla cargar_bloomberg_1min(const std::string& ticker, const std::string& fecha_inicio = "", const std::string& fecha_fin = "",
                                   bool estampa_cierre = true) {
    const std::string clean = limpiar_ticker_bbg(ticker);
    const std::string ruta = unir(dir_cache_bbg_1min(), clean + ".npz");
    if (!es_fichero(ruta))
        throw SimboloNoEncontradoError("INCOMPATIBILIDAD DE CONTRATO (Venue=BLOOMBERG_XLSX, Símbolo='" + ticker
                                       + "'): No se encontró el fichero intradía de Bloomberg en " + ruta + ".");
    const auto d = fmt::cargar_npz(ruta);
    auto it_t = d.find("t"), it_o = d.find("ohlcv");
    if (it_t == d.end() || it_o == d.end()) throw ValueError("KeyError: npz de Bloomberg sin 't' u 'ohlcv'");
    std::vector<int64_t> t = it_t->second.a_int64();
    const std::vector<double> ohlcv = it_o->second.a_double();
    const size_t n = t.size(), nc = it_o->second.forma.size() > 1 ? it_o->second.forma[1] : 1;
    if (nc < 5 || ohlcv.size() != n * nc) throw ValueError("ohlcv de Bloomberg con forma inesperada");
    for (size_t k = 1; k < n; ++k)
        if (t[k] - t[k - 1] <= 0)
            throw ErrorMonotoniaTiempo("VIOLACIÓN DE MONOTONÍA TEMPORAL EN CACHÉ BLOOMBERG 1M (" + ticker
                                       + "): Se detectaron timestamps invertidos o duplicados en la fuente cruda.");
    if (estampa_cierre) for (auto& x : t) x += 60;
    std::vector<int64_t> ts(n), imp(n);
    std::vector<double> col[5];
    for (auto& c : col) c.resize(n);
    for (size_t k = 0; k < n; ++k) {
        ts[k] = t[k] * NS_S;
        bool mal = false;
        for (size_t c = 0; c < 5; ++c) col[c][k] = ohlcv[k * nc + c];
        for (size_t c = 0; c < 4; ++c) { const double x = col[c][k]; if (std::isnan(x) || std::isinf(x) || x <= 0.0) mal = true; }
        imp[k] = mal ? 1 : 0;
    }
    Tabla df;
    df.poner_ts("fecha_utc", ts);
    df.poner_f64("open", col[0]); df.poner_f64("high", col[1]); df.poner_f64("low", col[2]); df.poner_f64("close", col[3]);
    df.poner_f64("volume", col[4]);
    df.poner_f64("precio", col[3]);
    df.poner_i64("es_imputado", imp);
    filtrar_rango(df, "fecha_utc", fecha_inicio, fecha_fin);
    return df;
}

inline std::string simbolo_fuente_macro_d1(const std::string& ticker) {
    try {
        std::ifstream f(unir(dir_cache_macro_d1(), "MANIFIESTO.json"), std::ios::binary);
        if (!f) return "DESCONOCIDO";
        std::stringstream ss;
        ss << f.rdbuf();
        const JSON man = parsear_json(ss.str());
        const JSON* e = man.get(ticker);
        if (!e || e->t != JSON::T::OBJ) return "DESCONOCIDO";
        const JSON* s = e->get("simbolo");
        if (!s) return "DESCONOCIDO";
        return s->t == JSON::T::STR ? s->s : "DESCONOCIDO";
    } catch (const std::exception&) {
        return "DESCONOCIDO";
    }
}

inline Tabla cargar_bloomberg_diario(const std::string& ticker, const std::string& fecha_inicio = "", const std::string& fecha_fin = "") {
    const std::string clean = limpiar_ticker_bbg(ticker);
    const std::string ruta_csv = unir(dir_cache_macro_d1(), clean + ".csv");
    Tabla df;
    if (es_fichero(ruta_csv)) {
        const CSV raw = leer_csv(ruta_csv);
        auto [ses, disp] = estampar_diario_disponible(raw.texto("fecha"));
        df.poner_ts("fecha_utc", disp);
        df.poner_ts("fecha_sesion", ses);
        std::vector<std::string> presentes;
        for (const char* c : {"open", "high", "low", "close", "volume"}) if (raw.indice(c) >= 0) { presentes.push_back(c); df.poner_f64(c, raw.numerica(c)); }
        if (raw.indice("close") >= 0) df.poner_f64("precio", df.col("close").f);
        else df.poner_f64("precio", raw.numerica("valor"));
        df.poner_str("fuente", std::vector<std::string>(df.filas(), "YAHOO_VIA_MACRO_D1"));
        df.poner_str("simbolo_fuente", std::vector<std::string>(df.filas(), simbolo_fuente_macro_d1(clean)));
    } else if (es_fichero(ruta_macro_btc_diario())) {
        const CSV raw = leer_csv(ruta_macro_btc_diario());
        const std::string col_last = clean + "_PX_LAST";
        if (raw.indice(col_last) < 0)
            throw SimboloNoEncontradoError("INCOMPATIBILIDAD DE CONTRATO (Venue=BLOOMBERG_XLSX, Símbolo='" + ticker
                                           + "'): No se encontró el ticker en cachas/macro_d1 ni en " + ruta_macro_btc_diario() + ".");
        auto [ses, disp] = estampar_diario_disponible(raw.texto("fecha_sesion"));
        std::vector<int64_t> fu;
        for (const auto& s : raw.texto("datetime")) fu.push_back(ns_desde_iso(recortar(s)));
        for (size_t k = 0; k < fu.size(); ++k) fu[k] = std::max(fu[k], disp[k]);
        const std::vector<double> precio = raw.numerica(col_last);
        auto o_def = [&](const std::string& c) { return raw.indice(c) >= 0 ? raw.numerica(c) : precio; };
        if (raw.indice(clean + "_PX_VOLUME") < 0)
            throw ValueError("AttributeError: 'float' object has no attribute 'fillna' (pd.to_numeric(0.0).fillna en el Python)");
        std::vector<double> vol = raw.numerica(clean + "_PX_VOLUME");
        for (auto& x : vol) if (std::isnan(x)) x = 0.0;
        df.poner_ts("fecha_utc", fu);
        df.poner_ts("fecha_sesion", ses);
        df.poner_f64("open", o_def(clean + "_PX_OPEN"));
        df.poner_f64("high", o_def(clean + "_PX_HIGH"));
        df.poner_f64("low", o_def(clean + "_PX_LOW"));
        df.poner_f64("close", precio);
        df.poner_f64("volume", vol);
        df.poner_f64("precio", precio);
        df.poner_str("fuente", std::vector<std::string>(df.filas(), "BLOOMBERG_CSV"));
        df.poner_str("simbolo_fuente", std::vector<std::string>(df.filas(), clean));
    } else {
        throw SimboloNoEncontradoError("No se encontraron repositorios locales de Bloomberg para '" + ticker + "'.");
    }
    filtrar_rango(df, "fecha_utc", fecha_inicio, fecha_fin);
    auto [mono, inc] = verificar_monotonicidad(df.col("fecha_utc").i);
    if (!mono)
        throw ErrorMonotoniaTiempo("VIOLACIÓN DE MONOTONÍA TEMPORAL EN CACHÉ BLOOMBERG D1 (" + ticker + "): "
                                   + (inc.empty() ? std::string("Índice temporal no monótono") : inc.front().descripcion) + ".");
    verificar_estampa_diaria(df);
    df.poner_i64("es_imputado", imputado_por_precios(df, {"open", "high", "low", "close", "precio"}));
    return df;
}

inline Tabla cargar_datos_bloomberg(const std::string& simbolo, const std::string& fecha_inicio, const std::string& fecha_fin,
                                    const std::string& resolucion = "1m", bool estampa_cierre = true) {
    if (resolucion == "1m" || resolucion == "intradia" || resolucion == "minuto") return cargar_bloomberg_1min(simbolo, fecha_inicio, fecha_fin, estampa_cierre);
    if (resolucion == "1d" || resolucion == "diario") return cargar_bloomberg_diario(simbolo, fecha_inicio, fecha_fin);
    throw ValueError("Resolución de Bloomberg '" + resolucion + "' no soportada (use '1m' o '1d').");
}

// ======================================================================================================================= YAHOO
inline std::string dir_cache_sinergia() { return unir(raiz_datos(), "cachas/sinergia_2_0"); }

inline const std::map<std::string, std::string>& mapa_archivos_yahoo() {
    static const std::map<std::string, std::string> M = {
        {"SPX", "SPX.csv"}, {"^GSPC", "SPX.csv"}, {"NDX", "NDX.csv"}, {"^NDX", "NDX.csv"}, {"RTY", "RTY.csv"}, {"^RUT", "RTY.csv"},
        {"DXY", "DXY.csv"}, {"DX-Y.NYB", "DXY.csv"}, {"HYG", "HYG.csv"}, {"IEF", "IEF.csv"}, {"KRE", "KRE.csv"}, {"EEM", "EEM.csv"},
        {"REAL10Y", "REAL10Y.csv"}, {"VIX", "VIX.csv"}, {"^VIX", "VIX.csv"}};
    return M;
}
inline bool sinergia_ajustado_retroactivo(const std::string& b) { return b == "EEM" || b == "HYG" || b == "IEF" || b == "KRE" || b == "TIP" || b == "VNQ"; }
inline constexpr double UMBRAL_DIVIDENDO_REL = 1e-4;

// Fila de trabajo de la tabla Yahoo (fecha de texto + columnas de valor).
struct TablaYahoo {
    std::vector<std::string> fecha;
    std::vector<double> valor, dividendo_rel;
    bool hay_dividendo = false, hay_fuente = false;
};

// merge(how='inner'/'left', on='fecha'): orden de la izquierda; cada clave de la izquierda se empareja con todas las de la
// derecha en su orden (las caches no traen fechas repetidas; si las trajeran, se replicaria igual que pandas).
inline std::multimap<std::string, size_t> indice_por_fecha(const std::vector<std::string>& f) {
    std::multimap<std::string, size_t> m;
    for (size_t k = 0; k < f.size(); ++k) m.emplace(f[k], k);
    return m;
}

inline TablaYahoo migrar_a_cierre_crudo(const TablaYahoo& sin, const std::string& base) {
    const CSV crudo = leer_csv(unir(dir_cache_macro_d1(), base + ".csv"));
    const auto fcr = crudo.texto("fecha");
    const auto close = crudo.numerica("close");
    const auto idx = indice_por_fecha(fcr);
    TablaYahoo o;
    std::vector<double> adj, raw;
    for (size_t k = 0; k < sin.fecha.size(); ++k) {
        auto [a, b] = idx.equal_range(sin.fecha[k]);
        for (auto it = a; it != b; ++it) { o.fecha.push_back(sin.fecha[k]); adj.push_back(sin.valor[k]); raw.push_back(close[it->second]); }
    }
    const size_t n = adj.size();
    std::vector<double> f(n), d(n, 0.0);
    for (size_t k = 0; k < n; ++k) f[k] = adj[k] / raw[k];
    for (size_t k = 1; k < n; ++k) d[k] = 1.0 - f[k - 1] / f[k];
    for (size_t k = 0; k < n; ++k) if (!(std::isfinite(d[k]) && std::fabs(d[k]) > UMBRAL_DIVIDENDO_REL)) d[k] = 0.0;
    o.valor = raw;
    o.dividendo_rel = d;
    o.hay_dividendo = o.hay_fuente = true;
    return o;
}

inline Tabla cargar_datos_yahoo_disco(const std::string& simbolo, const std::string& fecha_inicio, const std::string& fecha_fin) {
    const auto& M = mapa_archivos_yahoo();
    auto it = M.find(simbolo);
    const std::string archivo = it != M.end() ? it->second : simbolo + ".csv";
    const std::string ruta = unir(dir_cache_sinergia(), archivo);
    if (!es_fichero(ruta))
        throw SimboloNoEncontradoError("INCOMPATIBILIDAD DE CONTRATO (Venue=YAHOO, Símbolo='" + simbolo + "'): No se encontró la caché local en " + ruta + ".");
    const CSV raw = leer_csv(ruta);
    if (raw.indice("fecha") < 0 || raw.indice("valor") < 0)
        throw ValueError("Estructura inválida en archivo " + ruta + ": esperadas columnas 'fecha,valor'.");
    TablaYahoo y;
    y.fecha = raw.texto("fecha");
    y.valor = raw.numerica("valor");
    const std::string base = archivo.substr(0, archivo.rfind('.'));
    const std::string ruta_macro = unir(dir_cache_macro_d1(), base + ".csv");
    const bool migrada = sinergia_ajustado_retroactivo(base) && es_fichero(ruta_macro);
    if (migrada) y = migrar_a_cierre_crudo(y, base);
    auto [ses, disp] = estampar_diario_disponible(y.fecha);
    const int64_t t_ini = ns_desde_iso(fecha_inicio), t_fin = ns_desde_iso(fecha_fin);
    // Enriquecimiento OHLCV (merge left con macro_d1), dentro de un try/except: pass como en el Python.
    std::vector<std::pair<std::string, std::vector<double>>> extra_ohlcv;
    if (es_fichero(ruta_macro)) {
        try {
            const CSV mac = leer_csv(ruta_macro);
            const auto idx = indice_por_fecha(mac.texto("fecha"));
            std::vector<std::string> cols;
            for (const char* c : {"open", "high", "low", "volume"}) if (mac.indice(c) >= 0) cols.push_back(c);
            std::vector<std::vector<double>> src;
            for (const auto& c : cols) src.push_back(mac.numerica(c));
            // left merge: replica filas si hay fechas repetidas a la derecha
            TablaYahoo y2;
            std::vector<int64_t> s2, d2;
            std::vector<std::vector<double>> v2(cols.size());
            for (size_t k = 0; k < y.fecha.size(); ++k) {
                auto [a, b] = idx.equal_range(y.fecha[k]);
                auto empujar = [&](std::optional<size_t> j) {
                    y2.fecha.push_back(y.fecha[k]);
                    y2.valor.push_back(y.valor[k]);
                    if (y.hay_dividendo) y2.dividendo_rel.push_back(y.dividendo_rel[k]);
                    s2.push_back(ses[k]); d2.push_back(disp[k]);
                    for (size_t c = 0; c < cols.size(); ++c) v2[c].push_back(j ? src[c][*j] : NaN);
                };
                if (a == b) empujar(std::nullopt);
                else for (auto q = a; q != b; ++q) empujar(q->second);
            }
            y2.hay_dividendo = y.hay_dividendo;
            y2.hay_fuente = y.hay_fuente;
            y = y2; ses = s2; disp = d2;
            for (size_t c = 0; c < cols.size(); ++c) extra_ohlcv.emplace_back(cols[c], v2[c]);
        } catch (const std::exception&) {
            extra_ohlcv.clear();
        }
    }
    std::vector<uint8_t> m(disp.size());
    for (size_t k = 0; k < m.size(); ++k) m[k] = (disp[k] >= t_ini && disp[k] <= t_fin) ? 1 : 0;
    Tabla df;
    df.poner_ts("fecha_utc", disp);
    df.poner_ts("fecha_sesion", ses);
    df.poner_f64("precio", y.valor);
    for (const char* c : {"open", "high", "low", "volume"})
        for (const auto& [nom, v] : extra_ohlcv) if (nom == c) df.poner_f64(nom, v);
    if (y.hay_dividendo) df.poner_f64("dividendo_rel", y.dividendo_rel);
    if (y.hay_fuente) df.poner_str("fuente", std::vector<std::string>(y.fecha.size(), "YAHOO_VIA_MACRO_D1"));
    df.filtrar(m);
    auto [mono, inc] = verificar_monotonicidad(df.col("fecha_utc").i);
    if (!mono)
        throw ErrorMonotoniaTiempo("VIOLACIÓN DE MONOTONÍA TEMPORAL EN CACHÉ YAHOO (" + simbolo + "): "
                                   + (inc.empty() ? std::string("Índice temporal no monótono") : inc.front().descripcion) + ".");
    std::vector<int64_t> imp(df.filas());
    const auto& pv = df.col("precio").f;
    for (size_t k = 0; k < imp.size(); ++k) imp[k] = (std::isnan(pv[k]) || std::isinf(pv[k]) || pv[k] <= 0.0) ? 1 : 0;
    df.poner_i64("es_imputado", imp);
    df.poner_str("tipo_precio", std::vector<std::string>(df.filas(), (sinergia_ajustado_retroactivo(base) && !migrada) ? "AJUSTADO_RETROACTIVO" : "CIERRE"));
    return df;
}

}  // namespace ark::dataset::ingesta
