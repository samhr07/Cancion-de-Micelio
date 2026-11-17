// motor_ingesta.hpp -- CALCA DE recopilacion_unificada/motor_ingesta.py (ARK-DATA-02-CXX)
// Reactor ARK / Cancion del Micelio -- Jefe de Arquitectos de Codigo, 2026-11-02
//
// Porte tal cual (orden del Comandante): reloj intrinseco nu (calcular_incrementos_nu_causales, calibrar_umbral_nu_causal),
// CausalWindowEngine, alineacion multivariada (INTERSECCION_ESTRICTA / UNION_CON_MASCARA), carga de crudos por venue,
// rejilla diaria con anotacion de dias, lista separada de dias no regulares, procesamiento, cobertura diaria, testigo
// anti-lookahead sobre el pipeline real (auditar_causalidad_pipeline con el PCG64 de numpy) y de fuente, libro de hashes,
// manifiesto de integridad y ejecutar_extraccion. Tambien la parte de calendarios.py (advertir_auditoria_anual) y de
// ventanas_informacion.py (ventana_barra.cierre_hi y tabla_manifiesto) que el motor usa.
// Semantica de pandas 3 que se respeta: concat de una columna de texto, booleana o de fechas con NaN -> dtype object (que la
// rejilla rellena hacia delante); concat de enteros con NaN -> float; merge exterior: int -> float, bool -> object, texto y
// fechas conservan su tipo con nulos.
#pragma once
#include "ark/dataset/calendario.hpp"
#include "ark/dataset/control_calidad.hpp"
#include "ark/dataset/guardas_causales.hpp"

namespace ark::dataset::ingesta {

using ark::dataset::Calendario;
using ark::dataset::ClaseDia;
using ark::dataset::TipoDia;

inline const std::vector<std::string>& columnas_anotacion_dia() {
    static const std::vector<std::string> v{"tipo_dia", "nombre_dia", "fraccion_sesion", "hay_dato", "es_sesion_con_dato", "calendario",
                                            "dias_calendario_acumulados", "sesiones_acumuladas", "dias_acumulados", "estado_supervivencia",
                                            "cierre_real_utc"};
    return v;
}

inline double delta_reloj_segundos(TipoReloj r) {
    switch (r) {
        case TipoReloj::PARED_SEGUNDOS: return 1.0;
        case TipoReloj::PARED_MINUTOS: return 60.0;
        case TipoReloj::PARED_HORAS: return 3600.0;
        case TipoReloj::PARED_DIAS: return 86400.0;
        default: return 86400.0;   // MAPA_DELTA_RELOJ_SEGUNDOS.get(reloj, 86400.0)
    }
}

// ================================================================================================= calendarios (lote)
inline Calendario calendario_para(const std::string& sim, const std::optional<std::string>& fuente) {
    return ark::dataset::calendario_de(sim.c_str(), fuente ? fuente->c_str() : nullptr);
}
inline std::vector<int64_t> sesiones_cal(const Calendario& cal, int64_t z_ini, int64_t z_fin) {
    std::vector<int64_t> v;
    for (int64_t z = z_ini; z <= z_fin; ++z) if (cal.es_dia_operable(z) && !cal.es_feriado(z)) v.push_back(z);
    return v;
}
inline bool verificado_en_z(const Calendario& cal, int64_t z_ini, int64_t z_fin) {
    return cal.verificado_en(yyyymmdd_desde_dias(z_ini), yyyymmdd_desde_dias(z_fin));
}
inline std::optional<std::string> advertir_auditoria_anual(const Calendario& cal, int64_t z_fin) {
    const uint32_t h = ark::dataset::verificado_hasta(cal.id());
    if (h == 0 || yyyymmdd_desde_dias(z_fin) <= h) return std::nullopt;
    return std::string("AUDITORIA ANUAL REQUERIDA (") + cal.nombre() + "): el rango llega a " + texto_fecha_z(z_fin)
           + " y el calendario solo esta verificado hasta " + texto_fecha_z(dias_desde_yyyymmdd(h))
           + ". Sus reglas discrecionales (SIFMA: Viernes Santo con NFP) deben revisarse para los anos nuevos; una sesion en un dia que "
             "el calendario da por cerrado saldra como 'fuera de calendario', no como corrupcion.";
}
struct ClasLote { std::vector<std::string> tipo, nombre; std::vector<double> frac; };
inline ClasLote clasificar_lote(const Calendario& cal, const std::vector<int64_t>& ns) {
    ClasLote r;
    for (int64_t t : ns) {
        const ClaseDia c = cal.clasificar(piso_div(t, NS_DIA));
        r.tipo.push_back(ark::dataset::nombre_tipo_dia(c.tipo));
        r.nombre.push_back(c.nombre);
        r.frac.push_back(c.fraccion);
    }
    return r;
}
inline bool es_tipo_sesion_str(const std::string& t) { return t == "SESION_REGULAR" || t == "SESION_REDUCIDA"; }

// ======================================================================================= ventanas_informacion.py (manifiesto)
namespace detalle_vi {
struct TextoPerfil { const char* nombre; const char* sesion_reloj; const char* acta; const char* nota; };
inline const char* ACTA_COBERTURA = "boveda_fase1/actas/L5_mapas_red_macro/acta_bordes_diarios.txt (adendas f, g, h)";
inline const char* ACTA_FX = "JEFE DE ARQUITECTOS DE CODIGO/AUDITORIA_CAUSAL/COBERTURA_CALENDARIOS.md (JPY=X contra el minuto de Bloomberg) y "
                             "boveda_fase1/actas/L5_mapas_red_macro/acta_i10_hora_dxy.txt (Yahoo JPY=X 20:00 NY)";
inline const TextoPerfil& texto_perfil(const std::string& n) {
    static const TextoPerfil T[] = {
        {"NYSE_RV", "EEUU", ACTA_COBERTURA,
         "Renta variable y ETF de EE. UU.: subasta de cierre 16:00 ET (13:00 en media jornada). Coincide con la cobertura "
         "medida (EEUU 16.0, error 0.0 h)."},
        {"NYSE_INDICE", "EEUU", ACTA_COBERTURA,
         "Índice calculado (SPX, NDX, RTY): su apertura no es un precio negociable y en datos viejos de Yahoo suele ser "
         "el cierre anterior."},
        {"VIX", "EEUU", ACTA_COBERTURA,
         "Cboe publica a las 16:15 ET; la cobertura midió 17:00 (desfase de una hora no resuelto, I72). Horario extendido "
         "desde las 03:15 ET."},
        {"RENDIMIENTO_CBOE", "EEUU", "reloj_macro.SESION (asignado a EEUU sin medición propia)",
         "Índices de rendimiento de Cboe en Yahoo (^FVX, ^TNX, ^TYX). Hora no medida."},
        {"TESORO", "EEUU", "reloj_macro.SESION (REAL10Y y US02Y en EEUU, sin medición)",
         "Curva del Tesoro (fixing de fin de día; la metodología habla de ~15:30 ET, no leída aquí en fuente primaria). "
         "No tiene apertura."},
        {"BLOOMBERG_RENDIMIENTO", "EEUU", "acta_hora_fijacion.txt (USGG30YR fija a las 16:00 ET; USGG10YR sin medir)",
         "Rendimiento genérico de Bloomberg (USGG10YR)."},
        {"CME_ENERGIA_METAL", nullptr, ACTA_COBERTURA,
         "Futuros de COMEX, NYMEX y Brent de ICE en Yahoo: la barra va de las 18:00 NY de la víspera a las 17:00 NY. Su "
         "apertura NO es la mañana: empieza la tarde anterior."},
        {"CBOT_GRANO", "CBOT", ACTA_COBERTURA, "Granos de CBOT en Yahoo: la barra va de las 20:00 NY de la víspera a las 15:00 NY."},
        {"CME_RV_BLOOMBERG", nullptr, "sin medición (ES1 de Bloomberg)",
         "E-mini S&P de Globex: reabre a las 18:00 NY de la víspera; la hora de su PX_LAST diario no está medida."},
        {"DXY", "FX", "acta_i10_hora_dxy.txt (IC 15:03-18:26 NY) y reloj_macro (FX en [14, 17] NY)",
         "Índice del dólar: barras en las 24 horas, sin sesión; su cierre diario es una convención del proveedor. "
         "Intervalo = unión del de reloj_macro y del IC de I10."},
        {"FX_YAHOO", nullptr, ACTA_FX,
         "Cambio de Yahoo (=X): la barra etiquetada D cierra a las 00:00 UTC de D+1 (medido contra el minuto de "
         "Bloomberg, también en verano británico). Diverge de reloj_macro, que describe la convención de 17:00 NY de "
         "Bloomberg: es otro proveedor."},
        {"FX_YAHOO_ANALOGIA", nullptr, ACTA_FX,
         "Otro par de Yahoo (USDCNH = CNH=X): misma convención esperada que JPY=X, no medida; intervalo desde las 17:00 NY "
         "hasta las 00:00 UTC."},
        {"CRIPTO_UTC", nullptr, "adaptador_binance (close_time de la vela 1d)", "Vela diaria de Binance: [D 00:00 UTC, D+1 00:00 UTC)."},
        {"CRIPTO_NY16", nullptr, "FERIADOS_Y_PUNTOS_DEL_ECONOMISTA.md § 3.2 (cierre verificado al minuto)",
         "Corte NY_16: el día D cubre (D-1 16:00 NY, D 16:00 NY], horario de verano incluido."},
        {"CRIPTO_NO_VERIFICADO", nullptr, "sin medición", "Vela diaria de cripto de Yahoo o de Bitfinex: se espera el día UTC; no está medido."},
        {"CRIPTO_BLOOMBERG", nullptr, "sin medición (XBTUSD de Bloomberg)",
         "Bitcoin de Bloomberg: su convención diaria no está leída; intervalo desde las 16:00 NY hasta las 00:00 UTC."},
        {"EUROPA_INDICE", "EUROPA", ACTA_COBERTURA, "SX5E de Yahoo (calendario suizo, § 10 de la bitácora). Índice calculado."},
        {"DESCONOCIDO", nullptr, "", "Instrumento sin perfil: su cierre puede caer en cualquier momento del día UTC de la sesión."},
    };
    for (const auto& t : T) if (n == t.nombre) return t;
    return T[sizeof(T) / sizeof(T[0]) - 1];
}
inline const char* zona_iana(ark::dataset::Zona z) {
    switch (z) {
        case ark::dataset::Zona::NUEVA_YORK: return "America/New_York";
        case ark::dataset::Zona::BERLIN: return "Europe/Berlin";
        default: return "UTC";
    }
}
inline const char* proc_nombre(ark::dataset::Procedencia p) {
    switch (p) {
        case ark::dataset::Procedencia::REGLA: return "REGLA";
        case ark::dataset::Procedencia::MEDIDO: return "MEDIDO";
        default: return "NO_VERIFICADO";
    }
}
inline JV instante(const ark::dataset::Instante& i) {
    JV j = JV::objeto();
    j.poner("zona", JV::cadena(zona_iana(i.zona)));
    j.poner("dias", JV::entero(i.dias));
    j.poner("hora", JV::real(i.hora));
    return j;
}
inline JV par_instantes(const ark::dataset::Instante& a, const ark::dataset::Instante& b) {
    JV l = JV::lista();
    l.push(instante(a));
    l.push(instante(b));
    return l;
}
}  // namespace detalle_vi

// tabla_manifiesto: un registro por activo con su perfil, horas locales y procedencia.
inline JV tabla_manifiesto_ventanas(const std::vector<std::pair<std::string, std::optional<std::string>>>& sim_fuente, CorteDiario corte) {
    JV filas = JV::lista();
    for (const auto& [sim, fuente] : sim_fuente) {
        const auto& p = ark::dataset::perfil_de(sim.c_str(), fuente ? fuente->c_str() : nullptr, corte);
        const auto& tx = detalle_vi::texto_perfil(p.nombre);
        JV r = JV::objeto();
        r.poner("simbolo", JV::cadena(sim));
        r.poner("fuente", fuente ? JV::cadena(*fuente) : JV::nulo());
        r.poner("perfil", JV::cadena(p.nombre));
        r.poner("sesion_reloj", tx.sesion_reloj ? JV::cadena(tx.sesion_reloj) : JV::nulo());
        r.poner("cierre", detalle_vi::par_instantes(p.cierre_lo, p.cierre_hi));
        r.poner("proc_cierre", JV::cadena(detalle_vi::proc_nombre(p.proc_cierre)));
        r.poner("cierre_reducido", p.tiene_cierre_reducido ? detalle_vi::par_instantes(p.red_lo, p.red_hi) : JV::nulo());
        r.poner("proc_cierre_reducido", JV::cadena(detalle_vi::proc_nombre(p.proc_cierre_reducido)));
        if (p.inicio_es_cierre_previo) r.poner("inicio", JV::cadena("cierre de la sesion previa"));
        else r.poner("inicio", p.tiene_inicio ? detalle_vi::par_instantes(p.ini_lo, p.ini_hi) : JV::nulo());
        r.poner("proc_inicio", JV::cadena(detalle_vi::proc_nombre(p.proc_inicio)));
        r.poner("negociable", JV::booleano(p.negociable));
        r.poner("acta", JV::cadena(tx.acta));
        r.poner("nota", JV::cadena(tx.nota));
        r.poner("aviso", JV::cadena(simbolo_limpio(sim.c_str()) == "EEM"
                                        ? "EEM cotiza en NY pero su subyacente cierra en Asia: valor liquidativo rancio (reloj_macro.SUBYACENTE_RANCIO)."
                                        : ""));
        filas.push(std::move(r));
    }
    return filas;
}

// ventana_barra(simbolo, fecha_sesion, fuente, corte).cierre_hi en ns UTC (ValueError como en Python si no hay sesion previa).
inline int64_t cierre_real_utc_ns(const std::string& sim, int64_t fecha_sesion_ns, const std::optional<std::string>& fuente, CorteDiario corte) {
    const int64_t z = piso_div(fecha_sesion_ns, NS_DIA);
    const auto v = ark::dataset::ventana_barra(sim.c_str(), yyyymmdd_desde_dias(z), fuente ? fuente->c_str() : nullptr, corte);
    if (!v.ok) {
        const Calendario cal = ark::dataset::calendario_de(simbolo_limpio(sim.c_str()).c_str(), fuente ? fuente->c_str() : nullptr);
        int64_t zi = z;
        const bool sin = !ark::dataset::es_tipo_sesion(cal.clasificar(z).tipo);
        const int64_t zd = (sin && cal.sesion_previa(z, zi)) ? zi : z;
        throw ValueError("Sin sesión en las dos semanas anteriores a " + texto_fecha_z(zd) + " (" + cal.nombre() + ").");
    }
    return v.cierre_hi_ms * 1000000LL;
}

// ============================================================================================= reloj intrinseco (Axioma 4)
inline std::vector<double> calcular_incrementos_nu_causales(const Tabla& df, TipoReloj reloj, int64_t lookback_W = 60) {
    const size_t n = df.filas();
    if (n == 0) return {};
    std::vector<uint8_t> activo(n, 1);
    if (df.tiene("es_imputado")) {
        const auto v = a_dobles(df.col("es_imputado"));
        for (size_t k = 0; k < n; ++k) activo[k] = static_cast<int64_t>(std::isnan(v[k]) ? 0.0 : v[k]) == 0;
    }
    if (df.tiene("market_open_mask")) {
        const auto v = a_dobles(df.col("market_open_mask"));
        for (size_t k = 0; k < n; ++k) activo[k] = activo[k] && static_cast<int64_t>(std::isnan(v[k]) ? 1.0 : v[k]) == 1;
    }
    if (es_pared(reloj)) return std::vector<double>(n, 1.0);
    auto num0 = [&](const std::string& c) {   // to_numeric(...).fillna(0.0).astype(float)
        std::vector<double> v = a_dobles(df.col(c));
        for (auto& x : v) if (std::isnan(x)) x = 0.0;
        return v;
    };
    auto limpiar = [&](std::vector<double> v) {
        for (size_t k = 0; k < n; ++k) v[k] = (std::isfinite(v[k]) && v[k] > 0.0 && activo[k]) ? v[k] : 0.0;
        return v;
    };
    const std::string col_vol = df.tiene("volume") ? "volume" : (df.tiene("volumen") ? "volumen" : "");
    if (reloj == TipoReloj::NU_VOLUMEN_BASE) {
        if (col_vol.empty()) throw IncompatibilidadContratoError("Se requiere columna 'volume' o 'volumen' para el reloj intrínseco NU_VOLUMEN_BASE.");
        return limpiar(num0(col_vol));
    }
    if (reloj == TipoReloj::NU_VOLUMEN_NOCIONAL) {
        std::vector<double> noc;
        if (df.tiene("quote_volume")) noc = num0("quote_volume");
        else {
            if (col_vol.empty())
                throw IncompatibilidadContratoError("Se requiere columna 'quote_volume' o 'volume'/'volumen' para el reloj NU_VOLUMEN_NOCIONAL.");
            const auto px = num0("precio"), vol = num0(col_vol);
            noc.resize(n);
            for (size_t k = 0; k < n; ++k) noc[k] = px[k] * vol[k];
        }
        return limpiar(noc);
    }
    if (reloj == TipoReloj::NU_TRADES) {
        std::vector<double> tr;
        if (df.tiene("trades_count")) tr = num0("trades_count");
        else if (df.tiene("n")) tr = num0("n");
        else if (!col_vol.empty()) { const auto v = num0(col_vol); tr.resize(n); for (size_t k = 0; k < n; ++k) tr[k] = v[k] > 0.0 ? 1.0 : 0.0; }
        else tr.assign(n, 1.0);
        return limpiar(tr);
    }
    if (reloj == TipoReloj::NU_DELTA_INFO) {
        const std::vector<double> precios = a_dobles(df.col("precio"));
        std::vector<double> logp(n), lr(n, 0.0), di(n, 0.0);
        for (size_t k = 0; k < n; ++k) logp[k] = std::isnan(precios[k]) ? NaN : std::log(std::max(1e-12, precios[k]));
        for (size_t k = 1; k < n; ++k) lr[k] = logp[k] - logp[k - 1];
        const double inv_2ln2 = 1.0 / (2.0 * std::log(2.0));
        const int64_t w = std::max<int64_t>(2, lookback_W);
        for (size_t idx = 1; idx < n; ++idx) {
            if (!activo[idx]) continue;
            const size_t ini = static_cast<size_t>(std::max<int64_t>(1, static_cast<int64_t>(idx) - w));
            if (idx - ini < 2) continue;
            const std::vector<double> h(lr.begin() + static_cast<std::ptrdiff_t>(ini), lr.begin() + static_cast<std::ptrdiff_t>(idx));
            const double mu = np_mean(h), sd = np_std(h);
            const double med = np_median(h);
            std::vector<double> dev(h.size());
            for (size_t j = 0; j < h.size(); ++j) dev[j] = std::fabs(h[j] - med);
            const double mad = 1.4826 * np_median(dev);
            const double sigma = py_max({sd, mad, 1e-8});
            const double z = std::fabs(lr[idx] - mu) / sigma;
            di[idx] = py_min(50.0, std::pow(z, 2.0) * inv_2ln2);
        }
        return di;
    }
    throw IncompatibilidadContratoError(std::string("Tipo de reloj desconocido: ") + nombre(reloj));
}

inline double calibrar_umbral_nu_causal(const std::vector<double>& dnu, int64_t W, int64_t S, std::optional<double> umbral) {
    if (umbral) {
        if (!std::isfinite(*umbral) || *umbral <= 0.0) throw IncompatibilidadContratoError("umbral_nu (" + repr_py(*umbral) + ") debe ser finito y > 0.0.");
        return *umbral;
    }
    if (dnu.empty()) return 1.0;
    const size_t m = static_cast<size_t>(std::max<int64_t>(1, std::min<int64_t>(static_cast<int64_t>(dnu.size()), W)));
    std::vector<double> pos;
    for (size_t k = 0; k < m; ++k) if (dnu[k] > 0.0) pos.push_back(dnu[k]);
    const double media = pos.empty() ? 1.0 : np_mean(pos);
    return std::max(1e-12, static_cast<double>(std::max<int64_t>(1, S)) * media);
}

// CausalWindowEngine: en cada disparo entrega solo [t - W + 1, t].
struct CausalWindowEngine {
    Tabla df;
    int64_t W = 60, S = 1;
    TipoReloj tipo_reloj = TipoReloj::PARED_DIAS;
    std::vector<double> delta_nu;
    double umbral_nu = 1.0;
    std::vector<int64_t> indices_disparo;

    CausalWindowEngine(const Tabla& d, int64_t lookback_W, int64_t stride_S = 1, TipoReloj reloj = TipoReloj::PARED_DIAS,
                       std::optional<double> umbral = std::nullopt)
        : df(d), W(lookback_W), S(std::max<int64_t>(1, stride_S)), tipo_reloj(reloj) {
        delta_nu = calcular_incrementos_nu_causales(df, tipo_reloj, W);
        df.poner_f64("delta_nu", delta_nu);
        std::vector<double> acum(delta_nu.size());
        double s = 0.0;
        for (size_t k = 0; k < acum.size(); ++k) { s += delta_nu[k]; acum[k] = s; }
        df.poner_f64("nu_acumulado", std::move(acum));
        umbral_nu = calibrar_umbral_nu_causal(delta_nu, W, S, umbral);
        const int64_t N = static_cast<int64_t>(df.filas());
        if (N < W) return;
        if (es_pared(tipo_reloj)) { for (int64_t i = W - 1; i < N; i += S) indices_disparo.push_back(i); return; }
        double acu = 0.0;
        for (int64_t i = W - 1; i < N; ++i) {
            acu += delta_nu[static_cast<size_t>(i)];
            if (acu >= umbral_nu) { indices_disparo.push_back(i); acu = std::fmod(acu, umbral_nu); }
        }
    }
    size_t size() const { return indices_disparo.size(); }
    // stream(): (timestamp_actual, ventana_datos, mascara_imputacion) del disparo k.
    Tabla ventana(size_t k) const {
        const int64_t i = indices_disparo[k];
        std::vector<size_t> idx;
        for (int64_t j = i - W + 1; j <= i; ++j) idx.push_back(static_cast<size_t>(j));
        return tomar_filas(df, idx);
    }
    int64_t t_actual(size_t k) const { return df.col("fecha_utc").i[static_cast<size_t>(indices_disparo[k])]; }
};

// ================================================================================================ alineacion multivariada
namespace detalle_merge {
// Reune las filas `fila` (-1 = falta) de una columna con la semantica de pandas para las faltantes.
inline Columna reunir(const Columna& c, const std::vector<int64_t>& fila) {
    bool falta = false;
    for (int64_t r : fila) if (r < 0) { falta = true; break; }
    Columna n;
    n.nombre = c.nombre;
    n.tipo = c.tipo;
    n.objeto = c.objeto;
    const size_t m = fila.size();
    switch (c.tipo) {
        case TipoColumna::F64:
            n.f.resize(m);
            for (size_t k = 0; k < m; ++k) n.f[k] = fila[k] < 0 ? NaN : c.f[static_cast<size_t>(fila[k])];
            break;
        case TipoColumna::I64:
            if (falta) {
                n.tipo = TipoColumna::F64;
                n.f.resize(m);
                for (size_t k = 0; k < m; ++k) n.f[k] = fila[k] < 0 ? NaN : static_cast<double>(c.i[static_cast<size_t>(fila[k])]);
            } else {
                n.i.resize(m);
                for (size_t k = 0; k < m; ++k) n.i[k] = c.i[static_cast<size_t>(fila[k])];
            }
            break;
        case TipoColumna::BOOL:
            if (falta) {   // bool con NaN -> object
                n.tipo = TipoColumna::STR;
                n.objeto = true;
                n.s.resize(m);
                n.sn.assign(m, 0);
                for (size_t k = 0; k < m; ++k) {
                    if (fila[k] < 0) n.sn[k] = 1;
                    else n.s[k] = c.i[static_cast<size_t>(fila[k])] ? "True" : "False";
                }
            } else {
                n.i.resize(m);
                for (size_t k = 0; k < m; ++k) n.i[k] = c.i[static_cast<size_t>(fila[k])];
            }
            break;
        case TipoColumna::TS:
            n.i.resize(m);
            for (size_t k = 0; k < m; ++k) n.i[k] = fila[k] < 0 ? NAT : c.i[static_cast<size_t>(fila[k])];
            break;
        case TipoColumna::STR:
            n.s.resize(m);
            if (falta || !c.sn.empty()) n.sn.assign(m, 0);
            for (size_t k = 0; k < m; ++k) {
                if (fila[k] < 0) n.sn[k] = 1;
                else {
                    n.s[k] = c.s[static_cast<size_t>(fila[k])];
                    if (c.es_nulo_str(static_cast<size_t>(fila[k]))) n.sn[k] = 1;
                }
            }
            break;
    }
    return n;
}
}  // namespace detalle_merge

// pd.merge(L, R, on='fecha_utc', how=inner|outer) con claves unicas, ya ordenado por fecha_utc (como el sort_values posterior).
inline Tabla merge_fecha(const Tabla& L, const Tabla& R, bool outer) {
    const auto& kl = L.col("fecha_utc").i;
    const auto& kr = R.col("fecha_utc").i;
    std::map<int64_t, std::pair<int64_t, int64_t>> m;
    for (size_t k = 0; k < kl.size(); ++k) m[kl[k]].first = static_cast<int64_t>(k) + 1;
    for (size_t k = 0; k < kr.size(); ++k) m[kr[k]].second = static_cast<int64_t>(k) + 1;
    std::vector<int64_t> claves, fl, fr;
    for (const auto& [k, v] : m) {
        if (!outer && (v.first == 0 || v.second == 0)) continue;
        claves.push_back(k);
        fl.push_back(v.first - 1);
        fr.push_back(v.second - 1);
    }
    Tabla t;
    for (const auto& c : L.columnas) {
        if (c.nombre == "fecha_utc") { Columna k; k.nombre = "fecha_utc"; k.tipo = TipoColumna::TS; k.i = claves; t.columnas.push_back(std::move(k)); }
        else t.columnas.push_back(detalle_merge::reunir(c, fl));
    }
    for (const auto& c : R.columnas)
        if (c.nombre != "fecha_utc") t.columnas.push_back(detalle_merge::reunir(c, fr));
    return t;
}

// fillna(v).astype(int) de una columna numerica (o de un object de booleanos) -> I64.
inline std::vector<int64_t> llenar_entero(const Columna& c, int64_t v) {
    std::vector<int64_t> r(c.size());
    if (c.tipo == TipoColumna::STR) {
        for (size_t k = 0; k < r.size(); ++k) r[k] = c.es_nulo_str(k) ? v : (c.s[k] == "True" ? 1 : 0);
        return r;
    }
    const auto d = a_dobles(c);
    for (size_t k = 0; k < r.size(); ++k) r[k] = std::isnan(d[k]) ? v : static_cast<int64_t>(d[k]);
    return r;
}

struct Alineado { Tabla df; int64_t filas_excluidas = -1; };   // -1 = sin attrs["filas_excluidas"]

inline Alineado alinear_multivariado_causal(const Crudos& dfs, PoliticaCalendario politica) {
    Alineado out;
    if (dfs.size() == 1) {
        Tabla d = dfs.front().second;
        if (!d.tiene("market_open_mask")) {
            const size_t n = d.filas();
            std::vector<int64_t> m(n);
            if (d.tiene("es_sesion_con_dato")) { const auto& v = d.col("es_sesion_con_dato").i; for (size_t k = 0; k < n; ++k) m[k] = v[k] ? 1 : 0; }
            else {
                const auto v = a_dobles(d.col("es_imputado"));
                for (size_t k = 0; k < n; ++k) m[k] = std::clamp<int64_t>(1 - static_cast<int64_t>(v[k]), 0, 1);
            }
            d.poner_i64("market_open_mask", std::move(m));
        }
        out.df = std::move(d);
        return out;
    }
    const std::string sim0 = dfs.front().first;
    auto tabla = [&](const std::string& sim, const Tabla& d) {
        std::vector<std::string> cols{"fecha_utc", "precio", "es_imputado"};
        std::vector<std::string> anot;
        for (const auto& c : columnas_anotacion_dia()) if (d.tiene(c)) anot.push_back(c);
        cols.insert(cols.end(), anot.begin(), anot.end());
        Tabla s = d.seleccionar(cols);
        std::vector<int64_t> om(s.filas(), 1);
        if (d.tiene("es_sesion_con_dato")) { const auto& v = d.col("es_sesion_con_dato").i; for (size_t k = 0; k < om.size(); ++k) om[k] = v[k] ? 1 : 0; }
        s.poner_i64("open_mask_" + sim, std::move(om));
        s.renombrar("precio", "precio_" + sim);
        s.renombrar("es_imputado", "es_imputado_" + sim);
        for (const auto& c : anot) s.renombrar(c, c + "_" + sim);
        return s;
    };
    const bool inter = politica == PoliticaCalendario::INTERSECCION_ESTRICTA;
    Tabla m;
    bool primero = true;
    for (const auto& [sim, d] : dfs) {
        Tabla s = tabla(sim, d);
        m = primero ? std::move(s) : merge_fecha(m, s, !inter);
        primero = false;
    }
    const size_t n = m.filas();
    std::vector<std::string> cprecio, copen, cimp;
    for (const auto& [sim, d] : dfs) { cprecio.push_back("precio_" + sim); copen.push_back("open_mask_" + sim); cimp.push_back("es_imputado_" + sim); }
    if (inter) {
        std::vector<uint8_t> abiertos(n, 1);
        for (const auto& c : copen) { const auto v = a_dobles(m.col(c)); for (size_t k = 0; k < n; ++k) if (v[k] != 1.0) abiertos[k] = 0; }
        int64_t excl = 0;
        for (uint8_t a : abiertos) excl += a ? 0 : 1;
        m.filtrar(abiertos);
        const size_t r = m.filas();
        m.poner_f64("precio", a_dobles(m.col("precio_" + sim0)));
        std::vector<int64_t> imp(r, 0);
        for (const auto& c : cimp) { const auto v = a_dobles(m.col(c)); for (size_t k = 0; k < r; ++k) imp[k] = std::max(imp[k], static_cast<int64_t>(v[k])); }
        m.poner_i64("es_imputado", std::move(imp));
        m.poner_i64("market_open_mask", std::vector<int64_t>(r, 1));
        for (const auto& c : copen) m.quitar(c);
        out.df = std::move(m);
        out.filas_excluidas = excl;
        return out;
    }
    for (const auto& c : copen) m.poner_i64(c, llenar_entero(m.col(c), 0));
    for (const auto& c : cimp) m.poner_i64(c, llenar_entero(m.col(c), 1));
    for (const auto& [sim, d] : dfs)
        for (const char* c : {"hay_dato", "es_sesion_con_dato"}) {
            const std::string nom = std::string(c) + "_" + sim;
            if (!m.tiene(nom)) continue;
            m.poner_bool(nom, llenar_entero(m.col(nom), 0));
        }
    std::vector<uint8_t> con_hist(n, 1);
    for (const auto& c : cprecio) {
        auto& v = m.col(c).f;
        for (size_t k = 1; k < n; ++k) if (std::isnan(v[k])) v[k] = v[k - 1];
        for (size_t k = 0; k < n; ++k) if (std::isnan(v[k])) con_hist[k] = 0;
    }
    bool alguno = false;
    for (uint8_t x : con_hist) alguno = alguno || x;
    if (!alguno)
        throw ViolacionCausalidadError("CENSURA IZQUIERDA MULTIVARIADA: No existe ningún instante con historia causal válida para todos los activos.");
    m.filtrar(con_hist);
    const size_t r = m.filas();
    std::vector<int64_t> mom(r, std::numeric_limits<int64_t>::max()), imp(r, 0);
    for (const auto& c : copen) { const auto& v = m.col(c).i; for (size_t k = 0; k < r; ++k) mom[k] = std::min(mom[k], v[k]); }
    for (const auto& c : cimp) { const auto& v = m.col(c).i; for (size_t k = 0; k < r; ++k) imp[k] = std::max(imp[k], v[k]); }
    std::vector<int64_t> ei(r);
    for (size_t k = 0; k < r; ++k) ei[k] = (mom[k] == 0 || imp[k] == 1) ? 1 : 0;
    m.poner_i64("market_open_mask", std::move(mom));
    m.poner_i64("es_imputado", std::move(ei));
    m.poner_f64("precio", m.col("precio_" + sim0).f);
    out.df = std::move(m);
    return out;
}

// ========================================================================================================= carga de crudos
inline std::pair<Crudos, double> cargar_crudos(const ContratoExtraccion& c) {
    double delta = delta_reloj_segundos(c.tipo_reloj);
    Crudos crudos;
    for (const auto& sim : c.simbolos) {
        Tabla df;
        if (c.venue == Venue::YAHOO) {
            df = cargar_datos_yahoo_disco(sim, c.fecha_inicio, c.fecha_fin);
            delta = 86400.0;
        } else if (c.venue == Venue::BLOOMBERG_XLSX) {
            if (c.tipo_reloj == TipoReloj::PARED_MINUTOS || c.tipo_reloj == TipoReloj::PARED_SEGUNDOS) {
                delta = 60.0;
                df = cargar_datos_bloomberg(sim, c.fecha_inicio, c.fecha_fin, "1m");
            } else if (!es_pared(c.tipo_reloj)) {
                try {
                    df = cargar_datos_bloomberg(sim, c.fecha_inicio, c.fecha_fin, "1m");
                    delta = 60.0;
                } catch (const SimboloNoEncontradoError&) {
                    df = cargar_datos_bloomberg(sim, c.fecha_inicio, c.fecha_fin, "1d");
                    delta = 86400.0;
                }
            } else {
                delta = 86400.0;
                df = cargar_datos_bloomberg(sim, c.fecha_inicio, c.fecha_fin, "1d");
            }
        } else if (c.venue == Venue::BINANCE_UM) {
            std::string res;
            if (c.tipo_reloj == TipoReloj::PARED_MINUTOS || c.tipo_reloj == TipoReloj::PARED_SEGUNDOS || !es_pared(c.tipo_reloj)
                || c.nivel == NivelDato::L1_ENRIQUECIDO) {
                res = "1m";
                delta = 60.0;
            } else if (c.tipo_reloj == TipoReloj::PARED_HORAS) {
                res = "1h";
                delta = 3600.0;
            } else {
                res = c.corte_diario_cripto == CorteDiario::NY_16 ? "1d_ny16" : "1d";
                delta = 86400.0;
            }
            df = cargar_datos_binance(sim, c.fecha_inicio, c.fecha_fin, res);
        } else if (c.venue == Venue::BITFINEX) {
            delta = 1.0;
            df = cargar_datos_bitfinex(sim, c.fecha_inicio, c.fecha_fin, c.nivel);
        } else {
            throw IncompatibilidadContratoError(std::string("INCOMPATIBILIDAD DE SEDE (Axioma 2): El venue '") + nombre(c.venue)
                                                + "' no cuenta con repositorio local activo en este nodo del Reactor ARK.");
        }
        if (df.vacia())
            throw IncompatibilidadContratoError("INCOMPATIBILIDAD DE RANGO TEMPORAL (Axioma 2): El símbolo '" + sim + "' en el venue '" + nombre(c.venue)
                                                + "' no contiene muestras en el rango [" + c.fecha_inicio + ", " + c.fecha_fin + "].");
        crudos.emplace_back(sim, std::move(df));
    }
    return {crudos, delta};
}

// Sesiones cuya disponibilidad (sesion + 1 d) cae en [fecha_inicio, fecha_fin], redondeado hacia dentro a dias (ns).
inline std::pair<int64_t, int64_t> rango_sesiones_contrato(const ContratoExtraccion& c) {
    int64_t lo = ns_desde_iso(recortar(c.fecha_inicio)) - NS_DIA;
    const int64_t lo_n = normalizar_dia(lo);
    lo = lo_n + (lo != lo_n ? NS_DIA : 0);
    const int64_t hi = normalizar_dia(ns_desde_iso(recortar(c.fecha_fin)) - NS_DIA);
    return {lo, hi};
}

inline std::optional<std::string> fuente_de(const Tabla& df) {
    if (!df.tiene("fuente") || df.vacia()) return std::nullopt;
    return celda_texto(df.col("fuente"), 0);
}

// ===================================================================================================== rejilla diaria
namespace detalle_rejilla {
inline bool es_cero(const std::string& c) {
    static const std::set<std::string> s{"volume", "volumen", "quote_volume", "trades_count", "taker_buy_base", "taker_buy_quote", "dividendo_rel"};
    return s.count(c) > 0;
}
// Columna -> object (cadenas como str(x) de Python; nulo aparte).
inline void a_objeto(Columna& c) {
    if (c.tipo == TipoColumna::STR) { c.objeto = true; if (c.sn.empty()) c.sn.assign(c.s.size(), 0); return; }
    const size_t n = c.size();
    Columna o;
    o.nombre = c.nombre;
    o.tipo = TipoColumna::STR;
    o.objeto = true;
    o.s.resize(n);
    o.sn.assign(n, 0);
    for (size_t k = 0; k < n; ++k) {
        if (celda_nula(c, k)) o.sn[k] = 1;
        else o.s[k] = celda_texto(c, k);
    }
    c = std::move(o);
}
// pd.concat([d, filas[d.columns]]): filas con un valor escalar por columna (`nulo` = NaN; si no, `valor` real o entero).
inline void anadir(Columna& c, size_t nn, bool nulo, double valor, bool entero) {
    switch (c.tipo) {
        case TipoColumna::F64: for (size_t k = 0; k < nn; ++k) c.f.push_back(nulo ? NaN : valor); break;
        case TipoColumna::I64:
            if (entero && !nulo) { for (size_t k = 0; k < nn; ++k) c.i.push_back(static_cast<int64_t>(valor)); }
            else {
                std::vector<double> v = a_dobles(c);
                c.i.clear();
                c.tipo = TipoColumna::F64;
                c.f = std::move(v);
                for (size_t k = 0; k < nn; ++k) c.f.push_back(nulo ? NaN : valor);
            }
            break;
        default:   // BOOL, TS, STR con un escalar numerico o NaN -> object
            a_objeto(c);
            for (size_t k = 0; k < nn; ++k) {
                c.s.push_back(nulo ? std::string() : (entero ? std::to_string(static_cast<int64_t>(valor)) : repr_py(valor)));
                c.sn.push_back(nulo ? 1 : 0);
            }
            break;
    }
}
}  // namespace detalle_rejilla

inline Tabla aplicar_rejilla_diaria(const ContratoExtraccion& c, const std::string& sim, const Tabla& df) {
    const std::optional<std::string> fuente = fuente_de(df);
    const Calendario cal = calendario_para(sim, fuente);
    Tabla d = df;
    d.poner_bool("_original", std::vector<int64_t>(d.filas(), 1));
    if (c.rejilla_diaria != RejillaDiaria::SOLO_SESIONES && !d.vacia()) {
        const int64_t hi = rango_sesiones_contrato(c).second;
        const int64_t ayer = normalizar_dia(ahora_ns()) - NS_DIA;
        const auto& fs = d.col("fecha_sesion").i;
        const int64_t mn = *std::min_element(fs.begin(), fs.end()), mx = *std::max_element(fs.begin(), fs.end());
        const int64_t hi_g = std::max(std::min(hi, ayer), mx);
        const std::set<int64_t> exist(fs.begin(), fs.end());
        std::vector<int64_t> nuevos;
        for (int64_t t = mn; t <= hi_g; t += NS_DIA) {
            if (c.rejilla_diaria == RejillaDiaria::SESIONES_Y_FERIADOS && !cal.es_dia_operable(piso_div(t, NS_DIA))) continue;
            if (!exist.count(t)) nuevos.push_back(t);
        }
        if (!nuevos.empty()) {
            const size_t nn = nuevos.size();
            for (auto& col : d.columnas) {
                if (col.nombre == "fecha_sesion") { for (int64_t t : nuevos) col.i.push_back(t); }
                else if (col.nombre == "fecha_utc") { for (int64_t t : nuevos) col.i.push_back(t + NS_DIA); }
                else if (col.nombre == "_original") { for (size_t k = 0; k < nn; ++k) col.i.push_back(0); }
                else if (col.nombre == "es_imputado") detalle_rejilla::anadir(col, nn, false, 1.0, true);
                else if (detalle_rejilla::es_cero(col.nombre)) detalle_rejilla::anadir(col, nn, false, 0.0, false);
                else detalle_rejilla::anadir(col, nn, true, 0.0, false);
            }
            const auto& f2 = d.col("fecha_sesion").i;
            std::vector<size_t> orden(f2.size());
            std::iota(orden.begin(), orden.end(), size_t{0});
            std::stable_sort(orden.begin(), orden.end(), [&](size_t a, size_t b) { return f2[a] < f2[b]; });
            d = tomar_filas(d, orden);
            for (auto& col : d.columnas) {
                if (col.tipo != TipoColumna::STR || !col.objeto || col.nombre == "_original" || col.sn.empty()) continue;
                for (size_t k = 1; k < col.s.size(); ++k)
                    if (col.sn[k] && !col.sn[k - 1]) { col.s[k] = col.s[k - 1]; col.sn[k] = 0; }
            }
        }
    }
    const auto& fs = d.col("fecha_sesion").i;
    const size_t n = d.filas();
    ClasLote cl = clasificar_lote(cal, fs);
    d.poner_str("tipo_dia", cl.tipo);
    d.poner_str("nombre_dia", cl.nombre);
    d.poner_f64("fraccion_sesion", cl.frac);
    std::vector<int64_t> imp(n, 0);
    if (d.tiene("es_imputado")) imp = llenar_entero(d.col("es_imputado"), 1);
    const std::vector<int64_t> orig = d.col("_original").i;
    std::vector<int64_t> hay(n), ses(n), ei(n);
    for (size_t k = 0; k < n; ++k) {
        hay[k] = (orig[k] && imp[k] == 0) ? 1 : 0;
        ses[k] = (hay[k] && es_tipo_sesion_str(cl.tipo[k])) ? 1 : 0;
        ei[k] = (imp[k] | (orig[k] ? 0 : 1)) ? 1 : 0;
    }
    d.poner_bool("hay_dato", std::move(hay));
    d.poner_bool("es_sesion_con_dato", std::move(ses));
    d.poner_str("calendario", std::vector<std::string>(n, cal.nombre()));
    d.poner_i64("es_imputado", std::move(ei));
    const CatalogoUniversoHistorico& cat = c.catalogo_universo ? *c.catalogo_universo : catalogo_base();
    cat.auditar_dataframe_supervivencia(d, sim, "fecha_sesion");
    const std::vector<int64_t> fsd = d.col("fecha_sesion").i;        // copias: los poner_* de abajo realojan las columnas
    const std::vector<int64_t> valida = d.col("es_sesion_con_dato").i;
    std::vector<int64_t> dcal(n, 0), sacum(n, 0);
    int64_t prev = -1;
    for (size_t i = 0; i < n; ++i) {
        if (!valida[i]) continue;
        if (prev >= 0) {
            dcal[i] = std::max<int64_t>(1, piso_div(fsd[i] - fsd[static_cast<size_t>(prev)], NS_DIA));
            const int64_t zi = piso_div(fsd[static_cast<size_t>(prev)] + NS_DIA, NS_DIA), zf = piso_div(fsd[i], NS_DIA);
            sacum[i] = std::max<int64_t>(1, static_cast<int64_t>(sesiones_cal(cal, zi, zf).size()));
        } else {
            dcal[i] = 1;
            sacum[i] = 1;
        }
        prev = static_cast<int64_t>(i);
    }
    d.poner_i64("dias_calendario_acumulados", dcal);
    d.poner_i64("sesiones_acumuladas", std::move(sacum));
    d.poner_i64("dias_acumulados", std::move(dcal));
    std::vector<int64_t> cierre(n);
    for (size_t k = 0; k < n; ++k) cierre[k] = cierre_real_utc_ns(sim, fsd[k], fuente, c.corte_diario_cripto);
    d.poner_ts("cierre_real_utc", std::move(cierre));
    d.quitar("_original");
    return d;
}

// Dias operables del calendario en TODO el rango del contrato mas los que trae la fuente, anotados (solo para la lista).
struct TablaDias { std::vector<int64_t> fecha_sesion; ClasLote clas; std::vector<uint8_t> hay_dato; std::string calendario; };

inline TablaDias tabla_dias_contrato(const ContratoExtraccion& c, const std::string& sim, const Tabla& df) {
    const std::optional<std::string> fuente = fuente_de(df);
    const Calendario cal = calendario_para(sim, fuente);
    const auto [lo, hi] = rango_sesiones_contrato(c);
    const int64_t ayer = normalizar_dia(ahora_ns()) - NS_DIA;
    const auto& ses = df.col("fecha_sesion").i;
    const int64_t hi_g = ses.empty() ? std::min(hi, ayer) : std::max(std::min(hi, ayer), *std::max_element(ses.begin(), ses.end()));
    std::set<int64_t> dias;
    for (int64_t t = lo; t <= hi_g; t += NS_DIA) if (cal.es_dia_operable(piso_div(t, NS_DIA))) dias.insert(t);
    for (int64_t t : ses) dias.insert(t);
    std::set<int64_t> validos;
    if (df.tiene("es_imputado")) {
        const auto imp = llenar_entero(df.col("es_imputado"), 1);
        for (size_t k = 0; k < ses.size(); ++k) if (imp[k] == 0) validos.insert(ses[k]);
    } else {
        validos.insert(ses.begin(), ses.end());
    }
    TablaDias r;
    r.fecha_sesion.assign(dias.begin(), dias.end());
    r.clas = clasificar_lote(cal, r.fecha_sesion);
    for (int64_t t : r.fecha_sesion) r.hay_dato.push_back(validos.count(t) ? 1 : 0);
    r.calendario = cal.nombre();
    return r;
}

inline void dias_no_sesion(const std::string& sim, const TablaDias& d, JV& out) {
    for (size_t k = 0; k < d.fecha_sesion.size(); ++k) {
        const std::string& tipo = d.clas.tipo[k];
        const bool en_sesion = es_tipo_sesion_str(tipo);
        const bool hay = d.hay_dato[k] != 0;
        std::string cat;
        if (tipo == "FERIADO" || tipo == "CIERRE_EXTRAORDINARIO" || tipo == "SESION_REDUCIDA") cat = (en_sesion || !hay) ? tipo : "FUERA_DE_CALENDARIO";
        else if (tipo == "SESION_REGULAR" && !hay) cat = "SESION_SIN_DATO";
        else if (tipo == "FIN_DE_SEMANA" && hay) cat = "FUERA_DE_CALENDARIO";
        else continue;
        JV r = JV::objeto();
        r.poner("simbolo", JV::cadena(sim));
        r.poner("fecha_sesion", JV::cadena(texto_fecha(d.fecha_sesion[k])));
        r.poner("categoria", JV::cadena(cat));
        r.poner("tipo_dia", JV::cadena(tipo));
        r.poner("nombre_dia", JV::cadena(d.clas.nombre[k]));
        r.poner("fraccion_sesion", JV::real(d.clas.frac[k]));
        r.poner("hay_dato", JV::booleano(hay));
        r.poner("calendario", JV::cadena(d.calendario));
        out.push(std::move(r));
    }
}

// ========================================================================================================= procesamiento
struct ExtraProcesado {
    JV dias_no_sesion = JV::lista();
    std::vector<std::pair<std::string, Calendario>> calendarios;
    int64_t filas_excluidas = 0;
    std::vector<std::pair<std::string, std::optional<std::string>>> fuentes;
    JV supervivencia = JV::objeto();
};
struct ResultadoProcesado {
    Tabla df;
    std::vector<IncidenteHueco> incidentes;
    ConteoHuecos conteo;
    bool es_monotono = true;
    std::vector<std::string> advertencias;
    Crudos dfs_limpios;
    ExtraProcesado extra;
};

template <class V>
inline void poner_par(std::vector<std::pair<std::string, V>>& v, const std::string& k, V x) {
    for (auto& kv : v) if (kv.first == k) { kv.second = std::move(x); return; }
    v.emplace_back(k, std::move(x));
}

inline ResultadoProcesado procesar_crudos(const ContratoExtraccion& c, const Crudos& crudos, double delta) {
    ResultadoProcesado R;
    R.conteo = conteo_vacio();
    const bool diario = delta >= 86400.0;
    for (const auto& [sim, df0] : crudos) {
        Tabla df = df0;
        auto [mono, inc_m] = verificar_monotonicidad_incidentes(df.col("fecha_utc").i);
        if (!mono) {
            R.es_monotono = false;
            R.advertencias.push_back("VIOLACIÓN DE MONOTONICIDAD (" + sim + "): Detectadas " + std::to_string(inc_m.size())
                                     + " inversiones temporales o duplicados.");
        }
        R.incidentes.insert(R.incidentes.end(), inc_m.begin(), inc_m.end());
        verificar_estampa_diaria(df);
        if (df.tiene("tipo_precio")) {
            const Columna& tp = df.col("tipo_precio");
            bool aj = false;
            for (size_t k = 0; k < df.filas() && !aj; ++k) aj = !celda_nula(tp, k) && celda_texto(tp, k) == "AJUSTADO_RETROACTIVO";
            if (aj)
                R.advertencias.push_back("PRECIO AJUSTADO RETROACTIVO (" + sim + "): el NIVEL se reescribe hacia atrás con dividendos futuros; "
                                         "los retornos son causales (retorno total), los niveles no. No usar niveles, umbrales de precio ni "
                                         "cocientes de nivel con esta serie.");
        }
        const std::string col_h = (delta >= 86400.0 && df.tiene("fecha_sesion")) ? "fecha_sesion" : "fecha_utc";
        auto [inc_h, cnt] = clasificar_huecos_ontologicos(df, col_h, "precio", c.venue, sim, delta);
        R.incidentes.insert(R.incidentes.end(), inc_h.begin(), inc_h.end());
        for (const auto& kv : cnt)
            for (auto& g : R.conteo) if (g.first == kv.first) g.second += kv.second;
        if (diario && df.tiene("fecha_sesion")) {
            dias_no_sesion(sim, tabla_dias_contrato(c, sim, df), R.extra.dias_no_sesion);
            df = aplicar_rejilla_diaria(c, sim, df);
            std::optional<std::string> f = df.tiene("fuente") ? std::optional<std::string>(celda_texto(df.col("fuente"), 0)) : std::nullopt;
            poner_par(R.extra.fuentes, sim, f);
            poner_par(R.extra.calendarios, sim, calendario_para(sim, f));
        }
        std::vector<std::string> cprecio;
        for (const char* x : {"precio", "open", "high", "low", "close", "best_bid", "best_ask", "mid_price"}) if (df.tiene(x)) cprecio.push_back(x);
        Tabla bl = blindar_serie_con_mascara(df, "fecha_utc", cprecio, !c.crudo);
        const CatalogoUniversoHistorico& cat = c.catalogo_universo ? *c.catalogo_universo : catalogo_base();
        JV info = cat.auditar_dataframe_supervivencia(bl, sim, bl.tiene("fecha_sesion") ? "fecha_sesion" : "fecha_utc");
        const int64_t npre = info.get("filas_pre_listado")->i, npost = info.get("filas_deslistado")->i;
        R.extra.supervivencia.poner(sim, info);
        if (npre > 0)
            R.advertencias.push_back("SUPERVIVENCIA UNIVERSO (" + sim + "): " + std::to_string(npre)
                                     + " filas anteriores a fecha_listado marcadas como PRE_LISTADO (es_sesion_con_dato=False).");
        if (npost > 0)
            R.advertencias.push_back("SUPERVIVENCIA UNIVERSO (" + sim + "): " + std::to_string(npost)
                                     + " filas posteriores a fecha_deslistado marcadas como DESLISTADO (forward-fill cortado a NaN).");
        if (c.variables && !c.variables->empty() && !c.crudo) {
            std::vector<std::string> falt;
            for (const auto& v : *c.variables) if (!bl.tiene(v)) falt.push_back(v);
            if (!falt.empty())
                throw IncompatibilidadContratoError("INCOMPATIBILIDAD DE VARIABLES (Axioma 2): Las variables " + repr_lista_str(falt) + " no existen en el venue '"
                                                    + nombre(c.venue) + "' (nivel '" + nombre(c.nivel) + "'). Disponibles: "
                                                    + repr_lista_str(nombres_columnas(bl)) + ".");
            std::set<std::string> fijas{"fecha_utc", "fecha_sesion", "precio", "es_imputado"};
            for (const auto& a : columnas_anotacion_dia()) fijas.insert(a);
            std::vector<std::string> ret{"fecha_utc"};
            if (bl.tiene("fecha_sesion")) ret.push_back("fecha_sesion");
            for (const auto& v : *c.variables) if (!fijas.count(v)) ret.push_back(v);
            ret.push_back("precio");
            ret.push_back("es_imputado");
            for (const auto& a : columnas_anotacion_dia()) if (bl.tiene(a)) ret.push_back(a);
            Tabla s;
            for (const auto& x : ret) s.columnas.push_back(bl.col(x));
            bl = std::move(s);
        }
        R.dfs_limpios.emplace_back(sim, std::move(bl));
    }
    Alineado al = alinear_multivariado_causal(R.dfs_limpios, c.politica_calendario);
    Tabla& db = al.df;
    if (db.vacia()) {
        throw IncompatibilidadContratoError(std::string("INCOMPATIBILIDAD MULTIVARIADA: La alineación '") + nombre(c.politica_calendario) + "' entre "
                                            + repr_lista_str(c.simbolos) + " produjo una intersección vacía.");
    }
    const int64_t n_excl = al.filas_excluidas < 0 ? 0 : al.filas_excluidas;
    if (diario && R.dfs_limpios.size() > 1 && !R.extra.calendarios.empty()) {
        std::vector<int64_t> fs = db.col("fecha_utc").i;
        for (auto& t : fs) t -= NS_DIA;
        db.poner_ts("fecha_sesion", fs);
        for (const auto& [sim, cal] : R.extra.calendarios) {
            const ClasLote cl = clasificar_lote(cal, fs);
            for (const char* cc : {"tipo_dia", "nombre_dia", "fraccion_sesion"}) {
                const std::string col = std::string(cc) + "_" + sim;
                const bool frac = std::string(cc) == "fraccion_sesion";
                if (db.tiene(col)) {
                    Columna& x = db.col(col);
                    if (frac) { for (size_t k = 0; k < x.f.size(); ++k) if (std::isnan(x.f[k])) x.f[k] = cl.frac[k]; }
                    else if (!x.sn.empty()) {
                        const auto& src = std::string(cc) == "tipo_dia" ? cl.tipo : cl.nombre;
                        for (size_t k = 0; k < x.s.size(); ++k) if (x.sn[k]) { x.s[k] = src[k]; x.sn[k] = 0; }
                    }
                } else if (frac) db.poner_f64(col, cl.frac);
                else db.poner_str(col, std::string(cc) == "tipo_dia" ? cl.tipo : cl.nombre);
            }
            const std::string cn = "calendario_" + sim;
            if (db.tiene(cn)) db.poner_str(cn, std::vector<std::string>(db.filas(), cal.nombre()));
        }
    }
    if (n_excl)
        R.advertencias.push_back("INTERSECCION_ESTRICTA: " + std::to_string(n_excl) + " filas con algún mercado fuera de sesión se separaron de la serie; "
                                 "los días quedan en manifiesto.dias_no_sesion (no se descartan).");
    const std::vector<double> dnu = calcular_incrementos_nu_causales(db, c.tipo_reloj, c.lookback_W);
    db.poner_f64("delta_nu", dnu);
    std::vector<double> acum(dnu.size());
    double s = 0.0;
    for (size_t k = 0; k < dnu.size(); ++k) { s += dnu[k]; acum[k] = s; }
    db.poner_f64("nu_acumulado", std::move(acum));
    R.extra.filas_excluidas = n_excl;
    R.df = std::move(db);
    return R;
}

// ============================================================================================================ cobertura
inline JV cobertura_diaria(const ContratoExtraccion& c, const Crudos& dfs) {
    const auto [lo, hi] = rango_sesiones_contrato(c);
    const int64_t zlo = piso_div(lo, NS_DIA), zhi = piso_div(hi, NS_DIA);
    JV por = JV::objeto();
    std::vector<std::pair<std::string, double>> kappas;
    for (const auto& [sim, d] : dfs) {
        const std::optional<std::string> fuente = fuente_de(d);
        const Calendario cal = calendario_para(sim, fuente);
        const std::vector<int64_t> esp = sesiones_cal(cal, zlo, zhi);
        const std::set<int64_t> sesp(esp.begin(), esp.end());
        std::set<int64_t> val;
        const auto imp = a_dobles(d.col("es_imputado"));
        const auto& fs = d.col("fecha_sesion").i;
        for (size_t k = 0; k < fs.size(); ++k) if (imp[k] == 0.0) val.insert(piso_div(fs[k], NS_DIA));
        int64_t en = 0;
        std::vector<std::string> fuera, faltan;
        for (int64_t z : val) { if (sesp.count(z)) ++en; else fuera.push_back(texto_fecha_z(z)); }
        for (int64_t z : esp) if (!val.count(z)) faltan.push_back(texto_fecha_z(z));
        const double kappa = esp.empty() ? NaN : static_cast<double>(en) / static_cast<double>(esp.size());
        JV r = JV::objeto();
        r.poner("calendario", JV::cadena(cal.nombre()));
        r.poner("calendario_verificado", JV::booleano(verificado_en_z(cal, zlo, zhi)));
        r.poner("esperadas", JV::entero(static_cast<int64_t>(esp.size())));
        r.poner("observadas_en_calendario", JV::entero(en));
        r.poner("kappa", JV::real(kappa));
        r.poner("faltantes", JV::entero(static_cast<int64_t>(esp.size()) - en));
        r.poner("fuera_de_calendario", JV::entero(static_cast<int64_t>(fuera.size())));
        faltan.resize(std::min<size_t>(faltan.size(), 10));
        std::vector<std::string> f10(fuera.begin(), fuera.begin() + static_cast<std::ptrdiff_t>(std::min<size_t>(fuera.size(), 10)));
        r.poner("ejemplos_faltantes", JV::lista_str(faltan));
        r.poner("ejemplos_fuera", JV::lista_str(f10));
        por.poner(sim, std::move(r));
        if (std::isfinite(kappa)) poner_par(kappas, sim, kappa);
    }
    std::string peor = por.o.front().first;
    if (!kappas.empty()) {
        peor = kappas.front().first;
        double m = kappas.front().second;
        for (const auto& [s, k] : kappas) if (k < m) { m = k; peor = s; }
    }
    double kmin = 0.0;
    for (const auto& [s, k] : kappas) if (s == peor) kmin = k;
    JV out = JV::objeto();
    JV rango = JV::lista();
    rango.push(JV::cadena(texto_fecha(lo)));
    rango.push(JV::cadena(texto_fecha(hi)));
    out.poner("rango_sesiones", std::move(rango));
    out.poner("por_simbolo", std::move(por));
    out.poner("simbolo_peor", JV::cadena(peor));
    out.poner("kappa_min", JV::real(kmin));
    return out;
}

// ======================================================================================================== testigo y salida
inline Tabla salida_auditable(const ContratoExtraccion& c, const Crudos& crudos, double delta) {
    Tabla df = procesar_crudos(c, crudos, delta).df;
    if (c.modo == ModoEjecucion::MODO_RODANTE_CAUSAL) {
        const CausalWindowEngine m(df, c.lookback_W, c.stride_S, c.tipo_reloj, c.umbral_nu);
        std::vector<int64_t> marca(df.filas(), 0);
        for (int64_t i : m.indices_disparo) marca[static_cast<size_t>(i)] = 1;
        df.poner_bool("es_disparo", std::move(marca));
    }
    return df;
}

inline JV testigo_fuente(const ContratoExtraccion& c, const Crudos& crudos, int64_t corte) {
    JV rep = JV::objeto();
    rep.poner("corte", JV::cadena(texto_timestamp(corte)));
    rep.poner("max_diff", JV::real(0.0));
    rep.poner("simbolos_distintos", JV::lista());
    Crudos cc;
    try {
        ContratoExtraccion c2 = c;
        c2.fecha_fin = iso_timestamp(corte);
        cc = cargar_crudos(c2).first;
    } catch (const std::exception& e) {
        rep.poner("max_diff", JV::real(std::numeric_limits<double>::infinity()));
        rep.poner("excepcion", JV::cadena(clase_excepcion(e) + ": " + e.what()));
        return rep;
    }
    double md = 0.0;
    for (const auto& [sim, full] : crudos) {
        const Tabla* obs = nullptr;
        for (const auto& [s2, t2] : cc) if (s2 == sim) { obs = &t2; break; }
        Tabla vacia;
        if (!obs) { vacia = tomar_filas(full, {}); obs = &vacia; }
        const RepComparacion r = comparar_tablas_hasta(full, *obs, corte);
        if (r.max_diff > 0.0) {
            JV x = JV::objeto();
            x.poner("simbolo", JV::cadena(sim));
            for (auto& kv : r.a_jv().o) x.poner(kv.first, kv.second);
            rep.get("simbolos_distintos")->push(std::move(x));
            md = std::max(md, r.max_diff);
        }
    }
    rep.poner("max_diff", JV::real(md));
    return rep;
}

struct ResultadoExtraccion {
    Tabla df;                                    // payload (modo estatico) o engine.df (modo rodante)
    std::optional<CausalWindowEngine> motor;     // solo en MODO_RODANTE_CAUSAL
    JV manifiesto;                               // ManifiestoIntegridad.a_dict()
    std::vector<IncidenteHueco> inventario_huecos;
};

inline ResultadoExtraccion ejecutar_extraccion(const ContratoExtraccion& c) {
    validar_contrato(c);
    auto [crudos, delta] = cargar_crudos(c);
    ResultadoProcesado P = procesar_crudos(c, crudos, delta);
    std::vector<std::string>& adv = P.advertencias;
    const int64_t hi_ses = rango_sesiones_contrato(c).second;
    {
        std::vector<std::pair<std::string, Calendario>> unicos;
        for (const auto& [s, cal] : P.extra.calendarios) poner_par(unicos, std::string(cal.nombre()), cal);
        for (const auto& [n, cal] : unicos)
            if (auto a = advertir_auditoria_anual(cal, piso_div(hi_ses, NS_DIA))) adv.push_back(*a);
    }
    JV testigo = auditar_causalidad_pipeline([&](const Crudos& d) { return salida_auditable(c, d, delta); }, crudos, c.n_cortes_testigo,
                                             c.semilla_testigo);
    double score = testigo.get("score")->d;
    if (c.testigo_en_fuente && !testigo.get("cortes")->a.empty()) {
        std::vector<std::string> cs;
        for (const auto& x : testigo.get("cortes")->a) cs.push_back(x.s);
        std::sort(cs.begin(), cs.end());
        const int64_t corte = ns_desde_iso(cs[cs.size() / 2]);
        JV tf = testigo_fuente(c, crudos, corte);
        const double md = tf.get("max_diff")->d;
        testigo.poner("fuente", std::move(tf));
        if (!std::isnan(score)) score = std::max(score, md);
    }
    if (std::isnan(score)) {
        const JV* av = testigo.get("aviso");
        adv.push_back("TESTIGO NO CERTIFICABLE: " + (av ? av->s : std::string("sin cortes")));
    } else if (score > 0.0) {
        adv.push_back("FUGA DE FUTURO DETECTADA por el testigo (score=" + repr_py(score) + "): ver manifiesto.detalle_testigo.");
    }
    JV det_cob = JV::objeto();
    int64_t esperados = 0, reales = 0;
    double cob = 0.0;
    bool todos_ses = true;
    for (const auto& [s, d] : P.dfs_limpios) todos_ses = todos_ses && d.tiene("fecha_sesion");
    if (delta >= 86400.0 && todos_ses) {
        det_cob = cobertura_diaria(c, P.dfs_limpios);
        const std::string peor = det_cob.get("simbolo_peor")->s;
        const JV* ps = det_cob.get("por_simbolo")->get(peor);
        esperados = ps->get("esperadas")->i;
        reales = ps->get("observadas_en_calendario")->i;
        cob = 100.0 * det_cob.get("kappa_min")->d;
        for (const auto& [sim, d] : det_cob.get("por_simbolo")->o) {
            if (!d.get("calendario_verificado")->b)
                adv.push_back("CALENDARIO NO VERIFICADO (" + sim + "): '" + d.get("calendario")->s
                              + "' no está verificado en el rango pedido; la cobertura es aproximada.");
            if (d.get("fuera_de_calendario")->i > 0) {
                std::vector<std::string> ej;
                for (size_t k = 0; k < d.get("ejemplos_fuera")->a.size() && k < 3; ++k) ej.push_back(d.get("ejemplos_fuera")->a[k].s);
                adv.push_back("SESIONES FUERA DE CALENDARIO (" + sim + "): " + std::to_string(d.get("fuera_de_calendario")->i)
                              + " sesiones con dato en días que '" + d.get("calendario")->s + "' da por cerrados (ej. " + repr_lista_str(ej) + ").");
            }
        }
    } else {
        const int64_t ti = ns_desde_iso(recortar(c.fecha_inicio)), tf = ns_desde_iso(recortar(c.fecha_fin));
        const double seg = std::max(1.0, static_cast<double>(tf - ti) / 1e9);
        const double dias = seg / 86400.0;
        const double brutos = seg / std::max(1.0, delta);
        if ((c.venue == Venue::YAHOO || c.venue == Venue::BLOOMBERG_XLSX) && dias >= 5.0) {
            if (delta >= 86400.0) esperados = std::max<int64_t>(1, static_cast<int64_t>(brutos * (252.0 / 365.0)));
            else esperados = std::max<int64_t>(1, static_cast<int64_t>(brutos * (5.0 / 7.0)));
        } else {
            esperados = std::max<int64_t>(1, static_cast<int64_t>(brutos));
        }
        reales = static_cast<int64_t>(P.df.filas());
        cob = std::min(100.0, (static_cast<double>(reales) / static_cast<double>(std::max<int64_t>(1, esperados))) * 100.0);
    }
    JV resumen_h = JV::objeto(), admite = JV::objeto(), tipo_p = JV::objeto();
    for (const auto& sim : c.simbolos) {
        const std::string su = mayusculas(sim);
        bool aj = false;
        auto primero_aj = [&](const std::string& col) {
            return P.df.filas() > 0 && !celda_nula(P.df.col(col), 0) && celda_texto(P.df.col(col), 0) == "AJUSTADO_RETROACTIVO";
        };
        if (su == "VNQ") aj = true;
        else if (P.df.tiene("tipo_precio_" + sim)) aj = primero_aj("tipo_precio_" + sim);
        else if (P.df.tiene("tipo_precio")) aj = primero_aj("tipo_precio");
        if (aj) {
            admite.poner(su, JV::booleano(false));
            tipo_p.poner(su, JV::cadena("AJUSTADO_RETROACTIVO"));
            adv.push_back("Activo '" + su + "' tiene tipo_precio='AJUSTADO_RETROACTIVO'. Vetado para estrategias de nivel nominal.");
        } else {
            admite.poner(su, JV::booleano(true));
            tipo_p.poner(su, JV::cadena("CRUDO"));
        }
    }
    if (c.auditar_libro_hashes && !P.df.vacia()) {
        LibroHashesInmutable libro(c.ledger_hashes);
        for (const auto& sim : c.simbolos) {
            const std::string colp = P.df.tiene("precio_" + sim) ? "precio_" + sim : (P.df.tiene("precio") ? std::string("precio") : std::string());
            if (colp.empty()) continue;
            Tabla sub;
            Columna f = P.df.col(P.df.tiene("fecha_sesion") ? "fecha_sesion" : "fecha_utc");
            f.nombre = "fecha_sesion";
            sub.columnas.push_back(std::move(f));
            Columna p = P.df.col(colp);
            p.nombre = "precio";
            sub.columnas.push_back(std::move(p));
            JV r = libro.auditar_y_registrar_barras(sub, sim, nombre(c.venue));
            if (!r.get("es_inmutable")->b)
                adv.push_back("REVISIÓN RETROSPECTIVA DETECTADA en libro de hashes para " + sim + ": " + std::to_string(r.get("revisiones_detectadas")->i)
                              + " barras alteradas por proveedor.");
            resumen_h.poner(sim, std::move(r));
        }
    }
    const bool aprobado = P.es_monotono && score == 0.0 && (100.0 - cob) <= c.tolerancia_huecos_pct;
    JV m = JV::objeto();
    m.poner("puntos_esperados", JV::entero(esperados));
    m.poner("puntos_reales", JV::entero(reales));
    m.poner("tasa_cobertura_pct", JV::real(round_py(cob, 4)));
    m.poner("es_monotono", JV::booleano(P.es_monotono));
    m.poner("lookahead_score", JV::real(round_py(score, 6)));
    m.poner("total_huecos", JV::entero(static_cast<int64_t>(P.incidentes.size())));
    JV conteo = JV::objeto();
    for (const auto& [k, v] : P.conteo) conteo.poner(k, JV::entero(v));
    m.poner("conteo_por_tipo", std::move(conteo));
    m.poner("advertencias", JV::lista_str(adv));
    m.poner("aprobado_para_modelo", JV::booleano(aprobado));
    m.poner("detalle_testigo", std::move(testigo));
    m.poner("detalle_cobertura", std::move(det_cob));
    m.poner("dias_no_sesion", P.extra.dias_no_sesion);
    m.poner("ventanas_informacion", P.extra.fuentes.empty() ? JV::lista() : tabla_manifiesto_ventanas(P.extra.fuentes, c.corte_diario_cripto));
    m.poner("libro_hashes", std::move(resumen_h));
    m.poner("admite_estrategias_de_nivel", std::move(admite));
    m.poner("tipo_precio", std::move(tipo_p));
    m.poner("supervivencia_universo", P.extra.supervivencia);
    ResultadoExtraccion out;
    out.manifiesto = std::move(m);
    out.inventario_huecos = P.incidentes;
    if (c.modo == ModoEjecucion::MODO_RODANTE_CAUSAL) {
        out.motor.emplace(P.df, c.lookback_W, c.stride_S, c.tipo_reloj, c.umbral_nu);
        out.df = out.motor->df;
    } else {
        out.df = std::move(P.df);
    }
    return out;
}

}  // namespace ark::dataset::ingesta
