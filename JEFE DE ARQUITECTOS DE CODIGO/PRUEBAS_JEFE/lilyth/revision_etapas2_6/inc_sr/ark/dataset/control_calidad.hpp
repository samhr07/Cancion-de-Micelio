// control_calidad.hpp -- CALCA DE contratos.py, control_calidad.py, catalogo_universo.py y libro_hashes.py (ARK-DATA-02-CXX)
// Reactor ARK / Cancion del Micelio -- Jefe de Arquitectos de Codigo, 2026-11-02
//
// Porte tal cual del Python de ARQUITECTO DE DATOS/recopilacion_unificada (orden del Comandante: nada que discutir, solo
// migrar). Reglas, umbrales y mensajes son los del Python; donde pandas o numpy deciden un detalle numerico se copia su
// algoritmo (suma por pares de numpy, mediana, desviacion tipica, generador PCG64 en numpy_rng.hpp).
//   - Formato: repr de float de Python, repr de cadenas y listas, str(pd.Timestamp), isoformat, str(pd.Timedelta).
//   - JV: valor JSON ordenado y json.dumps de Python (separadores ', ' y ': ', ensure_ascii, NaN e Infinity).
//   - contratos.py: enumeraciones, MATRIZ_COMPATIBILIDAD, IncidenteHueco y ContratoExtraccion.
//   - control_calidad.py: RELOJES_PERMITIDOS_POR_VENUE, validar_contrato, clasificar_huecos_ontologicos,
//     blindar_serie_con_mascara, comparar_tablas_hasta, perturbar_futuro y auditar_causalidad_pipeline.
//   - catalogo_universo.py: CATALOGO_HISTORICO_BASE y auditar_dataframe_supervivencia.
//   - libro_hashes.py: computar_hash_barra y LibroHashesInmutable (mismo JSONL que el Python).
#pragma once
#include "ark/dataset/adaptadores_fuentes.hpp"
#include "ark/dataset/guardas_causales.hpp"
#include "ark/dataset/libro_hashes.hpp"
#include "ark/dataset/numpy_rng.hpp"
#include <charconv>
#include <chrono>
#include <functional>
#include <optional>
#include <set>
#include <utility>

namespace ark::dataset::ingesta {

// =========================================================================================================== formato Python
// repr(float) de Python (y str de np.float64): la cifra mas corta que vuelve al mismo double; notacion fija si el exponente
// decimal esta en [-4, 16), cientifica en otro caso; '.0' en los enteros de notacion fija.
inline std::string repr_py(double x) {
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    char b[64];
    const auto r = std::to_chars(b, b + sizeof b, x, std::chars_format::scientific);
    std::string t(b, r.ptr);
    const bool neg = !t.empty() && t[0] == '-';
    if (neg) t.erase(0, 1);
    const size_t pe = t.find('e');
    std::string mant = t.substr(0, pe);
    const int e10 = std::stoi(t.substr(pe + 1));
    std::string dig;
    for (char c : mant) if (c != '.') dig += c;
    while (dig.size() > 1 && dig.back() == '0') dig.pop_back();
    std::string out;
    if (e10 >= -4 && e10 < 16) {
        if (e10 >= 0) {
            const size_t ent = static_cast<size_t>(e10) + 1;
            if (dig.size() <= ent) out = dig + std::string(ent - dig.size(), '0') + ".0";
            else out = dig.substr(0, ent) + "." + dig.substr(ent);
        } else {
            out = "0." + std::string(static_cast<size_t>(-e10 - 1), '0') + dig;
        }
    } else {
        out = dig.substr(0, 1);
        if (dig.size() > 1) out += "." + dig.substr(1);
        char e[16];
        std::snprintf(e, sizeof e, "e%c%02d", e10 < 0 ? '-' : '+', e10 < 0 ? -e10 : e10);
        out += e;
    }
    return neg ? "-" + out : out;
}

// f"{x:.Nf}" de Python (redondeo correcto, como printf de la UCRT).
inline std::string fmt_f(double x, int dec) {
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    char b[512];
    std::snprintf(b, sizeof b, "%.*f", dec, x);
    return b;
}

// round(x, n) de Python sobre un float.
inline double round_py(double x, int n) {
    if (!std::isfinite(x)) return x;
    char b[512];
    std::snprintf(b, sizeof b, "%.*f", n, x);
    return std::strtod(b, nullptr);
}

// repr(str) de Python (comillas simples salvo que la cadena tenga ' y no ").
inline std::string repr_str_py(const std::string& s) {
    const bool usa_dobles = s.find('\'') != std::string::npos && s.find('"') == std::string::npos;
    const char q = usa_dobles ? '"' : '\'';
    std::string r(1, q);
    for (unsigned char c : s) {
        if (c == '\\') r += "\\\\";
        else if (c == static_cast<unsigned char>(q)) { r += '\\'; r += static_cast<char>(c); }
        else if (c == '\n') r += "\\n";
        else if (c == '\r') r += "\\r";
        else if (c == '\t') r += "\\t";
        else if (c < 0x20 || c == 0x7f) { char h[8]; std::snprintf(h, sizeof h, "\\x%02x", c); r += h; }
        else r += static_cast<char>(c);
    }
    return r + q;
}

inline std::string repr_lista_str(const std::vector<std::string>& v) {
    std::string r = "[";
    for (size_t k = 0; k < v.size(); ++k) { if (k) r += ", "; r += repr_str_py(v[k]); }
    return r + "]";
}

// str(pd.Timestamp(ns).date()) en UTC.
inline std::string texto_fecha(int64_t ns) {
    const FechaCivil c = civil_desde_dias(piso_div(ns, NS_DIA));
    char b[32];
    std::snprintf(b, sizeof b, "%04d-%02d-%02d", c.y, c.m, c.d);
    return b;
}
inline std::string texto_fecha_z(int64_t z) { return texto_fecha(z * NS_DIA); }

// pd.Timestamp.isoformat(): 'YYYY-MM-DDTHH:MM:SS[.ffffff|.fffffffff]+00:00'.
inline std::string iso_timestamp(int64_t ns) {
    std::string t = texto_timestamp(ns);
    if (t.size() > 10) t[10] = 'T';
    return t;
}

// str(pd.Timedelta(ns)) para un desplazamiento no negativo: 'D days HH:MM:SS[.fffffffff]'.
inline std::string texto_timedelta(int64_t ns) {
    const int64_t d = ns / NS_DIA, r = ns % NS_DIA;
    const int64_t s = r / NS_S, f = r % NS_S;
    char b[96];
    std::snprintf(b, sizeof b, "%lld days %02lld:%02lld:%02lld", static_cast<long long>(d), static_cast<long long>(s / 3600),
                  static_cast<long long>((s / 60) % 60), static_cast<long long>(s % 60));
    std::string out = b;
    if (f) {
        char x[16];
        if (f % 1000 == 0) std::snprintf(x, sizeof x, ".%06lld", static_cast<long long>(f / 1000));
        else std::snprintf(x, sizeof x, ".%09lld", static_cast<long long>(f));
        out += x;
    }
    return out;
}

inline int64_t ahora_ns() {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}
// pd.Timestamp.now(tz='UTC').isoformat(): microsegundos.
inline std::string ahora_iso() {
    const int64_t ns = ahora_ns() / 1000 * 1000;
    std::string t = iso_timestamp(ns);
    if (t.find('.') == std::string::npos) t.insert(19, ".000000");
    return t;
}

// ============================================================================================= JV: json.dumps de Python
struct JV {
    enum class T : uint8_t { NUL, BOOL, INT, DBL, STR, ARR, OBJ } t = T::NUL;
    bool b = false;
    int64_t i = 0;
    double d = 0.0;
    std::string s;
    std::vector<JV> a;
    std::vector<std::pair<std::string, JV>> o;

    static JV nulo() { return JV{}; }
    static JV booleano(bool v) { JV x; x.t = T::BOOL; x.b = v; return x; }
    static JV entero(int64_t v) { JV x; x.t = T::INT; x.i = v; return x; }
    static JV real(double v) { JV x; x.t = T::DBL; x.d = v; return x; }
    static JV cadena(std::string v) { JV x; x.t = T::STR; x.s = std::move(v); return x; }
    static JV lista() { JV x; x.t = T::ARR; return x; }
    static JV objeto() { JV x; x.t = T::OBJ; return x; }
    static JV lista_str(const std::vector<std::string>& v) { JV x = lista(); for (const auto& e : v) x.a.push_back(cadena(e)); return x; }

    JV& push(JV v) { a.push_back(std::move(v)); return a.back(); }
    // d[k] = v: sustituye en su sitio si existe (dict de Python) o anade al final.
    JV& poner(const std::string& k, JV v) {
        for (auto& kv : o) if (kv.first == k) { kv.second = std::move(v); return kv.second; }
        o.emplace_back(k, std::move(v));
        return o.back().second;
    }
    JV* get(const std::string& k) { for (auto& kv : o) if (kv.first == k) return &kv.second; return nullptr; }
    const JV* get(const std::string& k) const { for (const auto& kv : o) if (kv.first == k) return &kv.second; return nullptr; }
};

namespace detalle_jv {
inline void cadena(std::string& out, const std::string& s) {
    out += '"';
    for (size_t k = 0; k < s.size();) {
        const unsigned char c = static_cast<unsigned char>(s[k]);
        if (c < 0x80) {
            if (c == '"') out += "\\\"";
            else if (c == '\\') out += "\\\\";
            else if (c == '\n') out += "\\n";
            else if (c == '\r') out += "\\r";
            else if (c == '\t') out += "\\t";
            else if (c == '\b') out += "\\b";
            else if (c == '\f') out += "\\f";
            else if (c < 0x20) { char h[8]; std::snprintf(h, sizeof h, "\\u%04x", c); out += h; }
            else out += static_cast<char>(c);
            ++k;
            continue;
        }
        uint32_t cp = 0;
        size_t n = 0;
        if ((c & 0xE0) == 0xC0) { cp = c & 0x1F; n = 2; }
        else if ((c & 0xF0) == 0xE0) { cp = c & 0x0F; n = 3; }
        else { cp = c & 0x07; n = 4; }
        for (size_t j = 1; j < n && k + j < s.size(); ++j) cp = (cp << 6) | (static_cast<unsigned char>(s[k + j]) & 0x3F);
        k += n;
        char h[16];
        if (cp >= 0x10000) {
            const uint32_t v = cp - 0x10000;
            std::snprintf(h, sizeof h, "\\u%04x\\u%04x", 0xD800 + (v >> 10), 0xDC00 + (v & 0x3FF));
        } else {
            std::snprintf(h, sizeof h, "\\u%04x", cp);
        }
        out += h;
    }
    out += '"';
}
inline void volcar(std::string& out, const JV& v) {
    switch (v.t) {
        case JV::T::NUL: out += "null"; break;
        case JV::T::BOOL: out += v.b ? "true" : "false"; break;
        case JV::T::INT: out += std::to_string(v.i); break;
        case JV::T::DBL:
            if (std::isnan(v.d)) out += "NaN";
            else if (std::isinf(v.d)) out += v.d > 0 ? "Infinity" : "-Infinity";
            else out += repr_py(v.d);
            break;
        case JV::T::STR: cadena(out, v.s); break;
        case JV::T::ARR:
            out += '[';
            for (size_t k = 0; k < v.a.size(); ++k) { if (k) out += ", "; volcar(out, v.a[k]); }
            out += ']';
            break;
        case JV::T::OBJ:
            out += '{';
            for (size_t k = 0; k < v.o.size(); ++k) {
                if (k) out += ", ";
                cadena(out, v.o[k].first);
                out += ": ";
                volcar(out, v.o[k].second);
            }
            out += '}';
            break;
    }
}
}  // namespace detalle_jv

inline std::string json_dumps(const JV& v) { std::string s; detalle_jv::volcar(s, v); return s; }

// ================================================================================================================ contratos.py
enum class Venue : uint8_t { YAHOO, BLOOMBERG_XLSX, BINANCE_UM, BITFINEX, INTERACTIVE_BROKERS, IEX, DUKASCOPY };
enum class ModoEjecucion : uint8_t { MODO_ESTATICO_EXPLORATORIO, MODO_RODANTE_CAUSAL };
enum class TipoReloj : uint8_t { PARED_SEGUNDOS, PARED_MINUTOS, PARED_HORAS, PARED_DIAS, NU_TRADES, NU_VOLUMEN_BASE, NU_VOLUMEN_NOCIONAL,
                                 NU_DELTA_INFO };
enum class PoliticaCalendario : uint8_t { INTERSECCION_ESTRICTA, UNION_CON_MASCARA };
enum class RejillaDiaria : uint8_t { SOLO_SESIONES, SESIONES_Y_FERIADOS, CALENDARIO_COMPLETO };
enum class ModoUniverso : uint8_t { ESTATICO, PUNTO_EN_EL_TIEMPO };
enum class TipoHueco : uint8_t { ESTRUCTURAL_CALENDARIO, ILIQUIDEZ_TRANSITORIA, CAIDA_RED_SOCKET, CENSURA_IZQUIERDA, CORRUPCION_VALOR };
using CorteDiario = ark::dataset::CorteDiarioCripto;   // UTC_00 / NY_16

inline const char* nombre(Venue v) {
    static const char* n[] = {"YAHOO", "BLOOMBERG_XLSX", "BINANCE_UM", "BITFINEX", "INTERACTIVE_BROKERS", "IEX", "DUKASCOPY"};
    return n[static_cast<int>(v)];
}
inline const char* nombre(NivelDato v) {
    static const char* n[] = {"L1_OHLCV", "L1_ENRIQUECIDO", "L2_DEPTH", "L3_MBO", "TOPS"};
    return n[static_cast<int>(v)];
}
inline const char* nombre(ModoEjecucion v) { return v == ModoEjecucion::MODO_ESTATICO_EXPLORATORIO ? "MODO_ESTATICO_EXPLORATORIO" : "MODO_RODANTE_CAUSAL"; }
inline const char* nombre(TipoReloj v) {
    static const char* n[] = {"PARED_SEGUNDOS", "PARED_MINUTOS", "PARED_HORAS", "PARED_DIAS", "NU_TRADES", "NU_VOLUMEN_BASE",
                              "NU_VOLUMEN_NOCIONAL", "NU_DELTA_INFO"};
    return n[static_cast<int>(v)];
}
inline const char* nombre(PoliticaCalendario v) { return v == PoliticaCalendario::INTERSECCION_ESTRICTA ? "INTERSECCION_ESTRICTA" : "UNION_CON_MASCARA"; }
inline const char* nombre(RejillaDiaria v) {
    static const char* n[] = {"SOLO_SESIONES", "SESIONES_Y_FERIADOS", "CALENDARIO_COMPLETO"};
    return n[static_cast<int>(v)];
}
inline const char* nombre(TipoHueco v) {
    static const char* n[] = {"ESTRUCTURAL_CALENDARIO", "ILIQUIDEZ_TRANSITORIA", "CAIDA_RED_SOCKET", "CENSURA_IZQUIERDA", "CORRUPCION_VALOR"};
    return n[static_cast<int>(v)];
}
inline const char* nombre(CorteDiario v) { return v == CorteDiario::NY_16 ? "NY_16" : "UTC_00"; }
inline bool es_pared(TipoReloj r) { return static_cast<int>(r) <= static_cast<int>(TipoReloj::PARED_DIAS); }

// MATRIZ_COMPATIBILIDAD: Venue -> niveles permitidos.
inline std::vector<NivelDato> niveles_validos(Venue v) {
    using N = NivelDato;
    switch (v) {
        case Venue::YAHOO: return {N::L1_OHLCV};
        case Venue::BLOOMBERG_XLSX: return {N::L1_OHLCV};
        case Venue::BINANCE_UM: return {N::L1_OHLCV, N::L1_ENRIQUECIDO};
        case Venue::BITFINEX: return {N::L1_OHLCV, N::L2_DEPTH, N::L3_MBO};
        case Venue::INTERACTIVE_BROKERS: return {N::L1_OHLCV, N::L2_DEPTH};
        case Venue::IEX: return {N::TOPS};
        default: return {N::L1_OHLCV, N::TOPS};
    }
}

// control_calidad.RELOJES_PERMITIDOS_POR_VENUE (None para los venues sin entrada).
inline std::optional<std::vector<TipoReloj>> relojes_validos(Venue v) {
    using R = TipoReloj;
    switch (v) {
        case Venue::YAHOO: return std::vector<R>{R::PARED_DIAS, R::NU_VOLUMEN_BASE, R::NU_VOLUMEN_NOCIONAL, R::NU_DELTA_INFO};
        case Venue::BLOOMBERG_XLSX:
            return std::vector<R>{R::PARED_MINUTOS, R::PARED_HORAS, R::PARED_DIAS, R::NU_VOLUMEN_BASE, R::NU_VOLUMEN_NOCIONAL, R::NU_DELTA_INFO};
        case Venue::BINANCE_UM:
            return std::vector<R>{R::PARED_MINUTOS, R::PARED_HORAS, R::PARED_DIAS, R::NU_TRADES, R::NU_VOLUMEN_BASE, R::NU_VOLUMEN_NOCIONAL,
                                  R::NU_DELTA_INFO};
        case Venue::BITFINEX: return std::vector<R>{R::PARED_SEGUNDOS, R::PARED_MINUTOS, R::NU_DELTA_INFO};
        default: return std::nullopt;
    }
}

struct IncidenteHueco {
    std::string inicio_utc, fin_utc;
    double duracion_segundos = 0.0;
    TipoHueco tipo = TipoHueco::CORRUPCION_VALOR;
    std::string descripcion, activo;
    JV a_jv() const {
        JV j = JV::objeto();
        j.poner("inicio_utc", JV::cadena(inicio_utc));
        j.poner("fin_utc", JV::cadena(fin_utc));
        j.poner("duracion_segundos", JV::real(duracion_segundos));
        j.poner("tipo", JV::cadena(nombre(tipo)));
        j.poner("descripcion", JV::cadena(descripcion));
        j.poner("activo", JV::cadena(activo));
        return j;
    }
};

class CatalogoUniversoHistorico;   // catalogo_universo.py, mas abajo

struct ContratoExtraccion {
    Venue venue = Venue::YAHOO;
    NivelDato nivel = NivelDato::L1_OHLCV;
    std::vector<std::string> simbolos;
    std::string fecha_inicio, fecha_fin;
    ModoEjecucion modo = ModoEjecucion::MODO_ESTATICO_EXPLORATORIO;
    TipoReloj tipo_reloj = TipoReloj::PARED_DIAS;
    PoliticaCalendario politica_calendario = PoliticaCalendario::UNION_CON_MASCARA;
    std::optional<std::vector<std::string>> variables;
    int64_t lookback_W = 60;
    int64_t stride_S = 1;
    std::optional<double> umbral_nu;
    double tolerancia_huecos_pct = 5.0;
    bool crudo = false;
    RejillaDiaria rejilla_diaria = RejillaDiaria::SESIONES_Y_FERIADOS;
    CorteDiario corte_diario_cripto = CorteDiario::UTC_00;
    int64_t n_cortes_testigo = 3;
    uint64_t semilla_testigo = 20261031;
    bool testigo_en_fuente = true;
    bool auditar_libro_hashes = true;
    ModoUniverso modo_universo = ModoUniverso::ESTATICO;
    const CatalogoUniversoHistorico* catalogo_universo = nullptr;
    // Solo C++: el libro de hashes (vacio = <raiz>/telemetria/ledger_ingesta_hashes.jsonl, el LEDGER_DEFECTO de Python).
    std::string ledger_hashes;
};

// =========================================================================================== numpy: sumas y estadisticos
// pairwise_sum de numpy (bloque 128, desenrollado de 8); np.sum = 0.0 + pairwise_sum(todo) (medido: identico en 20 000
// sumas aleatorias).
inline double suma_pares(const double* a, size_t n) {
    if (n < 8) {
        double r = 0.0;
        for (size_t k = 0; k < n; ++k) r += a[k];
        return r;
    }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; ++j) r[j] = a[j];
        size_t k = 8;
        for (; k < n - (n % 8); k += 8)
            for (int j = 0; j < 8; ++j) r[j] += a[k + j];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; k < n; ++k) res += a[k];
        return res;
    }
    size_t n2 = n / 2;
    n2 -= n2 % 8;
    return suma_pares(a, n2) + suma_pares(a + n2, n - n2);
}
inline double np_sum(const std::vector<double>& v) { return 0.0 + suma_pares(v.data(), v.size()); }
inline double np_mean(const std::vector<double>& v) { return np_sum(v) / static_cast<double>(v.size()); }
// np.std(ddof=0) = sqrt(sum((x - media) * (x - media)) / n), con la misma suma por pares.
inline double np_std(const std::vector<double>& v) {
    const double mu = np_mean(v);
    std::vector<double> d(v.size());
    for (size_t k = 0; k < v.size(); ++k) d[k] = (v[k] - mu) * (v[k] - mu);
    return std::sqrt(np_sum(d) / static_cast<double>(v.size()));
}
// np.median: NaN si hay NaN; impar = el central; par = media de los dos centrales.
inline double np_median(std::vector<double> v) {
    for (double x : v) if (std::isnan(x)) return NaN;
    const size_t n = v.size();
    std::sort(v.begin(), v.end());
    if (n % 2) return v[n / 2];
    return (0.0 + v[n / 2 - 1] + v[n / 2]) / 2.0;
}
// max(a, b, c) y min(a, b) de Python sobre floats (con NaN el primero gana si nadie le es estrictamente mayor/menor).
inline double py_max(std::initializer_list<double> xs) { double m = *xs.begin(); for (double x : xs) if (x > m) m = x; return m; }
inline double py_min(double a, double b) { return b < a ? b : a; }

// ============================================================================================ columnas (pandas -> double)
// to_numpy(dtype=float) / pd.to_numeric(errors='coerce').values.astype(float).
inline std::vector<double> a_dobles(const Columna& c) {
    std::vector<double> v;
    if (c.tipo == TipoColumna::F64) return c.f;
    if (c.tipo == TipoColumna::I64 || c.tipo == TipoColumna::BOOL) {
        v.reserve(c.i.size());
        for (int64_t x : c.i) v.push_back(static_cast<double>(x));
        return v;
    }
    if (c.tipo == TipoColumna::TS) {
        for (int64_t x : c.i) v.push_back(x == NAT ? NaN : static_cast<double>(x));
        return v;
    }
    for (size_t k = 0; k < c.s.size(); ++k) {
        double x = NaN;
        if (!c.es_nulo_str(k) && !pandas_xstrtod(c.s[k], x)) x = NaN;
        v.push_back(x);
    }
    return v;
}
inline bool es_numerica(const Columna& c) { return c.tipo == TipoColumna::F64 || c.tipo == TipoColumna::I64 || c.tipo == TipoColumna::BOOL; }

// Valor de celda como lo ve `==` de pandas en columnas no numericas (texto canonico; nulo aparte).
inline bool celda_nula(const Columna& c, size_t k) {
    switch (c.tipo) {
        case TipoColumna::F64: return std::isnan(c.f[k]);
        case TipoColumna::TS: return c.i[k] == NAT;
        case TipoColumna::STR: return c.es_nulo_str(k);
        default: return false;
    }
}
inline std::string celda_texto(const Columna& c, size_t k) {
    switch (c.tipo) {
        case TipoColumna::F64: return repr_py(c.f[k]);
        case TipoColumna::I64: return std::to_string(c.i[k]);
        case TipoColumna::BOOL: return c.i[k] ? "True" : "False";
        case TipoColumna::TS: return texto_timestamp(c.i[k]);
        default: return c.es_nulo_str(k) ? std::string("nan") : c.s[k];
    }
}

inline std::vector<std::string> nombres_columnas(const Tabla& t) {
    std::vector<std::string> v;
    for (const auto& c : t.columnas) v.push_back(c.nombre);
    return v;
}

// Copia de las filas `idx` (en ese orden) de una tabla.
inline Tabla tomar_filas(const Tabla& t, const std::vector<size_t>& idx) {
    Tabla r;
    for (const auto& c : t.columnas) {
        Columna n;
        n.nombre = c.nombre;
        n.tipo = c.tipo;
        n.objeto = c.objeto;
        if (c.tipo == TipoColumna::F64) { n.f.reserve(idx.size()); for (size_t k : idx) n.f.push_back(c.f[k]); }
        else if (c.tipo == TipoColumna::STR) {
            n.s.reserve(idx.size());
            for (size_t k : idx) n.s.push_back(c.s[k]);
            if (!c.sn.empty()) { n.sn.reserve(idx.size()); for (size_t k : idx) n.sn.push_back(c.sn[k]); }
        } else { n.i.reserve(idx.size()); for (size_t k : idx) n.i.push_back(c.i[k]); }
        r.columnas.push_back(std::move(n));
    }
    return r;
}

// df[df[col] <= t] (columna TS).
inline Tabla filas_hasta(const Tabla& t, int64_t corte, const std::string& col = "fecha_utc") {
    std::vector<size_t> idx;
    const auto& v = t.col(col).i;
    for (size_t k = 0; k < v.size(); ++k) if (v[k] != NAT && v[k] <= corte) idx.push_back(k);
    return tomar_filas(t, idx);
}

// ====================================================================================================== control_calidad.py
inline void validar_contrato(const ContratoExtraccion& c) {
    const auto nv = niveles_validos(c.venue);
    if (std::find(nv.begin(), nv.end(), c.nivel) == nv.end()) {
        std::vector<std::string> n;
        for (auto x : nv) n.push_back(nombre(x));
        throw IncompatibilidadContratoError(std::string("INCOMPATIBILIDAD FÍSICA DE CONTRATO: El venue '") + nombre(c.venue)
                                            + "' NO soporta el nivel de datos '" + nombre(c.nivel) + "'. Niveles válidos soportados: "
                                            + repr_lista_str(n) + ".");
    }
    const auto rv = relojes_validos(c.venue);
    if (rv && std::find(rv->begin(), rv->end(), c.tipo_reloj) == rv->end()) {
        std::vector<std::string> n;
        for (auto x : *rv) n.push_back(nombre(x));
        throw IncompatibilidadContratoError(std::string("INCOMPATIBILIDAD DE RESOLUCIÓN TEMPORAL: El venue '") + nombre(c.venue)
                                            + "' NO soporta el reloj '" + nombre(c.tipo_reloj) + "'. Relojes soportados: " + repr_lista_str(n) + ".");
    }
    if (c.corte_diario_cripto == CorteDiario::NY_16 && (c.venue != Venue::BINANCE_UM || c.tipo_reloj != TipoReloj::PARED_DIAS))
        throw IncompatibilidadContratoError("INCOMPATIBILIDAD DE CORTE DIARIO: corte_diario_cripto=NY_16 exige venue BINANCE_UM y reloj PARED_DIAS "
                                            "(la barra se construye desde el minuto de Binance).");
    if (c.simbolos.empty()) throw IncompatibilidadContratoError("La lista de símbolos solicitados no puede estar vacía.");
    for (const auto& s : c.simbolos)
        if (recortar(s).empty()) throw IncompatibilidadContratoError("Símbolo inválido o vacío en contrato: " + repr_str_py(s));
    if (c.lookback_W < 1) throw IncompatibilidadContratoError("lookback_W (" + std::to_string(c.lookback_W) + ") debe ser >= 1.");
    if (c.stride_S < 1) throw IncompatibilidadContratoError("stride_S (" + std::to_string(c.stride_S) + ") debe ser >= 1.");
    if (c.umbral_nu && (!std::isfinite(*c.umbral_nu) || *c.umbral_nu <= 0.0))
        throw IncompatibilidadContratoError("umbral_nu (" + repr_py(*c.umbral_nu) + ") debe ser finito y > 0.0.");
    int64_t ti = 0, tf = 0;
    try {
        ti = ns_desde_iso(recortar(c.fecha_inicio));
        tf = ns_desde_iso(recortar(c.fecha_fin));
    } catch (const std::exception& e) {
        throw IncompatibilidadContratoError(std::string("Formato de fecha inválido (requiere ISO-8601): ") + e.what());
    }
    if (ti >= tf)
        throw IncompatibilidadContratoError("Rango temporal inválido: fecha_inicio (" + c.fecha_inicio + ") debe ser estrictamente anterior a fecha_fin ("
                                            + c.fecha_fin + ").");
}

// verificar_monotonicidad -> IncidenteHueco (tipo CORRUPCION_VALOR, activo GLOBAL_TIME_INDEX).
inline std::pair<bool, std::vector<IncidenteHueco>> verificar_monotonicidad_incidentes(const std::vector<int64_t>& ts) {
    auto [mono, inc] = verificar_monotonicidad(ts);
    std::vector<IncidenteHueco> out;
    for (const auto& m : inc) {
        IncidenteHueco h;
        h.inicio_utc = m.inicio_utc;
        h.fin_utc = m.fin_utc;
        h.duracion_segundos = m.duracion_segundos;
        h.tipo = TipoHueco::CORRUPCION_VALOR;
        h.descripcion = m.descripcion;
        h.activo = "GLOBAL_TIME_INDEX";
        out.push_back(h);
    }
    return {mono, out};
}

using ConteoHuecos = std::vector<std::pair<std::string, int64_t>>;   // {t.value: n} en el orden de TipoHueco
inline ConteoHuecos conteo_vacio() {
    ConteoHuecos c;
    for (int k = 0; k < 5; ++k) c.emplace_back(nombre(static_cast<TipoHueco>(k)), 0);
    return c;
}
inline void sumar_conteo(ConteoHuecos& a, TipoHueco t, int64_t n = 1) {
    for (auto& kv : a) if (kv.first == nombre(t)) { kv.second += n; return; }
}

inline std::pair<std::vector<IncidenteHueco>, ConteoHuecos> clasificar_huecos_ontologicos(const Tabla& df, const std::string& columna_fecha,
                                                                                         const std::string& columna_precio, Venue venue,
                                                                                         const std::string& activo,
                                                                                         double delta_esperado_segundos = 86400.0) {
    std::vector<IncidenteHueco> inc;
    ConteoHuecos conteo = conteo_vacio();
    if (df.vacia()) return {inc, conteo};
    const auto& fechas = df.col(columna_fecha).i;
    const std::vector<double> precios = a_dobles(df.col(columna_precio));
    for (size_t k = 0; k < precios.size(); ++k) {
        const double p = precios[k];
        if (std::isnan(p) || std::isinf(p) || p <= 0) {
            const std::string t = texto_timestamp(fechas[k]);
            IncidenteHueco h;
            h.inicio_utc = t;
            h.fin_utc = t;
            h.duracion_segundos = 0.0;
            h.tipo = TipoHueco::CORRUPCION_VALOR;
            h.descripcion = "Corrupción de valor: precio=" + repr_py(p) + " no físico o no finito en " + t;
            h.activo = activo;
            inc.push_back(h);
            sumar_conteo(conteo, TipoHueco::CORRUPCION_VALOR);
        }
    }
    if (fechas.size() > 1) {
        const double tolerancia = delta_esperado_segundos * 1.5;
        for (size_t k = 0; k + 1 < fechas.size(); ++k) {
            const int64_t salto = fechas[k + 1] / NS_S - fechas[k] / NS_S;
            if (!(static_cast<double>(salto) > tolerancia)) continue;
            const int64_t t_ini = fechas[k], t_fin = fechas[k + 1];
            const double sd = static_cast<double>(salto);
            TipoHueco tipo;
            std::string desc;
            if (venue == Venue::YAHOO || venue == Venue::BLOOMBERG_XLSX) {
                const int dia = dia_semana_lunes0(piso_div(t_ini, NS_DIA));
                if (delta_esperado_segundos >= 86400) {
                    if (dia == 4 && sd <= 3.5 * 86400) {
                        tipo = TipoHueco::ESTRUCTURAL_CALENDARIO;
                        desc = "Cierre de fin de semana TradFi (" + texto_fecha(t_ini) + " a " + texto_fecha(t_fin) + ")";
                    } else if (sd <= 4.5 * 86400 && (dia == 3 || dia == 4)) {
                        tipo = TipoHueco::ESTRUCTURAL_CALENDARIO;
                        desc = "Feriado / Fin de semana extendido TradFi (" + fmt_f(sd / 86400, 1) + " días)";
                    } else {
                        tipo = TipoHueco::CAIDA_RED_SOCKET;
                        desc = "HUECO ANÓMALO EN DÍA HÁBIL (" + fmt_f(sd / 86400, 1) + " días ausentes en feed " + nombre(venue) + ")";
                    }
                } else {
                    const int64_t hora = (piso_div(t_ini, NS_S) - piso_div(t_ini, NS_DIA) * 86400) / 3600;
                    if (dia == 4 && 40 * 3600 <= salto && salto <= 55 * 3600) {
                        tipo = TipoHueco::ESTRUCTURAL_CALENDARIO;
                        desc = "Cierre fin de semana intradía TradFi (" + fmt_f(sd / 3600, 1) + "h)";
                    } else if (3300 <= salto && salto <= 3900 && (hora == 21 || hora == 22)) {
                        tipo = TipoHueco::ESTRUCTURAL_CALENDARIO;
                        desc = "Pausa diaria de mantenimiento CME (" + fmt_f(sd / 60, 1) + " min)";
                    } else {
                        tipo = TipoHueco::CAIDA_RED_SOCKET;
                        desc = "BRECHA INTRADÍA EN SESIÓN (" + fmt_f(sd / 60, 1) + " min ausentes en feed " + nombre(venue) + ")";
                    }
                }
            } else {
                tipo = TipoHueco::CAIDA_RED_SOCKET;
                desc = "DESCONEXIÓN EN CRIPTO 24/7: Brecha de " + fmt_f(sd, 1) + " segundos sin datos";
            }
            IncidenteHueco h;
            h.inicio_utc = texto_timestamp(t_ini);
            h.fin_utc = texto_timestamp(t_fin);
            h.duracion_segundos = sd;
            h.tipo = tipo;
            h.descripcion = desc;
            h.activo = activo;
            inc.push_back(h);
            sumar_conteo(conteo, tipo);
        }
    }
    return {inc, conteo};
}

inline Tabla blindar_serie_con_mascara(const Tabla& df, const std::string& columna_fecha, const std::vector<std::string>& columnas_valor,
                                       bool imputar_ffill = true) {
    Tabla out = df;
    auto [mono, incm] = verificar_monotonicidad(out.col(columna_fecha).i);
    if (!mono) {
        const std::string detalle = incm.empty() ? std::string("Índice temporal no monótono") : incm.front().descripcion;
        throw ErrorMonotoniaTiempo("VIOLACIÓN DE FLECHA TEMPORAL EN BLINDAJE: " + detalle
                                   + ". Queda prohibido reordenar silenciosamente eventos fuera de secuencia.");
    }
    const size_t n = out.filas();
    std::vector<uint8_t> previo(n, 0);
    if (out.tiene("es_imputado")) {
        const auto v = a_dobles(out.col("es_imputado"));
        for (size_t k = 0; k < n; ++k) previo[k] = (!std::isnan(v[k]) && v[k] != 0.0) ? 1 : 0;
    }
    std::vector<std::string> presentes;
    for (const auto& c : columnas_valor) if (out.tiene(c)) presentes.push_back(c);
    if (presentes.empty())
        throw IncompatibilidadContratoError("Ninguna de las columnas de valor " + repr_lista_str(columnas_valor) + " existe en el DataFrame ("
                                            + repr_lista_str(nombres_columnas(out)) + ").");
    std::vector<uint8_t> invalido(n, 0);
    for (const auto& c : presentes) {
        std::vector<double> v = a_dobles(out.col(c));
        for (size_t k = 0; k < n; ++k)
            if (std::isnan(v[k]) || std::isinf(v[k]) || v[k] <= 0.0) { invalido[k] = 1; v[k] = NaN; }
        out.poner_f64(c, std::move(v));
    }
    if (n > 0 && invalido[0])
        throw ViolacionCausalidadError("CENSURA IZQUIERDA INICIAL: La primera muestra de la serie contiene valores no válidos. "
                                       "Está prohibido usar bfill (mirar hacia el futuro) para imputar el inicio de la serie.");
    if (imputar_ffill)
        for (const auto& c : presentes) {
            auto& v = out.col(c).f;
            for (size_t k = 1; k < v.size(); ++k) if (std::isnan(v[k])) v[k] = v[k - 1];
        }
    std::vector<int64_t> imp(n);
    for (size_t k = 0; k < n; ++k) imp[k] = (previo[k] || invalido[k]) ? 1 : 0;
    out.poner_i64("es_imputado", std::move(imp));
    return out;
}

inline const std::set<std::string>& columnas_no_perturbables() {
    static const std::set<std::string> s{"fecha_utc", "fecha_sesion", "es_imputado", "faltan", "market_open_mask"};
    return s;
}

// _diferencia_max_con_nan: NaN es un valor; patrones de NaN distintos -> inf.
inline double diferencia_max_con_nan(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < a.size(); ++k) if (std::isnan(a[k]) != std::isnan(b[k])) return std::numeric_limits<double>::infinity();
    double m = 0.0;
    bool hay = false;
    for (size_t k = 0; k < a.size(); ++k) {
        if (std::isnan(a[k]) || a[k] == b[k]) continue;
        const double d = std::fabs(a[k] - b[k]);
        if (!hay || d > m) { m = d; hay = true; }
    }
    return hay ? m : 0.0;
}

struct RepComparacion {
    int64_t filas_ref = 0, filas_obs = 0;
    std::vector<std::string> columnas_distintas;
    double max_diff = 0.0;
    std::string excepcion;   // solo testigo: "Clase: mensaje"
    JV a_jv(bool con_filas = true) const {
        JV j = JV::objeto();
        if (!excepcion.empty()) {
            j.poner("max_diff", JV::real(max_diff));
            j.poner("columnas_distintas", JV::lista_str(columnas_distintas));
            j.poner("excepcion", JV::cadena(excepcion));
            return j;
        }
        if (con_filas) { j.poner("filas_ref", JV::entero(filas_ref)); j.poner("filas_obs", JV::entero(filas_obs)); }
        j.poner("columnas_distintas", JV::lista_str(columnas_distintas));
        j.poner("max_diff", JV::real(max_diff));
        return j;
    }
};

inline RepComparacion comparar_tablas_hasta(const Tabla& ref, const Tabla& obs, int64_t t_corte, const std::string& col_fecha = "fecha_utc") {
    const Tabla a = filas_hasta(ref, t_corte, col_fecha), b = filas_hasta(obs, t_corte, col_fecha);
    RepComparacion r;
    r.filas_ref = static_cast<int64_t>(a.filas());
    r.filas_obs = static_cast<int64_t>(b.filas());
    if (a.filas() != b.filas() || a.col(col_fecha).i != b.col(col_fecha).i) {
        r.max_diff = std::numeric_limits<double>::infinity();
        r.columnas_distintas = {col_fecha};
        return r;
    }
    std::set<std::string> nom;
    for (const auto& c : a.columnas) nom.insert(c.nombre);
    for (const auto& c : b.columnas) nom.insert(c.nombre);
    for (const auto& n : nom) {
        double d;
        if (!a.tiene(n) || !b.tiene(n)) d = std::numeric_limits<double>::infinity();
        else {
            const Columna& ca = a.col(n);
            const Columna& cb = b.col(n);
            if (es_numerica(ca) && es_numerica(cb)) d = diferencia_max_con_nan(a_dobles(ca), a_dobles(cb));
            else {
                bool ig = true;
                for (size_t k = 0; k < a.filas() && ig; ++k) {
                    const bool na = celda_nula(ca, k), nb = celda_nula(cb, k);
                    if (na && nb) continue;
                    if (na != nb || celda_texto(ca, k) != celda_texto(cb, k)) ig = false;
                }
                d = ig ? 0.0 : std::numeric_limits<double>::infinity();
            }
        }
        if (d > 0.0) { r.columnas_distintas.push_back(n); r.max_diff = std::max(r.max_diff, d); }
    }
    return r;
}

using Crudos = std::vector<std::pair<std::string, Tabla>>;

inline Crudos perturbar_futuro(const Crudos& crudos, int64_t t_corte, numpy_rng::Generador& rng, const std::string& col_fecha = "fecha_utc",
                               double frac_nan = 0.03, double frac_baja = 0.03) {
    Crudos out;
    for (const auto& [sim, df] : crudos) {
        Tabla d = df;
        std::vector<size_t> idx;
        const auto& t = d.col(col_fecha).i;
        for (size_t k = 0; k < t.size(); ++k) if (t[k] != NAT && t[k] > t_corte) idx.push_back(k);
        if (!idx.empty()) {
            for (auto& c : d.columnas) {
                if (columnas_no_perturbables().count(c.nombre) || !es_numerica(c) || c.tipo == TipoColumna::BOOL) continue;
                std::vector<double> v = a_dobles(c);
                const auto ruido = rng.normal(0.0, 0.5, idx.size());
                for (size_t j = 0; j < idx.size(); ++j) v[idx[j]] = v[idx[j]] * std::exp(ruido[j]);
                const auto u = rng.random(idx.size());
                for (size_t j = 0; j < idx.size(); ++j) if (u[j] < frac_nan) v[idx[j]] = NaN;
                const std::string nom = c.nombre;
                c = Columna{};
                c.nombre = nom;
                c.tipo = TipoColumna::F64;
                c.f = std::move(v);
            }
            const auto u = rng.random(idx.size());
            std::vector<uint8_t> mant(d.filas(), 1);
            bool hay = false;
            for (size_t j = 0; j < idx.size(); ++j) if (u[j] < frac_baja) { mant[idx[j]] = 0; hay = true; }
            if (hay) d.filtrar(mant);
        }
        out.emplace_back(sim, std::move(d));
    }
    return out;
}

// Nombre de la clase Python de una excepcion C++ del porte (para "Clase: mensaje").
inline std::string clase_excepcion(const std::exception& e) {
    if (dynamic_cast<const SimboloNoEncontradoError*>(&e)) return "SimboloNoEncontradoError";
    if (dynamic_cast<const IncompatibilidadContratoError*>(&e)) return "IncompatibilidadContratoError";
    if (dynamic_cast<const ErrorMonotoniaTiempo*>(&e)) return "ErrorMonotoniaTiempo";
    if (dynamic_cast<const ViolacionCausalidadError*>(&e)) return "ViolacionCausalidadError";
    if (dynamic_cast<const ValueError*>(&e)) return "ValueError";
    return "Exception";
}

inline JV auditar_causalidad_pipeline(const std::function<Tabla(const Crudos&)>& procesar, const Crudos& crudos, int64_t n_cortes = 3,
                                      uint64_t semilla = 20261031, const std::string& col_fecha = "fecha_utc") {
    JV inf = JV::objeto();
    inf.poner("score", JV::real(NaN));
    inf.poner("n_cortes", JV::entero(0));
    inf.poner("cortes", JV::lista());
    inf.poner("fallos", JV::lista());
    const Tabla base = procesar(crudos);
    std::vector<int64_t> tiempos = base.col(col_fecha).i;
    std::sort(tiempos.begin(), tiempos.end());
    tiempos.erase(std::unique(tiempos.begin(), tiempos.end()), tiempos.end());
    const int64_t nt = static_cast<int64_t>(tiempos.size());
    if (nt < 10 || n_cortes < 1) {
        inf.poner("aviso", JV::cadena("Salida con " + std::to_string(nt) + " instantes: no hay cortes suficientes para certificar."));
        return inf;
    }
    numpy_rng::Generador rng(semilla);
    const int64_t lo = static_cast<int64_t>(0.2 * static_cast<double>(nt)), hi = static_cast<int64_t>(0.9 * static_cast<double>(nt));
    std::set<int64_t> pos;
    for (int64_t p : rng.integers(lo, std::max(lo + 1, hi), static_cast<size_t>(n_cortes))) pos.insert(p);
    std::vector<int64_t> invalidos;
    static const char* const PRECIOS[] = {"precio", "open", "high", "low", "close", "best_bid", "best_ask", "mid_price"};
    for (const auto& [sim, d] : crudos) {
        std::vector<std::vector<double>> cols;
        for (const char* c : PRECIOS) if (d.tiene(c)) cols.push_back(a_dobles(d.col(c)));
        if (cols.empty()) continue;
        const auto& t = d.col(col_fecha).i;
        for (size_t k = 0; k < d.filas(); ++k) {
            bool malo = false;
            for (const auto& v : cols) if (!std::isfinite(v[k]) || v[k] <= 0.0) { malo = true; break; }
            if (malo) invalidos.push_back(t[k]);
        }
    }
    if (base.tiene("es_imputado")) {
        const auto v = a_dobles(base.col("es_imputado"));
        const auto& t = base.col(col_fecha).i;
        for (size_t k = 0; k < v.size(); ++k) if (v[k] == 1.0) invalidos.push_back(t[k]);
    }
    int64_t n_dirigidos = 0;
    if (!invalidos.empty()) {
        const std::set<int64_t> conj(invalidos.begin(), invalidos.end());
        std::set<int64_t> ks;
        for (int64_t t : conj) ks.insert(static_cast<int64_t>(std::lower_bound(tiempos.begin(), tiempos.end(), t) - tiempos.begin()));
        std::vector<int64_t> kv;
        for (int64_t k : ks) if (lo <= k && k < hi && conj.count(tiempos[static_cast<size_t>(k)])) kv.push_back(k);
        const auto perm = rng.permutation(kv);
        for (size_t j = 0; j < perm.size() && j < 2; ++j) {
            n_dirigidos += pos.count(perm[j]) ? 0 : 1;
            pos.insert(perm[j]);
        }
    }
    inf.poner("cortes_dirigidos", JV::entero(n_dirigidos));
    double score = 0.0;
    JV& cortes = *inf.get("cortes");
    for (int64_t p : pos) {
        const int64_t c = tiempos[static_cast<size_t>(p)];
        cortes.push(JV::cadena(texto_timestamp(c)));
        Crudos trunc;
        for (const auto& [s, d] : crudos) trunc.emplace_back(s, filas_hasta(d, c, col_fecha));
        const Crudos pert = perturbar_futuro(crudos, c, rng, col_fecha);
        const std::pair<const char*, const Crudos*> variantes[2] = {{"truncado", &trunc}, {"perturbacion", &pert}};
        for (const auto& [nom, datos] : variantes) {
            RepComparacion rep;
            try {
                rep = comparar_tablas_hasta(base, procesar(*datos), c, col_fecha);
            } catch (const std::exception& e) {
                rep = RepComparacion{};
                rep.max_diff = std::numeric_limits<double>::infinity();
                rep.excepcion = clase_excepcion(e) + ": " + e.what();
            }
            if (rep.max_diff > 0.0) {
                JV f = JV::objeto();
                f.poner("corte", JV::cadena(texto_timestamp(c)));
                f.poner("prueba", JV::cadena(nom));
                for (auto& kv : rep.a_jv().o) f.poner(kv.first, kv.second);
                inf.get("fallos")->push(std::move(f));
            }
            score = std::max(score, rep.max_diff);
        }
    }
    inf.poner("score", JV::real(score));
    inf.poner("n_cortes", JV::entero(static_cast<int64_t>(pos.size())));
    return inf;
}

// ================================================================================================== catalogo_universo.py
struct ActivoPuntoEnElTiempo {
    std::string simbolo, fecha_listado;
    std::optional<std::string> fecha_deslistado;
    std::string motivo_desliste = "ACTIVO_VIVO";
    std::string venue, descripcion;
};

class CatalogoUniversoHistorico {
public:
    CatalogoUniversoHistorico() {
        const ActivoPuntoEnElTiempo base[] = {
            {"SPX", "1957-03-04", std::nullopt, "ACTIVO_VIVO", "YAHOO", "S&P 500 Index"},
            {"NDX", "1985-01-31", std::nullopt, "ACTIVO_VIVO", "YAHOO", "Nasdaq 100 Index"},
            {"VIX", "1993-01-02", std::nullopt, "ACTIVO_VIVO", "YAHOO", "Cboe Volatility Index"},
            {"DXY", "1973-03-01", std::nullopt, "ACTIVO_VIVO", "BLOOMBERG", "US Dollar Index"},
            {"ES1", "1997-09-09", std::nullopt, "ACTIVO_VIVO", "BLOOMBERG", "E-mini S&P 500 Future"},
            {"HYG", "2007-04-11", std::nullopt, "ACTIVO_VIVO", "YAHOO", "iShares iBoxx High Yield"},
            {"IEF", "2002-07-30", std::nullopt, "ACTIVO_VIVO", "YAHOO", "iShares 7-10 Year Treasury"},
            {"EEM", "2003-04-11", std::nullopt, "ACTIVO_VIVO", "YAHOO", "iShares MSCI Emerging Markets"},
            {"VNQ", "2004-09-29", std::nullopt, "ACTIVO_VIVO", "YAHOO", "Vanguard Real Estate ETF"},
            {"BTCUSDT", "2019-09-08", std::nullopt, "ACTIVO_VIVO", "BINANCE_UM", "Bitcoin USD-M Perp"},
            {"ETHUSDT", "2019-11-28", std::nullopt, "ACTIVO_VIVO", "BINANCE_UM", "Ethereum USD-M Perp"},
            {"SOLUSDT", "2020-09-24", std::nullopt, "ACTIVO_VIVO", "BINANCE_UM", "Solana USD-M Perp"},
            {"LUNAUSDT_OLD", "2020-08-27", std::string("2022-05-13"), "QUIEBRA_COLAPSO", "BINANCE_UM", "Terra Classic Colapso 2022"},
            {"FTTUSDT", "2021-01-15", std::string("2022-11-14"), "QUIEBRA_COLAPSO", "BINANCE_UM", "FTX Token Quiebra 2022"},
        };
        for (const auto& a : base) activos_.emplace_back(a.simbolo, a);
    }
    void registrar_activo(const ActivoPuntoEnElTiempo& a) {
        const std::string k = mayusculas(a.simbolo);
        for (auto& kv : activos_) if (kv.first == k) { kv.second = a; return; }
        activos_.emplace_back(k, a);
    }
    const ActivoPuntoEnElTiempo* obtener_activo(const std::string& s) const {
        const std::string k = mayusculas(s);
        for (const auto& kv : activos_) if (kv.first == k) return &kv.second;
        return nullptr;
    }

    // auditar_dataframe_supervivencia: anade/sustituye estado_supervivencia y anula las filas fuera del ciclo de vida.
    JV auditar_dataframe_supervivencia(Tabla& d, const std::string& simbolo, const std::string& columna_fecha = "fecha_sesion") const {
        const std::string su = mayusculas(simbolo);
        const ActivoPuntoEnElTiempo* f = obtener_activo(su);
        JV info = JV::objeto();
        info.poner("simbolo", JV::cadena(su));
        info.poner("en_catalogo", JV::booleano(f != nullptr));
        info.poner("filas_pre_listado", JV::entero(0));
        info.poner("filas_deslistado", JV::entero(0));
        info.poner("motivo_desliste", JV::cadena(f ? f->motivo_desliste : "NO_REGISTRADO"));
        info.poner("es_valido_causal", JV::booleano(true));
        const size_t n = d.filas();
        if (!f) {
            d.poner_str("estado_supervivencia", std::vector<std::string>(n, "COTIZANDO"));
            return info;
        }
        const std::string colf = d.tiene(columna_fecha) ? columna_fecha : "fecha_utc";
        const auto& t = d.col(colf).i;
        const int64_t zi = piso_div(ns_desde_iso(f->fecha_listado), NS_DIA);
        const std::optional<int64_t> zf = f->fecha_deslistado ? std::optional<int64_t>(piso_div(ns_desde_iso(*f->fecha_deslistado), NS_DIA))
                                                             : std::nullopt;
        std::vector<std::string> est(n);
        std::vector<uint8_t> pre(n, 0), post(n, 0);
        int64_t npre = 0, npost = 0;
        for (size_t k = 0; k < n; ++k) {
            const int64_t z = piso_div(t[k], NS_DIA);
            if (z < zi) { est[k] = "PRE_LISTADO"; pre[k] = 1; ++npre; }
            else if (zf && z > *zf) { est[k] = "DESLISTADO"; post[k] = 1; ++npost; }
            else est[k] = "COTIZANDO";
        }
        d.poner_str("estado_supervivencia", std::move(est));
        info.poner("filas_pre_listado", JV::entero(npre));
        info.poner("filas_deslistado", JV::entero(npost));
        info.poner("es_valido_causal", JV::booleano(npre == 0 && npost == 0));
        for (const auto* m : {&pre, &post}) {
            if ((m == &pre ? npre : npost) == 0) continue;
            if (d.tiene("es_sesion_con_dato")) { auto& v = d.col("es_sesion_con_dato").i; for (size_t k = 0; k < n; ++k) if ((*m)[k]) v[k] = 0; }
            if (d.tiene("hay_dato")) { auto& v = d.col("hay_dato").i; for (size_t k = 0; k < n; ++k) if ((*m)[k]) v[k] = 0; }
            if (d.tiene("es_imputado")) {
                Columna& c = d.col("es_imputado");
                if (c.tipo == TipoColumna::F64) { for (size_t k = 0; k < n; ++k) if ((*m)[k]) c.f[k] = 1.0; }
                else { for (size_t k = 0; k < n; ++k) if ((*m)[k]) c.i[k] = 1; }
            }
            for (const char* c : {"precio", "close", "open", "high", "low", "mid_price", "best_bid", "best_ask"}) {
                if (!d.tiene(c)) continue;
                Columna& col = d.col(c);
                if (col.tipo != TipoColumna::F64) { std::vector<double> v = a_dobles(col); d.poner_f64(c, std::move(v)); }
                auto& v = d.col(c).f;
                for (size_t k = 0; k < n; ++k) if ((*m)[k]) v[k] = NaN;
            }
        }
        return info;
    }

private:
    std::vector<std::pair<std::string, ActivoPuntoEnElTiempo>> activos_;
};

inline const CatalogoUniversoHistorico& catalogo_base() {
    static const CatalogoUniversoHistorico c;
    return c;
}

// ======================================================================================================== libro_hashes.py
inline std::string sha256_hex(const std::string& s) {
    uint8_t h[32];
    char x[65];
    sha256_bin(s.data(), s.size(), h);
    hex_de(h, x);
    return x;
}

inline std::string computar_hash_barra(const std::string& simbolo, const std::string& fecha_sesion, double precio,
                                       std::optional<double> o = std::nullopt, std::optional<double> h = std::nullopt,
                                       std::optional<double> l = std::nullopt, std::optional<double> v = std::nullopt) {
    const double pc = precio;
    const double po = (o && std::isfinite(*o)) ? *o : pc, ph = (h && std::isfinite(*h)) ? *h : pc, pl = (l && std::isfinite(*l)) ? *l : pc;
    const double vv = (v && std::isfinite(*v)) ? *v : 0.0;
    const std::string cad = mayusculas(recortar(simbolo)) + "|" + recortar(fecha_sesion) + "|" + fmt_f(po, 6) + "|" + fmt_f(ph, 6) + "|"
                            + fmt_f(pl, 6) + "|" + fmt_f(pc, 6) + "|" + fmt_f(vv, 2);
    return sha256_hex(cad);
}

inline std::string ledger_defecto() { return unir(raiz_datos(), "telemetria/ledger_ingesta_hashes.jsonl"); }

class LibroHashesInmutable {
public:
    explicit LibroHashesInmutable(std::string path = "") : path_(path.empty() ? ledger_defecto() : std::move(path)) { cargar(); }
    const std::string& path() const { return path_; }

    // auditar_y_registrar_barras sobre una tabla con fecha_sesion (o fecha_utc) y precio (o close) [y open/high/low/volume].
    JV auditar_y_registrar_barras(const Tabla& df, const std::string& simbolo, const std::string& venue = "DESCONOCIDO") {
        JV rep = JV::objeto();
        auto cerrar = [&](int64_t total, int64_t nuevas, int64_t verif, const JV& incs) {
            rep.poner("total_barras_evaluadas", JV::entero(total));
            rep.poner("nuevas_registradas", JV::entero(nuevas));
            rep.poner("verificadas_inmutables", JV::entero(verif));
            rep.poner("revisiones_detectadas", JV::entero(static_cast<int64_t>(incs.a.size())));
            rep.poner("es_inmutable", JV::booleano(incs.a.empty()));
            rep.poner("total_incidentes", JV::entero(static_cast<int64_t>(incs.a.size())));
            rep.poner("incidentes", incs);
            rep.poner("ledger_path", JV::cadena(path_));
        };
        if (df.vacia()) { cerrar(0, 0, 0, JV::lista()); return rep; }
        const std::string sk = mayusculas(simbolo);
        const std::string colf = df.tiene("fecha_sesion") ? "fecha_sesion" : "fecha_utc";
        const std::string colp = df.tiene("precio") ? "precio" : (df.tiene("close") ? "close" : "");
        if (colp.empty()) throw ValueError("No se encontró columna de precio ('precio' o 'close') para auditar hashes.");
        const auto& fechas = df.col(colf).i;
        const std::vector<double> pr = a_dobles(df.col(colp));
        auto opc = [&](const char* c) { return df.tiene(c) ? std::optional<std::vector<double>>(a_dobles(df.col(c))) : std::nullopt; };
        const auto vo = opc("open"), vh = opc("high"), vl = opc("low");
        const auto vv = df.tiene("volume") ? opc("volume") : opc("volumen");
        JV incs = JV::lista();
        std::vector<JV> nuevas;
        int64_t verif = 0;
        const std::string now_iso = ahora_iso();
        for (size_t k = 0; k < df.filas(); ++k) {
            const std::string fs = texto_fecha(fechas[k]);
            const double p = pr[k];
            if (!std::isfinite(p)) continue;
            const std::optional<double> o = vo ? std::optional<double>((*vo)[k]) : std::nullopt;
            const std::optional<double> h = vh ? std::optional<double>((*vh)[k]) : std::nullopt;
            const std::optional<double> l = vl ? std::optional<double>((*vl)[k]) : std::nullopt;
            const double vol = vv ? (*vv)[k] : 0.0;
            const std::string ha = computar_hash_barra(sk, fs, p, o, h, l, vol);
            const std::string clave = sk + "\x1f" + fs;
            auto it = indice_.find(clave);
            if (it != indice_.end()) {
                if (ha == it->second.first) { ++verif; continue; }
                const double pp = it->second.second;
                const double dif = std::fabs(p - pp) / std::max(1e-12, pp) * 100.0;
                JV inc = JV::objeto();
                inc.poner("simbolo", JV::cadena(sk));
                inc.poner("fecha_sesion", JV::cadena(fs));
                inc.poner("hash_previo", JV::cadena(it->second.first));
                inc.poner("hash_nuevo", JV::cadena(ha));
                inc.poner("precio_previo", JV::real(pp));
                inc.poner("precio_nuevo", JV::real(p));
                inc.poner("diferencia_relativa_pct", JV::real(round_py(dif, 4)));
                inc.poner("descripcion", JV::cadena("REVISIÓN RETROSPECTIVA DETECTADA en " + sk + " fecha " + fs + ": Precio previo congelado="
                                                    + fmt_f(pp, 4) + " vs Nuevo descargado=" + fmt_f(p, 4) + " (Variación: " + fmt_f(dif, 2)
                                                    + "%). Hash alterado."));
                incs.push(std::move(inc));
            } else {
                JV r = JV::objeto();
                r.poner("simbolo", JV::cadena(sk));
                r.poner("fecha_sesion", JV::cadena(fs));
                r.poner("fecha_utc", JV::cadena(texto_timestamp(fechas[k])));
                r.poner("venue", JV::cadena(venue));
                r.poner("hash_sha256", JV::cadena(ha));
                r.poner("precio_cierre", JV::real(p));
                r.poner("volumen", JV::real(vol));
                r.poner("timestamp_registro_utc", JV::cadena(now_iso));
                nuevas.push_back(std::move(r));
                indice_.emplace(clave, std::make_pair(ha, p));
            }
        }
        if (!nuevas.empty()) {
            std::filesystem::create_directories(std::filesystem::path(path_).parent_path());
            std::ofstream f(path_, std::ios::binary | std::ios::app);
            for (const auto& r : nuevas) f << json_dumps(r) << "\n";
        }
        cerrar(static_cast<int64_t>(df.filas()), static_cast<int64_t>(nuevas.size()), verif, incs);
        return rep;
    }

private:
    std::string path_;
    std::unordered_map<std::string, std::pair<std::string, double>> indice_;   // (simbolo, fecha) -> (hash, precio_cierre)

    void cargar() {
        if (!es_fichero(path_)) return;
        std::ifstream f(path_, std::ios::binary);
        std::string lin;
        while (std::getline(f, lin)) {
            const std::string t = recortar(lin);
            if (t.empty()) continue;
            try {
                const JSON r = parsear_json(t);
                const JSON* s = r.get("simbolo");
                const JSON* fs = r.get("fecha_sesion");
                const JSON* h = r.get("hash_sha256");
                if (!s || !fs || !h) continue;
                const std::string k = mayusculas(s->s) + "\x1f" + fs->s;
                if (indice_.count(k)) continue;
                const JSON* p = r.get("precio_cierre");
                indice_.emplace(k, std::make_pair(h->s, p ? p->n : 0.0));
            } catch (...) {
                continue;
            }
        }
    }
};

}  // namespace ark::dataset::ingesta
