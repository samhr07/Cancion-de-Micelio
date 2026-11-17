// tabla_ingesta.hpp -- BASE COMUN DEL PORTE DE LA INGESTA (ARK-DATA-02-CXX)
// Reactor ARK / Cancion del Micelio -- Jefe de Arquitectos de Codigo, 2026-11-02
//
// Calca de lo que los adaptadores Python de ARQUITECTO DE DATOS usan de pandas, contratos.py y control_calidad.py:
//   - Tabla: columnas con nombre y orden (como un DataFrame): F64, I64, TS (datetime64[ns, UTC] en ns), BOOL y STR.
//   - Excepciones con la misma jerarquia que contratos.py (SimboloNoEncontradoError es una IncompatibilidadContratoError;
//     ErrorMonotoniaTiempo es una ViolacionCausalidadError) y ValueError.
//   - pd.to_datetime(cadena ISO, utc=True) -> ns; pd.Timestamp -> texto como str() de pandas.
//   - pd.read_csv con conversion numerica (campo vacio o no numerico = NaN, como to_numeric(errors='coerce')).
//   - JSON minimo (manifiestos, metadatos de pandas en parquet y registros L3).
//   - control_calidad.verificar_monotonicidad, estampar_diario_disponible y verificar_estampa_diaria.
#pragma once
#include "ark/dataset/zona_horaria.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace ark::dataset::ingesta {

inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
inline constexpr int64_t NS_S = 1000000000LL;
inline constexpr int64_t NS_DIA = 86400LL * NS_S;
inline constexpr int64_t NAT = std::numeric_limits<int64_t>::min();   // pandas NaT

// ---------------------------------------------------------------------------------------------- excepciones (contratos.py)
struct IncompatibilidadContratoError : std::runtime_error { using std::runtime_error::runtime_error; };
struct SimboloNoEncontradoError : IncompatibilidadContratoError { using IncompatibilidadContratoError::IncompatibilidadContratoError; };
struct ViolacionCausalidadError : std::runtime_error { using std::runtime_error::runtime_error; };
struct ErrorMonotoniaTiempo : ViolacionCausalidadError { using ViolacionCausalidadError::ViolacionCausalidadError; };
struct ValueError : std::runtime_error { using std::runtime_error::runtime_error; };

// ------------------------------------------------------------------------------------------------------------------ tabla
enum class TipoColumna : uint8_t { F64, I64, TS, BOOL, STR };

struct Columna {
    std::string nombre;
    TipoColumna tipo = TipoColumna::F64;
    std::vector<double> f;          // F64
    std::vector<int64_t> i;         // I64, TS (ns UTC; NAT = nulo), BOOL (0/1)
    std::vector<std::string> s;     // STR
    // STR: nulo por celda (vacio = ninguna nula), como el NaN de la dtype 'str' de pandas 3. `objeto` = dtype object de pandas
    // (lo que deja un concat de una columna de texto, booleana o de fechas con NaN): sus celdas van como str(x) de Python.
    std::vector<uint8_t> sn;
    bool objeto = false;
    size_t size() const { return tipo == TipoColumna::F64 ? f.size() : (tipo == TipoColumna::STR ? s.size() : i.size()); }
    bool es_nulo_str(size_t k) const { return !sn.empty() && sn[k] != 0; }
};

struct Tabla {
    std::vector<Columna> columnas;

    size_t filas() const { return columnas.empty() ? 0 : columnas.front().size(); }
    bool vacia() const { return filas() == 0; }
    bool tiene(const std::string& n) const { for (const auto& c : columnas) if (c.nombre == n) return true; return false; }
    Columna& col(const std::string& n) {
        for (auto& c : columnas) if (c.nombre == n) return c;
        throw ValueError("columna inexistente: " + n);
    }
    const Columna& col(const std::string& n) const {
        for (const auto& c : columnas) if (c.nombre == n) return c;
        throw ValueError("columna inexistente: " + n);
    }
    Columna& poner(const std::string& n, TipoColumna t) {   // df[n] = ... : sustituye en su sitio o anade al final
        for (auto& c : columnas) if (c.nombre == n) { c = Columna{}; c.nombre = n; c.tipo = t; return c; }
        columnas.push_back(Columna{});
        columnas.back().nombre = n;
        columnas.back().tipo = t;
        return columnas.back();
    }
    void poner_f64(const std::string& n, std::vector<double> v) { poner(n, TipoColumna::F64).f = std::move(v); }
    void poner_i64(const std::string& n, std::vector<int64_t> v) { poner(n, TipoColumna::I64).i = std::move(v); }
    void poner_ts(const std::string& n, std::vector<int64_t> v) { poner(n, TipoColumna::TS).i = std::move(v); }
    void poner_bool(const std::string& n, std::vector<int64_t> v) { poner(n, TipoColumna::BOOL).i = std::move(v); }
    void poner_str(const std::string& n, std::vector<std::string> v) { poner(n, TipoColumna::STR).s = std::move(v); }
    void quitar(const std::string& n) {
        columnas.erase(std::remove_if(columnas.begin(), columnas.end(), [&](const Columna& c) { return c.nombre == n; }), columnas.end());
    }
    void renombrar(const std::string& a, const std::string& b) { for (auto& c : columnas) if (c.nombre == a) c.nombre = b; }
    // df[mask] (y reset_index): conserva las filas con m[k] != 0.
    void filtrar(const std::vector<uint8_t>& m) {
        for (auto& c : columnas) {
            size_t w = 0;
            if (c.tipo == TipoColumna::F64) { for (size_t k = 0; k < c.f.size(); ++k) if (m[k]) c.f[w++] = c.f[k]; c.f.resize(w); }
            else if (c.tipo == TipoColumna::STR) {
                const bool hay_n = !c.sn.empty();
                for (size_t k = 0; k < c.s.size(); ++k) if (m[k]) { if (hay_n) c.sn[w] = c.sn[k]; c.s[w++] = c.s[k]; }
                c.s.resize(w);
                if (hay_n) c.sn.resize(w);
            }
            else { for (size_t k = 0; k < c.i.size(); ++k) if (m[k]) c.i[w++] = c.i[k]; c.i.resize(w); }
        }
    }
    // df[cols] en ese orden (las que existan).
    Tabla seleccionar(const std::vector<std::string>& cols) const {
        Tabla t;
        for (const auto& n : cols) if (tiene(n)) t.columnas.push_back(col(n));
        return t;
    }
};

// -------------------------------------------------------------------------------------------------------- tiempo (pandas)
inline int64_t ns_desde_civil(int y, int m, int d, int hh = 0, int mi = 0, int ss = 0, int64_t frac_ns = 0) {
    return (dias_desde_civil(y, m, d) * 86400LL + hh * 3600LL + mi * 60LL + ss) * NS_S + frac_ns;
}

// pd.to_datetime(cadena, utc=True) para las formas que usa el proyecto: 'YYYY-MM-DD', 'YYYY-MM-DD HH:MM[:SS[.f...]]', con 'T'
// o espacio, y sufijo opcional 'Z', '+HH:MM', '-HH:MM' o '+HHMM'. Sin sufijo es UTC.
inline int64_t ns_desde_iso(const std::string& t) {
    auto num = [&](size_t a, size_t b) -> int {
        if (b > t.size()) throw ValueError("fecha ISO invalida: " + t);
        int v = 0;
        for (size_t k = a; k < b; ++k) { if (t[k] < '0' || t[k] > '9') throw ValueError("fecha ISO invalida: " + t); v = v * 10 + (t[k] - '0'); }
        return v;
    };
    if (t.size() < 10 || t[4] != '-' || t[7] != '-') throw ValueError("fecha ISO invalida: " + t);
    const int y = num(0, 4), m = num(5, 7), d = num(8, 10);
    if (!fecha_valida(static_cast<uint32_t>(y * 10000 + m * 100 + d))) throw ValueError("fecha inexistente: " + t);
    int hh = 0, mi = 0, ss = 0;
    int64_t frac = 0, off = 0;
    size_t p = 10;
    if (p < t.size() && (t[p] == 'T' || t[p] == ' ')) {
        hh = num(p + 1, p + 3);
        mi = num(p + 4, p + 6);
        p += 6;
        if (p < t.size() && t[p] == ':') { ss = num(p + 1, p + 3); p += 3; }
        if (p < t.size() && t[p] == '.') {
            ++p;
            int64_t esc = 100000000LL;
            while (p < t.size() && t[p] >= '0' && t[p] <= '9') { frac += (t[p] - '0') * esc; esc /= 10; ++p; }
        }
    }
    if (p < t.size()) {
        if (t[p] == 'Z') { ++p; }
        else if (t[p] == '+' || t[p] == '-') {
            const int sg = t[p] == '-' ? -1 : 1;
            const int oh = num(p + 1, p + 3);
            size_t q = p + 3;
            if (q < t.size() && t[q] == ':') ++q;
            const int om = num(q, q + 2);
            off = sg * (oh * 3600LL + om * 60LL);
            p = q + 2;
        }
    }
    while (p < t.size() && t[p] == ' ') ++p;
    if (p != t.size()) throw ValueError("fecha ISO invalida: " + t);
    return ns_desde_civil(y, m, d, hh, mi, ss, frac) - off * NS_S;
}

inline int64_t piso_div(int64_t a, int64_t b) { return a / b - ((a % b != 0) && ((a < 0) != (b < 0))); }
inline int64_t normalizar_dia(int64_t ns) { return piso_div(ns, NS_DIA) * NS_DIA; }   // .dt.normalize() en UTC

// str(pd.Timestamp(ns, tz='UTC')): 'YYYY-MM-DD HH:MM:SS[.ffffff | .fffffffff]+00:00'.
inline std::string texto_timestamp(int64_t ns) {
    if (ns == NAT) return "NaT";
    const int64_t s = piso_div(ns, NS_S), frac = ns - s * NS_S;
    const int64_t dia = piso_div(s, 86400), seg = s - dia * 86400;
    const FechaCivil c = civil_desde_dias(dia);
    char b[64];
    std::snprintf(b, sizeof b, "%04d-%02d-%02d %02lld:%02lld:%02lld", c.y, c.m, c.d, static_cast<long long>(seg / 3600),
                  static_cast<long long>((seg / 60) % 60), static_cast<long long>(seg % 60));
    std::string r = b;
    if (frac) {
        char f[16];
        if (frac % 1000 == 0) std::snprintf(f, sizeof f, ".%06lld", static_cast<long long>(frac / 1000));
        else std::snprintf(f, sizeof f, ".%09lld", static_cast<long long>(frac));
        r += f;
    }
    return r + "+00:00";
}

// UTC -> desplazamiento de America/New_York en ese instante (tz_convert). Verano: de 07:00 UTC (02:00 EST) del dia de inicio
// a 06:00 UTC (02:00 EDT) del dia de fin.
inline int64_t offset_ny_desde_utc_s(int64_t utc_s) {
    const int anio = civil_desde_dias(piso_div(utc_s - 5 * 3600, 86400)).y;
    int64_t ini = 0, fin = 0;
    if (!detalle_zona::verano_ny(anio, ini, fin)) return -5 * 3600;
    const int64_t a = ini * 86400 + 7 * 3600, b = fin * 86400 + 6 * 3600;
    return (utc_s >= a && utc_s < b) ? -4 * 3600 : -5 * 3600;
}

// ------------------------------------------------------------------------------------------------------- rutas y ficheros
inline std::string variable_entorno(const char* n) {
#ifdef _MSC_VER
    char* v = nullptr;
    size_t len = 0;
    std::string r;
    if (_dupenv_s(&v, &len, n) == 0 && v) { r = v; std::free(v); }
    return r;
#else
    const char* v = std::getenv(n);
    return v ? std::string(v) : std::string();
#endif
}

inline std::string raiz_datos() {
    if (const std::string e = variable_entorno("ARK_RAIZ_DATOS"); !e.empty()) return e;
    if (const std::string e = variable_entorno("ARK_RAIZ_REPO"); !e.empty()) return e;
    throw IncompatibilidadContratoError("ARK_RAIZ_DATOS no definido: el C++ no infiere la raiz del repositorio");
}
// os.path.join(raiz, "cachas", "x") de Python: separador nativo en toda la ruta (los mensajes de error la citan).
inline std::string unir(const std::string& a, const std::string& b) {
    return (std::filesystem::path(a) / std::filesystem::path(b)).make_preferred().string();
}
inline bool es_fichero(const std::string& r) { std::error_code ec; return std::filesystem::is_regular_file(r, ec); }

// sorted(glob(dir/prefijo*sufijo)) en un solo directorio.
inline std::vector<std::string> glob_simple(const std::string& dir, const std::string& prefijo, const std::string& sufijo) {
    std::vector<std::string> v;
    std::error_code ec;
    if (!std::filesystem::is_directory(dir, ec)) return v;
    for (const auto& e : std::filesystem::directory_iterator(dir, ec)) {
        const std::string n = e.path().filename().string();
        if (n.size() >= prefijo.size() + sufijo.size() && n.compare(0, prefijo.size(), prefijo) == 0
            && n.compare(n.size() - sufijo.size(), sufijo.size(), sufijo) == 0)
            v.push_back(unir(dir, n));
    }
    std::sort(v.begin(), v.end());
    return v;
}

inline std::string mayusculas(std::string s) { for (char& c : s) if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 32); return s; }
inline std::string minusculas(std::string s) { for (char& c : s) if (c >= 'A' && c <= 'Z') c = static_cast<char>(c + 32); return s; }
inline std::string reemplazar(std::string s, const std::string& a, const std::string& b) {
    if (a.empty()) return s;
    size_t p = 0;
    while ((p = s.find(a, p)) != std::string::npos) { s.replace(p, a.size(), b); p += b.size(); }
    return s;
}
inline std::string recortar(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

// ------------------------------------------------------------------------------------------- numeros de pandas (C engine)
// Calca de precise_xstrtod (pandas/_libs/src/parser/tokenizer.c), el conversor por omision de read_csv (float_precision=None)
// y de to_numeric (to_double): acumula hasta 17 digitos significativos en un double y escala con una tabla de potencias de 10.
// NO es un redondeo correcto: con 17 digitos puede quedar a 1 ulp del valor exacto, y eso es lo que el dataset Python ve.
inline const double* potencias_10_pandas() {
    static const std::vector<double> e = [] {
        std::vector<double> v(309);
        for (int k = 0; k <= 308; ++k) v[static_cast<size_t>(k)] = std::strtod(("1e" + std::to_string(k)).c_str(), nullptr);   // literal 1eK
        return v;
    }();
    return e.data();
}

// Devuelve true y el valor si la cadena entera es un numero para pandas (con espacios a los lados permitidos).
inline bool pandas_xstrtod(const std::string& s, double& out) {
    const double* e = potencias_10_pandas();
    const char* p = s.c_str();
    auto esp = [](char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f'; };
    auto dig = [](char c) { return c >= '0' && c <= '9'; };
    while (esp(*p)) ++p;
    bool neg = false;
    if (*p == '-') { neg = true; ++p; } else if (*p == '+') { ++p; }
    double number = 0.0;
    int exponent = 0, num_digits = 0, num_decimals = 0;
    const int max_digits = 17;
    while (dig(*p)) {
        if (num_digits < max_digits) { number = number * 10. + (*p - '0'); ++num_digits; }
        else ++exponent;
        ++p;
    }
    if (*p == '.') {
        ++p;
        while (num_digits < max_digits && dig(*p)) { number = number * 10. + (*p - '0'); ++p; ++num_digits; ++num_decimals; }
        if (num_digits >= max_digits) while (dig(*p)) ++p;
        exponent -= num_decimals;
    }
    if (num_digits == 0) return false;
    if (neg) number = -number;
    if (*p == 'e' || *p == 'E') {
        bool nexp = false;
        ++p;
        if (*p == '-') { nexp = true; ++p; } else if (*p == '+') { ++p; }
        int nd = 0, n = 0;
        while (nd < max_digits && dig(*p)) { n = n * 10 + (*p - '0'); ++nd; ++p; }
        exponent += nexp ? -n : n;
        if (nd == 0) --p;
    }
    if (exponent > 308) return false;   // ERANGE: el conversor lo da por no numerico
    if (exponent > 0) number *= e[exponent];
    else if (exponent < -308) {
        if (exponent < -616) number = 0.;
        else { number /= e[-308 - exponent]; number /= e[308]; }
    } else number /= e[-exponent];
    if (number == HUGE_VAL || number == -HUGE_VAL) return false;
    while (esp(*p)) ++p;
    if (*p != '\0') return false;
    out = number;
    return true;
}

// Suma compensada de pandas (group_sum / group_mean de _libs/groupby.pyx): Kahan con la compensacion a cero si sale NaN.
struct SumaKahanPandas {
    double sumx = 0.0, comp = 0.0;
    int64_t nobs = 0;
    void anadir(double v) {
        if (std::isnan(v)) return;
        ++nobs;
        const double y = v - comp, t = sumx + y;
        comp = t - sumx - y;
        if (comp != comp) comp = 0.0;
        sumx = t;
    }
};

// ---------------------------------------------------------------------------------------------------------- CSV (pandas)
struct CSV {
    std::vector<std::string> cabecera;
    std::vector<std::vector<std::string>> filas;
    int indice(const std::string& n) const {
        for (size_t k = 0; k < cabecera.size(); ++k) if (cabecera[k] == n) return static_cast<int>(k);
        return -1;
    }
    std::vector<std::string> texto(const std::string& n) const {
        const int k = indice(n);
        if (k < 0) throw ValueError("columna CSV inexistente: " + n);
        std::vector<std::string> v;
        v.reserve(filas.size());
        for (const auto& f : filas) v.push_back(static_cast<size_t>(k) < f.size() ? f[static_cast<size_t>(k)] : std::string());
        return v;
    }
    // pd.to_numeric(errors='coerce') de lo que read_csv dejo en la columna: numero completo o NaN.
    std::vector<double> numerica(const std::string& n) const {
        std::vector<double> v;
        for (const auto& s : texto(n)) {
            double x = NaN;
            v.push_back(pandas_xstrtod(s, x) ? x : NaN);
        }
        return v;
    }
};

inline CSV leer_csv(const std::string& ruta) {
    std::ifstream f(ruta, std::ios::binary);
    if (!f) throw SimboloNoEncontradoError("no se puede abrir " + ruta);
    CSV c;
    std::string lin;
    bool primera = true;
    while (std::getline(f, lin)) {
        if (!lin.empty() && lin.back() == '\r') lin.pop_back();
        if (primera && lin.size() >= 3 && static_cast<unsigned char>(lin[0]) == 0xEF) lin = lin.substr(3);   // BOM
        if (lin.empty()) continue;   // read_csv salta lineas vacias
        std::vector<std::string> campos;
        std::string cur;
        bool comillas = false;
        for (char ch : lin) {
            if (ch == '"') { comillas = !comillas; continue; }
            if (ch == ',' && !comillas) { campos.push_back(cur); cur.clear(); continue; }
            cur += ch;
        }
        campos.push_back(cur);
        if (primera) { c.cabecera = campos; primera = false; }
        else c.filas.push_back(std::move(campos));
    }
    return c;
}

// ---------------------------------------------------------------------------------------------------------- JSON minimo
struct JSON {
    enum class T : uint8_t { NULO, BOOL, NUM, STR, ARR, OBJ } t = T::NULO;
    bool b = false;
    double n = 0.0;
    bool entero = false;
    int64_t i = 0;
    std::string s;
    std::vector<JSON> a;
    std::vector<std::pair<std::string, JSON>> o;
    const JSON* get(const std::string& k) const { for (const auto& kv : o) if (kv.first == k) return &kv.second; return nullptr; }
};

namespace detalle_json {
struct P {
    const char* p;
    const char* e;
    void esp() { while (p < e && (*p == ' ' || *p == '\n' || *p == '\r' || *p == '\t')) ++p; }
    [[noreturn]] void error() { throw ValueError("JSON invalido"); }
    std::string cadena() {
        if (p >= e || *p != '"') error();
        ++p;
        std::string r;
        while (p < e && *p != '"') {
            if (*p == '\\') {
                ++p;
                if (p >= e) error();
                switch (*p) {
                    case 'n': r += '\n'; break; case 't': r += '\t'; break; case 'r': r += '\r'; break;
                    case 'b': r += '\b'; break; case 'f': r += '\f'; break;
                    case 'u': {
                        if (e - p < 5) error();
                        unsigned cp = static_cast<unsigned>(std::stoul(std::string(p + 1, p + 5), nullptr, 16));
                        p += 4;
                        if (cp >= 0xD800 && cp < 0xDC00 && e - p > 6 && p[1] == '\\' && p[2] == 'u') {
                            const unsigned lo = static_cast<unsigned>(std::stoul(std::string(p + 3, p + 7), nullptr, 16));
                            cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                            p += 6;
                        }
                        if (cp < 0x80) r += static_cast<char>(cp);
                        else if (cp < 0x800) { r += static_cast<char>(0xC0 | (cp >> 6)); r += static_cast<char>(0x80 | (cp & 0x3F)); }
                        else if (cp < 0x10000) { r += static_cast<char>(0xE0 | (cp >> 12)); r += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); r += static_cast<char>(0x80 | (cp & 0x3F)); }
                        else { r += static_cast<char>(0xF0 | (cp >> 18)); r += static_cast<char>(0x80 | ((cp >> 12) & 0x3F)); r += static_cast<char>(0x80 | ((cp >> 6) & 0x3F)); r += static_cast<char>(0x80 | (cp & 0x3F)); }
                        break;
                    }
                    default: r += *p;
                }
                ++p;
            } else {
                r += *p++;
            }
        }
        if (p >= e) error();
        ++p;
        return r;
    }
    JSON valor() {
        esp();
        if (p >= e) error();
        JSON v;
        if (*p == '{') {
            v.t = JSON::T::OBJ; ++p; esp();
            if (p < e && *p == '}') { ++p; return v; }
            for (;;) {
                esp(); std::string k = cadena(); esp();
                if (p >= e || *p != ':') error();
                ++p;
                v.o.emplace_back(std::move(k), valor());
                esp();
                if (p < e && *p == ',') { ++p; continue; }
                if (p < e && *p == '}') { ++p; return v; }
                error();
            }
        }
        if (*p == '[') {
            v.t = JSON::T::ARR; ++p; esp();
            if (p < e && *p == ']') { ++p; return v; }
            for (;;) {
                v.a.push_back(valor()); esp();
                if (p < e && *p == ',') { ++p; continue; }
                if (p < e && *p == ']') { ++p; return v; }
                error();
            }
        }
        if (*p == '"') { v.t = JSON::T::STR; v.s = cadena(); return v; }
        if (e - p >= 4 && std::string(p, p + 4) == "true") { v.t = JSON::T::BOOL; v.b = true; p += 4; return v; }
        if (e - p >= 5 && std::string(p, p + 5) == "false") { v.t = JSON::T::BOOL; p += 5; return v; }
        if (e - p >= 4 && std::string(p, p + 4) == "null") { p += 4; return v; }
        if (e - p >= 3 && std::string(p, p + 3) == "NaN") { v.t = JSON::T::NUM; v.n = NaN; p += 3; return v; }
        const char* q = p;
        bool ent = true;
        if (q < e && (*q == '-' || *q == '+')) ++q;
        while (q < e && ((*q >= '0' && *q <= '9') || *q == '.' || *q == 'e' || *q == 'E' || *q == '-' || *q == '+')) {
            if (*q == '.' || *q == 'e' || *q == 'E') ent = false;
            ++q;
        }
        if (q == p) error();
        const std::string num(p, q);
        v.t = JSON::T::NUM;
        v.n = std::strtod(num.c_str(), nullptr);
        v.entero = ent;
        if (ent) v.i = std::strtoll(num.c_str(), nullptr, 10);
        p = q;
        return v;
    }
};
}  // namespace detalle_json

inline JSON parsear_json(const std::string& t) {
    detalle_json::P q{t.data(), t.data() + t.size()};
    JSON v = q.valor();
    q.esp();
    if (q.p != q.e) throw ValueError("JSON con texto sobrante");
    return v;
}

// ------------------------------------------------------------------------------------------- control_calidad.py (guardas)
struct IncidenteMonotonia {
    std::string inicio_utc, fin_utc, descripcion;
    double duracion_segundos = 0.0;
};

// verificar_monotonicidad: compara en SEGUNDOS (astype('datetime64[s]')): dos marcas del mismo segundo son un duplicado.
inline std::pair<bool, std::vector<IncidenteMonotonia>> verificar_monotonicidad(const std::vector<int64_t>& ts_ns) {
    std::vector<IncidenteMonotonia> inc;
    if (ts_ns.size() <= 1) return {true, inc};
    auto seg = [](int64_t ns) {   // datetime64[ns] -> [s] trunca hacia cero (como numpy)
        return ns / NS_S;
    };
    for (size_t k = 0; k + 1 < ts_ns.size(); ++k) {
        const int64_t d = seg(ts_ns[k + 1]) - seg(ts_ns[k]);
        if (d <= 0) {
            IncidenteMonotonia x;
            x.inicio_utc = texto_timestamp(ts_ns[k]);
            x.fin_utc = texto_timestamp(ts_ns[k + 1]);
            x.duracion_segundos = static_cast<double>(d);
            x.descripcion = d == 0 ? ("Timestamp duplicado detectado en índice " + std::to_string(k) + " (" + x.inicio_utc + ")")
                                   : ("INVERSIÓN TEMPORAL (Viaje en el tiempo): " + x.fin_utc + " ocurre antes o igual a " + x.inicio_utc);
            inc.push_back(x);
        }
    }
    return {inc.empty(), inc};
}

inline constexpr int64_t DESFASE_DISPONIBILIDAD_DIARIA_NS = NS_DIA;

// estampar_diario_disponible: (fecha_sesion 00:00 UTC, fecha_utc = sesion + 1 dia); una sesion con hora es un error.
inline std::pair<std::vector<int64_t>, std::vector<int64_t>> estampar_diario_disponible(const std::vector<std::string>& fechas) {
    std::vector<int64_t> ses, disp;
    for (const auto& f : fechas) {
        const int64_t t = ns_desde_iso(recortar(f));
        if (t != normalizar_dia(t))
            throw ViolacionCausalidadError("ESTAMPADO DIARIO: la fecha de sesión trae hora; se esperaba un día puro (YYYY-MM-DD). "
                                           "Una hora dentro de la sesión no dice cuándo cerró la barra.");
        ses.push_back(t);
        disp.push_back(t + DESFASE_DISPONIBILIDAD_DIARIA_NS);
    }
    return {ses, disp};
}

inline void verificar_estampa_diaria(const Tabla& df, const std::string& col_fecha = "fecha_utc", const std::string& col_sesion = "fecha_sesion") {
    if (df.vacia() || !df.tiene(col_sesion)) return;
    const auto& t = df.col(col_fecha).i;
    const auto& s = df.col(col_sesion).i;
    size_t k = 0;
    int64_t mejor = std::numeric_limits<int64_t>::min();
    for (size_t r = 0; r < t.size(); ++r) {
        const int64_t adel = (s[r] + DESFASE_DISPONIBILIDAD_DIARIA_NS) - t[r];
        if (adel > mejor) { mejor = adel; k = r; }
    }
    if (mejor > 0) {
        // str(pd.Timedelta): 'D days HH:MM:SS[.ffffff|.fffffffff]'
        const int64_t dd = mejor / NS_DIA, r = mejor % NS_DIA, sg = r / NS_S, fr = r % NS_S;
        char td[96];
        std::snprintf(td, sizeof td, "%lld days %02lld:%02lld:%02lld", static_cast<long long>(dd), static_cast<long long>(sg / 3600),
                      static_cast<long long>((sg / 60) % 60), static_cast<long long>(sg % 60));
        std::string tds = td;
        if (fr) {
            char x[16];
            if (fr % 1000 == 0) std::snprintf(x, sizeof x, ".%06lld", static_cast<long long>(fr / 1000));
            else std::snprintf(x, sizeof x, ".%09lld", static_cast<long long>(fr));
            tds += x;
        }
        throw ViolacionCausalidadError("ESTAMPA DIARIA ANTICIPADA: la sesión " + texto_timestamp(s[k]).substr(0, 10) + " está rotulada en "
                                       + texto_timestamp(t[k]) + ", " + tds + " antes de existir.");
    }
}

// Filtro por rango [ini, fin] sobre una columna TS (df[col >= t_ini]; df[col <= t_fin]).
inline void filtrar_rango(Tabla& df, const std::string& col, const std::string& ini, const std::string& fin) {
    if (!ini.empty()) {
        const int64_t a = ns_desde_iso(ini);
        std::vector<uint8_t> m(df.filas());
        const auto& t = df.col(col).i;
        for (size_t k = 0; k < m.size(); ++k) m[k] = (t[k] != NAT && t[k] >= a) ? 1 : 0;
        df.filtrar(m);
    }
    if (!fin.empty()) {
        const int64_t b = ns_desde_iso(fin);
        std::vector<uint8_t> m(df.filas());
        const auto& t = df.col(col).i;
        for (size_t k = 0; k < m.size(); ++k) m[k] = (t[k] != NAT && t[k] <= b) ? 1 : 0;
        df.filtrar(m);
    }
}

}  // namespace ark::dataset::ingesta
