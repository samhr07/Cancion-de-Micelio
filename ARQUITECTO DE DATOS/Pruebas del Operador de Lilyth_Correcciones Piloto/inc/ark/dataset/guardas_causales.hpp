// guardas_causales.hpp -- VENTANA DE INFORMACION DE CADA BARRA DIARIA Y GUARDA DE PARES PREDICTIVOS (ARK-DATA-01-CXX)
// Reactor ARK / Cancion del Micelio -- correccion del Jefe 2026-11-01
//
// Porte EXACTO de ARQUITECTO DE DATOS/recopilacion_unificada/ventanas_informacion.py: perfiles por instrumento (cierre e
// inicio de la barra como INTERVALOS [lo, hi] en UTC, con su procedencia REGLA / MEDIDO / NO_VERIFICADO), ventana_barra() y
// la guarda con latencia L >= 0:
//     EJECUTABLE     si  c_X^hi + L <  s_Y^lo
//     REZAGO         si  c_X^lo + L >= s_Y^hi
//     NO_VERIFICADO  en otro caso (los intervalos se pisan)
// La desigualdad es ESTRICTA: observar un precio y operar a ese mismo precio no es ejecutable.
//
// Lo que corrige frente a la version anterior: la hora de cierre ya no es "NYSE 16:00 ET para todo TradFi" (DXY, VIX, la
// curva del Tesoro y los futuros son intervalos NO_VERIFICADO o MEDIDO; "SPX -> DXY +0 ms" era precision inventada), y la
// frontera es estricta (antes margen >= L daba EJECUTABLE con margen 0 y L = 0).
#pragma once

#include "tipos_dataset.hpp"
#include "calendario.hpp"
#include <cstdint>
#include <cmath>
#include <limits>
#include <string>

namespace ark::dataset {

enum class ClasificacionCausalPredictiva : uint8_t {
    EJECUTABLE = 0,
    REZAGO = 1,
    NO_VERIFICADO = 2,
    SIN_SESION = 3,    // el objetivo es un dia sin sesion: no tiene retorno
    NO_DEFINIDO = 4    // apertura a cierre sobre un objetivo sin apertura (un fixing)
};

[[nodiscard]] constexpr const char* nombre_clase(ClasificacionCausalPredictiva c) noexcept {
    switch (c) {
        case ClasificacionCausalPredictiva::EJECUTABLE: return "EJECUTABLE";
        case ClasificacionCausalPredictiva::REZAGO: return "REZAGO";
        case ClasificacionCausalPredictiva::NO_VERIFICADO: return "NO_VERIFICADO";
        case ClasificacionCausalPredictiva::SIN_SESION: return "SIN_SESION";
        default: return "NO_DEFINIDO";
    }
}

enum class Procedencia : uint8_t { REGLA = 0, MEDIDO = 1, NO_VERIFICADO = 2 };
[[nodiscard]] constexpr const char* nombre_procedencia(Procedencia p) noexcept {
    return p == Procedencia::REGLA ? "REGLA" : (p == Procedencia::MEDIDO ? "MEDIDO" : "NO_VERIFICADO");
}

enum class CorteDiarioCripto : uint8_t { UTC_00 = 0, NY_16 = 1 };

constexpr double LATENCIA_DIARIA_S = 0.0;
constexpr double LATENCIA_INTRADIA_S = 0.208;   // la medida por el ARK (ventanas_informacion.LATENCIA_INTRADIA_S)

struct Instante {
    Zona zona;
    int dias;
    double hora;
};

struct PerfilVentana {
    const char* nombre;
    Instante cierre_lo, cierre_hi;
    Procedencia proc_cierre;
    bool tiene_cierre_reducido;
    Instante red_lo, red_hi;
    Procedencia proc_cierre_reducido;
    bool tiene_inicio;
    Instante ini_lo, ini_hi;
    Procedencia proc_inicio;
    bool inicio_es_cierre_previo;
    bool negociable;
};

namespace detalle_ventana {
constexpr Zona NY = Zona::NUEVA_YORK, U = Zona::UTC, BER = Zona::BERLIN;
constexpr Procedencia RG = Procedencia::REGLA, ME = Procedencia::MEDIDO, NV = Procedencia::NO_VERIFICADO;
constexpr Instante X{U, 0, 0.0};   // relleno de los campos que no aplican

inline const PerfilVentana PERFILES[] = {
    {"NYSE_RV", {NY, 0, 16.0}, {NY, 0, 16.0}, RG, true, {NY, 0, 13.0}, {NY, 0, 13.0}, RG, true, {NY, 0, 9.5}, {NY, 0, 9.5}, RG, false, true},
    {"NYSE_INDICE", {NY, 0, 16.0}, {NY, 0, 16.0}, RG, true, {NY, 0, 13.0}, {NY, 0, 13.0}, RG, true, {NY, 0, 9.5}, {NY, 0, 9.5}, RG, false, false},
    {"VIX", {NY, 0, 16.0}, {NY, 0, 17.0}, NV, false, X, X, NV, true, {NY, 0, 3.0}, {NY, 0, 9.5}, NV, false, false},
    {"RENDIMIENTO_CBOE", {NY, 0, 15.0}, {NY, 0, 17.0}, NV, false, X, X, NV, true, {NY, 0, 0.0}, {NY, 0, 9.5}, NV, false, false},
    {"TESORO", {NY, 0, 15.0}, {NY, 0, 17.0}, NV, true, {NY, 0, 13.0}, {NY, 0, 17.0}, NV, false, X, X, NV, false, false},
    {"BLOOMBERG_RENDIMIENTO", {NY, 0, 15.0}, {NY, 0, 17.0}, NV, false, X, X, NV, false, X, X, NV, false, false},
    {"CME_ENERGIA_METAL", {NY, 0, 16.0}, {NY, 0, 17.0}, ME, true, {NY, 0, 12.0}, {NY, 0, 17.0}, NV, true, {NY, -1, 18.0}, {NY, -1, 19.0}, ME, false, true},
    {"CBOT_GRANO", {NY, 0, 14.0}, {NY, 0, 15.0}, ME, true, {NY, 0, 12.0}, {NY, 0, 15.0}, NV, true, {NY, -1, 20.0}, {NY, -1, 21.0}, ME, false, true},
    {"CME_RV_BLOOMBERG", {NY, 0, 16.0}, {NY, 0, 17.0}, NV, true, {NY, 0, 12.0}, {NY, 0, 17.0}, NV, true, {NY, -1, 18.0}, {NY, -1, 18.0}, RG, false, true},
    {"DXY", {NY, 0, 14.0}, {NY, 0, 18.5}, NV, false, X, X, NV, false, X, X, NV, true, false},
    {"FX_YAHOO", {U, 1, -1.0}, {U, 1, 0.0}, ME, false, X, X, ME, false, X, X, ME, true, true},
    {"FX_YAHOO_ANALOGIA", {NY, 0, 17.0}, {U, 1, 0.0}, NV, false, X, X, NV, false, X, X, NV, true, true},
    {"CRIPTO_UTC", {U, 1, 0.0}, {U, 1, 0.0}, RG, false, X, X, RG, true, {U, 0, 0.0}, {U, 0, 0.0}, RG, false, true},
    {"CRIPTO_NY16", {NY, 0, 16.0}, {NY, 0, 16.0}, RG, false, X, X, RG, true, {NY, -1, 16.0}, {NY, -1, 16.0}, RG, false, true},
    {"CRIPTO_NO_VERIFICADO", {U, 1, -1.0}, {U, 1, 0.0}, NV, false, X, X, NV, false, X, X, NV, true, true},
    {"CRIPTO_BLOOMBERG", {NY, 0, 16.0}, {U, 1, 0.0}, NV, false, X, X, NV, false, X, X, NV, true, true},
    {"EUROPA_INDICE", {NY, 0, 11.0}, {NY, 0, 12.0}, ME, false, X, X, ME, true, {BER, 0, 9.0}, {BER, 0, 9.0}, RG, false, false},
    {"DESCONOCIDO", {U, 0, 0.0}, {U, 1, 0.0}, NV, false, X, X, NV, false, X, X, NV, true, false},
};

[[nodiscard]] inline const PerfilVentana& perfil_por_nombre(const char* n) noexcept {
    for (const auto& p : PERFILES)
        if (std::strcmp(p.nombre, n) == 0) return p;
    return PERFILES[sizeof(PERFILES) / sizeof(PERFILES[0]) - 1];
}

struct Asig { const char* sim; const char* perfil; };
inline const Asig SIMBOLO_PERFIL[] = {
    {"HYG", "NYSE_RV"}, {"IEF", "NYSE_RV"}, {"EEM", "NYSE_RV"}, {"KRE", "NYSE_RV"}, {"TIP", "NYSE_RV"}, {"VNQ", "NYSE_RV"},
    {"XLE", "NYSE_RV"}, {"XLP", "NYSE_RV"}, {"SPY", "NYSE_RV"}, {"IBIT", "NYSE_RV"}, {"LQD", "NYSE_RV"},
    {"SPX", "NYSE_INDICE"}, {"^GSPC", "NYSE_INDICE"}, {"NDX", "NYSE_INDICE"}, {"^NDX", "NYSE_INDICE"}, {"RTY", "NYSE_INDICE"},
    {"^RUT", "NYSE_INDICE"},
    {"VIX", "VIX"}, {"^VIX", "VIX"},
    {"US05Y", "RENDIMIENTO_CBOE"}, {"US10Y", "RENDIMIENTO_CBOE"}, {"US30Y", "RENDIMIENTO_CBOE"}, {"^FVX", "RENDIMIENTO_CBOE"},
    {"^TNX", "RENDIMIENTO_CBOE"}, {"^TYX", "RENDIMIENTO_CBOE"},
    {"REAL10Y", "TESORO"}, {"US03M", "TESORO"},
    {"USGG10YR", "BLOOMBERG_RENDIMIENTO"},
    {"CL", "CME_ENERGIA_METAL"}, {"NG", "CME_ENERGIA_METAL"}, {"BZ", "CME_ENERGIA_METAL"}, {"GC", "CME_ENERGIA_METAL"},
    {"HG", "CME_ENERGIA_METAL"},
    {"ZC", "CBOT_GRANO"}, {"ZW", "CBOT_GRANO"}, {"ZS", "CBOT_GRANO"},
    {"ES1", "CME_RV_BLOOMBERG"},
    {"DXY", "DXY"}, {"DX-Y.NYB", "DXY"},
    {"JPY", "FX_YAHOO"},
    {"USDCNH", "FX_YAHOO_ANALOGIA"},
    {"XBTUSD", "CRIPTO_BLOOMBERG"},
    {"SX5E", "EUROPA_INDICE"},
};

[[nodiscard]] inline const char* asignado(const std::string& s) noexcept {
    for (const auto& a : SIMBOLO_PERFIL)
        if (s == a.sim) return a.perfil;
    return nullptr;
}
}  // namespace detalle_ventana

// Perfil de ventana de un instrumento. `fuente` desambigua (US02Y); `corte` es el corte diario de cripto de Binance.
[[nodiscard]] inline const PerfilVentana& perfil_de(const char* simbolo, const char* fuente = nullptr,
                                                    CorteDiarioCripto corte = CorteDiarioCripto::UTC_00) {
    using namespace detalle_ventana;
    const std::string s = simbolo_limpio(simbolo);
    std::string f = fuente ? fuente : "";
    for (auto& ch : f) ch = static_cast<char>((ch >= 'a' && ch <= 'z') ? ch - 32 : ch);
    if (s.size() >= 4 && s.compare(s.size() - 4, 4, "USDT") == 0)
        return perfil_por_nombre(corte == CorteDiarioCripto::NY_16 ? "CRIPTO_NY16" : "CRIPTO_UTC");
    if (!s.empty() && s[0] == 'T' && s.size() >= 3 && s.compare(s.size() - 3, 3, "USD") == 0 && s.size() <= 8 && !asignado(s))
        return perfil_por_nombre("CRIPTO_NO_VERIFICADO");
    if (s == "BTC" || s == "ETH") return perfil_por_nombre("CRIPTO_NO_VERIFICADO");
    if (s == "US02Y") return perfil_por_nombre(f == "YAHOO_VIA_MACRO_D1" ? "RENDIMIENTO_CBOE" : "TESORO");
    const char* p = asignado(s);
    return perfil_por_nombre(p ? p : "DESCONOCIDO");
}

[[nodiscard]] inline int64_t a_utc_ms(uint32_t fecha_sesion, const Instante& i) noexcept {
    return utc_ms_desde_local(i.zona, fecha_sesion, i.dias, i.hora);
}

struct VentanaBarra {
    bool ok = false;                 // false: no hubo sesion en las dos semanas previas (Python lanza ValueError)
    uint32_t fecha_sesion = 0;       // la fila pedida
    uint32_t sesion_informativa = 0; // la sesion cuya informacion lleva la fila (otra si la fila no es sesion)
    bool sin_sesion = false;
    TipoDia tipo_dia = TipoDia::SESION_REGULAR;
    const PerfilVentana* perfil = nullptr;
    int64_t cierre_lo_ms = 0, cierre_hi_ms = 0;
    Procedencia proc_cierre = Procedencia::NO_VERIFICADO;
    bool tiene_inicio = false;
    int64_t inicio_lo_ms = 0, inicio_hi_ms = 0;
    Procedencia proc_inicio = Procedencia::NO_VERIFICADO;
    bool negociable = false;
};

namespace detalle_ventana {
inline void cierre_sesion(const PerfilVentana& p, uint32_t d, TipoDia tipo, int64_t& lo, int64_t& hi, Procedencia& proc) noexcept {
    if (tipo == TipoDia::SESION_REDUCIDA && p.tiene_cierre_reducido) {
        lo = a_utc_ms(d, p.red_lo); hi = a_utc_ms(d, p.red_hi); proc = p.proc_cierre_reducido;
    } else {
        lo = a_utc_ms(d, p.cierre_lo); hi = a_utc_ms(d, p.cierre_hi); proc = p.proc_cierre;
    }
}
}  // namespace detalle_ventana

// Ventana de informacion (UTC) de la barra diaria de `simbolo` en `fecha_sesion` (ventanas_informacion.ventana_barra).
[[nodiscard]] inline VentanaBarra ventana_barra(const char* simbolo, uint32_t fecha_sesion, const char* fuente = nullptr,
                                                CorteDiarioCripto corte = CorteDiarioCripto::UTC_00) {
    VentanaBarra v;
    const PerfilVentana& p = perfil_de(simbolo, fuente, corte);
    const std::string s = simbolo_limpio(simbolo);
    const Calendario cal = calendario_de(s.c_str(), fuente);
    const int64_t z = dias_desde_yyyymmdd(fecha_sesion);
    v.perfil = &p;
    v.fecha_sesion = fecha_sesion;
    v.negociable = p.negociable;
    v.tipo_dia = cal.clasificar(z).tipo;
    v.sin_sesion = !es_tipo_sesion(v.tipo_dia);
    int64_t zi = z;
    if (v.sin_sesion && !cal.sesion_previa(z, zi)) return v;
    const TipoDia tipo_info = v.sin_sesion ? cal.clasificar(zi).tipo : v.tipo_dia;
    v.sesion_informativa = yyyymmdd_desde_dias(zi);
    detalle_ventana::cierre_sesion(p, v.sesion_informativa, tipo_info, v.cierre_lo_ms, v.cierre_hi_ms, v.proc_cierre);
    if (p.inicio_es_cierre_previo) {
        int64_t zp;
        if (!cal.sesion_previa(zi, zp)) return v;
        detalle_ventana::cierre_sesion(p, yyyymmdd_desde_dias(zp), cal.clasificar(zp).tipo, v.inicio_lo_ms, v.inicio_hi_ms, v.proc_inicio);
        v.tiene_inicio = true;
    } else if (p.tiene_inicio) {
        v.inicio_lo_ms = a_utc_ms(v.sesion_informativa, p.ini_lo);
        v.inicio_hi_ms = a_utc_ms(v.sesion_informativa, p.ini_hi);
        v.proc_inicio = p.proc_inicio;
        v.tiene_inicio = true;
    } else {
        v.proc_inicio = p.proc_inicio;
    }
    v.ok = true;
    return v;
}

struct ResultadoGuarda {
    ClasificacionCausalPredictiva clase = ClasificacionCausalPredictiva::NO_VERIFICADO;
    double margen_s = std::numeric_limits<double>::quiet_NaN();   // inicio_lo - (info_hi + L); > 0 es holgura garantizada
};

// Nucleo de la guarda (ventanas_informacion.clasificar_ventanas). Latencia en segundos; negativa no tiene sentido.
[[nodiscard]] inline ResultadoGuarda clasificar_ventanas(int64_t info_lo_ms, int64_t info_hi_ms, int64_t ini_lo_ms, int64_t ini_hi_ms,
                                                         double latencia_s = LATENCIA_DIARIA_S) noexcept {
    ResultadoGuarda r;
    if (!(latencia_s >= 0.0)) return r;   // NO_VERIFICADO, margen NaN: Python lanza ValueError
    // pd.Timedelta(seconds=L) redondea a nanosegundos; aqui se compara en ns enteros para la misma frontera.
    const int64_t lat_ns = static_cast<int64_t>(std::llround(latencia_s * 1e9));
    const int64_t ilo = info_lo_ms * 1000000, ihi = info_hi_ms * 1000000, slo = ini_lo_ms * 1000000, shi = ini_hi_ms * 1000000;
    if (ihi + lat_ns < slo) r.clase = ClasificacionCausalPredictiva::EJECUTABLE;
    else if (ilo + lat_ns >= shi) r.clase = ClasificacionCausalPredictiva::REZAGO;
    else r.clase = ClasificacionCausalPredictiva::NO_VERIFICADO;
    r.margen_s = static_cast<double>(slo - (ihi + lat_ns)) / 1e9;
    return r;
}

enum class TipoRetorno : uint8_t { CIERRE_A_CIERRE = 0, APERTURA_A_CIERRE = 1 };

// ?Puede la barra diaria de X en fecha_x predecir de forma ejecutable el retorno de Y en fecha_y? (clasificar_par de Python).
// fecha_anterior_y = 0: la sesion previa de Y segun su calendario.
[[nodiscard]] inline ResultadoGuarda clasificar_par(const char* x, uint32_t fecha_x, const char* y, uint32_t fecha_y,
                                                    TipoRetorno tipo = TipoRetorno::CIERRE_A_CIERRE,
                                                    double latencia_s = LATENCIA_DIARIA_S, const char* fuente_x = nullptr,
                                                    const char* fuente_y = nullptr,
                                                    CorteDiarioCripto corte_x = CorteDiarioCripto::UTC_00,
                                                    CorteDiarioCripto corte_y = CorteDiarioCripto::UTC_00,
                                                    uint32_t fecha_anterior_y = 0) {
    ResultadoGuarda r;
    const VentanaBarra vx = ventana_barra(x, fecha_x, fuente_x, corte_x);
    const VentanaBarra vy = ventana_barra(y, fecha_y, fuente_y, corte_y);
    if (!vx.ok || !vy.ok) return r;
    if (vy.sin_sesion) { r.clase = ClasificacionCausalPredictiva::SIN_SESION; return r; }
    int64_t slo, shi;
    if (tipo == TipoRetorno::CIERRE_A_CIERRE) {
        uint32_t d_ant = fecha_anterior_y;
        if (d_ant == 0) {
            int64_t zp;
            if (!calendario_de(simbolo_limpio(y).c_str(), fuente_y).sesion_previa(dias_desde_yyyymmdd(vy.fecha_sesion), zp)) return r;
            d_ant = yyyymmdd_desde_dias(zp);
        }
        const VentanaBarra va = ventana_barra(y, d_ant, fuente_y, corte_y);
        if (!va.ok) return r;
        slo = va.cierre_lo_ms; shi = va.cierre_hi_ms;
    } else {
        if (!vy.tiene_inicio) { r.clase = ClasificacionCausalPredictiva::NO_DEFINIDO; return r; }
        slo = vy.inicio_lo_ms; shi = vy.inicio_hi_ms;
    }
    return clasificar_ventanas(vx.cierre_lo_ms, vx.cierre_hi_ms, slo, shi, latencia_s);
}

// Guarda intradia: las barras ya van estampadas a su cierre, asi que la informacion de X esta completa en su sello y el
// retorno de Y empieza en el sello de su barra anterior (clasificar_par_intradia de Python).
class GuardaAsincroniaCxx {
public:
    // Instantes puntuales (lo = hi). 0 = sin hora -> NO_VERIFICADO. Frontera ESTRICTA: margen <= L no es ejecutable.
    [[nodiscard]] static ClasificacionCausalPredictiva clasificar_par(uint64_t cierre_x_utc_ms, uint64_t inicio_retorno_y_utc_ms,
                                                                      double latencia_s = LATENCIA_INTRADIA_S) noexcept {
        if (cierre_x_utc_ms == 0 || inicio_retorno_y_utc_ms == 0) return ClasificacionCausalPredictiva::NO_VERIFICADO;
        const auto x = static_cast<int64_t>(cierre_x_utc_ms), y = static_cast<int64_t>(inicio_retorno_y_utc_ms);
        return clasificar_ventanas(x, x, y, y, latencia_s).clase;
    }
};

}  // namespace ark::dataset
