// calendario.hpp -- CALENDARIOS DE SESIONES POR ACTIVO (ARK-DATA-01-CXX, correccion del Jefe 2026-11-01)
// Reactor ARK / Cancion del Micelio
//
// Porte EXACTO de ARQUITECTO DE DATOS/recopilacion_unificada/calendarios.py (la referencia canonica del dataset): NYSE,
// SIFMA_BONOS, CME_ICE_EEUU, SIX_SUIZA, CRIPTO_24_7, FX_YAHOO y HABILES_SIN_VERIFICAR, con la misma clasificacion de cada dia
// (tipo_dia, nombre_dia, fraccion_sesion) y el mismo mapa activo -> calendario. La identidad dia a dia contra Python en
// 1971-2030 la comprueba tests/contraste_python.py.
//
// Lo que corrige frente a la version anterior (auditoria del Jefe, 2026-11-01):
//   - horario de verano de EE. UU. con su historia real (antes, la regla de 2007 para todos los anos: 1 040 dias habiles de
//     1971-2006 con el cierre UTC errado en 1 h);
//   - los 10 cierres extraordinarios de NYSE desde 1999 (11-S, Reagan, Ford, Sandy, Bush, Carter) ;
//   - Juneteenth desde 2022 (no 2021: NYSE abrio el 2021-06-18);
//   - la media jornada del 3 de julio y del 24 de diciembre solo de lunes a jueves (en viernes son el feriado observado; antes
//     16 dias salian feriado Y media jornada a la vez, con un cierre a las 17-18 UTC en un dia cerrado);
//   - un calendario por activo (antes todo TradFi era NYSE);
//   - fechas anteriores a 1970 sin _mkgmtime64 (devolvia -1).
// La regla de MLK se aplica, como en Python, a todos los anos: NYSE no cerro ese dia hasta 1998, y el calendario solo esta
// verificado desde 1999 (`verificado_desde`). Antes de 1999 es una regla, no una medicion.
#pragma once

#include "tipos_dataset.hpp"
#include "zona_horaria.hpp"
#include <cstdint>
#include <cstring>
#include <cmath>
#include <limits>
#include <string>

namespace ark::dataset {

enum class TipoDia : uint8_t {
    SESION_REGULAR = 0,
    SESION_REDUCIDA = 1,
    FERIADO = 2,
    CIERRE_EXTRAORDINARIO = 3,
    FIN_DE_SEMANA = 4
};

[[nodiscard]] constexpr const char* nombre_tipo_dia(TipoDia t) noexcept {
    switch (t) {
        case TipoDia::SESION_REGULAR: return "SESION_REGULAR";
        case TipoDia::SESION_REDUCIDA: return "SESION_REDUCIDA";
        case TipoDia::FERIADO: return "FERIADO";
        case TipoDia::CIERRE_EXTRAORDINARIO: return "CIERRE_EXTRAORDINARIO";
        default: return "FIN_DE_SEMANA";
    }
}

[[nodiscard]] constexpr bool es_tipo_sesion(TipoDia t) noexcept {
    return t == TipoDia::SESION_REGULAR || t == TipoDia::SESION_REDUCIDA;
}

enum class CalendarioId : uint8_t {
    NYSE = 0,
    SIFMA_BONOS = 1,
    CME_ICE_EEUU = 2,
    SIX_SUIZA = 3,
    CRIPTO_24_7 = 4,
    FX_YAHOO = 5,
    HABILES_SIN_VERIFICAR = 6
};

[[nodiscard]] constexpr const char* nombre_calendario(CalendarioId c) noexcept {
    switch (c) {
        case CalendarioId::NYSE: return "NYSE";
        case CalendarioId::SIFMA_BONOS: return "SIFMA_BONOS";
        case CalendarioId::CME_ICE_EEUU: return "CME_ICE_EEUU";
        case CalendarioId::SIX_SUIZA: return "SIX_SUIZA";
        case CalendarioId::CRIPTO_24_7: return "CRIPTO_24_7";
        case CalendarioId::FX_YAHOO: return "FX_YAHOO";
        default: return "HABILES_SIN_VERIFICAR";
    }
}

// Rango verificado de cada calendario (YYYYMMDD; 0 = no verificado contra ninguna serie; hasta 0 = la regla no caduca).
[[nodiscard]] constexpr uint32_t verificado_desde(CalendarioId c) noexcept {
    switch (c) {
        case CalendarioId::NYSE: return 19990101;
        case CalendarioId::SIFMA_BONOS: return 20060101;
        case CalendarioId::CME_ICE_EEUU: return 20070101;
        case CalendarioId::SIX_SUIZA: return 20070330;
        case CalendarioId::CRIPTO_24_7: return 20140917;
        case CalendarioId::FX_YAHOO: return 19990101;
        default: return 0;
    }
}
[[nodiscard]] constexpr uint32_t verificado_hasta(CalendarioId c) noexcept {
    return c == CalendarioId::SIFMA_BONOS ? 20261001 : 0;
}

struct ClaseDia {
    TipoDia tipo = TipoDia::SESION_REGULAR;
    const char* nombre = "";
    double fraccion = 1.0;   // NaN: la reduccion existe pero su duracion no esta modelada
};

namespace detalle_cal {

constexpr double NAN_D = std::numeric_limits<double>::quiet_NaN();
constexpr int LUN = 0, JUE = 3, DOM = 6;

[[nodiscard]] constexpr int anio_de(int64_t z) noexcept { return civil_desde_dias(z).y; }
[[nodiscard]] constexpr int64_t D(int y, int m, int d) noexcept { return dias_desde_civil(y, m, d); }

[[nodiscard]] inline int64_t pascua(int anio) noexcept {
    const int a = anio % 19, b = anio / 100, c = anio % 100, d = b / 4, e = b % 4;
    const int f = (b + 8) / 25, g = (b - f + 1) / 3;
    const int h = (19 * a + b - d - g + 15) % 30;
    const int i = c / 4, k = c % 4;
    const int l = (32 + 2 * e + 2 * i - h - k) % 7;
    const int m = (a + 11 * h + 22 * l) / 451;
    const int x = h + l - 7 * m + 114;
    return D(anio, x / 31, x % 31 + 1);
}

[[nodiscard]] inline int64_t observado_sab_vie_dom_lun(int64_t z) noexcept {
    const int w = dia_semana_lunes0(z);
    return w == 5 ? z - 1 : (w == 6 ? z + 1 : z);
}

// Un mapa pequeno fecha -> (nombre, fraccion) con la semantica de un dict de Python: asignar una clave ya presente cambia
// el valor y conserva la posicion.
struct Mapa {
    int64_t z[40];
    const char* nom[40];
    double frac[40];
    int n = 0;
    void poner(int64_t k, const char* s, double f = 0.0) noexcept {
        for (int i = 0; i < n; ++i)
            if (z[i] == k) { nom[i] = s; frac[i] = f; return; }
        if (n < 40) { z[n] = k; nom[n] = s; frac[n] = f; ++n; }
    }
    void quitar(int64_t k) noexcept {
        for (int i = 0; i < n; ++i)
            if (z[i] == k) {
                for (int j = i + 1; j < n; ++j) { z[j - 1] = z[j]; nom[j - 1] = nom[j]; frac[j - 1] = frac[j]; }
                --n;
                return;
            }
    }
    [[nodiscard]] int buscar(int64_t k) const noexcept {
        for (int i = 0; i < n; ++i)
            if (z[i] == k) return i;
        return -1;
    }
};

// -------------------------------------------------------------------------------------------------------------------- NYSE
struct Extra { int64_t z; const char* nom; };
inline const Extra EXTRA_NYSE[] = {
    {D(2001, 9, 11), "11-S"}, {D(2001, 9, 12), "11-S"}, {D(2001, 9, 13), "11-S"}, {D(2001, 9, 14), "11-S"},
    {D(2004, 6, 11), "Duelo nacional: Reagan"},
    {D(2007, 1, 2), "Duelo nacional: Ford"},
    {D(2012, 10, 29), "Huracan Sandy"}, {D(2012, 10, 30), "Huracan Sandy"},
    {D(2018, 12, 5), "Duelo nacional: G. H. W. Bush"},
    {D(2025, 1, 9), "Duelo nacional: Carter"},
};
inline const Extra EXTRA_SIFMA[] = {
    {D(2001, 9, 11), "11-S"}, {D(2001, 9, 12), "11-S"},
    {D(2004, 6, 11), "Duelo nacional: Reagan"},
    {D(2012, 10, 30), "Huracan Sandy"},
    {D(2018, 12, 5), "Duelo nacional: G. H. W. Bush"},
};
inline const int64_t ABIERTOS_CME_SOBRE_NYSE[] = {D(2001, 9, 14), D(2007, 1, 2), D(2025, 1, 9)};
inline const int VS_ABIERTO_SIFMA[] = {2007, 2010, 2012, 2015, 2021, 2023, 2026};

constexpr double NYSE_FRAC_MEDIA = 3.5 / 6.5;
constexpr double SIFMA_FRAC_ANTICIPADO = 6.0 / 9.0;
constexpr double SIFMA_FRAC_VS_NFP = 4.0 / 9.0;

template <size_t K>
[[nodiscard]] inline const char* en_extra(const Extra (&t)[K], int64_t z) noexcept {
    for (const auto& e : t)
        if (e.z == z) return e.nom;
    return nullptr;
}
[[nodiscard]] inline bool abierto_cme(int64_t z) noexcept {
    for (int64_t a : ABIERTOS_CME_SOBRE_NYSE)
        if (a == z) return true;
    return false;
}
[[nodiscard]] inline bool vs_abierto_sifma(int anio) noexcept {
    for (int a : VS_ABIERTO_SIFMA)
        if (a == anio) return true;
    return false;
}

inline void nombres_nyse(int y, Mapa& f) noexcept {
    const int64_t an = D(y, 1, 1);
    const int w = dia_semana_lunes0(an);
    if (w == 6) f.poner(an + 1, "Ano Nuevo (observado)");
    else if (w < 5) f.poner(an, "Ano Nuevo");
    f.poner(enesimo_dia_semana(y, 1, LUN, 3), "Martin Luther King Jr.");
    f.poner(enesimo_dia_semana(y, 2, LUN, 3), "Washington's Birthday");
    f.poner(pascua(y) - 2, "Viernes Santo");
    f.poner(enesimo_dia_semana(y, 5, LUN, -1), "Memorial Day");
    if (y >= 2022) f.poner(observado_sab_vie_dom_lun(D(y, 6, 19)), "Juneteenth");
    f.poner(observado_sab_vie_dom_lun(D(y, 7, 4)), "Independence Day");
    f.poner(enesimo_dia_semana(y, 9, LUN, 1), "Labor Day");
    f.poner(enesimo_dia_semana(y, 11, JUE, 4), "Thanksgiving");
    f.poner(observado_sab_vie_dom_lun(D(y, 12, 25)), "Navidad");
}

[[nodiscard]] inline bool feriado_nyse(int64_t z) noexcept {
    Mapa f;
    nombres_nyse(anio_de(z), f);
    return f.buscar(z) >= 0 || en_extra(EXTRA_NYSE, z) != nullptr;
}

inline void reducidas_nyse(int y, Mapa& r) noexcept {
    r.poner(enesimo_dia_semana(y, 11, JUE, 4) + 1, "Viernes despues de Thanksgiving (13:00 ET)", NYSE_FRAC_MEDIA);
    const int64_t j3 = D(y, 7, 3), d24 = D(y, 12, 24);
    if (dia_semana_lunes0(j3) <= 3) r.poner(j3, "Vispera de Independence Day (13:00 ET)", NYSE_FRAC_MEDIA);
    if (dia_semana_lunes0(d24) <= 3) r.poner(d24, "Nochebuena (13:00 ET)", NYSE_FRAC_MEDIA);
}

// ------------------------------------------------------------------------------------------------------------------- SIFMA
inline void nombres_sifma(int y, Mapa& f) noexcept {
    nombres_nyse(y, f);
    if (vs_abierto_sifma(y)) f.quitar(pascua(y) - 2);
    f.poner(enesimo_dia_semana(y, 10, LUN, 2), "Columbus Day");
    const int64_t vet = D(y, 11, 11);
    const int w = dia_semana_lunes0(vet);
    if (w == 6) f.poner(vet + 1, "Veterans Day (observado)");
    else if (w < 5) f.poner(vet, "Veterans Day");
}

[[nodiscard]] inline bool feriado_sifma(int64_t z) noexcept {
    Mapa f;
    nombres_sifma(anio_de(z), f);
    return f.buscar(z) >= 0 || en_extra(EXTRA_SIFMA, z) != nullptr;
}

[[nodiscard]] inline int64_t habil_anterior_sifma(int64_t z) noexcept {
    int64_t x = z - 1;
    while (dia_semana_lunes0(x) >= 5 || feriado_sifma(x)) --x;
    return x;
}

// Los nombres de las visperas se construyen con el del feriado: se guardan en un almacen estatico por cadena (pocos y fijos).
[[nodiscard]] inline const char* nombre_vispera(const char* feriado) noexcept {
    static const char* const BASES[] = {"Ano Nuevo", "Ano Nuevo (observado)", "Martin Luther King Jr.", "Washington's Birthday",
                                        "Viernes Santo", "Memorial Day", "Juneteenth", "Independence Day", "Labor Day",
                                        "Thanksgiving", "Navidad", "Columbus Day", "Veterans Day", "Veterans Day (observado)"};
    static const char* const VISP[] = {"Vispera de Ano Nuevo (14:00 ET)", "Vispera de Ano Nuevo (observado) (14:00 ET)",
                                       "Vispera de Martin Luther King Jr. (14:00 ET)", "Vispera de Washington's Birthday (14:00 ET)",
                                       "Vispera de Viernes Santo (14:00 ET)", "Vispera de Memorial Day (14:00 ET)",
                                       "Vispera de Juneteenth (14:00 ET)", "Vispera de Independence Day (14:00 ET)",
                                       "Vispera de Labor Day (14:00 ET)", "Vispera de Thanksgiving (14:00 ET)",
                                       "Vispera de Navidad (14:00 ET)", "Vispera de Columbus Day (14:00 ET)",
                                       "Vispera de Veterans Day (14:00 ET)", "Vispera de Veterans Day (observado) (14:00 ET)"};
    for (size_t i = 0; i < sizeof(BASES) / sizeof(BASES[0]); ++i)
        if (std::strcmp(BASES[i], feriado) == 0) return VISP[i];
    return "Vispera de feriado (14:00 ET)";
}

inline void reducidas_sifma(int y, Mapa& r) noexcept {
    Mapa a, b;
    nombres_sifma(y, a);
    nombres_sifma(y + 1, b);
    for (const Mapa* m : {&a, &b})
        for (int i = 0; i < m->n; ++i) {
            const int64_t v = habil_anterior_sifma(m->z[i]);
            if (anio_de(v) == y) r.poner(v, nombre_vispera(m->nom[i]), SIFMA_FRAC_ANTICIPADO);
        }
    const int64_t fijas[3] = {enesimo_dia_semana(y, 11, JUE, 4) + 1, D(y, 12, 24), D(y, 12, 31)};
    const char* nfij[3] = {"Viernes despues de Thanksgiving (14:00 ET)", "Nochebuena (14:00 ET)", "Nochevieja (14:00 ET)"};
    for (int i = 0; i < 3; ++i)
        if (dia_semana_lunes0(fijas[i]) < 5 && !feriado_sifma(fijas[i])) r.poner(fijas[i], nfij[i], SIFMA_FRAC_ANTICIPADO);
    if (vs_abierto_sifma(y)) r.poner(pascua(y) - 2, "Viernes Santo con NFP (12:00 ET)", SIFMA_FRAC_VS_NFP);
}

// --------------------------------------------------------------------------------------------------------------------- SIX
inline void nombres_six(int y, Mapa& f) noexcept {
    const int64_t p = pascua(y);
    Mapa t;
    t.poner(D(y, 1, 1), "Ano Nuevo");
    t.poner(D(y, 1, 2), "Berchtoldstag");
    t.poner(p - 2, "Viernes Santo");
    t.poner(p + 1, "Lunes de Pascua");
    t.poner(D(y, 5, 1), "Dia del Trabajo");
    t.poner(p + 39, "Ascension");
    t.poner(p + 50, "Lunes de Pentecostes");
    t.poner(D(y, 8, 1), "Fiesta nacional suiza");
    t.poner(D(y, 12, 25), "Navidad");
    t.poner(D(y, 12, 26), "San Esteban");
    if (y >= 2010) {
        t.poner(D(y, 12, 24), "Nochebuena");
        t.poner(D(y, 12, 31), "Nochevieja");
    }
    for (int i = 0; i < t.n; ++i)
        if (dia_semana_lunes0(t.z[i]) < 5) f.poner(t.z[i], t.nom[i]);
}

// ---------------------------------------------------------------------------------------------------------------------- FX
[[nodiscard]] inline bool es_sesion_fx_yahoo(int64_t z) noexcept {
    const int y = anio_de(z);
    const int64_t ini = enesimo_dia_semana(y, 3, DOM, -1), fin = enesimo_dia_semana(y, 10, DOM, -1);
    const int64_t manana = z + 1;
    if (ini <= manana && manana < fin) return dia_semana_lunes0(manana) < 5;
    return dia_semana_lunes0(z) < 5;
}

}  // namespace detalle_cal

// ======================================================================================================== Calendario
class Calendario {
public:
    explicit constexpr Calendario(CalendarioId id = CalendarioId::HABILES_SIN_VERIFICAR) noexcept : id_(id) {}
    [[nodiscard]] constexpr CalendarioId id() const noexcept { return id_; }
    [[nodiscard]] constexpr const char* nombre() const noexcept { return nombre_calendario(id_); }

    [[nodiscard]] bool verificado_en(uint32_t ini, uint32_t fin) const noexcept {
        const uint32_t d = verificado_desde(id_), h = verificado_hasta(id_);
        if (d == 0 || ini < d) return false;
        return h == 0 || fin <= h;
    }

    [[nodiscard]] bool es_dia_operable(int64_t z) const noexcept {
        switch (id_) {
            case CalendarioId::CRIPTO_24_7: return true;
            case CalendarioId::FX_YAHOO: return detalle_cal::es_sesion_fx_yahoo(z);
            default: return dia_semana_lunes0(z) < 5;
        }
    }

    [[nodiscard]] bool es_feriado(int64_t z) const noexcept {
        using namespace detalle_cal;
        switch (id_) {
            case CalendarioId::NYSE: return feriado_nyse(z);
            case CalendarioId::SIFMA_BONOS: return feriado_sifma(z);
            case CalendarioId::CME_ICE_EEUU: return feriado_nyse(z) && !abierto_cme(z);
            case CalendarioId::SIX_SUIZA: {
                Mapa f;
                nombres_six(anio_de(z), f);
                return f.buscar(z) >= 0;
            }
            default: return false;
        }
    }

    // Igual que Calendario.clasificar de Python para un dia.
    [[nodiscard]] ClaseDia clasificar(int64_t z) const noexcept {
        using namespace detalle_cal;
        ClaseDia c;
        if (!es_dia_operable(z)) return {TipoDia::FIN_DE_SEMANA, "", 0.0};
        const int y = anio_de(z);
        if (es_feriado(z)) {
            Mapa nom;
            const char* extra = nullptr;
            switch (id_) {
                case CalendarioId::NYSE: nombres_nyse(y, nom); extra = en_extra(EXTRA_NYSE, z); break;
                case CalendarioId::SIFMA_BONOS: nombres_sifma(y, nom); extra = en_extra(EXTRA_SIFMA, z); break;
                case CalendarioId::CME_ICE_EEUU:
                    nombres_nyse(y, nom);
                    extra = abierto_cme(z) ? nullptr : en_extra(EXTRA_NYSE, z);
                    break;
                case CalendarioId::SIX_SUIZA: nombres_six(y, nom); break;
                default: break;
            }
            const int k = nom.buscar(z);
            if (extra != nullptr && k < 0) return {TipoDia::CIERRE_EXTRAORDINARIO, extra, 0.0};
            return {TipoDia::FERIADO, k >= 0 ? nom.nom[k] : "Feriado", 0.0};
        }
        Mapa red;
        switch (id_) {
            case CalendarioId::NYSE: reducidas_nyse(y, red); break;
            case CalendarioId::SIFMA_BONOS: reducidas_sifma(y, red); break;
            case CalendarioId::CME_ICE_EEUU: {
                Mapa r;
                reducidas_nyse(y, r);
                for (int i = 0; i < r.n; ++i)
                    if (!(feriado_nyse(r.z[i]) && !abierto_cme(r.z[i]))) red.poner(r.z[i], nombre_cme(r.nom[i]), NAN_D);
                break;
            }
            default: break;
        }
        const int k = red.buscar(z);
        if (k >= 0) return {TipoDia::SESION_REDUCIDA, red.nom[k], red.frac[k]};
        return c;
    }

    [[nodiscard]] ClaseDia clasificar_fecha(uint32_t f) const noexcept { return clasificar(dias_desde_yyyymmdd(f)); }

    // Sesiones del calendario en [ini, fin] (dias desde 1970).
    [[nodiscard]] int64_t contar_sesiones(int64_t ini, int64_t fin) const noexcept {
        int64_t n = 0;
        for (int64_t z = ini; z <= fin; ++z)
            if (es_dia_operable(z) && !es_feriado(z)) ++n;
        return n;
    }

    // Ultima sesion (dia operable y no feriado) estrictamente anterior a z, en las dos semanas previas (como en Python).
    [[nodiscard]] bool sesion_previa(int64_t z, int64_t& out) const noexcept {
        for (int k = 1; k < 15; ++k)
            if (es_tipo_sesion(clasificar(z - k).tipo)) { out = z - k; return true; }
        return false;
    }

private:
    CalendarioId id_;

    [[nodiscard]] static const char* nombre_cme(const char* nyse) noexcept {
        if (std::strcmp(nyse, "Viernes despues de Thanksgiving (13:00 ET)") == 0) return "Viernes despues de Thanksgiving (cierre anticipado)";
        if (std::strcmp(nyse, "Vispera de Independence Day (13:00 ET)") == 0) return "Vispera de Independence Day (cierre anticipado)";
        if (std::strcmp(nyse, "Nochebuena (13:00 ET)") == 0) return "Nochebuena (cierre anticipado)";
        return nyse;
    }
};

// Simbolo canonico: sin los sufijos de Bloomberg y en mayusculas (calendarios.calendario_de / ventanas_informacion._limpio).
[[nodiscard]] inline std::string simbolo_limpio(const char* s) {
    std::string x = s ? s : "";
    for (const char* suf : {" Index", " Comdty", " Curncy"}) {
        size_t p;
        while ((p = x.find(suf)) != std::string::npos) x.erase(p, std::strlen(suf));
    }
    size_t a = 0, b = x.size();
    while (a < b && (x[a] == ' ' || x[a] == '\t' || x[a] == '\r' || x[a] == '\n')) ++a;
    while (b > a && (x[b - 1] == ' ' || x[b - 1] == '\t' || x[b - 1] == '\r' || x[b - 1] == '\n')) --b;
    x = x.substr(a, b - a);
    for (auto& ch : x) ch = static_cast<char>((ch >= 'a' && ch <= 'z') ? ch - 32 : ch);
    return x;
}

[[nodiscard]] inline Calendario calendario_de(const char* simbolo, const char* fuente = nullptr) {
    const std::string s = simbolo_limpio(simbolo);
    const std::string f = fuente ? fuente : "";
    if (s == "US02Y") return Calendario(f == "YAHOO_VIA_MACRO_D1" ? CalendarioId::NYSE : CalendarioId::SIFMA_BONOS);
    const bool fin_usdt = s.size() >= 4 && s.compare(s.size() - 4, 4, "USDT") == 0;
    const bool fin_usd = s.size() >= 3 && s.compare(s.size() - 3, 3, "USD") == 0;
    if (fin_usdt || (!s.empty() && s[0] == 'T' && fin_usd && s.size() <= 8)) return Calendario(CalendarioId::CRIPTO_24_7);
    static const char* const NYSE_[] = {"SPX", "^GSPC", "NDX", "^NDX", "RTY", "^RUT", "VIX", "^VIX", "HYG", "IEF", "EEM", "KRE",
                                        "TIP", "VNQ", "XLE", "XLP", "SPY", "IBIT", "US05Y", "US10Y", "US30Y", "USGG10YR"};
    static const char* const SIFMA_[] = {"REAL10Y", "US03M"};
    static const char* const CME_[] = {"DXY", "DX-Y.NYB", "ES1", "CL", "GC", "HG", "NG", "BZ", "ZC", "ZW", "ZS"};
    static const char* const CRIPTO_[] = {"BTC", "ETH", "XBTUSD"};
    static const char* const FX_[] = {"JPY", "USDCNH"};
    for (const char* x : NYSE_) if (s == x) return Calendario(CalendarioId::NYSE);
    for (const char* x : SIFMA_) if (s == x) return Calendario(CalendarioId::SIFMA_BONOS);
    for (const char* x : CME_) if (s == x) return Calendario(CalendarioId::CME_ICE_EEUU);
    if (s == "SX5E") return Calendario(CalendarioId::SIX_SUIZA);
    for (const char* x : CRIPTO_) if (s == x) return Calendario(CalendarioId::CRIPTO_24_7);
    for (const char* x : FX_) if (s == x) return Calendario(CalendarioId::FX_YAHOO);
    return Calendario(CalendarioId::HABILES_SIN_VERIFICAR);
}

// ====================================================================================== API anterior (compatibilidad)
// La conservan ark.cpp y las pruebas DS2-DS3. Ahora son la regla de NYSE ya corregida; para un activo concreto, usar
// calendario_de() y ventana_barra() (guardas_causales.hpp), que llevan el horario de cada mercado.
class MotorCalendarioTradFi {
public:
    // 0 = domingo ... 6 = sabado (convencion anterior de esta clase).
    [[nodiscard]] static constexpr int dia_semana(uint32_t yyyymmdd) noexcept {
        return (dia_semana_lunes0(dias_desde_yyyymmdd(yyyymmdd)) + 1) % 7;
    }
    [[nodiscard]] static constexpr bool es_fin_de_semana(uint32_t yyyymmdd) noexcept {
        const int ds = dia_semana(yyyymmdd);
        return ds == 0 || ds == 6;
    }
    [[nodiscard]] static uint32_t viernes_santo(int anio) noexcept {
        return yyyymmdd_desde_dias(detalle_cal::pascua(anio) - 2);
    }
    // Verano de EE. UU. a la hora del cierre regular (16:00 ET), con la historia de America/New_York.
    [[nodiscard]] static bool es_horario_verano_eeuu(uint32_t yyyymmdd) noexcept {
        const int64_t pared = dias_desde_yyyymmdd(yyyymmdd) * 86400 + 16 * 3600;
        return offset_s(Zona::NUEVA_YORK, pared) == -4 * 3600;
    }
    // Feriado o cierre extraordinario de NYSE.
    [[nodiscard]] static bool es_feriado_nyse(uint32_t yyyymmdd, const char** out_nombre = nullptr) noexcept {
        const ClaseDia c = Calendario(CalendarioId::NYSE).clasificar_fecha(yyyymmdd);
        const bool fer = c.tipo == TipoDia::FERIADO || c.tipo == TipoDia::CIERRE_EXTRAORDINARIO;
        if (fer && out_nombre) *out_nombre = c.nombre;
        return fer;
    }
    [[nodiscard]] static bool es_media_jornada(uint32_t yyyymmdd) noexcept {
        return Calendario(CalendarioId::NYSE).clasificar_fecha(yyyymmdd).tipo == TipoDia::SESION_REDUCIDA;
    }
    // Cierre por la regla de NYSE (16:00 ET, 13:00 ET en media jornada) o 00:00 UTC de D+1 en cripto. Devuelve 0 si el
    // instante es anterior a 1970 (no cabe en uint64): 0 significa "sin hora" y la guarda lo trata como NO_VERIFICADO.
    [[nodiscard]] static uint64_t calcular_cierre_real_utc_ms(uint32_t yyyymmdd, Venue venue) noexcept {
        int64_t t;
        if (venue == Venue::BINANCE || venue == Venue::BITFINEX) t = utc_ms_medianoche(yyyymmdd) + 86400000LL;
        else t = utc_ms_desde_local(Zona::NUEVA_YORK, yyyymmdd, 0, es_media_jornada(yyyymmdd) ? 13.0 : 16.0);
        return t > 0 ? static_cast<uint64_t>(t) : 0;
    }
    [[nodiscard]] static uint64_t calcular_apertura_real_utc_ms(uint32_t yyyymmdd, Venue venue) noexcept {
        int64_t t;
        if (venue == Venue::BINANCE || venue == Venue::BITFINEX) t = utc_ms_medianoche(yyyymmdd);
        else t = utc_ms_desde_local(Zona::NUEVA_YORK, yyyymmdd, 0, 9.5);
        return t > 0 ? static_cast<uint64_t>(t) : 0;
    }
};

}  // namespace ark::dataset
