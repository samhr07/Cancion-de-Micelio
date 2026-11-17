// zona_horaria.hpp -- ARITMETICA CIVIL Y HUSOS HORARIOS SIN EL CRT (ARK-DATA-01-CXX, correccion del Jefe 2026-11-01)
// Reactor ARK / Cancion del Micelio
//
// Sustituye a _mkgmtime64 (que devuelve -1 antes de 1970 y, convertido a uint64, daba 18446744073698751616) y a la regla
// de verano de EE. UU. de 2007 aplicada a todos los anos (1 040 dias habiles de 1971-2006 con el cierre UTC errado en 1 h).
//
//   dias_desde_civil / civil_desde_dias : algoritmo de H. Hinnant (calendario gregoriano proleptico, enteros, exacto).
//   offset_s(zona, segundos_de_pared)   : desplazamiento UTC de una hora LOCAL de pared, con la semantica fold = 0 de
//                                         zoneinfo (hora inexistente -> desplazamiento previo al cambio; hora ambigua -> la
//                                         primera, la de verano). Es la que usa ventanas_informacion._a_utc en Python.
//
// Zonas: America/New_York (historia de tzdata desde 1967), Europe/Berlin (sin verano 1950-1979; reglas de la UE desde 1980) y
// UTC. Verificadas contra zoneinfo hora a hora en 1970-2040 por tests/contraste_python.py.
#pragma once

#include <cstdint>
#include <cmath>

namespace ark::dataset {

// ---------------------------------------------------------------------------------------------------------- fechas civiles
[[nodiscard]] constexpr int64_t dias_desde_civil(int64_t y, int64_t m, int64_t d) noexcept {
    y -= m <= 2;
    const int64_t era = (y >= 0 ? y : y - 399) / 400;
    const int64_t yoe = y - era * 400;
    const int64_t doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;
    const int64_t doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    return era * 146097 + doe - 719468;
}

struct FechaCivil {
    int y, m, d;
};

[[nodiscard]] constexpr FechaCivil civil_desde_dias(int64_t z) noexcept {
    z += 719468;
    const int64_t era = (z >= 0 ? z : z - 146096) / 146097;
    const int64_t doe = z - era * 146097;
    const int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    const int64_t y = yoe + era * 400;
    const int64_t doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    const int64_t mp = (5 * doy + 2) / 153;
    const int64_t d = doy - (153 * mp + 2) / 5 + 1;
    const int64_t m = mp + (mp < 10 ? 3 : -9);
    return {static_cast<int>(y + (m <= 2)), static_cast<int>(m), static_cast<int>(d)};
}

[[nodiscard]] constexpr int64_t dias_desde_yyyymmdd(uint32_t f) noexcept {
    return dias_desde_civil(f / 10000, (f / 100) % 100, f % 100);
}

[[nodiscard]] constexpr uint32_t yyyymmdd_desde_dias(int64_t z) noexcept {
    const FechaCivil c = civil_desde_dias(z);
    return static_cast<uint32_t>(c.y * 10000 + c.m * 100 + c.d);
}

// 0 = lunes ... 6 = domingo (convencion de datetime.weekday de Python).
[[nodiscard]] constexpr int dia_semana_lunes0(int64_t z) noexcept {
    const int64_t w = (z + 3) % 7;   // 1970-01-01 fue jueves (3)
    return static_cast<int>(w >= 0 ? w : w + 7);
}

[[nodiscard]] constexpr bool fecha_valida(uint32_t f) noexcept {
    const int y = static_cast<int>(f / 10000), m = static_cast<int>((f / 100) % 100), d = static_cast<int>(f % 100);
    if (y < 1 || m < 1 || m > 12 || d < 1) return false;
    const FechaCivil c = civil_desde_dias(dias_desde_civil(y, m, d));
    return c.y == y && c.m == m && c.d == d;
}

// n-esimo `dow` (0 = lunes) del mes; n = -1 es el ultimo. Igual que calendarios.enesimo_dia_semana.
[[nodiscard]] constexpr int64_t enesimo_dia_semana(int anio, int mes, int dow, int n) noexcept {
    if (n > 0) {
        const int64_t p = dias_desde_civil(anio, mes, 1);
        const int64_t d = p + ((dow - dia_semana_lunes0(p)) % 7 + 7) % 7;
        return d + 7 * (n - 1);
    }
    const int64_t sig = dias_desde_civil(anio + (mes == 12 ? 1 : 0), mes % 12 + 1, 1);
    const int64_t ult = sig - 1;
    return ult - ((dia_semana_lunes0(ult) - dow) % 7 + 7) % 7;
}

// ---------------------------------------------------------------------------------------------------------------- husos
enum class Zona : uint8_t { UTC = 0, NUEVA_YORK = 1, BERLIN = 2 };

namespace detalle_zona {

// Inicio y fin del horario de verano de America/New_York en hora LOCAL ESTANDAR de pared (segundos desde 1970-01-01
// 00:00 local): el cambio es a las 02:00 locales en los dos sentidos. Devuelve false si el ano no tiene verano (antes de 1967).
inline bool verano_ny(int anio, int64_t& ini_dia, int64_t& fin_dia) noexcept {
    constexpr int DOM = 6;
    if (anio < 1967) return false;
    if (anio == 1974) { ini_dia = dias_desde_civil(1974, 1, 6); fin_dia = enesimo_dia_semana(1974, 10, DOM, -1); return true; }
    if (anio == 1975) { ini_dia = dias_desde_civil(1975, 2, 23); fin_dia = enesimo_dia_semana(1975, 10, DOM, -1); return true; }
    if (anio <= 1986) { ini_dia = enesimo_dia_semana(anio, 4, DOM, -1); fin_dia = enesimo_dia_semana(anio, 10, DOM, -1); return true; }
    if (anio <= 2006) { ini_dia = enesimo_dia_semana(anio, 4, DOM, 1); fin_dia = enesimo_dia_semana(anio, 10, DOM, -1); return true; }
    ini_dia = enesimo_dia_semana(anio, 3, DOM, 2);
    fin_dia = enesimo_dia_semana(anio, 11, DOM, 1);
    return true;
}

// Europe/Berlin: sin verano entre 1950 y 1979; 1980 del primer domingo de abril al ultimo de septiembre; UE desde 1981: del
// ultimo domingo de marzo al ultimo de septiembre (hasta 1995) o de octubre (desde 1996). Cambio a las 01:00 UTC.
inline bool verano_berlin(int anio, int64_t& ini_dia, int64_t& fin_dia) noexcept {
    constexpr int DOM = 6;
    if (anio < 1980) return false;
    if (anio == 1980) { ini_dia = enesimo_dia_semana(1980, 4, DOM, 1); fin_dia = enesimo_dia_semana(1980, 9, DOM, -1); return true; }
    ini_dia = enesimo_dia_semana(anio, 3, DOM, -1);
    fin_dia = enesimo_dia_semana(anio, anio <= 1995 ? 9 : 10, DOM, -1);
    return true;
}

}  // namespace detalle_zona

// Desplazamiento UTC (segundos, hora local = UTC + offset) de la hora de pared `pared_s` (segundos desde 1970-01-01 00:00 de
// pared), con fold = 0. Floor de la division para fechas anteriores a 1970.
[[nodiscard]] inline int64_t offset_s(Zona z, int64_t pared_s) noexcept {
    if (z == Zona::UTC) return 0;
    const int64_t dia = (pared_s >= 0 ? pared_s : pared_s - 86399) / 86400;
    const int anio = civil_desde_dias(dia).y;
    int64_t ini = 0, fin = 0;
    if (z == Zona::NUEVA_YORK) {
        constexpr int64_t EST = -5 * 3600, EDT = -4 * 3600;
        if (!detalle_zona::verano_ny(anio, ini, fin)) return EST;
        // Primavera: de 02:00 se salta a 03:00 (las de 02:00-03:00 no existen; fold = 0 -> EST). Otono: 01:00-02:00 se repite
        // (fold = 0 -> la primera, EDT).
        const int64_t t_ini = ini * 86400 + 3 * 3600;   // primera hora de pared ya en verano
        const int64_t t_fin = fin * 86400 + 2 * 3600;   // primera hora de pared de nuevo en invierno
        return (pared_s >= t_ini && pared_s < t_fin) ? EDT : EST;
    }
    // Berlin: CET = +1 h, CEST = +2 h; cambio a las 01:00 UTC = 02:00 CET (primavera) y 03:00 CEST (otono).
    constexpr int64_t CET = 3600, CEST = 2 * 3600;
    if (!detalle_zona::verano_berlin(anio, ini, fin)) return CET;
    const int64_t t_ini = ini * 86400 + 3 * 3600;   // 02:00-03:00 CET no existe -> fold = 0 -> CET
    const int64_t t_fin = fin * 86400 + 3 * 3600;   // 02:00-03:00 CEST se repite -> fold = 0 -> CEST
    return (pared_s >= t_ini && pared_s < t_fin) ? CEST : CET;
}

// Instante UTC (ms) de la hora local `hora` (decimal; negativa = antes de la medianoche) del dia `fecha + dias`, con la
// aritmetica de ventanas_informacion._a_utc: medianoche de pared + round(hora * 3600) s, y despues a UTC con fold = 0.
[[nodiscard]] inline int64_t utc_ms_desde_local(Zona z, uint32_t fecha_sesion, int dias, double hora) noexcept {
    const int64_t base = dias_desde_yyyymmdd(fecha_sesion) + dias;
    const int64_t seg = static_cast<int64_t>(std::nearbyint(hora * 3600.0));   // round-half-even, como round() de Python
    const int64_t pared = base * 86400 + seg;
    return (pared - offset_s(z, pared)) * 1000;
}

// Medianoche UTC (ms) de una fecha YYYYMMDD; con signo (fechas anteriores a 1970 dan negativo, no -1).
[[nodiscard]] constexpr int64_t utc_ms_medianoche(uint32_t f) noexcept { return dias_desde_yyyymmdd(f) * 86400000LL; }

}  // namespace ark::dataset
