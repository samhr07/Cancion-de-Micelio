// fase_calendario.hpp -- fase de la purga diurna desde el calendario del activo (ARK-DATA-01-CXX), C++20, sin memoria dinamica.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04, por orden del Comandante ("conectelos"): la fase de la purga deja de armarla
// el llamador a mano ("dia >= 5 = fin de semana", sin feriados) y sale de calendario.hpp con la hora LOCAL del mercado.
//
// fase = tipo * bins_dia + bin_hora, con
//   tipo     0 = dia laborable del calendario (sesion regular o reducida, de lunes a viernes), 1 = no laborable (fin de semana,
//            feriado, cierre extraordinario; en CRIPTO_24_7 el sabado y el domingo tambien, aunque haya sesion: su perfil de
//            actividad es otro). Declarado: las sesiones reducidas van con las laborables (~3 al ano: un perfil propio no maduraria);
//   bin_hora floor(segundo local del dia * bins_dia / 86400), en la zona del mercado (zona_de).
// dia = numero de dia LOCAL (dias desde 1970-01-01 en la zona), monotono: es el que la purga usa para comprometer cada dia.
// La conversion UTC -> local es exacta (las fronteras del horario de verano se fijan en UTC: Nueva York cambia a las 02:00 locales,
// 07:00 UTC en primavera y 06:00 UTC en otono; Berlin a las 01:00 UTC), sin la ambiguedad de pared de la hora repetida.
#pragma once

#include "ark/dataset/calendario.hpp"
#include "ark/dataset/zona_horaria.hpp"

#include <cstddef>
#include <cstdint>

namespace ark::dataset {

// Zona del horario de cada calendario. SIX usa la regla de Berlin (Zurich comparte la de la UE). CRIPTO y los habiles sin verificar
// van en UTC. FX_YAHOO en Nueva York (el corte diario del mercado de divisas es a las 17:00 ET).
[[nodiscard]] constexpr Zona zona_de(CalendarioId c) noexcept {
    switch (c) {
        case CalendarioId::NYSE:
        case CalendarioId::SIFMA_BONOS:
        case CalendarioId::CME_ICE_EEUU:
        case CalendarioId::FX_YAHOO: return Zona::NUEVA_YORK;
        case CalendarioId::SIX_SUIZA: return Zona::BERLIN;
        default: return Zona::UTC;
    }
}

// Desplazamiento (s) de la zona en el INSTANTE utc_s (hora local = UTC + desplazamiento).
[[nodiscard]] inline int64_t offset_desde_utc_s(Zona z, int64_t utc_s) noexcept {
    if (z == Zona::UTC) return 0;
    const int64_t dia = (utc_s >= 0 ? utc_s : utc_s - 86399) / 86400;
    const int anio = civil_desde_dias(dia).y;
    int64_t ini = 0, fin = 0;
    if (z == Zona::NUEVA_YORK) {
        constexpr int64_t EST = -5 * 3600, EDT = -4 * 3600;
        if (!detalle_zona::verano_ny(anio, ini, fin)) return EST;
        const int64_t t_ini = ini * 86400 + 7 * 3600;   // 02:00 EST = 07:00 UTC
        const int64_t t_fin = fin * 86400 + 6 * 3600;   // 02:00 EDT = 06:00 UTC
        return (utc_s >= t_ini && utc_s < t_fin) ? EDT : EST;
    }
    constexpr int64_t CET = 3600, CEST = 2 * 3600;
    if (!detalle_zona::verano_berlin(anio, ini, fin)) return CET;
    const int64_t t_ini = ini * 86400 + 1 * 3600;       // 01:00 UTC en los dos sentidos
    const int64_t t_fin = fin * 86400 + 1 * 3600;
    return (utc_s >= t_ini && utc_s < t_fin) ? CEST : CET;
}

struct FaseCalendarioCxx {
    bool activa = false;                                     // false: la fase y el dia los entrega el llamador
    CalendarioId calendario = CalendarioId::HABILES_SIN_VERIFICAR;
    Zona zona = Zona::UTC;                                   // normalmente zona_de(calendario)
    std::size_t bins_dia = 24;                               // fases por dia y por tipo de dia
    [[nodiscard]] std::size_t bins_fase() const noexcept { return 2 * bins_dia; }   // lo que pide la purga (bins_fase)
};

// Configuracion por omision para un simbolo: su calendario (calendario_de) y la zona de ese calendario.
[[nodiscard]] inline FaseCalendarioCxx fase_calendario_de(const char* simbolo, std::size_t bins_dia = 24, const char* fuente = nullptr) {
    FaseCalendarioCxx f{};
    f.activa = true;
    f.calendario = calendario_de(simbolo, fuente).id();
    f.zona = zona_de(f.calendario);
    f.bins_dia = bins_dia ? bins_dia : 1;
    return f;
}

struct FaseDiaCxx {
    std::size_t fase = 0;
    int64_t dia = 0;          // dia local desde 1970-01-01
    bool laborable = true;
    TipoDia tipo = TipoDia::SESION_REGULAR;
};

// Fase del instante utc_ms (marca de cierre de la barra). Sin reservas.
[[nodiscard]] inline FaseDiaCxx fase_de(const FaseCalendarioCxx& f, int64_t utc_ms) noexcept {
    const int64_t utc_s = (utc_ms >= 0 ? utc_ms : utc_ms - 999) / 1000;
    const int64_t local_s = utc_s + offset_desde_utc_s(f.zona, utc_s);
    FaseDiaCxx r{};
    r.dia = (local_s >= 0 ? local_s : local_s - 86399) / 86400;
    const int64_t seg = local_s - r.dia * 86400;
    const std::size_t bins = f.bins_dia ? f.bins_dia : 1;
    std::size_t bin = static_cast<std::size_t>((static_cast<uint64_t>(seg) * bins) / 86400u);
    if (bin >= bins) bin = bins - 1;
    r.tipo = Calendario(f.calendario).clasificar(r.dia).tipo;
    r.laborable = es_tipo_sesion(r.tipo) && dia_semana_lunes0(r.dia) < 5;
    r.fase = (r.laborable ? 0u : 1u) * bins + bin;
    return r;
}

}  // namespace ark::dataset
