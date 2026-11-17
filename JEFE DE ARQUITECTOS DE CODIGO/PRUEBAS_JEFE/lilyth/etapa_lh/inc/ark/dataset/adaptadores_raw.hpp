// adaptadores_raw.hpp -- LECTOR DE CSV DIARIO Y REJILLA DIARIA CAUSAL (ARK-DATA-01-CXX)
// Reactor ARK / Cancion del Micelio -- ARQUITECTO DE DATOS (Flash) y correccion del Jefe 2026-11-01
//
// 1. LectorCSVDataset::leer_csv_macro: la fila CRUDA de la fuente (cachas/macro_d1). Columnas por NOMBRE de la cabecera (antes
//    por posicion), fecha validada, campo vacio o no numerico = NaN contado (antes 0.0 en silencio), estampa R1 (antes 0),
//    es_imputado si el cierre no es valido y huella SHA-256 binaria de la barra.
// 2. construir_rejilla_diaria: la tabla diaria del dataset, con la misma secuencia que motor_ingesta._procesar_crudos:
//      a. rejilla: inserta los dias operables del calendario que la fuente no trae (feriados, cierres, sesiones sin dato);
//      b. anotacion: tipo_dia, fraccion_sesion, hay_dato, es_sesion_con_dato, es_imputado;
//      c. supervivencia (catalogo de ciclo de vida);
//      d. dias y sesiones acumulados entre sesiones con dato (convencion viernes -> lunes) y cierre_real_utc por perfil;
//      e. blindaje: todo precio invalido (NaN, inf, <= 0) a NaN, prohibida la censura izquierda, arrastre causal (ffill) y
//         es_imputado |= precio invalido;
//      f. supervivencia otra vez (corta el arrastre despues del desliste).
//    La identidad celda a celda contra el Python real la comprueba tests/contraste_python.py sobre cachas/macro_d1.
// 3. SerializadorBinarioDataset: volcado plano de BarraCausal (version 2: la estructura cambio).
#pragma once

#include "tipos_dataset.hpp"
#include "zona_horaria.hpp"
#include "calendario.hpp"
#include "guardas_causales.hpp"
#include "libro_hashes.hpp"
#include "supervivencia.hpp"

#include <windows.h>
#include <string>
#include <cstdint>
#include <cmath>
#include <vector>
#include <cstdio>
#include <charconv>
#include <algorithm>

namespace ark::dataset {

// =======================================================================
// LectorMemoryMappedWin32: RAII para mapeo directo de archivos en Windows
// =======================================================================
class LectorMemoryMappedWin32 {
private:
    HANDLE h_file_{INVALID_HANDLE_VALUE};
    HANDLE h_map_{nullptr};
    const char* data_{nullptr};
    size_t size_{0};

public:
    LectorMemoryMappedWin32() = default;
    ~LectorMemoryMappedWin32() { cerrar(); }

    LectorMemoryMappedWin32(const LectorMemoryMappedWin32&) = delete;
    LectorMemoryMappedWin32& operator=(const LectorMemoryMappedWin32&) = delete;

    bool abrir(const std::string& ruta) {
        cerrar();
        h_file_ = CreateFileA(ruta.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
        if (h_file_ == INVALID_HANDLE_VALUE) return false;
        LARGE_INTEGER li_size;
        if (!GetFileSizeEx(h_file_, &li_size) || li_size.QuadPart == 0) { cerrar(); return false; }
        size_ = static_cast<size_t>(li_size.QuadPart);
        h_map_ = CreateFileMappingA(h_file_, nullptr, PAGE_READONLY, 0, 0, nullptr);
        if (!h_map_) { cerrar(); return false; }
        data_ = static_cast<const char*>(MapViewOfFile(h_map_, FILE_MAP_READ, 0, 0, 0));
        if (!data_) { cerrar(); return false; }
        return true;
    }

    void cerrar() {
        if (data_) { UnmapViewOfFile(data_); data_ = nullptr; }
        if (h_map_) { CloseHandle(h_map_); h_map_ = nullptr; }
        if (h_file_ != INVALID_HANDLE_VALUE) { CloseHandle(h_file_); h_file_ = INVALID_HANDLE_VALUE; }
        size_ = 0;
    }

    [[nodiscard]] const char* data() const noexcept { return data_; }
    [[nodiscard]] size_t size() const noexcept { return size_; }
    [[nodiscard]] bool ok() const noexcept { return data_ != nullptr && size_ > 0; }
};

// =======================================================================
// ParserFastFloat: conversion numerica sin copias (std::from_chars, IEEE-754 redondeo al mas cercano)
// =======================================================================
class ParserFastFloat {
public:
    // YYYY-MM-DD estricto (digitos y fecha que existe) a YYYYMMDD.
    static inline bool parse_fecha(const char*& p, const char* end, uint32_t& out_fecha) noexcept {
        if (end - p < 10 || p[4] != '-' || p[7] != '-') return false;
        for (int i : {0, 1, 2, 3, 5, 6, 8, 9})
            if (p[i] < '0' || p[i] > '9') return false;
        const uint32_t y = static_cast<uint32_t>((p[0] - '0') * 1000 + (p[1] - '0') * 100 + (p[2] - '0') * 10 + (p[3] - '0'));
        const uint32_t m = static_cast<uint32_t>((p[5] - '0') * 10 + (p[6] - '0'));
        const uint32_t d = static_cast<uint32_t>((p[8] - '0') * 10 + (p[9] - '0'));
        const uint32_t f = y * 10000 + m * 100 + d;
        if (!fecha_valida(f)) return false;
        out_fecha = f;
        p += 10;
        return true;
    }

    // Un campo [a, b): vacio -> NaN; numero valido ocupando todo el campo -> su valor; cualquier otra cosa -> NaN y malo = true.
    static inline double parse_campo(const char* a, const char* b, bool& malo) noexcept {
        while (a < b && (*a == ' ' || *a == '\t')) ++a;
        while (b > a && (b[-1] == ' ' || b[-1] == '\t' || b[-1] == '\r')) --b;
        if (a == b) return NAN_DATASET;
        double v = 0.0;
        auto r = std::from_chars(a, b, v);
        if (r.ec == std::errc{} && r.ptr == b) return v;
        malo = true;
        return NAN_DATASET;
    }

    // API anterior: lee un doble y avanza p hasta la coma o el fin de linea.
    static inline double parse_double(const char*& p, const char* end) noexcept {
        const char* a = p;
        while (p < end && *p != ',' && *p != '\n') ++p;
        bool malo = false;
        return parse_campo(a, p, malo);
    }
};

struct InformeLectura {
    size_t filas = 0;            // filas de datos leidas
    size_t lineas_malas = 0;     // fecha invalida o ilegible (se omiten)
    size_t campos_malos = 0;     // campos con texto no numerico (quedan NaN)
    size_t filas_vacias = 0;     // filas cuyo cierre no es valido (es_imputado = 1)
    bool cabecera_ok = false;
    bool capacidad_agotada = false;
};

// =======================================================================
// LectorCSVDataset: fila cruda de la fuente hacia BarraCausal
// =======================================================================
class LectorCSVDataset {
public:
    static size_t leer_csv_macro(const std::string& ruta, BarraCausal* out_barras, size_t max_barras, Venue venue, const char* simbolo,
                                 InformeLectura* informe = nullptr) {
        InformeLectura inf;
        const size_t n = leer(ruta, out_barras, max_barras, venue, simbolo, inf);
        if (informe) *informe = inf;
        return n;
    }

private:
    static size_t leer(const std::string& ruta, BarraCausal* out, size_t max_barras, Venue, const char* simbolo, InformeLectura& inf) {
        if (!out || max_barras == 0) return 0;
        LectorMemoryMappedWin32 mm;
        if (!mm.abrir(ruta)) return 0;
        const char* p = mm.data();
        const char* end = p + mm.size();
        if (end - p >= 3 && static_cast<unsigned char>(p[0]) == 0xEF && static_cast<unsigned char>(p[1]) == 0xBB &&
            static_cast<unsigned char>(p[2]) == 0xBF)
            p += 3;   // BOM de UTF-8

        // Cabecera: indices de columna por nombre.
        int i_f = -1, i_o = -1, i_h = -1, i_l = -1, i_c = -1, i_v = -1, ncol = 0;
        {
            const char* a = p;
            while (p < end && *p != '\n') ++p;
            const char* fin = p;
            if (p < end) ++p;
            const char* s = a;
            for (int k = 0;; ++k) {
                const char* e = s;
                while (e < fin && *e != ',') ++e;
                std::string nom(s, e);
                while (!nom.empty() && (nom.back() == '\r' || nom.back() == ' ')) nom.pop_back();
                while (!nom.empty() && nom.front() == ' ') nom.erase(nom.begin());
                for (auto& ch : nom) ch = static_cast<char>((ch >= 'A' && ch <= 'Z') ? ch + 32 : ch);
                if (nom == "fecha" || nom == "date") i_f = k;
                else if (nom == "open") i_o = k;
                else if (nom == "high") i_h = k;
                else if (nom == "low") i_l = k;
                else if (nom == "close" || nom == "precio") i_c = k;
                else if (nom == "volume" || nom == "volumen") i_v = k;
                ncol = k + 1;
                if (e >= fin) break;
                s = e + 1;
            }
        }
        inf.cabecera_ok = (i_f >= 0 && i_c >= 0);
        if (!inf.cabecera_ok) return 0;

        size_t count = 0;
        const char* campos[64];
        const char* fines[64];
        while (p < end) {
            const char* a = p;
            while (p < end && *p != '\n') ++p;
            const char* fin = p;
            if (p < end) ++p;
            if (fin > a && fin[-1] == '\r') --fin;
            if (fin == a) continue;   // linea vacia
            int k = 0;
            for (const char* s = a; k < 64;) {
                const char* e = s;
                while (e < fin && *e != ',') ++e;
                campos[k] = s; fines[k] = e; ++k;
                if (e >= fin) break;
                s = e + 1;
            }
            if (k <= i_f) { ++inf.lineas_malas; continue; }
            uint32_t fecha = 0;
            const char* pf = campos[i_f];
            if (!(fines[i_f] - pf == 10 && ParserFastFloat::parse_fecha(pf, fines[i_f], fecha))) { ++inf.lineas_malas; continue; }
            if (count >= max_barras) { inf.capacidad_agotada = true; break; }
            bool malo = false;
            auto campo = [&](int idx) { return (idx >= 0 && idx < k) ? ParserFastFloat::parse_campo(campos[idx], fines[idx], malo) : NAN_DATASET; };
            BarraCausal& b = out[count];
            b = BarraCausal{};
            b.fecha_sesion = fecha;
            b.open = campo(i_o);
            b.high = campo(i_h);
            b.low = campo(i_l);
            b.close = campo(i_c);
            b.volume = campo(i_v);
            if (malo) ++inf.campos_malos;
            const bool c_ok = std::isfinite(b.close) && b.close > 0.0;
            b.es_imputado = !c_ok;
            b.hay_dato = c_ok;
            b.es_sesion_con_dato = c_ok;   // provisional: la rejilla lo cruza con el calendario
            if (!c_ok) { ++inf.filas_vacias; b.flags_control |= FLAG_FILA_VACIA; }
            b.valor_nominal = b.close * b.volume;
            b.fraccion_sesion = 1.0;
            b.estado_supervivencia = EstadoSupervivencia::COTIZANDO;
            b.estampa_disponible_ms = static_cast<uint64_t>(std::max<int64_t>(0, utc_ms_medianoche(fecha) + 86400000LL));
            huella_barra(simbolo, fecha, b.open, b.high, b.low, b.close, b.volume, b.hash_sha256);
            ++count;
            (void)ncol;
        }
        inf.filas = count;
        return count;
    }
};

// =======================================================================
// Rejilla diaria causal
// =======================================================================
enum class RejillaDiaria : uint8_t { SOLO_SESIONES = 0, SESIONES_Y_FERIADOS = 1, CALENDARIO_COMPLETO = 2 };

enum class ErrorRejilla : uint8_t {
    OK = 0,
    VACIA = 1,
    NO_MONOTONA = 2,          // ErrorMonotoniaTiempo de Python: la fuente trae fechas repetidas o desordenadas
    CENSURA_IZQUIERDA = 3,    // ViolacionCausalidadError: la primera fila no tiene precio valido (prohibido el bfill)
    SIN_SESION_PREVIA = 4,    // ventana_barra no encuentra sesion en las dos semanas previas (ValueError de Python)
    CAPACIDAD = 5
};

[[nodiscard]] constexpr const char* nombre_error(ErrorRejilla e) noexcept {
    switch (e) {
        case ErrorRejilla::OK: return "OK";
        case ErrorRejilla::VACIA: return "VACIA";
        case ErrorRejilla::NO_MONOTONA: return "NO_MONOTONA";
        case ErrorRejilla::CENSURA_IZQUIERDA: return "CENSURA_IZQUIERDA";
        case ErrorRejilla::SIN_SESION_PREVIA: return "SIN_SESION_PREVIA";
        default: return "CAPACIDAD";
    }
}

struct OpcionesRejilla {
    const char* fuente = nullptr;                        // p. ej. "YAHOO_VIA_MACRO_D1" (desambigua US02Y)
    CorteDiarioCripto corte = CorteDiarioCripto::UTC_00;
    RejillaDiaria rejilla = RejillaDiaria::SESIONES_Y_FERIADOS;
    uint32_t fin = 0;                                    // ultima sesion de la rejilla; 0 = la ultima de la fuente
};

struct ResultadoRejilla {
    size_t n = 0;
    ErrorRejilla error = ErrorRejilla::OK;
    uint32_t fecha_error = 0;
    size_t insertadas = 0;
    size_t pre_listado = 0, deslistado = 0;
};

// Construye la rejilla desde las filas crudas (en orden de fecha). `out` necesita capacidad para las insertadas.
[[nodiscard]] inline ResultadoRejilla construir_rejilla_diaria(const char* simbolo, const BarraCausal* crudas, size_t n_crudas,
                                                               BarraCausal* out, size_t max_out, const OpcionesRejilla& op = {}) {
    ResultadoRejilla r;
    if (!crudas || n_crudas == 0 || !out) { r.error = ErrorRejilla::VACIA; return r; }
    for (size_t i = 1; i < n_crudas; ++i)
        if (crudas[i].fecha_sesion <= crudas[i - 1].fecha_sesion) { r.error = ErrorRejilla::NO_MONOTONA; r.fecha_error = crudas[i].fecha_sesion; return r; }

    const std::string sim = simbolo_limpio(simbolo);
    const Calendario cal = calendario_de(sim.c_str(), op.fuente);

    // a. rejilla: fusion de las crudas con los dias operables que faltan en [primera, hi]
    const int64_t z0 = dias_desde_yyyymmdd(crudas[0].fecha_sesion);
    int64_t z_hi = dias_desde_yyyymmdd(crudas[n_crudas - 1].fecha_sesion);
    if (op.fin != 0) z_hi = std::max(z_hi, dias_desde_yyyymmdd(op.fin));
    size_t k = 0, n = 0;
    auto poner = [&](const BarraCausal& b) -> bool {
        if (n >= max_out) { r.error = ErrorRejilla::CAPACIDAD; return false; }
        out[n++] = b;
        return true;
    };
    for (int64_t z = z0; z <= z_hi; ++z) {
        const uint32_t f = yyyymmdd_desde_dias(z);
        if (k < n_crudas && crudas[k].fecha_sesion == f) {
            if (!poner(crudas[k])) return r;
            ++k;
            continue;
        }
        if (op.rejilla == RejillaDiaria::SOLO_SESIONES) continue;
        if (op.rejilla == RejillaDiaria::SESIONES_Y_FERIADOS && !cal.es_dia_operable(z)) continue;
        BarraCausal b{};
        b.fecha_sesion = f;
        b.open = b.high = b.low = b.close = NAN_DATASET;
        b.volume = 0.0;   // un dia sin sesion no negocio: volumen 0, no desconocido
        b.es_imputado = true;
        b.flags_control = FLAG_FILA_INSERTADA;
        b.estado_supervivencia = EstadoSupervivencia::COTIZANDO;
        b.estampa_disponible_ms = static_cast<uint64_t>(std::max<int64_t>(0, utc_ms_medianoche(f) + 86400000LL));
        if (!poner(b)) return r;
        ++r.insertadas;
    }
    r.n = n;

    // b. anotacion por calendario
    for (size_t i = 0; i < n; ++i) {
        BarraCausal& b = out[i];
        const ClaseDia c = cal.clasificar_fecha(b.fecha_sesion);
        const bool original = (b.flags_control & FLAG_FILA_INSERTADA) == 0;
        b.tipo_dia = static_cast<uint8_t>(c.tipo);
        b.fraccion_sesion = c.fraccion;
        if (std::isnan(c.fraccion)) b.flags_control |= FLAG_FRACCION_NAN;
        b.hay_dato = original && !b.es_imputado;
        b.es_sesion_con_dato = b.hay_dato && es_tipo_sesion(c.tipo);
        b.es_imputado = b.es_imputado || !original;
        if (c.tipo == TipoDia::FIN_DE_SEMANA) b.tipo_hueco = TipoHueco::FIN_DE_SEMANA;
        else if (!es_tipo_sesion(c.tipo)) b.tipo_hueco = TipoHueco::CIERRE_PROGRAMADO;
        else b.tipo_hueco = b.hay_dato ? TipoHueco::NINGUNO : TipoHueco::CORTE_FEED;
    }

    // c. supervivencia (primera pasada, dentro de la rejilla, como en Python)
    CatalogoSupervivenciaCxx::aplicar_mascara(sim.c_str(), out, n);

    // d. dias y sesiones acumulados; cierre real por perfil
    int64_t z_prev = 0;
    bool hay_prev = false;
    for (size_t i = 0; i < n; ++i) {
        BarraCausal& b = out[i];
        const int64_t z = dias_desde_yyyymmdd(b.fecha_sesion);
        if (b.es_sesion_con_dato) {
            if (hay_prev) {
                b.dias_acumulados = static_cast<uint16_t>(std::min<int64_t>(65535, std::max<int64_t>(1, z - z_prev)));
                b.sesiones_acumuladas = static_cast<uint16_t>(std::min<int64_t>(65535, std::max<int64_t>(1, cal.contar_sesiones(z_prev + 1, z))));
            } else {
                b.dias_acumulados = 1;
                b.sesiones_acumuladas = 1;
            }
            z_prev = z;
            hay_prev = true;
        } else {
            b.dias_acumulados = 0;
            b.sesiones_acumuladas = 0;
        }
        const VentanaBarra v = ventana_barra(sim.c_str(), b.fecha_sesion, op.fuente, op.corte);
        if (!v.ok) { r.error = ErrorRejilla::SIN_SESION_PREVIA; r.fecha_error = b.fecha_sesion; return r; }
        b.cierre_lo_utc_ms = v.cierre_lo_ms;
        b.cierre_hi_utc_ms = v.cierre_hi_ms;
        b.cierre_real_utc_ms = static_cast<uint64_t>(std::max<int64_t>(0, v.cierre_hi_ms));
        b.apertura_real_utc_ms = v.tiene_inicio ? static_cast<uint64_t>(std::max<int64_t>(0, v.inicio_lo_ms)) : 0;
    }

    // e. blindaje: precio invalido -> NaN; censura izquierda prohibida; arrastre causal
    auto invalido = [](double x) { return !std::isfinite(x) || x <= 0.0; };
    std::vector<uint8_t> inval(n, 0);   // fuera del camino caliente: la rejilla se construye una vez
    for (size_t i = 0; i < n; ++i) {
        BarraCausal& b = out[i];
        for (double* x : {&b.close, &b.open, &b.high, &b.low})
            if (invalido(*x)) { *x = NAN_DATASET; inval[i] = 1; }
    }
    // DIVERGE DE PYTHON: la censura izquierda se mira en la primera fila DENTRO del ciclo de vida. Python la mira en la fila 0,
    // que la supervivencia de la rejilla ya puso en NaN si es anterior al listado: por eso el pipeline de Python no puede cargar
    // DXY (datos desde 1971, listado 1973-03-01) ni VIX (1990 contra 1993-01-02) de cachas/macro_d1. Las filas anteriores al
    // listado quedan PRE_LISTADO, sin precio, y no alimentan el arrastre.
    size_t i0 = 0;
    while (i0 < n && out[i0].estado_supervivencia == EstadoSupervivencia::PRE_LISTADO) ++i0;
    if (i0 == n) { r.error = ErrorRejilla::VACIA; return r; }
    if (inval[i0]) { r.error = ErrorRejilla::CENSURA_IZQUIERDA; r.fecha_error = out[i0].fecha_sesion; return r; }
    for (size_t i = i0 + 1; i < n; ++i) {
        BarraCausal& b = out[i];
        const BarraCausal& a = out[i - 1];
        if (std::isnan(b.close)) b.close = a.close;
        if (std::isnan(b.open)) b.open = a.open;
        if (std::isnan(b.high)) b.high = a.high;
        if (std::isnan(b.low)) b.low = a.low;
    }
    for (size_t i = 0; i < n; ++i) {
        out[i].es_imputado = out[i].es_imputado || inval[i] != 0;
        out[i].valor_nominal = out[i].close * out[i].volume;
    }

    // f. supervivencia otra vez: corta el arrastre despues del desliste
    CatalogoSupervivenciaCxx::aplicar_mascara(sim.c_str(), out, n);
    for (size_t i = 0; i < n; ++i) {
        if (out[i].estado_supervivencia == EstadoSupervivencia::PRE_LISTADO) ++r.pre_listado;
        if (out[i].estado_supervivencia == EstadoSupervivencia::DESLISTADO) ++r.deslistado;
        out[i].valor_nominal = out[i].close * out[i].volume;
    }
    return r;
}

// Conveniencia fuera del camino caliente: lee el CSV y construye la rejilla.
[[nodiscard]] inline ResultadoRejilla cargar_rejilla_diaria(const std::string& ruta, const char* simbolo, std::vector<BarraCausal>& out,
                                                            const OpcionesRejilla& op = {}, InformeLectura* informe = nullptr) {
    std::vector<BarraCausal> crudas(200000);
    InformeLectura inf;
    const size_t n = LectorCSVDataset::leer_csv_macro(ruta, crudas.data(), crudas.size(), Venue::YAHOO, simbolo, &inf);
    if (informe) *informe = inf;
    ResultadoRejilla r;
    if (n == 0) { r.error = ErrorRejilla::VACIA; out.clear(); return r; }
    const int64_t z0 = dias_desde_yyyymmdd(crudas[0].fecha_sesion);
    int64_t z1 = dias_desde_yyyymmdd(crudas[n - 1].fecha_sesion);
    if (op.fin != 0) z1 = std::max(z1, dias_desde_yyyymmdd(op.fin));
    out.resize(static_cast<size_t>(z1 - z0 + 1) + n);
    r = construir_rejilla_diaria(simbolo, crudas.data(), n, out.data(), out.size(), op);
    out.resize(r.error == ErrorRejilla::OK ? r.n : 0);
    return r;
}

// =======================================================================
// SerializadorBinarioDataset: volcado plano (version 2)
// =======================================================================
class SerializadorBinarioDataset {
public:
    static constexpr uint32_t MAGICO = 0x41524B44;   // "ARKD"
    static constexpr uint32_t VERSION = 2;

    static bool guardar_binario(const std::string& ruta, const BarraCausal* barras, size_t n) {
        if (!barras || n == 0) return false;
        FILE* f = nullptr;
        fopen_s(&f, ruta.c_str(), "wb");
        if (!f) return false;
        const uint32_t magico = MAGICO, version = VERSION;
        const uint64_t total = n, tam = sizeof(BarraCausal);
        bool ok = std::fwrite(&magico, sizeof(magico), 1, f) == 1 && std::fwrite(&version, sizeof(version), 1, f) == 1 &&
                  std::fwrite(&total, sizeof(total), 1, f) == 1 && std::fwrite(&tam, sizeof(tam), 1, f) == 1 &&
                  std::fwrite(barras, sizeof(BarraCausal), n, f) == n;
        ok = (std::fclose(f) == 0) && ok;
        return ok;
    }

    static size_t cargar_binario(const std::string& ruta, BarraCausal* out_barras, size_t max_barras) {
        if (!out_barras || max_barras == 0) return 0;
        FILE* f = nullptr;
        fopen_s(&f, ruta.c_str(), "rb");
        if (!f) return 0;
        uint32_t magico = 0, version = 0;
        uint64_t total = 0, tam = 0;
        if (std::fread(&magico, sizeof(magico), 1, f) != 1 || magico != MAGICO || std::fread(&version, sizeof(version), 1, f) != 1 ||
            version != VERSION || std::fread(&total, sizeof(total), 1, f) != 1 || std::fread(&tam, sizeof(tam), 1, f) != 1 ||
            tam != sizeof(BarraCausal)) {
            std::fclose(f);
            return 0;
        }
        const size_t leer = std::min(max_barras, static_cast<size_t>(total));
        const size_t leidos = std::fread(out_barras, sizeof(BarraCausal), leer, f);
        std::fclose(f);
        return leidos;
    }
};

}  // namespace ark::dataset
