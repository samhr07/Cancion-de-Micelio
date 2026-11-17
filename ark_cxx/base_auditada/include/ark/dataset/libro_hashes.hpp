// libro_hashes.hpp -- HUELLA SHA-256 DE CADA BARRA Y LIBRO APPEND-ONLY EN DISCO (ARK-DATA-01-CXX)
// Reactor ARK / Cancion del Micelio -- ARQUITECTO DE DATOS (Flash) y correccion del Jefe 2026-11-01
//
// Huella canonica: la MISMA cadena que libro_hashes.computar_hash_barra de Python,
//     SIMBOLO_EN_MAYUSCULAS|YYYY-MM-DD|open:.6f|high:.6f|low:.6f|close:.6f|volumen:.2f
// (open/high/low ausentes o no finitos -> el cierre; volumen ausente o no finito -> 0.0), asi que una barra da el mismo hash en
// los dos lenguajes (lo comprueba tests/contraste_python.py sobre las 140 000 barras de cachas/macro_d1). Antes era otra cadena
// (%u|%.8f), y la barra guardaba 32 caracteres hex en ASCII (media huella).
//
// Libro: con una ruta, cada barra nueva se ANADE a un fichero binario de registros fijos y se carga al abrir; la primera huella
// registrada de una clave (venue, simbolo, fecha) es la inmutable. Antes la ruta se guardaba y nunca se escribia: una
// revision retroactiva solo se detectaba dentro del mismo proceso. Sin ruta (constructor por omision) el libro vive en memoria.
//
// SHA-256 con Windows CNG (bcrypt.dll), con el proveedor abierto una sola vez y sin memoria dinamica por barra.
#pragma once

#include "tipos_dataset.hpp"
#include "zona_horaria.hpp"
#include <windows.h>
#include <bcrypt.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <unordered_map>

#pragma comment(lib, "bcrypt.lib")

namespace ark::dataset {

namespace detalle_hash {
// Proveedor SHA-256 de CNG abierto una vez (inicializacion de static local: segura entre hilos).
[[nodiscard]] inline BCRYPT_ALG_HANDLE proveedor_sha256() noexcept {
    static BCRYPT_ALG_HANDLE h = [] {
        BCRYPT_ALG_HANDLE a = nullptr;
        if (BCryptOpenAlgorithmProvider(&a, BCRYPT_SHA256_ALGORITHM, nullptr, 0) != 0) a = nullptr;
        return a;
    }();
    return h;
}
}  // namespace detalle_hash

// SHA-256 binario de `n` bytes. Devuelve false si CNG falla (la huella queda en ceros).
inline bool sha256_bin(const void* datos, size_t n, uint8_t out[32]) noexcept {
    std::memset(out, 0, 32);
    BCRYPT_ALG_HANDLE a = detalle_hash::proveedor_sha256();
    if (!a) return false;
    return BCryptHash(a, nullptr, 0, static_cast<PUCHAR>(const_cast<void*>(datos)), static_cast<ULONG>(n), out, 32) == 0;
}

inline void hex_de(const uint8_t h[32], char out[65]) noexcept {
    static const char* d = "0123456789abcdef";
    for (int i = 0; i < 32; ++i) { out[2 * i] = d[h[i] >> 4]; out[2 * i + 1] = d[h[i] & 15]; }
    out[64] = '\0';
}

// Cadena canonica de la barra (sin memoria dinamica). Devuelve su longitud, o 0 si no cabe o el cierre no es finito.
inline int cadena_canonica_barra(const char* simbolo, uint32_t fecha_sesion, double o, double h, double l, double c, double v,
                                 char* buf, size_t cap) noexcept {
    if (!std::isfinite(c) || !buf || cap == 0) return 0;
    char sim[64];
    size_t k = 0;
    // str(simbolo).upper().strip(): se quitan espacios de los extremos y se pasa a mayusculas (ASCII).
    const char* s = simbolo ? simbolo : "";
    while (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n') ++s;
    while (*s && k < sizeof(sim) - 1) { char ch = *s++; sim[k++] = (ch >= 'a' && ch <= 'z') ? static_cast<char>(ch - 32) : ch; }
    while (k > 0 && (sim[k - 1] == ' ' || sim[k - 1] == '\t' || sim[k - 1] == '\r' || sim[k - 1] == '\n')) --k;
    sim[k] = '\0';
    const double oo = std::isfinite(o) ? o : c, hh = std::isfinite(h) ? h : c, ll = std::isfinite(l) ? l : c;
    const double vv = std::isfinite(v) ? v : 0.0;
    const int len = std::snprintf(buf, cap, "%s|%04u-%02u-%02u|%.6f|%.6f|%.6f|%.6f|%.2f", sim, fecha_sesion / 10000,
                                  (fecha_sesion / 100) % 100, fecha_sesion % 100, oo, hh, ll, c, vv);
    return (len > 0 && static_cast<size_t>(len) < cap) ? len : 0;
}

// Huella binaria de la barra fisica. false si el cierre no es finito (la fila no tiene barra que firmar).
inline bool huella_barra(const char* simbolo, uint32_t fecha_sesion, double o, double h, double l, double c, double v,
                         uint8_t out[32]) noexcept {
    char buf[512];
    const int n = cadena_canonica_barra(simbolo, fecha_sesion, o, h, l, c, v, buf, sizeof(buf));
    if (n <= 0) { std::memset(out, 0, 32); return false; }
    return sha256_bin(buf, static_cast<size_t>(n), out);
}

// =======================================================================
// EntradaRegistroHash: registro binario fijo del libro en disco (64 bytes)
// =======================================================================
#pragma pack(push, 4)
struct EntradaRegistroHash {
    uint32_t fecha_sesion;     // YYYYMMDD
    char simbolo[20];          // terminado en '\0' (antes 12 bytes: "LUNAUSDT_OLD" no cabia con su terminador)
    uint8_t venue_id;
    uint8_t flags;
    uint16_t reservado;
    uint32_t reservado2;
    uint8_t hash_sha256[32];
};
#pragma pack(pop)
static_assert(sizeof(EntradaRegistroHash) == 64, "EntradaRegistroHash tiene que medir 64 bytes.");

class LibroHashesDataset {
public:
    enum class ResultadoAuditoriaHash : uint8_t {
        NUEVA_BARRA_REGISTRADA = 0,
        BARRA_IDENTICA_VALIDADA = 1,
        MUTACION_RETROACTIVA_CAZADA = 2,
        SIN_BARRA = 3   // el cierre no es finito: no hay barra que firmar
    };

    // Sin ruta: libro en memoria (no deja rastro en disco). Con ruta: carga el fichero y anade cada barra nueva al final.
    LibroHashesDataset() = default;
    explicit LibroHashesDataset(std::string ruta_ledger) : ruta_(std::move(ruta_ledger)) { cargar(); }

    LibroHashesDataset(const LibroHashesDataset&) = delete;
    LibroHashesDataset& operator=(const LibroHashesDataset&) = delete;
    ~LibroHashesDataset() { if (f_) std::fclose(f_); }

    // Huella en hex de la barra (API anterior; misma cadena que Python).
    static std::string calcular_hash_barra(const char* simbolo, uint32_t fecha_sesion, double o, double h, double l, double c, double v) {
        uint8_t b[32];
        if (!huella_barra(simbolo, fecha_sesion, o, h, l, c, v, b)) return "";
        char x[65];
        hex_de(b, x);
        return std::string(x, 64);
    }

    ResultadoAuditoriaHash auditar_y_registrar(const char* simbolo, Venue venue, uint32_t fecha_sesion, double o, double h, double l,
                                               double c, double v, std::string* out_hash = nullptr) {
        uint8_t b[32];
        if (!huella_barra(simbolo, fecha_sesion, o, h, l, c, v, b)) {
            if (out_hash) out_hash->clear();
            return ResultadoAuditoriaHash::SIN_BARRA;
        }
        if (out_hash) { char x[65]; hex_de(b, x); out_hash->assign(x, 64); }
        const std::string k = clave(venue, simbolo, fecha_sesion);
        auto it = idx_.find(k);
        if (it != idx_.end())
            return std::memcmp(it->second.data(), b, 32) == 0 ? ResultadoAuditoriaHash::BARRA_IDENTICA_VALIDADA
                                                              : ResultadoAuditoriaHash::MUTACION_RETROACTIVA_CAZADA;
        std::array<uint8_t, 32> hh{};
        std::memcpy(hh.data(), b, 32);
        idx_.emplace(k, hh);
        if (!ruta_.empty()) anadir(simbolo, venue, fecha_sesion, b);
        return ResultadoAuditoriaHash::NUEVA_BARRA_REGISTRADA;
    }

    [[nodiscard]] size_t total_barras_registradas() const noexcept { return idx_.size(); }
    [[nodiscard]] size_t registros_ilegibles() const noexcept { return ilegibles_; }
    [[nodiscard]] bool persistente() const noexcept { return !ruta_.empty(); }
    [[nodiscard]] bool error_disco() const noexcept { return error_disco_; }

private:
    std::string ruta_;
    std::unordered_map<std::string, std::array<uint8_t, 32>> idx_;
    FILE* f_ = nullptr;
    size_t ilegibles_ = 0;
    bool error_disco_ = false;

    static std::string clave(Venue v, const char* s, uint32_t f) {
        return std::to_string(static_cast<int>(v)) + ":" + simbolo_mayus(s) + ":" + std::to_string(f);
    }
    static std::string simbolo_mayus(const char* s) {
        std::string x = s ? s : "";
        for (auto& ch : x) ch = static_cast<char>((ch >= 'a' && ch <= 'z') ? ch - 32 : ch);
        return x;
    }

    void cargar() {
        FILE* r = nullptr;
        if (fopen_s(&r, ruta_.c_str(), "rb") == 0 && r) {
            EntradaRegistroHash e{};
            size_t n;
            while ((n = std::fread(&e, 1, sizeof(e), r)) == sizeof(e)) {
                e.simbolo[sizeof(e.simbolo) - 1] = '\0';
                if (!fecha_valida(e.fecha_sesion)) { ++ilegibles_; continue; }
                std::array<uint8_t, 32> hh{};
                std::memcpy(hh.data(), e.hash_sha256, 32);
                idx_.emplace(clave(static_cast<Venue>(e.venue_id), e.simbolo, e.fecha_sesion), hh);   // la primera gana
            }
            if (n != 0) ++ilegibles_;   // cola truncada (escritura interrumpida): no se reutiliza
            std::fclose(r);
        }
    }

    void anadir(const char* s, Venue v, uint32_t f, const uint8_t h[32]) {
        if (!f_ && fopen_s(&f_, ruta_.c_str(), "ab") != 0) { f_ = nullptr; error_disco_ = true; return; }
        EntradaRegistroHash e{};
        e.fecha_sesion = f;
        const std::string m = simbolo_mayus(s);
        std::memcpy(e.simbolo, m.c_str(), std::min(m.size(), sizeof(e.simbolo) - 1));
        e.venue_id = static_cast<uint8_t>(v);
        std::memcpy(e.hash_sha256, h, 32);
        if (std::fwrite(&e, sizeof(e), 1, f_) != 1 || std::fflush(f_) != 0) error_disco_ = true;
    }
};

} // namespace ark::dataset
