// formato_npz.hpp -- LECTOR DE .npz / .npy SIN DEPENDENCIAS (ARK-DATA-02-CXX, ingesta intradia)
// Reactor ARK / Cancion del Micelio -- Jefe de Arquitectos de Codigo, 2026-11-02
//
// Lo que hace numpy.load sobre las caches que leen los adaptadores de ARQUITECTO DE DATOS (cachas/cuadrivector_oos,
// vrp_dvol, bitfinex_1s, bbg_1min, calibracion_l3), sin zlib:
//   - ZIP: directorio central (con ZIP64: numpy escribe cada miembro con force_zip64), cabecera local, miembros STORED (0) y
//     DEFLATE (8); CRC-32 comprobado en cada miembro (un miembro corrupto es un error, nunca datos a medias).
//   - DEFLATE (RFC 1951): bloques almacenados, Huffman fijo y dinamico, con tablas de 15 bits.
//   - NPY (versiones 1, 2 y 3): descr, fortran_order y shape; tipos <f4 <f8 <i4 <i8 |b1 |u1 y <U<n> (UTF-32LE -> UTF-8).
// Conversiones: a_double() es astype(np.float64) (float -> double exacto), a_int64() es astype(np.int64) desde enteros,
// a_bool() es astype(bool).
#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace ark::dataset::formatos {

struct ErrorFormato : std::runtime_error {
    using std::runtime_error::runtime_error;
};

inline std::vector<uint8_t> leer_fichero(const std::string& ruta) {
    std::ifstream f(ruta, std::ios::binary | std::ios::ate);
    if (!f) throw ErrorFormato("no se puede abrir " + ruta);
    const std::streamsize n = f.tellg();
    f.seekg(0);
    std::vector<uint8_t> b(static_cast<size_t>(n));
    if (n > 0 && !f.read(reinterpret_cast<char*>(b.data()), n)) throw ErrorFormato("lectura incompleta de " + ruta);
    return b;
}

// ---------------------------------------------------------------------------------------------------------------- CRC-32
inline uint32_t crc32(const uint8_t* p, size_t n) {
    static const std::array<uint32_t, 256> T = [] {
        std::array<uint32_t, 256> t{};
        for (uint32_t i = 0; i < 256; ++i) {
            uint32_t c = i;
            for (int k = 0; k < 8; ++k) c = (c & 1u) ? (0xEDB88320u ^ (c >> 1)) : (c >> 1);
            t[i] = c;
        }
        return t;
    }();
    uint32_t c = 0xFFFFFFFFu;
    for (size_t i = 0; i < n; ++i) c = T[(c ^ p[i]) & 0xFFu] ^ (c >> 8);
    return c ^ 0xFFFFFFFFu;
}

// --------------------------------------------------------------------------------------------------------------- DEFLATE
namespace detalle_inflate {

constexpr int MAXBITS = 15;
constexpr uint32_t INVALIDO = 0xFFFFFFFFu;

// Tabla completa de 2^15 entradas: entrada = (simbolo << 4) | longitud. Codigos LSB-primero (indice = codigo invertido).
struct Tabla {
    std::vector<uint32_t> t = std::vector<uint32_t>(1u << MAXBITS, INVALIDO);
    void construir(const uint8_t* longitudes, int n) {
        std::fill(t.begin(), t.end(), INVALIDO);
        int cuenta[MAXBITS + 1] = {0};
        for (int i = 0; i < n; ++i) ++cuenta[longitudes[i]];
        cuenta[0] = 0;
        int siguiente[MAXBITS + 2] = {0};
        int codigo = 0;
        for (int b = 1; b <= MAXBITS; ++b) {
            codigo = (codigo + cuenta[b - 1]) << 1;
            siguiente[b] = codigo;
        }
        for (int s = 0; s < n; ++s) {
            const int L = longitudes[s];
            if (L == 0) continue;
            const uint32_t c = static_cast<uint32_t>(siguiente[L]++);
            uint32_t r = 0;
            for (int k = 0; k < L; ++k) r |= ((c >> k) & 1u) << (L - 1 - k);
            if (r >= (1u << L)) throw ErrorFormato("deflate: codigo Huffman sobresuscrito");
            for (uint32_t f = r; f < (1u << MAXBITS); f += (1u << L)) t[f] = (static_cast<uint32_t>(s) << 4) | static_cast<uint32_t>(L);
        }
    }
};

struct Lector {
    const uint8_t* p;
    size_t n, pos = 0;
    uint64_t buf = 0;
    int bits = 0;
    Lector(const uint8_t* d, size_t m) : p(d), n(m) {}
    void llenar() {
        while (bits <= 56) {
            const uint64_t b = (pos < n) ? p[pos] : 0u;   // relleno con ceros: el fin real se valida con 'consumido'
            ++pos;
            buf |= b << bits;
            bits += 8;
        }
    }
    uint32_t ver(int k) { if (bits < k) llenar(); return static_cast<uint32_t>(buf & ((1ull << k) - 1)); }
    void saltar(int k) { buf >>= k; bits -= k; }
    uint32_t leer(int k) { if (k == 0) return 0; const uint32_t v = ver(k); saltar(k); return v; }
    void alinear_byte() { const int r = bits & 7; saltar(r); }
    size_t consumido() const { return pos - static_cast<size_t>(bits / 8); }
    int decodificar(const Tabla& T) {
        const uint32_t e = T.t[ver(MAXBITS)];
        if (e == INVALIDO) throw ErrorFormato("deflate: codigo Huffman invalido");
        saltar(static_cast<int>(e & 15u));
        return static_cast<int>(e >> 4);
    }
};

inline const uint16_t BASE_LON[29] = {3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
inline const uint8_t EXTRA_LON[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
inline const uint16_t BASE_DIS[30] = {1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
inline const uint8_t EXTRA_DIS[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};

}  // namespace detalle_inflate

// Descomprime un flujo DEFLATE crudo (sin cabecera zlib/gzip). esperado: tamano de salida conocido (0 si no).
inline std::vector<uint8_t> inflar(const uint8_t* datos, size_t n, size_t esperado = 0, size_t* consumido = nullptr) {
    using namespace detalle_inflate;
    std::vector<uint8_t> out;
    out.reserve(esperado ? esperado : n * 4);
    Lector L(datos, n);
    Tabla lit, dis;
    bool final = false;
    while (!final) {
        final = L.leer(1) != 0;
        const uint32_t tipo = L.leer(2);
        if (tipo == 0) {
            L.alinear_byte();
            const uint32_t len = L.leer(16), nlen = L.leer(16);
            if ((len ^ 0xFFFFu) != nlen) throw ErrorFormato("deflate: bloque almacenado con LEN/NLEN incoherentes");
            for (uint32_t i = 0; i < len; ++i) out.push_back(static_cast<uint8_t>(L.leer(8)));
            continue;
        }
        if (tipo == 3) throw ErrorFormato("deflate: tipo de bloque reservado");
        if (tipo == 1) {
            uint8_t ll[288], dl[30];
            for (int i = 0; i < 144; ++i) ll[i] = 8;
            for (int i = 144; i < 256; ++i) ll[i] = 9;
            for (int i = 256; i < 280; ++i) ll[i] = 7;
            for (int i = 280; i < 288; ++i) ll[i] = 8;
            for (int i = 0; i < 30; ++i) dl[i] = 5;
            lit.construir(ll, 288);
            dis.construir(dl, 30);
        } else {
            const int hlit = static_cast<int>(L.leer(5)) + 257, hdist = static_cast<int>(L.leer(5)) + 1, hclen = static_cast<int>(L.leer(4)) + 4;
            static const int ORDEN[19] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};
            uint8_t cl[19] = {0};
            for (int i = 0; i < hclen; ++i) cl[ORDEN[i]] = static_cast<uint8_t>(L.leer(3));
            Tabla tcl;
            tcl.construir(cl, 19);
            uint8_t lens[320] = {0};
            int i = 0;
            while (i < hlit + hdist) {
                const int s = L.decodificar(tcl);
                if (s < 16) {
                    lens[i++] = static_cast<uint8_t>(s);
                } else {
                    int rep = 0;
                    uint8_t v = 0;
                    if (s == 16) {
                        if (i == 0) throw ErrorFormato("deflate: repeticion sin longitud previa");
                        v = lens[i - 1];
                        rep = 3 + static_cast<int>(L.leer(2));
                    } else if (s == 17) {
                        rep = 3 + static_cast<int>(L.leer(3));
                    } else {
                        rep = 11 + static_cast<int>(L.leer(7));
                    }
                    if (i + rep > hlit + hdist) throw ErrorFormato("deflate: repeticion fuera de rango");
                    while (rep--) lens[i++] = v;
                }
            }
            if (lens[256] == 0) throw ErrorFormato("deflate: falta el codigo de fin de bloque");
            lit.construir(lens, hlit);
            dis.construir(lens + hlit, hdist);
        }
        for (;;) {
            const int s = L.decodificar(lit);
            if (s < 256) {
                out.push_back(static_cast<uint8_t>(s));
            } else if (s == 256) {
                break;
            } else {
                const int li = s - 257;
                if (li >= 29) throw ErrorFormato("deflate: simbolo de longitud invalido");
                const size_t len = BASE_LON[li] + L.leer(EXTRA_LON[li]);
                const int di = L.decodificar(dis);
                if (di >= 30) throw ErrorFormato("deflate: simbolo de distancia invalido");
                const size_t d = BASE_DIS[di] + L.leer(EXTRA_DIS[di]);
                if (d > out.size()) throw ErrorFormato("deflate: distancia antes del inicio");
                const size_t ini = out.size() - d;
                for (size_t k = 0; k < len; ++k) out.push_back(out[ini + k]);
            }
        }
    }
    if (L.consumido() > n) throw ErrorFormato("deflate: flujo truncado");
    if (consumido) *consumido = L.consumido();
    return out;
}

// ------------------------------------------------------------------------------------------------------------------- ZIP
inline uint16_t le16(const uint8_t* p) { return static_cast<uint16_t>(p[0] | (p[1] << 8)); }
inline uint32_t le32(const uint8_t* p) { return static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8) | (static_cast<uint32_t>(p[2]) << 16) | (static_cast<uint32_t>(p[3]) << 24); }
inline uint64_t le64(const uint8_t* p) { return static_cast<uint64_t>(le32(p)) | (static_cast<uint64_t>(le32(p + 4)) << 32); }

struct MiembroZip {
    std::string nombre;
    uint16_t metodo = 0;
    uint32_t crc = 0;
    uint64_t comprimido = 0, descomprimido = 0, offset_local = 0;
};

inline std::vector<MiembroZip> directorio_zip(const std::vector<uint8_t>& z) {
    const size_t n = z.size();
    if (n < 22) throw ErrorFormato("zip: fichero demasiado corto");
    size_t eocd = n;
    for (size_t i = n - 22 + 1; i-- > (n > 22 + 65535 ? n - 22 - 65535 : 0);) {
        if (le32(&z[i]) == 0x06054b50u) { eocd = i; break; }
    }
    if (eocd == n) throw ErrorFormato("zip: sin fin de directorio central");
    uint64_t entradas = le16(&z[eocd + 10]), tam_cd = le32(&z[eocd + 12]), off_cd = le32(&z[eocd + 16]);
    if (eocd >= 20 && le32(&z[eocd - 20]) == 0x07064b50u) {   // localizador ZIP64
        const uint64_t off64 = le64(&z[eocd - 20 + 8]);
        if (off64 + 56 > n || le32(&z[off64]) != 0x06064b50u) throw ErrorFormato("zip: fin de directorio ZIP64 invalido");
        entradas = le64(&z[off64 + 32]);
        tam_cd = le64(&z[off64 + 40]);
        off_cd = le64(&z[off64 + 48]);
    }
    if (off_cd + tam_cd > n) throw ErrorFormato("zip: directorio central fuera del fichero");
    std::vector<MiembroZip> v;
    size_t p = static_cast<size_t>(off_cd);
    for (uint64_t e = 0; e < entradas; ++e) {
        if (p + 46 > n || le32(&z[p]) != 0x02014b50u) throw ErrorFormato("zip: entrada del directorio central invalida");
        MiembroZip m;
        m.metodo = le16(&z[p + 10]);
        m.crc = le32(&z[p + 16]);
        m.comprimido = le32(&z[p + 20]);
        m.descomprimido = le32(&z[p + 24]);
        const uint16_t ln = le16(&z[p + 28]), lx = le16(&z[p + 30]), lc = le16(&z[p + 32]);
        m.offset_local = le32(&z[p + 42]);
        m.nombre.assign(reinterpret_cast<const char*>(&z[p + 46]), ln);
        size_t x = p + 46 + ln;
        const size_t xfin = x + lx;
        while (x + 4 <= xfin) {
            const uint16_t id = le16(&z[x]), tam = le16(&z[x + 2]);
            if (id == 0x0001) {
                size_t q = x + 4;
                if (m.descomprimido == 0xFFFFFFFFu) { m.descomprimido = le64(&z[q]); q += 8; }
                if (m.comprimido == 0xFFFFFFFFu) { m.comprimido = le64(&z[q]); q += 8; }
                if (m.offset_local == 0xFFFFFFFFu) { m.offset_local = le64(&z[q]); q += 8; }
            }
            x += 4u + tam;
        }
        v.push_back(m);
        p = xfin + lc;
    }
    return v;
}

inline std::vector<uint8_t> extraer_miembro(const std::vector<uint8_t>& z, const MiembroZip& m) {
    const size_t o = static_cast<size_t>(m.offset_local);
    if (o + 30 > z.size() || le32(&z[o]) != 0x04034b50u) throw ErrorFormato("zip: cabecera local invalida en " + m.nombre);
    const size_t ini = o + 30 + le16(&z[o + 26]) + le16(&z[o + 28]);
    if (ini + m.comprimido > z.size()) throw ErrorFormato("zip: miembro truncado " + m.nombre);
    std::vector<uint8_t> d;
    if (m.metodo == 0) {
        d.assign(z.begin() + static_cast<std::ptrdiff_t>(ini), z.begin() + static_cast<std::ptrdiff_t>(ini + m.comprimido));
    } else if (m.metodo == 8) {
        d = inflar(&z[ini], static_cast<size_t>(m.comprimido), static_cast<size_t>(m.descomprimido));
    } else {
        throw ErrorFormato("zip: metodo de compresion no soportado en " + m.nombre);
    }
    if (d.size() != m.descomprimido) throw ErrorFormato("zip: tamano descomprimido distinto en " + m.nombre);
    if (crc32(d.data(), d.size()) != m.crc) throw ErrorFormato("zip: CRC-32 distinto en " + m.nombre);
    return d;
}

// ------------------------------------------------------------------------------------------------------------------- NPY
struct ArregloNpy {
    std::string descr;            // p. ej. "<f8", "<U10"
    std::vector<size_t> forma;
    bool fortran = false;
    std::vector<uint8_t> datos;   // bytes crudos en el orden del fichero

    size_t elementos() const { size_t k = 1; for (size_t s : forma) k *= s; return k; }
    size_t tam_item() const {
        if (descr.size() < 3) throw ErrorFormato("npy: descr invalido " + descr);
        const size_t b = static_cast<size_t>(std::stoul(descr.substr(2)));
        return descr[1] == 'U' ? 4 * b : b;
    }
    char tipo() const { return descr.size() > 1 ? descr[1] : '?'; }

    // Indice lineal en orden C de la posicion lineal i (para leer arreglos 2D con fortran_order True).
    size_t pos_c(size_t i) const {
        if (!fortran || forma.size() < 2) return i;
        // i en orden C -> (fila, col) -> offset Fortran
        size_t resto = i, off = 0, paso = 1;
        std::vector<size_t> idx(forma.size());
        for (size_t d = forma.size(); d-- > 0;) { idx[d] = resto % forma[d]; resto /= forma[d]; }
        for (size_t d = 0; d < forma.size(); ++d) { off += idx[d] * paso; paso *= forma[d]; }
        return off;
    }

    std::vector<double> a_double() const {
        const size_t n = elementos(), w = tam_item();
        std::vector<double> v(n);
        for (size_t i = 0; i < n; ++i) {
            const uint8_t* p = &datos[pos_c(i) * w];
            if (descr == "<f8") { double x; std::memcpy(&x, p, 8); v[i] = x; }
            else if (descr == "<f4") { float x; std::memcpy(&x, p, 4); v[i] = static_cast<double>(x); }
            else if (descr == "<i8") { int64_t x; std::memcpy(&x, p, 8); v[i] = static_cast<double>(x); }
            else if (descr == "<i4") { int32_t x; std::memcpy(&x, p, 4); v[i] = static_cast<double>(x); }
            else if (descr == "|b1" || descr == "|u1") { v[i] = p[0] ? 1.0 : 0.0; }
            else throw ErrorFormato("npy: a_double no soporta " + descr);
        }
        return v;
    }
    std::vector<int64_t> a_int64() const {
        const size_t n = elementos(), w = tam_item();
        std::vector<int64_t> v(n);
        for (size_t i = 0; i < n; ++i) {
            const uint8_t* p = &datos[pos_c(i) * w];
            if (descr == "<i8") { int64_t x; std::memcpy(&x, p, 8); v[i] = x; }
            else if (descr == "<i4") { int32_t x; std::memcpy(&x, p, 4); v[i] = x; }
            else if (descr == "|b1" || descr == "|u1") { v[i] = p[0]; }
            else if (descr == "<f8") { double x; std::memcpy(&x, p, 8); v[i] = static_cast<int64_t>(x); }   // astype(int64): trunca
            else throw ErrorFormato("npy: a_int64 no soporta " + descr);
        }
        return v;
    }
    std::vector<uint8_t> a_bool() const {
        const size_t n = elementos(), w = tam_item();
        std::vector<uint8_t> v(n);
        for (size_t i = 0; i < n; ++i) {
            const uint8_t* p = &datos[pos_c(i) * w];
            bool b = false;
            for (size_t k = 0; k < w; ++k) b = b || p[k] != 0;   // astype(bool): distinto de cero (incluye -0.0 -> False solo en float)
            if (descr == "<f8") { double x; std::memcpy(&x, p, 8); b = x != 0.0; }
            if (descr == "<f4") { float x; std::memcpy(&x, p, 4); b = x != 0.0f; }
            v[i] = b ? 1 : 0;
        }
        return v;
    }
    std::vector<std::string> a_cadenas() const {
        if (tipo() != 'U') throw ErrorFormato("npy: a_cadenas exige <U");
        const size_t n = elementos(), w = tam_item(), nch = w / 4;
        std::vector<std::string> v(n);
        for (size_t i = 0; i < n; ++i) {
            const uint8_t* p = &datos[pos_c(i) * w];
            std::string s;
            for (size_t k = 0; k < nch; ++k) {
                const uint32_t c = le32(p + 4 * k);
                if (c == 0) break;   // numpy rellena con NUL
                if (c < 0x80) s += static_cast<char>(c);
                else if (c < 0x800) { s += static_cast<char>(0xC0 | (c >> 6)); s += static_cast<char>(0x80 | (c & 0x3F)); }
                else if (c < 0x10000) { s += static_cast<char>(0xE0 | (c >> 12)); s += static_cast<char>(0x80 | ((c >> 6) & 0x3F)); s += static_cast<char>(0x80 | (c & 0x3F)); }
                else { s += static_cast<char>(0xF0 | (c >> 18)); s += static_cast<char>(0x80 | ((c >> 12) & 0x3F)); s += static_cast<char>(0x80 | ((c >> 6) & 0x3F)); s += static_cast<char>(0x80 | (c & 0x3F)); }
            }
            v[i] = s;
        }
        return v;
    }
};

inline ArregloNpy parsear_npy(const std::vector<uint8_t>& b) {
    if (b.size() < 10 || b[0] != 0x93 || std::memcmp(&b[1], "NUMPY", 5) != 0) throw ErrorFormato("npy: magia invalida");
    const uint8_t mayor = b[6];
    size_t hlen, hini;
    if (mayor == 1) { hlen = le16(&b[8]); hini = 10; }
    else if (mayor == 2 || mayor == 3) { hlen = le32(&b[8]); hini = 12; }
    else throw ErrorFormato("npy: version no soportada");
    if (hini + hlen > b.size()) throw ErrorFormato("npy: cabecera truncada");
    const std::string h(reinterpret_cast<const char*>(&b[hini]), hlen);
    ArregloNpy a;
    auto valor = [&](const std::string& clave) -> std::string {
        const size_t k = h.find("'" + clave + "'");
        if (k == std::string::npos) throw ErrorFormato("npy: falta " + clave);
        size_t p = h.find(':', k) + 1;
        while (p < h.size() && h[p] == ' ') ++p;
        if (h[p] == '\'') { const size_t q = h.find('\'', p + 1); return h.substr(p + 1, q - p - 1); }
        if (h[p] == '(') { const size_t q = h.find(')', p); return h.substr(p, q - p + 1); }
        size_t q = p;
        while (q < h.size() && h[q] != ',' && h[q] != '}') ++q;
        return h.substr(p, q - p);
    };
    a.descr = valor("descr");
    a.fortran = valor("fortran_order").find("True") != std::string::npos;
    const std::string f = valor("shape");
    size_t p = 1;
    while (p < f.size()) {
        while (p < f.size() && (f[p] == ' ' || f[p] == ',')) ++p;
        if (p >= f.size() || f[p] == ')') break;
        size_t q = p;
        while (q < f.size() && f[q] >= '0' && f[q] <= '9') ++q;
        a.forma.push_back(static_cast<size_t>(std::stoull(f.substr(p, q - p))));
        p = q;
    }
    const size_t ini = hini + hlen, bytes = a.elementos() * a.tam_item();
    if (ini + bytes > b.size()) throw ErrorFormato("npy: datos truncados");
    a.datos.assign(b.begin() + static_cast<std::ptrdiff_t>(ini), b.begin() + static_cast<std::ptrdiff_t>(ini + bytes));
    return a;
}

// np.load(npz): nombre (sin .npy) -> arreglo.
inline std::map<std::string, ArregloNpy> cargar_npz(const std::string& ruta) {
    const std::vector<uint8_t> z = leer_fichero(ruta);
    std::map<std::string, ArregloNpy> m;
    for (const MiembroZip& mi : directorio_zip(z)) {
        std::string nom = mi.nombre;
        if (nom.size() > 4 && nom.compare(nom.size() - 4, 4, ".npy") == 0) nom.resize(nom.size() - 4);
        m.emplace(nom, parsear_npy(extraer_miembro(z, mi)));
    }
    return m;
}

// Flujo gzip (RFC 1952) con uno o varios miembros concatenados: las capturas *.jsonl.gz.
inline std::vector<uint8_t> desgzip(const std::vector<uint8_t>& g) {
    std::vector<uint8_t> out;
    size_t p = 0;
    while (p + 18 <= g.size()) {
        if (g[p] != 0x1f || g[p + 1] != 0x8b || g[p + 2] != 8) throw ErrorFormato("gzip: cabecera invalida");
        const uint8_t flg = g[p + 3];
        size_t q = p + 10;
        if (flg & 4) { q += 2u + le16(&g[q]); }
        if (flg & 8) { while (q < g.size() && g[q]) ++q; ++q; }
        if (flg & 16) { while (q < g.size() && g[q]) ++q; ++q; }
        if (flg & 2) q += 2;
        if (q > g.size()) throw ErrorFormato("gzip: cabecera truncada");
        size_t usado = 0;
        std::vector<uint8_t> d = inflar(&g[q], g.size() - q, 0, &usado);
        q += usado;
        if (q + 8 > g.size()) throw ErrorFormato("gzip: falta el epilogo");
        if (crc32(d.data(), d.size()) != le32(&g[q])) throw ErrorFormato("gzip: CRC-32 distinto");
        if (static_cast<uint32_t>(d.size()) != le32(&g[q + 4])) throw ErrorFormato("gzip: ISIZE distinto");
        out.insert(out.end(), d.begin(), d.end());
        p = q + 8;
    }
    return out;
}

}  // namespace ark::dataset::formatos
