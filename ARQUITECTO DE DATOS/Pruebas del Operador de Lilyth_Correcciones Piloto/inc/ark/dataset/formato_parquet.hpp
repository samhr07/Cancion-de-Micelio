// formato_parquet.hpp -- LECTOR DE PARQUET SIN DEPENDENCIAS (ARK-DATA-02-CXX, ingesta intradia)
// Reactor ARK / Cancion del Micelio -- Jefe de Arquitectos de Codigo, 2026-11-02
//
// Lo que hace pandas.read_parquet (pyarrow) sobre data/*.parquet, el insumo de adaptador_binance.cargar_binance_parquet:
//   - pie FileMetaData en Thrift compacto (esquema, grupos de filas, columnas) y PageHeader de cada pagina;
//   - codecs UNCOMPRESSED y SNAPPY (descompresor propio);
//   - paginas DATA_PAGE (v1), DATA_PAGE_V2 y DICTIONARY_PAGE; codificaciones PLAIN, PLAIN_DICTIONARY y RLE_DICTIONARY;
//   - niveles de definicion (columnas OPTIONAL de un nivel) por RLE/bit-packing: un nulo es NaN en DOUBLE/FLOAT y se marca
//     aparte en INT64/INT32 (pandas daria NaT/NaN);
//   - tipos fisicos DOUBLE, FLOAT, INT64 e INT32; el tipo logico TIMESTAMP se expone con su unidad (ns/us/ms).
// Columnas anidadas o repetidas no se soportan (error explicito): las fuentes del dataset son tablas planas.
// Metadatos clave/valor (p. ej. 'pandas' con index_columns) se exponen tal cual para replicar el indice de pandas.
#pragma once
#include "ark/dataset/formato_npz.hpp"
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace ark::dataset::formatos {

// ---------------------------------------------------------------------------------------------------------------- SNAPPY
inline std::vector<uint8_t> dessnappy(const uint8_t* p, size_t n) {
    size_t i = 0;
    uint64_t largo = 0;
    for (int s = 0; i < n; s += 7) {
        const uint8_t b = p[i++];
        largo |= static_cast<uint64_t>(b & 0x7F) << s;
        if (!(b & 0x80)) break;
        if (s > 56) throw ErrorFormato("snappy: varint de longitud invalido");
    }
    std::vector<uint8_t> out;
    out.reserve(static_cast<size_t>(largo));
    while (i < n) {
        const uint8_t tag = p[i++];
        const uint8_t t = tag & 3u;
        if (t == 0) {
            size_t len = (tag >> 2) + 1u;
            if (len > 60) {
                const size_t nb = len - 60;
                if (i + nb > n) throw ErrorFormato("snappy: literal truncado");
                len = 0;
                for (size_t k = 0; k < nb; ++k) len |= static_cast<size_t>(p[i + k]) << (8 * k);
                len += 1;
                i += nb;
            }
            if (i + len > n) throw ErrorFormato("snappy: literal truncado");
            out.insert(out.end(), p + i, p + i + len);
            i += len;
        } else {
            size_t len, off;
            if (t == 1) {
                if (i + 1 > n) throw ErrorFormato("snappy: copia truncada");
                len = ((tag >> 2) & 7u) + 4u;
                off = (static_cast<size_t>(tag >> 5) << 8) | p[i];
                i += 1;
            } else if (t == 2) {
                if (i + 2 > n) throw ErrorFormato("snappy: copia truncada");
                len = (tag >> 2) + 1u;
                off = le16(p + i);
                i += 2;
            } else {
                if (i + 4 > n) throw ErrorFormato("snappy: copia truncada");
                len = (tag >> 2) + 1u;
                off = le32(p + i);
                i += 4;
            }
            if (off == 0 || off > out.size()) throw ErrorFormato("snappy: desplazamiento invalido");
            const size_t ini = out.size() - off;
            for (size_t k = 0; k < len; ++k) out.push_back(out[ini + k]);
        }
    }
    if (out.size() != largo) throw ErrorFormato("snappy: longitud final distinta");
    return out;
}

// ------------------------------------------------------------------------------------------------------- THRIFT COMPACTO
struct LectorThrift {
    const uint8_t* p;
    size_t n, i = 0;
    LectorThrift(const uint8_t* d, size_t m, size_t ini = 0) : p(d), n(m), i(ini) {}
    uint8_t byte() { if (i >= n) throw ErrorFormato("thrift: fin inesperado"); return p[i++]; }
    uint64_t varint() {
        uint64_t v = 0;
        for (int s = 0;; s += 7) {
            const uint8_t b = byte();
            v |= static_cast<uint64_t>(b & 0x7F) << s;
            if (!(b & 0x80)) return v;
            if (s > 63) throw ErrorFormato("thrift: varint invalido");
        }
    }
    int64_t zigzag() { const uint64_t v = varint(); return static_cast<int64_t>(v >> 1) ^ -static_cast<int64_t>(v & 1); }
    std::string cadena() {
        const size_t L = static_cast<size_t>(varint());
        if (i + L > n) throw ErrorFormato("thrift: cadena truncada");
        std::string s(reinterpret_cast<const char*>(p + i), L);
        i += L;
        return s;
    }
    // Cabecera de campo: devuelve (id, tipo); tipo 0 = STOP.
    std::pair<int, int> campo(int& ultimo) {
        const uint8_t b = byte();
        if (b == 0) return {0, 0};
        const int tipo = b & 0x0F, delta = b >> 4;
        const int id = delta ? ultimo + delta : static_cast<int>(zigzag());
        ultimo = id;
        return {id, tipo};
    }
    std::pair<int, size_t> lista() {   // (tipo de elemento, tamano)
        const uint8_t b = byte();
        size_t tam = b >> 4;
        if (tam == 15) tam = static_cast<size_t>(varint());
        return {b & 0x0F, tam};
    }
    void saltar(int tipo) {
        switch (tipo) {
            case 1: case 2: return;                         // BOOL true/false (valor en el tipo)
            case 3: byte(); return;                          // BYTE
            case 4: case 5: case 6: varint(); return;        // I16 I32 I64
            case 7: i += 8; return;                          // DOUBLE
            case 8: { const size_t L = static_cast<size_t>(varint()); i += L; return; }   // BINARY
            case 9: case 10: { auto [t, k] = lista(); for (size_t e = 0; e < k; ++e) saltar_elem(t); return; }   // LIST SET
            case 11: {                                       // MAP
                const size_t k = static_cast<size_t>(varint());
                if (k == 0) return;
                const uint8_t kv = byte();
                for (size_t e = 0; e < k; ++e) { saltar_elem(kv >> 4); saltar_elem(kv & 0x0F); }
                return;
            }
            case 12: { int u = 0; for (;;) { auto [id, t] = campo(u); if (t == 0) return; saltar(t); } }   // STRUCT
            default: throw ErrorFormato("thrift: tipo desconocido");
        }
    }
    void saltar_elem(int t) { if (t == 1 || t == 2) { byte(); return; } saltar(t); }   // BOOL en listas ocupa un byte
};

// --------------------------------------------------------------------------------------------------------------- PARQUET
enum class TipoFisico : int { BOOLEAN = 0, INT32 = 1, INT64 = 2, INT96 = 3, FLOAT = 4, DOUBLE = 5, BYTE_ARRAY = 6, FIXED = 7 };

struct ColumnaParquet {
    std::string nombre;
    TipoFisico tipo = TipoFisico::DOUBLE;
    bool opcional = false;
    std::string unidad_timestamp;   // "ns", "us", "ms" o "" si no es TIMESTAMP
    bool timestamp_ajustado_utc = false;
    std::vector<double> dobles;     // DOUBLE/FLOAT (nulo -> NaN)
    std::vector<int64_t> enteros;   // INT64/INT32
    std::vector<uint8_t> nulo;      // 1 si el valor es nulo
};

struct TablaParquet {
    int64_t filas = 0;
    std::vector<ColumnaParquet> columnas;
    std::vector<std::pair<std::string, std::string>> metadatos;
    const ColumnaParquet* columna(const std::string& n) const {
        for (const auto& c : columnas) if (c.nombre == n) return &c;
        return nullptr;
    }
};

namespace detalle_parquet {

struct MetaColumna {
    int tipo = -1, codec = 0;
    int64_t num_valores = 0, data_off = -1, dict_off = -1, comprimido = 0;
    std::vector<std::string> ruta;
};

// RLE / bit-packing hibrido de Parquet: decodifica 'cuantos' valores de ancho 'bw' desde p[0..n).
inline std::vector<uint32_t> rle_hibrido(const uint8_t* p, size_t n, int bw, size_t cuantos) {
    std::vector<uint32_t> v;
    v.reserve(cuantos);
    size_t i = 0;
    const size_t bytes_valor = static_cast<size_t>((bw + 7) / 8);
    while (v.size() < cuantos && i < n) {
        uint64_t h = 0;
        for (int s = 0;; s += 7) {
            if (i >= n) throw ErrorFormato("parquet: RLE truncado");
            const uint8_t b = p[i++];
            h |= static_cast<uint64_t>(b & 0x7F) << s;
            if (!(b & 0x80)) break;
        }
        if (h & 1u) {   // grupos de 8 bit-packed
            const size_t grupos = static_cast<size_t>(h >> 1), total = grupos * 8;
            const size_t nbytes = grupos * static_cast<size_t>(bw);
            if (i + nbytes > n) throw ErrorFormato("parquet: bit-packing truncado");
            for (size_t k = 0; k < total && v.size() < cuantos; ++k) {
                uint32_t x = 0;
                for (int b = 0; b < bw; ++b) {
                    const size_t bit = k * static_cast<size_t>(bw) + static_cast<size_t>(b);
                    x |= static_cast<uint32_t>((p[i + bit / 8] >> (bit % 8)) & 1u) << b;
                }
                v.push_back(x);
            }
            i += nbytes;
        } else {        // racha repetida
            const size_t cnt = static_cast<size_t>(h >> 1);
            uint32_t x = 0;
            for (size_t k = 0; k < bytes_valor; ++k) x |= static_cast<uint32_t>(p[i + k]) << (8 * k);
            i += bytes_valor;
            for (size_t k = 0; k < cnt && v.size() < cuantos; ++k) v.push_back(x);
        }
    }
    if (v.size() < cuantos) throw ErrorFormato("parquet: faltan valores RLE");
    return v;
}

inline std::vector<uint8_t> descomprimir(int codec, const uint8_t* p, size_t n) {
    if (codec == 0) return std::vector<uint8_t>(p, p + n);
    if (codec == 1) return dessnappy(p, n);
    throw ErrorFormato("parquet: codec no soportado (" + std::to_string(codec) + ")");
}

}  // namespace detalle_parquet

inline TablaParquet leer_parquet(const std::string& ruta) {
    using namespace detalle_parquet;
    const std::vector<uint8_t> f = leer_fichero(ruta);
    const size_t n = f.size();
    if (n < 12 || std::memcmp(&f[0], "PAR1", 4) != 0 || std::memcmp(&f[n - 4], "PAR1", 4) != 0) throw ErrorFormato("parquet: magia invalida");
    const size_t lm = le32(&f[n - 8]);
    if (lm + 8 > n) throw ErrorFormato("parquet: pie invalido");
    LectorThrift T(f.data(), n - 8, n - 8 - lm);

    TablaParquet tabla;
    struct Esquema { std::string nombre; int tipo = -1, repeticion = 0, hijos = 0; std::string unidad; bool utc = false; };
    std::vector<Esquema> esquema;
    std::vector<std::vector<MetaColumna>> grupos;

    auto leer_unidad = [&](LectorThrift& L, Esquema& e) {   // LogicalType union -> TIMESTAMP (campo 8)
        int u = 0;
        for (;;) {
            auto [id, t] = L.campo(u);
            if (t == 0) return;
            if (id == 8 && t == 12) {   // TimestampType {1 isAdjustedToUTC bool, 2 unit TimeUnit}
                int u2 = 0;
                for (;;) {
                    auto [id2, t2] = L.campo(u2);
                    if (t2 == 0) break;
                    if (id2 == 1) e.utc = (t2 == 1);
                    else if (id2 == 2 && t2 == 12) {
                        int u3 = 0;
                        for (;;) {
                            auto [id3, t3] = L.campo(u3);
                            if (t3 == 0) break;
                            e.unidad = id3 == 1 ? "ms" : (id3 == 2 ? "us" : (id3 == 3 ? "ns" : "?"));
                            L.saltar(t3);
                        }
                    } else L.saltar(t2);
                }
            } else L.saltar(t);
        }
    };

    int ultimo = 0;
    for (;;) {
        auto [id, t] = T.campo(ultimo);
        if (t == 0) break;
        if (id == 2 && t == 9) {
            auto [te, k] = T.lista();
            (void)te;
            for (size_t e = 0; e < k; ++e) {
                Esquema s;
                int u = 0;
                for (;;) {
                    auto [id2, t2] = T.campo(u);
                    if (t2 == 0) break;
                    if (id2 == 1) s.tipo = static_cast<int>(T.zigzag());
                    else if (id2 == 3) s.repeticion = static_cast<int>(T.zigzag());
                    else if (id2 == 4) s.nombre = T.cadena();
                    else if (id2 == 5) s.hijos = static_cast<int>(T.zigzag());
                    else if (id2 == 6) { const int conv = static_cast<int>(T.zigzag()); if (conv == 9) s.unidad = "ms"; else if (conv == 10) s.unidad = "us"; }
                    else if (id2 == 10 && t2 == 12) leer_unidad(T, s);
                    else T.saltar(t2);
                }
                esquema.push_back(s);
            }
        } else if (id == 3 && (t == 6 || t == 5)) {
            tabla.filas = T.zigzag();
        } else if (id == 4 && t == 9) {
            auto [te, k] = T.lista();
            (void)te;
            for (size_t g = 0; g < k; ++g) {
                std::vector<MetaColumna> cols;
                int u = 0;
                for (;;) {
                    auto [id2, t2] = T.campo(u);
                    if (t2 == 0) break;
                    if (id2 == 1 && t2 == 9) {
                        auto [tc, kc] = T.lista();
                        (void)tc;
                        for (size_t c = 0; c < kc; ++c) {
                            MetaColumna mc;
                            int u3 = 0;
                            for (;;) {
                                auto [id3, t3] = T.campo(u3);
                                if (t3 == 0) break;
                                if (id3 == 3 && t3 == 12) {
                                    int u4 = 0;
                                    for (;;) {
                                        auto [id4, t4] = T.campo(u4);
                                        if (t4 == 0) break;
                                        if (id4 == 1) mc.tipo = static_cast<int>(T.zigzag());
                                        else if (id4 == 3 && t4 == 9) { auto [ts, ks] = T.lista(); (void)ts; for (size_t q = 0; q < ks; ++q) mc.ruta.push_back(T.cadena()); }
                                        else if (id4 == 4) mc.codec = static_cast<int>(T.zigzag());
                                        else if (id4 == 5) mc.num_valores = T.zigzag();
                                        else if (id4 == 7) mc.comprimido = T.zigzag();
                                        else if (id4 == 9) mc.data_off = T.zigzag();
                                        else if (id4 == 11) mc.dict_off = T.zigzag();
                                        else T.saltar(t4);
                                    }
                                } else T.saltar(t3);
                            }
                            cols.push_back(mc);
                        }
                    } else T.saltar(t2);
                }
                grupos.push_back(cols);
            }
        } else if (id == 5 && t == 9) {
            auto [te, k] = T.lista();
            (void)te;
            for (size_t e = 0; e < k; ++e) {
                std::string kk, vv;
                int u = 0;
                for (;;) {
                    auto [id2, t2] = T.campo(u);
                    if (t2 == 0) break;
                    if (id2 == 1) kk = T.cadena();
                    else if (id2 == 2) vv = T.cadena();
                    else T.saltar(t2);
                }
                tabla.metadatos.emplace_back(kk, vv);
            }
        } else {
            T.saltar(t);
        }
    }
    if (esquema.empty()) throw ErrorFormato("parquet: esquema vacio");
    for (size_t e = 1; e < esquema.size(); ++e) {
        if (esquema[e].hijos > 0) throw ErrorFormato("parquet: columnas anidadas no soportadas");
        if (esquema[e].repeticion == 2) throw ErrorFormato("parquet: columnas repetidas no soportadas");
        ColumnaParquet c;
        c.nombre = esquema[e].nombre;
        c.tipo = static_cast<TipoFisico>(esquema[e].tipo);
        c.opcional = esquema[e].repeticion == 1;
        c.unidad_timestamp = esquema[e].unidad;
        c.timestamp_ajustado_utc = esquema[e].utc;
        tabla.columnas.push_back(c);
    }

    for (const auto& g : grupos) {
        if (g.size() != tabla.columnas.size()) throw ErrorFormato("parquet: columnas del grupo distintas del esquema");
        for (size_t ci = 0; ci < g.size(); ++ci) {
            const MetaColumna& mc = g[ci];
            ColumnaParquet& col = tabla.columnas[ci];
            const TipoFisico tf = col.tipo;
            const size_t ancho = (tf == TipoFisico::DOUBLE || tf == TipoFisico::INT64) ? 8u : 4u;
            if (!(tf == TipoFisico::DOUBLE || tf == TipoFisico::FLOAT || tf == TipoFisico::INT64 || tf == TipoFisico::INT32))
                throw ErrorFormato("parquet: tipo fisico no soportado en " + col.nombre);
            std::vector<uint8_t> diccionario;
            size_t n_dic = 0;
            int64_t leidos = 0;
            size_t pos = static_cast<size_t>(mc.dict_off >= 0 && mc.dict_off < mc.data_off ? mc.dict_off : mc.data_off);
            auto empujar = [&](const uint8_t* v, bool es_nulo) {
                col.nulo.push_back(es_nulo ? 1 : 0);
                if (tf == TipoFisico::DOUBLE) { double x = std::numeric_limits<double>::quiet_NaN(); if (!es_nulo) std::memcpy(&x, v, 8); col.dobles.push_back(x); }
                else if (tf == TipoFisico::FLOAT) { float x = std::numeric_limits<float>::quiet_NaN(); if (!es_nulo) std::memcpy(&x, v, 4); col.dobles.push_back(static_cast<double>(x)); }
                else if (tf == TipoFisico::INT64) { int64_t x = 0; if (!es_nulo) std::memcpy(&x, v, 8); col.enteros.push_back(x); }
                else { int32_t x = 0; if (!es_nulo) std::memcpy(&x, v, 4); col.enteros.push_back(x); }
            };
            while (leidos < mc.num_valores) {
                LectorThrift H(f.data(), n, pos);
                int tipo_pag = -1, tam_desc = 0, tam_comp = 0, nvals = 0, codif = 0, nnulos = 0, lon_def = 0, lon_rep = 0;
                bool v2_comprimido = true;
                int u = 0;
                for (;;) {
                    auto [id, t] = H.campo(u);
                    if (t == 0) break;
                    if (id == 1) tipo_pag = static_cast<int>(H.zigzag());
                    else if (id == 2) tam_desc = static_cast<int>(H.zigzag());
                    else if (id == 3) tam_comp = static_cast<int>(H.zigzag());
                    else if ((id == 5 || id == 7 || id == 8) && t == 12) {
                        int u2 = 0;
                        for (;;) {
                            auto [id2, t2] = H.campo(u2);
                            if (t2 == 0) break;
                            if (id2 == 1) nvals = static_cast<int>(H.zigzag());
                            else if (id == 8 && id2 == 2) nnulos = static_cast<int>(H.zigzag());
                            else if (id == 8 && id2 == 4) codif = static_cast<int>(H.zigzag());
                            else if (id == 8 && id2 == 5) lon_def = static_cast<int>(H.zigzag());
                            else if (id == 8 && id2 == 6) lon_rep = static_cast<int>(H.zigzag());
                            else if (id == 8 && id2 == 7) v2_comprimido = (t2 == 1);
                            else if (id != 8 && id2 == 2) codif = static_cast<int>(H.zigzag());
                            else H.saltar(t2);
                        }
                    } else H.saltar(t);
                }
                const size_t ini = H.i;
                if (ini + static_cast<size_t>(tam_comp) > n) throw ErrorFormato("parquet: pagina truncada");
                pos = ini + static_cast<size_t>(tam_comp);
                if (tipo_pag == 2) {   // DICTIONARY_PAGE
                    diccionario = descomprimir(mc.codec, &f[ini], static_cast<size_t>(tam_comp));
                    n_dic = static_cast<size_t>(nvals);
                    continue;
                }
                if (tipo_pag != 0 && tipo_pag != 3) continue;   // indices u otras: no aportan valores
                std::vector<uint8_t> cuerpo;
                std::vector<uint32_t> defs;
                size_t q = 0;
                if (tipo_pag == 0) {
                    cuerpo = descomprimir(mc.codec, &f[ini], static_cast<size_t>(tam_comp));
                    if (static_cast<int>(cuerpo.size()) != tam_desc) throw ErrorFormato("parquet: tamano de pagina distinto");
                    if (col.opcional) {
                        const size_t L = le32(&cuerpo[0]);
                        defs = rle_hibrido(&cuerpo[4], L, 1, static_cast<size_t>(nvals));
                        q = 4 + L;
                    }
                } else {
                    const size_t niv = static_cast<size_t>(lon_def + lon_rep);
                    if (lon_rep) throw ErrorFormato("parquet: niveles de repeticion no soportados");
                    if (col.opcional && lon_def) defs = rle_hibrido(&f[ini], static_cast<size_t>(lon_def), 1, static_cast<size_t>(nvals));
                    std::vector<uint8_t> vals = v2_comprimido ? descomprimir(mc.codec, &f[ini + niv], static_cast<size_t>(tam_comp) - niv)
                                                              : std::vector<uint8_t>(&f[ini + niv], &f[ini + static_cast<size_t>(tam_comp)]);
                    cuerpo = std::move(vals);
                    (void)nnulos;
                }
                const size_t presentes = defs.empty() ? static_cast<size_t>(nvals)
                                                      : static_cast<size_t>(std::count(defs.begin(), defs.end(), 1u));
                if (codif == 0) {   // PLAIN
                    if (q + presentes * ancho > cuerpo.size()) throw ErrorFormato("parquet: valores PLAIN truncados");
                    size_t k = 0;
                    for (int r = 0; r < nvals; ++r) {
                        const bool nul = !defs.empty() && defs[static_cast<size_t>(r)] == 0;
                        empujar(nul ? nullptr : &cuerpo[q + (k++) * ancho], nul);
                    }
                } else if (codif == 2 || codif == 8) {   // PLAIN_DICTIONARY / RLE_DICTIONARY
                    if (diccionario.empty()) throw ErrorFormato("parquet: pagina con diccionario sin diccionario");
                    const int bw = cuerpo[q];
                    const std::vector<uint32_t> idx = rle_hibrido(&cuerpo[q + 1], cuerpo.size() - q - 1, bw, presentes);
                    size_t k = 0;
                    for (int r = 0; r < nvals; ++r) {
                        const bool nul = !defs.empty() && defs[static_cast<size_t>(r)] == 0;
                        if (nul) { empujar(nullptr, true); continue; }
                        const uint32_t ix = idx[k++];
                        if (ix >= n_dic) throw ErrorFormato("parquet: indice de diccionario fuera de rango");
                        empujar(&diccionario[ix * ancho], false);
                    }
                } else {
                    throw ErrorFormato("parquet: codificacion no soportada (" + std::to_string(codif) + ") en " + col.nombre);
                }
                leidos += nvals;
            }
        }
    }
    for (const auto& c : tabla.columnas)
        if (static_cast<int64_t>(c.nulo.size()) != tabla.filas) throw ErrorFormato("parquet: filas leidas distintas en " + c.nombre);
    return tabla;
}

}  // namespace ark::dataset::formatos
