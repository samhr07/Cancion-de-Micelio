// motor_ventanas.hpp -- calca C++ de ARQUITECTO DE DATOS/filtros_de_ventana (contratos_ventanas.py y motor_ventanas.py, ARK-WIN-01):
// el motor causal de ventanas rodantes en memoria O(W) que consume las tablas del dataset (motor_ingesta.hpp) y arreglos sueltos.
// Mismas reglas, mismos mensajes y mismas sumas que el Python (np.sum por pares en 1-D, secuencial por columna en axis=0 con
// M > 1, cumsum secuencial, el lazo de Shannon de oraculo_lazo.hpp). Lo que en Python es un generador aqui es una retrollamada
// por ventana. Sin identidad bit a bit, declarado: la rama TV (ADMM con LAPACK) y el proyector aniquilador (inv/pinv/cond de
// LAPACK), que se contrastan con tolerancia; auditar_huella_memoria_streaming (tracemalloc) cuenta ventanas y no mide RAM.
// Jefe de Arquitectos, 2026-11-03.
// CORRECCION D1-D6 (2026-11-03, decisiones del Comandante; acta PORTE_INGESTA_Y_VENTANAS_CXX.md §4). Desde aqui la cabecera ya
// no es la calca del principal sino la del oraculo corregido PRUEBAS_JEFE/ventanas_corregido/filtros_de_ventana/:
//   D1 el testigo pasa a float las columnas enteras de su copia antes del choque;
//   D2 sin `timestamp_utc` ni `t`, el sello es `fecha_utc` en segundos de epoca;
//   D3 un DatetimeIndex se lleva a segundos sea cual sea su unidad;
//   D4 la mascara del vector de reloj se convierte a booleana al entrar;
//   D5 un indice no entero (texto incluido) no tiene sellos: se numera;
//   D6 un None de venue/nivel/simbolo es nulo tambien al principio del stream.
#pragma once

#include "ark/dataset/control_calidad.hpp"
#include "ark/dataset/deconvolucion_tv.hpp"
#include "ark/dataset/errores_py.hpp"
#include "ark/dataset/guardas_causales.hpp"
#include "ark/dataset/numpy_formato.hpp"
#include "ark/dataset/numpy_rng.hpp"
#include "ark/dataset/oraculo_lazo.hpp"
#include "ark/dataset/tabla_ingesta.hpp"
#include "ark/dataset/zona_horaria.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace ark::dataset::ventanas {

using ingesta::Columna;
using ingesta::JV;
using ingesta::Tabla;
using ingesta::TipoColumna;

// ============================================================================================== contratos_ventanas.py: excepciones
ARK_ERROR_PY(ErrorFiltroVentana, ErrorPy)
ARK_ERROR_PY(ErrorViolacionCausalidad, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorMatrizNoCausal, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorMonotoniaTiempo, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorContratoInvalido, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorSincronizacionInvalida, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorDimensionInsuficiente, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorParametroInvalido, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorMascaraImputacionExcedida, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorRelojIntrinsecoInvalido, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorAislamientoMemoria, ErrorFiltroVentana)
ARK_ERROR_PY(ErrorResiduoTemporalInvalido, ErrorFiltroVentana)

enum class TipoRelojVentana : uint8_t { PARED, NU_TRADES, NU_VOLUMEN, NU_NOCIONAL, NU_INFO_SHANNON };
enum class PoliticaSincronizacion : uint8_t { INTERSECCION_ESTRICTA, UNION_CON_MASCARA };

inline const char* nombre(TipoRelojVentana r) {
    switch (r) {
        case TipoRelojVentana::PARED: return "PARED";
        case TipoRelojVentana::NU_TRADES: return "NU_TRADES";
        case TipoRelojVentana::NU_VOLUMEN: return "NU_VOLUMEN";
        case TipoRelojVentana::NU_NOCIONAL: return "NU_NOCIONAL";
        default: return "NU_INFO_SHANNON";
    }
}
inline const char* nombre(PoliticaSincronizacion p) {
    return p == PoliticaSincronizacion::INTERSECCION_ESTRICTA ? "INTERSECCION_ESTRICTA" : "UNION_CON_MASCARA";
}

inline const std::vector<std::string>& metodos_imputacion() {
    static const std::vector<std::string> m = {"FORWARD_FILL", "CONGELAR_RELOJ", "PROHIBIR_INTERPOLACION"};
    return m;
}

inline const std::map<std::string, std::set<std::string>>& matriz_compatibilidad() {
    static const std::map<std::string, std::set<std::string>> m = {
        {"YAHOO", {"L1", "L1_OHLCV"}},
        {"BLOOMBERG", {"L1", "L1_OHLCV"}},
        {"BLOOMBERG_XLSX", {"L1", "L1_OHLCV"}},
        {"BINANCE", {"L1", "L1_OHLCV", "L1_ENRIQUECIDO", "L2", "L2_DEPTH"}},
        {"BINANCE_UM", {"L1", "L1_OHLCV", "L1_ENRIQUECIDO", "L2", "L2_DEPTH"}},
        {"BINANCE_SPOT", {"L1", "L1_OHLCV", "L1_ENRIQUECIDO", "L2", "L2_DEPTH"}},
        {"BITFINEX", {"L1", "L1_OHLCV", "L2", "L2_MID", "L2_DEPTH", "L3", "L3_MBO"}},
        {"INTERACTIVE_BROKERS", {"L1", "L1_OHLCV", "L2", "L2_DEPTH"}},
        {"IBKR", {"L1", "L1_OHLCV", "L2", "L2_DEPTH"}},
        {"IEX", {"TOPS", "L1", "L1_OHLCV"}},
        {"DUKASCOPY", {"L1", "L1_OHLCV", "TOPS"}},
    };
    return m;
}

namespace detalle_vent {

inline const std::set<std::string>& tokens_nulos() {
    static const std::set<std::string> t = {"", "NAN", "NONE", "NULL", "NAT", "<NA>"};
    return t;
}
inline bool es_espacio(unsigned char c) { return c == ' ' || (c >= 0x09 && c <= 0x0d) || (c >= 0x1c && c <= 0x1f); }
inline std::string strip(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && es_espacio(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && es_espacio(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}
inline std::string upper(std::string s) {
    for (auto& c : s) if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    return s;
}
inline bool es_nulo_token(const std::string& s) { return tokens_nulos().count(upper(strip(s))) > 0; }
inline std::string lista_ordenada(const std::set<std::string>& s) {
    std::vector<std::string> v(s.begin(), s.end());
    return ingesta::repr_lista_str(v);
}
inline std::string i2s(int64_t v) { return std::to_string(v); }
inline std::string f2s(double v) { return ingesta::repr_py(v); }   // str(float) de Python
inline std::string pct1(double v) { return ingesta::fmt_f(v, 1); }

// Division entera por piso de Python (int // int).
inline int64_t div_piso(int64_t a, int64_t b) {
    int64_t q = a / b;
    if ((a % b != 0) && ((a < 0) != (b < 0))) --q;
    return q;
}
// float // float de Python (float_floor_div de CPython).
inline double div_piso_f(double vx, double wx) {
    double mod = std::fmod(vx, wx);
    double div = (vx - mod) / wx;
    if (mod != 0.0) {
        if ((wx < 0) != (mod < 0)) { mod += wx; div -= 1.0; }
    }
    double fd;
    if (div != 0.0) {
        fd = std::floor(div);
        if (div - fd > 0.5) fd += 1.0;
    } else {
        fd = std::copysign(0.0, vx / wx);
    }
    return fd;
}
// np.sum de un tramo contiguo 1-D (0 + suma por pares).
inline double suma_np(const double* a, size_t n) { return 0.0 + ingesta::suma_pares(a, n); }

}  // namespace detalle_vent

// validar_compatibilidad_venue_nivel (Axioma 2).
inline void validar_compatibilidad_venue_nivel(const std::optional<std::string>& venue = std::nullopt,
                                               const std::optional<std::string>& nivel = std::nullopt,
                                               const std::optional<std::string>& simbolo = std::nullopt) {
    using namespace detalle_vent;
    if (venue && es_nulo_token(*venue))
        throw ErrorContratoInvalido("VIOLACIÓN AXIOMA 2: Identificador de 'venue' inválido, nulo o vacío: " + ingesta::repr_str_py(*venue) + ".");
    if (nivel && es_nulo_token(*nivel))
        throw ErrorContratoInvalido("VIOLACIÓN AXIOMA 2: Identificador de 'nivel' inválido, nulo o vacío: " + ingesta::repr_str_py(*nivel) + ".");
    if (simbolo && es_nulo_token(*simbolo))
        throw ErrorContratoInvalido("VIOLACIÓN AXIOMA 2: Identificador de 'simbolo' inválido, nulo o vacío: " + ingesta::repr_str_py(*simbolo) + ".");
    if (venue && nivel) {
        const std::string vn = upper(strip(*venue)), nn = upper(strip(*nivel));
        const auto& M = matriz_compatibilidad();
        const auto it = M.find(vn);
        if (it != M.end() && !it->second.count(nn))
            throw ErrorContratoInvalido("VIOLACIÓN AXIOMA 2 (Incompatibilidad Física Venue x Nivel): El venue '" + *venue +
                                        "' NO soporta físicamente el nivel de microestructura '" + *nivel + "'. Niveles permitidos para " + vn +
                                        ": " + lista_ordenada(it->second) + ".");
    }
}

// ConfiguracionVentanaCausal (dataclass congelada); validar() = __post_init__.
struct ConfiguracionVentanaCausal {
    int64_t lookback_W = 240;
    int64_t stride_S = 1;
    TipoRelojVentana reloj = TipoRelojVentana::PARED;
    bool forzar_divisores_exactos = true;
    bool bidireccional_kantelhardt = true;
    bool deconvolucion_tv_activa = false;
    double lambda_tv = 1.0;
    double umbral_imputacion_maximo = 0.20;
    int64_t min_muestras_requeridas = 32;
    int64_t s_min = 4;
    double max_ratio_escalas = 4.0;
    std::string metodo_imputacion = "FORWARD_FILL";   // el valor (str) del miembro o la cadena suelta
    bool metodo_es_enum = true;                        // MetodoImputacion(...) o str suelta
    bool exigir_monotonia_estricta = true;
    bool abortar_si_excede_imputacion = false;
    bool permitir_sort_automatico = false;
    PoliticaSincronizacion politica_sincronizacion = PoliticaSincronizacion::INTERSECCION_ESTRICTA;
    std::optional<std::string> venue_esperado, nivel_esperado, simbolo_esperado;
    std::optional<std::vector<std::string>> activos_esperados;
    bool congelar_en_huecos = true;
    bool estampar_al_cierre = false;
    std::optional<int64_t> delta_t_barra_segundos;
    double tolerancia_jitter_segundos = 0.0;
    bool proteger_escritura_vistas = true;
    int64_t max_elementos_materializables = 10000000;
    int64_t orden_polinomio_detrending = 1;
    std::optional<int64_t> ventana_innovacion_shannon;
    bool lazo_feedback_shannon = true;
    std::optional<double> tamano_tick, precio_referencia;
    double sigma_floor_minimo = 1e-4;

    void validar() const {
        using namespace detalle_vent;
        if (lookback_W < min_muestras_requeridas)
            throw ErrorParametroInvalido("lookback_W=" + i2s(lookback_W) + " no puede ser menor que min_muestras_requeridas=" +
                                         i2s(min_muestras_requeridas) + ".");
        if (stride_S < 1) throw ErrorParametroInvalido("stride_S=" + i2s(stride_S) + " debe ser estrictamente >= 1.");
        if (lambda_tv <= 0) throw ErrorParametroInvalido("lambda_tv=" + f2s(lambda_tv) + " debe ser estrictamente positivo (> 0).");
        if (!(0.0 <= umbral_imputacion_maximo && umbral_imputacion_maximo <= 1.0))
            throw ErrorParametroInvalido("umbral_imputacion_maximo=" + f2s(umbral_imputacion_maximo) + " debe estar en [0.0, 1.0].");
        if (s_min < 2) throw ErrorParametroInvalido("s_min=" + i2s(s_min) + " debe ser >= 2 para garantizar al menos 2 grados de libertad.");
        if (max_ratio_escalas <= 1.0)
            throw ErrorParametroInvalido("max_ratio_escalas=" + f2s(max_ratio_escalas) + " debe ser estrictamente > 1.0.");
        if (tamano_tick && *tamano_tick <= 0.0)
            throw ErrorParametroInvalido("tamano_tick=" + f2s(*tamano_tick) + " debe ser estrictamente positivo (> 0).");
        if (precio_referencia && *precio_referencia <= 0.0)
            throw ErrorParametroInvalido("precio_referencia=" + f2s(*precio_referencia) + " debe ser estrictamente positivo (> 0).");
        if (sigma_floor_minimo <= 0.0)
            throw ErrorParametroInvalido("sigma_floor_minimo=" + f2s(sigma_floor_minimo) + " debe ser estrictamente positivo (> 0).");
        if (ventana_innovacion_shannon && (*ventana_innovacion_shannon < 2 || *ventana_innovacion_shannon > lookback_W))
            throw ErrorParametroInvalido("ventana_innovacion_shannon=" + i2s(*ventana_innovacion_shannon) + " debe estar en [2, lookback_W=" +
                                         i2s(lookback_W) + "].");
        if (max_elementos_materializables <= 0)
            throw ErrorParametroInvalido("max_elementos_materializables=" + i2s(max_elementos_materializables) +
                                         " debe ser estrictamente positivo (> 0).");
        if (orden_polinomio_detrending < 0)
            throw ErrorParametroInvalido("orden_polinomio_detrending=" + i2s(orden_polinomio_detrending) + " debe ser >= 0.");
        if (permitir_sort_automatico)
            throw ErrorViolacionCausalidad("VIOLACIÓN AXIOMA 2: Prohibido activar 'permitir_sort_automatico'. "
                                           "Reordenar pasivamente un stream de feed en vivo encubre discontinuidades de red y produce lookahead bias.");
        if (!exigir_monotonia_estricta)
            throw ErrorMonotoniaTiempo("VIOLACIÓN AXIOMA 2: Prohibido desactivar 'exigir_monotonia_estricta'. "
                                       "La monotonía temporal estricta (t_k > t_{k-1}) es una ley física innegociable del Reactor ARK.");
        if (delta_t_barra_segundos && *delta_t_barra_segundos <= 0)
            throw ErrorParametroInvalido("delta_t_barra_segundos=" + i2s(*delta_t_barra_segundos) + " debe ser estrictamente positivo (> 0).");
        if (tolerancia_jitter_segundos < 0.0)
            throw ErrorParametroInvalido("tolerancia_jitter_segundos=" + f2s(tolerancia_jitter_segundos) + " no puede ser negativa.");
        validar_compatibilidad_venue_nivel(venue_esperado, nivel_esperado, simbolo_esperado);
        const std::string mv = metodo_es_enum ? metodo_imputacion : upper(metodo_imputacion);
        const auto& vm = metodos_imputacion();
        if (std::find(vm.begin(), vm.end(), mv) == vm.end()) {
            std::set<std::string> s(vm.begin(), vm.end());
            throw ErrorViolacionCausalidad("VIOLACIÓN AXIOMA 1: El método de imputación '" + metodo_imputacion +
                                           "' viola la causalidad temporal al propagar o interpolar datos del futuro hacia el pasado. Solo se admiten: " +
                                           lista_ordenada(s) + ".");
        }
    }

    double obtener_sigma_floor(std::optional<double> precio_actual = std::nullopt) const {
        if (tamano_tick && *tamano_tick > 0) {
            const std::optional<double> pref = (precio_actual && *precio_actual > 0) ? precio_actual : precio_referencia;
            if (pref && *pref > 0) {
                const double st = *tamano_tick / *pref;
                return st > sigma_floor_minimo ? st : sigma_floor_minimo;
            }
        }
        return sigma_floor_minimo;
    }
};

// ======================================================================================== entradas: arreglos, marcos y trozos
enum class Dt : uint8_t { F64, I64, BOOL, OBJ, USTR };

// Celda de un arreglo object de numpy (str, None o float) o de un arreglo de np.str_.
struct Obj {
    enum Tipo : uint8_t { STR, NONE, FLOAT } t = NONE;
    std::string s;
    double f = 0.0;
};

// Un np.ndarray 1-D (o con `forma` 2-D/0-D para los datos).
struct Arreglo {
    Dt dt = Dt::F64;
    std::vector<double> f;
    std::vector<int64_t> i;    // I64 y BOOL (0/1)
    std::vector<Obj> o;        // OBJ y USTR
    std::vector<int64_t> forma;  // vacio = 1-D; {} con escalar = true: 0-D
    bool escalar = false;
    size_t size() const { return dt == Dt::F64 ? f.size() : ((dt == Dt::OBJ || dt == Dt::USTR) ? o.size() : i.size()); }
    static Arreglo de_f64(std::vector<double> v) { Arreglo a; a.dt = Dt::F64; a.f = std::move(v); return a; }
    static Arreglo de_i64(std::vector<int64_t> v) { Arreglo a; a.dt = Dt::I64; a.i = std::move(v); return a; }
    static Arreglo de_bool(std::vector<int64_t> v) { Arreglo a; a.dt = Dt::BOOL; a.i = std::move(v); return a; }
    static Arreglo de_obj(std::vector<Obj> v, bool ustr = false) { Arreglo a; a.dt = ustr ? Dt::USTR : Dt::OBJ; a.o = std::move(v); return a; }
};

enum class TipoIndice : uint8_t { RANGO, ENTERO, FECHA, OTRO, STR };   // STR: indice de texto de pandas 3 (StringDtype)
struct Indice {
    TipoIndice tipo = TipoIndice::RANGO;
    int64_t inicio = 0;               // RANGO: inicio, inicio + 1, ...
    std::vector<int64_t> valores;     // ENTERO (valores) o FECHA (ns UTC)
    int64_t unidad_ns = 1;            // FECHA: ns por unidad del dtype (astype(int64) de un DatetimeIndex[us] da microsegundos)
};

// Un pd.DataFrame: tabla + indice. Las columnas TS de la tabla van en ns; `unidad_ts_ns[c]` = ns por unidad del dtype de pandas
// (to_numpy(int64) de un datetime64[us] da microsegundos).
struct Marco {
    Tabla tabla;
    Indice indice;
    std::map<std::string, int64_t> unidad_ts_ns;
    size_t filas() const { return tabla.filas(); }
};

// Un valor de un dict de trozo: None, un arreglo (o pd.Series, cuyo indice el motor ignora en ese camino) o un DataFrame.
struct Valor {
    enum Tipo : uint8_t { NADA, ARR, MARCO } t = NADA;
    Arreglo a;
    Marco m;
    static Valor nada() { return {}; }
    static Valor arr(Arreglo x) { Valor v; v.t = ARR; v.a = std::move(x); return v; }
    static Valor marco(Marco x) { Valor v; v.t = MARCO; v.m = std::move(x); return v; }
};

// Un trozo del generador: dict, DataFrame, Series o arreglo.
struct Trozo {
    enum Tipo : uint8_t { DICT, MARCO, SERIE, ARREGLO } t = ARREGLO;
    std::vector<std::pair<std::string, Valor>> dict;
    Marco m;
    Arreglo a;
    Indice indice;   // SERIE
    const Valor* buscar(const std::string& k) const {
        for (const auto& kv : dict) if (kv.first == k) return &kv.second;
        return nullptr;
    }
};

// Entrada de generar_ventanas: un DataFrame, una Series o un arreglo/lista.
struct Datos {
    enum Tipo : uint8_t { MARCO, SERIE, ARREGLO } t = ARREGLO;
    Marco m;
    Arreglo a;
    Indice indice;
    size_t len() const { return t == MARCO ? m.filas() : (a.forma.empty() ? a.size() : static_cast<size_t>(a.forma[0])); }
};

// Un activo de generar_ventanas_multiactivo: DataFrame, Series o arreglo.
struct ItemActivo {
    enum Tipo : uint8_t { MARCO, SERIE, ARREGLO } t = ARREGLO;
    Marco m;
    Arreglo a;
    Indice indice;
};
using SeriesDict = std::vector<std::pair<std::string, ItemActivo>>;
using ArreglosDict = std::vector<std::pair<std::string, Arreglo>>;

// ================================================================================================ VentanaCausalCertificada
struct VentanaCausalCertificada {
    int64_t t_actual = 0;
    std::vector<double> valores;           // fila mayor; 1-D si m == 0
    std::vector<uint8_t> mascara_huecos;
    int64_t w = 0, m = 0;                  // m = 0: arreglo 1-D de largo w; m >= 1: matriz (w, m)
    double reloj_acumulado = 0.0;
    std::vector<int64_t> divisores_exactos;
    int64_t residuo_temporal_maximo = 0;
    double porcentaje_descarte_causal = 0.0;
    double lookahead_score = 0.0;
    bool es_valida = true;
    std::string venue = "GENERICO", nivel = "L1", simbolo = "SYNTH";
    std::vector<std::string> activos = {"SYNTH"};
    std::optional<std::vector<uint8_t>> market_open_mask;
    TipoRelojVentana tipo_reloj = TipoRelojVentana::PARED;
    double delta_reloj_ventana = 0.0;
    int64_t barras_omitidas = 0;
    std::vector<std::string> advertencias;
    std::vector<std::pair<std::pair<std::string, std::string>, std::string>> clasificacion_pares;
    std::optional<std::vector<double>> delta_reloj_por_activo;

    bool tiene_rezago_causal() const { for (const auto& kv : clasificacion_pares) if (kv.second == "REZAGO") return true; return false; }
    int64_t num_pares_rezago() const { int64_t n = 0; for (const auto& kv : clasificacion_pares) if (kv.second == "REZAGO") ++n; return n; }
    std::string obtener_clase_par(const std::string& u, const std::string& v) const {
        for (const auto& kv : clasificacion_pares) if (kv.first.first == u && kv.first.second == v) return kv.second;
        return "NO_VERIFICADO";
    }
    int64_t longitud() const { return w; }
    bool tiene_discontinuidad_temporal() const { return barras_omitidas > 0; }
    int64_t num_puntos_imputados() const {
        if (w == 0) return 0;
        int64_t n = 0; for (auto b : mascara_huecos) n += b ? 1 : 0; return n;
    }
    double densidad_perdida_total() const {
        const int64_t tot = static_cast<int64_t>(mascara_huecos.size()) + barras_omitidas;
        if (tot == 0) return 0.0;
        return static_cast<double>(num_puntos_imputados() + barras_omitidas) / static_cast<double>(tot);
    }
    bool es_multiactivo() const { return m > 1; }
    int64_t num_activos() const { return m >= 1 ? m : 1; }
    double porcentaje_imputado() const {
        if (w == 0) return 0.0;
        return 100.0 * static_cast<double>(num_puntos_imputados()) / static_cast<double>(mascara_huecos.size());
    }
    double ratio_imputacion() const {
        if (w == 0) return 0.0;
        return static_cast<double>(num_puntos_imputados()) / static_cast<double>(mascara_huecos.size());
    }
    // valor_presente: escalar (1-D) o la ultima fila (matriz).
    std::vector<double> valor_presente() const {
        if (m == 0) return {valores.back()};
        return std::vector<double>(valores.end() - m, valores.end());
    }
    std::vector<double> normalizar_zscore_causal(double eps = 1e-12) const {
        std::vector<double> z(valores.size(), 0.0);
        if (m == 0) {
            const size_t n = valores.size();
            const double mu = detalle_vent::suma_np(valores.data(), n) / static_cast<double>(n);
            std::vector<double> d(n);
            for (size_t k = 0; k < n; ++k) d[k] = (valores[k] - mu) * (valores[k] - mu);
            const double sigma = std::sqrt(detalle_vent::suma_np(d.data(), n) / static_cast<double>(n));
            if (!(sigma < eps)) for (size_t k = 0; k < n; ++k) z[k] = (valores[k] - mu) / sigma;
            return z;
        }
        const size_t W = static_cast<size_t>(w), M = static_cast<size_t>(m);
        std::vector<double> mu(M), sg(M);
        for (size_t j = 0; j < M; ++j) mu[j] = suma_eje0(valores, W, M, j, false) / static_cast<double>(W);
        std::vector<double> d(valores.size());
        for (size_t r = 0; r < W; ++r) for (size_t j = 0; j < M; ++j) {
            const double x = valores[r * M + j] - mu[j];
            d[r * M + j] = x * x;
        }
        for (size_t j = 0; j < M; ++j) sg[j] = std::sqrt(suma_eje0(d, W, M, j, false) / static_cast<double>(W));
        for (size_t r = 0; r < W; ++r) for (size_t j = 0; j < M; ++j) {
            const double den = sg[j] > eps ? sg[j] : (std::isnan(sg[j]) ? sg[j] : eps);   // np.maximum(sigma, eps)
            z[r * M + j] = (sg[j] < eps) ? 0.0 : (valores[r * M + j] - mu[j]) / den;
        }
        return z;
    }
    // np.sum(x, axis=0) de una matriz (W, M) contigua: M == 1 por pares; M > 1 secuencial por columna.
    static double suma_eje0(const std::vector<double>& x, size_t W, size_t M, size_t j, bool) {
        if (M == 1) return detalle_vent::suma_np(x.data(), W);
        double s = 0.0;
        for (size_t r = 0; r < W; ++r) s += x[r * M + j];
        return s;
    }
};
using Emisor = std::function<void(const VentanaCausalCertificada&)>;

// ================================================================================================= funciones del Axioma 1, 3 y 6
// _resolver_fecha_sesion: fecha de sesion (yyyymmdd) del sello de una barra.
inline uint32_t resolver_fecha_sesion(int64_t t) {
    const double tn = static_cast<double>(t);
    constexpr int64_t NSD = 86400LL * 1000000000LL;
    auto fecha_ns = [&](double base_unidades, int64_t m, int p) -> uint32_t {
        const int64_t base = static_cast<int64_t>(base_unidades);
        double frac = base_unidades - static_cast<double>(base);
        if (p) frac = ingesta::round_py(frac, p);
        const int64_t ns = base * m + static_cast<int64_t>(frac * static_cast<double>(m));
        return yyyymmdd_desde_dias(detalle_vent::div_piso(ns, NSD));
    };
    if (tn > 1e14) return fecha_ns(tn, 1, 0);
    if (tn > 1e11) return fecha_ns(tn, 1000000, 6);
    if (tn > 1e8) return fecha_ns(tn, 1000000000, 9);
    return yyyymmdd_desde_dias(static_cast<int64_t>(detalle_vent::div_piso_f(tn, 86400.0)));
}
inline std::string texto_yyyymmdd(uint32_t f) {
    char b[16];
    std::snprintf(b, sizeof b, "%04u-%02u-%02u", f / 10000, (f / 100) % 100, f % 100);
    return b;
}

// validar_matriz_causal: A (fila mayor, filas x cols) estrictamente triangular inferior y finita.
inline bool validar_matriz_causal(const std::vector<double>& A, int64_t filas, int64_t cols, double tol = 1e-15) {
    if (filas != cols)
        throw ErrorParametroInvalido("A debe ser una matriz cuadrada, recibido shape=(" + std::to_string(filas) + ", " + std::to_string(cols) + ")");
    for (double v : A)
        if (!std::isfinite(v))
            throw ErrorMatrizNoCausal("VIOLACIÓN AXIOMA 1: La matriz A contiene valores no finitos (NaN o Inf), "
                                      "impidiendo certificar su estructura causal triangular inferior.");
    double mx = 0.0;
    bool hay = false;
    for (int64_t i = 0; i < filas; ++i)
        for (int64_t j = i + 1; j < cols; ++j) {
            const double a = std::fabs(A[static_cast<size_t>(i * cols + j)]);
            if (a > tol) hay = true;
            if (a > mx) mx = a;
        }
    if (hay) {
        char b[64];
        std::snprintf(b, sizeof b, "%.2e", mx);
        throw ErrorMatrizNoCausal(std::string("VIOLACIÓN AXIOMA 1: La matriz A no es causal triangular inferior. "
                                              "Fuga retrospectiva detectada en el triángulo superior: max|A[i, j > i]|=") + b + ".");
    }
    return true;
}

inline std::vector<int64_t> calcular_divisores_exactos(int64_t W, int64_t s_min = 4, double max_ratio = 4.0) {
    const int64_t s_max = std::max(s_min, static_cast<int64_t>(detalle_vent::div_piso_f(static_cast<double>(W), max_ratio)));
    std::vector<int64_t> d;
    for (int64_t s = s_min; s <= s_max; ++s) if (W % s == 0) d.push_back(s);
    return d;
}

inline int64_t validar_escala_multiescala(int64_t W, int64_t s, bool forzar_divisores_exactos = true) {
    if (s < 2) throw ErrorParametroInvalido("Escala s=" + std::to_string(s) + " debe ser >= 2.");
    if (s > W) throw ErrorParametroInvalido("Escala s=" + std::to_string(s) + " no puede ser mayor que la ventana W=" + std::to_string(W) + ".");
    const int64_t r = W % s;
    if (forzar_divisores_exactos && r != 0)
        throw ErrorResiduoTemporalInvalido("VIOLACIÓN AXIOMA 6: Escala s=" + std::to_string(s) + " no divide exactamente a W=" + std::to_string(W) +
                                           " (residuo=" + std::to_string(r) + "). Bajo forzar_divisores_exactos=True, el residuo temporal debe ser "
                                           "idénticamente 0 para evitar mutilación del instante presente.");
    return r;
}

namespace detalle_vent {
// Inversa de una matriz simetrica pequena por Gauss-Jordan con pivote parcial (sustituto de np.linalg.inv / pinv).
inline bool invertir(std::vector<double> a, size_t n, std::vector<double>& inv) {
    inv.assign(n * n, 0.0);
    for (size_t i = 0; i < n; ++i) inv[i * n + i] = 1.0;
    for (size_t c = 0; c < n; ++c) {
        size_t p = c;
        for (size_t r = c + 1; r < n; ++r) if (std::fabs(a[r * n + c]) > std::fabs(a[p * n + c])) p = r;
        if (a[p * n + c] == 0.0) return false;
        if (p != c) for (size_t k = 0; k < n; ++k) { std::swap(a[c * n + k], a[p * n + k]); std::swap(inv[c * n + k], inv[p * n + k]); }
        const double d = a[c * n + c];
        for (size_t k = 0; k < n; ++k) { a[c * n + k] /= d; inv[c * n + k] /= d; }
        for (size_t r = 0; r < n; ++r) {
            if (r == c) continue;
            const double f = a[r * n + c];
            if (f == 0.0) continue;
            for (size_t k = 0; k < n; ++k) { a[r * n + k] -= f * a[c * n + k]; inv[r * n + k] -= f * inv[c * n + k]; }
        }
    }
    return true;
}
// Autovalores de una matriz simetrica (Jacobi) para cond_2 = l_max / l_min.
inline std::vector<double> autovalores_sim(std::vector<double> a, size_t n) {
    for (int it = 0; it < 100; ++it) {
        double off = 0.0;
        for (size_t i = 0; i < n; ++i) for (size_t j = i + 1; j < n; ++j) off += a[i * n + j] * a[i * n + j];
        if (off < 1e-300) break;
        for (size_t p = 0; p < n; ++p)
            for (size_t q = p + 1; q < n; ++q) {
                if (a[p * n + q] == 0.0) continue;
                const double th = (a[q * n + q] - a[p * n + p]) / (2.0 * a[p * n + q]);
                const double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (size_t k = 0; k < n; ++k) {
                    const double akp = a[k * n + p], akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for (size_t k = 0; k < n; ++k) {
                    const double apk = a[p * n + k], aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
            }
    }
    std::vector<double> e(n);
    for (size_t i = 0; i < n; ++i) e[i] = a[i * n + i];
    return e;
}
}  // namespace detalle_vent

// construir_proyector_aniquilador: M_X = I - X (X^T X)^-1 X^T sobre t = linspace(-1, 1, W), X = [t^0 .. t^orden]. W x W fila mayor.
inline std::vector<double> construir_proyector_aniquilador(int64_t W, int64_t orden = 1) {
    if (W < orden + 1)
        throw ErrorParametroInvalido("W=" + std::to_string(W) + " debe ser >= orden + 1=" + std::to_string(orden + 1) + " para construir proyector.");
    if (orden < 0) throw ErrorParametroInvalido("orden=" + std::to_string(orden) + " debe ser >= 0.");
    const size_t n = static_cast<size_t>(W), p = static_cast<size_t>(orden + 1);
    std::vector<double> t(n);
    if (n == 1) t[0] = -1.0;
    else {
        const double paso = 2.0 / static_cast<double>(n - 1);
        for (size_t k = 0; k < n; ++k) t[k] = static_cast<double>(k) * paso + -1.0;
        t[n - 1] = 1.0;
    }
    std::vector<double> X(n * p);
    for (size_t r = 0; r < n; ++r)
        for (size_t k = 0; k < p; ++k) X[r * p + k] = k == 0 ? 1.0 : (k == 1 ? t[r] : (k == 2 ? t[r] * t[r] : std::pow(t[r], static_cast<double>(k))));
    std::vector<double> XtX(p * p, 0.0);
    for (size_t a = 0; a < p; ++a) for (size_t b = 0; b < p; ++b) for (size_t r = 0; r < n; ++r) XtX[a * p + b] += X[r * p + a] * X[r * p + b];
    std::vector<double> inv;
    const auto ev = detalle_vent::autovalores_sim(XtX, p);
    double emax = 0.0, emin = std::numeric_limits<double>::infinity();
    for (double e : ev) { emax = std::max(emax, std::fabs(e)); emin = std::min(emin, std::fabs(e)); }
    if (!detalle_vent::invertir(XtX, p, inv)) throw ErrorParametroInvalido("XtX singular");
    std::vector<double> XI(n * p, 0.0);
    for (size_t r = 0; r < n; ++r) for (size_t b = 0; b < p; ++b) for (size_t a = 0; a < p; ++a) XI[r * p + b] += X[r * p + a] * inv[a * p + b];
    std::vector<double> Mx(n * n);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (size_t b = 0; b < p; ++b) s += XI[i * p + b] * X[j * p + b];
            Mx[i * n + j] = (i == j ? 1.0 : 0.0) - s;
        }
    (void)emax; (void)emin;
    return Mx;
}

// aplicar_detrending_ortogonal: M_X @ valores (1-D de largo W, o matriz W x M fila mayor con m > 0).
inline std::vector<double> aplicar_detrending_ortogonal(const std::vector<double>& v, int64_t W, int64_t m = 0, int64_t orden = 1) {
    const auto Mx = construir_proyector_aniquilador(W, orden);
    const size_t n = static_cast<size_t>(W), M = m == 0 ? 1 : static_cast<size_t>(m);
    std::vector<double> o(n * M, 0.0);
    for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < M; ++j) {
        double s = 0.0;
        for (size_t k = 0; k < n; ++k) s += Mx[i * n + k] * v[k * M + j];
        o[i * M + j] = s;
    }
    return o;
}

// segmentar_kantelhardt_bidireccional: (segmentos fila mayor, numero de segmentos, descarte).
struct Segmentacion {
    std::vector<double> segmentos;
    int64_t n_segmentos = 0, s = 0;
    double descarte = 0.0;
};
inline Segmentacion segmentar_kantelhardt_bidireccional(const std::vector<double>& arr, int64_t s, bool bidireccional = true) {
    const int64_t W = static_cast<int64_t>(arr.size());
    if (s < 2 || s > W) throw ErrorParametroInvalido("Escala s=" + std::to_string(s) + " inválida para W=" + std::to_string(W) + ".");
    const int64_t Ns = W / s;
    if (Ns == 0) throw ErrorDimensionInsuficiente("W=" + std::to_string(W) + " insuficiente para escala s=" + std::to_string(s) + ".");
    Segmentacion r;
    r.s = s;
    for (int64_t k = 0; k < Ns; ++k) r.segmentos.insert(r.segmentos.end(), arr.begin() + k * s, arr.begin() + (k + 1) * s);
    const int64_t res = W % s;
    if (res == 0 || !bidireccional) {
        r.descarte = res == 0 ? 0.0 : static_cast<double>(res) / static_cast<double>(W);
        r.n_segmentos = Ns;
        return r;
    }
    for (int64_t k = 0; k < Ns; ++k) r.segmentos.insert(r.segmentos.end(), arr.begin() + (W - (k + 1) * s), arr.begin() + (W - k * s));
    r.n_segmentos = 2 * Ns;
    r.descarte = 0.0;
    return r;
}

// alinear_bloques_welch_al_presente: (x[r:], r).
inline std::pair<std::vector<double>, int64_t> alinear_bloques_welch_al_presente(const std::vector<double>& x, int64_t nperseg, int64_t noverlap) {
    const int64_t N = static_cast<int64_t>(x.size());
    if (nperseg > N) throw ErrorDimensionInsuficiente("nperseg=" + std::to_string(nperseg) + " no puede ser mayor que N=" + std::to_string(N) + ".");
    if (noverlap >= nperseg) throw ErrorParametroInvalido("noverlap=" + std::to_string(noverlap) + " debe ser < nperseg=" + std::to_string(nperseg) + ".");
    const int64_t step = nperseg - noverlap;
    const int64_t nb = detalle_vent::div_piso(N - noverlap, step);
    if (nb < 1) throw ErrorDimensionInsuficiente("Datos insuficientes para al menos un bloque de Welch.");
    const int64_t util = (nb - 1) * step + nperseg;
    const int64_t r = N - util;
    const int64_t desde = r < 0 ? std::max<int64_t>(0, N + r) : std::min(r, N);   // arr[r:] de Python
    return {std::vector<double>(x.begin() + desde, x.end()), r};
}

// desacoplar_rampa_convolucional_tv: ADMM de deconvolucion_tv (rama importable) o el umbral suave de respaldo.
inline bool& tv_admm_disponible() { static bool b = true; return b; }
inline std::vector<double> desacoplar_rampa_convolucional_tv(const std::vector<double>& arr, double lambda_tv = 1.0,
                                                             std::optional<int64_t> ancho_kernel_W = std::nullopt) {
    const int64_t N = static_cast<int64_t>(arr.size());
    int64_t w_ef;
    if (ancho_kernel_W) w_ef = std::min(*ancho_kernel_W, std::max<int64_t>(2, N - 1));
    else w_ef = std::max<int64_t>(4, std::min<int64_t>(20, std::max<int64_t>(2, N / 4)));
    if (tv_admm_disponible()) {
        ConfiguracionTV c;
        c.ancho_kernel_W = w_ef;
        c.lambda_tv = lambda_tv;
        c.min_muestras = std::min<int64_t>(32, w_ef + 1);
        c.max_iter = 500;
        return resolver_deconvolucion_tv_admm(arr, c).serie_deconvolucionada;
    }
    std::vector<double> rec(arr.size());
    if (arr.empty()) return rec;
    const double um = lambda_tv / static_cast<double>(N);
    rec[0] = arr[0];
    double acum = 0.0;
    for (size_t k = 1; k < arr.size(); ++k) {
        const double d = arr[k] - arr[k - 1];
        const double sg = d > 0 ? 1.0 : (d < 0 ? -1.0 : 0.0);
        acum += sg * std::max(std::fabs(d) - um, 0.0);
        rec[k] = arr[0] + acum;
    }
    return rec;
}

// =================================================================================================== _LazoFeedbackShannon
class LazoFeedbackShannon {
public:
    int64_t W;
    std::optional<int64_t> ventana_fija;
    double tick;
    int64_t s_min;
    LazoOraculo oraculo;
    PropuestaLazo ultimo_p;
    bool hay_ultimo = false;

    LazoFeedbackShannon(int64_t W_, const ConfiguracionVentanaCausal& c)
        : W(W_), ventana_fija(c.ventana_innovacion_shannon),
          tick((c.tamano_tick && *c.tamano_tick > 0) ? *c.tamano_tick : 0.01), s_min(c.s_min), oraculo(W_, c.s_min, tick, 4096) {
        if (ventana_fija && *ventana_fija > 0) { oraculo.w_innov = *ventana_fija; oraculo.w_max = *ventana_fija; }
    }
    double paso(double val_k, bool mask_k) {
        const PropuestaLazo p = oraculo.calcular(val_k, mask_k);
        ultimo_p = p;
        hay_ultimo = true;
        oraculo.confirmar(val_k, mask_k, p);
        if (ventana_fija && *ventana_fija > 0) { oraculo.w_innov = *ventana_fija; oraculo.w_max = *ventana_fija; }
        return p.d_nu;
    }
};

// ===================================================================================================== conversiones numpy/pandas
namespace detalle_vent {

inline std::string repr_obj(const Obj& o, bool ustr) {
    if (o.t == Obj::NONE) return "None";
    if (o.t == Obj::FLOAT) return ingesta::repr_py(o.f);
    return ustr ? "np.str_(" + ingesta::repr_str_py(o.s) + ")" : ingesta::repr_str_py(o.s);
}
inline std::string str_obj(const Obj& o) {
    if (o.t == Obj::NONE) return "None";
    if (o.t == Obj::FLOAT) return ingesta::repr_py(o.f);
    return o.s;
}

// float(str) de Python (lo que hace np.asarray(objetos, float)).
inline double float_de_str(const std::string& s0) {
    const std::string s = strip(s0);
    std::string l = s;
    for (auto& c : l) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (l == "nan" || l == "+nan" || l == "-nan") return std::numeric_limits<double>::quiet_NaN();
    if (l == "inf" || l == "+inf" || l == "infinity" || l == "+infinity") return std::numeric_limits<double>::infinity();
    if (l == "-inf" || l == "-infinity") return -std::numeric_limits<double>::infinity();
    double v = 0.0;
    const char* b = s.data();
    const char* e = s.data() + s.size();
    if (b != e && *b == '+') ++b;
    const auto r = std::from_chars(b, e, v);
    if (s.empty() || r.ec != std::errc() || r.ptr != e) throw ValueError("could not convert string to float: " + ingesta::repr_str_py(s0));
    return v;
}

// np.asarray(x, dtype=float)
inline std::vector<double> a_f64(const Arreglo& a) {
    std::vector<double> o(a.size());
    switch (a.dt) {
        case Dt::F64: return a.f;
        case Dt::I64: case Dt::BOOL: for (size_t k = 0; k < o.size(); ++k) o[k] = static_cast<double>(a.i[k]); return o;
        default:
            for (size_t k = 0; k < o.size(); ++k) {
                const Obj& x = a.o[k];
                if (x.t == Obj::FLOAT) o[k] = x.f;
                else if (x.t == Obj::NONE) throw TypeError("float() argument must be a string or a real number, not 'NoneType'");
                else o[k] = float_de_str(x.s);
            }
            return o;
    }
}
// np.asarray(x, dtype=np.int64) sobre un ndarray (astype: float -> trunc; NaN/inf -> INT64_MIN).
inline int64_t f64_a_i64(double v) {
    if (!std::isfinite(v) || v >= 9.2233720368547758e18 || v < -9.2233720368547758e18) return std::numeric_limits<int64_t>::min();
    return static_cast<int64_t>(v);
}
inline std::vector<int64_t> a_i64(const Arreglo& a) {
    std::vector<int64_t> o(a.size());
    switch (a.dt) {
        case Dt::F64: for (size_t k = 0; k < o.size(); ++k) o[k] = f64_a_i64(a.f[k]); return o;
        case Dt::I64: case Dt::BOOL: return a.i;
        default:
            for (size_t k = 0; k < o.size(); ++k) {
                const Obj& x = a.o[k];
                if (x.t == Obj::FLOAT) o[k] = f64_a_i64(x.f);
                else throw ValueError("invalid literal for int() with base 10: " + repr_obj(x, false));
            }
            return o;
    }
}
// np.asarray(x, dtype=bool)
inline std::vector<uint8_t> a_bool(const Arreglo& a) {
    std::vector<uint8_t> o(a.size());
    switch (a.dt) {
        case Dt::F64: for (size_t k = 0; k < o.size(); ++k) o[k] = a.f[k] != 0.0; return o;   // NaN -> True
        case Dt::I64: case Dt::BOOL: for (size_t k = 0; k < o.size(); ++k) o[k] = a.i[k] != 0; return o;
        default:
            for (size_t k = 0; k < o.size(); ++k) {
                const Obj& x = a.o[k];
                o[k] = x.t == Obj::NONE ? 0 : (x.t == Obj::FLOAT ? (x.f != 0.0) : !x.s.empty());
            }
            return o;
    }
}

// Columnas de pandas -> numpy
inline std::vector<double> col_f64(const Columna& c) {
    std::vector<double> o(c.size());
    switch (c.tipo) {
        case TipoColumna::F64: return c.f;
        case TipoColumna::I64: case TipoColumna::BOOL: for (size_t k = 0; k < o.size(); ++k) o[k] = static_cast<double>(c.i[k]); return o;
        case TipoColumna::STR:
            for (size_t k = 0; k < o.size(); ++k) o[k] = c.es_nulo_str(k) ? std::numeric_limits<double>::quiet_NaN() : float_de_str(c.s[k]);
            return o;
        default: throw TypeError("Cannot cast DatetimeArray to dtype float64");
    }
}
inline std::vector<int64_t> col_i64(const Columna& c, int64_t unidad_ns) {
    std::vector<int64_t> o(c.size());
    switch (c.tipo) {
        case TipoColumna::F64:
            for (double v : c.f) if (std::isnan(v)) throw ValueError("cannot convert float NaN to integer");
            for (size_t k = 0; k < o.size(); ++k) o[k] = f64_a_i64(c.f[k]);
            return o;
        case TipoColumna::I64: case TipoColumna::BOOL: return c.i;
        case TipoColumna::TS:
            for (size_t k = 0; k < o.size(); ++k) o[k] = c.i[k] == std::numeric_limits<int64_t>::min() ? c.i[k] : c.i[k] / unidad_ns;
            return o;
        default: throw ValueError("invalid literal for int() with base 10");
    }
}
inline std::vector<uint8_t> col_bool(const Columna& c) {
    std::vector<uint8_t> o(c.size());
    switch (c.tipo) {
        case TipoColumna::F64: for (size_t k = 0; k < o.size(); ++k) o[k] = c.f[k] != 0.0; return o;
        case TipoColumna::I64: case TipoColumna::BOOL: case TipoColumna::TS: for (size_t k = 0; k < o.size(); ++k) o[k] = c.i[k] != 0; return o;
        default: for (size_t k = 0; k < o.size(); ++k) o[k] = c.es_nulo_str(k) ? 1 : !c.s[k].empty(); return o;   // bool(nan) = True
    }
}
// df[col].to_numpy() de una columna de texto: arreglo object con str y NaN.
inline Arreglo col_obj(const Columna& c) {
    Arreglo a;
    a.dt = Dt::OBJ;
    a.o.resize(c.size());
    for (size_t k = 0; k < c.size(); ++k) {
        Obj& x = a.o[k];
        if (c.tipo == TipoColumna::STR) {
            if (c.es_nulo_str(k)) { x.t = Obj::FLOAT; x.f = std::numeric_limits<double>::quiet_NaN(); }
            else { x.t = Obj::STR; x.s = c.s[k]; }
        } else if (c.tipo == TipoColumna::F64) {
            x.t = Obj::FLOAT; x.f = c.f[k];
        } else {
            x.t = Obj::FLOAT; x.f = static_cast<double>(c.i[k]);
        }
    }
    return a;
}
inline bool numerica(const Columna& c) { return c.tipo == TipoColumna::F64 || c.tipo == TipoColumna::I64; }   // select_dtypes(np.number)
inline const Columna* primera_numerica(const Tabla& t) {
    for (const auto& c : t.columnas) if (numerica(c)) return &c;
    return nullptr;
}
inline int64_t unidad(const Marco& m, const std::string& c) {
    const auto it = m.unidad_ts_ns.find(c);
    return it == m.unidad_ts_ns.end() ? 1 : it->second;
}
// D2: `fecha_utc` (TS) en segundos de epoca, si la tabla la trae.
inline std::optional<std::vector<int64_t>> segundos_fecha_utc(const Tabla& t) {
    for (const auto& c : t.columnas)
        if (c.nombre == "fecha_utc" && c.tipo == TipoColumna::TS) {
            std::vector<int64_t> v(c.i.size());
            for (size_t k = 0; k < v.size(); ++k) v[k] = div_piso(c.i[k], 1000000000LL);
            return v;
        }
    return std::nullopt;
}
// Sellos desde el indice (DatetimeIndex en segundos o indice entero); nullopt si no hay.
inline std::optional<std::vector<int64_t>> ts_de_indice(const Indice& ix, size_t n) {
    if (ix.tipo == TipoIndice::RANGO) {
        std::vector<int64_t> v(n);
        for (size_t k = 0; k < n; ++k) v[k] = ix.inicio + static_cast<int64_t>(k);
        return v;
    }
    if (ix.tipo == TipoIndice::ENTERO) return ix.valores;
    if (ix.tipo == TipoIndice::FECHA) {
        std::vector<int64_t> v(ix.valores.size());
        for (size_t k = 0; k < v.size(); ++k) v[k] = div_piso(ix.valores[k], 1000000000LL);   // D3: ns -> s, sin importar la unidad
        return v;
    }
    return std::nullopt;   // D5: OTRO o STR (texto) -> sin sellos
}

}  // namespace detalle_vent

// ============================================================================================================= FiltroVentanaEngine
class FiltroVentanaEngine {
public:
    ConfiguracionVentanaCausal config;
    std::optional<double> ultimo_lookahead_score;
    std::map<std::string, std::vector<double>> multiactivo_aux;   // vol / trades / noc, fila mayor (N, M)
    int64_t aux_N = 0, aux_M = 0;

    explicit FiltroVentanaEngine(ConfiguracionVentanaCausal c = {}) : config(std::move(c)) {
        config.validar();
        validar_metodo_imputacion_causal();
    }

    void validar_metodo_imputacion_causal() const {
        if (config.metodo_es_enum) return;
        const auto& vm = metodos_imputacion();
        if (std::find(vm.begin(), vm.end(), config.metodo_imputacion) == vm.end()) {
            std::set<std::string> s(vm.begin(), vm.end());
            throw ErrorViolacionCausalidad("VIOLACIÓN AXIOMA 1: El método de imputación '" + config.metodo_imputacion +
                                           "' viola la causalidad temporal o utiliza interpolación retrospectiva/biespectral. Solo se permite la lista "
                                           "blanca causal: " + detalle_vent::lista_ordenada(s) + ".");
        }
    }

    struct Columnas {   // los nombres opcionales de _stream_ventanas_core
        std::optional<std::string> datos, timestamps, mascara, volumen, trades;
    };

    // ------------------------------------------------------------------------------------------- _stream_ventanas_core
    void stream_ventanas_core(const std::vector<Trozo>& trozos, const Columnas& col, const std::optional<std::string>& venue,
                              const std::optional<std::string>& nivel, const std::optional<std::string>& simbolo, const Emisor& emitir) {
        using namespace detalle_vent;
        auto o_ = [](const std::optional<std::string>& a, const std::optional<std::string>& b) -> std::optional<std::string> {
            return (a && !a->empty()) ? a : ((b && !b->empty()) ? b : std::optional<std::string>{});
        };
        auto o_b = [](const std::optional<std::string>& a, const std::optional<std::string>& b) -> std::optional<std::string> {
            return (a && !a->empty()) ? a : b;   // `a or b` de Python
        };
        (void)o_;
        validar_compatibilidad_venue_nivel(o_b(venue, config.venue_esperado), o_b(nivel, config.nivel_esperado),
                                           o_b(simbolo, config.simbolo_esperado));
        if (config.venue_esperado && venue && *venue != *config.venue_esperado)
            throw ErrorContratoInvalido("VIOLACIÓN AXIOMA 2: Venue '" + *venue + "' no coincide con el contractualmente esperado '" +
                                        *config.venue_esperado + "'.");
        if (config.nivel_esperado && nivel && *nivel != *config.nivel_esperado)
            throw ErrorContratoInvalido("VIOLACIÓN AXIOMA 2: Nivel '" + *nivel + "' no coincide con el contractualmente esperado '" +
                                        *config.nivel_esperado + "'.");
        if (config.simbolo_esperado && simbolo && *simbolo != *config.simbolo_esperado)
            throw ErrorContratoInvalido("VIOLACIÓN AXIOMA 2: Símbolo '" + *simbolo + "' no coincide con el contractualmente esperado '" +
                                        *config.simbolo_esperado + "'.");

        const int64_t W = config.lookback_W, S = config.stride_S;
        const std::optional<int64_t> dtn = config.delta_t_barra_segundos;
        const int64_t offset_cierre = (config.estampar_al_cierre && dtn) ? *dtn : 0;
        const auto divisores = calcular_divisores_exactos(W, config.s_min, config.max_ratio_escalas);
        int64_t residuo_max = 0;
        double pct_descarte = 0.0;
        calc_residuo(divisores, residuo_max, pct_descarte);

        const size_t W2 = static_cast<size_t>(2 * W), Wz = static_cast<size_t>(W);
        std::vector<double> ring_datos(W2), ring_nu(W2), ring_dnu(W2);
        std::vector<uint8_t> ring_mask(W2);
        std::vector<int64_t> ring_gaps(W2, 0);
        LazoFeedbackShannon lazo(W, config);

        std::optional<std::string> venue_det = venue, nivel_det = nivel, simbolo_det = simbolo;
        if (venue_det && venue_det->empty()) {}   // se conserva: `venue_detectado = venue` (cadena vacia incluida)
        bool compat_ok = false;
        // v_prev_raw de Python: nullopt = None inicial
        std::optional<Obj> v_prev, n_prev, s_prev;

        int64_t cursor = 0, siguiente = W - 1, ts_prev = 0;
        double ultimo_valido = 0.0;
        bool tiene_previo = false;
        double nu_acum = 0.0;
        const TipoRelojVentana reloj = config.reloj;
        const bool congelar = config.congelar_en_huecos;
        const bool prohibir = config.metodo_imputacion == "PROHIBIR_INTERPOLACION";
        const auto& venue_esp = config.venue_esperado;
        const auto& nivel_esp = config.nivel_esperado;
        const auto& simbolo_esp = config.simbolo_esperado;

        for (const Trozo& tz : trozos) {
            std::optional<std::vector<double>> datos, vol, trades, noc;
            std::optional<std::vector<int64_t>> tsv;
            std::optional<std::vector<uint8_t>> maskv;
            std::optional<Arreglo> ven, niv, sim;
            std::vector<int64_t> forma_datos;
            bool datos_escalar = false;

            auto de_marco = [&](const Marco& mk, bool camino_dict, const Valor* ts_dict, const Valor* mask_dict) {
                const Tabla& t = mk.tabla;
                // datos
                if (col.datos && t.tiene(*col.datos)) datos = col_f64(t.col(*col.datos));
                else {
                    for (const char* c : {"c", "close", "mid_price", "precio", "val", "price"})
                        if (t.tiene(c)) { datos = col_f64(t.col(c)); break; }
                    if (!datos) {
                        const Columna* pc = primera_numerica(t);
                        if (pc) datos = col_f64(*pc);
                        else throw ErrorParametroInvalido("El DataFrame del chunk no contiene columnas numéricas válidas.");
                    }
                }
                // sellos
                const bool ts_ya = camino_dict && ts_dict && ts_dict->t == Valor::ARR;
                if (!ts_ya) {
                    if (!camino_dict && col.timestamps && t.tiene(*col.timestamps)) tsv = col_i64(t.col(*col.timestamps), unidad(mk, *col.timestamps));
                    else if (t.tiene("timestamp_utc")) tsv = col_i64(t.col("timestamp_utc"), unidad(mk, "timestamp_utc"));
                    else if (t.tiene("t")) tsv = col_i64(t.col("t"), unidad(mk, "t"));
                    else if (auto fs = segundos_fecha_utc(t)) tsv = fs;
                    else tsv = ts_de_indice(mk.indice, t.filas());
                }
                // mascara
                const bool m_ya = camino_dict && mask_dict && mask_dict->t == Valor::ARR;
                if (!m_ya) {
                    if (!camino_dict && col.mascara && t.tiene(*col.mascara)) maskv = col_bool(t.col(*col.mascara));
                    else if (t.tiene("es_imputado")) maskv = col_bool(t.col("es_imputado"));
                    else if (t.tiene("faltan")) maskv = col_bool(t.col("faltan"));
                }
                // volumen / transacciones / nocional
                if (!vol) {
                    if (!camino_dict && col.volumen && t.tiene(*col.volumen)) vol = col_f64(t.col(*col.volumen));
                    else for (const char* c : {"v", "volume", "vol", "volumen"}) if (t.tiene(c)) { vol = col_f64(t.col(c)); break; }
                }
                if (!trades) {
                    if (!camino_dict && col.trades && t.tiene(*col.trades)) trades = col_f64(t.col(*col.trades));
                    else for (const char* c : {"n_trades", "trades", "count", "n"}) if (t.tiene(c)) { trades = col_f64(t.col(c)); break; }
                }
                if (!noc) for (const char* c : {"notional", "q", "quote_volume", "nocional"}) if (t.tiene(c)) { noc = col_f64(t.col(c)); break; }
                if (!ven && t.tiene("venue")) ven = col_obj(t.col("venue"));
                if (!niv && t.tiene("nivel")) niv = col_obj(t.col("nivel"));
                if (!sim && t.tiene("simbolo")) sim = col_obj(t.col("simbolo"));
            };

            if (tz.t == Trozo::DICT) {
                const Valor* vd = nullptr;
                if (col.datos && !col.datos->empty() && tz.buscar(*col.datos)) vd = tz.buscar(*col.datos);
                else {
                    for (const char* k : {"datos", "close", "c", "precio", "val", "price"}) if (tz.buscar(k)) { vd = tz.buscar(k); break; }
                    if ((!vd || vd->t == Valor::NADA) && !tz.dict.empty()) vd = &tz.dict.front().second;
                }
                const Valor* vt = nullptr;
                if (col.timestamps && !col.timestamps->empty() && tz.buscar(*col.timestamps)) vt = tz.buscar(*col.timestamps);
                else for (const char* k : {"timestamps", "timestamp_utc", "t"}) if (tz.buscar(k)) { vt = tz.buscar(k); break; }
                const Valor* vm = nullptr;
                if (col.mascara && !col.mascara->empty() && tz.buscar(*col.mascara)) vm = tz.buscar(*col.mascara);
                else for (const char* k : {"mascara_huecos", "es_imputado", "faltan"}) if (tz.buscar(k)) { vm = tz.buscar(k); break; }
                auto arr_de = [](const Valor* v) -> std::optional<Arreglo> { if (v && v->t == Valor::ARR) return v->a; return std::nullopt; };
                const Valor* vv = (col.volumen && !col.volumen->empty() && tz.buscar(*col.volumen)) ? tz.buscar(*col.volumen)
                                  : (tz.buscar("v") ? tz.buscar("v") : tz.buscar("volume"));
                const Valor* vtr = (col.trades && !col.trades->empty() && tz.buscar(*col.trades)) ? tz.buscar(*col.trades)
                                   : (tz.buscar("trades") ? tz.buscar("trades") : tz.buscar("n_trades"));
                if (auto a = arr_de(vt)) tsv = detalle_vent::a_i64(*a);
                if (auto a = arr_de(vm)) maskv = detalle_vent::a_bool(*a);
                if (auto a = arr_de(vv)) vol = detalle_vent::a_f64(*a);
                if (auto a = arr_de(vtr)) trades = detalle_vent::a_f64(*a);
                if (auto a = arr_de(tz.buscar("venue"))) ven = *a;
                if (auto a = arr_de(tz.buscar("nivel"))) niv = *a;
                if (auto a = arr_de(tz.buscar("simbolo"))) sim = *a;
                if (vd && vd->t == Valor::MARCO) {
                    de_marco(vd->m, true, vt, vm);
                } else if (vd && vd->t == Valor::ARR) {
                    datos = detalle_vent::a_f64(vd->a);
                    forma_datos = vd->a.forma;
                    datos_escalar = vd->a.escalar;
                }
            } else if (tz.t == Trozo::MARCO) {
                de_marco(tz.m, false, nullptr, nullptr);
            } else if (tz.t == Trozo::SERIE) {
                datos = detalle_vent::a_f64(tz.a);
                tsv = ts_de_indice(tz.indice, tz.a.size());
            } else {
                datos = detalle_vent::a_f64(tz.a);
                forma_datos = tz.a.forma;
                datos_escalar = tz.a.escalar;
            }
            if (!datos) continue;
            if (forma_datos.size() > 1) {
                bool uno = false;
                for (auto d : forma_datos) if (d == 1) uno = true;
                if (!uno) {
                    std::string sh = "(";
                    for (size_t q = 0; q < forma_datos.size(); ++q) sh += (q ? ", " : "") + std::to_string(forma_datos[q]);
                    throw ErrorParametroInvalido("Los datos deben ser unidimensionales, shape=" + sh + ")");
                }
            }
            (void)datos_escalar;
            const std::vector<double>& D = *datos;
            const size_t Nk = D.size();
            if (Nk == 0) continue;
            std::vector<int64_t> TS;
            if (tsv) {
                if (tsv->size() != Nk)
                    throw ErrorParametroInvalido("Longitud de timestamps (" + std::to_string(tsv->size()) + ") no coincide con datos (" +
                                                 std::to_string(Nk) + ") en chunk.");
                TS = *tsv;
            } else {
                TS.resize(Nk);
                for (size_t k = 0; k < Nk; ++k) TS[k] = cursor + static_cast<int64_t>(k);
            }
            std::vector<uint8_t> MK;
            if (maskv) {
                if (maskv->size() != Nk)
                    throw ErrorParametroInvalido("Longitud de mascara_huecos (" + std::to_string(maskv->size()) + ") no coincide con datos (" +
                                                 std::to_string(Nk) + ") en chunk.");
                MK = *maskv;
            } else {
                MK.assign(Nk, 0);
            }
            if (reloj == TipoRelojVentana::NU_VOLUMEN) {
                if (!vol) throw ErrorRelojIntrinsecoInvalido("Reloj NU_VOLUMEN requiere columna de volumen.");
                if (vol->size() != Nk)
                    throw ErrorParametroInvalido("Longitud de volumen (" + std::to_string(vol->size()) + ") no coincide con datos (" + std::to_string(Nk) + ").");
            } else if (reloj == TipoRelojVentana::NU_TRADES) {
                if (!trades) throw ErrorRelojIntrinsecoInvalido("Reloj NU_TRADES requiere columna de trades.");
                if (trades->size() != Nk)
                    throw ErrorParametroInvalido("Longitud de trades (" + std::to_string(trades->size()) + ") no coincide con datos (" + std::to_string(Nk) + ").");
            } else if (reloj == TipoRelojVentana::NU_NOCIONAL) {
                if (noc) {
                } else if (vol) {
                    const size_t nv = vol->size();
                    if (nv != Nk && nv != 1 && Nk != 1)
                        throw ValueError("operands could not be broadcast together with shapes (" + std::to_string(Nk) + ",) (" + std::to_string(nv) + ",) ");
                    const size_t nn = std::max(nv, Nk);
                    std::vector<double> x(nn);
                    for (size_t k = 0; k < nn; ++k) x[k] = D[Nk == 1 ? 0 : k] * (*vol)[nv == 1 ? 0 : k];
                    noc = std::move(x);
                } else {
                    throw ErrorRelojIntrinsecoInvalido("Reloj NU_NOCIONAL requiere datos de precio y volumen.");
                }
                if (noc->size() != Nk)
                    throw ErrorParametroInvalido("Longitud de nocional (" + std::to_string(noc->size()) + ") no coincide con datos (" + std::to_string(Nk) + ").");
            }

            for (size_t idx = 0; idx < Nk; ++idx) {
                const int64_t kg = cursor + static_cast<int64_t>(idx);
                auto revisar = [&](const std::optional<Arreglo>& arr, std::optional<Obj>& prev, std::optional<std::string>& det,
                                   const std::optional<std::string>& esp, const char* campo, const char* Campo, const char* plural,
                                   const char* extra) {
                    if (!arr) return;
                    const Obj& x = arr->o.empty() ? Obj{} : arr->o[idx];
                    Obj xv = x;
                    if (arr->dt == Dt::F64) { xv.t = Obj::FLOAT; xv.f = arr->f[idx]; }
                    else if (arr->dt == Dt::I64 || arr->dt == Dt::BOOL) { xv.t = Obj::FLOAT; xv.f = static_cast<double>(arr->i[idx]); }
                    const bool ustr = arr->dt == Dt::USTR;
                    // `x is None or (x is not prev and x != prev)` (D6: el None siempre se revisa)
                    bool igual;
                    if (xv.t == Obj::NONE) igual = false;
                    else if (!prev) igual = false;
                    else if (prev->t == Obj::STR && xv.t == Obj::STR) igual = prev->s == xv.s;
                    else if (prev->t == Obj::FLOAT && xv.t == Obj::FLOAT) igual = prev->f == xv.f;
                    else igual = prev->t == Obj::NONE && xv.t == Obj::NONE;
                    if (igual) return;
                    const bool nulo = xv.t == Obj::NONE || (xv.t == Obj::FLOAT && std::isnan(xv.f)) || es_nulo_token(str_obj(xv));
                    if (nulo)
                        throw ErrorContratoInvalido(std::string("VIOLACIÓN AXIOMA 2: Valor de '") + campo + "' nulo, NaN o vacío detectado en índice " +
                                                    std::to_string(kg) + ": " + repr_obj(xv, ustr) + ".");
                    const std::string vs = strip(str_obj(xv));
                    if (!det) det = vs;
                    else if (vs != *det) {
                        if (cursor > 0 && idx == 0)
                            throw ErrorContratoInvalido(std::string("VIOLACIÓN AXIOMA 2 (Contaminación en Frontera de Chunks): El chunk actual presenta ") +
                                                        campo + "='" + vs + "', pero los chunks previos establecieron '" + *det + "'.");
                        throw ErrorContratoInvalido(std::string("VIOLACIÓN AXIOMA 2: Stream heterogéneo detectado en índice ") + std::to_string(kg) +
                                                    ". Múltiples " + plural + " en el mismo stream: '" + *det + "' vs '" + vs + "'. " + extra);
                    }
                    if (esp && vs != *esp)
                        throw ErrorContratoInvalido(std::string("VIOLACIÓN AXIOMA 2: ") + Campo + " '" + vs + "' en índice " + std::to_string(kg) +
                                                    " no coincide con el contractualmente esperado '" + *esp + "'.");
                    prev = xv;
                };
                revisar(ven, v_prev, venue_det, venue_esp, "venue", "Venue", "venues", "Prohibido mezclar venues en la misma ventana.");
                revisar(niv, n_prev, nivel_det, nivel_esp, "nivel", "Nivel", "niveles", "Prohibido mezclar niveles L1/L2/L3 en la misma ventana.");
                revisar(sim, s_prev, simbolo_det, simbolo_esp, "simbolo", "Símbolo", "símbolos", "Prohibido mezclar símbolos en la misma ventana.");
                if (!compat_ok) {
                    validar_compatibilidad_venue_nivel(o_b(venue_det, venue_esp), o_b(nivel_det, nivel_esp), o_b(simbolo_det, simbolo_esp));
                    compat_ok = true;
                }
                // 2. monotonia y resolucion de barra
                const int64_t ts_k = TS[idx] + offset_cierre;
                int64_t gap_k = 0;
                if (kg > 0) {
                    const int64_t d_t = ts_k - ts_prev;
                    if (d_t < 0) {
                        if (cursor > 0 && idx == 0)
                            throw ErrorMonotoniaTiempo("VIOLACIÓN AXIOMA 2: Inversión temporal en frontera de chunks: anterior=" + i2s(ts_prev) +
                                                       " > nuevo=" + i2s(ts_k) + " (delta=" + i2s(d_t) + ").");
                        throw ErrorMonotoniaTiempo("VIOLACIÓN AXIOMA 2 (Inversión Temporal en índice " + i2s(kg - 1) + "): t[" + i2s(kg - 1) + "]=" +
                                                   i2s(ts_prev) + " > t[" + i2s(kg) + "]=" + i2s(ts_k) + " (delta=" + i2s(d_t) + "). "
                                                   "Paquete desordenado o reloj anacrónico detectado. Prohibido reordenar pasivamente con sort_values.");
                    } else if (d_t == 0) {
                        if (cursor > 0 && idx == 0)
                            throw ErrorMonotoniaTiempo("VIOLACIÓN AXIOMA 2: Timestamp duplicado en frontera de chunks: anterior=" + i2s(ts_prev) +
                                                       " == nuevo=" + i2s(ts_k) + ".");
                        throw ErrorMonotoniaTiempo("VIOLACIÓN AXIOMA 2 (Timestamp Duplicado en índice " + i2s(kg - 1) + "): t[" + i2s(kg - 1) + "]=" +
                                                   i2s(ts_prev) + " == t[" + i2s(kg) + "]=" + i2s(ts_k) + ". "
                                                   "Colisión temporal de eventos simultáneos no resueltos por reloj intrínseco.");
                    }
                    if (dtn && reloj == TipoRelojVentana::PARED) {
                        const double jit = config.tolerancia_jitter_segundos;
                        if (static_cast<double>(d_t) < static_cast<double>(*dtn) - jit)
                            throw ErrorMonotoniaTiempo("VIOLACIÓN AXIOMA 2 (Colisión Sub-Barra en índice " + i2s(kg - 1) + "): delta_t=" + i2s(d_t) +
                                                       "s < delta_t_barra_segundos=" + i2s(*dtn) + "s (tol=" + f2s(jit) + "s) entre t[" + i2s(kg - 1) +
                                                       "]=" + i2s(ts_prev) + " y t[" + i2s(kg) + "]=" + i2s(ts_k) + ".");
                        else if (d_t > *dtn)
                            gap_k = div_piso(d_t - 1, *dtn);
                    }
                }
                ts_prev = ts_k;
                // 3. extraccion e imputacion causal
                const double val_raw = D[idx];
                bool mask_k = MK[idx] != 0;
                double val_k;
                if (!std::isfinite(val_raw)) {
                    if (prohibir)
                        throw ErrorViolacionCausalidad("VIOLACIÓN AXIOMA 1: Valor no finito (NaN/Inf) detectado en índice " + i2s(kg) +
                                                       " bajo política PROHIBIR_INTERPOLACION.");
                    if (kg == 0 || !tiene_previo)
                        throw ErrorViolacionCausalidad("VIOLACIÓN AXIOMA 1: Valor inicial no finito (NaN/Inf) en índice " + i2s(kg) +
                                                       " sin dato pasado disponible. Imputar este punto requeriría backward-fill (bfill) desde el "
                                                       "futuro, lo cual está prohibido.");
                    val_k = ultimo_valido;
                    mask_k = true;
                } else {
                    if (mask_k && prohibir)
                        throw ErrorViolacionCausalidad("VIOLACIÓN AXIOMA 1: Muestra imputada detectada en índice " + i2s(kg) +
                                                       " bajo política PROHIBIR_INTERPOLACION.");
                    val_k = val_raw;
                    ultimo_valido = val_k;
                    tiene_previo = true;
                }
                // 4. reloj intrinseco
                double d_nu;
                if (reloj == TipoRelojVentana::PARED) d_nu = (congelar && mask_k) ? 0.0 : 1.0;
                else if (reloj == TipoRelojVentana::NU_VOLUMEN) {
                    const double v = (*vol)[idx];
                    if (v < 0.0)
                        throw ErrorRelojIntrinsecoInvalido("VIOLACIÓN AXIOMA 4: Volumen negativo detectado en índice " + i2s(kg) + ": " + f2s(v) + " < 0.");
                    d_nu = (congelar && mask_k) ? 0.0 : v;
                } else if (reloj == TipoRelojVentana::NU_TRADES) {
                    const double v = (*trades)[idx];
                    if (v < 0.0)
                        throw ErrorRelojIntrinsecoInvalido("VIOLACIÓN AXIOMA 4: Conteo de trades negativo detectado en índice " + i2s(kg) + ": " + f2s(v) +
                                                           " < 0.");
                    d_nu = (congelar && mask_k) ? 0.0 : v;
                } else if (reloj == TipoRelojVentana::NU_NOCIONAL) {
                    const double v = (*noc)[idx];
                    if (v < 0.0)
                        throw ErrorRelojIntrinsecoInvalido("VIOLACIÓN AXIOMA 4: Nocional negativo detectado en índice " + i2s(kg) + ": " + f2s(v) + " < 0.");
                    d_nu = (congelar && mask_k) ? 0.0 : v;
                } else {
                    d_nu = lazo.paso(val_k, mask_k);
                }
                nu_acum += d_nu;
                // 5. buferes espejo
                const size_t pos = static_cast<size_t>(kg % W);
                ring_datos[pos] = ring_datos[pos + Wz] = val_k;
                ring_mask[pos] = ring_mask[pos + Wz] = mask_k;
                ring_nu[pos] = ring_nu[pos + Wz] = nu_acum;
                ring_dnu[pos] = ring_dnu[pos + Wz] = d_nu;
                ring_gaps[pos] = ring_gaps[pos + Wz] = gap_k;
                // 6. emision
                if (kg == siguiente) {
                    const size_t st = static_cast<size_t>((kg + 1) % W);
                    VentanaCausalCertificada v;
                    v.t_actual = ts_k;
                    v.w = W;
                    std::vector<double> sl(ring_datos.begin() + static_cast<std::ptrdiff_t>(st), ring_datos.begin() + static_cast<std::ptrdiff_t>(st + Wz));
                    v.mascara_huecos.assign(ring_mask.begin() + static_cast<std::ptrdiff_t>(st), ring_mask.begin() + static_cast<std::ptrdiff_t>(st + Wz));
                    v.reloj_acumulado = ring_nu[st + Wz - 1];
                    v.delta_reloj_ventana = suma_np(&ring_dnu[st], Wz);
                    int64_t n_imp = 0;
                    for (auto b : v.mascara_huecos) n_imp += b ? 1 : 0;
                    int64_t omit = 0;
                    for (size_t q = st + 1; q < st + Wz; ++q) omit += ring_gaps[q];
                    double pct;
                    if (omit > 0) {
                        pct = static_cast<double>(n_imp + omit) / static_cast<double>(W + omit);
                        v.advertencias.push_back("Gap temporal camuflado en ventana t=" + i2s(ts_k) + ": " + i2s(omit) +
                                                 " barras omitidas sin máscara (densidad efectiva de pérdida=" + pct1(pct * 100) + "%).");
                    } else {
                        pct = static_cast<double>(n_imp) / static_cast<double>(W);
                    }
                    if (pct > config.umbral_imputacion_maximo) {
                        if (config.abortar_si_excede_imputacion)
                            throw ErrorMascaraImputacionExcedida("VIOLACIÓN AXIOMA 2: Porcentaje de imputación/pérdida " + pct1(pct * 100) + "% en t=" +
                                                                 i2s(ts_k) + " supera el umbral máximo permitido (" +
                                                                 pct1(config.umbral_imputacion_maximo * 100) + "%).");
                        v.es_valida = false;
                        v.advertencias.push_back("Ventana en t=" + i2s(ts_k) + " invalidada: Porcentaje de imputación/pérdida " + pct1(pct * 100) +
                                                 "% excede el umbral tolerado (" + pct1(config.umbral_imputacion_maximo * 100) + "%).");
                    }
                    if (config.deconvolucion_tv_activa) {
                        v.valores = desacoplar_rampa_convolucional_tv(sl, config.lambda_tv);
                        double mx = 0.0;
                        bool primero = true;
                        for (size_t q = 0; q < sl.size(); ++q) {
                            const double d = std::fabs(v.valores[q] - sl[q]);
                            if (primero || d > mx || std::isnan(d)) { mx = d; primero = false; }
                        }
                        v.lookahead_score = mx;
                    } else {
                        v.valores = std::move(sl);
                        v.lookahead_score = ultimo_lookahead_score ? *ultimo_lookahead_score : 0.0;
                    }
                    v.divisores_exactos = divisores;
                    v.residuo_temporal_maximo = residuo_max;
                    v.porcentaje_descarte_causal = pct_descarte;
                    v.venue = o_b(o_b(venue_det, venue_esp), std::string("GENERICO")).value();
                    v.nivel = o_b(o_b(nivel_det, nivel_esp), std::string("L1")).value();
                    v.simbolo = o_b(o_b(simbolo_det, simbolo_esp), std::string("SYNTH")).value();
                    v.tipo_reloj = reloj;
                    v.barras_omitidas = omit;
                    emitir(v);
                    siguiente += S;
                }
            }
            cursor += static_cast<int64_t>(Nk);
        }
    }

    // --------------------------------------------------------------------------------------------- generar_ventanas
    void generar_ventanas(const Datos& d, const std::optional<Arreglo>& timestamps, const std::optional<Arreglo>& mascara,
                          const std::optional<std::string>& venue, const std::optional<std::string>& nivel,
                          const std::optional<std::string>& simbolo, const std::optional<std::string>& columna_valor, const Emisor& emitir) {
        validar_metodo_imputacion_causal();
        const int64_t N = static_cast<int64_t>(d.len());
        const int64_t W = config.lookback_W;
        if (N < W)
            throw ErrorDimensionInsuficiente("Longitud de datos insuficiente: N=" + std::to_string(N) + " < lookback_W=" + std::to_string(W) + ".");
        if (d.t != Datos::MARCO) {
            if (config.reloj == TipoRelojVentana::NU_VOLUMEN)
                throw ErrorRelojIntrinsecoInvalido("Reloj NU_VOLUMEN en generar_ventanas requiere DataFrame con columna de volumen.");
            if (config.reloj == TipoRelojVentana::NU_TRADES)
                throw ErrorRelojIntrinsecoInvalido("Reloj NU_TRADES en generar_ventanas requiere DataFrame con columna de trades.");
            if (config.reloj == TipoRelojVentana::NU_NOCIONAL)
                throw ErrorRelojIntrinsecoInvalido("Reloj NU_NOCIONAL en generar_ventanas requiere DataFrame con columnas de precio y volumen.");
        }
        Trozo tz;
        if (d.t == Datos::MARCO && !timestamps && !mascara) {
            tz.t = Trozo::MARCO;
            tz.m = d.m;
        } else {
            tz.t = Trozo::DICT;
            tz.dict.push_back({"datos", d.t == Datos::MARCO ? Valor::marco(d.m) : Valor::arr(d.a)});
            tz.dict.push_back({"timestamps", timestamps ? Valor::arr(*timestamps) : Valor::nada()});
            tz.dict.push_back({"mascara_huecos", mascara ? Valor::arr(*mascara) : Valor::nada()});
        }
        Columnas c;
        c.datos = columna_valor;
        stream_ventanas_core({tz}, c, venue, nivel, simbolo, emitir);
    }
    std::vector<VentanaCausalCertificada> lista_ventanas(const Datos& d, const std::optional<Arreglo>& ts = std::nullopt,
                                                         const std::optional<Arreglo>& mk = std::nullopt) {
        std::vector<VentanaCausalCertificada> out;
        generar_ventanas(d, ts, mk, std::nullopt, std::nullopt, std::nullopt, std::nullopt, [&](const VentanaCausalCertificada& v) { out.push_back(v); });
        return out;
    }

    // ----------------------------------------------------------------------------------- calcular_vector_reloj_intrinseco
    std::pair<std::vector<double>, std::vector<double>> calcular_vector_reloj_intrinseco(const Datos& d, const Arreglo& mascara,
                                                                                       const std::optional<Arreglo>& precios = std::nullopt,
                                                                                       const std::optional<Arreglo>& columna_volumen = std::nullopt,
                                                                                       const std::optional<Arreglo>& columna_trades = std::nullopt) {
        using namespace detalle_vent;
        const size_t N = mascara.size();
        const TipoRelojVentana reloj = config.reloj;
        const std::vector<uint8_t> mascara_b = a_bool(mascara);   // D4: np.asarray(mascara, dtype=bool) al entrar
        auto congelar = [&](std::vector<double>& dn) {
            if (!config.congelar_en_huecos) return;
            if (mascara_b.size() != dn.size())
                throw IndexError("boolean index did not match indexed array along axis 0; size of axis is " + std::to_string(dn.size()) +
                                 " but size of corresponding boolean axis is " + std::to_string(mascara_b.size()));
            for (size_t k = 0; k < dn.size(); ++k) if (mascara_b[k]) dn[k] = 0.0;
        };
        auto cumsum = [](const std::vector<double>& v) { std::vector<double> c(v.size()); double s = 0.0; for (size_t k = 0; k < v.size(); ++k) { s += v[k]; c[k] = s; } return c; };
        auto revisar_neg = [&](const std::vector<double>& dn, const char* que) {
            for (size_t k = 0; k < dn.size(); ++k)
                if (dn[k] < 0)
                    throw ErrorRelojIntrinsecoInvalido(std::string("VIOLACIÓN AXIOMA 4: ") + que + " negativo detectado en índice " + std::to_string(k) + ": " +
                                                       f2s(dn[k]) + " < 0.");
        };
        if (reloj == TipoRelojVentana::PARED) {
            std::vector<double> dn(N, 1.0);
            congelar(dn);
            return {dn, cumsum(dn)};
        }
        if (reloj == TipoRelojVentana::NU_VOLUMEN || reloj == TipoRelojVentana::NU_TRADES) {
            const bool vol = reloj == TipoRelojVentana::NU_VOLUMEN;
            std::vector<double> dn;
            if (d.t == Datos::MARCO) {
                const char* cand_v[] = {"v", "volume", "vol", "volumen"};
                const char* cand_t[] = {"n_trades", "trades", "count", "n"};
                const char* hallada = nullptr;
                for (const char* c : (vol ? cand_v : cand_t)) if (d.m.tabla.tiene(c)) { hallada = c; break; }
                if (!hallada)
                    throw ErrorRelojIntrinsecoInvalido(vol ? "El DataFrame no contiene columna de volumen para el reloj NU_VOLUMEN."
                                                           : "El DataFrame no contiene columna de conteo de transacciones para NU_TRADES.");
                dn = col_f64(d.m.tabla.col(hallada));
            } else if (vol ? columna_volumen.has_value() : columna_trades.has_value()) {
                dn = a_f64(vol ? *columna_volumen : *columna_trades);
            } else {
                throw ErrorRelojIntrinsecoInvalido(vol ? "Reloj NU_VOLUMEN requiere suministrar datos de volumen."
                                                       : "Reloj NU_TRADES requiere suministrar conteo de transacciones.");
            }
            revisar_neg(dn, vol ? "Volumen" : "Conteo de trades");
            congelar(dn);
            return {dn, cumsum(dn)};
        }
        if (reloj == TipoRelojVentana::NU_NOCIONAL) {
            std::vector<double> dn;
            if (d.t == Datos::MARCO) {
                const Tabla& t = d.m.tabla;
                const char* cn = nullptr;
                for (const char* c : {"notional", "q", "quote_volume", "nocional"}) if (t.tiene(c)) { cn = c; break; }
                if (cn) dn = col_f64(t.col(cn));
                else {
                    const char* cp = nullptr;
                    for (const char* c : {"c", "close", "mid_price"}) if (t.tiene(c)) { cp = c; break; }
                    const char* cv = nullptr;
                    for (const char* c : {"v", "volume", "vol"}) if (t.tiene(c)) { cv = c; break; }
                    if (!(cp && cv)) throw ErrorRelojIntrinsecoInvalido("No se encontraron columnas para calcular el nocional (P * V).");
                    const auto p = col_f64(t.col(cp)), v = col_f64(t.col(cv));
                    dn.resize(p.size());
                    for (size_t k = 0; k < p.size(); ++k) dn[k] = p[k] * v[k];
                }
            } else if (precios && columna_volumen) {
                const auto p = a_f64(*precios), v = a_f64(*columna_volumen);
                if (p.size() != v.size() && p.size() != 1 && v.size() != 1)
                    throw ValueError("operands could not be broadcast together with shapes (" + std::to_string(p.size()) + ",) (" + std::to_string(v.size()) + ",) ");
                const size_t nn = std::max(p.size(), v.size());
                dn.resize(nn);
                for (size_t k = 0; k < nn; ++k) dn[k] = p[p.size() == 1 ? 0 : k] * v[v.size() == 1 ? 0 : k];
            } else {
                throw ErrorRelojIntrinsecoInvalido("Reloj NU_NOCIONAL requiere datos de precio y volumen.");
            }
            revisar_neg(dn, "Nocional");
            congelar(dn);
            return {dn, cumsum(dn)};
        }
        // NU_INFO_SHANNON
        std::vector<double> pr;
        if (precios) pr = a_f64(*precios);
        else if (d.t == Datos::MARCO) {
            const char* cp = nullptr;
            for (const char* c : {"c", "close", "mid_price", "precio"}) if (d.m.tabla.tiene(c)) { cp = c; break; }
            if (!cp) throw ErrorRelojIntrinsecoInvalido("No se encontraron precios para calcular retornos de Shannon.");
            pr = col_f64(d.m.tabla.col(cp));
        } else {
            pr = a_f64(d.a);
        }
        LazoFeedbackShannon lazo(config.lookback_W, config);
        const auto& mb = mascara_b;
        std::vector<double> dn(N, 0.0);
        for (size_t k = 0; k < N; ++k) {
            if (k >= pr.size()) throw IndexError("index " + std::to_string(k) + " is out of bounds for axis 0 with size " + std::to_string(pr.size()));
            dn[k] = lazo.paso(pr[k], mb[k] != 0);
        }
        return {dn, cumsum(dn)};
    }

    // ---------------------------------------------------------------------------------------- auditar_congelamiento_reloj
    JV auditar_congelamiento_reloj(const std::vector<double>& nu, const std::vector<uint8_t>& mask) const {
        int64_t fugas = 0, tot = 0;
        for (size_t k = 1; k < nu.size() && k < mask.size(); ++k) {
            if (mask[k]) ++tot;
            if ((nu[k] - nu[k - 1] > 1e-12) && mask[k]) ++fugas;
        }
        JV r = JV::objeto();
        r.poner("total_huecos", JV::entero(tot));
        r.poner("fugas_en_huecos", JV::entero(fugas));
        r.poner("es_congelado_puro", JV::booleano(fugas == 0));
        r.poner("reloj_acumulado_total", JV::real(nu.empty() ? 0.0 : nu.back()));
        return r;
    }

    // ------------------------------------------------------------------------------ auditar_integridad_contrato_y_reloj
    struct CampoContrato {   // str, arreglo o None
        std::optional<std::string> s;
        std::optional<Arreglo> a;
    };
    JV auditar_integridad_contrato_y_reloj(const std::vector<int64_t>& ts, std::optional<int64_t> delta_nominal = std::nullopt,
                                           const std::optional<std::vector<uint8_t>>& mascara = std::nullopt, const CampoContrato& venue = {},
                                           const CampoContrato& nivel = {}, const CampoContrato& simbolo = {}) const {
        using namespace detalle_vent;
        JV r = JV::objeto();
        const int64_t N = static_cast<int64_t>(ts.size());
        if (N < 2) {
            r.poner("es_monotono", JV::booleano(true));
            r.poner("inversiones", JV::entero(0));
            r.poner("duplicados", JV::entero(0));
            r.poner("total_muestras", JV::entero(N));
        } else {
            int64_t inv = 0, dup = 0;
            for (int64_t k = 1; k < N; ++k) {
                const int64_t d = ts[static_cast<size_t>(k)] - ts[static_cast<size_t>(k - 1)];
                if (d < 0) ++inv;
                if (d == 0) ++dup;
            }
            r.poner("es_monotono", JV::booleano(inv == 0 && dup == 0));
            r.poner("inversiones", JV::entero(inv));
            r.poner("duplicados", JV::entero(dup));
            r.poner("total_muestras", JV::entero(N));
            r.poner("t_inicio", JV::entero(ts.front()));
            r.poner("t_fin", JV::entero(ts.back()));
            const std::optional<int64_t> dtn = (delta_nominal && *delta_nominal != 0) ? delta_nominal : config.delta_t_barra_segundos;
            if (dtn && *dtn > 0) {
                int64_t gaps = 0, col = 0, falt = 0, mx = std::numeric_limits<int64_t>::min();
                for (int64_t k = 1; k < N; ++k) {
                    const int64_t d = ts[static_cast<size_t>(k)] - ts[static_cast<size_t>(k - 1)];
                    if (d > *dtn) ++gaps;
                    if (d > 0 && d < *dtn) ++col;
                    falt += std::max<int64_t>(0, div_piso(d - 1, *dtn));
                    mx = std::max(mx, d);
                }
                r.poner("delta_nominal", JV::entero(*dtn));
                r.poner("gaps_detectados", JV::entero(gaps));
                r.poner("colisiones_sub_barra", JV::entero(col));
                r.poner("barras_faltantes_no_declaradas", JV::entero(falt));
                r.poner("max_gap", JV::entero(mx));
            }
        }
        if (mascara) {
            int64_t tot = 0, raf = 0, act = 0;
            for (auto m : *mascara) {
                if (m) { ++tot; ++act; if (act > raf) raf = act; }
                else act = 0;
            }
            r.poner("total_imputados", JV::entero(tot));
            r.poner("ratio_imputacion_global", JV::real(mascara->empty() ? 0.0 : static_cast<double>(tot) / static_cast<double>(mascara->size())));
            r.poner("rafaga_maxima_imputada", JV::entero(raf));
        }
        auto resolver = [&](const CampoContrato& c, const char* nombre_c, const std::optional<std::string>& esp) -> std::optional<std::string> {
            if (!c.s && !c.a) return esp;
            if (c.s) {
                const std::string n = nombre_c;
                validar_compatibilidad_venue_nivel(n == "venue" ? c.s : std::nullopt, n == "nivel" ? c.s : std::nullopt,
                                                   n == "simbolo" ? c.s : std::nullopt);
                if (esp && strip(*c.s) != *esp)
                    throw ErrorContratoInvalido(std::string("VIOLACIÓN AXIOMA 2: ") + nombre_c + " '" + *c.s + "' != esperado '" + *esp + "'.");
                return strip(*c.s);
            }
            const Arreglo& a = *c.a;
            std::vector<Obj> el(a.size());
            for (size_t k = 0; k < a.size(); ++k) {
                if (a.dt == Dt::F64) { el[k].t = Obj::FLOAT; el[k].f = a.f[k]; }
                else if (a.dt == Dt::I64 || a.dt == Dt::BOOL) { el[k].t = Obj::FLOAT; el[k].f = static_cast<double>(a.i[k]); }
                else el[k] = a.o[k];
            }
            for (size_t k = 0; k < el.size(); ++k) {
                const Obj& x = el[k];
                if (x.t == Obj::NONE || (x.t == Obj::FLOAT && std::isnan(x.f)) || es_nulo_token(str_obj(x)))
                    throw ErrorContratoInvalido(std::string("VIOLACIÓN AXIOMA 2: ") + nombre_c + " nulo o vacío en índice " + std::to_string(k) + ".");
            }
            std::vector<std::string> unicos;
            for (const auto& x : el) {
                const std::string v = strip(str_obj(x));
                if (std::find(unicos.begin(), unicos.end(), v) == unicos.end()) unicos.push_back(v);
            }
            if (unicos.size() > 1)
                throw ErrorContratoInvalido(std::string("VIOLACIÓN AXIOMA 2: Stream heterogéneo en '") + nombre_c + "': " + ingesta::repr_lista_str(unicos) + ".");
            const std::optional<std::string> vu = unicos.size() == 1 ? std::optional<std::string>(unicos[0]) : esp;
            if (esp && vu != esp)
                throw ErrorContratoInvalido(std::string("VIOLACIÓN AXIOMA 2: ") + nombre_c + " '" + (vu ? *vu : std::string("None")) + "' != esperado '" + *esp + "'.");
            return vu;
        };
        const auto v = resolver(venue, "venue", config.venue_esperado);
        const auto n = resolver(nivel, "nivel", config.nivel_esperado);
        const auto s = resolver(simbolo, "simbolo", config.simbolo_esperado);
        validar_compatibilidad_venue_nivel(v, n, s);
        r.poner("contrato_valido", JV::booleano(true));
        r.poner("venue", JV::cadena((v && !v->empty()) ? *v : "GENERICO"));
        r.poner("nivel", JV::cadena((n && !n->empty()) ? *n : "L1"));
        r.poner("simbolo", JV::cadena((s && !s->empty()) ? *s : "SYNTH"));
        return r;
    }

    // --------------------------------------------------------------------------------- sincronizar_series_multiactivo
    struct Sincronizado {
        std::vector<int64_t> ts;
        std::vector<double> X;       // (N, M) fila mayor
        std::vector<uint8_t> Mask, Open;
        int64_t N = 0, M = 0;
        std::vector<std::string> activos;
    };
    Sincronizado sincronizar_series_multiactivo(const SeriesDict& series, const std::optional<ArreglosDict>& ts_dict = std::nullopt,
                                                const std::optional<ArreglosDict>& mask_dict = std::nullopt) {
        using namespace detalle_vent;
        if (series.empty()) throw ErrorParametroInvalido("series_dict debe ser un diccionario no vacío con los activos a sincronizar.");
        std::vector<std::string> activos;
        for (const auto& kv : series) activos.push_back(kv.first);
        auto tupla = [](const std::vector<std::string>& v) {
            std::string r = "(";
            for (size_t k = 0; k < v.size(); ++k) r += (k ? ", " : "") + ingesta::repr_str_py(v[k]);
            return r + (v.size() == 1 ? ",)" : ")");
        };
        if (config.activos_esperados) {
            const std::set<std::string> a(activos.begin(), activos.end()), b(config.activos_esperados->begin(), config.activos_esperados->end());
            if (a != b)
                throw ErrorContratoInvalido("VIOLACIÓN AXIOMA 3: Los activos suministrados " + tupla(activos) + " no coinciden con los esperados " +
                                            tupla(*config.activos_esperados) + ".");
        }
        auto buscar = [](const std::optional<ArreglosDict>& d, const std::string& k) -> const Arreglo* {
            if (!d) return nullptr;
            for (const auto& kv : *d) if (kv.first == k) return &kv.second;
            return nullptr;
        };
        struct P { std::vector<int64_t> ts; std::vector<double> y; std::vector<uint8_t> mask; std::optional<std::vector<double>> vol, trades, noc; };
        std::vector<P> pd(series.size());
        for (size_t q = 0; q < series.size(); ++q) {
            const std::string& k = series[q].first;
            const ItemActivo& it = series[q].second;
            P& p = pd[q];
            if (it.t == ItemActivo::MARCO) {
                const Tabla& t = it.m.tabla;
                const char* cy = nullptr;
                for (const char* c : {"c", "close", "mid_price", "precio", "val", "price"}) if (t.tiene(c)) { cy = c; break; }
                if (cy) p.y = col_f64(t.col(cy));
                else {
                    const Columna* pc = primera_numerica(t);
                    if (!pc) throw ErrorParametroInvalido("Activo '" + k + "' no contiene columnas numéricas válidas.");
                    p.y = col_f64(*pc);
                }
                if (t.tiene("timestamp_utc")) p.ts = col_i64(t.col("timestamp_utc"), unidad(it.m, "timestamp_utc"));
                else if (t.tiene("t")) p.ts = col_i64(t.col("t"), unidad(it.m, "t"));
                else if (auto fs = segundos_fecha_utc(t)) p.ts = *fs;
                else if (auto ix = ts_de_indice(it.m.indice, t.filas())) p.ts = *ix;
                else if (const Arreglo* a = buscar(ts_dict, k)) p.ts = a_i64(*a);
                else throw ErrorParametroInvalido("No se pudieron resolver timestamps para el activo '" + k + "'.");
                if (t.tiene("es_imputado")) p.mask = col_bool(t.col("es_imputado"));
                else if (t.tiene("faltan")) p.mask = col_bool(t.col("faltan"));
                else if (const Arreglo* a = buscar(mask_dict, k)) p.mask = a_bool(*a);
                else p.mask.assign(p.y.size(), 0);
                for (const char* c : {"v", "volume", "vol", "volumen"}) if (t.tiene(c)) { p.vol = col_f64(t.col(c)); break; }
                for (const char* c : {"trades", "num_trades", "n", "trade_count"}) if (t.tiene(c)) { p.trades = col_f64(t.col(c)); break; }
                for (const char* c : {"nocional", "notional", "quote_volume", "q"}) if (t.tiene(c)) { p.noc = col_f64(t.col(c)); break; }
            } else if (it.t == ItemActivo::SERIE) {
                p.y = a_f64(it.a);
                if (auto ix = ts_de_indice(it.indice, it.a.size())) p.ts = *ix;
                else if (const Arreglo* a = buscar(ts_dict, k)) p.ts = a_i64(*a);
                else throw ErrorParametroInvalido("Faltan timestamps para la serie '" + k + "'.");
                if (const Arreglo* a = buscar(mask_dict, k)) p.mask = a_bool(*a);
                else p.mask.assign(p.y.size(), 0);
            } else {
                p.y = a_f64(it.a);
                if (const Arreglo* a = buscar(ts_dict, k)) p.ts = a_i64(*a);
                else throw ErrorParametroInvalido("Faltan timestamps para el array '" + k + "'.");
                if (const Arreglo* a = buscar(mask_dict, k)) p.mask = a_bool(*a);
                else p.mask.assign(p.y.size(), 0);
            }
            for (size_t j = 0; j + 1 < p.ts.size(); ++j) {
                const int64_t d = p.ts[j + 1] - p.ts[j];
                if (d <= 0) {
                    if (d < 0)
                        throw ErrorMonotoniaTiempo("VIOLACIÓN AXIOMA 2/3: Inversión temporal en activo '" + k + "' en índice " + std::to_string(j) + ": t[" +
                                                   std::to_string(j) + "]=" + std::to_string(p.ts[j]) + " > t[" + std::to_string(j + 1) + "]=" +
                                                   std::to_string(p.ts[j + 1]) + ".");
                    throw ErrorMonotoniaTiempo("VIOLACIÓN AXIOMA 2/3: Timestamp duplicado en activo '" + k + "' en índice " + std::to_string(j) + ": t[" +
                                               std::to_string(j) + "]=" + std::to_string(p.ts[j]) + ".");
                }
            }
        }
        const size_t M = activos.size();
        const int64_t Wl = config.lookback_W;
        auto ss = [](const std::vector<int64_t>& a, int64_t v) { return static_cast<size_t>(std::lower_bound(a.begin(), a.end(), v) - a.begin()); };
        auto fuera = [](size_t i, size_t n) { return IndexError("index " + std::to_string(i) + " is out of bounds for axis 0 with size " + std::to_string(n)); };
        Sincronizado S;
        S.activos = activos;
        S.M = static_cast<int64_t>(M);
        multiactivo_aux.clear();
        if (config.politica_sincronizacion == PoliticaSincronizacion::INTERSECCION_ESTRICTA) {
            std::set<int64_t> comun(pd[0].ts.begin(), pd[0].ts.end());
            for (size_t q = 1; q < M; ++q) {
                std::set<int64_t> b(pd[q].ts.begin(), pd[q].ts.end()), c;
                std::set_intersection(comun.begin(), comun.end(), b.begin(), b.end(), std::inserter(c, c.begin()));
                comun.swap(c);
            }
            if (static_cast<int64_t>(comun.size()) < Wl)
                throw ErrorSincronizacionInvalida("VIOLACIÓN AXIOMA 3: Intersección temporal estricta insuficiente: " + std::to_string(comun.size()) +
                                                  " muestras comunes < lookback_W=" + std::to_string(Wl) +
                                                  ". Posible desfase de huso horario (timezone mismatch) o feeds desalineados.");
            S.ts.assign(comun.begin(), comun.end());
            const size_t N = S.ts.size();
            S.N = static_cast<int64_t>(N);
            S.X.assign(N * M, 0.0);
            S.Mask.assign(N * M, 0);
            S.Open.assign(N * M, 1);
            for (size_t q = 0; q < M; ++q) {
                for (size_t r = 0; r < N; ++r) {
                    const size_t ix = ss(pd[q].ts, S.ts[r]);
                    if (ix >= pd[q].y.size()) throw fuera(ix, pd[q].y.size());
                    S.X[r * M + q] = pd[q].y[ix];
                    if (ix >= pd[q].mask.size()) throw fuera(ix, pd[q].mask.size());
                    S.Mask[r * M + q] = pd[q].mask[ix];
                }
                for (size_t r = 0; r < N; ++r)
                    if (!std::isfinite(S.X[r * M + q]))
                        throw ErrorSincronizacionInvalida("VIOLACIÓN AXIOMA 3: Valor NaN o no finito detectado en activo '" + activos[q] +
                                                          "' bajo política INTERSECCION_ESTRICTA. Prohibido emitir ventanas con NaNs.");
            }
            aux_N = S.N; aux_M = S.M;
            auto llenar = [&](const char* clave, auto obtener) {
                std::vector<double> mat(N * M, 0.0);
                for (size_t q = 0; q < M; ++q)
                    for (size_t r = 0; r < N; ++r) mat[r * M + q] = obtener(q, ss(pd[q].ts, S.ts[r]), r);
                multiactivo_aux[clave] = std::move(mat);
            };
            if (config.reloj == TipoRelojVentana::NU_VOLUMEN) {
                for (size_t q = 0; q < M; ++q)
                    if (!pd[q].vol) throw ErrorRelojIntrinsecoInvalido("VIOLACIÓN AXIOMA 4: Activo '" + activos[q] + "' no contiene columna de volumen para reloj NU_VOLUMEN.");
                llenar("vol", [&](size_t q, size_t ix, size_t) { if (ix >= pd[q].vol->size()) throw fuera(ix, pd[q].vol->size()); return (*pd[q].vol)[ix]; });
            } else if (config.reloj == TipoRelojVentana::NU_TRADES) {
                for (size_t q = 0; q < M; ++q)
                    if (!pd[q].trades) throw ErrorRelojIntrinsecoInvalido("VIOLACIÓN AXIOMA 4: Activo '" + activos[q] + "' no contiene columna de trades para reloj NU_TRADES.");
                llenar("trades", [&](size_t q, size_t ix, size_t) { if (ix >= pd[q].trades->size()) throw fuera(ix, pd[q].trades->size()); return (*pd[q].trades)[ix]; });
            } else if (config.reloj == TipoRelojVentana::NU_NOCIONAL) {
                for (size_t q = 0; q < M; ++q)
                    if (!pd[q].vol && !pd[q].noc)
                        throw ErrorRelojIntrinsecoInvalido("VIOLACIÓN AXIOMA 4: Activo '" + activos[q] + "' requiere datos de precio y volumen para reloj NU_NOCIONAL.");
                llenar("noc", [&](size_t q, size_t ix, size_t r) {
                    if (pd[q].noc) { if (ix >= pd[q].noc->size()) throw fuera(ix, pd[q].noc->size()); return (*pd[q].noc)[ix]; }
                    if (ix >= pd[q].vol->size()) throw fuera(ix, pd[q].vol->size());
                    return S.X[r * M + q] * (*pd[q].vol)[ix];
                });
            }
            return S;
        }
        // UNION_CON_MASCARA
        std::set<int64_t> uni;
        for (const auto& p : pd) uni.insert(p.ts.begin(), p.ts.end());
        std::vector<int64_t> tu(uni.begin(), uni.end());
        const size_t Nu = tu.size();
        const double nanv = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> X(Nu * M, nanv);
        std::vector<uint8_t> Mk(Nu * M, 1), Op(Nu * M, 0);
        for (size_t q = 0; q < M; ++q) {
            if (pd[q].y.size() != pd[q].ts.size() || pd[q].mask.size() != pd[q].ts.size()) {
                const size_t a = pd[q].ts.size(), b = pd[q].y.size() != a ? pd[q].y.size() : pd[q].mask.size();
                throw ValueError("shape mismatch: value array of shape (" + std::to_string(b) + ",) could not be broadcast to indexing result of shape (" +
                                 std::to_string(a) + ",)");
            }
            for (size_t j = 0; j < pd[q].ts.size(); ++j) {
                const size_t ix = ss(tu, pd[q].ts[j]);
                X[ix * M + q] = pd[q].y[j];
                Mk[ix * M + q] = pd[q].mask[j];
                Op[ix * M + q] = 1;
            }
        }
        size_t ini = 0;
        for (size_t q = 0; q < M; ++q) {
            size_t pri = Nu;
            for (size_t r = 0; r < Nu; ++r) if (!std::isnan(X[r * M + q])) { pri = r; break; }
            if (pri == Nu) throw ErrorSincronizacionInvalida("Activo '" + activos[q] + "' no tiene ninguna observación válida.");
            ini = std::max(ini, pri);
        }
        tu.erase(tu.begin(), tu.begin() + static_cast<std::ptrdiff_t>(ini));
        X.erase(X.begin(), X.begin() + static_cast<std::ptrdiff_t>(ini * M));
        Mk.erase(Mk.begin(), Mk.begin() + static_cast<std::ptrdiff_t>(ini * M));
        Op.erase(Op.begin(), Op.begin() + static_cast<std::ptrdiff_t>(ini * M));
        const size_t Nr = tu.size();
        if (static_cast<int64_t>(Nr) < Wl)
            throw ErrorSincronizacionInvalida("VIOLACIÓN AXIOMA 3: Longitud post-inicio conjunto (" + std::to_string(Nr) + ") < lookback_W=" + std::to_string(Wl) + ".");
        for (size_t q = 0; q < M; ++q) {
            double last = X[q];
            for (size_t r = 0; r < Nr; ++r) {
                if (std::isnan(X[r * M + q])) {
                    X[r * M + q] = last;
                    Op[r * M + q] = 0;
                    Mk[r * M + q] = 1;
                } else {
                    last = X[r * M + q];
                }
            }
        }
        aux_N = static_cast<int64_t>(Nr); aux_M = S.M;
        auto llenar_u = [&](const char* clave, auto valor) {
            std::vector<double> mat(Nr * M, 0.0);
            for (size_t q = 0; q < M; ++q)
                for (size_t j = 0; j < pd[q].ts.size(); ++j) {
                    const size_t ix = ss(tu, pd[q].ts[j]);
                    const bool ok = ix < Nr && tu[std::min(ix, Nr - 1)] == pd[q].ts[j];
                    if (ok) mat[ix * M + q] = valor(q, j);
                }
            multiactivo_aux[clave] = std::move(mat);
        };
        if (config.reloj == TipoRelojVentana::NU_VOLUMEN) {
            for (size_t q = 0; q < M; ++q)
                if (!pd[q].vol) throw ErrorRelojIntrinsecoInvalido("VIOLACIÓN AXIOMA 4: Activo '" + activos[q] + "' no contiene columna de volumen para reloj NU_VOLUMEN.");
            llenar_u("vol", [&](size_t q, size_t j) { return (*pd[q].vol)[j]; });
        } else if (config.reloj == TipoRelojVentana::NU_TRADES) {
            for (size_t q = 0; q < M; ++q)
                if (!pd[q].trades) throw ErrorRelojIntrinsecoInvalido("VIOLACIÓN AXIOMA 4: Activo '" + activos[q] + "' no contiene columna de trades para reloj NU_TRADES.");
            llenar_u("trades", [&](size_t q, size_t j) { return (*pd[q].trades)[j]; });
        } else if (config.reloj == TipoRelojVentana::NU_NOCIONAL) {
            for (size_t q = 0; q < M; ++q)
                if (!pd[q].vol && !pd[q].noc)
                    throw ErrorRelojIntrinsecoInvalido("VIOLACIÓN AXIOMA 4: Activo '" + activos[q] + "' requiere datos de precio y volumen para reloj NU_NOCIONAL.");
            llenar_u("noc", [&](size_t q, size_t j) { return pd[q].noc ? (*pd[q].noc)[j] : pd[q].y[j] * (*pd[q].vol)[j]; });
        }
        S.ts = std::move(tu);
        S.N = static_cast<int64_t>(Nr);
        S.X = std::move(X);
        S.Mask = std::move(Mk);
        S.Open = std::move(Op);
        return S;
    }

    // ------------------------------------------------------------------------------------ generar_ventanas_multiactivo
    void generar_ventanas_multiactivo(const SeriesDict& series, const std::optional<ArreglosDict>& ts_dict,
                                      const std::optional<ArreglosDict>& mask_dict, const Emisor& emitir) {
        using namespace detalle_vent;
        const Sincronizado S = sincronizar_series_multiactivo(series, ts_dict, mask_dict);
        const size_t N = static_cast<size_t>(S.N), M = static_cast<size_t>(S.M);
        const int64_t W = config.lookback_W, St = config.stride_S;
        const size_t Wz = static_cast<size_t>(W);
        const TipoRelojVentana reloj = config.reloj;
        const bool congelar = config.congelar_en_huecos;
        const auto divisores = calcular_divisores_exactos(W, config.s_min, config.max_ratio_escalas);
        int64_t residuo_max = 0;
        double pct_descarte = 0.0;
        calc_residuo(divisores, residuo_max, pct_descarte);
        std::vector<double> dnu(N * M, 0.0);
        if (reloj == TipoRelojVentana::PARED) {
            for (size_t m = 0; m < M; ++m)
                for (size_t r = 0; r < N; ++r)
                    dnu[r * M + m] = !S.Open[r * M + m] ? 0.0 : ((congelar && S.Mask[r * M + m]) ? 0.0 : 1.0);
        } else if (reloj == TipoRelojVentana::NU_INFO_SHANNON) {
            std::vector<LazoFeedbackShannon> lazos;
            lazos.reserve(M);
            for (size_t m = 0; m < M; ++m) lazos.emplace_back(W, config);
            for (size_t m = 0; m < M; ++m)
                for (size_t r = 0; r < N; ++r) {
                    if (!S.Open[r * M + m]) dnu[r * M + m] = 0.0;
                    else if (congelar && S.Mask[r * M + m]) dnu[r * M + m] = lazos[m].paso(S.X[r * M + m], true);
                    else dnu[r * M + m] = lazos[m].paso(S.X[r * M + m], S.Mask[r * M + m] != 0);
                }
        } else {
            const char* clave = reloj == TipoRelojVentana::NU_VOLUMEN ? "vol" : (reloj == TipoRelojVentana::NU_TRADES ? "trades" : "noc");
            const auto it = multiactivo_aux.find(clave);
            if (it == multiactivo_aux.end())
                throw ErrorRelojIntrinsecoInvalido(reloj == TipoRelojVentana::NU_VOLUMEN
                                                       ? "Reloj NU_VOLUMEN en multiactivo requiere suministrar datos de volumen para todos los activos."
                                                       : (reloj == TipoRelojVentana::NU_TRADES ? "Reloj NU_TRADES en multiactivo requiere datos de conteo de transacciones."
                                                                                               : "Reloj NU_NOCIONAL en multiactivo requiere datos de precio y volumen."));
            const auto& mat = it->second;
            for (size_t m = 0; m < M; ++m)
                for (size_t r = 0; r < N; ++r) {
                    const double v = mat[r * M + m];
                    dnu[r * M + m] = (!S.Open[r * M + m] || (congelar && S.Mask[r * M + m])) ? 0.0 : (v > 0.0 ? v : 0.0);   // max(0.0, v)
                }
        }
        std::vector<double> nua(N * M);
        for (size_t m = 0; m < M; ++m) { double s = 0.0; for (size_t r = 0; r < N; ++r) { s += dnu[r * M + m]; nua[r * M + m] = s; } }
        std::map<std::string, std::string> cache;
        for (int64_t i = W - 1; i < static_cast<int64_t>(N); i += St) {
            const size_t a = static_cast<size_t>(i - W + 1), b = static_cast<size_t>(i + 1);
            VentanaCausalCertificada v;
            v.t_actual = S.ts[static_cast<size_t>(i)];
            v.w = W;
            v.m = static_cast<int64_t>(M);
            v.valores.assign(S.X.begin() + static_cast<std::ptrdiff_t>(a * M), S.X.begin() + static_cast<std::ptrdiff_t>(b * M));
            v.mascara_huecos.assign(S.Mask.begin() + static_cast<std::ptrdiff_t>(a * M), S.Mask.begin() + static_cast<std::ptrdiff_t>(b * M));
            std::vector<uint8_t> op(S.Open.begin() + static_cast<std::ptrdiff_t>(a * M), S.Open.begin() + static_cast<std::ptrdiff_t>(b * M));
            int64_t nimp = 0;
            for (auto x : v.mascara_huecos) nimp += x ? 1 : 0;
            const double pct = static_cast<double>(nimp) / static_cast<double>(W * static_cast<int64_t>(M));
            if (pct > config.umbral_imputacion_maximo) {
                if (config.abortar_si_excede_imputacion)
                    throw ErrorMascaraImputacionExcedida("VIOLACIÓN AXIOMA 3: Imputación multiactivo " + pct1(pct * 100) + "% > " +
                                                         pct1(config.umbral_imputacion_maximo * 100) + "% en t=" + i2s(v.t_actual) + ".");
                v.es_valida = false;
                v.advertencias.push_back("Ventana multiactivo en t=" + i2s(v.t_actual) + " invalidada: Imputación " + pct1(pct * 100) + "% > " +
                                         pct1(config.umbral_imputacion_maximo * 100) + "%.");
            }
            std::vector<double> dpa(M);
            for (size_t m = 0; m < M; ++m) {
                if (M == 1) dpa[m] = suma_np(&dnu[a], Wz);
                else { double s = 0.0; for (size_t r = a; r < b; ++r) s += dnu[r * M + m]; dpa[m] = s; }
            }
            v.delta_reloj_ventana = suma_np(dpa.data(), M) / static_cast<double>(M);
            v.reloj_acumulado = suma_np(&nua[static_cast<size_t>(i) * M], M) / static_cast<double>(M);
            bool cierre = false;
            for (auto x : op) if (!x) cierre = true;
            if (cierre) {
                std::vector<std::string> ac;
                for (size_t m = 0; m < M; ++m) {
                    bool c = false;
                    for (size_t r = 0; r < Wz; ++r) if (!op[r * M + m]) c = true;
                    if (c) ac.push_back(S.activos[m]);
                }
                v.advertencias.push_back("Asimetría de tiempo propio en UNION_CON_MASCARA: " + std::to_string(ac.size()) +
                                         " activos presentan sesiones inactivas/cerradas en la ventana: " + ingesta::repr_lista_str(ac) + ".");
            }
            const uint32_t fecha = resolver_fecha_sesion(v.t_actual);
            const std::string fs = texto_yyyymmdd(fecha);
            for (size_t u = 0; u < M; ++u)
                for (size_t w = 0; w < M; ++w) {
                    if (S.activos[u] == S.activos[w]) continue;
                    const bool uo = op[(Wz - 1) * M + u] != 0, vo = op[(Wz - 1) * M + w] != 0;
                    std::string clase;
                    if (!vo) clase = "SIN_SESION";
                    else if (!uo) clase = "REZAGO";
                    else {
                        const std::string clave = S.activos[u] + "\x1f" + fs + "\x1f" + S.activos[w] + "\x1f" + fs;
                        auto itc = cache.find(clave);
                        if (itc == cache.end()) {
                            std::string c;
                            try {
                                c = nombre_clase(ark::dataset::clasificar_par(S.activos[u].c_str(), fecha, S.activos[w].c_str(), fecha).clase);
                            } catch (...) {
                                c = "NO_VERIFICADO";
                            }
                            itc = cache.emplace(clave, c).first;
                        }
                        clase = itc->second;
                    }
                    v.clasificacion_pares.push_back({{S.activos[u], S.activos[w]}, clase});
                    if (clase == "REZAGO")
                        v.advertencias.push_back("Guarda de asincronía causal: Par (" + S.activos[u] + " -> " + S.activos[w] +
                                                 ") clasificado en REZAGO temporal en t=" + i2s(v.t_actual) + ".");
                }
            v.market_open_mask = std::move(op);
            v.delta_reloj_por_activo = std::move(dpa);
            v.lookahead_score = ultimo_lookahead_score ? *ultimo_lookahead_score : 0.0;
            v.divisores_exactos = divisores;
            v.residuo_temporal_maximo = residuo_max;
            v.porcentaje_descarte_causal = pct_descarte;
            std::set<std::string> distintos(S.activos.begin(), S.activos.end());
            v.venue = distintos.size() > 1 ? "MULTIVENUE" : S.activos[0];
            v.nivel = "MULTINIVEL";
            std::string sim;
            for (size_t m = 0; m < M; ++m) sim += (m ? "," : "") + S.activos[m];
            v.simbolo = sim;
            v.activos = S.activos;
            v.tipo_reloj = reloj;
            v.barras_omitidas = 0;
            emitir(v);
        }
    }

    // ----------------------------------------------------------------------------------- auditar_testigo_anti_lookahead
    struct ReporteTestigo {
        double max_diff = 0.0;
        int64_t cortes_evaluados = 0;
        std::vector<int64_t> lista_cortes;
        int64_t corte_peor_diff = 0;
        double max_diff_perturbacion = 0.0, max_diff_truncado = 0.0;
        bool es_causal_estricto = true;
    };
    ReporteTestigo auditar_testigo_anti_lookahead(const Datos& datos, const std::optional<std::vector<int64_t>>& t_eval = std::nullopt,
                                                  bool t_eval_es_lista = false, const std::optional<Arreglo>& timestamps = std::nullopt,
                                                  const std::optional<Arreglo>& mascara = std::nullopt, double magnitud_shock = 1e12,
                                                  bool modo_exhaustivo = true, int64_t max_cortes_muestreo = 50, bool verificar_truncado = true,
                                                  bool perturbacion_aleatoria = true) {
        const int64_t N = static_cast<int64_t>(datos.len());
        const int64_t W = config.lookback_W, S = config.stride_S;
        if (N < W + 1)
            throw ErrorDimensionInsuficiente("Longitud N=" + std::to_string(N) + " insuficiente para auditar testigo anti-lookahead con W=" + std::to_string(W) +
                                             " (se requiere N >= W + 1).");
        std::vector<int64_t> ind;
        for (int64_t i = W - 1; i < N; i += S) if (i < N - 1) ind.push_back(i);
        if (ind.empty())
            throw ErrorParametroInvalido("No existe ningún paso de emisión con muestras futuras disponibles para N=" + std::to_string(N) + ", W=" +
                                         std::to_string(W) + ", S=" + std::to_string(S) + ".");
        auto mas_cercano = [&](int64_t t) {
            size_t k = 0;
            int64_t mejor = std::llabs(ind[0] - t);
            for (size_t q = 1; q < ind.size(); ++q) if (std::llabs(ind[q] - t) < mejor) { mejor = std::llabs(ind[q] - t); k = q; }
            return ind[k];
        };
        std::vector<int64_t> cortes;
        if (t_eval && !t_eval_es_lista) {
            const int64_t t = (*t_eval)[0];
            if (t < W - 1 || t >= N - 1)
                throw ErrorParametroInvalido("t_eval_indice=" + std::to_string(t) + " debe ser un índice posicional en [" + std::to_string(W - 1) + ", " +
                                             std::to_string(N - 2) + "] para evaluar la ventana y disponer de al menos una muestra futura (idx_pos + 1 < N).");
            cortes = {mas_cercano(t)};
        } else if (t_eval) {
            std::set<int64_t> c;
            for (int64_t t : *t_eval) {
                if (t < W - 1 || t >= N - 1)
                    throw ErrorParametroInvalido("Índice de corte " + std::to_string(t) + " fuera de rango [" + std::to_string(W - 1) + ", " + std::to_string(N - 2) + "].");
                c.insert(mas_cercano(t));
            }
            cortes.assign(c.begin(), c.end());
        } else {
            if (modo_exhaustivo || static_cast<int64_t>(ind.size()) <= max_cortes_muestreo) cortes = ind;
            else {
                const double fr[] = {0.0, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0};
                std::set<int64_t> c;
                for (double f : fr) c.insert(ind[static_cast<size_t>(f * static_cast<double>(ind.size() - 1))]);
                const int64_t paso = std::max<int64_t>(1, static_cast<int64_t>(ind.size()) / max_cortes_muestreo);
                for (size_t q = 0; q < ind.size(); q += static_cast<size_t>(paso)) c.insert(ind[q]);
                cortes.assign(c.begin(), c.end());
                if (static_cast<int64_t>(cortes.size()) > max_cortes_muestreo) cortes.resize(static_cast<size_t>(max_cortes_muestreo));
            }
        }
        const auto orig = lista_ventanas(datos, timestamps, mascara);
        if (orig.size() < ind.size()) throw IndexError("list index out of range");
        double mx_g = 0.0, mx_p = 0.0, mx_t = 0.0;
        int64_t peor = cortes[0];
        numpy_rng::Generador rng(12345);
        auto maxabs = [](const std::vector<double>& a, const std::vector<double>& b) {
            double m = 0.0;
            bool p = true;
            for (size_t k = 0; k < a.size() && k < b.size(); ++k) {
                const double d = std::fabs(a[k] - b[k]);
                if (p || d > m || std::isnan(d)) { m = d; p = false; }
            }
            return m;
        };
        auto pymax = [](std::initializer_list<double> xs) { double m = *xs.begin(); for (double x : xs) if (x > m) m = x; return m; };
        for (int64_t idx : cortes) {
            const size_t ko = static_cast<size_t>(std::find(ind.begin(), ind.end(), idx) - ind.begin());
            const VentanaCausalCertificada& vo = orig[ko];
            Datos dp;
            if (datos.t == Datos::MARCO) {
                dp = datos;
                const size_t n_fut = static_cast<size_t>(N - (idx + 1));
                for (auto& c : dp.m.tabla.columnas) {
                    if (!detalle_vent::numerica(c)) continue;
                    if (c.nombre == "timestamp_utc" || c.nombre == "t" || c.nombre == "es_imputado" || c.nombre == "faltan") continue;
                    if (c.tipo == TipoColumna::I64) {   // D1: la copia pasa la columna entera a float64 antes del choque
                        c.f.resize(c.i.size());
                        for (size_t q = 0; q < c.i.size(); ++q) c.f[q] = static_cast<double>(c.i[q]);
                        c.i.clear();
                        c.tipo = TipoColumna::F64;
                    }
                    std::vector<double> ruido(n_fut, 0.0);
                    if (perturbacion_aleatoria) for (auto& r : ruido) r = -1e6 + (1e6 - -1e6) * rng.next_double();
                    std::vector<double> nv(n_fut);
                    for (size_t k = 0; k < n_fut; ++k) {
                        const size_t p = static_cast<size_t>(idx + 1) + k;
                        const double base = c.tipo == TipoColumna::F64 ? c.f[p] : static_cast<double>(c.i[p]);
                        nv[k] = (base + magnitud_shock) + ruido[k];
                    }
                    for (size_t k = 0; k < n_fut; ++k) c.f[static_cast<size_t>(idx + 1) + k] = nv[k];
                }
            } else {
                dp.t = Datos::ARREGLO;
                std::vector<double> base = detalle_vent::a_f64(datos.a);
                const size_t n_fut = base.size() - static_cast<size_t>(idx + 1);
                std::vector<double> ruido(n_fut, 0.0);
                if (perturbacion_aleatoria) for (auto& r : ruido) r = -1e6 + (1e6 - -1e6) * rng.next_double();
                for (size_t k = 0; k < n_fut; ++k) base[static_cast<size_t>(idx + 1) + k] += magnitud_shock + ruido[k];
                dp.a = Arreglo::de_f64(std::move(base));
                dp.a.forma = datos.a.forma;
            }
            const auto pert = lista_ventanas(dp, timestamps, mascara);
            const VentanaCausalCertificada& vp = pert.at(ko);
            const double dv = maxabs(vo.valores, vp.valores);
            const auto pa = vo.valor_presente(), pb = vp.valor_presente();
            const double dpr = maxabs(pa, pb);
            const double drel = std::fabs(vo.reloj_acumulado - vp.reloj_acumulado);
            const double ddr = std::fabs(vo.delta_reloj_ventana - vp.delta_reloj_ventana);
            std::vector<double> ma(vo.mascara_huecos.begin(), vo.mascara_huecos.end()), mb(vp.mascara_huecos.begin(), vp.mascara_huecos.end());
            const double dm = maxabs(ma, mb);
            const double dt = std::fabs(static_cast<double>(vo.t_actual - vp.t_actual));
            const double dva = vo.es_valida == vp.es_valida ? 0.0 : 1.0;
            const double dz = maxabs(vo.normalizar_zscore_causal(), vp.normalizar_zscore_causal());
            const double dpc = pymax({dv, dpr, drel, ddr, dm, dt, dva, dz});
            mx_p = pymax({mx_p, dpc});
            double dtc = 0.0;
            if (verificar_truncado) {
                Datos dtr = datos;
                const size_t hasta = static_cast<size_t>(idx + 1);
                if (datos.t == Datos::MARCO) {
                    std::vector<size_t> filas(hasta);
                    for (size_t k = 0; k < hasta; ++k) filas[k] = k;
                    dtr.m.tabla = ingesta::tomar_filas(datos.m.tabla, filas);
                    if (dtr.m.indice.tipo != TipoIndice::RANGO && dtr.m.indice.valores.size() > hasta) dtr.m.indice.valores.resize(hasta);
                } else {
                    cortar(dtr.a, hasta);
                }
                std::optional<Arreglo> tst, mkt;
                if (timestamps) { tst = *timestamps; cortar(*tst, hasta); }
                if (mascara) { mkt = *mascara; cortar(*mkt, hasta); }
                std::vector<VentanaCausalCertificada> tr;
                generar_ventanas(dtr, tst, mkt, std::nullopt, std::nullopt, std::nullopt, std::nullopt,
                                 [&](const VentanaCausalCertificada& v) { tr.push_back(v); });
                if (!tr.empty()) {
                    const auto& vt = tr.back();
                    dtc = pymax({maxabs(vo.valores, vt.valores), std::fabs(vo.reloj_acumulado - vt.reloj_acumulado),
                                 std::fabs(vo.delta_reloj_ventana - vt.delta_reloj_ventana), std::fabs(static_cast<double>(vo.t_actual - vt.t_actual))});
                    mx_t = pymax({mx_t, dtc});
                }
            }
            const double dtot = pymax({dpc, dtc});
            if (dtot > mx_g) { mx_g = dtot; peor = idx; }
        }
        ultimo_lookahead_score = mx_g;
        ReporteTestigo r;
        r.max_diff = mx_g;
        r.cortes_evaluados = static_cast<int64_t>(cortes.size());
        r.lista_cortes = cortes;
        r.corte_peor_diff = peor;
        r.max_diff_perturbacion = mx_p;
        r.max_diff_truncado = mx_t;
        r.es_causal_estricto = mx_g == 0.0;
        return r;
    }
    double medir_lookahead_score(const Datos& d, const std::optional<Arreglo>& ts = std::nullopt, const std::optional<Arreglo>& mk = std::nullopt,
                                 int64_t max_cortes = 10) {
        return auditar_testigo_anti_lookahead(d, std::nullopt, false, ts, mk, 1e12, false, max_cortes).max_diff;
    }

    // ----------------------------------------------------------------------------------------- generar_ventanas_chunked
    void generar_ventanas_chunked(const std::vector<Trozo>& trozos, const Columnas& col, const std::optional<std::string>& venue,
                                  const std::optional<std::string>& nivel, const std::optional<std::string>& simbolo, const Emisor& emitir) {
        validar_metodo_imputacion_causal();
        stream_ventanas_core(trozos, col, venue, nivel, simbolo, emitir);
    }

    // ------------------------------------------------------------------------------------- materializar_matriz_ventanas
    std::vector<double> materializar_matriz_ventanas(const Datos& d, const std::optional<Arreglo>& ts = std::nullopt,
                                                     const std::optional<Arreglo>& mk = std::nullopt, const std::optional<std::string>& venue = std::nullopt,
                                                     const std::optional<std::string>& nivel = std::nullopt,
                                                     const std::optional<std::string>& simbolo = std::nullopt, int64_t* n_ventanas_out = nullptr) {
        const int64_t N = static_cast<int64_t>(d.len());
        const int64_t W = config.lookback_W, S = config.stride_S;
        if (N < W)
            throw ErrorDimensionInsuficiente("Longitud de datos insuficiente: N=" + std::to_string(N) + " < lookback_W=" + std::to_string(W) + ".");
        const int64_t nv = (N - W) / S + 1;
        const int64_t tot = nv * W;
        if (tot > config.max_elementos_materializables)
            throw ErrorAislamientoMemoria("VIOLACIÓN AXIOMA 5: Intento de materialización densa de " + std::to_string(tot) + " elementos (" +
                                          std::to_string(nv) + " ventanas x " + std::to_string(W) + " muestras) excede el umbral de seguridad de " +
                                          std::to_string(config.max_elementos_materializables) +
                                          " elementos. El aislamiento de memoria O(W) exige consumir las ventanas en streaming mediante generar_ventanas() o "
                                          "generar_ventanas_chunked() para evitar desbordamiento de RAM.");
        std::vector<double> mat(static_cast<size_t>(nv * W), 0.0);
        int64_t k = 0;
        generar_ventanas(d, ts, mk, venue, nivel, simbolo, std::nullopt, [&](const VentanaCausalCertificada& v) {
            if (k >= nv) throw IndexError("index " + std::to_string(k) + " is out of bounds for axis 0 with size " + std::to_string(nv));
            std::copy(v.valores.begin(), v.valores.end(), mat.begin() + k * W);
            ++k;
        });
        if (n_ventanas_out) *n_ventanas_out = nv;
        return mat;
    }

    // --------------------------------------------------------------------------- auditar_huella_memoria_streaming (sin tracemalloc)
    int64_t contar_ventanas_streaming(const std::function<void(const Emisor&)>& generador, int64_t n_max = 100000) {
        int64_t n = 0;
        struct Alto {};
        try {
            generador([&](const VentanaCausalCertificada&) { if (++n >= n_max) throw Alto{}; });
        } catch (const Alto&) {
        }
        return n;
    }

    // ----------------------------------------------------------------------------------------------- envoltorios del Axioma 6
    int64_t validar_escala_multiescala(int64_t s) const { return ark::dataset::ventanas::validar_escala_multiescala(config.lookback_W, s, config.forzar_divisores_exactos); }
    std::vector<double> construir_proyector_aniquilador(std::optional<int64_t> W = std::nullopt, std::optional<int64_t> orden = std::nullopt) const {
        return ark::dataset::ventanas::construir_proyector_aniquilador((W && *W) ? *W : config.lookback_W, orden ? *orden : config.orden_polinomio_detrending);
    }
    Segmentacion segmentar_kantelhardt_bidireccional(const std::vector<double>& v, int64_t s) const {
        return ark::dataset::ventanas::segmentar_kantelhardt_bidireccional(v, s, config.bidireccional_kantelhardt);
    }
    std::vector<double> desacoplar_rampa_convolucional_tv(const std::vector<double>& v, std::optional<double> lambda_tv = std::nullopt,
                                                          std::optional<int64_t> ancho = std::nullopt) const {
        return ark::dataset::ventanas::desacoplar_rampa_convolucional_tv(v, lambda_tv ? *lambda_tv : config.lambda_tv, ancho);
    }

private:
    void calc_residuo(const std::vector<int64_t>& divisores, int64_t& residuo_max, double& pct_descarte) const {
        const int64_t W = config.lookback_W;
        if (config.forzar_divisores_exactos) {
            residuo_max = 0;
            for (int64_t s : divisores) residuo_max = std::max(residuo_max, W % s);
            pct_descarte = 0.0;
        } else {
            const int64_t s_max = std::max(config.s_min, static_cast<int64_t>(detalle_vent::div_piso_f(static_cast<double>(W), config.max_ratio_escalas)));
            residuo_max = 0;
            for (int64_t s = config.s_min; s <= s_max; ++s) residuo_max = std::max(residuo_max, W % s);
            pct_descarte = (config.bidireccional_kantelhardt || residuo_max == 0) ? 0.0 : 100.0 * static_cast<double>(residuo_max) / static_cast<double>(W);
        }
    }
    static void cortar(Arreglo& a, size_t hasta) {
        if (a.dt == Dt::F64) { if (a.f.size() > hasta) a.f.resize(hasta); }
        else if (a.dt == Dt::OBJ || a.dt == Dt::USTR) { if (a.o.size() > hasta) a.o.resize(hasta); }
        else if (a.i.size() > hasta) a.i.resize(hasta);
        if (!a.forma.empty()) a.forma[0] = std::min<int64_t>(a.forma[0], static_cast<int64_t>(hasta));
    }
};

// Nombre de la clase Python de una excepcion del porte.
inline std::string clase_error(const std::exception& e) {
    if (const auto* p = dynamic_cast<const ErrorPy*>(&e)) return p->clase();
    return ingesta::clase_excepcion(e);
}

}  // namespace ark::dataset::ventanas
