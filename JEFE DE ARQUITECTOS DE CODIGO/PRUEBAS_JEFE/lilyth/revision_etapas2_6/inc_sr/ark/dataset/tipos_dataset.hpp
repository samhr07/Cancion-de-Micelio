// tipos_dataset.hpp -- TIPOS DEL DATASET C++ (ARK-DATA-01-CXX)
// Reactor ARK / Cancion del Micelio -- ARQUITECTO DE DATOS (Flash) y correccion del Jefe 2026-11-01
//
// BarraCausal es una fila de la REJILLA DIARIA del dataset con la misma semantica que la tabla de motor_ingesta.py:
//   - estampa_disponible_ms : 00:00 UTC del dia siguiente a la sesion (regla R1). Antes no se rellenaba nunca (0 en el 100 %).
//   - es_imputado           : la fila no trae dato propio (fila vacia del proveedor, dia insertado por la rejilla, precio
//                             invalido o fuera del ciclo de vida); su precio va ARRASTRADO hacia delante (o NaN fuera del ciclo).
//   - hay_dato              : la fuente trajo la fila y su precio es valido.
//   - es_sesion_con_dato    : hay_dato y el calendario del activo dice que ese dia es sesion.
//   - dias_acumulados / sesiones_acumuladas : en una fila con sesion, dias naturales y sesiones del calendario desde la sesion
//                             con dato anterior (la convencion viernes -> lunes); 0 en las demas filas.
//   - cierre_lo/hi          : intervalo UTC del cierre de la sesion cuya informacion lleva la fila (ventanas_informacion).
// Precio NaN = no hay precio (antes 0.0, que se confundia con un precio y borraba el retorno siguiente).
#pragma once

#include <cstdint>
#include <cstddef>
#include <cmath>
#include <limits>
#include <span>
#include <array>
#include <algorithm>
#include <cstring>
#include <type_traits>

namespace ark::dataset {

enum class Venue : uint8_t {
    YAHOO = 0,
    BLOOMBERG = 1,
    BINANCE = 2,
    BITFINEX = 3
};

enum class EstadoSupervivencia : uint8_t {
    PRE_LISTADO = 0,
    COTIZANDO = 1,
    DESLISTADO = 2
};

[[nodiscard]] constexpr const char* nombre_estado(EstadoSupervivencia e) noexcept {
    return e == EstadoSupervivencia::PRE_LISTADO ? "PRE_LISTADO" : (e == EstadoSupervivencia::COTIZANDO ? "COTIZANDO" : "DESLISTADO");
}

enum class TipoHueco : uint8_t {
    NINGUNO = 0,
    CIERRE_PROGRAMADO = 1,   // feriado o cierre extraordinario del calendario
    FIN_DE_SEMANA = 2,       // dia no operable del calendario
    CORTE_FEED = 3,          // dia de sesion sin dato valido del proveedor
    DESCONEXION = 4,         // reservado para el intradia
    FUERA_DE_CICLO = 5       // antes del listado o despues del desliste
};

enum class ModoUniverso : uint8_t {
    ESTATICO = 0,
    PUNTO_EN_EL_TIEMPO = 1
};

enum class TipoObservableDataset : uint8_t {
    PRECIO_CIERRE = 0,
    RETORNO_LOG = 1,
    RETORNO_APERTURA_A_CIERRE = 2,
    VOLUMEN = 3,
    VOLUMEN_LOG = 4,
    FACTOR_NU = 5
};

constexpr double NAN_DATASET = std::numeric_limits<double>::quiet_NaN();

// Bits de flags_control
constexpr uint8_t FLAG_FILA_INSERTADA = 1;   // la fila la inserto la rejilla (la fuente no la trae)
constexpr uint8_t FLAG_FILA_VACIA = 2;       // la fuente trae la fila pero sin precio valido
constexpr uint8_t FLAG_FRACCION_NAN = 4;     // sesion reducida sin duracion modelada

// =======================================================================
// BarraCausal: estructura fija, trivialmente copiable, sin punteros ni memoria dinamica.
// =======================================================================
#pragma pack(push, 8)
struct BarraCausal {
    uint64_t estampa_disponible_ms;   // 00:00 UTC de D+1 (regla R1)
    uint32_t fecha_sesion;            // YYYYMMDD
    uint32_t reservado_padding;

    double open;
    double high;
    double low;
    double close;
    double volume;
    double valor_nominal;             // close * volume (nocional aproximado; ver extractor)
    double fraccion_sesion;           // 1.0 regular, 0.0 cerrado, p. ej. 3.5/6.5 en media jornada; NaN si no modelada

    uint64_t cierre_real_utc_ms;      // = cierre_hi_utc_ms (la hora de cierre que publica motor_ingesta: cierre_real_utc)
    uint64_t apertura_real_utc_ms;    // = inicio_lo de la barra; 0 si el instrumento no tiene apertura (un fixing)
    int64_t cierre_lo_utc_ms;         // intervalo del cierre (ventanas_informacion): lo
    int64_t cierre_hi_utc_ms;         //                                              hi

    uint16_t dias_acumulados;         // dias naturales desde la sesion con dato anterior (0 si la fila no es sesion con dato)
    uint16_t sesiones_acumuladas;     // sesiones del calendario en (anterior, actual]
    bool es_sesion_con_dato;
    bool es_imputado;
    bool hay_dato;
    EstadoSupervivencia estado_supervivencia;
    TipoHueco tipo_hueco;
    uint8_t tipo_dia;                 // TipoDia de calendario.hpp
    uint8_t flags_control;
    uint8_t reservado[5];

    uint8_t hash_sha256[32];          // SHA-256 binario de la barra fisica de la fuente (0 si la fila no trae dato)

    [[nodiscard]] double retorno_apertura_a_cierre() const noexcept {
        if (!es_sesion_con_dato || !(open > 0.0) || !(close > 0.0)) return 0.0;
        return std::log(close) - std::log(open);
    }

    [[nodiscard]] constexpr bool es_ejecutable() const noexcept {
        return es_sesion_con_dato && (estado_supervivencia == EstadoSupervivencia::COTIZANDO);
    }
};
#pragma pack(pop)

static_assert(sizeof(BarraCausal) == 152, "BarraCausal cambio de tamano: subir la version del serializador.");
static_assert(std::is_trivially_copyable_v<BarraCausal>, "BarraCausal tiene que ser trivialmente copiable.");

// =======================================================================
// ContratoExtraccionCxx: Parametros inmutables de consulta de datos
// =======================================================================
struct ContratoExtraccionCxx {
    Venue venue;
    char simbolo[16];
    uint32_t fecha_inicio; // YYYYMMDD
    uint32_t fecha_fin;    // YYYYMMDD
    ModoUniverso modo_universo;
    bool estricto_sin_huecos;

    constexpr ContratoExtraccionCxx()
        : venue(Venue::YAHOO), simbolo{0}, fecha_inicio(0), fecha_fin(0),
          modo_universo(ModoUniverso::ESTATICO), estricto_sin_huecos(false) {}

    constexpr ContratoExtraccionCxx(Venue v, const char* s, uint32_t fi, uint32_t ff, ModoUniverso mu = ModoUniverso::ESTATICO)
        : venue(v), simbolo{0}, fecha_inicio(fi), fecha_fin(ff),
          modo_universo(mu), estricto_sin_huecos(false) {
        if (s) {
            size_t i = 0;
            while (i < 15 && s[i] != '\0') {
                simbolo[i] = s[i];
                ++i;
            }
            simbolo[i] = '\0';
        }
    }
};

// =======================================================================
// VentanaCausalBuffer: Ring-buffer estatico en memoria continua (cero heap). Empujar es O(1).
// =======================================================================
template <size_t N>
class VentanaCausalBuffer {
    static_assert(N >= 1, "La ventana necesita al menos una barra.");
private:
    std::array<BarraCausal, N> buffer_{};
    size_t head_{0};
    size_t count_{0};

public:
    constexpr VentanaCausalBuffer() noexcept = default;

    void push(const BarraCausal& b) noexcept {
        buffer_[head_] = b;
        head_ = (head_ + 1) % N;
        if (count_ < N) ++count_;
    }

    [[nodiscard]] constexpr size_t size() const noexcept { return count_; }
    [[nodiscard]] constexpr size_t capacity() const noexcept { return N; }
    [[nodiscard]] constexpr bool full() const noexcept { return count_ == N; }
    [[nodiscard]] constexpr bool empty() const noexcept { return count_ == 0; }

    void clear() noexcept {
        head_ = 0;
        count_ = 0;
    }

    // Acceso causal: idx = 0 es la barra mas reciente (t), idx = 1 es (t-1)... Fuera de rango devuelve una barra vacia
    // (antes devolvia en silencio una barra vieja o sin escribir).
    [[nodiscard]] const BarraCausal& operator[](size_t idx) const noexcept {
        static const BarraCausal vacia{};
        if (idx >= count_) return vacia;
        size_t pos = (head_ + N - 1 - idx) % N;
        return buffer_[pos];
    }

    // Retornos log del cierre, mas reciente primero: out_r[i] = ln C(t-i) - ln C(t-i-1), con el cierre arrastrado de la
    // rejilla (una fila sin sesion repite el cierre y da 0; la sesion siguiente cobra el movimiento entero). NaN si falta un
    // precio (fuera del ciclo de vida) o si no es un precio: P <= 0 o infinito no tiene logaritmo (mandato P <= 0, 2026-11-02;
    // antes log(max(P, 1e-12))).
    size_t extraer_retornos_log(double* out_r, size_t max_out) const noexcept {
        if (!out_r || count_ < 2) return 0;
        size_t n = std::min(max_out, count_ - 1);
        auto valido = [](double p) { return std::isfinite(p) && p > 0.0; };
        for (size_t i = 0; i < n; ++i) {
            const double c1 = (*this)[i].close, c0 = (*this)[i + 1].close;
            out_r[i] = (valido(c1) && valido(c0)) ? std::log(c1) - std::log(c0) : NAN_DATASET;
        }
        return n;
    }
};

} // namespace ark::dataset
