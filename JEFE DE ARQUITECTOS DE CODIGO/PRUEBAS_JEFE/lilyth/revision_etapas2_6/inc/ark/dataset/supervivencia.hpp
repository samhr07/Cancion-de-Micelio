// supervivencia.hpp -- CICLO DE VIDA PUNTO EN EL TIEMPO DE CADA ACTIVO (ARK-DATA-01-CXX)
// Reactor ARK / Cancion del Micelio -- ARQUITECTO DE DATOS (Flash) y correccion del Jefe 2026-11-01
//
// Porte EXACTO de catalogo_universo.CATALOGO_HISTORICO_BASE y de auditar_dataframe_supervivencia (la referencia canonica):
// antes del listado la fila es PRE_LISTADO y despues del desliste DESLISTADO; en las dos, sin sesion, sin dato, imputada y
// con los precios en NaN (no 0.0: un 0.0 se lee como precio y corta el arrastre de forma silenciosa).
//
// La version anterior tenia OTRAS fechas que Python (SPX 1927 contra 1957, DXY 1971 contra 1973, VIX 1990 contra 1993, LUNA
// 2019-08-26/2022-05-12 contra 2020-08-27/2022-05-13, FTT 2019 contra 2021-01-15...) y BTC en vez de EEM. Ahora es la misma
// tabla. Las fechas siguen sin fuente primaria citada en el catalogo de Python: son las del Arquitecto de Datos.
// Un activo fuera del catalogo queda COTIZANDO (como en Python): marcar supervivencia exige catalogarlo.
#pragma once

#include "tipos_dataset.hpp"
#include <cstdint>
#include <cstring>
#include <cctype>

namespace ark::dataset {

enum class MotivoDesliste : uint8_t {
    ACTIVO_VIVO = 0,
    DESLISTADO_REGULATORIO = 1,
    QUIEBRA_COLAPSO = 2,
    FUSION_ADQUISICION = 3,
    FIN_CONTRATO_FUTURO = 4,
    DESLISTADO_VOLUMEN = 5
};

struct RegistroCicloVidaActivo {
    const char* simbolo;
    uint32_t fecha_listado;   // YYYYMMDD
    uint32_t fecha_desliste;  // YYYYMMDD; 0 si sigue vivo
    MotivoDesliste motivo;
    const char* venue;
    const char* descripcion;
};

inline constexpr RegistroCicloVidaActivo CATALOGO_HISTORICO_CXX[] = {
    {"SPX", 19570304, 0, MotivoDesliste::ACTIVO_VIVO, "YAHOO", "S&P 500 Index"},
    {"NDX", 19850131, 0, MotivoDesliste::ACTIVO_VIVO, "YAHOO", "Nasdaq 100 Index"},
    {"VIX", 19930102, 0, MotivoDesliste::ACTIVO_VIVO, "YAHOO", "Cboe Volatility Index"},
    {"DXY", 19730301, 0, MotivoDesliste::ACTIVO_VIVO, "BLOOMBERG", "US Dollar Index"},
    {"ES1", 19970909, 0, MotivoDesliste::ACTIVO_VIVO, "BLOOMBERG", "E-mini S&P 500 Future"},
    {"HYG", 20070411, 0, MotivoDesliste::ACTIVO_VIVO, "YAHOO", "iShares iBoxx High Yield"},
    {"IEF", 20020730, 0, MotivoDesliste::ACTIVO_VIVO, "YAHOO", "iShares 7-10 Year Treasury"},
    {"EEM", 20030411, 0, MotivoDesliste::ACTIVO_VIVO, "YAHOO", "iShares MSCI Emerging Markets"},
    {"VNQ", 20040929, 0, MotivoDesliste::ACTIVO_VIVO, "YAHOO", "Vanguard Real Estate ETF"},
    {"BTCUSDT", 20190908, 0, MotivoDesliste::ACTIVO_VIVO, "BINANCE_UM", "Bitcoin USD-M Perp"},
    {"ETHUSDT", 20191128, 0, MotivoDesliste::ACTIVO_VIVO, "BINANCE_UM", "Ethereum USD-M Perp"},
    {"SOLUSDT", 20200924, 0, MotivoDesliste::ACTIVO_VIVO, "BINANCE_UM", "Solana USD-M Perp"},
    {"LUNAUSDT_OLD", 20200827, 20220513, MotivoDesliste::QUIEBRA_COLAPSO, "BINANCE_UM", "Terra Classic Colapso 2022"},
    {"FTTUSDT", 20210115, 20221114, MotivoDesliste::QUIEBRA_COLAPSO, "BINANCE_UM", "FTX Token Quiebra 2022"},
};

class CatalogoSupervivenciaCxx {
public:
    [[nodiscard]] static const RegistroCicloVidaActivo* ficha(const char* simbolo) noexcept {
        if (!simbolo) return nullptr;
        char u[32];
        size_t k = 0;
        for (; simbolo[k] && k < sizeof(u) - 1; ++k) u[k] = static_cast<char>(std::toupper(static_cast<unsigned char>(simbolo[k])));
        u[k] = '\0';
        for (const auto& r : CATALOGO_HISTORICO_CXX)
            if (std::strcmp(r.simbolo, u) == 0) return &r;
        return nullptr;
    }

    [[nodiscard]] static bool en_catalogo(const char* simbolo) noexcept { return ficha(simbolo) != nullptr; }

    [[nodiscard]] static EstadoSupervivencia evaluar(const char* simbolo, uint32_t fecha_sesion) noexcept {
        const RegistroCicloVidaActivo* r = ficha(simbolo);
        if (!r) return EstadoSupervivencia::COTIZANDO;
        if (fecha_sesion < r->fecha_listado) return EstadoSupervivencia::PRE_LISTADO;
        if (r->fecha_desliste != 0 && fecha_sesion > r->fecha_desliste) return EstadoSupervivencia::DESLISTADO;
        return EstadoSupervivencia::COTIZANDO;
    }

    // auditar_dataframe_supervivencia: fuera del ciclo de vida, sin sesion, sin dato, imputada y precios NaN.
    static void aplicar_mascara(const char* simbolo, BarraCausal* barras, size_t n) noexcept {
        if (!barras || n == 0 || !simbolo) return;
        const bool cat = en_catalogo(simbolo);
        for (size_t i = 0; i < n; ++i) {
            BarraCausal& b = barras[i];
            const EstadoSupervivencia st = evaluar(simbolo, b.fecha_sesion);
            b.estado_supervivencia = st;
            if (!cat || st == EstadoSupervivencia::COTIZANDO) continue;
            b.es_sesion_con_dato = false;
            b.hay_dato = false;
            b.es_imputado = true;
            b.open = b.high = b.low = b.close = NAN_DATASET;
            b.tipo_hueco = TipoHueco::FUERA_DE_CICLO;
        }
    }
};

} // namespace ark::dataset
