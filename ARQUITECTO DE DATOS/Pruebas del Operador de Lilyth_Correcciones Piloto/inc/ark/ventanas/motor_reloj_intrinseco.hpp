// motor_reloj_intrinseco.hpp -- ventanas rodantes en RELOJ INTRINSECO con W de ejecucion (ARK-WIN-01-CXX), C++20, zero-heap.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04, por orden del Comandante: "actualiza el motor a reloj intrinseco y haz que su
// W no sea fija, sino que sea un input del operador basado en el pipeline ya hecho de las ventanas rodantes".
//
// Lo que reutiliza (no duplica) del pipeline ya hecho: la ADUANA fila a fila de MotorVentanasCausalesCxx (monotonia, colision de
// sub-barra, disponibilidad R1, imputacion causal sin bfill, sesiones, huecos por calendario o cadencia) y su definicion del reloj
// d_nu por fila (PARED 1, NU_VOLUMEN, NU_TRADES, NU_NOCIONAL, NU_INFO_SHANNON con su lazo; d_nu = 0 en filas imputadas si
// congelar_en_huecos, y en filas sin sesion). El motor de pared se usa como COMPUERTA: cada fila pasa por su paso() transaccional y
// se lee su FilaAdmitidaCxx. Si la aduana rechaza la fila, este motor tampoco cambia.
//
// Lo nuevo: barras INTRINSECAS. Las fronteras de barra son los multiplos k delta_nu del nu acumulado por las filas que entran en
// barras (CON ARRASTRE del exceso, orden del Comandante 2026-11-04): la barra se cierra en la primera fila en que el reloj cruza la
// frontera siguiente, y el exceso sobre esa frontera queda como credito de la barra siguiente. Una fila no se reparte entre barras
// (declarado: repartirla inventaria datos): si una sola fila cruza k > 1 fronteras, cierra UNA barra y las k - 1 fronteras restantes
// quedan sin barra; se cuentan (`fronteras_barra`, `fronteras_saltadas`). Asi las fronteras de un ancho r delta_nu (r entero) son
// siempre fronteras del ancho base. Las filas sin sesion no entran en ninguna barra (reloj congelado). Con PARED y delta_nu = 1,
// cada fila con sesion y sin imputar es una barra: el control de identidad exige reproducir el motor de pared bit a bit en ese caso.
// Valor de la barra: ULTIMO (niveles), SUMA (flujos: una fila imputada aporta 0) o MEDIA_NU (intensivos, ponderada por d_nu).
// Ventana: las ultimas W barras (W y el paso S en barras son ENTRADAS de ejecucion, W <= W_MAX). Presupuesto de imputacion contado
// en FILAS con la formula del motor de pared: (imputadas con sesion + omitidas) / (filas con sesion + omitidas), sin el hueco
// anterior a la primera fila de la ventana. Sello de la ventana = maxima disponibilidad de sus filas.
// DECLARADO: el lazo de Shannon de la compuerta va acotado por W_MAX (la capacidad), no por la W de ejecucion; abortar por
// presupuesto no se admite aqui (la compuerta ya confirmo la fila cuando se conoce el presupuesto de la barra).
#pragma once

#include "ark/ventanas/contratos_ventanas.hpp"
#include "ark/ventanas/motor_ventanas.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

namespace ark::ventanas {

enum class AgregacionBarra : uint8_t { ULTIMO = 0, SUMA = 1, MEDIA_NU = 2 };

enum class ErrorMotorIntrinseco : uint8_t { NINGUNO = 0, BASE_INVALIDA, W_FUERA_DE_RANGO, S_NULO, DELTA_NU_INVALIDO, ABORTAR_NO_SOPORTADO };

template <size_t W_MAX>
struct ConfiguracionMotorIntrinsecoCxx {
    ConfiguracionVentanaCausalCxx<W_MAX> base{};   // aduana, reloj, imputacion, umbral; su stride_S se ignora
    size_t W = 0;                                  // OBLIGATORIO: lo entrega el operador (barras), min_muestras_requeridas <= W <= W_MAX
    size_t S = 1;                                  // paso de emision, en barras intrinsecas
    double delta_nu = std::numeric_limits<double>::quiet_NaN();   // OBLIGATORIO: ancho de barra en unidades de nu
    AgregacionBarra agregacion = AgregacionBarra::ULTIMO;

    [[nodiscard]] ErrorMotorIntrinseco validar() const noexcept {
        if (base.validar() != ErrorConfiguracion::NINGUNO) return ErrorMotorIntrinseco::BASE_INVALIDA;
        if (W < base.min_muestras_requeridas || W > W_MAX) return ErrorMotorIntrinseco::W_FUERA_DE_RANGO;
        if (S < 1) return ErrorMotorIntrinseco::S_NULO;
        if (!(delta_nu > 0.0 && std::isfinite(delta_nu))) return ErrorMotorIntrinseco::DELTA_NU_INVALIDO;
        if (base.abortar_si_excede_imputacion) return ErrorMotorIntrinseco::ABORTAR_NO_SOPORTADO;
        return ErrorMotorIntrinseco::NINGUNO;
    }
};

// Ventana certificada en reloj intrinseco: las n = W barras mas recientes (la mas vieja en [0]).
template <size_t W_MAX>
struct VentanaIntrinsecaCxx {
    size_t n{0};
    uint64_t t_actual{0};                       // maxima disponibilidad de las filas de la ventana (B5)
    uint64_t t_cierre_ms{0};                    // cierre de la ultima fila de la ultima barra
    bool disponibilidad_declarada{false};
    std::array<double, W_MAX> valores{};
    std::array<bool, W_MAX> mascara{};          // la barra contiene alguna fila con sesion imputada
    std::array<double, W_MAX> nu_barra{};       // nu de las filas de cada barra (con arrastre puede ser < delta_nu o > delta_nu)
    std::array<uint32_t, W_MAX> filas_barra{};  // filas con sesion que forman cada barra
    std::array<uint32_t, W_MAX> fronteras_barra{};   // fronteras k delta_nu que cruza cada barra (>= 1)
    uint32_t fronteras_saltadas{0};             // en la ventana: suma de (fronteras_barra - 1)
    double reloj_acumulado{0.0};                // nu acumulado de la compuerta
    double delta_reloj_ventana{0.0};            // suma de nu_barra
    bool hay_imputadas{false};
    bool es_valida{true};
    double ratio_imputacion{0.0};
    uint32_t barras_omitidas{0};                // en FILAS (calendario o cadencia), dentro de la ventana
    uint32_t filas{0};                          // filas con sesion en la ventana
    TipoRelojVentana tipo_reloj{TipoRelojVentana::PARED};
    double delta_nu{0.0};
    std::array<uint8_t, W_MAX> flags_hueco{};       // OR de flags_hueco de las filas con sesion que forman cada barra
    std::array<bool, W_MAX> solo_mercado_quieto{};  // true si la barra es imputada pero todas sus filas imputadas son MERCADO_QUIETO sin corrupcion ni omitidas
};

struct PasoIntrinsecoCxx {
    ResultadoPaso resultado{ResultadoPaso::SIN_EMISION};   // error de la aduana, SIN_EMISION o VENTANA_EMITIDA
    bool barra_cerrada{false};
};

template <size_t W_MAX>
class MotorVentanasIntrinsecoCxx {
    static_assert(W_MAX >= 4, "W_MAX >= 4");

    ConfiguracionMotorIntrinsecoCxx<W_MAX> cfg_{};
    ErrorMotorIntrinseco error_{ErrorMotorIntrinseco::BASE_INVALIDA};
    MotorVentanasCausalesCxx<W_MAX> compuerta_{};
    VentanaCausalCertificadaCxx<W_MAX> descarte_{};   // ventana de la compuerta (no se usa: su presupuesto es en filas de pared)

    // barras cerradas (anillo espejo de capacidad W_MAX; la ventana son sus ultimas W)
    RingBufferEspejoCxx<double, W_MAX> r_val_{};
    RingBufferEspejoCxx<bool, W_MAX> r_mask_{};
    RingBufferEspejoCxx<double, W_MAX> r_nu_{};
    RingBufferEspejoCxx<uint32_t, W_MAX> r_filas_{}, r_imp_{}, r_omit_{}, r_gap0_{}, r_front_{};
    RingBufferEspejoCxx<uint64_t, W_MAX> r_disp_{}, r_cierre_{};
    RingBufferEspejoCxx<bool, W_MAX> r_decl_{}, r_quieto_{};
    RingBufferEspejoCxx<uint8_t, W_MAX> r_flags_{};

    // barra en curso
    double a_nu_ = 0.0, a_val_ = 0.0, a_pond_ = 0.0;
    bool a_mask_ = false, a_decl_ = true, a_vacia_ = true;
    uint32_t a_filas_ = 0, a_imp_ = 0, a_imp_quieto_ = 0, a_omit_ = 0, a_gap0_ = 0;
    uint8_t a_flags_ = 0;
    uint64_t a_disp_ = 0, a_cierre_ = 0;
    double credito_ = 0.0;   // exceso sobre la ultima frontera cruzada, en [0, delta_nu): nu que ya cuenta para la barra siguiente

    size_t barras_ = 0;
    size_t siguiente_ = 0;
    uint32_t emitidas_ = 0, invalidadas_ = 0;

    void reiniciar_barra() noexcept {
        a_nu_ = a_val_ = a_pond_ = 0.0;
        a_mask_ = false; a_decl_ = true; a_vacia_ = true;
        a_filas_ = a_imp_ = a_imp_quieto_ = a_omit_ = a_gap0_ = 0;
        a_flags_ = 0;
        a_disp_ = a_cierre_ = 0;
    }

public:
    MotorVentanasIntrinsecoCxx() noexcept = default;
    explicit MotorVentanasIntrinsecoCxx(const ConfiguracionMotorIntrinsecoCxx<W_MAX>& c) noexcept { configurar(c); }

    // Configura y reinicia. La W y el paso los entrega el operador (por ejemplo, la W derivada de su banda).
    ErrorMotorIntrinseco configurar(const ConfiguracionMotorIntrinsecoCxx<W_MAX>& c) noexcept {
        cfg_ = c;
        error_ = c.validar();
        ConfiguracionVentanaCausalCxx<W_MAX> b = c.base;
        b.stride_S = std::numeric_limits<size_t>::max() / 4;   // la compuerta casi nunca emite: solo se lee su fila admitida
        compuerta_ = MotorVentanasCausalesCxx<W_MAX>(b);
        reset();
        return error_;
    }
    [[nodiscard]] ErrorMotorIntrinseco error_configuracion() const noexcept { return error_; }
    [[nodiscard]] const ConfiguracionMotorIntrinsecoCxx<W_MAX>& configuracion() const noexcept { return cfg_; }
    [[nodiscard]] size_t barras() const noexcept { return barras_; }
    [[nodiscard]] uint32_t ventanas_emitidas() const noexcept { return emitidas_; }
    [[nodiscard]] uint32_t ventanas_invalidadas() const noexcept { return invalidadas_; }
    [[nodiscard]] double nu_acumulado() const noexcept { return compuerta_.nu_acumulado(); }
    [[nodiscard]] const MotorVentanasCausalesCxx<W_MAX>& compuerta() const noexcept { return compuerta_; }

    void reset() noexcept {
        compuerta_.reset();
        r_val_.clear(); r_mask_.clear(); r_nu_.clear(); r_filas_.clear(); r_imp_.clear(); r_omit_.clear(); r_gap0_.clear();
        r_front_.clear(); r_disp_.clear(); r_cierre_.clear(); r_decl_.clear(); r_quieto_.clear(); r_flags_.clear();
        reiniciar_barra();
        credito_ = 0.0;
        barras_ = 0;
        siguiente_ = cfg_.W;
        emitidas_ = invalidadas_ = 0;
    }

    // Una fila nativa. Transaccional: si la aduana la rechaza, ni la compuerta ni este motor cambian.
    [[nodiscard]] PasoIntrinsecoCxx paso(const EntradaPasoCxx& e, VentanaIntrinsecaCxx<W_MAX>& out) noexcept {
        PasoIntrinsecoCxx p{};
        if (error_ != ErrorMotorIntrinseco::NINGUNO) { p.resultado = ResultadoPaso::ERROR_CONFIGURACION_INVALIDA; return p; }
        const ResultadoPaso rc = compuerta_.paso(e, descarte_);
        if (es_error(rc) || rc == ResultadoPaso::FILA_DESCARTADA_INTERSECCION) { p.resultado = rc; return p; }
        const FilaAdmitidaCxx& f = compuerta_.ultima_fila();
        // huecos: cuentan aunque la fila no tenga sesion (es informacion que falta)
        if (a_vacia_) a_gap0_ += f.omitidas; else a_omit_ += f.omitidas;
        if (!f.con_sesion) return p;   // reloj congelado: la fila no entra en ninguna barra
        if (a_vacia_) { a_vacia_ = false; }
        ++a_filas_;
        a_flags_ |= f.flags_hueco;
        if (f.mascara) {
            ++a_imp_;
            a_mask_ = true;
            // Es mercado quieto puro si trae FLAG_MERCADO_QUIETO (8) sin FLAG_FILA_VACIA (2), FLAG_PAUSA_UMBRAL (64) ni FLAG_CORRUPCION_CAIDA (128)
            const bool fila_quieta = ((f.flags_hueco & 8u) != 0u) && ((f.flags_hueco & (2u | 64u | 128u)) == 0u) && std::isfinite(e.valor);
            if (fila_quieta) ++a_imp_quieto_;
        }
        a_nu_ += f.d_nu;
        switch (cfg_.agregacion) {
            case AgregacionBarra::ULTIMO: a_val_ = f.valor; break;
            case AgregacionBarra::SUMA: if (!f.mascara) a_val_ += f.valor; break;
            case AgregacionBarra::MEDIA_NU: a_val_ += f.valor * f.d_nu; a_pond_ += f.d_nu; break;
        }
        a_disp_ = std::max(a_disp_, f.disponible_ms);
        a_decl_ = a_decl_ && f.disponibilidad_declarada;
        a_cierre_ = f.ts_ms;
        const double reloj = credito_ + a_nu_;   // nu desde la ultima frontera cruzada
        if (!(reloj >= cfg_.delta_nu)) return p;
        // cierre de la barra: k fronteras cruzadas; el exceso sobre la ultima pasa a la barra siguiente
        double k = std::floor(reloj / cfg_.delta_nu);
        if (k < 1.0) k = 1.0;   // redondeo: reloj >= delta_nu garantiza al menos una
        credito_ = reloj - k * cfg_.delta_nu;
        if (credito_ < 0.0) credito_ = 0.0;
        if (credito_ >= cfg_.delta_nu) credito_ = 0.0;   // no ocurre salvo redondeo en el ultimo bit
        const uint32_t fronteras = k >= 4294967295.0 ? 4294967295u : static_cast<uint32_t>(k);
        const double valor = cfg_.agregacion == AgregacionBarra::MEDIA_NU ? (a_pond_ > 0.0 ? a_val_ / a_pond_ : std::numeric_limits<double>::quiet_NaN()) : a_val_;
        const bool bar_solo_quieto = a_mask_ && (a_imp_ > 0) && (a_imp_quieto_ == a_imp_) && (a_omit_ == 0) &&
                                     ((a_flags_ & (2u | 64u | 128u)) == 0u);
        r_val_.push(valor); r_mask_.push(a_mask_); r_nu_.push(a_nu_); r_filas_.push(a_filas_); r_imp_.push(a_imp_);
        r_omit_.push(a_omit_); r_gap0_.push(a_gap0_); r_front_.push(fronteras); r_disp_.push(a_disp_); r_cierre_.push(a_cierre_); r_decl_.push(a_decl_);
        r_flags_.push(a_flags_); r_quieto_.push(bar_solo_quieto);
        ++barras_;
        reiniciar_barra();
        p.barra_cerrada = true;
        if (barras_ < siguiente_) return p;
        siguiente_ += cfg_.S;
        // emision de las ultimas W barras
        const size_t W = cfg_.W, o = W_MAX - W;
        const double* v = r_val_.ventana_actual() + o;
        const bool* m = r_mask_.ventana_actual() + o;
        const double* nu = r_nu_.ventana_actual() + o;
        const uint32_t* fl = r_filas_.ventana_actual() + o;
        const uint32_t* im = r_imp_.ventana_actual() + o;
        const uint32_t* om = r_omit_.ventana_actual() + o;
        const uint32_t* g0 = r_gap0_.ventana_actual() + o;
        const uint32_t* fr = r_front_.ventana_actual() + o;
        const uint64_t* dp = r_disp_.ventana_actual() + o;
        const bool* dc = r_decl_.ventana_actual() + o;
        const uint8_t* fg = r_flags_.ventana_actual() + o;
        const bool* sq = r_quieto_.ventana_actual() + o;
        out.n = W;
        std::memcpy(out.valores.data(), v, W * sizeof(double));
        std::memcpy(out.mascara.data(), m, W * sizeof(bool));
        std::memcpy(out.nu_barra.data(), nu, W * sizeof(double));
        std::memcpy(out.filas_barra.data(), fl, W * sizeof(uint32_t));
        std::memcpy(out.fronteras_barra.data(), fr, W * sizeof(uint32_t));
        std::memcpy(out.flags_hueco.data(), fg, W * sizeof(uint8_t));
        std::memcpy(out.solo_mercado_quieto.data(), sq, W * sizeof(bool));
        uint32_t filas = 0, imp = 0, omit = 0, salt = 0;
        double dnu = 0.0;
        uint64_t td = 0;
        bool decl = true, hay = false;
        for (size_t i = 0; i < W; ++i) {
            filas += fl[i]; imp += im[i]; omit += om[i] + (i ? g0[i] : 0u);   // el hueco anterior a la primera fila queda fuera
            dnu += nu[i]; td = std::max(td, dp[i]); decl = decl && dc[i]; hay = hay || m[i]; salt += fr[i] - 1u;
        }
        out.fronteras_saltadas = salt;
        const double den = static_cast<double>(filas + omit);
        const double ratio = den > 0.0 ? static_cast<double>(imp + omit) / den : std::numeric_limits<double>::quiet_NaN();
        out.t_actual = td;
        out.t_cierre_ms = r_cierre_.ventana_actual()[W_MAX - 1];
        out.disponibilidad_declarada = decl;
        out.reloj_acumulado = compuerta_.nu_acumulado();
        out.delta_reloj_ventana = dnu;
        out.hay_imputadas = hay;
        out.ratio_imputacion = ratio;
        out.es_valida = ratio <= cfg_.base.umbral_imputacion_maximo;
        out.barras_omitidas = omit;
        out.filas = filas;
        out.tipo_reloj = cfg_.base.reloj;
        out.delta_nu = cfg_.delta_nu;
        if (!out.es_valida) ++invalidadas_;
        ++emitidas_;
        p.resultado = ResultadoPaso::VENTANA_EMITIDA;
        return p;
    }
};

// Ancho de barra derivado de la historia (calibracion; reserva nada): nu medio por fila con sesion y sin imputar, para que el
// motor intrinseco emita en promedio una barra por fila de pared. Es una entrada del operador; esta es la opcion por omision.
template <size_t W_MAX>
[[nodiscard]] inline double delta_nu_por_fila(const ConfiguracionVentanaCausalCxx<W_MAX>& base, const EntradaPasoCxx* filas, size_t n) noexcept {
    ConfiguracionVentanaCausalCxx<W_MAX> b = base;
    b.stride_S = std::numeric_limits<size_t>::max() / 4;
    b.abortar_si_excede_imputacion = false;
    MotorVentanasCausalesCxx<W_MAX> g(b);
    VentanaCausalCertificadaCxx<W_MAX> d{};
    double nu = 0.0;
    size_t k = 0;
    for (size_t i = 0; i < n; ++i) {
        const ResultadoPaso r = g.paso(filas[i], d);
        if (es_error(r) || r == ResultadoPaso::FILA_DESCARTADA_INTERSECCION) continue;
        const FilaAdmitidaCxx& f = g.ultima_fila();
        if (f.con_sesion && !f.mascara) { nu += f.d_nu; ++k; }
    }
    return k ? nu / static_cast<double>(k) : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace ark::ventanas
