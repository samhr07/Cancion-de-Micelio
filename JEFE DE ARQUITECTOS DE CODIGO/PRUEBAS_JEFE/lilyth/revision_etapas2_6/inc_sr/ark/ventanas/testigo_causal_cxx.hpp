// testigo_causal_cxx.hpp -- TESTIGO CAUSAL DEL PIPELINE (ARK-WIN-01-CXX)
// Reactor ARK / Canción del Micelio -- ARQUITECTO DE DATOS & COMANDANTE
// Reescrito por el Jefe de Arquitectos de Código (2026-11-01), B1 del acta AUDITORIA_CAUSAL/AUDITORIA_VENTANAS_CXX.md.
//
// El testigo anterior perturbaba x[c+1..] y solo alimentaba x[0..c]: el futuro perturbado no se leía nunca (con choque 1e12 y
// con choque 0 daba lo mismo), y auditaba paso(), que no puede ver el futuro por construcción. No veía cuatro de cuatro fugas.
//
// Definición de trabajo (CLAUDE.md del Jefe): X_t es causal si es F_t-medible, con F_t lo que llegó hasta t. Para cada corte c:
//   (a) truncado:     el pipeline sobre el CRUDO[0..c];
//   (b) perturbación: el pipeline sobre el crudo entero con las filas > c cambiadas (precios y huecos);
// y toda ventana emitida en una fila <= c tiene que ser bit a bit la misma en (a) y en (b), y emitirse en las dos.
// El pipeline es todo lo que va del crudo a la ventana (la preparación del que llama incluida): ahí vive la fuga.
// Mutante obligatorio: el mismo pipeline con el crudo adelantado una fila (la fila k ve la k+1). Si el testigo no lo caza, no
// tiene potencia sobre estos datos y su veredicto no vale (es_causal = false).
//
// El testigo reserva memoria (copias del crudo y una huella de 8 bytes por ventana); no corre en el ciclo caliente.
#pragma once

#include "ark/ventanas/contratos_ventanas.hpp"
#include "ark/ventanas/motor_ventanas.hpp"
#include "ark/ventanas/motor_multiactivo.hpp"
#include <cstdint>
#include <cstddef>
#include <cmath>
#include <cstring>
#include <vector>
#include <array>
#include <utility>
#include <memory>
#include <algorithm>
#include <limits>
#include <iostream>

namespace ark::ventanas {

// =======================================================================
// 1. Generador determinista (XorShift64) y huellas exactas
// =======================================================================

class RngDeterministicoCxx {
private:
    uint64_t state_{123456789ULL};

public:
    constexpr explicit RngDeterministicoCxx(uint64_t seed = 123456789ULL) noexcept
        : state_(seed != 0 ? seed : 123456789ULL) {}

    uint64_t siguiente_u64() noexcept {
        uint64_t x = state_;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        state_ = x;
        return x;
    }

    double siguiente_uniforme(double min_val, double max_val) noexcept {
        const double u01 = static_cast<double>(siguiente_u64() >> 11) * (1.0 / 9007199254740992.0);   // [0, 1)
        return min_val + u01 * (max_val - min_val);
    }
};

class HuellaFnv64 {
    uint64_t h_{1469598103934665603ULL};

public:
    void bytes(const void* p, size_t n) noexcept {
        const auto* b = static_cast<const unsigned char*>(p);
        for (size_t i = 0; i < n; ++i) {
            h_ ^= b[i];
            h_ *= 1099511628211ULL;
        }
    }
    template <class T>
    void valor(const T& v) noexcept { bytes(&v, sizeof(T)); }
    [[nodiscard]] uint64_t valor_final() const noexcept { return h_; }
};

// Huella de todo lo que una ventana mono-activo trae de los datos (bits exactos, NaN incluidos).
template <size_t W>
[[nodiscard]] uint64_t huella_ventana(const VentanaCausalCertificadaCxx<W>& v) noexcept {
    HuellaFnv64 h;
    h.valor(v.t_actual);
    h.valor(v.t_cierre_ms);
    h.valor(v.disponibilidad_declarada);
    h.bytes(v.valores.data(), sizeof(double) * W);
    h.bytes(v.mascara_huecos.data(), sizeof(bool) * W);
    h.bytes(v.con_sesion.data(), sizeof(bool) * W);
    h.valor(v.reloj_acumulado);
    h.valor(v.delta_reloj_ventana);
    h.valor(v.es_valida);
    h.valor(v.ratio_imputacion);
    h.valor(v.barras_omitidas);
    return h.valor_final();
}

template <size_t W, size_t M>
[[nodiscard]] uint64_t huella_ventana(const VentanaCausalCertificadaMultiactivoCxx<W, M>& v) noexcept {
    HuellaFnv64 h;
    h.valor(v.t_actual);
    h.valor(v.t_cierre_ms);
    h.valor(v.disponibilidad_declarada);
    for (size_t t = 0; t < W; ++t) {
        h.bytes(v.valores[t].data(), sizeof(double) * M);
        h.bytes(v.mascara_huecos[t].data(), sizeof(bool) * M);
        h.bytes(v.market_open_mask[t].data(), sizeof(bool) * M);
    }
    h.bytes(v.delta_reloj_por_activo.data(), sizeof(double) * M);
    h.bytes(v.reloj_acum_por_activo.data(), sizeof(double) * M);
    h.valor(v.es_valida);
    h.valor(v.ratio_imputacion_global);
    h.valor(v.barras_omitidas);
    for (size_t u = 0; u < M; ++u) {
        h.bytes(v.clasificacion_pares[u].data(), sizeof(ClasePar) * M);
        h.bytes(v.margen_pares_s[u].data(), sizeof(double) * M);
    }
    return h.valor_final();
}

// =======================================================================
// 2. Crudo auditable: filas que se pueden truncar, perturbar y adelantar
// =======================================================================

namespace detalle_testigo {
inline void perturbar_entrada(double& valor, bool& con_sesion_puede_nan, RngDeterministicoCxx& rng) noexcept {
    if (std::isfinite(valor)) valor *= 1.0 + rng.siguiente_uniforme(-0.5, 0.5);
    if (con_sesion_puede_nan && rng.siguiente_uniforme(0.0, 1.0) < 0.25) valor = std::numeric_limits<double>::quiet_NaN();
}
}  // namespace detalle_testigo

struct SerieCrudaMonoCxx {
    std::vector<EntradaPasoCxx> filas;

    [[nodiscard]] size_t n() const noexcept { return filas.size(); }
    [[nodiscard]] SerieCrudaMonoCxx truncado(size_t c) const {
        return {std::vector<EntradaPasoCxx>(filas.begin(), filas.begin() + static_cast<std::ptrdiff_t>(std::min(c + 1, n())))};
    }
    [[nodiscard]] SerieCrudaMonoCxx perturbado(size_t c, uint64_t semilla) const {
        SerieCrudaMonoCxx o = *this;
        RngDeterministicoCxx rng(semilla);
        for (size_t k = c + 1; k < o.n(); ++k) {
            bool puede = o.filas[k].con_sesion;
            detalle_testigo::perturbar_entrada(o.filas[k].valor, puede, rng);
            o.filas[k].volumen *= 1.0 + rng.siguiente_uniforme(0.0, 1.0);
        }
        return o;
    }
    // La fila k lleva los datos (no el sello) de la k+1: la fuga de un paso del mutante obligatorio.
    [[nodiscard]] SerieCrudaMonoCxx adelantado() const {
        SerieCrudaMonoCxx o = *this;
        for (size_t k = 0; k + 1 < o.n(); ++k) {
            o.filas[k].valor = filas[k + 1].valor;
            o.filas[k].volumen = filas[k + 1].volumen;
        }
        return o;
    }
};

template <size_t M>
struct SerieCrudaMultiCxx {
    std::vector<uint64_t> ts;
    std::vector<std::array<EntradaActivoCxx, M>> filas;
    std::vector<int32_t> sesiones_omitidas;   // vacío = -1 en todas

    [[nodiscard]] size_t n() const noexcept { return filas.size(); }
    [[nodiscard]] SerieCrudaMultiCxx truncado(size_t c) const {
        SerieCrudaMultiCxx o;
        const auto e = static_cast<std::ptrdiff_t>(std::min(c + 1, n()));
        o.ts.assign(ts.begin(), ts.begin() + e);
        o.filas.assign(filas.begin(), filas.begin() + e);
        if (!sesiones_omitidas.empty()) o.sesiones_omitidas.assign(sesiones_omitidas.begin(), sesiones_omitidas.begin() + e);
        return o;
    }
    [[nodiscard]] SerieCrudaMultiCxx perturbado(size_t c, uint64_t semilla) const {
        SerieCrudaMultiCxx o = *this;
        RngDeterministicoCxx rng(semilla);
        for (size_t k = c + 1; k < o.n(); ++k) {
            for (size_t m = 0; m < M; ++m) {
                bool puede = o.filas[k][m].con_sesion;
                detalle_testigo::perturbar_entrada(o.filas[k][m].valor, puede, rng);
            }
        }
        return o;
    }
    [[nodiscard]] SerieCrudaMultiCxx adelantado() const {
        SerieCrudaMultiCxx o = *this;
        for (size_t k = 0; k + 1 < o.n(); ++k)
            for (size_t m = 0; m < M; ++m) o.filas[k][m].valor = filas[k + 1][m].valor;
        return o;
    }
};

// =======================================================================
// 3. Pipelines estándar: preparación del que llama + motor
// =======================================================================

using HuellasEmitidas = std::vector<std::pair<size_t, uint64_t>>;   // (fila del crudo en la que se emitió, huella)

// preparar: SerieCrudaMonoCxx -> SerieCrudaMonoCxx con el mismo número de filas (lo que el que llama hace antes del motor).
template <size_t W, class Preparar>
[[nodiscard]] auto pipeline_motor_mono(const ConfiguracionVentanaCausalCxx<W>& cfg, Preparar preparar) {
    return [cfg, preparar](const SerieCrudaMonoCxx& crudo) {
        const SerieCrudaMonoCxx prep = preparar(crudo);
        HuellasEmitidas out;
        auto motor = std::make_unique<MotorVentanasCausalesCxx<W>>(cfg);
        auto win = std::make_unique<VentanaCausalCertificadaCxx<W>>();
        for (size_t k = 0; k < prep.n(); ++k)
            if (motor->paso(prep.filas[k], *win) == ResultadoPaso::VENTANA_EMITIDA) out.emplace_back(k, huella_ventana<W>(*win));
        return out;
    };
}

template <size_t W, size_t M, class Preparar>
[[nodiscard]] auto pipeline_motor_multi(const ConfiguracionVentanaCausalCxx<W>& cfg, const std::array<double, M>& ticks,
                                        Preparar preparar) {
    return [cfg, ticks, preparar](const SerieCrudaMultiCxx<M>& crudo) {
        const SerieCrudaMultiCxx<M> prep = preparar(crudo);
        HuellasEmitidas out;
        auto motor = std::make_unique<MotorVentanasMultiactivoCxx<W, M>>(cfg);
        for (size_t m = 0; m < M; ++m) motor->establecer_tick(m, ticks[m]);
        auto win = std::make_unique<VentanaCausalCertificadaMultiactivoCxx<W, M>>();
        for (size_t k = 0; k < prep.n(); ++k) {
            const int32_t so = prep.sesiones_omitidas.empty() ? -1 : prep.sesiones_omitidas[k];
            if (motor->paso(prep.ts[k], prep.filas[k].data(), *win, so) == ResultadoPaso::VENTANA_EMITIDA)
                out.emplace_back(k, huella_ventana<W, M>(*win));
        }
        return out;
    };
}

struct IdentidadCxx {
    template <class D>
    D operator()(const D& d) const { return d; }
};

// =======================================================================
// 4. Testigo
// =======================================================================

struct ReporteTestigoPipelineCxx {
    size_t cortes_evaluados{0};
    size_t ventanas_comparadas{0};
    size_t ventanas_distintas{0};        // misma fila, huella distinta
    size_t emisiones_desparejas{0};      // emitida en (a) y no en (b), o al revés, en filas <= c
    size_t primer_corte_con_fuga{std::numeric_limits<size_t>::max()};
    size_t primera_fila_con_fuga{std::numeric_limits<size_t>::max()};
    bool sin_discrepancias{true};        // el pipeline auditado pasó (a) == (b) en todos los cortes
    bool mutante_cazado{false};          // el testigo vio la fuga de un paso: tiene potencia
    bool es_causal{false};               // sin_discrepancias Y mutante_cazado

    void imprimir(std::ostream& os = std::cout) const {
        os << "  testigo del pipeline: cortes=" << cortes_evaluados << " ventanas comparadas=" << ventanas_comparadas
           << " distintas=" << ventanas_distintas << " desparejas=" << emisiones_desparejas
           << " mutante cazado=" << (mutante_cazado ? "si" : "NO") << " -> " << (es_causal ? "CAUSAL" : "FUGA o SIN POTENCIA");
        if (!sin_discrepancias) os << " (primer corte con fuga " << primer_corte_con_fuga << ", fila " << primera_fila_con_fuga << ")";
        os << "\n";
    }
};

class TestigoCausalPipelineCxx {
public:
    // Solo el ensayo (a) == (b), sin el mutante. Cortes c = primer_corte, primer_corte + paso, ..., n - 2.
    template <class Datos, class Pipeline>
    static ReporteTestigoPipelineCxx comparar(const Datos& crudo, Pipeline&& pipe, size_t primer_corte, size_t paso_cortes,
                                              uint64_t semilla = 20261101ULL) {
        ReporteTestigoPipelineCxx r{};
        const size_t n = crudo.n();
        if (paso_cortes == 0) paso_cortes = 1;
        for (size_t c = primer_corte; c + 1 < n; c += paso_cortes) {
            const HuellasEmitidas a = pipe(crudo.truncado(c));
            const HuellasEmitidas b = pipe(crudo.perturbado(c, semilla + c));
            ++r.cortes_evaluados;
            size_t ia = 0, ib = 0;
            while (ia < a.size() || (ib < b.size() && b[ib].first <= c)) {
                const bool hay_a = ia < a.size();
                const bool hay_b = ib < b.size() && b[ib].first <= c;
                size_t fila = 0;
                bool mal = false;
                if (hay_a && hay_b && a[ia].first == b[ib].first) {
                    fila = a[ia].first;
                    ++r.ventanas_comparadas;
                    if (a[ia].second != b[ib].second) { ++r.ventanas_distintas; mal = true; }
                    ++ia;
                    ++ib;
                } else if (hay_a && (!hay_b || a[ia].first < b[ib].first)) {
                    fila = a[ia].first;
                    ++r.emisiones_desparejas;
                    mal = true;
                    ++ia;
                } else {
                    fila = b[ib].first;
                    ++r.emisiones_desparejas;
                    mal = true;
                    ++ib;
                }
                if (mal && r.sin_discrepancias) {
                    r.sin_discrepancias = false;
                    r.primer_corte_con_fuga = c;
                    r.primera_fila_con_fuga = fila;
                }
            }
        }
        if (r.cortes_evaluados == 0) r.sin_discrepancias = false;   // sin cortes no se certifica nada
        return r;
    }

    // Auditoría completa: el pipeline y su mutante de un paso.
    template <class Datos, class Pipeline>
    static ReporteTestigoPipelineCxx auditar(const Datos& crudo, Pipeline&& pipe, size_t primer_corte, size_t paso_cortes,
                                             uint64_t semilla = 20261101ULL) {
        ReporteTestigoPipelineCxx r = comparar(crudo, pipe, primer_corte, paso_cortes, semilla);
        auto mutante = [&pipe](const Datos& d) { return pipe(d.adelantado()); };
        const ReporteTestigoPipelineCxx rm = comparar(crudo, mutante, primer_corte, paso_cortes, semilla);
        r.mutante_cazado = !rm.sin_discrepancias;
        r.es_causal = r.sin_discrepancias && r.mutante_cazado;
        return r;
    }
};

// =======================================================================
// 5. Compatibilidad: la interfaz anterior, con el ensayo real
// =======================================================================

struct ReporteTestigoCausalCxx {
    double max_diff{0.0};
    double max_diff_perturbacion{0.0};
    double max_diff_truncado{0.0};
    size_t cortes_evaluados{0};
    size_t corte_peor_diff{0};
    bool mutante_cazado{false};
    bool es_causal_estricto{false};

    void imprimir(std::ostream& os = std::cout) const {
        os << "  testigo (interfaz anterior, ensayo real): cortes=" << cortes_evaluados << " max_diff=" << max_diff
           << " mutante cazado=" << (mutante_cazado ? "si" : "NO") << " -> " << (es_causal_estricto ? "CAUSAL" : "FUGA o SIN POTENCIA")
           << "\n";
    }
};

// Audita el motor sobre arrays planos (sin preparación). Para auditar una preparación, TestigoCausalPipelineCxx.
template <size_t W>
class TestigoCausalAntiLookaheadCxx {
public:
    static ReporteTestigoCausalCxx auditar(const uint64_t* timestamps_ms, const double* precios, const bool* mascaras, size_t N,
                                           MotorVentanasCausalesCxx<W>& motor, double /*magnitud_shock*/ = 1e12,
                                           bool /*verificar_truncado*/ = true, bool /*perturbacion_aleatoria*/ = true,
                                           bool modo_exhaustivo = true, size_t max_cortes_muestreo = 50) {
        ReporteTestigoCausalCxx rep{};
        if (!timestamps_ms || !precios || N < W + 1) return rep;
        SerieCrudaMonoCxx crudo;
        crudo.filas.resize(N);
        for (size_t k = 0; k < N; ++k) {
            crudo.filas[k].ts_ms = timestamps_ms[k];
            crudo.filas[k].valor = precios[k];
            crudo.filas[k].mascara = mascaras ? mascaras[k] : false;
        }
        const size_t paso = modo_exhaustivo ? 1 : std::max<size_t>(1, (N - W) / std::max<size_t>(1, max_cortes_muestreo));
        const auto pipe = pipeline_motor_mono<W>(motor.configuracion(), IdentidadCxx{});
        const ReporteTestigoPipelineCxx r = TestigoCausalPipelineCxx::auditar(crudo, pipe, W - 1, paso);
        rep.cortes_evaluados = r.cortes_evaluados;
        rep.max_diff = r.sin_discrepancias ? 0.0 : 1.0;
        rep.max_diff_perturbacion = rep.max_diff;
        rep.max_diff_truncado = rep.max_diff;
        rep.corte_peor_diff = r.sin_discrepancias ? 0 : r.primer_corte_con_fuga;
        rep.mutante_cazado = r.mutante_cazado;
        rep.es_causal_estricto = r.es_causal;
        motor.registrar_auditoria_testigo(rep.max_diff);
        return rep;
    }
};

} // namespace ark::ventanas
