// ============================================================================
// test_fase1_etiquetas_ventanas.cpp
// Verificación de la Fase 1:
//   - Banderas en BarraCausal (FLAG_MERCADO_QUIETO, FLAG_PUENTE_ITO,
//     FLAG_CHOQUE_EN_HUECO, FLAG_PAUSA_UMBRAL, FLAG_CORRUPCION_CAIDA)
//   - Propagación sin pérdida por EntradaPasoCxx -> MotorVentanasCausalesCxx
//     -> MotorVentanasIntrinsecoCxx -> FlujoBarrasIntrinsecasCxx -> BarraIntrinsecaCxx
//   - Clasificación exacta en clasificar_barra_intrinseca() (distinguiendo
//     MERCADO_QUIETO puro de CORRUPCION_O_CAIDA y preservando etiquetas en
//     agregación multi-fila)
// ============================================================================

#include "protocolo_inputs_lilyth.hpp"

#include <cstdio>
#include <type_traits>
#include <vector>

namespace {

using namespace ark::dataset;
using namespace ark::ventanas;
using namespace ark::lilyth::piloto;

int g_fallas = 0;

void verificar(const char* codigo, const char* desc, bool cond) {
    std::printf("[%s] %s: %s\n", cond ? "PASA" : "FALLA", codigo, desc);
    if (!cond) ++g_fallas;
}

BarraCausal crear_barra_base(uint64_t ts_ms, double precio) {
    BarraCausal b{};
    b.cierre_real_utc_ms = ts_ms;
    b.estampa_disponible_ms = ts_ms;
    b.fecha_sesion = 20261118;
    b.open = precio;
    b.high = precio;
    b.low = precio;
    b.close = precio;
    b.volume = 10.0;
    b.valor_nominal = precio * 10.0;
    b.fraccion_sesion = 1.0;
    b.es_sesion_con_dato = true;
    b.es_imputado = false;
    b.hay_dato = true;
    b.estado_supervivencia = EstadoSupervivencia::COTIZANDO;
    b.tipo_hueco = TipoHueco::NINGUNO;
    b.flags_control = 0;
    return b;
}

} // namespace

int main() {
    std::printf("====================================================================================================\n");
    std::printf(" FASE 1: VERIFICACION DE ETIQUETAS Y PROPAGACION POR ARK/DATASET -> ARK/VENTANAS\n");
    std::printf("====================================================================================================\n");

    // CF1.1: Ortogonalidad de los 8 bits de flags_control e invariancia de layout de BarraCausal
    {
        const uint8_t todos[8] = {
            FLAG_FILA_INSERTADA,
            FLAG_FILA_VACIA,
            FLAG_FRACCION_NAN,
            FLAG_MERCADO_QUIETO,
            FLAG_PUENTE_ITO,
            FLAG_CHOQUE_EN_HUECO,
            FLAG_PAUSA_UMBRAL,
            FLAG_CORRUPCION_CAIDA
        };
        bool disjuntos = true;
        uint8_t union_bits = 0;
        for (int i = 0; i < 8; ++i) {
            if (todos[i] == 0) disjuntos = false;
            for (int j = i + 1; j < 8; ++j) {
                if ((todos[i] & todos[j]) != 0u) disjuntos = false;
            }
            union_bits |= todos[i];
        }
        const bool layout_ok = std::is_trivially_copyable_v<BarraCausal> && (union_bits == 0xFFu);
        verificar("CF1.1", "8 bits de flags_control ortogonales (0xFF) y BarraCausal trivialmente copiable intacta",
                  disjuntos && layout_ok);
    }

    // CF1.2: Propagación 1-a-1 por FlujoBarrasIntrinsecasCxx (PARED, delta_nu = 1.0, congelar_en_huecos = false)
    {
        FuenteFilasCxx<64> fuente{};
        fuente.base.reloj = TipoRelojVentana::PARED;
        fuente.base.congelar_en_huecos = false;
        fuente.base.metodo_imputacion = MetodoImputacion::FORWARD_FILL;
        fuente.base.umbral_imputacion_maximo = 1.0;
        fuente.agregacion = AgregacionBarra::ULTIMO;

        FlujoBarrasIntrinsecasCxx<64> flujo{};
        const auto err = flujo.configurar(fuente, 1.0);

        // 1) Barra válida
        BarraCausal b_val = crear_barra_base(1000, 100.0);
        auto out_val = flujo.paso(a_entrada_paso_cxx(b_val, b_val.close));

        // 2) Barra de mercado quieto (es_imputado = true, precio arrastrado 100.0, sin FLAG_FILA_VACIA)
        BarraCausal b_quieto = crear_barra_base(2000, 100.0);
        b_quieto.es_imputado = true;
        b_quieto.hay_dato = false;
        b_quieto.flags_control = FLAG_FILA_INSERTADA;
        auto out_quieto = flujo.paso(a_entrada_paso_cxx(b_quieto, b_quieto.close));

        // 3) Barra corrupta (FLAG_FILA_VACIA / NaN)
        BarraCausal b_corr = crear_barra_base(3000, NAN_DATASET);
        b_corr.es_imputado = true;
        b_corr.hay_dato = false;
        b_corr.flags_control = FLAG_FILA_VACIA;
        auto out_corr = flujo.paso(a_entrada_paso_cxx(b_corr, b_corr.close));

        // 4) Barra desconexión explícita (TipoHueco::DESCONEXION con precio arrastrado)
        BarraCausal b_desc = crear_barra_base(4000, 100.0);
        b_desc.es_imputado = true;
        b_desc.hay_dato = false;
        b_desc.tipo_hueco = TipoHueco::DESCONEXION;
        auto out_desc = flujo.paso(a_entrada_paso_cxx(b_desc, b_desc.close));

        const bool ok_1a1 = (err == ErrorMotorIntrinseco::NINGUNO) &&
                            out_val.has_value() && !out_val->imputada &&
                            (clasificar_barra_intrinseca(*out_val) == ClasePasoLilyth::VALIDO) &&
                            out_quieto.has_value() && out_quieto->imputada && out_quieto->solo_mercado_quieto &&
                            (clasificar_barra_intrinseca(*out_quieto) == ClasePasoLilyth::MERCADO_QUIETO) &&
                            out_corr.has_value() && out_corr->imputada && !out_corr->solo_mercado_quieto &&
                            (clasificar_barra_intrinseca(*out_corr) == ClasePasoLilyth::CORRUPCION_O_CAIDA) &&
                            out_desc.has_value() && out_desc->imputada && !out_desc->solo_mercado_quieto &&
                            (clasificar_barra_intrinseca(*out_desc) == ClasePasoLilyth::CORRUPCION_O_CAIDA);

        verificar("CF1.2", "FlujoBarrasIntrinsecasCxx distingue VALIDO, MERCADO_QUIETO (solo_mercado_quieto=true) y CORRUPCION_O_CAIDA",
                  ok_1a1);
    }

    // CF1.3: Agregación multi-fila (delta_nu = 2.0):
    //   - Barra A (1 fila VALIDO + 1 fila MERCADO_QUIETO) -> solo_mercado_quieto = true (MERCADO_QUIETO)
    //   - Barra B (1 fila MERCADO_QUIETO + 1 fila CORRUPCION_O_CAIDA) -> solo_mercado_quieto = false (CORRUPCION_O_CAIDA)
    {
        FuenteFilasCxx<64> fuente{};
        fuente.base.reloj = TipoRelojVentana::PARED;
        fuente.base.congelar_en_huecos = false;
        fuente.base.metodo_imputacion = MetodoImputacion::FORWARD_FILL;
        fuente.base.umbral_imputacion_maximo = 1.0;
        fuente.agregacion = AgregacionBarra::ULTIMO;

        FlujoBarrasIntrinsecasCxx<64> flujo{};
        (void)flujo.configurar(fuente, 2.0);

        BarraCausal f1 = crear_barra_base(1000, 100.0);
        BarraCausal f2 = crear_barra_base(2000, 100.0);
        f2.es_imputado = true;
        f2.hay_dato = false;
        f2.flags_control = FLAG_FILA_INSERTADA;

        auto p1 = flujo.paso(a_entrada_paso_cxx(f1, f1.close));
        auto barra_a = flujo.paso(a_entrada_paso_cxx(f2, f2.close));

        BarraCausal f3 = crear_barra_base(3000, 100.0);
        f3.es_imputado = true;
        f3.hay_dato = false;
        f3.flags_control = FLAG_FILA_INSERTADA;

        BarraCausal f4 = crear_barra_base(4000, NAN_DATASET);
        f4.es_imputado = true;
        f4.hay_dato = false;
        f4.flags_control = FLAG_FILA_VACIA;

        auto p3 = flujo.paso(a_entrada_paso_cxx(f3, f3.close));
        auto barra_b = flujo.paso(a_entrada_paso_cxx(f4, f4.close));

        const bool ok_agreg = !p1.has_value() && barra_a.has_value() &&
                              barra_a->imputada && barra_a->solo_mercado_quieto &&
                              (clasificar_barra_intrinseca(*barra_a) == ClasePasoLilyth::MERCADO_QUIETO) &&
                              !p3.has_value() && barra_b.has_value() &&
                              barra_b->imputada && !barra_b->solo_mercado_quieto &&
                              (clasificar_barra_intrinseca(*barra_b) == ClasePasoLilyth::CORRUPCION_O_CAIDA);

        verificar("CF1.3", "En agregacion multi-fila (delta_nu=2), VALIDO+QUIETO preserva MERCADO_QUIETO y QUIETO+CORRUPTO fuerza CORRUPCION_O_CAIDA",
                  ok_agreg);
    }

    // CF1.4: Supervivencia de las etiquetas FLAG_PUENTE_ITO, FLAG_CHOQUE_EN_HUECO y FLAG_PAUSA_UMBRAL
    {
        FuenteFilasCxx<64> fuente{};
        fuente.base.reloj = TipoRelojVentana::PARED;
        fuente.base.congelar_en_huecos = false;
        fuente.base.metodo_imputacion = MetodoImputacion::FORWARD_FILL;
        fuente.base.umbral_imputacion_maximo = 1.0;
        fuente.agregacion = AgregacionBarra::ULTIMO;

        FlujoBarrasIntrinsecasCxx<64> flujo{};
        (void)flujo.configurar(fuente, 1.0);

        BarraCausal b_ito = crear_barra_base(1000, 125.0);
        b_ito.es_imputado = true;
        b_ito.flags_control = FLAG_PUENTE_ITO;
        auto out_ito = flujo.paso(a_entrada_paso_cxx(b_ito, b_ito.close));

        BarraCausal b_chq = crear_barra_base(2000, 150.0);
        b_chq.flags_control = FLAG_CHOQUE_EN_HUECO;
        auto out_chq = flujo.paso(a_entrada_paso_cxx(b_chq, b_chq.close));

        BarraCausal b_pau = crear_barra_base(3000, 150.0);
        b_pau.es_imputado = true;
        b_pau.flags_control = FLAG_PAUSA_UMBRAL;
        auto out_pau = flujo.paso(a_entrada_paso_cxx(b_pau, b_pau.close));

        const bool ok_flags = out_ito.has_value() && barra_tiene_puente_ito(*out_ito) && !barra_tiene_choque_en_hueco(*out_ito) &&
                              out_chq.has_value() && barra_tiene_choque_en_hueco(*out_chq) && !barra_tiene_puente_ito(*out_chq) &&
                              out_pau.has_value() && barra_tiene_pausa_umbral(*out_pau) &&
                              (clasificar_barra_intrinseca(*out_pau) == ClasePasoLilyth::CORRUPCION_O_CAIDA);

        verificar("CF1.4", "Etiquetas FLAG_PUENTE_ITO, FLAG_CHOQUE_EN_HUECO y FLAG_PAUSA_UMBRAL viajan intactas hasta BarraIntrinsecaCxx",
                  ok_flags);
    }

    std::printf("----------------------------------------------------------------------------------------------------\n");
    std::printf(" FALLAS FASE 1: %d\n", g_fallas);
    std::printf("====================================================================================================\n");
    return g_fallas == 0 ? 0 : 1;
}
