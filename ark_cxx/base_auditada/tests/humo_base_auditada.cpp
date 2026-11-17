// humo_base_auditada.cpp -- toda la base auditada en UNA unidad de traduccion: compila sin choques de nombres entre bloques.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04. Uso: compilar_humo.bat (rc 0 = compila y corre).
#include "ark/dataset/adaptadores_fuentes.hpp"
#include "ark/dataset/adaptadores_raw.hpp"
#include "ark/dataset/calendario.hpp"
#include "ark/dataset/control_calidad.hpp"
#include "ark/dataset/deconvolucion_tv.hpp"
#include "ark/dataset/errores_py.hpp"
#include "ark/dataset/extractor_cxx.hpp"
#include "ark/dataset/formato_npz.hpp"
#include "ark/dataset/formato_parquet.hpp"
#include "ark/dataset/guardas_causales.hpp"
#include "ark/dataset/libro_hashes.hpp"
#include "ark/dataset/motor_ingesta.hpp"
#include "ark/dataset/motor_ventanas.hpp"
#include "ark/dataset/numpy_formato.hpp"
#include "ark/dataset/numpy_rng.hpp"
#include "ark/dataset/oraculo_lazo.hpp"
#include "ark/dataset/sincronizador_multiactivo.hpp"
#include "ark/dataset/supervivencia.hpp"
#include "ark/dataset/tabla_ingesta.hpp"
#include "ark/dataset/tablas_lazo_ventanas.hpp"
#include "ark/dataset/tipos_dataset.hpp"
#include "ark/dataset/zona_horaria.hpp"
#include "ark/ventanas/contratos_ventanas.hpp"
#include "ark/ventanas/desacoplamiento_artefactos.hpp"
#include "ark/ventanas/lazo_shannon.hpp"
#include "ark/ventanas/lazo_tablas.hpp"
#include "ark/ventanas/motor_multiactivo.hpp"
#include "ark/ventanas/motor_ventanas.hpp"
#include "ark/ventanas/motor_reloj_intrinseco.hpp"
#include "ark/ventanas/flujo_barras.hpp"
#include "ark/dataset/fase_calendario.hpp"
#include "ark/fractal/observador_fractal_filas.hpp"
#include "ark/metrologia/observador_saltos_filas.hpp"
#include "ark/metrologia/detector_deslizante.hpp"
#include "ark/metrologia/observador_delta_fractal.hpp"
#include "ark/ventanas/testigo_causal_cxx.hpp"
#include "ark/holm.hpp"
#include "ark/fractal/observador_fractal_multiescala.hpp"
#include "ark/metrologia/observador_memoria.hpp"
#include "ark/metrologia/superestadistica_inversa.hpp"
#include "ark/coordinacion/observador_coordinacion.hpp"
#include "ark/metrologia/calibrador_familia.hpp"
#include "ark/metrologia/observador_saltos_familia.hpp"
#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/lilyth/modelo_lineal_nulo.hpp"
#include "ark/lilyth/nulo_lilyth.hpp"
#include "ark/lilyth/armonico_lilyth.hpp"
#include "ark/lilyth/fin_calor_lilyth.hpp"
#include "ark/lilyth/reloj_marea_lilyth.hpp"

#include <cstdio>

int main() {
    ark::ventanas::MotorVentanasCausalesCxx<32> motor{};
    ark::fractal::SensorDFACausal<50, 2> dfa(ark::fractal::Integracion::FLUJO);
    ark::metrologia::SensorWhittleExacto<256> elw(ark::metrologia::ConfiguracionSensorELW{128, 30, 1, ark::metrologia::Integracion::FLUJO, -0.5, 2.0});
    (void)motor; (void)dfa; (void)elw;
    std::printf("humo base auditada: compila y enlaza\n");
    return 0;
}
// instanciacion de las lentes y el conjunto (2026-11-05)
template class ark::metrologia::ObservadorMemoriaDFA<50, 2>;
template class ark::metrologia::ObservadorDeltaFractal<ark::metrologia::ObservadorMemoria<1024>>;
template class ark::metrologia::ObservadorDeltaFractal<ark::metrologia::ObservadorMemoriaDFA<50, 2>>;
// operador de coordinacion de pares (2026-11-05)
template class ark::coordinacion::SensorCoordinacion<240>;
template class ark::coordinacion::ObservadorCoordinacion<240>;
template class ark::metrologia::ObservadorSaltosFamilia<>;
// operador de Lilyth, etapas 1-6 (2026-11-05; etapas 2-6 de Argon corregidas por Code el 2026-11-07)
template class ark::lilyth::ObservadorLilyth<ark::metrologia::ObservadorSaltos<>>;
template class ark::lilyth::detalle_nulo::EspejoPurgaNulo<ark::lilyth::RasgosDelta<ark::metrologia::ObservadorSaltos<>>::b_max>;
