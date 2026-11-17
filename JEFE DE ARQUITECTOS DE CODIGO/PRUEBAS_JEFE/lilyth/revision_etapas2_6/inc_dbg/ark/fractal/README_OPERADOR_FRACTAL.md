# Operador fractal -- `ObservadorFractalUniversal` y `PiramideFractal` (ARK-FRAC-01)

Codigo: `observador_fractal_multiescala.hpp`, `sensor_dfa_causal.hpp`, `regularizador_dual_fractal.hpp`, `tipos_fractal.hpp`
(`nucleo_observador_universal.hpp` es un borrador sin uso). Acta: `JEFE DE ARQUITECTOS DE CODIGO/AUDITORIA_CAUSAL/BLOQUE_FRACTAL_UNIVERSAL_CXX.md`;
viabilidad del sensor: `AUDITORIA_CAUSAL/MANUAL_OPERADOR_FRACTAL_VIABILIDAD.md`. Es telemetria (declaracion metrologica de EL ABISMO:
prohibido usar H en escalas < 1 h como gatillo taker).

## Que hace

Cualquier senal causal x_t -> sensor DFA-1 solapado (W = 50, paso S = 2) -> H^ -> innovacion nu = (H^ - mu0)/s0 -> CUSUM bilateral
(k = Delta/2, h por el ARL0 del nulo) -> regularizador dual del Fisico en el borde presente: L2 (Tikhonov, lambda de Morozov) o L1
(TV de Condat, lambda_TV = cuantil 1 - alfa del maximo del puente de sumas parciales en el nulo), elegida por el certificado del puente
(decision del Comandante del 2026-11-03). Exporta tau50 medido por inyeccion. `PiramideFractal` corre seis observadores 1m, 5m, 15m,
1h, 4h y D1, cada uno con su calibracion.

## Entradas obligatorias

`ConfiguracionFractal.arl0_objetivo_barras`; `integracion` (FLUJO: H = alfa del DFA; NIVEL: H = alfa - 1); `agregacion` de la piramide
(ULTIMO niveles, SUMA flujos, MEDIA intensivos).

## Nulos que acepta (catalogo en `../metrologia/README_NULOS.md`)

| `nulo` | que es | cuando |
|---|---|---|
| `HISTORICO` (= NINGUNO) | la historia causal del propio observable: mu0, s0, lambda_TV y h miden su reposo tal como es | observables cuyo reposo no es H 0.5 (volumen, profundidad) o cuando la pregunta es "cambio respecto de lo habitual" |
| `SIGNOS` | los incrementos del observable con signo aleatorio (`replicas_nulo_signos` mundos): H 0.5 con la marginal y la heterocedasticidad de X | precios y flujos simetricos (el nulo real de EL ABISMO) |

La ganancia del sensor (Delta) y tau50 se miden SIEMPRE en el mundo de signos, sea cual sea `nulo`. No ofrece PERMUTACION ni BLOQUES:
para un observable positivo (volumen, horquilla) el nulo sin memoria correcto seria PERMUTACION del log, que este bloque no implementa
(limitacion declarada; el operador de memoria si la tiene).

## Limites medidos

- **El sensor DFA W = 50 no es viable en ninguna escala** (19 series, dos nulos): s0 ~ 0.15-0.2 por barra; un cambio de H de 0.15 tarda
  ~12 h en verse a 1 min (`MEDIDO_DELTA_FRACTAL.md`). Para un sensor con menos ruido, ver `../metrologia/README_OPERADOR_MEMORIA.md`.
- Conexion (2026-11-04): vive en la base auditada y compila solo contra ella.
- **h del CUSUM (2026-11-05, NO medido):** `h_por_arl` lee el ARL en el nulo con >= 30 alarmas y, mas alla, extrapola en recta
  sobre log ARL (mitad superior del rango legible), como `buscar_h_con` del operador delta. En el operador delta esa
  extrapolacion resulto permisiva con colas pesadas y mucha distancia (x7-14 a ~140 veces). Para el CUSUM la forma teorica
  (ARL ~ e^{2kh}) es casi lineal, pero las colas de H^ no estan medidas. Ver
  `../ventanas/README_VENTANAS_Y_DETECTORES.md` seccion 5 y `mensajes/CONTEXTO_AUDITOR_EXTERNO_CALIBRACION_H.md`.

## Entrada por filas en reloj intrinseco (2026-11-04, orden del Comandante)

`ObservadorFractalFilas<W, S, N_REG, W_FLUJO>` (`observador_fractal_filas.hpp`): filas -> aduana del motor de pared auditado -> barra
intrinseca (`ark/ventanas/flujo_barras.hpp`, la fuente comun `FuenteFilasCxx`: reloj, delta_nu con arrastre, agregacion) ->
`ObservadorFractalUniversal`, que no se toca. `calibrar_filas(cfg, fuente, filas)` y `paso_fila(fila)`.

- Una barra con alguna fila imputada entra como NaN en el ciclo caliente: las W ventanas del DFA que la contienen no se estiman y el
  regularizador corta su ventana. En la calibracion entra con el valor que le dio la aduana (arrastre en ULTIMO, 0 de aportacion en
  SUMA), porque el nulo reconstruye el nivel sumando incrementos y un NaN lo anularia; se reporta `imputadas_calibracion()`.
- Controles (`PRUEBAS_JEFE/conexion_filas/`): G3 en PARED con delta_nu 1 la calibracion y todo el estado son bit a bit los del
  observador alimentado con los valores; G4 en NU_VOLUMEN con imputadas y cerradas = el observador sobre la serie de barras; G7a cero
  reservas; G8a causalidad (adulterar el futuro no cambia lo anterior y si lo posterior).
- La piramide 1m..D1 (`PiramideFractal`) sigue con su agregador de calendario en segundos y NO tiene camino por filas: unificar la
  agregacion multiescala con el reloj intrinseco esta por decidir (ver `ark_cxx/base_auditada/README.md`).

## El DFA como lente de memoria del operador delta (2026-11-05, orden del Comandante: "no se desecha")

`ark/metrologia/observador_memoria_dfa.hpp` (`ObservadorMemoriaDFA<W = 50, S = 2>`): el sensor DFA-1 de este bloque, sin tocar, alimenta el
operador delta (DISJUNTO o DESLIZANTE) sobre su H^, con la misma tuberia que el ELW (nulo_sensor y nulos_contraste de entrada; umbrales
autocalibrados sobre ese H^; filas por la fuente comun y el motor de reloj intrinseco con W = 50). Es intercambiable con el ELW en el
observador conjunto (`observador_delta_fractal.hpp`). Medido en sinteticos (`PRUEBAS_JEFE/delta_fractal/MEDIDO_DELTA_FRACTAL_CONJUNTO.md`):
detecta un cambio de H 0.5 -> 0.75 en 8/8 con ~380 barras de retardo (el ELW en 205), porque su H^ tiene el doble de ruido (sd 0.151 contra
0.077); a cambio, ventana de 50 contra 184. **Pendiente: calibrarlo sobre observables que no sean retornos** (volumen, horquilla, OFI), con
datos reales.
