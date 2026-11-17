# Base auditada C++20 del pipeline nuevo (`ark_cxx/base_auditada/`)

Orden del Comandante (2026-11-04): el pipeline nuevo (dataset, ventanas, Holm-Bonferroni, operador delta y operador fractal) se
conecta **solo** a la base corregida; la base se duplica para el Arquitecto de Datos y el Fisico de Cabecera; las versiones no
corregidas (Python o C++ defectuoso) se archivan en la carpeta de cada agente y no se usan.

**Esta carpeta es la fuente de verdad.** Las copias de los agentes (`ARQUITECTO DE DATOS/nucleo_cxx_auditado/`,
`FISICO DE CABECERA/nucleo_cxx_auditado/`) se verifican contra ella:

```
python -I ark_cxx/base_auditada/verificar_base.py --verificar "<carpeta de la copia>"
```

Una copia que no verifica no es la base: se recopia o se archiva. Nadie edita una copia; los cambios entran aqui, con su prueba, y se
regenera el manifiesto (`verificar_base.py --generar`).

**El ARK 1.0 (`ark_cxx/src/ark.cpp`, el binario que corre en papel) NO usa esta base**: es el pipeline antiguo y no se toca (orden
del Comandante). Sigue compilando contra `ark_cxx/include/`, donde quedan la capa de datos de Flash, las ventanas de Antigravity y la
coordinacion SIN corregir porque el ARK 1.0 depende de ellas; esas carpetas llevan la marca `NO_AUDITADO.md`.

## Contenido y procedencia

| bloque | ruta | de donde sale | acta y bateria |
|---|---|---|---|
| capa de datos (Flash, corregida) + porte de la ingesta y las ventanas del dataset (D1-D6 corregidos) | `include/ark/dataset/` | `JEFE.../PRUEBAS_JEFE/dataset_cxx/include/ark/dataset/` | `AUDITORIA_CAUSAL/CORRECCION_DATASET_CXX.md`, `AUDITORIA_CAUSAL/PORTE_INGESTA_Y_VENTANAS_CXX.md` |
| ventanas rodantes streaming (Antigravity, corregidas) + lazo de Shannon nuevo (prior, gauss/Laplace, phi) | `include/ark/ventanas/` | `PRUEBAS_JEFE/ventanas_cxx/corregido/include/ark/ventanas/` | `AUDITORIA_CAUSAL/CORRECCION_VENTANAS_CXX.md`, `CORRECCION_SHANNON_LAPLACE_CXX.md` |
| motor de ventanas en reloj intrinseco (W y S de ejecucion, arrastre del exceso) y fuente comun de filas / flujo de barras | `include/ark/ventanas/motor_reloj_intrinseco.hpp`, `flujo_barras.hpp` | Code (2026-11-04) | `PRUEBAS_JEFE/ventanas_intrinseco/` (16/16), `PRUEBAS_JEFE/conexion_filas/`, `PRUEBAS_JEFE/metrologia/prueba_whittle.cpp` modo `filas` |
| fase de la purga desde el calendario (hora local, feriados) | `include/ark/dataset/fase_calendario.hpp` | Code (2026-11-04) | `PRUEBAS_JEFE/conexion_filas/` (G1) |
| Holm-Bonferroni (ARK-STAT-01) | `include/ark/holm.hpp` | `PRUEBAS_JEFE/holm_cxx/include/ark/holm.hpp` | `AUDITORIA_CAUSAL/CORRECCION_HOLM_CXX.md` |
| operador fractal (ARK-FRAC-01) | `include/ark/fractal/` | Code | `AUDITORIA_CAUSAL/BLOQUE_FRACTAL_UNIVERSAL_CXX.md`; `PRUEBAS_JEFE/bloque_fractal_universal/` |
| operador delta y temperatura, superestadistica, sensor de memoria ELW (ARK-METRO-01) | `include/ark/metrologia/` | Code | `PRUEBAS_JEFE/metrologia/MEDIDO.md`, `MEDIDO_REAL.md`, `MEDIDO_WHITTLE.md` |
| operador de coordinacion de pares (TE, MI y sinergia neta aparte, PIT aleatorizado, nulos SIGNOS y ROTACION, operador delta por componente) | `include/ark/coordinacion/` | Code (2026-11-05) | `PRUEBAS_JEFE/coordinacion_sensor/` |
| operador de Lilyth, ETAPAS 1-6 (restitucion tras un evento del operador delta; nulo lineal de Onsager, vias (a) y (b) con guarda de solape, espejo exacto de la purga, armonico emparejado y Holm de 4, fin del calor, marea descriptiva; Etapa 7 en espera) | `include/ark/lilyth/` | Argon; Code reviso, corrigio y conecto (Etapa 1 2026-11-05, Etapas 2-6 2026-11-07); contraste del resumen bifurcado por Code (POOL_LOO por omision / INTERCAMBIABLE, 2026-11-15); pool derivado (M_MIN por omision / DERIVADO) y mascara de fuente de los nulos (2026-11-16) | `PRUEBAS_JEFE/lilyth/revision_etapa1/`, `PRUEBAS_JEFE/lilyth/revision_etapas2_6/`, `PRUEBAS_JEFE/lilyth/bifurcacion_intercambiable/`, `PRUEBAS_JEFE/lilyth/pool_y_fuente_nulos/` |

Fuera de la base, a proposito: el operador de coordinacion HEREDADO (`ark_cxx/include/ark/coordinacion/`, defectos medidos el
11-03; el nuevo, `include/ark/coordinacion/`, si es de la base), el borrador `nucleo_observador_universal.hpp` (sin uso) y todo el ARK 1.0 (`alertas`, `regla_e`, `kelly`, ...).

Comprobado el 2026-11-04: toda la base compila en una sola unidad con `/W4` sin avisos (`tests/compilar_humo.bat`), y las baterias
del fractal, de la metrologia y del ELW pasan compilando desde aqui.

## Quien hace que (para no ser redundantes)

| capa | responsable unico | lo que los operadores NO deben rehacer |
|---|---|---|
| lectura de fuentes, sellos R1 de disponibilidad, calendarios y feriados, supervivencia, guardas causales, fase de la purga | `ark/dataset/` (ingesta; `fase_calendario.hpp`) | los operadores por filas no arman la fase: sale del calendario del activo en hora local (`FuenteFilasCxx.fase`) |
| validacion causal fila a fila (monotonia, imputacion declarada, presupuesto, huecos por calendario), reloj acumulado nu (PARED, NU_TRADES, NU_VOLUMEN, NU_NOCIONAL, NU_INFO_SHANNON), barras intrinsecas | `ark/ventanas/` (motor streaming, sin memoria dinamica; `FuenteFilasCxx` y `FlujoBarrasIntrinsecasCxx` comunes a los tres operadores) | los operadores no validan ni imputan: consumen la barra ya certificada (valor, imputada, nu, cierre, fase) |
| familias de pruebas multiples | `ark/holm.hpp` | nadie reimplementa Holm |
| estimacion de H, saltos, temperatura, superestadistica | `ark/fractal/`, `ark/metrologia/` | su anillo interno de W muestras queda solo como ciclo caliente (W derivada en tiempo de ejecucion) |
| agregacion multiescala | UNA sola (pendiente de unificar: hoy existen la escalera 2^j de metrologia, la piramide 1m..D1 del fractal y las escalas de divisores del motor) | -- |

## Estado de la conexion (2026-11-04)

Hecho: los operadores ya compilan y se prueban SOLO contra esta base (no contra `ark_cxx/include/`).

Hecho el 2026-11-04 (orden del Comandante: "actualiza el motor a reloj intrinseco y haz que su W no sea fija, sino un input del
operador"):
- `ark/ventanas/motor_reloj_intrinseco.hpp` (nuevo): usa el motor de pared auditado como aduana fila a fila (cambio aditivo en
  `motor_ventanas.hpp`: el accesor `ultima_fila()`, sin tocar ninguna regla; bateria de ventanas 504/504, monticulo, acople y oraculo del lazo 26/26, mutantes 44/44, paridad con el
  Python canonico identica bit a bit) y cierra una barra cada vez que el nu acumulado llega a `delta_nu`. W y S son entradas de
  ejecucion (W <= W_MAX).
- `ObservadorMemoria` (ELW) consume filas: `calibrar_filas` y `paso_fila`, con la W del motor = la W que deriva su banda. Detalle,
  controles y limites en `include/ark/metrologia/README_OPERADOR_MEMORIA.md`.

Hecho el 2026-11-04 (orden del Comandante: arrastre del exceso y "conectelos"):
- Arrastre: las fronteras de barra son los multiplos de `delta_nu` del reloj acumulado; el exceso pasa a la barra siguiente y las
  fronteras que una sola fila cruza de mas se cuentan (`fronteras_saltadas`). Controles T2 (referencia a mano) y T2b (conservacion
  del reloj: sum nu - delta_nu * fronteras = credito en [0, delta_nu)).
- `ark/ventanas/flujo_barras.hpp`: `FuenteFilasCxx` (aduana, reloj, delta_nu, agregacion, fase) comun a los tres operadores y
  `FlujoBarrasIntrinsecasCxx` (el motor con W = 1: barra a barra).
- `ark/dataset/fase_calendario.hpp`: la fase de la purga desde el calendario del activo y su hora local.
- `ark/fractal/observador_fractal_filas.hpp` y `ark/metrologia/observador_saltos_filas.hpp`: el fractal y el ObservadorSaltos sobre
  el crudo por filas, sin tocar los observadores; `ObservadorMemoria` usa la misma fuente y la fase de calendario.
- Controles en `PRUEBAS_JEFE/conexion_filas/` (G1-G9b, identidad bit a bit contra los observadores alimentados a mano, cero reservas,
  causalidad) y su arnes de mutantes (`mutaciones_conexion.py`: 15/15 cazados, testigos pasan).

Regla comun anotada (Comandante, 2026-11-04): ningun operador estima sobre datos imputados. El sensor de memoria no estima una
ventana con alguna barra imputada; el fractal recibe NaN en la barra imputada (sus W ventanas DFA no se estiman); el ObservadorSaltos
la recibe como no valida (un nivel de su escalera solo es valido si lo son sus dos hijas). El coste es medible: con 0.2 % de filas
imputadas, el 36 % de las ventanas de W ~ 170 barras del sensor de memoria queda sin estimar (F2). Los huecos de calendario (sesion
cerrada) no invalidan: congelan el reloj.

Hecho el 2026-11-05 (orden del Comandante: "haga el detector deslizante, pero no descarte el que va por ventanas independientes"):
- `ark/metrologia/detector_deslizante.hpp` y `ConfiguracionObservador.modo` (DISJUNTO por omision, identico bit a bit al de antes;
  DESLIZANTE: nivel en cada paso nativo, temperatura en cada barra). Reutiliza la escalera, `RingBufferEspejoCxx` y la busqueda de h
  (`detalle_nivel::buscar_h_con`, generalizada sin cambiar sus operaciones). Controles 31/31, mutantes 16/16 y lo medido en
  `JEFE.../PRUEBAS_JEFE/deslizante/` (`MEDIDO_DESLIZANTE.md`): en estos mundos no compra latencia; decisiones abiertas en el acta.

Hecho el 2026-11-05 (especificacion aprobada por el Comandante: "reutilizando las ventanas deslizantes e independientes, nulos
autocalibrados y el pipeline ya hecho"): `ark/metrologia/observador_memoria_dfa.hpp` (el DFA del fractal como lente de memoria) y
`observador_delta_fractal.hpp` (crudo + memoria ELW o DFA, ARL0 repartido por canales, etiqueta MEMORIA_PURA / CON_VOLATILIDAD). Guia de
ventanas y detectores con ventajas y desventajas: `include/ark/ventanas/README_VENTANAS_Y_DETECTORES.md`. Controles y lo medido en
`JEFE.../PRUEBAS_JEFE/delta_fractal/`.

Las dos agregaciones multiescala quedan (decision del Comandante, 2026-11-05). Hoy, sobre las barras intrinsecas, la escalera 2^j de saltos agrupa
por cuenta de barras y la piramide 1m..D1 del fractal agrupa por cubetas de calendario en segundos (sin camino por filas).
