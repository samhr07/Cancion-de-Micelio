# Ventanas y detectores de la base auditada: que tipo usar, con sus ventajas y desventajas

Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-05, por orden del Comandante. Todo lo medido es SINTETICO (verificacion del instrumento);
con datos reales esta pendiente. Detalle y numeros: `JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/deslizante/MEDIDO_DESLIZANTE.md`.

Principio: no hay una ventana ni un detector "bala de plata". Cada combinacion responde mejor en una situacion; por eso se conservan
todas y se elige por configuracion.

## 1. Reloj de la barra (de que esta hecha cada barra)

| reloj | donde | ventajas | desventajas |
|---|---|---|---|
| PARED (tiempo de calendario) | `motor_ventanas.hpp` | barra = minuto/hora, comparable con cualquier fuente; la aduana (monotonia, imputacion declarada, huecos de calendario) es la base de todo | la actividad cambia dentro del dia: una barra de madrugada y una de apertura no pesan lo mismo |
| INTRINSECO (nu: trades, volumen, nocional, informacion) | `motor_reloj_intrinseco.hpp`, `flujo_barras.hpp` | cada barra lleva la misma cantidad de mercado; arrastre del exceso (las fronteras son multiplos de delta_nu y las saltadas se cuentan); W y S son entradas del operador | una fila grande salta fronteras (barras de ~1.3 filas con delta_nu = nu medio); con imputaciones, ningun operador estima la barra (regla comun) |

## 2. Ventanas del operador delta (saltos y temperatura): DISJUNTO o DESLIZANTE

`ConfiguracionObservador.modo` (`ark/metrologia/`, README_OPERADOR_DELTA.md). Por omision DISJUNTO.

| | DISJUNTO (ventanas independientes) | DESLIZANTE |
|---|---|---|
| que mira | la escala j al cerrar su barra 2^j; la temperatura al cerrar su ventana de T* barras | el nivel en CADA paso nativo en todas las escalas; la temperatura en cada barra de su escala |
| umbral h | el mas bajo para el mismo ARL0 (pocas miradas, independientes) | mas alto (miradas solapadas y correlacionadas): medido 10.37 contra 10.06 ... 6.78 contra 5.17 |
| saltos de nivel | iguales en saltos grandes (1, 2, 4 s0: 42.5 / 13 / 3.5 pasos) | iguales en grandes; 5-6 % mas lento en pequenos (0.25 y 0.5 s0), por el h mayor |
| primera alarma de temperatura | espera a cerrar la ventana, pero con h bajo: en los mundos medidos llega antes (x4: 32.5; x1.5: 257.5 pasos) | no espera la ventana pero paga el h mayor: llega despues (x4: 35; x1.5: 301.5). Ganaria si la escala que detecta tuviera T* m grande |
| confirmacion de temperatura | una segunda alarma con una ventana entera nueva | por acumulacion (h_c = ARL0/alfa): **confirma antes**, x4: 93 contra 137; x1.5: 493 contra 1081.5; x1.3: 1493.5 contra 1814.5 |
| falsas alarmas en calma | la referencia | familia: 0.85 de la referencia; confirmaciones falsas de temperatura: 2.3 veces en iid, la mitad en AR(1) 0.5 |
| coste | la referencia | 1.15-1.19 veces por paso; memoria de calibracion O(L 2^J) |
| superestadistica, regularizador, forma de DERIVA | sobre barras disjuntas | sobre barras disjuntas (las mismas: piden ventanas independientes) |
| cuando conviene | deteccion temprana del primer aviso; cuando las falsas confirmaciones cuestan | confirmar pronto un cambio de regimen de dispersion; cuando la decision espera la confirmacion |

## 3. Ventanas del sensor de memoria (H): solapadas, de W muestras cada S

| sensor | donde | ventajas | desventajas |
|---|---|---|---|
| Whittle local exacto (ELW) | `ark/metrologia/whittle_exacto.hpp`, `observador_memoria.hpp` | menos ruido; W y m derivados de la banda (suelo de ruido del nulo, techo interno y de localidad) | necesita historia para derivar la banda; banda no viable en series cortas |
| DFA-1 solapado (W = 50, S = 2) | `ark/fractal/sensor_dfa_causal.hpp`, `observador_fractal_multiescala.hpp` | barato, ventana corta (reacciona antes), sin banda que derivar | ruido alto sobre retornos (s0 0.15-0.2: no viable como gatillo sobre retornos); **pendiente calibrarlo sobre otros observables** (volumen, horquilla, OFI): no se descarta |

Las dos ventanas del sensor alimentan el mismo operador delta, en cualquiera de sus dos modos (seccion 2).

Medido juntos en el observador conjunto (`ObservadorDeltaFractal`, sinteticos): un cambio de H 0.5 -> 0.75 lo detectan las dos en 8/8, la
ELW en 205 barras y el DFA en ~380 (su H^ tiene sd 0.151 contra 0.077); el modo del operador sobre H^ no cambia el retardo (1.00 y 0.99).
La etiqueta conjunta usa la temperatura de la escala NATIVA del crudo, porque la de las escalas gruesas tambien cambia con la memoria.

## 4. Agregaciones multiescala (las dos quedan, decision del Comandante del 2026-11-05)

| | escalera 2^j (operador delta) | piramide 1m..D1 (fractal) |
|---|---|---|
| unidad | cuenta de barras (sirve en reloj intrinseco) | segundos de calendario |
| ventajas | escalas exactas y anidadas; camino por filas | escalas con significado humano (hora, dia) |
| desventajas | las escalas no son horas del reloj | sin camino por filas |

## 5. Calibracion del umbral h y presupuesto de la familia -- DEFECTO ABIERTO (2026-11-05); resuelto solo en coordinacion

**Estado (2026-11-05, tarde).** El operador de coordinacion tiene el modo `FAMILIA_EPISODIOS`: mide la tasa de episodios de la familia
sobre B realizaciones del nulo y biseca la rareza comun A. En calma da 8.2e-5 contra el limite 1.0e-4, donde antes daba 5.7e-4 (acta
`JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/calibracion_h/MEDIDO_FASE2.md`). **Los demas detectores** (ObservadorSaltos solo, memoria,
fractal, conjunto) **siguen con el reparto por canal y lo de abajo les aplica.**


Todos los detectores de la base fijan su umbral h por el ARL0 medido SOBRE EL NULO (`metrologia/canal_nivel.hpp`, `buscar_h_con`; el
fractal con su propia `h_por_arl`, del mismo tipo). El ARL(h) se lee con >= 30 alarmas en el nulo; si el ARL0 pedido no es legible, h se
EXTRAPOLA con una recta en log ARL(h) sobre [h_b/2, h_b] y la ficha lo marca (`h_extrapolado`). El presupuesto de una familia se reparte
por igual entre sus canales: `ARL0_canal = ARL0_familia x 2 x tot` (tot = suma de las J de todos sus observadores).

| | lo medido (sinteticos) |
|---|---|
| pocos canales, ruido gaussiano (operador delta solo; conjunto crudo + memoria, ~2 veces de extrapolacion) | conservador o en su ARL0 |
| muchos canales y volatilidad estocastica (operador de coordinacion, 90 componentes, ~140 veces de extrapolacion por canal) | **la familia emite 5.7 veces el limite**; el crudo solo, x14 |
| la misma familia con ARL0 menores (menos extrapolacion) | x0.6 cerca, x2 a ~10 veces, x7-14 a ~140 veces |
| prueba sobre el MISMO nulo de calibracion | el mismo exceso: el nulo no es el problema |

Hipotesis (no medida): con colas pesadas log ARL(h) es concava y la recta sobreestima el ARL al extrapolar, dando h demasiado bajo. Con
ruido gaussiano casi no pasa: log ARL ~ h + (1/2) log h. **Consecuencia practica, hasta resolverlo:** si una familia tiene sus canales
con `h_extrapolado` y el ARL0 por canal muy por encima de la historia, su tasa de falsas alarmas NO es la declarada. Opciones en
discusion y preguntas al auditor externo: `mensajes/CONTEXTO_AUDITOR_EXTERNO_CALIBRACION_H.md`; acta
`JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/coordinacion_sensor/MEDIDO_COORDINACION.md`.
