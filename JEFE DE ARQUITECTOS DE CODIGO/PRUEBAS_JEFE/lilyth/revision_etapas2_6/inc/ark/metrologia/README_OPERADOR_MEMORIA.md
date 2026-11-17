# Operador de memoria -- `ObservadorMemoria` (sensor ELW + operador delta y temperatura sobre H)

Codigo: `whittle_exacto.hpp` (estimador, `SensorWhittleExacto`, `derivar_banda`) y `observador_memoria.hpp`. Prerregistro y actas:
`JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/metrologia/EXPECTATIVAS_WHITTLE.md` (enmiendas 1-3) y `MEDIDO_WHITTLE.md`. Es telemetria.

## Que hace

observable crudo -> **Whittle local exacto** (Shimotsu y Phillips 2005, media desconocida de Shimotsu 2010; fuente verificada solo en
secundaria: Blevins 2025, arXiv:2511.15689) sobre ventanas de W muestras con m frecuencias -> H^(t) cada S muestras, con
H = d^ + 1/2 - orden (FLUJO 0, NIVEL 1: vale tanto para la serie de ruido como para su acumulado) -> `ObservadorSaltos` (ver
`README_OPERADOR_DELTA.md`) con FLUJO y agregacion MEDIA: el canal de nivel ve **saltos de H** y el de temperatura **regimenes de la
dispersion de H**.

Semiparametrico: solo supone la ley de potencias del espectro cerca de la frecuencia cero. No supone fGn ni nada sobre la memoria corta
(decision del Comandante: el bloque es agnostico al observable).

## Banda (W, m), toda derivada

1. Suelo de ruido (enmienda 3): parte de m = ceil(1/(4 s*^2)) y sube m hasta que la sd MEDIDA de d^ sobre el **nulo de ruido** quede
   <= s*.
2. Techo de la ventana interna: m <= ceil(W/w) - 1 (por encima de 2 pi / w una estadistica rodante de w muestras no tiene informacion).
3. Techo de localidad: el sesgo de d^(m; W), medido en sub-ventanas de W pareadas contra d^(m; W_ref) y sumado el sesgo de la propia
   referencia por la ley (m/n)^2, tiene que quedar <= s* (enmiendas 1 y 2). Se elige la menor W que lo cumple.
4. Si no hay W <= W_MAX: NO MEDIBLE con su motivo (no se inventa una W).

## Entradas obligatorias

| campo | que es |
|---|---|
| `banda.ruido_objetivo` (s*) | sd de H^ por ventana que exige el usuario |
| `nulo_ruido` | PERMUTACION o SIGNOS (de el se deriva el suelo; BLOQUES no mide ruido) |
| `delta.arl0_pasos_nativos` | ARL0 de la familia en pasos de H |
| `banda.ventana_interna` | w si el observable es una estadistica rodante (1 si no lo es) |
| `orden` | FLUJO o NIVEL del observable |
| `paso_sensor` (S) | resolucion temporal de H^ |

## Nulos que acepta (catalogo en `README_NULOS.md`)

| campo | valores | que hace |
|---|---|---|
| `nulo_ruido` | PERMUTACION, SIGNOS | historia sustituta (del CRUDO) sobre la que se mide el ruido del estimador y se fija m |
| `nulo_sensor` | NINGUNO, PERMUTACION, SIGNOS, BLOQUES (b = W) | con que H^ se calibra el operador delta: NINGUNO = H^ real (los regimenes lentos de H cuentan como reposo); un nulo = H^ de un sustituto del crudo (el reposo es "sin memoria") |
| `nulos_contraste` | mascara de los tres | H^ de cada nulo contra el real: media, sd, ACF1, tau_int y la diferencia con su se |
| `delta.nulo_nivel`, `delta.nulos_separacion` | los del operador delta | actuan sobre la serie H^ que recibe |

Todos los sustitutos se construyen sobre el observable CRUDO, no sobre H^: permutar H^ destruiria su propia autocorrelacion de solape
y daria un umbral demasiado bajo. Con NIVEL actuan sobre los incrementos y se reintegran.

Guia por variable: retornos y precios, `nulo_ruido` = SIGNOS (conserva la volatilidad, que tambien ensucia al estimador); volumen,
nu, horquilla, profundidad (en log), PERMUTACION; estadisticas rodantes (Shannon), el nulo hay que construirlo aguas arriba
(permutar la entrada y recalcular la estadistica) y declarar `ventana_interna`.

## Limites medidos (`MEDIDO_WHITTLE.md`)

- El sesgo de la W elegida se tolera hasta s* (MSE <= 2 varianzas): el nivel de H^ lleva ese sesgo y no se cita sin restarlo.
- La memoria corta que aparece DESPUES de calibrar se lee como salto de H (control de la trampa: 8 de 8).
- Una entropia de Shannon de 32 muestras no es medible a s* 0.08 con W <= 2048.
- Conexion (2026-11-04): vive en la base auditada y compila solo contra ella (`ark_cxx/base_auditada/README.md`, reparto de
  responsabilidades). La fase de la purga sigue saliendo del llamador, no de `ark/dataset/calendario.hpp` (pendiente).

## Entrada por filas en reloj intrinseco (2026-11-04, orden del Comandante)

Dos caminos con el mismo estimador:

| camino | llamadas | quien valida y remuestrea |
|---|---|---|
| serie ya validada (el de siempre) | `calibrar(cfg, x, ok, fase, dia)` y `paso(x, fase, dia)` | el llamador |
| filas crudas | `calibrar_filas(cfg, fuente, filas, fase, dia)` y `paso_fila(EntradaPasoCxx, fase, dia)` | `ark/ventanas/motor_reloj_intrinseco.hpp` |

`FuenteFilas` es la fuente comun de la base (`ark/ventanas/flujo_barras.hpp`, la misma del fractal y del ObservadorSaltos por filas):
`base` (la configuracion del motor de ventanas: reloj PARED / NU_TRADES / NU_VOLUMEN / NU_NOCIONAL / NU_INFO_SHANNON, tick,
presupuesto de imputacion, sesiones), `delta_nu` (cuanto nu forma una barra; NaN = la media de nu por fila con sesion y sin imputar
de la historia de calibracion), `agregacion` (ULTIMO para niveles, SUMA para flujos, MEDIA_NU) y `fase` (la fase de la purga desde el
calendario del activo: ver abajo).

- Las fronteras de barra son los multiplos de `delta_nu` del reloj acumulado, CON ARRASTRE del exceso (orden del Comandante,
  2026-11-04). Una fila no se reparte: si cruza k > 1 fronteras cierra una barra y las k - 1 restantes se cuentan como saltadas
  (`fronteras_barra`, `fronteras_saltadas`). Las filas de sesion cerrada congelan el reloj.
- Fase de la purga: con `fuente.fase.activa` (por ejemplo `ark::dataset::fase_calendario_de("BTCUSDT", 24)`), la fase y el dia de cada
  salida salen del cierre de la ultima barra de la ventana en hora LOCAL del mercado (`ark/dataset/fase_calendario.hpp`: tipo de dia
  laborable / no laborable con los feriados del calendario, bin horario local) y `delta.bins_fase` se fija a 2 bins_dia. Sin activar,
  los entrega el llamador por fila.
- **W del motor = la W que deriva el operador** (`banda().W`, Whittle local con su suelo y sus techos) y S = `paso_sensor`: no hay W de
  compilacion en el camino; W_MAX es solo el tope de memoria.
- Calibracion: la historia en filas se convierte en la serie de barras con el MISMO motor y la banda se deriva sobre esas barras.
- Una ventana con alguna fila imputada NO se estima (el sensor no recibe la ventana; la salida queda invalida y el operador delta la
  trata como hueco). Medido en F2 (0.2 % de filas imputadas, volumen 1..9 por fila, `delta_nu` = nu medio por fila): sin arrastre
  W = 170 barras eran ~255 filas (1.5 filas por barra) y quedaba sin estimar el 45 % de las salidas; CON arrastre, W = 168 barras y el
  36 % sin estimar, que corresponde a ~1.3 filas por barra ((1 - 0.002)^223 = 0.64). La barra sigue siendo algo mayor que una fila
  porque una fila con mas nu que `delta_nu` cruza varias fronteras y cierra una sola barra.
  Los huecos de calendario (sesion cerrada) NO invalidan: congelan el reloj.
- `abortar_por_presupuesto` no esta soportado en el motor intrinseco (`ErrorMotorIntrinseco::ABORTAR_NO_SOPORTADO`).
- El lazo de Shannon queda acotado por W_MAX.
- Controles (`PRUEBAS_JEFE/metrologia/prueba_whittle.cpp`, modo `filas`; `PRUEBAS_JEFE/ventanas_intrinseco/`;
  `PRUEBAS_JEFE/conexion_filas/`): F1 en PARED con delta_nu 1 el camino por filas es el de siempre bit a bit (banda y 14 919 salidas);
  F2 en NU_VOLUMEN con imputadas y cerradas el ciclo caliente reproduce la serie de barras (73 608 salidas, 0 distintas); F3 cero
  reservas en `paso_fila`; motor T1-T8 y T2b (arrastre) 16/16; G9 y G9b la fase de calendario en el ciclo caliente = la calibracion a
  mano con 48 fases.

## Modo del operador delta sobre H (2026-11-05)

`ConfiguracionMemoria.delta.modo` = DISJUNTO o DESLIZANTE: el operador delta sobre H^ elige su modo como sobre cualquier observable
(ventajas y desventajas en `README_OPERADOR_DELTA.md` y `ark/ventanas/README_VENTANAS_Y_DETECTORES.md`). La ventana del sensor (W, S) no
cambia con el modo: es la del ELW.

## Lente de memoria del observador conjunto (2026-11-05)

`ObservadorMemoria` expone la interfaz comun de las lentes (`Configuracion`, `ventana()`, `paso_sensor()`, `reiniciar_flujo()`,
`recalibrar_arl0()`), igual que `ObservadorMemoriaDFA`, para usarse en `ObservadorDeltaFractal` (`observador_delta_fractal.hpp`): un
observable, la lente del crudo (saltos de nivel y temperatura) y esta (cambios de memoria), con un solo ARL0 repartido por canales y la
etiqueta MEMORIA_PURA / MEMORIA_CON_VOLATILIDAD. Ver `README_OPERADOR_DELTA.md` (seccion del conjunto).

## Defecto abierto: el h del operador sobre H^ (2026-11-05)

El operador delta sobre H^ es el mismo `ObservadorSaltos` y hereda su defecto: si su ARL0 (en pasos de H^) queda muy por encima de lo que
la serie nula de H^ puede leer con 30 alarmas, el h se extrapola y la tasa real de falsas alarmas puede ser varias veces la declarada
(medido x2 a ~10 veces de extrapolacion y x7-14 a ~140 con colas pesadas). H^ estimado por ventanas hereda las colas del observable.
Ver `README_OPERADOR_DELTA.md` (seccion final) y `mensajes/CONTEXTO_AUDITOR_EXTERNO_CALIBRACION_H.md`.
