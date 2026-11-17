# Nulos del bloque ARK-METRO-01 (y del Bloque Fractal): catalogo comun

Decision del Comandante (2026-11-04): **los nulos son ENTRADAS del sistema y dependen de la variable a observar**. Ningun operador
fija un nulo por dentro; el usuario elige, como elige las familias de Holm. Este archivo es el catalogo; cada operador tiene su README
con que nulos acepta y para que (`README_OPERADOR_DELTA.md`, `README_OPERADOR_MEMORIA.md`, `../fractal/README_OPERADOR_FRACTAL.md`).

Codigo: `sustitutos.hpp` (`NuloSustituto`, `sustituto(segs, nulo, bloque, mu, rng)`). Cada nulo actua **dentro de cada tramo
contiguo** valido: nunca mezcla datos separados por un hueco. Con un observable NIVEL (precio, nivel acumulado) el sustituto se aplica
a los **incrementos** y se reintegra desde el primer valor del tramo (`observador_memoria.hpp`, `sustituto_crudo`).

| nulo | que conserva | que destruye | pregunta que contesta |
|---|---|---|---|
| `NULO_NINGUNO` (en el fractal: `HISTORICO`) | todo: es la historia tal cual | nada | "cambio respecto del reposo habitual de ESTA variable" (los regimenes habituales cuentan como reposo) |
| `NULO_PERMUTACION` | la marginal (cada valor) | todo orden temporal: memoria, agrupamiento de volatilidad, equilibrio local | "hay algo mas que la marginal" (memoria de cualquier tipo) |
| `NULO_SIGNOS` | `|y - mu|` en su orden (la volatilidad y su agrupamiento) | la direccion, y con ella toda autocorrelacion lineal | "hay memoria de la direccion mas alla de la volatilidad" |
| `NULO_BLOQUES` (b) | el contenido de cada bloque de b muestras | el orden entre bloques (lo lento, mas alla de b) | "algo persiste mas alla de b" (separacion de escalas) |

## Que nulo usar segun la variable

| variable | tipo | nulo de ruido / de memoria | por que |
|---|---|---|---|
| log retornos, retornos de apertura y cierre | FLUJO simetrico | **SIGNOS** | conserva la heterocedasticidad, que tambien ensucia a los estimadores de memoria; es el nulo real de EL ABISMO |
| precio, log precio | NIVEL | **SIGNOS** sobre los incrementos | igual que los retornos, reintegrado |
| flujo firmado (desequilibrio, signo agresor) | FLUJO simetrico | SIGNOS o PERMUTACION | SIGNOS si la pregunta es la direccion; PERMUTACION si es cualquier memoria |
| volumen, nocional, transacciones por minuto (nu), Q | FLUJO positivo | **PERMUTACION** (sobre el log si el operador recibe el log) | SIGNOS cambia de signo `y - mu` y puede dar valores imposibles (volumen negativo) |
| horquilla, profundidad, volatilidad realizada | intensivo positivo | **PERMUTACION** (sobre el log) | idem |
| razones (taker_buy_ratio, ticket) | intensivo acotado | PERMUTACION | SIGNOS sacaria valores fuera de [0, 1] |
| metrica de Shannon u otra estadistica rodante de w muestras | intensivo, rodante | PERMUTACION **de la serie de entrada, recalculando la estadistica** (aguas arriba) | permutar la estadistica ya construida destruye tambien la memoria de su propio solape y subestima el ruido; hoy los operadores reciben una sola serie, asi que este nulo hay que construirlo fuera y pasarlo como historia (limitacion declarada) |
| cualquier variable, pregunta "cambio respecto de lo habitual" | -- | NINGUNO | calibra con la propia historia |

BLOQUES no es nulo de ruido de un estimador de memoria (conserva la memoria dentro de la ventana): se usa para la separacion de
escalas del canal de temperatura (b = T* 2^k) y para preguntar si algo persiste mas alla de b.

**Nulo de CONTRASTE contra nulo de CALIBRACION (2026-11-11, decision del Comandante).**

- **La tabla de arriba es para contrastar** (estimadores de memoria, separacion de escalas): ahi destruir la memoria es la pregunta.
- **Como `nulo_nivel` del operador delta es otra cosa:** ese nulo CALIBRA el canal de nivel (mu0, s0, tau_int y h).
  - Con PERMUTACION o SIGNOS sobre una serie con memoria, tau_int sale el de la serie sin memoria (1.00 contra 15.6 en un AR(1) phi
    0.88) y el ARL0 no se cumple: 28-30 veces el presupuesto de alarmas, tambien sin choques.
  - El volumen, nu, la horquilla y la volatilidad realizada tienen memoria fuerte en este proyecto (Hurst del volumen 0.79-1.02). Para
    calibrar el operador delta sobre ellos, `nulo_nivel` = NINGUNO, no PERMUTACION.
  - Sigue siendo una entrada: se elige segun la serie. Detalle en `README_OPERADOR_DELTA.md`; acta en
    `JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/auditoria_delta_choques/`.

## Nulos del operador de coordinacion de pares (2026-11-05)

Viven en `../coordinacion/sensor_coordinacion.hpp`, no en `sustitutos.hpp`, porque actuan sobre DOS series:

| nulo | sobre que | conserva | destruye | para |
|---|---|---|---|---|
| SIGNOS de centro CERO | cada retorno base de X y de Y, por un +-1 independiente | \|r\| exacto, y con el la varianza realizada de todo bloque y su co-movimiento entre activos | la dependencia de signo (dentro y entre activos) | DIRECCION y DIRECCION_NORMALIZADA |
| ROTACION | Y desplazada n/2 contra X (valores y mascara) | la memoria de cada serie (volatilidad incluida) | el emparejamiento entre activos | MAGNITUD |

DIVERGE de `NULO_SIGNOS` de `sustitutos.hpp`, que centra en la media: con media no nula no conserva |r|.

## Preguntas abiertas sobre la validez del nulo (2026-11-05)

- La validez de un nulo para un observable se puede medir: la tasa de falsas alarmas sobre una prueba nueva tiene que ser la misma que
  sobre el propio nulo de calibracion. Solo se ha medido para retornos (y para las informaciones del operador de coordinacion, en
  sinteticos); no para volumen, nocional, Q ni las metricas propias. Medido el 2026-11-11 para `nulo_nivel` = PERMUTACION sobre un AR(1)
  phi 0.88: 1 062 falsas alarmas sobre una prueba nueva contra un presupuesto de 37.5. **No es valido con memoria.**
- Una familia que mezcle observables de tipos distintos (flujo simetrico, positivo, razon) necesita un nulo CONJUNTO valido para todos;
  ningun sustituto de esta tabla lo es para todos a la vez.
- Con el reloj de Shannon (barras cortadas por sorpresa) no esta decidido si el sustituto va sobre las filas (cambia el propio reloj) o
  sobre las barras ya cortadas (conserva un reloj definido con los datos reales).
Contexto: `mensajes/CONTEXTO_AUDITOR_EXTERNO_CALIBRACION_H.md`.
