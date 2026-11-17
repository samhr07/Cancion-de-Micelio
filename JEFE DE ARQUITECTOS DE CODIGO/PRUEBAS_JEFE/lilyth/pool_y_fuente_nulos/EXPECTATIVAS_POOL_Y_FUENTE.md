# Pool derivado y fuente de los nulos de Lilyth -- expectativas selladas (ultima fase del operador)

Code (Opus 5.5), 2026-11-16. Ordenes del Comandante: "prepara ambas" (pool y fuente), nivel de decision = Holm de 4 de Lilyth, prueba
bilateral; ramificar la calibracion (dias sin viento para el operador delta, dias de viento para Lilyth); D1, D2 y D3 aprobados ("es
bueno tener la comparacion": la recuperacion con su propia marea temporal, un offset). Se sella ANTES de escribir codigo. Base previa en
`base_previa/` (verificada: identica, 70 archivos). Borrador anterior en `borrador_previo.md`.

Correccion de una afirmacion mia en la conversacion: Lilyth NO toma los segmentos "de calma del operador delta". Toma los tramos validos de
la historia de calibracion `x_cal` que le pasa el llamador (`calibrar`, paso 2); hoy esa historia es, en la practica, la misma calma del
delta. Por eso la fuente se ramifica con una mascara, sin tocar al operador delta.

## Diseno

A. Pool derivado (regla 6: entrada nueva, la de hoy por omision)
- A1. `ConfiguracionLilyth::pool_nulo` = `PoolNuloLilyth::{M_MIN = 0, DERIVADO = 1}`, por omision M_MIN (bit a bit).
- A2. DERIVADO: al empezar `evaluar_nulos_representacion`, K = numero de eventos candidatos (RESTITUIDA o CENSURADA_HORIZONTE) de la rama;
  objetivo por evento n* = ceil((2 * 4 / alfa)^(1/K)) - 1 (alfa 0.05: K=1 -> 159, K=2 -> 12, K=3 -> 5, K=4 -> 3); umbral del evento =
  max(m_min, n*). La via (a) se usa si la excedencia tiene >= umbral; si no, la via (b) con objetivo = umbral (su tope de intentos sigue
  siendo el de hoy). El evento queda evaluado si llega a m_min (como hoy); si llega a m_min pero no a n*, sigue evaluado y lo dice
  `resoluble` (la condicion exacta prod (n_k + 1) >= 160 ya la calcula la rama INTERCAMBIABLE). K se fija con los candidatos, no se
  recalcula tras evaluar (sin bucle). DERIVADO solo AGREGA nulos: nunca evalua un evento que hoy no se evalua.
- A3. Ningun numero nuevo: 2 (bilateral), 4 (Holm de Lilyth) y alfa.

B. Fuente de los nulos (regla 6)
- B1. `calibrar(...)` gana un ultimo argumento opcional `std::span<const std::uint8_t> fuente_nulo = {}`. Vacio: hoy, bit a bit. Con
  mascara: los segmentos del nulo son los tramos validos Y dentro de la mascara (un 0 corta el segmento); la purga, el largo de cada
  realizacion (n_cal) y B no cambian. La ficha registra `fuente_externa` y `n_fuente` (muestras usadas).
- B2. Las tres lecturas son mascaras del llamador sobre la misma historia de calibracion con viento:
  CALMA = fuera de las ventanas de las tormentas del catalogo; VIENTO = dentro de ellas; COMPLETA = todo.
- B3. D2 (deja-uno-fuera): los eventos contrastados viven en el flujo de operacion, disjunto de la historia de calibracion; el arnes lo
  verifica. El README lo deja escrito como responsabilidad del llamador.
- B4. Las ventanas salen de la nota del bufer (dos pasadas): pasada 1 = Lilyth calibrado en calma corre sobre la historia de calibracion
  como flujo y mide tau_rest de las tormentas del catalogo; bufer posterior = cuantil 1 - alfa de tau_rest consolidado (censurado ->
  T_max efectivo); bufer previo = mayor rezago l <= T_max con tasa de apertura previa a la tormenta sobre 1/ARL0 (binomial unilateral,
  Holm sobre los rezagos), 0 si ninguno.

## Mundos sinteticos (verdad por construccion)

Proceso comun: flujo y_t con AR(1) de umbral (phi = 0.5 si |y_{t-1}| < 4, phi = 0.92 si no; innovacion N(0,1)): en calma decae deprisa,
en las excursiones grandes despacio. Es la "memoria normal" con agrupamiento: ninguna excursion grande es especial. Tormenta = innovacion
de amplitud 12 con signo alterno en t0. Mecanismo lento (solo donde se dice): componente aditiva z = 4 * 0.97^tau sumada al flujo
observado x = y + z (no entra en el umbral).
- M0: ninguna tormenta lleva z. Verdad: ninguna recomposicion es lenta mas alla de la memoria normal.
- M1: todas las tormentas (calibracion y operacion) llevan z.
- M2: una de cada cuatro tormentas lleva z, en calibracion y en operacion.
Delta calibrado en una serie de calma independiente del mismo proceso (3000 pasos, sin tormentas). Historia de calibracion de Lilyth:
6000 pasos con una tormenta cada 300 (catalogo conocido). Flujo de operacion: 3200 pasos, tormenta cada 200 (16). Replicas: 40 por mundo.
Lectura por serie y referencia: "rechaza lenta" = p_holm de esa referencia < alfa y S > 0 (como E7b), con la rama INTERCAMBIABLE.

## Expectativas

Condiciones (deciden):
- EC5 (identidad): con mascara vacia y M_MIN, las baterias de Lilyth (Etapas 2-6) y el volcado del acople de E7b en 3 mundos (IID,
  RAFAGAS, CAMBIO) son identicos bit a bit a la base previa (mismos instrumentos que XB.0; la excepcion esporadica conocida de XB.0 se
  reporta si reaparece y no se re-explica).
- EC6a (pool, mecanica): con DERIVADO, en series de un solo evento, cada evento evaluado tiene >= m_min nulos; los que llegan a 159 dan
  resoluble verdadero y p_min <= alfa/4, los que no dan resoluble falso; `resoluble` coincide con prod (n_k + 1) >= 160 en todas.
- EC6b (pool, tamano): mundo lineal (AR(1) phi 0.82, calibracion en calma, tormenta = innovacion grande, K = 1, 400 replicas):
  con DERIVADO, rechazo de Holm de Lilyth en la referencia Onsager con Wilson inferior <= alfa.
- EC8 (fuente, mecanica): con mascara, `n_fuente` = muestras validas dentro de la mascara; mascara toda en 1 = mascara vacia bit a bit;
  los eventos de operacion no se solapan con la historia de calibracion.

Expectativas (se reportan se cumplan o no):
- EC6c (potencia): mismo mundo lineal pero la tormenta lleva z: DERIVADO rechaza (Onsager, por lenta) en >= 50 %; M_MIN en 0 % (con K = 1
  su p minimo es 2/9).
- EC1 (M0, Onsager): VIENTO tamano cercano a nominal (Wilson inferior <= alfa).
- EC2 (M0, Onsager): CALMA rechaza "lenta" por encima de alfa (Wilson inferior > alfa): la calma no tiene excursiones grandes que
  recomponer, su modelo lineal decae como phi 0.5.
- EC1b (M0, Onsager): COMPLETA rechaza menos que CALMA y mas o igual que VIENTO (un modelo lineal no reproduce la memoria de umbral).
- EC3 (M1, Onsager): CALMA y COMPLETA rechazan por lenta en >= 80 %; VIENTO no (Wilson inferior <= alfa): ciega a lo universal.
- EC4 (M2, por evento, Onsager): con VIENTO la media de u_pit de las tormentas sin z esta en [0.4, 0.6] y la de las tormentas con z es
  mayor; con CALMA la media de u_pit de las tormentas sin z es > 0.6.
- EC9 (bufer): el bufer previo derivado es 0 en los tres mundos (las tormentas sinteticas no avisan); el posterior es mayor en M1 que en M0.
- EC7 (mutantes): >= 8 de 10 cazados; los equivalentes se declaran. Lista: mascara ignorada; mascara invertida; mascara que no corta el
  segmento; n* sin el "- 1"; n* con alfa en vez de 2*4/alfa; via (a) con m_min en vez del umbral; via (b) con m_min en vez del umbral;
  K = eventos evaluados en vez de candidatos; DERIVADO por omision; piso m_min retirado (evento con < m_min evaluado).
