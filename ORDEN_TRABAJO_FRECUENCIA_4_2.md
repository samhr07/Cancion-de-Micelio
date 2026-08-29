# Orden de Trabajo — La Frecuencia Óptima de Operación · v4.2

Sucede a `ORDEN_TRABAJO_HORIZONTE_4_1.md` y a su `ADENDA_C_REFORMULACION_S1_4_1.md`.

**El cambio de planteamiento.** Hasta ahora la pregunta era «¿existe señal a horizonte `H`?». Es la
pregunta equivocada, y su forma binaria ha producido cuatro veredictos que no decidían nada.

La comisión es un **impuesto a la frecuencia**, no al tamaño. Con `maker+maker` a 3.6 pb ida y
vuelta y ciclo de trabajo del 10 %, el coste anual solo en comisiones es el **946 % del capital a
2 minutos de posición, el 63 % a 30 minutos, el 8 % a 4 horas y el 1.3 % a un día.** No hay
capital que arregle una frecuencia mal elegida, y hay frecuencias donde el problema de coste
simplemente deja de existir.

**La pregunta correcta, y la que este documento contrasta:**

> ¿Existe un par (periodo de tenencia `H`, umbral de operación `θ`) donde la rentabilidad neta
> esperada por unidad de tiempo sea positiva y distinguible de cero?

Es un problema de optimización con óptimo **interior**: a `H` corta el coste fijo domina; a `H`
larga el movimiento crece pero se opera menos veces y la financiación se acumula. **El óptimo
existe aunque no sea rentable**, y localizarlo es el resultado, cruce o no cruce el cero.

`Micelio.py` no se toca. El conjunto de prueba de la v3.2 no se abre.

---

## §0. Registro de errores y deudas heredadas

| Origen | Estado | Dónde se salda |
|---|---|---|
| `c(u)` con tarifas **ASUMIDAS** VIP 0; criterio §8 de `PREREGISTRO_3_2.md` **no cumplido** desde la v3.1 | **ABIERTO.** Toda la curva de coste de este documento descansa en ese número | §1, compuerta |
| `C_respaldo` sin medir; sin él maker gana siempre porque no paga nada al fallar | **ABIERTO** | §3 |
| Selección adversa medida solo a 1 s, 5 s y 50 s; los horizontes de este documento son 10³–10⁵ s | **ABIERTO** | §2 |
| Tasa de llenado `p` medida solo a 60 s (41.3 %) | **ABIERTO** | §2 |
| Financiación del perpetuo nunca medida sobre dato propio; se usa un valor de referencia | **ABIERTO** | §2 |
| §1 de la v4.1: instrumento sin resolución más allá de 2 min | **Resuelto en la Adenda C**, pendiente de ejecutar | §1, compuerta |
| Micro-precio (§6 de la v4.1), `η̂` (§4), calibración contemporánea contra Cont–Kukanov–Stoikov, curva en U, nulo espectral de `ν` | **PENDIENTES**, ninguna bloquea | §7 |

⚠ **Ningún resultado de este documento es citable mientras `c(u)` siga siendo un número asumido.**

---

## §1. Compuertas de entrada — nada de lo demás se ejecuta antes

### 1.1 Tarifas reales

Consultar `/fapi/v1/commissionRate` con endpoint firmado **al arrancar y sin cachear**, junto con
el nivel VIP efectivo de la cuenta y el estado del descuento BNB. Reportar:

```
maker real, taker real, nivel VIP, descuento BNB activo (si/no)
c(u) maker+maker en pb  y  el factor del paso 3 recalculado
```

Referencia externa **a verificar, no a creer** (fuente secundaria, julio 2026): futuros USDⓈ-M
VIP 0 maker 0.0200 %, VIP 1 0.0180 %, VIP 2 0.0160 %, VIP 3 0.0120 %; descuento BNB en futuros
10 %; ruta «Holder VIP» por activos en 100 000 / 200 000 / 3 000 000 USD para los niveles 1/2/3.
**Si la tarifa real difiere, manda la real y esta orden se recalcula entera.**

### 1.2 El instrumento de la Adenda C

Ejecutar `§C.4` puntos 1, 2 y 3 **antes** de cualquier barrido: suelo de ruido con ≥ 200
barajados por horizonte, calibración de la identidad contra la medición directa a 60 s y 120 s, y
`R_cum(H)` y `σ²_r(H)` **por fragmento**.

⚠ **Compuerta dura:** un horizonte solo entra al barrido del §4 si
`R²_req(H) ≥ 3 · q95(|R²_nulo|)` en ese horizonte. Los que no la pasan **no se barren**: se
reportan como no resolubles con la cobertura que harían falta. Barrer un horizonte ciego es
fabricar el máximo del §5.

---

## §2. La curva de coste, medida en función de `H`

Hasta ahora el lastre se ha tratado como una constante. **No lo es**, y sus tres componentes se
mueven en direcciones distintas con el horizonte. Esto es medición, no modelo.

```
L(H)  =  c(u)  +  A(H)  +  F(H)  +  (1 − p(H))·C_respaldo
```

| término | qué es | cómo se mide | dirección esperada |
|---|---|---|---|
| `c(u)` | comisión ida y vuelta | §1.1, endpoint firmado | **constante en pb** |
| `A(H)` | selección adversa | markout del punto medio a `H` tras el llenado, sobre `captura_v32` | crece y **debería saturar** |
| `F(H)` | financiación del perpetuo | tasa real histórica de `captura_estacional`; **es con signo** | crece ~lineal con `H` |
| `p(H)` | tasa de llenado maker | misma definición que la v4.0 (con «NIVEL BARRIDO» contado como llenado) | crece con `H` |
| `C_respaldo` | coste de no llenar | §3 | constante |

⚠ **`F(H)` puede ser negativo y eso no es un error.** La financiación se cobra o se paga según el
lado y el signo de la tasa. Un sistema que opera condicionado al signo de la financiación tiene un
lastre menor que uno que no lo hace. **Reportar `F(H)` como distribución con su signo, no como
coste medio en valor absoluto** — es la primera vez que este proyecto mira la financiación y
tratarla como un coste sin signo destruiría información.

⚠ **La financiación es discreta, no continua.** Se liquida cada 8 h en instantes fijos. A `H` menor
que 8 h, una posición paga financiación **solo si cruza uno de esos instantes**, con probabilidad
≈ `H/28800`. Modelarlo como coste continuo `F·H/28800` es una aproximación al valor esperado y
**subestima la varianza**. Reportar las dos cosas: el valor esperado y la fracción de operaciones
que efectivamente cruzan una liquidación.

### 2.1 Salida forzada

`c(u)` tiene dos valores medidos en este proyecto: **26.03 USD/BTC** con salida maker y **37.65**
con salida forzada (v4.0). La diferencia son 2.7 pb, comparable a toda la señal medida. **Ambos
van en toda tabla, como dos columnas**, y `p(H)` es lo que determina la mezcla real.

---

## §3. `C_respaldo` — se salda aquí, y su definición se declara antes de medir

Es el término que falta desde que se identificó y sin él la ecuación de coste maker está sesgada a
favor de maker por construcción.

**Definición operativa, declarada antes de mirar dato:** para cada orden maker colocada y **no
llenada** dentro de la ventana de espera, `C_respaldo` es el desplazamiento adverso del punto medio
entre la colocación y el momento de cancelación, en pb, con signo respecto de la dirección
pretendida. Es lo que costó no haber cruzado al principio.

Medir sobre `captura_v32`, reportar la distribución completa (no solo la media) y **por horizonte
de espera** `H` ∈ la misma rejilla del §4, porque `p` depende de `H` y el término entra como
`(1−p(H))·C_respaldo`.

⚠ **Nota de prioridad, corregida respecto de la v4.1.** Allí se escribió que `C_respaldo` pierde
peso al crecer `H` porque `p` sube. Es cierto para el factor `(1−p)` y **falso para el otro
factor**: cuanto más esperas, más se puede haber escapado el precio, así que `C_respaldo` **crece
con `H`**. El producto no es monótono y **hay que medirlo, no razonarlo**. Ese error de
razonamiento se registra aquí.

---

## §4. El experimento — la superficie `R(H, θ)`

### 4.1 El objeto

Ganancia neta esperada por ida y vuelta, en pb, operando solo cuando la señal supera el umbral:

```
g(H, θ)  =  E[ movimiento capturado | |señal| > θ ]  −  L(H)
```

Rentabilidad anual sobre capital, apalancamiento 1:

```
R(H, θ)  =  d(H, θ) · (T_año / H) · g(H, θ) / 10⁴
```

con `d(H, θ)` la **fracción de tiempo medida** en que la señal supera el umbral — no un 10 %
supuesto.

⚠ **`E[movimiento capturado]` se mide empíricamente, no se deriva de `R²`.** El factor
`κ = 1.755` que se ha venido usando para pasar de `R²` a movimiento del decil superior supone
normalidad, y la distribución de retornos de este mercado tiene curtosis 1179.7. **Reportar el
`κ` empírico medido** junto al gaussiano, en toda fila. Si difieren en más del 20 %, todas las
curvas requeridas publicadas hasta hoy están mal por ese factor al cuadrado.

### 4.2 La rejilla, declarada y congelada antes de mirar

```
H  ∈  {15 min, 30 min, 1 h, 2 h, 4 h, 8 h, 24 h, 72 h}          8 valores
θ  ∈  {q50, q70, q80, q90, q95, q99} de |señal| en ENTRENAMIENTO  6 valores
```

48 celdas. Predictor: el mismo del §1.3 de la v4.1 —flujo firmado acumulado, sin `G(τ)`, sin
`τ₀`, sin `β`— y **además** una segunda columna con la financiación como predictor exógeno, que es
la única variable de este documento que no se disipa en segundos. Las dos columnas siempre, como
punto medio y precio de transacción.

⚠ **`H = 72 h` está fuera del alcance verificado (711 s) por dos órdenes de magnitud** y depende de
que `captura_estacional` acumule. Va en la rejilla porque el óptimo previsto cae cerca del borde
superior y truncar la rejilla antes del óptimo lo desplazaría artificialmente hacia dentro. Marcado
como extrapolación en la propia tabla.

### 4.3 ⚠ El modo de fallo que este documento introduce, y no existía antes

**Un barrido de 48 celdas produce un máximo aunque no haya nada.** Es la primera vez que este
proyecto elige parámetros sobre una rejilla, y es exactamente donde muere la disciplina de nulos
que ha sostenido todo lo anterior. Cuatro guardas, las cuatro obligatorias:

1. **La rejilla se congela con hash antes de tocar dato.** Ampliarla después es una enmienda con
   su motivo y su constancia de qué se había visto.
2. **Superficie nula.** Repetir el barrido completo sobre ≥ 200 realizaciones con los signos de
   flujo barajados dentro de cada ventana, y reportar la **distribución del máximo de la
   superficie** bajo el nulo. La cifra que decide no es `R(H*, θ*)`: es
   **`P(max R_nulo ≥ R(H*, θ*))`**. Sin ese número el máximo no se reporta.
3. **La celda ganadora se valida fuera de la muestra donde se eligió.** Elegir en entrenamiento,
   validar en validación. Si el proyecto quiere el conjunto de prueba, hace falta preregistro
   nuevo y autorización explícita de Samuel.
4. **Se reporta la superficie entera**, no el máximo. Un óptimo aislado en una superficie ruidosa
   es sobreajuste; un óptimo con vecinos que también son buenos es una propiedad. **La suavidad de
   la superficie es evidencia y hay que mirarla.**

### 4.4 Expectativa declarada ANTES de ejecutar

Con la aritmética disponible —lastre fijo 4.49 pb, financiación 0.968 pb por 8 h, `σ₁` medida,
ciclo de trabajo 10 %— y suponiendo `R²` **constante** con el horizonte (el caso más favorable
posible, porque la capacidad predictiva decae):

| `H` | `R²` = 0.4 % | 1.0 % | 2.0 % | 4.0 % |
|---|---|---|---|---|
| 15 min | −113 % | −83 % | −53 % | −9 % |
| 30 min | −48 % | −27 % | −5 % | +26 % |
| 1 h | −18 % | −3 % | +12 % | +35 % |
| 2 h | −5 % | +6 % | +16 % | +32 % |
| **4 h** | **+0.3 %** | **+8 %** | +15 % | +27 % |
| 8 h | +2 % | +7 % | +13 % | +21 % |
| 1 día | +2 % | +5 % | +8 % | +13 % |
| 3 días | +1 % | +3 % | +5 % | +7 % |

*(rentabilidad anual sobre capital, `H_p = 0.500`, apalancamiento 1, sin `C_respaldo`)*

**Tres predicciones que se declaran ahora y se contrastan después:**

- **P1 — El óptimo cae entre 2 y 8 horas.** Es el punto donde el coste fijo ya se diluyó y la
  financiación todavía no domina.
- **P2 — Todo lo que esté por debajo de 30 minutos sale negativo**, incluso con `R²` generoso.
- **P3 — La rentabilidad en el óptimo es de un dígito porcentual anual** salvo que el `R²` medido
  supere el 2 %, que es más de lo que este proyecto ha medido nunca fuera de muestra.

⚠ **Si el resultado no cumple P1 y P2, sospechar del cálculo antes que celebrar.** Esas dos
predicciones salen de aritmética de costes que no depende de ninguna hipótesis de mercado; si el
experimento las contradice, lo más probable es que haya un error en la contabilidad del lastre.

---

## §5. Regla de lectura, escrita antes

| condición | lectura |
|---|---|
| `P(max R_nulo ≥ R(H*, θ*)) > 0.05` | **no hay óptimo.** El máximo es lo que da una rejilla de 48 celdas por azar. Se reporta el `p` y se cierra |
| `p ≤ 0.05`, la celda ganadora **valida fuera de muestra**, y sus vecinas también son positivas | **hay banda viable.** Se preregistra el contraste económico completo, y solo entonces se habla de capital |
| `p ≤ 0.05` pero la celda no valida, o está aislada | **sobreajuste.** Se reporta y no se decide |
| El óptimo cae en el borde de la rejilla (`72 h`) | **rejilla truncada.** No se lee el máximo; se amplía la rejilla con enmienda declarada |
| La superficie es negativa en toda celda resoluble | **resultado.** «No hay frecuencia rentable con este predictor y este coste», y se reporta como tal |

⚠ **La última fila es un desenlace tan válido como las otras y probablemente el más informativo.**
No hay que evitarlo.

---

## §6. Lo que este documento NO hace

- **No toca `Micelio.py`.** Sigue sin cambios desde la v2.2.
- **No abre el conjunto de prueba de la v3.2.** Ya se abrió dos veces.
- **No implementa el coste lineal del NMPC** (v3.4). Va después de que exista `c(u, estado)` y de
  que alguna celda pase.
- **No ejecuta `PLAN_CAPITAL_5_0.md`.** Sigue condicional.
- **No optimiza el predictor.** El predictor tonto está deliberadamente fijado. Refinarlo antes de
  saber si alguna frecuencia funciona es optimizar dentro de una región que puede no existir.
- **No añade apalancamiento al cálculo.** Multiplica ganancia y pérdida por igual y no cambia el
  signo de ninguna celda. Es una decisión de riesgo posterior, no de investigación.
- **No resucita `ω_m` ni decide `ω_G`.**

---

## §7. Deudas que se saldan aquí, fuera de la ruta crítica

Ninguna bloquea el §4. Se ejecutan después, en este orden.

**7.1 Calibración contemporánea contra Cont–Kukanov–Stoikov (2014).** Regresión
`∆mid_k = α + β·OFI_k + ε_k` en rejilla de 10 s, `β` por submuestras de media hora, errores de
White. Referencias verificadas contra el artículo primario (arXiv:1011.6402v3): `R²` medio 65 %;
**35–60 % excluyendo del OFI los eventos que cambian el precio**, que es su propio control de
tautología y el número honesto; desequilibrio de transacciones solo 32 %; ambos juntos 67 %;
término cuadrático insignificante; `λ̂ ≈ 0.98` en `β = c/AD^λ`.
**Interpretación declarada antes:** un `R²` muy por debajo de ese rango **no es un hallazgo sobre
BTCUSDT**, es la medida de la degradación de nuestro OFI-L1 aproximado, cuyo residuo de
reconciliación es del 93.5 %. Si `λ` sale lejos de 1, entonces sí hay argumento para ingerir
`@depth`; hasta ahora no lo hay.

**7.2 Curva en U: ¿volatilidad o actividad?** Contrastar `σ_t` contra `ν_t` como predictores de
`σ_{t+1}` fuera de muestra en las cuatro capturas. Si `ν` domina, la pregunta abierta del §4.3 del
traspaso queda resuelta a favor de **actividad** y se retira de la lista, citando Clark (1973) y
Tauchen & Pitts (1983) **como marco, no como hallazgo**.

**7.3 Micro-precio (§6 de la v4.1).** Solo el punto 1: fracción de ceros contra el 96.67 % del
punto medio. Si no baja sustancialmente, se abandona ahí.

**7.4 `η̂` con barrido de colapso de ráfagas** en `{0, 10, 50, 200} ms`. Deuda de reporte: el
veredicto `η̂ = 0.62` se emitió sobre un contador que puede estar contaminado por barridos de
varios niveles contados como continuaciones de 1 tick.

**7.5 Nulo espectral sobre `ν_t`.** El nulo multitaper de la v2.2 se midió sobre dato contaminado
y sobre la serie de **precio**; la hipótesis de ciclo endógeno vive en la serie de **actividad** y
nunca se contrastó. Nulo de ruido rojo simulado, nunca tabla asintótica, con la guarda de banda de
la v2.1 §2. **Expectativa declarada: no habrá pico.**

---

## §8. Modos de fallo

| # | fallo | firma | detección |
|---|---|---|---|
| 1 | Máximo de una rejilla de 48 celdas leído como hallazgo | óptimo aislado, no valida fuera de muestra | §4.3, superficie nula y `P(max)` |
| 2 | Barrer horizontes donde el instrumento no resuelve | celdas con valores grandes y erráticos | §1.2, compuerta dura de resolución |
| 3 | `c(u)` asumido | toda la curva de coste desplazada | §1.1, endpoint firmado |
| 4 | Financiación tratada como coste sin signo | se pierde una fuente de ventaja | §2, distribución con signo |
| 5 | `κ = 1.755` gaussiano sobre datos de curtosis 1179.7 | error al cuadrado en toda curva requerida | §4.1, `κ` empírico al lado |
| 6 | Óptimo en el borde de la rejilla | truncamiento leído como óptimo | §5, fila 4 |
| 7 | `d(H, θ)` supuesto en 10 % | rentabilidad inflada o deflactada por un factor libre | §4.1, `d` medido |
| 8 | `C_respaldo` omitido «porque pesa poco a `H` grande» | maker sale siempre ganador | §3, incluida la corrección del razonamiento previo |
| 9 | Ventanas solapadas en el conteo de vueltas | vueltas por año infladas | test de no solapamiento |

---

## §9. Criterios de aceptación

- [ ] Tarifas reales del endpoint firmado; `c(u)` deja de ser asumido; factor del paso 3 recalculado.
- [ ] Compuerta de resolución de la Adenda C aplicada; horizontes no resolubles excluidos y listados.
- [ ] `L(H)` con sus cinco términos medidos, dos columnas de `c(u)` (maker y salida forzada).
- [ ] `F(H)` con signo y con la fracción de operaciones que cruzan liquidación de financiación.
- [ ] `C_respaldo` con definición declarada antes, distribución completa, y por `H`.
- [ ] Rejilla congelada con hash antes de tocar dato.
- [ ] Superficie `R(H, θ)` completa reportada, no solo el máximo.
- [ ] `P(max R_nulo ≥ R(H*, θ*))` con ≥ 200 realizaciones.
- [ ] Celda ganadora validada fuera de la muestra de selección.
- [ ] `κ` empírico junto al gaussiano en toda fila.
- [ ] `d(H, θ)` medido, nunca supuesto.
- [ ] Contraste de P1, P2 y P3 declarado y reportado, se cumplan o no.
- [ ] `Micelio.py` sin cambios (`git diff --stat` vacío en el reporte). Suites en verde.

---

## §10. Prioridad

1. **§1.1** — tarifas reales. Minutos, y sin ellas nada es citable.
2. **§1.2** — instrumento de la Adenda C. Horas.
3. **§2 y §3** — la curva de coste y `C_respaldo`. **Es la mitad del experimento y es la mitad que
   este proyecto nunca ha medido bien.**
4. **§4** — la superficie. El experimento que decide.
5. **§7** — las deudas, en su orden.

---

## §11. Por qué esto es una decisión y no otra iteración

Si el §4 devuelve una banda viable, la siguiente pregunta es de capital y ya está contestada: la
escalera de descuentos de Binance no empieza hasta los 100 000 USD, así que **por debajo de esa
cifra el capital no es una palanca** y la conversación es sobre si el rendimiento justifica las
horas.

Si el §4 devuelve la superficie entera negativa, es una **refutación limpia y publicable**: se
habrá medido que en este mercado, con este predictor y con estos costes, no existe frecuencia de
operación rentable para un participante que paga tarifa minorista. Eso cierra el proyecto con un
resultado, no con un abandono.

**Las dos salidas son aceptables. La que no lo es es seguir iterando sin haber medido la
superficie.**
