# Preregistro — ORDEN_TRABAJO_MIGRACION_3_2

**Fecha de escritura: 2026-08-08.** Este documento se commitea **antes** de ejecutar cualquier
medición del §3 en adelante de la v3.2. Su hash se referencia en el reporte final.

> §9 de la orden: "Commitear `PREREGISTRO_3_2.md` con fecha antes de ejecutar nada; referenciar
> el hash en el reporte. Es el tercer documento de este proyecto que lo exige y las dos veces
> anteriores salvó un veredicto."

**Estado de la captura al escribir esto:** `telemetria/captura_v32` lanzada el 2026-08-08 a las
20:13 con `--horas=48`, primera captura que persiste `b`, `B`, `a`, `A` y `u` (commit
`16bcec2`). **No se ha mirado ni un dato de ella.** Las capturas anteriores no sirven para esta
decisión: ninguna trae las cantidades del libro.

---

## 0. Qué se decide, y qué NO

Se decide **cuál** de estos cuatro modelos anidados sobrevive fuera de muestra:

```
Δp_t = Σ_{k≥0} h(k)·ε_{t−k}·f(v_{t−k}) + η_t     G(τ) = G₀(1 + τ/τ₀)^(−β)     f(v) = v^δ
```

| id | restricción | parámetros libres |
|---|---|---|
| **M0** — reloj de ticks sin capa predictiva | `G₀ = 0` | 0 |
| **M1** — MkII / OFI | `β = 0`, `δ = 0.5` | 1 |
| **M1′** — MkII con `c_t` adaptativo | `β = 0`, `δ = 0.5`, `G₀ → Y·σ_t/√V_best` | 1 + def. |
| **M2** — propagador completo | ninguna | 4 |
| **M2+L** — M2 más el lado de límite | — | 4 + núcleo `G_L` |

**No se decide** integrarlo en el filtro (eso es la v3.3), ni ingerir `@depth`, ni retirar la
cadena EMD/HHT, ni borrar `ω_m,max`, `γ_ω`, `C_ON`, `C_OFF`, `W` ni `K`. `Micelio.py` **no se
toca** en toda la tanda.

**M0 es el nulo, no un candidato.** Cualquier otro modelo tiene que batirlo fuera de muestra
para existir. Que gane M0 es un resultado, no un fracaso de la medición.

---

## 1. Compuerta de datos (§2.1) — se evalúa ANTES y se obedece

Se ejecuta `python captura_larga.py --resumen --dir=telemetria/captura_v32` y **se pega su
salida literal en la primera sección del reporte**, pase o no pase.

### 1.1 ⚠ El umbral va en TICKS, no en horas

El §2.1 de la orden pide "≥ 8 h". **Eso no acota nada en un proyecto que declaró el reloj de
transacciones primario**, con `ν` variando por un factor 20 entre capturas: 8 h a 102 tx/s son
2.9 M ticks y a 5.6 tx/s son 161 k.

El umbral correcto sale de encadenar los requisitos que el propio documento ya impone:

```
bloques efectivos ≥ 15   y   longitud ≥ 5·embargo    →   prueba UTILIZABLE ≥ 75·embargo
la banda de embargo se DESCARTA de la prueba         →   0.20·N − embargo ≥ 75·embargo
                                                     →   N ≥ 380·embargo
```

Con el piso `embargo = 1 950` ticks eso son **≥ 741 000 ticks continuos y limpios**. Se
**recalcula** si `H*_ticks` sale mayor que el piso, porque entonces el embargo crece y el
requisito con él. El número se **deriva** en `migracion_v32.ticks_de_compuerta()`, no se escribe
como literal.

⚠ **Corrección sobre la primera versión de este documento**, que decía `375·embargo = 731 250`.
Esa cuenta olvidaba que **el embargo se descarta del conjunto de prueba en lugar de repartirse**,
y daba **14 bloques, no 15**. Lo caza el control 6 de `migracion_v32.py`, que además comprueba la
contraprueba: con 375·embargo salen 14. Es exactamente el defecto que este documento obliga a
marcar —presentar 12 bloques como 15— y apareció en su propio umbral.

Lo que eso significa en horas depende del régimen, y por eso no es un umbral horario:

| `ν` | horas necesarias |
|---|---|
| 3.35 tx/s (el de la captura al lanzarla) | **61.4 h** |
| 12.4 tx/s (medido a los 30 min) | 16.6 h |
| 39 tx/s | 5.3 h |
| 102 tx/s | 2.0 h |

⚠ **Consecuencia operativa declarada ahora: la captura de 48 h puede no bastar.** Si al terminar
no hay **741 000** ticks limpios, la compuerta **no pasa** y hay que relanzar, por mucho que el
reloj de pared diga 48 h. No se rebajará el umbral para poder concluir.

### 1.2 La tabla

| requisito | umbral |
|---|---|
| tramo continuo limpio | **≥ 380 · embargo ticks** (≥ 741 000 con el piso) |
| hueco temporal interno máximo | **< 300 s** |
| `tr_maker` y `q` persistidos | obligatorio |
| `bookTicker` con `b`, `B`, `a`, `A` | obligatorio |
| filtro `tick_valido` aplicado en origen | obligatorio |

**Si no pasa, no hay decisión y no se ejecuta nada de §3 en adelante.** La única acción
admisible es capturar más. Un reporte que presente conclusiones sobre datos que no pasan la
compuerta es un entregable rechazado, y esto se escribe aquí para que no haya margen de
interpretación después.

**Se declara ahora que 12 bloques de bootstrap no son 15.** Si el número de bloques efectivos
queda por debajo de 15, todos los intervalos se marcan **subpotenciados** en el propio reporte,
no en una nota al pie. Es el defecto exacto que la Adenda A de la v3.1 dejó abierto.

---

## 2. Banco de pruebas — congelado aquí

### 2.1 Reloj y partición

- **Reloj primario: transacciones.** Todo `τ` en ticks; todo resultado **además** en segundos
  con la `ν` del tramo.
- **Corte temporal 60 / 20 / 20**: entrenamiento, validación, prueba. Sin validación cruzada
  aleatoria — mezclaría futuro con pasado.
- **El conjunto de prueba se abre UNA sola vez**, al final, con todos los modelos ya fijados.
  Si hace falta volver a él, el resultado deja de ser fuera de muestra y así se reporta.

### 2.2 ⚠ El embargo: `H*` en segundos NO es una constante

El §2.2 de la orden propone `H*·ν` con `H* = 50 s`. Pero la sesión 2026-08-08 (e) estableció
que `H*` pasa a `H*_ticks = (c/σ_tick)²` y que **su valor en segundos depende de `ν`**, que en
este proyecto varía por un factor 20. Tomar 50 s como constante reintroduciría el error que esa
sesión corrigió.

**Se declara:** el embargo es en **ticks**, calculado como

```
embargo = max( H*_ticks estimado SOLO sobre el conjunto de entrenamiento , 1950 ticks )
```

El piso de 1 950 ticks es el valor que da la orden (`50 s × 39 tx/s`) y se conserva como cota
inferior para que un `H*_ticks` pequeño no anule el embargo. **El número resultante se congela
antes de tocar validación**, se reporta explícitamente y se **verifica por test**, no por
comentario.

### 2.3 Métrica

- **Decide: `LL/N` fuera de muestra** sobre el conjunto de prueba.
- Intervalo por **bootstrap por bloques móviles**, longitud `≥ 5 · embargo`.
- `R²` se reporta como descriptivo y **no decide**.
- **Prohibido decidir con `LL` dentro de muestra.** M2 ⊃ M1 estrictamente, así que M2 no puede
  ajustar peor in-sample; usarlo para comparar es motivo de rechazo del entregable (§8, fallo 9).

### 2.4 Identidad de muestra

`assert` de que los índices del conjunto de evaluación son **idénticos** para M0, M1, M1′, M2 y
M2+L. El fallo clásico de esta comparación es evaluar M1 sobre la malla de `bookTicker` y M2
sobre la de `@trade` y atribuir a los modelos una diferencia que es de muestreo.

---

## 3. Observable primario: punto medio (§3.4)

**Todo el banco se corre sobre `p_mid = (bid + ask)/2`.** El precio de transacción se mantiene
como control secundario y **las dos columnas van siempre en el reporte**.

**Regla falsable declarada ahora:** si M1 (o cualquier modelo) **gana sobre precio de
transacción y pierde sobre punto medio**, lo medido es el rebote bid-ask y la conclusión es que
ese modelo **no aporta**. No se elegirá la columna que favorezca al modelo.

Justificación, con los números ya medidos: `R(1) = −0.0052` es rebote bid-ask (v3.1 §2) y
`s_eff = 0.2062 USD/BTC` ≈ 2 ticks (v3.1 §1). El spread mediano medido el 2026-08-08 sobre el
libro real es **0.10 USD = exactamente 1 tick**. Predecir que la siguiente transacción cae al
otro lado del spread es predecir algo real y **no operable**.

---

## 4. Alineación del libro (§3.1–§3.3)

### 4.1 La regla

Para cada transacción en `t`, el estado del libro utilizable es el del **último `bookTicker` con
timestamp estrictamente anterior** a `t`; `e_t` se construye con ése y el inmediatamente previo.
**Nunca interpolar. Nunca `ffill` hacia atrás.**

### 4.2 Test de fuga — dirección esperada declarada ANTES

Se ejecuta el ajuste de M1 **tres veces** y **las tres cifras van al reporte**:

Se ejecuta el ajuste de M1 **cuatro veces** y **las cuatro cifras van al reporte**:

| variante | expectativa declarada | qué prueba |
|---|---|---|
| alineación correcta | el valor a reportar | — |
| libro **adelantado** 1 actualización | **debe subir claramente** | prueba de humo |
| libro **retrasado** 10 actualizaciones | debe bajar | prueba de humo |
| **libro BARAJADO entre transacciones** | **debe caer al nivel de M0** | **el control con poder** |

- Si la fuga deliberada **no** mejora → el pipeline no está usando el libro y `e_t` es ruido.
  **Se para y se arregla el pipeline**; no se reporta el número del medio como si valiera.
- Si la correcta y la adelantada **empatan** → hay fuga en la correcta. Mismo tratamiento.
- **Si la barajada NO cae al nivel de M0**, la ganancia de M1 no viene del emparejamiento entre
  libro y transacción: viene de la marginal de `e_t` o de la propia especificación, y M1 no está
  midiendo lo que dice medir.

⚠ **Las dos primeras variantes son casi tautológicas y así se etiquetan.** Adelantar el libro una
actualización captura la actualización que causó **esa misma** transacción, así que va a subir
siempre; sólo demuestra que el código lee el campo. **La cuarta variante es el único control con
poder**, y es la que este proyecto ya ha usado con éxito dos veces: barajar destruyendo el
emparejamiento y conservando la marginal (el `C = 0.783` sobre barajados de la v2.2, y el
propagador con signos barajados de la v3.1 §2, factor 150–700).

### 4.3 ⚠ Esto es OFI-L1 APROXIMADO, y así se etiqueta en todo el reporte

`@bookTicker` es un *snapshot* con estrangulamiento: las cancelaciones entre actualizaciones son
invisibles, y son justo las que `e_t` pretende penalizar. El OFI real exige el flujo diferencial
de `@depth`, que hoy no se ingiere. **En ningún punto del reporte se escribirá "OFI" a secas.**

Diagnóstico obligatorio, con umbral declarado:

- **fracción de transacciones cuyo `q` excede la cantidad del mejor nivel del último snapshot**.
  Cada una prueba al menos una actualización perdida.
  **Si supera el 20 %, M1 se evalúa igualmente pero su resultado es una COTA INFERIOR**, no una
  medición, y así se rotula en la tabla.
- **reconciliación** de `ΔV_b` y `ΔV_a` con las transacciones observadas; el residuo es actividad
  de límite no observada. Se reporta como porcentaje.

---

## 5. Test A — el ajuste que responde la pregunta (§4)

### 5.1 Malla y elección

- `δ ∈ {0, 0.25, 0.5, 1}`, elegida **en validación**, nunca en prueba. La celda `δ = 0.5` se
  reporta **aunque no gane**: es la de MkII y es la pregunta del operador.
- `v` se normaliza por su **mediana del tramo** antes de elevar a `δ`. Sin eso `G₀` y `δ` quedan
  confundidos por la escala y sus errores estándar no son interpretables.
- Errores estándar por bootstrap por bloques, nunca analíticos.

### 5.2 Los dos contrastes que deciden la migración

⚠ **CORRECCIÓN: `β̂` NO es interpretable por sí sola, y el contraste del §4.2 de la orden tal
como está escrito no es válido.** Encontrado por el control positivo de `migracion_v32.py`
**antes** de tocar dato real: sobre una serie generada con `β = 0` verdadero, el ajuste devuelve
`β̂ = 5.0` (el borde) y **acierta igualmente la curva**. La razón es que con `τ₀` grande,
`(1 + τ/τ₀)^(−5) ≈ 1` en todo el rango observado, o sea **observacionalmente idéntico a
`β = 0`**. La pareja `(τ₀, β)` no está identificada por separado; lo que los datos determinan es
la **curva `G(τ)` sobre `[0, K]`**.

Se declara como cantidad identificada el **decaimiento relativo sobre el rango ajustado**:

```
D = G(K) / G(0)          D ≈ 1  ->  impacto permanente  ->  M1 / MkII
                         D < 1  ->  impacto transitorio ->  propagador
```

Verificado: con verdad `β = 0.35` (`D` verdadero 0.6808) el ajuste da `D = 0.6994`; con verdad
`β = 0` (`D = 1`) da `D = 0.9835`, pese a devolver `β̂ = 5.0` en ese segundo caso.

| contraste | regla declarada |
|---|---|
| **¿impacto permanente?** | IC bootstrap al 95 % de **`D = G(K)/G(0)`**. Si **contiene 1** → MkII está bien especificada y **se prefiere por parsimonia** |
| **respaldo obligatorio** | el contraste anidado **M1 contra M2 fuera de muestra** (§9 y §10). Es el que decide si los dos discrepan |
| **¿`δ = 0.5`?** | si `δ = 0.5` es el máximo en validación **o está dentro de su IC** → MkII acierta también aquí |

**`β̂` y `τ̂₀` se reportan, pero marcados como NO identificados individualmente.** `τ₀` sigue
siendo el sucesor honesto de `ω_m` como *escala*, pero su valor puntual sólo es interpretable
junto con `β̂`, y su IC conjunto es lo que hay que mirar.

**La respuesta a "¿migramos a OFI?" es el resultado de estos contrastes.** No se decide por
preferencia arquitectónica, ni por elegancia, ni por lo que diga este documento.

### 5.3 Sobrepaso: `ω_G` contra nulo simulado

El estadístico `R(final)/R(pico)` **queda retirado**: es degenerado bajo el nulo (p5 = −560,
mediana = −11 cuando el pico es ruido; v3.1 §2). Se sustituye por

```
G(τ) = G₀ · (1 + τ/τ₀)^(−β) · cos(ω_G·τ + φ)      (M2-osc, 5 parámetros)
```

y se contrasta `ω_G = 0` por **razón de verosimilitud con nulo SIMULADO**: generar bajo M2
ajustado con `ω_G = 0`, reajustar M2-osc, construir la distribución de `2ΔLL` y comparar.

⚠ **Prohibido usar la tabla χ².** Con observaciones dependientes `2ΔLL ~ χ²` no se cumple, y es
la **tercera aparición** del mismo patrón en este proyecto (χ² sobre NIS en v2.1; `(1−γ)/2` en
v3.1 §3.2; y aquí). Todo contraste de verosimilitud de esta tanda va contra nulo simulado.

⚠ **`ω_G` NO es `ω_m` resucitada**, y el reporte lo dirá con estas palabras. `ω_m` era una
frecuencia del *precio*, extraída por EMD, sin nulo — y la sesión 2026-08-08 (f) mostró que el
período que devuelve la descomposición es el de la ventana (`T = 2L/k`, error 0.0000). `ω_G` es
una frecuencia de la *respuesta al impulso del flujo*, con estimador propio, error estándar y
nulo simulado.

- `ω_G ≠ 0` → el sistema es un oscilador **forzado**, y la premisa de la v3.0 sobrevive en la
  única forma en que podía sobrevivir.
- `ω_G = 0` → el propagador es monótono, y `ω_m,max` y `γ_ω` se borran de
  `constantes_micelio.py`.

### 5.4 Residuo

Ljung-Box sobre el residuo de **cada** modelo, a rezagos hasta `2 · embargo`.

**Declarado ahora: un modelo con residuo estructurado no gana aunque tenga mejor `LL/N`.** Eso
significa que su ventaja es de ajuste local y no de especificación. Umbral: `|ρ₁| < 0.05` con
decaimiento a cero, el mismo que la v3.1.

---

## 6. Test B — ¿aporta el lado de límite? (§5)

La pregunta **no** es "OFI contra signo": `ε` ya es el lado taker del OFI y ya se midió con
factor 150–700 sobre su control. Lo que M1 añade es la actividad de **límite**.

- `e_t = e_t^(taker) + e_t^(límite)`, con el taker reconstruible desde `@trade`.
- Si `LL/N(M2+L) − LL/N(M2)` no supera el margen del §8, **el lado de límite no aporta**, y
  `@depth` y el anti-*spoofing* de MkII son coste sin retorno. **Ese resultado se publica igual.**

### 6.1 Sesgo de Jensen por `√|e|` — guarda declarada

MkII usa `sgn(e)·√|e|`. Con `e` medido con error —y por §4.3 se mide con error— la concavidad de
`√` sesga hacia **sobreestimar** el impacto de los `e` pequeños, y `e ≈ 0` es la mayoría de las
observaciones.

**Guarda:** simular con `e` real más ruido de magnitud calibrada por §4.3, ajustar, y reportar el
sesgo en `G₀`. **Si el sesgo es del orden de `G₀`, el estimador no vale** y se usa `f(v) = v^δ`
con `δ` libre — que es lo que M2 hace de todas formas.

### 6.2 `c_t` como hipótesis, no como especificación

| variante | pregunta |
|---|---|
| `G₀` constante | referencia |
| `G₀ → Y·σ_t` | ¿la ganancia viene sólo del escalado por volatilidad? |
| `G₀ → Y/√V_best` | ¿aporta la cantidad en el mejor nivel? |
| `G₀ → Y·σ_t/√V_best` | ¿aportan las dos? |

⚠ **Se etiqueta `V_best`, nunca "profundidad".** `bookTicker` da la cantidad en el mejor nivel;
la profundidad real necesita `@depth`. Llamarlo profundidad sería el deslizamiento de nombres que
ya produjo la colisión `C_max`/`I_max` en la v1.1.

**Declarado ahora:** si la variante `Y·σ_t` sola captura toda la ganancia, **M1′ no está midiendo
microestructura: está reescalando por volatilidad**, y eso ya lo hace `ρ_k` en el filtro. Se
reportará con esas palabras.

---

## 7. Test C — ¿en qué reloj vive el propagador? (§6)

`τ₀` es el número que sustituye a `ω_m`, así que en qué reloj vive es decisorio.

Regresión de `log τ₀` contra `log ν` por bloques:

| pendiente | lectura |
|---|---|
| ≈ 0 | el impacto decae en **tiempo de transacciones** (`τ₀` en ticks estable) |
| ≈ 1 | el impacto decae en **tiempo de pared** (`τ₀/ν` en segundos estable) |

**Se reporta la pendiente con su IC, no un veredicto.** `τ₀` va en ticks **y** en segundos.

⚠ **Confundidor declarado:** `ν` correlaciona con volatilidad y con sesión. Se incluye el
regresor de régimen del §8 como control. Sin él, la pendiente mide otra cosa.

### 7.1 ⚠ Sólo DENTRO de la captura de 48 h — nunca entre capturas

Este test se corre **exclusivamente sobre bloques de `captura_v32`**. Comparar `τ₀` entre
capturas distintas confundiría el estimador con los **cambios de pipeline**: la propia lista de
limitaciones del §12 reconoce dos omisiones de campos entre versiones (`tr_maker` hasta la v3.1,
`B`/`A` hasta la v3.2), y una captura sin esos campos no es la misma medición hecha en otro
régimen, es otra medición.

Las 48 h dan variación natural de `ν` **por sesión** —Asia, Europa, América— dentro de un mismo
código, que es el diseño limpio y el único admisible aquí.

**Consecuencia: el §7 NO se ejecuta con un tramo parcial.** Espera a la captura completa aunque
la compuerta del §1 ya haya pasado con menos.

---

## 8. Régimen — congelado antes de mirar

Idéntico al preregistro de la v3.1, para que las dos tandas sean comparables:

- **Primario: volatilidad realizada** por bloques.
- **Secundario 1: magnitud del desbalance de flujo** (`|Σ ε·f(v)|` por bloque).
- **Secundario 2: sesión** (Asia / Europa / América), exógena y sin circularidad.

Cortes por **terciles de la distribución del primario sobre todo el registro**, no por
inspección. **No se cambian después.**

Se añade como **exploratorio y no confirmatorio**, por venir de la sesión 2026-08-08 (f):
**φ′ = ticks/BTC**, que mostró asociación positiva con la volatilidad realizada en 16/16 ventanas
(p de signos 3.05e-05, +18 % a +41 % entre terciles). No entra en ninguna regla de decisión de
esta tanda.

---

### 8.1 ⚠ Estado de SSA, y por qué NO es un cuarto intento de resucitar `ω_m`

Entre la v3.1 y esta tanda hubo dos sesiones que no están en la orden de trabajo, y de ellas
salen tres cifras que este documento cita. Se declaran aquí para que nadie tenga que
reconstruirlas, y sobre todo para vallar SSA.

**Sesión 2026-08-08 (e)** — la firma de volatilidad debe medirse en tiempo de ticks; en reloj de
pared el estimador tiene sesgo propio (−0.0796 sobre barajados, contra −0.0054 en ticks) y los
muestreos solapado y no solapado discrepan **en signo**. De ahí `H*_ticks = (c/σ_tick)²` y el
§2.2 de este documento.

**Sesión 2026-08-08 (f)** — cambio de herramienta pedido por el operador: SSA en lugar de EMD.
El objetivo declarado era **una ventana de toma de datos con logeo completo**, no un modelo. Lo
que salió fue **negativo para `ω_m`, no un rescate**:

- los autovectores del precio real son **los armónicos de la ventana**, `T = 2L/k`, con error
  mediano **0.0000** al escalón, idéntico sobre un paseo aleatorio y sobre incrementos barajados,
  e invariante al cambiar `ν` por un factor 18;
- el primer autovector se lleva 84–99 % de la energía y los pares detectados **0.00–0.06 %**,
  contra **85.9 %** en un control positivo con períodos conocidos;
- el período y `φ′` son **independientes**: `L` explica el 79 % de la varianza de `T`.

O sea que SSA **no reabrió `ω_m`: mostró que el período que devuelve una descomposición es del
instrumento**. Lo hizo con el nulo que la EMD nunca tuvo (Monte Carlo SSA) y con dos límites de
ese nulo documentados —inservible sobre series casi blancas, demasiado severo con señales fuertes.

**Se declara para esta tanda:**

1. **Ninguna cantidad derivada de SSA entra en ninguna regla de decisión de la v3.2.** Ni `L`, ni
   `T`, ni la ortogonalidad, ni los pares.
2. `ssa.py` **no se importa desde `Micelio.py`** y no se importará en esta tanda.
3. `φ′` entra sólo como **regresor exploratorio** (§8), nunca confirmatorio.
4. Si algún resultado del §5.3 se leyera como "el ciclo vuelve", el §5.3 ya lo bloquea: `ω_G` es
   una frecuencia de la **respuesta al impulso del flujo**, con nulo simulado, no del precio.

La suite de `ssa.py --autotest` (11/11) se cita en el §13 sólo como comprobación de que el
código de análisis de la sesión (f) sigue en verde, no como evidencia sobre el mercado.

---

## 9. Regla de decisión (§9.1) — en este orden, y se para en el primer fallo

1. **Compuerta de datos** (§1 de este documento). Si no pasa → **no hay decisión**; capturar más.
2. **`LL/N`(mejor modelo) > `LL/N`(M0)** fuera de muestra, con IC bootstrap que **excluya 0**.
   Si falla → **gana M0**: la capa predictiva se retira y el proyecto pasa a la lista de
   supervivientes del §6.2 de la v3.0 (reloj de ticks, ingesta, riesgo, Loeper, NMPC como
   ejecutor).
3. **Criterio económico** — reformulado, ver §9.1. Si falla → **gana M0**.
   **Ésta es la única comparación con contenido económico y ningún `LL/N` la sustituye.**
4. **M2 contra M1** por el margen del §10. Si M2 no lo supera → **gana M1 (MkII)**, por parsimonia.
5. **M2+L contra M2** por el mismo margen. Decide si hace falta `@depth`.

**El paso en el que se detuvo se indica explícitamente en el reporte.**

### 9.1 ⚠ El paso 3, reformulado — la media incondicional NO es el criterio

Tal como lo escribe el §9.1 de la orden —"`E[Δp | ℱ_t]` al horizonte `H*` supera `c(u)`"— el
criterio es **inejecutable y provoca un abandono falso**. `E[Δp | ℱ_t]` es una variable
aleatoria, no un número. Siendo el precio aproximadamente una martingala, el valor
**incondicional** de ese pronóstico está cerca de cero **por construcción**, así que el paso 3
fallaría siempre, hubiera o no estructura explotable. Y el §11 abandona la hipótesis entera
sobre ese fallo, de forma irreversible por diseño. Es el error más caro disponible en todo el
procedimiento.

Lo que decide económicamente **no es la media del pronóstico, es la masa de su distribución más
allá del coste**. Se declara ahora, con todo estimado **fuera de muestra** sobre el conjunto de
prueba:

Sea `μ̂_t = E[Δp | ℱ_t]` al horizonte `H*`, y `c(u)` el coste de ida y vuelta del esquema de
ejecución declarado.

| # | condición | umbral |
|---|---|---|
| **3a** | **decil superior** de `\|μ̂_t\|` | `q90(\|μ̂_t\|) ≥ 1.5 · c(u)` |
| **3b** | fracción de ticks con `\|μ̂_t\| ≥ 1.5 · c(u)` | `frac ≥ 20 · embargo / N_prueba` |

- **3a** es la señal: cuando el modelo habla fuerte, ¿habla por encima del coste?
- **3b** es la potencia: esa fracción debe corresponder a **≥ 20 oportunidades NO solapadas** en
  el conjunto de prueba. Dos oportunidades separadas por menos que el embargo son la misma
  oportunidad contada dos veces. El umbral se escribe en función de `N_prueba` para que **se
  autoescale** con el tamaño de la captura en vez de ser un número que envejece.
- **1.5 × `c(u)`** es el mismo "con margen" de la v3.1, conservado para que las dos tandas usen
  la misma vara.

### 9.2 ⚠ Guarda contra el abandono por falta de potencia

**Si 3a pasa y 3b falla, el resultado es NO DECIDIBLE, no abandono.** Que no haya suficientes
oportunidades independientes para medir es una carencia de datos, no evidencia de ausencia de
estructura. En ese caso la única acción admisible es capturar más, exactamente como en la
compuerta del §1.

El §11 (abandono) **sólo se activa si 3a falla con 3b satisfecha**, es decir: hay potencia
suficiente para haber detectado la señal, y la señal no está.

Esta guarda se escribe aquí porque el abandono es irreversible y este proyecto ya anuló dos
veredictos por emitirlos sobre datos que no los sostenían.

---

## 10. El margen, declarado ahora (§9.2)

### 10.1 ⚠ El §9.2 de la orden dice BIC y escribe otra cosa

BIC por observación es `k·ln(N)/(2N)`; la orden escribe `k/(2N)`, que es **medio AIC**. Con
`N ≈ 10⁶` la diferencia es un factor `ln N ≈ 14`.

Y hay un problema mayor debajo: **el IC sale de 15–138 bloques mientras la penalización usaría un
`N` de ticks**. Mezclar el `N` ingenuo con la incertidumbre de `N_eff` es incoherente, y las
elecciones plausibles difieren aquí en **cinco órdenes de magnitud** (medido, con `k = 3`):

| penalización | `N = 2e5`, `N_eff = 15` | `N = 1e6`, `N_eff = 138` |
|---|---|---|
| `k/(2N)` — lo que escribe la orden | 1.00e-05 | 1.50e-06 |
| `k·ln(N)/(2N)` — BIC con `N` ingenuo | 1.19e-04 | 2.07e-05 |
| **`k·ln(N_eff)/(2N)`** — BIC coherente | **2.71e-05** | **7.39e-06** |
| `k·ln(N_eff)/(2N_eff)` | 0.271 | 0.054 |

Con `N` de ticks la parsimonia **no vota**; con `N_eff` de bloques **no vota nada más**.

### 10.2 Lo que se declara

Preferir el modelo más complejo exige **las dos** condiciones:

**(a) Estadística.** El límite inferior del IC bootstrap al 95 % de `ΔLL/N` debe superar

```
k · ln(N_eff) / (2N)
```

con `k` la diferencia de parámetros, `N` el número de observaciones de prueba y `N_eff` el
número de **bloques efectivos** de bootstrap. Es BIC con el tamaño muestral **efectivo** —el
único coherente con observaciones dependientes— expresado por observación para que sea
comparable con la métrica.

**(b) Económica, y es la que decide en la práctica.** El modelo complejo debe mejorar el
estadístico del §9.1 en al menos **`0.1 · c(u)`**:

```
q90(|μ̂_t|)_complejo  −  q90(|μ̂_t|)_simple  ≥  0.1 · c(u)
```

**Se declara explícitamente que (a) pasará casi siempre con `N` grande y que (b) es el criterio
vinculante.** Eso es deliberado: la decisión final de este proyecto es económica, así que el
desempate debe estar en escala económica y no en nats por observación, donde cualquier elección
de penalización es defendible y ninguna es interpretable.

---

## 11. Criterio de ABANDONO (§9.3)

Sin cambios respecto al §4.2 de la v3.1, que sigue vigente y **este documento no relaja**.

Si los pasos 2 o 3 del §9 fallan sobre datos que **sí** pasan la compuerta, la hipótesis de
**estructura explotable a `H*`** queda abandonada.

> "No hay otra escala a la que retirarse: ésta es donde está el dinero."

---

## 12. Limitaciones conocidas ANTES de medir

Se anotan ahora para que no se presenten después como matices:

1. **El escalón de comisiones no se puede leer de la cuenta.** `/fapi/v1/commissionRate` es
   firmado y el Modo LECTURA no tiene credenciales. Se usan las tarifas públicas VIP 0 **marcadas
   como asumidas**, y el criterio correspondiente queda **no cumplido**, no fingido. Idéntico a la
   v3.1.
2. **`bookTicker` no es OFI real** (§4.3). Todo resultado de M1, M1′ y M2+L es "OFI-L1 aproximado".
3. **Toda captura anterior al 2026-08-08 se descarta**, no se limpia a posteriori: el feed emitía
   trades con `p = 0` sin filtrar y el factor en el incremento máximo era de 5 500.
4. **Las capturas anteriores a `captura_v32` no tienen `B` ni `A`**, así que ninguna sirve para
   M1, M1′ ni M2+L. Es una omisión propia, la segunda de la misma familia tras `tr_maker`.
5. **`H*` en segundos no es una constante** (sesión 2026-08-08 e). De ahí el §2.2 de este
   documento.
6. **El nulo espectral de la v2.2 se midió sobre datos sucios** y no se cita como cerrado.
7. **Ningún resultado de esta tanda dice nada sobre la hipótesis (B)** —escala de decenas de
   minutos— que sigue sin decidir desde la v2.2.
8. **La captura de 48 h puede no alcanzar la compuerta.** A la `ν` con que arrancó (3.35 tx/s)
   harían falta 61.4 h para los 741 000 ticks. Si al terminar no llega, se relanza; no se rebaja
   el umbral.
9. **`e_t` se estima con `ε` que incluye transacciones de todo el mercado.** En Modo LECTURA no
   operamos, así que no hay impacto propio que descontar; la guarda del §7.3 de la orden
   (separación propio/ajeno por `orderId`) es de la v3.3, no de ésta, y se anota para que no se
   olvide al integrar.

---

## 13. Transversales que se verifican en el reporte

- [ ] Suite en verde: **56/56** de `tests_v13.py` + **11/11** de `ssa.py --autotest` + los nuevos.
- [ ] `Micelio.py` **sin cambios** (`git diff --stat` en el reporte).
- [ ] **Ningún literal nuevo** en la ruta del filtro; `grep` de `Δn` y de constantes mágicas
      incluido en el reporte. Si aparece un `Δn` fijo, se ha reabierto el agujero que el §7 de la
      v3.1 cerraba.
- [ ] Todo contraste de verosimilitud **contra nulo simulado**, ninguno contra tabla χ².
- [ ] Salida literal de `--resumen` de la compuerta, pase o no pase.
- [ ] Número de bloques de bootstrap **reportado explícitamente**.
