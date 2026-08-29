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

⚠ **El piso es un literal heredado, y hay que resolverlo cuando `H*_ticks` exista.** 1 950 sale
de `ν = 39 tx/s`; a la `ν ≈ 11 tx/s` de la captura eso son **174 s**, no 50. Si `H*_ticks`
estimado sobre entrenamiento sale del orden de 400, **el piso está haciendo todo el trabajo y
multiplica por ~5 el requisito de datos sin justificación viva** — y el requisito de datos es lo
que decide si hay decisión o no. No es un error de corrección, pero es justo la clase de
constante que el `grep` del §13 busca. **Declarado ahora: en cuanto `H*_ticks` esté medido, o se
justifica el piso con un argumento propio o se retira**, y el reporte dice cuál de las dos.

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

⚠ **SEGUNDA CORRECCIÓN: `D` tampoco tiene nulo, y el umbral 1 es inservible.** El estimador de
`D` está sesgado, y **el signo del sesgo no es estable**: con `N = 40 000`, `K = 120` y signos iid
la mediana bajo `D = 1` verdadero sale **0.9835**; con `N = 20 000`, `K = 60` y `γ = 0.3` sale
**1.0023**. Medida la tasa de error de la regla ingenua sobre impacto **permanente** verdadero:

```
"D_hat < 1  ->  transitorio"   declara transitorio en el 33 % de los sorteos
```

**33 % de falsos positivos contra un 5 % nominal**, y el sesgo empuja hacia el modelo **más
complejo** justo en la frontera de decisión. Sería la **cuarta** aparición del patrón que el
§5.3 de este mismo documento prohíbe.

**Se declara:** el umbral se toma del **cuantil 5 % de la distribución simulada de `D̂` bajo
`D = 1`**, con `N`, `K`, `γ` y `f(v)` **emparejados con la muestra real**, nunca del valor
asintótico 1. Con la configuración del control salió `q05 = 0.9817`. La **curva de sesgo** se
reporta sobre una malla de `D` verdadero (1.00, 0.90, 0.80, 0.70, 0.50, 0.30), no sobre dos
puntos. Implementado en `migracion_v32.nulo_de_D` y `curva_de_sesgo`.

⚠ **TERCERA CORRECCIÓN: `D` no separa "decae a un suelo" de "decae a cero", y `G(∞)` es
justamente lo que el §5 de la v3.1 designó referencia móvil del sistema.** El núcleo original
`(1 + τ/τ₀)^(−β)` tiene `G(∞) = 0` **por construcción**: todo el impacto es transitorio se mida
lo que se mida. La sospecha ya estaba anotada en la v3.1 §3 y se quedó sin arreglar.

**Se declara el núcleo con suelo permanente:**

```
G(τ)/G₀ = f_∞ + (1 − f_∞)·(1 + τ/τ₀)^(−β)        f_∞ = G(∞)/G(0)
```

`f_∞` es la **fracción permanente**, que es la cantidad con contenido económico y la que la v3.3
necesitará de todos modos. Verificado: con verdad `f_∞ = 0.60`, el ajuste **sin** suelo no puede
representarlo y el ajuste **con** suelo devuelve **0.572**. `M1` queda anidado como `f_∞ = 1`,
más limpio que por `β = 0`.

Se reporta además el **perfil `D` a `K/4`, `K/2` y `K`** para ver si la caída se aplana. ⚠ Los
intervalos no son igual de anchos —el segundo abarca el doble de rezago—, así que el criterio es
que caiga **menos** en el tramo ancho; exigir un factor 2 sobre las caídas crudas fallaba sobre
un suelo verdadero del 60 %.

⚠ **`K` se congela ANTES: `K = 2 · embargo`.** `D` depende del rango en que se evalúa, así que
"`D ≈ 1`" no significa nada hasta que `K` esté declarado.

### 5.2.bis Potencia, no sólo talla — y por qué eso cambia cómo se lee el paso 4

Una malla que salta de 1.00 a 0.90 da la **talla** del contraste y no su **potencia**, y la
asimetría no es inocua: con potencia baja el contraste diría "permanente" por defecto, y el paso
4 de la regla de decisión (*si M2 no supera el margen, gana M1*) convertiría esa falta de
potencia en **una victoria de MkII que nadie declaró**.

Medido sobre la malla fina (`N = 20 000`, `K = 60`, `γ = 0.3`, 24–40 sorteos por punto):

| `D` verdadero | 1.000 | 0.990 | 0.970 | 0.950 | 0.900 | 0.800 | 0.700 | 0.500 |
|---|---|---|---|---|---|---|---|---|
| `D̂` mediana | 0.9950 | 0.9848 | 0.9684 | 0.9487 | 0.9034 | 0.8032 | 0.6986 | 0.5006 |
| **potencia** | 0.00 | **0.00** | 0.62 | **0.96** | 1.00 | 1.00 | 1.00 | 1.00 |

`q05 = 0.9698`, `sd` del nulo `= 0.0158`, **talla = 0.05** (exactamente la nominal).

**Desviación mínima detectable con potencia ≥ 0.80: `D = 0.95`**, o sea un **5 %**. Una
extrapolación por normalidad desde `q05` y la mediana sugería `D ≈ 0.989`; la medición la
desmiente por un factor 5. Es la razón de medir la potencia en vez de inferirla.

⚠ **Consecuencia declarada ahora: el paso 4 NO es una prueba económica y no se leerá como tal.**
Con esta resolución, rechazar `D = 1` puede ocurrir con un transitorio económicamente
irrelevante. **El paso 3 sigue siendo la única compuerta con dinero detrás.**

### 5.2.ter `f_∞` se estima SECUENCIALMENTE, y tiene su propia malla de sesgo

`f_∞` hereda la degeneración por el otro lado: **con `β → 0` el núcleo tiende a 1 sea cual sea
`f_∞`**, así que `f_∞` no está identificado cuando no hay decaimiento.

**Se declara el orden:** `f_∞` se estima y se reporta **sólo condicionado a haber rechazado
`D = 1`**. Si no se rechaza, `f_∞` no se reporta, porque no hay contenido que reportar.

Su sesgo es **un orden de magnitud peor que el de `D`**, y el signo es el peligroso —
infravalorar el permanente sobrestima el transitorio, o sea predice más reversión de la que hay,
que es sobreoperación en la magnitud que el §5 de la v3.1 designó referencia móvil:

| `f_∞` verdadero | 0.00 | 0.20 | 0.40 | 0.60 | 0.80 |
|---|---|---|---|---|---|
| `f̂_∞` mediana | 0.0115 | 0.1782 | 0.3713 | 0.5820 | 0.8130 |
| sesgo relativo | — | **−10.9 %** | **−7.2 %** | −3.0 % | +1.6 % |

⚠ **Nota de método, y es la razón de fondo del nulo simulado.** Con M1 anidado como `f_∞ = 1`,
el nulo está en la **frontera del espacio de parámetros**, donde la asintótica de la razón de
verosimilitud falla **incluso con observaciones independientes** (Self & Liang 1987). Con
observaciones dependientes, doblemente. Que el nulo sea simulado deja de ser buena práctica y
pasa a ser **obligatorio**; no es la disciplina general del proyecto, es este caso concreto.

### 5.2.quater ⚠ `q05 = 0.9698` es PROVISIONAL — lo que se congela es el procedimiento

Circulan tres configuraciones con tres resultados distintos: `0.9835` (`N=40 k, K=120`, iid),
`1.0023` (`N=20 k, K=60, γ=0.3`) y `0.9950` (la malla fina). **Si el cuantil saliera de una y la
curva de sesgo de otra, el umbral y la corrección serían incoherentes entre sí.**

**Se declara:** lo congelado aquí es el **procedimiento**, no el número. El `q05` se **regenera
cuando la captura cierre**, emparejado a la `γ̂`, la `ν` y la `K` medidas sobre el conjunto de
entrenamiento. Y `γ` se estima del propio dato, así que su incertidumbre se propaga **simulando
sobre la distribución bootstrap de `γ̂`**; si no, el umbral queda condicionado a un puntual y su
cobertura no es la declarada.

| contraste | regla declarada |
|---|---|
| **¿impacto permanente?** | `D̂ ≥ q05` de la distribución simulada bajo `D = 1`, regenerada con `γ̂`, `ν` y `K` medidas → MkII está bien especificada y **se prefiere por parsimonia**. Nunca contra el umbral 1 |
| **¿cuánto es permanente?** | IC bootstrap de **`f_∞`**, **sólo si se rechazó `D = 1`**, más el perfil `D` a `K/4`, `K/2`, `K` |
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

⚠ **CUARTA CORRECCIÓN: la regresión de `log τ̂₀` contra `log ν` que pide el §6 de la orden es
INEJECUTABLE.** Si `τ̂₀` no está identificada individualmente —y no lo está, ver §5.2—, regresar
su logaritmo no mide nada. El test quedaba inservible sin que se notara.

**Se declara el sustituto, sobre el funcional que SÍ está identificado.** Se computa `D` dos
veces por bloque: una con el rezago tope **fijo en ticks** y otra **fijo en segundos**
(`K_b = round(T · ν_b)`), y se compara la dispersión entre bloques:

| lectura | firma |
|---|---|
| decae en **tiempo de transacciones** | `D` a `K` ticks estable; `D` a `T` segundos varía con `ν` |
| decae en **tiempo de pared** | al revés |

Sin `τ₀`, sin regresión. Implementado en `migracion_v32.reloj_del_propagador`, y validado **por
los dos lados**, que es lo que lo hace un diagnóstico y no un artefacto:

| control | generado con | `sd` de `D` a `K` ticks | `sd` de `D` a `T` segundos | lectura |
|---|---|---|---|---|
| directo | núcleo fijo **en ticks** | **0.0218** | 0.1284 | tiempo de transacciones ✔ |
| **espejo** | núcleo fijo **en segundos** | 0.1201 | **0.0123** | tiempo de pared ✔ |

⚠ **El control espejo es obligatorio y por poco no lo pongo.** Sin él no se distingue "el
diagnóstico detecta el reloj" de "`D` a rezago fijo en ticks es mecánicamente más estable porque
`ν` no entra en su definición". Es exactamente el patrón —control positivo sin su negativo— que
las cinco enmiendas anteriores vinieron a cerrar, repetido una vez más.

**Se reportan las dos dispersiones y las dos series de `D`, no un veredicto.**

Y hay palanca real: `ν` se movió de **3.35 a 12.4 tx/s dentro de la misma captura** en la primera
hora, así que el contraste no necesita comparar capturas de versiones distintas del pipeline.

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

### 11.1 ⚠ ENMIENDA — el abandono se acota a la ESCALA DE SEGUNDOS

**Redactada el 2026-08-09, ANTES de que exista ningún resultado del §9** y sin haber ajustado un
solo modelo sobre `captura_v33`. Se enmienda ahora precisamente porque hacerlo después sería
mover la portería.

La redacción anterior —"la hipótesis de estructura explotable a `H*` queda abandonada", con la
coletilla *"no hay otra escala a la que retirarse: ésta es donde está el dinero"*— **era correcta
mientras se creía que `H*` era el horizonte operativo. No lo es.**

`H*` está definido como el punto donde el movimiento típico **iguala** el coste. O sea: es el
horizonte en el que habría que acertar el movimiento **entero** para no perder. Eso no es un
punto de operación, es un **suelo por debajo del cual operar es imposible**. La v3.1 §1.3 lo
convirtió en "el número movible" que fija el horizonte del propagador y del NMPC, y ahí está el
error de lectura.

El R² necesario para que la ventaja supere el lastre cae con `1/H`. Con `σ₁ = 3.77 USD/BTC·s^½`
—el plateau de la firma de volatilidad de la v3.1, no el despeje `c/√H*`, que es tautológico
porque `H* := (c/σ₁)²`— y un lastre de comisiones (37.65) más selección adversa (5.78) medida en
la sesión de ejecución pasiva, operando sólo el decil superior de la señal:

| horizonte | `σ(H)` | lastre | **R² requerido** |
|---|---|---|---|
| **103 s = `H*`** | 38.3 | 43.4 | **41.8 %** |
| 10 min | 92.3 | 43.4 | **7.2 %** |
| 1 h | 226.2 | 43.4 | **1.2 %** |
| 4 h | 452.4 | 45.8 | 0.33 % |
| ~2.3 días | 1 680.6 | 86.1 | 0.09 % |

(financiación del perpetuo a partir de 1 h, ~6.3 USD/BTC por cada 8 h; es lo que crea un óptimo
interior en vez de "cuanto más largo, mejor")

**Un R² del 42 % sobre dirección de precio es fantasía.** Ésa es la razón *estructural* por la
que el paso 3 puede fallar, y **no tiene nada que ver con que haya o no señal**.

**Lo que se enmienda, y sólo esto:**

> Si los pasos 2 o 3 del §9 fallan sobre datos que sí pasan la compuerta, queda abandonada la
> hipótesis de **estructura explotable a escala de SEGUNDOS**. Eso **no** decide nada sobre
> escalas de minutos a horas, que pasan a ser la línea principal —la hipótesis (B) de la v2.2,
> sin decidir desde entonces y despriorizada justo cuando se derivó `H*`.

**Lo que NO se relaja:** los criterios de éxito, el margen, la compuerta, el orden de la regla de
decisión y la obligación de reportar el paso en que se detuvo siguen exactamente igual.

⚠ **Reservas de la aritmética de arriba, para que no se cite como más de lo que es:** supone
normalidad conjunta, difusividad sostenida, operar sólo el decil superior y tarifas VIP 0. **La
difusividad está medida sólo hasta 8 192 ticks (~12 min a `ν = 11`)** — más allá es
extrapolación y hay que verificarla. Y que el R² *requerido* caiga no dice nada sobre el R²
*alcanzable*: a escala de días el flujo de órdenes predice poco y el juego es otro. Además ~150
operaciones al año hacen que estimar un Sharpe tarde años, y eso es un coste real de irse largo.

### 11.2 El criterio, con la enmienda aplicada

Sin más cambios respecto al §4.2 de la v3.1, que sigue vigente y **este documento no relaja**.

Si los pasos 2 o 3 del §9 fallan sobre datos que **sí** pasan la compuerta, la hipótesis de
**estructura explotable a escala de segundos** queda abandonada.

### 11.3 Lo que la v3.2 decide, y lo que ya no

La v3.2 **no se para**: permanente contra transitorio es una pregunta bien planteada y el
estimador está validado con 18 controles. Pero **baja la apuesta que tiene encima**. Cualquiera
que sea el resultado, decide la **forma del núcleo de impacto a escala de segundos**, y eso es
una capa de **ejecución** —cómo colocar una orden— no la capa de **decisión** —cuándo operar y en
qué dirección.

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

## 12.bis ⚠ REGISTRO DE ENMIENDAS — obligatorio en el reporte

Un preregistro que se edita sin dejar rastro del motivo no protege de nada, y este proyecto sabe
cómo se lee un veredicto emitido sobre criterios movidos a posteriori. **Toda enmienda se acredita
con: hash antes, hash después, la regla retirada, la que la sustituye, el control que la invalidó,
y la constancia de que la captura no se había tocado.**

| # | commit antes → después | regla retirada | regla que la sustituye | qué la invalidó | ¿se miró dato real? |
|---|---|---|---|---|---|
| 1 | `190edda` → `2ad0700` | «el IC de `β̂` contiene 0» | IC de `D = G(K)/G(0)` contiene 1 | control positivo de `migracion_v32.py`: con `β = 0` verdadero el ajuste devuelve `β̂ = 5.0` y acierta la curva; `(τ₀, β)` no están identificados por separado | **NO** |
| 2 | `2ad0700` → (esta) | compuerta `375·embargo = 731 250` | `380·embargo = 741 000`, derivado en `ticks_de_compuerta()` | el embargo se **descarta** del conjunto de prueba en vez de repartirse; con 375 salen 14 bloques, no 15 | **NO** |
| 3 | `2ad0700` → (esta) | «`D < 1` → transitorio» | `D̂ < q05` de la distribución simulada bajo `D = 1` | el estimador de `D` está sesgado, con **signo no estable**, y la regla ingenua da **33 % de falsos positivos** sobre impacto permanente verdadero | **NO** |
| 4 | `2ad0700` → (esta) | núcleo `(1+τ/τ₀)^(−β)` con `G(∞) = 0` | núcleo con suelo `f_∞ + (1−f_∞)(1+τ/τ₀)^(−β)` | `D` no separa «decae a un suelo» de «decae a cero», y `G(∞)` es la referencia móvil del §5 de la v3.1 | **NO** |
| 5 | `2ad0700` → (esta) | regresión de `log τ̂₀` contra `log ν` | dispersión de `D` en reloj de ticks contra reloj de pared | `τ̂₀` no está identificada individualmente, así que la regresión no mide nada | **NO** |



| 6 | (esta) | malla de `D` sin puntos en `[0.95, 1.00]` | malla fina con potencia medida | la malla daba la **talla** y no la **potencia**; con potencia baja el paso 4 convertiría la falta de resolución en una victoria de MkII no declarada | **NO** |
| 7 | (esta) | control del reloj de un solo lado | control **espejo** con núcleo fijo en segundos | sin el espejo no se distingue el diagnóstico del artefacto de que `ν` no entra en la definición de `D` a ticks fijos | **NO** |
| 8 | (esta) | `f_∞` estimado siempre | `f_∞` **sólo si se rechazó `D = 1`** | con `β → 0` el núcleo tiende a 1 sea cual sea `f_∞`: no está identificado sin decaimiento | **NO** |

| 9 | (esta, 2026-08-09) | abandono de "estructura explotable a `H*`" | abandono acotado a **escala de SEGUNDOS** | `H*` es el horizonte donde habría que acertar el movimiento entero: un **suelo**, no un punto de operación. El R² requerido cae con `1/H` — 41.8 % en `H*`, 7.2 % a 10 min, 1.2 % a 1 h | **NO** |

**Las nueve enmiendas son anteriores a que exista ningún resultado del §9**, y ninguna se motivó
en un resultado: ocho por controles del propio estimador y la novena por una implicación
aritmética de cifras ya publicadas (`c(u)` con probabilidad de llenado, y la firma de
volatilidad de la v3.1).

⚠ **Sobre `captura_v33`: no se ha ajustado ningún modelo sobre ella, ni se ha formado la
partición de prueba.** Lo único ejecutado es `--resumen`, que es la comprobación de compuerta que
el §1 de este documento ordena hacer y pegar literalmente.

### El estadístico ha cambiado tres veces; la pregunta, ninguna

```
R(final)/R(pico)  →  IC de β̂  →  D = G(K)/G(0)  →  D con suelo f_∞
```

Cada sustitución la forzó un control que invalidaba a la anterior: el primero era degenerado bajo
el nulo, el segundo no estaba identificado, el tercero no separaba suelo de caída a cero. **La
pregunta es la misma desde la v3.1: impacto permanente contra impacto transitorio.** Enunciarla
independiente del estimador es lo que impide que una cuarta sustitución se lea desde fuera como
mover la portería — y se enuncia así en el reporte, con esta tabla al lado.

**Regla para el futuro: a partir del primer contacto con el conjunto de prueba, este documento
queda congelado.** Cualquier cambio posterior invalida el carácter fuera de muestra del resultado
y así se reportaría.

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
