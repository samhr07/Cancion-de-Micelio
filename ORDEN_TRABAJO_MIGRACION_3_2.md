# Orden de Trabajo — Decisión de Migración por Contraste Anidado v3.2

Sucede a `ORDEN_TRABAJO_PROPAGADOR_3_1.md`. Su §2 estableció que **el propagador con
forzamiento medido existe** (`R(τ) > 0` en 799/799 rezagos, factor 150–700 sobre el control con
signos barajados) y dejó **sin decidir su forma** (el estadístico `R(final)/R(pico)` es
degenerado bajo el nulo). El ajuste paramétrico de `G₀`, `τ₀`, `β`, `δ` sigue pendiente.

**Lo que este documento hace:** convierte una elección de arquitectura en una **comparación de
modelos anidados sobre el mismo dato, con el mismo corte y la misma métrica**, decidida fuera de
muestra y con un criterio económico declarado antes de mirar.

**La pregunta que NO se responde aquí:** "¿OFI o delta de ticks?" — está mal planteada. Son
puntos del mismo espacio de modelos, y uno de ellos es el nulo. Ver §1.

**Lo que NO cambia:** `Micelio.py` no se toca. Todo es código de análisis aparte, como en la
v2.2, la v3.0 y la v3.1. La integración en el filtro es la v3.3, y sólo si §9 la habilita.

---

## §0. Registro de errores previos y heredados

| Origen | Error | Estado |
|---|---|---|
| v3.1 §2 (ejecución) | Estadístico `R(final)/R(pico)` **degenerado bajo el nulo** (p5 = −560, mediana = −11 cuando el pico es ruido). El `p = 0.917` no significa nada | **Abierto.** §4.3 lo sustituye por un contraste sobre parámetros del núcleo, no sobre razones de estimadores puntuales |
| v3.1 §2 (ejecución) | `G(0) = 0` interpretado a punto de parada obligatoria. Es artefacto de definir `R(τ) = E[(p_{t+τ} − p_t)·ε_t]`, nulo en `τ = 0` por construcción | Resuelto y documentado. **No repetir el amago:** el test de signo se evalúa sobre `τ ≥ 2`, nunca sobre `τ ∈ {0, 1}` |
| v3.1 §2 (ejecución) | `R(1) = −0.0052 < 0` — es el rebote bid-ask, no signo invertido | Resuelto. ⚠ **Pero implica que todo `R(τ)` medido sobre precio de transacción arrastra el rebote.** Ver §3.4: el banco de pruebas se corre sobre punto medio, no sobre precio de transacción |
| v3.1 §2 (ejecución) | Tramo de 19 157 ticks (0.95 h, ν = 5.6 tx/s) para un procedimiento dimensionado a 48 h | **Abierto.** §2.1 fija el mínimo de datos como compuerta de entrada, no como nota al pie |
| MkII (documento) | `c_t = Y·σ_t/√L_t` introduce **dos constantes libres** (`Y` y la definición operativa de `L_t`) donde `G₀` viene estimado con error estándar | §5.3 lo trata como hipótesis contrastable, no como especificación |
| v2.2 y v3.0 (producción) | Trades con `p = 0` sin filtrar; factor 5 500 en el incremento máximo | Corregido en `mercado.tick_valido`. **Toda captura anterior al 2026-08-08 se descarta**, no se limpia a posteriori |

---

## §1. El objeto de decisión: tres modelos anidados, no tres alternativas

Los candidatos que se están sopesando **no son arquitecturas rivales**. Son restricciones
sucesivas del mismo núcleo:

```
Δp_t = Σ_{k≥0} h(k) · ε_{t−k} · f(v_{t−k}) + η_t,    G(τ) = G₀(1 + τ/τ₀)^(−β),   f(v) = v^δ
```

| id | nombre | restricción | parámetros libres |
|---|---|---|---|
| **M0** | reloj de ticks sin capa predictiva ("delta de ticks") | `G₀ = 0` | 0 |
| **M1** | MkII — OFI con impacto instantáneo y permanente | `β = 0`, `δ = 0.5` | 1 (`G₀`) |
| **M1′** | MkII con `c_t` adaptativo | `β = 0`, `δ = 0.5`, `G₀ → Y·σ_t/√L_t` | 1 (`Y`) + def. de `L_t` |
| **M2** | propagador completo | ninguna | 4 (`G₀`, `τ₀`, `β`, `δ`) |

**Consecuencias inmediatas, y son las que ordenan todo el documento:**

1. **M0 es el nulo, no un candidato.** "Reemplazar `ω_m` por un delta de ticks" deja `A` en
   velocidad constante en reloj de ticks — que es lo que la v2.0 ya implementa y lo que el AR(2)
   de la v3.0 agotó a siete decimales (`φ₁ + φ₂ = 1`, `a = −0.216061`). No es una migración: es
   la supresión de la capa predictiva. Cualquier otro modelo tiene que **batirlo fuera de
   muestra** para existir.
2. **M1 ⊂ M2 estrictamente.** M2 no puede ajustar peor dentro de muestra. La comparación por
   `R²` o log-verosimilitud **in-sample no decide nada** y su uso es motivo de rechazo del
   entregable.
3. **La decisión es fuera de muestra y por parsimonia**, con el margen declarado por adelantado
   (§9.2). Si M2 no bate a M1 por más del margen, **gana M1**, porque tiene tres parámetros
   menos.
4. `δ = 0.5` de MkII ya está en el barrido `δ ∈ {0, 0.25, 0.5, 1}` del §2.3 de la v3.1. **La
   pregunta "¿MkII?" es literalmente una celda de una tabla que aún no se ha llenado.**

⚠ **El sucesor honesto de `ω_m` es `τ₀`**, no un delta de ticks libre: es una escala
característica en tiempo de transacciones, estimada con error estándar en vez de elegida. Si
alguien escribe `Δn` como constante en el código, se ha reabierto el agujero que el §7 de la
v3.1 cerraba.

---

## §2. El banco de pruebas común — condición previa a cualquier ajuste

**Ningún resultado de §4–§7 es admisible si no sale de este banco.** El fallo clásico de esta
comparación es evaluar M1 sobre una muestra alineada a `bookTicker` y M2 sobre la de `@trade`, y
atribuir a los modelos una diferencia que es de muestreo.

### 2.1 Compuerta de datos

| requisito | valor | por qué |
|---|---|---|
| Tramo continuo limpio | **≥ 8 h**, objetivo 48 h | 12 bloques de bootstrap contra los ~15 mínimos en la v3.1 |
| Sin cortes de DNS internos | verificado por diferencias de timestamp | la captura larga de la v2.2 se partió por un corte de 36 442 s |
| `tr_maker` y `q` persistidos | obligatorio | sin ellos no hay `ε` ni `f(v)` |
| `@bookTicker` simultáneo | obligatorio para M1/M1′ | `captura_dual.py` ya lo hace |
| Filtro `tick_valido` aplicado | obligatorio | §0 |

Si el tramo no llega a 8 h, **no se ejecuta el §5 ni el §7** y el reporte lo dice en la primera
línea. Relanzar `captura_larga.py`.

### 2.2 Reloj, corte y embargo

- **Reloj primario: transacciones.** Todo `τ` en ticks. Todo resultado se reporta **además** en
  segundos con la `ν` del tramo, y §6 contrasta en cuál de los dos vive el decaimiento.
- **Corte temporal 60 / 20 / 20** — entrenamiento, validación (para elegir `δ` y la malla de
  §7), prueba (se toca **una sola vez**, al final).
- **Banda de embargo ≥ `H*`** entre bloques contiguos, en ticks: `H*·ν`. Con `H* = 50 s` y
  `ν = 39`, son ~1 950 ticks. Verificarlo por test, no por comentario.

### 2.3 Métrica única

**Log-verosimilitud fuera de muestra por observación**, `LL/N`, sobre el conjunto de prueba, con
intervalo por **bootstrap por bloques** de longitud `≥ 5·H*·ν`. `R²` se reporta como
descriptivo; no decide.

---

## §3. Alineación de `bookTicker` con `@trade` — donde vive el fallo más peligroso

M1 necesita `e_t` (OFI de nivel 1). M2 sólo necesita `ε` y `q`, que vienen en el mismo mensaje
que el precio. **Esa asimetría es una fuente de fuga de información unilateral**, y favorece a
M1 de forma invisible.

### 3.1 La regla

Para cada transacción en `t`, el estado del libro utilizable es el del **último `bookTicker`
con timestamp estrictamente anterior** a `t`, y `e_t` se construye con ése y el inmediatamente
previo. Nunca interpolar. Nunca `ffill` hacia atrás.

### 3.2 Test de fuga obligatorio

Ejecutar el ajuste de M1 **tres veces**:

| variante | `LL/N` esperado |
|---|---|
| alineación correcta (§3.1) | el valor a reportar |
| libro adelantado 1 actualización (fuga deliberada) | **debe subir claramente** |
| libro retrasado 10 actualizaciones | debe bajar |

Si la variante con fuga deliberada **no** mejora, el pipeline no está usando el libro en
absoluto y `e_t` es ruido. Si la correcta y la adelantada empatan, hay fuga en la correcta.
**Las tres cifras van en el reporte**, no sólo la del medio.

### 3.3 ⚠ `bookTicker` no da OFI de verdad — declararlo, y cuantificarlo

`@bookTicker` es un *snapshot* con estrangulamiento: las cancelaciones entre actualizaciones son
invisibles, y son precisamente las que la fórmula de `e_t` pretende penalizar. **Lo que se
contrasta es "OFI-L1 aproximado desde bookTicker", y así debe etiquetarse en todo el reporte.**
El OFI real exige el flujo diferencial `@depth`, que hoy no se ingiere.

Diagnóstico exigido, para saber cuán mala es la aproximación:

- Fracción de transacciones cuyo `q` **excede** la cantidad en el mejor nivel del último snapshot
  observado. Cada una es prueba de al menos una actualización perdida.
- Reconciliación: la parte de `ΔV_b`, `ΔV_a` atribuible a transacciones observadas frente a la
  variación total. El residuo es actividad de límite no observada.

Reportar ambas como porcentaje. **Si la primera supera el 20 %, M1 se evalúa igualmente pero su
resultado se marca como cota inferior**, no como medición.

### 3.4 ⚠ Todo sobre punto medio, no sobre precio de transacción

`R(1) = −0.0052` es rebote bid-ask. Un modelo que "predice" que la siguiente transacción cae al
otro lado del spread predice algo **real y no operable**: no se puede cruzar el spread y cobrar
el spread. Con `s_eff = 0.2062 USD/BTC` ≈ 2 ticks, esa cantidad es del orden de todo lo que hay
en juego a rezagos cortos.

**El banco de pruebas se corre sobre `p_mid = (bid + ask)/2`.** El precio de transacción se
mantiene como control secundario, y **si M1 gana sobre transacción y pierde sobre punto medio,
lo que se ha medido es el rebote y la conclusión es que M1 no aporta.** Reportar las dos
columnas siempre.

---

## §4. Test A — el ajuste del núcleo decide M1 contra M2

### 4.1 Estimación

Ajustar `G(τ) = G₀(1 + τ/τ₀)^(−β)` y `f(v) = v^δ` sobre el conjunto de entrenamiento, con
`δ ∈ {0, 0.25, 0.5, 1}` elegido en **validación** (no en prueba, §2.2). Errores estándar por
bootstrap por bloques.

Normalizar `v` por su mediana del tramo antes de elevar a `δ`: si no, `G₀` y `δ` quedan
confundidos por la escala y sus errores estándar no son interpretables.

### 4.2 La tabla que decide

| modelo | `LL/N` prueba | Δ contra M0 | IC bootstrap 95 % | nº parám. |
|---|---|---|---|---|
| M0 | — | 0 (referencia) | — | 0 |
| M1 (`β = 0`, `δ = 0.5`) | | | | 1 |
| M1′ (`c_t` adaptativo) | | | | 1 + def. |
| M2 (`δ` en validación) | | | | 4 |

Y, por separado, la que contesta la pregunta directamente:

| contraste | estadístico | regla |
|---|---|---|
| `β = 0`? | `β̂` con IC bootstrap | si el IC **contiene 0**, MkII está bien especificada y se prefiere por parsimonia |
| `δ = 0.5`? | `LL/N` de validación por celda | si `δ = 0.5` es el máximo o está dentro de su IC, MkII acierta también aquí |

**Esto es la respuesta a la pregunta de la migración, y sale de un ajuste, no de una preferencia.**

### 4.3 Sustituto del estadístico degenerado

El sobrepaso (raíces complejas) **no se contrasta con `R(final)/R(pico)`**. Se contrasta
ajustando un núcleo que **admita** oscilación y viendo si la admite:

```
G(τ) = G₀ · (1 + τ/τ₀)^(−β) · cos(ω_G·τ + φ)      (M2-osc, 5 parámetros)
```

Contraste `ω_G = 0` por razón de verosimilitud **con nulo simulado**: generar bajo M2 ajustado
(`ω_G = 0`), reajustar M2-osc, y construir la distribución de `2ΔLL`. Comparar el valor real
contra esa distribución.

⚠ **`ω_G` no es `ω_m` resucitada, y hay que decirlo en el reporte para que nadie lo lea así.**
`ω_m` era una frecuencia del *precio*, extraída por EMD, sin nulo. `ω_G` es una frecuencia de la
*respuesta al impulso del flujo*, con estimador propio, error estándar y nulo simulado. Si sale
`ω_G ≠ 0`, el sistema es un oscilador **forzado** y la premisa de la v3.0 sobrevive en la única
forma en que podía sobrevivir. Si sale `ω_G = 0`, el propagador es monótono y `ω_m,max` y `γ_ω`
se borran de `constantes_micelio.py` (§7 de la v3.1).

### 4.4 Residuo

Ljung-Box sobre el residuo de cada modelo, a rezagos hasta `2·H*·ν`. **Un modelo con residuo
estructurado no gana aunque tenga mejor `LL/N`**: significa que su ventaja es de ajuste local y
no de especificación.

---

## §5. Test B — ¿aporta el lado de límite algo sobre el flujo de transacciones?

Ésta es la pregunta empírica que separa realmente a MkII del propagador, y no es "OFI contra
signo": **el signo de transacción `ε` ya es el lado taker del OFI**, y ya se midió con factor
150–700 sobre su control. Lo que MkII añade es la actividad de límite.

### 5.1 Descomposición

```
e_t = e_t^(taker) + e_t^(limite)
```

con `e_t^(taker)` reconstruible desde `@trade` y `e_t^(limite)` el residuo de §3.3. Modelo
anidado:

| id | regresores |
|---|---|
| M2 | `ε·f(v)` con núcleo `G(τ)` |
| **M2+L** | M2 más `e_t^(limite)` con su propio núcleo `G_L(τ)` |

Si `LL/N(M2+L) − LL/N(M2)` no excede el margen del §9.2, **el lado de límite no aporta**, y toda
la maquinaria de `bookTicker`, `@depth` y anti-*spoofing* de MkII es coste sin retorno. Ése es un
resultado perfectamente publicable y ahorra una dependencia entera.

### 5.2 ⚠ La raíz cuadrada sesga bajo ruido de medición

MkII usa `sgn(e)·√|e|`. Si `e` se mide con error —y por §3.3 se mide con error—, entonces
`E[√|e + ξ|] ≠ √|E[e]|`: la concavidad de `√` sesga sistemáticamente hacia **sobreestimar** el
impacto de los `e` pequeños. Con `e ≈ 0` el sesgo es máximo, y `e ≈ 0` es la mayoría de las
observaciones.

**Guarda:** simular con `e` real más ruido de magnitud calibrada por §3.3, ajustar, y reportar
el sesgo en `G₀`. Si el sesgo es del orden de `G₀`, el estimador no vale y hay que usar
`f(v) = v^δ` con `δ` libre en vez de fijar `0.5` — que es lo que M2 hace de todas formas.

### 5.3 M1′ y las constantes libres de `c_t`

`c_t = Y·σ_t/√L_t` se contrasta como **hipótesis**, no se adopta como especificación:

| variante | pregunta |
|---|---|
| `G₀` constante | referencia |
| `G₀ → Y·σ_t` | ¿la ganancia viene sólo del escalado por volatilidad? |
| `G₀ → Y/√L_t` | ¿aporta la profundidad? |
| `G₀ → Y·σ_t/√L_t` | ¿aportan las dos? |

⚠ **`L_t` desde `bookTicker` es la cantidad en el mejor nivel, no la profundidad.** Llamarlo
profundidad en el reporte sería la clase de deslizamiento de nombres que ya produjo la colisión
`C_max`/`I_max`. Etiquetarlo `V_best`, y decir explícitamente que la profundidad real necesita
`@depth`.

Si la variante `Y·σ_t` sola captura toda la ganancia, **M1′ no está midiendo microestructura:
está reescalando por volatilidad**, y eso ya lo hace `ρ_k` en el filtro.

---

## §6. Test C — ¿en qué reloj vive el propagador?

Pregunta abierta declarada en `micelio_hallazgos.tex`, y aquí es decisoria porque `τ₀` es el
número que sustituye a `ω_m`.

Las capturas disponibles tienen `ν` muy distintas (5.6, 39 y ~102 tx/s medidas en distintas
sesiones). Ajustar `τ₀` por bloques y contrastar:

| hipótesis | predicción |
|---|---|
| el impacto decae en tiempo de transacciones | `τ₀` [ticks] estable, `τ₀/ν` [s] varía con `ν` |
| el impacto decae en tiempo de pared | `τ₀/ν` [s] estable, `τ₀` [ticks] varía con `ν` |

Regresión de `log τ₀` contra `log ν` por bloques: pendiente 0 favorece la primera; pendiente 1,
la segunda. **Reportar la pendiente con su IC, no el veredicto.**

⚠ **Confundidor:** `ν` correlaciona con volatilidad y con sesión. Incluir el regresor de régimen
del §6 de la v3.1 (volatilidad realizada, congelado por adelantado) como control, o el resultado
mide otra cosa.

---

## §7. Embebido markoviano — sólo si M2 gana

Un núcleo de ley de potencias **no es markoviano** y no cabe exacto en estado finito. Este
apartado especifica la aproximación, para que la v3.3 no la improvise.

### 7.1 La forma

```
G(τ) ≈ Σ_{j=1..J} w_j · exp(−τ/T_j),      J ∈ {2, 3}
```

Estado aumentado del EAKF: `x = [S, v, R_n, I₁, …, I_J]`, con

```
I_j(n+1) = exp(−1/T_j)·I_j(n) + w_j·ε_n·f(v_n)
```

`ε_n f(v_n)` es **entrada exógena medida**, no estado y no ruido. Ésa es toda la diferencia con
la v3.0, donde el forzamiento vivía en el residuo.

### 7.2 ⚠ Prony es mal condicionado — la guarda es estructural, no numérica

Ajustar simultáneamente `w_j` y `T_j` a una ley de potencias es un problema notoriamente mal
condicionado: infinitas combinaciones dan la misma curva a tolerancia. Un ajuste que converja
sin quejarse **no es evidencia de identificabilidad**.

**Guarda:** fijar las tasas en una malla geométrica declarada por adelantado, p. ej.
`T_j = τ₀·10^(j−1)` para `j = 1..J`, y **ajustar sólo los pesos `w_j`** — que entonces es
mínimos cuadrados lineal, bien condicionado y con errores estándar honestos.

**Y no reportar nunca `w_j` ni `T_j` individuales como magnitudes con significado.** Sólo:

- el error de reconstrucción `max_τ |G_ajustado(τ) − G(τ)|/G(τ)` sobre `τ ∈ [1, 2·H*·ν]`;
- el `LL/N` de la versión embebida contra el de M2 con núcleo exacto.

Si el segundo cae más que el margen del §9.2, `J` es insuficiente. Subir a `J = 3` antes de
concluir nada sobre la aproximación.

### 7.3 ⚠ Doble contabilidad del impacto propio

Si `ε` incluye **nuestras propias** transacciones y la capa de Loeper ya modela nuestro impacto,
el sistema cuenta su impacto dos veces: sobreestima el movimiento esperado y opera de más.

**Guarda:** separar `ε` en propio y ajeno. En demo tenemos nuestro `orderId`, así que es
identificable. El propagador se estima **sólo con flujo ajeno**; el propio entra por Loeper,
donde ya está. Test permanente: inyectar transacciones propias sintéticas y verificar que
`G(τ)` estimada no se mueve.

### 7.4 ⚠ `Q` absorbe la entrada si no se recalibra

Añadir `B·u` sin recalibrar `Q` cambia el NIS por una razón que no es la calidad del modelo. La
comparación entre EAKF con y sin propagador **exige recalibrar `Q` en ambos con el mismo
procedimiento**, o no compara modelos: compara calibraciones.

Recordatorio de la v1.2: la inflación de `Q` llevó el NIS de 19.35 a 0.53 dejando `ρ₁` en +0.44.
**El NIS sin `ρ₁` al lado no es evidencia de nada.**

---

## §8. Modos de fallo previstos

Tabla de vigilancia. Cada fila lleva su detección; una fila sin detección implementada es una fila
que no se está vigilando.

| # | fallo | firma | detección |
|---|---|---|---|
| 1 | Comparar M1 y M2 sobre muestras distintas | M1 gana por márgenes grandes y estables | §2: `assert` de identidad de índices entre conjuntos de evaluación |
| 2 | Fuga por alineación de libro | `R²` fuera de muestra sospechosamente alto | §3.2, las tres variantes |
| 3 | La ganancia es el rebote bid-ask | M1 gana sobre precio de transacción, pierde sobre punto medio | §3.4, dos columnas siempre |
| 4 | `δ` elegida en el conjunto de prueba | contraste múltiple encubierto sobre 4 celdas | §2.2: `δ` se fija en **validación**; prueba se abre una vez |
| 5 | Sesgo de Jensen por `√\|e\|` | `G₀` sobreestimada, peor cuanto más ruidoso `e` | §5.2, simulación con ruido calibrado |
| 6 | Prony no identificable | ajuste converge, parámetros absurdos, curva correcta | §7.2: malla fija, sólo pesos libres |
| 7 | Doble contabilidad del impacto propio | el sistema opera más al aumentar su propio tamaño | §7.3, separación por `orderId` |
| 8 | `Q` absorbiendo la entrada | NIS mejora, `ρ₁` intacta | §7.4, recalibración simétrica + `ρ₁` obligatoria |
| 9 | `LL` in-sample usada para decidir | M2 gana siempre (es superconjunto) | §1.2, rechazo del entregable |
| 10 | `ν` confundida con volatilidad en §6 | pendiente de `log τ₀` vs `log ν` significativa y espuria | §6, control de régimen congelado |
| 11 | Bootstrap con bloques cortos | IC demasiado estrechos, todo "significativo" | longitud ≥ `5·H*·ν`, y nº de bloques reportado |
| 12 | Signo de `ε` invertido | `R(τ) < 0` en **todo** el rango, no sólo en `τ = 1` | §4 de la v3.1, ya verificado: 799/799 positivos |
| 13 | Ganar a M0 sin ganar a los costes | `LL/N` mejor, `E[Δp] < c(u)` | §9.1, criterio 3 |
| 14 | Tramo insuficiente y reporte igualmente concluyente | 12 bloques presentados como 15 | §2.1, compuerta antes de ejecutar |
| 15 | `Δn` reintroducido como literal | constante mágica en el código | revisión: `grep` de literales en la ruta del filtro |

⚠ **Patrón, ya en su tercera aparición:** umbral asintótico usado donde la distribución finita es
otra (χ² sobre NIS en v2.1; `(1−γ)/2` en v3.1 §3.2; y ahora `2ΔLL ~ χ²` en §4.3, que con
observaciones dependientes **no se cumple**). Todo contraste de razón de verosimilitud de este
documento va contra **nulo simulado**, nunca contra la tabla χ².

---

## §9. Preregistro — criterios fijados ANTES de mirar

Commitear `PREREGISTRO_3_2.md` con fecha antes de ejecutar nada; referenciar el hash en el
reporte. Es el tercer documento de este proyecto que lo exige y las dos veces anteriores salvó
un veredicto.

### 9.1 Regla de decisión

Se aplica **en este orden**, y se para en el primer fallo:

1. **Compuerta de datos** (§2.1). Si no pasa, no hay decisión: hay que capturar más.
2. **`LL/N`(mejor modelo) > `LL/N`(M0)** fuera de muestra, con IC bootstrap que **excluya 0**.
   Si falla → **gana M0**: la capa predictiva se retira y el proyecto pasa a la lista de
   supervivientes del §6.2 de la v3.0 (reloj de ticks, ingesta, riesgo, Loeper, NMPC como
   ejecutor).
3. **`E[Δp | ℱ_t]` al horizonte `H*` supera `c(u)` con margen** en al menos un régimen
   congelado. Si falla → M0 igualmente, por §4.2 de la v3.1. **Esta es la única comparación con
   contenido económico y ningún `LL/N` la sustituye.**
4. **M2 contra M1 por §9.2.** Si M2 no supera el margen → **gana M1** (MkII), por parsimonia.
5. **M2+L contra M2** por el mismo margen. Decide si hace falta `@depth`.

### 9.2 El margen, declarado ahora

`ΔLL/N` mínimo para preferir el modelo más complejo: **el límite inferior del IC bootstrap al
95 % debe superar `k/(2N)`**, con `k` la diferencia de parámetros y `N` el tamaño de la muestra
de prueba — es el criterio BIC expresado por observación, con la incertidumbre del estimador
incorporada en vez de ignorada.

### 9.3 Criterio de abandono

Sin cambios respecto al §4.2 de la v3.1, que sigue vigente y **este documento no relaja**. Si el
paso 2 o el 3 del §9.1 fallan sobre datos que pasan la compuerta, la hipótesis de estructura
explotable a `H*` queda abandonada. No hay otra escala a la que retirarse: ésta es donde está el
dinero.

---

## §10. Criterios de aceptación

**§2 — Banco de pruebas**
- [ ] Tramo continuo limpio ≥ 8 h, sin cortes, con `tr_maker`, `q` y `bookTicker`.
- [ ] Identidad de índices entre conjuntos de evaluación verificada por `assert`.
- [ ] Embargo ≥ `H*·ν` verificado por test, no por comentario.
- [ ] Nº de bloques de bootstrap reportado explícitamente.

**§3 — Alineación**
- [ ] Las **tres** variantes del test de fuga reportadas, con la deliberada mejorando.
- [ ] Porcentaje de transacciones con `q` > cantidad en mejor nivel del último snapshot.
- [ ] Reconciliación de `ΔV_b`, `ΔV_a` con transacciones observadas.
- [ ] Todos los resultados en **dos columnas**: punto medio y precio de transacción.

**§4 — Test A**
- [ ] Tabla de `LL/N` fuera de muestra para M0, M1, M1′, M2 con IC bootstrap.
- [ ] `β̂` con IC: ¿contiene 0?
- [ ] `δ` elegida en validación; celda `δ = 0.5` reportada aunque no gane.
- [ ] `ω_G` contra **nulo simulado**, no contra χ².
- [ ] Ljung-Box del residuo a rezagos hasta `2·H*·ν`, por modelo.

**§5 — Test B**
- [ ] `LL/N`(M2+L) − `LL/N`(M2) con IC.
- [ ] Sesgo de Jensen sobre `G₀` cuantificado por simulación.
- [ ] Las cuatro variantes de `c_t` del §5.3, con `L_t` etiquetado `V_best`.

**§6 — Test C**
- [ ] Pendiente de `log τ₀` contra `log ν` con IC, controlada por régimen.
- [ ] `τ₀` reportada en ticks **y** en segundos.

**§7 — Embebido (sólo si M2 gana)**
- [ ] Malla `T_j` declarada antes del ajuste; sólo pesos libres.
- [ ] Error de reconstrucción de `G(τ)` reportado; `w_j`, `T_j` individuales **no** interpretados.
- [ ] Separación de `ε` propio/ajeno con test de inyección sintética.
- [ ] `Q` recalibrada con el mismo procedimiento en ambos modelos; NIS **y** `ρ₁` juntos.

**§9 — Preregistro**
- [ ] `PREREGISTRO_3_2.md` commiteado con fecha anterior a la primera ejecución; hash citado.
- [ ] Regla de decisión aplicada en orden, con el paso en que se detuvo indicado.

**Transversal**
- [ ] Suite en verde: 56/56 + los nuevos.
- [ ] `Micelio.py` **sin cambios**.
- [ ] Ningún literal nuevo en la ruta del filtro (`grep` incluido en el reporte).
- [ ] Todo contraste de verosimilitud contra nulo simulado.

---

## §11. Fuera de alcance

- **Implementar el ganador en el filtro.** Este documento decide **cuál**; la integración es la
  v3.3 y usa el §7 como especificación.
- **Ingestión de `@depth`.** Sólo si el §5.1 demuestra que el lado de límite aporta. Decidirlo
  antes sería añadir una dependencia por si acaso.
- **Retirar la cadena EMD/HHT del código.** Depende de `ω_G` (§4.3). Hasta entonces queda, sin
  alimentar decisiones.
- **`ω_m,max`, `γ_ω`, `C_ON`, `C_OFF`, `W`, `K`**: su borrado está condicionado a §4.3, no a este
  documento.
- **v2.1 §5 (log-precio) y §6 (cuantiles empíricos del NIS)**: siguen en cola y **no
  bloqueadas**. Con curtosis 1179.7, §6 mantiene la prioridad que la v3.1 le dio.
- Fase 2 (ALS), Fase 3, Testnet con credenciales, las 30 corridas.
- CMA-ES, IMPC, DeepONet.
- Reconciliación documental del PDF.
