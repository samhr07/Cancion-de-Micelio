# Orden de Trabajo — Las Dos Curvas: dónde (si es que hay dónde) la señal paga el peaje · v4.1

Sucede a `ORDEN_TRABAJO_TICK_GRANDE_3_3.md`, ejecutada el 2026-08-10 (b).

**Por qué v4.1 y no v3.4.** La etiqueta `v3.4` está reservada desde el 2026-08-09 para el cambio
de coste lineal en el NMPC, que toca `Micelio.py`. Este documento **no toca `Micelio.py`** y es la
continuación natural de la v4.0 (costes de ejecución): junta la curva de coste de la v4.0 con la
curva de señal de la v3.2. No se reusa el número reservado.

---

## §0. Registro de errores previos y heredados

| Origen | Error | Estado |
|---|---|---|
| v3.3 §1 (documento) | `γ̂ = 0.798` etiquetada «autocorrelación de signos» era `C(1)`, no el exponente de decaimiento | **Corregido por Code**, 2026-08-10 (b). El `γ` correcto es 0.5217. La conclusión «impacto permanente» venía del `γ` equivocado |
| v3.3 §1 + §6 (ejecución) | **`γ` se midió DOS VECES en la misma sesión con dos métodos y salieron 0.5217 y 0.3539.** La discrepancia (0.168) es **mayor que todo el intervalo de hipótesis** (`β` de permanencia 0 a `β` de difusividad 0.239) | **ABIERTO.** §3.1 |
| v3.3 §6 (ejecución) | La columna `N_eff` de la tabla usa `N_eff = N^γ` (sin la constante `c_var`) y la columna `inflación` usa `√(Var_LRD/Var_iid)`. **Son incoherentes entre sí**: ninguna fila cumple `inflación = √(N/N_eff)` | **ABIERTO.** §3.2 |
| v3.3 §1 (ejecución) | `γ` ajustada en rezagos `[10, 2000]`, `H` en `n ∈ [256, 16384]`. La relación `H = (2−γ)/2 − β` es asintótica y exige **el mismo régimen de escala**. Solapamiento real: `[256, 2000]`, menos de una década | **ABIERTO.** §3.3 |
| v3.2 (ejecución) | `β = −0.160` sobre punto medio, `D = 11.53`, núcleo **creciente** | Tratado hasta ahora como medición. `text_3.txt` afirma que `β < 0` **nunca se ha reportado** y es firma de mala especificación, no del mercado. §3.4 lo reclasifica |
| v3.2 (ejecución, limitación 1) | Verosimilitud gaussiana mal especificada: `y_mid` es **96.67 % ceros exactos** | **ABIERTO.** §6 propone el micro-precio como observable que no tiene ese problema |
| Traspaso §6 contra Preregistro §11.1 | El `R²` requerido a 10 min figura como **2.5 %** en un documento y **7.2 %** en el otro | **Explicado**: 2.5 % usa `H_p = 0.591`, 7.2 % usa `H_p = 0.5`. §2 lo cierra midiendo `H_p` |

---

## §1. La única hipótesis viva: el horizonte. Y cómo se contrasta

### 1.1 Por qué es la única

La v3.2 midió, fuera de muestra, sobre 203 998 observaciones de prueba:

```
q90(|mu|) = 11.30 USD/BTC = 1.736 pb
max(|mu|) = 17.12 USD/BTC = 2.630 pb        <- el MAXIMO, no el decil
comision ida y vuelta maker+maker = 4.000 pb
```

**Ni una sola observación del conjunto de prueba llega a pagar la comisión bruta**, y falta un
factor 1.52 incluso retirando por completo el margen de seguridad de 1.5. Esto no es una brecha
que un estimador mejor pueda cerrar: es categórico al horizonte de 2 233 ticks (~128 s).

Las palancas disponibles, y su tamaño real:

| palanca | techo alcanzable | ¿cierra 3.46×? |
|---|---|---|
| Descuento BNB sobre comisión | 4.0 → 3.6 pb (−10 %) | **No** |
| Nivel VIP 1 (exige volumen 30 d) | maker 2.0 → 1.6 pb por lado (−20 %) | **No**, y fuera de alcance con capital restringido |
| Reformulación en tick grande | el propio §3 de la v3.3 lo dijo antes de correrlo: **0 pb** | **No** |
| Micro-precio / mejor observable | corrige sesgo de atenuación; optimista, factor ~1.5 | **No** |
| **Horizonte** | comisión **fija en pb**, movimiento **crece con `H`** | **Es la única con el orden de magnitud correcto** |

⚠ **Se declara antes de medir:** el resultado «las dos curvas no se cruzan en ninguna banda» es
un desenlace perfectamente admisible y se reporta como tal. Este documento no busca rescatar el
proyecto; busca decidirlo.

### 1.2 Las dos curvas

Todo se reduce a dos funciones de `H` sobre los mismos ejes.

**Curva requerida** (derivada, no medida). Operando sólo el decil superior de una señal, con
`E[|μ̂| \| decil superior] ≈ 1.755·σ_μ`:

```
R²_req(H) = ( lastre(H) / (1.755 · σ(H)) )²

sigma(H)  = sigma_1 · H^H_p              sigma_1 = 3.77 USD/BTC·s^(1/2)
lastre(H) = c(u) + seleccion_adversa + financiacion(H)
```

Esta forma **reproduce exactamente** la tabla del §11.1 del `PREREGISTRO_3_2.md` (41.8 % a 103 s,
7.2 % a 10 min, 1.2 % a 1 h, 0.33 % a 4 h), así que está verificada contra lo ya publicado.

**Curva medida** (esto es lo que no existe). `R²_medido(H)` fuera de muestra, del mejor predictor
disponible, para `H` en la banda de minutos a horas.

**Se conoce un solo punto de la curva medida: `R² = 0.014` a `H ≈ 128 s`.** Todo lo demás del
argumento del horizonte es extrapolación con la curva medida supuesta **plana**, que es
precisamente lo que no puede ser: la capacidad predictiva del flujo de órdenes decae con el
horizonte.

Si se supone plana —el caso más favorable posible— el cruce está en:

| `H_p` | cruce con `R²_medido = 1.4 %` |
|---|---|
| 0.500 (difusivo) | **51.3 min** |
| 0.591 (medido en `captura_v33`) | **14.9 min** |

**Un factor 3.4 en el horizonte de operación colgando de un exponente que no ha replicado.** De
ahí que el §2 vaya primero.

### 1.3 Procedimiento

1. **Predictor deliberadamente simple, sin propagador.** Regresión de `r_{t→t+H}` (log-retorno del
   punto medio) sobre flujo firmado acumulado en ventanas retrospectivas
   `{H/4, H/2, H, 2H}`. Sin `G(τ)`, sin `τ₀`, sin `β`. **Motivo declarado:** la degeneración
   `β`–`τ₀` mató cuatro estadísticos y no debe contaminar la medición del horizonte. Refinar el
   predictor sólo tiene sentido si la curva cruza con el predictor tonto o queda cerca.
2. **Ventanas NO solapadas.** El `N` efectivo de esta medición es `T_captura / H`, y hay que
   reportarlo por cada `H`. Con `captura_v33` (16.38 h): 196 ventanas a 5 min, 98 a 10 min,
   16 a 1 h. **A partir de 30 min la medición es un piloto, no una decisión, y se etiqueta así.**
3. **Partición.** Entrenar sobre el 60 % de entrenamiento de la v3.2, medir `R²` sobre validación.
   **El conjunto de prueba de la v3.2 NO se abre.** Ya se abrió dos veces.
4. **Control negativo obligatorio:** el mismo procedimiento sobre flujo con signos barajados
   dentro de cada ventana. `R²_barajado` va en la misma tabla, siempre. Sin él, un `R²` de 1 %
   sobre 16 ventanas no significa nada.
5. **Reportar la tabla completa** de `H` ∈ {1, 2, 5, 10, 15, 30, 60} min con
   `R²_medido`, `R²_barajado`, `R²_req(H_p = 0.5)`, `R²_req(H_p = medido)`, `n_ventanas`.

### 1.4 Regla de lectura, declarada antes de mirar

| resultado | lectura |
|---|---|
| `R²_medido` supera `R²_req` en alguna `H` con `n_ventanas ≥ 30` y control barajado limpio | **hay banda viable**; se preregistra el contraste económico completo sobre esa banda |
| `R²_medido` decae más rápido que `R²_req` en todo el rango | **no hay banda viable con este predictor**; el paso 3 falla en todo horizonte medible |
| Cruce sólo con `n_ventanas < 30` | **no decidible**; espera a `captura_estacional`, no se decide con 16 ventanas |
| `R²_barajado` comparable al real | la medición no tiene potencia a esa `H`; se reporta y no se lee |

---

## §2. `H_p`: el exponente del que cuelga la curva requerida. Va PRIMERO

`H_p = 0.591` sobre `captura_v33` (ν = 17.5 tx/s) y `≈ 0.5` sobre `captura_larga` (ν = 39 tx/s).
**No ha replicado.** Y de él depende un factor 3.4 en el horizonte de cruce.

### 2.1 Lo que hay que hacer

1. **Firma de volatilidad sobre las cuatro capturas en disco** (`captura_larga`, `captura_v31b`,
   `captura_v32`, `captura_v33`), **con su control barajado al lado en todas**, que es lo que la
   sesión 2026-08-08 (e) dejó como obligatorio.
2. **En los dos relojes, y declarado cuál es cuál.** El hallazgo de subdifusión de la v2.2 fue un
   artefacto de mezclar relojes en el mismo estimador. La curva requerida vive en **segundos de
   pared** (la comisión y la financiación se pagan en tiempo de calendario), así que el `H_p` que
   entra en `σ(H)` debe ser el de **reloj de pared**, no el de ticks. El de ticks se reporta al
   lado, como control.
3. **`σ₁` re-medido por captura.** `σ₁ = 3.77 USD/BTC·s^(1/2)` viene del plateau de la v3.1 y
   **está en USD/BTC a un precio que ya envejeció**. Reexpresar en **pb·s^(−H_p)** y recalcular.
   Es el mismo error que el §3 de la v3.3 cazó en `c(u)`.
4. **Rango de validez explícito.** El alcance verificado es 711 s. Toda fila de 1 h en adelante es
   extrapolación y se marca como tal **en la propia tabla**, no en una nota al pie.

### 2.2 Compuerta

| resultado | consecuencia |
|---|---|
| `H_p` compatible entre las cuatro capturas | se usa el valor común; la curva requerida queda fijada |
| `H_p` varía y correlaciona con `ν` | **la superdifusión es un efecto de la tasa de actividad, no del mercado**; la curva requerida se construye por régimen de `ν` y el §1 se corre por régimen |
| `H_p` varía sin patrón | la curva requerida se reporta **como banda** entre `H_p = 0.5` y el máximo medido, y el §1 se lee contra la banda entera |

⚠ **Nunca se elige el `H_p` que hace cruzar las curvas.** Si el §1 cruza sólo con el `H_p` más
favorable, el resultado es «no decidible», no «hay banda viable».

---

## §3. Correcciones al reporte de la v3.3

Ninguna cambia un veredicto. Todas cambian un número que otros documentos ya están citando.

### 3.1 `γ` medida dos veces, dos respuestas

| método | `γ` | dónde |
|---|---|---|
| Regresión log-log de `C(ℓ)`, rezagos `[10, 2000]` | **0.5217** | v3.3 §1 |
| Escalado de `Var(media de bloque)` contra `b` | **0.3539** | v3.3 §6, fila «signos `ε`» |

Sobre la **misma serie de signos del mismo tramo**. Propagado a `β` implícita:

| `γ` usada | `β` implícita | `β` de difusividad | posición relativa |
|---|---|---|---|
| 0.5217 | +0.1484 | +0.2391 | 62 % |
| 0.3539 | **+0.2324** | +0.3231 | 72 % |

**El `β` implícito casi se dobla.** La posición relativa (62 % contra 72 % del camino hacia
difusividad) es más estable que el valor absoluto, y ésa es la única lectura que sobrevive a la
discrepancia.

**Qué hacer:** un tercer estimador que no comparta modo de fallo con ninguno de los dos —
**Whittle local o regresión log-periodograma (GPH)** — sobre la misma serie. Y declarar por qué
se esperan sesgos: el ACF muestral de un proceso de memoria larga está **sesgado a la baja** por
centrar con la media muestral, lo que empina el ajuste log-log y **sobreestima `γ`**. Eso apunta
a que 0.5217 es el sesgado y 0.3539 el más cercano — pero es una predicción, y se verifica, no se
asume.

⚠ **Y una identidad estructural que hay que escribir en el reporte, porque cambia cómo se lee
todo el §1 de la v3.3:**

```
beta_implicita = (2 - gamma)/2 - H = 0   <=>   gamma = 2 - 2H
```

que es exactamente la relación entre exponente de ACF y exponente de Hurst del **ruido gaussiano
fraccionario**. O sea: `β` implícita **no mide el núcleo**; mide **cuánto se desvía el par
(signos, precio) de la relación fGn**. Si signos y precio comparten el mismo proceso de memoria
larga, `β` sale 0 **mecánicamente**, sin que eso diga nada sobre permanencia. Con `corr(γ̂, Ĥ) =
+0.18` entre las 8 sub-series y `Cov` aportando el 14 % de `Var(β)`, los dos estimadores no son
independientes. **`β` implícita es un contraste de coherencia entre dos exponentes, no una tercera
medición del núcleo**, y presentarla como «tercera ruta al mismo parámetro» —que es como la
escribió la v3.3— sobrevende lo que es.

### 3.2 La tabla de `N_eff` es internamente incoherente

`tick_grande.ic_memoria_larga` devuelve dos cantidades que deberían satisfacer
`inflación = √(N/N_eff)` y no la satisfacen en ninguna fila:

| serie | `N` | inflación reportada | `N_eff` reportado (`= N^γ`) | **`N_eff` coherente (`= N/inflación²`)** |
|---|---|---|---|---|
| `ΔLL` de prueba | 199 532 | 3.15 | 8 574 | **20 109** |
| signos `ε` | 618 692 | 115.4 | 112 | **46.5** |
| incrementos de precio | 618 691 | 0.71 | 99 510 | **1 227 318** |

**La causa:** `n_eff = n**gamma_eff` descarta la constante `c_var` del escalado, así que no es un
número de observaciones independientes — es `N^γ` a secas. La columna con contenido es la de
inflación, que sí usa el intercepto ajustado.

Consecuencias que hay que propagar:

- **Signos: `N_eff ≈ 47`, no 112.** El diagnóstico de la v3.3 («esto explica por qué murieron
  cuatro estadísticos») es más fuerte de lo que se reportó, no más débil.
- **Incrementos de precio: `N_eff > N`.** Eso **no es una pérdida del 6.2×**: es
  **antipersistencia**, la firma del rebote bid-ask, y la tabla actual la esconde detrás de una
  columna de «pérdida». Es un resultado, y coherente con `ρ₁(retornos) = −0.216` de la v3.0.
- **La fila de `ΔLL` no cambia el veredicto del paso 2**, que se recalculó con la inflación (la
  columna correcta) y sobrevive.

⚠ **Y una reserva sobre el alcance del argumento.** `Var(media) ∝ N^(−γ)` es el resultado de Beran
**para la media muestral**. No se transporta sin más a estimadores de exponentes, a pendientes de
regresión ni a estadísticos de razón de verosimilitud. Que el `N_eff` de los signos sea 47 **no
implica** que el `ΔLL` de la v3.2 tenga 47 grados de libertad — y de hecho su propio escalado da
`γ_eff = 0.742`, muy distinto. **`N_eff` se mide por estadístico, midiendo el escalado de ese
estadístico. No se importa de la serie de signos.** Escribirlo como regla del proyecto.

### 3.3 `γ` y `H` deben ajustarse en la misma ventana de escala

Rangos actuales: `γ` en `[10, 2000]` ticks, `H` en `[256, 16384]` ticks. Solapamiento `[256,
2000]`: **menos de una década**. La relación `H = (2−γ)/2 − β` es asintótica en un régimen de
escala concreto y combinar exponentes de ventanas distintas no la satisface por construcción.

**Qué hacer:** reajustar ambos en la ventana común `[256, 2000]` y reportar `β` implícita como
**curva contra la ventana de ajuste**, no como punto. Si `β` implícita se mueve más que la
distancia entre las dos hipótesis al variar la ventana, **el §1 de la v3.3 no está midiendo nada**
y así se reporta.

### 3.4 El núcleo creciente se reclasifica como mala especificación

`β = −0.160` con `τ₀` pegada a su cota inferior y `D = G(K)/G(0) = 11.53`. `text_3.txt` es
categórico: `β` se reporta en la literatura en `(0, 1)`, típicamente `0.2–0.6`, y valores
negativos son **firma de mala especificación del estimador**, no de no estacionariedad legítima.
Y hay un argumento propio, que no depende del resumen: un núcleo que **crece sin cota** sobre
`[0, K]` implica impacto de mercado creciente indefinidamente, que es económicamente imposible
(arbitraje ilimitado).

**Consecuencia:** `D = 11.53` deja de citarse como «núcleo creciente medido» y pasa a citarse como
**«ajuste no identificado en el régimen `τ₀` en cota»**, que es lo que el propio §5.2 de la v3.2 ya
decía sobre la identificabilidad `β`–`τ₀`. No se abre un quinto estadístico para arreglarlo: se
degrada la afirmación.

---

## §4. `η̂`: un control que falta, y su prioridad real

`η̂ = 0.6200` cerró el marco de tick grande. La fórmula es correcta y se puede verificar sin el
artículo, por argumento de barreras:

> Tras un salto al alza, `X` queda en `(k+1/2+η)α`. Continuar exige recorrer `α`; alternar exige
> recorrer `2ηα`. Por parada óptima sobre martingala continua,
> `P(cont)/P(alt) = 2ηα/α = 2η`, de donde `η̂ = N_c/(2·N_a)`. Con `η = 1/2` las dos distancias se
> igualan y `N_c = N_a`.

La verificación contra Robert & Rosenbaum (2011) **sigue siendo obligatoria** para el rango de
validez, pero la fórmula no está en duda.

⚠ **Lo que sí está en duda es el conteo.** El argumento de barreras supone que **cada cambio de
precio es un evento independiente de cruce de barrera**. En Binance, una orden de mercado que
barre varios niveles genera **varios `aggTrade` consecutivos a precios distintos, mismo signo,
milisegundos aparte**. Eso entra al contador como una ristra de continuaciones cuando es **un solo
evento económico**, e infla `N_c` — es decir, infla `η̂` justo en la dirección que cerró la
compuerta. El modo de fallo 3 de la v3.3 cubrió los saltos de más de 1 tick; **no cubrió los
barridos disfrazados de continuaciones de 1 tick**.

**Control:** recomputar `η̂` colapsando cambios consecutivos que ocurran dentro de la misma ráfaga
de agresión (mismo signo de agresor y separación `< 50 ms`, umbral declarado antes de contar), y
reportar la sensibilidad de `η̂` al umbral en `{0, 10, 50, 200} ms`.

⚠ **Y su prioridad, sin adornos: BAJA.** El §3 de la v3.3 estableció **antes de correr nada** que
la reformulación en tick grande arregla la ciencia y no la economía. Aunque `η̂` bajara de 1/2 y el
marco se reabriera, el spread cruzado es el **0.26 %** del coste total y el spread implícito de
Dayri–Rosenbaum sería `2ηα ≈ 1.24` ticks: **0.019 pb contra 4 pb de comisión**. No mueve nada.
Se corre porque un veredicto publicado sobre un contador contaminado es deuda, no porque decida.

---

## §5. `C_respaldo`, re-priorizado

Sigue sin medir sobre `captura_v32` y sigue siendo el término que falta de la ecuación de coste
maker: sin él, maker gana siempre porque no paga nada al fallar.

⚠ **Pero su peso cae con el horizonte, y eso reordena su prioridad.** El término es
`(1−p)·C_respaldo`, con `p` la tasa de llenado. Medida a 60 s: `p = 41.3 %`, luego `1−p = 58.7 %`.
A un horizonte de 10–60 min, `p` sube mucho y `(1−p)` se hunde. **`C_respaldo` es crítico para la
ejecución a escala de segundos —que es la que el §1 puede estar a punto de enterrar— y secundario
para la banda de minutos a horas.**

**Decisión:** se mide de todas formas, porque `lastre(H)` del §1 lo necesita para no estar
subestimado, pero **después** del §1 y del §2, y con `p(H)` medida al horizonte que el §1 señale,
no a 60 s.

---

## §6. El micro-precio de Stoikov como tercer observable

`text_4.txt` lo trae y el proyecto no lo ha usado. Es la única idea nueva de los resúmenes que
ataca un defecto **declarado** de la v3.2, no una hipótesis nueva.

**El defecto:** `y_mid` es **96.67 % ceros exactos**. La verosimilitud gaussiana está mal
especificada en la marginal, la propia v3.2 lo declaró como limitación 1, y eso ataca
directamente a `μ̂` — un observable que casi siempre vale cero atenúa cualquier estimación de
media condicional. **Si `μ̂` está atenuado, `q90(|μ̂|) = 11.30` es un suelo, no una medición.**

**El micro-precio** se construye con lo que ya se persiste (`bookTicker`: precios y cantidades de
los dos mejores niveles), es martingala por construcción, y se mueve en **cada actualización de
libro**, no sólo cuando el precio cruza un tick.

**Lo que se pide, y sólo esto:**

1. Implementar el micro-precio y reportar **la fracción de ceros** frente al 96.67 % del punto
   medio. Si no baja sustancialmente, se abandona la idea aquí y no se sigue.
2. Si baja: reajustar **M0/M1/M2 sin cambiar nada más** sobre el micro-precio como observable, y
   reportar `q90(|μ̂|)` junto a los dos observables ya publicados. **Tres columnas, siempre**, que
   es la regla heredada.

⚠ **Declarado antes de correr:** esto **no puede** cerrar un factor 3.46. Un factor 1.5 por
corregir atenuación es optimista y aun así deja el paso 3 fallando a 128 s. Se hace porque la
verosimilitud mal especificada contamina **toda** cifra de `μ̂` que el §1 vaya a producir, no
porque rescate la economía.

---

## §7. Lo que NO se hace

- **No se busca un quinto estadístico** para permanente contra transitorio. El §3.4 **degrada** la
  afirmación existente en vez de intentar mejorarla.
- **No se toca `Micelio.py`.** Sigue sin cambios desde la v2.2.
- **No se reabre el conjunto de prueba de la v3.2.** Ya se abrió dos veces. Todo el §1 corre sobre
  entrenamiento y validación.
- **No se reajusta M2 con `K` de 10 min.** Exige ~62 h continuas y el §1 está diseñado
  explícitamente para no necesitarlo.
- **No se resucita `ω_m` ni se decide `ω_G`.** `ω_m,max` y `γ_ω` siguen en el limbo.
- **No se persigue la reducción de comisiones.** El techo realista es −20 % y hace falta −71 %.
  Escrito aquí para que no vuelva a la mesa.

---

## §8. Modos de fallo

| # | fallo | firma | detección |
|---|---|---|---|
| 1 | Se elige el `H_p` que hace cruzar las curvas | cruce que sólo existe con el valor favorable | §2.2, las dos curvas requeridas siempre juntas |
| 2 | `R²` medido sobre ventanas solapadas | `R²` que no baja al aumentar `H` | ventanas no solapadas, `n_ventanas` en toda fila |
| 3 | Cruce declarado con 16 ventanas | resultado que no replica en `captura_estacional` | umbral `n_ventanas ≥ 30` escrito antes |
| 4 | Control barajado ausente a `H` grande | `R²` de 1 % sobre pocas ventanas leído como señal | control barajado obligatorio en toda fila |
| 5 | `σ₁` en USD/BTC a precio envejecido | `R²_req` mal por el cuadrado del error de precio | §2.1 punto 3, todo en pb |
| 6 | `N_eff` importado de la serie de signos a otro estadístico | IC mal por órdenes de magnitud | §3.2, `N_eff` se mide por estadístico |
| 7 | `η̂` sobre barridos | veredicto de compuerta sobre contador contaminado | §4, barrido de umbral |
| 8 | Micro-precio adoptado como observable primario porque «sale mejor» | tres observables, se reporta el que gana | §6, tres columnas siempre, ninguna se retira |
| 9 | `β` implícita presentada como medición independiente del núcleo | coherencia con fGn leída como evidencia | §3.1, la identidad escrita en el reporte |

---

## §9. Criterios de aceptación

**§1 — Las dos curvas**
- [ ] Tabla completa de `H` con `R²_medido`, `R²_barajado`, las dos `R²_req`, y `n_ventanas`.
- [ ] Ventanas no solapadas, verificado por test.
- [ ] Conjunto de prueba de la v3.2 no abierto, verificado por hash de partición.
- [ ] Regla del §1.4 aplicada literalmente y el desenlace nombrado.

**§2 — `H_p`**
- [ ] Firma de volatilidad sobre las cuatro capturas, con control barajado en todas.
- [ ] Reloj declarado por fila; el de pared es el que entra en `σ(H)`.
- [ ] `σ₁` reexpresado en pb.
- [ ] Extrapolación más allá de 711 s marcada en la propia tabla.

**§3 — Correcciones**
- [ ] Tercer estimador de `γ` (Whittle local o GPH) con su predicción de sesgo declarada antes.
- [ ] `β` implícita como curva contra ventana de ajuste, no como punto.
- [ ] Tabla de `N_eff` recalculada con la definición coherente y las tres consecuencias del §3.2.
- [ ] Identidad `β = 0 ⟺ γ = 2−2H` escrita en el reporte.
- [ ] `D = 11.53` reclasificado en todo texto donde aparezca.

**§4/§5/§6**
- [ ] `η̂` con barrido de umbral de colapso `{0, 10, 50, 200} ms`.
- [ ] `C_respaldo` con `p(H)` al horizonte que señale el §1.
- [ ] Fracción de ceros del micro-precio contra el 96.67 %; tres columnas si continúa.

**Transversal**
- [ ] `Micelio.py` sin cambios (`git diff --stat` vacío en el reporte).
- [ ] Suites en verde.
- [ ] Toda fórmula tomada de `text*.txt` verificada contra fuente primaria, con cita.

---

## §10. Prioridad

1. **§2 — `H_p`.** Barato, sobre dato en disco, y la curva requerida del §1 no se puede dibujar
   sin él. Horas.
2. **§1 — Las dos curvas.** El experimento que decide el proyecto. Todo lo demás es subordinado.
3. **§3 — Correcciones.** Baratas, y hay documentos citando las cifras malas ahora mismo.
4. **§6 — Micro-precio.** Sólo el punto 1 (fracción de ceros) es barato; el punto 2 sólo si pasa.
5. **§5 — `C_respaldo`.** Después de que el §1 diga a qué horizonte medir `p`.
6. **§4 — `η̂`.** Deuda de reporte, no decisión.

---

## §11. La compuerta de infraestructura que bloquea todo lo demás

`captura_estacional` lleva corriendo desde el 2026-08-10 16:30:53 UTC, 21 días previstos. **Es el
único activo del proyecto que puede dar potencia al §1 por encima de 30 min.**

⚠ **Un solo hueco la parte y no hay reintento barato: son tres semanas.**

- **Desactivar la suspensión automática de la máquina.** Está declarado desde la v3.2 y sigue
  siendo el criterio que salió en rojo en la compuerta.
- **Verificar hoy, y cada día, `n_huecos` y `trades_perdidos`** contra el detector de huecos por
  `aggTradeId` de la v2.0 §3.4. Un hueco silencioso descubierto el día 21 cuesta 21 días.
- **Reportar el tramo continuo más largo acumulado** en cada revisión, no el total capturado. Es
  el tramo continuo el que fija el `H` máximo medible: `H_max ≈ T_continuo / 30`.

Con 21 días continuos, `H = 1 h` da ~500 ventanas no solapadas y `H = 4 h` da ~126. Con un hueco a
mitad, ambas cifras se parten.
