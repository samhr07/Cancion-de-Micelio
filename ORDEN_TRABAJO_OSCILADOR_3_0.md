# Orden de Trabajo — Oscilador Forzado y Presupuesto Energético v3.0

Sucede a `ORDEN_TRABAJO_EXISTE_OMEGA_2_2.md`. La v2.2 dio veredicto **(A) dentro de
2.5–300 s**: en esa banda no hay escala característica, y `ω_m` es salida del algoritmo, no del
mercado. La hipótesis **(B)** —escala más larga— sigue abierta y la captura de 48 h corre.

**Lo que este documento cambia:** deja de estimar `ω_m` como cantidad primaria. `ω₀ = √(k/m)`
es **derivada**; el artefacto de ventana vive entero en el paso "elegir ventana → extraer modo
dominante". Este documento mide las primitivas —`k`, `m`, `γ`— por regresión sobre todos los
datos, sin ventana espectral, sin banda de resolubilidad, con hipótesis nula.

**Lo que NO cambia:** la premisa física. El mercado como oscilador amortiguado bajo forzamiento
externo, con intercambio de energía, se conserva **estrictamente**. Lo que se conserva mejor
que antes, de hecho: aquí cada término tiene unidades y el balance energético se puede
**verificar**, cosa que hoy no existe en ninguna parte del diseño.

---

## §0. Registro de errores previos

| Origen | Error | Estado |
|---|---|---|
| v2.1 §3.4 | "Fracción de duplicados" como criterio de observable — un proxy que engañaba | Corregido en v2.2 §4.3 |
| v2.2 §5 (brazo IAAFT) | IAAFT aplicado al **nivel** de precio (paseo aleatorio), no converge: 11.6 s en un observable y 118 s en el otro | Rehacer **sobre incrementos**. La lectura de la v2.2 no se apoyaba en él, pero separar (A) de (A′) sigue pendiente |
| Diseño (todas las versiones) | `ω_m` tratada como cantidad primaria medible | **Es derivada.** `k`, `m`, `γ` son las primitivas |
| Diseño (todas las versiones) | `S_ref` elegida por preproceso (nodos de fase), y luego `ΔS` usada para medir | Circularidad. El §4 la elimina: tendencia y ciclo se estiman **conjuntamente** |

Nota metodológica: la forma de onda `e^{sin(ω τ)}` **no explica el nulo de la v2.2**. Su
descomposición en armónicos (Bessel modificadas, confirmado por FFT) concentra el **94.4 %** de
la potencia AC en la fundamental, con el 2º armónico **12.4 dB** por debajo. Para un estimador
espectral es prácticamente una sinusoide: el multitaper la habría visto. La no sinusoidalidad
queda descartada como explicación y no debe reabrirse.

---

## §1. La premisa, en forma falsable

### 1.1 El modelo

```
m ẍ  +  γ ẋ  +  k x  =  F(τ)              x = S − S_ref,  τ = tiempo de ticks
```

- `m` — inercia del precio: cuánto cuesta acelerarlo
- `γ` — disipación: pérdida por fricción de ejecución
- `k` — fuerza recuperadora hacia `S_ref`
- `F(τ)` — forzamiento externo

### 1.2 La condición que la hace falsable

Un oscilador forzado cuya `ω` varía libremente de ventana en ventana **no tiene contenido
empírico**: cualquier serie lo satisface, incluido el paseo aleatorio, que es exactamente lo que
hizo la EMD. La física pone la restricción:

- **Lineal:** `ω₀ = √(k/m)` es propiedad del sistema y **no la cambia el forzamiento**. El
  forzamiento cambia amplitud y fase. El espectro muestra resonancia de anchura `γ/m`.
- **No lineal** (Duffing, péndulo): `ω = ω(A)` varía, pero **como función de una cantidad
  observable**, no libremente.

La versión falsable de la premisa es la segunda, y su predicción concreta es que **`ω` debe
correlacionar con la amplitud**. Por eso el §7 de la v2.1 pasa aquí de bloqueado a **central**
(§6).

**Condición de coherencia, para el registro:** `ω` debe variar más despacio que la ventana
necesaria para medirla. Si cambia a la misma escala a la que se mide, es inobservable por
construcción, y el modelo vuelve a ser infalsable.

---

## §2. El presupuesto energético

### 2.1 El balance

```
E        = ½ m ẋ²  +  ½ k x²
dE/dτ    = ẋ(m ẍ + k x) = F·ẋ − γẋ²
```

```
P_ext  = F · ẋ        potencia inyectada por el forzamiento
P_diss = γ ẋ²         potencia disipada
```

**Los tres términos son medibles por separado.** Ésa es la propiedad que ninguna versión
anterior del diseño tenía: si el balance no cierra, el modelo está mal, y hay una prueba que lo
dice.

### 2.2 Análogos de mercado — candidatos, con su razonamiento

⚠ **Ninguna de estas identificaciones se acepta sin derivación dimensional.** Se listan con su
argumento; el §2.3 fija el procedimiento.

| término | candidato | razonamiento |
|---|---|---|
| `ẋ²` | **varianza realizada** por tick | `ẋ` es velocidad de precio; su cuadrado integrado es varianza realizada. Usar **núcleos realizados** (Barndorff-Nielsen/Hansen/Lunde/Shephard) o estimador de dos escalas: el ruido de microestructura aquí es severo (curtosis 601.8, retícula `tickSize`, rebote bid-ask con ρ₁ = −0.164) |
| `m` | **profundidad del libro**, candidato `∝ 1/λ` | Inercia = resistencia a la aceleración. Un libro profundo cuesta más mover: es literalmente inercia. `λ` ya existe en el sistema |
| `γ` | **impacto temporal / spread efectivo** | La disipación es pérdida irreversible por fricción de ejecución. El spread efectivo sale del modelo de Roll, y el ρ₁ = −0.164 ya medido **es** esa medición |
| `F` | **desbalance de flujo de órdenes** (volumen firmado) | Es la fuerza externa que empuja el precio. `es_maker` **ya se captura y no se usa** |
| `k` | **por medir** | Es el hueco. Sin `k` no hay pozo de potencial y no hay oscilador. §3 |

### 2.3 Procedimiento de derivación dimensional — obligatorio antes de implementar

1. Fijar las unidades de `x`, `ẋ`, `ẍ` en tiempo de ticks, con la nomenclatura de relojes del
   §2.4 de la v2.0.
2. Extraer de la Sec. 4 del PDF las **unidades reales de `λ`** — no asumirlas. La comprobación
   de que `λS²Γ` es adimensional fija `λ`, y de ahí sale si `m ∝ 1/λ` es sostenible o necesita
   un factor.
3. Verificar que los cuatro términos de `m ẍ + γẋ + kx = F` tienen **unidades idénticas**.
4. Verificar que `Q = √(km)/γ` sale **adimensional**. Si no, la identificación está mal.
5. Test automático que falla si cualquiera de las anteriores se rompe. Las trampas del 2π y del
   factor 125 fueron errores dimensionales silenciosos; aquí hay cuatro identificaciones nuevas
   a la vez.

---

## §3. Medir `k` y `Q` — el test con hipótesis nula

### 3.1 El oscilador discretizado es un AR(2)

Con `Δτ = 1 tick` (`Δn = 1`, condición ya vigente desde la v2.0):

```
x_t = φ₁ x_{t-1} + φ₂ x_{t-2} + ε_t
```

y, salvo una escala común `D = m + γ + k` que el forzamiento fija:

```
m ∝ −φ₂            γ ∝ φ₁ + 2φ₂            k ∝ 1 − φ₁ − φ₂
```

**`Q` es identificable sin conocer `D`**, que es lo que lo hace útil:

```
Q = √(k·m) / γ = √[(1 − φ₁ − φ₂)(−φ₂)] / (φ₁ + 2φ₂)
```

Condición de oscilación: raíces complejas, `φ₁² + 4φ₂ < 0`, equivalente a `Q > ½`.

### 3.2 Validación numérica — hecha, reproducir como test

| caso | Q real | Q recuperada por AR(2) | raíces | `k` recuperada (real 0.25) |
|---|---|---|---|---|
| subamortiguado | 5.000 | **5.094** | complejas | 0.2519 |
| subamortiguado | 1.515 | **1.500** | complejas | 0.2491 |
| crítico | 0.500 | **0.501** | complejas (frontera) | 0.2523 |
| sobreamortiguado | 0.167 | **0.166** | **reales** | 0.2529 |
| **paseo aleatorio** | — | **5.6e-05** | **reales** | **4.65e-06** |

En el paseo aleatorio: `φ₁ = 1.00066`, `φ₂ = −0.00066` — el paseo puro sería `(1, 0)`. **`k → 0`
⇒ sin fuerza recuperadora ⇒ sin oscilador.**

Es decir: **el nulo de este test es exactamente un paseo aleatorio**, que es la hipótesis que la
v2.2 no pudo rechazar. Por primera vez la pregunta se hace con un estadístico cuyo nulo es la
alternativa que nos preocupa, en lugar de con una herramienta que siempre devuelve un modo.

### 3.3 Lo que hay que añadir a la validación de arriba

- **Errores estándar de `φ₁`, `φ₂`** y propagación a `Q` y `k` (delta method o bootstrap por
  bloques, que respeta la dependencia serial).
- **Contraste explícito `H₀: k = 0`** con su p-valor. Es el resultado principal del documento.
- **Robustez a la cola.** Con curtosis 601.8, mínimos cuadrados es ineficiente y sus errores
  estándar son optimistas. Ajustar también por **IRLS/Huber** y reportar ambos. Si difieren de
  forma material, prevalece el robusto.
- **`γ < 0`** significa inyección neta de energía: o el sistema es inestable, o la
  discretización está mal. **Es una guarda, no un resultado**: si sale, parar y revisar antes de
  interpretar.
- **Estabilidad temporal:** repetir por bloques a lo largo de las 48 h. Si `k` y `Q` cambian de
  bloque en bloque tanto como cambiaba `ω_m` de ventana en ventana, no hemos avanzado y hay que
  decirlo.

---

## §4. El ciclo estocástico de Harvey dentro del EAKF

### 4.1 Por qué es la implementación de la premisa, no una conveniencia

```
⎡ψ_t ⎤       ⎡ cos λ   sin λ⎤ ⎡ψ_{t-1} ⎤   ⎡κ_t ⎤
⎢    ⎥ = ρ · ⎢              ⎥ ⎢        ⎥ + ⎢    ⎥        ρ ∈ (0,1)
⎣ψ*_t⎦       ⎣−sin λ   cos λ⎦ ⎣ψ*_{t-1}⎦   ⎣κ*_t⎦
```

Es la forma en espacio de estados del oscilador amortiguado con forzamiento aleatorio, término
por término: `ρ` es el amortiguamiento, `λ` la frecuencia natural, y **`κ_t` es el forzamiento
externo**. Correspondencia con el §3: `ρ = √(−φ₂)`, `λ = arccos(φ₁/(2√(−φ₂)))`.

Cuatro propiedades que atacan directamente lo hallado en las últimas cuatro sesiones:

1. **Nulo nativo.** `ρ → 0` o `σ_κ² → 0` ⇒ el ciclo no existe. La pregunta de la v2.2 se
   responde de forma **continua y sobre todos los datos**, no por experimento aparte.
2. **`σ_ω` sale de la covarianza de parámetros.** Sin la descomposición instrumental/genuina
   del §4.2 de la v2.1.
3. **Desaparece el traspaso externo `ω_m → A_arm`**, donde viven la trampa del 2π, la guarda de
   banda y el artefacto de rejilla. **Sin ventana de análisis no hay artefacto de ventana.**
4. **Elimina la circularidad de `S_ref`.** El modelo estructural estima **tendencia + ciclo +
   irregular conjuntamente**: la tendencia absorbe `S_ref` en vez de fijarla por preproceso.
   Los nodos de fase (Sec. 2.6) dejan de ser necesarios para definir `ΔS`.

### 4.2 Coste, que es real

`A` pasa a depender de parámetros a estimar (`λ`, `ρ`), así que el filtro deja de ser lineal en
parámetros. Tres vías, de menor a mayor ambición:

| vía | descripción | recomendación |
|---|---|---|
| **(a)** | `λ`, `ρ` fijos por ML fuera de línea sobre la captura de 48 h; el filtro corre lineal | **Empezar aquí.** No toca la arquitectura y da el veredicto |
| (b) | Rejilla de `(λ, ρ)` con banco de filtros y pesos por verosimilitud | Si (a) muestra que los parámetros varían de forma material |
| (c) | Estado aumentado con EKF sobre `(x, λ, ρ)` | Solo si (b) resulta insuficiente. Reintroduce no linealidad en el filtro que hoy no existe |

**Esta orden de trabajo implementa (a) y solo (a).** El filtro de producción no se toca hasta
tener veredicto.

### 4.3 Comparación obligatoria

Verosimilitud del modelo estructural **con** componente cíclica contra **sin** ella (solo
tendencia + irregular), sobre los mismos datos. Contraste de razón de verosimilitudes.

⚠ El contraste está en la frontera del espacio de parámetros (`σ_κ² = 0`), donde la
distribución asintótica **no es** χ² estándar sino una mezcla. Usar la corrección adecuada o
bootstrap paramétrico. No reportar un p-valor χ² ingenuo: sería exactamente el tipo de veredicto
sobre supuestos no verificados que este proyecto ya ha tenido que anular dos veces.

---

## §5. Rediseño de `Φ`, `Ψ`, `Ω` desde el balance energético

Se derivan del §2.1; no se postulan.

### 5.1 `Ω` → factor de calidad

```
Ω ≡ Q = √(k·m)/γ
```

Adimensional, físicamente interpretable, y con nulo en `Q = ½` (frontera subamortiguado /
sobreamortiguado). **Sustituye a `C`** (concentración espectral), que era un proxy indirecto de
la misma pregunta y que la v2.2 demostró inflado por artefacto (0.917 con escalera contra 0.624
limpia, y 0.783 sobre datos **barajados**, donde no hay estructura alguna).

### 5.2 `Ω_crit` → condición energética

```
Ω_crit  ≡  P_ext / P_diss = 1
```

Si el forzamiento supera a la disipación, el sistema gana energía y el régimen cambia. Es una
condición física, medible y con significado, no un umbral ajustado.

### 5.3 `Φ` y `Ψ`

Se rederivan del balance `dE/dτ = P_ext − P_diss`. **Esta orden de trabajo no fija su forma
final:** primero hay que saber si `k > 0`. Si `k = 0` no hay energía potencial, y `Φ`/`Ψ`
necesitan otra base — probablemente el desbalance de flujo directamente, sin pasar por el
oscilador.

**Restricción de diseño, no negociable:** cualquier forma propuesta debe tener unidades
verificables y entrar en el balance del §2.1. Nada que no se pueda contrastar contra la
ecuación de energía.

### 5.4 Migración

Nada de la Sec. 1.4 actual se retira hasta tener el veredicto del §3. Mientras tanto, `Ω` sigue
gobernada por la guarda de banda de la v2.1 §2 —que, con `omega_valida` al 12.5 %, la mantiene
congelada o en cero la mayor parte del tiempo. Es el comportamiento correcto: la capa no tiene
insumo y la guarda lo dice.

---

## §6. El test de amplitud — ahora central

De bloqueado en la v2.2 a resultado principal, porque **es la predicción que distingue la
premisa no lineal de la nula**.

- Correlación parcial `ρ(ω, A | ν)`, con `A` la amplitud de la componente cíclica estimada en
  el §4 y `ω` su frecuencia instantánea.
- ⚠ **Controlar por ν es obligatorio.** `ν` varía por factor 20 y arrastra a casi todo.
  Reportar también `ρ(ω, ν)` y `ρ(A, ν)` para poder leer la parcial.
- Si `k = 0` sale del §3, este test **no se ejecuta**: no hay oscilador cuya frecuencia
  correlacionar. Es la única dependencia de orden dentro del documento.

Además, con `k` y `Q` estimados por bloques: **¿covaría `Q` con el régimen de volatilidad?** Un
oscilador que se sobreamortigua cuando sube la volatilidad es una predicción física concreta y
contrastable, y es la versión medible de tu observación de que un armónico muy perturbado no lo
parecerá.

---

## §7. Qué decide y qué confirma

| ítem | papel |
|---|---|
| **§3 — `k` y `Q` por AR(2) con `H₀: k = 0`** | **DECIDE.** Todo lo demás depende de esto |
| §4 — Harvey vía (a), razón de verosimilitudes | **DECIDE**, por segunda vía independiente |
| §6 — amplitud–frecuencia | Decide entre premisa lineal y no lineal, **si** `k > 0` |
| v2.2 §4 multiescala sobre 48 h | **Confirma.** Sigue abierta la hipótesis (B) |
| Contraste de periodicidad de **calendario** (funding 8 h, sesiones, diario) | **Confirma**, y es barato: período conocido *a priori* ⇒ se contrasta, no se estima. Hipótesis única con mucha más potencia que una búsqueda libre. Verificar el intervalo de funding vigente de BTCUSDT contra la documentación |
| Wavelet con significancia Torrence-Compo | **Confirma.** Único hueco de la banda corta: el multitaper es global y promediaría un ciclo **intermitente** |
| IAAFT rehecho sobre incrementos | **Confirma.** Separa (A) de (A′): "no hay ciclo" contra "hay memoria, no ciclo" |

---

## §8. Criterios de aceptación

**§2 — Dimensiones**
- [ ] Los cinco pasos del §2.3, con las unidades de `λ` **extraídas del PDF**, no asumidas.
- [ ] Test automático: los cuatro términos de la ecuación con unidades idénticas; `Q`
      adimensional.
- [ ] Varianza realizada con núcleo realizado o estimador de dos escalas, **no** suma ingenua de
      cuadrados. Reportar el **gráfico de firma de volatilidad**: dice empíricamente a qué `Δτ`
      deja de dominar el ruido de microestructura, que es la misma pregunta que la elección de
      `K` y hoy se resuelve por conjetura.

**§3 — El test decisivo**
- [ ] AR(2) sobre la captura de 48 h, tres observables, con la nomenclatura de relojes vigente.
- [ ] **`H₀: k = 0` con p-valor y errores estándar por bootstrap por bloques.**
- [ ] Reproducir la tabla de validación del §3.2 como test permanente, **incluida la fila del
      paseo aleatorio**: es la demostración de que el nulo es el correcto.
- [ ] Ajuste robusto IRLS/Huber en paralelo a MCO; ambos reportados.
- [ ] Guarda `γ < 0`: si dispara, **parar y revisar**, no interpretar.
- [ ] `k` y `Q` por bloques a lo largo de las 48 h, con su dispersión. Si varían tanto como
      variaba `ω_m`, decirlo explícitamente.

**§4 — Harvey**
- [ ] Vía (a) únicamente. **El filtro de producción no se toca.**
- [ ] Razón de verosimilitudes con la corrección de frontera o bootstrap paramétrico. **Un
      p-valor χ² ingenuo es motivo de rechazo del entregable.**
- [ ] Coherencia con el §3: `ρ = √(−φ₂)`, `λ = arccos(φ₁/(2√(−φ₂)))` dentro de tolerancia. Dos
      caminos al mismo número; si discrepan, uno está mal.

**§5 — Rediseño**
- [ ] `Q` y `P_ext/P_diss` calculados y en telemetría, **sin sustituir todavía** a `Ω` ni a
      `Ω_crit`.
- [ ] **Verificación del balance energético**: `dE/dτ` medida contra `P_ext − P_diss` medida,
      con su residuo. Es la prueba que el diseño nunca tuvo, y su resultado vale por sí solo
      aunque `k = 0`.

**Transversal**
- [ ] Suite en verde: **53/53 + los nuevos**.
- [ ] `Micelio.py` **sin cambios de modelo**. Todo en código de análisis aparte, como en la v2.2.
- [ ] Cada identificación física del §2.2 marcada con `# NOTA DE INTERPRETACION:` y su
      derivación dimensional en el comentario.

---

## §9. Fuera de alcance

- **Vías (b) y (c) del §4.2.** Solo si (a) lo justifica.
- **Retirar la Sec. 1.4 actual.** No antes del veredicto del §3.
- **Forma final de `Φ` y `Ψ`** (§5.3). Depende de si `k > 0`.
- **Cambiar el observable de producción.** El §4.3 de la v2.2 lo evalúa.
- **v2.1 §5 (log-precio) y §6 (cuantiles empíricos del NIS).** Siguen en cola y **no
  bloqueadas**: no dependen de `ω`. Las compuertas del NIS están mal calibradas pase lo que pase
  —curtosis 601.8 contra una χ²— y son el mejor trabajo para la espera de las 48 h.
- **Fase 2 (ALS), Fase 3 (`Ω_crit`), Testnet con credenciales, las 30 corridas.**
- CMA-ES, IMPC, DeepONet.
- El silencio de `@aggTrade` en fstream.
- Reconciliación documental del PDF. La lista ya justifica sesión propia: Secs. 3.5 y 8.6.1
  superadas por ALS; contradicción L1/cuadrática (4.5 contra 6.1, prevalece 6.1); dimensiones
  del NMPC en 7.5.3; Sec. 2.6 como `θ ≡ π/2 (mód π)`; Sec. 4.6 con `n_ticks` por transacción;
  Sec. 2.2 con la banda de resolubilidad; **Sec. 2 entera** pendiente del veredicto del §3; y
  **Sec. 1.4** pendiente del §5.
