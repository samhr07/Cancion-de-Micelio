# Orden de Trabajo — Convergencia de la Cadena de Medición de ω_m v2.1

Sucede a `ORDEN_TRABAJO_RELOJES_2_0.md`, ejecutada y verificada (42/42, `CLAUDE.md` sesión
2026-08-05). Nada de lo que construyó se deshace.

**Tesis de esta orden de trabajo:** el sistema ya no descarta datos ni corrige dos veces con la
misma medición, pero `ω_m` **todavía no es una medida del mercado**. Es en buena parte una
medida de la rejilla de muestreo. Todo lo que cuelga de `ω_m` —`ρ_k`, `Φ`, `Ψ`, `Ω`, `κΩ²`,
`μΩ²`, `Ω_crit`, la rama de `A`— hereda esa contaminación. Esta tanda cierra la cadena de
medición antes de volver a decidir nada sobre el modelo.

**El A/B queda congelado.** Razón: depende de `ω_m` y de `C`, que son justo lo que aquí se
repara. Correr 24 h ahora daría una medida muy precisa de la rejilla. El §8 fija la condición
**medible** de descongelamiento, para que no sea una decisión de criterio.

---

## §0. Registro de errores previos

| Origen | Error | Estado |
|---|---|---|
| v1.3 §D.4 (documento) | `C_ON = 0.50`, `C_OFF = 0.35` propuestos "**[CALIBRAR]** contra la distribución real de `C`" — pero se usaron como si fueran operativos, y la distribución no existía | La v2.0 mide `C` limpia en **0.624** contra **0.917** con escalera. Los `C` de 0.94–1.00 que mantenían enganchado el armónico eran en buena parte artefacto. Umbrales **invalidados**, ver §4.4 |
| v2.0 §1.4 (documento) | Se afirmó que el desajuste de tasas "no es una violación de Nyquist" y que no había aliasing en la cadena de estimación | Cierto para el lazo de control. **Falso para la cadena EMD**: el §5.2 demuestra que la rejilla determinaba la respuesta. La afirmación era correcta en su alcance y se leyó más ancha de lo que era |

---

## §1. Diagnóstico

### 1.1 Tres rejillas, tres respuestas, sin orden

| vía | duplicados | período | `C` |
|---|---|---|---|
| A — pared 0.5 s sobre precio a 1 Hz | 74.2 % | 95.2 s | 0.917 |
| B — pared 0.5 s con precio siempre fresco | 50.1 % | 262.6 s | 0.737 |
| C — una muestra cada K transacciones | 44.1 % | **13.2 s** | 0.624 |

Factor 20× entre extremos. Y el detalle que lo delata: **no están ordenadas**. La vía B tiene
mejores datos que la A —precio fresco, la mitad de duplicados— y se aleja de C. Arreglar los
duplicados sin arreglar el espaciado empeoró la respuesta, que es lo que ocurre cuando dos
defectos se cancelaban parcialmente. Firma de una medición donde el aparato domina sobre el
fenómeno.

### 1.2 La banda de resolubilidad — las dos vías de pared eran inreportables

La EMD no puede resolver un período más largo que ~W/3 (necesita al menos 3 ciclos en la
ventana para una IMF fiable) ni más corto que ~5 espaciados (necesita extremos para interpolar
envolventes). Aplicado a las tres vías, con W = 384 en todas:

| vía | Δτ | banda válida | período medido | ciclos en ventana | veredicto |
|---|---|---|---|---|---|
| A | 0.500 s | 2.50 – 64.00 s | 95.2 s | 2.02 | **fuera de banda** |
| B | 0.500 s | 2.50 – 64.00 s | 262.6 s | **0.73** | **fuera de banda** |
| C | 0.519 s | 2.59 – 66.37 s | 13.2 s | 15.08 | dentro |

En la vía B **el "período" es más largo que la ventana de observación**. Eso no es un ciclo: es
una tendencia — precisamente lo que `R_n` debe absorber — devuelta como si fuera la señal.

**Consecuencia operativa:** una guarda de tres líneas habría rechazado ambas respuestas sin
necesidad de la comparación a tres bandas. Es la adición más valiosa de esta orden de trabajo,
y su ausencia es la razón de que un `ω_m` estructuralmente inválido circulara durante tres
sesiones sin que nada protestara.

### 1.3 El 44 % residual no es un defecto de muestreo

La vía C conserva 44.1 % de duplicados, y son genuinos: con `tickSize = 0.10 USD`, a 27 tx/s y
K = 14, el precio de BTC no se mueve en 14 transacciones.

Es decir: **a resolución de transacción el precio no es una señal continua, es un proceso de
saltos sobre una retícula.** La HHT presupone lo contrario — la EMD interpola envolventes con
splines cúbicos y la transformada de Hilbert de un escalón tiene contenido en todo el espectro.
No se arregla muestreando mejor, porque no es un problema de muestreo: es que el observable
elegido está cuantizado. Ver §3.

### 1.4 ν varía por un factor 20

18.8 tx/s (REST, tranquilo) · 26–32 tx/s (WS) · **517.9 tx/s** media hora después, el mismo día
y el mismo par. ν es el puente entre relojes, así que esa variación entra en `σ_rel_tick =
σ_rel_s/√ν` (factor 4.5), en `K`, en `ω_m = f/ν` y en `ΔQ*`.

**Cualquier constante medida en un solo instante es sospechosa por construcción**, incluido el
`ΔQ* = 386.40 USD` de la v2.0. Las constantes de esta capa se miden como **distribuciones**,
no como valores, y se reportan con percentiles.

### 1.5 La cola: qué rompe y qué no

`y0`: curtosis **+207.8**, asimetría −9.52, Shapiro W = 0.09993.

Precisión importante, porque suena peor de lo que es: **el filtro de Kalman sigue siendo el
mejor estimador lineal insesgado aunque el ruido no sea gaussiano**, si los dos primeros
momentos están bien. Las estimaciones de estado no son basura.

Lo que sí queda roto es la **cuantificación de incertidumbre**, y con ella todo lo que cuelga
del NIS: el criterio de burn-in, la ventana adaptativa del EMD y las guardas. `NIS ~ χ²_m`
presupone gaussianidad; con curtosis 207 la distribución real no se parece a una χ². Ver §6.

---

## §2. Guarda de banda para ω_m — PRIORIDAD 1

### 2.1 La guarda

```python
def omega_en_banda(periodo, dtau, W, ciclos_min=3.0, muestras_por_ciclo_min=5.0):
    """§1.2: la EMD solo puede resolver períodos dentro de una banda fijada por
    la ventana y el espaciado. Fuera de ella la respuesta es un artefacto de la
    rejilla, no una medida del mercado.

    Cota superior: se necesitan >= 3 ciclos en la ventana para una IMF fiable.
    Cota inferior: se necesitan >= 5 muestras por ciclo para localizar extremos;
    2 (Nyquist) no basta porque la EMD interpola envolventes, no reconstruye.
    """
    T_min = muestras_por_ciclo_min * dtau
    T_max = W * dtau / ciclos_min
    return T_min <= periodo <= T_max
```

`ciclos_min` y `muestras_por_ciclo_min` son **[CALIBRAR]** — los valores de arriba son el punto
de partida de la literatura de EMD, no una medición de este sistema. Barrerlos contra la
captura del §5.2 y elegir por evidencia.

### 2.2 Qué hacer cuando cae fuera — no fabricar un valor

Sigue el patrón de `hay_medicion` de la v2.0: **se propaga la invalidez, no se inventa un
número.** Campo nuevo `omega_valida` (uint8) en `MICELIO_DTYPE`.

⚠ **La guarda tiene que alcanzar a `Φ`, `Ψ` y `Ω`, no solo a la rama de `A`.** Un `ω_m` fuera
de banda que se cuele en `Φ = ΣQ·ω_m` propaga a `Ω`, y `Ω` entra **al cuadrado** en el costo
del NMPC (`κΩ²`, `μΩ²`) y en `Ω_crit`. Es la vía por la que un artefacto de muestreo llega a
mover dinero.

Comportamiento exigido con `omega_valida = 0`:

| consumidor | comportamiento |
|---|---|
| rama de `A` | velocidad constante, incondicionalmente |
| `ρ_k` | término `γ_ω·ω_m` a cero; los demás siguen |
| `Φ`, `Ψ`, `Ω` | conservar último valor válido **con contador de antigüedad** |
| `Ω` con antigüedad > `T_OMEGA_RANCIA` | degradar: el NMPC pasa al comportamiento de burn-in |
| telemetría | `omega_valida`, `edad_omega`, y `ω_m` a **NaN**, no a cero |

`NaN` y no cero por la misma razón que la v2.0 puso `y1` a NaN: cero es un valor con
significado —"sin ciclo"— y confundirlo con "sin medición" es exactamente el defecto que
persigue este proyecto desde la primera sesión.

---

## §3. Observable menos cuantizado para la cadena EMD

### 3.1 El cambio

La cadena EMD pasa a consumir el **precio medio de `@bookTicker`**, `mid = (bid + ask)/2`. El
reloj de ticks y `z₀` del filtro siguen con `@trade`.

Razones:

- Resolución de **medio tick** (0.05 USD) en vez de un tick, por ser promedio de dos valores
  de retícula.
- Se actualiza cuando se mueve **cualquiera** de los dos lados, no solo cuando hay transacción.
  Medido en la v2.0: 91–288 msg/s.
- Es el observable que la Sec. 2 necesita —una señal suave que tamizar— mientras `@trade` es el
  que la Sec. 4.6 necesita —eventos que contar.

### 3.2 La separación no introduce un segundo reloj

`mid` es un **estado que se lee**, no un evento que se cuenta. El muestreo sigue siendo "cada K
transacciones": el reloj lo marca `@trade`, y en cada marca se lee el `mid` vigente. `ν` y `Δn`
no cambian de definición. Esto **no** viola la regla de puente único del §2.3 de la v2.0.

### 3.3 Riesgo conocido

El `mid` se mueve por parpadeo de cotizaciones sin que haya transacción, y puede oscilar entre
dos niveles por actividad de creadores de mercado sin contenido informativo. Es un ruido de
naturaleza distinta al de `@trade`. **Medirlo, no suponerlo**: reportar la fracción de cambios
de `mid` sin transacción intermedia. Si domina, reconsiderar.

### 3.4 Criterio de éxito

Fracción de duplicados a la misma `K`, `mid` contra precio de transacción. Referencia actual:
**44.1 %**. Si `mid` no baja de forma material, el observable no era el problema y hay que
volver sobre §1.3 con otra hipótesis.

---

## §4. σ_ω: separar la anchura instrumental de la genuina

### 4.1 Por qué antes de modelar

Tratar `ω_m` como distribución es correcto. Pero hay dos fuentes de anchura, y si construimos
el modelo sobre "ω es ancha" mientras buena parte de la anchura es instrumental, estaríamos
**incorporando el aparato a la física**.

### 4.2 La descomposición

```
σ_instr²  : dispersión de ω_m entre rejillas — MISMOS datos, variando K
            (K/2, K, 2K) y el desplazamiento de inicio de ventana
σ_total²  : dispersión de ω_m en el tiempo, a rejilla fija
σ_genuina² = σ_total² − σ_instr²
```

Si `σ_genuina²` sale **negativa**, es un resultado con contenido, no un fallo: significa que
toda la anchura observada es del aparato. Registrarlo como tal.

Reportar siempre la **fracción instrumental** `σ_instr²/σ_total²`. Es el número que gobierna el
§8.

### 4.3 σ_ω dentro de ρ_k — `DIVERGE DEL PDF (Sec. 7.3.3)`

Hoy `ρ_k = 1 + γ_ω·|ω_m| + γ_Q·|ΣQ|` infla `Q` con la **magnitud** de ω: dice "los ciclos
rápidos son más inciertos". Lo que se quiere es que infle con la **incertidumbre**: "los ciclos
mal conocidos son más inciertos". Con `ω` promovida a determinante de la dinámica en `A_arm`,
lo segundo es lo físicamente correcto.

La propagación de incertidumbre paramétrica es derivable, **sin constantes libres**:

```
J   = ∂A_arm/∂ω                    evaluada en Δn = 1
Q_ω = (J · x)(J · x)ᵀ · σ_ω²
```

Con

```
∂cos(ω)/∂ω      = −sin(ω)
∂(sin(ω)/ω)/∂ω  = (ω·cos(ω) − sin(ω))/ω²
∂(−ω·sin(ω))/∂ω = −sin(ω) − ω·cos(ω)
```

**Rama de Taylor obligatoria**, misma disciplina que `A_arm`: para `|ω| < 1e-6` usar
`−sin(ω) ≈ −ω`, `(ω·cos ω − sin ω)/ω² ≈ −ω/3`, `−sin ω − ω·cos ω ≈ −2ω`. Los tres tienden a
cero, así que `Q_ω → 0` cuando `ω → 0`: sin ciclo, no hay incertidumbre de ciclo. Es la
comprobación de coherencia de la construcción.

Implementar **ambas** variantes (`γ_ω·|ω_m|` y `Q_ω`), registrar las dos en telemetría, y no
decidir aquí. Decide el dato.

### 4.4 Los umbrales de `C` quedan invalidados — y hay un sustituto mejor

`C_ON = 0.50` / `C_OFF = 0.35` se fijaron contra una distribución que no existía, y la v2.0
demuestra que `C` estaba inflada por la escalera (0.917 → 0.624). Además `C` volverá a cambiar
al cambiar el observable (§3). Siguen **[CALIBRAR]**, y ahora contra una **distribución**, no
un punto.

**Propuesta a evaluar, no a imponer:** `C` es un *proxy* de "¿hay un ciclo bien definido?".
`σ_ω/ω` es una **medida directa de lo mismo**. Conmutar la rama de `A` por `σ_ω/ω` en vez de
por `C` sustituye un indicador indirecto por uno directo, y compone de forma natural con la
guarda de banda del §2: fuera de banda ⇒ inválida; dentro de banda con `σ_ω/ω` alta ⇒
velocidad constante; dentro de banda con `σ_ω/ω` baja ⇒ armónico.

Implementar el criterio nuevo **en paralelo** al de `C`, registrar ambos, decidir con datos.

---

## §5. Espacio logarítmico — con la justificación correcta

### 5.1 Lo que NO hace: limpiar la cola

Verificado numéricamente sobre una innovación con la estructura medida (mayoría de
transacciones sin mover el precio, cuantizadas a `tickSize`, unas pocas saltando):

| | curtosis | asimetría |
|---|---|---|
| innovación en precio | 9461 | 36.6 |
| innovación en log | **9413** | 36.4 |

Reducción del 0.5 %. La razón es estructural: **la curtosis es invariante bajo transformaciones
afines**, y `ln` es casi exactamente afín en el rango relevante. Las innovaciones abarcan
±0.41 % del nivel de precio; en esa franja el término cuadrático de `ln` pesa **2.05e-3**
respecto al lineal.

La log-normalidad limpiaría colas si vinieran de que el **nivel** de precio es multiplicativo.
Aquí vienen de una mezcla con masa puntual en cero por la retícula de `tickSize`. Es otro
mecanismo, y `ln` no lo toca.

**Anotar esto en el código** donde se implemente el cambio, para que en tres sesiones nadie
concluya que la cola quedó resuelta.

### 5.2 Lo que SÍ hace, y por lo que se adopta

En log-precio la volatilidad relativa es **nativa**. El apaño `q_S = (σ_rel·S)²` de la v1.2
—que existe porque un movimiento de 45k a 90k dejaba `q_S` mal escalada por 4×— desaparece:
`q_ln = σ_ln²`, sin factor de nivel.

### 5.3 Superficie de cambio — auditar como en la forma afín

`x[0]` pasa de precio a log-precio. Todo consumidor necesita `S = exp(x[0])`:

- `γ_0(S)` y la condición terminal de Loeper
- el centro de la malla de Loeper
- el NMPC y las restricciones de caja de la Sec. 6.2
- la abrazadera de nocional y la de posición (§B.4 de la v1.3) — **estas son de seguridad**
- telemetría

La cadena EMD pasa también a log: se tamiza `ln(mid)`, y `R_n` queda en log-precio. `z₀ =
ln(P_spot)`, `z₁ = R_n^ln`. Coherente de extremo a extremo o no se hace.

### 5.4 ⚠ El hueco de Jensen — decisión explícita, no accidente

`E[exp(x)] ≠ exp(E[x])`. Al volver a espacio de precio para el NMPC y para las guardas de
riesgo:

- `exp(x̂)` da la **mediana** de la distribución de precio
- `exp(x̂ + P₀₀/2)` da la **media**

No son lo mismo, y la diferencia crece con `P₀₀`. Con la incertidumbre típica el sesgo es
despreciable, **pero no en burn-in ni tras un halt**, que es justo cuando `P` está inflada.

Elegir una, documentarla con `# NOTA DE INTERPRETACION:`, y **publicar `P₀₀` junto al precio
reconstruido** para que la magnitud del hueco sea auditable en telemetría en vez de invisible.

---

## §6. Cuantiles empíricos en lugar de umbrales χ²

Todo lo que hoy compara NIS contra un umbral χ² está mal calibrado (§1.5): burn-in, ventana
adaptativa del EMD, guardas.

**El arreglo no es transformar los datos** —§5.1 muestra por qué no funcionaría— **sino medir
la distribución que ya existe.** Sustituir el umbral por el **cuantil empírico** de la
distribución observada de NIS sobre ventana deslizante, conservando la *intención* de diseño
(el percentil que el PDF pretendía) y descartando el *supuesto* (que ese percentil corresponde
a una χ²).

Arranque en frío: durante las primeras `N_NIS_BOOTSTRAP` muestras no hay distribución que medir.
Usar χ² como prior, **registrar el cambio a empírico en telemetría**, y no dejarlo implícito.

Reportar la razón `cuantil_empírico / cuantil_χ²`. Si es de orden 1, el supuesto gaussiano era
inocuo aquí y se anota; si es grande, cuantifica cuánto llevaban desviadas las compuertas.

---

## §7. La hipótesis del péndulo — ω como función del estado

### 7.1 El planteamiento

Un péndulo con arrastre cuadrático tiene período dependiente de la amplitud. Si el ciclo
estructural se comporta así, `ω_m` no es un parámetro ruidoso: es **función del estado**.

Y el diseño ya lo contempla: **`Ψ = Φ/ΔS = ΣQ·ω_m/ΔS` es el acoplamiento frecuencia-amplitud.**
La analogía no contradice la Sec. 1.4, la describe.

### 7.2 El test

Correlación entre `ω_m` y `ΔS` sobre la corrida larga.

⚠ **Controlando por ν.** `ω_m = f/ν` es mecánicamente antiproporcional a ν, y ν varía 20×
(§1.4). Sin controlar, aparecería correlación espuria con cualquier magnitud que covaríe con la
actividad — y `ΔS` covaría con la actividad. Usar correlación parcial `ρ(ω_m, ΔS | ν)`, y
reportar también `ρ(ω_m, ν)` y `ρ(f, ΔS)` para poder distinguir los casos.

### 7.3 Alcance

Esta orden de trabajo **mide**, no modela. Si la correlación parcial resulta material, la
promoción de `ω` a función del estado es materia de la v2.2 y toca la Sec. 1.4 del PDF.

---

## §8. El A/B, congelado — con condición medible de descongelamiento

Congelado por dependencia: el veredicto depende de `ω_m` y de `C`, que son lo que esta orden
repara. Correr 24 h ahora daría una medida muy precisa de la rejilla.

**Se descongela cuando se cumplan las cuatro:**

1. `σ_instr²/σ_total² < 0.5` (§4.2) — la anchura de `ω_m` está dominada por el mercado, no por
   el aparato.
2. Fracción de duplicados con `mid` **materialmente** por debajo del 44.1 % actual (§3.4).
3. `omega_valida = 1` en más del 90 % de las muestras de una corrida larga (§2).
4. `§5.2` repetido en **régimen agitado** (ν > 300 tx/s) además del tranquilo, con conclusión
   coherente. Un solo tramo de mercado no generaliza, como el propio reporte de la v2.0 advierte.

Que las cuatro sean números y no juicios es deliberado: es lo que impide que el A/B se
descongele por impaciencia.

---

## §9. Criterios de aceptación

**§2 — Guarda de banda**
- [ ] `omega_en_banda` implementada; barrido de `ciclos_min` y `muestras_por_ciclo_min` sobre la
      captura del §5.2, con la elección justificada por evidencia.
- [ ] **Test de regresión histórica**: alimentar la guarda con las tres vías del §5.2 y verificar
      que rechaza A y B y acepta C. Es el test que demuestra que el defecto de tres sesiones
      queda cerrado.
- [ ] `omega_valida` propagada a `Φ`, `Ψ` y `Ω`, no solo a la rama de `A`. **Test que lo fuerza**:
      inyectar ω fuera de banda y verificar que `Ω` no cambia de valor.
- [ ] `edad_omega` y degradación por rancidez implementadas y probadas.
- [ ] `ω_m` a NaN —no a cero— cuando es inválida, en telemetría y en el volcado.

**§3 — Observable**
- [ ] Cadena EMD sobre `mid` de `@bookTicker`; reloj y `z₀` siguen en `@trade`.
- [ ] Fracción de duplicados a igual `K`, `mid` contra transacción, reportada como número.
- [ ] Fracción de cambios de `mid` sin transacción intermedia, reportada (§3.3).

**§4 — σ_ω**
- [ ] Descomposición instrumental/genuina ejecutada; `σ_instr²/σ_total²` en el reporte.
- [ ] `Q_ω` implementada con rama de Taylor; **test de que `Q_ω → 0` cuando `ω → 0`**.
- [ ] Ambas variantes de `ρ_k` registradas en paralelo, sin decidir.
- [ ] Criterio `σ_ω/ω` registrado en paralelo a `C`, sin decidir.

**§5 — Espacio log**
- [ ] Cadena completa en log: `x[0]`, `z₀`, `z₁`, EMD sobre `ln(mid)`.
- [ ] **Auditoría escrita** de consumidores de `x[0]`, con atención explícita a las abrazaderas
      de riesgo de la v1.3. Lista en el commit.
- [ ] Decisión de Jensen documentada; `P₀₀` publicada junto al precio reconstruido.
- [ ] Nota en el código de que `ln` **no** resuelve la cola, con los números del §5.1.
- [ ] Test de ida y vuelta: `exp(ln(S)) == S` dentro de tolerancia, en todo el rango operativo.

**§6 — Cuantiles**
- [ ] Cuantil empírico deslizante sustituye el umbral χ² en burn-in, ventana adaptativa y guardas.
- [ ] Cambio de prior χ² a empírico registrado en telemetría.
- [ ] Razón `cuantil_empírico/cuantil_χ²` reportada.

**§7 — Péndulo**
- [ ] `ρ(ω_m, ΔS | ν)` parcial, más `ρ(ω_m, ν)` y `ρ(f, ΔS)` como control.

**Transversal**
- [ ] Suite v1.3 + v2.0 en verde: **42/42 + los nuevos**.
- [ ] Toda constante nueva reportada como **distribución con percentiles**, no como valor (§1.4).
- [ ] Nomenclatura de relojes del §2.4 de la v2.0 mantenida sin excepciones.

---

## §10. Fuera de alcance

- **El A/B**, hasta las cuatro condiciones del §8.
- **Fase 2 (ALS)**. El terreno mejora —`y0` ya está por debajo del umbral material y lo que
  queda es microestructura, que es lo que ALS debe estimar— pero con curtosis 207 el ALS
  estándar necesitará mucha muestra. La nota de la v1.2 sobre ALS-IRLS (Huber) se refuerza;
  aun así, primero el estándar, que sin línea base no hay mejora que medir.
- **Promoción de `ω` a función del estado** (§7.3): materia de la v2.2 si el dato lo pide.
- **El silencio de `@aggTrade`** en fstream. Sin explicación, sin bloquear nada. `@trade` cubre
  el caso y además es lo que §3 necesita.
- **Cadencia del Hilo Rápido a 518 tx/s.** Subió a 41.1 ms/ciclo (3.7×) por drenar el lote.
  Vigilar; si el drenado necesita hilo propio, es decisión de la v2.2.
- Fase 3 (`Ω_crit`), Testnet con credenciales, las 30 corridas.
- CMA-ES, IMPC, DeepONet.
- Reconciliación documental del PDF: Secs. 3.5 y 8.6.1 superadas por ALS; contradicción
  L1/cuadrática (4.5 contra 6.1, prevalece 6.1); dimensiones del NMPC en 7.5.3; Sec. 2.6 como
  `θ ≡ π/2 (mód π)`; Sec. 4.6 con `n_ticks` por transacción; y ahora la Sec. 2.2, que debe
  incorporar la banda de resolubilidad como condición de validez de `ω_m`.
