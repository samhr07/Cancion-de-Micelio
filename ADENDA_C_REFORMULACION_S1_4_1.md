# Adenda C — Reformulación del §1 de `ORDEN_TRABAJO_HORIZONTE_4_1.md`

Motivada por la ejecución del 2026-08-10 (c). **No es una enmienda al criterio de decisión**: la
regla del §1.4 se mantiene palabra por palabra. Lo que se sustituye es el **instrumento** con el
que se mide `R²_medido(H)`, porque la ejecución demostró que el instrumento diseñado no tiene
resolución en la banda donde vive la pregunta.

`Micelio.py` no se toca. El conjunto de prueba de la v3.2 no se abre.

---

## §C.0. Qué estableció la ejecución del §1, y qué no

**Lo que estableció, y se conserva sin tocar:**

A `H = 60 s` y `H = 120 s` la medición **sí** tiene resolución, y el resultado es un fallo
categórico: `R²_medido ≈ 0` contra un requerido de 38–83 %. El margen instrumento/requisito es de
13.8× y 3.0× respectivamente (§C.1). **La hipótesis del horizonte queda refutada por debajo de
2 minutos**, y esa refutación no depende de ninguna mejora del predictor: ni un predictor perfecto
(`R² = 1`) paga a 60 s si el requisito es 83 %.

**Lo que NO estableció, y la ejecución lo etiquetó correctamente como piloto:**

De `300 s` en adelante, los `R²` negativos y decrecientes **no son mediciones de ausencia de
señal**. Son la firma del error de estimación de los coeficientes, y son cuantitativamente
compatibles con ella:

| `H` | `n_train` | `−k/n_train` (ruido puro, `k = 4`) | `R²_medido` | razón |
|---|---|---|---|---|
| 300 s | 118 | −0.0339 | −0.0264 | 0.78 |
| 600 s | 59 | −0.0678 | −0.0554 | 0.82 |
| 900 s | 39 | −0.1018 | −0.0763 | 0.75 |

Tres filas consecutivas con la misma razón. **La curva medida en esa banda es la curva del ruido
del estimador**, no del mercado.

Y la confirmación independiente ya está en el propio reporte: a `1800 s` el control barajado dio
**+0.1292**, mayor en valor absoluto que el real. Un control barajado grande no es una anomalía —
**es la resolución del instrumento, medida**.

---

## §C.1. El error de diseño: el control barajado nunca se leyó como lo que es

El §1 pidió el control barajado como comprobación cualitativa («¿es el real mayor que el
barajado?»). Es más que eso: **el barajado es el suelo de ruido del instrumento**, y con él se
puede decidir *antes de mirar el dato real* si el experimento puede responder la pregunta.

Dos defectos, los dos baratos de corregir:

1. **Un solo barajado por horizonte.** Un sorteo es una moneda, no una distribución. Hacen falta
   ≥ 200 para tener `q95(|R²_nulo|)`.
2. **No se comparó con el requisito.** Es la comparación que decide si la fila es informativa.

Simulado bajo el nulo con `k = 4` y los `n_train`/`n_val` reales de la ejecución (4 000 sorteos):

| `H` | `n_train` | `n_val` | **`q95(\|R²_nulo\|)` = resolución** | `R²_req` (borde favorable) | **razón req/resolución** | ¿informativa? |
|---|---|---|---|---|---|---|
| 60 s | 590 | 112 | **0.060** | 82.7 % | **13.8** | **sí** |
| 120 s | 295 | 56 | **0.127** | 38.5 % | **3.0** | sí, al límite |
| 300 s | 118 | 23 | **0.339** | 14.0 % | 0.41 | **no** |
| 600 s | 59 | 11 | **0.879** | 6.5 % | 0.074 | **no** |
| 900 s | 39 | 8 | **1.641** | 4.2 % | 0.025 | **no** |
| 1800 s | 20 | 4 | **8.006** | 1.9 % | 0.0024 | **no** |

*(Borde favorable = `H_p = 0.552` con `c(u) = 26.03`, la esquina más benévola de la banda. Es la
que hay que resolver para refutar la banda entera.)*

**A 30 minutos el instrumento tiene una resolución 400 veces más gruesa que la cantidad que debe
medir.** Cualquier número que devuelva ahí — positivo o negativo — es ruido.

⚠ **Regla que se incorpora al proyecto:** un control barajado se sortea ≥ 200 veces, se reporta su
`q95`, y **ninguna fila se lee como resultado si `R²_req < 3 · q95(|R²_nulo|)`**. Se reporta como
«el instrumento no resuelve», que es distinto de «no hay señal» y distinto de «faltan datos».

---

## §C.2. Y esperar a `captura_estacional` NO lo arregla. Es estructural

La resolución del método de ventanas escala como `q95 ≈ k / n_val`, con
`n_val = T_validación / H`. Despejando la cobertura necesaria para que
`R²_req ≥ 3 · q95`, con la partición 3:1 de la ejecución:

| `H` | `R²_req` favorable | `n_val` necesario (`k = 4`) | **cobertura de validación** | **cobertura total (con entrenamiento)** |
|---|---|---|---|---|
| 300 s | 14.0 % | 184 | 15 h | 61 h |
| 600 s | 6.5 % | 396 | 66 h | **11 días** |
| 900 s | 4.2 % | 620 | 155 h | **26 días** |
| 1800 s | 1.9 % | 1 333 | 667 h | **111 días** |
| 3600 s | 0.90 % | 2 866 | 2 866 h | **478 días** |

Reduciendo a un solo regresor (`k = 1`) se gana un factor 1.4 en `n`. No cambia el orden de
magnitud.

**`captura_estacional` son 21 días y su banda objetivo empieza en 15 minutos.** Aunque la
alimentación se arreglara hoy y el resto corriera sin un solo hueco, **no alcanza para `H = 15 min`
y no llega ni de lejos a `H = 30 min`.** El §11 de la v4.1 presenta la captura como «el único
activo que puede dar potencia al §1»; con el instrumento actual **no puede**, y eso no es una
carencia de datos: es una propiedad del método de trocear en ventanas no solapadas.

---

## §C.3. El instrumento nuevo: medir la covarianza, no trocear en ventanas

### C.3.1 La identidad

Para un predictor lineal construido sobre el flujo firmado `ε_t`, el `R²` a horizonte `H` es una
razón de momentos de segundo orden:

```
R²(H)  =  Corr( ε_t ,  p_{t+H} − p_t )²

        =        R_cum(H)²
           ─────────────────────
            Var(ε) · σ²_r(H)
```

con `R_cum(H) ≡ Cov(ε_t, p_{t+H} − p_t)`.

**Las dos piezas ya están medidas en este proyecto.** `R_cum(τ)` es **exactamente** la función de
respuesta que la v3.1 §2 definió como `R(τ) = E[(p_{t+τ} − p_t)·ε_t]` y midió en 799 rezagos, con
`R(τ) > 0` en 799/799 y factor 150–700 sobre el control barajado. `σ²_r(H)` es la firma de
volatilidad del §2 de la v4.1. **No hay estadístico nuevo**, como manda el §7: hay una razón entre
dos cantidades ya medidas y ya controladas.

### C.3.2 Por qué gana tanto

El método de ventanas usa **un dato por cada `H` segundos**. La covarianza usa **cada tick como
origen de una pareja `(ε_t, p_{t+H} − p_t)`**. Las parejas se solapan, así que no aportan
información independiente en proporción 1:1 — pero el factor `ε_t` decorrelaciona el sumando, y la
varianza del estimador escala como `H^(1−γ)/N` y **no** como `H/N`. Con `γ ≈ 0.5` la penalización
por horizonte es `√H`, no `H`.

Simulado bajo el nulo (precio de paseo aleatorio, signos con memoria larga independientes del
precio) sobre el tramo continuo que **ya está en disco** — `captura_estacional`, 1 829 242 ticks,
23.33 h:

| `H` | ventanas no solapadas equivalentes | **resolución del método de ventanas** | **resolución de la identidad** | `R²_req` favorable | razón req/resolución |
|---|---|---|---|---|---|
| 900 s | 93 | ~0.19 | **0.00016** | 4.16 % | **260** |
| 1800 s | 46 | ~0.75 | **0.00024** | 1.94 % | **81** |
| 3600 s | 23 | ~3.0 | **0.00033** | 0.90 % | **27** |
| 7200 s | 11 | ~12 | **0.00038** | 0.42 % | **11** |

**Cuatro órdenes de magnitud de resolución, sobre datos que ya existen.** Las cuatro filas superan
el margen de 3× exigido en el §C.1.

⚠ **Ese suelo simulado es una cota inferior y hay que medirlo sobre el dato real.** El nulo
simulado supone incrementos de precio iid; el dato real tiene agrupamiento de volatilidad y
superdifusión, que lo inflan. **El suelo operativo se mide con el control barajado sobre la serie
real, ≥ 200 sorteos, por horizonte.** Aunque saliera 100× peor que el simulado, quedaría margen
hasta `H = 1 h`.

### C.3.3 Demostración directa de por qué el método viejo no puede funcionar

Sobre serie sintética con señal **inyectada de magnitud creciente y conocida**, `N = 400 000`,
`H = 8 000` (49 ventanas no solapadas — más de las que la ejecución tuvo a 900 s):

| señal inyectada `a` | **identidad** | **ventanas no solapadas** |
|---|---|---|
| 0.00 | 0.00004 | 0.01335 |
| 0.05 | 0.00010 | 0.06639 |
| 0.15 | 0.00170 | 0.00582 |
| 0.40 | 0.00281 | 0.00210 |

**La identidad es monótona en la señal. Las ventanas no solapadas no tienen relación con ella** —
el valor más alto de toda la columna sale con señal nula. Es el mismo experimento que la ejecución
corrió sobre dato real, con la verdad conocida.

### C.3.4 Lo que la identidad NO da, declarado

- **Es una cota inferior.** Usa un solo regresor (`ε_t`). Un predictor más rico rinde más. Si la
  cota inferior cruza la curva requerida, es un resultado positivo; si no cruza, **no** demuestra
  que un predictor rico tampoco cruce. Para eso está C.3.5.
- **Es una estimación por sustitución de momentos, no un `R²` fuera de muestra.** No hay
  parámetros ajustados (es una correlación), así que la diferencia dentro/fuera de muestra es
  pequeña — pero no nula, y el punto 3 del §C.4 la controla.
- **No es robusta a la no estacionariedad.** Este proyecto ya documenta que `β` cambia de signo
  entre tramos, que `H_p` recorre [0.371, 0.552] y que la antipersistencia no replicó. Por eso el
  punto 2 del §C.4 es obligatorio.

### C.3.5 La versión multivariante, si hace falta

El `R²` del mejor predictor lineal sobre `k` agregados de flujo pasado
`X_j = Σ_{i<L_j} ε_{t−i}` también es una razón de momentos de segundo orden:

```
R²(H) = c' Σ⁻¹ c / σ²_r(H)        c_j = Cov(X_j, r_{t→t+H})   Σ_jk = Cov(X_j, X_k)
```

`c` sale del mismo correlograma cruzado y `Σ` de la autocovarianza del flujo — **las dos medibles a
resolución de tick sobre la muestra entera**. Se corre **sólo si la versión univariante queda
cerca del cruce**, y con corrección de muestra finita por `k` declarada antes.

---

## §C.4. Procedimiento

1. **Suelo de ruido primero, señal después.** Para cada `H ∈ {60, 120, 300, 600, 900, 1800, 3600,
   7200} s`: ≥ 200 barajados de `ε` (preservando el precio), `q95(|R²_nulo|)`, y la marca
   informativa/no informativa contra `R²_req/3`. **Esta tabla se produce y se congela antes de
   calcular un solo `R²` real.**
2. **Por fragmento, no agregado.** `R_cum(H)` y `σ²_r(H)` se calculan **en cada tramo continuo por
   separado** y se reporta la dispersión entre tramos junto a la media. **Si el signo de
   `R_cum(H)` cambia entre tramos, ése es el resultado** y ninguna media lo arregla — es la misma
   inestabilidad que mató a `β` y a `H_p`, y aquí se ve por 20 USD de cómputo.
3. **Calibración contra el instrumento viejo donde el viejo sí ve.** A `H = 60 s` y `H = 120 s` la
   medición directa por ventanas tiene resolución (§C.1). **La identidad debe reproducir allí lo
   que la ejecución ya midió.** Si no lo reproduce, la identidad está mal aplicada a esta serie y
   se abandona sin gastar nada más. Es calibrar el termómetro nuevo contra el viejo en el rango
   donde el viejo funciona.
4. **`N_eff` por estadístico.** Según la regla del §3.2 de la v4.1: se mide el escalado de **este**
   estadístico, no se importa el de la serie de signos. Reportar `N` y `N_eff` en toda fila con IC.
5. **Tabla final** con, por `H`: `R²_identidad` (media entre fragmentos ± dispersión),
   `q95(|R²_nulo|)`, `R²_req` en las **cuatro esquinas** de la banda del §C.5, `n_fragmentos`,
   `N_eff`, y la marca de informatividad.

---

## §C.5. La curva requerida también es una banda en COSTE, no sólo en `H_p`

La ejecución usó `lastre = c(u) + selección adversa = 26.03 + 5.78 = 31.81 USD/BTC = 4.888 pb`,
con `c(u)` la comisión maker+maker pura. **Pero la v4.0 midió `c(u) = 37.65 USD/BTC` con salida
forzada**, que es el coste realista de una estrategia que no puede garantizar salir como maker.

```
lastre optimista  =  4.888 pb   (c(u) = 26.03, todo maker)
lastre realista   =  6.674 pb   (c(u) = 37.65, salida forzada)
factor en R2      =  1.86
```

La curva requerida es por tanto un rectángulo de dos ejes —`H_p ∈ [0.371, 0.552]` y
`c(u) ∈ [26.03, 37.65]`— y hay que reportar **las dos esquinas**:

| `H` | esquina favorable (`H_p = 0.552`, `c = 26.03`) | esquina adversa (`H_p = 0.371`, `c = 37.65`) |
|---|---|---|
| 300 s | 13.99 % | 38.80 % |
| 600 s | 6.51 % | 23.20 % |
| 900 s | 4.16 % | 17.17 % |
| 1800 s | 1.94 % | 10.27 % |
| 3600 s | 0.90 % | 6.14 % |
| 7200 s | 0.42 % | 3.67 % |
| 14400 s | 0.20 % | 2.19 % |

⚠ **`C_respaldo` (§5 de la v4.1) sigue sin medir y va dentro del lastre.** Mientras no exista, la
esquina adversa está **subestimada** y hay que escribirlo en la tabla, no en una nota.

---

## §C.6. Regla de lectura — sustituye al §1.4

Cuatro desenlaces, no tres. Se aplica fila por fila.

| condición | lectura |
|---|---|
| `R²_req < 3·q95(\|R²_nulo\|)` | **el instrumento no resuelve.** No se lee el valor. Se reporta la resolución y la cobertura que haría falta |
| Resuelve, y `R²_identidad` supera la **esquina adversa**, con signo estable entre fragmentos | **hay banda viable.** Se preregistra el contraste económico completo sobre esa banda |
| Resuelve, y `R²_identidad` queda por debajo de la **esquina favorable** | **no hay banda viable a ese `H`** con predictor lineal sobre flujo. Es refutación, y se reporta como tal |
| Resuelve, y cae **entre las dos esquinas**, o el signo cambia entre fragmentos | **no decidible.** Falta cerrar `C_respaldo`, `H_p` o la estacionariedad — se nombra cuál |

⚠ **Nunca se elige la esquina que hace cruzar.** Si sólo cruza contra la favorable, el desenlace es
«no decidible», no «hay banda viable». Es la misma guarda del §2.2 de la v4.1, extendida al coste.

---

## §C.7. Consecuencia sobre el §11: los apagones dejan de ser catastróficos

El §11 de la v4.1 escribe `H_max ≈ T_continuo / 30` y concluye que un hueco parte la captura. **Eso
es cierto para el método de ventanas y falso para la identidad.**

- Método de ventanas: un fragmento de longitud `ℓ` aporta `floor(ℓ/H)` ventanas. Un fragmento de
  `2H` aporta **2**.
- Identidad: un fragmento de longitud `ℓ` aporta `ℓ − H` parejas. Un fragmento de `2H` aporta
  **`H` parejas**.

**Un apagón cuesta `H` segundos de parejas en cada costura, no el tramo entero.** La cantidad a
vigilar deja de ser el tramo continuo más largo y pasa a ser:

```
N_parejas(H) = SUMA sobre fragmentos de max(0, longitud_i - H)
```

Con eso, `captura_estacional` es útil **desde ya** en la banda de 15–30 min aunque siga
fragmentada, y su tramo de 23.33 h ya resuelve hasta `H = 2 h` según §C.3.2.

⚠ **Esto NO retira la instrucción de arreglar la alimentación.** El límite pasa a ser otro y es
peor: **la reproducibilidad entre fragmentos**, que es donde este proyecto ha perdido `β`, `H_p` y
la antipersistencia. Necesitas fragmentos **independientes y en horas del día distintas**, y eso
sólo lo da la captura corriendo. El portátil sigue teniendo que estar enchufado — pero ahora un
apagón cuesta minutos de dato, no tres semanas de espera.

---

## §C.8. Modos de fallo

| # | fallo | firma | detección |
|---|---|---|---|
| 1 | Suelo de ruido estimado con un solo barajado | control barajado mayor que el real y leído como anomalía | ≥ 200 sorteos, `q95` reportada |
| 2 | Fila no resoluble leída como refutación | `R²` negativo creciente en magnitud con `H` | §C.6 fila 1, marca obligatoria |
| 3 | Identidad aplicada sin calibrar | resuelve donde el método viejo veía, y discrepa | §C.4 punto 3, compuerta de abandono |
| 4 | Agregación entre fragmentos que oculta cambio de signo | media limpia sobre tramos que se cancelan | §C.4 punto 2, dispersión obligatoria |
| 5 | Elegir la esquina favorable de la banda | cruce que sólo existe con `c(u) = 26.03` | §C.6, las cuatro esquinas siempre |
| 6 | `N_eff` importado de la serie de signos | IC mal por órdenes de magnitud | §C.4 punto 4 |
| 7 | Versión multivariante corrida antes que la univariante | `k` parámetros donde 1 bastaba para refutar | §C.3.5, condicionada |
| 8 | Superdifusión y agrupamiento de volatilidad inflando el nulo real | suelo simulado optimista usado como operativo | §C.3.2, el suelo se mide sobre la serie real |

---

## §C.9. Criterios de aceptación

- [ ] Tabla de resolución (`q95(|R²_nulo|)`, ≥ 200 sorteos) producida y congelada **antes** de
      calcular `R²` real.
- [ ] Identidad calibrada contra la medición directa a 60 s y 120 s; discrepancia reportada.
- [ ] `R_cum(H)` y `σ²_r(H)` por fragmento, con dispersión y con el signo de cada fragmento.
- [ ] Las cuatro esquinas de la banda requerida en toda tabla; `C_respaldo` marcada como faltante.
- [ ] `N` y `N_eff` de **este** estadístico en toda fila con IC.
- [ ] Regla del §C.6 aplicada literalmente y el desenlace nombrado por fila.
- [ ] `N_parejas(H)` reportada para `captura_estacional` en su estado actual fragmentado.
- [ ] Verdicto de 60 s y 120 s de la ejecución anterior **conservado sin cambios**.
- [ ] `Micelio.py` sin cambios. Conjunto de prueba de la v3.2 no abierto. Suites en verde.

---

## §C.10. Prioridad

1. **§C.4 puntos 1 y 3** — suelo de ruido y calibración. Sin ellos nada de lo demás se lee. Horas.
2. **§C.4 punto 2** — la identidad por fragmento sobre `captura_v33` y sobre el tramo de 23.33 h de
   `captura_estacional`. **Ésta es la medición que decide el proyecto.**
3. **§C.5** — la banda de coste, aritmética pura.
4. **§C.3.5** — multivariante, sólo si la univariante queda cerca del cruce.
5. **§C.7** — `N_parejas(H)` sobre el estacional fragmentado, para saber qué hay ya.
