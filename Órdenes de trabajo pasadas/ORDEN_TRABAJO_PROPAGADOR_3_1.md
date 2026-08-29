# Orden de Trabajo — Horizonte Derivado y Propagador con Forzamiento Medido v3.1

Sucede a `ORDEN_TRABAJO_OSCILADOR_3_0.md`. La v3.0 dio veredicto **`k = 0`**: sobre 446 892
ticks, `φ₁ + φ₂ = 1` idénticamente, y el AR(2) queda agotado a siete decimales por paseo
aleatorio más rebote bid-ask (`a = ρ₁(retornos) = −0.216061`). Sin fuerza recuperadora, sin
oscilador libre, y Harvey sin parametrización válida (`φ₂ > 0`).

**Lo que este documento cambia — y es un cambio de método, no de modelo:**

1. **El horizonte de operación deja de ser una elección de diseño y pasa a derivarse de los
   costes.** Es un número medido, y **móvil**: cambia con el régimen de volatilidad.
2. **El forzamiento deja de ser una hipótesis y pasa a ser un dato.** `es_maker` y `q` llevan
   capturándose desde la v2.0 sin usarse. Con `F` observado, el test de la v3.0 se rehace sin el
   sesgo que lo lastraba.

**Lo que NO cambia:** `Micelio.py` no se toca. Todo es código de análisis aparte, como en la
v2.2 y la v3.0.

---

## §0. Registro de errores previos

| Origen | Error | Estado |
|---|---|---|
| v3.0 §3 (documento) | El AR(2) mete `F` **dentro del residuo** `ε`. El flujo de órdenes tiene **memoria larga** (ley de potencias, por fragmentación de órdenes), así que es un regresor omitido y autocorrelado: sesga y confunde los coeficientes | El veredicto `k = 0` sigue en pie para la parte **no forzada**. Pero no es un test de la premisa de la v3.0, que era un oscilador **forzado**. §2 lo rehace |
| Diseño (todas las versiones) | Horizonte de operación fijado por diseño ("entradas rápidas"), nunca derivado | §1. Con las comisiones medidas, la escala de tick es económicamente irrelevante |
| v3.0 §3.3 (documento) | Guarda especificada sobre `γ < 0`; el fallo real fue `m < 0`, y `γ` fue positiva en **38/38** bloques | Guarda corregida por Code (`oscilador.guarda_masa`). Patrón a recordar: se vigiló el parámetro equivocado |
| Producción, v2.0–v3.0 | Trades con `p = 0`, `q = 0` **sin filtrar** (~0.2 %): innovación espuria de −65 000 USD cada ~30 s | Corregido en `mercado.tick_valido`. ⚠ **Toda cifra de varianza previa está mal**: la curtosis SUBE de 601.8 a 1179.7 al limpiar, y el incremento máximo cae de 65 245 a 11.8 USD (**factor 5 500**) |

⚠ **Consecuencia del último punto que esta orden de trabajo asume:** el nulo espectral de la
v2.2 (multitaper sin exceso sobre ruido rojo, ×1.78) se midió sobre datos sucios. Unos pocos
impulsos de 5 500× la amplitud típica son deltas: levantan el suelo de ruido en **todas** las
frecuencias y entierran cualquier pico moderado. **Ese nulo no está cerrado hasta rehacerse
limpio**, y hasta entonces no se cita como establecido.

---

## §1. El horizonte se deriva de los costes

### 1.1 La medición que lo motiva

Con `tickSize = 0.10 USD` y BTC en ~63 000:

| esquema | coste ida y vuelta | en ticks |
|---|---|---|
| maker + maker (~0.02 % c/u) | 25.2 USD/BTC | **252** |
| maker + taker | 44.1 USD/BTC | 441 |
| taker + taker (~0.05 % c/u) | 63.0 USD/BTC | **630** |

⚠ **Leer el escalón de comisiones real de la cuenta**, no asumir estos porcentajes — mismo
principio que `exchangeInfo` en la v1.3. El orden de magnitud no cambia con el escalón, pero el
número sí.

Horizonte al que BTC se mueve lo que cuesta la ida y vuelta:

| vol anual | maker | taker |
|---|---|---|
| 25 % | 81 s | 505 s |
| 40 % | **32 s** | 197 s |
| 60 % | 14 s | 88 s |
| 90 % | 6 s | 39 s |

**Dos consecuencias que reordenan el proyecto:**

- Una predicción a escala de tick, aunque fuera perfecta, no paga las comisiones. El AR(2) de la
  v3.0 midió a `Δn = 1 tick` con rezagos 1–2: **midió bien, en una escala que no es la de este
  sistema.**
- **La latencia de 300 ms deja de ser una restricción.** Contra un horizonte de 32 s es el
  **0.95 %**. Lleva seis versiones condicionando decisiones de diseño y, al horizonte correcto,
  es despreciable.

### 1.2 La definición

```
H*  resuelve   σ(H*) = c(u)
```

- `σ(H)` — volatilidad realizada al horizonte `H`, medida con **núcleo realizado** o estimador
  de dos escalas (el ruido de microestructura aquí es severo: curtosis 1179.7 en limpio).
- `c(u)` — coste de ida y vuelta: comisiones + spread efectivo cruzado (si taker) + impacto.

⚠ **`c` depende del tamaño de orden `u` a través del impacto**, luego `H* = H*(u)`. Órdenes
mayores necesitan horizontes mayores. Eso acopla `H*` directamente con la variable de decisión
del NMPC y con la capa de Loeper, que ya existe. **No es un parámetro nuevo: es una relación
entre dos que ya están.**

Spread efectivo por el modelo de Roll, que ya está medido a medias: `s_eff = 2σ_r·√(−ρ₁)`, con
`ρ₁ = −0.216061` de la v3.0 y `σ_r` **medida sobre datos limpios** (la sucia no vale, §0).

### 1.3 Lo que `H*` sustituye

`H*` es el "número movible" que el proyecto no tenía. Pasa a fijar, en vez de constantes
elegidas:

- el horizonte de predicción del propagador (§2)
- el horizonte del NMPC
- la escala de muestreo `Δτ` relevante
- el umbral mínimo de ventaja para que una operación tenga sentido

Publicarlo en telemetría como magnitud viva, recalculado con la volatilidad del régimen. Si
`H*` se dispara —volatilidad baja, coste constante— es información operativa: el sistema debe
operar menos, no más.

---

## §2. El propagador con forzamiento medido

### 2.1 El modelo

```
Δp_t = Σ_{k≥0} h(k) · ε_{t−k} · f(v_{t−k})  +  η_t
```

con `G(τ)` el propagador (impacto de una transacción a rezago `τ`), `h(0) = G(0)`,
`h(k) = G(k) − G(k−1)`, `ε` el signo de la transacción y `f(v)` la dependencia con el volumen.

**Por qué esto sí contrasta la premisa de la v3.0:** en el AR(2), `F` estaba en el residuo. Aquí
está medido. Y **la respuesta al impulso `G(τ)` puede ser oscilatoria aunque el AR(2) libre no
lo sea** — un impacto que sobrepasa y revierte *es* un oscilador amortiguado, solo que forzado
por algo observado en vez de por ruido supuesto. Es la premisa del oscilador forzado con la
pieza que faltaba.

También contesta el punto de las "regiones con frecuencia" mejor que una wavelet: si `G(τ)`
tiene raíces complejas en unos regímenes y reales en otros, eso **es** el fenómeno, con
estimador propio y nulo propio.

### 2.2 ⚠ Convención de signo — verificar antes que nada

El campo `m` de Binance indica **si el comprador fue el maker**. Luego:

```
m = True   ->  comprador maker  ->  el taker vendía  ->  ε = −1  (venta agresora)
m = False  ->  comprador taker  ->  ε = +1           (compra agresora)
```

**Confirmarlo contra la documentación vigente y contra los datos** (el impacto inmediato `G(0)`
debe salir **positivo**; si sale negativo, el signo está invertido). Un signo al revés produce
un propagador espejo, plausible y completamente equivocado — la misma familia de fallo que el 2π
y el factor 125, y la tercera vez que este proyecto se juega un resultado en un signo.

### 2.3 Parametrizar el núcleo, no los rezagos

A `ν ≈ 102 tx/s` y `H* ≈ 32 s`, cubrir el horizonte con rezagos libres son ~3 300 coeficientes.
En su lugar, forma paramétrica:

```
G(τ) = G₀ · (1 + τ/τ₀)^(−β)          3 parámetros
f(v) = v^δ                            1 parámetro
```

**Cuatro parámetros estimados con errores estándar, frente a 3 300 libres.** Es la dirección que
pide el proyecto: menos constantes, y las que quedan medidas en vez de elegidas.

Barrer `δ ∈ {0, 0.25, 0.5, 1}` y elegir por verosimilitud fuera de muestra. El impacto
individual es cóncavo (`δ` pequeño) en la literatura; que salga otra cosa sería informativo.

---

## §3. El contraste de difusividad — un nulo derivado de teoría

### 3.1 La relación

Si la autocorrelación de signos decae como `C(τ) ~ τ^(−γ)` y el propagador como
`G(τ) ~ τ^(−β)`, la condición para que el precio sea difusivo (retornos sin autocorrelación) es

```
β = (1 − γ)/2
```

Es el resultado central de la línea Bouchaud–Gefen–Potters–Wyart / Lillo–Farmer: el decaimiento
del impacto compensa la memoria del flujo, y de ahí que el precio sea martingala **pese a que
el flujo sea muy predecible**.

Así que el contraste tiene tres lecturas, fijadas por teoría y no por umbral ajustado:

| resultado | lectura |
|---|---|
| `β ≈ (1−γ)/2` | difusivo — **sin ventaja explotable** |
| `β < (1−γ)/2` | el impacto decae demasiado despacio — **momentum** |
| `β > (1−γ)/2` | el impacto decae demasiado rápido — **reversión** |

**Se medirán `β_k ≠ 0` con enorme significancia. Eso no es la señal: es el impacto, y ya está
descontado.** La señal es la desviación respecto a `(1−γ)/2`.

### 3.2 ⚠ El umbral asintótico NO vale como umbral finito

Simulado con signos de memoria larga y propagador de ley de potencias, buscando el `β` donde
`ρ₁(retornos)` cruza cero:

| γ | β teórico `(1−γ)/2` | β con `ρ₁ = 0` medido | error |
|---|---|---|---|
| 0.30 | 0.350 | 0.446 | **0.096** |
| 0.50 | 0.250 | 0.314 | 0.064 |
| 0.70 | 0.150 | 0.181 | 0.031 |

El signo y el orden son correctos —el cruce se mueve como predice la teoría, y el error
converge al crecer `γ`— pero **hay sesgo sistemático de tamaño finito**, del orden de 0.03–0.10.
Con `K` finito y la transformación `sign()`, la relación asintótica no se alcanza.

**Por tanto: calibrar el umbral por simulación**, con `N`, `K`, `γ` y `f(v)` emparejados a los
datos reales, y comparar `β` medida contra esa distribución simulada — no contra `(1−γ)/2`
directamente. Es el mismo error que los umbrales χ² sobre el NIS (v2.1 §6): usar una asintótica
donde la distribución finita es otra. Que aparezca por segunda vez es motivo para tratarlo como
patrón, no como caso.

---

## §4. Preregistro — criterios fijados ANTES de mirar

Este proyecto ha anulado dos veredictos por emitirlos sobre datos que no los sostenían. Los
criterios se **escriben y se commitean con fecha antes de ejecutar nada**, y el commit se
referencia en el reporte.

### 4.1 Criterio de éxito

Las cuatro, conjuntamente:

1. **`R²` fuera de muestra > 0** en corte temporal, con **banda de embargo** de al menos `H*`
   entre entrenamiento y prueba. Sin embargo hay fuga por solapamiento de horizontes.
2. **`E[Δp | información en t]` al horizonte `H*` supera `c(u)`** con margen, no lo iguala.
   Es la única comparación que decide algo económico.
3. **`β` fuera del intervalo simulado** del §3.2, con su dirección declarada.
4. **Residuo del propagador sin autocorrelación remanente** (Ljung-Box a rezagos largos). Si
   queda estructura, el modelo está incompleto y el `R²` no es de fiar.

### 4.2 Criterio de abandono — declarado por adelantado

Si sobre datos limpios, con forzamiento medido, referencia móvil (§5) y régimen congelado (§6),
al horizonte `H*` derivado:

- `β` cae dentro del intervalo simulado, **y**
- `E[Δp]` no supera `c(u)` en ningún régimen, **y**
- el residuo es blanco

entonces **la hipótesis de estructura explotable a este horizonte queda abandonada**, y el
proyecto pasa a lo que la v3.0 §6.2 lista como superviviente: reloj de ticks, ingesta, capa de
riesgo, Loeper y NMPC como ejecutor de órdenes, sin capa predictiva.

Escribir esto ahora, y no después, es lo que impide que la hipótesis se vuelva infalsable por
retirada. `ω_m` lleva tres refutaciones de rigor creciente y cada una se respondió con "quizá a
otra escala". Ésta es la escala donde está el dinero, y no hay otra a la que retirarse.

---

## §5. Referencia móvil — endógena, sin ventana arbitraria

La v3.0 midió reversión a un **nivel fijo** con intercepto constante. Reversión a una tendencia
móvil aparecería como `k = 0` igualmente.

⚠ **Media móvil rezagada, nunca centrada.** Una MA que incluye el punto actual inyecta reversión
por construcción y daría `k > 0` espurio. Sería la sexta versión del mismo defecto.

**Pero hay una referencia mejor y sin ventana libre:** en el marco del propagador, `G(∞)` es el
**impacto permanente** y `G(τ) − G(∞)` el transitorio. El impacto permanente acumulado **es** la
referencia móvil, y sale del mismo ajuste. Elimina la ventana de la MA como parámetro elegido.

Implementar ambas —MA rezagada a varias escalas como control, permanente del propagador como
principal— y reportar si coinciden.

---

## §6. Régimen congelado por adelantado

"En ciertas regiones hay frecuencia y en otras no" solo es contrastable si las regiones se
definen **por adelantado desde algo observable**. Si se eligen mirando dónde apareció, siempre
aparecerán: la v2.2 midió `C = 0.783` sobre incrementos **barajados**, donde por construcción no
hay nada.

**Procedimiento:** declarar **un** regresor de régimen primario y como mucho dos secundarios,
commitearlos antes de ejecutar, y no cambiarlos después.

Candidatos, por orden de preferencia:

1. **Volatilidad realizada** — ya hace falta para `H*` (§1.2), así que no añade maquinaria.
2. **Magnitud del desbalance de flujo** — ya se captura.
3. **Sesión** (Asia / Europa / América) — exógena, sin riesgo de circularidad.

Profundidad del libro sería mejor teóricamente, pero exige `@depth` que hoy no se ingiere.

---

## §7. Constantes que este documento elimina

| desaparece | por qué |
|---|---|
| horizonte de operación (implícito) | derivado: `H*(u)` de coste y volatilidad medidas |
| `ΔS_ref` | candidato a rederivarse de la escala del impacto transitorio |
| `C_ON`, `C_OFF` | nunca se calibraron; `C` se demostró inflada por artefacto |
| `ciclos_min`, `muestras_por_ciclo_min` | de la guarda de banda, si la cadena EMD se retira |
| `W`, `K` de la EMD | ídem |
| `ω_m,max`, `γ_ω` | si el propagador sale monótono |

**Aparecen:** `G₀`, `τ₀`, `β`, `δ` — cuatro, todos **estimados con errores estándar**, ninguno
elegido. Balance neto claramente favorable, que es el objetivo.

---

## §8. Criterios de aceptación

**§1 — Horizonte**
- [ ] Escalón de comisiones **leído de la cuenta**, no asumido.
- [ ] `σ(H)` con núcleo realizado o dos escalas, sobre **datos limpios**; gráfico de firma de
      volatilidad incluido.
- [ ] `s_eff` de Roll con `σ_r` limpia.
- [ ] `H*(u)` publicado en telemetría como magnitud viva, con su dependencia del tamaño.

**§2 — Propagador**
- [ ] **Test de convención de signo**: `G(0) > 0`. Si sale negativo, parar. Documentar la
      correspondencia del campo `m` con evidencia, no con suposición.
- [ ] Núcleo paramétrico de 4 parámetros con errores estándar por bootstrap por bloques.
- [ ] `δ` barrida y elegida por verosimilitud **fuera de muestra**.
- [ ] **Raíces de `G(τ)`**: ¿complejas (oscilatorio, sobrepasa y revierte) o reales (monótono)?
      Es el contraste directo de la premisa del oscilador forzado.

**§3 — Difusividad**
- [ ] `γ` medida de la autocorrelación de signos, con su ajuste de ley de potencias y bondad.
- [ ] `β` del propagador.
- [ ] **Umbral calibrado por simulación** con `N`, `K`, `γ`, `f(v)` emparejados. Comparar contra
      `(1−γ)/2` directamente es motivo de rechazo del entregable (§3.2).
- [ ] Reproducir la tabla del §3.2 como test permanente, incluido el sesgo de tamaño finito.

**§4 — Preregistro**
- [ ] Criterios de §4.1 y §4.2 **commiteados con fecha antes** de ejecutar; hash referenciado en
      el reporte.
- [ ] Corte temporal con banda de embargo ≥ `H*`, verificada.

**§5 y §6**
- [ ] MA **rezagada** verificada por test (una MA centrada debe hacer fallar el test).
- [ ] Regresor de régimen commiteado antes de mirar.

**Transversal**
- [ ] Suite en verde: **56/56 + los nuevos**.
- [ ] `Micelio.py` sin cambios de modelo.
- [ ] Nulo espectral de la v2.2 **rehecho sobre datos limpios** antes de citarse como cerrado.

---

## §9. Fuera de alcance

- **Implementar el propagador en el filtro.** Este documento mide; si el resultado lo justifica,
  la integración es otra versión.
- **Retirar la Sec. 1.4 ni la Sec. 2 del PDF.** No antes de que §4.1 o §4.2 se resuelvan.
- **Hipótesis (B)** (escala de decenas de minutos) y §4 multiescala de la v2.2: la captura larga
  sufrió un corte de DNS de 10 h y hay que relanzarla. Con `H*` derivado, además, esa escala ya
  no es la prioritaria.
- **Wavelet Torrence-Compo** e **IAAFT sobre incrementos**: confirmatorios, no decisorios.
- **v2.1 §5 (log-precio) y §6 (cuantiles empíricos del NIS)**: siguen en cola, **no bloqueadas**.
  Con la curtosis real en 1179.7 —el doble de lo que se creía— las compuertas χ² están peor
  calibradas de lo que se pensaba, y §6 gana prioridad.
- Fase 2 (ALS), Fase 3, Testnet con credenciales, las 30 corridas.
- CMA-ES, IMPC, DeepONet.
- Reconciliación documental del PDF, que ya exige sesión propia.
