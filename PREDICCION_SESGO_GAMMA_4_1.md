# Predicción de sesgo del tercer estimador de `γ` — v4.1 §3.1

**Escrito ANTES de calcular ningún GPH ni ningún Whittle local.** El §3.1 de
`ORDEN_TRABAJO_HORIZONTE_4_1.md` lo exige con esas palabras: *«un tercer estimador que no
comparta modo de fallo con ninguno de los dos … Y declarar por qué se esperan sesgos … pero es
una predicción, y se verifica, no se asume.»*

Este archivo se congela con su commit. Cualquier cosa que se mida después se contrasta contra lo
que está aquí, no al revés.

---

## 1. El punto de partida

Dos estimadores sobre **la misma serie de signos del mismo tramo** (`captura_v33`,
entrenamiento, 618 692 ticks) dan dos respuestas:

| método | `γ` | qué ajusta |
|---|---|---|
| A — regresión log-log de `C(ℓ)`, rezagos `[10, 2000]` | **0.5217** | la **autocovarianza** en el dominio del tiempo |
| B — escalado de `Var(media de bloque)` contra `b` | **0.3539** | la **varianza agregada**, también en el tiempo |

Ninguno de los dos es espectral. Ése es el hueco que llena el tercero.

## 2. El tercer estimador y la equivalencia de parámetros

Se usan **los dos** estimadores semiparamétricos en frecuencia, porque no cuesta más y sus
modos de fallo tampoco coinciden entre sí:

- **GPH** — regresión log-periodograma (Geweke & Porter-Hudak 1983). Regresa
  `log I(λ_j)` sobre `log|2 sin(λ_j/2)|` en las `m` frecuencias de Fourier más bajas;
  la pendiente es `−2d`.
- **Whittle local** (Robinson 1995). Minimiza la verosimilitud de Whittle localizada en la
  misma banda `[1, m]`. Es asintóticamente más eficiente que GPH y no arrastra el sesgo de
  regresión de éste.

La traducción al vocabulario del proyecto, y hay que escribirla porque en este documento `γ`
significa algo muy concreto:

```
C(l) ~ l^(-gamma)     <=>     f(lambda) ~ lambda^(-2d)  cuando lambda -> 0
                      <=>     gamma = 1 - 2*d
                      <=>     d = (1 - gamma)/2
```

Con eso, las dos cifras en disputa se leen así:

| | `γ` | `d` equivalente |
|---|---|---|
| A (log-log de `C`) | 0.5217 | **0.2392** |
| B (`Var` de bloque) | 0.3539 | **0.3231** |

⚠ **Coincidencia algebraica que NO hay que malinterpretar:** `d = (1−γ)/2` es literalmente la
misma expresión que la `β` de difusividad del §1 de la v3.3. Son cantidades distintas que
comparten fórmula. No es evidencia de nada y no se citará como tal.

Régimen de validez: `d ∈ (0, 0.5)` es memoria larga estacionaria, que es donde caen las dos
cifras. Si algún estimador devuelve `d ≥ 0.5` la serie no es estacionaria en el sentido que
estos estimadores suponen, y eso invalidaría el estimador, no la serie.

## 3. Las predicciones, en orden de compromiso

### P1 — Dirección. **`γ_GPH < 0.5217` y `γ_LW < 0.5217`.**

Mecanismo: la ACF muestral se centra con la **media muestral**, no con la media verdadera. Bajo
memoria larga la media muestral converge lentísimo, así que absorbe parte de la propia
dependencia de largo alcance y deprime `Ĉ(ℓ)` **más cuanto mayor es `ℓ`**. El resultado es un
ajuste log-log **más empinado del que corresponde**, o sea `γ` sobreestimada. El estimador A es
el único de los tres que centra por la media muestral en el dominio del tiempo.

Los estimadores en frecuencia no comparten ese modo de fallo: el periodograma en las frecuencias
de Fourier `λ_j = 2πj/N` con `j ≥ 1` es **exactamente invariante** a la media muestral, porque
`Σ_t e^{−i λ_j t} = 0`. No es que el sesgo sea pequeño: no existe por esa vía.

**Falsación:** si cualquiera de los dos devuelve `γ ≥ 0.5217`, P1 queda refutada y la lectura es
que el sesgo por centrado no domina, con lo que la cifra de la v3.3 queda reivindicada.

### P2 — Magnitud. **`γ_GPH` y `γ_LW` caen en `[0.28, 0.46]`**, o sea más cerca de B que de A.

Es la predicción que la propia orden apunta («0.5217 es el sesgado y 0.3539 el más cercano»)
puesta en un intervalo que se puede fallar. Es la más comprometida de las cuatro y la que menos
mecanismo tiene detrás: el tamaño del sesgo de centrado depende de `d` y de `N` y no lo he
calculado, sólo su signo.

**Falsación:** cualquier valor fuera de `[0.28, 0.46]`.

### P3 — Dependencia de la banda. **`γ̂` CRECE con `m`** en el barrido `m = N^α`, `α ∈ {0.4, 0.5, 0.6, 0.7}`.

Mecanismo, y es el sesgo propio del tercer estimador —que hay que declarar igual que el del
primero—: GPH y Whittle local son válidos sólo en `λ → 0`. Al ensanchar la banda entran
frecuencias donde manda la estructura de **corto** alcance. El flujo de órdenes de este mercado
tiene `C(1) = 0.798`, una correlación positiva de corto alcance enorme, que aplana el espectro
lejos del origen; contaminar con esas frecuencias **reduce `d̂`**, y por `γ = 1 − 2d` eso
**aumenta `γ̂`**.

**Falsación:** que `γ̂` decrezca con `m`, o que no se mueva de forma monótona.

**Y es la predicción con más consecuencia práctica:** si `γ̂` depende fuerte de `m`, ninguno de
los tres números es «la» `γ`, y la discrepancia de §3.1 deja de ser un problema de estimador
para pasar a ser un problema de **ausencia de un régimen de escala único** — que es exactamente
lo que el §3.3 pregunta por la otra vía. Las dos secciones estarían midiendo la misma patología.

### P4 — Control de verdad conocida. **Sobre fGn con `d` conocida, ambos recuperan `d` con error < 0.05** en `m = N^0.5`.

Se corre **antes** de mirar el dato real. Si un estimador no recupera una verdad conocida, su
número sobre BTC no entra al reporte. Es el mismo criterio con el que en la v3.3 se tiró a la
basura mi generador de suma de AR(1) y se sustituyó por fGn.

## 4. Lo que este ejercicio NO puede decidir

Aunque los tres estimadores converjan, **`γ` no mide el núcleo del propagador**. Entra al
proyecto sólo a través de `β_implícita = (2−γ)/2 − H`, y el propio §3.1 establece que esa
cantidad es un **contraste de coherencia entre dos exponentes**, no una tercera medición de la
permanencia del impacto. Convergencia aquí no reabre permanente-contra-transitorio.

## 5. Registro

| | |
|---|---|
| escrito | 2026-08-12, antes de ejecutar `gph` o `whittle_local` sobre dato real |
| series a las que aplica | `captura_v33` entrenamiento (la que la v3.3 midió) y el tramo continuo de 23.33 h de `captura_estacional` como réplica independiente |
| `Micelio.py` | sin cambios |
