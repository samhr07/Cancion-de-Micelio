# Preregistro — Cribado de activos (v4.0)

**Fecha de redacción:** 2026-08-29. **Autor:** Samuel Hoyos R.
**Commiteado ANTES de cualquier medición**, como en `PREREGISTRO_3_1.md` (`66bed94`).
Proyecto anterior cerrado en `ACTA_CIERRE_MICELIO_1_0_A_3_1.md`.

> Regla que gobierna este documento: **ningún umbral declarado aquí se rebaja después.**
> Si el cribado no pasa la compuerta, el resultado del cribado es "no pasa", no un umbral
> nuevo. Si un número de este documento resulta estar mal, se corrige **con un commit
> fechado que diga qué cambió y por qué**, nunca en silencio.

---

## 1. Hipótesis, declarada antes de medir

**Condicionar por volatilidad DENTRO de un activo no mejora el margen** — medido, y
predicho antes de medirlo (§1 estratificado, 2026-08-28). **Elegir un activo
permanentemente más volátil sí puede hacerlo**, porque el mecanismo es distinto:

- Los tramos volátiles de BTC son **tramos de información**, donde la competencia por
  capturar el movimiento es máxima. Más σ y más competencia se cancelan.
- Un activo de menor capitalización es **estructuralmente** más volátil y **menos
  competido**. Ahí σ sube sin que suba la competencia.

La hipótesis es sobre el **mecanismo**, no sobre σ: predice que la relación entre σ y
capturabilidad tiene **signo distinto** entre estratos (dentro de un activo) y entre
activos.

## 2. Predicción falsable, declarada antes de medir

> **En el cribado, σ₁ y R² predictivo NO estarán negativamente correlacionados entre
> activos**, a diferencia de lo que ocurre entre estratos dentro de BTC.

Operacionalizada, para que no admita reinterpretación posterior:

- Estadístico: **ρ de Spearman** entre `sigma1_roll` y `r2_pred_oos` sobre los pares del
  cribado que superen el filtro de datos (§6.1), con su p-valor bilateral.
- **La hipótesis MUERE en la etapa 1** si `ρ_Spearman ≤ −0.30` con `p < 0.05`.
  En ese caso **no se captura nada** y el proyecto se cierra con su acta.
- **La hipótesis SOBREVIVE** (no se confirma: sobrevive) si `ρ_Spearman > −0.30` o
  si `p ≥ 0.05`. Sobrevivir no autoriza la etapa 2 por sí solo: hace falta además que
  algún par pase la **compuerta 1** del §4.
- Se declara de antemano que con ~30 pares el poder es bajo: `p < 0.05` bilateral en
  Spearman con n = 30 exige |ρ| ≳ 0.36. **Un "no rechaza" aquí significa "no se detectó",
  no "no lo hay"**, y así se reportará.

## 3. Criterio de diseño, derivado y fijo

### 3.1 La identidad económica

Predictor `ŷ` del retorno a horizonte `H`, con `R²` **fuera de muestra** contra el
retorno realizado `y_H`. Correlación `ρ = √R²`. Estrategia de signo, entrada y salida
**maker** (el §6.3 explica por qué maker y no taker).

Bajo normalidad bivariante, el margen bruto esperado por ida y vuelta es

```
E[y_H · signo(ŷ)] = kappa · rho · sigma_H ,      kappa = sqrt(2/pi) = 0.79788
```

y el coste es `2·c_lado`. El equilibrio da, con `sigma_H = sigma_1 · H^H_p`:

```
R2_req(H) = ( 2 * c_lado / ( kappa * sigma_1 * H**H_p ) )**2
c_lado_equilibrio = kappa * sqrt(R2_medido) * sigma_1 * H**H_p / 2
```

**`R²_req` es lineal en `c_lado²` e inversa en `σ₁²`.** Como la comisión de futuros
USDⓈ-M es **idéntica en todos los pares** (VIP 0: maker 0.0200 %, taker 0.0500 %), la
única variable libre entre pares es `σ₁`. Todo el cribado se ordena por ese número.

### 3.2 Las cifras que vienen de BTC

| magnitud | valor | origen |
|---|---|---|
| comisión de equilibrio con el `R²` medido en BTC | **0.254 pb por lado** | §1 estratificado, 2026-08-28 |
| comisión real (maker, VIP 0) | **2.000 pb por lado** | tarifa pública USDⓈ-M |
| razón, = factor necesario sobre `σ₁` de BTC agregado | **7.87×** | 2.000 / 0.254 |
| factor necesario sobre `estacional_0` | **3.4×** | §1 estratificado |
| factor necesario sobre `estacional_1` | **1.2×** | §1 estratificado |
| `R²` del modelo direccional | **0.002 – 0.008** fuera de muestra | proyecto BTC |
| `R²` del modelo de volatilidad | **0.615** fuera de muestra | proyecto BTC |

⚠ **Estas seis cifras se toman como DADAS y no se remiden aquí.** Provienen de sesiones
(2026-08-28, 2026-08-29) cuyo código y datos **no están en este repositorio** (§7.2).

### 3.3 Comprobación de consistencia interna, declarada como test

De la identidad del §3.1, el `σ₁` que hace `R²_req = R²_medido` a `H = 3600 s` con
`c_lado = 2.000 pb` y `H_p = 0.5` es **`σ₁ = 1.30 pb·s^(−0.5)`**, e implica

```
R2_referencia = ( 2*2.000 / (0.79788 * 1.30 * 3600**0.5) )**2 = 0.004131
```

**Predicción de consistencia:** `0.004131` cae dentro de la banda declarada `0.002–0.008`
del modelo direccional de BTC. Se fija como test permanente
(`test_v40_umbral_sigma1_es_consistente_con_R2_declarado`). **Si fallara, el que está mal
es este preregistro**, y se corrige aquí con commit fechado — no se ajusta el umbral para
que cuadre.

### 3.4 ⚠ Una inconsistencia entre las cifras dadas, registrada y no resuelta

Invirtiendo la identidad del §3.1 con `H_p = 0.5`, los tres factores del §3.2 implican:

| referencia | factor declarado | `σ₁` implicado [pb·s^(−0.5)] |
|---|---|---|
| agregado (2.000 / 0.254) | 7.87× | **0.165** |
| `estacional_0` | 3.4× | 0.382 |
| `estacional_1` | 1.2× | 1.083 |

El agregado sale **por debajo de los dos estratos**, cuando debería quedar entre ellos.
La explicación más probable es que la comisión de equilibrio de 0.254 pb se calculó con el
`R²` del predictor **agrupado**, que es menor que el de cada estrato ajustado por separado
—justo lo que hace que condicionar por volatilidad no pague—, mientras que 3.4× y 1.2×
usan el `R²` de su propio estrato. Con `R²` distintos por fila, los tres números son
compatibles y la tabla no es contradictoria.

**No se resuelve aquí y no se usa ninguno de los tres como umbral.** La compuerta C1.1 se
apoya únicamente en `σ₁ ≥ 2.60 pb·s^(−H_p)`, que es absoluto y no depende de qué `R²` de
BTC se tome como referencia. Los factores 3.4× y 1.2× se reportan como contexto.
Queda como **primera pregunta a verificar** cuando el material de las sesiones del §7.2
esté disponible.

## 4. Compuerta 1 — umbrales declarados ANTES de correr el cribado

Ningún par pasa a la etapa 2 si incumple **cualquiera** de estos. Se evalúan en este orden
y se reporta **la primera** que falla, junto con todas.

| # | criterio | umbral | razón |
|---|---|---|---|
| **C1.1** | `σ₁` corregido por Roll | **≥ 2.60 pb·s^(−H_p)** | 1.30 (equilibrio exacto) **con margen de 2×** |
| **C1.2** | horquilla relativa mediana | **≤ 1.00 pb** | por encima, el spread se come el margen antes de empezar |
| **C1.3** | granularidad `minQty · precio` | **≤ 5 USD** y ≥ 40 lotes en el nocional máximo de orden | Sec. A.3 de la v1.3: por debajo, el controlador continuo se vuelve interruptor |
| **C1.4** | fracción del movimiento en el 1 % de los minutos | **≤ 0.35** | volatilidad que es toda salto no se captura con órdenes maker |
| **C1.5** | curtosis de retornos de 1 min | **≤ 60** | ídem, por otra vía; se reporta siempre, no solo cuando falla |
| **C1.6** | cobertura de datos | ≥ 95 % de los minutos esperados en 90 días | un par con huecos no es un par comparable |

**Lectura de C1.1 declarada aquí para que no se reinterprete después:** "1.30 con margen de
2×" significa **`σ₁ ≥ 2.60`**, no `σ₁ ≥ 1.30`. El margen es un factor sobre el umbral, no
una tolerancia sobre el fallo.

**C1.1 se evalúa sobre `σ₁` CORREGIDO POR ROLL, no sobre el crudo.** σ de velas de 1 minuto
**sobreestima** por el rebote bid-ask, y sobreestima **más** cuanto más ancha es la
horquilla — que es exactamente la dirección que favorece a los pares malos. Se reportan las
dos y su diferencia; **la compuerta usa la corregida**.

## 5. Regla de parada, acordada antes

- **Etapa 1** termina cuando la tabla completa de ~30 pares está medida y reportada.
  Se reporta **la tabla entera**, ordenada por `σ₁`, no los tres primeros.
- **Etapa 2** se corre sobre **UN** par: el de mayor `σ₁` corregido entre los que pasan la
  compuerta 1. **No se cambia de par a mitad de la etapa 2 porque otro "se ve mejor".**
  Si ese par falla los cuatro invariantes de tubería, se para y se diagnostica la captura
  — no se prueba con el siguiente.
- **Etapa 3** para, con veredicto negativo, si `R²_max` de la identidad multivariante **no
  supera el suelo por rotación circular** en el par nuevo. Para con veredicto positivo si
  `R(H, θ)` tiene un máximo interior con `R > 0` **usando la horquilla real del par**.
- **El proyecto entero para** si la etapa 1 mata la predicción del §2.

## 6. Decisiones metodológicas fijadas antes de medir

### 6.1 Filtro de datos, antes de cualquier estadístico
Se descartan velas con `volume = 0` o `numTrades = 0`, y pares con cobertura < 95 %.
⚠ El feed emitió transacciones con `p = 0` en ~0.2 % de los casos (v3.0); en klines el
análogo es la vela vacía, y **enmascara la cola en vez de crearla** (la curtosis SUBE al
limpiar). Se limpia **antes** de medir varianza, no después.

### 6.2 Control barajado obligatorio
`σ(H)` y `H_p` se ajustan también sobre los **retornos barajados** del mismo par. Bajo
barajado los retornos son iid y `H_p` debe salir **0.5**. Se declara el criterio:
**si `|H_p_barajado − 0.5| > 0.05`, el estimador tiene sesgo propio y el `H_p` real de ese
par no se reporta como medida.** Es el control del A.3.4, que ya salvó dos análisis.

### 6.3 Maker, no taker
El proyecto viejo midió que con comisiones taker `σ` **no alcanza el coste dentro de los
300 s** (v3.1 §1). Todo el criterio económico se evalúa **maker+maker**. Por eso importa
C1.4: el decil superior del movimiento **no se puede capturar con órdenes maker**, luego un
activo cuya volatilidad vive en saltos discretos no sirve aunque su `σ₁` pase.

### 6.4 `R²` predictivo del cribado es un PROXY y se declara como tal
De klines de 1 minuto se construye un predictor con la misma familia de rasgos que el
proyecto BTC (desequilibrio de flujo desde `takerBuyBaseVolume`, retornos rezagados, rango,
`numTrades`), ajustado **fuera de muestra** por partición temporal 70/30 sin solapamiento.
**No es el modelo de BTC** y no se compara con su `R²` en valor absoluto: se usa **solo**
para la correlación entre pares del §2, que es un estadístico de **orden**.

## 7. Criterios de aceptación que se declaran NO CUMPLIDOS de antemano

Se declaran ahora para que no se presenten después como resultados.

### 7.1 El cribado NO puede ejecutarse desde esta sesión
`fapi.binance.com`, `api.binance.com` y `data-api.binance.vision` están **denegados por la
política de egreso** de esta sesión remota (403 al CONNECT; GitHub sí resuelve, luego es
política y no red). **Etapa 1 queda implementada y probada contra datos sintéticos con
verdad conocida, y sin ejecutar contra mercado.** La corrida real es local, en un solo
comando.

### 7.2 Las cifras de BTC del §3.2 no son reproducibles en este repositorio
El código y los datos de las sesiones 2026-08-11 a 2026-08-29 —el modelo de volatilidad
`R² = 0.615`, `captura_estacional.py`, el §1 estratificado, el §C.3.5, la corrección
falsable/no-falsable del 2026-08-29— **no están en este árbol**. `git log` termina en
`d591256` (v3.1 §2). Las seis cifras del §3.2 entran como **constantes declaradas**, no
como mediciones replicadas, y el código las exige como parámetros explícitos con
**centinela en cero**: olvidarlas falla ruidosamente en vez de operar con un literal.

### 7.3 El escalón de comisiones sigue asumido, no leído
`/fapi/v1/commissionRate` es firmado. Se usan tarifas públicas VIP 0. Mismo criterio no
cumplido que en la v3.1, arrastrado y declarado.

## 8. Lo que NO se hace

- **No se toca `Micelio.py`.**
- **No se refina el predictor** antes de saber si algún activo cruza.
- **No se rebaja ningún umbral** declarado en el §4.
- **No se cambia de par** a mitad de la etapa 2.
- No se adopta ninguna hipótesis abierta del proyecto viejo (§6 del acta de cierre).

---

**Firmado por commit.** Cualquier medición fechada después de este commit se contrasta
contra los umbrales de arriba tal como están escritos.
