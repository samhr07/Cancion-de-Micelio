# Orden de trabajo — Cribado de activos (v5.1)

**Fecha:** 2026-08-29. Hipótesis, predicción falsable y umbrales en
`PREREGISTRO_CRIBADO_5_1.md`, **commiteado antes que este documento**. Proyecto
anterior sigue registrado en `CLAUDE.md` (llega a la v4.2) y `PLAN_CAPITAL_5_0.md`.

Este documento dice **cómo se ejecuta**. No redefine ningún umbral: los umbrales
están en el preregistro y aquí solo se aplican.

---

## §0 — Tarea residual de BTC: YA ESTÁ HECHA en `master`

⚠ **Esta sección se escribió sobre una base obsoleta (`d591256`) y su premisa era falsa.**
La identidad multivariante del §C.3.5 **no está sin contrastar**: `identidad.py` la
implementa (`r2_identidad_multi`, `R² = c' Σ⁻¹ c / σ²_r(H)` sobre agregados de flujo
pasado `X_j = Σ_{i<L_j} ε_{t−i}`), con suelo por rotación circular, autotest 6/6, y la
sesión 2026-08-27 la ejecutó. **No se duplica.** El módulo que esta rama traía
(`identidad_multivariante.py`) se ha retirado por eso.

Lo que el registro de `master` establece, y que entra aquí como dato de partida:

| | valor |
|---|---|
| `σ(t)` desde 3 rezagos de (`ν`, `σ`) | `R² = +0.6153` contra suelo 0.0798 → **7.7×**, se lee |
| `ν(t) → σ(t+1)` | `R² = +0.5464` contra suelo 0.0508 → **11×**, se lee |
| suelo por rotación contra barajado | **14× más ancho** (0.0798 contra 0.0058) |
| §1 estratificado por `σ` pronosticada | **margen mínimo 9.3×** — LÍNEA CERRADA |
| v4.2 §4, superficie `R(H, θ)` | **NO HAY ÓPTIMO**, `P(max R_nulo ≥ R*) = 0.985` |

⚠ **Y una corrección de fondo a lo que esta rama afirmaba.** `identidad.py` advierte que
«la multivariante no puede rendir menos que la univariante» **NO ES UN TEOREMA**: sólo se
cumple si los rasgos univariantes están **anidados** en los multivariantes. Con agregados
de flujo `X_j` que no contienen `ε_t`, la multivariante **puede rendir menos**. La
redacción anterior de esta rama («superior por construcción») era demasiado fuerte, y el
test que la acompañaba sólo pasaba porque usaba el mismo conjunto de rasgos en los dos
lados, que es anidamiento por construcción.

**Consecuencia sobre el factor 3.4×:** con la multivariante ya medida y las filas
decisivas marcadas NO FALSABLES (`b8f3155`), el factor **no baja**. La relación sigue
siendo `factor_nuevo = factor_actual / √(R²_nuevo/R²_referencia)`, pero no hay ganancia
que aplicar. El cribado se ordena por `σ₁` contra el umbral absoluto de C1.1, que no
depende de ese factor.

## §1 — Etapa 1: cribado sin captura

```
python cribado_activos.py                        # 30 pares, 90 dias
python cribado_activos.py --pares=40 --dias=90 --horquilla=120
python cribado_activos.py --sin-cache            # ignora cribado_cache/
```

### Lo que mide por par

| magnitud | de dónde sale | para qué |
|---|---|---|
| `σ₁`, `H_p` | ajuste log-log de `σ(H)` en **[60, 3600] s** sobre klines de 1 min | **C1.1**, y el orden de la tabla |
| `σ₁` corregido por **Roll** | `s = 2√(−cov₁)`, se resta `s²/2` de la varianza de cada horizonte | **C1.1 se evalúa sobre este** |
| horquilla relativa mediana | `/fapi/v1/ticker/bookTicker`, muestreado N veces | **C1.2** |
| `minQty · precio` y `minNotional` | `/fapi/v1/exchangeInfo` de **Mainnet** | **C1.3** |
| fracción del movimiento en el 1 % de los minutos | `Σ|r|` de la cola / `Σ|r|` total | **C1.4** |
| curtosis de retornos de 1 min | momento cuarto estandarizado | **C1.5** |
| cobertura | minutos observados / esperados | **C1.6** |
| `R²` predictivo fuera de muestra | partición temporal 70/30, bloques de 60 min no solapados | la **predicción falsable** del §2 |
| `H_p` barajado | los mismos retornos, permutados | el **control** del §6.2 |

### ⚠ Roll no es un adorno: cambia veredictos

σ de velas de 1 min **sobreestima** por rebote bid-ask, y sobreestima **más** cuanto
más ancha es la horquilla — la dirección que favorece justo a los pares malos.
Verificado inyectando un rebote de tamaño conocido
(`test_v51_roll_recupera_el_spread_y_desinfla_sigma1`): `s_eff` recuperado 4.932 pb
contra 5.000 inyectados, y `σ₁` de 0.9092 → 0.6614 contra una verdad de 0.6455.

En el ensayo de humo sobre pares sintéticos, un par con `σ₁` crudo de **3.51** cayó a
**2.55** tras Roll y **pasó de superar la compuerta a no superarla**. Sin la
corrección, el cribado habría elegido un par cuya volatilidad era horquilla.

### ⚠ `σ₁` es una EXTRAPOLACION, y hay que leerla como tal

El ajuste vive en [60, 3600] s (rango declarado) y `σ₁` es su valor en **H = 1 s**:
sesenta veces por debajo del punto más corto ajustado. Un error `d` en `H_p` se
convierte en un factor `60^(−d)` sobre `σ₁` — `d = 0.012` ya son ~5 %.

Medido sobre paseos aleatorios con verdad conocida, 8 semillas:

| | valor |
|---|---|
| error de `σ₁` | **−1.0 % a −4.4 %**, siempre del mismo signo |
| error de `H_p` | +0.0005 a +0.0077 |
| error de `σ(3600)` (dentro del rango) | 4.2 %, y más fino que el de `σ₁` |

El test no se limita a tolerar el error: **verifica el mecanismo**, exigiendo que el
error de `σ₁` sea el que predice `60^(−d)` a partir del error de `H_p`. Si no lo
fuera, habría otra causa y la tolerancia sola no la vería.

**Consecuencia operativa:** el reporte publica la columna `+-%` (error típico
relativo de `σ₁`, propagado del ajuste, **cota inferior** porque los residuos de
horizontes anidados están correlacionados) y marca `[MARGINAL en C1.1]` a los pares
que caen a menos de dos errores típicos del umbral. **No cambia el veredicto** — el
umbral es el declarado — pero un par marginal no está decidido por el dato.

### El control barajado es bloqueante

Bajo barajado los retornos son iid y `H_p` **debe** salir 0.5. Si
`|H_p_barajado − 0.5| > 0.05`, el estimador tiene sesgo propio y **el `H_p` de ese par
no se reporta como medida**. Es el control del A.3.4, que ya salvó dos análisis en el
proyecto viejo (la pendiente de la firma, y `C = 0.783` sobre barajados en la v2.2).

### Salida

Tabla completa ordenada por `σ₁` — **no los tres primeros** —, motivos de caída
agregados (dice dónde está el cuello de botella), el veredicto de la predicción
falsable con su `|ρ| mínimo detectable`, y el ganador. Se escribe a
`cribado_etapa1.txt`.

### Presupuesto y reanudabilidad

~87 peticiones de `klines` por par (peso 10) × 30 pares ≈ 26 100 de peso contra
2 400/min: **~11 min mínimo**. Hay token bucket real al 85 % del límite, backoff
exponencial, y **cache en disco por par** (`cribado_cache/`), así que una corrida
interrumpida se reanuda sin volver a descargar. El proyecto ya perdió media hora de
captura por no tener esto (v2.2).

---

## §2 — Etapa 2: captura y réplica sobre UN par

El ganador del §1. **Reutilizar `captura_estacional.py` sin modificar la lógica.**

**Antes de cualquier resultado nuevo**, replicar los cuatro invariantes ya
establecidos en BTC, que sirven de control de tubería:

1. alineado libro-transacción,
2. `E[y_mid·ε] > 0`,
3. `R²` contemporáneo contra Cont-Kukanov-Stoikov,
4. `ν` prediciendo `σ` mejor que `σ`.

**Si esos cuatro no replican, el problema es la captura, no el activo, y no se
sigue.** No se prueba con el siguiente par: se diagnostica la captura.

---

## §3 — Etapa 3: el mismo criterio económico, sin rebajarlo

- Identidad de covarianza multivariante (§C.3.5), con techo `R²_max` **por fila**,
  suelo por rotación circular y la marca falsable/no-falsable, **con `identidad.py` de
  `master`**, que ya la implementa y tiene autotest 6/6.
- `L(H)` medida con sus cinco términos, **con la horquilla real del par nuevo**, que
  allí no será despreciable — a diferencia de BTC.
- Superficie `R(H, θ)` con cap definida como covarianza.
- Regla de parada: la del §5 del preregistro, **acordada antes**.

---

## §4 — Lo que NO se hace

- **No se toca `Micelio.py`.**
- **No se refina el predictor** antes de saber si algún activo cruza.
- **No se rebaja ningún umbral** declarado en el §4 del preregistro.
- **No se cambia de par** a mitad de la etapa 2 porque otro «se ve mejor».

---

## §5 — Estado de ejecución

| pieza | estado |
|---|---|
| identidad multivariante (§C.3.5) | **YA HECHA en `master`**: `identidad.py`, sesión 2026-08-27. No se duplica |
| `cribado_activos.py` | **implementado y probado**; **sin correr contra mercado** (egreso denegado, §7.1) |
| suite de aceptación | **12/12 criterios v5.1** (`python tests_v13.py`) |
| Etapa 2 | no iniciada; depende del §1 |
| Etapa 3 | no iniciada; depende del §2 |

⚠ **El cribado no se ha ejecutado contra mercado.** `fapi.binance.com`,
`api.binance.com` y `data-api.binance.vision` están denegados por la política de
egreso de la sesión remota donde se escribió este código (403 al CONNECT; GitHub sí
resuelve, luego es política y no red). **La corrida real es local, en un solo
comando.** Todo lo que se afirma arriba sobre los estimadores está medido contra
datos sintéticos con verdad conocida, y así está declarado.
