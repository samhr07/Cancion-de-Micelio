# Traspaso de Sesión — Micelio · agosto 2026

Documento de continuidad. Escrito para que otra conversación retome el proyecto sin repetir
trabajo ni reabrir preguntas cerradas.

---

## 0. Estado en una página

**El sistema mide estructura real y no la puede cobrar.**

La v3.2 se ejecutó de punta a punta. El propagador de impacto con memoria (M2) bate al nulo fuera
de muestra con holgura: `IC95 = [+6.86e−3, +8.08e−3]`, `R² = +0.014`, residuo blanco. La
estructura existe y está medida.

Pero el paso económico falla: `q90(|μ̂|) = 11.30 USD/BTC` contra un umbral de `39.05`. **Factor
3.46 en `μ`, factor ~12 en `R²`.** El obstáculo son **4 pb de comisión de ida y vuelta**, no la
microestructura: haría falta una comisión de 0.579 pb por lado contra los 2.000 de VIP 0.

Los pasos 4 y 5 no se ejecutaron. El criterio de abandono del §11 **no se activó**, porque su
condición (3a falla con 3b satisfecha) no se cumple: 3b falla de forma degenerada.

**Dos preguntas quedan abiertas** —permanente contra transitorio, y `ω_G`— y ninguna se arregla
capturando más datos.

**Lo que se aprendió en la sesión y reordena el proyecto:** `H*` no es el horizonte de operación,
es el suelo por debajo del cual operar es imposible. El `R²` requerido cae con el horizonte, así
que la banda viable está en decenas de minutos, no en segundos.

---

## 1. Decisiones arquitectónicas tomadas

Cada una con su evidencia, no solo su veredicto.

| decisión | motivo |
|---|---|
| **MkII no es rival del propagador: es su caso particular** (`β = 0`, `δ = 0.5`) | la elección de arquitectura se convirtió en una celda de una tabla de ajuste |
| **«Delta de ticks» es el nulo M0**, no un candidato | equivale a velocidad constante en reloj de ticks, que el AR(2) ya agotó (`φ₁+φ₂ = 1`) |
| **Todo se mide sobre punto medio**, con precio de transacción como columna de control | `R(1) = −0.0052` es rebote bid-ask; la horquilla mediana es 1 tick exacto |
| **MPC adversario: rechazado** | sin retroacción no hay juego; un oponente modelado es *menos* robusto que el minimax ya implementado |
| **Problema de la secretaria: rechazado** | no se cumple ninguna de sus cinco premisas; la regla de cuándo operar es la banda muerta `\|α\| > c` |
| **NMPC: coste lineal, no cuadrático** | `J = −α·u + c·(u⁺+u⁻) + ½R(u⁺−u⁻)²`. Sin el término lineal no hay banda muerta y el bot nunca se abstiene. Es **v3.4**, toca `Micelio.py` |
| **`γ` no se deriva del techo de riesgo** | hacerlo cancela `α` algebraicamente y el tamaño deja de responder a la señal |
| **`φ′` retirado como variable** | es `1/tamaño medio de operación` por definición; se escribe `q̄` |
| **Piso de embargo de 1 950: retirado** | `H*_ticks` medido es 2 233 y ata él |

---

## 2. Resultados de la v3.2 (commits `35d6c7e`, `0db4fcf`)

### 2.1 La decisión

```
PASO 1  compuerta (1 031 155 ticks continuos / 848 540)      PASA
PASO 2  M2 contra M0, IC95 [+6.86e-3, +8.08e-3]              PASA
PASO 3  q90(|mu|) = 11.30  vs  1.5*c(u) = 39.05              FALLA
PASOS 4 y 5                                                  no ejecutados
```

`μ̂` de M1 sale **idénticamente cero**, y no es un fallo: con impacto permanente el movimiento ya
ocurrió y no queda nada que capturar. Es el contenido económico de MkII.

### 2.2 Lo que salió y no estaba previsto

**El signo de `D` se invierte entre observables.** Misma M2, misma muestra:

| observable | `β` | `D` |
|---|---|---|
| punto medio | −0.160 | 11.53 (~~núcleo creciente~~ → **ajuste no identificado**, v4.1 §3.4) |
| precio de transacción | +4.69 | 0.107 |

⚠ **Reclasificado el 2026-08-12.** `β < 0` con `τ₀` en cota es firma de mala especificación del
estimador, no del mercado. `D = 11.53` no se cita como medición de la forma del núcleo.

El «impacto transitorio» de la v3.1 §2 era el rebote bid-ask. **El preregistro lo cazó porque
obligaba a llevar las dos columnas siempre.**

**`D > 1` no cabe en la dicotomía escrita.** La regla del §5.2 era de una cola; aplicada
literalmente habría dicho «impacto permanente, MkII bien especificada» sobre datos que dicen lo
contrario. Se leyó la misma simulación entera: `D = 11.53 > q95 = 7.27` → se rechaza `D = 1` por
arriba.

**El contraste se derrumba a la `γ` real.** Con `γ̂ = 0.798`, bajo `D = 1` verdadero el estimador
devuelve cualquier cosa entre 0.0 y 7.3. La malla de calibración se midió a `γ = 0.3` y no
transporta.

**Super-difusión.** `H_p = 0.591` contra 0.488 del control barajado. No replica la sesión (e),
que dio difusivo sobre `captura_larga`.

**`δ = 0` gana en validación** en los dos observables y los dos modelos; `δ = 0.5` (MkII) da 1.77×
menos. El impacto informa el **signo**, no el volumen.

**Control de fuga barajado:** cae al 0.1 % de la ganancia, factor 1 000.

### 2.3 `ω_G` — el catch más importante de la sesión

Primera ejecución: `2ΔLL = −149.078`, `p = 0.6250`. **Imposible en modelos anidados**: M2-osc
contiene a M2. Era un fallo del optimizador (Nelder-Mead, 5 parámetros, arranque fijo), y el 37.5 %
del nulo caía aún más abajo.

Leído sin mirar el signo, ese `p` decía «propagador monótono», y la consecuencia declarada del
§5.3 es **borrar `ω_m,max` y `γ_ω` de `constantes_micelio.py`**. Se habrían borrado sobre un fallo
numérico.

Corregido sembrando el modelo grande en la solución del anidado: `2ΔLL = +134.009`, `p = 0.0000`,
0 negativos en 40 sorteos. **Pero el rechazo no sostiene «oscilador forzado»:**

```
omega_G = +0.000009 rad/tick  ->  periodo 737 227 ticks = 11.7 h
K = 4466                      ->  razon 165.1
el coseno recorre el 0.606 % de un ciclo en [0, K]
```

No oscila: actúa como una **inclinación lenta** del núcleo. `ω_G` queda **sin decidir**, y
`ω_m,max` y `γ_ω` **ni se borran ni se resucitan**.

**El control 19** (`migracion_v32.py` pasa a 20/20) falla solo con núcleo creciente — justo el
régimen donde cayó el dato real. Los 18 controles anteriores usaban núcleos decrecientes.

> **Principio que se sigue de ahí:** los controles positivos deben cubrir **el complemento** de lo
> que crees, no la vecindad de tu hipótesis.

### 2.4 Dos lecturas post-hoc, declaradas

1. `H*_ticks` por punto fijo (la firma de volatilidad no tiene plateau; la fórmula da 6 347 o
   1 688 según qué `σ` se elija).
2. La lectura de dos colas de `D`.

**Ambas son correctas y ambas se tomaron después de ver el dato.** Van como lectura del
preregistro, no como enmienda. El reporte final debe listarlas en sección propia con la fecha y
qué se había visto ya, para que un lector externo las descuente.

---

## 3. Línea de ejecución (v4.0) — cuánto cuesta operar

Todo sobre `captura_v32`, dato en disco.

| cantidad | valor |
|---|---|
| tasa de llenado maker (60 s) | **41.3 %** |
| markout 1 s | media −3.13, mediana −1.15 |
| markout 50 s | media −5.78 |
| aporte del 5 % peor a la media | 30 % → **toxicidad general, no de cola** |
| `ρ`(flujo previo, markout) | +0.407 (1 s), +0.320 (5 s), +0.144 (50 s) |
| rotación de cola | 10.5× (**cota inferior**, `bookTicker` estrangulado) |
| llenados dependientes de la atribución | 1.2 % (era 30.7 %) |
| `c(u)` con salida forzada | 37.65 USD/BTC |
| **`H*`** | **50 s → ~103 s** |

`ρ` decae monótonamente con el horizonte: firma de efecto de microestructura real, no de
correlación espuria. **Los peores llenados son anticipables por el flujo previo.**

**`C_respaldo` sigue sin medir.** Es el término `(1−p)·C_respaldo` de la ecuación de coste maker,
no existe en ningún documento del proyecto, y sin él maker gana siempre porque no paga nada al
fallar. Punto de equilibrio estimado: si el precio se escapa más de ~9 USD/BTC mientras esperas,
la orden límite fue un error.

---

## 4. Hipótesis

### 4.1 Corroboradas

| hipótesis | evidencia |
|---|---|
| Propagador de impacto con forzamiento medido existe | `R(τ) > 0` en 799/799 rezagos, factor 150–700 sobre barajado; M2 bate a M0 con IC que excluye 0 |
| El impacto informa el signo, no el volumen | `δ = 0` gana en validación en los dos observables |
| Selección adversa medible y anticipable | markout −5.78 USD/BTC; `ρ` con flujo previo decayendo con el horizonte |
| Tamaño medio de operación asociado a volatilidad | test disjunto 14/16, `p = 0.0042`. **Es Jones–Kaul–Lipson (1994)**, no un hallazgo nuevo |
| `q̄` (ex-`φ′`) es proxy de `G₀` | `ρ = +0.6323`, 8/8; nulo por desplazamiento circular −0.0908, 2/8 |

### 4.2 Refutadas — no revisitar sin evidencia nueva

- **Frecuencia estructural modal `ω_m` vía EMD/HHT**: artefacto de ventana de análisis.
- **Fuerza recuperadora AR(2)**: `k ≈ 4.67e−7`, `p ≈ 0.62`, `φ₁+φ₂ = 1`.
- **Subdifusión a 30–300 s**: artefacto de mezclar dos relojes en el estimador de volatilidad.
- **`ω_G` como oscilador forzado**: fuera de banda, 0.606 % de un ciclo en el rango ajustado.
- **`φ′` como magnitud propia**: es `1/q̄` por definición (`corr = +1.000000` exacto).
- **`φ′` como modulador de `Q_k`**: devolvería al residuo una cantidad medible.

### 4.3 Abiertas

| pregunta | estado |
|---|---|
| **Permanente contra transitorio** | 4 estadísticos muertos. **Puede estar ya contestada** — ver §5 |
| **`ω_G`** | sin decidir; fuera de banda. `ω_m,max` y `γ_ω` en el limbo |
| **¿Es `β` no estacionario?** | −0.23, +0.04, +0.21, −0.16 en cuatro tramos. Media +0.0058 |
| **¿En qué reloj vive el propagador?** | §7 de la v3.2 no ejecutado; exige captura completa |
| **Difusividad más allá de 12 min** | alcance máximo verificado 711 s. Toda la fila de 1 h en adelante es extrapolación |
| **Curva en U: ¿volatilidad o actividad?** | `σ_pared = √ν·σ_tick`; el factor mecánico esperable es ~1.92× y la variación observada 2–3× |
| **¿Aporta el lado de límite sobre el flujo de transacciones?** | M2+L no ejecutado |
| **`C_respaldo`** | sin medir |

---

## 5. La ruta que puede cerrar la pregunta abierta

Salida de la revisión de literatura de esta sesión. **`H = (2−γ)/2 − β`**, con las cantidades ya
medidas:

| cantidad | valor |
|---|---|
| `γ̂` | 0.798 |
| `H_p` | 0.591 |
| `(2−γ̂)/2` (valor si el impacto fuera permanente) | 0.601 |
| **`β` implícita** | **+0.010** |
| `β` de difusividad, `(1−γ̂)/2` | +0.101 |

**Impacto esencialmente permanente**, y la difusividad queda lejos. Es una **tercera ruta al mismo
parámetro** que no hereda la degeneración `β`–`τ₀` que mató a los estadísticos 2, 3 y 4. Coherente
con que `μ̂` de M1 saliera idénticamente cero.

⚠ **Sin verificar contra fuente primaria.** La relación es asintótica y vale en un régimen de
escala concreto.

### 5.1 Y una corrección que afecta retroactivamente

Con memoria larga, `Var(media) ∝ N^(−γ)`, no `N^(−1)`:

```
N_eff ≈ N^0.798 ≈ 63 000   de 1 031 155 ticks   ->   factor ~16
```

Eso explica los cuatro estimadores fallidos mejor que ningún defecto de diseño. Y la longitud de
bloque debe crecer con `N` (entre `N^(1/3)` y `N^(1/2)`), no quedarse en `5·embargo`. **Puede
ensanchar retroactivamente los IC de la v3.2, incluido el del paso 2.**

---

## 6. La aritmética que ordena las expectativas

```
1 tick                 = 0.10 USD/BTC = 0.010 pb
comision ida y vuelta  = 4 pb         ~ 384 ticks
horquilla cruzada      = 1 tick       = 0.26 % del coste total
```

**Hay que predecir ~390 ticks de movimiento para pagar la ida y vuelta.** Un agotamiento de cola
da 1 tick.

Y el `R²` requerido cae con el horizonte, así que `H*` es un **suelo**, no un horizonte de
operación. Con `H_p = 0.591` (super-difusión), el `R²` requerido a 10 min es ~2.5 % y el medido a
2 min es 1.4 % — **la primera vez que un número medido y uno requerido están en el mismo orden**.
Con la reserva de que `H_p` no ha replicado.

⚠ **Verificar `c(u)`:** 26.03 USD/BTC con 4 pb implica BTC ≈ 65 000. Si el tramo corrió cerca de
96 000, el factor del paso 3 no es 3.46 sino ~5.1. **La tabla de costes debe estar en puntos
básicos**, no en USD/BTC fija desde la v3.1.

---

## 7. Errores encontrados y corregidos en esta sesión

El registro importa tanto como los resultados: nueve enmiendas del preregistro, **todas motivadas
por un control positivo y ninguna por un resultado**.

| # | error | quién lo encontró | consecuencia |
|---|---|---|---|
| 1 | `β` y `τ₀` no identificados por separado (`β̂ = 5.0` sobre verdad `β = 0`) | control positivo | invalidó la regla «¿el IC de `β̂` contiene 0?» |
| 2 | Compuerta `375·embargo` → **`380·embargo`** | Code, corrigiendo a Claude | 14 bloques presentados como 15 |
| 3 | `D` sin nulo: **33 % de falsos positivos** contra 5 % nominal | control positivo | umbral pasa al `q05` simulado |
| 4 | `G(∞) = 0` sin suelo `f_∞` | heredado de la v3.1 §3 | núcleo reformulado con `f_∞` |
| 5 | §7 muerto: `τ₀` no interpretable → regresión sin sentido | Claude | sustituido por dispersión de `D` en dos relojes |
| 6 | «NIVEL BARRIDO» clasificado como no llenado | Claude preguntando por la definición | llenado 16.2 % → **41.3 %**; markout 1 s −1.39 → −3.13 |
| 7 | «99.1 % cancelación», mal por factor 38 | Code, contrastando contra cantidad independiente | con L1 cancelación y reposición no son separables |
| 8 | `2ΔLL < 0` en modelos anidados | Code, mirando el signo | habría borrado `ω_m,max` y `γ_ω` sobre un fallo numérico |
| 9 | Argumento dimensional de `φ′` (ticks = transacciones, no ticks de precio) | Code, verificando antes de correr | todo el argumento de λ de Kyle era falso |
| 10 | Derivación de `γ` en la que `α` se cancela | Claude | el tamaño dejaba de responder a la señal |
| 11 | `minQty` (0.001 BTC ≈ 96 USD) manda sobre notional de 50 USDT | Claude | la tabla de granularidad estaba desplazada ×2 |

⚠ **Patrón recurrente, cinco apariciones:** umbral asintótico usado donde la distribución finita es
otra. Y **mezclar dos relojes en un estimador**, otras tantas.

---

## 8. Documentos producidos en esta sesión

| documento | qué es |
|---|---|
| `ORDEN_TRABAJO_MIGRACION_3_2.md` | comparación anidada M0/M1/M1′/M2; ejecutada |
| `PREREGISTRO_3_2.md` | 9 enmiendas acreditadas; `190edda` → `2ad0700` |
| `ORDEN_TRABAJO_EJECUCION_4_0.md` | captura estacional, dos relojes, `C_respaldo`, resistencia |
| `PLAN_CAPITAL_5_0.md` | **condicional**; granularidad, Merton con techo, validación, abandono |
| `ORDEN_TRABAJO_TICK_GRANDE_3_3.md` | reformulación en tick grande; ruta `H = (2−γ)/2 − β` |

**Capturas vivas:** `captura_v33` (dato de la v3.2, bloque 201+) y `captura_estacional`
(lanzada 16:30:53 UTC, 21 días).

**Suites:** `tests_v13` 56/56 · `ssa` 11/11 · `migracion_v32` **20/20** · `cola` 9/9.
`Micelio.py` **sin cambios desde la v2.2**.

---

## 9. Trabajo futuro, por prioridad

**Inmediato, sobre dato en disco, sin bloqueos:**

1. **`H = (2−γ)/2 − β`** — puede cerrar la pregunta abierta. Verificar fuente primaria primero.
2. **`N_eff` y longitud de bloque creciente** — barato, y puede invalidar IC ya publicados. Debe
   correr antes que nada que use IC.
3. **`C_respaldo` sobre `captura_v32`** — el que más cambia el diseño. Pendiente desde que se
   identificó.
4. **`c(u)` en puntos básicos** al precio real del tramo.
5. **`η̂ = N^(c)/(2N^(a)))`** — escalar de dos contadores; compuerta del marco de tick grande.
6. **Difusividad por klines de 1 minuto** — con `H_p` inestable, pasa de útil a decisivo.
7. **Cuantiles empíricos del NIS** (v2.1 §6) — en cola cuatro sesiones, curtosis 1179.7.
8. **Sellar la pendiente `log G₀` vs `log q̄`** como predicción de `δ` (pendiente = `1 − δ`).

**Cuando la captura estacional acumule:** curva en U en dos relojes con regla de colapso; azar con
riesgos competitivos; §7 de la v3.2 (reloj del propagador).

**Requiere autorización explícita** (preregistro nuevo): reajustar M2 a `K` de 10 min. ⚠ Necesita
**~62 h continuas**; un hueco de 2 h parte el tramo. La de 1 h (~15.6 días sin cortes) está fuera
de alcance con el montaje actual. **Desactivar la suspensión automática de la máquina.**

**Condicionado a que el paso 3 pase en alguna banda:** v3.4 (coste lineal en el NMPC),
`PLAN_CAPITAL_5_0.md`, demo, despliegue.

---

## 10. Reglas heredadas que no se renegocian

- **`Micelio.py` no se toca durante fases de prueba.** Código de análisis estrictamente aparte.
- **Preregistro antes de mirar.** Enmiendas con hash antes/después, motivo, y constancia de si
  había resultado a la vista.
- **Controles positivos antes de dato real**, y cubriendo el complemento de la hipótesis.
- **Todo contraste de verosimilitud contra nulo simulado**, nunca contra tabla χ².
- **NIS y `ρ₁` siempre juntos.** Inflar `Q` siempre «mejora» el NIS.
- **Punto medio y precio de transacción, las dos columnas siempre.**
- **No abrir un quinto estadístico** para permanente/transitorio: la pregunta necesita
  reformularse, y eso es decisión de Samuel.
- **Comprensión antes de despliegue.** Si no se puede explicar por qué el bot abrió una posición
  concreta, no se distingue un fallo de una mala racha.

---

## 11. Para quien retome esto

**La colaboración:** Samuel decide y es el único portador de continuidad entre sesiones. Claude
revisa arquitectura, metodología y matemáticas, y produce órdenes de trabajo. Code implementa y
ejecuta. El traspaso a Code es `CLAUDE.md` en la raíz del repositorio.

**Lo que más ha funcionado:** que Claude y Code se corrijan mutuamente con Samuel decidiendo cuál
versión sigue. En esta sesión, cada error de la tabla del §7 lo cazó el otro lado. Ninguna de las
nueve enmiendas habría salido de dos modelos de acuerdo entre sí.

**Lo que más ha costado:** presuponer comportamiento donde no lo hay. `ω_m` costó ocho versiones.
El patrón de rescate es siempre el mismo — una premisa necesita un nulo **antes** de construir
encima, y la literatura ayuda porque trae los nulos ya pensados, no porque exima del test.

**Y la reserva honesta sobre el estado:** el sistema mide estructura real. Todavía no hay
evidencia de que la pueda cobrar, y «no hay señal explotable» sigue siendo un desenlace vivo y
perfectamente reportable.
