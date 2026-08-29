# Orden de Trabajo — Reformulación en Tick Grande y Cierre de la v3.2 · v3.3

Sucede a la v3.2, que quedó **cerrada en el paso 3** (señal 3.46× por debajo del coste) y con dos
preguntas abiertas: permanente contra transitorio, y `ω_G`.

**Lo que este documento NO hace: abrir un quinto estadístico.** Cuatro han muerto y el patrón
apunta a que la pregunta está mal planteada sobre el observable equivocado, no a que falte un
estimador mejor.

**Lo que hace:** tres rutas nuevas, todas sobre dato en disco, ninguna con estadístico nuevo.

⚠ **Advertencia sobre las fuentes.** `text.txt` … `text_5.txt` son resúmenes generados por otro
modelo. Son **pistas, no citas**. Toda fórmula que se implemente debe verificarse contra el
artículo primario antes de escribir código. Este proyecto ya perdió ocho versiones por construir
sobre una premisa no verificada.

---

## §1. La ruta que puede cerrar la pregunta sin estadístico nuevo

Bajo el marco del propagador, el exponente de Hurst del precio, la memoria del flujo y el
decaimiento del núcleo están ligados:

```
H  =  (2 − γ)/2  −  β
```

Comprobación de coherencia: con `β = (1−γ)/2` (condición de difusividad de Bouchaud) sale
`H = 1/2` exacto; con `β = 0` (impacto permanente) sale `H = (2−γ)/2`.

**Con las cantidades ya medidas en la v3.2:**

| cantidad | valor | origen |
|---|---|---|
| `γ̂` | 0.798 | autocorrelación de signos |
| `H_p` | 0.591 | firma de volatilidad, control barajado 0.488 |
| `(2−γ̂)/2` | 0.601 | valor de `H` si el impacto fuera permanente |
| **`β` implícita** | **+0.010** | `0.601 − 0.591` |
| `β` de difusividad | +0.101 | `(1−γ̂)/2` |

**`β` implícita ≈ 0.010, es decir, impacto esencialmente permanente.** Y la difusividad
(`β = 0.101`) queda a varios errores estándar.

Esto es una **tercera ruta al mismo parámetro**, y usa dos cantidades mucho mejor estimadas que
`β̂` directa. `γ̂` sale de un promedio sobre millones de pares; `H_p` de una regresión sobre la
firma de volatilidad. Ninguna de las dos hereda el problema de identificabilidad `β`–`τ₀` que
mató a los estadísticos 2, 3 y 4.

### 1.1 Lo que hay que hacer

1. **Verificar la relación `H = (2−γ)/2 − β` contra fuente primaria.** Es asintótica y vale en un
   régimen de escala concreto; si no aplica al rango de rezagos medido, todo lo anterior cae.
2. **Errores estándar de `γ̂` y `H_p`** por bootstrap por bloques (con la corrección del §5), y
   propagarlos a `β` implícita.
3. **Contraste declarado:** ¿el IC de `β` implícita excluye `+0.101` (difusividad)? ¿Contiene 0
   (permanencia)?
4. **Reportar los tres valores de `β` juntos:** implícita (+0.010), ajuste sobre punto medio
   (−0.160), ajuste sobre precio de transacción (+4.69). La discrepancia **es** el resultado, no
   un problema a resolver eligiendo uno.

⚠ Si la ruta funciona, el desenlace probable es **«no se puede rechazar impacto permanente, y se
rechaza difusividad»**. Eso selecciona M1/MkII por la puerta de atrás — y es coherente con que
`μ̂` de M1 saliera idénticamente cero: con impacto permanente no queda nada que capturar.

---

## §2. La compuerta de tick grande: `η̂`

`η` es el parámetro de zona de incertidumbre de Robert–Rosenbaum. Con `N^(c)` continuaciones
(el precio salta en la misma dirección que el salto previo) y `N^(a)` alternancias:

```
η̂ = N^(c) / (2·N^(a))          [VERIFICAR contra Robert & Rosenbaum 2011]
```

Es un **escalar de dos contadores**. Se computa en minutos sobre dato en disco.

| resultado | consecuencia |
|---|---|
| `η̂ < 1/2` | régimen de tick grande confirmado; §3 y §4 aplican |
| `η̂ ≥ 1/2` | el marco no aplica; §3 y §4 **no se ejecutan** y se reporta así |

⚠ **Y una duda que hay que resolver antes de aceptar el marco.** BTCUSDT es de tick grande por el
criterio operativo —horquilla mediana clavada en 1 tick— pero su tick **relativo** es diminuto:

```
α/P = 0.1/96 000 = 1.04e−6  ≈  0.010 pb por tick
```

En la literatura de renta variable, «tick grande» significa `α/P` de 1e−4 a 1e−3: **dos o tres
órdenes de magnitud mayor**. Que las dos condiciones coexistan aquí no está en los artículos
citados. Reportarlo explícitamente como límite de la analogía, no darlo por bueno.

---

## §3. La aritmética del tick, y qué puede y no puede rescatar esto

```
1 tick                    = 0.10 USD/BTC  = 0.010 pb
comisión ida y vuelta     = 4 pb          ≈ 384 ticks
horquilla cruzada         = 1 tick        = 0.26 % del coste total
|μ̂| medida (q90)          = 11.30 USD/BTC = 113 ticks
umbral 1.5·c(u)           = 39.05         = 390 ticks
```

**Hay que predecir ~390 ticks de movimiento para pagar la ida y vuelta. Un agotamiento de cola da
1 tick.** La capa de microestructura opera a una escala ~400× por debajo del umbral de
rentabilidad.

⚠ **Consecuencia que ordena el resto del documento: la reformulación en tick grande arregla la
ciencia, no la economía.** Explica `β < 0`, `H_p = 0.591` y la inversión de signo entre
observables. No mueve el paso 3 ni un punto básico.

⚠ **Y verificar `c(u)`.** 26.03 USD/BTC con 4 pb implica BTC ≈ 65 000. Si el tramo corrió cerca de
96 000, `c(u)` son 38.4 y el factor pasa de 3.46 a **5.1** (de 12 a **26** en `R²`). La tabla de
costes debe expresarse en **puntos básicos** y recalcularse al precio del tramo, no en USD/BTC
fija desde la v3.1.

---

## §4. El propagador sobre el precio eficiente

El modelo de zonas de incertidumbre separa el precio eficiente `X_t` (difusivo, continuo) del
precio observado `P_t` (a saltos, anclado a la rejilla). `P_t` salta solo cuando `X_t` escapa de
la banda `±(α/2 + η̂α)`.

**Es sobre `X_t` donde la pregunta permanente/transitorio está bien planteada**, porque `X_t` no
tiene rebote bid-ask — y el rebote es lo que la v3.2 identificó como causa de la inversión de
signo entre columnas.

Reconstruir `X_t` con `η̂` del §2 y reajustar M2 sobre él. Reportar `β`, `D` y `H` sobre las
**tres** series: `X_t`, punto medio, precio de transacción.

⚠ La reconstrucción **introduce el modelo en el dato**. Si M2 sobre `X_t` sale limpio, parte de
esa limpieza es que `X_t` se construyó suponiendo lo que M2 mide. Declararlo, y no presentar el
ajuste sobre `X_t` como confirmación independiente del §1.

---

## §5. La predicción falsable — el test que de verdad discrimina

La literatura afirma algo contrastable y arriesgado:

> El propagador crece mientras la cola resiste, y **el decaimiento empieza en el rezago del
> agotamiento de la cola**, no en el rezago 0.

Es falsable con lo que ya tienes: `cola.py` mide agotamientos, y el propagador está implementado.

**Procedimiento:** estimar el rezago mediano hasta agotamiento `τ_dep` en la muestra, y medir
`G(τ)` **condicionado** a que la cola siga viva contra condicionado a que ya se haya agotado.

| resultado | lectura |
|---|---|
| `G` crece antes de `τ_dep` y decae después | la literatura aplica; el marco es correcto |
| `G` crece en ambos regímenes | el crecimiento no es por agotamiento; el marco no explica el dato |
| sin diferencia | `τ_dep` no es la escala relevante |

**Éste es el apartado con más valor de información del documento**, porque puede refutar el marco
recién adoptado. Un marco que explica tres anomalías a posteriori y no arriesga ninguna predicción
es exactamente lo que fue `ω_m`.

---

## §6. Corrección al bootstrap por bloques

Con memoria larga, `Var(media) ∝ N^(−γ)`, no `N^(−1)`. Con `γ̂ = 0.798`:

```
N_eff ≈ N^γ = (1 031 155)^0.798 ≈ 63 000        [VERIFICAR la forma exacta]
```

**Factor ~16 de pérdida de información.** Eso explica los cuatro estimadores fallidos mejor que
cualquier defecto de diseño, y **no se arregla capturando más**: con `β` cambiando de signo entre
tramos (−0.23, +0.04, +0.21, −0.16), promediar más días promedia algo que se mueve.

Y la longitud de bloque **debe crecer con `N`**, entre `O(N^(1/3))` y `O(N^(1/2))`, no quedarse en
la regla fija de `5·embargo`. Recalcular los IC del §1 con la longitud correcta.

⚠ **Esto puede ensanchar retroactivamente los IC de la v3.2**, incluido el del paso 2
(`[+6.86e−3, +8.08e−3]`). Rehacerlo y reportar si el paso 2 sigue pasando. Si no pasa, es un
resultado y se reporta como tal.

---

## §7. Lo que NO se hace

- **Ningún quinto estadístico** para permanente/transitorio. §1 y §5 son rutas, no estimadores.
- **No se reescribe el EKF** con difusión fraccionaria ni procesos de Hawkes. `text.txt` §8 lo
  propone; es una reescritura de arquitectura sobre un resumen no verificado, con el paso 3 ya
  fallado. `Micelio.py` no se toca.
- **No se implementa el QRM completo** (intensidades `λ^L`, `λ^C`, `λ^M` condicionales al estado).
  Exige `@depth`, que sigue sin ingerirse.
- **No se resucita `ω_m`.** `ω_G` quedó fuera de banda y sin decidir; `ω_m,max` y `γ_ω` **no se
  borran ni se restauran**.
- **No se reabre el conjunto de prueba de la v3.2.** §1, §2, §4 y §5 corren sobre entrenamiento y
  validación, o sobre `captura_v32`. Cualquier cosa que exija prueba requiere preregistro nuevo y
  autorización explícita.

---

## §8. Modos de fallo

| # | fallo | firma | detección |
|---|---|---|---|
| 1 | Fórmula tomada del resumen sin verificar | resultado limpio que no replica | cita primaria obligatoria en el reporte, por fórmula |
| 2 | `H = (2−γ)/2 − β` fuera de su régimen | `β` implícita incoherente con ambos ajustes | §1.1, verificar rango de validez |
| 3 | `η̂` calculada sobre saltos de más de 1 tick | `η̂ > 1/2` espurio | definir el evento de salto antes de contar |
| 4 | `X_t` reconstruido y usado como confirmación independiente | circularidad | §4, declarado |
| 5 | El marco explica todo y no predice nada | tres anomalías resueltas, cero refutable | §5 es obligatorio, no opcional |
| 6 | `c(u)` en USD/BTC a precio equivocado | el paso 3 falla por un factor distinto del real | §3, expresar en pb |
| 7 | Bloques fijos bajo memoria larga | IC demasiado estrechos en toda la v3.2 | §6, longitud creciente |
| 8 | `τ_dep` estimado sobre la misma muestra que `G` | el condicionamiento se autoconfirma | partición o validación cruzada |
| 9 | Reapertura silenciosa del conjunto de prueba | veredicto sobre dato ya visto | §7, autorización explícita |

---

## §9. Criterios de aceptación

**§1 — Relación de coherencia**
- [ ] Relación verificada contra fuente primaria, con su rango de validez citado.
- [ ] `γ̂` y `H_p` con IC por bootstrap corregido (§6).
- [ ] `β` implícita con IC propagado; contraste contra `+0.101` y contra `0`.
- [ ] Los tres valores de `β` reportados juntos, sin elegir uno.

**§2 — `η̂`**
- [ ] Fórmula verificada contra Robert & Rosenbaum.
- [ ] Definición del evento de salto declarada antes de contar.
- [ ] `α/P` reportado junto a `η̂`, con la reserva sobre el tick relativo.

**§3 — Aritmética**
- [ ] Tabla de costes en puntos básicos, al precio real del tramo.
- [ ] Paso 3 recalculado; factor actualizado.

**§4/§5 — Propagador y predicción falsable**
- [ ] `β`, `D`, `H` sobre las tres series.
- [ ] `G(τ)` condicionada a cola viva contra agotada, con `τ_dep` estimado fuera de la muestra de
      ajuste.
- [ ] Circularidad de `X_t` declarada en el reporte.

**§6 — Bootstrap**
- [ ] `N_eff` reportado junto a `N` en toda tabla con IC.
- [ ] IC del paso 2 de la v3.2 recalculados; se indica si sigue pasando.

**Transversal**
- [ ] `Micelio.py` sin cambios.
- [ ] Conjunto de prueba de la v3.2 no reabierto.
- [ ] Sección propia listando las lecturas post-hoc de la v3.2 (`H*_ticks` por punto fijo, `D` de
      dos colas) con la fecha y qué se había visto ya.
- [ ] Suites en verde.

---

## §10. Prioridad

1. **§1** — puede cerrar la pregunta abierta con cantidades ya medidas. Horas de trabajo.
2. **§6** — barato y puede invalidar IC ya publicados. Debe correr antes que nada que use IC.
3. **§2** — escalar de dos contadores, compuerta para §4 y §5.
4. **§3** — recálculo aritmético, cierra la duda del precio.
5. **§5** — el que puede refutar el marco. Después de §2.
6. **§4** — el más caro y el más circular. El último.

Nada de esto es urgente. Las dos capturas siguen vivas y no se rompe nada si esperan.
