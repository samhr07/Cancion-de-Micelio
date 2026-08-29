# Orden de Trabajo — Relojes y Reconstrucción del Acoplamiento con el Mercado v2.0

Sucede a `ORDEN_TRABAJO_RIESGO_1_3.md`, que queda como registro. La v1.3 está **ejecutada y
verificada** (33/33, ver `CLAUDE.md` sesión 2026-08-04); nada de lo que construyó se deshace
aquí.

Este documento es de **reorganización estructural**, no de calibración. Cambia dónde vive el
tiempo en el sistema.

---

## §0. Registro de errores de la v1.3 (del documento, no del código)

Ambos detectados por Code en implementación y ya corregidos. Se anotan para que el rastro
quede completo:

| Error de la v1.3 | Realidad |
|---|---|
| `minNotional = 100 USDT`, y de ahí "ninguna orden es legal" | Es **50 USDT**. La conclusión era falsa. El margen sigue estrecho: por debajo de BTC = 50 000 vuelve a serlo, de ahí que `cantidad_minima_legal` dependa del precio. |
| Clamp de posición `u ← min(u, nocional_max/S − |inv|)` sin signo | Al saturar anula **ambas** direcciones, incluida la que reduce posición: **rompe la ruta de cierre del halt.** Corregido direccionalmente. |

El segundo es el más grave que ha producido este flujo de trabajo hasta ahora: una capa de
seguridad que inhabilitaba el mecanismo de seguridad al que el propio documento asignaba el
criterio de aceptación más importante. Sirva de precedente para el §8.

---

## §1. Diagnóstico verificado

Cuatro hallazgos, todos comprobados contra el código. Los tres primeros son **estructurales**:
no dependen de la implementación, sino de la forma del bloque de entorno y de la cadencia de
los hilos.

### 1.1 El EMD tamiza una escalera

`Micelio.py:1461` muestrea con reloj de pared:

```python
if (t_now - t_ultima_muestra) >= 0.9 * PERIODO_MUESTREO_EMD:
```

Toma `env_arr[0]["P_spot"]` cada 0.5 s **independientemente de si cambió**. Con el feed
degradado a 0.76 Hz, cerca de la mitad de las muestras del buffer son duplicados literales.

El comentario inmediatamente encima se preocupa de que el espaciado **irregular** rompa los
splines de la EMD. Resolvió la uniformidad de la malla *temporal* justo mientras la malla de
*valores* se volvía escalera: uniforme en `t`, constante a trozos en `S`. Y la transformada de
Hilbert de un escalón tiene contenido en todo el espectro — es la singularidad que el propio
docstring de `hht.py` advierte.

**Consecuencia:** `ω_m` y `C` pueden estar contaminados en origen. Esto es una explicación del
resultado del A/B **independiente** de la de propagación entre correcciones, y hay que
descartarla antes de volver a decidir sobre `A`.

### 1.2 `n_ticks` cuenta paquetes, no transacciones

`Micelio.py:908`: `n_ticks += 1` una vez por iteración del bucle de red. Y `ENV_DTYPE:116`
tiene **una sola casilla escalar** de `P_spot`: el bloque de entorno estructuralmente no puede
transportar un lote.

Pero `aggTrades` **no es un muestreo, es entrega completa por lotes**: un sondeo devuelve todos
los aggTrades del intervalo, paginables por `fromId`. En BTCUSDT perp el ritmo es de orden
10–50 aggTrades/s. Las 0.76 mediciones/s del reporte son la **tasa de paquetes**, no la de
información.

Propagación del error de ν:

| Magnitud | Efecto |
|---|---|
| `ω_m = f/ν` | **sobreestimada** por el factor de lote |
| `ρ_k = 1 + γ_ω·ω_m + γ_Q·ΣQ` | **Q inflada** por esa vía |
| `OMEGA_M_MAX = 3e-3` | el techo puede estar saturando sin que se note |
| Piso de Nyquist `ω_m ≤ 0.5/Ticks` | se aplica sobre el tick equivocado |
| `Φ = ΣQ·ω_m` | hereda el sesgo, y con él Ψ y Ω |
| `c²_vol = k·ω_m·ν` | **inmune** — como `ω_m = f/ν`, ν se cancela exacto: `c²_vol = k·f` |

La última fila se deja escrita para que nadie la "corrija" dos veces.

### 1.3 `Q` no depende de Δt

`Micelio.py:1148`:

```python
Q_k = rho_k * np.diag([q_S, q_base * 1e-2, q_base])
```

No hay Δt. `Σ AⁱQAⁱᵀ` con `Q` constante crece con el **número de pasos**, no con el tiempo
transcurrido. Es decir: **hoy la tasa del bucle cambia el filtro.** El ruido de proceso
inyectado es proporcional a cuántas veces despertó el planificador, y `q_base` está calibrada
en silencio contra ~90 Hz. Cambiar `PERIODO_HILO_RAPIDO` recalibra el filtro sin que nada lo
diga.

Esto **corrige una afirmación previa**: se dijo que N pasos pequeños equivalen a uno grande.
Vale para `A`, que es `exp(FΔt)` exacta. **No vale para `Q` tal como está escrita.**

### 1.4 El bucle es ~118× más rápido que la medición

No es una violación de Nyquist: es sobremuestreo, no submuestreo, y no produce aliasing en el
lazo de control. Con períodos modales de 33–188 s y mediciones a 0.76 Hz hay entre 25 y 143
muestras por ciclo — holgadamente por encima de Nyquist para observar `ω_m`.

**Sí hay aliasing genuino en un sitio**: la microestructura por debajo de ~1.3 s se pliega
dentro de la banda observada, y contamina precisamente `r_S,base`. Se aborda en §4.3.

---

## §2. Principio rector: dos relojes, un puente

La estimación y la ejecución viven en tiempos físicamente distintos. Hasta ahora el sistema
usaba uno solo, y todas las conversiones ocurrían implícitamente en el punto de uso.

### 2.1 Reloj de transacciones — para la estimación

Los retornos son más estacionarios y más gaussianos en tiempo de transacción que en tiempo de
calendario. El proyecto ya tiene montada esa contracción: ν, `ω_m` en 1/Ticks, `ρ_k`, `c²_vol`.
El EAKF y la cadena EMD→Hilbert pasan a este reloj.

### 2.2 Reloj de pared — para la ejecución

Irreductible. `τ_d` y `τ_max` (Sec. 6.5), las 7 guardas de riesgo, los rate limits, el funding,
la expiración de órdenes, el horizonte de Loeper y el NMPC. Nada de esto tiene significado en
ticks.

### 2.3 ν es el único puente

`ν [Ticks/s]` ya existe y ya es la conversión. Lo que cambia es que deja de ser una constante
más y pasa a ser **el único punto de contacto entre los dos relojes**, con las obligaciones que
eso impone (§2.4).

### 2.4 Reglas no negociables de nomenclatura y conversión

Las dos trampas de la v1.3 (el 2π y el factor 125) fueron errores de conversión entre relojes,
y **ninguna produjo excepción, `nan` ni log**. Con dos relojes formales el riesgo se multiplica.

1. **Toda magnitud lleva su reloj en el nombre.** Sin excepciones, ni en variables locales:
   `dt_s`, `dn_ticks`, `dq_usd`; `omega_ang_rad_s`, `omega_ang_rad_tick`, `omega_m_ciclos_tick`;
   `sigma_rel_s`, `sigma_rel_tick`.
2. **Cada conversión ocurre en una única función**, en `constantes_micelio.py`, y esa función
   tiene un test que **invierte la construcción** — el patrón de `periodo_implicito`, que es el
   único que atrapa esta clase de fallo.
3. **Prohibido convertir en el punto de uso.** Si un consumidor necesita una magnitud en otro
   reloj, se publica ya convertida en el bloque compartido.

---

## §3. Ingesta por lotes — el cambio arquitectónico

Es el cambio de mayor superficie del documento. Sin él, nada de §4 y §5 es implementable.

### 3.1 `ENV_DTYPE` deja de ser escalar

La casilla única `P_spot` pasa a un **Ring Buffer de datos de mercado**, espejo del de
actuación (SPSC, un productor en el Motor de Red, consumidores en los Hilos Rápido y Lento).
Cada entrada:

```
aggTradeId  uint64   # 'a' — identidad, para deduplicar y detectar huecos
precio      float64  # 'p'
cantidad    float64  # 'q'  [BTC]
T_trade     float64  # 'T' — hora del trade EN EL EXCHANGE, no de recepción
es_maker    uint8    # 'm' — lado del taker, para el desbalance de flujo
```

Dimensionar el anillo para el peor caso de ráfaga observado, con un factor de holgura de 4×.
La medición del tamaño real de lote es la **primera tarea** del §8: hoy nadie la ha hecho.

### 3.2 `n_ticks` incrementa por transacción

`n_ticks += len(lote)`, no `+= 1`. Con eso ν pasa a ser la tasa real de transacciones y todas
las filas de la tabla de §1.2 se corrigen a la vez.

### 3.3 Deduplicación obligatoria por `aggTradeId`

El sondeo REST devuelve ventanas **solapadas**: sin deduplicar, los mismos trades entran al
filtro varias veces. Eso es exactamente el bug de las 90 correcciones con otro disfraz, y sería
la tercera aparición de la misma familia en tres sesiones.

Regla: se procesa un trade si y solo si `aggTradeId > ultimo_id_procesado`. El contador vive en
el Motor de Red y **no se reinicia entre reconexiones**.

### 3.4 Detección de huecos — un modo de fallo nuevo y necesario

Si `aggTradeId` salta más de lo que trae el lote, **hubo transacciones que no se vieron**. Eso
es distinto del socket mudo: el feed funciona y aun así falta información. Hoy el sistema no
puede detectarlo porque no tiene identidad de trade.

Publicar `n_huecos` y `trades_perdidos` en telemetría. Un hueco invalida el conteo de ticks
para ν en esa ventana, y **debe excluir el tramo del análisis de innovaciones** por la misma
razón que `hay_medicion`: un ρ₁ calculado sobre una serie con huecos silenciosos es basura con
apariencia de dato.

---

## §4. El filtro bajo Δn = 1

**Un tick = un paso de predicción + una corrección.** Es la condición prioritaria de esta orden
de trabajo.

### 4.1 `A_arm` en espacio de ticks — el 2π reaparece

```
omega_ang_rad_tick = 2π · omega_m_ciclos_tick        # [rad/Tick]
```

`omega_m_desde_hz` devuelve **ciclos/tick** (frecuencia ordinaria: `hht.frecuencia_instantanea`
ya dividió por 2π). `A_arm` sale de `s̈ = −ω²s` y necesita **angular**. Misma trampa, nuevo
espacio, mismo remedio: una sola función de conversión y `test_periodo_implicito` reescrito en
ticks.

Con Δn = 1 la matriz es:

```
        ⎡  cos(ω)      sin(ω)/ω    0 ⎤
A_arm = ⎢ -ω·sin(ω)    cos(ω)      0 ⎥        ω ≡ omega_ang_rad_tick
        ⎣  0           0           1 ⎦
```

Control de magnitud: con `ω_m ≈ 1.26e-3 ciclos/tick`, `ω ≈ 7.9e-3 rad/tick`, período ≈ 795
ticks. A 20 transacciones/s son ≈ 40 s — coherente con el ciclo estructural conocido.

**Rama de Taylor obligatoria** para `|ω| < 1e-6`, como en la v1.3.

**Premio:** con Δn = 1 fijo, `A_arm` es **constante entre cambios de régimen** y se cachea; solo
se reconstruye cuando `ω_m` cambia.

⚠ **`x[1]` cambia de unidades.** La velocidad pasa de USD/BTC por segundo a **USD/BTC por
transacción**. Auditar todo consumidor de `x[1]` y toda constante que la escale, `q_base`
incluida. Es el mismo riesgo que motivó la forma afín de la v1.3.

### 4.2 `Q` por tick — una constante física de verdad

```
q_S_tick = (sigma_rel_tick · S)²
sigma_rel_tick = sigma_rel_s / sqrt(nu_por_segundo)      # incrementos iid entre ticks
```

Esa conversión vive en el sitio único del §2.4, con su test inverso.

Lo que se gana no es cosmético: **`Q` por tick es varianza inyectada por transacción**, que es
una propiedad del mercado. `Q` por paso de bucle era varianza por despertar del planificador,
que es una propiedad del sistema operativo. El reloj de ticks saca al planificador del
presupuesto de ruido del filtro.

### 4.3 `R` y la microestructura — no se resuelve aquí, pero se acota

Procesar cada transacción exige que `R` refleje ruido de microestructura: transacciones
consecutivas **no son observaciones independientes** — el rebote bid-ask las correlaciona.
Ese es justo el papel de `r_S,base`, hoy medido mal por 8×.

**Esta orden de trabajo NO lo calibra.** La regla de no inflar `Q` ni `R` con constantes ad-hoc
sigue vigente y la Fase 2 (ALS) sigue siendo su dueña. Lo que hace §4 es darle a la Fase 2 un
objetivo bien definido: `r_S,base` como ruido de observación **por transacción**, que es una
cantidad con significado físico, a diferencia de "por ciclo de control".

Guarda mientras tanto: si `ρ₁` de la innovación de `y0` sube al procesar por transacción
respecto a procesar por paquete, es señal de correlación de microestructura no modelada —
registrarlo, no compensarlo.

### 4.4 Lo que Δn = 1 elimina estructuralmente

No lo parchea: lo vuelve imposible de expresar.

- **Corrección múltiple con la misma medición.** Un tick, un paso, una corrección.
- **El EMD tamizando duplicados** (§5).
- **El planificador dentro de `Q`** (§4.2).
- **La inconsistencia de Δt**, bug bloqueante arrastrado desde la v1.1: se disuelve, porque
  deja de haber Δt en el filtro.

### 4.5 Lo que Δn = 1 NO elimina

Escrito explícitamente para que no se asuma lo contrario:

- `τ_d`, `τ_max` y las 7 guardas de riesgo: reloj de pared.
- El horizonte de Loeper y su condición terminal: reloj de pared.
- El NMPC y los rate limits: reloj de pared.
- **La trampa del 2π**: reaparece en espacio de ticks (§4.1).
- La latencia del feed: sigue siendo un problema, solo que menor de lo que parecía.

---

## §5. El EMD deja de tamizar una escalera

### 5.1 Muestreo por transacción

El buffer del EMD se llena **cada K transacciones**, no cada `PERIODO_MUESTREO_EMD` segundos.
`K` se fija para conservar la cobertura actual de la ventana: con `W = 384` muestras cubriendo
varios ciclos estructurales, y períodos de 795 ticks, `K ≈ ν · PERIODO_MUESTREO_EMD` como punto
de partida — **[CALIBRAR]** contra la distribución real de ν.

Esto satisface a la vez las dos exigencias que hoy están en conflicto: malla uniforme (en
ticks) y ausencia de duplicados. El comentario de `Micelio.py:1455` sobre el jitter del
planificador queda obsoleto y debe reescribirse: en reloj de ticks el jitter no existe.

⚠ **`ω_m` sale ya en ciclos/tick de forma nativa**, sin pasar por `f_hz/ν`. Conservar
`omega_m_desde_hz` para el camino de pared, pero el consumidor del filtro usa la vía directa.
Verificar que ambas coinciden — es un test de consistencia entre relojes, y es barato.

### 5.2 Verificación de la contaminación previa

Test explícito, porque de él depende si el veredicto del A/B tiene que reabrirse:

Tamizar la **misma** serie de mercado por dos vías — muestreada por reloj de pared con
retención de orden cero (como hoy) y muestreada por transacción — y comparar `f_hz` y `C`. Si
difieren de forma material, queda demostrado que `ω_m` estaba contaminado y **el resultado del
A/B de la v1.3 se anula por una segunda causa independiente**.

---

## §6. El reloj de volumen y el factor φ'

### 6.1 El punto de partida es correcto

`Φ = ΣQ · ω_m` (`Micelio.py:1574`) **ya es el acoplamiento volumen–frecuencia**. El reloj de
volumen no es una idea externa al diseño: es la misma intuición física reapareciendo en la capa
de la Sec. 1.4. Eso es lo que justifica tratarlo aquí y no como una variante exótica.

### 6.2 Por qué `Δt · ΣQ` no puede ser el reloj

Dos obstrucciones, ambas fatales para un reloj y ninguna para un peso:

1. **`ΣQ` ya es una integral temporal.** Se acumula desde el nodo de fase (Sec. 6.3).
   Multiplicarla por `Δt` integra el tiempo dos veces: a actividad constante crece como `t²`
   dentro del ciclo.
2. **`ΣQ` se reinicia en cada nodo de fase** (`Q_ref`). Un reloj debe ser monótono; ése
   retrocedería a cero varias veces por hora, y con él toda magnitud derivada.

### 6.3 La formulación correcta — y la coincidencia útil

Lo que hace el trabajo es el **incremento**, no el acumulado:

```
dφ' = dQ                    volumen transado DURANTE el intervalo   [USD]
τ_Q = Σ dQ / ΔQ*            reloj de volumen, monótono               [cuantos]
```

`ΔQ*` es el cuanto de volumen. **[CALIBRAR]**, con una elección recomendada: fijarlo en la
**mediana de volumen por aggTrade**, de modo que a actividad típica el reloj de volumen y el de
ticks avancen a la misma tasa media. Así la única diferencia entre ambos es la **ponderación
por tamaño**, que es exactamente la propiedad que se quiere poner a prueba, y la degradación de
uno al otro es continua en vez de un salto de escala.

Requiere un acumulador monótono nuevo, `Q_acumulado_total`, **distinto de `ΣQ_T̄`** y que no se
reinicia en los nodos. Nombrarlos de forma que no puedan confundirse.

**La observación que unifica:** el peso que se buscaba y el reloj que se buscaba son el mismo
objeto. `dφ' = dQ` sirve de incremento de reloj y de factor de ponderación sin cambiar de
forma. No hacen falta dos construcciones.

### 6.4 φ' como peso — dónde entra

El principio —"para ponderar no uses el transcurso del tiempo, usa el mercado que pasó"— aplica
a todo sitio donde hoy hay un olvido o un decaimiento clavado al reloj de pared:

| Consumidor actual | Hoy | Con φ' |
|---|---|---|
| Refractario del detector de nodos | `5 · PERIODO_HILO_LENTO` | volumen mínimo entre nodos |
| OU de λ (`dt_sim`) | segundos | volumen transado |
| EMA de `f_hz` y `dt_muestreo` | por muestra de reloj | por cuanto de volumen |
| Ventana adaptativa `W_k` del EMD | muestras temporales | cuantos de volumen |
| Burn-in (Sec. 7.1) | ciclos | transacciones o cuantos |

El refractario es el caso más claro: no quieres dos nodos de fase separados por poco *tiempo*,
quieres que estén separados por poco *mercado*. En un mercado muerto, el refractario temporal
dispara nodos espurios; en una ráfaga, los suprime justo cuando son reales.

### 6.5 ⚠ El riesgo de doble conteo en Ω — y la regla que lo cierra

Si `ω_m` se midiera en reloj de volumen, **ya llevaría dentro la información de volumen que
`Φ = ΣQ·ω_m` vuelve a multiplicar.** Ω se construye sobre Φ y Ψ, y Ω determina `κ` y `μ` vía
`Ω_crit` — o sea que el error se propagaría al costo entero del NMPC, escalado al cuadrado
(`κΩ²`, `μΩ²`), sin producir ningún síntoma local.

**Regla:** `Φ`, `Ψ` y `Ω` se calculan **siempre con `ω_m` en reloj de transacciones**, aunque el
filtro corra en reloj de volumen. Publicar `omega_m_ciclos_tick` y `omega_m_ciclos_cuanto` como
campos separados y fijar por escrito qué consumidor lee cuál. Misma disciplina que
`ω_m` / `ω_ang`, y por la misma razón: la conversión implícita en el punto de uso es lo que ha
producido todos los fallos silenciosos de este proyecto.

### 6.6 Alcance de esta versión

El reloj de **transacciones** es el primario y es lo que se implementa. El reloj de **volumen**
se implementa como **selector configurable** sobre la misma maquinaria (`RELOJ ∈ {TICKS,
VOLUMEN}`), no como una rama paralela: si `ΔQ*` se fija según §6.3, la diferencia se reduce a
cómo se cuenta el incremento. φ' como peso (§6.4) se implementa completo en ambos casos, porque
no depende del selector.

Ninguna decisión sobre cuál gana. Se decide con datos, y esos datos no existen todavía.

---

## §7. Lo que no cambia

- Toda la capa de riesgo de la v1.3: las 7 guardas, la máquina de episodios, las 6 compuertas,
  la ruta de cierre, el disparo forzado. Reloj de pared, sin tocar.
- `I_max = 0.50 BTC` y el principio de separación (Sec. B.1 de la v1.3).
- La forma afín de la predicción y el manejo de nodos de fase.
- `c²_vol = k·ω_m·ν`, inmune al error de ν (§1.2).
- Las convenciones del proyecto: español, referencia al PDF, `NOTA DE INTERPRETACION`,
  `DIVERGE DEL PDF`, constantes derivadas, ASCII en texto impreso.

---

## §8. Criterios de aceptación

**§3 — Ingesta**
- [ ] **Medición previa, antes de implementar nada**: tamaño de lote real de `aggTrades` a 1 Hz
      —mediana, p95, máximo— y tasa real de transacciones por segundo. Es el dato que dimensiona
      el anillo y que hoy nadie tiene.
- [ ] Ring Buffer de mercado con los cinco campos; dimensionado a 4× la ráfaga p99.
- [ ] `n_ticks` incrementa por transacción; ν comparada contra la medición anterior.
- [ ] **Test de deduplicación**: inyectar lotes solapados y verificar que ningún `aggTradeId`
      se procesa dos veces.
- [ ] **Test de huecos**: inyectar un salto de `aggTradeId` y verificar que se detecta, se
      publica y **excluye el tramo** del análisis de innovaciones.

**§4 — Filtro**
- [ ] `test_periodo_implicito_ticks`: el período que `A_arm` codifica de verdad coincide con el
      real dentro del 5 %, en espacio de ticks. Es la única defensa contra el 2π.
- [ ] Rama de Taylor: `ω = 0` exacto reproduce `[[1,1],[0,1]]` bit a bit.
- [ ] **Test de invariancia a la tasa**: con `Q(Δt)` bien derivada, N pasos de `Δn=1` y un paso
      de `Δn=N` dan el mismo `x` y la misma `P` dentro de tolerancia numérica. Este test es el
      que demuestra que §1.3 quedó cerrado.
- [ ] Auditoría de consumidores de `x[1]` tras el cambio de unidades, con la lista en el commit.
- [ ] Consistencia entre relojes: `ω_m` nativa en ticks contra `f_hz/ν` coinciden.

**§5 — EMD**
- [ ] Buffer muestreado por transacción; cero duplicados verificado.
- [ ] **Test de contaminación de §5.2**: misma serie por ambas vías, comparación de `f_hz` y
      `C`. Si difieren de forma material, anotar en `CLAUDE.md` que el A/B de la v1.3 queda
      anulado por segunda causa independiente.

**§6 — Volumen**
- [ ] `Q_acumulado_total` monótono, distinto de `ΣQ_T̄`, sin reinicio en nodos.
- [ ] `ΔQ*` fijado como mediana medida de volumen por aggTrade, con el dato en el commit.
- [ ] **Test de la regla §6.5**: verificar que `Φ`, `Ψ` y `Ω` leen `omega_m_ciclos_tick` y no la
      variante de volumen, con el selector en `VOLUMEN`. Un test que lo fuerza a leer la
      equivocada y comprueba que el resultado difiere — para que el test no sea vacuo.
- [ ] φ' aplicado en los cinco consumidores de §6.4; refractario de nodos verificado en mercado
      muerto y en ráfaga.

**Transversal**
- [ ] Nomenclatura del §2.4 aplicada sin excepciones; `grep` de nombres sin sufijo de reloj
      como parte de la suite.
- [ ] Toda conversión entre relojes en una única función, cada una con test inverso.
- [ ] La suite de la v1.3 sigue en verde: **33/33 + los nuevos**.

---

## §9. Fuera de alcance

- **Fase 2 (ALS)** sobre `r_S,base` y `r_EMD`. §4.3 le prepara el terreno y le define el
  objetivo; no la ejecuta.
- **Veredicto del A/B.** Requiere ≥ 24 h continuas y, ahora, que §5.2 haya descartado la
  contaminación de `ω_m`. Antes de eso cualquier veredicto es prematuro por partida doble.
- **Desatascar el WebSocket de futuros.** Sigue pendiente, pero **baja al segundo puesto**: el
  sistema no está escaso de datos, está descartando el 95–98 % de los que ya recibe. El lote se
  arregla sin depender del ISP.
- **Fase 3 (`Ω_crit`)**, Testnet con credenciales, las 30 corridas.
- `C_ON` / `C_OFF`, `ω_m,max`, `ΔS_ref` — siguen `[CALIBRAR]`, y sus datos vendrán de la corrida
  larga posterior a esta orden de trabajo.
- CMA-ES, IMPC, DeepONet.
- Reconciliación documental del PDF: Secs. 3.5 y 8.6.1 superadas por ALS; contradicción
  L1/cuadrática (4.5 contra 6.1, prevalece 6.1); dimensiones del NMPC en 7.5.3; Sec. 2.6
  formalizada como `θ ≡ π/2 (mód π)`; y ahora también la Sec. 4.6, cuyo `n_ticks` hay que
  redefinir explícitamente como transacciones y no como paquetes.
