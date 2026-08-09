# Preregistro — ORDEN_TRABAJO_EJECUCION_4_0

**Fecha de escritura: 2026-08-09.** Se commitea **antes** de ejecutar el §3 en adelante. Su hash
se referencia en el reporte.

**Estado al escribir:** `captura_estacional` lanzada el **2026-08-09 a las 16:30:53 UTC** con
`--dias=21` (§1, lo único que el calendario no perdona). **No se ha mirado ni un dato de ella.**
`captura_v33` sigue viva y en un solo tramo; esta línea **no la toca**.

> §2 de la orden: "con la misma disciplina que la v3.2: enmiendas acreditadas con hash
> antes/después, motivo, y constancia de si había resultado a la vista."

---

## 0. Qué línea es ésta y qué NO decide

Línea **paralela**, no sucesora. La predictiva (v3.2 → v3.3 → v3.4) sigue su curso sin tocarse.
Aquí se mide **cuánto cuesta operar**, que entraba en el sistema como tarifa fija y resultó no
serlo: la tasa de llenado maker medida (**41.3 %**) y el markout (**−5.78 USD/BTC** a 50 s) mueven
`c(u)` de 25.2 a **~37.7** y `H*` de 50 s a **~103 s**.

**No decide** nada de la v3.2, ni toca `Micelio.py`, ni modifica la función de coste del NMPC
(eso es v3.4 y va después de que `c(u, estado)` exista).

⚠ **Fuera de alcance por decisión conjunta, y se anota para que no vuelva:** MPC adversario y
teoría de juegos. **Sin retroacción no hay juego, y un oponente modelado es MENOS robusto que el
minimax que ya está.** La pregunta de "¿aguanta un flujo hostil?" se responde en el §6 como
**prueba de resistencia** —supervivencia, no rendimiento—, no como partida.

---

## 1. Lo que se congela AHORA (§2 de la orden)

### 1.1 Estimadores del §3 y su malla

Por casilla de **hora UTC**, sobre todo el registro:

| símbolo | definición |
|---|---|
| `σ_pared(h)` | volatilidad realizada **por segundo** |
| `σ_tick(h)` | volatilidad realizada **por transacción** |
| `ν(h)` | transacciones por segundo |

- **Estimador de dos escalas o núcleo realizado**, nunca varianza muestral cruda: la curtosis
  limpia es 1179.7 y el ruido de microestructura domina la escala corta.
- **Tres columnas siempre: media, media recortada al 5 %, y mediana.** Un solo dato macro puede
  fabricar una U donde no la hay.
- Malla: 24 casillas horarias. El perfil **semanal** (168 casillas) es **exploratorio** hasta
  ≥ 4 semanas, y así se rotula pase lo que pase.

### 1.2 Regla de colapso del §3.3 — se aplica SIN DISCUSIÓN

| resultado | consecuencia |
|---|---|
| `σ_tick(h)` plana (**rango < 1.3×**) y U sólo en `σ_pared` | **el eje de volatilidad se ELIMINA** del §4. El regresor pasa a ser `ν`. `H*_ticks` se declara constante y se verifica |
| U en **ambos** relojes | el eje se queda: hay estacionalidad de volatilidad real |
| U **sólo** en `σ_tick` | inesperado: **PARAR y reportar**, no seguir |

Predicción mecánica declarada: si `σ_tick` fuera constante, `σ_pared = √ν · σ_tick`, así que una U
en `ν` produce una U en `σ_pared`. Con el rango de `ν` ya medido (3.35 → 16.6 tx/s) el factor
mecánico es `√4.96 ≈ 2.2×`, y la variación intradía típica de `σ` en BTC es de orden 2–3×.
**Entre la mitad y la totalidad de la U podría ser estacionalidad de actividad.**

Entregable: regresión de `log σ_pared(h)` contra `log ν(h)` sobre las 24 casillas, con pendiente
e IC. **Pendiente 0.5 con `σ_tick` plana = la U es enteramente mecánica.**

⚠ Si sale el primer caso, es evidencia independiente para el §7 de la v3.2 (en qué reloj vive el
propagador). **Se reporta como convergencia y NO se usa para modificar la v3.2**, que está
congelada.

### 1.3 Desenlaces del §4.1 y su taxonomía

| desenlace | definición | signo económico |
|---|---|---|
| **LLENADO BENIGNO** | volumen transado al nivel ≥ posición en cola **y el nivel sobrevive** | el bueno |
| **LLENADO POR BARRIDO** | ídem, pero el nivel desaparece en el mismo evento | el caro |
| **SUPERADO** | el mejor bid mejora por encima de tu precio; quedas descolgado | beneficio perdido |
| **CENSURADO** | se alcanza el horizonte (120 s) | sin información |

Se reportan `P(barrido | llenado)` y **los dos markouts por separado**. Esa diferencia es la
medición de selección adversa y es la cantidad que decide si merece la pena retirarse.

⚠ **La censura no se descarta ni se funde.** Modelo: azar en tiempo discreto con **riesgos
competitivos**. Prohibida la regresión logística sobre "llenado sí/no", que fundiría desenlaces de
signo opuesto.

⚠ **Nada de LightGBM ni árboles.** Con ~10³ observaciones dependientes y 4 variables, un modelo
flexible ajusta el ruido de la ventana.

### 1.4 Variable de estado del §4.2

```
T_despeje    = Q_pos / (ritmo de consumo del nivel)   [s]
T_residencia = duración esperada del nivel de precio   [s]
variable primaria = T_despeje / T_residencia           [adimensional]
```

Se congela **antes de ver su poder predictivo**.

⚠ **Requisito de entrada, no apéndice: resolver primero la inconsistencia declarada en el §4.2.**
2 058.1 BTC transados al mejor bid dan ~0.53 BTC/s, con lo que una cola de 22.7 BTC se despejaría
en ~43 s; pero se midió `P(sin llenar a 60 s) = 92.8 %` en ese tercil (**cifra corregida**: el
99.5 % es anterior al arreglo del llenado adverso, commit `0e3b9e0`, y toda cifra de ejecución
previa a ese commit se descarta). **Los dos números siguen sin ser compatibles.**

Primero verificar que ambos salen del **mismo tramo**. Si se sostiene, la lectura es que el cuello
de botella no es despejar la cola sino que **el nivel muere antes**, y entonces manda
`T_residencia`, no `Q₀`.

### 1.5 Eventos programados del §7 — la lista, ANTES de mirar si tienen efecto

| evento | hora | nota |
|---|---|---|
| liquidación de financiación | 00:00 / 08:00 / 16:00 UTC | **verificar contra documentación vigente de Binance**, no por memoria |
| apertura de renta variable EE.UU. | 13:30 UTC | ⚠ **se desplaza con el horario de verano** |
| publicaciones macro (IPC, FOMC) | calendario publicado | son eventos, no estacionalidad |
| apertura/cierre CME | fin de semana | el hueco es estructural |

Regresor: **distancia en minutos al próximo evento programado**. Compite con el regresor cíclico
(`sin`, `cos` de la hora) y **se elige en validación, no por preferencia**. No se usan los dos.

⚠ **Husos por biblioteca, nunca por desplazamiento fijo.** El horario de verano desplaza efectos
una hora a mitad de registro y no se ve.

---

## 2. Compuerta de datos — y NO es la de la v3.2

```
compuerta = min sobre las 168 casillas (hora UTC × día) de los minutos
            observados en esa casilla  ≥  30
```

Es de **cobertura**, no de continuidad: un corte de 40 min no invalida la estacionalidad; que
falte sistemáticamente la franja 03:00–05:00 UTC sí. Se reporta la **matriz 24×7 completa**, no un
porcentaje agregado.

**Si no pasa, no se ejecuta el §3 en adelante.** La única acción admisible es seguir capturando.

---

## 3. Tamaño muestral efectivo — declarado antes de que infle nada

⚠ Inserciones hipotéticas cada 10 s con horizonte 120 s **se solapan 12 veces**. El número
efectivo de observaciones independientes es del orden de las ventanas no solapadas, **no** el de
inserciones.

- **Bootstrap por bloques de DÍA COMPLETO.**
- Se reporta el **número de días**, nunca el de inserciones, como tamaño muestral.

## 4. Contrafactual — límite de validez declarado

Las inserciones son hipotéticas: una orden real **añadiría** a la cola y sería consumida antes que
lo que entre después. Se declara el **tamaño de orden supuesto** y se verifica que es despreciable
frente al `Q₀` mediano. **A tamaño grande el análisis no vale, y se dice** en vez de extrapolar.

## 5. `C_respaldo` (§5) — el término que no existe en ningún documento del proyecto

```
C_maker = p·(comisión_m + markout)  +  (1−p)·C_respaldo
C_taker = comisión_t + media horquilla + impacto
```

Sin `C_respaldo`, **maker gana siempre porque no paga nada por fallar**.

Se mide para cada inserción **no llenada** al alcanzar `τ ∈ {10, 30, 60, 120} s`: qué habría
costado cruzar en ese instante contra el precio de referencia del momento de la inserción.

⚠ **Sesgo de selección declarado ahora:** no llenarse correlaciona con que el precio se fue en
contra, así que `C_respaldo` condicionado a no-llenado **no es** el coste incondicional. Se reporta
la **distribución completa**, no la media, y **SUPERADO se separa de CENSURADO** — el primero es
justamente el caso en que el precio huyó.

⚠ **Markout contra MICROPRECIO de Stoikov, no contra punto medio.** El punto medio arrastra la
deriva por desbalance y atribuiría a selección adversa algo que era predecible. La primera vez se
reportan **ambas columnas** para cuantificar la diferencia.

## 6. Prueba de resistencia (§6) — criterio de éxito declarado

Inyección sintética de **500 BTC de venta** con cuatro calendarios: de golpe, TWAP 5 min, TWAP 1 h,
y ráfagas irregulares. Se mide: convergencia del EAKF (**NIS y `ρ₁` JUNTOS**), si la capa de riesgo
corta y en cuánto, si el impacto permanente estimado arrastra la referencia móvil, si el NMPC
intenta operar **contra** el flujo, y si alguna guarda salta por la razón equivocada.

**El criterio de éxito es SUPERVIVENCIA, no rendimiento:** el sistema no diverge, no acumula
inventario no pedido, y la capa de riesgo actúa. Que además gane dinero no es la pregunta.

⚠ **Dos honestidades obligatorias en el reporte:**
1. **Parámetros congelados ANTES de la inyección.** Si el propagador se reestima sobre el dato que
   contiene la inyección, el test es circular y no mide nada.
2. **`f(v) = v^δ` extrapolada a 500 BTC está dos o tres órdenes de magnitud fuera del rango
   ajustado** (el volumen mediano por transacción es 0.0060 BTC). La trayectoria de precio
   generada es **ilustrativa, no una predicción**, y se reporta el **percentil de tamaño de
   transacción** al que corresponde el ajuste para que se vea cuánto se extrapola.

---

## 7. Limitaciones conocidas ANTES de medir

1. **Perfil semanal exploratorio** hasta ≥ 4 semanas. Con 14 días son 2 observaciones por casilla
   hora×día.
2. **La rotación de cola (10.5×) es cota inferior**: `bookTicker` está estrangulado y las altas y
   bajas entre snapshots se ven netas.
3. **Con L1, cancelación y reposición no son separables**, sólo su neto.
4. **Toda cifra de ejecución anterior al commit `0e3b9e0` se descarta**: el llenado adverso estaba
   clasificado como no llenado y el sesgo ocultaba justamente los llenados malos.
5. **Riesgo de interferencia entre capturas.** Dos procesos contra el mismo endpoint pueden chocar
   con los límites de conexión de Binance. Implementado: esta captura **se detiene sola** si
   `captura_v33` lleva más de 1 800 s sin datos. La v3.2 tiene prioridad.
6. `pyarrow` **no estaba instalado** en este entorno; se instaló (25.0.0) para cumplir el formato
   columnar que el §1 exige. Es una dependencia nueva del entorno de análisis, no de `Micelio.py`.

---

## 8. Registro de enmiendas

Mismo formato que la v3.2: hash antes → después, regla retirada, la que la sustituye, qué la
invalidó, y **constancia de si había resultado a la vista**.

| # | commit antes → después | regla retirada | regla que la sustituye | qué la invalidó | ¿resultado a la vista? |
|---|---|---|---|---|---|
| — | — | (ninguna todavía) | — | — | — |

**Regla: a partir del primer contacto con el conjunto de evaluación del §4, este documento queda
congelado.**
