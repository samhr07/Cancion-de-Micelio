# Orden de Trabajo — Coste de Ejecución Dependiente del Estado, Estacionalidad y Resistencia v4.0

Abre una **línea paralela**, no sucede a la v3.2. La línea predictiva (v3.2 → v3.3 integración en
el filtro, v3.4 función de coste del NMPC) sigue su curso sin tocarse. Esta línea mide la otra
mitad del problema: **cuánto cuesta operar**, que hasta ahora entraba en el sistema como tarifa
fija y resultó no serlo.

**Lo que la motiva, en una línea:** la tasa de llenado maker medida (41.3 %) y el markout medido
(−5.78 USD/BTC a 50 s) mueven `c(u)` de 25.2 a ~37.7 USD/BTC y `H*` de 50 s a ~103 s. `c(u)` no
es una constante, y todo lo que este proyecto ha derivado de ella hereda esa dependencia.

**Lo único de este documento que el calendario no perdona:** §1. La estacionalidad exige semanas
de reloj de pared y no se acelera con cómputo ni con mejor método. **Arranca hoy.** Todo lo demás
se ejecuta sobre datos que se van acumulando.

**Lo que NO cambia:** `Micelio.py` no se toca. `captura_v33` no se toca ni se reinicia. La
decisión de la v3.2 sigue congelada y se hace sobre su propio dato.

---

## §0. Registro de errores previos de esta línea

| Origen | Error | Estado |
|---|---|---|
| `captura_v32`, análisis de cola | «NIVEL BARRIDO» clasificado como **no llenado**. Si el nivel se transa entero, la orden se ejecutó | Corregido: 16.2 % → **41.3 %** de llenado; markout 1 s de −1.39 a **−3.13**. ⚠ El sesgo ocultaba precisamente los llenados malos: **toda cifra de ejecución anterior al commit `0e3b9e0` se descarta** |
| `captura_v32`, flujos de cola | «99.1 % cancelación», mal por factor 38: atribución capada por la bajada neta de `B` entre snapshots | Corregido. Y el fondo no era de código: **con L1, cancelación y reposición no son separables**, sólo su neto. Lo identificable es la rotación (**10.5×**, cota inferior por estrangulamiento) |
| Propuesta externa (matriz 3×3×3) | Los tres ejes —sesión, volatilidad, `Q₀`— son **colineales por la propia tesis que los define**. 27 celdas, ~6–8 pobladas, y las de mayor interés vacías por construcción | §4 sustituye la matriz por una variable de estado continua |
| Propuesta externa (regresión logística) | Descarta la censura (20.1 %) y funde desenlaces de **signo económico opuesto** | §4.3: azar en tiempo discreto con riesgos competitivos |
| Todo el proyecto hasta ahora | `c(u)` tratada como tarifa. El término `(1−p)·C_respaldo` **no existe en ningún documento** | §5. Es el término que decide maker contra taker, y nunca se ha medido |

⚠ **Patrón, cuarta aparición:** mezclar dos relojes en un estimador. Ya invalidó la subdifusión de
30–300 s y ya obligó a redefinir `H*` en ticks. §3 lo vigila en la curva en U **antes** de que la
U entre en ninguna decisión.

---

## §1. La captura de estacionalidad — se lanza HOY

### 1.1 Especificación

| campo | valor |
|---|---|
| proceso | **independiente** de `captura_v33`; conexión WebSocket propia |
| streams | `@trade` y `@bookTicker` de BTCUSDT perp |
| duración mínima | **14 días**; objetivo 21 |
| rotación de fichero | **diaria**, con marca UTC en el nombre |
| formato | columnar comprimido (Parquet/Feather). **CSV no** |
| filtro | `tick_valido` en origen, como `captura_v33` |
| husos | **UTC exclusivamente**, en persistencia y en análisis |

Estimación de volumen a `ν ≈ 11 tx/s`: ~13.5 M transacciones y ~25 M actualizaciones de libro en
14 días → **2–4 GB comprimido**. Verificar el espacio libre **antes de lanzar**, no el día 9.

### 1.2 ⚠ La compuerta aquí NO es la de la v3.2

La v3.2 exige **un tramo continuo** de `380·embargo` ticks, porque su estimador necesita bloques
contiguos. La estacionalidad necesita otra cosa: **cobertura**. Un corte de 40 min no la invalida;
que falte sistemáticamente la franja 03:00–05:00 UTC sí.

```
compuerta = min sobre las 168 casillas (hora UTC × día de la semana) de
            los minutos observados en esa casilla  ≥  30
```

Con 14 días son 2 observaciones por casilla de hora×día y 14 por casilla de hora. **Eso basta para
el perfil horario marginal y NO basta para la interacción hora×día.** Declararlo ahora: el perfil
semanal es exploratorio hasta que haya ≥ 4 semanas.

Reportar la **matriz de cobertura completa** en el reporte, no un porcentaje agregado.

### 1.3 Vigilancia diaria

Un resumen automático cada 24 h con: minutos cubiertos por casilla, número y duración de
desconexiones, `ν` media diaria, espacio en disco restante. Si el proceso muere, la pérdida debe
detectarse en horas, no al final.

⚠ **Riesgo real de interferencia:** dos procesos contra el mismo endpoint pueden chocar con los
límites de conexión de Binance. Vigilar reconexiones en **ambas** capturas. Si `captura_v33` se
degrada, la que se para es ésta, no aquélla — la v3.2 tiene prioridad y su compuerta es continua.

---

## §2. Preregistro de esta línea

Commitear `PREREGISTRO_4_0.md` antes de ejecutar §3 en adelante, con la misma disciplina que la
v3.2: enmiendas acreditadas con hash antes/después, motivo, y constancia de si había resultado a
la vista. Lo que se congela:

1. Los estimadores de §3.1 y su malla horaria.
2. La regla de colapso de §3.3 (qué se elimina si `σ_tick` sale plana).
3. La definición de la variable de estado de §4.2, **antes** de ver su poder predictivo.
4. Los desenlaces de §4.1 y su taxonomía.
5. La lista de eventos programados de §7, **antes** de mirar si tienen efecto.

---

## §3. Test A — la curva en U en los dos relojes

### 3.1 Los dos estimadores

Por casilla de hora UTC, sobre todo el registro:

```
σ_pared(h)  =  volatilidad realizada por SEGUNDO
σ_tick(h)   =  volatilidad realizada por TRANSACCIÓN
ν(h)        =  transacciones por segundo
```

⚠ **Ambos con estimador de dos escalas o núcleo realizado**, nunca varianza muestral cruda: la
curtosis limpia es 1179.7 y el ruido de microestructura domina la escala corta. Es el mismo
requisito que el §1.2 de la v3.1.

⚠ **Mediana y media recortada al 5 % además de la media.** Un solo dato macro puede crear una U
donde no la hay. Las tres columnas van al reporte.

### 3.2 La relación que se contrasta

Si la volatilidad por transacción fuera constante:

```
σ_pared(h) = √ν(h) · σ_tick(h)     →     una U en ν produce una U en σ_pared
```

Con el rango de `ν` ya medido en este proyecto (3.35 → 12.4 tx/s, factor 3.7), el factor mecánico
esperable es `√3.7 ≈ 1.92×`. La variación intradía típica de `σ` en BTC es de orden 2–3×.
**Entre la mitad y la totalidad de la U podría ser estacionalidad de actividad.**

Entregable: regresión de `log σ_pared(h)` contra `log ν(h)` sobre las 24 casillas, con su
pendiente e IC. Pendiente 0.5 con `σ_tick` plana = la U es enteramente mecánica.

### 3.3 Regla de colapso — declarada ANTES

| resultado | consecuencia, aplicada sin discusión |
|---|---|
| `σ_tick(h)` plana (rango < 1.3×) y U sólo en `σ_pared` | **el eje de volatilidad se elimina** de §4. El regresor pasa a ser `ν`. `H*_ticks` se declara constante y se verifica |
| U presente en ambos relojes | el eje se queda; hay estacionalidad de volatilidad real |
| U sólo en `σ_tick` | resultado inesperado: **parar y reportar**, no seguir |

⚠ Si sale el primer caso, es evidencia independiente para el §7 de la v3.2 (en qué reloj vive el
propagador). **Se reporta como convergencia, no se usa para modificar la v3.2**, que está
congelada.

---

## §4. Test B — la variable de estado de ejecución

### 4.1 Taxonomía de desenlaces — corregida y ampliada

Sobre inserciones hipotéticas al mejor bid (y simétrico al ask), horizonte 120 s:

| desenlace | definición | signo económico |
|---|---|---|
| **LLENADO BENIGNO** | volumen transado al nivel ≥ posición en cola, **y el nivel sobrevive** | el bueno |
| **LLENADO POR BARRIDO** | ídem, pero el nivel desaparece en el mismo evento | el caro |
| **SUPERADO** | el mejor bid mejora por encima de tu precio; quedas descolgado | beneficio perdido |
| **CENSURADO** | se alcanza el horizonte | sin información |

**La subdivisión del llenado es la medición de selección adversa**, y es nueva: hasta ahora los
dos casos iban juntos. La diferencia de markout entre ambos es la cantidad que decide si merece
la pena retirarse. Reportar `P(barrido | llenado)` y los dos markouts por separado.

### 4.2 La variable de estado — y una inconsistencia que hay que resolver primero

La matriz 3×3×3 se descarta (§0). En su lugar, dos cantidades continuas:

```
T_despeje  =  Q_pos / (ritmo de consumo del nivel)      [s]
T_residencia = duración esperada del nivel de precio     [s]
```

y la variable primaria es su cociente, adimensional.

⚠ **Antes de usarla, resolver esta inconsistencia de los datos actuales:** 2 058.1 BTC transados
al mejor bid sobre el tramo analizado dan ~0.53 BTC/s; una cola de 22.7 BTC se despejaría en
**~43 s**. Pero se midió `P(sin llenar a 60 s) = 99.5 %` en ese tercil. **Los dos números no son
compatibles.**

Primero verificar que ambos salen del **mismo tramo** — si no, la comparación no vale y se rehace.
Si se sostiene, la lectura es que el cuello de botella no es despejar la cola sino que el nivel
muere antes, y entonces `T_residencia` es la variable que manda, no `Q₀`. **Esa reconciliación es
requisito de entrada al §4.3**, no un apéndice.

### 4.3 El modelo

**Azar en tiempo discreto con riesgos competitivos**, no regresión logística: hay censura (20.1 %)
y desenlaces con signo opuesto que no se pueden fundir.

- Variables: `T_despeje`, `T_residencia`, `ν`, y el eje de volatilidad **sólo si §3.3 lo conserva**.
- Estacionalidad como regresor cíclico (`sin`, `cos` de la hora) o como distancia al evento
  programado más cercano (§7). **No las dos**, y se elige en validación.
- **Nada de LightGBM ni árboles.** Con ~10³ observaciones dependientes y 4 variables, un modelo
  flexible ajusta el ruido de la ventana.

### 4.4 ⚠ El tamaño muestral efectivo no es el número de inserciones

Inserciones hipotéticas cada 10 s con horizonte 120 s **se solapan 12 veces**. El número efectivo
de observaciones independientes es del orden de ventanas no solapadas, no de inserciones.
**Bootstrap por bloques de día completo**, y reportar el número de días, no el de inserciones.

### 4.5 Contrafactual

Las inserciones son hipotéticas: tu orden real **añadiría** a la cola y sería consumida antes que
lo que entre después. Válido a tamaño pequeño; declarar el tamaño supuesto y verificar que es
despreciable frente a `Q₀` mediano. A tamaño grande el análisis no vale y hay que decirlo.

---

## §5. Test C — el término que falta: `C_respaldo`

### 5.1 La comparación bien planteada

```
C_maker = p·(comisión_m + markout)  +  (1−p)·C_respaldo
C_taker = comisión_t + media horquilla + impacto
```

`C_respaldo` es lo que cuesta haber esperado y no haberte llenado. **No existe en ningún documento
de este proyecto.** Sin él, maker gana siempre en la tabla porque no paga nada por fallar.

### 5.2 Cómo se mide con lo que ya hay

Para cada inserción **no llenada** al alcanzar `τ ∈ {10, 30, 60, 120} s`: qué habría costado
cruzar en ese instante, contra el precio de referencia del momento de la inserción.

⚠ **Sesgo de selección declarado:** no llenarse correlaciona con que el precio se fue en tu
contra, así que `C_respaldo` condicionado a no-llenado **no es** el coste incondicional. Reportar
la distribución completa, no sólo la media, y separar SUPERADO de CENSURADO — el primero es
justamente el caso en que el precio huyó.

### 5.3 Markout contra microprecio

`markout(h) = ε_fill · (microprecio(t+h) − p_fill)`, con el microprecio de Stoikov, **no** contra
el punto medio. El punto medio arrastra la deriva por desbalance y atribuiría a selección adversa
algo que era predecible. Reportar ambas columnas la primera vez, para cuantificar la diferencia.

### 5.4 El entregable

Una función `c(u, estado)` tabulada o paramétrica, con su incertidumbre, que sustituya a la tarifa
fija en el cálculo de `H*`. **Esto no es un módulo nuevo: es una constante del proyecto que pasa a
estar medida.** Balance de constantes favorable, que es el criterio de siempre.

---

## §6. Test D — pruebas de resistencia

No es teoría de juegos y no modela a ningún adversario. Comprueba que el sistema **sobrevive** a
un flujo grande, que es una pregunta que hoy no tiene respuesta.

### 6.1 Protocolo

Sobre reproducción de datos capturados, inyectar una venta sintética de 500 BTC con cuatro
calendarios: de golpe, TWAP 5 min, TWAP 1 h, y en ráfagas irregulares. En cada uno, medir:

- ¿converge el EAKF o diverge? NIS y `ρ₁` de la innovación **juntos**.
- ¿corta la capa de riesgo, y en cuánto tiempo?
- ¿se dispara el impacto permanente estimado y arrastra la referencia móvil del §5 de la v3.1?
- ¿el NMPC intenta operar **contra** el flujo?
- ¿alguna guarda se dispara por la razón equivocada? (patrón v3.0 §3.3: se vigiló el parámetro
  que no era)

### 6.2 ⚠ Dos honestidades que van en el reporte

1. **Parámetros congelados antes de la inyección.** Si el propagador se reestima sobre el dato que
   contiene la inyección, el test es circular y no mide nada.
2. **`f(v) = v^δ` extrapolada a 500 BTC está muy fuera del rango ajustado** — dos o tres órdenes
   de magnitud sobre el tamaño típico de transacción. La trayectoria de precio generada es
   ilustrativa, no una predicción. **Reportar el percentil del tamaño de transacción al que
   corresponde el ajuste**, para que el lector vea cuánto se está extrapolando.

El criterio de éxito es **supervivencia**, no rendimiento: el sistema no diverge, no acumula
inventario no pedido, y la capa de riesgo actúa. Que además gane dinero no es la pregunta.

---

## §7. Eventos programados, en vez de «al mediodía»

«Los movimientos grandes ocurren al mediodía» no está medido en este proyecto y sería una segunda
afirmación temporal apilada sobre una primera (§3) que puede ser artefacto. Se sustituye por
eventos **exógenos y con hora publicada**, sin riesgo de circularidad:

| evento | hora | verificar |
|---|---|---|
| liquidación de financiación | 00:00 / 08:00 / 16:00 UTC | **contra documentación vigente de Binance**, no por memoria |
| apertura de renta variable EE.UU. | 13:30 UTC | ⚠ **se desplaza con el horario de verano de EE.UU.** |
| publicaciones macro (IPC, FOMC) | calendario publicado | son eventos, no estacionalidad |
| apertura/cierre CME | fin de semana | el hueco es estructural |

Regresor: **distancia en minutos al próximo evento programado**. Compite con el regresor cíclico
de §4.3 y se elige en validación, no por preferencia.

⚠ **El horario de verano es un fallo silencioso clásico.** Persistir en UTC y calcular la hora del
evento con biblioteca de husos, nunca con desplazamiento fijo.

---

## §8. Modos de fallo previstos

| # | fallo | firma | detección |
|---|---|---|---|
| 1 | La captura nueva degrada `captura_v33` | reconexiones en la v33 | vigilancia en ambas; se para la nueva, no la v33 |
| 2 | Compuerta de continuidad aplicada a estacionalidad | se descarta dato útil por un corte de 40 min | §1.2: la compuerta es de **cobertura**, no de continuidad |
| 3 | U fabricada por un solo evento macro | media con U, mediana sin ella | §3.1, tres estimadores en columnas |
| 4 | Ruido de microestructura inflando `σ` a escala corta | `σ_tick` decreciente con el rezago | dos escalas o núcleo realizado |
| 5 | Horario de verano | efectos que se desplazan una hora a mitad de registro | §7, husos por biblioteca |
| 6 | Hora×día confundidas | perfil semanal con 2 observaciones por casilla | §1.2, marcado exploratorio |
| 7 | Solapamiento de inserciones inflando `N` | IC absurdamente estrechos | §4.4, bloques de día |
| 8 | `C_respaldo` leído como incondicional | maker parece peor de lo que es | §5.2, sesgo declarado y desenlaces separados |
| 9 | Markout contra punto medio | selección adversa sobreestimada | §5.3, ambas columnas |
| 10 | Test de resistencia circular | el sistema «resiste» lo que se ajustó a él | §6.2, parámetros congelados |
| 11 | Extrapolación de `f(v)` presentada como predicción | trayectoria de 500 BTC citada como cifra | §6.2, percentil reportado |
| 12 | Disco lleno el día 9 | captura muerta sin aviso | §1.3, resumen diario con espacio libre |
| 13 | Modelo flexible ajustando la ventana | rendimiento excelente dentro, nulo fuera | §4.3, prohibición de árboles |
| 14 | La U usada para decidir antes de §3.3 | eje de volatilidad en §4 sin haber pasado el colapso | §2, orden de congelado |

---

## §9. Criterios de aceptación

**§1 — Captura**
- [ ] Lanzada, con hora UTC de arranque en el reporte.
- [ ] Espacio en disco verificado antes de lanzar.
- [ ] Resumen diario funcionando; primera salida pegada al reporte.
- [ ] Matriz de cobertura 24×7 completa, no un porcentaje.

**§3 — Dos relojes**
- [ ] `σ_pared`, `σ_tick` y `ν` por hora UTC, con estimador robusto al ruido.
- [ ] Media, media recortada y mediana, en columnas separadas.
- [ ] Pendiente de `log σ_pared` contra `log ν` con IC.
- [ ] Regla de colapso de §3.3 **aplicada sin discusión** al resultado que salga.

**§4 — Ejecución**
- [ ] Inconsistencia de §4.2 resuelta **antes** de ajustar nada.
- [ ] Taxonomía de cuatro desenlaces, con llenado subdividido y sus dos markouts.
- [ ] Azar con riesgos competitivos; censura tratada, no descartada.
- [ ] Bootstrap por bloques de día; **número de días** reportado, no de inserciones.
- [ ] Tamaño de orden supuesto declarado y comparado con `Q₀` mediano.

**§5 — Coste**
- [ ] `C_respaldo` medida a los cuatro horizontes, con distribución completa.
- [ ] SUPERADO y CENSURADO separados.
- [ ] Markout contra microprecio, con la columna de punto medio como comparación.
- [ ] `c(u, estado)` entregada con incertidumbre.

**§6 — Resistencia**
- [ ] Cuatro calendarios, con parámetros congelados previamente.
- [ ] NIS **y** `ρ₁` juntos en cada caso.
- [ ] Percentil de tamaño de transacción del ajuste de `f(v)` reportado.

**Transversal**
- [ ] `PREREGISTRO_4_0.md` commiteado antes de §3; hash en el reporte.
- [ ] `Micelio.py` sin cambios (`git diff --stat`).
- [ ] `captura_v33` intacta y en un solo tramo.
- [ ] Suite en verde, incluidos los nuevos.
- [ ] UTC en toda la persistencia y todo el análisis.

---

## §10. Fuera de alcance

- **Modificar la función de coste del NMPC** (`R_maker`, penalizaciones de ejecución). Toca
  `Micelio.py`; es v3.4 y va después de que `c(u, estado)` exista.
- **MPC adversario / teoría de juegos.** Retirado por decisión conjunta: sin retroacción no hay
  juego, y un oponente modelado es **menos** robusto que el minimax que ya está.
- **Ingerir `@depth`.** El horquillado cayó a 1.2 % y el argumento de datos se debilitó. Sigue
  siendo la única vía para la rotación real de cola, pero no es prioritario.
- **Cualquier cosa que toque la v3.2.** Su decisión está congelada y se hace sobre su dato.
- **Perfil semanal como confirmatorio.** Exploratorio hasta ≥ 4 semanas.
- **Terminal Bloomberg.** Si se retoma, es para regresores exógenos (base CME, financiación entre
  venues, DXY), nunca para microestructura, y **con la cuestión de licencia resuelta con el
  representante de Bloomberg, no con TI.**
- Fase 2 (ALS), Fase 3, Testnet con credenciales, las 30 corridas.
