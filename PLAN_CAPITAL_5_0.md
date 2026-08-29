# Plan de Despliegue de Capital v5.0 — condicional

**Este documento no es ejecutable hoy.** Describe qué hacer *si* las pruebas pendientes salen
favorables. Se redacta ahora por la misma razón que se redactan los preregistros: escribir las
reglas antes de que haya resultados es lo que impide ajustarlas a ellos.

**Condiciones de activación.** Ninguna parte de este plan se ejecuta hasta que las cuatro se
cumplan:

| # | condición | estado |
|---|---|---|
| 1 | v3.2 cierra con un modelo de impacto seleccionado | captura en curso |
| 2 | Paso 3 de la regla de decisión pasa: `E[Δp]` supera `c(u)` con margen | **no tocado** |
| 3 | `c(u, estado)` medida, con `C_respaldo` incluida (v4.0 §5) | pendiente |
| 4 | Demo con dinero ficticio, semanas, con desconexiones y llenados parciales | pendiente |

Si la 2 falla, este documento se archiva sin ejecutarse. Ése es el desenlace declarado en la
enmienda 9 del preregistro 3.2 y **no se renegocia aquí**.

**Supuestos numéricos** (verificar antes de usar): precio BTC ≈ 96 000 USD, `minQty` = `stepSize`
= 0.001 BTC, notional mínimo 50 USDT, `σ` diaria de BTC ≈ 3 %, `c` ida y vuelta ≈ 37.65 USD/BTC,
inyección 300 000 COP/mes ≈ 80 USD.

---

## §1. El problema de la secretaria: por qué no aplica

Se documenta el rechazo con su motivo, no sólo el veredicto.

El problema de la secretaria resuelve: **una** elección, **irrevocable**, **sin recuerdo**, sobre
`n` candidatos de los que sólo se observa el **rango relativo**, maximizando la probabilidad de
acertar **el mejor**. La regla `1/e` (37 %, no 2/3) sale de esa estructura.

Ninguna de las cinco premisas se cumple aquí:

| premisa | tu caso |
|---|---|
| una sola elección | operas cientos de veces al año |
| irrevocable | puedes cerrar y reabrir |
| sin recuerdo | la oportunidad rechazada vuelve |
| sólo rango relativo | observas `α` y `c` **en dólares** |
| objetivo: acertar el mejor | objetivo: valor esperado positivo repetido |

La cuarta es la más costosa: renunciar a información cardinal para usar sólo ordinal es tirar
justo lo que tanto trabajo ha costado medir.

**Y la respuesta que buscas ya existe en el plan.** «¿Cuándo invertir?» está contestado por la
banda muerta del NMPC: se opera cuando `|α| > c_efectivo`, y no se opera cuando no. Ésa **es** la
regla de umbral óptima bajo el objetivo real. Añadir un muestreo de 2/3 encima no la mejora: la
retrasa.

### 1.1 Dónde sí hay un problema de parada óptima

Tu intuición apunta a algo real, pero está en otro sitio: **cuándo dejar de esperar un llenado
maker y cruzar**. Ahí sí hay parada óptima con información que llega — observas la cola
avanzando, el precio moviéndose, `C_respaldo` creciendo — y decides cancelar. La v4.0 §5 lo
resuelve hoy en **lazo abierto** (elegir `τ*` de una tabla). La versión de lazo cerrado
—cancelar cuando el estado se deteriora— es estrictamente mejor y es el sucesor natural.

Ése es el sitio donde vale la pena invertir el esfuerzo de parada óptima.

### 1.2 Y la fase de muestreo ya está en el plan

Se llama demo. No dura 2/3 de nada: dura **hasta que los estimadores convergen**, que es un
criterio medible y no una fracción arbitraria. Observar el mercado es gratis; sólo la parte de
ejecución —tasas de llenado, markout real— exige poner órdenes, y para eso está la cuenta demo.

---

## §2. Granularidad: el primer año no hay nada que dimensionar

```
orden mínima = max(0.001 BTC, 50 USDT/precio) = 0.001 BTC ≈ 96 USD
```

El filtro que manda es `minQty`, **no** el notional de 50 USDT. A 1× de apalancamiento, el número
de lotes disponibles es `floor(W/96)`:

| mes | capital | lotes a 1× | granularidad |
|---|---|---|---|
| 1 | 80 | **0** | no se puede operar |
| 3 | 240 | 2 | 50 % |
| 6 | 480 | 5 | 20 % |
| 12 | 960 | 10 | 10 % |
| 24 | 1 920 | 20 | 5 % |

**Consecuencia que simplifica todo el plan:** una fórmula continua de dimensionamiento no
significa nada sobre una rejilla de 2 o 5 puntos. El error de cuantización supera cualquier
refinamiento que Merton pueda aportar.

| fase | capital | regla de tamaño |
|---|---|---|
| **F1** | < 240 USD | **no operar.** Acumular inyecciones |
| **F2** | 240 – 960 | **un lote o ninguno.** Sólo la banda muerta decide |
| **F3** | > 960 | Merton con techo (§3) empieza a tener sentido |

Todo el aparato de Kelly/Merton queda **dormido el primer año**. Lo único que tiene que
funcionar en F2 es el umbral `|α| > c_efectivo`. Menos piezas, menos cosas que pueden fallar, y
una prueba mucho más limpia de si el modelo acierta.

⚠ **No apalancarse para alcanzar el mínimo.** En el mes 1 harían falta ~1.2× sólo para poder
poner una orden: el tamaño te lo impondría el exchange, no tu criterio de riesgo, y el techo del
§3 quedaría anulado antes de aplicarse. Si no cabe una orden a 1×, no hay capital suficiente.

---

## §3. Dimensionamiento en F3 — dos cosas separadas

```
u* = min(  α·W/(γ·σ_h²)  ,   0.02·W/(3·σ_H)  )
            └── señal ──┘     └── techo ────┘
```

**No se deriva `γ` del techo.** Hacerlo cancela `α` algebraicamente y convierte toda la fórmula
en una constante:

```
γ = 3α/(0.02σ)   ⟹   u* = 0.02W/(3σ)   ← α desaparece
```

`γ` es una preferencia tuya, se elige, y se documenta como elección. El techo hace su trabajo
aparte.

### 3.1 ⚠ El horizonte del techo no es el de tenencia

`σ_H` debe incluir el **peor caso de salida**, no la tenencia prevista. Si piensas salir en
30 min pero en malas condiciones tardas 2 h, el riesgo es el de 2 h:

| `σ` usada | horizonte | techo resultante |
|---|---|---|
| 3 % diaria | 1 día | 0.22·W ← demasiado restrictivo |
| 0.43 % | 30 min (previsto) | 1.55·W ← **falsamente permisivo** |
| 0.87 % | 2 h (con retraso de salida) | **0.77·W** |

La diferencia entre 1.55 y 0.77 es exactamente `C_respaldo` reapareciendo en la capa de riesgo.
**Usar la fila de 30 min autorizaría el doble de exposición de la que es honesta.**

Con el techo en 0.77·W, `u*` queda por debajo de 1× — o sea, **el apalancamiento nunca está
justificado por esta cuenta**, y lo que limita en F2 y F3 sigue siendo el tamaño del lote.

### 3.2 Kelly fraccional

Si en F3 se adopta la rama de señal: **un cuarto de Kelly**, no Kelly completo. No por timidez —
`μ` estimado con error hace que Kelly completo sobreapueste sistemáticamente, y con este número
de operaciones el error en `μ` es grande.

---

## §4. Presupuesto de operaciones

```
coste ida y vuelta = 0.039 % del nocional
```

Para que las comisiones no superen el 20 % del capital al año:

| apalancamiento | idas y vueltas/año |
|---|---|
| 1× | ~500 |
| 3× | ~170 |

A 1× son ~2 al día. Es compatible con horizontes de 10 min – 1 h **si eres selectivo**. Es un
límite duro que **no depende de si la señal funciona**: si el diseño quiere operar más a menudo,
las comisiones se lo comen antes de que la señal importe.

Instrumentar un contador: si el ritmo proyectado supera el presupuesto, el umbral de la banda
muerta sube hasta que encaje.

---

## §5. El P&L no puede validar nada, y eso cambia el diseño

Años necesarios para distinguir un Sharpe `S` de cero con potencia 0.8:

| Sharpe verdadero | años |
|---|---|
| 0.5 | ~31 |
| 1.0 | ~8 |
| 1.5 | ~3.5 |
| 2.0 | ~2 |

**Ni siquiera una estrategia excelente se valida por su P&L en un plazo útil.** Y depende del
tiempo transcurrido, no del número de operaciones: operar más no acelera la respuesta.

Por tanto la fase con dinero real **no es una prueba de rentabilidad**. Es una prueba de si las
predicciones del modelo sobre cantidades observables se cumplen fuera de muestra:

| predicción | se contrasta contra | converge en |
|---|---|---|
| `P(llenado \| estado)` | llenados reales | ~decenas de órdenes |
| markout esperado | markout realizado | ~decenas de llenados |
| `C_respaldo(τ)` | coste real de cruzar tarde | ~decenas de no-llenados |
| impacto de nuestras órdenes | movimiento tras ejecutar | ~decenas |
| deslizamiento contra el previsto | precio de llenado | inmediato |

Cada orden da una observación de cada una. Eso converge en **semanas**, no en años.

**Regla:** el sistema se considera validado si estas cinco coinciden con lo predicho dentro de su
incertidumbre. Si coinciden y el P&L es negativo, el modelo es correcto y no hay `α` — resultado
limpio. Si no coinciden, el modelo está mal, **independientemente de lo que diga el P&L**. Un
P&L positivo con predicciones fallidas es suerte, y se trata como tal.

---

## §6. Frenos y criterio de abandono

### 6.1 Frenos

| freno | umbral | acción |
|---|---|---|
| diario | 8 % del saldo | apagado hasta el día siguiente |
| mensual | 25 % del saldo | **modo lectura** hasta la siguiente inyección y recalibración |
| interruptor manual | — | apagado inmediato, **sin depender de que el bot decida apagarse** |

Con una posición de un lote sobre 240 USD, el freno diario salta con un movimiento del ~6.7 % en
BTC: poco frecuente, bien calibrado. A 3× saltaría con un 2.2 %, que ocurre varias veces al mes —
otro argumento contra el apalancamiento. **Si el freno salta cada semana no es un freno, es el
modo normal de operación**, y la tentación será ensancharlo.

### 6.2 ⚠ La inyección esconde la pérdida

> *«Si pierde el 25 %, la inyección del mes siguiente recupera el 100 %.»*

Aritméticamente cierto, y es el razonamiento por el que un sistema perdedor se sostiene años. El
saldo nunca parece roto porque se repone; la pérdida acumulada crece sin que nada la muestre.

**El P&L acumulado se lleva en una cuenta separada del saldo.** El criterio de abandono se
escribe **en pesos, antes de desplegar**, se commitea con hash, y no se renegocia:

```
si  pérdida acumulada > ____ COP   →   apagado definitivo
```

*(El número lo pones tú. Sugerencia de orden de magnitud: seis inyecciones, 1 800 000 COP. Lo
que no vale es dejarlo en blanco «hasta ver cómo va».)*

---

## §7. Fases de despliegue

| fase | duración | capital | qué se prueba |
|---|---|---|---|
| **D1** demo | ≥ 6 semanas | ficticio | supervivencia operativa: desconexiones, errores de API, llenados parciales, caídas del exchange |
| **D2** real mínimo | ≥ 3 meses | 240 USD, 1 lote | las cinco predicciones del §5 |
| **D3** F2 completa | ≥ 6 meses | hasta 960 | banda muerta con varios lotes |
| **D4** F3 | — | > 960 | Merton con techo |

Cada fase exige que las predicciones del §5 se cumplan en la anterior. **No se salta ninguna por
buenos resultados**: un buen resultado en D2 es exactamente el momento de mayor tentación y de
menor evidencia.

---

## §8. Modos de fallo

| # | fallo | firma | detección |
|---|---|---|---|
| 1 | `minQty` cambia | órdenes rechazadas | consultar `/fapi/v1/exchangeInfo` al arrancar, no cachear |
| 2 | Apalancarse para alcanzar el mínimo | el techo del §3 nunca se aplica | `assert` de que `u* ≤ W/precio` |
| 3 | `γ` derivado del techo | el tamaño no responde a `α` | test: dos `α` distintas deben dar `u*` distintos |
| 4 | `σ` del horizonte previsto, no del peor caso | exposición al doble de lo honesto | §3.1, la fila de 2 h |
| 5 | Presupuesto de operaciones excedido | comisiones > 20 % anual | contador con umbral adaptativo |
| 6 | P&L usado como validación | se concluye antes de tiempo | §5, las cinco predicciones |
| 7 | Pérdida acumulada oculta por inyecciones | saldo sano, pérdida creciente | contabilidad separada |
| 8 | Freno ensanchado tras saltar | umbral distinto del preregistrado | umbral en fichero versionado |
| 9 | Liquidación en lugar de drawdown | pérdida ≫ freno diario | sin apalancamiento no aplica; con él, vigilar margen |
| 10 | Fase saltada por buenos resultados | D3 sin completar D2 | requisito explícito |

---

## §9. Criterios de aceptación antes de D2

- [ ] Las cuatro condiciones de activación cumplidas y documentadas.
- [ ] `minQty`, `stepSize`, notional mínimo y **comisiones reales de la cuenta** leídos del
      exchange, no supuestos.
- [ ] `PREREGISTRO_5_0.md` con: `γ` elegida y justificada, techo, frenos, presupuesto de
      operaciones y **criterio de abandono en COP**. Commiteado con hash.
- [ ] Contabilidad de P&L acumulado separada del saldo, con test.
- [ ] Interruptor manual probado: apaga el bot sin que el bot participe.
- [ ] D1 completada ≥ 6 semanas incluyendo al menos una desconexión real.
- [ ] Las cinco predicciones del §5 instrumentadas y comparándose automáticamente.
- [ ] Samuel puede explicar por qué el bot abrió una posición concreta, elegida al azar del
      registro de D1.

Ese último punto no es formalidad. Sin él no se distingue un fallo de una mala racha, y ésa es la
diferencia entre parar a tiempo y no parar.

---

## §10. Fuera de alcance

- Todo lo anterior a la condición 2. **Este documento no adelanta que exista `α`.**
- Parada óptima en lazo cerrado para `τ*` (§1.1): sucesor natural de la v4.0 §5, no de esto.
- Apalancamiento por encima de 1×: no lo justifica ninguna cuenta de este documento.
- Aumentar la inyección mensual: no cambia ninguna conclusión; el límite es la rejilla de lotes y
  el tiempo de validación, no el capital.
- Optimizar rentabilidad en D2. El objetivo de D2 es información, no dinero.
