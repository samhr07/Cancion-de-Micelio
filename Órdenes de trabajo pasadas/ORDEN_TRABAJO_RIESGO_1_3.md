# Orden de Trabajo — Riesgo de Cuenta y Modelo Oscilatorio v1.3

Sucede a `ORDEN_TRABAJO_CALIBRACION_1_2.md`, que queda como registro.

**Precondición.** Este documento asume implementado el **Modo LECTURA** (feed público real de
Mainnet, sin credenciales, ejecución físicamente imposible). Si no está, hacerlo primero: las
Secciones D y E de aquí abajo se calibran sobre datos reales de Mainnet, no sobre Testnet ni
sobre mocks.

**Alcance.**

- **Sección A** — hallazgo bloqueante de granularidad. Leer antes que nada.
- **Sección B** — capa de riesgo de cuenta.
- **Sección C** — máquina de episodios y reaprovisionamiento automático del faucet.
- **Sección D** — matriz `A` oscilatoria (promueve la "Sección E" aparcada en la v1.2).
- **Sección E** — protocolo A/B y criterios de decisión.

**Convenciones heredadas:** comentarios en español, referencia a la sección del PDF al
implementar una fórmula, `# NOTA DE INTERPRETACION:` para huecos del PDF, `# DIVERGE DEL PDF
(Sec X.Y):` para contradicciones, constantes derivadas en `constantes_micelio.py` y nunca como
literales, texto impreso en ASCII.

---

## A. Hallazgo bloqueante — el tope por orden colisiona con la granularidad del instrumento

### A.1 La medición

BTCUSDT perpetuo (USDⓈ-M) tiene `stepSize = minQty = 0.001 BTC`. Con BTC en el orden de
63 000 USD, **un lote mínimo vale ~63 USD**. El mínimo nocional del perp ha estado en 100 USDT.

| minNotional real | Orden legal más pequeña | Con un tope de 100 USD por orden |
|---|---|---|
| 100 USDT | 0.002 BTC ≈ 126 USD | **Ninguna orden es legal.** El bot corre y nunca opera. |
| 10 USDT | 0.001 BTC ≈ 63 USD | `u ∈ {0, 0.001}` — control **binario** |

En el segundo caso el `floor` de `apply_filters` aniquila la salida del NMPC. Un controlador
continuo con 1.6 lotes de rango efectivo no es un controlador: es un interruptor, y no
distinguirías un NMPC bien calibrado de uno roto porque ambos emiten la misma secuencia de
ceros y unos.

### A.2 Tarea previa obligatoria

1. Leer `exchangeInfo` de Testnet **y** de Mainnet y registrar `stepSize`, `minQty`, `tickSize`
   y el nocional mínimo vigente de BTCUSDT. **No hardcodear.** Hoy `Micelio.py:1347` pasa
   `apply_filters(u_c, 1e-5, 1e-5, 10.0, P_spot)` — los cuatro valores son inventados.
2. Publicar los filtros reales en el bloque de hot-reloading y refrescarlos periódicamente.
3. Guarda de arranque `verificar_resolucion_control` (ver B.4).

### A.3 Decisión de dimensionamiento

Se exige **al menos 40 lotes de resolución** entre cero y el tope por orden, para que la
cuantización sea un redondeo y no una decisión binaria:

```
NOCIONAL_MAX_ORDEN / (stepSize · S) >= 40
```

Con `stepSize·S ≈ 63 USD`, eso fija `NOCIONAL_MAX_ORDEN ≈ 3 000 USD` (≈47 lotes). En Testnet
el dinero es de juguete: no hay razón para operar por debajo de la resolución del instrumento.

---

## B. Capa de riesgo de cuenta

### B.1 Principio de separación — `I_max` NO es un parámetro de riesgo

`I_max` es el numerador de γ_0, y `Γ = γ_0` exactamente:

```
γ_0 = I_max / (S · ΔS_ref)          [BTC³/USD²]
```

Bajar `I_max` para "que quepa" en un presupuesto de cuenta baja γ_0 en el mismo factor, baja
`λS²Γ` en el mismo factor, y **deja el freno de singularidad de Loeper (Sec. 4.4.3)
estructuralmente inalcanzable**. Con `I_max = 0.5 BTC` hoy se tiene `λS²Γ ≈ 0.075`; recortar
`I_max` a 0.0016 BTC lo llevaría a ~2.4e-4 contra un umbral de 1. Correrías 30 episodios
validando un sistema al que le desconectaste uno de los dos mecanismos de seguridad, y nada lo
reportaría.

**Regla:** `I_max` se queda en 0.50 BTC. El límite de cuenta es una **abrazadera aguas abajo**,
en el Motor de Red, después del Ring Buffer y antes de firmar. Dos cosas distintas, dos sitios
distintos, y la abrazadera no entra en ninguna fórmula del PDF.

### B.2 Constantes

Van en `constantes_micelio.py`, sección nueva `1.bis LÍMITES DE CUENTA`, con la marca explícita
de que **no son parámetros del modelo**.

```python
# --- Capa de riesgo de cuenta -------------------------------------------------
# NO son parámetros del modelo: no aparecen en ninguna ecuación del PDF. Son una
# abrazadera externa que acota lo que el sistema puede hacerle a la cuenta, y
# debe seguir funcionando aunque el modelo esté completamente equivocado.

NOCIONAL_MAX_ORDEN     = 3_000.0   # [USD] techo por orden individual (ver A.3)
PERDIDA_MAX_EPISODIO   = 3_000.0   # [USD] drawdown de EQUITY desde el inicio
APALANCAMIENTO         = 5         # [adimensional] ver B.3
EPS_HOLGURA_POSICION   = 0.02      # [adimensional] margen de la abrazadera sobre I_max
```

El tope de **posición** se deriva, no se declara:

```python
def nocional_max_posicion(S: float, i_max: float = I_MAX,
                          eps: float = EPS_HOLGURA_POSICION) -> float:
    """Abrazadera de exposición abierta, en USD.

    Se DERIVA de I_max con una holgura pequeña, deliberadamente. Si se declarara
    como constante independiente y quedara por debajo de I_max·S, la abrazadera
    mordería ANTES que las restricciones de caja de la Sec. 6.2 y el NMPC operaría
    permanentemente contra un límite que no sabe que existe: dejarías de medir el
    controlador para medir el clamp.

    La holgura del 2 % hace que sea un backstop contra bugs, no un controlador
    competidor.
    """
    return i_max * S * (1.0 + eps)
```

Con S ≈ 63 000: `I_max·S ≈ 31 500 USD` de nocional, abrazadera ≈ 32 100 USD.

### B.3 Apalancamiento y la guarda que lo ata al tope de pérdida

Mantener `I_max = 0.50 BTC` exige ~31 500 USD de nocional. A 1x eso requiere 31 500 USD de
margen, que probablemente excede lo que entrega el faucet. Se resuelve con apalancamiento,
**pero eso introduce liquidación, que el modelo no contempla en ninguna ecuación.** La única
defensa es que el kill switch dispare siempre antes:

```python
def verificar_cap_antes_de_liquidacion(S, equity_inicio, apalancamiento,
                                       perdida_max, mmr=0.004):
    """El tope de pérdida debe morder ANTES de que la posición sea liquidable.

    Si no se cumple, el kill switch es decorativo: el exchange cierra la posición
    por su cuenta y el episodio termina por una vía que el sistema no controla ni
    registra. `mmr` es el maintenance margin rate del primer tramo de BTCUSDT;
    LEERLO de `leverageBracket`, no asumir el valor por defecto.
    """
    nocional = nocional_max_posicion(S)
    margen   = nocional / apalancamiento
    colchon_liquidacion = margen - nocional * mmr
    if perdida_max >= colchon_liquidacion:
        raise ErrorDimensional(
            f"PERDIDA_MAX_EPISODIO={perdida_max} no muerde antes de la liquidacion "
            f"(colchon={colchon_liquidacion:.0f} USD). Bajar el cap, bajar el "
            f"apalancamiento, o subir el equity inicial."
        )
```

Con 5x: margen ≈ 6 420 USD, colchón ≈ 6 292 USD contra un cap de 3 000 USD. Dispara con
holgura de más de 2×. **Con 20x no se cumpliría** — por eso el apalancamiento es una constante
sujeta a guarda y no una preferencia.

Equity mínimo para abrir episodio:

```
EQUITY_MIN_EPISODIO = nocional_max_posicion(S)/APALANCAMIENTO + PERDIDA_MAX_EPISODIO
                    ≈ 6 420 + 3 000 ≈ 9 420 USD
```

Si el faucet entrega menos, **no arrancar**: registrar el saldo real y reportarlo. Las opciones
entonces son bajar `APALANCAMIENTO` (no sirve, empeora el margen), bajar `PERDIDA_MAX_EPISODIO`,
o bajar `I_max` **a sabiendas** recalculando γ_0 y verificando que `λS²Γ` siga siendo alcanzable
(Sec. B.1). Nunca lo último en silencio.

### B.4 Ubicación y contrato

La capa vive en el **Motor de Red**, como último filtro antes de firmar, y debe poder detener
el sistema **sin cooperación del Hilo Rápido**. Si viviera dentro del NMPC no sería una capa de
seguridad: sería parte de lo que debe vigilar.

Orden de operaciones por orden candidata, estrictamente:

```
1. leer u del Ring Buffer
2. clamp de nocional        -> u <- min(u, NOCIONAL_MAX_ORDEN / S)
3. clamp de posición        -> u <- min(u, (nocional_max_posicion(S)/S) - |inv|)
4. apply_filters (floor a stepSize, minQty, minNotional)
5. RE-validar nocional y posición DESPUÉS de cuantizar
6. si u == 0 tras el floor -> no enviar, contabilizar como "bajo resolución"
7. firmar y enviar
```

El paso 5 no es redundante: `apply_filters` ya arrastraba el bug de validar `minNotional` antes
del floor. Aquí la cuantización solo puede reducir, pero la re-validación es lo que documenta
que se pensó.

Guarda de arranque, junto a las de la Sec. 0.5:

```python
def verificar_resolucion_control(S, step_size, nocional_max_orden, min_lotes=40):
    """A.3: la cuantización debe ser un redondeo, no una decisión binaria."""
    lotes = nocional_max_orden / (step_size * S)
    if lotes < min_lotes:
        raise ErrorDimensional(
            f"Resolucion de control insuficiente: {lotes:.1f} lotes entre 0 y el "
            f"tope por orden (minimo {min_lotes}). El NMPC seria un interruptor."
        )
```

### B.5 Guardas obligatorias

Cada una es un motivo de halt distinto y debe registrarse como tal (`causa_halt`).

| # | Guarda | Condición de disparo |
|---|---|---|
| 1 | **Drawdown de equity** | `equity_actual − equity_inicio_episodio <= −PERDIDA_MAX_EPISODIO` |
| 2 | **Exposición** | `|inventario|·S > nocional_max_posicion(S)` tras reconciliación |
| 3 | **Tasa de órdenes** | `> N_ORDENES_MAX_MIN` en ventana de 60 s |
| 4 | **Rechazos** | `M_RECHAZOS_MAX` consecutivos del exchange |
| 5 | **Estado desconocido** | órdenes en vuelo sin confirmación durante `> T_CONFIRM_MAX` |
| 6 | **Desincronía** | `|inv_local − inv_exchange| > TOL_INV` en la reconciliación periódica |
| 7 | **Reloj** | `|offset| > 500 ms` contra `/fapi/v1/time` |

Notas de implementación:

- **La pérdida se mide sobre EQUITY, no sobre PnL realizado.** Si solo cuentas realizado, una
  posición abierta sangra indefinidamente más allá del tope sin disparar nada. El equity viene
  del `ACCOUNT_UPDATE` del User Data Stream, **no** de un cálculo local: el cálculo local es
  precisamente lo que la guarda 6 existe para desconfiar.
- **Guardas 3 y 4 son protección contra bugs, no contra el modelo.** En Testnet el fallo
  probable no es un NMPC malo: es un bucle que dispara órdenes. Son independientes del Token
  Bucket, que regula pesos de API, no riesgo.
- **Semántica del halt: cerrar y parar, no congelar.** Congelar deja exposición abierta sin
  supervisión. Cerrar tiene además la virtud de ejercitar la ruta de cierre, que de otro modo
  es código que nunca se prueba. Si el cierre falla, reintentar con backoff y escalar la alerta;
  jamás dar el halt por completado sin posición plana confirmada por el exchange.

### B.6 Loeper permanece simulado — y lo que eso implica

`IS_TESTNET` se sustituye por un tri-estado `MODO ∈ {LECTURA, TESTNET, MAINNET}`, porque hoy
ese booleano decide a la vez tres cosas que deben moverse por separado: el generador de precios,
el modelo de `λ` y el modelo de fills.

En `TESTNET`, `λ` sigue viniendo del proceso OU (Sec. 8.1.1). Es correcto: no hay impacto de
mercado sobre un libro simulado. Consecuencias que hay que dejar escritas en el código:

- **Testnet NO puede calibrar `μ_OU`, `θ_OU`, `σ_OU`, `η` ni `λ_min`.** Esas salen del feed de
  solo lectura de Mainnet. Marcar con `# NOTA DE INTERPRETACION:` para que ninguna corrida de
  Testnet las "confirme" — confirmaría el propio proceso OU.
- **Las estadísticas de fill de Testnet no son evidencia sobre calidad de ejecución.** El libro
  es delgado y errático. Testnet valida plomería: firma, filtros, reconciliación, kill switch.
- Con `I_max` intacto (B.1), `λS²Γ` se mantiene en rango realista y el freno de singularidad
  **sí** se ejercita. Es una razón adicional para no tocar `I_max`.

---

## C. Máquina de episodios y faucet automático

### C.1 Por qué el faucet automático necesita compuertas

Un reaprovisionamiento automático sin más convierte el kill switch en un bucle: ante un bug,
el sistema quema los 30 episodios en una hora y produce 30 datasets idénticos e inútiles. El
valor de un halt está en que **alguien mire por qué ocurrió**. Las compuertas de abajo conservan
el automatismo pero preservan esa propiedad.

### C.2 Estados

```
ARRANQUE -> [verificar equity >= EQUITY_MIN_EPISODIO] -> OPERANDO
OPERANDO -> [guarda B.5 dispara]                      -> CERRANDO
CERRANDO -> [posición plana confirmada]               -> CERRADO
CERRADO  -> volcar telemetría, correr diagnóstico, escribir resumen
         -> [compuertas C.3]  --pasa-->  REAPROVISIONANDO -> ARRANQUE
                              --falla--> DETENIDO (requiere intervención)
```

`DETENIDO` es terminal y ruidoso. Nunca se sale de él por software.

### C.3 Compuertas antes de reaprovisionar

Todas deben cumplirse:

| Compuerta | Umbral sugerido | Motivo |
|---|---|---|
| Episodios automáticos consecutivos | `< 5` | Revisión humana obligatoria cada 5 |
| Muestras de telemetría del episodio | `>= MUESTRAS_MIN_EPISODIO` | Un episodio que muere en 30 s es un bug, no una pérdida. No recargar. |
| Causa de halt repetida | no dos veces seguidas la misma | Si repite, es determinista: no lo arregla más dinero |
| Causa de halt | debe ser la guarda 1 (drawdown) | Las guardas 2-7 son fallos de sistema. Recargar sobre ellas es tapar el bug. |
| `diagnostico.py` | ejecutado, resumen escrito | Convierte el bucle en un pipeline de datos y no en una tragamonedas |
| Enfriamiento | `>= 60 s` | Evita bucles apretados si algo falla temprano |

### C.4 El adaptador de faucet

**El faucet de Binance Futures Testnet es una función de la interfaz web, no un endpoint
documentado de la API pública.** No hay que fingir lo contrario. Implementarlo como adaptador
con dos realizaciones:

```python
class Reaprovisionador:
    """Interfaz. La ausencia de una implementación automática NUNCA debe bloquear
    la ruta de halt: el bot para igual, y espera."""
    def solicitar(self, monto_objetivo: float) -> bool: ...

class ReaprovisionadorManual(Reaprovisionador):
    """Por defecto. Emite alerta, escribe el resumen y espera a que el equity de la
    cuenta supere EQUITY_MIN_EPISODIO por sondeo periódico. Sin intervención en el
    proceso: Samuel recarga por la web y el bot lo detecta solo."""

class ReaprovisionadorAutomatico(Reaprovisionador):
    """Opcional. Si se automatiza, hacerlo contra el endpoint que use la web de
    Testnet, aislado en esta clase y solo aquí. Debe:
      - tolerar rate limits del faucet (son estrictos) sin reintentar en bucle
      - fallar de forma limpia degradando a ReaprovisionadorManual
      - no ejecutarse jamás si MODO == MAINNET (aserción dura, no if)
    """
```

La detección por sondeo de equity es la parte que de verdad importa: hace que el modo manual y
el automático sean **el mismo camino de código**, con la única diferencia de quién provoca la
recarga. Eso hace que el modo automático no sea una ruta sin probar.

### C.5 Telemetría por episodio

- `id_episodio` (entero monótono) como campo nuevo en `TELEM_DTYPE`, para que los ~30 episodios
  sean separables en el análisis.
- Al cerrar: volcado del bloque, más un `resumen_episodio_NNN.json` con equity inicial y final,
  `causa_halt`, duración, número de órdenes enviadas / rechazadas / bajo resolución, y las
  estadísticas del diagnóstico (NIS mediana, ρ₁ por canal, n efectivo).
- **Reset limpio al abrir episodio**: ΣQ, racha de burn-in, `S_ref`, ventana del EMD, inventario
  local, contadores de las guardas. El estado que sobrevive a un halt es fuente segura de
  confusión al analizar.
- La alerta va **fuera de la ruta crítica**: su fallo jamás debe tumbar el proceso, mismo
  principio que `log()`.

---

## D. Matriz `A` oscilatoria — promoción de la Sección E de la v1.2

Aparcada en la v1.2 "hasta datos de cuenta demo" porque validarla contra un mock sinusoidal que
nosotros pusimos no probaba nada. Con Modo LECTURA disponible, se levanta el aparcamiento.

Recordatorio de la medición que la motiva (v1.2 §1.6), con 10 % de error en ω — el caso
realista:

| Modelo | NIS | ρ₁ |
|---|---|---|
| Velocidad constante (actual) | 5.24 | +0.841 |
| Sinusoide prescrita, ω −10 % | 2.57 | +0.591 |
| **Oscilador armónico, ω −10 %** | **1.43** | **−0.050** |

### D.1 ⚠ Trampa del 2π — verificar ANTES de escribir nada

`hht.frecuencia_instantanea` divide por 2π: devuelve `f` en **Hz (ciclos/s)**. Y
`omega_m_desde_hz` devuelve `ω_m = f/ν` en **ciclos/tick**. Ambas son frecuencias
**ordinarias**, no angulares.

Pero `A_arm` sale de `s̈ = −ω²s`, cuya solución es `cos(ωt)`: ahí `ω` es **angular**. Usar
`ω_m` directamente introduce un factor 2π ≈ 6.28 de error, es decir un **628 %**. La tabla de
arriba muestra que ya con −30 % el modelo se degrada a ρ₁ = +0.544; con 2π el armónico sería
peor que velocidad constante y concluirías, equivocadamente, que la propuesta no sirve.

### D.2 ⚠ Segunda trampa: `1/Ticks` contra `Δt` en segundos

`ω_m` está en 1/Ticks; el `Δt` del Hilo Rápido está en segundos. `ω_m·Δt` **no es
adimensional**. Con `ω_m ≈ 1.26e-3` y `Δt = 0.01 s` daría `1.26e-5` en lugar del `1.57e-3`
correcto: un factor 125. `A_arm` degeneraría a velocidad constante y el fallo sería
**silencioso** — parecería que el armónico "no aporta".

**Resolución: dos variables de frecuencia distintas para dos usos distintos.** No convertir en
el punto de uso; publicar ambas.

```
ω_m   [1/Ticks]  -> ρ_k (7.3.3) y c²_vol = k·ω_m·ν (4.5).  SIN CAMBIOS.
ω_ang [rad/s] = 2π · f_hz  ->  A_arm, y SOLO A_arm.
```

Añadir `w_ang` a `MICELIO_DTYPE` (escrito por el Hilo Lento bajo el mismo seqlock). Magnitud de
control: con un ciclo de 40 s, `f = 0.025 Hz`, `ω_ang = 0.157 rad/s`, `ω_ang·Δt = 1.57e-3` con
`Δt = 0.01 s`.

### D.3 Forma afín — evita reescribir el significado del estado

La v1.2 formula `A_arm` sobre `s = S − S_ref`. Hacerlo literalmente cambia el significado de
`x[0]` de precio absoluto a desviación, y **todos** los consumidores hay que auditarlos:
`q_S = (σ_rel·S)²`, `γ_0(S)`, el centro de la malla de Loeper, `z₀ = P_spot`, la telemetría, el
NMPC. Es una superficie de bug grande y silenciosa.

Se implementa en su lugar la **forma afín equivalente**, que conserva `x[0] = S` absoluto:

```
x_ref  = [S_ref, 0, 0]ᵀ
x_pred = x_ref + A_arm · (x_k − x_ref)
P_pred = A_arm · P_k · A_armᵀ + Q_k        (sin cambios: un offset no afecta la covarianza)
```

La tercera componente de `x_ref` es 0, así que `R_n` pasa intacto por la fila `[0,0,1]`.

**Beneficio adicional: el manejo de nodos de fase sale gratis.** La v1.2 exige
`s ← s + (S_ref_viejo − S_ref_nuevo)` con `P` sin tocar. En coordenadas absolutas eso es
exactamente **no hacer nada**: `S_ref` salta, el offset se aplica en la predicción siguiente, y
no hay discontinuidad en `x` ni en `P`. Se elimina toda una clase de bugs, incluida la cascada
de innovación espuria → pico de NIS → ventana del EMD a `W_min` → racha de burn-in rota.

Matriz:

```
        ⎡  cos(ωΔt)      sin(ωΔt)/ω    0 ⎤
A_arm = ⎢ -ω·sin(ωΔt)    cos(ωΔt)      0 ⎥          ω ≡ ω_ang [rad/s], Δt [s]
        ⎣  0             0             1 ⎦
```

**Rama de Taylor obligatoria.** `sin(ωΔt)/ω` es 0/0 cuando ω→0. Con `|ωΔt| < 1e-6` usar
`sin(ωΔt)/ω ≈ Δt·(1 − (ωΔt)²/6)` y `−ω·sin(ωΔt) ≈ −ω²Δt`. Sin esto, un régimen sin ciclo
dominante produce `nan` y contamina `P` de forma irreversible.

### D.4 Mitigaciones del riesgo que introduce

`A_arm` asciende `ω_m` de modulador de `Q` a **determinante de la dinámica del estado**. Hoy un
`ω_m` ruidoso solo ensancha la incertidumbre; con el armónico haría ruidosa la transición misma.

1. **EMA sobre `f_hz`** en el Hilo Lento antes de derivar `ω_ang`. La mediana sobre las últimas
   12 muestras ya existe en `hht.py`; la EMA va encima.
2. **Conmutación por concentración espectral, con histéresis.** Definir `C ∈ [0,1]` como la
   fracción de energía de la IMF dominante (ya se calcula en el colapso espectral de la
   Sec. 2.5). Si `C ≥ C_ON` usar `A_arm(ω_ang)`; si `C ≤ C_OFF` usar velocidad constante.
   `C_OFF < C_ON` para matar el chatter. Sugerido: `C_ON = 0.5`, `C_OFF = 0.35` — **[CALIBRAR]**
   sobre la distribución real de `C` en Mainnet, que es un dato que aún no tenemos.
   No escalar `ω` por `C`: eso sesga la frecuencia a la baja incluso cuando el ciclo es nítido.
3. **Registrar la rama activa en telemetría** (`rama_A`, 0 = velocidad constante, 1 = armónico),
   para poder condicionar el análisis del A/B sobre ella.

---

## E. Protocolo A/B y regla de decisión

### E.1 Ejecución en paralelo

Instanciar **dos EAKF sobre el mismo flujo de mediciones**: uno con velocidad constante, otro
con `A_arm`. Solo uno alimenta el control; el otro corre en sombra. Son matrices 3×3 — el costo
es despreciable contra los 2.2 ms/ciclo medidos de Loeper+NMPC.

Es la única comparación honesta: dos corridas distintas verían mercados distintos.

### E.2 Datos

**Modo LECTURA sobre Mainnet, no Testnet.** El libro de Testnet es simulado y su serie de
precios no es un mercado real; la comparación de modelos de precio ahí no significa nada.
Mínimo 24 h continuas, cubriendo al menos un ciclo completo de sesiones asiática / europea /
americana.

### E.3 Regla de decisión

Sobre las innovaciones del canal `y0`, con recorte de transitorio, leyendo **mediana** de NIS
(la media está inflada por cola pesada, medido en la v1.2: media 147 contra mediana 4.69):

| Condición | Decisión |
|---|---|
| `ρ₁(y0)` del armónico `< 0.20` **y** menor que el de velocidad constante | **Adoptar `A_arm`** |
| Ambos con `ρ₁ ≥ 0.20` | No adoptar. El mismatch no es el que creíamos: reabrir el diagnóstico |
| Armónico peor | No adoptar. **Verificar antes D.1 y D.2** — un 2π o un factor 125 se ven exactamente así |

**No usar `y1` para esta compuerta.** Su ρ₁ = 0.87 es artefacto del solapamiento de ventanas del
EMD (retención de orden cero: el Hilo Lento publica a 2 Hz, el Rápido corre a ~89 Hz) y no es
evidencia sobre `A`. Esto resuelve la pregunta que la v1.2 dejó abierta para este documento.

### E.4 Consecuencia sobre la Fase 2 (ALS)

Si se adopta `A_arm` y `ρ₁(y0)` baja de 0.20, **la compuerta de Ljung-Box queda abierta y la
Fase 2 se desbloquea**. Es la vía prevista para atacar `r_S,base` y `r_EMD`, medidos mal por 8×
y 161× respectivamente y candidatos número uno a explicar el NIS residual. La Fase 2 no se
ejecuta en esta orden de trabajo; solo se registra la condición que la habilita.

---

## F. Criterios de aceptación

**Sección A**
- [ ] `exchangeInfo` leído de ambos entornos; filtros publicados en hot-reload; ningún literal.
- [ ] `verificar_resolucion_control` falla con `NOCIONAL_MAX_ORDEN = 100` y pasa con 3 000.

**Sección B**
- [ ] `I_max = 0.50 BTC` sin cambios; `λS²Γ` reportado al arranque y del orden de 0.075.
- [ ] `nocional_max_posicion` derivada de `I_max`, no declarada.
- [ ] `verificar_cap_antes_de_liquidacion` pasa con 5x y **falla con 20x** (test explícito).
- [ ] `MODO` tri-estado sustituye a `IS_TESTNET` en los tres usos.
- [ ] Las 7 guardas implementadas, cada una con su `causa_halt` distinguible.
- [ ] **Test de disparo forzado**: inyectar una pérdida ficticia y verificar la cadena completa
      — detección → cierre → posición plana confirmada → halt → alerta → volcado. Un freno que
      nunca se probó no es un freno. **Este es el criterio de aceptación más importante del
      documento.**
- [ ] Test de la ruta de cierre con el cierre fallando: reintenta, escala, no reporta éxito.

**Sección C**
- [ ] Máquina de estados completa; `DETENIDO` inalcanzable-de-salida por software.
- [ ] Las 6 compuertas de C.3, cada una probada por separado forzando su condición.
- [ ] `ReaprovisionadorAutomatico` con aserción dura contra `MODO == MAINNET`.
- [ ] Modo manual y automático comparten el camino de sondeo de equity.
- [ ] `id_episodio` en telemetría; reset limpio verificado (ΣQ, burn-in, `S_ref`, EMD, inventario).

**Sección D**
- [ ] `w_ang = 2π·f_hz` publicado aparte de `ω_m`; `ω_m` sin cambios en ρ_k y c²_vol.
- [ ] **Test de equivalencia**: con `ω_ang → 0`, `A_arm` reproduce `[[1,Δt],[0,1]]` bit a bit.
- [ ] **Test de la trampa 2π**: con un ciclo sintético de período conocido, verificar que el
      período implícito en `A_arm` coincide con el real dentro del 5 %. Este test es la única
      defensa contra D.1 y D.2, que fallan en silencio.
- [ ] Rama de Taylor probada: `ω = 0` exacto no produce `nan`.
- [ ] Forma afín: nodo de fase con salto de `S_ref` no produce discontinuidad en `x` ni en `P`,
      y no rompe la racha del burn-in (contrastar contra la forma en `s`).
- [ ] Histéresis de conmutación sin chatter; `rama_A` en telemetría.

**Sección E**
- [ ] Dos EAKF en paralelo; sobrecosto medido `< 0.3 ms/ciclo`.
- [ ] Reporte comparativo con NIS **mediana**, ρ₁/ρ₂/ρ₃ por canal, n efectivo y ρ₁ mínimo
      detectable, condicionado por `rama_A`.

---

## G. Fuera de alcance de la v1.3

- **Fase 2 (ALS)** — sigue tras su compuerta. E.4 registra qué la habilita.
- **Fase 3 (`Ω_crit`)** — requiere Testnet operativo, que es lo que esta orden construye.
- **Mainnet con dinero real** — no antes de que las 30 corridas de Testnet estén cerradas y el
  test de disparo forzado se haya ejecutado en cada una.
- **CMA-ES, IMPC, DeepONet** — sin cambios respecto a la v1.2.
- **Actualizaciones documentales pendientes** del PDF: Secs. 3.5 y 8.6.1 superadas por ALS;
  contradicción L1/cuadrática (4.5 contra 6.1, prevalece 6.1); dimensiones del NMPC en 7.5.3;
  Sec. 2.6 formalizada como `θ ≡ π/2 (mód π)`.
- **`PRECIO_REFERENCIA = 45 000`** está obsoleto (BTC ~63 000), lo que deja `K_USD` corto en
  ~40 %. No rompe `γ_Q` porque `SUMA_Q_MAX = 1.3e5` domina el `max()`, pero corregirlo y dejar
  constancia.
