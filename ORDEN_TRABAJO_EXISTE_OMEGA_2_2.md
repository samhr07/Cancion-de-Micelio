# Orden de Trabajo — ¿Existe ω_m? v2.2

Sucede a `ORDEN_TRABAJO_OMEGA_2_1.md`, ejecutada parcialmente (§2 y §4 completos, mediciones
del §3, suite 51/51). Nada de lo construido se deshace.

**Este documento no implementa nada del modelo.** Es un experimento con una sola pregunta:

> ¿`ω_m` es una propiedad del mercado, o un artefacto de la ventana de análisis?

De la respuesta depende si `ω_m`, `Φ`, `Ψ`, `Ω`, `Ω_crit`, `A_arm` y el detector de nodos de
fase siguen en el diseño. Por eso el §7 de la v2.1 (correlación `ω_m`–`ΔS`) queda **bloqueado**:
mide si `ω_m` covaría con la amplitud, y todavía no sabemos si `ω_m` es algo.

---

## §0. Registro de errores previos

| Origen | Error | Estado |
|---|---|---|
| v2.1 §3.4 (documento) | Se fijó "fracción de duplicados" como criterio para elegir observable. Era un **proxy**, y el proxy engañaba: el precio de transacción se mueve más que `mid` porque **rebota entre bid y ask** — movimiento sin información, el mismo ρ₁ = −0.164 de la v2.0. | **Criterio equivocado.** Medido contra lo que de verdad importa, `mid` gana: fracción instrumental 62 % contra 195 %, validez 25.0 % contra 12.5 %. §3 **se reabre** (ver §4) |
| v2.1 §3 (documento) | Se justificó `mid` por "se actualiza cuando se mueve cualquiera de los dos lados" | Cierto para la tasa de **mensajes** (674 msg/s), falso para el movimiento: solo el **1.1 %** mueve el mid, o sea 7.7 cambios/s contra 102 tx/s. Con spread de un tick, el mid queda clavado |
| v2.0 §5.2 (medición) | El 44.1 % de duplicados se leyó como propiedad del observable | Es propiedad de **ν**: con ν = 102 y K = 51, el mismo observable da 9.1 %. Confirma el §1.4 de la v2.1 |

Los tres son el mismo error de método: **un indicador indirecto usado como si fuera la
cantidad**. Es la razón de que este documento mida directamente lo que quiere saber.

---

## §1. Estado de la pregunta

### 1.1 La evidencia que la plantea

| | precio de transacción | `mid` |
|---|---|---|
| fracción instrumental `σ_instr²/σ_total²` | **195 %** | 62 % |
| `omega_valida` | **12.5 %** | 25.0 % |

`σ_genuina² < 0` significa: **cambiar la rejilla mueve `ω_m` más que cambiar de tramo de
mercado**. Y tres cuartas partes del tiempo el período medido cae fuera de la banda de
resolubilidad — es decir, no es una medida, y hasta la v2.1 se usaba igual.

### 1.2 Por qué la EMD no puede responder esto sola

**La EMD no tiene hipótesis nula.** Siempre devuelve modos: con ruido blanco devuelve IMFs, con
un paseo aleatorio devuelve un "ciclo dominante" con su frecuencia instantánea y su `C`, todo
bien formado. Estructuralmente no puede emitir "aquí no hay ciclo".

Y hay un resultado establecido —EMD como banco de filtros sobre ruido fraccionario, línea
Flandrin/Rilling/Gonçalvès y Wu/Huang— según el cual sobre paseos aleatorios la EMD se comporta
como un banco diádico y **el modo dominante lo fija la ventana, no la señal**. Un paseo
aleatorio tiene espectro 1/f²: energía creciente hacia baja frecuencia, sin escala
característica. El modo dominante acaba siendo el más largo que la ventana sostenga.

`σ_instr² > σ_total²` es exactamente esa firma.

### 1.3 La explicación alternativa, igual de compatible

La banda medida es **[2.5 s, 64 s]**: estrecha y toda en el extremo rápido. Si el ciclo
estructural vive en decenas de minutos —plausible para flujo de órdenes y regímenes de
mercado— la configuración actual **no puede verlo por construcción**, y el 87.5 % de rechazos
sería la guarda funcionando bien.

Con `Δτ = 0.5 s`, ver un ciclo de 30 min exige `W ≥ 10 800`. Estamos en `W = 384`: factor 28×.

**Dos lecturas, no distinguibles con lo medido:**

- **(A)** No hay escala característica. `ω_m` no existe como cantidad física.
- **(B)** Existe, pero fuera de la ventana donde llevamos tres sesiones mirando.

---

## §2. Experimento 0 — el histograma que falta (hacer primero, cuesta una tarde)

El reporte da 12.5 % de validez pero **no dice por qué lado falla**. Es un histograma sobre
capturas que ya existen, y mueve la prioridad de todo lo demás.

- Distribución del período medido, con las cotas de banda superpuestas.
- Fracción de rechazos **por arriba** (`T > W·Δτ/ciclos_min`) contra **por abajo**
  (`T < m_min·Δτ`), por observable.
- Para los rechazos por arriba: cociente `T / (W·Δτ)`. Si se concentra cerca de un valor fijo,
  es ya un indicio fuerte de que el período lo fija la ventana.

**Lectura:** rechazo dominante por arriba apunta a (B) o a que la EMD devuelve la tendencia;
por abajo, a que la energía dominante está en microestructura.

---

## §3. Experimento 1 — nulo por sustitutos

Es el test decisivo y es esencialmente gratis: no necesita datos nuevos.

### 3.1 Las tres vías

Mismo esquema de tres vías que funcionó en el §5.2 de la v2.0, y por la misma razón: con dos no
se puede separar la causa.

| vía | qué conserva | qué destruye |
|---|---|---|
| **REAL** | todo | — |
| **IAAFT** (fase aleatorizada) | espectro de potencia y distribución marginal | estructura de fase / no lineal |
| **BARAJADO** (incrementos permutados en el tiempo) | distribución marginal exacta —curtosis 207, retícula de `tickSize`, masa en cero— y los tiempos de llegada reales | toda estructura temporal |

El barajado es el nulo más limpio posible aquí: al permutar los **incrementos observados** se
conserva por construcción todo lo que hace rara a esta serie, y solo se destruye el orden.

### 3.2 Qué comparar

Pasar las tres vías por la cadena `ω_m` **completa y sin modificar**, y comparar:

- distribución de `ω_m` (no solo la media: la anchura es el objeto de estudio)
- `σ_instr²/σ_total²`
- tasa de `omega_valida`
- distribución de `C`

### 3.3 Lectura

| resultado | lectura |
|---|---|
| REAL ≈ BARAJADO | **La cadena no extrae información temporal.** `ω_m` es salida del algoritmo, no del mercado |
| REAL ≈ IAAFT ≠ BARAJADO | Lo que se extrae es solo el espectro 1/f². Hay memoria, no hay ciclo |
| REAL ≠ ambos | Hay estructura genuina. Pasar al §4 para saber a qué escala |

---

## §4. Experimento 2 — barrido multiescala, y reapertura del §3

### 4.1 Barrer Δτ, no W

La banda es `[m_min·Δτ, W·Δτ/ciclos_min]`. Subir `W` solo mueve el techo y cuesta cómputo
cuadrático en el tamizado; subir `Δτ` (o sea `K`) **mueve la banda entera** y es barato.

| config | Δτ | banda con W = 384 |
|---|---|---|
| C0 (actual) | 0.5 s | 2.5 – 64 s |
| C1 | 2 s | 10 – 256 s |
| C2 | 8 s | 40 – 1 024 s |
| C3 | 30 s | 150 – 3 840 s |
| C4 | 120 s | 600 – 15 360 s |

Cinco configuraciones cubren de 2.5 s a 4.3 h.

### 4.2 El estadístico que discrimina

Para cada configuración, `T_dominante / (W·Δτ)`.

- **Constante entre configuraciones** ⇒ el período lo fija la ventana ⇒ **(A)**.
- **Alguna configuración donde `T` deja de seguir a la ventana y se ancla** ⇒ **(B)**, y esa
  configuración marca dónde vive la escala.

Reportar además, por configuración y por observable, `σ_instr²/σ_total²` y la tasa de
`omega_valida`. Puede haber una escala donde la fracción instrumental caiga por debajo del 50 %
aunque en C0 sea del 195 %.

### 4.3 §3 se reabre aquí, con el criterio corregido

El observable no se elige por duplicados (§0). Se elige por **`σ_instr²/σ_total²` y tasa de
validez**, y se evalúa **por configuración**, porque la respuesta puede depender de la escala:
`mid` es más pegajoso a escala fina, pero a `Δτ = 30 s` la pegajosidad de microestructura es
irrelevante y podría ganar por otra razón.

Candidatos a evaluar en las cinco configuraciones: precio de transacción, `mid`, y **precio
medio ponderado por volumen sobre el intervalo Δτ** — que a escalas gruesas es el observable
natural y que además usa el volumen, coherente con `dφ' = dQ` de la v2.0.

### 4.4 Requisito de datos — el poste largo

Probar hasta `T ≈ 1 h` con 3 ventanas independientes pide ≥ 9 h continuas; con holgura para
independencia y para el régimen agitado, **48 h de captura continua**.

A 102 tx/s son ~17.6 M transacciones. Con `MERCADO_DTYPE` a ~40 B/trade, del orden de 700 MB.
Planificar almacenamiento y volcado por bloques **antes** de lanzar: el defecto operativo de la
v2.0 —volcado solo al llenar, con pérdida en parada no limpia— costó varias corridas enteras.

**Lanzar la captura primero.** Los §2 y §3 se ejecutan sobre datos existentes mientras corre.

---

## §5. Experimento 3 — árbitro independiente de la cadena EMD

Si la pregunta es "¿hay exceso de potencia a alguna escala?", conviene un estimador que **sí
tenga nulo**:

- **Multitaper** (Thomson) sobre log-retornos, con nulo de **ruido rojo AR(1)** ajustado a los
  propios datos y bandas de confianza al 95 %.
- Opcional si el multitaper insinúa algo: **wavelet de Morlet** con significancia contra AR(1)
  (esquema Torrence-Compo), que da además localización temporal del pico.

Ejecutar sobre las mismas tres vías del §3. Un pico significativo sobre REAL que no aparezca en
IAAFT ni en BARAJADO es evidencia positiva de escala característica **independiente de la EMD**.

Caso especialmente informativo: **el multitaper encuentra un pico que la EMD no recupera**. Eso
no es (A) ni (B): es **(D)**, la herramienta es la equivocada, y se trata en el §6.

---

## §6. Matriz de decisión — y las alternativas al enfoque de frecuencia

### 6.1 Qué hacer con cada resultado

| resultado | veredicto | acción sobre el diseño |
|---|---|---|
| REAL ≈ BARAJADO, `T/(WΔτ)` constante | **(A)** No hay escala característica | `ω_m` **sale**. Ver §6.2 |
| REAL ≈ IAAFT ≠ BARAJADO | **(A′)** Hay memoria, no ciclo | `ω_m` sale; considerar enfoque de **memoria larga** (Hurst, volatilidad realizada) en lugar de ciclo |
| REAL ≠ ambos, `T` se ancla en Ci | **(B)** La banda estaba mal puesta | Mover `K` y `W` a Ci. **El diseño se conserva entero.** Recalibrar todo lo que dependa de `ω_m` |
| Pico multitaper que la EMD no recupera | **(D)** Herramienta equivocada | Sustituir la cadena de estimación. Ver §6.3 |
| Mezclado o sin resolución | **(C)** No decidible aún | Extender captura. **No decidir.** Es el desenlace más probable y hay que tenerlo previsto |

### 6.2 Si sale (A): qué sobrevive y qué hay que reemplazar

No es el final del proyecto. Sobreviven intactos:

- Toda la v2.0: reloj de ticks, ingesta por lotes, `Q(Δt)`, deduplicación, detección de huecos.
- Toda la capa de riesgo de la v1.3.
- Loeper y el NMPC — necesitan `λ` y `Γ`, no `ω`.
- **`R_n`**, con alta probabilidad: estimar una tendencia es mucho más robusto que extraer un
  ciclo, y el residuo de la EMD no depende de que las IMFs signifiquen algo.

Queda por rediseñar, y hay que decirlo sin adornos:

| pieza | depende de ω_m | reemplazo a explorar |
|---|---|---|
| `A_arm` | sí | velocidad constante permanente (ya es la rama por defecto) |
| `ρ_k`, término `γ_ω·ω_m` | sí | ponderar por volatilidad realizada o por `σ_ω` |
| `c²_vol = k·ω_m·ν = k·f` | sí | volatilidad realizada por tick, que sí es medible |
| `Φ`, `Ψ`, `Ω`, `Ω_crit` | sí | **redefinición completa de la Sec. 1.4** |
| Nodos de fase (Sec. 2.6), `S_ref`, reinicio de `ΣQ` | sí | `S_ref` como media móvil o residuo EMD; reinicio de `ΣQ` por umbral de volumen |

### 6.3 Si sale (D), o si se quiere un enfoque más robusto en cualquier caso

Cuatro alternativas, de menor a mayor cambio arquitectónico:

1. **Espectro → banda → filtro pasabanda → Hilbert.** Multitaper identifica la banda
   significativa; un pasabanda de fase cero extrae la componente; Hilbert da fase y frecuencia
   instantánea. Determinista, con nulo, **sin mezcla de modos**. Es el reemplazo más directo de
   la cadena EMD→HHT y conserva la forma de `ω_m` y `R_n`.

2. **SSA con nulo Monte Carlo.** Adaptativa como la EMD, pero con test de significancia
   establecido (Monte Carlo SSA contra ruido rojo). Devuelve pares oscilatorios explícitos. Es
   el sustituto natural si se quiere conservar el carácter data-adaptive.

3. **Wavelet con significancia Torrence-Compo.** Si el pico existe pero es intermitente,
   localiza *cuándo*. Útil como diagnóstico aunque no sea el estimador de producción.

4. **⭐ Ciclo estocástico dentro del propio EAKF** (componente cíclica de Harvey):

   ```
   ⎡ψ_t ⎤       ⎡ cos λ   sin λ⎤ ⎡ψ_{t-1} ⎤   ⎡κ_t ⎤
   ⎢    ⎥ = ρ · ⎢              ⎥ ⎢        ⎥ + ⎢    ⎥      ρ ∈ (0,1)
   ⎣ψ*_t⎦       ⎣−sin λ   cos λ⎦ ⎣ψ*_{t-1}⎦   ⎣κ*_t⎦
   ```

   Es **`A_arm` con amortiguamiento**, y `λ` se estima **conjuntamente con el estado** en vez de
   inyectarse desde fuera. Ventajas que atacan directamente todo lo hallado en las últimas
   sesiones:

   - **Tiene nulo:** si `ρ → 0` o `σ_κ² → 0` en la estimación, el ciclo no existe. La pregunta
     de este documento pasaría a responderse **dentro** del filtro, de forma continua.
   - **`σ_ω` sale nativa** de la matriz de covarianza de los parámetros: no hace falta la
     descomposición instrumental/genuina del §4.2 de la v2.1.
   - **Desaparece el traspaso externo `ω_m → A_arm`**, que es exactamente donde viven la trampa
     del 2π, la guarda de banda y el problema de rejilla. No hay ventana de análisis, luego no
     hay artefacto de ventana.
   - **`ρ < 1` es el péndulo amortiguado** de tu analogía, ahora como parámetro estimado en vez
     de como intuición: amplitud y fase erran, que es lo que haría un ciclo de mercado real.

   Coste: `A` pasa a depender de parámetros a estimar, así que el filtro deja de ser lineal en
   parámetros. Exige estado aumentado con EKF, o máxima verosimilitud sobre rejilla de `(λ, ρ)`,
   o estimación paramétrica en línea. **Es un cambio de arquitectura de la Sec. 2 entera**, y
   por eso es candidato a v3.0, no un reemplazo directo.

**Nada de esto se implementa en la v2.2.** Se cataloga para que la decisión, cuando llegue, no
se tome bajo presión de tiempo ni con una sola opción sobre la mesa.

---

## §7. Lo que queda bloqueado o diferido

| ítem | estado | motivo |
|---|---|---|
| **A/B de `A`** | **congelado** | Ninguna de las 4 condiciones del §8 de la v2.1 se cumple. Y si sale (A), el A/B deja de tener objeto |
| **v2.1 §7** (correlación `ω_m`–`ΔS`) | **bloqueado** | Mide si `ω_m` covaría con la amplitud. Prematuro |
| **v2.1 §4.4** (`σ_ω/ω` en paralelo a `C`) | **diferido** | Ambos son criterios sobre `ω_m` |
| **v2.1 §5** (log-precio) | **en cola, NO bloqueado** | Independiente de `ω_m`. Puede correr en paralelo a la captura de 48 h |
| **v2.1 §6** (cuantiles empíricos del NIS) | **en cola, NO bloqueado** | Independiente, y las compuertas NIS están mal calibradas pase lo que pase con `ω_m`. Buen candidato para la espera |
| **Fase 2 (ALS)** | sin cambios | Compuerta propia |
| **Fase 3, Testnet, 30 corridas** | sin cambios | — |

---

## §8. Criterios de aceptación

**§2 — Histograma**
- [ ] Distribución de períodos con cotas de banda, por observable.
- [ ] Fracción de rechazos por arriba contra por abajo, como número.
- [ ] Distribución de `T/(W·Δτ)` para los rechazos por arriba.

**§3 — Sustitutos**
- [ ] Tres vías (REAL, IAAFT, BARAJADO) por la cadena **sin modificar**.
- [ ] Comparación de las cuatro magnitudes del §3.2, con distribuciones y no solo medias.
- [ ] **Test de que el nulo es un nulo**: verificar que BARAJADO conserva la distribución
      marginal de incrementos —curtosis, masa en cero, retícula— dentro de tolerancia. Si el
      barajado alteró la marginal, no es el nulo que se pretendía.

**§4 — Multiescala**
- [ ] Captura de **48 h continuas** lanzada **antes** que nada, con volcado por bloques
      verificado contra parada no limpia.
- [ ] Cinco configuraciones C0–C4 × tres observables.
- [ ] `T/(W·Δτ)` por configuración: el estadístico que decide (A) contra (B).
- [ ] `σ_instr²/σ_total²` y validez por configuración y observable.
- [ ] Número de ventanas independientes por configuración **reportado explícitamente**. La
      reserva de la v2.1 —8 ventanas contra 9 rejillas— no debe repetirse en silencio.

**§5 — Árbitro**
- [ ] Multitaper con nulo AR(1) y bandas al 95 %, sobre las tres vías.
- [ ] Comparación explícita: ¿coincide el pico del multitaper con el modo dominante de la EMD?

**Transversal**
- [ ] Suite en verde: **51/51 + los nuevos**.
- [ ] **Ningún cambio al modelo en esta tanda.** Si algo del pipeline hay que tocar para medir,
      va en código de análisis aparte, no en `Micelio.py`. Un experimento que modifica lo que
      mide no mide nada.
- [ ] Veredicto explícito contra la matriz del §6.1, incluida la opción **(C) no decidible**,
      que es un desenlace legítimo y probable.

---

## §9. Fuera de alcance

- **Implementar cualquiera de las alternativas del §6.3.** Se catalogan, no se construyen.
- **Rediseñar la Sec. 1.4** aunque salga (A). Eso sería la v3.0 y necesita su propio documento.
- **Cambiar el observable de producción.** El §4.3 lo evalúa; el cambio, si procede, va después
  del veredicto.
- CMA-ES, IMPC, DeepONet.
- El silencio de `@aggTrade` en fstream. Sin explicación, sin bloquear nada.
- Cadencia del Hilo Rápido en régimen de 518 tx/s — vigilar durante la captura de 48 h, que es
  la primera ocasión de verlo sostenido.
- Reconciliación documental del PDF. La lista sigue creciendo y merece su propia sesión: Secs.
  3.5 y 8.6.1 superadas por ALS; contradicción L1/cuadrática (4.5 contra 6.1, prevalece 6.1);
  dimensiones del NMPC en 7.5.3; Sec. 2.6 como `θ ≡ π/2 (mód π)`; Sec. 4.6 con `n_ticks` por
  transacción; Sec. 2.2 con la banda de resolubilidad como condición de validez; y, pendiente
  del veredicto de este documento, **posiblemente la Sec. 2 entera**.
