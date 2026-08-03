# Canción del Micelio — Contexto del Proyecto

Bot de trading algorítmico. Autor: Samuel Hoyos R. Documento de diseño: `Trading_Bot.pdf` (Julio 2026).
Fuente principal: `Micelio.py` (orquestador consolidado).

**Este archivo recoge un diagnóstico previo hecho sobre `Micelio.py` contrastándolo con el PDF.
Los tests citados se corrieron en un contenedor aislado sin GPU (NumPy/SciPy/CasADi 3.7.2 con
IPOPT y qpOASES). No se han validado contra CUDA ni acados.**

---

## ESTADO ACTUAL (2026-08-03)

Dos tandas de trabajo aplicadas, en este orden:

1. **Correcciones estructurales** — el orden de trabajo de abajo, ejecutado en su totalidad.
   Detalle en "Sesión 2026-08-02".
2. **Fase de calibración** — Sección 0 y Fase 1 de `ORDEN_TRABAJO_CALIBRACION_1.1.md`.
   Detalle en "Sesión 2026-08-02 (b)" al final. **Esta tanda deja obsoletas varias
   afirmaciones de la primera**, señaladas donde corresponde.

Las secciones intermedias se conservan como registro del diagnóstico original — describen
el estado *anterior* del código.

3. **Calibración v1.2** — Secciones A, B, C y la medición de la D de
   `ORDEN_TRABAJO_CALIBRACION_1.2.md`. Detalle en "Sesión 2026-08-03". Integra la cadena
   EMD → Hilbert, que era la precondición que bloqueaba las Fases 2 y 3.

**Las Fases 2 y 3 siguen sin ejecutar**, cada una por su compuerta: la 2 porque Ljung-Box
sigue rechazando blancura (ρ₁ = 0.87 en `y1`), la 3 porque necesita datos de Testnet.
Ninguna de las dos está bloqueada ya por la ausencia del HHT.

### Archivos

- `Micelio.py` — orquestador (3 procesos).
- `constantes_micelio.py` — **única** definición de las constantes de acoplamiento.
- `hht.py` — cadena EMD → Hilbert (Sección 2 del PDF).
- `diagnostico.py` — reporte offline de consistencia del filtro (`python diagnostico.py`).

### Nomenclatura (v1.2 Sec. A.1) — dos símbolos renombrados

El orden de trabajo v1.1 introdujo dos colisiones con el PDF. El PDF llegó primero y es la
fuente, así que se renombró en el orden de trabajo:

| v1.1 | Choca con | Nombre definitivo |
|---|---|---|
| `C_max` (USD) | `C_max` de la Sec. 6.2, que está en **BTC** | **`K_USD`** |
| `ΔS_max` (USD/BTC) | `ΔS_max` de la Sec. 7.4.1, que es un **margen %** | **`ΔS_ref`** |

`C_max` y `ΔS_max` quedan **reservados** para su significado del PDF. La guarda
`verificar_dominio_malla` sigue vigilando que el margen de malla cubra a `ΔS_ref`.

El entorno de esta sesión solo tenía NumPy: **numba, CasADi y pyarrow no están instalados**.
Los núcleos que el PDF asigna a CUDA (Sec. 7.4) y a acados (Sec. 7.5) están implementados
como referencia CPU con el mismo layout de arrays, para que el porte sea mecánico.

---

## Arquitectura (según Sec. 7.6 del PDF)

Tres procesos aislados vía `multiprocessing` + memoria compartida, para evadir el GIL:

- **Proc 1 — Motor de Red** (I/O bound, `asyncio`): WebSocket, Watchdog, Token Bucket,
  cuantización de órdenes (Sec. 8).
- **Proc 2 — Hilo Rápido** (CPU bound): EAKF, disparo de la malla de Loeper a GPU,
  solver NMPC (CasADi/acados), telemetría (Sec. 7).
- **Proc 3 — Hilo Lento** (CPU bound): EMD de ventana deslizante, Transformada de Hilbert,
  proceso Ornstein-Uhlenbeck del Micelio (Sec. 2).

IPC: `shared_memory` con lecturas lock-free vía buffers NumPy. Start method obligatorio `spawn`
(contexto CUDA).

---

## Hallazgo crítico: la EDP de Loeper está mal planteada

**Prioridad máxima. Todo lo demás se construye sobre esta malla.**

El esquema discreto de la Sec. 7.4.3 marcha **hacia adelante** en `t` desde `t0`:

```
U[i,j+1] = U[i,j] - dt * ( 0.5*sigma2*S_i^2*Gamma[i,j] / D[i,j] - costos )
```

Pero de la Ec. 4, `∂U/∂t = -½σ²S²Γ - ...`, que es una parábola **backward**: se resuelve
retrocediendo desde una condición terminal, como todo pricing de opciones. Marchar hacia
adelante la vuelve antidifusiva y por tanto mal planteada.

### Evidencia experimental

Con `λ = 1e-9` (fricción despreciable, singularidad de Loeper imposible), Γ hizo esto:

```
j=10..13:  9.8e-07 → 9.8e-07 → 9.8e-07 → 1.3e-06
j=14..17:  3.1e-06 → 1.3e-05 → 7.4e-05 → 4.3e-04
j=18..21:  2.5e-03 → 1.5e-02 → 9.2e-02 → 8.5e-01   ← D cruza a negativo
```

Multiplicación por ~5 por paso. **No es la Singularidad de Loeper**: `D` se mantuvo en
0.999998 hasta j=20. Es blowup de ecuación de calor invertida amplificando ruido de escala
de malla.

### Por qué CFL no salva esto

Se intentaron tres estrategias, todas fallidas:

1. **CFL literal** (`dt ≤ ds²/2·D_max`): `n_t` explotó a **7057 pasos** por ciclo de control.
   Inviable contra el presupuesto de latencia del Hilo Rápido.
2. **Fijar `n_t`, ensanchar `ds`**: hay una **dependencia circular no documentada** — `ds`
   mayor → dominio más ancho → `S_max` mayor → como el coeficiente difusivo va con `S²`, la
   cota CFL se endurece → `ds` aún mayor. El punto fijo divergió a `ds = 5.9e6`, con el
   dominio cayendo en precios negativos (−119M USD).
3. **Dominio fijo, reducir resolución `n_s`**: "sobrevive" solo porque quedan 5 puntos de malla.

La nota sobre CFL de la Sec. 7.4.1 da una falsa sensación de seguridad: CFL presupone un
problema bien planteado de entrada.

### Confirmación por horizonte corto

Con `tau_pred` = 30 s expresado en años (9.5e-07), 61 puntos de malla: Γ plano en 1.0000e-06
durante los 20 pasos, sin singularidad. El esquema solo "funciona" cuando el horizonte es tan
corto que no alcanza a explotar.

### Corrección

Invertir la marcha: definir condición terminal en `t0 + Δτ_pred` (payoff de la cobertura) e
integrar hacia `t0`, que es donde el NMPC interpola. Alternativa: Crank-Nicolson implícito,
que además libera de CFL — pero en `@cuda.jit` implica solver tridiagonal, más caro que el
explícito.

### Ambigüedad de unidades del horizonte

Sec. 4.5 define `Δτ_pred = T̄/N` con `T̄` en **ticks**. Sec. 7.4.1 lo usa como eje `t` mientras
`r_USD` y `q_BTC` se integran en **años**. No es cosmético: es exactamente la diferencia entre
el caso que explotó y el que quedó estable. La arquitectura de tiempo dual (Sec. 4.5) resuelve
el problema conceptualmente, pero la discretización de 7.4.3 no la refleja — solo tiene un eje
temporal.

---

## Lo que sí está validado

Cadena EMD → Hilbert (Sec. 2), implementada con extensión por reflexión + splines cúbicos:

- Ventanas de 64–256 muestras: **7–12 ms por llamada**. El Hilo Lento corre a 0.5 s → dos
  órdenes de magnitud de margen. **El EMD no es el cuello de botella.**
- Colapso espectral de la Sec. 2.5 (frecuencia ponderada por energía, excluyendo IMF1 y
  residuo): dio 0.4697 Hz sobre señal sintética con modo dominante en 0.4 Hz. Correcto.

Análisis dimensional de los acoplamientos, verificado:

- `λ = ω_m/Ψ` → `(1/Ticks)·(Ticks/BTC) = 1/BTC` ✓
- `c²_vol = k·ω_m·ν` → `(1/Ticks)·(Ticks/Años) = 1/Años` ✓
- **Falta declarar κ** en `R(Ω) = R_base + κΩ²`. Con `[Ω] = BTC/Ticks²`, κ necesita
  `Ticks⁴/BTC²`.

---

## Fallas bloqueantes en `Micelio.py`

El sistema no puede hacer lo que declara el PDF hasta que se corrijan:

1. **`fast_thread_process` nunca abre `shm_mic_name` ni `shm_param_name`** (los recibe como
   argumentos y los ignora). Cascada:
   - Los 11 parámetros de hot-reloading documentados en el encabezado **no se leen en ninguna
     parte del programa**.
   - `Q_k` clavado en `np.eye(3)*0.01` en vez de `ρ_k·diag(q_S,q_v,q_Rn)` con
     `ρ_k = 1+γ_ω|ω_m|+γ_Q|ΣQ|` (Sec. 7.3.3).
   - `R_k` clavado en `diag([0.005,0.02])` en vez de `r_S,base·e^(βτ_d)`.
   - **Toda la modulación endógena del Micelio sobre el filtro — la tesis central del sistema —
     está desconectada.**

2. **Burn-in con umbral estático** `tr_Pk < 1.0`. Sec. 7.1 exige derivada de la traza en banda
   muerta: `|Tr(P_k)−Tr(P_{k-1})|/Δt ≤ ε_burn`, sostenida N ticks continuos.

3. **`Δt` inconsistente**: `A` tiene `0.001` hardcodeado, el loop corre `time.sleep(0.01)`.
   Factor 10. Y Sec. 8 dice que desde Colombia la señal tarda ~300 ms, así que ninguno de los
   dos es defendible. Debe medirse por ciclo (`t_inicio_ciclo` ya se calcula y se descarta).

4. **Watchdog inexistente**: el `except asyncio.TimeoutError` está **fuera** del `while True`.
   Al primer timeout marca dropout y la corrutina termina. Sec. 8.4.2 exige reconexión infinita
   con backoff exponencial.

5. **Token Bucket es un `pass`** (Sec. 8.2, completamente especificado en el PDF).

6. **No hay Ring Buffer**: `ACTUATOR_DTYPE.itemsize` reserva un solo slot. Sec. 7.6.2 exige
   SPSC lock-free con puntero de escritura (el campo `id_accion` ya lleva el comentario).
   Sin spinlocks → las escrituras multi-campo del Hilo Lento son carreras de datos.

---

## Fallas no bloqueantes (el código corre, el resultado es incorrecto)

- **`omega` nunca se escribe**: siempre 0, así que `λ_sim = lam_sim + eta*abs(omega)` deja
  muerto el acoplamiento con Ω de la Sec. 8.1.2.
- **Falta el piso `λ_min`**: la fórmula del PDF es `max(λ_min, λ_ruido + η|Ω|)` precisamente
  "para evitar fricciones negativas computacionales". El OU con σ=0.2 sí cruza a negativo.
- **`dt_sim = 0.1` contra `time.sleep(0.5)`**: el proceso OU revierte a la media 5× más lento
  en tiempo de pared de lo calibrado.
- **`mu_OU = 0.5` hardcodeado**: con `[λ] = 1/BTC` (Sec. 4.4.1), λ≈0.5/BTC es enorme —
  `λS²Γ` dispararía la singularidad de inmediato a precios de BTC. Debe venir de métricas
  reales de Mainnet.
- **`apply_filters` nunca se llama**, y tiene bug de orden: valida `minQty`/`minNotional`
  **antes** del floor, así que la cuantización puede devolver una orden que ya no cumple
  `minNotional`. Revalidar después de discretizar.
- **Telemetría descarta datos**: `telem_buffer` es lista de Python y al llenarse se vacía sin
  escribir. Sec. 8.6.2 pide buffer NumPy preasignado con sobrescritura O(1) y volcado a
  parquet. Además solo guarda `Tr(P)`; Sec. 8.6.1 exige también `x_k` y la innovación `ỹ_k`
  — **sin `ỹ_k` no se puede hacer el Covariance Matching offline, que es el propósito
  declarado del logging.**

### Menores

- `np.linalg.inv(S)` sin guarda → usar `np.linalg.solve`.
- `asyncio.get_event_loop()` deprecado en Python 3.12 → `asyncio.run()`.
- Sin pinning a núcleos (`os.sched_setaffinity`) pese a Sec. 7.6.1.
- `main()` sin `try...finally`: depende solo del signal handler para el `unlink` de
  `/dev/shm`. Cualquier excepción no-señal filtra memoria compartida.

---

## Contradicciones internas del PDF (corregir el documento, no solo el código)

1. **L1 vs cuadrática**: Sec. 4.5 dice "los costos de exchange aplican como penalizaciones
   continuas de Norma L1". Sec. 6.1 dice "se descarta el uso directo de normas L1" y justifica
   la cuadrática por no diferenciabilidad en SQP. **6.1 tiene la razón**; borrar la frase de 4.5.

2. **Dimensiones del NMPC**: Sec. 7.5.1 define `x^c_k ∈ R^(1×1) = [I_k]`. Sec. 7.5.2 define
   `Q(Ω) ∈ R^(2×2)`. Sec. 7.5.3 escribe `(x^c_k)^T Q(Ω) x^c_k` — un escalar no se multiplica
   por una 2×2. Y Sec. 6.1 formula el costo sobre `e_k`, no sobre `x^c_k`.
   **Reconciliación propuesta**: el vector penalizado es `[e_k; I_k] ∈ R²`, con
   `Q = diag(q_Δ, q_inv(Ω))` y `e_k = I_k − Δ_k` (relación algebraica ya definida en 7.5.1).

3. **Condición inicial/terminal de U ausente**: Sec. 7.4 nunca especifica `U(S,·)` en el borde
   temporal ni la forma explícita de `costos`. En los tests hubo que **asumir**
   `U = ½γ₀(S−S_k)²` y `costos = q_BTC·U − (r_USD−q_BTC)·S·∂U/∂S` (derivado de Ec. 4).
   Tal como está documentado, **el solver de GPU es irreproducible**.

---

## Orden de trabajo recomendado — TODO EJECUTADO

1. ~~Inversión temporal de la EDP de Loeper~~ ✔ `resolver_malla_loeper`
2. ~~Conectar `shm_param` y `shm_mic` al Hilo Rápido; restaurar modulación de Q y R~~ ✔
3. ~~Watchdog con backoff + Token Bucket~~ ✔ clase `TokenBucket`
4. ~~Burn-in por derivada de la traza~~ ✔ (ver advertencia sobre ε_burn más abajo)
5. ~~Ring Buffer SPSC + spinlocks~~ ✔ seqlock + secuencia monótona global
6. ~~Telemetría con `ỹ_k` y volcado real~~ ✔ (a `.npy`; parquet requiere pyarrow)

---

## Sesión 2026-08-02 — Correcciones aplicadas y hallazgos nuevos

### La inversión de Loeper, verificada por contraste

Se corrió el esquema **forward literal de la Sec. 7.4.3** contra el **backward corregido**,
sobre la misma malla, el mismo `dt` y con `λ = 1e-9` (singularidad imposible):

| | Γ inicial | Γ final | D |
|---|---|---|---|
| Forward (7.4.3) | 2e-8 | **1.07e+3** | cruza a negativo en **j=21** |
| Backward (corregido) | 2e-8 | **2e-8** | 0.99999996 |

Reproduce el hallazgo original, **incluido el cruce exacto en j=21**. Con el problema ya
bien planteado, CFL recupera su sentido: 150 pasos para 1 h de horizonte, **3 pasos** para
los 30 s operativos. Latencia medida Loeper+NMPC: **2.2 ms/ciclo** contra 10 ms de
presupuesto del Hilo Rápido.

Sobre la **ambigüedad de unidades del horizonte**: se resolvió documentando que un único eje
en años es suficiente, porque `ν` entra exclusivamente a través de `c²_vol = k·ω_m·ν`
([1/Ticks]·[Ticks/Años] = [1/Años]). No hacen falta dos ejes que sincronizar.

### Huecos NUEVOS del PDF (corregir el documento)

Todos son de la misma familia que el κ faltante: **fórmulas correctas cuyas constantes no
tienen magnitud ni unidades declaradas**. Cada una rompió el sistema en ejecución.

1. **`τ_max` de la Sec. 6.5 nunca se implementó.** El PDF dice "si τ_d supera un umbral de
   tolerancia máximo τ_max, el paquete se descarta", pero no había tal umbral en el código:
   en el primer ciclo `timestamp` vale 0, `τ_d ≈ 1.75e9 s` y `e^(β·τ_d)` desborda con
   `OverflowError`. Añadido como `P_TAU_MAX`.
   **Sutileza:** `τ_d` puede salir *levemente negativo* (el Hilo Rápido muestrea `t_ahora` al
   inicio del ciclo y el Motor de Red publica un timestamp más nuevo mientras trabaja). Eso
   es granularidad de reloj, no un paquete inválido; hay que saturarlo a 0. Rechazarlo anula
   K, impide que P se contraiga y **rompe la racha del burn-in indefinidamente**.

2. **Sec. 7.3.3 no declara unidades de `γ_ω` ni `γ_Q`.** Con `[ω_m] = 1/Ticks` y `[ΣQ] = USD`
   (Sec. 1.1), para que ρ_k sea adimensional hace falta `[γ_ω] = Ticks` y `[γ_Q] = 1/USD`.
   Con `γ_Q = 0.2` y ΣQ ~ 1e5 USD salía **ρ_k ≈ 2×10⁴** y la traza nunca entraba en banda
   muerta.

3. **Sec. 8.1.2 no declara la magnitud de `η`.** En `λ_sim = máx(λ_min, λ_ruido + η|Ω|)` el
   término `η|Ω|` debe ser *comparable* a `λ_ruido ~ μ_OU`, no dominarlo. Con `η = 0.15` y
   `|Ω| ~ O(1)` salía `λ_sim ≈ 0.15/BTC` → `λS²Γ ≈ 6` → **singularidad permanente**.

4. **Sec. 8.1.1 no acota `σ_OU`.** La desviación estacionaria de un OU es `σ/√(2θ)`; para que
   λ fluctúe *alrededor* de μ_OU hay que exigir `σ/√(2θ) ≪ μ_OU`. Con `σ = 0.05, θ = 0.5`
   salía `σ_est = 0.05`, **cien veces μ_OU**. (Esto generaliza la nota previa de que "el OU
   con σ=0.2 sí cruza a negativo": el problema no es solo el signo, es la escala.)

5. **La Sec. 7.4.4 solo devuelve `M_Γ`, pero la 7.5.1 exige `Δ_k` como TVP.** El kernel de
   GPU debe devolver **ambas** superficies, `M_Γ` y `M_Δ`. Tal como está documentado, el NMPC
   no tiene de dónde sacar su objetivo de cobertura.

6. **`Δ` y `I_k` no son dimensionalmente comparables.** La 7.5.1 escribe `e_k = I_k − Δ_k`
   con `I_k` en BTC, pero de la 4.4.1 `[Δ] = ∂U/∂S = BTC²/USD`. La cantidad de cobertura
   denominada en BTC es `S·∂U/∂S`. Se usa esa conversión.

7. **Condición terminal de U** (contradicción #3 previa): se ancló el payoff
   `U = ½γ₀(S − S_ref)²` al **nodo de fase `S_ref` de la Sec. 2.6**, no al precio actual.
   Anclarlo en `S_k` daría `Δ(S_k) ≡ 0` por simetría y el objetivo del NMPC sería
   trivialmente nulo. Anclarlo en `S_ref` hace que la cobertura siga el desplazamiento
   estructural `ΔS`, que es la variable que el modelo declara significativa.

8. **`f_Loeper` de la Sec. 7.5.2 se nombra pero no se define.** Se asumió
   `f_Loeper = R_base·(1/D_k − 1)`: homogénea con `R_base`, nula sin impacto (Γ→0) y
   divergente en la singularidad — que es el comportamiento que la propia 7.5.2 describe.

### ⚠ ε_burn NO es una constante estática — RESUELTO en la sesión (b) por el criterio NIS

La Sec. 7.1 da el criterio pero no el valor. Medido contra el generador sintético:
`dTr/dt` vive en ~0.01–0.03 en régimen estacionario y el burn-in cierra en ~215–240 ciclos
(≈2.4 s) con `ε_burn = 1.0`.

**Pero tratarlo como constante es incorrecto al entrar a mercado real.** La derivada de la
traza depende de `Δt` (en producción es la latencia real medida, ~300 ms desde Colombia, no
la cadencia nominal), de `ρ_k` —que el Micelio modula tick a tick vía `ω_m` y `ΣQ`— y del
régimen de volatilidad vigente. Un umbral fijo calibrado en régimen tranquilo **bloqueará**
el burn-in en uno agitado; calibrado en régimen agitado **dejará pasar** un filtro que aún no
convergió. En Mainnet debe hacerse **adaptativo**: normalizarlo contra la escala de `Tr(P)`
(criterio relativo) o contra un percentil móvil de la propia derivada.

### Otros defectos encontrados en ejecución

- **El Ring Buffer SPSC se atascaba en silencio.** El Hilo Rápido publica a ~100 Hz pero el
  presupuesto de pesos de Binance (Sec. 8.2) sostiene del orden de **1 orden/s**: el
  productor lapea al consumidor y sobrescribe slots. Sin detección de sobrepaso el consumidor
  espera para siempre una secuencia que ya no existe y **el lazo de control se rompe sin
  emitir ningún error** (el inventario simplemente deja de actualizarse). Al resincronizar se
  salta al **último** publicado, no al más antiguo: Sec. 8.2.3 purga órdenes obsoletas y por
  horizonte recedente (Sec. 6.4) solo el `U_0` más fresco es válido.
- **El mock de precio era ruido blanco ±1 USD.** El detector de nodos de fase disparaba en
  casi todos los ticks, ΣQ se reiniciaba sin parar y ρ_k saltaba. Sustituido por una señal
  con ciclo estructural real (150 USD, 40 s), que es lo que el modelo de la Sec. 1 supone.

### Verificación end-to-end

Los tres procesos corren estables; burn-in cierra a los ~238 ciclos; se emiten ~460 órdenes
en 20 s; el **inventario sigue la cobertura** (sube a 0.123 BTC al alejarse el precio del
nodo, baja a −0.008 al revertir); el freno de singularidad entra y se despeja correctamente;
telemetría con los 7 campos incluida `ỹ_k`; sin fugas en memoria compartida al cerrar.

### Pendiente

- `ω_m` y `R_n` siguen siendo mocks de Testnet en `Micelio.py`. La cadena EMD → Hilbert está
  validada (7–12 ms) pero **vive fuera de este archivo**; falta integrarla.
- Todos los parámetros marcados `[CALIBRAR]` requieren métricas reales de Mainnet.
- Sin validar contra CUDA ni acados.

---

## Sesión 2026-08-02 (b) — Fase de calibración

Ejecuta la **Sección 0** y la **Fase 1** de `ORDEN_TRABAJO_CALIBRACION_1.1.md`.
Las Fases 2 y 3 quedan sin hacer, y por razones del propio documento (ver abajo).

### Sección 0 — las constantes de acoplamiento son fórmulas, no valores

`constantes_micelio.py` es ahora la **única** definición de γ_0, γ_ω, γ_Q, κ y μ. Ninguna
aparece como literal en `Micelio.py`: el bloque de hot-reloading almacena los **límites
estructurales** (I_max, C_max, ω_m,max, ΔS_max, Ω_crit, ΣQ_max, c, c') y las derivadas se
evalúan llamando al módulo.

- `γ_0 = I_max/(S·ΔS_max)` **depende de S**, así que se reevalúa cada ciclo con el precio
  filtrado. Verificado: la cobertura implícita en ΔS_max da exactamente I_max = 0.5 BTC.
- `κ` y `μ` se evalúan dentro de `resolver_nmpc` a partir de Ω_crit. Comprobado que la
  cascada funciona: con Ω_crit = 1, 2 o 4, `R(Ω_crit)` da **101× R_base en los tres casos**.
  Esto es lo que hace segura la búsqueda de orden cero de la Fase 3.
- Guardas de la Sec. 0.5 al arranque del Hilo Rápido, ambas probadas contra su fallo.
- Constantes de ruido (Sec. 0.6) **sin tocar**. La guarda de σ_OU advierte sin corregir.

**Dos colisiones de nombres entre documentos**, ambas documentadas en el módulo:
`C_max` (BTC en la Sec. 6.2 del PDF, USD en la Sec. 0.2 del orden de trabajo) y `ΔS_max`
(margen relativo de malla en 7.4.1, desplazamiento absoluto en USD/BTC en 0.2). La segunda
es peligrosa y tiene guarda propia (`verificar_dominio_malla`).

### Fase 1 — el NIS cambió el veredicto del sistema

`ε_k = ỹᵀS⁻¹ỹ` se calcula en línea y tiene los dos consumos que pide el documento:
telemetría (el escalar, no la matriz) y la **ventana adaptativa W_k de la Sec. 2.2.1**, que
nunca se había implementado.

**El burn-in pasa a criterio NIS** (`DIVERGE DEL PDF (Sec. 7.1)` en el código). Esto resuelve
la advertencia sobre ε_burn de la sesión anterior: el NIS es adimensional y auto-normalizado,
su banda χ²_m no depende de Δt, ni de ρ_k, ni del régimen.

Y lo primero que hizo fue delatar al criterio viejo. Con 115 s de telemetría real:

```
NIS medio = 13.93   (teórico 2)        -> SOBRECONFIADO
legacy dTr/dt: mediana 0.0053          -> racha de 1391 ticks "convergido"
Ljung-Box y0 / y1 / multivariante      -> RECHAZA BLANCURA en los tres
```

El criterio de la traza daba el filtro por sano mientras subestimaba su incertidumbre ~7×.

**Fase 2 (ALS) queda bloqueada por su propia compuerta, y es correcto.** Los mocks son
deterministas: `R_n` es un coseno puro muestreado a 0.5 s y el precio una sinusoide de 40 s
que un modelo de velocidad constante no puede seguir. Hay *model mismatch* genuino, así que
correr ALS absorbería el error de modelo dentro de Q. Concuerda con lo que el propio orden de
trabajo advierte: integrar la cadena EMD → Hilbert es **precondición** para que el NIS
signifique algo, porque `R_n` entra directo en el vector de medición.

⚠ **Al leer el reporte de `diagnostico.py`, mirar siempre el recorte de transitorio.** Sin
descartar el arranque en frío, el NIS daba 3333 y la curtosis 3568 por la innovación del
primer ciclo (`x2 = 0` contra `R_n ≈ 45000`, perfectamente legítima por la Sec. 7.3.4). Con
recorte automático el veredicto sobrevive, así que es real — pero los tres tests de la Fase 1
dan falsos positivos sin él.

### Defectos corregidos de paso

- El alias del módulo de constantes colisionaba con `K`, la ganancia de Kalman: Python trata
  `K` como local en toda la función y las guardas del arranque reventaban con
  `UnboundLocalError`. Renombrado a `CTE`.
- **Un `print` fallido podía matar el Hilo Rápido en silencio.** Bajo `spawn` los hijos
  heredan el handle de stdout del padre; si el consumidor de ese pipe se cierra, el siguiente
  print lanza `BrokenPipeError` y el supervisor solo ve morir un hijo sin causa. Todo el
  diagnóstico pasa ahora por `log()`, que nunca propaga. El supervisor además ya dice **qué**
  proceso murió y con qué exitcode.
- Los mensajes que se imprimen van en ASCII: la consola de Windows es cp1252 y no puede
  codificar letras griegas.

### Pendiente tras esta fase

- **Integrar la cadena EMD → Hilbert.** Es la precondición de todo lo demás. Al hacerlo,
  aplicar la conversión de la Sec. 0.4: el HHT entrega f en **Hz** y el modelo exige
  **1/Ticks** (`CTE.omega_m_desde_hz`); ojo con que ν se almacena en Ticks/**Años**.
- Fase 2 (ALS) — bloqueada hasta que Ljung-Box pase.
- Fase 3 (Ω_crit por búsqueda de orden cero) — necesita datos de Testnet.
- Marcar en el PDF las Secs. 3.5 y 8.6.1 como superadas por ALS (la telemetría de `ỹ_k` sigue
  siendo correcta; cambia el método que la consume, no el dato).

---

## Sesión 2026-08-03 — Integración de la cadena EMD → Hilbert (v1.2 A, B, C)

### `hht.py` — la Sección 2 del PDF, por fin dentro del sistema

`ω_m` y `R_n` dejan de ser mocks. El Hilo Lento mantiene una ventana deslizante de precios,
ejecuta EMD → Hilbert y consume el resultado completo: frecuencia consolidada por colapso
espectral (Sec. 2.5), residuo real como `R_n`, y nodos de fase (Sec. 2.6) que actualizan
`S_ref` y reinician ΣQ.

Cuatro cosas que hubo que descubrir midiendo, porque el PDF no las cubre:

1. **La Transformada de Hilbert va sobre el vector EXTENDIDO, no sobre las IMFs truncadas.**
   El paso 2 de la Sec. 2.3 lo dice, pero es fácil truncar primero. Hacerlo mal reintroduce un
   borde duro exactamente en t0 —el único instante que el control usa— y daba errores de
   frecuencia del **200-600 %**. Corregido: 2.1 % sobre el ciclo de 40 s.

2. **ω_m no se puede evaluar en una sola muestra.** La Sec. 2.5 la define como instantánea,
   pero `f` sale de derivar numéricamente la fase desenrollada y una muestra suelta es puro
   ruido: cambiar la semilla movía el error del 2 % al 45 %. Se toma la **mediana** sobre las
   últimas 12 muestras (mediana y no media: los saltos del desenrollado son atípicos aislados).
   Encima, EMA temporal en el Hilo Lento, porque ω_m alimenta ρ_k y c²_vol.

3. **Dimensionado de la ventana: manda muestras/ciclo, NO ciclos/ventana.** Barrido de
   dt ∈ {0.25, 0.5, 1.0} s × W ∈ {128…512}, 12 semillas por punto:

   | muestras/ciclo | error mediano |
   |---|---|
   | 40 (dt=1.0 s) | 30–58 % |
   | **80 (dt=0.5 s)** | **5–17 %** |
   | 160 (dt=0.25 s) | 9–48 % (pierde ciclos) |

   Un intento previo de muestrear a 2 s "para abarcar más ciclos" degradó el error al 61 %.
   Config final: muestreo 0.5 s, `W_min=192`, `W_max=384`. Latencia 6–10 ms contra 500 ms.
   **Warm-up de ~96 s** antes del primer tamizado: es real, la EMD necesita esa historia.

4. **El nodo de fase se detecta por el signo de la IMF dominante en t0 ENTRE iteraciones**,
   no comparando `imf[-2]` contra `imf[-1]` dentro de una misma ventana. El borde se recalcula
   en cada llamada y la comparación intra-ventana no disparaba casi nunca. Con el cambio:
   **13 nodos detectados contra 14 reales**, desviación mediana 1.0 s. Más período refractario
   de 0.8 medios ciclos para matar el chatter.

### ⚠ Otra inconsistencia del PDF: la Sec. 2.6 se contradice

La prosa define el nodo como "cruza el eje horizontal (**amplitud cero**)" y lo refuerza con
"El retorno del precio exactamente al nodo fuerza a ΔS → 0". Pero lo formaliza como
`θ_unwrap ≡ 0 (mód π)`, y con la convención estándar `Re(Z) = IMF = A·cos θ`:

```
θ ≡ 0   (mód π)  ⟺  IMF = ±A   -> EXTREMO de la oscilación
θ ≡ π/2 (mód π)  ⟺  IMF = 0    -> CRUCE POR CERO
```

La fórmula detecta los picos, justo lo contrario de lo que pide la prosa. Se implementó el
cruce por cero (lo que exige la física del modelo). **Corregir en el PDF**: la condición debe
ser `θ ≡ π/2 (mód π)`, o reescribirse directamente como cruce por cero de la IMF dominante.

### Sección C — `q_S` relativa

`q_S = (σ_rel·S_k)²`, mismas unidades (USD/BTC)² que ya tenía, así que el álgebra del filtro
no cambia; lo que cambia es que deja de depender del nivel de precio (con `q_S` fijo, BTC de
45k a 90k lo desescala por 4×). `σ_rel = 2.22e-6` reproduce exactamente el `q_S = 0.01`
anterior a S = 45 000, así que no altera el comportamiento actual. **No es una
adimensionalización del filtro** — esa se evaluó y se descartó, ver la nota en `q_S_relativa`.

### A.2 — Ljung-Box ahora se lee por magnitud

`diagnostico.py` imprime ρ₁, ρ₂, ρ₃ junto al p-valor, el `n` efectivo tras el recorte, y el
ρ₁ mínimo detectable con ese `n`. La función `_rho_minimo_detectable` reproduce **exactamente**
la tabla de la v1.2 (n=500 → 0.250, n=2 000 → 0.125, n=11 500 → 0.052, n=50 000 → 0.025).
El veredicto distingue "RECHAZA BLANCURA" (ρ₁ ≥ 0.20, model mismatch) de "rechazo NO MATERIAL"
(ρ₁ < 0.20, con n grande Ljung-Box rechaza por correlaciones irrelevantes).

⚠ **PROHIBIDO inflar Q o R con una constante ad-hoc.** Está medido en la v1.2: Q×7 lleva el
NIS de 19.35 a 3.22 pero solo baja ρ₁ de +0.79 a +0.57, y Q×50 lo deja en +0.44. Un error
determinista no lo describe ninguna matriz de covarianza — la inflación no corrige, oculta.

### Sección D — la medición, y una expectativa del documento que NO se cumplió

La v1.2 predice: "buena parte del exceso de NIS debería evaporarse sola, porque `R_n` real es
lentamente variable". **No ocurrió.** Medido sobre dos bloques de 10 000 muestras:

| | NIS medio | NIS mediana | var(y0) vs `r_S,base`=0.5 | var(y1) vs `r_EMD`=1.0 |
|---|---|---|---|---|
| Bloque 0 (pre-HHT, `R_n` de respaldo) | 848 | 4.11 | 446 | 0.40 |
| Bloque 1 (**HHT activo**) | 147 | 4.69 | 4.03 (**8×**) | 161 (**161×**) |

Lecturas, en orden de importancia:

1. **La mediana del NIS es ~4.7 contra un teórico de 2** — solo 2.3× de exceso. La *media* de
   147 está inflada por cola pesada (|y1| máximo 280 USD). Mirar la mediana.
2. **`R` está gruesamente mal especificada, y es medible directamente.** `var(y0) = 4.03`
   contra `r_S,base = 0.5` es exactamente la varianza del ruido que el mock inyecta
   (`gauss(0, 2)` → 4.0): el filtro cree tener 8× más precisión de la que tiene. Y `r_EMD`
   está 161× por debajo. Eso es *justo* lo que la Fase 2 (ALS) existe para estimar.
3. **`R_n` real NO es lentamente variable.** Sale de una EMD de ventana deslizante recalculada
   a 2 Hz; al desplazarse la ventana el residuo salta, y salta más cuando dispara un nodo de
   fase. La premisa de la v1.2 sobre este punto era optimista.

### ⚠ Cómo leer ρ₁ aquí: no todo rechazo es evidencia sobre `A`

`diagnostico.py` ahora reporta, junto a ρ₁, con qué frecuencia cambia cada componente y cómo
decae la autocorrelación. Con eso, los dos rechazos significan cosas distintas:

- **`y1` (ρ₁ = 0.87)**: `R_n` proviene de ventanas EMD **solapadas**, así que valores
  consecutivos están correlacionados *por construcción*. En el bloque 0 la componente llegaba
  a repetirse ~8 ciclos idénticos (retención de orden cero: el Hilo Lento publica a 2 Hz, el
  Rápido corre a 89 Hz). Su ρ₁ **no es evidencia sobre la matriz `A`**.
- **`y0` (ρ₁ = 0.73, decae a 0.12 en el rezago 10)**: error de modelo suave y de corto
  alcance. Ese sí apunta a lo que describe la Sección E — `A` asume velocidad constante contra
  una oscilación.

Un mismatch estructural del tipo de la Sección E dejaría correlación a rezagos del orden del
ciclo estructural (miles de muestras), no solo en ρ₁. Aquí no aparece.

### En Windows `unlink()` no destruye la memoria compartida — el bot no arrancaba dos veces

Defecto real de arranque, encontrado en ejecución. En POSIX `shm_unlink` retira el nombre de
inmediato; en Windows el bloque es un objeto del kernel con conteo de referencias y **vive
mientras algún proceso lo tenga mapeado**. `SharedMemory.unlink()` es, de hecho, un no-op.

Consecuencia: tras una caída sucia (o con un proceso hijo rezagado que el padre no alcanzó a
matar en cascada), el bloque sobrevive y el arranque siguiente moría con
`FileExistsError: [WinError 183]` al intentar destruir-y-recrear.

`allocate_shared_memory` ahora **se adjunta al bloque huérfano y lo recicla** si su tamaño
alcanza; solo si es más chico intenta liberarlo y recrearlo. Probado contra el fallo real
(bloque huérfano de 57 bytes retenido por un PID vivo).

### Higiene de secretos

`.gitignore` excluye ahora `*Credenciales*`, `*.key`, `*.pem`, `.env*` y `secrets.*`.
El repo `origin` es **público**: un secreto commiteado sobrevive en el historial de git
aunque después se borre el archivo, así que la única defensa barata es que nunca entre.
Las credenciales de la cuenta demo de Binance viven **fuera** del árbol de trabajo.

### Pendiente tras esta fase

- **Fase 2 sigue BLOQUEADA** por la compuerta (ρ₁ = 0.87 ≥ 0.20), y así se deja: el documento
  manda no correr ALS con mismatch. Pero conviene decidir en la v1.3 si el ρ₁ de `y1` debe
  contar para la compuerta, dado que es un artefacto del solapamiento de ventanas.
- **`r_S,base` y `r_EMD` están medidos y mal por 8× y 161×.** No se tocaron: la Fase 2 es su
  dueña. Es el candidato número uno a explicar el NIS residual.
- Fase 3 (Ω_crit) — necesita Testnet.
- Sección E (modelo de oscilador armónico en vez de velocidad constante) queda **aparcada por
  decisión explícita** hasta que el bot opere en cuenta demo. `A` se queda como está.
- `ω_m,max` en `constantes_micelio` sigue siendo un TODO(HHT): ahora que la cadena está
  integrada, debe fijarse con el máximo empírico observado.

## Convenciones

- Comentarios y nombres de variables en español, consistente con el código y el PDF existentes.
- Referenciar la sección del PDF en los comentarios al implementar una fórmula.
- Cualquier suposición que rellene un hueco del PDF debe marcarse explícitamente con
  `# NOTA DE INTERPRETACION:` y describir qué se asumió.
- Cuando una decisión **contradiga** al PDF (como el NIS sobre ε_burn, o ALS sobre Covariance
  Matching), no aplicarla en silencio: marcar `# DIVERGE DEL PDF (Sec X.Y):` con la
  justificación, para poder reconciliar el documento después.
- Las constantes de acoplamiento se **derivan** en `constantes_micelio.py`, nunca se escriben
  como literales. Un número mágico reabre el agujero que esas fórmulas cierran.
- Todo texto que se **imprime** va en ASCII (la consola es cp1252); los comentarios y
  docstrings sí llevan acentos y símbolos.
