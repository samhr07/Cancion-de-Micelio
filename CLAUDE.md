# Canción del Micelio — Contexto del Proyecto

Bot de trading algorítmico. Autor: Samuel Hoyos R. Documento de diseño: `Trading_Bot.pdf` (Julio 2026).
Fuente principal: `Micelio.py` (orquestador consolidado).

**Este archivo recoge un diagnóstico previo hecho sobre `Micelio.py` contrastándolo con el PDF.
Los tests citados se corrieron en un contenedor aislado sin GPU (NumPy/SciPy/CasADi 3.7.2 con
IPOPT y qpOASES). No se han validado contra CUDA ni acados.**

---

## ESTADO ACTUAL (2026-08-05)

Cinco tandas de trabajo aplicadas, en este orden:

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

4. **Riesgo de cuenta y modelo oscilatorio v1.3** — `ORDEN_TRABAJO_RIESGO_1_3.md`, completo:
   precondición (Modo LECTURA) y Secciones A, B, C, D, E y F. Detalle en
   "Sesión 2026-08-04".
5. **Relojes y acoplamiento con el mercado v2.0** — `ORDEN_TRABAJO_RELOJES_2_0.md`, completo:
   §2 a §8. Detalle en "Sesión 2026-08-05". Reorganización estructural: el filtro y el
   EMD pasan a **reloj de transacciones** (Δn = 1) y el sistema deja de descartar el
   **96 %** de los datos que ya recibía. **42/42** criterios (`python tests_v13.py`).

**Las Fases 2 y 3 siguen sin ejecutar.** La 3 necesita Testnet. La 2 sigue tras su
compuerta, pero la v1.3 cambió el terreno: el ρ₁ = 0.87 de `y1` resultó ser artefacto de
corregir el filtro 90 veces con la misma medición, y eso está corregido en origen (ver
"El filtro corregía 90 veces por medición").

### Archivos

- `Micelio.py` — orquestador (3 procesos).
- `constantes_micelio.py` — **única** definición de las constantes de acoplamiento y de
  los límites de cuenta (Sec. 1.bis).
- `hht.py` — cadena EMD → Hilbert (Sección 2 del PDF) + concentración espectral `C`.
- `mercado.py` — **v1.3**: tri-estado `MODO`, lectura de `exchangeInfo`, feed público real.
- `riesgo.py` — **v1.3**: capa de riesgo de cuenta, 7 guardas, ruta de cierre.
- `episodios.py` — **v1.3**: máquina de episodios y adaptador de faucet.
- `dinamica.py` — **v1.3**: matrices `A`, conmutador de rama, EAKF sombra.
- `diagnostico.py` — reporte offline (`python diagnostico.py [--episodio=N]`).
- `tests_v13.py` — criterios de aceptación de la Sec. F de la v1.3 **y del §8 de la v2.0**
  (`python tests_v13.py`). El nombre se conserva para no romper referencias.
- `test_contaminacion_emd.py` — **v2.0 §5.2** sobre mercado real
  (`python test_contaminacion_emd.py --capturar=440`).

### Cómo se arranca

```
python Micelio.py                     # arranca en LECTURA (feed real, sin ejecucion)
MICELIO_MODO=TESTNET python Micelio.py

python tests_v13.py                   # 33 criterios de aceptacion de la Sec. F
python tests_v13.py --sin-red         # omite los que consultan exchangeInfo

python diagnostico.py                 # reporte de consistencia + A/B de la Sec. E
python diagnostico.py --episodio=3    # un solo episodio (Sec. C.5)
```

`MODO` es tri-estado y su valor por omisión es **LECTURA**. Elevarlo a TESTNET o MAINNET
tiene que ser un acto explícito del operador: el modo que puede tocar la cuenta nunca es
el que sale por descuido.

⚠ **No arranques dos bots a la vez.** `verificar_instancia_unica` lo impide desde la v1.3,
pero conviene saber por qué existe: antes se adjuntaban en silencio a la misma memoria
compartida y los datos de ambos quedaban inservibles sin ningún error. Si acabas de matar
uno, espera 5 s (tolerancia del latido) antes de arrancar el siguiente.

⚠ **El HHT tarde ~192 s en dar su primera estimación** (384 muestras a 0.5 s). Hasta
entonces `C = 0`, `ω_ang = 0` y la rama de `A` es velocidad constante. Una corrida más corta
que eso no dice nada sobre las Secciones D y E.

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
  → **RESUELTO en la v1.3, y no como se esperaba.** No era el solapamiento de ventanas: era
  que el filtro corregía ~45 veces con la misma medición. La Sec. E.3 además excluye `y1` de
  la compuerta por decisión del documento. Ver "El filtro corregía 90 veces por medición".
- **`r_S,base` y `r_EMD` están medidos y mal por 8× y 161×.** No se tocaron: la Fase 2 es su
  dueña. Es el candidato número uno a explicar el NIS residual. → **Sigue vigente.**
- Fase 3 (Ω_crit) — necesita Testnet. → **Sigue vigente.**
- Sección E (modelo de oscilador armónico en vez de velocidad constante) queda **aparcada por
  decisión explícita** hasta que el bot opere en cuenta demo. `A` se queda como está.
  → **DESAPARCADA por la v1.3**, que la promueve a su Sección D. `A_arm` está implementada,
  con conmutación por concentración espectral. El veredicto sobre adoptarla sigue pendiente
  del A/B.
- `ω_m,max` en `constantes_micelio` sigue siendo un TODO(HHT): ahora que la cadena está
  integrada, debe fijarse con el máximo empírico observado. → **Sigue vigente**, pero ya hay
  material: sobre mercado real se observaron períodos de 33–188 s.

## Sesión 2026-08-04 — Riesgo de cuenta y modelo oscilatorio (v1.3 completa)

Ejecuta `ORDEN_TRABAJO_RIESGO_1_3.md` entero: precondición y Secciones A–F. Es la primera
tanda en la que el bot corre contra **datos reales de Mainnet**.

### RESULTADOS DE UN VISTAZO

Todo lo de abajo está **medido en ejecución**, no razonado. Cada fila tiene su subsección
más abajo con el detalle y el porqué.

**Estado de las secciones del orden de trabajo**

| Sección | Estado | Nota |
|---|---|---|
| Precondición — Modo LECTURA | ✔ operativo | feed público real de Mainnet, sin credenciales |
| A — filtros y resolución | ✔ completa | dos cifras del documento corregidas |
| B — riesgo de cuenta | ✔ completa | las 7 guardas + disparo forzado probado |
| C — episodios y faucet | ✔ completa | 6 compuertas, `DETENIDO` terminal |
| D — matriz `A` oscilatoria | ✔ completa | ambas trampas verificadas por test |
| E — protocolo A/B | ⚠ infraestructura sí, **veredicto no** | el feed degradado lo invalida; y la v2.0 §5.2 lo anula por segunda causa |
| F — criterios de aceptación | ✔ **33/33** | `python tests_v13.py` |

**Constantes y guardas al arranque, con S real de 63 920 USD**

| Magnitud | Valor medido | Criterio |
|---|---|---|
| `nocional_max_orden` | 3 000 USD = **46.9 lotes** | ≥ 40 exigidos (Sec. A.3) |
| `nocional_max_posicion` | 32 599 USD | derivado de `I_max`, no declarado |
| colchón de liquidación a 5× | 6 259 USD contra cap de 3 000 | **2.1×** de holgura |
| ídem a 20× | 1 366 USD | **rechazado**, como debe |
| `equity_min_episodio` | 9 520 USD | por debajo, no se arranca |
| `λS²γ₀` | **0.1065** | freno de Loeper ALCANZABLE |
| ídem con `I_max` recortado a 0.0016 | 3.4e-4 | freno **desconectado** (contraprueba B.1) |
| offset de reloj vs `/fapi/v1/time` | **0.279 s** | guarda 7 en 0.500 s |

**Sección D — las dos trampas, medidas con `periodo_implicito`**

| | período que `A_arm` codifica de verdad |
|---|---|
| correcto (`ω_ang = 2π·f`) | **40.00 s** (real 40 s, error 0.00 %) |
| sin el 2π | 251 s |
| con `ω_m` en 1/Ticks | 5 027 s |

**Efecto de los cuatro defectos corregidos**

| Defecto | Antes | Después |
|---|---|---|
| WS mudo + `ticker/price` rezagado | `Tr(P)` → **6e9**, NIS `nan`, 0 correcciones | `Tr(P)` ~10, acotado |
| corrección por ciclo en vez de por paquete | `y1` con ρ₁ = 0.87 (artefacto) | corrige 1 vez por paquete |
| dos bots sobre la misma memoria | **294** conmutaciones / 10 000 ciclos | **1** / 40 000 ciclos |
| `A_arm` sin rama de Taylor | `nan` irreversible en `P` | identidad exacta con `ω=0` |

**Sobrecosto del EAKF sombra:** 0.046–0.10 ms/ciclo contra un presupuesto de 0.3 y los
2.2 ms de Loeper+NMPC. No es un problema y no lo será.

### Modo LECTURA — la precondición, y lo que destapó

`IS_TESTNET` se sustituye por `MODO ∈ {LECTURA, TESTNET, MAINNET}`. El booleano decidía a
la vez tres cosas que debían moverse por separado (generador de precios, modelo de λ,
modelo de fills), y con una sola bandera no se podía pedir "precios reales de Mainnet pero
sin ejecutar", que es justo lo que las Secciones D y E necesitan.

En LECTURA no hay credenciales cargadas y `assert_ejecucion_permitida` es una **aserción
dura**, no un `if`: la ejecución tiene que ser imposible, no improbable. La contabilidad es
de papel (`riesgo.CuentaPapel`) sobre precios reales, y está marcada como tal — **no es el
`ACCOUNT_UPDATE`** que la Sec. B.5 exige en TESTNET y MAINNET.

### ⚠ El WebSocket de futuros conecta y NO entrega datos

Fallo real, medido, y del peor tipo. Desde esta red `wss://fstream.binance.com` (futuros)
**completa el handshake y luego no manda ni un solo mensaje**, mientras
`wss://stream.binance.com:9443` (spot) funciona con normalidad. Sin excepción, sin cierre,
sin error.

Eso **derrota al watchdog de la Sec. 8.4.2 tal como estaba escrito**, porque ese watchdog
solo reacciona a excepciones y un socket sano que no habla no levanta ninguna. El efecto
completo se vio en la primera corrida de verificación: `valido=False` en todos los ciclos,
NIS en `nan`, burn-in eternamente en racha 0 y **Tr(P) creciendo hasta 6e9** sin una sola
corrección. El bot se creía en un mercado en calma.

Defensa en dos capas, ambas en `mercado.FeedPublico`:
1. **Detector de estancamiento**: 12 s sin mensaje sobre un socket conectado levanta
   `EstancamientoFeed`, que el watchdog sí ve.
2. **Degradación a sondeo REST** tras dos estancamientos. El objeto del feed se conserva
   entre reconexiones; si se reconstruyera, el contador volvería a cero y el bot
   reintentaría eternamente un socket ya demostrado mudo.

### ⚠ El endpoint REST obvio era el equivocado, por 7.6 segundos

Al degradar hay que elegir endpoint. Medido el 2026-08-04, con `τ_d` = (reloj local al
recibir) − (timestamp del propio dato), que es el retardo del sensor de la Sec. 6.5:

| endpoint | peso | τ_d mediana | volumen |
|---|---|---|---|
| `ticker/price` | 1 | **7.64 s** | no |
| `ticker/bookTicker` | 2 | 1.54 s | no |
| **`aggTrades`** | 20 | **0.70 s** | **sí** |

`ticker/price` era la elección obvia y es la equivocada: su campo `time` no es la hora del
servidor sino la del último cambio de precio publicado, y llega con ~7.6 s de rezago.
Contra el `τ_max = 2 s` de la Sec. 6.5 eso significa **rechazar el 100 % de los paquetes**,
con el mismo cuadro clínico que el socket mudo. Y el síntoma no apuntaba al endpoint por
ningún lado.

Se usa `aggTrades`, que además trae volumen — sin él ΣQ queda en cero y con ello mueren Φ,
Ψ y Ω, o sea todo el acoplamiento endógeno de la Sec. 1.4. Peso 20 a 1 Hz = 1 200/min
contra 2 400 disponibles.

Detalle que decide: se usa `aiohttp` y no `urllib` porque la sesión **reutiliza la conexión
TLS** — 0.27 s por llamada contra ~0.9 s abriendo socket nuevo cada vez. A 1 Hz de sondeo,
esa diferencia es la que decide si `τ_d` cabe bajo `τ_max`.

Offset de reloj medido contra `/fapi/v1/time`: **0.279 s**, coherente con los ~300 ms desde
Colombia que anticipa la Sec. 8, y holgado contra la guarda 7 (500 ms).

### ⚠ El filtro corregía 90 veces con la misma medición — y eso explica el ρ₁ de `y1`

`DIVERGE DEL PDF (Sec. 7.3)`: el filtro corrige **una vez por paquete nuevo**, no una vez
por ciclo de control.

El Hilo Rápido corre a ~90 Hz; el feed entrega a 2 Hz (R_n desde el Hilo Lento) y hasta
1 Hz (precio, con el feed degradado). Corrigiendo cada ciclo, la misma medición entra al
filtro decenas de veces: P se contrae como si hubiera decenas de observaciones
independientes, el filtro se declara mucho más seguro de lo que está, y **la innovación
queda autocorrelacionada por construcción**.

Eso no es una hipótesis: **es la explicación del ρ₁ = 0.87 de `y1`** que la v1.2 dejó como
pregunta abierta para este documento. La Sec. E.3 lo resuelve excluyendo `y1` de la
compuerta; esto lo resuelve en el origen, y además impide que `y0` heredara el mismo
artefacto al degradar el feed — lo que habría contaminado justo la serie sobre la que se
decide `A_arm`.

Efecto medido al aplicarlo, con el mismo feed: **Tr(P) pasa de 6e9 a ~10**, acotado, y
`valido=True`. Con el paquete repetido el NIS ni siquiera existía.

**Consecuencia sobre la telemetría, y hay que tenerla presente al analizar.** Se sigue
registrando una fila por ciclo de control (~90 Hz), porque `tr_P`, `x_k` y `rama_A`
evolucionan en cada uno. Pero la innovación **solo existe cuando llegó un paquete**, y en
los demás ciclos vale cero por relleno. Con el feed degradado eso es el 98.7 % de las filas.

Sin marcar esa distinción, Ljung-Box y el NIS se calculan sobre una serie que es casi toda
ceros: **ρ₁ sale 0.0000 y el veredicto del A/B es basura con apariencia de dato** — se vio
tal cual antes de corregirlo. De ahí el campo `hay_medicion` y la función
`diagnostico.solo_observaciones`, por la que pasa todo análisis basado en la innovación.
`diagnostico.py` reporta ahora explícitamente cuántas observaciones reales hay y avisa si
bajan del 5 % de los ciclos, porque en ese régimen un NIS bajo no significa "filtro
conservador" sino "feed lento".

### ⚠ Dos bots compartiendo memoria en silencio — el chatter que no era chatter

Defecto real, y de los que enseñan a desconfiar del síntoma. Al medir la conmutación de la
rama de `A` sobre datos reales salieron **294 cambios en 10 000 ciclos**, con permanencia
mediana de 0.27 s. Eso es chatter de manual, y el sospechoso obvio era la histéresis de la
Sec. D.4.2. Pero la histéresis estaba bien: aguanta 200 oscilaciones dentro de la zona
muerta sin conmutar una sola vez.

La causa era que **había dos instancias del bot corriendo a la vez, y ambas escribían la
misma memoria compartida**. La recuperación de bloques huérfanos que la v1.2 añadió para
Windows (donde `unlink()` es un no-op) tiene un efecto de segundo orden que entonces no se
vio: si el bloque existe porque hay otro Micelio **vivo**, la segunda instancia se adjunta
a él en vez de fallar. Dos Hilos Lentos escribiendo el mismo seqlock, dos Hilos Rápidos
publicando en el mismo Ring Buffer.

Con una sola instancia, mismo mercado y misma configuración: **12 cambios en 20 000 ciclos**.

`verificar_instancia_unica` lo corta antes de reservar nada. La detección es por **latido**
y no por PID: un PID se recicla, y en Windows no hay forma barata y portable de preguntar si
un PID sigue siendo el mismo proceso; un timestamp refrescado cada segundo no tiene esa
ambigüedad. Si aborta, **no limpia la memoria compartida** — es de la otra instancia, que
sigue trabajando, y liberarla sería exactamente el daño que se evita.

### Sección A — los filtros del instrumento, medidos

`exchangeInfo` leído de ambos entornos. Los cuatro slots del bloque de hot-reloading nacen
en **cero como centinela**, no con un valor por defecto: olvidarse de leerlos falla
ruidosamente en vez de operar con un literal.

Antes se llamaba `apply_filters(u_c, 1e-5, 1e-5, 10.0, P_spot)` — los cuatro inventados y
los cuatro equivocados por órdenes de magnitud. Un `stepSize` de 1e-5 contra el real de
1e-3 hace que el `floor` sea casi la identidad, y entonces nada en las pruebas revela que
el controlador continuo se convierte en interruptor al llegar al exchange.

**Dos correcciones a las cifras de la Sec. A.1**, ambas materiales:

| entorno | stepSize | minQty | tickSize | minNotional | 1 lote a 63 800 |
|---|---|---|---|---|---|
| MAINNET | 0.001 | 0.001 | 0.10 | **50 USDT** | 63.8 USD |
| TESTNET | **0.0001** | 0.0001 | 0.10 | 50 USDT | 6.4 USD |

1. **minNotional es 50, no 100.** La primera fila de la tabla de A.1 ("ninguna orden es
   legal") no se cumple hoy. Pero el margen es estrecho: por debajo de BTC = 50 000 un lote
   deja de alcanzar el nocional mínimo y la orden legal más pequeña pasa a 0.002 BTC. Por
   eso `cantidad_minima_legal` depende del precio en vez de fijarse una vez.
2. ⚠ **Testnet es 10× MÁS FINO que Mainnet**, y es una trampa silenciosa: la guarda de
   resolución evaluada contra Testnet da 476 lotes y pasa cómodamente, contra Mainnet da 47
   y va mucho más justa. Calibrar contra Testnet produciría un sistema que funciona en
   pruebas y se degrada a interruptor en producción. **`verificar_resolucion_control` se
   evalúa siempre contra el stepSize de Mainnet**, y el orquestador publica el más grueso
   de los dos.

### Sección B — la abrazadera va aguas abajo, e `I_max` no se toca

Los números al arranque, con S real de 63 920:

```
nocional_max_orden    = 3000 USD (46.9 lotes; minimo exigido 40)
nocional_max_posicion = 32599 USD (DERIVADO de I_max con holgura 2%)
apalancamiento 5x -> margen 6520, colchon liquidacion 6259 contra cap 3000 (2.1x)
lambda*S^2*gamma_0    = 0.1065  -> freno de Loeper ALCANZABLE
```

El principio de separación de la Sec. B.1 está verificado por contraprueba en los tests:
con `I_max = 0.0016 BTC` el producto cae a 3.4e-4 y el freno de singularidad queda
**estructuralmente inalcanzable** — correrías 30 episodios validando un sistema con un
mecanismo de seguridad desconectado y nada lo reportaría. Por eso `I_max` se queda en 0.50
y el límite de cuenta vive en el Motor de Red, después del Ring Buffer y antes de firmar.

**El `mmr` no se puede leer sin credenciales**: `leverageBracket` es un endpoint firmado.
`mercado.leer_mmr` devuelve el valor por defecto **inflado por un factor de seguridad de
2×** y marca que fue asumido, no leído. La guarda queda conservadora ante la duda.

⚠ **Una desviación del clamp de la Sec. B.4, deliberada.** El documento escribe
`u ← min(u, nocional_max_posicion(S)/S − |inv|)` sin distinguir compra de venta. Tomado
literalmente, al llegar al tope ese término vale 0 y anula **ambas** componentes, incluida
la que REDUCE la posición: el sistema quedaría atrapado en el límite, incapaz de deshacer,
y **rompería la ruta de cierre del halt**, que necesita emitir exactamente esas órdenes. Se
aplica de forma direccional sobre el inventario resultante. Marcado
`# NOTA DE INTERPRETACION:` y con test propio.

Las 7 guardas están implementadas con `causa_halt` distinguible y detalle propio. La
semántica es **cerrar y parar, no congelar**; `RutaDeCierre` reintenta con backoff, escala
la alerta y **jamás reporta éxito sin posición plana confirmada** — con test que la fuerza
a fallar.

### Sección C — episodios, y por qué el faucet automático necesita compuertas

`episodios.py`. La máquina vive en el **Motor de Red**, no en el orquestador: la capa de
riesgo, la cuenta y la ruta de cierre ya están en ese proceso, y moverla fuera obligaría a
sincronizar el cierre por memoria compartida justo en el momento en que menos se puede
confiar en el estado.

`ARRANQUE → OPERANDO → CERRANDO → CERRADO → REAPROVISIONANDO → ARRANQUE`, con `DETENIDO`
terminal. Verificado que **no existe transición de salida de `DETENIDO` por software**: la
única asignación a `ARRANQUE` desde otro estado está en `sondear_equity`, que exige estar en
`REAPROVISIONANDO`.

Las **seis compuertas de la Sec. C.3** están probadas una por una, forzando cada condición
por separado y comprobando que bloquea *la suya y solo la suya*. La que más importa:
**solo se recarga sobre la guarda 1 (drawdown)**. Las guardas 2–7 son fallos de sistema, y
recargar sobre ellas es tapar el bug con dinero — hay un test que recorre las seis causas de
sistema y verifica que ninguna pasa.

Sobre el adaptador de faucet, el documento tiene razón en no fingir: **el faucet de Testnet
es una función de la interfaz web, no un endpoint de la API**. Lo que hace que el modo
automático no sea una ruta sin probar es que **manual y automático comparten el sondeo de
equity**; la única diferencia es quién provoca la recarga. `ReaprovisionadorAutomatico`
lleva aserción dura contra `MODO == MAINNET` **en el constructor**, no un `if` que devuelve
`False`: un `if` dejaría el sistema corriendo con un reaprovisionador que sobre dinero real
es un sinsentido peligroso.

**Reset limpio (Sec. C.5), y cómo cruza procesos.** Se hace al ABRIR episodio, no al cerrar,
para que el estado residual de un cierre fallido tampoco se herede. El Motor de Red reinicia
lo suyo (capa de riesgo, cuenta) y publica `id_episodio` en memoria compartida; el Hilo
Rápido y el Hilo Lento **detectan el cambio y reinician lo suyo por su cuenta** — burn-in,
NIS, `P`, inventario y rama de `A` en uno; ΣQ, `S_ref`, ventana del EMD y `ω_m` en el otro.
No hace falta señalización adicional porque el campo ya es monótono y lo escribe un único
productor.

`id_episodio` va también en `TELEM_DTYPE`, y `diagnostico.py --episodio=N` filtra por él:
mezclar episodios falsearía Ljung-Box y el NIS, porque cada uno arranca con su propio
transitorio de burn-in.

### Sección D — las dos trampas, y por qué hacen falta tests

`A_arm` sale de `s̈ = −ω²s`, cuya solución es `cos(ωt)`: ahí ω es **angular**. Pero
`hht.frecuencia_instantanea` divide por 2π y `omega_m_desde_hz` devuelve ciclos/tick — las
dos son frecuencias **ordinarias**. Y `ω_m` está en 1/Ticks mientras `Δt` está en segundos.

Ninguno de los dos errores produce excepción, `nan` ni log. Producen un modelo que "no
aporta", y la conclusión equivocada sería que la propuesta no sirve. Medido con
`periodo_implicito`, que invierte la construcción de la matriz:

| | período que `A_arm` codifica de verdad |
|---|---|
| correcto (`ω_ang = 2π·f`) | **40.00 s** (real 40 s, error 0.00 %) |
| sin el 2π | 251 s |
| con `ω_m` en 1/Ticks | 5 027 s |

Resolución: **dos variables publicadas por separado**, sin conversión en el punto de uso.
`ω_m` [1/Ticks] sigue alimentando ρ_k y c²_vol sin cambios; `ω_ang` [rad/s] alimenta
`A_arm` y solo `A_arm`.

La predicción pasa a **forma afín** (`x_pred = x_ref + A·(x_k − x_ref)`), que conserva
`x[0] = S` absoluto y evita auditar a todos los consumidores. Verificado: al saltar `S_ref`
100 USD en un nodo de fase, **P queda idéntica** y `x[0]` se mueve 0.00012 USD, contra los
100 USD que saltaría la formulación sobre la desviación.

Concentración espectral `C` medida: **0.967** en un ciclo nítido contra **0.487** con
energía repartida, así que `C_ON = 0.50` discrimina. La histéresis aguanta 200 oscilaciones
dentro de la zona muerta con **cero conmutaciones**.

Sobre mercado real, el conmutador entra en la rama armónica con `C` de 0.81–0.99 y períodos
de 47–188 s, y vuelve a velocidad constante al caer `C` por debajo de 0.35. Ojo con el
warm-up: la ventana del EMD son 384 muestras a 0.5 s, o sea **192 s** antes del primer
tamizado; hasta entonces `C = 0` y la rama es velocidad constante, que es lo correcto.

### Sección E — dos EAKF en paralelo

El sombra usa siempre la rama contraria a la de control y comparte `z_k`, `R_k` y `Q_k`:
cualquier otra diferencia contaminaría la comparación. Sobrecosto medido **0.046–0.10
ms/ciclo** contra un presupuesto de 0.3 y los 2.2 ms de Loeper+NMPC.

`diagnostico.py` reensambla las innovaciones de cada modelo cruzando por `rama_A` y aplica
la regla de decisión de la Sec. E.3 sobre `y0`, leyendo **mediana** de NIS. Si el armónico
sale peor, el reporte remite explícitamente a `test_D_trampa_2pi_periodo_implicito`, porque
un 2π o un factor 125 se ven exactamente así.

**Primera medición end-to-end (8 min de LECTURA, 355 observaciones):**

| modelo | NIS mediana | ρ₁ | ρ₂ | ρ₃ |
|---|---|---|---|---|
| velocidad constante | 1.166 | +0.268 | −0.037 | −0.089 |
| oscilador armónico | 1.237 | +0.483 | +0.360 | +0.326 |

⚠ **Esto NO es el veredicto, y el reporte ahora lo dice él mismo.** Se añadieron tres
condiciones de "no decidible" que se evalúan antes de aplicar la regla, porque un veredicto
emitido sobre datos que no lo sostienen se lee igual que uno bueno:

1. La corrida dura 0.13 h contra las ≥ 24 h que exige la Sec. E.2.
2. **0.76 mediciones/s contra ~90 Hz de ciclo de control.** Con el feed degradado, `A_arm`
   propaga el oscilador ~1.3 s entre correcciones, así que un error del 10 % en ω se
   acumula mucho más que en velocidad constante. Con estos datos el A/B mide el **feed**, no
   el modelo — y eso explica por sí solo que el armónico salga peor.
3. ρ₁ mínimo detectable con n = 355 es 0.296, o sea que el 0.268 de velocidad constante
   ni siquiera es distinguible del ruido muestral.

Lo que sí queda demostrado es que la infraestructura del A/B funciona de punta a punta y
que el veredicto, cuando llegue, será legible.

### Sección F — la suite de aceptación, y para qué sirve de verdad

`python tests_v13.py` → **33/33**. Sin dependencia de pytest a propósito: el runner son
veinte líneas, así que la suite corre en cualquier entorno donde corra el bot.
`--sin-red` omite los dos que consultan `exchangeInfo`.

Cada test lleva en el nombre la casilla de la Sección F que cubre, y **cada uno imprime el
número que midió**, no un "OK" pelado — la salida es en sí misma el registro de calibración.

Los tres que más valen:

- `test_B_disparo_forzado_cadena_completa` — el criterio que el propio documento marca como
  el más importante. Recorre **detección → cierre → posición plana confirmada → halt →
  alerta → volcado del resumen** inyectando una pérdida ficticia. "Un freno que nunca se
  probó no es un freno."
- `test_B_cierre_fallido_reintenta_escala_y_no_miente` — fuerza el cierre a fallar siempre y
  verifica que reintenta, escala y **propaga la excepción**. Reportar un halt como
  completado con exposición abierta es el fallo más caro del sistema, porque a partir de ahí
  nadie está mirando.
- `test_D_trampa_2pi_periodo_implicito` — la única defensa contra D.1 y D.2, que no producen
  excepción, ni `nan`, ni log. Sin este test, un factor 2π o un factor 125 se leerían como
  "el oscilador armónico no aporta".

⚠ **Cuatro de los cinco fallos iniciales de la suite eran defectos de los tests, no del
código.** El más instructivo: en `d[k] = v` Python evalúa `v` **antes** que `k`, así que
`disparadas[capa.evaluar(S)] = capa.detalle_halt` leía el detalle antes de que la guarda
disparara. Un test que falla por su propia culpa gasta el mismo tiempo que uno real.

### Pendiente tras esta fase

- **El A/B todavía no tiene veredicto.** La Sec. E.2 exige ≥ 24 h continuas de Modo LECTURA
  cubriendo las sesiones asiática, europea y americana. Lo ejecutado es la infraestructura y
  su verificación, no la corrida de decisión.
- **`C_ON` y `C_OFF` siguen `[CALIBRAR]`** contra la distribución real de `C` en Mainnet.
  Los 0.50/0.35 son los sugeridos por el documento; el 0.967/0.487 medido es sobre señal
  sintética, no sobre mercado.
- **El WebSocket de futuros sigue mudo desde esta red, y es hoy el cuello de botella
  principal.** El bot opera degradado a 0.76 mediciones/s contra ~90 Hz de ciclo de control,
  y eso basta para invalidar el A/B por sí solo (ver arriba). El spot sí funciona, así que
  parece filtrado de `fstream.binance.com` y no un problema de código. **Resolver esto es la
  tarea de mayor rendimiento pendiente**: sin feed rápido, ni la Sección E ni la Fase 2
  pueden concluir nada.
- **Fase 2 (ALS)**: sigue tras su compuerta, pero con el artefacto de la medición repetida
  eliminado, el ρ₁ que se mida ahora sí es evidencia sobre `A`.
- **`r_S,base` y `r_EMD`** siguen medidos mal por 8× y 161× (v1.2). No se tocaron: la Fase 2
  es su dueña, y la regla de no inflar Q ni R con una constante ad-hoc sigue en pie.
- Fase 3 (`Ω_crit`) — necesita Testnet operativo.
- `ω_m,max` en `constantes_micelio` sigue siendo un `TODO(HHT)`. Ahora hay material para
  fijarlo: sobre mercado real se observaron períodos de 33–188 s.
- **`PRECIO_REFERENCIA` se actualizó a 63 000** (estaba en 45 000 desde la v1.1). Se
  comprobó que **no altera γ_Q**, porque su denominador es `max(ΣQ_max, K_USD)` y ΣQ_max =
  1.3e5 sigue dominando. Queda marcado como ancla que se desactualiza sola con el precio.
- `.gitignore` excluye ahora `episodios/` — los resúmenes y diagnósticos por episodio se
  regeneran en cada corrida.
- Hay un `Addendum_I_Cancion_del_Micelio.pdf` sin rastrear en el árbol. **No lo he tocado**;
  decide tú si entra al repo.

### Orden sugerido para la v1.4

1. **Desatascar el feed de futuros.** Todo lo demás depende de ello y hoy nada más lo hace.
2. Corrida de LECTURA de ≥ 24 h y veredicto del A/B (Sec. E.3).
3. Según ese veredicto: Fase 2 (ALS) sobre `r_S,base` y `r_EMD`, o reabrir el diagnóstico
   de `A` si ambos modelos siguen con ρ₁ ≥ 0.20.
4. Testnet con credenciales → Fase 3 (`Ω_crit`) y las 30 corridas, con el test de disparo
   forzado ejecutado en cada una (Sec. G).

## Sesión 2026-08-05 — Relojes y reconstrucción del acoplamiento (v2.0)

Ejecuta `ORDEN_TRABAJO_RELOJES_2_0.md`: precondición de medición y §2 a §8. Es una tanda de
**reorganización estructural**, no de calibración: cambia dónde vive el tiempo en el sistema.
Suite: **42/42** (`python tests_v13.py`), los 33 de la v1.3 más 9 nuevos.

### RESULTADOS DE UN VISTAZO

**La tesis del documento, confirmada por medición**

| | valor |
|---|---|
| lotes de `aggTrades` a 1 Hz | p50 = 3, p95 = 6, máx = 1000 (el límite del API saturando) |
| factor de lote | **24.9×** |
| **transacciones descartadas por tener `P_spot` escalar** | **96.0 %** |

Las "0.76 mediciones/s" de la v1.3 eran la tasa de **paquetes**, no la de información. El
sistema no estaba escaso de datos: los estaba tirando.

**Observaciones que llegan al filtro, por etapa**

| etapa | obs/s | % de ticks con medición |
|---|---|---|
| v1.3 | 0.76 | 0.90 % |
| anillo publicando, sin consumir | 2.95 | 3.2 % |
| **§4 completo (Δn = 1)** | **24–94** (sigue al mercado) | **100 %** |

**Autocorrelación de la innovación, antes y después**

| canal | v1.3 | v2.0 §4 | v2.0 + multi-tasa |
|---|---|---|---|
| `y0` (precio) | +0.268 | **−0.164** | −0.164 (*rechazo no material*) |
| `y1` (residuo) | +0.819 | +0.819 | **+0.476** |
| repetición de `y1` | 11.6× | 11.6× | **1.5×** |

### ⚠ El WebSocket de futuros NO estaba filtrado — la v1.3 concluyó mal

La v1.3 dio por filtrada la ruta porque `btcusdt@aggTrade` conectaba y callaba mientras spot
funcionaba. Sondeado en serio, en el **mismo host y el mismo socket**:

| stream en `fstream.binance.com` | |
|---|---|
| `@bookTicker` | 91–288 msg/s |
| `@trade` | 26–518 msg/s |
| `@depth@100ms` | 9.8 msg/s |
| `@aggTrade`, `@markPrice`, `@kline`, `@ticker` | **mudos** |

Y con suscripción explícita por mensaje el servidor **confirma**
`{"result":["btcusdt@aggTrade"],"id":99}` y aun así no manda un dato. El socket es
bidireccional, Binance nos oye y nos contesta; los certificados son legítimos (DigiCert) y no
hay proxy. **Descartado el filtrado de ISP, y con él la necesidad de montar un VPS.**

La causa última del silencio de `@aggTrade` **queda sin explicar**, y se deja anotada como
tal. Lo operativo es que `@trade` da lo mismo sin agregar —precio, cantidad, hora del
exchange, identidad y lado del taker— con τ_d = 0.198 s (p90 0.274) y **cero huecos en 1170
ids consecutivos**. Es ahora el stream primario (`mercado.STREAM_TRANSACCIONES`).

### ⚠ La tasa de transacciones varía por un factor 20 — cuidado con "lo típico"

Medido el mismo día sobre el mismo par: **18.8 tx/s** (REST, momento tranquilo), **26–32 tx/s**
(WS), **517.9 tx/s** media hora después. Cualquier constante calibrada contra "la tasa típica"
es sospechosa de origen. Por eso `K` (muestreo del EMD) y el avance del reloj se derivan de ν
en ejecución, y el anillo se dimensiona contra el **pico**: 16 384 entradas ≈ 31 s de colchón
a 518 tx/s.

### §3 — La ingesta por lotes

`MERCADO_DTYPE` es un Ring Buffer SPSC con **identidad de trade**, que es lo que hasta ahora
faltaba y lo que habilita las dos cosas siguientes:

- **Deduplicación por `aggTradeId`** (§3.3). El sondeo REST devuelve ventanas solapadas; sin
  deduplicar, los mismos trades entran al filtro varias veces — el bug de las 90 correcciones
  con otro disfraz. El contador **no se reinicia entre reconexiones**, a propósito.
- **Detección de huecos** (§3.4). Un modo de fallo nuevo y distinto del socket mudo: el feed
  funciona y aun así falta información. Se publican `n_huecos` y `trades_perdidos`, y se marca
  el instante para poder excluir el tramo del análisis.

### §4 — El filtro bajo Δn = 1, y lo que elimina estructuralmente

Un tick = un paso de predicción + una corrección. No parchea: vuelve **imposible de
expresar** corregir dos veces con la misma medición, descartar el resto del lote, y —lo más
sutil— meter al planificador dentro de `Q`.

Ese último era un defecto real: `Q` no dependía de Δt, así que `Σ AⁱQAⁱᵀ` crecía con el
**número de pasos** y no con el tiempo. El ruido de proceso era proporcional a cuántas veces
despertó el planificador, y `q_base` estaba calibrada en silencio contra ~90 Hz. Ahora
`q_S_tick = (σ_rel_tick·S)²` con `σ_rel_tick = σ_rel_s/√ν` — varianza por **transacción**,
que es una propiedad del mercado. El test de invariancia a la tasa lo cierra: N pasos de
Δn=1 dan el mismo `x` y la misma `P` que un paso de Δn=N.

⚠ Ese test exige `Q_N = Σ AⁱQAⁱᵀ`, **no `N·Q`**. Solo coinciden si `A = I`; en cuanto `A`
propaga, suponer `N·Q` sobreestima la certeza en la posición. `dinamica.acumular_Q` existe
para eso, y el test comprueba explícitamente que la versión ingenua **no** pasa — si pasara,
el test no discriminaría y no valdría.

**La trampa del 2π reaparece en espacio de ticks**, intacta. Medida con
`periodo_implicito_ticks` sobre un ciclo de 795 ticks: correcto **795.0**, sin el 2π **4995**.

**Auditoría de `x[1]` tras el cambio de unidades** (USD/BTC por segundo → por transacción):
el único consumidor es la telemetría (`reg["x1"]`). Ninguna fórmula la consume, así que el
cambio es seguro. La constante que la escala, `q_base·1e-2`, sí queda mal calibrada y pasa a
ser `[CALIBRAR]` de la Fase 2.

### ⚠ El §4 arregló medio problema — y medirlo lo delató

Con Δn = 1, `y0` bajó a ρ₁ = −0.164. Pero `y1` seguía en **+0.819**, y el propio reporte
explicaba por qué: `R_n` cambiaba en el 8.6 % de los pasos y se repetía **11.6×**. Era el bug
de las 90 correcciones **desplazado de canal** — el filtro pasó a decenas de ticks/s mientras
el Hilo Lento sigue publicando `R_n` a 2 Hz. Tercera aparición de la misma familia en tres
sesiones, que es justo contra lo que advierte el §3.3.

Cerrado con **actualización secuencial multi-tasa**: si `R_n` es fresco, `H` completa y m = 2;
si no, solo la fila del precio y m = 1. Asimilar un `R_n` retenido contrae `P` con información
que no existe. La frescura se detecta con el contador del seqlock, no comparando valores —
comparar valores confundiría "no cambió" con "no llegó".

En telemetría, `y1` va a **NaN** cuando no se observó, no a cero: cero significaría
"innovación nula", que es lo contrario. Resultado: ρ₁ de 0.819 → **0.476** y repetición de
11.6× → 1.5×. Lo que queda es la correlación genuina de `R_n` por ventanas EMD solapadas, que
la Sec. E.3 de la v1.3 ya excluye de la compuerta.

### El signo de ρ₁ cambió, y eso dice algo

`y0` pasó de **+0.268 a −0.164**, con decaimiento rápido (ρ₁₀ = −0.022). Un ρ₁ negativo a
rezago 1 es la firma del **rebote bid-ask**: microestructura, no desajuste de `A`. Es
exactamente lo que el §4.3 anticipa al procesar por transacción, y —como manda— **se registra,
no se compensa**. La Fase 2 (ALS) sigue siendo la dueña de `r_S,base`, ahora con un objetivo
bien definido: ruido de observación **por transacción**, que sí es una cantidad con
significado físico.

Dato relacionado: la curtosis de `y0` es **+207**. Con datos por transacción la mayoría de los
trades no mueven el precio y unos pocos saltan, así que la cola es enorme. Refuerza la nota de
la v1.2 sobre ALS-IRLS (Huber) — pero primero el ALS estándar, que sin baseline no hay mejora
que medir.

### §5 — El EMD dejaba de tamizar una escalera

El buffer se llenaba con reloj de pared: `P_spot` cada 0.5 s "hubiera cambiado o no". Con el
feed lento, media ventana eran **duplicados literales**. El comentario que había ahí se
preocupaba del espaciado irregular y resolvió la uniformidad de la malla *temporal* justo
mientras la de *valores* se volvía escalera — y la transformada de Hilbert de un escalón tiene
contenido en todo el espectro. Ahora se muestrea **cada K transacciones**, lo que da malla
uniforme en ticks y cero duplicados a la vez; en reloj de ticks el jitter del planificador no
existe, así que la preocupación original desaparece en vez de resolverse.

### ⚠ §5.2 — CONTAMINACIÓN CONFIRMADA: la sincronización era LOCAL, no global

Ejecutado sobre **11 870 transacciones reales** capturadas en 439 s (ν = 27.0 tx/s,
K = 14 ticks/muestra), tamizando la misma serie por tres vías con W = 384 en las tres.
Se añadió una tercera vía a las dos que pide el documento, precisamente para poder
**separar las dos causas posibles**:

| vía | duplicados | f_hz | período | C |
|---|---|---|---|---|
| **A** escalera (v1.3): pared 0.5 s sobre precio a 1 Hz | 74.2 % | 0.01051 | 95.2 s | **0.917** |
| **B** pared 0.5 s con precio siempre fresco | 50.1 % | 0.00381 | 262.6 s | 0.737 |
| **C** ticks (v2.0): una muestra cada K transacciones | 44.1 % | 0.07587 | **13.2 s** | 0.624 |

```
A vs C (v1.3 contra v2.0):  f_hz difiere  86.2 %    C difiere 0.294
A vs B (efecto DUPLICADOS): f_hz difiere  63.8 %
B vs C (efecto ESPACIO):    f_hz difiere  95.0 %    C difiere 0.114
```

**Veredicto: contaminación material, y la causa dominante es el ESPACIO de muestreo**
(95.0 %), no los duplicados (63.8 %). Sin la vía B esto no se podría afirmar — un
resultado positivo A-vs-C no diría cuál de las dos cosas hay que arreglar.

Consecuencias, en orden de gravedad:

1. **El A/B de la v1.3 queda anulado por una segunda causa independiente.** No solo
   corría con 0.76 mediciones/s: además `ω_m` venía de tamizar una escalera. Un período
   de 95 s donde el tiempo de transacción dice 13 s es un factor **7**.
2. **`C` estaba inflada por la escalera, y `C` es lo que decide la rama de `A`.** La
   escalera suprime las IMFs rápidas donde el valor se repite, así que la energía se
   concentra en el modo lento y `C` sube: 0.917 con 74 % de duplicados contra 0.624 con
   44 %. Los `C` de 0.94–1.00 que en la v1.3 mantenían enganchado el oscilador armónico
   **eran en buena parte artefacto**. Con muestreo limpio `C` sigue por encima de
   `C_ON = 0.50`, pero con mucho menos margen — y alimentando una frecuencia 7× distinta.
3. **Respuesta a la pregunta de fondo: la sincronización era LOCAL.** El §4 puso el
   filtro en reloj de transacciones, pero el EMD —el otro consumidor de datos de mercado,
   y el productor de `ω_m` y `C`— seguía en reloj de pared. La cadena estaba sincronizada
   a medias. Ahora lo está de extremo a extremo.

⚠ **Lo que el muestreo por ticks NO arregla:** la vía C conserva un **44.1 %** de
duplicados. No son artificiales: a 27 tx/s y K = 14, el precio de BTC genuinamente no se
mueve en 14 transacciones (tickSize 0.10 USD). El muestreo por ticks elimina los
duplicados por *sobremuestrear una variable rancia*, no los que trae el mercado. Queda
como límite conocido del estimador, no como defecto pendiente.

Reproducible con `python test_contaminacion_emd.py --capturar=440`; la captura queda en
`telemetria/captura_trades.npz`. La versión determinista, con verdad conocida y sin red,
está en la suite (`test_v20_contaminacion_emd_reloj_de_pared`).

**Un solo tramo de mercado no generaliza.** Conviene repetirlo en régimen agitado antes
de dar por cerrada la magnitud del efecto.

### §6 — El reloj de volumen

`ΔQ*` **medido**: mediana de volumen por transacción = 0.0060 BTC ≈ **386 USD** a S ≈ 64 400.
Con eso, a actividad típica el reloj de volumen y el de ticks avanzan a la misma tasa media y
lo único que los separa es la ponderación por tamaño — que es la propiedad a poner a prueba.

`Q_acumulado_total` es monótono y **distinto de ΣQ**, que se reinicia en cada nodo de fase.
φ' se aplica al refractario de nodos, el caso que el documento llama más claro: no quieres dos
nodos separados por poco *tiempo*, quieres que lo estén por poco *mercado*.

La regla §6.5 está codificada en `CONSUMIDOR_DE_OMEGA` y tiene test propio que **no es vacuo**:
fuerza la variante equivocada y comprueba que el resultado difiere (170 %).

### Defecto operativo que la v2.0 introdujo, y su arreglo

Con una fila de telemetría por transacción y la tasa variando 20×, el volcado solo-al-llenar
hacía impredecible cuándo hay datos en disco (20 s a 518 tx/s, 6 min a 26 tx/s) y una parada
no limpia se llevaba todo lo acumulado. Costó varias corridas de verificación enteras. Añadido
volcado **también por tiempo** (60 s), con bloque parcial.

### TABLA DE MEDICIONES — todos los números crudos de la v2.0

Recogidos aquí para poder analizarlos sin releer la narrativa. Todos son de ejecución real
salvo donde diga sintético.

**Feed y latencia** (WebSocket `btcusdt@trade`, Mainnet, 45 s)

| magnitud | valor |
|---|---|
| τ_d | p10 0.189 s · **p50 0.198 s** · p90 0.274 s · max 0.522 s |
| fracción que supera τ_max = 2.0 s | **0.0 %** |
| huecos en ids consecutivos | **0 de 1170** |
| cobertura temporal | 44.8 s de mercado en 45 s de reloj = **0.99×** (tiempo real, no replay) |
| llegada en ráfagas | 85.7 % de los mensajes a < 2 ms del anterior; ráfaga p90 = 35, max = 203 |
| separación real entre trades | p50 0.0 ms · p90 74.0 ms |
| offset de reloj vs `/fapi/v1/time` | 0.279 s |

**Endpoints REST comparados** (sesión reutilizada, `τ_d` = local al recibir − timestamp del dato)

| endpoint | peso | τ_d mediana | trae volumen |
|---|---|---|---|
| `ticker/price` | 1 | **7.64 s** ← el obvio, y el equivocado | no |
| `ticker/bookTicker` | 2 | 1.54 s | no |
| **`aggTrades`** | 20 | **0.70 s** | **sí** |

Latencia por llamada: **0.27 s** con sesión `aiohttp` reutilizada contra **~0.9 s** abriendo
socket nuevo con `urllib` en cada sondeo.

**Lote y volumen** (`aggTrades` a 1 Hz, 46 llamadas en 61 s)

| magnitud | valor |
|---|---|
| transacciones nuevas por lote | p50 = 3 · p95 = 6 · **max = 1000** (el límite del API saturando) |
| media por lote | 24.9 |
| tasa real | 18.8 tx/s |
| **descarte por `P_spot` escalar** | **96.0 %** |
| volumen por transacción | mediana **0.0060 BTC** · p95 1.0050 · max 30.1420 |
| ΔQ* resultante a S ≈ 64 400 | **386.40 USD** |

**Ljung-Box tras Δn = 1 + multi-tasa** (5857 ticks, 240.9 s, 24.3 ticks/s)

| canal | ρ₁ | ρ₂ | ρ₃ | ρ₁₀ | ρ₄₅ | cambia en | veredicto |
|---|---|---|---|---|---|---|---|
| `y0` precio | **−0.164** | −0.077 | −0.097 | −0.022 | −0.001 | 100.0 % | rechazo NO material |
| `y1` residuo | +0.476 | +0.220 | −0.100 | +0.111 | +0.004 | 66.4 % | rechaza blancura |

Antes de la multi-tasa, `y1` daba ρ₁ = +0.819, ρ₂ = +0.723, ρ₃ = +0.555, cambiando en el
8.6 % de los pasos (repetición **11.6×**).

**Shapiro-Wilk** — colas, y por qué importan para la Fase 2

| canal | W | curtosis exceso | asimetría |
|---|---|---|---|
| `y0` precio | 0.09993 | **+207.8** | −9.52 |
| `y1` residuo | 0.54108 | +14.09 | +0.79 |

Con datos por transacción la mayoría de los trades no mueven el precio y unos pocos saltan:
esa cola es real, no un outlier a recortar.

**Latencia del lazo**

| magnitud | v1.3 | v2.0 |
|---|---|---|
| cadencia del Hilo Rápido | 11.0–11.9 ms/ciclo | **41.1 ms/ciclo** |
| sobrecosto del EAKF sombra | 0.046–0.11 ms/ciclo | 0.10–0.11 ms/ciclo (presupuesto 0.3) |

⚠ La cadencia subió 3.7× porque cada ciclo drena el lote de transacciones. Sigue por debajo
del presupuesto para el trabajo útil, pero **hay que vigilarlo en régimen de 518 tx/s**: es el
número que decide si el drenado necesita su propio hilo.

**Salidas de los tests que llevan número**

| test | resultado |
|---|---|
| `periodo_implicito_ticks` | correcto **795.0** ticks · sin 2π **4995** (ciclo real de 795) |
| invariancia a la tasa | 12 pasos ≡ 1 paso de Δn=12; la versión ingenua `N·Q` se desvía 5.02e-03 |
| regla §6.5 | Φ con la ω equivocada se desvía **170 %** |
| sobrepaso del anillo | lapeo de 500 detectado y contabilizado en 1 sobrepaso |
| §5.2 sintético | ticks recupera 100 muestras/ciclo (verdad 75, error 25 %); pared da 194 y difiere **48 %** |
| σ_rel por tick | σ_s/√ν = 4.354e-07 con ν = 26 tx/s |

### Pendiente tras esta fase

- **El A/B sigue sin veredicto.** El test de contaminación del §5.2 ya está ejecutado y salió
  POSITIVO (ver arriba), así que el resultado de la v1.3 queda anulado por partida doble.
  Falta la corrida de ≥ 24 h (§E.2), ahora sobre la cadena ya sincronizada de extremo a extremo.
- **Fase 2 (ALS)** sigue tras su compuerta, pero el terreno cambió: `y0` ya está por debajo del
  umbral material y lo que queda es microestructura, que es lo que ALS debe estimar.
- El silencio de `@aggTrade` en fstream sigue **sin explicación**. No bloquea nada.
- `C_ON`/`C_OFF`, `ω_m,max`, `ΔS_ref`, `ΔQ*` y `q_base·1e-2` siguen `[CALIBRAR]`.
- Fase 3 (`Ω_crit`) y las 30 corridas de Testnet, sin cambios.

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
