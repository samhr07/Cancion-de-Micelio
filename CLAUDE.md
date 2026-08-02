# Canción del Micelio — Contexto del Proyecto

Bot de trading algorítmico. Autor: Samuel Hoyos R. Documento de diseño: `Trading_Bot.pdf` (Julio 2026).
Fuente principal: `Micelio.py` (orquestador consolidado).

**Este archivo recoge un diagnóstico previo hecho sobre `Micelio.py` contrastándolo con el PDF.
Los tests citados se corrieron en un contenedor aislado sin GPU (NumPy/SciPy/CasADi 3.7.2 con
IPOPT y qpOASES). No se han validado contra CUDA ni acados.**

---

## ESTADO ACTUAL (2026-08-02)

El orden de trabajo de abajo está **ejecutado en su totalidad**. Ver la sección
"Sesión 2026-08-02" al final para el detalle de qué se corrigió, cómo se verificó, y los
huecos NUEVOS del PDF que aparecieron al hacerlo. Las secciones intermedias se conservan
como registro del diagnóstico original — describen el estado *anterior* del código.

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

### ⚠ ε_burn NO es una constante estática

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

## Convenciones

- Comentarios y nombres de variables en español, consistente con el código y el PDF existentes.
- Referenciar la sección del PDF en los comentarios al implementar una fórmula.
- Cualquier suposición que rellene un hueco del PDF debe marcarse explícitamente con
  `# NOTA DE INTERPRETACION:` y describir qué se asumió.
