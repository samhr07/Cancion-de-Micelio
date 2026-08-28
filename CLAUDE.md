# Canción del Micelio — Contexto del Proyecto

Bot de trading algorítmico. Autor: Samuel Hoyos R. Documento de diseño: `Trading_Bot.pdf` (Julio 2026).
Fuente principal: `Micelio.py` (orquestador consolidado).

**Este archivo recoge un diagnóstico previo hecho sobre `Micelio.py` contrastándolo con el PDF.
Los tests citados se corrieron en un contenedor aislado sin GPU (NumPy/SciPy/CasADi 3.7.2 con
IPOPT y qpOASES). No se han validado contra CUDA ni acados.**

---

## ESTADO ACTUAL (2026-08-09)

### Lo primero que hay que saber

**La v3.2 ya se ejecutó y su regla de decisión se paró en el paso 3.** Queda una línea abierta.
No se abre una tercera.

| línea | qué decide | estado |
|---|---|---|
| **v3.2** `ORDEN_TRABAJO_MIGRACION_3_2` | impacto permanente contra transitorio (M0/M1/M1′/M2) | **EJECUTADA** el 2026-08-10 sobre 1 031 155 ticks continuos. Paso 2 pasa, **paso 3 falla**: `q90(\|μ̂\|) = 11.30` contra `1.5·c(u) = 39.05`. Falta el §7 (exige captura completa) |
| **v4.0** `ORDEN_TRABAJO_EJECUCION_4_0` | cuánto cuesta operar | `captura_estacional` corriendo, 21 días. §4.2 y §5 ya resueltos |
| **v4.1** `ORDEN_TRABAJO_HORIZONTE_4_1` | a qué horizonte (si a alguno) la señal paga el peaje | §11, §2, §3 hechos. El §1 se rehizo el 2026-08-23 y su veredicto quedó **RETRACTADO el mismo día**: `σ₁` y `ν` se agruparon sobre las 24 h y el requisito varía **8.4×** según la hora. **Pendiente: §1 estratificado por casilla horaria.** Quedan §6, §5, §4 |

✅ **Estado de la captura estacional (2026-08-22): SU COMPUERTA PASA.** Es la primera vez.
**32 429 468 transacciones y 223 546 502 snapshots de libro, 0 partes ilegibles** de 3 207, 2.0 GB.
Cobertura **15 522 min = 10.78 días equivalentes**, **mínimo por casilla 60** contra los ≥ 30 que
pide la compuerta, **0 casillas vacías de 168**. Calidad: 0 precios ≤ 0, 0 ids duplicados, 0
tiempos no monótonos, `maker` = 49.29 % (no degenerado). El proceso lleva **240 h sin reiniciarse**
(PID 17332, desde el 2026-08-12) y el vigilante no ha tenido que reponerlo ni una vez.

**Dos tramos continuos largos, y son de REGÍMENES OPUESTOS** — que es lo que este proyecto
llevaba pidiendo desde que la v3.2 encontró que nada replica entre capturas con `ν` distinta:

| tramo | horas | ticks | `ν` | recorrido de precio | ventanas limpias 15/30/60/120/240 min |
|---|---|---|---|---|---|
| 2026-08-12 11:46 → 08-17 23:36 UTC | **131.83** | 8 898 312 | 18.75 tx/s | 62 484 – 64 601 (**3.4 %**, lateral) | 513 / 251 / 120 / 55 / 23 |
| 2026-08-18 01:26 → 08-21 11:54 UTC | **82.47** | 19 985 140 | **67.31 tx/s** | 63 979 – **79 555** (**24.3 %**, tendencia) | 325 / 160 / 78 / 37 / 18 |

`ν` difiere **3.6×** entre los dos y el recorrido de precio **7×**. Ventana limpia = sin ningún
hueco > 10 s. Ids ausentes 1.04 % y 1.56 % respectivamente.

⚠ **Hueco de 24.26 h el 2026-08-21 (06:54 → 07:09 local del día siguiente), y la causa es NUEVA.**
Reconstruido del registro de eventos: 06:54:02 `Kernel-Power 105` **cambio de fuente de energía**
(se desenchufó) + `506` entrada en espera moderna → 06:54:55 `172` **«Conectividad: Disconnected.
Motivo: Policy Setting»** → 07:27:40 `42` suspensión → hueco en el propio registro → 08-22
07:09:26 eventos `Kernel-Boot` y salto del reloj de `08-21T12:27:44` a `08-22T12:09:26`. Eso es
**hibernación (S4)**, no apagado: el PID sobrevivió intacto y volvió a volcar solo.
**Corrige la regla que teníamos escrita:** la espera moderna es inofensiva **enchufado** — probado
en 8 ciclos — pero **en batería Windows desconecta la red por directiva**. Segunda vez que la
batería cuesta datos (la primera, el 2026-08-10, por umbral crítico). **Mantener el portátil
enchufado es la única medida que ha hecho falta y la única que faltaba.**

⚠ **La cobertura ya no bloquea nada. Rehacer el §1 de la v4.1 sobre estos dos tramos es lo
siguiente**, y por primera vez con potencia: contra las 8/4/2 ventanas con las que el §1 se declaró
NO DECIDIBLE, ahora hay **838 a 15 min, 411 a 30 min y 198 a 1 h** sumando los dos tramos. Y el
tramo en tendencia es justo el régimen donde un `R²` no nulo es plausible.

⚠ **Lo que la v3.2 dejó sin decidir NO es por falta de datos**: el estimador de `D` **no tiene
potencia** a la autocorrelación de signos real (`γ̂ = 0.798`) — bajo `D = 1` verdadero devuelve
cualquier cosa entre 0.0 y 7.3. Capturar más no lo arregla; hace falta otro estadístico.

**`Micelio.py` no se toca desde la v2.2.** Todo lo posterior es código de análisis aparte.

**Preregistros vigentes y congelados:** `PREREGISTRO_3_2.md` (9 enmiendas, todas anteriores a
mirar dato) y `PREREGISTRO_4_0.md` (1 enmienda). Los dos llevan registro de enmiendas con hash
antes/después y constancia de si había resultado a la vista.

**Suites:** `tests_v13.py` 56/56 · `ssa.py` 11/11 · `migracion_v32.py` **20/20** · `cola.py` 9/9 ·
`difusividad.py` 5/5 · `tick_grande.py` **28/28** · `horizonte.py` **14/14**.

⚠ **Retractaciones vigentes — no citar lo retirado:**
- ⚠ **No existe «la» `γ` de este mercado** (v4.1 §3.1, 2026-08-12). Estimada por GPH y Whittle
  local, `γ̂` recorre de **+0.69 a −0.26** según el ancho de banda `m = N^α`, en las dos capturas.
  El mismo barrido sobre fGn de `γ` conocida es **plano** (salto máximo 0.06–0.13 contra 0.45–0.48
  sobre el dato), así que el fallo es de la serie: **el flujo de órdenes no tiene un régimen de
  escala único.** Cualquier fórmula que use `γ` hereda esa indeterminación.
- ⚠ **`β` implícita NO mide el núcleo, y en el §1 de la v3.3 no medía nada** (v4.1 §3.3). Es la
  identidad `β = 0 ⟺ γ = 2−2H` del ruido gaussiano fraccionario: un contraste de coherencia entre
  dos exponentes. Ajustando `γ` y `H` **en la misma ventana**, `β` recorre 0.40 (v33) y 0.49
  (estacional) contra una distancia entre hipótesis de 0.21 y 0.16. El **+0.148** que la v3.3
  publicó es un punto arbitrario de esa curva.
- **El impacto «transitorio» de la v3.1 §2 era el REBOTE BID-ASK** (sesión 2026-08-10). Sobre el
  precio de transacción `D = 0.107`; sobre el punto medio, que no tiene rebote, `D = 11.53` y el
  signo se invierte. El spread mediano de la captura es 1 tick exacto.
- ⚠ **`D = 11.53` NO es «núcleo creciente medido»** (v4.1 §3.4, sesión 2026-08-12). Se cita como
  **«ajuste no identificado en el régimen `τ₀` en cota»** y nada más. `β = −0.160` con `τ₀` pegada
  a su cota inferior es firma de **mala especificación del estimador**: `β` se reporta en la
  literatura en `(0, 1)`, y un núcleo que crece sin cota sobre `[0, K]` implica impacto de mercado
  creciente indefinidamente, que es económicamente imposible (arbitraje ilimitado). El propio
  §5.2 de la v3.2 ya decía que `β` y `τ₀` no están identificados por separado. **No se abre un
  quinto estadístico para arreglarlo: se degrada la afirmación.**
- **La difusividad de la sesión (e) NO replica.** Sobre `captura_v33` la pendiente de la firma en
  ticks es **+0.091** (`H_p = 0.591`) con control barajado plano (−0.012), contra el +0.007 que
  dio `captura_larga`. **Este tramo es super-difusivo.**
- **Los «tres regímenes» de difusividad NO existen** (commit `432f459`). El rango de ajuste
  estaba fijo en ticks y la banda en segundos difería por factor 7; a banda común la reversión
  desaparece, y sobre banda común la pendiente **no es estimable** (cambia de +0.02 a +1.52 solo
  con la densidad de la rejilla). La difusividad más allá de ~12 min **sigue sin verificarse**.
- **φ′ se retira del vocabulario; se escribe `q̄`** (commit `9b2267e`). Es `1/q̄` y el hallazgo es
  Jones-Kaul-Lipson (1994), no un descubrimiento. Ni la estimación de `δ` por esa vía ni el apoyo
  a M1′ sobreviven a sus nulos.
- **Toda cifra de ejecución anterior al commit `0e3b9e0` se descarta**: el llenado adverso estaba
  clasificado como no-llenado y el sesgo ocultaba justo los llenados malos.

### Historial

Ocho tandas de trabajo aplicadas, en este orden:

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
6. **Convergencia de la cadena de medición de ω_m v2.1** — `ORDEN_TRABAJO_OMEGA_2_1.md`,
   **PARCIAL**: §2 (prioridad 1) y §4 completos, más las mediciones del §3. Detalle en
   "Sesión 2026-08-07". Quedan §5, §6 y §7.
7. **¿Existe ω_m? v2.2** — `ORDEN_TRABAJO_EXISTE_OMEGA_2_2.md`, experimentos §2, §3 y §5.
   Detalle en "Sesión 2026-08-07 (b)". **Ningún cambio al modelo.** Veredicto: **(A) dentro de
   la banda medida** — `ω_m` es salida del algoritmo, no del mercado — y **(C) fuera de ella**,
   pendiente de la captura de 48 h. **53/53** criterios.

8. **Oscilador forzado v3.0** — `ORDEN_TRABAJO_OSCILADOR_3_0.md`, §2.3 y §3 (el que decide).
   Detalle en "Sesión 2026-08-08". **Ningún cambio de modelo.** Veredicto: **k = 0, m < 0,
   raíces reales — NO HAY OSCILADOR**, y todo el AR(2) queda explicado a 7 decimales por paseo
   aleatorio + rebote bid-ask. **56/56** criterios.

9. **Cambio de herramienta: SSA en lugar de EMD, y φ′ contra el precio** — sesión
   2026-08-08 (f). **Ningún cambio a `Micelio.py`.** La búsqueda de `ω_m` por atractor de
   frecuencias queda **abandonada por decisión del operador**. Dos resultados: la
   descomposición SSA del precio real es **la escalera armónica de la ventana** (`T = 2L/k`,
   error 0.0000, idéntica sobre un paseo aleatorio), y **φ′ = ticks/BTC se asocia
   positivamente con la volatilidad realizada en 16 de 16 ventanas** (p de signos 3.05e-05).
   **El período y φ′ son independientes** (`L` explica el 79 % de la varianza de `T`; φ′ no
   añade nada significativo): se apuntaba a la **frecuencia** y lo que resultó medible es la
   **amplitud**.

⚠ **Lee las sesiones 2026-08-07 (b) y 2026-08-08 antes de tocar nada que dependa de `ω_m`.**
`ω_m`, `Φ`, `Ψ`, `Ω`, `Ω_crit`, `A_arm` y los nodos de fase **no tienen sustento empírico**.
Lo que sobrevive está en "Qué sobreviviría si (A) se confirma".

⚠ **El feed emite transacciones con `p = 0`** (~0.2 %). Corregido en `mercado.tick_valido`, pero
**toda medición de varianza anterior al 2026-08-08 está contaminada** — ver el hallazgo de datos
en la sesión 2026-08-08.

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
  (`python test_contaminacion_emd.py --capturar=440`), con `--barrido` para el §2.1 de la v2.1.
- `captura_dual.py` / `analisis_v21.py` — **v2.1**: captura simultánea de `@trade` y
  `@bookTicker` y las mediciones de §3, §4.2 y las condiciones del §8.
- `captura_larga.py` — **v2.2**: captura continua por bloques (48 h) con escritura atómica.
- `experimento_v22.py` — **v2.2**: nulos por sustitutos, histograma de banda y árbitro
  multitaper. `python experimento_v22.py`.
- `ssa.py` — **v3.2**: Análisis Espectral Singular. Descomposición exacta, w-correlación,
  barrido de `L` con parada por **mínimo local**, detección de pares, **Monte Carlo SSA**
  contra nulos AR(1) / ARIMA(1,1,0) / barajado, escalera de ventana y color de ruido.
  `python ssa.py --autotest` → **11/11**. **No se importa desde `Micelio.py`.**
- `barrido_ssa.py` / `graficar_ssa.py` — **v3.2**: toma de datos por ventanas (logea `t`,
  precio, id de tick, cantidad, signo y volumen neto, más autovalores, autovectores y
  w-correlación de **cada** `L`) y las figuras de matplotlib sobre ese log.
- `phi_precio.py` — **v3.2**: φ′ = ticks por volumen inyectado contra el precio, en bloques
  disjuntos, con nulo por desplazamiento circular y por barajado.
- `horizonte.py` — **v4.1**: `H_p` en los dos relojes con control barajado, `σ₁` en pb, y las
  dos curvas del §1 (`R²` medido contra `R²` requerido). `--autotest` → **14/14**.
- `tick_grande.py` — **v3.3 + v4.1 §3**: `γ` como exponente (no `C(1)`), `H`, `β` implícita, `η̂`
  de Robert-Rosenbaum, costes en pb y `N_eff` bajo memoria larga. La v4.1 añade los dos
  estimadores **espectrales** de `γ` (`gph`, `whittle_local`), el barrido de banda, `N_eff`
  coherente y `curva_beta_implicita`. `--autotest` → **28/28**. **No se importa desde `Micelio.py`.**
- `correcciones_v41.py` — **v4.1 §3**: el ejecutor de las correcciones a la v3.3, sobre las **dos**
  series (entrenamiento de `captura_v33` y el tramo continuo de 23.33 h de `captura_estacional`).
  `PREDICCION_SESGO_GAMMA_4_1.md` lleva las cuatro predicciones congeladas antes de medir.
- `experimento_v32.py` — **v3.2**: **el ejecutor del preregistro**. Etapas `muestra`, `fuga`,
  `delta`, `A`, `osc`, `B`, `decision`; abre el conjunto de prueba sólo en la última. `log()`
  transcribe a ASCII por sí sola, para que la consola cp1252 deje de ser una regla que recordar.
- `migracion_v32.py` — **v3.2**: el estimador M0/M1/M2/M2-osc y sus 18 controles.
- `oscilador.py` / `experimento_v30.py` — **v3.0**: primitivas `k`, `m`, `γ`, `Q` por AR(2) con
  hipótesis nula, verificación dimensional y descomposición del rebote bid-ask.
  `python experimento_v30.py`. **No se importa desde `Micelio.py`.**

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

## Sesión 2026-08-07 — Convergencia de la cadena de medición de ω_m (v2.1, §2 y §4)

Ejecuta `ORDEN_TRABAJO_OMEGA_2_1.md`. **Parcial**: §2 (prioridad 1) y §4 completos, más
las mediciones del §3. Quedan §5, §6 y §7. Suite: **51/51**.

### ⚠ LO QUE HAY QUE LEER PRIMERO: ω_m no tiene todavía anchura de mercado

Medido sobre 91 931 transacciones y 606 082 actualizaciones de libro en 899 s (ν = 102 tx/s):

| observable | σ_total² | σ_instr² | σ_genuina² | fracción instrumental |
|---|---|---|---|---|
| precio de transacción | 3.94e-05 | 7.70e-05 | **−3.76e-05** | **195 %** |
| `mid` de bookTicker | 1.88e-04 | 1.16e-04 | +7.14e-05 | 62 % |

Con el precio de transacción, **σ_genuina² sale NEGATIVA**. El §4.2 anticipa exactamente este
caso: "es un resultado con contenido, no un fallo — significa que toda la anchura observada es
del aparato". Dicho de otro modo: **la dispersión de `ω_m` entre rejillas es MAYOR que su
dispersión en el tiempo.** Cambiar la rejilla mueve `ω_m` más que cambiar de tramo de mercado.

Eso justifica por sí solo que el A/B siga congelado, y lo justifica con un número.

⚠ **Con la reserva del tamaño de muestra.** Son 8 ventanas temporales contra 9 rejillas: una
razón de varianzas con esos grados de libertad tiene un intervalo de confianza enorme. Lo
robusto es el signo y el orden de magnitud (instrumental ≥ total), no el 195 % exacto.

### ⚠ §3 REFUTADO POR MEDICIÓN: `mid` es MÁS pegajoso, no menos

El §3 propone pasar la cadena EMD al `mid` de `@bookTicker` por tener "resolución de medio
tick" y actualizarse "cuando se mueve cualquiera de los dos lados". La segunda mitad es cierta
para la tasa de MENSAJES y falsa para el movimiento real:

| | valor |
|---|---|
| actualizaciones de libro | 674.2 msg/s (**6.6×** la tasa de transacciones) |
| updates que **mueven** el mid | **1.1 %** → solo **7.7 cambios/s** |
| transacciones | 102.3 tx/s |
| duplicados a igual K — precio de transacción | p10 2.5 % · **mediana 9.1 %** · p90 15.8 % |
| duplicados a igual K — `mid` | p10 18.0 % · **mediana 23.8 %** · p90 29.0 % |

**`mid` da 2.6× MÁS duplicados que el precio de transacción.** La explicación es la
microestructura: con `tickSize = 0.10` el spread suele ser de un tick exacto, así que el mid se
queda clavado mientras las transacciones alternan entre bid y ask — el mismo rebote bid-ask que
la v2.0 detectó como ρ₁ negativo hace que el precio de transacción se mueva MÁS.

El propio §3.4 fija el criterio: "Si `mid` no baja de forma material, el observable no era el
problema y hay que volver sobre §1.3 con otra hipótesis." **Así que NO se cambia el observable.**
La cadena EMD sigue sobre el precio de transacción.

`§3.3`: el **49.9 %** de los cambios de `mid` no tiene transacción detrás — justo en la frontera
del "si domina, reconsiderar".

Dato que obliga a releer la v2.0: el 44.1 % de duplicados que se midió entonces era con
ν = 27 tx/s y K = 14. Con ν = 102 y K = 51 el mismo observable da **9.1 %**. La fracción de
duplicados **no es una propiedad del observable, es una propiedad de ν** — otra confirmación
del §1.4.

### §2 — La guarda de banda, y cuánto tiempo `ω_m` no es medible

`hht.omega_en_banda(periodo, dtau, W)` con la banda `[m_min·Δτ, W·Δτ/ciclos_min]`.

**Test de regresión histórica** (el que cierra el defecto de tres sesiones), con las tres vías
medidas en la v2.0 y W = 384:

| vía | Δτ | banda | período medido | veredicto |
|---|---|---|---|---|
| A escalera | 0.500 s | [2.5, 64.0] s | 95.2 s | **RECHAZA** |
| B pared fresca | 0.500 s | [2.5, 64.0] s | 262.6 s | **RECHAZA** (0.73 ciclos en ventana) |
| C ticks | 0.519 s | [2.6, 66.4] s | 13.2 s | acepta |

Y el dato que duele: sobre mercado real, `omega_valida` sale cierta en **12.5 %** de las
ventanas con precio de transacción y **25.0 %** con `mid`. La condición 3 del §8 pide > 90 %.
O sea que **tres cuartas partes del tiempo `ω_m` no es una medida**, y hasta ahora se usaba
igual.

La invalidez **alcanza a Φ, Ψ y Ω**, que es lo que el §2.2 marca como lo que más importa: con
`omega_valida = 0` se congela `Ω` en su último valor válido en vez de recalcularlo con un
artefacto. Test que lo fuerza: con los datos de prueba, `Ω` habría saltado de 1.5 a 0.13 y el
costo `κΩ²` del NMPC con ella. Más `edad_omega`, degradación a `Ω = 0` pasados
`T_OMEGA_RANCIA = 120 s`, y `w_ang` a **NaN**, nunca a cero.

⚠ **Lo que la guarda NO hace, y conviene tenerlo escrito:** comprueba RESOLUBILIDAD, no
EXISTENCIA. Una rampa con ruido produce una IMF dominante de ~20 s que cae dentro de banda y se
acepta, aunque no sea un ciclo de mercado. Detectar "no hay ciclo" es el papel de `C` y, según
el §4.4, de `σ_ω/ω`.

**El barrido del §2.1 no discrimina**, y el reporte lo dice en vez de fingir. Las 25
combinaciones que rechazan A y B aceptando C empatan todas en ventanas aceptadas: con 3
ventanas independientes no hay resolución para separarlas. Se conservan los valores del
documento (3.0 y 5.0) por estar dentro de la región válida, en vez de elegir un extremo del
empate y llamarlo evidencia.

### §4 — `Q_ω`, y una corrección numérica encontrada al testear

`Q_ω = (J·x)(J·x)ᵀ·σ_ω²` con `J = ∂A_arm/∂ω`, implementada sobre la desviación `x − x_ref`
coherente con la forma afín. Verificada contra derivada numérica a **4.5e-10**, simétrica,
semidefinida positiva, escalando con σ², y **exactamente cero con ω = 0** — la comprobación de
coherencia del §4.3: sin ciclo, no hay incertidumbre de ciclo.

⚠ **El jacobiano necesita su propio umbral de Taylor, y el valor estaba mal por razonarlo en
vez de medirlo.** `(ω·cos ω − sin ω)/ω²` resta dos cantidades de orden ω para dar una de orden
ω³: cancelación catastrófica. La rama de Taylor tiene el error opuesto. Medido:

| ω | error relativo de la fórmula exacta | error relativo de Taylor |
|---|---|---|
| 1e-6 | 7.8e-05 | 1.0e-13 |
| 1e-4 | 1.1e-08 | 1.0e-09 |
| **3e-4** | **~1e-08** | **~3e-08** ← cruce |
| 1e-3 | 5.6e-11 | 1.0e-07 |
| 1e-2 | 8.5e-13 | 1.0e-05 |

`UMBRAL_TAYLOR_JACOBIANO = 3e-4`. Un primer intento puso 1e-3 "porque ahí se cruzan", sin
medir, y dejaba una discrepancia del 0.2 %. Con 3e-4 queda en **3.0e-08**.

⚠ **Esto NO es un límite de convergencia del modelo armónico.** Ni `sin(ω)/ω` ni
`(ω·cos ω − sin ω)/ω²` divergen: tienen singularidades REMOVIBLES con límite finito (Δt y 0).
El umbral es una conmutación NUMÉRICA entre dos formas de evaluar la misma función, elegida
donde se cruzan sus curvas de error en float64. El único límite de validez real es la guarda
de banda del §2, y depende de la rejilla, no de la dinámica.

### Tres tests míos que estaban mal

Se anotan porque el patrón se repite: un test que falla por su propia culpa cuesta lo mismo
que uno real.

1. Le pedía a la guarda de banda que detectara **ausencia** de ciclo, cuando solo promete
   resolubilidad.
2. Medía el salto entre ramas de Taylor evaluando en **dos ω distintos**, así que medía la
   pendiente propia de la función (0.2 %) y no la discrepancia entre ramas (3e-8).
3. El sobrecosto del EAKF sombra fallaba por **carga de máquina** con una captura corriendo en
   paralelo. Corregido a **mínimo de varios bloques**, que es el estimador correcto de un coste
   intrínseco: el planificador solo puede añadir tiempo, nunca quitarlo. Medido 0.267 ms contra
   0.506 del peor bloque.

### Pendiente de la v2.1

- **§5 (log-precio), §6 (cuantiles empíricos) y §7 (péndulo)** sin ejecutar.
- **Ninguna de las cuatro condiciones del §8 se cumple**, así que el A/B sigue congelado:
  fracción instrumental 195 % / 62 % contra < 50 %; `omega_valida` 12.5 % / 25 % contra > 90 %;
  duplicados de `mid` peores, no mejores; §5.2 en régimen agitado sin repetir.
- El criterio `σ_ω/ω` del §4.4 en paralelo a `C` está **sin implementar**.
- `ciclos_min` y `muestras_por_ciclo_min` siguen `[CALIBRAR]`: hace falta una captura mucho
  más larga para que el barrido discrimine.

## Sesión 2026-08-07 (b) — ¿Existe ω_m? (v2.2, experimentos §2, §3 y §5)

Ejecuta `ORDEN_TRABAJO_EXISTE_OMEGA_2_2.md`. **Ningún cambio al modelo**, como exige su §8:
todo es código de análisis aparte. Suite **53/53**. Captura de 48 h lanzada y corriendo.

### ⚠ VEREDICTO: (A) dentro de la banda medida, (C) fuera de ella

**Dentro de ~2.5 s – 300 s no hay escala característica.** `ω_m` es salida del algoritmo, no
del mercado. Tres líneas de evidencia independientes, todas en la misma dirección:

| evidencia | resultado |
|---|---|
| §3 REAL contra BARAJADO, precio de transacción | KS **p = 0.98** → indistinguibles |
| §3 REAL contra BARAJADO, `mid` | KS p = 0.28 → indistinguibles |
| §3 ídem en la otra captura (ν = 27 contra 102) | KS p = 0.60 → indistinguibles |
| §3 medianas de período, REAL / BARAJADO | **112.0 s / 112.7 s** |
| §5 multitaper contra AR(1), precio | 6.6 % de frecuencias sobre el umbral 95 % (azar ≈ 5 %) |
| §5 multitaper contra AR(1), `mid` | 7.4 % |
| §5 pico multitaper / modo EMD | **184×** de separación |
| §2 rechazos de banda | **100 % por arriba**, 0 % por abajo |
| §2 `T/(W·Δτ)` de los rechazados | concentrado en **~0.60** (dispersión 24 % en `mid`) |
| §4.2 fracción instrumental (precio) | **195 %** |

**Barajar los incrementos destruye todo el orden temporal y `ω_m` no se entera.** Y el nulo es
un nulo de verdad: curtosis 601.8 contra 601.8, masa en cero 65.6 % contra 65.6 % — la marginal
se conserva exacta, así que la única diferencia posible era la estructura temporal, y no la hay.

### La demostración que lo cierra: la EMD no tiene hipótesis nula

Fijado como test permanente (`test_v22_emd_no_tiene_hipotesis_nula`). Sobre un **paseo
aleatorio puro** —sin escala característica por construcción— la cadena devuelve:

```
periodo 118.1 s con C = 0.690      <- ventana completa (384 muestras)
periodo  37.5 s                    <- media ventana
```

El período **escala con la ventana**, no con la señal. Eso es exactamente lo que el §1.2
anticipa citando la línea Flandrin/Rilling/Gonçalvès y Wu/Huang: sobre ruido fraccionario la
EMD se comporta como un banco diádico y el modo dominante lo fija la ventana.

Y explica de golpe tres sesiones de síntomas: los `C` de 0.94–1.00 de la v1.3, el factor 7 del
§5.2 de la v2.0, la `σ_genuina² < 0` de la v2.1. No eran defectos distintos: **eran la misma
cosa, medida desde tres sitios.**

### Lo que NO queda decidido, y por qué importa

La hipótesis **(B)** —ciclo en decenas de minutos— **no queda probada ni refutada**. El
multitaper alcanza ~300 s y la EMD ~65 s con la configuración actual. Ver un ciclo de 30 min con
`Δτ = 0.5 s` exige `W ≥ 10 800` contra los 384 actuales: factor 28×.

Por eso está corriendo la **captura de 48 h** (`captura_larga.py`, bloques cada 5 min con
escritura atómica a temporal y `rename`, para no repetir el defecto de volcado de la v2.0). El
§4 multiescala se ejecuta sobre ella, y hasta entonces **no se decide sobre el diseño**.

### Reservas honestas de esta medición

- **Potencia estadística baja.** 4 ventanas independientes. Un KS que no rechaza significa "no
  se detectó diferencia", NO "son iguales". Lo que sostiene la lectura no es el p-valor sino que
  las medianas coinciden y que las tres líneas convergen.
- **El brazo IAAFT no es fiable aquí.** Da 11.6 s en un observable y 118 s en el otro: aplicado
  al NIVEL de precio (paseo aleatorio) no converge. La lectura **no se apoya en él**, se apoya en
  BARAJADO. Con eso, separar (A) de (A′) —"no hay ciclo" contra "hay memoria, no ciclo"— queda
  pendiente.
- **Un solo par de capturas**, de 899 s y 439 s.

### TABLA DE MEDICIONES — todos los números crudos de la v2.2

Captura principal: **91 931 transacciones en 899 s** (ν = 102.3 tx/s, K = 51, W = 384,
Δτ ≈ 0.51 s). Captura secundaria: 11 870 transacciones en 439 s (ν = 27.0, K = 14).
Reproducible con `python experimento_v22.py`.

**§2 — Por qué lado falla la banda** (banda [2.5, 65.0] s; ventana W·Δτ = 195 s)

| observable | períodos p10 / MED / p90 / max | dentro | **rechazo por arriba** | por abajo | ventanas indep. |
|---|---|---|---|---|---|
| precio de transacción | 52.2 / **112.0** / 263.0 / 305.1 s | 2/8 (25.0 %) | **6/8 (75.0 %)** | 0 (0.0 %) | 4 |
| `mid` | 36.9 / **104.3** / 153.9 / 187.9 s | 3/8 (37.5 %) | **5/8 (62.5 %)** | 0 (0.0 %) | 4 |

`T/(W·Δτ)` de los rechazados por arriba:

| observable | p10 | MED | p90 | dispersión |
|---|---|---|---|---|
| precio de transacción | 0.490 | **0.596** | 1.412 | 50.9 % |
| `mid` | 0.525 | **0.604** | 0.865 | **24.4 %** ← concentrado |

**§3 — Nulo por sustitutos, precio de transacción**

| vía | f_hz p10 / MED / p90 | C MED | `omega_valida` | T MED | σ_instr²/σ_total² |
|---|---|---|---|---|---|
| **REAL** | 0.00384 / **0.00893** / 0.01921 | 0.702 | 12.5 % | **112.0 s** | **195.3 %** |
| IAAFT | 0.04939 / 0.08711 / 0.15974 | 0.654 | 87.5 % | 11.6 s | 92.8 % |
| **BARAJADO** | 0.00535 / **0.00888** / 0.01223 | 0.783 | 0.0 % | **112.7 s** | 482.8 % |

**§3 — Nulo por sustitutos, `mid`**

| vía | f_hz p10 / MED / p90 | C MED | `omega_valida` | T MED | σ_instr²/σ_total² |
|---|---|---|---|---|---|
| **REAL** | 0.00662 / **0.00967** / 0.03004 | 0.698 | 25.0 % | 104.3 s | **62.0 %** |
| IAAFT | 0.00664 / 0.00847 / 0.01944 | 0.741 | 25.0 % | 118.2 s | 122.1 % |
| **BARAJADO** | 0.00418 / 0.00638 / 0.01274 | 0.815 | 0.0 % | 158.6 s | 143.0 % |

**§3 — Kolmogorov-Smirnov sobre la distribución de f_hz**

| comparación | p | lectura |
|---|---|---|
| REAL vs BARAJADO, precio | **0.9801** | indistinguibles |
| REAL vs BARAJADO, `mid` | 0.2827 | indistinguibles |
| REAL vs BARAJADO, captura ν = 27 | 0.6000 | indistinguibles |
| REAL vs IAAFT, precio | 0.0025 | distintas ⚠ ver reserva sobre IAAFT |
| REAL vs IAAFT, `mid` | 0.9801 | indistinguibles ⚠ inconsistente con la fila anterior |

**§3 — El nulo ES un nulo** (marginal de incrementos, REAL contra BARAJADO)

| observable | curtosis | masa en cero |
|---|---|---|
| precio de transacción | 601.8 → **601.8** | 65.6 % → **65.6 %** |
| `mid` | 165.7 → **165.7** | 95.9 % → **95.9 %** |
| captura ν = 27 | 216.8 → **216.8** | 74.3 % → **74.3 %** |

**§5 — Multitaper (Thomson) contra nulo AR(1)**, NW = 4, 7 tapers, bandas al 95 %

| observable | n log-retornos | Δτ | AR(1) `a` | frecuencias sobre umbral | pico | exceso | modo EMD | separación |
|---|---|---|---|---|---|---|---|---|
| precio de transacción | 1800 | 0.171 s | 0.4816 | 59/900 = **6.6 %** | 0.6 s | ×1.78 | 112.0 s | **184×** |
| `mid` | 1802 | 0.171 s | 0.3125 | 67/901 = **7.4 %** | 0.5 s | ×1.64 | 103.4 s | **190×** |

Por azar se esperaría ~5 % de frecuencias por encima del umbral del 95 %. **Ninguno de los dos
observables muestra exceso significativo sobre ruido rojo.** El multitaper cubre períodos de
0.34 s a ~308 s, o sea que **sí incluye** el modo de 112 s que la EMD devuelve — y ahí no ve nada.

**El test que lo cierra** (`test_v22_emd_no_tiene_hipotesis_nula`, sintético y determinista)

| serie | ventana | período que la EMD devuelve | C |
|---|---|---|---|
| paseo aleatorio puro | 384 muestras | **118.1 s** | 0.690 |
| el mismo | 192 muestras | **37.5 s** | — |

**§4.2 — Fracción instrumental por observable** (8 ventanas temporales, 9 rejillas)

| observable | σ_total² | σ_instr² | σ_genuina² | fracción |
|---|---|---|---|---|
| precio de transacción | 3.94e-05 | 7.70e-05 | **−3.76e-05** | **195 %** |
| `mid` | 1.88e-04 | 1.16e-04 | +7.14e-05 | 62 % |

### Qué sobreviviría si (A) se confirma a todas las escalas

Del §6.2, y conviene tenerlo escrito para que el veredicto no se lea como catástrofe: sobreviven
**toda la v2.0** (reloj de ticks, ingesta por lotes, `Q(Δt)`, deduplicación, huecos), **toda la
capa de riesgo de la v1.3**, **Loeper y el NMPC** (necesitan `λ` y `Γ`, no `ω`) y con alta
probabilidad **`R_n`** — estimar una tendencia es mucho más robusto que extraer un ciclo.

Habría que rediseñar `A_arm`, el término `γ_ω·ω_m` de `ρ_k`, `c²_vol`, los nodos de fase y la
Sec. 1.4 entera (`Φ`, `Ψ`, `Ω`, `Ω_crit`). Eso sería una v3.0 con su propio documento.

El §6.3 cataloga cuatro alternativas sin implementar ninguna. La cuarta —**ciclo estocástico de
Harvey dentro del propio EAKF**— merece atención cuando llegue el momento: es `A_arm` con
amortiguamiento y `λ` estimada CONJUNTAMENTE con el estado, así que **tiene nulo** (si `ρ → 0`
el ciclo no existe), `σ_ω` sale nativa de la covarianza de parámetros, y desaparece el traspaso
externo `ω_m → A_arm` donde viven la trampa del 2π, la guarda de banda y el problema de rejilla.

### Estado de la captura de 48 h y qué falta para cerrar la v2.2

Lanzada con `captura_larga.py --horas=48`, bloques cada 300 s. Primer bloque verificado:
22 623 transacciones a 75.4 tx/s.

⚠ **Defecto propio corregido durante el lanzamiento.** `np.savez_compressed` **añade `.npz`** al
nombre si no lo lleva, así que escribiendo en `bloque_00000.npz.tmp` acababa creando
`bloque_00000.npz.tmp.npz` y el `os.replace` fallaba con `WinError 2`. La primera captura corrió
**media hora sin guardar un solo bloque**. Se arregla pasando un descriptor de archivo abierto en
vez de la ruta. Es exactamente el modo de fallo que el §4.4 pedía evitar, y casi se repite.

Pendiente para cerrar la v2.2:

- **§4 multiescala** (C0–C4, Δτ de 0.5 s a 120 s, cinco configuraciones × tres observables) sobre
  la captura larga. Es el estadístico `T/(W·Δτ)` por configuración lo que decide entre (A) y (B).
- **§4.3**: reevaluar el observable por configuración, ahora con el criterio corregido
  —`σ_instr²/σ_total²` y validez, no duplicados— e incluyendo el precio medio ponderado por
  volumen, que a escalas gruesas es el observable natural.
- Repetir el §5.2 de la v2.0 en **régimen agitado** (ν > 300 tx/s), condición 4 del §8 de la v2.1.
- Un IAAFT fiable —sobre incrementos y no sobre el nivel— para poder separar (A) de (A′).

Y en cola, **no bloqueadas** por el veredicto porque no dependen de `ω_m` (§7 de la v2.2):
**v2.1 §5** (log-precio) y **v2.1 §6** (cuantiles empíricos del NIS). Las compuertas del NIS
están mal calibradas pase lo que pase con `ω_m`, así que son el mejor candidato para la espera.

## Sesión 2026-08-08 — Oscilador forzado: el veredicto (v3.0, §2.3 y §3)

Ejecuta `ORDEN_TRABAJO_OSCILADOR_3_0.md`. **Ningún cambio de modelo en `Micelio.py`**, como
exige su §8. Suite **56/56**.

### ⚠ VEREDICTO: k = 0. NO HAY OSCILADOR. Y ahora se sabe qué había en su lugar

El AR(2) sobre **446 892 ticks reales** (3.19 h continuas, tras descartar los ceros del feed):

| magnitud | valor |
|---|---|
| φ₁ | **+0.783939** ± 0.000087 |
| φ₂ | **+0.216061** ± 0.000089 |
| `k = 1 − φ₁ − φ₂` | **+4.67e-07** ± 1.7e-04 |
| `m = −φ₂` | **−0.216061** ← **NEGATIVA** |
| raíces | **REALES** (sin oscilación) |
| `Q` | **indefinida** (`k·m < 0`) |
| **`H₀: k = 0`** | **p = 0.6175 → NO se rechaza** |

Y el ajuste robusto de Huber, que por el §3.3 prevalece si difiere: `k = 6.7e-16`, tres órdenes
de magnitud aún más cerca de cero.

### La explicación mecánica, exacta a 7 decimales

`m < 0` no es "poca inercia". Es la firma de que el segundo rezago viene de la
**autocorrelación de los RETORNOS**, no de una inercia del precio. Si

```
x_t = x_{t-1} + r_t        con    r_t = a·r_{t-1} + ε
```

entonces, sin aproximar nada:

```
x_t = (1+a)·x_{t-1} − a·x_{t-2} + ε      ⟹   φ₁ = 1+a,  φ₂ = −a
                                          ⟹   φ₁ + φ₂ = 1  IDÉNTICAMENTE
                                          ⟹   k = 0  POR CONSTRUCCIÓN,  m = a
```

Medido sobre los datos reales:

| | valor |
|---|---|
| `a = ρ₁(retornos)` — el rebote bid-ask | **−0.216061** |
| φ₁ predicho `1+a` / medido | +0.783939 / **+0.783939** (error **2.0e-07**) |
| φ₂ predicho `−a` / medido | +0.216061 / **+0.216061** (error **2.7e-07**) |

**Todo el AR(2) del mercado real es un paseo aleatorio más rebote bid-ask.** No queda residuo
que atribuir a una fuerza recuperadora. Y encaja con lo ya medido: la v2.0 reportó ρ₁ = −0.164
en la innovación del filtro y lo identificó como rebote bid-ask; es la misma cantidad.

### Segunda vía independiente, y gratis: Harvey tampoco existe

El §4.1 da la correspondencia `ρ = √(−φ₂)`, `λ = arccos(φ₁/(2√(−φ₂)))`. Con **φ₂ = +0.216 > 0**,
`√(−φ₂)` no es real: **el ciclo estocástico de Harvey no tiene parametrización válida sobre
estos datos**. No hace falta ajustar nada por máxima verosimilitud para saberlo — el §4 queda
respondido por el mismo número que respondió el §3.

### ⚠ La guarda que faltaba: el §3.3 pide vigilar `γ < 0` y el fallo real fue `m < 0`

`γ` salió positiva (+1.216) en el ajuste global y en **0 de 38** bloques fue negativa. La guarda
que el documento especifica no habría disparado nunca. La que hacía falta es sobre la masa, y
está añadida (`oscilador.guarda_masa`) con su test.

### ⚠ HALLAZGO DE DATOS QUE OBLIGA A RELEER TRES SESIONES

**El feed emite transacciones con `p = 0` y `q = 0`**, con id de trade válido y monótono. En las
tres rutas de captura:

| captura | ceros | fracción |
|---|---|---|
| `captura_trades` (stream simple `/ws/`) | 29 / 11 870 | **0.244 %** |
| `captura_dual` (stream combinado) | 76 / 91 931 | 0.083 % |
| `captura_larga` (combinado) | 990 / 467 326 | 0.212 % |

**En producción no se filtraban.** `publicar_trade` habría escrito `P_spot = 0` y el EAKF habría
tomado una medición de 0 USD/BTC: a ~30 tx/s, una innovación espuria de −65 000 USD **cada
~30 s**. Corregido en `mercado.tick_valido`, en el borde de la ingesta.

Y la corrección de las cifras del registro es contraintuitiva:

| sobre `captura_dual` | con ceros | sin ceros |
|---|---|---|
| curtosis de incrementos | **601.8** | **1179.7** |
| \|incremento\| máximo | 65 245.5 USD | **11.8 USD** |

⚠ **La curtosis SUBE al limpiar.** Los ceros inflaban tanto la varianza que el denominador σ⁴
aplastaba el cociente: estaban **enmascarando** la cola real, no creándola. Toda cifra de
varianza medida sobre la serie sucia —incluida la curtosis 601.8 que la v2.1 §1.5 y la v2.2
citan— está mal, y el factor 5 500 en el incremento máximo dice cuánto.

### Estabilidad por bloques (§3.3)

38 bloques de 5 min sobre el tramo continuo:

| magnitud | p10 | mediana | p90 | dispersión |
|---|---|---|---|---|
| `k` | −6.44e-05 | +1.66e-05 | +1.88e-04 | **250 %** |
| `Q` | 0.0018 | 0.0035 | 0.0103 | (finito en 21/38) |

`k` cambia de signo entre bloques y su dispersión es del 250 %. El §3.3 pregunta explícitamente
si `k` y `Q` varían tanto como variaba `ω_m` de ventana en ventana: **sí**. Con la diferencia de
que ahora sabemos por qué — están fluctuando alrededor de cero.

### Lo que esto cierra y lo que no

**Cierra** la pregunta del §1.2 de la v2.2 por una vía con hipótesis nula correcta, sobre todos
los datos, sin ventana espectral y sin banda de resolubilidad. El resultado coincide con el de
la v2.2 (`REAL ≈ BARAJADO`) y con el multitaper, que no veía exceso sobre ruido rojo.

**No cierra** la hipótesis (B) —escala de decenas de minutos—: 3.19 h de tramo continuo dan
poca resolución ahí, y el AR(2) a Δn = 1 tick es un instrumento de escala fina. La captura
larga sufrió un **corte de DNS de 36 442 s** (10 h) que partió los datos en dos tramos; hay que
relanzarla para el §4 multiescala de la v2.2.

⚠ **Reserva sobre el alcance del AR(2):** se ajusta sobre el NIVEL con intercepto constante, así
que mide reversión a un nivel FIJO. Una reversión a una tendencia móvil aparecería como `k = 0`
igualmente. Los bloques de 5 min mitigan esto —ahí el nivel es localmente constante— y también
dan `k ≈ 0`, pero no lo eliminan del todo. El modelo estructural de Harvey con tendencia +
ciclo lo resolvería, y es justo el que no tiene parametrización válida.

## Sesión 2026-08-08 (b) — Horizonte derivado (v3.1, §1 completo)

Ejecuta `ORDEN_TRABAJO_PROPAGADOR_3_1.md`. **Preregistro commiteado ANTES de medir** en
`66bed94` (`PREREGISTRO_3_1.md`), como exige su §4. `Micelio.py` sin cambios.

### §1 — El horizonte, derivado por primera vez

Sobre el tramo continuo limpio de 446 892 ticks (3.19 h, ν = 39 tx/s):

| magnitud | valor |
|---|---|
| ρ₁ de retornos (rebote bid-ask) | −0.216061 |
| σ_r sobre datos **limpios** | 0.2218 USD |
| **`s_eff` de Roll** | **0.2062 USD/BTC** ≈ 2 ticks |

Coste de ida y vuelta y horizonte al que la volatilidad lo iguala:

| esquema | `c(u)` | en ticks | **`H*`** |
|---|---|---|---|
| maker + maker | 25.97 USD/BTC | 260 | **50.0 s** |
| maker + taker | 45.66 USD/BTC | 457 | ~200 s |
| taker + taker | 65.35 USD/BTC | 653 | **fuera del rango medido** |

**Firma de volatilidad** (`σ(H)` contra `H`), el gráfico que el §8 exige y que dice
empíricamente a qué escala deja de dominar la microestructura:

| H [s] | 0.5 | 1 | 2 | 5 | 10 | 30 | 60 | 120 | 300 |
|---|---|---|---|---|---|---|---|---|---|
| σ [USD] | 2.22 | 3.34 | 5.11 | 8.30 | 11.91 | 20.58 | 28.24 | 35.47 | 51.66 |
| σ/√H | 3.15 | 3.34 | 3.61 | 3.71 | **3.77** | **3.76** | 3.65 | 3.24 | 2.98 |

La meseta está en **H = 10–30 s**: por debajo, σ/√H cae; por encima, también. Ese es el rango
donde la difusión es más limpia, y coincide en orden de magnitud con el `H*` derivado.

**Dos consecuencias que reordenan el proyecto**, y ahora con números propios:

1. Con comisiones taker, `σ` **no alcanza el coste dentro de los 300 s medidos**. Al régimen de
   volatilidad de esta captura, una ida y vuelta taker no se paga. El §1.3 lo anticipa: si `H*`
   se dispara, **el sistema debe operar menos, no más**.
2. **La latencia de 300 ms deja de ser una restricción.** Contra `H* = 50 s` es el **0.6 %**.
   Lleva seis versiones condicionando decisiones de diseño.

⚠ **Criterio NO CUMPLIDO, declarado en el preregistro:** el escalón de comisiones **no se pudo
leer de la cuenta** — `/fapi/v1/commissionRate` es firmado y el Modo LECTURA no tiene
credenciales. Se usan las tarifas públicas VIP 0, marcadas como asumidas en el código.

### §2 — Bloqueado por una omisión mía, ya corregida

El §2 necesita `ε` del campo `m` de Binance. **Mis capturas nunca lo persistieron**: el bot lo
lleva en `MERCADO_DTYPE` desde la v2.0, pero `captura_dual.py` y `captura_larga.py` guardaban
precio, tiempo, id y cantidad, no `m`. Sin él no hay forzamiento medido y el §2 es inejecutable.

Corregido: `captura_larga.py` guarda ahora `tr_maker` y además **filtra los ceros del feed en
origen**. Captura nueva lanzada; primer bloque verificado (42.5 % de trades con `m = True`).

### ⚠ §3 — La calibración por simulación NO reproduce la tabla del documento

El §3.2 mide un sesgo de tamaño finito de 0.03–0.10 entre `β` teórico `(1−γ)/2` y el `β` al que
`ρ₁(retornos)` cruza cero, **decreciente** con γ. Mi simulación da sesgos de **0.40–0.55** y
**crecientes**. Dos defectos propios encontrados y corregidos por el camino:

1. El núcleo `h(k) = G(k) − G(k−1)` se construía desde τ = 1, así que `h(1)` salía exactamente
   cero — el propagador no decaía en el primer rezago, que es donde está casi todo el efecto.
2. La convolución era `mode="same"`, que **centra el núcleo y mete el futuro en el presente**.
   Sobre un test de predictibilidad eso fabrica la señal que se pretende medir. Es la misma
   familia que la media móvil centrada contra la que advierte el §5.
3. El generador de signos no producía la γ pedida (0.30 → 0.542 medido). Exponente del núcleo
   corregido a `α = (1+γ)/2`, y la calibración pasa a evaluarse contra la **γ medida**.

Aun así el sesgo no converge al del documento. Sospecha principal: mi `G(τ) = (1+τ)^(−β)` tiene
**impacto permanente `G(∞) = 0`**, mientras que el propagador de Bouchaud tiene `G(∞) > 0` — con
todo el impacto transitorio, el precio revierte por construcción y el cruce se desplaza a `β`
altas. **No se ajusta el simulador hasta que cuadre**: eso sería ajustar el instrumento al
resultado esperado. Se deja como no reproducido y el §3 no se ejecuta sobre datos reales hasta
resolverlo.

## Sesión 2026-08-08 (c) — Adenda A: la pendiente de la firma. El hallazgo NO sobrevive

Ejecuta `ADENDA_A_PROPAGADOR_3_1.md`, que sustituye el §3 entero. **Simulador eliminado** como
exige el A.5. `Micelio.py` sin cambios.

### ⚠ EL CONTROL OBLIGATORIO DEL A.3.4 FALLA — y por eso existe

Sobre incrementos **barajados**, donde por construcción no hay estructura temporal alguna, la
pendiente **no sale cero**:

| serie | pendiente [30, 180] s | pendiente [0.5, 180] s |
|---|---|---|
| REAL | −0.1530 | −0.1194 |
| **BARAJADA (control)** | **−0.0403** | **−0.0796** |

El A.3.4 es explícito: "Si sale distinta [de cero], el estimador tiene sesgo propio y el
resultado sobre datos reales no vale." **El estimador tiene sesgo negativo propio.** La
comparación honesta es −0.153 contra −0.040, no contra cero.

### Y el bootstrap tampoco separa

`IC 95 % = [−0.1522, +0.0699]`, mediana −0.0530, **p contra cero = 0.4167**. Bloques móviles de
900 s, 12 bloques efectivos — por debajo de los ~15 que el A.3.3 pide como mínimo.

### Los dos estimadores de la firma discrepan EN SIGNO

| rango | solapada | no solapada |
|---|---|---|
| [0.5, 180] s | −0.1194 (reversión) | **+0.0111 (momentum)** |
| [0.5, 30] s | −0.1235 (reversión) | **+0.0421 (momentum)** |
| [30, 180] s | −0.1530 (reversión) | −0.0580 (reversión) |

El A.3.2 pide reportar el acuerdo de signo entre rangos; aquí falla algo peor, el acuerdo entre
**estimadores**. La ponderación es la causa: la firma solapada pone un punto de partida por
tick, así que **sobrepondera los tramos de alta actividad**, que son también los de alta
volatilidad. Solapada y no solapada no estiman la misma cantidad.

Por eso `H_lo` sale distinto según el estimador: la firma solapada decrece monótonamente desde
0.5 s —no hay región de microestructura creciente— mientras la no solapada sí la tiene y pone la
meseta en 10–30 s, como decía el A.2.

### Veredicto sobre el hallazgo preliminar del A.2

La razón de varianzas **sí** es monótona (`VR(45)=0.936`, `VR(60)=0.845`, `VR(90)=0.785`,
`VR(120)=0.703`, `VR(180)=0.567`), que es lo que hacía atractivo el hallazgo. Pero **el control
barajado reproduce el mismo comportamiento cualitativo**, el bootstrap da p = 0.42, y los dos
estimadores discrepan en signo.

**No se puede rechazar difusividad con estos datos.** La pista del A.2 no sobrevive a su propio
control — que es exactamente para lo que el A.3.4 lo puso, y la segunda vez que este control
salva un análisis (la primera fue `C = 0.783` sobre barajados en la v2.2).

⚠ Esto **no cierra** el criterio de abandono del §4.2 del preregistro: ese exige además que
`E[Δp]` no supere `c(u)` en ningún régimen y que el residuo sea blanco, y ninguna de las dos se
ha medido todavía. Falta también el §2 (`β`, `γ`, forma de `G(τ)`), ahora ya posible porque la
captura persiste `tr_maker`.

### Lo que hace falta antes de volver a mirar la pendiente

1. **Resolver la discrepancia entre estimadores.** Mientras solapada y no solapada den signos
   opuestos, ninguna pendiente es reportable. La ponderación por actividad es la sospecha.
2. **Más datos.** 12 bloques de bootstrap contra los ~15 mínimos; y el A.3.3 dimensiona el
   procedimiento para 48 h, no para 3.19 h.
3. **Corregir el sesgo del estimador** o contrastar siempre contra el barajado en vez de contra
   cero.

## Sesión 2026-08-08 (d) — §2: el propagador existe; el sobrepaso NO está establecido

Primera ejecución del §2 con **forzamiento medido**, ahora que la captura persiste `tr_maker`.
Tramo continuo de 19 157 ticks (0.95 h, ν = 5.6 tx/s — tramo tranquilo).

### El test de signo del §2.2: NO falla, y por poco lo doy por fallado

`G(0) = 0` **exactamente**. Estuve a punto de aplicar la parada obligatoria, pero el cero es
**artefacto de mi estimador**: uso `R(τ) = E[(p_{t+τ} − p_t)·ε_t]`, que en `τ = 0` es
idénticamente nulo por construcción. No es un signo invertido, es una definición. El impacto
inmediato en sentido de propagador es `R(1)`.

Y `R(1) = −0.0052 < 0`, que tampoco es signo invertido: **es el rebote bid-ask**. Una compra
agresora ejecuta al ask y el trade siguiente tiende al bid. Es la misma cantidad que el
`ρ₁ = −0.216` de la v3.0.

**La convención `m = True → ε = −1` es CORRECTA**, y el criterio honesto lo demuestra:

| | valor |
|---|---|
| `R(τ) > 0` para τ ≥ 2 | **799 de 799 rezagos (100 %)** |
| con `ε` invertido | espejo exacto, negativa en todo el rango |

Si el signo estuviera invertido, `R` sería negativa en TODO el rango, no solo en τ = 1.

### El propagador existe, y con margen enorme

| τ (ticks) | 1 | 2 | 10 | 50 | 100 | 200 | 300 | 500 | 800 |
|---|---|---|---|---|---|---|---|---|---|
| `R(τ)` [USD/BTC] | −0.005 | +0.005 | +0.088 | +0.510 | +0.943 | +1.525 | +1.839 | +1.834 | +1.621 |

**Control con signos barajados** (precios intactos, relación destruida): el pico se desploma de
**1.9488 a 0.0028–0.0155**, o sea un factor **~150–700**. La relación entre signo de transacción
y movimiento posterior del precio es real y masiva. **El propagador con forzamiento medido
existe**, que es lo que el §2 quería establecer.

### ⚠ Pero el sobrepaso NO está establecido, y el fallo es de mi estadístico

`R` alcanza su máximo en τ = 398 ticks (**71.1 s**, cerca del `H* = 50 s`) y cae a 1.621 en
τ = 800 — un descenso del 17 %. Eso *parece* sobrepaso y reversión, que es la firma de raíces
complejas, o sea de oscilador.

**No se puede afirmar.** Elegí como estadístico la razón `R(final)/R(pico)`, y bajo el nulo esa
razón es **degenerada**: cuando el pico es ruido (~0.003), la razón explota a valores absurdos
(p5 = −560, mediana = −11). El `p = 0.917` que devuelve no significa nada. Es un estadístico mal
elegido, no un resultado.

El control B (incrementos barajados, signos intactos) da mediana 0.722 y p95 = 1.000, con el
real en 0.832: **dentro del rango del nulo**. Con 24 ventanas independientes a τ = 800, un
descenso del 17 % es perfectamente compatible con ruido de muestra.

**Conclusión honesta:** el propagador existe y es grande; su forma —monótona o con sobrepaso—
**queda sin decidir**, y hace falta (a) un estadístico con nulo bien definido y (b) mucho más
dato. La captura de 8 h con `tr_maker` está corriendo.

### Pendiente del §2

- Ajuste del núcleo paramétrico de 4 parámetros (`G₀`, `τ₀`, `β`, `δ`) con bootstrap por bloques.
- `δ` barrida y elegida por verosimilitud **fuera de muestra**.
- `γ` de la autocorrelación de signos, y la comprobación de coherencia
  `pendiente ≈ (1−γ)/2 − β` del A.4, que es la validación del mecanismo por dos vías.
- Raíces de `G(τ)` con un estadístico cuyo nulo no sea degenerado.

## Sesión 2026-08-08 (e) — La firma en tiempo de ticks resuelve la discrepancia

Tres tareas pedidas: firma en reloj de ticks, control positivo, `τ_pico` en ambos relojes.

### ⚠ CONFIRMADO: el sesgo era mezclar dos relojes, no el agrupamiento de volatilidad

La confirmación estaba en mis propios datos y no la señalé: el barajado destruye el
**agrupamiento** de volatilidad pero conserva los **tiempos de llegada**. Con incrementos iid,
los tramos con más ticks siguen teniendo más varianza por segundo, mecánicamente. Que la
pendiente barajada saliera −0.0403 y no cero era exactamente esa firma.

**Firma en tiempo de ticks** (`σ(n)/√n` contra `n` en ticks), 446 892 ticks, ν = 39 tx/s:

| n [ticks] | 2 | 8 | 32 | 128 | 512 | 1024 | 4096 | 8192 |
|---|---|---|---|---|---|---|---|---|
| solapada | 0.1964 | 0.2088 | 0.3045 | 0.4579 | 0.5631 | 0.5839 | 0.5721 | 0.5393 |
| **no solapada** | **0.1984** | **0.2127** | **0.3044** | **0.4584** | **0.5640** | **0.5835** | 0.5466 | 0.5474 |

**Los dos estimadores coinciden a tres decimales.** La discrepancia de signo desaparece **por
construcción**, no por corrección de sesgo.

Y el control negativo se limpia solo:

| pendiente sobre barajados | reloj de pared | **reloj de ticks** |
|---|---|---|
| rango completo | −0.0796 | **−0.0054** |
| rango largo | −0.0403 | **−0.0028** |

**El sesgo del estimador desaparece.** Es la confirmación de que venía de los dos relojes.

**Resultado en tiempo de ticks:**

| rango | pendiente | lectura |
|---|---|---|
| [2, 256] ticks | +0.228 | microestructura: el rebote bid-ask muriéndose |
| **[256, 8192] ticks** | **+0.007 / +0.003** | **DIFUSIVO** (control barajado: −0.003) |

Más allá de la microestructura, **el precio es difusivo en tiempo de ticks**. La "reversión" del
A.2 era íntegramente el artefacto de los dos relojes.

Consecuencia aceptada: `H*` pasa a `H*_ticks = (c/σ_tick)²`, y su valor en segundos depende de
ν. Con ν variando por factor 20, **`H*` en segundos no es una constante** — más honesto, no menos.

### ⚠ El control positivo funciona, pero su nulo NO es conservador

Sustituyendo la razón degenerada por un bootstrap paramétrico bajo núcleo monótono
`G(τ) = G∞·τ/(τ₀+τ)`, con σ calibrada a la volatilidad real (0.1806 USD):

| | valor |
|---|---|
| descenso medido pico→final | **16.8 %** |
| descenso espurio bajo el nulo monótono | MED 0.0 %, p90 0.1 %, **máx 0.4 %** |
| sorteos ≥ 16.8 % | **0 de 60** |

Parecía cerrado. **No lo está**, por dos razones que encontré después:

1. **`R` no tiene un solo pico.** Extendiendo el rango: 1.95 en τ=398, baja a 1.55 en 1000,
   **sube a 2.09 en 1800**, vuelve a bajar. El "descenso" depende de dónde se trunca:

   | max_rezago | 600 | 800 | 1000 | 1500 | 2000 | 2300 |
   |---|---|---|---|---|---|---|
   | pico | 398 | 398 | 398 | 1500 (borde) | 1642 | 1642 |
   | descenso | 9.9 % | 16.8 % | 20.4 % | 0.0 % | 12.1 % | 16.9 % |

2. **Mi nulo usa signos iid.** `rng.permutation(eps_real)` destruye la memoria larga del flujo
   de órdenes — que es precisamente lo que produce estas ondulaciones a rezagos largos. El nulo
   es **demasiado estrecho**, y por eso da máx 0.4 % donde el real da 16.8 %.

**Conclusión:** el sobrepaso sigue sin establecerse, y ahora se sabe qué haría falta — un nulo
monótono con **flujo de memoria larga**, no con signos iid. Es la misma lección otra vez: el
control positivo solo vale si su nulo reproduce las propiedades del dato que importan.

⚠ Nota de calibración: con σ ajustada a la volatilidad real el nulo da máx 0.4 %; con otra SNR
da 15.4 %. **El resultado depende críticamente de la SNR**, así que el número aislado no
significa nada sin declarar cómo se calibró.

### `τ_pico` en ambos relojes — sin conclusión

| captura | ν [tx/s] | pico [ticks] | pico [s] |
|---|---|---|---|
| v31 | 5.6 | 398 (o 1642 según ventana) | 71 (o 293) |
| larga (regla de tick) | 39.0 | 668 | 17.1 |

**No es estable en ninguno de los dos relojes**, pero el propio `τ_pico` no es una cantidad bien
definida mientras `R` tenga varios máximos locales. La pregunta —¿en qué reloj vive el
decaimiento del impacto?— sigue abierta y es la que sostiene la arquitectura de dos relojes
desde la v2.0. Hace falta más dato antes de responderla.

### Lo que sí queda establecido de esta sesión

- La firma de volatilidad **debe medirse en tiempo de ticks**; en reloj de pared el estimador
  tiene sesgo propio y los dos muestreos discrepan en signo.
- Más allá de la microestructura (n > 256 ticks) el precio es **difusivo**.
- El propagador **existe** (sesión d): pico 1.95 contra 0.003–0.016 con signos barajados.

## Sesión 2026-08-08 (f) — SSA en lugar de EMD, y φ′ contra el precio (v3.2)

Cambio de estrategia pedido por el operador: **se abandona la búsqueda de `ω_m` como atractor
de frecuencias** y se pasa a SSA (Análisis Espectral Singular). Objetivo declarado: una
**ventana de toma de datos** con logeo completo, no un modelo. **`Micelio.py` sin cambios.**

### Por qué SSA y no EMD: la EMD no tiene hipótesis nula, SSA sí

La v2.2 dejó demostrado que la EMD devuelve un "ciclo" de 118.1 s sobre un paseo aleatorio
puro, y de 37.5 s al partir la ventana en dos. SSA **no arregla eso por sí solo** — se
comprobó y también fabrica pares oscilatorios sobre ruido. Lo que sí tiene es el test que
falta: **Monte Carlo SSA** (Allen & Smith 1996), que contrasta los autovalores observados
contra los de un nulo ajustado a los mismos datos.

### La condición de parada del barrido de `L`

Buscar `ortogonalidad == 0` no termina nunca: el ruido de medición —el mismo que alimenta `R`
en el EAKF— siempre filtra energía entre componentes, así que el mínimo alcanzable es
estrictamente positivo y desconocido. Implementado como **rejilla finita** (no hay bucle) con
elección por **mínimo local interior**, y dos salvaguardas:

- si el mínimo cae en un extremo de la rejilla, devuelve `hay_minimo_local = False` en lugar
  de entregar el borde como si fuera una elección;
- cada punto del barrido lleva al lado **el mismo estadístico medido sobre sustitutos
  barajados**. Una métrica que baja igual sin estructura no mide separación, mide la rejilla.

### Controles del propio estimador (`python ssa.py --autotest`) — 11/11

Este proyecto llevaba cinco sesiones ejecutando controles **negativos** con disciplina. Estos
son la mitad que faltaba, los **positivos**: dada una verdad conocida, ¿el estimador la
recupera?

| control | resultado |
|---|---|
| reconstrucción exacta (suma de elementales = serie) | error **1.1e-12** |
| dos períodos conocidos, 120 y 37 muestras | **0.63 %** y **0.02 %** de error |
| ruido blanco → β | **+0.002 ± 0.017** |
| paseo aleatorio → β | **+1.977 ± 0.018** (ROJO) |
| señal enterrada en AR(1), MC-SSA | la marca, **p = 0.000** |
| AR(1) puro, MC-SSA | 6/20 sobre p95 |
| curva monótona → mínimo local | **no discrimina** (correcto) |
| paseo aleatorio → escalera `2L/k` | error **0.0317** (cae en ella) |
| señal con ciclos → escalera `2L/k` | error **0.6020** (la rompe) |

⚠ **Sobre un paseo aleatorio puro, SSA también fabrica pares oscilatorios**, con período que
escala con la ventana: 152.4 muestras con L=400 (0.38·L) y 60.4 con L=200 (0.30·L). Es el
defecto de la EMD intacto. Lo que cambia es que **MC-SSA lo detecta**: 1/20 componentes sobre
el p95, o sea el 5 % del azar.

### ⚠ HALLAZGO: los autovectores del precio real son los armónicos de la VENTANA

Medido sobre 8 ventanas de 8 192 ticks de `captura_v31b` (97 742 tx, 4.33 h, ν = 6.0 tx/s) y
8 de `captura_larga` (480 757 tx, 14.39 h, ν por ventana de **21 a 109 tx/s**).

`T/L` de los primeros autovectores, **idéntico en las 16 ventanas**, con `L` entre 96 y 768 y
en los dos observables (precio y volumen neto acumulado):

```
1.600   0.889   0.667   0.500   0.400   0.333   0.286
```

Que es `T = 2L/k`. Error mediano al escalón, y energía del primer autovector:

| serie | error a `2L/k` | E(EOF1) | E del par dominante |
|---|---|---|---|
| **8 ventanas reales, v31b** | **0.0000** (6/8 exacto) | 97.6–99.8 % | **0.00–0.01 %** |
| **8 ventanas reales, larga** | **0.0000–0.0335** | 84.2–99.4 % | 0.00–0.06 % |
| ctrl_paseo (paseo aleatorio) | 0.0347 | 90.8 % | 0.16 % |
| ctrl_barajado (incrementos reales barajados) | 0.0264 | 98.1 % | 0.03 % |
| **ctrl_positivo (T = 120 y 37 conocidos)** | **0.1250 / 0.4235** | **45.8 %** | **85.93 %** |

El control positivo recupera su verdad: 128 ticks (verdad 120) con el 85.9 % de la energía y
37 ticks (verdad 37) con el 11.8 %, y MC-SSA marca **exactamente las componentes 1, 2, 3 y 4**
—los dos pares verdaderos, ni una más—. Sobre BTC real, la descomposición es
**indistinguible de la de un paseo aleatorio**, y los pares que el detector encuentra llevan
0.00–0.06 % de la energía.

**La escalera no se mueve al cambiar ν por un factor 18** (6 → 109 tx/s). Si el período
viniera del mercado, cambiar la tasa de transacciones lo movería.

### Ortogonalidad por ventana

| serie | media \|w-corr\| | razón real/barajado |
|---|---|---|
| 8 ventanas reales v31b, precio | 0.1425–0.1456 | **0.882–0.993** |
| 8 ventanas reales larga, precio | 0.1392–0.1532 | 0.944–1.043 |
| ctrl_paseo | 0.1459 | 0.998 |
| ctrl_barajado | 0.1448 | 0.987 |
| **ctrl_positivo** | **0.0933** | **0.629** |

La métrica discrimina —el control positivo se separa con claridad— y sobre datos reales dice
que la separación es la de un sustituto barajado. La curva contra `L` es además **plana** en
las ventanas reales: los mínimos locales que se eligen son poco profundos, así que **la `L`
elegida sobre datos reales no es una elección informada**, y así queda anotado.

### Sobre retornos la escalera se rompe, pero no aparece ciclo

Aplicar SSA al **nivel** garantiza que el primer EOF sea la tendencia y se lleve ~98 %. Se
repitió sobre `retorno` y `vol_neto` (v31b, 8 ventanas):

| observable | escalera | E(EOF1) | E del par dominante | razón real/nulo |
|---|---|---|---|---|
| retorno | 0.125–6.16 (med **1.16**) | 0.3–11.5 % | 0.65–3.00 % | 0.63–0.83 |
| vol_neto | 0.007–0.125 (med **0.032**) | 1.0–5.6 % | 0.56–5.54 % | 0.81–1.16 |
| ctrl_positivo (retorno) | 5.53 | 29.7 % | **58.58 %** | 0.862 |
| ctrl_paseo (retorno) | 11.28 | 0.28 % | 0.56 % | 0.750 |

En retornos ya no hay escalera, pero tampoco energía oscilatoria: 0.65–3.00 % contra 58.58 %
del control. Y `ctrl_paseo` da razón 0.750, dentro del rango de los datos reales, así que esa
métrica **no discrimina** en este observable.

### ⚠ Dos límites de MC-SSA encontrados por los controles

1. **Sobre series casi blancas el test es inservible.** Con `nulo=ar1` sobre retornos marca
   **23–30 de 30** componentes en los datos reales… y **30/30 en `ctrl_paseo`**, que es ruido
   iid puro. El sesgo anticonservador conocido (proyectar sobre las EOF de los propios datos)
   se vuelve fatal cuando el espectro es plano y todos los autovalores son parecidos.
2. **`ar1_incrementos` es demasiado severo con señales fuertes.** En `ctrl_positivo` sobre el
   nivel da **0/30**: pierde un ciclo verdadero y grande, porque el AR(1) ajustado a los
   incrementos absorbe la propia sinusoide y ensancha el nulo. `ar1` sobre el nivel lo detecta
   con 4/4.

**Ningún nulo es el bueno por sí solo**, y por eso el log guarda los tres. Es la misma lección
de la sesión (e): el control solo vale si su nulo reproduce las propiedades del dato que
importan.

### ⚠ φ′ CONTRA EL PRECIO — y la separación que el operador pidió

`φ′ = ticks por volumen inyectado [ticks/BTC]`, la derivada del reloj de transacciones
respecto del de volumen. Sobre bloques **disjuntos** de m = 64 ticks: `φ′_b = m / Σq`. Es el
inverso del tamaño medio de operación: φ′ alta = muchas operaciones pequeñas.

El criterio del operador: si `T/L` es constante en todas las ventanas pero la asociación
φ′–precio **varía**, entonces la escalera es del aparato y la asociación es del mercado.
Medido:

| magnitud (v31b) | rango entre ventanas | sd | signo | p de signos |
|---|---|---|---|---|
| **T/L (período SSA)** | **0.0908** | **0.0314** | — | — |
| **Spearman(log φ′, V)** | 0.2486 | 0.0942 | **8/8** | **0.0078** |
| Pearson(log φ′, log V) | 0.4679 | 0.1411 | 7/8 | 0.0703 |
| ρ(log φ′, \|ΔP\|) | 0.3441 | 0.1037 | 7/8 | 0.0703 |
| ρ(log φ′, flujo neto) | 0.7449 | 0.3111 | 4/8 | 1.0000 |
| ρ(log φ′, nivel de precio) ESPURIA | 0.4707 | 0.1406 | 2/8 | 0.2891 |
| ρ predictiva (φ′ → volatilidad siguiente) | 0.2237 | 0.0796 | 5/8 | 0.7266 |
| ρ predictiva (φ′ → retorno siguiente) | 0.2695 | 0.0830 | 2/8 | 0.2891 |

**El período es la magnitud más estable de todas** (rango 0.09) y la asociación φ′–volatilidad
varía 3× más. Son dos cosas distintas, que es exactamente lo que había que decidir.

**Y la asociación es real.** Replicada en las dos capturas, con nulo por **desplazamiento
circular** del volumen contra el precio (conserva intacto el agrupamiento de operaciones
grandes, que barajar destruye; los precios no se tocan, así que la volatilidad del bloque es
idéntica bajo el nulo y lo único que se rompe es el emparejamiento):

| captura | ν | ρ Spearman mediana | p < 0.05 individual | V(φ′ alto)/V(φ′ bajo) |
|---|---|---|---|---|
| v31b | 4.6–9.6 tx/s | +0.204 (0.033 a 0.282) | 4/8 | **1.181** (1.07–1.25) |
| larga | 21–109 tx/s | **+0.338** (0.114 a 0.484) | **6/8**, z hasta **+5.27** | **1.412** (1.12–1.87) |

**16 de 16 ventanas con signo positivo → p de signos = 3.05e-05.** En terciles: la volatilidad
realizada es un **18 % mayor** (v31b) y un **41 % mayor** (larga) en el tercil alto de φ′ que
en el bajo, y la razón supera 1 en las 16 ventanas.

Lecturas, en orden:

1. **Cuando el mercado se fragmenta en operaciones pequeñas, la volatilidad realizada por
   unidad de volumen sube.** Es una propiedad del mercado, no del estimador: el nulo no la
   reproduce y el efecto es más fuerte donde ν es mayor.
2. ⚠ **Pearson no era ciego, estaba roto — y el defecto era mío.** Hay bloques con volatilidad
   **exactamente cero** (64 ticks sin que el precio se mueva): **2.82 %** en `captura_v31b` y
   **0.62 %** en `captura_larga`. Metidos en un logaritmo se van a −27.6, unas **30
   desviaciones** fuera, y dominaban la correlación ellos solos. Con ellos dentro, Pearson daba
   −0.093 a +0.158 en v31b y −0.558 a +0.391 en larga, errático y sin significación.
   Excluyéndolos —y solo de las correlaciones sobre `log V`— pasa a **7/8 y 8/8 positivo**, con
   medianas **+0.193** y **+0.363**, de acuerdo con Spearman. Spearman no los necesita: el
   rango de un cero está bien definido y es el más bajo, que es justo lo que son.
   **Sin la comprobación gráfica no se habría visto**: los ceros aparecieron como una fila de
   puntos aplastada contra el borde del eje logarítmico.
3. **No hay componente direccional ni predictiva.** ρ(φ′, retorno) da 2/8 y 3/8; ρ(φ′, retorno
   siguiente) da 2/8 y 3/8; ρ(φ′, volatilidad del bloque **siguiente**) da 5/8 y 3/8. φ′ dice
   **cuánto** se mueve el precio en el mismo bloque, no **hacia dónde**, y no lo anticipa.
4. La correlación con el **nivel** de precio es espuria por construcción (dos series no
   estacionarias) y sale 2/8: se calcula porque se preguntó por ella, y se reporta marcada.

### ⚠ El período y φ′ son INDEPENDIENTES — se apuntaba a la frecuencia y lo medible era la amplitud

Medido cruzando las dos magnitudes sobre las **16 ventanas** (8 de `captura_v31b` y 8 de
`captura_larga`), que es la comprobación que faltaba para cerrar la separación:

| par | Pearson | Spearman | t (n = 16) |
|---|---|---|---|
| **`T` contra `L`** (la ventana de análisis) | **+0.891** | **+0.938** | **+7.36** |
| `T/L` contra φ′ mediana | −0.347 | −0.162 | −1.38 |
| `T` contra φ′ mediana | −0.469 | −0.306 | −1.99 |
| `T/L` contra ν | −0.105 | −0.165 | — |
| `T/L` contra ρ(φ′, V) | −0.153 | −0.103 | — |
| φ′ mediana contra ν | −0.082 | −0.035 | — |

Con n = 16 el umbral es |t| > 2.14. **El período lo determina `L` y nada más**: `L` explica el
**79 %** de la varianza de `T`, y la dependencia residual de φ′ —una vez quitada `L`— no es
significativa. No se puede descartar una dependencia débil con 16 ventanas; sí se puede
descartar que sea el mecanismo dominante, porque apenas queda varianza por repartir.

**Lo que esto significa para el diseño del sistema.** Desde la v1.1 el proyecto apunta a una
**frecuencia**: `ω_m`, `ω_ang`, `A_arm`, los nodos de fase, `Φ`, `Ψ`, `Ω`, `Ω_crit`, `C`,
`c²_vol`, el término `γ_ω·ω_m` de `ρ_k`. Todas esas cantidades son una frecuencia o una fase.
La v3.0 ya había dictado `k = 0, raíces reales, NO HAY OSCILADOR`; esta sesión añade que el
período que devuelve la descomposición es el de la ventana y no el del mercado.

**Lo único que sobrevive a la medición es una relación de amplitud**: φ′ contra la volatilidad
realizada, 16/16 ventanas, p de signos 3.05e-05, con un efecto del +18 % al +41 % entre
terciles. φ′ no mueve el período: mueve **cuánto** se mueve el precio.

**Consecuencia, y queda marcada como hipótesis de diseño SIN PROBAR:** una variable que
modula amplitud y no fase no entra en la matriz de transición `A` —que es donde vive la
dinámica— sino en las **covarianzas de ruido**. En términos del EAKF, φ′ es candidata natural
a entrar en `Q_k` y `R_k`, que es exactamente donde la Sec. 7.3.3 del PDF pone `ρ_k`. La forma
`ρ_k = 1 + γ_ω|ω_m| + γ_Q|ΣQ|` tiene el término equivocado: `γ_ω|ω_m|` es la frecuencia, que
no tiene sustento empírico, mientras que un término en φ′ sí lo tendría. **Nada de esto está
medido todavía** —haría falta comprobar que modular `Q` con φ′ mejora el NIS y la blancura de
la innovación— y no se ha tocado `Micelio.py`.

### Reservas de esta sesión

- **La `L` elegida sobre datos reales no significa nada.** La curva de ortogonalidad contra
  `L` es plana; los mínimos son de baja prominencia. El barrido discrimina en el control
  positivo y no en el mercado.
- **Todo esto mide DENTRO de la ventana**: 8 192 ticks son ~22 min a ν = 6 y ~2 min a ν = 109,
  y `L ≤ 1024`. **No dice nada sobre la hipótesis (B) de la v2.2** (escala de decenas de
  minutos), que sigue sin decidir.
- El efecto de φ′ es **modesto y variable**: individualmente significativo en 10 de 16
  ventanas. Lo que lo sostiene es la consistencia de signo, no la magnitud de ninguna.
- `captura_larga` no persiste `tr_maker`, así que en esa réplica el flujo firmado no es
  interpretable. `captura_v31b` sí lo trae.
- m = 64 ticks por bloque es una elección **sin calibrar**. Habría que barrer m y comprobar a
  qué escala vive la asociación — que es la misma pregunta de los dos relojes de la sesión (e).

### Pendiente

- ⚠ **φ′ NO va en `Q_k`, y la propuesta de arriba está mal planteada.** Se corrige aquí en vez
  de borrarla, porque el error es instructivo. Ver la nota siguiente.
- Barrer `m` en φ′ y ver si la asociación vive en el reloj de ticks o en el de volumen.

### ⚠ Corrección sobre φ′ (2026-08-09): qué es, qué no es, y dónde va

**Qué es.** `φ′ = m/Σq` con `m = 64` **fijo**, así que `corr(log φ′, −log volumen) = +1.000000`
exactamente: **φ′ es función solo del volumen**, y sus unidades son transacciones/BTC = 1/BTC.
**No** es precio/volumen, así que no es dimensionalmente la λ de Kyle ni `G₀`. En este proyecto
«tick» significa *transacción* desde la v2.0, y confundirlo con el tick de precio lleva justo a
esa lectura.

El hallazgo se reformula sin misterio: **a número de transacciones fijo, menos volumen va con más
volatilidad.**

**No es artefacto de construcción.** Test disjunto —φ′ de los ticks **pares**, σ de los
**impares**, sin nada compartido—: mediana **+0.1883**, **14/16** ventanas, p de signos
**0.0042**. Sobrevive.

**Pero sí es un proxy de `G₀`, y eso importa más.** Estimando `G₀` por bloque con la regresión
`Δp = G₀·ε·q`:

| | mediana | signo | valores |
|---|---|---|---|
| **REAL** `ρ(log φ′, G₀)` | **+0.6323** | **8/8** | 0.618 0.600 0.691 0.602 0.734 0.549 0.670 0.646 |
| **NULO** (desplazamiento circular de `q, ε` contra precios) | −0.0908 | 2/8 | −0.183 … +0.164 |

El nulo era obligatorio porque `G₀ = Σ(Δp·εq)/Σ(εq)²` y `φ′ ∝ 1/Σq` **comparten `q`**. No lo
reproduce: la asociación es real.

### ⚠ Segunda corrección, del mismo día: el enlace con `G₀` NO sobrevive

**1. Retirar φ′ del vocabulario. Escribir `q̄` (tamaño medio de operación).** El
`corr = +1.000000` no es un hallazgo: es una identidad algebraica, porque `log φ′ = log 64 −
log Σq`. φ′ no es una variable, es un **nombre para `1/q̄`**, y el nombre hacía daño — invitaba a
tratar como descubrimiento una reexpresión. Con `q̄` cada frase se lee sola: «φ′ alta se asocia a
volatilidad» → «**operaciones pequeñas se asocian a volatilidad**».

**2. Y así reformulado, el hallazgo tiene nombre en la literatura.** A igual número de
operaciones, operaciones más pequeñas van con más volatilidad: eso es **Jones, Kaul & Lipson
(1994)** — el número de transacciones, no su tamaño, es lo que porta la información de
volatilidad. Uno de los hechos estilizados mejor replicados de la microestructura. **No es un
descubrimiento.** Lo bueno es que la tubería mide cosas reales y el test disjunto lo confirma
limpiamente; lo útil es que ahora la literatura dice qué esperar.

**3. Mi nulo de φ′–`G₀` no contrastaba lo que yo creía.** Bajo desplazamiento circular el
numerador de `G₀ = Σ(Δp·εq)/Σ(εq)²` se vuelve una suma de signos aleatorios, así que lo que se
destruye es la **alineación precio-flujo** — no el acoplamiento por **escala de volumen**, que
sobrevive porque la magnitud sigue yendo como `1/√Σq²`. El nulo acreditaba algo cierto pero
distinto de lo que hacía falta.

**El test que sí discrimina, y falla.** Con `Δp = G·ε·q^δ` y `q ≈ q̄` dentro del bloque sale
`log G₀ = const + (1−δ)·log φ′`, así que la pendiente da `δ` por una ruta independiente. Medido:

```
pendiente = +1.6286   IC95 bootstrap [+1.4896, +1.7562]   ->   delta = -0.63
```

**`δ = −0.63` está fuera del rango admisible** (`δ ∈ [0,1]`): implicaría que operaciones más
grandes mueven **menos** el precio en términos absolutos. La derivación supone `q` poco disperso
dentro del bloque, y el **coeficiente de variación de `q` por bloque es 2.60 (p90 3.94)** — la
aproximación `q ≈ q̄` está gruesamente violada. **El caveat resultó ser la restricción
vinculante, no una nota al pie.** No hay `δ` que sellar como predicción.

**4. Y la predicción sobre M1′ tampoco se sostiene.** La correlación parcial controlando por el
tamaño medio:

| | ρ(σ, G₀ \| q̄) |
|---|---|
| REAL | +0.2563 (976 bloques) |
| **NULO — signo barajado dentro del bloque** | **+0.2625** (666 bloques) |

El nulo conserva `Δp` y `q` intactos —y con ellos todo acoplamiento por escala— y rompe **solo**
la alineación precio-flujo. **Lo reproduce exactamente**, así que la asociación es **mecánica**:
`σ` y `G₀` comparten `Δp`, y condicionar a `G₀ > 0` selecciona bloques donde el precio se movió
*con* el flujo. **M1′ no recibe apoyo por esta vía.**

Lo que sí queda medido y no es mecánico: `ρ(q̄, G₀) = −0.7022`, o sea que operaciones más grandes
tienen menor impacto **por unidad**, que es la concavidad del impacto (`δ < 1`). Pero cuantificar
`δ` desde aquí exige bloques con `q` mucho menos disperso.

**Dónde va en la arquitectura, y no es en `Q`.** Descomponiendo:

```
Δp  =  φ′·(flujo OBSERVADO)  +  φ′·(flujo NO observado)
```

El primer término es **medido** y pertenece a la ecuación de estado, con φ′ como coeficiente. Solo
el segundo va a `Q`. Meter φ′ entero en `Q` es devolver al residuo una cantidad que se puede
medir — **el mismo error que el §0 de la v3.1 diagnosticó en el AR(2)**, reaparecido en otro sitio.

Y si algo de φ′ acaba en `Q`, dos correcciones más:
1. **La dependencia es cuadrática**: `Var(Δp) = φ′²·Var(flujo)`. Un término `γ_φ|φ′|` es
   incoherente con su propia derivación.
2. **No es aditivo con lo que ya está.** `ρ_k = 1 + γ_ω|ω_m| + γ_Q|ΣQ|` trata flujo e impacto
   como contribuciones independientes y no lo son: φ′ es impacto **por** flujo, así que la
   cantidad física es el producto `(φ′·|ΣQ|)²`, no dos sumandos. Sustituir `γ_ω|ω_m|` por un
   tercer sumando conserva la incoherencia y solo cambia la etiqueta.
- Barrer `m` en φ′ y ver si la asociación vive en el reloj de ticks o en el de volumen.
- Nulo de MC-SSA que funcione sobre series casi blancas (el sesgo anticonservador lo invalida).
- SSA multivariante (M-SSA) sobre `[precio, φ′]`, que es el paso natural ahora que hay una
  asociación establecida entre los dos canales.
- Sigue en pie todo lo de la v3.1: núcleo paramétrico del propagador, `γ` de la
  autocorrelación de signos y la comprobación `pendiente ≈ (1−γ)/2 − β`.

## Sesión 2026-08-10 — v3.2 EJECUTADA: la compuerta pasó y la decisión se para en el paso 3

Ejecuta `PREREGISTRO_3_2.md` de punta a punta sobre `captura_v33`, con el preregistro
**congelado desde el 2026-08-09** y sus 9 enmiendas todas anteriores a mirar dato.
`Micelio.py` sin cambios (`git diff --stat` vacío). Suite **56/56**.
Ejecutor nuevo: `experimento_v32.py`.

### ⚠ VEREDICTO: el paso 2 pasa, el paso 3 falla. Hay señal medible y está 3.5× por debajo del coste

```
PASO 1  compuerta          1 031 154 ticks continuos contra 741 000     PASA
PASO 2  LL/N contra M0     M2 gana; IC95 [+6.857e-03, +8.084e-03]       PASA
PASO 3  criterio economico q90(|mu|) = 11.30  contra  1.5*c(u) = 39.05  FALLA
PASOS 4 y 5                no se ejecutan: la regla se para en el primero que falla
```

**No es «no hay señal».** M2 bate a M0 fuera de muestra con un IC que excluye 0 con holgura,
`R²` fuera de muestra +0.014, y el residuo es blanco (`|ρ₁| < 0.007` en los tres modelos). La
estructura está ahí y está medida. Lo que falla es que **no llega a pagar las comisiones**:

| | valor |
|---|---|
| `\|μ̂\|` de M2 en prueba | q50 **5.16** · q90 **11.30** · q99 15.39 · **máx 17.12** USD/BTC |
| umbral `1.5·c(u)` | **39.05** USD/BTC |
| `c(u)` que el paso 3a exigiría | 7.53 USD/BTC = **0.579 pb** por lado |
| comisión maker VIP 0 **asumida** | **2.000 pb** por lado → falta un **factor 3.46** |

`c(u)` es aquí **comisión pura**: maker+maker no cruza el spread, así que
`c = 2 × 0.0002 × 65 100 = 26.03 USD/BTC`. El obstáculo no es la microestructura, son 4 pb.

⚠ **`μ̂` de M1 es CERO IDÉNTICAMENTE, y no es un fallo del cálculo.** Con impacto permanente
(`f_∞ = 1`) el núcleo es constante, así que `μ̂_t(H) = G₀·Σ_j x_{t−j}·[G(H+j) − G(j)] ≡ 0`: el
movimiento ya ocurrió y no queda nada que capturar. Es el contenido económico de la hipótesis de
MkII, y por eso el paso 3 se evalúa sobre `μ̂` y no sobre el ajuste.

⚠ **El §11 (abandono) NO se activa**, y hay que decir por qué con precisión. El §9.2 lo
condiciona a «3a falla **con 3b satisfecha**». Aquí 3b también falla, pero **3b es degenerado en
este caso**: `frac(|μ̂| ≥ 1.5c) = 0` exactamente porque `máx|μ̂| = 17.12 < 39.05`, o sea que su
fallo es una *consecuencia* de 3a y no una medida de potencia. La potencia sí está acreditada por
otro lado —el paso 2 detecta la estructura con un IC que excluye 0—, así que el caso real no es
ninguno de los dos que el §9 contempla: **señal medible, por debajo del coste.**

### ⚠ El resultado que reordena la lectura: el signo de `D` se INVIERTE entre observables

El §3 del preregistro declara el punto medio como observable primario y el precio de transacción
como control secundario, con regla falsable escrita antes. Ajustando **la misma M2 sobre las dos
columnas**, con la misma muestra y la misma malla:

| observable | `β` | `τ₀` | `D = G(K)/G(0)` | lectura |
|---|---|---|---|---|
| **punto medio** (primario) | **−0.160** | pegada al límite | **11.53** | ~~núcleo creciente~~ → **ajuste no identificado** (v4.1 §3.4) |
| precio de transacción (control) | **+4.69** | pegada al límite | **0.107** | decae al 11 % |

**El impacto «transitorio» que la v3.1 §2 midió sobre el precio de transacción era el rebote
bid-ask.** El spread mediano de esta captura es **0.1000 USD = exactamente 1 tick** (p90 igual),
así que una compra imprime en el ask y la siguiente vuelve al bid: el 89 % de la «reversión» es
Roll, no un propagador. Sobre el punto medio, que no tiene rebote, **la reversión desaparece y el
signo se invierte**.

Es la misma familia de hallazgo que el §5.2 de la v2.0 y la Adenda A: una cantidad que parecía del
mercado y era del instrumento. Aquí el preregistro lo cazó porque obligaba a llevar **las dos
columnas siempre**.

### ⚠ `D > 1` no es «permanente»: la dicotomía del §5.2 no contiene el resultado

⚠ **RECLASIFICADO por la v4.1 §3.4 (2026-08-12).** Todo este apartado se leyó en su día como
«núcleo creciente medido». **No lo es.** `β < 0` con `τ₀` en cota es firma de mala especificación
del estimador, no del mercado, y la afirmación queda degradada a **«ajuste no identificado en el
régimen `τ₀` en cota»**. Lo que sigue se conserva como registro de lo que se midió y de por qué
la regla de una cola del §5.2 no contenía el caso — no como evidencia sobre la forma del núcleo.

`D = 11.53` con `β < 0` significa núcleo que **crece** con el rezago: ni permanente (`D = 1`) ni
transitorio (`D < 1`). El perfil lo confirma: `D(K/4) = 9.24 → D(K/2) = 10.32 → D(K) = 11.53`.

La regla escrita del §5.2 es **de una cola** —«`D̂ ≥ q05` → no se rechaza `D = 1` → permanente»— y
aplicada literalmente habría declarado **impacto permanente y MkII bien especificada** sobre datos
que dicen lo contrario. Se lee la distribución simulada **completa**, que el propio §5.2 manda
generar, y por arriba también:

```
D_hat = 11.5252     q05 = 0.0000     q95 = 7.2704     ->  se rechaza D = 1 POR ARRIBA
```

Esto **no es una enmienda**: es leer entera la misma simulación. Callarlo sería reportar
«permanente» sobre evidencia de lo contrario.

### ⚠ Y el contraste de `D` se derrumba a la γ real — el nulo emparejado lo delata

El §5.2.quater obliga a **regenerar `q05` con la `γ̂` medida**, y ahí está el hallazgo:

| | γ = 0.3 (§5.2.bis, antes de datos) | **γ̂ = 0.798 (medido)** |
|---|---|---|
| `D̂` mediana bajo `D = 1` verdadero | 0.9950 | **0.2854** |
| `q05` | 0.9698 | **0.0000** |
| `q95` | — | **7.2704** |
| sd del nulo | 0.0158 | enorme |

**Con la autocorrelación de signos real, el estimador de `D` devuelve cualquier cosa entre 0 y 7.3
cuando la verdad es exactamente 1.** El contraste permanente-contra-transitorio **no tiene
potencia a la γ de este mercado**, y eso lo dice el propio procedimiento que el preregistro
congeló. La malla de potencia del §5.2.bis se midió a γ = 0.3 y no transporta.

Dicho de otro modo: `D̂ = 11.53` cae fuera del nulo, pero el nulo es tan ancho que la cifra no
sostiene una lectura fina. Lo que sí sostiene es que **el signo del efecto no es el que la
dicotomía esperaba**.

⚠ Reserva: el nulo se regeneró con **24 sorteos**, no más, por coste (206 s). `q95` con 24
muestras es una estimación pobre.

### ⚠ Este tramo es SUPER-DIFUSIVO — la sesión 2026-08-08 (e) NO replica

Firma de volatilidad en **tiempo de ticks**, no solapada, sobre entrenamiento, **con su control de
incrementos barajados al lado** (que es lo que la sesión (e) estableció como obligatorio):

| n [ticks] | 2 | 8 | 32 | 64 | 256 | 1024 | 4096 | 16384 |
|---|---|---|---|---|---|---|---|---|
| σ/√n **real** | 0.799 | 0.520 | 0.343 | **0.327** | 0.433 | 0.505 | 0.588 | **0.634** |
| σ/√n **barajado** | 1.102 | 1.102 | 1.100 | 1.097 | 1.099 | 1.127 | 1.138 | 0.993 |

```
pendiente [256, 16384]   real +0.0907     barajado -0.0123
H_p                      real  0.591      barajado  0.488     (0.5 = difusivo)
```

**El control barajado sale plano y el real no.** La sesión 2026-08-08 (e) concluyó «más allá de la
microestructura el precio es difusivo» con pendiente +0.007 sobre `captura_larga` (ν = 39 tx/s);
sobre `captura_v33` (ν = 17.5 tx/s) la pendiente es **+0.091**, trece veces mayor y con el control
limpio. **La difusividad no es una propiedad estable del mercado**, o al menos no lo es entre estas
dos capturas. ~~Y encaja por dos vías con el núcleo creciente: un `G` que crece y un `H_p > 0.5`
son el mismo hecho medido desde dos sitios.~~ ⚠ **Esa última frase se retira** (v4.1, 2026-08-12):
el «núcleo creciente» está reclasificado como ajuste no identificado (§3.4) y `H_p = 0.591` no
sobrevivió al reloj de pared (§2, banda [0.371, 0.552]). Dos cifras retiradas no se corroboran
entre sí.

⚠ **No hay plateau, y eso obligó a resolver una ambigüedad del preregistro.** El §2.2 escribe
`H*_ticks = (c/σ_tick)²` sin decir **cuál** `σ_tick`. Esa fórmula es la solución de `σ(H) = c`
*bajo difusión exacta*, y aquí la curva es en **U**: `(c/σ)²` da **6347** con el σ del mínimo y
**1688** con el del extremo, y **ninguno cumple `σ(H*) = c`**. Se resuelve por **punto fijo**
—interpolar `log σ` contra `log n` y despejar—, que es lo que la fórmula quiere decir y que bajo
difusión coincide exactamente con ella. **La elección se tomó después de ver la firma y se declara
como lectura, no como enmienda.**

### El §2.2 queda resuelto: el piso de 1 950 se puede retirar

Era una pregunta abierta declarada en el preregistro («o se justifica el piso con un argumento
propio o se retira, y el reporte dice cuál de las dos»):

```
H*_ticks medido (punto fijo)  =  2233 ticks  =  127.7 s a nu = 17.49
piso heredado                 =  1950 ticks
```

**El `H*_ticks` medido ata, el piso no.** El piso no está multiplicando el requisito de datos, y
como es un literal heredado sin argumento vivo, **se retira**. Con el embargo medido la compuerta
exige 848 540 ticks y hay 1 031 154 → sigue pasando.

### §4.2 — el test de fuga pasa, y el control con poder pasa por un factor 1 000

Cuatro variantes, las cuatro al reporte como exige el preregistro. Métrica: `ΔLL/N` de M1 contra
**su propio** M0 en validación (los niveles absolutos no son comparables entre variantes porque
cambia el observable).

| variante | `ΔLL/N` vs su M0 | `G₀` | expectativa declarada | resultado |
|---|---|---|---|---|
| **correcta** | +4.290e-04 | +0.001101 | el valor a reportar | — |
| adelantada +1 | **+1.841e-03** | +0.002285 | debe subir claramente | **4.3×** ✔ |
| retrasada −10 | +4.098e-04 | +0.001083 | debe bajar | baja ✔ |
| **BARAJADA** | **+4.359e-07** | +0.010056 | debe caer al nivel de M0 | **0.1 %** ✔ |

La barajada —el único control con poder— cae a la milésima parte de la ganancia. **La ganancia de
M1 viene del emparejamiento libro-transacción**, no de la marginal de `e_t` ni de la
especificación.

### §4.3 — esto es OFI-L1 APROXIMADO, y los dos diagnósticos obligatorios

| diagnóstico | valor | umbral declarado |
|---|---|---|
| `q` excede la cantidad del mejor nivel del último snapshot | **2.77 %** | 20 % → **M1 es medición, no cota inferior** |
| residuo de reconciliación `ΔV` contra transacciones | **93.5 %** | — (es actividad de límite no observada) |

El 93.5 % no es un error: dice que el libro se mueve casi todo por actividad de límite invisible
entre snapshots, que es exactamente lo que el §4.3 advierte. La etiqueta «OFI-L1 aproximado» se
mantiene en todo el reporte.

Y la predicción falsable del signo se sostiene sobre el observable primario:
`E[y_mid·ε] = +0.017275 > 0`.

### §5.1 — `δ = 0` gana en validación, y `δ = 0.5` (MkII) NO es el máximo

| δ | M1 `ΔLL/N` | M2 `ΔLL/N` | `β` de M2 | `D` de M2 |
|---|---|---|---|---|
| **0.00** | **+3.173e-03** | **+7.763e-03** | −0.160 | 11.53 |
| 0.25 | +1.641e-03 | +6.356e-03 | −0.166 | 12.70 |
| **0.50** (MkII) | +4.244e-04 | **+4.387e-03** | −0.191 | 18.62 |
| 1.00 | +1.677e-05 | +1.571e-03 | −0.312 | 118.07 |

`δ = 0.5` da **1.77× menos** ganancia que `δ = 0` sobre el punto medio, y el orden es monótono en
los dos observables y los dos modelos. **El volumen no ayuda: lo que informa es el signo.** La
celda `δ = 0.5` se reporta aunque no gane, como manda el §5.1.

### §5.4 — el residuo es blanco en los tres modelos

| modelo | ρ₁ | ρ₂ | ρ₃ | ρ₁₀ | ρ₁₀₀ | ρ a rezago `embargo` |
|---|---|---|---|---|---|---|
| M0 | +0.0070 | +0.0086 | +0.0111 | +0.0145 | +0.0081 | +0.0022 |
| M1 | −0.0017 | −0.0003 | +0.0021 | +0.0056 | +0.0037 | +0.0020 |
| M2 | −0.0068 | −0.0053 | −0.0027 | +0.0010 | +0.0015 | +0.0019 |

Todos por debajo del umbral declarado `|ρ₁| < 0.05`, así que **ningún modelo queda descalificado
por residuo estructurado**. Nótese que M2 no deja residuo pese a `D = 11.5`.

### ⚠ §5.3 — `ω_G`: el estadístico salió NEGATIVO, y eso no es un resultado

Primera ejecución del contraste sobre `captura_v33`:

```
2*dLL observado = -149.078     nulo p95 = 1.159     p simulado = 0.6250
```

**M2-osc contiene a M2, así que `2·ΔLL` no puede ser negativo.** Y el 37.5 % de los sorteos del
nulo caían aún más abajo. Eso no es evidencia de nada: es que M2 y M2-osc arrancaban del **mismo
punto fijo** y aterrizaban en cuencas distintas, con Nelder-Mead sobre 5 parámetros.

Leído sin mirar el signo, el `p = 0.6250` decía «`ω_G = 0`, el propagador es monótono» y la
consecuencia declarada del §5.3 es borrar `ω_m,max` y `γ_ω` de `constantes_micelio.py`. **Se
habría borrado sobre un fallo del optimizador.**

**Corrección numérica, no de la regla:** el modelo grande se **siembra en la solución del
anidado** (`migracion_v32.contraste_omega_G`). No cambia la definición del estadístico; hace que
sea el que dice ser. Con eso:

```
2*dLL observado = +134.009     nulo p95 = 5.336     p simulado = 0.0000
anidamiento: obs >= 0 y nulo >= 0 en los 40 sorteos
```

⚠ **Pero el rechazo NO sostiene «oscilador forzado», y la guarda de banda lo dice.**

| | valor |
|---|---|
| `ω_G` ajustado | **+0.000009 rad/tick** |
| período implícito | **737 227 ticks = 42 151 s ≈ 11.7 h** a ν = 17.49 |
| `K` (rango ajustado) | 4 466 ticks |
| razón período/`K` | **165.1×** |

Dentro de `[0, K]` el coseno recorre el **0.606 % de un ciclo**: no oscila. Para `ω·K ≪ 1`,
`cos(ωτ+φ) ≈ cos φ − ωτ·sin φ`, o sea que el término está actuando como una **inclinación lenta
del núcleo**, no como un ciclo. Es la misma situación que la guarda de banda de la v2.1 §2 cazó
para `ω_m`, y por eso se implementó aquí la equivalente.

**`ω_G` queda SIN DECIDIR como frecuencia. No se borra ni se resucita nada en
`constantes_micelio.py`** — la consecuencia del §5.3 estaba condicionada a `ω_G = 0`, y eso no es
lo que salió; y el rechazo tampoco apoya un ciclo.

**Control nuevo, el 19, y lleva su propio negativo.** `migracion_v32.py` pasa a **20 controles**.
Y hay una lección en por qué los 18 anteriores no lo vieron: **el fallo sólo aparece con núcleo
CRECIENTE** (`β < 0`, `τ₀` pegada a su cota), que es justo el régimen donde cayó el mercado real y
que ningún control sintético anterior visitaba. Medido en el test:

```
2*dLL sin sembrar = -866.363      sembrado = +0.417      -> el control DISCRIMINA
```

Con `β = +0.5` los dos caminos coinciden y el test sería vacuo — por eso el control comprueba las
dos ramas y no sólo la buena.

### TABLA DE MEDICIONES — v3.2

**Compuerta (§1), salida literal de `--resumen`**

| magnitud | valor |
|---|---|
| transacciones | 1 031 155 en 16.38 h → **ν = 17.49 tx/s** |
| **tramo continuo más largo** | **1 031 155 ticks**, uno solo (corte a 300 s) |
| hueco máximo | 8.8 s |
| snapshots de libro | 8 097 860 (137.4 msg/s guardados) |
| OFI-L1 | 8 097 859 valores, mediana \|e\| = 0.1290, 100 % no nulos |
| `updateId` | **0 retrocesos, 0 repetidos** |
| `tr_maker`, `q`, cantidades de `bookTicker` | los tres persistidos |
| criterio en rojo | **sólo la vitalidad** (suspensión de la máquina), no la calidad del dato |

**Muestra y partición**

| magnitud | valor |
|---|---|
| observaciones tras alinear | 1 031 154 |
| precio | 64 794.4 – 65 482.7 USD |
| spread mediano | **0.1000 USD = 1 tick exacto** (p90 igual) |
| `y_mid` nulos | **96.67 %** |
| `y_tr` nulos | 72.61 % |
| `s_eff` de Roll (limpio) | 1.1280 USD/BTC (ρ₁ = −0.4260) |
| `c(u)` maker+maker | **26.0306 USD/BTC** ⚠ tarifas ASUMIDAS VIP 0 |
| embargo / `K` | **2233 / 4466** ticks |
| partición 60/20/20 | 618 692 / 203 998 / 203 998 (embargo verificado por test) |
| bloques de bootstrap en prueba | **19** en el diseño, **18** efectivos en el paso 2 (> 15) |
| `γ̂` de signos (entrenamiento) | **+0.7981** |

**Ajuste M2 elegido (δ = 0, punto medio, entrenamiento)**

```
G0 = +0.003541    tau0 pegada al limite    beta = -0.1596    f_inf = 0    sigma = 0.2583
```

`β` y `τ₀` **no están identificados individualmente** (§5.2): `τ₀` se pega a su cota inferior y
lo que los datos determinan es la curva `G(τ)` sobre `[0, K]`, que aquí es una potencia creciente
sin escala.

**LL/N fuera de muestra (conjunto de prueba, N = 203 998, abierto una sola vez)**

| modelo | LL/N |
|---|---|
| M0 | −0.09905503 |
| M1 | −0.09857358 |
| **M2** | **−0.09150251** |

`ΔLL/N`(M2 − M0) = **+7.553e-03**, IC95 bootstrap por bloques móviles
**[+6.857e-03, +8.084e-03]**, 18 bloques efectivos.

### Limitaciones de esta ejecución, declaradas

1. **La verosimilitud gaussiana está mal especificada en la marginal**: `y_mid` es **96.7 % ceros
   exactos**. Los `ΔLL/N` son comparables **entre modelos** (misma familia, misma muestra) pero su
   nivel no es interpretable como bondad de ajuste. El preregistro no anticipó esto.
2. **`c(u)` usa tarifas asumidas VIP 0**, criterio del §8 **no cumplido** (endpoint firmado, Modo
   LECTURA sin credenciales). Igual que en la v3.1.
3. **El nulo de `D` se regeneró con 24 sorteos** y el contraste de `ω_G` con 40; no más, por coste.
4. **Un solo tramo, 16.38 h, un solo régimen de ν.** La comparación con la sesión (e) ya muestra
   que la firma de volatilidad no es estable entre capturas.
5. **El §7 (Test C, en qué reloj vive el propagador) NO se ejecuta**: su §7.1 exige la captura
   completa y prohíbe expresamente correrlo con un tramo parcial.
6. **Nada de esto dice nada sobre la hipótesis (B)** —escala de decenas de minutos—, que sigue sin
   decidir desde la v2.2. El §11.1 lo deja escrito: el abandono estaba acotado a segundos y ni
   siquiera se activa.
7. **Se volvió al conjunto de prueba una segunda vez** para añadir descriptivos de `|μ̂|` (q50,
   q99, máx y el `c(u)` equivalente). **Ningún modelo, umbral ni regla cambió**, y ninguna decisión
   depende de esas cifras; se declara por la regla del §2.1.

### Qué queda decidido y qué no

**Decidido:**
- La regla de decisión se para en el **paso 3**, y así se reporta.
- El **piso de embargo de 1 950 se retira**: `H*_ticks` medido es 2233 y ata él.
- **`δ = 0`**: el volumen no aporta al forzamiento; informa el signo.
- El impacto transitorio de la v3.1 §2 **era rebote bid-ask**; sobre el punto medio no está.
- Este tramo es **super-difusivo** con control barajado limpio.

**No decidido, y por qué:**
- **Permanente contra transitorio sigue sin resolverse.** No por falta de datos: porque el
  estimador de `D` no tiene potencia a `γ̂ = 0.798`. Hace falta un estadístico cuyo nulo no se
  ensanche con la memoria del flujo — que es la misma lección de la sesión (e) por cuarta vez.
- **`ω_G` sin decidir.** El contraste corregido rechaza `ω_G = 0`, pero el `ω_G` ajustado tiene
  un período **165×** mayor que el rango ajustado: es una inclinación del núcleo, no un ciclo.
  **`ω_m,max` y `γ_ω` NO se borran de `constantes_micelio.py`, y tampoco se resucitan.**
- **La forma del núcleo** (creciente contra permanente contra transitorio) queda como medición
  puntual sin contraste con potencia detrás.
- **El §7**, por diseño del propio preregistro.

## Sesión 2026-08-10 (b) — v3.3 §6, §1, §2 y §3. La compuerta de tick grande CIERRA el marco

Ejecuta `ORDEN_TRABAJO_TICK_GRANDE_3_3.md` en su orden de prioridad. **Ningún estadístico nuevo**,
como manda su §7. `Micelio.py` sin cambios. Módulo nuevo: `tick_grande.py`, **19 controles**.

### ⚠ LO PRIMERO: `γ` significaba dos cosas, y la orden mezcló las dos

La tabla del §1 escribe `γ̂ = 0.798` y lo etiqueta «autocorrelación de signos». Ese número es el
que midió la v3.2 y es **`C(1) = corr(ε_t, ε_{t+1})`**, la autocorrelación **a rezago 1**.

Pero `H = (2−γ)/2 − β` viene del marco del propagador (Bouchaud, Gefen, Potters & Wyart 2004),
donde `γ` es el **exponente de decaimiento** de `C(ℓ) ~ ℓ^(−γ)`. Son dos cantidades distintas y no
hay razón para que coincidan. Medidas las dos sobre el mismo tramo:

| | valor |
|---|---|
| `C(1)` — lo que la v3.2 llamó `γ̂` | **+0.7981** |
| **exponente `γ` de `C(ℓ) ~ ℓ^(−γ)`** | **+0.5217** |

`C(ℓ)`: 0.798 (ℓ=1) · 0.762 (2) · 0.666 (10) · 0.586 (25) · 0.061 (2000). Ajuste log-log sobre
rezagos [10, 2000], 42 puntos, R² = 0.939.

**Con el `γ` correcto la predicción del §1 no se cumple.** La orden anticipaba
`β implícita ≈ +0.010`, «impacto esencialmente permanente». Sale:

| | valor |
|---|---|
| `H_p` (firma en ticks, n ∈ [256, 16384]) | 0.5907 |
| **`β` implícita = (2−γ)/2 − H** | **+0.1484** |
| `β` de difusividad = (1−γ)/2 | +0.2392 |

`β` implícita no es 0 (permanencia) ni 0.239 (difusividad): **cae en medio**, al 62 % del camino.
La conclusión que la orden esperaba —«no se puede rechazar impacto permanente»— **venía del `γ`
equivocado**.

### ⚠ Y el bootstrap por bloques NO SIRVE para estos estimadores. Se exhibe el fallo

Con bloque `= N^(1/2) = 1015`, que es lo que el §6 de la orden pide:

```
gamma  0.5217  ->  IC de bootstrap [0.7551, 0.9900]
H      0.5907  ->  IC de bootstrap [0.0904, 0.1826]
```

**Los intervalos no contienen el punto estimado.** No es ruido: remuestrear bloques de longitud
`b` produce una serie cuya dependencia **muere en `b`**, y los dos estimadores ajustan sobre
rezagos **mayores** que `b` (hasta 2 000 y 16 384). Miden la longitud del bloque, no el mercado.

Sustituto correcto: **submuestreo contiguo**, que conserva la estructura temporal de cada réplica.
8 sub-series de 77 336 ticks:

| sub-serie | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| `γ` | 0.384 | 0.508 | 0.624 | 0.616 | 0.438 | 0.560 | 0.724 | 0.585 |
| `H` | 0.665 | 0.725 | 0.689 | 0.542 | 0.432 | 0.601 | 0.639 | 0.627 |
| **`β`** | +0.143 | +0.021 | **−0.000** | +0.150 | **+0.349** | +0.119 | **−0.001** | +0.081 |

```
beta media +0.1076   sd 0.1152   intervalo t OPTIMISTA [+0.0135, +0.2017]
```

**Tres sub-series dan `β ≈ 0` exacto (permanencia) y una da +0.349 (más allá de difusividad).**
El intervalo excluye tanto 0 como 0.239, pero **por muy poco y con un intervalo que es cota
inferior de la incertidumbre** (supone sub-series independientes, y con memoria larga no lo son).

**Lectura honesta: el §1 tampoco cierra la pregunta.** Da un `β` central positivo y pequeño, con
una heterogeneidad entre tramos que abarca desde permanencia hasta pasada la difusividad. Es la
misma inestabilidad que la v3.0 encontró en `k` y que la propia orden anota (`β` cambiando de signo
entre tramos). Sólo que ahora está cuantificada y con el `γ` correcto. Sólo 8 réplicas:
**subpotenciado**, declarado.

**Los tres `β` juntos, sin elegir** (§1.1 punto 4):

```
implicita por H y gamma       +0.1484
ajuste M2 sobre punto medio   -0.1600
ajuste M2 sobre precio trans. +4.6900
```

La discrepancia es el resultado.

### ⚠ §2 — `η̂ = 0.62 ≥ 1/2`: el marco de tick grande NO APLICA. §4 y §5 no se ejecutan

La compuerta que la orden puso para poder rechazar barato ha hecho justo eso.

| variante | `η̂` | `N_c` | `N_a` |
|---|---|---|---|
| **sólo saltos de 1 tick** (la del artículo) | **0.6200** | 98 718 | 79 616 |
| todos los saltos (control del modo de fallo 3) | 1.1614 | 197 468 | 85 011 |

El control del modo de fallo 3 funciona: contar saltos de más de un tick casi **duplica** `η̂`. Con
la definición correcta sale 0.62, por encima de 1/2.

**Consecuencia declarada por la propia orden: `§4` (propagador sobre el precio eficiente `X_t`) y
`§5` (la predicción falsable del agotamiento de cola) NO se ejecutan.** Y ese `§5` era «el apartado
con más valor de información del documento». Se cierra sin gastarlo, que es exactamente para lo que
sirve una compuerta barata.

⚠ **Y la reserva del §2 sobre el tick relativo queda confirmada por el número:**

```
alpha/P = 0.10 / 65 076 = 1.54e-6 = 0.0154 pb por tick
```

En renta variable «tick grande» significa `α/P` de 1e−4 a 1e−3: **dos o tres órdenes de magnitud
mayor**. BTCUSDT tiene la horquilla clavada en 1 tick **y** un tick relativo diminuto, y esa
combinación no está en los artículos. Con `η̂ = 0.62` la cuestión es académica aquí, pero se anota.

27.4 % de los ticks mueven el precio de transacción (282 480 de 1 031 153).

### §6 — `N_eff` medido, y el paso 2 de la v3.2 SOBREVIVE

El §6 propone `N_eff = N^γ` con la `γ` de los signos. **Se mide en vez de suponerse**, y sobre la
serie que toca — la de `ΔLL`, no la de signos, que no tienen por qué compartir exponente. Método:
ajustar `Var(media de bloque)` contra el tamaño de bloque y evaluar en `b = N`.

| serie | `γ_eff` | R² | `N` | `N_eff` | pérdida | inflación del error estándar |
|---|---|---|---|---|---|---|
| **`ΔLL` de prueba** | **0.7421** | 0.993 | 199 532 | **8 574** | **23.3×** | **3.15×** |
| signos `ε` | 0.3539 | 0.927 | 618 692 | **112** | **5 518×** | 115.4× |
| incrementos de precio | 0.8630 | 0.981 | 618 691 | 99 510 | 6.2× | 0.71× |

**Recálculo del IC del paso 2:**

```
publicado en la v3.2 (bloques de 5*embargo) : [+6.857e-03, +8.084e-03]
corregido por memoria larga                 : [+6.649e-03, +8.456e-03]
```

**El paso 2 sigue pasando.** El error estándar se ensancha 3.15× pero la media está a ~15 errores
estándar de cero, así que el intervalo sigue excluyendo 0 con holgura. La preocupación del §6 era
legítima y el resultado es que no muerde **en este IC**.

⚠ **Donde sí muerde es en los estimadores basados en signos: `N_eff = 112` sobre 618 692 ticks,
pérdida de 5 518×.** Eso explica mejor que cualquier defecto de diseño por qué han muerto cuatro
estadísticos: todos dependían de la relación entre signo y precio, y ahí el tamaño muestral
efectivo es de tres cifras. **No se arregla capturando más** — con `γ_eff = 0.354`, multiplicar
`N` por 10 sólo multiplica `N_eff` por 2.3.

⚠ **Nota de procedimiento:** este recálculo toca el conjunto de prueba. El §7 prohíbe **reabrirlo
para un veredicto nuevo**; el §6 **ordena** recalcular un IC ya publicado. Se recomputó el mismo
estadístico sobre los mismos modelos, sin ajustar ni elegir nada, y no se deriva ninguna decisión
nueva de él.

### §3 — La sospecha del precio NO se cumple, y la aritmética en pb

La orden sospecha que el tramo pudiera haber corrido cerca de 96 000, con lo que `c(u)` serían 38.4
y el factor pasaría de 3.46 a 5.1. **Medido: el tramo va de 64 794.4 a 65 482.7, mediana 65 076.4.**
`c(u) = 26.03 USD/BTC` es correcta para este tramo y **el factor 3.46 se mantiene**.

| magnitud | USD/BTC | pb | ticks |
|---|---|---|---|
| 1 tick | 0.100 | **0.0154** | 1.0 |
| comisión ida y vuelta maker | 26.031 | **4.0000** | 260.3 |
| `\|μ̂\|` q90 medida (v3.2) | 11.297 | 1.7360 | 113.0 |
| **umbral `1.5·c(u)`** | **39.046** | **6.0000** | **390.5** |

**Hay que predecir 390 ticks de movimiento para pagar la ida y vuelta. Un agotamiento de cola da
1 tick.** La capa de microestructura opera **390× por debajo** del umbral de rentabilidad. Y en pb
la tabla ya no envejece con el precio: la comisión son 4 pb sea cual sea el nivel de BTC.

### Un control que falló y lo que enseñó

El control 2 de `tick_grande.py` —¿recupera el estimador un `γ` conocido?— falló al primer intento
con **R² = 0.87**, por debajo del umbral de 0.95 que yo mismo había puesto. **El fallo era del
generador, no del estimador**: una suma finita de seis AR(1) es ley de potencias sólo a trozos, con
ondulaciones entre las `τ`. Bajar el umbral para que pasara habría sido ajustar el control al
resultado. Se sustituyó por **ruido gaussiano fraccionario**, cuyo exponente sí es verdad conocida
(`C(ℓ) ~ ℓ^(2H−2)`, o sea `γ = 2−2H`).

Y así apareció algo que cambia cómo hay que leer el R² sobre dato real:

| `H` verdadera | `γ` verdadera | `γ` medida | **R²** | `C(1)` |
|---|---|---|---|---|
| 0.75 | 0.50 | 0.5436 | **0.977** | 0.243 |
| 0.65 | 0.70 | 0.6647 | **0.660** | 0.129 |

**Con memoria débil el R² se desploma mientras `γ` sigue bien recuperada**, porque `C(ℓ)` es más
pequeña y el ruido muestral domina el logaritmo. Yo había escrito una compuerta «R² < 0.90 → la
relación no aplica»: **habría descartado un caso donde el estimador funciona perfectamente.** El R²
pasa a reportarse como diagnóstico, nunca como compuerta.

### Qué queda de la v3.3

**Cerrado:**
- **§2**: `η̂ = 0.62` → el marco de tick grande no aplica. **§4 y §5 quedan cancelados por la
  compuerta**, no pendientes.
- **§3**: precio del tramo verificado, `c(u)` correcta, factor 3.46 confirmado, tabla en pb.
- **§6**: `N_eff` medido en tres series; el paso 2 de la v3.2 sobrevive al IC corregido.

**Ejecutado pero sin cerrar:**
- **§1**: la relación se verificó en sus dos límites y `γ` se midió bien, pero `β` implícita sale
  +0.148 con sub-series entre −0.001 y +0.349. **No separa permanencia de difusividad**, y el
  intervalo disponible es una cota inferior de la incertidumbre con sólo 8 réplicas.

**Lo que esto añade a la lectura de fondo:** las tres rutas de la v3.3 no rescatan la economía —el
§3 de la orden ya lo decía— y la que podía cerrar la ciencia (§1) tropieza con la misma
heterogeneidad entre tramos que viene apareciendo desde la v3.0. El `N_eff = 112` de los signos es
la explicación cuantitativa de por qué.

## Sesión 2026-08-10 (c) — v4.1 §2 y §1. Las dos curvas NO se cruzan donde hay potencia

Ejecuta `ORDEN_TRABAJO_HORIZONTE_4_1.md` en su orden de prioridad: §11 (compuerta de
infraestructura) → §2 (`H_p`) → §1 (las dos curvas). `Micelio.py` sin cambios. Módulo nuevo:
`horizonte.py`, **14 controles**.

### ⚠ §11 PRIMERO: la captura estacional muere por BATERÍA, no por suspensión por inactividad

`captura_estacional` es «el único activo del proyecto que puede dar potencia al §1 por encima de
30 min», y lleva **dos apagones sólo hoy**. Los dos capturadores tienen **exactamente los mismos
huecos**, así que es la máquina y no los procesos:

```
Evento 524  10:52:05  "Critical desencadenador de bateria cumplido"
Evento 42   10:52:05  "El sistema esta entrando en suspension"
Bateria: estado 1 = DESCARGANDO      nivel critico configurado: 2 %
```

**El portátil está con batería, no enchufado.** `STANDBYIDLE = 0` —que sí está configurado
correctamente, verificado— **no protege de esto**: la suspensión por batería crítica pasa por
encima de todo, y `impedir_suspension()` (`SetThreadExecutionState`) tampoco la evita. Es una causa
**distinta** de la que se corrigió el 2026-08-09.

Cobertura acumulada del estacional: **1 080 min = 0.75 días equivalentes**, 147 casillas vacías de
168, mínimo por casilla 0 contra los ≥ 30 que pide su compuerta. Apagones de hoy: 01:00–06:00 y
10:52–16:47 local.

⚠ **Con 21 días continuos, `H = 1 h` da ~500 ventanas y `H = 4 h` da ~126. Cada apagón parte el
tramo continuo, y es el tramo continuo el que fija `H_max ≈ T_continuo/30`.** Sin resolver la
alimentación, el §1 no puede pasar de piloto en la banda que importa.

### ⚠ §2 — `H_p = 0.591` NO sobrevive, y la compuerta resuelve por BANDA

Firma de volatilidad sobre las **cuatro** capturas, **reloj de pared** (que es donde vive la curva
requerida: comisión y financiación se pagan en tiempo de calendario), con **control barajado en
todas**, ajuste en `H ∈ [60, 600] s`:

| captura | `ν` | `H_p` real | barajado | **corregido** | \|sesgo\| | efecto | fiable |
|---|---|---|---|---|---|---|---|
| `captura_larga` | 38.97 | 0.407 | 0.401 | 0.506 | 0.099 | 0.093 | **NO** |
| `captura_v31b` | 5.76 | 0.414 | 0.543 | **0.371** | 0.043 | 0.086 | sí |
| `captura_v32` | 10.13 | 0.527 | 0.520 | **0.508** | 0.020 | 0.027 | sí |
| `captura_v33` | 17.49 | 0.548 | 0.496 | **0.552** | 0.004 | 0.048 | sí |

**El 0.591 de la v3.3 era del reloj de TICKS.** En reloj de pared, sobre la misma captura y con su
control, `captura_v33` da **0.552**. Y entre capturas la banda es **[0.371, 0.552]**, que **cruza
0.5**: la superdifusión no es una propiedad estable.

⚠ **Dos cosas que hay que declarar sobre cómo se llegó a esa tabla.**

1. **La rejilla de horizontes se densificó a mitad de camino**, y no es cosmético. Con la rejilla
   original el ajuste de `[60, 600]` tenía **4 puntos** y el propio control barajado —que debe dar
   0.5 exacto— salía entre **0.254 y 0.528**. Ahí la corrección mueve más ruido que señal. Con 11
   puntos los controles quedan en 0.401–0.543.
2. **El criterio de fiabilidad se escribió DESPUÉS de ver la primera corrida**, y la fila que
   descarta es justo la que rompe la monotonía con `ν`. Eso es selección post-hoc y se reporta con
   las dos correlaciones, sin elegir:

```
corr(log nu, H_p)  con las 4 filas : +0.6821   (n = 4)
corr(log nu, H_p)  solo fiables    : +0.9616   (n = 3)
```

**Con `n ≤ 4` una correlación no es evidencia**: para `n = 3` hace falta `|r| > 0.997` para
`p < 0.05`. La rama de correlación del §2.2 **no se puede invocar**, y se aplica la **tercera**:
la curva requerida va como **banda `[0.371, 0.552]`** y el §1 se lee contra la banda entera.

`σ₁` reexpresada en pb (§2.1 punto 3), recalculada por cada `H_p` sobre la firma de `captura_v33`:

| `H_p` | 0.371 | 0.500 | 0.552 |
|---|---|---|---|
| `σ₁` [pb·s^(−H_p)] | 0.7356 | 0.4061 | 0.3196 |

### ⚠ §1 — LAS DOS CURVAS NO SE CRUZAN, y donde podrían cruzarse no hay potencia

Predictor deliberadamente tonto (§1.3): regresión de `r_{t→t+H}` del punto medio sobre flujo
firmado acumulado en `{H/4, H/2, H, 2H}`. Sin `G(τ)`, sin `τ₀`, sin `β`. Ventanas **no solapadas**
(verificado por test). Entrenado en el 60 % de entrenamiento, medido en validación. **El conjunto
de prueba de la v3.2 NO se abrió.**

| `H` | `n_vent` | **`R²_medido`** | `R²_barajado` | `R²_req` (0.371) | (0.500) | (0.552) | |
|---|---|---|---|---|---|---|---|
| 60 s | 112 | **−0.0008** | −0.0070 | 68.72 % | 78.40 % | 82.67 % | |
| 120 s | 56 | **+0.0019** | **+0.0146** | 41.09 % | 39.20 % | 38.46 % | |
| 300 s | 23 | −0.0264 | −0.0051 | 20.82 % | 15.68 % | 13.99 % | piloto |
| 600 s | 11 | −0.0554 | −0.0195 | 12.45 % | 7.84 % | 6.51 % | piloto |
| 900 s | 8 | −0.0763 | +0.0062 | 9.21 % | 5.23 % | 4.16 % | piloto · extrapola |
| 1800 s | 4 | −0.0461 | **+0.1292** | 5.51 % | 2.61 % | 1.93 % | piloto · extrapola |
| 3600 s | 2 | insuficiente | | | | | |

**Aplicando el §1.4 literalmente, el desenlace es doble:**

1. **Donde hay potencia (60 y 120 s, `n ≥ 30`): fallo categórico.** `R²_medido ≈ 0` contra un
   requerido de **38–83 %**. No es una brecha que un predictor mejor cierre: aunque un predictor
   perfecto alcanzara `R² = 0.5`, seguiría sin pagar a 60 s. Y en las dos filas el barajado es
   comparable o **mayor** que el real, que por el §1.4 se lee como «la medición no tiene potencia
   a esa `H`» — o, dicho de otro modo, no hay nada que medir.
2. **Donde el requerido se vuelve plausible (≥ 15 min, 1.9–5.5 %): sólo 4–8 ventanas.** Por el
   §1.4 eso es **NO DECIDIBLE**, no un fallo — y la única acción admisible es esperar a
   `captura_estacional`. Que es justo lo que el §11 está bloqueando.

⚠ **Y una reconciliación que hay que escribir para que las dos cifras no se lean como
contradicción.** La v3.2 reportó `R² = 0.014` a 128 s y aquí sale ~0. **No es la misma cantidad.**
El `R²` de la v3.2 predice `Δp_t` con un núcleo que incluye `h(0)·x_t`, o sea **el signo de la
transacción que está ocurriendo**: es ajuste contemporáneo, y no es operable. El `R²` de aquí
predice el retorno **futuro** a `H` desde flujo **pasado**, que es lo único que se puede negociar.
Que uno sea 1.4 % y el otro 0 no es inconsistencia: es la diferencia entre explicar y predecir.
(El `μ̂` del paso 3 de la v3.2 ya usaba sólo rezagos ≥ 1, así que aquella parte estaba bien.)

### Dos controles míos que fallaron, y lo que enseñaron

**El control de la curva requerida pasa a la primera** y reproduce la tabla del §11.1 del
`PREREGISTRO_3_2.md` con tres decimales: 41.77 % contra 41.8 a 103 s, 7.17 contra 7.2 a 10 min,
1.20 contra 1.2 a 1 h. La forma está verificada contra lo ya publicado.

⚠ **El control de superdifusión estaba mal planteado y era mío.** Fabricaba superdifusión con una
**tendencia lineal**, y una tendencia constante **no aparece en la firma de volatilidad**: añade la
misma constante a todos los incrementos y `std(x + c) = std(x)`. Superdifusión es incrementos
**positivamente autocorrelacionados**, no deriva. Sustituido por movimiento browniano fraccionario,
cuyo `H` es verdad conocida:

| fBm | real | barajado | corregido |
|---|---|---|---|
| `H = 0.65` | 0.6833 | 0.5372 | **0.6461** |
| `H = 0.35` | 0.3977 | 0.5224 | **0.3754** |

Recupera la verdad por los dos lados y el barajado vuelve a 0.5 en los dos, que es lo que hace del
control un discriminador y no un adorno.

### Qué queda de la v4.1

**Hecho:**
- **§11**: causa de los apagones identificada (**batería**, no inactividad) y cuantificada.
- **§2**: `H_p` sobre las cuatro capturas, dos relojes, control en todas, `σ₁` en pb. Resuelve por
  **banda [0.371, 0.552]**. `H_p = 0.591` retirado.
- **§1**: las dos curvas, con control barajado y ventanas no solapadas. **No se cruzan donde hay
  potencia; donde podrían, no hay potencia.**

**Pendiente, en el orden del §10:**
- **§3** (correcciones a la v3.3): tercer estimador de `γ` por Whittle/GPH, `β` implícita como
  curva contra la ventana de ajuste, tabla de `N_eff` con la definición coherente, la identidad
  `β = 0 ⟺ γ = 2−2H`, y reclasificar `D = 11.53`.
- **§6** (micro-precio): sólo la fracción de ceros contra el 96.67 %.
- **§5** (`C_respaldo`) y **§4** (`η̂` con barrido de colapso), en ese orden.

## Sesión 2026-08-12 — v4.1 §3. No hay «la» `γ`, y `β` implícita no mide nada

Ejecuta el §3 completo de `ORDEN_TRABAJO_HORIZONTE_4_1.md`. **Ningún quinto estadístico**, como
manda su §7; el §3.4 **degrada** una afirmación en vez de intentar mejorarla. `Micelio.py` sin
cambios. Módulos: `tick_grande.py` pasa a **28 controles**, ejecutor nuevo `correcciones_v41.py`.

Todo se corre sobre **dos series**, que es lo que la v3.3 no pudo hacer:

| serie | `N` | `ν` | duración |
|---|---|---|---|
| `captura_v33` entrenamiento — la que la v3.3 midió | 618 692 | 17.49 tx/s | — |
| **`captura_estacional`, tramo continuo más largo** — réplica independiente | **1 829 242** | 21.78 tx/s | **23.33 h** |

### ⚠ LO PRIMERO: la orden esperaba que el §3 fuera cosmético. No lo es

El §3 se escribió como «correcciones baratas, ninguna cambia un veredicto». Dos de las cuatro
**retiran una cantidad del vocabulario del proyecto**, y las dos replican en la captura nueva.

### §3.1 — El tercer estimador no arbitra entre los dos primeros: los explica

`gph()` y `whittle_local()` en `tick_grande.py`, con la traducción `γ = 1 − 2d` escrita en
`PREDICCION_SESGO_GAMMA_4_1.md` (commit `f1fd84f`, **congelado antes de calcular ninguna cifra**).

| método | dominio | `γ` en `captura_v33` | `γ` en `estacional` |
|---|---|---|---|
| A — regresión log-log de `C(ℓ)`, `[10, 2000]` | tiempo | +0.5217 | +0.6336 |
| B — escalado de `Var(media de bloque)` | tiempo | +0.3539 | +0.4201 |
| **C1 — GPH**, `m = N^0.5` | frecuencia | **+0.4187** ±0.046 | **+0.6263** ±0.035 |
| **C2 — Whittle local**, `m = N^0.5` | frecuencia | **+0.4216** ±0.036 | **+0.6482** ±0.027 |
| `C(1)` (lo que la v3.2 llamó `γ̂`) | — | +0.7981 | +0.8093 |

**Las cuatro predicciones congeladas, y cómo salieron:**

| | predicción | `captura_v33` | `estacional` |
|---|---|---|---|
| **P1** | `γ_espectral < γ_A` (sesgo por centrado con la media muestral) | **se cumple** | **REFUTADA** (LW 0.6482 > 0.6336) |
| **P2** | `γ_espectral ∈ [0.28, 0.46]`, o sea más cerca de B | **se cumple** | **REFUTADA** (0.63 y 0.65) |
| **P3** | `γ̂` **crece** con la banda `m` | **REFUTADA** | **REFUTADA** |
| **P4** | los dos recuperan `d` conocida sobre fGn con error < 0.05 | **se cumple** (peor 0.0117) | ídem |

**P1 y P2 no replican.** Sobre la captura nueva los dos estimadores espectrales caen **encima de
A**, no de B. La hipótesis de la orden —«0.5217 es el sesgado y 0.3539 el más cercano»— **no
recibe apoyo**: se cumple en un tramo y falla en el otro.

⚠ **P3 falló, y mi mecanismo declarado estaba al revés.** Escribí que la estructura de corto
alcance «aplana el espectro lejos del origen», lo que reduciría `d̂` y subiría `γ̂`. Es falso para
autocorrelación **positiva**: un proceso con `C(1) = 0.80` tiene espectro que **decae** al alejarse
del origen, así que ampliar la banda **empina** la pendiente, sube `d̂` y **baja** `γ̂`. Eso es lo
que se midió, en las dos series y en los dos estimadores.

### ⚠ Y la consecuencia de P3 es el hallazgo de la sesión: **no existe «la» `γ`**

| `α` (`m = N^α`) | `m` | GPH v33 | LW v33 | GPH estacional | LW estacional |
|---|---|---|---|---|---|
| 0.4 | 207 / 319 | +0.5116 | +0.5440 | +0.6513 | +0.6917 |
| 0.5 | 786 / 1 352 | +0.4187 | +0.4216 | +0.6263 | +0.6482 |
| 0.6 | 2 984 / 5 719 | +0.0484 | −0.0258 | +0.2646 | +0.2240 |
| 0.7 | 11 324 / 24 187 | **−0.1796** | **−0.2441** | **−0.1653** | **−0.2577** |

`γ̂` recorre **0.69 a 0.95** según la banda, sobre una disputa cuyo rango entero es 0.17. Y `γ < 0`
significa `d > 0.5`, o sea **fuera del rango de estacionariedad** que los propios estimadores
suponen.

**El control que decide cómo se lee esto**, y sin él el resultado sería ambiguo: el **mismo
barrido sobre fGn de `γ` conocida (0.50) y la misma longitud**:

| | salto máximo entre bandas |
|---|---|
| fGn, verdad conocida (v33 / estacional) | **0.1274 / 0.0578** |
| **dato real** | **0.4474 / 0.4817** |
| razón | **3.5× / 8.3×** |

Sobre verdad conocida el barrido es **plano** (0.44–0.56 en las cuatro bandas); sobre el dato se
derrumba. **El fallo es de la serie, no del estimador: el flujo de órdenes de este mercado no
tiene un régimen de escala único.** La discrepancia del §3.1 nunca fue un problema de estimador.

### §3.2 — `N_eff` coherente, y las tres consecuencias

`inflación = √(N/N_eff)` debe satisfacerse; `N^γ_eff` descarta la constante del escalado y no la
satisface en ninguna fila. Corregido en `ic_memoria_larga`, que ahora devuelve las dos.

| serie | `N` | inflación | **`N_eff` = `N`/infl²** | `N^γ` (v3.3) |
|---|---|---|---|---|
| `ΔLL` de prueba (v3.3, aritmética) | 199 532 | 3.15 | **20 109** | 8 574 |
| **signos `ε`** v33 | 618 692 | 115.40 | **46.5** | 112 |
| **signos `ε`** estacional | 1 829 242 | 115.99 | **136.0** | 427 |
| incrementos de precio v33 | 618 691 | **0.71** | **1 238 283** | 99 510 |
| incrementos de precio estacional | 1 829 241 | 1.14 | 1 415 781 | 167 153 |

1. **Signos: `N_eff = 46.5`, no 112.** El diagnóstico de la v3.3 sobre por qué murieron cuatro
   estadísticos es **más fuerte** de lo que se reportó, no más débil. Y la **inflación es
   notablemente estable entre capturas** (115.40 contra 115.99) pese a triplicarse `N`.
2. ⚠ **La antipersistencia NO replica.** En v33 `infl = 0.71 < 1` y `N_eff > N`, que es la firma
   del rebote bid-ask (`ρ₁(retornos) = −0.216`, v3.0) y que la tabla de la v3.3 escondía tras una
   columna de «pérdida». En el tramo estacional `infl = 1.14 > 1`: pérdida de 1.3×, **no**
   antipersistencia. Es otra cantidad que cambia entre tramos.
3. **La fila de `ΔLL` no se recalculó**: sale por aritmética sobre la inflación 3.15 ya publicada.
   El §7 prohíbe reabrir el conjunto de prueba y aquí no hace falta. El paso 2 sigue pasando.

⚠ **Regla del proyecto, escrita:** `Var(media) ∝ N^(−γ)` es el resultado de Beran **para la media
muestral**. No se transporta a estimadores de exponentes, a pendientes de regresión ni a razones
de verosimilitud. **`N_eff` se mide por estadístico, midiendo el escalado de ese estadístico. No
se importa de la serie de signos.**

### ⚠ §3.3 — `β` implícita se mueve más que la distancia entre las hipótesis. En las dos series

La identidad que el §3.1 obliga a escribir, y que cambia cómo se lee todo el §1 de la v3.3:

```
beta_implicita = (2 - gamma)/2 - H = 0   <=>   gamma = 2 - 2H
```

que es exactamente la relación entre exponente de ACF y exponente de Hurst del **ruido gaussiano
fraccionario**. **`β` implícita no mide el núcleo**: mide cuánto se desvía el par (signos, precio)
de la relación fGn. Si los dos comparten proceso de memoria larga, sale 0 **mecánicamente**. El
control 7 de `tick_grande.py` lo demuestra sobre fGn puro —donde no hay propagador ninguno— y da
`β = −0.039`.

`γ` y `H` **ajustadas en la misma ventana de escala** (la v3.3 las tomó de ventanas con menos de
una década de solapamiento):

| ventana [ticks] | `γ` v33 | `H` v33 | **`β` v33** | `γ` est. | `H` est. | **`β` est.** |
|---|---|---|---|---|---|---|
| [16, 128] | +0.263 | 0.498 | **+0.370** | +0.322 | 0.603 | **+0.236** |
| [32, 256] | +0.433 | 0.584 | +0.200 | +0.512 | 0.679 | +0.065 |
| [64, 512] | +0.644 | 0.656 | +0.022 | +0.789 | 0.714 | −0.108 |
| [128, 1024] | +0.789 | 0.639 | −0.033 | +1.067 | 0.665 | **−0.198** |
| **[256, 2000]** común | +0.726 | 0.606 | **+0.031** | +0.839 | 0.588 | **−0.007** |
| [512, 4096] | +0.653 | 0.591 | +0.083 | +0.776 | 0.558 | +0.054 |
| [1024, 8192] | +0.536 | 0.576 | +0.156 | +0.485 | 0.528 | +0.229 |
| [2048, 16384] | +0.570 | 0.623 | +0.092 | +0.487 | 0.467 | **+0.290** |
| *mezcla de la v3.3* | +0.522 | 0.598 | *+0.142* | +0.634 | 0.537 | *+0.146* |

```
                       captura_v33        estacional
recorrido de beta        0.4030             0.4881
distancia entre hipotesis 0.2050            0.1583      (0 contra beta_dif medio)
```

**En las dos series el recorrido dobla o triplica la distancia entre las hipótesis.** La regla de
lectura del §3.3, escrita antes: *«si `β` implícita se mueve más que la distancia entre las dos
hipótesis al variar la ventana, el §1 de la v3.3 no está midiendo nada y así se reporta.»*

**Se reporta: el §1 de la v3.3 no está midiendo nada.** La `β` implícita de +0.148 que aquella
sesión publicó es un punto arbitrario de una curva que cruza el cero y llega a ±0.29 sin salir
del rango de escalas del propio experimento.

Nótese además que en la ventana común `β` sale ≈ 0 en las dos series (+0.031 y −0.007) — que es
justo lo que la identidad fGn predice **sin que haya propagador**, así que ni siquiera eso apoya
la permanencia.

### §3.4 — `D = 11.53` reclasificado en los tres documentos donde se citaba

`β = −0.160` con `τ₀` pegada a su cota es firma de **mala especificación del estimador**: `β` se
reporta en la literatura en `(0, 1)`, y un núcleo que crece sin cota sobre `[0, K]` implica
impacto de mercado creciente indefinidamente, económicamente imposible (arbitraje ilimitado). El
propio §5.2 de la v3.2 ya decía que `β` y `τ₀` no están identificados por separado.

`D = 11.53` deja de citarse como «núcleo creciente medido» y pasa a **«ajuste no identificado en
el régimen `τ₀` en cota»**. Aplicado en `CLAUDE.md` (retractaciones, tabla de la v3.2, apartado
`D > 1`, y la frase de la firma de volatilidad que lo usaba como corroboración) y en
`TRASPASO_SESION_2026-08-10.md`. **No se abre un quinto estadístico.**

### Controles nuevos (`python tick_grande.py --autotest` → 28/28)

| # | control | resultado |
|---|---|---|
| 6 | GPH y LW recuperan `d` conocida sobre fGn (`H` = 0.85…0.55) | peor error **0.0117** |
| 6b | sin memoria, `d = 0` dentro de 2 errores estándar | GPH −0.0003, LW −0.0064 |
| 7 | `β = 0 ⟺ γ = 2−2H`, y sobre fGn puro `β ≈ 0` sin propagador | **−0.0392** |
| 8 | `N_eff` coherente satisface `infl = √(N/N_eff)`; serie antipersistente marcada | exacto a 1e−9 |

⚠ **Un control mío falló y el fallo era del umbral.** Puse `|d| < 0.06` sobre ruido blanco sin
calcular el error estándar; con `N = 2^16` y `m = 256` el error estándar de GPH es **0.0401**, así
que un sorteo suelto en 0.0684 está a 1.7 σ de cero. **Un umbral por debajo del propio error del
estimador rechaza estimadores correctos.** Sustituido por media de 8 sorteos contra `2·se/√k`,
que es **más exigente** que el original: si el estimador estuviera sesgado, promediar lo dejaría
lejos de cero mientras el umbral se estrecha. Es el tercer umbral de este proyecto puesto «a ojo»
que hubo que medir (los anteriores: `UMBRAL_TAYLOR_JACOBIANO` en la v2.1 y la compuerta de `R²`
en la v3.3).

### Qué queda de la v4.1 tras el §3

**Hecho:** §11, §2, §1 (sesión anterior) y **§3 completo**, con los cinco criterios de aceptación
del §9 cubiertos y replicado en una segunda captura.

**Lo que el §3 retira del vocabulario:**
- **`γ` no es un número de este mercado.** Depende de la banda por 0.7–0.95. Cualquier fórmula que
  la use —incluida `H = (2−γ)/2 − β`— hereda esa indeterminación.
- **`β` implícita no es una medición del núcleo.** Es un contraste de coherencia con fGn, y en el
  rango de escalas del experimento no discrimina.
- **`D = 11.53` no describe la forma del núcleo.**

**Pendiente, en el orden del §10:** §6 (micro-precio, sólo la fracción de ceros es barato), §5
(`C_respaldo`, después de que el §1 diga a qué horizonte), §4 (`η̂` con barrido de umbral de
colapso — deuda de reporte, no decisión).

⚠ **Y lo que de verdad desatasca el proyecto:** el tramo de **23.33 h continuas** da **93
ventanas a 15 min, 46 a 30 min y 23 a 1 h**, contra las 8/4/2 con las que el §1 se declaró NO
DECIDIBLE. La banda donde el `R²` requerido baja a 1.9–5.5 % ya tiene potencia. **Rehacer el §1
sobre ese tramo es lo siguiente.**

## Sesión 2026-08-27 — Adenda C: el §1 medido por covarianza. Módulo `identidad.py`

Ejecuta `ADENDA_C_REFORMULACION_S1_4_1.md`. **No es una enmienda al criterio**: el §1.4 se
mantiene palabra por palabra y lo que se sustituye es el **instrumento**. `Micelio.py` sin
cambios. El conjunto de prueba de la v3.2 no se abre.

### ⚠ LO PRIMERO, porque es lo que el operador preguntó: `σ` con `ν` SIGUE EN PIE

La adenda argumenta que un `R²` fuera de muestra tiene suelo de ruido `≈ k/n` y que ninguna
fila se lee si la cantidad a medir no supera **3× `q95(|R²_nulo|)`**. Aplicada la misma regla,
con **200 sorteos**, al hallazgo de predecir volatilidad con actividad:

| | `R²` medido | `q95(\|R²_nulo\|)` | razón | ¿se lee? |
|---|---|---|---|---|
| `σ(t)` desde 3 rezagos de (`ν`, `σ`) | **+0.6153** | 0.0798 | **7.7×** | **sí** |
| **`ν(t)` → `σ(t+1)`** | **+0.5464** | 0.0508 | **11×** | **sí** |

El nulo es **rotación circular** (conserva la memoria de la serie), que da un suelo **14× más
ancho** que el barajado (0.0798 contra 0.0058). Aun contra ese, pasa con holgura. La crítica de
la adenda apunta al régimen de `n` pequeño del método de trocear en ventanas; el trabajo de
`σ`/`ν` tiene `n = 3 081` con `k = 7` y está dos órdenes por encima del suelo. **No hay nada
que retractar ahí.**

### La identidad, y sus controles con verdad conocida (6/6)

```
R2(H) = Corr(eps_t, p_{t+H} - p_t)^2 = R_cum(H)^2 / (Var(eps) * sigma_r(H)^2)
```

Reproduce el §C.3.3 sobre serie sintética con señal inyectada de magnitud conocida
(`N = 400 000`):

| señal inyectada | **identidad** | ventanas no solapadas, con ajuste fuera de muestra |
|---|---|---|
| 0.00 (nulo) | 2.6e-07 | **−0.1941** |
| 0.20 | 3.0e-06 | −0.1661 |
| 0.50 | 2.4e-05 | −0.1319 |
| 1.00 | **6.9e-05** | −0.1018 |

La identidad es **monótona por encima de su propio suelo** y separa señal de nulo por **263×**.
Las ventanas quedan **dominadas por el sesgo `−k/n`**: las cuatro salen negativas y su recorrido
(0.092) es menor que el desplazamiento (0.194). Es la tesis de la adenda, medida.

### ⚠ Desviación declarada sobre el nulo del §C.4.1

El §C.4.1 pide «barajados de `ε` preservando el precio». Un barajado destruye también la
**memoria larga de `ε`**, que es real (`N_eff` de los signos = 46.5) y que infla la varianza del
estimador por el solapamiento de las parejas — **el modo de fallo 8 de la propia adenda**. Se usa
**rotación circular**, que conserva toda la autocorrelación de `ε` y rompe sólo el
emparejamiento, y se reportan los dos. Sobre el control sintético con `ε` de memoria, la rotación
da un suelo **48× más ancho** (0.000432 contra 0.000009): la diferencia no es cosmética.

### ⚠ Tres defectos propios, los tres de familias que este proyecto ya conoce

1. **Alineación de rezagos, en mi propio control.** Inyecté la señal como `a·cumsum(ε)`, pero
   `cumsum[j] − cumsum[i] = Σ_{s=i+1..j} ε_s`, que **excluye `ε_i`**: la «señal» era predecible
   desde el futuro, no desde el presente, y la identidad devolvía 0 **correctamente**. Misma
   familia que el `mode="same"` de la v3.1 §3 y el `h(1) = 0` del núcleo.
2. **Dos aserciones mal formuladas.** Exigía monotonía incluso donde la señal inyectada cae por
   debajo del suelo de muestreo, y afirmaba que las ventanas «no son monótonas» — medido, **sí**
   lo son, pero eso no las salva: lo cierto es peor, el sesgo las domina.
3. **Memoria.** Cargar los cinco fragmentos a la vez son ~1.1 GB con ~1 GB libre, y el proceso
   murió. Convertido en generador. **Tercera vez que la RAM de este portátil decide la
   arquitectura del análisis.**

### Orden codificado, no confiado a la memoria

`--etapa=medir` **se niega a correr** si no existe el archivo de suelo congelado. El §C.4.1 exige
producir la tabla de resolución antes de calcular un solo `R²` real, y un orden que depende de
que alguien se acuerde no es un orden.

---

## Sesion 2026-08-23 -- v4.1 §1 rehecho, y RETRACTADO el mismo dia por estacionalidad

⚠⚠ **EL VEREDICTO DE ESTA SESION NO SE SOSTIENE. NO CITARLO.** Lo que sigue se conserva entero
porque los datos y los controles valen; **la lectura no**. Objecion del operador, verificada con
medicion el mismo dia:

1. **Los tramos estan definidos por cuando se cayo la conexion, no por calendario**, y el
   confundido es severo: el tramo 1 es **100 % dias habiles**, el tramo 2 es **100 % fin de
   semana**, el 0 mezcla ambos. Etiquetarlos por `ν` y llamar a eso «regimen medido» es
   incorrecto: buena parte de esa `ν` es dia de la semana.
2. **`σ₁` se estimo agrupando las 24 horas**, y `R²_req ∝ (lastre/σ)²`. Medido sobre las cuatro
   muestras (16 340 ventanas de 60 s):

   | | `σ₆₀` [pb] | **`R²_req(300 s)`** |
   |---|---|---|
   | agrupando todas las horas (lo publicado) | 5.301 | **5.52 %** |
   | **mejor hora, UTC 15 (mediodia EE.UU.)** | **9.536** | **1.71 %** |
   | peor hora, UTC 23 | 3.289 | 14.34 % |
   | habil (lun-vie) | 5.937 | 4.40 % |
   | fin de semana | 3.418 | 13.28 % |

   **El requisito varia 8.4× segun la hora.** El maximo de `σ` cae en UTC 15, que es el mediodia
   de EE.UU.: es estructura diurna, no ruido.
3. **`ν` entra en el predictor como ESCALAR.** `rasgos_flujo` construye las ventanas
   retrospectivas como `f·H·ν` **ticks** con la `ν` media del tramo. Medida por hora, `ν` recorre
   **0.55× a 2.17×** de su media (factor **3.96**, max UTC 15, min UTC 4): en la hora punta la
   ventana cubre el **46 %** de los `H` segundos que dice cubrir. No es solo una descripcion
   pobre de `ν` -- esta dentro del estimador.

**Por que esto no es un matiz.** Los dos lados de la desigualdad se mueven al condicionar por
hora, y solo se ha medido que se mueve uno. Si la relacion predictiva es ella misma estacional,
agrupar la **atenua hacia cero**: una senal viva en una fracción de las horas aparece diluida por
su cuota de varianza. El `+0.0076` medido a 60 s con 1 341 ventanas y controles limpios es
compatible con un `R²` de un dígito alto concentrado en las horas activas -- y el requisito en
esas horas baja a **1.71 %**. **La medicion agrupada no puede descartar el cruce**, que es
justo lo que el veredicto afirmaba.

**Lo que SI sobrevive de la sesion:** el alineado libro-transaccion (100 % en los cuatro tramos,
spread de 1 tick, `E[y_mid·ε] > 0` en los cuatro), el control positivo de potencia y su
maquinaria, y la banda de `H_p` **ajustada en ventanas cortas** (`[10,100]`, `[30,300]`,
`[60,600]` s), donde el ciclo diurno de 24 h no alcanza a inclinar la pendiente. Los ajustes en
`[60, 3600]` y `[300, 14400]` **si** estan contaminados por el ciclo y no se leen.

**Lo que hay que rehacer, y en este orden:** (a) ventanas retrospectivas definidas en
**segundos**, no en ticks via `ν` escalar; (b) `σ₁` y por tanto `R²_req` **por casilla horaria**;
(c) el `R²` medido **estratificado** por la misma casilla, con su control barajado y su control
de potencia dentro de cada estrato; (d) recomponer el veredicto por estrato, no agrupando.


### Medicion de la estacionalidad (mismo dia) -- modulo `estacionalidad.py`, 8/8 controles

Modelo pedido por el operador: armonicos de 24 h + aperturas de sesion como **pulsos** con
respuesta exponencial + fin de semana con **amplitud y fase propias**. Ajuste lineal (un par
cos/sin con coeficientes libres YA es amplitud y fase libres, asi que no hay minimos locales --
que es lo que hundio el estadistico de `omega_G` en la v3.2). Unico parametro no lineal: `tau`,
barrido en rejilla y reportado como curva.

⚠ **Esto NO resucita `omega_m`.** Aquella era una frecuencia **endogena** sin ancla y sin nulo
(la EMD devolvia 118 s sobre un paseo aleatorio y el periodo escalaba con la ventana). Aqui el
periodo **no se estima**: es 24 h y 168 h, conocido a priori, con causa exogena verificable. Se
estiman amplitud y fase, con nulo propio (barajar hora y finde) y validacion **por dias enteros**.

**Perfil descriptivo** (3 199 casillas de 5 min, 13 dias UTC):

| | pico | valle | razon |
|---|---|---|---|
| habil | **UTC 13-15** (apertura NY): `ν` 61.8 tx/s, `σ` 3.92 pb/30 s | UTC 4: 15.9, 1.43 | `ν` 3.9× · `σ` 2.7× |
| finde | **UTC 21-22**: `ν` 20.4, `σ` 1.95 | UTC 3: 8.2, 0.45 | — |

habil/finde: `ν` **2.13×**, `σ` **1.98×**. El pico de fin de semana llega **~7 h mas tarde**.

⚠ **Fuera de muestra por dias enteros, el modelo sobre el NIVEL es peor que una constante**
(`log σ`: M1 −0.013, M3 −0.179; `log ν`: M1 +0.008 contra nulo −0.008). Con 13 dias, el **nivel
del dia** domina y no se predice desde la hora. Hay que separar nivel y forma:

| | var entre dias | var dentro del dia | recorrido del nivel diario |
|---|---|---|---|
| `log σ` | 46 % | 54 % | **13.75×** |
| `log ν` | 51 % | 49 % | **20.78×** |

**Forma diurna con el nivel del dia retirado a ambos lados, `R²` fuera de muestra:**

| modelo | `log σ` | `log ν` |
|---|---|---|
| M1 armonicos 24 h | +0.0294 | +0.1205 |
| M3 + amplitud y **fase** de finde | +0.0073 | **+0.1359** |
| M4 + **pulsos**, `tau` = 2 h | −0.0033 | **+0.1394** |
| NULO (hora y finde barajados) | −0.0110 | −0.0103 |

**Conclusiones, y una corrige una cifra mia de esta misma sesion:**

1. **`ν` escalar queda refutado con dato.** Su forma diurna replica fuera de muestra a **13× el
   nulo**, y el modelo completo del operador —pulsos con decaimiento + desfase de finde— es el
   mejor de la escalera. `rasgos_flujo` debe definir sus ventanas en **segundos**, no en
   `f·H·ν` ticks.
2. ⚠ **El «el requisito varia 8.4× segun la hora» que se midio antes esta INFLADO.** Salia de
   medias horarias agrupadas, con ~13 observaciones por hora y dias desbalanceados entre horas.
   La amplitud diurna de `σ` que replica fuera de muestra corresponde a **~1.13×**, no a 2.7×.
   Anadir finde y pulsos a `σ` **empeora** el ajuste: con ~4 dias de fin de semana, sobreajusta.
3. **La palanca grande no es la hora, es el DIA.** El nivel diario de `σ` recorre **13.75×** y
   `R²_req ∝ σ⁻²`, o sea **~190×** de recorrido en el requisito entre dias. Y a diferencia de
   `omega_m`, la volatilidad diaria es persistente y pronosticable por vias establecidas
   (HAR, GARCH), asi que **esa** es la estratificacion con contenido.
4. **13 dias son pocos** para el perfil semanal, y la propia captura lo avisa en su informe de
   cobertura: «con 14 dias hay 2 observaciones por casilla hora x dia; el perfil semanal es
   EXPLORATORIO hasta >= 4 semanas». Sigue corriendo hasta el 2026-09-02.

⚠ **Un control mio fallo, y era un SIGNO — la cuarta vez en este proyecto** (tras el 2π de la
v1.3, el factor 125 y la convencion de `ε` del propagador). Con `cos(2π(h−φ)/24)` la fase es
`atan2(c_sin, c_cos)` y estaba escrito `atan2(−c_sin, c_cos)`, que devuelve `24 − φ`. Ademas lo
contrastaba contra una verdad que mezclaba habil, finde y pulso. Corregido: recupera **15.00 h y
20.00 h exactas** sobre verdad conocida.

---

### H1 y H2 (2026-08-23) -- `hipotesis_liquidez.py`, 7/7 controles

Predicciones **congeladas en el docstring del modulo antes de calcular nada** (P1-P6).

#### ⚠ H2 CONFIRMADA, y con una desviacion cuantitativa que importa

`ν` = numero de perturbaciones, `σ` = respuesta. Sobre 3 199 casillas de 5 min, 13 dias:

| relacion `log σ` contra `log ν` | pendiente | corr |
|---|---|---|
| todo | +0.8557 | +0.8775 |
| **entre dias** | **+0.8830** | **+0.9741** |
| dentro del dia | +0.8045 | +0.7866 |

**P4 REFUTADA en su forma literal: la pendiente NO es 0.5, es ~0.8.** Y no es artefacto de
discretizacion — barriendo el paso de submuestreo, con el tamano tipico del movimiento pasando de
0.68 a 2.55 ticks, la pendiente va de 0.881 a **0.786** y ahi converge:

| paso [s] | 10 | 30 | 60 | 120 | 300 |
|---|---|---|---|---|---|
| ticks por paso | 0.68 | 1.21 | 1.67 | 2.13 | 2.55 |
| pendiente | 0.881 | 0.856 | 0.839 | 0.801 | **0.786** |

Pendiente > 0.5 significa que **las perturbaciones NO son independientes**: se agrupan y se
refuerzan. Coherente con la memoria larga del flujo de ordenes que este proyecto ya midio.

⚠ **NO CONFUNDIR ESE 0.8 CON `H_p`.** Es un exponente **transversal** entre ventanas
(`σ` de una ventana contra `ν` de esa ventana), no el exponente de escala temporal de la firma de
volatilidad, que sigue en la banda `[0.372, 0.554]`. Es exactamente la confusion que costo el §1
de la v3.3, donde `γ` significaba dos cosas.

**P5 CONFIRMADA:** normalizar por `ν` reduce la varianza de `log(varianza)` un **63.7 %**, y
**entre dias un 77 %** (1.8895 -> 0.4340). `ν` contiene la mayor parte de `σ`, sobre todo su nivel
diario. No la contiene entera: queda un 36 % sin explicar, casi todo dentro del dia.

**P6 CONFIRMADA:** tras quitar `log ν`, al residuo no le queda estructura diurna (`R²` fuera de
muestra del modelo diurno: −0.039, contra −0.179 sobre `log σ` crudo).

**Y lo que decide en la practica** -- `R²` fuera de muestra por dias enteros, prediciendo la `σ`
de la casilla **siguiente**:

| predictor | `R²` fuera de muestra |
|---|---|
| constante | +0.0000 |
| **`log ν(t)`** | **+0.5433** |
| `log σ(t)` | +0.5101 |
| `log ν(t)` + `log σ(t)` | +0.5575 |
| `log ν(t)` barajado dentro del dia (NULO) | +0.1611 |

**`ν` predice la volatilidad futura MEJOR que la propia volatilidad presente.** Ese es el
contenido operativo de H2: da un estado observable, adelantado y pronosticable.

⚠ **Esto tiene nombre y se declaro antes de medirlo:** subordinacion del precio al reloj de
transacciones, Clark (1973); Hipotesis de Mezcla de Distribuciones, Tauchen & Pitts (1983). El
proyecto ya tropezo con la version de 1994 al medir `q̄`. **Confirmarla valida la tuberia y da
un modelo con el que trabajar; no es un descubrimiento.**

#### H1 -- la profundidad de nivel 1 NO aporta sobre el flujo

| | `b` (`log ν`) | `c` (`log D`) | parcial(`σ`, `D` \| `ν`) | nulo IC95 |
|---|---|---|---|---|
| **dentro del dia** (n = 1 622) | +0.8209 | **+0.0649** | **+0.0425** | [−0.0438, +0.0418] |
| entre dias (n = 13) | +0.8337 | −0.1358 | −0.3000 | \|r\| critico 0.632 |

**Dentro del dia gana P2 (la intuicion del operador, `c > 0`) y no P1 (Kyle, `c < 0`)** — pero por
un pelo: +0.0425 contra un techo de nulo de +0.0418, o sea **0.18 % de la varianza residual**.
Entre dias el signo se invierte y con 13 dias no hay potencia para nada.

**P3 CONFIRMADA y es la lectura principal: `|c|` es despreciable frente a `b`.** Sabiendo el flujo,
la cola del mejor precio no anade practicamente nada.

Dato colateral que si es solido: **`corr(log D, log ν) = −0.485`** (−0.255 dentro del dia, −0.658
entre dias). **El libro esta MAS FINO cuando hay mas actividad**, no mas saturado.

Perfil de profundidad: en dias habiles es plano (22–31 BTC); en fin de semana **se desploma a
12–15 BTC entre las 12 y las 20 UTC** contra 22–27 el resto. Con 4 dias de fin de semana, es
exploratorio.

⚠ **Limitacion que no se arregla con estos datos:** `@bookTicker` da solo el **nivel 1**. La
«saturacion del libro» con ordenes en reposo a varios niveles **no es observable aqui**. Para
contrastar H1 de verdad haria falta capturar `@depth`.

⚠ **Un nulo mio estaba mal especificado y lo delato el signo.** La primera version barajaba la
profundidad dentro del dia sobre las series SIN desmediar, y su IC95 salia
`[−0.1065, −0.0554]`: descentrado, y con el valor real (−0.0456) **fuera del nulo por el lado
contrario**. La causa: barajar dentro del dia **conserva intacto el termino entre dias**, asi que
el nulo no destruia lo que yo creia. Desmediando por dia, el nulo queda centrado en **+0.0004**.
Es la misma leccion que el control de potencia del §1 de esta misma sesion, y la tercera vez en
el dia: **hay que separar nivel diario y forma intradia antes de construir cualquier nulo.**

---

### Retroalimentacion `σ`-`ν` y prediccion por etapas (2026-08-23) -- `retroalimentacion.py`, 8/8

#### A -- El VAR: la retroalimentacion es ASIMETRICA

`x_t = [log ν_t, log σ_t]`, 3 rezagos de 5 min, 3 081 observaciones contiguas, 13 dias.

```
log ν(t)      <-  +0.6053 * log ν(t-1)   -0.0543 * log σ(t-1)
log σ(t)      <-  +0.2428 * log ν(t-1)   +0.2025 * log σ(t-1)

autovalores 0.5693 y 0.2384   ->  vida media 1.23 casillas = 6.2 min
```

**`ν` empuja a `σ` (+0.243); `σ` casi no empuja a `ν` (−0.054), y encima con signo
negativo** — mas volatilidad ahora va con algo MENOS actividad despues. `ν` es ademas
mucho mas persistente (0.605 contra 0.203): la actividad tiene memoria, la volatilidad
casi no.

**Direccion, fuera de muestra por dias enteros** (que es lo que decide, no un test F
dentro de muestra):

| objetivo | solo rezagos de `ν` | solo rezagos de `σ` | ambos |
|---|---|---|---|
| `log ν(t)` | **+0.7719** | +0.6488 | +0.7721 |
| `log σ(t)` | +0.5864 | +0.6039 | **+0.6153** |

**Anadir `σ` a la prediccion de `ν` aporta +0.0002: nada.** Anadir `ν` a la de `σ` aporta
+0.011, y `σ` aporta +0.029 sobre `ν` sola. O sea: **no es un lazo simetrico, es un motor
(`ν`) y una respuesta (`σ`)** — pero para PREDECIR `σ` conviene llevar las dos, porque son
~95 % redundantes y el 5 % restante es real.

⚠ Matiz frente al hallazgo de H2: con UN rezago `ν(t)` batia a `σ(t)` (0.543 contra 0.510);
con TRES, la historia propia de `σ` adelanta a la de `ν` (0.604 contra 0.586). No es
contradiccion: los rezagos de `σ` acaban conteniendo la informacion de `ν` de forma
indirecta.

**Respuesta al impulso** (choque unitario propagado con `A_1`; es una aproximacion, el
modelo tiene 3 rezagos):

| pasos (5 min c/u) | 1 | 2 | 3 | 4 | 6 | 8 | 12 |
|---|---|---|---|---|---|---|---|
| choque en `ν` → `σ` | **+0.243** | +0.196 | +0.125 | +0.075 | +0.025 | +0.008 | +0.001 |
| choque en `σ` → `ν` | −0.054 | −0.044 | −0.028 | −0.017 | −0.006 | −0.002 | −0.000 |

Un choque de actividad se agota en `σ` en **~1 hora**.

#### B -- Prediccion por etapas: la barra de error VIVA

⚠ **Precision sobre la propuesta, porque sin ella el test es vacuo.** "Predecir `σ_t` desde
`t−1`" y "predecir `σ_{t+1}` desde `t`" son **la misma operacion desplazada un paso**: bajo
estacionariedad sus errores coinciden por construccion. Medido: RMSE 0.5446 contra 0.5631,
**razon 1.034**, y errores consecutivos con correlacion **+0.0018** (o sea, el error de un
paso es blanco: no queda estructura que exprimir en el nivel). Eso **no es un hallazgo**.

Lo que si tiene contenido es la version util de la misma idea, y ahi funciona:

**(1) El error reciente PREDICE el error siguiente, y mejor cuanto mas larga la ventana:**

| ventana | 15 min | 30 min | 60 min | 120 min |
|---|---|---|---|---|
| n | 1 005 | 492 | 234 | 106 |
| corr(log RMS pasado, log RMS futuro) | +0.2029 | +0.3219 | +0.3961 | **+0.5940** |
| `R²` fuera de muestra | +0.033 | +0.088 | +0.117 | **+0.292** |

**Hay barra de error viva y autocalibrada a escala de 1-2 h.** Es exactamente lo que el
operador buscaba con las dos etapas, formulado de modo que no sea tautologico.

**(2) NIS -- ?es honesta la incertidumbre declarada?**

```
var(z) con sd CONSTANTE  = 1.0000   (por construccion; no informa)
var(z) con sd PREDICHA   = 1.3089   (1.0 = honesta)
curtosis: 12.71 con sd constante  ->  8.77 con sd predicha
fraccion |z| > 3 : 1.266 %  ->  1.915 %   (normal: 0.270 %)
```

La barra viva **acerca la varianza a 1 y reduce la curtosis**, pero **empeora la cola
lejana**: en los tramos tranquilos la `sd` predicha es demasiado pequena y una sorpresa da
un `z` enorme. Sigue siendo sobreconfiada, igual que el NIS del EAKF en la Fase 1.

**(3) `ν` predice el TAMANO del error, debilmente:** `R²` fuera de muestra sobre
`log|error|` = +0.0140 con `log ν(t−1)` y +0.0272 con los rezagos completos.

⚠ **Defecto propio corregido, y sin corregirlo el veredicto del (2) era otro.** La `sd`
predicha se estima modelando `log|e|`, y hay que deshacer el sesgo de Jensen: para
`e ~ N(0,s)`, `E[log|e|] = log s − 0.63518`, asi que `s = exp(E[log|e|] + 0.63518)`. Yo usaba
`exp(E[log|e|])·√(π/2)` — que es la conversion valida para `E|e|`, **no** para
`exp(E[log|e|])`. El error era de **1.506× en `s`, o sea 2.27× en `var(z)`**: daba 2.97 y
habria declarado deshonesta una barra que esta en 1.31.

#### Lo que esto entrega para el §1

`σ` a un paso se predice con `R²` fuera de muestra **+0.6153** y **error relativo mediano
del 30 %**. Como `R²_req ∝ σ⁻²`, un 30 % en `σ` son **~60 % en el requisito**: suficiente
para una **compuerta gruesa** de operar / no operar, insuficiente para dimensionar posicion
con ella. Y con la barra viva del punto (1), esa compuerta puede llevar su propia confianza
adjunta en vez de un umbral fijo.

---

### Registro de lo ejecutado (lectura retractada, cifras validas)


Rehace el §1 de `ORDEN_TRABAJO_HORIZONTE_4_1.md` sobre `captura_estacional`, que el 2026-08-22
pasó su compuerta. `Micelio.py` sin cambios. Módulo nuevo: `curvas_estacional.py`.
Acta completa en `telemetria/acta_v41_sec1_estacional.txt`.

### ⚠ VEREDICTO: segundo desenlace del §1.4 — NO HAY BANDA VIABLE

```
filas con n_val >= 30 (decidibles)              : 15 de 26
filas decidibles que CRUZAN con control limpio  :  0
mejor razon R2_medido / R2_req                  : 0.0758   (tramo 2, H = 300 s)
```

El §1.1 declaró este desenlace admisible **antes de medir**: «este documento no busca rescatar el
proyecto; busca decidirlo».

### Cuatro tramos, y el régimen pasa a ser una variable medida

| tramo | horas | ticks | `ν` | precio | recorrido | `σ₁` difusivo [pb·s^−½] |
|---|---|---|---|---|---|---|
| 0 | 131.83 | 8 898 311 | 18.75 | 62 484 – 64 601 | 3.4 % | 0.3827 |
| **1** | **82.47** | **19 985 140** | **67.31** | **63 979 – 79 555** | **24.3 %** | **1.0622** |
| 2 | 34.72 | 4 099 830 | 32.80 | 75 588 – 78 058 | 3.3 % | 0.6333 |
| 3 | 23.33 | 1 829 242 | 21.78 | 63 212 – 64 470 | 2.0 % | 0.3862 |

Alineado libro-transacción **100 % en los cuatro**, spread mediano **1 tick exacto** en los cuatro,
y `E[y_mid·ε] > 0` en los cuatro — la convención `m = True → ε = −1` se sostiene en cada réplica.

### §2 — La banda de `H_p` REPLICA, y la superdifusión no

24 ajustes (4 tramos × 6 ventanas de escala), reloj de pared, control barajado en todos:

```
banda medida ahora : [0.3719, 0.5535]   mediana 0.4932   |sesgo| del control <= 0.0438
banda publicada    : [0.371 , 0.552 ]
```

Reproduce la banda del 2026-08-10 con tres decimales, y la **mediana es 0.493: difusivo**. El
control barajado vuelve a 0.50 en las 24 filas, así que el estimador está limpio a estas
longitudes. **`H_p = 0.591` de la v3.3 queda fuera de la banda entera** — la retractación se
confirma con 24 ajustes en vez de 4, y **la superdifusión no replica en ningún tramo**.

### §1 — Las dos curvas. Tabla por tramo

Predictor tonto del §1.3 (flujo firmado acumulado en `{H/4, H/2, H, 2H}`, sin `G(τ)`), ventanas
no solapadas, partición 60/20/20 **en tiempo** con embargo de 4 h, y el último 20 % de cada tramo
**sin abrir**. El conjunto de prueba de la v3.2 tampoco se abrió.

| tramo | `H` | n_ent | n_val | **R²_medido** | bar_global | bar_ventana | **R²_req mín** | potencia |
|---|---|---|---|---|---|---|---|---|
| 0 | 60 s | 4 745 | 1 341 | **+0.0076** | −0.0052 | −0.0100 | 77.38 % | 0.99 |
| 0 | 120 s | 2 372 | 670 | +0.0041 | −0.0484 | −0.0179 | 43.31 % | 0.99 |
| 0 | 300 s | 949 | 267 | −0.0227 | +0.0132 | −0.0470 | 15.75 % | 0.98 |
| 0 | 600 s | 474 | 133 | −0.0777 | −0.0226 | −0.0914 | 7.33 % | 0.58 |
| **1** | **60 s** | 2 969 | **748** | **−0.0032** | −0.0021 | −0.0009 | **10.05 %** | 0.92 |
| **1** | **120 s** | 1 484 | **374** | **−0.0058** | −0.0028 | +0.0038 | **5.62 %** | 0.81 |
| **1** | **300 s** | 593 | **149** | **−0.0001** | +0.0014 | +0.0073 | **2.04 %** | **1.18** |
| 2 | 60 s | 1 249 | 176 | −0.0046 | −0.0029 | −0.0099 | 28.26 % | 0.63 |
| 2 | 120 s | 624 | 88 | −0.0079 | −0.0160 | −0.0085 | 15.82 % | 0.73 |
| 2 | 300 s | 249 | 35 | +0.0044 | −0.0203 | −0.0195 | **5.75 %** | 0.91 |
| 3 | 60 s | 840 | 39 | −0.0187 | −0.0098 | +0.0091 | 76.00 % | 0.86 |

**La fila que decide es `tramo 1, H = 300 s`**, y hay que leerla entera: es el régimen **más
favorable** de los cuatro —24 % de recorrido de precio, `σ₁` 2.8× la del tramo 3, así que el
requisito se desploma al **2.04 %**—, tiene **149 ventanas** no solapadas, el control positivo
demuestra que una señal del 2 % **se recupera** ahí (devuelve 2.4 %), y el `R²` medido fuera de
muestra es **−0.0001**. Cero exacto donde el peaje era más barato que nunca.

### ⚠ El control positivo de potencia — sin él «R² = 0» no habría sido un resultado

La cuarta fila del §1.4 manda no leer una `H` donde el barajado es comparable al real. Aquí **los
dos salen ≈ 0**, así que la regla no discrimina entre «no hay señal» y «no hay potencia» — y sólo
una de las dos decide el proyecto. Se inyecta en el objetivo una señal del tamaño **exactamente
igual al `R²` requerido** y se mide con el mismo procedimiento:

| n_val | ≥ 267 | 133–176 | 88 | ≤ 49 |
|---|---|---|---|---|
| razón recuperado/inyectado | **0.98–1.18** | 0.58–0.63 | 0.73 | negativa |

**11 de las 15 filas decidibles tienen potencia demostrada**, y son las que sostienen el veredicto.
Las 4 que no (H ≥ 600 s en los tramos 0 y 1) se marcan y **no se leen**, aunque su `R²` medido sea
negativo y «favorezca» la conclusión.

⚠ **Tres versiones del control fueron mías y estaban mal**, y el patrón es el mismo de siempre:

1. Normalizando con estadísticos de **entrenamiento**, la señal inyectada llegaba encogida a
   validación —los rasgos de flujo no son estacionarios en escala— y el control declaraba «sin
   potencia» con 1 249 puntos de ajuste.
2. Normalizando **globalmente**, seguía contaminado: la varianza del retorno real cambia hasta
   **4.6×** entre bloques (columna `sd(y)val/ent`: 0.38 en el tramo 0, 1.85 en el 2), así que el
   control medía agrupamiento de volatilidad y no tamaño muestral.
3. Con **una sola dirección** de señal el resultado saltaba de 0.32 a 1.6 entre celdas vecinas: los
   cuatro rasgos son sumas acumuladas anidadas y hay direcciones casi degeneradas. Se promedia
   sobre **25 sorteos** y la razón pasa a ser monótona en `n`, que es lo que debe ser.

La versión correcta normaliza **dentro de cada bloque** —la pregunta es si una fracción `R²_req`
de la varianza *de validación* sería visible— y promedia sobre direcciones.

### Lo que esto cierra, y lo que deja abierto

**Cierra.** La hipótesis del horizonte era, según el §1.1, **la única palanca con el orden de
magnitud correcto** para el factor 3.46 que dejó parada la v3.2. Medida en cuatro regímenes de `ν`
que van de 18.8 a 67.3 tx/s, con hasta 1 341 ventanas y con potencia acreditada, **no cruza en
ninguno**. El mejor margen es **0.076**, o sea un factor **13** de defecto, y el §1.2 dice que
refinar el predictor sólo tiene sentido «si la curva cruza con el predictor tonto o queda cerca».
No queda cerca.

**No cierra:**
- El `+0.0076` del tramo 0 a 60 s (1 341 ventanas, controles en −0.005 y −0.010) es **señal real y
  fuera de muestra**. Es el mismo cuadro que el paso 2 de la v3.2: **hay estructura medible y es
  ~100× menor que el peaje**. Lo que se refuta no es que exista señal, es que pague.
- El §1 se midió con el predictor que el propio documento impone. Un predictor mejor no está
  descartado *en principio*; está descartado *por el criterio del documento*, que es distinto y hay
  que decirlo así.
- Nada de esto toca `H > 4 h`: a 4 h el tramo más largo da 5 ventanas de validación. La banda de
  decenas de minutos sí queda cubierta, y es la que importaba.

### Reservas declaradas

1. **`c(u)` sigue usando tarifas asumidas VIP 0.** Mismo criterio no cumplido que en la v3.1 y la
   v3.2: `/fapi/v1/commissionRate` es firmado y el Modo LECTURA no tiene credenciales.
2. **La partición es 60/20/20 en tiempo**, así que entrenamiento y validación son tramos
   *distintos* de mercado. Con la volatilidad cambiando hasta 4.6× entre bloques, eso penaliza al
   `R²` fuera de muestra — y es deliberado: es la situación real de operar.
3. **Cuatro tramos no son cuatro muestras independientes** del mercado: son cuatro trozos de doce
   días consecutivos de BTCUSDT.
4. `curvas_estacional.py` **no tiene suite propia**; importa las funciones de `horizonte.py`, que
   sigue en **14/14**, así que la aritmética de las dos curvas es la ya verificada. Lo nuevo y no
   cubierto por controles es la carga por tramos, que sí se verificó por sus invariantes
   (cobertura 100 %, spread de 1 tick, `G(0) > 0` en los cuatro).

---

### `c(u)` real (2026-08-23) -- `coste.py`, 13/13 controles

**`c(u)` es el COSTE DE TRANSACCION de ida y vuelta por unidad, en pb**, y la `u` esta ahi
porque depende de la accion (maker o taker). **No** es el coeficiente difusivo de Loeper
-- ese es `sigma^2`, la misma `sigma` que ahora se pronostica -- ni `c2_vol = k*omega_m*nu`.
Tres `c` distintas, y `coste.py` pasa a ser la fuente unica.

**Descomposicion, con procedencia por componente:**

| componente | procedencia | valor |
|---|---|---|
| comision | endpoint **FIRMADO** | ⛔ **ASUMIDA** VIP 0. **81.8 %** de `c(u)` |
| cruce de spread | la captura, 34.8 M ticks | ✅ **MEDIDO: 0.0146 pb** |
| seleccion adversa | v4.0 §5 | ✅ medida: 0.888 pb |
| financiacion | endpoint **PUBLICO** | ✅ **LEIDA: 0.4261 pb / 8 h** |

**Lo que se cerro hoy:**

1. **La financiacion se LEE**, ya no se asume. 500 periodos, 2026-03-10 → 08-24:
   `|tasa|` mediana **0.4261 pb/8 h** (p90 0.8928, ultimos 30 dias 0.6147), media con signo
   **+0.2333** — y con signo importa, porque **es una transferencia, no un coste: un corto
   la COBRA**. La constante asumida era 0.9681 pb/8 h, o sea **2.27x demasiado alta**. El
   efecto sobre el §1 es nulo por debajo de 1 h y < 2.5 % a 2 h, asi que no cambia ningun
   veredicto — pero deja de ser un numero inventado.
2. ⚠ **Defecto corregido en `horizonte.py`**: alli la financiacion es
   `max(H_s − 3600, 0)/(8·3600)`, o sea **cero por debajo de una hora**, sin justificacion.
   Una tenencia de `H` cruza una marca de financiacion con probabilidad `H/(8 h)`, asi que
   su coste esperado es **lineal en `H` desde 0**. Numericamente irrelevante (0.0044 pb a
   300 s contra 4 pb de comision) pero corregido en `coste.py`.
3. **El spread es despreciable, y por mucho mas de lo que se creia.** Mediana **0.0146 pb**
   sobre 34.8 M de ticks alineados, contra 4.00 pb de comision maker ida y vuelta: **274x**.
   Y la diferencia taker−maker (6 pb) es **411x** el spread. **Cruzar el spread no es lo que
   cuesta; cuesta el escalon de comision.**
4. ⚠ **Inconsistencia del repo, encontrada al unificar.** `propagador.py` y `cola.py` usan
   taker = `0.0005` (correcto para futuros USD-M VIP 0); `horizonte.py` **ignora la taker por
   completo** y su `LASTRE_IDA_VUELTA_PB` sólo lleva maker.

**Y eso ultimo importa mas que todo lo demas de este apartado:**

| esquema | `c(u)` a H = 300 s | `R²_req` × | tramo 1 a 300 s |
|---|---|---|---|
| **maker + maker** (lo que el §1 uso) | 4.8926 pb | ×1.00 | 2.04 % |
| maker + taker | 7.4558 pb | ×2.32 | 4.74 % |
| **taker + taker** | **10.0190 pb** | **×4.19** | **8.55 %** |

**Todo el §1 se midio suponiendo llenado maker perfecto en las dos patas — el mejor caso
posible.** Con ejecucion taker el requisito se multiplica por 4.2.

**Lo que NO se puede cerrar, y por que:**

- ⛔ **El criterio del §8 sigue incumplido.** `/fapi/v1/commissionRate` es firmado y el
  **81.8 %** de `c(u)` descansa en ese numero. `coste.py` trae el lector firmado listo
  (`leer_comision_firmado`) y una bandera `COMISIONES_LEIDAS` que **no se puede poner a
  cierto sin pasar tarifas explicitas**, con test que lo comprueba.
- ⚠ **NO valen las credenciales de Testnet.** El escalon es propiedad de la CUENTA (nivel
  VIP, descuento BNB, referido) y la de Testnet no es la de Mainnet. La v1.3 ya midio que el
  `stepSize` de Testnet es 10x mas fino y dejo escrito que calibrar contra el entorno
  equivocado produce un sistema que funciona en pruebas y se degrada en produccion.
- **Basta una clave de MAINNET de SOLO LECTURA.** Este endpoint no necesita permiso de
  trading ni de retiro, y pedir mas permisos de los necesarios es el riesgo que no hay que
  correr.
- ⚠ **La tercera pata: `c(u)` no es constante.** Con llenado maker incierto,
  `c_efectivo = p·c_maker + (1−p)·C_respaldo`. `cola.py` ya tiene la estructura; `p` y
  `C_respaldo` son el §5 de la v4.1, pendiente y a su vez a la espera de que el §1 diga a
  que horizonte.

#### Lectura firmada de TESTNET (2026-08-23) -- verifica el codigo, NO cierra el criterio

El operador paso las credenciales de la cuenta **demo**. Se advirtio antes de usarlas que
Testnet no cierra el §8 y se uso igual, con la procedencia marcada. **Las claves no se
escriben en ningun archivo del repo ni se imprimen**; el lector las toma en memoria.

| endpoint firmado | resultado [TESTNET] | contra lo asumido |
|---|---|---|
| `/fapi/v1/commissionRate` maker | **0.000200** | **COINCIDE exacto** |
| `/fapi/v1/commissionRate` taker | **0.000400** | ⚠ **el repo asume 0.000500: 25 % mas alta** |
| `/fapi/v2/account` | `feeTier = 0` | confirma VIP 0 |
| `/fapi/v1/leverageBracket` tramo 1 | **mmr = 0.0040** (nocional ≤ 50 000) | **COINCIDE exacto** |

**Tres cosas que esto si establece:**

1. **El codigo de firma funciona de punta a punta.** `leer_comision_firmado` esta verificado
   contra un servidor real; cuando haya credenciales de Mainnet es correr y ya.
2. ⚠ **La comision taker del repo esta desactualizada en 4 modulos.** `propagador.py`,
   `cola.py`, `tick_grande.py` y `coste.py` asumen `0.0005`; la lectura da `0.0004`
   -- Binance bajo la taker de futuros de 0.0500 % a 0.0400 % en algun momento y el proyecto
   no se entero. **Se conserva 0.0005 por omision porque es la CONSERVADORA** (mas coste =
   requisito mas duro = conclusion negativa mas robusta) y se expone la leida al lado. Con la
   leida, `c(u)` taker+taker baja de 10.02 a **8.02 pb** y `R²_req` se multiplica por **0.64**.
3. ✅ **El hueco de `mmr` de la v1.3 queda cerrado.** Aquella sesion dejo escrito que
   `leverageBracket` es firmado, que `mercado.leer_mmr` devuelve el valor asumido **inflado
   por un factor de seguridad de 2x** y que la guarda queda conservadora ante la duda. La
   lectura confirma `mmr = 0.0040` para el primer tramo, **exactamente el valor asumido**: el
   factor 2x era conservadurismo, no ignorancia.

⚠ **Lo que NO establece, y por que se mantiene el criterio como incumplido:** el escalon es
propiedad de la **cuenta de Mainnet**. Una cuenta de Testnet nace en `feeTier = 0` por
construccion y no sabe nada del descuento BNB (−10 %), del nivel VIP real ni de un referido.
`COMISIONES_LEIDAS` sigue en `False` y hay test que lo comprueba. **Basta una clave de
Mainnet de SOLO LECTURA** -- este endpoint no necesita permiso de trading ni de retiro.

#### Lectura de MAINNET (2026-08-23): el criterio del §8 pasa a PARCIAL

El operador paso tambien las credenciales de Mainnet. **Las claves no se escriben en ningun
archivo del repo ni se imprimen**, y se verifico con `git grep` que no quedan en el arbol ni
en el indice -- el `origin` es publico.

⚠ **`/fapi/*` devolvio `-2015` desde la IP 191.104.7.185.** El diagnostico por
`/sapi/v1/account/apiRestrictions` lo separa sin ambiguedad:

```
enableReading    True      <- la clave es valida y llega
enableFutures    False     <- ESTE es el bloqueo
ipRestrict       False
enableWithdrawals / trading / margin ...  todos False
```

La clave es exactamente la de solo lectura que se pidio, pero **los endpoints de futuros
exigen `enableFutures` aparte**, y Binance solo deja activarlo sobre una clave con
**lista blanca de IP**.

**Lo que si se leyo, y cambia el estado del criterio:** `/api/v3/account` devolvio
`commissionRates` de spot = **0.00100000** maker y taker, que es exactamente el escalon
**VIP 0** de spot. El nivel VIP de Binance es **unificado** entre spot y futuros -- lo fija
el volumen a 30 dias y la tenencia de BNB, no el producto -- asi que **la cuenta esta en
VIP 0 tambien en futuros**, y de ahi salen `maker 0.0200 %` y `taker 0.0400 %`.

| pieza | procedencia |
|---|---|
| nivel VIP = 0 | ✅ **LEIDO de la cuenta de Mainnet** |
| tarifas de futuros para VIP 0 | tabla publica |
| descuento BNB en futuros | ⛔ **sin leer** (`/fapi/v1/feeBurn` necesita `enableFutures`) |

**El §8 pasa de INCUMPLIDO a PARCIAL.** Y lo que falta solo puede mover la comision hacia
ABAJO (−10 % si el BNB burn esta activo: `c(u)` maker de 4.8926 a 4.4926 pb, `R²_req` ×0.84).

⚠ **Eso importa para leer el §1, y en la direccion buena:** la cifra actual es una **cota
superior** del coste. Para una conclusion **negativa** como la del §1 es justo lo que se
quiere -- si no cruza con el coste maximo, tampoco cruzaria con el real. Solo morderia si la
conclusion fuera positiva.

⚠ **Y la taker se CORRIGE de 0.000500 a 0.000400**, por dos vias independientes que coinciden:
la lectura de Testnet y la derivacion del VIP 0 leido de Mainnet. El `0.0005` que arrastraban
`propagador.py`, `cola.py` y `tick_grande.py` es una tarifa **vieja**. Consecuencia:
`c(u)` taker+taker baja de 10.02 a **8.02 pb** y `R²_req` ×0.64. **Atencion al sentido: es un
cambio que AFLOJA un criterio**, y esos hay que mirarlos dos veces; aqui no toca ningun
veredicto porque todo el §1 se midio con maker+maker, que no cambia.

⚠ **Un control propio fallo al hacer el cambio, y estuvo bien que fallara.** El test tenia
clavado `taker+taker = 10.00 pb` y detecto que la constante se habia movido. Es exactamente
lo que impide que una tarifa cambie en silencio -- que es como el `0.0005` viejo sobrevivio
en cuatro modulos sin que nadie lo notara.

**Para cerrar el §8 del todo** hace falta, sobre la clave de Mainnet: activar
*«Restringir el acceso solo a IP de confianza»*, anadir la IP de la maquina de captura, y
entonces marcar **«Habilitar Futuros»** (Binance no lo ofrece sin lista blanca). Sin permisos
de trading ni de retiro. Con eso, `coste.py --leer` lee `commissionRate` y `feeBurn` reales.

---

## HOJA DE RUTA tras la sesión 2026-08-23 — qué falta, y el dimensionamiento de posición

Escrita a partir de lo medido el 2026-08-23, no de lo planeado antes. Todo lo que se
afirma aquí tiene su número en las secciones de esa sesión.

### 0. El estado en una frase

**Hay un modelo del RÉGIMEN y no hay una señal.** `ν` y `σ` están medidas, tienen ciclo con
período conocido, se pronostican a un paso (`R²` = 0.615) y traen barra de error viva. Lo
que no existe es `α`: el §1 midió `R² ≈ 0` con el predictor que su propio documento impone,
y su veredicto está **retractado** por estacionalidad, no confirmado. Sin `α` no hay nada
que dimensionar. **El dimensionamiento es aguas abajo de una señal que todavía no existe**,
y por eso lo de abajo se escribe como maquinaria condicional, no como plan de ejecución.

---

### 1. El dimensionamiento de posición

#### 1.1 La fórmula, y de dónde sale cada término

Con el coste **lineal** de la decisión de diseño 2 (2026-08-09) y penalización cuadrática
de riesgo, el objetivo por ciclo es

```
J(u) = -alpha*u + c*|u| + (1/2)*R*u^2
```

cuyo óptimo es una **banda muerta** seguida de una rampa:

```
|alpha| <= c        ->  u* = 0                       <- el bot se abstiene
|alpha| >  c        ->  u* = sign(alpha) * (|alpha| - c) / R
```

y con `R = gamma * sigma_H^2` (varianza al horizonte de tenencia `H`):

```
u* = sign(alpha) * (|alpha| - c) / (gamma * sigma_H^2)
u_final = sign(u*) * min( |u*|, techo_de_riesgo )        <- decisión de diseño 3: `min`, NO derivación
```

| término | qué es | estado hoy |
|---|---|---|
| `alpha` | ventaja esperada por unidad, en pb | ⛔ **NO EXISTE.** §1 retractado; el mejor `R²` con potencia fue ≈ 0 |
| `c` | coste ida y vuelta, en pb | ⚠ **4.00 pb con tarifas ASUMIDAS VIP 0**. Criterio incumplido desde la v3.1 |
| `sigma_H` | volatilidad al horizonte de tenencia | ✅ **pronosticable**, `R²` = 0.615 a un paso de 5 min |
| `gamma` | aversión al riesgo | decisión de política. **No se deriva del techo** (cancela `alpha`) |
| techo | `I_max`, `nocional_max_posicion` | ✅ vivo en la capa de riesgo de la v1.3 |

#### 1.2 Lo que la sesión de hoy cambia, y es lo más importante del apartado

**Un tamaño fijo está mal por dos órdenes de magnitud.** El nivel diario de `sigma` recorre
**13.75×**, así que `sigma^2` recorre **~190×**. Como `u* ∝ 1/sigma^2`, un tamaño constante
estaría mal dimensionado por ese factor entre el día más tranquilo y el más agitado. **El
dimensionamiento TIENE que ser condicional al pronóstico de `sigma`**, y hoy por primera
vez ese pronóstico existe.

#### 1.3 ⚠ Y el pronóstico puntual NO se puede enchufar tal cual: sobredimensiona 1.9×

El error de pronóstico de `log sigma` tiene `sd = 0.5676` y es aproximadamente normal en
logaritmos. Entonces, con `sigma_real = sigma_hat * exp(e)` y `e ~ N(0, 0.5676^2)`:

```
E[ 1/sigma_real^2 ]  =  (1/sigma_hat^2) * exp(2 * 0.5676^2)  =  1.905 / sigma_hat^2
```

**Enchufar `sigma_hat` en `1/sigma^2` sobredimensiona por un factor 1.9.** Es la misma
desigualdad de Jensen que ya mordió hoy en el NIS de la barra de error, y por tercera vez
en el proyecto (tras `q̄`/`G_0` y el propio NIS). No es un matiz: es casi el doble de
posición.

**Y la dispersión del pronóstico es mayor que su sesgo.** El error relativo mediano en
`sigma` es del **30 %**, que en `sigma^2` son **~69 %**; y el cuantil 90 del pronóstico está
en `2.07 * sigma_hat`, o sea **4.29× en `u`** entre dimensionar con la mediana y dimensionar
con el decil superior.

**Regla que sale de esto, y se escribe como decisión:** dimensionar con un **cuantil
superior** de la distribución predictiva de `sigma`, no con su mediana. Cuál cuantil es una
elección de política —el decil superior cuesta 4.3× de tamaño frente a la mediana— pero
usar la mediana **no** es una opción neutral: es una elección que sobredimensiona.

#### 1.4 La barra de error es viva, así que el cuantil también debe serlo

El error de pronóstico reciente predice el siguiente (corr **+0.594** a 2 h, `R²` fuera de
muestra **+0.292**). O sea que la anchura de la distribución predictiva **no es constante** y
se puede estimar en línea. Consecuencia operativa: en tramos donde el error reciente es
grande, el cuantil superior se aleja más y el tamaño baja **solo**, sin umbral fijo.

⚠ Con la reserva medida: la barra viva lleva `var(z) = 1.31` contra el 1.0 honesto, y
**empeora la cola lejana** (`|z| > 3` pasa del 1.27 % al 1.92 %). Es sobreconfiada justo en
los tramos tranquilos, que son donde un salto sorprende más. Cualquier dimensionamiento que
cuelgue de ella necesita **suelo** además de cuantil.

#### 1.5 El horizonte de tenencia entra al cuadrado, y su exponente es una banda

`sigma_H = sigma_1 * H^H_p` con `H_p` en **[0.372, 0.554]** (24 ajustes, 4 tramos). Anclando
en `H_ref = 100 s`, la banda contribuye a `sigma_H^2`:

| `H` | factor de incertidumbre en `sigma_H^2` por la banda de `H_p` |
|---|---|
| 300 s | **1.49×** |
| 1 h | 3.69× |
| 4 h | **6.10×** |

**Por debajo de ~10 min la banda de `H_p` no es la restricción vinculante** (manda el error
de pronóstico de `sigma`, 69 % en `sigma^2`); **por encima de una hora sí lo es**. Eso
importa porque el §1, si alguna vez cruza, cruzará en la banda de decenas de minutos.

#### 1.6 La granularidad convierte la rampa en escalera

`minQty = 0.001 BTC ≈ 96 USD` (decisión de diseño 4) es el paso mínimo expresable. Cerca de
la banda muerta, `u*` redondea a cero — lo cual es **correcto** (es abstenerse), pero
significa que **la banda muerta efectiva es más ancha que `c`** y hay que medirla, no
suponerla. Con `I_max = 0.5 BTC` quedan 500 pasos de resolución sobre el rango completo.

---

### 2. Lo que falta por hacer, en orden

#### 2.1 BLOQUEANTE — rehacer el §1 con lo aprendido

Sin esto no hay `alpha` y todo lo demás es maquinaria sin motor. Cuatro cambios, los cuatro
justificados por medición de hoy:

1. **Ventanas retrospectivas en SEGUNDOS**, no en `f*H*nu` ticks con `nu` escalar. `nu`
   recorre 0.55×–2.17× de su media, así que en la hora punta la ventana cubre el **46 %** de
   los segundos que dice cubrir.
2. **`sigma_1` y por tanto `R2_req` condicionales al PRONÓSTICO de `sigma`**, no a la hora de
   reloj. Es donde está el recorrido de 190×, contra el ~8× de la hora — y ese 8× además
   salió inflado por medias horarias confundidas con el nivel del día.
3. **Estratos declarados antes de mirar el `R²`.** Recomendación: terciles del `sigma`
   pronosticado, que es conocido *ex ante* y captura la palanca grande.
4. **Control positivo de potencia dentro de cada estrato.** Ya existe y está probado; hoy
   acreditó 11 de 15 filas y descartó las otras 4.

**Criterio de lectura, sin cambios:** el §1.4 se aplica literalmente y «no hay banda viable»
sigue siendo un desenlace admisible.

#### 2.2 `c(u)` REAL — criterio incumplido desde la v3.1

`/fapi/v1/commissionRate` es firmado. Con tarifas asumidas, todo el criterio económico
descansa en un número que nadie ha leído de la cuenta. Es barato de arreglar y lleva tres
versiones pendiente.

#### 2.3 Extender el pronóstico de `sigma` al horizonte de tenencia

Hoy `sigma` se pronostica a **un paso de 5 min**. El dimensionamiento necesita `sigma_H`
con `H` = el horizonte que el §1 señale. Hay que medir directamente a esa escala en vez de
extrapolar con `H_p`, precisamente porque `H_p` es una banda y su contribución crece con `H`
(§1.5).

#### 2.4 Cerrar la v4.1

- **§6** (micro-precio de Stoikov): sólo la fracción de ceros es barata. Hoy se midió que
  `y_mid` es nulo el **95–97 %** de las veces en los cuatro tramos, así que el margen que
  §6 persigue sigue ahí.
- **§5** (`C_respaldo`) — después de que el §1 diga a qué horizonte.
- **§4** (`η̂` con barrido de colapso) — deuda de reporte, no decisión.

#### 2.5 Datos

- La captura corre hasta el **2026-09-02**. Con **≥ 4 semanas** el perfil semanal deja de ser
  exploratorio por su propio criterio; hoy son 13 días y ~4 de fin de semana, y eso ya
  bastó para que los bloques de finde sobreajustaran en `sigma`.
- **`@depth` no se está capturando.** H1 (presión del libro) sólo se pudo contrastar con
  nivel 1 y salió que no aporta; para contrastarla de verdad haría falta profundidad a
  varios niveles. Es una decisión de captura, no de análisis.

#### 2.6 Lo que NO hay que hacer

- **No tocar `Micelio.py`** hasta que el §1 decida. Su auditoría de 2026-08-09 sigue vigente:
  el código está sano y su modelo no.
- **No refinar el predictor del §1** antes de rehacerlo bien. El §1.2 lo condiciona a que la
  curva «cruce o quede cerca», y con el predictor tonto el mejor margen fue **0.076**.
- **No derivar `gamma` del techo de riesgo** (decisión 3): cancela `alpha` y el tamaño deja
  de responder a la señal.
- **No construir ningún nulo sin separar antes nivel diario y forma intradía.** Tres nulos
  propios fallaron hoy por exactamente eso.

---

## AUDITORÍA DE `Micelio.py` (2026-08-09) — qué pasaría si se arrancara hoy

`Micelio.py` no se toca desde la v2.2. Desde entonces se han refutado varias de las cantidades
que lo gobiernan, y **el código no lo sabe**. Esta sección dice qué sigue en pie, qué está muerto
y qué haría el bot si alguien lo arrancara. **No es una lista de tareas** — no se toca nada hasta
que la v3.2 decida.

### ⚠ Lo que más importa: el objetivo de posición NO CONTIENE NINGUNA SEÑAL

La condición terminal de Loeper es `U = ½γ₀(S − S_ref)²` con `γ₀ = I_max/(S·ΔS_max)`. La
cobertura objetivo denominada en BTC sale de ahí por derivación directa:

```
S·∂U/∂S  =  S·γ₀·(S − S_ref)  =  I_max · (S − S_ref) / ΔS_max
```

**El inventario objetivo es una función lineal del desplazamiento respecto del nodo de fase,
escalada por el techo de riesgo. No hay `α` en ninguna parte.** El bot no compra porque espere
que el precio suba: compra porque el precio se ha alejado de `S_ref`, y compra exactamente
`I_max` cuando se aleja `ΔS_max`.

Eso es la decisión 3 de abajo llevada al límite: no es que `γ` derivado del techo de riesgo
*cancele* `α` — es que **`α` nunca entró en la formulación**. Y el ancla `S_ref` es un nodo de
fase de la EMD, que es justo lo que la v2.2, la v3.0 y la v3.2 refutaron.

### ⚠ Segundo: el bot no puede abstenerse

El coste del NMPC es **puramente cuadrático** (línea 1104):

```
J = Σ [ q_Δ·e_k² + q_inv(Ω)·I_k² + R_eff,k·(u_c,k² + u_v,k²) ]  +  p_Δ·e_N² + p_inv·I_N²
```

Sin término lineal en `|u|` **no hay banda muerta**: el óptimo de una cuadrática con objetivo no
nulo es siempre `u ≠ 0`. El bot opera *siempre* que haya desviación, por pequeña que sea, y solo
lo frenan las restricciones de caja y el freno de singularidad. Es exactamente lo que la decisión
2 corrige, y es la razón de que esa decisión sea v3.4 y no un detalle.

### Qué está muerto pero conectado, y qué lo salva

| cantidad | estado empírico | qué hace hoy en el código |
|---|---|---|
| `ω_m`, `ω_ang` | **sin sustento** (v2.2, v3.0, v3.2) | alimenta `A_arm` y `c²_vol = k·ω_m·ν` |
| `A_arm` / rama armónica | el oscilador **no existe** (`k = 0`, raíces reales) | se conmuta por `C` con histéresis |
| `C` (concentración espectral) | inflada por la escalera de ventana | decide la rama de `A` |
| nodos de fase → `S_ref` | son armónicos de la ventana (`T = 2L/k`) | **compuerta de todo el lazo** y ancla del objetivo |
| `Ω`, `Φ`, `Ψ` | **sin sustento** | `q_inv(Ω)` y `R_eff = R_base + κΩ²` |

**Lo que lo salva de hacer daño, y es un accidente afortunado:** la guarda de banda de la v2.1
declara `omega_valida` cierta solo el **12.5 %** de las ventanas. Cuando es falsa, `ω_ang` va a
NaN, `Ω` se congela y degrada a 0 pasados 120 s, y la rama de `A` cae a velocidad constante. O
sea que **el 87.5 % del tiempo el acoplamiento endógeno está efectivamente desconectado** y el
filtro corre como un EAKF de velocidad constante en reloj de ticks — que es justo lo que la v3.0
§6.2 lista como superviviente.

Dicho de otro modo: **el sistema funciona hoy porque su parte refutada casi nunca se activa.**

Y hay una compuerta más, en la línea 2486: `if is_burnt_in and dropout == 0 and nu > 0.0 and
S_ref > 0.0`. Sin un nodo de fase detectado **el bot no opera en absoluto**. La cadena entera
cuelga de un detector cuya base empírica cayó.

### Lo que sí sobrevive intacto

Coincide con la lista del §6.2 de la v3.0, y la auditoría lo confirma leyendo el código:

- **Reloj de transacciones** (Δn = 1), ingesta por lotes, deduplicación por `aggTradeId`,
  detección de huecos, `Q(Δt)` acumulada correctamente.
- **Capa de riesgo entera** (v1.3): 7 guardas con `causa_halt` distinguible, ruta de cierre que
  no reporta éxito sin posición plana confirmada, máquina de episodios con `DETENIDO` terminal.
- **EAKF** con corrección una vez por paquete, actualización multi-tasa, NIS y burn-in.
- **Loeper backward** y su condición CFL — el esquema es correcto; lo discutible es su
  condición terminal, no su integración.
- **Infraestructura**: instancia única por latido, seqlock, ring buffer SPSC con detección de
  sobrepaso, telemetría con `ỹ_k`, recuperación de memoria compartida huérfana en Windows.

### Defectos de código encontrados en esta pasada

Ninguno nuevo de corrección. Los que había siguen documentados en las secciones históricas.

⚠ **Lo que sí hay que anotar como riesgo latente:** `apply_filters` valida contra `minNotional`
y el orden de trabajo v4.0 (decisión 4) establece que **el filtro que ata es `minQty` = 0.001 BTC
≈ 96 USD**, no los 50 USDT del nocional. A precios actuales `minQty` es el doble de restrictivo.
No es un error de la v1.3 —entonces se midieron los dos— pero sí una cifra que envejeció.

### Conclusión de la auditoría

**El código está sano; su modelo no.** Lo que hay que cambiar cuando la v3.2 decida no son bugs:
son tres decisiones de diseño —el objetivo de posición sin señal, el coste sin término lineal, y
el ancla en un nodo de fase refutado— y las tres están ya identificadas y fuera del alcance de
esta tanda.

## Decisiones de diseño tomadas fuera de sesión (2026-08-09)

Acordadas en conversación entre el operador y Claude. **No son tareas**: son decisiones que
cierran discusiones abiertas y que hay que conocer antes de tocar lo que afectan.

**1. El problema de la secretaria queda RECHAZADO.** No se cumple ninguna de sus cinco premisas
—elección única, irrevocable, sin recuerdo, solo rango ordinal, objetivo «el mejor»— y el mercado
no es ninguna de esas cosas. La regla de **cuándo operar** es la banda muerta `|α| > c_efectivo`.
Donde sí hay un problema de parada óptima genuino es en **`τ*`: cuándo cancelar la orden maker y
cruzar**, que es el sucesor natural del §5 de la v4.0 en lazo cerrado.

**2. El coste del NMPC es LINEAL, no cuadrático.**

```
J = −α·u + c·(u⁺ + u⁻) + ½R(u⁺ − u⁻)²      con  u⁺, u⁻ ≥ 0
```

La complementariedad `u⁺·u⁻ = 0` sale **gratis** porque `c > 0`. **Sin el término lineal no hay
banda muerta y el bot nunca se abstiene** — opera siempre, aunque `α` sea ruido. Es **v3.4**:
toca `Micelio.py` y va **después** de que `c(u, estado)` exista.

⚠ Esto reabre la contradicción #1 del PDF por el otro lado. La Sec. 6.1 descartaba la norma L1
por no diferenciable en SQP; la formulación de arriba la recupera **sin** perder
diferenciabilidad, separando `u` en parte positiva y negativa. La Sec. 4.5 tenía razón en pedir
L1 y la 6.1 en rechazar la formulación ingenua.

**3. `γ` NO se deriva del techo de riesgo.** Hacerlo **cancela `α` algebraicamente** y el tamaño
deja de responder a la señal: el bot operaría el mismo tamaño con señal fuerte y con señal
nula. Son **dos términos separados con un `min`**, no uno derivado del otro.

**4. El filtro que manda es `minQty = 0.001 BTC ≈ 96 USD, no el `minNotional` de 50 USDT.** La
v1.3 §A midió los dos y se quedó con el análisis del nocional; a los precios actuales el que ata
es `minQty`. Se consulta `/fapi/v1/exchangeInfo` **al arrancar y sin cachear** — el propio
proyecto ya documentó que Testnet es 10× más fino que Mainnet y que calibrar contra el entorno
equivocado produce un sistema que funciona en pruebas y se degrada en producción.

**5. `PLAN_CAPITAL_5_0.md` es CONDICIONAL** y no se ejecuta hasta que pase el **paso 3** de la
regla de decisión de la v3.2 — el criterio económico, el único con dinero detrás.

### Deudas §7.2 y §7.5 saldadas (2026-08-27) — `actividad.py`, 5/5

#### §7.2 — La curva en U se resuelve a favor de ACTIVIDAD

`R²` fuera de muestra prediciendo `log σ(t+1)`, por bloques contiguos, con nulo por rotación:

| fragmento | n | `σ(t)` | **`ν(t)`** | ambos | nulo q95 |
|---|---|---|---|---|---|
| estacional_0 | 1 516 | +0.4978 | **+0.4711** | +0.5119 | 0.1663 |
| estacional_1 | 987 | +0.7029 | **+0.7540** | +0.7550 | 0.3286 |
| estacional_2 | 415 | +0.1594 | **+0.2485** | +0.2456 | 0.0401 |
| estacional_3 | 277 | +0.1395 | **+0.2322** | +0.2013 | 0.1297 |
| captura_v33 | 189 | +0.4358 | **+0.4427** | +0.4363 | 0.3853 |

**`ν` gana en 4 de 5 fragmentos.** La pregunta «curva en U: ¿volatilidad o actividad?» del §4.3
del traspaso **se retira de la lista, resuelta a favor de actividad**. Marco: Clark (1973) y
Tauchen & Pitts (1983) — **no es hallazgo**.

⚠ Reserva de potencia: por fragmento el `n` es de 189 a 1 516 casillas y los márgenes contra el
nulo van de **1.15×** (captura_v33) a 6.2× (estacional_2). El resultado fuerte es el **agrupado**
(`n = 3 081`, `R² = 0.6153`, **7.7×** el suelo); por fragmento la comparación `σ` contra `ν` es
pareada y sufre menos, pero los niveles absolutos de las dos últimas filas no se leen solos.

#### ⚠ §7.5 — NO hay pico en `ν`. Y el pico que apareció primero era MI NULO

Multitaper sobre el residuo de `log ν` tras retirar el ciclo diurno (que explica el 23.3 %),
tramo contiguo de 82.33 h, **nulo de ruido rojo SIMULADO** — nunca tabla asintótica.

**Primera corrida, con AR(1) como pedía el guion:** pico a **904 s**, exceso **×4.10**,
`p global = 0.0000`, y dentro de la banda de la v2.1 §2. Es decir: **hallazgo**.

**No lo era.** La firma no cuadraba: **63 de 494 frecuencias (12.8 %) sobre el q95 puntual**,
cuando un ciclo genuino da un pico *estrecho*, no un exceso repartido por el 13 % del espectro.
Barriendo el orden del nulo:

| nulo | frecuencias sobre q95 | pico | `p` global |
|---|---|---|---|
| AR(1) | 12.3 % | ×4.63 | 0.000 |
| AR(5) | 9.1 % | ×2.99 | 0.000 |
| **AR(20)** (elegido por AIC) | **3.8 %** | **×1.48** | **0.335** |

Y el exceso estaba concentrado en **600–900 s (26 % de esa banda sobre q95)**, que es donde un
AR(1) desajusta más: **banda ancha, no pico**. Un AR(1) no puede imitar la persistencia real de
`ν` (`ρ₁ = 0.8845`), y su desajuste se lee como señal.

**Con el orden elegido por AIC: NO HAY PICO.** La expectativa declarada antes de correr se
cumple, y la pregunta de ciclo endógeno abierta desde la v2.2 **queda cerrada sobre la serie de
actividad**, que es donde nunca se había contrastado.

⚠ **El patrón de la sesión, y conviene tenerlo escrito:** es la **cuarta vez en un día** que un
nulo demasiado estrecho estuvo a punto de producir un hallazgo falso — el barajado de `ε` en la
Adenda C (75–1538× demasiado estrecho), los tres nulos del 2026-08-23 que no separaban nivel
diario de forma intradía, y ahora el AR(1) del espectro. **La regla que sale: el nulo tiene que
reproducir la propiedad del dato que hace ancho al estadístico, y si no se sabe cuál es, se
barre el parámetro que la controla y se mira si el veredicto se mueve.**

---

### ⚠ ADENDA C CERRADA (2026-08-27) — no hay banda viable donde el instrumento resuelve

`identidad.py`, 7/7 controles. Suelo congelado (5 fragmentos × 9 horizontes) **antes** de
calcular ningún `R²` real, como exige el §C.4.1. Acta en `telemetria/acta_adendaC_medir.txt`.

#### El desenlace del §C.6, fila por fila

| `H` | desenlace | fragmentos sobre su propio suelo | falta un factor |
|---|---|---|---|
| 60 s | **NO HAY BANDA VIABLE** | **5 de 5** | **577×** |
| 120 s | NO HAY BANDA VIABLE | 2 de 5 | 359× |
| 300 s | NO HAY BANDA VIABLE | 1 de 5 | 203× |
| 600 s | NO HAY BANDA VIABLE | **0 de 5** | 127× |
| 900 s | NO HAY BANDA VIABLE | 0 de 5 | 127× |
| 1800 s | NO HAY BANDA VIABLE | 0 de 5 | 96× |
| 3600 s | NO HAY BANDA VIABLE | 0 de 5 | 43× |
| 7200 s | **el instrumento no resuelve** | — | *no se lee* |
| 14400 s | **el instrumento no resuelve** | — | *no se lee* |

**Tres lecturas, en orden de importancia:**

1. **A 60 s la señal es REAL y está medida**: los 5 fragmentos superan su propio suelo, con
   **signo positivo en los 5**. Coherente con el propagador de la v3.1 §2. **Y falta un factor
   577.** No es que no haya señal: es que no paga.
2. **El margen se cierra monótonamente con `H`** — 577 → 359 → 203 → 127 → 96 → **43** — y
   **exactamente donde se cerraría (≥ 2 h) el instrumento deja de resolver.** Es un dato de
   entrada directo para el §4 de la v4.2, cuyo óptimo previsto (P1) cae entre 2 y 8 h.
3. **De 10 min a 1 h el `R²` medido no supera su propio suelo en ningún fragmento** (0 de 5).
   En esa banda no hay medición, hay ruido — y el `R²` que se reporte ahí no significa nada.

#### ⚠ La corrección que hizo falta: la identidad hay que ponderarla por TIEMPO, no por ticks

**La compuerta de calibración del §C.4.3 lo cazó.** La identidad tal como el §C.3.2 la propone
—cada tick como origen— no reproducía el método de ventanas: daba **0.026** contra **0.0076** a
60 s, factor 3.4. Diagnóstico, con el mismo predictor y el mismo bloque de validación:

| orígenes | `R²` a 60 s | `R²` a 120 s | n |
|---|---|---|---|
| **cada tick** (lo que propone la adenda) | 0.028261 | 0.016233 | 2 223 663 |
| **uno por segundo de reloj** | **0.001634** | **0.000708** | 422 527 |
| ídem, sólo validación | 0.004515 | 0.003221 | 76 644 |
| método de ventanas (OOS, publicado) | 0.0076 | 0.0041 | 1 341 |

**Factor 17 entre ponderar por ticks y por tiempo.** Usar cada tick importa una **ponderación
por actividad**: los tramos con más transacciones aportan más parejas y son también los de más
volatilidad. **La curva requerida vive en tiempo de calendario** —la comisión se paga por ida y
vuelta, la volatilidad se acumula en segundos—, así que el estimando tiene que ser el ponderado
por tiempo. Es la misma lección que la sesión 2026-08-08 (e), donde la firma solapada
sobreponderaba los tramos activos y los dos estimadores discrepaban **en signo**.

Con orígenes uniformes, la identidad **sí** reproduce: 0.0045 contra 0.0076, y el propio §C.1 da
`q95 = 0.060` para las ventanas a 60 s, o sea que su 0.0076 está muy dentro de su propio ruido.
**La compuerta pasa.**

⚠ **Y la ventaja de resolución SOBREVIVE**: 422 527 orígenes uniformes contra 1 341 ventanas no
solapadas, **315×**. El §C.3.2 se sostiene; lo que había que corregir era el estimando, no el
método.

#### Otras dos correcciones propias

- ⚠ **`ε_t` no estaba en la multivariante**, así que el control «la multivariante no puede rendir
  menos que la univariante» **no era un teorema**: los agregados más cortos promedian ~10³ ticks
  y diluyen el signo suelto. Con `ε_t` como primera columna, la cota inferior del §C.3.4 se
  cumple por construcción.
- ⚠ **La regla del §C.6 exigía signo estable a magnitudes por debajo de su propio suelo.** El
  signo de una covarianza dominada por ruido es aleatorio: pedirle consistencia convierte una
  refutación limpia en «no decidible» y esconde el resultado. La estabilidad de signo **solo se
  exige a un resultado positivo**. Corregido el orden; sin eso, cinco de los nueve horizontes se
  habrían reportado como «no decidible» en vez de como refutación.

#### El suelo, con el estimador correcto

| fragmento | orígenes a 60 s | q95 a 60 s | q95 a 4 h |
|---|---|---|---|
| estacional_1 (20.0 M ticks) | 284 565 | **0.000049** | 0.000424 |
| estacional_0 (8.9 M) | 422 587 | 0.000129 | 0.000530 |
| estacional_2 (4.1 M) | 117 209 | 0.000222 | 0.002211 |
| estacional_3 (1.8 M) | 79 098 | 0.000582 | 0.003872 |
| captura_v33 (1.0 M) | 52 245 | 0.001352 | 0.007219 |

Razón rotación/barajado: **4× a 81×, mediana 25×**. Sigue justificando la desviación declarada
sobre el nulo del §C.4.1 — con el suelo de barajado, filas que no resuelven se leerían como si
resolvieran.

---

### v4.2 `ORDEN_TRABAJO_FRECUENCIA_4_2` — arranque (2026-08-27)

**El cambio de planteamiento:** de «¿existe señal a horizonte `H`?» a «¿existe un par
(tenencia `H`, umbral `θ`) con rentabilidad neta por unidad de tiempo positiva?». La comisión es
un **impuesto a la frecuencia**, y el óptimo es interior: existe aunque no sea rentable.

#### ⚠ §1.1 — COMPUERTA NO PASA. Nada de la v4.2 es citable todavía

```
/fapi/v1/commissionRate  ->  HTTP 401  -2015
permisos de la clave: enableReading=True   enableFutures=FALSE   ipRestrict=False
```

Lo que **sí** está leído de la cuenta de Mainnet (2026-08-23): **nivel VIP = 0**, vía
`commissionRates` de spot = 0.001. De ahí, por tabla pública, futuros VIP 0 = maker 0.0200 % /
taker 0.0400 %. Lo que **no** se puede leer: el **descuento BNB** (`/fapi/v1/feeBurn` también
exige `enableFutures`), que valdría −10 %.

⚠ **Nótese que el §4.4 de la orden construye su tabla de expectativas con `lastre = 4.49 pb`, que
es maker CON descuento BNB** (3.60 + 0.888). Si el BNB no está activo, el lastre real es
**4.89 pb** y toda esa tabla se desplaza en contra. Es exactamente el número que la compuerta
existe para fijar.

**Para abrir la compuerta**, sobre la clave de Mainnet: activar *«Restringir el acceso solo a IP
de confianza»*, añadir la IP de la máquina de captura, y entonces marcar **«Habilitar Futuros»**
— Binance no ofrece esa casilla sin lista blanca. Sin permisos de trading ni de retiro.

**Decisión de ejecución, declarada:** se ejecutan §1.2, §2, §3 y las deudas del §7, que **no
dependen** de `c(u)`. **NO se ejecuta el §4 (la superficie)** hasta que la compuerta pase: es la
etapa que produce un titular citable y su curva de coste entera cuelga del número asumido. `c(u)`
entra de forma **aditiva**, así que la superficie se recalcula en minutos cuando lleguen las
tarifas.

#### Estado de las deudas del §7 al arrancar

| deuda | estado |
|---|---|
| **7.1** calibración contra Cont–Kukanov–Stoikov | ✅ **HECHA** el 2026-08-27, ver apartado propio. El OFI-L1 **no** está degradado |
| **7.2** curva en U | implementada en `actividad.py`, pendiente de correr |
| 7.3 micro-precio | pendiente |
| 7.4 `η̂` con barrido de colapso | pendiente |
| **7.5** nulo espectral sobre `ν` | implementado en `actividad.py` con nulo **simulado**, pendiente de correr |

#### §1.2 — suelo de ruido de la Adenda C, medido sobre `estacional_0` (131.83 h)

| `H` | 60 s | 120 s | 300 s | 600 s | 900 s | 1800 s | 3600 s | 7200 s | 14400 s |
|---|---|---|---|---|---|---|---|---|---|
| **q95 rotación** | 0.000182 | 0.000192 | 0.000283 | 0.000337 | 0.000317 | 0.000313 | 0.000501 | 0.000553 | **0.000648** |
| q95 barajado | 0.000002 | 0.000003 | 0.000002 | 0.000002 | 0.000003 | 0.000002 | 0.000002 | 0.000002 | 0.000002 |
| razón | 82× | 75× | 131× | 151× | 125× | 167× | 227× | 275× | **292×** |

**Dos cosas quedan establecidas:**

1. **El suelo crece 3.6× mientras `H` crece 240×.** El método de ventanas iba como `k/n_val`, que
   en ese mismo rango empeora ~100×. **La tesis del §C.3.2 se sostiene sobre dato real**, y el
   orden de magnitud coincide con lo que la adenda había simulado (1.6e-4 a 3.8e-4).
2. ⚠ **Barajar `ε` subestima el suelo entre 75× y 292×.** El §C.4.1 pedía barajar; la desviación
   a **rotación circular** —declarada por conservar la memoria larga de `ε`, que es el modo de
   fallo 8 de la propia adenda— resultó **necesaria**: con el suelo de barajado, filas que **no**
   resuelven se habrían leído como si resolvieran.

⚠ **Decisión con consecuencia declarada: tope de 2 M de parejas**, con submuestreo **sistemático**
(no aleatorio: el aleatorio rompería el solapamiento, que es lo que da anchura al nulo) y **el
mismo conjunto para la estimación y para el suelo**. Se pierde resolución respecto a usar las
20 M y se reporta como tal; siguen siendo ~10⁵ veces más parejas que ventanas tenía el método
viejo.

#### Infraestructura: el techo de RAM, resuelto

Los cinco fragmentos pasan a `.npy` sueltos abiertos con `mmap_mode='r'` — **1.4 GB de disco**
contra 7.7 GB de RAM con ~1 GB libre. Tres procesos habían muerto por eso. **Lección
complementaria, que costó dos muertes más: un solo trabajo pesado a la vez** — el suelo y la
Tarea 1 corriendo en paralelo se mataron entre ellos.

---

### Tarea 1 (2026-08-27) — Calibración contra Cont, Kukanov & Stoikov (2014). `cont2014.py`, 6/6

⚠ **TODO ESTE APARTADO ES CONTEMPORÁNEO** (`Δmid_k` contra `OFI_k` del **mismo** intervalo de
10 s). No es comparable con la Adenda C ni con el §1 de la v4.1, que son predictivos. Ver la
convención de la Tarea 2.

#### ⚠ EL RESULTADO INVIERTE LA EXPECTATIVA DECLARADA: el instrumento NO está degradado

| serie | n_sub | OFI | OFI **sin** eventos que mueven el precio | transac. | ambos | cuadrát. | sig(TI) |
|---|---|---|---|---|---|---|---|
| **captura_v33** | 45 | **75.0 %** | **72.4 %** | 43.3 % | 84.3 % | 76.7 % | **91.1 %** |
| **estacional_tramo3** | 46 | **75.1 %** | **72.6 %** | 43.0 % | 84.3 % | 76.8 % | **91.3 %** |
| publicado (50 acciones NYSE) | 50 | 65.0 % | 35–60 % | 32.0 % | 67.0 % | 68.0 % | 31.0 % |

La interpretación se declaró antes de correr: *«un `R²` sustancialmente por debajo de 35–65 % no
es un hallazgo sobre BTCUSDT; es una medida de la degradación de nuestro OFI-L1 aproximado»*.
**Sale al revés.** Nuestro OFI-L1 desde `bookTicker` **reproduce y supera** las cifras publicadas.

**Consecuencias, y la primera es la que importa:**

1. ⚠ **El déficit predictivo NO se puede achacar al instrumento.** El residuo de reconciliación
   del 93.5 % (v3.2 §4.3) **no** se traduce en un OFI degradado. Explicar contemporáneamente
   funciona al **75 %**; predecir da **≈ 0**. Es la distinción de la Tarea 2 con número puesto,
   y **refuerza** la conclusión negativa del §1 en vez de debilitarla.
2. **Replica a tres cifras entre dos capturas separadas por semanas** (75.0/75.1, 72.4/72.6,
   84.3/84.3, 91.1/91.3). En un proyecto donde `H_p`, `β` y la antipersistencia **no**
   replicaron, esta relación es la primera que sí. Es una propiedad estable del mercado.
3. **La relación es LINEAL.** El término cuadrático `OFI·|OFI|` sube 75.0 → 76.8 %, en línea con
   su 65 → 68 %. **No introducir no linealidad**, como manda el artículo.
4. ⚠ **Dos diferencias reales con NYSE, no artefactos:**
   - **El desequilibrio de transacciones pesa mucho más aquí**: significativo en el **91 %** de
     submuestras contra su **31 %**, y por sí solo da 43 % contra su 32 %. En BTCUSDT el flujo
     de transacciones lleva información propia que en renta variable no lleva.
   - **Su control de tautología casi no muerde aquí**: 75.0 → 72.4 % contra su 65 → 35–60 %.
     Hipótesis, y hay que tratarla como tal: con el **spread clavado en 1 tick** (medido: mediana
     y p90 = 0.1000 USD en los cuatro tramos), el punto medio se mueve exactamente cuando una
     cola se vacía, y el vaciado es visible como caída de cantidad **a precio fijo** antes de que
     el precio salte. En NYSE, con spreads de varios ticks, ese canal es mucho más débil.

#### `λ` de `β_i = c / AD_i^λ` — replica, pero NO es la misma regresión que la suya

| serie | n | `λ̂` | `z` contra `λ = 1` | recorrido de `AD` (p10→p90) |
|---|---|---|---|---|
| captura_v33 | 45 | **+0.1180 ± 0.0519** | −16.98 | 1.92× |
| estacional_tramo3 | 46 | **+0.1172 ± 0.0503** | −17.57 | 1.95× |
| publicado | 50 acciones | ≈ 0.98 | no rechazable en 35/50 | órdenes de magnitud |

⚠ **NO se concluye «λ ≠ 1, luego la profundidad más allá del primer nivel domina», y la razón es
metodológica:** la suya es una regresión **TRANSVERSAL entre 50 activos distintos** con
profundidades que difieren en órdenes de magnitud; la nuestra es **TEMPORAL dentro de un solo
activo** sobre un recorrido de profundidad de **1.9×**. Es la misma forma funcional ajustada
sobre un eje de variación distinto, y no tienen por qué compartir exponente — sobre todo cuando
dentro del instrumento la profundidad y la actividad co-varían (`corr(log D, log ν) = −0.485`,
medido el 2026-08-23).

**Lo que sí queda establecido:** *dentro* de BTCUSDT e intradía, **el coeficiente de impacto es
casi independiente de la profundidad de nivel 1** (`λ = 0.12`, replicado). Y eso coincide con
H1, que por otra vía midió que la profundidad de L1 **no aporta nada** sobre el flujo (parcial
+0.0425 contra un techo de nulo de +0.0418). **Dos medidas independientes dicen que la liquidez
que gobierna no está en L1.**

**Eso sí es un argumento para ingerir `@depth`** — el primero que existe en el proyecto — pero
formulado así, no como «λ ≠ 1». Y sigue siendo una decisión de captura, no de análisis.

---

## ⚠ CONVENCIÓN OBLIGATORIA — todo `R²` se cita como CONTEMPORÁNEO o PREDICTIVO

**Ninguna cifra de `R²` de este proyecto se escribe sin una de esas dos etiquetas.** Son
cantidades distintas y confundirlas ha estado a punto de costar una lectura equivocada más de
una vez.

| etiqueta | qué regresa | ejemplo del proyecto | ejemplo de la literatura |
|---|---|---|---|
| **CONTEMPORÁNEO** | `Δp_k` contra el flujo del **mismo** intervalo | el `R² = 0.014` de la v3.2 (núcleo con `h(0)·x_t` dentro) | **65 %** de Cont, Kukanov & Stoikov (2014) |
| **PREDICTIVO** | `r_{t→t+H}` futuro contra flujo **pasado** | el `R² ≈ 0` del §1 de la v4.1; toda la Adenda C | — |

**El 65 % de la literatura es CONTEMPORÁNEO y no es comparable con nada de la Adenda C.**
Explicar el movimiento que ya ocurrió no es predecir el que viene, y sólo lo segundo se puede
negociar. Es la misma reconciliación que la sesión 2026-08-10 (c) escribió para el `R² = 0.014`
de la v3.2 contra el `R² ≈ 0` del §1, elevada aquí a regla.

Corolario práctico: **un `R²` contemporáneo alto no rebaja el requisito del paso 3.** El peaje
se paga por entrar y salir, y para eso hace falta saber antes.

---

## Decisión pendiente para Samuel — de dónde saldría el predictor de la banda abierta

**No se implementa nada de esto sin orden de trabajo y preregistro nuevos.** Queda escrito para
que la decisión exista cuando toque tomarla.

**El problema de escalas.** El desequilibrio del libro es un predictor de **segundos**: el propio
Cont-Kukanov-Stoikov reporta que sus autocorrelaciones se desvanecen hacia los 10 s. La banda de
horizonte que la v4.1 dejó abierta es de **15 a 60 minutos**. Son escalas incompatibles por dos
o tres órdenes de magnitud. La Adenda C lo va a confirmar o refutar con la curva `R²(H)` por
identidad de covarianza — es exactamente lo que esa curva mide.

**Si se confirma que la señal de libro no sobrevive más allá de unos minutos**, el predictor
para la banda abierta tiene que venir de variables **exógenas** —tasa de financiación del
perpetuo, base contra CME, DXY— y **no de microestructura**. Sería un cambio de familia de
datos, no un refinamiento del estimador.

### ⚠ El argumento estructural, que no requiere medición y conviene tener escrito

**El nivel de equilibrio de la predictibilidad del flujo lo fija el coste del participante
marginal, no el nuestro.**

Un participante de VIP alto paga una comisión maker cercana a cero; nosotros pagamos **4 pb de
ida y vuelta a VIP 0** — y ese VIP 0 está **leído de la cuenta** (2026-08-23). Mientras exista
alguien capaz de operar rentablemente con un margen mucho menor que el nuestro, la
predictibilidad se compite **hasta su suelo, no hasta el nuestro**: entre uno y dos órdenes de
magnitud por debajo de nuestro umbral.

**Consecuencia dura: mejorar el estimador no mueve ese suelo.** No es un problema de método.

Y es coherente con las tres cosas que este proyecto ya midió por separado:

- el **paso 3 de la v3.2**, que se paró por un factor 3.46 en comisiones;
- **`max|μ̂| = 2.63 pb` contra 4 pb de comisión** — ni el máximo de la señal medida llega al
  peaje;
- el **§1 de la v4.1 a 60 y 120 s**, donde la medición sí resuelve y da `R² ≈ 0` contra un
  requisito del 38–83 %.

Tres mediciones independientes, con instrumentos distintos, apuntando al mismo sitio.

---

## Convenciones

- **Todo `R²` se etiqueta CONTEMPORÁNEO o PREDICTIVO.** Son cantidades distintas y no
  comparables; ver la sección propia más arriba. El 65 % de Cont-Kukanov-Stoikov es
  contemporáneo.
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
