# Orden de Trabajo — Fase de Calibración y Diagnóstico

**Contexto:** las correcciones estructurales de la sesión 2026-08-02 están aplicadas y verificadas
(ver `CLAUDE.md`). Este documento cubre la **fase siguiente**: instrumentar el sistema para saber
si el filtro está bien calibrado, y calibrarlo con métodos que no sean prueba y error.

Basado en una revisión de literatura externa (documento de calibración HFT) **filtrada por
aplicabilidad**. Varias de las técnicas que ese documento recomienda NO aplican todavía y están
listadas explícitamente en la sección "No hacer aún" — leerla antes de empezar.

---

## 0. Constantes de acoplamiento — RESUELTAS, implementar como fórmulas

El análisis dimensional de las constantes de interconexión del modelo está **cerrado y aprobado por
el autor**. Ya no son valores libres: cada una se **deriva** de un límite físico o de capital del
sistema.

**Regla dura: implementarlas como fórmulas evaluadas a partir de los límites estructurales, nunca
como literales numéricos.** Un número mágico en el código vuelve a abrir el agujero que esto cierra.

Consolidar todo en un módulo nuevo, `constantes_micelio.py`, con una entrada por constante que
incluya símbolo, sección del PDF, unidades y la fórmula. El resto del código importa desde ahí;
nadie más las define.

### 0.1 Constantes derivadas (aprobadas)

| Constante | Fórmula | Unidades | Sec. PDF |
|---|---|---|---|
| `γ_0` | `I_max / (S · ΔS_max)` | BTC³/USD² | 7.4 (payoff terminal) |
| `γ_ω` | `1 / ω_m,max` | Ticks | 7.3.3 |
| `γ_Q` | `1 / max(max(ΣQ_T̄), C_max)` | 1/USD | 7.3.3 |
| `κ` | `c · R_base / Ω_crit²` | Ticks⁴/BTC² | 6.1 |
| `μ` | `c' · q_base / Ω_crit²` | Ticks⁴/BTC² | 6.1 |

### 0.2 Límites estructurales de los que dependen

- `I_max` — exposición máxima absoluta de inventario [BTC]
- `C_max` — capital operativo total [USD]
- `ω_m,max` — mayor frecuencia modal vía HHT [1/Ticks]
- `ΔS_max` — desplazamiento máximo esperado respecto al nodo de fase [USD/BTC]
- `Ω_crit` — umbral de crisis estructural [BTC/Ticks²] (ver Fase 3)
- `c`, `c'` — factores de escala adimensionales (sugerido `c = c' = 100`: en `Ω_crit` el costo de
  ejecución vale 101× el fee base, así el solver difiere órdenes *antes* de que las restricciones
  de caja se activen — la "aversión preventiva" de la Sec. 6.2)

### 0.3 Notas de derivación que deben quedar en comentarios

**`γ_0` es un Gamma, no un Delta.** Como `Γ = ∂²U/∂S² = γ_0` exactamente, sus unidades son las de
`Γ` en la Sec. 4.4.1: BTC³/USD². La fórmula sale de exigir que la cobertura en BTC (`S·∂U/∂S`)
sature `I_max` cuando el precio se aleja `ΔS_max` del nodo de fase. Una versión anterior con
unidades de Delta (BTC²/USD) daba `λS²Γ ≈ 25` → **singularidad de Loeper permanente desde el primer
ciclo**, y una cobertura implícita de 169 BTC contra un límite de 0.5 BTC.

**`γ_ω` tiene un piso teórico.** `ω_m` es la inversa de los ticks para completar un ciclo, luego el
mínimo son 2 ticks (Nyquist en espacio de ticks) → `ω_m,max ≤ 0.5 /Ticks` → `γ_ω ≥ 2 Ticks`. Usarlo
como sanity check si el estimador HHT devuelve algo fuera de rango.

**`γ_Q` usa `max`, no `min`.** `ΣQ` es volumen de **mercado** (Secs. 1.1 y 6.3), órdenes de magnitud
mayor que `C_max`. Tomar el `min` del denominador *maximiza* `γ_Q`, que es lo contrario de acotar.
Con `max`, el término queda ≤1 mientras `ΣQ` no supere su máximo histórico, y `C_max` actúa como
piso del denominador durante el arranque en frío cuando aún no hay historial.

**`κ` y `μ` se anclan a `Ω_crit`, no a un `Ω_max`.** `Ω` no tiene cota superior natural: por la
Sec. 2.6, en el nodo de fase `ΔS → 0` → `Ψ → ∞` → `Ω → ∞`. Y eso es intencional (Sec. 6.1: `R(Ω)`
debe divergir para forzar la inacción). Por tanto el papel de `κ` no es normalizar a ≤1 sino fijar
**dónde muerde** la divergencia — y ese punto ya lo define `Ω_crit`.

### 0.4 Conversión de unidades de `ω_m` — pendiente de implementación

La cadena EMD → Hilbert calcula frecuencia en **Hz** (tiempo físico), pero el modelo exige
**1/Ticks**. La conversión es `ω_m = f / ν`, con `ν` la tasa de llegada de ticks de la Sec. 4.5.

Sin esto, `γ_ω·ω_m` no es adimensional en ejecución aunque lo sea en el papel. Implementar junto con
la integración de la cadena EMD.

### 0.5 Guardas de ejecución (implementar)

Dos asserts baratos que delatan un error dimensional en el primer ciclo en vez de a los veinte
minutos:

1. **Cobertura acotada:** `S · γ_0 · ΔS_max ≤ I_max` (con tolerancia). Si falla, `γ_0` está mal
   escalado.
2. **Fricción sub-crítica al arranque:** `λ · S² · γ_0 < 1`. Si falla en el primer ciclo, hay
   singularidad de Loeper por construcción, no por condiciones de mercado.

Ambas con mensaje explícito nombrando la constante sospechosa.

### 0.6 NO tocar — constantes de ruido, pendientes de datos de cuenta demo

Estas quedan **deliberadamente sin revisar** hasta tener datos reales de la cuenta demo, que es de
donde saldrá su calibración:

`σ_OU`, `θ_OU`, `μ_OU`, `λ_min`, `η`, `β` (sensibilidad al jitter)

Claude Code no debe cambiar sus valores. Si una tarea parece exigirlo, dejar un
`# TODO(calibración demo):` describiendo la restricción encontrada y seguir adelante.

Sí mantener la guarda ya conocida sobre `σ_OU`: la desviación estacionaria de un OU es `σ/√(2θ)`, y
para que `λ` fluctúe *alrededor* de `μ_OU` hace falta `σ/√(2θ) ≪ μ_OU`. Verificarlo al arranque y
advertir, sin corregir el valor.

### 0.7 Adimensionales por declaración

`q_Δ`, `q_base` y `R_base` son **adimensionales**. Esto no es cosmético: es lo que hace que los tres
términos de `J` sean conmensurables (todos quedan en BTC²·adimensional una vez `e_k` está en BTC).
Declararlo explícitamente en el módulo.

---

## Principio de secuencia (importante)

Diagnosticar → calibrar → meta-optimizar. **En ese orden, sin saltos.**

La razón no es estilística. ALS (Fase 2) es matemáticamente válido **sólo si `A` y `H` son
correctas y las covarianzas de ruido son las únicas incógnitas** (Odelson et al. 2006, y toda la
literatura posterior asume lo mismo). Si `A` tiene estructura no modelada, ALS la absorberá dentro
de `Q` y devolverá una respuesta confiadamente equivocada. El test de Ljung-Box de la Fase 1 es
precisamente lo que detecta ese caso. Por eso va primero.

Y la meta-optimización (Fase 3+) sobre mocks sintéticos ajustaría hiperparámetros al generador de
mocks, no al mercado.

---

## FASE 1 — Diagnóstico de consistencia del filtro

**Prioridad máxima. Barato, sin dependencias nuevas más allá de `scipy.stats`, y desbloquea todo
lo demás.**

### 1.1 NIS online — y la unificación con la Sec. 2.2.1

La Secuencia de Innovación Normalizada al Cuadrado es:

```
ε_k = ỹ_kᵀ S_k⁻¹ ỹ_k
```

**Esto ya está en el PDF.** La Sec. 2.2.1 lo define con ese mismo símbolo `ε_k` y esa misma fórmula
para modular la ventana adaptativa del EMD (`α = e^{−γε_k}`, `W_k = ⌊W_min + α(W_max−W_min)⌉`).
Es el mismo escalar con dos usos: diagnóstico de consistencia del filtro **y** driver de la ventana.

Implementar:

- Calcular `ε_k` **en línea en el Hilo Rápido**, donde `S_k` ya existe. Usar `np.linalg.solve`,
  no invertir.
- Publicarlo en memoria compartida (campo nuevo en `MICELIO_DTYPE` o en el bloque de entorno).
- **Consumirlo en dos sitios:**
  - Telemetría: registrar el escalar `ε_k` (no la matriz `S_k` — logear matrices es caro e
    innecesario si el escalar se computa en línea).
  - Hilo Lento: implementar la ventana adaptativa `W_k` de la Sec. 2.2.1, que hasta ahora nunca
    se implementó (el EMD usa ventana fija).

**Interpretación:** con `m = 2` (rango del vector de medición), si el filtro es consistente
`ε_k ~ χ²` con 2 grados de libertad → media teórica = 2. Promedio persistentemente por debajo de
2 → filtro conservador (sobreestima su incertidumbre). Por encima → sobreconfiado, divergencia
inminente.

### 1.2 ⭐ Reemplazar `ε_burn` por el criterio NIS

**Este es el resultado más valioso de toda la revisión de literatura. Léelo con atención.**

`CLAUDE.md` deja abierto que `ε_burn` no puede ser constante: depende de `Δt`, de `ρ_k` (que el
Micelio modula tick a tick) y del régimen de volatilidad. Un umbral fijo bloquea el burn-in en
mercado agitado y lo deja pasar prematuramente en mercado tranquilo.

**El NIS resuelve esto de raíz porque es adimensional y auto-normalizado.** Su distribución de
referencia (χ² con `m` g.l.) no depende de `Δt`, ni de `ρ_k`, ni del régimen. Es la misma en todo
momento.

Criterio propuesto:

```
burn-in completo  ⟺  media móvil de ε_k dentro de la banda de confianza χ²_m
                     durante N ticks continuos
```

Implementación sugerida: ventana móvil de ~100 muestras, banda de dos colas al 95% de
`χ²_m / m` (usar `scipy.stats.chi2.ppf`). Mantener la derivada de la traza como criterio
**secundario** y registrar ambos en telemetría durante un tiempo, para poder contrastarlos.

Documentar en el código que esto sustituye el criterio de la Sec. 7.1 del PDF, con la razón.
Es un cambio al diseño, no solo a la implementación.

### 1.3 Ljung-Box multivariante sobre la innovación

Test de blancura sobre la secuencia de innovación normalizada. Estadístico:

```
Q_LB = n(n+2) · Σ_{j=1..h} ρ̂_j² / (n−j)
```

Comparar contra `χ²` con `h` grados de libertad. Ejecutar **offline**, sobre los volcados de
telemetría, no en el lazo de control.

**Qué significa el rechazo, y por qué importa más que el resto:** si `Q_LB` excede el crítico, hay
predictibilidad remanente que la matriz cinemática `A` no capturó. Eso es *model mismatch*, no
mala calibración de ruido. **Recalibrar `Q` y `R` en ese caso es tapar el síntoma:** ALS absorbería
el error de modelo dentro de `Q`. La respuesta correcta es aumento de estado (Sec. 3.6 del PDF —
derivadas de orden superior, sesgos de fricción).

**Ljung-Box es la compuerta de entrada a la Fase 2.** Si no pasa, no correr ALS.

### 1.4 Shapiro-Wilk sobre los residuales

Test de normalidad (`scipy.stats.shapiro`). Errores asimétricos o de colas pesadas justifican pasar
a ALS-IRLS con umbrales de Huber en vez de ALS estándar. Es el criterio para decidir si la
Fase 2.3 vale la pena; no un fin en sí mismo.

### 1.5 Entregable

Un módulo `diagnostico.py` + un script que consuma los volcados de telemetría y emita un reporte:
NIS medio con su banda, resultado de Ljung-Box a varios rezagos, Shapiro-Wilk, y series temporales
de los tres. Que se pueda correr después de cada sesión de Testnet.

---

## FASE 2 — Calibración offline de `Q_base` y `R_base` por ALS

**Precondición dura: Ljung-Box debe pasar. Si no, volver a 1.3.**

### 2.1 Contexto: esto contradice el PDF, y el PDF está desactualizado

Las Secs. 3.5 y 8.6.1 mandan **Covariance Matching** (Mehra, 1970) como método de calibración
offline — de hecho es el propósito declarado de logear `ỹ_k`. La literatura moderna es clara en que
Covariance Matching tiene tres defectos para tiempo real: no garantiza convergencia global, produce
estimaciones de varianza alta, y **frecuentemente devuelve matrices no semidefinidas positivas**.
Lo último es fatal: valores propios negativos en `Q` o `R` revientan la descomposición de Cholesky.

ALS (Odelson, Rajamani & Rawlings, *Automatica* 42(2):303–308, 2006) colapsa el procedimiento de
tres pasos de Mehra en una sola optimización lineal, con varianza menor y condiciones de unicidad
demostradas.

**Acción documental:** marcar las Secs. 3.5 y 8.6.1 del PDF como superadas. La telemetría de `ỹ_k`
sigue siendo correcta y necesaria — cambia el método que la consume, no el dato.

### 2.2 Implementación de ALS

Formulación: definir la autocovarianza de la innovación con rezago `j`, `C_j = E[ỹ_{k+j} ỹ_kᵀ]`, y
mediante álgebra de `vec` y producto de Kronecker armar el sistema lineal `𝒜χ = b`, donde `χ`
agrupa las incógnitas de `Q` y `R`.

Tres advertencias, todas verificadas en la literatura:

1. **Verificar la condición de unicidad antes de confiar en el resultado.** Hay resultados
   publicados de que las covarianzas de ruido de sistemas con entradas desconocidas **no son
   identificables de forma única** vía ALS. Tu `R_n` aumentado es efectivamente un canal de
   entrada desconocida / sesgo. Comprobar el rango de `𝒜` y reportarlo explícitamente. Si el
   sistema es rank-deficient, decirlo en vez de devolver una solución arbitraria.

2. **ALS calibra `Q_base`, no `ρ_k`.** La modulación endógena del Micelio
   (`ρ_k = 1 + γ_ω|ω_m| + γ_Q|ΣQ|`) se mantiene exactamente como está. Que no compitan por la misma
   escala: ALS fija la base térmica, el Micelio la expande dinámicamente. Documentarlo, porque es
   fácil de confundir.

3. **Positividad.** La versión completa inyecta ALS en un programa semidefinido (SDP) con
   restricciones `Q ⪰ 0`, `R ⪰ 0`. Eso implica dependencia nueva (`cvxpy` + solver de punto
   interior). **Empezar sin SDP:** ALS sin restricción + proyección de valores propios negativos a
   cero, y verificar cuánto se corrige. Añadir SDP sólo si la proyección resulta ser agresiva.
   Menos dependencias hasta que se demuestre que hacen falta.

### 2.3 ALS-IRLS — sólo si Shapiro-Wilk lo justifica

Si los residuales muestran colas pesadas (que es lo esperable en cripto), existe una variante
robusta: **ALS-IRLS** (Li & Deng, arXiv:2603.08158, marzo 2026; código en
`github.com/jiahongljh/als-irls`). Sustituye el criterio cuadrático por la función de costo de
Huber y añade umbralización adaptativa sobre las innovaciones crudas.

**Corre ALS estándar primero.** Sin baseline no se puede medir la mejora. El paper reporta reducción
de RMSE de más de dos órdenes de magnitud, pero eso es en simulación controlada con un modelo
ε-contamination específico — no diseñes esperando ese número.

---

## FASE 3 — `Ω_crit` por búsqueda de orden cero

`Ω_crit` (el umbral de crisis estructural que dispara las restricciones de caja asimétricas de la
Sec. 6.2) **no se puede ajustar por gradiente**: rompe la diferenciabilidad de la trayectoria KKT.
Es genuinamente el único hiperparámetro que necesita un método sin derivadas.

Pero eso **no** justifica traer CMA-ES. Es un escalar. Barrido 1-D o bisección sobre datos de
Testnet, con criterio explícito (p.ej. maximizar detección de transiciones de fase reales
minimizando falsos bloqueos). Documentar el criterio antes de correr la búsqueda, no después.

**`Ω_crit` ahora tiene efecto en cascada.** Tras la Sección 0, `κ` y `μ` se derivan de él
(`κ = c·R_base/Ω_crit²`, `μ = c'·q_base/Ω_crit²`). Consecuencias para la implementación:

- Los tres deben recalcularse juntos. `Ω_crit` no puede quedar como valor suelto en un sitio y
  `κ` hardcodeado en otro.
- La búsqueda de orden cero sobre `Ω_crit` calibra los tres de una vez — no hace falta una búsqueda
  aparte para `κ` y `μ`.
- El criterio de búsqueda debe evaluar el sistema **completo** (costo de ejecución + aversión de
  inventario + restricciones de caja), no sólo el punto de disparo de las cajas, porque mover
  `Ω_crit` mueve simultáneamente la pendiente de `R(Ω)` y de `q_inv(Ω)`.

---

## NO HACER AÚN — y por qué

El documento de investigación recomienda estas técnicas. Todas son reales y bien fundamentadas.
Ninguna aplica todavía.

### CMA-ES / meta-optimización global
**Bloqueado por falta de datos reales.** Los mocks son sintéticos — el generador de precio es una
señal con ciclo estructural de 150 USD / 40 s que *nosotros* elegimos. Correr CMA-ES contra eso
ajusta hiperparámetros al generador, no al mercado. Se necesita historial tick de Mainnet primero.

Nota: el documento atribuye CMA-ES al marco *AutoQuant* (arXiv:2512.22476). **AutoQuant usa
optimización bayesiana, no evolutiva.** Y excluye explícitamente impacto de mercado y restricciones
de capacidad — justo lo que el modelo de Loeper existe para capturar. Es un buen template de
*screening con costos realistas*, no de validación de este sistema.

### IMPC / Control Óptimo Inverso
**Estructuralmente inaplicable.** El IOC recupera la función de costo subyacente a partir de
**trayectorias expertas observadas**. No hay dataset de demostraciones expertas de cobertura. No hay
de qué invertir. Si algún día existe un historial de operación manual bien ejecutada, se reconsidera.

### DeepONet / transporte óptimo de semimartingalas para `c²_vol`
**Contradice tu propio modelo.** Ya tienes `c²_vol = k·ω_m·ν` derivado de la Canción del Micelio, con
consistencia dimensional verificada (`[1/Ticks]·[Ticks/Años] = [1/Años]`). Sustituirlo por un
operador neuronal aprendido descarta la tesis central del sistema. No hacer.

### PBO/CSCV, Ratio de Sharpe Deflactado, Walk-Forward, Stationary Bootstrap
Correctos e importantes — pero son la **compuerta de salida de la meta-optimización**. No hay nada
que deflactar hasta que existan backtests reales sobre datos reales. Cuando llegue ese momento, son
obligatorios, no opcionales.

---

## Pendientes heredados que esta fase no cubre

- `ω_m` y `R_n` siguen siendo mocks. La cadena EMD → Hilbert está validada (7–12 ms) pero vive fuera
  de `Micelio.py`. **Integrarla es precondición para que el NIS y la ventana adaptativa signifiquen
  algo**, porque `R_n` entra directamente en el vector de medición.
- Sin validar contra CUDA ni acados. Los núcleos son referencia CPU con el mismo layout.
- Las contradicciones documentales #1 (L1 vs cuadrática) y #2 (dimensiones del NMPC) siguen abiertas
  en el PDF.

---

## Convenciones

Las de `CLAUDE.md`. Recordatorio de las dos que más aplican aquí:

- Comentarios y variables en español; referenciar la sección del PDF al implementar una fórmula.
- Toda suposición que rellene un hueco del PDF va marcada con `# NOTA DE INTERPRETACION:`.

Y una específica de esta fase: **cuando una decisión contradiga el PDF** (como ALS sobre Covariance
Matching, o NIS sobre `ε_burn`), no la apliques en silencio. Marca
`# DIVERGE DEL PDF (Sec X.Y):` con la justificación. El PDF es el documento de diseño y hay que
poder reconciliarlo después.
