# Orden de Trabajo — Calibración v1.2

Sustituye a `ORDEN_TRABAJO_CALIBRACION_1.1.md`. Ese documento se conserva como registro.

**Estado:** Sección 0 y Fase 1 de la v1.1 están **ejecutadas y verificadas** (ver `CLAUDE.md`,
sesión 2026-08-02 (b)). Este documento cubre lo que sigue, y corrige dos errores del anterior.

---

## A. Correcciones a la v1.1 (errores del documento, no del código)

### A.1 Colisiones de nombres que introdujo el orden de trabajo

Ambas fueron detectadas en implementación y están bien documentadas en `constantes_micelio.py`.
Se resuelven **renombrando en el orden de trabajo**, porque el PDF es la fuente y llegó primero:

| Símbolo v1.1 | Conflicto | Nombre nuevo |
|---|---|---|
| `C_max` (USD, en `γ_Q`) | El PDF Sec. 6.2 usa `C_max` en **BTC** (capacidad de billetera para las cajas) | **`K_USD`** — capital operativo total en USD |
| `ΔS_max` (USD/BTC, en `γ_0`) | El PDF Sec. 7.4.1 usa `ΔS_max` como **margen porcentual** de la malla | **`ΔS_ref`** — desplazamiento absoluto esperado respecto al nodo de fase [USD/BTC] |

Fórmulas corregidas:

```
γ_0 = I_max / (S · ΔS_ref)              [BTC³/USD²]
γ_Q = 1 / max(max(ΣQ_T̄), K_USD)         [1/USD]
```

`C_max` y `ΔS_max` quedan reservados **exclusivamente** para su significado del PDF. Mantener
`verificar_dominio_malla` como guarda.

### A.2 Ljung-Box: reportar magnitud, no solo p-valor

El test es muy potente con n grande y rechaza por correlaciones sin relevancia práctica:

| n | ρ₁ mínimo que rechaza (h=20, α=0.05) |
|---|---|
| 500 | 0.250 |
| 2 000 | 0.125 |
| 11 500 (≈115 s a 100 Hz) | **0.052** |
| 50 000 | 0.025 |

**`diagnostico.py` debe imprimir `ρ₁` (y `ρ₂`, `ρ₃`) junto al p-valor**, y el veredicto debe
leerse sobre la magnitud. Con ρ₁ ≈ 0.8 no hay ambigüedad; con ρ₁ ≈ 0.03 el p-valor diría lo
mismo y significaría algo muy distinto.

Añadir a la salida el `n` efectivo tras el recorte de transitorio.

---

## B. Prioridad 1 — Integrar la cadena EMD → Hilbert

**Es la precondición de todo lo demás y bloquea las Fases 2 y 3.**

Mientras `R_n` sea un coseno puro y el precio una sinusoide de 40 s, el NIS y el Ljung-Box están
midiendo el **generador de mocks**, no el filtro. El `NIS = 13.93` reportado es real como medida,
pero su causa es artificial.

Tareas:

1. Traer la implementación validada (extensión por reflexión + splines cúbicos, 7–12 ms para
   ventanas de 64–256) al Hilo Lento.
2. Aplicar la conversión de la Sec. 0.4: el HHT entrega `f` en **Hz**, el modelo exige
   **1/Ticks**. `ω_m = f/ν`, con `ν` en Ticks/Años (usar `CTE.omega_m_desde_hz`).
3. Colapso espectral de la Sec. 2.5 (ponderado por energía, excluyendo IMF1 y residuo).
4. `R_n` pasa a ser el residuo real del EMD, no un coseno.
5. Detección de nodos de fase (Sec. 2.6) sobre la fase desenrollada real.
6. La ventana adaptativa `W_k` ya implementada pasa a alimentarse del NIS real.

**Expectativa:** buena parte del exceso de NIS debería evaporarse sola, porque `R_n` real es
lentamente variable — que es justo lo que `A` asume para ese estado.

---

## C. Prioridad 2 — `q_S` como varianza relativa

Sustituir la varianza absoluta de precio del ruido de proceso por una parametrización relativa:

```
q_S = (σ_rel · S_k)²          [σ_rel] adimensional
```

Unidades idénticas a las que `q_S` ya debe tener — `(USD/BTC)²` — así que el álgebra del filtro
no cambia. Ventaja: `q_S` deja de depender del nivel de precio. Con `q_S` fijo, un movimiento de
BTC de 45k a 90k lo deja mal por un factor de 4; con `σ_rel` fijo se corrige solo.

`σ_rel` entra en `constantes_micelio.py` como constante adimensional anclada a un observable
(volatilidad realizada del par). Marcar `[CALIBRAR]`.

**No es una adimensionalización del filtro.** Se evaluó y descartó normalizar el EAKF completo:
escalar `Q` y `R` sin escalar `x`, `P`, `z` y `H` no es un cambio de variables sino un cambio de
modelo, y lleva el NIS a ~1e11 con `Tr(P)` colapsando a 1e-10 (el filtro declara certeza casi
perfecta mientras está equivocado). La transformación completa sí es exacta pero no mejora el
condicionamiento — `cond(P₀)` queda igual, porque el mal condicionamiento viene de la asimetría
de confianza deliberada de la Sec. 7.3.4, que es información y no escala.

---

## D. Prioridad 3 — Volver a medir, y regla de decisión

Correr `diagnostico.py` sobre telemetría **posterior** a la integración del EMD. Interpretación:

| NIS medio | ρ₁ | Diagnóstico | Acción |
|---|---|---|---|
| ≈ m (=2) | pequeño (<0.05) | Filtro consistente | Desbloquear Fase 2 |
| > m | **grande** (>0.2) | *Model mismatch*: `A` no captura la dinámica | **NO** recalibrar ruido. Ver E. |
| > m | pequeño | Ruido mal calibrado de verdad | Fase 2 (ALS) lo estima |
| < m | pequeño | Filtro conservador | Fase 2 (ALS) lo estima |

### ⚠ Prohibido: constante de inflación ad-hoc sobre `Q` o `R`

Se evaluó y se descartó con números. Sobre el mock actual:

| | NIS | ρ₁ |
|---|---|---|
| `Q` sin tocar | 19.35 | +0.79 |
| `Q` × 7 | 3.22 | **+0.57** |
| `Q` × 50 | 0.53 | **+0.44** |

Inflar `Q` lleva el NIS al objetivo y **deja la autocorrelación casi intacta**. Un error
determinista no lo describe ninguna matriz de covarianza: la inflación no corrige, oculta. Si el
NIS queda alto con ρ₁ pequeño, el estimador correcto es ALS, no una constante a ojo.

---

## E. Aparcado hasta que el bot opere — modelo de velocidad constante

**No implementar en esta fase. Se documenta para que la decisión quede registrada, no para
ejecutarla.**

Si tras integrar el EMD el Ljung-Box sigue rechazando con ρ₁ grande, la causa sería que `A`
asume **velocidad constante** mientras la tesis del sistema (Canción del Micelio) declara que el
mercado es **oscilatorio**, y el Hilo Lento ya estima su frecuencia dominante `ω_m`. La dirección
natural sería un modelo de oscilador armónico alrededor de `S_ref` alimentado por `ω_m`, que
atacaría el mismatch en su origen en vez de ensanchar la incertidumbre para tolerarlo.

Es un cambio de diseño y exige datos con estructura real, no mocks. **Decisión diferida hasta que
el bot opere en cuenta demo.** Mientras tanto, `A` se queda como está.

---

## F. Fase 2 — ALS (bloqueada, sin cambios respecto a v1.1)

Compuerta de entrada: **Ljung-Box debe pasar, leído por magnitud de ρ₁, no solo por p-valor.**

ALS es válido sólo si `A` y `H` son correctas y las covarianzas son las únicas incógnitas. Con
mismatch, ALS absorbe el error de modelo dentro de `Q` y devuelve una respuesta confiadamente
equivocada.

Al implementarla, mantener las tres advertencias de la v1.1:

1. Verificar la condición de unicidad (rango de `𝒜`) antes de confiar en el resultado. Hay
   resultados publicados de no-identificabilidad con entradas desconocidas, y `R_n` aumentado es
   efectivamente un canal de ese tipo. Reportarlo explícitamente si es rank-deficient.
2. ALS calibra `Q_base`, **no** `ρ_k`. La modulación del Micelio se mantiene intacta; que no
   compitan por la misma escala.
3. Empezar sin SDP (ALS + proyección de valores propios negativos) y medir cuánto corrige la
   proyección antes de añadir `cvxpy`.

ALS-IRLS (Huber) sólo si Shapiro-Wilk lo justifica, y **después** de tener el baseline de ALS
estándar.

---

## G. Fase 3 — `Ω_crit` (bloqueada, necesita Testnet)

Sin cambios respecto a v1.1. Recordatorio: `κ` y `μ` se derivan de `Ω_crit`, así que la búsqueda
de orden cero calibra los tres a la vez y el criterio debe evaluar el sistema completo (costo de
ejecución + aversión de inventario + cajas), no sólo el punto de disparo.

Ya está verificado que la cascada funciona: con `Ω_crit` = 1, 2 o 4, `R(Ω_crit)` da 101× `R_base`
en los tres casos.

---

## H. NO HACER AÚN (sin cambios)

- **CMA-ES / meta-optimización global** — bloqueada por falta de datos reales. Optimizar contra
  mocks sintéticos ajusta hiperparámetros al generador. (Nota: *AutoQuant* usa optimización
  bayesiana, no evolutiva, y excluye impacto de mercado — no es template de este sistema.)
- **IMPC / Control Óptimo Inverso** — estructuralmente inaplicable: no hay trayectorias expertas
  que invertir.
- **DeepONet / transporte óptimo para `c²_vol`** — contradice el modelo propio, que ya deriva
  `c²_vol = k·ω_m·ν` con consistencia dimensional verificada.
- **PBO/CSCV, DSR, Walk-Forward, Stationary Bootstrap** — obligatorios, pero son la compuerta de
  salida de la meta-optimización. No hay nada que deflactar todavía.

---

## I. Pendientes heredados

- Sin validar contra CUDA ni acados. Los núcleos son referencia CPU con el mismo layout.
- Contradicciones del PDF #1 (L1 vs cuadrática) y #2 (dimensiones del NMPC) siguen abiertas.
- Marcar en el PDF las Secs. 3.5 y 8.6.1 como superadas por ALS.
- Constantes de ruido (`σ_OU`, `θ_OU`, `μ_OU`, `λ_min`, `η`, `β`) **sin tocar** hasta datos de
  cuenta demo. La guarda de `σ_OU` advierte sin corregir.
- Volcado de telemetría a `.npy`; parquet requiere `pyarrow`.

---

## J. Convenciones

Las de `CLAUDE.md`. Las tres que más aplican aquí:

- `# DIVERGE DEL PDF (Sec X.Y):` cuando una decisión contradiga el documento de diseño.
- `# NOTA DE INTERPRETACION:` cuando se rellene un hueco del PDF.
- Constantes de acoplamiento **derivadas** en `constantes_micelio.py`, nunca literales.
