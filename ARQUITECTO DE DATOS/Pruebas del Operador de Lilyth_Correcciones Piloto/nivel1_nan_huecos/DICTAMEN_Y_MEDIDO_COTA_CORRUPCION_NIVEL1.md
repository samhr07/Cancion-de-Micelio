# Dictamen y Medido Empírico Final — Nivel 1: Protocolo de Inputs y Cotas de Corrupción (`ObservadorLilyth`)

**Fecha de cierre:** 2026-11-17  
**Expectativas selladas:** `EXPECTATIVAS_COTA_NIVEL1.md` (`sha256: 54aa457e9dc69ddec1a982f4904e2d71e354f4ae9ea2dada4aba77c7f52eafa3`, `SELLO_EXPECTATIVAS.txt`) + anexo declarado `ANEXO_DESVIOS_EXPECTATIVAS_NIVEL1.md` (`D1`–`D5`).  
**Regímenes comparados:** **Versión Estricta (`ESTRICTO`, $\alpha = 0.05$)** y **Versión Dulce (`DULCE`, error $< 50\%$)**.  
**Regla de reporte:** Ninguna cota se extrapola más allá de la malla efectivamente corrida; todas las cotas $(q_{\max}, r_{\max})$ se calculan con parada monótona estricta (`break` al primer fallo).

---

## 1. Verificación de las 6 Condiciones Selladas (`CN1.1` a `CN1.6`)

| Condición | Qué verifica | Resultado medido (`test_cota_corrupcion_nivel1.log` / `mutaciones_nivel1.log`) |
|---|---|---|
| **`CN1.1`** | Identidad bit a bit en formato hexadecimal IEEE-754 (`%a`) en `M1`, `M2` y `M3` (limpios y con huecos) contra `observador_lilyth_original.hpp` (`volcado_identidad_nivel1.cpp`) | **`[PASA]`**: `114` líneas hexadecimales idénticas byte a byte (`SHA256 = 8BCED50D6877D3832D1481C3F28ECCECD779C766C577474900A337220DB208F9`, `IDENTICOS = True`) |
| **`CN1.2`** | Tamaño nominal bajo $H_0$ en $q = 0$ ($R = 200$ réplicas, $a = 8.0\,\sigma_\varepsilon$, `CUANTIL_RESOLUCION` + `DERIVADO`, Wilson inferior $\le 0.05$) y potencia base bajo $H_1$ ($z_\tau = 4.0\,\sigma_\varepsilon\,0.97^\tau \ge 50\%$) | **`[PASA]`**:<br>• `M1_H0`: `4/200 = 0.0200 [0.0078, 0.0503]`<br>• `M2_H0`: `7/199 = 0.0352 [0.0171, 0.0708]`<br>• `M3_H0`: `7/200 = 0.0350 [0.0171, 0.0705]`<br>• `M1_H1` (Potencia): **`123/200 = 0.6150 [0.5460, 0.6797]`** ($61.5\% \ge 50\%$) |
| **`CN1.3`** | Supervivencia en `M3` (`NIVEL + PurgaDiurna`) bajo $q = 2\%$ de `MERCADO_QUIETO` (`BarraCausal`) y corte limpio en `FRONTERA_ESTRUCTURAL` (`FIN_DE_SEMANA`) | **`[PASA]`**: `PROTOCOLO` tiene `0` cortes (`dCens = 0.00` en trayectoria; `200/200` series válidas en $R=200$), mientras `ACTUAL` pierde `100%` en trayectoria fija (`54/200` series válidas en $R=200$, `dCens = 0.96`) |
| **`CN1.4`** | Escenario con $25\%$ de choques con `Z2` corrupto ($r=3$) y $50\%$ con corte en `Z3` ($r=1$) sobre $R = 200$ réplicas $H_0$ | **`[PASA]`**: `eval_z2 = 0`, `eval_z3 = 0`, y rechazo de Holm bajo $H_0$ con **`Wilson lo = 0.002954 <= 0.05`** |
| **`CN1.5`** | Derivación algorítmica monótona de cotas e invariantes estructurales a priori ($\text{cota}_{\text{ESTRICTO}} \le \text{cota}_{\text{DULCE}}$ y $\|\text{sesgo}_{\text{reloj}}\| < \|\text{sesgo}_{\text{comprime}}\|$ en `Z3a`) | **`[PASA]`** |
| **`CN1.6`** | Caza de mutantes reales por reemplazo de texto en copia de cabecera (`mutaciones_nivel1.py`, condición $\ge 5/6$) | **`[PASA]`**: **6 de 6 mutantes cazados** (`MN1.0` testigo pasa con `0` fallas; `MN1.1` a `MN1.6` todos cazados) |

---

## 2. Detalle de la Caza de Mutantes Reales de Cabecera (`mutaciones_nivel1.log`)

| ID | Archivo mutado | Cambio exacto en la cabecera C++ | Estado | Condición que lo caza y lectura medida |
|---|---|---|---|---|
| **`MN1.0`** | — | Testigo sin cambios | **PASA** | `TESTIGO PASA (0 fallas)` |
| **`MN1.1`** | `protocolo_inputs_lilyth.hpp` | `MERCADO_QUIETO` en `NIVEL` puesto como `valido_para_operador = false` | **CAZADO** | `CN1.3` |
| **`MN1.2`** | `observador_lilyth.hpp` | `tau_ef = rama.pasos_validos_z3` en `RELOJ_FISICO_TOPE_R` (comprime $\tau$) | **CAZADO** | `CN1.5` (sesgo negativo en $\tau_{\text{rest}}$ de `Z3a`) |
| **`MN1.3`** | `observador_lilyth.hpp` | No excluir `CENSURADA_DATO_INVALIDO` en `evaluar_nulos_representacion` | **CAZADO** | `CN1.1` y `CN1.4` (`eval_z3 = 134` ranuras inválidas metidas al estrato e inflación estadística del falso rechazo bajo $H_0$ hasta **`Wilson lo = 0.170033 > 0.05`**) |
| **`MN1.4`** | `protocolo_inputs_lilyth.hpp` | `FRONTERA_ESTRUCTURAL` puesto como `valido_para_operador = true` | **CAZADO** | `CN1.3` |
| **`MN1.5`** | `observador_lilyth.hpp` | Omisión invertida en `PoliticaZ3Piloto` (`RELOJ_FISICO_TOPE_R` y `r_max_z3_piloto = 3` por omisión) | **CAZADO** | `CN1.1` (rompe el volcado hexadecimal bit a bit contra `observador_lilyth_original.hpp`) |
| **`MN1.6`** | `observador_lilyth.hpp` | Ignorar `if (z2_corrupto)` en la vía `(b)` (`evaluar_una_referencia`) | **CAZADO** | `CN1.4` (`eval_z2 = 31` perfiles corruptos inyectados en vía `(b)`) |

---

## 3. Resultados Medidos por Zona y Experimento

### 3.1. Zona `Z1` (`hist_prev_est`, $t < t_{\text{cambio}}$) con `ark::lilyth::innovaciones()` y `ejecutar_via_b_evento()` Reales (`test_punto2_z1.log` y `test_puntos3_4_operador.log`)

- **Borde inmediato `Z1a` ($\Delta t = 1$, hueco termina en $t_{\text{cambio}} - 1$, malla medida $r \in \{1, 2, 3, 5\}$ en `innovaciones()` y $r \in \{1, 2, 3\}$ en vía `(b)`):**
  - **Versión Estricta ($\le 0.15\sigma_\varepsilon$):** **FALLA en las 3 políticas** (`RMSE = 0.4345–0.4707` en `M1`, `0.2742–0.3858` en `M2`). En el primer paso del perfil (`e_iny[0]`), el error es de `1.21` a `2.35` $\sigma_\varepsilon$. Por tanto, **$r_{\max, Z1a}^{\text{ESTRICTO}} = 0$**.
  - **Versión Dulce ($< 0.50\sigma_\varepsilon$):** En `M1`, `COMPRESION_BASE` falla desde $r = 2$ (`0.5538`), mientras que **`PUENTE_CAUSAL` pasa en toda la malla medida** (`r <= 3` en vía `(b)` completa con $\text{RMSE} = 0.3274\text{–}0.3317$ y $|dS|/\text{ee}_S \le 0.2188$; `r <= 5` en `innovaciones()` aislado con $\text{RMSE} = 0.4387$ en `M1` y `0.2907` en `M2`). Por tanto, en la malla de vía `(b)`, **$r_{\max, Z1a}^{\text{DULCE}} = 3$**.
- **Interior `Z1b` ($\Delta t \ge 2$, hueco termina en $\le t_{\text{cambio}} - 2$):**
  - En `M1` (`AR(1)`, $p = 1$): las 3 políticas dan **`0.0000` exacto** para todo $r \in \{1, 2, 3, 5\}$.
  - En `M2` (`ARFIMA_WHITTLE`, $\hat{d} = 0.3133$): en $\Delta t = 2$, **`CORTE_SUFIJO` falla `ESTRICTO`** en vía `(b)` porque amputar la cola fraccionaria previa desvía el estadístico $S$ en **`0.2701` $\text{ee}_S > 0.25$** (`RMSE = 0.1286–0.1358`). En cambio, **`PUENTE_CAUSAL` pasa `ESTRICTO` y `DULCE` en toda la malla medida** (`r <= 3` en vía `(b)` con $\text{RMSE} = 0.0099\text{–}0.0199$, $|dS|/\text{ee}_S \le 0.2188$; `r <= 5` en `innovaciones()` aislado con $\text{RMSE} = 0.0577$).

### 3.2. Zona `Z2` (`perfil_est` `[t_cambio, t_apertura]`) con `ejecutar_via_b_evento()` Real (`test_puntos3_4_operador.log`)
- Para cualquier $r \in \{1, 2, 3\}$, tanto `COMPRESION_BASE` como `PUENTE_CAUSAL` **fallan `ESTRICTO` y `DULCE`** ($\text{RMSE}(e_{\text{iny}})/\sigma_\varepsilon = \mathbf{1.1992\text{–}1.2817} > 0.50$).
- **Cota en `Z2`:** **$r_{\max, Z2}^{\text{ESTRICTO}} = 0$ y $r_{\max, Z2}^{\text{DULCE}} = 0$**. Si una ranura necesita la vía `(b)` y tiene $r \ge 1$ inválido en `[t_cambio, t_apertura]`, se bloquea la vía `(b)` (`sin_nulo_b = true`), sin cerrar la ranura cuando la vía `(a)` tiene excedencias suficientes para evaluarla.

### 3.3. Zona `Z3` (`ObservadorLilyth::avanzar()` real en $t > t_{\text{apertura}}$, malla $r \in \{1, 2, 3, 5\}$, `test_puntos3_4_operador.log`)
- **Subcaso `Z3a` (Hueco en $\tau = 1\dots r$ antes del cruce de $\eta \le 0.5$):**
  - `CORTA_AL_PRIMERO` (actual): falla `ESTRICTO` y `DULCE` desde $r = 1$ (`dCens = 0.76`, $|dS|/\text{ee}_S = 7.7135$).
  - `MUTANTE_COMPRIME_TAU`: falla `ESTRICTO` y `DULCE` en $r \in \{1, 2, 3\}$ (`Sesgo(tau_r) = -0.7368` a `-0.9474`, $|dS|/\text{ee}_S = 1.7506$ a $3.5012$).
  - `RELOJ_FISICO_SIN_COMPR`: **pasa `ESTRICTO` y `DULCE` en toda la malla medida** ($r \in \{1, 2, 3, 5\}$, con `dCens = 0.00`, `Sesgo(tau_r) = 0.0000`, $|dS|/\text{ee}_S \le 0.1094$).
- **Subcaso `Z3b` (Adversarial: hueco exactamente encima de $\tau = \tau_{\text{rest}}$ en el $100\%$ de los eventos):**
  - Falla `ESTRICTO` y `DULCE` desde $r = 1$ (`Sesgo(tau_r) = +1.1818` pasos, $|dS|/\text{ee}_S = 1.5865$), por lo que **$r_{\max, Z3b}^{\text{ESTRICTO}} = r_{\max, Z3b}^{\text{DULCE}} = 0$**.

---

### 3.4. Flujo Completo en $R = 200$ Réplicas con Cota Monótona (`test_punto6_wilson_h0_h1.log`)

#### A. Mercado Quieto (`ILIQUIDEZ_TRANSITORIA`, `[6A]`)
1. **En `M3_H0` (`NIVEL + PurgaDiurna`, donde la barra quieta mantiene $\Delta x_t = 0$):**
   - `ACTUAL` (`valido = !b.es_imputado`): pierde el $82\%$ de los eventos en $q = 1\%$ (`160/200` series válidas), el $96\%$ en $q = 2\%$ (`54/200`) y el **$100\%$ (`0/200` series válidas)** en $q \ge 5\%$. Falla `ESTRICTO` y `DULCE` desde $q = 1\%$.
   - `PROTOCOLO` (`MERCADO_QUIETO -> valido = true`): conserva **`200/200` series válidas en todas las tasas $q \in \{1\%, 2\%, 5\%, 10\%\}$**, con `dCens` de `0.01` a `0.07` ($\le 0.10$) y falso rechazo bajo $H_0$ entre `0.015` y `0.045` (`Wilson lo <= 0.024 <= 0.05`). **En `NIVEL` (`M3_H0`), pasa `ESTRICTO` y `DULCE` hasta $q_{\max} = 10.0\%$**.
2. **En `M1_H1` (`FLUJO`, donde inyectar $x_t = 0$ reemplaza físicamente el flujo activo y amortigua la cola lenta $z_\tau$):**
   - La retención de potencia `RetPot_H1` cae monótonamente con $q$: `0.85` en $q=1.0\%$ (`0.8455 < 0.85`, por lo que falla por $0.0045$ el corte `ESTRICTO` de `0.85`), `0.76` en $q=2.0\%$ (pasa `DULCE` $\ge 0.50$), `0.48` en $q=5.0\%$ (falla `DULCE` $< 0.50$) y `0.27` en $q=10.0\%$.
   - Por tanto, si se fuerza $x_t = 0$ sobre una serie de `FLUJO` con choques lentos, la cota monótona combinada es **$q_{\max}^{\text{ESTRICTO}} = 0.0\%$ ($1.0\%$ al $84.6\%$ de potencia) y $q_{\max}^{\text{DULCE}} = 2.0\%$**.

#### B. Corrupción Física Extremo a Extremo (`ObservadorSaltos` + `ObservadorLilyth` ven ambos los `NaN`, `[6B.1]`, $R = 200$)
| $q$ (%) | Racha $r$ | Política | `dCens` | $\|\Delta\bar{S}\|/\overline{\text{ee}}_S$ | Falso rechazo $H_0$ [Wilson 95 %] | `RetPot_H1` | **`ESTRICTO`** | **`DULCE`** |
|---:|---:|---|---:|---:|---|---:|---|---|
| **0.5%** | 1, 2, 3 | `ACTUAL` | 0.04–0.06 | 0.0007–0.0137 | `0.020–0.030 [lo <= 0.014]` | 1.00–1.07 | **PASA** | **PASA** |
| **0.5%** | 1, 2, 3 | `PROTEGIDO` | **0.03–0.06** | 0.0121–0.0313 | `0.020 [lo = 0.008]` | 0.99–1.10 | **PASA** | **PASA** |
| **1.0%** | 1 | `ACTUAL` / `PROTEGIDO` | 0.11 / 0.10 | 0.0069 / 0.0357 | `0.045 / 0.030` | 0.98 / 1.08 | **FALLA** *(dCens > 0.10)* | **PASA** |
| **1.0%** | 2, 3 | `PROTEGIDO` | **0.04** | 0.0044–0.0187 | `0.035 [lo = 0.017]` | 0.95–1.02 | **PASA** | **PASA** |
| **2.0%** | 1 | `ACTUAL` / `PROTEGIDO` | 0.18 / 0.16 | 0.0078 / 0.0575 | `0.040 / 0.035` | 0.97 / 1.17 | **FALLA** *(dCens > 0.10)* | **PASA** |
| **2.0%** | 2 | `ACTUAL` | **0.10** *(0.104)* | 0.0405 | `0.010 [lo = 0.003]` | 1.00 | **FALLA** | **PASA** |
| **2.0%** | 2 | `PROTEGIDO` | **0.09** | 0.0575 | `0.010 [lo = 0.003]` | 1.09 | **PASA** | **PASA** |
| **5.0%** | 1, 2, 3 | `PROTEGIDO` | **0.16–0.29** | 0.0448–0.2371 | `0.020–0.040 [lo <= 0.020]` | 0.92–1.25 | **FALLA** *(dCens > 0.10)* | **PASA** |

#### C. Corrupción Física con `forzar_via_b = true` (`Z1+Z2+Z3` activas en el $100\%$ de los eventos, `[6B.2]`, $R = 200$)
- En **`DULCE`**, `PROTEGIDO` **pasa en todas las celdas medidas hasta $q = 5.0\%$ ($r \in \{1, 2\}$)** (`dCens = 0.05–0.34 < 0.50`, `|dS|/ee <= 0.3244 <= 1.00`, `Wilson lo <= 0.014 <= 0.10`, `RetPot_H1 >= 0.94`).
- En **`ESTRICTO`**, al forzar que el $100\%$ de los eventos usen vía `(b)` y rechazar cualquier evento cuyo perfil `Z2` tenga un paso corrupto, `dCens` supera `0.10` en $r = 1$ desde $q = 1.0\%$ (`dCens = 0.11`), mientras que para $r = 2$ pasa `ESTRICTO` hasta **$q = 2.0\%$** (`dCens = 0.05` en $q=1\%$, `0.10` en $q=2\%$).

---

## 4. Tabla Resumen Final de Cotas Monótonas Medidas (Sin Extrapolaciones)

| Ámbito / Zona (Malla efectivamente corrida) | **Cota Monótona `ESTRICTO`** | **Cota Monótona `DULCE`** |
|---|---|---|
| **1. `MERCADO_QUIETO` en `NIVEL` (`M3_H0`, malla $q \in \{0.5\%\dots 10\%\}$)** | **$q_{\max} = 10.0\%$** | **$q_{\max} = 10.0\%$** |
| **2. `MERCADO_QUIETO` forzado como $x_t=0$ sobre `FLUJO` activo (`M1_H1`, malla $q \in \{1\%\dots 10\%\}$)** | **$q_{\max} = 0.0\%$** *(en $q=1\%$ `RetPot = 0.846 < 0.85`)* | **$q_{\max} = 2.0\%$** *(en $q=5\%$ `RetPot = 0.48 < 0.50`)* |
| **3. Zona `Z1a` (Borde $t_{\text{cambio}}-1$, vía `(b)`, malla $r \in \{1, 2, 3\}$)** | **$r_{\max} = 0$** | **$r_{\max} = 3$ pasos** *(con `PUENTE_CAUSAL`)* |
| **4. Zona `Z1b` (Interior $t \le t_{\text{cambio}}-2$, vía `(b)`, malla $r \in \{1, 2, 3\}$)** | **$r_{\max} = 3$ pasos** *(con `PUENTE_CAUSAL`; `CORTE_SUFIJO` falla)* | **$r_{\max} = 3$ pasos** *(medido hasta $r=5$ en `innovaciones()` aislado)* |
| **5. Zona `Z2` (Perfil `[t_cambio, t_a]`, vía `(b)`, malla $r \in \{1, 2, 3\}$)** | **$r_{\max} = 0$** | **$r_{\max} = 0$** |
| **6. Zona `Z3a` (`avanzar` fuera de $\tau_{\text{rest}}$, malla $r \in \{1, 2, 3\}$ en unificado / $\{1..5\}$ en `P3_4`)** | **$r_{\max} = 3$ pasos** *(reloj físico sin comprimir, sesgo `0.0000`)* | **$r_{\max} = 3$ pasos** *(pasa también en $r=5$ en `test_puntos3_4`)* |
| **7. Zona `Z3b` (`avanzar` encima de $\tau_{\text{rest}}$ en el $100\%$ de eventos, malla $r \in \{1, 2, 3\}$)** | **$r_{\max} = 0$** | **$r_{\max} = 0$** |
| **8. Flujo Extremo a Extremo con `ObservadorSaltos` viendo `NaN` (`6B.1`, malla $q \in \{0.5\%\dots 5\%\}$)** | **$q_{\max} = 0.5\%$** *(en $r=1$)* / **$2.0\%$** *(en $r=2, 3$)* | **$q_{\max} = 5.0\%$** *(en toda la malla $r \in \{1, 2, 3\}$)* |
