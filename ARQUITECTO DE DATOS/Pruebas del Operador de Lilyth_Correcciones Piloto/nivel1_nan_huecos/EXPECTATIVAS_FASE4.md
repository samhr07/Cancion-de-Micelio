# Expectativas Selladas — Fase 4 (Fusión Limpia): Puente Histórico de Itô, Llave `Z_hueco` Derivada e Imputación Múltiple de Rubin (`ark/lilyth/`)

**Fecha de sello:** 2026-11-18  
**Autores:** Arquitecto de Datos (Argon) + Auditoría Externa + Revisión de Código (Opus 5.5) bajo orden del Comandante (Samuel Hoyos R.).  
**Regla de inmutabilidad:** Este documento se sella con SHA-256 en `SELLO_EXPECTATIVAS_FASE4.txt` **antes** de compilar y ejecutar el arnés ciego `test_fase4_ciego_rubin.cpp` y `mutaciones_fase4_puente.py`. Ningún umbral, fórmula, tamaño muestral ($R, M$), ni bloque de semillas puede modificarse después del sello. Todo criterio que falle sobre el bloque ciego se reporta textualmente como `[FALLA]` sin ajustar umbrales ni cambiar semillas.

---

## 1. Bloque Ciego de Semillas, Parámetros Derivados y Reglas de Combinación

1. **Separación estricta calibración vs. confirmación ciega:**
   - Quedan **prohibidas** en confirmación las semillas usadas en desarrollo/calibración (`seed = 202`, `seed = 303`, `202 + 17*s`, y `20261117ULL + 10007ULL * rep`).
   - **Bloque ciego `CF4.P2a` / `CF4.P2b` / `CF4.P2c` (`NIVEL M3` + `PurgaDiurna`, $\phi = 0.65$, $n_{\text{op}} = 6000$):** semilla `20261119ULL`.
   - **Bloque ciego `CF4.P3` (`M2` ARFIMA, $d = 0.25, \phi_1 = 0.40, n_{\text{op}} = 6000$, `forzar_via_b = true`):** $16$ realizaciones independientes con semillas `20261119ULL + 31337ULL * s` ($s = 0, \dots, 15$).
   - **Bloque ciego `CF4.P5` (`M1_H0`, `M2_H0`, `M3_H0`, `M1_H1`):** $R = 200$ réplicas por mundo con `ParametrosReplicaH0H1::semilla_base = 20261119ULL` e índices `rep = 0..199` (`20261119ULL + 10007ULL * rep`, y `rep + 50000` interno en `M3_H0`).
   - **Bloque ciego `CF4.P7` (`FLUJO` con deriva $\mu = +5.0\,\sigma_e$):** semilla `20261119ULL ^ 0xF00DCAFEULL`.

2. **Derivación de la llave de amplitud `z_max_efectivo` desde $\alpha_{\text{hueco}}$ (cero literales a ojo):**
   - En `ConfiguracionPuenteHistoricoIto` y `ConfiguracionLilyth`, el umbral por omisión (`z_max_hueco = NaN`) se deriva automáticamente de $\alpha_{\text{hueco}} = 0.01$:
     $$z_{\max}(\alpha_{\text{hueco}}) = \Phi^{-1}\!\left(1 - \frac{\alpha_{\text{hueco}}}{2}\right) = \Phi^{-1}(0.995) = 2.5758293.$$
   - Salto mínimo detectable con potencia $1 - \beta = 0.80$ en un hueco de longitud $r$ ($\Delta t = r + 1$):
     $$\Delta_{\min}(r; \alpha_{\text{hueco}}, 0.80) = \left(z_{\max}(\alpha_{\text{hueco}}) + \Phi^{-1}(0.80)\right)\frac{\sigma_{\Delta t}}{\sigma_e}.$$

3. **Contrato explícito de `valor_z0_imputado()` en `Z0` (`ObservadorSaltos`):**
   - Si el hueco **no** pasa la doble llave ($r > r_{\max}$ o $|Z_{\text{hueco}}| > z_{\max}$): `NaN` en `Z0` y `Z1..Z3`.
   - En **`NIVEL`**, si reconstruye: inyecta `res.x_rng[m][k]` en `Z0` y `Z1` preservando continuidad de nivel sin teletransporte (`+Delta`).
   - En **`FLUJO`**, si reconstruye: mantiene `Z0` causal usando la media estacionaria calibrada `mod_calma.mu` (**nunca `0.0` literal**) para no desplazar `t_cambio` hacia dentro del hueco ni crear un escalón falso cuando $\mu \neq 0$, reservando la trayectoria estocástica de Itô para `Z1b` + Rubin + blindaje `Z2`/`Z3`.

4. **Imputación Múltiple por Reglas de Rubin (1987) + Ajuste de Holm (`ark::holm`):**
   - Número de imputaciones preregistrado: **$M = 16$** en `CF4.P2a` y `CF4.P3`; **$M = 5$** en el Monte Carlo masivo `CF4.P5` ($99.0\%$ de eficiencia relativa bajo $\gamma_M \le 0.05$).
   - Para cada unidad evaluada (una serie individual, nunca promediando series antes de aplicar Rubin):
     $$\bar{S}_M = \frac{1}{M}\sum_{m=0}^{M-1} S^{(m)}, \qquad \bar{W}_M = \frac{1}{M}\sum_{m=0}^{M-1} \left(\text{ee}_S^{(m)}\right)^2, \qquad B_M = \frac{1}{M-1}\sum_{m=0}^{M-1} \left(S^{(m)} - \bar{S}_M\right)^2,$$
     $$T_M = \bar{W}_M + \left(1 + \frac{1}{M}\right) B_M, \qquad \text{ee}_{\text{Rubin}} = \sqrt{T_M}, \qquad \gamma_M = \frac{(1 + 1/M) B_M}{T_M},$$
     $$\nu_M = (M-1)\left(1 + \frac{\bar{W}_M}{(1 + 1/M)B_M}\right)^2, \qquad t_M = \frac{\bar{S}_M}{\text{ee}_{\text{Rubin}}}, \qquad p_{\text{Rubin}} = P(|T_{\nu_M}| \ge |t_M|).$$
   - En `CF4.P5`, el rechazo de `B4_RUBIN` inserta $p_{\text{Rubin}}$ en el índice `1` (contraste Onsager) de la familia armónica de 4 pruebas (tomando la mediana entre las $M$ imputaciones para los índices auxiliares `0, 2, 3`) y aplica **`ark::holm()`** (`p_holm_rubin[1] < 0.05`).

---

## 2. Controles Falsables Sellados (`CF4.P1` – `CF4.P7`)

| Control | Escenario y Diseño | Criterio de Aceptación Sellado (Estricto) |
|---|---|---|
| **`CF4.P1`** | No regresión Nivel 1 con políticas piloto en modo base (`usar_protocolo_huecos = false`) en `M1`, `M2`, `M3` (limpios y con huecos). | Volcado hexadecimal `%a` real generado por `ObservadorLilyth` **idéntico bit a bit** a `volcado_orig.txt` (`0` diferencias). |
| **`CF4.P2a`** | `NIVEL M3` + `PurgaDiurna` (semilla ciega `20261119ULL`), micro-cortes limpios $r = 4$ fuera de choques reales; tabla comparativa completa `B0`, `B1`, `B1b`, `B2_LINEAL`, `B3`, `B4_RUBIN` ($M = 16$). | $\|\bar{S}_M - S_0\| / \text{ee}_{\text{Rubin}} \le 0.25$; $\text{Alm\_Tmp\_Esp}(\text{B4}) \le \text{Alm\_Tmp\_Esp}(\text{B2})$; $\max_m \text{err\_anc}^{(m)} \le 10^{-10}$. |
| **`CF4.P2b`** | Matriz completa de choque oculto en `NIVEL M3` (semilla ciega `20261119ULL`): $\Delta/\sigma_e \in \{4, 6, 8, 10, 12\} \times r \in \{1, 2, 3, 4\}$ con historia causal real (`hist = 64`), evaluada en $K = 20$ posiciones tranquilas independientes. | 1) Reporte íntegro de la matriz $(r, \Delta_{\min}(80\%), \bar{Z}_{\text{hueco}}, \text{tasa\_bloqueo})$.<br>2) En toda celda con $\Delta/\sigma_e \ge \Delta_{\min}(r; 0.01, 0.80)$ la tasa de bloqueo es $\ge 80\%$ (y $100\%$ en $+12\,\sigma_e, r \in \{1, 2\}$ con `FLAG_CHOQUE_EN_HUECO` activo y `reconstruido == false`). |
| **`CF4.P2c`** | Discriminación efectiva de `B4` frente a `B2_LINEAL` (calmas inventadas) y frente a `B3_ITO_SIN_FILTRO` (choques ocultos) en `NIVEL M3` (semilla ciega `20261119ULL`):<br>(i) $K = 24$ micro-cortes limpios $r = 4$: razón de varianza **intra-hueco** $\text{RazónVar} = \frac{\overline{\text{Var}}_{s=0..r}(\Delta X^{\text{método}}_s)}{\overline{\text{Var}}_{s=0..r}(\Delta X^{\text{física}}_s)}$ + lectura informativa del sesgo continuo en `sd.temp[0].s2` de `ObservadorSaltos`.<br>(ii) $K = 10$ micro-cortes $r = 2$ con choques ocultos $+10\,\sigma_e$ pasando por `ObservadorSaltos` + `ObservadorLilyth`. | 1) En micro-cortes limpios: $\text{RazónVar}(\text{B2\_LINEAL}) < 0.10$ mientras $\text{RazónVar}(\text{B4\_ITO\_RNG}) \in [0.60, 1.40]$.<br>2) En los $10$ micro-cortes con choque oculto $+10\,\sigma_e$: `B3_ITO_SIN_FILTRO` bloquea $0/10$ (reconstruye todos), mientras `B4` bloquea $\ge 8/10$ ($\ge 80\%$) con `FLAG_CHOQUE_EN_HUECO`. |
| **`CF4.P3`** | `M2` (`ARFIMA`, $16$ semillas ciegas, `forzar_via_b = true`, `INTERCAMBIABLE`, `Z1b` en `t_cambio - 2`, $r = 3$, Rubin $M = 16$ calculado **por serie**). Tabla completa `B0`, `B0b_CORTE_SUFIJO`, `B1`, `B2`, `B4_MEDIA`, `B4_RUBIN16`. | 1) $\text{RMSE}_{\text{rel}}(e_{\text{iny}})$ real vía `ark::lilyth::innovaciones()` para `B4_RUBIN16` $\le 0.15$ y $\le \text{RMSE}(\text{B0b\_CORTE\_SUFIJO})$.<br>2) $\|\frac{1}{16}\sum_s (\bar{S}_{M,s} - S_{0,s})\| / \overline{\text{ee}}_{\text{Rubin}} \le 0.25$.<br>3) Coherencia Rubin: $10^{-6} < \overline{\gamma}_M \le 0.25$ y $\text{err\_anc} \le 10^{-10}$.<br>4) Reporte honesto de `B1` frente a `B4` sin afirmar ventaja inexistente. |
| **`CF4.P4`** | Máquina de pausa dura en corte largo $r = 12 > r_{\max} = 3$. | Exactamente $9$ pasos en `EN_PAUSA_DURA` (`FLAG_PAUSA_UMBRAL`) y `cortar_segmento_al_reenganchar == true` al reconectar. |
| **`CF4.P5`** | Monte Carlo ciego $R = 200$ (`q = 2%`, `r = 2`) en `M1_H0`, `M2_H0`, `M3_H0` y `M1_H1`. Tabla comparativa `LIMPIO q0`, `B0`, `B1`, `B2`, `B3`, `B4_RUBIN` ($M = 5$ + `ark::holm()`). | 1) En `M1_H0`, `M2_H0` y `M3_H0`: límite inferior de Wilson al $95\%$ del rechazo de `B4_RUBIN` (`p_holm_rubin[1] < 0.05`) $\le \alpha = 0.05$.<br>2) Retención de potencia en `M1_H1` $\ge 0.85$.<br>3) $\|\bar{S}_M - \bar{S}_0\| / \overline{\text{ee}}_{\text{Rubin}} \le 0.25$ en `M1_H0`, `M2_H0` y `M3_H0`.<br>4) $\Delta\text{Censura}(\text{M1\_H0}) \le p_{\text{teor}}(\bar{l}) + 1.96\sqrt{\frac{p_{\text{teor}}(1-p_{\text{teor}})}{R\cdot \bar{K}}}$, con $p_{\text{teor}}(\bar{l}) = 1 - (1 - q/r)^{\bar{l} + r - 1}$ derivado del ancho medio $\bar{l}$ real de `Z2` en `LIMPIO`. |
| **`CF4.P6`** | Blindaje obligatorio de `FLAG_PUENTE_ITO` en `Z2` y `Z3` dentro de `ObservadorLilyth`. | **Ambos** (`veto_z2 && veto_z3`): en `Z2` activa `perfil_corrupto_z2_cruda && sin_nulo_b && !evaluado`; en `Z3` deja `eta[0] = NaN` y conserva `tau_rest_consolidado` real ($< 10^{-12}$). |
| **`CF4.P7`** | Deriva no nula en `FLUJO` ($\mu = +5.0\,\sigma_e$, AR(1) $\phi = 0.65$, semilla ciega) con $12$ micro-cortes limpios $r = 2$ alimentando `ObservadorSaltos` vía `valor_z0_imputado(..., mod.mu)`. | `0` alarmas espurias de salto negativo en los huecos (frente a disparo inmediato si `valor_z0_imputado` inyecta `0.0` literal). |

---

## 3. Mutantes Obligatorios Sellados (`mutaciones_fase4_puente.py`)

El script debe compilar con `/W4 /permissive-`, pasar el testigo y cazar **5 de 5** mutantes reales con $R = 200$ fijo desde este preregistro (cada patrón con `count == 1` exacto):

1. **`MF4.0` (Testigo limpio):** sin cambios $\rightarrow$ `TESTIGO PASA (0 fallas)`.
2. **`MN1.7` (`observador_lilyth.hpp`):** anular `veto_z3_ito` y `veto_ito_z2` $\rightarrow$ cazado por **`CF4.P6`**.
3. **`MN1.8` (`puente_historico_ito.hpp`):** anular `if (std::abs(out.z_hueco) > out.z_max_efectivo)` $\rightarrow$ cazado por **`CF4.P2b`** y **`CF4.P2c`**.
4. **`MN1.9` (`puente_historico_ito.hpp`):** colapsar `out.b_entre = 0.0;` en `combinar_rubin()` $\rightarrow$ cazado por `gamma_M > 1e-6` en **`CF4.P3`**.
5. **`MN1.10` (`puente_historico_ito.hpp`):** sustituir `reposo_flujo` por `0.0` literal en `valor_z0_imputado()` $\rightarrow$ cazado por **`CF4.P7`**.
