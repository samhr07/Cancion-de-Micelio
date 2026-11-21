# Expectativas Selladas — Fase 4, Adenda v2: Escalado Diurno Causal en `puente_historico_ito.hpp`, Tope `r_max = 2` en `NIVEL` y Caza Limpia 5/5 de Mutantes

**Fecha de sello:** 2026-11-19  
**Predecesor inmutable:** `EXPECTATIVAS_FASE4.md` (SHA-256 `49fe152d00248465489ca2c561647d61033197f2543e819f179bec25e4806ef1`) y sus bitácoras ciegas `test_fase4_ciego_rubin.log` / `mutaciones_fase4_puente.log`, que permanecen intactas como registro histórico de v1.  
**Regla de inmutabilidad:** Este documento se sella con SHA-256 en `SELLO_EXPECTATIVAS_FASE4_V2.txt` **antes** de compilar y ejecutar el arnés ciego `test_fase4_ciego_rubin_v2.cpp` y `mutaciones_fase4_puente_v2.py`. Ningún umbral, fórmula, tamaño muestral ($R, M, K$), ni bloque de semillas puede modificarse después del sello. Todo criterio que falle sobre el tercer bloque ciego se reporta textualmente como `[FALLA]` sin ajustes post-hoc.

---

## 1. Motivación Exacta de la Adenda v2 (Derivada de la Corrida Ciega v1)

La corrida ciega v1 (`20261119ULL`) aprobó `CF4.P1`, `CF4.P2a`, `CF4.P2c`, `CF4.P3`, `CF4.P4`, `CF4.P5`, `CF4.P6` y `CF4.P7`, y dejó exactamente dos hallazgos científicos abiertos que esta Adenda v2 resuelve de raíz:

1. **Hallazgo 1 (`CF4.P2b` en v1 — Heterocedasticidad diurna en `NIVEL + PurgaDiurna` + tope operativo `r_max = 2` + error de diseño estadístico del criterio v1):**
   - *Causa física en el código:* en v1, `evaluar_y_reconstruir_hueco_ito()` asumía incrementos homocedásticos ($\kappa_s \equiv 1$) sobre `M3` (`PurgaDiurna`, donde la amplitud física oscila $g(f) \in [0.5, 1.5]$), y `pos_calmas` muestreaba con paso `240` (múltiplo de `24`, cayendo siempre en las mismas fases ascendentes $\kappa > 1$). Además, la matriz $4 \times 5$ demostró que para $r \ge 3$ la cota $\Delta_{\min}(80\%)$ se degrada a $12.64\,\sigma_e\text{–}15.35\,\sigma_e$, acotando el dominio seguro de reconstrucción en `NIVEL` a **$r_{\max} = 2$**.
   - *Causa estadística en el criterio v1:* incluso para un proceso gaussiano estacionario perfecto sin marea diurna, en $(r = 2, +12\,\sigma_e)$ la potencia teórica exacta es $\Phi(12 / 2.8316 - 2.5758) = 95.17\%$ (no $100\%$), por lo que exigir `20/20 = 100%` producía una tasa de fallo falso del $63\%$ bajo código exacto.
   - *Solución en v2:*
     - `ConfiguracionPuenteHistoricoIto::r_max` pasa a **`2`** por omisión (tope operativo en `NIVEL`; en `ARFIMA Z1b` puede configurarse explícitamente).
     - `evaluar_y_reconstruir_hueco_ito()` y `delta_minima_detectable_ito()` incorporan el parámetro opcional **`escalas_fase`** ($\kappa_0, \dots, \kappa_r$), alimentado estrictamente de forma **causal** desde `EspejoPurgaNulo::factor(fase, dia)` ($\kappa_s = 1 / f_s = \sqrt{v_{\text{fase}}(f_s) / v_{\text{media}}}$, usando solo días previos comprometidos, cero fuga futura):
       $$c_j = \sum_{i=j}^{\Delta t - 1} \kappa_i \psi_{i-j}, \qquad v_{\text{borde}} = \sigma_e^2 \sum_{j=0}^{\Delta t - 1} c_j^2, \qquad w_i = \kappa_i \sigma_e^2 \sum_{j=0}^{i} \psi_{i-j} c_j,$$
       $$\Delta X_{\text{esp}} = \sum_{s=0}^{\Delta t - 1} \kappa_s \mu_s, \qquad Z_{\text{hueco}} = \frac{(X_{\text{fin}} - X_{\text{ini}}) - \Delta X_{\text{esp}}}{\sqrt{v_{\text{borde}}}},$$
       y las realizaciones de Hoffman–Ribak proyectan $\Delta\tilde{X}_s = \kappa_s \tilde{z}_s$ con pesos $w_s$, preservando anclaje exacto ($\sum \Delta X_s^{(m)} = X_{\text{fin}} - X_{\text{ini}}$), homogeneidad $Z_{\text{hueco}} \sim \mathcal{N}(0, 1)$ en las 24 horas del día y varianza local exacta al volver a pasar por `PurgaDiurna`.
     - `CF4.P2b` pasa a muestreo **estratificado uniforme sobre las 24 fases diurnas**, separando calibración bajo $H_0$ (incluyendo razón pico/valle que caza `MN1.11`) y coherencia binomial exacta con $p_{\text{teor}}(r, A)$ en la matriz $4 \times 5$.

2. **Hallazgo 2 (`MN1.10` escapó de `CF4.P7` en v1):**
   - Con $\mu = +5.0\,\sigma_e$, $\phi = 0.65$ ($\tau_{\text{int}} = 4.714$) y $r = 2$, inyectar `0.0` producía $g_{\text{GLR}} = 3.06 < h_{\text{GLR}}(\text{ARL}_0 = 450) \approx 5.8$, por lo que `ObservadorSaltos` no disparaba alarma ni bajo el mutante.
   - *Solución en v2:* `CF4.P7` evalúa una serie `FLUJO` con deriva positiva $\mu = +10.0\,\sigma_e$ (donde inyectar `0.0` durante $r = 2$ produce $g_{\text{GLR}} \approx 12.25 \gg h_{\text{GLR}}$, disparando alarmas negativas masivas en `ObservadorSaltos`) y verifica tanto `0` alarmas negativas espurias como el invariante contractual $|x_{\Delta}(t) - \hat{\mu}_{\text{calma}}| \le 10^{-12}$.

---

## 2. Tercer Bloque Ciego de Semillas (Estrictamente Disjunto de Calibración y de v1)

Quedan prohibidas todas las semillas de calibración (`101, 202, 303, 20261117ULL`), de Nivel 1 (`20261118ULL`) y de Fase 4 v1 (`20261119ULL`). Toda la Adenda v2 corre sobre el tercer bloque ciego **`SEMILLA_CIEGA_V2 = 20261120ULL`**:

- `CF4.P2a`, `CF4.P2b`, `CF4.P2c` (`NIVEL M3` + `PurgaDiurna`, $\phi = 0.65$, $n_{\text{op}} = 6000$): `20261120ULL`.
- `CF4.P3` (`M2` ARFIMA, $d = 0.25, \phi_1 = 0.40$, $n_{\text{op}} = 6000$, `forzar_via_b = true`): $16$ realizaciones `20261120ULL + 31337ULL * s` ($s = 0 \dots 15$).
- `CF4.P5` (`M1_H0`, `M2_H0`, `M3_H0`, `M1_H1`): $R = 200$ réplicas con `ParametrosReplicaH0H1::semilla_base = 20261120ULL` (`rep = 0..199`).
- `CF4.P6` (blindaje `Z2 && Z3`): `20261120ULL + 7777ULL`.
- `CF4.P7` (`FLUJO` con deriva $\mu = +10.0\,\sigma_e$, AR(1) $\phi = 0.65$): `20261120ULL ^ 0xF00DCAFEULL`.

---

## 3. Controles Falsables Sellados (`CF4.P1` – `CF4.P7`, Adenda v2)

| Control | Escenario y Diseño (Bloque Ciego `20261120ULL`) | Criterio de Aceptación Sellado (Estricto) |
|---|---|---|
| **`CF4.P1`** | No regresión Nivel 1 con modo base (`usar_protocolo_huecos = false`) en `M1`, `M2`, `M3` (limpios y con huecos). | Volcado hexadecimal `%a` real generado por `ObservadorLilyth` **idéntico bit a bit** a `volcado_orig.txt` (`0` diferencias). |
| **`CF4.P2a`** | `NIVEL M3` + `PurgaDiurna`, micro-cortes limpios **$r = 2$** (tope operativo `NIVEL`) fuera de choques reales; `B4_RUBIN16` con escalado diurno causal (`EspejoPurgaNulo`) + Rubin $M = 16$. Tabla `B0`, `B1`, `B1b`, `B2`, `B3`, `B4_RUBIN16`. | $\|\bar{S}_M - S_0\| / \text{ee}_{\text{Rubin}} \le 0.25$; $\text{Alm\_Tmp\_Esp}(\text{B4}) \le \text{Alm\_Tmp\_Esp}(\text{B2})$; $\max_m \text{err\_anc}^{(m)} \le 10^{-10}$. |
| **`CF4.P2b`** | **Subcontrol H0 (Homogeneidad diurna, $r = 2$, $\Delta = 0$, $K = 192 = 8$ días tranquilos por cada una de las 24 fases):** evalúa `Z_hueco` con `escalas_fase` causal (y reporta A/B sin escalado).<br>**Subcontrol H1 (Matriz completa $4 \times 5$ de choque oculto en $K = 96 = 4$ posiciones tranquilas por fase):** $r \in \{1, 2, 3, 4\} \times A \in \{4, 6, 8, 10, 12\}\,\sigma_{\text{local}}$ ($\Delta_{\text{fis}} = A \cdot \bar{\kappa}_{\text{ef}}(t_h, r)\,\sigma_{e,\text{purg}}$). | **H0:** (1) $\text{SD}(Z_{\text{hueco}}) \in [0.80, 1.20]$; (2) razón pico/valle $\text{SD}(Z_{\text{pico}}) / \text{SD}(Z_{\text{valle}}) \in [0.70, 1.40]$ (caza `MN1.11`); (3) límite inferior de Wilson al 99% de la tasa de falsos bloqueos $\le \alpha_{\text{hueco}} = 0.01$.<br>**H1:** (4) En toda celda con $A \ge \Delta_{\min}(r; 0.01, 0.80)$, el límite superior de Wilson al 95% de $\hat{p}$ es $\ge 0.80$ **y** $\hat{p} \ge p_{\text{teor}} - 3\sqrt{p_{\text{teor}}(1 - p_{\text{teor}})/K}$; (5) en régimen autorizado ($r \le 2$) a $+12\,\sigma_{\text{local}}$: $\hat{p}(r=1, +12\sigma) \ge 95\%$ y $\hat{p}(r=2, +12\sigma) \ge 85\%$ con `FLAG_CHOQUE_EN_HUECO` activo y `reconstruido == false`. |
| **`CF4.P2c`** | Discriminación efectiva de `B4` (con escalado diurno causal) frente a `B2_LINEAL` y `B3_ITO_SIN_FILTRO` en `NIVEL M3`:<br>(i) $K = 48$ micro-cortes limpios $r = 2$ (2 por fase en las 24 fases): razón de varianza intra-hueco $\text{RazónVar}$ + telemetría `sd.temp[0].s2`.<br>(ii) $K = 12$ micro-cortes $r = 2$ espaciados a través del ciclo diurno con choques ocultos $+10\,\sigma_{\text{local}}$ pasando por `ObservadorSaltos` + `ObservadorLilyth`. | 1) En micro-cortes limpios ($r = 2$): $\text{RazónVar}(\text{B2\_LINEAL}) < 0.10$ y $\text{RazónVar}(\text{B4\_ITO\_RNG}) \in [0.60, 1.40]$.<br>2) En los $12$ micro-cortes con choque oculto $+10\,\sigma_{\text{local}}$: `B3_ITO_SIN_FILTRO` bloquea $0/12$, mientras `B4` bloquea $\ge 8/12$ con `FLAG_CHOQUE_EN_HUECO`. |
| **`CF4.P3`** | `M2` (`ARFIMA`, $16$ semillas ciegas, `forzar_via_b = true`, `INTERCAMBIABLE`, `Z1b` en `t_cambio - 2`, $r = 3$, Rubin $M = 16$ por serie). | $\text{RMSE}_{\text{rel}}(e_{\text{iny}}) \le 0.15$ y $\le \text{RMSE}(\text{B0b\_CORTE\_SUFIJO})$; $\|\overline{\Delta S}\|/\overline{\text{ee}}_{\text{Rubin}} \le 0.25$; $10^{-6} < \overline{\gamma}_M \le 0.25$; $\text{err\_anc} \le 10^{-10}$. |
| **`CF4.P4`** | Máquina de pausa dura en corte largo $r = 12 > r_{\max} = 3$. | Exactamente $9$ pasos en `EN_PAUSA_DURA` (`FLAG_PAUSA_UMBRAL`) y `cortar_segmento_al_reenganchar == true`. |
| **`CF4.P5`** | Monte Carlo ciego $R = 200$ (`q = 2%`, `r = 2`) en `M1_H0`, `M2_H0`, `M3_H0` (con escalado diurno causal para `B4`) y `M1_H1`. Rubin $M = 5$ + `ark::holm()`. | 1) Wilson 95% lo $\le 0.05$ en `M1_H0`, `M2_H0`, `M3_H0`.<br>2) Retención de potencia en `M1_H1` $\ge 0.85$.<br>3) $\|\bar{S}_M - \bar{S}_0\| / \overline{\text{ee}}_{\text{Rubin}} \le 0.25$ en `M1_H0`, `M2_H0`, `M3_H0`.<br>4) $\Delta\text{Censura}(\text{M1\_H0}) \le p_{\text{teor}}(\bar{l}) + 1.96\sqrt{\frac{p_{\text{teor}}(1 - p_{\text{teor}})}{R \cdot \bar{K}}}$. |
| **`CF4.P6`** | Blindaje obligatorio de `FLAG_PUENTE_ITO` en `Z2 && Z3` dentro de `ObservadorLilyth`. | `veto_z2 && veto_z3` exacto (`perfil_corrupto_z2_cruda && sin_nulo_b && !evaluado`; `eta[0] = NaN` y `tau_rest_consolidado` intacto $< 10^{-12}$). |
| **`CF4.P7`** | Deriva no nula en `FLUJO` ($\mu = +10.0\,\sigma_e$, AR(1) $\phi = 0.65$) con $12$ micro-cortes limpios $r = 2$ alimentando `ObservadorSaltos` vía `valor_z0_imputado(..., mod.mu)`. | `nrec_d >= 10`, **`0` alarmas negativas espurias** en las ventanas de hueco en `ObservadorSaltos`, y error contractual $\max |x_{\Delta}(t) - \hat{\mu}_{\text{calma}}| \le 10^{-12}$ dentro de los huecos reconstruidos. |

---

## 4. Mutantes Obligatorios Sellados (`mutaciones_fase4_puente_v2.py`, $R = 200$ fijo)

Compilación limpia con `/W4 /permissive-` (0 warnings), cada patrón con `count == 1` exacto, testigo limpio con `0` fallas y **5 de 5** mutantes cazados por su propia prueba específica:

1. **`MF4.0` (Testigo limpio):** sin cambios $\rightarrow$ `TESTIGO PASA (0 fallas)`.
2. **`MN1.7` (`observador_lilyth.hpp`):** anular `veto_z3_ito` y `veto_ito_z2` $\rightarrow$ cazado por **`CF4.P6`**.
3. **`MN1.8` (`puente_historico_ito.hpp`):** anular `if (std::abs(out.z_hueco) > out.z_max_efectivo)` $\rightarrow$ cazado por **`CF4.P2b`** y **`CF4.P2c`**.
4. **`MN1.9` (`puente_historico_ito.hpp`):** colapsar `out.b_entre = 0.0` en `combinar_rubin()` $\rightarrow$ cazado por **`CF4.P3`**.
5. **`MN1.10` (`puente_historico_ito.hpp`):** sustituir `reposo_flujo` por `0.0` literal en `valor_z0_imputado()` $\rightarrow$ cazado por **`CF4.P7`** (con disparo real de `ObservadorSaltos`).
6. **`MN1.11` (`puente_historico_ito.hpp`):** anular `usar_escalado_diurno` ($\kappa \equiv 1.0$) en `evaluar_y_reconstruir_hueco_ito()` $\rightarrow$ cazado por **`CF4.P2b`** (homogeneidad pico/valle bajo $H_0$).
