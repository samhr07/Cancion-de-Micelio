# Expectativas Selladas — Nivel 2: Horizonte Causal Dinámico de Hosking, Memoria Adaptable Sin Heap y Protocolo Híbrido `Z1b`

**Fecha de sello:** 2026-11-21  
**Base de referencia inmutable:** `lilyth-auditoria-cota @ dda2db0c` (`ARQUITECTO DE DATOS/Pruebas del Operador de Lilyth_Correcciones Piloto/`)  
**Predecesores inmutables:**
- Nivel 1 (`EXPECTATIVAS_COTA_NIVEL1.md`, SHA-256 `54aa457e9dc69ddec1a982f4904e2d71e354f4ae9ea2dada4aba77c7f52eafa3`).
- Fase 4 `v1`–`v2.2` (`EXPECTATIVAS_FASE4.md` `49fe152d...`, `V2` `32ad43d2...`, `V2_1` `f9b0aed3...`, `V2_2` `7c63f1dc...`).
- Contrato consensuado (`contrato_nivel2_y_nivel3_via_b_lilyth.md`) con las **4 precisiones de ingeniería** aprobadas por el Comandante:
  1. **Memoria adaptable al rango calibrado:** `cap_perfil_adapt_ = min(HIST_MAX, Rasgos::l_glr * max_m_calibrado)` en arena pre-reservada durante `configurar()`/`calibrar()`, con cero reservas en `paso()`, sin inflar `sizeof(RanuraLilyth)` y sin copiar la arena histórica dentro de `lil_local`.
  2. **Inversión algebraica exacta en `simular()`:** centrar `mu_serie` antes de `(1-L)^d` en `ajustar()` e `innovaciones()`, e invertir exactamente en `simular()` mediante `x_out[t] = m.mu_serie + sum_k pi_k(-d) * (m.mu + u_c)` (con `tamano_trabajo_simular(m_lin, n)` = `2 * tot` también en `construir_estacionaria_nula()`).
  3. **Banderas explícitas opt-in (Regla 6):** `usar_enmienda_nivel2 = false` por omisión en piloto para preservar no-regresión bit a bit (`CN1.1`/`CF4.P1`); veto `Z2` exclusivo a vía `(b)` sin cerrar la ranura si vía `(a)` la cubre; y prohibición estricta en `Z1a` (`t_cambio - 1`) de usar `t_cambio` como borde derecho del puente de Itô.
  4. **Ejecución aislada del DAG:** completar e independizar el **Bloque Nivel 2 (`T2.1`–`T2.9`)** antes de abrir el **Bloque Nivel 3 (`T3.1`–`T3.8`)**, sin tocar `master` ni `ark_cxx/base_auditada/`.

**Regla de inmutabilidad:** Este documento se sella con SHA-256 en `SELLO_EXPECTATIVAS_NIVEL2.txt` **antes** de modificar `modelo_lineal_nulo.hpp`, `observador_lilyth.hpp` o `nulo_lilyth.hpp`. Ningún umbral, fórmula, tamaño muestral ni semilla puede alterarse después del sello. Cualquier incumplimiento se reporta textualmente como `[FALLA]`.

---

## 1. Contrato Matemático y Semillas Ciegas del Nivel 2

Quedan prohibidas las semillas de calibración (`101, 202, 303, 20261117ULL`) y de bloques anteriores (`20261118ULL`–`20261122ULL`). Toda la batería corre bajo el bloque **`SEMILLA_CIEGA_NIVEL2 = 20261123ULL`**:

- **`NV2.1` (`AR(p)`, $d = 0$, $p \in \{1, 4, 8\}$):** `20261123ULL + 101ULL * p`.
- **`NV2.2` (`ARFIMA`, $d \in \{-0.35, +0.25, +0.35\}$):** `20261123ULL + 1009ULL * idx_d`.
- **`NV2.3a` (supervivencia `Z1b` con `q = 2%`, ráfagas `r = 2`, $d = +0.35$):** `16` series con `20261123ULL + 7919ULL * s` ($s = 0 \dots 15$).
- **`NV2.3b` (cortes duros $r=4$, choque $+12\sigma$, corte `hist_seg_purgada_`, `Z1a`, `Z2`, $n_{\text{hist}} < p$):** `20261123ULL + 5003ULL`.
- **`NV2.4` (`K10` cero-heap con $l > 256$):** `20261123ULL + 9001ULL`.
- **`NV2.5` (no-regresión `CN1.1`/`CF4.P1` con `usar_enmienda_nivel2 = false`):** semillas canónicas de `volcado_orig.txt` (`20261117ULL`).

### Fórmulas exactas obligatorias (`eps_Q^2(alfa=0.05, n=4096) = z_{0.975}^2 / 4096 = 0.0009378602...`)
- **`AR(p)` ($d = 0$):** `m_inv_arfima = 0`, `m_sim_arfima = 0`, `n_hist_req = p`.
- **`ARFIMA(p, d, 0)` ($d \neq 0$):** mediante la autocovarianza cerrada de Hosking (1981),
  $$\gamma_z(0) = \frac{\Gamma(1 - 2d)}{\Gamma(1 - d)^2}, \qquad \gamma_z(h) = \gamma_z(h - 1)\frac{h - 1 + d}{h - d},$$
  $$\text{MSE}_\pi(d, M) = \gamma_z(0)\sum_{k=0}^M \pi_k(d)^2 + 2\sum_{h=1}^M \gamma_z(h)\sum_{j=0}^{M-h} \pi_j(d)\pi_{j+h}(d) - 1,$$
  se fija `m_inv_arfima = min{ M in [1, cap] : MSE_pi(d, M) <= eps_Q^2 }` y `n_hist_req = p + m_inv_arfima`. Para `simular()`, `m_sim_arfima = M_gauss(-d; eps_Q^2, max_seg)` con `burn_in = min(max_seg, B_AR + m_sim_arfima)` y `burn_in_topado = (M_gauss_sin_tope > max_seg)`.

---

## 2. Controles Falsables Sellados (`NV2.1` – `NV2.5`)

| Control | Escenario y Diseño (`SEMILLA_CIEGA_NIVEL2 = 20261123ULL`) | Criterio de Aceptación Sellado (Estricto) |
|---|---|---|
| **`NV2.1`** | **Exactitud Markov en `AR(p)` ($d = 0$):** modelos estables con $p \in \{1, 4, 8\}$ ($n = 4096$), perfil $l = 32$, `innovaciones()` con `n_hist = p` frente a `n_hist = 2048`. | `m.n_hist_req == p`, `m.m_inv_arfima == 0`, `m.m_sim_arfima == 0`, y $\max_j \|e_j^{(p)} - e_j^{(2048)}\| \le 10^{-15}$. |
| **`NV2.2`** | **Horizonte exacto de Hosking + Pre-centrado `mu_serie` + Ida y vuelta en `ARFIMA(p, d, 0)`:**<br>a) Verificación analítica de `horizonte_inv_arfima_hosking(d, eps_Q^2, 4096)` para $d \in \{-0.40, -0.35, -0.25, +0.25, +0.35, +0.40\}$ ($n = 4096, \alpha = 0.05$).<br>b) Invariancia de traslación $\mu_z = 50.0 \mapsto 1050.0$ e ida y vuelta `innovaciones -> simular` (`t_iny = 0, l = n = 512`) para $d \in \{-0.35, +0.25, +0.35\}$, $\phi_1 = 0.40$.<br>c) Precisión de cola fraccionaria sobre $K = 512$ realizaciones estacionarias independientes ($d \in \{-0.35, +0.25, +0.35\}$): comparar `M_pi_exact` frente al truncado `M = 8` (fórmula de ruido blanco para $d > 0$). | a) Coincidencia exacta con `{338, 211, 82, 85, 233, 418}`, `n_hist_req == p + m_inv_arfima`, `burn_in_topado == true` en $d \ge +0.25$ y `false` en $d \le -0.25$.<br>b) $\max_j \|e_j(1050) - e_j(50)\| \le 10^{-11}$ y $\max_t \|x_{\text{rec}}[t] - z[t]\| \le 10^{-12}$.<br>c) $\text{RMSE}(\hat{e}^{(M_\pi^{\text{exact}})} - \varepsilon) / \sigma_\varepsilon \le 1.15\,\varepsilon_Q$ ($\varepsilon_Q = 0.0306245$) y estrictamente menor que `RMSE(M = 8)`. |
| **`NV2.3`** | **Protocolo Híbrido `Z1b` + Cortes Duros + `Z1a`/`Z2` + `n_hist < p`:**<br>a) Supervivencia `Z1b` bajo `q = 2.0%` en ráfagas $r = 2$ sobre `ARFIMA` ($d = +0.35$, `n_hist_req = 234`), `64` ventanas causales en `16` series independientes: comparar `PUENTE_ITO_HISTORICO` (`r_max = 2, alfa_hueco = 0.01`) frente a `CORTE_SUFIJO`.<br>b) Cortes duros en `ObservadorLilyth` (`usar_enmienda_nivel2 = true`): (1) hueco largo $r = 4 > r_{\max}=2$ en `Z1b` a distancia `d_cut` de `t_cambio`; (2) choque oculto $+12\,\sigma$ en micro-hueco $r = 2$ en `Z1b`; (3) cambio de `hist_seg_purgada_` a distancia `d_seg` en rama `PURGADA`; (4) hueco en borde inmediato `Z1a` (`t_cambio - 1`); (5) `r = 1`, `FLAG_PUENTE_ITO` o corte de segmento dentro de `Z2`; (6) `n_hist_prev_disp < mod.p`. | a) `PUENTE_ITO_HISTORICO` logra `n_hist_prev_disp == n_hist_req` en $\ge 85\%$ de ventanas (con `n_huecos_puenteados_z1b > 0`) y $\text{RMSE}(e)/\sigma_\varepsilon \le 0.10$, frente a supervivencia $< 35\%$ de `CORTE_SUFIJO`.<br>b) (1) `break` exacto en `n_hist_prev_disp == d_cut` con `hist_prev_truncada_en_corte == true`; (2) `break` exacto en borde del choque por `\|Z_hueco\| > z_max`; (3) `break` exacto en `n_hist_prev_disp == d_seg` con `hist_prev_truncada_en_corte == true`; (4) `Z1a` jamás usa `t_cambio` como borde derecho del puente; (5) `Z2` marca `perfil_corrupto_z2 == true` y veta solo vía `(b)` (`sin_nulo_b == true`) manteniendo vía `(a)` intacta si procede; (6) `n_hist_prev_disp < p` marca `sin_nulo_b == true`. |
| **`NV2.4`** | **`K10` Cero-Heap en `paso()` y Memoria Adaptable con `l > 256`:** contador global de `operator new` activo durante todo `paso()` con `ObservadorSaltos<24, 120, 120>`, `escalas_max = 3` (`cap_perfil_adapt_ >= 480`), y eventos cerrados con perfil real de longitud $l > 256$. | `reservas_en_paso == 0`, `perfil_evento_n2().size() == l` (sin truncar en `256`), `e_iny_evento_n2().size() == l`, y ausencia de arreglos estáticos de tamaño `HIST_MAX` dentro de `RanuraLilyth`. |
| **`NV2.5`** | **No Regresión Bit a Bit (`CN1.1` / `CF4.P1`):** modo por omisión (`usar_enmienda_nivel2 = false`) sobre `M1`, `M2` y `M3` (limpios y con huecos). | Volcado `%a` idéntico byte a byte a `nivel1_nan_huecos/volcado_orig.txt` (`0` diferencias). |

---

## 3. Mutantes Obligatorios Sellados (`mutaciones_nivel2.py`)

Compilación limpia con `/W4 /permissive-` (0 warnings), cada sustitución con `count == 1` exacto, testigo limpio `MN2.0` con `0` fallas y **3 de 3** mutantes cazados por su propia prueba específica:

1. **`MN2.0` (Testigo limpio):** sin cambios $\rightarrow$ `TESTIGO PASA (0 fallas)`.
2. **`M_N2_1` (`modelo_lineal_nulo.hpp`):** omitir el restado de `m.mu_serie` antes de `(1-L)^d` en `innovaciones()` $\rightarrow$ cazado por **`NV2.2`**.
3. **`M_N2_2` (`modelo_lineal_nulo.hpp`):** sustituir `horizonte_inv_arfima_hosking` por la fórmula de ruido blanco (`horizonte_sim_gauss(+d)`) en `ajustar()` $\rightarrow$ cazado por **`NV2.2`**.
4. **`M_N2_3` (`observador_lilyth.hpp`):** comprimir rezagos o ignorar el corte por cambio de `hist_seg_purgada_` en el escaneo `Z1b` $\rightarrow$ cazado por **`NV2.3`**.
