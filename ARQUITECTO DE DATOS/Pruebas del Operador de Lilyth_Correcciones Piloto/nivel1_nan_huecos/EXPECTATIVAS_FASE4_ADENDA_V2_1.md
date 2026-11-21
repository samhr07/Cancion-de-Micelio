# Expectativas Selladas — Fase 4, Adenda v2.1: Dimensionamiento Muestral en `CF4.P2c`/`CF4.P7`, Cuarto Bloque Ciego (`20261121ULL`) y Caza 5/5 de Mutantes

**Fecha de sello:** 2026-11-19  
**Predecesores inmutables:**  
- `EXPECTATIVAS_FASE4.md` (SHA-256 `49fe152d00248465489ca2c561647d61033197f2543e819f179bec25e4806ef1`) y sus bitácoras `test_fase4_ciego_rubin.log` / `mutaciones_fase4_puente.log`.  
- `EXPECTATIVAS_FASE4_ADENDA_V2.md` (SHA-256 `32ad43d217d9a792df5c80a23c1a968c1cdb5f9ed80f13814ce4479c5d332e1c`) y su bitácora `test_fase4_ciego_rubin_v2.log`.  
Todos los documentos anteriores y las cabeceras del operador (`puente_historico_ito.hpp`, `observador_lilyth.hpp`) permanecen **100% intactos** sin modificación alguna.  

**Regla de inmutabilidad:** Este documento se sella con SHA-256 en `SELLO_EXPECTATIVAS_FASE4_V2_1.txt` **antes** de compilar y ejecutar el arnés ciego `test_fase4_ciego_rubin_v2_1.cpp` y `mutaciones_fase4_puente_v2_1.py`. Ningún umbral, fórmula, tamaño muestral ($R, M, K$), ni bloque de semillas puede modificarse después del sello. Todo criterio que falle sobre el cuarto bloque ciego se reporta textualmente como `[FALLA]` sin ajustes post-hoc.

---

## 1. Motivo Exacto de la Adenda v2.1 (Cero Cambios al Operador)

En la corrida ciega de la Adenda v2 (`20261120ULL`), el algoritmo de `puente_historico_ito.hpp` pasó `CF4.P1`, `CF4.P2a`, `CF4.P2b` (con homogeneidad diurna pico/valle = `1.093`), `CF4.P3`, `CF4.P4`, `CF4.P5` y `CF4.P6`, y dejó en evidencia **3 defectos puros de dimensionamiento muestral en el arnés de pruebas**, sin ningún defecto en `puente_historico_ito.hpp` ni en `observador_lilyth.hpp`:

1. **`CF4.P2c` (i) — Denominador muestral con `r = 2`:** evaluar la razón de varianza intra-hueco sobre apenas $K = 48$ micro-huecos cortos ($r = 2$, $2$ grados de libertad por hueco) introduce una fluctuación muestral excesiva en el denominador físico $\text{Var}_{\text{física}}$. Usando las **$K = 192$** posiciones tranquilas estratificadas ($8$ por fase diurna), con guarda `n_rec_p2c >= 180`, se estabiliza el denominador físico preservando exactamente los mismos umbrales ($\text{RazónVar}(\text{B2}) < 0.10$, $\text{RazónVar}(\text{B4}) \in [0.60, 1.40]$).
2. **`CF4.P2c` (ii) — Recorrido incompleto para insertar `12` choques ocultos:** al escanear solo `muestras_h0` (que terminaba alrededor de $t \approx 3200$) con separación $\ge 320$ barras, solo cabían $9$ choques ocultos (`9/9` bloqueados por `B4` vs `0/9` por `B3`), fallando la aserción `(n_iny_ocultos == 12)`. Recorriendo todo el rango `[200, nop - 40)` con separación $\ge 300$ barras sobre las $12$ fases pares distintas $\{0, 2, \dots, 22\}$, se colocan efectivamente los $12$ choques ocultos $+10\,\sigma_{\text{local}}$ (`n_iny_ocultos == 12`, `chq_o3 == 0`, `chq_o4 >= 8`).
3. **`CF4.P7` — Guarda muestral sobre `12` huecos bajo llave $\alpha = 0.01$:** en una serie AR(1) limpia de $12$ huecos, la probabilidad de que $3$ o más huecos muestren $|Z_{\text{hueco}}| > 2.5758$ por cola normal o estimación es no despreciable (en `20261120ULL`, `B4` tuvo `alarmas_neg_en_hueco = 0` y error contractual `0.00e+00`, pero `huecos_rec = 9 < 10`). Aumentando a **$n_{\text{op}} = 5600$** con **$24$** micro-cortes limpios `r = 2` (`th = 400 + 200 * k`, $k = 0 \dots 23$) y guarda `nrec_d >= 18`, se preservan las exigencias estrictas (`alarmas_neg_hueco == 0` y `|xd - mu| <= 1e-12`) mientras el mutante `MN1.10` dispara alarmas negativas masivas en `ObservadorSaltos`.

---

## 2. Cuarto Bloque Ciego de Semillas (Estrictamente Disjunto de Calibración, v1 y v2)

Quedan prohibidas todas las semillas de calibración (`101, 202, 303, 20261117ULL`), de Nivel 1 (`20261118ULL`), de Fase 4 v1 (`20261119ULL`) y de Fase 4 v2 (`20261120ULL`). Toda la Adenda v2.1 corre sobre el cuarto bloque ciego **`SEMILLA_CIEGA_V2_1 = 20261121ULL`**:

- `CF4.P2a`, `CF4.P2b`, `CF4.P2c` (`NIVEL M3` + `PurgaDiurna`, $\phi = 0.65$, $n_{\text{op}} = 6000$): `20261121ULL`.
- `CF4.P3` (`M2` ARFIMA, $d = 0.25, \phi_1 = 0.40$, $n_{\text{op}} = 6000$, `forzar_via_b = true`): $16$ realizaciones `20261121ULL + 31337ULL * s` ($s = 0 \dots 15$).
- `CF4.P5` (`M1_H0`, `M2_H0`, `M3_H0`, `M1_H1`): $R = 200$ réplicas con `ParametrosReplicaH0H1::semilla_base = 20261121ULL` (`rep = 0..199`).
- `CF4.P6` (blindaje `Z2 && Z3`): `20261121ULL + 7777ULL`.
- `CF4.P7` (`FLUJO` con deriva $\mu = +10.0\,\sigma_e$, AR(1) $\phi = 0.65$, $n_{\text{op}} = 5600$, $24$ huecos): `20261121ULL ^ 0xF00DCAFEULL`.

---

## 3. Controles Falsables Sellados (`CF4.P1` – `CF4.P7`, Adenda v2.1)

| Control | Escenario y Diseño (Bloque Ciego `20261121ULL`) | Criterio de Aceptación Sellado (Estricto) |
|---|---|---|
| **`CF4.P1`** | No regresión Nivel 1 con modo base (`usar_protocolo_huecos = false`) en `M1`, `M2`, `M3` (limpios y con huecos). | Volcado hexadecimal `%a` real generado por `ObservadorLilyth` **idéntico bit a bit** a `volcado_orig.txt` (`0` diferencias). |
| **`CF4.P2a`** | `NIVEL M3` + `PurgaDiurna`, micro-cortes limpios **$r = 2$** (tope operativo `NIVEL`) fuera de choques reales; `B4_RUBIN16` con escalado diurno causal (`EspejoPurgaNulo`) + Rubin $M = 16$. Tabla `B0`, `B1`, `B1b`, `B2`, `B3`, `B4_RUBIN16`. | $\|\bar{S}_M - S_0\| / \text{ee}_{\text{Rubin}} \le 0.25$; $\text{Alm\_Tmp\_Esp}(\text{B4}) \le \text{Alm\_Tmp\_Esp}(\text{B2})$; $\max_m \text{err\_anc}^{(m)} \le 10^{-10}$. |
| **`CF4.P2b`** | **Subcontrol H0 (Homogeneidad diurna, $r = 2$, $\Delta = 0$, $K = 192 = 8$ días tranquilos por cada una de las 24 fases):** evalúa `Z_hueco` con `escalas_fase` causal (y reporta A/B sin escalado).<br>**Subcontrol H1 (Matriz completa $4 \times 5$ de choque oculto en $K = 96 = 4$ posiciones tranquilas por fase):** $r \in \{1, 2, 3, 4\} \times A \in \{4, 6, 8, 10, 12\}\,\sigma_{\text{local}}$ ($\Delta_{\text{fis}} = A \cdot \bar{\kappa}_{\text{ef}}(t_h, r)\,\sigma_{e,\text{purg}}$). | **H0:** (1) $\text{SD}(Z_{\text{hueco}}) \in [0.80, 1.20]$; (2) razón pico/valle $\text{SD}(Z_{\text{pico}}) / \text{SD}(Z_{\text{valle}}) \in [0.70, 1.40]$ (caza `MN1.11`); (3) límite inferior de Wilson al 99% de la tasa de falsos bloqueos $\le \alpha_{\text{hueco}} = 0.01$.<br>**H1:** (4) En toda celda con $A \ge \Delta_{\min}(r; 0.01, 0.80)$, el límite superior de Wilson al 95% de $\hat{p}$ es $\ge 0.80$ **y** $\hat{p} \ge p_{\text{teor}} - 3\sqrt{p_{\text{teor}}(1 - p_{\text{teor}})/K}$; (5) en régimen autorizado ($r \le 2$) a $+12\,\sigma_{\text{local}}$: $\hat{p}(r=1, +12\sigma) \ge 95\%$ y $\hat{p}(r=2, +12\sigma) \ge 85\%$ con `FLAG_CHOQUE_EN_HUECO` activo y `reconstruido == false`. |
| **`CF4.P2c`** | Discriminación efectiva de `B4` (con escalado diurno causal) frente a `B2_LINEAL` y `B3_ITO_SIN_FILTRO` en `NIVEL M3`:<br>(i) **$K = 192$** micro-cortes limpios $r = 2$ ($8$ por fase en las $24$ fases): razón de varianza intra-hueco $\text{RazónVar}$ + telemetría `sd.temp[0].s2`.<br>(ii) **$K = 12$** micro-cortes $r = 2$ en `[200, nop - 40)` separados $\ge 300$ barras sobre las $12$ fases pares distintas con choques ocultos $+10\,\sigma_{\text{local}}$ pasando por `ObservadorSaltos` + `ObservadorLilyth`. | 1) En micro-cortes limpios ($r = 2$, `n_rec_p2c >= 180`): $\text{RazónVar}(\text{B2\_LINEAL}) < 0.10$ y $\text{RazónVar}(\text{B4\_ITO\_RNG}) \in [0.60, 1.40]$.<br>2) En los $12$ micro-cortes con choque oculto $+10\,\sigma_{\text{local}}$ (`n_iny_ocultos == 12`): `B3_ITO_SIN_FILTRO` bloquea $0/12$, mientras `B4` bloquea $\ge 8/12$ con `FLAG_CHOQUE_EN_HUECO`. |
| **`CF4.P3`** | `M2` (`ARFIMA`, $16$ semillas ciegas, `forzar_via_b = true`, `INTERCAMBIABLE`, `Z1b` en `t_cambio - 2`, $r = 3$, Rubin $M = 16$ por serie). | $\text{RMSE}_{\text{rel}}(e_{\text{iny}}) \le 0.15$ y $\le \text{RMSE}(\text{B0b\_CORTE\_SUFIJO})$; $\|\overline{\Delta S}\|/\overline{\text{ee}}_{\text{Rubin}} \le 0.25$; $10^{-6} < \overline{\gamma}_M \le 0.25$; $\text{err\_anc} \le 10^{-10}$. |
| **`CF4.P4`** | Máquina de pausa dura en corte largo $r = 12 > r_{\max} = 3$. | Exactamente $9$ pasos en `EN_PAUSA_DURA` (`FLAG_PAUSA_UMBRAL`) y `cortar_segmento_al_reenganchar == true`. |
| **`CF4.P5`** | Monte Carlo ciego $R = 200$ (`q = 2%`, `r = 2`) en `M1_H0`, `M2_H0`, `M3_H0` (con escalado diurno causal para `B4`) y `M1_H1`. Rubin $M = 5$ + `ark::holm()`. | 1) Wilson 95% lo $\le 0.05$ en `M1_H0`, `M2_H0`, `M3_H0`.<br>2) Retención de potencia en `M1_H1` $\ge 0.85$.<br>3) $\|\bar{S}_M - \bar{S}_0\| / \overline{\text{ee}}_{\text{Rubin}} \le 0.25$ en `M1_H0`, `M2_H0`, `M3_H0`.<br>4) $\Delta\text{Censura}(\text{M1\_H0}) \le p_{\text{teor}}(\bar{l}) + 1.96\sqrt{\frac{p_{\text{teor}}(1 - p_{\text{teor}})}{R \cdot \bar{K}}}$. |
| **`CF4.P6`** | Blindaje obligatorio de `FLAG_PUENTE_ITO` en `Z2 && Z3` dentro de `ObservadorLilyth`. | `veto_z2 && veto_z3` exacto (`perfil_corrupto_z2_cruda && sin_nulo_b && !evaluado`; `eta[0] = NaN` y `tau_rest_consolidado` intacto $< 10^{-12}$). |
| **`CF4.P7`** | Deriva no nula en `FLUJO` ($\mu = +10.0\,\sigma_e$, AR(1) $\phi = 0.65$, $n_{\text{op}} = 5600$) con **$24$** micro-cortes limpios $r = 2$ (`th = 400 + 200*k`) alimentando `ObservadorSaltos` vía `valor_z0_imputado(..., mod.mu)`. | `nrec_d >= 18`, **`0` alarmas negativas espurias** en las ventanas de hueco en `ObservadorSaltos`, y error contractual $\max |x_{\Delta}(t) - \hat{\mu}_{\text{calma}}| \le 10^{-12}$ dentro de los huecos reconstruidos. |

---

## 4. Mutantes Obligatorios Sellados (`mutaciones_fase4_puente_v2_1.py`, $R = 200$ fijo)

Compilación limpia con `/W4 /permissive-` (0 warnings), cada patrón con `count == 1` exacto, testigo limpio con `0` fallas y **5 de 5** mutantes cazados por su propia prueba específica:

1. **`MF4.0` (Testigo limpio):** sin cambios $\rightarrow$ `TESTIGO PASA (0 fallas)`.
2. **`MN1.7` (`observador_lilyth.hpp`):** anular `veto_z3_ito` y `veto_ito_z2` $\rightarrow$ cazado por **`CF4.P6`**.
3. **`MN1.8` (`puente_historico_ito.hpp`):** anular `if (std::abs(out.z_hueco) > out.z_max_efectivo)` $\rightarrow$ cazado por **`CF4.P2b`** y **`CF4.P2c`**.
4. **`MN1.9` (`puente_historico_ito.hpp`):** colapsar `out.b_entre = 0.0` en `combinar_rubin()` $\rightarrow$ cazado por **`CF4.P3`**.
5. **`MN1.10` (`puente_historico_ito.hpp`):** sustituir `reposo_flujo` por `0.0` literal en `valor_z0_imputado()` $\rightarrow$ cazado por **`CF4.P7`** (con disparo real de `ObservadorSaltos`).
6. **`MN1.11` (`puente_historico_ito.hpp`):** anular `usar_escalado_diurno` ($\kappa \equiv 1.0$) en `evaluar_y_reconstruir_hueco_ito()` $\rightarrow$ cazado por **`CF4.P2b`** (homogeneidad pico/valle bajo $H_0$).
