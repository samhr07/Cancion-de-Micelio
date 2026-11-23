# Expectativas Selladas — Nivel 3: Simetría de Entorno, Ventana Derivada y Parada Temprana V4 en Vía (b)

**Fecha de sello:** 2026-11-22  
**Base de referencia inmutable:** `lilyth-auditoria-cota @ 6621797f7deaa5ecf6f795d24dcc038ef1540aab` (`ARQUITECTO DE DATOS/Pruebas del Operador de Lilyth_Correcciones Piloto/`)  
**Predecesores inmutables (SHA-256):**
- Nivel 1 (`EXPECTATIVAS_COTA_NIVEL1.md`): `54aa457e9dc69ddec1a982f4904e2d71e354f4ae9ea2dada4aba77c7f52eafa3`
- Fase 4 `v1` (`EXPECTATIVAS_FASE4.md`): `49fe152d00248465489ca2c561647d61033197f2543e819f179bec25e4806ef1`
- Fase 4 `v2` (`EXPECTATIVAS_FASE4_ADENDA_V2.md`): `32ad43d217d9a792df5c80a23c1a968c1cdb5f9ed80f13814ce4479c5d332e1c`
- Fase 4 `v2.1` (`EXPECTATIVAS_FASE4_ADENDA_V2_1.md`): `f9b0aed3ca2f3117789e672d7770647bcba3f1c53c474a03fc3653d9245f682d`
- Fase 4 `v2.2` (`EXPECTATIVAS_FASE4_ADENDA_V2_2.md`): `7c63f1dcbdfdd14e5ec11f65a1df08ef9cabf55b2b339b86ce39092bd1b21242`
- Nivel 2 (`EXPECTATIVAS_NIVEL2.md`, commit `6621797f7deaa5ecf6f795d24dcc038ef1540aab`): `2387e846df90c7b87bb90d70cbd7d32b8d54ca85ea269386464a260fe6187d41`
- Sonda empírica de diseño (`nivel2_3_auditoria_argon/sonda_glr_reinicio_via_b.log @ dda2db0c`) y contrato cerrado (`contrato_nivel2_y_nivel3_via_b_lilyth.md`).

**Regla de inmutabilidad y Regla 6:**
- Este documento se sella con SHA-256 en `SELLO_EXPECTATIVAS_NIVEL3.txt` **antes** de editar `nulo_lilyth.hpp` o `observador_lilyth.hpp`. Ningún umbral, fórmula, tamaño muestral ni semilla puede alterarse después del sello; cualquier incumplimiento se reporta textualmente como `[FALLA]`.
- Bandera explícita opt-in **`ConfiguracionLilyth::usar_enmienda_nivel3 = false`** por omisión en piloto para preservar no-regresión bit a bit (`CN1.1`, `CF4.P1`, Nivel 1, Fase 4 y Nivel 2).

---

## 1. Contrato Matemático y Semillas Ciegas del Nivel 3

Queda prohibido reusar las semillas de calibración (`101, 202, 303, 20261117ULL`), Fase 4 (`20261118ULL`–`20261122ULL`), Nivel 2 (`20261123ULL`), sonda (`2026112100ULL`/`2026120000ULL`) y `etapa_lh` (`11000`–`41000`). Toda la batería corre bajo el bloque ciego exclusivo **`SEMILLA_CIEGA_NIVEL3 = 20261124ULL`**:

- **`NV3.1` (Banco `R1`, `40` semillas $s = 0 \dots 39$):**
  - `sem_cal(s) = splitmix64(20261124ULL ^ (0x9E3779B97F4A7C15ULL * (s + 1)))`
  - `sem_op(s) = splitmix64(20261124ULL ^ (0xBF58476D1CE4E5B9ULL * (s + 101)))`
  - `sem_nulo(s) = splitmix64(20261124ULL ^ (0x94D049BB133111EBULL * (s + 211)))`
- **`NV3.2` (Banco `R2`, `24` semillas: `R2a` $s = 0 \dots 11$ y `R2b` $s = 12 \dots 23$):**
  - `sem_cal(s) = splitmix64(20261124ULL ^ (0x517CC1B727220A95ULL * (s + 1)))`
  - `sem_op(s) = splitmix64(20261124ULL ^ (0xD6E8FEB86659FD93ULL * (s + 307)))`
  - `sem_nulo(s) = splitmix64(20261124ULL ^ (0xA0761D6478BD642FULL * (s + 409)))`
- **`NV3.3a` (`S2_N3` y `n_sim`):** `20261124ULL + 1009ULL`
- **`NV3.3b` (`ARFIMA` `mu_serie` + `V4`):** `20261124ULL + 3001ULL`
- **`NV3.3c` (`K10` con $l > 256$ y `K1`):** `20261124ULL + 7001ULL`
- **`NV3.3d` (No-regresión Regla 6):** semillas canónicas de `volcado_orig.txt` (`20261117ULL`), Fase 4 (`20261122ULL`) y Nivel 2 (`20261123ULL`).

### Ecuaciones exactas cerradas bajo `usar_enmienda_nivel3 = true`

1. **Ventana de exclusión derivada (sin `T_MAX` fijo):**
   ```cpp
   ventanas[i] = construir_ventana_inyeccion(tiempos_iny[i], l, w_ancla, guarda, t_max_ef);
   ```
   El solape aplica únicamente entre pares de inyecciones deterministas $i \neq j$.
2. **Longitud `n_sim` y soporte de sorteo sin literales `4/8/16`:**
   ```cpp
   const std::size_t w_pre = w_ancla + guarda + (integracion == Integracion::NIVEL ? 1u : 0u);
   const std::size_t w_sorteo = std::max<std::size_t>(1, w_ancla);
   const std::size_t w_post = l + t_max_ef;
   const std::size_t n_sim = w_pre + w_sorteo + w_post;
   ```
   Sorteo uniforme $t_{\text{iny}} \sim \mathcal{U}[w_{\text{pre}},\, w_{\text{pre}} + w_{\text{sorteo}} - 1]$.
3. **Retirada del veto asimétrico `hay_previa` en inyección simple:**
   - Ni `hay_previa` ni `hay_posterior` marcan `solapada[i] = true` bajo `usar_enmienda_nivel3 = true`.
   - `con_ajena_previa` y `con_reapertura_posterior` se conservan como telemetría pura en `ResultadoInyeccionesRealizacion` y se agregan en `ReferenciaNuloEvento` (`con_ajena_previa_b`).
4. **Parada temprana `V4_OPCION_B_PARADA_TEMP` y telemetría de escala:**
   - Al abrir `lil_local` una ranura con `signo == signo_real` en $t_{\text{ap}} \in [t_0, t_0 + l - 2]$, cesa inmediatamente la inyección restante `e_iny[t_ap - t0 + 1 .. l - 1]` (o `perfil_est` en `sin_memoria`), restaurando las innovaciones del nulo fresco `eps_nulo` desde $t_{\text{ap}} + 1$ y recalculando `u`, `x_est` (con `mu_serie` + `mu` en `ARFIMA`) y `x_sim` (`FLUJO`/`NIVEL`/`PURGADA`) para todo $t > t_{\text{ap}}$, de modo que $\tau \ge 1$ evoluciona $100\%$ bajo proceso nulo puro.
   - La admisión sigue exigiendo `signo == signo_real` en $[t_0, t_0 + l]$ y registra `paradas_tempranas_b` y `escala_coincide_b` (`ev.escala_origen == escala_real`) como telemetría diagnóstica para Nivel 4.
   - Se admite un gancho opcional genérico `GanchoHijas` sin costo en tiempo de ejecución (por omisión `SinHijasPaso{}` sin efecto) en `ejecutar_inyecciones_en_realizacion` / `ejecutar_via_b_evento` para inyectar réplicas endógenas en línea en pruebas de procesos autoexcitantes (`R2`) sin acoplar `replicas_lilyth.hpp` ni modificar `etapa_lh`.

---

## 2. Controles Falsables Sellados (`NV3.1` – `NV3.3`)

| Control | Escenario y Diseño (`SEMILLA_CIEGA_NIVEL3 = 20261124ULL`) | Criterio de Aceptación Sellado (Estricto) |
|---|---|---|
| **`NV3.1`** | **Banco `R1` (`AR(1)` $\phi = 0.80$, `40` semillas ciegas, `forzar_via_b = true`):** `ObservadorSaltos<24, 64, 120>`, `ObservadorLilyth<ObsDelta, 8, 64>`, `FLUJO`, `SUMA`, `arl0_pasos_nativos = 450.0`, `escalas_max = 3`, `alfa = ALFA_PROYECTO (0.05)`, `n_cal = 4096`, `n_op = 6144`, `b_max_nulo = 16`, `b_mc_resumen = 1999`, `contraste_resumen = INTERCAMBIABLE`, `forzar_via_b = true`, `usar_enmienda_nivel2 = true`, `usar_enmienda_nivel3 = true`. | - $\|S_{\text{medio}}\| \le 2.0\,\text{ee}(S_{\text{medio}})$ sobre las 40 semillas ($\|z_S\| \le 2.0$, con $\ge 38/40$ semillas evaluables).<br>- Tasa de rechazo bilateral a $\alpha = 0.05$ (`p_bilateral < 0.05` bajo `INTERCAMBIABLE`) $\in [0.0, 0.10]$ ($\le 4/40$).<br>- $p_{\det} = \sum n_{\text{nulos},b} / \sum \text{intentos}_b \ge 0.40$ y `paradas_tempranas_b_total > 0`.<br>- Supervivencia al corte `S3` (`p_det_inicial > 0`) $\ge 95\%$ de eventos reales medibles y cobertura completa `n_nulos_b >= m_min` bajo el tope `S3` (`ceil(m_min / p_det_inicial)`) $\ge 55\%$ *(nota: truncamiento binomial esperado ~25–30% con $p_{\det} \approx 0.5$; no se relaja `S3`)*.<br>- `solapes_b_total == 0` en inyecciones simples y `con_ajena_previa_b_total > 0`; registro informativo de `escala_coincide_b` para Nivel 4. |
| **`NV3.2`** | **Banco `R2` (`AR(1) + Hawkes` 2 signos bajo $H_0$ exacto, `24` semillas ciegas):** `ObservadorSaltos<24, 4, 120>`, `ObservadorLilyth<ObsH, 8, 128, 2048>`, `FLUJO`, `SUMA`, `arl0_pasos_nativos = 800.0`, `escalas_max = 1`, `alfa = 0.05`, $\phi = 0.88$, $\sigma_e = \sqrt{1 - \phi^2}$, `n_cal = 8192`, `n_test = 6000`, `b_max_nulo = 8`, `b_mc_resumen = 1999`, `INTERCAMBIABLE`, `forzar_via_b = true`, `usar_enmienda_nivel2 = true`, `usar_enmienda_nivel3 = true`. Madres `tasa_madre = 1/300` ($\pm 1$ equiprobable, $a_1 = 9.5\,\sigma_e$) + hijas en cada `o.abrio_ranura` con retardo $d_0 + \text{Geom}(1 - e^{-1/\tau_{\text{ret}}})$ (`d0 = 6`, `tau_ret = 4.0`, $a_1 = 9.5\,\sigma_e$) vía `GanchoHijas`:<br>- `R2a` ($s = 0 \dots 11$): `n_mismo = 0.30`, `n_contrario = 0.15` ($\rho = 0.45$).<br>- `R2b` ($s = 12 \dots 23$): `n_mismo = 0.45`, `n_contrario = 0.25` ($\rho = 0.70$). | - Todas las `24/24` semillas con contraste `INTERCAMBIABLE` válido.<br>- Semillas con $S < 0$ dentro del intervalo binomial $[7/24,\, 17/24]$ (corrigiendo el $20/24$ previo).<br>- Falsos rechazos (`p_bilateral < 0.05`) $\le 3/24$.<br>- Actividad de 2 signos verificada (`cortadas_b_total > 0` y `con_reapertura_posterior_total > 0`), sin descartes por `hay_previa` (`solapes_b_total == 0` y `con_ajena_previa_b_total > 0`). |
| **`NV3.3`** | **Integridad Estructural (`NV3.3a`–`NV3.3d`):**<br>a) **`S2_N3` y `n_sim` derivado:** `AR(1)` $\phi = 0.85$, `ARL0 = 1e6`, `T_MAX = 128` (`t_max_ef < T_MAX`): par en borde `w_ancla + guarda + l + t_max_ef` vs par en `w_ancla + guarda + l + t_max_ef + 16`; exactitud de `n_sim` en `FLUJO` y `NIVEL`.<br>b) **`ARFIMA` `mu_serie` + `V4`:** `ARFIMA(1, d=+0.35, 0)` con `paradas_tempranas > 0`, traslación $\mu_z = 50.0 \mapsto 1050.0$, y coherencia en `NIVEL` + `PURGADA`.<br>c) **`K10` ($l > 256$) + `K1`:** `ObservadorSaltos<24, 120, 120>`, `escalas_max = 3`, ambas enmiendas `true`, perfil real $l > 256$ evaluado por vía `(b)`, contador global `operator new` e identidad `K1`.<br>d) **No-regresión completa Regla 6:** `ConfiguracionLilyth{}.usar_enmienda_nivel3 == false` por omisión; verificación externa `CN1.1`, Nivel 1, Fases 1/2/3, Fase 4 `v2.2` (`CF4.P1`), Nivel 2 y sellos SHA-256 previos. | a) Borde: `solapes == 2` y `0` admitidas; intermedio (`< ... + T_MAX`): `solapes == 0` y `2` admitidas; `n_sim == w_pre + max(1, w_ancla) + l + t_max_ef` exacto.<br>b) Invariancia $\max \|x_{\text{sim}}^{(1050)} - (x_{\text{sim}}^{(50)} + 1000)\| \le 10^{-10}$ y finitud/consistencia en `NIVEL` + `PURGADA`.<br>c) `reservas_en_paso == 0`, `perfil.size() == l > 256`, `e_iny.size() == l`, `n_nulos_b > 0`, y `K1` bit a bit (`0` diferencias).<br>d) `usar_enmienda_nivel3 == false` por omisión, `CN1.1` (`114/114` líneas `%a` idénticas a `volcado_orig.txt`), baterías previas con `0` fallas y sellos intactos. |

---

## 3. Mutantes Obligatorios Sellados (`mutaciones_nivel3.py`)

Compilación limpia con `/W4 /permissive-` (0 warnings), cada sustitución con `count == 1` exacto, testigo limpio `MN3.0` con `0` fallas y **3 de 3** mutantes cazados por su propia prueba específica:

1. **`MN3.0` (Testigo limpio):** sin cambios $\rightarrow$ `TESTIGO PASA (0 fallas)`.
2. **`M_N3_1` (`nulo_lilyth.hpp`):** desactivar parada temprana `V4_OPCION_B_PARADA_TEMP` (seguir inyectando `e_iny[t - t0]` después de que `lil_local` abra en $t_{\text{ap}} < t_0 + l - 1$) $\rightarrow$ cazado por **`NV3.1`** (y reportado en `NV3.2`).
3. **`M_N3_2a` (`nulo_lilyth.hpp`):** reactivar el veto `if (hay_previa) solapada[i] = true;` bajo `usar_enmienda_nivel3 = true` $\rightarrow$ cazado por **`NV3.1`** y **`NV3.2`**.
4. **`M_N3_2b` (`nulo_lilyth.hpp`):** usar `T_MAX` fijo en lugar de `t_max_ef` en `construir_ventana_inyeccion(...)` bajo `usar_enmienda_nivel3 = true` $\rightarrow$ cazado por **`NV3.3`** (`NV3.3a / S2_N3`).
