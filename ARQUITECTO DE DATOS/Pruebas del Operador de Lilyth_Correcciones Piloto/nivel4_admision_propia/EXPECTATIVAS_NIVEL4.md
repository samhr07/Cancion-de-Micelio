# Expectativas Pre-Registradas y Sello — Nivel 4 (`No-Op Demostrado`: `P0_NIVEL3_BASE`)

**Rama de auditoría:** `origin/lilyth-auditoria-cota`
**Commit base (cierre de Nivel 3):** `241ae8eec6616c582ef6c861359ec01062a90788`
**Rama `master` / `HEAD` intacta:** `abf6ce41197888b9312999fc4686ca58c608edcd`
**Directorio exclusivo del bloque:** `PILOT/nivel4_admision_propia/` (sin ediciones en `PILOT/inc/ark/lilyth/*.hpp`, `ark_cxx/base_auditada/`, `etapa_lh/`, ni sellos previos)

---

## 1. Cadena de sellos predecesores intactos (SHA-256)

1. `nivel1_nan_huecos/EXPECTATIVAS_COTA_NIVEL1.md` — `54aa457e9dc69ddec1a982f4904e2d71e354f4ae9ea2dada4aba77c7f52eafa3`
2. `nivel1_nan_huecos/EXPECTATIVAS_FASE4.md` — `49fe152d00248465489ca2c561647d61033197f2543e819f179bec25e4806ef1`
3. `nivel1_nan_huecos/EXPECTATIVAS_FASE4_ADENDA_V2.md` — `32ad43d217d9a792df5c80a23c1a968c1cdb5f9ed80f13814ce4479c5d332e1c`
4. `nivel1_nan_huecos/EXPECTATIVAS_FASE4_ADENDA_V2_1.md` — `f9b0aed3ca2f3117789e672d7770647bcba3f1c53c474a03fc3653d9245f682d`
5. `nivel1_nan_huecos/EXPECTATIVAS_FASE4_ADENDA_V2_2.md` — `7c63f1dcbdfdd14e5ec11f65a1df08ef9cabf55b2b339b86ce39092bd1b21242`
6. `nivel2_horizonte_dinamico/EXPECTATIVAS_NIVEL2.md` — `2387e846df90c7b87bb90d70cbd7d32b8d54ca85ea269386464a260fe6187d41`
7. `nivel3_via_b_parada/EXPECTATIVAS_NIVEL3.md` — `47fd5587cf6ded1682b1c4dd4887500f7a84678412e240f5a71dce18a90d9a94`

---

## 2. Evidencia previa de la Sonda Empírica Multiescala (Fase 4.0)

Antes de este sello se ejecutó la sonda `nivel4_admision_propia/sonda_admision_propia_via_b.cpp` (`sonda_admision_propia_via_b.log`) con semilla exclusiva `SEMILLA_SONDA_NIVEL4 = 2026112500ULL` (`32` hilos, `19.52 s`, `min_escalas=3` en los 3 bancos), comparando la línea base `P0` (`P0_NIVEL3_BASE`) frente a las variantes candidatas de admisión de `propia` (`P1a`, `P1b`: exigir coincidencia dura de `escala_origen == escala_real`; `P2`: exigencia dura de excedencia `g_disparo >= g_real` (doble condicionamiento); `P3`: ambas condiciones simultáneas):

1. **`M1_AR1_3ESC`** (`24` semillas, `158` eventos reales: `j=0: 109`, `j=1: 23`, `j=2: 26`, ventana `t_max_ef = [3, 13]`):
   - **`P0`**: `z_S = +0.589` (`S_medio = +0.01812`, `ee = 0.03077`), `p_det = 0.4742`, `survS3 = 0.9810`, `cobS3 = 0.7595`, `Eg_dif = -0.011`, `esc== = 89.0%`.
   - **`P1a` / `P1b`**: empeoraron la calibración (`z_S = +0.772` / `+0.910`) y recortaron la cobertura en escala gruesa `j=2` de `0.8846` a `0.5769`.
   - **`P2` / `P3`**: colapsaron `p_det` (`0.2121` / `0.1969`) y `cobS3` (`0.3291` / `0.3038`).
2. **`M2_ARFIMA_3ESC`** (`16` semillas, `165` eventos reales: `j=0: 66`, `j=1: 52`, `j=2: 47`, ventana `t_max_ef = [3, 33]`):
   - **`P0`**: `z_S = +0.066` (`S_medio = +0.00292`, `ee = 0.04414`), `0/15` rechazos, `p_det = 0.2107`, `cobS3 = 0.3576`, `Eg_dif = -0.052`, `esc== = 72.9%`.
   - **`P1a` / `P1b`**: recortaron severamente la cobertura en `j=2` de `0.3191` a `0.1064` / `0.1277`.
   - **`P2` / `P3`**: colapsaron `cobS3` (`0.0970` / `0.0909`).
3. **`M3_DIAG_GRUESAS`** (`N_op = 120000`, `N_mc = 200` realizaciones por evento, `t_max_ef = 5`):
   - Eventos con `j_real = 2` dispararon legítimamente en escala fina/intermedia `j'=0 (29.3%)`, `j'=1 (29.5%)` y `j'=2 (41.2%)`.
   - **`P1a` / `P1b`**: redujeron `p_det` a menos de la mitad (`0.5856 -> 0.2411 / 0.2428`) sin mejorar `S` (`+0.0530` frente a `+0.0534 / +0.0530`).
   - **`P2` / `P3`**: colapsaron `p_det` a `0.3928` / `0.1411`.

---

## 3. Decisión Formal de Nivel 4: Cierre por "No-Op Demostrado"

Con base en la evidencia empírica anterior, el **Nivel 4 se cierra sin modificar `PILOT/inc/ark/lilyth/nulo_lilyth.hpp`** (`P0_NIVEL3_BASE` permanece intacta):

1. **No se introduce filtrado duro por `escala_origen == escala_real` (`P1a`/`P1b`)** en la Vía (b), dado que destruye la detección y cobertura en escalas gruesas (`j=1,2`) cuando las excursiones nulas disparan causalmente en escalas más finas (`j' < j`), sin mejorar ni la calibración ni la simetría de $w_{\text{ancla}}$.
2. **No se introduce la exigencia de excedencia `g_disparo >= g_real` (`P2`/`P3`)**, dado que introduce sesgo por doble condicionamiento sobre la magnitud del evento real y colapsa severamente `p_det` y `cobS3`.
3. **Regla 6 estricta (Opt-in):** `ConfiguracionLilyth{}` mantiene `usar_enmienda_nivel2 == false` y `usar_enmienda_nivel3 == false` por omisión.

---

## 4. Protocolo Falsable Pre-Registrado (`SEMILLA_CIEGA_NIVEL4 = 20261125ULL`)

La batería confirmatoria de cierre (`nivel4_admision_propia/test_nv4_admision_propia.cpp`) utiliza la semilla ciega inédita **`SEMILLA_CIEGA_NIVEL4 = 20261125ULL`**, disjunta de la usada en Nivel 2 (`20261123ULL`), Nivel 3 (`20261124ULL`) y la sonda Fase 4.0 (`2026112500ULL`).

### `NV4.1` — Prueba integral E2E (`ObservadorLilyth::completar_nulos()`) bajo `AR(1)` (`phi = 0.80`) con 3 escalas
- **Configuración:** `24` semillas (`s = 0..23`), `n_cal = 12000`, `n_op = 12000`, `ObsDelta = ObservadorSaltos<24,64,120>`, `ObsLil = ObservadorLilyth<ObsDelta,8,64>`, `FLUJO`, `SUMA`, `arl0 = 450.0`, `escalas_max = 3`, `b_max_nulo = 8`, `hilos_nulo = 1` (paralelismo externo por semilla), `b_mc_resumen = 1999`, `INTERCAMBIABLE`, `forzar_via_b = true`, `usar_enmienda_nivel2 = true`, `usar_enmienda_nivel3 = true`.
- **Derivación determinista de semillas:**
  - `sem_cal  = splitmix64(20261125ULL ^ (0x9E3779B97F4A7C15ULL * (s + 1)))`
  - `sem_op   = splitmix64(20261125ULL ^ (0xBF58476D1CE4E5B9ULL * (s + 101)))`
  - `sem_nulo = splitmix64(20261125ULL ^ (0x94D049BB133111EBULL * (s + 211)))`
- **Criterios de aceptación (`FALLAS NV4.1 == 0`):**
  1. `ok_esc`: `cal == 24/24`, `min_esc == 3`, `k_med_j[0..2] > 0`, `n_sem_eval >= 23`.
  2. `ok_cal`: `|z_S| <= 2.0`, `rech(p < 0.05) <= 3/24`.
  3. `ok_cov`: `p_det >= 0.38`, `surv_s3 >= 0.90`, `cob_s3 >= 0.60`, `solapes == 0`.
  4. `ok_j12`: `p_det_j12 >= 0.42` y `escala_cruzada_j12 = n_nulos_b_j12 - escala_coincide_b_j12 > 0`.

### `NV4.2` — Prueba integral E2E (`ObservadorLilyth::completar_nulos()`) bajo `ARFIMA(1, d=+0.35, 0)` con 3 escalas
- **Configuración:** `phi1 = 0.40`, `d = +0.35`, `mu = 50.0`, `16` semillas (`s = 0..15`), `n_cal = 12000`, `n_op = 8000`, `ObsLil128 = ObservadorLilyth<ObsDelta,8,128>`, `ARFIMA_WHITTLE`, `d_whittle_cruda = 0.35`, `arl0 = 450.0`, `escalas_max = 3`, `b_max_nulo = 4`, `hilos_nulo = 1`, `b_mc_resumen = 1999`, `INTERCAMBIABLE`, `forzar_via_b = true`, `usar_enmienda_nivel2 = true`, `usar_enmienda_nivel3 = true`.
- **Derivación determinista de semillas:**
  - `sem_cal  = splitmix64(20261125ULL ^ (0x517CC1B727220A95ULL * (s + 1)))`
  - `sem_op   = splitmix64(20261125ULL ^ (0xD6E8FEB86659FD93ULL * (s + 307)))`
  - `sem_nulo = splitmix64(20261125ULL ^ (0xA0761D6478BD642FULL * (s + 409)))`
- **Criterios de aceptación (`FALLAS NV4.2 == 0`):**
  1. `ok_esc`: `cal == 16/16`, `min_esc == 3`, `k_med_j[0..2] > 0`, `n_sem_eval >= 14`.
  2. `ok_cal`: `|z_S| <= 2.0`, `rech(p < 0.05) <= 2/16`.
  3. `ok_cov`: `p_det >= 0.16`, `cob_s3 >= 0.25`, `solapes == 0`.
  4. `ok_j12`: `p_det_j12 >= 0.15` y `escala_cruzada_j12 = n_nulos_b_j12 - escala_coincide_b_j12 > 0`.

### `NV4.3` — Contraste pareado ciego (`P0..P3`) y verificación de Regla 6
- **Configuración:** `AR(1)` con `phi = 0.80`, `n_cal = 12000`, `n_op = 120000`, `J_ = 3`, hasta `16` eventos cerrados medibles por cada escala `j = 0, 1, 2`, `N_mc = 150` realizaciones acopladas con base `SEMILLA_CIEGA_NIVEL4 + 9001ULL`, evaluando las definiciones pareadas de sonda (`P0`: sin filtro; `P1a`: `escala_origen == escala_real`; `P2`: excedencia `g_disparo >= g_real`; `P3`: `escala_origen == escala_real && g_disparo >= g_real`; `Eg_dif = E[g_disparo - g_real]`).
- **Criterios de aceptación (`FALLAS NV4.3 == 0`):**
  1. `ok_p2`: `p_det_global(P2) <= 0.70 * p_det_global(P0)` y `Eg_dif(P2) > Eg_dif(P0) + 0.05`.
  2. `ok_p1`: `p_det_j2(P1a) <= 0.75 * p_det_j2(P0)` y `p_det_j2(P3) <= 0.50 * p_det_j2(P0)`.
  3. `ok_r6`: `!ConfiguracionLilyth{}.usar_enmienda_nivel2 && !ConfiguracionLilyth{}.usar_enmienda_nivel3`.

---

## 5. Campaña de Mutación de Regresión (`mutaciones_nivel4.py`)

El arnés aislado `nivel4_admision_propia/mutaciones_nivel4.py` opera sobre copias temporales de `PILOT/` en `%TEMP%` (sin tocar el árbol real) ejecutando el testigo completo y las dos regresiones prohibidas sobre el operador real (`test_nv4_admision_propia.exe --solo-operador`):

1. **`MN4.0` (Testigo intacto):** ejecuta `test_nv4_admision_propia.exe` completo (`NV4.1`–`NV4.3`); debe finalizar con `RC == 0` (`FALLAS NV4: 0`).
2. **`M_N4_1` (Reintroducir filtro duro `escala_origen == escala_real`):** modifica `ejecutar_inyecciones_en_realizacion` en copia temporal de `nulo_lilyth.hpp` para exigir coincidencia de escala cuando `usar_enmienda_nivel3` está activa. Debe **FALLAR** en `test_nv4_admision_propia.exe --solo-operador` (`RC != 0`, colapso de `escala_cruzada_j12` / caída de `p_det_j12` y cobertura).
3. **`M_N4_2` (Reintroducir exigencia de excedencia `g_disparo >= g_real`):** inyecta `g_real` en la cadena `ObservadorLilyth::completar_nulos()` -> `ejecutar_via_b_evento()` -> `ejecutar_inyecciones_en_realizacion()` en copias temporales y exige `propia->ev.g_disparo >= g_real` para admitir bajo `usar_enmienda_nivel3`. Debe **FALLAR** en `test_nv4_admision_propia.exe --solo-operador` (`RC != 0`, caída severa de `p_det` y cobertura en `NV4.1`/`NV4.2`).

Cualquier falla en `NV4.1`–`NV4.3`, cualquier mutante sobreviviente (`VIVO`), o cualquier alteración de identidad/sellos bloquea el cierre y la subida a `origin/lilyth-auditoria-cota`.
