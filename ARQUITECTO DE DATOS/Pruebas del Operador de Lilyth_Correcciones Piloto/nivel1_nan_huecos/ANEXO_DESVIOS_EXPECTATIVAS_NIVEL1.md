# Anexo Declarado de Desvíos al Sello y Especificación de Cierre — Nivel 1 (`ObservadorLilyth`)

**Fecha:** 2026-11-17  
**Documento sellado original (intacto):** `EXPECTATIVAS_COTA_NIVEL1.md` (`sha256: 54aa457e9dc69ddec1a982f4904e2d71e354f4ae9ea2dada4aba77c7f52eafa3`, `SELLO_EXPECTATIVAS.txt`).  
**Propósito de este anexo:** Declarar con transparencia total los desvíos entre el documento sellado y la implementación de los Puntos 2 a 7 señalados en la segunda auditoría de Opus 5.5, sin modificar ni borrar el sello original, y fijar las reglas exactas para cerrar los 6 pendientes restantes sin redundancia.

---

## 1. Desvíos Declarados respecto a `EXPECTATIVAS_COTA_NIVEL1.md`

### Desvío D1 — Parámetros de los mundos de réplicas $H_0$ / $H_1$ (`ParametrosReplicaH0H1`)
- **Lo sellado en Sección 1.3:** Choques puntuales con $a = 8.0\,\sigma_\varepsilon$, componente lenta $z_\tau = 4.0 \times 0.97^\tau$ y $R = 200$ réplicas.
- **Lo ocurrido en el código:**
  1. En la primera corrida de `test_punto6_wilson_h0_h1.cpp` se puso `n_cal = 2048`, que era menor al mínimo de calibración de `ObservadorSaltos` ($\lceil 1/\alpha \rceil N_{\text{REG}} = 2400$, más $480$ pasos de maduración de `PurgaDiurna` en `M3`, total $2880$), dando `0` series válidas. Subir `n_cal = 3200` fue una corrección mecánica necesaria.
  2. Sin embargo, en el mismo cambio se modificaron `amp_choque_sigma` de `8.5` a `12.0` y `amp_lenta_sigma` de `4.0` a `4.5` sin declararlo, y en `test_cota_corrupcion_nivel1.cpp` se corrió con `80` réplicas por omisión (`40` en `mutaciones_nivel1.py`) en vez de `200`.
- **Corrección para el cierre:** En `ParametrosReplicaH0H1` se mantiene `n_cal = 3200` (requisito de `ObservadorSaltos` + `PurgaDiurna`) y se restauran exactamente los parámetros sellados del choque ($a = 8.0\,\sigma_\varepsilon$) y de la cola lenta ($z_\tau = 4.0\,\sigma_\varepsilon \times 0.97^\tau$), o bien se reportan lado a lado con `HorizonteLilyth::CUANTIL_RESOLUCION` + `PoolNuloLilyth::DERIVADO` (más `b_max_nulo = 128` para que la vía `(a)` tenga población suficiente para el cuantil $1-\alpha/8$, como documentó Code en `MEDIDO_HORIZONTE.md`) de modo que la potencia base en $q = 0$ supere el $50\%$ y `RetPot_H1` no sea ruido.

### Desvío D2 — Condición `CN1.1` (Identidad bit a bit contra la cabecera original)
- **Lo sellado en Sección 4 (`CN1.1`):** Con la política piloto en sus valores por omisión (`0`), los volcados de `ObservadorLilyth` en `M1`, `M2` y `M3` deben ser idénticos bit a bit a la cabecera original sin tocar (`observador_lilyth_original.hpp`).
- **Lo ocurrido en el código:** `test_cota_corrupcion_nivel1.cpp` solo comprobó que los tres `enum` de `ConfiguracionLilyth` valieran `0` por omisión y que hubiera al menos una ranura censurada ante un hueco, sin compilar ni comparar contra `observador_lilyth_original.hpp`.
- **Corrección para el cierre:** Un ejecutable de volcado hexadecimal (`volcado_identidad_nivel1.cpp`, análogo a `volcado_ab.cpp` de Code) se compila dos veces:
  1. Contra `observador_lilyth_original.hpp` (la copia exacta de `ark_cxx/base_auditada/`), y
  2. Contra `observador_lilyth.hpp` (la cabecera piloto con valores por omisión),
  volcando en formato `%a` (hexadecimal IEEE-754 exacto) todos los campos de las ranuras cerradas, `ficha_nulo_{cruda,purgada}` y `resumen_armonico_{cruda,purgada}` sobre `M1`, `M2` y `M3` (tanto con `q = 0` como con huecos inyectados). `CN1.1` pasa **sí y solo si ambos volcados son idénticos byte a byte (`0` líneas distintas)**.

### Desvío D3 — Condición `CN1.4` y Mutante `MN1.3` (Exclusión de ranuras corruptas/censuradas)
- **Lo sellado en Sección 4 (`CN1.4` y `MN1.3`):** Con un $25\%$ de los choques corrompidos en `Z2` ($r = 3$), excluir esas ranuras de $K$ mantiene el límite inferior de Wilson del rechazo de Holm $\le \alpha = 0.05$, mientras que el mutante `MN1.3` (incluir ranuras `CENSURADA_DATO_INVALIDO` / `NO_MEDIDA` codificadas en $T_{\max}+1$ dentro de los estratos de `INTERCAMBIABLE`) debe ser cazado por `CN1.4`.
- **Lo ocurrido en el código:**
  1. En `CN1.4` solo se corrompió 1 evento en `Z2` ($r=1$) y 2 en `Z3` por serie sobre 40–80 réplicas.
  2. Además, en `observador_lilyth.hpp`, cuando `MN1.3` quitó `CENSURADA_DATO_INVALIDO` del `continue;` en `evaluar_nulos_representacion`, `evaluar_una_referencia` evaluaba esa ranura cuyo `rama.tau_rest_consolidado` es `NaN` (que `codificar_tau` convierte en `T_max + 1`). Con solo 2 eventos cortados en `Z3` de 6, y con nulos de vía `(b)` que también tenían censuras en `T_max + 1`, el desplazamiento de $S$ no alcanzó a subir el `Wilson lo` de 40 réplicas por encima de $0.05$, dejando a `MN1.3` **VIVO**.
- **Corrección para el cierre:**
  - En `CN1.4` se aplica exactamente el escenario sellado: en un mundo $H_0$ limpio (donde las trayectorias reales restituyen rápido, $\tau_{\text{rest}} \ll T_{\max}$, con `0` censuras en el nulo de calma), se corrompe el $25\%$ de los eventos en `Z2` ($r = 3$) y otro $50\%$ de los eventos en `Z3` ($r = 1$ tras la apertura, antes de restituir).
  - En la cabecera sin mutar, las ranuras con `Z2` corrupto en vía `(b)` quedan con `evaluado = false` y las de `Z3` con `CENSURADA_DATO_INVALIDO` se saltan en `evaluar_nulos_representacion`, manteniendo el rechazo de Holm bajo $H_0$ con `Wilson lo <= 0.05`.
  - En el mutante `MN1.3` (que no excluye `CENSURADA_DATO_INVALIDO` en `evaluar_nulos_representacion`), esos eventos reales cortados en `Z3` entran al estrato codificados en $T_{\max} + 1$ contra nulos que restituyeron en $\tau \ll T_{\max}$, dando $u_{\text{PIT}} \approx 1.0$ y disparando el falso rechazo bajo $H_0$ muy por encima de $\alpha = 0.05$ (`Wilson lo > 0.05`), cazando a `MN1.3` de forma determinista (`6/6` mutantes cazados).

### Desvío D4 — Condición `CN1.5` (Derivación algorítmica de cotas y monotonía)
- **Lo sellado en Sección 4 (`CN1.5`):** Las cotas $(q_{\max}, r_{\max})$ se calculan algorítmicamente recorriendo las tablas medidas bajo `ESTRICTO` y `DULCE`.
- **Lo ocurrido en el código:**
  1. Se metieron en `cn1_5_ok` comparaciones contra los números que habían salido en los Puntos 2 a 4 (`r_max_z1_borde_estricto == 0 && r_max_z1_borde_dulce >= 2...`).
  2. El bucle que calculaba `q_max` y `r_max` hacía `if (cmp.pasa_estricto) q_max = q;` **sin exigir monotonía** (por eso en mercado quieto con $R=200$, donde $q=1\%$ pasaba, $q=2\%$ fallaba con `RetPot = 0.82` y $q=5\%$ volvía a pasar con `0.89`, el bucle sin freno habría reportado $5\%$ en vez de detenerse en el primer fallo en $1\%$).
- **Corrección para el cierre:**
  - Toda cota $q_{\max}$ o $r_{\max}$ se calcula con **regla de parada monótona estricta (`break` al primer fallo)**: es el mayor $q$ (o $r$) tal que **todos los niveles desde $0$ hasta $q$ (o $r$) inclusive** pasan el régimen (`ESTRICTO` o `DULCE`).
  - `CN1.5` verifica únicamente invariantes estructurales a priori:
    1. Monotonía de regímenes: para toda zona, $\text{cota}_{\text{ESTRICTO}} \le \text{cota}_{\text{DULCE}}$;
    2. Coherencia de `Z3a` bajo reloj físico: en $r = 1$, el sesgo de $\tau_{\text{rest}}$ con `RELOJ_FISICO_TOPE_R` es menor en valor absoluto que con `MUTANTE_COMPRIME_TAU`;
    3. Ninguna cota se extrapola más allá de la malla efectivamente corrida en el binario.

### Desvío D5 — Acoplamiento con el Detector Delta (`ObservadorSaltos`) y Vía `(b)` en `6B`
- **Lo observado por Opus (Defectos técnicos 1 y 6):**
  1. En `6B`, `ObservadorSaltos` corrió sobre la serie limpia mientras solo `ObservadorLilyth` veía los `NaN`. Eso mide el efecto aislado sobre Lilyth pero oculta cuántos eventos dejan de abrir (o se abren en falso) cuando `ObservadorSaltos` también recibe los `NaN`.
  2. Además, con `forzar_via_b = false` en `6B`, casi todos los eventos iban por vía `(a)`, por lo que `Z1` y `Z2` apenas intervenían en `6B`.
- **Corrección para el cierre:** En `6B` se corren y reportan explícitamente los dos modos:
  - **`6B.1` (Flujo extremo a extremo real: `ObservadorSaltos` + `ObservadorLilyth` ven ambos los `NaN`):** mide la pérdida real de aperturas del sistema completo y el falso rechazo bajo $H_0$.
  - **`6B.2` (Flujo con `forzar_via_b = true`):** obliga a todos los eventos sobrevivientes a pasar por `Z1` (`hist_prev_est`) y `Z2` (`perfil_est`) en la vía `(b)`, ejercitando simultáneamente las políticas de `Z1`, `Z2` y `Z3`.
