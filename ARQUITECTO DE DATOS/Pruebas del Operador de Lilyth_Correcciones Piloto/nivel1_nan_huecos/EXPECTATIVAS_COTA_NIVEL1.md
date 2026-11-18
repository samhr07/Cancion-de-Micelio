# Expectativas Selladas — Nivel 1: Protocolo de Inputs y Cota de Corrupción (`ObservadorLilyth`)

**Fecha de sello:** 2026-11-17  
**Autores:** Arquitecto de Datos (Argon) bajo auditoría de Jefe de Arquitectos de Código (Opus 5.5) y orden del Comandante.  
**Regla de inmutabilidad:** Este documento se sella con SHA-256 en `SELLO_EXPECTATIVAS.txt` **antes** de escribir el código corregido de los puntos 2 al 7. Ningún umbral, fórmula, semilla ni criterio puede modificarse después del sello. Cualquier condición o expectativa que falle en la corrida real se reporta como `[FALLA]` o `[NO SE CUMPLE]` sin alterar el código ni los topes.

---

## 1. Diseño Experimental Conectado al Operador Base y `ark/dataset/`

1. **Cero fórmulas sustitutas:**
   - Todas las innovaciones se extraen exclusivamente con `ark::lilyth::innovaciones()` y la desviación típica residual `sd_innov` del `ModeloLinealNulo` real (`ark/lilyth/modelo_lineal_nulo.hpp`).
   - La vía `(b)` se ejecuta exclusivamente con `ark::lilyth::detalle_nulo::ejecutar_via_b_evento()` y `ObservadorLilyth::completar_nulos()` (`ark/lilyth/nulo_lilyth.hpp`, `ark/lilyth/observador_lilyth.hpp`).
   - La trayectoria post-apertura $\eta(\tau)$ se evalúa exclusivamente a través de `ObservadorLilyth::avanzar()` instrumentado en la copia piloto `inc/ark/lilyth/observador_lilyth.hpp`.
   - El contraste estadístico y el ajuste de multiplicidad se leen directamente de `ficha_nulo_{cruda,purgada}()` (`contraste_onsager`, `contraste_onsager_pool_loo`) y `resumen_armonico_{cruda,purgada}()` (`p_holm[1]`, `familia_holm_valida`).
2. **Emparejamiento estricto por evento físico:**
   - Un evento de la corrida con huecos/corrupción se compara contra la corrida limpia ($q = 0, r = 0$) **únicamente si coincide su clave exacta `(t_apertura, escala_origen, signo)`**.
   - Los eventos limpios que no aparecen o quedan en `CENSURADA_DATO_INVALIDO` / `NO_MEDIDA_*` se contabilizan en la tasa de pérdida/censura $\Delta\text{Censura}$, nunca emparejándolos por índice de vector contra otro choque distinto.
3. **Mundos sintéticos con tamaño nominal bajo $H_0$ y potencia bajo $H_1$ ($R = 200$ réplicas por mundo):**
   - **`M1_H0` (`FLUJO AR(1)`, $\phi = 0.80$, `ModoNuloLineal::AR_LEVINSON_DURBIN`):** Choques puntuales de 1 paso ($a = 8\,\sigma_\varepsilon$) que decaen con la propia ley $\text{AR}(1)$ del mundo ($z_\tau = 0$).
   - **`M2_H0` (`FLUJO ARFIMA(1, d=0.25, 0)`, $\phi_1 = 0.40$, `ModoNuloLineal::ARFIMA_WHITTLE`):** Choques puntuales de 1 paso ($a = 8\,\sigma_\varepsilon$) que decaen con la propia memoria fraccionaria $(1-L)^{-d}(1-\phi_1 L)^{-1}$ del mundo ($z_\tau = 0$).
   - **`M3_H0` (`NIVEL AR(1)` $\phi = 0.65$ + `PurgaDiurna` 24 fases, `ModoNuloLineal::AR_LEVINSON_DURBIN`):** Choques puntuales de 1 paso modulados por la marea diurna de 24 fases ($z_\tau = 0$).
   - **`M1_H1`, `M2_H1`, `M3_H1` (Mundos con mecanismo lento para medir potencia):** Mismos mundos base añadiendo tras cada choque la componente lenta aditiva $z_\tau = 4.0 \times 0.97^\tau$.
   - **Semillas deterministas:** Para la réplica $\text{rep} \in \{0, \dots, R-1\}$, semilla `20261117ULL + 10007ULL * rep`.

---

## 2. Los Dos Regímenes de Umbrales Sellados: Versión Estricta vs. Versión Dulce

Cada celda $(q, r, \text{Zona}, \text{Política})$ será evaluada y reportada simultáneamente bajo dos regímenes de tolerancia:

| Métrica (medida sobre el operador real y eventos emparejados por `t_apertura`) | **Versión Estricta (`ESTRICTO`)** *(Metrología fina al nivel $\alpha = 0.05$)* | **Versión Dulce (`DULCE`)** *(Tolerancia operativa $< 50\%$ de error)* |
|---|---|---|
| **1. Error relativo de innovación en el perfil** $\text{RMSE}_{\text{rel}}(e_{\text{iny}}) = \frac{\text{RMSE}(e_{\text{iny}}^{\text{cont}} - e_{\text{iny}}^{\text{limp}})}{\text{sd\_innov}}$ | **$\le 0.15$** ($15\%$ de $\sigma_\varepsilon$) | **$< 0.50$** ($< 50\%$ de $\sigma_\varepsilon$) |
| **2. Caída de probabilidad de detección inicial en vía `(b)`** $\Delta p_{\text{det}} = p_{\text{det,ini}}(0) - p_{\text{det,ini}}(q,r)$ | **$\le 0.05$** ($\le 5$ puntos porcentuales) | **$\le 0.25$** (conserva $\ge 75\%$ de la detección inicial) |
| **3. Desvío del estadístico de resumen $S$** $\frac{|\bar{S}(q,r) - \bar{S}(0)|}{\overline{\text{ee}}_S}$ | **$\le 0.25$** ($\le \frac{1}{4}$ error estándar) | **$\le 1.00$** ($\le 1$ error estándar) |
| **4. Sesgo medio en tiempo de restitución consolidado** $|\overline{\tau_{\text{rest}}^{\text{cont}} - \tau_{\text{rest}}^{\text{limp}}}|$ | **$\le 0.25$ pasos** | **$\le 1.50$ pasos** |
| **5. Tasa de pérdida/censura adicional de eventos reales** $\Delta\text{Censura} = 1 - \frac{K_{\text{evaluados}}(q,r)}{K_{\text{evaluados}}(0)}$ | **$\le 0.10$** ($\le 10\%$ de eventos perdidos) | **$< 0.50$** (conserva $> 50\%$ de los eventos) |
| **6. Falso rechazo bajo $H_0$ (`p_holm[1] < 0.05`)** en las $R$ réplicas | Límite inferior de Wilson al $95\%$ **$\le \alpha = 0.05$** | Límite inferior de Wilson al $95\%$ **$\le 2\alpha = 0.10$** |
| **7. Retención de potencia bajo $H_1$** $\frac{\text{Potencia}(q,r)}{\text{Potencia}(0)}$ | **$\ge 0.85$** ($\ge 85\%$ de la potencia limpia) | **$\ge 0.50$** ($\ge 50\%$ de la potencia limpia) |

---

## 3. Batería de Experimentos por Zona (Malla de Tasas $q$ y Rachas $r$)

### Experimento E1 — `ILIQUIDEZ_TRANSITORIA` (Mercado Quieto) desde `ark/dataset/`
- **Entrada:** Flujo de `ark::dataset::BarraCausal` donde una fracción $q \in \{0\%, 0.5\%, 1\%, 2\%, 5\%, 10\%\}$ de los pasos fuera del choque son barras de mercado quieto (`es_imputado = true`, sin `FLAG_FILA_VACIA`, con el mismo nivel $x_t = x_{t-1}$ en `NIVEL`).
- **Comparación A/B sobre la misma trayectoria física:**
  - `ACTUAL`: `valido = !b.es_imputado` (comportamiento actual de `observador_saltos_filas.hpp`).
  - `PROTOCOLO`: `valido = true` cuando `clasificar_barra_dataset(b) == ClasePasoLilyth::MERCADO_QUIETO`.
- **Qué se mide:** Las 7 métricas de la Sección 2 bajo `ESTRICTO` y `DULCE` en `M1`, `M2` y `M3`.

### Experimento E2.Z1 — Corrupción Física en `hist_prev_est` (`t < t_cambio`)
- **Malla:** Distancias al inicio del choque $\Delta t = t_{\text{cambio}} - t_{\text{fin\_hueco}} \in \{1, 2, 3, 5, 10, 20\}$ y longitudes de racha $r \in \{1, 2, 3, 5, 8\}$, tanto en `M1` (`AR_LEVINSON_DURBIN`) como en `M2` (`ARFIMA_WHITTLE`).
- **Políticas comparadas:**
  1. `COMPRESION_BASE`: El código real actual de `copiar_perfil_estacionario` (salta `hist_val_est == 0` con `continue`, sin inyectar valores artificiales).
  2. `CORTE_SUFIJO`: Trunca `hist_prev_est` para conservar únicamente los pasos contiguos válidos inmediatamente anteriores a `t_cambio`.
  3. `PUENTE_CAUSAL`: Sustituye el paso inválido por la predicción causal a 1 paso del `ModeloLinealNulo` real hasta $r \le r_{\max}$.
- **Expectativa declarada de antemano (basada en el hallazgo honesto de la 1ª corrida):**
  - En $\Delta t = 1$ (hueco pegado a `t_cambio - 1`), **se espera que las tres políticas fallen el régimen `ESTRICTO` ($\text{RMSE}_{\text{rel}} > 0.15$)** porque la innovación $\varepsilon_{t_{\text{cambio}}-1}$ es irrecuperable, pero que para $r \le 2$ entren en el régimen **`DULCE` ($\text{RMSE}_{\text{rel}} < 0.50$)**.
  - Para $\Delta t > p$ en `AR(p)` (`M1`), `CORTE_SUFIJO` y `PUENTE_CAUSAL` deben dar $\text{RMSE}_{\text{rel}} = 0.0000$ exacto (cumpliendo `ESTRICTO` y `DULCE`).
  - En `M2` (`ARFIMA_WHITTLE`, memoria larga $d = 0.25$), para $\Delta t \ge 3$ y $r \le 2$, `PUENTE_CAUSAL` debe dar menor $\text{RMSE}_{\text{rel}}$ que `COMPRESION_BASE` (que desfasa la convolución fraccionaria $(1-L)^d$) y que `CORTE_SUFIJO` (que trunca la cola fraccionaria).

### Experimento E2.Z2 — Corrupción Física en `perfil_est` (`[t_cambio, t_apertura]`)
- **Malla:** Rachas $r \in \{0, 1, 2, 3\}$ dentro de `[t_cambio, t_apertura]` evaluadas con `ejecutar_via_b_evento()` real.
- **Políticas comparadas:** `COMPRESION_BASE` vs `PUENTE_CAUSAL` vs declarar `NO_MEDIDA_PERFIL_CORRUPTO` (excluida de $K$).
- **Expectativa declarada:** Con $r \ge 1$ en `[t_cambio, t_apertura]`, el choque inyectado pierde energía real ($\text{RMSE}_{\text{rel}} \ge 0.50$), fallando `ESTRICTO` (y en $r \ge 1\text{–}2$ también `DULCE`), por lo que $r_{\max, Z2}^{\text{ESTRICTO}} = 0$.

### Experimento E2.Z3 — Corrupción Física en `avanzar()` (`t > t_apertura`)
- **Malla:** Rachas $r \in \{0, 1, 2, 3, 5, 8\}$ posicionadas en dos regímenes:
  - `(Z3a)` Antes de la ventana de cruce de $1/e$ ($\tau \in [2, \tau_{\text{rest}} - 2]$).
  - `(Z3b)` Directamente sobre el inicio del cruce de $1/e$ ($\tau = \tau_{\text{rest}}$ limpio).
- **Políticas comparadas en `ObservadorLilyth::avanzar()`:**
  1. `CORTA_AL_PRIMERO` (actual: `CENSURADA_DATO_INVALIDO` en $r=1$).
  2. `MUTANTE_COMPRIME_TAU` (salta el paso inválido pero incrementa $\tau$ solo en válidos).
  3. `RELOJ_FISICO_SIN_COMPRIMIR` (mantiene $\tau = t - t_a$, escribe `eta[tau-1] = NaN` en el paso inválido y evalúa la consolidación de 3 pasos sobre pasos válidos mientras $r \le r_{\max}$).
- **Expectativa declarada:**
  - En `(Z3a)` (antes del cruce), `RELOJ_FISICO_SIN_COMPRIMIR` preserva $\tau_{\text{rest}}$ con sesgo $0.00$ (**cumple `ESTRICTO` y `DULCE`**), mientras `MUTANTE_COMPRIME_TAU` sesga $\tau_{\text{rest}}$ en $-r$ pasos.
  - En `(Z3b)` (justo encima del cruce de $1/e$), un hueco de longitud $r$ retrasa la consolidación en $+r$ pasos: para $r = 1$, el sesgo medio ponderado por la probabilidad de caer en el cruce debe cumplir `DULCE` ($\le 1.50$ pasos), mientras que para $r \ge 3$ rompe la ventana de consolidación de 3 pasos.

---

## 4. Condiciones Falsables (`CN1.1`–`CN1.6`) y Mutantes de Cabecera (`MN1.1`–`MN1.6`)

### Condiciones (el binario devuelve `0` si y solo si todas pasan)
- **`CN1.1` (Identidad bit a bit en $q = 0, r = 0$):** Con la política piloto en sus valores por omisión, los volcados de `ObservadorLilyth` en `M1`, `M2` y `M3` son idénticos bit a bit a la cabecera original sin tocar.
- **`CN1.2` (Tamaño nominal de los mundos $H_0$ en $q = 0$):** En `M1_H0`, `M2_H0` y `M3_H0` ($R = 200$ réplicas), el límite inferior de Wilson del rechazo de Holm (`p_holm[1] < 0.05`) es $\le \alpha = 0.05$.
- **`CN1.3` (Supervivencia en `M3` bajo `MERCADO_QUIETO`):** Con $q = 2\%$ en `M3_H0`, `PROTOCOLO` produce $0$ `CENSURADA_DATO_INVALIDO` por mercado quieto y conserva $\Delta\text{Censura} = 0$, mientras `ACTUAL` pierde $> 50\%$ de las ranuras.
- **`CN1.4` (Coherencia de exclusión en `NO_MEDIDA` / `CENSURADA_DATO_INVALIDO`):** En `M1_H0` con un $25\%$ de choques corrompidos en `Z2` ($r = 3$), marcar esas ranuras como `NO_MEDIDA_PERFIL_CORRUPTO` (excluidas de $K$) mantiene el límite inferior de Wilson del rechazo de Holm $\le \alpha = 0.05$.
- **`CN1.5` (Derivación algorítmica de cotas):** Todas las cotas $(q_{\max}, r_{\max})$ reportadas al final del log se calculan automáticamente recorriendo las tablas medidas bajo `ESTRICTO` y `DULCE`, sin constantes escritas a mano en el `printf` final.
- **`CN1.6` (Caza de mutantes reales de cabecera $\ge 5/6$):** El script `mutaciones_nivel1.py` aplica reemplazos de texto sobre una copia de `inc/ark/`:
  - `MN1.1`: En `protocolo_inputs_lilyth.hpp`, tratar `MERCADO_QUIETO` como `valido = false` $\rightarrow$ cazado por `CN1.3`.
  - `MN1.2`: En `observador_lilyth.hpp` (`avanzar`), incrementar el índice de `eta` solo cuando `ok == true` (comprimir $\tau$) $\rightarrow$ cazado por sesgo negativo de $\tau_{\text{rest}}$ en `E2.Z3`.
  - `MN1.3`: En `observador_lilyth.hpp` (`evaluar_nulos_representacion`), incluir ranuras `CENSURADA_DATO_INVALIDO` / `NO_MEDIDA` codificadas en `T_max + 1` dentro de los estratos de `INTERCAMBIABLE` $\rightarrow$ cazado por `CN1.4` (inflación de rechazo bajo $H_0$).
  - `MN1.4`: En `protocolo_inputs_lilyth.hpp`, no cortar segmento (`valido = true`) en `FRONTERA_ESTRUCTURAL` (`FIN_DE_SEMANA`) $\rightarrow$ cazado por control de salto de sesión en `NIVEL`.
  - `MN1.5`: En `observador_lilyth.hpp`, activar por omisión el puente/relajación de `Z3` cuando la configuración pide el modo base $\rightarrow$ cazado por `CN1.1` (identidad bit a bit).
  - `MN1.6`: En `observador_lilyth.hpp` (`copiar_perfil_estacionario`), admitir perfiles con `CORRUPCION_O_CAIDA` ($r \ge 2$) en `Z2` sin marcar `NO_MEDIDA_PERFIL_CORRUPTO` bajo el modo protegido $\rightarrow$ cazado por el control de error de innovación en `E2.Z2`.
