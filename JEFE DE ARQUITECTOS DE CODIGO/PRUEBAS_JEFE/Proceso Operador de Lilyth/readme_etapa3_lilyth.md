# Operador de Lilyth -- Prerregistro y Preexpectativas de la Etapa 3 (`nulo_lilyth.hpp`)

Auditoria Externa (Argon) para el Jefe de Arquitectos de Codigo (Code, Opus 5.5) y el Comandante, 2026-11-06.
Documento escrito y entregado **ANTES** de escribir y correr `nulo_lilyth.hpp` y los controles `S1`-`S3`, para su sellado con SHA-256 conforme a la seccion 3.7 del plano maestro (`Planos de las Etapas 3 a 8 del Operador de Lilyth`).

---

## 1. Cierre de la Etapa 3.0 (E1.21: continuacion contra replica en `observador_lilyth.hpp`)

- **Regla implementada**: cuando el espejo del operador delta se abre (`abre == true`) existiendo al menos una ranura viva (`r.ocupada == true` tras avanzar las ranuras vivas en el paso $t$), la nueva apertura se evalua contra la **ranura viva mas reciente** (aquella con mayor `id_evento` / `t_apertura`):
  - Si `a_abre.signo == ref_viva->signo` **y** `!ref_viva->ha_restituido()` (ninguna de sus ramas activas tiene `tau_rest_consolidado` finito) $\implies$ **continuacion**: no corta ninguna ranura viva, no abre ranura nueva, activa `SalidaPaso.continuacion = true` e incrementa `Contadores.continuaciones`.
  - Si `a_abre.signo != ref_viva->signo` **o** `ref_viva->ha_restituido()` $\implies$ **replica**: activa `SalidaPaso.replica = true`, incrementa `Contadores.replicas`, corta las ranuras vivas (si `cortar_abiertas_con_replica == true`) y abre una ranura nueva.
- **Declaracion de por que se evalua contra la ranura viva mas reciente**: con `cortar_abiertas_con_replica == true` (valor por omision), toda replica cierra las ranuras abiertas previas y a lo sumo hay una ranura viva. Cuando `cortar_abiertas_con_replica == false` pueden coexistir varias ranuras vivas de choques sucesivos; en ese caso, el regimen local del observable inmediatamente anterior a la nueva alarma esta determinado por el ultimo choque admitido (la ranura viva mas reciente). Comparar contra una ranura mas antigua ignoraria un choque intermedio de signo contrario o ya restituido que ya cambio el estado del sistema.
- **Verificacion de la Etapa 3.0**:
  - `etapa3_0_lilyth.cpp`: 18/18 controles en verde (`K1`, `K2`, `K2b`, `K3`, `K4`, `K4b`, `K5`, `K5b`, `K5c`, `K6`, `K6b`, `K7`, `K8`, `K9`, `K10`, `K11`) contra `ObservadorSaltos<24, 64, 120>` real.
  - `mutaciones_etapa3_0.py`: 22/22 mutantes compilados y cazados (`M01`..`M18` de Etapa 1 mas `M19_continuacion_corta`, `M20_continuacion_ignora_signo`, `M21_continuacion_ignora_si_ya_restituyo` y `M22_compara_con_ranura_mas_antigua`).

---

## 2. Diseno y parametros derivados de la Etapa 3 (`nulo_lilyth.hpp`)

1. **Dos poblaciones nulas por representacion activa (`PURGADA`, `CRUDA`)**:
   - **Nulo sin memoria (`iid`)**: sustituto `PERMUTACION` o `SIGNOS` (`ark::metrologia::sustituto(segs, nulo, bloque, mu, rng)`) sobre la historia de calibracion de la representacion; en `NIVEL` se sustituyen los incrementos y se reintegra desde el origen de cada segmento (E1.22).
   - **Nulo lineal (`Onsager`)**: realizaciones del `ModeloLinealNulo` de la Etapa 2 ajustado sobre los segmentos de calma de la representacion (en `NIVEL`, ajustado sobre incrementos y reintegrado).
2. **Via (a), espontaneos**:
   - $B$ realizaciones independientes de la longitud de la historia de calibracion (`familia::realizaciones(m_min, arl0, n, B_max, b_pedido)` con $m_{\min} = \text{familia::m\_min}(\alpha) = 8$ para $\alpha = 0.05$). Si $B$ alcanza `B_max`, se marca `b_topado = true`.
   - Ejecucion mediante `familia::en_paralelo(B, hilos, f)` con un generador `std::mt19937_64` propio por indice de realizacion `b` (sembrado deterministicamente desde `semilla_base ^ mezclar(b)`), garantizando identidad bit a bit entre 1 hilo y $N$ hilos (`N2`).
   - En cada realizacion corre una **copia nueva** del operador delta ya calibrado (`Obs copia_delta = delta_calibrado; copia_delta.reiniciar_flujo();`, Regla 7, control `N3`) y un `ObservadorLilyth` con la misma configuracion y regla E1.21.
   - **Emparejamiento por excedencia (sin bins)**: para cada evento real $k$ con escala de origen $j_k$, signo $s_k$ y estadistico de disparo $g_k$, sus eventos nulos espontaneos son aquellos con $j_{\text{nulo}} = j_k$, $s_{\text{nulo}} = s_k$ y $g_{\text{nulo}} \ge g_k$.
3. **Derivacion de $T_{\max}$ (E1.9)**:
   - Se mide la poblacion de eventos espontaneos del nulo lineal con horizonte `T_MAX` (el tope de plantilla).
   - Si el modelo lineal es invalido por raiz unitaria (`radio_max >= 1.0`) o la integracion es `NIVEL` con incrementos sin memoria ($p = 0, d = 0$), se fija $T_{\max} = W_{\text{ancla}}(j=0) = \text{N\_REG}$ y se marca `t_max_raiz_unitaria = true` (`P10`).
   - En caso estacionario, sea $f_{\text{cens}}$ la fraccion de eventos del nulo lineal que no restituyen dentro de `T_MAX`. Si $f_{\text{cens}} > \alpha$, el cuantil $1 - \alpha$ cae fuera del tope de memoria: se fija $T_{\max} = \text{T\_MAX}$ y se marca `t_max_censurado = true`. Si $f_{\text{cens}} \le \alpha$, se asigna a los censurados el valor `T_MAX` y se toma $T_{\max} = \max(2, \lceil Q_{1-\alpha}(\tau_{\text{rest}}^{\text{primer\_toque}}) \rceil)$ mediante `est::cuantil`.
4. **Via (b), inyeccion generativa (`completar_nulos()`, fuera del ciclo caliente)**:
   - Se activa para todo evento real $k$ cuya poblacion nula lineal en (a) tenga menos de $m_{\min}(\alpha) = 8$ eventos con $g \ge g_k$.
   - Si el modo configurado era `FASE_ALEATORIZADA_FFT`, `ModeloLinealNulo::ajustar` ya cayo a `AR_LEVINSON_DURBIN` con `cayo_a_generativo = true` (`P12`).
   - Se extrae el perfil del evento real $k$ en $[t_{\text{cambio}}, t_a]$ (longitud $l = t_a - t_{\text{cambio}} + 1$) junto con su historia previa y se obtienen sus innovaciones `e_iny` con `innovaciones(...)`.
   - **Ventana de exclusion de una inyeccion en $t_{\text{iny}}$**:
     $$\mathcal{V}(t_{\text{iny}}, l) = \bigl[\,t_{\text{iny}} - (W_{\text{ancla}} + \text{guarda}),\; t_{\text{iny}} + l + \text{T\_MAX}\,\bigr].$$
     Varias inyecciones en una misma realizacion solo se admiten si sus ventanas $\mathcal{V}$ son disjuntas; si dos ventanas candidatas se solapan o si el operador delta dispara un evento espontaneo dentro de $\mathcal{V}(t_{\text{iny}}, l)$ fuera del propio choque inyectado, la inyeccion se descarta y se incrementa `Contadores.solapes`.
   - **Exigencia del mismo disparo**: sobre la realizacion inyectada corre una copia limpia del operador delta y de `ObservadorLilyth`. Solo entra al nulo si `ObservadorLilyth` abre una ranura en $[t_{\text{iny}}, t_{\text{iny}} + l]$ con el mismo signo $s_k$. Si no dispara, se incrementa `no_detectados` (exportando `tasa_no_detectados`).
   - **Regla de parada [DECIDIDO]**: se ejecutan primero $m_{\min}$ inyecciones de sondeo para estimar $\hat{p}_{\text{det}} = k_0 / m_{\min}$ (fijada al cerrar la tanda inicial de $m_{\min}$ inyecciones, sin actualizarla con las inyecciones posteriores). Si $\hat{p}_{\text{det}} = 0$, la busqueda se detiene de inmediato, el evento queda con `sin_nulo_b = true` y se declara en la ficha. Si $\hat{p}_{\text{det}} > 0$, el tope de intentos queda fijado en
     $$M_{\text{tope}} = \min\!\left(M_{\text{max\_mem}},\; \left\lceil \frac{m_{\min}}{\hat{p}_{\text{det}}} \right\rceil\right),$$
     repitiendo hasta reunir $m_{\min}$ eventos detectados o agotar $M_{\text{tope}}$ intentos.
5. **PIT de punto medio y pruebas del resumen [DECIDIDO]**:
   - Para cada evento real $k$ y cada nulo (sin memoria y lineal), con valores nulos $v_1, \dots, v_m$ (donde los no restituidos en $T_{\max}$ toman el valor tope $T_{\max} + 1$):
     $$u_k = \frac{\#\{i : v_i < \tau_k\} + \frac{1}{2}\,\#\{i : v_i = \tau_k\}}{m}.$$
   - Estadistico global sobre los $K$ eventos reales medidos:
     $$S_{\text{real}} = \frac{1}{K} \sum_{k=1}^{K} \left(u_k - \frac{1}{2}\right).$$
   - Distribucion nula de $S$ por Monte Carlo ($B_{\text{MC}} = 1999$ conjuntos de tamano $K$ muestreados de la poblacion nula correspondiente, evaluando el PIT de punto medio de cada elemento muestreado contra el resto de su nulo): valor $p$ bilateral
     $$p = \frac{1 + \#\bigl\{b : |S^{(b)}| \ge |S_{\text{real}}|\bigr\}}{1 + B_{\text{MC}}}.$$

---

## 3. Preexpectativas selladas de los controles S1, S2 y S3 (escritas ANTES de correr)

### S1 -- Replicas de Hawkes sobre AR(1) de calma conocida

- **Mundo sintetico**: serie de `FLUJO` con fondo estacionario $\text{AR}(1)$, $x_t = \phi x_{t-1} + \varepsilon_t$, $\phi = 0.8$ ($\tau_{\text{int}} = \frac{1+\phi}{1-\phi} = 9.0$, semivida teorica de la media $\tau_{1/2} = \frac{\ln 2}{-\ln 0.8} \approx 3.106$ pasos), $\varepsilon_t \sim \mathcal{N}(0, 1 - \phi^2)$ ($\text{sd}(x) = 1$).
- **Llegadas de choques**: proceso de Hawkes univariado en tiempo discreto con intensidad por paso
  $$\lambda(t) = \lambda_0 + \sum_{t_i < t} n \, (1 - e^{-\beta}) e^{-\beta (t - 1 - t_i)},$$
  tasa base $\lambda_0 = 1 / 600$ pasos$^{-1}$, decaimiento del nucleo $\beta = 1 / 12$ pasos$^{-1}$ (vida media del nucleo $\approx 12$ pasos $< T_{\max}$), evaluado en las tres razones de ramificacion $n \in \{0.0, 0.3, 0.6\}$. Cada choque inyecta en las innovaciones un pulso de amplitud $A = 5.5\,\sigma_\varepsilon$ durante $l_{\text{pulso}} = 4$ pasos, alternando el signo de las hijas o asignando signos independientes $\pm 1$ con probabilidad $1/2$ (conocidos en la verdad del generador junto con el identificador de su madre).
- **Que debe salir en cada razon de ramificacion $n \in \{0.0, 0.3, 0.6\}$**:
  1. **(i) Continuaciones y replicas (E1.21) evento por evento**:
     - En **todas** las $n \in \{0.0, 0.3, 0.6\}$, el vector de decisiones `(abrio_ranura, continuacion, replica)` y el desenlace de cada ranura cerrada coinciden **100 % bit a bit** con el motor de referencia independiente alimentado con la verdad del flujo.
     - En $n = 0.0$ (Poisson sin autoexcitacion), la fraccion de ranuras cortadas por replica es baja ($\le 0.15$); en $n = 0.3$ sube a un regimen intermedio ($0.10 \le f_{\text{rep}} \le 0.45$); en $n = 0.6$ es mayor que en $n = 0.0$ ($f_{\text{rep}}(0.6) > f_{\text{rep}}(0.0)$), reflejando el agrupamiento en cascada de las hijas.
  2. **(ii) Curva `CIF_replica(tau)` de Aalen-Johansen contra la incidencia simulada**:
     - En cada $n \in \{0.0, 0.3, 0.6\}$, el estimador `aalen_johansen_cruda().total.cif_replica[tau - 1]` permanece dentro de la banda de $\pm 2.5\,\text{ee}$ binomial/bootstrap respecto de la incidencia acumulada verdadera de corte por replica medida sobre 50 realizaciones independientes del mismo mundo de Hawkes con el mismo $n$.
  3. **(iii) Invariancia de $\Delta\tau_{\text{Onsager}}$ en las tres razones $n \in \{0.0, 0.3, 0.6\}$**:
     - Como entre choques el medio relaja exactamente con el $\text{AR}(1)$ lineal de $\phi = 0.8$ y las replicas solo censuran trayectorias por riesgo competitivo (sin cambiar la fisica de restitucion), en **las tres** razones $n = 0.0$, $n = 0.3$ y $n = 0.6$ el contraste bilateral de Onsager no rechaza la hipotesis nula lineal al nivel $\alpha = 0.05$ ($p_{\Delta\tau,\text{Onsager}} \ge 0.05$ en cada $n$, y $|S_{\text{Onsager}}| \le 2.5\,\text{ee}_{\text{nulo}}$), mientras que $\Delta\tau_{\text{iid}} > 0$ rechaza el nulo sin memoria ($p_{\Delta\tau,\text{iid}} < 0.05$).
  4. **(iv) Estabilidad de la parada de la via (b) entre las tres $n$**:
     - La probabilidad media de deteccion en el nulo generativo $\hat{p}_{\text{det}}(n)$ medida en la via (b) es alta ($\hat{p}_{\text{det}} \ge 0.70$) y estable entre $n = 0.0$, $0.3$ y $0.6$ ($|\hat{p}_{\text{det}}(n) - \hat{p}_{\text{det}}(0)| \le 0.15$), alcanzando $m_{\min} = 8$ nulos detectados por evento sin agotar el tope en los eventos de SNR nominal.

### S2 -- Ventanas solapadas y contaminacion del ancla

- **Longitud exacta de la ventana de guarda de una inyeccion**:
  $$D_{\text{min}} = W_{\text{ancla}} + \text{guarda} + l + \text{T\_MAX}.$$
  Con $j_0 = 0$ ($m = 1$), $W_{\text{ancla}} = 120$, $\text{guarda} = \lceil \tau_{\text{int}} \rceil$, perfil de longitud $l$ y horizonte `T_MAX`, dos inyecciones en $t_1 < t_2$ se solapan si y solo si $\Delta t = t_2 - t_1 \le D_{\text{min}}$.
- **Que debe salir en cada separacion $\Delta t$**:
  1. **Separacion subcritica ($\Delta t = \lfloor D_{\text{min}} / 2 \rfloor < D_{\text{min}}$)**:
     - Al intentar insertar $K_{\text{par}}$ segundas inyecciones a distancia $\Delta t < D_{\text{min}}$ de una primera inyeccion activa en la misma realizacion, el verificador de ventana rechaza el 100 % de las segundas inyecciones: `Contadores.solapes == K_par` exacto, y **0** inyecciones solapadas entran a la poblacion del nulo.
  2. **Separacion en el borde exacto ($\Delta t = D_{\text{min}}$) vs. supercritica ($\Delta t = D_{\text{min}} + 10 > D_{\text{min}}$)**:
     - En $\Delta t = D_{\text{min}}$ (los extremos $t_1 + l + \text{T\_MAX}$ y $t_2 - (W_{\text{ancla}} + \text{guarda})$ coinciden en el mismo paso), se cuenta como solape (`solapes == K_par`).
     - En $\Delta t = D_{\text{min}} + 10 > D_{\text{min}}$ (con umbral alto o sin disparos espontaneos en la ventana), `solapes == 0` por colision entre inyecciones y todas las inyecciones detectadas entran al nulo.
  3. **Efecto fisico de la contaminacion (con guarda vs. sin guarda)**:
     - En un $\text{AR}(1)$ con $\phi = 0.85$ donde se inyecta un choque precedido a distancia corta ($\Delta t = 40 < W_{\text{ancla}} = 120$) por otro choque grande del mismo signo $+1$ cuya cola eleva la ventana del ancla $Y^*$ del segundo choque:
       - **Con la guarda activa**: la inyeccion contaminada se rechaza y la mediana de $\tau_{\text{rest}}$ del nulo limpio coincide con la semivida verdadera de las inyecciones aisladas (diferencia $\le 1.0$ paso).
       - **Sin la guarda (mutante / modo sin filtro de solape)**: el segundo choque entra con $Y^*$ sesgado al alza (menor $dY = Y(t_a) - Y^*$), acortando artificialmente su $\tau_{\text{rest}}$ medido en mas de $1.5$ pasos respecto del nulo limpio.

### S3 -- Integridad de la parada de la via (b)

- **Que debe salir en cada regimen de $p_{\text{det}}$**:
  1. **Regimen intermedio ($0.25 \le p_{\text{det}} \le 0.85$, amplitud calibrada cerca del umbral del GLR)**:
     - Con $p_{\text{det}}$ verdadera medida de forma independiente sobre 500 inyecciones, al ejecutar la parada secuencial sobre $R = 100$ eventos sin tope truncado (o verificando la ley condicional):
       - El numero medio de intentos hasta obtener $m_{\min} = 8$ detecciones es compatible con la media teorica de la Binomial Negativa $m_{\min} / p_{\text{det}}$ dentro de $2.5\,\text{ee}$ ($\text{ee} = \sqrt{m_{\min}(1 - p_{\text{det}}) / (R \, p_{\text{det}}^2)}$).
       - Cuando se aplica el tope derivado $M_{\text{tope}} = \lceil m_{\min} / \hat{p}_{\text{det}} \rceil$ calculado **exclusivamente** con las primeras $m_{\min} = 8$ inyecciones, ningun evento supera $M_{\text{tope}}$ intentos (`intentos_realizados <= M_tope` en el 100 % de los casos).
  2. **Regimen sub-umbral ($p_{\text{det}} = 0$, perfil de amplitud $0.05\,\sigma_0$ invisible para el GLR)**:
     - En las primeras $m_{\min} = 8$ inyecciones se obtiene $\hat{p}_{\text{det}} = 0$.
     - El procedimiento realiza **exactamente** $m_{\min} = 8$ intentos (no sigue intentando ni divide por cero), marca `sin_nulo_b = true`, deja `n_nulos_lineal < m_min` y registra el motivo en la ficha del evento.

---

## 4. Preexpectativas del resto de controles de la Etapa 3 (`P1`-`P12`, `N1`-`N3`)

| Control | Mundo sintetico | Preexpectativa cuantitativa |
|---|---|---|
| `P1` | Ruido iid $\mathcal{N}(0,1)$ con disparo real | Mediana de $\tau_{\text{rest}}$ crudo $= 1$ paso; $p_{\Delta\tau,\text{iid}} \ge \alpha$ y $p_{\Delta\tau,\text{Onsager}} \ge \alpha$ (compatible con 0) |
| `P2` | $\text{AR}(1)$ con $\phi = 0.8$ y choques | Mediana de $\tau_{\text{rest}}$ en el nulo lineal dentro de $\pm 1.5$ pasos de la semivida con adelanto de primer cruce; exceso sobre el nulo iid $\Delta\tau_{\text{iid}} > 1.5$ pasos con $p_{\Delta\tau,\text{iid}} < \alpha$; $\Delta\tau_{\text{Onsager}}$ compatible con 0 ($p_{\Delta\tau,\text{Onsager}} \ge \alpha$) |
| `P2b` | Mismo $\text{AR}(1)$ $\phi = 0.8$ forzando la via (b) | La mediana de $\tau_{\text{rest}}$ del nulo generativo por inyeccion (b) coincide con la del nulo espontaneo (a) dentro de $\pm 1.0$ paso |
| `P3a` | Proceso TAR (autorregresivo de umbral: lento en extremos $\phi_{\text{ext}} = 0.92$, rapido cerca del origen para igualar la ACF de calma) | $\Delta\tau_{\text{Onsager}} > 0$ (con el signo del umbral) y $p_{\Delta\tau,\text{Onsager}} < \alpha$ |
| `P3b` | $\text{AR}(1)$ ($\phi = 0.75$) y $\text{AR}(2)$ ($\phi_1 = 0.6, \phi_2 = 0.2$) lineales | Curva media $\bar{\eta}(\tau)$ dentro de la banda de $\pm 2.5\,\text{ee}$ (bootstrap de eventos) respecto de la ACF teorica/muestral de la calma para $\tau \in [1, 10]$ |
| `P10` | Paseo aleatorio en `NIVEL` (incrementos iid) | `t_max_efectivo == W_ancla`, `t_max_raiz_unitaria == true`, y $\Delta\tau_{\text{Onsager}}$ compatible con 0 |
| `P12` | Configuracion de nulos | `LILYTH_NULO_BLOQUES` se rechaza en `calibrar`; `FASE_ALEATORIZADA_FFT` en la via (b) cae al modelo generativo y marca `cayo_a_generativo = true` en la ficha |
| `N1` | Control de tamano bajo el nulo ($\ge 200$ semillas independientes) | Tasa empirica de rechazo a nivel $\alpha = 0.05$ de $p_{\Delta\tau,\text{iid}}$ y $p_{\Delta\tau,\text{Onsager}}$ menor o igual a $\alpha + 2.5\sqrt{\alpha(1-\alpha)/200} \approx 0.088$ |
| `N2` | Independencia de realizaciones (1 hilo vs. 4 hilos) | Poblaciones de eventos nulos, $T_{\max}$ y estadisticos identicos **bit a bit** entre `hilos = 1` y `hilos = 4` |
| `N3` | Copia limpia del operador delta por realizacion (Regla 7) | Correr con una copia nueva por realizacion vs. arrastrar el estado de `PurgaDiurna`/escalera de una realizacion a la siguiente produce salidas distintas, demostrando la necesidad de la Regla 7 |
