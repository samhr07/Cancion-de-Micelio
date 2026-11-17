# Planos de las Etapas 3 a 8 del Operador de Lilyth -- Code para Argon, 2026-11-06

Continuan el plano de la Etapa 2 (`PLANO_ETAPA2_MODELO_LINEAL_NULO.md`). La especificacion manda: el README canonico
(`PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md`, Enmienda 1 E1.1-E1.22 y su seccion 9). Donde este plano y el README discrepen, manda
el README y me lo dices. Los puntos marcados **[DECIDIDO]** los cerro el Comandante el 2026-11-06 (tabla del final); lo demas lo
decides tu y lo declaras. **La Etapa 7 queda en espera** (decision del Comandante: su bloque exige modificar el dataset).

## 0. Reglas comunes a todas las etapas

1. **Entrega** en `PRUEBAS_JEFE/lilyth/etapaN/` (cabecera, bateria `etapaN_lilyth.cpp`, `mutaciones_etapaN.py`, sus logs). Nada en
   `ark_cxx/`: yo verifico y conecto. Compila con `/std:c++20 /W4` sin avisos y `/MD /fp:precise`.
2. **Se prueba contra el `ObservadorSaltos` real** (`metrologia/observador_saltos.hpp`), nunca contra un simulado. Cada etapa vuelve a
   pasar K1 (identidad del operador delta con y sin Lilyth) de la Etapa 1.
3. **Referencia independiente**: cada cantidad nueva se recalcula con un motor escrito aparte (en la bateria, sin compartir codigo con
   la cabecera) y se compara bit a bit o a la tolerancia que declares con su razon.
4. **Mutantes en TU codigo**, con testigo; un mutante que no compila no esta cazado. Al menos los que lista cada etapa.
5. **Sin numeros a ojo.** Todo tope, ventana o umbral se deriva (de alfa, de una ficha, de un ARL) o es un tope de memoria contado y
   declarado. Si no ves como derivarlo, pregunta antes de escribirlo. (Los parametros de los MUNDOS SINTETICOS de los controles si se
   eligen: son la verdad conocida, y se declaran en la bateria.)
6. **Ciclo caliente sin reservas** (`paso`, `paso_filas`): lo controla un contador de `operator new` como K10. La calibracion si puede
   reservar.
7. **Lilyth no recalibra al operador delta.** Para correr nulos, se copia el operador delta ya calibrado (es copiable), se llama a
   `reiniciar_flujo()` y se le da la serie del nulo. **Cada realizacion parte de una copia NUEVA del original**: la `PurgaDiurna` de la
   copia sigue acumulando y no puede arrastrar estado de una realizacion a otra.
8. Archivos que vas a leer, todos en `ark_cxx/base_auditada/include/ark/`: `lilyth/observador_lilyth.hpp` (la Etapa 1 conectada),
   `metrologia/{observador_saltos,tipos_metrologia,canal_nivel,canal_temperatura,escalera_diadica,estadistica_metrologia,sustitutos,
   calibrador_familia,whittle_exacto}.hpp`, `holm.hpp`, `ventanas/flujo_barras.hpp`, `metrologia/observador_saltos_filas.hpp`.

## Etapa 3.0 -- Continuacion (E1.21) en `observador_lilyth.hpp`

Es un cambio de la cabecera de la Etapa 1, y va primero porque los nulos se miden con la misma regla de aperturas.

- Trabaja sobre una COPIA en `etapa3/observador_lilyth.hpp`.
- Regla: al abrir el espejo con una ranura viva,
  - mismo signo que esa ranura **y** la ranura no ha restituido (en ninguna de sus ramas activas) -> **continuacion**: no corta, no
    abre ranura, `Contadores.continuaciones++`;
  - signo contrario, o la ranura ya restituyo -> **replica**, como hoy.
  Con varias ranuras vivas, la regla se evalua contra la mas reciente; declara por que.
- `SalidaPaso` gana `continuacion` (bool).
- Controles: K5 se rehace (en el mundo de pulsos, las reaperturas dentro del mismo pulso pasan a continuaciones; las replicas que
  quedan son las de signo contrario o tras restituir); un caso a mano de cada una de las cuatro combinaciones (signo igual o contrario,
  restituida o no); el `rama_ref` de la bateria aplica la misma regla y casa bit a bit.
- Mutantes: la continuacion corta; ignora el signo; ignora si ya restituyo; compara con la ranura mas antigua.

## Etapa 3 -- Nulos (a) y (b) conectados al operador delta (`nulo_lilyth.hpp`)

### 3.1 Que produce

Para cada representacion medida (PURGADA, CRUDA), dos poblaciones de eventos nulos, cada uno con la misma ficha que un evento real
(signo, escala, g de la alarma, eta(tau), los dos tau_rest, desenlace):

- **sin memoria**: sustituto `PERMUTACION` o `SIGNOS` (`ark::metrologia::sustituto(segs, nulo, bloque, mu, rng)`, en `sustitutos.hpp`) de la historia de
  calibracion. En NIVEL se sustituyen los incrementos y se reintegra (como en E1.22);
- **lineal**: realizaciones del `ModeloLinealNulo` de la Etapa 2 (uno por representacion, E1.22).

Y con ellas: T_max (E1.9), las dos referencias de cada evento real y las dos pruebas del resumen que tocan aqui (seccion 3.5).

### 3.2 Via (a), espontaneos (README 8.2)

- B realizaciones de la longitud de la historia de calibracion; en cada una corre una copia nueva del operador delta (regla 7) y un
  `ObservadorLilyth` con la configuracion de produccion. Se guardan las ranuras cerradas.
- **Excedencia**: para el evento real k, sus nulos son los eventos nulos de la misma escala de origen y el mismo signo con
  `g >= g_k` (el `g` de `EventoAlarma`). Sin bins.
- Vale (a) para k si hay al menos `familia::m_min(alfa)` nulos (8 con alfa 0.05). B con `familia::realizaciones(m_min, arl0, n, B_max,
  b_pedido)`; si B topa en B_max, la ficha lo marca (`b_topado`).
- Las realizaciones son independientes: `familia::en_paralelo(B, hilos, f)`, un generador por realizacion sembrado por indice
  (nunca un generador compartido, ni reanudado por avance).

### 3.3 Via (b), inyeccion generativa (cuando (a) no llega a m_min)

- Solo el modelo lineal (E1.19). En el nulo sin memoria no se simula: tau_iid = el medido sobre su sustituto, como en (a).
- El perfil del evento real k son los pasos `[t_cambio, t_a]` de la serie de la representacion. Se pasa por `innovaciones(...)` con su
  historia previa, y las innovaciones se inyectan con `simular(...)` en un tiempo al azar (uniforme en la parte de la realizacion que
  deja W_ancla + guarda de historia antes y T_MAX despues).
- **Ventana de una inyeccion** = [t_iny - (W_ancla + guarda), t_iny + l + T_MAX]: la historia del ancla, el perfil y el horizonte de
  medida. Varias inyecciones en la misma realizacion solo si sus ventanas **no se solapan** (separacion derivada de la ventana, no
  elegida). Si dos ventanas se solapan (o un evento espontaneo del nulo cae dentro de la ventana de una inyeccion), la inyeccion
  afectada no entra al nulo y se cuenta en `Contadores.solapes`: su ancla o su trayectoria estarian contaminadas por otro choque.
- **El disparo es el mismo**: la copia del operador delta corre sobre la realizacion; el evento nulo es la primera apertura de Lilyth
  en `[t_iny, t_iny + l]` con el signo de k. Si no dispara, se cuenta como **no detectado** y no entra (la tasa de no detectados va a la
  ficha). Sin esto se mediria una restitucion que el detector no habria visto.
- **Parada [DECIDIDO]:** se repite hasta m_min eventos detectados o un tope de intentos derivado, ceil(m_min / p_det), con p_det
  estimada en las primeras m_min inyecciones y su tope de memoria contado. Si p_det = 0 en esas inyecciones, el evento queda sin nulo
  (b) y la ficha lo dice.
- **Matiz del Comandante:** la parada y la integridad del instrumento se verifican sobre sinteticos con verdad conocida, con replicas
  de Hawkes y con ventanas solapadas (controles S1-S3 de 3.7).
- `nulo_por_inyeccion = true` en la ficha del evento.

### 3.4 T_max (E1.9)

- Los eventos del nulo lineal se miden con horizonte T_MAX (el tope de memoria). T_max = cuantil 1 - alfa de su `tau_rest_primer_toque`
  (`est::cuantil`); si mas de una fraccion alfa queda censurada en T_MAX, el cuantil no existe dentro del tope: T_max = T_MAX y la ficha
  marca `t_max_censurado`.
- Raiz unitaria (radio_max >= 1 en la Etapa 2, o NIVEL con incrementos sin memoria): T_max = W_ancla y se marca.
- Con T_max fijado, `ConfiguracionLilyth.t_max_efectivo` deja de ser una entrada y lo escribe `calibrar`.

### 3.5 Las referencias y sus pruebas

- Por evento k y por nulo: `u_k = PIT de punto medio` de su tau_rest en la distribucion de sus nulos (censurados en T_max como el
  mayor valor; empates a la mitad). Exceso por evento = tau_rest_k - mediana del nulo, con la mediana y su recorrido al lado.
- **Prueba del resumen [DECIDIDO]**: estadistico S = media de (u_k - 1/2) sobre los eventos; su distribucion nula,
  por Monte Carlo, sacando conjuntos del mismo tamano de eventos nulos (cada uno con su PIT contra el resto del nulo); p bilateral.
  Asi salen p_Delta_tau_iid (nulo sin memoria) y p_Delta_tau_Onsager (nulo lineal). Los otros dos p de la familia de E1.16 son de la
  Etapa 4.
- P3b: eta_barra(tau) del real contra la ACF de la calma de la representacion (rho_tau en FLUJO; en NIVEL, la respuesta acumulada
  del AR de incrementos), con banda de 2 ee por bootstrap de eventos.

### 3.6 Orquestador

`ObservadorLilyth::calibrar(const Obs& delta, cfg, x_cal, valido_cal, fase_cal, dia_cal, nu_cal)`: llama al `configurar` de la Etapa 1,
ajusta los modelos de la Etapa 2, corre (a), fija T_max, guarda los nulos. La via (b) se dispara por evento al cerrarse una ranura real
que no tiene m_min nulos en (a): **fuera del ciclo caliente** (`completar_nulos()`, que el llamador invoca cuando quiera; el ciclo
caliente nunca simula).

### 3.7 Controles

| id | control |
|---|---|
| P1 | iid con disparo real: tau_rest crudo ~1 paso y Delta_tau_iid compatible con 0 (p no rechaza a alfa en >= 1 - alfa de las semillas) |
| P2 | AR(1) phi 0.8 con choques: Delta_tau_iid recupera ln 2 / (-ln phi) corregido por el adelanto de primer cruce, medido |
| P2b | el mismo mundo forzando la via (b): la semivida sale dentro del error de la de (a) |
| P3a | TAR con la misma ACF de calma: Delta_tau_Onsager != 0 con el signo del umbral |
| P3b | AR(1) y AR(2): eta_barra dentro de la banda de la ACF |
| P10 | paseo aleatorio en NIVEL: T_max = W_ancla marcado y Delta_tau = 0 |
| P12 | BLOQUES se rechaza; FASE_ALEATORIZADA_FFT en (b) cae al generativo y la ficha lo marca |
| N1 | tamano: en un mundo nulo, las dos pruebas rechazan a tasa <= alfa (+ su error de Monte Carlo) sobre >= 200 semillas |
| N2 | independencia de realizaciones: con 1 hilo y con N hilos los nulos son identicos bit a bit |
| N3 | cada realizacion parte de una copia nueva: sembrar una realizacion con estado previo cambia el resultado (control de la regla 7) |
| S1 | **replicas de Hawkes** (matiz del Comandante): AR(1) de calma conocida; los choques llegan por un Hawkes exponencial con razon de ramificacion n en {0, 0.3, 0.6} (n = 0 es Poisson, el control) y nucleo conocido; cada choque es un perfil inyectado por innovaciones. Se conoce cada choque, su madre y su tiempo. Se exige: (i) continuaciones y replicas (E1.21) = las de la verdad, por evento; (ii) CIF_replica de Aalen-Johansen dentro de su banda contra la incidencia de replicas simulada; (iii) Delta_tau_Onsager compatible con 0 en las tres n (el proceso de calma es lineal: las replicas cortan trayectorias, no cambian la restitucion); (iv) la parada de (b) llega a m_min con p_det estable entre las tres n |
| S2 | **ventanas solapadas** (matiz del Comandante): inyecciones a separacion conocida, por debajo y por encima de la ventana de 3.3. Por debajo: `solapes` = las inyecciones solapadas de la verdad y ninguna entra al nulo; por encima: 0 solapes y todas entran. En un mundo donde la contaminacion SI cambiaria el resultado (el choque anterior todavia sin restituir), el nulo con la guarda reproduce la semivida verdadera y sin la guarda no |
| S3 | **integridad de la parada**: con p_det conocido por construccion (choques cerca del umbral del detector, de amplitud fijada), el numero de intentos hasta m_min sigue su distribucion negativa binomial; con p_det = 0 (choques por debajo de lo detectable) el evento queda sin nulo (b) y se declara, sin bucle infinito |

Mutantes: excedencia `g > g_k` por `g <= g_k`; nulos de cualquier escala o signo; (b) sin exigir que el detector dispare; T_max con
cuantil alfa en vez de 1 - alfa; generador compartido entre hilos; reutilizar la copia del operador delta; PIT sin punto medio;
inyeccion en una serie fija (P8 otra vez, aqui conectado); **sin la guarda de solapamiento**; ventana de solapamiento sin T_MAX;
**p_det estimada con todas las inyecciones** (se mira el resultado para decidir cuando parar); parada sin tope.

Antes de correr S1-S3, escribe sus expectativas (que debe salir en cada n y cada separacion) en el README de la etapa y entregalas ANTES del codigo; yo las
sello con su sha256 antes de que los corras: es la preexpectativa del instrumento.

## Etapa 4 -- Armonico (`armonico_lilyth.hpp`)

- Entrada: eta_barra(tau) por representacion (ya la suma la Etapa 1) y las curvas de los eventos nulos de la Etapa 3.
- **`PRONY_CURVA_MEDIA`** (E1.8): minimos cuadrados de eta_barra(tau) = a1 eta_barra(tau-1) + a2 eta_barra(tau-2) para tau en [2, T_max].
  Oscilador si a1^2 + 4 a2 < 0; periodo 2 pi / arccos(a1 / (2 sqrt(-a2))), **sin recortar**. Estadistico de sobrepaso: min eta_barra en
  [2, T_max].
- **`PRONY_FASE_REESCALADA`** (E1.20): cada eta_k se lee en tau / tau_rest_k; rejilla comun con paso 1 / mediana(tau_rest) (una barra a
  la escala tipica), interpolacion lineal, sin extrapolar; excluye los que no restituyen y los cuenta.
- **`SOBREPASO_POR_EVENTO`**: min eta_k en [tau_rest_k, T_max] por evento, PIT de punto medio contra sus nulos, estadistico como en 3.5.
- **Contrastes** (E1.8): el min de eta_barra contra su distribucion nula, por Monte Carlo con conjuntos de eventos nulos del mismo tamano,
  una vez contra el nulo sin memoria (hay oscilacion?) y otra contra el lineal (hay mas que en calma?). Son los p 3 y 4 de E1.16.
- **Etiqueta `ARMONICO_SOSPECHOSO`** por evento: raices complejas en su propia eta_k (Prony por evento, con al menos 4 puntos) y min por
  debajo del cuantil alfa del nulo. Su tasa en el nulo va al lado. Nunca veta.
- **Holm** (E1.16): `ark::holm(p)` sobre los cuatro p (dos de la Etapa 3, dos de esta), en `ResumenLilyth`. Un p no valido hace fallar
  la familia entera (falla dura de `holm.hpp`), no se omite.
- Controles: P5 (AR(2) subamortiguado: raices complejas y rechazo contra los dos nulos), P5b (AR(1) phi < 0: raiz real negativa, nunca
  oscilador), Prony contra una recurrencia exacta a 1e-12, periodo de un coseno amortiguado conocido, tasa de la etiqueta en un AR(1)
  <= alfa.
- Mutantes: signo del discriminante; recortar el periodo; Prony desde tau = 1; rejilla reescalada incluyendo los no restituidos; Holm
  sobre tres p; omitir un p NaN.

## Etapa 5 -- Fin del calor (`fin_calor_lilyth.hpp`)

Sustituye la entrada `fijar_calor_activo` de la Etapa 1 (que se queda para los controles).

- **Inicio del calor**: una alarma de `temp[j].alarma` (`SalidaTemperatura`) con el signo que aumenta la temperatura. Lee en
  `canal_temperatura.hpp` cual es: comprueba en el codigo si temperatura alta es log beta bajo, no lo supongas.
- **Fin** (E1.7): detector secuencial unilateral hacia mu_0,beta (la ficha de temperatura de la escala: `ficha_temperatura(j).nulo`, con
  mu0, s0 y tau_int de log beta), calibrado sobre el nulo de log beta **desplazado a mu_calor**.
- **mu_calor [DECIDIDO]**: media de log beta en las l* ventanas posteriores a la alarma de temperatura (su propia memoria, ningun
  numero nuevo). Mientras no hay l* ventanas cerradas, el fin del calor no puede declararse (el calor sigue activo).
- **h** con ARL_0,caliente = T_max / alfa (aproximacion geometrica, declarada), con la busqueda de `canal_nivel.hpp`
  (`detalle_nivel::buscar_h_con`). ⚠ Esa busqueda extrapola cuando no hay alarmas legibles, y la extrapolacion ya dio una tasa x5.7 en
  el operador de coordinacion: si `extrapolado`, la ficha lo marca y el fin del calor de esa escala se declara **no calibrado** (no se
  usa).
- Con varias escalas de temperatura calientes a la vez, el calor termina cuando terminan todas; declaralo.
- Controles: P9 (calor apagado en un t conocido: fin en falso antes de t con tasa <= alfa en T_max; el retraso se mide y se reporta);
  consolidado contra primer toque en un mundo con calor (K4b de la Etapa 1 rehecho con el calor detectado).
- Mutantes: hipotesis nula invertida (el defecto F4); mu_calor desde la amplitud de la alarma; ARL con T_max en vez de T_max / alfa;
  usar un h extrapolado; terminar el calor con la primera escala que termina.

## Etapa 6 -- Reloj doble y marea

- **Reloj de actividad** (E1.11, D7): la Etapa 1 ya acumula nu y guarda la pared. Aqui: nu_punto* = (nu en W_ancla) / W_ancla congelado
  al abrir; `tau_act = Delta_nu_rest / nu_punto*` en barras de calma, por rama. NU_TRADES, NU_VOLUMEN o NU_NOCIONAL; **nunca**
  NU_INFO_SHANNON (se rechaza en `configurar`).
- Por filas: `paso_filas` con las barras de `flujo_barras.hpp` / `observador_saltos_filas.hpp`; en barras intrinsecas la pared no
  avanza a paso fijo (P4).
- **Marea** (E1.14): por evento, Delta_tau_marea = tau_rest^cruda - tau_rest^purgada y Delta_RMST_marea. Agregado: mediana de los
  pareados y su IC por bootstrap de eventos. **[DECIDIDO]: descriptivo**, fuera de la familia de Holm (que sigue con los cuatro p de
  E1.16).
- Controles: P4 (actividad congelada: tau grande en pared y normal en actividad; perdida de rigidez: grande en los dos), P11 (AR(1) por
  perfil diurno determinista: CRUDA depende de la fase del choque; PURGADA recupera la semivida; Delta_tau_marea mide el sesgo).
- Mutantes: nu_punto* sin congelar; NU_INFO_SHANNON aceptado; tau_act en pasos de pared; marea con la rama cruzada (cruda contra
  purgada de otro evento).

## Etapa 7 -- Capa 1 (`sensor_eka_lilyth.hpp`): EN ESPERA

**En espera por decision del Comandante (2026-11-06):** el bloque de la Capa 1 exige modificar el dataset, y se hara entonces. En
`ark_cxx/base_auditada/` no hay hoy ningun sensor Eka ni L_Eka de medianas moviles (la unica "Capa 1" de la base es el DFA del bloque
fractal), aunque el §9.3 del README lo da por existente. No se construye nada de esta etapa. Lo que ya queda fijado para entonces:

- Interfaz: el sensor expone su estado de medianas (`std::array<double, N_CAPA1>`) y una funcion `L(D, s, medianas)`; Lilyth ya copia
  las medianas por valor en cada ranura (E1.12, hecho en la Etapa 1).
- Control P6: sin congelar las medianas, tau_100 <= W por construccion (el control demuestra el artefacto); congeladas, no.
- Mutante: medianas leidas del anillo vivo en vez de la copia.

## Etapa 8 -- Paso a la base (la hago yo)

Mover el README canonico a `ark_cxx/base_auditada/include/ark/lilyth/`, conectar las cabeceras verificadas, la prueba de humo, el
manifiesto, `verificar_base.py` y las copias del Arquitecto de Datos y del Fisico de Cabecera. Argon no toca `ark_cxx/`.

## Decisiones del Comandante (2026-11-06)

| etapa | punto | decision |
|---|---|---|
| 3 | estadistico de la prueba del resumen | **aprobado**: media de PIT de punto medio, nulo por Monte Carlo con conjuntos de eventos nulos |
| 3 | parada de la via (b) | **aprobado**: ceil(m_min / p_det) con p_det de las primeras m_min inyecciones; **matiz**: se verifica sobre sinteticos con replicas de Hawkes y con ventanas solapadas (S1-S3) |
| 5 | ventanas de mu_calor | **aprobado**: las l* de la propia alarma de temperatura |
| 6 | el contraste de marea | **aprobado**: descriptivo, fuera de la familia de Holm |
| 7 | que sensor es la Capa 1 | **en espera**: el bloque exige modificar el dataset |

Orden de trabajo para Argon: 3.0, 3, 4, 5, 6. La 7 no se empieza; la 8 la hago yo.
