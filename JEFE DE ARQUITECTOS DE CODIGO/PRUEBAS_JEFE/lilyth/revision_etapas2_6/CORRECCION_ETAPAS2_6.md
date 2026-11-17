# Correccion de las Etapas 2-6 del Operador de Lilyth -- Code, 2026-11-07

Orden del Comandante (2026-11-07): "Procedelo a corregir tu" y autorizacion de la opcion (b) para el nulo sin memoria. Defectos de
partida en `REVISION_ETAPAS2_6.md`. Las cabeceras corregidas estan en `inc/ark/lilyth/` de esta carpeta (la entrega de Argon, intacta, en
`argon_tal_cual/`); lo que cambia lleva "Code (revision ...)". Cada correccion tiene su control y su mutante.

## 1. Que se corrigio

| id | defecto | correccion | control | mutante |
|---|---|---|---|---|
| D1 | relleno del nulo sin memoria: eventos fabricados sin detector cuando (a) daba < 16 | quitado. Si la excedencia (a) no llega a m_min: **via (b) sin memoria** (aprobada): el perfil real, en el espacio de la representacion, se escribe en un sustituto sin memoria de la calma y corre el mismo detector; misma parada y misma guarda. Si (b) tampoco llega a m_min, el evento queda **sin nulo** y la ficha lo dice | D1, D1b; P2 y P5 pasan sin relleno | C01, C02, C04 |
| D1b | `construir_muestra_nula_emparejada` completaba la muestra con g < g_k, con eventos de cualquier escala o signo y, si no habia nada, con valores 1.0 | devuelve solo la excedencia; el pool del Monte Carlo usa la misma regla en la funcion libre `pool_pit_loo` (PIT de cada nulo contra su excedencia sin el, solo con m_min) | D1 (excedencia exacta, ningun evento con < m_min); D1c (el pool de `pool_pit_loo` = referencia escrita a mano) | C03, C10 |
| D1c | el armonico filtraba los nulos con \|dY\| >= 0.5 mediana (numero escrito) y los "enriquecia" con eventos fabricados hasta 64 | usa los nulos EMPAREJADOS de la Etapa 3 de cada evento real; Monte Carlo emparejado (un nulo al azar de cada conjunto por sorteo); el min de eta_barra real se calcula sobre los mismos eventos que tienen nulo | P5 (Etapa 4) | C11, C12 |
| D2 | la guarda tomaba las replicas del propio choque como eventos ajenos (646 descartes en la bateria) | regla final: **solo contamina una apertura en la ventana previa [ini, t_iny)**. Una apertura posterior al perfil no contamina, este viva o cerrada la ranura propia (las realarmas llegaban tambien con la ranura ya cerrada; la primera version solo cubria la viva). El destino de cada inyeccion es la funcion explicita `destino_inyeccion` (ADMITIDA, CORTADA, NO_MEDIDA); las CORTADAS quedan fuera, simetrico con los eventos reales cortados. Contadores `con_ajena_previa` y `con_reapertura_posterior` | S2 (aisladas: solapes == previas y posteriores > 0); S2c (prueba unitaria de `destino_inyeccion`: el corte por una replica de signo contrario no se puede forzar en un mundo, porque el espejo no abre con el evento del operador delta abierto) | C05, C06; M09, M10 de Argon ahora se cazan |
| D5 | la parte fisica de S2 movia el ancla a mano | S2 corre el instrumento con y sin la guarda (`aplicar_guarda`); mundo con ARL0 = 1e6 para que no haya alarmas espontaneas (con 2000 sobre 1200 pasos la probabilidad era ~0.45 y la guarda descartaba con razon) | S2: sin guarda la segunda inyeccion contaminada da tau 19.5 contra 27.0 aislada | M09, M10 |
| D6 | numeros a ojo | B = familia::realizaciones con el ARL0 de la familia (sin el factor 2 ni el 1000 de respaldo); alarmas del fin del calor = m_min(alfa) (no 8); B de Monte Carlo y bootstrap = `b_mc_resumen`, por omision 1999 (el prerregistro; el codigo usaba 499) y con la condicion 1/(B+1) <= alfa/4; banda de P3b z_{1 - alfa/(2H)} (no 2.5) sin suelo escrito, hasta el K de la calma (no 10); NIVEL parte de x0 (no 100); P3a con el signo y p < alfa (sin 2.0 ni 0.15). Se retiro la sobrecarga de `calibrar` que deducia la configuracion del operador delta (NIVEL por la longitud de la historia) | N3 (B = 5), P3a | C09 |
| D7 | los nulos entraban al operador delta con fase 0 y un dia por paso; la rama purgada se purgaba dos veces | `EspejoPurgaNulo`: copia invertible de `PurgaDiurna` (misma aritmetica) con `factor(fase, dia)`. La serie nula se construye en el espacio de la representacion y `llevar_a_crudo` la lleva al crudo (y = p / factor en la rama purgada) con las fases de la calibracion y los dias desplazados | D7: espejo = PurgaDiurna bit a bit (6 600 pasos), ida y vuelta 1.9e-16, la rama PURGADA calibra sus nulos | C07, C08 |
| -- | historico de ranuras cerradas: lo que no cabia se perdia en silencio | `Contadores.historico_desbordado` | -- | -- |
| D4 | K4b sobrescribia su resultado con la ultima ranura del pulso; K5 esperaba algo que el mundo no da | K4b sobre la primera ranura del pulso; K5 en un mundo de pulsos de +8 sd (la ranura no restituye durante el pulso): una ranura por pulso, 8 continuaciones, = referencia. El mundo de +3 sd queda informativo | K4b, K5 | -- |

Ademas S1 separa lo que es del instrumento (E1.21 = verdad, CIF de replicas, estabilidad de p_det) de las **preexpectativas de Argon**
(Delta_tau_Onsager ~ 0 en las tres n; p_det >= 0.70), que se informan aparte y no hacen fallar la bateria. La batería 3 gana D1, D1b y D7.

## 2. Resultados (`etapaN.log`)

| bateria | resultado | cifras |
|---|---|---|
| 2 | PASA | -- |
| 3.0 | PASA | K5: una ranura por pulso, 8 continuaciones |
| 3 | PASA | P2: Delta_tau_iid = 4 (p 0.0005), Onsager p 0.13; N1: rechazo 0.005 y 0.005 (cota 0.088); N3: B = 5; S2 sin guarda 19.5 contra 27.0; D1: los 19 eventos usan (b) sin memoria, 6 quedan sin nulo lineal (declarado) |
| 4 | PASA | P5 sin relleno; "P5 emparejado" recalcula min, p y mediana del nulo bit a bit (la mediana hace falta: los dos p saturan en 1/(B+1)) |
| 5 | PASA | P9: fin en falso 0.030 (alfa 0.05) |
| 6 | PASA | P4, P11 |

Mutantes, con la bateria completa y los testigos en verde:
- de Argon: **68/68** cazados (`mutantes_argon_corregido.log`);
- de Code (`mutaciones_revision.py`): **12/12** (`mutaciones_revision.log`). En la primera pasada C12 no se cazaba (11/12); lo caza
  "P5 emparejado" tras anadirle la mediana del nulo, y C11 sigue cazado.

Bateria de la Etapa 1 contra este orquestador (`etapa1_contra_nuevo.log`): fallan **K2, K3, K5 y K2b**, los cuatro porque cuentan como
replicas lo que con E1.21 son continuaciones (la apertura con la ranura viva del mismo signo ya no corta ni abre). Es el cambio de
semantica decidido, no una regresion; la bateria 3.0 cubre la semantica nueva (K4b, K5).

## 2b. Conexion a la base (2026-11-07)

Las seis cabeceras de `inc/ark/lilyth/` estan en `ark_cxx/base_auditada/include/ark/lilyth/` del worktree; el arbol `inc/` entero es
identico a la base (`diff -rq` vacio), asi que las seis baterias pasan desde ella. Humo: la base compila en una unidad con /W4 sin
avisos, con las seis cabeceras incluidas e instanciados `ObservadorLilyth<ObservadorSaltos<>>` y `EspejoPurgaNulo`. README de la base
reescrito para las Etapas 1-6; manifiesto regenerado (70 archivos); las copias de `ARQUITECTO DE DATOS/` y `FÍSICO DE CABECERA/`
verifican identicas. La copia principal NO se toco: la pone al dia `despliegue\actualizar_base_auditada.ps1 -Aplicar`.

## 3. Hallazgos para el Comandante

1. **S1(iii) no se cumple con n = 0.6**: el contraste de Onsager rechaza (p = 0.005 en la corrida final; 0.011 en una intermedia). Con replicas de Hawkes del mismo signo, E1.21 hace
   de las hijas continuaciones del mismo choque, y la restitucion es mas lenta que la respuesta lineal al perfil observado (las hijas llegan
   despues de t_a). No es un defecto: Onsager mide justo esa desviacion. Pero significa que **Delta_tau_Onsager mezcla dos cosas** (no
   linealidad y forzamiento posterior al disparo) cuando hay cascadas; para separarlas haria falta un nulo con replicas.
2. **En los mundos AR, el nulo sin memoria nunca llega a m_min por (a)**: el detector, calibrado sobre la serie con memoria, casi no alarma
   sobre el sustituto permutado. Toda la referencia sin memoria sale de la via (b) sin memoria. Es la opcion aprobada, y conviene saberlo.
3. Con dY de una sola observacion, la "restitucion" de primer toque llega por regresion a la media del disparo (K5 informativo).
