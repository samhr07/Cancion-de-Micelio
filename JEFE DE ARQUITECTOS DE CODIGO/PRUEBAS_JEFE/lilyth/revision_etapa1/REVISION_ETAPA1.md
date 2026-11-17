# Revision de la Etapa 1 del Operador de Lilyth (entrega de Argon) y conexion a la base -- Code, 2026-11-05

Entrega revisada: `PRUEBAS_JEFE/operador de lilyth/observador_lilyth.cpp` (cabecera), `prueba_lilyth_etapa1.cpp` y
`readme_operador_lilyth.md`. Copia tal cual en `argon_tal_cual/` (los originales de Argon no se tocan). Version corregida y conectada:
`ark_cxx/base_auditada/include/ark/lilyth/observador_lilyth.hpp`, con lo que cambia marcado "Code:".

## 1. La entrega tal cual

- Compila (con la ruta corta: el include relativo superaba 260 caracteres) y **su bateria pasa 21/21**.
- Lo que esta bien y se conserva sin cambios: los tipos (`RamaTrayectoriaLilyth`, `RanuraLilyth`, los desenlaces), Aalen-Johansen y
  su referencia independiente, el desglose por signo, eta con signo, primer toque y consolidado, el reloj normalizado, la copia de las
  medianas de la Capa 1, el rechazo de NINGUNO y BLOQUES.
- **Pero su bateria prueba contra un operador delta simulado** (`OperadorDeltaPrueba`, con una `Salida` que el real no tiene), asi que no
  ve como se comporta con el operador real. Y dos controles no pueden fallar: la identidad compara dos simulados que Lilyth nunca toca, y
  el "mutante de un paso al futuro" altera la ENTRADA del arnes, no el codigo de Lilyth.

## 2. Con el operador delta real (`argon_tal_cual/sonda_real.cpp`)

FLUJO iid N(0,1), ARL0 5 000, pulsos de +3 sd durante 200 pasos en t = 1000, 2000 y 3000:

| evento | alarma (t_a real) | Lilyth de Argon abre en | por que |
|---|---|---|---|
| 1 | 1000 | **1120** | abre con el evento ya CLASIFICADO (t_cierre), no con la alarma |
| 2 | 2001 | **2241** | idem, 240 pasos tarde |
| 3 | 3002 | **3122** | idem |

Abre 120-240 pasos despues del choque, en mitad del pulso, y mide desde ahi: los `tau_rest` que da (14, 1, 47) no tienen relacion con
la vuelta real (el pulso termina a los 200 pasos). Es el punto 9.1 del README (apertura = alarma) sin cumplir.

## 3. Defectos y su correccion

| id | defecto (entrega de Argon) | correccion de Code |
|---|---|---|
| R1 | apertura con el evento clasificado (`Salida.eventos`, que el operador delta emite al cerrar el evento) | espejo del estado "evento de nivel abierto" con la misma regla de `ObservadorSaltos::paso`: abre la primera alarma de nivel (de fina a gruesa) con el evento cerrado; se cierra cuando el operador delta emite el evento de canal 0, y se comprueba que casa con la apertura capturada |
| R2 | Lilyth POSEE un operador delta (`Obs obs_`) y su `paso(y, valido)` le pasa `valido` como `den` (firma real: `paso(x, den, valido, fase, dia)`) | no posee nada: recibe la `Salida` por referencia constante; se retira esa sobrecarga |
| R3 | reloj propio (`++t_actual_`) | el reloj es `Salida.t`; un paso saltado se cuenta (`pasos_desalineados`) y no se mide |
| R4 | `W_ancla = det.L * m`: `det.L` es la memoria del GLR, no N_REG (los dos valen 120 por omision y el error no se veia); 120 escrito de respaldo | `N_REG` del tipo del operador delta (`RasgosDelta`), sin respaldo |
| R5 | la ventana del ancla podia salirse del anillo (HIST_MAX) y la mediana se calculaba con lo que quedara, en silencio | la ventana entera tiene que estar en el anillo: si no, `NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA` |
| R6 | el factor de la purga lo entregaba el llamador, sin relacion con la purga del operador delta | espejo de la `PurgaDiurna` del operador delta (misma clase, configuracion y entradas), alimentado tambien con la historia de calibracion, porque la del operador delta la acumula y la conserva en el ciclo caliente (`reiniciar_flujo` no la toca) |
| R7 | NIVEL: tras un incremento purgado invalido la rama purgada quedaba muerta para siempre | el corte abre un segmento nuevo de la reintegracion; las ramas vivas se censuran y un ancla que cruza un corte da `NO_MEDIDA_ANCLA_INSUFICIENTE` |
| R8 | una replica no cerraba las ramas ya consolidadas: seguian sumando a la curva media eta_barra la trayectoria del choque nuevo | la replica las cierra como RESTITUIDA en su tau |
| R9 | el calor se deducia de los eventos de temperatura clasificados (y un evento de NIVEL cerrado por la temperatura, `estado = CAMBIO_TEMPERATURA` con canal 0, se leia como de temperatura) | Etapa 1: el calor es una entrada (`fijar_calor_activo`); el fin del calor calibrado es la Etapa 5 |
| R10 | `t_max_efectivo = 128` por omision | 0 = sin fijar, y sin fijar no se configura (lo deriva la Etapa 3) |
| R11 | `SalidaPaso` llevaba RANURAS ranuras por valor y se ponia a cero en cada paso (~70 KB por paso con T_MAX 512) | salida pequena; las cerradas se copian a `cerrada(i)` solo al cerrar |
| R12 | signo de dY inconsistente contado como "ancla insuficiente" | motivo propio, `NO_MEDIDA_SIGNO_INCONSISTENTE`; cada motivo se cuenta |
| R13 | la integracion se configuraba en Lilyth aparte de la del operador delta | sale de la `ConfiguracionObservador` del operador delta (una sola fuente) |

Tambien: include `ark/metrologia/...` (antes por nombre suelto), codigos de nulo de la base (antes una copia), y el camino por filas
(`paso_filas`, con nu y pared de la barra intrinseca).

## 4. Controles contra el operador delta real (`etapa1_lilyth.cpp`): 17/17 (`etapa1_lilyth.log`)

Un motor de referencia escrito aparte reconstruye cada ranura desde la serie grabada, la ficha de la escala y la lista de aperturas, y se
compara bit a bit. Operador delta `ObservadorSaltos<24, 64, 120>` (memoria del GLR 64 != N_REG 120, para separar R4).

| id | resultado |
|---|---|
| K1 | identidad del operador delta con y sin Lilyth, 6 000 pasos |
| K2 | 23 eventos de canal 0 = 23 aperturas, 0 inconsistentes |
| K2b | con eventos de temperatura en medio (historia 2^20: con 2^14 el canal de temperatura no calibra, medido; varianza x9): 322 eventos de canal 0 y 31 de canal 1, 0 inconsistentes |
| K3 | el primer pulso abre en la alarma (t = 1000-1010), no en su clasificacion |
| K4 | FLUJO: las 23 ranuras = referencia (ancla, dY, eta, primer toque, consolidado, desenlace, tau_fin) |
| K4b | el calor separa primer toque y consolidado |
| K5, K5b | replicas (16 de 23 aperturas) cortan; eta_barra = referencia a 1e-12 |
| K6, K6b | NIVEL con purga (4 fases, varianza diurna x4): 11 ranuras = referencia, espejo de la purga y segmentos; 3 anclas cruzan un corte |
| K7 | HIST_MAX 128: ancla fuera de memoria = referencia |
| K8 | Aalen-Johansen de la corrida = la referencia independiente de Argon (22 eventos) |
| K9 | causalidad por truncado (estado en t < 3000 igual con la serie cortada) |
| K10 | cero reservas en `paso` |
| K11 | configuracion y reloj desalineado |

Sobre el "mutante de un paso al futuro" pedido para la causalidad: en streaming Lilyth no tiene acceso al futuro (recibe una barra por
paso), asi que esa fuga no es expresable en su codigo; la unica via es la configuracion (historia de calibracion que incluya la prueba),
y esa es responsabilidad del llamador. Los mutantes de abajo cubren lo que si puede fallar.

## 5. Mutantes (`mutaciones_etapa1.py`)

**18/18 cazados, testigo pasa** (`mutaciones_etapa1.log`). Un mutante que no compila no cuenta como cazado.

| id | mutante | lo caza |
|---|---|---|
| M1 | abre con cualquier alarma (sin espejo) | K2, K2b |
| M2, M3 | guarda 0; ancla hasta t_cambio | K4, K5b, K6, K7 |
| M4 | W_ancla con det.L (el defecto R4) | K4, K5b, K6, K7 |
| M5 | dY a la mitad | K4, K5b, K6, K7 |
| M6 | sin la guarda del signo | K4, K5b, K7 |
| M7 | el consolidado ignora el calor | K4, K4b, K7 |
| M8 | la replica no cierra las ramas consolidadas (el defecto R8) | K4, K5b, K6, K7 |
| M9 | espejo de la purga sin la historia de calibracion (R6) | K6 |
| M10, M11 | el corte no abre segmento (R7); origen del nivel purgado en x(t) | K6 |
| M12 | Aalen-Johansen sin S(tau - 1) | K8 |
| M13, M15 | tau_fin de la replica y horizonte uno de mas | K4, K6, K7 |
| M14 | sin la guarda del reloj (R3) | K11 |
| M16 | sin la guarda de memoria del ancla (R5) | K7 |
| M17 | PURGADA y CRUDA invertidas | K4, K6, K7 |
| M18 | el espejo se cierra con eventos de temperatura | K2b (la primera pasada no lo cazaba: ningun mundo tenia eventos de temperatura; de ahi K2b) |

## 6. Hallazgo para el Comandante: replicas del mismo choque

Con un desplazamiento sostenido el GLR vuelve a alarmar respecto de mu0 cada pocos pasos. El operador delta cierra su evento (SALTO o
TRANSITORIO) mucho antes de que el observable vuelva y reabre otro sobre el MISMO choque. En el mundo de pulsos, de 14 aperturas 7 son
asi (en la corrida con la varianza x9 añadida despues, 16 de 23): el pulso de 1000-1199 da aperturas en 1000, 1065 y 1195; la de 1195 se ancla dentro del pulso (Y* = 2.79) y "restituye" en 2 pasos
cuando el pulso termina. Hoy Lilyth las trata como replicas (cortan y abren ranura). Para la Etapa 3 hace falta decidir si una apertura
del mismo signo con una ranura viva sin restituir es una "continuacion" del choque (no corta ni abre) o una replica.

## 7. Conexion con la base

- `ark_cxx/base_auditada/include/ark/lilyth/observador_lilyth.hpp` y su `README_OPERADOR_LILYTH.md`.
- Prueba de humo: incluye la cabecera e instancia `ObservadorLilyth<ObservadorSaltos<>>`; compila con /W4 sin avisos.
- Fila nueva en `ark_cxx/base_auditada/README.md`; manifiesto regenerado (65 archivos) y las copias del Arquitecto de Datos y del
  Fisico de Cabecera (`nucleo_cxx_auditado`) verificadas IDENTICAS. La principal se pone al dia con
  `desplieguectualizar_base_auditada.ps1 -Aplicar` (lo corre el Comandante; copia los arboles enteros, `lilyth/` incluido).
- El README canonico del operador sigue en `PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md`. ⚠ La copia que Argon dejo en
  `PRUEBAS_JEFE/operador de lilyth/readme_operador_lilyth.md` es una variante (recoge las correcciones de Code dentro de E1.2, E1.13 y
  E1.14 en vez de E1.16-E1.20 y un seccion 9 mas corta): no es la canonica.
- La bateria de Argon (`prueba_lilyth_etapa1.cpp`) no compila contra la version conectada (su simulado no es un `ObservadorSaltos`); la
  sustituye `etapa1_lilyth.cpp`, que conserva su referencia de Aalen-Johansen.
