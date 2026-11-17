# Operador delta (saltos) y temperatura -- `ObservadorSaltos` (ARK-METRO-01)

Codigo: `observador_saltos.hpp` (orquestador), `canal_nivel.hpp`, `canal_temperatura.hpp`, `escalera_diadica.hpp`,
`superestadistica_inversa.hpp`, `sustitutos.hpp`. Pruebas y actas: `JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/metrologia/`
(`EXPECTATIVAS.md`, `MEDIDO.md`, `MEDIDO_REAL.md`, `MEDIDO_DELTA_FRACTAL.md`). Es telemetria: no decide operar.

## Que hace

Ingiere cualquier observable theta(t) a su resolucion nativa: representacion estacionaria (NIVEL -> incrementos), purga diurna causal
opcional, escalera diadica 2^j y, por escala, dos canales:
- **Canal de nivel** (saltos de la media): GLR de Lorden con salto desconocido, h por el ARL0 del nulo, curva t* de Lorden y tau50
  medido por inyeccion. Clasifica en cascada entre escalas: SALTO (persiste hasta la escala mas gruesa), TRANSITORIO (no confirma a
  tiempo), DERIVA (pendiente significativa tras la alarma).
- **Canal de temperatura** (regimenes de la dispersion): superestadistica directa (T* por la curtosis local contra la gaussiana del
  mismo T, momentos de log beta corregidos, q de Tsallis, T1 equilibrio local), separacion de escalas contra los nulos elegidos y GLR
  sobre log beta^. CAMBIO_TEMPERATURA si se confirma en el tiempo.
- f(beta) por inversion de Fredholm (Chambolle-Pock, TV o Tikhonov sobre Df, lambda por discrepancia) en `superestadistica_inversa.hpp`.

Cada evento sale con su ficha formal (`FichaEvento`): escala de origen y de confirmacion, alarma (signo, punto de cambio, amplitud y su
se), plazo, kappa y la ficha completa de la escala (nulo y detector).

## Entradas obligatorias

`ConfiguracionObservador.arl0_pasos_nativos` (una falsa alarma de la familia cada ARL0 pasos nativos; se reparte entre los 2J canales).

## Nulos que acepta (catalogo en `README_NULOS.md`)

| campo | valores | que calibra | por omision |
|---|---|---|---|
| `nulo_nivel` | NINGUNO, PERMUTACION, SIGNOS, BLOQUES (b = N_REG) | la historia con que el canal de nivel mide mu0, s0, tau_int y h | NINGUNO (la historia tal cual) |
| `nulos_separacion` | mascara de PERMUTACION, BLOQUES (b = T* 2^k), SIGNOS | contra que se contrasta la separacion de escalas (Holm entre los elegidos) | ninguno: la separacion queda NO contrastada y `es_superestadistica` = false |
| `replicas_nulo` | B | sustitutos por nulo | 199 |

Guia por variable: `nulo_nivel` = NINGUNO busca cambios respecto del reposo habitual de la variable. SIGNOS (flujos simetricos) o
PERMUTACION (positivos, acotados) fijan el reposo en "sin memoria". Para la separacion de escalas sobre retornos se usaron BLOQUES y
PERMUTACION (SIGNOS conservaria la memoria de la volatilidad, que es justo lo que la separacion mide).

**Que significa "evento" con cada nulo (auditoria del 2026-11-11, decision del Comandante: opciones 1 y 4).**
Acta: `JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/auditoria_delta_choques/AUDITORIA_DELTA_CALIBRACION_CON_CHOQUES.md`.

- **Con NINGUNO el evento es RELATIVO A LA HISTORIA.**
  - El umbral h se busca para que la propia historia de calibracion de una alarma cada ARL0 pasos de la escala. Las alarmas dentro de
    la calibracion quedan fijadas por ese presupuesto, haya los choques que haya: medido 34-38 contra 37.5, y 150 contra 150.
  - "Evento" significa "mas raro que una vez cada ARL0 en la historia". Un tipo de choque mas frecuente que eso es parte de lo
    normal y no se ve: con AR(1) phi 0.88, ARL0 800 y un choque cada 300 pasos se vieron 29 de 202, contra 198 de 202 calibrando
    sobre la misma serie sin choques. Lo que sube es h (1.09 -> 3.07); tau_int no cambia.
  - **Para que un tipo de evento cuente, su frecuencia tiene que ser menor que una alarma cada ARL0 pasos de la escala. El ARL0 se
    elige con eso en mente.**
- **PERMUTACION y SIGNOS son para series SIN memoria.**
  - El canal de nivel mide mu0, s0, tau_int y h sobre la historia SUSTITUIDA. Con PERMUTACION, tau_int sale 1 aunque la serie tenga
    memoria: en un AR(1) phi 0.88, 1.00 contra 15.6.
  - El GLR normaliza entonces la serie real con la memoria equivocada y el ARL0 no se cumple: 1 133 alarmas contra un presupuesto de
    37.5 (28-30 veces), tambien sin choques. El operador se convierte en un detector de memoria.
  - El nulo es una ENTRADA y no se prohibe dentro del operador: depende de la naturaleza de la serie. **Si se sabe que la serie tiene
    memoria, no se elige PERMUTACION ni SIGNOS para `nulo_nivel`.**

## Limites medidos

- Sobre retornos de 1 min el equilibrio local gaussiano no se cumple en ninguna escala (T1 falla; `MEDIDO_REAL.md`).
- El canal de nivel sobre retornos crudos esta casi ciego (el agrupamiento de volatilidad engorda el nulo).
- El ARL0 realizado depende de cuanto se extrapola el h: con ruido gaussiano y poca extrapolacion es conservador en las escalas
  finas; con colas pesadas y mucha extrapolacion es PERMISIVO (x2 a ~10 veces, x14 a ~140 veces, medido el 2026-11-05).
  Ver la seccion "Defecto abierto" al final.
- En el camino de serie (`calibrar`/`paso`) la fase de la purga la sigue entregando el llamador. En el camino por filas sale del
  calendario (abajo).
- Con `nulo_nivel` = NINGUNO las alarmas dentro de la calibracion son las del presupuesto (n / ARL0 por escala), haya los choques que
  haya. Los eventos mas frecuentes que el presupuesto no se ven (2026-11-11, guia de nulos arriba).
- Con `nulo_nivel` = PERMUTACION o SIGNOS sobre una serie con memoria, el ARL0 no se cumple: 28-30 veces el presupuesto en un AR(1)
  phi 0.88 (2026-11-11).
- **Las alarmas no son los choques: quien estime la excitacion entre choques sobre las alarmas tiene que invertir el detector, no
  trasladar sus cifras 1 a 1** (2026-11-11, al construir el nulo con replicas de Lilyth; `PRUEBAS_JEFE/lilyth/etapa_lh/`).
  Medido en un AR(1) phi 0.88 con choques de 9.5 sd, una escala, L = 4, ARL0 800:
  - **Realarmas.** Un choque aislado da en media ~1.76 alarmas, el 89 % en los 10 pasos siguientes. Las alarmas se agrupan aunque los
    choques no esten agrupados.
  - **Pendiente de la ligadura.** Un Hawkes ajustado sobre las alarmas responde a la excitacion verdadera entre choques con pendiente
    0.66-0.88, no 1. Un intervalo calculado sobre las alarmas y trasladado tal cual sale demasiado estrecho y puede no contener su
    propio punto (se vio: 0.342 con [0.181, 0.263]). La cota correcta divide el recorrido por la pendiente medida.
  - **Signo contrario invisible.** La excitacion de signo contrario (un choque que dispara otro del signo opuesto) no mueve el ajuste
    sobre las alarmas: su pendiente no se resuelve a 2 ee y su masa queda en [0, inf). Este detector no la identifica.
  - Con choques de DOS pasos (dos innovaciones grandes seguidas) las alarmas se agrupan sin que haya ninguna excitacion: es una
    desviacion legitima de un nulo iid, no replicas.

## Entrada por filas sobre el observable crudo (2026-11-04, orden del Comandante)

`ObservadorSaltosFilas<J_MAX, L, N_REG, B_MAX, EV_MAX, W_FLUJO>` (`observador_saltos_filas.hpp`): filas -> aduana del motor de pared
auditado -> barra intrinseca (`ark/ventanas/flujo_barras.hpp`, fuente comun) -> `ObservadorSaltos::paso`, que no se toca (su escalera
2^j agrega barras intrinsecas). Por barra: `valido` = sin filas imputadas; `den` = nu de la barra (con RAZON_DE_SUMAS la escalera da
sum(valor) / sum(nu), el observable por unidad de reloj; con las demas no se usa); fase y dia del calendario del activo si
`fuente.fase.activa` (y entonces `bins_fase` = 2 bins_dia), o los del llamador para la fila que cierra la barra.

Fase de calendario (`ark/dataset/fase_calendario.hpp`): fase = tipo x bins_dia + bin horario LOCAL, con tipo 0 laborable (sesion
regular o reducida de lunes a viernes) y 1 no laborable (fin de semana, feriado, cierre extraordinario; en cripto el sabado y el domingo
aunque haya sesion). Zona por calendario: NYSE, SIFMA, CME y FX en Nueva York; SIX en Berlin; cripto en UTC. Conversion UTC -> local
exacta (G1a: 967 680 horas de 1971-2030 en NY y Berlin). Declarado: las sesiones reducidas van con las laborables (unas 3 al ano: un
perfil propio no maduraria con 20 dias por fase).

Controles (`PRUEBAS_JEFE/conexion_filas/`): G5 en PARED con delta_nu 1 = ObservadorSaltos sobre los valores, bit a bit; G6a NU_VOLUMEN
con imputadas y purga con el calendario de BTC (48 fases) = ObservadorSaltos sobre la serie de barras calculada a mano; G6b
RAZON_DE_SUMAS con den = nu; G7b cero reservas; G8b causalidad.

## Modo de deteccion: DISJUNTO o DESLIZANTE (2026-11-05, orden del Comandante: los dos quedan)

### Ventajas y desventajas (medido en sinteticos; la tabla completa con los otros tipos de ventana esta en
`ark/ventanas/README_VENTANAS_Y_DETECTORES.md`)

| | DISJUNTO | DESLIZANTE |
|---|---|---|
| ventaja | h mas bajo al mismo ARL0; primera alarma de temperatura antes (x4: 32.5 contra 35 pasos; x1.5: 257.5 contra 301.5); menos confirmaciones falsas en iid | confirma la temperatura antes (por acumulacion: x4 93 contra 137; x1.5 493 contra 1081.5); evalua todas las escalas en cada paso (no espera a que cierren) |
| desventaja | espera a que cierre la barra o la ventana; confirmar exige una ventana entera nueva | h mayor (miradas correlacionadas): saltos pequenos 5-6 % mas lentos y primera alarma de temperatura despues; 2.3 veces mas confirmaciones falsas en iid |
| coste | referencia | 1.15-1.19 veces |
| usar cuando | importa el primer aviso o el coste de una confirmacion falsa | importa confirmar pronto un cambio de regimen |

`ConfiguracionObservador.modo` (`tipos_metrologia.hpp`, `ModoDeteccion`). Por omision DISJUNTO, que es el operador de siempre: la
copia congelada antes del cambio (`PRUEBAS_JEFE/deslizante/ref/`) y el modo DISJUNTO dan lo mismo bit a bit (D1, cinco mundos).

DESLIZANTE (`detector_deslizante.hpp`):
- **Nivel**: en CADA paso nativo, en todas las escalas. La barra de la escala j es la de las ultimas 2^j barras nativas (escalera
  deslizante: el mismo arbol y el mismo orden de suma que la disjunta, con cualquier agregacion; en los instantes alineados coincide bit a
  bit). El GLR suma l barras separadas 2^j que terminan en el paso actual. Nulo (mu0, s0, tau_int) = el de la escala disjunta; h se
  recalibra sobre la misma construccion deslizante con el mismo reparto del ARL0 (cada canal: ARL0 2J pasos nativos). tau50 por inyeccion
  en una fase al azar (rampa de 2^j pasos; exacta para SUMA y MEDIA, linealizada para las demas, entera para ULTIMO).
- **Temperatura**: la ventana de T* barras avanza UNA barra de su escala (no T*); el GLR suma ventanas separadas T* barras. Nulo de log
  beta = el disjunto; h recalibrado. Declarado: no baja a paso nativo en las escalas gruesas (costaria L T* 2^j de memoria por escala).
  Una ventana constante (s^2 = 0) no se estima y reinicia el detector (el disjunto la salta sin reiniciar).
- Tras una alarma el detector espera D - 1 entradas (D = 2^j en nivel, T* en temperatura): el siguiente elemento se forma solo con datos
  posteriores, como tras el reinicio del disjunto. Una entrada no valida reinicia.
- Siguen sobre barras disjuntas, porque piden ventanas independientes: la superestadistica (T*, momentos, separacion de escalas), el
  regularizador dual (norma, borde, derivada; su GLR disjunto corre pero no se usa) y la ventana de forma de DERIVA. El clasificador es
  el mismo; en DESLIZANTE sus plazos y las fichas de los eventos son las deslizantes (`ficha_nivel(j)`, `ficha_temperatura(j)`).
- Con `nulo_nivel` distinto de NINGUNO, el sustituto se aplica a los pasos NATIVOS por tramos (BLOQUES con b = N_REG 2^j; numerador y
  denominador juntos) y despues se construye la barra deslizante (DIVERGE del disjunto, que sustituye las barras de la escala).
- Reutilizado, sin rehacer: `Acumulador`/`combinar`/`valor` (escalera disjunta), `RingBufferEspejoCxx` (`ark/ventanas`, la ventana de T*
  barras), `detalle_nivel::buscar_h_con` (la busqueda de h de siempre, generalizada sin cambiar sus operaciones). La entrada por filas y
  el sensor de memoria no cambian: el modo viaja en la configuracion (`delta.modo` en `ObservadorMemoria`).
- Coste: cero reservas en el ciclo caliente; memoria de calibracion O(L 2^J) dobles.

Controles (`PRUEBAS_JEFE/deslizante/`, `prueba_deslizante.exe identidad`, 31/31; `mutaciones_deslizante.py`, 16/16 cazados): D1 disjunto
= congelado; D2 escalera deslizante = disjunta alineada y = arbol a mano siempre (6 agregaciones); D3 detector espaciado (D = 1 =
`DetectorGLR` bit a bit; D > 1 = formula directa); D4 temperatura alineada = disjunta; D5 escala 0 deslizante = disjunta; D6
causalidad; D7 cero reservas; D8 entrada por filas; D9 ARL(h) en su objetivo sobre el nulo; D10 fichas de los eventos.

Medido (`PRUEBAS_JEFE/deslizante/MEDIDO_DESLIZANTE.md`, iid y AR(1) 0.5, ARL0 2000): en calma cumple (falsas alarmas de la familia 0.85
del disjunto; nivel en su ARL0 en 16/16 canales; temperatura 10/12, dos con 25 % mas por h extrapolado) y cuesta 1.15-1.19 veces el
disjunto. Pero **no compra latencia**: iguala en saltos grandes, pierde 5-6 % en los pequenos (h mayor por las miradas correlacionadas) y
pierde en temperatura (x4: 35 contra 32.5 pasos; confirmado 537 contra 137, porque la escala con T* = 256 alarma primero y su confirmacion
exige una ventana entera nueva). Decisiones abiertas (a)-(c) en el acta.

Enmienda 1 (2026-11-05, decision del Comandante): en DESLIZANTE la temperatura se confirma POR ACUMULACION. La apertura (G >= h) no
reinicia el GLR; el evento se confirma cuando la evidencia llega a h_c (`FichaDetector.h_confirmacion`, ARL(h_c) = ARL0 / alfa sobre el
mismo nulo) con el mismo signo, y vence en el plazo t* = 2 h_c tau / delta^2 (`plazo_elementos_temperatura`, la misma funcion en el canal
y el observador). La alarma de confirmacion sale marcada (`SalidaTemperatura.confirmacion`). El DISJUNTO no cambia. Medido: confirma antes
que el disjunto en x4, x2, x1.5 y x1.3 (93 contra 137 ... 1493.5 contra 1814.5); la primera alarma sigue sin mejorar; en el mundo iid da
2.3 veces mas confirmaciones falsas (en AR(1) 0.5, la mitad). Controles D9b, D11, D11b; mutantes 21/21.

## Un observable, dos lentes: `ObservadorDeltaFractal` (2026-11-05, especificacion aprobada por el Comandante)

`observador_delta_fractal.hpp`. Lente del crudo (este operador, DISJUNTO o DESLIZANTE) + lente de memoria (`ObservadorMemoria` ELW u
`ObservadorMemoriaDFA`, con su operador delta sobre H^ en cualquiera de los dos modos), sobre las MISMAS barras (misma fuente, misma aduana,
mismo delta_nu). Entrada obligatoria: `arl0_familia_barras` (una falsa alarma de las dos lentes juntas cada ARL0 barras). Reparto sin
numeros a mano: todos los canales (2 J_crudo + 2 J_memoria) con la misma tasa por barra; se calibra cada lente, se leen sus J y se
recalibran los umbrales con `recalibrar_arl0` (misma historia, otro ARL0: identico a calibrar de nuevo). Cada evento de memoria sale con su
etiqueta: MEMORIA_CON_VOLATILIDAD si la temperatura de la escala NATIVA del crudo abrio una alarma dentro de la evidencia del evento
(enmienda 1: la dispersion de las escalas gruesas cambia tambien con la memoria, porque la varianza de una suma de 2^j barras va como
2^{2Hj}); MEMORIA_PURA si no.

Medido en sinteticos (`PRUEBAS_JEFE/delta_fractal/MEDIDO_DELTA_FRACTAL_CONJUNTO.md`): calma en el 15-30 % del limite de falsas alarmas; un
cambio de varianza lo ve el crudo (8/8) y no la memoria; un cambio de memoria lo ve la memoria (8/8; ELW 205 barras, DFA ~380); pureza de
la etiqueta en W2 66-72 % (cumple en DESLIZANTE, falla por poco en DISJUNTO); con los dos cambios, 100 % con volatilidad. Controles F2-F8
24/24; mutantes 10/10 (`mutaciones_delta_fractal.py`).

## Defecto abierto: el h extrapolado y el presupuesto de la familia (2026-11-05)

**Estado (2026-11-05, tarde):** resuelto para el operador de coordinacion por calibracion de familia (`FAMILIA_EPISODIOS`, acta
`PRUEBAS_JEFE/calibracion_h/MEDIDO_FASE2.md`). `buscar_h_con` no se toco; sigue extrapolando y afecta a quien reparte por canal.

`buscar_h_con` lee el ARL(h) en el nulo con >= 30 alarmas y, mas alla, extrapola en linea recta sobre log ARL(h) ([h_b/2, h_b], 8
puntos). El reparto del presupuesto entre muchos canales pide a cada uno un ARL0 que su historia no puede leer. Medido con el operador de
coordinacion (que usa este operador para cada componente y para el crudo de cada activo): con ~140 veces de extrapolacion por canal y
volatilidad estocastica, **este mismo operador sobre retornos da 14 veces sus falsas alarmas nominales**, tambien dentro de muestra.
Afecta a toda familia grande: crudo con muchas escalas, `ObservadorDeltaFractal` (alli midio bien porque la extrapolacion era ~2 veces) y
`ObservadorCoordinacion`. Detalle comun: `../ventanas/README_VENTANAS_Y_DETECTORES.md` seccion 5. Dilema y opciones (calibrar a nivel de
familia, cola con fundamento, nulo mas largo, menos canales): `mensajes/CONTEXTO_AUDITOR_EXTERNO_CALIBRACION_H.md`. Hasta decidir: leer la
ficha (`det.h_extrapolado`, `det.alarmas_nulo`) antes de confiar en el ARL0 declarado.

## Presupuesto de familia: `calibrador_familia.hpp` y `ObservadorSaltosFamilia` (2026-11-05)

Orden del Comandante: la calibracion de familia de la coordinacion para todos los operadores, sin eliminar la arquitectura antigua.
Prerregistro y actas: `JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/calibracion_familia/`.

- **`calibrador_familia.hpp`** (lo generico):
  - `FichaSubfamilia`, `ParametrosFamilia` y `PresupuestoFamilia`;
  - m_min y B;
  - `episodios_coalescencia`;
  - la busqueda de A* (acotamiento x4 y biseccion en log A);
  - la permutacion por bloques y el reparto en hilos.

  El operador de coordinacion se rehizo sobre el y da sus mismas cifras.
- **`ObservadorSaltosFamilia`**: el operador delta sobre un crudo con presupuesto de familia.
  - `BONFERRONI_CANALES` (omision) es el ObservadorSaltos de siempre, bit a bit.
  - `FAMILIA_EPISODIOS`: la familia son sus 2J canales. El nulo es la historia mas B-1 sustitutos por BLOQUES (decision del
    Comandante), con b = N_REG 2^(J-1) y J las escalas de la historia sola; con NIVEL se permutan los incrementos. El operador conserva
    esas J escalas (`ConfiguracionObservador.escalas_max`, nuevo; 0 = sin tope, lo de siempre).
- **Medido en calma** (crudo solo, ARL0 20 000, 2^16 barras):

  | | BONFERRONI | FAMILIA |
  |---|---|---|
  | GAUSS | 2.7e-5 | 7.3e-5 |
  | SV | 8.2e-5 | 5.9e-5 |

  El limite es 1.0e-4. Los dos metodos cumplen; la familia acerca la tasa al nominal 5e-5 desde los dos lados y cuesta x85 en
  calibracion. Con pocos canales no compra seguridad: compra exactitud.

## Estado de cada canal y episodio por canal (2026-11-14, ordenes del Comandante: "se bifurca, no se desmonta")

Actas: `JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/episodio_canal/` (`MEDIDO_PASO1.md`, `MEDIDO_PASO2.md`). Registro de decisiones y
de como se propagan: `JEFE DE ARQUITECTOS DE CODIGO/DECISIONES_DE_DISENO_Y_PROPAGACION.md`. Todo es aditivo: con las entradas por
omision la salida es la de siempre, bit a bit (regresion A/B de 14 baterias en cada paso).

- **Paso 1, `Salida.nivel_abierto` y `Salida.temperatura_abierto`:** el estado real de cada canal al final del paso; `evento_abierto`
  es su o. Motivo: la familia y la coordinacion usaban la junta en su calibracion, y Lilyth habia construido su propio espejo del
  canal de nivel.
- **Paso 2, rama de dos ARL0:** `ConfiguracionObservador.arl0_temperatura_pasos_nativos` (0 = el de nivel, lo de siempre) y
  `recalibrar_arl0(nivel, temperatura)` (el de un argumento vuelve a la rama de un ARL0). Cada canal reparte su ARL0 entre escalas igual
  que hoy (arl0 2J / 2^j).
- **Paso 2, `DefinicionEpisodio` en `ObservadorSaltosFamilia` y `ObservadorCoordinacion`:** CONJUNTO (omision) o POR_CANAL. La salida
  gana en las dos ramas el episodio de NIVEL y el de TEMPERATURA (su union es el episodio de siempre). Con FAMILIA_EPISODIOS y POR_CANAL
  cada canal es su propia subfamilia con presupuesto 2 ARL0 (Bonferroni sobre 2; en la coordinacion, 2 ARL0_sub dentro de cada
  subfamilia de nulo) y su rareza A*, calibradas en alterna (un delta con dos ARL0: conserva la precedencia del CAMBIO_TEMPERATURA sobre el
  evento de nivel) hasta que ninguna cambie, con tope de rondas declarado. Con BONFERRONI_CANALES los umbrales no cambian.
- **Por que:** en reposo, una alarma de temperatura sin confirmar mantiene abierto el episodio CONJUNTO durante todo su plazo de
  confirmacion, m (ceil(t_Lorden) + 1) pasos nativos, cientos a miles (Red E6b, `PRUEBAS_JEFE/red/MEDIDO_E6B.md`): con varios canales el
  episodio casi nunca se cierra y un cambio real no abre episodio nuevo.
