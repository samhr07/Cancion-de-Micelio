# Enmienda de E1.9 (horizonte T_max de Lilyth) -- MEDIDO

Code (Opus 5.5), 2026-11-16. Expectativas selladas en `EXPECTATIVAS_HORIZONTE.md` (sha256 269b1f53..., `SELLO_EXPECTATIVAS.txt`) antes
del codigo. Orden del Comandante: modificar la regla de Argon y aclarar el por que. Todo sintetico.

## La regla y su por que

E1.9 (Argon): T_max = cuantil 1 - alfa de `tau_rest` (primer toque) del nulo lineal de la via (a), guarda de censura en alfa. Con esa
regla queda ~alfa de los nulos sin restituir dentro de T_max, codificados en T_max + 1. Un evento real lento tambien sale censurado y
empata con todos ellos: su p minimo queda en el orden de esa fraccion (c en una cola, 2c en dos) por grande que sea el pool. En el pool
derivado, 41 de 47 tormentas lentas censuradas empataban.

Enmienda (E1.9b): T_max = cuantil 1 - alfa/8 y guarda en alfa/8. alfa/8 = alfa/(2 * 4): 2 por la prueba bilateral y 4 por la familia
de Holm de Lilyth, los mismos del pool derivado. Ningun numero nuevo.

## Lo que entro en la base (regla 6: E1.9 por omision, bit a bit)

- `nulo_lilyth.hpp`: `HorizonteLilyth::{CUANTIL_ALFA = 0, CUANTIL_RESOLUCION = 1}`, `nivel_horizonte(h, alfa)` (alfa o alfa/8), ficha con
  `nivel_t_max`.
- `observador_lilyth.hpp`: `ConfiguracionLilyth::horizonte` (CUANTIL_ALFA por omision, con su validez); `derivar_t_max` recibe el nivel
  de la ficha en vez de `cfg_.alfa`. Nada mas cambia.
- README de la base (bullet nuevo) y README canonico (fila E1.9b en la tabla de enmiendas, con el por que). Manifiesto regenerado
  (70), copias de los dos agentes identicas, humo compila y enlaza.

## Condiciones

| | resultado |
|---|---|
| **ET0** identidad (CUANTIL_ALFA) | **CUMPLE.** 6 baterias de Lilyth identicas (`identidad_baterias.log`); volcado del acople de E7b identico bit a bit en IID, RAFAGAS y CAMBIO (4289 lineas). |
| **ET1** mecanica | **PASA** en los 4 casos x 400 series: el T_max coincide con un recalculo independiente desde `nulos_lineal` (cuantil 1 - alfa/8 con interpolacion lineal, censura en 64 y guarda en alfa/8, recorte [2, 64]); nunca es menor que el de E1.9 sobre la misma calibracion; ninguna raiz unitaria; la omision sigue siendo E1.9 (anadido antes de los mutantes). |
| **ET2** tamano (DERIVADO + CUANTIL_RESOLUCION, mundo lineal, K = 1) | **PASA.** Rechazo Onsager 4/346 = 0.0116 [0.0045, 0.0293]. |

## Expectativas

- **ET3 NO SE CUMPLE.** Tormentas lentas censuradas: 42 de 370 (11.4 %) contra 47 de 371 con E1.9 (12.7 %); baja, pero poco. Potencia por
  lenta: 28/390 = 7.2 % [5.0, 10.2] contra 6.7 % con E1.9 (espere >= 13.4 %). T_max p10/p50/p90 = 6/10/18 contra 6/9/14.
- **ET4 SE CUMPLE.** Con CUANTIL_RESOLUCION el patron de las fuentes se mantiene: M0 CALMA 37/38 rechaza "lenta", VIENTO 0/40, COMPLETA
  10/40; M2 VIENTO separa las tormentas con z (u medio 0.684 contra 0.586 sin z). Los T_max suben (VIENTO 32-41, COMPLETA 20-34, CALMA
  6-11).
- **ET5: 4 de 5 cazados.** H1 (nivel alfa), H2 (alfa/4), H4 (cuantil alfa en vez de 1 - alfa) y H5 (CUANTIL_RESOLUCION por omision) los
  caza ET1. H3 (guarda de censura en el alfa del proyecto) es equivalente en estos mundos: ningun nulo lineal queda censurado en el tope
  de 64 (fraccion mediana 0.0000), asi que la guarda nunca actua.

## Por que la potencia casi no se mueve (diagnostico posterior declarado, `diag_horizonte_*.log`)

La via (a) deja **11 nulos lineales por serie** (p10/p50/p90 = 5/11/17). Con 11 puntos, el cuantil 1 - alfa/8 (posicion 0.994 * 10) y
el 1 - alfa (0.95 * 10) caen los dos entre los dos ultimos valores ordenados: casi el maximo de la muestra en los dos casos. La regla es
correcta pero no se puede estimar con esa poblacion. Estimar el cuantil 1 - alfa/8 pide al menos 1/(alfa/8) = 8/alfa = **160** nulos
lineales.

Esa poblacion la fija B, el numero de realizaciones de la via (a): `familia::realizaciones` apunta hoy a m_min eventos
(B = ceil(m_min * ARL0 / n), topado en `b_max_nulo` = 16). Con el horizonte nuevo, el objetivo deberia salir del mismo nivel:
B = ceil((8/alfa) * ARL0 / n). En este mundo (ARL0 2000, n 3000) eso es B = 107 contra 6 hoy, por encima del tope por omision.

## Lo que queda listo y lo abierto

Listo: la regla enmendada como rama seleccionable, con E1.9 por omision, verificada contra un recalculo independiente y sin cambiar nada
cuando no se pide.

Abierto, a decision del Comandante (no se toca sin orden):
- Con CUANTIL_RESOLUCION, derivar B del mismo nivel (objetivo 8/alfa nulos lineales en vez de m_min) y decidir que hacer con el tope
  `b_max_nulo` (subirlo, o derivarlo tambien). Sin eso, la enmienda no cambia casi nada en la practica.
