# Enmienda de E1.9 (horizonte T_max de Lilyth) -- expectativas selladas

Code (Opus 5.5), 2026-11-16. Orden del Comandante: "ve con la primera, modificando la regla de Argon y aclarando el por que" (la primera
decision abierta de `pool_y_fuente_nulos/MEDIDO_POOL_Y_FUENTE.md`). Se sella ANTES de escribir codigo. Base previa en `base_previa/`
(verificada: identica, 70 archivos).

## La regla de Argon y por que no alcanza

E1.9 (Argon): T_max = cuantil 1 - alfa de `tau_rest` (primer toque) en el nulo lineal de la via (a); si la fraccion censurada en el tope
T_MAX pasa de alfa, T_max = T_MAX y se marca censurado; con raiz unitaria, T_max = W_ancla.

Con esa regla queda por construccion ~alfa de los nulos sin restituir dentro de T_max (censurados, codificados en T_max + 1). Un evento real
lento tambien sale censurado y EMPATA con todos ellos. Con una fraccion c de nulos censurados, el p minimo de un evento censurado es del
orden de c (una cola) o 2c (dos colas), aunque el pool tenga miles de nulos: medido en el pool derivado, 41 de 47 tormentas lentas
censuradas empataban y la potencia de un evento quedo en 6.7 %. Ningun tamano de pool lo arregla; lo arregla el horizonte.

## La regla nueva (sin numero nuevo)

Para que un evento censurado pueda resolver el primer escalon del Holm de 4 de Lilyth (alfa/4) en la prueba bilateral, la fraccion de
nulos censurados tiene que ser c <= alfa/(2 * 4) = alfa/8. Por tanto: T_max = cuantil 1 - alfa/8 del mismo nulo lineal, con la guarda de
censura en el tope tambien en alfa/8. 2 = bilateral y 4 = tamano de la familia: los mismos de n* del pool derivado.

## Diseno (regla 6: la de Argon por omision, bit a bit)

- H1. `ConfiguracionLilyth::horizonte` = `HorizonteLilyth::{CUANTIL_ALFA = 0, CUANTIL_RESOLUCION = 1}`, por omision CUANTIL_ALFA (E1.9 de
  Argon). CUANTIL_RESOLUCION: `derivar_t_max` recibe nivel alfa/8 en vez de alfa (cuantil y guarda). Nada mas cambia.
- H2. La ficha registra `nivel_t_max` (el nivel usado). README del operador (base) y tabla de enmiendas del README canonico (fila E1.9b
  con el por que).

## Expectativas

Condiciones:
- ET0 (identidad): con CUANTIL_ALFA, las 6 baterias de Lilyth y el volcado del acople de E7b (IID, RAFAGAS, CAMBIO) identicos bit a bit
  a la base previa.
- ET1 (mecanica): con CUANTIL_RESOLUCION, el T_max derivado de cada serie = un recalculo independiente desde `nulos_lineal` de la ficha
  (cuantil 1 - alfa/8 con interpolacion lineal entre ordenados, censura en T_MAX + guarda en alfa/8, recorte [2, T_MAX], raiz unitaria
  igual que antes); y >= el T_max de CUANTIL_ALFA sobre la misma calibracion, en todas las series del mundo lineal (400).
- ET2 (tamano): mundo lineal de EC6 (AR(1) phi 0.82, tormenta 12, K = 1, 400 replicas), DERIVADO + CUANTIL_RESOLUCION: rechazo Onsager
  con Wilson inferior <= alfa.

Expectativas (se reportan se cumplan o no):
- ET3 (potencia): mismo mundo con la tormenta lenta (z = 4 * 0.97^tau), DERIVADO + CUANTIL_RESOLUCION: la fraccion de tormentas lentas
  censuradas baja de 12.7 % y la potencia por lenta al menos se duplica respecto a 6.7 % (>= 13.4 %).
- ET4 (fuentes): con CUANTIL_RESOLUCION el patron de las fuentes se mantiene en M0 (CALMA rechaza por encima de alfa, VIENTO con Wilson
  inferior <= alfa) y en M2 VIENTO sigue separando las tormentas con z (u medio mayor que sin z).
- ET5 (mutantes, >= 4 de 5): nivel alfa en vez de alfa/8; nivel alfa/4; guarda de censura en alfa con el cuantil en alfa/8; cuantil
  alfa/8 en vez de 1 - alfa/8; CUANTIL_RESOLUCION por omision.
