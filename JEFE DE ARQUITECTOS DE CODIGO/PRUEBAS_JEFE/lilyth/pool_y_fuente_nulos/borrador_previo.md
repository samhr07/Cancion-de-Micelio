# Pool derivado y fuente de los nulos de Lilyth -- BORRADOR de expectativas (NO sellado)

Code (Opus 5.5), 2026-11-16. Ordenes del Comandante: "prepara ambas" (la condicion del pool y el control de la fuente de nulos), nivel de
decision = el Holm de 4 de Lilyth, prueba bilateral; y su propuesta de ramificar: dias sin viento para el nulo del operador delta, dias de
viento (los eventos excluidos) para una calibracion exclusiva de Lilyth. Se sella cuando el Comandante apruebe los puntos marcados (D1-D3).

## A. Pool derivado (entrada nueva, regla 6)

A1. `ConfiguracionLilyth::pool_nulo` = `PoolNuloLilyth::{M_MIN = 0, DERIVADO = 1}`; por omision M_MIN (el codigo de hoy, bit a bit).
A2. Con DERIVADO, la condicion es la que ya calcula `resoluble`: prod_k (n_k + 1) >= 2 m / alfa, con m = 4 (el Holm de Lilyth) -> 160 con
    alfa = 0.05. Sin ningun numero nuevo: sale de alfa y del tamano de la familia.
A3. Al evaluar, con K eventos evaluados en la serie, el objetivo por evento es n* = ceil(160^(1/K)) - 1 (K=1: 159; K=2: 12; K=3: 5;
    K=4: 3). Si un evento tiene menos de n* nulos emparejados, se completan por la via que ya existe (mas realizaciones de la via (a); si
    no alcanza, inyeccion de la via (b)) hasta n* o hasta el tope de coste declarado; si no se llega, NO RESOLUBLE con motivo. m_min(alfa)
    sigue como piso de la estabilidad del PIT.
A4. El tope de coste no se elige: es el b_max_nulo que ya existe. Se mide el coste real y se reporta.

## B. Fuente de los nulos de Lilyth (entrada nueva, regla 6)

B1. `ConfiguracionLilyth::fuente_nulo` = `FuenteNuloLilyth::{CALMA_DELTA = 0, VIENTO = 1, COMPLETA = 2}`; por omision CALMA_DELTA (hoy).
    - CALMA_DELTA: los segmentos de calma del operador delta (hoy). Pregunta: "mas lenta que la calma".
    - VIENTO (propuesta del Comandante): las ventanas de los eventos excluidos con su bufer derivado (nota del bufer), DEJANDO FUERA la
      ventana del evento que se contrasta (deja-uno-fuera). Pregunta: "mas lenta que las otras tormentas de su tamano".
    - COMPLETA: toda la serie salvo la ventana del evento contrastado, sin quitar el post-evento. Pregunta: "mas lenta que su memoria
      normal, viento incluido".
B2. Las tres alimentan la misma maquinaria (sustitutos sin memoria, ajuste Onsager, vias (a) y (b)); solo cambia de donde salen los
    segmentos. El operador delta no se toca: sigue calibrado en calma.

## C. Control con verdad conocida (sinteticos, antes de cualquier dato real)

Mundos (semillas nuevas), con el catalogo = las tormentas inyectadas, conocidas por construccion:
- M0 "solo agrupamiento": volatilidad con agrupamiento (tipo HETERO/GARCH), tormentas que son solo agrupamiento, sin mecanismo extra de
  recomposicion. Verdad: no hay recomposicion lenta mas alla de la memoria.
- M1 "mecanismo universal": todas las tormentas llevan ademas una recomposicion lenta extra (semivida mayor que la del agrupamiento).
- M2 "mecanismo en una minoria": solo 1 de cada 4 tormentas lo lleva.

Expectativas (se reportan se cumplan o no; las marcadas "condicion" deciden):
- EC1 (condicion de tamano) M0: COMPLETA y VIENTO rechazan "lenta" con Wilson inferior <= alfa.
- EC2 (expectativa) M0: CALMA_DELTA rechaza por encima de alfa (anticonservadora: la calma no tiene agrupamiento que recomponer).
- EC3 (expectativa) M1: CALMA_DELTA y COMPLETA detectan (>= 80 %); VIENTO NO detecta (el mecanismo comun a todas las tormentas queda
  dentro del nulo: ciega a lo universal).
- EC4 (expectativa) M2: VIENTO detecta las tormentas con mecanismo (las distingue de las otras) con potencia >= la de COMPLETA.
- EC5 (condicion) identidad: con CALMA_DELTA y M_MIN, todo bit a bit igual a la base actual (baterias y volcado de E7b, como XB.0).
- EC6 (condicion) pool: con DERIVADO, toda serie que llega a n* da resoluble verdadero y su p_min <= alfa/4; las que no, NO RESOLUBLE
  con motivo; tamano nominal en M0 con pools de 159 (K = 1).
- EC7 (mutantes): >= 8 de 10, equivalentes declarados.

## D. Decisiones del Comandante antes de sellar

D1. VIENTO cambia la pregunta: un mecanismo que tengan TODAS las tormentas queda dentro del nulo y Lilyth no lo ve (EC3). Propongo
    conservar CALMA_DELTA o COMPLETA como segunda referencia al lado de VIENTO, igual que hoy conviven sin memoria y Onsager, en vez de
    sustituir. ¿Aprobado?
D2. Deja-uno-fuera en VIENTO y COMPLETA (el evento contrastado nunca entra en su propio nulo). ¿Aprobado?
D3. El catalogo del Economista entra solo con datos reales, despues del control sintetico y de la nota del bufer (pasadas 1 y 2). ¿Aprobado?
