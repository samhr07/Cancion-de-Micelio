# Pool derivado y fuente de los nulos de Lilyth -- MEDIDO (ultima fase del operador)

Code (Opus 5.5), 2026-11-16. Expectativas selladas en `EXPECTATIVAS_POOL_Y_FUENTE.md` (sha256 9ece06cd..., `SELLO_EXPECTATIVAS.txt`) antes
de escribir codigo. Todo sintetico; nada es para operar.

## Lo que entro en la base (regla 6: lo de hoy por omision, bit a bit)

- `nulo_lilyth.hpp`: `PoolNuloLilyth::{M_MIN = 0, DERIVADO = 1}`, `objetivo_pool_derivado(K, alfa)` = ceil((8/alfa)^(1/K)) - 1, y en la
  ficha `fuente_externa`, `n_fuente`, `pool_derivado`, `k_candidatos_pool`, `objetivo_pool`.
- `observador_lilyth.hpp`: `ConfiguracionLilyth::pool_nulo` (M_MIN por omision, con su validez); `calibrar(..., nu_cal, fuente_nulo)`
  con mascara opcional (un 0 corta el segmento del nulo; vacia = hoy; largo distinto de `x_cal` = no calibra); en
  `evaluar_nulos_representacion` K = candidatos (RESTITUIDA o CENSURADA_HORIZONTE) y umbral = max(m_min, n*); `evaluar_una_referencia`
  recibe el umbral (elige la via (a) y es el objetivo de la via (b)); el piso para evaluar sigue en m_min; si la via (b) no llega y la
  excedencia ya tenia m_min (y no es menor), se queda la via (a): DERIVADO solo agrega nulos.
- README del operador y README de la base actualizados; manifiesto regenerado (70 archivos), copias de los dos agentes identicas, humo
  de la base compila y enlaza.

## Condiciones

| | resultado |
|---|---|
| **EC5** identidad (mascara vacia, M_MIN) | **CUMPLE.** Las 6 baterias de Lilyth (Etapas 2-6) identicas contra `base_previa/` (`identidad_baterias.log`); el volcado del acople de E7b en IID, RAFAGAS y CAMBIO identico bit a bit: 4289 lineas (`volcado_previa.txt` = `volcado_nueva.txt`). La excepcion esporadica de XB.0 no reaparecio. |
| **EC6a** sellado (resoluble <=> prod (n_k+1) >= 160) | **FALLA, y el error es mio al sellar.** 104 de 346 (tamano) y 122 de 390 (potencia) discrepan. Las dos causas, medidas: empates del real en su extremo (censura en T_max + 1), y extremos asimetricos con la midpoint PIT, donde solo una cola alcanza el max \|S\| y el p minimo es 1/prod en vez de 2/prod (caso real: n = 139, p_min = 1/140). El producto no es ni suficiente ni necesario. Mi segundo atajo, declarado tras sellar (EC6a': "resoluble implica prod >= 160"), tambien es falso (46 y 16 series) y se retira. |
| **EC6a''** (declarado tras sellar) | **PASA.** El `p_min_alcanzable` de la cabecera coincide con una enumeracion completa independiente de los estratos (<= 1e-12) y `resoluble` = (p_min <= alfa/4) en 346 de 346 y 390 de 390 series. Piso m_min: 0 eventos evaluados por debajo. Resolubles con p_min > alfa/4: 0. |
| **EC6m** mecanica (anadida antes de los mutantes) | **PASA.** Objetivo n* = el menor n con (n+1)^K >= 160 (calculado aparte), K = conteo independiente de candidatos, ningun evento bajo el umbral sin intentar la via (b), ninguna via (b) "completada" bajo el umbral, M_MIN por omision. |
| **EC6b** tamano de DERIVADO (mundo lineal, K = 1, 400) | **PASA.** Rechazo Onsager 5/346 = 0.0145 [0.0062, 0.0334]. |
| **EC8** mecanica de la fuente | **PASA.** `n_fuente` = muestras validas en la mascara; mascara toda en 1 = mascara vacia bit a bit (p crudos, p de Holm, S, B, T_max); historia de calibracion y flujo de operacion disjuntos (D2); y, anadido antes de los mutantes, numero de segmentos = tramos de la mascara. |

## Expectativas

**Pool (mundo lineal AR(1) phi 0.82; tormenta de 12 en t = 200; 400 replicas):**

| | M_MIN | DERIVADO |
|---|---|---|
| K (candidatos) | -- | 1: 313, 2: 35, 0: 52 (tamano); 1: 340, 2: 55 (potencia) |
| series resolubles | 11 / 17 | 206 / 242 |
| tamano (rechazo Onsager) | 0/330 | 5/346 = 0.0145 |
| potencia por lenta (z = 4 * 0.97^tau) | 1/383 = 0.0026 | **26/390 = 0.0667 [0.046, 0.096]** |

**EC6c NO se cumple:** espere >= 50 % y sale 6.7 %. DERIVADO convierte un evento unico de imposible (M_MIN: p minimo 2/9) en posible con
tamano nominal, pero la potencia es baja. Diagnostico posterior declarado (`diag_potencia.log`):
- T_max (cuantil 1 - alfa del nulo lineal de la calma) = 6 / 9 / 14 (p10 / p50 / p90).
- De 371 tormentas lentas evaluadas, 47 (12.7 %) salen censuradas; 41 de esas empatan con nulos censurados. Ese empate impide resolver:
  con 3 nulos censurados de 159, p minimo 0.025 > 0.0125.
- Las otras se recomponen dentro de T_max, y el real tiene que superar a sus 159 nulos: 170 de 371 tienen u > 0.9, pocas superan a todos.
  Ademas la via (b) inyecta las innovaciones del perfil del propio evento, que ya llevan parte de z.
- Lo que haria falta, y no esta en esta fase: un T_max con menos censura en el nulo. Con un evento, resolver el primer escalon exige que la
  fraccion de nulos censurados sea menor que alfa/8 (bilateral) o alfa/4 (una cola); con T_max en el cuantil 1 - alfa queda ~alfa. La
  regla que sale sin numero nuevo es T_max = cuantil 1 - alfa/8 del nulo lineal. Es un cambio de E1.9 (decision de Argon y del
  Comandante), no se toca sin orden.

**Fuente (mundos de umbral; delta calibrado en calma; 40 replicas por mundo; rama INTERCAMBIABLE):**

Bufer de la pasada 1 (nota del bufer, Lilyth de la Etapa 1 con horizonte 64 sobre la historia de calibracion como flujo; 100 tormentas
por mundo, 87-91 medidas): posterior 17 (M0), 19 (M1), 16 (M2); previo 0 en los tres (tasa de apertura en calma 0.0023-0.0027 por paso,
1/ARL0 = 0.004). **EC9 SE CUMPLE.**

| mundo | fuente | T_max medio | Onsager rechaza "lenta" | u_pit Onsager medio (sin z / con z) |
|---|---|---|---|---|
| M0 solo agrupamiento | CALMA | 4.7 | 36/37 = 0.973 | 0.9995 / -- |
| | VIENTO | 20.2 | 0/40 [0, 0.088] | 0.650 / -- |
| | COMPLETA | 12.6 | 12/40 = 0.300 | 0.946 / -- |
| M1 mecanismo en todas | CALMA | 6.8 | 35/40 = 0.875 | -- / 0.9997 |
| | VIENTO | 23.5 | 0/40 | -- / 0.673 |
| | COMPLETA | 18.6 | 1/40 = 0.025 | -- / 0.904 |
| M2 mecanismo en 1 de 4 | CALMA | 5.6 | 39/40 = 0.975 | 0.9989 / 1.0000 |
| | VIENTO | 20.2 | 0/40 | **0.583 / 0.704** |
| | COMPLETA | 14.6 | 6/40 = 0.150 | 0.896 / 0.970 |

- **EC1 SE CUMPLE:** en M0 VIENTO tiene tamano (0/40).
- **EC2 SE CUMPLE:** en M0 CALMA rechaza "lenta" en el 97 % sin que haya ningun mecanismo: la calma no tiene excursiones grandes que
  recomponer.
- **EC1b SE CUMPLE:** COMPLETA queda en medio (30 %).
- **EC3 NO SE CUMPLE, y el motivo cambia la lectura.** Espere que COMPLETA detectara el mecanismo universal y sale 2.5 %: la historia
  completa contiene las mismas tormentas con z, asi que su memoria lineal tambien lo absorbe. Solo CALMA "detecta" en M1, pero CALMA
  tambien rechaza en M0 (97 % contra 87.5 %): **no hay fuente que separe un mecanismo comun a todas las tormentas de la memoria normal del
  viento.** Con estos datos son observacionalmente lo mismo; separarlos pediria un modelo estructural.
- **EC4 SE CUMPLE:** en M2 VIENTO calibra (las tormentas sin z quedan en 0.58) y separa las que llevan z (0.70); CALMA da ~1 a todas y no
  separa nada.
- La referencia sin memoria rechaza "lenta" casi siempre con CALMA y COMPLETA (memoria), y menos con VIENTO (15/40 en M0, 0 en M1).

**EC7 (mutantes, `mutaciones_pool_fuente.log`):** el script dice 9 de 10. M1-M3 (mascara ignorada, invertida, sin cortar) los caza EC8;
M4-M9 (n* sin -1, con 1/alfa, via (a) o (b) con m_min, K sin censurados, DERIVADO por omision) los caza EC6m. M10 (piso m_min retirado)
sale VIVO en el informe porque el script excluia la linea entera de EC6a (que falla en la base por la condicion mal sellada); su
subcontrol de piso SI lo caza: 2 y 5 eventos evaluados bajo m_min contra 0 en la base (`build/mutante_10.log`). Con esa lectura, 10 de 10.
M8 es un proxy declarado de "K = evaluados" (K sin los censurados).

## Errores mios, declarados

1. EC6a sellado: "resoluble <=> prod >= 160" ignora los empates y la asimetria de los extremos (ver arriba). Es la misma familia que el
   XB.1f de la bifurcacion: un atajo de cuenta donde la enumeracion exacta ya existia.
2. EC6a' (declarado tras sellar) tambien era falso; se retira.
3. EC3: COMPLETA incluye las tormentas, no podia detectar un mecanismo que todas llevan.
4. EC6c: potencia muy optimista; no considere que con un evento el real tiene que superar a todos sus nulos.
5. El script de mutantes enmascaro M10 por excluir la linea de EC6a entera.

## Lo que queda listo y lo abierto

Listo en la base: el pool derivado y la mascara de fuente, con lo de hoy por omision. La comparacion pedida (D1) queda medida: CALMA =
"mas lenta que la calma" (marea basal, rechaza siempre que haya viento), VIENTO = "que otras tormentas" (con su propia marea temporal de
recuperacion restada, el offset del Comandante), COMPLETA = en medio.

Abierto, a decision del Comandante:
- T_max = cuantil 1 - alfa/8 del nulo lineal para que un evento unico pueda resolver (cambia E1.9).
- El catalogo real del Economista entra solo con datos reales, despues de las dos pasadas del bufer (D3).
