# Bifurcacion del contraste del resumen de Lilyth: POOL_LOO (por omision) e INTERCAMBIABLE -- expectativas selladas

Code (Opus 5.5), 2026-11-15. Orden del Comandante: "bifurca el resumen de Lilyth, y la recomposicion se lee con ambas" (las dos
referencias, sin memoria y Onsager, cada una con su lectura y su Holm de red), con los puntos 1-6 aprobados. Se sella ANTES de escribir
codigo. Base previa copiada en `base_previa/` (verificada contra el manifiesto: identica, 70 archivos).

## Diseno (puntos 1-6, aprobados)

1. Entrada nueva `ConfiguracionLilyth::contraste_resumen` = `ContrasteResumenLilyth::{POOL_LOO = 0, INTERCAMBIABLE = 1}`, por omision
   POOL_LOO (el codigo de hoy, linea por linea). INTERCAMBIABLE: por evento real evaluado, estrato = {el real} U {sus nulos emparejados};
   u de cada miembro = PIT de punto medio contra el resto del estrato; S = media(u_real - 0.5); nulo por Monte Carlo sorteando un miembro
   uniforme por estrato; mismo B (`b_mc_resumen`) y misma semilla; p = (1 + extremos)/(1 + B).
2. Cubre los dos usos de `evaluar_contraste_resumen`: Delta tau (iid y Onsager) y SOBREPASO_POR_EVENTO (estrato {s_k} U su conjunto).
3. Resolucion declarada, no forzada: `ContrasteResumenNulo` gana `intercambiable`, `pool_min_estrato`, `p_min_alcanzable` (p EXACTO de la
   configuracion mas extrema, con el piso 1/(B+1)) y `resoluble` (p_min_alcanzable <= alfa/4, el primer escalon del Holm de 4). m_min no
   se toca.
4. Con INTERCAMBIABLE la ficha lleva tambien la rama vieja al lado (`contraste_iid_pool_loo`, `contraste_onsager_pool_loo` y los dos del
   sobrepaso por evento), para el A/B en la misma corrida. El Holm de 4 usa la rama elegida.
5. Pruebas: identidad bit a bit de POOL_LOO contra la base previa, controles con verdad conocida, mutantes; manifiesto regenerado y copias
   de los agentes recopiadas al final.
6. E7b otra vez, en mundos nuevos (otras semillas), con las dos ramas y leyendo la recomposicion con AMBAS referencias.

Supuesto declarado: los estratos de eventos distintos pueden compartir nulos (via (a): las excedencias espontaneas son un pool comun); el
Monte Carlo los sortea como independientes, igual que la rama vieja.

## Expectativas

**XB.0 (identidad, condicion).** Con POOL_LOO, todo campo que ya existia es identico bit a bit entre la base previa y la nueva:
(i) las baterias de Lilyth de las Etapas 2-6 (`PRUEBAS_JEFE/lilyth/etapa_lh/`, compiladas contra cada base) dan logs identicos salvo las
lineas de tiempo; (ii) un volcado del acople de E7b en 3 mundos (IID, RAFAGAS, CAMBIO; semillas de E7b) con los contrastes del resumen
(todos sus campos), el resumen armonico (p crudos, p de Holm, sobrepaso) y, por ranura cerrada, u_pit y delta_tau de las dos referencias,
impresos en hexadecimal: identicos. Falla con cualquier diferencia.

**XB.0b (A/B, condicion).** Con INTERCAMBIABLE en la misma corrida: los campos `*_pool_loo` son identicos bit a bit a los del modo
POOL_LOO, y todo lo que no es el contraste del resumen ni depende de el (ranuras, u_pit, delta_tau, P3b, sobrepaso de la curva) tambien.

**XB.0c (coherencia, condicion).** En INTERCAMBIABLE: el u del real en su estrato coincide con el u_pit de su referencia (|dif| <= 1e-12),
cada estrato tiene n_nulos_usados + 1 miembros, y los p crudos 1 y 2 del Holm de 4 son los p_bilateral de los contrastes nuevos.

**XB.1 (verdad conocida, unidad, condicion).**
- a) Enumeracion exacta: en 20 configuraciones al azar de estratos pequenos (K <= 4, tamanos 2-5, con empates) `p_min_alcanzable`
  coincide con la enumeracion completa (<= 1e-12) y el p de Monte Carlo (B 1999) cae a <= 4 sd binomiales del p exacto enumerado.
- b) Tamano: eventos intercambiables (el real es un miembro al azar de su estrato), tau enteras con empates y censura codificada en
  T_max + 1, pools de 8, K en {1, 3, 10}; 2000 replicas por K: tasa de rechazo a alfa = 0.05 <= 0.05 + 2.58 * sqrt(0.05*0.95/2000)
  (0.0626) en cada K.
- c) Evento unico censurado: estrato {T_max + 1} U 8 nulos <= T_max: `p_min_alcanzable` >= 1/9, `resoluble` falso, y p de Monte Carlo
  >= 1/9 - 3 sd.
- d) Potencia: K = 10, nulos U{1..8}, real U{4..11}: rechazo >= 0.8 en 500 replicas.
- e) Degenerado: todos los miembros iguales: S = 0, p = 1, p_min = 1.
- f) Resolucion: K = 1 con 99 nulos distintos: resoluble verdadero (p_min = 0.01); K = 1 con 50: falso (1/51).

**EXB.v (expectativa, se reporta se cumpla o no).** La rama vieja, en las mismas replicas de XB.1b, rechaza mas que 0.0626 con K = 1.

**XB.2 (Lilyth sobre E7b, mundos nuevos; IID 8101+, HETERO 8201+, RAFAGAS 8301+, 10 mundos cada uno; CAMBIO 8401+, 5, solo lectura).**
Por referencia, Holm de red sobre las series de su p de Holm de Lilyth (Onsager: p_holm[1]; sin memoria: p_holm[0]), cada rama.
- Condicion (rama INTERCAMBIABLE, Onsager, reposo IID + HETERO, 20 mundos): tasa de rechazo de red con Wilson inferior <= alfa.
- Condicion (rama INTERCAMBIABLE, Onsager, RAFAGAS): >= 8 de 10 mundos rechazan, y en ellos por lenta (S > 0).
- Expectativa (sin memoria, reposo): rechaza en >= 5 de 20 (la memoria del solape de ventanas); se lee como memoria, no como fallo.
- Expectativa (rama POOL_LOO al lado): contra Onsager en reposo rechaza al menos tanto como la intercambiable.

**XB.3 (coste, reportado).** La rama nueva anade menos del 5 % al tiempo de Lilyth por mundo.

**XB.4 (mutantes).** Diez mutantes sobre la rama nueva y el acople: real fuera del estrato, Monte Carlo del pool global, p sin el +1,
p_min sin la cola inferior, resoluble con alfa en vez de alfa/4, omision invertida, sobrepaso sin bifurcar, `*_pool_loo` sin llenar,
u del real contra el estrato que lo incluye, Holm de 4 con el p viejo. Se espera cazar >= 8; los equivalentes se declaran.
