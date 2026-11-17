# Estado de la bifurcacion al reinicio del equipo (Code, 2026-11-15, ~21:00 local)

Nada se pierde: el codigo y los resultados ya medidos estan en disco. Lo parado eran corridas de diagnostico.

## Hecho
- Expectativas selladas: `EXPECTATIVAS_BIFURCACION.md`, sha256 e9cfa2c9... (`SELLO_EXPECTATIVAS.txt`). Base previa en `base_previa/`.
- Rama en la base (`ark_cxx/base_auditada/include/ark/lilyth/`): `nulo_lilyth.hpp` (enum, campos, `evaluar_contraste_intercambiable`),
  `armonico_lilyth.hpp` (campos `*_pool_loo` del sobrepaso), `observador_lilyth.hpp` (entrada `contraste_resumen`, validez, Delta tau y
  sobrepaso por evento), README del operador. Humo de la base: compila. **Falta regenerar el manifiesto y recopiar a los agentes.**
- `red/e7b_red.cpp`: gancho de configuracion nulo por omision (copia previa en `red/build/e7b_red_antes_gancho_config.cpp`).
- XB.0 (i) CUMPLE (`identidad_baterias.log`: las 6 baterias de Lilyth identicas; fuente revision_etapas2_6, declarado).
- XB.1 (`prueba_intercambiable.log`): a-e PASAN; f FALLA por mi error al sellar (p bilateral: p_min = 2/(n+1), hacen falta 159 nulos);
  EXB.v no se cumple.
- Volcados completos (`volcado_previa/nueva/inter.txt`): IID y RAFAGAS identicos bit a bit; XB.0b identico salvo lineas H en IID y
  RAFAGAS. **Abierto:** en CAMBIO, la corrida completa 'nueva' difiere de la previa en UNA serie (aperturas 752), solo en el sobrepaso
  contra Onsager (SE ons, SC ons y su Holm). Aislado (solo CAMBIO, `determinismo/`), previa x2 y nueva dan lo mismo y coinciden con la
  completa previa. Sospecha: no determinismo del acople o de Lilyth que depende de lo corrido antes en el proceso. Las tres corridas
  completas de repeticion (`determinismo/completa_*.txt`) se pararon en el primer mundo: **truncadas, no valen**.
- XB.0c: estratos 4314, u y tamanos 0 distintos; los 11 "p distintos" son series con T_max = 2, donde Lilyth no calcula el Holm de 4
  en ninguna rama (H < 3); con Holm calculado, 0 de 19. Control mal planteado por mi: se declara.

## Pendiente, en orden
1. Repetir las tres corridas completas (`determinismo/`: nueva x2, previa x1) y aislar el no determinismo del mundo CAMBIO.
2. Mutantes: `python -I mutaciones_bifurcacion.py <scratchpad> 3 10` (necesita `volcado_previa.txt`).
3. XB.2: `build\e7b_bif.exe 32` (ya compilado; ~1 h).
4. Manifiesto (`verificar_base.py --generar`), copias de los agentes, acta `MEDIDO_BIFURCACION.md`, README de Red, bitacora, tablero y
   memoria.

## Cerrado (2026-11-16)
Los cuatro pendientes se hicieron: determinismo repetido (5 de 6 corridas completas identicas; la excepcion no se reproduce), mutantes 8/10,
XB.2 (0 fallas), manifiesto y copias de agentes. Todo en `MEDIDO_BIFURCACION.md`.
