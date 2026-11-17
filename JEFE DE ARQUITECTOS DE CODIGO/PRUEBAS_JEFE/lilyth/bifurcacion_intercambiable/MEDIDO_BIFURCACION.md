# Bifurcacion del contraste del resumen de Lilyth (POOL_LOO / INTERCAMBIABLE) -- MEDIDO

Code (Opus 5.5), 2026-11-15/16. Orden del Comandante: bifurcar el resumen de Lilyth y leer la recomposicion de la Red con AMBAS
referencias (sin memoria y Onsager), puntos 1-6 aprobados. Expectativas selladas antes del codigo: `EXPECTATIVAS_BIFURCACION.md`
(sha256 e9cfa2c9ddb2cea568b9b9d199f9cf30dd35256b4dfb4e42b41dab2995521beb, `SELLO_EXPECTATIVAS.txt`). Solo sinteticos. Nada de esto es
para operar.

## 1. Lo construido (en `ark_cxx/base_auditada/include/ark/lilyth/` del worktree, SIN aplicar a la principal)

- `nulo_lilyth.hpp`: `enum class ContrasteResumenLilyth { POOL_LOO = 0, INTERCAMBIABLE = 1 }`; `ContrasteResumenNulo` gana
  `intercambiable`, `pool_min_estrato`, `p_min_alcanzable`, `resoluble`; `FichaNuloRepresentacion` gana `contraste_iid_pool_loo` y
  `contraste_onsager_pool_loo`; funcion nueva `detalle_nulo::evaluar_contraste_intercambiable(estratos, B, semilla, umbral)`.
- `armonico_lilyth.hpp`: `ContrasteSobrepasoLilyth` gana `sobrepaso_por_evento_{iid,onsager}_pool_loo`.
- `observador_lilyth.hpp`: entrada `ConfiguracionLilyth::contraste_resumen` (POOL_LOO por omision, al final de la estructura) y su
  validez; en `evaluar_nulos_representacion` los estratos {tau del real} U {tau de sus nulos emparejados} y la rama nueva en el Delta tau;
  en `contraste_sob` los estratos {s_k} U sn y la rama nueva en el sobrepaso por evento. Con POOL_LOO no se arma ningun estrato: el
  camino es el de antes.
- `README_OPERADOR_LILYTH.md`: la decision documentada.
- Manifiesto regenerado (70 archivos); las dos copias de los agentes recopiadas enteras y las tres verifican identicas; humo de la base compila.
- Fuera de la base: `PRUEBAS_JEFE/red/e7b_red.cpp` gana un gancho de configuracion de Lilyth nulo por omision (la bateria no lo usa; la
  copia previa queda en `red/build/e7b_red_antes_gancho_config.cpp`).

## 2. Declarado al concretar (despues del sello, antes de leer resultados salvo donde se dice)

1. **XB.0 (i) usa las baterias de `revision_etapas2_6/`, no las de `etapa_lh/`.** La de la Etapa 4 de `etapa_lh/` es de la rama L+H
   (Holm de 6), que no esta en la base, y no compila contra ella; las de `revision_etapas2_6/` son las que acompanaron a la base.
2. **XB.1f estaba mal sellada (error mio, aritmetico).** El p es BILATERAL: con un solo evento por encima de todos sus n nulos, tambien
   el miembro de abajo alcanza |S| maxima, asi que p_min = 2/(n + 1), no 1/(n + 1). Con 99 nulos p_min = 0.02 y no resuelve. Resolver
   alfa/4 pide n >= 159 (medido: 0.0125 con 159, 0.0126 con 158). La cifra "n >= 79" de `PRUEBAS_JEFE/red/MEDIDO_E7B.md` tenia el mismo
   error y se corrige alli.
3. **XB.0c estaba mal planteado en un caso (error mio).** Con T_max efectivo < 3, `evaluar_armonico_representacion` sale antes de la
   familia de Holm de 4 (H < 3) en LAS DOS ramas: p_crudos queda NaN aunque el contraste tenga p. Eso no es un defecto de la rama nueva.
4. **El script de mutantes tuvo un defecto propio en su primera pasada** (`mutaciones_bifurcacion_primera_defectuosa.log`): al leer los
   contadores de XB.0c recogia tambien el "4" de "Holm de 4" y daba a todo mutante por cazado. Corregido (compara contra los contadores
   del testigo) y corrido otra vez.

## 3. Resultados por la regla sellada

| | resultado |
|---|---|
| **XB.0 (i)** baterias de Lilyth 2-6, base previa contra nueva | **CUMPLE**: las 6 identicas (`identidad_baterias.log`) |
| **XB.0 (ii)** volcado del acople de E7b en IID, RAFAGAS y CAMBIO | **CUMPLE en 5 de 6 corridas completas, con una excepcion esporadica** (abajo) |
| **XB.0b** INTERCAMBIABLE con la vieja al lado | **CUMPLE**: lo que no es Holm, identico bit a bit al volcado de la base previa en los 3 mundos |
| **XB.0c** coherencia | u del real contra u_pit: **0 distintos de 4 314 estratos**; tamano del estrato: **0 distintos**; p del Holm de 4 contra el contraste nuevo: 0 distintos en las 19 series con Holm (T_max >= 3); los 11 restantes son T_max = 2 (punto 2.3) |
| **XB.1a** enumeracion exacta | **PASA**: p_min 0/20 distinto; Monte Carlo a menos de 4 sd en 20/20 (peor 2.5) |
| **XB.1b** tamano con pools de 8, tau empatadas y censura | **PASA**: K = 1: 0.000, K = 3: 0.038, K = 10: 0.046 (techo 0.0626) |
| **EXB.v** la vieja anticonservadora con K = 1 | **NO SE CUMPLE**: 0.0285 / 0.040 / 0.048 en ese mundo (ver 4) |
| **XB.1c** evento unico censurado | **PASA**: p_min 2/9 = 0.222, no resoluble, p 0.235 (vieja: 0.135) |
| **XB.1d** potencia | **PASA**: 0.928 |
| **XB.1e** degenerado | **PASA**: S 0, p 1, p_min 1 |
| **XB.1f** resolucion | **FALLA** (sellada mal, punto 2.2) |
| **XB.4** mutantes | **8 de 10** (ver 6) |

**La excepcion de XB.0 (ii).** La primera corrida completa con la base nueva (`corrida` de las 19:23) dio, en el mundo CAMBIO, UNA serie
distinta de la previa (la de 752 aperturas), y solo en el sobrepaso contra Onsager (SE ons, SC ons y su Holm): Delta tau, ranuras, u_pit
y P3b identicos. Repetido: el mundo CAMBIO aislado (previa dos veces y nueva) identico; tres corridas completas mas (nueva dos veces y
previa) identicas a la previa en los 3 mundos; y la corrida INTERCAMBIABLE que corrio a la vez que la excepcion tambien coincide con la
previa en lo que no es Holm. Es decir: 1 de 6 corridas completas, no reproducible, en una cantidad que con POOL_LOO calcula el codigo
viejo linea por linea. Sin estado estatico, relojes ni miembros `mutable` en Lilyth ni en el operador delta (revisado). **No aislado.**
Se declara como no determinismo esporadico del acople o de Lilyth; no se puede excluir del todo que lo cause la rama.

## 4. Lecturas

- **La rama nueva tiene tamano nominal en verdad conocida, con pools pequenos, empates y censura**, y potencia (0.93 con K = 10). Con
  K = 1 es conservadora (0.000): con 9 miembros por estrato el p bilateral no baja de 2/9.
- **La vieja no salio anticonservadora en el mundo de XB.1b** (EXB.v). El mecanismo de E7b (el real por encima de todo su nulo, con el
  nulo empatado en sus DOS extremos) no se da en un mundo de tau exponenciales: el miembro mas bajo del pool suele tener u = 0, y entonces
  el pool si alcanza |S| = 0.5. En el caso de XB.1c la vieja da 0.135 por eso. El defecto existe (E7b lo midio: p = 0.0005 con un evento)
  pero exige empates en los dos extremos; no es general.
- **Resolucion:** con un evento hace falta n >= 159 para que un rechazo pueda pasar el primer escalon del Holm de 4. Con m_min = 8,
  las familias de un evento no pueden rechazar con la rama nueva. Se reporta en `resoluble`.

## 5. XB.2 y XB.3 (E7b en mundos nuevos, las dos ramas, las dos referencias)

`e7b_bif.exe 32` -> `e7b_bif.log` (4 037 s, 0 fallas). Mundos nuevos: IID 8101-8110, HETERO 8201-8210, RAFAGAS 8301-8310, CAMBIO
8401-8405. Una sola corrida de Lilyth en INTERCAMBIABLE por mundo; la rama vieja sale de los `*_pool_loo` con su Holm de 4 rehecho por
`aplicar_familia_holm_cuatro`. Por referencia y por rama, Holm de red sobre el p de Holm de Lilyth de cada serie valida.

**Holm de red, mundos que rechazan (lenta = alguna serie rechazada con S > 0):**

| mundo | Onsager, nueva | Onsager, vieja | sin memoria, nueva | sin memoria, vieja |
|---|---|---|---|---|
| IID (10) | **0/10** | 0/10 | 9/10 (9 lentas) | 10/10 (10) |
| HETERO (10) | **1/10** (1 lenta) | **4/10** (3) | 8/10 (8) | 10/10 (10) |
| RAFAGAS (10) | **10/10** (10 lentas) | 10/10 (10) | 10/10 (10) | 10/10 (10) |
| CAMBIO (5, lectura) | 5/5 (5) | 5/5 (5) | 5/5 (5) | 5/5 (5) |

**Series que rechazan su Holm de 4** (Onsager): IID nueva 1/50, vieja 2/50; HETERO nueva 1/57, vieja 5/57; RAFAGAS 55/63 las dos.
S Onsager de la nueva, mediana: IID 0.000, HETERO +0.035, RAFAGAS **+0.263**, CAMBIO +0.063. Series con eventos no resolubles
(p_min > alfa/4): 55 % en IID, 54 % en HETERO (pools de 8 y familias de pocos eventos), 2 % en RAFAGAS, 0 % en CAMBIO.

| expectativa | resultado |
|---|---|
| **XB.2 condicion:** nueva contra Onsager en reposo, Wilson inferior <= alfa | **CUMPLE**: 1/20 = 0.05 [0.009, 0.236] |
| **XB.2 condicion:** nueva contra Onsager en rafagas, >= 8/10 y todos por lenta | **CUMPLE**: 10/10, 10 por lenta |
| **XB.2e** sin memoria en reposo rechaza en >= 5/20 | **SE CUMPLE**: 17/20 |
| **XB.2e'** la vieja contra Onsager en reposo rechaza al menos tanto como la nueva | **SE CUMPLE**: vieja 4/20, nueva 1/20 |
| **XB.3** coste de la rama nueva | **0.004 %** del tiempo de la serie (los dos contrastes rehechos aparte) |

**Lectura (las dos referencias, como pidio el Comandante):**
- **Contra Onsager**, la red se recompone en reposo como su propia memoria lineal (tamano nominal con la rama nueva) y tras las rafagas
  se recompone MAS LENTO que ella (10/10, S +0.26). Es la respuesta de E7b, ahora confirmada en mundos nuevos y por la regla sellada.
- **Contra el nulo sin memoria** rechaza tambien en reposo (17/20 nueva, 20/20 vieja): mide la memoria del solape de ventanas de las
  series de E6b (12 emisiones), no la recomposicion. Las dos lecturas juntas separan "la serie tiene memoria" (sin memoria) de "tras
  el choque tarda mas que su memoria" (Onsager).
- **La rama vieja es anticonservadora en reposo HETERO** (4/10 mundos, 5/57 series contra 1/10 y 1/57 de la nueva). En IID no se ve
  (0/10 las dos). Con potencia (rafagas, cambio) las dos ramas coinciden.
- Mas de la mitad de las series en reposo no pueden rechazar con la rama nueva aunque quisieran (no resolubles): la decision sobre
  m_min (subirlo a >= 159 eventos en el nulo emparejado, o no) sigue abierta.

## 6. Mutantes (`mutaciones_bifurcacion.py`, `mutaciones_bifurcacion.log`)

Base copiada al scratchpad (ruta corta) con un cambio; prueba de unidad + volcado IID en `nueva` e `inter`.

| | mutante | resultado |
|---|---|---|
| M1 | real fuera del estrato en el Monte Carlo | CAZADO (XB.1a, 1b) |
| M2 | Monte Carlo del pool global | CAZADO (XB.1a) |
| M3 | p sin el +1 | **NO CAZADO**: corre el p en <= 1/2000, por debajo de la tolerancia de XB.1a; debilidad declarada |
| M4 | p_min sin la cola inferior | CAZADO (XB.1a) |
| M5 | resoluble con alfa en vez de alfa/4 | **NO CAZADO, equivalente en el mundo IID** (ningun p_min cae entre alfa/4 y alfa) |
| M6 | omision invertida | CAZADO (XB.0, XB.0b) |
| M7 | sobrepaso por evento sin bifurcar | CAZADO (XB.0b) |
| M8 | `contraste_iid_pool_loo` sin llenar | CAZADO (XB.0b) |
| M9 | u del real contra el estrato que lo incluye | CAZADO (XB.1a) |
| M10 | Holm de 4 con el p viejo en INTERCAMBIABLE | CAZADO (XB.0c) |

## 7. Limitaciones

- Solo sinteticos. La rama nueva supone estratos independientes entre eventos (comparten nulos en la via (a)), igual que la vieja.
- El no determinismo esporadico de 3 queda sin aislar.
- M3 y M5 no los caza ningun control.
