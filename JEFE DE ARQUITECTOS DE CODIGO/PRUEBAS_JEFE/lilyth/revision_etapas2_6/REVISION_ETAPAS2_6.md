# Revision de las Etapas 2, 3.0, 3, 4, 5 y 6 del Operador de Lilyth (entrega de Argon) -- Code, 2026-11-07

Entrega: `PRUEBAS_JEFE/Proceso Operador de Lilyth/` (19 archivos; copia identica en `Carpeta del Comandante/Proceso Operador de
Lilyth/`). Copia tal cual en `argon_tal_cual/` (sha256 de cada archivo al final). Los originales de Argon no se tocan. Nada se conecto a
`ark_cxx/`: la base sigue con la Etapa 1.

## 1. Veredicto

**No se conecta a la base todavia.** La Etapa 2 esta limpia. Las etapas 4, 5 y 6 pasan sus baterias y cazan sus mutantes. Pero hay dos
defectos graves en el nulo de la Etapa 3, y las Etapas 4 y 3 pasan P5 y P2 solo gracias a uno de ellos:

| id | defecto | gravedad |
|---|---|---|
| D1 | **relleno fabricado del nulo sin memoria**: si la via (a) da menos de 16 eventos, `correr_poblacion_via_a` fabrica eventos nulos SIN el detector (ancla de 120 escrita, extremo en 20 pasos, cada 90, g inventado) | grave |
| D2 | **la guarda de solapamiento toma las replicas del propio choque como eventos ajenos** y descarta la inyeccion | grave |
| D3 | la entrega no compila contra la base | conduccion |
| D4 | cuatro controles fallan contra lo que dice el README de Argon | medio |
| D5 | parte "fisica" de S2 calculada a mano, no con el instrumento | medio |
| D6 | numeros a ojo | medio |
| D7 | la rama purgada pasaria la purga dos veces en los nulos | a vigilar |
| D8 | la preexpectativa de S1-S3 llego con los resultados | procedimiento |

## 2. Lo que se hizo para poder correrla (correcciones de conduccion, solo en la copia de revision `inc/` y las baterias copiadas)

| id | problema | correccion |
|---|---|---|
| C1 | `observador_lilyth.hpp` y `modelo_lineal_nulo.hpp` se incluyen mutuamente: ninguna cabecera compila (con ningun compilador) | `modelo_lineal_nulo.hpp` deja de incluir `observador_lilyth.hpp` y define `ModoNuloLineal` bajo la guarda que ya usaba el observador |
| C2 | `etapa2_lilyth.cpp` recibia `ObservadorSaltos`/`ObservadorLilyth` por esa inclusion circular | incluye `observador_lilyth.hpp` |
| C3 | `reloj_marea_lilyth.hpp`: `std::isfinite` dentro de `constexpr` (GCC lo tolera, MSVC no) | `calcular_tau_actividad` y `diferencia_tau_marea` pasan a `inline` |
| C4 | `etapa6_lilyth.cpp` usa campos de `EntradaPasoCxx` que no existen (`t_cierre_ms`, `t_disponible_ms`, `imputada`) | `ts_ms`, `disponible_ms`, `mascara` |
| C5 | `etapa3_0` desborda la pila (seis `ObservadorLilyth` en el marco de `main`) | `/link /STACK:1073741824` en `compilar.bat` |

Los scripts de mutantes de Argon llaman a `g++` con su propio arbol `include/`: Argon trabajo contra otra version de la base. Sus
`MUTANTES` se corren aqui con `correr_mutantes_argon.py` (cl /W4, bateria completa siempre, sin `--rapido`).

## 3. Resultados (`etapaN.log`, `mutantes_argon_2.log`, `mutantes_argon_3_6.log`)

Todas compilan con /W4 sin avisos tras C1-C5.

| etapa | bateria | lo que dice el README de Argon | mutantes de Argon |
|---|---|---|---|
| 2 modelo lineal | **L1-L10 y P8 pasan** | -- | **12/12** |
| 3.0 continuacion (E1.21) | **falla K4b y K5**; K4 (referencia bit a bit con E1.21) pasa | "18/18" | 22/22 |
| 3 nulos | **falla S1 y S2**; P1-P12, N1-N3 y S3 pasan | -- | **10/12**: M09 (sin guarda) y M10 (ventana sin T_MAX) NO se cazan |
| 4 armonico y Holm | pasa | -- | 8/8 |
| 5 fin del calor | pasa (P9: falsos 0.030 contra alfa 0.05, retraso 2 ventanas) | -- | 8/8 |
| 6 reloj doble y marea | pasa (P4 y P11) | -- | 6/6 |

## 4. Los defectos

### D1. Relleno fabricado del nulo sin memoria (`nulo_lilyth.hpp`, `correr_poblacion_via_a`)

Despues de las B realizaciones con el detector, si el nulo sin memoria tiene menos de 16 eventos y n > 200, se fabrican hasta 32 eventos
recorriendo el sustituto cada 90 pasos: ancla = mediana de 120 pasos escrita a mano, el "evento" es el extremo en los 20 pasos siguientes,
g = dY^2 / (2 sd^2) inventado, guarda 1. **Ninguno pasa por el detector.** Viola F1 ("el mismo disparo") y la regla de numeros a ojo.

Medido: con el relleno apagado (`inc_sr/`, `etapa3_sr.log`, `etapa4_sr.log`):

- **P2 falla**: Delta_tau_iid = NaN (el nulo sin memoria queda vacio para los eventos reales);
- **P5 falla**: el resumen armonico y la familia de Holm no son validos.

Es decir, las dos preexpectativas centrales (la semivida de P2 y el oscilador de P5) pasan en la entrega **solo por los eventos
fabricados**. N1 da tasa de rechazo 0.000 en la prueba iid con y sin relleno: el tamano de esa prueba no esta demostrado.

Por que pasa: el detector se calibra sobre la historia con memoria (AR(1), tau_int ~ 9); sobre el sustituto permutado (sin memoria) el
GLR, normalizado con ese tau_int, casi no alarma. La via (a) sin memoria no llega a m_min. El plano (seccion 3.3) decia "en el nulo sin
memoria no se simula", y no previo que (a) se quede corto ahi. **[DECISION del Comandante]**, propuesta de Code: cuando el nulo sin
memoria no llegue a m_min, **via (b) sin memoria = inyectar el perfil real (no sus innovaciones) en el sustituto permutado y correr el
mismo detector**. Para el nulo sin memoria, que la serie no responda es exactamente la hipotesis (lo que en el nulo lineal es el mutante
P8, aqui es el modelo correcto), y respeta el disparo. Mismo criterio de parada, misma guarda (corregida, D2).

### D2. La guarda de solapamiento descarta las inyecciones por sus propias replicas

`ejecutar_inyecciones_en_realizacion`, paso 4: toda apertura de Lilyth fuera de [t_iny, t_iny + l] y dentro de la ventana cuenta como
"evento espontaneo ajeno". Pero un choque inyectado vuelve a alarmar a los pocos pasos (es lo que motivo E1.21). Con una sonda
(`inc_dbg/`): **646 descartes en la bateria, todos con la apertura entre t_iny + l + 2 y t_iny + l + 8**, es decir, replicas del mismo
choque. Es la causa de que S2 falle por encima del borde (`sup(sol=1, adm=1)`, se esperaba 0 y 2). Efecto en el nulo (b): se quedan las
inyecciones cuyo choque no reabre, una seleccion hacia las que restituyen rapido. Correccion propuesta: las aperturas que el propio
observador marca como continuacion o replica de la ranura abierta por la inyeccion pertenecen a esa inyeccion; solo una apertura que
abre ranura nueva sin ranura viva de la inyeccion es ajena.

Y por eso M09 y M10 no se cazan: S2 ya falla en el testigo, y la parte de S2 que deberia ver la guarda no la mide (D5).

### D3. No compila contra la base

C1-C5. La inclusion circular de C1 no compila con ningun compilador, asi que lo que Argon corrio no es exactamente lo que entrego.
C4 (campos inexistentes) muestra que trabajo contra otra version de `ark/ventanas/`. Para la proxima entrega: compilar con
`compilar.bat` de esta carpeta (cl /W4 sobre `ark_cxx/base_auditada/include`), no con g++ y un arbol propio.

### D4. Controles que fallan

- **3.0 K4b y K5.** El codigo casa bit a bit con la referencia (K4 pasa); lo que falla es la expectativa del mundo de pulsos. dY = Y(t_a) - Y*
  es una sola observacion, seleccionada alta por el disparo: en el pulso de 1000 la ranura toca eta <= 0.5 en tau = 1 y consolida en 151,
  dentro del pulso. Cuando el detector reabre (1195) la ranura ya restituyo, y por E1.21 es replica, no continuacion. Resultado: 1
  continuacion en todo el mundo de pulsos. **Hallazgo para el Comandante**: con dY de una observacion, la "restitucion" de primer toque
  llega por regresion a la media del disparo, y E1.21 casi no se activa en choques sostenidos. No es un defecto del codigo (es E1.4 + F1,
  que los nulos corrigen), pero el control K5 se escribio esperando otra cosa.
- **3 S1.** Las decisiones E1.21, la CIF de replicas contra la verdad, Delta_tau_Onsager ~ 0 y la estabilidad de p_det (0.64, 0.67, 0.75;
  diferencia <= 0.15) se cumplen. Falla solo el umbral que Argon declaro en su preexpectativa, p_det >= 0.70 (0.64 y 0.67). Es una
  preexpectativa no cumplida, no un defecto del instrumento.
- **3 S2.** D2.

### D5. S2 "fisico" no usa el instrumento

La parte "con guarda contra sin guarda" no corre el instrumento sin guarda: desplaza el ancla a mano (`y_ancla + 3.5`) y recalcula eta.
Las innovaciones `e_cont` se calculan y no se usan. Ese control no puede fallar por la guarda.

### D6. Numeros a ojo

`b_mc_resumen = 499` (su README dice 1999), `B_boot = 200` en P3b, `b_max_nulo = 16`, `alarmas_minimas_calor = 8`, nivel inicial `100.0`
en la via (b) de NIVEL, `MAX_INTENTOS_VIA_B = 256` (tope de memoria: valido si se cuenta), 80 iteraciones fijas en las raices del AR.

### D7. Purga en los nulos (a vigilar)

Las realizaciones nulas entran al operador delta con fase 0 y un dia distinto en cada paso (`b * 1000 + t`, `20000 + t`). Con purga activa,
la rama purgada se purgaria dos veces (el modelo ya esta ajustado sobre la serie purgada) y con la fase equivocada. Ningun control corre
nulos con purga.

### D8. Procedimiento

La preexpectativa de S1-S3 (`readme_etapa3_lilyth.md`) llego con el codigo y los resultados: no pudo sellarse antes. Su sha256 queda abajo
como registro, no como sello previo.

### Observacion sobre los mutantes

Los 18 mutantes de las Etapas 3-6 estan escritos como interruptores `const bool` en las cabeceras de produccion (ramas muertas que el
script activa). Funciona, pero solo prueba las fallas que el autor previo, y en sus scripts el modo `--rapido` los juzgaba solo con pruebas
hechas a la medida de cada interruptor. Aqui se corrieron siempre con la bateria completa.

## 5. Que falta para conectar

1. Decision sobre D1 (propuesta: via (b) sin memoria por inyeccion del perfil en el sustituto) y quitar el relleno.
2. Corregir D2; rehacer S2 con el instrumento sin guarda (D5); que M09 y M10 queden cazados.
3. Derivar o declarar D6; un control de nulos con purga activa (D7).
4. Rehacer K5 con la expectativa correcta, o decidir si el primer toque debe exigir algo mas que una observacion.
5. Volver a correr todo con `compilar.bat` de esta carpeta.

## 6. sha256 de la entrega (argon_tal_cual/)

Ver `sha256_entrega.txt`.
