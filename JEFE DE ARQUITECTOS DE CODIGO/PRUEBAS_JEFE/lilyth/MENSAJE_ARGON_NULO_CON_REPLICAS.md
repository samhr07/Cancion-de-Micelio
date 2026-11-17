# Para Argon: un tercer nulo de Lilyth, lineal con replicas (Hawkes por signo)

Code (Jefe de Arquitectos de Codigo, Opus 5.5), 2026-11-07, por orden del Comandante. **Directriz de trabajo: las cuatro decisiones
(H.1-H.4) las cerro el Comandante el 2026-11-08** y estan en la seccion 3. Contexto: `revision_etapas2_6/CORRECCION_ETAPAS2_6.md`
(seccion 3) y `README_OPERADOR_LILYTH.md` (E1.24).

## 1. El problema, medido

En el control S1 de la Etapa 3 (replicas de Hawkes sobre un AR(1) de calma conocido), el contraste contra el nulo de Onsager rechaza
con razon de ramificacion n = 0.6 (p 0.005). No es un defecto del instrumento:

- por E1.21, una alarma del mismo signo con la ranura viva sin restituir es una continuacion: la ranura no se corta y sigue midiendo;
- las hijas de Hawkes llegan despues de t_a, fuera del perfil de l* pasos que la via (b) inyecta;
- en el nulo de Onsager la serie solo responde a ese perfil, sin hijas, y restituye antes.

Asi que Delta_tau_Onsager suma dos cosas: la autoexcitacion (choques nuevos que llegan despues) y lo que no explica ni la memoria
lineal ni la autoexcitacion. Para separarlas hace falta un nulo que tenga replicas.

Una continuacion puede ser una hija (un choque nuevo) o el mismo choque que persiste mas alla del perfil. El nulo con replicas modela
lo primero; lo segundo es justo lo que debe quedar en el exceso.

## 2. El nulo L+H

**El modelo lineal de la Etapa 2 con un proceso de replicas encima.** En la via (b), para cada inyeccion:

1. Se simula la calma con el `ModeloLinealNulo` de la representacion (E1.22) y se inyecta el perfil real, como hoy.
2. Se generan hijas en t_a + s con un proceso de Hawkes de DOS TIPOS, mismo signo y signo contrario (H.2). Cada hija puede tener las
   suyas: un cumulo de ramificacion, con tope de memoria y su desborde contado.
3. Cada hija es un perfil de SU signo, tomado al azar de los perfiles de calibracion de ese signo, y entra por las innovaciones
   (filtro inverso A(L)), como la madre.
4. Corre el mismo detector, con la misma guarda de solape y la misma parada ceil(m_min / p_det). Las hijas son parte de la
   inyeccion: sus aperturas no son "ajenas". Una hija del mismo signo con la ranura viva sin restituir es continuacion; una de signo
   contrario es replica y corta, exactamente igual que en los eventos reales (E1.21). Las inyecciones cortadas quedan fuera, simetrico
   con los reales cortados.

Con eso:

```
Delta_tau_Onsager = [Delta_tau_Onsager - Delta_tau_LH]  +  Delta_tau_LH
                      autoexcitacion                       lo que ni la memoria lineal ni las replicas explican
```

**Ajuste del proceso de replicas, sin numeros a mano:**

- **Datos.** La secuencia de alarmas del operador delta CONGELADO (umbrales de produccion) reproducido sobre la historia de
  calibracion. Nunca los eventos que se evaluan: eso es fuga de validacion.
- **Dos tipos por signo con marca de escala (H.2, H.3).** El tipo de cada alarma es su signo (+ o -); su escala es una marca, no un
  proceso aparte. Nucleo matricial por pares (madre, hija) de tipos:

  ```
  phi_{s,s'}(t) = sum_j a_{s,s',j} exp(-t / tau_j),  a >= 0,  s, s' en {+, -}
  ```

  Las tau_j son las constantes de la escalera del propio operador delta (tau_j = 2^j pasos de la barra base, j en las escalas
  calibradas). Con las tau fijas la log-verosimilitud es concava en (mu, a), con un maximo unico. Se reportan por separado la excitacion
  del mismo signo (a_{+,+}, a_{-,-}) y la de signo contrario (a_{+,-}, a_{-,+}).
- **Tasa base por bin de fase de la purga** (los `bins_fase` de la `ConfiguracionObservador`), no constante, para cada tipo. Con base
  constante, una mezcla de Poisson de tasa variable se lee como autoexcitacion: Filimonov y Sornette (2015). Ya lo medimos aqui: en la
  calibracion L3 contra 1 minuto, la mitad de las ventanas con eta >= 1 las fabricaba la base constante.
- **Ramificacion.** La matriz de ramificacion es N_{s,s'} = sum_j a_{s,s',j} tau_j, y lo que decide es su radio espectral rho(N). Con
  rho(N) >= 1, o con su cota superior del IC >= 1, el proceso es NO IDENTIFICADO.
- **Biblioteca de perfiles** (marcas de las hijas), por signo: los perfiles de las alarmas de calibracion de cada signo.

## 3. Decisiones del Comandante (2026-11-08)

| # | decision | consecuencia para el codigo |
|---|---|---|
| H.1 | **la referencia L+H entra en la familia de Holm**: de 4 a 6 (p_Delta_tau_LH y p_sobrepaso_LH) | `ResumenArmonicoYHolm` pasa a 6 p. b_mc tiene que cumplir 1/(B+1) <= alfa/6 (1999 lo cumple con alfa = 0.05); `configuracion_lilyth_nulo_valida` lo comprueba con 6 |
| H.2 | **las hijas se clasifican de mismo signo y de signo contrario** | Hawkes de dos tipos con nucleo cruzado (seccion 2). Las hijas de signo contrario entran en el nulo como replicas que cortan, igual que en los eventos reales |
| H.3 | **ajuste agrupado por signo, con la escala como marca** | un solo ajuste de dos tipos; no un Hawkes por escala |
| H.4 | **NO DECIDIBLE si no pasa de sus nulos** | ver abajo |

**H.4, como se aplica.** El ajuste de replicas tiene que superar su propio nulo antes de usarse:

- **Nulo del ajuste:** Poisson con la misma tasa base por fase y por tipo, sin excitacion. Estadistico: la razon de verosimilitud
  Hawkes contra Poisson-por-fase. Su distribucion se simula bajo el nulo ajustado (bootstrap parametrico), con B derivado de alfa y no
  escrito.
- Si la razon de verosimilitud no supera su nulo a alfa, o si la calibracion no tiene alarmas suficientes para calcularla (que es la
  misma situacion: no puede superarlo), **la referencia L+H es NO DECIDIBLE** para toda la representacion. No se sustituye por n = 0,
  porque eso la volveria Onsager.
- Lo mismo con rho(N) >= 1: NO IDENTIFICADO.
- La ficha reporta el estadistico, su nulo y el recuento de alarmas.

**Consecuencia para Holm, aplicada por Code salvo que el Comandante diga otra cosa.** `holm.hpp` falla de forma dura con NaN, y la
familia se fija ANTES de ver el resultado (si encogiera de 6 a 4 al no decidirse, el alfa dependeria del dato). Por eso, con L+H NO
DECIDIBLE, sus dos p entran con **p = 1, marcados NO DECIDIBLE**. Eso es conservador: no infla el alfa y la familia sigue siendo de 6.

## 4. Controles que la bateria tiene que traer (sinteticos con verdad conocida)

| id | mundo | lo que se exige |
|---|---|---|
| R1 | AR(1) sin replicas | el ajuste NO pasa su nulo en una fraccion >= 1 - alfa de las corridas (tamano del contraste del ajuste); L+H queda NO DECIDIBLE con sus p = 1; Onsager con tamano <= alfa |
| R2 | el mundo S1 con n = 0.3 y 0.6, solo del mismo signo | el ajuste pasa su nulo; rho(N) recuperado dentro de su error; los terminos de signo contrario compatibles con 0; L+H con tamano <= alfa; Onsager sigue rechazando en n = 0.6 |
| R2b | replicas de signo contrario con verdad conocida | los a_{+,-} y a_{-,+} recuperados; en el nulo esas hijas cortan la inyeccion; los cortes del nulo y de los reales estan en la misma proporcion dentro de su banda |
| R3 | TAR no lineal con replicas | L+H rechaza: la no linealidad no se la come la autoexcitacion |
| R4 | mezcla de Poisson con tasa por fase, sin excitacion | con base por fase el ajuste no pasa su nulo (NO DECIDIBLE); con base constante sale excitacion falsa (el control se reporta: es la razon de la decision) |
| R5 | critico, rho(N) = 1 | NO IDENTIFICADO, contado; ningun evento con referencia L+H |
| R6 | causalidad | el ajuste solo ve la calibracion: el mutante que ajusta con los eventos evaluados TIENE que fallar |
| R7 | solape | la guarda trata a las hijas como parte de su inyeccion; S2 sigue pasando |
| R8 | identidad | con a = 0 forzado (sin la compuerta de H.4), L+H = Onsager bit a bit |
| R9 | Holm | con L+H NO DECIDIBLE, la familia es de 6 con p = 1 en esos dos puestos y las otras cuatro decisiones de Holm no cambian respecto de una corrida de referencia escrita aparte |

Mutantes minimos:
- ajuste con los eventos evaluados;
- base constante;
- un solo tipo, que funde los dos signos;
- hijas de signo contrario tratadas como continuacion;
- hijas sin descendencia;
- marcas tomadas de los eventos evaluados;
- rho(N) >= 1 aceptado;
- tau_j escritos a mano en vez de la escalera;
- L+H usado sin pasar su nulo;
- familia encogida a 4 cuando L+H no se decide.

## 5. Como entregar (lo que aprendimos en la revision de las Etapas 2-6)

- Compilar con MSVC contra `ark_cxx/base_auditada/include` (la copia de `ARQUITECTO DE DATOS/nucleo_cxx_auditado` es identica): cl
  /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX. La entrega anterior no compilaba contra la base (C1-C5).
- Ningun relleno ni respaldo silencioso: lo que no alcanza m_min o no pasa su nulo queda sin referencia y se cuenta.
- Ningun numero escrito: lo que no sale de la teoria o de la calibracion es una entrada declarada.
- Los mutantes, como hasta ahora (interruptores `const bool`), corridos con la bateria completa.
- La implementacion vive en `nulo_lilyth.hpp` (via (b)), `armonico_lilyth.hpp` (familia de 6) y una cabecera nueva para el ajuste:
  `ark/lilyth/replicas_lilyth.hpp`. Partir de las cabeceras CORREGIDAS de la base, no de la entrega anterior.
