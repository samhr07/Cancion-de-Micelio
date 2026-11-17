# Operador de Lilyth (`ark/lilyth/`) -- restitucion tras un evento del operador delta. ETAPAS 1-6

Argon (Auditoria Externa) escribio las Etapas 1-6. Code (Jefe de Arquitectos de Codigo, Opus 5.5) las reviso contra el operador delta
real, corrigio lo que no pasaba y las conecto a la base: la Etapa 1 el 2026-11-05 y las Etapas 2-6 el 2026-11-07.

- Diseno, decisiones del Comandante (E1.1-E1.22) y preexpectativas: `JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md`
  (canonico; este README lo sustituye al cerrar la Etapa 8).
- Revision de la Etapa 1: `PRUEBAS_JEFE/lilyth/revision_etapa1/REVISION_ETAPA1.md`.
- Revision y correccion de las Etapas 2-6 (defectos D1-D8, conduccion C1-C5, controles y mutantes):
  `PRUEBAS_JEFE/lilyth/revision_etapas2_6/REVISION_ETAPAS2_6.md` y `CORRECCION_ETAPAS2_6.md`.

Es telemetria: no decide operar.

## Cabeceras

| cabecera | etapa | que hace |
|---|---|---|
| `observador_lilyth.hpp` | 1, 3, 4, 5, 6 | el orquestador: ranuras, ancla, eta(tau), tau_rest, continuacion E1.21; `calibrar` y `completar_nulos` |
| `modelo_lineal_nulo.hpp` | 2 | nulo lineal de Onsager: AR por Levinson-Durbin o ARFIMA(p, d, 0), filtro inverso A(L), simulacion (E1.22) |
| `nulo_lilyth.hpp` | 3 | vias (a) por excedencia emparejada y (b) por inyeccion, guarda de solape, espejo de la purga, PIT LOO, P3b |
| `armonico_lilyth.hpp` | 4 | Prony sobre la curva media, sobrepaso armonico emparejado, familia de Holm de 4 (E1.16) |
| `fin_calor_lilyth.hpp` | 5 | fin del calor calibrado desde las alarmas de temperatura del operador delta (mu_calor de las l* ventanas) |
| `reloj_marea_lilyth.hpp` | 6 | relojes de actividad y de pared, contraste de marea DESCRIPTIVO (fuera de la familia de Holm) |

## Que hace

```
ObservadorSaltos (calibrado)  --Salida por referencia constante-->  ObservadorLilyth
   apertura en t_a = alarma.t (espejo de la regla del operador delta); ancla Y* congelada antes del cambio
   eta(tau) = (Y(t_a + tau) - Y*) / (Y(t_a) - Y*) en las ramas PURGADA y CRUDA
   tau_rest_primer_toque, tau_rest_consolidado, replicas y continuaciones (E1.21), Aalen-Johansen por signo
calibrar:        modelos lineales por representacion (E1.22), via (a), T_max derivado (E1.9), poblaciones nulas
completar_nulos: cada ranura cerrada contra el nulo sin memoria y el de Onsager (u = PIT), resumen S = media(u - 1/2) con Monte Carlo,
                 armonico emparejado, Holm de 4, marea descriptiva
```

## Uso minimo

```cpp
#include "ark/lilyth/observador_lilyth.hpp"
using Obs = ark::metrologia::ObservadorSaltos<>;
auto delta = std::make_unique<Obs>();
delta->calibrar(cd, x_cal, {}, v_cal, f_cal, d_cal);          // cd.arl0_pasos_nativos FIJADO: sin el, Lilyth no calibra
auto lil = std::make_unique<ark::lilyth::ObservadorLilyth<Obs>>();   // grande: en el monticulo
ark::lilyth::ConfiguracionLilyth cl{};                        // t_max_efectivo lo deriva calibrar
lil->calibrar(*delta, cl, cd, x_cal, v_cal, f_cal, d_cal);    // la MISMA historia y la MISMA configuracion que el operador delta
for (...) {
    const auto s = delta->paso(x, NaN, valido, fase, dia);
    const auto o = lil->paso(s, x, valido, d_nu, fase, dia);  // las MISMAS entradas que el operador delta
    for (std::size_t i = 0; i < o.ranuras_cerradas; ++i) usar(lil->cerrada(i));
}
lil->completar_nulos();                                        // fuera del ciclo caliente
const auto& h = lil->resumen_armonico_purgada();               // h.p_crudos, h.p_holm, h.rechaza_holm (familia de 4)
const auto& m = lil->resumen_marea();                          // descriptivo, sin p
```

Por filas: `lil->paso_filas(salida_de_ObservadorSaltosFilas::paso_fila, fase_de_barra, fase, dia)` (nu y pared salen de la barra).
`configurar` (Etapa 1, sin nulos) sigue disponible para medir solo la maquinaria, con `t_max_efectivo` fijado a mano.

## Contrato

- No posee ni toca al operador delta. Su reloj es el `Salida.t` del operador delta: un paso saltado se cuenta (`pasos_desalineados`) y
  no se mide.
- N_REG, B_MAX y J_MAX salen del tipo del operador delta (`RasgosDelta`). La integracion y la purga salen de su
  `ConfiguracionObservador`, que se pasa explicita (se retiro la sobrecarga que la deducia).
- La purga de la rama PURGADA es un espejo de la `PurgaDiurna`. El espejo de los nulos (`EspejoPurgaNulo`) la reproduce bit a bit
  sobre la historia de calibracion (`espejo_purga_discrepancias()` tiene que dar 0) y lleva cada nulo al espacio crudo con las fases y
  los dias de la calibracion (`llevar_a_crudo`).
- Sin numeros a mano: B de la via (a) = `familia::realizaciones` con el ARL0 del operador delta; m_min = `familia::m_min(alfa)` en la
  via (a), en la (b), en el pool LOO y en el fin del calor; `b_mc_resumen` (1999 por omision) tiene que cumplir 1/(B+1) <= alfa/4 o
  la configuracion no es valida. RANURAS, T_MAX, HIST_MAX y `b_max_nulo` son topes de memoria con su desborde contado.
- Cada perdida tiene su motivo (`DesenlaceRanura`) y cada inyeccion su destino (`DestinoInyeccion`: ADMITIDA, CORTADA, NO_MEDIDA).

## Decisiones del Comandante implementadas

- **E1.21 continuacion:** una alarma del mismo signo con la ranura viva sin restituir es continuacion: ni corta ni abre ranura.
- **E1.22:** el modelo lineal se ajusta sobre la serie que Lilyth mide (pasos de la barra base), con su propio K por
  `est::tau_integrado` (k_max = n/4); un modelo por representacion; en NIVEL modela dY. Innovaciones: residuos remuestreados por
  omision, gaussianas como alternativa.
- **Prueba del resumen:** S = media(u_k - 1/2) con PIT LOO del pool nulo y Monte Carlo; familia de Holm de 4.
- **Bifurcacion de la prueba del resumen (2026-11-15, regla 6):** `ConfiguracionLilyth::contraste_resumen` =
  `ContrasteResumenLilyth::POOL_LOO` (por omision, la de arriba, identica bit a bit) o `INTERCAMBIABLE`: cada evento real forma un
  estrato con sus nulos emparejados, el u de cada miembro es su PIT contra el resto del estrato y el Monte Carlo sortea un miembro por
  estrato. Cubre el Delta tau y el sobrepaso por evento. POOL_LOO es anticonservadora con pools de m_min eventos y tau empatadas (el
  real por encima de todo su nulo cae fuera del soporte del pool: `PRUEBAS_JEFE/red/MEDIDO_E7B.md`). Con INTERCAMBIABLE la ficha
  lleva la rama vieja al lado (`contraste_*_pool_loo`, `sobrepaso_por_evento_*_pool_loo`) y cada contraste declara
  `pool_min_estrato`, `p_min_alcanzable` (exacto) y `resoluble` (p_min <= alfa/4). Con un solo evento y extremos simetricos el p
  bilateral minimo es 2/(n + 1) (n >= 159 nulos para el primer escalon del Holm de 4); con empates en un extremo (censura en T_max + 1)
  los dos extremos dan |S| distinto y el minimo es 1/(n + 1) o mas, segun los empates: `p_min_alcanzable` es la enumeracion exacta y es
  lo unico que decide. Acta: `PRUEBAS_JEFE/lilyth/bifurcacion_intercambiable/MEDIDO_BIFURCACION.md`.
- **Pool derivado (2026-11-16, regla 6):** `ConfiguracionLilyth::pool_nulo` = `PoolNuloLilyth::M_MIN` (por omision, identica bit a bit)
  o `DERIVADO`: con K eventos candidatos (RESTITUIDA o CENSURADA_HORIZONTE) el objetivo por evento es
  n* = ceil((2 * 4 / alfa)^(1/K)) - 1 (K = 1: 159; 2: 12; 3: 5; 4: 3), sin ningun numero nuevo; umbral = max(m_min, n*): la via (a) se usa
  si la excedencia llega al umbral y la via (b) persigue el umbral. Solo agrega nulos: si la via (b) no llega y la excedencia ya tenia
  m_min, se queda la via (a). La ficha dice `pool_derivado`, `k_candidatos_pool` y `objetivo_pool`. Limite medido: con un solo evento la
  potencia sigue baja porque el real tiene que superar a todos sus nulos y, si sale censurado, empata con los nulos censurados (con
  T_max = cuantil 1 - alfa del nulo lineal quedan ~alfa de ellos). Acta: `PRUEBAS_JEFE/lilyth/pool_y_fuente_nulos/MEDIDO_POOL_Y_FUENTE.md`.
- **Fuente de los nulos (2026-11-16, regla 6):** `calibrar(..., nu_cal, fuente_nulo)` acepta una mascara opcional del largo de `x_cal`.
  Vacia: hoy, bit a bit. Con mascara, los segmentos del nulo son los tramos validos Y dentro de la mascara (un 0 corta el segmento); la
  purga, el largo de cada realizacion y B no cambian. La ficha dice `fuente_externa` y `n_fuente`. Las lecturas del Comandante son
  mascaras del llamador sobre una historia de calibracion con viento: CALMA (fuera de las ventanas del catalogo), VIENTO (dentro) y
  COMPLETA (todo). Cada una contesta otra pregunta: "mas lenta que la calma", "que otras tormentas" o "que su memoria normal". El evento
  contrastado nunca entra en su propio nulo si la historia de calibracion y el flujo de operacion son disjuntos: responsabilidad del
  llamador. El operador delta no se toca (se calibra en calma). Las ventanas salen del bufer derivado (NOTA_LILYTH_BUFER_10_SESIONES.md).
- **Horizonte T_max, enmienda de E1.9 (2026-11-16, regla 6):** `ConfiguracionLilyth::horizonte` = `HorizonteLilyth::CUANTIL_ALFA` (por
  omision, E1.9 de Argon, bit a bit) o `CUANTIL_RESOLUCION`: T_max = cuantil 1 - alfa/8 del nulo lineal y guarda de censura en alfa/8
  (la ficha dice `nivel_t_max`). Por que: con E1.9 queda ~alfa de los nulos censurados en T_max + 1, un evento real lento censurado empata
  con todos y su p minimo queda en el orden de esa fraccion aunque el pool sea enorme; resolver alfa/4 en dos colas pide c <= alfa/8.
  Requisito medido: el cuantil 1 - alfa/8 solo se estima con >= 8/alfa = 160 nulos lineales en la via (a); con la B de hoy (objetivo
  m_min) salen ~11 y el T_max apenas cambia. Acta: `PRUEBAS_JEFE/lilyth/horizonte_resolucion/MEDIDO_HORIZONTE.md`.
- **Via (b) hasta ceil(m_min / p_det)** con guarda de solape: solo una apertura en la ventana previa [ini, t_iny) contamina; una
  reapertura despues del perfil es replica (CORTADA) y se excluye, simetrico con los eventos reales cortados.
- **Nulo sin memoria por la via (b) (2026-11-07):** cuando la excedencia emparejada no alcanza m_min, el perfil real se escribe en un
  sustituto permutado y corre el mismo detector. No se rellena con eventos fabricados; con menos de m_min el evento queda "sin nulo".
- **mu_calor** de las l* ventanas de la alarma; **marea** descriptiva.

## Limites declarados

- **Etapa 7 en espera:** necesita modificar el dataset. **Etapa 8** (cierre y acta) es de Code.
- **Onsager mezcla dos cosas en mundos de Hawkes:** las hijas del mismo signo pasan a continuacion (E1.21), de modo que Delta_tau
  contra Onsager recoge no linealidad y forzamiento posterior al disparo a la vez (control S1(iii): rechaza con n = 0.6, p 0.005).
  Hace falta un nulo con replicas: E1.24 decidida el 2026-11-08 (lineal + Hawkes de dos signos, en la familia de Holm de 6),
  pendiente de construir; directriz para Argon en
  `PRUEBAS_JEFE/lilyth/MENSAJE_ARGON_NULO_CON_REPLICAS.md`.
- **Estado de E1.24 (2026-11-11): construido en `PRUEBAS_JEFE/lilyth/etapa_lh/`, NO conectado a esta base.**
  - Funciona:
    - la inferencia indirecta del estructural (converge en 2-3 iteraciones);
    - la cota de rho, que invierte la ligadura con pendiente 0.66-0.88;
    - pasan R1, R3, R6, R7, R8 y D-ARFIMA.
  - **No pasa el tamano en mundos con replicas:** L+H rechaza 5/24 y Onsager 4/24, con S < 0 en 20 de 24 (la recuperacion real es mas
    rapida que la del nulo).
  - Dos sondas:
    - el doble conteo de replicas en el perfil inyectado queda descartado (perfiles de 1-4 pasos, replicas a >= 6);
    - sin replicas ya hay un sesgo pequeno y real (S medio -0.023, ee 0.0097) con tamano correcto (2/40). Las replicas lo agrandan y
      las hijas de L+H no lo corrigen.
  - **Queda como CURIOSIDAD, por decision del Comandante:** probablemente faltan conceptos por conectar, no un error de codigo.
  - Detalle en la seccion 8.6 del README canonico (`PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md`).
  - Lo que el detector hace con las replicas (realarmas, pendiente, signo contrario invisible) esta en los limites medidos de
    `../metrologia/README_OPERADOR_DELTA.md`.
- **En mundos AR el nulo sin memoria casi nunca alcanza m_min por la via (a):** sale entero de la via (b). Con un observable de
  memoria basal (Hurst ~0.7 en actividad, volumen y flujo; ~0.5 en retornos) en dato real tambien: la ficha reporta la fraccion por
  cada via, y para esos observables el nulo lineal que corresponde es `ARFIMA_WHITTLE` con la d del sensor ELW.
- **dY de una sola observacion** restituye al primer toque por regresion a la media: el tau_rest crudo no es publicable, el exceso
  sobre los nulos si.
- Los eventos de los primeros W_ancla + guarda pasos del ciclo caliente no tienen ancla.
- La bateria de la Etapa 1 contra este orquestador falla K2, K3, K5 y K2b: median replicas que con E1.21 son continuaciones. La
  bateria 3.0 cubre la semantica nueva.

## Verificado (2026-11-07, desde esta base)

Seis baterias (Etapas 2, 3.0, 3, 4, 5, 6) pasan; mutantes de Argon 68/68 y de las correcciones de Code 12/12 cazados; la base entera
compila en una unidad con /W4 sin avisos (`tests/compilar_humo.bat`).
