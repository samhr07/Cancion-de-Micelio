# Operador de coordinacion de pares (`ark/coordinacion/`)

Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-05, diseno acordado con el Comandante. Prerregistro y pruebas:
`JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/coordinacion_sensor/` (`EXPECTATIVAS_COORDINACION.md`, `MEDIDO_COORDINACION.md`).

**No es el modulo heredado.** El de `ark_cxx/include/ark/coordinacion/` (KSG con copula de rangos exactos, reticulo PhiID de 16
atomos, regularizador dual, CUSUM con h de Siegmund, compuerta Q, constantes de EL ABISMO) queda NO AUDITADO y solo lo usa el
ARK 1.0. Su auditoria (11-03) midio: alarma el 53-57 % del reposo, Siegmund 38 veces corto, compuerta Q = moneda, la copula sesga
+0.036 nats con datos continuos.

## Que hace

```
retornos de X e Y (barras alineadas, reloj de pared)
  -> SensorCoordinacion<W>   theta^ cada S barras sobre las ultimas W barras (todas las escalas y observables)
  -> ObservadorSaltos        uno por componente de theta^, calibrado sobre su nulo (operador delta: nivel y temperatura)
  +  ObservadorSaltos        sobre el crudo de X y el de Y (para la etiqueta)
  =  ObservadorCoordinacion<W>   un solo presupuesto de falsas alarmas, etiqueta y guarda de pares en cada evento
```

| cabecera | contenido |
|---|---|
| `ksg_directo.hpp` | KSG-1 (MI), Frenzel-Pompe (informacion condicional), PIT aleatorizado, splitmix64 |
| `sensor_coordinacion.hpp` | escalas admisibles, observables, theta^ por ventana, nulos SIGNOS de centro cero y ROTACION |
| `observador_coordinacion.hpp` | calibracion sobre los nulos, reparto del ARL0, etiqueta, ciclo caliente |

## Decisiones (todas del Comandante, 2026-11-05)

- **Cantidades aparte, cada una con su estimador directo** (4 busquedas de vecinos por observable y escala):
  TE(A->B) = I(A_p; B_f | B_p) por Frenzel-Pompe; MI(A->B) = I(A_p; B_f) por KSG-1; SIN = TE - MI (sinergia neta). Los 16 atomos no.
- **PIT aleatorizado**, sin recorte en 0.
- **Pasado de una barra por escala.** Escalas s = 2^j que dividen W, con s <= W / max_ratio (el del modulo de ventanas) y n_j = W/s
  >= k + 3. **W y k son entradas**; el conjunto de escalas admisibles sale de ellos.
- **Dos nulos**: DIRECCION con SIGNOS de centro cero en los retornos base (conserva el RV de cada bloque y su co-movimiento);
  MAGNITUD con ROTACION de Y n/2 (conserva la memoria de cada serie, rompe el emparejamiento). DIRECCION_NORMALIZADA (R_j/sqrt(RV_j))
  es una herramienta opcional para el confundido de la volatilidad sincronizada.
- **El operador delta es el unico detector** (sin regularizador propio). Solo pares. Reloj de pared. La guarda de pares
  (EJECUTABLE, REZAGO, NO_VERIFICADO) va sellada en cada evento; fuera de EJECUTABLE no se lee como transmision.
- La informacion mutua del mismo instante I(X_t; Y_t) queda para un modulo de reporte posterior.

## Uso minimo

```cpp
#include "ark/coordinacion/observador_coordinacion.hpp"
using Obs = ark::coordinacion::ObservadorCoordinacion<240>;   // W = 240
auto ob = std::make_unique<Obs>();                              // grande: siempre en el monticulo
Obs::Configuracion c{};
c.arl0_familia_barras = 20000;          // una falsa alarma de TODA la familia cada 20 000 barras
c.sensor.paso_sensor = 8;               // S
c.sensor.k_vecinos = 4;                 // k
c.clase_par = ark::ventanas::ClasePar::EJECUTABLE;
ob->calibrar(c, rx_historia, ry_historia, okx, oky);
for (...) { auto s = ob->paso(rx, ry); for (std::size_t i = 0; i < s.n_eventos; ++i) usar(s.eventos[i]); }
```

## Dos presupuestos (`Configuracion.presupuesto`)

| | `BONFERRONI_CANALES` (omision) | `FAMILIA_EPISODIOS` (fase 2 de la calibracion de h) |
|---|---|---|
| crudos de X e Y | dentro del presupuesto, con su parte | solo etiqueta: cada uno con ARL0 = ARL0_familia, fuera del presupuesto |
| unidad de falsa alarma | alarma por canal | **episodio** de coalescencia por subfamilia: abre cuando un componente abre, cierra cuando todos cerraron |
| reparto | ARL0 x 2 tot por canal, extrapolado canal a canal | dos subfamilias por nulo, S1 = DIRECCION + NORMALIZADA (SIGNOS) y S2 = MAGNITUD (ROTACION), cada una con ARL0_sub = 2 ARL0_familia; igual rareza A por canal dentro de la subfamilia |
| como se fija el umbral | `buscar_h_con` por canal (recta en log ARL fuera de lo legible) | **se mide**: la maquina de estados de cada componente corre sobre B realizaciones del sustituto juntadas como tramos y se biseca log A hasta tasa de episodios <= 1/ARL0_sub |
| B y certificacion | -- | m_min = ceil((z_{1-alfa/2}/ln 2)^2) = 8; B = min(B_max, ceil(m_min ARL0_sub / N)); si B queda topado, el objetivo es `arl0_max_certificable` = B N / m_min y la ficha lo marca (`acotado_por_historia`) |
| salida extra | `episodio_abierto`, `episodio_nuevo` por subfamilia | lo mismo, y `ficha_subfamilia(q)` |

Lo generico (m_min, B, episodios, busqueda de A*) vive desde el 2026-11-05 en `metrologia/calibrador_familia.hpp`, compartido con
los demas operadores; las cifras de este modo no cambiaron (control G1). Reglas aprobadas por el Comandante sobre el dictamen del auditor externo Argon (`mensajes/RESPUESTA_ARGON_A_REVISION_CODE_CALIBRACION_H.md`).
En `FAMILIA_EPISODIOS`, `recalibrar_arl0` lanza `logic_error`: B depende del ARL0, se calibra de nuevo. `Configuracion.hilos` (0 =
todos) reparte la calibracion; el ciclo caliente sigue sin reservas. Controles Q1-Q8 y mutantes:
`PRUEBAS_JEFE/calibracion_h/` (`EXPECTATIVAS_FASE2.md`, `fase2.cpp`, `mutaciones_fase2.py`).

## Limites declarados

- Solo el camino de serie (barras ya alineadas). La entrada por filas y el motor multiactivo vienen despues.
- En BONFERRONI_CANALES, con muchos componentes el ARL0 por canal queda muy por encima de la historia y el h de cada componente sale
  extrapolado (log-lineal, declarado en su ficha). FAMILIA_EPISODIOS no depende de esa extrapolacion para la tasa de la familia.
- En FAMILIA_EPISODIOS cada subfamilia queda fuera de muestra en x1.6-1.8 de su objetivo: esta dentro del factor 2 de m_min = 8 y lo
  explica el sesgo de elegir el A minimo. Para tener margen se calibraria a la cota superior del intervalo (decision abierta).
- Coste: O(n^2) por busqueda de vecinos; con W = 240 y los tres observables, unos milisegundos por theta^.

## Lo medido (2026-11-05, sinteticos; acta `PRUEBAS_JEFE/coordinacion_sensor/MEDIDO_COORDINACION.md`)

- Controles de identidad 18/18, a la primera: PIT = referencia; MI = KSG heredado bit a bit donde aquel no recorta; CMI = referencia
  ingenua; verdad gaussiana conocida (MI, TE y sinergia nula); sensor = bloques agregados a mano; nulos (RV exacto con SIGNOS, rotacion y
  su mascara); cada componente = ObservadorSaltos solo con su parte del ARL0; reparto; etiqueta; causalidad; cero reservas. Mutantes 19/19.
- **Calma, con BONFERRONI_CANALES: FALLA.** La familia (90 componentes, ARL0 20 000) emite 5.7 veces su limite de falsas alarmas.
  Causa: la extrapolacion del h del operador delta (seccion 5 de `../ventanas/README_VENTANAS_Y_DETECTORES.md`), no este sensor ni sus
  nulos: dentro de muestra da lo mismo y el crudo solo se pasa x14.
- **Calma, con FAMILIA_EPISODIOS: PASA** (fase 3, acta `PRUEBAS_JEFE/calibracion_h/MEDIDO_FASE2.md`):
  - 8.2e-5 episodios por barra contra el limite 1.0e-4; las 80 subfamilias certificadas (B 5, E* 8);
  - P3-P7 pasan, coste x6.9;
  - P2 falla por una replica: DIRECCION TE X->Y 6/8, porque los umbrales mas altos retrasan la deteccion de 1 900 a 7 750 barras.

  **Para leer eventos como detecciones, usar este modo.**
- Como contraste (con el mismo sesgo), cada fenomeno lo detecta el grupo que corresponde, en 8/8 y mucho antes que en calma: acoplamiento
  lineal -> TE X->Y de DIRECCION (~1 900 barras, escala nativa); volatilidad sincronizada -> MAGNITUD (~1 600); acoplamiento de producto
  -> SIN X->Y (~3 500). La DIRECCION cruda se contamina con la volatilidad sincronizada y la NORMALIZADA no. Con la varianza x4, el 97 % de
  los eventos sale COORDINACION_CON_VOLATILIDAD.
- Coste medido: W 240, S 8, tres observables: unos 6 minutos de calibracion por par con 2^16 barras (tres series de theta^).
