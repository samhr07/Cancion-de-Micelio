# Nota de metrica — `phi` del flujo de mercado y su derivada `Omega`

Linea nueva abierta el **2026-09-05** a peticion del operador. Modulo:
`flujo_omega.py`, **18/18** controles. `Micelio.py` sin cambios.

> **Esto es una METRICA, no una senal.** No reabre la linea de microestructura,
> que la v4.2 §4 cerro con 78 de 79 celdas negativas y la Adenda C con un
> deficit de 43x a 577x. El bloque predictivo del modulo va marcado
> EXPLORATORIO y con su nulo al lado, por la convencion
> CONTEMPORANEO / PREDICTIVO del proyecto.

---

## 1. La definicion

```
phi(t)   = Q_neto(t) * P(t) / tau_0(t)        [USD/s]
Omega(t) = dphi/dt                            [USD/s^2]
```

| simbolo | que es | unidades | de donde sale |
|---|---|---|---|
| `Q_neto` | flujo neto de volumen de mercado firmado | BTC | `sum eps*q` por casilla, con `eps` de `propagador.signo_transaccion` |
| `P` | precio | USD/BTC | punto medio de cierre de la casilla (respaldo: precio de transaccion) |
| `tau_0` | tiempo de restablecimiento del libro | s | **medido**, tres estimadores (§3) |

Dimensionalmente `phi` es un **caudal de nocional neto por unidad de tiempo de
recuperacion del libro**, y `Omega` su aceleracion.

### ⚠ Colision de nombres, y hay que tenerla presente

Este proyecto **ya tiene** un `Omega` y un `phi`, y los dos estan **retirados**:

- el `Omega` de la Sec. 1.4 del PDF (`[Omega] = BTC/Ticks^2`), que alimenta
  `q_inv(Omega)` y `R_base + kappa*Omega^2`, y que cuelga de `omega_m` —
  refutada en la v2.2, la v3.0 y la v3.2;
- `phi'` = ticks/BTC, retirada en el commit `9b2267e` por ser `1/q_bar`.

**El `phi` y el `Omega` de esta nota no son ninguno de los dos**: otras
unidades, otra definicion, otro proposito. En actas se escriben `phi_F` y
`Omega_F` si hay riesgo de confusion. `flujo_omega.py` **no se importa desde
`Micelio.py`** ni toca `constantes_micelio.py`.

---

## 2. La "EDP" es la regla de la cadena, y se mide EXACTA

Como las tres entradas son variables, la derivada total se abre en tres canales
aditivos:

```
Omega = (P/tau_0)*dQ/dt  +  (Q/tau_0)*dP/dt  -  (Q*P/tau_0^2)*dtau_0/dt
        \____ FLUJO ____/    \___ PRECIO ___/    \______ LIBRO ________/
```

y donde `phi != 0`, en forma logaritmica:

```
Omega/phi = Q'/Q + P'/P - tau_0'/tau_0
```

**La descomposicion no es una aproximacion.** Con diferencias centradas y medias
de los extremos, la regla del producto es una identidad algebraica exacta sobre
la rejilla: `T_FLUJO + T_PRECIO + T_LIBRO = Omega` a `3e-16` relativo
(control 1). Importa porque lo que se va a leer es el **reparto de varianza**
entre los tres canales, y un reparto que no suma al total no es un reparto.

El detalle que lo hace exacto: para `f = g*h`,

```
g[k+1]h[k+1] - g[k-1]h[k-1] = mean(g)*diff(h) + mean(h)*diff(g)
```

con `mean` y `diff` sobre los extremos `k-1` y `k+1`. Aplicada dos veces reparte
la diferencia en tres sumandos sin residuo. Ojo: la media del producto `Y*Z` **no
es** el producto de las medias; usar esta ultima rompe la identidad.

### Las dos formas

El operador enuncio primero `phi = Q*tau_0/P` y despues `phi = Q*P/tau_0`. Se
implementan **las dos** (`--forma=P/tau`, por omision, y `--forma=tau/P`) porque
la unica diferencia observable entre ellas es el **signo del canal LIBRO**, que
es justo lo que la descomposicion mide (control 4). Cambiar de forma cuesta una
bandera; leerlas al reves, no.

---

## 3. `tau_0` — el punto delicado, y no se toma del ajuste del nucleo

⚠ **El `tau0` de `migracion_v32.py` NO sirve aqui**, por dos razones
independientes:

1. Esta **retractado**. Sale pegado a su cota inferior y la v4.1 §3.4 lo
   degrado a *"ajuste no identificado en el regimen `tau_0` en cota"*; el §5.2
   de la v3.2 ya decia que `beta` y `tau_0` no estan identificados por separado.
2. Es **un numero por ajuste**, y `phi(t)` necesita una **serie** `tau_0(t)`.

Se mide directamente del libro, con tres estimadores, y se reporta la
concordancia entre ellos en vez de elegir uno a ciegas — la leccion de la v4.1
§3.1, donde `gamma` recorria de +0.69 a −0.26 segun la banda y por eso "no
existe la gamma de este mercado".

| estimador | definicion | parametro libre | papel |
|---|---|---|---|
| **`tau_agot`** | profundidad L1 / caudal negociado [BTC / (BTC/s)] | **ninguno** | **PRIMARIO** |
| `tau_recup(theta)` | semi-recuperacion tras una caida de profundidad `>= theta`, dividida por `ln 2` | `theta`, **barrido** | literal |
| `tau_upd` | 1 / tasa de actualizacion del libro | ninguno | diagnostico |

### Por que SEMI-recuperacion y no recuperacion total

La primera version media el tiempo hasta volver al nivel de profundidad
**exacto** previo. Es la lectura literal de "restablecimiento" y esta mal:
bajo relajacion exponencial hacia un equilibrio, la profundidad se **acerca** a
su nivel previo pero no lo alcanza en tiempo finito, asi que el retorno al nivel
exacto es un **tiempo de primer paso que lo dispara el ruido de cotizacion**, no
la reposicion. Medido sobre la captura sintetica de verificacion, esa version
daba **0.8 s** con `theta = 0.30` — parpadeo — contra una constante verdadera de
3.4 a 24 s.

La semi-recuperacion vale `tau*ln 2` **exactamente** bajo relajacion
exponencial, asi que tiene verdad conocida y se controla: el control 7 recupera
una constante de 5.0 s con **0.41 %** de error.

### ⚠ Y `tau_recup` es DISPERSO

Solo esta definido en las casillas que tuvieron un evento de caida `>= theta`.
Con `theta = 0.50` sobre la captura de verificacion quedan **314 casillas de
17 041, en 226 rachas**: como entrada de `phi` fragmenta la serie y la derivada
casi no tiene donde apoyarse. Es la segunda razon —tras el parametro libre— para
que el primario sea `tau_agot`, que esta definido en todas las casillas validas.

`theta` se **barre** (0.30 / 0.50 / 0.70 / 0.90) y se reporta como curva. Este
proyecto lleva **cinco** umbrales puestos a ojo que hubo que corregir midiendo
(`UMBRAL_TAYLOR_JACOBIANO`, la compuerta de `R^2` de la v3.3, el `|d| < 0.06`
del §3, la fraccion de ceros del §7.3 y la razon fija del §4); este no es el
sexto.

### ⚠ Limitacion heredada, y no se arregla con analisis

`@bookTicker` da **solo el nivel 1**. La reposicion en niveles mas profundos es
invisible, asi que `tau_0` medido aqui es el de la **cola visible**, no el del
libro entero. Ya hay dos medidas independientes de que la liquidez que gobierna
no vive en L1: H1 (2026-08-23, parcial +0.0425 contra techo de nulo +0.0418) y
`lambda = 0.12` en `cont2014.py`. **Este es el tercer argumento para ingerir
`@depth`**, y sigue siendo una decision de captura, no de analisis.

---

## 4. Como se corre

```
python flujo_omega.py --autotest                    18 controles con verdad conocida
python flujo_omega.py --etapa=serie                 construye y cachea la rejilla
python flujo_omega.py --etapa=dia                   los Omega por dia y por bloque
python flujo_omega.py --etapa=relacion              relacion entre variables
```

Banderas: `--fuente=estacional|v33`, `--tau=tau_agot|tau_upd|tau_recup_0.50`,
`--forma=P/tau|tau/P`, `--dt=10` (casilla de `phi`), `--bloque=3600` (bloque de
reporte; 3600 s → **24 Omega por dia**).

Salidas: `telemetria/rejilla_omega_<fuente>.npz`,
`telemetria/omega_bloques_<fuente>.csv` (una fila por bloque, 20 columnas) y
`telemetria/omega_resumen_<fuente>.json`.

### Memoria

Una sola pasada por parte, sobre una rejilla global. Con `dt = 10 s`, tres
semanas son 181 k casillas y ~30 MB: lo que no cabe es el dato crudo (223 M de
filas de libro), no la rejilla. Nunca hay mas de una parte en memoria — la misma
disciplina que `curvas_estacional.limites_tramos` tuvo que adoptar tras morir
tres procesos.

### Huecos

La captura tiene huecos (el de 24.26 h del 2026-08-21, entre otros) y casillas
sin transacciones. Una diferencia centrada que cruce un hueco **no es una
derivada**, es la pendiente entre dos regimenes distintos. La rejilla se parte
en **rachas** de casillas validas consecutivas y se deriva dentro de cada racha
(control 10).

---

## 5. Dos cosas que hay que saber antes de leer la primera tabla

### 5.1 El canal PRECIO va a salir ~0, y es estructural

`dP/P` por casilla de 10 s es del orden de `1e-5`; `dQ/Q` es de orden 1 porque
`Q_neto` cambia de signo continuamente. El canal PRECIO entra en `Var(Omega)`
con peso **~1e-7**. No es un fallo del estimador: es que a esta escala el precio
es, comparado con el flujo, una constante. Si se quiere que el canal PRECIO
pese, hay que subir `--bloque` y `--dt` — y entonces el que se degrada es
`tau_0`.

### 5.2 El nulo por rotacion se satura si no se quita el ciclo diurno

Casi todas estas variables (`nu`, `tau_0`, `prof`, `rv`) llevan el mismo ciclo
de 24 h, medido el 2026-08-23: `nu` recorre **3.9x** entre UTC 13-15 y UTC 4.
Dos series con el mismo periodo correlacionan a **+0.99 por tener el mismo
periodo**, y una rotacion circular las vuelve a alinear en la siguiente vuelta —
asi que el suelo del nulo sale **tambien** en 0.99 y la razon en 1.00. Medido en
la corrida de verificacion, columna CRUDO.

Por eso `informe_relacion` da **dos columnas**: CRUDO y RESIDUO (sin el nivel
del dia ni la forma intradia, `quitar_ciclo`). **Solo se lee la columna
RESIDUO**, y solo si `razon >= 3`. Aquel dia fallaron **tres** nulos propios por
exactamente esto, y la regla que quedo escrita es: separar nivel diario y forma
intradia **antes** de construir cualquier nulo. Ademas el nivel del **dia** es la
palanca grande (`log sigma` recorre 13.75x entre dias), asi que retirarlo no es
quitar una molestia, es quitar el efecto dominante para poder ver el resto.

El modulo **avisa** cuando hay menos de 4 bloques por parametro de ciclo: con
pocos dias la columna RESIDUO esta sobre-restada. Con `--bloque=3600`, 24
casillas horarias y `D` dias hacen falta `>= 4*(D+24)` bloques, o sea unos
**5 dias** para empezar a leerla y las **3 semanas** de `captura_estacional`
para leerla comoda.

---

## 6. Estado de la verificacion

**Lo que esta verificado (18/18, `python flujo_omega.py --autotest`):**

| # | control | resultado |
|---|---|---|
| 1 | regla de la cadena exacta, las dos formas | residuo `2.9e-16` relativo |
| 2 | `Omega` recupera una derivada analitica | error `8.4e-06` con `dt = 0.01 s` |
| 2b | cada canal por separado, no solo la suma | FLUJO `8.4e-06`, LIBRO `6.1e-07` |
| 3 | `tau_0` constante → canal LIBRO exactamente 0 | `0.00e+00` |
| 3b | `phi` constante → `Omega` exactamente 0 | exacto |
| 4 | las dos formas dan el canal LIBRO con signo opuesto | `-1.03` contra `+2.6e-09` |
| 5 | el reparto de `Var(Omega)` suma 1 | `1.000000000000` |
| 6 | busqueda por bloques == fuerza bruta `O(n^2)` | 400 consultas |
| 6b | el resultado no depende del tamano de bloque | identico con `B = 7` y `B = 64` |
| 7 | `tau_recup` recupera una constante de relajacion conocida | 5.0206 s contra 5.0 (**0.41 %**) |
| 7b | evento sin recuperar → CENSURADO, no rapido | exacto |
| 8 | `tau_agot` = profundidad / caudal | exacto a `1e-9` |
| 9 | flujo neto firmado por casilla | exacto |
| 10 | un hueco parte la serie en vez de derivar a traves de el | 2 rachas |
| 11 | la rotacion rompe el emparejamiento y conserva la marginal | `0.999 -> 0.363` |
| 12 | Spearman con empates promediados | `-1` y `+1` exactos |
| 13 | censurados = escalar acotado por el numero de filas | exacto |

**Lo que NO esta hecho: no hay ni una cifra de mercado real.** El contenedor de
esta sesion no tiene las capturas — `telemetria/` esta en `.gitignore` y los
datos viven en la maquina del operador. La ruta de ingesta (parquet → indice →
rejilla → `Omega` → bloques → CSV) se verifico de punta a punta contra una
**captura sintetica** en el formato real, con ciclo diurno, memoria de flujo,
profundidad que baja cuando sube `nu` y un hueco de 40 min. **Ninguna cifra de
esa corrida es una medicion** y ninguna aparece en esta nota como tal.

Para tener numeros de verdad basta correr las tres etapas en la maquina de la
captura, en ese orden.

### Dos defectos propios encontrados en la verificacion

1. **Los censurados se sumaban a todo el array de casillas.** `+= cens` sobre un
   array de `nb` posiciones daba a cada casilla el total, y el recuento final
   salia multiplicado por `nb`: **41 869 440 censurados sobre 1.7 M de filas de
   libro**. Lo delato que el numero era imposible, no una asercion — por eso
   ahora hay una (control 13).
2. **El reparto de varianza sumaba 1.00025.** `np.cov` normaliza con `ddof=1` y
   `np.var` con `ddof=0`; mezclarlas da `1 + 1/(n-1)`, que con `n = 4000` son
   `2.5e-4` — invisible a ojo. Lo cazo el control 5 porque pide la suma a
   `1e-10`. Es el argumento para poner las tolerancias donde el estimador puede
   llegar y no donde uno se conforma.

---

## 7. Lo que falta, en orden

1. **Correr las tres etapas sobre `captura_estacional`** y sobre `captura_v33`.
   Es lo unico que convierte esto en mediciones.
2. **Decidir `dt` y `--bloque` midiendo**, no a ojo. `dt = 10 s` viene de
   `cont2014.py` y `--bloque = 3600` de querer 24 `Omega` por dia; ninguno esta
   calibrado. El barrido natural es `dt` contra la estabilidad de `tau_agot` y
   contra el peso del canal PRECIO (§5.1).
3. **Comprobar si `tau_0` es una cantidad bien definida**, que es la pregunta que
   decide si la metrica se sostiene: si `tau_agot` y `tau_recup(theta)` no se
   mueven juntas sobre dato real, `tau_0` no esta bien definida y hay que decirlo
   antes de construir nada encima. La tabla de concordancia de `--etapa=serie`
   es exactamente ese contraste.
4. **`@depth`**, si se quiere un `tau_0` del libro y no de la cola visible (§3).
