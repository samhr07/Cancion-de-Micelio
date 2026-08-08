# Adenda A al `ORDEN_TRABAJO_PROPAGADOR_3_1.md` — Sustitución del §3

**Reemplaza íntegramente el §3 de la v3.1** (contraste de difusividad por calibración de
simulador). El resto del documento —§1, §2, §4, §5, §6, §7— queda **sin cambios y vigente**.

**Motivo: el §3 estaba mal especificado por mí.** Code hizo bien en no ajustar el simulador
hasta que cuadrase; el simulador no hay que arreglarlo, hay que **eliminarlo**.

---

## A.0 Auditoría de mi propia simulación

Antes de reemplazar nada, verifiqué la tabla del §3.2. Tres resultados:

| comprobación | resultado |
|---|---|
| Sensibilidad al truncamiento `K` | **Ninguna.** `K = 500`, `3000`, `20000` dan idéntico a 3 decimales. No era el problema |
| Generador de signos | Degrada, aunque menos de lo hallado por Code: γ pedida 0.30 → **0.330** medida, 0.50 → **0.508**, pero 0.70 → **0.552** |
| Sesgo contra la γ **medida** | **+0.113, +0.063, −0.042** — **cambia de signo** |

La tercera fila es la que cierra el asunto. **Un sesgo de tamaño finito no cambia de signo.** El
cruce no está siguiendo a la teoría; sigue a otra cosa.

**El error es conceptual, no numérico.** La relación viene del espectro:

```
Ĝ(f) ~ f^(β−1) ,  S_ε(f) ~ f^(γ−1)   ⟹   S_p(f) = |Ĝ(f)|²·S_ε(f) ~ f^(2β+γ−3)
difusivo exige  S_p(f) ~ f^(−2)      ⟹   β = (1−γ)/2
```

Es una condición **asintótica de baja frecuencia** (`f → 0`, rezagos largos). Yo la puse a
prueba con `ρ₁`, un estadístico de rezago 1 dominado por las frecuencias **altas**. Estadístico
de un régimen para contrastar una condición del otro. La teoría está bien; el contraste que
especifiqué no la mide.

---

## A.1 Lo que la reemplaza: la pendiente de la firma de volatilidad

El §1 ya construyó la máquina correcta. Del mismo desarrollo espectral, el exponente de Hurst
del precio es `H_p = 1 − β − γ/2`, y como `σ(H) ~ H^(H_p)`:

```
  d log[ σ(H)/√H ]                              (1−γ)
  ─────────────────  =  H_p − ½  =  ─────  −  β
     d log H                                      2
```

**La pendiente de la firma de volatilidad ES la desviación respecto de la condición de
difusividad.** Directamente, sobre datos reales, sin simulador y sin umbral calibrado.

| pendiente | `H_p` | lectura |
|---|---|---|
| `= 0` | 0.50 | difusivo — `β = (1−γ)/2` — **sin ventaja explotable** |
| `< 0` | < 0.50 | subdifusivo — `β > (1−γ)/2` — el impacto decae demasiado rápido — **REVERSIÓN** |
| `> 0` | > 0.50 | superdifusivo — `β < (1−γ)/2` — **MOMENTUM** |

**Consecuencia que conviene ver:** para obtener la **dirección** no hace falta estimar ni `β` ni
`γ`. La pendiente sola basta. `β` y `γ` siguen haciendo falta para saber si el propagador es el
**mecanismo**, pero no para el veredicto. Dos niveles, y el que decide es el barato.

Constantes que esto elimina respecto al §3 original: el simulador entero, su `K`, su generador
de signos, su exponente `α`, y el umbral calibrado.

---

## A.2 Lo que los datos del §1 ya dicen

Aplicando lo anterior a la tabla de firma de volatilidad ya medida (3.19 h, ν = 39 tx/s):

| región | rango | pendiente | `H_p` | lectura |
|---|---|---|---|---|
| microestructura | 0.5–10 s | +0.0612 | 0.561 | rebote bid-ask muriéndose |
| **meseta** | 10–30 s | **−0.0022** | **0.498** | **difusivo, limpio** |
| **larga** | **30–300 s** | **−0.1067** | **0.393** | **REVERSIÓN** |
| por encima del ruido | 10–300 s | −0.0731 | 0.427 | reversión |

Y en razón de varianzas contra la meseta: `VR(60 s) = 0.937`, `VR(120 s) = 0.739`,
`VR(300 s) = 0.627`. **Monótona en tres puntos consecutivos.**

Traducido: `β − (1−γ)/2 = +0.107` en la región de 30–300 s. El impacto decaería más rápido de
lo que exige la difusividad.

**Y cae exactamente donde `H* = 50 s`**, que se derivó de comisiones y volatilidad por un camino
completamente independiente.

### A.2.1 Por qué esto NO es todavía un resultado

- Son **38 ventanas no solapadas** a 300 s en 3.19 h.
- Mi error estándar (`z ≈ −2.0` a 300 s, `−2.2` a 120 s) supone **momentos gaussianos**. Con
  curtosis 1179.7 el error real es **bastante mayor**: ese `z = −2` es probablemente `z ≈ −1`.
- La región de microestructura da pendiente **positiva** (+0.061), lo que recuerda que la firma
  tiene estructura propia a escalas cortas y que el rango de ajuste no es inocuo.

Lo que sí tiene peso es la **monotonía en tres puntos** y la **coincidencia con `H*`**. Es la
primera pista cuantitativa a favor de que haya reversión a este horizonte, y la captura de 48 h
la resuelve o la mata.

---

## A.3 Especificación del contraste

### A.3.1 Estadístico

Regresión de `log[σ(H)/√H]` sobre `log H`, con `σ(H)` estimada por **ventanas solapadas** (más
eficiente; la dependencia la absorbe el bootstrap del A.3.3).

### A.3.2 Rangos — declarados ANTES de ajustar

El rango de ajuste no puede elegirse mirando dónde sale la pendiente que gusta. Se fija así:

1. **`H_lo` se determina por la firma misma**, no por criterio: es donde `σ(H)/√H` deja de
   crecer, o sea el final de la región de microestructura. En la captura actual son **10 s**.
   Recalcularlo sobre los datos nuevos y **commitearlo antes** del ajuste.
2. **`H_hi`** = el mayor horizonte con al menos **50 ventanas no solapadas**. Con 48 h eso da
   `H_hi ≈ 3 456 s`; se redondea a la baja a **1 800 s** para dejar margen.
3. Reportar la pendiente en **los tres rangos**: `[H_lo, H_hi]`, `[H_lo, 300 s]` y
   `[300 s, H_hi]`. Si las tres coinciden en signo, el resultado es robusto al rango; si no,
   decirlo y no elegir.

### A.3.3 Errores estándar — bootstrap por bloques, no fórmula

La fórmula de Lo–MacKinlay supone momentos gaussianos. **Con curtosis 1179.7 es inservible**, y
usarla sería el mismo error que los umbrales χ² sobre el NIS: asintótica gaussiana donde la
distribución real no lo es. Tercera aparición del patrón en este proyecto.

- **Bootstrap por bloques móviles**, longitud de bloque ≥ `5 · H_hi` para preservar la memoria
  larga del flujo. Con `H_hi = 1 800 s` son bloques de 9 000 s, o **19 bloques en 48 h** —
  ajustado. Reportar el número de bloques efectivo; si baja de ~15, reducir `H_hi`.
- Intervalo de confianza al 95 % para la pendiente, **por percentiles**, no por `±1.96·EE`.
- **Contraste de la pendiente contra cero**, con su p-valor bootstrap.

### A.3.4 Control obligatorio

Repetir todo el procedimiento sobre **incrementos barajados** de los mismos datos. La pendiente
debe salir estadísticamente indistinguible de cero. Si sale distinta, el estimador tiene sesgo
propio y el resultado sobre datos reales no vale.

Es el control que la v2.2 introdujo y que salvó aquel análisis (`C = 0.783` sobre datos
barajados, donde no hay nada). Aquí es igual de obligatorio.

---

## A.4 Relación con el §2, que sigue vigente

`β` y `γ` **siguen estimándose** según el §2 y el propagador con forzamiento medido, ahora que
`captura_larga.py` persiste `tr_maker`. Pero cambian de papel:

| magnitud | papel nuevo |
|---|---|
| pendiente de la firma | **DECIDE** la dirección — difusivo, reversión o momentum |
| `γ` de la autocorrelación de signos | **explica** — ¿es la memoria del flujo? |
| `β` del propagador | **explica** — ¿es el decaimiento del impacto? |
| coherencia `pendiente ≈ (1−γ)/2 − β` | **valida el mecanismo**: dos vías independientes al mismo número |

Esa última fila es la comprobación más valiosa que gana el documento. Si la pendiente medida
directamente y la calculada desde `β` y `γ` coinciden, el modelo del propagador está confirmado
como mecanismo. Si discrepan, hay algo más y se sabe cuánto.

---

## A.5 Criterios de aceptación que sustituyen a los del §3

- [ ] **Ningún simulador.** Si aparece código que calibra un umbral por simulación, es motivo de
      rechazo del entregable.
- [ ] `H_lo` determinado por la firma y **commiteado antes** del ajuste, con el gráfico.
- [ ] Pendiente en los tres rangos del A.3.2, con acuerdo de signo reportado explícitamente.
- [ ] IC al 95 % **por bootstrap de bloques móviles**, con la longitud de bloque y el número
      efectivo de bloques en el reporte. Fórmula gaussiana: motivo de rechazo.
- [ ] **Control sobre incrementos barajados**: pendiente indistinguible de cero.
- [ ] `β` y `γ` del §2, y la comprobación de coherencia `pendiente ≈ (1−γ)/2 − β` con su residuo.
- [ ] Reproducir la tabla del A.2 sobre los datos nuevos, para poder comparar con la captura
      corta y ver si el hallazgo preliminar sobrevive.

---

## A.6 Lo que esta adenda NO cambia

- **§1** — completo y aceptado. Nota adicional: que `/fapi/v1/commissionRate` no sea legible sin
  credenciales **no invalida el resultado**. Las tarifas VIP 0 públicas son una **cota superior**
  del coste real (descuento BNB y volumen solo lo bajan), luego `H* = 50 s` es una **cota
  superior del horizonte necesario**. El error va en dirección conservadora. Registrarlo así en
  lugar de como criterio incumplido a secas.
- **§2** — vigente sin cambios. Sigue siendo la vía para `β`, `γ` y la forma de `G(τ)`, incluida
  la comprobación de si sus raíces son complejas.
- **§4 (preregistro), §5 (referencia móvil), §6 (régimen congelado), §7 (constantes)** — sin
  cambios. El criterio de abandono del §4.2 sigue en pie y ahora se evalúa con la pendiente en
  lugar de con `β` contra umbral simulado.
- **`Micelio.py`** — sin cambios de modelo, como en toda la serie desde la v2.2.

---

## A.7 Nota de método, para el registro

Van tres veces que este proyecto tropieza con lo mismo: **una distribución asintótica usada como
si fuera la finita**.

1. Umbrales χ² sobre el NIS, con curtosis 1179.7 (v2.1 §6, aún en cola).
2. Razón de verosimilitudes de Harvey en la frontera del espacio de parámetros (v3.0 §4.3).
3. Ésta: `β = (1−γ)/2` contrastada con un estadístico de frecuencia alta, y luego un error
   estándar gaussiano sobre una cola de curtosis 1179.7.

Conviene tratarlo como patrón y no como incidentes. Regla de trabajo propuesta: **todo umbral o
error estándar que provenga de una asintótica se calibra por remuestreo sobre los propios datos,
o no se usa.**
