# Preregistro — ORDEN_TRABAJO_PROPAGADOR_3_1

**Fecha de escritura: 2026-08-08.** Este documento se commitea **antes** de ejecutar
cualquier medición del §1, §2 o §3 de la v3.1. Su hash se referencia en el reporte final.

> §4 de la orden: "Este proyecto ha anulado dos veredictos por emitirlos sobre datos que no
> los sostenían. Los criterios se escriben y se commitean con fecha antes de ejecutar nada."

---

## 1. Criterio de ÉXITO (§4.1) — las cuatro conjuntamente

1. **`R²` fuera de muestra > 0** en corte temporal, con **banda de embargo ≥ `H*`** entre
   entrenamiento y prueba.
2. **`E[Δp | información en t]` al horizonte `H*` supera `c(u)` con margen**, no lo iguala.
   Se declara "con margen" como **≥ 1.5 × c(u)**, fijado aquí para que no se ajuste después.
3. **`β` fuera del intervalo simulado** del §3.2, con dirección declarada:
   - `β < ` intervalo ⇒ momentum
   - `β > ` intervalo ⇒ reversión
   Se usa el intervalo **al 95 %** de la distribución simulada.
4. **Residuo del propagador sin autocorrelación remanente**: Ljung-Box sobre el residuo a
   rezagos largos con `|ρ₁| < 0.05` y decaimiento a cero.

## 2. Criterio de ABANDONO (§4.2)

Si sobre datos limpios, con forzamiento medido, referencia móvil (§5) y régimen congelado
(§6), al horizonte `H*` derivado se cumplen **las tres**:

- `β` cae **dentro** del intervalo simulado al 95 %, **y**
- `E[Δp]` no supera `c(u)` en **ningún** régimen declarado, **y**
- el residuo es blanco (`|ρ₁| < 0.05`)

entonces **la hipótesis de estructura explotable a este horizonte queda abandonada**. El
proyecto pasa a lo que la v3.0 §6.2 lista como superviviente: reloj de ticks, ingesta, capa de
riesgo, Loeper y NMPC como ejecutor, **sin capa predictiva**.

> "`ω_m` lleva tres refutaciones de rigor creciente y cada una se respondió con 'quizá a otra
> escala'. Ésta es la escala donde está el dinero, y no hay otra a la que retirarse."

## 3. Regresor de RÉGIMEN (§6) — declarado antes de mirar

- **Primario: volatilidad realizada** por bloques, con núcleo realizado. Se elige por la razón
  que da el §6 —ya hace falta para `H*`, así que no añade maquinaria— y no por su desempeño.
- **Secundario 1: magnitud del desbalance de flujo** (|suma de ε·f(v)| por bloque).
- **Secundario 2: sesión** (Asia / Europa / América), exógena y sin riesgo de circularidad.

Los cortes de régimen se fijan por **terciles de la distribución del regresor primario sobre
todo el registro**, no por inspección. No se cambian después.

## 4. Convención de signo (§2.2) — predicción declarada

Campo `m` de Binance = "el comprador fue el maker". Por tanto:

```
m = True   -> taker vendía  -> eps = -1
m = False  -> taker compraba -> eps = +1
```

**Predicción falsable: `G(0) > 0`.** Si sale negativo, el signo está invertido y **se para**;
no se le da la vuelta al signo y se sigue como si nada.

## 5. Partición temporal

- Corte temporal simple: primer 70 % entrenamiento, último 30 % prueba.
- **Banda de embargo entre ambos: `≥ H*`**, verificada explícitamente.
- Sin validación cruzada aleatoria: mezclaría futuro con pasado.

## 6. Lo que NO se decide con esta tanda

- Integrar el propagador en el filtro (§9 de la orden).
- Retirar la Sec. 1.4 ni la Sec. 2 del PDF.
- La hipótesis (B) de escala larga.

## 7. Limitaciones conocidas ANTES de medir

Se anotan ahora para que no se presenten después como matices:

- **El escalón de comisiones no se puede leer de la cuenta**: el endpoint es firmado y el Modo
  LECTURA no tiene credenciales. Se usarán las tarifas públicas VIP 0 de USDⓈ-M **marcadas
  como asumidas**, y el criterio §8 correspondiente quedará **no cumplido**, no fingido.
- La captura larga sufrió un corte de DNS de 10 h; el tramo continuo utilizable es de ~3.2 h.
- El nulo espectral de la v2.2 se midió sobre datos sucios y **no se cita como cerrado** hasta
  rehacerse.
