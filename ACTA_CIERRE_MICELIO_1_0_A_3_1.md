# Acta de cierre — Canción del Micelio, v1.0 a v3.1

**Fecha:** 2026-08-29. **Autor:** Samuel Hoyos R.
**Alcance:** cierra el proyecto que va del PDF `Trading_Bot.pdf` (Julio 2026) hasta la
sesión 2026-08-08 (d), última registrada en `CLAUDE.md` de este repositorio.

Este documento no continúa el proyecto: lo liquida. Deja escrito qué quedó
**establecido**, qué quedó **refutado**, qué **sobrevive como infraestructura** y qué
queda **abierto y sin dueño**. El proyecto nuevo (cribado de activos, v4.0) arranca con
preregistro y orden de trabajo propios y **no hereda ninguna hipótesis** de este.

---

## 1. Lo que el proyecto se propuso

El PDF postula que el precio tiene un **modo oscilatorio endógeno** `ω_m`, extraíble por
EMD → Hilbert, del que cuelga toda la Sec. 1.4: los acoplamientos `Φ`, `Ψ`, `Ω`, el umbral
`Ω_crit`, los nodos de fase, la matriz de estado armónica `A_arm` y la modulación del filtro
de Kalman por `ρ_k = 1 + γ_ω|ω_m| + γ_Q|ΣQ|`. La tesis central del sistema era que ese modo
existe, se mide y **modula el filtro**.

## 2. Veredicto: el modo oscilatorio no existe

Cuatro líneas de evidencia independientes, todas convergentes, ninguna refutada:

| vía | sesión | resultado |
|---|---|---|
| Sustitutos por barajado (KS sobre la distribución de `f_hz`) | v2.2 | **p = 0.98** — REAL y BARAJADO indistinguibles; medianas 112.0 s / 112.7 s |
| Multitaper (Thomson) contra nulo AR(1) | v2.2 | 6.6 % / 7.4 % de frecuencias sobre el umbral 95 % (azar ≈ 5 %) — **sin exceso** |
| AR(2) sobre 446 892 ticks reales | v3.0 | `k = +4.7e-07 ± 1.7e-04`, **H₀: k = 0 no se rechaza (p = 0.6175)**; `m < 0`, raíces reales |
| Ciclo estocástico de Harvey | v3.0 | `φ₂ = +0.216 > 0` ⟹ `√(−φ₂)` no es real: **sin parametrización válida** |

Y la demostración que lo cierra, fijada como test permanente
(`test_v22_emd_no_tiene_hipotesis_nula`): sobre un **paseo aleatorio puro** la cadena
EMD → Hilbert devuelve un período de 118.1 s con `C = 0.690` en ventana de 384 muestras y
37.5 s en media ventana. **El período escala con la ventana, no con la señal.**
`ω_m` era salida del algoritmo, no del mercado.

La explicación mecánica quedó exacta a 7 decimales: todo el AR(2) del mercado real es
**paseo aleatorio + rebote bid-ask**. Con `x_t = x_{t−1} + r_t` y `r_t = a·r_{t−1} + ε`
sale `φ₁ = 1+a`, `φ₂ = −a`, luego `k = 0` **por construcción**. Medido:
`a = ρ₁(retornos) = −0.216061`, `φ₁` predicho/medido +0.783939 / +0.783939 (error 2.0e-07).
No queda residuo que atribuir a una fuerza recuperadora.

**Muere con `ω_m`:** `A_arm`, el término `γ_ω·ω_m` de `ρ_k`, `c²_vol`, los nodos de fase,
`S_ref` como ancla de la condición terminal de `U`, y la Sec. 1.4 entera
(`Φ`, `Ψ`, `Ω`, `Ω_crit`, `A_arm`, `κΩ²`).

## 3. Lo que sobrevive, y es mucho

Nada de lo siguiente depende de `ω_m`:

- **Reloj de transacciones (v2.0).** El sistema descartaba el **96 %** de los datos que ya
  recibía por tener `P_spot` escalar. Con Δn = 1 pasa de 0.76 a 24–94 observaciones/s y el
  100 % de los ticks lleva medición. `Q(Δt)` con `Q_N = Σ AⁱQAⁱᵀ` (no `N·Q`), ingesta por
  lotes, deduplicación por id de trade, detección de huecos.
- **Capa de riesgo de cuenta (v1.3).** Tri-estado `MODO`, 7 guardas con `causa_halt`
  distinguible, ruta de cierre que jamás reporta éxito sin posición plana confirmada,
  máquina de episodios con `DETENIDO` terminal, `verificar_instancia_unica`.
- **Loeper y el NMPC.** Necesitan `λ` y `Γ`, no `ω`. La inversión temporal de la EDP
  (backward desde condición terminal) sigue siendo correcta y verificada por contraste.
- **`R_n`** con alta probabilidad: estimar una tendencia es mucho más robusto que extraer
  un ciclo.
- **La disciplina de medición**, que es el activo más valioso que deja el proyecto: nulo
  por barajado, preregistro commiteado antes de medir, control obligatorio del estimador
  contra sí mismo, y la regla de no ajustar el instrumento hasta que dé el resultado
  esperado.

## 4. Lo que el proyecto aprendió a costa de sesiones enteras

Se registra porque el patrón se repitió y volverá a repetirse:

1. **La misma medición entrando N veces al filtro.** Apareció tres veces con tres disfraces
   (90 correcciones por paquete en v1.3; `R_n` retenido 11.6× en v2.0; ventanas REST
   solapadas sin deduplicar). Siempre se manifiesta como autocorrelación de la innovación y
   siempre se lee mal como error de modelo.
2. **Un estimador sin control es un generador de hallazgos.** La pendiente de la firma de
   volatilidad (Adenda A) daba −0.153 y parecía reversión; sobre **incrementos barajados**,
   donde no hay estructura alguna, daba −0.040. El sesgo era del estimador. Segunda vez que
   el control barajado salva un análisis.
3. **Un estadístico puede ser degenerado bajo su propio nulo.** La razón `R(final)/R(pico)`
   explota cuando el pico es ruido; el `p = 0.917` que devolvía no significaba nada.
4. **Los datos sucios enmascaran, no solo ensucian.** El feed emite trades con `p = 0`
   (~0.2 %); al limpiarlos la curtosis **sube** de 601.8 a 1179.7 y el |incremento| máximo
   cae de 65 245 a 11.8 USD. Los ceros inflaban σ⁴ y aplastaban el cociente.
5. **Un test que falla por su propia culpa cuesta lo mismo que uno real.** Al menos siete
   casos documentados (evaluación de `d[k] = v`, dos ω distintos para medir un salto entre
   ramas, carga de máquina en una medición de coste intrínseco, `mode="same"` metiendo el
   futuro en el presente).

## 5. Resultados que el proyecto nuevo hereda como MEDICIONES, no como hipótesis

Son números, no creencias, y entran al proyecto nuevo como constantes de diseño:

| magnitud | valor | sesión |
|---|---|---|
| `ρ₁` de retornos (rebote bid-ask) | −0.216061 | v3.0 |
| `s_eff` de Roll en BTCUSDT | 0.2062 USD/BTC ≈ 2 ticks | v3.1 §1 |
| `H*` maker+maker | **50.0 s** | v3.1 §1 |
| `H*` maker+taker | ~200 s | v3.1 §1 |
| `H*` taker+taker | **fuera del rango medido (300 s)** | v3.1 §1 |
| meseta de `σ(H)/√H` | H = 10–30 s | v3.1 §1 |
| latencia contra `H*` | 300 ms / 50 s = **0.6 %** | v3.1 §1 |
| propagador con forzamiento medido `R(τ)` | pico 1.95 USD/BTC a τ = 398 ticks (71 s) | v3.1 §2 |
| control de signos barajados | pico se desploma a 0.003–0.016 (**factor 150–700×**) | v3.1 §2 |

**La consecuencia de diseño más importante del proyecto viejo:** con comisiones taker, `σ`
no alcanza el coste dentro de los 300 s medidos. Una ida y vuelta taker **no se paga** al
régimen de volatilidad de BTC. Eso, y no `ω_m`, es lo que empuja al proyecto nuevo.

## 6. Lo que queda abierto y SIN DUEÑO

El proyecto nuevo no los adopta. Se listan para que nadie los dé por cerrados:

- La hipótesis **(B)** de la v2.2 —ciclo en decenas de minutos— **no queda probada ni
  refutada**. Requiere `W ≥ 10 800` contra los 384 usados; la captura de 48 h sufrió un
  corte de DNS de 10 h y quedó partida.
- **Separar (A) de (A′)** —"no hay ciclo" contra "hay memoria, no ciclo"— quedó pendiente:
  el brazo IAAFT no era fiable aplicado al nivel de precio.
- **La forma del propagador** (monótona o con sobrepaso) sin decidir: hace falta un
  estadístico con nulo no degenerado.
- **La discrepancia de signo entre firma solapada y no solapada** (Adenda A) sin resolver.
  Sospecha: la solapada sobrepondera los tramos de alta actividad.
- **Fase 2 (ALS)** nunca se corrió. `r_S,base` y `r_EMD` siguen medidos mal por 8× y 161×.
- **Fase 3 (`Ω_crit`)** y las 30 corridas de Testnet: sin ejecutar, y ahora sin sentido,
  porque `Ω_crit` muere con `ω_m`.
- El silencio de `@aggTrade` en `fstream.binance.com` sigue **sin explicación**.
- El escalón de comisiones nunca se leyó de la cuenta (`/fapi/v1/commissionRate` es firmado
  y el Modo LECTURA no tiene credenciales). Criterio de aceptación del §8 de la v3.1
  **NO CUMPLIDO** y declarado como tal.

## 7. Estado del código al cerrar

`Micelio.py` queda **como está** y no se toca. Los módulos de análisis
(`oscilador.py`, `propagador.py`, `experimento_v22.py`, `experimento_v30.py`) no se importan
desde el orquestador y quedan como registro reproducible de las mediciones citadas arriba.
Suite de aceptación: **56/56** (`python tests_v13.py`).

---

**Cerrado.** Lo que siga vive en `PREREGISTRO_CRIBADO_4_0.md` y
`ORDEN_TRABAJO_CRIBADO_4_0.md`, con hipótesis, umbrales y regla de parada propios.
