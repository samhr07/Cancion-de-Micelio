# Pruebas del Operador de Lilyth — Correcciones Piloto

Banco de pruebas aislado en `ARQUITECTO DE DATOS/Pruebas del Operador de Lilyth_Correcciones Piloto/` sobre una copia local de cabeceras (`inc/`).  
**Regla de seguridad:** `ark_cxx/base_auditada/` y `ARQUITECTO DE DATOS/nucleo_cxx_auditado/` NO se modifican durante los pilotos.

---

## 1. Estructura de la carpeta piloto

```text
Pruebas del Operador de Lilyth_Correcciones Piloto/
├── inc/                         # Copia aislada de ark_cxx/base_auditada/include (70 archivos)
├── compilar.bat                 # Compilador MSVC C++20 (/W4 /O2 /fp:precise /I"inc")
├── build/                       # Binarios (.exe) y objetos (.obj) de las pruebas piloto
├── PLAN_Y_ESTRUCTURA_PILOTO.md  # Este plan y diseño de intervención por niveles
├── nivel1_nan_huecos/           # Nivel 1: separación de causas de NaN y corte sin compresión de rezagos
├── nivel2_arfima_memoria/       # Nivel 2: historia previa dinámica (d, p, radio_max) y simetría ARFIMA
├── nivel3_solape_entorno/       # Nivel 3: simetría de hay_previa y ventana t_max_ef en vía (b)
└── nivel4_admision_propia/      # Nivel 4: admisión en [t0, t0+l], escala_origen, g >= g_real y p_det
```

---

## 2. Plan de intervención — Nivel 1 (`nivel1_nan_huecos/`)

### Diagnóstico actual en código
1. En `ObservadorLilyth::paso()` y `guardar()` (`inc/ark/lilyth/observador_lilyth.hpp`), `hist_val_cruda_`, `hist_val_purgada_`, `hist_val_est_cruda_` y `hist_val_est_purgada_` son banderas binarias (`0`/`1`) que mezclan en el mismo estado:
   - dato crudo corrupto / no finito (`!valido || !std::isfinite(x)`),
   - barra imputada sin trades en `paso_filas()` (`b.imputada`),
   - primer paso válido en modo `NIVEL` (`!hay_prev_`),
   - fase todavía inmadura en `PurgaDiurna` (`!vp`).
2. En `copiar_perfil_estacionario()`, cuando aparece cualquier paso inválido dentro de `[t_ini_perf, t_apertura]` o `[t_prev_ini, t_ini_perf)`, el bucle simplemente lo salta (`if (... && std::isfinite(...))`) y comprime la secuencia, desalineando los rezagos de $A(L)(1-L)^d$ y mezclando segmentos separados.

### Cómo se ejecutará el cambio en `inc/ark/lilyth/observador_lilyth.hpp`
1. **Clasificación explícita de causa en el anillo (cero reservas en `paso()`):**
   - Incorporar `enum class CausaHuecoPaso : std::uint8_t`:
     - `VALIDO = 0`
     - `DATO_CORRUPTO = 1` (`!valido || !std::isfinite(x)`)
     - `BARRA_IMPUTADA = 2` (en `paso_filas()` cuando `b.imputada == true`, en lugar de confundirlo con `DATO_CORRUPTO`)
     - `INICIO_NIVEL = 3` (primer paso válido en `Integracion::NIVEL`)
     - `PURGA_INMADURA = 4` (dato válido cuya fase diurna aún no cumple `dias_minimos_purga`)
   - Añadir al anillo `hist_causa_cruda_` y `hist_causa_purgada_` (`std::array<CausaHuecoPaso, HIST_MAX>`) y exponer en `Contadores` el conteo separado por causa, sin penalizar en `Holm` ni inflar `corruptos` con barras vacías o purga inmadura.
2. **Corte estricto en `copiar_perfil_estacionario()` sin compresión temporal:**
   - **Para `hist_prev_*` (`[t_prev_ini .. t_ini_perf)`):** recorrer **hacia atrás desde `t_ini_perf - 1`** y detenerse (`break`) en el primer paso no válido (`hist_val_est_* == 0`, `!isfinite`) o cambio de `hist_seg_purgada_`, copiando únicamente el tramo continuo inmediatamente previo al evento en orden cronológico.
   - **Para `perfil_*` (`[t_ini_perf .. t_apertura]`):** recorrer **hacia atrás desde `t_apertura`** y detenerse ante cualquier hueco o cambio de segmento purgado; conservar la cola continua hasta `t_apertura` y registrar en la ranura si el perfil o `hist_prev` quedó truncado por hueco (`perfil_cortado_hueco`, `hist_prev_cortado_hueco`), en lugar de comprimir rezagos alrededor del hueco.

### Batería de verificación piloto (`nivel1_nan_huecos/`)
1. **Control C1 (Identidad en series limpias):** sin `NaN` ni cortes de segmento, el piloto debe coincidir bit a bit con `ark_cxx/base_auditada`.
2. **Control C2 (Corte limpio de `hist_prev`):** inyectar un `NaN` (y en otro caso un cambio de `seg_purgado_`) a $m$ pasos antes de `t_cambio`; verificar que `n_hist_prev_* == m` exactos (no salta el hueco ni copia muestras anteriores).
3. **Control C3 (No compresión en `perfil_est`):** inyectar un `NaN` dentro de `[t_cambio, t_apertura]`; verificar que no comprime el vector ni mezcla rezagos previos al hueco.
4. **Control C4 (Discriminación de causas):** alimentar `DATO_CORRUPTO`, `BARRA_IMPUTADA`, `INICIO_NIVEL` y `PURGA_INMADURA`; comprobar conteos separados en `contadores()` y cero penalización falsa en la familia de `Holm`.
5. **Regresión sobre baterías oficiales:** compilar y correr las baterías de `JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/lilyth/revision_etapas2_6/` apuntando a `inc/` del piloto (`/I"inc"`) verificando `/W4` limpio y `0` fallas.
