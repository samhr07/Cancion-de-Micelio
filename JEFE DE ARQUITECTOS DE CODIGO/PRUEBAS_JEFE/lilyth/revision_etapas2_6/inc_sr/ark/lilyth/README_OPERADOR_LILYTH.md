# Operador de Lilyth (`ark/lilyth/`) -- restitucion tras un evento del operador delta. ETAPA 1

Argon (Auditoria Externa) escribio la Etapa 1 y Code (Jefe de Arquitectos de Codigo, Opus 5.5) la reviso contra el operador delta real
y la conecto a la base, 2026-11-05. Diseno, decisiones del Comandante y preexpectativas (P1-P12):
`JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/lilyth/README_OPERADOR_LILYTH.md` (canonico; este README pasa a sustituirlo al cerrar la
Etapa 8). Revision, controles y mutantes: `PRUEBAS_JEFE/lilyth/revision_etapa1/` (`REVISION_ETAPA1.md`). Es telemetria: no decide operar.

## Que hace en la Etapa 1

```
ObservadorSaltos (calibrado)  --Salida por referencia constante-->  ObservadorLilyth
   alarma de nivel con el evento cerrado (espejo de la regla del operador delta) = apertura en t_a = alarma.t
   ancla Y* = mediana causal de [t_cambio - guarda - W_ancla, t_cambio - guarda - 1], W_ancla = N_REG 2^j0, guarda = ceil(tau_int) 2^j0
   eta(tau) = (Y(t_a + tau) - Y*) / (Y(t_a) - Y*) en las ramas PURGADA y CRUDA
   tau_rest_primer_toque, tau_rest_consolidado (sin calor), reloj de actividad y de pared, replicas, Aalen-Johansen por signo
```

## Uso minimo

```cpp
#include "ark/lilyth/observador_lilyth.hpp"
using Obs = ark::metrologia::ObservadorSaltos<>;
auto delta = std::make_unique<Obs>();
delta->calibrar(cd, x_cal, {}, v_cal, f_cal, d_cal);
auto lil = std::make_unique<ark::lilyth::ObservadorLilyth<Obs>>();   // grande: en el monticulo
ark::lilyth::ConfiguracionLilyth cl{};
cl.t_max_efectivo = T_max;                     // obligatorio en la Etapa 1 (lo deriva la Etapa 3)
lil->configurar(cl, cd, *delta, x_cal, v_cal, f_cal, d_cal);   // la MISMA historia que recibio calibrar (espejo de la purga)
for (...) {
    const auto s = delta->paso(x, NaN, valido, fase, dia);
    lil->fijar_calor_activo(calor);            // Etapa 1: entrada; el fin del calor calibrado es la Etapa 5
    const auto o = lil->paso(s, x, valido, d_nu, fase, dia);   // las MISMAS entradas que el operador delta
    for (std::size_t i = 0; i < o.ranuras_cerradas; ++i) usar(lil->cerrada(i));
}
```

Por filas: `lil->paso_filas(salida_de_ObservadorSaltosFilas::paso_fila, fase_de_barra, fase, dia)` (nu y pared salen de la barra).

## Contrato

- No posee ni toca al operador delta. Su reloj es el `Salida.t` del operador delta: un paso saltado se cuenta (`pasos_desalineados`) y
  no se mide.
- N_REG y B_MAX salen del tipo del operador delta (`RasgosDelta`). La integracion y la purga salen de su `ConfiguracionObservador`.
- La purga de la rama PURGADA es un espejo de la `PurgaDiurna` del operador delta: necesita la historia de calibracion porque la del
  operador delta la acumula y la conserva en el ciclo caliente.
- Sin numeros a mano: T_max sin fijar (0) no configura; RANURAS, T_MAX y HIST_MAX son topes de memoria con su desborde contado.
- Cada perdida tiene su motivo (`DesenlaceRanura`): purga inmadura, ancla insuficiente (historia corta o corte de la reintegracion en
  NIVEL), ancla fuera de memoria, signo inconsistente.

## Limites declarados (Etapa 1)

- **Sin nulos todavia:** `tau_rest` crudo lleva la regresion a la media (F1). Lo publicable es el exceso sobre los dos nulos de la
  Etapa 3; hasta entonces esto es la maquinaria.
- **Replicas del mismo choque.** Con un desplazamiento sostenido el GLR vuelve a alarmar respecto de mu0, el operador delta cierra su
  evento antes de que el observable vuelva y reabre otro sobre el mismo choque: Lilyth lo cuenta como replica (corta y abre una
  ranura con el ancla ya dentro del choque). Medido en `REVISION_ETAPA1.md`; la regla para distinguir "continuacion" de "replica" es
  una decision pendiente del Comandante.
- Los eventos de los primeros W_ancla + guarda pasos del ciclo caliente no tienen ancla (la historia de Lilyth empieza en el ciclo
  caliente).
- `armonico_sospechoso` por evento es telemetria provisional: cualquier eta < 0 en tau >= 2 la enciende, incluida la vuelta por
  ruido. El contraste es de la Etapa 4.
