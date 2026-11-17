# Nota de Lilyth sobre el bufer de +-10 sesiones del muestreo de nulos

Code (Opus 5.5), 2026-11-16. Orden del Comandante: "anade la nota de Lilyth para quitar el buffer de +10 sesiones a ojo".
Destino: `ECONOMISTA DE CABECERA/CALIBRACION_DE_NULOS.md` (secciones 1-2) y `Carpeta del Comandante/CHECKLIST_MAESTRO_MUESTREO_NULOS.md`
(instruccion 4). Esta nota vive en el worktree de Code; no se edito ningun archivo del checkout principal. Copiarla alli es decision del
Comandante.

## Que se objeta

La directriz principal extirpa cada evento del catalogo con `M(E_k) = [t_ini - 10, t_fin + 10]` sesiones. El 10 se justifica con
mecanismos (fuga de informacion antes, relajacion despues), pero no se mide: es un numero a ojo, el mismo para un indice en 1d, una tasa,
una divisa y cripto a 1 minuto, y en sesiones aunque el muestreo sea intradia o en reloj intrinseco.

## Como se deriva, con los operadores que ya existen

Los dos lados del bufer responden a preguntas distintas y cada uno tiene su instrumento:

1. **Bufer posterior (relajacion).** Es el tiempo que el sistema tarda en volver a su calma tras el evento, y eso es exactamente lo que
   Lilyth mide: el tiempo de recomposicion `tau_rest` por evento. Regla: bufer posterior = cuantil `1 - alfa` de `tau_rest` de los
   eventos del catalogo del mismo activo, reloj y Delta t. Si ese cuantil esta censurado (eventos que no se recomponen dentro de T_max),
   el bufer es al menos T_max y se declara censurado; no se rellena.
2. **Bufer previo (fuga de informacion).** Es el retardo antes de `t_ini` desde el cual el comportamiento ya no se distingue de la calma.
   El instrumento es el operador delta: la tasa de alarmas en cada retardo previo al evento, contra su tasa de calma (1/ARL0). Regla:
   bufer previo = el retardo mas lejano en que la tasa de alarmas previa al evento supera a la de calma con la familia corregida (Holm).
   Si ningun retardo la supera, el bufer previo es cero: no hubo fuga medible.

Los dos salen en las unidades del reloj de la serie (pasos nativos o intrinsecos) y se traducen a sesiones solo para mostrarlos.

## La circularidad, y como se rompe

Lilyth necesita una calma para calibrarse y la calma necesita el bufer que Lilyth mide. Se rompe en dos pasadas declaradas:

- Pasada 1 (descriptiva, sin contraste): Lilyth corre con la calma de hoy (la del operador delta) solo para medir `tau_rest` de los
  eventos del catalogo; no se lee ningun p.
- Pasada 2: con el bufer derivado de la pasada 1 se fijan las ventanas de calma (para el operador delta) y de viento (para Lilyth), y
  entonces se calibra y se contrasta.

Si el bufer de la pasada 2 recalculado difiere del de la pasada 1 en mas de la resolucion del reloj, se reporta; no se itera hasta que
cuadre (eso seria ajustar el instrumento al resultado).

## Lo que no cambia

El catalogo de eventos y la jerarquia H0^Basal / H0^Expansion / H1 del Economista se conservan. Solo cambia de donde sale el ancho de la
mascara. Y por la ramificacion que propuso el Comandante, la calma resultante calibra el operador delta y las ventanas de viento
(eventos con su bufer derivado) quedan para Lilyth.
