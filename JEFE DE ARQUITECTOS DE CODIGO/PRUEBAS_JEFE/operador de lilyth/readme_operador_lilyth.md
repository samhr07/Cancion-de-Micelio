# Operador de Lilyth -- restitucion tras un evento (PREEXPECTATIVA, Etapa 1 en construccion)

Jefe de Arquitectos de Codigo (Opus 5.5) y Auditoria Externa (Argon), 2026-11-05. Diseno conversado con el Comandante sobre las fuentes:
- Fisico de Cabecera, `BITACORA DE VIAJE DEL FISICO.md`, entradas 38 (formulacion) y 40 (las cuatro sutilezas);
- Arquitecto de Datos, `BITACORA DEL ARQUITECTO.md`, capitulo XXVII (especificacion C++).

**Estado: preexpectativa cerrada con las correcciones de Code del 2026-11-05.** Fija que se construye, que NO se construye porque ya
existe, y que se espera medir antes de escribir cada etapa. Cuando haya codigo completo, el README pasa a
`ark_cxx/base_auditada/include/ark/lilyth/` con su manifiesto. Es telemetria: no decide operar. Datos: solo sinteticos por ahora.

## 1. Que pregunta responde, y cual no

- **Lilyth:** despues de un evento detectado, el observable vuelve a donde estaba? Cuanto tarda? Con que forma?
- **No es el operador de memoria.** `ObservadorMemoria` pregunta si cambio H. Lilyth pregunta por la trayectoria despues del cambio.
- Lo mas cercano en la base es la clasificacion SALTO / TRANSITORIO / DERIVA del operador delta. Lilyth no la sustituye: le anade el
  tiempo y la forma de la vuelta.

**Notacion (aprobada):** `tau_rest` = semivida de restitucion. El `tau50` de `FichaEscala` es OTRA cosa (la latencia del detector,
medida por inyeccion) y no se toca.

## 2. Cadena

```
observable Y (serie o filas, ya validadas por ark/ventanas)
  -> ObservadorSaltos (operador delta)          emite FichaEvento: t_cambio, signo, amplitud y su se, escala de origen
  -> Lilyth copia la ficha en una ranura        canal_nivel reinicia como siempre; Lilyth guarda historial
  -> Y* congelado en t_cambio^-                 nunca se actualiza mientras la ranura viva (separado para purgada y cruda)
  -> trayectoria eta(tau) = (Y - Y*) / dY       en reloj de pared y en reloj de actividad; sobre PURGADA, CRUDA o AMBAS
  -> tau_rest, v_rel, sobrepaso, censura
  -> EXCESO sobre los nulos configurados (sin memoria y lineal de calma)   (lo unico que se publica)
```

## 3. Decisiones del Comandante (2026-11-05)

| id | decision |
|---|---|
| D1 | **Disparo = evento del operador delta.** Sin umbral propio (el "\|dY\| > theta_det" de los apuntes era un numero a ojo). El presupuesto de falsas alarmas es el del operador delta (ARL0 o FAMILIA_EPISODIOS). |
| D2 | **Ancla Y\*:** mediana causal previa al evento, congelada. El "cuantil del vacio" ya existe: es el `mu0` del nulo del canal de nivel y se reutiliza. La moda solo para observables con reticulo real (spread en ticks). *(Sustituido en E1.3).* |
| D3 | **Calor:** lo detecta `canal_temperatura` (no se duplica). Lilyth suspende la medicion mientras dure. *(Sustituido en E1.6 y E1.7).* |
| D4 | **`canal_nivel` no se toca:** reinicia tras cada alarma (sin eso su ARL0 deja de valer). Lilyth copia la ficha en su ranura y guarda el historial para compararlo despues, que es lo que pedia la sutileza 1.2. |
| D5 | **Sutileza 2 del Fisico (denominador que colapsa): no aplica.** El GLR termico normaliza con el `s0` total de la serie nula de log beta; `var_log_beta` solo se reporta. Queda como guarda para quien use el exceso superestadistico. |
| D6 | **Cadencia:** la nativa de cada operador; S = 1 para los observables simples. Ya existe (`paso_sensor`). |
| D7 | **Reloj doble (preexpectativa, para no duplicar):** la infraestructura ya existe (`motor_reloj_intrinseco`, `FuenteFilas`). Lilyth solo anota el reloj en t0 y en el cruce. El reloj de actividad es NU_TRADES o NU_VOLUMEN, **no** el nu de Shannon: este acumula z^2, grande justo en el choque, y haria circular la comparacion. |
| D8 | **Replicas:** dentro de una cascada el ancla no se re-ancla (idea base aprobada). Hay que purgar los numeros a ojo: el "8 ranuras" pasa a ser un tope de memoria (como W_MAX), con el desborde contado y declarado; se retira el "coalescer la replica de menor energia". |
| D9 | **Elastico contra plastico:** se difiere al ultimo operador, el de red (tipo hipergrafo). Alli se decidira contra que se mide (la propia variable, el precio o la marea, que es la propuesta del Comandante). Aqui no se decide nada de eso. |
| D10 | **Armonico (hipotesis del Comandante):** una accion de control sobre una variable (cuentas institucionales absorbiendo una cascada) se traduce en una respuesta oscilatoria. Se conserva y solo se descarta con datos. Contraste propuesto en la seccion 5 (P5) y precisado en E1.8 y E1.13. |
| D11 | **Holm:** Lilyth mide y no emite pruebas propias en esta version; si un dia declara algo, pasa por `holm.hpp`. *(Sustituido en E1.16: el resumen contrasta cuatro hipotesis por `holm.hpp`).* |
| D12 | **Primero lo lineal:** la primera medicion es la prueba de respuesta lineal (Onsager, P3). Lo no lineal solo si lo lineal no lo explica. |

## 4. Flancos que los apuntes no cubrian (aprobados)

| id | flanco | resolucion |
|---|---|---|
| F1 | **Regresion a la media.** Se dispara en un extremo, y hasta un ruido iid "se relaja un 50 %" en un paso. | **Nulo obligatorio:** el mismo disparo y la misma medida sobre la serie nula del operador delta. Se publica el exceso de `tau_rest` sobre su nulo, nunca el `tau_rest` crudo. *(Precisado en E1.1 y 8.2).* |
| F2 | **Primer paso con ruido:** el cruce de 0.5 \|dY\| ocurre pronto por puro ruido. | dY = la amplitud del GLR (media sobre l pasos, con su se). El cruce se exige con una banda relativa a ese se. *(Sustituido por E1.4).* |
| F3 | **La Capa 1 (`SensorEkaLilyth`) reintroduce el artefacto:** sus medianas moviles de W absorben el choque y L_Eka vuelve a 1 por construccion antes de W. | Congelar tambien las medianas de la Capa 1 en t0. *(Precisado en E1.12).* |
| F4 | **Fin del calor:** una Page-CUSUM con max(0, .) toca 0 a menudo por ruido y declararia "condensado" demasiado pronto. | El fin del calor se calibra como la alarma: una racha fria de longitud m con su propio ARL sobre el nulo. *(Sustituido por E1.7).* |
| F5 | **Lo medido antes, leido con su ventana.** | Ver la seccion 6. |

## 5. Preexpectativas (escritas antes del codigo; todas sobre sinteticos con verdad conocida)

| id | mundo | se espera |
|---|---|---|
| P1 | ruido iid con el disparo del operador delta | `tau_rest` crudo de ~1 paso (regresion a la media) y exceso sobre el nulo compatible con 0. Si el crudo no sale corto, el control no es tal. |
| P2 | AR(1) con phi conocido y choques inyectados | el exceso recupera la semivida ln 2 / (-ln phi) dentro de su error. |
| P3 | Onsager | en el AR(1) lineal, el `tau_rest` de los choques coincide con el que predice la ACF de la calma (el `tau_int` del nulo, que ya calcula el canal de nivel; en tiempo discreto tau_int bilateral = (1 + rho)/(1 - rho), no 1/lambda). En un mundo no lineal (umbral o TAR) difieren con el signo conocido. |
| P4 | reloj doble | con la actividad congelada: `tau_rest` grande en pared y normal en el reloj de actividad. Con perdida de rigidez: grande en los dos. |
| P5 | armonico (D10) | un AR(2) subamortiguado (verdad conocida) produce **sobrepaso**: Y - Y* cambia de signo mas alla de su banda despues de t0. El contraste de sobrepaso contra su nulo lo detecta; en el AR(1) no. Se contrasta sobre la trayectoria posterior al evento, no sobre la PSD de la serie estacionaria: una accion de control aparece en la respuesta al choque, no tiene por que aparecer en reposo. |
| P6 | Capa 1 | sin congelar sus medianas, `tau_100` <= W por construccion (el control demuestra el artefacto); con las medianas congeladas, no. |
| P7 | cascadas | la ranura no re-ancla; el desborde se cuenta; la censura deja ver la cola lenta. |

No decide: datos reales, operar, ni elastico contra plastico (D9).

## 6. Lo medido antes, con su alcance (correcciones del Comandante)

- **El Lilyth macro era sobre todo la ventana:** tau proporcional a W^0.75 (orden 23, `acta_lilyth_recovery.txt`). Es el artefacto que
  el ancla congelada evita.
- **La elasticidad cambia en el tiempo:** 2024 fue mas elastico que 2026 en mediana, y hace unas semanas el mercado volvio a serlo
  (semana del 2026-10-19 ELASTICA; ago-oct 2026 parecido a 2024; `acta_termometro_semanal_marea.txt`). Promediar un ano entero y
  concluir que "no vuelve" seria decirse mentiras.
- **`e < 0` de HIFAS** (el medio no rebota tras los barridos, a escala de segundos) es de W1 y W1P: una propiedad de esas semanas y de
  esa escala, no del mercado.
- **Semividas por evento en L3** (HIFAS, T_c 2.1-2.6 s) son sustancia, y la ley de Lilyth cambio de signo con dD causal: se toman en
  cuenta al disenar la Capa 1.
- **No se juzga el radar por si caza cisnes.** Sera general; ese juicio, si llega, es posterior y aparte.
- Las cifras del capitulo XXVII ("1.59 M snapshots", "SNR 0.57 sigma", "epsilon_calma = 0.000") no tienen acta y no se citan hasta
  tenerla.

## 7. Decisiones cerradas por el Comandante (2026-11-05, tras la consulta a Argon)

1. D10 mientras no este verificado: **etiqueta** (`ARMONICO_SOSPECHOSO`): el evento se mide igual y se marca; nunca veto.
2. Replicas: **riesgos competitivos** (incidencia acumulada de Aalen-Johansen, con "restituye" y "llega una replica" como
   desenlaces), no Kaplan-Meier.
3. W del ancla: **derivada de la escala de origen del evento** (N_REG 2^j0 barras), conectada al resto del pipeline.

## 8. Enmienda 1 (2026-11-05, antes de cualquier codigo): dictamen de Argon, correcciones de Code y directivas del Comandante

Fuentes: dos dictamenes de Argon (`Carpeta del Comandante/Conversaciones Argon-Code/dictamen_auditor_externo_operador_lilyth.md` y su
cierre consolidado), la revision de Code contra `escalera_diadica.hpp` y `tipos_metrologia.hpp`, y las directivas del Comandante.
**Donde esta seccion contradice a las secciones 3-5, manda esta.**

### 8.1 Que cambia y por que

| id | antes | ahora | quien | por que |
|---|---|---|---|---|
| E1.1 | un solo nulo | **dos referencias**: exceso sobre el nulo **sin memoria** (PERMUTACION o SIGNOS: Delta_tau_iid) y exceso de **Onsager** sobre un nulo **lineal con la ACF de la calma** (Delta_tau_Onsager) | Argon | con un solo nulo, P2 y P3 se contradecian (si se resta el AR(1) de calma, el exceso de P2 sale 0) |
| E1.2 | el nulo podia ser NULO_NINGUNO | Lilyth **exige un sustituto sin memoria valido**; con `NULO_NINGUNO` o `NULO_BLOQUES` en `nulo_sin_memoria` se niega a calibrar | Argon, Code | restar la serie real de si misma da exceso 0; `BLOQUES` conserva la memoria dentro del bloque y no es un nulo sin memoria |
| E1.3 | Y* = mediana; "cuantil del vacio = mu0" | Y* = **mediana causal de Y en W_ancla = N_REG 2^j0 barras nativas**, terminada en t_cambio - ceil(tau_int,j0) 2^j0 - 1 (guarda contra el inicio de una rampa), copiada en la ranura al abrir. `mu0` **no** es ancla | Argon (W: Comandante) | con NIVEL, `mu0` es la deriva por barra (~0), no el nivel previo |
| E1.4 | dY = amplitud del GLR; cruce con banda en se | **dY = Y(t_a) - Y\***, con signo(dY) = signo de la alarma; `amplitud_glr` se guarda aparte. **eta(tau) = (Y(t_a + tau) - Y\*) / dY**, eta(0) = 1. `tau_rest` = primer tau >= 1 con eta <= 0.5, **sin banda** | Argon | con NIVEL la amplitud es una pendiente por barra; la banda en unidades de se vale 1/sqrt(2g) y fabrica exceso en los choques fuertes. eta con signo ve el cruce aunque salte al otro lado |
| E1.5 | Onsager contra tau_int | Onsager por **P3a** (el nulo lineal pasado por el mismo detector, a igual SNR) y **P3b** (la curva media eta_barra(tau) sin detenerse en 0.5, que en un proceso lineal sigue la ACF). Semivida de la media: ln 2 / (-ln phi), phi = (tau_int - 1)/(tau_int + 1) | Argon; precision de Code | el primer cruce se adelanta por la difusion; con memoria mayor que la de un AR(1) el GLR condiciona en la trayectoria y la formula ya no es exacta |
| E1.6 | "suspender durante el calor" | el reloj corre desde la alarma siempre; se exportan **`tau_rest_primer_toque`** y **`tau_rest_consolidado`** (eta sigue <= 0.5 al terminar el calor) | Argon | muchos choques vuelven estando calientes |
| E1.7 | fin del calor por "racha fria con ARL sobre el nulo de calma" (F4) | **deteccion secuencial unilateral hacia mu_0,beta**, calibrada sobre el nulo de log beta^ **desplazado a mu_calor**, con ARL_0,caliente = T_max / alfa (aproximacion geometrica, declarada). **mu_calor sale de las ventanas posteriores a la alarma**, no de su amplitud | Argon; mu_calor de Code | F4 invertia la hipotesis nula. La racha de 5 ventanas suponia independencia y tardaba ~62 ventanas. La amplitud de la alarma esta sesgada por la seleccion del GLR |
| E1.8 | sobrepaso en [tau_rest, 4 tau_rest] | sobrepaso en la **curva media**: min de eta_barra(tau) para tau en [2, T_max], contra su distribucion en el nulo. AR(2) por **Prony** (minimos cuadrados de la recurrencia) sobre eta_barra; oscilador = **raices complejas** (a1^2 + 4 a2 < 0) y se reporta su periodo **sin recortarlo**. Dos contrastes: contra el nulo iid (hay oscilacion?) y contra el lineal de calma (hay mas que en calma? Es el que contrasta la hipotesis del Comandante). Etiqueta por evento con su tasa en el nulo al lado | Code (ventana y raices); Argon (estadistico) | el 4 no es universal cerca del amortiguamiento critico; el ruido browniano se promedia en la media. **El rebote bid-ask (phi < 0) es una raiz real negativa y nunca da raices complejas**. Se descarta el "a1 > 0, periodo > 4" del cierre de Argon: elimina osciladores de 2-4 barras sin razon; Yule-Walker no aplica a una curva de relajacion |
| E1.9 | T_max sin definir | **T_max = cuantil 1 - alfa de `tau_rest` en el nulo lineal**; si el nulo lineal no restituye (raiz unitaria, p. ej. precio en nivel con nulo de paseo aleatorio), **T_max = W_ancla** y se marca | Argon; caso raiz unitaria de Code | sin T_max el RMST no esta definido; en un paseo aleatorio el cuantil no existe |
| E1.10 | Kaplan-Meier o riesgos competitivos | **Aalen-Johansen**: incidencias CIF_rest y CIF_replica, RMST(T_max) = suma (1 - CIF_rest), desglosado por signo del choque | Comandante, Argon | una replica no es censura no informativa |
| E1.11 | reloj de actividad | NU_TRADES, NU_VOLUMEN o NU_NOCIONAL (nunca NU_INFO_SHANNON); se guarda nu_punto* = (nu en W_ancla) / W_ancla y se exporta tau_act = Delta_nu_rest / nu_punto* en barras de calma | Argon | P4 compara tau_pared con tau_act en las mismas unidades |
| E1.12 | Capa 1 congelada en t0 | las medianas se **copian por valor en cada ranura** en t_cambio^-; el anillo del sensor sigue vivo | Argon | varias ranuras a la vez; el anillo sirve a eventos futuros |
| E1.13 | nulos y detector armonico fijos por dentro | **ENTRADAS de configuracion del pipeline** (`ConfiguracionLilyth`): `nulo_sin_memoria` (**solo `NULO_PERMUTACION` o `NULO_SIGNOS`**; `BLOQUES` y `NINGUNO` se rechazan), `modo_nulo_lineal` y `modo_armonico`. En `PRONY_FASE_REESCALADA` ($\tilde{\tau} = \tau/\tau_{\text{rest},k}$), los eventos que no restituyen ($\tau_{\text{rest}} = \infty$) quedan excluidos por construccion y **se declara en la ficha cuantos quedaron fuera** | Comandante; correcciones de Code | `BLOQUES` conserva memoria interna; excluir $\tau = \infty$ al dividir por $\tau_{\text{rest}}$ es una seleccion que debe declararse |
| E1.14 | medir solo sobre la purgada o solo sobre el crudo | **ENTRADA de representacion** (`ModoRepresentacionLilyth`: `PURGADA`, `CRUDA`, `AMBAS_CONTRASTE_MAREA`). `PurgaDiurna` (`escalera_diadica.hpp`) aplica $y_{\text{purgado}} = y \cdot \sqrt{v_{\text{media}} / v_{\text{fase}}}$ con el segundo momento crudo por fase (no resta $\mu_0$ ni corrige la media aparte; en observables positivos como volumen arrastra tambien la media diurna; en retornos, solo la varianza). Necesita `dias_minimos` (20) dias por fase antes de ser valida: hasta entonces la rama `PURGADA` **no mide**. En `NIVEL`, un incremento purgado invalido rompe la reintegracion: la ranura queda **censurada por fin de dato y se declara**, nunca se rellena | Comandante, Code, Argon | refleja la operacion exacta de `PurgaDiurna` en `escalera_diadica.hpp` y prohibe imputar en la reintegracion de `NIVEL` |
| E1.15 | duda sobre si disparar tambien con eventos de temperatura | En Etapa 1 Lilyth se dispara **solo con los eventos del canal de nivel** (los de temperatura gobiernan el estado de calor y su fin). La **Segunda Familia (Restitucion de Temperatura: cuanto tarda $\log \hat{\beta}$ en volver tras un `CAMBIO_TEMPERATURA`)** queda registrada como pendiente para la version posterior | Comandante, Code | separa limpiamente la restitucion de nivel (con el calor como modulador) de la relajacion termica pura |
| E1.16 | D11 decia que Lilyth no emite pruebas propias | El resumen de Lilyth contrasta **cuatro hipotesis por `ark/holm.hpp`**: $\Delta\tau_{\text{iid}}$, $\Delta\tau_{\text{Onsager}}$, sobrepaso contra el nulo iid y sobrepaso contra el nulo lineal de calma | Code | actualiza D11 a la estructura real del resumen |

### 8.2 Como se generan los eventos del nulo (aprobado: (a) + (b))

- **(a) Espontaneos.** B realizaciones del sustituto, con el detector de **produccion** (umbrales congelados, sin recalibrar). Para cada
  evento real k se usan los eventos nulos con **g >= g_k** (excedencia: no hay bins). Vale (a) si hay al menos m_min(alfa) = 8; B por la
  regla de `calibrador_familia.hpp` con tope B_max.
- **(b) Inyeccion generativa** (cuando (a) no llega a m_min). **Correccion de Code al cierre de Argon:** superponer el perfil del choque
  en un sustituto ya calculado no se propaga (la serie no responde): en FLUJO daria restitucion instantanea y en NIVEL un desplazamiento
  permanente, sea cual sea la memoria. Por eso:
  - nulo sin memoria: no hace falta simular; tau_iid = 1 por definicion (o el medido sobre el sustituto permutado/signos);
  - nulo lineal: **modelo generativo** configurable (`ModoNuloLineal`, E1.13). AR(p) ajustado a la ACF de la calma por Levinson-Durbin
    (con el orden $p \le K_{\text{corte}}$ por Quenouille, el ultimo $|\kappa_p| > z_{1-\alpha/2}/\sqrt{n}$ dentro de `FichaNulo.rezago_corte`;
    se declara que probar $K$ rezagos a nivel $\alpha$ tiende a elegir un $p$ ligeramente mayor, que es el error conservador al sobreajustar
    un poco el nulo lineal) o ARFIMA(p, d, 0) con la $d$ del Whittle exacto. **Nota sobre `FASE_ALEATORIZADA_FFT`:** solo sirve para el
    regimen espontaneo (a); si se eligio como nulo lineal y un evento cae en la via (b), Lilyth cae automaticamente al modelo generativo
    (una serie FFT fija no se propaga, P8) y lo marca en la ficha (`nulo_fft_degradado_a_generativo = true`).

### 8.3 Preexpectativas nuevas o cambiadas (sustituyen a las de la seccion 5 donde coinciden)

| id | mundo | se espera |
|---|---|---|
| P1 | iid | Delta_tau_iid ~ 0 y Delta_tau_Onsager ~ 0 |
| P2 | AR(1) phi | Delta_tau_iid ~ ln 2 / (-ln phi) (corregido por el adelanto de primer cruce, medido); Delta_tau_Onsager ~ 0 |
| P2b | AR(1) phi, solo nulo (b) | el modelo generativo inyectado recupera la misma semivida que (a) dentro de su error |
| P3a | TAR (umbral) con la misma ACF de calma | Delta_tau_Onsager != 0 con el signo del umbral |
| P3b | AR(1) y AR(2) lineales | eta_barra(tau) sigue la ACF de calma dentro de 2 ee |
| P5 | AR(2) subamortiguado | raices complejas en eta_barra y min eta_barra < 0 contra el nulo |
| P5b | AR(1) con phi < 0 (rebote) | raiz real negativa: **nunca** oscilador |
| P8 | mutante: inyeccion en serie fija | el control (b) tiene que fallar (restitucion en tau = 1 en FLUJO) |
| P9 | fin del calor | con el calor apagado en un t conocido, el fin declarado en falso antes de t tiene tasa <= alfa en T_max, y el retraso se mide |
| P10 | raiz unitaria | paseo aleatorio con nulo de paseo: T_max = W_ancla marcado; Delta_tau = 0 |
| P11 | contraste de marea (E1.14) | con estacionalidad multiplicativa $y\sqrt{v_{\text{media}}/v_{\text{fase}}}$: tras madurar los 20 dias por fase, `PURGADA` recupera la semivida verdadera; antes de los 20 dias `PURGADA` no mide; en `NIVEL` un incremento purgado invalido censura la ranura sin rellenar |

## 9. Estructura del codigo y Etapa 1

* **Cabeceras de la Etapa 1 (`ark/lilyth/`):**
  - `tipos_lilyth.hpp`: contratos, enumeraciones (`ModoRepresentacionLilyth`, `ModoNuloLineal`, `ModoArmonico`, `DesenlaceRanura`),
    `ConfiguracionLilyth`, `RamaTrayectoriaLilyth<T_MAX>`, `RanuraLilyth<T_MAX>`, `TablaAalenJohansen<T_MAX>`, `FichaResumenLilyth<T_MAX>`.
  - `aalen_johansen.hpp`: acumulador causal y fuera de linea de riesgos competitivos (`CIF_rest`, `CIF_replica`, `S`, `RMST(T_max)`)
    desglosado en total, signo $+$ y signo $-$.
  - `observador_lilyth.hpp`: orquestador `ObservadorLilyth<Obs, RANURAS, T_MAX>` que envuelve o lee la salida de `Obs` (`ObservadorSaltos`)
    sin modificar `ObservadorSaltos` bit a bit, gestiona el anillo causal pre-evento para congelar $Y^*$ (cruda y purgada reintegrada),
    abre y avanza las `RANURAS` en ciclo caliente con cero reservas de memoria, registra $\eta(\tau)$, `tau_rest_primer_toque`,
    `tau_rest_consolidado`, réplicas que cortan trayectorias y desborde de ranuras.