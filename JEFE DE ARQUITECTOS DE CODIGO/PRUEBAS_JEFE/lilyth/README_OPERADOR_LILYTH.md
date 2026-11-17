# Operador de Lilyth -- restitucion tras un evento (PREEXPECTATIVA, sin codigo)

Jefe de Arquitectos de Codigo (Opus 5.5) y Auditoria Externa (Argon), 2026-11-05. Diseno conversado con el Comandante sobre las fuentes:
- Fisico de Cabecera, `BITACORA DE VIAJE DEL FISICO.md`, entradas 38 (formulacion) y 40 (las cuatro sutilezas);
- Arquitecto de Datos, `BITACORA DEL ARQUITECTO.md`, capitulo XXVII (especificacion C++).

**Estado (2026-11-07): Etapas 1-6 construidas (Argon), corregidas y conectadas a `ark_cxx/base_auditada/include/ark/lilyth/` (Code; `revision_etapa1/`, `revision_etapas2_6/`). E1.21 implementada. Etapa 7 en espera (necesita modificar el dataset); Etapa 8 de Code.** Lo que sigue es la preexpectativa original, cerrada antes del codigo. Fija que se construye, que NO se construye porque ya existe, y que se espera medir antes de escribir una
linea. Cuando haya codigo, el README pasa a `ark_cxx/base_auditada/include/ark/lilyth/` con su manifiesto. Es telemetria: no decide
operar. Datos: solo sinteticos por ahora.

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
| D2 | **Ancla Y\*:** mediana causal previa al evento, congelada. El "cuantil del vacio" ya existe: es el `mu0` del nulo del canal de nivel y se reutiliza. La moda solo para observables con reticulo real (spread en ticks). *(Precisado en E1.3).* |
| D3 | **Calor:** lo detecta `canal_temperatura` (no se duplica). Lilyth suspende la medicion mientras dure. *(Precisado en E1.6 y E1.7).* |
| D4 | **`canal_nivel` no se toca:** reinicia tras cada alarma (sin eso su ARL0 deja de valer). Lilyth copia la ficha en su ranura y guarda el historial para compararlo despues, que es lo que pedia la sutileza 1.2. |
| D5 | **Sutileza 2 del Fisico (denominador que colapsa): no aplica.** El GLR termico normaliza con el `s0` total de la serie nula de log beta; `var_log_beta` solo se reporta. Queda como guarda para quien use el exceso superestadistico. |
| D6 | **Cadencia:** la nativa de cada operador; S = 1 para los observables simples. Ya existe (`paso_sensor`). |
| D7 | **Reloj doble (preexpectativa, para no duplicar):** la infraestructura ya existe (`motor_reloj_intrinseco`, `FuenteFilas`). Lilyth solo anota el reloj en t0 y en el cruce. El reloj de actividad es NU_TRADES o NU_VOLUMEN, **no** el nu de Shannon: este acumula z^2, grande justo en el choque, y haria circular la comparacion. |
| D8 | **Replicas:** dentro de una cascada el ancla no se re-ancla (idea base aprobada). Hay que purgar los numeros a ojo: el "8 ranuras" pasa a ser un tope de memoria (como W_MAX), con el desborde contado y declarado; se retira el "coalescer la replica de menor energia". |
| D9 | **Elastico contra plastico:** se difiere al ultimo operador, el de red (tipo hipergrafo). Alli se decidira contra que se mide (la propia variable, el precio o la marea, que es la propuesta del Comandante). Aqui no se decide nada de eso. |
| D10 | **Armonico (hipotesis del Comandante):** una accion de control sobre una variable (cuentas institucionales absorbiendo una cascada) se traduce en una respuesta oscilatoria. Se conserva y solo se descarta con datos. Contraste propuesto en la seccion 5 (P5) y precisado en E1.8 y E1.13. |
| D11 | **Holm:** Lilyth mide y no emite pruebas propias en esta version; si un dia declara algo, pasa por `holm.hpp`. *(Sustituido por E1.16).* |
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

## 8. Enmienda 1 (2026-11-05, antes de cualquier codigo): dictamen de Argon, correcciones de Code y directivas finales del Comandante

Fuentes: dos dictamenes de Argon (`Carpeta del Comandante/Conversaciones Argon-Code/dictamen_auditor_externo_operador_lilyth.md` y su
cierre consolidado), la lectura de Code contra el codigo y las tres directivas adicionales del Comandante previas a la Etapa 1.
**Donde esta seccion contradice a las secciones 3-5, manda esta.**

### 8.1 Que cambia y por que

| id | antes | ahora | quien | por que |
|---|---|---|---|---|
| E1.1 | un solo nulo | **dos referencias**: exceso sobre el nulo **sin memoria** (PERMUTACION o SIGNOS: Delta_tau_iid) y exceso de **Onsager** sobre un nulo **lineal con la ACF de la calma** (Delta_tau_Onsager) | Argon | con un solo nulo, P2 y P3 se contradecian (si se resta el AR(1) de calma, el exceso de P2 sale 0) |
| E1.2 | el nulo podia ser NULO_NINGUNO | Lilyth **exige un sustituto**; con NINGUNO se niega a calibrar | Argon | restar la serie real de si misma da exceso 0 por construccion |
| E1.3 | Y* = mediana; "cuantil del vacio = mu0" | Y* = **mediana causal de Y en W_ancla = N_REG 2^j0 barras nativas**, terminada en t_cambio - ceil(tau_int,j0) 2^j0 - 1 (guarda contra el inicio de una rampa), copiada en la ranura al abrir. `mu0` **no** es ancla | Argon (W: Comandante) | con NIVEL, `mu0` es la deriva por barra (~0), no el nivel previo |
| E1.4 | dY = amplitud del GLR; cruce con banda en se | **dY = Y(t_a) - Y\***, con signo(dY) = signo de la alarma; `amplitud_glr` se guarda aparte. **eta(tau) = (Y(t_a + tau) - Y\*) / dY**, eta(0) = 1. `tau_rest` = primer tau >= 1 con eta <= 0.5, **sin banda** | Argon | con NIVEL la amplitud es una pendiente por barra; la banda en unidades de se vale 1/sqrt(2g) y fabrica exceso en los choques fuertes. eta con signo ve el cruce aunque salte al otro lado |
| E1.5 | Onsager contra tau_int | Onsager por **P3a** (el nulo lineal pasado por el mismo detector, a igual SNR) y **P3b** (la curva media eta_barra(tau) sin detenerse en 0.5, que en un proceso lineal sigue la ACF). Semivida de la media: ln 2 / (-ln phi), phi = (tau_int - 1)/(tau_int + 1) | Argon; precision de Code | el primer cruce se adelanta por la difusion; con memoria mayor que la de un AR(1) el GLR condiciona en la trayectoria y la formula ya no es exacta |
| E1.6 | "suspender durante el calor" | el reloj corre desde la alarma siempre; se exportan **`tau_rest_primer_toque`** y **`tau_rest_consolidado`** (eta sigue <= 0.5 al terminar el calor) | Argon | muchos choques vuelven estando calientes |
| E1.7 | fin del calor por "racha fria con ARL sobre el nulo de calma" (F4) | **deteccion secuencial unilateral hacia mu_0,beta**, calibrada sobre el nulo de log beta^ **desplazado a mu_calor**, con ARL_0,caliente = T_max / alfa (aproximacion geometrica, declarada). **mu_calor sale de las ventanas posteriores a la alarma**, no de su amplitud | Argon; mu_calor de Code | F4 invertia la hipotesis nula. La racha de 5 ventanas suponia independencia y tardaba ~62 ventanas. La amplitud de la alarma esta sesgada por la seleccion del GLR |
| E1.8 | sobrepaso en [tau_rest, 4 tau_rest] | sobrepaso en la **curva media**: min de eta_barra(tau) para tau en [2, T_max], contra su distribucion en el nulo. AR(2) por **Prony** (minimos cuadrados de la recurrencia) sobre eta_barra; oscilador = **raices complejas** (a1^2 + 4 a2 < 0) y se reporta su periodo **sin recortarlo**. Dos contrastes: contra el nulo iid (hay oscilacion?) y contra el lineal de calma (hay mas que en calma? Es el que contrasta la hipotesis del Comandante). Etiqueta por evento con su tasa en el nulo al lado | Code (ventana y raices); Argon (estadistico) | el 4 no es universal cerca del amortiguamiento critico; el ruido browniano se promedia en la media. **El rebote bid-ask (phi < 0) es una raiz real negativa y nunca da raices complejas**. Se descarta el "a1 > 0, periodo > 4" del cierre de Argon: elimina osciladores de 2-4 barras sin razon; Yule-Walker no aplica a una curva de relajacion |
| E1.9 | T_max sin definir | **T_max = cuantil 1 - alfa de `tau_rest` en el nulo lineal**; si el nulo lineal no restituye (raiz unitaria, p. ej. precio en nivel con nulo de paseo aleatorio), **T_max = W_ancla** y se marca | Argon; caso raiz unitaria de Code | sin T_max el RMST no esta definido; en un paseo aleatorio el cuantil no existe |
| E1.9b (2026-11-16) | E1.9 deja ~alfa de los nulos censurados en T_max + 1 | `ConfiguracionLilyth::horizonte` = CUANTIL_ALFA (E1.9, por omision) o **CUANTIL_RESOLUCION: T_max = cuantil 1 - alfa/8 del nulo lineal, guarda de censura en alfa/8** | Comandante (orden), Code | un evento real lento sale censurado y empata con todos los nulos censurados: su p minimo es del orden de esa fraccion por grande que sea el pool; resolver el primer escalon del Holm de 4 (alfa/4) en la prueba bilateral pide c <= alfa/(2 * 4) = alfa/8. Requisito medido: estimar ese cuantil pide >= 8/alfa = 160 nulos lineales en la via (a) (hoy salen ~11): ver `PRUEBAS_JEFE/lilyth/horizonte_resolucion/MEDIDO_HORIZONTE.md` |
| E1.10 | Kaplan-Meier o riesgos competitivos | **Aalen-Johansen**: incidencias CIF_rest y CIF_replica, RMST(T_max) = suma (1 - CIF_rest), desglosado por signo del choque | Comandante, Argon | una replica no es censura no informativa |
| E1.11 | reloj de actividad | NU_TRADES, NU_VOLUMEN o NU_NOCIONAL (nunca NU_INFO_SHANNON); se guarda nu_punto* = (nu en W_ancla) / W_ancla y se exporta tau_act = Delta_nu_rest / nu_punto* en barras de calma | Argon | P4 compara tau_pared con tau_act en las mismas unidades |
| E1.12 | Capa 1 congelada en t0 | las medianas se **copian por valor en cada ranura** en t_cambio^-; el anillo del sensor sigue vivo | Argon | varias ranuras a la vez; el anillo sirve a eventos futuros |
| E1.13 | nulos y detector armonico fijos por dentro | **ENTRADAS de configuracion del pipeline** (`ConfiguracionLilyth`): tanto la pareja de nulos (`nulo_sin_memoria`, `modo_nulo_lineal`) como la politica del detector de armonicos (`modo_armonico`) son entradas explicitas para poder adaptarlas segun la naturaleza de cada variable | Comandante | igual que en `README_NULOS.md`, ningun operador clava un nulo ni un detector de forma rigida por dentro |
| E1.14 | medir solo sobre la purgada o solo sobre el crudo | **ENTRADA de representacion** (`ModoRepresentacionLilyth`: `PURGADA`, `CRUDA`, `AMBAS_CONTRASTE_MAREA`). Cuando se eligen ambas, Lilyth mide en paralelo sobre la serie purgada (en NIVEL, reintegrando causalmente los incrementos purgados) y sobre la serie cruda, y publica ademas la diferencia $\Delta\tau_{\text{marea}} = \tau_{\text{rest}}^{\text{cruda}} - \tau_{\text{rest}}^{\text{purgada}}$ | Comandante, Code, Argon | permite contrastar directamente que parte de la restitucion se debe a la marea diurna exogena y que parte es invariante endogena |
| E1.15 | duda sobre si disparar tambien con eventos de temperatura | En Etapa 1 Lilyth se dispara **solo con los eventos del canal de nivel** (los de temperatura gobiernan el estado de calor y su fin). La **Segunda Familia (Restitucion de Temperatura: cuanto tarda $\log \hat{\beta}$ en volver tras un `CAMBIO_TEMPERATURA`)** queda registrada como pendiente para la version posterior | Comandante, Code | separa limpiamente la restitucion de nivel (con el calor como modulador) de la relajacion termica pura |
| E1.16 | D11: Lilyth no emite pruebas | el **resumen** de Lilyth contrasta cuatro cosas como una familia por `holm.hpp`: Delta_tau_iid, Delta_tau_Onsager, sobrepaso contra el nulo sin memoria y sobrepaso contra el nulo lineal. Las fichas por evento siguen siendo medidas, no pruebas | Code | la estructura del codigo (seccion 9) ya los contrasta: sin Holm serian cuatro pruebas sueltas |
| E1.17 | la purga escrita como (y - mu0) sqrt(v_media / v_fase) + mu0 | **la purga real** (`PurgaDiurna`, `escalera_diadica.hpp`) es **y sqrt(v_media / v_fase)**, con el **segundo momento crudo** por fase y sin restar mu0. Por tanto la "marea" que se contrasta en E1.14 es la estacionalidad **multiplicativa** (en observables positivos, como el volumen, arrastra tambien el patron de la media; en retornos, solo el de la varianza). La purga necesita `dias_minimos` (20) dias por fase antes de ser valida: hasta entonces la rama PURGADA no mide. En NIVEL, un incremento purgado invalido rompe la reintegracion: la ranura queda censurada por fin de dato y se declara, nunca se rellena | Code (leido en el codigo) | el README de la Enmienda 1 describia una purga que no es la del operador delta |
| E1.18 | `nulo_sin_memoria` admitia BLOQUES | solo **PERMUTACION** o **SIGNOS**; BLOQUES se rechaza como NINGUNO (E1.2) | Code | BLOQUES conserva la memoria dentro de cada bloque: no es un nulo sin memoria |
| E1.19 | `FASE_ALEATORIZADA_FFT` como nulo lineal | vale **solo para el regimen (a)**. Si se elige, la via (b) no puede usarlo (una serie fija no se propaga, P8): cae al modelo generativo (`AR_LEVINSON_DURBIN`, o `ARFIMA_WHITTLE` si se eligio memoria larga) y la ficha lo marca. El orden p por Quenouille (ultimo kappa_p significativo dentro del K de la serie que mide Lilyth; corregido en E1.22: no el `rezago_corte` de `FichaNulo`) prueba K rezagos a alfa cada uno y tiende a elegir un p de mas: es el error conservador para un nulo lineal y queda declarado | Code | sin la caida explicita, (b) correria sobre una serie que no responde |
| E1.20 | `PRONY_FASE_REESCALADA` sin restriccion | divide por el `tau_rest` de cada evento, asi que excluye por construccion los eventos que no restituyen (tau = inf): la ficha declara cuantos quedaron fuera | Code | es una seleccion; sin declararla, la curva reescalada se lee como de toda la poblacion |
| E1.21 | toda apertura del operador delta con una ranura viva era una replica | una apertura **del mismo signo** con una ranura viva que **no ha restituido** es una **continuacion** del mismo choque: no corta ni abre ranura, y se cuenta. Es replica (corta, competing risk) solo si el signo es contrario o si la ranura ya restituyo | Comandante (2026-11-05), propuesta de Code | con un desplazamiento sostenido el GLR vuelve a alarmar respecto de mu0 y el operador delta reabre eventos sobre el MISMO choque (16 de 23 aperturas en el mundo de pulsos de la revision de la Etapa 1): tratarlas como replicas crea artefactos (anclas dentro del choque, restituciones de 1-3 pasos) |
| E1.22 | el orden p del AR del nulo lineal acotado por "el `rezago_corte` ya medido en `FichaNulo`" (E1.19, 8.2, texto de Code) | (1) **el modelo se ajusta sobre la serie que Lilyth mide, en pasos de la barra base**, con su propio K por `est::tau_integrado` sobre esa serie. El `rezago_corte` de `FichaNulo` no sirve: se mide sobre la serie de la escala j0, en pasos de esa escala (`canal_nivel.hpp`, `calibrar`). (2) **Un modelo por representacion**: con `AMBAS_CONTRASTE_MAREA`, uno sobre la PURGADA y otro sobre la CRUDA. (3) **FLUJO modela y tal cual; NIVEL modela los incrementos dY y reintegra**: si los incrementos son estacionarios el nivel tiene raiz unitaria y aplica P10 (T_max = W_ancla). (4) Las innovaciones de la simulacion son un input, `InnovacionesNulo` (`GAUSS` o `RESIDUOS_REMUESTREADOS`, este por omision: con colas pesadas las gaussianas estrechan el nulo) | Code (2026-11-06), aprobado por el Comandante | un K medido en pasos de la escala j0 usado como tope de un AR en pasos de la barra base es un error de unidades, de la familia del R4 de la Etapa 1 (`det.L` por N_REG); y con un solo modelo, (b) contrastaria una de las dos ramas contra un nulo que no es el suyo |
| E1.23 | el nulo sin memoria, cuando la excedencia (a) no llegaba a m_min, se rellenaba con eventos fabricados sin pasar por el detector (relleno con g < g_k, de cualquier escala o signo, o con 1.0; y el armonico los "enriquecia" hasta 64) | **via (b) sin memoria**: el perfil real, en el espacio de la representacion, se escribe en un sustituto sin memoria (permutado o de signos) de la calma y corre el MISMO detector, con la misma parada ceil(m_min/p_det) y la misma guarda. Si tampoco llega a m_min, el evento queda **sin nulo** y la ficha lo dice. Nada se fabrica | Comandante (2026-11-07, opcion (b)), propuesta de Code en `revision_etapas2_6/REVISION_ETAPAS2_6.md` (D1) | un evento fabricado no pasa por el detector, asi que no tiene el sesgo de seleccion de los reales: el contraste se volvia una comparacion contra un nulo que no es el suyo |
| E1.24 | Dos referencias: sin memoria y Onsager | un **tercer nulo, lineal con replicas (L+H)**: el modelo de la Etapa 2 con un Hawkes de **dos tipos (mismo signo y signo contrario)**, nucleo cruzado con las constantes de tiempo de la escalera del operador delta, escala como marca, tasa base por bin de fase; hijas de signo contrario cortan como replicas. El ajuste tiene que **superar su nulo** (Poisson por fase, razon de verosimilitud con bootstrap parametrico): si no, L+H es **NO DECIDIBLE**; rho(N) >= 1 es NO IDENTIFICADO. **Entra en la familia de Holm (de 4 a 6)**; si no se decide, sus dos p entran con p = 1 marcados (la familia no encoge). Delta_tau_Onsager - Delta_tau_LH mide la autoexcitacion | Comandante (2026-11-08: H.1 si, H.2 los dos signos, H.3 agrupado, H.4 NO DECIDIBLE si no pasa sus nulos), propuesta de Code; directriz para Argon en `MENSAJE_ARGON_NULO_CON_REPLICAS.md` | con replicas de Hawkes (S1, n = 0.6) Onsager rechaza con p 0.005: las hijas del mismo signo pasan a continuacion (E1.21), llegan despues del perfil inyectado y el nulo de Onsager no las tiene |

### 8.2 Como se generan los eventos del nulo (aprobado: (a) + (b))

- **(a) Espontaneos.** B realizaciones del sustituto, con el detector de **produccion** (umbrales congelados, sin recalibrar). Para cada
  evento real k se usan los eventos nulos con **g >= g_k** (excedencia: no hay bins). Vale (a) si hay al menos m_min(alfa) = 8; B por la
  regla de `calibrador_familia.hpp` con tope B_max.
- **(b) Inyeccion generativa** (cuando (a) no llega a m_min). **Correccion de Code al cierre de Argon:** superponer el perfil del choque
  en un sustituto ya calculado no se propaga (la serie no responde): en FLUJO daria restitucion instantanea y en NIVEL un desplazamiento
  permanente, sea cual sea la memoria. Por eso:
  - nulo sin memoria: el perfil real se escribe en un sustituto permutado o de signos y corre el mismo detector (E1.23, 2026-11-07); con menos de m_min, el evento queda sin nulo;
  - nulo lineal: **modelo generativo** configurable (`ModoNuloLineal`, E1.13). AR(p) ajustado a la ACF de la calma por Levinson-Durbin
    (con el orden $p \le K_{\text{corte}}$ seleccionado automaticamente sin numeros a mano por el ultimo coeficiente de reflexion parcial
    significativo $|\kappa_p| > z_{1-\alpha/2}/\sqrt{n}$ de Quenouille dentro del K medido con `est::tau_integrado` sobre la misma serie
    que Lilyth mide, en pasos de la barra base (E1.22; no el `rezago_corte` de `FichaNulo`, que va en pasos de la escala j0), o
    ARFIMA(p, d, 0) con la $d$ del Whittle exacto si el observable tiene memoria larga). Un modelo por representacion; en NIVEL se
    modelan los incrementos y se reintegra (E1.22). Se convierte el perfil real de los $l^*_k$ pasos en sus
    innovaciones (filtro inverso $A(L)$), se inyectan en una realizacion del modelo en tiempos al azar y se deja correr el filtro sin
    inyectar nada mas. Se marca `nulo_por_inyeccion = true`.
  - (a) puede usar el sustituto de fases aleatorizadas o el mismo modelo; (b) solo el modelo (E1.19).

### 8.3 Especificacion de las tres directivas nuevas del Comandante (E1.13, E1.14, E1.15)

#### 1. Nulos y detector de armonicos como inputs del pipeline (E1.13)

Ningun nulo ni criterio armonico queda cableado como constante interna:
- `ConfiguracionLilyth.nulo_sin_memoria`: `NuloSustituto` (`NULO_PERMUTACION` o `NULO_SIGNOS`, segun la tabla de `README_NULOS.md`;
  `NULO_NINGUNO` se rechaza por E1.2 y `NULO_BLOQUES` por E1.18).
- `ConfiguracionLilyth.modo_nulo_lineal`: `AR_LEVINSON_DURBIN` (por omision, orden $p$ derivado por Quenouille hasta el K de la serie medida, E1.22),
  `ARFIMA_WHITTLE` (para observables con memoria larga fraccionaria) o `FASE_ALEATORIZADA_FFT` (solo para el regimen espontaneo (a);
  en (b) cae al modelo generativo, E1.19).
- `ConfiguracionLilyth.modo_armonico`:
  - `PRONY_CURVA_MEDIA`: recurrencia AR(2) por Prony sobre $\bar{\eta}(\tau)$ en tiempo de barra $\tau \in [2, T_{\max}]$ (`a1^2 + 4 a2 < 0` y $\min \bar{\eta} < 0$ contra el nulo).
  - `PRONY_FASE_REESCALADA`: Prony sobre la curva media alineada por la semivida de cada evento ($\tilde{\tau} = \tau / \tau_{\text{rest},k}$), util cuando la frecuencia de oscilacion cambia de un episodio a otro y promediar en $\tau$ absoluto cancelaria fases. Excluye los eventos sin restitucion y la ficha los cuenta (E1.20).
  - `SOBREPASO_POR_EVENTO`: contraste de la distribucion de sobrepasos individuales contra el nulo (para observables discretos o altamente no estacionarios).
  - `TODOS` / `DESACTIVADO`.

#### 2. Medicion dual sobre serie PURGADA y serie CRUDA para contrastar la marea (E1.14)

El operador delta detecta sobre la serie purgada de la estacionalidad diurna (`y_purgado = y * sqrt(v_media / v_fase)`, con v el segundo
momento crudo por fase; corregido en E1.17: la formula anterior con mu0 no es la del codigo).
`ConfiguracionLilyth.representacion` acepta `PURGADA`, `CRUDA` o `AMBAS_CONTRASTE_MAREA` (por omision `AMBAS_CONTRASTE_MAREA`):
- **Rama `PURGADA`:** En `FLUJO` mide directamente sobre $y_{\text{purgado}}(t)$. En `NIVEL` mide sobre el nivel purgado reintegrado
  causalmente paso a paso, $Y_{\text{purgado}}(t) = Y(0) + \sum_{s=1}^t \Delta Y_{\text{purgado}}(s)$, calculando su propia ancla
  $Y^*_{\text{purgada}}$ y su propio $\Delta Y_{\text{purgado}} = Y_{\text{purgado}}(t_a) - Y^*_{\text{purgada}}$. *(Regla de higiene:
  nunca mezclar $Y(t)$ crudo con $Y^*$ purgado).*
- **Rama `CRUDA`:** Para el mismo disparo $t_a$ del detector, congela el ancla $Y^*_{\text{cruda}}$ sobre la serie original sin purgar y
  mide $\eta_{\text{cruda}}(\tau) = (Y_{\text{cruda}}(t_a + \tau) - Y^*_{\text{cruda}}) / \Delta Y_{\text{cruda}}$.
- **Contraste de marea publicado en la ficha:**
  $$\Delta\tau_{\text{marea}} = \tau_{\text{rest}}^{\text{cruda}} - \tau_{\text{rest}}^{\text{purgada}}, \qquad \Delta\text{RMST}_{\text{marea}} = \text{RMST}_{\text{cruda}}(T_{\max}) - \text{RMST}_{\text{purgada}}(T_{\max})$$
  Si $\Delta\tau_{\text{marea}} \approx 0$, la restitucion del observable es **invariante a la marea diurna**; si difieren significativamente,
  cuantifica cuanto acelera o frena la marea horaria el retorno aparente hacia el ancla.

#### 3. Pendiente declarado para version posterior: Segunda Familia -- Restitucion de Temperatura (E1.15)

En la Etapa 1, Lilyth abre ranuras **unicamente ante eventos del canal de nivel** (`canal == 0`), mientras que los eventos de
`canal_temperatura` (`canal == 1`) actuan como moduladores del estado termico (`E1.6` y `E1.7`).
Queda anotado como **pendiente aprobado para una version posterior** instanciar una **Segunda Familia de Lilyth sobre la Temperatura**:
- Disparada por las aperturas/confirmaciones de `CAMBIO_TEMPERATURA`.
- Congelando como ancla $\log \beta^*$ la mediana causal de $\log \hat{\beta}$ previa al choque termico.
- Midiendo el tiempo de relajacion termica $\tau_{\text{rest},\beta}$, su contraste de Onsager sobre la ACF de $\log \hat{\beta}$ y si
  existen oscilaciones de volatilidad post-choque.

### 8.4 Preexpectativas nuevas o cambiadas (sustituyen a las de la seccion 5 donde coinciden)

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
| P11 | contraste de marea (E1.14) | en un AR(1) multiplicado por perfil diurno determinista: sin purga (`CRUDA`) `tau_rest` depende de la fase horaria del choque; en `PURGADA` recupera la semivida verdadera y $\Delta\tau_{\text{marea}}$ mide el sesgo horario |

| P12 | E1.18 y E1.19 | `nulo_sin_memoria` = BLOQUES se rechaza; `FASE_ALEATORIZADA_FFT` en (b) cae al modelo generativo y la ficha lo marca |

### 8.5 Hallazgos de la Etapa 3 sobre sinteticos (2026-11-07, revision y correccion de Code)

1. **Onsager mezcla no linealidad y autoexcitacion.** Con replicas de Hawkes del mismo signo (S1, n = 0.6) el contraste de Onsager
   rechaza (p 0.005). El Comandante lo confirma: hace falta un nulo con replicas. Decidido el 2026-11-08 como E1.24 (directriz para Argon).
2. **Con memoria basal, el nulo sin memoria sale casi entero de la via (b).** En los mundos AR, el detector calibrado sobre la serie con
   memoria casi no alarma sobre el sustituto permutado, y la via (a) nunca alcanza m_min.
   - El Comandante senala que la literatura mide un Hurst basal de ~0.7. En este proyecto ese nivel aparece en la actividad, el volumen
     y el flujo: volumen 0.79-1.02 (`acta_hurst_oinfo_volumen.txt`) y flujo L3 0.75-1.06 (`Comportamiento de manada.md`). En los
     retornos el Hurst es ~0.5 (ELW de BTC/ETH ~0.47, `MEDIDO_WHITTLE.md`).
   - Consecuencia: si el operador delta mira un observable con memoria (volatilidad, actividad, flujo), en dato real el nulo sin memoria
     se construira sobre todo por la via (b). Es lo esperado, no un fallo. La ficha tiene que reportar que fraccion de los eventos usa
     la via (a) y cual la (b).
   - Con H ~ 0.7 (d ~ 0.2) un AR truncado en K por Quenouille puede quedarse corto de memoria; para esos observables el nulo lineal que
     corresponde es `ARFIMA_WHITTLE` (E1.13), con la d del sensor ELW.
3. **Con dY de una sola observacion, el primer toque restituye por regresion a la media.** El tau_rest crudo no se publica; el exceso
   sobre los nulos si.

### 8.6 Estado de E1.24 (2026-11-11, Code construyendo; NO conectado a la base)

Lo que funciona, en `PRUEBAS_JEFE/lilyth/etapa_lh/`:

- **Inferencia indirecta de las cuatro masas.** Converge en 2-3 iteraciones en todos los mundos (antes, 2 de 10 semillas no convergian
  en 64).
  - La tolerancia de cada masa es 2 x la mayor de dos sd de una realizacion: bajo el estructural y la del propio auxiliar real
    (bootstrap parametrico de si mismo).
  - Un contraste de cada masa contra Onsager se probo y se descarto: Onsager da casi siempre 0 exacto en la cruzada y cualquier valor
    positivo salia con p = 1/(B+1).
- **La cota de rho invierte la ligadura.** La pendiente del auxiliar respecto del estructural es 0.66-0.88, no 1. Se mide por la secante
  entre el estructural y el mismo a media masa, con numeros aleatorios comunes y diferencia pareada resuelta a 2 ee. El punto se corrige
  con un paso de Newton.
  - **Este detector no identifica la masa cruzada** (pendiente no resuelta): queda en [0, inf).
- **Bateria:** pasan R3, R6, R7, R8 y D-ARFIMA; R1 queda en su banda con 12 semillas (2/12).

Lo que NO funciona, y bloquea conectar:

- **El contraste L+H tiene tamano de mas.** Con 12 semillas por mundo de replicas del mismo signo, L+H rechaza 5 de 24 y Onsager 4 de 24
  (esperado ~1.2).
- **El sesgo es comun a los dos nulos y de signo fijo.** S < 0 en 20 de 24: la recuperacion real es mas rapida que la de los nulos. Ya
  asoma sin replicas (R1: S < 0 en 10 de 12, sin rechazos de mas) y crece con n (S medio ~ -0.04 en n = 0.3 y ~ -0.08 en n = 0.6).
- **Las hijas de L+H no lo corrigen.** Sospechas, sin medir:
  1. el perfil inyectado es el del evento real y ya lleva dentro sus replicas, asi que L+H las cuenta dos veces;
  2. el sesgo de base del contraste de Onsager en estos mundos.
- **Sondas (2026-11-11, por orden del Comandante): queda como CURIOSIDAD, no se sigue escarbando.**
  - **(a) Doble conteo: descartado.** Los perfiles inyectados miden 1-4 pasos (118 ranuras, mundo n = 0.6) y las replicas llegan >= 6
    pasos despues de la madre: un perfil no puede llevar madre e hija a la vez (`sonda_a_perfil.cpp` y su log).
  - **(b) Sesgo de base: existe y es pequeno.** Sin replicas, con 40 semillas, Onsager rechaza 2/40 (tamano correcto) y R1 pasa, pero
    S medio = -0.023, ee 0.0097 (t ~ -2.4). Las replicas lo agrandan y las hijas de L+H no lo corrigen (`sonda_b_r1_s40.log`).
  - Lectura del Comandante: probablemente faltan conceptos por conectar, no un error de codigo.

## 9. Estructura del codigo (propuesta por Code, aprobada por el Comandante el 2026-11-05)

### 9.1 Como se engancha sin tocar nada auditado

La `Salida` de `ObservadorSaltos` ya expone en cada paso las alarmas de cada escala (`nivel[j].alarma`: t_cambio, signo, l*, g) y las
de temperatura (`temp[j].alarma`, `log_beta`). Lilyth lee todo desde fuera:
- **apertura:** una alarma de nivel en un paso en que el operador delta no tenia evento abierto (`evento_abierto` falso en t - 1): es t_a,
  con su ficha;
- **replica:** una alarma de nivel mientras una ranura sigue viva;
- **calor:** las alarmas de `temp[j]` y su `log_beta`, que alimenta el fin del calor (E1.7).

`ObservadorSaltos` no se modifica. Control obligatorio: con Lilyth y sin Lilyth, las salidas del operador delta son identicas bit a bit.

### 9.2 Archivos (en `ark_cxx/base_auditada/include/ark/lilyth/`)

| archivo | que hace | de que se alimenta |
|---|---|---|
| `tipos_lilyth.hpp` | `ConfiguracionLilyth` (E1.13, E1.14), `RanuraLilyth` (POD: copia de la apertura, Y* y dY por representacion, t_a, j0, W_ancla, guarda, nu(t_a), nu_punto*, eta hasta T_MAX, medianas de la Capa 1), `FichaLilyth` por evento, `ResumenLilyth` | `tipos_metrologia.hpp` |
| `trayectoria_lilyth.hpp` | eta(tau) con signo, `tau_rest_primer_toque` y `_consolidado`, v_rel, tau_act | -- |
| `riesgos_competitivos.hpp` | contadores en_riesgo, restituidos y replicas por tau; CIF de Aalen-Johansen; RMST(T_max) por signo; cero reservas | -- |
| `modelo_lineal_nulo.hpp` | Levinson-Durbin sobre la ACF de calma de la serie que mide Lilyth (un modelo por representacion; NIVEL sobre dY), orden por Quenouille hasta su propio K (E1.19, E1.22), ARFIMA(p, d, 0) opcional, filtro inverso A(L) para la inyeccion y simulacion. Firmas en `PLANO_ETAPA2_MODELO_LINEAL_NULO.md` | `whittle_exacto.hpp` (su `coeficientes(d)` ya da los pesos fraccionarios), `estadistica_metrologia.hpp` |
| `nulo_lilyth.hpp` | (a) B realizaciones con el detector de produccion congelado y excedencia g >= g_k; (b) inyeccion generativa; las dos referencias; T_max (o W_ancla si hay raiz unitaria) | `sustitutos.hpp` (PERMUTACION, SIGNOS), `calibrador_familia.hpp` (m_min, B, `en_paralelo`), una copia del `ObservadorSaltos` con los umbrales de produccion |
| `armonico_lilyth.hpp` | eta_barra(tau), Prony AR(2), raices complejas y su periodo, min eta_barra contra los dos nulos, etiqueta `ARMONICO_SOSPECHOSO`, segun `modo_armonico` | -- |
| `fin_calor_lilyth.hpp` | detector unilateral hacia mu_0,beta sobre el nulo caliente (mu_calor de las ventanas posteriores), ARL_0,caliente = T_max / alfa | la ficha de `canal_temperatura` (mu0, s0 y tau_int de log beta); la busqueda de h de `canal_nivel` |
| `observador_lilyth.hpp` | el orquestador (9.3) | todo lo anterior y `holm.hpp` (E1.16) |
| `sensor_eka_lilyth.hpp` | Capa 1, la ultima | flujo de filas |

### 9.3 El orquestador

`ObservadorLilyth<Obs, RANURAS, T_MAX>`, con `Obs` el `ObservadorSaltos` concreto:
- `calibrar(const Obs& delta, cfg, historia, den, valido, fase, dia, nu)` construye los dos nulos, T_max y el modelo lineal a partir de las
  fichas del operador delta ya calibrado. No recalibra al operador delta.
- `paso(const Obs::Salida& s, double y, double nu)`: ciclo caliente, cero reservas.
- RANURAS y T_MAX son **topes de memoria** (como W_MAX), no parametros: el desborde se cuenta y se declara; si el T_max derivado supera
  T_MAX, la ficha lo marca.
- El resumen pasa sus cuatro contrastes por `holm.hpp` (E1.16).
- **Por filas:** el camino de `ObservadorSaltosFilas` (`flujo_barras.hpp`) ya entrega el nu de cada barra (NU_TRADES / NU_VOLUMEN).
  Lilyth guarda tambien la marca de pared de cada barra: en barras intrinsecas la pared no avanza a paso fijo (P4).
- **Capa 1:** el operador delta detecta con L_Eka de medianas moviles; Lilyth mide cada ranura con sus medianas congeladas. El sensor
  expone su estado de medianas y una funcion L(D, s, medianas); la ranura recalcula Y con su copia (E1.12).

### 9.4 Orden de construccion (cada etapa con sus controles y mutantes en `PRUEBAS_JEFE/lilyth/`)

1. **Tipos, trayectoria, ranuras y riesgos competitivos** (sin nulos). Controles: identidad del operador delta con y sin Lilyth;
   causalidad por truncado y perturbacion del futuro (con el mutante de un paso, que tiene que fallar); eta y los dos tau sobre casos a
   mano; Aalen-Johansen contra una referencia escrita aparte, con replicas que cortan trayectorias; desborde de ranuras contado; cero
   reservas en `paso`; las dos ramas de representacion (E1.14, E1.17) con la purga real.
2. **Modelo lineal nulo:** Levinson contra un AR conocido; el filtro inverso reconstruye las innovaciones; P8.
3. **Nulos (a) y (b) conectados al operador delta:** P1, P2, P2b, P3a, P3b, P10, P12.
4. **Armonico:** P5, P5b.
5. **Fin del calor:** P9.
6. **Reloj doble y marea:** P4, P11.
7. **Capa 1:** P6.
8. Mover el README a la base, regenerar el manifiesto y las copias de los dos agentes.

Reparto acordado: Argon construye cada etapa; Code la verifica (controles propios, mutantes sobre su codigo, `verificar_base.py`) y la
expande antes de pasar a la siguiente.
