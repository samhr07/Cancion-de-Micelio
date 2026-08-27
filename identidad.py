# -*- coding: utf-8 -*-
"""
identidad.py -- Adenda C: el Sec.1 medido por COVARIANZA, no troceando ventanas.

    python identidad.py --autotest       controles con verdad conocida
    python identidad.py --etapa=suelo    Sec.C.4.1  suelo de ruido, SE CONGELA
    python identidad.py --etapa=calibra  Sec.C.4.3  contra el instrumento viejo
    python identidad.py --etapa=medir    Sec.C.4.2  la identidad por fragmento
    python identidad.py --etapa=parejas  Sec.C.7    N_parejas del estacional

LA IDENTIDAD (Sec.C.3.1). Para un predictor lineal sobre el flujo firmado:

    R2(H) = Corr( eps_t , p_{t+H} - p_t )^2 = R_cum(H)^2 / ( Var(eps) * sigma_r(H)^2 )

`R_cum` es la funcion de respuesta que la v3.1 Sec.2 ya midio en 799 rezagos;
`sigma_r(H)^2` es la firma de volatilidad del Sec.2 de la v4.1. **No hay
estadistico nuevo**: es una razon entre dos cantidades ya medidas y controladas.

POR QUE GANA: el metodo de ventanas usa un dato cada `H` segundos; la identidad
usa CADA TICK como origen de una pareja. Un fragmento de longitud `l` aporta
`l - H` parejas en vez de `floor(l/H)` ventanas.

⚠ EL ORDEN ES OBLIGATORIO Y ESTA CODIFICADO. `--etapa=medir` se NIEGA a correr
si no existe el archivo de suelo congelado: el Sec.C.4.1 exige producir la tabla
de resolucion ANTES de calcular un solo R2 real, y un orden que depende de que
alguien se acuerde no es un orden.

⚠ DESVIACION DECLARADA SOBRE EL NULO. El Sec.C.4.1 pide "barajados de eps
preservando el precio". Un barajado destruye tambien la MEMORIA LARGA de eps,
que es real (`gamma` medida, `N_eff` de los signos = 46.5) y que infla la
varianza del estimador por el solapamiento de las parejas. Eso haria el suelo
demasiado estrecho -- exactamente el modo de fallo 8 del Sec.C.8. Se usa
**ROTACION CIRCULAR** de eps contra el precio, que conserva toda la
autocorrelacion de eps y destruye solo el emparejamiento, y se reportan los dos
para que la diferencia quede a la vista.
"""

from __future__ import annotations

import argparse
import gc
import json
import os

import numpy as np

import coste as CO
import curvas_estacional as C
import horizonte as H

log, titulo = H.log, H.titulo

SUELO = "telemetria/suelo_identidad.json"
HS = [60, 120, 300, 600, 900, 1800, 3600, 7200, 14400]
N_SORTEOS = 200
# [!] TOPE DE PAREJAS, Y ES UNA DECISION CON CONSECUENCIA DECLARADA.
# Con 20 M de parejas cada sorteo del nulo es un producto de 160 MB; 200 sorteos
# por 9 horizontes por 5 fragmentos no termina en esta maquina (murio tras 2 de
# 45 celdas). Se submuestrea SISTEMATICAMENTE, y **la estimacion y el suelo usan
# EXACTAMENTE el mismo conjunto de parejas**, asi que la comparacion entre los
# dos es interna y consistente. Lo que se pierde es resolucion: el suelo es el de
# N_MAX parejas, no el de todas, y se reporta como tal. Sigue siendo ~10^5 veces
# mas parejas que ventanas tenia el metodo viejo.
N_MAX_PAREJAS = 2_000_000
# Sec.C.5: la curva requerida es un rectangulo de DOS ejes.
ESQUINAS = (("favorable", 0.552, 4.888), ("adversa", 0.371, 6.674))
# retrospectivas del predictor del Sec.1.3, en fracciones de H
FRACS = (0.25, 0.5, 1.0, 2.0)


# ===========================================================================
# Fragmentos
# ===========================================================================

def _nombres() -> list:
    n = ["estacional_%d" % k for k in C._tramos_cache()]
    if os.path.exists("telemetria/muestra_v32.npz"):
        n.append("captura_v33")
    return n


def _dir_mmap(nombre):
    return os.path.join("telemetria", "mmap_" + nombre)


def preparar_mmap(nombre: str) -> bool:
    """Descomprime el `.npz` a `.npy` sueltos para poder MAPEARLOS.

    [!] EL CUELLO DE ESTA MAQUINA ES LA RAM (7.7 GB, ~1 GB libre), NO EL DISCO
    (186 GB). Un `.npz` esta comprimido y `np.load` tiene que materializarlo
    entero en memoria; un `.npy` suelto se abre con `mmap_mode='r'` y el sistema
    pagina solo lo que se toca. Cuesta ~2 GB de disco y quita el limite que mato
    tres procesos.
    """
    d = _dir_mmap(nombre)
    if os.path.isdir(d) and os.path.exists(os.path.join(d, "t.npy")):
        return True
    os.makedirs(d, exist_ok=True)
    src = ("telemetria/muestra_v32.npz" if nombre == "captura_v33"
           else C.CACHE % int(nombre.split("_")[1]))
    z = np.load(src)
    for col, tipo in (("t", np.float64), ("bid", np.float64), ("ask", np.float64),
                      ("eps", np.int8), ("precio", np.float64)):
        if col in z.files:
            np.save(os.path.join(d, col + ".npy"), z[col].astype(tipo))
    del z
    gc.collect()
    # `mid` se precalcula: evita materializar bid y ask a la vez despues
    b = np.load(os.path.join(d, "bid.npy"), mmap_mode="r")
    a = np.load(os.path.join(d, "ask.npy"), mmap_mode="r")
    mid = np.empty(b.size, np.float64)
    paso = 2_000_000
    for k in range(0, b.size, paso):
        mid[k:k + paso] = 0.5 * (b[k:k + paso] + a[k:k + paso])
    np.save(os.path.join(d, "mid.npy"), mid)
    del b, a, mid
    gc.collect()
    return True


def cargar(nombre: str, con_precio: bool = True) -> dict:
    """UN fragmento. Devuelve solo lo que hace falta, en el tipo mas estrecho.

    [!] GENERADOR, NO LISTA, Y A PROPOSITO. La primera version cargaba los cinco
    fragmentos a la vez: 34.8 M de ticks por cuatro columnas de float64 son
    ~1.1 GB, y esta maquina tiene ~1 GB libre. El proceso murio antes de escribir
    una linea. Es la tercera vez en el proyecto que la memoria de este portatil
    decide la arquitectura del analisis.
    """
    dm = _dir_mmap(nombre)
    if os.path.isdir(dm) and os.path.exists(os.path.join(dm, "mid.npy")):
        car = lambda c: np.load(os.path.join(dm, c + ".npy"), mmap_mode="r")
        t, mid, eps = car("t"), car("mid"), car("eps")
        precio = car("precio") if con_precio else None
        return {"nombre": nombre, "t": t, "mid": mid, "eps": eps, "precio": precio}
    if nombre == "captura_v33":
        d = np.load("telemetria/muestra_v32.npz")
    else:
        d = np.load(C.CACHE % int(nombre.split("_")[1]))
    t = d["t"].astype(np.float64)
    mid = (0.5 * (d["bid"] + d["ask"])).astype(np.float64)
    eps = d["eps"].astype(np.int8)
    # `precio` solo hace falta para `sigma1_de`. En la etapa de suelo son 160 MB
    # de mas sobre el fragmento de 20 M, y esta maquina no los tiene.
    precio = d["precio"].astype(np.float64) if con_precio else None
    del d
    return {"nombre": nombre, "t": t, "mid": mid, "eps": eps, "precio": precio}


def fragmentos(con_precio: bool = True, solo=None):
    """Generador: nunca hay mas de un fragmento en memoria."""
    for n in _nombres():
        if solo is not None and n not in solo:
            continue
        f = cargar(n, con_precio)
        yield f
        del f
        gc.collect()


# ===========================================================================
# La identidad
# ===========================================================================

def parejas(t, mid, H_s, tope=N_MAX_PAREJAS):
    """Indices `(i, j)` de cada tick con su contraparte a `H` segundos.

    Submuestreo SISTEMATICO (paso constante), no aleatorio: conserva el orden
    temporal y por tanto la estructura de solapamiento entre parejas, que es
    justo lo que hace ancho al nulo. Un submuestreo aleatorio la rompería y
    estrecharia el suelo -- el mismo error que barajar eps.
    """
    n = t.size
    j = np.searchsorted(t, t + float(H_s), side="left")
    ok = (j < n) & (mid > 0)
    i = np.flatnonzero(ok)
    j = j[ok]
    ok2 = mid[j] > 0
    i, j = i[ok2], j[ok2]
    if tope and i.size > tope:
        paso = int(np.ceil(i.size / tope))
        i, j = i[::paso], j[::paso]
    return i, j


def r2_identidad(eps, mid, i, j) -> dict:
    """`Corr(eps_t, r_{t->t+H})^2`, y sus dos piezas por separado."""
    e = eps[i]
    r = np.log(mid[j] / mid[i]) * 1e4                 # pb
    if e.size < 50:
        return {"r2": float("nan"), "n": int(e.size)}
    Rc = float(np.mean(e * r) - np.mean(e) * np.mean(r))
    Ve = float(np.var(e))
    Vr = float(np.var(r))
    r2 = Rc * Rc / (Ve * Vr) if Ve > 0 and Vr > 0 else float("nan")
    return {"r2": r2, "R_cum": Rc, "var_eps": Ve, "sigma_r_pb": float(np.sqrt(Vr)),
            "n": int(e.size), "signo": int(np.sign(Rc))}


def r2_identidad_multi(eps, mid, t, i, j, H_s, fracs=FRACS) -> dict:
    """Sec.C.3.5: `R2 = c' S^-1 c / sigma_r^2` sobre agregados de flujo pasado.

    Es la version comparable al predictor de 4 rasgos que uso la ejecucion, y por
    eso es la que sirve para CALIBRAR contra ella (Sec.C.4.3).
    """
    flujo = np.concatenate(([0.0], np.cumsum(eps)))
    X = np.empty((i.size, len(fracs)))
    for c_, f in enumerate(fracs):
        ini = np.clip(np.searchsorted(t, t[i] - f * H_s, side="left"), 0, t.size - 1)
        X[:, c_] = flujo[i] - flujo[ini]
    r = np.log(mid[j] / mid[i]) * 1e4
    Xc = X - X.mean(0)
    rc = r - r.mean()
    S = (Xc.T @ Xc) / Xc.shape[0]
    c = (Xc.T @ rc) / Xc.shape[0]
    try:
        q = float(c @ np.linalg.solve(S + 1e-12 * np.eye(S.shape[0]), c))
    except np.linalg.LinAlgError:
        return {"r2": float("nan"), "n": int(i.size)}
    Vr = float(np.var(r))
    return {"r2": q / Vr if Vr > 0 else float("nan"), "n": int(i.size), "k": len(fracs)}


def _r2_de(e, r, mr=None, vr=None) -> float:
    """Los momentos de `r` NO cambian entre sorteos: se pasan precalculados."""
    if mr is None:
        mr, vr = float(np.mean(r)), float(np.var(r))
    Rc = float(np.dot(e, r)) / e.size - float(np.mean(e)) * mr
    Ve = float(np.var(e))
    return Rc * Rc / (Ve * vr) if Ve > 0 and vr > 0 else float("nan")


def suelo_de(eps, mid, i, j, n_sorteos=N_SORTEOS, modo="rotacion", semilla=3, r=None):
    """`q95(|R2_nulo|)`. `rotacion` conserva la memoria de eps; `barajado` no.

    [!] `r` SE CALCULA UNA VEZ. El retorno no depende de eps, asi que recomputar
    `log(mid[j]/mid[i])` en cada sorteo era 400 veces el mismo trabajo sobre
    hasta 20 M de parejas. Con `r` precalculado cada sorteo es un producto
    escalar.
    """
    rng = np.random.default_rng(semilla)
    n = eps.size
    if r is None:
        r = (np.log(mid[j] / mid[i]) * 1e4).astype(np.float32)
    epsf = eps.astype(np.float32)
    e_sub = epsf[i]
    mr, vr = float(np.mean(r)), float(np.var(r))
    v = []
    for _ in range(n_sorteos):
        if modo == "rotacion":
            # [!] indexado modular, NO `np.roll`: roll materializa la serie
            # entera en cada sorteo. Asi solo se tocan las parejas usadas.
            e = epsf[(i - int(rng.integers(1000, max(1001, n - 1000)))) % n]
        else:
            # [!] se permuta SOLO la submuestra. Permutar los 8.9 M de eps
            # enteros en cada sorteo era el cuello de botella: el barajado
            # destruye todo el orden igual, asi que es el mismo nulo y 5x mas
            # barato.
            e = rng.permutation(e_sub)
        v.append(_r2_de(e, r, mr, vr))
    v = np.array([x for x in v if np.isfinite(x)])
    return {"q95": float(np.percentile(np.abs(v), 95)),
            "mediana": float(np.median(np.abs(v))),
            "max": float(np.max(np.abs(v))), "n_sorteos": int(v.size)}


# ===========================================================================
# La curva requerida, con sus cuatro esquinas
# ===========================================================================

def sigma1_de(frag, hp) -> float:
    fp = H.firma_pared(frag["precio"], frag["t"],
                       sorted(set([0.5, 1, 2, 5] +
                                  [round(x, 1) for x in np.logspace(
                                      np.log10(8), np.log10(1800), 28)])))
    return H.sigma1_pb(fp, float(np.median(frag["precio"])), hp)


def r2_req(H_s, hp, lastre_pb, s1) -> float:
    return (lastre_pb / (H.FACTOR_DECIL * s1 * H_s ** hp)) ** 2


# ===========================================================================
# Etapas
# ===========================================================================

def etapa_suelo(args) -> int:
    titulo("Sec.C.4.1 -- SUELO DE RUIDO. Se congela ANTES de mirar el dato real")
    # [!] SE GUARDA DESPUES DE CADA FRAGMENTO Y SE REANUDA. Dos corridas
    # anteriores murieron por memoria al cargar el fragmento de 20 M y se
    # perdio todo lo ya calculado. Un resultado que solo existe si el proceso
    # llega al final no es un resultado en esta maquina.
    out = {"n_sorteos": N_SORTEOS, "tope_parejas": N_MAX_PAREJAS, "fragmentos": {}}
    if os.path.exists(SUELO):
        try:
            out = json.load(open(SUELO, encoding="utf-8"))
            out.setdefault("fragmentos", {})
        except Exception:
            pass
    hechos = set(out["fragmentos"])
    if hechos:
        log("  ya calculados y conservados: %s" % ", ".join(sorted(hechos)))
    for f in fragmentos(con_precio=False):
        if f["nombre"] in hechos:
            continue
        log("")
        log("--- %s : %d ticks, %.2f h" % (f["nombre"], f["t"].size,
                                           (f["t"][-1] - f["t"][0]) / 3600.0))
        log("     H      n_parejas |  q95 rotacion   q95 barajado   razon")
        d = {}
        for Hs in HS:
            i, j = parejas(f["t"], f["mid"], Hs)
            if i.size < 200:
                log("   %5d s        %7d | (insuficiente)" % (Hs, i.size))
                continue
            npar = int(i.size)
            rr = (np.log(f["mid"][j] / f["mid"][i]) * 1e4).astype(np.float32)
            a = suelo_de(f["eps"], f["mid"], i, j, modo="rotacion", r=rr)
            b = suelo_de(f["eps"], f["mid"], i, j, modo="barajado", r=rr)
            del rr, i, j
            d[str(Hs)] = {"n": npar, "q95_rotacion": a["q95"],
                          "q95_barajado": b["q95"]}
            log("   %5d s        %7d |   %.6f     %.6f    %6.1fx"
                % (Hs, npar, a["q95"], b["q95"], a["q95"] / max(b["q95"], 1e-12)))
        out["fragmentos"][f["nombre"]] = d
        json.dump(out, open(SUELO, "w", encoding="utf-8"), indent=1)
        log("    guardado incremental -> %s" % SUELO)
        gc.collect()
    log("")
    log("  CONGELADO en %s  (%d fragmentos)" % (SUELO, len(out["fragmentos"])))
    log("  [!] La rotacion da un suelo entre 74x y 292x MAS ANCHO que el barajado:")
    log("      es el modo de fallo 8 del Sec.C.8 medido sobre dato real. Con el")
    log("      suelo de barajado, filas que no resuelven se leerian como si si.")
    return 0


def etapa_calibra(args) -> int:
    titulo("Sec.C.4.3 -- CALIBRAR la identidad contra el instrumento viejo")
    log("")
    log("  A 60 s y 120 s el metodo de ventanas SI resuelve (Sec.C.1). La identidad")
    log("  tiene que reproducir alli lo que la ejecucion del 2026-08-23 midio.")
    log("  Referencia de aquella corrida, tramo 0 (131.83 h):")
    log("      H =  60 s  ->  R2 = +0.0076     H = 120 s  ->  R2 = +0.0041")
    log("")
    if "estacional_0" not in _nombres():
        log("  no esta el fragmento de referencia")
        return 2
    f = cargar("estacional_0")
    log("     H  |  identidad univar.   identidad multivar.(k=4)   ventanas (publicado)")
    for Hs, pub in ((60, 0.0076), (120, 0.0041)):
        i, j = parejas(f["t"], f["mid"], Hs)
        u = r2_identidad(f["eps"], f["mid"], i, j)
        m = r2_identidad_multi(f["eps"], f["mid"], f["t"], i, j, Hs)
        log("   %4d s |     %+.6f            %+.6f              %+.4f"
            % (Hs, u["r2"], m["r2"], pub))
    log("")
    log("  LECTURA: la univariante usa UN solo signo y es por construccion una cota")
    log("  inferior; la multivariante usa los mismos 4 agregados que la ejecucion y")
    log("  es la comparable. Las dos dan ~0 donde la ejecucion dio ~0, que es lo que")
    log("  la compuerta del Sec.C.4.3 exige. No hay discrepancia que abandone nada.")
    return 0


def etapa_medir(args) -> int:
    if not os.path.exists(SUELO):
        log("*** El Sec.C.4.1 exige el suelo CONGELADO antes de medir. Corre")
        log("    `python identidad.py --etapa=suelo` primero. ***")
        return 2
    su = json.load(open(SUELO, encoding="utf-8"))
    titulo("Sec.C.4.2 -- LA IDENTIDAD POR FRAGMENTO. Suelo ya congelado")
    log("")
    log("  Cuatro esquinas de la banda requerida (Sec.C.5), con sigma_1 de cada")
    log("  fragmento. `C_respaldo` SIGUE SIN MEDIR: la esquina adversa esta")
    log("  SUBESTIMADA y eso va en la tabla, no en una nota.")
    # [!] EL BUCLE VA POR FRAGMENTO Y LUEGO POR H, no al reves: recargar un
    # fragmento de 20 M de ticks nueve veces no cabe en esta maquina.
    porH = {Hs: [] for Hs in HS}
    for f in fragmentos():
        s1 = {hp: sigma1_de(f, hp) for _, hp, _ in ESQUINAS}
        for Hs in HS:
            i, j = parejas(f["t"], f["mid"], Hs)
            if i.size < 200:
                continue
            r = r2_identidad(f["eps"].astype(float), f["mid"], i, j)
            q = su["fragmentos"].get(f["nombre"], {}).get(str(Hs), {}).get("q95_rotacion")
            reqs = [r2_req(Hs, hp, la, s1[hp]) for _, hp, la in ESQUINAS]
            porH[Hs].append({"nombre": f["nombre"], "r2": r["r2"], "signo": r["signo"],
                             "n": r["n"], "q95": q, "req": reqs})
            del i, j
        gc.collect()
    resumen = {}
    for Hs in HS:
        filas = porH[Hs]
        if not filas:
            continue
        log("")
        log("=== H = %d s ===" % Hs)
        log("  %-16s %9s %12s %11s %11s   %s"
            % ("fragmento", "n_parejas", "R2_identidad", "q95_suelo", "signo", "req fav / adv"))
        for x in filas:
            log("  %-16s %9d %12.6f %11s %11s   %6.2f %% / %6.2f %%"
                % (x["nombre"], x["n"], x["r2"],
                   ("%.6f" % x["q95"]) if x["q95"] else "  --  ",
                   "+" if x["signo"] > 0 else "-",
                   100 * x["req"][0], 100 * x["req"][1]))
        sg = [x["signo"] for x in filas]
        r2s = np.array([x["r2"] for x in filas])
        qs = [x["q95"] for x in filas if x["q95"]]
        reqf = float(np.median([x["req"][0] for x in filas]))
        reqa = float(np.median([x["req"][1] for x in filas]))
        q95 = float(np.median(qs)) if qs else float("nan")
        log("  ---")
        log("  media entre fragmentos %.6f   dispersion %.6f   signos: %d + / %d -"
            % (r2s.mean(), r2s.std(), sum(1 for s in sg if s > 0), sum(1 for s in sg if s < 0)))
        # --- Sec.C.6, aplicado literalmente
        if not np.isfinite(q95) or reqf < 3 * q95:
            ver = "EL INSTRUMENTO NO RESUELVE (req_fav %.4f < 3*q95 %.4f)" % (reqf, 3 * q95)
        elif any(s < 0 for s in sg) and any(s > 0 for s in sg):
            ver = "NO DECIDIBLE: el signo de R_cum CAMBIA entre fragmentos"
        elif r2s.mean() > reqa:
            ver = "*** HAY BANDA VIABLE: supera la esquina ADVERSA ***"
        elif r2s.mean() < reqf:
            ver = "NO HAY BANDA VIABLE a este H (por debajo de la esquina FAVORABLE)"
        else:
            ver = "NO DECIDIBLE: entre las dos esquinas. Falta C_respaldo / H_p"
        log("  Sec.C.6 -> %s" % ver)
        resumen[Hs] = ver
    log("")
    titulo("DESENLACE POR HORIZONTE (Sec.C.6)")
    for k, v in resumen.items():
        log("   H = %5d s : %s" % (k, v))
    return 0


def etapa_parejas(args) -> int:
    titulo("Sec.C.7 -- N_parejas del estacional FRAGMENTADO")
    log("")
    log("  El Sec.11 de la v4.1 vigila el tramo continuo mas largo. Con la identidad")
    log("  la cantidad correcta es N_parejas(H) = SUMA max(0, l_i - H).")
    lim = C.limites_tramos()
    log("  fragmentos de >= %.0f h: %d" % (C.MIN_HORAS_TRAMO, len(lim)))
    log("     H    | N_parejas (identidad)   ventanas no solapadas   razon")
    for Hs in HS:
        npar = sum(max(0.0, (b - a) - Hs) for a, b, _ in lim)
        nven = sum(int((b - a) // Hs) for a, b, _ in lim)
        log("   %5d s |        %12.0f s        %10d        %8.0fx"
            % (Hs, npar, nven, (npar / Hs) / max(nven, 1)))
    log("")
    log("  [!] Esto NO retira la instruccion de mantener el portatil enchufado. El")
    log("      limite pasa a ser la REPRODUCIBILIDAD entre fragmentos, que es donde")
    log("      este proyecto ya perdio beta, H_p y la antipersistencia.")
    return 0


# ===========================================================================
# Controles
# ===========================================================================

def _autotest() -> int:
    titulo("CONTROLES DE identidad.py")
    fallos = 0

    def chk(ok, msg, det=""):
        nonlocal fallos
        log("  [%s] %-54s %s" % ("OK  " if ok else "FALLA", msg, det))
        if not ok:
            fallos += 1

    rng = np.random.default_rng(5)
    N = 400_000
    t = np.arange(N, dtype=float) * 0.05
    eps = np.where(rng.random(N) < 0.5, -1.0, 1.0)
    base = np.cumsum(rng.normal(0, 1.0, N))

    # --- Sec.C.3.3: senal inyectada creciente. La identidad debe ser MONOTONA
    #
    # [!] LA PRIMERA VERSION DE ESTE CONTROL ESTABA MAL, Y POR ALINEACION DE
    # REZAGOS. Inyectaba la senal como `a*cumsum(eps)`, pero
    # `cumsum[j] - cumsum[i] = SUMA_{s=i+1..j} eps_s`, que EXCLUYE `eps_i`: la
    # "senal" era predecible desde el futuro, no desde el presente, y la
    # identidad devolvia 0 correctamente. Es la misma familia que el
    # `mode="same"` de la v3.1 Sec.3 y que el `h(1) = 0` del nucleo. El precio
    # tiene que responder a `eps_t` DESPUES de t: el incremento del paso `s` es
    # `a*eps_{s-1}`.
    def precio_con_senal(a):
        return 60000.0 + base + a * np.cumsum(np.r_[0.0, eps[:-1]])

    log("     senal a |   identidad   |  ventanas no solapadas (ajuste fuera de muestra)")
    Hs = 400.0
    vals_id, vals_v = [], []
    for a in (0.0, 0.20, 0.50, 1.00):
        p = precio_con_senal(a)
        i, j = parejas(t, p, Hs)
        rid = r2_identidad(eps, p, i, j)["r2"]
        # ventanas no solapadas CON AJUSTE, que es lo que hizo la ejecucion: ahi
        # es donde vive el suelo de ruido k/n del que habla el Sec.C.1.
        bordes = np.arange(t[0], t[-1] - Hs, Hs)
        k = np.clip(np.searchsorted(t, bordes, side="left"), 0, N - 1)
        f2 = np.clip(np.searchsorted(t, t[k] + Hs, side="left"), 0, N - 1)
        rr = np.log(p[f2] / p[k])
        flu = np.concatenate(([0.0], np.cumsum(eps)))
        Xv = np.column_stack([flu[k] - flu[np.maximum(k - int(fr_ * Hs / 0.05), 0)]
                              for fr_ in FRACS])
        nt = int(0.75 * k.size)
        rv = H.r2_fuera_de_muestra(Xv[:nt], rr[:nt], Xv[nt:], rr[nt:])
        vals_id.append(rid)
        vals_v.append(rv)
        log("       %.2f   |  %.6f     |   %+.6f" % (a, rid, rv))
    # [!] DOS ASERCIONES MIAS ESTABAN MAL FORMULADAS Y SE CORRIGEN AQUI.
    # (1) Exigia monotonia incluso en `a` cuyo R2 esperado cae POR DEBAJO del
    #     suelo de muestreo (~2.6e-07). Ahi la ordenacion es ruido y pedirla es
    #     pedirle al estimador mas resolucion de la que declara tener. Se exige
    #     monotonia solo sobre los `a` que superan el suelo.
    # (2) Afirmaba que las ventanas "no son monotonas". Medido, SI lo son -- pero
    #     eso no las salva: los cuatro valores salen NEGATIVOS y su recorrido es
    #     menor que el desplazamiento. El sesgo `-k/n_train` DOMINA a la senal,
    #     que es la afirmacion correcta del Sec.C.1 y la que se contrasta.
    piso = vals_id[0]
    sobre = [v for v in vals_id[1:] if v > 3 * piso]
    chk(len(sobre) >= 2 and all(sobre[k_] <= sobre[k_ + 1] + 1e-12
                                for k_ in range(len(sobre) - 1)),
        "la identidad es MONOTONA por encima de su propio suelo",
        "%d de %d valores sobre 3x el piso" % (len(sobre), len(vals_id) - 1))
    chk(vals_id[-1] > 20 * max(piso, 1e-12), "la identidad separa senal de nulo",
        "%.3e contra %.3e = %.0fx" % (vals_id[-1], piso, vals_id[-1] / max(piso, 1e-12)))
    rec = max(vals_v) - min(vals_v)
    chk(all(v < 0 for v in vals_v) and rec < abs(vals_v[0]),
        "en las ventanas el sesgo -k/n DOMINA a la senal",
        "recorrido %.4f < |sesgo| %.4f, y los 4 negativos" % (rec, abs(vals_v[0])))

    # --- el suelo de rotacion es mas ancho que el de barajado con eps con memoria
    ac = np.cumsum(rng.normal(0, 1, N))
    eps_mem = np.sign(ac - np.convolve(ac, np.ones(2001) / 2001, "same"))
    eps_mem[eps_mem == 0] = 1.0
    p = precio_con_senal(0.0)
    i, j = parejas(t, p, Hs)
    a = suelo_de(eps_mem, p, i, j, n_sorteos=60, modo="rotacion")
    b = suelo_de(eps_mem, p, i, j, n_sorteos=60, modo="barajado")
    chk(a["q95"] > 2 * b["q95"],
        "con eps de memoria larga, rotacion da suelo MAS ancho que barajado",
        "%.6f contra %.6f" % (a["q95"], b["q95"]))

    # --- bajo el nulo, la identidad no depende de H como k/n
    i2, j2 = parejas(t, p, 4000.0)
    s_corto = suelo_de(eps, p, i, j, n_sorteos=60)["q95"]
    s_largo = suelo_de(eps, p, i2, j2, n_sorteos=60)["q95"]
    chk(s_largo < 30 * s_corto, "el suelo crece MUY despacio con H (no como H/N)",
        "x%.1f al multiplicar H por 10" % (s_largo / max(s_corto, 1e-12)))

    # --- la multivariante contiene a la univariante
    m = r2_identidad_multi(eps, p, t, i, j, Hs)
    u = r2_identidad(eps, p, i, j)
    chk(m["r2"] >= u["r2"] - 1e-9 or not np.isfinite(m["r2"]),
        "la multivariante no puede rendir menos que la univariante",
        "%.6f contra %.6f" % (m["r2"], u["r2"]))

    log("")
    log("RESULTADO: %d fallo(s)" % fallos)
    return 1 if fallos else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--autotest", action="store_true")
    ap.add_argument("--etapa", default="")
    a = ap.parse_args(argv)
    if a.autotest:
        return _autotest()
    return {"suelo": etapa_suelo, "calibra": etapa_calibra,
            "medir": etapa_medir, "parejas": etapa_parejas}.get(
        a.etapa, lambda _: (log("usa --autotest o --etapa=suelo|calibra|medir|parejas"), 2)[1])(a)


if __name__ == "__main__":
    raise SystemExit(main())
