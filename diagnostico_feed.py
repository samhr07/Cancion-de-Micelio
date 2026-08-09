# -*- coding: utf-8 -*-
"""
diagnostico_feed.py -- catalogo de perdidas de datos de una captura.

    python diagnostico_feed.py --dir=telemetria/captura_v32
    python diagnostico_feed.py --dir=telemetria/captura_v32 --log=captura_v32.log

No modifica nada ni toca la captura en curso: solo lee lo que hay en disco. Se
puede correr tantas veces como se quiera mientras la captura crece.

POR QUE HACE FALTA
------------------
La compuerta de la v3.2 es sobre el tramo continuo MAS LARGO, asi que una
desconexion no cuesta el tiempo que dura: cuesta **todo lo acumulado antes**. Con
741 000 ticks de requisito y cortes cada pocas horas, la pregunta operativa no es
"cuantos datos llevo" sino "cual es la tasa de cortes y me deja llegar".

TRES CLASES DE PERDIDA, Y SOLO UNA SE VE EN EL LOG
---------------------------------------------------
1. **Hueco temporal con desconexion.** El watchdog lo registra. Es el visible.
2. **Hueco temporal SIN error registrado.** El proceso estuvo congelado -- la
   maquina suspendida, el planificador expulsandolo, o el bucle de eventos
   detenido. No hay error porque nunca se ejecuto el codigo que lo detectaria.
   Se reconoce porque el hueco NO va acompanado de lineas `[!]` proporcionales.
3. **Perdida SILENCIOSA: hueco de `tr_id` sin hueco temporal.** El feed sigue
   conectado y aun asi faltan transacciones. Ninguna logica de reconexion la ve,
   y es la mas peligrosa para el estimador porque no deja rastro en el eje de
   tiempos: la convolucion cruza el hueco como si no existiera.

⚠ Los ids tambien saltan por el filtro `tick_valido` (~0.2 % de trades con p=0),
asi que un salto de 1 no es prueba de perdida. Se reportan por tamano.
"""

from __future__ import annotations

import argparse
import os
import re

import numpy as np

import captura_larga as cl


HUECO_MIN = 60.0     # [s] por debajo de esto es cadencia normal de mercado lento


def eventos_de_hueco(t: np.ndarray, ident: np.ndarray,
                     umbral: float = HUECO_MIN) -> list:
    """Huecos temporales, con cuantas transacciones se perdieron en cada uno."""
    dt = np.diff(t)
    idx = np.flatnonzero(dt > umbral)
    ev = []
    for i in idx:
        salto_id = int(ident[i + 1] - ident[i] - 1) if ident is not None else -1
        ev.append({
            "i": int(i),
            "t_ini": float(t[i]), "t_fin": float(t[i + 1]),
            "duracion_s": float(dt[i]),
            "trades_perdidos": salto_id,
            # Tasa implicita durante el hueco, si los ids son de fiar.
            "tasa_implicita": (salto_id / dt[i]) if salto_id > 0 else float("nan"),
        })
    return ev


def perdidas_silenciosas(t: np.ndarray, ident: np.ndarray,
                         umbral_t: float = HUECO_MIN) -> dict:
    """Saltos de `tr_id` SIN hueco temporal: datos perdidos estando conectado.

    Es la clase que ninguna reconexion detecta. Se separa por tamano porque el
    filtro de ceros del feed produce saltos de 1-2 que no son perdida.
    """
    dt = np.diff(t)
    salto = np.diff(ident) - 1
    sin_hueco = (dt <= umbral_t) & (salto > 0)
    tam = salto[sin_hueco]
    return {
        "n_eventos": int(np.count_nonzero(sin_hueco)),
        "trades_perdidos": int(tam.sum()) if tam.size else 0,
        "de_1": int(np.count_nonzero(tam == 1)),
        "de_2": int(np.count_nonzero(tam == 2)),
        "de_3_a_10": int(np.count_nonzero((tam >= 3) & (tam <= 10))),
        "mayores_de_10": int(np.count_nonzero(tam > 10)),
        "mayor": int(tam.max()) if tam.size else 0,
        "fraccion": float(tam.sum() / max(1, ident[-1] - ident[0])),
    }


def parsear_log(ruta: str) -> dict:
    """Cuenta y clasifica las lineas de error del log de captura.

    El log no lleva marca de tiempo en las lineas de error, asi que solo se usa
    para CONTAR y CLASIFICAR; la cronologia sale de los datos, que si estan
    fechados por el exchange.
    """
    if not ruta or not os.path.exists(ruta):
        return {}
    tipos, bloques = {}, []
    with open(ruta, "r", encoding="utf-8", errors="replace") as fh:
        for linea in fh:
            m = re.search(r"\[!\]\s+(\w+):", linea)
            if m:
                tipos[m.group(1)] = tipos.get(m.group(1), 0) + 1
                continue
            m = re.search(r"\[\s*([\d.]+) h\]", linea)
            if m:
                bloques.append(float(m.group(1)))
    return {"tipos": tipos, "n_errores": sum(tipos.values()),
            "horas_bloque": bloques}


def probabilidad_de_llegar(tramos_ticks, objetivo: int, nu: float) -> dict:
    """Con la tasa de cortes observada, ¿se llega al objetivo en UN tramo?

    Modelo de renovacion simple: los cortes llegan a tasa constante `lambda`
    estimada como (n_cortes / tiempo_total). Un tramo continuo sobrevive hasta
    el objetivo con probabilidad exp(-lambda * T_objetivo).

    ⚠ Es un modelo grosero -- supone cortes independientes y tasa constante, y
    con dos o tres eventos la estimacion de lambda es pesima. Sirve para ordenar
    magnitudes y decidir si hace falta cambiar algo, no para prometer nada.
    """
    n_cortes = max(0, len(tramos_ticks) - 1)
    total_ticks = int(sum(tramos_ticks))
    if n_cortes == 0 or nu <= 0:
        return {"n_cortes": 0, "p": float("nan"),
                "ticks_medios_por_tramo": total_ticks,
                "nota": "sin cortes observados: no hay tasa que estimar"}
    ticks_medios = total_ticks / (n_cortes + 1)
    p = float(np.exp(-objetivo / ticks_medios))
    return {"n_cortes": n_cortes, "ticks_medios_por_tramo": ticks_medios,
            "p": p, "horas_por_corte": (ticks_medios / nu) / 3600.0,
            "nota": "modelo de renovacion con tasa constante; grosero"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="telemetria/captura_v32")
    ap.add_argument("--log", default=None)
    ap.add_argument("--umbral", type=float, default=HUECO_MIN)
    a = ap.parse_args()

    d = cl.cargar_larga(a.dir)
    t = d["tr_t"]
    ident = d.get("tr_id")
    n = len(t)
    dur = float(t[-1] - t[0])
    nu = n / dur if dur > 0 else float("nan")

    print("=" * 74)
    print("CATALOGO DE PERDIDAS DE DATOS")
    print("=" * 74)
    print("directorio   : %s" % a.dir)
    print("transacciones: %d en %.2f h (nu media %.2f tx/s)" % (n, dur / 3600.0, nu))
    print("")

    # --- 1. tramos continuos ---
    tramos = cl.tramos_continuos(t, cl.HUECO_MAX)
    largos = [b - a_ for a_, b in tramos]
    print("-" * 74)
    print("TRAMOS CONTINUOS (corte a %.0f s)" % cl.HUECO_MAX)
    print("-" * 74)
    print("%-4s %10s %10s %12s %10s" % ("#", "ticks", "horas", "nu [tx/s]", "% del total"))
    for k, (i0, i1) in enumerate(tramos):
        dh = float(t[i1 - 1] - t[i0])
        print("%-4d %10d %10.2f %12.2f %9.1f %%"
              % (k, i1 - i0, dh / 3600.0,
                 (i1 - i0) / dh if dh > 0 else float("nan"),
                 100.0 * (i1 - i0) / n))
    k_mejor = int(np.argmax(largos))
    print("tramo utilizable: #%d con %d ticks  (la compuerta pide %d)"
          % (k_mejor, largos[k_mejor], cl.TICKS_COMPUERTA))

    # --- 2. huecos temporales ---
    ev = eventos_de_hueco(t, ident, a.umbral)
    print("")
    print("-" * 74)
    print("HUECOS TEMPORALES > %.0f s" % a.umbral)
    print("-" * 74)
    if not ev:
        print("ninguno")
    else:
        print("%-6s %10s %12s %14s %12s"
              % ("#", "dur [s]", "dur [min]", "trades perdidos", "tasa impl."))
        for k, e in enumerate(ev):
            print("%-6d %10.1f %12.2f %14d %12.2f"
                  % (k, e["duracion_s"], e["duracion_s"] / 60.0,
                     e["trades_perdidos"], e["tasa_implicita"]))
        total_perdido = sum(e["duracion_s"] for e in ev)
        print("tiempo total en hueco: %.2f h (%.1f %% de la captura)"
              % (total_perdido / 3600.0, 100.0 * total_perdido / dur))

    # --- 3. perdidas silenciosas ---
    print("")
    print("-" * 74)
    print("PERDIDA SILENCIOSA (salto de id SIN hueco temporal)")
    print("-" * 74)
    if ident is None:
        print("sin tr_id: no se puede evaluar")
    else:
        s = perdidas_silenciosas(t, ident, a.umbral)
        print("eventos            : %d" % s["n_eventos"])
        print("trades perdidos    : %d (%.3f %% del rango de ids)"
              % (s["trades_perdidos"], 100.0 * s["fraccion"]))
        print("por tamano         : de 1 = %d | de 2 = %d | 3-10 = %d | >10 = %d | mayor = %d"
              % (s["de_1"], s["de_2"], s["de_3_a_10"], s["mayores_de_10"], s["mayor"]))
        print("AVISO: los saltos de 1-2 son compatibles con el filtro de ceros del")
        print("       feed (~0.2 %% de trades con p=0). Los mayores de 10 no lo son.")

    # --- 4. log ---
    ruta_log = a.log
    if ruta_log is None:
        cand = os.path.basename(os.path.normpath(a.dir)) + ".log"
        ruta_log = cand if os.path.exists(cand) else None
    pl = parsear_log(ruta_log)
    print("")
    print("-" * 74)
    print("LOG DE LA CAPTURA")
    print("-" * 74)
    if not pl:
        print("sin log (%s)" % ruta_log)
    else:
        print("fichero: %s | %d errores registrados" % (ruta_log, pl["n_errores"]))
        for k, v in sorted(pl["tipos"].items(), key=lambda kv: -kv[1]):
            print("   %-32s %d" % (k, v))
        # Contraste clave: un hueco largo con pocos errores NO es un fallo de red.
        if ev and pl["n_errores"] >= 0:
            mayor = max(e["duracion_s"] for e in ev)
            # El watchdog reintenta con backoff acotado a 60 s, asi que un corte
            # de red de D segundos deberia dejar del orden de D/60 lineas.
            esperadas = mayor / 60.0
            print("")
            print("   hueco mayor: %.0f s -> un corte de RED dejaria del orden de %.0f"
                  % (mayor, esperadas))
            print("   lineas de error (el backoff se satura en 60 s). Registradas: %d."
                  % pl["n_errores"])
            if pl["n_errores"] < 0.25 * esperadas:
                print("   VEREDICTO: demasiado pocas. El proceso estuvo CONGELADO, no")
                print("              desconectado -- suspension de la maquina o el bucle")
                print("              de eventos detenido. El watchdog no fallo: nunca")
                print("              llego a ejecutarse.")
            else:
                print("   VEREDICTO: compatible con un corte de red sostenido.")

    # --- 5. proyeccion ---
    print("")
    print("-" * 74)
    print("PROYECCION")
    print("-" * 74)
    pr = probabilidad_de_llegar(largos, cl.TICKS_COMPUERTA, nu)
    if pr["n_cortes"] == 0:
        print("sin cortes observados todavia: no hay tasa que estimar")
    else:
        print("cortes observados      : %d en %.2f h" % (pr["n_cortes"], dur / 3600.0))
        print("ticks medios por tramo : %.0f  (%.2f h a la nu media)"
              % (pr["ticks_medios_por_tramo"], pr["horas_por_corte"]))
        if pr["n_cortes"] < 3:
            # ⚠ Con uno o dos eventos la tasa no se puede estimar, y publicar
            # "P = 0.0000" o "harian falta 115 millones de intentos" es dar
            # apariencia de medicion a un dato que no la tiene.
            print("P(un tramo alcance %d): NO ESTIMABLE con %d corte(s)"
                  % (cl.TICKS_COMPUERTA, pr["n_cortes"]))
            print("       hacen falta al menos 3 eventos para que la tasa signifique")
            print("       algo. Lo unico defendible hoy: el tramo medio observado es")
            print("       %.0f ticks contra los %d que pide la compuerta, o sea un")
            print("       factor %.1f." % (pr["ticks_medios_por_tramo"],
                                           cl.TICKS_COMPUERTA,
                                           cl.TICKS_COMPUERTA / pr["ticks_medios_por_tramo"]))
        else:
            print("P(un tramo alcance %d): %.4f" % (cl.TICKS_COMPUERTA, pr["p"]))
            print("AVISO: %s" % pr["nota"])
            if pr["p"] < 0.5:
                n_int = int(np.ceil(1.0 / max(pr["p"], 1e-9)))
                print("       con esta tasa harian falta del orden de %d intentos" % n_int)

    falta = cl.TICKS_COMPUERTA - largos[k_mejor]
    if falta > 0:
        nu_m = largos[k_mejor] / max(1.0, float(t[tramos[k_mejor][1] - 1]
                                                - t[tramos[k_mejor][0]]))
        print("")
        print("faltan %d ticks en el tramo actual = %.1f h a su nu (%.2f tx/s)"
              % (falta, falta / nu_m / 3600.0, nu_m))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
