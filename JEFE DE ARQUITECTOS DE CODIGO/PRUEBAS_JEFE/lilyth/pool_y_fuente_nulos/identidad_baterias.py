# -*- coding: utf-8 -*-
"""EC5 (baterias) de la ultima fase de Lilyth, pool derivado y fuente (Code, 2026-11-16); copia de XB.0 (i) de la bifurcacion.

Compila las baterias de Lilyth de las Etapas 2-6 (fuentes de PRUEBAS_JEFE/lilyth/revision_etapas2_6/, las de la base; las de etapa_lh/ son de la rama L+H, que no esta en la base) contra la base PREVIA (base_previa/) y contra
la base NUEVA (ark_cxx/base_auditada/include), las corre y compara los logs linea a linea quitando las lineas de tiempo. Con POOL_LOO por
omision tienen que ser identicos. Uso: python -I identidad_baterias.py
"""
import concurrent.futures as cf
import os
import re
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
FUENTES = os.path.normpath(os.path.join(AQUI, "..", "revision_etapas2_6"))
PREVIA = os.path.join(AQUI, "base_previa", "include")
NUEVA = os.path.normpath(os.path.join(AQUI, "..", "..", "..", "..", "ark_cxx", "base_auditada", "include"))
BATERIAS = ["etapa2_lilyth", "etapa3_0_lilyth", "etapa3_lilyth", "etapa4_lilyth", "etapa5_lilyth", "etapa6_lilyth"]
TIEMPO = re.compile(r"(\d+(\.\d+)?\s*(s|ms|seg|segundos)\b)|(\bt(iempo)?\s*[=:])", re.I)


def compilar_y_correr(nombre, inc, etiqueta):
    d = os.path.join(AQUI, "build", "identidad", etiqueta)
    os.makedirs(os.path.join(d, nombre), exist_ok=True)
    exe = os.path.join(d, nombre + ".exe")
    bat = os.path.join(d, nombre + ".bat")
    with open(bat, "w") as f:
        f.write('@echo off\r\ncall "C:\\Program Files\\Microsoft Visual Studio\\18\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat" >nul 2>&1\r\n')
        f.write('cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%s" "%s" /Fo"%s\\\\" /Fe"%s" /link /STACK:1073741824\r\n'
                % (inc, os.path.join(FUENTES, nombre + ".cpp"), os.path.join(d, nombre), exe))
    c = subprocess.run(["cmd", "/c", bat], capture_output=True, text=True, encoding="utf-8", errors="replace")
    avisos = [l for l in (c.stdout + c.stderr).splitlines() if "warning" in l or "error" in l]
    if not os.path.exists(exe):
        return None, avisos
    r = subprocess.run([exe], cwd=FUENTES, capture_output=True, text=True, encoding="utf-8", errors="replace")
    salida = r.stdout + r.stderr
    with open(os.path.join(d, nombre + ".log"), "w", encoding="utf-8") as f:
        f.write(salida)
    return (r.returncode, salida), avisos


def una(nombre):
    a, av_a = compilar_y_correr(nombre, PREVIA, "previa")
    b, av_b = compilar_y_correr(nombre, NUEVA, "nueva")
    if a is None or b is None:
        return nombre, "NO COMPILA", av_a + av_b
    la = [l for l in a[1].splitlines() if not TIEMPO.search(l)]
    lb = [l for l in b[1].splitlines() if not TIEMPO.search(l)]
    quitadas = len(a[1].splitlines()) - len(la)
    if a[0] == b[0] and la == lb:
        return nombre, "IDENTICA (rc %s, %d lineas comparadas, %d de tiempo quitadas)" % (a[0], len(la), quitadas), av_a + av_b
    dif = [(i, x, y) for i, (x, y) in enumerate(zip(la, lb)) if x != y][:5]
    return nombre, "DISTINTA (rc %s/%s, lineas %d/%d) primeras: %r" % (a[0], b[0], len(la), len(lb), dif), av_a + av_b


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    malas = 0
    with cf.ThreadPoolExecutor(len(BATERIAS)) as ex:
        for nombre, estado, avisos in ex.map(una, BATERIAS):
            print("%-16s %s" % (nombre, estado))
            for v in avisos[:6]:
                print("    " + v)
            malas += not estado.startswith("IDENTICA")
    print("EC5 baterias: %s" % ("CUMPLE" if malas == 0 else "FALLA (%d baterias)" % malas))
    return 1 if malas else 0


if __name__ == "__main__":
    sys.exit(main())
