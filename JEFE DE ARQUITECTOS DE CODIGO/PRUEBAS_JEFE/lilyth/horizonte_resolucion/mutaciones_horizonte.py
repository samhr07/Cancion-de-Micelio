# -*- coding: utf-8 -*-
"""ET5 de la enmienda de E1.9 (horizonte; Code, 2026-11-16): cinco mutantes sellados sobre la base nueva.

Cada mutante: copia de ark_cxx/base_auditada/include en <scratch>/mN, un reemplazo que tiene que aparecer exactamente una vez, compila
control_horizonte.cpp contra la copia y lo corre (pool 400, fuente 40). Cazado = el control sale con codigo != 0 (alguna condicion falla)
o no compila o revienta. Uso: python -I mutaciones_pool_fuente.py <scratch> [paralelo] [hilos]
"""
import concurrent.futures as cf
import os
import re
import shutil
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(AQUI, "..", "..", "..", "..", "ark_cxx", "base_auditada", "include"))
OBS = os.path.join("ark", "lilyth", "observador_lilyth.hpp")
NUL = os.path.join("ark", "lilyth", "nulo_lilyth.hpp")
MUT = [
    ("H1 nivel alfa en vez de alfa/8", NUL, "return h == HorizonteLilyth::CUANTIL_RESOLUCION ? alfa / (2.0 * 4.0) : alfa;",
     "return h == HorizonteLilyth::CUANTIL_RESOLUCION ? alfa : alfa;"),
    ("H2 nivel alfa/4", NUL, "return h == HorizonteLilyth::CUANTIL_RESOLUCION ? alfa / (2.0 * 4.0) : alfa;",
     "return h == HorizonteLilyth::CUANTIL_RESOLUCION ? alfa / 4.0 : alfa;"),
    ("H3 guarda de censura en el alfa del proyecto", NUL, "    if (out_frac_cens > alfa) {", "    if (out_frac_cens > 0.05) {"),
    ("H4 cuantil alfa en vez de 1 - alfa", NUL, "const double q = ark::metrologia::est::cuantil(taus, 1.0 - alfa);",
     "const double q = ark::metrologia::est::cuantil(taus, alfa);"),
    ("H5 CUANTIL_RESOLUCION por omision", OBS, "HorizonteLilyth horizonte = HorizonteLilyth::CUANTIL_ALFA;",
     "HorizonteLilyth horizonte = HorizonteLilyth::CUANTIL_RESOLUCION;"),
]
VCVARS = r"C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"


def uno(i, nombre, rel, viejo, nuevo, scratch, hilos):
    d = os.path.join(scratch, "m%d" % (i + 1))
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(BASE, os.path.join(d, "include"))
    f = os.path.join(d, "include", rel)
    s = open(f, encoding="utf-8").read()
    n = s.count(viejo)
    if n != 1:
        return nombre, "REEMPLAZO INVALIDO (%d apariciones)" % n, False
    open(f, "w", encoding="utf-8", newline="").write(s.replace(viejo, nuevo))
    src = os.path.join(d, "control_m.cpp")
    shutil.copy(os.path.join(AQUI, "control_horizonte.cpp"), src)
    exe = os.path.join(d, "control_m.exe")
    bat = os.path.join(d, "c.bat")
    with open(bat, "w") as b:
        b.write('@echo off\r\ncall "%s" >nul 2>&1\r\ncd /d "%s"\r\n' % (VCVARS, d))
        b.write('cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%s" "%s" /Fe"%s" /link /STACK:1073741824\r\n'
                % (os.path.join(d, "include"), src, exe))
    c = subprocess.run(["cmd", "/c", bat], capture_output=True, text=True, encoding="utf-8", errors="replace")
    if not os.path.exists(exe):
        err = [l for l in (c.stdout + c.stderr).splitlines() if "error" in l][:2]
        return nombre, "NO COMPILA: %s" % err, False
    r = subprocess.run([exe, str(hilos), "todo", "400", "40", "resolucion"], capture_output=True, text=True, encoding="utf-8", errors="replace")
    with open(os.path.join(AQUI, "build", "mutante_%d.log" % (i + 1)), "w", encoding="utf-8") as g:
        g.write(r.stdout + r.stderr)
    fallas = [l for l in r.stdout.splitlines() if l.startswith("[FALLA]") and not l.startswith("[FALLA] condicion EC6a (")]
    cazado = bool(fallas) or (r.returncode not in (0, 1)) or "FALLAS DE CONDICIONES" not in r.stdout
    detalle = "; ".join(re.sub(r"^\[FALLA\] condicion ", "", l)[:110] for l in fallas[:3]) or ("rc %d" % r.returncode)
    return nombre, ("CAZADO: " if cazado else "VIVO: ") + detalle, cazado


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    scratch = sys.argv[1]
    par = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    hilos = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    os.makedirs(os.path.join(AQUI, "build"), exist_ok=True)
    # EC6a sellado falla en la base sin mutar (condicion mal especificada, declarada): no cuenta para cazar.
    with cf.ThreadPoolExecutor(par) as ex:
        futs = [ex.submit(uno, i, *m, scratch, hilos) for i, m in enumerate(MUT)]
        res = [f.result() for f in futs]
    cazados = 0
    for nombre, estado, ok in res:
        print("%-52s %s" % (nombre, estado))
        cazados += ok
    print("ET5: %d de %d cazados (>= 4 esperados)" % (cazados, len(MUT)))


if __name__ == "__main__":
    main()
