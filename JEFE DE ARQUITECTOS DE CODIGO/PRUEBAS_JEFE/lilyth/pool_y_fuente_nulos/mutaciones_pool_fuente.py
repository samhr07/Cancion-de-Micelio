# -*- coding: utf-8 -*-
"""EC7 de la ultima fase de Lilyth (pool derivado y fuente; Code, 2026-11-16): diez mutantes sellados sobre la base nueva.

Cada mutante: copia de ark_cxx/base_auditada/include en <scratch>/mN, un reemplazo que tiene que aparecer exactamente una vez, compila
control_fuente_pool.cpp contra la copia y lo corre (pool 400, fuente 40). Cazado = el control sale con codigo != 0 (alguna condicion falla)
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
    ("M1 mascara ignorada", OBS, "const bool en_fuente = fuente_nulo.empty() || fuente_nulo[i] != 0;", "const bool en_fuente = true;"),
    ("M2 mascara invertida", OBS, "const bool en_fuente = fuente_nulo.empty() || fuente_nulo[i] != 0;", "const bool en_fuente = fuente_nulo.empty() || fuente_nulo[i] == 0;"),
    ("M3 mascara que no corta el segmento", OBS, "if (v && std::isfinite(y) && en_fuente) {\n                    cur_c.push_back(y);\n                    ++nf_c;",
     "if (v && std::isfinite(y)) {\n                    if (en_fuente) cur_c.push_back(y);\n                    nf_c += en_fuente;"),
    ("M4 n* sin el - 1", NUL, "return r > 1.0 ? static_cast<std::size_t>(r) - 1 : 0;", "return r > 1.0 ? static_cast<std::size_t>(r) : 0;"),
    ("M5 n* con 1/alfa en vez de 2*4/alfa", NUL, "const double base = 2.0 * 4.0 / alfa;", "const double base = 1.0 / alfa;"),
    ("M6 via (a) con m_min en vez del umbral", OBS, "if (!cfg_.forzar_via_b && esp.size() >= umbral) {", "if (!cfg_.forzar_via_b && esp.size() >= fn.m_min) {"),
    ("M7 via (b) con m_min en vez del umbral", OBS, "r.signo, r.w_ancla, r.guarda_ancla, umbral, cfg_.t_max_efectivo,",
     "r.signo, r.w_ancla, r.guarda_ancla, fn.m_min, cfg_.t_max_efectivo,"),
    ("M8 K sin los censurados (proxy de 'K = evaluados')", OBS,
     "if (ra.desenlace == DesenlaceRanura::RESTITUIDA || ra.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE) ++fn.k_candidatos_pool;",
     "if (ra.desenlace == DesenlaceRanura::RESTITUIDA) ++fn.k_candidatos_pool;"),
    ("M9 DERIVADO por omision", OBS, "PoolNuloLilyth pool_nulo = PoolNuloLilyth::M_MIN;", "PoolNuloLilyth pool_nulo = PoolNuloLilyth::DERIVADO;"),
    ("M10 piso m_min retirado", OBS, "            if (taus.size() < fn.m_min) {", "            if (taus.empty()) {"),
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
    shutil.copy(os.path.join(AQUI, "control_fuente_pool.cpp"), src)
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
    r = subprocess.run([exe, str(hilos), "todo", "400", "40"], capture_output=True, text=True, encoding="utf-8", errors="replace")
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
    print("EC7: %d de %d cazados (>= 8 esperados)" % (cazados, len(MUT)))


if __name__ == "__main__":
    main()
