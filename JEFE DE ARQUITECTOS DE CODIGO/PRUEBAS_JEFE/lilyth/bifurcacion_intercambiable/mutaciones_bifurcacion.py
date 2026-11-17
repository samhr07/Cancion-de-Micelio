# -*- coding: utf-8 -*-
"""XB.4: mutantes de la bifurcacion del contraste del resumen de Lilyth (Code, 2026-11-15).

Cada mutante es una copia de la base nueva (ark/ entero) en una carpeta corta del scratchpad (las rutas largas pasaban de 260 caracteres)
con un cambio de texto. Se compilan contra ella prueba_intercambiable.cpp y volcado_ab.cpp (/DBASE_NUEVA) y se corren:
  U   la prueba de unidad (XB.1)                    -> cazado si cambia el conjunto de [FALLA] o el codigo de salida contra el testigo
  X0  volcado 'nueva' del mundo IID                 -> cazado si no es identico al bloque IID del volcado de la base previa (XB.0 ii)
  X0b volcado 'inter' del mundo IID                 -> cazado si sus lineas no H no son identicas a las de su propio 'nueva' (XB.0b)
  X0c contadores de XB.0c del volcado 'inter'       -> cazado si alguno no es 0
Aparte (no es un control): si las lineas HI (la rama nueva) cambian contra el testigo. Un mutante que no compila NO esta cazado.
Uso: python -I mutaciones_bifurcacion.py SCRATCH [paralelo] [hilos] [ids]
"""
import concurrent.futures as cf
import os
import re
import shutil
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(AQUI, "..", "..", "..", "..", "ark_cxx", "base_auditada", "include", "ark"))
RED = os.path.normpath(os.path.join(AQUI, "..", "..", "red"))
NUL = "ark/lilyth/nulo_lilyth.hpp"
OBS = "ark/lilyth/observador_lilyth.hpp"

MUTANTES = [
    ("M0", "testigo: sin cambio", []),
    ("M1", "real fuera del estrato en el Monte Carlo", [(NUL, "for (const auto& uk : u) dist.emplace_back(0, uk.size() - 1);",
                                                         "for (const auto& uk : u) dist.emplace_back(1, uk.size() - 1);", False)]),
    ("M2", "Monte Carlo del pool global (estrato al azar)", [(NUL, "for (std::size_t k = 0; k < K; ++k) sb += u[k][dist[k](rng)] - 0.5;",
                                                               "for (std::size_t k = 0; k < K; ++k) { const std::size_t kk = rng() % K; sb += u[kk][std::uniform_int_distribution<std::size_t>(0, u[kk].size() - 1)(rng)] - 0.5; }", False)]),
    ("M3", "p sin el +1", [(NUL, "    res.p_bilateral = static_cast<double>(1 + extremos) / static_cast<double>(1 + res.b_mc);\n\n    // p exacto",
                            "    res.p_bilateral = static_cast<double>(extremos) / static_cast<double>(res.b_mc);\n\n    // p exacto", False)]),
    ("M4", "p_min sin la cola inferior", [(NUL, "        if (c >= s_max - 1e-15) p_ext += pc;\n", "", False)]),
    ("M5", "resoluble con alfa en vez de alfa/4", [(OBS, "cfg_.alfa / 4.0)", "cfg_.alfa)", True)]),
    ("M6", "omision invertida (INTERCAMBIABLE por omision)", [(OBS, "ContrasteResumenLilyth contraste_resumen = ContrasteResumenLilyth::POOL_LOO;",
                                                               "ContrasteResumenLilyth contraste_resumen = ContrasteResumenLilyth::INTERCAMBIABLE;", False)]),
    ("M7", "sobrepaso por evento sin bifurcar", [(OBS, "        if (cfg_.contraste_resumen == ContrasteResumenLilyth::INTERCAMBIABLE) {\n            arm.sobrepaso.sobrepaso_por_evento_iid_pool_loo",
                                                  "        if (false) {\n            arm.sobrepaso.sobrepaso_por_evento_iid_pool_loo", False)]),
    ("M8", "contraste_iid_pool_loo sin llenar", [(OBS, "            fn.contraste_iid_pool_loo = fn.contraste_iid;\n", "", False)]),
    ("M9", "u del real contra el estrato que lo incluye", [(NUL, "            for (std::size_t q = 0; q < e.size(); ++q) if (q != m) resto.push_back(e[q]);\n            uk[m]",
                                                            "            for (std::size_t q = 0; q < e.size(); ++q) resto.push_back(e[q]);\n            uk[m]", False)]),
    ("M10", "Holm de 4 con el p viejo en INTERCAMBIABLE", [(OBS, "const double p_dt_iid = fn.contraste_iid.p_bilateral;",
                                                             "const double p_dt_iid = cfg_.contraste_resumen == ContrasteResumenLilyth::INTERCAMBIABLE ? fn.contraste_iid_pool_loo.p_bilateral : fn.contraste_iid.p_bilateral;", False)]),
]

SCRATCH = None
HILOS = "10"


def bloque_iid(txt):
    out, dentro = [], False
    for l in txt.splitlines():
        if l.startswith("MUNDO "):
            dentro = l.strip() == "MUNDO IID"
            continue
        if dentro:
            out.append(l)
    return out


def compilar(d, fuente, nombre, extra):
    bat = os.path.join(d, nombre + ".bat")
    with open(bat, "w") as f:
        f.write('@echo off\r\ncall "C:\\Program Files\\Microsoft Visual Studio\\18\\Community\\VC\\Auxiliary\\Build\\vcvars64.bat" >nul 2>&1\r\n')
        f.write('cl /nologo %s /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%s" /I"%s" /I"%s" "%s" /Fo"%s\\\\" /Fe"%s" /link /STACK:1073741824\r\n'
                % (extra, d, os.path.join(RED, "inc"), RED, fuente, os.path.join(d, "o_" + nombre), os.path.join(d, nombre + ".exe")))
    os.makedirs(os.path.join(d, "o_" + nombre), exist_ok=True)
    subprocess.run(["cmd", "/c", bat], capture_output=True, text=True, encoding="utf-8", errors="replace")
    exe = os.path.join(d, nombre + ".exe")
    return exe if os.path.exists(exe) else None


def correr(m):
    mid, desc, cambios = m
    d = os.path.join(SCRATCH, "mb", mid)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(BASE, os.path.join(d, "ark"))
    for rel, viejo, nuevo, todos in cambios:
        p = os.path.join(d, rel)
        s = open(p, encoding="utf-8").read()
        n = s.count(viejo)
        if n == 0 or (n != 1 and not todos):
            return mid, desc, "MAL DEFINIDO (%d apariciones)" % n, None
        s = s.replace(viejo, nuevo)
        open(p, "w", encoding="utf-8", newline="").write(s)
    u = compilar(d, os.path.join(AQUI, "prueba_intercambiable.cpp"), "u", "")
    v = compilar(d, os.path.join(AQUI, "volcado_ab.cpp"), "v", "/DBASE_NUEVA")
    if u is None or v is None:
        return mid, desc, "NO COMPILA (no esta cazado)", None
    ru = subprocess.run([u], capture_output=True, text=True, encoding="utf-8", errors="replace")
    fallas_u = tuple(sorted(l.split()[1] for l in ru.stdout.splitlines() if l.startswith("[FALLA]")))
    rn = subprocess.run([v, HILOS, "nueva", "IID"], capture_output=True, text=True, encoding="utf-8", errors="replace")
    ri = subprocess.run([v, HILOS, "inter", "IID"], capture_output=True, text=True, encoding="utf-8", errors="replace")
    n_iid, i_iid = bloque_iid(rn.stdout), bloque_iid(ri.stdout)
    # (corregido tras la primera pasada: el findall de todos los numeros recogia el "4" de "Holm de 4")
    m = re.search(r"u_pit (\d+); estrato distinto de n_nulos_usados (\d+); series con p del Holm de 4 distinto del contraste nuevo (\d+)", ri.stderr)
    xb0c = tuple(int(x) for x in m.groups()) if m else None
    return mid, desc, None, {"u": (ru.returncode, fallas_u), "nueva": n_iid, "inter": i_iid, "xb0c": xb0c,
                             "hi": [l for l in i_iid if l.startswith("HI")], "rc": (rn.returncode, ri.returncode)}


def main():
    global SCRATCH, HILOS
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    SCRATCH = sys.argv[1]
    par = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    HILOS = sys.argv[3] if len(sys.argv) > 3 else "10"
    solo = sys.argv[4].split(",") if len(sys.argv) > 4 else None
    previa = bloque_iid(open(os.path.join(AQUI, "volcado_previa.txt"), encoding="utf-8").read())
    t = correr(MUTANTES[0])
    T = t[3]
    print("M0   testigo: unidad rc %s fallas %s; nueva == previa %s; XB.0c %s" % (T["u"][0], T["u"][1], T["nueva"] == previa, T["xb0c"]))
    sys.stdout.flush()
    with cf.ThreadPoolExecutor(par) as ex:
        res = list(ex.map(correr, [m for m in MUTANTES[1:] if solo is None or m[0] in solo]))
    cazados = 0
    for mid, desc, est, R in res:
        if est is not None:
            print("%-4s %-50s %s" % (mid, desc, est))
            continue
        por = []
        if R["u"] != T["u"]: por.append("U %s %s" % R["u"])
        if R["nueva"] != previa: por.append("X0")
        no_h = lambda b: [l for l in b if not l.startswith("H")]
        if no_h(R["inter"]) != no_h(R["nueva"]): por.append("X0b")
        # XB.0c contra el testigo (el testigo tiene p distintos en las series con T_max = 2, sin Holm de 4 en ninguna rama)
        if R["xb0c"] is None or R["xb0c"] != T["xb0c"]: por.append("X0c %s" % (R["xb0c"],))
        if R["rc"] != (0, 0): por.append("rc %s" % (R["rc"],))
        cambia = R["hi"] != T["hi"]
        if por:
            cazados += 1
            print("%-4s %-50s CAZADO (%s)" % (mid, desc, "; ".join(por)))
        else:
            print("%-4s %-50s NO CAZADO%s" % (mid, desc, " (cambia la salida de la rama nueva)" if cambia else " (salida identica: equivalente en estos mundos)"))
        sys.stdout.flush()
    print("%d/%d cazados" % (cazados, len(res)))


if __name__ == "__main__":
    main()
