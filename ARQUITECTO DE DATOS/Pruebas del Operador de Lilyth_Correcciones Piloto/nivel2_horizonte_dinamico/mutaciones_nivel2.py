# -*- coding: utf-8 -*-
"""Nivel 2 Sellado (EXPECTATIVAS_NIVEL2.md §3):
Caza de mutantes reales (MN2.0 + M_N2_1..M_N2_3) sobre horizonte de Hosking,
pre-centrado mu_serie en innovaciones() y cortes duros en Z1b (hist_seg_purgada_).

Soporta --solo-compilar para validar en seco unicidad de patrones (count == 1) y compilacion /W4.
"""
import concurrent.futures as cf
import os
import re
import shutil
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
RAIZ_PILOTO = os.path.normpath(os.path.join(AQUI, ".."))
INC_DIR = os.path.join(RAIZ_PILOTO, "inc")
VCVARS = r"C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat"

MOD = os.path.join("inc", "ark", "lilyth", "modelo_lineal_nulo.hpp")
OBS = os.path.join("inc", "ark", "lilyth", "observador_lilyth.hpp")

MUTANTES = [
    ("MN2.0", "testigo: sin cambio (Nivel 2 sellado)", None, None, None),
    ("M_N2_1", "omitir restado de mu_serie en innovaciones() (modelo_lineal_nulo.hpp)", MOD, "NV2.2",
     [
         ("    const bool centrar_n2 = m.enmienda_nivel2 && std::isfinite(m.mu_serie);\n"
          "    const double mu_pre = centrar_n2 ? m.mu_serie : 0.0;",
          "    const bool centrar_n2 = m.enmienda_nivel2 && std::isfinite(m.mu_serie);\n"
          "    const double mu_pre = 0.0; (void)centrar_n2;")
     ]),
    ("M_N2_2", "sustituir horizonte_inv_arfima_hosking por horizonte_sim_gauss(+d) en ajustar()", MOD, "NV2.2",
     [
         ("m.m_inv_arfima = detalle_lineal::horizonte_inv_arfima_hosking(m.d, eps_q2, cap_inv);",
          "bool top_mut_inv = false; m.m_inv_arfima = detalle_lineal::horizonte_sim_gauss(m.d, eps_q2, cap_inv, top_mut_inv); (void)top_mut_inv;")
     ]),
    ("M_N2_3", "ignorar corte por hist_seg_purgada_ en el escaneo Z1b (observador_lilyth.hpp)", OBS, "NV2.3b",
     [
         ("if (es_purgada && (!tiene_ap || hist_seg_purgada_[idx_out] != seg_apertura)) {",
          "if (es_purgada && !tiene_ap && hist_seg_purgada_[idx_out] != seg_apertura) {")
     ]),
]


def compilar_uno(idx, mid, desc, rel, killer_esp, cambios, scratch):
    d = os.path.join(scratch, "mn2_%d" % idx)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(INC_DIR, os.path.join(d, "inc"))
    os.makedirs(os.path.join(d, "nivel2_horizonte_dinamico"), exist_ok=True)
    shutil.copy(os.path.join(AQUI, "test_nivel2_horizonte.cpp"),
                os.path.join(d, "nivel2_horizonte_dinamico", "test_nivel2_horizonte.cpp"))

    if rel is not None:
        f = os.path.join(d, rel)
        s = open(f, encoding="utf-8").read()
        for viejo, nuevo in cambios:
            n = s.count(viejo)
            if n != 1:
                return mid, desc, killer_esp, None, "MAL DEFINIDO (%d apariciones)" % n
            s = s.replace(viejo, nuevo)
        open(f, "w", encoding="utf-8", newline="\n").write(s)

    src = os.path.join(d, "nivel2_horizonte_dinamico", "test_nivel2_horizonte.cpp")
    exe = os.path.join(d, "test_nivel2_horizonte.exe")
    bat = os.path.join(d, "c.bat")
    with open(bat, "w") as b:
        b.write('@echo off\r\ncall "%s" >nul 2>&1\r\ncd /d "%s"\r\n' % (VCVARS, d))
        b.write('cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%s" "%s" /Fe"%s" /link /STACK:1073741824\r\n'
                % (os.path.join(d, "inc"), src, exe))
    c = subprocess.run(["cmd", "/c", bat], capture_output=True, text=True, encoding="utf-8", errors="replace")
    warns = [l for l in (c.stdout + c.stderr).splitlines() if "warning C" in l]
    if not os.path.exists(exe):
        err = [l for l in (c.stdout + c.stderr).splitlines() if "error" in l][:2]
        return mid, desc, killer_esp, None, "NO COMPILA: %s" % err
    if warns:
        return mid, desc, killer_esp, None, "COMPILA CON WARNINGS (/W4, %d): %s" % (len(warns), warns[0])
    return mid, desc, killer_esp, exe, "COMPILA LIMPIO (/W4, 0 warnings)"


def ejecutar_uno(mid, desc, killer_esp, exe, estado_comp):
    if exe is None:
        return mid, desc, estado_comp, False
    r = subprocess.run([exe], capture_output=True, text=True, encoding="utf-8", errors="replace")
    fallas = [l for l in r.stdout.splitlines() if l.startswith("[FALLA]")]

    if mid == "MN2.0":
        ok_testigo = (r.returncode == 0) and (len(fallas) == 0) and ("FALLAS NIVEL 2: 0" in r.stdout)
        return mid, desc, ("TESTIGO PASA (0 fallas)" if ok_testigo else "TESTIGO FALLA: %s" % fallas), ok_testigo

    cazado_por_control = any(l.startswith("[FALLA] %s:" % killer_esp) for l in fallas)
    cazado = cazado_por_control and (r.returncode != 0)
    detalle = "; ".join(re.sub(r"^\[FALLA\] ", "", l) for l in fallas) or ("rc %d" % r.returncode)
    if not cazado:
        if r.returncode != 0 and not fallas:
            mot = "CRASH/SIN SALIDA (rc=%d)" % r.returncode
        else:
            mot = "VIVO" if not fallas else ("NO CAZADO POR %s (%s)" % (killer_esp, detalle))
        return mid, desc, mot, False
    return mid, desc, "CAZADO (%s): %s" % (killer_esp, detalle), True


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    solo_compilar = "--solo-compilar" in sys.argv
    args_pos = [a for a in sys.argv[1:] if a != "--solo-compilar"]
    scratch = args_pos[0] if args_pos else os.path.join(os.environ.get("TEMP", r"C:\Temp"), "ark_mut_nivel2_sellado")
    os.makedirs(scratch, exist_ok=True)
    print("====================================================================================================")
    print(" NIVEL 2 SELLADO: CAZA DE MUTANTES (MN2.0 + M_N2_1..M_N2_3 | 0 warnings /W4)")
    print("====================================================================================================")

    with cf.ThreadPoolExecutor(4) as ex:
        futs = [ex.submit(compilar_uno, i, *m, scratch) for i, m in enumerate(MUTANTES)]
        compilados = [f.result() for f in futs]

    if solo_compilar:
        todos_ok = True
        for mid, desc, killer_esp, exe, estado in compilados:
            print("%-8s %-72s %s" % (mid, desc, estado))
            if exe is None:
                todos_ok = False
        print("====================================================================================================")
        sys.exit(0 if todos_ok else 1)

    res = [ejecutar_uno(mid, desc, killer_esp, exe, estado) for (mid, desc, killer_esp, exe, estado) in compilados]

    cazados = 0
    for mid, desc, estado, ok in res:
        print("%-8s %-72s %s" % (mid, desc, estado))
        if mid != "MN2.0" and ok:
            cazados += 1
    total_mut = len(MUTANTES) - 1
    testigo_ok = res[0][3]
    print("----------------------------------------------------------------------------------------------------")
    print("RESULTADO NIVEL 2 MUTANTES: %d de %d mutantes cazados | Testigo MN2.0: %s (%d de %d controles)"
          % (cazados, total_mut, "PASA" if testigo_ok else "FALLA",
             cazados + (1 if testigo_ok else 0), len(MUTANTES)))
    print("====================================================================================================")
    sys.exit(0 if (testigo_ok and cazados == total_mut) else 1)


if __name__ == "__main__":
    main()
