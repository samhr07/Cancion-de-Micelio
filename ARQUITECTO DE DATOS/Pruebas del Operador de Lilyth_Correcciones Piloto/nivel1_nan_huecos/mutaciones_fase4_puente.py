# -*- coding: utf-8 -*-
"""Fase 4 Sellada (EXPECTATIVAS_FASE4.md): Caza de mutantes reales sobre puente de Ito,
etiquetas FLAG_PUENTE_ITO, imputacion multiple de Rubin y deriva en FLUJO.

R = 200 fijo desde el preregistro sellado (sin reescalado post-hoc).
Soporta --solo-compilar para validar en seco unicidad de patrones y compilacion /W4.
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

OBS = os.path.join("inc", "ark", "lilyth", "observador_lilyth.hpp")
ITO = os.path.join("inc", "ark", "lilyth", "puente_historico_ito.hpp")

MUTANTES = [
    ("MF4.0", "testigo: sin cambio (Fase 4 sellada)", None, None),
    ("MN1.7", "ignorar FLAG_PUENTE_ITO en Z2 y Z3 dentro de observador_lilyth.hpp", OBS,
     [
         ("const bool veto_z3_ito = ((flags_hueco & ark::dataset::FLAG_PUENTE_ITO) != 0u);",
          "const bool veto_z3_ito = false;"),
         ("const bool veto_ito_z2 = tiene_idx && ((hist_flags_hueco_[idx] & (ark::dataset::FLAG_PUENTE_ITO | ark::dataset::FLAG_CHOQUE_EN_HUECO | ark::dataset::FLAG_PAUSA_UMBRAL | ark::dataset::FLAG_CORRUPCION_CAIDA)) != 0u);",
          "const bool veto_ito_z2 = false;")
     ]),
    ("MN1.8", "ignorar llave de amplitud |Z_hueco| > z_max_ef en puente_historico_ito.hpp", ITO,
     [
         ("if (std::abs(out.z_hueco) > out.z_max_efectivo) {",
          "if (false && std::abs(out.z_hueco) > out.z_max_efectivo) {")
     ]),
    ("MN1.9", "colapsar varianza entre imputaciones B_M = 0 en combinar_rubin (imputacion unica)", ITO,
     [
         ("out.b_entre = ss_entre / static_cast<double>(validas - 1);",
          "(void)ss_entre; out.b_entre = 0.0;")
     ]),
    ("MN1.10", "inyectar deriva cero (0.0) en valor_z0_imputado bajo FLUJO cuando mu != 0", ITO,
     [
         ("const double reposo_flujo = std::isfinite(mu_calma) ? mu_calma : 0.0;",
          "(void)mu_calma; const double reposo_flujo = 0.0;")
     ]),
]


def compilar_uno(idx, mid, desc, rel, cambios, scratch):
    d = os.path.join(scratch, "mf4_%d" % idx)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(INC_DIR, os.path.join(d, "inc"))
    os.makedirs(os.path.join(d, "nivel1_nan_huecos"), exist_ok=True)
    for fn in ("protocolo_inputs_lilyth.hpp", "banco_sintetico_comun.hpp", "test_fase4_ciego_rubin.cpp"):
        shutil.copy(os.path.join(AQUI, fn), os.path.join(d, "nivel1_nan_huecos", fn))

    if rel is not None:
        f = os.path.join(d, rel)
        s = open(f, encoding="utf-8").read()
        for viejo, nuevo in cambios:
            n = s.count(viejo)
            if n != 1:
                return mid, desc, None, "MAL DEFINIDO (%d apariciones)" % n
            s = s.replace(viejo, nuevo)
        open(f, "w", encoding="utf-8", newline="\n").write(s)

    src = os.path.join(d, "nivel1_nan_huecos", "test_fase4_ciego_rubin.cpp")
    exe = os.path.join(d, "test_f4_ciego.exe")
    bat = os.path.join(d, "c.bat")
    with open(bat, "w") as b:
        b.write('@echo off\r\ncall "%s" >nul 2>&1\r\ncd /d "%s"\r\n' % (VCVARS, d))
        b.write('cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%s" "%s" /Fe"%s" /link /STACK:1073741824\r\n'
                % (os.path.join(d, "inc"), src, exe))
    c = subprocess.run(["cmd", "/c", bat], capture_output=True, text=True, encoding="utf-8", errors="replace")
    if not os.path.exists(exe):
        err = [l for l in (c.stdout + c.stderr).splitlines() if "error" in l][:2]
        return mid, desc, None, "NO COMPILA: %s" % err
    warns = [l for l in (c.stdout + c.stderr).splitlines() if "warning C" in l]
    return mid, desc, exe, ("COMPILA LIMPIO (/W4, 0 warnings)" if not warns else "COMPILA (%d warnings)" % len(warns))


def ejecutar_uno(mid, desc, exe, estado_comp):
    if exe is None:
        return mid, desc, estado_comp, False
    r = subprocess.run([exe, "200", os.path.join(AQUI, "volcado_orig.txt")],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    fallas = [l for l in r.stdout.splitlines() if l.startswith("[FALLA]")]
    if mid == "MF4.0":
        ok_testigo = (r.returncode == 0) and (len(fallas) == 0) and ("FALLAS FASE 4: 0" in r.stdout)
        return mid, desc, ("TESTIGO PASA (0 fallas)" if ok_testigo else "TESTIGO FALLA: %s" % fallas), ok_testigo

    cazado = bool(fallas) or (r.returncode != 0) or ("FALLAS FASE 4: 0" not in r.stdout)
    detalle = "; ".join(re.sub(r"^\[FALLA\] ", "", l) for l in fallas) or ("rc %d" % r.returncode)
    return mid, desc, ("CAZADO: " + detalle if cazado else "VIVO"), cazado


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    solo_compilar = "--solo-compilar" in sys.argv
    args_pos = [a for a in sys.argv[1:] if a != "--solo-compilar"]
    scratch = args_pos[0] if args_pos else os.path.join(os.environ.get("TEMP", r"C:\Temp"), "ark_mut_f4_sellado")
    os.makedirs(scratch, exist_ok=True)
    print("====================================================================================================")
    print(" FASE 4 SELLADA: CAZA DE MUTANTES (MF4.0 + MN1.7-MN1.10, R=200 fijo)")
    print("====================================================================================================")

    # Compilacion en paralelo (los 5 ejecutables en ~5s)
    with cf.ThreadPoolExecutor(5) as ex:
        futs = [ex.submit(compilar_uno, i, *m, scratch) for i, m in enumerate(MUTANTES)]
        compilados = [f.result() for f in futs]

    if solo_compilar:
        todos_ok = True
        for mid, desc, exe, estado in compilados:
            print("%-6s %-76s %s" % (mid, desc, estado))
            if exe is None:
                todos_ok = False
        print("====================================================================================================")
        sys.exit(0 if todos_ok else 1)

    # Ejecucion secuencial para que cada corrida R=200 aproveche 100% de los 32 hilos sin contencion
    res = [ejecutar_uno(mid, desc, exe, estado) for (mid, desc, exe, estado) in compilados]

    cazados = 0
    for mid, desc, estado, ok in res:
        print("%-6s %-76s %s" % (mid, desc, estado))
        if mid != "MF4.0" and ok:
            cazados += 1
    total_mut = len(MUTANTES) - 1
    testigo_ok = res[0][3]
    print("----------------------------------------------------------------------------------------------------")
    print("RESULTADO FASE 4 MUTANTES: %d de %d mutantes cazados | Testigo MF4.0: %s (%d de %d controles)"
          % (cazados, total_mut, "PASA" if testigo_ok else "FALLA",
             cazados + (1 if testigo_ok else 0), len(MUTANTES)))
    print("====================================================================================================")
    sys.exit(0 if (testigo_ok and cazados == total_mut) else 1)


if __name__ == "__main__":
    main()
