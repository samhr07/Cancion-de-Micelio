# -*- coding: utf-8 -*-
"""Fase 4 (Caza de mutantes nuevos sobre puente de Ito y lectura de etiquetas FLAG_PUENTE_ITO).

Complementa a mutaciones_nivel1.py (MN1.0-MN1.6) sin alterar el contrato sellado de CN1.6,
verificando en copia aislada de inc/ y nivel1_nan_huecos/ que MN1.7 y MN1.8 quedan cazados.
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
    ("MF4.0", "testigo: sin cambio (Fase 4)", None, None),
    ("MN1.7", "ignorar FLAG_PUENTE_ITO en Z2 y Z3 dentro de observador_lilyth.hpp", OBS,
     [
         ("const bool veto_z3_ito = ((flags_hueco & ark::dataset::FLAG_PUENTE_ITO) != 0u);",
          "const bool veto_z3_ito = false;"),
         ("const bool veto_ito_z2 = tiene_idx && ((hist_flags_hueco_[idx] & (ark::dataset::FLAG_PUENTE_ITO | ark::dataset::FLAG_CHOQUE_EN_HUECO | ark::dataset::FLAG_PAUSA_UMBRAL | ark::dataset::FLAG_CORRUPCION_CAIDA)) != 0u);",
          "const bool veto_ito_z2 = false;")
     ]),
    ("MN1.8", "ignorar llave de amplitud |Z_hueco| > z_max en puente_historico_ito.hpp", ITO,
     [
         ("if (std::abs(out.z_hueco) > cfg.z_max_hueco) {",
          "if (false && std::abs(out.z_hueco) > cfg.z_max_hueco) {")
     ]),
]


def correr_uno(idx, mid, desc, rel, cambios, scratch):
    d = os.path.join(scratch, "mf4_%d" % idx)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(INC_DIR, os.path.join(d, "inc"))
    os.makedirs(os.path.join(d, "nivel1_nan_huecos"), exist_ok=True)
    for fn in ("protocolo_inputs_lilyth.hpp", "banco_sintetico_comun.hpp", "test_comparacion_puentes_nivel1.cpp"):
        shutil.copy(os.path.join(AQUI, fn), os.path.join(d, "nivel1_nan_huecos", fn))

    if rel is not None:
        f = os.path.join(d, rel)
        s = open(f, encoding="utf-8").read()
        for viejo, nuevo in cambios:
            n = s.count(viejo)
            if n != 1:
                return mid, desc, "MAL DEFINIDO (%d apariciones)" % n, False
            s = s.replace(viejo, nuevo)
        open(f, "w", encoding="utf-8", newline="\n").write(s)

    src = os.path.join(d, "nivel1_nan_huecos", "test_comparacion_puentes_nivel1.cpp")
    exe = os.path.join(d, "test_f4.exe")
    bat = os.path.join(d, "c.bat")
    with open(bat, "w") as b:
        b.write('@echo off\r\ncall "%s" >nul 2>&1\r\ncd /d "%s"\r\n' % (VCVARS, d))
        b.write('cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%s" "%s" /Fe"%s" /link /STACK:1073741824\r\n'
                % (os.path.join(d, "inc"), src, exe))
    c = subprocess.run(["cmd", "/c", bat], capture_output=True, text=True, encoding="utf-8", errors="replace")
    if not os.path.exists(exe):
        err = [l for l in (c.stdout + c.stderr).splitlines() if "error" in l][:2]
        return mid, desc, "NO COMPILA: %s" % err, False

    r = subprocess.run([exe, "200", os.path.join(AQUI, "volcado_orig.txt")], capture_output=True, text=True, encoding="utf-8", errors="replace")
    fallas = [l for l in r.stdout.splitlines() if l.startswith("[FALLA]")]
    if mid == "MF4.0":
        ok_testigo = (r.returncode == 0) and (len(fallas) == 0) and ("FALLAS FASE 4: 0" in r.stdout)
        return mid, desc, ("TESTIGO PASA (0 fallas)" if ok_testigo else "TESTIGO FALLA: %s" % fallas), ok_testigo

    cazado = bool(fallas) or (r.returncode != 0) or ("FALLAS FASE 4: 0" not in r.stdout)
    detalle = "; ".join(re.sub(r"^\[FALLA\] ", "", l) for l in fallas) or ("rc %d" % r.returncode)
    return mid, desc, ("CAZADO: " + detalle if cazado else "VIVO"), cazado


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    scratch = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ.get("TEMP", r"C:\Temp"), "ark_mut_f4")
    os.makedirs(scratch, exist_ok=True)
    print("====================================================================================================")
    print(" FASE 4: CAZA DE MUTANTES DE ETIQUETAS E ITO (mutaciones_fase4_puente.py)")
    print("====================================================================================================")
    with cf.ThreadPoolExecutor(3) as ex:
        futs = [ex.submit(correr_uno, i, *m, scratch) for i, m in enumerate(MUTANTES)]
        res = [f.result() for f in futs]

    cazados = 0
    for mid, desc, estado, ok in res:
        print("%-6s %-66s %s" % (mid, desc, estado))
        if mid != "MF4.0" and ok:
            cazados += 1
    total = len(MUTANTES) - 1
    print("----------------------------------------------------------------------------------------------------")
    print("RESULTADO FASE 4 MUTANTES: %d de %d cazados" % (cazados, total))
    print("====================================================================================================")
    sys.exit(0 if (res[0][3] and cazados == total) else 1)


if __name__ == "__main__":
    main()
