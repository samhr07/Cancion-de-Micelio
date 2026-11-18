# -*- coding: utf-8 -*-
"""Punto 7 (CN1.6): Caza de mutantes reales por reemplazo de texto sobre copia de cabecera (Nivel 1).

Sigue exactamente la arquitectura de mutaciones_horizonte.py y mutaciones_bifurcacion.py:
copia inc/ y nivel1_nan_huecos/ a una ruta corta temporal, aplica un unico reemplazo exacto
en la cabecera C++, compila test_cota_corrupcion_nivel1.cpp y verifica que las condiciones
CN1.1-CN1.5 cacen al mutante.
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
PRO = os.path.join("nivel1_nan_huecos", "protocolo_inputs_lilyth.hpp")

MUTANTES = [
    ("MN1.0", "testigo: sin cambio", None, None, None),
    ("MN1.1", "MERCADO_QUIETO tratado como invalido en NIVEL", PRO,
     "out.x_para_operador = ultimo_x_valido; // dY = 0 exacto en mercado quieto\n                    out.valido_para_operador = true;",
     "out.x_para_operador = ultimo_x_valido;\n                    out.valido_para_operador = false;"),
    ("MN1.2", "comprimir reloj tau en RELOJ_FISICO_TOPE_R (avanzar)", OBS,
     "const std::size_t tau_ef = (cfg_.politica_z3_piloto == ConfiguracionLilyth::PoliticaZ3Piloto::MUTANTE_COMPRIME_TAU)\n                                       ? rama.pasos_validos_z3\n                                       : tau;",
     "const std::size_t tau_ef = rama.pasos_validos_z3;"),
    ("MN1.3", "no excluir CENSURADA_DATO_INVALIDO en evaluar_nulos_representacion", OBS,
     "if (rama.desenlace != DesenlaceRanura::RESTITUIDA && rama.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE) continue;",
     "if (rama.desenlace != DesenlaceRanura::RESTITUIDA && rama.desenlace != DesenlaceRanura::CENSURADA_HORIZONTE && rama.desenlace != DesenlaceRanura::CENSURADA_DATO_INVALIDO) continue;"),
    ("MN1.4", "ignorar FRONTERA_ESTRUCTURAL (valido=true en cierre de sesion)", PRO,
     "if (clase == ClasePasoLilyth::FRONTERA_ESTRUCTURAL) {\n            ++n_frontera_estructural;\n            racha_corrupta_actual = 0;\n            out.x_para_operador = ark::metrologia::NaN;\n            out.valido_para_operador = false;\n            return out;\n        }",
     "if (clase == ClasePasoLilyth::FRONTERA_ESTRUCTURAL) {\n            ++n_frontera_estructural;\n            racha_corrupta_actual = 0;\n            out.x_para_operador = 0.0;\n            out.valido_para_operador = true;\n            return out;\n        }"),
    ("MN1.5", "omision invertida en PoliticaZ3Piloto (RELOJ_FISICO por omision)", OBS,
     "PoliticaZ3Piloto politica_z3_piloto = PoliticaZ3Piloto::CORTA_AL_PRIMERO;\n    std::size_t r_max_z3_piloto = 0;",
     "PoliticaZ3Piloto politica_z3_piloto = PoliticaZ3Piloto::RELOJ_FISICO_TOPE_R;\n    std::size_t r_max_z3_piloto = 3;"),
    ("MN1.6", "ignorar bloqueo de Z2 corrupto en la via (b)", OBS,
     "if (z2_corrupto) {",
     "if (false && z2_corrupto) {"),
]


def correr_uno(idx, mid, desc, rel, viejo, nuevo, scratch):
    d = os.path.join(scratch, "mn1_%d" % idx)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(INC_DIR, os.path.join(d, "inc"))
    os.makedirs(os.path.join(d, "nivel1_nan_huecos"), exist_ok=True)
    for fn in ("protocolo_inputs_lilyth.hpp", "banco_sintetico_comun.hpp", "test_cota_corrupcion_nivel1.cpp"):
        shutil.copy(os.path.join(AQUI, fn), os.path.join(d, "nivel1_nan_huecos", fn))

    if rel is not None:
        f = os.path.join(d, rel)
        s = open(f, encoding="utf-8").read()
        n = s.count(viejo)
        if n != 1:
            return mid, desc, "MAL DEFINIDO (%d apariciones)" % n, False
        open(f, "w", encoding="utf-8", newline="\n").write(s.replace(viejo, nuevo))

    src = os.path.join(d, "nivel1_nan_huecos", "test_cota_corrupcion_nivel1.cpp")
    exe = os.path.join(d, "test_m.exe")
    bat = os.path.join(d, "c.bat")
    with open(bat, "w") as b:
        b.write('@echo off\r\ncall "%s" >nul 2>&1\r\ncd /d "%s"\r\n' % (VCVARS, d))
        b.write('cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%s" "%s" /Fe"%s" /link /STACK:1073741824\r\n'
                % (os.path.join(d, "inc"), src, exe))
    c = subprocess.run(["cmd", "/c", bat], capture_output=True, text=True, encoding="utf-8", errors="replace")
    if not os.path.exists(exe):
        err = [l for l in (c.stdout + c.stderr).splitlines() if "error" in l][:2]
        return mid, desc, "NO COMPILA: %s" % err, False

    r = subprocess.run([exe, "40", os.path.join(AQUI, "volcado_orig.txt")], capture_output=True, text=True, encoding="utf-8", errors="replace")
    fallas = [l for l in r.stdout.splitlines() if l.startswith("[FALLA]")]
    if mid == "MN1.0":
        ok_testigo = (r.returncode == 0) and (len(fallas) == 0) and ("FALLAS DE CONDICIONES: 0" in r.stdout)
        return mid, desc, ("TESTIGO PASA (0 fallas)" if ok_testigo else "TESTIGO FALLA: %s" % fallas), ok_testigo

    cazado = bool(fallas) or (r.returncode != 0) or ("FALLAS DE CONDICIONES: 0" not in r.stdout)
    detalle = "; ".join(re.sub(r"^\[FALLA\] condicion ", "", l) for l in fallas) or ("rc %d" % r.returncode)
    return mid, desc, ("CAZADO: " + detalle if cazado else "VIVO"), cazado


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    scratch = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ.get("TEMP", r"C:\Temp"), "ark_mut_n1")
    os.makedirs(scratch, exist_ok=True)
    print("====================================================================================================")
    print(" PUNTO 7 (CN1.6): CAZA DE MUTANTES REALES SOBRE COPIA DE CABECERA (mutaciones_nivel1.py)")
    print("====================================================================================================")
    with cf.ThreadPoolExecutor(3) as ex:
        futs = [ex.submit(correr_uno, i, *m, scratch) for i, m in enumerate(MUTANTES)]
        res = [f.result() for f in futs]

    cazados = 0
    for mid, desc, estado, ok in res:
        print("%-6s %-58s %s" % (mid, desc, estado))
        if mid != "MN1.0" and ok:
            cazados += 1
    total_mut = len(MUTANTES) - 1
    print("----------------------------------------------------------------------------------------------------")
    print("RESULTADO CN1.6: %d de %d mutantes cazados (condicion sellada: >= 5 de 6)" % (cazados, total_mut))
    print("====================================================================================================")
    sys.exit(0 if cazados >= 5 else 1)


if __name__ == "__main__":
    main()
