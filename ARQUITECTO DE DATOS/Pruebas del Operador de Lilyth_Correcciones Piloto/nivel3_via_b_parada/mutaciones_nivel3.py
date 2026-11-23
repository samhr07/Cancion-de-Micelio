# -*- coding: utf-8 -*-
"""Nivel 3 Sellado (EXPECTATIVAS_NIVEL3.md §3):
Caza de mutantes reales (MN3.0 + M_N3_1, M_N3_2a, M_N3_2b) sobre parada temprana V4,
retirada del veto asimetrico hay_previa y ventana derivada t_max_ef en nulo_lilyth.hpp.

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

NUL = os.path.join("inc", "ark", "lilyth", "nulo_lilyth.hpp")

CONTROLES = {
    "NV3.1": ("test_nv3_1_banco_r1.cpp", "test_nv3_1_banco_r1.exe", "FALLAS NV3.1: 0"),
    "NV3.2": ("test_nv3_2_banco_r2.cpp", "test_nv3_2_banco_r2.exe", "FALLAS NV3.2: 0"),
    "NV3.3": ("test_nv3_3_integridad.cpp", "test_nv3_3_integridad.exe", "FALLAS NV3.3: 0"),
}

MUTANTES = [
    ("MN3.0", "testigo: sin cambio (Nivel 3 sellado)", None,
     ["NV3.1", "NV3.2", "NV3.3"], [], []),
    ("M_N3_1", "desactivar parada temprana V4 en nulo_lilyth.hpp", NUL,
     ["NV3.1"], ["NV3.2"],
     [
         ("    const bool parada_temprana_v4 = usar_enmienda_nivel3;",
          "    bool parada_temprana_v4 = false;")
     ]),
    ("M_N3_2a", "reactivar veto hay_previa en nulo_lilyth.hpp", NUL,
     ["NV3.1", "NV3.2"], [],
     [
         ("        const bool vetar_por_ajena = !usar_enmienda_nivel3;",
          "        bool vetar_por_ajena = true;")
     ]),
    ("M_N3_2b", "usar T_MAX fijo en lugar de t_max_ef en construir_ventana_inyeccion", NUL,
     ["NV3.3"], [],
     [
         ("    const std::size_t t_max_ventana = usar_enmienda_nivel3 ? t_max_ef : T_MAX;",
          "    const std::size_t t_max_ventana = T_MAX;")
     ]),
]


def _extraer_metricas(ctrl, stdout):
    """Extrae metricas clave impresas por cada control para diagnostico veridico."""
    tags = []
    if ctrl in ("NV3.1", "NV3.2"):
        m_z = re.search(r"z_S\s*=\s*([+\-0-9.]+)", stdout)
        m_par = re.search(r"paradas_tempranas_b_total\s*=\s*(\d+)", stdout)
        m_sol = re.search(r"solapes_b_total\s*=\s*(\d+)", stdout)
        m_sneg = re.search(r"semillas con S < 0\s*=\s*(\d+\s*/\s*\d+)", stdout)
        if m_z:
            tags.append("z_S=%s" % m_z.group(1))
        if m_sneg:
            tags.append("S<0=%s" % m_sneg.group(1).replace(" ", ""))
        if m_par:
            tags.append("parB=%s" % m_par.group(1))
        if m_sol:
            tags.append("solB=%s" % m_sol.group(1))
    elif ctrl == "NV3.3":
        for line in stdout.splitlines():
            if "[NV3.3a]" in line:
                tags.append(line.strip())
                break
    return ", ".join(tags)


def compilar_caso(idx, mid, desc, rel, req_ctrls, rep_ctrls, cambios, scratch):
    d = os.path.join(scratch, "mn3_%d_%s" % (idx, mid.lower()))
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(INC_DIR, os.path.join(d, "inc"))
    sub_n3 = os.path.join(d, "nivel3_via_b_parada")
    os.makedirs(sub_n3, exist_ok=True)

    ctrls = req_ctrls + [c for c in rep_ctrls if c not in req_ctrls]
    for ctrl in ctrls:
        cpp_name = CONTROLES[ctrl][0]
        shutil.copy(os.path.join(AQUI, cpp_name), os.path.join(sub_n3, cpp_name))

    if rel is not None:
        f = os.path.join(d, rel)
        s = open(f, encoding="utf-8").read()
        for viejo, nuevo in cambios:
            n = s.count(viejo)
            if n != 1:
                return mid, desc, rel, req_ctrls, rep_ctrls, {}, "MAL DEFINIDO (%d apariciones)" % n
            s = s.replace(viejo, nuevo)
        open(f, "w", encoding="utf-8", newline="\n").write(s)

    bat = os.path.join(d, "c.bat")
    with open(bat, "w", encoding="utf-8") as b:
        b.write('@echo off\r\ncall "%s" >nul 2>&1\r\ncd /d "%s"\r\n' % (VCVARS, d))
        for ctrl in ctrls:
            cpp_name, exe_name, _ = CONTROLES[ctrl]
            src = os.path.join(sub_n3, cpp_name)
            exe = os.path.join(d, exe_name)
            b.write('cl /nologo /std:c++20 /O2 /fp:precise /MD /W4 /permissive- /EHsc /utf-8 /DNOMINMAX /I"%s" "%s" /Fo"%s\\\\" /Fe"%s" /link /STACK:1073741824 || exit /b 1\r\n'
                    % (os.path.join(d, "inc"), src, d, exe))

    c = subprocess.run(["cmd", "/c", bat], capture_output=True, text=True, encoding="utf-8", errors="replace")
    txt_c = (c.stdout or "") + "\n" + (c.stderr or "")
    warns = [l for l in txt_c.splitlines() if "warning C" in l]
    exes = {ctrl: os.path.join(d, CONTROLES[ctrl][1]) for ctrl in ctrls}
    faltan = [ctrl for ctrl, p in exes.items() if not os.path.exists(p)]

    if c.returncode != 0 or faltan:
        err = [l for l in txt_c.splitlines() if "error" in l][:2]
        return mid, desc, rel, req_ctrls, rep_ctrls, {}, "NO COMPILA (%s): %s" % (faltan, err)
    if warns:
        return mid, desc, rel, req_ctrls, rep_ctrls, {}, "COMPILA CON WARNINGS (/W4, %d): %s" % (len(warns), warns[0])
    return mid, desc, rel, req_ctrls, rep_ctrls, exes, "COMPILA LIMPIO (/W4, 0 warnings; %s)" % ",".join(ctrls)


def _cazado_por(ctrl, rc, stdout):
    fallas = [l.strip() for l in stdout.splitlines() if "[FALLA]" in l]
    if rc == 0:
        return False, fallas
    if ctrl == "NV3.1":
        return any("NV3.1" in l for l in fallas), fallas
    if ctrl == "NV3.2":
        return any("NV3.2" in l for l in fallas), fallas
    if ctrl == "NV3.3":
        return any("NV3.3a" in l for l in fallas) and ("FALLAS NV3.3: 0" not in stdout), fallas
    return False, fallas


def ejecutar_caso(comp_res):
    mid, desc, _rel, req_ctrls, rep_ctrls, exes, estado_comp = comp_res
    if not exes:
        return mid, desc, estado_comp, False, []

    salidas = {}
    ctrls = req_ctrls + [c for c in rep_ctrls if c not in req_ctrls]
    for ctrl in ctrls:
        r = subprocess.run([exes[ctrl]], capture_output=True, text=True, encoding="utf-8", errors="replace")
        salidas[ctrl] = (r.returncode, (r.stdout or "").replace("\r\n", "\n"))

    detalles = []
    if mid == "MN3.0":
        ok_all = True
        for ctrl in req_ctrls:
            rc, out = salidas[ctrl]
            pass_tok = CONTROLES[ctrl][2]
            fallas = [l.strip() for l in out.splitlines() if "[FALLA]" in l]
            ok_c = (rc == 0) and (len(fallas) == 0) and (pass_tok in out)
            ok_all = ok_all and ok_c
            detalles.append("%s: %s (%s)" % (ctrl, "PASA" if ok_c else "FALLA", _extraer_metricas(ctrl, out) or pass_tok))
        msg = "TESTIGO PASA (0 fallas)" if ok_all else ("TESTIGO FALLA: " + " | ".join(detalles))
        return mid, desc, msg, ok_all, detalles

    ok_kill_all = True
    cazadores_ok = []
    fallos_motivo = []
    for ctrl in req_ctrls:
        rc, out = salidas[ctrl]
        killed, fallas = _cazado_por(ctrl, rc, out)
        met = _extraer_metricas(ctrl, out)
        subs = "; ".join(re.sub(r"^\[FALLA\]\s*", "", l) for l in fallas) or ("rc=%d" % rc)
        detalles.append("%s[REQ]: killed=%s | %s | fallas=[%s]" % (ctrl, killed, met, subs))
        if killed:
            cazadores_ok.append("%s (%s; %s)" % (ctrl, met, subs))
        else:
            ok_kill_all = False
            fallos_motivo.append("%s=%s" % (ctrl, "VIVO" if rc == 0 and not fallas else subs))

    for ctrl in rep_ctrls:
        rc, out = salidas[ctrl]
        killed, fallas = _cazado_por(ctrl, rc, out)
        met = _extraer_metricas(ctrl, out)
        subs = "; ".join(re.sub(r"^\[FALLA\]\s*", "", l) for l in fallas) or "PASA"
        detalles.append("%s[REP]: %s | %s | [%s]" % (ctrl, "FALLA" if killed else "PASA", met, subs))

    if ok_kill_all:
        msg = "CAZADO (%s): %s" % ("+".join(req_ctrls), " || ".join(cazadores_ok))
        return mid, desc, msg, True, detalles
    msg = "NO CAZADO (%s): %s" % ("+".join(req_ctrls), " || ".join(fallos_motivo))
    return mid, desc, msg, False, detalles


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    solo_compilar = "--solo-compilar" in sys.argv
    args_pos = [a for a in sys.argv[1:] if a != "--solo-compilar"]
    scratch = args_pos[0] if args_pos else os.path.join(os.environ.get("TEMP", r"C:\Temp"), "ark_mut_nivel3_sellado")
    os.makedirs(scratch, exist_ok=True)

    print("====================================================================================================")
    print(" NIVEL 3 SELLADO: CAZA DE MUTANTES (MN3.0 + M_N3_1, M_N3_2a, M_N3_2b | 0 warnings /W4)")
    print("====================================================================================================")

    with cf.ThreadPoolExecutor(4) as ex:
        futs = [ex.submit(compilar_caso, i, *m, scratch) for i, m in enumerate(MUTANTES)]
        compilados = [f.result() for f in futs]

    if solo_compilar:
        todos_ok = True
        for mid, desc, _rel, _req, _rep, exes, estado in compilados:
            print("%-8s %-64s %s" % (mid, desc, estado))
            if not exes:
                todos_ok = False
        print("====================================================================================================")
        sys.exit(0 if todos_ok else 1)

    res = [ejecutar_caso(c) for c in compilados]

    cazados = 0
    for mid, desc, msg, ok, detalles in res:
        print("%-8s %-64s %s" % (mid, desc, msg))
        for d in detalles:
            print("         -> %s" % d)
        if mid != "MN3.0" and ok:
            cazados += 1

    total_mut = len(MUTANTES) - 1
    testigo_ok = res[0][3]
    print("----------------------------------------------------------------------------------------------------")
    print("RESULTADO NIVEL 3 MUTANTES: %d de %d mutantes cazados | Testigo MN3.0: %s (%d de %d controles)"
          % (cazados, total_mut, "PASA" if testigo_ok else "FALLA",
             cazados + (1 if testigo_ok else 0), len(MUTANTES)))
    print("====================================================================================================")
    sys.exit(0 if (testigo_ok and cazados == total_mut) else 1)


if __name__ == "__main__":
    main()
