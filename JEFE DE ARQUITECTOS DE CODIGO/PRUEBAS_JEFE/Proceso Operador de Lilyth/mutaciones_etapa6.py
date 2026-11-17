#!/usr/bin/env python3
"""
mutaciones_etapa6.py -- Bateria de mutantes de la Etapa 6 del Operador de Lilyth (`reloj_marea_lilyth.hpp`).
Regla 4 del plano de Code (2026-11-06):
  - Primero verifica el testigo sin mutar (debe compilar sin avisos y pasar con codigo 0).
  - Cada mutante modifica `include/ark/lilyth/reloj_marea_lilyth.hpp`, DEBE compilar limpio bajo `-Werror`,
    y la bateria `etapa6_lilyth.cpp` debe fallar (codigo != 0). Un mutante que no compila NO esta cazado.
"""

import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
HEADER = ROOT / "include" / "ark" / "lilyth" / "reloj_marea_lilyth.hpp"
TEST_CPP = ROOT / "etapa6_lilyth.cpp"
BIN = ROOT / "etapa6_lilyth_mut_bin"

CXX_CMD = [
    "g++",
    "-std=c++20",
    "-O2",
    "-Wall",
    "-Wextra",
    "-Wpedantic",
    "-Wno-mismatched-new-delete",
    "-Werror",
    "-pthread",
    f"-I{ROOT / 'include'}",
    str(TEST_CPP),
    "-o",
    str(BIN),
]

MUTANTES = [
    (
        "M01_nu_punto_sin_congelar",
        "const bool congelar_nu_punto_en_apertura = true;",
        "const bool congelar_nu_punto_en_apertura = false;",
        "No congela nu_punto* en W_ancla al abrir (recalcula con la ventana movil posterior al choque)",
    ),
    (
        "M02_nu_info_shannon_aceptado",
        "const bool aceptar_shannon = false;",
        "const bool aceptar_shannon = true;",
        "Acepta NU_INFO_SHANNON como reloj de actividad en configurar/calibrar (violando D7 y E1.11)",
    ),
    (
        "M03_tau_act_en_pasos_de_pared",
        "const bool usar_pasos_de_pared_en_tau_act = false;",
        "const bool usar_pasos_de_pared_en_tau_act = true;",
        "Calcula tau_act en pasos de pared en vez de Delta_nu_rest / nu_punto* (violando E1.11)",
    ),
    (
        "M04_marea_rama_cruzada",
        "const bool parear_mismo_evento = true;",
        "const bool parear_mismo_evento = false;",
        "Calcula Delta_tau_marea con la rama cruzada (cruda de k contra purgada de otro evento k' != k)",
    ),
    (
        "M05_marea_signo_invertido",
        "const bool invertir_signo_marea = false;",
        "const bool invertir_signo_marea = true;",
        "Invierte el signo de Delta_tau_marea y Delta_RMST_marea (purgada - cruda en vez de cruda - purgada)",
    ),
    (
        "M06_marea_no_descriptiva",
        "res.descriptivo_fuera_de_holm = true;",
        "res.descriptivo_fuera_de_holm = false;",
        "Marca el contraste de marea como no descriptivo / dentro de Holm (violando la decision del Comandante)",
    ),
]


def compilar_y_correr(modo_rapido: bool) -> tuple[bool, bool, str]:
    comp = subprocess.run(CXX_CMD, capture_output=True, text=True)
    if comp.returncode != 0:
        return False, False, comp.stderr
    cmd = [str(BIN)] + (["--rapido"] if modo_rapido else [])
    run = subprocess.run(cmd, capture_output=True, text=True)
    return True, (run.returncode == 0), (run.stdout + "\n" + run.stderr)


def main() -> int:
    original = HEADER.read_text(encoding="utf-8")
    try:
        print("====================================================================")
        print(" AUDITORIA DE MUTANTES ETAPA 6 -- OPERADOR DE LILYTH")
        print("====================================================================")
        comp_ok, test_ok, out = compilar_y_correr(modo_rapido=False)
        if not comp_ok or not test_ok:
            print("[FALLA FATAL] El testigo sin mutar no paso:")
            print(out)
            return 1
        print("[OK] Testigo sin mutar compila limpio y pasa todos los controles (modo completo).")

        cazados = 0
        for mid, target, replacement, desc in MUTANTES:
            if original.count(target) != 1:
                print(f"[FALLA] {mid}: patron objetivo no es unico ({original.count(target)} ocurrencias)")
                return 1
            mutado = original.replace(target, replacement, 1)
            HEADER.write_text(mutado, encoding="utf-8")

            comp_m, pasa_m, _ = compilar_y_correr(modo_rapido=True)
            if not comp_m:
                print(f"[FALLA] {mid}: NO COMPILO (un mutante que no compila no esta cazado)")
                return 1
            if not pasa_m:
                cazados += 1
                print(f"  [CAZADO] {mid}: {desc}")
            else:
                comp_f, pasa_f, _ = compilar_y_correr(modo_rapido=False)
                if comp_f and not pasa_f:
                    cazados += 1
                    print(f"  [CAZADO] {mid}: {desc}")
                else:
                    print(f"  [ VIVO ] {mid}: SOBREVIVIO A LA BATERIA ({desc})")
                    return 1

        print("====================================================================")
        print(f" RESULTADO FINAL: {cazados}/{len(MUTANTES)} MUTANTES COMPILARON Y FUERON CAZADOS")
        print("====================================================================")
        return 0
    finally:
        HEADER.write_text(original, encoding="utf-8")
        if BIN.exists():
            BIN.unlink()


if __name__ == "__main__":
    sys.exit(main())