#!/usr/bin/env python3
"""
mutaciones_etapa4.py -- Bateria de mutantes de la Etapa 4 del Operador de Lilyth (`armonico_lilyth.hpp`).
Regla 4 del plano de Code (2026-11-06):
  - Primero verifica el testigo sin mutar (debe compilar sin avisos y pasar con codigo 0).
  - Cada mutante modifica `include/ark/lilyth/armonico_lilyth.hpp`, DEBE compilar limpio bajo `-Werror`,
    y la bateria `etapa4_lilyth.cpp` debe fallar (codigo != 0). Un mutante que no compila NO esta cazado.
"""

import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
HEADER = ROOT / "include" / "ark" / "lilyth" / "armonico_lilyth.hpp"
TEST_CPP = ROOT / "etapa4_lilyth.cpp"
BIN = ROOT / "etapa4_lilyth_mut_bin"

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
        "M01_signo_discriminante",
        "if (disc < 0.0 && a2 < 0.0) {",
        "if (disc > 0.0 && a2 < 0.0) {",
        "Invierte el signo del discriminante (a1^2 + 4*a2 > 0 en vez de < 0) al clasificar raices complejas",
    ),
    (
        "M02_recortar_periodo",
        "const bool recortar_periodo_corto = false;",
        "const bool recortar_periodo_corto = true;",
        "Recorta y descarta periodos <= 4 barras en Prony (violando E1.8)",
    ),
    (
        "M03_prony_desde_tau_1",
        "const std::size_t tau_ini = 2;",
        "const std::size_t tau_ini = 1;",
        "Ajusta la recurrencia de Prony desde tau = 1 en vez de tau = 2",
    ),
    (
        "M04_reescalada_incluye_no_restituidos",
        "const bool excluir_no_restituidos = true;",
        "const bool excluir_no_restituidos = false;",
        "Incluye los eventos no restituidos (tau_rest = inf) en PRONY_FASE_REESCALADA (violando E1.20)",
    ),
    (
        "M05_holm_sobre_tres_p",
        "const std::size_t n_familia_holm = 4;",
        "const std::size_t n_familia_holm = 3;",
        "Aplica ark::holm sobre 3 valores p en vez de la familia de 4 hipotesis (violando E1.16)",
    ),
    (
        "M06_omitir_p_nan_en_holm",
        "const bool sanear_nan_antes_de_holm = false;",
        "const bool sanear_nan_antes_de_holm = true;",
        "Sanea/omite un p = NaN a 1.0 antes de ark::holm en vez de propagar la falla dura",
    ),
    (
        "M07_etiqueta_ignora_raices_complejas",
        "return pr.valido && pr.raices_complejas && (min_eta_desde_2 < cuantil_alfa_nulo);",
        "return pr.valido && (min_eta_desde_2 < cuantil_alfa_nulo);",
        "Etiqueta ARMONICO_SOSPECHOSO sin exigir raices complejas en el Prony individual del evento",
    ),
    (
        "M08_sobrepaso_desde_tau_1",
        "// Estadistico de sobrepaso: min de eta(tau) para tau en [2, H]\n    for (std::size_t tau = 2; tau <= H; ++tau) {",
        "// Estadistico de sobrepaso: min de eta(tau) para tau en [2, H]\n    for (std::size_t tau = 1; tau <= H; ++tau) {",
        "Calcula min_eta_desde_2 incluyendo tau = 1 en vez de [2, T_max]",
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
        print(" AUDITORIA DE MUTANTES ETAPA 4 -- OPERADOR DE LILYTH")
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

            # Los primeros 6 mutantes se cazan con el modo rapido (<1s); M07 requiere el control estadistico AR(1) completo
            usar_rapido = True
            comp_m, pasa_m, _ = compilar_y_correr(modo_rapido=usar_rapido)
            if not comp_m:
                print(f"[FALLA] {mid}: NO COMPILO (un mutante que no compila no esta cazado)")
                return 1
            if not pasa_m:
                cazados += 1
                print(f"  [CAZADO] {mid}: {desc}")
            else:
                # Si paso en modo rapido, correr la bateria completa por seguridad
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