#!/usr/bin/env python3
"""
mutaciones_etapa5.py -- Bateria de mutantes de la Etapa 5 del Operador de Lilyth (`fin_calor_lilyth.hpp`).
Regla 4 del plano de Code (2026-11-06):
  - Primero verifica el testigo sin mutar (debe compilar sin avisos y pasar con codigo 0).
  - Cada mutante modifica `include/ark/lilyth/fin_calor_lilyth.hpp`, DEBE compilar limpio bajo `-Werror`,
    y la bateria `etapa5_lilyth.cpp` debe fallar (codigo != 0). Un mutante que no compila NO esta cazado.
"""

import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
HEADER = ROOT / "include" / "ark" / "lilyth" / "fin_calor_lilyth.hpp"
TEST_CPP = ROOT / "etapa5_lilyth.cpp"
BIN = ROOT / "etapa5_lilyth_mut_bin"

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
        "M01_hipotesis_nula_invertida_f4",
        "const bool hipotesis_nula_invertida_f4 = false;",
        "const bool hipotesis_nula_invertida_f4 = true;",
        "Invierte la hipotesis nula (defecto F4: evalua contra mu_0,beta en vez de mu_calor)",
    ),
    (
        "M02_mu_calor_desde_amplitud_alarma",
        "const bool mu_calor_desde_amplitud = false;",
        "const bool mu_calor_desde_amplitud = true;",
        "Fija mu_calor desde la amplitud de la alarma en vez de las l* ventanas posteriores (violando E1.7)",
    ),
    (
        "M03_arl_sin_dividir_por_alfa",
        "const bool dividir_arl_por_alfa = true;",
        "const bool dividir_arl_por_alfa = false;",
        "Calibra h con ARL_0,caliente = T_max en vez de T_max / alfa (violando E1.7)",
    ),
    (
        "M04_usar_h_extrapolado",
        "const bool permitir_h_extrapolado = false;",
        "const bool permitir_h_extrapolado = true;",
        "Acepta y usa un umbral h extrapolado por buscar_h_con en vez de declararlo no calibrado",
    ),
    (
        "M05_terminar_calor_con_primera_escala",
        "const bool terminar_con_primera_escala = false;",
        "const bool terminar_con_primera_escala = true;",
        "Termina el calor global cuando termina la primera escala en vez de esperar a que terminen todas",
    ),
    (
        "M06_signo_inicio_calor_invertido",
        "if (!est.caliente && st.alarma.hay && st.alarma.signo < 0) {",
        "if (!est.caliente && st.alarma.hay && st.alarma.signo > 0) {",
        "Usa signo > 0 (enfriamiento) en vez de signo < 0 (subida de temperatura = caida de log_beta) para iniciar calor",
    ),
    (
        "M07_declarar_fin_sin_esperar_l_star",
        "if (est.ventanas_mu_calor >= est.l_star) {",
        "if (est.ventanas_mu_calor >= 1) {",
        "Cierra mu_calor en 1 sola ventana sin esperar las l* ventanas posteriores a la alarma",
    ),
    (
        "M08_estadistico_fin_calor_hacia_abajo",
        "const double g = (S > 0.0) ? (S * S) / (2.0 * static_cast<double>(l) * tau_int_ * factor_mu) : 0.0;",
        "const double g = (S < 0.0) ? (S * S) / (2.0 * static_cast<double>(l) * tau_int_ * factor_mu) : 0.0;",
        "Evalua el GLR unilateral en direccion opuesta (S < 0 en vez de S > 0 hacia mu_0,beta)",
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
        print(" AUDITORIA DE MUTANTES ETAPA 5 -- OPERADOR DE LILYTH")
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