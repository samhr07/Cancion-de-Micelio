#!/usr/bin/env python3
# mutaciones_etapa2.py -- Arnes de mutaciones sobre `modelo_lineal_nulo.hpp` (Etapa 2 del Operador de Lilyth).
# Regla de Code:
#   - Los mutantes van en el codigo de produccion (`modelo_lineal_nulo.hpp`), NUNCA en el arnes de pruebas.
#   - Primero corre el testigo sin mutar (debe compilar y pasar con rc == 0).
#   - Un mutante que no compila NO cuenta como cazado (falla el arnes).

import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
HDR = ROOT / "include" / "ark" / "lilyth" / "modelo_lineal_nulo.hpp"
SRC = ROOT / "etapa2_lilyth.cpp"
BIN = ROOT / "etapa2_mut"

MUTANTES = [
    (
        "M1_signo_kappa",
        "const double kappa_k = (rho[k] - num) / v;",
        "const double kappa_k = -(rho[k] - num) / v;",
    ),
    (
        "M2_indice_recursion_desplazado",
        "phi_k[j - 1] = phi_prev[j - 1] - kappa_k * phi_prev[k - j - 1];",
        "phi_k[j - 1] = phi_prev[j - 1] - kappa_k * phi_prev[j - 1];",
    ),
    (
        "M3_sin_factor_1_menos_kappa2",
        "v *= (1.0 - kappa_k * kappa_k);",
        "v *= 1.0;",
    ),
    (
        "M4_filtro_inverso_sin_condiciones_iniciales",
        "return (idx < n_hist) ? historia[idx] : perfil[idx - n_hist];",
        "return (idx < n_hist) ? m.mu : perfil[idx - n_hist];",
    ),
    (
        "M5_coeficientes_d_en_lugar_de_menos_d",
        "ark::metrologia::elw::coeficientes(-m.d, cap_buf, pi_ptr);",
        "ark::metrologia::elw::coeficientes(m.d, cap_buf, pi_ptr);",
    ),
    (
        "M6a_burn_in_cero_en_ajustar",
        "m.burn_in = static_cast<std::size_t>(std::ceil(std::log(eps) / std::log(m.radio_max)));",
        "(void)eps; m.burn_in = 0;",
    ),
    (
        "M6b_burn_in_cero_en_simular",
        "const std::size_t burn = (e_iny != nullptr && t_iny == 0 && l >= n) ? 0 : m.burn_in;",
        "const std::size_t burn = 0;",
    ),
    (
        "M7_residuos_a_traves_de_cortes",
        """    for (const auto& g : segs_modelo) {
        if (g.size() <= m.p) continue;
        for (std::size_t t = m.p; t < g.size(); ++t) {
            double e = g[t] - m.mu;
            for (std::size_t i = 1; i <= m.p; ++i) {
                e -= m.phi[i - 1] * (g[t - i] - m.mu);
            }
            m.residuos.push_back(e);
            suma_res += e;
        }
    }""",
        """    {
        std::vector<double> todo;
        for (const auto& g : segs_modelo) todo.insert(todo.end(), g.begin(), g.end());
        for (std::size_t t = m.p; t < todo.size(); ++t) {
            double e = todo[t] - m.mu;
            for (std::size_t i = 1; i <= m.p; ++i) e -= m.phi[i - 1] * (todo[t - i] - m.mu);
            m.residuos.push_back(e);
            suma_res += e;
        }
    }""",
    ),
    (
        "M8_quenouille_para_en_primer_no_significativo",
        """        if (std::abs(kappa_k) > umbral_q) {
            m.p = k;
            m.phi = phi_k;
            v_en_p = v;
        }""",
        """        if (std::abs(kappa_k) > umbral_q) {
            m.p = k;
            m.phi = phi_k;
            v_en_p = v;
        } else {
            break;
        }""",
    ),
    (
        "M9_sin_caida_generativo_en_fft",
        "m.cayo_a_generativo = true;",
        "m.cayo_a_generativo = false;",
    ),
    (
        "M10_acepta_d_mayor_igual_medio",
        "if (!(std::abs(cfg.d_externa) < 0.5)) {",
        "if (!(std::abs(cfg.d_externa) < 1.5)) {",
    ),
    (
        "M11_sin_guarda_raiz_unitaria_K",
        "if (ti.K > k_max) {",
        "if (false && ti.K > k_max) {",
    ),
]


def compilar_y_correr():
    cmd_comp = [
        "g++",
        "-std=c++20",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-Wno-mismatched-new-delete",
        "-Werror",
        f"-I{ROOT / 'include'}",
        str(SRC),
        "-o",
        str(BIN),
    ]
    c = subprocess.run(cmd_comp, capture_output=True, text=True)
    if c.returncode != 0:
        return False, -1, c.stderr
    r = subprocess.run([str(BIN)], capture_output=True, text=True)
    return True, r.returncode, r.stdout


def main():
    original = HDR.read_text(encoding="utf-8")
    try:
        comp_ok, rc, out = compilar_y_correr()
        if not comp_ok or rc != 0:
            print("[FALLA] El testigo sin mutar no pasa:")
            print(out)
            sys.exit(1)
        print("[OK   ] Testigo sin mutar compila a /Werror y pasa todos los controles (0 fallas)")

        cazados = 0
        for mid, aguja, reemplazo in MUTANTES:
            if aguja not in original:
                print(f"[FALLA] {mid}: no se encontro el patron en modelo_lineal_nulo.hpp")
                sys.exit(1)
            mutado = original.replace(aguja, reemplazo, 1)
            HDR.write_text(mutado, encoding="utf-8")
            comp_ok, rc, out = compilar_y_correr()
            if not comp_ok:
                print(f"[FALLA] {mid}: NO COMPILA (un mutante que no compila no cuenta como cazado):\n{out}")
                sys.exit(1)
            fallidas = [
                line.split()[1]
                for line in out.splitlines()
                if line.startswith("[FALLA]")
            ]
            if rc != 0 and fallidas:
                cazados += 1
                print(f"[CAZADO] {mid:<48s} -> cazado por {', '.join(fallidas)}")
            else:
                print(f"[VIVO  ] {mid:<48s} -> SOBREVIVIO!")
                sys.exit(1)

        print(f"\nPASA: {cazados}/{len(MUTANTES)} mutantes compilados y cazados (testigo pasa).")
    finally:
        HDR.write_text(original, encoding="utf-8")
        if BIN.exists():
            BIN.unlink()


if __name__ == "__main__":
    main()