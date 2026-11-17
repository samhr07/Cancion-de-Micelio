#!/usr/bin/env python3
# mutaciones_etapa3_0.py -- Arnes de mutaciones sobre `observador_lilyth.hpp` (Etapa 1 + Etapa 3.0 E1.21).
# Incluye los 18 mutantes de la Etapa 1 (M1..M18 de REVISION_ETAPA1.md) y los 4 mutantes de la Etapa 3.0 (M19..M22):
#   M19: la continuacion corta las ranuras vivas
#   M20: la continuacion ignora el signo
#   M21: la continuacion ignora si ya restituyo
#   M22: con varias ranuras vivas, compara con la ranura mas antigua en lugar de la mas reciente

import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
HDR = ROOT / "include" / "ark" / "lilyth" / "observador_lilyth.hpp"
SRC = ROOT / "etapa3_0_lilyth.cpp"
BIN = ROOT / "etapa3_0_mut"

MUTANTES = [
    (
        "M01_abre_sin_espejo",
        "if (!espejo_abierto_) {",
        "if (true) {",
    ),
    (
        "M02_guarda_ancla_cero",
        "r.guarda_ancla = static_cast<std::size_t>(std::ceil(ti)) * m;",
        "(void)ti; r.guarda_ancla = 0;",
    ),
    (
        "M03_ancla_hasta_t_cambio",
        "r.t_fin_ancla = r.t_cambio - static_cast<std::int64_t>(r.guarda_ancla) - 1;",
        "r.t_fin_ancla = r.t_cambio;",
    ),
    (
        "M04_w_ancla_con_det_L",
        "r.w_ancla = Rasgos::n_reg * m;",
        "r.w_ancla = f.det.L * m;",
    ),
    (
        "M05_dy_a_la_mitad",
        "const double dy = y_a - med;",
        "const double dy = 0.5 * (y_a - med);",
    ),
    (
        "M06_sin_guarda_signo_dy",
        "if (!(std::isfinite(dy) && dy * static_cast<double>(signo) > 0.0)) {",
        "if (!std::isfinite(dy)) { (void)signo;",
    ),
    (
        "M07_consolidado_ignora_calor",
        "if (!calor_activo_ && !std::isfinite(rama.tau_rest_consolidado)) {",
        "if (!std::isfinite(rama.tau_rest_consolidado)) {",
    ),
    (
        "M08_replica_no_cierra_consolidadas",
        "if (std::isfinite(rama->tau_rest_consolidado)) {\n                rama->desenlace = DesenlaceRanura::RESTITUIDA;\n                rama->tau_fin = static_cast<std::size_t>(rama->tau_rest_consolidado);\n            }",
        "if (std::isfinite(rama->tau_rest_consolidado)) {\n                continue;\n            }",
    ),
    (
        "M09_espejo_purga_sin_historia_calibracion",
        "if (purga_activa_) {\n            bool hp = false;",
        "if (false && purga_activa_) {\n            bool hp = false;",
    ),
    (
        "M10_corte_no_abre_segmento",
        "if (purgado_vivo_) ++seg_purgado_;",
        "if (purgado_vivo_) {}",
    ),
    (
        "M11_origen_nivel_purgado_en_x_t",
        "if (!purgado_vivo_) { nivel_purgado_ = x - y_rep; purgado_vivo_ = true; }",
        "if (!purgado_vivo_) { nivel_purgado_ = x; purgado_vivo_ = true; }",
    ),
    (
        "M12_aalen_johansen_sin_supervivientes_prev",
        "c_rest += s_prev * (static_cast<double>(dr) / denom);",
        "c_rest += (static_cast<double>(dr) / denom);",
    ),
    (
        "M13_tau_fin_replica_uno_de_mas",
        "rama->tau_fin = tau;",
        "rama->tau_fin = tau + 1;",
    ),
    (
        "M14_sin_guarda_reloj_desalineado",
        "if (s.t != t_ + 1) {",
        "if (false && s.t != t_ + 1) {",
    ),
    (
        "M15_horizonte_uno_de_mas",
        "if (tau >= cfg_.t_max_efectivo) {",
        "if (tau > cfg_.t_max_efectivo) {",
    ),
    (
        "M16_sin_guarda_memoria_ancla",
        "if (!indice(t_ini, i0)) return DesenlaceRanura::NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA;",
        "(void)i0;",
    ),
    (
        "M17_purgada_y_cruda_invertidas",
        "preparar_rama(r.purgada, true, medir_purgada(), t_ini, r.t_fin_ancla, yp, okp, r.signo);\n        preparar_rama(r.cruda, false, medir_cruda(), t_ini, r.t_fin_ancla, yc, okc, r.signo);",
        "preparar_rama(r.purgada, false, medir_purgada(), t_ini, r.t_fin_ancla, yc, okc, r.signo);\n        preparar_rama(r.cruda, true, medir_cruda(), t_ini, r.t_fin_ancla, yp, okp, r.signo);",
    ),
    (
        "M18_espejo_cierra_con_eventos_temperatura",
        "if (e.canal != 0) continue;",
        "if (false && e.canal != 0) continue;",
    ),
    (
        "M19_continuacion_corta",
        """                if (mismo_signo && sin_restituir) {
                    out.continuacion = true;
                    ++cont_.continuaciones;
                }""",
        """                if (mismo_signo && sin_restituir) {
                    out.continuacion = true;
                    ++cont_.continuaciones;
                    for (auto& r : ranuras_) {
                        if (!r.ocupada || t_ <= r.t_apertura) continue;
                        cortar_por_replica(r, static_cast<std::size_t>(t_ - r.t_apertura), a_abre.signo);
                        if (!r.alguna_rama_abierta()) cerrar(r, out);
                    }
                }""",
    ),
    (
        "M20_continuacion_ignora_signo",
        "const bool mismo_signo = (a_abre.signo == ref_viva->signo);",
        "const bool mismo_signo = true;",
    ),
    (
        "M21_continuacion_ignora_si_ya_restituyo",
        "const bool sin_restituir = !ref_viva->ha_restituido();",
        "const bool sin_restituir = true;",
    ),
    (
        "M22_compara_con_ranura_mas_antigua",
        "if (ref_viva == nullptr || r.id_evento > ref_viva->id_evento) ref_viva = &r;",
        "if (ref_viva == nullptr || r.id_evento < ref_viva->id_evento) ref_viva = &r;",
    ),
]


import concurrent.futures
import os
import tempfile


def compilar_y_correr(include_override=None, bin_path=BIN):
    cmd_comp = [
        "g++",
        "-std=c++20",
        "-O2",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        "-Wno-mismatched-new-delete",
        "-Werror",
    ]
    if include_override is not None:
        cmd_comp.append(f"-I{include_override}")
    cmd_comp.extend([
        f"-I{ROOT / 'include'}",
        str(SRC),
        "-o",
        str(bin_path),
    ])
    c = subprocess.run(cmd_comp, capture_output=True, text=True)
    if c.returncode != 0:
        return False, -1, c.stderr
    r = subprocess.run([str(bin_path)], capture_output=True, text=True)
    return True, r.returncode, r.stdout


def evaluar_mutante(args):
    idx, mid, aguja, reemplazo, original = args
    if aguja not in original:
        return idx, mid, "NO_PATRON", []
    mutado = original.replace(aguja, reemplazo, 1)
    with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
        tmp = pathlib.Path(tmpdir)
        hdr_dir = tmp / "ark" / "lilyth"
        hdr_dir.mkdir(parents=True, exist_ok=True)
        (hdr_dir / "observador_lilyth.hpp").write_text(mutado, encoding="utf-8")
        bin_mut = tmp / f"etapa3_0_mut_{idx}"
        comp_ok, rc, out = compilar_y_correr(include_override=tmp, bin_path=bin_mut)
        if not comp_ok:
            return idx, mid, f"NO_COMPILA:\n{out}", []
        fallidas = [
            line.split()[1]
            for line in out.splitlines()
            if line.startswith("[FALLA]")
        ]
        if rc != 0 and fallidas:
            return idx, mid, "CAZADO", fallidas
        return idx, mid, "VIVO", []


def main():
    original = HDR.read_text(encoding="utf-8")
    try:
        comp_ok, rc, out = compilar_y_correr()
        if not comp_ok or rc != 0:
            print("[FALLA] El testigo sin mutar no pasa:")
            print(out)
            sys.exit(1)
        print("[OK   ] Testigo sin mutar compila a /Werror y pasa todos los controles (0 fallas)")

        tareas = [
            (i, mid, aguja, reemplazo, original)
            for i, (mid, aguja, reemplazo) in enumerate(MUTANTES)
        ]
        max_w = min(len(tareas), os.cpu_count() or 4)
        resultados = [None] * len(tareas)
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_w) as ex:
            for idx, mid, estado, fallidas in ex.map(evaluar_mutante, tareas):
                resultados[idx] = (mid, estado, fallidas)

        cazados = 0
        for mid, estado, fallidas in resultados:
            if estado == "CAZADO":
                cazados += 1
                print(f"[CAZADO] {mid:<48s} -> cazado por {', '.join(fallidas)}")
            elif estado == "NO_PATRON":
                print(f"[FALLA] {mid}: no se encontro el patron en observador_lilyth.hpp")
                sys.exit(1)
            elif estado.startswith("NO_COMPILA"):
                print(f"[FALLA] {mid}: {estado}")
                sys.exit(1)
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