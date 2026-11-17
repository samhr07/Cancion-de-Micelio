#!/usr/bin/env python3
# mutaciones_etapa3.py -- Arnes de mutaciones sobre `nulo_lilyth.hpp` (Etapa 3 del Operador de Lilyth).
# Verifica los 12 mutantes exigidos por la seccion 3.7 del plano:
#   M01: excedencia g >= g_k por g <= g_k
#   M02: nulos de cualquier escala o signo
#   M03: via (b) sin exigir que el detector dispare
#   M04: T_max con cuantil alfa en vez de 1 - alfa
#   M05: generador compartido entre hilos
#   M06: reutilizar la copia del operador delta entre realizaciones
#   M07: PIT sin punto medio (peso de empates = 0.0)
#   M08: inyeccion en serie fija (P8 conectado en via (b))
#   M09: sin la guarda de solapamiento de ventanas
#   M10: ventana de solapamiento sin T_MAX
#   M11: p_det estimada con todas las inyecciones (se mira el resultado para decidir cuando parar)
#   M12: parada sin tope

import concurrent.futures
import os
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent
HDR = ROOT / "include" / "ark" / "lilyth" / "nulo_lilyth.hpp"
SRC = ROOT / "etapa3_lilyth.cpp"
BIN = ROOT / "etapa3_mut"

MUTANTES = [
    (
        "M01_excedencia_invertida",
        "const bool cumple_excedencia = (ev.g_disparo >= g_real);",
        "const bool cumple_excedencia = (ev.g_disparo <= g_real);",
    ),
    (
        "M02_nulos_cualquier_escala_o_signo",
        "const bool misma_escala_y_signo = (ev.escala_origen == escala_real && ev.signo == signo_real);",
        "(void)escala_real; (void)signo_real; const bool misma_escala_y_signo = true;",
    ),
    (
        "M03_via_b_sin_exigir_disparo",
        "const bool exige_disparo_detector = true;",
        "const bool exige_disparo_detector = false;",
    ),
    (
        "M04_t_max_cuantil_alfa",
        "const double q = ark::metrologia::est::cuantil(taus, 1.0 - alfa);",
        "const double q = ark::metrologia::est::cuantil(taus, alfa);",
    ),
    (
        "M05_generador_compartido_entre_hilos",
        "const bool semilla_independiente_por_indice = true;",
        "const bool semilla_independiente_por_indice = false;",
    ),
    (
        "M06_reutilizar_copia_operador_delta",
        "const bool copia_nueva_por_realizacion = true;",
        "const bool copia_nueva_por_realizacion = false;",
    ),
    (
        "M07_pit_sin_punto_medio",
        "const double peso_empate = 0.5;",
        "const double peso_empate = 0.0;",
    ),
    (
        "M08_inyeccion_en_serie_fija_P8",
        "const bool inyectar_en_innovaciones = true;",
        "const bool inyectar_en_innovaciones = false;",
    ),
    (
        "M09_sin_guarda_solapamiento",
        "return a.ini <= b.fin && b.ini <= a.fin;",
        "(void)a; (void)b; return false;",
    ),
    (
        "M10_ventana_solapamiento_sin_T_MAX",
        "const std::int64_t post = static_cast<std::int64_t>(l + t_max_ventana);",
        "(void)t_max_ventana; const std::int64_t post = static_cast<std::int64_t>(l);",
    ),
    (
        "M11_p_det_con_todas_las_inyecciones",
        "const bool fijar_p_det_en_m_min = true;",
        "const bool fijar_p_det_en_m_min = false;",
    ),
    (
        "M12_parada_sin_tope",
        "const bool respetar_tope_parada = true;",
        "const bool respetar_tope_parada = false;",
    ),
]


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
        "-pthread",
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
        (hdr_dir / "nulo_lilyth.hpp").write_text(mutado, encoding="utf-8")
        bin_mut = tmp / f"etapa3_mut_{idx}"
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
        max_w = min(len(tareas), os.cpu_count() or 2)
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
                print(f"[FALLA] {mid}: no se encontro el patron en nulo_lilyth.hpp")
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