#!/usr/bin/env python3
"""
mutaciones_nivel4.py -- Campana de mutacion de regresion sellada para Nivel 4.

Opera exclusivamente sobre copias temporales aisladas en %TEMP%/ark_mut_nivel4_sellado,
sin modificar los encabezados reales del piloto. Verifica:
  MN4.0  : testigo sin mutar sobre test_nv4_admision_propia.exe completo (debe PASAR, 0 fallas).
  M_N4_1 : reintroduce filtro duro escala_origen == escala_real en ejecutar_inyecciones_en_realizacion
           cuando usar_enmienda_nivel3 esta activa (debe FALLAR --solo-operador -> CAZADO).
  M_N4_2 : reintroduce exigencia de excedencia g_disparo >= g_real propagada desde completar_nulos()
           cuando usar_enmienda_nivel3 esta activa (debe FALLAR --solo-operador -> CAZADO).
"""

from __future__ import annotations

import argparse
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile


ROOT_PILOT = pathlib.Path(__file__).resolve().parents[1]
RUTA_TEST_SRC = ROOT_PILOT / "nivel4_admision_propia" / "test_nv4_admision_propia.cpp"
RUTA_NULO = ROOT_PILOT / "inc" / "ark" / "lilyth" / "nulo_lilyth.hpp"
RUTA_OBS = ROOT_PILOT / "inc" / "ark" / "lilyth" / "observador_lilyth.hpp"
RUTA_COMPILAR = ROOT_PILOT / "compilar.bat"


def reemplazar_exacto(texto: str, objetivo: str, reemplazo: str, esperado: int, etiqueta: str) -> str:
    cuenta = texto.count(objetivo)
    if cuenta != esperado:
        raise RuntimeError(
            f"[{etiqueta}] coincidencia exacta esperada={esperado}, observada={cuenta} para:\n{objetivo}"
        )
    return texto.replace(objetivo, reemplazo)


def aplicar_mutacion(tag: str, inc_tmp: pathlib.Path) -> None:
    if tag == "MN4.0":
        return

    nulo_path = inc_tmp / "ark" / "lilyth" / "nulo_lilyth.hpp"
    obs_path = inc_tmp / "ark" / "lilyth" / "observador_lilyth.hpp"

    nulo_txt = nulo_path.read_text(encoding="utf-8")
    obs_txt = obs_path.read_text(encoding="utf-8")

    if tag == "M_N4_1":
        nulo_txt = reemplazar_exacto(
            nulo_txt,
            "            if (dest == DestinoInyeccion::ADMITIDA && propia->extraible) {",
            "            if (dest == DestinoInyeccion::ADMITIDA && propia->extraible && (!usar_enmienda_nivel3 || escala_real < 0 || propia->ev.escala_origen == escala_real)) {",
            1,
            "M_N4_1 filtro escala",
        )
        nulo_path.write_text(nulo_txt, encoding="utf-8", newline="\n")
        return

    if tag == "M_N4_2":
        nulo_txt = reemplazar_exacto(
            nulo_txt,
            "    GanchoHijas gancho_hijas = SinHijasPaso{}) {",
            "    GanchoHijas gancho_hijas = SinHijasPaso{}, double g_real_mut = NaN) {",
            2,
            "M_N4_2 firmas ejecutar_inyecciones y ejecutar_via_b_evento",
        )
        nulo_txt = reemplazar_exacto(
            nulo_txt,
            "            perfil_est, e_iny, signo_real, w_ancla, guarda, t_max_ef, t_arr, n_sim, rng, true, usar_enmienda_nivel3, escala_real, gancho_hijas);",
            "            perfil_est, e_iny, signo_real, w_ancla, guarda, t_max_ef, t_arr, n_sim, rng, true, usar_enmienda_nivel3, escala_real, gancho_hijas, g_real_mut);",
            1,
            "M_N4_2 reenvio g_real_mut",
        )
        nulo_txt = reemplazar_exacto(
            nulo_txt,
            "            if (dest == DestinoInyeccion::ADMITIDA && propia->extraible) {",
            "            if (dest == DestinoInyeccion::ADMITIDA && propia->extraible && (!usar_enmienda_nivel3 || !std::isfinite(g_real_mut) || (std::isfinite(propia->ev.g_disparo) && propia->ev.g_disparo >= g_real_mut))) {",
            1,
            "M_N4_2 excedencia g_disparo >= g_real",
        )
        obs_txt = reemplazar_exacto(
            obs_txt,
            "                cfg_.semilla_nulo ^ (sal * (idx + 1)), ref, ev_b, cont_.solapes, einy_n2, cfg_.usar_enmienda_nivel3, r.escala_origen);",
            "                cfg_.semilla_nulo ^ (sal * (idx + 1)), ref, ev_b, cont_.solapes, einy_n2, cfg_.usar_enmienda_nivel3, r.escala_origen, detalle_nulo::SinHijasPaso{}, r.g_disparo);",
            1,
            "M_N4_2 llamada completar_nulos",
        )
        nulo_path.write_text(nulo_txt, encoding="utf-8", newline="\n")
        obs_path.write_text(obs_txt, encoding="utf-8", newline="\n")
        return

    raise ValueError(f"Mutante desconocido: {tag}")


def compilar_en_temp(tmp_root: pathlib.Path, tag: str, src_rel: str, exe_name: str, inc_tmp: pathlib.Path) -> pathlib.Path:
    cmd = [
        "cmd.exe",
        "/d",
        "/c",
        str(tmp_root / "compilar.bat"),
        src_rel,
        exe_name,
        str(inc_tmp),
    ]
    proc = subprocess.run(cmd, cwd=str(tmp_root), capture_output=True, text=True)
    sal = (proc.stdout or "") + (proc.stderr or "")
    avisos = [ln for ln in sal.splitlines() if "warning C" in ln]
    if proc.returncode != 0 or avisos:
        raise RuntimeError(
            f"[{tag}] error de compilacion (rc={proc.returncode}, avisos_w4={len(avisos)}):\n{sal}"
        )
    exe_path = tmp_root / "build" / exe_name
    if not exe_path.is_file():
        raise RuntimeError(f"[{tag}] binario esperado no generado: {exe_path}")
    return exe_path


def ejecutar_campana(solo_compilar: bool, hilos: int) -> int:
    for p in [RUTA_TEST_SRC, RUTA_NULO, RUTA_OBS, RUTA_COMPILAR]:
        if not p.is_file():
            raise FileNotFoundError(f"Archivo requerido no encontrado: {p}")

    work_root = pathlib.Path(tempfile.gettempdir()) / "ark_mut_nivel4_sellado"
    if work_root.exists():
        shutil.rmtree(work_root)
    work_root.mkdir(parents=True, exist_ok=True)

    casos = [
        ("MN4.0", "testigo intacto (P0_NIVEL3_BASE)", False, False),
        ("M_N4_1", "filtro duro escala_origen == escala_real", True, True),
        ("M_N4_2", "exigencia de excedencia g_disparo >= g_real", True, True),
    ]

    try:
        shutil.copy2(RUTA_COMPILAR, work_root / "compilar.bat")
        (work_root / "nivel4_admision_propia").mkdir(parents=True, exist_ok=True)
        shutil.copy2(RUTA_TEST_SRC, work_root / "nivel4_admision_propia" / "test_nv4_admision_propia.cpp")

        print("================================================================================")
        print(f" CAMPANA DE MUTACION NIVEL 4 (hilos={hilos}, solo_compilar={int(solo_compilar)})")
        print(f"  Directorio aislado: {work_root}")
        print("================================================================================")

        cazados = 0
        mutantes_tot = 0
        testigo_ok = False

        for tag, desc, es_mutante, solo_op in casos:
            inc_tmp = work_root / f"inc_{tag}"
            shutil.copytree(ROOT_PILOT / "inc", inc_tmp)
            aplicar_mutacion(tag, inc_tmp)
            exe_name = f"nv4_{tag}.exe"
            exe_path = compilar_en_temp(
                work_root,
                tag,
                "nivel4_admision_propia\\test_nv4_admision_propia.cpp",
                exe_name,
                inc_tmp,
            )
            print(f"[{tag}] compilacion limpia (/W4 0 avisos): {desc}")
            if solo_compilar:
                continue

            args = [str(exe_path), f"--hilos={hilos}"]
            if solo_op:
                args.append("--solo-operador")
            run = subprocess.run(args, cwd=str(work_root), capture_output=True, text=True)
            sal = run.stdout or ""
            fallas_lineas = [ln.strip() for ln in sal.splitlines() if "FALLA" in ln]

            if not es_mutante:
                pasa_mn40 = (run.returncode == 0) and ("FALLAS NV4: 0" in sal)
                testigo_ok = pasa_mn40
                estado = "TESTIGO PASA (0 fallas)" if pasa_mn40 else f"TESTIGO FALLA (rc={run.returncode})"
                print(f"[{tag}] {estado}")
                for ln in fallas_lineas:
                    print(f"        {ln}")
                if not pasa_mn40:
                    print(sal)
            else:
                mutantes_tot += 1
                cazado = (run.returncode != 0)
                if cazado:
                    cazados += 1
                estado = f"CAZADO (rc={run.returncode})" if cazado else "VIVO (NO CAZADO)"
                print(f"[{tag}] {estado}")
                for ln in fallas_lineas[:12]:
                    print(f"        {ln}")

        print("--------------------------------------------------------------------------------")
        if solo_compilar:
            print("MODO --solo-compilar COMPLETADO: 3/3 compilaciones limpias (/W4 0 avisos)")
            return 0

        print(f"RESUMEN MUTACIONES NIVEL 4: MN4.0={'PASA' if testigo_ok else 'FALLA'} | {cazados} de {mutantes_tot} mutantes cazados")
        ok = testigo_ok and (cazados == mutantes_tot)
        print(f"FALLAS MUTACIONES NV4: {0 if ok else 1}")
        return 0 if ok else 1
    finally:
        if work_root.exists():
            shutil.rmtree(work_root, ignore_errors=True)


def main() -> int:
    parser = argparse.ArgumentParser(description="Campana de mutacion Nivel 4 en copias temporales.")
    parser.add_argument("--solo-compilar", action="store_true", help="Verifica anclas y compilacion /W4 sin ejecutar")
    parser.add_argument("--hilos", type=int, default=max(1, os.cpu_count() or 4), help="Hilos para test_nv4_admision_propia")
    args = parser.parse_args()
    return ejecutar_campana(args.solo_compilar, max(1, args.hilos))


if __name__ == "__main__":
    sys.exit(main())
