#!/usr/bin/env python3
# correr_baterias.py -- compila y corre las baterias de Lilyth contra inc/ de esta carpeta (o el que se pase), en paralelo.
# Uso: python correr_baterias.py [nombre ...]   (sin nombres: las de las Etapas 2-6 y la de L+H si existe)
import concurrent.futures as cf
import pathlib
import subprocess
import sys
import time

AQUI = pathlib.Path(__file__).resolve().parent
TODAS = ["etapa2_lilyth", "etapa3_0_lilyth", "etapa3_lilyth", "etapa4_lilyth", "etapa5_lilyth", "etapa6_lilyth", "etapa_lh_lilyth"]


def compilar(nombre):
    r = subprocess.run(["cmd", "/c", str(AQUI / "compilar.bat"), f"{nombre}.cpp", f"{nombre}.exe"], cwd=AQUI,
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    salida = (r.stdout + r.stderr).strip().splitlines()
    avisos = [l for l in salida if "warning" in l or "error" in l]
    return r.returncode == 0, avisos


def correr(nombre):
    t0 = time.time()
    r = subprocess.run([str(AQUI / "build" / f"{nombre}.exe")], cwd=AQUI, capture_output=True, text=True, encoding="utf-8", errors="replace")
    (AQUI / f"{nombre}.log").write_text(r.stdout + r.stderr, encoding="utf-8")
    fallas = [l for l in r.stdout.splitlines() if "FALLA" in l or "FAIL" in l]
    return r.returncode, time.time() - t0, fallas


def una(nombre):
    ok, avisos = compilar(nombre)
    if not ok or avisos:
        return nombre, "NO COMPILA" if not ok else "AVISOS", avisos[:6], 0.0
    rc, seg, fallas = correr(nombre)
    return nombre, "PASA" if rc == 0 else f"rc={rc}", fallas[:12], seg


def main():
    nombres = sys.argv[1:] or [n for n in TODAS if (AQUI / f"{n}.cpp").exists()]
    with cf.ThreadPoolExecutor(max_workers=len(nombres)) as ex:
        for nombre, estado, det, seg in ex.map(una, nombres):
            print(f"{nombre}: {estado} ({seg:.0f} s)", flush=True)
            for d in det:
                print(f"    {d}", flush=True)


if __name__ == "__main__":
    main()
