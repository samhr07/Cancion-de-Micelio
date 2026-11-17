#!/usr/bin/env python3
# correr_mutantes_argon.py -- Code (revision de las Etapas 2-6 de Lilyth, 2026-11-07).
# Corre los MUTANTES de los scripts de Argon (argon_tal_cual/mutaciones_etapaN.py, que llaman a g++ con su propio arbol) con el
# compilador del proyecto (cl /W4) sobre la copia de revision `inc/`, y SIEMPRE con la bateria completa (sin --rapido).
# Cazado = la bateria con el mutante falla en un control que el testigo pasa, o se cae. Un mutante que no compila o cuyo texto no
# esta en la cabecera NO esta cazado.
import importlib.util, pathlib, re, subprocess, sys

AQUI = pathlib.Path(__file__).resolve().parent
ETAPAS = {"2": "modelo_lineal_nulo", "3_0": "observador_lilyth", "3": "nulo_lilyth",
          "4": "armonico_lilyth", "5": "fin_calor_lilyth", "6": "reloj_marea_lilyth"}

def mutantes(n):
    s = importlib.util.spec_from_file_location("m" + n, AQUI / "argon_tal_cual" / f"mutaciones_etapa{n}.py")
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
    return m.MUTANTES

def compilar(n):
    r = subprocess.run(["cmd", "/c", str(AQUI / "compilar.bat"), f"etapa{n}_lilyth.cpp", f"mut_{n}.exe"], cwd=AQUI, capture_output=True,
                       text=True, encoding="mbcs", errors="replace")
    salida = (r.stdout + r.stderr)
    errores = [l for l in salida.splitlines() if " error " in l or "warning C" in l]
    return r.returncode == 0 and not errores, errores

def correr(n):
    try:
        r = subprocess.run([str(AQUI / "build" / f"mut_{n}.exe")], cwd=AQUI, capture_output=True, text=True,
                           encoding="mbcs", errors="replace", timeout=900)
    except subprocess.TimeoutExpired:
        return -999, {"TIMEOUT"}
    txt = r.stdout + r.stderr
    fallas = set(re.findall(r"\[FALLA\]\s+(\S+)", txt))
    fallas |= {l.strip()[:90] for l in txt.splitlines() if "FALLA:" in l and "[FALLA]" not in l and "RESULTADO" not in l}
    return r.returncode, fallas

def main(etapas):
    for n in etapas:
        hdr = AQUI / "inc" / "ark" / "lilyth" / (ETAPAS[n] + ".hpp")
        original = hdr.read_text(encoding="utf-8")
        ok, err = compilar(n)
        if not ok:
            print(f"[{n}] TESTIGO NO COMPILA: {err[:3]}"); continue
        rc0, f0 = correr(n)
        print(f"[{n}] testigo rc={rc0} fallas={sorted(f0)}", flush=True)
        cazados = 0; lista = mutantes(n)
        try:
            for mu in lista:
                nombre, viejo, nuevo = mu[0], mu[1], mu[2]
                if viejo not in original:
                    print(f"[{n}] {nombre}: NO CAZADO (el texto no esta en la cabecera)", flush=True); continue
                hdr.write_text(original.replace(viejo, nuevo, 1), encoding="utf-8")
                ok, err = compilar(n)
                if not ok:
                    print(f"[{n}] {nombre}: NO CAZADO (no compila: {err[:1]})", flush=True); continue
                rc, f = correr(n)
                nuevas = f - f0
                cazado = (rc != 0 and rc0 == 0) or bool(nuevas) or rc not in (0, 1)
                cazados += cazado
                print(f"[{n}] {nombre}: {'CAZADO' if cazado else 'NO CAZADO'} rc={rc} nuevas={sorted(nuevas)}", flush=True)
        finally:
            hdr.write_text(original, encoding="utf-8")
        print(f"[{n}] {cazados}/{len(lista)} cazados", flush=True)
    ok, _ = compilar(etapas[-1]) if etapas else (True, [])

if __name__ == "__main__":
    main(sys.argv[1:] or list(ETAPAS))
