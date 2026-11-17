# -*- coding: utf-8 -*-
"""verificar_base.py -- manifiesto sha256 de la base auditada y verificacion de copias.

  python -I verificar_base.py --generar            escribe MANIFIESTO_SHA256.txt de include/ (solo la fuente: ark_cxx/base_auditada)
  python -I verificar_base.py --verificar DIR      compara DIR (una copia: DIR/ark/...) contra el manifiesto; rc 0 si es identica

Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04. Una copia que no verifica NO es la base auditada: se recopia o se archiva.
"""
import hashlib
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
INCLUDE = os.path.join(AQUI, "include")
MANIFIESTO = os.path.join(AQUI, "MANIFIESTO_SHA256.txt")


def huellas(raiz):
    out = {}
    for d, _, fs in os.walk(raiz):
        for f in fs:
            if f.endswith((".hpp", ".md")):
                p = os.path.join(d, f)
                rel = os.path.relpath(p, raiz).replace(os.sep, "/")
                if not rel.startswith("ark/"):   # solo el arbol ark/ (el README propio de una copia no cuenta)
                    continue
                with open(p, "rb") as fh:
                    out[rel] = hashlib.sha256(fh.read()).hexdigest()
    return out


def leer_manifiesto():
    m = {}
    with open(MANIFIESTO, encoding="utf-8") as fh:
        for linea in fh:
            if linea.strip() and not linea.startswith("#"):
                h, rel = linea.strip().split("  ", 1)
                m[rel] = h
    return m


def main():
    if len(sys.argv) >= 2 and sys.argv[1] == "--generar":
        h = huellas(INCLUDE)
        with open(MANIFIESTO, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# sha256  ruta (relativa a include/) -- fuente: ark_cxx/base_auditada/include\n")
            for rel in sorted(h):
                fh.write(f"{h[rel]}  {rel}\n")
        print(f"manifiesto: {len(h)} archivos")
        return 0
    if len(sys.argv) >= 3 and sys.argv[1] == "--verificar":
        m, c = leer_manifiesto(), huellas(sys.argv[2])
        falta = sorted(set(m) - set(c))
        sobra = sorted(set(c) - set(m))
        distinto = sorted(r for r in set(m) & set(c) if m[r] != c[r])
        for r in falta:
            print("FALTA    ", r)
        for r in sobra:
            print("SOBRA    ", r)
        for r in distinto:
            print("DISTINTO ", r)
        ok = not (falta or sobra or distinto)
        print(f"copia {'IDENTICA' if ok else 'NO VERIFICA'}: {len(c)} archivos, {len(falta)} faltan, {len(sobra)} sobran, {len(distinto)} distintos")
        return 0 if ok else 1
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main())
