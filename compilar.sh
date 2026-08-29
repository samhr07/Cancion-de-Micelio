#!/usr/bin/env bash
# Compila micelio_hallazgos.tex -> micelio_hallazgos.pdf
# Dos pasadas: la primera genera el .toc, la segunda lo resuelve.
set -e
TEX="micelio_hallazgos"

command -v pdflatex >/dev/null || {
  echo "Falta pdflatex. En Debian/Ubuntu:"
  echo "  sudo apt install texlive-latex-recommended texlive-latex-extra texlive-lang-spanish"
  exit 1; }

echo "[1/2] primera pasada..."
pdflatex -interaction=nonstopmode -halt-on-error "$TEX.tex" > /dev/null
echo "[2/2] segunda pasada (indice)..."
pdflatex -interaction=nonstopmode -halt-on-error "$TEX.tex" > /dev/null

rm -f "$TEX".{aux,log,out,toc}
echo "OK -> $TEX.pdf  ($(du -h "$TEX.pdf" | cut -f1))"

# Verificacion visual opcional: primera pagina a PNG
if command -v pdftoppm >/dev/null; then
  pdftoppm -png -r 100 -f 1 -l 1 "$TEX.pdf" verificacion
  echo "Verificacion -> verificacion-1.png"
fi
