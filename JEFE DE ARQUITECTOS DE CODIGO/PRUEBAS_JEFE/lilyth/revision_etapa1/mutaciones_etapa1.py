# -*- coding: utf-8 -*-
"""Mutantes de la Etapa 1 del Operador de Lilyth (Code, 2026-11-05).

Cada mutante cambia UNA cosa de ark/lilyth/observador_lilyth.hpp en una copia de la base (build/mut/<id>/include), compila
etapa1_lilyth.cpp contra ella y exige que FALLE (rc != 0). Un mutante que no compila NO esta cazado. El testigo tiene que PASAR.
Uso: python -I mutaciones_etapa1.py [hilos]
"""
import concurrent.futures as cf
import os
import shutil
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(AQUI, "..", "..", "..", "..", "ark_cxx", "base_auditada", "include"))
H = "ark/lilyth/observador_lilyth.hpp"
MUTANTES = [
    ("M0", "testigo: sin cambio", None, None),
    ("M1", "abre con cualquier alarma (sin espejo del evento abierto)", "        if (!espejo_abierto_) {\n            for (std::size_t j = 0;",
     "        if (true) {\n            for (std::size_t j = 0;"),
    ("M2", "guarda del ancla 0", "        r.guarda_ancla = static_cast<std::size_t>(std::ceil(ti)) * m;", "        r.guarda_ancla = 0 * m;"),
    ("M3", "el ancla llega hasta t_cambio", "        r.t_fin_ancla = r.t_cambio - static_cast<std::int64_t>(r.guarda_ancla) - 1;", "        r.t_fin_ancla = r.t_cambio;"),
    ("M4", "W_ancla con la memoria del GLR (det.L) en vez de N_REG", "        r.w_ancla = Rasgos::n_reg * m;", "        r.w_ancla = f.det.L * m;"),
    ("M5", "dY a la mitad", "        const double dy = y_a - med;", "        const double dy = 0.5 * (y_a - med);"),
    ("M6", "sin la guarda del signo de dY", "        if (!(std::isfinite(dy) && dy * static_cast<double>(signo) > 0.0)) {", "        if (!std::isfinite(dy)) {"),
    ("M7", "el consolidado ignora el calor", "            if (!calor_activo_ && !std::isfinite(rama.tau_rest_consolidado)) {",
     "            if (!std::isfinite(rama.tau_rest_consolidado)) {"),
    ("M8", "la replica no cierra las ramas ya consolidadas (Argon)",
     "            if (std::isfinite(rama->tau_rest_consolidado)) {\n                rama->desenlace = DesenlaceRanura::RESTITUIDA;",
     "            if (std::isfinite(rama->tau_rest_consolidado)) {\n                continue;\n                rama->desenlace = DesenlaceRanura::RESTITUIDA;"),
    ("M9", "espejo de la purga sin la historia de calibracion", "        if (purga_activa_) {\n            bool hp = false;",
     "        if (false) {\n            bool hp = false;"),
    ("M10", "el corte de la reintegracion no abre segmento", "                    if (purgado_vivo_) ++seg_purgado_;", "                    (void)0;"),
    ("M11", "origen del nivel purgado en x(t) en vez de x(t - 1)", "nivel_purgado_ = x - y_rep;", "nivel_purgado_ = x;"),
    ("M12", "Aalen-Johansen sin S(tau - 1)", "                c_rest += s_prev * (static_cast<double>(dr) / denom);",
     "                c_rest += (static_cast<double>(dr) / denom);"),
    ("M13", "tau_fin de la replica uno de mas", "                rama->tau_fin = tau;\n            }\n        }\n        if (corto)",
     "                rama->tau_fin = tau + 1;\n            }\n        }\n        if (corto)"),
    ("M14", "sin la guarda del reloj del operador delta", "        if (s.t != t_ + 1) {", "        if (false) {"),
    ("M15", "horizonte uno de mas", "        if (tau >= cfg_.t_max_efectivo) {", "        if (tau > cfg_.t_max_efectivo) {"),
    ("M16", "sin la guarda de memoria del ancla", "        if (!indice(t_ini, i0)) return DesenlaceRanura::NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA;", "        (void)i0;"),
    ("M17", "PURGADA y CRUDA invertidas", "        return cfg_.representacion != ModoRepresentacionLilyth::CRUDA;",
     "        return cfg_.representacion == ModoRepresentacionLilyth::CRUDA;"),
    ("M18", "el espejo se cierra con eventos de temperatura", "            if (e.canal != 0) continue;\n            ++cont_.cierres_espejo;",
     "            ++cont_.cierres_espejo;"),
]


def correr(m):
    mid, desc, viejo, nuevo = m
    d = os.path.join(AQUI, "build", "mut", mid)
    shutil.rmtree(d, ignore_errors=True)
    inc = os.path.join(d, "include")
    shutil.copytree(BASE, inc)
    if viejo is not None:
        ruta = os.path.join(inc, H)
        s = open(ruta, encoding="utf-8").read()
        n = s.count(viejo)
        if n != 1:
            return mid, desc, "MAL DEFINIDO (%d apariciones)" % n, False
        open(ruta, "w", encoding="utf-8", newline="\n").write(s.replace(viejo, nuevo))
    exe = "mut_%s.exe" % mid
    r = subprocess.run(["cmd", "/c", r".\compilar.bat", "etapa1_lilyth.cpp", exe, inc], capture_output=True, text=True, cwd=AQUI)
    ruta_exe = os.path.join(AQUI, "build", exe)
    if not os.path.exists(ruta_exe):
        return mid, desc, "NO COMPILA (no esta cazado)", False
    p = subprocess.run([ruta_exe], capture_output=True, text=True, cwd=AQUI)
    fallas = [l.split()[2] for l in p.stdout.splitlines() if l.startswith("[FALLA]")]
    if mid == "M0":
        return mid, desc, "testigo PASA" if p.returncode == 0 else "TESTIGO FALLA: " + ",".join(fallas), p.returncode == 0
    return mid, desc, ("CAZADO por " + ",".join(fallas)) if p.returncode != 0 else "NO CAZADO", p.returncode != 0


def main():
    hilos = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    # cada compilacion escribe su .obj en build\: en paralelo se pisarian; se compila en serie y se corre en paralelo
    res = [correr(m) for m in MUTANTES] if hilos <= 1 else list(cf.ThreadPoolExecutor(1).map(correr, MUTANTES))
    ok = all(r[3] for r in res)
    for mid, desc, txt, bien in res:
        print("[%s] %-4s %s: %s" % ("OK " if bien else "MAL", mid, desc, txt))
    cazados = sum(1 for r in res[1:] if r[3])
    print("\n%d/%d mutantes cazados; testigo %s" % (cazados, len(res) - 1, "PASA" if res[0][3] else "FALLA"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
