#!/usr/bin/env python3
# mutaciones_revision.py -- Code, 2026-11-07. Mutantes de las CORRECCIONES de Code sobre la entrega de Argon (D1, D2, D6, D7). Cada uno
# deshace una correccion en la cabecera (inc/) y la bateria completa tiene que fallar en un control que el testigo pasa. Un mutante que no
# compila o cuyo texto no esta en la cabecera NO esta cazado. Usa las funciones de correr_mutantes_argon.py.
import pathlib, sys
import correr_mutantes_argon as R

AQUI = pathlib.Path(__file__).resolve().parent
L = AQUI / "inc" / "ark" / "lilyth"

MUTANTES = [
    # (nombre, cabecera, etapa que lo debe cazar, viejo, nuevo)
    ("C01_via_a_acepta_menos_de_m_min", "observador_lilyth.hpp", "3",
     "if (!cfg_.forzar_via_b && esp.size() >= fn.m_min) {", "if (!cfg_.forzar_via_b && esp.size() >= 1) {"),
    ("C02_via_b_acepta_menos_de_m_min", "observador_lilyth.hpp", "3",
     "            if (taus.size() < fn.m_min) {\n                ref.sin_nulo_b = true;", "            if (taus.empty()) {\n                ref.sin_nulo_b = true;"),
    ("C03_pool_loo_sin_m_min", "nulo_lilyth.hpp", "3",
     "if (taus_i.size() < m_min) continue;", "if (taus_i.empty()) continue;"),
    ("C04_iid_por_innovaciones_lineales", "observador_lilyth.hpp", "3",
     "evaluar_una_referencia(fn, es_purgada, true, r, idx, ref_iid, u_pool_iid, emp_iid[idx]);",
     "evaluar_una_referencia(fn, es_purgada, fn.modelo_lineal.valido ? false : true, r, idx, ref_iid, u_pool_iid, emp_iid[idx]);"),
    ("C05_aperturas_posteriores_son_ajenas", "nulo_lilyth.hpp", "3",
     "const bool solo_contamina_lo_previo = true;", "const bool solo_contamina_lo_previo = false;"),
    ("C06_cortada_cuenta_como_admitida", "nulo_lilyth.hpp", "3",
     "    if (desenlace == 3u) return DestinoInyeccion::CORTADA;", "    if (desenlace == 3u) return DestinoInyeccion::ADMITIDA;"),
    ("C07_crudo_sin_despurgar", "nulo_lilyth.hpp", "3",
     "        if (hay_y && ctx.purga_activa && en_purgado) {", "        if (hay_y && ctx.purga_activa && en_purgado && false) {"),
    ("C08_espejo_sin_comprometer_el_dia", "nulo_lilyth.hpp", "3",
     "        if (dia != dia_ && dia_ >= 0) comprometido(s, n, db, vm);", "        (void)dia;"),
    ("C09_B_con_factor_2", "observador_lilyth.hpp", "3",
     "mm, cfg_delta.arl0_pasos_nativos, x_cal.size(),", "mm, 2.0 * cfg_delta.arl0_pasos_nativos, x_cal.size(),"),
    ("C10_excedencia_en_pool_sin_omitir", "nulo_lilyth.hpp", "3",
     "(ev.escala_origen, ev.signo, ev.g_disparo, pool, m_min, t_max_ef, i, n_a);",
     "(ev.escala_origen, ev.signo, ev.g_disparo, pool, m_min, t_max_ef, std::numeric_limits<std::size_t>::max(), n_a);"),
    ("C11_armonico_real_sobre_todas_las_ranuras", "observador_lilyth.hpp", "4",
     "arm.sobrepaso.min_eta_barra_real_emparejado_onsager = min_eta_barra_de(idx_ons);",
     "arm.sobrepaso.min_eta_barra_real_emparejado_onsager = arm.prony_curva_media.min_eta_desde_2 + 0.5;"),
    ("C12_armonico_sobrepaso_con_pool_iid", "observador_lilyth.hpp", "4",
     "arm.sobrepaso.min_eta_barra_real_emparejado_onsager, sets_ons, H, cfg_.b_mc_resumen,",
     "arm.sobrepaso.min_eta_barra_real_emparejado_onsager, sets_iid, H, cfg_.b_mc_resumen,"),
]

def main():
    originales = {}
    testigos = {}
    cazados = 0
    try:
        for nombre, cab, etapa, viejo, nuevo in MUTANTES:
            hdr = L / cab
            if cab not in originales:
                originales[cab] = hdr.read_text(encoding="utf-8")
            if etapa not in testigos:
                ok, err = R.compilar(etapa)
                if not ok:
                    print(f"[{etapa}] TESTIGO NO COMPILA: {err[:2]}"); return
                testigos[etapa] = R.correr(etapa)
                print(f"[{etapa}] testigo rc={testigos[etapa][0]} fallas={sorted(testigos[etapa][1])}", flush=True)
            orig = originales[cab]
            if orig.count(viejo) != 1:
                print(f"{nombre}: NO CAZADO (el texto no esta, o no es unico: {orig.count(viejo)})", flush=True); continue
            hdr.write_text(orig.replace(viejo, nuevo, 1), encoding="utf-8")
            ok, err = R.compilar(etapa)
            if not ok:
                print(f"{nombre}: NO CAZADO (no compila: {err[:1]})", flush=True)
            else:
                rc, f = R.correr(etapa)
                rc0, f0 = testigos[etapa]
                nuevas = f - f0
                cazado = bool(nuevas) or (rc != 0 and rc0 == 0)
                cazados += cazado
                print(f"{nombre}: {'CAZADO' if cazado else 'NO CAZADO'} rc={rc} nuevas={sorted(nuevas)}", flush=True)
            hdr.write_text(orig, encoding="utf-8")
    finally:
        for cab, txt in originales.items():
            (L / cab).write_text(txt, encoding="utf-8")
    print(f"{cazados}/{len(MUTANTES)} cazados", flush=True)

if __name__ == "__main__":
    main()
