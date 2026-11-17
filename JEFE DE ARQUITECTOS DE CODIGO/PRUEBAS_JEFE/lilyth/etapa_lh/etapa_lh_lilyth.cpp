// etapa_lh_lilyth.cpp -- controles del nulo lineal con replicas (L+H, E1.24) CONTRA EL OPERADOR DELTA REAL. C++20.
// Code (Jefe de Arquitectos de Codigo, Opus 5.5), 2026-11-11, por la directriz MENSAJE_ARGON_NULO_CON_REPLICAS.md.
// Los controles de unidad del ajuste (U1-U4, R1u, R4, R5, R2bu) estan en replicas_unidad.cpp; la familia de 6 y R9 en etapa4_lilyth.cpp.
//
// Opcion (b) del Comandante (2026-11-11): el nucleo es ESTRUCTURAL (hijas que disparan las alarmas de la simulacion, por inferencia
// indirecta), no el de los choques: R2 y R2b ya no comparan con la n verdadera de los choques (se reporta al lado) sino con lo que el
// estructural tiene que cumplir. El operador delta se calibra sobre una historia LIMPIA del mismo mundo y Lilyth sobre la que tiene los
// choques: el operador delta calibrado sobre choques frecuentes deja de verlos (segundo hallazgo, auditado aparte).
//   R1   choques SIN replicas: las alarmas no se agrupan mas que en el nulo de Onsager con su detector (las realarmas del detector no
//        son replicas): NO DECIDIBLE en una fraccion >= 1 - alfa de las semillas (banda binomial); L+H entra en Holm con p = 1
//   R2   replicas del MISMO signo (n = 0.3 y 0.6): pasa su nulo, la inferencia indirecta converge (la simulacion reproduce el Hawkes de
//        las alarmas reales), la cruzada estructural compatible con 0, y L+H con tamano <= alfa sobre las semillas; Onsager en n = 0.6
//        se reporta (el hallazgo de S1 sobre semillas)
//   R2b  replicas de signo CONTRARIO: la cruzada estructural excluye el 0 (la detecta), y la fraccion de inyecciones cortadas por una
//        replica en el nulo L+H es la de los eventos reales dentro de su banda
//   R3   TAR no lineal con replicas: L+H rechaza (la no linealidad no se la come la autoexcitacion)
//   R6   causalidad: la ficha de replicas y la biblioteca de marcas son las de calibrar despues de completar_nulos, y no cambian con
//        el tramo medido (caza M1 y M6)
//   R7   solape: las hijas son parte de su inyeccion; solapes = aperturas previas, y las reaperturas posteriores (las hijas) no
//        contaminan
//   R8   identidad: con a = 0 (generador de hijas presente pero sin excitacion) la via (b) L+H es la de Onsager bit a bit

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <random>
#include <thread>
#include <vector>

#include "ark/lilyth/observador_lilyth.hpp"

using namespace ark::metrologia;
using namespace ark::lilyth;

namespace {

using ObsH = ObservadorSaltos<24, 4, 120>;
constexpr std::size_t TMAXC = 128;
using LilH = ObservadorLilyth<ObsH, 8, TMAXC, 2048>;

int fallos = 0;
void chk(bool ok, const char* id, const char* txt) {
    std::printf("[%s] %-4s %s\n", ok ? "OK   " : "FALLA", id, txt);
    if (!ok) ++fallos;
}
bool igual(double a, double b) { return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof a) == 0; }
double banda(std::size_t n, double alfa) {
    return alfa * double(n) + est::norm_ppf(1.0 - alfa) * std::sqrt(double(n) * alfa * (1.0 - alfa));
}

// Mundo: AR(1) (o TAR en el tramo medido) con choques que son un Hawkes de dos tipos en el tiempo de los choques: madres Bernoulli
// por paso; cada choque tiene Poisson(n_mismo) hijas de su signo y Poisson(n_contrario) del signo contrario, con retardo
// d0 + Geom(1 - exp(-1 / tau_ret)) (d0 >= el tiempo muerto del espejo: ninguna hija se pierde). La verdad del proceso de ALARMAS es
// la del proceso de choques mientras el detector dispare en cada choque.
struct Parametros {
    double phi = 0.88;
    bool tar = false;             // tramo medido TAR: phi 0.96 si |x| > umbral_tar, 0.35 si no (la calibracion es AR(phi_cal_tar))
    double umbral_tar = 2.2;
    double sd_ruido_tar = 0.9;
    double tasa_madre = 1.0 / 300.0;
    double n_mismo = 0.0, n_contrario = 0.0;
    std::size_t d0 = 6;
    double tau_ret = 4.0;
    // Amplitud del choque en sd de la innovacion. UN paso (a2 = 0): un choque de dos pasos es una desviacion legitima del nulo de
    // Onsager iid (dos innovaciones grandes seguidas) y el ajuste la ve como agrupamiento sin que haya replicas (sonda_ii.log: 2 pasos
    // sin replicas DECIDIBLE con rho 0.225). R3 (TAR) sigue con dos pasos a proposito: ahi se busca que L+H rechace.
    double a1 = 9.5, a2 = 0.0;
    std::size_t n_cal = 40000, n_test = 8000;
    double arl0 = 800.0;
};

std::vector<double> generar(const Parametros& P, std::size_t n, bool tramo_medido, std::uint64_t sem, std::size_t* choques = nullptr) {
    std::mt19937_64 rng(sem);
    const bool usar_tar = P.tar && tramo_medido;
    const double sd_e = P.tar ? P.sd_ruido_tar : std::sqrt(1.0 - P.phi * P.phi);
    std::normal_distribution<double> nd(0.0, sd_e);
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::vector<double> eps(n);
    for (double& e : eps) e = nd(rng);
    std::vector<std::pair<std::size_t, int>> cola;
    for (std::size_t t = 200; t + 200 < n; ++t)
        if (u(rng) < P.tasa_madre) cola.emplace_back(t, u(rng) < 0.5 ? +1 : -1);
    std::geometric_distribution<int> geo(1.0 - std::exp(-1.0 / P.tau_ret));
    std::poisson_distribution<int> pm(std::max(P.n_mismo, 1e-300)), pc(std::max(P.n_contrario, 1e-300));
    std::size_t cab = 0;
    while (cab < cola.size() && cola.size() < n) {
        const auto [t, s] = cola[cab++];
        const int km = P.n_mismo > 0.0 ? pm(rng) : 0, kc = P.n_contrario > 0.0 ? pc(rng) : 0;
        for (int k = 0; k < km + kc; ++k) {
            const std::size_t th = t + P.d0 + static_cast<std::size_t>(geo(rng));
            if (th + 2 < n) cola.emplace_back(th, k < km ? s : -s);
        }
    }
    for (const auto& [t, s] : cola) {
        eps[t] += s * P.a1 * sd_e;
        if (t + 1 < n) eps[t + 1] += s * P.a2 * sd_e;
    }
    if (choques) *choques = cola.size();
    std::vector<double> x(n);
    double cur = 0.0;
    for (std::size_t t = 0; t < n; ++t) {
        const double ph = P.tar ? (usar_tar && std::abs(cur) > P.umbral_tar ? 0.96 : 0.35) : P.phi;
        cur = ph * cur + eps[t];
        x[t] = cur;
    }
    return x;
}

struct Resultado {
    EstadoReplicas estado = EstadoReplicas::SIN_AJUSTAR;
    double rho = NaN, rho_inf = NaN, rho_sup = NaN, p_lr = NaN;
    double nc = NaN, nc_inf = NaN, nc_sup = NaN, nd = NaN, nd_inf = NaN, nd_sup = NaN;
    std::size_t alarmas = 0, bib0 = 0, bib1 = 0, iter = 0;
    bool convergio = false;
    double rho_aux = NaN;
    bool lh_utilizable = false;
    double p_iid = NaN, p_ons = NaN, p_lh = NaN, s_lh = NaN, s_ons = NaN;
    std::size_t ev_lh = 0, ev_ons = 0;
    // R2b: cortes
    std::size_t reales_registrados = 0, reales_cortados = 0, nulo_detectadas = 0, nulo_cortadas = 0;
    bool holm_lh_uno = false, familia_valida = false;
};

Resultado correr(const Parametros& P, std::uint64_t sem, std::unique_ptr<LilH>* guardar = nullptr, std::unique_ptr<ObsH>* guardar_d = nullptr) {
    Resultado R{};
    const auto x_cal = generar(P, P.n_cal, false, sem);
    const auto x_test = generar(P, P.n_test, true, sem + 7777);
    Parametros L = P;
    L.tasa_madre = 0.0;
    const auto x_limpia = generar(L, P.n_cal, false, sem + 3333);
    ConfiguracionObservador cd{};
    cd.arl0_pasos_nativos = P.arl0;
    cd.medir_tau50 = false;
    cd.escalas_max = 1;
    auto d = std::make_unique<ObsH>();
    d->calibrar(cd, x_limpia);
    auto lil = std::make_unique<LilH>();
    ConfiguracionLilyth cl{};
    cl.representacion = ModoRepresentacionLilyth::CRUDA;
    cl.b_max_nulo = 6;
    cl.hilos_nulo = 1;
    cl.semilla_nulo = sem * 31 + 5;
    if (!lil->calibrar(*d, cl, cd, x_cal)) return R;
    const auto& fr = lil->ficha_replicas();
    R.estado = fr.estado; R.rho = fr.rho; R.rho_inf = fr.rho_inf; R.rho_sup = fr.rho_sup; R.p_lr = fr.p_lr;
    R.nc = fr.n_signo_contrario; R.nc_inf = fr.n_contrario_inf; R.nc_sup = fr.n_contrario_sup;
    R.nd = fr.n_mismo_signo; R.nd_inf = fr.n_mismo_signo_inf; R.nd_sup = fr.n_mismo_signo_sup;
    R.alarmas = fr.n_alarmas;
    R.iter = fr.iter_indirecta; R.convergio = fr.convergio; R.rho_aux = fr.auxiliar.rho;
    auto d_run = std::make_unique<ObsH>(*d);
    d_run->reiniciar_flujo();
    for (double v : x_test) {
        const auto s = d_run->paso(v, NaN, true);
        (void)lil->paso(s, v, true);
    }
    lil->completar_nulos();
    const auto& fn = lil->ficha_nulo_cruda();
    R.bib0 = fn.biblioteca_innov[0].size(); R.bib1 = fn.biblioteca_innov[1].size();
    R.lh_utilizable = fn.lh_utilizable;
    R.p_iid = fn.contraste_iid.p_bilateral; R.p_ons = fn.contraste_onsager.p_bilateral; R.p_lh = fn.contraste_lh.p_bilateral;
    R.s_lh = fn.contraste_lh.estadistico_s; R.s_ons = fn.contraste_onsager.estadistico_s;
    R.ev_lh = fn.contraste_lh.n_eventos; R.ev_ons = fn.contraste_onsager.n_eventos;
    const auto& aj = lil->aalen_johansen_cruda().total;
    R.reales_registrados = aj.n_eventos_registrados; R.reales_cortados = aj.n_cortados_replica;
    for (const auto& r : lil->cerradas_historico()) {
        R.nulo_detectadas += r.cruda.nulo_lh.n_nulos_b + r.cruda.nulo_lh.cortadas_b;
        R.nulo_cortadas += r.cruda.nulo_lh.cortadas_b;
    }
    const auto& arm = lil->resumen_armonico_cruda();
    R.familia_valida = arm.familia_holm_valida;
    R.holm_lh_uno = arm.lh_no_decidible && arm.p_crudos[4] == 1.0 && arm.p_crudos[5] == 1.0 && arm.n_familia == 6;
    if (guardar) *guardar = std::move(lil);
    if (guardar_d) *guardar_d = std::move(d);
    return R;
}

template <class F>
std::vector<Resultado> semillas(const Parametros& P, std::size_t S, std::uint64_t base, F&& informe) {
    std::vector<Resultado> out(S);
    ark::metrologia::familia::en_paralelo(S, 0, [&](std::size_t i) { out[i] = correr(P, base + 101 * i); });
    for (std::size_t i = 0; i < S; ++i) informe(i, out[i]);
    return out;
}

void imprimir(std::size_t i, const Resultado& r) {
    std::printf("         [%2zu] est=%d alarmas=%zu rho_aux=%.3f it=%zu conv=%d rho=%.3f [%.3f, %.3f] p_lr=%.4f N_mismo=%.3f [%.3f, %.3f] N_contr=%.3f [%.3f, %.3f] "
                "bib=%zu/%zu LH=%d p_iid=%.3f p_ons=%.4f (S %+.3f) p_lh=%.4f (S %+.3f) ev=%zu/%zu cortes real %zu/%zu nulo %zu/%zu\n",
                i, int(r.estado), r.alarmas, r.rho_aux, r.iter, int(r.convergio), r.rho, r.rho_inf, r.rho_sup, r.p_lr, r.nd, r.nd_inf, r.nd_sup, r.nc, r.nc_inf, r.nc_sup, r.bib0,
                r.bib1, int(r.lh_utilizable), r.p_iid, r.p_ons, r.s_ons, r.p_lh, r.s_lh, r.ev_ons, r.ev_lh, r.reales_cortados, r.reales_registrados,
                r.nulo_cortadas, r.nulo_detectadas);
}

}  // namespace

int main(int argc, char** argv) {
    const double alfa = 0.05;
    // Numero de semillas de los controles de tamano (diseno del control, declarado): con S = 20 la banda binomial superior de
    // rechazos a alfa = 0.05 es 2.6.
    std::size_t S = 20;
    bool solo_r2 = false;   // --solo-r2: diagnostico, solo los dos mundos de R2
    bool solo_r1 = false;   // --solo-r1: diagnostico, solo el mundo sin replicas
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--rapido") == 0) S = 4;
        else if (std::strcmp(argv[i], "--solo-r2") == 0) solo_r2 = true;
        else if (std::strcmp(argv[i], "--solo-r1") == 0) solo_r1 = true;
        else if (std::strncmp(argv[i], "--S=", 4) == 0) S = static_cast<std::size_t>(std::atoi(argv[i] + 4));
    }

    // R1: AR(1) sin replicas
    if (!solo_r2) {
        Parametros P{};
        P.n_mismo = 0.0; P.n_contrario = 0.0;
        std::printf("      R1: AR(1) phi 0.88 con choques Poisson SIN replicas, %zu semillas\n", S);
        const auto res = semillas(P, S, 11000, imprimir);
        std::size_t decid = 0, uno = 0, onsr = 0;
        for (const auto& r : res) { decid += r.estado == EstadoReplicas::DECIDIBLE; uno += r.holm_lh_uno; onsr += r.p_ons < alfa; }
        std::printf("      R1: DECIDIBLE %zu/%zu (banda %.1f); familia de 6 con p_LH = 1 en %zu/%zu; Onsager rechaza %zu/%zu\n", decid, S,
                    banda(S, alfa), uno, S, onsr, S);
        chk(double(decid) <= banda(S, alfa) && uno + decid == S,
            "R1", "sin replicas el ajuste no pasa su nulo (<= alfa + banda) y L+H entra en Holm con p = 1");
    }
    if (solo_r1) return fallos ? 1 : 0;

    // R2: mismo signo, n = 0.3 y 0.6
    for (double n_m : {0.3, 0.6}) {
        Parametros P{};
        P.n_mismo = n_m; P.n_contrario = 0.0;
        std::printf("      R2: replicas de mismo signo n = %.1f, %zu semillas\n", n_m, S);
        const auto res = semillas(P, S, n_m < 0.5 ? 21000 : 22000, imprimir);
        std::size_t dec = 0, cubre = 0, contr0 = 0, lh_rech = 0, lh_ev = 0, ons_rech = 0, ons_ev = 0;
        for (const auto& r : res) {
            dec += r.estado == EstadoReplicas::DECIDIBLE;
            cubre += r.convergio;
            contr0 += (r.nc_inf <= 0.0);
            if (std::isfinite(r.p_lh)) { ++lh_ev; lh_rech += r.p_lh < alfa; }
            if (std::isfinite(r.p_ons)) { ++ons_ev; ons_rech += r.p_ons < alfa; }
        }
        std::printf("      R2 n=%.1f (verdad por choque): DECIDIBLE %zu/%zu; converge %zu/%zu; contrario compatible con 0 %zu/%zu; "
                    "L+H rechaza %zu/%zu (banda %.1f); Onsager rechaza %zu/%zu\n",
                    n_m, dec, S, cubre, S, contr0, S, lh_rech, lh_ev, banda(lh_ev, alfa), ons_rech, ons_ev);
        chk(dec == S && cubre == S && double(S - contr0) <= banda(S, alfa) && lh_ev > 0 && double(lh_rech) <= banda(lh_ev, alfa),
            n_m < 0.5 ? "R2a" : "R2b'",
            "replicas de mismo signo: pasa su nulo, la inferencia converge, cruzada ~ 0 y L+H con tamano <= alfa");
        if (n_m > 0.5) std::printf("      R2 (reporte): Onsager rechaza en n = 0.6 en %zu de %zu semillas (S1 lo vio en una)\n", ons_rech, ons_ev);
    }

    if (solo_r2) {
        if (fallos) std::printf("FALLA: %d fallas (solo R2)\n", fallos);
        else std::printf("PASA (solo R2)\n");
        return fallos ? 1 : 0;
    }

    // R2b: signo contrario
    {
        Parametros P{};
        P.n_mismo = 0.0; P.n_contrario = 0.4;
        std::printf("      R2b: replicas de signo contrario n = 0.4, %zu semillas\n", S);
        const auto res = semillas(P, S, 23000, imprimir);
        std::size_t cubre = 0, dec = 0;
        std::size_t rr = 0, rc = 0, nd = 0, nc = 0;
        for (const auto& r : res) {
            dec += r.estado == EstadoReplicas::DECIDIBLE;
            cubre += (r.nc_inf > 0.0);
            rr += r.reales_registrados; rc += r.reales_cortados; nd += r.nulo_detectadas; nc += r.nulo_cortadas;
        }
        const double p1 = rr ? double(rc) / double(rr) : NaN, p2 = nd ? double(nc) / double(nd) : NaN;
        const double pp = (rr + nd) ? double(rc + nc) / double(rr + nd) : NaN;
        const double ee = std::sqrt(pp * (1 - pp) * (1.0 / double(std::max<std::size_t>(rr, 1)) + 1.0 / double(std::max<std::size_t>(nd, 1))));
        const double z = est::norm_ppf(1.0 - alfa / 2.0);
        std::printf("      R2b: DECIDIBLE %zu/%zu; cruzada estructural que excluye el 0 %zu/%zu; cortados reales %zu/%zu = %.3f, "
                    "en el nulo L+H %zu/%zu = %.3f (banda +-%.3f)\n",
                    dec, S, cubre, S, rc, rr, p1, nc, nd, p2, z * ee);
        chk(dec == S && double(S - cubre) <= banda(S, alfa) && std::abs(p1 - p2) <= z * ee,
            "R2b", "replicas de signo contrario: la cruzada estructural excluye el 0 y cortes del nulo L+H = cortes reales en su banda");
    }

    // R3: TAR no lineal con replicas
    {
        Parametros P{};
        P.tar = true; P.n_mismo = 0.5; P.n_contrario = 0.0; P.a1 = 5.0; P.a2 = 5.0; P.arl0 = 400.0;
        std::printf("      R3: TAR con replicas de mismo signo n = 0.5, %zu semillas\n", S);
        const auto res = semillas(P, S, 24000, imprimir);
        std::size_t rech = 0, ev = 0, pos = 0;
        for (const auto& r : res) if (std::isfinite(r.p_lh)) { ++ev; rech += r.p_lh < alfa; pos += r.s_lh > 0.0; }
        std::printf("      R3: L+H rechaza %zu/%zu, con S > 0 en %zu\n", rech, ev, pos);
        chk(ev > 0 && 2 * rech > ev, "R3", "TAR con replicas: L+H rechaza en la mayoria de semillas (la no linealidad no se la come la autoexcitacion)");
    }

    // R6, R7, R8 sobre un mundo con replicas
    {
        Parametros P{};
        P.n_mismo = 0.6;
        std::unique_ptr<LilH> lil;
        std::unique_ptr<ObsH> d;
        const auto x_cal = generar(P, P.n_cal, false, 31000);
        Parametros L = P;
        L.tasa_madre = 0.0;
        const auto x_limpia = generar(L, P.n_cal, false, 31333);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = P.arl0; cd.medir_tau50 = false; cd.escalas_max = 1;
        d = std::make_unique<ObsH>();
        d->calibrar(cd, x_limpia);
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA; cl.b_max_nulo = 6; cl.hilos_nulo = 0; cl.semilla_nulo = 77;
        auto l1 = std::make_unique<LilH>();
        auto l2 = std::make_unique<LilH>();
        l1->calibrar(*d, cl, cd, x_cal);
        l2->calibrar(*d, cl, cd, x_cal);
        const FichaReplicas f_cal = l1->ficha_replicas();
        const auto bib_cal = l1->ficha_nulo_cruda().biblioteca_innov;
        // dos tramos medidos distintos
        for (int k = 0; k < 2; ++k) {
            auto& lk = k == 0 ? l1 : l2;
            auto dd = std::make_unique<ObsH>(*d);
            dd->reiniciar_flujo();
            Parametros Q = P;
            if (k == 1) { Q.n_mismo = 0.0; Q.tasa_madre = 1.0 / 150.0; }
            const auto xt = generar(Q, P.n_test, true, 32000 + k);
            for (double v : xt) { const auto s = dd->paso(v, NaN, true); (void)lk->paso(s, v, true); }
            lk->completar_nulos();
        }
        auto misma_ficha = [](const FichaReplicas& a, const FichaReplicas& b) {
            bool ok = a.estado == b.estado && igual(a.lr, b.lr) && igual(a.p_lr, b.p_lr) && igual(a.rho, b.rho) && igual(a.rho_sup, b.rho_sup);
            ok = ok && a.hawkes.mu.size() == b.hawkes.mu.size();
            for (std::size_t i = 0; ok && i < a.hawkes.mu.size(); ++i) ok = igual(a.hawkes.mu[i], b.hawkes.mu[i]);
            for (std::size_t m = 0; m < 2; ++m)
                for (std::size_t s = 0; s < 2; ++s)
                    for (std::size_t j = 0; j < REPLICAS_J_MAX; ++j) ok = ok && igual(a.hawkes.a[m][s][j], b.hawkes.a[m][s][j]);
            return ok;
        };
        auto misma_bib = [](const auto& a, const auto& b) {
            bool ok = true;
            for (std::size_t s = 0; s < 2; ++s) {
                ok = ok && a[s].size() == b[s].size();
                for (std::size_t i = 0; ok && i < a[s].size(); ++i) {
                    ok = a[s][i].size() == b[s][i].size();
                    for (std::size_t k = 0; ok && k < a[s][i].size(); ++k) ok = igual(a[s][i][k], b[s][i][k]);
                }
            }
            return ok;
        };
        const bool r6 = misma_ficha(f_cal, l1->ficha_replicas()) && misma_ficha(f_cal, l2->ficha_replicas()) &&
                        misma_bib(bib_cal, l1->ficha_nulo_cruda().biblioteca_innov) && misma_bib(bib_cal, l2->ficha_nulo_cruda().biblioteca_innov);
        std::printf("      R6: estado %d, rho %.3f, biblioteca %zu/%zu; cerradas en los dos tramos %zu y %zu\n", int(f_cal.estado), f_cal.rho,
                    bib_cal[0].size(), bib_cal[1].size(), l1->cerradas_historico().size(), l2->cerradas_historico().size());
        chk(r6 && f_cal.estado == EstadoReplicas::DECIDIBLE, "R6", "la ficha de replicas y las marcas son las de calibrar y no dependen del tramo medido");

        // R7: inyecciones aisladas con hijas en el mundo de replicas
        const LilH::Ranura* ref = nullptr;
        for (const auto& r : l1->cerradas_historico())
            if (r.cruda.desenlace == DesenlaceRanura::RESTITUIDA && r.n_perfil_cruda > 0) { ref = &r; break; }
        bool r7 = ref != nullptr;
        if (ref) {
            std::size_t sol = 0, prev = 0, post = 0, hijas = 0, det = 0;
            for (int rep = 0; rep < 40; ++rep) {
                const std::array<std::int64_t, 1> tiempos{400};
                auto res = l1->evaluar_inyecciones_realizacion(
                    false, std::span<const double>(ref->hist_prev_cruda_est.data(), ref->n_hist_prev_cruda),
                    std::span<const double>(ref->perfil_cruda_est.data(), ref->n_perfil_cruda), ref->signo, ref->w_ancla,
                    ref->guarda_ancla, tiempos, 900, 5000 + rep, false, true, true);
                sol += res.solapes; prev += res.con_ajena_previa; post += res.con_reapertura_posterior; hijas += res.hijas; det += res.detectadas;
            }
            std::printf("      R7: 40 inyecciones aisladas con hijas: hijas %zu, solapes %zu, previas %zu, reaperturas posteriores %zu, admitidas %zu\n",
                        hijas, sol, prev, post, det);
            r7 = hijas > 0 && sol == prev && post > 0 && det > 0;
        }
        chk(r7, "R7", "las hijas son parte de su inyeccion: solapes = aperturas previas y las reaperturas posteriores no contaminan");

        // R8: identidad con a = 0
        bool r8 = ref != nullptr;
        if (ref) {
            const auto& fn = l1->ficha_nulo_cruda();
            AjusteHawkesReplicas cero = l1->ficha_replicas().hawkes;
            for (auto& x : cero.a) for (auto& y : x) y.fill(0.0);
            const detalle_nulo::GeneradorHijasNulo gen0{&cero, &fn.biblioteca_innov};
            ReferenciaNuloEvento ra{}, rb{};
            std::vector<EventoNuloLilyth<TMAXC>> ea, eb;
            std::size_t sa = 0, sb = 0;
            const std::span<const double> hp(ref->hist_prev_cruda_est.data(), ref->n_hist_prev_cruda);
            const std::span<const double> pf(ref->perfil_cruda_est.data(), ref->n_perfil_cruda);
            detalle_nulo::ejecutar_via_b_evento<TMAXC, RasgosDelta<ObsH>::b_max>(
                *d, *l1, Integracion::FLUJO, false, false, cl.nulo_sin_memoria, fn.segmentos_calma, fn.modelo_lineal,
                l1->contexto_purga_nulo(), NaN, hp, pf, ref->signo, ref->w_ancla, ref->guarda_ancla, fn.m_min,
                l1->configuracion().t_max_efectivo, 99, ra, ea, sa);
            detalle_nulo::ejecutar_via_b_evento<TMAXC, RasgosDelta<ObsH>::b_max>(
                *d, *l1, Integracion::FLUJO, false, false, cl.nulo_sin_memoria, fn.segmentos_calma, fn.modelo_lineal,
                l1->contexto_purga_nulo(), NaN, hp, pf, ref->signo, ref->w_ancla, ref->guarda_ancla, fn.m_min,
                l1->configuracion().t_max_efectivo, 99, rb, eb, sb, &gen0);
            r8 = ea.size() == eb.size() && ra.intentos_b == rb.intentos_b && ra.n_nulos_b == rb.n_nulos_b && sa == sb && ea.size() > 0;
            for (std::size_t i = 0; r8 && i < ea.size(); ++i) {
                r8 = igual(ea[i].tau_rest_consolidado, eb[i].tau_rest_consolidado) && igual(ea[i].dy, eb[i].dy);
                for (std::size_t k = 0; r8 && k < TMAXC; ++k) r8 = igual(ea[i].eta[k], eb[i].eta[k]);
            }
            std::printf("      R8: via (b) Onsager %zu eventos en %zu intentos; L+H con a = 0 %zu en %zu\n", ea.size(), ra.intentos_b, eb.size(),
                        rb.intentos_b);
        }
        chk(r8, "R8", "con a = 0 la via (b) L+H es la de Onsager bit a bit");
    }

    // D-ARFIMA: la via (a) con una calibracion larga y nulo ARFIMA no puede dejar ceros despues de la pila de 16 384 de `simular`.
    {
        const std::size_t n = 40000;
        const auto x = generar(Parametros{}, n, false, 41000);
        ConfiguracionModeloLineal cm{};
        cm.alfa = alfa; cm.modo = ModoNuloLineal::ARFIMA_WHITTLE; cm.d_externa = 0.3;
        const auto m = ajustar({x}, cm);
        std::mt19937_64 rng(42);
        std::vector<double> est;
        detalle_nulo::construir_estacionaria_nula(true, LILYTH_NULO_PERMUTACION, {x}, m, n, rng, est);
        std::size_t ceros = 0;
        for (std::size_t t = 16384; t < n; ++t) ceros += (est[t] == 0.0);
        std::printf("      D-ARFIMA: modelo %s d=%.2f, ceros exactos despues de 16 384: %zu de %zu\n", m.valido ? "valido" : "INVALIDO", m.d, ceros,
                    n - 16384);
        chk(m.valido && ceros == 0, "D-AR", "la via (a) ARFIMA con calibracion larga no deja la cola en cero");
    }

    if (fallos) std::printf("FALLA: %d fallas\n", fallos);
    else std::printf("PASA\n");
    return fallos ? 1 : 0;
}
