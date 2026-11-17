// sonda_a_perfil.cpp -- sonda (a) del sesgo comun de los nulos (2026-11-11, Code, orden del Comandante).
// Pregunta: el perfil inyectado en los nulos lleva dentro las replicas reales (doble conteo en L+H)? El perfil va del punto de cambio
// estimado a la alarma; las replicas del mundo llegan d0 = 6 + Geom pasos despues de la madre. Se mide la longitud de los perfiles de
// las ranuras cerradas en el mundo de replicas n = 0.6 (el de R2) y cuantas replicas reales caen dentro de algun perfil.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <random>
#include <vector>

#include "ark/lilyth/observador_lilyth.hpp"

using namespace ark::metrologia;
using namespace ark::lilyth;
using ObsH = ObservadorSaltos<24, 4, 120>;
using LilH = ObservadorLilyth<ObsH, 8, 128, 2048>;

// el mismo mundo que etapa_lh_lilyth.cpp (un paso), devolviendo ademas los tiempos de madres e hijas
std::vector<double> generar(double tasa, double n_mismo, std::size_t n, std::uint64_t sem, std::vector<std::size_t>& madres,
                            std::vector<std::size_t>& hijas) {
    const double phi = 0.88, sd = std::sqrt(1 - phi * phi);
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> nd(0, sd);
    std::uniform_real_distribution<double> u(0, 1);
    std::vector<double> e(n);
    for (auto& v : e) v = nd(rng);
    std::vector<std::pair<std::size_t, int>> cola;
    for (std::size_t t = 200; t + 200 < n; ++t)
        if (u(rng) < tasa) cola.emplace_back(t, u(rng) < 0.5 ? 1 : -1);
    const std::size_t n_madres = cola.size();
    std::geometric_distribution<int> geo(1.0 - std::exp(-1.0 / 4.0));
    std::poisson_distribution<int> pm(std::max(n_mismo, 1e-300));
    for (std::size_t c = 0; c < cola.size() && cola.size() < n; ++c) {
        const auto [t, s] = cola[c];
        const int k = n_mismo > 0 ? pm(rng) : 0;
        for (int i = 0; i < k; ++i) { const std::size_t th = t + 6 + geo(rng); if (th + 2 < n) cola.emplace_back(th, s); }
    }
    for (std::size_t c = 0; c < cola.size(); ++c) (c < n_madres ? madres : hijas).push_back(cola[c].first);
    for (const auto& [t, s] : cola) e[t] += s * 9.5 * sd;
    std::vector<double> x(n);
    double c = 0;
    for (std::size_t t = 0; t < n; ++t) { c = phi * c + e[t]; x[t] = c; }
    return x;
}

int main() {
    std::map<std::size_t, std::size_t> hist;
    std::size_t ranuras = 0, con_hija_dentro = 0, perfiles = 0, con_choque_dentro = 0;
    for (std::uint64_t sem : {22000ull, 22101ull, 22202ull}) {
        std::vector<std::size_t> mc, hc, mt, ht;
        const auto x_cal = generar(1.0 / 300, 0.6, 40000, sem, mc, hc);
        const auto x_test = generar(1.0 / 300, 0.6, 8000, sem + 7777, mt, ht);
        std::vector<std::size_t> m0, h0;
        const auto x_limpia = generar(0.0, 0.0, 40000, sem + 3333, m0, h0);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 800; cd.medir_tau50 = false; cd.escalas_max = 1;
        auto d = std::make_unique<ObsH>();
        d->calibrar(cd, x_limpia);
        auto lil = std::make_unique<LilH>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA; cl.b_max_nulo = 6; cl.hilos_nulo = 0; cl.semilla_nulo = sem * 31 + 5;
        lil->calibrar(*d, cl, cd, x_cal);
        auto dr = std::make_unique<ObsH>(*d);
        dr->reiniciar_flujo();
        for (double v : x_test) { const auto s = dr->paso(v, NaN, true); (void)lil->paso(s, v, true); }
        for (const auto& r : lil->cerradas_historico()) {
            ++ranuras;
            if (r.n_perfil_cruda == 0) continue;
            ++perfiles;
            ++hist[r.n_perfil_cruda];
            const std::int64_t fin = r.t_apertura, ini = fin - static_cast<std::int64_t>(r.n_perfil_cruda) + 1;
            bool dentro = false;
            for (std::size_t th : ht) dentro = dentro || (static_cast<std::int64_t>(th) >= ini && static_cast<std::int64_t>(th) <= fin);
            // la hija dentro del perfil solo cuenta si su madre (o una hija anterior) tambien cae en el perfil o antes: replica DENTRO
            con_hija_dentro += dentro;
            // control de alineacion: el perfil tiene que contener ALGUN choque real (madre o hija) si los indices son los del tramo
            bool choque = dentro;
            for (std::size_t tm : mt) choque = choque || (static_cast<std::int64_t>(tm) >= ini && static_cast<std::int64_t>(tm) <= fin);
            con_choque_dentro += choque;
        }
    }
    std::printf("ranuras cerradas %zu, con perfil %zu\n", ranuras, perfiles);
    std::printf("longitud del perfil (pasos): ");
    for (const auto& [l, c] : hist) std::printf("%zu:%zu ", l, c);
    std::printf("\nperfiles que contienen algun choque real (control de alineacion): %zu de %zu", con_choque_dentro, perfiles);
    std::printf("\nperfiles que contienen el instante de alguna replica real: %zu de %zu\n", con_hija_dentro, perfiles);
}
