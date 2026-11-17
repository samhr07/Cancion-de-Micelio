// sonda_ii.cpp -- diagnostico de la inferencia indirecta del nulo L+H: una semilla por mundo, con la traza del punto fijo.
// Code, 2026-11-11. Compilar con /DARK_LILYTH_DEPURAR_II. Mundos: choques de UN paso (a2 = 0) sin replicas y con n = 0.6, y el de dos
// pasos sin replicas (el de R1 que rechazo).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

#include "ark/lilyth/observador_lilyth.hpp"

using namespace ark::metrologia;
using namespace ark::lilyth;
using ObsH = ObservadorSaltos<24, 4, 120>;
using LilH = ObservadorLilyth<ObsH, 8, 128, 2048>;

std::vector<double> generar(double tasa, double n_mismo, double a2, std::size_t n, std::uint64_t sem) {
    const double phi = 0.88, sd = std::sqrt(1 - phi * phi);
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> nd(0, sd);
    std::uniform_real_distribution<double> u(0, 1);
    std::vector<double> e(n);
    for (auto& v : e) v = nd(rng);
    std::vector<std::pair<std::size_t, int>> cola;
    for (std::size_t t = 200; t + 200 < n; ++t) if (u(rng) < tasa) cola.emplace_back(t, u(rng) < 0.5 ? 1 : -1);
    std::geometric_distribution<int> geo(1.0 - std::exp(-1.0 / 4.0));
    std::poisson_distribution<int> pm(std::max(n_mismo, 1e-300));
    for (std::size_t c = 0; c < cola.size() && cola.size() < n; ++c) {
        const auto [t, s] = cola[c];
        const int k = n_mismo > 0 ? pm(rng) : 0;
        for (int i = 0; i < k; ++i) { const std::size_t th = t + 6 + geo(rng); if (th + 2 < n) cola.emplace_back(th, s); }
    }
    for (const auto& [t, s] : cola) { e[t] += s * 9.5 * sd; if (t + 1 < n) e[t + 1] += s * a2 * sd; }
    std::vector<double> x(n);
    double c = 0;
    for (std::size_t t = 0; t < n; ++t) { c = phi * c + e[t]; x[t] = c; }
    return x;
}

int main() {
    struct M { const char* nombre; double n_mismo, a2; };
    for (const M m : {M{"un paso, sin replicas", 0.0, 0.0}, M{"un paso, n = 0.6", 0.6, 0.0}, M{"dos pasos, sin replicas", 0.0, 4.5}}) {
        const auto t0 = std::chrono::steady_clock::now();
        const std::size_t n = 40000;
        const auto x = generar(1.0 / 300, m.n_mismo, m.a2, n, 5);
        const auto xl = generar(0.0, 0.0, 0.0, n, 6);
        ConfiguracionObservador cd{};
        cd.arl0_pasos_nativos = 800; cd.medir_tau50 = false; cd.escalas_max = 1;
        auto d = std::make_unique<ObsH>();
        d->calibrar(cd, xl);
        auto lil = std::make_unique<LilH>();
        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA; cl.b_max_nulo = 6; cl.hilos_nulo = 0; cl.semilla_nulo = 9;
        std::printf("== %s ==\n", m.nombre);
        lil->calibrar(*d, cl, cd, x);
        const auto& f = lil->ficha_replicas();
        const double seg = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("   estado %d (%s) alarmas %zu p_lr %.4f rho_aux %.3f rho_est %.3f [%.3f, %.3f] iter %zu  %.0f s\n", int(f.estado), f.motivo,
                    f.n_alarmas, f.p_lr, f.auxiliar.rho, f.rho, f.rho_inf, f.rho_sup, f.iter_indirecta, seg);
        std::printf("   sd real por masa %.3f %.3f %.3f %.3f | N contrario %.3f [%.3f, %.3f]\n", f.sd_masa_real[0][0], f.sd_masa_real[0][1], f.sd_masa_real[1][0],
                    f.sd_masa_real[1][1], f.n_signo_contrario, f.n_contrario_inf, f.n_contrario_sup);
    }
}
