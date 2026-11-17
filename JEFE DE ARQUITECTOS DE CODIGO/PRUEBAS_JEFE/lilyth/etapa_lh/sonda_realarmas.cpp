// sonda_realarmas.cpp -- con choques AISLADOS (sin replicas) y el operador delta calibrado en limpio: cuantas alarmas trae cada
// choque, a que distancia, de que signo, y que hace Lilyth con ellas (apertura, continuacion, replica). Code, 2026-11-11. Diagnostico.
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

int main() {
    const double phi = 0.88, sd_e = std::sqrt(1 - phi * phi);
    const std::size_t n = 200000;
    std::mt19937_64 rng(7);
    std::normal_distribution<double> nd(0, sd_e);
    std::vector<double> el(n), e(n);
    for (std::size_t t = 0; t < n; ++t) el[t] = e[t] = nd(rng);
    std::vector<std::pair<std::size_t, int>> ch;
    for (std::size_t t = 1000; t + 1000 < n; t += 1000) { const int s = (t / 1000) % 2 ? 1 : -1; ch.emplace_back(t, s); e[t] += s * 9.5 * sd_e; e[t + 1] += s * 4.5 * sd_e; }
    std::vector<double> x(n), xl(n);
    double c = 0, cl = 0;
    for (std::size_t t = 0; t < n; ++t) { c = phi * c + e[t]; x[t] = c; cl = phi * cl + el[t]; xl[t] = cl; }
    ConfiguracionObservador cd{};
    cd.arl0_pasos_nativos = 800; cd.medir_tau50 = false; cd.escalas_max = 1;
    auto d = std::make_unique<ObsH>();
    d->calibrar(cd, xl);
    auto lil = std::make_unique<LilH>();
    ConfiguracionLilyth cfl{};
    cfl.representacion = ModoRepresentacionLilyth::CRUDA; cfl.t_max_efectivo = 128;
    lil->configurar(cfl, cd, *d, xl);
    auto dd = std::make_unique<ObsH>(*d);
    dd->reiniciar_flujo();
    struct A { std::size_t t; int s; bool ap, co, re; };
    std::vector<A> al;
    for (std::size_t t = 0; t < n; ++t) {
        const auto s = dd->paso(x[t], NaN, true);
        const auto o = lil->paso(s, x[t], true);
        if (o.alarma) al.push_back({t, o.signo_alarma, o.abrio_ranura, o.continuacion, o.replica});
    }
    // por choque: alarmas en [t, t + 200), por desfase
    std::map<std::size_t, std::size_t> hist;
    std::size_t mismo = 0, contrario = 0, lejos = 0, ap = 0, co = 0, re = 0, primera_mismo = 0;
    std::vector<std::size_t> por_choque(ch.size(), 0);
    for (const auto& a : al) {
        bool cerca = false;
        for (std::size_t k = 0; k < ch.size(); ++k) {
            if (a.t >= ch[k].first && a.t < ch[k].first + 200) {
                cerca = true;
                const std::size_t dt = a.t - ch[k].first;
                ++hist[dt / 10];
                ++por_choque[k];
                if (a.s == ch[k].second) ++mismo; else ++contrario;
                if (por_choque[k] == 1 && a.s == ch[k].second) ++primera_mismo;
                ap += a.ap; co += a.co; re += a.re;
                break;
            }
        }
        if (!cerca) ++lejos;
    }
    std::map<std::size_t, std::size_t> dist;
    for (auto v : por_choque) ++dist[v];
    std::printf("choques %zu, alarmas %zu; en [t, t+200): mismo signo %zu, contrario %zu; lejos de choques %zu (ARL0 800 -> esperadas ~%.0f en todo n)\n",
                ch.size(), al.size(), mismo, contrario, lejos, double(n) / 800.0);
    std::printf("cerca: aperturas %zu, continuaciones %zu, replicas %zu; primera alarma del choque con su signo: %zu\n", ap, co, re, primera_mismo);
    std::printf("alarmas por choque:"); for (auto& [k, v] : dist) std::printf(" %zu:%zu", k, v); std::printf("\n");
    std::printf("desfase (decenas de pasos):"); for (auto& [k, v] : hist) std::printf(" %zu:%zu", k * 10, v); std::printf("\n");
}
