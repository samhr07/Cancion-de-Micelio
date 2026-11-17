// sonda_deteccion.cpp -- cuanto detecta el operador delta calibrado SOBRE una historia con choques, segun la densidad de choques.
// Code, 2026-11-11. Diagnostico del mundo de la bateria L+H (no es un control).
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

int main() {
    const double phi = 0.88, sd_e = std::sqrt(1 - phi * phi);
    for (double tasa : {1.0 / 300, 1.0 / 1000, 1.0 / 3000}) {
        for (int cal_limpia = 0; cal_limpia < 2; ++cal_limpia) {
            const std::size_t n = 60000;
            std::mt19937_64 rng(5);
            std::normal_distribution<double> nd(0, sd_e);
            std::uniform_real_distribution<double> u(0, 1);
            std::vector<double> eps(n), eps_l(n);
            for (std::size_t t = 0; t < n; ++t) eps_l[t] = eps[t] = nd(rng);
            std::size_t choques = 0;
            for (std::size_t t = 200; t + 200 < n; ++t)
                if (u(rng) < tasa) { const int s = u(rng) < 0.5 ? 1 : -1; eps[t] += s * 9.5 * sd_e; eps[t + 1] += s * 4.5 * sd_e; ++choques; }
            std::vector<double> x(n), xl(n);
            double c = 0, cl = 0;
            for (std::size_t t = 0; t < n; ++t) { c = phi * c + eps[t]; x[t] = c; cl = phi * cl + eps_l[t]; xl[t] = cl; }
            ConfiguracionObservador cd{};
            cd.arl0_pasos_nativos = 800; cd.medir_tau50 = false; cd.escalas_max = 1;
            auto d = std::make_unique<ObsH>();
            d->calibrar(cd, cal_limpia ? xl : x);
            auto lil = std::make_unique<LilH>();
            ConfiguracionLilyth cfl{};
            cfl.representacion = ModoRepresentacionLilyth::CRUDA; cfl.t_max_efectivo = 128;
            lil->configurar(cfl, cd, *d, cal_limpia ? xl : x);
            auto dd = std::make_unique<ObsH>(*d);
            dd->reiniciar_flujo();
            std::size_t alarmas = 0;
            for (std::size_t t = 0; t < n; ++t) { const auto s = dd->paso(x[t], NaN, true); if (lil->paso(s, x[t], true).alarma) ++alarmas; }
            std::printf("tasa 1/%.0f  calibracion %s: choques %zu, alarmas %zu (%.2f por choque)\n", 1 / tasa, cal_limpia ? "SIN choques" : "CON choques",
                        choques, alarmas, double(alarmas) / double(choques));
        }
    }
}
