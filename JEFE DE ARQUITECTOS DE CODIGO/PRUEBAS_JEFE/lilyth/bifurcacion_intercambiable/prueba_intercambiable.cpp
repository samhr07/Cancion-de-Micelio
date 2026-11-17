// prueba_intercambiable.cpp -- XB.1 y EXB.v de la bifurcacion del contraste del resumen de Lilyth (Code, 2026-11-15).
// Verdad conocida, sin Lilyth: llama a detalle_nulo::evaluar_contraste_intercambiable (la rama nueva) y, para el lado a lado, arma la rama
// vieja como la arma el sobrepaso por evento (u del real contra sus nulos, pool = PIT de cada nulo contra los otros de su conjunto) y llama
// a detalle_nulo::evaluar_contraste_resumen. Expectativas selladas en EXPECTATIVAS_BIFURCACION.md (sha256 e9cfa2c9...).
// Uso: prueba_intercambiable.exe
#include "ark/lilyth/nulo_lilyth.hpp"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace dn = ark::lilyth::detalle_nulo;
using Estratos = std::vector<std::vector<double>>;

static int FALLAS = 0;
static void control(const char* id, bool ok, const std::string& det) {
    std::printf("[%s] %s %s\n", ok ? "PASA" : "FALLA", id, det.c_str());
    FALLAS += !ok;
}
static std::string g(double x) { char b[64]; std::snprintf(b, sizeof b, "%.6g", x); return b; }

static constexpr std::size_t B = 1999;
static constexpr double ALFA = 0.05;

// u de cada miembro contra el resto de su estrato (la misma cuenta que la rama nueva)
static std::vector<std::vector<double>> us(const Estratos& E) {
    std::vector<std::vector<double>> u;
    for (const auto& e : E) {
        std::vector<double> uk(e.size()), r;
        for (std::size_t m = 0; m < e.size(); ++m) {
            r.clear();
            for (std::size_t q = 0; q < e.size(); ++q) if (q != m) r.push_back(e[q]);
            uk[m] = dn::pit_punto_medio(e[m], r);
        }
        u.push_back(uk);
    }
    return u;
}

// rama vieja (POOL_LOO) tal como la arma contraste_sob
static ark::lilyth::ContrasteResumenNulo vieja(const Estratos& E, std::uint64_t sem) {
    std::vector<double> uk, dk, pool, r, sn;
    for (const auto& e : E) {
        sn.assign(e.begin() + 1, e.end());
        for (std::size_t i = 0; i < sn.size(); ++i) {
            r.clear();
            for (std::size_t j = 0; j < sn.size(); ++j) if (j != i) r.push_back(sn[j]);
            pool.push_back(dn::pit_punto_medio(sn[i], r));
        }
        uk.push_back(dn::pit_punto_medio(e[0], sn));
        dk.push_back(e[0] - ark::metrologia::est::cuantil(sn, 0.5));
    }
    return dn::evaluar_contraste_resumen(uk, dk, pool, B, sem);
}
static ark::lilyth::ContrasteResumenNulo nueva(const Estratos& E, std::uint64_t sem) {
    return dn::evaluar_contraste_intercambiable(E, B, sem, ALFA / 4.0);
}

// enumeracion completa de las configuraciones: p exacto del S observado y p de la configuracion mas extrema
static void enumerar(const Estratos& E, double& p_obs, double& p_ext) {
    const auto u = us(E);
    const std::size_t K = u.size();
    double s_obs = 0.0;
    for (const auto& uk : u) s_obs += uk[0] - 0.5;
    s_obs /= double(K);
    std::vector<std::size_t> idx(K, 0);
    std::vector<double> todos;
    double peso_total = 0.0;
    std::vector<double> pesos;
    for (;;) {
        double s = 0.0, w = 1.0;
        for (std::size_t k = 0; k < K; ++k) { s += u[k][idx[k]] - 0.5; w /= double(u[k].size()); }
        todos.push_back(std::abs(s / double(K)));
        pesos.push_back(w);
        peso_total += w;
        std::size_t k = 0;
        while (k < K && ++idx[k] == u[k].size()) { idx[k] = 0; ++k; }
        if (k == K) break;
    }
    double smax = 0.0;
    for (double v : todos) smax = std::max(smax, v);
    p_obs = 0.0; p_ext = 0.0;
    for (std::size_t i = 0; i < todos.size(); ++i) {
        if (todos[i] >= std::abs(s_obs) - 1e-15) p_obs += pesos[i];
        if (todos[i] >= smax - 1e-15) p_ext += pesos[i];
    }
    p_obs /= peso_total; p_ext /= peso_total;
    if (smax <= 1e-15) p_ext = 1.0;
}

// tau entera con empates y censura: 1 + floor(Exp(media 3)); > T_max = 6 se codifica en 7
static double tau_censurada(std::mt19937_64& rng) {
    std::exponential_distribution<double> ex(1.0 / 3.0);
    const double t = 1.0 + std::floor(ex(rng));
    return t > 6.0 ? 7.0 : t;
}

int main() {
    std::printf("Bifurcacion del contraste del resumen de Lilyth: controles con verdad conocida (XB.1, EXB.v), B = %zu, alfa = %g\n", B, ALFA);

    // a) enumeracion exacta
    {
        std::mt19937_64 rng(0xB1F0A0ULL);
        std::size_t mal_min = 0, mal_mc = 0;
        double peor_z = 0.0;
        for (int c = 0; c < 20; ++c) {
            const std::size_t K = 1 + rng() % 4;
            Estratos E(K);
            for (auto& e : E) { const std::size_t n = 2 + rng() % 4; for (std::size_t m = 0; m < n; ++m) e.push_back(double(1 + rng() % 4)); }
            double p_obs, p_ext;
            enumerar(E, p_obs, p_ext);
            const auto R = nueva(E, 1000 + c);
            mal_min += std::abs(R.p_min_alcanzable - std::max(p_ext, 1.0 / double(B + 1))) > 1e-12;
            const double sd = std::sqrt(std::max(p_obs * (1 - p_obs), 1e-300) / double(B));
            const double z = std::abs(R.p_bilateral - p_obs) / (sd + 1.0 / double(B + 1));
            peor_z = std::max(peor_z, z);
            mal_mc += z > 4.0;
        }
        control("XB.1a", mal_min == 0 && mal_mc == 0, "p_min distinto de la enumeracion " + std::to_string(mal_min) + "/20; p de Monte Carlo a > 4 sd " +
                                                       std::to_string(mal_mc) + "/20 (peor " + g(peor_z) + " sd)");
    }

    // b) tamano, y EXB.v (rama vieja en las mismas replicas)
    {
        const double techo = ALFA + 2.58 * std::sqrt(ALFA * (1 - ALFA) / 2000.0);
        bool ok = true;
        std::string det, det_v;
        bool vieja_anticons_k1 = false;
        for (std::size_t K : {std::size_t(1), std::size_t(3), std::size_t(10)}) {
            std::mt19937_64 rng(0xB1F0B0ULL + K);
            std::size_t rn = 0, rv = 0;
            for (int rep = 0; rep < 2000; ++rep) {
                Estratos E(K);
                for (auto& e : E) for (int m = 0; m < 9; ++m) e.push_back(tau_censurada(rng));
                rn += nueva(E, 50000 + rep).p_bilateral < ALFA;
                rv += vieja(E, 50000 + rep).p_bilateral < ALFA;
            }
            const double tn = rn / 2000.0, tv = rv / 2000.0;
            ok = ok && tn <= techo;
            det += "K " + std::to_string(K) + ": " + g(tn) + "; ";
            det_v += "K " + std::to_string(K) + ": " + g(tv) + "; ";
            if (K == 1) vieja_anticons_k1 = tv > techo;
        }
        control("XB.1b", ok, "tasa de rechazo de la intercambiable (techo " + g(techo) + "): " + det);
        std::printf("[%s] EXB.v rama vieja en las mismas replicas: %s(se esperaba > %s con K = 1)\n", vieja_anticons_k1 ? "SE CUMPLE" : "NO SE CUMPLE",
                    det_v.c_str(), g(techo).c_str());
    }

    // c) evento unico censurado con 8 nulos <= T_max
    {
        const Estratos E{{7.0, 1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 6.0, 6.0}};
        const auto R = nueva(E, 77);
        const auto V = vieja(E, 77);
        const double sd = std::sqrt((1.0 / 9) * (8.0 / 9) / double(B));
        control("XB.1c", R.p_min_alcanzable >= 1.0 / 9 - 1e-12 && !R.resoluble && R.p_bilateral >= 1.0 / 9 - 3 * sd,
                "p_min " + g(R.p_min_alcanzable) + ", resoluble " + std::to_string(int(R.resoluble)) + ", p " + g(R.p_bilateral) + " (rama vieja: p " +
                    g(V.p_bilateral) + ")");
    }

    // d) potencia
    {
        std::mt19937_64 rng(0xB1F0D0ULL);
        std::uniform_int_distribution<int> nul(1, 8), real(4, 11);
        std::size_t r = 0;
        for (int rep = 0; rep < 500; ++rep) {
            Estratos E(10);
            for (auto& e : E) { e.push_back(real(rng)); for (int m = 0; m < 8; ++m) e.push_back(nul(rng)); }
            r += nueva(E, 90000 + rep).p_bilateral < ALFA;
        }
        control("XB.1d", r >= 400, "rechazo " + g(r / 500.0) + " (>= 0.8)");
    }

    // e) degenerado
    {
        const Estratos E{{5, 5, 5}, {5, 5, 5, 5}, {5, 5}};
        const auto R = nueva(E, 5);
        control("XB.1e", R.estadistico_s == 0.0 && R.p_bilateral == 1.0 && R.p_min_alcanzable == 1.0,
                "S " + g(R.estadistico_s) + ", p " + g(R.p_bilateral) + ", p_min " + g(R.p_min_alcanzable));
    }

    // f) resolucion con un evento: 99 nulos (sellado: resoluble) y 50 (sellado: no)
    {
        auto con = [](int n) { std::vector<double> e{double(n + 1)}; for (int i = 1; i <= n; ++i) e.push_back(i); return Estratos{e}; };
        const auto R99 = nueva(con(99), 9), R50 = nueva(con(50), 9);
        control("XB.1f", R99.resoluble && !R50.resoluble, "99 nulos: p_min " + g(R99.p_min_alcanzable) + " resoluble " + std::to_string(int(R99.resoluble)) +
                                                           "; 50 nulos: p_min " + g(R50.p_min_alcanzable) + " resoluble " + std::to_string(int(R50.resoluble)));
        // declarado despues del sello: con p bilateral la configuracion mas extrema tiene DOS colas, p_min = 2/(n+1); hace falta n >= 159
        const auto R199 = nueva(con(199), 9), R159 = nueva(con(159), 9), R158 = nueva(con(158), 9);
        std::printf("[DATO] XB.1f' (posterior): p_min con 199 / 159 / 158 nulos = %s / %s / %s; resoluble %d / %d / %d (umbral alfa/4 = %s)\n",
                    g(R199.p_min_alcanzable).c_str(), g(R159.p_min_alcanzable).c_str(), g(R158.p_min_alcanzable).c_str(), int(R199.resoluble),
                    int(R159.resoluble), int(R158.resoluble), g(ALFA / 4).c_str());
    }

    std::printf("%d fallas\n", FALLAS);
    return FALLAS ? 1 : 0;
}
