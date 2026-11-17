// sustitutos.hpp -- nulos por sustitutos del bloque ARK-METRO-01, elegibles por el usuario (decision del Comandante, 2026-11-04).
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-04.
//
// Un nulo destruye una propiedad para contrastar si los datos la tienen; por eso NO se fija dentro del operador: el usuario escoge
// cuales aplica (mascara de bits), como escoge las familias de Holm. Cada nulo actua DENTRO de cada segmento contiguo (no mezcla
// tramos separados por un hueco).
//   PERMUTACION  conserva la marginal; destruye todo orden temporal (memoria, equilibrio local, memoria de beta).
//   BLOQUES      conserva el contenido de bloques de longitud b (equilibrio local dentro del bloque); destruye el orden entre
//                bloques (la memoria lenta). b la fija el contraste que lo use (en la separacion de escalas, b = T*).
//   SIGNOS       conserva |y - mu| en su orden (memoria de la volatilidad); destruye la direccion (y con ella toda
//                autocorrelacion lineal: no sirve de nulo de nivel para un flujo con memoria lineal).
// Fuera del ciclo caliente (reserva).
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace ark::metrologia {

enum NuloSustituto : std::uint32_t { NULO_NINGUNO = 0, NULO_PERMUTACION = 1u, NULO_BLOQUES = 2u, NULO_SIGNOS = 4u };
inline constexpr std::uint32_t NULOS_TODOS[3] = {NULO_PERMUTACION, NULO_BLOQUES, NULO_SIGNOS};
inline const char* nombre_nulo(std::uint32_t n) {
    return n == NULO_PERMUTACION ? "PERMUTACION" : n == NULO_BLOQUES ? "BLOQUES" : n == NULO_SIGNOS ? "SIGNOS" : "NINGUNO";
}

// Sustituto de unos segmentos. `bloque` solo lo usa BLOQUES (>= 1). `mu` es el centro de SIGNOS.
inline std::vector<std::vector<double>> sustituto(const std::vector<std::vector<double>>& segs, std::uint32_t nulo, std::size_t bloque,
                                                  double mu, std::mt19937_64& rng) {
    std::vector<std::vector<double>> out = segs;
    for (auto& g : out) {
        if (nulo == NULO_PERMUTACION) {
            std::shuffle(g.begin(), g.end(), rng);
        } else if (nulo == NULO_BLOQUES) {
            const std::size_t b = std::max<std::size_t>(1, bloque), nb = g.size() / b;
            if (nb < 2) continue;
            std::vector<std::size_t> orden(nb);
            for (std::size_t i = 0; i < nb; ++i) orden[i] = i;
            std::shuffle(orden.begin(), orden.end(), rng);
            std::vector<double> h(g.size());
            for (std::size_t i = 0; i < nb; ++i) std::copy_n(g.begin() + static_cast<std::ptrdiff_t>(orden[i] * b), b, h.begin() + static_cast<std::ptrdiff_t>(i * b));
            std::copy(g.begin() + static_cast<std::ptrdiff_t>(nb * b), g.end(), h.begin() + static_cast<std::ptrdiff_t>(nb * b));   // resto al final
            g.swap(h);
        } else if (nulo == NULO_SIGNOS) {
            std::bernoulli_distribution B(0.5);
            for (auto& v : g) v = mu + (B(rng) ? 1.0 : -1.0) * (v - mu);
        }
    }
    return out;
}

}  // namespace ark::metrologia
