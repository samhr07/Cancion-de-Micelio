// ksg_directo.hpp -- estimadores directos de informacion (KSG-1 y Frenzel-Pompe) y PIT aleatorizado, zero-heap. Jefe de Arquitectos
// de Codigo (Opus 5.5), 2026-11-05, por orden del Comandante (operador de coordinacion, cantidades APARTE con 4 busquedas de vecinos).
// Prerregistro: JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/coordinacion_sensor/EXPECTATIVAS_COORDINACION.md.
//
//   MI  (KSG-1, Kraskov-Stoegbauer-Grassberger 2004)   I(A;B)   = psi(k) + psi(n) - < psi(n_a + 1) + psi(n_b + 1) >
//   CMI (Frenzel-Pompe 2007)                            I(X;Y|Z) = psi(k) - < psi(n_xz + 1) + psi(n_yz + 1) - psi(n_z + 1) >
// con eps_i la distancia (norma del maximo) al k-esimo vecino en el espacio CONJUNTO y los conteos estrictamente menores que eps_i
// en cada subespacio: las mismas distancias valen para todos los terminos, y por eso los sesgos se cancelan (la TE como resta de dos
// KSG hechos por separado no los cancela: 1.5 veces el ruido en la sonda del 11-05).
// DIVERGE del modulo heredado (ark/coordinacion/estimador_ksg.hpp, NO AUDITADO): sin recorte en 0 (sesgaba el nulo hacia arriba),
// sin la copula de rangos exactos (rejilla (r + 0.5)/n: empates de DISTANCIA en L_inf, +0.036 nats con datos continuos). El MI es el
// mismo calculo que el heredado (digamma por la misma recurrencia, mismo orden de sumas): la prueba K2 lo exige bit a bit.
//
// PIT aleatorizado (decision del Comandante, 2026-11-03): u_i = (#{v < v_i} + U_i #{v = v_i}) / n, U_i ~ U(0, 1) independiente
// por punto. Conserva la dependencia (I(U_x; U_y) = I(X; Y) tambien con datos discretos) y no deja rejilla. El generador es
// splitmix64 (portable: no depende de la libreria estandar), sembrado por el que llama.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace ark::coordinacion {

inline constexpr double NaN_C = std::numeric_limits<double>::quiet_NaN();
inline constexpr double GAMMA_EULER_C = 0.577215664901532860606512;   // el mismo literal que el modulo heredado

// splitmix64 (Steele, Lea y Flood 2014): un paso y un generador de uniformes en [0, 1).
[[nodiscard]] constexpr std::uint64_t splitmix64(std::uint64_t x) noexcept {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

// Semilla de una ventana: mezcla encadenada de sus identificadores.
[[nodiscard]] constexpr std::uint64_t semilla_mezcla(std::uint64_t s, std::uint64_t a, std::uint64_t b, std::uint64_t c, std::uint64_t d) noexcept {
    std::uint64_t h = splitmix64(s);
    h = splitmix64(h ^ a);
    h = splitmix64(h ^ b);
    h = splitmix64(h ^ c);
    return splitmix64(h ^ d);
}

class GeneradorU01 {
    std::uint64_t e_ = 0;
public:
    constexpr explicit GeneradorU01(std::uint64_t semilla = 0) noexcept : e_(semilla) {}
    [[nodiscard]] constexpr std::uint64_t siguiente() noexcept { e_ += 0x9E3779B97F4A7C15ULL; std::uint64_t z = e_;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; return z ^ (z >> 31); }
    [[nodiscard]] constexpr double u01() noexcept { return static_cast<double>(siguiente() >> 11) * 0x1.0p-53; }
};

// PIT aleatorizado de v[0..n-1] en u[0..n-1]; idx es memoria de trabajo de n indices. Orden: por valor y, entre iguales, por indice
// (el mismo orden que un stable_sort; std::sort con desempate no reserva memoria). Las U se sacan en ese orden.
inline void pit_aleatorizado(const double* v, std::size_t n, double* u, std::size_t* idx, GeneradorU01& g) noexcept {
    for (std::size_t i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx, idx + n, [v](std::size_t a, std::size_t b) { return v[a] < v[b] || (v[a] == v[b] && a < b); });
    const double dn = static_cast<double>(n);
    for (std::size_t a = 0; a < n;) {
        std::size_t b = a;
        while (b < n && v[idx[b]] == v[idx[a]]) ++b;
        for (std::size_t j = a; j < b; ++j) u[idx[j]] = (static_cast<double>(a) + g.u01() * static_cast<double>(b - a)) / dn;
        a = b;
    }
}

// Estimadores directos sobre variables escalares (todas las de este operador lo son: pasado de una barra por escala).
template <std::size_t N_MAX>
class KSGDirecto {
    std::array<double, N_MAX> dist_{};
    std::array<double, N_MAX> eps_{};
    std::array<double, N_MAX + 2> psi_{};

public:
    KSGDirecto() noexcept {
        psi_[0] = NaN_C;
        psi_[1] = -GAMMA_EULER_C;
        for (std::size_t n = 2; n < psi_.size(); ++n) psi_[n] = psi_[n - 1] + 1.0 / static_cast<double>(n - 1);
    }
    [[nodiscard]] double psi(std::size_t n) const noexcept { return psi_[n]; }

    // I(a; b), KSG-1, norma del maximo. NaN si n < k + 2 (el estimador no esta definido).
    [[nodiscard]] double mi(const double* a, const double* b, std::size_t n, std::size_t k) noexcept {
        if (k < 1 || n < k + 2 || n > N_MAX) return NaN_C;
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t m = 0;
            for (std::size_t j = 0; j < n; ++j) {
                if (j == i) continue;
                const double da = std::abs(a[i] - a[j]), db = std::abs(b[i] - b[j]);
                dist_[m++] = da > db ? da : db;
            }
            std::nth_element(dist_.begin(), dist_.begin() + static_cast<std::ptrdiff_t>(k - 1), dist_.begin() + static_cast<std::ptrdiff_t>(n - 1));
            eps_[i] = dist_[k - 1];
        }
        double suma = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            std::uint32_t na = 0, nb = 0;
            for (std::size_t j = 0; j < n; ++j) {
                if (j == i) continue;
                if (std::abs(a[i] - a[j]) < eps_[i]) ++na;
                if (std::abs(b[i] - b[j]) < eps_[i]) ++nb;
            }
            suma += (psi_[na + 1] + psi_[nb + 1]);
        }
        return psi_[k] - suma / static_cast<double>(n) + psi_[n];
    }

    // I(x; y | z), Frenzel-Pompe, norma del maximo, eps del espacio conjunto (x, y, z). NaN si n < k + 2.
    [[nodiscard]] double cmi(const double* x, const double* y, const double* z, std::size_t n, std::size_t k) noexcept {
        if (k < 1 || n < k + 2 || n > N_MAX) return NaN_C;
        for (std::size_t i = 0; i < n; ++i) {
            std::size_t m = 0;
            for (std::size_t j = 0; j < n; ++j) {
                if (j == i) continue;
                const double dx = std::abs(x[i] - x[j]), dy = std::abs(y[i] - y[j]), dz = std::abs(z[i] - z[j]);
                const double dxy = dx > dy ? dx : dy;
                dist_[m++] = dxy > dz ? dxy : dz;
            }
            std::nth_element(dist_.begin(), dist_.begin() + static_cast<std::ptrdiff_t>(k - 1), dist_.begin() + static_cast<std::ptrdiff_t>(n - 1));
            eps_[i] = dist_[k - 1];
        }
        double suma = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            std::uint32_t nxz = 0, nyz = 0, nz = 0;
            for (std::size_t j = 0; j < n; ++j) {
                if (j == i) continue;
                const double dz = std::abs(z[i] - z[j]);
                if (!(dz < eps_[i])) continue;
                ++nz;
                if (std::abs(x[i] - x[j]) < eps_[i]) ++nxz;
                if (std::abs(y[i] - y[j]) < eps_[i]) ++nyz;
            }
            suma += (psi_[nxz + 1] + psi_[nyz + 1] - psi_[nz + 1]);
        }
        return psi_[k] - suma / static_cast<double>(n);
    }
};

}  // namespace ark::coordinacion
