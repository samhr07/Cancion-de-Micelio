// holm.hpp -- ARK 1.0: la UNICA definicion de Holm-Bonferroni del ARK (ARK-STAT-01-CXX, orden del Comandante, 2026-11-01).
//
// p_holm,(i) = min(1, max_{k <= i} (M - k + 1) p_(k)), con p_(1) <= ... <= p_(M), devuelto en la posicion ORIGINAL de cada p.
// Controla la tasa de error de la familia (FWER) bajo cualquier dependencia entre las M pruebas (Holm 1979).
//
// Contrato (falla dura): la familia entera vale NaN, con el mismo tamano M, si UN solo p no es finito (NaN, +Inf, -Inf) o cae
// fuera del intervalo cerrado [0, 1]. Nunca se sanea a 1 ni se recorta a 0: un p viciado corta en seco toda adopcion que
// dependa de la familia (cualquier comparacion con NaN da falso). Se valida ANTES de ordenar, asi que el comparador de
// std::sort nunca ve un NaN (con NaN no es un orden debil estricto y el resultado dependia de la posicion: acta
// AUDITORIA_HOLM_BONFERRONI.md, §A.1).
//
// Autocontenida (solo la biblioteca estandar) para que la prueba aislada y los mutantes la compilen sola.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace ark {

// true si p es un p valido: finito y en [0, 1]
inline bool p_valido(double p) { return std::isfinite(p) && p >= 0.0 && p <= 1.0; }

inline std::vector<double> holm(const std::vector<double>& p) {
    const size_t M = p.size();
    for (double x : p)
        if (!p_valido(x)) return std::vector<double>(M, std::numeric_limits<double>::quiet_NaN());
    std::vector<size_t> o(M);
    for (size_t k = 0; k < M; ++k) o[k] = k;
    std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return p[a] < p[b]; });
    std::vector<double> h(M);
    double mx = 0.0;
    for (size_t r = 0; r < M; ++r) {
        double v = std::min(1.0, static_cast<double>(M - r) * p[o[r]]);
        mx = std::max(mx, v);
        h[o[r]] = mx;
    }
    return h;
}

}  // namespace ark
