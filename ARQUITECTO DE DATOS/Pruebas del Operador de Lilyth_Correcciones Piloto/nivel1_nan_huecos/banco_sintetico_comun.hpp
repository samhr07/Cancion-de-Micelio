// banco_sintetico_comun.hpp -- Base sintetica determinista compartida para las pruebas piloto de Lilyth (Niveles 1 a 4).
// Semillas, generadores y mundos (M1: AR(1) FLUJO, M2: ARFIMA(1,d,0) FLUJO, M3: NIVEL AR(1) con PurgaDiurna) coordinados.
#pragma once

#ifndef USAR_CABECERA_ORIGINAL
#include "ark/lilyth/observador_lilyth.hpp"
#endif

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

namespace ark::lilyth::piloto {

struct ParametrosBancoComun {
    std::size_t n_cal = 4096;
    std::size_t n_op = 12000;
    std::size_t separacion_choques = 450;
    std::size_t duracion_pulso = 4;
    double amp_pulso_sigma = 6.0;
    std::uint64_t semilla_base = 20261117ULL;
};

struct SerieSinteticaCoordinada {
    std::vector<double> x_cal;
    std::vector<std::size_t> fase_cal;
    std::vector<std::int64_t> dia_cal;

    std::vector<double> x_op_limpia;
    std::vector<double> eps_op_limpia;
    std::vector<std::size_t> fase_op;
    std::vector<std::int64_t> dia_op;
    std::vector<std::size_t> tiempos_choque;
    std::vector<int> signos_choque;
};

// Mundo M1: FLUJO AR(1) estacionario con choques de pulso conocidos en las innovaciones
inline SerieSinteticaCoordinada generar_mundo_m1_ar1_flujo(
    double phi = 0.80,
    const ParametrosBancoComun& p = {},
    std::uint64_t sal_semilla = 0) {
    SerieSinteticaCoordinada out{};
    const double sd_e = std::sqrt(std::max(1e-12, 1.0 - phi * phi));
    std::mt19937_64 rng(ark::lilyth::detalle_nulo::splitmix64(p.semilla_base ^ sal_semilla ^ 0xA1B2C3D4E5F60718ULL));
    std::normal_distribution<double> nd(0.0, sd_e);

    out.x_cal.resize(p.n_cal);
    out.fase_cal.resize(p.n_cal, 0);
    out.dia_cal.resize(p.n_cal, 0);
    double u = 0.0;
    for (std::size_t s = 0; s < 200; ++s) u = phi * u + nd(rng); // burn-in
    for (std::size_t t = 0; t < p.n_cal; ++t) {
        u = phi * u + nd(rng);
        out.x_cal[t] = u;
        out.fase_cal[t] = t % 24;
        out.dia_cal[t] = static_cast<std::int64_t>(t / 24);
    }

    out.x_op_limpia.resize(p.n_op);
    out.eps_op_limpia.resize(p.n_op);
    out.fase_op.resize(p.n_op, 0);
    out.dia_op.resize(p.n_op, 0);
    for (std::size_t t = 0; t < p.n_op; ++t) {
        out.eps_op_limpia[t] = nd(rng);
        out.fase_op[t] = (p.n_cal + t) % 24;
        out.dia_op[t] = static_cast<std::int64_t>((p.n_cal + t) / 24);
    }

    int signo = 1;
    for (std::size_t t0 = 350; t0 + 300 < p.n_op; t0 += p.separacion_choques) {
        out.tiempos_choque.push_back(t0);
        out.signos_choque.push_back(signo);
        for (std::size_t k = 0; k < p.duracion_pulso && (t0 + k) < p.n_op; ++k) {
            out.eps_op_limpia[t0 + k] += static_cast<double>(signo) * p.amp_pulso_sigma * sd_e;
        }
        signo = -signo;
    }

    for (std::size_t t = 0; t < p.n_op; ++t) {
        u = phi * u + out.eps_op_limpia[t];
        out.x_op_limpia[t] = u;
    }
    return out;
}

// Mundo M2: FLUJO ARFIMA(1, d, 0) con memoria larga fraccionaria d in (-0.5, 0.5)
inline SerieSinteticaCoordinada generar_mundo_m2_arfima_flujo(
    double d = 0.25,
    double phi1 = 0.40,
    const ParametrosBancoComun& p = {},
    std::uint64_t sal_semilla = 0) {
    SerieSinteticaCoordinada out{};
    const double sd_e = std::sqrt(std::max(1e-12, 1.0 - phi1 * phi1));
    std::mt19937_64 rng(ark::lilyth::detalle_nulo::splitmix64(p.semilla_base ^ sal_semilla ^ 0x55AA55AA12345678ULL));
    std::normal_distribution<double> nd(0.0, sd_e);

    const std::size_t burn = 1024;
    const std::size_t tot = burn + p.n_cal + p.n_op;
    std::vector<double> eps(tot, 0.0);
    for (std::size_t t = 0; t < tot; ++t) eps[t] = nd(rng);

    out.eps_op_limpia.resize(p.n_op);
    for (std::size_t t = 0; t < p.n_op; ++t) out.eps_op_limpia[t] = eps[burn + p.n_cal + t];

    int signo = 1;
    for (std::size_t t0 = 350; t0 + 300 < p.n_op; t0 += p.separacion_choques) {
        out.tiempos_choque.push_back(t0);
        out.signos_choque.push_back(signo);
        for (std::size_t k = 0; k < p.duracion_pulso && (t0 + k) < p.n_op; ++k) {
            const double imp = static_cast<double>(signo) * p.amp_pulso_sigma * sd_e;
            eps[burn + p.n_cal + t0 + k] += imp;
            out.eps_op_limpia[t0 + k] += imp;
        }
        signo = -signo;
    }

    std::vector<double> u(tot, 0.0);
    for (std::size_t t = 0; t < tot; ++t) {
        u[t] = eps[t] + (t >= 1 ? phi1 * u[t - 1] : 0.0);
    }

    const std::size_t max_k = 2048;
    std::vector<double> pi_inv(max_k, 0.0);
    ark::metrologia::elw::coeficientes(-d, max_k, pi_inv.data());

    auto int_frac = [&](std::size_t idx) noexcept -> double {
        double s = 0.0;
        const std::size_t lim = std::min(idx + 1, max_k);
        for (std::size_t k = 0; k < lim; ++k) s += pi_inv[k] * u[idx - k];
        return s;
    };

    out.x_cal.resize(p.n_cal);
    out.fase_cal.resize(p.n_cal, 0);
    out.dia_cal.resize(p.n_cal, 0);
    for (std::size_t t = 0; t < p.n_cal; ++t) {
        out.x_cal[t] = int_frac(burn + t);
        out.fase_cal[t] = t % 24;
        out.dia_cal[t] = static_cast<std::int64_t>(t / 24);
    }

    out.x_op_limpia.resize(p.n_op);
    out.fase_op.resize(p.n_op, 0);
    out.dia_op.resize(p.n_op, 0);
    for (std::size_t t = 0; t < p.n_op; ++t) {
        out.x_op_limpia[t] = int_frac(burn + p.n_cal + t);
        out.fase_op[t] = (p.n_cal + t) % 24;
        out.dia_op[t] = static_cast<std::int64_t>((p.n_cal + t) / 24);
    }
    return out;
}

// Mundo M3: NIVEL (incrementos AR(1) modulados por perfil diurno de 24 fases) para probar PurgaDiurna y cortes de segmento
inline SerieSinteticaCoordinada generar_mundo_m3_nivel_purga(
    double phi = 0.65,
    const ParametrosBancoComun& p = {},
    std::uint64_t sal_semilla = 0) {
    auto base = generar_mundo_m1_ar1_flujo(phi, p, sal_semilla ^ 0x9988776655443322ULL);
    auto factor_fase = [](std::size_t f) noexcept -> double {
        // Marea diurna determinista en [0.6, 1.6]
        return 1.0 + 0.5 * std::sin(2.0 * 3.14159265358979323846 * static_cast<double>(f % 24) / 24.0);
    };

    double nivel = 100.0;
    for (std::size_t t = 0; t < base.x_cal.size(); ++t) {
        nivel += base.x_cal[t] * factor_fase(base.fase_cal[t]);
        base.x_cal[t] = nivel;
    }
    for (std::size_t t = 0; t < base.x_op_limpia.size(); ++t) {
        nivel += base.x_op_limpia[t] * factor_fase(base.fase_op[t]);
        base.x_op_limpia[t] = nivel;
    }
    return base;
}

// ============================================================================
// Mundos H0 (Tamaño nominal: choque puntual de 1 paso que decae con la ley lineal
// exacta del mundo, z_tau = 0) y H1 (Potencia: añade cola lenta z_tau = amp_lenta * rho_lenta^tau)
// para las R = 200 réplicas con intervalo de Wilson (Punto 6 de la Auditoría).
// ============================================================================

struct ParametrosReplicaH0H1 {
    std::size_t n_cal = 3200;
    std::size_t n_op = 3200;
    std::size_t separacion_choques = 220;
    double amp_choque_sigma = 8.0;
    bool lenta_h1 = false;
    double amp_lenta_sigma = 4.0;
    double rho_lenta = 0.97;
    std::uint64_t semilla_base = 20261117ULL;
};

inline SerieSinteticaCoordinada generar_mundo_m1_h0_h1(
    double phi,
    const ParametrosReplicaH0H1& p,
    std::size_t idx_replica) {

    const std::uint64_t sem = p.semilla_base + 10007ULL * idx_replica;
    SerieSinteticaCoordinada out{};
    const double sd_e = std::sqrt(std::max(1e-12, 1.0 - phi * phi));
    std::mt19937_64 rng(ark::lilyth::detalle_nulo::splitmix64(sem ^ 0xA1B2C3D4E5F60718ULL));
    std::normal_distribution<double> nd(0.0, sd_e);

    out.x_cal.resize(p.n_cal);
    out.fase_cal.resize(p.n_cal, 0);
    out.dia_cal.resize(p.n_cal, 0);
    double u = 0.0;
    for (std::size_t s = 0; s < 200; ++s) u = phi * u + nd(rng);
    for (std::size_t t = 0; t < p.n_cal; ++t) {
        u = phi * u + nd(rng);
        out.x_cal[t] = u;
        out.fase_cal[t] = t % 24;
        out.dia_cal[t] = static_cast<std::int64_t>(t / 24);
    }

    out.x_op_limpia.resize(p.n_op);
    out.eps_op_limpia.resize(p.n_op);
    out.fase_op.resize(p.n_op, 0);
    out.dia_op.resize(p.n_op, 0);
    for (std::size_t t = 0; t < p.n_op; ++t) {
        out.eps_op_limpia[t] = nd(rng);
        out.fase_op[t] = (p.n_cal + t) % 24;
        out.dia_op[t] = static_cast<std::int64_t>((p.n_cal + t) / 24);
    }

    // Choque puntual de 1 solo paso en la innovación (bajo H0 evoluciona con la ley AR(1) exacta)
    int signo = 1;
    for (std::size_t t0 = 300; t0 + 250 < p.n_op; t0 += p.separacion_choques) {
        out.tiempos_choque.push_back(t0);
        out.signos_choque.push_back(signo);
        out.eps_op_limpia[t0] += static_cast<double>(signo) * p.amp_choque_sigma * sd_e;
        signo = -signo;
    }

    for (std::size_t t = 0; t < p.n_op; ++t) {
        u = phi * u + out.eps_op_limpia[t];
        out.x_op_limpia[t] = u;
    }

    // Si es H1, sumar componente lenta aditiva z_tau = signo * amp_lenta * sd_e * rho_lenta^tau
    if (p.lenta_h1) {
        for (std::size_t k = 0; k < out.tiempos_choque.size(); ++k) {
            const std::size_t t0 = out.tiempos_choque[k];
            const double sgn = static_cast<double>(out.signos_choque[k]);
            double z = sgn * p.amp_lenta_sigma * sd_e;
            for (std::size_t t = t0; t < p.n_op && (t - t0) < 160; ++t) {
                out.x_op_limpia[t] += z;
                z *= p.rho_lenta;
            }
        }
    }
    return out;
}

inline SerieSinteticaCoordinada generar_mundo_m2_h0_h1(
    double d,
    double phi1,
    const ParametrosReplicaH0H1& p,
    std::size_t idx_replica) {

    const std::uint64_t sem = p.semilla_base + 10007ULL * idx_replica;
    SerieSinteticaCoordinada out{};
    const double sd_e = std::sqrt(std::max(1e-12, 1.0 - phi1 * phi1));
    std::mt19937_64 rng(ark::lilyth::detalle_nulo::splitmix64(sem ^ 0x55AA55AA12345678ULL));
    std::normal_distribution<double> nd(0.0, sd_e);

    const std::size_t burn = 512;
    const std::size_t tot = burn + p.n_cal + p.n_op;
    std::vector<double> eps(tot, 0.0);
    for (std::size_t t = 0; t < tot; ++t) eps[t] = nd(rng);

    out.eps_op_limpia.resize(p.n_op);
    for (std::size_t t = 0; t < p.n_op; ++t) out.eps_op_limpia[t] = eps[burn + p.n_cal + t];

    int signo = 1;
    for (std::size_t t0 = 300; t0 + 250 < p.n_op; t0 += p.separacion_choques) {
        out.tiempos_choque.push_back(t0);
        out.signos_choque.push_back(signo);
        const double imp = static_cast<double>(signo) * p.amp_choque_sigma * sd_e;
        eps[burn + p.n_cal + t0] += imp;
        out.eps_op_limpia[t0] += imp;
        signo = -signo;
    }

    std::vector<double> u(tot, 0.0);
    for (std::size_t t = 0; t < tot; ++t) {
        u[t] = eps[t] + (t >= 1 ? phi1 * u[t - 1] : 0.0);
    }

    const std::size_t max_k = 1024;
    std::vector<double> pi_inv(max_k, 0.0);
    ark::metrologia::elw::coeficientes(-d, max_k, pi_inv.data());

    auto int_frac = [&](std::size_t idx) noexcept -> double {
        double s = 0.0;
        const std::size_t lim = std::min(idx + 1, max_k);
        for (std::size_t k = 0; k < lim; ++k) s += pi_inv[k] * u[idx - k];
        return s;
    };

    out.x_cal.resize(p.n_cal);
    out.fase_cal.resize(p.n_cal, 0);
    out.dia_cal.resize(p.n_cal, 0);
    for (std::size_t t = 0; t < p.n_cal; ++t) {
        out.x_cal[t] = int_frac(burn + t);
        out.fase_cal[t] = t % 24;
        out.dia_cal[t] = static_cast<std::int64_t>(t / 24);
    }

    out.x_op_limpia.resize(p.n_op);
    out.fase_op.resize(p.n_op, 0);
    out.dia_op.resize(p.n_op, 0);
    for (std::size_t t = 0; t < p.n_op; ++t) {
        out.x_op_limpia[t] = int_frac(burn + p.n_cal + t);
        out.fase_op[t] = (p.n_cal + t) % 24;
        out.dia_op[t] = static_cast<std::int64_t>((p.n_cal + t) / 24);
    }

    if (p.lenta_h1) {
        for (std::size_t k = 0; k < out.tiempos_choque.size(); ++k) {
            const std::size_t t0 = out.tiempos_choque[k];
            const double sgn = static_cast<double>(out.signos_choque[k]);
            double z = sgn * p.amp_lenta_sigma * sd_e;
            for (std::size_t t = t0; t < p.n_op && (t - t0) < 160; ++t) {
                out.x_op_limpia[t] += z;
                z *= p.rho_lenta;
            }
        }
    }
    return out;
}

inline SerieSinteticaCoordinada generar_mundo_m3_h0_h1(
    double phi,
    const ParametrosReplicaH0H1& p,
    std::size_t idx_replica) {

    auto base = generar_mundo_m1_h0_h1(phi, p, idx_replica + 50000);
    auto factor_fase = [](std::size_t f) noexcept -> double {
        return 1.0 + 0.5 * std::sin(2.0 * 3.14159265358979323846 * static_cast<double>(f % 24) / 24.0);
    };

    double nivel = 100.0;
    for (std::size_t t = 0; t < base.x_cal.size(); ++t) {
        nivel += base.x_cal[t] * factor_fase(base.fase_cal[t]);
        base.x_cal[t] = nivel;
    }
    for (std::size_t t = 0; t < base.x_op_limpia.size(); ++t) {
        nivel += base.x_op_limpia[t] * factor_fase(base.fase_op[t]);
        base.x_op_limpia[t] = nivel;
    }
    return base;
}

} // namespace ark::lilyth::piloto
