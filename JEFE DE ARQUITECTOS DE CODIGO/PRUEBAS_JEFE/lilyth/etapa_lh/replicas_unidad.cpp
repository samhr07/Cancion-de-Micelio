// replicas_unidad.cpp -- controles de unidad del proceso de replicas del nulo L+H (replicas_lilyth.hpp), con verdad conocida. C++20.
// Code, 2026-11-11.
//   U1   log-verosimilitud = referencia de fuerza bruta O(N^2) escrita aparte
//   U2   recuperacion: Hawkes de dos tipos con base por fase conocida; rho verdadero dentro de [cuantil alfa, cuantil 1 - alfa] del
//        bootstrap; KKT del optimo (gradiente nulo donde a > 0, no positivo donde a = 0)
//   U3   cumulo: media de descendientes de cada tipo = N (I - N)^-1 (teoria) dentro de 4 errores de Monte Carlo (caza M4 y M5)
//   U4   escalera: tau_k = 2^k hasta 2^K >= T_max (caza M8)
//   R1u  tamano del contraste del ajuste: Poisson homogeneo, la fraccion DECIDIBLE <= alfa + banda binomial (caza M9)
//   R4   Poisson modulado por fase SIN excitacion: con base por fase la fraccion DECIDIBLE <= alfa + banda; con base constante se
//        reporta la excitacion falsa (caza M2)
//   R5   ramificacion critica (rho = 1): NO_IDENTIFICADO (caza M7)
//   R2bu cruzada conocida: la masa de signo contrario verdadera cae en su intervalo basico (caza M3)
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "ark/lilyth/replicas_lilyth.hpp"

using namespace ark::lilyth;
namespace dr = ark::lilyth::detalle_replicas;

namespace {
int fallos = 0;
void chk(bool ok, const char* id, const char* txt) {
    std::printf("[%s] %-4s %s\n", ok ? "OK   " : "FALLA", id, txt);
    if (!ok) ++fallos;
}

// Hawkes verdadero de la bateria: un nucleo de una sola tau por par (madre, hija), base por bin.
AjusteHawkesReplicas verdad(std::size_t J_esc, const std::array<double, REPLICAS_J_MAX>& tau, std::size_t bins,
                            std::vector<double> mu, double n_mm, double n_mc, std::size_t j_kern) {
    AjusteHawkesReplicas h{};
    h.J = J_esc; h.tau = tau; h.bins = bins; h.mu = std::move(mu);
    const double q = std::exp(-1.0 / tau[j_kern]);
    for (std::size_t m = 0; m < 2; ++m) {
        h.a[m][m][j_kern] = n_mm * (1.0 - q);
        h.a[m][1 - m][j_kern] = n_mc * (1.0 - q);
    }
    dr::calcular_masas(h);
    return h;
}

// Referencia de fuerza bruta de la log-verosimilitud (definicion directa).
double ll_bruta(const std::vector<AlarmaReplica>& ev, std::size_t n_pasos, const std::vector<double>& pasos_por_bin,
                const AjusteHawkesReplicas& h) {
    double ll = 0.0;
    for (std::size_t i = 0; i < ev.size(); ++i) {
        const std::size_t s = ev[i].tipo;
        double lam = h.mu[s * h.bins + ev[i].bin];
        for (std::size_t k = 0; k < ev.size(); ++k) {
            if (!(ev[k].t < ev[i].t)) continue;
            for (std::size_t j = 0; j < h.J; ++j)
                lam += h.a[ev[k].tipo][s][j] * std::pow(std::exp(-1.0 / h.tau[j]), double(ev[i].t - ev[k].t - 1));
        }
        ll += std::log(lam);
    }
    for (std::size_t s = 0; s < 2; ++s)
        for (std::size_t b = 0; b < h.bins; ++b) ll -= h.mu[s * h.bins + b] * pasos_por_bin[b];
    for (const auto& e : ev) {
        for (std::size_t s = 0; s < 2; ++s)
            for (std::size_t j = 0; j < h.J; ++j) {
                double c = 0.0;
                for (std::int64_t t = e.t + 1; t < static_cast<std::int64_t>(n_pasos); ++t)
                    c += std::pow(std::exp(-1.0 / h.tau[j]), double(t - e.t - 1));
                ll -= h.a[e.tipo][s][j] * c;
            }
    }
    return ll;
}

struct Mundo {
    std::vector<std::uint32_t> bin_de_paso;
    std::vector<std::vector<std::int64_t>> pasos_de_bin;
    std::vector<double> pasos_por_bin;
};
Mundo mundo(std::size_t n, std::size_t bins, std::size_t periodo) {
    Mundo w{};
    w.bin_de_paso.resize(n);
    w.pasos_de_bin.assign(bins, {});
    w.pasos_por_bin.assign(bins, 0.0);
    for (std::size_t t = 0; t < n; ++t) {
        const std::uint32_t b = static_cast<std::uint32_t>((t / periodo) % bins);
        w.bin_de_paso[t] = b;
        w.pasos_de_bin[b].push_back(static_cast<std::int64_t>(t));
        w.pasos_por_bin[b] += 1.0;
    }
    return w;
}

// banda binomial superior de rechazos: alfa n + z_{1 - alfa} sqrt(n alfa (1 - alfa))
double banda(std::size_t n, double alfa) {
    return alfa * double(n) + ark::metrologia::est::norm_ppf(1.0 - alfa) * std::sqrt(double(n) * alfa * (1.0 - alfa));
}
}  // namespace

int main() {
    const double alfa = 0.05;
    std::array<double, REPLICAS_J_MAX> tau{};
    const std::size_t J = dr::escalera_tau(64, tau);

    // U4: escalera
    {
        std::array<double, REPLICAS_J_MAX> t1{}, t2{};
        const std::size_t j1 = dr::escalera_tau(128, t1), j2 = dr::escalera_tau(129, t2);
        bool ok = j1 == 8 && j2 == 9 && t1[0] == 1.0 && t1[7] == 128.0 && t2[8] == 256.0 && J == 7;
        for (std::size_t k = 0; k < j1; ++k) ok = ok && t1[k] == double(std::size_t{1} << k);
        chk(ok, "U4", "escalera tau_k = 2^k hasta 2^K >= T_max (128 -> 8 tau, 129 -> 9)");
    }

    // U1: log-verosimilitud contra la definicion
    {
        const auto w = mundo(3000, 2, 100);
        auto h = verdad(J, tau, 2, {0.004, 0.012, 0.006, 0.010}, 0.4, 0.15, 3);
        std::mt19937_64 rng(11);
        std::vector<AlarmaReplica> ev;
        dr::simular_historia(h, true, 3000, w.pasos_de_bin, 100000, rng, ev, w.bin_de_paso);
        const auto d = dr::preparar_datos(ev, 3000, w.pasos_por_bin, J, tau, false);
        const double a = dr::log_verosimilitud(d, h), b = ll_bruta(ev, 3000, w.pasos_por_bin, h);
        std::printf("      U1: %zu eventos, LL recursiva %.10f, bruta %.10f\n", ev.size(), a, b);
        chk(ev.size() > 20 && std::abs(a - b) <= 1e-9 * (1.0 + std::abs(b)), "U1", "log-verosimilitud = referencia de fuerza bruta");
    }

    // U2: recuperacion con base por fase y KKT
    {
        const std::size_t n = 200000;
        const auto w = mundo(n, 2, 500);
        auto h = verdad(J, tau, 2, {0.002, 0.008, 0.003, 0.006}, 0.45, 0.2, 3);
        std::mt19937_64 rng(21);
        std::vector<AlarmaReplica> ev;
        dr::simular_historia(h, true, n, w.pasos_de_bin, 1000000, rng, ev, w.bin_de_paso);
        ConfiguracionReplicas cr{};
        cr.alfa = alfa; cr.t_max = 64; cr.bins = 2; cr.semilla = 22; cr.hilos = 0;
        const auto f = ajustar_replicas(ev, w.bin_de_paso, cr);
        // KKT
        const auto d = dr::preparar_datos(ev, n, w.pasos_por_bin, J, tau, false);
        double kkt = 0.0;
        for (std::size_t m = 0; m < 2; ++m)
            for (std::size_t s = 0; s < 2; ++s)
                for (std::size_t j = 0; j < J; ++j) {
                    double g = -d.C[m][j];
                    for (std::size_t i = 0; i < d.n; ++i) {
                        if (d.tipo[i] != s) continue;
                        double lam = f.hawkes.mu[s * 2 + d.bin[i]];
                        for (std::size_t mm = 0; mm < 2; ++mm)
                            for (std::size_t jj = 0; jj < J; ++jj) lam += f.hawkes.a[mm][s][jj] * d.R[i * 2 * J + mm * J + jj];
                        g += d.R[i * 2 * J + m * J + j] / lam;
                    }
                    const double gr = g / d.C[m][j];
                    // interior: |g| ~ 0; borde (a ~ 0): g <= 0
                    if (f.hawkes.a[m][s][j] > 1e-7) kkt = std::max(kkt, std::abs(gr));
                    else kkt = std::max(kkt, std::max(0.0, gr));
                }
        std::printf("      U2: %zu eventos, rho verdad %.3f, ajuste %.3f, [%.3f, %.3f]; N(+,+) %.3f N(+,-) %.3f; mu %.4f %.4f %.4f %.4f; "
                    "KKT %.2e; iter %zu; p_lr %.4f; estado %d\n",
                    ev.size(), h.rho, f.rho, f.rho_inf, f.rho_sup, f.hawkes.N[0][0], f.hawkes.N[0][1], f.hawkes.mu[0], f.hawkes.mu[1],
                    f.hawkes.mu[2], f.hawkes.mu[3], kkt, f.hawkes.iteraciones, f.p_lr, int(f.estado));
        chk(f.estado == EstadoReplicas::DECIDIBLE && h.rho >= f.rho_inf && h.rho <= f.rho_sup && kkt <= 1e-3 && !f.hawkes.tope_iteraciones,
            "U2", "recupera rho dentro de su intervalo bootstrap, DECIDIBLE y optimo KKT");
    }

    // U3: cumulo contra la teoria N (I - N)^-1
    {
        auto h = verdad(J, tau, 1, {0.0, 0.0}, 0.4, 0.25, 2);
        // teoria
        const double a = h.N[0][0], b = h.N[0][1], c = h.N[1][0], dd = h.N[1][1];
        const double det = (1 - a) * (1 - dd) - b * c;
        const double inv00 = (1 - dd) / det, inv01 = b / det;
        const double e_pp = a * inv00 + b * (c / det), e_pm = a * inv01 + b * ((1 - a) / det);
        std::mt19937_64 rng(31);
        const std::size_t M = 200000;
        double sp = 0, sm = 0, sp2 = 0, sm2 = 0;
        std::vector<std::pair<std::int64_t, std::uint8_t>> hijas;
        for (std::size_t i = 0; i < M; ++i) {
            dr::generar_cumulo(h, 0, 0, std::numeric_limits<std::int64_t>::max() / 4, rng, hijas);
            double p = 0, m = 0;
            for (const auto& x : hijas) (x.second == 0 ? p : m) += 1.0;
            sp += p; sm += m; sp2 += p * p; sm2 += m * m;
        }
        const double mp = sp / M, mm = sm / M;
        const double ep = std::sqrt((sp2 / M - mp * mp) / M), em = std::sqrt((sm2 / M - mm * mm) / M);
        std::printf("      U3: hijas + %.4f (teoria %.4f, ee %.4f)  hijas - %.4f (teoria %.4f, ee %.4f)\n", mp, e_pp, ep, mm, e_pm, em);
        chk(std::abs(mp - e_pp) <= 4 * ep && std::abs(mm - e_pm) <= 4 * em, "U3", "cumulo de dos tipos = N (I - N)^-1 dentro de 4 ee");
    }

    // R1u y R4: tamano con Poisson homogeneo y con Poisson modulado por fase
    {
        const std::size_t semillas = 40, n = 40000;
        const auto w = mundo(n, 2, 400);
        std::size_t dec_hom = 0, dec_fase = 0, dec_fase_cte = 0;
        for (std::size_t r = 0; r < semillas; ++r) {
            std::mt19937_64 rng(1000 + r);
            AjusteHawkesReplicas hom{}; hom.J = J; hom.tau = tau; hom.bins = 2; hom.mu = {0.004, 0.004, 0.004, 0.004};
            AjusteHawkesReplicas fas{}; fas.J = J; fas.tau = tau; fas.bins = 2; fas.mu = {0.001, 0.008, 0.001, 0.008};
            std::vector<AlarmaReplica> e1, e2;
            dr::simular_historia(hom, false, n, w.pasos_de_bin, 1000000, rng, e1, w.bin_de_paso);
            dr::simular_historia(fas, false, n, w.pasos_de_bin, 1000000, rng, e2, w.bin_de_paso);
            ConfiguracionReplicas cr{}; cr.alfa = alfa; cr.t_max = 64; cr.bins = 2; cr.semilla = 2000 + r;
            dec_hom += ajustar_replicas(e1, w.bin_de_paso, cr).estado == EstadoReplicas::DECIDIBLE;
            dec_fase += ajustar_replicas(e2, w.bin_de_paso, cr).estado == EstadoReplicas::DECIDIBLE;
            ConfiguracionReplicas c1 = cr; c1.bins = 1;
            std::vector<std::uint32_t> uno(n, 0);
            dec_fase_cte += ajustar_replicas(e2, uno, c1).estado == EstadoReplicas::DECIDIBLE;
        }
        const double bd = banda(semillas, alfa);
        std::printf("      R1u/R4: DECIDIBLE sin excitacion: homogeneo %zu/%zu, por fase con base por fase %zu/%zu, "
                    "por fase con base CONSTANTE %zu/%zu (banda %.1f)\n",
                    dec_hom, semillas, dec_fase, semillas, dec_fase_cte, semillas, bd);
        chk(double(dec_hom) <= bd, "R1u", "tamano del contraste del ajuste en un Poisson homogeneo <= alfa + banda");
        chk(double(dec_fase) <= bd, "R4", "Poisson por fase sin excitacion: con base por fase no pasa su nulo (<= alfa + banda)");
        std::printf("      R4 (reporte): con base constante la excitacion falsa pasa su nulo en %zu de %zu\n", dec_fase_cte, semillas);
    }

    // R5: ramificacion critica
    {
        const std::size_t n = 30000;
        const auto w = mundo(n, 1, n);
        auto h = verdad(J, tau, 1, {0.0005, 0.0005}, 1.0, 0.0, 2);
        std::mt19937_64 rng(51);
        // generador propio con desborde ignorado (truncado por el horizonte)
        std::vector<std::pair<std::int64_t, std::uint8_t>> todos;
        std::poisson_distribution<int> pm(h.mu[0] * double(n));
        for (std::uint8_t s = 0; s < 2; ++s) {
            const int k = pm(rng);
            std::uniform_int_distribution<std::int64_t> u(0, std::int64_t(n) - 1);
            for (int c = 0; c < k; ++c) {
                const std::int64_t t0 = u(rng);
                todos.emplace_back(t0, s);
                std::vector<std::pair<std::int64_t, std::uint8_t>> cola{{t0, s}};
                std::size_t cab = 0;
                while (cab < cola.size() && todos.size() < 20000) {
                    auto [tp, mp] = cola[cab++];
                    const double q = std::exp(-1.0 / tau[2]);
                    std::poisson_distribution<int> pk(h.a[mp][mp][2] / (1.0 - q));
                    const int kk = pk(rng);
                    for (int z = 0; z < kk; ++z) {
                        std::geometric_distribution<std::int64_t> g(1.0 - q);
                        const std::int64_t th = tp + 1 + g(rng);
                        if (th >= std::int64_t(n)) continue;
                        cola.emplace_back(th, mp);
                        todos.emplace_back(th, mp);
                    }
                }
            }
        }
        std::sort(todos.begin(), todos.end());
        std::vector<AlarmaReplica> ev;
        for (auto& [t, s] : todos) { AlarmaReplica a{}; a.t = t; a.tipo = s; ev.push_back(a); }
        ConfiguracionReplicas cr{}; cr.alfa = alfa; cr.t_max = 64; cr.bins = 1; cr.semilla = 52;
        const auto f = ajustar_replicas(ev, w.bin_de_paso, cr);
        std::printf("      R5: %zu eventos, rho %.3f, cota %.3f, explosivos %zu, estado %d\n", ev.size(), f.rho, f.rho_sup, f.boot_explosivos,
                    int(f.estado));
        chk(f.estado == EstadoReplicas::NO_IDENTIFICADO, "R5", "ramificacion critica: NO_IDENTIFICADO");
    }

    // R2bu: cruzada conocida y cruzada nula
    {
        const std::size_t n = 200000;
        const auto w = mundo(n, 1, n);
        bool ok = true;
        for (int caso = 0; caso < 2; ++caso) {
            const double nc = caso == 0 ? 0.3 : 0.0;
            auto h = verdad(J, tau, 1, {0.004, 0.004}, 0.3, nc, 3);
            std::mt19937_64 rng(61 + caso);
            std::vector<AlarmaReplica> ev;
            dr::simular_historia(h, true, n, w.pasos_de_bin, 1000000, rng, ev, w.bin_de_paso);
            ConfiguracionReplicas cr{}; cr.alfa = alfa; cr.t_max = 64; cr.bins = 1; cr.semilla = 63 + caso;
            const auto f = ajustar_replicas(ev, w.bin_de_paso, cr);
            std::printf("      R2bu caso %d: N contrario verdad %.3f, ajuste %.3f, basico [%.3f, %.3f]; mismo signo %.3f [%.3f, %.3f]\n",
                        caso, nc, f.n_signo_contrario, f.n_contrario_inf, f.n_contrario_sup, f.n_mismo_signo, f.n_mismo_signo_inf,
                        f.n_mismo_signo_sup);
            ok = ok && nc >= f.n_contrario_inf && nc <= f.n_contrario_sup && 0.3 >= f.n_mismo_signo_inf && 0.3 <= f.n_mismo_signo_sup;
        }
        chk(ok, "R2bu", "excitacion de mismo signo y de signo contrario recuperadas en su intervalo basico (0.3 y 0.0)");
    }

    if (fallos) std::printf("FALLA: %d fallas\n", fallos);
    else std::printf("PASA\n");
    return fallos ? 1 : 0;
}
