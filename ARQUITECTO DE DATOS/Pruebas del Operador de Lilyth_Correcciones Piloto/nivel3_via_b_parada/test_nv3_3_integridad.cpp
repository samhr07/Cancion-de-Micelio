// ============================================================================
// test_nv3_3_integridad.cpp -- Control sellado NV3.3 (T3.7):
//   NV3.3a: S2_N3 (t_max_ef < T_MAX) y n_sim derivado sin literales.
//   NV3.3b: ARFIMA(1, d=+0.35, 0) mu_serie=50->1050 + V4, y NIVEL+Purga + V4.
//   NV3.3c: K10 (l = 320 > 256 evaluado por via (b), 0 reservas en paso())
//           + K1 (identidad bit a bit del operador delta con y sin Lilyth).
//   NV3.3d_cfg: Opt-in Regla 6 por omision == false y arena_cap == 0.
// Semilla ciega exclusiva: SEMILLA_CIEGA_NIVEL3 = 20261124ULL.
// ============================================================================

#include <atomic>
#include <cstdlib>
#include <new>

namespace {
std::atomic<bool> g_contar_new{false};
std::atomic<std::size_t> g_news{0};
}  // namespace

void* operator new(std::size_t sz) {
    if (g_contar_new.load(std::memory_order_relaxed)) {
        g_news.fetch_add(1, std::memory_order_relaxed);
    }
    void* p = std::malloc(sz == 0 ? 1 : sz);
    if (!p) throw std::bad_alloc();
    return p;
}

void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <memory>
#include <random>
#include <span>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;

constexpr std::uint64_t SEMILLA_CIEGA_NIVEL3 = 20261124ULL;
constexpr std::uint64_t SEM_3A = SEMILLA_CIEGA_NIVEL3 + 1009ULL;
constexpr std::uint64_t SEM_3B = SEMILLA_CIEGA_NIVEL3 + 3001ULL;
constexpr std::uint64_t SEM_3C = SEMILLA_CIEGA_NIVEL3 + 7001ULL;

bool bit_igual_double(double a, double b) noexcept {
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

bool alarma_bit_igual(const EventoAlarma& a, const EventoAlarma& b) noexcept {
    return a.hay == b.hay &&
           a.signo == b.signo &&
           a.t == b.t &&
           a.t_cambio == b.t_cambio &&
           a.l == b.l &&
           bit_igual_double(a.g, b.g) &&
           bit_igual_double(a.amplitud, b.amplitud) &&
           bit_igual_double(a.se_amplitud, b.se_amplitud);
}

std::vector<double> generar_ar1_unidad(std::size_t n, double phi, double mu, std::uint64_t sem) {
    const std::size_t burn = 512;
    const double sd_e = std::sqrt(std::max(1e-12, 1.0 - phi * phi));
    std::mt19937_64 rng(detalle_nulo::splitmix64(sem));
    std::normal_distribution<double> nd(0.0, sd_e);
    std::vector<double> out(n, 0.0);
    double u = 0.0;
    for (std::size_t t = 0; t < burn + n; ++t) {
        u = phi * u + nd(rng);
        if (t >= burn) out[t - burn] = mu + u;
    }
    return out;
}

std::vector<double> generar_arfima1(
    std::size_t n, double phi1, double d, double mu_serie, double sigma_e, std::uint64_t sem) {
    const std::size_t burn = 1024;
    const std::size_t tot = burn + n;
    std::vector<double> pi_inv(tot, 0.0);
    elw::coeficientes(-d, tot, pi_inv.data());
    std::vector<double> u(tot, 0.0);
    std::vector<double> x(n, 0.0);
    std::mt19937_64 rng(detalle_nulo::splitmix64(sem));
    std::normal_distribution<double> nd(0.0, sigma_e);
    double u_prev = 0.0;
    for (std::size_t s = 0; s < tot; ++s) {
        const double u_cur = phi1 * u_prev + nd(rng);
        u[s] = u_cur;
        u_prev = u_cur;
        if (s >= burn) {
            double suma = 0.0;
            for (std::size_t k = 0; k <= s; ++k) suma += pi_inv[k] * u[s - k];
            x[s - burn] = mu_serie + suma;
        }
    }
    return x;
}

}  // namespace

int main() {
    std::size_t fallas = 0;
    auto check = [&](const char* tag, bool cond, const char* msg) {
        std::printf("  [%s] %-10s : %s\n", cond ? "PASA" : "FALLA", tag, msg);
        if (!cond) ++fallas;
    };

    std::printf("================================================================================\n");
    std::printf("NV3.3 (T3.7) -- Integridad Estructural (S2_N3, ARFIMA+V4, K10 l>256, K1, R6)\n");
    std::printf("  SEMILLA_CIEGA_NIVEL3 = %lluULL\n", static_cast<unsigned long long>(SEMILLA_CIEGA_NIVEL3));
    std::printf("================================================================================\n\n");

    // -------------------------------------------------------------------------
    // NV3.3a: S2_N3 (ventana t_max_ef < T_MAX) y n_sim derivado sin literales
    // -------------------------------------------------------------------------
    {
        const auto d_fl = dimensiones_simulacion_via_b(120, 14, 5, 22, Integracion::FLUJO);
        const auto d_nv = dimensiones_simulacion_via_b(120, 14, 5, 22, Integracion::NIVEL);
        const auto d_0  = dimensiones_simulacion_via_b(0, 8, 3, 10, Integracion::NIVEL);

        const bool ok_dim =
            d_fl.w_pre == 134 && d_fl.w_sorteo == 120 && d_fl.w_post == 27 &&
            d_fl.n_sim == 281 && d_fl.t_iny_min == 134 && d_fl.t_iny_max == 253 &&
            static_cast<std::size_t>(d_fl.t_iny_max) + 5u + 22u == d_fl.n_sim - 1u &&
            d_nv.w_pre == 135 && d_nv.w_sorteo == 120 && d_nv.w_post == 27 &&
            d_nv.n_sim == 282 && d_nv.t_iny_min == 135 && d_nv.t_iny_max == 254 &&
            static_cast<std::size_t>(d_nv.t_iny_max) + 5u + 22u == d_nv.n_sim - 1u &&
            d_0.w_pre == 9 && d_0.w_sorteo == 1 && d_0.w_post == 13 &&
            d_0.n_sim == 23 && d_0.t_iny_min == 9 && d_0.t_iny_max == 9 &&
            static_cast<std::size_t>(d_0.t_iny_max) + 3u + 10u == d_0.n_sim - 1u;

        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 128>;
        const auto x_cal = generar_ar1_unidad(4096, 0.85, 0.0, SEM_3A + 1ULL);

        ConfiguracionObservador cd_ref{};
        cd_ref.integracion = Integracion::FLUJO;
        cd_ref.agregacion = AgregacionEscala::SUMA;
        cd_ref.arl0_pasos_nativos = 450.0;
        cd_ref.medir_tau50 = false;
        cd_ref.escalas_max = 1;
        cd_ref.semilla = SEM_3A + 1ULL;
        auto od_ref = std::make_unique<Obs>();
        od_ref->calibrar(cd_ref, x_cal);

        ConfiguracionLilyth cl_ref{};
        cl_ref.representacion = ModoRepresentacionLilyth::CRUDA;
        cl_ref.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl_ref.b_max_nulo = 8;
        cl_ref.hilos_nulo = 1;
        cl_ref.semilla_nulo = SEM_3A + 2ULL;
        cl_ref.usar_enmienda_nivel2 = true;
        cl_ref.usar_enmienda_nivel3 = true;

        auto lil_ref = std::make_unique<Lil>();
        const bool cal_ref = lil_ref->calibrar(*od_ref, cl_ref, cd_ref, x_cal);

        ConfiguracionLilyth cl_base = cl_ref;
        cl_base.usar_enmienda_nivel3 = false;
        auto lil_base = std::make_unique<Lil>();
        const bool cal_base = lil_base->calibrar(*od_ref, cl_base, cd_ref, x_cal);

        const std::size_t t_max_ef = lil_ref->configuracion().t_max_efectivo;
        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od_ref->ficha_nivel(0).nulo.tau_int)));
        const std::size_t l = 4;
        const std::size_t dist_borde = w_ancla + guarda + l + t_max_ef;
        const std::size_t dist_entre = w_ancla + guarda + l + t_max_ef + 16u;
        const std::size_t dist_tmax = w_ancla + guarda + l + 128u;
        const bool ok_horiz = (t_max_ef + 16u < 128u) && (dist_entre < dist_tmax);

        ConfiguracionObservador cd_hi = cd_ref;
        cd_hi.arl0_pasos_nativos = 1.0e6;
        auto od_1e6 = std::make_unique<Obs>();
        od_1e6->calibrar(cd_hi, x_cal);

        ConfiguracionLilyth cl_1e6 = cl_ref;
        cl_1e6.t_max_efectivo = t_max_ef;
        auto lil_1e6 = std::make_unique<Lil>();
        const bool cfg_hi = lil_1e6->configurar(cl_1e6, cd_hi, *od_1e6, x_cal);

        const auto& fn = lil_ref->ficha_nulo_cruda();
        const std::vector<double> hist_prev(fn.modelo_lineal.p + 8u, 0.0);
        const std::vector<double> perfil{4.0, 6.0, 7.5, 8.5};
        std::vector<double> e_iny(l, 0.0);
        innovaciones(fn.modelo_lineal, hist_prev.data(), hist_prev.size(), perfil.data(), l, e_iny.data());

        const std::array<std::int64_t, 2> t_borde{250, static_cast<std::int64_t>(250u + dist_borde)};
        const std::array<std::int64_t, 2> t_entre{250, static_cast<std::int64_t>(250u + dist_entre)};
        const std::size_t n_sim = 1200;

        auto correr_s2 = [&](std::span<const std::int64_t> ts, bool usar_n3, std::uint64_t sem) {
            std::mt19937_64 rng(detalle_nulo::splitmix64(sem));
            return detalle_nulo::ejecutar_inyecciones_en_realizacion<128, RasgosDelta<Obs>::b_max>(
                *od_1e6, *lil_1e6, Integracion::FLUJO, false, false, cl_1e6.nulo_sin_memoria,
                fn.segmentos_calma, fn.modelo_lineal, lil_1e6->contexto_purga_nulo(), NaN,
                perfil, e_iny, +1, w_ancla, guarda, t_max_ef, ts, n_sim, rng, true, usar_n3);
        };

        const auto ra_borde = correr_s2(t_borde, true,  SEM_3A + 10ULL);
        const auto ra_ent_3 = correr_s2(t_entre, true,  SEM_3A + 11ULL);
        const auto ra_ent_0 = correr_s2(t_entre, false, SEM_3A + 11ULL);

        bool hallada = false;
        ResultadoInyeccionesRealizacion<128> rb_3{}, rb_0{};
        for (std::uint64_t s = 0; s < 32; ++s) {
            const std::uint64_t sem = SEM_3A + 100ULL + s;
            const auto probe = lil_ref->evaluar_inyecciones_realizacion(
                false, hist_prev, perfil, +1, w_ancla, guarda, t_entre, n_sim, sem, false, true);
            if (probe.con_ajena_previa == 0 && probe.solapes == 0 && probe.eventos_admitidos.size() == 2) {
                hallada = true;
                rb_3 = probe;
                rb_0 = lil_base->evaluar_inyecciones_realizacion(
                    false, hist_prev, perfil, +1, w_ancla, guarda, t_entre, n_sim, sem, false, true);
                break;
            }
        }

        const bool ok_3a = ok_dim && cal_ref && cal_base && cfg_hi && ok_horiz &&
                           ra_borde.solapes == 2 && ra_borde.eventos_admitidos.empty() &&
                           ra_ent_3.solapes == 0 && ra_ent_3.eventos_admitidos.size() == 2 &&
                           ra_ent_0.solapes == 2 && ra_ent_0.eventos_admitidos.empty() &&
                           hallada && rb_3.solapes == 0 && rb_3.eventos_admitidos.size() == 2 &&
                           rb_0.solapes == 2 && rb_0.eventos_admitidos.empty();

        std::printf("[NV3.3a] t_max_ef=%zu (<T_MAX=128) n_sim(fl=%zu,nv=%zu) | borde(sol=%zu adm=%zu) | entre n3=1(sol=%zu adm=%zu) n3=0(sol=%zu adm=%zu) | pass(n3=1 sol=%zu adm=%zu, n3=0 sol=%zu)\n",
                    t_max_ef, d_fl.n_sim, d_nv.n_sim,
                    ra_borde.solapes, ra_borde.eventos_admitidos.size(),
                    ra_ent_3.solapes, ra_ent_3.eventos_admitidos.size(),
                    ra_ent_0.solapes, ra_ent_0.eventos_admitidos.size(),
                    rb_3.solapes, rb_3.eventos_admitidos.size(), rb_0.solapes);
        check("NV3.3a", ok_3a, "S2_N3 exacto en borde e intermedio, n_sim derivado en FLUJO/NIVEL y passthrough");
    }

    // -------------------------------------------------------------------------
    // NV3.3b: ARFIMA(1, d=+0.35, 0) mu_serie=50->1050 + V4, y NIVEL+Purga + V4
    // -------------------------------------------------------------------------
    {
        using Obs = ObservadorSaltos<24, 64, 120>;
        using Lil = ObservadorLilyth<Obs, 8, 64>;

        const auto z50 = generar_arfima1(4096, 0.40, 0.35, 50.0, 1.0, SEM_3B + 1ULL);
        std::vector<double> z1050 = z50;
        for (double& v : z1050) v += 1000.0;
        const auto mom50 = est::momentos(std::vector<std::vector<double>>{z50});

        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.agregacion = AgregacionEscala::SUMA;
        cd.arl0_pasos_nativos = 1.0e6;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        cd.semilla = SEM_3B + 1ULL;
        auto od50 = std::make_unique<Obs>();
        auto od1050 = std::make_unique<Obs>();
        od50->calibrar(cd, z50);
        od1050->calibrar(cd, z1050);

        ConfiguracionLilyth cl{};
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.modo_nulo_lineal = ModoNuloLineal::ARFIMA_WHITTLE;
        cl.d_whittle_cruda = 0.35;
        cl.b_max_nulo = 2;
        cl.hilos_nulo = 1;
        cl.semilla_nulo = SEM_3B + 2ULL;
        cl.usar_enmienda_nivel2 = true;
        cl.usar_enmienda_nivel3 = true;

        auto lil50 = std::make_unique<Lil>();
        auto lil1050 = std::make_unique<Lil>();
        const bool cal50 = lil50->calibrar(*od50, cl, cd, z50);
        const bool cal1050 = lil1050->calibrar(*od1050, cl, cd, z1050);

        const std::size_t nh = std::min<std::size_t>(lil50->ficha_nulo_cruda().modelo_lineal.n_hist_req, 1500u);
        const std::vector<double> hp50(z50.begin() + static_cast<std::ptrdiff_t>(2000u - nh), z50.begin() + 2000);
        const std::vector<double> hp1050(z1050.begin() + static_cast<std::ptrdiff_t>(2000u - nh), z1050.begin() + 2000);

        constexpr std::size_t L = 12;
        std::vector<double> pf50(L, 0.0), pf1050(L, 0.0);
        for (std::size_t k = 0; k < L; ++k) {
            const double bump = (k < 4 ? (6.0 + 3.0 * static_cast<double>(k)) : 15.0) * mom50.sd;
            pf50[k] = z50[2000u + k] + bump;
            pf1050[k] = pf50[k] + 1000.0;
        }

        const std::size_t w_ancla = 120;
        const std::size_t guarda = static_cast<std::size_t>(std::ceil(std::max(1.0, od50->ficha_nivel(0).nulo.tau_int)));
        const std::array<std::int64_t, 1> t_arr{450};
        const std::size_t n_sim = 900;

        const auto r50 = lil50->evaluar_inyecciones_realizacion(
            false, hp50, pf50, +1, w_ancla, guarda, t_arr, n_sim, SEM_3B + 3ULL, false, true, 0);
        const auto r1050 = lil1050->evaluar_inyecciones_realizacion(
            false, hp1050, pf1050, +1, w_ancla, guarda, t_arr, n_sim, SEM_3B + 3ULL, false, true, 0);

        bool ok_arf = cal50 && cal1050 &&
                      r50.paradas_tempranas > 0 && r1050.paradas_tempranas > 0 &&
                      r50.detectadas == 1 && r1050.detectadas == 1 &&
                      r50.eventos_admitidos.size() == 1 && r1050.eventos_admitidos.size() == 1;
        double dif_anc = NaN, dif_dy = NaN, max_dif_eta = NaN, max_dif_xsim = NaN;
        if (ok_arf) {
            const auto& e50 = r50.eventos_admitidos.front();
            const auto& e1050 = r1050.eventos_admitidos.front();
            dif_anc = std::abs(e1050.y_ancla - (e50.y_ancla + 1000.0));
            dif_dy = std::abs(e1050.dy - e50.dy);
            max_dif_eta = 0.0;
            const double x0_50 = e50.y_ancla + e50.dy;
            const double x0_1050 = e1050.y_ancla + e1050.dy;
            max_dif_xsim = std::abs(x0_1050 - (x0_50 + 1000.0));
            ok_arf = ok_arf && (e50.pasos_observados == e1050.pasos_observados) && (e50.tau_fin == e1050.tau_fin);
            for (std::size_t tau = 0; tau < e50.pasos_observados && tau < e50.eta.size(); ++tau) {
                max_dif_eta = std::max(max_dif_eta, std::abs(e1050.eta[tau] - e50.eta[tau]));
                const double xt50 = e50.y_ancla + e50.dy * e50.eta[tau];
                const double xt1050 = e1050.y_ancla + e1050.dy * e1050.eta[tau];
                max_dif_xsim = std::max(max_dif_xsim, std::abs(xt1050 - (xt50 + 1000.0)));
            }
            const bool ambos_nan = std::isnan(e50.tau_rest_consolidado) && std::isnan(e1050.tau_rest_consolidado);
            ok_arf = ok_arf && (dif_anc <= 1e-10) && (dif_dy <= 1e-10) &&
                     (max_dif_eta <= 1e-10) && (max_dif_xsim <= 1e-10) &&
                     (ambos_nan || e50.tau_rest_consolidado == e1050.tau_rest_consolidado);
        }

        // NIVEL + PURGADA + V4
        constexpr std::size_t BINS = 4;
        using ObsP = ObservadorSaltos<24, 64, 120, 8>;
        using LilP = ObservadorLilyth<ObsP, 8, 64>;

        const auto inc_cal = generar_ar1_unidad(4096, 0.70, 0.0, SEM_3B + 101ULL);
        const std::array<double, BINS> escala_f{0.8, 1.3, 1.0, 1.6};
        std::vector<double> x_cal(4096, 0.0);
        std::vector<std::size_t> f_cal(4096, 0);
        std::vector<std::int64_t> d_cal(4096, 0);
        double lvl = 100.0;
        for (std::size_t t = 0; t < 4096; ++t) {
            f_cal[t] = t % BINS;
            d_cal[t] = static_cast<std::int64_t>(t / BINS);
            if (t > 0) lvl += inc_cal[t] * escala_f[f_cal[t]];
            x_cal[t] = lvl;
        }

        ConfiguracionObservador cd_np{};
        cd_np.integracion = Integracion::NIVEL;
        cd_np.purga = true;
        cd_np.bins_fase = BINS;
        cd_np.dias_minimos_purga = 5;
        cd_np.arl0_pasos_nativos = 450.0;
        cd_np.medir_tau50 = false;
        cd_np.escalas_max = 1;
        cd_np.semilla = SEM_3B + 102ULL;
        auto od_np = std::make_unique<ObsP>();
        od_np->calibrar(cd_np, x_cal, {}, {}, f_cal, d_cal);

        ConfiguracionLilyth cl_np{};
        cl_np.representacion = ModoRepresentacionLilyth::PURGADA;
        cl_np.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl_np.b_max_nulo = 2;
        cl_np.hilos_nulo = 1;
        cl_np.semilla_nulo = SEM_3B + 103ULL;
        cl_np.usar_enmienda_nivel2 = true;
        cl_np.usar_enmienda_nivel3 = true;

        auto lil_np = std::make_unique<LilP>();
        const bool cal_np = lil_np->calibrar(*od_np, cl_np, cd_np, x_cal, {}, f_cal, d_cal);

        const auto est = generar_ar1_unidad(512, 0.70, 0.0, SEM_3B + 104ULL);
        std::vector<double> x_lac;
        std::vector<std::size_t> f_lac;
        std::vector<std::int64_t> d_lac;
        std::size_t sf_lac = 0;
        detalle_nulo::llevar_a_crudo<8>(
            Integracion::NIVEL, true, lil_np->contexto_purga_nulo(), est, 100.0, 3, 0, x_lac, f_lac, d_lac, sf_lac);

        detalle_nulo::PasoCrudoNulo<8> pc(Integracion::NIVEL, true, lil_np->contexto_purga_nulo(), 100.0, 3, 0);
        std::size_t sf_pc = 0;
        double max_dif_pc = 0.0;
        bool fd_ok = true;
        for (std::size_t t = 0; t < est.size(); ++t) {
            std::size_t ft = 0;
            std::int64_t dt = 0;
            const double xp = pc.paso(t, est[t], ft, dt, sf_pc);
            if (ft != f_lac[t] || dt != d_lac[t]) fd_ok = false;
            const double d = std::abs(xp - x_lac[t]);
            if (!std::isfinite(d) || d > max_dif_pc) max_dif_pc = d;
        }

        const std::vector<double> hp_np(lil_np->ficha_nulo_purgada().modelo_lineal.p + 8u, 0.0);
        const std::vector<double> pf_np{5.0, 7.0, 9.0, 11.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0};
        const std::size_t guarda_np = static_cast<std::size_t>(std::ceil(std::max(1.0, od_np->ficha_nivel(0).nulo.tau_int)));
        const std::array<std::int64_t, 1> t_np{450};
        const auto r_np = lil_np->evaluar_inyecciones_realizacion(
            true, hp_np, pf_np, +1, w_ancla, guarda_np, t_np, 900, SEM_3B + 105ULL, false, true, 0);

        bool fin_np = cal_np && fd_ok && (sf_pc == sf_lac) && (max_dif_pc == 0.0) &&
                      (r_np.paradas_tempranas > 0) && (r_np.detectadas == 1) &&
                      (r_np.sin_factor_purga == 0) && (r_np.eventos_admitidos.size() == 1);
        if (fin_np) {
            const auto& ev = r_np.eventos_admitidos.front();
            fin_np = std::isfinite(ev.y_ancla) && std::isfinite(ev.dy) && ev.pasos_observados > 0;
            for (std::size_t tau = 0; tau < ev.pasos_observados && fin_np; ++tau) {
                fin_np = std::isfinite(ev.eta[tau]);
            }
        }

        std::printf("[NV3.3b] ARFIMA+V4(par=%zu,%zu dif_anc=%.2e dif_dy=%.2e max|dif_eta|=%.2e max|dif_xsim|=%.2e) | NIVEL+Purga(dif_pc=%.1e par=%zu sf=%zu)\n",
                    r50.paradas_tempranas, r1050.paradas_tempranas, dif_anc, dif_dy, max_dif_eta, max_dif_xsim,
                    max_dif_pc, r_np.paradas_tempranas, r_np.sin_factor_purga);
        check("NV3.3b", ok_arf && fin_np, "ARFIMA mu_serie + V4 invariante <= 1e-10 y coherencia NIVEL+PURGADA+V4");
    }

    // -------------------------------------------------------------------------
    // NV3.3c: K1 bit a bit + K10 (l=320 > 256 evaluado por via (b))
    // -------------------------------------------------------------------------
    {
        using ObsS3 = ObservadorSaltos<24, 120, 120>;
        using ObsL3 = ObservadorLilyth<ObsS3, 8, 64>;

        const auto x_cal = generar_ar1_unidad(10000, 0.70, 10.0, SEM_3C + 10ULL);
        ConfiguracionObservador cd3{};
        cd3.integracion = Integracion::FLUJO;
        cd3.agregacion = AgregacionEscala::SUMA;
        cd3.arl0_pasos_nativos = 450.0;
        cd3.medir_tau50 = false;
        cd3.escalas_max = 3;
        cd3.semilla = SEM_3C + 1ULL;

        auto od_solo = std::make_unique<ObsS3>();
        auto od_con_lil = std::make_unique<ObsS3>();
        od_solo->calibrar(cd3, x_cal);
        od_con_lil->calibrar(cd3, x_cal);

        ConfiguracionLilyth cl3{};
        cl3.representacion = ModoRepresentacionLilyth::CRUDA;
        cl3.modo_nulo_lineal = ModoNuloLineal::AR_LEVINSON_DURBIN;
        cl3.b_max_nulo = 8;
        cl3.hilos_nulo = 1;
        cl3.semilla_nulo = SEM_3C + 2ULL;
        cl3.forzar_via_b = true;
        cl3.usar_enmienda_nivel2 = true;
        cl3.usar_enmienda_nivel3 = true;

        auto lil3 = std::make_unique<ObsL3>();
        const bool cal3 = lil3->calibrar(*od_con_lil, cl3, cd3, x_cal);

        // 1. K1 bit a bit
        auto x_k1 = generar_ar1_unidad(2000, 0.70, 10.0, SEM_3C + 20ULL);
        const auto mom_k1 = est::momentos(std::vector<std::vector<double>>{x_k1});
        x_k1[700] += 6.0 * mom_k1.sd;
        x_k1[1300] -= 6.0 * mom_k1.sd;

        od_solo->reiniciar_flujo();
        od_con_lil->reiniciar_flujo();
        lil3->reiniciar_flujo_conservando_purga();

        std::size_t diffs_k1 = 0;
        g_news.store(0, std::memory_order_relaxed);
        g_contar_new.store(true, std::memory_order_release);
        for (std::size_t t = 0; t < x_k1.size(); ++t) {
            const auto s1 = od_solo->paso(x_k1[t], NaN, true, 0, 0);
            const auto s2 = od_con_lil->paso(x_k1[t], NaN, true, 0, 0);
            (void)lil3->paso(s2, x_k1[t], true, NaN, 0, 0);

            bool eq = (s1.t == s2.t) &&
                      (s1.n_eventos == s2.n_eventos) &&
                      (s1.evento_abierto == s2.evento_abierto) &&
                      (s1.niveles == s2.niveles);
            for (std::size_t j = 0; eq && j < s1.niveles; ++j) {
                eq = bit_igual_double(s1.nivel[j].g, s2.nivel[j].g) &&
                     alarma_bit_igual(s1.nivel[j].alarma, s2.nivel[j].alarma);
            }
            for (std::size_t e = 0; eq && e < s1.n_eventos; ++e) {
                eq = (s1.eventos[e].canal == s2.eventos[e].canal) &&
                     (s1.eventos[e].escala_origen == s2.eventos[e].escala_origen) &&
                     alarma_bit_igual(s1.eventos[e].alarma, s2.eventos[e].alarma);
            }
            if (!eq) ++diffs_k1;
        }
        g_contar_new.store(false, std::memory_order_release);
        const std::size_t news_k1 = g_news.load(std::memory_order_relaxed);
        const std::size_t ap_k1 = lil3->contadores().aperturas;
        const std::size_t esp_inc_k1 = lil3->contadores().espejo_inconsistente;

        // 2. K10 con l = 320 > 256 evaluado por via (b)
        auto x_k10 = generar_ar1_unidad(1400, 0.70, 10.0, SEM_3C + 30ULL);
        const auto mom_k10 = est::momentos(std::vector<std::vector<double>>{x_k10});
        for (std::size_t t = 800; t <= 1119; ++t) {
            const std::size_t k = t - 800u;
            const double bump = (k < 300u)
                ? ((1.0 + 2.0 * static_cast<double>(k) / 300.0) * mom_k10.sd)
                : ((3.0 + 6.0 * static_cast<double>(k - 300u) / 19.0) * mom_k10.sd);
            x_k10[t] += bump;
        }

        od_con_lil->reiniciar_flujo();
        lil3->reiniciar_flujo_conservando_purga();

        g_news.store(0, std::memory_order_relaxed);
        g_contar_new.store(true, std::memory_order_release);
        for (std::size_t t = 0; t < x_k10.size(); ++t) {
            auto s = od_con_lil->paso(x_k10[t], NaN, true, 0, 0);
            for (std::size_t j = 0; j < s.niveles; ++j) s.nivel[j].alarma.hay = false;
            if (t == 1119) {
                s.nivel[2].alarma.hay = true;
                s.nivel[2].alarma.signo = 1;
                s.nivel[2].alarma.t = 1119;
                s.nivel[2].alarma.t_cambio = 800;
                s.nivel[2].alarma.l = 80;
                s.nivel[2].alarma.g = 20.0;
                s.nivel[2].alarma.amplitud = 9.0 * mom_k10.sd;
                s.nivel[2].alarma.se_amplitud = mom_k10.sd;
            }
            (void)lil3->paso(s, x_k10[t], true, NaN, 0, 0);
        }
        g_contar_new.store(false, std::memory_order_release);
        const std::size_t news_k10 = g_news.load(std::memory_order_relaxed);

        lil3->completar_nulos();
        const auto& fn3 = lil3->ficha_nulo_cruda();
        bool enc_1119 = false;
        ObsL3::Ranura r1119{};
        for (const auto& rc : lil3->cerradas_historico()) {
            if (rc.t_apertura == 1119) {
                r1119 = rc;
                enc_1119 = true;
                break;
            }
        }
        const auto pf3 = enc_1119 ? lil3->perfil_evento_n2(r1119, false) : std::span<const double>{};
        const auto ei3 = enc_1119 ? lil3->e_iny_evento_n2(r1119, false) : std::span<const double>{};
        const auto& ref1119 = r1119.cruda.nulo_onsager;

        const bool des_ok = enc_1119 &&
                            (r1119.cruda.desenlace == DesenlaceRanura::RESTITUIDA ||
                             r1119.cruda.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE);
        const bool ok_3c = cal3 &&
                           od_con_lil->escalas() == 3 &&
                           lil3->cap_perfil_adaptable() == 480 &&
                           diffs_k1 == 0 && news_k1 == 0 && esp_inc_k1 == 0 && ap_k1 > 0 &&
                           news_k10 == 0 &&
                           lil3->contadores().arena_perfiles_desbordada == 0 &&
                           sizeof(ObsL3::Ranura) < 16384 &&
                           des_ok &&
                           pf3.size() == 320 && pf3.size() > 256 &&
                           ei3.size() == 320 &&
                           ref1119.nulo_por_inyeccion &&
                           ref1119.n_nulos_b > 0 &&
                           ref1119.n_nulos_b >= fn3.m_min &&
                           ref1119.solapes_b == 0;

        std::printf("[NV3.3c] K1(diffs=%zu news=%zu ap=%zu esp_inc=%zu) | K10(news=%zu des=%d pf=%zu ei=%zu nulB=%zu/%zu parB=%zu solB=%zu sizeof=%zu)\n",
                    diffs_k1, news_k1, ap_k1, esp_inc_k1,
                    news_k10, enc_1119 ? static_cast<int>(r1119.cruda.desenlace) : -1,
                    pf3.size(), ei3.size(), ref1119.n_nulos_b, fn3.m_min,
                    ref1119.paradas_tempranas_b, ref1119.solapes_b, sizeof(ObsL3::Ranura));
        check("NV3.3c", ok_3c, "K1 bit a bit + K10 l=320 > 256 por via (b) con reservas_en_paso==0");
    }

    // -------------------------------------------------------------------------
    // NV3.3d_cfg: Regla 6 por omision == false y arena_cap == 0
    // -------------------------------------------------------------------------
    {
        ConfiguracionModeloLineal cm_def{};
        ConfiguracionLilyth cl_def{};
        cl_def.t_max_efectivo = 32;

        using ObsS = ObservadorSaltos<24, 64, 60>;
        using ObsL = ObservadorLilyth<ObsS, 8, 64>;
        const auto x_def = generar_ar1_unidad(2000, 0.40, 0.0, SEMILLA_CIEGA_NIVEL3 + 9500ULL);
        ConfiguracionObservador cd_def{};
        cd_def.integracion = Integracion::FLUJO;
        cd_def.agregacion = AgregacionEscala::SUMA;
        cd_def.arl0_pasos_nativos = 450.0;
        cd_def.medir_tau50 = false;
        cd_def.escalas_max = 1;
        auto od_def = std::make_unique<ObsS>();
        od_def->calibrar(cd_def, x_def);
        auto lil_def = std::make_unique<ObsL>();
        const bool cfg_ok = lil_def->configurar(cl_def, cd_def, *od_def, x_def);

        const bool ok_cfg = !cm_def.usar_enmienda_nivel2 &&
                            !cl_def.usar_enmienda_nivel2 &&
                            !cl_def.usar_enmienda_nivel3 &&
                            cfg_ok &&
                            lil_def->arena_perfiles_capacidad() == 0;
        check("NV3.3d_cfg", ok_cfg, "usar_enmienda_nivel2/3 == false por omision y arena_cap == 0");
    }

    std::printf("\n--------------------------------------------------------------------------------\n");
    std::printf("FALLAS NV3.3: %zu\n", fallas);
    std::fflush(stdout);
    return fallas == 0 ? 0 : 1;
}
