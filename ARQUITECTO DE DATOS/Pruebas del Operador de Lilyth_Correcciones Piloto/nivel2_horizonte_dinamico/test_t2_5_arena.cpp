// ============================================================================
// test_t2_5_arena.cpp -- Micro-test de verificacion de T2.5:
//   1. Adaptabilidad de cap_perfil_adaptable() segun el rango calibrado.
//   2. Pre-reserva opt-in en arena_perfiles_n2_ (0 cuando usar_enmienda_nivel2=false)
//      y cero reservas dinamicas en paso() (contador global operator new = 0).
//   3. Cero inflado de sizeof(RanuraLilyth) (< 16384 bytes) y cero copia de
//      arena al clonar lil_local tras reiniciar_flujo_conservando_purga().
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

void operator delete(void* p) noexcept {
    std::free(p);
}

void operator delete(void* p, std::size_t) noexcept {
    std::free(p);
}

#include "ark/lilyth/observador_lilyth.hpp"
#include "ark/metrologia/observador_saltos.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

namespace {

using namespace ark::metrologia;
using namespace ark::lilyth;

std::vector<double> generar_ar1(std::size_t n, double phi, std::uint64_t sem) {
    std::mt19937_64 rng(sem);
    std::normal_distribution<double> n01(0.0, 1.0);
    std::vector<double> x(n, 0.0);
    double cur = 0.0;
    for (std::size_t t = 0; t < n; ++t) {
        cur = phi * cur + n01(rng);
        x[t] = cur;
    }
    return x;
}

}  // namespace

int main() {
    std::size_t fallas = 0;
    auto check = [&](bool cond, const char* msg) {
        if (!cond) {
            ++fallas;
            std::printf("  [FALLA] %s\n", msg);
        }
    };

    std::printf("================================================================================\n");
    std::printf(" TEST T2.5: MEMORIA ADAPTABLE Y ARENA CERO-HEAP EN ObservadorLilyth\n");
    std::printf("================================================================================\n\n");

    // Check 1a: L_GLR = 64, escalas_max = 1 -> cap_perfil_adaptable() == 64
    {
        using ObsS = ObservadorSaltos<24, 64, 120>;
        using ObsL = ObservadorLilyth<ObsS, 8, 64>;
        const auto x_cal = generar_ar1(4096, 0.70, 2026112301ULL);
        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 450.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 1;
        auto od = std::make_unique<ObsS>();
        od->calibrar(cd, x_cal);

        ConfiguracionLilyth cl{};
        cl.t_max_efectivo = 32;
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.usar_enmienda_nivel2 = false;
        auto lil_off = std::make_unique<ObsL>();
        lil_off->configurar(cl, cd, *od, x_cal);

        std::printf("  Caso 1a (L=64, J=1, N2=false): cap=%zu, arena_cap=%zu, sizeof(Ranura)=%zu\n",
                    lil_off->cap_perfil_adaptable(), lil_off->arena_perfiles_capacidad(), sizeof(ObsL::Ranura));
        check(lil_off->cap_perfil_adaptable() == 64, "L=64 J=1 -> cap == 64");
        check(lil_off->arena_perfiles_capacidad() == 0, "usar_enmienda_nivel2=false -> arena_cap == 0");
        check(sizeof(ObsL::Ranura) < 16384, "sizeof(Ranura) libre de arreglos HIST_MAX=4096");
    }

    // Check 1b + 2 + 3: L_GLR = 120, escalas_max = 3 (J=3, m_max=4) -> cap == 480
    {
        using ObsS = ObservadorSaltos<24, 120, 120>;
        using ObsL = ObservadorLilyth<ObsS, 8, 64>;
        const auto x_cal = generar_ar1(10000, 0.70, 2026112302ULL);
        ConfiguracionObservador cd{};
        cd.integracion = Integracion::FLUJO;
        cd.arl0_pasos_nativos = 450.0;
        cd.medir_tau50 = false;
        cd.escalas_max = 3;
        auto od = std::make_unique<ObsS>();
        od->calibrar(cd, x_cal);

        ConfiguracionLilyth cl{};
        cl.t_max_efectivo = 32;
        cl.representacion = ModoRepresentacionLilyth::CRUDA;
        cl.usar_enmienda_nivel2 = true;
        auto lil = std::make_unique<ObsL>();
        lil->configurar(cl, cd, *od, x_cal);

        std::printf("  Caso 1b (L=120, J=%zu, N2=true): cap=%zu, arena_cap=%zu (min_esp=%zu), sizeof(Ranura)=%zu\n",
                    od->escalas(), lil->cap_perfil_adaptable(), lil->arena_perfiles_capacidad(),
                    std::size_t{1024} * 2 * lil->cap_perfil_adaptable(), sizeof(ObsL::Ranura));

        check(od->escalas() == 3, "Escalas calibradas == 3 con 10000 barras");
        check(lil->cap_perfil_adaptable() == 480, "L=120 J=3 -> cap_perfil_adaptable == 480 (> 256)");
        check(lil->arena_perfiles_capacidad() >= std::size_t{1024} * 2 * 480, "Arena pre-reservada >= 1024*2*480");
        check(sizeof(ObsL::Ranura) < 16384, "sizeof(Ranura) < 16384");

        // Cero reservas en paso()
        const auto x_op = generar_ar1(2000, 0.70, 2026112303ULL);
        od->reiniciar_flujo();
        g_news.store(0, std::memory_order_relaxed);
        g_contar_new.store(true, std::memory_order_release);
        for (std::size_t t = 0; t < x_op.size(); ++t) {
            const auto s = od->paso(x_op[t], NaN, true, 0, 0);
            (void)lil->paso(s, x_op[t], true, NaN, 0, 0);
        }
        g_contar_new.store(false, std::memory_order_release);
        const std::size_t news_paso = g_news.load(std::memory_order_relaxed);
        std::printf("  K10 durante 2000 pasos caliente: operator new = %zu\n", news_paso);
        check(news_paso == 0, "Cero reservas dinamicas (new == 0) en paso()");

        // Cero copia de arena en lil_local tras reiniciar_flujo_conservando_purga()
        lil->reiniciar_flujo_conservando_purga();
        const auto cap_main_post_reset = lil->arena_perfiles_capacidad();
        auto lil_copia = std::make_unique<ObsL>(*lil);
        std::printf("  Post reiniciar_flujo: main arena_cap=%zu | lil_local copia arena_cap=%zu usados=%zu\n",
                    cap_main_post_reset, lil_copia->arena_perfiles_capacidad(), lil_copia->arena_perfiles_usados());
        check(cap_main_post_reset >= std::size_t{1024} * 2 * 480, "Main conserva pre-reserva tras reiniciar");
        check(lil_copia->arena_perfiles_capacidad() == 0 && lil_copia->arena_perfiles_usados() == 0,
              "lil_local copia 0 bytes de arena");
    }

    std::printf("\nFALLAS T2.5: %zu\n", fallas);
    return fallas == 0 ? 0 : 1;
}
