// Sonda de revision (Code, 2026-11-05): el ObservadorLilyth de Argon TAL CUAL alimentado por el ObservadorSaltos REAL de la base.
// Mundo: FLUJO, ruido iid N(0,1); en la prueba, pulsos de +3 sd durante 200 pasos en t = 1000, 2000 y 3000 (vuelven a 0: verdad
// conocida, tau_rest ~ el fin del pulso). Se mide cuando abre cada ranura contra la alarma que abrio el evento en el operador delta.
#include <cstdio>
#include <memory>
#include <random>
#include <vector>

#include "ark/metrologia/observador_saltos.hpp"
#include "observador_lilyth.hpp"

using namespace ark::metrologia;
using Obs = ObservadorSaltos<>;
using Lil = ark::lilyth::ObservadorLilyth<Obs, 8, 512, 4096>;

int main() {
    std::mt19937_64 rng(20261105);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<double> hist(1 << 14);
    for (auto& v : hist) v = nd(rng);

    auto delta = std::make_unique<Obs>();
    ConfiguracionObservador cd{};
    cd.arl0_pasos_nativos = 5000.0;
    delta->calibrar(cd, hist);

    auto lil = std::make_unique<Lil>();
    ark::lilyth::ConfiguracionLilyth cl{};
    cl.representacion = ark::lilyth::ModoRepresentacionLilyth::CRUDA;
    cl.t_max_efectivo = 400;
    lil->configurar(cl);

    const int N = 4000;
    for (int t = 0; t < N; ++t) {
        double y = nd(rng);
        if ((t >= 1000 && t < 1200) || (t >= 2000 && t < 2200) || (t >= 3000 && t < 3200)) y += 3.0;
        const auto s = delta->paso(y);
        for (std::size_t j = 0; j < s.niveles; ++j)
            if (s.nivel[j].alarma.hay) std::printf("t=%5d  ALARMA nivel j=%zu signo=%+d t_cambio=%lld g=%.1f\n", t, j, s.nivel[j].alarma.signo,
                                                    static_cast<long long>(s.nivel[j].alarma.t_cambio), s.nivel[j].alarma.g);
        for (std::size_t i = 0; i < s.n_eventos; ++i) {
            const auto& e = s.eventos[i];
            std::printf("t=%5d  EVENTO emitido canal=%d estado=%s alarma.t=%lld t_cierre=%lld\n", t, e.canal, NOMBRE_ESTADO[static_cast<int>(e.estado)],
                        static_cast<long long>(e.alarma.t), static_cast<long long>(e.t_cierre));
        }
        Lil::SalidaPaso o{};
        lil->paso_con_salida_delta(y, true, 1.0, true, 1.0, s, o);
        if (o.abrio_ranura) {
            const auto& r = lil->ranuras()[o.indice_ranura_abierta];
            std::printf("t=%5d  LILYTH abre ranura: t_apertura=%lld t_cambio=%lld dY=%.3f Y*=%.3f\n", t, static_cast<long long>(r.t_apertura),
                        static_cast<long long>(r.t_cambio), r.cruda.dy, r.cruda.y_ancla);
        }
        for (std::size_t i = 0; i < o.ranuras_cerradas_en_paso; ++i) {
            const auto& r = o.cerradas[i];
            std::printf("t=%5d  LILYTH cierra: t_apertura=%lld desenlace=%d tau_rest=%.0f w_ancla=%zu\n", t, static_cast<long long>(r.t_apertura),
                        static_cast<int>(r.cruda.desenlace), r.cruda.tau_rest_consolidado, r.w_ancla);
        }
    }
    std::printf("eventos de nivel recibidos por Lilyth: %zu, desbordes: %zu\n", lil->eventos_nivel_recibidos(), lil->desbordes_ranura());
    return 0;
}
