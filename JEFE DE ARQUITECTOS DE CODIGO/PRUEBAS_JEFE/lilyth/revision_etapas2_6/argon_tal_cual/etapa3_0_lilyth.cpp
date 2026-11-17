// etapa3_0_lilyth.cpp -- Controles de la Etapa 3.0 (E1.21: continuacion vs replica) del Operador de Lilyth
// CONTRA EL OPERADOR DELTA REAL.
// Cada ranura que cierra el modulo, asi como la decision de si cada apertura del espejo es ranura nueva, continuacion o replica,
// se compara bit a bit con un motor de referencia independiente escrito aparte que aplica la regla E1.21 desde la serie grabada y
// las alarmas brutas del operador delta.
//
//   K1   identidad: el operador delta con y sin Lilyth da la misma salida bit a bit (Lilyth no lo posee ni lo toca)
//   K2   espejo: cada evento de canal 0 que emite el operador delta casa con una apertura del espejo (aperturas de ranura + continuaciones)
//   K2b  espejo con eventos de temperatura en medio (historia 2^20, varianza x9): 0 inconsistentes
//   K3   disparo en la alarma: t_apertura = alarma.t del evento (no su t_cierre)
//   K4   FLUJO, rama cruda: aperturas, continuaciones, replicas, ancla, dY, eta, primer toque, consolidado, desenlaces y tau_fin = ref (bit a bit)
//   K4b  el calor separa primer toque y consolidado
//   K5   E1.21 en el mundo de pulsos: las reaperturas del mismo signo dentro del mismo pulso pasan a continuaciones; las replicas que quedan
//        son las de signo contrario o tras restituir (tanto CORTADA_POR_REPLICA como RESTITUIDA al cortar tras consolidar)
//   K5b  eta_barra(tau) = media de las trayectorias de referencia (1e-12)
//   K5c  cuatro combinaciones a mano (signo igual/contrario x restituida/no restituida) + varias ranuras vivas (evalua contra la mas reciente)
//   K6   NIVEL + purga: espejo de la purga, reintegracion por segmentos y regla E1.21 = referencia independiente (bit a bit)
//   K6b  corte de la reintegracion dentro de la ventana del ancla da NO_MEDIDA_ANCLA_INSUFICIENTE
//   K7   memoria: con HIST_MAX chico, W_ancla + guarda fuera del anillo da NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA (= referencia)
//   K8   Aalen-Johansen contra la referencia independiente de Argon, con cortes por replica
//   K9   causalidad por truncado: el estado de Lilyth en t < 3000 no cambia si la serie sigue o se corta en 3000
//   K10  cero reservas de memoria en paso
//   K11  configuracion: rechaza NINGUNO, BLOQUES, t_max sin fijar; reloj desalineado se cuenta y no mide

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <random>
#include <vector>

#include "ark/lilyth/observador_lilyth.hpp"

static std::size_t g_reservas = 0;
static bool g_contar = false;
void* operator new(std::size_t n) {
    if (g_contar) ++g_reservas;
    if (void* p = std::malloc(n ? n : 1)) return p;
    throw std::bad_alloc{};
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

using namespace ark::metrologia;
using namespace ark::lilyth;

namespace {

using Obs = ObservadorSaltos<24, 64, 120>;   // L (memoria del GLR) = 64 != N_REG = 120
constexpr std::size_t NREG = 120, TMAXC = 512;
using Lil = ObservadorLilyth<Obs, 8, TMAXC, 4096>;
using LilChico = ObservadorLilyth<Obs, 8, TMAXC, 128>;

int fallos = 0;
void chk(bool ok, const char* id, const char* txt) {
    std::printf("[%s] %-4s %s\n", ok ? "OK   " : "FALLA", id, txt);
    if (!ok) ++fallos;
}
bool igual(double a, double b) { return (std::isnan(a) && std::isnan(b)) || std::memcmp(&a, &b, sizeof a) == 0; }

// ------------------------------------------------------------------------------------------------- mundos

struct Mundo {
    ConfiguracionObservador cd{};
    std::vector<double> xc;               // historia de calibracion
    std::vector<std::uint8_t> vc;
    std::vector<std::size_t> fc;
    std::vector<std::int64_t> dc;
    std::vector<double> x;                // prueba
    std::vector<std::uint8_t> v;
    std::vector<std::size_t> f;
    std::vector<std::int64_t> d;
    std::vector<std::uint8_t> calor;      // calor activo en el paso (entrada de la Etapa 1)
};

// FLUJO iid con pulsos de +3 sd durante 200 pasos; dos pulsos seguidos (5000 y 5300) para forzar replica del mismo signo TRAS restituir,
// y cambio de temperatura (3300..3700) que genera aperturas de signo contrario con ranura viva sin restituir.
Mundo mundo_flujo(std::uint64_t semilla, std::size_t n) {
    Mundo m;
    std::mt19937_64 rng(semilla);
    std::normal_distribution<double> nd(0.0, 1.0);
    m.cd.arl0_pasos_nativos = 5000.0;
    m.xc.resize(1 << 14);
    for (auto& z : m.xc) z = nd(rng);
    m.x.resize(n); m.v.assign(n, 1); m.calor.assign(n, 0);
    for (std::size_t t = 0; t < n; ++t) {
        double y = nd(rng);
        for (std::size_t p : {1000u, 2500u, 4000u, 5000u, 5300u}) if (t >= p && t < p + 200) y += 3.0;
        if (t >= 3300 && t < 3700) y *= 3.0;   // cambio de temperatura (varianza x9)
        m.x[t] = y;
        if ((t >= 1000 && t < 1150) || (t >= 3300 && t < 3500)) m.calor[t] = 1;
        if (t == 2600 || t == 4050) m.v[t] = 0;   // barras invalidas dentro de trayectorias
    }
    return m;
}

// NIVEL con purga (4 fases por dia de 200 pasos, varianza diurna x4), escalones de nivel de 10 sd durante 200 pasos y barras invalidas.
Mundo mundo_nivel(std::uint64_t semilla, std::size_t n) {
    Mundo m;
    std::mt19937_64 rng(semilla);
    std::normal_distribution<double> nd(0.0, 1.0);
    m.cd.arl0_pasos_nativos = 5000.0;
    m.cd.integracion = Integracion::NIVEL;
    m.cd.purga = true;
    m.cd.bins_fase = 4;
    m.cd.dias_minimos_purga = 20;
    const std::size_t dia_pasos = 200;
    auto sd_fase = [](std::size_t ff) { return ff == 1 ? 2.0 : 1.0; };
    const std::size_t nc = 1 << 14;
    double lvl = 100.0;
    for (std::size_t t = 0; t < nc; ++t) {
        const std::size_t ff = (t % dia_pasos) / (dia_pasos / 4);
        lvl += sd_fase(ff) * nd(rng);
        m.xc.push_back(lvl); m.vc.push_back(1); m.fc.push_back(ff); m.dc.push_back(static_cast<std::int64_t>(t / dia_pasos));
    }
    double base = lvl;
    m.calor.assign(n, 0);
    for (std::size_t t = 0; t < n; ++t) {
        const std::size_t tt = nc + t;
        const std::size_t ff = (tt % dia_pasos) / (dia_pasos / 4);
        base += sd_fase(ff) * nd(rng);
        double y = base;
        for (std::size_t p : {1000u, 2500u, 4000u}) if (t >= p && t < p + 200) y += 10.0 * sd_fase(ff) + 10.0;
        m.x.push_back(y);
        m.v.push_back((t == 2450 || t == 3000 || t == 4100) ? 0 : 1);
        m.f.push_back(ff); m.d.push_back(static_cast<std::int64_t>(tt / dia_pasos));
    }
    return m;
}

// ------------------------------------------------------------------------------------------------- corrida y referencia independiente

struct Apertura { std::int64_t t, t_cambio; int signo, j; std::size_t l; };

template <class L_>
struct Corrida {
    std::unique_ptr<Obs> delta = std::make_unique<Obs>();
    std::unique_ptr<L_> lil = std::make_unique<L_>();
    std::vector<typename L_::Ranura> cerradas;
    std::vector<Apertura> emitidas;           // aperturas segun los eventos de canal 0 del operador delta (alarma.t...)
    std::vector<Apertura> candidatas_delta;   // aperturas del espejo leidas directamente de las alarmas brutas del operador delta
    std::vector<std::int64_t> t_espejo_lil;   // pasos en que Lilyth abrio ranura, conto apertura o marco continuacion
    std::vector<std::int64_t> t_aperturas;    // pasos en que Lilyth efectivamente intento abrir ranura (no continuacion)
    std::vector<std::uint64_t> huella;        // por paso (K9)
    std::size_t temperatura = 0;              // eventos de canal 1 emitidos
};

std::uint64_t mezclar(std::uint64_t h, double v) { std::uint64_t b; std::memcpy(&b, &v, 8); return (h ^ b) * 1099511628211ull; }
std::uint64_t mezclar(std::uint64_t h, std::int64_t v) { return (h ^ static_cast<std::uint64_t>(v)) * 1099511628211ull; }

template <class L_>
void correr(Corrida<L_>& c, const Mundo& m, std::size_t n, std::size_t tmax, bool contar_reservas) {
    c.delta->calibrar(m.cd, m.xc, {}, m.vc, m.fc, m.dc);
    ConfiguracionLilyth cl{};
    cl.t_max_efectivo = tmax;
    c.lil->configurar(cl, m.cd, *c.delta, m.xc, m.vc, m.fc, m.dc);
    bool espejo_ref = false;
    std::size_t ap_prev = 0, cont_prev = 0;
    for (std::size_t t = 0; t < n; ++t) {
        const std::size_t fa = m.f.empty() ? 0 : m.f[t];
        const std::int64_t di = m.d.empty() ? 0 : m.d[t];
        const auto s = c.delta->paso(m.x[t], NaN, m.v[t] != 0, fa, di);
        if (!espejo_ref) {
            for (std::size_t j = 0; j < s.niveles && j < c.delta->escalas(); ++j) {
                if (!s.nivel[j].alarma.hay) continue;
                const auto& a = s.nivel[j].alarma;
                c.candidatas_delta.push_back({a.t, a.t_cambio, a.signo, static_cast<int>(j), a.l});
                espejo_ref = true;
                break;
            }
        }
        c.lil->fijar_calor_activo(m.calor[t] != 0);
        if (contar_reservas) g_contar = true;
        const auto o = c.lil->paso(s, m.x[t], m.v[t] != 0, NaN, fa, di);
        g_contar = false;
        const auto& k = c.lil->contadores();
        if (o.abrio_ranura || k.aperturas > ap_prev) c.t_aperturas.push_back(o.t);
        if (o.abrio_ranura || o.continuacion || k.aperturas > ap_prev || k.continuaciones > cont_prev) c.t_espejo_lil.push_back(o.t);
        ap_prev = k.aperturas;
        cont_prev = k.continuaciones;
        for (std::size_t i = 0; i < o.ranuras_cerradas; ++i) c.cerradas.push_back(c.lil->cerrada(i));
        for (std::size_t i = 0; i < s.n_eventos; ++i) {
            const auto& e = s.eventos[i];
            if (e.canal == 0) {
                c.emitidas.push_back({e.alarma.t, e.alarma.t_cambio, e.alarma.signo, e.escala_origen, e.alarma.l});
                espejo_ref = false;
            } else {
                ++c.temperatura;
            }
        }
        std::uint64_t h = 1469598103934665603ull;
        for (const auto& r : c.lil->ranuras()) {
            h = mezclar(h, static_cast<std::int64_t>(r.ocupada)); h = mezclar(h, r.t_apertura);
            for (const auto* b : {&r.purgada, &r.cruda}) {
                h = mezclar(h, static_cast<std::int64_t>(b->desenlace)); h = mezclar(h, static_cast<std::int64_t>(b->pasos_observados));
                h = mezclar(h, b->y_ancla); h = mezclar(h, b->dy); h = mezclar(h, b->tau_rest_primer_toque); h = mezclar(h, b->tau_rest_consolidado);
            }
        }
        h = mezclar(h, static_cast<std::int64_t>(k.continuaciones));
        h = mezclar(h, static_cast<std::int64_t>(k.replicas));
        h = mezclar(h, c.lil->aalen_johansen_cruda().total.rmst);
        c.huella.push_back(h);
    }
}

struct RamaRef {
    DesenlaceRanura des = DesenlaceRanura::INACTIVA;
    double y_ancla = NaN, dy = NaN, t1 = NaN, tc = NaN, min2 = NaN;
    std::size_t tau_fin = 0, pasos = 0;
    std::vector<double> eta{};
};

// Inicializacion independiente de una rama en su apertura ta.
RamaRef iniciar_rama_ref(const std::vector<double>& y, const std::vector<std::uint8_t>& val, const std::vector<std::uint32_t>& seg,
                         std::int64_t ta, std::int64_t tcam, int signo, std::size_t m, double tau_int, std::size_t hist_max, bool sin_purga) {
    RamaRef r;
    if (sin_purga) { r.des = DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA; return r; }
    if (!val[static_cast<std::size_t>(ta)]) {
        r.des = seg.empty() ? DesenlaceRanura::CENSURADA_DATO_INVALIDO : DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA;
        return r;
    }
    const std::int64_t guarda = static_cast<std::int64_t>(std::ceil(std::max(1.0, tau_int))) * static_cast<std::int64_t>(m);
    const std::int64_t tf = tcam - guarda - 1, ti = tf - static_cast<std::int64_t>(NREG * m) + 1;
    if (ti < 0 || tf < ti) { r.des = DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE; return r; }
    if (static_cast<std::size_t>(ta - ti) >= std::min<std::size_t>(hist_max, static_cast<std::size_t>(ta) + 1)) {
        r.des = DesenlaceRanura::NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA; return r;
    }
    std::vector<double> w;
    for (std::int64_t t = ti; t <= tf; ++t) {
        if (!seg.empty() && seg[static_cast<std::size_t>(t)] != seg[static_cast<std::size_t>(ta)]) { r.des = DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE; return r; }
        if (val[static_cast<std::size_t>(t)]) w.push_back(y[static_cast<std::size_t>(t)]);
    }
    if (w.empty()) { r.des = DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE; return r; }
    std::sort(w.begin(), w.end());
    const double med = w.size() % 2 ? w[w.size() / 2] : 0.5 * (w[w.size() / 2 - 1] + w[w.size() / 2]);
    const double dy = y[static_cast<std::size_t>(ta)] - med;
    if (!(dy * signo > 0.0)) { r.des = DesenlaceRanura::NO_MEDIDA_SIGNO_INCONSISTENTE; return r; }
    r.des = DesenlaceRanura::ABIERTA; r.y_ancla = med; r.dy = dy;
    return r;
}

void avanzar_rama_ref(RamaRef& r, const std::vector<double>& y, const std::vector<std::uint8_t>& val, const std::vector<std::uint32_t>& seg,
                      const std::vector<std::uint8_t>& calor, std::int64_t ta, std::size_t t, std::size_t tmax) {
    if (r.des != DesenlaceRanura::ABIERTA) return;
    const std::size_t tau = t - static_cast<std::size_t>(ta);
    const bool ok = val[t] && (seg.empty() || seg[t] == seg[static_cast<std::size_t>(ta)]);
    if (!ok) {
        const bool rest = std::isfinite(r.tc);
        r.des = rest ? DesenlaceRanura::RESTITUIDA : DesenlaceRanura::CENSURADA_DATO_INVALIDO;
        r.tau_fin = rest ? static_cast<std::size_t>(r.tc) : r.pasos;
        return;
    }
    const double e = (y[t] - r.y_ancla) / r.dy;
    r.eta.push_back(e); r.pasos = tau;
    if (tau >= 2 && (!std::isfinite(r.min2) || e < r.min2)) r.min2 = e;
    if (e <= 0.5) {
        if (!std::isfinite(r.t1)) r.t1 = static_cast<double>(tau);
        if (!calor[t] && !std::isfinite(r.tc)) r.tc = static_cast<double>(tau);
    }
    if (tau >= tmax) {
        const bool rest = std::isfinite(r.tc);
        r.des = rest ? DesenlaceRanura::RESTITUIDA : DesenlaceRanura::CENSURADA_HORIZONTE;
        r.tau_fin = rest ? static_cast<std::size_t>(r.tc) : tau;
    }
}

void cortar_rama_ref(RamaRef& r, std::size_t tau) {
    if (r.des != DesenlaceRanura::ABIERTA) return;
    const bool rest = std::isfinite(r.tc);
    r.des = rest ? DesenlaceRanura::RESTITUIDA : DesenlaceRanura::CORTADA_POR_REPLICA;
    r.tau_fin = rest ? static_cast<std::size_t>(r.tc) : tau;
}

// Referencia independiente de una rama dada la lista de aperturas efectivas de ranura (que cortan por replica).
RamaRef rama_ref(const std::vector<double>& y, const std::vector<std::uint8_t>& val, const std::vector<std::uint32_t>& seg,
                 const std::vector<std::uint8_t>& calor, const std::vector<std::int64_t>& aperturas, std::int64_t ta, std::int64_t tcam,
                 int signo, std::size_t m, double tau_int, std::size_t tmax, std::size_t hist_max, bool sin_purga) {
    RamaRef r = iniciar_rama_ref(y, val, seg, ta, tcam, signo, m, tau_int, hist_max, sin_purga);
    if (r.des != DesenlaceRanura::ABIERTA) return r;
    for (std::size_t tau = 1;; ++tau) {
        const std::size_t t = static_cast<std::size_t>(ta) + tau;
        if (t >= y.size()) break;
        avanzar_rama_ref(r, y, val, seg, calor, ta, t, tmax);
        if (r.des != DesenlaceRanura::ABIERTA) break;
        if (std::find(aperturas.begin(), aperturas.end(), static_cast<std::int64_t>(t)) != aperturas.end()) {
            cortar_rama_ref(r, tau);
            break;
        }
    }
    return r;
}

template <class R>
bool rama_igual(const R& b, const RamaRef& r) {
    if (b.desenlace != r.des) return false;
    if (es_no_medida(r.des)) return true;
    if (!igual(b.y_ancla, r.y_ancla) || !igual(b.dy, r.dy) || !igual(b.tau_rest_primer_toque, r.t1) || !igual(b.tau_rest_consolidado, r.tc) ||
        !igual(b.min_eta_desde_2, r.min2) || b.tau_fin != r.tau_fin || b.pasos_observados != r.pasos)
        return false;
    for (std::size_t i = 0; i < r.eta.size(); ++i) if (!igual(b.eta[i], r.eta[i])) return false;
    return true;
}

// Purga independiente (PurgaDiurna nueva) sobre calibracion + prueba, con la reintegracion por segmentos del contrato.
void purgada_ref(const Mundo& m, std::vector<double>& yp, std::vector<std::uint8_t>& vp, std::vector<std::uint32_t>& sp) {
    PurgaDiurna<1440> pu(m.cd.bins_fase, m.cd.dias_minimos_purga);
    bool hp = false; double xp = NaN;
    for (std::size_t i = 0; i < m.xc.size(); ++i) {
        const bool ok = m.vc[i] && std::isfinite(m.xc[i]);
        const bool v = ok && hp;
        const double y = v ? m.xc[i] - xp : NaN;
        hp = ok; if (ok) xp = m.xc[i];
        bool q = false; (void)pu.paso(v ? y : NaN, m.fc[i], m.dc[i], q);
    }
    hp = false; xp = NaN;
    double nivel = NaN; bool vivo = false; std::uint32_t s = 0;
    for (std::size_t t = 0; t < m.x.size(); ++t) {
        const bool ok = m.v[t] && std::isfinite(m.x[t]);
        const bool v = ok && hp;
        const double y = v ? m.x[t] - xp : NaN;
        hp = ok; if (ok) xp = m.x[t];
        bool q = false;
        const double o = pu.paso(v ? y : NaN, m.f[t], m.d[t], q);
        if (v && q) {
            if (!vivo) { nivel = m.x[t] - y; vivo = true; }
            nivel += o;
            yp.push_back(nivel); vp.push_back(1);
        } else {
            if (vivo) ++s;
            vivo = false;
            yp.push_back(NaN); vp.push_back(0);
        }
        sp.push_back(s);
    }
}

// Motor de referencia independiente para el flujo completo de ranuras, continuaciones (E1.21) y replicas.
struct RanuraRefCompleta {
    std::uint64_t id = 0;
    std::int64_t ta = -1, tcam = -1;
    int signo = 0, j = 0;
    RamaRef purgada{}, cruda{};
    bool viva() const { return purgada.des == DesenlaceRanura::ABIERTA || cruda.des == DesenlaceRanura::ABIERTA; }
    bool restituida() const {
        return (purgada.des == DesenlaceRanura::ABIERTA && std::isfinite(purgada.tc)) ||
               (cruda.des == DesenlaceRanura::ABIERTA && std::isfinite(cruda.tc));
    }
};

template <class L_>
bool comparar_contra_ref(const Corrida<L_>& c, const Mundo& m, std::size_t tmax, std::size_t hist_max, bool con_purgada, std::size_t& comparadas,
                         std::size_t& con_desenlace, DesenlaceRanura buscado) {
    std::vector<double> yp; std::vector<std::uint8_t> vp; std::vector<std::uint32_t> sp;
    if (con_purgada) purgada_ref(m, yp, vp, sp);
    std::vector<std::uint8_t> vc(m.v);
    for (std::size_t t = 0; t < m.x.size(); ++t) if (!std::isfinite(m.x[t])) vc[t] = 0;

    // 1. Reconstruir de manera independiente que aperturas candidatas del espejo son continuaciones, replicas o nuevas ranuras.
    std::vector<RanuraRefCompleta> vivas;
    std::vector<std::int64_t> aperturas_ref;
    std::size_t cont_ref = 0, rep_ref = 0, ap_count_ref = 0;
    std::uint64_t seq_ref = 0;
    std::size_t idx_cand = 0;

    for (std::size_t t = 0; t < m.x.size(); ++t) {
        // Avanzar vivas en t
        for (auto& rv : vivas) {
            if (!rv.viva() || static_cast<std::int64_t>(t) <= rv.ta) continue;
            if (con_purgada) avanzar_rama_ref(rv.purgada, yp, vp, sp, m.calor, rv.ta, t, tmax);
            avanzar_rama_ref(rv.cruda, m.x, vc, {}, m.calor, rv.ta, t, tmax);
        }
        vivas.erase(std::remove_if(vivas.begin(), vivas.end(), [](const RanuraRefCompleta& r) { return !r.viva(); }), vivas.end());

        if (idx_cand < c.candidatas_delta.size() && c.candidatas_delta[idx_cand].t == static_cast<std::int64_t>(t)) {
            const Apertura ap = c.candidatas_delta[idx_cand++];
            const RanuraRefCompleta* mas_reciente = nullptr;
            for (const auto& rv : vivas) {
                if (mas_reciente == nullptr || rv.id > mas_reciente->id) mas_reciente = &rv;
            }
            bool es_cont = false;
            if (mas_reciente != nullptr) {
                if (ap.signo == mas_reciente->signo && !mas_reciente->restituida()) {
                    es_cont = true;
                    ++cont_ref;
                } else {
                    ++rep_ref;
                    for (auto& rv : vivas) {
                        const std::size_t tau = t - static_cast<std::size_t>(rv.ta);
                        cortar_rama_ref(rv.purgada, tau);
                        cortar_rama_ref(rv.cruda, tau);
                    }
                    vivas.clear();
                }
            }
            if (!es_cont) {
                ++ap_count_ref;
                aperturas_ref.push_back(ap.t);
                const auto& f = c.delta->ficha_nivel(static_cast<std::size_t>(ap.j));
                RanuraRefCompleta nr;
                nr.id = ++seq_ref;
                nr.ta = ap.t; nr.tcam = ap.t_cambio; nr.signo = ap.signo; nr.j = ap.j;
                nr.purgada = con_purgada ? iniciar_rama_ref(yp, vp, sp, ap.t, ap.t_cambio, ap.signo, static_cast<std::size_t>(f.m), f.nulo.tau_int, hist_max, false)
                                         : RamaRef{DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA};
                nr.cruda = iniciar_rama_ref(m.x, vc, {}, ap.t, ap.t_cambio, ap.signo, static_cast<std::size_t>(f.m), f.nulo.tau_int, hist_max, false);
                if (nr.viva()) vivas.push_back(std::move(nr));
            }
        }
    }

    const auto& k = c.lil->contadores();
    bool ok = (k.continuaciones == cont_ref && k.replicas == rep_ref && k.aperturas == ap_count_ref && c.t_aperturas == aperturas_ref);
    if (!ok) {
        std::printf("      contadores distintos: lil(ap=%zu, cont=%zu, rep=%zu) vs ref(ap=%zu, cont=%zu, rep=%zu)\n",
                    k.aperturas, k.continuaciones, k.replicas, ap_count_ref, cont_ref, rep_ref);
    }

    comparadas = con_desenlace = 0;
    for (const auto& r : c.cerradas) {
        const auto& f = c.delta->ficha_nivel(static_cast<std::size_t>(r.escala_origen));
        const auto rc = rama_ref(m.x, vc, {}, m.calor, aperturas_ref, r.t_apertura, r.t_cambio, r.signo, static_cast<std::size_t>(f.m),
                                 f.nulo.tau_int, tmax, hist_max, false);
        if (!rama_igual(r.cruda, rc)) {
            ok = false;
            std::printf("      cruda distinta en t_a=%lld: modulo %s, ref %s\n", static_cast<long long>(r.t_apertura),
                        NOMBRE_DESENLACE[static_cast<int>(r.cruda.desenlace)], NOMBRE_DESENLACE[static_cast<int>(rc.des)]);
        }
        if (rc.des == buscado) ++con_desenlace;
        const auto rp = con_purgada ? rama_ref(yp, vp, sp, m.calor, aperturas_ref, r.t_apertura, r.t_cambio, r.signo, static_cast<std::size_t>(f.m),
                                               f.nulo.tau_int, tmax, hist_max, false)
                                    : RamaRef{DesenlaceRanura::NO_MEDIDA_PURGA_INMADURA};
        if (!rama_igual(r.purgada, rp)) {
            ok = false;
            std::printf("      purgada distinta en t_a=%lld: modulo %s, ref %s\n", static_cast<long long>(r.t_apertura),
                        NOMBRE_DESENLACE[static_cast<int>(r.purgada.desenlace)], NOMBRE_DESENLACE[static_cast<int>(rp.des)]);
        }
        if (rp.des == buscado) ++con_desenlace;
        ++comparadas;
    }
    return ok;
}

// Referencia independiente de Aalen-Johansen (la de Argon, sin cambios).
struct EvAJ { int tipo; std::size_t t; };
void aj_ref(const std::vector<EvAJ>& ev, std::size_t tmax, std::vector<double>& cr, std::vector<double>& cp, double& rmst) {
    cr.assign(tmax, 0.0); cp.assign(tmax, 0.0);
    double s = 1.0, a = 0.0, b = 0.0; rmst = 0.0;
    for (std::size_t tau = 1; tau <= tmax; ++tau) {
        std::size_t n = 0, d1 = 0, d2 = 0;
        for (const auto& e : ev) {
            const std::size_t tf = std::min(e.t, tmax);
            if (tf >= tau) ++n;
            if (tf == tau && e.tipo == 1) ++d1;
            if (tf == tau && e.tipo == 2) ++d2;
        }
        if (n) { a += s * double(d1) / double(n); b += s * double(d2) / double(n); s *= double(n - d1 - d2) / double(n); }
        cr[tau - 1] = a; cp[tau - 1] = b; rmst += 1.0 - a;
    }
}

}  // namespace

int main() {
    const std::size_t N = 6000, TM = 400;

    // ---------------- FLUJO
    const Mundo mf = mundo_flujo(20261105, N);
    Corrida<Lil> cf;
    correr(cf, mf, N, TM, true);

    // K1: otro operador delta, mismas entradas, sin Lilyth
    {
        auto solo = std::make_unique<Obs>();
        solo->calibrar(mf.cd, mf.xc);
        auto con = std::make_unique<Obs>();
        con->calibrar(mf.cd, mf.xc);
        auto lil = std::make_unique<Lil>();
        ConfiguracionLilyth cl{}; cl.t_max_efectivo = TM;
        lil->configurar(cl, mf.cd, *con, mf.xc);
        bool id = true;
        for (std::size_t t = 0; t < N && id; ++t) {
            const auto a = solo->paso(mf.x[t], NaN, mf.v[t] != 0);
            const auto b = con->paso(mf.x[t], NaN, mf.v[t] != 0);
            (void)lil->paso(b, mf.x[t], mf.v[t] != 0);
            id = a.t == b.t && a.n_eventos == b.n_eventos && a.evento_abierto == b.evento_abierto && a.niveles == b.niveles;
            for (std::size_t j = 0; id && j < a.niveles; ++j)
                id = a.nivel[j].alarma.hay == b.nivel[j].alarma.hay && igual(a.nivel[j].g, b.nivel[j].g) && igual(a.temp[j].log_beta, b.temp[j].log_beta);
        }
        chk(id, "K1", "identidad del operador delta con y sin Lilyth (6000 pasos, salida por escala)");
    }
    {
        const auto& k = cf.lil->contadores();
        bool casa = k.espejo_inconsistente == 0 && k.cierres_espejo == cf.emitidas.size() && !cf.emitidas.empty() &&
                    (k.aperturas + k.continuaciones) == cf.t_espejo_lil.size();
        for (const auto& e : cf.emitidas)
            casa = casa && std::find(cf.t_espejo_lil.begin(), cf.t_espejo_lil.end(), e.t) != cf.t_espejo_lil.end();
        std::printf("      eventos de canal 0: %zu, de canal 1: %zu, espejo Lilyth: %zu (aperturas=%zu, continuaciones=%zu, replicas=%zu)\n",
                    cf.emitidas.size(), cf.temperatura, cf.t_espejo_lil.size(), k.aperturas, k.continuaciones, k.replicas);

        chk(casa, "K2", "espejo: cada evento de canal 0 casa con una apertura del espejo (aperturas + continuaciones), 0 inconsistentes");
        bool alarma = true;
        for (const auto& r : cf.cerradas) alarma = alarma && r.t_apertura == r.ficha_apertura.alarma.t;
        bool antes_de_cierre = true;
        for (const auto& e : cf.emitidas) antes_de_cierre = antes_de_cierre && std::find(cf.t_espejo_lil.begin(), cf.t_espejo_lil.end(), e.t) != cf.t_espejo_lil.end();
        chk(alarma && antes_de_cierre && cf.cerradas.front().t_apertura >= 1000 && cf.cerradas.front().t_apertura <= 1010,
            "K3", "disparo en la alarma: t_apertura = alarma.t (el primer pulso abre en 1000-1010, no en su clasificacion)");
    }
    {
        std::size_t n = 0, nb = 0;
        const bool ok = comparar_contra_ref(cf, mf, TM, 4096, false, n, nb, DesenlaceRanura::CORTADA_POR_REPLICA);
        bool calor_sep = false;
        bool pulso_1000_unico = false;
        std::size_t ranuras_en_pulso_1000 = 0;
        bool replica_tras_restituir_5300 = false;
        for (const auto& r : cf.cerradas) {
            std::printf("      t_a=%5lld t_cambio=%5lld j=%d signo=%+d dY=%7.3f Y*=%7.3f primer=%4.0f consolidado=%4.0f %s tau_fin=%zu\n",
                        static_cast<long long>(r.t_apertura), static_cast<long long>(r.t_cambio), r.escala_origen, r.signo, r.cruda.dy, r.cruda.y_ancla,
                        r.cruda.tau_rest_primer_toque, r.cruda.tau_rest_consolidado, NOMBRE_DESENLACE[static_cast<int>(r.cruda.desenlace)], r.cruda.tau_fin);
            if (r.t_apertura >= 1000 && r.t_apertura < 1200) {
                ++ranuras_en_pulso_1000;
                // Con E1.21 las reaperturas de 1065 y 1195 son continuaciones: el pulso de 1000 no se corta y consolida al terminar el pulso (~1200)
                calor_sep = std::isfinite(r.cruda.tau_rest_primer_toque) && r.t_apertura + r.cruda.tau_rest_primer_toque < 1150 &&
                            !(r.cruda.tau_rest_consolidado <= r.cruda.tau_rest_primer_toque);
                pulso_1000_unico = (r.cruda.desenlace == DesenlaceRanura::RESTITUIDA);
            }
            if (r.t_apertura >= 5000 && r.t_apertura < 5200 && r.signo_replica_corte == +1 && r.cruda.desenlace == DesenlaceRanura::RESTITUIDA) {
                replica_tras_restituir_5300 = true;
            }
        }
        pulso_1000_unico = pulso_1000_unico && (ranuras_en_pulso_1000 == 1);
        std::printf("      ranuras comparadas %zu, cortadas por replica %zu, continuaciones %zu, replicas %zu\n",
                    n, nb, cf.lil->contadores().continuaciones, cf.lil->contadores().replicas);
        chk(ok && n > 0, "K4", "FLUJO, rama cruda: aperturas, continuaciones, replicas, ancla, dY, eta, desenlace y tau_fin = referencia bit a bit");
        chk(calor_sep, "K4b", "el calor separa: primer toque dentro del calor y consolidado despues (pulso de 1000)");
        chk(cf.lil->contadores().continuaciones >= 2 && cf.lil->contadores().replicas >= 1 && nb >= 1 && pulso_1000_unico && replica_tras_restituir_5300,
            "K5", "E1.21 en mundo de pulsos: reaperturas del mismo pulso son continuaciones; quedan replicas de signo contrario y tras restituir");
    }
    // K5b: eta_barra = media de las trayectorias de referencia mientras la rama estuvo abierta
    {
        std::vector<std::uint8_t> vc(mf.v);
        std::vector<double> suma(TM, 0.0); std::vector<std::size_t> cuenta(TM, 0);
        for (const auto& r : cf.cerradas) {
            const auto& f = cf.delta->ficha_nivel(static_cast<std::size_t>(r.escala_origen));
            const auto rc = rama_ref(mf.x, vc, {}, mf.calor, cf.t_aperturas, r.t_apertura, r.t_cambio, r.signo, static_cast<std::size_t>(f.m),
                                     f.nulo.tau_int, TM, 4096, false);
            if (es_no_medida(rc.des)) continue;
            for (std::size_t i = 0; i < rc.eta.size(); ++i) { suma[i] += rc.eta[i]; ++cuenta[i]; }
        }
        bool ok = true;
        for (std::size_t tau = 1; tau <= TM; ++tau) {
            const double ref = cuenta[tau - 1] ? suma[tau - 1] / double(cuenta[tau - 1]) : NaN;
            if (!(std::abs(cf.lil->eta_barra_cruda(tau) - ref) <= 1e-12 * std::max(1.0, std::abs(ref))) &&
                !(std::isnan(ref) && std::isnan(cf.lil->eta_barra_cruda(tau))))
                ok = false;
        }
        chk(ok, "K5b", "eta_barra(tau) = media de las trayectorias de referencia (1e-12)");
    }

    // ---------------- K5c: caso a mano de las 4 combinaciones (signo igual/contrario x restituida/no restituida)
    //                       + caso con varias ranuras vivas (evalua contra la mas reciente, no contra la mas antigua)
    {
        auto delta = std::make_unique<Obs>();
        delta->calibrar(mf.cd, mf.xc);

        auto inyectar_paso = [](Lil& lil, std::int64_t t, double x, bool disparar, int signo, std::int64_t t_cam, bool cerrar_espejo) {
            Obs::Salida s{};
            s.t = t;
            s.niveles = 1;
            if (disparar) {
                s.nivel[0].alarma.hay = true;
                s.nivel[0].alarma.t = t;
                s.nivel[0].alarma.t_cambio = t_cam;
                s.nivel[0].alarma.signo = signo;
                s.nivel[0].alarma.l = 5;
                s.nivel[0].alarma.amplitud = signo * 4.0;
                s.nivel[0].alarma.se_amplitud = 0.5;
                s.nivel[0].alarma.g = 25.0;
            }
            if (cerrar_espejo) {
                s.n_eventos = 1;
                s.eventos[0].canal = 0;
                s.eventos[0].escala_origen = 0;
                s.eventos[0].alarma = s.nivel[0].alarma;
            }
            return lil.paso(s, x, true);
        };

        auto cerrar_espejo_en = [](Lil& lil, std::int64_t t, double x, std::int64_t t_ap, std::int64_t t_cam, int signo) {
            Obs::Salida s{};
            s.t = t;
            s.niveles = 1;
            s.n_eventos = 1;
            s.eventos[0].canal = 0;
            s.eventos[0].escala_origen = 0;
            s.eventos[0].alarma.hay = true;
            s.eventos[0].alarma.t = t_ap;
            s.eventos[0].alarma.t_cambio = t_cam;
            s.eventos[0].alarma.signo = signo;
            s.eventos[0].alarma.l = 5;
            return lil.paso(s, x, true);
        };

        // Preparar base: t = 0..199 en 0.0; en t = 200 abre ranura +1 con x = +10.0 (Y* = 0, dY = +10);
        // en t = 201..330 sigue en +10.0 (eta = 1.0 > 0.5, no restituida, y la ventana [209, 328] queda en +10.0).
        // Si restituir_en_331 == true, en t = 331 baja a +2.0 (eta = 0.2 <= 0.5, ya restituida); si no, sigue en +10.0.
        auto preparar_ranura_positiva = [&](Lil& lil, bool cortar_con_replica, bool restituir_en_331) {
            ConfiguracionLilyth cl{};
            cl.representacion = ModoRepresentacionLilyth::CRUDA;
            cl.t_max_efectivo = 200;
            cl.cortar_abiertas_con_replica = cortar_con_replica;
            lil.configurar(cl, mf.cd, *delta, mf.xc);
            for (std::int64_t t = 0; t < 200; ++t) inyectar_paso(lil, t, 0.0, false, 0, 0, false);
            inyectar_paso(lil, 200, 10.0, true, +1, 196, false);
            cerrar_espejo_en(lil, 201, 10.0, 200, 196, +1);
            for (std::int64_t t = 202; t <= 330; ++t) inyectar_paso(lil, t, 10.0, false, 0, 0, false);
            inyectar_paso(lil, 331, restituir_en_331 ? 2.0 : 10.0, false, 0, 0, false);
        };

        // Combinacion 1: mismo signo (+1) y NO restituida -> CONTINUACION (no corta, no abre ranura)
        Lil l1;
        preparar_ranura_positiva(l1, true, false);
        const auto o1 = inyectar_paso(l1, 332, 12.0, true, +1, 330, false);
        const bool c1_ok = o1.continuacion && !o1.replica && !o1.abrio_ranura && o1.ranuras_cerradas == 0 &&
                           l1.contadores().continuaciones == 1 && l1.contadores().replicas == 0 && l1.contadores().aperturas == 1 &&
                           l1.ranuras()[0].ocupada;

        // Combinacion 2: mismo signo (+1) y YA restituida -> REPLICA (cierra la anterior como RESTITUIDA en tau=131 y abre nueva)
        Lil l2;
        preparar_ranura_positiva(l2, true, true);
        const auto o2 = inyectar_paso(l2, 332, 12.0, true, +1, 330, false);
        const bool c2_ok = !o2.continuacion && o2.replica && o2.abrio_ranura && o2.ranuras_cerradas == 1 &&
                           l2.cerrada(0).cruda.desenlace == DesenlaceRanura::RESTITUIDA && l2.cerrada(0).cruda.tau_fin == 131 &&
                           l2.contadores().continuaciones == 0 && l2.contadores().replicas == 1 && l2.contadores().aperturas == 2;

        // Combinacion 3: signo contrario (-1) y NO restituida -> REPLICA (en x = +6.0, eta de la ranura +1 es 0.6 > 0.5: sigue sin restituir,
        // y para la nueva ranura -1 su Y* en [209, 328] es +10.0, dando dY = -4.0 < 0; corta la anterior como CORTADA_POR_REPLICA en tau=132)
        Lil l3;
        preparar_ranura_positiva(l3, true, false);
        const auto o3 = inyectar_paso(l3, 332, 6.0, true, -1, 330, false);
        const bool c3_ok = !o3.continuacion && o3.replica && o3.abrio_ranura && o3.ranuras_cerradas == 1 &&
                           l3.cerrada(0).cruda.desenlace == DesenlaceRanura::CORTADA_POR_REPLICA && l3.cerrada(0).cruda.tau_fin == 132 &&
                           l3.cerrada(0).signo_replica_corte == -1 &&
                           l3.contadores().continuaciones == 0 && l3.contadores().replicas == 1 && l3.contadores().aperturas == 2;

        // Combinacion 4: signo contrario (-1) y YA restituida -> REPLICA (cierra la anterior como RESTITUIDA en tau=131 y abre nueva)
        Lil l4;
        preparar_ranura_positiva(l4, true, true);
        const auto o4 = inyectar_paso(l4, 332, 6.0, true, -1, 330, false);
        const bool c4_ok = !o4.continuacion && o4.replica && o4.abrio_ranura && o4.ranuras_cerradas == 1 &&
                           l4.cerrada(0).cruda.desenlace == DesenlaceRanura::RESTITUIDA && l4.cerrada(0).cruda.tau_fin == 131 &&
                           l4.cerrada(0).signo_replica_corte == -1 &&
                           l4.contadores().continuaciones == 0 && l4.contadores().replicas == 1 && l4.contadores().aperturas == 2;

        // Caso 5: varias ranuras vivas (cortar_abiertas_con_replica = false).
        // Ranura 0 (mas antigua, t=200, signo +1, Y*=0, dY=+10) y Ranura 1 (mas reciente, t=332, signo -1, Y*=+10, dY=-4):
        // en x = +6.0 ninguna de las dos ha restituido (eta_0 = 0.6 > 0.5, eta_1 = 1.0 > 0.5).
        // En t = 334 llega otra apertura de signo -1 en x = +6.0: contra la mas reciente (Ranura 1, signo -1) es CONTINUACION;
        // si comparase contra la mas antigua (Ranura 0, signo +1) daria replica.
        Lil l5;
        preparar_ranura_positiva(l5, false, false);
        const auto o5_rep = inyectar_paso(l5, 332, 6.0, true, -1, 330, false);
        cerrar_espejo_en(l5, 333, 6.0, 332, 330, -1);
        const auto o5_cont = inyectar_paso(l5, 334, 6.0, true, -1, 333, false);
        const bool c5_ok = o5_rep.replica && o5_rep.abrio_ranura && l5.ranuras()[0].ocupada && l5.ranuras()[1].ocupada &&
                           o5_cont.continuacion && !o5_cont.replica && !o5_cont.abrio_ranura &&
                           l5.contadores().continuaciones == 1 && l5.contadores().replicas == 1 && l5.contadores().aperturas == 2;

        // Caso 6: signo inconsistente (alarma +1 pero dY = -2.0 < 0 -> NO_MEDIDA_SIGNO_INCONSISTENTE)
        Lil l6;
        preparar_ranura_positiva(l6, true, true);
        const auto o6 = inyectar_paso(l6, 332, -2.0, true, +1, 330, false);
        const bool c6_ok = !o6.abrio_ranura && l6.contadores().no_medida[static_cast<int>(DesenlaceRanura::NO_MEDIDA_SIGNO_INCONSISTENTE)] == 1;

        chk(c1_ok && c2_ok && c3_ok && c4_ok && c5_ok && c6_ok, "K5c",
            "4 combinaciones a mano (+/no rest, +/rest, -/no rest, -/rest), ranura mas reciente y signo inconsistente");
    }

    // ---------------- K2b: el espejo con eventos de temperatura (historia 2^20, varianza x9)
    {
        Mundo mt;
        std::mt19937_64 rng(4242);
        std::normal_distribution<double> nd(0.0, 1.0);
        mt.cd.arl0_pasos_nativos = 5000.0;
        mt.cd.escalas_max = 1;
        mt.cd.medir_tau50 = false;
        mt.xc.resize(350000);
        for (auto& z : mt.xc) z = nd(rng);
        const std::size_t nt = 30000;
        mt.x.resize(nt); mt.v.assign(nt, 1); mt.calor.assign(nt, 0);
        for (std::size_t t = 0; t < nt; ++t) mt.x[t] = nd(rng) * (t >= 5000 && t < 20000 ? 3.0 : 1.0);
        Corrida<Lil> ct;
        correr(ct, mt, nt, TM, false);
        const auto& k = ct.lil->contadores();
        bool casa = k.espejo_inconsistente == 0 && k.cierres_espejo == ct.emitidas.size() && ct.temperatura > 0 &&
                    (k.aperturas + k.continuaciones) == ct.t_espejo_lil.size();
        for (const auto& e : ct.emitidas) casa = casa && std::find(ct.t_espejo_lil.begin(), ct.t_espejo_lil.end(), e.t) != ct.t_espejo_lil.end();
        std::printf("      eventos de canal 0: %zu, de canal 1: %zu, espejo Lilyth: %zu (aperturas=%zu, continuaciones=%zu)\n",
                    ct.emitidas.size(), ct.temperatura, ct.t_espejo_lil.size(), k.aperturas, k.continuaciones);
        chk(casa, "K2b", "espejo con eventos de temperatura en medio (historia 2^20, varianza x9): 0 inconsistentes, cierres = eventos de canal 0");
    }

    // ---------------- NIVEL + purga
    const Mundo mn = mundo_nivel(77, N);
    Corrida<Lil> cn;
    correr(cn, mn, N, TM, false);
    {
        std::size_t n = 0, ninsuf = 0;
        const bool ok = comparar_contra_ref(cn, mn, TM, 4096, true, n, ninsuf, DesenlaceRanura::NO_MEDIDA_ANCLA_INSUFICIENTE);
        bool purgada_medida = false;
        for (const auto& r : cn.cerradas) purgada_medida = purgada_medida || r.purgada.desenlace == DesenlaceRanura::RESTITUIDA || r.purgada.desenlace == DesenlaceRanura::CENSURADA_HORIZONTE;
        std::printf("      ranuras comparadas %zu, con ancla insuficiente por corte de segmento %zu, continuaciones %zu\n",
                    n, ninsuf, cn.lil->contadores().continuaciones);
        chk(ok && n > 0 && purgada_medida, "K6", "NIVEL + purga: espejo de la purga, reintegracion por segmentos y E1.21 = referencia independiente (bit a bit)");
        chk(ninsuf >= 1, "K6b", "un corte de la reintegracion dentro de la ventana del ancla da NO_MEDIDA_ANCLA_INSUFICIENTE en la purgada");
    }

    // ---------------- memoria chica
    {
        Corrida<LilChico> cc;
        correr(cc, mf, N, TM, false);
        std::size_t n = 0, nfuera = 0;
        const bool ok = comparar_contra_ref(cc, mf, TM, 128, false, n, nfuera, DesenlaceRanura::NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA);
        chk(ok && nfuera >= 1, "K7", "HIST_MAX 128: W_ancla + guarda fuera del anillo da NO_MEDIDA_ANCLA_FUERA_DE_MEMORIA (= referencia)");
    }

    // ---------------- Aalen-Johansen
    {
        std::vector<EvAJ> ev;
        for (const auto& r : cf.cerradas) {
            const auto d = r.cruda.desenlace;
            if (d == DesenlaceRanura::INACTIVA || es_no_medida(d)) continue;
            ev.push_back({d == DesenlaceRanura::RESTITUIDA ? 1 : d == DesenlaceRanura::CORTADA_POR_REPLICA ? 2 : 0, r.cruda.tau_fin});
        }
        std::vector<double> cr, cp; double rm = 0.0;
        aj_ref(ev, TM, cr, cp, rm);
        const auto& aj = cf.lil->aalen_johansen_cruda().total;
        bool ok = std::abs(aj.rmst - rm) < 1e-9;
        for (std::size_t i = 0; i < TM; ++i)
            ok = ok && std::abs(aj.cif_rest[i] - cr[i]) < 1e-12 && std::abs(aj.cif_replica[i] - cp[i]) < 1e-12 &&
                 std::abs(aj.supervivientes[i] + aj.cif_rest[i] + aj.cif_replica[i] - 1.0) < 1e-12;
        std::printf("      eventos en Aalen-Johansen %zu, RMST %.3f\n", ev.size(), aj.rmst);
        chk(ok && !ev.empty(), "K8", "Aalen-Johansen de la corrida = referencia independiente de Argon (CIF, S + CIF = 1, RMST)");
    }

    // ---------------- causalidad por truncado
    {
        Corrida<Lil> ct;
        correr(ct, mf, N / 2, TM, false);
        bool ok = true;
        for (std::size_t t = 0; t < N / 2; ++t) ok = ok && ct.huella[t] == cf.huella[t];
        chk(ok, "K9", "causalidad por truncado: el estado de Lilyth en t < 3000 es el mismo con la serie cortada en 3000");
    }

    chk(g_reservas == 0, "K10", "cero reservas de memoria en paso (6000 pasos)");

    // ---------------- configuracion
    {
        auto d = std::make_unique<Obs>();
        d->calibrar(mf.cd, mf.xc);
        auto l = std::make_unique<Lil>();
        ConfiguracionLilyth c{};
        c.t_max_efectivo = TM;
        c.nulo_sin_memoria = LILYTH_NULO_NINGUNO;
        const bool a = !l->configurar(c, mf.cd, *d, mf.xc);
        c.nulo_sin_memoria = LILYTH_NULO_BLOQUES;
        const bool b = !l->configurar(c, mf.cd, *d, mf.xc);
        c.nulo_sin_memoria = LILYTH_NULO_SIGNOS; c.t_max_efectivo = 0;
        const bool e = !l->configurar(c, mf.cd, *d, mf.xc);
        c.t_max_efectivo = TM;
        const bool g = l->configurar(c, mf.cd, *d, mf.xc);
        (void)d->paso(0.0); const auto s1 = d->paso(0.0);
        const auto o = l->paso(s1, 0.0, true);
        chk(a && b && e && g && !o.emitida && l->contadores().pasos_desalineados == 1, "K11",
            "configuracion: rechaza NINGUNO, BLOQUES y T_max sin fijar; un paso desalineado se cuenta y no se mide");
    }

    std::printf("\n%s: %d fallas\n", fallos ? "FALLA" : "PASA", fallos);
    return fallos ? 1 : 0;
}