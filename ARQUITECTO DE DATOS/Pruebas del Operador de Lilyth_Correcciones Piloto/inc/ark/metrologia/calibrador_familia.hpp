// calibrador_familia.hpp -- calibracion de FAMILIA de cualquier operador delta (fase 2 de la calibracion de h), C++20.
// Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-05. Reglas aprobadas por el Comandante sobre el dictamen del auditor externo Argon
// (mensajes/RESPUESTA_ARGON_A_REVISION_CODE_CALIBRACION_H.md); prerregistro de la generalizacion:
// JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/calibracion_familia/EXPECTATIVAS_FAMILIA_GENERICA.md.
//
// Lo generico, extraido sin cambiar una operacion del operador de coordinacion (que se rehace sobre esto y da sus mismos numeros):
//   - la falsa alarma de la familia es un EPISODIO de coalescencia de los miembros de una subfamilia: abre cuando alguno abre, cierra
//     cuando todos cerraron, y la union vuelve a inactivo al empezar cada realizacion (tramo) del nulo juntado;
//   - m_min = ceil((z_{1-alfa/2}/ln 2)^2) (intervalo de factor 2); B = min(B_max, max(1, ceil(m_min ARL0_sub / N)));
//     si B queda topado, el objetivo es arl0_max_certificable = B N / m_min y la ficha lo marca;
//   - A* (rareza por canal, en barras) se MIDE: acotamiento por factor 4 desde el A de Bonferroni y biseccion en log A hasta que los
//     extremos difieran en <= 1 episodio (tope de iteraciones); A* = el extremo que cumple.
// Cada operador pone lo suyo: como construye sus B realizaciones, como fija la rareza A en sus miembros y como corre su maquina de
// estados sobre el nulo juntado (la funcion evaluar(A) -> episodios). Fuera del ciclo caliente.
#pragma once

#include "ark/metrologia/estadistica_metrologia.hpp"
#include "ark/metrologia/tipos_metrologia.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <mutex>
#include <span>
#include <thread>
#include <vector>

namespace ark::metrologia {

enum class PresupuestoFamilia : std::uint8_t { BONFERRONI_CANALES = 0, FAMILIA_EPISODIOS = 1 };
// RAMA (2026-11-14, paso 2; decisiones del Comandante: se bifurca, presupuesto SEPARADO, un delta con dos ARL0):
//   CONJUNTO   (omision) el episodio de siempre: abre con un evento de nivel o de temperatura y cierra cuando no queda ninguno abierto.
//   POR_CANAL  con FAMILIA_EPISODIOS, cada canal (nivel, temperatura) es su propia subfamilia con presupuesto 2 ARL0 (Bonferroni sobre 2)
//              y su propia rareza A*, calibradas en alterna hasta que ninguna cambie. Con BONFERRONI_CANALES los umbrales son los mismos.
enum class DefinicionEpisodio : std::uint8_t { CONJUNTO = 0, POR_CANAL = 1 };
// Canal de un episodio: el de nivel (abierto o con un evento de canal 0 emitido en el paso), el de temperatura, o los dos.
enum class CanalEpisodio : std::uint8_t { CONJUNTO = 0, NIVEL = 1, TEMPERATURA = 2 };
template <class SalidaDelta>
inline bool actividad_canal(const SalidaDelta& s, CanalEpisodio k) noexcept {
    if (k == CanalEpisodio::CONJUNTO) return s.evento_abierto || s.n_eventos > 0;
    const int canal = k == CanalEpisodio::NIVEL ? 0 : 1;
    if (canal == 0 ? s.nivel_abierto : s.temperatura_abierto) return true;
    for (std::size_t i = 0; i < s.n_eventos; ++i)
        if (s.eventos[i].canal == canal) return true;
    return false;
}

struct FichaSubfamilia {
    bool presente = false;                       // tiene al menos un miembro calibrado
    std::size_t componentes = 0, canales = 0;    // miembros; canales = 2 suma J_c
    std::size_t B = 0, m_min = 0;
    bool b_topado = false;                       // B quedo en B_max
    std::size_t barras_nulo = 0, pasos_nulo = 0; // B N barras; pasos evaluados (sin las uniones)
    double arl0_sub_barras = NaN;                // (subfamilias) x ARL0_familia
    double arl0_max_certificable = NaN;          // B N / m_min
    double arl0_objetivo_barras = NaN;           // min(ARL0_sub, arl0_max_certificable)
    double episodios_maximos = NaN;              // B N / arl0_objetivo
    bool certificado = false;                    // B N >= m_min ARL0_sub y la biseccion encontro un extremo que cumple
    bool acotado_por_historia = false;           // el objetivo es arl0_max_certificable, no ARL0_sub
    bool sin_cota = false;                       // ni con el A mas alto probado se cumplio
    bool en_suelo = false;                       // se cumplio aun con el A mas bajo admisible
    double a_bonferroni = NaN;                   // ARL0_sub x canales (punto de partida)
    double a_estrella = NaN;                     // A por canal en barras (el extremo que cumple)
    double a_otro_extremo = NaN;
    std::size_t episodios_estrella = 0, episodios_otro_extremo = 0;
    std::size_t evaluaciones = 0, iteraciones_biseccion = 0;
    double separacion_rotacion_barras = NaN;     // N / (B + 1), si el nulo es una ROTACION
};

struct ParametrosFamilia {
    std::size_t realizaciones_max = 10;     // B_max (presupuesto de computo declarado)
    double alfa_certificacion = 0.05;       // alfa del proyecto, fija m_min
    std::size_t iteraciones_max = 12;       // tope de la biseccion
    std::size_t expansiones_max = 10;       // tope del acotamiento (factor 4 cada una)
    std::size_t hilos = 0;                  // 0: hardware_concurrency
};

namespace familia {

// m_min = ceil((z_{1-alfa/2} / ln 2)^2)
inline std::size_t m_min(double alfa) {
    const double z = est::norm_ppf(1.0 - alfa / 2.0);
    const double r = z / std::log(2.0);
    return static_cast<std::size_t>(std::ceil(r * r));
}

// B = min(B_max, max(1, ceil(m_min ARL0_sub / N))); b_pedido sale sin topar (para la marca b_topado).
inline std::size_t realizaciones(std::size_t mm, double arl0_sub, std::size_t n, std::size_t b_max, double& b_pedido) {
    b_pedido = std::ceil(static_cast<double>(mm) * arl0_sub / static_cast<double>(n));
    return std::min(b_max, static_cast<std::size_t>(std::max(1.0, b_pedido)));
}

// La ficha de una subfamilia antes de la busqueda (todo lo que no depende de A).
inline void preparar_ficha(FichaSubfamilia& f, std::size_t B, std::size_t mm, double b_pedido, std::size_t b_max, std::size_t n,
                           std::size_t pasos_por_realizacion, double arl0_sub) {
    f.B = B;
    f.m_min = mm;
    f.b_topado = b_pedido > static_cast<double>(b_max);
    f.barras_nulo = B * n;
    f.pasos_nulo = B * pasos_por_realizacion;
    f.arl0_sub_barras = arl0_sub;
    f.arl0_max_certificable = static_cast<double>(B * n) / static_cast<double>(mm);
    f.arl0_objetivo_barras = std::min(arl0_sub, f.arl0_max_certificable);
    f.acotado_por_historia = f.arl0_max_certificable < arl0_sub;
    f.episodios_maximos = static_cast<double>(B * n) / f.arl0_objetivo_barras;
}

inline std::size_t n_hilos(std::size_t pedidos, std::size_t tareas) noexcept {
    std::size_t h = pedidos ? pedidos : static_cast<std::size_t>(std::thread::hardware_concurrency());
    if (h == 0) h = 1;
    return std::min(h, std::max<std::size_t>(tareas, 1));
}

// f(i) para i en [0, n), repartido entre hilos; la primera excepcion se relanza al final.
template <class F>
void en_paralelo(std::size_t n, std::size_t hilos, F&& f) {
    const std::size_t h = n_hilos(hilos, n);
    if (h <= 1) { for (std::size_t i = 0; i < n; ++i) f(i); return; }
    std::atomic<std::size_t> sig{0};
    std::exception_ptr err;
    std::mutex m;
    auto trab = [&] {
        for (;;) {
            const std::size_t i = sig.fetch_add(1);
            if (i >= n) return;
            try { f(i); } catch (...) { std::lock_guard<std::mutex> l(m); if (!err) err = std::current_exception(); }
        }
    };
    std::vector<std::thread> th;
    th.reserve(h - 1);
    for (std::size_t k = 1; k < h; ++k) th.emplace_back(trab);
    trab();
    for (auto& t : th) t.join();
    if (err) std::rethrow_exception(err);
}

}  // namespace familia

// Episodios de coalescencia: union de la actividad de los componentes miembro (act[c][p], miembro[c] != 0), un episodio por cada
// paso de inactivo a activo dentro de cada tramo [ini[r], fin[r]); la union vuelve a inactivo al empezar cada tramo.
inline std::size_t episodios_coalescencia(const std::vector<std::vector<std::uint8_t>>& act, std::span<const std::uint8_t> miembro,
                                          std::span<const std::size_t> ini, std::span<const std::size_t> fin) {
    std::size_t e = 0;
    for (std::size_t r = 0; r < ini.size(); ++r) {
        bool abierto = false;
        for (std::size_t p = ini[r]; p < fin[r]; ++p) {
            bool u = false;
            for (std::size_t c = 0; c < miembro.size() && !u; ++c)
                if (miembro[c] && act[c][p]) u = true;
            if (u && !abierto) ++e;
            abierto = u;
        }
    }
    return e;
}

namespace familia {

// Busca A* de una subfamilia con la ficha ya preparada (a_bonferroni y episodios_maximos). evaluar(A) fija la rareza A en los miembros,
// corre el nulo juntado y devuelve los episodios. suelo: el A mas bajo admisible. No deja los miembros en A*: eso lo hace el operador.
template <class Evaluar>
void buscar_a_estrella(FichaSubfamilia& f, double suelo, const ParametrosFamilia& p, Evaluar&& evaluar) {
    auto eval = [&](double A) { ++f.evaluaciones; return evaluar(A); };
    const double emax = f.episodios_maximos;
    double lo = NaN, hi = NaN;
    std::size_t elo = 0, ehi = 0;
    const double A0 = std::max(f.a_bonferroni, suelo);
    const std::size_t e0 = eval(A0);
    if (static_cast<double>(e0) <= emax) {
        hi = A0; ehi = e0;
        for (std::size_t k = 0; k < p.expansiones_max; ++k) {
            if (hi <= suelo) { f.en_suelo = true; break; }
            const double A = std::max(hi / 4.0, suelo);
            const std::size_t e = eval(A);
            if (static_cast<double>(e) > emax) { lo = A; elo = e; break; }
            hi = A; ehi = e;
        }
        if (!std::isfinite(lo) && hi <= suelo) f.en_suelo = true;
    } else {
        lo = A0; elo = e0;
        for (std::size_t k = 0; k < p.expansiones_max; ++k) {
            const double A = lo * 4.0;
            const std::size_t e = eval(A);
            if (static_cast<double>(e) <= emax) { hi = A; ehi = e; break; }
            lo = A; elo = e;
        }
        if (!std::isfinite(hi)) f.sin_cota = true;
    }
    if (std::isfinite(lo) && std::isfinite(hi)) {
        while (elo > ehi + 1 && f.iteraciones_biseccion < p.iteraciones_max) {
            const double mid = std::sqrt(lo * hi);
            const std::size_t e = eval(mid);
            if (static_cast<double>(e) <= emax) { hi = mid; ehi = e; } else { lo = mid; elo = e; }
            ++f.iteraciones_biseccion;
        }
    }
    if (f.sin_cota) { f.a_estrella = lo; f.episodios_estrella = elo; }
    else {
        f.a_estrella = hi; f.episodios_estrella = ehi;
        if (std::isfinite(lo)) { f.a_otro_extremo = lo; f.episodios_otro_extremo = elo; }
    }
    f.certificado = !f.sin_cota && !f.acotado_por_historia;
}

// Indices de una permutacion por BLOQUES de longitud b (el mismo algoritmo que sustitutos.hpp: orden de los bloques barajado, el
// resto al final). Con menos de 2 bloques, la identidad.
template <class Rng>
std::vector<std::size_t> permutacion_bloques(std::size_t n, std::size_t b, Rng& rng) {
    std::vector<std::size_t> idx(n);
    for (std::size_t i = 0; i < n; ++i) idx[i] = i;
    b = std::max<std::size_t>(1, b);
    const std::size_t nb = n / b;
    if (nb < 2) return idx;
    std::vector<std::size_t> orden(nb);
    for (std::size_t i = 0; i < nb; ++i) orden[i] = i;
    std::shuffle(orden.begin(), orden.end(), rng);
    for (std::size_t i = 0; i < nb; ++i)
        for (std::size_t k = 0; k < b; ++k) idx[i * b + k] = orden[i] * b + k;
    return idx;
}

}  // namespace familia

}  // namespace ark::metrologia
