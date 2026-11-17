// sensor_coordinacion.hpp -- sensor de coordinacion de un PAR: theta^ cada S barras sobre una ventana de W barras de los dos
// activos. Jefe de Arquitectos de Codigo (Opus 5.5), 2026-11-05, diseno acordado con el Comandante. Prerregistro:
// JEFE DE ARQUITECTOS DE CODIGO/PRUEBAS_JEFE/coordinacion_sensor/EXPECTATIVAS_COORDINACION.md.
//
// Escalas: s = 2^j que dividen W, con s <= W / max_ratio_escalas para j >= 1 (el maximo de subdivisiones del modulo de ventanas,
//   ConfiguracionVentanaCausalCxx::max_ratio_escalas) y n_j = W / s >= k + 3 (KSG definido: pares n_j - 1 >= k + 2). k es ENTRADA;
//   el conjunto admisible sale de W y k. Bloques de s barras alineados con el borde presente (W divisible: el primero empieza en 0).
// Observables por escala (de los retornos base r del bloque):
//   DIRECCION             R_j = suma de r
//   MAGNITUD              RV_j = suma de r^2
//   DIRECCION_NORMALIZADA R_j / sqrt(RV_j) (0 si RV_j = 0): herramienta opcional para el confundido de la volatilidad
// Cantidades por escala, observable y sentido A -> B, con PIT aleatorizado por ventana y pasado de una barra:
//   TE  = I(A_p; B_f | B_p)  (Frenzel-Pompe)
//   MI  = I(A_p; B_f)        (KSG-1)
//   SIN = TE - MI            (sinergia neta = - informacion de interaccion de A_p, B_p, B_f)
// Cuatro busquedas de vecinos por observable y escala. Componente c = ((o J + j) 2 + sentido) 3 + cantidad.
// Nulos (fuera de este sensor, en el observador): SIGNOS de centro CERO sobre los retornos base (conserva RV_j exacto en todas las
//   escalas) para DIRECCION y NORMALIZADA; ROTACION de Y n/2 contra X para MAGNITUD.
// Causalidad: cada theta^ usa solo las ultimas W barras recibidas; el primero sale con W barras y despues cada S. Ciclo caliente
// (paso): cero reservas.
#pragma once

#include "ark/coordinacion/ksg_directo.hpp"

#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace ark::coordinacion {

enum ObservableCoord : std::uint32_t { OBS_DIRECCION = 1u, OBS_MAGNITUD = 2u, OBS_DIRECCION_NORMALIZADA = 4u };
inline constexpr std::array<std::uint32_t, 3> OBS_TODOS{OBS_DIRECCION, OBS_MAGNITUD, OBS_DIRECCION_NORMALIZADA};
inline const char* nombre_observable(std::uint32_t o) {
    return o == OBS_DIRECCION ? "DIRECCION" : o == OBS_MAGNITUD ? "MAGNITUD" : o == OBS_DIRECCION_NORMALIZADA ? "DIRECCION_NORMALIZADA" : "?";
}
enum class Cantidad : std::uint8_t { TE = 0, MI = 1, SIN = 2 };
inline constexpr std::array<const char*, 3> NOMBRE_CANTIDAD{"TE", "MI", "SIN"};
enum class Sentido : std::uint8_t { X_A_Y = 0, Y_A_X = 1 };
inline constexpr std::array<const char*, 2> NOMBRE_SENTIDO{"X->Y", "Y->X"};

struct ConfiguracionSensorCoordinacion {
    std::size_t paso_sensor = 0;          // OBLIGATORIO >= 1: S, un theta^ cada S barras
    std::size_t k_vecinos = 0;            // OBLIGATORIO >= 1: k del KSG (entrada; las escalas admisibles salen de W y k)
    double max_ratio_escalas = 4.0;       // el del modulo de ventanas (s <= W / max_ratio)
    std::uint32_t observables = OBS_DIRECCION | OBS_MAGNITUD;
    std::uint64_t semilla = 20261105;     // semilla del PIT aleatorizado
};

struct DescriptorComponente { std::uint32_t observable = 0; std::size_t j = 0, s = 1; Sentido sentido = Sentido::X_A_Y; Cantidad cantidad = Cantidad::TE; };

// Escalas admisibles (s = 2^j) para W, k y max_ratio. Vacio si ni la nativa cabe.
inline std::vector<std::size_t> escalas_admisibles_coord(std::size_t W, std::size_t k, double max_ratio) {
    std::vector<std::size_t> s;
    for (std::size_t j = 0; j < 63; ++j) {
        const std::size_t b = std::size_t{1} << j;
        if (b > W || W % b != 0) break;
        if (j >= 1 && max_ratio > 0.0 && static_cast<double>(b) > static_cast<double>(W) / max_ratio) break;
        if (W / b < k + 3) break;
        s.push_back(b);
    }
    return s;
}

template <std::size_t W>
class SensorCoordinacion {
    static_assert(W >= 8, "W >= 8");
public:
    static constexpr std::size_t J_CAP = static_cast<std::size_t>(std::countr_zero(W)) + 1;
    static constexpr std::size_t C_MAX = 3 * J_CAP * 6;
    struct Resultado {
        bool valido = false;
        std::int64_t paso = -1;           // indice del theta^ desde el reinicio
        std::int64_t barra = -1;          // ultima barra de la ventana (indice desde el reinicio)
        std::array<double, C_MAX> theta{};
    };

private:
    ConfiguracionSensorCoordinacion cfg_{};
    std::array<std::size_t, J_CAP> s_{};
    std::size_t J_ = 0;
    std::array<std::uint32_t, 3> obs_{};
    std::size_t n_obs_ = 0;
    bool listo_ = false;
    // anillo de las ultimas W barras
    std::array<double, W> ax_{}, ay_{};
    std::array<std::uint8_t, W> aok_{};
    std::size_t cab_ = 0;
    std::int64_t barras_ = 0, pasos_ = 0;
    // trabajo
    std::array<double, W> x_{}, y_{}, vx_{}, vy_{}, ux_{}, uy_{};
    std::array<std::size_t, W> idx_{};
    KSGDirecto<W> ksg_{};

public:
    void configurar(const ConfiguracionSensorCoordinacion& cfg) {
        if (cfg.paso_sensor < 1) throw std::invalid_argument("ConfiguracionSensorCoordinacion.paso_sensor >= 1 es obligatorio");
        if (cfg.k_vecinos < 1) throw std::invalid_argument("ConfiguracionSensorCoordinacion.k_vecinos >= 1 es obligatorio");
        cfg_ = cfg;
        const auto s = escalas_admisibles_coord(W, cfg.k_vecinos, cfg.max_ratio_escalas);
        if (s.empty()) throw std::invalid_argument("sensor de coordinacion: ninguna escala admisible para W y k");
        J_ = s.size() < J_CAP ? s.size() : J_CAP;
        for (std::size_t j = 0; j < J_; ++j) s_[j] = s[j];
        n_obs_ = 0;
        for (std::uint32_t o : OBS_TODOS) if (cfg.observables & o) obs_[n_obs_++] = o;
        if (n_obs_ == 0) throw std::invalid_argument("sensor de coordinacion: ningun observable");
        listo_ = true;
        reiniciar();
    }
    void reiniciar() noexcept { cab_ = 0; barras_ = 0; pasos_ = 0; }

    [[nodiscard]] bool listo() const noexcept { return listo_; }
    [[nodiscard]] const ConfiguracionSensorCoordinacion& configuracion() const noexcept { return cfg_; }
    [[nodiscard]] std::size_t escalas() const noexcept { return J_; }
    [[nodiscard]] std::size_t escala(std::size_t j) const noexcept { return s_[j]; }
    [[nodiscard]] std::size_t observables() const noexcept { return n_obs_; }
    [[nodiscard]] std::uint32_t observable(std::size_t o) const noexcept { return obs_[o]; }
    [[nodiscard]] std::size_t componentes() const noexcept { return n_obs_ * J_ * 6; }
    [[nodiscard]] std::size_t ventana() const noexcept { return W; }
    [[nodiscard]] std::size_t paso_sensor() const noexcept { return cfg_.paso_sensor; }
    [[nodiscard]] static constexpr std::size_t indice(std::size_t o, std::size_t J, std::size_t j, Sentido d, Cantidad q) noexcept {
        return ((o * J + j) * 2 + static_cast<std::size_t>(d)) * 3 + static_cast<std::size_t>(q);
    }
    [[nodiscard]] DescriptorComponente descriptor(std::size_t c) const noexcept {
        DescriptorComponente d{};
        d.cantidad = static_cast<Cantidad>(c % 3);
        d.sentido = static_cast<Sentido>((c / 3) % 2);
        d.j = (c / 6) % J_;
        d.s = s_[d.j];
        d.observable = obs_[c / (6 * J_)];
        return d;
    }
    [[nodiscard]] std::string nombre(std::size_t c) const {
        const auto d = descriptor(c);
        return std::string(nombre_observable(d.observable)) + " s" + std::to_string(d.s) + " " + NOMBRE_SENTIDO[static_cast<int>(d.sentido)] + " " +
               NOMBRE_CANTIDAD[static_cast<int>(d.cantidad)];
    }

    // theta^ de una ventana cronologica x[0..W-1], y[0..W-1] (ya validada). `paso` siembra el PIT. Cero reservas.
    void theta_ventana(const double* x, const double* y, std::int64_t paso, Resultado& r) noexcept {
        const std::size_t k = cfg_.k_vecinos;
        for (std::size_t j = 0; j < J_; ++j) {
            const std::size_t s = s_[j], n = W / s, m = n - 1;
            for (std::size_t o = 0; o < n_obs_; ++o) {
                for (std::size_t b = 0; b < n; ++b) {
                    double sx = 0.0, sy = 0.0, qx = 0.0, qy = 0.0;
                    for (std::size_t i = b * s; i < (b + 1) * s; ++i) { sx += x[i]; sy += y[i]; qx += x[i] * x[i]; qy += y[i] * y[i]; }
                    if (obs_[o] == OBS_DIRECCION) { vx_[b] = sx; vy_[b] = sy; }
                    else if (obs_[o] == OBS_MAGNITUD) { vx_[b] = qx; vy_[b] = qy; }
                    else { vx_[b] = qx > 0.0 ? sx / std::sqrt(qx) : 0.0; vy_[b] = qy > 0.0 ? sy / std::sqrt(qy) : 0.0; }
                }
                GeneradorU01 gx(semilla_mezcla(cfg_.semilla, static_cast<std::uint64_t>(paso), j, obs_[o], 0));
                GeneradorU01 gy(semilla_mezcla(cfg_.semilla, static_cast<std::uint64_t>(paso), j, obs_[o], 1));
                pit_aleatorizado(vx_.data(), n, ux_.data(), idx_.data(), gx);
                pit_aleatorizado(vy_.data(), n, uy_.data(), idx_.data(), gy);
                const double* xp = ux_.data();
                const double* xf = ux_.data() + 1;
                const double* yp = uy_.data();
                const double* yf = uy_.data() + 1;
                const double te_xy = ksg_.cmi(xp, yf, yp, m, k), mi_xy = ksg_.mi(xp, yf, m, k);
                const double te_yx = ksg_.cmi(yp, xf, xp, m, k), mi_yx = ksg_.mi(yp, xf, m, k);
                r.theta[indice(o, J_, j, Sentido::X_A_Y, Cantidad::TE)] = te_xy;
                r.theta[indice(o, J_, j, Sentido::X_A_Y, Cantidad::MI)] = mi_xy;
                r.theta[indice(o, J_, j, Sentido::X_A_Y, Cantidad::SIN)] = te_xy - mi_xy;
                r.theta[indice(o, J_, j, Sentido::Y_A_X, Cantidad::TE)] = te_yx;
                r.theta[indice(o, J_, j, Sentido::Y_A_X, Cantidad::MI)] = mi_yx;
                r.theta[indice(o, J_, j, Sentido::Y_A_X, Cantidad::SIN)] = te_yx - mi_yx;
            }
        }
    }

    // Ciclo caliente: una barra del par. Devuelve theta^ cuando toca (valido = false si la ventana tiene barras invalidas).
    std::optional<Resultado> paso(double rx, double ry, bool valido) noexcept {
        if (!listo_) return std::nullopt;
        const bool ok = valido && std::isfinite(rx) && std::isfinite(ry);
        ax_[cab_] = ok ? rx : 0.0;
        ay_[cab_] = ok ? ry : 0.0;
        aok_[cab_] = ok ? 1 : 0;
        cab_ = (cab_ + 1) % W;
        ++barras_;
        if (barras_ < static_cast<std::int64_t>(W) || (barras_ - static_cast<std::int64_t>(W)) % static_cast<std::int64_t>(cfg_.paso_sensor) != 0)
            return std::nullopt;
        Resultado r{};
        r.paso = pasos_++;
        r.barra = barras_ - 1;
        bool todo = true;
        for (std::size_t i = 0; i < W; ++i) {
            const std::size_t q = (cab_ + i) % W;   // la mas antigua primero
            x_[i] = ax_[q];
            y_[i] = ay_[q];
            todo = todo && aok_[q];
        }
        if (!todo) { r.theta.fill(NaN_C); return r; }
        r.valido = true;
        theta_ventana(x_.data(), y_.data(), r.paso, r);
        for (std::size_t c = componentes(); c < C_MAX; ++c) r.theta[c] = NaN_C;
        return r;
    }

    // Serie de theta^ sobre dos series de barras (reinicia el sensor). Fuera del ciclo caliente.
    struct Serie {
        std::vector<std::vector<double>> theta;   // [componente][paso]
        std::vector<std::uint8_t> ok;              // [paso]
        std::vector<std::int64_t> barra;           // [paso]: ultima barra de la ventana
    };
    Serie serie(std::span<const double> x, std::span<const double> y, std::span<const std::uint8_t> ok) {
        if (x.size() != y.size() || (!ok.empty() && ok.size() != x.size())) throw std::invalid_argument("sensor de coordinacion: longitudes");
        Serie s;
        s.theta.assign(componentes(), {});
        reiniciar();
        for (std::size_t i = 0; i < x.size(); ++i) {
            const auto r = paso(x[i], y[i], ok.empty() || ok[i]);
            if (!r) continue;
            for (std::size_t c = 0; c < componentes(); ++c) s.theta[c].push_back(r->theta[c]);
            s.ok.push_back(r->valido ? 1 : 0);
            s.barra.push_back(r->barra);
        }
        reiniciar();
        return s;
    }
};

// ---------------- nulos de la coordinacion (fuera del ciclo caliente) ----------------

// SIGNOS de centro CERO: cada retorno por un +-1 independiente (U < 1/2). Conserva |r| (y con el RV_j de todas las escalas) y el
// co-movimiento de las magnitudes; destruye la dependencia de signo. DIVERGE de sustitutos.hpp (SIGNOS centrado en la media), que no
// conserva |r| cuando la media no es 0.
inline std::vector<double> nulo_signos_cero(std::span<const double> r, GeneradorU01& g) {
    std::vector<double> y(r.begin(), r.end());
    for (auto& v : y) v = g.u01() < 0.5 ? -v : v;
    return y;
}

// ROTACION: y_rot[i] = y[(i + d) mod n], con d = n / 2 por omision (la mayor distancia circular a los dos extremos; declarado, no
// ajustado). Con B realizaciones (calibracion de familia) d_b = floor(b n / (B + 1)), b = 1..B.
template <class T>
inline std::vector<T> nulo_rotacion(std::span<const T> y, std::size_t d) {
    const std::size_t n = y.size();
    std::vector<T> o(n);
    for (std::size_t i = 0; i < n; ++i) o[i] = y[(i + d) % n];
    return o;
}
template <class T>
inline std::vector<T> nulo_rotacion(std::span<const T> y) { return nulo_rotacion(y, y.size() / 2); }

}  // namespace ark::coordinacion
