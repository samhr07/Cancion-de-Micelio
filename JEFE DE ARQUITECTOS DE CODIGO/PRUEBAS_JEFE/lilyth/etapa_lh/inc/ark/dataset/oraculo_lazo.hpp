// oraculo_lazo.hpp -- calca C++ de ARQUITECTO DE DATOS/filtros_de_ventana/oraculo_lazo.py (LazoOraculo), el lazo de Shannon que
// usa motor_ventanas._LazoFeedbackShannon. Misma especificacion, mismo orden de operaciones (Neumaier, bisect, mediana por
// seleccion de numpy, pow/log/exp de la CRT), W en tiempo de ejecucion. Tablas: tablas_lazo_ventanas.hpp (generadas del JSON del
// principal). Jefe de Arquitectos, 2026-11-03.
#pragma once

#include "ark/dataset/tablas_lazo_ventanas.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <vector>

namespace ark::dataset::ventanas {

namespace oraculo {

inline constexpr double LN2 = 0.693147180559945309417;
inline constexpr double Z1A = 1.6448536269514722;
inline constexpr double C_MAD = 0.6744897501960817;
inline constexpr int B = 8;
inline constexpr int H = 64;
inline constexpr size_t EV_MIN_OBS = 40;
inline constexpr int64_t EV_MIN_COMPLETOS = 20;
inline constexpr double NAN_ = std::numeric_limits<double>::quiet_NaN();
inline constexpr double INF = std::numeric_limits<double>::infinity();

inline double psi1(double x) {
    double s = 0.0;
    while (x < 10.0) {
        s += 1.0 / (x * x);
        x += 1.0;
    }
    const double x2 = 1.0 / (x * x);
    return s + 1.0 / x + x2 / 2.0 +
           (x2 / x) * (1.0 / 6.0 + x2 * (-1.0 / 30.0 + x2 * (1.0 / 42.0 + x2 * (-1.0 / 30.0 + x2 * (5.0 / 66.0 + x2 * (-691.0 / 2730.0))))));
}

inline double psi1_inv(double y) {
    double lo = -30.0, hi = 30.0;
    for (int i = 0; i < 100; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (psi1(std::exp(mid)) > y) lo = mid;
        else hi = mid;
    }
    return std::exp(0.5 * (lo + hi));
}

struct Neumaier {
    double s = 0.0, c = 0.0;
    void add(double x) {
        const double t = s + x;
        if (std::fabs(s) >= std::fabs(x)) c += (s - t) + x;
        else c += (x - t) + s;
        s = t;
    }
    double val() const { return s + c; }
};

inline double mediana_ordenada(const std::vector<double>& v) {
    const size_t k = v.size();
    return (k % 2) ? v[k / 2] : 0.5 * (v[k / 2 - 1] + v[k / 2]);
}

template <class It>
inline double suma_seq(It a, It b) {
    Neumaier n;
    for (; a != b; ++a) n.add(*a);
    return n.val();
}
inline double suma_seq(const std::vector<double>& v) { return suma_seq(v.begin(), v.end()); }

inline double q95_nulo(int64_t N) {
    const auto& QN = tablas_lazo::REJILLA;
    const auto& Q95 = tablas_lazo::QTAB;
    if (N > QN.back()) return std::exp(tablas_lazo::MU_INF + Z1A * std::sqrt(tablas_lazo::V_ASINTOTICA / static_cast<double>(N)));
    const size_t j = static_cast<size_t>(std::lower_bound(QN.begin(), QN.end(), N) - QN.begin());
    if (QN[j] == N) return Q95[j];
    const size_t a = j - 1, b = j;
    const double w = (std::log(static_cast<double>(N)) - std::log(static_cast<double>(QN[a]))) /
                     (std::log(static_cast<double>(QN[b])) - std::log(static_cast<double>(QN[a])));
    return std::exp(std::log(Q95[a]) + w * (std::log(Q95[b]) - std::log(Q95[a])));
}

// np.median sobre valores no NaN: seleccion (el k-esimo; par = media de los dos centrales como np.mean de dos).
inline double np_median_sel(std::vector<double> v) {
    const size_t n = v.size();
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(n / 2), v.end());
    const double a = v[n / 2];
    if (n % 2) return a;
    const double b = *std::max_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(n / 2));
    return (0.0 + b + a) / 2.0;
}

}  // namespace oraculo

// El dict `p` de calcular(): claves presentes y sus valores.
struct PropuestaLazo {
    double d_nu = 0.0;
    bool invalido = false, primera = false;
    bool tiene_r = false, tiene_sh = false, tiene_sigma = false;
    double r = 0.0, sh = 0.0, sigma = 0.0;
};

class LazoOraculo {
public:
    // Atributos con el mismo nombre que en Python (los lee _LazoFeedbackShannon).
    int64_t W, s_min;
    double tick;
    size_t CAP;
    bool fuga;
    bool ini = false;
    double val_prev = oraculo::NAN_;
    std::deque<double> win;
    int64_t len = 0;
    int64_t N = 0;
    double mean = 0.0, M2 = 0.0;
    std::deque<double> hist;
    std::vector<double> ord, ord_abs;
    std::vector<std::pair<int64_t, double>> activos;
    int64_t k = 0;
    std::vector<oraculo::Neumaier> sum_rat;
    std::vector<int64_t> cnt;
    oraculo::Neumaier sum_d0;
    int64_t cnt_d0 = 0, completos = 0;
    std::vector<double> bloque;
    oraculo::Neumaier SG, SG2, SL, SL2, SQ;
    int64_t mG = 0, mL = 0, mQ = 0;
    double nu0 = oraculo::INF, alfa0 = oraculo::INF;
    int64_t uN = 0;
    oraculo::Neumaier S1, S2, P1, P10;
    std::vector<double> primeros;
    std::deque<double> ultimos;
    bool phi_ok = false;
    double w_phi = oraculo::NAN_, phi = oraculo::NAN_;
    int64_t cnt_w_max = -1;
    bool laplace = false;
    int64_t w_max = 8, w_innov = 8;
    double sh = oraculo::NAN_;
    int64_t tau_v = 0;
    double umbral = oraculo::NAN_;

    explicit LazoOraculo(int64_t W_ = 32, int64_t s_min_ = 4, double tick_ = 0.01, size_t CAP_ = 4096, bool fuga_ = false)
        : W(W_), s_min(s_min_), tick(tick_), CAP(CAP_), fuga(fuga_), sum_rat(oraculo::H + 1), cnt(oraculo::H + 1, 0) {}

    PropuestaLazo calcular(double val, bool mask) const {
        using namespace oraculo;
        PropuestaLazo p;
        if (!(std::isfinite(val) && val > 0.0)) { p.invalido = true; return p; }
        if (!ini) { p.primera = true; return p; }
        const double r = std::log(val) - std::log(val_prev);
        if (mask) { p.tiene_r = true; p.r = 0.0; return p; }
        const double qs = tick / val_prev;
        const double shh = qs * qs / 6.0;
        p.tiene_r = true; p.r = r; p.tiene_sh = true; p.sh = shh;
        const int64_t n = len;
        const int64_t kmin = std::min(std::max<int64_t>(2 * s_min, 8), w_innov);
        if (std::fabs(r) < 1e-12 || n < kmin) return p;
        std::vector<double> x(win.end() - n, win.end());
        if (fuga) x.push_back(r);
        const double nn = static_cast<double>(x.size());
        if (!laplace) {
            const double mu = suma_seq(x) / nn;
            std::vector<double> d(x.size());
            for (size_t i = 0; i < x.size(); ++i) d[i] = (x[i] - mu) * (x[i] - mu);
            const double v = suma_seq(d) / nn;
            const double s02 = std::max(M2 / static_cast<double>(N), shh);   // max(a, b) de Python: a salvo b > a
            const double s2 = std::isinf(nu0) ? s02 : (nu0 * s02 + nn * v) / (nu0 + nn - 1.0);
            const double s = std::sqrt(py_max(s2, shh));
            const double z = (r - mu) / s;
            p.d_nu = z * z / (2.0 * LN2);
            p.tiene_sigma = true; p.sigma = s;
        } else {
            std::vector<double> xs = x;
            std::sort(xs.begin(), xs.end());
            const double m = mediana_ordenada(xs);
            std::vector<double> a(x.size());
            for (size_t i = 0; i < x.size(); ++i) a[i] = std::fabs(x[i] - m);
            const double bh = suma_seq(a) / nn;
            const double bfl = std::sqrt(shh / 2.0);
            const double b0v = py_max(b0(), bfl);
            double b = std::isinf(alfa0) ? b0v : (alfa0 * b0v + nn * bh) / (alfa0 + nn);
            b = py_max(b, bfl);
            p.d_nu = std::fabs((r - m) / b) / (2.0 * LN2);
            p.tiene_sigma = true; p.sigma = std::sqrt(2.0) * b;
        }
        return p;
    }

    double b0() const {
        const double m = oraculo::mediana_ordenada(ord);
        std::vector<double> a(ord.size());
        for (size_t i = 0; i < ord.size(); ++i) a[i] = std::fabs(ord[i] - m);
        return oraculo::suma_seq(a) / static_cast<double>(ord.size());
    }

    void confirmar(double val, bool mask, const PropuestaLazo& p) {
        using namespace oraculo;
        if (p.invalido) return;
        if (p.primera) { val_prev = val; ini = true; return; }
        if (mask) { val_prev = val; return; }
        const double r = p.r, shh = p.sh;
        sh = shh;
        const double u = std::fabs(r);
        const double s0_pre = N > 0 ? M2 / static_cast<double>(N) : 0.0;
        // a. evento
        bool evento = false;
        if (ord_abs.size() >= EV_MIN_OBS && s0_pre > 0.0) {
            const int64_t kq = static_cast<int64_t>(ord_abs.size());
            const int64_t c = static_cast<int64_t>(std::ceil(0.975 * static_cast<double>(kq)));
            const double q = ord_abs[static_cast<size_t>(std::min(kq - 1, c - 1))];
            evento = u >= q;
        }
        // b. perfil
        std::vector<std::pair<int64_t, double>> quedan;
        for (const auto& [i0, s0e] : activos) {
            const int64_t h = k - i0;
            sum_rat[static_cast<size_t>(h)].add(r * r / s0e);
            cnt[static_cast<size_t>(h)] += 1;
            if (h == H) completos += 1;
            else quedan.push_back({i0, s0e});
        }
        activos = std::move(quedan);
        // c. alta
        if (evento) {
            activos.push_back({k, s0_pre});
            sum_d0.add(r * r / s0_pre);
            cnt_d0 += 1;
        }
        // d. Welford
        N += 1;
        const double d = r - mean;
        mean += d / static_cast<double>(N);
        M2 += d * (r - mean);
        // e. historia larga
        hist.push_back(r);
        ord.insert(std::upper_bound(ord.begin(), ord.end(), r), r);
        ord_abs.insert(std::upper_bound(ord_abs.begin(), ord_abs.end(), u), u);
        if (hist.size() > CAP) {
            const double viejo = hist.front();
            hist.pop_front();
            ord.erase(std::lower_bound(ord.begin(), ord.end(), viejo));
            const double av = std::fabs(viejo);
            ord_abs.erase(std::lower_bound(ord_abs.begin(), ord_abs.end(), av));
        }
        // f. bloques
        bloque.push_back(r);
        bool bloque_nuevo = false;
        if (static_cast<int>(bloque.size()) == B) {
            bloque_nuevo = true;
            const std::vector<double>& x = bloque;
            const double bm = suma_seq(x) / B;
            std::vector<double> t(x.size());
            for (size_t i = 0; i < x.size(); ++i) t[i] = (x[i] - bm) * (x[i] - bm);
            const double bv = suma_seq(t) / B;
            const double lg = std::log(py_max(bv, shh));
            SG.add(lg);
            SG2.add(lg * lg);
            mG += 1;
            std::vector<double> xs = x;
            std::sort(xs.begin(), xs.end());
            const double med = mediana_ordenada(xs);
            std::vector<double> ab(x.size());
            for (size_t i = 0; i < x.size(); ++i) ab[i] = std::fabs(x[i] - med);
            const double bb = suma_seq(ab) / B;
            const double lb = std::log(py_max(bb, std::sqrt(shh / 2.0)));
            SL.add(lb);
            SL2.add(lb * lb);
            mL += 1;
            std::vector<double> abs_ord = ab;
            std::sort(abs_ord.begin(), abs_ord.end());
            const double mad = mediana_ordenada(abs_ord) / C_MAD;
            if (mad > 0.0) {
                SQ.add(std::log(std::sqrt(bv) / mad));
                mQ += 1;
            }
            bloque.clear();
            if (mG >= 3) {
                const double sg = SG.val();
                const double V = (SG2.val() - sg * sg / static_cast<double>(mG)) / static_cast<double>(mG - 1);
                const double E = V - psi1((B - 1) / 2.0);
                nu0 = E > Z1A * V * std::sqrt(2.0 / static_cast<double>(mG - 1)) ? 2.0 * psi1_inv(E) : INF;
                const double sl = SL.val();
                const double VL = (SL2.val() - sl * sl / static_cast<double>(mL)) / static_cast<double>(mL - 1);
                const double EL = VL - tablas_lazo::VLAP_MEDIANA;
                alfa0 = EL > Z1A * VL * std::sqrt(2.0 / static_cast<double>(mL - 1)) ? psi1_inv(EL) : INF;
            }
        }
        // g. autocorrelacion de |r|
        if (uN >= 1) P1.add(u * ultimos.back());
        if (uN >= 10) P10.add(u * ultimos[ultimos.size() - 10]);
        uN += 1;
        S1.add(u);
        S2.add(u * u);
        if (primeros.size() < 10) primeros.push_back(u);
        ultimos.push_back(u);
        if (ultimos.size() > 10) ultimos.pop_front();
        phi_ok = false;
        if (uN >= 11) {
            const double Nn = static_cast<double>(uN);
            const double s1 = S1.val();
            const double ub = s1 / Nn;
            const double c0 = S2.val() / Nn - ub * ub;
            double ck[2];
            const int ks[2] = {1, 10};
            const Neumaier* Ps[2] = {&P1, &P10};
            for (int q = 0; q < 2; ++q) {
                const int k_ = ks[q];
                const double A = s1 - suma_seq(primeros.begin(), primeros.begin() + std::min<ptrdiff_t>(k_, static_cast<ptrdiff_t>(primeros.size())));
                const double Bk = s1 - suma_seq(ultimos.end() - std::min<ptrdiff_t>(k_, static_cast<ptrdiff_t>(ultimos.size())), ultimos.end());
                ck[q] = (Ps[q]->val() - ub * (A + Bk) + (Nn - k_) * ub * ub) / Nn;
            }
            if (c0 > 0.0) {
                const double r1 = ck[0] / c0, r10 = ck[1] / c0;
                if (r1 > 0.0 && r10 > 0.0 && r10 < r1) {
                    const double ph = std::pow(r10 / r1, 1.0 / 9.0);
                    double acc = 1.0;
                    for (int j = 1; j < 10; ++j) {
                        const double rj = r1 * std::pow(ph, static_cast<double>(j - 1));
                        acc += 2.0 * rj * rj;
                    }
                    if (r10 > Z1A * std::sqrt(acc / Nn)) {
                        phi_ok = true;
                        phi = ph;
                        w_phi = 2.5 * (-LN2 / std::log(ph));
                    }
                }
            }
        }
        // h. regimen
        if (N <= static_cast<int64_t>(CAP)) {
            if (N >= 8) {
                const double med = mediana_ordenada(ord);
                std::vector<double> a(ord.size());
                for (size_t i = 0; i < ord.size(); ++i) a[i] = std::fabs(ord[i] - med);
                const double mad = np_median_sel(std::move(a)) / C_MAD;
                if (mad > 0.0) {
                    const double Q = std::sqrt(M2 / static_cast<double>(N)) / mad;
                    laplace = Q >= q95_nulo(N);
                }
            }
        } else if (mQ >= 2) {
            const double z = (SQ.val() / static_cast<double>(mQ) - tablas_lazo::QBLOQ_MEDIA) /
                             (tablas_lazo::QBLOQ_SD / std::sqrt(static_cast<double>(mQ)));
            laplace = z >= Z1A;
        }
        // i. W_max
        if (cnt[H] >= EV_MIN_COMPLETOS && (bloque_nuevo || cnt[H] != cnt_w_max)) {
            cnt_w_max = cnt[H];
            recalcular_w_max();
        }
        // j. W_innov
        const int64_t lo = std::min(std::max<int64_t>(2 * s_min, 8), w_max);
        if (phi_ok) w_innov = static_cast<int64_t>(std::floor(py_min(py_max(w_phi, static_cast<double>(lo)), static_cast<double>(w_max)) + 0.5));
        else w_innov = w_max;
        // k. ventana
        win.push_back(r);
        if (static_cast<int64_t>(win.size()) > W) win.pop_front();
        len = std::min(len + 1, w_innov);
        k += 1;
        val_prev = val;
    }

    void recalcular_w_max() {
        using namespace oraculo;
        if (std::isinf(nu0)) { w_max = W; return; }
        std::vector<double> ex(H), pos(H);
        for (int h = 1; h <= H; ++h) ex[h - 1] = sum_rat[h].val() / static_cast<double>(cnt[h]) - 1.0;
        for (int h = 0; h < H; ++h) pos[h] = py_max(ex[h], 0.0);
        const double tot = suma_seq(pos);
        double exmax = ex[0];
        for (double e : ex) if (e > exmax) exmax = e;
        if (tot <= 0.0 || exmax <= 0.0) { w_max = W; return; }
        Neumaier acu;
        int64_t tv = H;
        for (int h = 0; h < H; ++h) {
            acu.add(pos[h]);
            if (acu.val() >= 0.5 * tot) { tv = h + 1; break; }
        }
        const double d0 = sum_d0.val() / static_cast<double>(cnt_d0);
        const double umb = 1.0 + exmax / 2.0;
        int64_t mejor = 8;
        for (int64_t n = 8; n <= W; ++n) {
            int64_t lat = H + 1;
            for (int64_t h = 1; h <= H; ++h) {
                Neumaier acc;
                for (int64_t idx = h - n; idx < h; ++idx)
                    acc.add(idx < 0 ? 1.0 : (idx == 0 ? d0 : py_max(1.0 + ex[static_cast<size_t>(idx - 1)], 0.0)));
                const double vw = acc.val() / static_cast<double>(n);
                const double est = std::isinf(nu0) ? 1.0 : (nu0 + static_cast<double>(n) * vw) / (nu0 + static_cast<double>(n) - 1.0);
                if (est >= umb) { lat = h; break; }
            }
            if (lat < tv) mejor = n;
        }
        w_max = mejor;
        tau_v = tv;
        umbral = umb;
    }

    double paso(double val, bool mask = false) {
        const PropuestaLazo p = calcular(val, mask);
        confirmar(val, mask, p);
        return p.d_nu;
    }

private:
    // max(a, b) / min(a, b) de Python: el primero salvo que el segundo sea estrictamente mayor / menor.
    static double py_max(double a, double b) { return b > a ? b : a; }
    static double py_min(double a, double b) { return b < a ? b : a; }
};

}  // namespace ark::dataset::ventanas
