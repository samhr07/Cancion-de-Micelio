// numpy_formato.hpp -- str(ndarray) de numpy 2 para un vector float64 1-D (array2string por defecto: precision 8, floatmode
// 'maxprec', suppress=False, linewidth 75, umbral de resumen 1000, edgeitems 3, separador ' '). Lo usa el mensaje
// "Invalid value '<array>' for dtype 'int64'" de pandas 3 que el testigo de motor_ventanas provoca al sumar el choque a una
// columna entera. Dragon4 de numpy = la cifra unica mas corta si cabe en la precision; si no, el redondeo exacto (printf de la
// UCRT). Jefe de Arquitectos, 2026-11-03.
#pragma once

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace ark::dataset::ventanas::npfmt {

// Digitos significativos y exponente decimal de la cifra unica mas corta (to_chars) de |x| > 0.
inline void mas_corta(double x, std::string& dig, int& e10) {
    char b[64];
    const auto r = std::to_chars(b, b + sizeof b, std::fabs(x), std::chars_format::scientific);
    std::string t(b, r.ptr);
    const size_t pe = t.find('e');
    dig.clear();
    for (size_t k = 0; k < pe; ++k) if (t[k] != '.') dig += t[k];
    e10 = std::stoi(t.substr(pe + 1));
    while (dig.size() > 1 && dig.back() == '0') dig.pop_back();
}

// Cientifica: mantisa con a lo sumo `prec` decimales (unica si cabe; si no, redondeada) -> (entero, fraccion, exponente).
inline void cientifica(double x, int prec, std::string& ent, std::string& frac, int& e10) {
    std::string dig;
    mas_corta(x, dig, e10);
    if (static_cast<int>(dig.size()) <= prec + 1) {
        ent = dig.substr(0, 1);
        frac = dig.substr(1);
    } else {
        char b[128];
        std::snprintf(b, sizeof b, "%.*e", prec, std::fabs(x));
        std::string t(b);
        const size_t pe = t.find('e');
        e10 = std::stoi(t.substr(pe + 1));
        ent = t.substr(0, 1);
        frac = pe > 2 ? t.substr(2, pe - 2) : "";
        while (!frac.empty() && frac.back() == '0') frac.pop_back();
    }
    if (x < 0 || (x == 0 && std::signbit(x))) ent = "-" + ent;
}

// Posicional: a lo sumo `prec` decimales (unica si cabe; si no, redondeada) -> (entero, fraccion) sin ceros finales.
inline void posicional(double x, int prec, std::string& ent, std::string& frac) {
    std::string dig;
    int e10 = 0;
    const bool neg = x < 0 || (x == 0 && std::signbit(x));
    if (x == 0) { ent = neg ? "-0" : "0"; frac = ""; return; }
    mas_corta(x, dig, e10);
    const int dec_unica = std::max(0, static_cast<int>(dig.size()) - 1 - e10);
    if (dec_unica <= prec) {
        if (e10 >= 0) {
            const size_t ne = static_cast<size_t>(e10) + 1;
            if (dig.size() <= ne) { ent = dig + std::string(ne - dig.size(), '0'); frac = ""; }
            else { ent = dig.substr(0, ne); frac = dig.substr(ne); }
        } else {
            ent = "0";
            frac = std::string(static_cast<size_t>(-e10 - 1), '0') + dig;
        }
    } else {
        char b[512];
        std::snprintf(b, sizeof b, "%.*f", prec, std::fabs(x));
        std::string t(b);
        const size_t pp = t.find('.');
        ent = t.substr(0, pp);
        frac = pp == std::string::npos ? "" : t.substr(pp + 1);
        while (!frac.empty() && frac.back() == '0') frac.pop_back();
    }
    if (neg) ent = "-" + ent;
}

inline std::string str_array_f64(const std::vector<double>& a) {
    const size_t n = a.size();
    if (n == 0) return "[]";
    const bool resumen = n > 1000;
    std::vector<double> datos;   // los que fijan el formato (los de los bordes si hay resumen)
    if (resumen) {
        datos.assign(a.begin(), a.begin() + 3);
        datos.insert(datos.end(), a.end() - 3, a.end());
    } else {
        datos = a;
    }
    std::vector<double> fin;
    for (double v : datos) if (std::isfinite(v)) fin.push_back(v);
    bool expo = false;
    {
        double mx = -1.0, mn = -1.0;
        for (double v : fin) {
            if (v == 0) continue;
            const double av = std::fabs(v);
            if (mx < 0 || av > mx) mx = av;
            if (mn < 0 || av < mn) mn = av;
        }
        if (mx >= 0 && (mx >= 1e8 || mn < 0.0001 || mx / mn > 1000.0)) expo = true;
    }
    int prec = 8, pad_left = 0, pad_right = 0, exp_size = -1;
    bool min_dig_eq = false;
    if (fin.empty()) {
        pad_left = 0; pad_right = 0; exp_size = -1;
    } else if (expo) {
        int maxf = 0, maxe = 0, maxi = 0;
        for (double v : fin) {
            std::string e, f; int e10;
            cientifica(v, 8, e, f, e10);
            const int ed = std::max(2, static_cast<int>(std::to_string(std::abs(e10)).size()));
            maxe = std::max(maxe, ed + 1);
            maxf = std::max(maxf, static_cast<int>(f.size()));
            maxi = std::max(maxi, static_cast<int>(e.size()));
        }
        exp_size = maxe - 1;
        prec = maxf;
        min_dig_eq = true;
        pad_left = maxi;
        pad_right = exp_size + 2 + prec;
    } else {
        int maxf = 0, maxi = 0;
        for (double v : fin) {
            std::string e, f;
            posicional(v, 8, e, f);
            maxi = std::max(maxi, static_cast<int>(e.size()));
            maxf = std::max(maxf, static_cast<int>(f.size()));
        }
        pad_left = maxi;
        pad_right = maxf;
    }
    if (fin.size() != datos.size()) {
        bool neginf = false;
        for (double v : datos) if (std::isinf(v) && v < 0) neginf = true;
        const int offset = pad_right + 1;
        pad_left = std::max({pad_left, 3 - offset, 3 + (neginf ? 1 : 0) - offset});
    }
    auto formatear = [&](double v) -> std::string {
        if (!std::isfinite(v)) {
            std::string ret = std::isnan(v) ? "nan" : (v < 0 ? "-inf" : "inf");
            const int w = pad_left + pad_right + 1 - static_cast<int>(ret.size());
            return std::string(static_cast<size_t>(std::max(0, w)), ' ') + ret;
        }
        std::string e, f;
        if (expo) {
            int e10;
            cientifica(v, prec, e, f, e10);
            if (min_dig_eq && static_cast<int>(f.size()) < prec) f += std::string(static_cast<size_t>(prec - static_cast<int>(f.size())), '0');
            std::string ex = std::to_string(std::abs(e10));
            if (static_cast<int>(ex.size()) < exp_size) ex = std::string(static_cast<size_t>(exp_size) - ex.size(), '0') + ex;
            std::string s = e + "." + f + "e" + (e10 < 0 ? "-" : "+") + ex;
            if (static_cast<int>(e.size()) < pad_left) s = std::string(static_cast<size_t>(pad_left) - e.size(), ' ') + s;
            return s;
        }
        posicional(v, prec, e, f);
        std::string s = e + "." + f;
        if (static_cast<int>(f.size()) < pad_right) s += std::string(static_cast<size_t>(pad_right) - f.size(), ' ');
        if (static_cast<int>(e.size()) < pad_left) s = std::string(static_cast<size_t>(pad_left) - e.size(), ' ') + s;
        return s;
    };
    // _formatArray (1-D): sangria " ", ancho 75, elementos hasta 74 columnas
    const std::string sangria = " ";
    const size_t ancho = 74;
    std::string s, linea = sangria;
    auto extender = [&](const std::string& palabra) {
        bool envolver = linea.size() + palabra.size() > ancho;
        if (linea.size() <= sangria.size()) envolver = false;
        if (envolver) {
            size_t fin_l = linea.size();
            while (fin_l > 0 && std::isspace(static_cast<unsigned char>(linea[fin_l - 1]))) --fin_l;
            s += linea.substr(0, fin_l) + "\n";
            linea = sangria;
        }
        linea += palabra;
    };
    if (resumen) {
        for (size_t i = 0; i < 3; ++i) { extender(formatear(a[i])); linea += " "; }
        extender("...");
        linea += " ";
        for (size_t i = 3; i > 1; --i) { extender(formatear(a[n - i])); linea += " "; }
    } else {
        for (size_t i = 0; i + 1 < n; ++i) { extender(formatear(a[i])); linea += " "; }
    }
    extender(formatear(a[n - 1]));
    s += linea;
    return "[" + s.substr(sangria.size()) + "]";
}

}  // namespace ark::dataset::ventanas::npfmt
