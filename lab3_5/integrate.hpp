#ifndef INTEGRATE_HPP
#define INTEGRATE_HPP

#include <cmath>

#include "../lab3_1/interpolator.hpp"

using func = double(double);

double integrate_rect(double l, double r, double h, func f) {
    double x1 = l;
    double x2 = l + h;
    double res = 0;
    while (x1 < r) {
        res += h * f((x1 + x2) * 0.5);
        x1 = x2;
        x2 += h;
    }
    return res;
}

double integrate_trap(double l, double r, double h, func f) {
    double x1 = l;
    double x2 = l + h;
    double res = 0;
    while (x1 < r) {
        res += h * (f(x1) + f(x2));
        x1 = x2;
        x2 += h;
    }
    return res * 0.5;
}

using vec = std::vector<double>;

double integrate_simp(double l, double r, double h, func f) {
    double x1 = l;
    double x2 = l + h;
    double res = 0;
    while (x1 < r) {
        vec x = {x1, (x1 + x2) * 0.5, x2};
        vec y = {f(x[0]), f(x[1]), f(x[2])};
        inter_lagrange lagr(x, y);
        res += lagr().integrate(x1, x2);
        x1 = x2;
        x2 += h;
    }
    return res;
}

inline double runge_romberg(double Fh, double Fkh, double k, double p) {
    return (Fh - Fkh) / (std::pow(k, p) - 1.0);
}

#endif /* INTEGRATE_HPP */
