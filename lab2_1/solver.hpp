#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cmath>

int iter_count = 0;

double f(double x) {
    return std::sin(x) - 2 * x * x + 0.5;
}

double f_s(double x) {
    return std::cos(x) - 4 * x;
}

double phi(double x) {
    return std::sqrt(0.5 * std::sin(x) + 0.25);
}

double phi_s(double x) {
    return std::cos(x) / (4.0 * phi(x));
}

double iter_solve(double l, double r, double eps) {
    iter_count = 0;
    double x_k = r;
    double dx = 1;
    double q = std::max(std::abs(phi_s(l)), std::abs(phi_s(r)));
    double eps_coef = q / (1 - q);
    do {
        double x_k1 = phi(x_k);
        dx = eps_coef * std::abs(x_k1 - x_k);
        ++iter_count;
        x_k = x_k1;
    } while (dx > eps);
    return x_k;
}

double newton_solve(double x0, double eps) {
    iter_count = 0;
    double x_k = x0;
    double dx = 1;
    do {
        double x_k1 = x_k - f(x_k) / f_s(x_k);
        dx = std::abs(x_k1 - x_k);
        ++iter_count;
        x_k = x_k1;
    } while (dx > eps);
    return x_k;
}

#endif /* SOLVER_HPP */
