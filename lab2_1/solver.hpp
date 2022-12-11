#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <cmath>

int iter_count = 0;

double f(double x) { return std::sin(x) - 2.0 * x * x + 0.5; }

double f_s(double x) { return std::cos(x) - 4.0 * x; }

double f_ss(double x) { return -std::sin(x) - 4.0; }

double phi(double x) { return std::sqrt(0.5 * std::sin(x) + 0.25); }

double phi_s(double x) { return std::cos(x) / (4.0 * phi(x)); }

double iter_solve(double l, double r, double eps) {
    iter_count = 0;
    double x_k = r;
    double dx = 1.0;
    double q = std::max(std::abs(phi_s(l)), std::abs(phi_s(r)));
    double eps_coef = q / (1.0 - q);
    do {
        double x_k1 = phi(x_k);
        dx = eps_coef * std::abs(x_k1 - x_k);
        ++iter_count;
        x_k = x_k1;
    } while (dx > eps);
    return x_k;
}

double newton_solve(double l, double r, double eps) {
    double x0 = l;
    if (!(f(x0) * f_ss(x0) > eps)) {
        x0 = r;
    }
    iter_count = 0;
    double x_k = x0;
    double dx = 1.0;
    do {
        double x_k1 = x_k - f(x_k) / f_s(x_k);
        dx = std::abs(x_k1 - x_k);
        ++iter_count;
        x_k = x_k1;
    } while (dx > eps);
    return x_k;
}

#endif /* SOLVER_HPP */
