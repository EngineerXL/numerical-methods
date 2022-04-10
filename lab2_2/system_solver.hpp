#ifndef SYSTEM_SOLVER_HPP
#define SYSTEM_SOLVER_HPP

#include "../lab1_1/lu.hpp"

int iter_count = 0;

const double a = 1;

double phi1(double x1, double x2) {
    (void)x1;
    return a + std::cos(x2);
}

double phi1_s(double x1, double x2) {
    (void)x1;
    return -std::sin(x2);
}

double phi2(double x1, double x2) {
    (void)x2;
    return a + std::sin(x1);
}

double phi2_s(double x1, double x2) {
    (void)x2;
    return std::cos(x1);
}

double phi(double x1, double x2) {
    return phi1_s(x1, x2) * phi2_s(x1, x2);
}

using pdd = std::pair<double, double>;

pdd iter_solve(double l1, double r1, double l2, double r2, double eps) {
    iter_count = 0;
    double x_1_k = r1;
    double x_2_k = r2;
    double q = -1;
    q = std::max(q, std::abs(phi(l1, r1)));
    q = std::max(q, std::abs(phi(l1, r2)));
    q = std::max(q, std::abs(phi(l2, r1)));
    q = std::max(q, std::abs(phi(l2, r2)));
    double eps_coef = q / (1 - q);
    double dx = 1;
    do {
        double x_1_k1 = phi1(x_1_k, x_2_k);
        double x_2_k1 = phi2(x_1_k, x_2_k);
        dx = eps_coef * (std::abs(x_1_k1 - x_1_k) + std::abs(x_2_k1 - x_2_k));
        ++iter_count;
        x_1_k = x_1_k1;
        x_2_k = x_2_k1;
    } while (dx > eps);
    return std::make_pair(x_1_k, x_2_k);
}

using matrix = matrix_t<double>;
using lu = lu_t<double>;
using vec = std::vector<double>;

double f1(double x1, double x2) {
    return x1 - std::cos(x2) - a;
}

double f2(double x1, double x2) {
    return x2 - std::sin(x1) - a;
}

matrix j(double x1, double x2) {
    matrix res(2);
    res[0][0] = 1.0;
    res[0][1] = std::sin(x2);
    res[1][0] = -std::cos(x1);
    res[1][1] = 1.0;
    return res;
}

double norm(const vec & v) {
    double res = 0;
    for (double elem : v) {
        res = std::max(res, std::abs(elem));
    }
    return res;
}

pdd newton_solve(double x1_0, double x2_0, double eps) {
    iter_count = 0;
    vec x_k = {x1_0, x2_0};
    double dx = 1;
    do {
        double x1 = x_k[0];
        double x2 = x_k[1];
        lu jacobi(j(x1, x2));
        vec f_k = {f1(x1, x2), f2(x1, x2)};
        vec x_k1 = x_k - jacobi.inv_matrix() * f_k;
        dx = norm(x_k1 - x_k);
        ++iter_count;
        x_k = x_k1;
    } while (dx > eps);
    return std::make_pair(x_k[0], x_k[1]);
}

#endif /* SYSTEM_SOLVER_HPP */
