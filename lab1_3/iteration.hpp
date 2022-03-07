#ifndef ITERATION_HPP
#define ITERATION_HPP

#include <cmath>
#include "../matrix.hpp"

class iter_solver {
private:
    using matrix = matrix_t<double>;
    using vec = std::vector<double>;

    matrix a;
    size_t n;
    double eps;

    static constexpr double INF = 1e18;
public:
    int iter_count;

    iter_solver(const matrix & _a, double _eps = 1e-6) : a(_a), n(_a.size()), eps(_eps) {}

    static double norm(matrix & m) {
        double res = -INF;
        for (size_t i = 0; i < m.size(); ++i) {
            res = std::max(res, norm(m[i]));
        }
        return res;
    }

    static double norm(const vec & v) {
        double res = 0;
        for (double elem : v) {
            res += std::abs(elem);
        }
        return res;
    }

    std::pair<matrix, vec> precalc_ab(const vec & b, matrix & alpha, vec & beta) {
        for (size_t i = 0; i < n; ++i) {
            beta[i] = b[i] / a[i][i];
            for (size_t j = 0; j < n; ++j) {
                if (i != j) {
                    alpha[i][j] = -a[i][j] / a[i][i];
                }
            }
        }
        return std::make_pair(alpha, beta);
    }

    vec solve_simple(const vec & b) {
        matrix alpha(n);
        vec beta(n);
        precalc_ab(b, alpha, beta);
        double eps_coef = 1.0;
        if (norm(alpha) - 1.0 < eps) {
            eps_coef = norm(alpha) / (1.0 - norm(alpha));
        }
        double eps_k = 1.0;
        vec x(beta);
        iter_count = 0;
        while (eps_k > eps) {
            vec x_k = beta + alpha * x;
            eps_k = eps_coef * norm(x_k - x);
            x = x_k;
            ++iter_count;
        }
        return x;
    }

    vec zeidel(const vec & x, const matrix & alpha, const vec & beta) {
        vec x_k(beta);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                x_k[i] += x_k[j] * alpha[i][j];
            }
            for (size_t j = i; j < n; ++j) {
                x_k[i] += x[j] * alpha[i][j];
            }
        }
        return x_k;
    }

    vec solve_zeidel(const vec & b) {
        matrix alpha(n);
        vec beta(n);
        precalc_ab(b, alpha, beta);
        matrix c(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i; j < n; ++j) {
                c[i][j] = alpha[i][j];
            }
        }
        double eps_coef = 1.0;
        if (norm(alpha) - 1.0 < eps) {
            eps_coef = norm(c) / (1.0 - norm(alpha));
        }
        double eps_k = 1.0;
        vec x(beta);
        iter_count = 0;
        while (eps_k > eps) {
            vec x_k = zeidel(x, alpha, beta);
            eps_k = eps_coef * norm(x_k - x);
            x = x_k;
            ++iter_count;
        }
        return x;
    }

    ~iter_solver() = default;
};

#endif /* ITERATION_HPP */
