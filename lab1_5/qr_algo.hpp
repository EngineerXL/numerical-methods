#ifndef QR_ALGO_HPP
#define QR_ALGO_HPP

#include <cmath>
#include <complex>
#include "../matrix.hpp"

class qr_algo {
private:
    using matrix = matrix_t<double>;
    using vec = std::vector<double>;
    using complex = std::complex<double>;
    using pcc = std::pair<complex, complex>;
    using vec_complex = std::vector<complex>;

    static constexpr double INF = 1e18;
    static constexpr complex COMPLEX_INF = complex(INF, INF);

    size_t n;
    matrix a;
    double eps;
    vec_complex eigen;

    double vtv(const vec & v) {
        double res = 0;
        for (double elem : v) {
            res += elem * elem;
        }
        return res;
    }

    double norm(const vec & v) {
        return std::sqrt(vtv(v));
    }

    matrix vvt(const vec & b) {
        size_t n_b = b.size();
        matrix res(n_b);
        for (size_t i = 0; i < n_b; ++i) {
            for (size_t j = 0; j < n_b; ++j) {
                res[i][j] = b[i] * b[j];
            }
        }
        return res;
    }

    double sign(double x) {
        if (x < eps) {
            return -1.0;
        } else if (x > eps) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    matrix householder(const vec & b, int id) {
        vec v(b);
        v[id] += sign(b[id]) * norm(b);
        return matrix(n, true) - (2.0 / vtv(v)) * vvt(v);
    }

    pcc solve_sq(double a11, double a12, double a21, double a22) {
        double a = 1.0;
        double b = -(a11 + a22);
        double c = a11 * a22 - a12 * a21;
        double d_sq = b * b - 4.0 * a * c;
        if (d_sq > eps) {
            complex bad(NAN, NAN);
            return std::make_pair(bad, bad);
        }
        complex d(0.0, std::sqrt(-d_sq));
        complex x1 = (-b + d) / (2.0 * a);
        complex x2 = (-b - d) / (2.0 * a);
        return std::make_pair(x1, x2);
    }

    bool check_diag() {
        for (size_t i = 0; i < n; ++i) {
            double col_sum = 0;
            for (size_t j = i + 2; j < n; ++j) {
                col_sum += a[j][i] * a[j][i];
            }
            double norm = std::sqrt(col_sum);
            if (!(norm < eps)) {
                return false;
            }
        }
        return true;
    }

    void calc_eigen() {
        for (size_t i = 0; i < n; ++i) {
            if (i < n - 1 and !(abs(a[i + 1][i]) < eps)) {
                auto [l1, l2] = solve_sq(a[i][i], a[i][i + 1], a[i + 1][i], a[i + 1][i + 1]);
                if (std::isnan(l1.real())) {
                    ++i;
                    continue;
                }
                eigen[i] = l1;
                eigen[++i] = l2;
            } else {
                eigen[i] = a[i][i];
            }
        }
    }

    bool check_eps() {
        if (!check_diag()) {
            return false;
        }
        vec_complex prev_eigen(eigen);
        calc_eigen();
        for (size_t i = 0; i < n; ++i) {
            double delta = std::norm(eigen[i] - prev_eigen[i]);
            if (delta > eps) {
                return false;
            }
        }
        return true;
    }

    void build() {
        iter_count = 0;
        while (!check_eps()) {
            ++iter_count;
            matrix q(n, true);
            matrix r(a);
            for (size_t i = 0; i < n - 1; ++i) {
                vec b(n);
                for (size_t j = i; j < n; ++j) {
                    b[j] = r[j][i];
                }
                matrix h = householder(b, i);
                q = q * h;
                r = h * r;
            }
            a = r * q;
        }
    }

public:
    int iter_count;

    qr_algo(const matrix & _a, double _eps) :
            n(_a.size()), a(_a), eps(_eps), eigen(n, COMPLEX_INF) {
        build();
    };

    vec_complex get_eigen_values() {
        calc_eigen();
        return eigen;
    }

    ~qr_algo() = default;
};

#endif /* QR_ALGO_HPP */
