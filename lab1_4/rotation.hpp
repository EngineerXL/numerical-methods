#ifndef ROTATION_HPP
#define ROTATION_HPP

#include <cmath>
#include "../matrix.hpp"

class rotation {
private:
    using matrix = matrix_t<double>;
    using vec = std::vector<double>;

    static constexpr double GLOBAL_EPS = 1e-9;

    size_t n;
    matrix a;
    double eps;
    matrix v;

    static double norm(const matrix & m) {
        double res = 0;
        for (size_t i = 0; i < m.rows(); ++i) {
            for (size_t j = 0; j < m.cols(); ++j) {
                if (i == j) {
                    continue;
                }
                res += m[i][j] * m[i][j];
            }
        }
        return std::sqrt(res);
    }

    double calc_phi(size_t i, size_t j) {
        if (std::abs(a[i][i] - a[j][j]) < GLOBAL_EPS) {
            return std::atan2(1.0, 1.0);
        } else {
            return 0.5 * std::atan2(2 * a[i][j], a[i][i] - a[j][j]);
        }
    }

    matrix create_rotation(size_t i, size_t j, double phi) {
        matrix u = matrix::identity(n);
        u[i][i] = std::cos(phi);
        u[i][j] = -std::sin(phi);
        u[j][i] = std::sin(phi);
        u[j][j] = std::cos(phi);
        return u;
    }

    void build() {
        iter_count = 0;
        while (norm(a) > eps) {
            ++iter_count;
            size_t i = 0, j = 1;
            for (size_t ii = 0; ii < n; ++ii) {
                for (size_t jj = 0; jj < n; ++jj) {
                    if (ii == jj) {
                        continue;
                    }
                    if (std::abs(a[ii][jj]) > std::abs(a[i][j])) {
                        i = ii;
                        j = jj;
                    }
                }
            }
            double phi = calc_phi(i, j);
            matrix u = create_rotation(i, j, phi);
            v = v * u;
            a = u.t() * a * u;
        }
    }

public:
    int iter_count;

    rotation(const matrix & _a, double _eps) {
        if (_a.rows() != _a.cols()) {
            throw std::invalid_argument("Matrix is not square");
        }
        a = matrix(_a);
        n = a.rows();
        eps = _eps;
        v = matrix::identity(n);
        build();
    };

    matrix get_eigen_vectors() {
        return v;
    }

    vec get_eigen_values() {
        vec res(n);
        for (size_t i = 0; i < n; ++i) {
            res[i] = a[i][i];
        }
        return res;
    }

    ~rotation() = default;
};

#endif /* ROTATION_HPP */
