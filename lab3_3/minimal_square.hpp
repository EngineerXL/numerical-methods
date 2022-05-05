#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include "../lab1_1/lu.hpp"
#include "../polynom.hpp"

class minimal_square_t {
    using vec = std::vector<double>;
    using matrix = matrix_t<double>;
    using lu = lu_t<double>;
    size_t n;
    vec x;
    vec y;
    polynom f;

    void build(size_t degree) {
        matrix lhs(degree);
        vec rhs(degree);
        for (size_t k = 0; k < degree; ++k) {
            for (size_t i = 0; i < degree; ++i) {
                double sum_x = 0;
                for (size_t j = 0; j < n; ++j) {
                    sum_x += bpow(x[j], k + i);
                }
                lhs[k][i] = sum_x;
            }
        }
        for (size_t k = 0; k < degree; ++k) {
            double sum_yx = 0;
            for (size_t j = 0; j < n; ++j) {
                sum_yx += y[j] * bpow(x[j], k);
            }
            rhs[k] = sum_yx;
        }
        lu ms_lu(lhs);
        vec solved = ms_lu.solve(rhs);
        f = polynom(solved);
    }

public:
    minimal_square_t(const vec & _x, const vec & _y, size_t degree) {
        if (_x.size() != _y.size()) {
            throw std::invalid_argument("Sizes does not match");
        }
        x = _x;
        y = _y;
        n = x.size();
        ++degree;
        build(degree);
    }

    friend std::ostream & operator << (std::ostream & out, const minimal_square_t & item) {
        out << item.f;
        return out;
    }

    double mse() {
        double res = 0;
        for (size_t i = 0; i < n; ++i) {
            res += bpow(f(x[i]) - y[i], 2);
        }
        return res;
    }

    double operator () (double x0) {
        return f(x0);
    }

    static double bpow(double x, int64_t y) {
        double z = 1.0;
        while (y) {
            if (y & 1) {
                z = z * x;
            }
            x = x * x;
            y = y >> 1;
        }
        return z;
    }
};

#endif /* CUBIC_SPLINE_HPP */
