#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include "../lab1_1/lu.hpp"
#include "../polynom.hpp"

class minimal_square_t {
    using vec = std::vector<double>;
    using matrix = matrix_t<double>;
    using lu = lu_t<double>;

    using func = std::function<double(double)>;
    using vf = std::vector<func>;

    size_t n;
    vec x;
    vec y;
    size_t m;
    vec a;
    vf phi;

    void build() {
        matrix lhs(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                lhs[i][j] = phi[j](x[i]);
            }
        }
        matrix lhs_t = lhs.t();
        lu lhs_lu(lhs_t * lhs);
        vec rhs = lhs_t * y;
        a = lhs_lu.solve(rhs);
    }

    double get(double x0) {
        double res = 0.0;
        for (size_t i = 0; i < m; ++i) {
            res += a[i] * phi[i](x0);
        }
        return res;
    }

public:
    minimal_square_t(const vec & _x, const vec & _y, const vf & _phi) {
        if (_x.size() != _y.size()) {
            throw std::invalid_argument("Sizes does not match");
        }
        n = _x.size();
        x = _x;
        y = _y;
        m = _phi.size();
        a.resize(m);
        phi = _phi;
        build();
    }

    friend std::ostream & operator << (std::ostream & out, const minimal_square_t & item) {
        for (size_t i = 0; i < item.m; ++i) {
            if (i) {
                out << ' ';
            }
            out << item.a[i];
        }
        return out;
    }

    double mse() {
        double res = 0;
        for (size_t i = 0; i < n; ++i) {
            res += std::pow(get(x[i]) - y[i], 2.0);
        }
        return res;
    }

    double operator () (double x0) {
        return get(x0);
    }
};

#endif /* CUBIC_SPLINE_HPP */
