#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

#include "../lab1_2/tridiag.hpp"

class cubic_spline_t {
    using vec = std::vector<double>;
    using tridiag = tridiag_t<double>;
    size_t n;
    vec x;
    vec y;
    vec a, b, c, d;

    void build_spline() {
        vec h(n + 1);
        h[0] = NAN;
        for (size_t i = 1; i <= n; ++i) {
            h[i] = x[i] - x[i - 1];
        }
        vec eq_a(n - 1);
        vec eq_b(n - 1);
        vec eq_c(n - 1);
        vec eq_d(n - 1);
        for (size_t i = 2; i <= n; ++i) {
            eq_a[i - 2] = h[i - 1];
            eq_b[i - 2] = 2.0 * (h[i - 1] + h[i]);
            eq_c[i - 2] = h[i];
            eq_d[i - 2] = 3.0 * ((y[i] - y[i - 1]) / h[i] -
                                 (y[i - 1] - y[i - 2]) / h[i - 1]);
        }
        eq_a[0] = 0.0;
        eq_c.back() = 0.0;
        // for (size_t i = 0; i < n - 1; ++i) {
        //     printf("%lf %lf %lf %lf\n", eq_a[i], eq_b[i], eq_c[i], eq_d[i]);
        // }
        tridiag system_of_eq(eq_a, eq_b, eq_c);
        vec c_solved = system_of_eq.solve(eq_d);
        for (size_t i = 2; i <= n; ++i) {
            c[i] = c_solved[i - 2];
        }
        for (size_t i = 1; i <= n; ++i) {
            a[i] = y[i - 1];
        }
        for (size_t i = 1; i < n; ++i) {
            b[i] =
                (y[i] - y[i - 1]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        }
        c[1] = 0.0;
        b[n] = (y[n] - y[n - 1]) / h[n] - (2.0 / 3.0) * h[n] * c[n];
        d[n] = -c[n] / (3.0 * h[n]);
    }

   public:
    cubic_spline_t(const vec& _x, const vec& _y) {
        if (_x.size() != _y.size()) {
            throw std::invalid_argument("Sizes does not match");
        }
        x = _x;
        y = _y;
        n = x.size() - 1;
        a.resize(n + 1);
        b.resize(n + 1);
        c.resize(n + 1);
        d.resize(n + 1);
        build_spline();
    }

    friend std::ostream& operator<<(std::ostream& out,
                                    const cubic_spline_t& spline) {
        for (size_t i = 1; i <= spline.n; ++i) {
            out << "i = " << i << ", a = " << spline.a[i]
                << ", b = " << spline.b[i] << ", c = " << spline.c[i]
                << ", d = " << spline.d[i] << '\n';
        }
        return out;
    }

    double operator()(double x0) {
        for (size_t i = 1; i <= n; ++i) {
            if (x[i - 1] <= x0 and x0 <= x[i]) {
                double x1 = x0 - x[i - 1];
                double x2 = x1 * x1;
                double x3 = x2 * x1;
                return a[i] + b[i] * x1 + c[i] * x2 + d[i] * x3;
            }
        }
        return NAN;
    }
};

#endif /* CUBIC_SPLINE_HPP */
