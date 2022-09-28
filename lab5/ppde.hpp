#ifndef PPDE_HPP
#define PPDE_HPP

#include "../lab1_2/tridiag.hpp"

/* Parabolic Partial Differential Equation Solver */
class ppde_t {
    const double EPS = 1e-9;

    int n, K;
    double l, T, h, tau, theta;
    double a, b, c;

    double h2;

    using vec = std::vector<double>;
    vec u0, ul;
    vec uk, uk1;

    using tridiag = tridiag_t<double>;

    void get_explicit(int i, double & uik, double & dudx, double & d2udx2) {
        uik = uk[i];
        dudx = (uk[i + 1] - uk[i - 1]) / (2 * h);
        d2udx2 = (uk[i - 1] - 2 * uik + uk[i + 1]) / h2;
    }

    void get_implicit(vec & trd_a, vec & trd_b, vec & trd_c) {
        trd_b[0] = 1;
        trd_c[0] = 0;
        for (int i = 1; i < n - 1; ++i) {
            trd_a[i] = (theta * tau * (b / (2 * h) - a / h2));
            trd_b[i] = (theta * tau * ((2 * a) / h2 - c)) + 1;
            trd_c[i] = (theta * tau * -(a / h2 + b / (2 * h)));
        }
        trd_a.back() = 0;
        trd_b.back() = 1;
    }

    static void print_vec(std::ostream & out, const vec & v) {
        size_t n = v.size();
        for (size_t i = 0; i < n; ++i) {
            if (i) {
                out << ' ';
            }
            out << v[i];
        }
        out << '\n';
    }

public:
    friend std::istream & operator >> (std::istream & in, ppde_t & item) {
        in >> item.n >> item.K;
        in >> item.l >> item.T >> item.h >> item.tau >> item.theta;
        if (item.theta < 0 or item.theta > 1) {
            throw std::invalid_argument("Theta is invalid!");
        }
        item.h2 = item.h * item.h;
        in >> item.a >> item.b >> item.c;
        item.uk.resize(item.n);
        item.u0.resize(item.K);
        item.ul.resize(item.K);
        for (int i = 0; i < item.n; ++i) {
            in >> item.uk[i];
        }
        for (int k = 0; k < item.K; ++k) {
            in >> item.u0[k];
        }
        for (int k = 0; k < item.K; ++k) {
            in >> item.ul[k];
        }
        return in;
    }

    void solve(std::ostream & out) {
        if (theta < EPS) {
            solve_explicit(out);
        } else if (theta < 1) {
            solve_combo(out);
        } else {
            solve_implicit(out);
        }
    }

    void solve_implicit(std::ostream & out) {
        print_vec(out, uk);
        vec trd_a(n), trd_b(n), trd_c(n), trd_d(n);
        for (int k = 0; k < K - 1; ++k) {
            get_implicit(trd_a, trd_b, trd_c);
            trd_d[0] = u0[k + 1];
            for (int i = 1; i < n - 1; ++i) {
                trd_d[i] = uk[i];
            }
            trd_d.back() = ul[k + 1];
            tridiag trd(trd_a, trd_b, trd_c);
            uk1 = trd.solve(trd_d);
            print_vec(out, uk1);
            swap(uk1, uk);
        }
    }

    void solve_combo(std::ostream & out) {
        print_vec(out, uk);
        vec trd_a(n), trd_b(n), trd_c(n), trd_d(n);
        for (int k = 0; k < K - 1; ++k) {
            get_implicit(trd_a, trd_b, trd_c);
            trd_d[0] = u0[k + 1];
            for (int i = 1; i < n - 1; ++i) {
                double uik, dudx, d2udx2;
                get_explicit(i, uik, dudx, d2udx2);
                double rhs = a * d2udx2 + b * dudx + c * uik;
                trd_d[i] = uik + (1 - theta) * tau * rhs;
            }
            trd_d.back() = ul[k + 1];
            tridiag trd(trd_a, trd_b, trd_c);
            uk1 = trd.solve(trd_d);
            print_vec(out, uk1);
            swap(uk1, uk);
        }
    }

    void solve_explicit(std::ostream & out) {
        print_vec(out, uk);
        for (int k = 0; k < K - 1; ++k) {
            uk1.assign(n, 0);
            uk1[0] = u0[k + 1];
            uk1.back() = ul[k + 1];
            for (int i = 1; i < n - 1; ++i) {
                double uik, dudx, d2udx2;
                get_explicit(i, uik, dudx, d2udx2);
                double rhs = a * d2udx2 + b * dudx + c * uik;
                uk1[i] = uik + tau * rhs;
            }
            print_vec(out, uk1);
            swap(uk1, uk);
        }
    }
};

#endif /* PPDE_HPP */
