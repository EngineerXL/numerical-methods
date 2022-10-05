#ifndef PPDE_HPP
#define PPDE_HPP

#include "../lab1_2/tridiag.hpp"

/* Parabolic Partial Differential Equation Solver */
class ppde_t {
    const double EPS = 1e-9;

    int n, K, boundary;
    double l, T, h, tau, theta;
    double a, b, c;

    double alpha_0, beta_0, alpha_l, beta_l;

    double h2;

    using vec = std::vector<double>;
    vec gamma_0, gamma_l;
    vec uk, uk1;

    using tridiag = tridiag_t<double>;
    vec trd_a, trd_b, trd_c, trd_d;

    void gen_explicit(int i, double & uik, double & dudx, double & d2udx2) {
        uik = uk[i];
        dudx = (uk[i + 1] - uk[i - 1]) / (2 * h);
        d2udx2 = (uk[i - 1] - 2 * uik + uk[i + 1]) / h2;
    }

    vec c_3t2p_0, c_3t2p_l;

    void prepare_3t2p() {
        c_3t2p_0.resize(4, 0);
        c_3t2p_0[0] = 2 * h * beta_0 - 3 * alpha_0;
        c_3t2p_0[1] = 4 * alpha_0;
        c_3t2p_0[2] = -alpha_0;
        c_3t2p_0[3] = 2 * h;

        c_3t2p_l.resize(4, 0);
        c_3t2p_l[0] = alpha_l;
        c_3t2p_l[1] = -4 * alpha_l;
        c_3t2p_l[2] = 2 * h * beta_l + 3 * alpha_l;
        c_3t2p_l[3] = 2 * h;
    }

    vec c_2t2p_0, c_2t2p_l;

    void prepare_2t2p() {
        c_2t2p_0.resize(4, 0);
        c_2t2p_0[0] = alpha_0 * (1 / tau + 2 * a / h2 - c) + beta_0 * (b - 2 * a / h);
        c_2t2p_0[1] = -2 * alpha_0 * a / h2;
        c_2t2p_0[2] = b - 2 * a / h;
        c_2t2p_0[3] = alpha_0 / tau;

        c_2t2p_l.resize(4, 0);
        c_2t2p_l[0] = -2 * alpha_l * a / h2;
        c_2t2p_l[1] = alpha_l * (1 / tau + 2 * a / h2 - c) + beta_l * (b + 2 * a / h);
        c_2t2p_l[2] = b + 2 * a / h;
        c_2t2p_l[3] = alpha_l / tau;
    }

    void prepare_trd() {
        trd_a.resize(n, 0);
        trd_b.resize(n, 0);
        trd_c.resize(n, 0);
        trd_d.resize(n, 0);
    }

    void gen_boundary_0_2t1p(int k) {
        trd_b[0] = beta_0 - alpha_0 / h;
        trd_c[0] = alpha_0 / h;
        trd_d[0] = gamma_0[k + 1];
    }

    void gen_boundary_0_3t2p(int k) {
        double coef_row = c_3t2p_0[2] / trd_c[1];
        trd_b[0] = c_3t2p_0[0] - trd_a[1] * coef_row;
        trd_c[0] = c_3t2p_0[1] - trd_b[1] * coef_row;
        trd_d[0] = c_3t2p_0[3] * gamma_0[k + 1] - trd_d[1] * coef_row;
    }

    void gen_boundary_0_2t2p(int k) {
        trd_b[0] = c_2t2p_0[0];
        trd_c[0] = c_2t2p_0[1];
        trd_d[0] = c_2t2p_0[2] * gamma_0[k + 1] + c_2t2p_0[3] * uk[0];
    }

    void gen_boundary_0(int k) {
        if (boundary == 1) {
            gen_boundary_0_2t1p(k);
        } else if (boundary == 2) {
            gen_boundary_0_3t2p(k);
        } else {
            gen_boundary_0_2t2p(k);
        }
    }

    void gen_boundary_l_2t1p(int k) {
        trd_a.back() = -alpha_l / h;
        trd_b.back() = alpha_l / h + beta_l;
        trd_d.back() = gamma_l[k + 1];
    }

    void gen_boundary_l_3t2p(int k) {
        double coef_row = c_3t2p_l[0] / trd_a[n - 2];
        trd_a.back() = c_3t2p_l[1] - trd_b[n - 2] * coef_row;
        trd_b.back() = c_3t2p_l[2] - trd_c[n - 2] * coef_row;
        trd_d.back() = c_3t2p_l[3] * gamma_l[k + 1] - trd_d[n - 2] * coef_row;
    }

    void gen_boundary_l_2t2p(int k) {
        trd_a.back() = c_2t2p_l[0];
        trd_b.back() = c_2t2p_l[1];
        trd_d.back() = c_2t2p_l[2] * gamma_l[k + 1] + c_2t2p_l[3] * uk[n - 1];
    }

    void gen_boundary_l(int k) {
        if (boundary == 1) {
            gen_boundary_l_2t1p(k);
        } else if (boundary == 2) {
            gen_boundary_l_3t2p(k);
        } else {
            gen_boundary_l_2t2p(k);
        }
    }

    void gen_implicit(int k, bool mode_combo) {
        prepare_trd();
        for (int i = 1; i < n - 1; ++i) {
            trd_a[i] = (theta * tau * (b / (2 * h) - a / h2));
            trd_b[i] = (theta * tau * ((2 * a) / h2 - c)) + 1;
            trd_c[i] = (theta * tau * -(a / h2 + b / (2 * h)));
            if (mode_combo) {
                double uik, dudx, d2udx2;
                gen_explicit(i, uik, dudx, d2udx2);
                double rhs = a * d2udx2 + b * dudx + c * uik;
                trd_d[i] = uik + (1 - theta) * tau * rhs;
            } else {
                trd_d[i] = uk[i];
            }
        }
        gen_boundary_0(k);
        gen_boundary_l(k);
    }

    void boundary_explicit_2t1p(int k) {
        uk1[0] = (gamma_0[k + 1] - alpha_0 * uk1[1] / h) / (beta_0 - alpha_0 / h);
        uk1[n - 1] = (gamma_l[k + 1] + alpha_l * uk1[n - 2] / h) / (alpha_l / h + beta_l);
    }

    void boundary_explicit_3t2p(int k) {
        double rhs0 = c_3t2p_0[3] * gamma_0[k + 1] - c_3t2p_0[2] * uk1[2] - c_3t2p_0[1] * uk1[1];
        double lhs0 = c_3t2p_0[0];
        uk1[0] = rhs0 / lhs0;
        double rhsl = c_3t2p_l[3] * gamma_l[k + 1] - c_3t2p_l[0] * uk1[n - 3] - c_3t2p_l[1] * uk1[n - 2];
        double lhsl = c_3t2p_l[2];
        uk1[n - 1] = rhsl / lhsl;
    }

    void boundary_explicit_2t2p(int k) {
        double rhs0 = c_2t2p_0[2] * gamma_0[k + 1] + c_2t2p_0[3] * uk[0] - c_2t2p_0[1] * uk1[1];
        double lhs0 = c_2t2p_0[0];
        uk1[0] = rhs0 / lhs0;
        double rhsl = c_2t2p_l[2] * gamma_l[k + 1] + c_2t2p_l[3] * uk[n - 1] - c_2t2p_l[0] * uk1[n - 2];
        double lhsl = c_2t2p_l[1];
        uk1[n - 1] = rhsl / lhsl;
    }

    void boundary_explicit(int k) {
        if (boundary == 1) {
            boundary_explicit_2t1p(k);
        } else if (boundary == 2) {
            boundary_explicit_3t2p(k);
        } else {
            boundary_explicit_2t2p(k);
        }
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
        in >> item.n >> item.K >> item.boundary;
        in >> item.l >> item.T >> item.h >> item.tau >> item.theta;
        if (item.theta < 0 or item.theta > 1) {
            throw std::invalid_argument("Theta is invalid!");
        }
        item.h2 = item.h * item.h;
        in >> item.a >> item.b >> item.c;
        in >> item.alpha_0 >> item.beta_0 >> item.alpha_l >> item.beta_l;
        item.uk.resize(item.n);
        item.gamma_0.resize(item.K);
        item.gamma_l.resize(item.K);
        for (int i = 0; i < item.n; ++i) {
            in >> item.uk[i];
        }
        for (int k = 0; k < item.K; ++k) {
            in >> item.gamma_0[k];
        }
        for (int k = 0; k < item.K; ++k) {
            in >> item.gamma_l[k];
        }
        return in;
    }

    void solve(std::ostream & out) {
        prepare_3t2p();
        prepare_2t2p();
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
        for (int k = 0; k < K - 1; ++k) {
            gen_implicit(k, false);
            tridiag trd(trd_a, trd_b, trd_c);
            uk1 = trd.solve(trd_d);
            print_vec(out, uk1);
            swap(uk1, uk);
        }
    }

    void solve_combo(std::ostream & out) {
        print_vec(out, uk);
        for (int k = 0; k < K - 1; ++k) {
            gen_implicit(k, true);
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
            for (int i = 1; i < n - 1; ++i) {
                double uik, dudx, d2udx2;
                gen_explicit(i, uik, dudx, d2udx2);
                double rhs = a * d2udx2 + b * dudx + c * uik;
                uk1[i] = uik + tau * rhs;
            }
            boundary_explicit(k);
            print_vec(out, uk1);
            swap(uk1, uk);
        }
    }
};

#endif /* PPDE_HPP */
