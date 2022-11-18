#ifndef HPDE_HPP
#define HPDE_HPP

#include "../lab1_2/tridiag.hpp"

/* Hyperbolic Partial Differential Equation Solver */
class hpde_t {
    const double EPS = 1e-9;

    using vec = std::vector<double>;
    using vecvec = std::vector<vec>;

    int n, K, u1_degree, boundary;
    double l, T, h, tau;
    int theta;

    double a, b, c, d;
    vecvec f;

    /*
     * u_prev = u_{k-1}
     * u_k = u_k
     * u_next = u_{k + 1}
     */
    vec u_prev, u_k, u_next;

    double alpha_0, beta_0, alpha_l, beta_l;
    vec gamma_0, gamma_l;

    vec ddx_psi1, d2dx2_psi1;

    double h2, tau2;

    void prepare_u1_1() {
        for (int i = 0; i < n; ++i) {
            u_k[i] = u_prev[i] + u_k[i] * tau;
        }
    }

    void prepare_u1_2() {
        for (int i = 0; i < n; ++i) {
            u_k[i] = u_prev[i] + tau2 / 2 * (a * d2dx2_psi1[i] + b * ddx_psi1[i] + c * u_prev[i] + f[i][0]) + u_k[i] * (tau + d * tau2 / 2);
        }
    }

    void prepare_u1() {
        if (u1_degree == 1) {
            prepare_u1_1();
        } else {
            prepare_u1_2();
        }
    }

    vec c_expl;
    const int EXPLICIT_COEFS = 6;

    void precalc_explicit() {
        c_expl.resize(EXPLICIT_COEFS);
        c_expl[0] = h2 * (2 - d * tau);
        c_expl[1] = tau2 * (2 * a - h * b);
        c_expl[2] = 4 * h2 + tau2 * (-4 * a + 2 * h2 * c);
        c_expl[3] = tau2 * (2 * a + b * h);
        c_expl[4] = h2 * (-2 - d * tau);
        c_expl[5] = 2 * h2 * tau2;
    }

    vec c_impl;
    const int IMPLICIT_COEFS = 6;

    void precalc_implicit() {
        c_impl.resize(IMPLICIT_COEFS);
        c_impl[0] = tau2 * (-2 * a + b * h);
        c_impl[1] = 2 * h2 + tau2 * (4 * a - 2 * h2 * c) - d * h2 * tau;
        c_impl[2] = tau2 * (-2 * a - b * h);
        c_impl[3] = 4 * h2;
        c_impl[4] = h2 * (-2 - d * tau);
        c_impl[5] = 2 * h2 * tau2;
    }

    vec c_3t2p_0, c_3t2p_l;
    const int COEFS_3T2P = 4;

    void prepare_3t2p() {
        c_3t2p_0.resize(COEFS_3T2P);
        c_3t2p_0[0] = 2 * h * beta_0 - 3 * alpha_0;
        c_3t2p_0[1] = 4 * alpha_0;
        c_3t2p_0[2] = -alpha_0;
        c_3t2p_0[3] = 2 * h;

        c_3t2p_l.resize(COEFS_3T2P);
        c_3t2p_l[0] = alpha_l;
        c_3t2p_l[1] = -4 * alpha_l;
        c_3t2p_l[2] = 2 * h * beta_l + 3 * alpha_l;
        c_3t2p_l[3] = 2 * h;
    }

    vec c_2t2p_0, c_2t2p_l;
    const int COEFS_2T2P = 6;

    void prepare_2t2p() {
        c_2t2p_0.resize(COEFS_2T2P);
        c_2t2p_0[0] = -alpha_0 * a - alpha_0 * h2 / 2 * (1 / tau2 - c - d / (2 * tau)) + beta_0 * (a * h - b * h2 / 2);
        c_2t2p_0[1] = alpha_0 * a;
        c_2t2p_0[2] = alpha_0 * h2 / tau2;
        c_2t2p_0[3] = -alpha_0 * h2 / (2 * tau2) + alpha_0 * h2 * d / (4 * tau);
        c_2t2p_0[4] = a * h - b * h2 / 2;
        c_2t2p_0[5] = alpha_0 * h2 / 2;

        c_2t2p_l.resize(COEFS_2T2P);
        c_2t2p_l[0] = -alpha_l * a;
        c_2t2p_l[1] = alpha_l * a + alpha_l * h2 / 2 * (1 / tau2 - c - d / (2 * tau)) + beta_l * (a * h + b * h2 / 2);
        c_2t2p_l[2] = -alpha_l * h2 / tau2;
        c_2t2p_l[3] = alpha_l * h2 / (2 * tau2) + alpha_l * h2 * d / (4 * tau);
        c_2t2p_l[4] = a * h + b * h2 / 2;
        c_2t2p_l[5] = -alpha_l * h2 / 2;
    }

    using tridiag = tridiag_t<double>;
    vec trd_a, trd_b, trd_c, trd_d;

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
        trd_d[0] = c_2t2p_0[4] * gamma_0[k + 1] - c_2t2p_0[5] * f[0][k + 1] - c_2t2p_0[2] * u_k[0] - c_2t2p_0[3] * u_prev[0];
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
        trd_d.back() = c_2t2p_0[4] * gamma_l[k + 1] - c_2t2p_0[5] * f[n - 1][k + 1] - c_2t2p_0[2] * u_k[n - 1] + c_2t2p_0[3] * u_prev[n - 1];
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

    void gen_implicit(int k) {
        prepare_trd();
        for (int i = 1; i < n - 1; ++i) {
            trd_a[i] = c_impl[0];
            trd_b[i] = c_impl[1];
            trd_c[i] = c_impl[2];
            trd_d[i] = c_impl[3] * u_k[i] + c_impl[4] * u_prev[i] + c_impl[5] * f[i][k + 1];
        }
        gen_boundary_0(k);
        gen_boundary_l(k);
    }

    void boundary_explicit_2t1p(int k) {
        u_next[0] = (gamma_0[k + 1] - alpha_0 * u_next[1] / h) / (beta_0 - alpha_0 / h);
        u_next[n - 1] = (gamma_l[k + 1] + alpha_l * u_next[n - 2] / h) / (alpha_l / h + beta_l);
    }

    void boundary_explicit_3t2p(int k) {
        double rhs0 = c_3t2p_0[3] * gamma_0[k + 1] - c_3t2p_0[2] * u_next[2] - c_3t2p_0[1] * u_next[1];
        double lhs0 = c_3t2p_0[0];
        u_next[0] = rhs0 / lhs0;
        double rhsl = c_3t2p_l[3] * gamma_l[k + 1] - c_3t2p_l[0] * u_next[n - 3] - c_3t2p_l[1] * u_next[n - 2];
        double lhsl = c_3t2p_l[2];
        u_next[n - 1] = rhsl / lhsl;
    }

    void boundary_explicit_2t2p(int k) {
        double rhs0 = c_2t2p_0[4] * gamma_0[k + 1] - c_2t2p_0[5] * f[0][k + 1] - c_2t2p_0[2] * u_k[0] - c_2t2p_0[3] * u_prev[0] - c_2t2p_0[1] * u_next[1];
        double lhs0 = c_2t2p_0[0];
        u_next[0] = rhs0 / lhs0;
        double rhsl = c_2t2p_l[4] * gamma_l[k + 1] - c_2t2p_l[5] * f[n - 1][k + 1] - c_2t2p_l[2] * u_k[n - 1] - c_2t2p_l[3] * u_prev[n - 1] - c_2t2p_l[0] * u_next[n - 2];
        double lhsl = c_2t2p_l[1];
        u_next[n - 1] = rhsl / lhsl;
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
    friend std::istream & operator >> (std::istream & in, hpde_t & item) {
        in >> item.n >> item.K >> item.u1_degree >> item.boundary;
        in >> item.l >> item.T >> item.h >> item.tau >> item.theta;
        item.h2 = item.h * item.h;
        item.tau2 = item.tau * item.tau;
        in >> item.a >> item.b >> item.c >> item.d;
        if (item.a < 0) {
            throw std::invalid_argument("a is invalid!");
        }
        item.f.resize(item.n, vec(item.K));
        for (int i = 0; i < item.n; ++i) {
            for (int j = 0; j < item.K; ++j) {
                in >> item.f[i][j];
            }
        }
        item.u_prev.resize(item.n);
        for (int i = 0; i < item.n; ++i) {
            in >> item.u_prev[i];
        }
        item.u_k.resize(item.n);
        for (int i = 0; i < item.n; ++i) {
            in >> item.u_k[i];
        }

        in >> item.alpha_0 >> item.beta_0 >> item.alpha_l >> item.beta_l;
        item.gamma_0.resize(item.K);
        item.gamma_l.resize(item.K);
        for (int k = 0; k < item.K; ++k) {
            in >> item.gamma_0[k];
        }
        for (int k = 0; k < item.K; ++k) {
            in >> item.gamma_l[k];
        }
        if (item.u1_degree == 2) {
            item.ddx_psi1.resize(item.n);
            item.d2dx2_psi1.resize(item.n);
            for (int i = 0; i < item.n; ++i) {
                in >> item.ddx_psi1[i];
            }
            for (int i = 0; i < item.n; ++i) {
                in >> item.d2dx2_psi1[i];
            }
        }
        return in;
    }

    void solve(std::ostream & out) {
        print_vec(out, u_prev);
        prepare_u1();
        prepare_3t2p();
        prepare_2t2p();
        if (theta == 0) {
            solve_explicit(out);
        } else {
            solve_implicit(out);
        }
    }

    void solve_implicit(std::ostream & out) {
        precalc_implicit();
        print_vec(out, u_k);
        for (int k = 1; k < K - 1; ++k) {
            gen_implicit(k);
            tridiag trd(trd_a, trd_b, trd_c);
            u_next = trd.solve(trd_d);
            print_vec(out, u_next);
            u_prev = u_k;
            u_k = u_next;
        }
    }

    void solve_explicit(std::ostream & out) {
        precalc_explicit();
        print_vec(out, u_k);
        for (int k = 1; k < K - 1; ++k) {
            u_next.assign(n, 0);
            for (int i = 1; i < n - 1; ++i) {
                double rhs = c_expl[1] * u_k[i - 1] + c_expl[2] * u_k[i] + c_expl[3] * u_k[i + 1] + c_expl[4] * u_prev[i] + c_expl[5] * f[i][k];
                double lhs = c_expl[0];
                u_next[i] = rhs / lhs;
            }
            boundary_explicit(k);
            print_vec(out, u_next);
            u_prev = u_k;
            u_k = u_next;
        }
    }
};

#endif /* HPDE_HPP */
