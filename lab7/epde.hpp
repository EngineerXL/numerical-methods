#ifndef EPDE_HPP
#define EPDE_HPP

#include "../matrix.hpp"

/* Elliptic Partial Differential Equation Solver */
class epde_t {
    using vec = std::vector<double>;
    using matr = matrix_t<double>;

    static constexpr double EPS = 1e-6;
    static const int N_COEFS = 7;

    int n, m, mu, iters;
    double omega, hx, hy, hx2, hy2;
    matr u;

    double coef[N_COEFS];
    double a, b, c, d, e, f, g;

    double alpha_0y, alpha_Ly, alpha_x0, alpha_xL;
    double beta_0y, beta_Ly, beta_x0, beta_xL;
    vec gamma_0y, gamma_Ly, gamma_x0, gamma_xL;

    void init() {
        a = coef[0];
        b = coef[1];
        c = coef[2];
        d = coef[3];
        e = coef[4];
        f = coef[5];
        g = coef[6];
        hx2 = hx * hx;
        hy2 = hy * hy;
        u = matr::uniform(n, m);
        gamma_x0.resize(n);
        gamma_xL.resize(n);
        gamma_0y.resize(m);
        gamma_Ly.resize(m);
        prepare_scheme();
    }

    static const int N_SCHEME = 5;
    double coef_sch[N_SCHEME];

    void prepare_scheme() {
        coef_sch[0] = 2 * (a / hx2 + c / hy2) - f;
        coef_sch[1] = c / hy2 - e / (2 * hy);
        coef_sch[2] = c / hy2 + e / (2 * hy);
        coef_sch[3] = a / hx2 - d / (2 * hx);
        coef_sch[4] = a / hx2 + d / (2 * hx);
    }

    void iter(const matr& u_prev, matr& u_next) {
        for (int i = 1; i < n - 1; ++i) {
            for (int j = 1; j < m - 1; ++j) {
                double rhs = coef_sch[1] * u_prev[i][j - 1] +
                             coef_sch[2] * u_prev[i][j + 1] +
                             coef_sch[3] * u_prev[i - 1][j] +
                             coef_sch[4] * u_prev[i + 1][j] + g;
                double lhs = coef_sch[0];
                u_next[i][j] = rhs / lhs;
            }
        }
        calc_boundary(u_next);
    }

    void calc_boundary(matr& u_) {
        for (int i = 0; i < n; ++i) {
            u_[i][0] = 0;
            u_[i][m - 1] = 0;
        }
        for (int j = 0; j < m; ++j) {
            u_[0][j] = 0;
            u_[n - 1][j] = 0;
        }

        for (int i = 0; i < n; ++i) {
            double rhs = gamma_x0[i] - alpha_x0 / hy * u_[i][1];
            double lhs = beta_x0 - alpha_x0 / hy;
            u_[i][0] += rhs / lhs;
        }
        for (int i = 0; i < n; ++i) {
            double rhs = gamma_xL[i] - alpha_xL / hy * u_[i][m - 2];
            double lhs = beta_xL - alpha_xL / hy;
            u_[i][m - 1] += rhs / lhs;
        }

        for (int j = 0; j < m; ++j) {
            double rhs = gamma_0y[j] - alpha_0y / hx * u_[1][j];
            double lhs = beta_0y - alpha_0y / hx;
            u_[0][j] += rhs / lhs;
        }
        for (int j = 0; j < m; ++j) {
            double rhs = gamma_Ly[j] + alpha_Ly / hx * u_[n - 2][j];
            double lhs = beta_Ly + alpha_Ly / hx;
            u_[n - 1][j] += rhs / lhs;
        }

        u_[n - 1][m - 1] /= 2;
        u_[n - 1][0] /= 2;
        u_[0][m - 1] /= 2;
        u_[0][0] /= 2;
    }

    double calc_delta(const matr& u_prev, const matr& u_next) {
        double delta = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                delta = std::max(delta, std::abs(u_prev[i][j] - u_next[i][j]));
            }
        }
        return delta;
    }

   public:
    friend std::istream& operator>>(std::istream& in, epde_t& item) {
        in >> item.n >> item.m >> item.mu >> item.omega;
        in >> item.hx >> item.hy;
        for (int i = 0; i < N_COEFS; ++i) {
            in >> item.coef[i];
        }
        item.init();
        in >> item.alpha_0y >> item.alpha_Ly >> item.alpha_x0 >> item.alpha_xL;
        in >> item.beta_0y >> item.beta_Ly >> item.beta_x0 >> item.beta_xL;

        for (int i = 0; i < item.n; ++i) {
            in >> item.gamma_x0[i];
        }
        for (int i = 0; i < item.n; ++i) {
            in >> item.gamma_xL[i];
        }

        for (int i = 0; i < item.m; ++i) {
            in >> item.gamma_0y[i];
        }
        for (int i = 0; i < item.m; ++i) {
            in >> item.gamma_Ly[i];
        }

        return in;
    }

    void solve_simple() {
        double delta = 1;
        while (delta > EPS) {
            matr u_next(n, m);
            iter(u, u_next);
            delta = calc_delta(u, u_next);
            u = u_next;
            ++iters;
        }
    }

    void solve_zeidel() {
        double delta = 1;
        while (delta > EPS) {
            matr u_next(u);
            iter(u_next, u_next);
            delta = calc_delta(u, u_next);
            u = u_next;
            ++iters;
        }
    }

    void solve_relax() {
        double delta = 1;
        while (delta > EPS) {
            matr u_next(n, m);
            iter(u, u_next);
            u_next = omega * u_next + (1 - omega) * u;
            // u_next = omega * (u_next - u) + u;
            delta = calc_delta(u, u_next);
            // std::cout << delta << std::endl;
            u = u_next;
            ++iters;
        }
    }

    void solve(std::ostream& out) {
        iters = 0;
        if (mu == 1) {
            solve_simple();
        } else if (mu == 2) {
            solve_zeidel();
        } else if (mu == 3) {
            solve_relax();
        } else {
            throw std::runtime_error("Invalid mu");
        }
        u.comma_out(false);
        out << iters << '\n';
        out << u;
    }
};

#endif /* PPDE_HPP */
