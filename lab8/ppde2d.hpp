#ifndef PPDE2D_HPP
#define PPDE2D_HPP

#include "../lab1_2/tridiag.hpp"
#include "tensor3.hpp"

/* 2-D Parabolic Partial Differential Equation Solver */
class ppde2d_t {
    using vec = std::vector<double>;
    using tensor3d = tensor3_t<double>;
    using tridiag = tridiag_t<double>;

    int nx, ny, K, psi;
    double lx, ly, T, hx, hy, tau;
    double ax, ay, bx, by, c;
    double alpha_0_y, beta_0_y, alpha_x_0, beta_x_0, alpha_lx_y, beta_lx_y,
        alpha_x_ly, beta_x_ly;
    double hx2, hy2;

    tensor3d f, u, u_abab, u0, gamma_0_y, gamma_x_0, gamma_lx_y, gamma_x_ly;
    vec trd_a_x, trd_b_x, trd_c_x, trd_d_x;
    vec trd_a_y, trd_b_y, trd_c_y, trd_d_y;

    void init() {
        hx2 = hx * hx;
        hy2 = hy * hy;
        u0 = tensor3d(nx, ny, 1, 0.5);
        f = tensor3d(nx, ny, K);
        gamma_0_y = tensor3d(1, ny, K);
        gamma_x_0 = tensor3d(nx, 1, K);
        gamma_lx_y = tensor3d(1, ny, K);
        gamma_x_ly = tensor3d(nx, 1, K);
        trd_a_x.resize(nx);
        trd_b_x.resize(nx);
        trd_c_x.resize(nx);
        trd_d_x.resize(nx);
        trd_a_y.resize(ny);
        trd_b_y.resize(ny);
        trd_c_y.resize(ny);
        trd_d_y.resize(ny);
    }

    vec coefs_x_1, coefs_y_1;
    void prepare_coefs_1() {
        coefs_x_1.resize(6, 0);
        coefs_x_1[0] = ax / hx2 - bx / (2 * hx);
        coefs_x_1[1] = -2 * ax / hx2 - 2 / tau + c;
        coefs_x_1[2] = ax / hx2 + bx / (2 * hx);
        coefs_x_1[3] = -ay / hy2 + by / (2 * hy);
        coefs_x_1[4] = -2 / tau + 2 * ay / hy2;
        coefs_x_1[5] = -ay / hy2 - by / (2 * hy);

        coefs_y_1.resize(6, 0);
        coefs_y_1[0] = ay / hy2 - by / (2 * hy);
        coefs_y_1[1] = -2 * ay / hy2 - 2 / tau + c;
        coefs_y_1[2] = ay / hy2 + by / (2 * hy);
        coefs_y_1[3] = -ax / hx2 + bx / (2 * hx);
        coefs_y_1[4] = 2 * ax / hx2 - 2 / tau;
        coefs_y_1[5] = -ax / hx2 - bx / (2 * hx);
    }

    vec coefs_x_2, coefs_y_2;
    void prepare_coefs_2() {
        coefs_x_2.resize(4, 0);
        coefs_x_2[0] = ax / hx2 - bx / (2 * hx);
        coefs_x_2[1] = -2 * ax / hx2 - 1 / tau + c;
        coefs_x_2[2] = ax / hx2 + bx / (2 * hx);
        coefs_x_2[3] = -1 / tau;

        coefs_y_2.resize(4, 0);
        coefs_y_2[0] = ay / hy2 - by / (2 * hy);
        coefs_y_2[1] = -2 * ay / hy2 - 1 / tau + c;
        coefs_y_2[2] = ay / hy2 + by / (2 * hy);
        coefs_y_2[3] = -1 / tau;
    }

    vec bound_coefs_x_0;
    void prepare_bound_x_0() {
        bound_coefs_x_0.resize(2);
        bound_coefs_x_0[0] = beta_x_0 - alpha_x_0 / hy;
        bound_coefs_x_0[1] = alpha_x_0 / hy;
    }

    vec bound_coefs_x_ly;
    void prepare_bound_x_ly() {
        bound_coefs_x_ly.resize(2);
        bound_coefs_x_ly[0] = -alpha_x_ly / hy;
        bound_coefs_x_ly[1] = beta_x_ly + alpha_x_ly / hy;
    }

    vec bound_coefs_0_y;
    void prepare_bound_0_y() {
        bound_coefs_0_y.resize(2);
        bound_coefs_0_y[0] = beta_0_y - alpha_0_y / hx;
        bound_coefs_0_y[1] = alpha_0_y / hx;
    }

    vec bound_coefs_lx_y;
    void prepare_bound_lx_y() {
        bound_coefs_lx_y.resize(2);
        bound_coefs_lx_y[0] = -alpha_lx_y / hx;
        bound_coefs_lx_y[1] = beta_lx_y + alpha_lx_y / hx;
    }

    void prepare_bound() {
        prepare_bound_0_y();
        prepare_bound_lx_y();
        prepare_bound_x_0();
        prepare_bound_x_ly();
    }

   public:
    friend std::istream& operator>>(std::istream& in, ppde2d_t& item) {
        in >> item.nx >> item.ny >> item.K >> item.psi;
        in >> item.lx >> item.ly >> item.T >> item.hx >> item.hy >> item.tau;
        in >> item.ax >> item.ay >> item.bx >> item.by >> item.c;
        item.init();
        in >> item.f >> item.u0;
        in >> item.alpha_0_y >> item.alpha_x_0 >> item.beta_0_y >>
            item.beta_x_0;
        in >> item.alpha_lx_y >> item.alpha_x_ly >> item.beta_lx_y >>
            item.beta_x_ly;
        in >> item.gamma_0_y >> item.gamma_x_0;
        in >> item.gamma_lx_y >> item.gamma_x_ly;
        return in;
    }

    void trd_x_set_bounds(const int j, const int k) {
        trd_b_x[0] = bound_coefs_0_y[0];
        trd_c_x[0] = bound_coefs_0_y[1];
        /* Нужно брать среднее, так как шагаем на половину шага */
        trd_d_x[0] = (gamma_0_y(0, j, k) + gamma_0_y(0, j, k + 1)) / 2;
        trd_a_x[nx - 1] = bound_coefs_lx_y[0];
        trd_b_x[nx - 1] = bound_coefs_lx_y[1];
        /* Нужно брать среднее, так как шагаем на половину шага */
        trd_d_x[nx - 1] = (gamma_lx_y(0, j, k) + gamma_lx_y(0, j, k + 1)) / 2;
    }

    void trd_y_set_bounds(const int i, const int k) {
        trd_b_y[0] = bound_coefs_x_0[0];
        trd_c_y[0] = bound_coefs_x_0[1];
        trd_d_y[0] = gamma_x_0(i, 0, k + 1);
        trd_a_y[ny - 1] = bound_coefs_x_ly[0];
        trd_b_y[ny - 1] = bound_coefs_x_ly[1];
        trd_d_y[ny - 1] = gamma_x_ly(i, 0, k + 1);
    }

    void calc_bounds_x(tensor3d& u_next, const int k) {
        for (int i = 0; i < nx; ++i) {
            /* Нужно брать среднее, так как шагаем на половину шага */
            double rhs_x_0 = (gamma_x_0(i, 0, k) + gamma_x_0(i, 0, k + 1)) / 2 -
                             bound_coefs_x_0[1] * u_next(i, 1, 0);
            double lhs_x_0 = bound_coefs_x_0[0];
            u_next(i, 0, 0) = rhs_x_0 / lhs_x_0;
            /* Нужно брать среднее, так как шагаем на половину шага */
            double rhs_x_ly =
                (gamma_x_ly(i, 0, k) + gamma_x_ly(i, 0, k + 1)) / 2 -
                bound_coefs_x_ly[0] * u_next(i, ny - 2, 0);
            double lhs_x_ly = bound_coefs_x_ly[1];
            u_next(i, ny - 1, 0) = rhs_x_ly / lhs_x_ly;
        }
    }

    void calc_bounds_y(tensor3d& u_next, const int k) {
        for (int j = 0; j < ny; ++j) {
            double rhs_0_y =
                gamma_0_y(0, j, k + 1) - bound_coefs_0_y[1] * u_next(1, j, 0);
            double lhs_0_y = bound_coefs_0_y[0];
            u_next(0, j, 0) = rhs_0_y / lhs_0_y;

            double rhs_lx_y = gamma_lx_y(0, j, k + 1) -
                              bound_coefs_lx_y[0] * u_next(nx - 2, j, 0);
            double lhs_lx_y = bound_coefs_lx_y[1];
            u_next(nx - 1, j, 0) = rhs_lx_y / lhs_lx_y;
        }
    }

    void solve(std::ostream& out) {
        prepare_bound();
        prepare_coefs_1();
        prepare_coefs_2();
        u = u0;
        u_abab = u0;
        out << u << '\n';
        for (int k = 0; k < K - 1; ++k) {
            if (psi == 1) {
                solve_1_x(u, u_abab, k);
                solve_1_y(u_abab, u, k);
            } else {
                solve_2_x(u, u_abab, k);
                solve_2_y(u_abab, u, k);
            }
            out << u << '\n';
        }
    }

    /* k -> k + 1/2 */
    void solve_1_x(const tensor3d& u_prev, tensor3d& u_next, const int k) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                trd_a_x[i] = coefs_x_1[0];
                trd_b_x[i] = coefs_x_1[1];
                trd_c_x[i] = coefs_x_1[2];
                /* Нужно брать среднее, так как шагаем на половину шага */
                trd_d_x[i] = coefs_x_1[3] * u_prev(i, j - 1, 0) +
                             coefs_x_1[4] * u_prev(i, j, 0) +
                             coefs_x_1[5] * u_prev(i, j + 1, 0) -
                             (f(i, j, k) + f(i, j, k + 1)) / 2;
            }
            trd_x_set_bounds(j, k);
            tridiag trd(trd_a_x, trd_b_x, trd_c_x);
            vec res = trd.solve(trd_d_x);
            for (int i = 0; i < nx; ++i) {
                u_next(i, j, 0) = res[i];
            }
        }
        calc_bounds_x(u_next, k);
    }

    /* k + 1/2 -> k + 1 */
    void solve_1_y(const tensor3d& u_prev, tensor3d& u_next, const int k) {
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                trd_a_y[j] = coefs_y_1[0];
                trd_b_y[j] = coefs_y_1[1];
                trd_c_y[j] = coefs_y_1[2];
                trd_d_y[j] = coefs_y_1[3] * u_prev(i - 1, j, 0) +
                             coefs_y_1[4] * u_prev(i, j, 0) +
                             coefs_y_1[5] * u_prev(i + 1, j, 0) -
                             f(i, j, k + 1);
            }
            trd_y_set_bounds(i, k);
            tridiag trd(trd_a_y, trd_b_y, trd_c_y);
            vec res = trd.solve(trd_d_y);
            for (int j = 0; j < ny; ++j) {
                u_next(i, j, 0) = res[j];
            }
        }
        calc_bounds_y(u_next, k);
    }

    /* k -> k + 1/2 */
    void solve_2_x(const tensor3d& u_prev, tensor3d& u_next, const int k) {
        for (int j = 1; j < ny - 1; ++j) {
            for (int i = 1; i < nx - 1; ++i) {
                trd_a_x[i] = coefs_x_2[0];
                trd_b_x[i] = coefs_x_2[1];
                trd_c_x[i] = coefs_x_2[2];
                /* Нужно брать среднее, так как шагаем на половину шага */
                trd_d_x[i] = coefs_x_2[3] * u_prev(i, j, 0) -
                             0.5 * (f(i, j, k) + f(i, j, k + 1)) / 2;
            }
            trd_x_set_bounds(j, k);
            tridiag trd(trd_a_x, trd_b_x, trd_c_x);
            vec res = trd.solve(trd_d_x);
            for (int i = 0; i < nx; ++i) {
                u_next(i, j, 0) = res[i];
            }
        }
        calc_bounds_x(u_next, k);
    }

    /* k + 1/2 -> k + 1 */
    void solve_2_y(const tensor3d& u_prev, tensor3d& u_next, const int k) {
        for (int i = 1; i < nx - 1; ++i) {
            for (int j = 1; j < ny - 1; ++j) {
                trd_a_y[j] = coefs_y_2[0];
                trd_b_y[j] = coefs_y_2[1];
                trd_c_y[j] = coefs_y_2[2];
                trd_d_y[j] =
                    coefs_y_2[3] * u_prev(i, j, 0) - 0.5 * f(i, j, k + 1);
            }
            trd_y_set_bounds(i, k);
            tridiag trd(trd_a_y, trd_b_y, trd_c_y);
            vec res = trd.solve(trd_d_y);
            for (int j = 0; j < ny; ++j) {
                u_next(i, j, 0) = res[j];
            }
        }
        calc_bounds_y(u_next, k);
    }
};

#endif /* PPDE2D_HPP */
