#ifndef SIMPLE_DESOLVE_HPP
#define SIMPLE_DESOLVE_HPP

#include "../de_utils.hpp"
#include <functional>

/* f(x, y, z) */
using func = std::function<double(double, double, double)>;
using vect = std::vector<tddd>;
using vec = std::vector<double>;

const double EPS = 1e-9;

bool leq(double a, double b) {
    return (a < b) or (std::abs(b - a) < EPS);
}

class euler {
private:
    double l, r;
    func f, g;
    double y0, z0;

public:
    euler(const double _l, const double _r,
        const func _f, const func _g,
        const double _y0, const double _z0) : l(_l), r(_r), f(_f), g(_g), y0(_y0), z0(_z0) {}

    vect solve(double h) {
        vect res;
        double xk = l;
        double yk = y0;
        double zk = z0;
        res.push_back(std::make_tuple(xk, yk, zk));
        while (leq(xk + h, r)) {
            double dy = h * f(xk, yk, zk);
            double dz = h * g(xk, yk, zk);
            xk += h;
            yk += dy;
            zk += dz;
            res.push_back(std::make_tuple(xk, yk, zk));
        }
        return res;
    }
};

class runge {
private:
    double l, r;
    func f, g;
    double y0, z0;

public:
    runge(const double _l, const double _r,
        const func _f, const func _g,
        const double _y0, const double _z0) : l(_l), r(_r), f(_f), g(_g), y0(_y0), z0(_z0) {}

    vect solve(double h) {
        vect res;
        double xk = l;
        double yk = y0;
        double zk = z0;
        res.push_back(std::make_tuple(xk, yk, zk));
        while (leq(xk + h, r)) {
            double K1 = h * f(xk, yk, zk);
            double L1 = h * g(xk, yk, zk);
            double K2 = h * f(xk + 0.5 * h, yk + 0.5 * K1, zk + 0.5 * L1);
            double L2 = h * g(xk + 0.5 * h, yk + 0.5 * K1, zk + 0.5 * L1);
            double K3 = h * f(xk + 0.5 * h, yk + 0.5 * K2, zk + 0.5 * L2);
            double L3 = h * g(xk + 0.5 * h, yk + 0.5 * K2, zk + 0.5 * L2);
            double K4 = h * f(xk + h, yk + K3, zk + L3);
            double L4 = h * g(xk + h, yk + K3, zk + L3);
            double dy = (K1 + 2.0 * K2 + 2.0 * K3 + K4) / 6.0;
            double dz = (L1 + 2.0 * L2 + 2.0 * L3 + L4) / 6.0;
            xk += h;
            yk += dy;
            zk += dz;
            res.push_back(std::make_tuple(xk, yk, zk));
        }
        return res;
    }
};

class adams {
private:
    double l, r;
    func f, g;
    double y0, z0;

public:
    adams(const double _l, const double _r,
        const func _f, const func _g,
        const double _y0, const double _z0) : l(_l), r(_r), f(_f), g(_g), y0(_y0), z0(_z0) {}

    double calc_tuple(func f, tddd xyz) {
        return f(std::get<0>(xyz), std::get<1>(xyz), std::get<2>(xyz));
    }

    vect solve(double h) {
        if (l + 3.0 * h > r) {
            throw std::invalid_argument("h is too big");
        }
        runge first_points(l, l + 3.0 * h, f, g, y0, z0);
        vect res = first_points.solve(h);
        size_t cnt = res.size();
        double xk = std::get<0>(res.back());
        double yk = std::get<1>(res.back());
        double zk = std::get<2>(res.back());
        while (leq(xk + h, r)) {
            /* Predictor */
            double dy = (h / 24.0) * (55.0 * calc_tuple(f, res[cnt - 1])
                                    - 59.0 * calc_tuple(f, res[cnt - 2])
                                    + 37.0 * calc_tuple(f, res[cnt - 3])
                                    - 9.0 * calc_tuple(f, res[cnt - 4]));
            double dz = (h / 24.0) * (55.0 * calc_tuple(g, res[cnt - 1])
                                    - 59.0 * calc_tuple(g, res[cnt - 2])
                                    + 37.0 * calc_tuple(g, res[cnt - 3])
                                    - 9.0 * calc_tuple(g, res[cnt - 4]));
            double xk1 = xk + h;
            double yk1 = yk + dy;
            double zk1 = zk + dz;
            res.push_back(std::make_tuple(xk1, yk1, zk1));
            ++cnt;
            /* Corrector */
            dy = (h / 24.0) * (9.0 * calc_tuple(f, res[cnt - 1])
                            + 19.0 * calc_tuple(f, res[cnt - 2])
                            - 5.0 * calc_tuple(f, res[cnt - 3])
                            + 1.0 * calc_tuple(f, res[cnt - 4]));
            dz = (h / 24.0) * (9.0 * calc_tuple(g, res[cnt - 1])
                            + 19.0 * calc_tuple(g, res[cnt - 2])
                            - 5.0 * calc_tuple(g, res[cnt - 3])
                            + 1.0 * calc_tuple(g, res[cnt - 4]));
            xk += h;
            yk += dy;
            zk += dz;
            res.pop_back();
            res.push_back(std::make_tuple(xk, yk, zk));
        }
        return res;
    }
};

double runge_romberg(const vect & y_2h, const vect & y_h, double p) {
    double coef = 1.0 / (std::pow(2, p) - 1.0);
    double res = 0.0;
    for (size_t i = 0; i < y_2h.size(); ++i) {
        res = std::max(res, coef * std::abs(std::get<1>(y_2h[i]) - std::get<1>(y_h[2 * i])));
    }
    return res;
}

#endif /* SIMPLE_DESOLVE_HPP */
