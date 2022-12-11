#ifndef POLYNOM_HPP
#define POLYNOM_HPP

#include <iostream>
#include <vector>

class polynom {
   private:
    using vec = std::vector<double>;

    vec data;
    size_t n;

    constexpr static double EPS = 1e-9;

   public:
    polynom() : data(1), n(1) {}

    polynom(int _n) : data(_n), n(_n) {}

    polynom(const vec& coef) : data(coef), n(data.size()) {}

    size_t size() const { return n; }

    double& operator[](size_t id) { return data[id]; }

    const double& operator[](size_t id) const { return data[id]; }

    friend polynom operator+(const polynom& lhs, const polynom& rhs) {
        polynom res(std::max(lhs.size(), rhs.size()));
        for (size_t i = 0; i < lhs.size(); ++i) {
            res[i] += lhs[i];
        }
        for (size_t i = 0; i < rhs.size(); ++i) {
            res[i] += rhs[i];
        }
        return res;
    }

    friend polynom operator-(const polynom& lhs, const polynom& rhs) {
        polynom res(std::max(lhs.size(), rhs.size()));
        for (size_t i = 0; i < lhs.size(); ++i) {
            res[i] += lhs[i];
        }
        for (size_t i = 0; i < rhs.size(); ++i) {
            res[i] -= rhs[i];
        }
        return res;
    }

    friend polynom operator*(double lambda, const polynom& p) {
        polynom res(p);
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] *= lambda;
        }
        return res;
    }

    friend polynom operator/(const polynom& p, double lambda) {
        polynom res(p);
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] /= lambda;
        }
        return res;
    }

    /* FFT? */
    friend polynom operator*(const polynom& lhs, const polynom& rhs) {
        polynom res(lhs.size() + rhs.size());
        for (size_t i = 0; i < lhs.size(); ++i) {
            for (size_t j = 0; j < rhs.size(); ++j) {
                res[i + j] += lhs[i] * rhs[j];
            }
        }
        while (res.n > 1 and std::abs(res.data.back()) < EPS) {
            res.data.pop_back();
            --res.n;
        }
        return res;
    }

    polynom derivative() {
        polynom res(n - 1);
        for (size_t i = 1; i < n; ++i) {
            res[i - 1] = data[i] * i;
        }
        return res;
    }

    polynom integrate() {
        polynom res(n + 1);
        for (size_t i = 1; i < n + 1; ++i) {
            res.data[i] = data[i - 1] / (double)i;
        }
        return res;
    }

    double integrate(double l, double r) {
        polynom F = integrate();
        return F(r) - F(l);
    }

    double operator()(double x) {
        double res = 0.0;
        double xi = 1.0;
        for (double elem : data) {
            res += elem * xi;
            xi *= x;
        }
        return res;
    }

    friend std::ostream& operator<<(std::ostream& out, const polynom& poly) {
        bool flag = false;
        int deg = 0;
        for (double elem : poly.data) {
            if (!(std::abs(elem) < EPS)) {
                if (flag and deg) {
                    out << (elem > EPS ? " + " : " - ");
                }
                out << std::abs(elem);
                flag = true;
                if (deg) {
                    out << " * x";
                    if (deg > 1) {
                        out << " ^ " << deg;
                    }
                }
            }
            ++deg;
        }
        if (!flag) {
            out << 0;
        }
        return out;
    }

    ~polynom() = default;
};

#endif /* POLYNOM_HPP */
