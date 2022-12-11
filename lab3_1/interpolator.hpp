#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP

#include "../polynom.hpp"

using vec = std::vector<double>;

class inter_lagrange {
    vec x;
    vec y;
    size_t n;

   public:
    inter_lagrange(const vec& _x, const vec& _y) : x(_x), y(_y), n(x.size()){};

    polynom operator()() {
        polynom res(vec({0}));
        for (size_t i = 0; i < n; ++i) {
            polynom li(vec({1}));
            for (size_t j = 0; j < n; ++j) {
                if (i == j) {
                    continue;
                }
                polynom xij(vec({-x[j], 1}));
                li = li * xij;
                li = li / (x[i] - x[j]);
            }
            res = res + y[i] * li;
        }
        return res;
    }
};

class inter_newton {
   private:
    using vvd = std::vector<std::vector<double> >;
    using vvb = std::vector<std::vector<bool> >;

    vec x;
    vec y;
    size_t n;

    vvd memo;
    vvb calc;

    double f(int l, int r) {
        if (calc[l][r]) {
            return memo[l][r];
        }
        calc[l][r] = true;
        double res;
        if (l + 1 == r) {
            res = (y[l] - y[r]) / (x[l] - x[r]);
        } else {
            res = (f(l, r - 1) - f(l + 1, r)) / (x[l] - x[r]);
        }
        return memo[l][r] = res;
    }

   public:
    inter_newton(const vec& _x, const vec& _y) : x(_x), y(_y), n(x.size()) {
        memo.resize(n, std::vector<double>(n));
        calc.resize(n, std::vector<bool>(n));
    };

    polynom operator()() {
        polynom res(vec({y[0]}));
        polynom li(vec({-x[0], 1}));
        int r = 0;
        for (size_t i = 1; i < n; ++i) {
            res = res + f(0, ++r) * li;
            li = li * polynom(vec({-x[i], 1}));
        }
        return res;
    }
};

#endif /* INTERPOLATOR_HPP */
