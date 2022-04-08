#ifndef INTERPOLATOR_HPP
#define INTERPOLATOR_HPP

#include "../polynom.hpp"

using vec = std::vector<double>;

class inter_lagrange {
    vec x;
    vec y;
    size_t n;

public:
    inter_lagrange(const vec & _x, const vec & _y) : x(_x), y(_y), n(x.size()) {};

    polynom operator () () {
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
    vec x;
    vec y;
    size_t n;
    using vint = std::vector<int>;

    double f(vint bits) {
        std::vector<size_t> ids;
        for (size_t i = 0; i < n; ++i) {
            if (bits[i]) {
                ids.push_back(i);
            }
        }
        size_t id0 = ids[0], idn = ids.back();
        if (ids.size() > 2) {
            vint l_bits(bits), r_bits(bits);
            l_bits[id0] ^= 1;
            r_bits[idn] ^= 1;
            return (f(r_bits) - f(l_bits)) / (x[id0] - x[idn]);
        } else {
            return (y[id0] - y[idn]) / (x[id0] - x[idn]);
        }
    }

public:
    inter_newton(const vec & _x, const vec & _y) : x(_x), y(_y), n(x.size()) {};

    polynom operator () () {
        polynom res(vec({y[0]}));
        polynom li(vec({-x[0], 1}));
        vint mask(n);
        mask[0] ^= 1;
        for (size_t i = 1; i < n; ++i) {
            mask[i] ^= 1;
            res = res + f(mask) * li;
            li = li * polynom(vec({-x[i], 1}));
        }
        return res;
    }
};

#endif /* INTERPOLATOR_HPP */
