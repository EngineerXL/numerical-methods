#ifndef TRIDIAG_HPP
#define TRIDIAG_HPP

#include <exception>
#include <iostream>
#include <vector>

template <class T>
class tridiag_t {
   private:
    using vec = std::vector<T>;

    const double EPS = 1e-9;

    int n;
    vec a;
    vec b;
    vec c;

   public:
    tridiag_t(const int& _n) : n(_n), a(n), b(n), c(n) {}

    tridiag_t(const vec& _a, const vec& _b, const vec& _c) {
        if (!(_a.size() == _b.size() and _a.size() == _c.size())) {
            throw std::invalid_argument("Sizes of a, b, c are invalid");
        }
        n = _a.size();
        a = _a;
        b = _b;
        c = _c;
    }

    friend std::istream& operator>>(std::istream& in, tridiag_t<T>& tridiag) {
        in >> tridiag.b[0] >> tridiag.c[0];
        for (int i = 1; i < tridiag.n - 1; ++i) {
            in >> tridiag.a[i] >> tridiag.b[i] >> tridiag.c[i];
        }
        in >> tridiag.a.back() >> tridiag.b.back();
        return in;
    }

    vec solve(const vec& d) {
        int m = d.size();
        if (n != m) {
            throw std::invalid_argument("Size of vector d is invalid");
        }
        vec p(n);
        p[0] = -c[0] / b[0];
        vec q(n);
        q[0] = d[0] / b[0];
        for (int i = 1; i < n; ++i) {
            p[i] = -c[i] / (b[i] + a[i] * p[i - 1]);
            q[i] = (d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]);
        }
        vec x(n);
        x.back() = q.back();
        for (int i = n - 2; i >= 0; --i) {
            x[i] = p[i] * x[i + 1] + q[i];
        }
        return x;
    }

    ~tridiag_t() = default;
};

#endif /* TRIDIAG_HPP */
