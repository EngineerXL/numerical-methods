#ifndef FFT_MULT_HPP
#define FFT_MULT_HPP

#include <algorithm>
#include <complex>
#include <string>
#include <vector>

using complex = std::complex<double>;
using vc = std::vector<complex>;
using vi = std::vector<int64_t>;

const int64_t BASE = 10;
const double PI = std::acos(-1);

void fft(vc & a, bool invert) {
    size_t n = a.size();
    if (n == 1) {
        return;
    }
    vc a0(n / 2), a1(n / 2);
    for (size_t i = 0, j = 0; j < n; ++i, j += 2) {
        a0[i] = a[j];
        a1[i] = a[j + 1];
    }
    fft(a0, invert);
    fft(a1, invert);
    double phi = 2.0 * PI / n;
    if (invert) {
        phi = -phi;
    }
    complex w(1), wn(std::cos(phi), std::sin(phi));
    for (size_t i = 0; i < n / 2; ++i) {
        a[i] = a0[i] + w * a1[i];
        a[n / 2 + i] = a0[i] - w * a1[i];
        if (invert) {
            a[i] /= 2.0;
            a[n / 2 + i] /= 2.0;
        }
        w *= wn;
    }
}

std::string fft_mult(const std::string & a, const std::string & b) {
    size_t max_size = std::max(a.size(), b.size());
    size_t n = 1;
    while (n < max_size) {
        n *= 2;
    }
    n *= 2;
    vc fa(n), fb(n);
    for (size_t i = 0; i < a.size(); ++i) {
        fa[a.size() - i - 1] = complex(a[i] - '0');
    }
    for (size_t i = 0; i < b.size(); ++i) {
        fb[b.size() - i - 1] = complex(b[i] - '0');
    }
    fft(fa, false);
    fft(fb, false);
    for (size_t i = 0; i < n; ++i) {
        fa[i] = fa[i] * fb[i];
    }
    fft(fa, true);
    vi res(n);
    for (size_t i = 0; i < n; ++i) {
        res[i] = (int64_t)(fa[i].real() + 0.5);
    }
    for (size_t i = 0; i < n - 1; ++i) {
        res[i + 1] += res[i] / BASE;
        res[i] %= BASE;
    }
    while (res.size() > 1 and res.back() == 0) {
        res.pop_back();
    }
    reverse(res.begin(), res.end());
    std::string ab;
    for (int64_t elem : res) {
        ab.push_back('0' + elem);
    }
    return ab;
}

std::string slow_mult(const std::string & a, const std::string & b) {
    size_t n = 4 * std::max(a.size(), b.size());
    vi res(n);
    for (size_t i = 0; i < a.size(); ++i) {
        int64_t ai = a[i] - '0';
        for (size_t j = 0; j < b.size(); ++j) {
            int64_t bj = b[j] - '0';
            res[(b.size() - 1 - j) + (a.size() - 1 - i)] += ai * bj;
        }
    }
    for (size_t i = 0; i < n - 1; ++i) {
        res[i + 1] += res[i] / BASE;
        res[i] %= BASE;
    }
    while (res.size() > 1 and res.back() == 0) {
        res.pop_back();
    }
    reverse(res.begin(), res.end());
    std::string ab;
    for (int64_t elem : res) {
        ab.push_back('0' + elem);
    }
    return ab;
}

#endif /* FFT_MULT_HPP */
