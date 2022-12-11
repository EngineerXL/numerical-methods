#ifndef FFT_MULT_HPP
#define FFT_MULT_HPP

#include <algorithm>
#include <complex>
#include <string>
#include <vector>

using complex = std::complex<double>;
using vc = std::vector<complex>;
const double PI = std::acos(-1);

using vi = std::vector<int64_t>;
const int64_t BASE = 10;

int rev_bits(int x, int lg_n) {
    int y = 0;
    for (int i = 0; i < lg_n; ++i) {
        y = y << 1;
        y ^= (x & 1);
        x = x >> 1;
    }
    return y;
}

void fft(vc& a, bool invert) {
    int n = a.size();
    int lg_n = 0;
    while ((1 << lg_n) < n) {
        ++lg_n;
    }
    for (int i = 0; i < n; ++i) {
        if (i < rev_bits(i, lg_n)) {
            swap(a[i], a[rev_bits(i, lg_n)]);
        }
    }
    for (int layer = 1; layer <= lg_n; ++layer) {
        int cluster = 1 << layer;
        double phi = (2.0 * PI) / cluster;
        if (invert) {
            phi *= -1;
        }
        complex wn = complex(std::cos(phi), std::sin(phi));
        for (int i = 0; i < n; i += cluster) {
            complex w(1);
            for (int j = 0; j < cluster / 2; ++j) {
                complex u = a[i + j];
                complex v = a[i + j + cluster / 2] * w;
                a[i + j] = u + v;
                a[i + j + cluster / 2] = u - v;
                w *= wn;
            }
        }
    }
    if (invert) {
        for (int i = 0; i < n; ++i) {
            a[i] /= n;
        }
    }
}

std::string fft_mult(const std::string& a, const std::string& b) {
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
        res[i] = (int64_t)round(fa[i].real());
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

std::string slow_mult(const std::string& a, const std::string& b) {
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
