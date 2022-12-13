#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <exception>
#include <iostream>
#include <vector>

template <class T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b) {
    size_t n = a.size();
    std::vector<T> c(n);
    for (size_t i = 0; i < n; ++i) {
        c[i] = a[i] + b[i];
    }
    return c;
}

template <class T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b) {
    size_t n = a.size();
    std::vector<T> c(n);
    for (size_t i = 0; i < n; ++i) {
        c[i] = a[i] - b[i];
    }
    return c;
}

template <class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
    size_t n = v.size();
    for (size_t i = 0; i < n; ++i) {
        if (i) {
            out << ' ';
        }
        out << v[i];
    }
    return out;
}

template <class T>
class matrix_t {
   private:
    using matrix = matrix_t<T>;
    using vec = std::vector<T>;

    size_t n, m;
    std::vector<vec> data;
    bool comma = true;

   public:
    matrix_t() : n(1), m(1), data(1) {}

    matrix_t(size_t _n) : n(_n), m(_n) { data.resize(n, vec(n)); }

    matrix_t(size_t _n, size_t _m) : n(_n), m(_m) { data.resize(n, vec(m)); }

    matrix_t(const matrix& other) {
        n = other.n;
        m = other.m;
        data = other.data;
        comma = other.comma;
    }

    matrix& operator=(const matrix& other) {
        if (this == &other) {
            return *this;
        }
        n = other.n;
        m = other.m;
        data = other.data;
        comma = other.comma;
        return *this;
    }

    void comma_out(bool flag) { comma = flag; }

    static matrix identity(size_t n) {
        matrix res(n, n);
        for (size_t i = 0; i < n; ++i) {
            res[i][i] = T(1);
        }
        return res;
    }

    size_t rows() const { return n; }

    size_t cols() const { return m; }

    void swap_rows(size_t i, size_t j) {
        if (i == j) {
            return;
        }
        for (size_t k = 0; k < m; ++k) {
            std::swap(data[i][k], data[j][k]);
        }
    }

    void swap_cols(size_t i, size_t j) {
        if (i == j) {
            return;
        }
        for (size_t k = 0; k < n; ++k) {
            std::swap(data[k][i], data[k][j]);
        }
    }

    matrix t() const {
        matrix_t<T> res(m, n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[j][i] = data[i][j];
            }
        }
        return res;
    }

    static matrix uniform(size_t n, size_t m) {
        matrix res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[i][j] = 0.1f;
            }
        }
        return res;
    }

    friend matrix operator+(const matrix& a, const matrix& b) {
        if (a.rows() != b.rows() or a.cols() != b.cols()) {
            throw std::invalid_argument("Sizes does not match");
        }
        size_t n = a.rows();
        size_t m = a.cols();
        matrix res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[i][j] = a[i][j] + b[i][j];
            }
        }
        return res;
    }

    friend matrix operator-(const matrix& a, const matrix& b) {
        if (a.rows() != b.rows() or a.cols() != b.cols()) {
            throw std::invalid_argument("Sizes does not match");
        }
        size_t n = a.rows();
        size_t m = a.cols();
        matrix res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[i][j] = a[i][j] - b[i][j];
            }
        }
        return res;
    }

    friend matrix operator*(T lambda, const matrix& a) {
        size_t n = a.rows();
        size_t m = a.cols();
        matrix res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                res[i][j] = lambda * a[i][j];
            }
        }
        return res;
    }

    friend vec operator*(const matrix& a, const vec& b) {
        if (a.cols() != b.size()) {
            throw std::invalid_argument("Sizes does not match");
        }
        size_t n = a.rows();
        size_t m = a.cols();
        vec c(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                c[i] += a[i][j] * b[j];
            }
        }
        return c;
    }

    friend matrix operator*(const matrix& a, const matrix& b) {
        if (a.cols() != b.rows()) {
            throw std::invalid_argument("Sizes does not match");
        }
        size_t n = a.rows();
        size_t k = a.cols();
        size_t m = b.cols();
        matrix res(n, m);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                for (size_t ii = 0; ii < k; ++ii) {
                    res[i][j] += a[i][ii] * b[ii][j];
                }
            }
        }
        return res;
    }

    vec& operator[](size_t i) { return data[i]; }

    const vec& operator[](size_t i) const { return data[i]; }

    friend std::ostream& operator<<(std::ostream& out, const matrix& matr) {
        for (size_t i = 0; i < matr.rows(); ++i) {
            for (size_t j = 0; j < matr.cols(); ++j) {
                if (j) {
                    out << ' ';
                    if (matr.comma) {
                        out << ',';
                    }
                }
                out << matr[i][j];
            }
            out << '\n';
        }
        return out;
    }

    friend std::istream& operator>>(std::istream& in, matrix& matr) {
        for (size_t i = 0; i < matr.rows(); ++i) {
            for (size_t j = 0; j < matr.cols(); ++j) {
                in >> matr[i][j];
            }
        }
        return in;
    }

    ~matrix_t() = default;
};

#endif /* MATRIX_HPP */
