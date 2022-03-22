#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>

template<class T>
std::vector<T> operator + (const std::vector<T> & a, const std::vector<T> & b) {
    size_t n = a.size();
    std::vector<T> c(n);
    for (size_t i = 0; i < n; ++i) {
        c[i] = a[i] + b[i];
    }
    return c;
}

template<class T>
std::vector<T> operator - (const std::vector<T> & a, const std::vector<T> & b) {
    size_t n = a.size();
    std::vector<T> c(n);
    for (size_t i = 0; i < n; ++i) {
        c[i] = a[i] - b[i];
    }
    return c;
}

template<class T>
class matrix_t {
private:
    using matrix = matrix_t<T>;
    using vec = std::vector<T>;

    std::vector<vec> _data;
    size_t _size;
public:
    matrix_t(size_t n, bool identity = false) : _size(n) {
        _data.resize(_size, vec(_size));
        if (identity) {
            for (size_t i = 0; i < _size; ++i) {
                _data[i][i] = 1;
            }
        }
    }

    size_t size() const {
        return _size;
    }

    void swap_rows(size_t i, size_t j) {
        if (i == j) {
            return;
        }
        for (size_t k = 0; k < _size; ++k) {
            std::swap(_data[i][k], _data[j][k]);
        }
    }

    void swap_cols(size_t i, size_t j) {
        if (i == j) {
            return;
        }
        for (size_t k = 0; k < _size; ++k) {
            std::swap(_data[k][i], _data[k][j]);
        }
    }

    matrix t() {
        matrix_t<T> res(*this);
        for (size_t i = 0; i < _size; ++i) {
            for (size_t j = i + 1; j < _size; ++j) {
                std::swap(res[i][j], res[j][i]);
            }
        }
        return res;
    }

    friend matrix operator + (const matrix & a, const matrix & b) {
        size_t n = a.size();
        matrix res(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                res[i][j] = a[i][j] + b[i][j];
            }
        }
        return res;
    }

    friend matrix operator - (const matrix & a, const matrix & b) {
        size_t n = a.size();
        matrix res(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                res[i][j] = a[i][j] - b[i][j];
            }
        }
        return res;
    }

    friend matrix operator * (T lambda, const matrix & a) {
        size_t n = a.size();
        matrix res(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                res[i][j] = lambda * a[i][j];
            }
        }
        return res;
    }

    friend vec operator * (const matrix & a, const vec & b) {
        size_t n = a.size();
        vec c(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                c[i] += a[i][j] * b[j];
            }
        }
        return c;
    }

    friend matrix operator * (const matrix & a, const matrix & b) {
        size_t n = a.size();
        matrix res(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                for (size_t k = 0; k < n; ++k) {
                    res[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return res;
    }

    vec & operator [] (size_t i) {
        return _data[i];
    }

    const vec & operator [] (size_t i) const {
        return _data[i];
    }

    friend std::ostream & operator << (std::ostream & out, const matrix & matr) {
        for (size_t i = 0; i < matr.size(); ++i) {
            for (size_t j = 0; j < matr.size(); ++j) {
                if (j) {
                    out << ", ";
                }
                out << matr[i][j];
            }
            out << '\n';
        }
        return out;
    }

    friend std::istream & operator >> (std::istream & in, matrix & matr) {
        for (size_t i = 0; i < matr.size(); ++i) {
            for (size_t j = 0; j < matr.size(); ++j) {
                in >> matr[i][j];
            }
        }
        return in;
    }

    ~matrix_t() = default;
};

#endif /* MATRIX_HPP */
