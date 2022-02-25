#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>

template<class T>
class matrix_t {
private:
    std::vector< std::vector<T> > _data;
    size_t _size;
public:
    matrix_t(size_t n, bool identity = false) : _size(n) {
        _data.resize(_size, std::vector<T>(_size));
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

    friend matrix_t<T> operator * (const matrix_t<T> & a, const matrix_t<T> & b) {
        size_t n = a.size();
        matrix_t<T> res(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                for (size_t k = 0; k < n; ++k) {
                    res._data[i][j] += a._data[i][k] * b._data[k][j];
                }
            }
        }
        return res;
    }

    std::vector<T> & operator [] (size_t i) {
        return _data[i];
    }

    friend std::ostream & operator << (std::ostream & out, const matrix_t<T> & matr) {
        for (size_t i = 0; i < matr.size(); ++i) {
            for (size_t j = 0; j < matr.size(); ++j) {
                if (j) {
                    out << ", ";
                }
                out << matr._data[i][j];
            }
            out << '\n';
        }
        return out;
    }

    friend std::istream & operator >> (std::istream & in, matrix_t<T> & matr) {
        for (size_t i = 0; i < matr.size(); ++i) {
            for (size_t j = 0; j < matr.size(); ++j) {
                in >> matr._data[i][j];
            }
        }
        return in;
    }

    ~matrix_t() = default;
};

#endif /* MATRIX_HPP */
