#ifndef LU_HPP
#define LU_HPP

#include <algorithm>
#include <cmath>
#include <utility>

#include "../matrix.hpp"

template <class T>
class lu_t {
   private:
    using matrix = matrix_t<T>;
    using vec = std::vector<T>;
    using pii = std::pair<size_t, size_t>;

    const T EPS = 1e-6;

    matrix l;
    matrix u;
    T det;
    std::vector<pii> swaps;

    void decompose() {
        size_t n = u.rows();
        for (size_t i = 0; i < n; ++i) {
            size_t max_el_ind = i;
            for (size_t j = i + 1; j < n; ++j) {
                if (abs(u[j][i]) > abs(u[max_el_ind][i])) {
                    max_el_ind = j;
                }
            }
            if (max_el_ind != i) {
                pii perm = std::make_pair(i, max_el_ind);
                swaps.push_back(perm);
                u.swap_rows(i, max_el_ind);
                l.swap_rows(i, max_el_ind);
                l.swap_cols(i, max_el_ind);
            }
            for (size_t j = i + 1; j < n; ++j) {
                if (abs(u[i][i]) < EPS) {
                    continue;
                }
                T mu = u[j][i] / u[i][i];
                l[j][i] = mu;
                for (size_t k = 0; k < n; ++k) {
                    u[j][k] -= mu * u[i][k];
                }
            }
        }
        det = (swaps.size() & 1 ? -1 : 1);
        for (size_t i = 0; i < n; ++i) {
            det *= u[i][i];
        }
    }

    void do_swaps(vec& x) {
        for (pii elem : swaps) {
            std::swap(x[elem.first], x[elem.second]);
        }
    }

   public:
    lu_t(const matrix& matr) {
        if (matr.rows() != matr.cols()) {
            throw std::invalid_argument("Matrix is not square");
        }
        l = matrix::identity(matr.rows());
        u = matrix(matr);
        decompose();
    }

    friend std::ostream& operator<<(std::ostream& out, const lu_t<T>& lu) {
        out << "Matrix L:\n" << lu.l << "Matrix U:\n" << lu.u;
        return out;
    }

    T get_det() { return det; }

    vec solve(vec b) {
        int n = b.size();
        do_swaps(b);
        vec z(n);
        for (int i = 0; i < n; ++i) {
            T summary = b[i];
            for (int j = 0; j < i; ++j) {
                summary -= z[j] * l[i][j];
            }
            z[i] = summary;
        }
        vec x(n);
        for (int i = n - 1; i >= 0; --i) {
            if (abs(u[i][i]) < EPS) {
                continue;
            }
            T summary = z[i];
            for (int j = n - 1; j > i; --j) {
                summary -= x[j] * u[i][j];
            }
            x[i] = summary / u[i][i];
        }
        return x;
    }

    matrix inv_matrix() {
        size_t n = l.rows();
        matrix res(n);
        for (size_t i = 0; i < n; ++i) {
            vec b(n);
            b[i] = T(1);
            vec x = solve(b);
            for (size_t j = 0; j < n; ++j) {
                res[j][i] = x[j];
            }
        }
        return res;
    }

    ~lu_t() = default;
};

#endif /* LU_HPP */
