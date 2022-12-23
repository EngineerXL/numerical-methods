#ifndef TENSOR3_HPP
#define TENSOR3_HPP

#include <cassert>
#include <iostream>

template <class T>
class tensor3_t {
   private:
    size_t nx, ny, nz;
    T* data;

    using tensor = tensor3_t<T>;

   public:
    tensor3_t() {
        nx = 1;
        ny = 1;
        nz = 1;
        data = new T[nx * ny * nz];
        for (size_t i = 0; i < nx * ny * nz; ++i) {
            data[i] = T();
        }
    }

    tensor3_t(size_t _nx, size_t _ny, size_t _nz, const T& elem = T())
        : nx(_nx), ny(_ny), nz(_nz) {
        data = new T[nx * ny * nz];
        for (size_t i = 0; i < nx * ny * nz; ++i) {
            data[i] = elem;
        }
    }

    tensor3_t(const tensor& t) : nx(t.nx), ny(t.ny), nz(t.nz) {
        delete[] data;
        data = new T[nx * ny * nz];
        for (size_t i = 0; i < nx * ny * nz; ++i) {
            data[i] = t.data[i];
        }
    }

    tensor& operator=(const tensor& t) {
        if (this == &t) {
            return *this;
        }
        nx = t.nx;
        ny = t.ny;
        nz = t.nz;
        delete[] data;
        data = new T[nx * ny * nz];
        for (size_t i = 0; i < nx * ny * nz; ++i) {
            data[i] = t.data[i];
        }
        return *this;
    }

    inline size_t ind(size_t i, size_t j, size_t k) const {
        assert(i < nx);
        assert(j < ny);
        assert(k < nz);
        return k * (nx * ny) + j + i * ny;
    }

    inline size_t ind_safe(size_t i, size_t j, size_t k) const {
        return k * (nx * ny) + j + i * ny;
    }

    T& operator()(size_t i, size_t j, size_t k) { return data[ind(i, j, k)]; }

    const T& operator()(size_t i, size_t j, size_t k) const {
        return data[ind(i, j, k)];
    }

    friend std::ostream& operator<<(std::ostream& out, const tensor& t) {
        for (size_t k = 0; k < t.nz; ++k) {
            if (k) {
                out << "\n\n";
            }
            for (size_t j = 0; j < t.ny; ++j) {
                if (j) {
                    out << '\n';
                }
                for (size_t i = 0; i < t.nx; ++i) {
                    if (i) {
                        out << ' ';
                    }
                    out << t.data[t.ind_safe(i, j, k)];
                }
            }
        }
        return out;
    }

    friend std::istream& operator>>(std::istream& in, tensor& t) {
        for (size_t k = 0; k < t.nz; ++k) {
            for (size_t j = 0; j < t.ny; ++j) {
                for (size_t i = 0; i < t.nx; ++i) {
                    in >> t.data[t.ind_safe(i, j, k)];
                }
            }
        }
        return in;
    }

    ~tensor3_t() { delete[] data; }
};

#endif /* TENSOR3_HPP */