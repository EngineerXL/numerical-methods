#ifndef TABLE_FUNCTION_HPP
#define TABLE_FUNCTION_HPP

#include <vector>
#include <exception>

class table_function_t {
    using vec = std::vector<double>;
    size_t n;
    vec x;
    vec y;

public:
    table_function_t(const vec & _x, const vec & _y) {
        if (_x.size() != _y.size()) {
            throw std::invalid_argument("Sizes does not match");
        }
        x = _x;
        y = _y;
        n = x.size();
    }

    double derivative1(double x0) {
        for (size_t i = 0; i < n - 2; ++i) {
            if (x[i] <= x0 and x0 <= x[i + 1]) {
                double dydx1 = (y[i + 1] - y[i + 0]) / (x[i + 1] - x[i + 0]);
                double dydx2 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
                double res = dydx1 + (dydx2 - dydx1) * (2.0 * x0 - x[i] - x[i + 1]) / (x[i + 2] - x[i]);
                return res;
            }
        }
        return NAN;
    }

    double derivative2(double x0) {
        for (size_t i = 0; i < n - 2; ++i) {
            if (x[i] <= x0 and x0 <= x[i + 1]) {
                double dydx1 = (y[i + 1] - y[i + 0]) / (x[i + 1] - x[i + 0]);
                double dydx2 = (y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]);
                double res = 2.0 * (dydx2 - dydx1) / (x[i + 2] - x[i]);
                return res;
            }
        }
        return NAN;
    }
};

#endif /* TABLE_FUNCTION_HPP */
