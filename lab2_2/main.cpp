#include <iostream>

#include "system_solver.hpp"

using namespace std;

using pdd = pair<double, double>;

int main() {
    cout.precision(9);
    cout << fixed;
    double l1, r1, l2, r2, eps;
    cin >> l1 >> r1 >> l2 >> r2 >> eps;
    auto [x0, y0] = iter_solve(l1, r1, l2, r2, eps);
    cout << "x_0 = " << x0 << ", y0 = " << y0 << endl;
    cout << "Решение методом простой итерации получено за " << iter_count
         << " итераций" << endl;
    auto [x0_n, y0_n] = newton_solve(r1, r2, eps);
    cout << "x_0 = " << x0_n << ", y0 = " << y0_n << endl;
    cout << "Решение методом Ньютона получена за " << iter_count << " итераций"
         << endl;
}
