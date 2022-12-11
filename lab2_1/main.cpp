#include <iostream>

#include "solver.hpp"

using namespace std;

int main() {
    cout.precision(9);
    cout << fixed;
    double l, r, eps;
    cin >> l >> r >> eps;
    double root;
    root = iter_solve(l, r, eps);
    cout << "x_0 = " << root << endl;
    cout << "Решение методом простой итерации получено за " << iter_count
         << " итераций" << endl;
    root = newton_solve(l, r, eps);
    cout << "x_0 = " << root << endl;
    cout << "Решение методом Ньютона получена за " << iter_count << " итераций"
         << endl;
}
