#include "iteration.hpp"

using namespace std;

using matrix = matrix_t<double>;
using vec = vector<double>;

int main() {
    cout.precision(6);
    cout << fixed;
    int n;
    double eps;
    cin >> n >> eps;
    matrix a(n);
    cin >> a;
    iter_solver my_solver(a, eps);
    vec b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }
    vec x = my_solver.solve_simple(b);
    cout << "Метод простых итераций" << endl;
    cout << "Решени получено за " << my_solver.iter_count << " итераций"
         << endl;
    cout << "Решение системы:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
    vec z = my_solver.solve_zeidel(b);
    cout << "Метод Зейделя" << endl;
    cout << "Решени получено за " << my_solver.iter_count << " итераций"
         << endl;
    cout << "Решение системы:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << z[i] << endl;
    }
}
