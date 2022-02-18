#include "lu.hpp"

using namespace std;

using matrix = matrix_t<double>;
using vec = vector<double>;
using lu = lu_t<double>;

matrix obr(lu lu_a, size_t n) {
    matrix res(n);
    for (size_t i = 0; i < n; ++i) {
        vec b(n);
        b[i] = 1;
        vec x = lu_a.solve(b);
        for (size_t j = 0; j < n; ++j) {
            res[j][i] = x[j];
        }
    }
    return res;
}

int main() {
    cout.precision(12);
    cout << fixed;
    int n;
    cin >> n;
    matrix a(n);
    cin >> a;
    lu lu_a(a);
    cout << "Обратная матрица:" << endl;
    cout << obr(lu_a, n) << endl;
    vec b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }
    vec x = lu_a.solve(b);
    cout << "Решение системы:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}
