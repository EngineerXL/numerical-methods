#include "rotation.hpp"

using namespace std;

using matrix = matrix_t<double>;

int main() {
    cout.precision(6);
    cout << fixed;
    int n;
    double eps;
    cin >> n >> eps;
    matrix a(n);
    cin >> a;
    rotation rot(a, eps);
    vector<double> lambda = rot.get_eigen_values();
    cout << "Собственные значения:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "l_" << i + 1 << " = " << lambda[i] << endl;
    }
    matrix s = rot.get_eigen_vectors();
    cout << "Собственные векторы:" << endl << s;
    cout << "Решение получено за " << rot.iter_count << " итераций" << endl;
}
