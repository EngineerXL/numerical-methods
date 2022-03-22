#include "qr_algo.hpp"

using namespace std;

using matrix = matrix_t<double>;
using complex_t = complex<double>;
using vec_complex = vector<complex_t>;

int main() {
    cout.precision(6);
    cout << fixed;
    int n;
    double eps;
    cin >> n >> eps;
    matrix a(n);
    cin >> a;
    qr_algo qr(a, eps);
    vec_complex lambda = qr.get_eigen_values();
    cout << "Собственные значения:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "l_" << i + 1 << " = " << lambda[i] << endl;
    }
}
