#include "tridiag.hpp"

using namespace std;

using tridiag = tridiag_t<double>;
using vec = vector<double>;

int main() {
    cout.precision(6);
    cout << fixed;
    int n;
    cin >> n;
    tridiag tridiag_a(n);
    cin >> tridiag_a;
    vec b(n);
    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }
    vec x = tridiag_a.solve(b);
    cout << "Решение системы:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << x[i] << endl;
    }
}
