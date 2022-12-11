#include <cmath>
#include <iostream>

#include "minimal_square.hpp"

using namespace std;

using vec = vector<double>;
using func = std::function<double(double)>;
using vf = vector<func>;

double f0(double x0) {
    (void)x0;
    return 1.0;
}

double f1(double x0) { return x0; }

double f2(double x0) { return x0 * x0; }

int main() {
    int n;
    cin >> n;
    vec x(n), y(n);
    for (int i = 0; i < n; ++i) {
        cin >> x[i];
    }
    for (int i = 0; i < n; ++i) {
        cin >> y[i];
    }

    cout.precision(4);
    cout << fixed;
    vf phi1 = {f0, f1};
    minimal_square_t ms1(x, y, phi1);
    cout << "Полученная функция первого порядка: " << ms1 << endl;
    cout << "Значение суммы квадратов ошибков: " << ms1.mse() << endl;

    vf phi2 = {f0, f1, f2};
    minimal_square_t ms2(x, y, phi2);
    cout << "Полученная функция второго порядка: " << ms2 << endl;
    cout << "Значение суммы квадратов ошибков: " << ms2.mse() << endl;
}
