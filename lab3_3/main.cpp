#include <cmath>
#include <iostream>
#include "minimal_square.hpp"

using namespace std;

using vec = vector<double>;

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
    minimal_square_t ms1(x, y, 1);
    cout << "Полученная функция первого порядка: " << ms1 << endl;
    cout << "Значение суммы квадратов ошибков: " << ms1.mse() << endl;

    minimal_square_t ms2(x, y, 2);
    cout << "Полученная функция второго порядка: " << ms2 << endl;
    cout << "Значение суммы квадратов ошибков: " << ms2.mse() << endl;
}
