#include <cmath>
#include <iostream>

#include "table_function.hpp"

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
    double x0;
    cin >> x0;

    cout.precision(4);
    cout << fixed;
    table_function_t f(x, y);
    cout << "Первая производная функции в точке x0 = " << x0
         << ", f'(x0) = " << f.derivative1(x0) << endl;
    cout << "Вторая производная функции в точке x0 = " << x0
         << ", f''(x0) = " << f.derivative2(x0) << endl;
}
