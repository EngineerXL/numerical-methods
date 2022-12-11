#include <cmath>
#include <iostream>

#include "cubic_spline.hpp"

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
    cubic_spline_t f(x, y);
    cout << "Полученные сплайны:\n" << f << endl;
    cout << "Значение функции в точке x0 = " << x0 << ", f(x0) = " << f(x0)
         << endl;
}
