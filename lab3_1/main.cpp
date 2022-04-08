#include <cmath>
#include "interpolator.hpp"

using namespace std;

using vec = vector<double>;

int main() {
    int n;
    cin >> n;
    vec x(n), y(n);
    for (int i = 0; i < n; ++i) {
        cin >> x[i];
        y[i] = atan(x[i]);
    }
    double x_star;
    cin >> x_star;

    inter_lagrange my_lagrange(x, y);
    polynom lagrange = my_lagrange();
    cout << "Интерполяционный многочлен Лагранжа: " << lagrange << endl;
    cout << "Погрешность в точке X*: " << abs(lagrange(x_star) - atan(x_star)) << endl;

    inter_newton my_newton(x, y);
    polynom newton = my_newton();
    cout << "Интерполяционный многочлен Ньютона: " << newton << endl;
    cout << "Погрешность в точке X*: " << abs(newton(x_star) - atan(x_star)) << endl;
}
