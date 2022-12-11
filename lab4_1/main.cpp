#include <cmath>
#include <iostream>

#include "simple_desolve.hpp"

using namespace std;

double g(double x, double y, double z) {
    return (z * x * (2 * x + 1) - y * (2 * x + 1)) / (x * x * (x + 1));
}

double f(double x, double y, double z) {
    (void)x;
    (void)y;
    return z;
}

using tddd = tuple<double, double, double>;

int main() {
    cout.precision(6);
    cout << fixed;
    double l, r, y0, z0, h;
    cin >> l >> r;
    cin >> y0 >> z0 >> h;

    runge de_euler(l, r, f, g, y0, z0);
    vector<tddd> sol_euler = de_euler.solve(h);
    cout << "Метод Эйлера:" << endl;
    print_data(sol_euler);
    cout << "Погрешность вычислений:" << endl;
    double euler_err =
        runge_romberg(de_euler.solve(h), de_euler.solve(h / 2), 1);
    cout << euler_err << endl;

    runge de_runge(l, r, f, g, y0, z0);
    vector<tddd> sol_runge = de_runge.solve(h);
    cout << "Метод Рунге-Кутты:" << endl;
    print_data(sol_runge);
    cout << "Погрешность вычислений:" << endl;
    double runge_err =
        runge_romberg(de_runge.solve(h), de_runge.solve(h / 2), 4);
    cout << runge_err << endl;

    adams de_adams(l, r, f, g, y0, z0);
    vector<tddd> sol_adams = de_adams.solve(h);
    cout << "Метод Адамса:" << endl;
    print_data(sol_adams);
    cout << "Погрешность вычислений:" << endl;
    double adams_err =
        runge_romberg(de_adams.solve(h), de_adams.solve(h / 2), 4);
    cout << adams_err << endl;
}
