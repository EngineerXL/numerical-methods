#include <cmath>
#include <iostream>
#include "simple_desolve.hpp"

using namespace std;

double g(double x, double y, double z) {
    return -y * pow(cos(x), 2) - z * tan(x);
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
    vector<double> euler_err = runge_romberg(de_euler.solve(h), de_euler.solve(h / 2), 1);
    print_err(euler_err);

    runge de_runge(l, r, f, g, y0, z0);
    vector<tddd> sol_runge = de_runge.solve(h);
    cout << "Метод Рунге-Кутта:" << endl;
    print_data(sol_runge);
    cout << "Погрешность вычислений:" << endl;
    vector<double> runge_err = runge_romberg(de_runge.solve(h), de_runge.solve(h / 2), 4);
    print_err(runge_err);

    adams de_adams(l, r, f, g, y0, z0);
    vector<tddd> sol_adams = de_adams.solve(h);
    cout << "Метод Адамса:" << endl;
    print_data(sol_adams);
    cout << "Погрешность вычислений:" << endl;
    vector<double> adams_err = runge_romberg(de_adams.solve(h), de_adams.solve(h / 2), 4);
    print_err(adams_err);

}
