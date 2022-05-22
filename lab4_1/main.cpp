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

void print_data(const vector<tddd> & v) {
    cout << "x = [";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) {
            cout << ", ";
        }
        cout << get<0>(v[i]);
    }
    cout << "]\n";
    cout << "y = [";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) {
            cout << ", ";
        }
        cout << get<1>(v[i]);
    }
    cout << "]\n";
}

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

    runge de_runge(l, r, f, g, y0, z0);
    vector<tddd> sol_runge = de_runge.solve(h);
    cout << "Метод Рунге-Кутта:" << endl;
    print_data(sol_runge);

    adams de_adams(l, r, f, g, y0, z0);
    vector<tddd> sol_adams = de_adams.solve(h);
    cout << "Метод Адамса:" << endl;
    print_data(sol_adams);

}
