#include "boundary_solver.hpp"

using namespace std;

double g(double x, double y, double z) {
    return ((2 * x + 1) * z - (x + 1) * y) / x;
}

double f(double x, double y, double z) {
    (void)x;
    (void)y;
    return z;
}

double px(double x) {
    return -(2 * x + 1) / x;
}

double qx(double x) {
    return (x + 1) / x;
}

double fx(double x) {
    (void)x;
    return 0.0;
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
    double h, eps;
    cin >> h >> eps;
    /*
    Краевые условия 3 рода
    alpha * y(a) + beta * y'(a) = y0
    delta * y(b) + gamma * y'(b) = y1
    */
    double a = 1, b = 2;
    double alpha = 0, beta = 1, y0 = 3 * exp(1);
    double delta = -2, gamma = 1, y1 = 0;

    shooting de_shooting(a, b, f, g, alpha, beta, y0, delta, gamma, y1);
    vector<tddd> sol_shooting = de_shooting.solve(h, eps);
    cout << "Метод стрельбы:" << endl;
    print_data(sol_shooting);

    fin_dif de_fin_dif(a, b, px, qx, fx, alpha, beta, y0, delta, gamma, y1);
    vector<tddd> sol_fin_dif = de_fin_dif.solve(h);
    cout << "Конечно-разностный метод:" << endl;
    print_data(sol_fin_dif);
}
