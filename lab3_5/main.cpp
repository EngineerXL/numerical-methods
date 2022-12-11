#include <iostream>

#include "integrate.hpp"

using namespace std;

double f(double x) { return (x * x) / (16.0 + x * x); }

int main() {
    double l, r;
    cin >> l >> r;
    double h1, h2;
    cin >> h1 >> h2;

    double rect1 = integrate_rect(l, r, h1, f);
    cout << "Метод прямоугольников с шагом " << h1 << ": " << rect1 << endl;
    double trap1 = integrate_trap(l, r, h1, f);
    cout << "Метод трапеций с шагом " << h1 << ": " << trap1 << endl;
    double simp1 = integrate_simp(l, r, h1, f);
    cout << "Метод Симпсона с шагом " << h1 << ": " << simp1 << endl;

    cout << endl;

    double rect2 = integrate_rect(l, r, h2, f);
    cout << "Метод прямоугольников с шагом " << h2 << ": " << rect2 << endl;
    double trap2 = integrate_trap(l, r, h2, f);
    cout << "Метод трапеций с шагом " << h2 << ": " << trap2 << endl;
    double simp2 = integrate_simp(l, r, h2, f);
    cout << "Метод Симпсона с шагом " << h2 << ": " << simp2 << endl;

    cout << endl;

    double rect_rr = runge_romberg(rect1, rect2, h2 / h1, 2);
    cout << "Погрешность вычислений методом прямоугольников: " << rect_rr
         << endl;
    double trap_rr = runge_romberg(trap1, trap2, h2 / h1, 2);
    cout << "Погрешность вычислений методом трапеций: " << trap_rr << endl;
    double simp_rr = runge_romberg(simp1, simp2, h2 / h1, 2);
    cout << "Погрешность вычислений методом Симпсона: " << simp_rr << endl;
}
