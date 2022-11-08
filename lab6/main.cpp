#include "hpde.hpp"

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    hpde_t hpdes;
    cin >> hpdes;
    cout.precision(12);
    cout << fixed;
    hpdes.solve(cout);
}
