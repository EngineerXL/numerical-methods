#include "epde.hpp"

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    epde_t epdes;
    cin >> epdes;
    cout.precision(12);
    cout << fixed;
    epdes.solve(cout);
}
