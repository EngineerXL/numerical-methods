#include "ppde.hpp"

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    ppde_t ppdes;
    cin >> ppdes;
    cout.precision(12);
    cout << fixed;
    ppdes.solve(cout);
}
