#include "fft_mult.hpp"
#include <iostream>

using namespace std;

int main() {
    string a, b;
    cin >> a >> b;
    string fft_res = fft_mult(a, b);
    cout << fft_res << endl;
}
