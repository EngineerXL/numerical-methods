#include <iostream>

#include "fft_mult.hpp"

using namespace std;

int main() {
    string a, b;
    cin >> a >> b;
    string fft_res = fft_mult(a, b);
    cout << fft_res << endl;
}
