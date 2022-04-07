#include "fft_mult.hpp"
#include <iostream>

using namespace std;

int main() {
    string a, b;
    cin >> a >> b;
    int64_t start, end;
    
    start = clock();
    string fft_res = fft_mult(a, b);
    end = clock();
    std::cout << "FFT time is " << (end - start) / 1000 << " ms" << endl;
    // cout << fft_res << std::endl;

    start = clock();
    string slow_res = slow_mult(a, b);
    end = clock();
    std::cout << "Naive time is " << (end - start) / 1000 << " ms" << endl;
    // cout << slow_res << std::endl;
}
