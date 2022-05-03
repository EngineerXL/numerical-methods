#include "fft_mult.hpp"
#include <iostream>

using namespace std;

int main() {
    string a, b;
    cin >> a >> b;
    std::cout.precision(3);
    int64_t start, end;
    double time;

    start = clock();
    string fft_res = fft_mult(a, b);
    end = clock();
    time = (end - start) / 1000.0;
    std::cout << std::fixed << "FFT time is " << time << " ms" << endl;

    start = clock();
    string slow_res = slow_mult(a, b);
    end = clock();
    time = (end - start) / 1000.0;
    std::cout << std::fixed << "Naive time is " << time << " ms" << endl;
}
