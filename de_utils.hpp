#ifndef DE_UTILS_HPP
#define DE_UTILS_HPP

#include <iostream>
#include <vector>
#include <tuple>

using tddd = std::tuple<double, double, double>;

void print_data(const std::vector<tddd> & v) {
    std::cout << "x = [";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) {
            std::cout << ", ";
        }
        std::cout << std::get<0>(v[i]);
    }
    std::cout << "]\n";
    std::cout << "y = [";
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) {
            std::cout << ", ";
        }
        std::cout << std::get<1>(v[i]);
    }
    std::cout << "]\n";
}

void print_err(const std::vector<double> & v) {
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) {
            std::cout << ", ";
        }
        std::cout << v[i];
    }
    std::cout << "\n";
}

#endif /* DE_UTILS_HPP */
