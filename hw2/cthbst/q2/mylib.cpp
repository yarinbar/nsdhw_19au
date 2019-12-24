#include <pybind11/pybind11.h>
#include <cmath>
#include <iostream>

using Point = std::pair<float, float>;

float dot(const Point &a, const Point &b) {
    return a.first * b.first + a.second * b.second;
}

float angle(Point a, Point b) {
    return acos(dot(a, b) / sqrt(dot(a, a) * dot(b, b)));
}

PYBIND11_MODULE(mylib, m) {
    m.def("angle",                    // function name
          &angle,                     // function pointer
          "Get angle of two vector."  // function doc string
    );
}
