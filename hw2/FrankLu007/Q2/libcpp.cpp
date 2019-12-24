#include <pybind11/pybind11.h>
#include <cmath>
#include <iostream>

typedef std::pair<float, float> point;

float dot(point &n1, point &n2) {
    return n1.first * n2.first + n1.second * n2.second;
}

float angle(point n1, point n2) {
    return acos(dot(n1, n2) / sqrt(dot(n1, n1) * dot(n2, n2)));
}

PYBIND11_MODULE(libcpp, n) {
    n.def("angle",                    // function name
          &angle,                     // function pointer
          "Calculate the angle between two vectors in 2D system."  // function doc string
    );
}