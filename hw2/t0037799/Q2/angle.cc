#include <pybind11/pybind11.h>

#include <cmath>
#include <utility>

double angle(std::pair<double, double> x, std::pair<double, double> y) {
  auto normalize = [](std::pair<double, double> x) {
    double length = sqrt(x.first * x.first + x.second * x.second);
    x.first /= length;
    x.second /= length;
    return x;
  };
  x = normalize(x);
  y = normalize(y);
  double inner_product = x.first * y.first + x.second * y.second;
  return acos(inner_product);
}

PYBIND11_MODULE(angle, m) {
  m.doc() = "calculate angle (radian)";
  m.def("angle", &angle, "angle of two vector");
}
