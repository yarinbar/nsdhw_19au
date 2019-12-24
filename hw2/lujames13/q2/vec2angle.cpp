#include <pybind11/pybind11.h>
#include <cmath>
namespace py = pybind11;

double vec2angle(std::pair<double, double> v1, std::pair<double, double> v2){
	double norm = sqrt(pow(v1.first, 2) + pow(v1.second, 2)) * sqrt(pow(v2.first, 2) + pow(v2.second, 2));
	double dot = v1.first * v2.first + v1.second * v2.second;
	double angle = acos(dot/norm);
	return angle;
}

PYBIND11_MODULE(vec2angle, m){
	m.doc() = "pybin11 plugin to calculate the angle between 2 vectors";
	m.def("vec2angle", &vec2angle, "A function which calculate the angle between 2 vectors", py::arg("v1"), py::arg("v2"));

}