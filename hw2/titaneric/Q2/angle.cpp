#include <pybind11/pybind11.h>

#include <cmath>
#include <utility>

namespace py = pybind11;
typedef std::pair<double, double> Vector;
double dot(Vector vec1, Vector vec2)
{
    return vec1.first * vec2.first + vec1.second * vec2.second;
}
double getAngle(Vector vec1, Vector vec2)
{
    double dot_result = dot(vec1, vec2);
    double vec1_len = sqrt(dot(vec1, vec1));
    double vec2_len = sqrt(dot(vec2, vec2));
    return (double)acos(dot_result / (vec1_len * vec2_len));
}

PYBIND11_MODULE(angle, m)
{
    m.doc() = "pybind11 angle plugin"; // optional module docstring
    m.def("getAngle", &getAngle, "A function which calculates angle between two vectors");
}