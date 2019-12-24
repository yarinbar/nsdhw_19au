#include <pybind11/pybind11.h>

#include <cmath>

#define PI (3.14159265358979323846)

namespace py = pybind11;

class Point
{
public:
    Point(float p_x, float p_y): x(p_x), y(p_y)
    {
    }
    float x, y;
};

float calc_angle(const Point &p1, const Point &p2)
{
    if(p1.x==p2.x && p1.y==p2.y)
        return -1.0;

    float len1 = sqrt(p1.x * p1.x + p1.y * p1.y);
    float len2 = sqrt(p2.x * p2.x + p2.y * p2.y);

    float dot = p1.x * p2.x + p1.y * p2.y;

    float a = dot / (len1 * len2);

    if(a >= 1.0)
        return 0.0;
    else if(a <= -1.0)
        return PI;
    else
        return acos(a);
}

PYBIND11_MODULE(angle, m) {
  m.doc() = "pybind11 angle calculator";      // module doc string
  py::class_<Point>(m, "Point")
    .def(py::init<float, float>());
  m.def("calc_angle",                              // function name
        &calc_angle,                               // function pointer
        "A function which calculate angle between two points." //function doc string
       );
}