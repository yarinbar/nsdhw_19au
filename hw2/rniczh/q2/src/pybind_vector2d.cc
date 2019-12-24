#include "vector2d.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

namespace py = pybind11;

// Wrapping the C++ function for python
namespace {
  void set_x(Vector2D &lhs, float x) { lhs.x() = x; }
  void set_y(Vector2D &lhs, float y) { lhs.y() = y; }

  float get_x(Vector2D &v) { return v.x(); }
  float get_y(Vector2D &v) { return v.y(); }
}

// Interface for the python
PYBIND11_MODULE(pybind_vector2d, m) {
  m.doc() = "vector2d plugin";

  m.def("angle", &Vector2D::angle, "Compute the angle of a vector (in radians)");

  py::class_<Vector2D>(m, "Vector2D")
    .def(py::init<float, float>()) // usage: v1 = vector2d.Vector2D(1,2)
    .def("setX", &set_x, "Set the x of a vector with argument: float")
    .def("setY", &set_y, "Set the y of a vector with argument: float")
    .def("getX", &get_x, "Get the x of a vector")
    .def("getY", &get_y, "Get the y of a vector")
    .def("norm", &Vector2D::norm, "Compute the normalize of a vector")
    .def(py::self + py::self, "Add of two vector")
    .def(py::self - py::self, "Sub of two vector")
    .def(py::self * py::self, "Dot product of two vectr")
    .def(py::self * float(), "Times the vector with scalar")
    .def(float() * py::self, "Times the vector with scalar")
    .def(py::self *= float(), "Lengthen the vector with scalar")
    .def(py::self += py::self, "Add to the lhs of the vector")
    .def(py::self -= py::self, "Sub to the lhs of the vector")
    .def(-py::self, "Negtaive of the vector")
    .def(+py::self, "Positive of the vector")
    .def("__repr__", &Vector2D::show);
}


