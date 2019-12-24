#include "Vec2.cc"
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

namespace py = pybind11;

PYBIND11_MODULE(vec2, mod) {
  mod.doc() = "two dimension vector implementation";
  py::class_<Vec2>(mod, "vec2", py::dynamic_attr())
    .def(py::init<>())
    .def(py::init<double, double>())
    .def(py::init<Vec2 const &>())
    .def("copy", &Vec2::copy)
    .def("size", &Vec2::size)
    .def("length", &Vec2::size)
    .def("dot", &Vec2::dot)
    .def("angle", &Vec2::angle)
    .def("x", (double &(Vec2::*)(double)) &Vec2::x)
    .def("x", (double const &(Vec2::*)() const) &Vec2::x)
    //.def_readwrite("x", &Vec2::x_)
    .def("y", (double &(Vec2::*)(double)) &Vec2::y)
    .def("y", (double const &(Vec2::*)() const) &Vec2::y)
    //.def_readwrite("y", &Vec2::y_)
    .def("__repr__", &Vec2::repr)
    .def("__eq__", [](const Vec2 &self, const Vec2 &other){ return self == other; })
    //.def("__ne__", [](const Vec2 &self, const Vec2 &other){ return self != other; })
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self += py::self)
    .def(py::self *= double())
    .def(double() * py::self)
    .def(py::self * double())
    .def(+py::self)
    .def(-py::self);
}
