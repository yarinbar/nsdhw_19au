#include "matrix.hpp"
#include <utility>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_matrix, mod) {
  mod.doc() = "matrix implementation";
  mod.def("multiply_naive", &multiply_naive, "naive version of matrix multiplication");
  mod.def("multiply_tile", &multiply_tile, "tiled version of matrix multiplication");
  mod.def("multiply_mkl", &multiply_mkl, "cblas version of matrix multiplication");
  py::class_<Matrix>(mod, "Matrix", py::dynamic_attr())
    .def(py::init<size_t, size_t>())
    .def(py::init<size_t, size_t, std::vector<double>>())
    .def(py::self * py::self)
    .def_property_readonly("nrow", &Matrix::nrow)
    .def_property_readonly("ncol", &Matrix::ncol)
    .def("__repr__", &Matrix::show)
    .def_property("array", &Matrix::ndarray, nullptr)
    .def("__eq__", [](const Matrix &self, const Matrix &other)
        { return self == other; })
    .def("__getitem__", [](const Matrix &self, std::pair<size_t, size_t>idx)
        { return self(idx.first, idx.second); })
    .def("__setitem__", [](Matrix &self, std::pair<size_t, size_t>idx, double val)
        { self(idx.first, idx.second) = val; });
}
