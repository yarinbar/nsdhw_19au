#include "Matrix.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

PYBIND11_MODULE(_matrix, m)
{
      m.doc() = "Matrix library";
      m.def("multiply_naive", &multiply_naive,
            "use naive method to multiply matrices.");
      m.def("multiply_mkl", &multiply_mkl,
            "use mkl to multiply matrices.");
      m.def("multiply_tile", &multiply_tile,
            "multiply matrices with tile optimization");

      py::class_<Matrix>(m, "Matrix")
          .def(py::init<int, int>())
          .def_property_readonly("nrow", &Matrix::nrow)
          .def_property_readonly("ncol", &Matrix::ncol)
          .def_property("array", &Matrix::array, nullptr)
          .def("__eq__", [](const Matrix &self, const Matrix &other) { return self == other; })
          .def("__setitem__", [](Matrix &self, std::pair<size_t, size_t> idx, double val) { return self(idx.first, idx.second) = val; })
          .def("__getitem__", [](const Matrix &self, std::pair<size_t, size_t> idx) { return self(idx.first, idx.second); });
}