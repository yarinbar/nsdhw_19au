#include <utility>
#include <vector>
using namespace std;

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;

#include <mkl.h>

class Matrix {
   public:
    Matrix(size_t n, size_t m) : n(n), m(m) {
        a = vector<double>(n * m, 0);
    }
    size_t nrow() const {
        return n;
    }
    size_t ncol() const {
        return m;
    }
    double *operator[](size_t idx) {
        return a.data() + idx * m;
    }
    const double *operator[](size_t idx) const {
        return a.data() + idx * m;
    }
    py::array_t<double> array() {
        auto capsule = py::capsule(a.data(), [](void *p) {});
        return py::array_t<double>({n, m}, a.data(), capsule);
    }

   private:
    size_t n;
    size_t m;
    vector<double> a;
};

bool isequal(const Matrix &A, const Matrix &B) {
    if (A.nrow() != B.nrow()) return 0;
    if (A.ncol() != B.ncol()) return 0;
    for (size_t i = 0; i < A.nrow(); i++) {
        for (size_t j = 0; j < A.ncol(); j++) {
            if (A[i][j] != B[i][j]) {
                return 0;
            }
        }
    }
    return 1;
}

Matrix multiply_naive(const Matrix &A, const Matrix &B) {
    if (A.ncol() != B.nrow()) throw out_of_range("size error");
    Matrix C(A.nrow(), B.ncol());
    for (size_t i = 0; i < A.nrow(); i++) {
        for (size_t j = 0; j < B.nrow(); j++) {
            for (size_t k = 0; k < B.ncol(); k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

vector<double> to1D(const Matrix &A) {
    vector<double> a(A.nrow() * A.ncol());
    for (size_t i = 0, ptr = 0; i < A.nrow(); i++) {
        for (size_t j = 0; j < A.ncol(); j++) {
            a[ptr++] = A[i][j];
        }
    }
    return a;
}

Matrix multiply_mkl(const Matrix &A, const Matrix &B) {
    if (A.ncol() != B.nrow()) throw out_of_range("size error");
    Matrix C(A.nrow(), B.ncol());
    vector<double> a = to1D(A);
    vector<double> b = to1D(B);
    vector<double> c(A.nrow() * B.ncol());
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A.nrow(), B.ncol(),
                A.ncol(), 1, a.data(), A.ncol(), b.data(), B.ncol(), 0,
                c.data(), C.ncol());
    for (size_t i = 0, ptr = 0; i < A.nrow(); i++) {
        for (size_t j = 0; j < B.ncol(); j++) {
            C[i][j] = c[ptr++];
        }
    }
    return C;
}

PYBIND11_MODULE(_matrix, m) {
    py::class_<Matrix>(m, "Matrix")
        .def(py::init<size_t, size_t>())
        .def_property_readonly("nrow", &Matrix::nrow)
        .def_property_readonly("ncol", &Matrix::ncol)
        .def("__eq__", &isequal)
        .def("__getitem__",
             [](const Matrix &A, const pair<size_t, size_t> idx) {
                 return A[idx.first][idx.second];
             })
        .def("__setitem__",
             [](Matrix &A, const pair<size_t, size_t> idx, double v) {
                 return A[idx.first][idx.second] = v;
             })
        .def_property_readonly("array", &Matrix::array);
    m.def("isequal", &isequal);
    m.def("multiply_naive", &multiply_naive);
    m.def("multiply_mkl", &multiply_mkl);
}
