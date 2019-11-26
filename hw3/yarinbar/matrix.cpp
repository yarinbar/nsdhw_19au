
#include <iostream>
#include <iomanip>
#include <mkl.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

class Matrix {

public:

    Matrix(size_t nrow, size_t ncol)
            : m_nrow(nrow), m_ncol(ncol)
    {
        size_t nelement = nrow * ncol;
        m_buffer = new double[nelement];
    }


    ~Matrix()
    {
        delete[] m_buffer;
    }

    // No bound check.
    double   operator() (size_t row, size_t col) const { return m_buffer[row*m_ncol + col]; }
    double & operator() (size_t row, size_t col)       { return m_buffer[row*m_ncol + col]; }

    size_t nrow() const { return m_nrow; }
    size_t ncol() const { return m_ncol; }


    size_t m_nrow;
    size_t m_ncol;
    double * m_buffer = nullptr;

};

void work(Matrix & matrix)
{
    for (size_t i=0; i<matrix.nrow(); ++i) // the i-th row
    {
        for (size_t j=0; j<matrix.ncol(); ++j) // the j-th column
        {
            matrix(i, j) = i*10 + j;
        }
    }
}

Matrix multiply_naive(Matrix const& A, Matrix const& B){

    // check for dim compatibility
    if(A.ncol() != B.nrow())
        throw "sizes dont match!";

    Matrix C(A.nrow(), B.ncol());

    for(unsigned int i = 0; i < A.nrow(); ++i)
        for(unsigned int j = 0; j < B.ncol(); ++j) {
            C(i, j) = 0;
            for(unsigned int k = 0; k < A.ncol(); k++)
                C(i, j) += A(i, k) * B(k, j);
        }
    return C;
}

Matrix multiply_mkl(const Matrix& A, const Matrix& B){

	// check for dim compatibility
    if(A.ncol() != B.nrow()){
        throw "sizes dont match!";
    }

	Matrix C(a.nrow(), b.ncol());

    mkl_set_num_threads(1);
	
    cblas_dgemm(
        CblasRowMajor /* const CBLAS_LAYOUT Layout */
      , CblasNoTrans /* const CBLAS_TRANSPOSE transa */
      , CblasNoTrans /* const CBLAS_TRANSPOSE transb */
      , A.nrow() /* const MKL_INT m */
      , B.ncol() /* const MKL_INT n */
      , A.ncol() /* const MKL_INT k */
      , 1.0 /* const double alpha */
      , A.m_buffer /* const double *a */
      , A.ncol() /* const MKL_INT lda */
      , B.m_buffer /* const double *b */
      , B.ncol() /* const MKL_INT ldb */
      , 0.0 /* const double beta */
      , C.m_buffer /* double * c */
      , C.ncol() /* const MKL_INT ldc */
    );
	
	return C;
}

PYBIND11_MODULE(_matrix, m) {
    py::class_<Matrix>(m, "Matrix")
      .def(py::init<size_t, size_t>())
      .def_property_readonly("nrow", &Matrix::nrow)
      .def_property_readonly("ncol", &Matrix::ncol)
      .def("__eq__", [](Matrix &a, Matrix &b) { return a == b; })
      .def("__getitem__",
           [](Matrix &m, std::pair<size_t, size_t> i) {
             return m(i.first, i.second);
           })
      .def("__setitem__", [](Matrix &m, std::pair<size_t, size_t> i,
                             double v) { m(i.first, i.second) = v; });

  m.def("multiply_naive", &multiply_naive, "naive");
  m.def("multiply_mkl", &multiply_mkl, "mkl");
}

// vim: set ff=unix fenc=utf8 et sw=4 ts=4 sts=4:

