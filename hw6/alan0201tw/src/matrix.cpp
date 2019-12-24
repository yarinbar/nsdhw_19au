#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "matrix.hpp"

#include <sstream>

#include <mkl.h>

namespace py = pybind11;

/*
 * Naive matrix matrix multiplication.
 */
Matrix operator*(Matrix const & mat1, Matrix const & mat2)
{
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range(
            "the number of first matrix column "
            "differs from that of second matrix row");
    }

    Matrix ret(mat1.nrow(), mat2.ncol());

    for (size_t i = 0; i < ret.nrow(); ++i)
    {
        for (size_t k = 0; k < ret.ncol(); ++k)
        {
            double v = 0;
            for (size_t j = 0; j < mat1.ncol(); ++j)
            {
                v += mat1(i,j) * mat2(j,k);
            }
            ret(i,k) = v;
        }
    }

    return ret;
}

Matrix multiply_naive(Matrix const & mat1, Matrix const & mat2)
{
    // redirect to operator *
    return mat1 * mat2;
}

// This example computes real matrix C=alpha*A*B+beta*C
//
// Reference : MKL Get started
Matrix multiply_mkl(Matrix const & mat1, Matrix const & mat2)
{
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range(
            "the number of first matrix column "
            "differs from that of second matrix row");
    }

    double* result = new double [mat1.nrow() * mat2.ncol()];

    cblas_dgemm(
        CblasRowMajor, CblasNoTrans, CblasNoTrans,
        mat1.nrow(), mat2.ncol(), mat1.ncol(),
        1.0, // alpha
        mat1.m_buffer, mat1.ncol(), mat2.m_buffer, mat2.ncol(),
        0.0, // beta
        result, mat2.ncol()
    );

    Matrix ret(mat1.nrow(), mat2.ncol(), result);

    return ret;
}

std::string Matrix::ToString()
{
    std::stringstream result("");
    for (size_t i=0; i<nrow(); ++i)
    {
        result << std::endl << " ";
        for (size_t j=0; j<ncol(); ++j)
        {
            result << " " << std::setw(2) << this->operator()(i, j);
        }
    }

    return result.str();
}

std::ostream & operator << (std::ostream & ostr, Matrix const & mat)
{
    for (size_t i=0; i<mat.nrow(); ++i)
    {
        ostr << std::endl << " ";
        for (size_t j=0; j<mat.ncol(); ++j)
        {
            ostr << " " << std::setw(2) << mat(i, j);
        }
    }

    return ostr;
}

//
// Reference :
// https://stackoverflow.com/questions/15829223/loop-tiling-blocking-for-large-dense-matrix-multiplication
//
Matrix multiply_tile(Matrix const & mat1, Matrix const & mat2, const unsigned int tsize)
{
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range(
            "the number of first matrix column "
            "differs from that of second matrix row");
    }

    Matrix ret(mat1.nrow(), mat2.ncol());

    const size_t tileSize = static_cast<size_t>(tsize);
    const size_t N = mat1.nrow();
    const size_t M = mat1.ncol();
    const size_t K = mat2.ncol();

    for(size_t i=0; i<N; i++)
    {
        for(size_t j=0; j<K; j++)
        {
            ret(i, j) = 0;
        }
    }

    // iterate by i,k,j instead of i,j,k gives another performance boost
    // i k j for row major, i j k for column major

    for (size_t i0 = 0; i0 < N; i0 += tileSize)
    {
        const size_t imax = i0 + tileSize > N ? N : i0 + tileSize;
        for (size_t k0 = 0; k0 < K; k0 += tileSize)
        {
            const size_t kmax = k0 + tileSize > K ? K : k0 + tileSize;
            for (size_t j0 = 0; j0 < M; j0 += tileSize)
            {
                const size_t jmax = j0 + tileSize > M ? M : j0 + tileSize;
                for (size_t i1 = i0; i1 < imax; ++i1)
                {
                    for (size_t k1 = k0; k1 < kmax; ++k1)
                    {
                        for (size_t j1 = j0; j1 < jmax; ++j1)
                        {
                            ret(i1, j1) += mat1(i1, k1) * mat2(k1, j1);
                        }
                    }
                }
            }
        }
    }

    return ret;
}


/**
 *  Reference :
 *  
 *  Basics : https://pybind11.readthedocs.io/en/stable/classes.html
 *  Operators : https://pybind11.readthedocs.io/en/stable/advanced/classes.html#operator-overloading
 *  Indexing operator : 
 *  https://pybind11.readthedocs.io/en/stable/reference.html#_CPPv4NK10object_apiixE6handle
 *  https://github.com/pybind/pybind11/blob/master/tests/test_sequences_and_iterators.cpp#L182
 */

PYBIND11_MODULE(_matrix, mod)
{
    mod.doc() = "Matrix class implementing naive multiplication and mkl delegation";
    mod.def("multiply_naive", &multiply_naive, "Naive implementation");
    mod.def("multiply_mkl", &multiply_mkl, "MKL implementation delegate");

    mod.def("multiply_tile", &multiply_tile, "Matrix Tiling implementation");

    // Reference : 
    // https://pybind11.readthedocs.io/en/stable/advanced/pycpp/numpy.html

    py::class_<Matrix>(mod, "Matrix", py::buffer_protocol())
        //
        .def_buffer([](Matrix &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),
                sizeof(double),
                py::format_descriptor<double>::format(),
                2,
                { m.nrow(), m.ncol() },
                { sizeof(double) * m.ncol(),
                  sizeof(double) }
            );
        })
        .def_property("array", &Matrix::numpy_data, nullptr)
        //
        .def( py::init<size_t, size_t>() )
        .def( py::init<size_t, size_t, std::vector<double>>() )
        .def_property_readonly("nrow", &Matrix::nrow)
        .def_property_readonly("ncol", &Matrix::ncol)
        .def("__str__", &Matrix::ToString )
        .def("__repr__", &Matrix::ToString )
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def("__getitem__",
            [](const Matrix &matrix, std::pair<size_t, size_t> index)
            {
                if (index.first >= matrix.nrow() || index.second >= matrix.ncol()) throw py::index_error();
                return matrix(index.first, index.second);
            }, py::is_operator())
        .def("__setitem__",
            [](Matrix &matrix, std::pair<size_t, size_t> index, double value)
            {
                if (index.first >= matrix.nrow() || index.second >= matrix.ncol()) throw py::index_error();
                matrix(index.first, index.second) = value;
            }, py::is_operator());
}