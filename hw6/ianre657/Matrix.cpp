#include <stdexcept>
#include <algorithm>
#include <cstring>
#include "Matrix.hpp"

namespace py = pybind11;

py::array_t<double> Matrix::array()
{
    auto capsule = py::capsule(m_buffer, [](void *p) { /*delete reinterpret_cast<double *>(p);*/ });
    return py::array_t<double>({m_nrow, m_ncol}, m_buffer, capsule);
}

Matrix multiply_naive(Matrix const &mat1, Matrix const &mat2)
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
            for (size_t j = 0; j < mat1.ncol(); ++j)
            {
                ret(i, k) += mat1(i, j) * mat2(j, k);
            }
        }
    }

    return ret;
}

Matrix &multiply_mkl(Matrix const &mat1, Matrix const &mat2)
{
    double alpha = 1.0;
    double beta = 0.0;

    Matrix *mat = new Matrix(mat1.nrow(), mat2.ncol());
    cblas_dgemm(
        CblasRowMajor,
        CblasNoTrans, CblasNoTrans,
        mat->nrow(), mat->ncol(), mat1.ncol(),
        alpha,
        mat1.raw_data(), mat1.ncol(),
        mat2.raw_data(), mat2.ncol(),
        beta,
        mat->raw_data(), mat->ncol());
    return *mat;
}

Matrix multiply_tile(Matrix const &mat1, Matrix const &mat2, size_t t_size)
{
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range(
            "the number of first matrix column "
            "differs from that of second matrix row");
    }

    size_t M = mat1.nrow();
    size_t N = mat1.ncol();
    size_t K = mat2.ncol();

    Matrix ret(M, K);
    std::fill(ret.raw_data(), ret.raw_data() + M * K, 0.0);

    for (size_t i = 0; i < M; i += t_size)
    {
        for (size_t k = 0; k < K; k += t_size)
        {
            for (size_t j = 0; j < N; j += t_size)
            {
                size_t it_max = std::min(i + t_size, M);
                for (size_t it = i; it < it_max; ++it)
                {
                    size_t kt_max = std::min(k + t_size, K);
                    for (size_t kt = k; kt < kt_max; ++kt)
                    {
                        size_t jt_max = std::min(j + t_size, N);
                        for (size_t jt = j; jt < jt_max; ++jt)
                        {
                            ret(it, kt) += mat1(it, jt) * mat2(jt, kt);
                        }
                    }
                }
            }
        }
    }
    return ret;
}