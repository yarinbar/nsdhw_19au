#include "StopWatch.hpp"

#include <mkl.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <functional>

struct Matrix {

public:
	
	Matrix() = delete;
	
    Matrix(size_t nrow, size_t ncol)
      : m_nrow(nrow), m_ncol(ncol)
    {
        reset_buffer(nrow, ncol);
    }

    Matrix(size_t nrow, size_t ncol, std::vector<double> const & vec)
      : m_nrow(nrow), m_ncol(ncol)
    {
        reset_buffer(nrow, ncol);
        (*this) = vec;
    }

    Matrix & operator=(std::vector<double> const & vec){
        if (size() != vec.size())
        {
            throw std::out_of_range("number of elements mismatch");
        }

        size_t k = 0;
        for (size_t i=0; i<m_nrow; ++i)
        {
            for (size_t j=0; j<m_ncol; ++j)
            {
                (*this)(i,j) = vec[k];
                ++k;
            }
        }

        return *this;
    }

    Matrix(Matrix const & other)
      : m_nrow(other.m_nrow), m_ncol(other.m_ncol)
      , m_elapsed(other.m_elapsed), m_nflo(other.m_nflo){
        reset_buffer(other.m_nrow, other.m_ncol);
        for (size_t i=0; i<m_nrow; ++i)
        {
            for (size_t j=0; j<m_ncol; ++j)
            {
                (*this)(i,j) = other(i,j);
            }
        }
    }

    Matrix & operator=(Matrix const & other){
        if (this == &other) { return *this; }
        if (m_nrow != other.m_nrow || m_ncol != other.m_ncol)
        {
            reset_buffer(other.m_nrow, other.m_ncol);
        }
        for (size_t i=0; i<m_nrow; ++i)
        {
            for (size_t j=0; j<m_ncol; ++j)
            {
                (*this)(i,j) = other(i,j);
            }
        }
        m_elapsed = other.m_elapsed;
        m_nflo = other.m_nflo;
        return *this;
    }

    Matrix(Matrix && other)
      : m_nrow(other.m_nrow), m_ncol(other.m_ncol)
      , m_elapsed(other.m_elapsed), m_nflo(other.m_nflo){
        reset_buffer(0, 0);
        std::swap(m_nrow, other.m_nrow);
        std::swap(m_ncol, other.m_ncol);
        std::swap(m_buffer, other.m_buffer);
    }

    Matrix & operator=(Matrix && other){
        if (this == &other) { return *this; }
        reset_buffer(0, 0);
        std::swap(m_nrow, other.m_nrow);
        std::swap(m_ncol, other.m_ncol);
        std::swap(m_buffer, other.m_buffer);
        std::swap(m_elapsed, other.m_elapsed);
        std::swap(m_nflo, other.m_nflo);
        return *this;
    }

    ~Matrix(){
        reset_buffer(0, 0);
    }

    double   operator() (size_t row, size_t col) const { return m_buffer[index(row, col)]; }
    double & operator() (size_t row, size_t col)       { return m_buffer[index(row, col)]; }

    double   operator[] (size_t idx) const { return m_buffer[idx]; }
    double & operator[] (size_t idx)       { return m_buffer[idx]; }

    size_t nrow() const { return m_nrow; }
    size_t ncol() const { return m_ncol; }

    size_t size() const { return m_nrow * m_ncol; }
    double buffer(size_t i) const { return m_buffer[i]; }
    std::vector<double> buffer_vector() const { return std::vector<double>(m_buffer, m_buffer+size()); }

    double   elapsed() const { return m_elapsed; }
    double & elapsed()       { return m_elapsed; }

    size_t   nflo() const { return m_nflo; }
    size_t & nflo()       { return m_nflo; }

    double gflops() const { return m_nflo / m_elapsed / 1.e9; }

    Matrix transpose() const;

public:

    size_t index(size_t row, size_t col) const
    {
        return row * m_ncol + col;
    }

    void reset_buffer(size_t nrow, size_t ncol)
    {
        if (m_buffer) { delete[] m_buffer; }
        const size_t nelement = nrow * ncol;
        if (nelement) { m_buffer = new double[nelement]; }
        else          { m_buffer = nullptr; }
        m_nrow = nrow;
        m_ncol = ncol;
    }

    size_t m_nrow = 0;
    size_t m_ncol = 0;
    double * m_buffer = nullptr;
    double m_elapsed = 0;
    size_t m_nflo = 0; // number of floating-point operations.

};

Matrix Matrix::transpose() const{
    Matrix ret(nrow(), ncol());

    for (size_t i=0; i<ret.nrow(); ++i)
    {
        for (size_t j=0; j<ret.ncol(); ++j)
        {
            ret(j, i) = (*this)(i, j);
        }
    }

    return ret;
}

bool operator== (Matrix const & mat1, Matrix const & mat2){
    if ((mat1.ncol() != mat2.ncol()) && (mat1.nrow() != mat2.ncol()))
    {
        return false;
    }

    for (size_t i=0; i<mat1.nrow(); ++i)
    {
        for (size_t j=0; j<mat1.ncol(); ++j)
        {
            if (mat1(i, j) != mat2(i, j))
            {
                return false;
            }
        }
    }

    return true;
}

bool operator!= (Matrix const & mat1, Matrix const & mat2){
    return !(mat1 == mat2);
}

/*
 * Throw an exception if the shapes of the two matrices don't support
 * multiplication.
 */
void validate_multiplication(Matrix const & mat1, Matrix const & mat2){
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range(
            "the number of first matrix column "
            "differs from that of second matrix row");
    }
}

/*
 * Get the number of floating-point operations.
 */
size_t calc_nflo(Matrix const & mat1, Matrix const & mat2)
{
    return mat1.nrow() * mat1.ncol() * mat2.ncol();
}

/*
 * Use MKL for the matrix matrix multiplication.
 */
Matrix multiply_mkl(Matrix const & mat1, Matrix const & mat2){
    mkl_set_num_threads(1);

    Matrix ret(mat1.nrow(), mat2.ncol());

    StopWatch sw;

    cblas_dgemm(
        CblasRowMajor /* const CBLAS_LAYOUT Layout */
      , CblasNoTrans /* const CBLAS_TRANSPOSE transa */
      , CblasNoTrans /* const CBLAS_TRANSPOSE transb */
      , mat1.nrow() /* const MKL_INT m */
      , mat2.ncol() /* const MKL_INT n */
      , mat1.ncol() /* const MKL_INT k */
      , 1.0 /* const double alpha */
      , mat1.m_buffer /* const double *a */
      , mat1.ncol() /* const MKL_INT lda */
      , mat2.m_buffer /* const double *b */
      , mat2.ncol() /* const MKL_INT ldb */
      , 0.0 /* const double beta */
      , ret.m_buffer /* double * c */
      , ret.ncol() /* const MKL_INT ldc */
    );

    ret.elapsed() = sw.lap();
    ret.nflo() = calc_nflo(mat1, mat2);

    return ret;
}

/*
 * Indirect naive matrix matrix multiplication.
 */
Matrix multiply_naive(Matrix const & mat1, Matrix const & mat2){
    validate_multiplication(mat1, mat2);

    Matrix ret(mat1.nrow(), mat2.ncol());

    StopWatch sw;

    for (size_t i=0; i<mat1.nrow(); ++i)
    {
        for (size_t k=0; k<mat2.ncol(); ++k)
        {
            double v = 0;
            for (size_t j=0; j<mat1.ncol(); ++j)
            {
                v += mat1(i,j) * mat2(j,k);
            }
            ret(i,k) = v;
        }
    }

    ret.elapsed() = sw.lap();
    ret.nflo() = calc_nflo(mat1, mat2);

    return ret;
}


void initialize(Matrix & mat)
{
    for (size_t i=0; i<mat.nrow(); ++i)
    {
        for (size_t j=0; j<mat.ncol(); ++j)
        {
            if (0 == i%2)
            {
                mat(i, j) = j;
            }
            else
            {
                mat(i, j) = mat.ncol() - 1 - j;
            }
        }
    }
}
