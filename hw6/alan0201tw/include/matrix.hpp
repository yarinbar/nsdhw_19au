#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

namespace py = pybind11;

class Matrix
{
public:

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

    Matrix(size_t nrow, size_t ncol, double* vec)
      : m_nrow(nrow), m_ncol(ncol)
    {
        reset_buffer(nrow, ncol);
        m_buffer = vec;
    }

    Matrix & operator=(std::vector<double> const & vec)
    {
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
    {
        reset_buffer(other.m_nrow, other.m_ncol);
        for (size_t i=0; i<m_nrow; ++i)
        {
            for (size_t j=0; j<m_ncol; ++j)
            {
                (*this)(i,j) = other(i,j);
            }
        }
    }

    Matrix & operator=(Matrix const & other)
    {
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
        return *this;
    }

    Matrix(Matrix && other)
      : m_nrow(other.m_nrow), m_ncol(other.m_ncol)
    {
        reset_buffer(0, 0);
        std::swap(m_nrow, other.m_nrow);
        std::swap(m_ncol, other.m_ncol);
        std::swap(m_buffer, other.m_buffer);
    }

    Matrix & operator=(Matrix && other)
    {
        if (this == &other) { return *this; }
        reset_buffer(0, 0);
        std::swap(m_nrow, other.m_nrow);
        std::swap(m_ncol, other.m_ncol);
        std::swap(m_buffer, other.m_buffer);
        return *this;
    }

    ~Matrix()
    {
        reset_buffer(0, 0);
    }

    double   operator() (size_t row, size_t col) const { return m_buffer[index(row, col)]; }
    double & operator() (size_t row, size_t col)       { return m_buffer[index(row, col)]; }

    size_t nrow() const { return m_nrow; }
    size_t ncol() const { return m_ncol; }

    size_t size() const { return m_nrow * m_ncol; }
    double* data() const { return m_buffer; }
    
    // Reference :
    // https://github.com/pybind/pybind11/issues/1042
    py::array_t<double> numpy_data() const
    {
        auto capsule = py::capsule(m_buffer, [](void *v) 
        {
            // delete reinterpret_cast<std::vector<int>*>(v); 
        });
        return py::array_t<double>({m_nrow, m_ncol}, m_buffer, capsule);
    }

    double buffer(size_t i) const { return m_buffer[i]; }
    std::vector<double> buffer_vector() const { return std::vector<double>(m_buffer, m_buffer+size()); }

    std::string ToString();

    bool operator==(const Matrix &other) const
    {
        if(nrow() == other.nrow() && ncol() == other.ncol())
        {
            for (size_t i = 0; i < nrow(); ++i)
            {
                for (size_t k = 0; k < ncol(); ++k)
                {
                    if(this->operator()(i, k) != other(i, k))
                    {
                        return false;
                    }
                }
            }

            return true;
        }

        return false;
    }

    bool operator!=(const Matrix &other) const
    {
        return !(*this == other);
    }

private:

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
    // py::array_t<double> m_buffer;
    double * m_buffer = nullptr;

    // friend free function for private field access
    friend Matrix multiply_mkl(Matrix const & mat1, Matrix const & mat2);
    friend Matrix multiply_tile(Matrix const & mat1, Matrix const & mat2, const unsigned int tsize);
    template<size_t N> friend struct Block;
    template<size_t N> friend struct Tiler;
};
