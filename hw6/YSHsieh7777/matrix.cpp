#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

#include <vector>
#include <utility>
#include <stddef.h>

#include <mkl.h>

namespace py = pybind11;

class Matrix {
public:

    Matrix(size_t nrow, size_t ncol)
      : m_nrow(nrow), m_ncol(ncol)
    {
        reset_buffer(nrow, ncol);
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

    ~Matrix()
    {
        reset_buffer(0, 0);
    }

    double   operator() (size_t row, size_t col) const { return m_buffer[index(row, col)]; }
    double & operator() (size_t row, size_t col)       { return m_buffer[index(row, col)]; }

    bool operator== (const Matrix &other) const
    {
        if(!((m_nrow == other.nrow()) && (m_ncol == other.ncol())))
        {
            return false;
        }      
        else
        {
            size_t s = this->size();
            for(size_t i=0; i<s; i++)
            {
                if(this->buffer(i) != other.buffer(i))
                    return false;
            }
            return true;
        }        
    }
    bool operator!= (const Matrix &other) const
    {
        return !(*this == other);
    }

    size_t nrow() const { return m_nrow; }
    size_t ncol() const { return m_ncol; }
    double * buffer() const { return m_buffer; }

    py::array_t<double> array()
    {
        auto capsule = py::capsule(m_buffer, [](void *p) { /*delete reinterpret_cast<double *>(p);*/ });
        return py::array_t<double>({m_nrow, m_ncol}, m_buffer, capsule);
    }

    size_t size() const { return m_nrow * m_ncol; }
    double buffer(size_t i) const { return m_buffer[i]; }
    std::vector<double> buffer_vector() const { return std::vector<double>(m_buffer, m_buffer+size()); }

private:

    size_t index(size_t row, size_t col) const
    {
        return m_ncol*row + col;
    }

    void reset_buffer(size_t nrow, size_t ncol)
    {
        if (m_buffer) { mkl_free(m_buffer); }
        const size_t nelement = nrow * ncol;
        if (nelement) { m_buffer = (double *)mkl_malloc( nrow*ncol*sizeof( double ), 64 ); }
        else          { m_buffer = nullptr; }
        // if (m_buffer) { delete[] m_buffer; }
        // const size_t nelement = nrow * ncol;
        // if (nelement) { m_buffer = new double[nelement]; }
        // else          { m_buffer = nullptr; }
        m_nrow = nrow;
        m_ncol = ncol;
    }

    size_t m_nrow = 0;
    size_t m_ncol = 0;
    double * m_buffer = nullptr;

};

/*
 * Naive matrix matrix multiplication.
 */
Matrix multiply_naive(Matrix const & mat1, Matrix const & mat2)
{
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range(
            "the number of first matrix column "
            "differs from that of second matrix row");
    }

    Matrix ret(mat1.nrow(), mat2.ncol());

    for (size_t i=0; i<ret.nrow(); ++i)
    {
        for (size_t k=0; k<ret.ncol(); ++k)
        {
            double v = 0;
            for (size_t j=0; j<mat1.ncol(); ++j)
            {
                v += mat1(i,j) * mat2(j,k);
            }
            ret(i,k) = v;
        }
    }

    return ret;
}

Matrix multiply_tile(Matrix const & mat1, Matrix const & mat2, size_t t_size)
{
    if (mat1.ncol() != mat2.nrow())
    {
        throw std::out_of_range(
            "the number of first matrix column "
            "differs from that of second matrix row");
    }

    size_t i_size = mat1.nrow();
    size_t j_size = mat1.nrow();
    size_t k_size = mat2.ncol();

    Matrix ret(i_size, k_size);

    for(size_t i=0; i<i_size; i++)
    {
        for(size_t k=0; k<k_size; k++)
        {
            ret(i, k) = 0;
        }
    }

    for(size_t i=0; i<i_size; i+=t_size)
    {
        for(size_t k=0; k<k_size; k+=t_size)
        {
            for(size_t j=0; j<j_size; j+=t_size)
            {
                size_t it_max = i+t_size < i_size ? i+t_size : i_size;
                for(size_t it=i; it<it_max; it++)
                {
                    size_t kt_max = k+t_size < k_size ? k+t_size : k_size;
                    for(size_t kt=k; kt<kt_max; kt++)
                    {
                        size_t jt_max = j+t_size < j_size ? j+t_size : j_size;
                        double v = 0;
                        for(size_t jt=j; jt<jt_max; jt++)
                        {
                            v += mat1(it,jt) * mat2(jt,kt);
                        }
                        ret(it,kt) += v;
                    }
                }
            }
        }
    }
    return ret;
}

Matrix multiply_mkl(Matrix const & mat1, Matrix const & mat2)
{
    float alpha = 1.0; 
    float beta = 0.0;
    Matrix ret(mat1.nrow(), mat2.ncol());
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                mat1.nrow(), mat2.ncol(), mat1.ncol(), alpha, 
                mat1.buffer(), mat1.ncol(), mat2.buffer(), mat2.ncol(), 
                beta, ret.buffer(), mat2.ncol());

    return ret;
}

PYBIND11_MODULE(_matrix, m) {
    m.doc() = "pybind11 matrix class";      // module doc string

    py::class_<Matrix>(m, "Matrix")
        .def(py::init<size_t, size_t>())
        .def("__getitem__", [](const Matrix &m_matrix, std::pair<size_t, size_t> index) {
            return m_matrix(index.first, index.second);
        })
        .def("__setitem__", [](Matrix &m_matrix, std::pair<size_t, size_t> index, double value) {
            m_matrix(index.first, index.second) = value;
        })
        .def_property_readonly("nrow", &Matrix::nrow)
        .def_property_readonly("ncol", &Matrix::ncol)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def_property_readonly("array", &Matrix::array);

    m.def("multiply_naive",                              // function name
        &multiply_naive,                               // function pointer
        "A function which multiply two matrices in naive way." //function doc string
        );
    m.def("multiply_mkl",                              // function name
        &multiply_mkl,                               // function pointer
        "A function which multiply two matrices with mkl library." //function doc string
        );
    m.def("multiply_tile",                              // function name
        &multiply_tile,                               // function pointer
        "A function which multiply two matrices with tile size." //function doc string
        );
}
