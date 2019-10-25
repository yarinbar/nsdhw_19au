
#include <iostream>
#include <iomanip>
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

private:

    size_t m_nrow;
    size_t m_ncol;
    double * m_buffer;

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

Matrix multiply_naive(const Matrix& a, const Matrix& b){

    // check for dim compatibility
    if(a.ncol() != b.nrow())
        throw "sizes dont match!";

    Matrix c(a.nrow(), b.ncol());

    for(unsigned int i = 0; i < a.ncol(); ++i)
        for(unsigned int j = 0; j < a.ncol(); ++j) {
            c(i, j) = 0;
            for(unsigned int k = 0; k < a.ncol(); k++)
                c(i, j) += a(i, k) * b(k, j);
        }
    return c;
}

PYBIND11_MODULE(_matrix, m) {
    py::class_<Matrix>(m, "Matrix")
        .def(py::init<size_t, size_t>());

    m.def("multiply_naive", &multiply_naive, "");
}




int main(int argc, char ** argv)
{
    size_t width = 5;

    Matrix matrix(width, width);

    work(matrix);

    std::cout << "matrix:";
    for (size_t i=0; i<matrix.nrow(); ++i) // the i-th row
    {
        std::cout << std::endl << " ";
        for (size_t j=0; j<matrix.ncol(); ++j) // the j-th column
        {
            std::cout << " " << std::setfill('0') << std::setw(2)
                      << matrix(i, j);
        }
    }
    std::cout << std::endl;

    return 0;
}

// vim: set ff=unix fenc=utf8 et sw=4 ts=4 sts=4:
