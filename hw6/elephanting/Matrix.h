#include <iostream>
#include <string>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

class Matrix
{
public:
    Matrix(size_t, size_t);
    ~Matrix();
    Matrix& operator=(Matrix const &);
    bool operator==(Matrix const &);
    Matrix(const Matrix &);
    double& operator()(size_t const, size_t const);
    double operator()(size_t const, size_t const) const;
    double* buffer() const;
    double* get_buffer();
    std::string repr() const;
    size_t nrow() const;
    size_t ncol() const;
    py::array_t<double> array();
private:
    size_t m_ncol;
    size_t m_nrow;
    double* m_buffer = nullptr;
    void resetbuf(size_t, size_t);
};

Matrix mkl_MM(Matrix const &, Matrix const &);

bool isequal(Matrix &, Matrix &);

Matrix multiply_naive(Matrix const &, Matrix const &);

struct Block;

void block_assign(Matrix const &, size_t, size_t, Block &);

void matrix_assign(Matrix &, size_t, size_t, Block const &);

Matrix multiply_tile(Matrix const &, Matrix const &, size_t);

Matrix operator*(Matrix const &, Matrix const &);

/* broadcasting */
Matrix operator*(double const &, Matrix const &);

Matrix operator*(Matrix const &, double const &);
