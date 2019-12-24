#ifdef PYTHON_LIB
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/embed.h>
#endif
#include <mkl.h>

#include <iostream>
#include <utility>
#include <cmath>
#include <map>
#include <vector>
#include <stdexcept>
#include <numeric>
#include <algorithm>

using namespace std;

#ifdef PYTHON_LIB
namespace py = pybind11;
#endif

class Matrix
{
public:
    Matrix() = default;
    Matrix(size_t r, size_t c) : m_nrow(r), m_ncol(c)
    {
        allocate_buffer(nrow(), ncol());
    };
    Matrix(size_t r, size_t c, vector<double> const &vec) : m_nrow(r), m_ncol(c)
    {
        allocate_buffer(nrow(), ncol());

        size_t vec_i = 0;
        for (size_t i = 0; i < nrow(); i++)
        {
            for (size_t j = 0; j < ncol(); j++)
            {
                // data[i][j] = vec[vec_i++];
                this->operator()(i, j) = vec[vec_i++];
            }
        }
    };
    Matrix(size_t r, size_t c, double *flat) : m_nrow(r), m_ncol(c)
    {
        allocate_buffer(nrow(), ncol());
        size_t vec_i = 0;
        for (size_t i = 0; i < nrow(); i++)
        {
            for (size_t j = 0; j < ncol(); j++)
            {
                // data[i][j] = flat[vec_i++];
                this->operator()(i, j) = flat[vec_i++];
            }
        }
    };
    Matrix(vector<vector<double>> const &d) : m_nrow(d.size()), m_ncol((d.size()) ? d[0].size() : 0)
    {
        allocate_buffer(nrow(), ncol());
        for (size_t i = 0; i < nrow(); i++)
        {
            for (size_t j = 0; j < ncol(); j++)
            {
                this->operator()(i, j) = d[i][j];
            }
        }
    }
    Matrix(const Matrix &target) : m_nrow(target.m_nrow), m_ncol(target.m_ncol)
    {
        allocate_buffer(nrow(), ncol());
        for (size_t i = 0; i < nrow(); i++)
        {
            for (size_t j = 0; j < ncol(); j++)
            {
                this->operator()(i, j) = target(i, j);
            }
        }
    };
    Matrix(Matrix &&target) : m_nrow(target.m_nrow), m_ncol(target.m_ncol)
    {
        allocate_buffer(nrow(), ncol());
        swap(m_nrow, target.m_nrow);
        swap(m_ncol, target.m_ncol);
        swap(m_data, target.m_data);
    }
    ~Matrix()
    {
        allocate_buffer(0, 0);
    };
    Matrix &operator=(const Matrix &target)
    {
        if (this == &target)
        {
            return *this;
        }

        if (m_nrow != target.m_nrow || m_ncol != target.m_ncol)
        {
            allocate_buffer(target.m_nrow, target.m_ncol);
        }
        for (size_t i = 0; i < nrow(); i++)
        {
            for (size_t j = 0; j < ncol(); j++)
            {
                this->operator()(i, j) = target(i, j);
            }
        }

        return *this;
    }
    Matrix &operator=(Matrix &&target)
    {
        if (this == &target)
        {
            return *this;
        }
        allocate_buffer(0, 0);
        std::swap(m_nrow, target.m_nrow);
        std::swap(m_ncol, target.m_ncol);
        std::swap(m_data, target.m_data);
        return *this;
    }
    double operator()(size_t row, size_t col) const { return m_data[flatten_index(row, col)]; }
    double &operator()(size_t row, size_t col)
    {
        return m_data[flatten_index(row, col)];
    }
    bool operator==(const Matrix &other) const
    {
        if (is_same_size(other))
        {
            for (size_t i = 0; i < nrow(); i++)
            {
                for (size_t j = 0; j < ncol(); j++)
                {
                    if (this->operator()(i, j) != other(i, j))
                        return false;
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
    bool is_same_size(const Matrix &other) const
    {
        return (nrow() == other.nrow()) && (ncol() == other.ncol());
    }
    Matrix &operator+=(const Matrix &other)
    {
        if (!nrow() && !ncol())
        {
            *this = other;
            return *this;
        }
        if (!is_same_size(other))
        {
            throw out_of_range("Matrices can not be added!!!");
        }
        for (size_t i = 0; i < nrow(); i++)
        {
            for (size_t j = 0; j < ncol(); j++)
            {
                // data[i][j] += other(i, j);
                this->operator()(i, j) += other(i, j);
            }
        }
        return *this;
    }

    double *getFlatData() const
    {
        double *flat = new double[nrow() * ncol()];
        size_t k = 0;
        for (size_t i = 0; i < nrow(); i++)
        {
            for (size_t j = 0; j < ncol(); j++)
            {
                flat[k++] = this->operator()(i, j);
            }
        }
        return flat;
    }
    vector<vector<double>> loadBlock(size_t u, size_t v,
                                     size_t u_tile, size_t v_tile,
                                     bool column_major = false) const
    {
        if (column_major)
        {
            vector<vector<double>> block(v_tile, vector<double>(u_tile, 0));
            for (size_t i = 0; i < u_tile; i++)
            {
                for (size_t j = 0; j < v_tile; j++)
                {
                    // block[j][i] = data[u + i][v + j];
                    block[j][i] = this->operator()(u + i, v + j);
                }
            }
            return block;
        }
        else
        {
            vector<vector<double>> block(u_tile, vector<double>(v_tile, 0));
            for (size_t i = 0; i < u_tile; i++)
            {
                for (size_t j = 0; j < v_tile; j++)
                {
                    // block[i][j] = data[u + i][v + j];
                    block[i][j] = this->operator()(u + i, v + j);
                }
            }
            return block;
        }
    }

    void saveBlock(size_t u, size_t v,
                   Matrix const &block)
    {
        for (size_t i = 0; i < block.nrow(); i++)
        {
            for (size_t j = 0; j < block.ncol(); j++)
            {
                // data[u + i][v + j] = block(i, j);
                this->operator()(u + i, v + j) = block(i, j);
            }
        }
    }

#ifndef PYTHON_LIB
    friend ostream &operator<<(ostream &ostr, Matrix const &mat)
    {
        for (size_t i = 0; i < mat.nrow(); ++i)
        {
            ostr << endl
                 << " ";
            for (size_t j = 0; j < mat.ncol(); ++j)
            {
                ostr << mat(i, j);
            }
        }

        return ostr;
    }
#endif

public:
    size_t nrow() const { return m_nrow; };
    size_t ncol() const { return m_ncol; };

    double *data() const
    {
        return m_data;
    }
#ifdef PYTHON_LIB
public:
    py::array_t<double> get_array()
    {
        /* We MUST provide base to prevent the data copy
            BUT we don't need to delete the original buffer
            because we've already have destructor for it.
            Thanks to https://github.com/pybind/pybind11/issues/1042
            array source: https://github.com/pybind/pybind11/blob/master/include/pybind11/numpy.h
        */
        auto capsule = py::capsule(data(), [](void *v) { /*delete reinterpret_cast<double *>(v);*/ });
        py::array_t<double> numpy_data(
            {nrow(), ncol()},
            {sizeof(double) * ncol(),
             sizeof(double)},
            data(),
            capsule);
        return numpy_data;
    }

#endif

private:
    const size_t flatten_index(size_t i, size_t j) const
    {
        return i * ncol() + j;
    }
    void allocate_buffer(size_t nr, size_t nc)
    {
        if (m_data)
        {
            delete[] m_data;
        }
        m_ncol = nc;
        m_nrow = nr;
        m_data = new double[m_ncol * m_nrow];
    }

private:
    size_t m_nrow = 0;
    size_t m_ncol = 0;
    double *m_data = nullptr;
};

Matrix multiply_naive(Matrix const &A, Matrix const &B, bool column_major = false)
{

    double value = 0;
    if (!column_major)
    {
        if (A.ncol() != B.nrow())
        {
            throw out_of_range("Matrices can not be multiplied!!!");
        }
        Matrix result(A.nrow(), B.ncol());
        for (size_t i = 0; i < result.nrow(); i++)
        {
            for (size_t j = 0; j < result.ncol(); j++)
            {
                value = 0;
                for (size_t k = 0; k < A.ncol(); k++)
                {
                    value += A(i, k) * B(k, j);
                }
                result(i, j) = value;
                // cout << result(i, j) << endl;
            }
        }
        return result;
    }
    else
    {
        if (A.ncol() != B.ncol())
        {
            throw out_of_range("Matrices can not be multiplied!!!");
        }
        Matrix result(A.nrow(), B.nrow());
        for (size_t i = 0; i < result.nrow(); i++)
        {
            for (size_t j = 0; j < result.ncol(); j++)
            {
                value = 0;
                for (size_t k = 0; k < A.ncol(); k++)
                {
                    value += A(i, k) * B(j, k);
                }
                result(i, j) = value;
            }
        }
        return result;
    }
}

// https://software.intel.com/en-us/mkl-tutorial-c-multiplying-matrices-using-dgemm
Matrix multiply_mkl(Matrix const &A, Matrix const &B)
{
    size_t m = A.nrow(); // A: m x k
    size_t k = A.ncol(); // B: k x n
    size_t n = B.ncol(); // C: m x n
    double alpha = 1.0;
    double beta = 0.0;

    double *A_flat = A.getFlatData();
    double *B_flat = B.getFlatData();
    double *C_flat = new double[m * n];

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                m, n, k, alpha, A_flat, k, B_flat, n, beta, C_flat, n);

    Matrix result(m, n, C_flat);

    delete A_flat;
    delete B_flat;
    delete C_flat;

    return result;
}

bool is_on_edge(size_t it, size_t i_limit)
{
    return it == (i_limit - 1);
}

Matrix determine_source_of_block(size_t i, size_t j,
                                 size_t it, size_t jt,
                                 size_t i_limit, size_t j_limit,
                                 size_t rem_row, size_t rem_col,
                                 size_t tile_size, const Matrix &matrix,
                                 map<pair<size_t, size_t>, Matrix> &matrix_map,
                                 bool column_major = false)
{
    auto block_index = make_pair(i, j);
    auto find_block = matrix_map.find(block_index);

    if (find_block != matrix_map.end())
    {
        return find_block->second;
    }
    size_t load_row = (is_on_edge(it, i_limit) && rem_row) ? rem_row : tile_size;
    size_t load_col = (is_on_edge(jt, j_limit) && rem_col) ? rem_col : tile_size;

    matrix_map[block_index] = Matrix(matrix.loadBlock(i, j, load_row, load_col, column_major));
    return matrix_map[block_index];
}
size_t power_of_2_tile_size(size_t tile_size)
{
    double pos = floor(log2((double)tile_size));
    return (size_t)pow(2, pos);
}
Matrix multiply_tile(Matrix const &A, Matrix const &B, size_t tile_size)
{
    if (A.ncol() != B.nrow())
    {
        throw out_of_range("Matrices can not be multiplied!!!");
    }
    auto matrices_sizes = {A.nrow(), A.ncol(), B.ncol()};
    size_t min_size = *min_element(matrices_sizes.begin(), matrices_sizes.end());
    if (min_size < tile_size)
    {
        throw out_of_range("Tile size is too big!!!");
    }

#if defined(PYTHON_LIB) && !defined(PYTHON_TEST)
    tile_size = 64;
#else
    // find suitable tile size which is the power of 2
    size_t new_tile_size = power_of_2_tile_size(tile_size);
    tile_size = new_tile_size;
#endif

    div_t tile_rowA_result = div((int)A.nrow(), (int)tile_size);
    const size_t num_tile_rowA = (tile_rowA_result.rem) ? tile_rowA_result.quot + 1 : tile_rowA_result.quot;

    div_t tile_colA_result = div((int)A.ncol(), (int)tile_size);
    const size_t num_tile_colA = (tile_colA_result.rem) ? tile_colA_result.quot + 1 : tile_colA_result.quot;

    div_t tile_colB_result = div((int)B.ncol(), (int)tile_size);
    const size_t num_tile_colB = (tile_colB_result.rem) ? tile_colB_result.quot + 1 : tile_colB_result.quot;

    map<pair<size_t, size_t>, Matrix> matrix_mapA;
    map<pair<size_t, size_t>, Matrix> matrix_mapB;

    size_t rem_rowA = tile_rowA_result.rem;
    size_t rem_colA = tile_colA_result.rem;
    size_t rem_colB = tile_colB_result.rem;

    Matrix result(A.nrow(), B.ncol());
    size_t augA_col = num_tile_colA * tile_size;

    // empty block
    Matrix empty_block_source;
    // for blocking MM result
    Matrix block_result;

    Matrix blockA;
    Matrix blockB;

    for (size_t i = 0, it = 0; i < result.nrow(); i += tile_size, it++)
    {
        for (size_t j = 0, jt = 0; j < result.ncol(); j += tile_size, jt++)
        {
            block_result = empty_block_source;
            for (size_t k = 0, kt = 0; k < augA_col; k += tile_size, kt++)
            {
                blockA = determine_source_of_block(i, k,
                                                   it, kt,
                                                   num_tile_rowA, num_tile_colA,
                                                   rem_rowA, rem_colA,
                                                   tile_size, A, matrix_mapA);

                blockB = determine_source_of_block(k, j,
                                                   kt, jt,
                                                   num_tile_colA, num_tile_colB,
                                                   rem_colA, rem_colB,
                                                   tile_size, B, matrix_mapB, true);

                block_result += multiply_naive(blockA, blockB, true);
            }
            result.saveBlock(i, j, block_result);
        }
    }

    return result;
}

#ifdef PYTHON_LIB
PYBIND11_MODULE(_matrix, m)
{
    m.doc() = "pybind11 matrix plugin"; // optional module docstring
    m.def("multiply_mkl", &multiply_mkl, "A function which calculates product of two matrices using MKL");
    m.def("multiply_naive", &multiply_naive, "A function which calculates product of two matrices using naive method",
          py::arg("A"), py::arg("B"), py::arg("column_major") = false);
    m.def("multiply_tile", &multiply_tile, "A function which calculates product of two matrices using tile method");
    py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
        .def_buffer([](Matrix &m) -> py::buffer_info {
            return py::buffer_info(
                m.data(),                                /* Pointer to buffer */
                sizeof(double),                          /* Size of one scalar */
                py::format_descriptor<double>::format(), /* Python struct-style format descriptor */
                2,                                       /* Number of dimensions */
                {m.nrow(), m.ncol()},                    /* Buffer dimensions */
                {sizeof(double) * m.ncol(),              /* Strides (in bytes) for each index */
                 sizeof(double)});
        })
        .def(py::init<size_t, size_t>())
        .def(py::init<size_t, size_t, vector<double>>())
        .def(py::init<vector<vector<double>>>())
        .def_property_readonly("nrow", &Matrix::nrow)
        .def_property_readonly("ncol", &Matrix::ncol)
        .def_property_readonly("array", &Matrix::get_array)
        .def("__eq__", [](const Matrix &self, const Matrix &other) { return self == other; })
        .def("__ne__", [](const Matrix &self, const Matrix &other) { return self != other; })
        .def("__getitem__", [](const Matrix &self, pair<size_t, size_t> idx) { return self(idx.first, idx.second); })
        .def("__setitem__", [](Matrix &self, pair<size_t, size_t> idx, double val) { self(idx.first, idx.second) = val; });
}
#else
int main()
{
    Matrix A({{1, 0},
              {0, 1}});
    Matrix B({{1, 2},
              {3, 4}});

    cout << A << endl;
    cout << B << endl;
    cout << multiply_mkl(A, B) << endl;

    cout << "------\n";

    A = Matrix({{1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}});
    B = Matrix({{1, 2, 3},
                {4, 5, 6},
                {7, 8, 9}});
    cout << A << endl;
    cout << B << endl;
    cout << multiply_tile(A, B, 2) << endl;
    return 0;
}
#endif