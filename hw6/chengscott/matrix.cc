#include <mkl.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <vector>
namespace py = pybind11;

class Matrix {
public:
  Matrix(size_t nrow, size_t ncol)
      : nrow_(nrow), ncol_(ncol), data_(nrow * ncol) {}
  Matrix(const Matrix &) = default;
  Matrix(Matrix &&) noexcept = default;
  ~Matrix() = default;
  Matrix &operator=(const Matrix &) = default;
  Matrix &operator=(Matrix &&) noexcept = default;

  double operator()(size_t i, size_t j) const {
    if (i >= nrow_ || j >= ncol_) {
      throw py::index_error();
    }
    return data_[i * ncol_ + j];
  }
  double &operator()(size_t i, size_t j) {
    if (i >= nrow_ || j >= ncol_) {
      throw py::index_error();
    }
    return data_[i * ncol_ + j];
  }

  size_t nrow() const { return nrow_; }
  size_t ncol() const { return ncol_; }

  bool operator==(const Matrix &rhs) const noexcept {
    return nrow_ == rhs.nrow_ && ncol_ == rhs.ncol_ && data_ == rhs.data_;
  }
  bool operator!=(const Matrix &rhs) const noexcept { return !(*this == rhs); }

  double *get_raw_data() noexcept { return data_.data(); }
  const double *get_raw_data() const noexcept { return data_.data(); }

  friend Matrix multiply_naive(const Matrix &, const Matrix &);
  friend Matrix multiply_mkl(const Matrix &, const Matrix &);

private:
  size_t nrow_, ncol_;
  std::vector<double> data_;
};

Matrix multiply_naive(const Matrix &lhs, const Matrix &rhs) {
  const size_t m = lhs.nrow_, n = lhs.ncol_, l = rhs.ncol_;
  Matrix ret(m, l);
  std::fill(std::begin(ret.data_), std::end(ret.data_), 0);
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < l; ++j)
      for (size_t k = 0; k < n; ++k)
        ret(i, j) += lhs(i, k) * rhs(k, j);
  return ret;
}

Matrix multiply_mkl(const Matrix &lhs, const Matrix &rhs) {
  // assert(lhs.ncol() == rhs.nrow());
  Matrix ret(lhs.nrow(), rhs.ncol());
  std::fill(std::begin(ret.data_), std::end(ret.data_), 0);
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, lhs.nrow(), rhs.ncol(),
              lhs.ncol(), 1., lhs.data_.data(), lhs.ncol(), rhs.data_.data(),
              rhs.ncol(), 1., ret.data_.data(), ret.ncol());
  return ret;
}

PYBIND11_MODULE(_matrix, m) {
  m.doc() = "class Matrix";
  m.def("multiply_naive", &multiply_naive);
  m.def("multiply_mkl", &multiply_mkl);
  py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
      .def(py::init<size_t, size_t>())
      .def_property_readonly("nrow", &Matrix::nrow)
      .def_property_readonly("ncol", &Matrix::ncol)
      .def("__eq__",
           [](const Matrix &lhs, const Matrix &rhs) { return lhs == rhs; })
      .def("__ne__",
           [](const Matrix &lhs, const Matrix &rhs) { return lhs != rhs; })
      .def("__getitem__",
           [](const Matrix &m, std::pair<size_t, size_t> i) {
             return m(i.first, i.second);
           })
      .def("__setitem__", [](Matrix &m, std::pair<size_t, size_t> i,
                             double v) { m(i.first, i.second) = v; })
      .def_buffer([](Matrix &m) -> py::buffer_info {
        return py::buffer_info(m.get_raw_data(), sizeof(double),
                               py::format_descriptor<double>::format(), 2,
                               {m.nrow(), m.ncol()},
                               {sizeof(double) * m.ncol(), sizeof(double)});
      })
      .def_property("array",
                    [](Matrix &m) {
                      return py::array(
                          {m.nrow(), m.ncol()}, m.get_raw_data(),
                          py::capsule(m.get_raw_data(), [](void *) {}));
                    },
                    nullptr);
  ;
}