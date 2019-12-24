#ifndef MATRIX_MUL
#define MATRIX_MUL

#include <ostream>
#include <stdexcept>
#include <vector>


class Matrix {
public:
  Matrix() = delete;
  ~Matrix() = default;
  Matrix(size_t row, size_t col) : nrow_(row), ncol_(col) {
    buf_.resize(row * col);
  }
  Matrix(Matrix const &rhs) {
    this->ncol_ = rhs.ncol_;
    this->nrow_ = rhs.nrow_;
    this->buf_ = rhs.buf_;
  }
  Matrix(Matrix &&rhs) {
    this->ncol_ = rhs.ncol_;
    this->nrow_ = rhs.nrow_;
    this->buf_ = std::move(rhs.buf_);
  }
  Matrix &operator=(Matrix const &rhs) {
    this->ncol_ = rhs.ncol_;
    this->nrow_ = rhs.nrow_;
    this->buf_ = rhs.buf_;
    return *this;
  }
  Matrix &operator=(Matrix &&rhs) {
    this->ncol_ = rhs.ncol_;
    this->nrow_ = rhs.nrow_;
    this->buf_ = std::move(rhs.buf_);
    return *this;
  }
  const double &get(size_t i, size_t j) const {
    if (i >= nrow_ || j >= ncol_)
      throw std::invalid_argument("index exceed dimension");
    return buf_.at(i * ncol_ + j);
  }
  double &set(size_t i, size_t j) {
    if (i >= nrow_ || j >= ncol_)
      throw std::invalid_argument("index exceed dimension");
    return buf_.at(i * ncol_ + j);
  }
  const double &operator()(size_t i, size_t j) const {
    return buf_[i * ncol_ + j];
  }
  double &operator()(size_t i, size_t j) { return buf_[i * ncol_ + j]; }
  bool operator==(Matrix const &n) const {
    return (ncol_ == n.ncol_) && (nrow_ == n.nrow_) && (buf_ == n.buf_);
  }
  double *data() { return buf_.data(); }
  const size_t &nrow() const { return nrow_; }
  const size_t &ncol() const { return ncol_; }

private:
  size_t nrow_;
  size_t ncol_;
  std::vector<double> buf_;
};
Matrix multiply_naive(Matrix &lhs, Matrix &rhs);
Matrix multiply_mkl(Matrix &lhs, Matrix &rhs);
Matrix multiply_tile(Matrix &lhs, Matrix &rhs, size_t tsize);
std::ostream &operator<<(std::ostream &os, Matrix const &m);


#endif
