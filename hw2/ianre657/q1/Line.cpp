#include "Line.hpp"
#include <stdexcept>
#include <vector>

size_t Line::check_range(size_t it) const {
  if (it >= data_ptr->size()) {
    throw std::out_of_range("Line index out of range");
  }
  return it;
}

size_t get_x_shift(size_t it) { return it << 1; }

size_t get_y_shift(size_t it) { return (it << 1) + 1; }

Line::Line() : m_size(0) { data_ptr = new std::vector<float>{}; }

// copy constructor
Line::Line(Line const &other) : m_size(other.m_size) {
  data_ptr =
      new std::vector<float>{other.data_ptr->begin(), other.data_ptr->end()};
}

// move constructor
Line::Line(Line &&other) : m_size(other.m_size) { data_ptr = other.data_ptr; }

Line &Line::operator=(Line const &other) {
  data_ptr =
      new std::vector<float>{other.data_ptr->begin(), other.data_ptr->end()};
  m_size = other.m_size;
  return *this;
}

Line &Line::operator=(Line &&other) {
  data_ptr = other.data_ptr;
  m_size = other.m_size;
  return *this;
}

Line::Line(size_t size) {
  data_ptr = new std::vector<float>{};
  data_ptr->resize(2 * size);
  m_size = size;
}

Line::~Line() { delete data_ptr; }

size_t Line::size() const { return m_size; }

float &Line::x(size_t it) const {
  size_t idx = get_x_shift(it);
  check_range(idx);
  return (*data_ptr)[idx];
}

float &Line::x(size_t it) { return (*data_ptr)[check_range(get_x_shift(it))]; }

float &Line::y(size_t it) const {
  return (*data_ptr)[check_range(get_y_shift(it))];
}

float &Line::y(size_t it) { return (*data_ptr)[check_range(get_y_shift(it))]; }
