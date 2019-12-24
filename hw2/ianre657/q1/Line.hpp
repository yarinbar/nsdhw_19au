#pragma once

#include <vector>

class Line {
public:
  Line();
  Line(Line const &);
  Line(Line &&);
  Line &operator=(Line const &);
  Line &operator=(Line &&);
  Line(size_t size);
  ~Line();
  size_t size() const;

  float &x(size_t it) const;
  float &x(size_t it);
  float &y(size_t it) const;
  float &y(size_t it);

private:
  size_t m_size;
  size_t check_range(size_t it) const;
  std::vector<float> *data_ptr;
};