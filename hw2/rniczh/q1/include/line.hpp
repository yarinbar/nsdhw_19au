#pragma once
#include <cstddef>
#include <vector>
#include <iostream>

class Line
{
public:
  Line() = default;
  Line(Line const &);
  Line(Line &&);
  Line &operator=(Line const &);
  Line &operator=(Line &&);
  Line(size_t size) :sz_ {size}, location_{std::vector<std::pair<float, float>>(size)} {};
  virtual ~Line() = default;
  size_t size() const { return sz_; };
  float const &x(size_t it) const { return location_.at(it).first; };
  float       &x(size_t it)       { return location_.at(it).first; };
  float const &y(size_t it) const { return location_.at(it).second; };
  float       &y(size_t it)       { return location_.at(it).second; };
  
private:
  size_t sz_ = 0;
  std::vector<std::pair<float, float>> location_;
}; /* end class Line */
