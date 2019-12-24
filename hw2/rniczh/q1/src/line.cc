#include "line.hpp"

#include <stdexcept>
#include <iostream>

// copy constructor
Line::Line(Line const &line) : sz_{line.sz_}, location_{line.location_} {}

// move constructor
Line::Line(Line &&line) {
  std::swap(sz_, line.sz_);
  location_.swap(line.location_);
}

// copy assignment
Line &Line::operator=(Line const &rhs) {
  if (&rhs == this) return *this; // check the self copy assignment
  Line tmp {rhs};
  std::swap(tmp, *this);
  return *this;
}

// move assignment
Line &Line::operator=(Line &&rhs) {
  if (&rhs == this) return *this;
  std::swap(sz_, rhs.sz_);
  location_.swap(rhs.location_);
  return *this;
}
