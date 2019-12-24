#include "Vec2.h"
#include <cmath>
#include <assert.h>
#include <sstream>

Vec2::Vec2(){
  x_ = y_ = 0;
}

Vec2::Vec2(Vec2 const &v){
  x_ = v.x(), y_ = v.y();
}

Vec2::Vec2(Vec2 &&v){
  x_ = v.x(), y_ = v.y();
}

Vec2 Vec2::copy(Vec2 const &v){
  return Vec2(v);
}

Vec2 &Vec2::operator=(Vec2 const &){
  return *this;
}

Vec2 &Vec2::operator=(Vec2 &&){
  return *this;
}

Vec2::Vec2(double x, double y){
  x_ = x, y_ = y;
}

Vec2::~Vec2(){ }

double Vec2::size() const{
  return sqrt((x_ * x_) + (y_ * y_));
}

double Vec2::dot(Vec2 v) const{
  return x_ * v.x() + y_ * v.y();
}

double Vec2::angle(Vec2 v) const{
  double mul_length = this->size() * v.size();

  if(mul_length == 0)
    throw std::runtime_error("cannot not calculate angle of zero length vec2");
  return acos(this->dot(v) / mul_length);
}

double const &Vec2::x() const{
  return x_;
}

double &Vec2::x(double v){
  return x_ = v;
}

double const &Vec2::y() const{
  return y_;
}

double &Vec2::y(double v){
  return y_ = v;
}

std::string Vec2::repr(){
  std::ostringstream ss;
  ss << "<vec2.vec2 (" << x_ << ", " << y_ << ") at "
     << std::hex << (void*)this << ">";
  return ss.str();
}
