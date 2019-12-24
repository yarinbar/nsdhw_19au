#include "vector2d.hpp"

Vector2D Vector2D::operator-() const {
  return Vector2D(-x(), -y());
}

Vector2D Vector2D::operator+() const {
  return Vector2D(+x(), +y());
}

Vector2D &Vector2D::operator+=(Vector2D const& rhs) {
  x() += rhs.x();
  y() += rhs.y();
  return *this;
}

Vector2D &Vector2D::operator-=(Vector2D const& rhs) {
  x() -= rhs.x();
  y() -= rhs.y();
  return *this;
}

Vector2D &Vector2D::operator*=(float const& scalar) {
  x() *= scalar;
  y() *= scalar;
  return *this;
}

Vector2D Vector2D::operator+(Vector2D const& rhs) const {
  return Vector2D(*this) += rhs;
}

Vector2D Vector2D::operator-(Vector2D const& rhs) const {
  return Vector2D(*this) -= rhs;
}


Vector2D Vector2D::operator*(float const& scalar) const {
  return Vector2D(*this) *= scalar;
}

float Vector2D::operator*(Vector2D const& vec) const {
  return x() * vec.x() + y() * vec.y();
}

Vector2D operator*(float scalar, const Vector2D &rhs) {
  return rhs * scalar;
}

float Vector2D::norm() const {
  return sqrt(x() * x() + y() * y());
}

float Vector2D::angle(Vector2D const &rhs) {
  return acos((*this * rhs)/(this->norm() * rhs.norm()));
}

std::string Vector2D::show() const {
  return "[" + std::to_string(x()) + ", " + std::to_string(y()) + "]";
}
