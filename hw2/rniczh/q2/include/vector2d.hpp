#pragma once
#include <cstddef>
#include <string>
#include <cmath>

class Vector2D
{
public:
  Vector2D()                            = default;
  Vector2D(float x, float y) :x_ {x}, y_ {y} {};
  Vector2D(Vector2D const &)            = default;
  Vector2D(Vector2D &&)                 = default;
  Vector2D &operator=(Vector2D const &) = default;
  Vector2D &operator=(Vector2D &&)      = default;
  virtual ~Vector2D()                   = default;
  float const &x() const { return x_; };
  float       &x()       { return x_; };
  float const &y() const { return y_; };
  float       &y()       { return y_; };
  
  // vector operations
  Vector2D operator-() const;
  Vector2D operator+() const;
  Vector2D operator+(Vector2D const& rhs) const;
  Vector2D operator-(Vector2D const& rhs) const;
  Vector2D operator*(float const& scalar) const;
  friend Vector2D operator*(float scalar, const Vector2D &rhs);
  Vector2D &operator+=(Vector2D const& rhs);
  Vector2D &operator-=(Vector2D const& rhs);
  Vector2D &operator*=(float const& scalar);
  float operator*(Vector2D const& vec) const;
  float norm() const;
  float angle(Vector2D const &);

  std::string show() const;

private:
  float x_ = 0.0f, y_ = 0.0f;
}; /* end class Vector2D */
