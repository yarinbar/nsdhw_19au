#include<stddef.h>
#include<string>
#include<limits>

bool double_equals(double a, double b)
{
    return std::abs(a - b) < std::numeric_limits<double>::epsilon();
}

class Vec2{
  public:
    Vec2();
    Vec2(Vec2 const &);
    Vec2(Vec2       &&);
    Vec2 &operator=(Vec2 const &);
    Vec2 &operator=(Vec2       &&);
    Vec2(double , double);
    ~Vec2();
    static Vec2 copy(Vec2 const &);
    double size() const;
    double dot(Vec2) const;
    double angle(Vec2) const;
    double const &x() const;
    double &x(double);
    double const &y() const;
    double &y(double);
    std::string repr();
    Vec2 operator+(const Vec2 &v) const { return Vec2(x_ + v.x_, y_ + v.y_); }
    Vec2 operator-(const Vec2 &v) const { return Vec2(x_ - v.x_, y_ - v.y_); }
    Vec2 operator*(const double value) const { return Vec2(x_ * value, y_ * value); }
    Vec2& operator+=(const Vec2 &v) { x_ += v.x_; y_ += v.y_; return *this; }
    Vec2& operator*=(const double v) { x_ *= v; y_ *= v; return *this; }
    bool operator==(const Vec2 &v) const
    { return double_equals(x_, v.x_) && double_equals(y_, v.y_); }
    bool operator!=(const Vec2 &v) const { return !((*this) == v); }
  private:
    double x_, y_;
};

Vec2 operator*(const double value, const Vec2 &v){
  return Vec2(v.x() * value, v.y() * value);
}

Vec2 operator-(const Vec2 &v){
  return Vec2(-v.x(), -v.y());
}

Vec2 operator+(const Vec2 &v){
  return Vec2(v.x(), v.y());
}
