#include "Line.h"

Line::Line() {}
Line::Line(Line const &rhs) { points = rhs.points; }
Line::Line(Line &&rhs) { std::swap(points, rhs.points); }
Line &Line::operator=(Line const &rhs) {
    points = rhs.points;
    return *this;
}
Line &Line::operator=(Line &&rhs) {
    std::swap(points, rhs.points);
    return *this;
}
Line::Line(size_t size) { points = std::vector<Point>(size, {0.0, 0.0}); }
Line::~Line() {}
size_t Line::size() const { return points.size(); };
const float &Line::x(size_t it) const { return points[it].first; }
float &Line::x(size_t it) { return points[it].first; }
const float &Line::y(size_t it) const { return points[it].second; }
float &Line::y(size_t it) { return points[it].second; }
