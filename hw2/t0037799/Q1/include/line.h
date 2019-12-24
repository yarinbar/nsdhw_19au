#ifndef LINE_H
#define LINE_H

#include <vector>

class Line {
public:
  Line() = default;
  Line(Line const &);
  Line(Line &&);
  Line &operator=(Line const &);
  Line &operator=(Line &&);
  Line(size_t size);
  ~Line() = default;
  size_t size() const { return datas_.size(); };
  const float &x(size_t it) const { return datas_.at(it).first; }
  float &x(size_t it) { return datas_.at(it).first; }
  const float &y(size_t it) const { return datas_.at(it).second; }
  float &y(size_t it) { return datas_.at(it).second; }

private:
  std::vector<std::pair<float, float>> datas_;
};

#endif
