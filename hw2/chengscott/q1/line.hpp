#pragma once
#include <array>
#include <vector>

class Line
{
public:
    using ValueType = float;
    Line() = default;
    explicit Line(size_t size) { line_.resize(size); }
    Line(Line const &line)
    {
        auto tmp{line.line_};
        line_.swap(tmp);
    }
    Line(Line &&line) : Line() { swap(line); }
    Line &operator=(Line const &line)
    {
        Line(line).swap(*this);
        return *this;
    }
    Line &operator=(Line &&line)
    {
        swap(line);
        return *this;
    }
    ~Line() = default;
    size_t size() const { return line_.size(); }
    const ValueType &x(size_t it) const { return line_.at(it).at(0); }
    ValueType &x(size_t it) { return line_.at(it).at(0); }
    const ValueType &y(size_t it) const { return line_.at(it).at(1); }
    ValueType &y(size_t it) { return line_.at(it).at(1); }
    void swap(Line &line) { std::swap(line_, line.line_); }

private:
    std::vector<std::array<ValueType, 2>> line_;
};
