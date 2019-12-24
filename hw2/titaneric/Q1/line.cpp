#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iostream>
#include <utility>

#include "line.h"

Line::Line(size_t size) :list(size, make_pair(0,0)){}
/* Define the copy constructor */
Line::Line(Line const &other)
{
    list = vector<pair<float, float>>(other.list);
}

/* Define the move constructor */
Line::Line(Line &&other)
{
    std::swap(other.list, list);
}

/* Define the copy assignment operator */
Line &Line::operator=(Line const &other)
{
    if (this == &other)
    {
        return *this;
    } // don't move to self.
    // This part is the same as what's in the copy constructor.
    list = vector<pair<float, float>>(other.list);
    return *this;
}

/* Define the move assignment operator */
Line &Line::operator=(Line &&other)
{
    if (this == &other)
    {
        return *this;
    } // don't move to self.
    std::swap(other.list, list);
    return *this;
}

void Line::check_range(size_t it) const
{
    if (it >= list.size())
    {
        throw std::out_of_range("Line index out of range");
    }
}
size_t Line::size() const
{
    return list.size();
}
float Line::x(const size_t it) const
{
    check_range(it);
    return list[it].first;
}

float &Line::x(const size_t it)
{
    check_range(it);
    return list[it].first;
}

float Line::y(const size_t it) const
{
    check_range(it);
    return list[it].second;
}

float &Line::y(const size_t it)
{
    check_range(it);
    return list[it].second;
}