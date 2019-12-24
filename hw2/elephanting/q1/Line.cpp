#include <iostream>
#include <vector>
#include "Line.h"
using namespace std;

/* define the default constructor */
Line::Line()
{
    m_size = 0;
    ptr = nullptr;
}

/* define the custom constructor */
Line::Line(size_t size)
{
    m_size = size;
    ptr = new vector<float>;
    (*ptr).resize(size * 2);    
}

void Line::check_range(size_t it) const
{ 
    if (it >= m_size)
        throw std::out_of_range("Line index out of range");
}

/* setter x */
float &Line::x(size_t it)
{
    check_range(it);
    return (*ptr)[2 * it];
}

/* setter y */
float &Line::y(size_t it)
{
    check_range(it);
    return (*ptr)[2 * it + 1];
}

/* getter x */
float &Line::x(size_t it) const
{
    check_range(it);
    return (*ptr)[2 * it];
}

/* getter y */
float &Line::y(size_t it) const
{
    check_range(it);
    return (*ptr)[2 * it + 1];
}

/* define the copy constructor */
Line::Line(Line const &other)
{
    m_size = other.m_size;
    ptr = new vector<float>(other.ptr->begin(), other.ptr->end());
}

size_t Line::size() const
{
    return m_size;
}

Line::~Line()
{
    delete ptr;
}