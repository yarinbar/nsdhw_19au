#include "line.hpp"

#include <algorithm>

/**
 *   constructor and copy/move operators
 */

Line::Line(size_t size)
    : m_size(size), m_points(std::make_shared<std::vector<float>>())
{
    m_points->resize(m_size * 2);
}

Line::Line(Line const & other)
{
    // do init and field copying if the other Line is not null
    if(other.m_points != nullptr)
    {
        m_points = std::make_shared<std::vector<float>>(*(other.m_points));
        m_size = other.m_size;
    }
    else
    {
        m_points = nullptr;
        m_size = 0;   
    }
}
Line::Line(Line && other)
{
    std::swap(other.m_size, m_size);
    std::swap(other.m_points, m_points);
}

Line & Line::operator=(Line const & other)
{
    if(this == &other) { return *this; } // don't copy myself.

    if(other.m_points != nullptr)
    {
        m_points = std::make_shared<std::vector<float>>(*(other.m_points));
        m_size = other.m_size;
    }
    else
    {
        m_points = nullptr;
        m_size = 0;   
    }
    return *this;
}

Line & Line::operator=(Line && other)
{
    if (this == &other) { return *this; } // don't move to self.
    std::swap(other.m_size, m_size);
    std::swap(other.m_points, m_points);
    return *this;
}

/**
 *   destructor
 */

// since we used smart pointer here, 
// we do not need to manually release the memory.
Line::~Line() = default;

/**
 *   public methods
 */

float & Line::x(size_t it) const
{
    bound_check(it);
    return m_points->at(it * 2);
}
float & Line::x(size_t it)
{
    bound_check(it);
    return m_points->at(it * 2);
}
float & Line::y(size_t it) const
{
    bound_check(it);
    return m_points->at(it * 2 + 1);
}
float & Line::y(size_t it)
{
    bound_check(it);
    return m_points->at(it * 2 + 1);
}

/**
 *   private methods
 */

void Line::bound_check(size_t it) const
{
    if (it >= m_size)
    {
        throw std::out_of_range("Line index out of range");
    }
}
