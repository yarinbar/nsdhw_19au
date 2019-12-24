
#include "line.hpp"

Line::Line(size_t size)
{
    m_points = new std::vector<point_t>(size);
}

Line::~Line()
{
    delete m_points;
}

size_t Line::size() const
{
    return m_points->size();
}

Line::Line(Line const &other)
{
    if(other.m_points)
    {
        if(m_points)
            delete m_points;

        m_points = new std::vector<point_t>(other.size());
        std::copy(other.m_points->begin(), other.m_points->end(), m_points->begin());
    }
    else
    {
        if(m_points)
        {
            delete m_points;
            m_points = nullptr;
        }
    }
}

Line::Line(Line &&other)
{
    if(m_points)
    {
        delete m_points;
        m_points = nullptr;
    }

    m_points = other.m_points;
    other.m_points = nullptr;
}

Line & Line::operator=(Line const &other)
{
    if(this == &other)
        return *this;

    if(other.m_points)
    {
        if(m_points)
            delete m_points;
        m_points = new std::vector<point_t>(other.size());
        std::copy(other.m_points->begin(), other.m_points->end(), m_points->begin());
    }
    else
    {
        if(m_points)
        {
            delete m_points;
            m_points = nullptr;
        }
    }
    return *this;
}

Line & Line::operator=(Line &&other)
{
    if(this == &other)
        return *this;
    
    if(m_points)
    {
        delete m_points;
        m_points = nullptr;
    }
    m_points = other.m_points;
    other.m_points = nullptr;
    return *this;
}

float & Line::x(size_t it) const
{
    return m_points->at(it).first;
}

float & Line::x(size_t it)
{
    return m_points->at(it).first;
}

float & Line::y(size_t it) const
{
    return m_points->at(it).second;
}

float & Line::y(size_t it)
{
    return m_points->at(it).second;
}