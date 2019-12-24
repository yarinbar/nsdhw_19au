#pragma once

#include <vector>
#include <memory>

class Line
{
    typedef std::shared_ptr< std::vector<float> > data_ref;

private:
    size_t m_size = 0; // number of points, should match the size of m_points
    data_ref m_points;

    void bound_check(size_t it) const;

public:

    // constructors
    Line() = default;
    Line(size_t size);

    // destructor
    ~Line();

    // copy / move operators
    Line(Line const & );
    Line(Line       &&);
    Line & operator=(Line const & );
    Line & operator=(Line       &&);

    // inline methods
    inline size_t size() const { return m_size; }

    // methods
    float & x(size_t it) const;
    float & x(size_t it);
    float & y(size_t it) const;
    float & y(size_t it);
};