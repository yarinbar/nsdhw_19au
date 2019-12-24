#include <vector>

class Line
{
public:
    Line();
    Line(Line const & );
    //Line(Line       &&);
    //Line & operator=(Line const & );
    //Line & operator=(Line       &&);
    Line(size_t size);
    ~Line();
    size_t size() const;
    float & x(size_t it) const;
    float & x(size_t it);
    float & y(size_t it) const;
    float & y(size_t it);
private:
    void check_range(size_t it) const;
    size_t m_size; // number of points.
    std::vector<float> *ptr;
};