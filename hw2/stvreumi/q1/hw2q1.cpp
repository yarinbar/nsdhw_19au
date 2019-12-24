#include <iostream>
#include <vector>
#include <tuple>

class Line
{
public:
    Line() = default;
    Line(Line const & l){
        coord_vec = l.coord_vec;
    } // copy constructor
    Line(size_t size) { coord_vec.resize(size); }
    size_t size() const { return coord_vec.size(); }
    ~Line() {};
    float   x(size_t it) const { return std::get<0>(coord_vec.at(it)); }
    float & x(size_t it)       { return std::get<0>(coord_vec.at(it)); }
    float   y(size_t it) const { return std::get<1>(coord_vec.at(it)); }
    float & y(size_t it)       { return std::get<1>(coord_vec.at(it)); }
private:
    // Member data.
    std::vector< std::tuple<float,float> > coord_vec;
}; /* end class Line */

int main(int, char **)
{
    Line line(3);
    line.x(0) = 0; line.y(0) = 1;
    line.x(1) = 1; line.y(1) = 3;
    line.x(2) = 2; line.y(2) = 5;

    Line line2(line);
    line2.x(0) = 9;

    std::cout << "line: number of points = " << line.size() << std::endl;
    for (size_t it=0; it<line.size(); ++it)
    {
        std::cout << "point " << it << ":"
                  << " x = " << line.x(it)
                  << " y = " << line.y(it) << std::endl;
    }

    std::cout << "line2: number of points = " << line.size() << std::endl;
    for (size_t it=0; it<line.size(); ++it)
    {
        std::cout << "point " << it << ":"
                  << " x = " << line2.x(it)
                  << " y = " << line2.y(it) << std::endl;
    }

    return 0;
}