#include <iostream>
#include <vector>

class Line
{
    std::vector <size_t> data1, data2;
public :
    Line()
    {
        data1.clear();
        data2.clear();
    }
    Line(size_t n)
    {
        data1.resize(n);
        data2.resize(n);
    }
    Line(Line & line)
    {
        data1 = line.data1;
        data2 = line.data2;
    }
    size_t & x(size_t n) {return data1[n];}
    size_t & y(size_t n) {return data2[n];}
    size_t size() {return data1.size();}
};

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