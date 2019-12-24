#include<iostream>
#include "line.hpp"

int main(int, char **)
{
    Line line(3);
    line.x(0) = 0; line.y(0) = 1;
    line.x(1) = 1; line.y(1) = 3;
    line.x(2) = 2; line.y(2) = 5;

    Line line2(line);
    Line line4;
    line2.x(0) = 9;
    
    std::cout << "line: number of points = " << line.size() << std::endl;
    for (size_t it=0; it<line.size(); ++it)
    {
        std::cout << "point " << it << ":"
                  << " x = " << line.x(it)
                  << " y = " << line.y(it) << std::endl;
    }

    std::cout << "line2: number of points = " << line2.size() << std::endl;
    for (size_t it=0; it<line2.size(); ++it)
    {
        std::cout << "point " << it << ":"
                  << " x = " << line2.x(it)
                  << " y = " << line2.y(it) << std::endl;
    }
    Line line3(std::move(line));
    
    std::cout << "line3: number of points = " << line3.size() << std::endl;
    for (size_t it=0; it<line3.size(); ++it)
    {
        std::cout << "point " << it << ":"
                  << " x = " << line3.x(it)
                  << " y = " << line3.y(it) << std::endl;
    }

    return 0;
}