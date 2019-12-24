#include <iostream>
#include <cmath>
#include "Vector.h"
#include <pybind11/pybind11.h>

Vector::Vector(float x, float y)
{
    m_x = x;
    m_y = y;
}

float Vector::len()
{
    return sqrt(pow(m_x, 2) + pow(m_y, 2));
}

float Vector::getX()
{
    return m_x;
}

float Vector::getY()
{
    return m_y;
}

float dot_product(Vector v1, Vector v2)
{
    return v1.getX() * v2.getX() + v1.getY() * v2.getY();
}

/* Law of cosines */
float angle(Vector v1, Vector v2)
{
    return acos(dot_product(v1, v2) / (v1.len() * v2.len()));
}

PYBIND11_MODULE(Vector, m) {
    m.doc() = "angle (in radians) between two vectors";
    py::class_<Vector> Vector(m, "Vector");
        .def(py::init<float, float>());
        .def_readwrite("m_x", &Vector::m_x);
        .def_readwrite("m_y", &Vector::m_y);
        .def("getX", &Vector::getX());
        .def("getY", &Vector::getY());
        .def("get length", &Vector::len());        
 
    m.def("dot_product", &dot_product, "vector dot product")
    m.def("angle", &angle, "angle between teo vectors")
}


/*
int main(int argc, char *argv[])
{
    Vector v1(1, 0);
    Vector v2(0, 1);
    float rad;
    rad = angle(v1, v2);
    std::cout << rad << std::endl;
    return 0;
}
*/