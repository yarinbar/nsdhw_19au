#include <pybind11/pybind11.h>

#include <cmath>
#include <tuple>


double calculate_angle(
    std::tuple<double, double> const & vec1
  , std::tuple<double, double> const & vec2
)
{

    const double vec1_length = 
        sqrt(pow(std::get<0>(vec1), 2.0)+pow(std::get<1>(vec1), 2.0));
    const double vec2_length = 
        sqrt(pow(std::get<0>(vec2), 2.0)+pow(std::get<1>(vec2), 2.0));
    if(vec1_length == 0 || vec2_length == 0) {
        // if length == 0, return -100, invalid input
        return -100.0;
    }
    const double dot_product = std::get<0>(vec1)*std::get<0>(vec2) +
        std::get<1>(vec1)*std::get<1>(vec2);
    
    return acos(dot_product/(vec1_length*vec2_length));

}

PYBIND11_MODULE(hw2q2, mod)
{

    mod.doc() = 
        "nsd hw2 q2 - calculate angles between 2-dimensional vectors";

    mod.def("calculate_angle", &calculate_angle, 
        "calculate angles between 2-dimensional vectors");

}