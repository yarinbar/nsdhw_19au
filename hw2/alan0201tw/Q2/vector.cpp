#include <pybind11/pybind11.h>

#include <cmath>
#include <tuple>

double compute_angle(
    std::tuple<double, double> const & vec0,
    std::tuple<double, double> const & vec1
    )
{
    const double x0 = std::get<0>(vec0);
    const double y0 = std::get<1>(vec0);
    const double x1 = std::get<0>(vec1);
    const double y1 = std::get<1>(vec1);

    const double dot = x0*x1 + y0*y1;
    const double det = x0*y1 - y0*x1; 

    // Reference : http://www.cplusplus.com/reference/cmath/atan2/
    //
    // If both arguments passed are zero, a domain error occurs.
    //
    // If a domain error occurs:
    // - And math_errhandling has MATH_ERRNO set: the global variable errno is set to EDOM.
    // - And math_errhandling has MATH_ERREXCEPT set: FE_INVALID is raised.
    //
    // Since we only use atan2 from the C++ standard math library, instead
    // of doing some proper error handling, I just return a nan value when
    // the arguments are not proper for atan2, and leave it to the python
    // side to deal with it.

    if(det == 0.0 && dot == 0.0)
        return nan("");

    // Principal arc tangent of y/x, in the interval [-pi,+pi] radians.
    // One radian is equivalent to 180/PI degrees.
    // 
    // When returning, the receiver should be expecting a positive angle.
    // In order to deal with negative radians ( range from -pi to 0 ).
    // We should map [-pi, 0) to [-pi + 2pi, 0 + 2pi) = [pi, 2pi)
    // which leads to the bigger angle between the two vectors.
    // So return 2*pi - value = [2*pi - pi, 2*pi- 2*pi)
    // = [ pi, 0 ). and the result is ( 2pi - (2pi + result) )
    // which is negative result, when result is negative.
    //
    // In conclusion, a simple abs will do the trick.

    return std::abs(std::atan2(det, dot));
}

PYBIND11_MODULE(vector, mod)
{
    mod.doc() = "example C extension module";
    mod.def("compute_angle", &compute_angle, 
        "calculates the angle (in radians) between two vectors in the "
        "2-dimensional Cartesian coordinate system");
}

/* vim: set et ts=4 sw=4: */
