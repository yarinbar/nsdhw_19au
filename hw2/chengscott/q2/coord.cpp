#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <cmath>
#include <numeric>
#include <vector>

namespace py = pybind11;

template <class ReturnType, class VectorType>
ReturnType angle(const VectorType &lhs, const VectorType &rhs)
{
    auto inner_product_ = [](const VectorType &lhs, const VectorType &rhs) -> ReturnType {
        return std::inner_product(std::begin(lhs), std::end(lhs), std::begin(rhs), 0);
    };
    ReturnType l2 = inner_product_(lhs, lhs), r2 = inner_product_(rhs, rhs);
    return std::acos(inner_product_(lhs, rhs) / std::sqrt(l2) / std::sqrt(r2));
}

PYBIND11_MODULE(coord, m)
{
    m.doc() = "Vector Arithmetics";
    m.def("angle", &angle<float, std::vector<float>>, "angle between two vectors");
}
