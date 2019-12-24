#include "PointVec.hpp"
#include <pybind11/pybind11.h>
namespace py = pybind11;

PYBIND11_MODULE(PointVec, m){
    m.doc() = "PointVec library";
    py::class_<PointVec> PointVec(m, "PointVec");
    
    PointVec.def(py::init<>());
    PointVec.def(py::init<double, double>());
    PointVec.def_readwrite("m_x", &PointVec::m_x)
        .def_readwrite("m_y", &PointVec::m_y);

    m.def("angle_btw",
        &angle_btw,
        "calculate angle between two Point Vectors"
    );
}