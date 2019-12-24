#include<pybind11/pybind11.h>
#include<cmath>
using namespace std;
double calculateAngle(pair<double,double> vec1,pair<double,double> vec2){
	double rawAngle = atan2(vec2.second, vec2.first) - atan2(vec1.second, vec1.first);
	return (rawAngle<0)?rawAngle+2*M_PI:rawAngle;
}
PYBIND11_MODULE(angle, m) {
    m.doc() = "pybind11 calculate Angle from two vectors";
    m.def("calculateAngle", &calculateAngle, "calaulate Angle.");
}
