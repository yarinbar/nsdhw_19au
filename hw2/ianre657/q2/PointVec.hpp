#ifndef __POINTVEC_LIB__
#define __POINTVEC_LIB__

struct PointVec{
    PointVec() = default;
    PointVec(double x, double y): m_x(x), m_y(y){};
    double m_x, m_y;
};

double angle_btw( const PointVec&, const PointVec& );
#endif