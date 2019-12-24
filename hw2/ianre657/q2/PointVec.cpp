#include <iostream>
#include <cmath>
#include "PointVec.hpp"


double angle_btw( const PointVec& p1, const PointVec& p2 ){
    
    double len1 = sqrt( double(p1.m_x*p1.m_x + p1.m_y*p1.m_y) );
    double len2 = sqrt( double(p2.m_x*p2.m_x + p2.m_y*p2.m_y) );

    double ux = p1.m_x/len1, uy = p1.m_y/len1;
    double vx = p2.m_x/len2, vy = p2.m_y/len2;

    double dotval = ux*vx + uy*vy;
    if (dotval <= -1.0){
        return M_PI;
    }else if(dotval >= 1.0){
        return .0;
    }
    return acos(dotval);
}