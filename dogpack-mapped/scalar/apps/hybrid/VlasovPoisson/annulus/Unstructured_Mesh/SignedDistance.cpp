#include "meshdefs.h"

//  Signed distance function: 
//
//      SignedDistance(x,y) < 0 inside the region
//      SignedDistance(x,y) = 0 on the boundary
//      SignedDistance(x,y) > 0 outside of the region
//
double SignedDistance(point pt)
{
    double xin = pt.x;
    double yin = pt.y;
    double Min(double,double);
    double Max(double,double);

    double rad = sqrt(xin*xin+yin*yin);
    double d1 = 1.0 - rad;      
    double d2 = rad - 10.0;       
    return Max(d1,d2);

}
