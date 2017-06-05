#include "meshdefs.h"
#include "dog_math.h"
#include "constants.h"  // for access to Pi

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

    double d1 =  2.0*pi  - yin;
    double d2 =  0.0     + xin;
    double d3 =  4.0*pi  - xin;
    double d4 =  2.0*pi  + yin;

    double dist = -Min(Min(d1,d2),Min(d3,d4));

    return dist;
}
