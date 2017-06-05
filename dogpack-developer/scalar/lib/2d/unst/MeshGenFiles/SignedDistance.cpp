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
    
  double d1 = 1.0e0 - yin;
  double d2 = 0.0e0 + xin;
  double d3 = 1.0e0 - xin;
  double d4 = 0.0e0 + yin;
  
  double dist = -Min(Min(d1,d2),Min(d3,d4));
    
  return dist;
}
