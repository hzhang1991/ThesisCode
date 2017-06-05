#include "meshdefs.h"
#include "dog_math.h"

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
    
  double d1 = 11.0e0 - yin;
  double d2 = 0.0e0 + xin;
  double d3 = 13.0e0 - xin;
  double d4 = 0.0e0 + yin;

  double d11 = 6.0e0 - yin;
  double d12 = 0.0e0 + xin;
  double d13 = 3.4e0/6.0e0*yin - xin;
  double d14 = 0.0e0 + yin;
  
  double dist1 = -Min(Min(d1,d2),Min(d3,d4));
  double dist2 = Min(Min(d11,d12),Min(d13,d14));
  double dist  = Max(dist1,dist2);  
  return dist;
}
