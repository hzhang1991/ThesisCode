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
    
  double d1 = 1.0e0 - yin;
  double d2 = 0.0e0 + xin;
  double d3 = 3.0e0 - xin;
  double d4 = 0.0e0 + yin;

  double d11 = 0.2e0 - yin;
  double d12 = xin   - 0.6;
  double d13 = 3.0   - xin;
  double d14 = 0.0e0 + yin;
  
  double dist1 = -Min(Min(d1,d2),Min(d3,d4));
  double dist2 = Min(Min(d11,d12),Min(d13,d14));
  double dist  = Max(dist1,dist2);  
  return dist;
}
