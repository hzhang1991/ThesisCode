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
    
  double d1 = 1.1e0 - yin;
  double d2 = 1.1e0 + xin;
  double d3 = 1.1e0 - xin;
  double d4 = 1.1e0 + yin;

  double dist2 = 1.1e0*1.1e0 - xin*xin-yin*yin;
  
  double dist = -Min(Min(d1,d2),Min(d3,d4));
    
  return -dist2;
}
