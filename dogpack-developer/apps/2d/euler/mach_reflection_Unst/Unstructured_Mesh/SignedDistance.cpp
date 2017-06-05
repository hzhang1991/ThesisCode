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
    
  double d1 = 2.0e0 - yin;
  double d2 = 0.0e0 + xin;
  double d3 = 3.0e0 - xin;
  double d4 = 0.0e0 + yin;

 //0.1666666667      0.0
 //3.2      1.75129581654


  double d12 = -xin+3.0;
  double d13 = -yin + 1.44337567297*(xin-0.5)/2.5;
  double d14 = 0.0e0 + yin;
  
  double dist1 = -Min(Min(d1,d2),Min(d3,d4));
  double dist2 = Min(d12,Min(d13,d14));
  double dist  = Max(dist1,dist2);  
  return dist;
}
