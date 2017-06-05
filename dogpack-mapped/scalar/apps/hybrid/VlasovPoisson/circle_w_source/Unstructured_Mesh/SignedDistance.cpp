#include "meshdefs.h"


//  Signed distance function: 
//
//      SignedDistance(x,y) < 0 inside the region
//      SignedDistance(x,y) = 0 on the boundary
//      SignedDistance(x,y) > 0 outside of the region
//
double SignedDistance(point pt)
{
  double Min(double,double);
  double xin = pt.x;
  double yin = pt.y;
  //double rad = sqrt(pow(xin-0.5,2)+pow(yin-0.5,2));
  double rad = sqrt(pow(xin,2)+pow(yin,2));
   
  //double dist = rad - 0.5e0;
  double dist = rad - 1.0;
    
  return dist;
}
