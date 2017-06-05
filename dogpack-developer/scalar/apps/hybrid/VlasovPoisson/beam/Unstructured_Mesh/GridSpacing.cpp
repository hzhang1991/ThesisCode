#include "meshdefs.h"

//  Grid spacing function: 
//
//     This functions sets a relative grid spacing of points
//
//     The input is a point (i.e., x and y coordinate), and the
//     output is a real positive number called "hdist". Large (small)
//     "hdist" means relatively large (small) grid spacing (relative to
//     the maximum and minimum values of GridSpacing.cpp).
//
//     For example, GridSpacing = 1   for all input x and y and
//                  GridSpacing = 55  for all input x and y 
//
//           will produce the same nearly uniform mesh since 
//           the relative variation of GridSpacing in both
//           examples is zero.
//
double GridSpacing(point pt)
{
  double Min(double,double);
  double xin = pt.x;
  double yin = pt.y;
  double rad = sqrt(pow(xin-0.5,2)+pow(yin-0.5,2));
  double offrad = sqrt(pow(xin-0.25,2)+pow(yin-0.5,2));

  double hdist = 0.10;//0.15+offrad;//1.0;
  
  return hdist;
}
