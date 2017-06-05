#ifndef NDIMS
#define NDIMS 0
#endif
#include <dogdefs.h>
// #include <stdlib.h>
// #include <stdio.h>
// #include "tensors.h"   // includes dimdefs.h
#include <cmath>

int test_constants()
{
  const double Pi = M_PI; // 4.0*atan(1.0);
  // This is not expensive since sqrt is implemented in hardware
  // but it requires the inclusion of <cmath>
  const double Sq2 = sqrt(2.0);
  const double Sq3 = sqrt(3.0);
  const double Sq5 = sqrt(5.0);
  const double Sq7 = sqrt(7.0);
  const double Sq10 = sqrt(10.0);
  const double Sq13 = sqrt(13.0);
  const double Sq19 = sqrt(19.0);
  const double Sq23 = sqrt(23.0);
  const double Sq71 = sqrt(71.0);

  // const double pi   = 3.141592653589793;
  // const double sq2  = 1.4142135623730951;
  // const double sq3  = 1.7320508075688772;
  // const double sq5  = 2.2360679774997898;
  // const double sq7  = 2.6457513110645907;
  // const double sq10 = 3.1622776601683795;
  // const double sq13 = 3.6055512754639891;
  // const double sq19 = 4.3588989435406740;
  // const double sq23 = 4.7958315233127191;
  // const double sq71 = 8.4261497731763590;

  printf("Pi   = %24.16f\n",Pi);
  printf("Sq2  = %24.16f\n",Sq2);
  printf("Sq3  = %24.16f\n",Sq3 );
  printf("Sq5  = %24.16f\n",Sq5 );
  printf("Sq7  = %24.16f\n",Sq7 );
  printf("Sq10 = %24.16f\n",Sq10);
  printf("Sq13 = %24.16f\n",Sq13);
  printf("Sq19 = %24.16f\n",Sq19);
  printf("Sq23 = %24.16f\n",Sq23);
  printf("Sq71 = %24.16f\n",Sq71);

  // verify exact agreement with computed values
  assert(Pi   == pi  );
  assert(Sq2  == sq2 );
  assert(Sq3  == sq3 );
  assert(Sq5  == sq5 );
  assert(Sq7  == sq7 );
  assert(Sq10 == sq10);
  assert(Sq13 == sq13);
  assert(Sq19 == sq19);
  assert(Sq23 == sq23);
  assert(Sq71 == sq71);
}

int main()
{
  test_constants();
}
