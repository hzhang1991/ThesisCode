#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
  const int numpts=xpts.getsize(1);

  for (int i=1; i<=numpts; i++)
    {
      const double x = xpts.get(i,1);
      const double y = xpts.get(i,2);
      
      const double r2 = abs(x-0.2);//pow(x-0.20,2)+pow(y-0.5,2);
      const double r  = sqrt(r2);

       if( x > (0.05-1.0e-15) && x < (0.25-1.0e-15) )// && y > 0.25 && y<0.75)
        { qvals.set(i,1,pow(1.0/2.0*(1.0+cos(pi*(x-0.15)*10.0)),6));}
      else
        {  qvals.set(i,1, 0.0 ); }

    }
}
