#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
  const int numpts=xpts.getsize(1);

  // Loop over grid points
  for (int i=1; i<=numpts; i++)
    {
      // Variables
      double x = xpts.get(i,1);
      double y = xpts.get(i,2);

      double p = 0.0;
      double u = 0.0;
      double v = 0.0;

      if (x>0.2-1.0e-14 && x<0.4-1.0e-14)
	{
	  //printf("here 2\n");
	  p = 0.0;
	}      
      
      qvals.set(i,1, p );
      qvals.set(i,2, u );
      qvals.set(i,3, v );
      
    }

}
