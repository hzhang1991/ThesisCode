#include "dogdefs.h"
#include "dog_math.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     2D Acoustics
//
void FluxFunc(const dTensor2& xpts,
	      const dTensor2& Q,
	      const dTensor2& Aux,
	      dTensor3& flux)
{
  const int numpts=xpts.getsize(1);
  
  for (int i=1; i<=numpts; i++)
    {
      // Variables
      double p = Q.get(i,1);
      double u = Q.get(i,2);
      double v = Q.get(i,3);

      double cs  = Aux.get(i,1);
      double cs2 = cs*cs;
      double rho = Aux.get(i,2);
      
      // 1-component of flux function
      flux.set(i,1,1, cs2*rho*u );
      flux.set(i,2,1, p/rho     );
      flux.set(i,3,1, 0.0       );
      
      // 2-component of flux function
      flux.set(i,1,2, cs2*rho*v );
      flux.set(i,2,2, 0.0       );
      flux.set(i,3,2, p/rho     );
    }
  
}
