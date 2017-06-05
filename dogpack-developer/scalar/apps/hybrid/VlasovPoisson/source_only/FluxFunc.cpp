#include "dogdefs.h"
#include "DogParams.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Simple advection equation
//
void FluxFunc(
    const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& Q, 
    const dTensor2& Aux, dTensor3& flux)
{
    const int numpts = Q.getsize(1);
    const int meqn   = Q.getsize(2);

    for (int i=1; i<=numpts; i++)
    for (int m=1; m<=meqn;   m++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double qc = Q.get(i,m);

        // Flux function
        if( dogParams.get_maux() > 0 )
        {
            // Variables
            double u  = Aux.get(i,1);
            double v  = Aux.get(i,2);
            flux.set(i,m,1, 0. );
            flux.set(i,m,2, 0. );
        }
        else
        {
            flux.set(i, m, 1, 0. );
            flux.set(i, m, 2, 0. );
        }
    }

}
