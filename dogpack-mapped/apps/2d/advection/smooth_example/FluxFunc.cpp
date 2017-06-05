#include "dogdefs.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Simple advection equation
//
void FluxFunc(const dTensor2& xpts, const dTensor2& Q, 
        const dTensor2& Aux, dTensor3& flux)
{
    int i;
    int numpts=xpts.getsize(1);
    double x,y,qc,u,v;

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i,1);
        y = xpts.get(i,2);

        // Variables
        qc = Q.get(i,1);
        u  = Aux.get(i,1);
        v  = Aux.get(i,2);

        // Flux function
        flux.set(i,1,1, u*qc );
        flux.set(i,1,2, v*qc );
    }

}
