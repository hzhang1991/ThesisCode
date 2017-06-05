#include "tensors.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// The expected format is Dflux.get(:, i, j) = \partial f_i, \partial q_j.
//
//     Simple advection equation, f'(q) = u
//     Burger's equation, f'(q) = q
//     Acoustics equation, f'(q) = [0 1; 1 0]
//
void DFluxFunc(const dTensor1& xpts, 
	       const dTensor2& Q,
	       const dTensor2& Aux,
	       dTensor3& Dflux)
{

    const int numpts=xpts.getsize();
    for(int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        Dflux.set(i,1,1, 0. );
    }

}
