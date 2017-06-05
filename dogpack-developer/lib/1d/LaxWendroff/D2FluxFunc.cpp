#include "tensors.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
// In 1d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k) = \partial^2 f_i / \partial q_j \partial q_k
//
//     Dummy Function Call if not using LaxWendroff
//
void D2FluxFunc(const dTensor1& xpts, 
		const dTensor2& Q,
		const dTensor2& Aux,
		dTensor4& D2flux)
{

    const int numpts=xpts.getsize();
    D2flux.setall(0.);

}
