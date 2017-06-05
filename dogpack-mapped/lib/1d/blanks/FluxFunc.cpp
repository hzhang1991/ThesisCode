#include "tensors.h"

// *REQUIRED* 
//
// Flux function f in the hyperbolic balance law:
//
//            q_t + f_x = psi,
//
// Input:
//
//    xpts( 1:numpts )           - The x-coordinates for a list of points
//    Q   ( 1:numpts, 1:meqn )   - The solution at each of these points
//    Aux ( 1:numpts, 1:maux )   - The auxilary function at each point.
//
// Output:
//
//    flux( 1:numpts, 1:meqn )  - The flux function f(q) defined at each point
//
// See also: AuxFunc.
void FluxFunc(const dTensor1& xpts, const dTensor2& Q,const dTensor2& Aux, dTensor2& flux)
{

    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {

        // Flux function f(q) in q_t + f(q)_x = psi.
        double x = xpts.get(i);

    }

}
