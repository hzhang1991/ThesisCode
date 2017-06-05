#include "tensors.h"

// *REQUIRED* 
//
// This is a user-required routine that defines the initial conditions for the
// problem.
//
// Each application is REQUIRED to define one of these.
//
// Input:
//
//    xpts( 1:numpts )           - The x-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: AuxFunc.
void QinitFunc(const dTensor1& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {

        // Spatial location.
        double x = xpts.get(i);

    }

}
