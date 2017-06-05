#include "dogdefs.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(
    const dTensor2* vel_vec,
    const dTensor2& xpts, 
    dTensor2& auxvals)
{
    const int numpts=xpts.getsize(1);

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        // u:  1-component of the advection velocity
        auxvals.set(i,1,  2.0*pi*(y-0.5) );

        // v:  2-component of the advection velocity
        auxvals.set(i,2, -2.0*pi*(x-0.5) );

    }

}
