#include "dogdefs.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
void AuxFunc(
    const dTensor2* vel_vec,
    const dTensor2& xpts, dTensor2& auxvals)
{

    const int numpts  = xpts.getsize(1);
    const int meqn    = xpts.getsize(2);

    for (int i=1; i<=numpts; i++)
    for (int m=1; m<=meqn;   m++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double vx = vel_vec->get(m,1);
        double vy = vel_vec->get(m,2);

        // u:  1-component of the advection velocity
        auxvals.set(i, 1, vx );

        // v:  2-component of the advection velocity
        auxvals.set(i, 2, vy  );

    }

}
