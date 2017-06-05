#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(
    const dTensor2* vel_vec,
    const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts  = qvals.getsize(1);
    const int meqn    = qvals.getsize(2);

    for (int i=1; i<=numpts; i++)
    for (int m=1; m<=meqn; m++)
    {

        const double x  = xpts.get(i,1);
        const double y  = xpts.get(i,2);

        const double vx = vel_vec->get(m,1);
        const double vy = vel_vec->get(m,2);

        double r  = sqrt( pow(x,2) + pow(y,2) );

        qvals.set(i,m, sin(2.0*pi*vx)*sin(2.0*pi*vy) );

    }
}
