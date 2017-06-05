#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc(
    const double t, 
    const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& qvals, 
    const dTensor2& auxvals, dTensor2& source)
{

    const int numpts = xpts.getsize(1);
    const int meqn  = source.getsize(2);

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        for (int m=1; m<=meqn; m++)
        {
            const double vx = vel_vec->get(m,1);
            const double vy = vel_vec->get(m,2);

            double tmp = -qvals.get(i,m);
            source.set(i,m, tmp );

        }
    }

}
