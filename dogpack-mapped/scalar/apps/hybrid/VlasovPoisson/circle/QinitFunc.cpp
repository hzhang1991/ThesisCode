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

    // f(t=0,x,v) = rho_0 * exp(-0.5*v^2 / temp ) / sqrt(2*pi*temp) 
    const double temp = 5.526350206e-4;
    const double rho0 = 1.0;

    const double one_over_pi = 1.0/pi;
    for (int i=1; i<=numpts; i++)
    for (int m=1; m<=meqn; m++)
    {

        const double x  = xpts.get(i,1);
        const double y  = xpts.get(i,2);

        const double vx = vel_vec->get(m,1);
        const double vy = vel_vec->get(m,2);

        // -- Constant initial density -- //
        qvals.set(i,m, rho0 * exp( -0.5*(vx*vx+vy*vy)/temp ) / ( 2.0*pi*temp ) );

//      // multiply initial conditions by a normalized maxwellian:
//      qvals.set(i,m, qvals.get(i,m) * exp(-vx*vx-vy*vy) * one_over_pi );

    }
}
