#include <cmath>
#include "dogdefs.h"
#include "VlasovParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(
    const dTensor2* vel_vec,
    const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts  = qvals.getsize(1);
    const int meqn    = qvals.getsize(2);

    const double alpha = 0.01;
    const double kx    = 0.5;
    const double ky    = 0.5;

    for (int i=1; i<=numpts; i++)
    for( int m=1; m<=meqn;   m++)
    {

        // configuration space coordinates
        const double x  = xpts.get(i,1);
        const double y  = xpts.get(i,2);

        // velocity space coordinates
        const double vx = vel_vec->get(m,1);
        const double vy = vel_vec->get(m,2);

        double tmp = 0.5/pi * exp( -0.5*(vx*vx+vy*vy) ) * ( 1.0 + alpha*cos(kx*x)*cos(ky*y) );
        qvals.set(i,m, tmp); 

    }

}
