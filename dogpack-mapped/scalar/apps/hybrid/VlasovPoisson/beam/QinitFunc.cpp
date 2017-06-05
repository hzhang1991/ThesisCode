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

    // Beam parameters:
    const double a      = 1.0;          // Parameter used in definition of beam
    const double vth2   = 1.0;          // Thermal velocity
    const double R2     = 0.25;         // Square of R  ( R^2 = a^2 / 4 )
    const double eta    = 0.25;         // Tune depression
    const double omega0 = 8.0;          // Derived by solving v_th = a*eta*omega0 / 2
    const double omega  = eta*omega0;   // definition of tune depression
    const double n0     = 2.0*pi*a*a*(1.0-eta*eta)*omega0*omega0;

    const double one_over_pi = 1.0/pi;
    for (int i=1; i<=numpts; i++)
    for (int m=1; m<=meqn; m++)
    {

        const double x  = xpts.get(i,1);
        const double y  = xpts.get(i,2);
        const double r2 = x*x+y*y;

        const double vx = vel_vec->get(m,1);
        const double vy = vel_vec->get(m,2);

        // -- Beam parameters -- //
        double tmp = n0 / ( 4.0*pi*pi*vth2*R2 ) * 
            exp( -(vx*vx+vy*vy)/(2.0*vth2) ) *
            exp( -( x* x+ y* y)/(2.0*  R2) );
        qvals.set(i,m, tmp );

    }
}
