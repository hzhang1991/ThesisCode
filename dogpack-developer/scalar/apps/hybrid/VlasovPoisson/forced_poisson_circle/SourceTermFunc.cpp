#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
//
// Exact solution: f = A(t) * (1-x^2-y^2) * exp( -vx^2 ) / sqrt(pi)
//
// Where, A(t) = ( 0.75 + 0.25*cos(2*pi*t) )
//
void SourceTermFunc(
    const double t, 
    const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& qvals, 
    const dTensor2& auxvals, dTensor2& source)
{

    const int numpts    = source.getsize(1);
    const int meqn      = source.getsize(2);

    // Time dependent part:
    const double A  = (0.75 + 0.25 * cos(2.0*pi*t) );
    const double At = -0.5*pi*sin(2.0*pi*t);

    const double one_over_sqpi = 1.0/sqrt(pi);
    for (int i=1; i<=numpts; i++)
    {

        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        double r2 = pow(x,2) + pow(y,2);
        double r  = sqrt( r2 );

        for (int m=1; m<=meqn; m++)
        {

            const double vx = vel_vec->get(m,1);
            const double vy = vel_vec->get(m,2);

            // spatial and velocity dependant part:
            double qt  = At * ( 1.0 - r2 ) * exp(-vx*vx ) * one_over_sqpi;
            double qx  = A  * ( -2.0*x   ) * exp(-vx*vx ) * one_over_sqpi;
            double qy  = A  * ( -2.0*y   ) * exp(-vx*vx ) * one_over_sqpi;
            double qvx = A  * ( 1.0 - r2 ) * ( -2.0*vx*exp(-vx*vx) ) * one_over_sqpi;
            double qvy = 0.;

            // Electric field:
            double E1 = A*( -0.5*x + x*x*x/3. );
            E1 = 0.;

            double tmp = qt + vx*qx + vy*qy - E1*qvx;
            source.set(i,m, tmp );

        }
    }

}
