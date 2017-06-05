#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// The initial conditions for this problem define a (non-physical) Gaussian
// pulse of the form:
//
//      f(x,v, vx, vy) = (1-(x^2+y^2)) * G( vx, sigma1, v1c ) * G( vy, sigma2, v2c ),
//
// where
//
//      G( v, sigma, vc ) := 1/(sqrt(2*pi)) * exp( -(v-xc)^2 / (2*sigma^2) ).
//
// Note that \int G(v) = 1, \int v*G(v) = xc, and 
//           \int v^2 G(v) = sigma^2 + vc^2
//
void QinitFunc(
    const dTensor2* vel_vec,
    const dTensor2& xpts, 
    dTensor2& qvals)
{

    const int numpts  = qvals.getsize(1);
    const int meqn    = qvals.getsize(2);

    const double one_over_pi = 1.0/pi;

    const double sigma1 = 1.0;
    const double sigma2 = 2.0;

    const double v1c =  1.2;
    const double v2c = -1.1;

    for (int i=1; i<=numpts; i++)
    for (int m=1; m<=meqn; m++)
    {
        double x  = xpts.get(i,1);
        double y  = xpts.get(i,2);

        double r2 = pow(x,2)+pow(y,2);
        double r  = sqrt(r2);

        double vx = vel_vec->get(m,1) - v1c;
        double vy = vel_vec->get(m,2) - v2c;

        qvals.set(i,m, (1.0-r2)*exp( 
            -(vx*vx)/(2.0*sigma1*sigma1) 
            -(vy*vy)/(2.0*sigma2*sigma2) ) / (2.0*pi*sigma1*sigma2) );

    }
}
