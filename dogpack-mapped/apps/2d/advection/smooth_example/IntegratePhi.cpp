///////////////////////////////////////////////////////////////////////////////
// Function that integrates the source term exactly
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"

// This function is ONLY appropriate for solving f_t + v f_x = psi, since
// the characteristics for v can be solved exactly.
void IntegratePhi( double tnew, const dTensor2 &xpts, const dTensor1 &tpts, dTensor3 &psi)
{

    const int numpts  = xpts.getsize(1);
    const int meqn    = psi.getsize(2);

    double a,x,v,t,E,qv;

    for( int i=1; i <= numpts; i++ )
    for( int me=1; me <= meqn; me++ )
    {

        a = xpts.get(i,1);  // x-coordinate of quadrature point //
        v = xpts.get(i,2);  // v-coordinate of quadrature point //

        for( int n=1; n <= tpts.getsize(); n++ )
        {

            t   = tpts.get(n);
            x = a+v*(t-tnew); // value along the characteristic at time t

            E    = 0.25 * (t -   t*t ) * cos( pi * x );
            qv = ( -8.0*v ) * sin(2.0*pi*(x-pi*t)) * exp(-4.0*v*v);

//          qt  = 2.0*pi*cos(2.0*pi*(x-pi*t))*(-pi)*exp(-4.0*v*v);
//          vqx = 2.0*pi*cos(2.0*pi*(x-pi*t))*( v )*exp(-4.0*v*v);

            double ps =
            2.0*pi*cos(2.0*pi*(x-pi*t))*(v-pi)*exp(-4.0*v*v) + E*qv;
            psi.set(i, me, n, ps );

        }
    }

}
///////////////////////////////////////////////////////////////////////////////
