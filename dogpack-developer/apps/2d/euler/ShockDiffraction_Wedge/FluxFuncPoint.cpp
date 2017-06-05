#include "dogdefs.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     2D Euler equations
//
void FluxFuncPoint(const dTensor2& xpts,
        const dTensor1& Q,
        const dTensor1& Aux,
        dTensor1& xflux,dTensor1& yflux)
{
    const int numpts=xpts.getsize(1);

    // Gas constant
    double const gamma = eulerParams.gamma;

    {
        // Variables
        double rho    = Q.get(1);
        double u1     = Q.get(2)/rho;
        double u2     = Q.get(3)/rho;
        double u3     = Q.get(4)/rho;
        double energy = Q.get(5);
        double press  = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

        // 1-component of flux function
        xflux.set(1, rho*u1 );
        xflux.set(2, rho*u1*u1 + press );
        xflux.set(3, rho*u1*u2 );
        xflux.set(4, rho*u1*u3 );
        xflux.set(5, u1*(energy+press) );

        // 2-component of flux function
        yflux.set(1, rho*u2 );
        yflux.set(2, rho*u2*u1 );
        yflux.set(3, rho*u2*u2 + press );
        yflux.set(4, rho*u2*u3 );
        yflux.set(5, u2*(energy+press) );
    }
}
