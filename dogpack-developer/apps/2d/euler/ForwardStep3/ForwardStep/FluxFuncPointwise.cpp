#include "dogdefs.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     2D Euler equations
//
void FluxFuncPointwise(const dTensor1& nvec,
        const dTensor1& Q,
        const dTensor1& Aux,
        dTensor1& flux)
{

    // Gas constant
    double const gamma = eulerParams.gamma;

        // Variables
        double rho    = Q.get(1);
        double u1     = Q.get(2)/rho;
        double u2     = Q.get(3)/rho;
        double u3     = Q.get(4)/rho;
        double energy = Q.get(5);
        double press  = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));
        
        flux.set(1,0.0);
        flux.set(2,0.0);
        flux.set(3,0.0);
        flux.set(4,0.0);
        flux.set(5,0.0);  

        // 1-component of flux function
        flux.set(1, flux.get(1)+nvec.get(1)*rho*u1 );
        flux.set(2, flux.get(2)+nvec.get(1)*(rho*u1*u1 + press ));
        flux.set(3, flux.get(3)+nvec.get(1)* rho*u1*u2 );
        flux.set(4, flux.get(4)+nvec.get(1)*rho*u1*u3 );
        flux.set(5, flux.get(5)+nvec.get(1)*u1*(energy+press) );

        // 2-component of flux function
        flux.set(1, flux.get(1)+nvec.get(2)*rho*u2 );
        flux.set(2, flux.get(2)+nvec.get(2)*rho*u2*u1 );
        flux.set(3, flux.get(3)+nvec.get(2)*(rho*u2*u2 + press ));
        flux.set(4, flux.get(4)+nvec.get(2)*rho*u2*u3 );
        flux.set(5, flux.get(5)+nvec.get(2)*u2*(energy+press) );
}
