#include "dogdefs.h"
#include <cmath>
#include "dog_math.h"
#include "EulerParams.h"
#include <iostream>
using namespace std;

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     2D Euler equations
//
void FluxFunc(const dTensor2& xpts,
        const dTensor2& Q,
        const dTensor2& Aux,
        dTensor3& flux)
{
    const int numpts=xpts.getsize(1);

    // Gas constant
    double const gamma = 1.4;//eulerParams.gamma;

    for (int i=1; i<=numpts; i++)
    {
        // Variables
        double rho    = Q.get(i,1);
        double u1     = Q.get(i,2)/rho;
        double u2     = Q.get(i,3)/rho;
        double u3     = Q.get(i,4)/rho;
        double energy = Q.get(i,5);
        double press  = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

        // 1-component of flux function
        flux.set(i,1,1, rho*u1 );
        flux.set(i,2,1, rho*u1*u1 + press );
        flux.set(i,3,1, rho*u1*u2 );
        flux.set(i,4,1, rho*u1*u3 );
        flux.set(i,5,1, u1*(energy+press) );

        // 2-component of flux function
        flux.set(i,1,2, rho*u2 );
        flux.set(i,2,2, rho*u2*u1 );
        flux.set(i,3,2, rho*u2*u2 + press );
        flux.set(i,4,2, rho*u2*u3 );
        flux.set(i,5,2, u2*(energy+press) );
    }
}
