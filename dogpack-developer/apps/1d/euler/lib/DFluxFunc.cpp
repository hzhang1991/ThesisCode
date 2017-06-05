#include <cmath>
#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// Euler Equations
//
void DFluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor3& Dflux)
{

    const int numpts=xpts.getsize();

    // Gas constant
    const double gamma = eulerParams.gamma;

    Dflux.setall(0.);
    for (int i=1; i<=numpts; i++)
    {

        // Not used: u2 and u3 ... 
        double q1 = Q.get(i,1);
        double q2 = Q.get(i,2);
        double q3 = Q.get(i,5);

        // pd{ f_1 }{ q_j }
        Dflux.set(i,1,1, 0.0e0);
        Dflux.set(i,1,2, 1.0e0);
        Dflux.set(i,1,5, 0.0e0);

        // pd{ f_2 }{ q_j }
        Dflux.set(i,2,1, 0.5*(gamma-3.0)*pow((q2/q1),2) );
        Dflux.set(i,2,2, -(q2/q1)*(gamma-3.0)           );
        Dflux.set(i,2,5, gamma-1.0                      );

        // pd{ f_3 }{ q_j }
        Dflux.set(i,5,1, 
            -q2*( gamma*(q1*q3-pow(q2,2)) + pow(q2,2))*pow(q1,-3) );
        Dflux.set(i,5,2, 
            0.5*(2.0*gamma*q1*q3+3.0*pow(q2,2)*(1.0-gamma))*pow(q1,-2) );
        Dflux.set(i,5,5, q2*gamma/q1 );

    }

}
