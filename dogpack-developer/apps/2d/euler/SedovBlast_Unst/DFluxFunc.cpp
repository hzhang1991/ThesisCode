#include <cmath>
#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// Jacobian of the Flux Function 
//
// Euler Equations
//
void DFluxFunc(const dTensor2& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux, 
        dTensor4& Dflux)
{

    const int numpts = xpts.getsize(1);
    const int meqn   = Q.getsize(2);

    // Gas constant
    const double gamma = 1.4;//eulerParams.gamma;

    Dflux.setall(0.);
    for (int i=1; i<=numpts; i++)
    {

        // Not used: u3 ... 
        double q1 = Q.get(i,1);
        double q2 = Q.get(i,2);
        double q3 = Q.get(i,3);
        double q4 = Q.get(i,4);
        double q5 = Q.get(i,5);

        //////////////////////
        // Jacobian, f'(q): //
        //////////////////////

        // pd{ f_1 }{ q_j }
        Dflux.set(i,1,2, 1, 1.0e0);

        // pd{ f_2 }{ q_j }
        Dflux.set(i,2,1, 1, 
            0.5*(gamma-3.0)*( pow((q2/q1),2) )
           +0.5*(gamma-1.0)*( pow((q3/q1),2) + 2.0*pow((q4/q1),3) ) );
        Dflux.set(i,2,2, 1, -(q2/q1)*(gamma-3.0)             );
        Dflux.set(i,2,3, 1, -(q3/q1)*(gamma-1.0)             );
        Dflux.set(i,2,4, 1, 1.5*( gamma-1.0 )*pow((q4/q1),2) );
        Dflux.set(i,2,5, 1, gamma-1.0                        );

        // pd{ f_3 }{ q_j }
        Dflux.set(i,3,1, 1, -q2*q3*pow(q1,-2) );
        Dflux.set(i,3,2, 1, q3 / q1 );
        Dflux.set(i,3,3, 1, q2 / q1 );

        // pd{ f_4 }{ q_j }
        Dflux.set(i,4,1, 1, -q2*q4*pow(q1,-2) );
        Dflux.set(i,4,2, 1, q4 / q1 );
        Dflux.set(i,4,4, 1, q2 / q1 );

        // pd{ f_5 }{ q_j }
        Dflux.set(i,5,1, 1, -q5*(gamma*q2*pow(q1,-2))
            +(gamma-1.0)*( pow(q2/q1,3) + (q2/q1)*pow(q3/q1,2) )
            +1.5*(gamma-1.0)*( (q2/q1)*pow(q4/q1,3) )
        );
        Dflux.set(i,5,2, 1, gamma*q5/q1
            + 1.5*(1.0-gamma)*pow(q2/q1,2)
            + 0.5*(1.0-gamma)*( pow(q3/q1,2) + pow(q4/q1,3 ) ) );
        Dflux.set(i,5,3, 1, q2*q3*pow(q1,-2)*(1.0-gamma)       );
        Dflux.set(i,5,4, 1, 1.5*(1.0-gamma)*q2*pow(q4/q1,3)*q4 );
        Dflux.set(i,5,5, 1, q2*gamma/q1                        );

        //////////////////////
        // Jacobian, g'(q): //
        //////////////////////
        // pd{ g_1 }{ q_j }
        Dflux.set(i,1,3, 2, 1.0e0);

        // pd{ g_2 }{ q_j }
        Dflux.set(i,2,1, 2, Dflux.get(i,3,1,1 ) );
        Dflux.set(i,2,2, 2, Dflux.get(i,3,2,1 ) );
        Dflux.set(i,2,3, 2, Dflux.get(i,3,3,1 ) );

        // pd{ g_3 }{ q_j }
        Dflux.set(i,3,1, 2, 
            0.5*(gamma-3.0)*( pow((q3/q1),2) ) +
            0.5*(gamma-1.0)*( pow((q2/q1),2) + 2.0*pow((q4/q1),3) ) );
        Dflux.set(i,3,2, 2, -(q2/q1)*(gamma-1.0)             );
        Dflux.set(i,3,3, 2, -(q3/q1)*(gamma-3.0)             );
        Dflux.set(i,3,4, 2, 1.5*( gamma-1.0 )*pow((q4/q1),2) );
        Dflux.set(i,3,5, 2, gamma-1.0                        );

        // pd{ g_4 }{ q_j }
        Dflux.set(i,4,1, 2, -q3*q4*pow(q1,-2) );
        Dflux.set(i,4,3, 2, q4 / q1 );
        Dflux.set(i,4,4, 2, q3 / q1 );

        // pd{ g_5 }{ q_j }
        Dflux.set(i,5,1, 2, -q5*(gamma*q3*pow(q1,-2))
            +(gamma-1.0)*( pow(q3/q1,3) + (q3/q1)*pow(q2/q1,2) )
            +1.5*(gamma-1.0)*( (q3/q1)*pow(q4/q1,3) )
        );
        Dflux.set(i,5,2, 2, Dflux.get(i,5,3,1)                 );
        Dflux.set(i,5,3, 2, gamma*q5/q1
            + 1.5*(1.0-gamma)*pow(q3/q1,2)
            + 0.5*(1.0-gamma)*( pow(q2/q1,2) + pow(q4/q1,3 ) ) );
        Dflux.set(i,5,4, 2, 1.5*(1.0-gamma)*q3*pow(q4/q1,3)*q4 );
        Dflux.set(i,5,5, 2, q3*gamma/q1                        );

    }

}
