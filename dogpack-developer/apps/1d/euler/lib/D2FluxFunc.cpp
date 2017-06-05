#include <cmath>
#include "tensors.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
//     ** Euler Equations **
//
// In 1d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k) = \partial^2 f_i / \partial q_j \partial q_k
//
// Inputs:
//
//     xpts( 1:numpts )         - a list of x-values at various spatial points.
//        Q( 1:numpts, 1:meqn ) - a vector of conserved variables
//      Aux( 1:numpts, 1:maux ) - vector of auxilary values
//   
// Output:
//
//    D2flux( 1:numpts, 1:meqn, 1:meqn, 1:meqn ) - f''(q) at each point.
//
// See also: FluxFunc and DFluxFunc.
void D2FluxFunc(const dTensor1& xpts, 
        const dTensor2& Q,
        const dTensor2& Aux,
        dTensor4& D2flux)
{

    // Gas constant
    const double gamma = eulerParams.gamma;
    const double gm1   = gamma - 1.0;

    D2flux.setall(0.);  // most terms are zero

    const int numpts=xpts.getsize();
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);      

        // (Characteristic) Variables
//      const double rho    = Q.get(i,1);
//      const double u1     = Q.get(i,2)/rho;
//      const double u2     = Q.get(i,3)/rho;
//      const double u3     = Q.get(i,4)/rho;
//      const double energy = Q.get(i,5);
//      const double press  = (gamma-1.0e0)*(energy-0.5e0*rho*(u1*u1+u2*u2+u3*u3));

        // (Conserved) Variables
        const double q1 = Q.get(i,1);
        const double q2 = Q.get(i,2);
        const double q3 = Q.get(i,3);
        const double q4 = Q.get(i,4);
        const double q5 = Q.get(i,5);

        // Computing the Hessian of the flux function: f''(q)
        D2flux.set(i,1,1,1, 0);
        D2flux.set(i,1,1,2, 0);
        D2flux.set(i,1,1,3, 0);
        D2flux.set(i,1,1,4, 0);
        D2flux.set(i,1,1,5, 0);
        D2flux.set(i,1,2,1, 0);
        D2flux.set(i,1,2,2, 0);
        D2flux.set(i,1,2,3, 0);
        D2flux.set(i,1,2,4, 0);
        D2flux.set(i,1,2,5, 0);
        D2flux.set(i,1,3,1, 0);
        D2flux.set(i,1,3,2, 0);
        D2flux.set(i,1,3,3, 0);
        D2flux.set(i,1,3,4, 0);
        D2flux.set(i,1,3,5, 0);
        D2flux.set(i,1,4,1, 0);
        D2flux.set(i,1,4,2, 0);
        D2flux.set(i,1,4,3, 0);
        D2flux.set(i,1,4,4, 0);
        D2flux.set(i,1,4,5, 0);
        D2flux.set(i,1,5,1, 0);
        D2flux.set(i,1,5,2, 0);
        D2flux.set(i,1,5,3, 0);
        D2flux.set(i,1,5,4, 0);
        D2flux.set(i,1,5,5, 0);
         
        D2flux.set(i,2,1,1, -gamma*pow( q2, 2 )/pow( q1, 3 ) - gamma*pow( q3, 2 )/pow( q1, 3 ) - gamma*pow( q4, 2 )/pow( q1, 3 ) + 3.0*pow( q2, 2 )/pow( q1, 3 ) + pow( q3, 2 )/pow( q1, 3 ) + pow( q4, 2 )/pow( q1, 3 ));
        D2flux.set(i,2,1,2, gamma*q2/pow( q1, 2 ) - 3.0*q2/pow( q1, 2 ));
        D2flux.set(i,2,1,3, gamma*q3/pow( q1, 2 ) - q3/pow( q1, 2 ));
        D2flux.set(i,2,1,4, gamma*q4/pow( q1, 2 ) - q4/pow( q1, 2 ));
        D2flux.set(i,2,1,5, 0);
        D2flux.set(i,2,2,1, gamma*q2/pow( q1, 2 ) - 3.0*q2/pow( q1, 2 ));
        D2flux.set(i,2,2,2, -gamma/q1 + 3.0/q1);
        D2flux.set(i,2,2,3, 0);
        D2flux.set(i,2,2,4, 0);
        D2flux.set(i,2,2,5, 0);
        D2flux.set(i,2,3,1, gamma*q3/pow( q1, 2 ) - q3/pow( q1, 2 ));
        D2flux.set(i,2,3,2, 0);
        D2flux.set(i,2,3,3, -gamma/q1 + 1/q1);
        D2flux.set(i,2,3,4, 0);
        D2flux.set(i,2,3,5, 0);
        D2flux.set(i,2,4,1, gamma*q4/pow( q1, 2 ) - q4/pow( q1, 2 ));
        D2flux.set(i,2,4,2, 0);
        D2flux.set(i,2,4,3, 0);
        D2flux.set(i,2,4,4, -gamma/q1 + 1/q1);
        D2flux.set(i,2,4,5, 0);
        D2flux.set(i,2,5,1, 0);
        D2flux.set(i,2,5,2, 0);
        D2flux.set(i,2,5,3, 0);
        D2flux.set(i,2,5,4, 0);
        D2flux.set(i,2,5,5, 0);
         
        D2flux.set(i,3,1,1, 2.0*q2*q3/pow( q1, 3 ));
        D2flux.set(i,3,1,2, -q3/pow( q1, 2 ));
        D2flux.set(i,3,1,3, -q2/pow( q1, 2 ));
        D2flux.set(i,3,1,4, 0);
        D2flux.set(i,3,1,5, 0);
        D2flux.set(i,3,2,1, -q3/pow( q1, 2 ));
        D2flux.set(i,3,2,2, 0);
        D2flux.set(i,3,2,3, 1/q1);
        D2flux.set(i,3,2,4, 0);
        D2flux.set(i,3,2,5, 0);
        D2flux.set(i,3,3,1, -q2/pow( q1, 2 ));
        D2flux.set(i,3,3,2, 1/q1);
        D2flux.set(i,3,3,3, 0);
        D2flux.set(i,3,3,4, 0);
        D2flux.set(i,3,3,5, 0);
        D2flux.set(i,3,4,1, 0);
        D2flux.set(i,3,4,2, 0);
        D2flux.set(i,3,4,3, 0);
        D2flux.set(i,3,4,4, 0);
        D2flux.set(i,3,4,5, 0);
        D2flux.set(i,3,5,1, 0);
        D2flux.set(i,3,5,2, 0);
        D2flux.set(i,3,5,3, 0);
        D2flux.set(i,3,5,4, 0);
        D2flux.set(i,3,5,5, 0);
         
        D2flux.set(i,4,1,1, 2.0*q2*q4/pow( q1, 3 ));
        D2flux.set(i,4,1,2, -q4/pow( q1, 2 ));
        D2flux.set(i,4,1,3, 0);
        D2flux.set(i,4,1,4, -q2/pow( q1, 2 ));
        D2flux.set(i,4,1,5, 0);
        D2flux.set(i,4,2,1, -q4/pow( q1, 2 ));
        D2flux.set(i,4,2,2, 0);
        D2flux.set(i,4,2,3, 0);
        D2flux.set(i,4,2,4, 1/q1);
        D2flux.set(i,4,2,5, 0);
        D2flux.set(i,4,3,1, 0);
        D2flux.set(i,4,3,2, 0);
        D2flux.set(i,4,3,3, 0);
        D2flux.set(i,4,3,4, 0);
        D2flux.set(i,4,3,5, 0);
        D2flux.set(i,4,4,1, -q2/pow( q1, 2 ));
        D2flux.set(i,4,4,2, 1/q1);
        D2flux.set(i,4,4,3, 0);
        D2flux.set(i,4,4,4, 0);
        D2flux.set(i,4,4,5, 0);
        D2flux.set(i,4,5,1, 0);
        D2flux.set(i,4,5,2, 0);
        D2flux.set(i,4,5,3, 0);
        D2flux.set(i,4,5,4, 0);
        D2flux.set(i,4,5,5, 0);
         
        D2flux.set(i,5,1,1, 2.0*q5*gamma*q2/pow( q1, 3 ) - 2.0*q5*q2/pow( q1, 3 ) - 3.0*gamma*pow( q2, 3 )/pow( q1, 4 ) - 3.0*gamma*q2*pow( q3, 2 )/pow( q1, 4 ) - 3.0*gamma*q2*pow( q4, 2 )/pow( q1, 4 ) + 2.0*q2*q5/pow( q1, 3 ) + 3.0*pow( q2, 3 )/pow( q1, 4 ) + 3.0*q2*pow( q3, 2 )/pow( q1, 4 ) + 3.0*q2*pow( q4, 2 )/pow( q1, 4 ));
        D2flux.set(i,5,1,2, -q5*gamma/pow( q1, 2 ) + q5/pow( q1, 2 ) + 3.0*gamma*pow( q2, 2 )/pow( q1, 3 ) + gamma*pow( q3, 2 )/pow( q1, 3 ) + gamma*pow( q4, 2 )/pow( q1, 3 ) - q5/pow( q1, 2 ) - 3.0*pow( q2, 2 )/pow( q1, 3 ) - pow( q3, 2 )/pow( q1, 3 ) - pow( q4, 2 )/pow( q1, 3 ));
        D2flux.set(i,5,1,3, 2.0*gamma*q2*q3/pow( q1, 3 ) - 2.0*q2*q3/pow( q1, 3 ));
        D2flux.set(i,5,1,4, 2.0*gamma*q2*q4/pow( q1, 3 ) - 2.0*q2*q4/pow( q1, 3 ));
        D2flux.set(i,5,1,5, -q2/pow( q1, 2 ));
        D2flux.set(i,5,2,1, -q5*gamma/pow( q1, 2 ) + q5/pow( q1, 2 ) + 3.0*gamma*pow( q2, 2 )/pow( q1, 3 ) + gamma*pow( q3, 2 )/pow( q1, 3 ) + gamma*pow( q4, 2 )/pow( q1, 3 ) - q5/pow( q1, 2 ) - 3.0*pow( q2, 2 )/pow( q1, 3 ) - pow( q3, 2 )/pow( q1, 3 ) - pow( q4, 2 )/pow( q1, 3 ));
        D2flux.set(i,5,2,2, -3.0*gamma*q2/pow( q1, 2 ) + 3.0*q2/pow( q1, 2 ));
        D2flux.set(i,5,2,3, -gamma*q3/pow( q1, 2 ) + q3/pow( q1, 2 ));
        D2flux.set(i,5,2,4, -gamma*q4/pow( q1, 2 ) + q4/pow( q1, 2 ));
        D2flux.set(i,5,2,5, 1/q1);
        D2flux.set(i,5,3,1, 2.0*gamma*q2*q3/pow( q1, 3 ) - 2.0*q2*q3/pow( q1, 3 ));
        D2flux.set(i,5,3,2, -gamma*q3/pow( q1, 2 ) + q3/pow( q1, 2 ));
        D2flux.set(i,5,3,3, -gamma*q2/pow( q1, 2 ) + q2/pow( q1, 2 ));
        D2flux.set(i,5,3,4, 0);
        D2flux.set(i,5,3,5, 0);
        D2flux.set(i,5,4,1, 2.0*gamma*q2*q4/pow( q1, 3 ) - 2.0*q2*q4/pow( q1, 3 ));
        D2flux.set(i,5,4,2, -gamma*q4/pow( q1, 2 ) + q4/pow( q1, 2 ));
        D2flux.set(i,5,4,3, 0);
        D2flux.set(i,5,4,4, -gamma*q2/pow( q1, 2 ) + q2/pow( q1, 2 ));
        D2flux.set(i,5,4,5, 0);
        D2flux.set(i,5,5,1, -q2/pow( q1, 2 ));
        D2flux.set(i,5,5,2, 1/q1);
        D2flux.set(i,5,5,3, 0);
        D2flux.set(i,5,5,4, 0);
        D2flux.set(i,5,5,5, 0);
 
    }

}
