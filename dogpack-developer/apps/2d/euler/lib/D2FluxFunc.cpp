#include <cmath>
#include "tensors.h"
#include "EulerParams.h"



// This is a user-supplied routine that sets the
// Derivative of the Jacobian of the Flux Function.  
//
// In 2d this is a 3-tensor of size numpts x meqn x meqn x meqn
//
//     Euler Equations
//
void D2FluxFunc(const dTensor2& xpts, 
		const dTensor2& Q,
		const dTensor2& Aux,
		dTensor5& D2flux)
{

    const int numpts = xpts.getsize(1);
    const int meqn   = Q.getsize(2);

    D2flux.setall(0.);


    // Gas constant
    double const gamma = 1.4;//eulerParams.gamma;
    for (int i=1; i<=numpts; i++)
    {

        // Not used: u3 ... 
        double q1 = Q.get(i,1);
        double q2 = Q.get(i,2);
        double q3 = Q.get(i,3);
        double q4 = Q.get(i,4);
        double q5 = Q.get(i,5);

        double en = q5;


        //////////////////////
        // Hessian, f''(q): //
        //////////////////////

////////D2flux.set(i,1,1,1,1, 0);
////////D2flux.set(i,1,1,2,1, 0);
////////D2flux.set(i,1,1,3,1, 0);
////////D2flux.set(i,1,1,4,1, 0);
////////D2flux.set(i,1,1,5,1, 0);
////////D2flux.set(i,1,2,1,1, 0);
////////D2flux.set(i,1,2,2,1, 0);
////////D2flux.set(i,1,2,3,1, 0);
////////D2flux.set(i,1,2,4,1, 0);
////////D2flux.set(i,1,2,5,1, 0);
////////D2flux.set(i,1,3,1,1, 0);
////////D2flux.set(i,1,3,2,1, 0);
////////D2flux.set(i,1,3,3,1, 0);
////////D2flux.set(i,1,3,4,1, 0);
////////D2flux.set(i,1,3,5,1, 0);
////////D2flux.set(i,1,4,1,1, 0);
////////D2flux.set(i,1,4,2,1, 0);
////////D2flux.set(i,1,4,3,1, 0);
////////D2flux.set(i,1,4,4,1, 0);
////////D2flux.set(i,1,4,5,1, 0);
////////D2flux.set(i,1,5,1,1, 0);
////////D2flux.set(i,1,5,2,1, 0);
////////D2flux.set(i,1,5,3,1, 0);
////////D2flux.set(i,1,5,4,1, 0);
////////D2flux.set(i,1,5,5,1, 0);
        D2flux.set(i,2,1,1,1, (-gamma*pow(q2,2) - gamma*pow(q3,2) - gamma*pow(q4,2) + 3*pow(q2,2) + pow(q3,2) + pow(q4,2))/pow(q1,3));
        D2flux.set(i,2,1,2,1, q2*(gamma - 3)/pow(q1,2));
        D2flux.set(i,2,1,3,1, q3*(gamma - 1)/pow(q1,2));
        D2flux.set(i,2,1,4,1, q4*(gamma - 1)/pow(q1,2));
//      D2flux.set(i,2,1,5,1, 0);
        D2flux.set(i,2,2,1,1, q2*(gamma - 3)/pow(q1,2));
        D2flux.set(i,2,2,2,1, (-gamma + 3)/q1);
//      D2flux.set(i,2,2,3,1, 0);
//      D2flux.set(i,2,2,4,1, 0);
//      D2flux.set(i,2,2,5,1, 0);
        D2flux.set(i,2,3,1,1, q3*(gamma - 1)/pow(q1,2));
//      D2flux.set(i,2,3,2,1, 0);
        D2flux.set(i,2,3,3,1, (-gamma + 1.0)/q1);
//      D2flux.set(i,2,3,4,1, 0);
//      D2flux.set(i,2,3,5,1, 0);
        D2flux.set(i,2,4,1,1, q4*(gamma - 1)/pow(q1,2));
//      D2flux.set(i,2,4,2,1, 0);
//      D2flux.set(i,2,4,3,1, 0);
        D2flux.set(i,2,4,4,1, (-gamma + 1.0)/q1);
//      D2flux.set(i,2,4,5,1, 0);
//      D2flux.set(i,2,5,1,1, 0);
//      D2flux.set(i,2,5,2,1, 0);
//      D2flux.set(i,2,5,3,1, 0);
//      D2flux.set(i,2,5,4,1, 0);
//      D2flux.set(i,2,5,5,1, 0);
        D2flux.set(i,3,1,1,1, 2*q2*q3/pow(q1,3));
        D2flux.set(i,3,1,2,1, -q3/pow(q1,2));
        D2flux.set(i,3,1,3,1, -q2/pow(q1,2));
//      D2flux.set(i,3,1,4,1, 0);
//      D2flux.set(i,3,1,5,1, 0);
        D2flux.set(i,3,2,1,1, -q3/pow(q1,2));
//      D2flux.set(i,3,2,2,1, 0);
        D2flux.set(i,3,2,3,1, 1/q1);
//      D2flux.set(i,3,2,4,1, 0);
//      D2flux.set(i,3,2,5,1, 0);
        D2flux.set(i,3,3,1,1, -q2/pow(q1,2));
        D2flux.set(i,3,3,2,1, 1/q1);
////////D2flux.set(i,3,3,3,1, 0);
////////D2flux.set(i,3,3,4,1, 0);
////////D2flux.set(i,3,3,5,1, 0);
////////D2flux.set(i,3,4,1,1, 0);
////////D2flux.set(i,3,4,2,1, 0);
////////D2flux.set(i,3,4,3,1, 0);
////////D2flux.set(i,3,4,4,1, 0);
////////D2flux.set(i,3,4,5,1, 0);
////////D2flux.set(i,3,5,1,1, 0);
////////D2flux.set(i,3,5,2,1, 0);
////////D2flux.set(i,3,5,3,1, 0);
////////D2flux.set(i,3,5,4,1, 0);
////////D2flux.set(i,3,5,5,1, 0);
        D2flux.set(i,4,1,1,1, 2*q2*q4/pow(q1,3));
        D2flux.set(i,4,1,2,1, -q4/pow(q1,2));
//      D2flux.set(i,4,1,3,1, 0);
        D2flux.set(i,4,1,4,1, -q2/pow(q1,2));
//      D2flux.set(i,4,1,5,1, 0);
        D2flux.set(i,4,2,1,1, -q4/pow(q1,2));
//      D2flux.set(i,4,2,2,1, 0);
//      D2flux.set(i,4,2,3,1, 0);
        D2flux.set(i,4,2,4,1, 1/q1);
//      D2flux.set(i,4,2,5,1, 0);
//      D2flux.set(i,4,3,1,1, 0);
//      D2flux.set(i,4,3,2,1, 0);
//      D2flux.set(i,4,3,3,1, 0);
//      D2flux.set(i,4,3,4,1, 0);
//      D2flux.set(i,4,3,5,1, 0);
        D2flux.set(i,4,4,1,1, -q2/pow(q1,2));
        D2flux.set(i,4,4,2,1, 1/q1);
//      D2flux.set(i,4,4,3,1, 0);
//      D2flux.set(i,4,4,4,1, 0);
//      D2flux.set(i,4,4,5,1, 0);
//      D2flux.set(i,4,5,1,1, 0);
//      D2flux.set(i,4,5,2,1, 0);
//      D2flux.set(i,4,5,3,1, 0);
//      D2flux.set(i,4,5,4,1, 0);
//      D2flux.set(i,4,5,5,1, 0);
        D2flux.set(i,5,1,1,1, q2*(2*en*gamma*q1 - 2*en*q1 - 3*gamma*pow(q2,2) - 3*gamma*pow(q3,2) - 3*gamma*pow(q4,2) + 2*q1*q5 + 3*pow(q2,2) + 3*pow(q3,2) + 3*pow(q4,2))/pow(q1,4));
        D2flux.set(i,5,1,2,1, (-en*gamma*q1 + en*q1 + 3*gamma*pow(q2,2) + gamma*pow(q3,2) + gamma*pow(q4,2) - q1*q5 - 3*pow(q2,2) - pow(q3,2) - pow(q4,2))/pow(q1,3));
        D2flux.set(i,5,1,3,1, 2*q2*q3*(gamma - 1)/pow(q1,3));
        D2flux.set(i,5,1,4,1, 2*q2*q4*(gamma - 1)/pow(q1,3));
        D2flux.set(i,5,1,5,1, -q2/pow(q1,2));
        D2flux.set(i,5,2,1,1, (-en*gamma*q1 + en*q1 + 3*gamma*pow(q2,2) + gamma*pow(q3,2) + gamma*pow(q4,2) - q1*q5 - 3*pow(q2,2) - pow(q3,2) - pow(q4,2))/pow(q1,3));
        D2flux.set(i,5,2,2,1, 3*q2*(-gamma + 1.0)/pow(q1,2));
        D2flux.set(i,5,2,3,1, q3*(-gamma + 1.0)/pow(q1,2));
        D2flux.set(i,5,2,4,1, q4*(-gamma + 1.0)/pow(q1,2));
        D2flux.set(i,5,2,5,1, 1/q1);
        D2flux.set(i,5,3,1,1, 2*q2*q3*(gamma - 1)/pow(q1,3));
        D2flux.set(i,5,3,2,1, q3*(-gamma + 1.0)/pow(q1,2));
        D2flux.set(i,5,3,3,1, q2*(-gamma + 1.0)/pow(q1,2));
//      D2flux.set(i,5,3,4,1, 0);
//      D2flux.set(i,5,3,5,1, 0);
        D2flux.set(i,5,4,1,1, 2.0*q2*q4*(gamma - 1)/pow(q1,3));
        D2flux.set(i,5,4,2,1, q4*(-gamma + 1.0)/pow(q1,2));
//      D2flux.set(i,5,4,3,1, 0);
        D2flux.set(i,5,4,4,1, q2*(-gamma + 1.0)/pow(q1,2));
//      D2flux.set(i,5,4,5,1, 0);
        D2flux.set(i,5,5,1,1, -q2/pow(q1,2));
        D2flux.set(i,5,5,2,1, 1.0/q1);
//      D2flux.set(i,5,5,3,1, 0);
//      D2flux.set(i,5,5,4,1, 0);
//      D2flux.set(i,5,5,5,1, 0);


        //////////////////////
        // Hessian, g''(q): //
        //////////////////////
///     D2flux.set(i,1,1,1,2, 0);
///     D2flux.set(i,1,1,2,2, 0);
///     D2flux.set(i,1,1,3,2, 0);
///     D2flux.set(i,1,1,4,2, 0);
///     D2flux.set(i,1,1,5,2, 0);
///     D2flux.set(i,1,2,1,2, 0);
///     D2flux.set(i,1,2,2,2, 0);
///     D2flux.set(i,1,2,3,2, 0);
///     D2flux.set(i,1,2,4,2, 0);
///     D2flux.set(i,1,2,5,2, 0);
///     D2flux.set(i,1,3,1,2, 0);
///     D2flux.set(i,1,3,2,2, 0);
///     D2flux.set(i,1,3,3,2, 0);
///     D2flux.set(i,1,3,4,2, 0);
///     D2flux.set(i,1,3,5,2, 0);
///     D2flux.set(i,1,4,1,2, 0);
///     D2flux.set(i,1,4,2,2, 0);
///     D2flux.set(i,1,4,3,2, 0);
///     D2flux.set(i,1,4,4,2, 0);
///     D2flux.set(i,1,4,5,2, 0);
///     D2flux.set(i,1,5,1,2, 0);
///     D2flux.set(i,1,5,2,2, 0);
///     D2flux.set(i,1,5,3,2, 0);
///     D2flux.set(i,1,5,4,2, 0);
///     D2flux.set(i,1,5,5,2, 0);
         
        D2flux.set(i,2,1,1,2, 2*q2*q3/pow(q1,3));
        D2flux.set(i,2,1,2,2, -q3/pow(q1,2));
        D2flux.set(i,2,1,3,2, -q2/pow(q1,2));
//      D2flux.set(i,2,1,4,2, 0);
//      D2flux.set(i,2,1,5,2, 0);
        D2flux.set(i,2,2,1,2, -q3/pow(q1,2));
//      D2flux.set(i,2,2,2,2, 0);
        D2flux.set(i,2,2,3,2, 1/q1);
//      D2flux.set(i,2,2,4,2, 0);
//      D2flux.set(i,2,2,5,2, 0);
        D2flux.set(i,2,3,1,2, -q2/pow(q1,2));
        D2flux.set(i,2,3,2,2, 1/q1);
//      D2flux.set(i,2,3,3,2, 0);
//      D2flux.set(i,2,3,4,2, 0);
//      D2flux.set(i,2,3,5,2, 0);
//      D2flux.set(i,2,4,1,2, 0);
//      D2flux.set(i,2,4,2,2, 0);
//      D2flux.set(i,2,4,3,2, 0);
//      D2flux.set(i,2,4,4,2, 0);
//      D2flux.set(i,2,4,5,2, 0);
//      D2flux.set(i,2,5,1,2, 0);
//      D2flux.set(i,2,5,2,2, 0);
//      D2flux.set(i,2,5,3,2, 0);
//      D2flux.set(i,2,5,4,2, 0);
//      D2flux.set(i,2,5,5,2, 0);
         
        D2flux.set(i,3,1,1,2, (-gamma*pow(q2,2) - gamma*pow(q3,2) - gamma*pow(q4,2) + pow(q2,2) + 3*pow(q3,2) + pow(q4,2))/pow(q1,3));
        D2flux.set(i,3,1,2,2, q2*(gamma - 1)/pow(q1,2));
        D2flux.set(i,3,1,3,2, q3*(gamma - 3)/pow(q1,2));
        D2flux.set(i,3,1,4,2, q4*(gamma - 1)/pow(q1,2));
//      D2flux.set(i,3,1,5,2, 0);
        D2flux.set(i,3,2,1,2, q2*(gamma - 1)/pow(q1,2));
        D2flux.set(i,3,2,2,2, (-gamma + 1)/q1);
//      D2flux.set(i,3,2,3,2, 0);
//      D2flux.set(i,3,2,4,2, 0);
//      D2flux.set(i,3,2,5,2, 0);
        D2flux.set(i,3,3,1,2, q3*(gamma - 3)/pow(q1,2));
//      D2flux.set(i,3,3,2,2, 0);
        D2flux.set(i,3,3,3,2, (-gamma + 3)/q1);
//      D2flux.set(i,3,3,4,2, 0);
//      D2flux.set(i,3,3,5,2, 0);
        D2flux.set(i,3,4,1,2, q4*(gamma - 1)/pow(q1,2));
//      D2flux.set(i,3,4,2,2, 0);
//      D2flux.set(i,3,4,3,2, 0);
        D2flux.set(i,3,4,4,2, (-gamma + 1)/q1);
//      D2flux.set(i,3,4,5,2, 0);
//      D2flux.set(i,3,5,1,2, 0);
//      D2flux.set(i,3,5,2,2, 0);
//      D2flux.set(i,3,5,3,2, 0);
//      D2flux.set(i,3,5,4,2, 0);
//      D2flux.set(i,3,5,5,2, 0);
         
        D2flux.set(i,4,1,1,2, 2*q3*q4/pow(q1,3));
//      D2flux.set(i,4,1,2,2, 0);
        D2flux.set(i,4,1,3,2, -q4/pow(q1,2));
        D2flux.set(i,4,1,4,2, -q3/pow(q1,2));
//      D2flux.set(i,4,1,5,2, 0);
//      D2flux.set(i,4,2,1,2, 0);
//      D2flux.set(i,4,2,2,2, 0);
//      D2flux.set(i,4,2,3,2, 0);
//      D2flux.set(i,4,2,4,2, 0);
//      D2flux.set(i,4,2,5,2, 0);
        D2flux.set(i,4,3,1,2, -q4/pow(q1,2));
//      D2flux.set(i,4,3,2,2, 0);
//      D2flux.set(i,4,3,3,2, 0);
        D2flux.set(i,4,3,4,2, 1/q1);
//      D2flux.set(i,4,3,5,2, 0);
        D2flux.set(i,4,4,1,2, -q3/pow(q1,2));
//      D2flux.set(i,4,4,2,2, 0);
        D2flux.set(i,4,4,3,2, 1/q1);
//      D2flux.set(i,4,4,4,2, 0);
//      D2flux.set(i,4,4,5,2, 0);
//      D2flux.set(i,4,5,1,2, 0);
//      D2flux.set(i,4,5,2,2, 0);
//      D2flux.set(i,4,5,3,2, 0);
//      D2flux.set(i,4,5,4,2, 0);
//      D2flux.set(i,4,5,5,2, 0);
         
        D2flux.set(i,5,1,1,2, q3*(2*en*gamma*q1 - 2*en*q1 - 3*gamma*pow(q2,2) - 3*gamma*pow(q3,2) - 3*gamma*pow(q4,2) + 2*q1*q5 + 3*pow(q2,2) + 3*pow(q3,2) + 3*pow(q4,2))/pow(q1,4));
        D2flux.set(i,5,1,2,2, 2*q2*q3*(gamma - 1)/pow(q1,3));
        D2flux.set(i,5,1,3,2, (-en*gamma*q1 + en*q1 + gamma*pow(q2,2) + 3*gamma*pow(q3,2) + gamma*pow(q4,2) - q1*q5 - pow(q2,2) - 3*pow(q3,2) - pow(q4,2))/pow(q1,3));
        D2flux.set(i,5,1,4,2, 2*q3*q4*(gamma - 1)/pow(q1,3));
        D2flux.set(i,5,1,5,2, -q3/pow(q1,2));
        D2flux.set(i,5,2,1,2, 2*q2*q3*(gamma - 1)/pow(q1,3));
        D2flux.set(i,5,2,2,2, q3*(-gamma + 1)/pow(q1,2));
        D2flux.set(i,5,2,3,2, q2*(-gamma + 1)/pow(q1,2));
//      D2flux.set(i,5,2,4,2, 0);
//      D2flux.set(i,5,2,5,2, 0);
        D2flux.set(i,5,3,1,2, (-en*gamma*q1 + en*q1 + gamma*pow(q2,2) + 3*gamma*pow(q3,2) + gamma*pow(q4,2) - q1*q5 - pow(q2,2) - 3*pow(q3,2) - pow(q4,2))/pow(q1,3));
        D2flux.set(i,5,3,2,2, q2*(-gamma + 1)/pow(q1,2));
        D2flux.set(i,5,3,3,2, 3*q3*(-gamma + 1)/pow(q1,2));
        D2flux.set(i,5,3,4,2, q4*(-gamma + 1)/pow(q1,2));
        D2flux.set(i,5,3,5,2, 1/q1);
        D2flux.set(i,5,4,1,2, 2*q3*q4*(gamma - 1)/pow(q1,3));
//      D2flux.set(i,5,4,2,2, 0);
        D2flux.set(i,5,4,3,2, q4*(-gamma + 1)/pow(q1,2));
        D2flux.set(i,5,4,4,2, q3*(-gamma + 1)/pow(q1,2));
//      D2flux.set(i,5,4,5,2, 0);
        D2flux.set(i,5,5,1,2, -q3/pow(q1,2));
//      D2flux.set(i,5,5,2,2, 0);
        D2flux.set(i,5,5,3,2, 1/q1);
//      D2flux.set(i,5,5,4,2, 0);
//      D2flux.set(i,5,5,5,2, 0);

    }

}
