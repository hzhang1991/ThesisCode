#include "dogdefs.h"

// This is a user-supplied routine that sets the
// Jacobian for each of the two Flux functions.
//
// For the 2D hyperbolic conservation law, q_t + f_x + g_y = 0,
// this function defines f'(q) and g'(q).
//
// The expected format is Dflux.get(:, i, j, 1) = \partial f_i, \partial q_j,
// and                    Dflux.get(:, i, j, 2) = \partial g_i, \partial q_j.
//
// The first index is an aribrary list of points needed to perform L2-Project.
//
//     Simple advection equation, f'(q) = u, g'(q) = v
//     Burger's equation,         f'(q) = q, g'(q) = q
//     Acoustics equation,        f'(q) = [0 1; 1 0], g'(q) = [0 1; 1 0].
//
//
void DFluxFunc(const dTensor2& xpts, 
               const dTensor2& Q,
               const dTensor2& Aux, 
               dTensor4& Dflux )
{

    const int numpts = xpts.getsize(1);
    const int meqn   = Q.getsize(2);

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get( i, 1);
        double y = xpts.get( i, 2);

        // Variables
        double q = Q.get(i, 1);

        // f'(q):
        Dflux.set( i, 1, 1, 1, Aux.get(i,1) );

        // g'(q):
        Dflux.set( i, 1, 1, 2, Aux.get(i,2) );

    }

}
