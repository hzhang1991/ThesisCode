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

    for (int n=1; n<=numpts; n++)
    {
        double x = xpts.get( n, 1);
        double y = xpts.get( n, 2);

        for( int m1=1; m1 <= meqn; m1++ )
        for( int m2=1; m2 <= meqn; m2++ )
        {
            // Variables
            double q = Q.get(n, m2);

            // f'(q):
            Dflux.set( n, m1, m2, 1, 0. );

            // g'(q):
            Dflux.set( n, m1, m2, 2, 0. );
        }

    }

    // We will now print an error, because users need to implement this
    // function for each problem they want to solve using the Lax-Wendroff
    // method.
    printf("Error in DFluxFunc.cpp: " );
    printf("  In order to use Lax-Wendroff, you need to implement this function.\n");
    exit(1);

}
