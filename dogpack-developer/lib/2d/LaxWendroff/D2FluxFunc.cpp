#include "dogdefs.h"

// This is a user-supplied routine that sets the Hessian
// of the flux function.
//
// In 2d these are two 3-tensors of size numpts x meqn x meqn x meqn
//
//     D2flux(:, i, j, k, 1) = \partial^2 f_i / \partial q_j \partial q_k
//     D2flux(:, i, j, k, 2) = \partial^2 g_i / \partial q_j \partial q_k
//
// Template for defining this in an application's folder.
//
void D2FluxFunc(const dTensor2& xpts, 
                const dTensor2& Q,
                const dTensor2& Aux, 
                dTensor5& D2flux )
{

    const int numpts = xpts.getsize(1);
    const int meqn   = Q.getsize(2);

    for (int n=1; n<=numpts; n++)
    {
        double x = xpts.get( n, 1);
        double y = xpts.get( n, 2);

        for( int m1=1; m1 <= meqn; m1++ )
        for( int m2=1; m2 <= meqn; m2++ )
        for( int m3=1; m3 <= meqn; m3++ )
        {
            // Variables
            double qj = Q.get(n, m2);
            double qk = Q.get(n, m3);

            // f'(q):
            D2flux.set( n, m1, m2, m3, 1, 0. );

            // g'(q):
            D2flux.set( n, m1, m2, m3, 2, 0. );
        }

    }

    // We will now print an error, because users need to implement this
    // function for each problem they want to solve using the Lax-Wendroff
    // method.
    printf("Error in D2FluxFunc.cpp: " );
    printf("  In order to use third-order Lax-Wendroff, you need to implement this function.\n");
    exit(1);

}
