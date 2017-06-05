#include "tensors.h"

// Template for flux function that is defined by each application.
//
// All lengths, ndims, meqn, maux are defined in the parameters.ini file.
//
// Parameters:
// ----------
//
//     xpts =    xpts( 1:numpts, 1:NDIMS )
//        Q =       Q( 1:numpts, 1:MEQN  )
//      Aux =     Aux( 1:numpts, 1:MAUX  )
//
// Each problem knows what dimension it's in.  Users need not concern
// themselves with what points are used here.
//
// Returns:
// --------
//
//  flux = flux( 1:numpts, 1:MEQN, 1:NDIMS )
//
void FluxFunc(const dTensor2& xpts, const dTensor2& Q, 
        const dTensor2& Aux, dTensor3& flux)
{

    const int numpts = xpts.getsize(1);
    const int meqn   =    Q.getsize(2);
    for (int i=1; i<=numpts; i++)
    {

        for( int me=1; me <= meqn; me++ )
        {
            // 1-component of flux function
            flux.set(i, me, 1, 0. );

            // 2-component of flux function
            flux.set(i, me, 2, 0. );
        }
    }

}
