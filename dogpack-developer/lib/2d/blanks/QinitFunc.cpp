#include "tensors.h"

// Template for the initial conditions that are defined by each application.
//
// Parameters:
// ----------
//
//     xpts =    xpts( 1:numpts, 1:NDIMS )
//
// Each problem knows what dimension it's in.  Users need not concern
// themselves with what points are used here for the initial projection onto
// the basis functions.
//
// Returns:
// --------
//
//  qvals = qvals( 1:numpts, 1:MAUX  )
//
// See also: AuxFunc.cpp

void QinitFunc(const dTensor2& xpts, 
	       dTensor2& qvals)
{

    const int numpts = xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        // Set any values into qvals
        //
        // Each application knows meqn, so no real need for this for loop:
        //
        const int meqn = qvals.getsize(2);
        for( int me=1; me <= meqn; me++ )
        {
            qvals.set(i, me, 1.0 );
        }
    }

}
