#include "tensors.h"

// Template for Auxilary function that is defined by each application.
//
// Parameters:
// ----------
//
//     xpts =    xpts( 1:numpts, 1:NDIMS )
//
// Each problem knows what dimension it's in.  Users need not concern
// themselves with what points are used here.
//
// For example, in dimension 2, xpts( :, 1 ) = x and xpts( :, 2 ) = y.
//
// Returns:
// --------
//
//  auxvals = auxvals( 1:numpts, 1:MAUX  )
//
// See also: QinitFunc.cpp

void AuxFunc(const dTensor2& xpts, 
	     dTensor2& auxvals)
{

    const int numpts = xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        // Set any values into auxvals (user defines maux in the
        // parameters.ini file)
        const int maux = auxvals.getsize(2);
        for( int ma=1; ma <= auxvals.getsize(2); ma++ )
        {
            auxvals.set(i,ma, 1.0 );
        }
    }

}
