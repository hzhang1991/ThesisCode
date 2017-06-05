#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// right-hand side function at all the points "xpts"
//
void RhsFunc(const dTensor2& xpts, 
        const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2, 
        dTensor2& rhs)
{

    const int numpts=xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        rhs.set(i,1, 1.0 );
    }

}
