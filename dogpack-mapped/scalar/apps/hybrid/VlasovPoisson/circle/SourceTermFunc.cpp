#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc(
    const double t, 
    const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& qvals, 
    const dTensor2& auxvals, dTensor2& source)
{

    const int numpts=xpts.getsize(1);
    const int meqn=source.getsize(2);

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        for (int m=1; m<=meqn; m++)
        {
            source.set(i,m, 0.0e0 );
        }
    }

}
