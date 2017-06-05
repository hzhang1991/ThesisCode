#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize(1);

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {

        // Variables
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double r=sqrt(x*x+y*y);

        // exact solution for p: cos(2*pi*t)*sin(2*pi*x)*sin(4*pi*y):
        //double p = exp(-200.0*pow(r-0.15,2.0))*sin(150.0*r);
        double p = exp(-200.0*pow(r-0.15,2.0));//*sin(150.0*r);
        double u = 0.0;
        double v = 0.0;

        qvals.set(i,1, p );
        qvals.set(i,2, u );
        qvals.set(i,3, v );

    }

}
