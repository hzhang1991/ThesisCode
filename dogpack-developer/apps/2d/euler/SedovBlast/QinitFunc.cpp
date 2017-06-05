#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "EulerParams.h"
#include "DogParamsCart2.h"
#include <iostream>
using namespace std;

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    //cout<<eulerParams.gamma<<" "<<dy<<endl;
    // Gas consant
    const double gamma = 1.4;

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double rad = sqrt(pow(x,2) + pow(y,2));

        double rho,u1,u2,u3,press,energy;
        if(fabs(x)<=dx && fabs(y)<=dy)
        {
            rho   =  1.0;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            energy =  0.244816/(dx*dy);
        }
        else
        {
            rho   =  1.0;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            energy =  1.0e-12;
        }

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
