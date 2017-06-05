#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);

    // Gas consant
    const double gamma = 1.4;//eulerParams.gamma;
    const double x0 = 1.0/6.0;//eulerParams.x0;

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double xstar = x0 + y*osq3;

        double rho,u1,u2,u3,press;

        if(x < 0.5)//xstar)
        {
            rho   =  8.0;
            u1    =  8.25;
            u2    =  0.0;//-8.25*sin(pi/6.0);
            u3    =  0.0;
            press =  116.5;
        }
        else
        {
            rho   =  1.4;
            u1    =  0.0;
            u2    =  0.0;
            u3    =  0.0;
            press =  1.0;
        } 

        double energy = press/(gamma-1.0e0) 
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
