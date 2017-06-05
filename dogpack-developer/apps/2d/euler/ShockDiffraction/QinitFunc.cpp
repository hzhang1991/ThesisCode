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


void qinitfuncs1(double x,double y,double &rho,double &press,double &u1,double &u2,double &u3,double &energy)
{

    const double gamma = 1.4;
    const double gm1   = 1.4 - 1.0;
    const double gp1   = 1.4 + 1.0;

    const double Mach  = 5.09;
    const double M2    = Mach*Mach;

    const double vs = Mach*sqrt( gamma*press/rho );

    // Correct with the Mach number for incoming shock
    {
        press = press*(2.0*gamma*M2- gm1)/gp1;
        rho   = rho*gp1*M2 / (gm1*M2 + 2.0);
        u1    = vs*(1.0 - 1.4 / rho );
        energy = press/(gamma-1.0e0)
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);
    }
}

void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    void qinitfuncs1(double x,double y,double &rho,double &press,double &u1,double &u2,double &u3,double &energy);

    const int numpts=xpts.getsize(1);
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double gamma = 1.4;

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double rad = sqrt(pow(x,2) + pow(y,2));

        double  rho   =  1.4;
        double  u1    =  0.0;
        double  u2    =  0.0;
        double  u3    =  0.0;
        double  press =  1.0;
        double energy = press/(gamma-1.0e0)
            + 0.5e0*rho*(u1*u1 + u2*u2 + u3*u3);

        if(x<0.5 && y>6.0)
        {
            qinitfuncs1(x,y,rho,press,u1,u2,u3,energy);
        }

        qvals.set(i,1, rho    );
        qvals.set(i,2, rho*u1 );
        qvals.set(i,3, rho*u2 );
        qvals.set(i,4, rho*u3 );
        qvals.set(i,5, energy );
    }
}
