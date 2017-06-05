#include <cmath>
#include <iostream>
#include <fstream>
#include "dogdefs.h"
#include "EulerParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//


//In this file we specify initial conditions for a number of
//2d Riemann problems that have been studied in many papers,
//see for example "Solution of Two-Dimensional Riemann Problems
//for Gas Dynamics without Riemann Problem
//Solvers" by Kurganov and Tadmor. We use
// the indicator riemann to indicate which Riemann problem
// we wish to solve.
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    const int numpts=xpts.getsize(1);

    // Gas consant
    const double gamma = 1.4;eulerParams.gamma;
    const double x0 = eulerParams.x0;

    // Loop over grid points
    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double y = xpts.get(i,2);

        double xstar = x0 + y*osq3;
       
        double rho,u1,u2,u3,press;

        double rho1,rho2,rho3,rho4;
        double u11,u12,u13,u14;
        double u21,u22,u23,u24;
        double press1,press2,press3,press4;
  
        double xcorner=0.0;double ycorner=0.0;

        int center=0;       

        if(center==0)
        {
            rho1=0.532258064516129;
            u11 = 1.206045378311055;
            u21 = 0.0;
            u3 = 0.0;
            press1=0.3;

            rho2   = 0.137992831541219;
            u12    = 1.206045378311055;
            u22    = 1.206045378311055;
            u3    =  0.0;
            press2 =  0.029032258064516;

            rho3   =  0.532258064516129;
            u13    =  0.0;
            u23    =  0.532258064516129;
            u3    =  0.0;
            press3 =  0.3;

            rho4   = 1.5;
            u14    = 0.0;
            u24    = 0.0;
            u3    = 0.0;
            press4 = 1.5;

            xcorner=0.3;
            ycorner=0.3;
        }

        const int riemann=1;
        if(riemann==1)
        {

            rho1=0.5323;
            u11 = 1.206;
            u21 = 0.0;
            u3 = 0.0;
            press1=0.3;

            rho2   = 0.138;
            u12    = 1.206;
            u22    = 1.206;
            u3    =  0.0;
            press2 =  0.029;

            rho3   =  0.5323;
            u13    =  0.0;
            u23    =  1.206;
            u3    =  0.0;
            press3 =  0.3;

            rho4   = 1.5;
            u14    = 0.0;
            u24    = 0.0;
            u3    = 0.0;
            press4 = 1.5;
        }
        if(riemann==2)
        {

            rho1=0.5065;
            u11 = 0.8939;
            u21 = 0.0;
            u3 = 0.0;
            press1=0.35;

            rho2   = 1.1;
            u12    =  0.8939;
            u22    =  0.8939;
            u3    =  0.0;
            press2 =  1.1;

            rho3   =  0.5065;
            u13    =  0.0;
            u23    =  0.8939;
            u3    =  0.0;
            press3 =  0.35;

            rho4   = 1.1;
            u14    = 0.0;
            u24    = 0.0;
            u3    = 0.0;
            press4 = 1.1;
        }
        if(riemann==3)
        {

            rho1=2.0;
            u11 = 0.75;
            u21 = 0.5;
            u3 = 0.0;
            press1=1.0;

            rho2   = 1.0;
            u12    =  -0.75;
            u22    =  0.5;
            u3    =  0.0;
            press2 =  1.0;

            rho3   =  3.0;
            u13    =  -0.75;
            u23    =  -0.5;
            u3    =  0.0;
            press3 =  1.0;

            rho4   = 1.0;
            u14    = 0.75;
            u24    = -0.5;
            u3    = 0.0;
            press4 = 1.0;
        }
        if(riemann==4)
        {

            rho1=1.0;
            u11 = -0.6259;
            u21 = 0.1;
            u3 = 0.0;
            press1=1.0;

            rho2   = 0.8;
            u12    =  0.1;
            u22    =  0.1;
            u3    =  0.0;
            press2 =  1.0;

            rho3   =  1.0;
            u13    =  0.1;
            u23    =  -0.6259;
            u3    =  0.0;
            press3 =  1.0;

            rho4   = 0.5917;
            u14    = 0.1;
            u24    = 0.1;
            u3    = 0.0;
            press4 = 0.4;
        }
        if(riemann==5)
        {

            rho1=1.0;
            u11 = 0.7276;
            u21 = 0.0;
            u3 = 0.0;
            press1=1.0;

            rho2   = 0.8;
            u12    =  0.0;
            u22    =  0.0;
            u3    =  0.0;
            press2 =  1.0;

            rho3   =  1.0;
            u13    =  0.0;
            u23    =  0.7276;
            u3    =  0.0;
            press3 =  1.0;

            rho4   = 0.5313;
            u14    = 0.0;
            u24    = 0.0;
            u3    = 0.0;
            press4 = 0.4;
        }




        if(x < xcorner)
        {
            if(y<ycorner)
            {rho   = rho2;
            u1    =  u12;
            u2    =  u22;
            u3    =  0.0;
            press =  press2;}
            else
            {
              rho=rho1;
              u1 = u11;
              u2 = u21;
              u3 = 0.0;
              press=press1;
            }
        }
        else
        {
            if(y<ycorner)
            {
              rho   =  rho3;
              u1    =  u13;
              u2    =  u23;
              u3    =  0.0;
              press =  press3;
            }
            else
            {
              rho   = rho4;
              u1    = u14;
              u2    = u24;
              u3    = 0.0;
              press = press4;
            }
             
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
