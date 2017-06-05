#include "tensors.h"
#include "dog_math.h"
#include <cmath>
#include "dogdefs.h"
#include <iostream>
using namespace std;

// This is the HLL(E) method, written specifically for LaxWendroff time
// stepping.
//
// In the Lax-Wendroff method, we already have access to flux function values on
// the "left" and "right" side of the interface, so there is no need to call the
// flux function to set those.
//
// See also: $DOGPACK/lib/2d/RiemannSolve.
double RiemannSolveLxW(const dTensor1& nvec,
    const dTensor1& xedge,
    const dTensor1& Ql,   const dTensor1& Qr,
    const dTensor1& Auxl, const dTensor1& Auxr,
    const dTensor1& ffl,  const dTensor1& ffr,
    dTensor1& Fl, dTensor1& Fr)
{

    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();

    // HLLE parameters
    double s1    = 0.;
    double s2    = 0.;
    double Qstar = 0.;

    // Calculate minimum and maximum HLLE speeds
    void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
                    const dTensor1& Ql, const dTensor1& Qr,
                    const dTensor1& Auxl, const dTensor1& Auxr,
                    double& s1, double& s2);
    SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

    int mcase = -1;

    // Calculate Fluxes
    if (fabs(s1)<=1.0e-12 && fabs(s2)<=1.0e-12)
    { mcase = 1; }
    else
    {
        if (s1*s2 > 0.0)
        {
            if (s1>0)
            { mcase = 1; }
            else
            { mcase = 2; }
        }
        else
        { mcase = 3; }
    }

    switch ( mcase )
    {
        case 1:  // both s1 and s2 are positive (or both zero)

            for (int m=1; m<=meqn; m++)
            { 
                Fl.set(m, ffl.get(m) );
                Fr.set(m, ffl.get(m) );
            }	
            break;

        case 2:  // both s1 and s2 are negative

            for (int m=1; m<=meqn; m++)
            { 
                Fl.set(m, ffr.get(m) );
                Fr.set(m, ffr.get(m) );
            }	
            break;

        case 3:  // s1 is negative and s2 is positive

            for (int m=1; m<=meqn; m++)
            {
//printf("you should never get here\n");
                //cout<<"s1= "<<s1<<" "<<s2<<endl;
                Qstar = (ffr.get(m)-ffl.get(m) 
                        + s1*Ql.get(m) - s2*Qr.get(m))/(s1-s2);
                Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
                            + (s1+s2)*Qstar - s1*Ql.get(m) - s2*Qr.get(m)) );
                Fr.set(m, Fl.get(m) );
            }	
            break;
        default:
            printf("you did not find a valid case\n");
            exit(1);
    }

    double smax_edge = Max(fabs(s1), fabs(s2));
    return smax_edge;
}
