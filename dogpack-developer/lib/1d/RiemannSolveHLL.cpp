#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"

// This is the HLL(E) method.  The HLLE method defines left and right
// speeds by
//
// s1 := minimum speed of waves for left state and intermediate state and
// s2 := maximum speed of waves for right state and intermediate state,
//       where intermediate state is theoretically the Roe average of
//       the left and right states and in practice is often the simple
//       average.
//
// (The original HLL method defines s1 and s2 as extreme characteristic
// speeds of Riemann solutions.)
//
// In case s1 < 0 < s2 approximate the Riemann solution by two
// shockwaves of speeds s1 and s2 separating a single constant
// intermediate state Q* from the left state Ql and the right state Qr.
// Define Fl := F(Ql) and Fr := F(Qr).
//
// Conservation (Rankine-Hugoniot) at the shock fronts
// gives a system of two equations which can be solved
// for the two unknowns Q* and F*:
//
// left  shock: (Fr - F*) = (Qr - Q*) s2
// right shock: (F* - Fl) = (Q* - Ql) s1
//
// Solving this system gives:
//
//   Q* = (Fr - Fl + s1 Ql - s2 Qr)/(s1-s2)  
//
// and
//
//   F* = (Fl + Fr - s1 Ql - s2 Qr + (s1+s2)Q*)/2.
//
double RiemannSolve(const dTensor1& xedge,
        const dTensor1& Ql,
        const dTensor1& Qr,
        const dTensor1& Auxl,
        const dTensor1& Auxr,
        dTensor1& Fl,
        dTensor1& Fr,
        void (*FluxFunc)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&),
        void (*SetWaveSpd)(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,double&,double&))
{
    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();
    dTensor1 ffl(meqn), ffr(meqn);
    dTensor2 Ql_tmp(1,meqn),Qr_tmp(1,meqn);
    dTensor2 Auxl_tmp(1,maux),Auxr_tmp(1,maux);
    dTensor2 ftmp(1,meqn);

    // Reshape Ql,Qr,Auxl,Auxr,xedge to conform with FluxFunc.cpp
    for (int m=1; m<=meqn; m++)
    {  
        Ql_tmp.set(1,m, Ql.get(m) );  
        Qr_tmp.set(1,m, Qr.get(m) );
    }
    for (int m=1; m<=maux; m++)
    {  
        Auxl_tmp.set(1,m, Auxl.get(m) );  
        Auxr_tmp.set(1,m, Auxr.get(m) );
    }

    // Evaluate flux function at Ql and Qr
    FluxFunc(xedge,Ql_tmp,Auxl_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffl.set(m, ftmp.get(1,m) );  }

    FluxFunc(xedge,Qr_tmp,Auxr_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffr.set(m, ftmp.get(1,m) );  }

    // Calculate minimum and maximum HLLE speeds
    double s1,s2;
    SetWaveSpd(xedge, Ql, Qr, Auxl, Auxr, s1, s2);

    // Calculate Fluxes
    int mcase;
    if (fabs(s1)<=1.0e-12 && fabs(s2)<=1.0e-12)
    { mcase = 1; }
    else
    {
        if (s1*s2 > 1.0e-12)
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
                double Qstar = (ffr.get(m)-ffl.get(m) + s1*Ql.get(m) 
                        - s2*Qr.get(m))/(s1-s2);
                Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
                            + (s1+s2)*Qstar - s1*Ql.get(m) - s2*Qr.get(m)) );
                Fr.set(m, Fl.get(m) );
            }

            break;
    }

    double smax_edge = Max(fabs(s1),fabs(s2));
    return smax_edge;
}
