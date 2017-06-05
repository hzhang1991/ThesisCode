#include <cmath>
#include "tensors.h"
#include "dog_math.h"

// HLL(E) Riemman solver option for Lax-Wendroff (and two-derivative) code.
//
// The multiderivative methods tuck information about time derivatives of the
// function into the Riemann solver.  This information has already been
// computed, so there's no need to recompute it at the cell interfaces.
//
// To use this solver, one needs to replace the linker in order to link to
// this.  Add the following to your Makefile:
//
//     RiemannSolveLxW  = $(LIB1D)/LaxWendroff/RiemannSolveLxW_HLL
//
// ( c.f. RiemannSolveLxW_LLF.cpp )
//
double RiemannSolveLxW(const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr, const dTensor1& Auxl,
        const dTensor1& Auxr, const dTensor1& ffl, const dTensor1& ffr,
        dTensor1& Fl, dTensor1& Fr,
        void (*SetWaveSpd)(const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,const dTensor1&,double&,double&))
{

    ///////////////////////////////////////////////////////////////////////////
    // HLL(E) Riemann Solver:
    ///////////////////////////////////////////////////////////////////////////

    const int meqn = Ql.getsize();

    // Calculate minimum and maximum HLLE speeds
    double s1,s2;
    SetWaveSpd(xedge,Ql,Qr,Auxl,Auxr,s1,s2);

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
                // Approximate single intermediate state value:
                double Qstar = (
                    ffr.get(m)-ffl.get(m) 
                    + s1*Ql.get(m) - s2*Qr.get(m) ) / (s1-s2);

                Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
                            + (s1+s2)*Qstar - s1*Ql.get(m) - s2*Qr.get(m)) );
                Fr.set(m, Fl.get(m) );
            }

            break;
    }

    double smax_edge = Max(fabs(s1),fabs(s2));
    return smax_edge;

}
