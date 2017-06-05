#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"

// Local Lax-Friedrichs solver
//
// This is the default Riemann solver for all of the 1D problems.  To use a
// different solver, one needs to link to RiemannSolveHLL in place of this
// one.
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
    double smax_edge = 0.0e0;
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
    SetWaveSpd(xedge,Ql,Qr,Auxl,Auxr,s1,s2);
    smax_edge = Max(fabs(s1),fabs(s2));

    // Calculate Fluxes
    for (int m=1; m<=meqn; m++)
    {
        Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
                    + smax_edge*(Ql.get(m) - Qr.get(m)) ) );
        Fr.set(m, Fl.get(m) );
    }

    return smax_edge;

}
