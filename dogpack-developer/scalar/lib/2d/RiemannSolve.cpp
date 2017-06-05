#include "tensors.h"
#include "dog_math.h"
#include <cmath>
#include "RiemannSolve.h"
#include <stdio.h>

double RiemannSolver::solve(const dTensor2* vel_vec, 
        const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr)
{

    dTensor1& ffl   = fetch_ffl  ();
    dTensor1& ffr   = fetch_ffr  ();
    dTensor2& xedge_tmp = fetch_xedge();
    dTensor2& Ql_tmp    = fetch_Ql   ();
    dTensor2& Qr_tmp    = fetch_Qr   ();
    dTensor2& Auxl_tmp  = fetch_Auxl ();
    dTensor2& Auxr_tmp  = fetch_Auxr ();
    dTensor3& ftmp  = fetch_flux2 ();

    int m,mcase;
    double smax_edge = 0.0e0;
    int meqn = Ql.getsize();
    int maux = Auxl.getsize();
    double s1,s2,Qstar;

    // Reshape Ql,Qr,Auxl,Auxr,xedge to conform with FluxFunc.cpp
    for (m=1; m<=meqn; m++)
    {  
        Ql_tmp.set(1,m, Ql.get(m) );  
        Qr_tmp.set(1,m, Qr.get(m) );
    }
    for (m=1; m<=maux; m++)
    {  
        Auxl_tmp.set(1,m, Auxl.get(m) );  
        Auxr_tmp.set(1,m, Auxr.get(m) );
    }
    for (m=1; m<=2; m++)
    {  xedge_tmp.set(1,m, xedge.get(m) );  }

    // Evaluate flux function at Ql and Qr
    void FluxFunc(
            const dTensor2* vel_vec,
            const dTensor2& xpts, const dTensor2& Q,
            const dTensor2& Aux, dTensor3& flux);
    FluxFunc(vel_vec, xedge_tmp,Ql_tmp,Auxl_tmp,ftmp);
    for (m=1; m<=meqn; m++)
    {  ffl.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) );  }

    FluxFunc(vel_vec, xedge_tmp,Qr_tmp,Auxr_tmp,ftmp);
    for (m=1; m<=meqn; m++)
    {  ffr.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) );  }

    // Calculate minimum and maximum HLLE speeds
    void SetWaveSpd(const dTensor2* vel_vec,
            const dTensor1& nvec, const dTensor1& xedge,
            const dTensor1& Ql, const dTensor1& Qr,
            const dTensor1& Auxl, const dTensor1& Auxr,
            double& s1,double& s2);
    SetWaveSpd(vel_vec, nvec,xedge,Ql,Qr,Auxl,Auxr,s1,s2);

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

            for (m=1; m<=meqn; m++)
            { 
                Fl.set(m, ffl.get(m) );
                Fr.set(m, ffl.get(m) );
            }	
            break;

        case 2:  // both s1 and s2 are negative

            for (m=1; m<=meqn; m++)
            { 
                Fl.set(m, ffr.get(m) );
                Fr.set(m, ffr.get(m) );
            }	
            break;

        case 3:  // s1 is negative and s2 is positive

            for (m=1; m<=meqn; m++)
            {
                Qstar = (ffr.get(m)-ffl.get(m) 
                        + s1*Ql.get(m) - s2*Qr.get(m))/(s1-s2);
                Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
                            + (s1+s2)*Qstar - s1*Ql.get(m) - s2*Qr.get(m)) );
                Fr.set(m, Fl.get(m) );
            }	
            break;
    }

    smax_edge = Max(fabs(s1),fabs(s2));
    return smax_edge;
}

// support the old RiemannSolve interface for backwards compatibility
//
double RiemannSolve(const dTensor2* vel_vec, const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr)
{
    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();
    RiemannSolver riemannSolver(meqn,maux);
    return riemannSolver.solve(vel_vec, nvec, xedge, Ql, Qr, Auxl, Auxr, Fl, Fr);
}
