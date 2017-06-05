#include "tensors.h"
#include "dog_math.h"
#include <cmath>
#include "RiemannSolve.h"
#include "dogdefs.h"
#include <iostream>
using namespace std;

// to use this RiemannSolver the user must implement the following
// "linking callbacks":
// * FluxFunc
// * SetWaveSpd
//
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
// Conservation (Rankine-Hugoniot) at the shock fronts
// gives a system of two equations which can be solved
// for the two unknowns Q* and F*:
//
// left  shock: (Fr - F*) = (Qr - Q*) s2
// right shock: (F* - Fl) = (Q* - Ql) s1
//
// Solving gives:
//
//   Q* = (Fr - Fl + s1 Ql - s2 Qr)/(s1-s2)  and
//   F* = (Fl + Fr - s1 Ql - s2 Qr + (s1+s2)Q*)/2.
//

double RiemannSolver::solve(const dTensor1& nvec,
                            const dTensor1& xedge,
                            const dTensor1& Ql,
                            const dTensor1& Qr,
                            const dTensor1& Auxl,
                            const dTensor1& Auxr,
                            dTensor1& Fl,
                            dTensor1& Fr)
{

    // Retrieve tensors local to this Riemann solver object:
    dTensor1& ffll       = fetch_ffl  ();
    dTensor1& ffrr       = fetch_ffr  ();
    dTensor2& xedge_tmp = fetch_xedge();
    dTensor2& Ql_tmp    = fetch_Ql   ();
    dTensor2& Qr_tmp    = fetch_Qr   ();
    dTensor2& Auxl_tmp  = fetch_Auxl ();
    dTensor2& Auxr_tmp  = fetch_Auxr ();
    dTensor1& fflr(ffll); 
    dTensor1& ffrl(ffrr); 

    // 2D problems require the "2D" flux function, which has indices,
    //    ftmp( 1:numpts, 1:meqn, 1:2 )
    dTensor3& ftmp      = fetch_flux2 ();  

    int mcase;
    double smax_edge = 0.0e0;
    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();
    double s1,s2,Qstar;

    // Reshape Ql,Qr,Auxl,Auxr,xedge to conform with FluxFunc.cpp
    //
    // Note that for a Riemann problem, mpoints = 1 always, which is the first
    // index in the tensors that get passed into FluxFunc.
    //
    for (int m=1; m<=meqn; m++)
    {  
        Ql_tmp.set(1, m, Ql.get(m) );  
        Qr_tmp.set(1, m, Qr.get(m) );
    }
    for (int m=1; m<=maux; m++)
    {  
        Auxl_tmp.set(1, m, Auxl.get(m) );  
        Auxr_tmp.set(1, m, Auxr.get(m) );
    }
    for (int m=1; m<=2; m++)
    {  xedge_tmp.set(1, m, xedge.get(m) );  }

    // Evaluate flux function at Ql and Qr
    void FluxFunc(const dTensor2& xpts,
                  const dTensor2& Q,
                  const dTensor2& Aux,
                  dTensor3& flux);
    FluxFunc(xedge_tmp,Ql_tmp,Auxl_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffll.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) );  }

    FluxFunc(xedge_tmp,Ql_tmp,Auxr_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  fflr.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) );  }

    FluxFunc(xedge_tmp,Qr_tmp,Auxl_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffrl.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) );  }

    FluxFunc(xedge_tmp,Qr_tmp,Auxr_tmp,ftmp);
    for (int m=1; m<=meqn; m++)
    {  ffrr.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) );  }

    // Calculate minimum and maximum HLLE speeds
    void SetWaveSpd(const dTensor1& nvec,
                    const dTensor1& xedge,
                    const dTensor1& Ql,
                    const dTensor1& Qr,
                    const dTensor1& Auxl,
                    const dTensor1& Auxr,
                    double& s1,
                    double& s2);
    SetWaveSpd(nvec, xedge, Ql, Qr, Auxl, Auxr, s1, s2);

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
                Fl.set(m, ffll.get(m) );
                Fr.set(m, fflr.get(m) );
            }	
            break;

        case 2:  // both s1 and s2 are negative

            for (int m=1; m<=meqn; m++)
            { 
                Fl.set(m, ffrl.get(m) );
                Fr.set(m, ffrr.get(m) );
            }	
            break;

        case 3:  // s1 is negative and s2 is positive

            for (int m=1; m<=meqn; m++)
            {
                Qstar = (ffrr.get(m)-ffll.get(m) 
                        + s1*Ql.get(m) - s2*Qr.get(m))/(s1-s2);
                Fl.set(m, 0.5e0*(ffll.get(m) + ffrr.get(m) 
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
double RiemannSolve(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr)
{
    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();
    RiemannSolver riemannSolver(meqn,maux);
    return riemannSolver.solve(nvec, xedge, Ql, Qr, Auxl, Auxr, Fl, Fr);
}
