#include "tensors.h"
#include "dog_math.h"
#include <cmath>

// This Riemann Solver is more efficient than the solver defined
// in RiemannSolve.cpp in the case of meshes with (many)
// coordinate-aligned edges (e.g. Cartesian meshes), but it
// requires the user to implement FluxFunc1 and FluxFunc2 in
// addition to FluxFunc and SetWaveSpd.
//
double RiemannSolve(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr)
{
    const int meqn = Ql.getsize();
    const int maux = Auxl.getsize();
    dTensor1 ffl(meqn), ffr(meqn);
    dTensor2 xedge_tmp(1,2),Ql_tmp(1,meqn),Qr_tmp(1,meqn);
    dTensor2 Auxl_tmp(1,maux),Auxr_tmp(1,maux);

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
    for (int m=1; m<=2; m++)
    {  xedge_tmp.set(1,m, xedge.get(m) );  }

    // Evaluate flux function at Ql and Qr
    const double n1 = nvec.get(1);
    const double n2 = nvec.get(2);
    assert_almost_eq(n1*n1+n2*n2,1.);
    // check for the cheap cases first
    if(n2==0.)
    {
        dTensor2 ftmp(1,meqn);

        void FluxFunc1(
            const dTensor2& xpts,
            const dTensor2& Q,
            const dTensor2& Aux,
            dTensor2& flux);
        FluxFunc1(xedge_tmp,Ql_tmp,Auxl_tmp,ftmp);
        for (int m=1; m<=meqn; m++)
        {  ffl.set(m, ftmp.get(1,m) );  }

        FluxFunc1(xedge_tmp,Qr_tmp,Auxr_tmp,ftmp);
        for (int m=1; m<=meqn; m++)
        {  ffr.set(m, ftmp.get(1,m) );  }
    }
    else if(n1==0.)
    {
        dTensor2 ftmp(1,meqn);

        void FluxFunc2(
            const dTensor2& xpts,
            const dTensor2& Q,
            const dTensor2& Aux,
            dTensor2& flux);
        FluxFunc2(xedge_tmp,Ql_tmp,Auxl_tmp,ftmp);
        for (int m=1; m<=meqn; m++)
        {  ffl.set(m, ftmp.get(1,m) );  }

        FluxFunc2(xedge_tmp,Qr_tmp,Auxr_tmp,ftmp);
        for (int m=1; m<=meqn; m++)
        {  ffr.set(m, ftmp.get(1,m) );  }
    }
    else
    {
        dTensor3 ftmp(1,meqn,2);

        void FluxFunc(const dTensor2& xpts, const dTensor2& Q,
            const dTensor2& Aux, dTensor3& flux);
        FluxFunc(xedge_tmp,Ql_tmp,Auxl_tmp,ftmp);
        for (int m=1; m<=meqn; m++)
        {  ffl.set(m, ftmp.get(1,m,1)*n1 + ftmp.get(1,m,2)*n2 );  }

        FluxFunc(xedge_tmp,Qr_tmp,Auxr_tmp,ftmp);
        for (int m=1; m<=meqn; m++)
        {  ffr.set(m, ftmp.get(1,m,1)*n1 + ftmp.get(1,m,2)*n2 );  }
    }

    // Calculate minimum and maximum HLLE speeds
    void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        double& s1,double& s2);
    double s1,s2;
    SetWaveSpd(nvec,xedge,Ql,Qr,Auxl,Auxr,s1,s2);

    // Calculate Fluxes
    int mcase = 0;
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
                const double Qstar = (ffr.get(m)-ffl.get(m) 
                        + s1*Ql.get(m) - s2*Qr.get(m))/(s1-s2);
                Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
                            + (s1+s2)*Qstar - s1*Ql.get(m) - s2*Qr.get(m)) );
                Fr.set(m, Fl.get(m) );
            }	
            break;
    }

    double smax_edge = Max(fabs(s1),fabs(s2));
    return smax_edge;
}
