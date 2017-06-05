#include "tensors.h"
#include "dog_math.h"
#include <cmath>
#include "RiemannSolveGLF.h"
#include <iostream>
using namespace std;

// to use this RiemannSolver the user must implement the following
// "linking callbacks":
// * FluxFunc
// * SetWaveSpd

double RiemannSolverGLF::solve(const dTensor1& nvec, const dTensor1& xedge,
		    const dTensor1& Ql, const dTensor1& Qr,
		    const dTensor1& Auxl, const dTensor1& Auxr,
		    dTensor1& Fl, dTensor1& Fr,const double smax)
{
    int m,mcase;
    double smax_edge = 0.0e0;
    int meqn = Ql.getsize();
    int maux = Auxl.getsize();
    dTensor1 ffl(meqn), ffr(meqn);
    dTensor3 ftmp(1,meqn,2);
    double s1,s2,Qstar;
    dTensor2 xedge_tmp(1,2),Ql_tmp(1,meqn),Qr_tmp(1,meqn);
    dTensor2 Auxl_tmp(1,maux),Auxr_tmp(1,maux);

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
    
    void FluxFunc(const dTensor2& xpts,
                  const dTensor2& Q,
                  const dTensor2& Aux,
                  dTensor3& flux);
    // Evaluate flux function at Ql and Qr
    FluxFunc(xedge_tmp,Ql_tmp,Auxl_tmp,ftmp);
    for (m=1; m<=meqn; m++)
      {  ffl.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) );  }
    
    FluxFunc(xedge_tmp,Qr_tmp,Auxr_tmp,ftmp);
    for (m=1; m<=meqn; m++)
      {  ffr.set(m, ftmp.get(1,m,1)*nvec.get(1) + ftmp.get(1,m,2)*nvec.get(2) );  }
    
    void SetWaveSpd(const dTensor1& nvec, 
                    const dTensor1& xedge,
                    const dTensor1& Ql, 
                    const dTensor1& Qr,
                    const dTensor1& Auxl, 
                    const dTensor1& Auxr,
                    double& s1,double& s2);
    // Calculate minimum and maximum HLLE speeds
    SetWaveSpd(nvec,xedge,Ql,Qr,Auxl,Auxr,s1,s2);

    // Local Lax-Friedrichs method
    smax_edge = smax;//Max(fabs(s1),fabs(s2));
    s1 = -smax_edge;
    s2 =  smax_edge;

    for (m=1; m<=meqn; m++)
      {
	Fl.set(m, 0.5e0*(ffl.get(m) + ffr.get(m) 
			 - smax_edge*(Qr.get(m)-Ql.get(m))) );
	Fr.set(m, Fl.get(m) );
      }    
    
    return smax_edge;
}

// support the old RiemannSolve interface for backwards compatibility
//
double RiemannSolveGLF(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr,const double smax)
{
  const int meqn = Ql.getsize();
  const int maux = Auxl.getsize();
  RiemannSolverGLF riemannSolver(meqn,maux);
  return riemannSolver.solve(nvec, xedge, Ql, Qr, Auxl, Auxr, Fl, Fr,smax);
}
