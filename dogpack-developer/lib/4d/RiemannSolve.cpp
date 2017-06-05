#include "tensors.h"
#include "dog_math.h"
#include <cmath>
#include "RiemannSolve.h"

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
			    const dTensor1& xface,
			    const dTensor1& Ql, 
			    const dTensor1& Qr,
			    const dTensor1& Auxl, 
			    const dTensor1& Auxr,
			    dTensor1& Fl, 
			    dTensor1& Fr)
{
  dTensor1& ffl   = fetch_ffl  ();
  dTensor1& ffr   = fetch_ffr  ();
  dTensor2& xface_tmp = fetch_xface();
  dTensor2& Ql_tmp    = fetch_Ql   ();
  dTensor2& Qr_tmp    = fetch_Qr   ();
  dTensor2& Auxl_tmp  = fetch_Auxl ();
  dTensor2& Auxr_tmp  = fetch_Auxr ();
  dTensor3& ftmp  = fetch_flux3 ();
  
  int m,mcase;
  double smax_edge = 0.0e0;
  int meqn = Ql.getsize();
  int maux = Auxl.getsize();
  double s1,s2,Qstar;
  
  // Reshape Ql,Qr,Auxl,Auxr,xface to conform with FluxFunc.cpp
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
  for (m=1; m<=3; m++)
    {  xface_tmp.set(1,m, xface.get(m) );  }
  
  // Evaluate flux function at Ql and Qr
  void FluxFunc(const dTensor2& xpts, const dTensor2& Q,
		const dTensor2& Aux, dTensor3& flux);
  FluxFunc(xface_tmp,Ql_tmp,Auxl_tmp,ftmp);
  for (m=1; m<=meqn; m++)
    {  ffl.set(m, ftmp.get(1,m,1)*nvec.get(1) 
	       + ftmp.get(1,m,2)*nvec.get(2) 
	       + ftmp.get(1,m,3)*nvec.get(3) );  }
  
  FluxFunc(xface_tmp,Qr_tmp,Auxr_tmp,ftmp);
  for (m=1; m<=meqn; m++)
    {  ffr.set(m, ftmp.get(1,m,1)*nvec.get(1) 
	       + ftmp.get(1,m,2)*nvec.get(2) 
	       + ftmp.get(1,m,3)*nvec.get(3) );  }
  
  // Calculate minimum and maximum HLLE speeds
  void SetWaveSpd(const dTensor1& nvec, const dTensor1& xface,
		  const dTensor1& Ql, const dTensor1& Qr,
		  const dTensor1& Auxl, const dTensor1& Auxr,
		  double& s1,double& s2);
  SetWaveSpd(nvec,xface,Ql,Qr,Auxl,Auxr,s1,s2);
  
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
