#ifndef RiemannSolve_h
#define RiemannSolve_h
#include "tensors.h"

class RiemannSolver
{
  dTensor1* ffl;   // flux function evaluated at left state
  dTensor1* ffr;   // flux function evaluated at right state
  dTensor2* xedge; // spatial location of Riemann solve (2D)
  dTensor2* xface; // spatial location of Riemann solve (3D)
  dTensor2* Ql;    // left state 
  dTensor2* Qr;    // right state
  dTensor2* Auxl;  // left aux array
  dTensor2* Auxr;  // right aux array
  dTensor2* flux1; 
  dTensor3* flux2;
  dTensor3* flux3;
  
 public:
  RiemannSolver(int meqn, int maux)
    {
      ffl   = new dTensor1(meqn);
      ffr   = new dTensor1(meqn);
      xedge = new dTensor2(1,2);
      xface = new dTensor2(1,3);
      Ql    = new dTensor2(1,meqn);
      Qr    = new dTensor2(1,meqn);
      Auxl  = new dTensor2(1,maux);
      Auxr  = new dTensor2(1,maux);
      flux1 = new dTensor2(1,meqn);
      flux2 = new dTensor3(1,meqn,2);
      flux3 = new dTensor3(1,meqn,3);
    }
  ~RiemannSolver()
    {
      delete ffl;
      delete ffr;
      delete xedge;
      delete xface;
      delete Ql;
      delete Qr;
      delete Auxl;
      delete Auxr;
      delete flux1;
      delete flux2;
      delete flux3;
    }
  double solve(const dTensor1& nvec, 
	       const dTensor1& xedge,
	       const dTensor1& Ql, 
	       const dTensor1& Qr,
	       const dTensor1& Auxl, 
	       const dTensor1& Auxr,
	       dTensor1& Fl, 
	       dTensor1& Fr);
 public:
  dTensor1& fetch_ffl  (){return *ffl  ;}
  dTensor1& fetch_ffr  (){return *ffr  ;}
  dTensor2& fetch_xedge(){return *xedge;}
  dTensor2& fetch_xface(){return *xface;}
  dTensor2& fetch_Ql   (){return *Ql   ;}
  dTensor2& fetch_Qr   (){return *Qr   ;}
  dTensor2& fetch_Auxl (){return *Auxl ;}
  dTensor2& fetch_Auxr (){return *Auxr ;}
  dTensor2& fetch_flux1(){return *flux1;}
  dTensor3& fetch_flux2(){return *flux2;}
  dTensor3& fetch_flux3(){return *flux3;}
};

#endif // RiemannSolve_h
