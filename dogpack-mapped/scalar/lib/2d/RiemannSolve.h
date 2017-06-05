#ifndef RiemannSolve_h
#define RiemannSolve_h
#include "tensors.h"

class RiemannSolver
{
    dTensor1* ffl;
    dTensor1* ffr;
    dTensor2* xedge;
    dTensor2* Ql;
    dTensor2* Qr;
    dTensor2* Auxl;
    dTensor2* Auxr;
    dTensor2* flux1;
    dTensor3* flux2;

    public:
    RiemannSolver(int meqn, int maux)
    {
        ffl   = new dTensor1(meqn);
        ffr   = new dTensor1(meqn);
        xedge = new dTensor2(1,2);
        Ql    = new dTensor2(1,meqn);
        Qr    = new dTensor2(1,meqn);
        Auxl  = new dTensor2(1,maux);
        Auxr  = new dTensor2(1,maux);
        flux1 = new dTensor2(1,meqn);
        flux2 = new dTensor3(1,meqn,2);
    }
    ~RiemannSolver()
    {
        delete ffl;
        delete ffr;
        delete xedge;
        delete Ql;
        delete Qr;
        delete Auxl;
        delete Auxr;
        delete flux1;
        delete flux2;
    }
    double solve(const dTensor2* vel_vec, const dTensor1& nvec, const dTensor1& xedge,
            const dTensor1& Ql, const dTensor1& Qr,
            const dTensor1& Auxl, const dTensor1& Auxr,
            dTensor1& Fl, dTensor1& Fr);
    public:
    dTensor1& fetch_ffl  (){return *ffl  ;}
    dTensor1& fetch_ffr  (){return *ffr  ;}
    dTensor2& fetch_xedge(){return *xedge;}
    dTensor2& fetch_Ql   (){return *Ql   ;}
    dTensor2& fetch_Qr   (){return *Qr   ;}
    dTensor2& fetch_Auxl (){return *Auxl ;}
    dTensor2& fetch_Auxr (){return *Auxr ;}
    dTensor2& fetch_flux1(){return *flux1;}
    dTensor3& fetch_flux2(){return *flux2;}
};

double RiemannSolve(const dTensor2* vel_vec, 
    const dTensor1& nvec, const dTensor1& xedge,
    const dTensor1& Ql, const dTensor1& Qr,
    const dTensor1& Auxl, const dTensor1& Auxr,
    dTensor1& Fl, dTensor1& Fr);

#endif // RiemannSolve_h