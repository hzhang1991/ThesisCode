#ifndef RiemannSolve_h
#define RiemannSolve_h
#include "tensors.h"

void FluxFunc(const dTensor2& xpts,
              const dTensor2& Q,
              const dTensor2& Aux,
              dTensor3& flux);

void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
                const dTensor1& Ql, const dTensor1& Qr,
                const dTensor1& Auxl, const dTensor1& Auxr,
                double& s1, double& s2);

// -------------------------------------------------------------------------- //
class RiemannSolverGLF
{

    // Left and right hand values of flux function.
    // ffl = ffl( 1:MEQN ) and ffr = ffr( 1:MEQN )
    dTensor1* ffl;
    dTensor1* ffr;

    // The indices for each of these tensors have been chosen to conform with
    // the user supplied function, FluxFunc.
    //
    // That is, we have the following indices:
    //
    //    xedge( 1:numpts, 1:NDIMS )
    //       Ql( 1:numpts, 1:MEQN  )
    //       Qr( 1:numpts, 1:MEQN  )
    //     Auxl( 1:numpts, 1:MAUX  )
    //     Auxr( 1:numpts, 1:MAUX  )
    //
    // Note that numpts = 1 for any Riemann problem.
    //
    dTensor2* xedge;
    dTensor2* Ql;
    dTensor2* Qr;
    dTensor2* Auxl;
    dTensor2* Auxr;

    // These fluxes can be used for 1D, and 2D Riemann solvers.
    // This class has yet to be pushed one level higher, although given the
    // format of DoGPack, it would make sense to make all codes call this
    // routine for solving their Riemann problems.
    //
    // 1D flux function in q_t = f(q)_x :
    //    flux1( 1:numpts, 1:MEQN          )
    //
    // 2D flux function in q_t = f(q)_x + g(q)_y :
    //
    //    flux2( 1:numpts, 1:MEQN, 1:NDIMS )
    //
    // In 2D, flux2( :, :, 1 ) = f(q), and flux2( :, :, 2 ) = g(q)
    //
    dTensor2* flux1;
    dTensor3* flux2;

    public:
    RiemannSolverGLF(int meqn, int maux)
    {
        ffl   = new dTensor1(meqn);
        ffr   = new dTensor1(meqn);
        xedge = new dTensor2(1, 2);
        Ql    = new dTensor2(1, meqn);
        Qr    = new dTensor2(1, meqn);
        Auxl  = new dTensor2(1, maux);
        Auxr  = new dTensor2(1, maux);
        flux1 = new dTensor2(1, meqn);
        flux2 = new dTensor3(1, meqn, 2);
    }
    ~RiemannSolverGLF()
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
    double solve(const dTensor1& nvec, const dTensor1& xedge,
            const dTensor1& Ql, const dTensor1& Qr,
            const dTensor1& Auxl, const dTensor1& Auxr,
            dTensor1& Fl, dTensor1& Fr,const double);
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
// -------------------------------------------------------------------------- //


// -------------------------------------------------------------------------- //
//
// Single function call used to provide backwards compatability
//
// -------------------------------------------------------------------------- //
double RiemannSolveGLF(const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr);
// -------------------------------------------------------------------------- //

#endif // RiemannSolve_h
