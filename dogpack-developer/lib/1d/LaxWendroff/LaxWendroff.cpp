#include <iostream>
#include <cmath>
#include "DogParams.h"
#include "tensors.h"
#include "dog_math.h"
#include "assert.h"

using namespace std;

// Right-hand side of the Lax-Wendroff discretization:
//
//      ( q(t+dt) - q(t) )/dt = -F_{x}
//
// where the the flux is given by
//
//      F = F1 + dt/2 * F2 + dt^2/6 * F3
//
//        F1 =  f(q)
//        F2 = -f'(q) * f(q)_{x}
//        F3 =  f''(q) * f(q)_{x} * f(q)_{x} - f'(q)*f(q)_{xx} 
//
// This is the equivalent of the ConstructL procedure, but in this case, we
// include dt.
// 
void LaxWendroff(const double dt, const int method[], const dTensor2& node,
        dTensorBC3& aux, dTensorBC3& q, 
        dTensorBC3& Lstar, dTensorBC1& smax)
{

    void EvalRiemmanData(const int& i, const dTensorBC3& q, const dTensorBC3& aux,  
            const dTensorBC3& LxWF, dTensor1& Ql, dTensor1& Qr, 
            dTensor1& LxFFluxl, dTensor1& LxFFluxr,
            dTensor1& Auxl, dTensor1& Auxr);

    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   kmax = q.getsize(3);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    assert_eq( kmax, 3 );
    dTensorBC2 Fm (melems, meqn, mbc);
    dTensorBC2 Fp (melems, meqn, mbc);
    dTensorBC3 N  (melems, meqn, kmax, mbc);
    dTensorBC3 Psi(melems, meqn, kmax, mbc);

    // Lax-Wendroff flux function
    dTensorBC3  LxWF(melems, meqn, kmax, mbc);

    // The Lax-Wendroff flux function broken up into terms involing
    // dt^0, dt^1, and dt^2.  These were previously used for 'stability
    // correction' terms that are no longer used anywhere in the code.
    // It is not necessary to break up the function like this, but at the time
    // made coding it up easier to do.  In the future, we should probably not
    // use as much storage as is used here.
    //
    dTensorBC3 LxWF1(melems, meqn, kmax, mbc);
    dTensorBC3 LxWF2(melems, meqn, kmax, mbc);
    dTensorBC3 LxWF3(melems, meqn, kmax, mbc);
    dTensorBC3 IntF1(melems, meqn, kmax, mbc);
    dTensorBC3 IntF2(melems, meqn, kmax, mbc);
    dTensorBC3 IntF3(melems, meqn, kmax, mbc);

    double RiemannSolveLxW(const dTensor1& xedge,
            const dTensor1& Ql, const dTensor1& Qr, const dTensor1& Auxl,
            const dTensor1& Auxr, const dTensor1& ffl, const dTensor1& ffr,
            dTensor1& Fl, dTensor1& Fr,
            void (*SetWaveSpd)(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
                const dTensor1&,double&,double&));
    void SetBndValues(const dTensor2&,dTensorBC3&,dTensorBC3&);
    void L2Project(int,int,int,const dTensor2&,const dTensorBC3&,
            const dTensorBC3&,dTensorBC3&,
            void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));
    void L2ProjectLxW(const int istart, const int iend, 
            const dTensor2& node, const dTensorBC3& qin, const dTensorBC3& auxin,  
            dTensorBC3& F1, dTensorBC3& F2, dTensorBC3& F3,
            dTensorBC3& LxWFlux1, dTensorBC3& LxWFlux2, dTensorBC3& LxWFlux3,
            void (*Func)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&),
            void (*DFunc)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor3&),
            void (*D2Func)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor4&));
    void FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void DFluxFunc (const dTensor1&,const dTensor2&,const dTensor2&,dTensor3&);
    void D2FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor4&);
    void SetWaveSpd(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,double&,double&);
    void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void LstarExtra(const dTensor2&,dTensorBC3&,dTensorBC3&,dTensorBC3&);

    // quick error check
    if( method[7] > 0 )
    {
        cout << "    error: have not implemented source term for LxW solver "
            << endl;
        exit(1);
    }

    // Grid spacing
    const double xlower = node.get(1,1);
    const double dx     = node.get(2,1) - node.get(1,1);

    // double check to make sure these were initialized to zero...
    LxWF1.setall(0.);
    LxWF2.setall(0.);
    LxWF3.setall(0.);

    IntF1.setall(0.);
    IntF2.setall(0.);
    IntF3.setall(0.);

    // Boundary conditions
    SetBndValues(node,aux,q);

    // Compute the Lax-Wendroff "flux" function:
    L2ProjectLxW(1-mbc, melems+mbc, node, q, aux,  
            IntF1, IntF2, IntF3,  LxWF1,  LxWF2,  LxWF3,
            &FluxFunc, &DFluxFunc, &D2FluxFunc);
#pragma omp parallel for
    for(int i=1-mbc; i<= melems+mbc; i++ )
    for(int k=1; k<= kmax; k++)	    
    for(int m=1; m<=meqn; m++)
    {
        double tmp = 
              LxWF1.get(i,m,k)
            + LxWF2.get(i,m,k) * 0.5 * dt
            + LxWF3.get(i,m,k) * pow(dt,2) / 6.0;
        LxWF.set(i,m,k, tmp );
    }

    // ---------------------------------------------------------
    // Part I: compute inter-element interaction fluxes
    // ---------------------------------------------------------

    // Boundary conditions
    SetBndValues(node,aux,q);

    // Loop over interior edges and solve Riemann problems
// TODO: this pragma statement is "buggy".  Probably because we're not setting
// the correct max wave speed.  At any rate, this doesn't agree with the
// serial code.  (-DS).
#pragma omp parallel for
    for (int i=(2-mbc); i<=(melems+mbc); i++)
    {

        // Local storage
        dTensor1   Ql(meqn),   Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);

        dTensor1 LxFFluxl(meqn),  LxFFluxr(meqn);
        dTensor1       Fl(meqn),        Fr(meqn);

        dTensor1 xedge(1);
        xedge.set(1, xlower + (double(i)-1.0)*dx );

        // Evaluate the Riemman data
        //
        // This access q(i) and q(i-1), as well as aux(i) and aux(i-1).
        // 
        // This routine simply uses what has already been computed to define
        // left and right hand states.
        //
        EvalRiemmanData(i, q, aux, LxWF, Ql, Qr, LxFFluxl, LxFFluxr, Auxl, Auxr); 
        double smax_edge = RiemannSolveLxW(xedge,Ql,Qr,Auxl,Auxr,
                LxFFluxl, LxFFluxr, Fl, Fr, &SetWaveSpd);

        // TODO: this *should* be a problem if the openmp flags are turned
        // on.
        // 
        // For some reason, we don't see this as a problem in the 2D code ...
        // (-DS)
        //
        smax.set(i-1, Max(smax_edge,smax.get(i-1)) );
        smax.set(i,   Max(smax_edge,smax.get(i)) );

        // Construct fluxes
        for (int m=1; m<=meqn; m++)
        {
            Fm.set(i,  m, Fr.get(m) );
            Fp.set(i-1,m, Fl.get(m) );
        }

    }
    SetBndValues(node,aux,q);  // clean up boundary values
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part II: compute intra-element contributions
    // ---------------------------------------------------------
    //
    //   N = int( F(q,x,t) * phi_x, x )/dx
    //
    // Compute ``N'' by projecting flux function onto the 
    // gradient of Legendre polynomials
    if (method[1]>1)
    {
#pragma omp parallel for
        for(int i=1-mbc; i<=(melems+mbc); i++)           
        for(int k=1; k<=kmax; k++)
        for(int m=1; m<=meqn; m++)
        {
            double tmp = IntF1.get(i,m,k) 
                + IntF2.get(i,m,k) * 0.5*dt
                + IntF3.get(i,m,k) * pow(dt,2)/6.0;
            N.set(i,m,k, tmp );
        }
    }
    else // case kmax == 1 - there is no integral term in this case
    {
        N.setall(0.);
    }
    SetBndValues(node,aux,q);   // clean up boundary values
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part III: compute source term
    // --------------------------------------------------------- 
    if ( method[7]>0 )
    // TODO this needs some work - source terms need to be expanded in taylor series as
    // well in order to attain high order:
    {   
        cout << "    Warning: LxW solver needs to expand source term." << endl;
        cout << " method will be first order for source term function" << endl;

        // Set source term on computational grid
        // Set values and apply L2-projection
        L2Project(0,1-mbc,melems+mbc,node,q,aux,Psi,&SourceTermFunc);
    }
    // ---------------------------------------------------------

    // Apply MPP Flux-limiter, if requested
    if( dogParams.using_positive_limiter() )
    {
        void ApplyPosMPPLimiter(const double dt, 
            const int method[], const dTensor2& node,
            const dTensorBC1& smax,
            const dTensorBC3& aux, const dTensorBC3& q, 
            dTensorBC2& Fm, dTensorBC2& Fp );
        ApplyPosMPPLimiter( dt, method, node, smax, aux, q, Fm, Fp );
    }

    // ---------------------------------------------------------
    // Part IV: construct Lstar
    // ---------------------------------------------------------
    if (method[7]==0)  // Without Source Term
    { 
#pragma omp parallel for
        for (int i=(2-mbc); i<=(melems+mbc-1); i++)
        for (int m=1; m <= meqn; m++)
        for (int k=1; k <= kmax; k++)
        {
            double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx;
            Lstar.set(i,m,k, tmp );
        }
    }
    else  // With Source Term
    {
#pragma omp parallel for
        for (int i=(2-mbc); i<=(melems+mbc-1); i++)
        for (int m=1; m <=meqn; m++)
        for (int k=1; k <= kmax; k++)
        {
            double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx
                + Psi.get(i,m,k);

            Lstar.set(i,m,k, tmp );
        }
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part V: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    LstarExtra(node,aux,q,Lstar);
    // ---------------------------------------------------------

}
