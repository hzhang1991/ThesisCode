#ifndef _LAX_WENDROFF_H__
#define _LAX_WENDROFF_H__

// Flux function
void FluxFunc(const dTensor2& xpts,
              const dTensor2& Q,
              const dTensor2& Aux,
              dTensor3& flux);

// Jacobian of the flux function:
void DFluxFunc(const dTensor2& xpts, 
               const dTensor2& Q,
               const dTensor2& Aux, 
               dTensor4& Dflux );

// Hessian of the flux function:
void D2FluxFunc(const dTensor2& xpts, 
               const dTensor2& Q,
               const dTensor2& Aux, 
               dTensor5& D2flux );

// Riemann solver that relies on the fact that we already have 
// f(ql) and f(qr) already computed:
//
// This is not the case for the generic Riemann solver, but for the
// Lax-Wendroff code, we compute time-averaged values from information
// obtained on the interior of the cell
//
double RiemannSolveLxW(const dTensor1& nvec,
    const dTensor1& xedge,
    const dTensor1& Ql,   const dTensor1& Qr,
    const dTensor1& Auxl, const dTensor1& Auxr,
    const dTensor1& ffl,  const dTensor1& ffr,
    dTensor1& Fl, dTensor1& Fr);

void LstarExtra(const dTensorBC4*,
        const dTensorBC4*,
        dTensorBC4*);
void ArtificialViscosity(const dTensorBC4* aux, 
        const dTensorBC4* q, 
        dTensorBC4* Lstar);

// This is the primary workhorse of the routine: it computes time-averages
// flux values, using only information from the inside of each element.
void L2ProjectLxW( const int mterms, 
    const double alpha, const double beta_dt, const double charlie_dt,
    const int istart, const int iend, 
    const int jstart, const int jend,
    const int QuadOrder,
    const int BasisOrder_auxin,
    const int BasisOrder_fout,    
    const dTensorBC4* qin,
    const dTensorBC4* auxin,
    dTensorBC4* F,
    dTensorBC4* G,
    void FluxFunc (const dTensor2& xpts, 
        const dTensor2& Q, const dTensor2& Aux, dTensor3& flux),
    void DFluxFunc (const dTensor2& xpts, 
        const dTensor2& Q, const dTensor2& aux, dTensor4& Dflux),
    void D2FluxFunc (const dTensor2& xpts, 
        const dTensor2& Q, const dTensor2& aux, dTensor5& D2flux) );

void L2ProjectGradAddLegendre(
    const int istart, const int iend, 
    const int jstart, const int jend,
    const int QuadOrder, const dTensorBC4* F, const dTensorBC4* G, 
    dTensorBC4* fout );

// MPP limiter (vacuous call unless an application overrides this)
void ApplyPosMPPLimiter( double dt, const dTensorBC3& smax, 
    const dTensorBC3& aux, const dTensorBC3& q, 
    dTensorBC3& FmI, dTensorBC3& GmI );

#endif
