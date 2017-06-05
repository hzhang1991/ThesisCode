#ifndef _LAXWENDROFF_UNST_H_
#define _LAXWENDROFF_UNST_H_

// Functions this routine needs to have access to:
//  double RiemannSolve(const dTensor1& nvec, const dTensor1& xedge,
//          const dTensor1& Ql, const dTensor1& Qr,
//          const dTensor1& Auxl, const dTensor1& Auxr,
//          dTensor1& Fl, dTensor1& Fr);

// Riemann solver that relies on the fact that we already have 
// f(ql) and f(qr) already computed:
double RiemannSolveLxW(const dTensor1& nvec,
    const dTensor1& xedge,
    const dTensor1& Ql,   const dTensor1& Qr,
    const dTensor1& Auxl, const dTensor1& Auxr,
    const dTensor1& ffl,  const dTensor1& ffr,
    dTensor1& Fl, dTensor1& Fr);

void SetBndValues_Unst(const mesh&,dTensor3*,dTensor3*);


void L2ProjectLxW_Unst( const int mterms,
        const double alpha, const double beta_dt, const double charlie_dt,
        const int istart, const int iend,               // Start-stop indices
        const int QuadOrder,
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,
        const mesh& Mesh, 
        const dTensor3* qin, const dTensor3* auxin,     // state vector
        dTensor3* F, dTensor3* G,                       // time-averaged Flux function
        void FluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& Aux, dTensor3& flux),
        void DFluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& aux, dTensor4& Dflux),
        void D2FluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& aux, dTensor5& D2flux) );

void L2ProjectGradAddLegendre_Unst(const int istart, const int iend, 
        const int QuadOrder, 
        const mesh& Mesh, 
        const dTensor3* F, 
        const dTensor3* G, 
        dTensor3* fout );

// Flux function and its (analytical) derivatives
void FluxFunc(const dTensor2& xpts, const dTensor2& Q,
        const dTensor2& Aux, dTensor3& flux);
void  DFluxFunc(const dTensor2& xpts, const dTensor2& Q, 
        const dTensor2& Aux, dTensor4& Dflux );
void D2FluxFunc(const dTensor2& xpts, const dTensor2& Q, 
        const dTensor2& Aux, dTensor5& D2flux );

void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, const dTensor2& auxvals,
        dTensor2& source);
void LstarExtra_Unst(const mesh&,const dTensor3*,const dTensor3*,dTensor3*);

#endif
