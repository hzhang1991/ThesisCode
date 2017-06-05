#ifndef _CONSTRUCTL_UNST_H_
#define _CONSTRUCTL_UNST_H_

// Functions this routine needs to have access to:
double RiemannSolve(const dTensor2* vel_vec, 
    const dTensor1& nvec, const dTensor1& xedge,
        const dTensor1& Ql, const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        dTensor1& Fl, dTensor1& Fr);
void SetBndValues_Unst(const mesh&,dTensor3*,dTensor3*);

void L2Project_Unst(
    const double t,
    const dTensor2* vel_vec,
    const int istart, 
    const int iend, 
    const int QuadOrder, 
    const int BasisOrder_qin,
    const int BasisOrder_auxin,
    const int BasisOrder_fout,
    const mesh& Mesh, 
    const dTensor3* qin, 
    const dTensor3* auxin, 
    dTensor3* fout, 
    void (*Func)(const double t, const dTensor2* vel_vec,
        const dTensor2&,const dTensor2&,
            const dTensor2&,dTensor2&));

void L2ProjectGrad_Unst(
    const dTensor2* vel_vec,
    const int istart, 
    const int iend, 
    const int QuadOrder, 
    const int BasisOrder_qin,
    const int BasisOrder_auxin,
    const int BasisOrder_fout,
    const mesh& Mesh, 
    const dTensor3* qin, 
    const dTensor3* auxin, 
    dTensor3* fout, 
    void (*Func)(const dTensor2* vel_vec,
            const dTensor2&,const dTensor2&,
            const dTensor2&,dTensor3&));

void LstarExtra_Unst(const mesh&,const dTensor3*,const dTensor3*,dTensor3*);

// User supplied functions
void FluxFunc(const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& Q,
    const dTensor2& Aux, dTensor3& flux);

void SourceTermFunc(
    const double t, 
    const dTensor2* vel_vec,
    const dTensor2& xpts, const dTensor2& qvals, const dTensor2& auxvals,
        dTensor2& source);

#endif
