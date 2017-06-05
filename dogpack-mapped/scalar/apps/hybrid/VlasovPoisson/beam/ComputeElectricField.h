#ifndef _COMPUTE_ELECTRIC_FIELD_H__
#define _COMPUTE_ELECTRIC_FIELD_H__

void ComputeEfield(const int space_order, const mesh& Mesh, const dTensor1& phi,
    dTensor2& E1, dTensor2& E2);

// Used for the applied electric field
void L2Project_Unst(
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
    void (*Func)(const dTensor2* vel_vec, const dTensor2&,const dTensor2&,
        const dTensor2&,dTensor2&));

#endif
