#ifndef L2Project_h
#define L2Project_h
#include "tensors.h"

void L2ProjectGradAdd(int istart, int iend, int jstart, int jend,
                      const dTensorBC4& qin,
                      const dTensorBC4& auxin, dTensorBC4& Fout,
                      void (*FluxFunc)(const dTensor2& xpts, const dTensor2& Q,
                          const dTensor2& Aux, dTensor3& flux));

void L2ProjectGradAdd(const int istart, 
                      const int iend, 
                      const int jstart, 
                      const int jend,
                      const int QuadOrder, 
                      const int BasisOrder_qin,
                      const int BasisOrder_auxin,
                      const int BasisOrder_fout,
                      const dTensorBC4* qin,
                      const dTensorBC4* auxin, 
                      dTensorBC4* fout,
                      void Func(const dTensor2&,const dTensor2&,
                                const dTensor2&,dTensor3&));

void L2Project(int istart, int iend, int jstart, int jend,
               const dTensorBC4& q,
               const dTensorBC4& aux, dTensorBC4& Fout,
               void (*Func)(const dTensor2& xpts,
                            const dTensor2& qvals,
                            const dTensor2& auxvals,
                            dTensor2& source));

void L2Project(const int istart,
               const int iend,
               const int jstart,
               const int jend,
               const int QuadOrder,
               const int BasisOrder_qin,
               const int BasisOrder_auxin,
               const int BasisOrder_fout,    
               const dTensorBC4* qin,
               const dTensorBC4* auxin,
               dTensorBC4* fout,
               void (*Func)(const dTensor2&,const dTensor2&,
                            const dTensor2&,dTensor2&));

void L2Project_extra(int istart, int iend, int jstart, int jend,
               const dTensorBC4& q,
               const dTensorBC4& aux, dTensorBC4& Fout,
               void (*Func)(const dTensor2&,const dTensor2&,
                            const dTensor2&,dTensor2&, void* data),
               void* data);

void L2Project_extra(const int istart,
                     const int iend,
                     const int jstart,
                     const int jend,
                     const int QuadOrder,
                     const int BasisOrder_qin,
                     const int BasisOrder_auxin,
                     const int BasisOrder_fout,    
                     const dTensorBC4* qin,
                     const dTensorBC4* auxin,
                     dTensorBC4* fout,
                     void (*Func)(const dTensor2&,const dTensor2&,
                                  const dTensor2&,dTensor2&, void* data),
                     void* data);

#endif // L2Project_h
