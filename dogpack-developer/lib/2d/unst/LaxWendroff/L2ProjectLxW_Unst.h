#ifndef _L2PROJECT_LXW_UNST_H_
#define _L2PROJECT_LXW_UNST_H_

// Quadrature rules (see: lib/2d/unst/Legendre2d_Unst)
void setQuadPoints_Unst(int QuadOrder, dTensor1& wgts, dTensor2& spts);
void SetLegendreAtPoints_Unst(const dTensor2& spts, dTensor2& phi);
void SetLegendreGrad_Unst( const dTensor2& spts, dTensor2& phi_xi, dTensor2& phi_eta );
void LegendreDiff2_Unst(const dTensor2& spts, 
    dTensor2* phi_xi2, dTensor2* phi_xieta, dTensor2* phi_eta2 );

#endif
