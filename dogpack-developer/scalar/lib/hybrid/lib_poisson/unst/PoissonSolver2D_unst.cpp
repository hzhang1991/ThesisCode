#include "dogdefs.h"
#include "edge_data_Unst.h"
#include "mesh.h"
#include "SparseCholesky.h"

void PoissonSolver2D_unst(const int space_order,
        const mesh& Mesh,			  
        const SparseCholesky& R,
        const dTensor1& rhs,
        dTensor1& phi,
        dTensor2& E1,
        dTensor2& E2)
{
    // Solve linear system using cholesky factorization + forward/backward substitution
    const int N = rhs.getsize();
    dTensor1 rhs_tmp(N);
    R.ForwardSubs(rhs,rhs_tmp);
    R.BackwardSubs(rhs_tmp,phi);

    // Compute electric field as gradient of the electric potential
    void ComputeEfield(const int space_order,
            const mesh& Mesh,
            const dTensor1& phi,
            dTensor2& E1,
            dTensor2& E2);
    ComputeEfield(space_order,Mesh,phi,E1,E2);
}
