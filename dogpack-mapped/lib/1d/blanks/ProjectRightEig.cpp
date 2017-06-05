#include "tensors.h"

// *REQUIRED*
//
// This is a user-supplied routine that projects
// Qvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
// For the 1D balance law: 
//
//            q_t + f_x = psi
//
// The diagonalized matrix is f'(q) = R \Lambda R^{-1}.
//
// This routine computes Q := R * W.
//
// Input:
//
// Aux_ave( 1:maux )
//   Q_ave( 1:meqn )
//   Qvals( 1:meqn, 1:numpts )
//
// Output:
//
//   Wvals( 1:meqn, 1:numpts )
//
// See also: ProjectRightEig.
void ProjectRightEig(const dTensor1& Aux_ave, const dTensor1& Q_ave, const dTensor2& Wvals,
		     dTensor2& Qvals)
{

    const int meqn = Qvals.getsize(1);
    const int mpts = Qvals.getsize(2);

}
