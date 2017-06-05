#include "tensors.h"

// *REQUIRED*
//
// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
// For the 1D balance law: 
//
//            q_t + f_x = psi,
//
// The diagonalized matrix is f'(q) = R \Lambda R^{-1}.
//
// This routine computes W := R^{-1} * Q.
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
// See also: ProjectRightEig, SetWaveSpd.
void ProjectLeftEig(const dTensor1& Aux_ave, const dTensor1& Q_ave, const dTensor2& Qvals,
		    dTensor2& Wvals)
{

    const int meqn = Qvals.getsize(1);
    const int mpts = Qvals.getsize(2);

}
