#include "dogdefs.h"

// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig(int ixy, const dTensor1& Aux_ave, const dTensor1& Q_ave, 
		    const dTensor2& Qvals, dTensor2& Wvals)
{    
  const int meqn = Qvals.getsize(1);
  const int kmax = Qvals.getsize(2);
    
  // Project onto left eigenvectors
  for (int k=1; k<=kmax; k++)
    {
      Wvals.set(1,k, Qvals.get(1,k) );
    }
}
