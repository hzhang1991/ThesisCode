#include "dogdefs.h"
#include "DogSolverCart3.h"

// Function that is called if time step is rejected due to CFL violation
void AfterReject(double dt, 
		 dTensorBC5& aux, 
		 dTensorBC5& q, 
		 DogSolverCart3& solver)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int mz   = q.getsize(3);
  const int meqn = q.getsize(4);
  const int kmax = q.getsize(5);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(4);
}
