#include "dogdefs.h"
#include "DogSolverCart2.h"

// Function that is called from ConstructL
// in order to enforce boundary conditions
// through fluxes
void SetBndFluxes(const dTensorBC4& q, 
		  const dTensorBC4& aux,
		  dTensorBC4& Fm,
		  dTensorBC4& Fp,
		  dTensorBC4& Gm,
		  dTensorBC4& Gp)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(3);
}
