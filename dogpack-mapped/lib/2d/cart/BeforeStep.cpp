#include "dogdefs.h"
#include "DogSolverCart2.h"

// Function that is called before each time stage
void BeforeStep(double dt, dTensorBC4& aux, dTensorBC4& q, DogSolverCart2& solver)
{
  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int kmax = q.getsize(4);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(3);
}
