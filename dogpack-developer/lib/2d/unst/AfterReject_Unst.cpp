#include "dogdefs.h"
#include "mesh.h"

// Function that is called if time step is rejected due to CFL violation
void AfterReject_Unst(const mesh& Mesh, const double dt, dTensor3& aux, dTensor3& q)
{
  const int NumElems = q.getsize(1);
  const int meqn     = q.getsize(2);
  const int kmax     = q.getsize(3);
  const int maux     = aux.getsize(2);
}
