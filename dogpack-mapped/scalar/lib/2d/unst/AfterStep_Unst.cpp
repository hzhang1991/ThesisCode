#include "dogdefs.h"
#include "mesh.h"

// Function that is called after each time step
void AfterStep_Unst(const double dt, const mesh& Mesh, dTensor3& aux, dTensor3& q)
{
    const int NumElems = q.getsize(1);
    const int meqn     = q.getsize(2);
    const int kmax     = q.getsize(3);
    const int maux     = aux.getsize(2);
}
