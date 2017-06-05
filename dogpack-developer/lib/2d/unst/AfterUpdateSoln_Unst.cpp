#include "dogdefs.h"
#include "mesh.h"

// Optional call to modify updated solution
// This is done *before* limiters are applied.
void AfterUpdateSoln_Unst(const mesh& Mesh, 
        dTensor3& aux, 
        dTensor3& q,
        double dt,
        double beta)
{
    const int NumElems = q.getsize(1);
    const int meqn     = q.getsize(2);
    const int kmax     = q.getsize(3);
    const int maux     = aux.getsize(2);
}
