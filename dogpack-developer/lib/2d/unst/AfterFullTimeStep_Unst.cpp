#include "dogdefs.h"
#include "mesh.h"

// Function that is called after a full time step (i.e., after all stages are complete)
void AfterFullTimeStep_Unst(const double dt, const mesh& Mesh,
			    const dTensor3& auxold, const dTensor3& qold,
			    const dTensor3& Lold, dTensor3& aux, dTensor3& q)
{
}
