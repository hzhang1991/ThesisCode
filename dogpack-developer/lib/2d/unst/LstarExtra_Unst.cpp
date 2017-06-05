#include "dogdefs.h"
#include "mesh.h"

// Optional call to function that adds extra piece into
// solution right-hand side (i.e., q_t = Lstar)
void LstarExtra_Unst(const mesh& Mesh, 
		     const dTensor3* q, 
		     const dTensor3* aux,
		     dTensor3* Lstar)
{
}
