#include "dogdefs.h"
#include "mesh.h"

// Function that is called after the intial conditions are set
void AfterQinit_Unst(const mesh& Mesh, dTensor3& aux, dTensor3& q)
{

    void ApplyPosLimiter_Unst(const mesh& Mesh, const dTensor3& aux, dTensor3& q);
    ApplyPosLimiter_Unst(Mesh,aux,q);
}
