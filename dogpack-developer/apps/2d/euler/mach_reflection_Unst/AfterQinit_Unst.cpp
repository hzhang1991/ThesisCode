#include "dogdefs.h"
#include "mesh.h"

// Function that is called after the intial conditions are set
void AfterQinit_Unst(const mesh& Mesh, dTensor3& aux, dTensor3& q)
{

    void ApplyPosLimiter_Unst(const mesh& Mesh, const dTensor3& aux, dTensor3& q);
    ApplyPosLimiter_Unst(Mesh,aux,q);

    const int NumElems = Mesh.get_NumElems();
    const int NumPhysElems = Mesh.get_NumPhysElems();
    const int NumEdges = Mesh.get_NumEdges();
    const int NumBndEdges = Mesh.get_NumBndEdges();
    const int meqn = q.getsize(2);
    const int kmax = q.getsize(3);
    for(int i =1;i<=NumElems;i++)
    {
       for(int m=1;m<=meqn;m++)
       {
           for(int k=2;k<=kmax;k++)
           {
              q.set(i,m,k,0.0);
           }
       }
    }
}
