#include "dogdefs.h"

void CopyQ_Unst(const dTensor3& qin, dTensor3& qout)
{
    const int NumElems = qin.getsize(1);
    const int meqn     = qin.getsize(2);
    const int kmax     = qin.getsize(3);

#if( NDIMS == 2 )
#pragma omp parallel for
#endif
    for (int i=1; i<=NumElems; i++)
    for (int k=1; k<=kmax; k++)
    for (int m=1; m<=meqn; m++)
    {
        double tmp1 = qin.get(i,m,k);
        qout.set(i,m,k, tmp1);
    }

}
