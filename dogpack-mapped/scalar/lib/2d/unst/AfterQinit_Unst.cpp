#include "dogdefs.h"
#include "DogParams.h"
#include "mesh.h"

// The 'do nothing' replacable function
void AfterQinit_Unst(const mesh& Mesh, dTensorBC5& qnew)
{

    // initialize the vlasov parameters:

    const int mx       = qnew.getsize(1);
    const int my       = qnew.getsize(2);
    const int NumElems = qnew.getsize(3);
    const int meqn     = qnew.getsize(4);
    const int kmax     = qnew.getsize(5);
    const int mbc      = qnew.getmbc();

}
