#include "tensors.h"
#include "assert.h"

// deprecated
void CopyQ(const dTensorBC4& qin,dTensorBC4& qout)
{
    qout.copyfrom(qin);
}
