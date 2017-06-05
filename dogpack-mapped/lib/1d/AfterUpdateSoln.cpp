#include "tensors.h"

// Optional call to modify updated solution
// This is done *before* limiters are applied.
void AfterUpdateSoln(const dTensor2& node,
		     const dTensorBC3& aux,
		     dTensorBC3& q,
		     double dt,
		     double beta)
{
}
