#include "dogdefs.h"

// Optional call to modify updated solution
// This is done *before* limiters are applied.
void AfterUpdateSoln(const dTensorBC4& aux,
		     dTensorBC4& q,
		     double dt,
		     double beta)
{
}
