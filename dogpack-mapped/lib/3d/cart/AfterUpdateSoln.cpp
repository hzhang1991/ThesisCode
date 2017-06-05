#include "dogdefs.h"

// Optional call to modify updated solution
// This is done *before* limiters are applied.
void AfterUpdateSoln(const dTensorBC5& aux,
		     dTensorBC5& q,
		     double dt,
		     double beta)
{
}
