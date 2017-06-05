#include "tensors.h"

// Function that is called before a full time step
void BeforeFullTimeStep(double dt, const dTensor2& node, const dTensor1& prim_vol,
		       dTensorBC3& auxold, dTensorBC3& aux, 
		       dTensorBC3& qold, dTensorBC3& q)
{
  const int melems = q.getsize(1);
  const int meqn   = q.getsize(2);
  const int kmax   = q.getsize(3);
  const int mbc    = q.getmbc();
  const int maux   = aux.getsize(4);
}
