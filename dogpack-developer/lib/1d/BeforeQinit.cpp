#include "tensors.h"
#include "dogdefs.h"
#include <cmath>

// Function that is called before initial conditions are set
void BeforeQinit(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q)
{
  const int melems = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   kmax = q.getsize(3);
  const int   mbc  = q.getmbc();
  const int   maux = aux.getsize(2);

 
}
