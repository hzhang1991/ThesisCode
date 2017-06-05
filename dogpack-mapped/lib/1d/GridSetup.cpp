#include "tensors.h"

void GridSetup(int method1, double xlow, double dx, dTensor2& node,
	       dTensor1& prim_vol)
{
  const int mnodes = node.getsize(1);
  const int melems = prim_vol.getsize();
  const double zmeth = double(method1);
  
  // Set variable "node"
#pragma omp parallel for
  for (int i=1; i<=mnodes; i++)
    {  node.set(i,1, xlow+(double(i)-1.0e0)*dx );  }
  
  // Set grid cell volumes
#pragma omp parallel for
  for (int j=1; j<=melems; j++)
    {  prim_vol.set(j, dx);  }
    
}
