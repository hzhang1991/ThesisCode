#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc(const dTensor2& xpts, 
		    const dTensor2& qvals,
		    const dTensor2& auxvals,
		    dTensor2& source)
{
  const int numpts=xpts.getsize(1);
  const int meqn=source.getsize(2);
  
  for (int i=1; i<=numpts; i++)
    {
      const double x = xpts.get(i,1);
      const double y = xpts.get(i,2);

      source.set(i,1,  0.0  );
      source.set(i,2,  0.0 );
      source.set(i,3,  0.0 );
    }
}
