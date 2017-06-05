#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, 
		    const dTensor2& auxvals, dTensor2& source)
{
  int i,m;
  int numpts=xpts.getsize(1);
  int meqn=source.getsize(2);
  double x,y;
  
  for (i=1; i<=numpts; i++)
    {
      x = xpts.get(i,1);
      y = xpts.get(i,2);
      
      for (m=1; m<=meqn; m++)
	{
	  source.set(i,m, 0.0e0 );
	}
    }
}
