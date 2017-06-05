#include "tensors.h"

// Euler time stepping for the spectral deferred correction method
void StepSDC(double dt, const dTensorBC4& aux, const dTensorBC4& qold, 
	     const dTensorBC4& Lrhs, dTensorBC4& qnew)
{
    const int mx     = qnew.getsize(1);
    const int my     = qnew.getsize(2);
    const int meqn   = qnew.getsize(3);
    const int kmax   = qnew.getsize(4);
    const int mbc    = qnew.getmbc();
    const int maux   = aux.getsize(3);

    #pragma omp parallel for
    // Euler time step loop
    for (int k=1; k<=kmax; k++)
      for (int m=1; m<=meqn; m++)
	for (int i=(1-mbc); i<=(mx+mbc); i++)
	  for (int j=(1-mbc); j<=(my+mbc); j++)	  
	    {
	      double tmp = qold.get(i,j,m,k) + dt*Lrhs.get(i,j,m,k);
	      qnew.set(i,j,m,k, tmp );
            }
}
