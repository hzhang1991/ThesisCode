#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "constants.h"
#include "tensors.h"
using namespace std;

void RelaxLimiter(const dTensor2& node, 
		  dTensorBC3& aux, 
		  dTensorBC3& q)
{
  const int melems = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   kmax = q.getsize(3);
  const int    mbc = q.getmbc();
  const int   maux = aux.getsize(2);
  const double dx  = dogParamsCart1.get_dx();
  void Relax(const dTensorBC1& eps, 
	     const dTensorBC3& aux,
	     dTensorBC3& q);

  // Indicator for artificial viscosity limiter
  int m = 1;
  double s0   = 0.05e0/pow(double(kmax),4);
  dTensorBC1 eps(melems,mbc);
  int mz = 0;

  printf("why am i here.\n");

#pragma omp parallel for
  for (int i=(1-mbc); i<=(melems+mbc); i++)
    {
      double Q2sum = 0.0;
      
      for (int k=1; k<=(kmax-1); k++)
	{  Q2sum = Q2sum + pow(q.get(i,m,k),2);  }
      
      double Q2kmax = pow(q.get(i,m,kmax),2);
      Q2sum = Q2sum + Q2kmax;

      double se;
      if (Q2kmax>1.0e-10)
	{  se =  (Q2kmax/Q2sum);  }
      else
	{  se = 0.0e0;  }

      if (se>=s0)
	{  
	  mz = mz + 1;  	  

	  // Set the artificial diffusion parameter
	  // in this element
	  eps.set(i, 2.0*dx/double(kmax) );
	}
      else
	{  eps.set(i, 0.0 );  }

      eps.set(i, 2.0*dx/double(kmax) );
      //eps.set(i, 0.0 );
    }

  // Apply relaxation
  Relax(eps,aux,q);
  
}


// Compute explicit time-stepping diffusion operator
void Relax(const dTensorBC1& eps, 
	   const dTensorBC3& aux,
	   dTensorBC3& q)
{
  const int    mx = q.getsize(1);
  const int  meqn = q.getsize(2);
  const int  kmax = q.getsize(3);
  const int   mbc = q.getmbc();
  const double dt = dogParams.get_dt();

#pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      double alpha2 = pow(aux.get(i,2,1),2);
      double kappa = eps.get(i)/alpha2;

      for (int k=1; k<=kmax; k++)
	{
	  double An = q.get(i,1,k);
	  double Bn = q.get(i,2,k);

	  if (kappa>1.0e-14)
	    {
	      q.set(i,2,k, (Bn-An)*exp(-dt/kappa) + An );
	    }
	  else
	    {
	      q.set(i,2,k, An );		    
	    }
	}
    }

}
