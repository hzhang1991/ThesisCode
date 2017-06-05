#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "DogParamsCart1.h"
#include "constants.h"
#include "tensors.h"
#include "diffusion.h"
using namespace std;

void ArtificialViscosity(const dTensor2& node, dTensorBC3& aux, 
			 dTensorBC3& q, dTensorBC3& Ldiff)
{
  const int melems = q.getsize(1);
  const int   meqn = q.getsize(2);
  const int   kmax = q.getsize(3);
  const int    mbc = q.getmbc();
  const int   maux = aux.getsize(2);
  const double dx  = node.get(2,1)-node.get(1,1);
  void CreateLdiff(const dTensorBC1& eps, 
		   const dTensorBC3& q, dTensorBC3& Ldiff);
  
  if (meqn>1)
    {
      cout << endl;
      cout << " ERROR in ArtificialViscosity.cpp: " << endl;
      cout << "     Currently only supported for scalar equations. " << endl;
      cout << "     meqn = " << meqn << endl;
      cout << endl;
      exit(1);
    }

  // Indicator for artificial viscosity limiter
  int m = 1;
  double s0   = 0.05e0/pow(double(kmax),4);
  double eps0;
  dTensorBC1 eps(melems,mbc);
  int mz = 0;

  //#pragma omp parallel for
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
    }

  if (mz>0)
    {  CreateLdiff(eps,q,Ldiff);  }

}


// Compute explicit time-stepping diffusion operator
void CreateLdiff(const dTensorBC1& eps, 
		 const dTensorBC3& q, dTensorBC3& Lstar)
{
  int i,k,ell;
  double x; 
  int    mx = q.getsize(1);
  int  meqn = q.getsize(2);
  int  kmax = q.getsize(3);
  int   mbc = q.getmbc();
  double dx = dogParamsCart1.get_dx();
  double Stmp;
  dTensorBC3  V(mx,meqn,kmax,mbc);
  dTensorBC2 FU(mx,1,mbc);
  dTensorBC2 FV(mx,1,mbc);
  
  // fluxes to evaluate gradient
  //#pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
      FU.set(i,1, -signv[0]*phiv[0]*q.get(i,1,1) );
      
      for (int k=2; k<=kmax; k++)
        {  FU.set(i,1, FU.get(i,1) - signv[k-1]*phiv[k-1]*q.get(i,1,k) );  }
    }
  
  // gradient
  //#pragma omp parallel for
  for (int i=(1-mbc); i<=(mx+mbc-1); i++)
    for (int k=1; k<=kmax; k++)
      {
        Stmp = Smat[k-1][0]*q.get(i,1,1);
        for (ell=2; ell<=kmax; ell++)
          {  Stmp = Stmp + Smat[k-1][ell-1]*q.get(i,1,ell);  }
        
        V.set(i,1,k, (phiv[k-1]*(FU.get(i+1,1) + 
                                 signv[k-1]*FU.get(i,1)) - Stmp)/dx );
      }
  
  // fluxes to evaluate solution
  //#pragma omp parallel for
  for (int i=(2-mbc); i<=(mx+mbc); i++)
    {
      FV.set(i,1, phiv[0]*V.get(i-1,1,1) );
      
      for (int k=2; k<=kmax; k++)
        {
          FV.set(i,1, FV.get(i,1) + phiv[k-1]*V.get(i-1,1,k) );   
        }
      
      FV.set(i,1, FV.get(i,1)*0.5*(eps.get(i)+eps.get(i-1)) );
      
    }
  
  // update solution
  //#pragma omp parallel for
  for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    for (int k=1; k<=kmax; k++)
      {
        Stmp = Smat[k-1][0]*V.get(i,1,1);
        for (ell=2; ell<=kmax; ell++)
          {  Stmp = Stmp + Smat[k-1][ell-1]*V.get(i,1,ell);  }
	Stmp = Stmp * eps.get(i);
        
        Lstar.set(i,1,k, Lstar.get(i,1,k) 
                  + (phiv[k-1]*(FV.get(i+1,1) + 
                                signv[k-1]*FV.get(i,1)) - Stmp)/dx );
      }
}
