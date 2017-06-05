#include "dogdefs.h"
#include "FullMatrix.h"

void ConjugateGradient(const int MaxIters,
		       const double TOL,
		       const FullMatrix& A,
		       const dTensor1& rhs,
		       dTensor1& phi)
{
  double DotProd(const dTensor1& avec,
		 const dTensor1& bvec);
  void MatMult(const FullMatrix& A,
	       const dTensor1& invec,
	       dTensor1& outvec);
  const int size = A.get_NumRows();
  int NumIters  = 0;
  int mflag = 0;
  
  dTensor1 r(size);
  dTensor1 d(size);
  dTensor1 Ad(size);

  for (int i=1; i<=size; i++)
    {  phi.set(i, 0.0 );  }

  MatMult(A,phi,r);

  for (int i=1; i<=size; i++)
    {  r.set(i, rhs.get(i)-r.get(i) );  }

  for (int i=1; i<=size; i++)
    {  d.set(i, r.get(i) );  }
  
  double rhs_norm = sqrt(DotProd(rhs,rhs));
  if (fabs(rhs_norm)<=1.0e-12)
    { rhs_norm = 1.0; }
  //printf(" rhs_norm = %e\n",rhs_norm);

  double rdotr = DotProd(r,r);
  //printf(" sqrt(rdotr) = %e\n",sqrt(rdotr));

  double rel_res_norm = sqrt(rdotr)/rhs_norm;
  if(rel_res_norm < TOL)
    {  mflag = 1; }

  NumIters = 1;
  while(mflag==0)
    {
      MatMult(A,d,Ad);
      double alpha = rdotr/DotProd(d,Ad);
      
      for (int i=1; i<=size; i++)
	{  phi.set(i, phi.get(i)+alpha*d.get(i) );  }

      for (int i=1; i<=size; i++)
	{  r.set(i, r.get(i)-alpha*Ad.get(i) );  }
      
      double rdotr_old = rdotr;
      rdotr = DotProd(r,r);

      rel_res_norm = sqrt(rdotr)/rhs_norm;
      if(rel_res_norm < TOL)
	{  mflag = 1;  }
	  
      //printf("   NumIters = %i,   rel_res_norm = %e\n",NumIters,rel_res_norm);
	  
      if (NumIters==MaxIters)
	{  mflag = 1;  }
	  
      if (mflag==0)
	{
	  double beta = rdotr/rdotr_old;

	  for (int i=1; i<=size; i++)
	    {  d.set(i, r.get(i)+beta*d.get(i) );  }

	  NumIters = NumIters+1;
	}
    }

  printf("  |----------------------------\n");
  printf("  | Conjugate Gradient Results:\n");
  printf("  |----------------------------\n");
  printf("  |  MaxIters = %i\n",MaxIters);
  printf("  |       TOL = %e\n",TOL);
  printf("  |  NumIters = %i\n",NumIters);
  printf("  |  residual = %e\n",rel_res_norm);
  printf("  |----------------------------\n");
  printf("\n");
}


double DotProd(const dTensor1& avec,
	       const dTensor1& bvec)
{
  const int length = avec.getsize();
  double prod = 0.0;

  for (int i=1; i<=length; i++)
    {
      prod = prod + avec.get(i)*bvec.get(i);
    }

  return prod;
}

void MatMult(const FullMatrix& A,
	     const dTensor1& invec,
	     dTensor1& outvec)
{


  for (int i=1; i<=A.get_NumRows(); i++)
    {
      int jmax = A.get_NZrow(i);
      double tmp = 0.0;

      for (int j=1; j<=jmax; j++)
	{
	  int ind    = A.get_Index(i,j);
	  double val = A.get_Value(i,j);

	  tmp = tmp + val*invec.get(ind);
	}
      outvec.set(i, tmp );
    }

}
