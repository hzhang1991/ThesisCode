#include "Limiters.h"
#include <cmath>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "Legendre2d.h"

void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
		  void (*ProjectRightEig)(int,const dTensor1&,
					  const dTensor1&,const dTensor2&,
					  dTensor2&),
		  void (*ProjectLeftEig)(int,const dTensor1&,
					 const dTensor1&,const dTensor2&,
					 dTensor2&))
{
  ApplyLimiterKrivodonova(aux, q, ProjectRightEig, ProjectLeftEig);

  if (dogParams.get_use_divfree()==1)
    project_onto_locally_divergence_free_subspace(q);
}

#if 0
  const int   mx = q.getsize(1);
  const int   my = q.getsize(2);
  const int mbc  = q.getmbc();
  
  // Divergence-free fix-up (moment-limiting does not preserve div(B)=0, this
  //                         fix-up will restore div(B)=0 preservation)
  if (dogParams.get_use_divfree())
    {
      double minmod(double,double);
      const double dx = dogParamsCart2.get_dx();
      const double dy = dogParamsCart2.get_dy();
      const double length = sqrt(dx*dx + dy*dy);
      const double lx5    = sqrt(5.0*dx*dx + dy*dy);
      const double ly5    = sqrt(dx*dx + 5.0*dy*dy);
      const int how_many = dogParams.get_how_many_vectors_divfree();

      #pragma omp parallel for
      for (int k=1; k<=how_many; k++)
	{
	  int mcomp = dogParams.get_which_compnt_divfree()[k];

	  switch ( dogParams.get_space_order() )
	    {
	    case 2:  // 2nd order in space
	      
              #pragma omp parallel for
	      for (int i=(3-mbc); i<=(mx+mbc-2); i++)	
		for (int j=(3-mbc); j<=(my+mbc-2); j++)
		  {
		    const double B3 = minmod(q.get(i,j,mcomp,2)/dx,
                                            -q.get(i,j,mcomp+1,3)/dy);
		    
		    q.set(i,j,mcomp,  2,  dx*B3 );
		    q.set(i,j,mcomp+1,3, -dy*B3 );
		  }
	      
	      break;
	      
	    case 3:  // 3rd order in space
	      
              #pragma omp parallel for
	      for (int i=(3-mbc); i<=(mx+mbc-2); i++)	
		for (int j=(3-mbc); j<=(my+mbc-2); j++)
		  {
		    double B3 = minmod(q.get(i,j,mcomp,2)/dx,       -q.get(i,j,mcomp+1,3)/dy);
		    double B6 = minmod(q.get(i,j,mcomp,5)/dx,       -q.get(i,j,mcomp+1,4)/(sq5*dy));
		    double B7 = minmod(q.get(i,j,mcomp,4)/(sq5*dx), -q.get(i,j,mcomp+1,6)/dy);

		    q.set(i,j,mcomp,  2,      dx*B3 );
		    q.set(i,j,mcomp+1,3,     -dy*B3 );
		    
		    q.set(i,j,mcomp,  5,      dx*B6 );
		    q.set(i,j,mcomp+1,4, -sq5*dy*B6 );
		    
		    q.set(i,j,mcomp,  4,  sq5*dx*B7 );
		    q.set(i,j,mcomp+1,6,     -dy*B7 );		
		  }

	      break;
	    }
	}
    }
#endif

