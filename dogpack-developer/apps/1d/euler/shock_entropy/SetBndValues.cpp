#include "tensors.h"

// Zeroth order extrapolation boundary conditions
void SetBndValues(const dTensor2& node, 
		  dTensorBC3& aux, 
		  dTensorBC3& q)
{
  const int melems = q.getsize(1);
  const int meqn   = q.getsize(2);
  const int kmax   = q.getsize(3);
  const int maux   = aux.getsize(2);
  const int mbc    = q.getmbc();
    
  for (int ell=1; ell<=kmax; ell++)
    { 
      // ***********************************************
      // LEFT BOUNDARY
      // ***********************************************
      for (int i=0; i>=(1-mbc); i--)
        {        
	  // q values
	  for (int m=1; m<=meqn; m++)
	    {
	      double tmp = q.get(1,m,ell);
                    
	      q.set(i,m,ell, tmp );
            }
	  
	  // aux values
	  for (int m=1; m<=maux; m++)
            {
	      double tmp = aux.get(1,m,ell);
                    
	      aux.set(i,m,ell, tmp );
            }
        }
      // ***********************************************  

      
      // ***********************************************
      // RIGHT BOUNDARY
      // ***********************************************
      for (int i=(melems+1); i<=(melems+mbc); i++)
        {        
	  // q values
	  for (int m=1; m<=meqn; m++)
	    {
	      double tmp = q.get(melems,m,ell);
                    
	      q.set(i,m,ell, tmp );
            }
	  
	  // aux values
	  for (int m=1; m<=maux; m++)
            {
	      double tmp = aux.get(melems,m,ell);
              
	      aux.set(i,m,ell, tmp );
            }
        }
      // ***********************************************
      
    }
  
}
