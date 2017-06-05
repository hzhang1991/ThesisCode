#include "tensors.h"

// Convert characteristic variables to conserved variables
void ConvertWtoQ(const dTensorBC3& aux,
		 const dTensorBC3& qold, 
		 const dTensorBC3& win, 
		 dTensorBC3& qout,
		 void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{
  const int melems = qold.getsize(1);
  const int   meqn = qold.getsize(2);
  const int   kmax = qold.getsize(3);
  const int   maux = aux.getsize(2);
  const int mbc = qold.getmbc();
  dTensor1 Aux_ave(maux),Q_ave(meqn);
  dTensor2 Qvals(meqn,kmax-1),Wvals(meqn,kmax-1);
  
  //
  // Loop over all interior grid cells
  //
  for (int i=(3-mbc); i<=(melems+mbc-2); i++)
    {
      // Store w and aux values in temporary arrays
      for (int m=1; m<=meqn; m++)
	{  Q_ave.set(m, qold.get(i,m,1) );  }
      
      for (int m=1; m<=maux; m++)
	{  Aux_ave.set(m, aux.get(i,m,1) );  }
      
      for (int m=1; m<=meqn; m++)
	{  
	  for (int k=1; k<=(kmax-1); k++)
	    {  Wvals.set(m,k, win.get(i,m,k) );  }
	}
      
      // Project Qvals onto left eigenvectors in cell i,
      //    the result is stored in Wvals
      ProjectRightEig(Aux_ave,Q_ave,Wvals,Qvals);
      
      // Place results of Qvals into larger array qout
      for (int k=2; k<=kmax; k++)	
	for (int m=1; m<=meqn; m++)
	  {
	    qout.set(i,m,k, Qvals.get(m,k-1) );
	  }
    }

}
