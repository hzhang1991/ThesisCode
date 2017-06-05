#include "tensors.h"

// Convert conserved variables to characteristic variables
void ConvertQtoW(const dTensorBC3& aux, 
		 const dTensorBC3& qold, 
		 dTensorBC3& dw_right, 
		 dTensorBC3& dw_left,
		 dTensorBC3& w_cent,
		 void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{
  const int melems = qold.getsize(1);
  const int   meqn = qold.getsize(2);
  const int   kmax = qold.getsize(3);
  const int   maux = aux.getsize(2);
  const int mbc = qold.getmbc();
  dTensor1 Aux_ave(maux),Q_ave(meqn);
  dTensor2 Qin(meqn,kmax-1),Wout(meqn,kmax-1);
  
  //
  // Loop over all interior grid cells
  //
  for (int i=(3-mbc); i<=(melems+mbc-2); i++)
    {
      // Store q and aux values in temporary arrays
      for (int m=1; m<=meqn; m++)
	{  Q_ave.set(m, qold.get(i,m,1) );  }
      
      for (int m=1; m<=maux; m++)
	{  Aux_ave.set(m, aux.get(i,m,1) );  }
      
      // ---------------------------------------
      // Part I:  dw_right = L*( q(i+1) - q(i) )
      // ---------------------------------------
      for (int m=1; m<=meqn; m++)       
	for (int k=1; k<=(kmax-1); k++)
	  {  Qin.set(m,k, qold.get(i+1,m,k) - qold.get(i,m,k) );  }
	
      // Project Qin onto left eigenvectors in cell i
      ProjectLeftEig(Aux_ave,Q_ave,Qin,Wout);
      
      // Store results in dw_right
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=(kmax-1); k++)
	  {  dw_right.set(i,m,k, Wout.get(m,k) );  }
      // ---------------------------------------
      
      // ---------------------------------------
      // Part II:  dw_left = L*( q(i) - q(i-1) )
      // ---------------------------------------
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=(kmax-1); k++)
	  {  Qin.set(m,k, qold.get(i,m,k) - qold.get(i-1,m,k) );  }
  
      // Project Qin onto left eigenvectors in cell i
      ProjectLeftEig(Aux_ave,Q_ave,Qin,Wout);
      
      // Store results in dw_left
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=(kmax-1); k++)
	  {  dw_left.set(i,m,k, Wout.get(m,k) );  }    
      // ---------------------------------------
      
      // ---------------------------------------
      // Part III:  w_cent = L*q(i)
      // ---------------------------------------
      for (int m=1; m<=meqn; m++)
	for (int k=1; k<=(kmax-1); k++)
	  {  Qin.set(m,k, qold.get(i,m,k+1) );  }	
      
      // Project Qin onto left eigenvectors in cell i
      ProjectLeftEig(Aux_ave,Q_ave,Qin,Wout);
      
      // Store results in w_cent
      for (int m=1; m<=meqn; m++)	
	for (int k=1; k<=(kmax-1); k++)
	  {  w_cent.set(i,m,k, Wout.get(m,k) );  }
      // ---------------------------------------
    }
}
