#include "dogdefs.h"
#include "dog_math.h"

// Convert characteristic variables to conserved variables
void ConvertWtoQ(int ixy, 
		 const dTensorBC4* aux, 
		 const dTensorBC4& qold,
                 const dTensorBC4& win, 
		 dTensorBC4& qout,
		 void (*ProjectRightEig)(int,
					 const dTensor1&,
					 const dTensor1&,
					 const dTensor2&,
					 dTensor2&))
{
  const int mx    = qold.getsize(1);
  const int my    = qold.getsize(2);
  const int meqn  = qold.getsize(3);
  const int mbc   = qold.getmbc();
  const int maux  = aux->getsize(3);
  const int ksize = win.getsize(4);
  
  // ----------------------------------------------------------
  // Key: storage of characteristic variables
  //       
  //      wx_cent(i,j,me,1)  = Lx * q(i,j,me,2)
  //      wx_cent(i,j,me,2)  = Lx * q(i,j,me,4)
  //      wx_cent(i,j,me,3)  = Lx * q(i,j,me,5)
  //
  //      wy_cent(i,j,me,1)  = Ly * q(i,j,me,3)
  //      wy_cent(i,j,me,2)  = Ly * q(i,j,me,4)
  //      wy_cent(i,j,me,3)  = Ly * q(i,j,me,6)
  //
  // ----------------------------------------------------------
  
  //
  // Loop over all interior grid cells
  //
  switch ( ixy )
    {
    case 1:
      // ############################################
      //   X-DIRECTION
      // ############################################
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)
	for (int j=(3-mbc); j<=(my+mbc-2); j++)
	  {
	    dTensor1 Aux_ave(maux),Q_ave(meqn);
	    dTensor2 Qxvals(meqn,ksize);
	    dTensor2 Wxvals(meqn,ksize);
	    
	    // Store w and aux values in temporary arrays
	    for (int m=1; m<=meqn; m++)
	      {  Q_ave.set(m, qold.get(i,j,m,1) );  }
	    
	    for (int m=1; m<=maux; m++)
	      {  Aux_ave.set(m, aux->get(i,j,m,1) );  }
	    
	    for (int m=1; m<=meqn; m++)
	      for (int k=1; k<=ksize; k++)
		{  Wxvals.set(m,k, win.get(i,j,m,k) );  }
	                
	    // X-Direction: Project onto right eigenvectors
	    ProjectRightEig(1,Aux_ave,Q_ave,Wxvals,Qxvals);
              
	    // Place results of Qvals into larger array qout
	    for (int m=1; m<=meqn; m++)
	      {		
		qout.set(i,j,m,2, Qxvals.get(m,1) );
		
		if (ksize>1)
		  {                      
		    qout.set(i,j,m,4, Qxvals.get(m,2) );
		    qout.set(i,j,m,5, Qxvals.get(m,3) );
		  }
	      }
	  }    
      break;

    case 2:
      // ############################################
      //   Y-DIRECTION
      // ############################################
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)
	for (int j=(3-mbc); j<=(my+mbc-2); j++)
	  {
	    dTensor1 Aux_ave(maux),Q_ave(meqn);
	    dTensor2 Qyvals(meqn,ksize);
	    dTensor2 Wyvals(meqn,ksize);
	    
	    // Store w and aux values in temporary arrays
	    for (int m=1; m<=meqn; m++)
	      {  Q_ave.set(m, qold.get(i,j,m,1) );  }
	    
	    for (int m=1; m<=maux; m++)
	      {  Aux_ave.set(m, aux->get(i,j,m,1) );  }
	    
	    for (int m=1; m<=meqn; m++)
	      for (int k=1; k<=ksize; k++)
		{  Wyvals.set(m,k, win.get(i,j,m,k) );  }	      
	    
	    // Y-Direction: Project onto right eigenvectors
	    ProjectRightEig(2,Aux_ave,Q_ave,Wyvals,Qyvals);
            
	    // Place results of Qvals into larger array qout
	    for (int m=1; m<=meqn; m++)
	      {		
		qout.set(i,j,m,3, Qyvals.get(m,1) );
		
		if (ksize>1)
		  {
		    // Here we assume that this method is called first
		    // with ixy=1 and then with ixy=2.
		    // This is to guarantee that we take the minmod of
		    // the limiting in each direction for the mixed term.		      
		    double currval = qout.get(i,j,m,4);
		    double  newval = Qyvals.get(m,2);
		    newval = minmod(currval,newval);
		    
		    qout.set(i,j,m,4, newval );
		    qout.set(i,j,m,6, Qyvals.get(m,3) );
		  }
	      }
	  }
      break;
    }

}

// Convert characteristic variables to conserved variables
void ConvertWtoQ(int ixy, const dTensorBC4& aux, const dTensorBC4& qold,
                 const dTensorBC4& win, dTensorBC4& qout,
		 void (*ProjectRightEig)(int,const dTensor1&,
					 const dTensor1&,const dTensor2&,
					 dTensor2&))
{
  void ConvertWtoQ(int ixy, const dTensorBC4* aux, const dTensorBC4& qold,
                   const dTensorBC4& win, dTensorBC4& qout,
                   void (*ProjectRightEig)(int,const dTensor1&,
                                           const dTensor1&,const dTensor2&,
                                           dTensor2&));
  
  ConvertWtoQ(ixy,&aux,qold,win,qout,ProjectRightEig);
}
