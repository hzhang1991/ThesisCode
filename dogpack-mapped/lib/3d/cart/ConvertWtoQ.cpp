#include "dogdefs.h"
#include "dog_math.h"

// Convert characteristic variables to conserved variables
void ConvertWtoQ(int ixy, 
		 const dTensorBC5& aux, 
		 const dTensorBC5& qold,
                 const dTensorBC5& win, 
		 dTensorBC5& qout,
		 void (*ProjectRightEig)(int,
					 const dTensor1&,
					 const dTensor1&,
					 const dTensor2&,
					 dTensor2&))
{
  const int mx    = qold.getsize(1);
  const int my    = qold.getsize(2);
  const int mz    = qold.getsize(3);
  const int meqn  = qold.getsize(4);
  const int mbc   = qold.getmbc();
  const int maux  = aux.getsize(4);
  const int ksize = win.getsize(5);

  // ----------------------------------------------------------
  //
  // Key: storage of characteristic variables
  //       
  //      wx_cent(i,j,k,me,1)  = Lx * q(i,j,k,me,2)
  //      wx_cent(i,j,k,me,2)  = Lx * q(i,j,k,me,5)
  //      wx_cent(i,j,k,me,3)  = Lx * q(i,j,k,me,6)
  //      wx_cent(i,j,k,me,4)  = Lx * q(i,j,k,me,8)
  //
  //      wy_cent(i,j,k,me,1)  = Ly * q(i,j,k,me,3)
  //      wy_cent(i,j,k,me,2)  = Ly * q(i,j,k,me,5)
  //      wy_cent(i,j,k,me,3)  = Ly * q(i,j,k,me,7)
  //      wy_cent(i,j,k,me,4)  = Ly * q(i,j,k,me,9)
  //
  //      wz_cent(i,j,k,me,1)  = Lz * q(i,j,k,me,4)
  //      wz_cent(i,j,k,me,2)  = Lz * q(i,j,k,me,6)
  //      wz_cent(i,j,k,me,3)  = Lz * q(i,j,k,me,7)
  //      wz_cent(i,j,k,me,4)  = Lz * q(i,j,k,me,10)
  //
  // ---------------------------------------------------------- 
  
  //
  // Loop over all interior grid cells
  //
  switch( ixy )
    {
    case 1:      
      // ############################################
      //   X-DIRECTION
      // ############################################
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)
	for (int j=(3-mbc); j<=(my+mbc-2); j++)
	  for (int k=(3-mbc); k<=(mz+mbc-2); k++)
	    {
	      dTensor1 Aux_ave(maux),Q_ave(meqn);
	      dTensor2 Qxvals(meqn,ksize);
	      dTensor2 Wxvals(meqn,ksize);
	      
              // Store w and aux values in temporary arrays
              for (int m=1; m<=meqn; m++)
		{  Q_ave.set(m, qold.get(i,j,k,m,1) );  }
	                    
              for (int m=1; m<=maux; m++)
		{  Aux_ave.set(m, aux.get(i,j,k,m,1) );  }
	                    
              for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  Wxvals.set(m,ell, win.get(i,j,k,m,ell) );  }		
              
	      // X-Direction: Project onto right eigenvectors
              ProjectRightEig(1,Aux_ave,Q_ave,Wxvals,Qxvals);
              
              // Place results of Qvals into larger array qout
              for (int m=1; m<=meqn; m++)
		{		
                  qout.set(i,j,k,m,2, Qxvals.get(m,1) );
		  
                  if (ksize>1)
		    {
                      qout.set(i,j,k,m,5, Qxvals.get(m,2) );
                      qout.set(i,j,k,m,6, Qxvals.get(m,3) );
		      qout.set(i,j,k,m,8, Qxvals.get(m,4) );
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
	  for (int k=(3-mbc); k<=(mz+mbc-2); k++)
	    {
	      dTensor1 Aux_ave(maux),Q_ave(meqn);
	      dTensor2 Qyvals(meqn,ksize);
	      dTensor2 Wyvals(meqn,ksize);

              // Store w and aux values in temporary arrays
              for (int m=1; m<=meqn; m++)
		{  Q_ave.set(m, qold.get(i,j,k,m,1) );  }
              
              for (int m=1; m<=maux; m++)
		{  Aux_ave.set(m, aux.get(i,j,k,m,1) );  }
              
              for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  Wyvals.set(m,ell, win.get(i,j,k,m,ell) );  }		
              
              // Y-Direction: Project onto right eigenvectors
              ProjectRightEig(2,Aux_ave,Q_ave,Wyvals,Qyvals);

              // Place results of Qvals into larger array qout
              for (int m=1; m<=meqn; m++)
		{		
                  qout.set(i,j,k,m,3, Qyvals.get(m,1) );
		  
                  if (ksize>1)
		    {
		      // Here we assume that this method is called first
		      // with ixy=1, then with ixy=2, and finally with ixy=3.
		      // This is to guarantee that we take the minmod of
		      // the limiting in each direction for the mixed term.	
		      double tmp5 = minmod( Qyvals.get(m,2), qout.get(i,j,k,m,5) );

                      qout.set(i,j,k,m,5, tmp5            );
                      qout.set(i,j,k,m,7, Qyvals.get(m,3) );
		      qout.set(i,j,k,m,9, Qyvals.get(m,4) );
		    }
		}
	    }
      break;

    case 3:
      // ############################################
      //   Z-DIRECTION
      // ############################################
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)
	for (int j=(3-mbc); j<=(my+mbc-2); j++)
	  for (int k=(3-mbc); k<=(mz+mbc-2); k++)
	    {
	      dTensor1 Aux_ave(maux),Q_ave(meqn);
	      dTensor2 Qzvals(meqn,ksize);
	      dTensor2 Wzvals(meqn,ksize);

              // Store w and aux values in temporary arrays
              for (int m=1; m<=meqn; m++)
		{  Q_ave.set(m, qold.get(i,j,k,m,1) );  }
              
              for (int m=1; m<=maux; m++)
		{  Aux_ave.set(m, aux.get(i,j,k,m,1) );  }
              
              for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  Wzvals.set(m,ell, win.get(i,j,k,m,ell) );  }		
              
              // Y-Direction: Project onto right eigenvectors
              ProjectRightEig(3,Aux_ave,Q_ave,Wzvals,Qzvals);	      

              // Place results of Qvals into larger array qout
              for (int m=1; m<=meqn; m++)
		{		
                  qout.set(i,j,k,m,4, Qzvals.get(m,1) );
		  
                  if (ksize>1)
		    {
		      // Here we assume that this method is called first
		      // with ixy=1, then with ixy=2, and finally with ixy=3.
		      // This is to guarantee that we take the minmod of
		      // the limiting in each direction for the mixed term.	
		      double tmp6 = minmod( Qzvals.get(m,2), qout.get(i,j,k,m,6) );
		      double tmp7 = minmod( Qzvals.get(m,3), qout.get(i,j,k,m,7) );

                      qout.set(i,j,k,m,6,  tmp6            );
                      qout.set(i,j,k,m,7,  tmp7            );
		      qout.set(i,j,k,m,10, Qzvals.get(m,4) );
		    }
		}
	    }
      break;
    }
  
}
