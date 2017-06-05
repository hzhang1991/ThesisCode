#include "tensors.h"
#include "DogParams.h"

// Convert conserved variables to characteristic variables
void ConvertQtoW(int ixy, 
		 const dTensorBC5& aux, 
		 const dTensorBC5& qold, 
                 dTensorBC5& dwp, 
		 dTensorBC5& dwm, 
		 dTensorBC5& w_cent,
		 void (*ProjectLeftEig)(int,
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
  const int ksize = w_cent.getsize(5);  
  const int space_order = dogParams.get_space_order();

  // ----------------------------------------------------------
  //
  // Key: storage of characteristic variables
  //       
  //      wx_cent(i,j,k,me,1)  = Lx * q(i,j,k,me,2)
  //      wx_cent(i,j,k,me,2)  = Lx * q(i,j,k,me,5)
  //      wx_cent(i,j,k,me,3)  = Lx * q(i,j,k,me,6)
  //      wx_cent(i,j,k,me,4)  = Lx * q(i,j,k,me,8)
  //
  //      dw_right(i,j,k,me,1) = Lx * (q(i+1,j,k,me,1)-q(i,j,k,me,1))
  //      dw_right(i,j,k,me,2) = Lx * (q(i,j+1,k,me,2)-q(i,j,k,me,2))
  //      dw_right(i,j,k,me,3) = Lx * (q(i,j,k+1,me,2)-q(i,j,k,me,2))
  //      dw_right(i,j,k,me,4) = Lx * (q(i+1,j,k,me,2)-q(i,j,k,me,2))
  //
  //      dw_left(i,j,k,me,1)  = Lx * (q(i,j,k,me,1)-q(i-1,j,k,me,1))
  //      dw_left(i,j,k,me,2)  = Lx * (q(i,j,k,me,2)-q(i,j-1,k,me,2))
  //      dw_left(i,j,k,me,3)  = Lx * (q(i,j,k,me,2)-q(i,j,k-1,me,2))
  //      dw_left(i,j,k,me,4)  = Lx * (q(i,j,k,me,2)-q(i-1,j,k,me,2))
  //
  //      wy_cent(i,j,k,me,1)  = Ly * q(i,j,k,me,3)
  //      wy_cent(i,j,k,me,2)  = Ly * q(i,j,k,me,5)
  //      wy_cent(i,j,k,me,3)  = Ly * q(i,j,k,me,7)
  //      wy_cent(i,j,k,me,4)  = Ly * q(i,j,k,me,9)
  //
  //      dw_back(i,j,k,me,1)  = Ly * (q(i,j+1,k,me,1)-q(i,j,k,me,1))
  //      dw_back(i,j,k,me,2)  = Ly * (q(i+1,j,k,me,3)-q(i,j,k,me,3))
  //      dw_back(i,j,k,me,3)  = Ly * (q(i,j,k+1,me,3)-q(i,j,k,me,3))
  //      dw_back(i,j,k,me,4)  = Ly * (q(i,j+1,k,me,3)-q(i,j,k,me,3))
  //
  //      dw_front(i,j,k,me,1) = Ly * (q(i,j,k,me,1)-q(i,j-1,k,me,1))
  //      dw_front(i,j,k,me,2) = Ly * (q(i,j,k,me,3)-q(i-1,j,k,me,3))
  //      dw_front(i,j,k,me,3) = Ly * (q(i,j,k,me,3)-q(i,j,k-1,me,3))
  //      dw_front(i,j,k,me,4) = Ly * (q(i,j,k,me,3)-q(i,j-1,k,me,3))
  //
  //      wz_cent(i,j,k,me,1)  = Lz * q(i,j,k,me,4)
  //      wz_cent(i,j,k,me,2)  = Lz * q(i,j,k,me,6)
  //      wz_cent(i,j,k,me,3)  = Lz * q(i,j,k,me,7)
  //      wz_cent(i,j,k,me,4)  = Lz * q(i,j,k,me,10)
  //
  //      dw_up(i,j,k,me,1)    = Lz * (q(i,j,k+1,me,1)-q(i,j,k,me,1))
  //      dw_up(i,j,k,me,2)    = Lz * (q(i+1,j,k,me,4)-q(i,j,k,me,4))
  //      dw_up(i,j,k,me,3)    = Lz * (q(i,j+1,k,me,4)-q(i,j,k,me,4))
  //      dw_up(i,j,k,me,4)    = Lz * (q(i,j,k+1,me,4)-q(i,j,k,me,4))
  //
  //      dw_down(i,j,k,me,1)  = Lz * (q(i,j,k,me,1)-q(i,j,k-1,me,1))
  //      dw_down(i,j,k,me,2)  = Lz * (q(i,j,k,me,4)-q(i-1,j,k,me,4))
  //      dw_down(i,j,k,me,3)  = Lz * (q(i,j,k,me,4)-q(i,j-1,k,me,4))
  //      dw_down(i,j,k,me,4)  = Lz * (q(i,j,k,me,4)-q(i,j,k-1,me,4))
  //
  // ---------------------------------------------------------- 
   
  //
  // Loop over all interior grid cells
  //
  switch( ixy )
    {
    case 1:
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)
	for (int j=(3-mbc); j<=(my+mbc-2); j++)
	  for (int k=(3-mbc); k<=(mz+mbc-2); k++)
	    {
	      iTensor1 nn(4);
	      dTensor1 Aux_ave(maux),Q_ave(meqn);
	      dTensor2 Qin(meqn,ksize),Wout(meqn,ksize);
	      	      
	      // index for x-direction
	      nn.set(1, 2 );
	      nn.set(2, 5 );
	      nn.set(3, 6 );
	      nn.set(4, 8 );
	      
	      // Store q and aux values in temporary arrays
	      for (int m=1; m<=meqn; m++)
		{  Q_ave.set(m, qold.get(i,j,k,m,1) );   }
	      
	      for (int m=1; m<=maux; m++)
		{  Aux_ave.set(m, aux.get(i,j,k,m,1) );  }
	      
	      // ------------------------------------------------
	      // Part I:  dw_right = Lx*( q(i+1,j,k) - q(i,j,k) )
	      // ------------------------------------------------
	      switch( space_order )
		{
		case 3:
		  for (int m=1; m<=meqn; m++)
		    {  
		      Qin.set(m,1, qold.get(i+1,j,k,m,1) - qold.get(i,j,k,m,1) );
		      Qin.set(m,2, qold.get(i,j+1,k,m,2) - qold.get(i,j,k,m,2) );
		      Qin.set(m,3, qold.get(i,j,k+1,m,2) - qold.get(i,j,k,m,2) );
		      Qin.set(m,4, qold.get(i+1,j,k,m,2) - qold.get(i,j,k,m,2) );
		    }
		  break;
		  
		case 2:
		  for (int m=1; m<=meqn; m++)
		    {  
		      Qin.set(m,1, qold.get(i+1,j,k,m,1) - qold.get(i,j,k,m,1) );
		    }
		  break;
		}
	      
	      // Project Qin onto left eigenvectors in cell (i,j,k)
	      ProjectLeftEig(1,Aux_ave,Q_ave,Qin,Wout);
	      
	      // Store results in dw_right
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  dwp.set(i,j,k,m,ell, Wout.get(m,ell) );  }		  
	      // --------------------------------------------
	      
	      // ------------------------------------------------
	      // Part II:  dw_left = Lx*( q(i,j,k) - q(i-1,j,k) )
	      // ------------------------------------------------
	      switch( space_order )		  
		{
		case 3:
		  for (int m=1; m<=meqn; m++)
		    {
		      Qin.set(m,1, qold.get(i,j,k,m,1) - qold.get(i-1,j,k,m,1) );
		      Qin.set(m,2, qold.get(i,j,k,m,2) - qold.get(i,j-1,k,m,2) );
		      Qin.set(m,3, qold.get(i,j,k,m,2) - qold.get(i,j,k-1,m,2) );
		      Qin.set(m,4, qold.get(i,j,k,m,2) - qold.get(i-1,j,k,m,2) );
		    }
		  break;
		  
		case 2:
		  for (int m=1; m<=meqn; m++)
		    {  
		      Qin.set(m,1, qold.get(i,j,k,m,1) - qold.get(i-1,j,k,m,1) );
		    }
		  break;
		}
	      
	      // Project Qin onto left eigenvectors in cell (i,j,k)
	      ProjectLeftEig(1,Aux_ave,Q_ave,Qin,Wout);
		
	      // Store results in dw_left
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  dwm.set(i,j,k,m,ell, Wout.get(m,ell) );  }
	      // --------------------------------------------
		
	      // --------------------------------------------
	      // Part III:  wx_cent = Lx*q(i,j,k)
	      // --------------------------------------------
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  Qin.set(m,ell, qold.get(i,j,k,m,nn.get(ell)) );  }	      
	      
	      // Project Qin onto left eigenvectors in cell (i,j,k)
	      ProjectLeftEig(1,Aux_ave,Q_ave,Qin,Wout);
	
	      // Store results in wx_cent
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  w_cent.set(i,j,k,m,ell, Wout.get(m,ell) );  }	      
	      // --------------------------------------------
	    }      
      break;
      
    case 2:
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)
	for (int j=(3-mbc); j<=(my+mbc-2); j++)
	  for (int k=(3-mbc); k<=(mz+mbc-2); k++)
	    {
	      iTensor1 nn(4);
	      dTensor1 Aux_ave(maux),Q_ave(meqn);
	      dTensor2 Qin(meqn,ksize),Wout(meqn,ksize);
	      	      
	      // index for x-direction
	      nn.set(1, 3 );
	      nn.set(2, 5 );
	      nn.set(3, 7 );
	      nn.set(4, 9 );
	      
	      // Store q and aux values in temporary arrays
	      for (int m=1; m<=meqn; m++)
		{  Q_ave.set(m, qold.get(i,j,k,m,1) );   }
	      
	      for (int m=1; m<=maux; m++)
		{  Aux_ave.set(m, aux.get(i,j,k,m,1) );  }
	      
	      // ------------------------------------------------
	      // Part I:  dw_right = Ly*( q(i,j+1,k) - q(i,j,k) )
	      // ------------------------------------------------
	      switch( space_order )
		{
		case 3:
		  for (int m=1; m<=meqn; m++)
		    {  
		      Qin.set(m,1, qold.get(i,j+1,k,m,1) - qold.get(i,j,k,m,1) );
		      Qin.set(m,2, qold.get(i+1,j,k,m,3) - qold.get(i,j,k,m,3) );
		      Qin.set(m,3, qold.get(i,j,k+1,m,3) - qold.get(i,j,k,m,3) );
		      Qin.set(m,4, qold.get(i,j+1,k,m,3) - qold.get(i,j,k,m,3) );
		    }
		  break;
		  
		case 2:
		  for (int m=1; m<=meqn; m++)
		    {  
		      Qin.set(m,1, qold.get(i,j+1,k,m,1) - qold.get(i,j,k,m,1) );
		    }
		  break;
		}
	      
	      // Project Qin onto left eigenvectors in cell (i,j,k)
	      ProjectLeftEig(2,Aux_ave,Q_ave,Qin,Wout);
	      
	      // Store results in dw_right
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  dwp.set(i,j,k,m,ell, Wout.get(m,ell) );  }		  
	      // --------------------------------------------
	      
	      // ------------------------------------------------
	      // Part II:  dw_left = Ly*( q(i,j,k) - q(i,j-1,k) )
	      // ------------------------------------------------
	      switch( space_order )		  
		{
		case 3:
		  for (int m=1; m<=meqn; m++)
		    {
		      Qin.set(m,1, qold.get(i,j,k,m,1) - qold.get(i,j-1,k,m,1) );
		      Qin.set(m,2, qold.get(i,j,k,m,3) - qold.get(i-1,j,k,m,3) );
		      Qin.set(m,3, qold.get(i,j,k,m,3) - qold.get(i,j,k-1,m,3) );
		      Qin.set(m,4, qold.get(i,j,k,m,3) - qold.get(i,j-1,k,m,3) );
		    }
		  break;
		  
		case 2:
		  for (int m=1; m<=meqn; m++)
		    {  
		      Qin.set(m,1, qold.get(i,j,k,m,1) - qold.get(i,j-1,k,m,1) );
		    }
		  break;
		}
	      
	      // Project Qin onto left eigenvectors in cell (i,j,k)
	      ProjectLeftEig(2,Aux_ave,Q_ave,Qin,Wout);
		
	      // Store results in dw_left
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  dwm.set(i,j,k,m,ell, Wout.get(m,ell) );  }
	      // --------------------------------------------
		
	      // --------------------------------------------
	      // Part III:  wy_cent = Ly*q(i,j,k)
	      // --------------------------------------------
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  Qin.set(m,ell, qold.get(i,j,k,m,nn.get(ell)) );  }	      
	      
	      // Project Qin onto left eigenvectors in cell (i,j,k)
	      ProjectLeftEig(2,Aux_ave,Q_ave,Qin,Wout);
	
	      // Store results in wx_cent
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  w_cent.set(i,j,k,m,ell, Wout.get(m,ell) );  }	      
	      // --------------------------------------------
	    }
      break;
      
    case 3:
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)
	for (int j=(3-mbc); j<=(my+mbc-2); j++)
	  for (int k=(3-mbc); k<=(mz+mbc-2); k++)
	    {
	      iTensor1 nn(4);
	      dTensor1 Aux_ave(maux),Q_ave(meqn);
	      dTensor2 Qin(meqn,ksize),Wout(meqn,ksize);
	      	      
	      // index for x-direction
	      nn.set(1, 4  );
	      nn.set(2, 6  );
	      nn.set(3, 7  );
	      nn.set(4, 10 );
	      
	      // Store q and aux values in temporary arrays
	      for (int m=1; m<=meqn; m++)
		{  Q_ave.set(m, qold.get(i,j,k,m,1) );   }
	      
	      for (int m=1; m<=maux; m++)
		{  Aux_ave.set(m, aux.get(i,j,k,m,1) );  }
	      
	      // ------------------------------------------------
	      // Part I:  dw_right = Lz*( q(i,j,k+1) - q(i,j,k) )
	      // ------------------------------------------------
	      switch( space_order )
		{
		case 3:
		  for (int m=1; m<=meqn; m++)
		    {  
		      Qin.set(m,1, qold.get(i,j,k+1,m,1) - qold.get(i,j,k,m,1) );
		      Qin.set(m,2, qold.get(i+1,j,k,m,4) - qold.get(i,j,k,m,4) );
		      Qin.set(m,3, qold.get(i,j+1,k,m,4) - qold.get(i,j,k,m,4) );
		      Qin.set(m,4, qold.get(i,j,k+1,m,4) - qold.get(i,j,k,m,4) );
		    }
		  break;
		  
		case 2:
		  for (int m=1; m<=meqn; m++)
		    {  
		      Qin.set(m,1, qold.get(i,j,k+1,m,1) - qold.get(i,j,k,m,1) );
		    }
		  break;
		}
	      
	      // Project Qin onto left eigenvectors in cell (i,j,k)
	      ProjectLeftEig(3,Aux_ave,Q_ave,Qin,Wout);
	      
	      // Store results in dw_right
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  dwp.set(i,j,k,m,ell, Wout.get(m,ell) );  }		  
	      // --------------------------------------------
	      
	      // ------------------------------------------------
	      // Part II:  dw_left = Lz*( q(i,j,k) - q(i,j,k-1) )
	      // ------------------------------------------------
	      switch( space_order )		  
		{
		case 3:
		  for (int m=1; m<=meqn; m++)
		    {
		      Qin.set(m,1, qold.get(i,j,k,m,1) - qold.get(i,j,k-1,m,1) );
		      Qin.set(m,2, qold.get(i,j,k,m,4) - qold.get(i-1,j,k,m,4) );
		      Qin.set(m,3, qold.get(i,j,k,m,4) - qold.get(i,j-1,k,m,4) );
		      Qin.set(m,4, qold.get(i,j,k,m,4) - qold.get(i,j,k-1,m,4) );
		    }
		  break;
		  
		case 2:
		  for (int m=1; m<=meqn; m++)
		    {  
		      Qin.set(m,1, qold.get(i,j,k,m,1) - qold.get(i,j,k-1,m,1) );
		    }
		  break;
		}
	      
	      // Project Qin onto left eigenvectors in cell (i,j,k)
	      ProjectLeftEig(3,Aux_ave,Q_ave,Qin,Wout);
		
	      // Store results in dw_left
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  dwm.set(i,j,k,m,ell, Wout.get(m,ell) );  }
	      // --------------------------------------------
		
	      // --------------------------------------------
	      // Part III:  wz_cent = Lz*q(i,j,k)
	      // --------------------------------------------
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  Qin.set(m,ell, qold.get(i,j,k,m,nn.get(ell)) );  }	      
	      
	      // Project Qin onto left eigenvectors in cell (i,j,k)
	      ProjectLeftEig(3,Aux_ave,Q_ave,Qin,Wout);
	
	      // Store results in wx_cent
	      for (int m=1; m<=meqn; m++)
		for (int ell=1; ell<=ksize; ell++)
		  {  w_cent.set(i,j,k,m,ell, Wout.get(m,ell) );  }	      
	      // --------------------------------------------
	    }
      break;
    }
  
}
