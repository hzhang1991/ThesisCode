#include "tensors.h"
#include "DogParams.h"

// Convert conserved variables to characteristic variables
void ConvertQtoW(int ixy, 
		 const dTensorBC4* aux,
		 const dTensorBC4& qold, 
                 dTensorBC4& dwp, 
		 dTensorBC4& dwm, 
		 dTensorBC4& w_cent,
		 void (*ProjectLeftEig)(int,
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
  const int ksize = w_cent.getsize(4);
  const int space_order = dogParams.get_space_order();
  
  // ----------------------------------------------------------
  // Key: storage of characteristic variables
  //       
  //      wx_cent(i,j,me,1)  = Lx * q(i,j,me,2)
  //      wx_cent(i,j,me,2)  = Lx * q(i,j,me,4)
  //      wx_cent(i,j,me,3)  = Lx * q(i,j,me,5)
  //
  //      dw_right(i,j,me,1) = Lx * (q(i+1,j,me,1)-q(i,j,me,1))
  //      dw_right(i,j,me,2) = Lx * (q(i+1,j,me,2)-q(i,j,me,2))
  //      dw_right(i,j,me,3) = Lx * (q(i,j+1,me,2)-q(i,j,me,2))
  //
  //      dw_left(i,j,me,1)  = Lx * (q(i,j,me,1)-q(i-1,j,me,1))
  //      dw_left(i,j,me,2)  = Lx * (q(i,j,me,2)-q(i-1,j,me,2))
  //      dw_left(i,j,me,3)  = Lx * (q(i,j,me,2)-q(i,j-1,me,2))
  //
  //      wy_cent(i,j,me,1)  = Ly * q(i,j,me,3)
  //      wy_cent(i,j,me,2)  = Ly * q(i,j,me,4)
  //      wy_cent(i,j,me,3)  = Ly * q(i,j,me,6)
  //
  //      dw_up(i,j,me,1)    = Ly * (q(i,j+1,me,1)-q(i,j,me,1))
  //      dw_up(i,j,me,2)    = Ly * (q(i+1,j,me,3)-q(i,j,me,3))
  //      dw_up(i,j,me,3)    = Ly * (q(i,j+1,me,3)-q(i,j,me,3))
  //
  //      dw_down(i,j,me,1)  = Ly * (q(i,j,me,1)-q(i,j-1,me,1))
  //      dw_down(i,j,me,2)  = Ly * (q(i,j,me,3)-q(i-1,j,me,3))
  //      dw_down(i,j,me,3)  = Ly * (q(i,j,me,3)-q(i,j-1,me,3))
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
	  {
	    iTensor1 nx(3);
            dTensor1 Aux_ave(maux),Q_ave(meqn);
            dTensor2 Qin(meqn,ksize),Wout(meqn,ksize);

	    // index for x-direction
	    nx.set(1, 2 );
	    nx.set(2, 4 );
	    nx.set(3, 5 );
	    
	        // Store q and aux values in temporary arrays
	    for (int m=1; m<=meqn; m++)
	      {  Q_ave.set(m, qold.get(i,j,m,1) );   }
 	    
	    for (int m=1; m<=maux; m++)
	      {  Aux_ave.set(m, aux->get(i,j,m,1) );  }
	    	    
	    // ############################################
	    //   X-DIRECTION
	    // ############################################
	    
	    // --------------------------------------------
	    // Part I:  dw_right = Lx*( q(i+1,j) - q(i,j) )
	    // --------------------------------------------
	    switch( space_order )
	      {
	      case 3:
		for (int m=1; m<=meqn; m++)
		  {  
		    Qin.set(m,1, qold.get(i+1,j,m,1) - qold.get(i,j,m,1) );
		    Qin.set(m,2, qold.get(i+1,j,m,2) - qold.get(i,j,m,2) );
		    Qin.set(m,3, qold.get(i,j+1,m,2) - qold.get(i,j,m,2) );
		  }
		break;

	      case 2:
		for (int m=1; m<=meqn; m++)
		  {  
		    Qin.set(m,1, qold.get(i+1,j,m,1) - qold.get(i,j,m,1) );
		  }
		break;
	      }

	    // Project Qin onto left eigenvectors in cell (i,j)
	    ProjectLeftEig(1,Aux_ave,Q_ave,Qin,Wout);
	    
	    // Store results in dw_right
	    for (int m=1; m<=meqn; m++)
	      for (int k=1; k<=ksize; k++)
		{  dwp.set(i,j,m,k, Wout.get(m,k) );  }	    
	    // --------------------------------------------
	    
	    // --------------------------------------------
	    // Part II:  dw_left = Lx*( q(i,j) - q(i-1,j) )
	    // --------------------------------------------
	    switch( space_order )
	      {
	      case 3:
		for (int m=1; m<=meqn; m++)
		  {
		    Qin.set(m,1, qold.get(i,j,m,1) - qold.get(i-1,j,m,1) );
		    Qin.set(m,2, qold.get(i,j,m,2) - qold.get(i-1,j,m,2) );
		    Qin.set(m,3, qold.get(i,j,m,2) - qold.get(i,j-1,m,2) );
		  }
		break;

	      case 2:
		for (int m=1; m<=meqn; m++)
		  {  
		    Qin.set(m,1, qold.get(i,j,m,1) - qold.get(i-1,j,m,1) );
		  }
		break;
	      }

	    // Project Qin onto left eigenvectors in cell (i,j)
	    ProjectLeftEig(1,Aux_ave,Q_ave,Qin,Wout);
	    
	    // Store results in dw_left
	    for (int m=1; m<=meqn; m++)
	      for (int k=1; k<=ksize; k++)
		{  dwm.set(i,j,m,k, Wout.get(m,k) );  }		
	    // --------------------------------------------
	    
	    // --------------------------------------------
	    // Part III:  wx_cent = Lx*q(i,j)
	    // --------------------------------------------
	    for (int m=1; m<=meqn; m++) 
	      for (int k=1; k<=ksize; k++)
		{  Qin.set(m,k, qold.get(i,j,m,nx.get(k)) );  }		

	    // Project Qin onto left eigenvectors in cell (i,j)
	    ProjectLeftEig(1,Aux_ave,Q_ave,Qin,Wout);
	
	    // Store results in wx_cent
	    for (int m=1; m<=meqn; m++)
	      for (int k=1; k<=ksize; k++)
		{  w_cent.set(i,j,m,k, Wout.get(m,k) );  }
	    // --------------------------------------------
	  }
      break;

    case 2:
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)
	for (int j=(3-mbc); j<=(my+mbc-2); j++)
	  {	        
	    iTensor1 ny(3);
            dTensor1 Aux_ave(maux),Q_ave(meqn);
            dTensor2 Qin(meqn,ksize),Wout(meqn,ksize);

	    // index for y-direction
	    ny.set(1, 3 );
	    ny.set(2, 4 );
	    ny.set(3, 6 );
	    
	    // Store q and aux values in temporary arrays
	    for (int m=1; m<=meqn; m++)
	      {  Q_ave.set(m, qold.get(i,j,m,1) );   }
 	    
	    for (int m=1; m<=maux; m++)
	      {  Aux_ave.set(m, aux->get(i,j,m,1) );  }
	    
	    // ############################################
	    //   Y-DIRECTION
	    // ############################################
	    
	    // --------------------------------------------
	    // Part I:  dw_up = Ly*( q(i,j+1) - q(i,j) )
	    // --------------------------------------------
	    switch( space_order )
	      {
	      case 3:
		for (int m=1; m<=meqn; m++)
		  {
		    Qin.set(m,1, qold.get(i,j+1,m,1) - qold.get(i,j,m,1) );
		    Qin.set(m,2, qold.get(i+1,j,m,3) - qold.get(i,j,m,3) );
		    Qin.set(m,3, qold.get(i,j+1,m,3) - qold.get(i,j,m,3) );
		  }
		break;
		
	      case 2:
		for (int m=1; m<=meqn; m++)
		  {
		    Qin.set(m,1, qold.get(i,j+1,m,1) - qold.get(i,j,m,1) );			
		  }
		break;
	      }

	    // Project Qin onto left eigenvectors in cell (i,j)
	    ProjectLeftEig(2,Aux_ave,Q_ave,Qin,Wout);
	    
	    // Store results in dw_up
	    for (int m=1; m<=meqn; m++)
	      for (int k=1; k<=ksize; k++)
		{  dwp.set(i,j,m,k, Wout.get(m,k) );  }		
	    // --------------------------------------------
	    
	    // --------------------------------------------
	    // Part II:  dw_down = Ly*( q(i,j) - q(i,j-1) )
	    // --------------------------------------------
	    switch( space_order )
	      {
	      case 3:
		for (int m=1; m<=meqn; m++)
		  {  
		    Qin.set(m,1, qold.get(i,j,m,1) - qold.get(i,j-1,m,1) );
		    Qin.set(m,2, qold.get(i,j,m,3) - qold.get(i-1,j,m,3) );
		    Qin.set(m,3, qold.get(i,j,m,3) - qold.get(i,j-1,m,3) );
		  }
		break;

	      case 2:
		for (int m=1; m<=meqn; m++)
		  {  
		    Qin.set(m,1, qold.get(i,j,m,1) - qold.get(i,j-1,m,1) );
		  }
		break;
	      }	      
	    
	    // Project Qin onto left eigenvectors in cell (i,j)
	    ProjectLeftEig(2,Aux_ave,Q_ave,Qin,Wout);
	    
	    // Store results in dw_down
	    for (int m=1; m<=meqn; m++)
	      for (int k=1; k<=ksize; k++)
		{  dwm.set(i,j,m,k, Wout.get(m,k) );  }
	    // --------------------------------------------
	    
	    // --------------------------------------------
	    // Part III:  wy_cent = Ly*q(i,j)
	    // --------------------------------------------
	    for (int m=1; m<=meqn; m++)
	      for (int k=1; k<=ksize; k++)
		{  Qin.set(m,k, qold.get(i,j,m,ny.get(k)) );  }	      
	    
	    // Project Qin onto left eigenvectors in cell (i,j)
	    ProjectLeftEig(2,Aux_ave,Q_ave,Qin,Wout);
	    
	    // Store results in wy_cent
	    for (int m=1; m<=meqn; m++)
	      for (int k=1; k<=ksize; k++)
		{  w_cent.set(i,j,m,k, Wout.get(m,k) );  }
	    // --------------------------------------------
	  }
      break;
    }
}


// Convert conserved variables to characteristic variables
void ConvertQtoW(int ixy, const dTensorBC4& aux, const dTensorBC4& qold, 
                 dTensorBC4& dwp, dTensorBC4& dwm, dTensorBC4& w_cent,
		 void (*ProjectLeftEig)(int,const dTensor1&,
					const dTensor1&,const dTensor2&,
					dTensor2&))
{
  void ConvertQtoW(int ixy, const dTensorBC4* aux, const dTensorBC4& qold, 
		   dTensorBC4& dwp, dTensorBC4& dwm, dTensorBC4& w_cent,
		   void (*ProjectLeftEig)(int,const dTensor1&,
					  const dTensor1&,const dTensor2&,
					  dTensor2&));

  ConvertQtoW(ixy,&aux,qold,dwp,dwm,w_cent,ProjectLeftEig);  
}
