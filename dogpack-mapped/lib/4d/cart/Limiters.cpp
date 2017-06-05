#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"

// -------------------------------------------------------------------------- //
//
// Limiter based on a modification of the following paper:
//      L. Krivodonova. "Limiters for high-order discontinuous
//      Galerkin methods." J. Comp. Phys., Vol. 226, pg 879-896.
//
// For systems, we have orders 1--3 implemented.  For scalar problems, we have
// orders 1--5 implemented.
//
// -------------------------------------------------------------------------- //
void ApplyLimiterKrivodonova(dTensorBC5& aux, 
                             dTensorBC5& q,
                             void (*ProjectRightEig)(int,
                                                     const dTensor1&,
                                                     const dTensor1&,
                                                     const dTensor2&,
                                                     dTensor2&),
                             void (*ProjectLeftEig)(int,
                                                    const dTensor1&,
                                                    const dTensor1&,
                                                    const dTensor2&,
                                                    dTensor2&))
{
  const int   mx = q.getsize(1);
  const int   my = q.getsize(2);
  const int   mz = q.getsize(3);
  const int meqn = q.getsize(4);
  const int kmax = q.getsize(5);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(4);
  const int space_order = dogParams.get_space_order();
  
    
  if(space_order==1)
    {  return;  }
  
  int msize;
  switch (space_order)
    { 
    case 2:
      msize = 1;
      break;
    case 3:
      msize = 4;
      break;
    }
  
  dTensorBC5  wx_cent(mx,my,mz,meqn,msize,mbc);
  dTensorBC5 dw_right(mx,my,mz,meqn,msize,mbc);
  dTensorBC5  dw_left(mx,my,mz,meqn,msize,mbc);

  dTensorBC5  wy_cent(mx,my,mz,meqn,msize,mbc);
  dTensorBC5  dw_back(mx,my,mz,meqn,msize,mbc);
  dTensorBC5 dw_front(mx,my,mz,meqn,msize,mbc);

  dTensorBC5  wz_cent(mx,my,mz,meqn,msize,mbc);
  dTensorBC5    dw_up(mx,my,mz,meqn,msize,mbc);
  dTensorBC5  dw_down(mx,my,mz,meqn,msize,mbc);

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
					  dTensor2&));

  void ConvertWtoQ(int ixy, 
		   const dTensorBC5& aux, 
		   const dTensorBC5& qold,
		   const dTensorBC5& win, 
		   dTensorBC5& qout,
		   void (*ProjectRightEig)(int,
					   const dTensor1&,
					   const dTensor1&,
					   const dTensor2&,
					   dTensor2&));

  // ----------------------------------------------------------
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

  // Moment limiter (highest to lowest)
  switch ( space_order )
    {
    default:
      unsupported_value_error(space_order);
      break;
      
    case 2:  // 2nd order in space 

      // Convert to characteristic variables in x-direction
      ConvertQtoW(1,aux,q,dw_right,dw_left,wx_cent,ProjectLeftEig);
      
      // Convert to characteristic variables in y-direction
      ConvertQtoW(2,aux,q,dw_back,dw_front,wy_cent,ProjectLeftEig);

      // Convert to characteristic variables in z-direction
      ConvertQtoW(3,aux,q,dw_up,dw_down,wz_cent,ProjectLeftEig);
      
      // Limit in x, y, and z-directions
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)   
	for (int j=(3-mbc); j<=(my+mbc-2); j++)   
	  for (int k=(3-mbc); k<=(mz+mbc-2); k++)   
	    for (int m=1; m<=meqn; m++)
	      {		
		{
		  // ---------------------------------------------------        
		  // x-direction: limit linear term
		  const double dwp        = dw_right.get(i,j,k,m,1);
		  const double dwm        =  dw_left.get(i,j,k,m,1);
		  const double wx_now     =  wx_cent.get(i,j,k,m,1);
		  const double wx_limited = minmod(wx_now, dwp*osq3, dwm*osq3);
		  		  
		  wx_cent.set(i,j,k,m,1, wx_limited );
		  // ---------------------------------------------------        
		}		
		
		{
		  // ---------------------------------------------------        
		  // y-direction: limit linear term
		  const double dwp        =  dw_back.get(i,j,k,m,1);
		  const double dwm        = dw_front.get(i,j,k,m,1);
		  const double wy_now     =  wy_cent.get(i,j,k,m,1);
		  const double wy_limited = minmod(wy_now, dwp*osq3, dwm*osq3);
		  
		  wy_cent.set(i,j,k,m,1, wy_limited );
		  // ---------------------------------------------------        
		}

		{
		  // ---------------------------------------------------        
		  // z-direction: limit linear term
		  const double dwp        =    dw_up.get(i,j,k,m,1);
		  const double dwm        =  dw_down.get(i,j,k,m,1);
		  const double wz_now     =  wz_cent.get(i,j,k,m,1);
		  const double wz_limited = minmod(wz_now, dwp*osq3, dwm*osq3);
		  
		  wz_cent.set(i,j,k,m,1, wz_limited );
		  // ---------------------------------------------------        
		}
	      }
      
      
      // Convert back to conserved variables in x-direction
      ConvertWtoQ(1,aux,q,wx_cent,q,ProjectRightEig);
      
      // Convert back to conserved variables in y-direction
      ConvertWtoQ(2,aux,q,wy_cent,q,ProjectRightEig);
      
      // Convert back to conserved variables in y-direction
      ConvertWtoQ(3,aux,q,wz_cent,q,ProjectRightEig);      
      
      break;
      
    case 3:  // 3rd order in space
      
      // Convert to characteristic variables in x-direction
      ConvertQtoW(1,aux,q,dw_right,dw_left,wx_cent,ProjectLeftEig);

      // Convert to characteristic variables in y-direction
      ConvertQtoW(2,aux,q,dw_back,dw_front,wy_cent,ProjectLeftEig);
      
      // Convert to characteristic variables in z-direction
      ConvertQtoW(3,aux,q,dw_up,dw_down,wz_cent,ProjectLeftEig);
      
      // Limit
      dTensor1 fac(4);
      fac.set(1, 0.6*osq3    );
      fac.set(2, 0.6*sq3*osq5 );
      fac.set(3, 0.6*sq3*osq5 );
      fac.set(4, 0.6*sq3*osq5 );
      
#pragma omp parallel for
      for (int i=(3-mbc); i<=(mx+mbc-2); i++)   
	for (int j=(3-mbc); j<=(my+mbc-2); j++) 
	  for (int k=(3-mbc); k<=(mz+mbc-2); k++)   
	    for (int m=1; m<=meqn; m++)
	      {
		int mflag = 0;
		
		{
		  // ---------------------------------------------------
		  // x-direction: limit quadratic term
		  const double dwp        = dw_right.get(i,j,k,m,4);
		  const double dwm        =  dw_left.get(i,j,k,m,4);
		  const double wx_now     =  wx_cent.get(i,j,k,m,4);
		  const double wx_limited = minmod(wx_now, fac.get(4)*dwp, fac.get(4)*dwm);
		  
		  wx_cent.set(i,j,k,m,4, wx_limited );
		  
		  if (fabs(wx_now - wx_limited)>=1.0e-14 || 
		      fabs(wx_limited)<1.0e-14)
		    {  mflag = 1;  }
		  // ---------------------------------------------------              
		}

		{
		  // ---------------------------------------------------
		  // y-direction: limit quadratic term
		  const double dwp        =  dw_back.get(i,j,k,m,4);
		  const double dwm        = dw_front.get(i,j,k,m,4);
		  const double wy_now     =  wy_cent.get(i,j,k,m,4);
		  const double wy_limited = minmod(wy_now, fac.get(4)*dwp, fac.get(4)*dwm);
		  
		  wy_cent.set(i,j,k,m,4, wy_limited );
		  
		  if (fabs(wy_now - wy_limited)>=1.0e-14 ||
		      fabs(wy_limited)<1.0e-14)
		    {  mflag = 1;  }
		  // ---------------------------------------------------        
		}
		
		{
		  // ---------------------------------------------------
		  // z-direction: limit quadratic term
		  const double dwp        =    dw_up.get(i,j,k,m,4);
		  const double dwm        =  dw_down.get(i,j,k,m,4);
		  const double wz_now     =  wz_cent.get(i,j,k,m,4);
		  const double wz_limited = minmod(wz_now, fac.get(4)*dwp, fac.get(4)*dwm);
		  
		  wz_cent.set(i,j,k,m,4, wz_limited );
		  
		  if (fabs(wz_now - wz_limited)>=1.0e-14 ||
		      fabs(wz_limited)<1.0e-14)
		    {  mflag = 1;  }
		  // ---------------------------------------------------        
		}

		if (mflag==1)
		  {  
		    for (int ell=1; ell<=3; ell++)
		      {
			// ---------------------------------------------------
			// x-direction: limit linear and mixed terms
			const double dwp        = dw_right.get(i,j,k,m,ell);
			const double dwm        =  dw_left.get(i,j,k,m,ell);
			const double wx_now     =  wx_cent.get(i,j,k,m,ell);
			const double wx_limited = minmod(wx_now, fac.get(ell)*dwp, fac.get(ell)*dwm);
			
			wx_cent.set(i,j,k,m,ell, wx_limited );
			// ---------------------------------------------------              
		      }
		    
		    for (int ell=1; ell<=3; ell++)
		      {			
			// ---------------------------------------------------
			// y-direction: limit linear and mixed terms
			const double dwp        =  dw_back.get(i,j,k,m,ell);
			const double dwm        = dw_front.get(i,j,k,m,ell);
			const double wy_now     =  wy_cent.get(i,j,k,m,ell);
			const double wy_limited = minmod(wy_now, fac.get(ell)*dwp, fac.get(ell)*dwm);
			
			wy_cent.set(i,j,k,m,ell, wy_limited );			
			// ---------------------------------------------------        
		      }

		    for (int ell=1; ell<=3; ell++)
		      {
			// ---------------------------------------------------
			// z-direction: limit mixed term
			const double dwp        =    dw_up.get(i,j,k,m,ell);
			const double dwm        =  dw_down.get(i,j,k,m,ell);
			const double wz_now     =  wz_cent.get(i,j,k,m,ell);
			const double wz_limited = minmod(wz_now, fac.get(ell)*dwp, fac.get(ell)*dwm);
			
			wz_cent.set(i,j,k,m,ell, wz_limited );
			// ---------------------------------------------------        
		      }
		  }
	      }
      
      // Convert back to conserved variables in x-direction
      ConvertWtoQ(1,aux,q,wx_cent,q,ProjectRightEig);
      
      // Convert back to conserved variables in y-direction
      ConvertWtoQ(2,aux,q,wy_cent,q,ProjectRightEig);
      
      // Convert back to conserved variables in y-direction
      ConvertWtoQ(3,aux,q,wz_cent,q,ProjectRightEig);

      break;          
    }
  
}
// -------------------------------------------------------------------------- //
