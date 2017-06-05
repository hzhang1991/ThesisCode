#include "dogdefs.h"
#include "dog_math.h"

void ApplyLimiter(double t, double dt, dTensorBC4& aux, dTensorBC4& q)
{
    // -------------------------------------------------------------  
    // Limiter based on a modification of the following paper:
    //      L. Krivodonova. "Limiters for high-order discontinuous
    //      Galerkin methods." J. Comp. Phys., Vol. 226, pg 879-896.
    // ------------------------------------------------------------- 
    int i,j,m,ksize;
    int   mx = q.getsize(1);
    int   my = q.getsize(2);
    int meqn = q.getsize(3);
    int kmax = q.getsize(4);
    int mbc  = q.getmbc();
    int maux = aux.getsize(2);
    int method1 = int((sqrt(1+8*kmax)-1)/2);
    if (method1==2)
    {  ksize = 1;  }
    else
    {  ksize = 3;  }
    dTensorBC4  wx_cent(mx,my,meqn,ksize,mbc);
    dTensorBC4  wy_cent(mx,my,meqn,ksize,mbc);
    dTensorBC4 dw_right(mx,my,meqn,ksize,mbc);
    dTensorBC4  dw_left(mx,my,meqn,ksize,mbc);
    dTensorBC4    dw_up(mx,my,meqn,ksize,mbc);
    dTensorBC4  dw_down(mx,my,meqn,ksize,mbc);
    void ConvertQtoW(int ixy, const dTensorBC4& aux, const dTensorBC4& qold, 
                     dTensorBC4& dwp, dTensorBC4& dwm, dTensorBC4& w_cent);
    void ConvertWtoQ(int ixy, const dTensorBC4& aux, const dTensorBC4& qold,
                     const dTensorBC4& win, dTensorBC4& qout);
    double dwm,dwp;
    double wx_now, wx_limited, wy_now, wy_limited;
    int mflag;

    // Moment limiter (highest to lowest)
    switch ( method1 )
    {
        case 2:  // 2nd order in space   

	  // Limit in both the x-direction and the y-direction	  
	  wx_cent.set(i,j,1,1, q.get(i,j,1,2) );
	  wy_cent.set(i,j,1,1, q.get(i,j,1,2) );
	  
          for (i=(3-mbc); i<=(mx+mbc-2); i++)
	  {
	      for (j=(3-mbc); j<=(my+mbc-2); j++)
	      {
		  // x-direction: limit linear term
		  dwp        = q.get(i+1,j,1,1) - 2.0*q.get(i,j,1,1)
		    + q.get(i-1,j,1,1);
		  dwm        = q.get(i,j,1,1) - 2.0*q.get(i-1,j,1,1) 
		    + q.get(i-2,j,1,1);
		  wx_now     = q.get(i,j,1,2) - q.get(i-1,j,1,2);
		  wx_limited = minmod(wx_now, dwp/sq3, dwm/sq3);
		
		  wx_cent.set(i,j,1,1, wx_cent.get(i-1,j,1,1) + wx_limited );
		  
		  // y-direction: limit linear term
		  dwp        = q.get(i,j+1,1,1) - 2.0*q.get(i,j,1,1)
		    + q.get(i,j-1,1,1);
		  dwm        = q.get(i,j,1,1) - 2.0*q.get(i,j-1,1,1)
		    + q.get(i,j-2,1,1);
		  wy_now     = q.get(i,j,1,3) - q.get(i,j-1,1,3);
		  wy_limited = minmod(wy_now, dwp/sq3, dwm/sq3);
		  
		  wy_cent.set(i,j,1,1, wy_cent.get(i,j-1,1,1) + wy_limited );
              }
          }

          break;

        case 3:  // 3rd order in space
	  
	  break;          
    }
 
}
