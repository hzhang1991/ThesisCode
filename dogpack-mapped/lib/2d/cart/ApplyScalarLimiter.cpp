#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"

void ApplyScalarLimiter(const dTensorBC4& aux, dTensorBC4& q)
{
    // -------------------------------------------------------------  
    // Limiter based on a modification of the following paper:
    //      L. Krivodonova. "Limiters for high-order discontinuous
    //      Galerkin methods." J. Comp. Phys., Vol. 226, pg 879-896.
    // ------------------------------------------------------------- 
    const int   mx = q.getsize(1);
    const int   my = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(2);
    const int method1 = int((sqrt(1+8*kmax)-1)/2);
//  const int ksize;
    const int m = 1;   // equation number is ALWAYS 1 here
    const double sq9 = 3.0e0;

    assert_printf(meqn==1, "Scalar Limiter only supports for scalar equations: "
            "meqn = %d", meqn);

    // TODO - play around with this indicator!
    const int morder = dogParams.get_space_order();
    const int kless = int((morder*(morder-1))/2);
    const double s0   = 0.5e0/pow(double(morder),4);
    const double Qtol = 1e-2;

    // Copy info into qold
    dTensorBC4 qold(mx,my,1,kmax,2,2); 

#pragma omp parallel for
    for (int i=(1-mbc); i<=(mx+mbc); i++)
      for (int j=(1-mbc); j<=(my+mbc); j++)
	for (int k=1; k<=kmax; k++)
	  {
	    qold.set(i,j,1,k, q.get(i,j,1,k) );
	  }

    // Moment limiter (highest to lowest)
#pragma omp parallel for 
    for (int i=(3-mbc); i<=(mx+mbc-2); i++)
      for (int j=(3-mbc); j<=(my+mbc-2); j++)      
        {

//////// Smoothness indicator /////////////
//          double Q2sum = 0.0;

//          for (int k=1; k<=kless; k++)
//          {  Q2sum = Q2sum + pow(q.get(i,j,m,k),2);  }

//          double Q2kmax = 0.0;

//          for (int k=(kless+1); k<=kmax; k++)
//          {  Q2kmax = Q2kmax + pow(q.get(i,j,m,k),2);  }

//          Q2sum = Q2sum + Q2kmax;

//          double se;
//          if (Q2kmax>1.0e-10)
//          {  se =  (Q2kmax/Q2sum);  }
//          else
//          {  se = 0.0e0;  }

//          if( se > s0 )
            if( q.get(i,j,m,1) < Qtol )
            { // apply the limiter 

                int mflag = 0; // flag to indicate if we need to continue limiting
                double qxp, qxm, qx_now, qx_limited;
                double qyp, qym, qy_now, qy_limited;

                switch( method1 )
                {
                    case 5:  // limit "5th order" terms: Q^{(11)} -- Q^{(15)}

                        // reset the flag to 0
                        mflag = 0;

                        /////////////////////////////////////////////////////////
                        //                      x^4 term                       //
                        /////////////////////////////////////////////////////////

                        // x-direction
                        qxp        = qold.get(i+1,j,m, 9) - qold.get(i  ,j,m,9);
                        qxm        = qold.get(i,  j,m, 9) - qold.get(i-1,j,m,9);
                        qx_now     = q.get(i,j,m,14);
                        //                    cout << "i = " << i << "   j = " << j << endl;
                        qx_limited = minmod(qx_now, qxp*sq7/sq9, qxm*sq7/sq9);

                        if( fabs(qx_limited - qx_now ) >= EPSILON ||
                                fabs(qx_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,14, qx_limited );

                        // y-direction
                        //
                        //  -- does not exist for this edge case!


                        /////////////////////////////////////////////////////////
                        //                    x^3 * y term                     //
                        /////////////////////////////////////////////////////////

                        // x-direction
                        qxp        = qold.get(i+1,j,m,7) - qold.get(i  ,j,m,7);
                        qxm        = qold.get(i,  j,m,7) - qold.get(i-1,j,m,7);
                        qx_now     = q.get(i,j,m,11);
                        qx_limited = minmod(qx_now, qxp*sq7/sq9, qxm*sq7/sq9);

                        if( fabs(qx_limited - qx_now ) >= EPSILON ||
                                fabs(qx_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,11, qx_limited );

                        // y-direction 
                        qyp        = qold.get(i,j+1,m,9) - qold.get(i,j,  m,9);
                        qym        = qold.get(i,j,  m,9) - qold.get(i,j-1,m,9);
                        qy_now     = q.get(i,j,m,11);
                        qy_limited = minmod(qy_now, qyp*sq7/sq9, qym*sq7/sq9);

                        if( fabs(qy_limited - qy_now ) >= EPSILON ||
                                fabs(qy_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,11, qy_limited );

                        /////////////////////////////////////////////////////////
                        //                    x^2 y^2 term                     //
                        /////////////////////////////////////////////////////////

                        // x-direction
                        qxp        = qold.get(i+1,j,m,8) - qold.get(i  ,j,m,8);
                        qxm        = qold.get(i,  j,m,8) - qold.get(i-1,j,m,8);
                        qx_now     = q.get(i,j,m,13);
                        qx_limited = minmod(qx_now, qxp*sq7/sq9, qxm*sq7/sq9);

                        if( fabs(qx_limited - qx_now ) >= EPSILON ||
                                fabs(qx_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,13, qx_limited );

                        // y-direction 
                        qyp        = qold.get(i,j+1,m,7) - qold.get(i,j,  m,7);
                        qym        = qold.get(i,j,  m,7) - qold.get(i,j-1,m,7);
                        qy_now     = q.get(i,j,m,13);
                        qy_limited = minmod(qy_now, qyp*sq7/sq9, qym*sq7/sq9);

                        if( fabs(qy_limited - qy_now ) >= EPSILON ||
                                fabs(qy_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,13, qy_limited );

                        /////////////////////////////////////////////////////////
                        //                    x y^3 term                       //
                        /////////////////////////////////////////////////////////

                        // x-direction
                        qxp        = qold.get(i+1,j,m,10) - qold.get(i  ,j,m,10);
                        qxm        = qold.get(i,  j,m,10) - qold.get(i-1,j,m,10);
                        qx_now     = q.get(i,j,m,12);
                        qx_limited = minmod(qx_now, qxp*sq7/sq9, qxm*sq7/sq9);

                        if( fabs(qx_limited - qx_now ) >= EPSILON ||
                                fabs(qx_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,12, qx_limited );

                        // y-direction 
                        qyp        = qold.get(i,j+1,m,8) - qold.get(i,j,  m,4);
                        qym        = qold.get(i,j,  m,8) - qold.get(i,j-1,m,4);
                        qy_now     = q.get(i,j,m,12);
                        qy_limited = minmod(qy_now, qyp*sq7/sq9, qym*sq7/sq9);

                        if( fabs(qy_limited - qy_now ) >= EPSILON ||
                                fabs(qy_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,12, qy_limited );


                        /////////////////////////////////////////////////////////
                        //                      y^4 term                       //
                        /////////////////////////////////////////////////////////

                        // y-direction 
                        qyp        = qold.get(i,j+1,m,10) - qold.get(i,j,  m,10);
                        qym        = qold.get(i,j,  m,10) - qold.get(i,j-1,m,10);
                        qy_now     = q.get(i,j,m,15);
                        qy_limited = minmod(qy_now, qyp*sq7/sq9, qym*sq7/sq9);

                        if( fabs(qy_limited - qy_now ) >= EPSILON ||
                                fabs(qy_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,15, qy_limited );

                        // no need to check x-direction for this edge case

                        // do we need to continue?
                        if( mflag == 0 ){ break; }


                    case 4:  // limit "4th order" terms: Q^{(7)} -- Q^{(10)}

                        // reset the flag to 0
                        mflag = 0;

                        /////////////////////////////////////////////////////////
                        //                      x^3 term                       //
                        /////////////////////////////////////////////////////////

                        // x-direction
                        qxp        = qold.get(i+1,j,m,5) - qold.get(i  ,j,m,5);
                        qxm        = qold.get(i,  j,m,5) - qold.get(i-1,j,m,5);
                        qx_now     = q.get(i,j,m,9);
                        qx_limited = minmod(qx_now, qxp*sq5/sq7, qxm*sq5/sq7);

                        if( fabs(qx_limited - qx_now ) >= EPSILON ||
                                fabs(qx_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,9, qx_limited );

                        // y-direction
                        //
                        //  -- does not exist for this edge case!


                        /////////////////////////////////////////////////////////
                        //                      xxy term                       //
                        /////////////////////////////////////////////////////////

                        // x-direction
                        qxp        = qold.get(i+1,j,m,4) - qold.get(i  ,j,m,4);
                        qxm        = qold.get(i,  j,m,4) - qold.get(i-1,j,m,4);
                        qx_now     = q.get(i,j,m,7);
                        qx_limited = minmod(qx_now, qxp*sq5/sq7, qxm*sq5/sq7);

                        if( fabs(qx_limited - qx_now ) >= EPSILON ||
                                fabs(qx_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,7, qx_limited );

                        // y-direction 
                        qyp        = qold.get(i,j+1,m,5) - qold.get(i,j,  m,5);
                        qym        = qold.get(i,j,  m,5) - qold.get(i,j-1,m,5);
                        qy_now     = q.get(i,j,m,7);
                        qy_limited = minmod(qy_now, qyp*sq5/sq7, qym*sq5/sq7);

                        if( fabs(qy_limited - qy_now ) >= EPSILON ||
                                fabs(qy_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,7, qy_limited );

                        /////////////////////////////////////////////////////////
                        //                      xyy term                       //
                        /////////////////////////////////////////////////////////

                        // x-direction
                        qxp        = qold.get(i+1,j,m,6) - qold.get(i  ,j,m,6);
                        qxm        = qold.get(i,  j,m,6) - qold.get(i-1,j,m,6);
                        qx_now     = q.get(i,j,m,8);
                        qx_limited = minmod(qx_now, qxp*sq5/sq7, qxm*sq5/sq7);

                        if( fabs(qx_limited - qx_now ) >= EPSILON ||
                                fabs(qx_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,8, qx_limited );

                        // y-direction 
                        qyp        = qold.get(i,j+1,m,4) - qold.get(i,j,  m,4);
                        qym        = qold.get(i,j,  m,4) - qold.get(i,j-1,m,4);
                        qy_now     = q.get(i,j,m,8);
                        qy_limited = minmod(qy_now, qyp*sq5/sq7, qym*sq5/sq7);

                        if( fabs(qy_limited - qy_now ) >= EPSILON ||
                                fabs(qy_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,8, qy_limited );

                        /////////////////////////////////////////////////////////
                        //                      yyy term                       //
                        /////////////////////////////////////////////////////////

                        // y-direction 
                        qyp        = qold.get(i,j+1,m,6) - qold.get(i,j,  m,6);
                        qym        = qold.get(i,j,  m,6) - qold.get(i,j-1,m,6);
                        qy_now     = q.get(i,j,m,10);
                        qy_limited = minmod(qy_now, qyp*sq5/sq7, qym*sq5/sq7);

                        if( fabs(qy_limited - qy_now ) >= EPSILON ||
                                fabs(qy_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,10, qy_limited );

                        // no need to check x-direction for this edge case

                        // do we need to continue?
                        if( mflag == 0 ){ break; }

                    case 3:  // limit "3rd order" terms: Q^{(4)} -- Q^{(6)}

                        // reset the flag to 0
                        mflag = 0;

                        /////////////////////////////////////////////////////////
                        //                      x^2 term                       //
                        /////////////////////////////////////////////////////////

                        // x-direction
                        qxp        = qold.get(i+1,j,m,2) - qold.get(i  ,j,m,2);
                        qxm        = qold.get(i,  j,m,2) - qold.get(i-1,j,m,2);
                        qx_now     = q.get(i,j,m,5);
                        qx_limited = minmod(qx_now, qxp*sq3/sq5, qxm*sq3/sq5);

                        if( fabs(qx_limited - qx_now ) >= EPSILON ||
                                fabs(qx_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,5, qx_limited );

                        // y-direction
                        //
                        //  -- does not exist for this edge case!


                        /////////////////////////////////////////////////////////
                        //                      xy term                        //
                        /////////////////////////////////////////////////////////

                        // x-direction
                        qxp        = qold.get(i+1,j,m,3) - qold.get(i  ,j,m,3);
                        qxm        = qold.get(i,  j,m,3) - qold.get(i-1,j,m,3);
                        qx_now     = q.get(i,j,m,4);
                        qx_limited = minmod(qx_now, qxp*sq3/sq5, qxm*sq3/sq5);

                        if( fabs(qx_limited - qx_now ) >= EPSILON ||
                                fabs(qx_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,4, qx_limited );

                        // y-direction 
                        qyp        = qold.get(i,j+1,m,2) - qold.get(i,j,  m,2);
                        qym        = qold.get(i,j,  m,2) - qold.get(i,j-1,m,2);
                        qy_now     = q.get(i,j,m,4);
                        qy_limited = minmod(qy_now, qyp*sq3/sq5, qym*sq3/sq5);

                        if( fabs(qy_limited - qy_now ) >= EPSILON ||
                                fabs(qy_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,4, qy_limited );

                        /////////////////////////////////////////////////////////
                        //                      yy term                        //
                        /////////////////////////////////////////////////////////
                        // y-direction 
                        qyp        = qold.get(i,j+1,m,3) - qold.get(i,j,  m,3);
                        qym        = qold.get(i,j,  m,3) - qold.get(i,j-1,m,3);
                        qy_now     = q.get(i,j,m,6);
                        qy_limited = minmod(qy_now, qyp*sq3/sq5, qym*sq3/sq5);

                        if( fabs(qy_limited - qy_now ) >= EPSILON ||
                                fabs(qy_limited) < EPSILON ) { mflag = 1; }
                        q.set(i,j,m,6, qy_limited );

                        // no need to check x-direction for this edge case

                        // do we need to continue?
                        if( mflag == 0 ){ break; }

                    case 2:  // limit "2nd order" terms:  Q^{(2)} -- Q^{(3)}

                        // x-direction: limit linear term
                        qxp        = qold.get(i+1,j,m,1) - qold.get(i  ,j,m,1);
                        qxm        = qold.get(i,  j,m,1) - qold.get(i-1,j,m,1);
                        qx_now     = q.get(i,j,m,2);
                        qx_limited = minmod(qx_now, qxp/sq3, qxm/sq3);

                        q.set(i,j,m,2, qx_limited );

                        // y-direction: limit linear term
                        qyp        = qold.get(i,j+1,m,1) - qold.get(i,j,  m,1);
                        qym        = qold.get(i,j,  m,1) - qold.get(i,j-1,m,1);
                        qy_now     = q.get(i,j,m,3);
                        qy_limited = minmod(qy_now, qyp/sq3, qym/sq3);

                        q.set(i,j,m,3, qy_limited );

                        break;

                }
            }

        }
}

