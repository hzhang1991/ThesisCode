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
// orders 1--5 implemented.  This routine calls ApplyScalarLimiter in the scalar
// case.
//
// -------------------------------------------------------------------------- //
void ApplyLimiterKrivodonova(dTensorBC4& aux, dTensorBC4& q,
     void (*ProjectRightEig)(int ixy, // Component of flux function (1 or 2)
            const dTensor1& Aux_ave, const dTensor1& Q_ave,   // Riemann states
            const dTensor2& Wvals,    // Characteristic variables
            dTensor2& Qvals           // Conserved variables
            ),
     void (*ProjectLeftEig)(int ixy,  // Component of flux function (1 or 2)
            const dTensor1& Aux_ave, const dTensor1& Q_ave,   // Riemann states
            const dTensor2& Qvals,    // Conserved variables
            dTensor2& Wvals)          // Characteristic variables
    )
{

    // Problem dimensions
    const int   mx = q.getsize(1);
    const int   my = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    // Number of basis functions
    const int space_order = dogParams.get_space_order();
    assert_eq(int((sqrt(1+8*kmax)-1)/2),space_order);
    const int ksize = (space_order==2)?1:3;

    // Nothing to limit in the first-order case
    if(space_order==1)
      { return; }

    // Scalar limiter is implemented for all orders up to 5.  Call this routine
    // instead if this is a scalar problem.
    if( meqn == 1 ) 
    {
        void ApplyScalarLimiter(const dTensorBC4& aux, dTensorBC4& q);
        ApplyScalarLimiter(aux, q);
        return;
    }

    // Storage
    dTensorBC4  wx_cent(mx,my,meqn,ksize,mbc);
    dTensorBC4  wy_cent(mx,my,meqn,ksize,mbc);
    dTensorBC4 dw_right(mx,my,meqn,ksize,mbc);
    dTensorBC4  dw_left(mx,my,meqn,ksize,mbc);
    dTensorBC4    dw_up(mx,my,meqn,ksize,mbc);
    dTensorBC4  dw_down(mx,my,meqn,ksize,mbc);
    void ConvertQtoW(int ixy, 
            const dTensorBC4& aux, 
            const dTensorBC4& qold, 
            dTensorBC4& dwp, 
            dTensorBC4& dwm, 
            dTensorBC4& w_cent,
            void (*ProjectLeftEig)(int,
                const dTensor1&,
                const dTensor1&,
                const dTensor2&,
                dTensor2&));
    void ConvertWtoQ(int ixy, 
            const dTensorBC4& aux, 
            const dTensorBC4& qold,
            const dTensorBC4& win, 
            dTensorBC4& qout,
            void (*ProjectRightEig)(int,
                const dTensor1&,
                const dTensor1&,
                const dTensor2&,
                dTensor2&));

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

    // Moment limiter (highest to lowest)
    // Reminder: this is for 3rd and 2nd-order case only.
    switch ( space_order )
    {

        default:
            unsupported_value_error(space_order);
            break;

        case 2:  // 2nd order in space   

            // Convert to characteristic variables in x-direction
            ConvertQtoW(1,aux,q,dw_right,dw_left,wx_cent,ProjectLeftEig);

            // Convert to characteristic variables in y-direction
            ConvertQtoW(2,aux,q,dw_up,dw_down,wy_cent,ProjectLeftEig);

            // Limit in both the x-direction and the y-direction
#pragma omp parallel for
            for (int i=(3-mbc); i<=(mx+mbc-2); i++)   
            for (int j=(3-mbc); j<=(my+mbc-2); j++)   
            for (int m=1; m<=meqn; m++)
            {
                // x-direction: limit linear term
                {
                    const double        dwp = dw_right.get(i,j,m,1);
                    const double        dwm =  dw_left.get(i,j,m,1);
                    const double     wx_now =  wx_cent.get(i,j,m,1);
                    const double wx_limited = minmod(wx_now, dwp/sq3, dwm/sq3);

                    wx_cent.set(i,j,m,1, wx_limited );
                }

                // y-direction: limit linear term
                {
                    const double dwp        =    dw_up.get(i,j,m,1);
                    const double dwm        =  dw_down.get(i,j,m,1);
                    const double wy_now     =  wy_cent.get(i,j,m,1);
                    const double wy_limited = minmod(wy_now, dwp/sq3, dwm/sq3);

                    wy_cent.set(i,j,m,1, wy_limited );
                }
            }

            // Convert back to conserved variables in x-direction
            ConvertWtoQ(1,aux,q,wx_cent,q,ProjectRightEig);

            // Convert back to conserved variables in y-direction
            ConvertWtoQ(2,aux,q,wy_cent,q,ProjectRightEig);

            break;

        case 3:  // 3rd order in space

            // Convert to characteristic variables in x-direction
            ConvertQtoW(1, aux, q, dw_right, dw_left, wx_cent, ProjectLeftEig);

            // Convert to characteristic variables in y-direction
            ConvertQtoW(2, aux, q, dw_up, dw_down, wy_cent, ProjectLeftEig);

            // Limit
            const double molifier = 0.6;
            const double fac1 = molifier*sq3*osq5; // osq3*osq5
            const double fac2 = molifier*osq3;     // osq3

#pragma omp parallel for
            for (int i=(3-mbc); i<=(mx+mbc-2); i++)   
            for (int j=(3-mbc); j<=(my+mbc-2); j++)       
            for (int m=1; m<=meqn; m++)
            {
                int mflag = 0;
                double dwp,dwm,wx_limited,wx_now,wy_limited,wy_now;

                // ---------------------------------------------------
                // x-direction: limit quadratic term
                dwp        = dw_right.get(i,j,m,2);
                dwm        =  dw_left.get(i,j,m,2);
                wx_now     =  wx_cent.get(i,j,m,3);
                wx_limited = minmod(wx_now, fac1*dwp, fac1*dwm);

                wx_cent.set(i,j,m,3, wx_limited );

                if (fabs(wx_now - wx_limited)>=1.0e-14 || 
                        fabs(wx_limited)<1.0e-14)
                {  mflag = 1;  }
                // ---------------------------------------------------	      

                // ---------------------------------------------------
                // y-direction: limit quadratic term
                dwp        =    dw_up.get(i,j,m,3);
                dwm        =  dw_down.get(i,j,m,3);
                wy_now     =  wy_cent.get(i,j,m,3);
                wy_limited = minmod(wy_now, fac1*dwp, fac1*dwm);

                wy_cent.set(i,j,m,3, wy_limited );

                if (fabs(wy_now - wy_limited)>=1.0e-14 ||
                        fabs(wy_limited)<1.0e-14)
                {  mflag = 1;  }
                // ---------------------------------------------------        

                if (mflag==1)
                {       
                    // ---------------------------------------------------
                    // x-direction: limit mixed term
                    dwp        = dw_right.get(i,j,m,3);
                    dwm        =  dw_left.get(i,j,m,3);
                    wx_now     =  wx_cent.get(i,j,m,2);
                    wx_limited = minmod(wx_now, fac1*dwp, fac1*dwm);

                    wx_cent.set(i,j,m,2, wx_limited );
                    // ---------------------------------------------------

                    // ---------------------------------------------------
                    // x-direction: limit linear term
                    dwp        = dw_right.get(i,j,m,1);
                    dwm        =  dw_left.get(i,j,m,1);
                    wx_now     =  wx_cent.get(i,j,m,1);
                    wx_limited = minmod(wx_now, fac2*dwp, fac2*dwm);

                    wx_cent.set(i,j,m,1, wx_limited );
                    // ---------------------------------------------------

                    // ---------------------------------------------------
                    // y-direction: limit mixed term
                    dwp        =    dw_up.get(i,j,m,2);
                    dwm        =  dw_down.get(i,j,m,2);
                    wy_now     =  wy_cent.get(i,j,m,2);
                    wy_limited = minmod(wy_now, fac1*dwp, fac1*dwm);

                    wy_cent.set(i,j,m,2, wy_limited );
                    // ---------------------------------------------------

                    // ---------------------------------------------------
                    // y-direction: limit linear term
                    dwp        =    dw_up.get(i,j,m,1);
                    dwm        =  dw_down.get(i,j,m,1);
                    wy_now     =  wy_cent.get(i,j,m,1);
                    wy_limited = minmod(wy_now, fac2*dwp, fac2*dwm);

                    wy_cent.set(i,j,m,1, wy_limited );
                    // ---------------------------------------------------
                }
            }

            // Convert back to conserved variables in x-direction
            ConvertWtoQ(1, aux, q, wx_cent, q, ProjectRightEig);

            // Convert back to conserved variables in y-direction
            ConvertWtoQ(2, aux, q, wy_cent, q, ProjectRightEig);

            break;          

    }

}
// -------------------------------------------------------------------------- //
