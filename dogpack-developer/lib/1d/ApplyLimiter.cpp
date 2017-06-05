#include <iostream>
#include "dogdefs.h"
#include "stdio.h"
#include "dog_math.h"
#include "constants.h"
#include "tensors.h"
#include "DogParamsCart1.h"
#include "LegendrePolys1d.h"
using namespace std;

// This is a wrapper function designed to select which limiter should be
// applied
void ApplyLimiter(const dTensor2& node, const dTensorBC3& aux, dTensorBC3& q, 
    void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
    void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{

    // TODO - introduce a switch that selects between the various options of
    // limiting that we have! We will likely need to dump each option into its
    // own local library (e.g. shallow water vs. Euler vs. MHD) because each
    // set of applications will likely need (or want) to override this.  -DS

    // This is the new "theta" limiter.
    void ApplyThetaLimiter(const dTensor2& node, const dTensorBC3& aux, dTensorBC3& q, 
        void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
        void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&));
    ApplyThetaLimiter( node, aux, q, ProjectLeftEig, ProjectRightEig );

}

// -------------------------------------------------------------  
// Moment limiter based on the following paper:
//
//      L. Krivodonova. "Limiters for high-order discontinuous
//      Galerkin methods." J. Comp. Phys., Vol. 226, pg 879-896.
//
// See also: ApplyThetaLimiter
// -------------------------------------------------------------
void ApplyLimiterKrivodonova(const dTensor2& node, const dTensorBC3& aux, dTensorBC3& q, 
        void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
        void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{

    void ConvertQtoW(const dTensorBC3&,const dTensorBC3&,dTensorBC3&,dTensorBC3&,dTensorBC3&,
            void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&));
    void ConvertWtoQ(const dTensorBC3&,const dTensorBC3&,const dTensorBC3&,dTensorBC3&,
            void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&));

    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   kmax = q.getsize(3);
    const int    mbc = q.getmbc();
    const int   maux = aux.getsize(2);
    dTensorBC3   w_cent(melems,meqn,kmax-1,mbc);
    dTensorBC3 dw_right(melems,meqn,kmax-1,mbc);
    dTensorBC3  dw_left(melems,meqn,kmax-1,mbc);

    ///////////////////////////////////////////////////////////////////////////
    // If these are placed inside the for loops, then we can run this is
    // parallel - DS
    ////double dw1_right, dw1_left, w2_now, w2_limited;
    ////double dw2_right, dw2_left, w3_now, w3_limited;
    ////double dw3_right, dw3_left, w4_now, w4_limited;
    ////double dw4_right, dw4_left, w5_now, w5_limited;
    ///////////////////////////////////////////////////////////////////////////

    const double dx = node.get(2,1)-node.get(1,1);

    // Convert to characteristic variables
    ConvertQtoW(aux,q,dw_right,dw_left,w_cent,ProjectLeftEig);

    // Moment limiter (highest to lowest)
    switch ( kmax )
    {
        // ****************************************************************
        case 2:  // 2nd order in space     

#pragma omp parallel for
            for (int i=(3-mbc); i<=(melems+mbc-2); i++)
                for (int m=1; m<=meqn; m++)
                {
                    // limit highest moment
                    double dw1_right  = dw_right.get(i,m,1);
                    double dw1_left   = dw_left.get(i,m,1);
                    double w2_now     = w_cent.get(i,m,1);
                    double w2_limited = minmod(w2_now, dw1_right/sq3, dw1_left/sq3);

                    w_cent.set(i,m,1, w2_limited );       
                }

            break;

            // ****************************************************************
        case 3:  // 3rd order in space     

#pragma omp parallel for    
            for (int i=(3-mbc); i<=(melems+mbc-2); i++)      
                for (int m=1; m<=meqn; m++)
                {
                    // limit highest moment
                    double dw2_right  = dw_right.get(i,m,2);
                    double dw2_left   = dw_left.get(i,m,2);
                    double w3_now     = w_cent.get(i,m,2);
                    double w3_limited =minmod(w3_now, (sq3/sq5)*dw2_right, (sq3/sq5)*dw2_left);

                    w_cent.set(i,m,2, w3_limited );

                    // only limit next order if highest moment was limited
                    if (fabs(w3_now - w3_limited)>1.0e-14 || 
                            fabs(w3_limited)<1.0e-14)
                    {
                        double dw1_right  = dw_right.get(i,m,1);
                        double dw1_left   = dw_left.get(i,m,1);
                        double w2_now     = w_cent.get(i,m,1);
                        double w2_limited = minmod(w2_now, dw1_right/sq3, dw1_left/sq3);

                        w_cent.set(i,m,1, w2_limited );
                    }          
                }
            break;

            // ****************************************************************
        case 4:  // 4th order in space     

#pragma omp parallel for
            for (int i=(3-mbc); i<=(melems+mbc-2); i++)     
                for (int m=1; m<=meqn; m++)
                {
                    double dw3_right  = dw_right.get(i,m,3);
                    double dw3_left   = dw_left.get(i,m,3);
                    double w4_now     = w_cent.get(i,m,3);
                    double w4_limited =minmod(w4_now, (sq5/sq7)*dw3_right, 
                            (sq5/sq7)*dw3_left);

                    w_cent.set(i,m,3, w4_limited );

                    // only limit next order if previous moment was limited
                    if (fabs(w4_now - w4_limited)>1.0e-14 || 
                            fabs(w4_limited)<1.0e-14)
                    {              
                        double dw2_right  = dw_right.get(i,m,2);
                        double dw2_left   = dw_left.get(i,m,2);
                        double w3_now     = w_cent.get(i,m,2);
                        double w3_limited =minmod(w3_now, (sq3/sq5)*dw2_right, 
                                (sq3/sq5)*dw2_left);

                        w_cent.set(i,m,2, w3_limited );

                        // only limit next order if previous moment was limited
                        if (fabs(w3_now - w3_limited)>1.0e-14 || 
                                fabs(w3_limited)<1.0e-14)
                        {
                            double dw1_right  = dw_right.get(i,m,1);
                            double dw1_left   = dw_left.get(i,m,1);
                            double w2_now     = w_cent.get(i,m,1);
                            double w2_limited = minmod(w2_now, dw1_right/sq3, 
                                    dw1_left/sq3);

                            w_cent.set(i,m,1, w2_limited );
                        }        
                    }  
                }      
            break;

            // ****************************************************************
        case 5:  // 5th order in space     

#pragma omp parallel for
            for (int i=(3-mbc); i<=(melems+mbc-2); i++)      
                for (int m=1; m<=meqn; m++)
                {
                    double dw4_right  = dw_right.get(i,m,4);
                    double dw4_left   = dw_left.get(i,m,4);
                    double w5_now     = w_cent.get(i,m,4);
                    double w5_limited =minmod(w5_now, (sq7/3.0)*dw4_right, 
                            (sq7/3.0)*dw4_left);

                    w_cent.set(i,m,4, w5_limited );

                    // only limit next order if previous moment was limited
                    if (fabs(w5_now - w5_limited)>1.0e-14 || 
                            fabs(w5_limited)<1.0e-14)
                    {
                        double dw3_right  = dw_right.get(i,m,3);
                        double dw3_left   = dw_left.get(i,m,3);
                        double w4_now     = w_cent.get(i,m,3);
                        double w4_limited =minmod(w4_now, (sq5/sq7)*dw3_right, 
                                (sq5/sq7)*dw3_left);

                        w_cent.set(i,m,3, w4_limited );

                        // only limit next order if previous moment was limited
                        if (fabs(w4_now - w4_limited)>1.0e-14 || 
                                fabs(w4_limited)<1.0e-14)
                        {    
                            double dw2_right  = dw_right.get(i,m,2);
                            double dw2_left   = dw_left.get(i,m,2);
                            double w3_now     = w_cent.get(i,m,2);
                            double w3_limited =minmod(w3_now, (sq3/sq5)*dw2_right, 
                                    (sq3/sq5)*dw2_left);

                            w_cent.set(i,m,2, w3_limited );

                            // only limit next order if previous moment was limited
                            if (fabs(w3_now - w3_limited)>1.0e-14 || 
                                    fabs(w3_limited)<1.0e-14)
                            {
                                double dw1_right  = dw_right.get(i,m,1);
                                double dw1_left   = dw_left.get(i,m,1);
                                double w2_now     = w_cent.get(i,m,1);
                                double w2_limited = minmod(w2_now, dw1_right/sq3, 
                                        dw1_left/sq3);

                                w_cent.set(i,m,1, w2_limited );
                            }        
                        }  
                    }
                }      
            break;

            // ****************************************************************
        case 6:  // 6th order in space     

#pragma omp parallel for
            for (int i=(3-mbc); i<=(melems+mbc-2); i++)      
                for (int m=1; m<=meqn; m++)
                {      
                    double dw5_right  = dw_right.get(i,m,5);
                    double dw5_left   = dw_left.get(i,m,5);
                    double w6_now     = w_cent.get(i,m,5);
                    double w6_limited =minmod(w6_now, (3.0/sq11)*dw5_right, 
                            (3.0*sq11)*dw5_left);

                    w_cent.set(i,m,5, w6_limited );

                    if (fabs(w6_now - w6_limited)>1.0e-14 || 
                            fabs(w6_limited)<1.0e-14)
                    {    
                        double dw4_right  = dw_right.get(i,m,4);
                        double dw4_left   = dw_left.get(i,m,4);
                        double w5_now     = w_cent.get(i,m,4);
                        double w5_limited =minmod(w5_now, (sq7/3.0)*dw4_right, 
                                (sq7/3.0)*dw4_left);

                        w_cent.set(i,m,4, w5_limited );

                        // only limit next order if previous moment was limited
                        if (fabs(w5_now - w5_limited)>1.0e-14 || 
                                fabs(w5_limited)<1.0e-14)
                        {
                            double dw3_right  = dw_right.get(i,m,3);
                            double dw3_left   = dw_left.get(i,m,3);
                            double w4_now     = w_cent.get(i,m,3);
                            double w4_limited =minmod(w4_now, (sq5/sq7)*dw3_right, 
                                    (sq5/sq7)*dw3_left);

                            w_cent.set(i,m,3, w4_limited );

                            // only limit next order if previous moment was limited
                            if (fabs(w4_now - w4_limited)>1.0e-14 || 
                                    fabs(w4_limited)<1.0e-14)
                            {      
                                double dw2_right  = dw_right.get(i,m,2);
                                double dw2_left   = dw_left.get(i,m,2);
                                double w3_now     = w_cent.get(i,m,2);
                                double w3_limited =minmod(w3_now, (sq3/sq5)*dw2_right, 
                                        (sq3/sq5)*dw2_left);

                                w_cent.set(i,m,2, w3_limited );

                                // only limit next order if previous moment was limited
                                if (fabs(w3_now - w3_limited)>1.0e-14 || 
                                        fabs(w3_limited)<1.0e-14)
                                {
                                    double dw1_right  = dw_right.get(i,m,1);
                                    double dw1_left   = dw_left.get(i,m,1);
                                    double w2_now     = w_cent.get(i,m,1);
                                    double w2_limited = minmod(w2_now, dw1_right/sq3, 
                                            dw1_left/sq3);

                                    w_cent.set(i,m,1, w2_limited );
                                }        
                            }  
                        }
                    }      
                }
            break;
    }

    // Convert characteristic variables to conserved variables
    ConvertWtoQ(aux,q,w_cent,q,ProjectRightEig);

}

// This is an experimental limiter based upon "Limiters for Unstructured
// Higher-Order Accurate Solutions of the Euler Equations", Krzysztof Michalak
// and Carl Ollivier-Gooch published in AIAA.
// 
// It is very similar to the Zhang and Shu limiter, but is applied locally,
// where estimates for the min and max of the function is computed using
// deviations in cell averages of the solution.

inline double phi_func( double x )
{ 
    // return (x<1.1)*Min(1.0/1.1*x,1.0)+(x>=1.1)*1.0;//*(x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
    return Min( 1.0, x/1.1 );
}

// This is the new "theta" limiter.
void ApplyThetaLimiter(const dTensor2& node, const dTensorBC3& aux, dTensorBC3& q, 
    void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
    void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(4);

    // Grid information, and tunable parameter
    const double dx = dogParamsCart1.get_dx();
    const double alpha = 30.0*pow(dx,1.1);
 
    //Only used for Shock-Detecting schemes
    const double alpha2 = 0.0;  // 30.0*pow(dx,0.6);

    // Set up the "normal" quadrature points and evaluate the basis functions
    // there.  This will be used to expedite evaluating the solution at each
    // of these points.  
    const int mpoints = kmax;
    const int MAX_KMAX =   6;
    dTensor2 phi (mpoints, MAX_KMAX);
    dTensor1 spts(mpoints);

    // Compare at the GLL Nodes//
    void setGaussLobattoPoints1d(dTensor1& x1d);
    setGaussLobattoPoints1d( spts  );
    evaluateLegendrePolys( spts, phi );

    // Storage for the max (and min) values on each cell
    dTensorBC2 maxcell (mx, meqn,mbc);
    dTensorBC2 mincell (mx, meqn,mbc);
    dTensorBC2 edgeval (mx+1, meqn,mbc);

    // Compute max and min on each cell...this will be done anyway
    // but we will store it for future computation
    for(int i=1-mbc; i <= mx+mbc; i++)
    {
        for(int me=1; me <= meqn; me++)
        { 
            double max1 = -1.0e6;
            double min1 =  1.0e6;
            for(int mp=1; mp <= mpoints; mp++)
            {
                double qnow = 0.0;   // (q this quad pts)
                for( int k=1; k <= kmax; k++ )
                {
                    qnow += q.get(i,me,k) * phi.get(mp,k);
                }
                max1 = Max(max1,qnow);
                min1 = Min(min1,qnow);

                if(mp==1)
                {edgeval.set(i,me,edgeval.get(i,me)-qnow);}

                if(mp==mpoints)
                {edgeval.set(i+1,me,qnow);}

            }
            maxcell.set(i,me,max1);
            mincell.set(i,me,min1);
        }
    }

    // Shock detecting parameter 
    const double cutoff = pow(dx, 0.5);

    //Now compare to neighboring cells....
    int icheck = 0;
    for(int i=1; i <= mx; i++)
    {
        double thetae = 1.0;
        int detector  = 0;

        //For Edge Detector from krivodonova 2004
        for(int me=1; me <= meqn; me++)
        { 
            if(fabs(edgeval.get(i,me))>cutoff || fabs(edgeval.get(i+1,me))>cutoff)
            {detector=1;}
        }
        if(detector==1)
        {
            for(int me=1; me <= meqn; me++)
            { 
                double Q1=q.get(i,me,1);
                double diff1M=maxcell.get(i-1,me)-Q1;
                double diff1m=mincell.get(i-1,me)-Q1;
         
                double diff2M=maxcell.get(i,me)-Q1;
                double diff2m=mincell.get(i,me)-Q1;

                double diff3M=maxcell.get(i+1,me)-Q1;
                double diff3m=mincell.get(i+1,me)-Q1;

                //Find the minimum and maximum deviations from the average value of cell i
                //Note that 0 is from comparing with the ith cell average itself.

                //double diffM=Max(0.0,Max(diff3M,diff1M));
                //double diffm=Min(0.0,Min(diff3m,diff1m));

                //double diffM=Max(3.0*pow(dx,1.5),Max(diff3M,diff1M));
                //double diffm=Min(-3.0*pow(dx,1.5),Min(diff3m,diff1m));
                
                //Add a tolerance alpha
                double diffM=Max(alpha,Max(diff3M,diff1M));
                double diffm=Min(-alpha,Min(diff3m,diff1m));
                
                double thetam1,thetaM1;
                //Attempt to compute theta bounding our data within bounds. 
                if( fabs( diff2m ) < 1e-14 )
                { thetam1 = 1.0; }
                else
                { thetam1=phi_func(diffm/diff2m); } 
                if( fabs( diff2M ) < 1e-14 )
                { thetaM1 = 1.0; }
                else
                { thetaM1=phi_func(diffM/diff2M); }

                if( Max( fabs(diff2M), fabs(diff2m) ) > alpha2 )
                { detector = 1; }        

                double theta=Min(1.0,Min(thetam1,thetaM1));
                assert_ge( theta, 0.0 );
                assert_le( theta, 1.0 );

                //if(theta<0.1){cout<<i<<" "<<me<<" theta= "<<theta<<endl;exit(1);}
                thetae=Min(theta,thetae);
                
            }

        }

        // Limit high-order entries
        //
        // Here, we use the same value of theta for all terms!
        //
        if( thetae < 1.0 && detector == 1 )
        {
            for(int me=1; me <= meqn; me++)
            { 
                for( int k=2; k <= kmax; k++ )
                {
                     q.set(i,me,k, q.get(i,me,k) * thetae );
                }  
            }
           icheck+=1;
        }
    }

}

/*
//Ignore this next limiter, This is the old version
void ApplyLimiter(const dTensor2& node, 
    const dTensorBC3& aux, const dTensorBC3& qold, dTensorBC3& q, 
    void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
    void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(4);

    // Set up the "normal" quadrature points
    const int mpoints = kmax;
    const int MAX_KMAX =   6;
    dTensor2 phi (mpoints, MAX_KMAX);
    dTensor1 spts(mpoints);

    // Use Legendre Points for plotting purposes //
    void setGaussLobattoPoints1d(dTensor1& x1d);
    setGaussLobattoPoints1d( spts  );
    evaluateLegendrePolys( spts, phi );

    dTensorBC2 maxcell (mx, meqn,mbc);
    dTensorBC2 mincell (mx, meqn,mbc);

    //Compute max and min on each cell...this will be done anyway
    //but we will store it for future computation
    for(int i=1-mbc; i <= mx+mbc; i++)
    for(int me=1; me <= meqn; me++)
    { 
        double max1=-1.0e16;
        double min1=1.0e16;
        for(int mp=1; mp <= mpoints; mp++)
        {
            double qnow = 0.0;   // (q at quad value
            for( int k=1; k <= kmax; k++ )
            {
                qnow += q.get(i,me,k) * phi.get(mp,k);
            }
            max1=Max(max1,qnow);
            min1=Min(min1,qnow);
        }
        maxcell.set(i,me,max1);
        mincell.set(i,me,min1);
    }
    //Now compare to neighboring cells....
    for(int i=1; i <= mx; i++)
    for(int me=1; me <= meqn; me++)
    { 
        double Q1=q.get(i,me,1);
        double diff1M=maxcell.get(i-1,me)-Q1;
        double diff1m=mincell.get(i-1,me)-Q1;
 
        double diff2M=maxcell.get(i,me)-Q1;
        double diff2m=mincell.get(i,me)-Q1;

        double diff3M=maxcell.get(i+1,me)-Q1;
        double diff3m=mincell.get(i+1,me)-Q1;

        //Find the minimum and maximum deviations from the average value of cell i
        //Note that 0 is from comparing with the ith cell average itself.

        double diffM=Max(0.0,Max(diff3M,diff1M));
        double diffm=Min(0.0,Min(diff3m,diff1m));
        double thetam1,thetaM1; 
        if( fabs( diff2m ) < 1e-12 )
        {
           thetam1 = 1.0;
        }
        else
        {
           thetam1=phi_func(diffm/diff2m);
        } 
        if( fabs( diff2M ) < 1e-12 )
        {
           thetaM1 = 1.0;
        }
        else
        {
           thetaM1=phi_func(diffM/diff2M);
        }

        double theta=Min(thetam1,thetaM1);
        assert_ge( theta, 0.0 );
        assert_le( theta, 1.0 );
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,me,k, q.get(i,me,k) * theta );
        }        
    }

}
*/
