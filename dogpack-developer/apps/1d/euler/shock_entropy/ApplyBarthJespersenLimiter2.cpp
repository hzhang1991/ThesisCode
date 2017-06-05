#include "stdio.h"
#include "dog_math.h"
#include "constants.h"
#include "tensors.h"
#include "dogdefs.h"

#include "DogParamsCart1.h"
#include "LegendrePolys1d.h"
#include <iostream>
using namespace std;


// This is an experimental limiter based upon "Limiters for Unstructured
// Higher-Order Accurate Solutions of the Euler Equations", Krzysztof Michalak
// and Carl Ollivier-Gooch published in AIAA.
// 
// It is very similar to the Zhang and Shu limiter, but is applied locally,
// where estimates for the min and max of the function is computed using
// deviations in cell averages of the solution.

inline double phi_func( double x )
{ 
    return Min(x,1.0);//(x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
}

void ApplyLimiter(const dTensor2& node, 
		  const dTensorBC3& aux,
		  dTensorBC3& q, 
		  void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
		  void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(4);

    int icheck=0;

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
    dTensorBC2 maxcellAve (mx, meqn,mbc);
    dTensorBC2 mincellAve (mx, meqn,mbc);

    dTensorBC2 edgeval (mx+1, meqn,mbc);

    double dx=1.0/mx;
    double alpha=30.0*pow(dx,1.5);

    double alpha2=0.0;//30.0*pow(dx,1.01);

    double cutoff=pow(dx,0.5);


    //Compute max and min on each cell...this will be done anyway
    //but we will store it for future computation
    for(int i=1-mbc; i <= mx+mbc; i++)
    {for(int me=1; me <= meqn; me++)
    { 
        double max1=-1.0e6;
        double min1=1.0e6;
        for(int mp=1; mp <= mpoints; mp++)
        {
            double qnow = 0.0;   // (q at quad value
            for( int k=1; k <= kmax; k++ )
            {
                qnow += q.get(i,me,k) * phi.get(mp,k);
            }
            max1=Max(max1,qnow);
            min1=Min(min1,qnow);
            if(mp==1){edgeval.set(i,me,edgeval.get(i,me)-qnow);}
            if(mp==mpoints){edgeval.set(i+1,me,qnow);}
        }
        maxcell.set(i,me,max1);
        maxcellAve.set(i,me,q.get(i,me,1));
        mincellAve.set(i,me,q.get(i,me,1));
        mincell.set(i,me,min1);
    }
   }
    //Now compare to neighboring cells....
    for(int i=1; i <= mx; i++)
    {double thetae=1.0;int detector=0;

    for(int me=1; me <= meqn; me++)
    { if(fabs(edgeval.get(i,me))>cutoff || fabs(edgeval.get(i+1,me))>cutoff)
      {detector=1;}
    }
    if(1)//detector==1)
    {for(int me=1; me <= meqn; me++)
    { 
        double Q1=q.get(i,me,1);
        double diff1M=maxcellAve.get(i-1,me)-Q1;
        double diff1m=mincellAve.get(i-1,me)-Q1;
 
        double diff2M=maxcell.get(i,me)-Q1;
        double diff2m=mincell.get(i,me)-Q1;

        double diff3M=maxcellAve.get(i+1,me)-Q1;
        double diff3m=mincellAve.get(i+1,me)-Q1;

        //Find the minimum and maximum deviations from the average value of cell i
        //Note that 0 is from comparing with the ith cell average itself.

        //double diffM=Max(0.0,Max(diff3M,diff1M));
        //double diffm=Min(0.0,Min(diff3m,diff1m));


        //double diffM=Max(3.0*pow(dx,1.5),Max(diff3M,diff1M));
        //double diffm=Min(-3.0*pow(dx,1.5),Min(diff3m,diff1m));
        //double diffM=Max(alpha,Max(diff3M,diff1M));
        //double diffm=Min(-alpha,Min(diff3m,diff1m));
        double diffM=Max(0.0,Max(diff3M,diff1M));
        double diffm=Min(0.0,Min(diff3m,diff1m));
        
        double thetam1,thetaM1; 
        if( fabs( diff2m ) < 1e-14 )
        {
           thetam1 = 1.0;
        }
        else
        {
           thetam1=phi_func(diffm/diff2m);
        } 
        if( fabs( diff2M ) < 1e-14 )
        {
           thetaM1 = 1.0;
        }
        else
        {
           thetaM1=phi_func(diffM/diff2M);
        }
        if (Max(fabs(diff2M),fabs(diff2m))>alpha2){detector=1;}        
        double theta=Min(1.0,Min(thetam1,thetaM1));
        assert_ge( theta, 0.0 );
        assert_le( theta, 1.0 );
        //if(theta<0.1){cout<<i<<" "<<me<<" theta= "<<theta<<endl;exit(1);}
        thetae=Min(theta,thetae);
        
        /*
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,me,k, q.get(i,me,k) * theta );
        }  */
     } 
    }
    if(thetae<1.0 && detector==1)//phi_func(1.0))
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
   //cout<<icheck<<endl;

}

void ApplyLimiter(const dTensor2& node, 
		  const dTensorBC3& aux,
          const dTensorBC3& qold,
		  dTensorBC3& q, 
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
