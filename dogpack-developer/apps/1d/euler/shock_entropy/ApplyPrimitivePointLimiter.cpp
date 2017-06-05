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
    return (x<1.1)*Min(1.0/1.1*x,1.0)+(x>=1.1)*1.0;//*(x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
}

void ApplyLimiter(const dTensor2& node, const dTensorBC3& aux, dTensorBC3& q, 
		  void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
		  void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(4);

    double gamma=1.4;

    int icheck=0;

    // Set up the "normal" quadrature points
    const int mpoints = kmax;
    const int MAX_KMAX =   6;
    dTensor2 phi (mpoints, MAX_KMAX);
    dTensor1 spts(mpoints);

    //Compare at the GLL Nodes//
    void setGaussLobattoPoints1d(dTensor1& x1d);
    setGaussLobattoPoints1d( spts  );
    evaluateLegendrePolys( spts, phi );

    //Store the max values on each cell
    dTensorBC2 maxcell (mx, meqn,mbc);
    dTensorBC2 mincell (mx, meqn,mbc);
    dTensorBC2 avecell (mx, meqn,mbc);


    double dx=1.0/mx;
    double alpha=200.0*pow(dx,1.0);
    //double alpha=1.0*pow(dx,1.5);
 
    //Only used for Shock-Detecting schemes
    double alpha2=0.0;//30.0*pow(dx,0.6);

    //Compute max and min on each cell...this will be done anyway
    //but we will store it for future computation
    for(int i=1-mbc; i <= mx+mbc; i++)
    { 
        double maxrho=-1.0e6;
        double minrho=1.0e6;
        double maxu=-1.0e6;
        double minu=1.0e6;
        double maxv=-1.0e6;
        double minv=1.0e6;
        double maxw=-1.0e6;
        double minw=1.0e6;
        double maxP=-1.0e6;
        double minP=1.0e6;
         
        for(int mp=1; mp <= mpoints; mp++)
        {
            double rhonow = 0.0;
            double unow   = 0.0;
            double vnow   = 0.0;
            double wnow   = 0.0;
            double Pnow   = 0.0;

            double mxnow  = 0.0;
            double mynow  = 0.0;
            double mznow  = 0.0;
            double enow   = 0.0;

            for( int k=1; k <= kmax; k++ )
            {
                rhonow+=q.get(i,1,k) * phi.get(mp,k);
                mxnow +=q.get(i,2,k) * phi.get(mp,k);
                mynow +=q.get(i,3,k) * phi.get(mp,k);
                mznow +=q.get(i,4,k) * phi.get(mp,k);
                enow  +=q.get(i,5,k) * phi.get(mp,k);  
            }

            unow=mxnow/q.get(i,1,1);
            vnow=mynow/q.get(i,1,1);
            wnow=mznow/q.get(i,1,1);
            Pnow=(gamma-1.0e0)*(enow/q.get(i,1,1)-0.5*(q.get(i,2,1)/q.get(i,1,1)*unow));
       
            maxrho=Max(maxrho,rhonow);
            minrho=Min(minrho,rhonow);
            maxu=Max(maxu,unow);
            minu=Min(minu,unow);
            maxv=Max(maxv,vnow);
            minv=Min(minv,vnow);
            maxw=Max(maxw,wnow);
            minw=Min(minw,wnow);
            maxP=Max(maxP,Pnow);
            minP=Min(minP,Pnow);

        }
        maxcell.set(i,1,maxrho);
        mincell.set(i,1,minrho);
        avecell.set(i,1,q.get(i,1,1));
        maxcell.set(i,2,maxu);
        mincell.set(i,2,minu);
        avecell.set(i,2,q.get(i,2,1)/q.get(i,1,1));
        maxcell.set(i,3,maxv);
        mincell.set(i,3,minv);
        avecell.set(i,3,q.get(i,3,1)/q.get(i,1,1));
        maxcell.set(i,4,maxw);
        mincell.set(i,4,minw);
        avecell.set(i,4,q.get(i,4,1)/q.get(i,1,1));
        maxcell.set(i,5,maxP);
        mincell.set(i,5,minP);
        avecell.set(i,5,(gamma-1.0e0)*(q.get(i,5,1)/q.get(i,1,1)-0.5*q.get(i,2,1)/q.get(i,1,1)*q.get(i,2,1)/q.get(i,1,1)));
    }


    //Now compare to neighboring cells....
    for(int i=1; i <= mx; i++)
    {double thetae=1.0;
    
    {
        double Q1rho=avecell.get(i,1);//q.get(i,me,1);
        double diff1Mrho=maxcell.get(i-1,1)-Q1rho;
        double diff1mrho=mincell.get(i-1,1)-Q1rho;
 
        double diff2Mrho=maxcell.get(i,1)-Q1rho;
        double diff2mrho=mincell.get(i,1)-Q1rho;

        double diff3Mrho=maxcell.get(i+1,1)-Q1rho;
        double diff3mrho=mincell.get(i+1,1)-Q1rho;
        double diffMrho=Max(alpha,Max(diff3Mrho,diff1Mrho));
        double diffmrho=Min(-alpha,Min(diff3mrho,diff1mrho));
        
        double thetam1rho,thetaM1rho;
        //Attempt to compute theta bounding our data within bounds. 
        if( fabs( diff2mrho ) < 1e-14 )
        {
           thetam1rho = 1.0;
        }
        else
        {
           thetam1rho=phi_func(diffmrho/diff2mrho);
        } 
        if( fabs( diff2Mrho ) < 1e-14 )
        {
           thetaM1rho = 1.0;
        }
        else
        {
           thetaM1rho=phi_func(diffMrho/diff2Mrho);
        }

        double thetar=Min(1.0,Min(thetam1rho,thetaM1rho));
        assert_ge( thetar, 0.0 );
        assert_le( thetar, 1.0 );
        
        double Q1=avecell.get(i,2);//q.get(i,me,1);
        double diff1M=maxcell.get(i-1,2)-Q1;
        double diff1m=mincell.get(i-1,2)-Q1;
 
        double diff2M=maxcell.get(i,2)-Q1;
        double diff2m=mincell.get(i,2)-Q1;

        double diff3M=maxcell.get(i+1,2)-Q1;
        double diff3m=mincell.get(i+1,2)-Q1;
        double diffM=Max(alpha,Max(diff3M,diff1M));
        double diffm=Min(-alpha,Min(diff3m,diff1m));
        
        double thetam1,thetaM1;
        //Attempt to compute theta bounding our data within bounds. 
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
        double thetamx=Min(1.0,Min(thetam1,thetaM1));
        assert_ge( thetamx, 0.0 );
        assert_le( thetamx, 1.0 );
        
        
        Q1=avecell.get(i,5);//q.get(i,me,1);
        diff1M=maxcell.get(i-1,5)-Q1;
        diff1m=mincell.get(i-1,5)-Q1;
 
        diff2M=maxcell.get(i,5)-Q1;
        diff2m=mincell.get(i,5)-Q1;

        diff3M=maxcell.get(i+1,5)-Q1;
        diff3m=mincell.get(i+1,5)-Q1;
        diffM=Max(alpha,Max(diff3M,diff1M));
        diffm=Min(-alpha,Min(diff3m,diff1m));
        
        thetam1,thetaM1;
        //Attempt to compute theta bounding our data within bounds. 
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
        double thetamp=Min(1.0,Min(thetam1,thetaM1));
        assert_ge( thetamp, 0.0 );
        assert_le( thetamp, 1.0 );

        thetae=Min(thetar,thetae);
        thetae=Min(thetamx,thetae);
        thetae=Min(thetamp,thetae);
 
        /*for( int k=2; k <= kmax; k++ )
        {
            q.set(i,1,k, q.get(i,1,k) * thetar );
        }
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,5,k,(q.get(i,5,k)-0.5*(q.get(i,2,1)/q.get(i,1,1))*q.get(i,2,k))* thetamp + 0.5*(q.get(i,2,1)/q.get(i,1,1))*q.get(i,2,k) * thetamx );
        }
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,2,k, q.get(i,2,k) * thetamx );
        }*/

    for(int me=1; me <= meqn; me++)
    { 
        Q1=avecell.get(i,me);//q.get(i,me,1);
        diff1M=maxcell.get(i-1,me)-Q1;
        diff1m=mincell.get(i-1,me)-Q1;
 
        diff2M=maxcell.get(i,me)-Q1;
        diff2m=mincell.get(i,me)-Q1;

        diff3M=maxcell.get(i+1,me)-Q1;
        diff3m=mincell.get(i+1,me)-Q1;
        

        /*
        double diff1M=q.get(i-1,me,1)-Q1;
        double diff1m=q.get(i-1,me,1)-Q1;

        double diff2M=maxcell.get(i,me)-Q1;
        double diff2m=mincell.get(i,me)-Q1;

        double diff3M=q.get(i+1,me,1)-Q1;
        double diff3m=q.get(i+1,me,1)-Q1;
        */       

        //Find the minimum and maximum deviations from the average value of cell i
        //Note that 0 is from comparing with the ith cell average itself.

        //double diffM=Max(0.0,Max(diff3M,diff1M));
        //double diffm=Min(0.0,Min(diff3m,diff1m));


        //double diffM=Max(3.0*pow(dx,1.5),Max(diff3M,diff1M));
        //double diffm=Min(-3.0*pow(dx,1.5),Min(diff3m,diff1m));

        
        //Add a tolerance alpha
        diffM=Max(alpha,Max(diff3M,diff1M));
        diffm=Min(-alpha,Min(diff3m,diff1m));
        
        thetam1,thetaM1;
        //Attempt to compute theta bounding our data within bounds. 
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
        double theta=Min(1.0,Min(thetam1,thetaM1));
        assert_ge( theta, 0.0 );
        assert_le( theta, 1.0 );
        //if(theta<0.1){cout<<i<<" "<<me<<" theta= "<<theta<<endl;exit(1);}
        thetae=Min(theta,thetae);
        
        
        /*for( int k=2; k <= kmax; k++ )
        {
            q.set(i,me,k, q.get(i,me,k) * theta );
        } */
      
    }
   }
    //Limit high order entries
   if(thetae<1.0)//phi_func(1.0))
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

//Ignore this next limiter, old version

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
