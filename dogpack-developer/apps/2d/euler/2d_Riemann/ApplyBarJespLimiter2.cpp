#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "Legendre2d.h"
#include "ApplyPosLimiter.h"
#include <iostream>
using namespace std;

// Routine that applies a modified Barth-Jespersen limiter to higher order moments of
// the conserved variables.  This limiter is locally applied, and conserves
// total mass by not adjusting cell averages.
//
// See: K. Michalak and C. Ollivier-Gooch, "Limiters for Unstructured Higher-Order
// Accurate Solutions of the Euler Equations"
inline double phi_func(double x)
{
   return Min(1.0,x);//(x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
}
void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
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

    double Min(double, double);
    const int   mx = q.getsize(1);
    const int   my = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(2);
    const int space_order = dogParams.get_space_order();

    // ------------------------------------------------ //
    // Number of points where we want to check solution
    //
    // Change this quantity as needed to check more points
    //
    // ------------------------------------------------ //
    const int mpoints   =  4+4*space_order+space_order*space_order;

    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //

    // Quadrature points
    dTensor2 spts(mpoints,2);
    dTensor1 w1d(space_order);
    void SetPositivePointsGauss(const int& space_order, dTensor2& spts,dTensor1& w1d);
    SetPositivePointsGauss(space_order,spts,w1d);

    // Basis functions evaluated at each point
    dTensor2 phi(mpoints,kmax);
    SetLegendreAtPoints(space_order,spts,phi);
    double dx=1.0/mx;

    //A matrix to store max and min values on each cell. 
    //We need this to see if a cell's point values fall
    //within the range of the values of its neighbors. 
    //This will determine bounds on each cell.

    dTensorBC3 MaxVal(mx,my,meqn,mbc); MaxVal.setall(-1.0e18);
    dTensorBC3 MinVal(mx,my,meqn,mbc); MinVal.setall(1.0e18);
    dTensorBC4 edgeval(mx+1,my+1,meqn,4,mbc);edgeval.setall(0.0e0);

    dTensorBC3 MaxValAve(mx,my,meqn,mbc); MaxValAve.setall(-1.0e18);
    dTensorBC3 MinValAve(mx,my,meqn,mbc); MinValAve.setall(1.0e18);

    // -------------------------------------------------------------- //
    //                                                                //
    // Loop over each cell and every quadrature point on each cell    //
    // I do not see a clever way of avoiding this right now 
    // unless we want to compute approximate max and min...           //
    // -------------------------------------------------------------- //
#pragma omp parallel for
    for(int i=1-mbc; i <= mx+mbc; i++)
    for(int j=1-mbc; j <= my+mbc; j++)
    for(int me=1; me <= meqn; me++)
    {
        double max1 = MaxVal.get(i,j,me);
        double min1 = MinVal.get(i,j,me);
        int s=1;
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double qnow = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                qnow += q.get(i,j,me,k) * phi.get(mp,k);
            }
            max1=Max(max1,qnow);
            min1=Min(min1,qnow);

            if (mp>4 && mp<=4+space_order)
            {                             

                edgeval.set(i,j,me,1,edgeval.get(i,j,me,1)+qnow*w1d.get(s));
                s=s+1;
                if(s==space_order){s=1;}
            }
            else if(mp>4+space_order && mp<=4+2*space_order)
            {
                edgeval.set(i,j,me,2,edgeval.get(i,j,me,2)+qnow*w1d.get(s));
                s=s+1;

                if(s==space_order){s=1;}
            }
            else if(mp>4+2*space_order && mp<=4+3*space_order)
            {
                edgeval.set(i,j,me,3,edgeval.get(i,j,me,3)+qnow*w1d.get(s));
                s=s+1;
                if(s==space_order){s=1;}
            } 
            else if (mp>4+3*space_order && mp<= 4*space_order)
            {
                edgeval.set(i,j,me,4,edgeval.get(i,j,me,4)+qnow*w1d.get(s));
                s=s+1;
                if(s==space_order){s=1;}
            } 

        }

        MaxVal.set(i,j,me,max1);
        MinVal.set(i,j,me,min1);

        MaxValAve.set(i,j,me,q.get(i,j,me,1));
        MinValAve.set(i,j,me,q.get(i,j,me,1));


    }
    
    double alpha=10.0*pow(dx,1.5);
    double cutoff=pow(dx,0.5);

#pragma omp parallel for
    for(int i=1; i <= mx; i++)
    for(int j=1; j <= my; j++)
    { double thetae=1.0;int detector=0;
    for(int me=1; me <= meqn; me++)
    { if(fabs(edgeval.get(i-1,j,3,me)-edgeval.get(i,j,1,me))>cutoff || fabs(edgeval.get(i,j,1,me)-edgeval.get(i-1,j,3,me))>cutoff || fabs(edgeval.get(i,j+1,2,me)-edgeval.get(i,j,4,me))>cutoff || fabs(edgeval.get(i,j-1,4,me)-edgeval.get(i,j,2,me))>cutoff)
      {detector=1;}
    }
    if(detector==1)
    {
    for(int me=1; me <= meqn; me++)
    {
        double Q1=q.get(i,j,me,1);

       //Find the deviation from the max and min values
       //on neighbouring cells from the average value
       //on our current cell.


        /*
        double difflM=MaxVal.get(i-1,j,me)-Q1;
        double difflm=MinVal.get(i-1,j,me)-Q1;

        double diffrM=MaxVal.get(i+1,j,me)-Q1;
        double diffrm=MinVal.get(i+1,j,me)-Q1;

        double diffuM=MaxVal.get(i,j+1,me)-Q1;
        double diffum=MinVal.get(i,j+1,me)-Q1;

        double diffdM=MaxVal.get(i,j-1,me)-Q1;
        double diffdm=MinVal.get(i,j-1,me)-Q1; 
    
        double diffM=Max(alpha,Max(Max(difflM,diffrM),Max(diffuM,diffdM)));
        double diffm=Min(-alpha,Min(Min(difflm,diffrm),Min(diffum,diffdm)));
 
        //cout<<"here "<<diffM<<endl;

        double diffdlM=MaxVal.get(i-1,j-1,me)-Q1;
        double diffdlm=MinVal.get(i-1,j-1,me)-Q1;

        double diffdrM=MaxVal.get(i+1,j-1,me)-Q1;
        double diffdrm=MinVal.get(i+1,j-1,me)-Q1;

        double diffurM=MaxVal.get(i+1,j+1,me)-Q1;
        double diffurm=MinVal.get(i+1,j+1,me)-Q1;

        double diffulM=MaxVal.get(i-1,j+1,me)-Q1;
        double diffulm=MinVal.get(i-1,j+1,me)-Q1;
        */

        
        double difflM=MaxValAve.get(i-1,j,me)-Q1;
        double difflm=MinValAve.get(i-1,j,me)-Q1;

        double diffrM=MaxValAve.get(i+1,j,me)-Q1;
        double diffrm=MinValAve.get(i+1,j,me)-Q1;

        double diffuM=MaxValAve.get(i,j+1,me)-Q1;
        double diffum=MinValAve.get(i,j+1,me)-Q1;

        double diffdM=MaxValAve.get(i,j-1,me)-Q1;
        double diffdm=MinValAve.get(i,j-1,me)-Q1;

        double diffM=Max(alpha,Max(Max(difflM,diffrM),Max(diffuM,diffdM)));
        double diffm=Min(-alpha,Min(Min(difflm,diffrm),Min(diffum,diffdm)));

        //double diffM=Max(0.0,Max(Max(difflM,diffrM),Max(diffuM,diffdM)));
        //double diffm=Min(0.0,Min(Min(difflm,diffrm),Min(diffum,diffdm)));
        //cout<<"here "<<diffM<<endl;

        double diffdlM=MaxValAve.get(i-1,j-1,me)-Q1;
        double diffdlm=MinValAve.get(i-1,j-1,me)-Q1;

        double diffdrM=MaxValAve.get(i+1,j-1,me)-Q1;
        double diffdrm=MinValAve.get(i+1,j-1,me)-Q1;

        double diffurM=MaxValAve.get(i+1,j+1,me)-Q1;
        double diffurm=MinValAve.get(i+1,j+1,me)-Q1;

        double diffulM=MaxValAve.get(i-1,j+1,me)-Q1;
        double diffulm=MinValAve.get(i-1,j+1,me)-Q1;
         

        diffM=Max(alpha,Max(Max(diffdlM,diffdrM),Max(diffurM,diffulM)));
        diffm=Min(-alpha,Min(Min(diffdlm,diffdrm),Min(diffurm,diffulm)));
        //if(me==1){diffm=Max(diffm,-Q1);}
        //cout<<"here2 "<<diffM<<" "<<diffdlM<<" "<<diffdrM<<" "<<diffurM<<" "<<diffulM<<endl;


        double diffcM=MaxVal.get(i,j,me)-Q1;
        double diffcm=MinVal.get(i,j,me)-Q1;
        //Compute the minimum theta value we need to bound
        //diffcM and diffcm between all of the neighbouring cell
        //differencesl
        double thetam1,thetaM1;
        if (fabs(diffcm)<1.0e-15)
        {
            thetam1=1.0;
        }
        else
        {
            thetam1=phi_func(diffm/diffcm);
        }
 

        if (fabs(diffcM)<1.0e-15)
        {
            thetaM1=1.0;
        }
        else
        {
            thetaM1=phi_func(diffM/diffcM);
        }

         

        double theta=Min(thetam1,thetaM1);
        theta=Min(1.0,theta);
        thetae=Min(theta,thetae);
        //if(1)//theta <1.0)
        //{cout<<"theta= "<<theta<<" "<<diffM<<" "<<diffm<<" "<<diffcM<<" "<<diffcm<<" "<<phi_func(diffm/diffcm)<<endl;}
        if(isnan(theta))
        {
         double rho=q.get(i,j,1,1);double u1=q.get(i,j,2,1)/q.get(i,j,1,1);double u2=q.get(i,j,3,1)/q.get(i,j,1,1);
         double eget=q.get(i,j,5,1);
         double Press=eget-0.5*rho*(u1*u1+u2*u2);
         cout<<MaxVal.get(i-1,j,me)<<endl;
         cout<<MaxVal.get(i+1,j,me)<<endl;
         cout<<MaxVal.get(i,j-1,me)<<endl;
         cout<<MaxVal.get(i,j+1,me)<<endl;
         cout<<i<<" "<<j<<" "<<me<<" "<<rho<<" "<<Press<<" here2 "<<diffM<<" "<<diffdlM<<" "<<diffdrM<<" "<<diffurM<<" "<<diffulM<<endl;}

        assert_ge(theta,0.0);
        assert_le(theta,1.0);
       /* for( int k=2;k<=kmax;k++)
        {
            q.set(i,j,me,k,q.get(i,j,me,k)*theta);
        } */ 
       }

       for(int k=2;k<=kmax;k++)
       {
          for(int me=1;me<=meqn;me++)
          {q.set(i,j,me,k,q.get(i,j,me,k)*thetae);}
       }
     }
    }

}


// Gaussian quadrature points.  However, an application will likely want to add in the edge
// points that are used for Riemann solves.
void SetPositivePointsGauss(const int& space_order, dTensor2& spts, dTensor1& w1d)
{

    dTensor1 s1d(space_order);

    // 1D Gaussian quadrature points
    switch ( space_order )
    {
        case 1:
            w1d.set(1, 2.0 );

            s1d.set(1, 0.0e0 );
            break;

        case 2:

            w1d.set(1,  1.0 );
            w1d.set(2,  1.0 );

            s1d.set(1, -1.0/sq3 );
            s1d.set(2,  1.0/sq3 );      
            break;

        case 3:

            w1d.set(1,  5.0/9.0 );
            w1d.set(2,  8.0/9.0 );
            w1d.set(3,  5.0/9.0 );

            s1d.set(1, -sq3/sq5 );
            s1d.set(2,  0.0e0 );
            s1d.set(3,  sq3/sq5 );
            break;

        case 4:

            w1d.set(1, (18.0 - sq3*sq10)/36.0 );
            w1d.set(2, (18.0 + sq3*sq10)/36.0 );
            w1d.set(3, w1d.get(2) );
            w1d.set(4, w1d.get(1) );


            s1d.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
            s1d.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );      
            break;

        case 5:

            w1d.set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
            w1d.set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
            w1d.set(3, 128.0/225.0 );
            w1d.set(4, w1d.get(2) );
            w1d.set(5, w1d.get(1) );

            s1d.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
            s1d.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            s1d.set(3,  0.0 );
            s1d.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            s1d.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );      
            break;
    }

    // This region has been commented out because we are no longer applying
    // the limiter at the corner points

    // 2D points -- corners
  spts.set(1,1, -1.0e0 );
  spts.set(1,2, -1.0e0 );

  spts.set(2,1,  1.0e0 );
  spts.set(2,2, -1.0e0 );

  spts.set(3,1, -1.0e0 );
  spts.set(3,2,  1.0e0 );

  spts.set(4,1,  1.0e0 );
  spts.set(4,2,  1.0e0 );

    // 2D points -- left, right, bottom and top edges
  for (int m=1; m<=space_order; m++)
  {
      double s = s1d.get(m);

      // left edge
      spts.set(4+m,1, -1.0e0 );
      spts.set(4+m,2,  s     );

      // right edge
      spts.set(4+space_order+m,1,  1.0e0 );
      spts.set(4+space_order+m,2,  s     );

      // bottom edge
      spts.set(4+2*space_order+m,1,  s     );
      spts.set(4+2*space_order+m,2, -1.0e0 );

      // top edge
      spts.set(4+3*space_order+m,1,  s     );
      spts.set(4+3*space_order+m,2,  1.0e0 );
  }

    // 2D points -- all interior points
    int z = 4+4*space_order;
    //z = 0;
    for (int m=1; m<=space_order; m++)
        for (int k=1; k<=space_order; k++)
        {
            double s1 = s1d.get(m);
            double s2 = s1d.get(k);

            z = z+1;

            spts.set(z,1, s1d.get(m) );
            spts.set(z,2, s1d.get(k) );
        }

}
