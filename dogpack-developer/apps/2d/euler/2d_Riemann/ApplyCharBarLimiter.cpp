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
   return (x<1.1)*Min(1.0/1.1*x,1.0)+(x>=1.1)*1.0;//(x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
}

void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q);


void ProjectLeftEig(int ixy,
                    const dTensor1& Aux_ave,
                    const dTensor1& Q_ave,
                    const dTensor2& Qvals,
                    dTensor2& Wvals);

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

    // Use Gauss Quadrature points augmented with edge values
    dTensor2 spts(mpoints,2);
    void SetPositivePointsGauss(const int& space_order, dTensor2& spts);
    SetPositivePointsGauss(space_order,spts);

    // Basis functions evaluated at each point
    dTensor2 phi(mpoints,kmax);
    SetLegendreAtPoints(space_order,spts,phi);


    //A matrix to store max and min values on each cell. 
    //We need this to see if a cell's point values fall
    //within the range of the values of its neighbors. 
    //This will determine bounds on each cell.

    dTensorBC3 MaxVal(mx,my,meqn,mbc); MaxVal.setall(-1.0e18);
    dTensorBC3 MinVal(mx,my,meqn,mbc); MinVal.setall(1.0e18);

    dTensorBC3 MaxValAve(mx,my,meqn,mbc); MaxValAve.setall(-1.0e18);
    dTensorBC3 MinValAve(mx,my,meqn,mbc); MinValAve.setall(1.0e18);



    //Store the max values on each cell
    dTensorBC3 maxcell (mx,my, meqn,mbc);
    dTensorBC3 mincell (mx,my, meqn,mbc);
    dTensorBC3 avecell (mx,my, meqn,mbc);

    dTensorBC3 maxW1c (mx,my, meqn,mbc);
    dTensorBC3 maxW1s (mx,my, meqn,mbc);
    dTensorBC3 maxW1e (mx,my, meqn,mbc);
    dTensorBC3 maxW1n (mx,my, meqn,mbc);
    dTensorBC3 maxW1w (mx,my, meqn,mbc);

    dTensorBC3 minW1c (mx,my, meqn,mbc);
    dTensorBC3 minW1s (mx,my, meqn,mbc);
    dTensorBC3 minW1e (mx,my, meqn,mbc);
    dTensorBC3 minW1n (mx,my, meqn,mbc);
    dTensorBC3 minW1w (mx,my, meqn,mbc);

    dTensorBC3 aveW1 (mx,my, meqn,mbc);

    dTensorBC3 maxW2c (mx,my, meqn,mbc);
    dTensorBC3 maxW2s (mx,my, meqn,mbc);
    dTensorBC3 maxW2e (mx,my, meqn,mbc);
    dTensorBC3 maxW2n (mx,my, meqn,mbc);
    dTensorBC3 maxW2w (mx,my, meqn,mbc);

    dTensorBC3 minW2c (mx,my, meqn,mbc);
    dTensorBC3 minW2s (mx,my, meqn,mbc);
    dTensorBC3 minW2e (mx,my, meqn,mbc);
    dTensorBC3 minW2n (mx,my, meqn,mbc);
    dTensorBC3 minW2w (mx,my, meqn,mbc);

    dTensorBC3 aveW2 (mx,my, meqn,mbc);


    dTensor2 WScell1(meqn,kmax);WScell1.setall(0.0) ;
    dTensor2 WScell2(meqn,kmax);WScell2.setall(0.0) ;

    dTensorBC4 Wcell1(mx,my,meqn,kmax,mbc);Wcell1.setall(0.0) ;
    dTensorBC4 Wcell2(mx,my,meqn,kmax,mbc);Wcell2.setall(0.0) ;
    dTensor1 Auxavs(maux);
    dTensor1 qavs(meqn);
    dTensor2 qsingle(meqn,kmax);
    dTensor2 Wsingle(meqn,kmax);
    // -------------------------------------------------------------- //
    //                                                                //
    // Loop over each cell and every quadrature point on each cell    //
    // I do not see a clever way of avoiding this right now 
    // unless we want to compute approximate max and min...           //
    // -------------------------------------------------------------- //
#pragma omp parallel for
    for(int i=1-mbc; i <= mx+mbc; i++)
    for(int j=1-mbc; j <= my+mbc; j++)
    {
     for(int me=1; me <= meqn; me++)
     {
        qavs.set(me,q.get(i,j,me,1));
        Auxavs.set(me,aux.get(i,j,me,1));
        double max1 = MaxVal.get(i,j,me);
        double min1 = MinVal.get(i,j,me);
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double qnow = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                qnow += q.get(i,j,me,k) * phi.get(mp,k);
                qsingle.set(me,k,q.get(i,j,me,k));
            }
            max1=Max(max1,qnow);
            min1=Min(min1,qnow);
        }

        MaxVal.set(i,j,me,max1);
        MinVal.set(i,j,me,min1);
        //Compute extrema values
        MaxValAve.set(i,j,me,q.get(i,j,me,1));
        MinValAve.set(i,j,me,q.get(i,j,me,1));


     }
     Wsingle.setall(0.0); 
     ProjectLeftEig(1,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              Wcell1.set(i,j,me,k,Wsingle.get(me,k));
            }
     }

     Wsingle.setall(0.0);
     ProjectLeftEig(2,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              Wcell2.set(i,j,me,k,Wsingle.get(me,k));
            }
     }


    }


    for(int i=1; i <= mx; i++)
    for(int j=1; j <= my; j++)
    {

     for(int me=1; me <= meqn; me++)
     {
        qavs.set(me,q.get(i,j,me,1));
        Auxavs.set(me,aux.get(i,j,me,1));
     }
     /////////////////////////////////////////////////////////////////
     for(int me=1; me <= meqn; me++)
     {
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            for( int k=1; k <= kmax; k++ )
            {
                qsingle.set(me,k,q.get(i,j,me,k));
            }
        }

     }
     Wsingle.setall(0.0);
     ProjectLeftEig(1,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            { 
              WScell1.set(me,k,Wsingle.get(me,k));
            }
     }

     Wsingle.setall(0.0);
     ProjectLeftEig(2,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            { 
              WScell2.set(me,k,Wsingle.get(me,k));
            }
     }

     for(int me=1; me <= meqn; me++)
     {

        double max1 = -1.0e6;
        double min1 = 1.0e6;
        double max2 = -1.0e6;
        double min2 = 1.0e6;
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double w1now = 0.0; 
            double w2now = 0.0; 
            for( int k=1; k <= kmax; k++ )
            {   
                w1now += WScell1.get(me,k) * phi.get(mp,k);
                w2now += WScell2.get(me,k) * phi.get(mp,k);
            }
            max1=Max(max1,w1now);
            min1=Min(min1,w1now);
            max2=Max(max2,w2now);
            min2=Min(min2,w2now);
        }
      
        maxW1c.set(i,j,me,max1);
        minW1c.set(i,j,me,min1);
        aveW1.set(i,j,me,WScell1.get(me,1));
        maxW2c.set(i,j,me,max2);
        minW2c.set(i,j,me,min2);
        aveW2.set(i,j,me,WScell2.get(me,1));

      }
     /////////////////////////////////////////////////////
     for(int me=1; me <= meqn; me++)
     {
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            for( int k=1; k <= kmax; k++ )
            {
                qsingle.set(me,k,q.get(i,j-1,me,k));
            }
        }

     }
     Wsingle.setall(0.0);
     ProjectLeftEig(1,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              WScell1.set(me,k,Wsingle.get(me,k));
            }
     }

     Wsingle.setall(0.0);
     ProjectLeftEig(2,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              WScell2.set(me,k,Wsingle.get(me,k));
            }
     }




     for(int me=1; me <= meqn; me++)
     {

        /////////////////////////////////////////////////
        double max1 = -1.0e6;
        double min1 = 1.0e6;
        double max2 = -1.0e6;
        double min2 = 1.0e6;
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double w1now = 0.0;
            double w2now = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                w1now += WScell1.get(me,k) * phi.get(mp,k);
                w2now += WScell2.get(me,k) * phi.get(mp,k);
            }
            max1=Max(max1,w1now);
            min1=Min(min1,w1now);
            max2=Max(max2,w2now);
            min2=Min(min2,w2now);
        }

        maxW1s.set(i,j,me,max1);
        minW1s.set(i,j,me,min1);
        maxW2s.set(i,j,me,max2);
        minW2s.set(i,j,me,min2);
     }
        /////////////////////////////////////////////////


     for(int me=1; me <= meqn; me++)
     {
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            for( int k=1; k <= kmax; k++ )
            {
                qsingle.set(me,k,q.get(i+1,j,me,k));
            }
        }

     }
     Wsingle.setall(0.0);
     ProjectLeftEig(1,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              WScell1.set(me,k,Wsingle.get(me,k));
            }
     }

     Wsingle.setall(0.0);
     ProjectLeftEig(2,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              WScell2.set(me,k,Wsingle.get(me,k));
            }
     }

     for(int me=1; me <= meqn; me++)
     {

        double max1 = -1.0e6;
        double min1 = 1.0e6; 
        double max2 = -1.0e6;
        double min2 = 1.0e6; 
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double w1now = 0.0;
            double w2now = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                w1now += WScell1.get(me,k) * phi.get(mp,k);
                w2now += WScell2.get(me,k) * phi.get(mp,k);
            }
            max1=Max(max1,w1now);
            min1=Min(min1,w1now);
            max2=Max(max2,w2now);
            min2=Min(min2,w2now);
        }

        maxW1e.set(i,j,me,max1);
        minW1e.set(i,j,me,min1);
        maxW2e.set(i,j,me,max2);
        minW2e.set(i,j,me,min2);
      }
        /////////////////////////////////////////////////
     for(int me=1; me <= meqn; me++)
     {
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            for( int k=1; k <= kmax; k++ )
            {
                qsingle.set(me,k,q.get(i,j+1,me,k));
            }
        }

     }
     Wsingle.setall(0.0);
     ProjectLeftEig(1,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              WScell1.set(me,k,Wsingle.get(me,k));
            }
     }

     Wsingle.setall(0.0);
     ProjectLeftEig(2,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              WScell2.set(me,k,Wsingle.get(me,k));
            }
     }

     for(int me=1; me <= meqn; me++)
     {
        double max1 = -1.0e6;
        double min1 = 1.0e6; 
        double max2 = -1.0e6;
        double min2 = 1.0e6; 
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double w1now = 0.0;
            double w2now = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                w1now += WScell1.get(me,k) * phi.get(mp,k);
                w2now += WScell2.get(me,k) * phi.get(mp,k);
            }
            max1=Max(max1,w1now);
            min1=Min(min1,w1now);
            max2=Max(max2,w2now);
            min2=Min(min2,w2now);
        }

        maxW1n.set(i,j,me,max1);
        minW1n.set(i,j,me,min1);
        maxW2n.set(i,j,me,max2);
        minW2n.set(i,j,me,min2);
     }
        /////////////////////////////////////////////////
     for(int me=1; me <= meqn; me++)
     {
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            for( int k=1; k <= kmax; k++ )
            {
                qsingle.set(me,k,q.get(i-1,j,me,k));
            }
        }

     }
     Wsingle.setall(0.0);
     ProjectLeftEig(1,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              WScell1.set(me,k,Wsingle.get(me,k));
            }
     }

     Wsingle.setall(0.0);
     ProjectLeftEig(2,Auxavs,qavs,qsingle,Wsingle);
     for(int me=1; me <= meqn; me++)
     {
            for( int k=1; k <= kmax; k++ )
            {
              WScell2.set(me,k,Wsingle.get(me,k));
            }
     }
     for(int me=1; me <= meqn; me++)
     {

        double max1 = -1.0e6;
        double min1 = 1.0e6; 
        double max2 = -1.0e6;
        double min2 = 1.0e6; 
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double w1now = 0.0;
            double w2now = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                w1now += WScell1.get(me,k) * phi.get(mp,k);
                w2now += WScell2.get(me,k) * phi.get(mp,k);
            }
            max1=Max(max1,w1now);
            min1=Min(min1,w1now);
            max2=Max(max2,w2now);
            min2=Min(min2,w2now);
        }

        maxW1w.set(i,j,me,max1);
        minW1w.set(i,j,me,min1);
        maxW2w.set(i,j,me,max2);
        minW2w.set(i,j,me,min2);
     }
        /////////////////////////////////////////////////////////////////////


    }
    //Add a tolerance.
    double dx=1.0/mx;
    double alpha=1000.0*pow(dx,1.5);
#pragma omp parallel for
    for(int i=1; i <= mx; i++)
    for(int j=1; j <= my; j++)
    { double thetae=1.0;
    for(int me=1; me <= meqn; me++)
    {

       //Find the deviation from the max and min values
       //on neighbouring cells from the average value
       //on our current cell.


        ////////////////////////////////////////////////W1 theta
        double Q1=aveW1.get(i,j,me);
        double difflM1=maxW1w.get(i,j,me)-Q1;
        double difflm1=minW1w.get(i,j,me)-Q1;

        double diffrM1=maxW1e.get(i,j,me)-Q1;
        double diffrm1=minW1e.get(i,j,me)-Q1;

        double diffuM1=maxW1n.get(i,j,me)-Q1;
        double diffum1=minW1n.get(i,j,me)-Q1;

        double diffdM1=maxW1s.get(i,j,me)-Q1;
        double diffdm1=minW1s.get(i,j,me)-Q1; 
    
        //Compute with a tolerance

        //double diffM1=Max(alpha,Max(Max(difflM1,diffrM1),Max(diffuM1,diffdM1)));
        //double diffm1=Min(-alpha,Min(Min(difflm1,diffrm1),Min(diffum1,diffdm1)));
        double diffM1=Max(alpha,Max(difflM1,diffrM1));
        double diffm1=Min(-alpha,Min(difflm1,diffrm1));
 
        double diffcM1=maxW1c.get(i,j,me)-Q1;
        double diffcm1=minW1c.get(i,j,me)-Q1;

        double thetam1,thetaM1;
        if (fabs(diffcm1)<1.0e-15)
        {
            thetam1=1.0;
        }
        else
        {
            thetam1=phi_func(diffm1/diffcm1);
        }
 

        if (fabs(diffcM1)<1.0e-15)
        {
            thetaM1=1.0;
        }
        else
        {
            thetaM1=phi_func(diffM1/diffcM1);
        }

         

        double theta1=Min(thetam1,thetaM1);

        ////////////////////////////////////////////////W2 theta
        double Q2=aveW2.get(i,j,me);
        double difflM2=maxW2w.get(i,j,me)-Q2;
        double difflm2=minW2w.get(i,j,me)-Q2;

        double diffrM2=maxW2e.get(i,j,me)-Q2;
        double diffrm2=minW2e.get(i,j,me)-Q2;

        double diffuM2=maxW2n.get(i,j,me)-Q2;
        double diffum2=minW2n.get(i,j,me)-Q2;

        double diffdM2=maxW2s.get(i,j,me)-Q2;
        double diffdm2=minW2s.get(i,j,me)-Q2;

        //Compute with a tolerance

        //double diffM2=Max(alpha,Max(Max(difflM2,diffrM2),Max(diffuM2,diffdM2)));
        //double diffm2=Min(-alpha,Min(Min(difflm2,diffrm2),Min(diffum2,diffdm2)));
        double diffM2=Max(alpha,Max(diffuM2,diffdM2));
        double diffm2=Min(-alpha,Min(diffum2,diffdm2));
 
        double diffcM2=maxW2c.get(i,j,me)-Q2;
        double diffcm2=minW2c.get(i,j,me)-Q2;

        double thetam2,thetaM2;
        if (fabs(diffcm2)<1.0e-15)
        {
            thetam2=1.0;
        }
        else
        {
            thetam2=phi_func(diffm2/diffcm2);
        }


        if (fabs(diffcM2)<1.0e-15)
        {
            thetaM2=1.0;
        }
        else
        {
            thetaM2=phi_func(diffM2/diffcM2);
        }



        double theta2=Min(thetam2,thetaM2);

 //////////////////////////////////////////////////////////////////////////////////
        //Compute a theta for our whole system. This is done because shocks are in characteristic variables not the physical variables.
        double theta11=Max(0.0,Min(1.0,theta1));
        double theta21=Max(0.0,Min(1.0,theta2));
 
        double theta=Max(theta11,theta21);
 
        //if(theta<1.0)
        //{printf("theta=%f\n",theta);}

        ////////////////////////////////
        thetae=Min(theta,thetae);
        }
       for(int k=2;k<=kmax;k++)
       {
          for(int me=1;me<=meqn;me++)
          {q.set(i,j,me,k,q.get(i,j,me,k)*thetae);}
       }
    }
    ApplyPosLimiter(aux, q);

}


// Gaussian quadrature points.  However, an application will likely want to add in the edge
// points that are used for Riemann solves.
void SetPositivePointsGauss(const int& space_order, dTensor2& spts)
{

    dTensor1 s1d(space_order);

    // 1D Gaussian quadrature points
    switch ( space_order )
    {
        case 1:
            s1d.set(1, 0.0e0 );
            break;

        case 2:
            s1d.set(1, -1.0/sq3 );
            s1d.set(2,  1.0/sq3 );      
            break;

        case 3:
            s1d.set(1, -sq3/sq5 );
            s1d.set(2,  0.0e0 );
            s1d.set(3,  sq3/sq5 );
            break;

        case 4:
            s1d.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
            s1d.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );      
            break;

        case 5:
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
