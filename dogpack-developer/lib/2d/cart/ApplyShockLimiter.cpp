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
//
// TODO - see our fantastic arxived paper on this work too! -DS (3/3/2015)

/*  Phi function used for limiting.
 * 
 *  An example of one of the phi functions we have experimented with. 
 *  This phi function penalizes non-monotonic distributions, which is important for
 *  limiting solutions where the neighboring value indicates oscillations. 
 *  It is used in place of what would typically be max(1,x). When we have shocks in several 
 *  neighboring cells we need to use a phi function such that phi(1)<1 and
 *  phi(x)=1 for x>alpha for 2 > alpha >1. 
 *
 *  The reason for this is that if neighboring cells have shocks the bounds we would normally
 *  use will be polluted and it is possible that x=1 even while the current
 *  cell is experiencing nonphysical oscillations. However if we use one of these
 *  phi functions over many steps we will penalize cells such that x=1. Eventually
 *  this will eliminate the oscillations because a smooth 
 *  monotonic function should actually satisfy x \ge 2.
 *
 */
inline double phi_func(double x)
{
    //return Min(1.0,x);    // This method does not work well in 2D!  (As described earlier)
    return Min(x/1.1,1.0);  // many other options exist!  For example, (x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
}

void ApplyLimiter(dTensorBC4& aux, dTensorBC4& q,
    void (*ProjectRightEig)(int, const dTensor1&, const dTensor1&, const dTensor2&, dTensor2&),
    void (*ProjectLeftEig)(int, const dTensor1&, const dTensor1&, const dTensor2&, dTensor2&))
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
    const int mpoints = 4+4*space_order+space_order*space_order;

    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //

    // Use Gauss Quadrature points augmented with edge values
    dTensor2 spts(mpoints,2);
    dTensor1 w1d(space_order);
    void SetPositivePointsGauss(const int& space_order, dTensor2& spts,dTensor1& w1d);
    SetPositivePointsGauss(space_order,spts,w1d);

    // Basis functions evaluated at each point
    dTensor2 phi(mpoints,kmax);
    SetLegendreAtPoints(space_order,spts,phi);
  
    // A characteristic length...something more complicated is required on unstructured grds 
    double dx=1.0/mx;

    //A matrix to store max and min values on each cell. 
    //We need this to see if a cell's point values fall
    //within the range of the values of its neighbors. 
    //This will determine bounds on each cell.

    dTensorBC3 MaxVal(mx,my,meqn,mbc); MaxVal.setall(-1.0e18);
    dTensorBC3 MinVal(mx,my,meqn,mbc); MinVal.setall(1.0e18);
    
    //Tensor to hold edge differences in the DG solution on different
    //Sides of faces. This is a shock detecting method introduced in the 
    //"Shock detection and limiting with discontinuous Galerkin methods
    //for hyperbolic conservation laws"-Krivodonova
    dTensorBC4 edgeval(mx+1,my+1,meqn,4,mbc);edgeval.setall(0.0e0);


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
        }
        
        //Compute Extrema Values
        MaxVal.set(i,j,me,max1);
        MinVal.set(i,j,me,min1);
    }

    // Add a tolerance to avoid smearing extrema
    double alpha = 1.0*pow(dx,1.5);  // TODO - do we want to include something involving dy??? -DS

    // Use this to determine when to turn on the the limiter 
    // double cutoff=pow(dx,2.5);
    // use the shock detector. Basically checking if any jumps accross edges are over the limit for a given cell.
    // If they are set detector=1
    //
#pragma omp parallel for
    for(int i=1; i <= mx; i++)
    for(int j=1; j <= my; j++)
    { 

        // Compute a single value of theta for this current element
        double thetae = 1.0;
        for(int me=1; me <= meqn; me++)
        {
            double Q1=q.get(i,j,me,1);

            // Find the (ratio of the) deviation from the max and min values
            // on neighbouring cells from the average value
            // on our current cell.

            double difflM=MaxVal.get(i-1,j,me)-Q1;
            double difflm=MinVal.get(i-1,j,me)-Q1;

            double diffrM=MaxVal.get(i+1,j,me)-Q1;
            double diffrm=MinVal.get(i+1,j,me)-Q1;

            double diffuM=MaxVal.get(i,j+1,me)-Q1;
            double diffum=MinVal.get(i,j+1,me)-Q1;

            double diffdM=MaxVal.get(i,j-1,me)-Q1;
            double diffdm=MinVal.get(i,j-1,me)-Q1; 

            //Compare the max and min computed with a tolerance alpha,
            //this is designed to prevent smearing extrema

            double diffM = Max(alpha,Max(Max(difflM,diffrM),Max(diffuM,diffdM)));
            double diffm = Min(-alpha,Min(Min(difflm,diffrm),Min(diffum,diffdm)));

            double diffcM=MaxVal.get(i,j,me)-Q1;
            double diffcm=MinVal.get(i,j,me)-Q1;

            // Compute the minimum theta value we need to bound
            // diffcM and diffcm between all of the neighbouring cell
            // differencesl
            double thetam1,thetaM1;
            if (fabs(diffcm)<1.0e-15)
            { thetam1=1.0; }
            else
            { thetam1=phi_func(diffm/diffcm); }

            if (fabs(diffcM)<1.0e-15)
            { thetaM1=1.0; }
            else
            { thetaM1=phi_func(diffM/diffcM); }

            double theta=Min(thetam1,thetaM1);

            // Compute a theta for our whole system. This is done because shocks 
            // are in characteristic variables not the physical variables.
            theta  = Min(1.0,theta);
            thetae = Min(theta,thetae);
            if(isnan(theta)) { cout<<"a nan theta has been encountered"<<endl;}

            // Last second chance to make sure we don't introduce extra
            // oscillations
            assert_ge(theta,0.0);
            assert_le(theta,1.0);

        }

        // Limit every physical variable using the computed value of theta
        for(int k=2;k<=kmax;k++)
        for(int me=1;me<=meqn;me++) {q.set(i,j,me,k,q.get(i,j,me,k)*thetae);}

    }

}

// Points that are used for Riemann solves.
//
// TODO - we should be able to use something from main library instead of this
// routine! -DS
//
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
