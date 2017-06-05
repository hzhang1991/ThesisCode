#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "mesh.h"                   // Required for all unstructured stuff
#include "Quadrature.h"             // 1D Quadrature rules
#include "MonomialsToLegendre.h"    // For evaluating basis functions
#include <iostream>
using namespace std;

// Routine that applies a modified Barth-Jespersen limiter to higher order moments of
// the conserved variables.  This limiter is locally applied, and conserves
// total mass by not adjusting cell averages.
//
// See: K. Michalak and C. Ollivier-Gooch, "Limiters for Unstructured Higher-Order
// Accurate Solutions of the Euler Equations"

/*
An example of one of the phi functions we have experimented with. This
phi function penalizes non-monotonic distributions. It is used in place
of what would typically be max(1,x). When we have shocks in several 
neighboring cells we need to use a phi function such that phi(1)<1 and
phi(x)=1 for x>alpha for 2 > alpha >1. The reason for this is
that if neighboring cells have shocks the bounds we would normally
use will be polluted and it is possible that x=1 even while the current
cell is experiencing nonphysical oscillations. However if we use one of these
phi functions over many steps we will penalize cells such that x=1. Eventually
this will eliminate the oscillations because a smooth 
monotonic function should actually satisfy x \ge 2.
*/


inline double phi_func(double x)
{
   return (x<1.5)*x*(1.0-4.0/27.0*x*x)+(x>=1.5)*1.0;
}


void ApplyLimiter_Unst(const mesh& Mesh, const dTensor3& aux, dTensor3& q)
{

    const int NumElems      = Mesh.get_NumElems();
    const int NumPhysElems  = Mesh.get_NumPhysElems();
    const int NumEdges      = Mesh.get_NumEdges();
    const int meqn          = q.getsize(2);
    const int kmax          = q.getsize(3);
    const int maux          = aux.getsize(2);
    const int space_order   = dogParams.get_space_order();

    // Do nothing in the case of piecewise constants
    if( space_order == 1 )
    { return; }

    // ------------------------------------------------ //
    // number of points where we want to check solution //
    // ------------------------------------------------ //
    const int space_order_sq = space_order*space_order;
    const int mpts_vec[] = {0, 3*space_order_sq, 18, 3*space_order_sq, 3*space_order_sq };  // TODO - FILL IN 2ND-ORDER CASE
    const int mpoints    = mpts_vec[space_order-1];

 
    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //
    dTensor2 spts(mpoints, 2);
    void SetPositivePointsBarJesp_Unst(const int& space_order, dTensor2& spts);
    SetPositivePointsBarJesp_Unst(space_order, spts);

    void SamplePhiAtPositivePointsBarJesp_Unst(const int& space_order, 
            const dTensor2& spts, dTensor2& phi);
    dTensor2 phi(mpoints, kmax);
    SamplePhiAtPositivePointsBarJesp_Unst(space_order, spts, phi);

    //--------------------------------------------------------------//
    // We will store the max and min values on each cell....//
    //--------------------------------------------------------------//
    dTensor2 MaxVal(NumElems,meqn);
    dTensor2 MinVal(NumElems,meqn);

//#pragma omp parallel for
    for(int  i=1;  i <= NumElems; i++)
    for(int me=1; me <= meqn; me++)
    {

        double minval = 1.0e16;
        double maxval = -1.0e16;
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double qnow = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                qnow += q.get(i,me,k) * phi.get(mp,k);
            }
            minval=Min(minval,qnow);
            maxval=Max(maxval,qnow);
           
        }

        MaxVal.set(i,me,maxval);
        MinVal.set(i,me,minval);
    }

//#pragma omp parallel for
    for(int  i=1;  i <= NumPhysElems; i++)
    {double thetae=1.0;
    for(int me=1; me <= meqn; me++)
    {
       int ileft  = Mesh.get_eelem(i,1);
       int iright = Mesh.get_eelem(i,2);
       double Q1=q.get(i,me,1);
       double diffM=0.0;
       double diffm=0.0;

       //Find the deviation from the max and min values
       //on neighbouring cells from the average value
       //on our current cell.
       for(int e=1;e<=3;e++)
       {
          int edge=Mesh.get_tedge(i,e);
          int index;
          int ileft  = Mesh.get_eelem(edge,1);
          int iright = Mesh.get_eelem(edge,2);
          if(ileft==i){index=iright;}
          else{index=ileft;}
          double diffM1=Max(diffM,MaxVal.get(index,me)-Q1);
          double diffm1=Min(diffm,MinVal.get(index,me)-Q1);
          if(isnan(diffM1) || isnan(diffm1)){cout<<NumPhysElems<<" problem here "<<index<<endl;exit(1);}
          diffM=diffM1;
       }
       double diffcM=MaxVal.get(i,me)-Q1;
       double diffcm=MinVal.get(i,me)-Q1;
        //Compute the minimum theta value we need to bound
        //diffcM and diffcm between all of the neighbouring cell
        //differencesl

       double thetam1,thetaM1;
       if (fabs(diffcm)<1.0e-14)
       {
           thetam1=1.0;
       }
       else
       {
           thetam1=phi_func(diffm/diffcm);
       }


       if (fabs(diffcM)<1.0e-14)
       {
           thetaM1=1.0;
       }
       else
       {
           thetaM1=phi_func(diffM/diffcM);
       }

       double theta=Min(thetam1,thetaM1);
       //cout<<" HERE "<<diffM<<" "<<diffm<<" "<<diffcM<<" "<<diffcm<<" "<<i<<endl;
       assert_ge(theta,0.0);
       assert_le(theta,1.0);
       //compute the single theta value to use for all entries in a system
       thetae=Min(theta,thetae);
       /*for( int k=2;k<=kmax;k++)
       {
            q.set(i,me,k,q.get(i,me,k)*theta);
       }*/
      }
      if(thetae<1.0)//phi_func(1.0))
      {//Limit all of the physical variables with thetae
        for(int me=1; me <= meqn; me++)
        { 
            for( int k=2; k <= kmax; k++ )
            {
                 q.set(i,me,k, q.get(i,me,k) * thetae );
            }  
        }

     }
    }

}


void SetPositivePointsBarJesp_Unst(const int& space_order, dTensor2& spts)
{

    // Positivity points with duplicates removed.  I am unaware of a clever
    // solution to construct these without saving duplicates
    //
    // TODO - hard code the points for second, fourth and fifth-order case
    //
    if( space_order == 3 )
    {
        spts.set(1,1, -2.206316679540750e-01 );
        spts.set(1,2, -3.333333333333334e-01 );
        spts.set(2,1, 5.539650012874083e-01  );
        spts.set(2,2, -2.206316679540750e-01 );
        spts.set(3,1, -3.333333333333334e-01 );
        spts.set(3,2, 5.539650012874083e-01  );
        spts.set(4,1, 1.666666666666667e-01  );
        spts.set(4,2, -3.333333333333333e-01 );
        spts.set(5,1, 1.666666666666667e-01  );
        spts.set(5,2, 1.666666666666667e-01  );
        spts.set(6,1, -3.333333333333333e-01 );
        spts.set(6,2, 1.666666666666667e-01  );
        spts.set(7,1, 5.539650012874083e-01  );
        spts.set(7,2, -3.333333333333334e-01 );
        spts.set(8,1, -2.206316679540750e-01 );
        spts.set(8,2, 5.539650012874083e-01  );
        spts.set(9,1, -3.333333333333334e-01 );
        spts.set(9,2, -2.206316679540750e-01);
        spts.set(10,1, -2.769825006437042e-01 );
        spts.set(10,2, -2.769825006437042e-01);
        spts.set(11,1, 5.539650012874084e-01 );
        spts.set(11,2, -2.769825006437042e-01);
        spts.set(12,1, -2.769825006437042e-01 );
        spts.set(12,2, 5.539650012874084e-01);
        spts.set(13,1, -8.333333333333334e-02 );
        spts.set(13,2, -8.333333333333333e-02);
        spts.set(14,1, 1.666666666666667e-01 );
        spts.set(14,2, -8.333333333333334e-02);
        spts.set(15,1, -8.333333333333333e-02 );
        spts.set(15,2, 1.666666666666667e-01);
        spts.set(16,1, 1.103158339770375e-01 );
        spts.set(16,2, 1.103158339770375e-01);
        spts.set(17,1, -2.206316679540750e-01 );
        spts.set(17,2, 1.103158339770375e-01);
        spts.set(18,1, 1.103158339770375e-01 );
        spts.set(18,2, -2.206316679540750e-01);
        return;

    } 

    // -------------------------------------------------------------- //
    // Positivity points for the general case.  There are many duplicates that
    // will be saved if you enter this part of the code!
    // -------------------------------------------------------------- //

    // Gaussian Quadrature points
    dTensor1 s1d_ga(space_order);
    setGaussPoints1d( s1d_ga );

    // Gauss lobatto points
    dTensor1 s1d_gl(space_order);  
    setGaussLobattoPoints1d( s1d_gl );

    // values for defining quadrature points (TODO - remove these)
    double A, B, C;

    // Vertices of Canonical triangle (in Clockwise order)
    double v1[] = {-1./3., -1./3.};
    double v2[] = {-1./3.,  2./3.};
    double v3[] = { 2./3., -1./3.};

// TESTING VERTICES
//  double v1[] = {-1., 0.};
//  double v2[] = { 1., 0.};
//  double v3[] = { 0., 1.};

    // 2D points (on the square) //
    dTensor3 uu( spts.getsize(1), spts.getsize(1), 2 );
    dTensor3 vv( spts.getsize(1), spts.getsize(1), 2 );
    for (int m=1; m<=space_order; m++)
    for (int k=1; k<=space_order; k++)
    {
        uu.set(m,k, 1, s1d_gl.get(m) );
        vv.set(m,k, 2, s1d_ga.get(k) );
    }

    // Weird mapping from unit square to triangles
    int z = 0;
    for( int m=1; m<=space_order; m++)
    for( int k=1; k<=space_order; k++)
    {

        double ua = uu.get(m,k,1);
        double vb = vv.get(m,k,2);

        A = 0.5*( 1.0 + vb );
        B = 0.25*( 1. + ua )*( 1. - vb );
        C = 0.25*( 1. - ua )*( 1. - vb );

        z++;
        spts.set( z, 1, A*v1[0] + B*v2[0] + C*v3[0] );
        spts.set( z, 2, A*v1[1] + B*v2[1] + C*v3[1] );
//      printf("spts(%d,:) = [%2.15e %2.15e];\n", z, spts.get(z,1), spts.get(z,2) );

        z++;
        spts.set( z, 1, A*v3[0] + B*v1[0] + C*v2[0] );
        spts.set( z, 2, A*v3[1] + B*v1[1] + C*v2[1] );
//      printf("spts(%d,:) = [%2.15e %2.15e];\n", z, spts.get(z,1), spts.get(z,2) );

        z++;
        spts.set( z, 1, A*v2[0] + B*v3[0] + C*v1[0] );
        spts.set( z, 2, A*v2[1] + B*v3[1] + C*v1[1] );
//      printf("spts(%d,:) = [%2.15e %2.15e];\n", z, spts.get(z,1), spts.get(z,2) );

    }

}


// Construct vector phi, such that Q \cdot phi = qvals( 1:numpts ).  That is,
// given a Galerkin representation:
//
//      q(x,y) = \sum_{k=1}^kmax Q^{(k)}_i \varphi^{(k)}( \xi, \eta ),
//
// we want phi such that
//
//   qvals( 1:numpts, 1:meqn ) = \sum_k  Q( i, 1:meqn, k ) * phi(1:numpts, k )
//
// TODO - OR SOME VARIATION OF THIS!
void SamplePhiAtPositivePointsBarJesp_Unst(const int& space_order, 
        const dTensor2& spts, dTensor2& phi)
{

    const int mpoints = phi.getsize(1);  assert_eq( mpoints, spts.getsize(1) );
    const int kmax    = phi.getsize(2);

    // Loop over each quadrature point and construct monomial polys
    //
    // TODO - consolidate this with L2Project_Unst.cpp
    //
    dTensor2   mu(mpoints, kmax); // monomial basis (non-orthogonal)
    for (int m=1; m<=mpoints; m++)
    {
        // coordinates
        const double xi   = spts.get(m,1);      
        const double xi2  = xi*xi;
        const double xi3  = xi2*xi;
        const double xi4  = xi3*xi;
        const double eta  = spts.get(m,2);
        const double eta2 = eta*eta;
        const double eta3 = eta2*eta;
        const double eta4 = eta3*eta;      

        // monomials basis (non-orthogonal)
        switch( kmax )
        {
            case 15:  // fifth order		    		    
                mu.set(m, 15, eta4     );
                mu.set(m, 14, xi4      );
                mu.set(m, 13, xi2*eta2 );
                mu.set(m, 12, eta3*xi  );
                mu.set(m, 11, xi3*eta  );

            case 10:  // fourth order
                mu.set(m, 10, eta3     );
                mu.set(m, 9,  xi3      );
                mu.set(m, 8,  xi*eta2  );
                mu.set(m, 7,  eta*xi2  );

            case 6:  // third order
                mu.set(m, 6,  eta2     );
                mu.set(m, 5,  xi2      );
                mu.set(m, 4,  xi*eta   );		    

            case 3:  // second order		    
                mu.set(m, 3, eta       );
                mu.set(m, 2, xi        );

            case 1:  // first order
                mu.set(m, 1, 1.0       );

                break;		    
        }
    }

    // Construct pointwise values of the Legendre polynomials
    for (int m=1; m<=mpoints; m++)    
    for (int k=1; k<=kmax; k++)
    {
        double tmp = 0.0;
        for (int j=1; j<=k; j++)
        {  tmp = tmp + Mmat[k-1][j-1]*mu.get(m,j);  }
        phi.set(m,k, tmp );      
    }

}
