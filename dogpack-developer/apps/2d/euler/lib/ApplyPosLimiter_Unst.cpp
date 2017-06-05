#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "mesh.h"                   // Required for all unstructured stuff
#include "Quadrature.h"             // 1D Quadrature rules
#include "MonomialsToLegendre.h"    // For evaluating basis functions
#include<iostream>
using namespace std;

// This is a positivity preserving limiter similar to that proposed in 
// "Maximum-Principle-Satisfying and Positivity-Preserving
// High Order Discontinuous Galerkin Schemes
// for Conservation Laws on Triangular Meshes", Zhang, Xia and Shu
// J. Sci. Comput. (2012).
//
// THIS METHOD IS WRITTEN FOR THE EULER EQUATIONS SPECIFICALLY
//
// In order to implement this for a different scheme, one should rewrite, or
// redefine what components should remain positiive.  This will require
// reworking the control flow logic for how time step lengths are chosen.
void ApplyPosLimiter_Unst(const mesh& Mesh, const dTensor3& aux, dTensor3& q)
{

    const int NumElems      = Mesh.get_NumElems();
    const int NumPhysElems  = Mesh.get_NumPhysElems();
    const int NumEdges      = Mesh.get_NumEdges();
    const int meqn          = q.getsize(2);
    const int kmax          = q.getsize(3);
    const int maux          = aux.getsize(2);
    const int space_order   = dogParams.get_space_order();
    //eps will be treated as zero
    double eps=1.0e-10;

    // Do nothing in the case of piecewise constants
    if( space_order == 1 )
    { return; }

    // ------------------------------------------------ //
    // number of points where we want to check solution //
    // ------------------------------------------------ //
    const int space_order_sq = space_order*space_order;
    const int mpts_vec[] = {0, 3*space_order_sq, 18/*18*/, 3*space_order_sq, 3*space_order_sq };  // TODO - FILL IN 2ND-ORDER CASE
    const int mpoints    = mpts_vec[space_order-1];

    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //
    dTensor2 spts(mpoints, 2);
    void SetPositivePoints_Unst(const int& space_order, dTensor2& spts);
    SetPositivePoints_Unst(space_order, spts);

    void SamplePhiAtPositivePoints_Unst(const int& space_order, 
            const dTensor2& spts, dTensor2& phi);
    dTensor2 phi(mpoints, kmax);
    SamplePhiAtPositivePoints_Unst(space_order, spts, phi);

    // -------------------------------------------------------------- //
    // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
    // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
    // -------------------------------------------------------------- //
    #pragma omp parallel for
    for(int  i=1;  i <= NumPhysElems; i++)
    {
                 // ---------------------------------- //
                 // Limit the Density and Energy
                 // ---------------------------------- //
                 double mine = eps;
                 double minrho = eps;
                 double thetae = 1.0;
                 double thetarho = 1.0;
                 if(q.get(i,1,1)<0.0 ){printf(" Density average is negative ");}
                 //Check our set of points
                 for(int mp=1; mp <= mpoints; mp++)
                 {
                    // evaluate e and rho at spts(mp)
                    double enow = 0.0;
                    double rhonow = 0.0;

                    for( int k=1; k <= kmax; k++ )
                    {
                       enow += q.get(i,5,k) * phi.get(mp,k);
                       rhonow += q.get(i,1,k) * phi.get(mp,k);
                    }
                    mine=Min(mine,enow);
                    minrho=Min(minrho,rhonow);
                 }
                 //Shrink the higher moments if necessary for positivity at any of the checked points
                 double Q1 = q.get(i,1,1);
                 if( fabs( Q1 - minrho ) < eps|| Q1 < eps ){ thetarho = 0.0; }
                 else if(minrho<0.0){ thetarho = Min( thetarho, fabs( (eps-Q1) / (minrho-Q1) ) ); }

                 Q1 = q.get(i,5,1);
                 if( fabs( Q1 - mine ) < eps|| Q1 < eps ){ thetae = 0.0; }
                 else if (mine<0.0){ thetae = Min( thetae, fabs( (eps-Q1) / (mine-Q1) ) ); }
                 //---------------------------------------------------
                 // limit e and rho //
                 // Shrink by multiplying by theta 0 < theta <1, also energy 
                 // doesn't need to be adjusted here seperately.
                 //---------------------------------------------------
                  for( int k=2; k <= kmax; k++ )
                 {
                    q.set(i,1,k, q.get(i,1,k) * thetarho );
                    //q.set(i,5,k, q.get(i,5,k) * thetae );
                 }

                 // ------------------------------------------------ //
                 // Limit the Momentum (to obtain positive pressure)
                 //
                 // This step works by considering an expansion of density, momentum
                 // and energy in terms of a single parameter, theta.
                 //
                 // ------------------------------------------------ //

                 double m = eps;
                 double ave=q.get(i,1,1)*q.get(i,5,1)-0.5*(q.get(i,2,1)*q.get(i,2,1)+q.get(i,3,1)*q.get(i,3,1)+q.get(i,4,1)*q.get(i,4,1));
                 double thetam = 1.0;

                 if(ave<0.0 || isnan(ave)){cout<<q.get(i,1,1)<<" "<<q.get(i,2,1)<<" "<<q.get(i,3,1)<<" "<<q.get(i,4,1)<<" "<<q.get(i,5,1)<<" Avg Pressure PROBLEM! "<<ave<<endl;exit(1);}


                 double mxa = q.get(i,2,1);
                 double mya = q.get(i,3,1);
                 double mza = q.get(i,4,1);
                 double ea = q.get(i,5,1);
                 double rhoa = q.get(i,1,1);
                 for(int mp=1; mp <= mpoints; mp++)
                 {
                    // evaluate q at spts(mp) //
                    double mxnow = 0.0;
                    double mynow = 0.0;
                    double mznow = 0.0;
                    double enow = 0.0;
                    double rhonow = 0.0;


                    for( int k=1; k <= kmax; k++ )
                    {
                       mxnow += q.get(i,2,k) * phi.get(mp,k);
                       mynow += q.get(i,3,k) * phi.get(mp,k);
                       mznow += q.get(i,4,k) * phi.get(mp,k);
                       enow += q.get(i,5,k) * phi.get(mp,k);
                       rhonow += q.get(i,1,k) * phi.get(mp,k);
                    }

                    double mxc = mxnow-mxa;
                    double myc = mynow-mya;
                    double mzc = mznow-mza;
                    double ec = enow-ea;
                    double rhoc = rhonow-rhoa;

                    double slope=ec*rhoa+rhoc*ea-(mxc*mxa+myc*mya+mzc*mza);
                    double press=rhonow*enow-0.5*(mxnow*mxnow+mynow*mynow+mznow*mznow);
                 // Here we explicitly compute the slope of the quadratic pressure
                 // (in terms of theta). We then Compute the intersection of the tangent 
                 // line with zero.
                 if( press<0.0 && fabs( slope ) < eps ){ thetam = thetam; }
                 else if (press<0.0){ thetam = Min( thetam, fabs( (eps-ave) / slope ) ); }

                    m = Min(m, press);
                 }
                //cout<<ave<<endl;
                //Here we compare the minimum tangent line zero intersection with the zero 
                // intersection of the secant from the most negative point and the point theta=0.
                double theta = 1.0;
                Q1 = ave;
                if( fabs( Q1 - m ) < eps || Q1 < eps){ theta = 0.0; }
                else if( m < eps){ theta = Min( thetam, fabs( (eps-Q1) / (m-Q1) ) ); }

                // limit all solution variables //
                for (int me=1;me<=meqn;me++)
                {for( int k=2; k <= kmax; k++ )
                {
                    q.set(i,me,k, q.get(i,me,k) * theta );
                }}
                


                 /*for(int mp=1; mp <= mpoints; mp++)
                 {
                    // evaluate q at spts(mp) //
                    double mxnow = 0.0;
                    double mynow = 0.0;
                    double mznow = 0.0;
                    double enow = 0.0;
                    double rhonow = 0.0;


                    for( int k=1; k <= kmax; k++ )
                    {
                       mxnow += q.get(i,2,k) * phi.get(mp,k);
                       mynow += q.get(i,3,k) * phi.get(mp,k);
                       mznow += q.get(i,4,k) * phi.get(mp,k);
                       enow += q.get(i,5,k) * phi.get(mp,k);
                       rhonow += q.get(i,1,k) * phi.get(mp,k);
                    }

                    double press=rhonow*enow-0.5*(mxnow*mxnow+mynow*mynow+mznow*mznow);

                 if( press<0.0){cout<<theta<<" pressure problem "<<press<<" "<<m<<" "<< Q1<<endl;exit(1);}
                 }*/

    }

}

/*
void SetPositivePoints_Unst(const int& space_order, dTensor2& spts)
{

    // Positivity points with duplicates removed.  I am unaware of a clever
    // solution to construct these without saving duplicates
    //
    // TODO - hard code the points for second, fourth and fifth-order case
    //

    // -------------------------------------------------------------- //
    // Positivity points for the general case.  There are many duplicates that
    // will be saved if you enter this part of the code!
    // -------------------------------------------------------------- //

    // Gaussian Quadrature points
    dTensor1 s1d_ga(space_order);
    setGaussPoints1d( s1d_ga );

    // Gauss lobatto points
    dTensor1 s1d_gl(2);  
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
    dTensor3 uu( 2, spts.getsize(1), 2 );
    dTensor3 vv( 2, spts.getsize(1), 2 );
    for (int m=1; m<=2; m++)
    for (int k=1; k<=space_order; k++)
    {
        uu.set(m,k, 1, s1d_gl.get(m) );
        vv.set(m,k, 2, s1d_ga.get(k) );
    }

    // Weird mapping from unit square to triangles
    int z = 0;
    for( int m=1; m<=2; m++)
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
        
        spts.set(19,1,  0.112615157582632 );
        spts.set(19,2,  0.112615157582632 );
        spts.set(20,1, -0.225230315165263 );
        spts.set(20,2,  0.112615157582632 );
        spts.set(21,1,  0.112615157582632 );
        spts.set(21,2, -0.225230315165263 );
        spts.set(22,1, -0.241757119823562 );
        spts.set(22,2, -0.241757119823562 );
        spts.set(23,1,  0.483514239647126 );
        spts.set(23,2, -0.241757119823562 );
        spts.set(24,1, -0.241757119823562 );
        spts.set(24,2,  0.483514239647126 );

}*/


void SetPositivePoints_Unst(const int& space_order, dTensor2& spts)
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

        /*spts.set(19,1,  0.112615157582632 );
        spts.set(19,2,  0.112615157582632 );
        spts.set(20,1, -0.225230315165263 );
        spts.set(20,2,  0.112615157582632 );
        spts.set(21,1,  0.112615157582632 );
        spts.set(21,2, -0.225230315165263 );
        spts.set(22,1, -0.241757119823562 );
        spts.set(22,2, -0.241757119823562 );
        spts.set(23,1,  0.483514239647126 );
        spts.set(23,2, -0.241757119823562 );
        spts.set(24,1, -0.241757119823562 );
        spts.set(24,2,  0.483514239647126 );*/
/*
      spts.set(1,1,  0.112615157582632 );
        spts.set(1,2,  0.112615157582632 );
        spts.set(2,1, -0.225230315165263 );
        spts.set(2,2,  0.112615157582632 );
        spts.set(3,1,  0.112615157582632 );
        spts.set(3,2, -0.225230315165263 );
        spts.set(4,1, -0.241757119823562 );
        spts.set(4,2, -0.241757119823562 );
        spts.set(5,1,  0.483514239647126 );
        spts.set(5,2, -0.241757119823562 );
        spts.set(6,1, -0.241757119823562 );
        spts.set(6,2,  0.483514239647126 );*/

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
void SamplePhiAtPositivePoints_Unst(const int& space_order, 
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
