#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "mesh.h"                   // Required for all unstructured stuff
#include "Quadrature.h"             // 1D Quadrature rules
#include "MonomialsToLegendre.h"    // For evaluating basis functions

// This is the positivity preserving limiter proposed in 
// "Maximum-Principle-Satisfying and Positivity-Preserving
// High Order Discontinuous Galerkin Schemes
// for Conservation Laws on Triangular Meshes", Zhang, Xia and Shu
// J. Sci. Comput. (2012).
//
// THIS METHOD ASSUMES THAT EVERY COMPONENT OF CONSERVED VARIABLES SHOULD STAY
// POSITIVE.
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

    // Do nothing in the case of piecewise constants
    if( space_order == 1 )
    { return; }

    // ------------------------------------------------ //
    // number of points where we want to check solution //
    // ------------------------------------------------ //
//  const int mpoints_on_edges   = 4*space_order;
//  const int mpoints_L2Proj     = (space_order)*(space_order);
    const int mpoints = 3*space_order*space_order;  // TODO - HARD CODE THESE INSTEAD!

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
    for(int me=1; me <= meqn; me++)
    {

        double m = 0.0;
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double qnow = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                qnow += q.get(i,me,k) * phi.get(mp,k);
            }
            m = Min(m, qnow);
        }

        double theta = 0.0;
        double Q1 = q.get(i,me,1);  assert_ge( Q1, -1e-13 );
        if( fabs( Q1 - m ) < 1.0e-14 ){ theta = 1.0; }
        else{ theta = Min( 1.0, fabs( Q1 / (Q1 - m) ) ); }

        // limit q //
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,me,k, q.get(i,me,k) * theta );
        }

    }

}


void SetPositivePoints_Unst(const int& space_order, dTensor2& spts)
{

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
