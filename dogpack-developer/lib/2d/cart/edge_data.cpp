#include <new>
#include <cmath>
#include "constants.h"
#include "tensors.h"
#include "assert.h"
#include "debug.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "edge_data.h"

// Constructor
edge_data::edge_data()
{

    const int mpoints1d = dogParams.get_space_order();
    const int kmax = dogParams.get_kmax();
    const int kmax_divfree = dogParams.get_kmax_divfree();

    // MOVED TO HEADER FILE -DS
//  MAX_KMAX = 15;

    // 1D Gaussian quadrature weights and points
    wgts1d      = new dTensor1(mpoints1d); wgts1d->setall(0.);
    xpts1d      = new dTensor1(mpoints1d); xpts1d->setall(0.);

    // Legendre basis functions
    phi_xl      = new dTensor2(mpoints1d, MAX_KMAX); phi_xl->setall(0.);
    phi_xr      = new dTensor2(mpoints1d, MAX_KMAX); phi_xr->setall(0.);
    phi_yl      = new dTensor2(mpoints1d, MAX_KMAX); phi_yl->setall(0.);
    phi_yr      = new dTensor2(mpoints1d, MAX_KMAX); phi_yr->setall(0.);

    // weights times Legendre basis functions
    wght_phi_xl = new dTensor2(mpoints1d, MAX_KMAX); wght_phi_xl->setall(0.);
    wght_phi_xr = new dTensor2(mpoints1d, MAX_KMAX); wght_phi_xr->setall(0.);
    wght_phi_yl = new dTensor2(mpoints1d, MAX_KMAX); wght_phi_yl->setall(0.);
    wght_phi_yr = new dTensor2(mpoints1d, MAX_KMAX); wght_phi_yr->setall(0.);

}

// Destructor
edge_data::~edge_data()
{

    // 1D weights and points
    delete wgts1d;
    delete xpts1d;

    // Legendre basis functions
    delete phi_xl;
    delete phi_xr;
    delete phi_yl;
    delete phi_yr;

    // weights times Legendre basis functions
    delete wght_phi_xl;
    delete wght_phi_xr;
    delete wght_phi_yl;
    delete wght_phi_yr;

};

// this code was moved from SetEdgeData
void edge_data::init()
{

    edge_data& Edge_Data = *this;
    const int morder     = dogParams.get_space_order();
    const int kmax       = dogParams.get_kmax();

    // Grid information.  (In general, this could depend on local values, for
    //                     example in an AMR code.)
    const double dx      = dogParamsCart2.get_dx();
    const double dy      = dogParamsCart2.get_dy();

// I don't see where these are used ... (6/9/2014 -DS)
//
//  const double length = sqrt(dx*dx + dy*dy);
//  const double lx5    = sqrt(5.0*dx*dx + dy*dy);
//  const double ly5    = sqrt(dx*dx + 5.0*dy*dy);

    // ---------------------------------
    // Quick error check
    // ---------------------------------
    assert_printf( kmax==(morder*(morder+1))/2, "\n"
            "       morder = %d\n"
            "         kmax = %d\n",
            morder, kmax);

    // Compute the 1D weights and points
    void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);
    setGaussPoints1d(*wgts1d, *xpts1d);

    // ----------------------------------------------
    // Evaluate Legendre basis functions on the edges
    // ----------------------------------------------
    double xi,xi2,xi3,xi4;
    double eta,eta2,eta3,eta4;
    for( int m=1; m<=morder; m++ )
    {

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Right edge (will be Left state in Riemann problem)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        xi  = 1.0;
        eta = xpts1d->get(m);

        xi2 = xi*xi;
        xi3 = xi2*xi;
        xi4 = xi3*xi;

        eta2 = eta*eta;
        eta3 = eta2*eta;
        eta4 = eta3*eta;

        // fifth-order terms
        phi_xl->set( m, 15,  105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
        phi_xl->set( m, 14,  105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
        phi_xl->set( m, 13,  5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
        phi_xl->set( m, 12,  sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
        phi_xl->set( m, 11,  sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

        // fourth-order terms
        phi_xl->set( m, 10,  sq7*(2.5*eta3 - 1.5*eta) );
        phi_xl->set( m, 9,   sq7*(2.5*xi3 - 1.5*xi) );
        phi_xl->set( m, 8,   sq3*sq5*xi*(1.5*eta2 - 0.5) );
        phi_xl->set( m, 7,   sq3*sq5*eta*(1.5*xi2 - 0.5) );

        // third-order terms
        phi_xl->set( m, 6,   sq5*(1.5*eta2 - 0.5) );
        phi_xl->set( m, 5,   sq5*(1.5*xi2 - 0.5) );
        phi_xl->set( m, 4,   3.0*xi*eta );                  

        // second-order terms
        phi_xl->set( m, 3,   sq3*eta );
        phi_xl->set( m, 2,   sq3*xi  );

        // first-order terms
        phi_xl->set( m, 1,   1.0 );

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Left edge (will be Right state in Riemann problem)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        xi  = -1.0;
        eta = xpts1d->get(m);

        xi2 = xi*xi;
        xi3 = xi2*xi;
        xi4 = xi3*xi;

        eta2 = eta*eta;
        eta3 = eta2*eta;
        eta4 = eta3*eta;

        // fifth-order terms
        phi_xr->set( m,15, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
        phi_xr->set( m,14, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
        phi_xr->set( m,13, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
        phi_xr->set( m,12, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
        phi_xr->set( m,11, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

        // fourth-order terms
        phi_xr->set( m,10, sq7*(2.5*eta3 - 1.5*eta) );
        phi_xr->set( m,9,  sq7*(2.5*xi3 - 1.5*xi) );
        phi_xr->set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
        phi_xr->set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );        

        // third-order terms
        phi_xr->set( m,6,  sq5*(1.5*eta2 - 0.5) );
        phi_xr->set( m,5,  sq5*(1.5*xi2 - 0.5) );
        phi_xr->set( m,4,  3.0*xi*eta );                  

        // second-order terms
        phi_xr->set( m,3,  sq3*eta );
        phi_xr->set( m,2,  sq3*xi  );

        // first-order terms
        phi_xr->set( m,1,  1.0 );

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Top edge (will be Left state in Riemann problem)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        xi  = xpts1d->get(m);
        eta = 1.0;

        xi2 = xi*xi;
        xi3 = xi2*xi;
        xi4 = xi3*xi;

        eta2 = eta*eta;
        eta3 = eta2*eta;
        eta4 = eta3*eta;

        // fifth-order terms
        phi_yl->set( m,15, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
        phi_yl->set( m,14, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
        phi_yl->set( m,13, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
        phi_yl->set( m,12, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
        phi_yl->set( m,11, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

        // fourth-order terms
        phi_yl->set( m,10, sq7*(2.5*eta3 - 1.5*eta) );
        phi_yl->set( m,9,  sq7*(2.5*xi3 - 1.5*xi) );
        phi_yl->set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
        phi_yl->set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );

        // third-order terms
        phi_yl->set( m,6,  sq5*(1.5*eta2 - 0.5) );
        phi_yl->set( m,5,  sq5*(1.5*xi2 - 0.5) );
        phi_yl->set( m,4,  3.0*xi*eta );                  

        // second-order terms
        phi_yl->set( m,3,  sq3*eta );
        phi_yl->set( m,2,  sq3*xi  );

        // first-order terms
        phi_yl->set( m,1,  1.0 );

        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Bottom edge (will be Right state in Riemann problem)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        xi  = xpts1d->get(m);
        eta = -1.0;

        xi2 = xi*xi;
        xi3 = xi2*xi;
        xi4 = xi3*xi;

        eta2 = eta*eta;
        eta3 = eta2*eta;
        eta4 = eta3*eta;

        // fifth-order terms
        phi_yr->set( m,15, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
        phi_yr->set( m,14, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
        phi_yr->set( m,13, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
        phi_yr->set( m,12, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
        phi_yr->set( m,11, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

        // fourth-order terms
        phi_yr->set( m,10, sq7*(2.5*eta3 - 1.5*eta) );
        phi_yr->set( m,9,  sq7*(2.5*xi3 - 1.5*xi) );
        phi_yr->set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
        phi_yr->set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );

        // third-order terms
        phi_yr->set( m,6,  sq5*(1.5*eta2 - 0.5) );
        phi_yr->set( m,5,  sq5*(1.5*xi2 - 0.5) );
        phi_yr->set( m,4,  3.0*xi*eta );                  

        // second-order terms
        phi_yr->set( m,3, sq3*eta );
        phi_yr->set( m,2, sq3*xi  );

        // first-order terms
        phi_yr->set( m,1, 1.0 );

    }

    // Save the product of the weights times the basis functions.  This is the
    // part that is required to perform integration along each edge.
    for( int m=1; m <= morder; m++)
    for( int k=1; k <= MAX_KMAX; k++)
    {
        wght_phi_xl->set(m,k, wgts1d->get(m)*phi_xl->get(m,k) );
        wght_phi_xr->set(m,k, wgts1d->get(m)*phi_xr->get(m,k) );
        wght_phi_yl->set(m,k, wgts1d->get(m)*phi_yl->get(m,k) );
        wght_phi_yr->set(m,k, wgts1d->get(m)*phi_yr->get(m,k) );
    }

    // Debugging portion.  Print points and phi values to standard output
    // (compare with Legendre2d::setPositivityPoints)
    if(debug3)
    {
        printf("\n phi_xl:");
        for (int m=1; m<=morder; m++)
        {
            printf("\n  phi_xl(%d,*) at (%24.16e, %24.16e)",m,1.,xpts1d->get(m));
            for (int k=1; k<=kmax; k++)
            {
                printf("\n   phi_xl(%d,%d) = %24.16e", m,k,phi_xl->get(m, k));
            }
        }
        printf("\n phi_xr:");
        for (int m=1; m<=morder; m++)
        {
            printf("\n  phi_xr(%d,*) at (%24.16e,  %24.16e)",m,-1.,xpts1d->get(m));
            for (int k=1; k<=kmax; k++)
            {
                printf("\n   phi_xr(%d,%d) = %24.16e", m,k,phi_xr->get(m, k));
            }
        }
        printf("\n phi_yl:");
        for (int m=1; m<=morder; m++)
        {
            printf("\n  phi_yl(%d,*) at (%24.16e, %24.16e)",m,xpts1d->get(m),1.);
            for (int k=1; k<=kmax; k++)
            {
                printf("\n   phi_yl(%d,%d) = %24.16e", m,k,phi_yl->get(m, k));
            }
        }
        printf("\n phi_yr:");
        for (int m=1; m<=morder; m++)
        {
            printf("\n  phi_yr(%d,*) at (%24.16e, %24.16e)",m,xpts1d->get(m),-1.);
            for (int k=1; k<=kmax; k++)
            {
                printf("\n   phi_yr(%d,%d) = %24.16e", m,k,phi_yr->get(m, k));
            }
        }
        printf("\n");
    }

}
