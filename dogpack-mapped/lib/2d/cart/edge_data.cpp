#include "edge_data.h"

edge_data::edge_data()
    // Constructor
{
    const int mpoints1d = dogParams.get_space_order();
    const int kmax = dogParams.get_kmax();
    const int kmax_divfree = dogParams.get_kmax_divfree();

    KMAX_MAX = 25;

    // 1D Gaussian quadrature weights and points
    wgts1d = new dTensor1(mpoints1d); wgts1d->setall(0.);
    xpts1d = new dTensor1(mpoints1d); xpts1d->setall(0.);

    // Legendre basis functions
    phi_xl = new dTensor2(mpoints1d,KMAX_MAX); phi_xl->setall(0.);
    phi_xr = new dTensor2(mpoints1d,KMAX_MAX); phi_xr->setall(0.);
    phi_yl = new dTensor2(mpoints1d,KMAX_MAX); phi_yl->setall(0.);
    phi_yr = new dTensor2(mpoints1d,KMAX_MAX); phi_yr->setall(0.);

    // weights times Legendre basis functions
    wght_phi_xl = new dTensor2(mpoints1d,KMAX_MAX); wght_phi_xl->setall(0.);
    wght_phi_xr = new dTensor2(mpoints1d,KMAX_MAX); wght_phi_xr->setall(0.);
    wght_phi_yl = new dTensor2(mpoints1d,KMAX_MAX); wght_phi_yl->setall(0.);
    wght_phi_yr = new dTensor2(mpoints1d,KMAX_MAX); wght_phi_yr->setall(0.);
}

edge_data::~edge_data()
    // Destructor
{
    delete wgts1d;
    delete xpts1d;

    // Legendre basis functions
    delete phi_xl;
    delete phi_xr;
    delete phi_yl;
    delete phi_yr;
}

// this code was moved from SetEdgeData
void edge_data::init()
{
    edge_data& Edge_Data = *this;
    const int morder = dogParams.get_space_order();
    const int kmax = dogParams.get_kmax();

    double xi,xi2,xi3,xi4;
    double eta,eta2,eta3,eta4;
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double length = sqrt(dx*dx + dy*dy);
    const double lx5    = sqrt(5.0*dx*dx + dy*dy);
    const double ly5    = sqrt(dx*dx + 5.0*dy*dy);

    // ---------------------------------
    // Quick error check
    // ---------------------------------
    assert_printf(kmax==(morder*morder), "\n"
            "       morder = %d\n"
            "         kmax = %d\n",
            morder,kmax);

    // ---------------------------------
    // Set quadrature weights and points
    // ---------------------------------
    switch ( morder )
    {
        case 1:
            wgts1d->set(1, 2.0 );

            xpts1d->set(1, 0.0 );
            break;

        case 2:
            wgts1d->set(1,  1.0 );
            wgts1d->set(2,  1.0 );

            xpts1d->set(1, -1.0/sq3 );
            xpts1d->set(2,  1.0/sq3 );
            break;

        case 3:
            wgts1d->set(1,  5.0/9.0 );
            wgts1d->set(2,  8.0/9.0 );
            wgts1d->set(3,  5.0/9.0 );

            xpts1d->set(1,  sq3/sq5 );
            xpts1d->set(2,  0.0 );
            xpts1d->set(3, -sq3/sq5 );
            break;

        case 4:
            wgts1d->set(1, (18.0 - sq3*sq10)/36.0 );
            wgts1d->set(2, (18.0 + sq3*sq10)/36.0 );
            wgts1d->set(3, wgts1d->get(2) );
            wgts1d->set(4, wgts1d->get(1) );

            xpts1d->set(1,  sqrt(3.0 + sqrt(4.8))/(-sq7) );
            xpts1d->set(2,  sqrt(3.0 - sqrt(4.8))/(-sq7) );
            xpts1d->set(3, -xpts1d->get(2) );
            xpts1d->set(4, -xpts1d->get(1) );           
            break;

        case 5:      
            wgts1d->set(1, (322.0 - 13.0*sq7*sq10)/900.0 );
            wgts1d->set(2, (322.0 + 13.0*sq7*sq10)/900.0 );
            wgts1d->set(3, 128.0/225.0 );
            wgts1d->set(4, wgts1d->get(2) );
            wgts1d->set(5, wgts1d->get(1) );

            xpts1d->set(1,  sqrt(5.0 + 2.0*sq10/sq7)/(-3.0) );
            xpts1d->set(2,  sqrt(5.0 - 2.0*sq10/sq7)/(-3.0) );
            xpts1d->set(3,  0.0 );
            xpts1d->set(4, -xpts1d->get(2) );
            xpts1d->set(5, -xpts1d->get(1) );
            break;
    }

    // ----------------------------------------------
    // Evaluate Legendre basis functions on the edges
    // ----------------------------------------------
    for (int m=1; m<=morder; m++)
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

        phi_xl->set( m,25, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*(105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0) );
        phi_xl->set( m,24, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq7*(2.5*xi3 - 1.5*xi) );
        phi_xl->set( m,23, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq7*(2.5*eta3 - 1.5*eta) );
        phi_xl->set( m,22, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq5*(1.5*xi2 - 0.5) );
        phi_xl->set( m,21, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq5*(1.5*eta2 - 0.5) );
        phi_xl->set( m,20, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq3*xi );
        phi_xl->set( m,19, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq3*eta );
        phi_xl->set( m,18, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
        phi_xl->set( m,17, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
        phi_xl->set( m,16, sq7*(2.5*eta3 - 1.5*eta)*sq7*(2.5*xi3 - 1.5*xi) );
        phi_xl->set( m,15, sq7*(2.5*eta3 - 1.5*eta)*sq5*(1.5*xi2 - 0.5) );
        phi_xl->set( m,14, sq7*(2.5*xi3 - 1.5*xi)*sq5*(1.5*eta2 - 0.5) );
        phi_xl->set( m,13, sq7*(2.5*eta3 - 1.5*eta)*sq3*xi );
        phi_xl->set( m,12, sq7*(2.5*xi3 - 1.5*xi)*sq3*eta );
        phi_xl->set( m,11, sq7*(2.5*eta3 - 1.5*eta) );
        phi_xl->set( m,10,  sq7*(2.5*xi3 - 1.5*xi) );
        phi_xl->set( m,9, 5.0*(1.5*xi2 - 0.5)*(1.5*eta2 - 0.5) );
        phi_xl->set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
        phi_xl->set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );
        phi_xl->set( m,6,  sq5*(1.5*eta2 - 0.5) );
        phi_xl->set( m,5,  sq5*(1.5*xi2 - 0.5) );
        phi_xl->set( m,4,  3.0*xi*eta );
        phi_xl->set( m,3, sq3*eta );
        phi_xl->set( m,2, sq3*xi  );
        phi_xl->set( m,1, 1.0 );


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

        phi_xr->set( m,25, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*(105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0) );
        phi_xr->set( m,24, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq7*(2.5*xi3 - 1.5*xi) );
        phi_xr->set( m,23, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq7*(2.5*eta3 - 1.5*eta) );
        phi_xr->set( m,22, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq5*(1.5*xi2 - 0.5) );
        phi_xr->set( m,21, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq5*(1.5*eta2 - 0.5) );
        phi_xr->set( m,20, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq3*xi );
        phi_xr->set( m,19, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq3*eta );
        phi_xr->set( m,18, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
        phi_xr->set( m,17, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
        phi_xr->set( m,16, sq7*(2.5*eta3 - 1.5*eta)*sq7*(2.5*xi3 - 1.5*xi) );
        phi_xr->set( m,15, sq7*(2.5*eta3 - 1.5*eta)*sq5*(1.5*xi2 - 0.5) );
        phi_xr->set( m,14, sq7*(2.5*xi3 - 1.5*xi)*sq5*(1.5*eta2 - 0.5) );
        phi_xr->set( m,13, sq7*(2.5*eta3 - 1.5*eta)*sq3*xi );
        phi_xr->set( m,12, sq7*(2.5*xi3 - 1.5*xi)*sq3*eta );
        phi_xr->set( m,11, sq7*(2.5*eta3 - 1.5*eta) );
        phi_xr->set( m,10,  sq7*(2.5*xi3 - 1.5*xi) );
        phi_xr->set( m,9, 5.0*(1.5*xi2 - 0.5)*(1.5*eta2 - 0.5) );
        phi_xr->set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
        phi_xr->set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );
        phi_xr->set( m,6,  sq5*(1.5*eta2 - 0.5) );
        phi_xr->set( m,5,  sq5*(1.5*xi2 - 0.5) );
        phi_xr->set( m,4,  3.0*xi*eta );
        phi_xr->set( m,3, sq3*eta );
        phi_xr->set( m,2, sq3*xi  );
        phi_xr->set( m,1, 1.0 );

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

        phi_yl->set( m,25, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*(105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0) );
        phi_yl->set( m,24, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq7*(2.5*xi3 - 1.5*xi) );
        phi_yl->set( m,23, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq7*(2.5*eta3 - 1.5*eta) );
        phi_yl->set( m,22, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq5*(1.5*xi2 - 0.5) );
        phi_yl->set( m,21, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq5*(1.5*eta2 - 0.5) );
        phi_yl->set( m,20, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq3*xi );
        phi_yl->set( m,19, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq3*eta );
        phi_yl->set( m,18, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
        phi_yl->set( m,17, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
        phi_yl->set( m,16, sq7*(2.5*eta3 - 1.5*eta)*sq7*(2.5*xi3 - 1.5*xi) );
        phi_yl->set( m,15, sq7*(2.5*eta3 - 1.5*eta)*sq5*(1.5*xi2 - 0.5) );
        phi_yl->set( m,14, sq7*(2.5*xi3 - 1.5*xi)*sq5*(1.5*eta2 - 0.5) );
        phi_yl->set( m,13, sq7*(2.5*eta3 - 1.5*eta)*sq3*xi );
        phi_yl->set( m,12, sq7*(2.5*xi3 - 1.5*xi)*sq3*eta );
        phi_yl->set( m,11, sq7*(2.5*eta3 - 1.5*eta) );
        phi_yl->set( m,10,  sq7*(2.5*xi3 - 1.5*xi) );
        phi_yl->set( m,9, 5.0*(1.5*xi2 - 0.5)*(1.5*eta2 - 0.5) );
        phi_yl->set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
        phi_yl->set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );
        phi_yl->set( m,6,  sq5*(1.5*eta2 - 0.5) );
        phi_yl->set( m,5,  sq5*(1.5*xi2 - 0.5) );
        phi_yl->set( m,4,  3.0*xi*eta );
        phi_yl->set( m,3, sq3*eta );
        phi_yl->set( m,2, sq3*xi  );
        phi_yl->set( m,1, 1.0 );


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

        phi_yr->set( m,25, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*(105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0) );
        phi_yr->set( m,24, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq7*(2.5*xi3 - 1.5*xi) );
        phi_yr->set( m,23, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq7*(2.5*eta3 - 1.5*eta) );
        phi_yr->set( m,22, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq5*(1.5*xi2 - 0.5) );
        phi_yr->set( m,21, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq5*(1.5*eta2 - 0.5) );
        phi_yr->set( m,20, (105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0)*sq3*xi );
        phi_yr->set( m,19, (105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0)*sq3*eta );
        phi_yr->set( m,18, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
        phi_yr->set( m,17, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
        phi_yr->set( m,16, sq7*(2.5*eta3 - 1.5*eta)*sq7*(2.5*xi3 - 1.5*xi) );
        phi_yr->set( m,15, sq7*(2.5*eta3 - 1.5*eta)*sq5*(1.5*xi2 - 0.5) );
        phi_yr->set( m,14, sq7*(2.5*xi3 - 1.5*xi)*sq5*(1.5*eta2 - 0.5) );
        phi_yr->set( m,13, sq7*(2.5*eta3 - 1.5*eta)*sq3*xi );
        phi_yr->set( m,12, sq7*(2.5*xi3 - 1.5*xi)*sq3*eta );
        phi_yr->set( m,11, sq7*(2.5*eta3 - 1.5*eta) );
        phi_yr->set( m,10,  sq7*(2.5*xi3 - 1.5*xi) );
        phi_yr->set( m,9, 5.0*(1.5*xi2 - 0.5)*(1.5*eta2 - 0.5) );
        phi_yr->set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
        phi_yr->set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );
        phi_yr->set( m,6,  sq5*(1.5*eta2 - 0.5) );
        phi_yr->set( m,5,  sq5*(1.5*xi2 - 0.5) );
        phi_yr->set( m,4,  3.0*xi*eta );
        phi_yr->set( m,3, sq3*eta );
        phi_yr->set( m,2, sq3*xi  );
        phi_yr->set( m,1, 1.0 );


    }
    for (int m=1; m<=morder; m++)
        for (int k=1; k<=KMAX_MAX; k++)
        {
            wght_phi_xl->set(m,k, wgts1d->get(m)*phi_xl->get(m,k));
            wght_phi_xr->set(m,k, wgts1d->get(m)*phi_xr->get(m,k));
            wght_phi_yl->set(m,k, wgts1d->get(m)*phi_yl->get(m,k));
            wght_phi_yr->set(m,k, wgts1d->get(m)*phi_yr->get(m,k));
        }

    // print points and phi values
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
