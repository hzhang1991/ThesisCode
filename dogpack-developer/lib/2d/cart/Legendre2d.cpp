// Legendre2d.cpp.
//
// This file describes common operations performed one the 2D (Cartesian) basis
// functions.
//
// See also: edge_data.cpp

#include <cmath>
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "Interval.h"
#include "dog_math.h"
#include "edge_data.h"
#include "Quadrature.h"             // 1D quadrature points/weights
#include "Legendre2d.h"

// -------------------------------------------------------------------------- //
// Routine for evaluating DG basis function at a list of points.
//
// Input:
// ------
//
//      spts( 1:numpts, 1:2 )  - points (in canonical coordinates)
//
//          xi  = xpts( :, 1 ) - mapped from (x-x_i)/dx = \xi  / 2
//          eta = xpts( :, 2 ) - mapped from (y-y_j)/dy = \eta / 2 
//
// Returns:
// --------
//
//      phi( 1:numpts, 1:kmax ) - basis functions evaluated at each point.  It
//                                is expected that kmax is consisten with 
//                                dogParamsCart2.get_kmax()
//
// See also: SetLegendreGrad.
// -------------------------------------------------------------------------- //
void SetLegendreAtPoints(const dTensor2& spts, dTensor2& phi)
{

    const int space_order   = dogParams.get_space_order();
    const int kmax          = phi.getsize(2);            // assert_eq(kmax, dogParams.get_kmax());
    const int mpoints       = spts.getsize(1);           // assert_eq(mpoints, phi.getsize(1));

    for (int m=1; m<=mpoints; m++)
    {

        // grid point (x, y) <-> (xi, eta )
        const double xi  = spts.get(m,1);
        const double eta = spts.get(m,2);
        const double xi2 = xi*xi;
        const double xi3 = xi*xi2;
        const double xi4 = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;     


        // Legendre basis functions at each gaussian quadrature point in the
        // interval [-1,1]x[-1,1].
        switch( kmax )
        {
            case 15:  // fifth order                                 
                phi.set( m,15, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
                phi.set( m,14, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
                phi.set( m,13, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
                phi.set( m,12, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
                phi.set( m,11, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

            case 10:  // fourth order
                phi.set( m,10, sq7*(2.5*eta3 - 1.5*eta) );
                phi.set( m,9,  sq7*(2.5*xi3 - 1.5*xi) );
                phi.set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
                phi.set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );

            case 6:  // third order
                phi.set( m,6,  sq5*(1.5*eta2 - 0.5) );
                phi.set( m,5,  sq5*(1.5*xi2 - 0.5) );
                phi.set( m,4,  3.0*xi*eta );                  

            case 3:  // second order                
                phi.set( m,3, sq3*eta );
                phi.set( m,2, sq3*xi  );

            case 1:  // first order
                phi.set( m,1, 1.0 );

                break;                
            default:
                unsupported_value_error(space_order);
        }

    }
}

// -------------------------------------------------------------------------- //
// Routine for evaluating DG basis function at a list of points.  This routine
// is identical to the previous one, but this one allows you to specifically
// assign the polynomial order before going to evaluate the basis functions.
// That is, it is assumed that space_order <= sphi.getsize(2), but this does
// not need to be the case.
//
// Input:
// ------
//
//      space_order            - Order of the method.
//
//      spts( 1:numpts, 1:2 )  - points (in canonical coordinates)
//
//          xi  = xpts( :, 1 ) - mapped from (x-x_i)/dx = \xi  / 2
//          eta = xpts( :, 2 ) - mapped from (y-y_j)/dy = \eta / 2 
//
// Returns:
// --------
//
//      phi( 1:numpts, 1:kmax ) - basis functions evaluated at each point.  It
//                                is expected that kmax is large enough to
//                                accomodate whatever space_order is fed in.
//
// See also: SetLegendreGrad.
// -------------------------------------------------------------------------- //
void SetLegendreAtPoints(const int& space_order, const dTensor2& spts, dTensor2& phi)
{
    // Number of points to sample
    const int mpoints = spts.getsize(1);

    for (int m=1; m<=mpoints; m++)
    {
        // grid point (x,y)
        const double xi   = spts.get(m,1);
        const double eta  = spts.get(m,2);
        const double xi2  = xi*xi;
        const double xi3  = xi*xi2;
        const double xi4  = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;     

        // Legendre basis functions at each gaussian quadrature point in the
        // interval [-1,1]x[-1,1].
        switch( space_order )
        {
            case 5:  // fifth order                                 
                phi.set( m,15, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
                phi.set( m,14, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
                phi.set( m,13, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
                phi.set( m,12, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
                phi.set( m,11, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

            case 4:  // fourth order
                phi.set( m,10, sq7*(2.5*eta3 - 1.5*eta) );
                phi.set( m,9,  sq7*(2.5*xi3 - 1.5*xi) );
                phi.set( m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
                phi.set( m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );

            case 3:  // third order
                phi.set( m,6,  sq5*(1.5*eta2 - 0.5) );
                phi.set( m,5,  sq5*(1.5*xi2 - 0.5) );
                phi.set( m,4,  3.0*xi*eta );                  

            case 2:  // second order                
                phi.set( m,3, sq3*eta );
                phi.set( m,2, sq3*xi  );

            case 1:  // first order
                phi.set( m,1, 1.0 );

                break;                

            default:
                unsupported_value_error(space_order);
        }
    }
}

// -------------------------------------------------------------------------- //
// Routine for evaluating the derivatives of the DG basis function at a list 
// of points.
//
// Input:
// ------
//
//   dx, dy - size of the element considered.  These are required to compute
//            derivatives.
//
//   
//   spts( 1:numpts, 1:2 )  - points (in canonical coordinates)
//
//          xi  = xpts( :, 1 ) 
//          eta = xpts( :, 2 )
//
// Returns:
// --------
//
// Derivatives of the basis functions:
//
//      phi_x( 1:numpts, 1:kmax ) - \pd{ phi^k }{ x }
//      phi_y( 1:numpts, 1:kmay ) - \pd{ phi^k }{ y }
//
// See also: SetLegendreAtPoints.
// -------------------------------------------------------------------------- //
void SetLegendreGrad(const double dx, const double dy,
        const dTensor2& spts, dTensor2& phi_x, dTensor2& phi_y)
{

    const int space_order = dogParams.get_space_order();
    const int mpoints     = spts.getsize(1);
    const int kmax        = dogParams.get_kmax();

    // Quick error check.  
    // TODO - do we always want this turned on? -DS
//  if( debug2 )
//  {
//      assert_eq(phi_x.getsize(1), mpoints);
//      assert_eq(kmax, phi_x.getsize(2) );
//      assert_eq(kmax, phi_y.getsize(2) );
//      const int mpoints1d = int(sqrt(mpoints));//space_order-1;
//      assert_eq(mpoints, mpoints1d*mpoints1d);
//  }

    // Common term that can be pulled out of the evaluation below
    const double tmpx = 2.0/dx;
    const double tmpy = 2.0/dy;

    for (int m=1; m<=(mpoints); m++)
    {

        // grid point (x,y) <-> (xi, eta)
        const double xi   = spts.get(m,1);
        const double eta  = spts.get(m,2);
        const double xi2  = xi*xi;
        const double xi3  = xi*xi2;
        const double xi4  = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;     

        // Gradient of Legendre basis functions at each gaussian quadrature point
        switch( kmax )
        {
            case 15:  // fifth order
                phi_x.set( m,15, 0.0 );
                phi_x.set( m,14, tmpx*(105.0/2.0*xi3  - 45.0/2.0*xi) );
                phi_x.set( m,13, tmpx*15.0*(1.5*eta2 - 0.5)*xi );
                phi_x.set( m,12, tmpx*sq3*sq7*(2.5*eta3 - 1.5*eta) );
                phi_x.set( m,11, tmpx*sq3*sq7*(7.5*xi2 - 1.5)*eta );

                phi_y.set( m,15, tmpy*(105.0/2.0*eta3  - 45.0/2.0*eta) );
                phi_y.set( m,14, 0.0 );
                phi_y.set( m,13, tmpy*15.0*(1.5*xi2 - 0.5)*eta );
                phi_y.set( m,12, tmpy*sq3*sq7*(7.5*eta2 - 1.5)*xi );
                phi_y.set( m,11, tmpy*sq3*sq7*(2.5*xi3 - 1.5*xi) );

            case 10:  // fourth order
                phi_x.set( m,10, 0.0 );
                phi_x.set( m,9,  tmpx*sq7*(7.5*xi2 - 1.5) );                  
                phi_x.set( m,8,  tmpx*sq3*sq5*(1.5*eta2 - 0.5) );
                phi_x.set( m,7,  tmpx*sq3*sq5*eta*3.0*xi );

                phi_y.set( m,10, tmpy*sq7*(7.5*eta2 - 1.5) );
                phi_y.set( m,9,  0.0 );
                phi_y.set( m,8,  tmpy*sq3*sq5*eta*3.0*xi );
                phi_y.set( m,7,  tmpy*sq3*sq5*(1.5*xi2 - 0.5) );

            case 6:  // third order
                phi_x.set( m,6,  0.0 );
                phi_x.set( m,5,  tmpx*sq5*3.0*xi );                   
                phi_x.set( m,4,  tmpx*3.0*eta );

                phi_y.set( m,6,  tmpy*sq5*3.0*eta );                  
                phi_y.set( m,5,  0.0 );
                phi_y.set( m,4,  tmpy*3.0*xi );

            case 3:  // second order
                phi_x.set( m,3,  0.0 );
                phi_x.set( m,2,  tmpx*sq3 );

                phi_y.set( m,3,  tmpy*sq3 );
                phi_y.set( m,2,  0.0 );

            case 1:  // first order
                phi_x.set( m,1,  0.0 );

                phi_y.set( m,1,  0.0 );
                break;

            default:
                unsupported_value_error(space_order);

        }
    }
}

// -------------------------------------------------------------------------- //
// Routine for evaluating the second derivative of the DG basis function.
//
// Input:
// ------
//
//   dx, dy - size of the element considered.  These are required to compute
//            derivatives.
//
//   
//   spts( 1:numpts, 1:2 )  - points (in canonical coordinates)
//
//          xi  = xpts( :, 1 ) 
//          eta = xpts( :, 2 )
//
// Returns:
// --------
//
// Derivatives of the basis functions:
//
//      phi_xx( 1:numpts, 1:kmax ) - \pd2{ phi^k }{ xx }
//      phi_xy( 1:numpts, 1:kmay ) - \pd2{ phi^k }{ xy }
//      phi_yy( 1:numpts, 1:kmay ) - \pd2{ phi^k }{ yy }
//
// See also: SetLegendreAtPoints.
// -------------------------------------------------------------------------- //
void LegendreDiff2(const double dx, const double dy, const dTensor2& spts, 
    dTensor2* phi_xx, dTensor2* phi_xy, dTensor2* phi_yy )
{

    const int space_order = dogParams.get_space_order();
    const int mpoints     = spts.getsize(1);
    const int kmax        = dogParams.get_kmax();

    // Derivatives that pop out from the change of variables
    const double tmpx = 2.0/dx;
    const double tmpy = 2.0/dy;

    const double tmpx2 = tmpx*tmpx;
    const double tmpxy = tmpx*tmpy;
    const double tmpy2 = tmpy*tmpy;

    for (int m=1; m<=(mpoints); m++)
    {

        // grid point (x,y) <-> (xi, eta)
        const double xi   = spts.get(m,1);
        const double eta  = spts.get(m,2);
        const double xi2  = xi*xi;
        const double xi3  = xi*xi2;
        const double xi4  = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;     

        // Gradient of Legendre basis functions at each gaussian quadrature point
        switch( kmax )
        {

            case 15:  // fifth order

                phi_xx->set( m, 15, 0. );
                phi_xx->set( m, 14, tmpx2*0.5*(315.*xi2 - 45.0) );
                phi_xx->set( m, 13, tmpx2*0.5*(45.*eta2 - 15.0) );
                phi_xx->set( m, 12, 0 );
                phi_xx->set( m, 11, tmpx2*15.*sq3*sq7*eta*xi );

                phi_xy->set( m, 15, 0. );
                phi_xy->set( m, 14, 0. );
                phi_xy->set( m, 13, tmpxy*45.*eta*xi );
                phi_xy->set( m, 12, tmpxy*sq3*sq7*0.5*(15.0*eta2 - 3.0) );
                phi_xy->set( m, 11, tmpxy*sq3*sq7*0.5*(15.0*xi2  - 3.0) );

                phi_yy->set( m, 15, tmpy2*0.5*(315.0*eta2 - 45.0) );
                phi_yy->set( m, 14, 0.0 );
                phi_yy->set( m, 13, tmpy2*0.5*(45.*xi2 - 15.0) );
                phi_yy->set( m, 12, tmpy2*15.0*sq3*sq7*eta*xi );
                phi_yy->set( m, 11, 0.0 );

            case 10:  // fourth order

                phi_xx->set( m, 10, 0. );
                phi_xx->set( m, 9, tmpx2*15.*sq7*xi );
                phi_xx->set( m, 8, 0 );
                phi_xx->set( m, 7, tmpx2*3.*sq3*sq5*eta );

                phi_xy->set( m, 10, 0. );
                phi_xy->set( m, 9,  0. );
                phi_xy->set( m, 8, tmpxy*3.*sq3*sq5*eta );
                phi_xy->set( m, 7, tmpxy*3.*sq3*sq5*xi );

                phi_yy->set( m, 10, tmpy2*15.*sq7*eta );
                phi_yy->set( m, 9,   0. );
                phi_yy->set( m, 8,  tmpy2*  3.*sq3*sq5*xi );
                phi_yy->set( m, 7,   0. );

            case 6:  // third order

                phi_xx->set( m, 6, 0.             );
                phi_xx->set( m, 5, tmpx2*3.0*sq5  );
                phi_xx->set( m, 4, 0.             );

                phi_xy->set( m, 6, 0. );
                phi_xy->set( m, 5, 0. );
                phi_xy->set( m, 4, tmpxy*3.       );

                phi_yy->set( m, 6, tmpy2*3.0*sq5  );
                phi_yy->set( m, 5, 0.             );
                phi_yy->set( m, 4, 0.             );

            case 3:  // second order (all second-derivatives vanish)

                phi_xx->set( m, 3, 0. );
                phi_xx->set( m, 2, 0. );
                 
                phi_xy->set( m, 3, 0. );
                phi_xy->set( m, 2, 0. );

                phi_yy->set( m, 3, 0. );
                phi_yy->set( m, 2, 0. );


            case 1:  // first order

                phi_xx->set( m, 1, 0. );
                phi_xy->set( m, 1, 0. );
                phi_yy->set( m, 1, 0. );
                break;

            default:
                unsupported_value_error(space_order);

        }
    }
}

// -------------------------------------------------------------------------- //
// Set (2D) quadrature weights and points.  These are constructed as a tensor
// product of 1D quadrature weights and points:
//
// Input:
// ------
//
//     mpoints1d - number of quadrature points in 1d.
//
// Returns:
// --------
//
//     wgt(1:numpts      ) - 2D quadrature weights.
//    spts(1:numpts, 1:2 ) - 2D quadrature points (in canonical variables).
//                           xi = spts(:,1) and eta = spts(:,2).
//
// See also: Quadrature.cpp.
// -------------------------------------------------------------------------- //
void SetQuadrature(const int mpoints1d, dTensor1& wgt, dTensor2& spts)
{

    const int mpoints = mpoints1d*mpoints1d;
    assert_eq(mpoints, wgt.getsize());
    assert_eq(mpoints, spts.getsize(1));

    // Construct the 1D quadrature weights and points
    dTensor1 w1d(mpoints1d), x1d(mpoints1d);
    setGaussPoints1d(w1d, x1d);

    // Tensor product Gaussian Quadrature
    int k=0;
    for (int m1=1; m1<=(mpoints1d); m1++)
    for (int m2=1; m2<=(mpoints1d); m2++)
    {
        k = k+1;
        wgt.set(k,  w1d.get(m1)*w1d.get(m2) );

        spts.set(k,1, x1d.get(m1) );
        spts.set(k,2, x1d.get(m2) );
    }

}

// -------------------------------------------------------------------------- //
// Compute the extrema of each basis function in the cell.  This is based on an
// L^\infty norm, so the min/max value of the element prescribed by
//
//      q^h = \sum_k Q^{(k)} \phi^{(k)}
//
// Can be estimated by
//
//      q^h >= -\sum_k abs( Q^{(k)} ) phi_max(k)
//
// or 
//
//      q^h <= \sum_k abs( Q^{(k)} ) phi_max(k).
//
// This routine is useful, for example applying the Positivity preserving
// limiter to determine rather quickly whether or not to sample the function at
// each quadrature point.
//
// See also: ???
// -------------------------------------------------------------------------- //
void Legendre2d::SetLegendreIntervals()
{
    const int kmax      = dogParams.get_kmax();
    phi_min_ptr         = new dTensor1(kmax);
    phi_max_ptr         = new dTensor1(kmax);

    // ???
    dTensor1& phi_min   = *phi_min_ptr;
    dTensor1& phi_max   = *phi_max_ptr;

    switch( kmax )
    {
        case 15:  // fifth order                                 
            phi_min.set(15, -9./7. );   phi_max.set(15, 3. );
            phi_min.set(14, -9./7. );   phi_max.set(14, 3. );
            phi_min.set(13, -5./2. );   phi_max.set(13, 5. );
            phi_min.set(12, -sq3*sq7 ); phi_max.set(12, sq3*sq7 );
            phi_min.set(11, -sq3*sq7 ); phi_max.set(11, sq3*sq7 );

        case 10:  // fourth order
            phi_min.set(10, -sq7 );     phi_max.set(10, sq7 );
            phi_min.set( 9, -sq7 );     phi_max.set( 9, sq7 );
            phi_min.set(8,  -sq3*sq5 ); phi_max.set(8,  sq3*sq5 );
            phi_min.set(7,  -sq3*sq5 ); phi_max.set(7,  sq3*sq5 );

        case 6:  // third order
            phi_min.set(6, -sq5/2. );   phi_max.set(6,  sq5 );
            phi_min.set(5, -sq5/2. );   phi_max.set(5,  sq5 );
            phi_min.set(4, -3.0 );      phi_max.set(4,  3.0 );

        case 3:  // second order                
            phi_min.set(3, -sq3 );      phi_max.set(3, sq3 );
            phi_min.set(2, -sq3 );      phi_max.set(2, sq3 );

        case 1:  // first order
            phi_min.set(1, 1.0 );       phi_max.set(1, 1.0 );

            break;
        default:
            unsupported_value_error(kmax);
    }

    phi_interval = new IntervalArray(kmax+1);
    for(int k=1;k<=kmax;k++)
    { phi_interval->fetch(k) = Interval(phi_min.get(k), phi_max.get(k)); }

}
// -------------------------------------------------------------------------- //

// -------------------------------------------------------------------------- //
// TODO - add a description of this funciton.
// -------------------------------------------------------------------------- //
void set_divfree_polys(int mpoints1d, const dTensor2& spts, dTensor3& phi_divfree)
{

    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double dx2 = dx*dx;
    const double dy2 = dy*dy;
    const double dx3 = dx*dx2;
    const double dy3 = dy*dy2;
    const double dx4 = dx*dx3;
    const double dy4 = dy*dy3;
    const double dx5 = dx*dx4;
    const double dy5 = dy*dy4;
    const double dx6 = dx*dx5;
    const double dy6 = dy*dy5;
    const double length = sqrt(dx*dx + dy*dy);
    const double lx5 = sqrt(5.0*dx*dx + dy*dy);
    const double ly5 = sqrt(dx*dx + 5.0*dy*dy);
    const double l4  = sqrt(3.0*dx4 + 45.0*dx2*dy2 + 35.0*dy4);
    const double l6  = sqrt(7.0*dx6 + 107.0*dx4*dy2 + 107.0*dy4*dx2 + 7.0*dy6);

    assert_eq(mpoints1d,dogParams.get_space_order());
    for (int m=1; m<=phi_divfree.getsize(1); m++)
    {
        // grid point (x,y)
        const double xi  = spts.get(m,1);
        const double eta = spts.get(m,2);
        const double xi2 = xi*xi;
        const double xi3 = xi*xi2;
        const double xi4 = xi*xi3;
        const double eta2 = eta*eta;
        const double eta3 = eta*eta2;
        const double eta4 = eta*eta3;

        // Divergence-free Legendre basis functions at each
        // gaussian quadrature point in the interval [-1,1]x[-1,1].
        switch( mpoints1d )
        {
            default:
                unsupported_value_error(mpoints1d);
            case 5:  // fifth order (6 basis vectors)
                eprintf("L2Project_divfree has not been implemented"
                        " for mpoints1d = %d", mpoints1d);
            case 4:  // fourth order (5 basis vectors)
                phi_divfree.set(m,14,1,  0.5*sq7*(5.0*eta3-3.0*eta)  );
                phi_divfree.set(m,14,2,  0.0  );

                phi_divfree.set(m,13,1,  0.0  );
                phi_divfree.set(m,13,2, -0.5*sq7*(5.0*xi3-3.0*xi) );

                phi_divfree.set(m,12,1,  (0.5*sq3*sq5)*(dx/length)*eta*(3.0*xi2-1.0) );
                phi_divfree.set(m,12,2, -(0.5*sq3*sq5)*(dy/length)*xi*(3.0*eta2-1.0) );

                phi_divfree.set(m,11,1,  (0.5*(sq5*sq7*dx*xi)/(l4*l6))
                        *((-21.0+105.0*eta2)*dy4+(-48.0+7.0*xi2+135.0*eta2)*dx2*dy2
                            +(9.0*eta2-3.0)*dx4) );
                phi_divfree.set(m,11,2, -(0.5*(sq5*sq7*dy*eta)/(l4*l6))
                        *((-21.0+35.0*eta2)*dy4+(-48.0+21.0*xi2+45.0*eta2)*dx2*dy2
                            +(-3.0+3.0*eta2)*dx4) );

                phi_divfree.set(m,10,1,  0.5*((sq3*sq7*dx)/(length*l4))*xi
                        *(dy2*(5.0*xi2-5.0)+dx2*(5.0*xi2-3.0)) );
                phi_divfree.set(m,10,2, -0.5*((sq3*sq7*dy)/(length*l4))*eta
                        *(dy2*(15.0*xi2-5.0)+dx2*(15.0*xi2-3.0)) );

            case 3:  // third order (4 basis vectors)
                phi_divfree.set(m,9,1,  0.5*sq5*(3.0*eta2-1.0) );
                phi_divfree.set(m,9,2,  0.0 );

                phi_divfree.set(m,8,1,  0.0 );
                phi_divfree.set(m,8,2, -0.5*sq5*(3.0*xi2-1.0) );

                phi_divfree.set(m,7,1,  3.0*sq5*xi*eta*(dx/lx5) );
                phi_divfree.set(m,7,2, -sq5/2.0*(3.0*eta2-1.0)*(dy/lx5) );

                phi_divfree.set(m,6,1,  sq5/2.0*(3.0*xi2-1.0)*(dx/ly5) );
                phi_divfree.set(m,6,2, -3.0*sq5*xi*eta*(dy/ly5) );

            case 2:  // second order (3 basis vectors)
                phi_divfree.set(m,5,1,  sq3*eta );
                phi_divfree.set(m,5,2,  0.0 );

                phi_divfree.set(m,4,1,  0.0 );
                phi_divfree.set(m,4,2, -sq3*xi );

                phi_divfree.set(m,3,1,  dx/length * sq3*xi  );
                phi_divfree.set(m,3,2, -dy/length * sq3*eta );

            case 1:  // first order (2 basis vectors)
                phi_divfree.set(m,2,1,  1.0 );
                phi_divfree.set(m,2,2,  0.0 );

                phi_divfree.set(m,1,1,  0.0 );
                phi_divfree.set(m,1,2, -1.0 );

                break;
        }
    }
}


// -------------------------------------------------------------------------- //
// Set values of the Legendre polys on each quadrature point
//
// TODO - what is this function used for?  (6/9/2014) -DS
// -------------------------------------------------------------------------- //
void SetLegendre(const int dummy_space_order, const dTensor2& spts, 
    dTensor2& phi)
{
    const int space_order = dogParams.get_space_order();
    assert_eq(dummy_space_order, space_order);
    SetLegendreAtPoints(spts, phi);
}

// -------------------------------------------------------------------------- //
// Constructor.
//
// TODO - document this function.  Does this allocate memory for every private
// variable?
// -------------------------------------------------------------------------- //
Legendre2d::Legendre2d():
    edgeData(0),
    quadratureWeights(0),
    quadraturePoints(0),
    phi(0),
    wgt_phi_transpose(0),
    phi_min_ptr(0),
    phi_max_ptr(0),
    phi_interval(0),
    quadratureWeightsCoarse(0), quadraturePointsCoarse(0), phiCoarse(0),
    phi_x(0), phi_y(0),
    numPositivityPoints(0), positivityPoints(0), phiAtPositivityPoints(0),
    numRiemannPoints(0), riemannPoints(0), phiAtRiemannPoints(0)    
{

    // Parameters defined by the problem.  
    // (TODO - we might as well pass these
    // into this constructor rather than rely on global variables.  What if we
    // want to persue p-refinement? -DS).
    const int space_order = dogParams.get_space_order();
    const int kmax        = dogParams.get_kmax();

    // --------------------------------------------------------------------- //
    // 1. Set quadrature rules necessary for integration accurate up to order 
    // space_order.  Additionally, compute the product of wgt*phi for efficient
    // integration.
    // --------------------------------------------------------------------- //
    const int mpoints1d = space_order;
    const int mpoints   = mpoints1d*mpoints1d;
    quadratureWeights   = new dTensor1(mpoints      );
    quadraturePoints    = new dTensor2(mpoints, 2   );
    phi                 = new dTensor2(mpoints, kmax);
    SetQuadrature(mpoints1d, *quadratureWeights, *quadraturePoints);
    SetLegendre(space_order, *quadraturePoints, *phi);

    // Compute the transpose of the product of wgts * phi.  This allow for 
    // greater efficiency when computing integrals because this order leverages
    // the order defined in the tensor class.
    wgt_phi_transpose   = new dTensor2(kmax, mpoints);
    for(int mp=1;mp<=mpoints;mp++)
    for(int k=1;k<=kmax;k++)
    {
        wgt_phi_transpose->set(k, mp, quadratureWeights->get(mp)*phi->get(mp,k));
    }

    // --------------------------------------------------------------------- //
    // 2. Set quadrature for integration accurate to order space_order-1. This
    // is used when projecting onto the gradient of the basis functions.
    // --------------------------------------------------------------------- //
    const int mpoints1dCoarse = space_order-1;
    const int mpointsCoarse   = mpoints1dCoarse*mpoints1dCoarse;
    quadratureWeightsCoarse   = new dTensor1( mpointsCoarse         );
    quadraturePointsCoarse    = new dTensor2( mpointsCoarse, 2      );
    phiCoarse                 = new dTensor2( mpointsCoarse, kmax   );
    SetQuadrature( mpoints1dCoarse, *quadratureWeightsCoarse, *quadraturePointsCoarse);
    SetLegendre( space_order, *quadraturePointsCoarse, *phiCoarse);

    // Grid information (This assumes constant grid size, but is needed to
    //                   compute derivatives ... )
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    phi_x           = new dTensor2(mpointsCoarse,kmax);
    phi_y           = new dTensor2(mpointsCoarse,kmax);
    SetLegendreGrad(dx, dy, *quadraturePointsCoarse, *phi_x, *phi_y);

    // --------------------------------------------------------------------- //
    // 3. Construct the edge data.  (1D quadrature points along each edge.)
    // This is required for the boundary integral in the weak formulation.
    // --------------------------------------------------------------------- //
    edgeData = new edge_data();
    edgeData->init();

    // --------------------------------------------------------------------- //
    // 4. Positivity points.  Additionally, we compute the min/max of each 
    // legendre polynomial in an L^\infty norm as an estimator of whether or not
    // to investigate whether or not to apply the limiter.
    // --------------------------------------------------------------------- //
    SetLegendreIntervals();
    setPositivityPoints();

    // --------------------------------------------------------------------- //
    // 5. ??? 
    // --------------------------------------------------------------------- //
    setRiemannPoints();

    // --------------------------------------------------------------------- //
    // 6. ???.  Do we need this?  (5/14/2014) -DS.
    // Check that positivity points concur with edgeData
    // --------------------------------------------------------------------- //
    check_edgeData();

};


// Destructor
Legendre2d::~Legendre2d()
{

    // --------------------------------------------------------------------- //
    // This memory is allocated in the constructor.
    // --------------------------------------------------------------------- //

    delete quadratureWeights;
    delete quadraturePoints;
    delete phi;
    delete wgt_phi_transpose;

    // Coarse points, weights, and legendre polynomials at these points.
    delete quadratureWeightsCoarse;
    delete quadraturePointsCoarse;
    delete phiCoarse;

    // Gradient of basis functions
    delete phi_x;
    delete phi_y;

    // Edge data
    delete edgeData;

    // --------------------------------------------------------------------- //
    // This memory is allocated in SetLegendreIntervals.
    // --------------------------------------------------------------------- //
    delete phi_min_ptr;
    delete phi_max_ptr;
    delete phi_interval;

    // --------------------------------------------------------------------- //
    // These are constructed in setPositivityPoints.  -DS (6/11/2014).
    // --------------------------------------------------------------------- //
    delete phiAtPositivityPoints;
    delete positivityPoints;

    // --------------------------------------------------------------------- //
    // This memory is allocated in setRiemannPoints. As of (6/11/2014) was
    // not deallocated.
    // --------------------------------------------------------------------- //
    delete riemannPoints;
    delete phiAtRiemannPoints;

}

// -------------------------------------------------------------------------- //
// setPositivityPoints.  Set the positivity points.
//
// This routine computes the union of of two tensor products that are defined by
// Gauss-Lobatto as well as Gauss-Legendre quadrature points.  Specifically, we
// compute
//
//     ( x1d_lobatto, x1d_gauss ) \cup ( x1d_gauss, x1d_lobatto ),
//
// where x1d_lobatto are the 1D Lobatto points (that include the endpoints), 
// and x1d_gauss are the Gaussian quadrature points (where the Riemann problems
// are solved).
//
// Because it is not necessary to use the Gauss-Lobatto points in any
// integrals we do not compute the weights.  See Theorem 2.1 in
// [ZhangShu10]:
//
// [ZhangShu10] X. Zhang, C-W. Shu, On positivity preserving high
// order discontinuous Galerkin schemes for compressible Euler
// equations on rectangular meshes, Journal of Computational
// Physics (2010), doi: 10.1016/j.jcp.2010.08.016
//
// See also: ...
// -------------------------------------------------------------------------- //
void Legendre2d::setPositivityPoints()
{

    const int polynomial_order          = dogParams.get_space_order()-1;
    const int mGaussLobattoPoints1d     = int(ceil((polynomial_order+3)/2.));
    const int mGaussPoints1d            = dogParams.get_space_order();

    // Grab the 1D points
    dTensor1 gaussLobattoPoints1d   ( mGaussLobattoPoints1d  );
    dTensor1 gaussPoints1d          ( mGaussPoints1d         );
    setGaussLobattoPoints1d         ( gaussLobattoPoints1d   );
    setGaussPoints1d                ( gaussPoints1d          );

    const int mpoints = mGaussLobattoPoints1d*mGaussPoints1d;

    // We will enforce positivity on a total of two meshes.
    numPositivityPoints = 2*mpoints;
    positivityPoints    = new dTensor2(numPositivityPoints,2);

    // -------------------------------------------------------- //
    // 1. Construct points of the form ( x1d_gauss, x1d_lobatto ).
    // -------------------------------------------------------- //
    int k1=0;
    for(int m1=1;m1<=mGaussPoints1d;m1++)
    for(int m2=1;m2<=mGaussLobattoPoints1d;m2++)
    {
        k1++;
        positivityPoints->set(k1,1,gaussPoints1d.get(m1));
        positivityPoints->set(k1,2,gaussLobattoPoints1d.get(m2));
    }

    // -------------------------------------------------------- //
    // 2. Construct points of the form ( x1d_lobatto, x1d_gauss ).
    // -------------------------------------------------------- //
    assert_eq(k1,mpoints);
    int k2=mpoints;
    for(int m1=1;m1<=mGaussPoints1d;m1++)
    for(int m2=1;m2<=mGaussLobattoPoints1d;m2++)
    {
        k2++;
        positivityPoints->set(k2,1,gaussLobattoPoints1d.get(m2));
        positivityPoints->set(k2,2,gaussPoints1d.get(m1));
    }
    assert_eq(k2,2*mpoints);

    // -------------------------------------------------------- //
    // 3. Evaluate the basis functions at the positivity points
    // -------------------------------------------------------- //
    phiAtPositivityPoints = new dTensor2(
        positivityPoints->getsize(1), dogParams.get_kmax());
    SetLegendreAtPoints(*positivityPoints, *phiAtPositivityPoints);

    // Finsished with this routine.

    // -------------------------------------------------------- //
    // 4. Debugging options
    // -------------------------------------------------------- //

    // print out coordinates and phi values at Riemann points
    // (compare with setRiemannPoints() and edge_data::init()).
    //
    // TODO - how do we insert these debug options? Do we always want this
    // turned on? (-DS)
    //
    if(debug3)
    {
        const int kmax = dogParams.get_kmax();
        int k1=0;
        printf("\n points at y boundaries:");
        for( int m1=1; m1 <= mGaussPoints1d;        m1++)
        for( int m2=1; m2 <= mGaussLobattoPoints1d; m2++)
        {
            k1++;
            if(m2==1||m2==mGaussLobattoPoints1d)
            {
                printf("\n   phi(%d,*) at (%24.16e, %24.16e ):",k1,
                        positivityPoints->get(k1,1),
                        positivityPoints->get(k1,2));
                for(int k=1;k<=kmax;k++)
                {
                    printf("\n    phi(%d,%d) = %24.16e",
                            k1,k,phiAtPositivityPoints->get(k1,k));
                }
                printf("\n");
            }
        }
        assert_eq(k1,mpoints);
        int k2=mpoints;
        printf("\n points at x boundaries:");
        for( int m1=1; m1 <= mGaussPoints1d;        m1++)
        for( int m2=1; m2 <= mGaussLobattoPoints1d; m2++)
        {
            k2++;
            if( m2==1 || m2 == mGaussLobattoPoints1d )
            {
                printf("\n   phi(%d,*) at (%24.16e, %24.16e ):",k2,
                        positivityPoints->get(k2,1),
                        positivityPoints->get(k2,2));
                for(int k=1;k<=kmax;k++)
                {
                    printf("\n    phi(%d,%d) = %24.16e",
                            k2,k,phiAtPositivityPoints->get(k2,k));
                }
                printf("\n");
            }
        }
        assert_eq(k2,2*mpoints);
    }

}

// -------------------------------------------------------------------------- //
// TODO - document this "function"
// If we have edge_data, why do we need this routine? -DS
// -------------------------------------------------------------------------- //
void Legendre2d::setRiemannPoints()
{

    const int mGaussPoints1d = dogParams.get_space_order();

    // Construct the (1D) Gaussian quadrature points
    dTensor1 gaussPoints1d(mGaussPoints1d);
    setGaussPoints1d(gaussPoints1d);

    numRiemannPoints = mGaussPoints1d*4;

    // TODO - where does this get deallocated? -DS (6/9/2014)
    riemannPoints = new dTensor2(numRiemannPoints, 2);

    int k1=0;
    for(int m1=1;m1<=mGaussPoints1d;m1++) // yr
    {
        k1++;
        riemannPoints->set(k1,1, gaussPoints1d.get(m1));
        riemannPoints->set(k1,2, -1.);
    }
    for(int m1=1;m1<=mGaussPoints1d;m1++) // yl
    {
        k1++;
        riemannPoints->set(k1,1, gaussPoints1d.get(m1));
        riemannPoints->set(k1,2, 1.);
    }
    for(int m1=1;m1<=mGaussPoints1d;m1++) // xr
    {
        k1++;
        riemannPoints->set(k1,1, -1.);
        riemannPoints->set(k1,2, gaussPoints1d.get(m1));
    }
    for(int m1=1;m1<=mGaussPoints1d;m1++) // xl
    {
        k1++;
        riemannPoints->set(k1,1, 1.);
        riemannPoints->set(k1,2, gaussPoints1d.get(m1));
    }
    assert_eq(k1, numRiemannPoints);

    // TODO - where does this get deallocated? -DS (6/9/2014)
    phiAtRiemannPoints = new dTensor2(riemannPoints->getsize(1),
            dogParams.get_kmax());
    SetLegendreAtPoints(*riemannPoints, *phiAtRiemannPoints);

    // print out coordinates and phi values at Riemann points
    // ( compare with setPositivityPoints() )
    if(debug3)
    {
        const int kmax = dogParams.get_kmax();
        for(int m1=1;m1<=numRiemannPoints;m1++)
        {
            printf("\n   phi(%d,*) at (%24.16e, %24.16e ):",m1,
                    riemannPoints->get(m1,1),
                    riemannPoints->get(m1,2));
            for(int k=1;k<=kmax;k++)
            {
                printf("\n    phi(%d,%d) = %24.16e",
                        m1,k,phiAtRiemannPoints->get(m1,k));
            }
            printf("\n");
        }
    }
}

// -------------------------------------------------------------------------- //
// For use after applying limiters
//
// should also create a linear projection
// project_linearly_onto_locally_divergence_free_subspace
// which uses averages instead of minmod.
// This could be used to project the electric field onto the
// affine subspace of functions whose divergence agrees with the
// charge density over permittivity (with highest coefficient
// zeroed).
//
// (moved from ApplyLimiter_divfree.cpp)
// 
// Input:
// ------
//
//     TODO
//
// Returns:
// --------
//
//     TODO
//
// -------------------------------------------------------------------------- //
void project_onto_locally_divergence_free_subspace(dTensorBC4& q)
{
    const int mx = q.getsize(1);
    const int my = q.getsize(2);
    const int mbc  = q.getmbc();
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    const double length = sqrt(dx*dx + dy*dy);
    const double lx5    = sqrt(5.0*dx*dx + dy*dy);
    const double ly5    = sqrt(dx*dx + 5.0*dy*dy);
    const int how_many = dogParams.get_how_many_vectors_divfree();

    //#pragma omp parallel for
    for (int k=1; k<=how_many; k++)
    {
        int mcomp = dogParams.get_which_compnt_divfree()[k];

        switch ( dogParams.get_space_order() )
        {
            default:
                unsupported_value_error(dogParams.get_space_order());
            case 1:
                break;
            case 2:  // 2nd order in space

#pragma omp parallel for
                for (int i=(3-mbc); i<=(mx+mbc-2); i++)   
                    for (int j=(3-mbc); j<=(my+mbc-2); j++)
                    {
                        const double B3 = minmod(q.get(i,j,mcomp,2)/dx,
                                -q.get(i,j,mcomp+1,3)/dy);

                        q.set(i,j,mcomp,  2,  dx*B3 );
                        q.set(i,j,mcomp+1,3, -dy*B3 );
                    }

                break;

            case 3:  // 3rd order in space

#pragma omp parallel for
                for (int i=(3-mbc); i<=(mx+mbc-2); i++)   
                    for (int j=(3-mbc); j<=(my+mbc-2); j++)
                    {
                        double B3 = minmod(q.get(i,j,mcomp,2)/dx,       -q.get(i,j,mcomp+1,3)/dy);
                        double B6 = minmod(q.get(i,j,mcomp,5)/dx,       -q.get(i,j,mcomp+1,4)/(sq5*dy));
                        double B7 = minmod(q.get(i,j,mcomp,4)/(sq5*dx), -q.get(i,j,mcomp+1,6)/dy);

                        q.set(i,j,mcomp,  2,      dx*B3 );
                        q.set(i,j,mcomp+1,3,     -dy*B3 );

                        q.set(i,j,mcomp,  5,      dx*B6 );
                        q.set(i,j,mcomp+1,4, -sq5*dy*B6 );

                        q.set(i,j,mcomp,  4,  sq5*dx*B7 );
                        q.set(i,j,mcomp+1,6,     -dy*B7 );          
                    }

                break;
        }
    }
}


// Debugging routine for checking edgeData.
void Legendre2d::check_edgeData() const
{

    const int mGaussPoints1d    = dogParams.get_space_order();
    const int kmax              = dogParams.get_kmax();

    const dTensor2& phi         = *phiAtRiemannPoints;


    int k1=0;
    for(int m1=1; m1 <= mGaussPoints1d; m1++) // yr
    {   
        k1++;
        for(int k=1;k<=kmax;k++)
            assert_eq(phi.get(k1,k), edgeData->phi_yr->get(m1,k));
    }

    for(int m1=1; m1 <= mGaussPoints1d; m1++) // yl
    { 
        k1++;
        for(int k=1;k<=kmax;k++)
            assert_eq(phi.get(k1,k), edgeData->phi_yl->get(m1,k));
    }
    for(int m1=1; m1 <= mGaussPoints1d; m1++) // xr
    {   
        k1++;
        for(int k=1;k<=kmax;k++)
            assert_eq(phi.get(k1,k), edgeData->phi_xr->get(m1,k));
    }
    for(int m1=1; m1 <= mGaussPoints1d; m1++) // xl
    { 
        k1++;
        for(int k=1;k<=kmax;k++)
            assert_eq(phi.get(k1,k), edgeData->phi_xl->get(m1,k));
    }
    assert_eq(k1,numRiemannPoints);

}

Legendre2d& Legendre2d::instance() // const
{

    // TODO - what's the correct way to dealloate this memory? (5/11/2014) -DS
    static Legendre2d* legendre2d = new Legendre2d;
    return *legendre2d;

}
