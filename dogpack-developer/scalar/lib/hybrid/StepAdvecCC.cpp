#include <math.h>
#include <iostream>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"

///////////////////////////////////////////////////////////////////////////////
// Function for stepping advection equation:
//
//            q_t + u q_x + v q_y = 0 
//
// Parameters
//
// dt - time step taken in each equation
// qold(mx,my,meqn,kmax) - q before translation
// qnew - q after being translated and projected onto the legendre basis
// speeds(1:meqn, 2 );  advection velocities for each equation.
// 
///////////////////////////////////////////////////////////////////////////////
void StepAdvecCC(double dt, 
    dTensorBC4& qold,            // Positiivty limiter affects qold
    const dTensor2& speeds,
    dTensorBC4& qnew )
{

    //-local parameters -----------------------------------------------------//
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);  assert_eq( meqn, speeds.getsize(1) );
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();
    const int space_order = dogParams.get_space_order();

    // number of points used for integration on each cell
    const int numpts = 4*space_order*space_order;

    // cell widths are uniform
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();
    //-----------------------------------------------------------------------//

    const double s_area = 4.0;

    // Fractional eta values - 
    //      frac_eta( 1:meqn, 1 ) = x-direction, 
    //      frac_eta( 1:meqn, 2 ) = y direction
    dTensor2 frac_eta  ( meqn, 2 );
    dTensor2 large_eta ( meqn, 2 );
    for( int m=1; m <= meqn; m++ )
    {

        large_eta.set ( m, 1, speeds.get(m,1) * dt / dx );
        frac_eta.set  ( m, 1, (double)( large_eta.get(m,1) - floor( large_eta.get(m,1) ) ) );

        large_eta.set ( m, 2, speeds.get(m,2) * dt / dy );
        frac_eta.set  ( m, 2, (double)( large_eta.get(m,2) - floor( large_eta.get(m,2) ) ) );
    }

    // integration points and weights
    // Because each 'equation' passed into this routine has its own speed, we
    // save a list of all the necessary weights for each equation.
    dTensor2  wgt(meqn, numpts);
    dTensor3 spts(meqn, numpts, 2);
    void setIntegrationPoints(const dTensor2& frac_eta, dTensor3& spts, dTensor2& wgt);
    setIntegrationPoints(frac_eta, spts, wgt);

    // 'old' integration points:
    dTensor3 spts_old ( meqn, numpts, 2 );
    iTensor2 ishift   ( meqn, numpts    );
    iTensor2 jshift   ( meqn, numpts    );
    void translateXi2D( 
        const dTensor2& large_eta, const dTensor3& xi, 
        dTensor3& xinew, iTensor2& ishift, iTensor2& jshift );
    translateXi2D( large_eta, spts, spts_old, ishift, jshift );

    // Legendre Basis functions, evaluated at all the
    // necessary points
    //
    dTensor3 phi    (meqn, numpts, kmax);
    dTensor3 phi_old(meqn, numpts, kmax);
    
    // Stolen and copied from lib/2d/cart/Legendre2d:
    // this version allows for meqn to be passed in ...
    void evaluateLegendrePolys( const dTensor3& spts, dTensor3& phi );
    evaluateLegendrePolys( spts,     phi      );
    evaluateLegendrePolys( spts_old, phi_old  );

    // -------------------------------------------------------------- //
    // Positivity preserving limiter                                  //
    // -------------------------------------------------------------- //
    if( dogParams.using_moment_limiter() )
    {

        for(int i=1; i<=mx; i++)
        for(int j=1; j<=my; j++)
        for(int me=1; me <= meqn; me++ )
        {

            // -------------------------------------------------------------- //
            // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
            // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
            // -------------------------------------------------------------- //
            double m = qold.get(i,j,me,1);
            for(int mp=1; mp <= phi_old.getsize(2); mp++)
            {
                // evaluate q at spts(mp) //
                double qnow = 0.0;
                for( int k=1; k <= kmax; k++ )
                {
                    qnow += qold.get(i, j, me, k) * phi_old.get(me, mp, k);
                }
                m = Min(m, qnow );
            }

            double Q1 = qold.get(i,j,me,1);

            double theta = 0.0;
            if( fabs( Q1 - m ) < EPSILON ) { theta = 0.0; }
            else{ theta = Min( 1.0, (Q1 - EPSILON) / (Q1 - m) ); }

            // limit q //
            for( int k=2; k <= kmax; k++ )
            {
                qold.set(i,j,me,k, qold.get(i,j,me,k) * theta );
            }

        }

    }
    // End of limiter //
    // //////////////////////////////////////////// //

    for(int i=1; i<=mx; i++)
    for(int j=1; j<=my; j++)
    for(int me=1; me <= meqn; me++ )
    {

        for( int k=1; k <= kmax; k++)
        {
            double sum = 0.0;
            for( int m=1; m <= numpts; m++ )
            {

                // hard coded periodicity is enforced here!
                int io = (int)( i + ishift.get(me,m) );
                io = iMod((io-1), mx)+1;

                int jo = (int)( j + jshift.get(me,m) );
                jo = iMod((jo-1), my)+1;

                // evaluate q at the old point:
                double qval = 0.0;
                for( int kq=1; kq <= kmax; kq++ )
                {
                    qval += qold.get(io, jo, me, kq) * phi_old.get(me, m, kq);
                }
                sum += wgt.get(me, m) * qval * phi.get(me, m, k);
            }
            qnew.set(i, j, me, k, sum / s_area );
        }

    }

//  void SetBndValues(const dTensor2& node, dTensorBC3& aux, dTensorBC3& q);
//  SetBndValues(node, auxvals, qnew);

}

void setIntegrationPoints(const dTensor2& frac_eta, dTensor3& spts, dTensor2& wgt)
{

    const int morder = dogParams.get_space_order();

    const int meqn   = wgt.getsize(1);
    const int numpts = wgt.getsize(2);

    //////////////////////////////////////////////////////////////////////////
    // Set quadrature weights and points (standard before a transformation)
    //////////////////////////////////////////////////////////////////////////
    dTensor1 w1d(morder), x1d(morder);
    switch ( morder )
    {
        case 1:
            w1d.set(1,  2.0e0 );
            x1d.set(1, 0.0e0 );

            break;

        case 2:
            w1d.set(1,   1.0 );
            w1d.set(2,   1.0 );

            x1d.set(1, -1.0/sq3 );
            x1d.set(2,  1.0/sq3 );

            break;

        case 3:
            w1d.set(1, 5.0e0/9.0e0 );
            w1d.set(2, 8.0e0/9.0e0 );
            w1d.set(3, 5.0e0/9.0e0 );

            x1d.set(1, -sq3/sq5 );
            x1d.set(2,  0.0e0 );
            x1d.set(3,  sq3/sq5 );

            break;

        case 4:
            w1d.set(1, (18.0 - sqrt(30.0))/36.0 );
            w1d.set(2, (18.0 + sqrt(30.0))/36.0 );
            w1d.set(3, (18.0 + sqrt(30.0))/36.0 );
            w1d.set(4, (18.0 - sqrt(30.0))/36.0 );

            x1d.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
            x1d.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
            x1d.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
            x1d.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );

            break;

        case 5:
            w1d.set(1, (322.0 - 13.0*sqrt(70.0))/900.0 );
            w1d.set(2, (322.0 + 13.0*sqrt(70.0))/900.0 );
            w1d.set(3,  128.0/225.0 );
            w1d.set(4, (322.0 + 13.0*sqrt(70.0))/900.0 );
            w1d.set(5, (322.0 - 13.0*sqrt(70.0))/900.0 );

            x1d.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
            x1d.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            x1d.set(3,  0.0 );
            x1d.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            x1d.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );

            break;
    }

    void setIntegrationPoints1D(
        const dTensor1& w1d,
        const dTensor1& x1d,
        double frac_eta, dTensor1& spts, dTensor1& wgt);

    // quadrature weights and points for integration in the x-direction:
    dTensor1 wx( 2*morder );
    dTensor1 sx( 2*morder );

    // quadrature weights and points for integration in the y-direction:
    dTensor1 wy( 2*morder );
    dTensor1 sy( 2*morder );

    for( int me=1; me <= meqn; me++ )
    {

        // Set the 1D integration points
        setIntegrationPoints1D(w1d, x1d, frac_eta.get(me,1), sx, wx);
        setIntegrationPoints1D(w1d, x1d, frac_eta.get(me,2), sy, wy);

        // Tensor product on all of the weights and points:
        int k=0;
        for( int i=1; i <= 2*morder; i++ )
        for( int j=1; j <= 2*morder; j++ )
        {
            k++;
            wgt.set (me, k,    wx.get(i)*wy.get( j )   );
            spts.set(me, k, 1, sx.get(i)               );
            spts.set(me, k, 2, sy.get(j)               );
        }

    }

    // check that the weights sum to 4:
//  double tmp = 0.;
//  for( int n=1; n <= wgt.getsize(2); n++ )
//  {
//      tmp += wgt.get(1,n);
//  }
//  printf("tmp = %f\n", tmp );
//  printf("frac_eta(1) = %f\n", frac_eta.get(1,1) );
//  printf("frac_eta(2) = %f\n", frac_eta.get(1,2) );


}

// Integration points for a 1D problem, where we take into account the exact
// location of the discontinuity
void setIntegrationPoints1D( const dTensor1& w1d, const dTensor1& x1d,
    double frac_eta, dTensor1& spts, dTensor1& wgt)
{

    const int mpts   = wgt.getsize();               
    const int morder = dogParams.get_space_order();
    assert_eq( mpts, 2*morder);

    // location of discontinuity and width of left and right hand side
    const double scut      = -1.0 + 2.0*frac_eta;
    const double left_len  = scut - (-1.0);
    const double right_len = 1.0  - scut;

    // left half of the integral
    const double bl = scut;
    const double al = -1.0;

    // right half of the integral
    const double br = 1.0;
    const double ar = scut;

    // modify weights and find modified quadrature points to evaluate q
    for( int m=1; m <= morder; m++ )
    {

        // weights
        wgt.set(m,           left_len / 2.0 * w1d.get(m) );
        wgt.set(m + morder, right_len / 2.0 * w1d.get(m) );

        // points
        double x = x1d.get(m);
        spts.set(m,         0.5*( (bl-al)*x + (al+bl) ) );
        spts.set(morder+m,  0.5*( (br-ar)*x + (ar+br) ) );
    }

}

// Function that traces a quadrature point back in time to its old value, and
// computes the index of its shift.  It is expected that xi \in [-1,1]
//
//    xi    ( 1:meqn, 1:numpts (==4*M^2), 1:2 )
//    xinew ( 1:meqn, 1:numpts (==4*M^2), 1:2 )
//  ishift( 1:meqn, 1:numpts )
//  jshift( 1:meqn, 1:numpts )
void translateXi2D( 
        const dTensor2& large_eta, const dTensor3& xi, 
        dTensor3& xinew, iTensor2& ishift, iTensor2& jshift )
{

    const int meqn = large_eta.getsize(1);
    assert_eq( meqn, ishift.getsize(1) );
    assert_eq( meqn, jshift.getsize(1) );
    assert_eq( meqn,     xi.getsize(1) );

    const int numpts = xi.getsize(2);
    assert_eq( numpts, ishift.getsize(2) );
    assert_eq( numpts, jshift.getsize(2) );

    assert_eq( xi.getsize(3), xinew.getsize(3) );
    assert_eq( xi.getsize(3), 2 );

    for( int me=1; me <= meqn; me++ )
    {
        double eta_x = large_eta.get(me,1);
        double eta_y = large_eta.get(me,2);

        for( int m=1; m <= numpts; m++ )
        {
            double xi_old = xi.get(me,m,1) - 2.0 * eta_x;
            double nu_old = xi.get(me,m,2) - 2.0 * eta_y;
            ishift.set( me, m, (int) floor( 0.5 * (xi_old+1.0) ) );
            jshift.set( me, m, (int) floor( 0.5 * (nu_old+1.0) ) );
            xinew.set(  me, m, 1, xi_old - 2.0 * ( ishift.get(me,m) ) );
            xinew.set(  me, m, 2, nu_old - 2.0 * ( jshift.get(me,m) ) );
        }
    }

}

void evaluateLegendrePolys(const dTensor3& spts, dTensor3& phi)
{

    const int meqn  = phi.getsize(1);
    assert_eq( meqn, spts.getsize(1) );

    const int space_order = dogParams.get_space_order();
    const int kmax = phi.getsize(3);

    // This assertion is only valid is kmax = kmax_cart ...
    // assert_eq(kmax,dogParams.get_kmax());
    switch( space_order )
    {

        case 3:
            assert_eq( kmax, 6 );
            break;

        case 2:
            assert_eq( kmax, 3 );
            break;

        case 1:
            assert_eq( kmax, 1 );
            assert_eq( dogParams.get_kmax(), 1 );
            break;
        default:
            unsupported_value_error( space_order );
    }


    const int mpoints = spts.getsize(2);
    assert_eq(mpoints,phi.getsize(2));
    for(int me=1; me <= meqn; me++ )
    for(int m=1; m<=mpoints; m++)
    {
        // grid point (x,y)
        const double xi  = spts.get(me, m,1);
        const double eta = spts.get(me, m,2);
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
                phi.set(me, m,15, 105.0/8.0*eta4 - 45.0/4.0*eta2 + 9.0/8.0 );
                phi.set(me, m,14, 105.0/8.0*xi4  - 45.0/4.0*xi2  + 9.0/8.0 );
                phi.set(me, m,13, 5.0/4.0*(3.0*xi2 - 1.0)*(3.0*eta2 - 1.0) );
                phi.set(me, m,12, sq3*sq7*(2.5*eta3 - 1.5*eta)*xi );
                phi.set(me, m,11, sq3*sq7*(2.5*xi3 - 1.5*xi)*eta );

            case 10:  // fourth order
                phi.set(me, m,10, sq7*(2.5*eta3 - 1.5*eta) );
                phi.set(me, m,9,  sq7*(2.5*xi3 - 1.5*xi) );
                phi.set(me, m,8,  sq3*sq5*xi*(1.5*eta2 - 0.5) );
                phi.set(me, m,7,  sq3*sq5*eta*(1.5*xi2 - 0.5) );

            case 6:  // third order
                phi.set(me, m,6,  sq5*(1.5*eta2 - 0.5) );
                phi.set(me, m,5,  sq5*(1.5*xi2 - 0.5) );
                phi.set(me, m,4,  3.0*xi*eta );                  

            case 3:  // second order                
                phi.set(me, m,3, sq3*eta );
                phi.set(me, m,2, sq3*xi  );

            case 1:  // first order
                phi.set(me, m,1, 1.0 );

                break;                
            default:
                unsupported_value_error(space_order);
        }

    }
}


