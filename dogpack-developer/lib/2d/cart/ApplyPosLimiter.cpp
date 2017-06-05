#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "Legendre2d.h"
#include "ApplyPosLimiter.h"

// Routine that applies the Zhang and Shu limiter to higher order moments of
// the conserved variables.  This limiter is locally applied, and conserves
// total mass by not adjusting cell averages.
//
// See: X. Zhang and C.-W. Shu, "Maximum-principle-satisfying and
// positivity-preserving high order schemes for conservation laws: Survey and
// new developments", Proceedings of the Royal Society A, v467 (2011),
// pp.2752-2776. for more details on the limiter
void ApplyPosLimiter(const dTensorBC4& aux, dTensorBC4& q)
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
    // const int mpoints   =  space_order*space_order;
    const int mpoints   =  4+4*space_order+space_order*space_order;
    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //

    // Quadrature points
    dTensor2 spts(mpoints,2);
    SetPositivePoints(space_order,spts);

    // Basis functions evaluated at each point
    dTensor2 phi(mpoints,kmax);
    SetLegendreAtPoints(space_order,spts,phi);

    // -------------------------------------------------------------- //
    //                                                                //
    // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
    // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
    //                                                                //
    // -------------------------------------------------------------- //
#pragma omp parallel for
    for(int i=1-mbc; i <= mx+mbc; i++)
    for(int j=1-mbc; j <= my+mbc; j++)
    for(int me=1; me <= meqn; me++)
    {
        double m = 0.0;
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double qnow = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                qnow += q.get(i,j,me,k) * phi.get(mp,k);
            }
            m = Min(m, qnow);
        }

        double theta = 0.0;
        double Q1 = q.get(i,j,me,1);
        if( fabs( Q1 - m ) < 1.0e-14 ){ theta = 1.0; }
        else{ theta = Min( 1.0, fabs( Q1 / (Q1 - m) ) ); }

        // limit q //
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,j,me,k, q.get(i,j,me,k) * theta );
        }

    }

}

// TODO - we can likely remove this routine because these (used to be) the
// Gaussian quadrature points.  However, an application will likely want to add in the edge
// points that are used for Riemann solves.
void SetPositivePoints(const int& space_order, dTensor2& spts)
{

    // 1D Gaussian quadrature points
    dTensor1 s1d(space_order);
    void setGaussPoints1d(dTensor1& x1d);
    setGaussPoints1d( s1d );

    // This region has been commented out because we are no longer applying
    // the limiter at the corner points

    // 2D points -- corners
//  spts.set(1,1, -1.0e0 );
//  spts.set(1,2, -1.0e0 );

//  spts.set(2,1,  1.0e0 );
//  spts.set(2,2, -1.0e0 );

//  spts.set(3,1, -1.0e0 );
//  spts.set(3,2,  1.0e0 );

//  spts.set(4,1,  1.0e0 );
//  spts.set(4,2,  1.0e0 );

    // 2D points -- left, right, bottom and top edges
//  for (int m=1; m<=space_order; m++)
//  {
//      double s = s1d.get(m);

//      // left edge
//      spts.set(4+m,1, -1.0e0 );
//      spts.set(4+m,2,  s     );

//      // right edge
//      spts.set(4+space_order+m,1,  1.0e0 );
//      spts.set(4+space_order+m,2,  s     );

//      // bottom edge
//      spts.set(4+2*space_order+m,1,  s     );
//      spts.set(4+2*space_order+m,2, -1.0e0 );

//      // top edge
//      spts.set(4+3*space_order+m,1,  s     );
//      spts.set(4+3*space_order+m,2,  1.0e0 );
//  }

    // 2D points -- all interior points
//  int z = 4+4*space_order;
    int z = 0;
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

// This routine is used in SL libraries
void SamplePhiAtPositivePoints(const int& space_order, 
        const dTensor2& spts, dTensor2& phi)
{
    const int mpoints = spts.getsize(1);

    for (int m=1; m<=mpoints; m++)
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
