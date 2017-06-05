#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "LegendrePolys1d.h"

// Apply positivity limiter.  1D Euler equations.
//
// This routine limits higher-order moments for Mass, momentum and Energy in
// order to gaurantee that pointwise values of the density and pressure remain
// positive.  This is constructed in a ``frozen-in-time'' setting, without
// regard to time discretizations, and is performed locally in each element.
//
// TODO - link to a useful reference on this.
//
// This is part I of the Zhang and Shu limiter for Euler equations (with a
// minor modification), and provably does not break the order of accuracy.
//
// See also: ApplyMPPLimiter.
void ApplyPosLimiter(const dTensorBC3& aux, dTensorBC3& q)
{

    double Min(double, double);
    const int   mx = q.getsize(1);
    const int meqn = q.getsize(2);
    const int kmax = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(2);
    const int space_order = dogParams.get_space_order();

    const double eps = 1.0e-14;
    

    // ------------------------------------------------ //
    // number of points where we want to check solution //
    // ------------------------------------------------ //
    const int mpoints =  kmax+2;

    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //
    const int MAX_KMAX = 6;
    dTensor1 spts(mpoints);
    dTensor2 phi(mpoints,MAX_KMAX);

    // Use Legendre Points for plotting purposes //
    setGaussPoints1dwBdry( spts );
    evaluateLegendrePolys( spts, phi );


    // -------------------------------------------------------------- //
    // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
    // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
    //
    //  We limit Density, Energy, and pressure only for Euler equations
    // -------------------------------------------------------------- //
#pragma omp parallel for
    for(int i=1-mbc; i <= mx+mbc; i++)
    {

        // ------------------------------------------ //
        // Limit the Density and Energy, first pass
        // ------------------------------------------ //
        double mrho  = 0.0;  // minimum point value (of density)
        double me    = 0.0;  // minimum point value (of energy)
        for(int mp=1; mp <= mpoints; mp++)
        {
            // evaluate q at spts(mp) //
            double rhonow = 0.0;
            double enow   = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                rhonow += q.get(i,1,k) * phi.get(mp,k);
                enow += q.get(i,5,k) * phi.get(mp,k);
            }
            mrho = Min(mrho, rhonow);
            me = Min(me, rhonow);
        }

        double thetarho = 0.0;
        double Q1    = q.get(i,1,1);
        if( fabs( Q1 - mrho ) < 1.0e-14 ){ thetarho = 0.0; }
        else{ thetarho = Min( 1.0, fabs( (eps-Q1) / (Q1 - mrho) ) ); }

        double thetae = 0.0;
        Q1    = q.get(i,5,1);
        if( fabs( Q1 - me ) < 1.0e-14 ){ thetae = 0.0; }
        else{ thetae = Min( 1.0, fabs( (eps-Q1) / (Q1 - me) ) ); }

        // TODO - This looks like we actually limit the solution twice! -DS
        //        Do we want one, or two values for theta?

        // limit q //
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,1,k, q.get(i,1,k) * thetarho );
            q.set(i,5,k, q.get(i,5,k) * thetae );
        }


        // ------------------------------------------------------------------ //
        //
        // Limit the momentum (to obtain positive pressure).
        //
        // This step works by considering an expansion of density, momentum
        // and energy in terms of a single parameter, theta.
        //
        // ------------------------------------------------------------------ //

        // Average momentum
        Q1 = q.get(i,2,1);
        double theta  = 1.0;
        double thetam = 1.0;
        double m      = 0.0;
        double minr   = 1.0;
        double maxr   = 1.0;
        double ave    = q.get(i,1,1)*q.get(i,5,1)-0.5*q.get(i,2,1)*q.get(i,2,1);

        // Average values for rho, momentum, and energy.
        double rhoa = q.get(i,1,1);
        double mxa  = q.get(i,2,1);
        double ea   = q.get(i,5,1);

        for(int mp=1; mp <= mpoints; mp++)
        {

            // evaluate q at spts(mp) //
            double mxnow  = 0.0;
            double enow   = 0.0;
            double rhonow = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
                rhonow  += q.get(i,1,k) * phi.get(mp,k);
                mxnow   += q.get(i,2,k) * phi.get(mp,k);
                enow    += q.get(i,5,k) * phi.get(mp,k);
            }

            double mxc = mxnow-mxa;

            double ec   = enow-ea;
            double rhoc = rhonow-rhoa;

            double slope=ec*rhoa+rhoc*ea-(mxc*mxa);
            double press=rhonow*enow-0.5*mxnow*mxnow;

            if( fabs( ave - press ) < 1.0e-14 ){ thetam = thetam; }
            else{ thetam = Min( thetam, fabs( (eps-ave) / slope ) ); }

            m = Min( m, press );

        }

        // ------------------------------------------------------------------ //
        // Final limiting step (limit Density, Momemtum, and Energy)
        // ------------------------------------------------------------------ //
        if( fabs( ave - m ) < eps ){ theta = Min(theta,1.0); }
        else{ theta = Min( thetam, fabs( (eps-ave) / (ave-m) ) );}
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,1,k, q.get(i,1,k) * theta );
            q.set(i,2,k, q.get(i,2,k) * theta );
            q.set(i,5,k, q.get(i,5,k) * theta );
        }

    }

}
