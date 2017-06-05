#include<cmath>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"

#include "LegendrePolys1d.h"

void ApplyPosLimiter(const dTensorBC3& aux, dTensorBC3& q)
{
    double Min(double, double);
    const int   mx = q.getsize(1);
    const int meqn = q.getsize(2);
    const int kmax = q.getsize(3);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(2);
    const int space_order = dogParams.get_space_order();

    // ------------------------------------------------ //
    // number of points where we want to check solution //
    // ------------------------------------------------ //
    const int mpoints =  kmax;

    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //
    const int MAX_KMAX = 6;
    dTensor1 spts(mpoints);
    dTensor2 phi(mpoints,MAX_KMAX);

    // Use Legendre Points for plotting purposes //
    setGaussPoints1d( spts );
    evaluateLegendrePolys( spts, phi );

    // -------------------------------------------------------------- //
    // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
    // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
    // -------------------------------------------------------------- //
#pragma omp parallel for
    for(int i=1-mbc; i <= mx+mbc; i++)
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
        const double Q1 = q.get(i,me,1);
        if( fabs( Q1 - m ) < 1.0e-14 ){ theta = 1.0; }
        else{ theta = Min( 1.0, fabs( Q1 / (Q1 - m) ) ); }

        // limit q //
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,me,k, q.get(i,me,k) * theta );
        }

    }

}
