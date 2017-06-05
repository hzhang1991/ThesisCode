#include <cmath>
#include <iostream>
#include "DogParams.h"
#include "dogdefs.h"
#include "dog_math.h"
#include "Legendre2d.h"
using namespace std;

// Positivity preserving limiter: 2D Euler equations.
//
// Routine that applies the Zhang and Shu limiter to higher order moments of
// the conserved variables.  This limiter is locally applied, and conserves
// total mass by not adjusting cell averages.
//
// See: X. Zhang and C.-W. Shu, "Maximum-principle-satisfying and
// positivity-preserving high order schemes for conservation laws: Survey and
// new developments", Proceedings of the Royal Society A, v467 (2011),
// pp.2752-2776. for more details on the limiter
//
// See also: ApplyPosMPPLimiter (Flux Limiter).
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

    double eps=1.0e-10;

    // ------------------------------------------------ //
    // number of points where we want to check solution //
    // ------------------------------------------------ //
    const int mpoints_on_corners = 4;
    const int mpoints_on_edges   = 4*space_order;   // Included for Riemann solves
    const int mpoints_L2Proj     = (space_order)*(space_order);

    // Change this quantity if needed to check more points //
    const int mpoints   =  mpoints_on_corners+mpoints_on_edges+mpoints_L2Proj;

    // ---------------------------------------------------------- //
    // sample basis at all points where we want to check solution //
    // ---------------------------------------------------------- //
    dTensor2 spts(mpoints,2);
    void SetPositivePoints(const int& space_order, dTensor2& spts);
    SetPositivePoints(space_order,spts);

    // ---------------------------------------------------------- //
    // Sample phi at the points we picked out
    // ---------------------------------------------------------- //
    void SetLegendreAtPoints(const int& space_order, 
        const dTensor2& spts, dTensor2& phi);
    dTensor2 phi(mpoints,kmax);
    SetLegendreAtPoints(space_order, spts, phi);

    const int istep = (mx / 5)+1;
    const int jstep = ((my / 5))+1;


    // -------------------------------------------------------------- //
    // q_limited = Q1 + \theta ( q(xi,eta) - Q1 )                     //
    // where theta = min(1, |Q1| / |Q1-m|; m = min_{i} q(xi_i, eta_i) //
    // -------------------------------------------------------------- //
#pragma omp parallel for
//    for(int i=1; i <= mx; i++)
//    for(int j=1; j <= my; j++)
    for(int i=1-mbc; i <= mx+mbc; i++)
    for(int j=1-mbc; j <= my+mbc; j++)
    {

       //if(i<istep || j>=jstep)
         {
         // ---------------------------------- //
         // Limit the Density and Energy
         // ---------------------------------- //
         double mine     = 0.0;
         double minrho   = 0.0;
         double thetae   = 1.0;
         double thetarho = 1.0;
         for(int mp=1; mp <= mpoints; mp++)
         {
            // evaluate e and rho at spts(mp)
            double enow   = 0.0;
            double rhonow = 0.0;
            for( int k=1; k <= kmax; k++ )
            {
               rhonow += q.get(i,j,1,k) * phi.get(mp,k);
               enow   += q.get(i,j,5,k) * phi.get(mp,k);
            }
            minrho = Min(minrho,rhonow);
            mine   = Min(mine,enow);
         }

         double Q1 = q.get(i,j,1,1);
        if (Q1 < 0.0 || isnan(Q1))
        {    cout << " Negative dens "<<Q1<<" "<<i<<" "<<j << endl;exit(1);}
         if( fabs( Q1 - minrho ) < eps || Q1 < eps ){ thetarho = 0.0; }
         else if(  minrho <eps) { thetarho = Min( thetarho, fabs( (eps-Q1) / (minrho-Q1) ) ); }

         Q1 = q.get(i,j,5,1);
         if( fabs( Q1 - mine ) < 1.0e-14 || isnan(Q1) ){ thetae = 0.0; }
         else{ thetae = Min( thetae, fabs( (eps-Q1) / (mine-Q1) ) ); }


         // limit e and rho //
         for( int k=2; k <= kmax; k++ )
         {
            q.set(i,j,1,k, q.get(i,j,1,k) * thetarho );
            //q.set(i,j,5,k, q.get(i,j,5,k) * thetae );
         }


         // ------------------------------------------------ //
         // Limit the Momentum (to obtain positive pressure)
         //
         // This step works by considering an expansion of density, momentum
         // and energy in terms of a single parameter, theta.
         //
         // ------------------------------------------------ //

         double m = eps;
         double ave=q.get(i,j,1,1)*q.get(i,j,5,1)-0.5*(q.get(i,j,2,1)*q.get(i,j,2,1)+q.get(i,j,3,1)*q.get(i,j,3,1)+q.get(i,j,4,1)*q.get(i,j,4,1));
         double thetam = 1.0;

         
        if (ave < 0.0)
        {    
            {cout << q.get(i,j,1,1)<<" "<<q.get(i,j,5,1)<<" Negative press "<<ave<<" "<<i<<" "<<j << endl;exit(1);}
        }
         double mxa  = q.get(i,j,2,1);
         double mya  = q.get(i,j,3,1);
         double mza  = q.get(i,j,4,1);
         double ea   = q.get(i,j,5,1);
         double rhoa = q.get(i,j,1,1);
         for(int mp=1; mp <= mpoints; mp++)
         {

            // evaluate q at spts(mp) //
            double mxnow  = 0.0;
            double mynow  = 0.0;
            double mznow  = 0.0;
            double enow   = 0.0;
            double rhonow = 0.0;

            for( int k=1; k <= kmax; k++ )
            {
               mxnow  += q.get(i,j,2,k) * phi.get(mp,k);
               mynow  += q.get(i,j,3,k) * phi.get(mp,k);
               mznow  += q.get(i,j,4,k) * phi.get(mp,k);
               enow   += q.get(i,j,5,k) * phi.get(mp,k);
               rhonow += q.get(i,j,1,k) * phi.get(mp,k);
            }

            double mxc  = mxnow-mxa;
            double myc  = mynow-mya;
            double mzc  = mznow-mza;
            double ec   = enow-ea;
            double rhoc = rhonow-rhoa;

            double slope = ec*rhoa+rhoc*ea-(mxc*mxa+myc*mya+mzc*mza);
            double press = rhonow*enow-0.5*(mxnow*mxnow+mynow*mynow+mznow*mznow);

            if( fabs( slope ) < eps ){ thetam = thetam; }
            else{ thetam = Min( thetam, fabs( (eps-ave) / slope ) ); }

            m = Min(m, press);
        }

        double theta = 1.0;
        Q1 = ave;
        if( fabs( Q1 - m ) < eps ){ theta = 0.0; }
        else if( m < eps){ theta = Min( thetam, fabs( (eps-Q1) / (m-Q1) ) ); }

        // limit q //
        for (int me=1;me<=5;me++)
        for( int k=2; k <= kmax; k++ )
        {
            q.set(i,j,me,k, q.get(i,j,me,k) * theta );
        }


         for(int mp=1; mp <= mpoints; mp++)
         {

            // evaluate q at spts(mp) //
            double mxnow  = 0.0;
            double mynow  = 0.0;
            double mznow  = 0.0;
            double enow   = 0.0;
            double rhonow = 0.0;

            for( int k=1; k <= kmax; k++ )
            {
               mxnow  += q.get(i,j,2,k) * phi.get(mp,k);
               mynow  += q.get(i,j,3,k) * phi.get(mp,k);
               mznow  += q.get(i,j,4,k) * phi.get(mp,k);
               enow   += q.get(i,j,5,k) * phi.get(mp,k);
               rhonow += q.get(i,j,1,k) * phi.get(mp,k);
            }

  
            double press = enow-0.5/rhonow*(mxnow*mxnow+mynow*mynow+mznow*mznow);
         //   if(i==11 && j==3){cout<<"point "<<mp<<" rho "<<rhonow<<" u "<<mxnow/rhonow<<" v "<<mynow/rhonow<<" e "<<enow<<endl;}

            if( rhonow<0.0 || press<0.0 || isnan(rhonow)|| isnan(press) ){ cout<<i<<" "<<j<<" "<<rhonow<<" pressure probs "<<press<<endl;exit(1); }
        }
        //if(i==11 && j==3){cout<<"THIS ONE!"<<endl;
        //cout<<q.get(i,j,3,1)<<" "<<q.get(i,j,3,2)<<" "<<q.get(i,j,3,3)<<" "<<q.get(i,j,3,4)<<" "<<q.get(i,j,3,5)<<" "<<q.get(i,j,3,6)<<endl;
       // }

     }

    }

}


// This routine sets quadrature points at the following locations:
//
//   corner points   : a total of 4 of them.
//   edge points     : a total of 4*space_order quadrature of them (space_order
//                     on each edge)
//   interior points : a total of space_order^2 of them.
//
// This routine differs from what is saved in the main library.
void SetPositivePoints(const int& space_order, dTensor2& spts)
{

    dTensor1 s1d(space_order);

    // 1D Gaussian quadrature points
    switch ( space_order )
    {
        case 1:
            s1d.set(1, 0.0e0 );
            break;

        case 2:
            s1d.set(1, -1.0/sq3 );
            s1d.set(2,  1.0/sq3 );      
            break;

        case 3:
            s1d.set(1, -sq3/sq5 );
            s1d.set(2,  0.0e0 );
            s1d.set(3,  sq3/sq5 );
            break;

        case 4:
            s1d.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
            s1d.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
            s1d.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );      
            break;

        case 5:
            s1d.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
            s1d.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            s1d.set(3,  0.0 );
            s1d.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            s1d.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );      
            break;
    }

    // This region has been commented out because we are no longer applying
    // the limiter at the corner points

    // 2D points -- corners
    spts.set(1,1, -1.0e0 );
    spts.set(1,2, -1.0e0 );

    spts.set(2,1,  1.0e0 );
    spts.set(2,2, -1.0e0 );

    spts.set(3,1, -1.0e0 );
    spts.set(3,2,  1.0e0 );

    spts.set(4,1,  1.0e0 );
    spts.set(4,2,  1.0e0 );

    // 2D points -- left, right, bottom and top edges
    for (int m=1; m<=space_order; m++)
    {
        double s = s1d.get(m);

        // left edge
        spts.set(4+m,1, -1.0e0 );
        spts.set(4+m,2,  s     );

        // right edge
        spts.set(4+space_order+m,1,  1.0e0 );
        spts.set(4+space_order+m,2,  s     );

        // bottom edge
        spts.set(4+2*space_order+m,1,  s     );
        spts.set(4+2*space_order+m,2, -1.0e0 );

        // top edge
        spts.set(4+3*space_order+m,1,  s     );
        spts.set(4+3*space_order+m,2,  1.0e0 );
    }

    // 2D points -- all interior points
    int z = 4+4*space_order;
    //z = 0;
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
