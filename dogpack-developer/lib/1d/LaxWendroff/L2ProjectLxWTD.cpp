#include <cmath>
#include <iostream>
#include "dogdefs.h"
#include "stdlib.h"
using namespace std;

// Modified version of the all purpose routine L2Project specifically written
// for projecting FluxFuncLxW onto the derivatives of the legendre basis
// functions.
//
// This routine also returns the coefficients of the Lax Wendroff Flux
// Function when expanded with legendre basis functions.
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           dTensor2    node(mnodes,1)
//           dTensorBC3 auxin(1-mbc:mnodes+mbc,maux,mpoints)
//           dTensorBC3   qin(1-mbc:mnodes+mbc,meqn,mpoints)
// ---------------------------------------------------------------------
//
//        F1 =  f(q)
//        F2 =  -f'(q) * f(q)_{x} * dt
//
// Parameters are given by:
//
// Integral terms:
//      IntF1 = \int( F1 \phi_x )/dx
//      IntF2 = \int( F2 \phi_x )/dx
//
//  Coefficients for boundary terms:
//      Flux     = \int( F1 \phi) / dx
//      LxWFlux2 = \int( F2 \phi) / dx
//
// See also L2projectLxW for a version of this function that also includes the
// third derivative.
//
void L2ProjectLxWTD(const int istart, const int iend, 
        const dTensor2& node, const dTensorBC3& qin, const dTensorBC3& auxin,  
        dTensorBC3& Flux, dTensorBC3& LxWFlux2, 
        dTensorBC3& IntF1, dTensorBC3& IntF2, 
        void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&),
        void (*DFunc)(const dTensor1&, const dTensor2&, const dTensor2&,
        dTensor3&) )
{    

    const int maux    = auxin.getsize(2);
    const int meqn    = qin.getsize(2);   // number of equations
    const int mpoints = qin.getsize(3);   // ( == kmax as well)

    dTensor1 wgt(mpoints), spts(mpoints);
    dTensor2 phi(mpoints,5), phi_x(mpoints,5);

    // -----------------
    // Quick error check
    // -----------------
    if( meqn > 1 && mpoints > 5)
    {
        cout << " Errror in L2projectLxW.cpp ... " << endl;
        cout << "     not implemented for: " << endl;
        cout << "       meqn = " << meqn << endl;
        cout << "      order = " << mpoints << endl;
        exit(1);
    }

    // ---------------------------------------------
    // Check for trivial case in the case 
    // ---------------------------------------------
    if ( mpoints == 1 )
    {
        IntF1.setall(0.);
        IntF2.setall(0.);
    }
    else
    {
        // ---------------------------------
        // Set quadrature weights and points
        // ---------------------------------
        switch ( mpoints )
        {
            case 1:
                wgt.set(1,  2.0e0 );

                spts.set(1, 0.0e0 );

                break;

            case 2:
                wgt.set(1,   1.0 );
                wgt.set(2,   1.0 );

                spts.set(1, -1.0/sq3 );
                spts.set(2,  1.0/sq3 );

                break;

            case 3:
                wgt.set(1, 5.0e0/9.0e0 );
                wgt.set(2, 8.0e0/9.0e0 );
                wgt.set(3, 5.0e0/9.0e0 );

                spts.set(1, -sq3/sq5 );
                spts.set(2,  0.0e0 );
                spts.set(3,  sq3/sq5 );

                break;

            case 4:
                wgt.set(1, (18.0 - sqrt(30.0))/36.0 );
                wgt.set(2, (18.0 + sqrt(30.0))/36.0 );
                wgt.set(3, (18.0 + sqrt(30.0))/36.0 );
                wgt.set(4, (18.0 - sqrt(30.0))/36.0 );

                spts.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
                spts.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
                spts.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
                spts.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );

                break;

            case 5:
                wgt.set(1, (322.0 - 13.0*sqrt(70.0))/900.0 );
                wgt.set(2, (322.0 + 13.0*sqrt(70.0))/900.0 );
                wgt.set(3,  128.0/225.0 );
                wgt.set(4, (322.0 + 13.0*sqrt(70.0))/900.0 );
                wgt.set(5, (322.0 - 13.0*sqrt(70.0))/900.0 );

                spts.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
                spts.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
                spts.set(3,  0.0 );
                spts.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
                spts.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );

                break;
        }

        // Loop over each quadrature point to construct Legendre polys
        const double dx = node.get(2,1) - node.get(1,1);
        for (int m=1; m<=(mpoints); m++)
        {
            // Legendre basis functions at each grid point
            phi.set( m,1, 1.0 );
            phi.set( m,2, sq3*spts.get(m) );
            phi.set( m,3, 0.5*sq5*( 3.0*pow(spts.get(m),2) - 1.0 ) );
            phi.set( m,4, 0.5*sq7*spts.get(m)
                    *(5.0*pow(spts.get(m),2) - 3.0) );
            phi.set( m,5, (105.0/8.0)*pow(spts.get(m),4) 
                    - (45.0/4.0)*pow(spts.get(m),2) + (9.0/8.0) );

            // 1st derivative of Legendre basis functions at each grid point
            phi_x.set( m,1, 0.0 );
            phi_x.set( m,2, 2.0*sq3/dx );
            phi_x.set( m,3, 6.0*sq5*spts.get(m)/dx );
            phi_x.set( m,4, 3.0*sq7*(5.0*pow(spts.get(m),2)-1.0)/dx );
            phi_x.set( m,5, 15.0*spts.get(m)*
                    (7.0*pow(spts.get(m),2)-3.0)/dx );

        }


        // ----------------------------------
        // Loop over all elements of interest
        // ----------------------------------    
        const double s_area = 2.0;
//#pragma omp paralell for
        for (int i=istart; i<=iend; i++)
        {

            // Local storage (so loop can be parallelized)
            dTensor2 qvals(mpoints,meqn), auxvals(mpoints,maux);
            dTensor2   fvals(mpoints,meqn);       // flux function, f(q)
            dTensor3  Dfvals(mpoints,meqn,meqn);  // deriv of flux function, Df(q)

            // placeholders for f'(q)f_x, f'(q)f_xx, ...
            dTensor2   tmpF  (mpoints,meqn);
            dTensor2   tmpF1 (mpoints,meqn);
            dTensor2   tmpF2 (mpoints,meqn);

            // placeholders for f_x( x_k ) and A*f_x( x_k )
            dTensor2   f_x  (mpoints, meqn);
            dTensor2   Af_x (mpoints, meqn);

            dTensor1 xpts(mpoints);

            const double xc = node.get(1,1) + (double(i)-0.5)*dx;

            // Loop over each quadrature point
            for (int m=1; m<=(mpoints); m++)
            {
                // grid point x
                xpts.set( m, xc + 0.5*dx*spts.get(m) );

                // Solution values (q) at each grid point
                for (int me=1; me<=meqn; me++)
                {
                    qvals.set(m,me, 0.0 );

                    for (int k=1; k<=mpoints; k++)
                    {
                        qvals.set(m,me, qvals.get(m,me) 
                                + phi.get(m,k) * qin.get(i,me,k) );
                    }
                }

                // Auxiliary values (aux) at each grid point
                for (int ma=1; ma<=maux; ma++)
                {
                    auxvals.set(m,ma, 0.0 );

                    for (int k=1; k<=mpoints; k++)
                    {
                        auxvals.set(m,ma, auxvals.get(m,ma) 
                                + phi.get(m,k) * auxin.get(i,ma,k) );
                    }
                }

                // Call user-supplied function to set fvals, Dfvals
                Func  (xpts,qvals,auxvals,fvals);
                DFunc (xpts,qvals,auxvals,Dfvals);

            }

            // Evaluate integrals
            // \int f(q) \phi_x dx
            //      and project f(q) onto legendre basis:
            // f(q) = \sum_k Flux(:,:,k) \phi^{(k)}
            for (int m1=1; m1<=meqn; m1++)
            {
                for (int m2=1; m2<=mpoints; m2++)
                {
                    double tmp1 = 0.0;
                    double tmp2 = 0.0;
                    for (int k=1; k<=mpoints; k++)
                    {
                        tmp1 += wgt.get(k)*fvals.get(k,m1)
                            *phi.get(k,m2);
                        tmp2 += wgt.get(k)*fvals.get(k,m1)
                            *phi_x.get(k,m2);
                    }
                    Flux.set  (i,m1,m2, tmp1/s_area );  // Legendre basis
                    IntF1.set (i,m1,m2, tmp2/s_area );  // projection onto derivative
                }
            }

            //////////////////////////////////////////////////////////////////////
            // Compute the higher order terms in Taylor Expansion
            //////////////////////////////////////////////////////////////////////
            if( mpoints > 1 ) // 2nd order terms
            {

                /////////////////////////////////////////////////////////////////
                // Project f'(q) f_x onto legendre basis
                /////////////////////////////////////////////////////////////////

                // evaluate pointwise values of 
                // f'(q) f(q)_x ( x_k ) 
                for(int m1=1; m1<=meqn; m1++)
                {
                    for(int k=1; k<=mpoints; k++)
                    {

                        double tmp = 0.0;
                        for(int m2=1; m2<=meqn; m2++)
                            for(int l=1;l<=mpoints;l++)
                            {
                                tmp += Dfvals.get(k,m1,m2) * Flux.get(i,m2,l) *
                                    phi_x.get(k,l); 
                            }

                        Af_x.set(k,m1,tmp);

                    }
                }

                // project f'(q)*f(q)_x onto legendre polynomials
                // by computing the integration
                for (int m1=1; m1<=meqn; m1++)
                {
                    for (int m2=1; m2<=mpoints; m2++)
                    {
                        double tmp1 = 0.0;
                        double tmp2 = 0.0;
                        for (int k=1; k<=mpoints; k++)
                        {
                            tmp1 += wgt.get(k)*Af_x.get(k,m1)
                                *phi.get(k,m2);
                            tmp2 += wgt.get(k)*Af_x.get(k,m1)
                                *phi_x.get(k,m2);
                        }

                        LxWFlux2.set(i,m1,m2, -tmp1 / s_area );
                        IntF2.set   (i,m1,m2, -tmp2 / s_area );
                    }
                }
                // Finished projecting f'(q) f_x onto legendre basis
            }

        }

    }

}
