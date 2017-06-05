#include <cmath>
#include <iostream>
#include "dogdefs.h"
#include "stdlib.h"
using namespace std;

// Modified version of the all purpose routine L2Project specifically written
// for projecting FluxFuncLxW onto the derivatives of the legendre basis
// functions.
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
//      F = F1 - dt/2 * F2 + dt^2/6 * F3
//
//        F1 =  f(q)
//        F2 =  f'(q) * f(q)_{x}
//        F3 =  (A_q f(q)_x f(q)_x) + A(A f(q)_x)_x
//
// stab_cor == 0 - stability corrections turned off.
// stab_cor != 0 - stability corrections turned on
//
// Parameters are given by:
//
// Integral terms:
//      IntF1 = \int( F1 \phi_x )/dx
//      IntF2 = \int( F2 \phi_x )/dx
//      IntF3 = \int( F3 \phi_x )/dx
//
//  Coefficients for boundary terms:
//      Flux     = \int( F1 \phi) / dx
//      LxWFlux2 = \int( F2 \phi) / dx
//      LxWFlux3 = \int( F3 \phi) / dx
//
void L2ProjectLxW(const int stab_cor, const int istart, const int iend, 
                  const dTensor2& node, const dTensorBC3& qin, const dTensorBC3& auxin,  
                  dTensorBC3& IntF1, dTensorBC3& IntF2, dTensorBC3& IntF3,
                  dTensorBC3& Flux, dTensorBC3& LxWFlux2, dTensorBC3& LxWFlux3,
                  void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&),
                  void (*DFunc)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor3&),
                  void (*D2Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor4&))
{    
    int l,i,k,m,me,ma,m1,m2,m3;
    int maux = auxin.getsize(2);
    int meqn = qin.getsize(2);   // number of equations
    int mpoints = qin.getsize(3);   // ( == kmax as well)
    double dx,xc,s_area,tmp,tmp1,tmp2;
    dTensor1 wgt(mpoints), spts(mpoints), xpts(mpoints);
    dTensor2 phi(mpoints,5), phi_x(mpoints,5), phi_xx(mpoints,5);
    dTensor2 qvals(mpoints,meqn), auxvals(mpoints,maux);
    dTensor2   fvals(mpoints,meqn);  // flux function, f(q)
    dTensor3  Dfvals(mpoints,meqn,meqn);  // deriv of flux function, Df(q)
    dTensor4 D2fvals(mpoints,meqn,meqn,meqn);  // D2f(q)

    dTensor2   tmpF(mpoints,meqn);// placeholder for f'(q)f_x, f'(q)f_xx,...
    dTensor2   tmpF1(mpoints,meqn);// placeholder for f'(q)f_x, f'(q)f_xx,...
    dTensor2   tmpF2(mpoints,meqn);// placeholder for f'(q)f_x, f'(q)f_xx,...

    dTensor2   f_x(mpoints,meqn); // placeholder for f_x(x_k)
    dTensor2   Af_x(mpoints,meqn); // placeholder for A*f_x(x_k)

    // -----------------
    // Quick error check
    // -----------------
    if (meqn<1 || maux <1 || mpoints<1 || mpoints>3 )
    {
        cout << " Error in L2projectLxW.cpp ... " << endl;
        cout << "         meqn = " << meqn << endl;
        cout << "         maux = " << maux << endl;
        cout << "      mpoints = " << mpoints << endl;
        cout << "      meqn = " << meqn << endl;
        cout << "       istart = " << istart << endl;
        cout << "         iend = " << iend << endl;
        cout << endl;
        exit(1);
    }

    if( meqn > 1 && mpoints > 3)
    {
        cout << " LxW method not implemented for: " << endl;
        cout << "       meqn = " << meqn << endl;
        cout << "      order = " << mpoints << endl;
        exit(1);
    }

    // ---------------------------------------------
    // Check for trivial case in the case 
    // ---------------------------------------------
    if ( mpoints == 1 )
    {
        for (i=istart; i<=iend; i++)
        for (m=1; m<=meqn; m++)
        {
            IntF1.set(i,m,1,0.0);
            IntF2.set(i,m,1,0.0);
            IntF3.set(i,m,1,0.0);
      }
  }
//  else
    {
        // ---------------------------------
        // Set quadrature weights and points
        // ---------------------------------
        //switch ( (mpoints-mopt) )
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
        dx = node.get(2,1) - node.get(1,1);
        //for (m=1; m<=(mpoints-mopt); m++)
        for (m=1; m<=(mpoints); m++)
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

            // 2st derivative of Legendre basis functions at each grid point
            phi_xx.set( m,1, 0.0 );
            phi_xx.set( m,2, 0.0 );
            phi_xx.set( m,3, 12.0*sq5*pow(dx,-2) );
            //        phi_xx.set( m,4, 60.0*sq7*pow(dx,-2)*spts.get(m) );
        }


        // ----------------------------------
        // Loop over all elements of interest
        // ----------------------------------    
        s_area = 2.0;
        for (i=istart; i<=iend; i++)
        {
            xc = node.get(1,1) + (double(i)-0.5)*dx;

            // Loop over each quadrature point
            //for (m=1; m<=(mpoints-mopt); m++)
            for (m=1; m<=(mpoints); m++)
            {
                // grid point x
                xpts.set( m, xc + 0.5*dx*spts.get(m) );

                // Solution values (q) at each grid point
                for (me=1; me<=meqn; me++)
                {
                    qvals.set(m,me, 0.0 );

                    for (k=1; k<=mpoints; k++)
                    {
                        qvals.set(m,me, qvals.get(m,me) 
                                + phi.get(m,k) * qin.get(i,me,k) );
                    }
                }

                // Auxiliary values (aux) at each grid point
                for (ma=1; ma<=maux; ma++)
                {
                    auxvals.set(m,ma, 0.0 );

                    for (k=1; k<=mpoints; k++)
                    {
                        auxvals.set(m,ma, auxvals.get(m,ma) 
                                + phi.get(m,k) * auxin.get(i,ma,k) );
                    }
                }

                // Call user-supplied function to set fvals, Dfvals, D2fvals
                Func  (xpts,qvals,auxvals,fvals);
                DFunc (xpts,qvals,auxvals,Dfvals);
                D2Func(xpts,qvals,auxvals,D2fvals);
            }

            // Evaluate integrals
            // \int f(q) \phi_x dx
            //      and project f(q) onto legendre basis:
            // f(q) = \sum_k Flux(:,:,k) \phi^{(k)}
            for (m1=1; m1<=meqn; m1++)
            {
                for (m2=1; m2<=mpoints; m2++)
                {
                    tmp1 = 0.0;
                    tmp2 = 0.0;
                    for (k=1; k<=mpoints; k++)
                    {
                        tmp1 += wgt.get(k)*fvals.get(k,m1)
                            *phi.get(k,m2);
                        tmp2 += wgt.get(k)*fvals.get(k,m1)
                            *phi_x.get(k,m2);
                    }
                    Flux.set(i,m1,m2, tmp1/s_area ); // Legendre basis
                    IntF1.set(i,m1,m2, tmp2/s_area );  // projection onto derivative
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
                for(m1=1; m1<=meqn; m1++)
                {
                    for(k=1; k<=mpoints; k++)
                    {

                        tmp = 0.0;
                        for(m2=1; m2<=meqn; m2++)
                            for(l=1;l<=mpoints;l++)
                            {
                                tmp += Dfvals.get(k,m1,m2) * Flux.get(i,m2,l) *
                                    phi_x.get(k,l); 
                            }

                        Af_x.set(k,m1,tmp);

                    }
                }

                // project f'(q)*f(q)_x onto legendre polynomials
                // by computing the integration
                for (m1=1; m1<=meqn; m1++)
                {
                    for (m2=1; m2<=mpoints; m2++)
                    {
                        tmp1 = 0.0;
                        tmp2 = 0.0;
                        for (k=1; k<=mpoints; k++)
                        {
                            tmp1 += wgt.get(k)*Af_x.get(k,m1)
                                *phi.get(k,m2);
                            tmp2 += wgt.get(k)*Af_x.get(k,m1)
                                *phi_x.get(k,m2);
                        }

                        LxWFlux2.set(i,m1,m2, - tmp1 / s_area );
                        IntF2.set   (i,m1,m2, - tmp2 / s_area );
                    }
                }
                // Finished projecting f'(q) f_x onto legendre basis
            }

            if( mpoints > 2) // third order terms
            {

                /////////////////////////////////////////////////////////////////
                // Compute term A(A*f_x)_x  where A = f'(q)
                /////////////////////////////////////////////////////////////////

                // compute pointwise values of A*(A*f_x)_x (x_k)
                // save result into tmpF(1:mpoints, 1:meqn)
                for(k=1;  k <=mpoints; k++)
                    for(m1=1; m1<=meqn; m1++)
                    {
                        // compute pointwise values of (A*f_x)_x (x_k)
                        tmp = 0.0;
                        for(m2=1; m2<=meqn; m2++)
                        {    
                            // [(A f_x)_x]_{m2} - the m2 equation evaluated at x_k
                            tmp1 = 0.0;
                            for(l=1;l<=mpoints;l++)
                            {
                                tmp1 += -LxWFlux2.get(i,m2,l) * phi_x.get(k,l);
                            }
                            tmp += Dfvals.get(k, m1, m2) * tmp1;
                        }

                        tmpF.set(k,m1,tmp);

                    }

                // project (A (A f_x)_x ) onto legendre basis
                for (m1=1; m1<=meqn; m1++)
                {
                    for (m2=1; m2<=mpoints; m2++)
                    {
                        tmp1 = 0.0;
                        tmp2 = 0.0;
                        for (k=1; k<=mpoints; k++)
                        {
                            tmp1 += wgt.get(k)*tmpF.get(k,m1)
                                *phi.get(k,m2);
                            tmp2 += wgt.get(k)*tmpF.get(k,m1)
                                *phi_x.get(k,m2);
                        }
                        LxWFlux3.set(i,m1,m2, tmp1 / s_area );
                        IntF3.set   (i,m1,m2, tmp2 / s_area );
                    }
                }

                /////////////////////////////////////////////////////////////////
                // Compute ( A_q f_x f_x )
                /////////////////////////////////////////////////////////////////

                // compute f_x(x_k)
                for(m1=1; m1<=meqn; m1++)
                    for(k=1;  k <=mpoints; k++)
                    {
                        tmp = 0.0;
                        for(l=1; l<=mpoints; l++)
                        {
                            tmp += Flux.get(i,m1,l) * phi_x.get(k,l);
                        }
                        f_x.set(k,m1, tmp);
                    }


                // compute A_q f_x f_x
                for(m1=1; m1<=meqn; m1++)
                    for(k=1; k<=mpoints; k++)
                    {
                        tmp = 0.0;
                        for(m2=1;m2<=meqn;m2++)
                            for(m3=1; m3<=meqn; m3++)
                            {
                                tmp += D2fvals.get(k, m1, m2, m3 ) *
                                    f_x.get(k,m2) * f_x.get(k,m3);
                            }
                        tmpF1.set(k,m1,tmp);
                    }

                // add in contribution from this term into LxWFlux3 ...
                for (m1=1; m1<=meqn; m1++)
                {
                    for (m2=1; m2<=mpoints; m2++)
                    {
                        tmp1 = 0.0;
                        tmp2 = 0.0;
                        for (k=1; k<=mpoints; k++)
                        {
                            tmp1 += wgt.get(k)*tmpF1.get(k,m1)
                                *phi.get(k,m2);
                            tmp2 += wgt.get(k)*tmpF1.get(k,m1)
                                *phi_x.get(k,m2);
                        }

                        LxWFlux3.set(i,m1,m2, LxWFlux3.get(i,m1,m2) + tmp1 / s_area );
                        IntF3.set   (i,m1,m2, IntF3.get(i,m1,m2) + tmp2 / s_area );
                    }
                }

            }// end mpoints > 2 case

        }

    }

}
