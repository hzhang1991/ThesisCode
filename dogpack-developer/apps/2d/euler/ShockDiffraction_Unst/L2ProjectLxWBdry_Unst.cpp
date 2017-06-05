#include "dogdefs.h"
#include "dog_math.h"
#include "mesh.h"
#include "DogParamsCart2.h"
#include "L2ProjectInline2d.h"
#include "assert.h"
#include "L2ProjectLxW_Unst.h"

// Modified version of the all purpose routine L2Project specifically written
// for projecting the "time-averaged" flux function onto the basis function.
//
// This routine also returns the coefficients of the Lax Wendroff Flux
// Function when expanded with legendre basis functions, and therefore the
// basis expansions produced by this routine can be used for all of the
// Riemann solves.
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           TODO - document the inputs here
// ---------------------------------------------------------------------
void L2ProjectLxWBdry_Unst( const int mterms,
        const double alpha, const double beta_dt, const double charlie_dt,
        const int istart, const int iend,               // Start-stop indices
        const int QuadOrder,
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,
        const mesh& Mesh, 
        const dTensor3* qin, const dTensor3* auxin,     // state vector
        dTensor3* F, dTensor3* G,                       // time-averaged Flux function
        dTensor3*  qtbdry,dTensor3* qttbdry, dTensor3* qtttbdry,iTensor1* indexbdry,
        void FluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& Aux, dTensor3& flux),
        void DFluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& aux, dTensor4& Dflux),
        void D2FluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& aux, dTensor5& D2flux) )
{    

    if( fabs( alpha ) < 1e-14 && fabs( beta_dt ) < 1e-14 && fabs( charlie_dt ) < 1e-14 )
    {
        F->setall(0.);
        G->setall(0.);
        return;
    }

    // starting and ending indices 
    const int   NumElems = Mesh.get_NumElems();
    const int NumPhysElems  = Mesh.get_NumPhysElems();
    assert_ge(istart,1);
    assert_le(iend,NumElems);

    // qin variable
    assert_eq(NumElems,qin->getsize(1));
    const int     meqn = qin->getsize(2);
    const int kmax_qin = qin->getsize(3);
    assert_eq(kmax_qin,(BasisOrder_qin*(BasisOrder_qin+1))/2);

    // auxin variable
    assert_eq(NumElems,auxin->getsize(1));
    const int       maux = auxin->getsize(2);
    const int kmax_auxin = auxin->getsize(3);
    assert_eq(kmax_auxin,(BasisOrder_auxin*(BasisOrder_auxin+1))/2);

    // fout variables
    assert_eq(NumElems,    F->getsize(1));
    const int mcomps_out = F->getsize(2);
    const int  kmax_fout = F->getsize(3);
    assert_eq(kmax_fout, (BasisOrder_fout*(BasisOrder_fout+1))/2 );

    // number of quadrature points
    assert_ge(QuadOrder, 1);
    assert_le(QuadOrder, 5);

    // Number of quadrature points
    int mpoints;
    switch( QuadOrder )
    {
        case 1:
            mpoints = 1;
            break;

        case 2:
            mpoints = 3;
            break;

        case 3:
            mpoints = 6;
            break;

        case 4:
            mpoints = 12;
            break;

        case 5:	     
            mpoints = 16;
            break;
    }

    const int kmax = iMax(iMax(kmax_qin, kmax_auxin), kmax_fout);
    dTensor2  phi(mpoints, kmax); // Legendre basis (orthogonal)
    dTensor2 spts(mpoints, 2);    // List of quadrature points
    dTensor1 wgts(mpoints);       // List of quadrature weights

    setQuadPoints_Unst( QuadOrder, wgts, spts );

    // ---------------------------------------------------------------------- //
    // Evaluate the basis functions at each point
    SetLegendreAtPoints_Unst(spts, phi);
    // ---------------------------------------------------------------------- //

    // ---------------------------------------------------------------------- //
    // First-order derivatives
    dTensor2 phi_xi (mpoints, kmax );
    dTensor2 phi_eta(mpoints, kmax );
    SetLegendreGrad_Unst( spts, phi_xi, phi_eta );
    // ---------------------------------------------------------------------- //

    // ---------------------------------------------------------------------- //
    // Second-order derivatives
    dTensor2 phi_xi2  (mpoints, kmax );
    dTensor2 phi_xieta(mpoints, kmax );
    dTensor2 phi_eta2 (mpoints, kmax );
    LegendreDiff2_Unst(spts, &phi_xi2, &phi_xieta, &phi_eta2 );
    // ---------------------------------------------------------------------- //

    // ------------------------------------------------------------- //
    // Loop over every grid cell indexed by user supplied parameters //
    // described by istart...iend, jstart...jend                     // 
    // ------------------------------------------------------------- //

    int indicator=0;

#pragma omp parallel for
    for (int i=istart; i<=iend; i++)
    {

        // These need to be defined locally.  Each mesh element carries its
        // own change of basis matrix, so these need to be recomputed for
        // each element.  The canonical derivatives, phi_xi, and phi_eta can
        // be computed and shared for each element.

        // First-order derivatives
        dTensor2   phi_x(mpoints, kmax_fout);   //   x-derivative of Legendre basis (orthogonal)
        dTensor2   phi_y(mpoints, kmax_fout);   //   y-derivative of Legendre basis (orthogonal)

        // Second-order derivatives
        dTensor2   phi_xx(mpoints, kmax_fout);   //   xx-derivative of Legendre basis (orthogonal)
        dTensor2   phi_xy(mpoints, kmax_fout);   //   xy-derivative of Legendre basis (orthogonal)
        dTensor2   phi_yy(mpoints, kmax_fout);   //   yy-derivative of Legendre basis (orthogonal)

        //find center of current cell
        const int    i1 = Mesh.get_tnode(i,1);
        const int    i2 = Mesh.get_tnode(i,2);
        const int    i3 = Mesh.get_tnode(i,3);

        // Corners:
        const double x1 = Mesh.get_node(i1,1);
        const double y1 = Mesh.get_node(i1,2);
        const double x2 = Mesh.get_node(i2,1);
        const double y2 = Mesh.get_node(i2,2);
        const double x3 = Mesh.get_node(i3,1);
        const double y3 = Mesh.get_node(i3,2);

        // Center of current cell:
        const double xc = (x1+x2+x3)/3.0;
        const double yc = (y1+y2+y3)/3.0;

        // Variables that need to be written to, and therefore are 
        // created for each thread
        dTensor2 xpts   (mpoints, 2);
        dTensor2 qvals  (mpoints, meqn);
        dTensor2 auxvals(mpoints, maux);

        // local storage for Flux function its Jacobian, and the Hessian:
        dTensor3    fvals(mpoints,             meqn, 2);  // flux function (vector)
        dTensor4        A(mpoints,       meqn, meqn, 2);  // Jacobian of flux
        dTensor5        H(mpoints, meqn, meqn, meqn, 2);  // Hessian of flux

        // Compute q, aux and fvals at each Gaussian Quadrature point
        // for this current cell indexed by (i,j)
        // Save results into dTensor2 qvals, auxvals and fvals.
        for (int m=1; m<= mpoints; m++)
        {

            // convert phi_xi and phi_eta derivatives
            // to phi_x and phi_y derivatives through Jacobian
            //
            // Note that: 
            //
            //     pd_x = J11 pd_xi + J12 pd_eta and
            //     pd_y = J21 pd_xi + J22 pd_eta.
            //
            // Squaring these operators yields the second derivatives.
            for (int k=1; k<=kmax_fout; k++)
            {
                phi_x.set(m,k, Mesh.get_jmat(i,1,1)*phi_xi.get(m,k)
                             + Mesh.get_jmat(i,1,2)*phi_eta.get(m,k) );
                phi_y.set(m,k, Mesh.get_jmat(i,2,1)*phi_xi.get(m,k)
                             + Mesh.get_jmat(i,2,2)*phi_eta.get(m,k) );

                phi_xx.set(m,k, Mesh.get_jmat(i,1,1)*Mesh.get_jmat(i,1,1)*phi_xi2.get(m,k)
                              + Mesh.get_jmat(i,1,1)*Mesh.get_jmat(i,1,2)*phi_xieta.get(m,k)
                              + Mesh.get_jmat(i,1,2)*Mesh.get_jmat(i,1,2)*phi_eta2.get(m,k)
                           );

                phi_xy.set(m,k, Mesh.get_jmat(i,1,1)*Mesh.get_jmat(i,2,1)*phi_xi2.get(m,k)
                             +(Mesh.get_jmat(i,1,2)*Mesh.get_jmat(i,2,1)
                             + Mesh.get_jmat(i,1,1)*Mesh.get_jmat(i,2,2))*phi_xieta.get(m,k)
                             + Mesh.get_jmat(i,1,2)*Mesh.get_jmat(i,2,2)*phi_eta2.get(m,k)
                           );

                phi_yy.set(m,k, Mesh.get_jmat(i,2,1)*Mesh.get_jmat(i,2,1)*phi_xi2.get(m,k)
                              + Mesh.get_jmat(i,2,1)*Mesh.get_jmat(i,2,2)*phi_xieta.get(m,k)
                              + Mesh.get_jmat(i,2,2)*Mesh.get_jmat(i,2,2)*phi_eta2.get(m,k)
                           );
            }

            // point on the unit triangle
            const double s = spts.get(m,1);
            const double t = spts.get(m,2);

            // point on the physical triangle
            xpts.set(m,1, xc + (x2-x1)*s + (x3-x1)*t );
            xpts.set(m,2, yc + (y2-y1)*s + (y3-y1)*t );

            // Solution values (q) at each grid point
            for (int me=1; me<=meqn; me++)
            {
                qvals.set(m,me, 0.0 );
                for (int k=1; k<=kmax_qin; k++)
                {
                    qvals.set(m,me, qvals.get(m,me) 
                            + phi.get(m,k) * qin->get(i,me,k) );
                }
            }

            // Auxiliary values (aux) at each grid point
            for (int ma=1; ma<=maux; ma++)
            {
                auxvals.set(m,ma, 0.0 );
                for (int k=1; k<=kmax_auxin; k++)
                {
                    auxvals.set(m,ma, auxvals.get(m,ma) 
                            + phi.get(m,k) * auxin->get(i,ma,k) );
                }
            } 
        }

        // ----------------------------------------------------------------- //
        //
        // Part I:
        //
        // Project the flux function onto the basis 
        // functions.  This is the term of order O( 1 ) in the
        // "time-averaged" Taylor expansion of f and g.
        //
        // ----------------------------------------------------------------- //

        // Call user-supplied function to set fvals
        FluxFunc(xpts, qvals, auxvals, fvals);

        // Evaluate integral on current cell (project onto Legendre basis) 
        // using Gaussian Quadrature for the integration
        //
        // TODO - do we want to optimize this by looking into using transposes,
        // as has been done in the 2d/cart code? (5/14/2014) -DS
        for (int me=1; me<=mcomps_out; me++)		
        for (int k=1; k<=kmax; k++)
        {
            double tmp1 = 0.0;
            double tmp2 = 0.0;
            for (int mp=1; mp <= mpoints; mp++)
            {
                tmp1 += wgts.get(mp)*fvals.get(mp, me, 1)*phi.get(mp, k);
                tmp2 += wgts.get(mp)*fvals.get(mp, me, 2)*phi.get(mp, k);
            }
            F->set(i, me, k,  2.0*tmp1 );
            G->set(i, me, k,  2.0*tmp2 );
        }

        // ----------------------------------------------------------------- //
        //
        // Part II:
        //
        // Project the derivative of the flux function onto the basis 
        // functions.  This is the term of order O( \dt ) in the
        // "time-averaged" Taylor expansion of f and g.
        //
        // ----------------------------------------------------------------- //

        // ----------------------------------------------------------------- //
        // Compute pointwise values for fx+gy:
        //
        // We can't multiply fvals of f, and g,
        // by alpha, otherwise we compute the wrong derivative here!
        //
        dTensor2 fx_plus_gy( mpoints, meqn ); fx_plus_gy.setall(0.);
        for( int mp=1; mp <= mpoints; mp++ )
        for( int me=1; me <= meqn; me++ )
        {
            double tmp = 0.;
            for( int k=2; k <= kmax; k++ )                
            {
                tmp += F->get( i, me, k ) * phi_x.get( mp, k );
                tmp += G->get( i, me, k ) * phi_y.get( mp, k );
            }
            fx_plus_gy.set( mp, me, tmp );
        }

        // Call user-supplied Jacobian to set f'(q) and g'(q):
        DFluxFunc( xpts, qvals, auxvals, A );

        // place-holders for point values of
        // f'(q)( fx + gy ) and g'(q)( fx + gy ):
        dTensor2 dt_times_fdot( mpoints, meqn );
        dTensor2 dt_times_gdot( mpoints, meqn );

        // Compute point values for f'(q) * (fx+gy) and g'(q) * (fx+gy):
        for( int mp=1; mp <= mpoints; mp++ )
        for( int m1=1; m1 <= meqn; m1++ )
        {
            double tmp1 = 0.;
            double tmp2 = 0.;
            for( int m2=1; m2 <= meqn; m2++ )
            {
                tmp1 += A.get(mp, m1, m2, 1 ) * fx_plus_gy.get(mp, m2);
                tmp2 += A.get(mp, m1, m2, 2 ) * fx_plus_gy.get(mp, m2);
            }
            dt_times_fdot.set( mp, m1, -beta_dt*tmp1 );
            dt_times_gdot.set( mp, m1, -beta_dt*tmp2 );
        }

        // ---  Third-order terms --- //
        //
        // These are the terms that are O( \dt^2 ) in the "time-averaged"
        // flux function.
        dTensor2 f_tt( mpoints, meqn );   f_tt.setall(0.);
        dTensor2 g_tt( mpoints, meqn );   g_tt.setall(0.);
        dTensor2 fx_plus_gy_t( mpoints, meqn );
        if( mterms > 2 )
        {

            // Construct the Hessian at each (quadrature) point
            D2FluxFunc( xpts, qvals, auxvals, H );

            // Second-order derivative terms
            dTensor2 qx_vals (mpoints, meqn);   qx_vals.setall(0.);
            dTensor2 qy_vals (mpoints, meqn);   qy_vals.setall(0.);

            dTensor2 fxx_vals(mpoints, meqn);   fxx_vals.setall(0.);
            dTensor2 gxx_vals(mpoints, meqn);   gxx_vals.setall(0.);

            dTensor2 fxy_vals(mpoints, meqn);   fxy_vals.setall(0.);
            dTensor2 gxy_vals(mpoints, meqn);   gxy_vals.setall(0.);

            dTensor2 fyy_vals(mpoints, meqn);   fyy_vals.setall(0.);
            dTensor2 gyy_vals(mpoints, meqn);   gyy_vals.setall(0.);

            for( int m=1; m <= mpoints; m++ )
            for( int me=1; me <= meqn; me++ )
            {
                // Can start at k=1, because derivative of a constant is
                // zero.
                double tmp_qx = 0.;
                double tmp_qy = 0.;
                for( int  k=2; k <= kmax; k++   )
                {
                    tmp_qx += phi_x.get(m,k) * qin->get(i,me,k);
                    tmp_qy += phi_y.get(m,k) * qin->get(i,me,k);
                }
                qx_vals.set(m,me, tmp_qx );
                qy_vals.set(m,me, tmp_qy );

                // First non-zero terms start at third-order.
                for( int  k=4; k <= kmax; k++   )
                {
                    fxx_vals.set(m,me, fxx_vals.get(m,me) + phi_xx.get(m,k)*F->get(i,me,k) );
                    gxx_vals.set(m,me, gxx_vals.get(m,me) + phi_xx.get(m,k)*G->get(i,me,k) );

                    fxy_vals.set(m,me, fxy_vals.get(m,me) + phi_xy.get(m,k)*F->get(i,me,k) );
                    gxy_vals.set(m,me, gxy_vals.get(m,me) + phi_xy.get(m,k)*G->get(i,me,k) );

                    fyy_vals.set(m,me, fyy_vals.get(m,me) + phi_yy.get(m,k)*F->get(i,me,k) );
                    gyy_vals.set(m,me, gyy_vals.get(m,me) + phi_yy.get(m,k)*G->get(i,me,k) );
                }

            }

            // ----------------------------------- //
            // Part I: Compute (f_x + g_y)_{,t}
            // ----------------------------------- //

            // Compute terms that get multiplied by \pd2{ f }{ q } and \pd2{ g }{ q }.
            for( int  m = 1;  m <= mpoints; m++ )
            for( int me = 1; me <= meqn; me++   )
            {

                double tmp = 0.;

                // Terms that get multiplied by the Hessian:
                for( int m1=1; m1 <= meqn; m1++ )
                for( int m2=1; m2 <= meqn; m2++ )
                {

                    tmp += H.get(m,me,m1,m2,1)*qx_vals.get(m,m1)*fx_plus_gy.get(m,m2);
                    tmp += H.get(m,me,m1,m2,2)*qy_vals.get(m,m1)*fx_plus_gy.get(m,m2);
                }

                // Terms that get multiplied by f'(q) and g'(q):
                for( int m1=1; m1 <= meqn; m1++ )
                {

                    tmp += A.get(m,me,m1,1)*( fxx_vals.get(m,m1)+gxy_vals.get(m,m1) );
                    tmp += A.get(m,me,m1,2)*( fxy_vals.get(m,m1)+gyy_vals.get(m,m1) );
                }

                fx_plus_gy_t.set( m, me, tmp );
            }

            // ----------------------------------- //
            // Part II: Compute 
            //      f'(q) * fx_plus_gy_t and 
            //      g'(q) * fx_plus_gy_t
            // ----------------------------------- //

            // Add in the third term that gets multiplied by A:
            for( int m=1; m <= mpoints; m++ )
            for( int m1=1; m1 <= meqn; m1++ )
            {
                double tmp1 = 0.;
                double tmp2 = 0.;
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp1 += A.get(m,m1,m2,1)*fx_plus_gy_t.get(m,m2);
                    tmp2 += A.get(m,m1,m2,2)*fx_plus_gy_t.get(m,m2);
                }
                f_tt.set( m, m1, tmp1 );
                g_tt.set( m, m1, tmp2 );
            }

            // ----------------------------------------------- //
            // Part III: Add in contributions from
            //      f''(q) * (fx_plus_gy, fx_plus_gy ) and 
            //      g''(q) * (fx_plus_gy, fx_plus_gy ).
            // ----------------------------------------------- //
            for( int m =1; m <= mpoints; m++ )
            for( int me =1; me <= meqn; me++ )
            {
                double tmp1 = 0.;
                double tmp2 = 0.;

                // Terms that get multiplied by the Hessian:
                for( int m1=1; m1 <= meqn; m1++ )
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp1 += H.get(m,me,m1,m2,1)*fx_plus_gy.get(m,m1)*fx_plus_gy.get(m,m2);
                    tmp2 += H.get(m,me,m1,m2,2)*fx_plus_gy.get(m,m1)*fx_plus_gy.get(m,m2);
                }

                f_tt.set( m, me, f_tt.get(m,me) + tmp1 );
                g_tt.set( m, me, g_tt.get(m,me) + tmp2 );
            }

        } // End of computing "third"-order terms

        // ---------------------------------------------------------- //
        // 
        // Construct basis coefficients (integrate_on_current_cell)
        //
        // ---------------------------------------------------------- //
        for (int me=1; me<=mcomps_out; me++)		
        for (int k=1; k<=kmax; k++)
        {

            double tmp1 = 0.0;
            double tmp2 = 0.0;
            for (int mp=1; mp<=mpoints; mp++)
            {
                tmp1 += wgts.get(mp)*phi.get(mp,k)*(
                    dt_times_fdot.get(mp, me) + charlie_dt*f_tt.get(mp, me) );
                tmp2 += wgts.get(mp)*phi.get(mp,k)*(
                    dt_times_gdot.get(mp, me) + charlie_dt*g_tt.get(mp, me) );
            }
            F->set(i,me,k,  F->get(i,me,k) + 2.0*tmp1 );
            G->set(i,me,k,  G->get(i,me,k) + 2.0*tmp2 );

        }

        int efound1=0;
        int bdryindic=0;
        if(i<=NumPhysElems){
        for(int e=1;e<=3;e++)
        {
          int iedge=Mesh.get_tedge(i,e);
        // Elements on either side of edge
          int ileft  = Mesh.get_eelem(iedge,1);
          int iright = Mesh.get_eelem(iedge,2);
          if (ileft>NumPhysElems || iright>NumPhysElems){bdryindic=1;}
        }
        }

           if(bdryindic)
           {
                dTensor2 ftt(meqn,kmax);
                dTensor2 gtt(meqn,kmax);
                dTensor2 fttx_plus_gtty(meqn,mpoints);fttx_plus_gtty.setall(0.0);
                dTensor2 qt(meqn,kmax);
                dTensor2 qtt(meqn,kmax);
                dTensor2 qttt(meqn,kmax);

                for( int me=1; me <= meqn; me++ )
                for( int  k=1;  k <= kmax;  k++ )
                {
                   double tmp1 = 0.;
                   double tmp2 = 0.;
                   for( int m=1; m <= mpoints; m++ )
                   {
                       tmp1 += wgts.get(m)*phi.get(m,k)*(
                        f_tt.get(m,me) );
                       tmp2 += wgts.get(m)*phi.get(m,k)*(
                        g_tt.get(m,me) );

//                  tmp1 += wgt_phi_transpose.get(k,m)*( dt_times_fdot.get(m, me ) );
//                  tmp2 += wgt_phi_transpose.get(k,m)*( dt_times_gdot.get(m, me ) );
                    }
                    ftt.set(me,k, 0.25*tmp1 );
                    gtt.set(me,k, 0.25*tmp2 );
                 }

                 for( int m=1; m <= mpoints; m++ )
                 for( int me=1; me <= meqn; me++ )
                 {
                    double tmp = 0.;
                    for( int k=2; k <= kmax; k++ )
                    {
                        tmp += ftt.get( me, k ) * phi_x.get( m, k );
                        tmp += gtt.get( me, k ) * phi_y.get( m, k );
                     }
                     fttx_plus_gtty.set( me,m, tmp );
                 }

                for( int me=1; me <= meqn; me++ )
                for( int  k=1;  k <= kmax;  k++ )
                {
                   double tmp1 = 0.;
                   double tmp2 = 0.;
                   double tmp3 = 0.;
                   for( int m=1; m <= mpoints; m++ )
                   {
                       tmp3 += wgts.get(m)*phi.get(m,k)*(
                         fttx_plus_gtty.get(me,m) );
                       tmp1 += wgts.get(m)*phi.get(m,k)*(
                         fx_plus_gy.get(m,me) );
                       tmp2 += wgts.get(m)*phi.get(m,k)*(
                         fx_plus_gy_t.get(m,me) );
                    }
                    qt.set(me,k, 0.25*tmp1 );
                    qtt.set(me,k, 0.25*tmp2 );
                    qttt.set(me,k, 0.25*tmp3 );
                 }


                //int j1=i-NumPhysElems;
                indicator=indicator+1;
                for( int me=1; me <= meqn; me++ )
                for( int  k=1;  k <= kmax;  k++ )
                {
                        qtbdry->set(indicator,me,k, qt.get(me,k) );
                        qttbdry->set(indicator,me,k, qtt.get(me,k) );
                        qtttbdry->set(indicator,me,k, qttt.get(me,k) );
                        indexbdry->set(i,indicator);
                }
    }
  }
}
