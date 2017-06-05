#include "dogdefs.h"
#include "dog_math.h"
#include "DogParamsCart2.h"
#include "L2ProjectInline2d.h"

// Modified version of the all purpose routine L2Project specifically written
// for projecting FluxFuncLxW onto the derivatives of the legendre basis
// functions.
//
// This routine also returns the coefficients of the Lax Wendroff Flux
// Function when expanded with legendre basis functions, and therefore the
// basis expansions produced by this routine can be used for all of the
// Riemann solves.
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           TODO
// ---------------------------------------------------------------------
void L2ProjectLxW( const int mterms,
        const double alpha, const double beta_dt, const double charlie_dt,
        const int istart, const int iend, 
        const int jstart, const int jend,
        const int QuadOrder,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,    
        const dTensorBC4* qin,
        const dTensorBC4* auxin,
        dTensorBC4* F,
        dTensorBC4* G,
        void FluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& Aux, dTensor3& flux),
        void DFluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& aux, dTensor4& Dflux),
        void D2FluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& aux, dTensor5& D2flux) )
{    

    if( fabs( alpha ) < 1e-14 && fabs( beta_dt ) < 1e-14 )
    {
        F->setall(0.);
        G->setall(0.);
        return;
    }

    // qin variable
    const int       mx = qin->getsize(1);
    const int       my = qin->getsize(2);
    const int     meqn = qin->getsize(3);
    const int     kmax = qin->getsize(4);

    // auxin variable
    assert_eq(mx,   auxin->getsize(1) );
    assert_eq(my,   auxin->getsize(2) );
    assert_eq(kmax, auxin->getsize(4) );
    const int       maux = auxin->getsize(3);

    // F variables
    assert_eq(mx,   F->getsize(1) );
    assert_eq(my,   F->getsize(2) );
    assert_eq(meqn, F->getsize(3) );
    assert_eq(kmax, F->getsize(4) );

    // mbc
    const int mbc = qin->getmbc();
    assert_eq(mbc,auxin->getmbc());
    assert_eq(mbc,F->getmbc());

    // starting and ending indeces
    assert_ge(istart, 1-mbc  );
    assert_le(iend,   mx+mbc );
    assert_ge(jstart, 1-mbc  );
    assert_le(jend,   my+mbc );

    // Grid information:
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    // number of quadrature points
    assert_ge( QuadOrder,  1 );
    assert_le( QuadOrder, 20 );
    int QuadOrder_MOD = QuadOrder;
    switch(QuadOrder)
    {
        case 7:
            QuadOrder_MOD = 8;
            break;
        case 9:
            QuadOrder_MOD = 10;
            break;
        case 11:
            QuadOrder_MOD = 12;
            break;
        case 13:
            QuadOrder_MOD = 14;
            break;
        case 15:
            QuadOrder_MOD = 16;
            break;
        case 17:
            QuadOrder_MOD = 18;
            break;
        case 19:
            QuadOrder_MOD = 20;
            break;
    }
    const int mpoints = QuadOrder_MOD*QuadOrder_MOD;

    // ---------------------------------------------------------------------- //
    // evaluate the basis functions, phi and its derivatives at all quadrature
    // points
    // ---------------------------------------------------------------------- //

    // set quadrature point and weight information
    dTensor1  wgt( mpoints    );
    dTensor2 spts( mpoints, 2 );
    void SetQuadWgtsPts(const int QuadOrder, dTensor1& wgt, dTensor2& spts );
    SetQuadWgtsPts(QuadOrder_MOD, wgt,spts);

    // Basis functions
    dTensor2 phi  (mpoints, kmax );

    // First-order derivatives
    dTensor2 phi_x(mpoints, kmax );
    dTensor2 phi_y(mpoints, kmax );

    // Second-order derivatives
    dTensor2 phi_xy(mpoints, kmax );
    dTensor2 phi_xx(mpoints, kmax );
    dTensor2 phi_yy(mpoints, kmax );

    // phi itself:
    void SetLegendrePolys(const int,const int,const dTensor2&,dTensor2&);
    SetLegendrePolys(mpoints, kmax, spts, phi);

    // phi's derivatives, phi_x and phi_y:
    void SetLegendrePolysGrad(const double dx, const double dy,
            const int mpoints, const int kmax, const dTensor2& spts, 
            dTensor2& phi_x, dTensor2& phi_y);
    SetLegendrePolysGrad(dx, dy, mpoints, kmax, spts, phi_x, phi_y);

    // phi's second derivatives, phi_xx, phi_xy, and phi_yy
    void LegendreDiff2(const double dx, const double dy, const dTensor2& spts, 
        dTensor2* phi_xx, dTensor2* phi_xy, dTensor2* phi_yy );
    LegendreDiff2(dx, dy, spts, &phi_xx, &phi_xy, &phi_yy );

    // ---------------------------------------------------------------------- //
    // For efficiency compute weight*phi and then take transpose
    dTensor2 wgt_phi_transpose    (kmax, mpoints);
//  dTensor2 wgt_phi_x_transpose  (kmax, mpoints);
//  dTensor2 wgt_phi_y_transpose  (kmax, mpoints);
//  dTensor2 wgt_phi_xx_transpose (kmax, mpoints);
//  dTensor2 wgt_phi_xy_transpose (kmax, mpoints);
//  dTensor2 wgt_phi_yy_transpose (kmax, mpoints);
    for(int mp=1; mp<=mpoints; mp++)
    for(int k=1;  k<=kmax; k++)
    {
        wgt_phi_transpose.set    (k, mp, wgt.get(mp) *    phi.get(mp,k) );
//      wgt_phi_x_transpose.set  (k, mp, wgt.get(mp) *  phi_x.get(mp,k) );
//      wgt_phi_y_transpose.set  (k, mp, wgt.get(mp) *  phi_y.get(mp,k) );
//      wgt_phi_xx_transpose.set (k, mp, wgt.get(mp) * phi_xx.get(mp,k) );
//      wgt_phi_xy_transpose.set (k, mp, wgt.get(mp) * phi_xy.get(mp,k) );
//      wgt_phi_yy_transpose.set (k, mp, wgt.get(mp) * phi_yy.get(mp,k) );
    }
    // ---------------------------------------------------------------------- //

    // ------------------------------------------------------------- //
    // Loop over every grid cell indexed by user supplied parameters //
    // described by istart...iend, jstart...jend                     // 
    // ------------------------------------------------------------- //
#pragma omp parallel for
    for (int i=istart; i<=iend; i++)
    {
        // Local storage for q, aux and xpts (all passed into user supplied
        // function)
        dTensor2   qvals(mpoints, meqn);
        dTensor2 auxvals(mpoints, maux);
        dTensor2    xpts(mpoints, 2);

        // local storage for Flux function its Jacobian, and the Hessian:
        dTensor3    fvals(mpoints,             meqn, 2);  // flux function (vector)
        dTensor4        A(mpoints,       meqn, meqn, 2);  // Jacobian of flux
        dTensor5        H(mpoints, meqn, meqn, meqn, 2);  // Hessian of flux

        for (int j=jstart; j<=jend; j++)
        {
            //find center of current cell
            const double xc = dogParamsCart2.get_xc(i);
            const double yc = dogParamsCart2.get_yc(j);

            // ------------------------------------------------------------- //
            //
            // Part I:
            //
            // Project the flux function onto the basis 
            // functions.  This is the term of order O( 1 ) in the
            // "time-averaged" Taylor expansion of f and g.
            //
            // ------------------------------------------------------------- //

            // Compute q, and aux at each Gaussian quadrature point
            // for this current cell indexed by (i,j)
            // Save results into qvals, and auxvals.
            L2ProjectInline2d::set_vals_at_each_Gaussian_quadrature_point(i, j, 
                    mpoints, meqn, maux, kmax, kmax, 
                    xc, yc, dx, dy, 
                    spts, phi, qin, auxin, xpts, qvals, auxvals);

            // ------------------------------------------------------------- //
            // Project the flux function onto the basis functions.  This is
            // the term of order O( 1 ) in the Taylor expansions of f_x + g_y.
            // ------------------------------------------------------------- //

            // Evaluate the Flux function at Gaussian quadrature point
            FluxFunc( xpts, qvals, auxvals, fvals );

            // copy onto smaller arrays (so we can call
            // integrate_on_current_cell):
            //
            // This is a bit of overkill, but we need to do something to
            // extra f, without the alpha in there.
            dTensor2 f( mpoints, meqn );
            dTensor2 g( mpoints, meqn );
            for( int m=1; m <= mpoints; m++ )
            for( int me=1; me <= meqn ; me++ )
            {
                f.set(m, me, fvals.get(m, me, 1 ) );
                g.set(m, me, fvals.get(m, me, 2 ) );
            }

            // Evaluate integral on current cell (project onto Legendre basis)
            // using Gaussian Quadrature for the integration
            //
            // Note: this does not have the factor of ALPHA in here yet!
            //
            L2ProjectInline2d::integrate_on_current_cell(false, i, j, 
                    meqn, kmax, mpoints, wgt_phi_transpose, f, F);
            L2ProjectInline2d::integrate_on_current_cell(false, i, j, 
                    meqn, kmax, mpoints, wgt_phi_transpose, g, G);

            // ------------------------------------------------------------- //
            //
            // Part II:
            //
            // Project the derivative of the flux function onto the basis 
            // functions.  This is the term of order O( \dt ) in the
            // "time-averaged" Taylor expansion of f and g.
            //
            // ------------------------------------------------------------- //


            // ------------------------------------------------------------- //
            // Compute pointwise values for fx+gy:
            //
            // We can't multiply fvals of f, and g,
            // by alpha, otherwise we compute the wrong derivative here!
            //
            dTensor2 fx_plus_gy( mpoints, meqn ); fx_plus_gy.setall(0.);
            for( int m=1; m <= mpoints; m++ )
            for( int me=1; me <= meqn; me++ )
            {
                double tmp = 0.;
                for( int k=2; k <= kmax; k++ )                
                {
                    tmp += F->get( i, j, me, k ) * phi_x.get( m, k );
                    tmp += G->get( i, j, me, k ) * phi_y.get( m, k );
                }
                fx_plus_gy.set( m, me, tmp );
            }

            // Call user-supplied Jacobian to set f'(q) and g'(q):
            DFluxFunc( xpts, qvals, auxvals, A );

            // place-holders for point values of
            // f'(q)( fx + gy ) and g'(q)( fx + gy ):
            dTensor2 dt_times_fdot( mpoints, meqn );
            dTensor2 dt_times_gdot( mpoints, meqn );

            // Compute point values for f'(q) * (fx+gy) and g'(q) * (fx+gy):
            for( int m=1; m <= mpoints; m++ )
            for( int m1=1; m1 <= meqn; m1++ )
            {
                double tmp1 = 0.;
                double tmp2 = 0.;
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp1 += A.get(m, m1, m2, 1 ) * fx_plus_gy.get(m, m2);
                    tmp2 += A.get(m, m1, m2, 2 ) * fx_plus_gy.get(m, m2);
                }
                dt_times_fdot.set( m, m1, -beta_dt*tmp1 );
                dt_times_gdot.set( m, m1, -beta_dt*tmp2 );
            }

            // ---  Third-order terms --- //
            //
            // These are the terms that ar O( \dt^2 ) in the "time-averaged"
            // flux function.
            dTensor2 f_tt( mpoints, meqn );   f_tt.setall(0.);
            dTensor2 g_tt( mpoints, meqn );   g_tt.setall(0.);
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
                        tmp_qx += phi_x.get(m,k) * qin->get(i,j,me,k);
                        tmp_qy += phi_y.get(m,k) * qin->get(i,j,me,k);
                    }
                    qx_vals.set(m,me, tmp_qx );
                    qy_vals.set(m,me, tmp_qy );

                    // First non-zero terms start at third-order.
                    for( int  k=4; k <= kmax; k++   )
                    {
                        fxx_vals.set(m,me, fxx_vals.get(m,me) + phi_xx.get(m,k)*F->get(i,j,me,k) );
                        gxx_vals.set(m,me, gxx_vals.get(m,me) + phi_xx.get(m,k)*G->get(i,j,me,k) );

                        fxy_vals.set(m,me, fxy_vals.get(m,me) + phi_xy.get(m,k)*F->get(i,j,me,k) );
                        gxy_vals.set(m,me, gxy_vals.get(m,me) + phi_xy.get(m,k)*G->get(i,j,me,k) );

                        fyy_vals.set(m,me, fyy_vals.get(m,me) + phi_yy.get(m,k)*F->get(i,j,me,k) );
                        gyy_vals.set(m,me, gyy_vals.get(m,me) + phi_yy.get(m,k)*G->get(i,j,me,k) );
                    }

                }

                // ----------------------------------- //
                // Part I: Compute (f_x + g_y)_{,t}
                // ----------------------------------- //

                // Compute terms that get multiplied by \pd2{ f }{ q } and \pd2{ g }{ q }.
                dTensor2 fx_plus_gy_t( mpoints, meqn );
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

            // ------------------------------------------ //
            // 
            // Construct basis coefficients
            //
            // << Copied from integrate_on_current_cell >>
            //
            // ------------------------------------------ //
            for( int me=1; me <= meqn; me++ )
            for( int  k=1;  k <= kmax;  k++ )
            {
                double tmp1 = 0.;
                double tmp2 = 0.;
                for( int m=1; m <= mpoints; m++ )
                {
                    tmp1 += wgt_phi_transpose.get(k,m)*( 
                        dt_times_fdot.get(m, me ) + charlie_dt*f_tt.get(m,me) );
                    tmp2 += wgt_phi_transpose.get(k,m)*(
                        dt_times_gdot.get(m, me ) + charlie_dt*g_tt.get(m,me) );

//                  tmp1 += wgt_phi_transpose.get(k,m)*( dt_times_fdot.get(m, me ) );
//                  tmp2 += wgt_phi_transpose.get(k,m)*( dt_times_gdot.get(m, me ) );
                }
                F->set(i,j,me,k, alpha*F->get(i,j,me,k) + 0.25*tmp1 );
                G->set(i,j,me,k, alpha*G->get(i,j,me,k) + 0.25*tmp2 );
            }
            // ------------------------------------------ //

            // second order and higher terms
//          else  
//          {
//              // first order case, insert the factor of alpha back in here:
//              for( int me=1; me <= meqn; me++ )
//              for( int  k=1;  k <= kmax;  k++ )
//              {
//                  F->fetch(i,j,me,k) *= alpha;
//                  G->fetch(i,j,me,k) *= alpha;
//              }
//          }
        }
    }

}
