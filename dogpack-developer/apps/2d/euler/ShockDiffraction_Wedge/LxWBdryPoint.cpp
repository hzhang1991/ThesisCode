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
void LxWBdryPoint( const int mterms,
        const double alpha, const double beta_dt, const double charlie_dt,
        const int QuadOrder,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,    
        dTensor1 qvals,
        dTensor1 auxvals,
        dTensor1* F,
        dTensor1* G,
        dTensor1* qtbdry,
        dTensor1* qttbdry,
        dTensor1* qtttbdry,
        void FluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& Aux, dTensor3& flux),
        void DFluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& aux, dTensor4& Dflux),
        void D2FluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& aux, dTensor5& D2flux) )
{    
void FluxFuncPoint(const dTensor2& xpts,
        const dTensor1& Q,
        const dTensor1& Aux,
        dTensor1& xflux,dTensor1& yflux);
void DFluxFuncPoint(const dTensor2& xpts,
        const dTensor1& Q,
        const dTensor1& Aux,
        dTensor4& Dflux);
void D2FluxFuncPoint(const dTensor2& xpts,
                const dTensor1& Q,
                const dTensor1& Aux,
                dTensor5& D2flux);

    if( fabs( alpha ) < 1e-14 && fabs( beta_dt ) < 1e-14 )
    {
        F->setall(0.);
        G->setall(0.);
        return;
    }

    // qin variable
    const int     meqn = qvals.getsize();

    // auxin variable
    const int       maux = auxvals.getsize();

    // F variables

    // mbc


    // Grid information:
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    // number of quadrature points
    const int mpoints = 1;

    // ---------------------------------------------------------------------- //
    // evaluate the basis functions, phi and its derivatives at all quadrature
    // points
    // ---------------------------------------------------------------------- //

    // ---------------------------------------------------------------------- //
    // For efficiency compute weight*phi and then take transpose
//  dTensor2 wgt_phi_x_transpose  (kmax, mpoints);
//  dTensor2 wgt_phi_y_transpose  (kmax, mpoints);
//  dTensor2 wgt_phi_xx_transpose (kmax, mpoints);
//  dTensor2 wgt_phi_xy_transpose (kmax, mpoints);
//  dTensor2 wgt_phi_yy_transpose (kmax, mpoints);
    // ---------------------------------------------------------------------- //

    // ------------------------------------------------------------- //
    // Loop over every grid cell indexed by user supplied parameters //
    // described by istart...iend, jstart...jend                     // 
    // ------------------------------------------------------------- //
        // Local storage for q, aux and xpts (all passed into user supplied
        // function)
        dTensor2    xpts(mpoints, 2);
        xpts.set(1,1,0.0);xpts.set(1,2,0.0);

        // local storage for Flux function its Jacobian, and the Hessian:
        dTensor1    xfvals(meqn);  // flux function (vector)
        dTensor1    yfvals(meqn);  // flux function (vector)
        dTensor4        A(mpoints,       meqn, meqn, 2);  // Jacobian of flux
        dTensor5        H(mpoints, meqn, meqn, meqn, 2);  // Hessian of flux
  

            //find center of current cell

            // ------------------------------------------------------------- //
            //
            // Part I:
            //
            // Project the flux function onto the basis 
            // functions.  This is the term of order O( 1 ) in the
            // "time-averaged" Taylor expansion of f and g.
            //
            // ------------------------------------------------------------- //

            // ------------------------------------------------------------- //
            // Project the flux function onto the basis functions.  This is
            // the term of order O( 1 ) in the Taylor expansions of f_x + g_y.
            // ------------------------------------------------------------- //

            // Evaluate the Flux function at Gaussian quadrature point
            FluxFuncPoint( xpts, qvals, auxvals, xfvals, yfvals );

            // copy onto smaller arrays (so we can call
            // integrate_on_current_cell):
            //
            // This is a bit of overkill, but we need to do something to
            // extra f, without the alpha in there.
            dTensor2 f( mpoints, meqn );
            dTensor2 g( mpoints, meqn );
            for( int me=1; me <= meqn ; me++ )
            {
                int m=1;
                f.set(m, me, xfvals.get(me) );
                g.set(m, me, yfvals.get(me) );
            }

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

            // Call user-supplied Jacobian to set f'(q) and g'(q):
            DFluxFuncPoint( xpts, qvals, auxvals, A );

            // place-holders for point values of
            // f'(q)( fx + gy ) and g'(q)( fx + gy ):
            dTensor2 dt_times_fdot( mpoints, meqn );
            dTensor2 dt_times_gdot( mpoints, meqn );

            // Compute point values for f'(q) * (fx+gy) and g'(q) * (fx+gy):
            for( int m1=1; m1 <= meqn; m1++ )
            {   int m=1;
                double tmp1 = 0.;
                double tmp2 = 0.;
                for( int m2=1; m2 <= meqn; m2++ )
                {
                    tmp1 += A.get(m, m1, m2, 1 ) * qtbdry->get(m2);
                    tmp2 += A.get(m, m1, m2, 2 ) * qtbdry->get(m2);
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
            dTensor2 fx_plus_gy_t( mpoints, meqn );
            if( mterms > 2 )
            {

                // Construct the Hessian at each (quadrature) point
                D2FluxFuncPoint( xpts, qvals, auxvals, H );

                // ----------------------------------- //
                // Part I: Compute 
                //      f'(q) * fx_plus_gy_t and 
                //      g'(q) * fx_plus_gy_t
                // ----------------------------------- //

                // Add in the third term that gets multiplied by A:
                for( int m1=1; m1 <= meqn; m1++ )
                {
                    int m=1;
                    double tmp1 = 0.;
                    double tmp2 = 0.;
                    for( int m2=1; m2 <= meqn; m2++ )
                    {
                        tmp1 += A.get(m,m1,m2,1)*qttbdry->get(m2);
                        tmp2 += A.get(m,m1,m2,2)*qttbdry->get(m2);
                    }
                    f_tt.set( m, m1, tmp1 );
                    g_tt.set( m, m1, tmp2 );
                }

                // ----------------------------------------------- //
                // Part III: Add in contributions from
                //      f''(q) * (fx_plus_gy, fx_plus_gy ) and 
                //      g''(q) * (fx_plus_gy, fx_plus_gy ).
                // ----------------------------------------------- //
                for( int me =1; me <= meqn; me++ )
                {
                    int m=1;
                    double tmp1 = 0.;
                    double tmp2 = 0.;

                    // Terms that get multiplied by the Hessian:
                    for( int m1=1; m1 <= meqn; m1++ )
                    for( int m2=1; m2 <= meqn; m2++ )
                    {
                        tmp1 += H.get(m,me,m1,m2,1)*qtbdry->get(m1)*qtbdry->get(m2);
                        tmp2 += H.get(m,me,m1,m2,2)*qtbdry->get(m1)*qtbdry->get(m2);
                    }

                    f_tt.set( m, me, f_tt.get(m,me) + tmp1 );
                    g_tt.set( m, me, g_tt.get(m,me) + tmp2 );
                }

            } // End of computing "third"-order terms

            for( int me =1; me <= meqn; me++ )
            {
                int m=1;
                F->set(me,alpha*F->get(me)+dt_times_fdot.get(m, me ) + charlie_dt*f_tt.get(m,me));
                G->set(me,alpha*G->get(me)+dt_times_gdot.get(m, me ) + charlie_dt*g_tt.get(m,me));
            }

}
