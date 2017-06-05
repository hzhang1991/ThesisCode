#include <stdlib.h>
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "edge_data.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"
#include "L2Project.h"
#include "Legendre2d.h"
#include "float.h" // for debugging

// Right-hand side of the Lax-Wendroff discretization:
//
//      ( q(t+dt) - q(t) )/dt = -F_{x} - G_{y}
//
// where the the flux is given by
//
//      F = F1 + dt/2 * F2 + dt^2/6 * F3
//
//        F1 =  f(q) + g(q)
//        F2 = -f'(q) * ( f(q)_x + g(q)_y ) - g'(q)*( f(q)_x + g(q)_y )
//        F3 =  f''(q) * f(q)_{x} * f(q)_{x} - f'(q)*f(q)_{xx} + ...
//
// This routine constructs the modified flux functions F and G using the
// Cauchy-Kowalewski procedure.
//
// The functions F and G are created by looking at time derivatives of q.  For
// 'second order' accuracy, we set alpha1 = 1, and beta1 = 1/2, and then the
// Taylor series looks like:
//
//    q^{n+1} = q^n + dt * (alpha1 q1_t + dt * beta1 q1_tt ).
//
// Which means that
//
//    ( q^{n+1} - q^n ) / dt =   (alpha1 q1_t + dt * beta1 q1_tt ) 
//
// Each time derivative leaves a single spatial derivative that can be pulled 
// out, so we write
//
//   q_{t}  = ( -f(q) )_x + ( -g(q) )_y, 
//   q_{tt} = ( f'(q) * ( f_x + g_y )_x + ( g'(q) * ( f_x + g_y )_y
//
// And then work with the weak formulation for F and G:
//
//    F: = ( alpha1*f - beta1*dt * ( f'(q)*( f_x+g_y ) )
//    G: = ( alpha1*g - beta1*dt * ( g'(q)*( f_x+g_y ) )
//
// The final recipe sets the following information in Lstar:
//
//    L(q) := -F_x - G_y.
//
// For the analogous RK time integrator, see ConstructL
//
void LaxWendroffTD(double dt, 
    double alpha1,  double beta1,      // parameteres used for each derivative.
    dTensorBC4& aux1, dTensorBC4& q1,    // set bndy values modifies these
    double alpha2,  double beta2,      // parameteres used for each derivative.
    dTensorBC4& aux2, dTensorBC4& q2,    // set bndy values modifies these
    dTensorBC4& Lstar, dTensorBC3& smax)
{

    if ( !dogParams.get_flux_term() )
    {  return;  }

    const edge_data& edgeData = Legendre2d::get_edgeData();
    const int space_order = dogParams.get_space_order();
    const int mx   = q1.getsize(1);
    const int my   = q1.getsize(2);
    const int meqn = q1.getsize(3);
    const int kmax = q1.getsize(4);
    const int mbc  = q1.getmbc();
    const int maux = aux1.getsize(3);

    // Flux values
    //
    // Space-order = number of quadrature points needed for 1D integration
    // along cell edges.
    //
    dTensorBC4 Fm(mx, my, meqn, space_order, mbc);
    dTensorBC4 Fp(mx, my, meqn, space_order, mbc);
    dTensorBC4 Gm(mx, my, meqn, space_order, mbc);
    dTensorBC4 Gp(mx, my, meqn, space_order, mbc);

    // Flux function
    void FluxFunc(const dTensor2& xpts,
                  const dTensor2& Q,
                  const dTensor2& Aux,
                  dTensor3& flux);

    // Jacobian of the flux function:
    void DFluxFunc(const dTensor2& xpts, 
                   const dTensor2& Q,
                   const dTensor2& Aux, 
                   dTensor4& Dflux );

    // Riemann solver that relies on the fact that we already have 
    // f(ql) and f(qr) already computed:
    double RiemannSolveLxW(const dTensor1& nvec,
        const dTensor1& xedge,
        const dTensor1& Ql,   const dTensor1& Qr,
        const dTensor1& Auxl, const dTensor1& Auxr,
        const dTensor1& ffl,  const dTensor1& ffr,
        dTensor1& Fl, dTensor1& Fr);

    void LstarExtra(const dTensorBC4*,
            const dTensorBC4*,
            dTensorBC4*);
    void ArtificialViscosity(const dTensorBC4* aux, 
            const dTensorBC4* q, 
            dTensorBC4* Lstar);

    // Grid information
    const double xlower = dogParamsCart2.get_xlow();
    const double ylower = dogParamsCart2.get_ylow();
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    // --------------------------------------------------------------------- //
    // Boundary data:
    // --------------------------------------------------------------------- //
    void SetBndValues(dTensorBC4& q, dTensorBC4& aux);
    SetBndValues(q1, aux1);
    SetBndValues(q2, aux2);
    // ---------------------------------------------------------

    // --------------------------------------------------------------------- //
    // Part 0: Compute the Lax-Wendroff "flux" function:
    //
    // Here, we include the extra information about time derivatives.
    // --------------------------------------------------------------------- //
    void L2ProjectLxW( double alpha, double beta_dt,
        const int istart, const int iend, 
        const int jstart, const int jend,
        const int QuadOrder,
        const int BasisOrder_qin,
        const int BasisOrder_auxin,
        const int BasisOrder_fout,    
        const dTensorBC4* qin,
        const dTensorBC4* auxin,
        dTensorBC4* F,
        dTensorBC4* G,
        void FluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& Aux, dTensor3& flux),
        void DFluxFunc (const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& aux, dTensor4& Dflux) );

    dTensorBC4 F(mx, my, meqn, kmax, mbc);
    dTensorBC4 G(mx, my, meqn, kmax, mbc);
    L2ProjectLxW( alpha1, dt*beta1, 
        1-mbc, mx+mbc, 1-mbc, my+mbc,
        space_order, space_order, space_order, space_order,
        &q1, &aux1, &F, &G, &FluxFunc, &DFluxFunc );

    // Add in contributions from the second term
    dTensorBC4 Ftmp(mx, my, meqn, kmax, mbc);
    dTensorBC4 Gtmp(mx, my, meqn, kmax, mbc);
    L2ProjectLxW( alpha2, dt*beta2, 
        1-mbc, mx+mbc, 1-mbc, my+mbc,
        space_order, space_order, space_order, space_order,
        &q2, &aux2, &Ftmp, &Gtmp, &FluxFunc, &DFluxFunc );

    // ---------------------------------------------------------------------- //
    // Add the result from the first projection onto the second:
    // ---------------------------------------------------------------------- //
    const int numel = Ftmp.numel();
assert_eq( numel, (mx+2*mbc)*(my+2*mbc)*meqn*kmax );
#pragma omp parallel for
    for(int v=0; v < numel; v++)
    {
//      double tmpf = F.vget(v);
//      double tmpg = G.vget(v);
        F.vfetch(v)  += Ftmp.vget(v);
        G.vfetch(v)  += Gtmp.vget(v);
//      assert_almost_eq( tmpf, F.vget(v) - Ftmp.vget(v) );
//      assert_almost_eq( tmpg, G.vget(v) - Gtmp.vget(v) );
    }
    // ---------------------------------------------------------------------- //

/*
if( fabs( alpha2 ) < 1e-14 && fabs( beta2 ) < 1e-14 )
{

    printf("we think alpha2 = beta2 = 0\n");
    printf("numel = %d\n", numel );

    for(int v=0; v < numel; v++)
    {
        assert_almost_eq( Ftmp.vget(v), 0. );
        assert_almost_eq( Gtmp.vget(v), 0. );
    }
}


    // Check to make sure L2ProjectLxW is a linear function:
    const double alpha = alpha1+alpha2;
    const double  beta = beta1+beta2;

    dTensorBC4    q( mx, my, meqn, kmax, mbc);
    dTensorBC4  aux( mx, my, maux, kmax, mbc);

    void CopyQ(const dTensorBC4& qin, dTensorBC4& qout);
    CopyQ(   q1,   q );
    CopyQ( aux1, aux );

    dTensorBC4 Fsum( mx, my, meqn, kmax, mbc );
    dTensorBC4 Gsum( mx, my, meqn, kmax, mbc );
    L2ProjectLxW( alpha, beta*dt,
        1-mbc, mx+mbc, 1-mbc, my+mbc,
        space_order, space_order, space_order, space_order,
        &q, &aux, &Fsum, &Gsum, &FluxFunc, &DFluxFunc );

    // ---------------------------------------------------------------------- //
    printf("alpha,   beta = %f, %f\n", alpha,beta);
    printf("alpha1, beta1 = %f, %f\n", alpha1,beta1);
    printf("alpha2, beta2 = %f, %f\n", alpha2,beta2);
    for(int i=1; i <= mx; i++ )
    for(int j=1; j <= my; j++ )
    for(int k=1; k <= kmax; k++ )
    {
        assert_almost_eq( Fsum.get(i,j,1,k), F.get(i,j,1,k) );
        assert_almost_eq( Gsum.get(i,j,1,k), G.get(i,j,1,k) );
    }
*/


    // ---------------------------------------------------------
    // Part I: compute source term
    // --------------------------------------------------------- 
    if( dogParams.get_source_term() > 0 )
    {
        // eprintf("error: have not implemented source term for LxW solver.");
        printf("Source term has not been implemented for LxW solver.  Terminating program.");
        exit(1);
    }
    Lstar.setall( 0. );

    // ---------------------------------------------------------
    // Part II: compute inter-element interaction fluxes
    //
    //   N = int( F(q,x,t) * phi_x, x ) / dA
    //
    // ---------------------------------------------------------

    // 1-direction: loop over interior edges and solve Riemann problems
    dTensor1 nvec(2);
    nvec.set(1, 1.0e0 );
    nvec.set(2, 0.0e0 );

#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc); i++)
    {

        dTensor1 Ql(meqn),   Qr(meqn);
        dTensor1 ffl(meqn),  ffr(meqn);
        dTensor1 Fl(meqn),   Fr(meqn);
        dTensor1 DFl(meqn),  DFr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);

        dTensor1 xedge(2);

        for (int j=(2-mbc); j<=(my+mbc-1); j++)
        {
            // ell indexes Riemann point along the edge
            for (int ell=1; ell<=space_order; ell++)
            {
                // Riemann data - q and f (from basis functions/q)
                for (int m=1; m<=meqn; m++)
                {
                    Ql.set (m, 0.0 );
                    Qr.set (m, 0.0 );
                    ffl.set(m, 0.0 );
                    ffr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        // phi_xl( xi=1.0, eta ), phi_xr( xi=-1.0, eta )
                        Ql.fetch(m)  += edgeData.phi_xl->get(ell,k)*q1.get(i-1, j, m, k );
                        Qr.fetch(m)  += edgeData.phi_xr->get(ell,k)*q1.get(i,   j, m, k );
                        ffl.fetch(m) += edgeData.phi_xl->get(ell,k)*F.get(i-1,  j, m, k );
                        ffr.fetch(m) += edgeData.phi_xr->get(ell,k)*F.get(i,    j, m, k );
                    }

                }

                // Riemann data - aux
                for (int m=1; m<=maux; m++)
                {
                    Auxl.set(m, 0.0 );
                    Auxr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Auxl.fetch(m) += edgeData.phi_xl->get(ell,k)*aux1.get(i-1, j, m, k);
                        Auxr.fetch(m) += edgeData.phi_xr->get(ell,k)*aux1.get(i,   j, m, k);
                    }
                }

                // Solve Riemann problem
                xedge.set(1, xlower + (double(i)-1.0)*dx );
                xedge.set(2, ylower + (double(j)-0.5)*dy );

                const double smax_edge = RiemannSolveLxW(
                    nvec, xedge, Ql, Qr, Auxl, Auxr, ffl, ffr, Fl, Fr);

                smax.set(i-1, j, 1, Max(dy*smax_edge,smax.get(i-1, j, 1)) );
                smax.set(i,   j, 1, Max(dy*smax_edge,smax.get(i,   j, 1)) );

                // Construct fluxes
                for (int m=1; m<=meqn; m++)
                {
                    Fm.set(i  , j, m, ell,  Fr.get(m) );
                    Fp.set(i-1, j, m, ell,  Fl.get(m) );
                }
            }
        }
    }


    // 2-direction: loop over interior edges and solve Riemann problems
    nvec.set(1, 0.0e0 );
    nvec.set(2, 1.0e0 );

#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    {
        dTensor1  Ql(meqn),   Qr(meqn);
        dTensor1  Fl(meqn),   Fr(meqn);
        dTensor1 ffl(meqn),  ffr(meqn);
        dTensor1 Auxl(maux),Auxr(maux);
        dTensor1 xedge(2);

        for (int j=(2-mbc); j<=(my+mbc); j++)
        for (int ell=1; ell<=space_order; ell++)
        {
            // Riemann data - q
            for (int m=1; m<=meqn; m++)
            {

                Ql.set  (m, 0.0 );
                Qr.set  (m, 0.0 );
                ffl.set (m, 0.0 );
                ffr.set (m, 0.0 );

                for (int k=1; k<=kmax; k++)
                {
                    Ql.fetch(m)  += edgeData.phi_yl->get(ell, k)*q1.get(i, j-1, m, k );
                    Qr.fetch(m)  += edgeData.phi_yr->get(ell, k)*q1.get(i, j,   m, k );
                    ffl.fetch(m) += edgeData.phi_yl->get(ell, k)*G.get(i, j-1, m, k );
                    ffr.fetch(m) += edgeData.phi_yr->get(ell, k)*G.get(i,   j, m, k );
                }
            }

            // Riemann data - aux
            for (int m=1; m<=maux; m++)
            {
                Auxl.set(m, 0.0 );
                Auxr.set(m, 0.0 );

                for (int k=1; k<=kmax; k++)
                {
                    Auxl.fetch(m) += edgeData.phi_yl->get(ell,k)*aux1.get(i,j-1,m,k);
                    Auxr.fetch(m) += edgeData.phi_yr->get(ell,k)*aux1.get(i,j,m,k);
                }
            }

            // Solve Riemann problem
            xedge.set(1, xlower + (double(i)-0.5)*dx );
            xedge.set(2, ylower + (double(j)-1.0)*dy );

            const double smax_edge = RiemannSolveLxW(
                nvec, xedge, Ql, Qr, Auxl, Auxr, ffl, ffr, Fl, Fr);

            smax.set(i, j-1, 2,  Max(dx*smax_edge, smax.get(i, j-1, 2)) );
            smax.set(i, j,   2,  Max(dx*smax_edge, smax.get(i, j,   2)) );

            // Construct fluxes
            for (int m=1; m<=meqn; m++)
            {
                Gm.set(i, j,   m, ell, Fr.get(m) );
                Gp.set(i, j-1, m, ell, Fl.get(m) );
            }
        }
    }

    // Compute ``flux differences'' dF and dG    
    const double half_dx_inv = 0.5/dx;
    const double half_dy_inv = 0.5/dy;
    const int mlength = Lstar.getsize(3);   assert_eq( meqn, mlength );

    // Use the four values, Gm, Gp, Fm, Fp to construct the boundary integral:
#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)    
    for (int j=(2-mbc); j<=(my+mbc-1); j++)
    for (int m=1; m<=mlength; m++)
    for (int k=1; k<=kmax; k++)
    {
        // 1-direction: dF
        double F1 = 0.0;
        double F2 = 0.0;
        for (int ell=1; ell<=space_order; ell++)
        {
            F1 = F1 + edgeData.wght_phi_xr->get(ell,k)*Fm.get(i,j,m,ell);
            F2 = F2 + edgeData.wght_phi_xl->get(ell,k)*Fp.get(i,j,m,ell);
        }

        // 2-direction: dG
        double G1 = 0.0;
        double G2 = 0.0;
        for (int ell=1; ell<=space_order; ell++)
        {
            G1 = G1 + edgeData.wght_phi_yr->get(ell,k)*Gm.get(i,j,m,ell);
            G2 = G2 + edgeData.wght_phi_yl->get(ell,k)*Gp.get(i,j,m,ell);
        }

        Lstar.fetch(i,j,m,k) -= (half_dx_inv*(F2-F1) + half_dy_inv*(G2-G1));

    }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part III: compute intra-element contributions
    // ---------------------------------------------------------
    // No need to call this if first-order in space
    if( dogParams.get_space_order() > 1 )
    {
        void L2ProjectGradAddLegendre(const int istart, 
                const int iend, 
                const int jstart, 
                const int jend,
                const int QuadOrder, 
                const dTensorBC4* F, 
                const dTensorBC4* G, 
                dTensorBC4* fout );
        L2ProjectGradAddLegendre( 1-mbc, mx+mbc, 1-mbc, my+mbc,
            space_order, &F, &G, &Lstar );

    }
    // ---------------------------------------------------------  

    // ---------------------------------------------------------
    // Part IV: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    LstarExtra(&q1, &aux1, &Lstar);
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part V: artificial viscosity limiter
    // ---------------------------------------------------------  
    if (dogParams.get_space_order()>1  &&
            dogParams.using_viscosity_limiter())
    {  ArtificialViscosity(&aux1, &q1, &Lstar);  }
    // ---------------------------------------------------------

}
