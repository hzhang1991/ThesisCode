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
#include "LaxWendroff.h"
#include <iostream>
using namespace std;

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
void LaxWendroff(double dt, 
    double alpha1,  double beta1,      // parameteres used for each derivative.
    dTensorBC4& aux, dTensorBC4& q,    // set bndy values modifies these
    dTensorBC4& Lstar, dTensorBC3& smax)
{

    if ( !dogParams.get_flux_term() )
    {  return;  }

    const edge_data& edgeData = Legendre2d::get_edgeData();
    const int space_order = dogParams.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    // Flux values
    //
    // Space-order = number of quadrature points needed for 1D integration
    // along cell edges.
    //
    dTensorBC4 Fm(mx, my, meqn, space_order, mbc);
    dTensorBC4 Fp(mx, my, meqn, space_order, mbc);
    dTensorBC4 Gm(mx, my, meqn, space_order, mbc);
    dTensorBC4 Gp(mx, my, meqn, space_order, mbc);

    dTensorBC3 FmI(mx, my, meqn, mbc);
    dTensorBC3 GmI(mx, my, meqn, mbc);
    dTensorBC3 qI(mx, my, meqn, mbc);
    dTensorBC3 auxI(mx, my, meqn, mbc);

    // Grid information
    const double xlower = dogParamsCart2.get_xlow();
    const double ylower = dogParamsCart2.get_ylow();
    const double dx = dogParamsCart2.get_dx();
    const double dy = dogParamsCart2.get_dy();

    // --------------------------------------------------------------------- //
    // Boundary data:
    //
    // TODO - this call is inconsistent with ConstructL, and the other
    // LaxWendroff routine.  This should be modified, but before doing so,
    // make sure SetBndValues is called before calling this routine. -DS
    // (12/16/2014).
    // --------------------------------------------------------------------- //
    void SetBndValues(dTensorBC4& q, dTensorBC4& aux);
    void SetBndValuesX(dTensorBC4& q, dTensorBC4& aux);
    void SetBndValuesY(dTensorBC4& q, dTensorBC4& aux);
    SetBndValues(q, aux);
    // ---------------------------------------------------------

    // --------------------------------------------------------------------- //
    // Part 0: Compute the Lax-Wendroff "flux" function:
    //
    // Here, we include the extra information about time derivatives.
    // --------------------------------------------------------------------- //
    dTensorBC4 F(mx, my, meqn, kmax, mbc);  F.setall(0.);
    dTensorBC4 G(mx, my, meqn, kmax, mbc);  G.setall(0.);
    L2ProjectLxW( 3, alpha1, beta1*dt, dt*dt/6.0,
        1-mbc, mx+mbc, 1-mbc, my+mbc,
        space_order, space_order, space_order,
        &q, &aux, &F, &G, &FluxFunc, &DFluxFunc, &D2FluxFunc );

    // ---------------------------------------------------------
    // Part I: compute source term
    // --------------------------------------------------------- 
    if( dogParams.get_source_term() > 0 )
    {
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

SetBndValuesX(q, aux);
#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc); i++)
    {

        // Information for Riemann states
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
                        Ql.fetch(m)  += edgeData.phi_xl->get(ell,k)*q.get(i-1, j, m, k );
                        Qr.fetch(m)  += edgeData.phi_xr->get(ell,k)*q.get(i,   j, m, k );
                        ffl.fetch(m) += edgeData.phi_xl->get(ell,k)*F.get(i-1, j, m, k );
                        ffr.fetch(m) += edgeData.phi_xr->get(ell,k)*F.get(i,   j, m, k );
                    }

                }

                // Riemann data - aux
                for (int m=1; m<=maux; m++)
                {
                    Auxl.set(m, 0.0 );
                    Auxr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Auxl.fetch(m) += edgeData.phi_xl->get(ell,k)*aux.get(i-1, j, m, k);
                        Auxr.fetch(m) += edgeData.phi_xr->get(ell,k)*aux.get(i,   j, m, k);
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

SetBndValuesY(q, aux);
#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    {

        // Information for Riemann states
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
                    Ql.fetch(m)  += edgeData.phi_yl->get(ell, k)*q.get(i, j-1, m, k );
                    Qr.fetch(m)  += edgeData.phi_yr->get(ell, k)*q.get(i, j,   m, k );
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
                    Auxl.fetch(m) += edgeData.phi_yl->get(ell,k)*aux.get(i,j-1,m,k);
                    Auxr.fetch(m) += edgeData.phi_yr->get(ell,k)*aux.get(i,j,m,k);
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

    // Compute Fluxes used for averages ... We will limit these fluxes to guarentee positive qunatity averages.
#pragma omp parallel for
    for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    for (int j=(2-mbc); j<=(my+mbc-1); j++)
    {
      for (int m=1; m<=meqn; m++)
    {

        // In case the MPP limiters are turned on, we compute the cell
        // averages separately from the higher-order terms.
        const int k=1;

        // 1-direction: dF
        double F1 = 0.0;
        double F2 = 0.0;
        for (int ell=1; ell<=space_order; ell++)
        {
            F1 = F1 + edgeData.wght_phi_xr->get(ell,k)*Fm.get(i,j,m,ell);
            F2 = F2 + edgeData.wght_phi_xl->get(ell,k)*Fm.get(i+1,j,m,ell);
        }

        // 2-direction: dG
        double G1 = 0.0;
        double G2 = 0.0;
        for (int ell=1; ell<=space_order; ell++)
        {
            G1 = G1 + edgeData.wght_phi_yr->get(ell,k)*Gm.get(i,j,m,ell);
            G2 = G2 + edgeData.wght_phi_yr->get(ell,k)*Gm.get(i,j+1,m,ell);
        }
        FmI.set (i,j,m, 0.5*F1);
        GmI.set (i,j,m, 0.5*G1);

        FmI.set (i+1,j,m, 0.5*F2);
        GmI.set (i,j+1,m, 0.5*G2);
        qI.set  (i,j,m, q.get(i,j,m,1));
     }
     for(int m=1;m<=maux;m++)
        {auxI.set(i,j,m, aux.get(i,j,m,1));}

    }

    // Apply flux limiter to that gaurantees positive averages.
    //
    // The default behaviour of this function is to do nothing.

    if( dogParams.using_positive_limiter() )
    { 
        void ApplyPosMPPLimiter( double dt, dTensorBC3& smax,
        dTensorBC3& aux, dTensorBC3& q, 
        dTensorBC3& fHat, dTensorBC3& gHat );
        ApplyPosMPPLimiter( dt, smax, auxI, qI, FmI, GmI ); 
    }

    // Use the four values, Gm, Gp, Fm, Fp to construct the boundary integral:
#pragma omp parallel for
    //for (int i=(2-mbc); i<=(mx+mbc-1); i++)    
    //for (int j=(2-mbc); j<=(my+mbc-1); j++)
    for (int i=(1); i<=(mx); i++)    
    for (int j=(1); j<=(my); j++)
    for (int m=1; m<=mlength; m++)
    {

        // TODO - this adds an extra loop for the non-positivity preserving code ...
        // Moreover, this is not the routine that DogSolveLxW calls.
        // -DS (12/16/2014).

        // Deal with the cell averages first
        
        int k = 1;

        // 1-direction: dF
        double F1 = 0.0;
        double F2 = 0.0;
        /*for (int ell=1; ell<=space_order; ell++)
        {
            F1 = F1 + edgeData.wght_phi_xr->get(ell,k)*FmI.get(i,j,m);
            F2 = F2 + edgeData.wght_phi_xl->get(ell,k)*FmI.get(i+1,j,m);
        }*/
        F1=2.0*FmI.get(i,j,m);
        F2=2.0*FmI.get(i+1,j,m);
        // 2-direction: dG
        double G1 = 0.0;
        double G2 = 0.0;
        /*for (int ell=1; ell<=space_order; ell++)
        {
            G1 = G1 + edgeData.wght_phi_yr->get(ell,k)*GmI.get(i,j,m);
            G2 = G2 + edgeData.wght_phi_yl->get(ell,k)*GmI.get(i,j+1,m);
        }*/
        G1=2.0*GmI.get(i,j,m);
        G2=2.0*GmI.get(i,j+1,m);
        Lstar.fetch(i,j,m,k) -= (half_dx_inv*(F2-F1) + half_dy_inv*(G2-G1));
        

    for (int k=2; k<=kmax; k++)
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
    }

    // ---------------------------------------------------------------------- //
    // Part III: compute intra-element contributions
    // ---------------------------------------------------------------------- //
    // No need to call this if first-order in space
    if(dogParams.get_space_order()>1)
    {
        L2ProjectGradAddLegendre( 1-mbc, mx+mbc, 1-mbc, my+mbc,
            space_order, &F, &G, &Lstar );
    }
    // ---------------------------------------------------------------------- //
    // ---------------------------------------------------------  

    // ---------------------------------------------------------------------- //
    // Part IV: add extra contributions to Lstar
    // ---------------------------------------------------------------------- //
    // Call LstarExtra
    LstarExtra(&q,&aux,&Lstar);
    // ---------------------------------------------------------------------- //

    // ---------------------------------------------------------------------- //
    // Part V: artificial viscosity limiter
    // ---------------------------------------------------------------------- //
    if (dogParams.get_space_order()>1  && dogParams.using_viscosity_limiter())
    {  ArtificialViscosity(&aux,&q,&Lstar);  }
    // ---------------------------------------------------------

}
