//#undef CHECK_BOUNDS // this file is okay, so omit tensor bounds check
#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "edge_data.h"
#include "DogParams.h"
#include "DogParamsCart2.h"
#include "L2Project.h"
#include "RiemannSolve.h"
#include "Legendre2d.h"
#include "float.h" // for debugging
#include <iostream>
using namespace std;

// ------------------------------------------------------------------------- //
// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t = RHS =  -[ f(q,x,y,t)_x + g(q,x,y,t)_y ] + Psi(q,x,y,t)
//
// It is expected that the user sets the correct boundary conditions before
// calling this routine.
//
// Input:
// ------
//
//     aux(1:mx, 1:my, 1:maux, 1:kmax ) - auxiliary function
//       q(1:mx, 1:my, 1:meqn, 1:kmax ) - vector of conserved variables
//
// Returns:
// --------
//
//   Lstar(1:mx, 1:my, 1:meqn, 1:kmax ) - right hand side of MOL formulation
//    smax(1:mx, 1:my, 1:2 )            - maximum speed observed 
//
// ------------------------------------------------------------------------- //
void ConstructL(const dTensorBC4& aux, 
        const dTensorBC4& q,
        dTensorBC4& Lstar, 
        dTensorBC3& smax)
{
void SetBndFlux(const dTensor1& nvec, const dTensor1& xedge,
                    dTensor1& Ql, dTensor1& Qr,
                    const dTensor1& Auxl, const dTensor1& Auxr,
                    int side);
    const edge_data& edgeData = Legendre2d::get_edgeData();
    const int space_order = dogParams.get_space_order();
    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    // Data for Riemann problems
    dTensorBC4 Fm(mx, my, meqn, space_order, mbc);
    dTensorBC4 Fp(mx, my, meqn, space_order, mbc);
    dTensorBC4 Gm(mx, my, meqn, space_order, mbc);
    dTensorBC4 Gp(mx, my, meqn, space_order, mbc);
    Fm.setall(0.);Fp.setall(0.);Gm.setall(0.);Gp.setall(0.);
    dTensorBC3 Smax1(mx, my, 4, mbc);
    Smax1.setall(0.);
    // If you need access elsewhere to one of these declarations
    // move it into DogSolverCart2.h rather than copying the
    // declaration. (By using a single declaration we can append
    // arguments with default values without modifying existing
    // calls.) -eaj
    //
    void FluxFunc(const dTensor2& xpts, 
            const dTensor2& Q, const dTensor2& Aux, dTensor3& flux);
    void LstarExtra(const dTensorBC4*, const dTensorBC4*, dTensorBC4*);
    void ArtificialViscosity(const dTensorBC4* aux, const dTensorBC4* q, 
            dTensorBC4* Lstar);

    // Grid information
    const double xlower = dogParamsCart2.get_xlow();
    const double ylower = dogParamsCart2.get_ylow();
    const double dx     = dogParamsCart2.get_dx();
    const double dy     = dogParamsCart2.get_dy();

    // ---------------------------------------------------------
    // Part I: compute source term
    // --------------------------------------------------------- 
    if ( dogParams.get_source_term()>0 )
    {        
        // Set source term on computational grid
        // Set values and apply L2-projection
        void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals,
                const dTensor2& auxvals, dTensor2& source);
        L2Project(1-mbc,mx+mbc,1-mbc,my+mbc,
                space_order,space_order,space_order,space_order,
                &q,&aux,&Lstar,&SourceTermFunc);
    }
    else
    {
        Lstar.setall(0.);
    }

    if ( !dogParams.get_flux_term() )
    {  return;  }
    // ---------------------------------------------------------

    const int istep = (mx / 5)+1;
    const int jstep = (my / 5)+1;
    // ---------------------------------------------------------
    // Part II: compute inter-element interaction fluxes
    // ---------------------------------------------------------

    // 1-direction: loop over interior edges and solve Riemann problems
    dTensor1 nvec(2);
    nvec.set(1, 1.0e0 );
    nvec.set(2, 0.0e0 );

#pragma omp parallel for
    for (int i=(2-1); i<=(mx+1); i++)
    //for (int i=(2-mbc); i<=(mx+mbc); i++)
    {
        dTensor1 Ql(meqn),Qr(meqn);
        dTensor1 Auxl(maux),Auxr(maux);
        dTensor1 Fl(meqn),Fr(meqn);
        dTensor1 xedge(2);

        RiemannSolver riemannSolver(meqn,maux);
        //for (int j=(2-mbc); j<=(my+mbc-1); j++)
        for (int j=(1); j<=(my+1); j++)
        {

           if(i<=istep || j>=jstep)
           {

            // "ell" indexes Riemann point along the edge
            for (int ell=1; ell<=space_order; ell++)
            {
                // Riemann data - q
                for (int m=1; m<=meqn; m++)
                {
                    Ql.set(m, 0.0 );
                    Qr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Ql.fetch(m) += edgeData.phi_xl->get(ell,k)*q.get(i-1,j,m,k);
                        Qr.fetch(m) += edgeData.phi_xr->get(ell,k)*q.get(i,j,m,k);
                    }
                }

                // Riemann data - aux
                for (int m=1; m<=maux; m++)
                {
                    Auxl.set(m, 0.0 );
                    Auxr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Auxl.fetch(m) += edgeData.phi_xl->get(ell,k)*aux.get(i-1,j,m,k);
                        Auxr.fetch(m) += edgeData.phi_xr->get(ell,k)*aux.get(i,j,m,k);
                    }
                }
                // Solve Riemann problem
                xedge.set(1, xlower + (double(i)-1.0)*dx );
                xedge.set(2, ylower + (double(j)-0.5)*dy );
                int iside=0;
                if(i==1){iside=1;SetBndFlux(nvec,xedge,Ql,Qr,Auxl,Auxr,iside);}
                else if (i==istep && j < jstep){

                /*
                for (int m=1; m<=meqn; m++)
                {
                    //Ql.set(m, q.get(i-1,j,m,1) );
                    Qr.set(m, q.get(i-1,j,m,1) );
                } 
                Qr.set(2,-Qr.get(2));*/
                iside=2;SetBndFlux(nvec,xedge,Ql,Qr,Auxl,Auxr,iside);
                }
                else if (j >=jstep && i==mx+1){
                /*
                for (int m=1; m<=meqn; m++)
                {
                    //Ql.set(m, q.get(i-1,j,m,1) );
                    Qr.set(m, q.get(i-1,j,m,1) );
                }*/   
                iside=2;SetBndFlux(nvec,xedge,Ql,Qr,Auxl,Auxr,iside);
                /* for (int m=1; m<=meqn; m++)
                {
                    Ql.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Ql.fetch(m) += edgeData.phi_xl->get(ell,k)*q.get(i-1,j,m,k);
                    }
                }*/

                }



                const double smax_edge
                    = riemannSolver.solve(nvec,xedge,Ql,Qr,Auxl,Auxr,Fl,Fr);

                // TODO: See comments in all ConstructL methods.  This *should* be a 
                // bug when OpenMP loops are turned on. (-DS)
                //
                // If we reverse the role of i and j in the above for loop,
                // then it will be OK, but the tensors don't run as fast when
                // that's the case.
                //
                //smax.set(i-1, j, 1, Max(dy*smax_edge, smax.get(i-1, j, 1) ) );
                //smax.set(i,   j, 1, Max(dy*smax_edge, smax.get(i,   j, 1) ) );

                Smax1.set(i,j,1,Max(dy*smax_edge, smax.get(i,   j, 1)));
                Smax1.set(i,j,2,Max(dy*smax_edge, smax.get(i-1, j, 1)));
              

                // Construct fluxes
                for (int m=1; m<=meqn; m++)
                {
                    Fm.set(i  ,j,m,ell, Fr.get(m) );
                    Fp.set(i-1,j,m,ell, Fl.get(m) );
                }
            }
        }
     }
    }


    // ---------------------------------------------------------------------- //
    // 2-direction: loop over interior edges and solve Riemann problems
    // ---------------------------------------------------------------------- //
    nvec.set(1, 0.0e0 );
    nvec.set(2, 1.0e0 );

#pragma omp parallel for
    //for (int i=(2-mbc); i<=(mx+mbc-1); i++)
    for (int i=1; i<=(mx+1); i++)
    {

        dTensor1 Ql(meqn),Qr(meqn);
        dTensor1 Auxl(maux),Auxr(maux);
        dTensor1 xedge(2);
        dTensor1 Fl(meqn),Fr(meqn);

        //for (int j=(2-mbc); j<=(my+mbc); j++)
        for (int j=1; j<=(my+1); j++)
         {
           if(i<=istep || j>=jstep)
           {

             for (int ell=1; ell<=space_order; ell++)
            {
                // Riemann data - q
                for (int m=1; m<=meqn; m++)
                {
                    Ql.set(m, 0.0 );
                    Qr.set(m, 0.0 );

                    for (int k=1; k<=kmax; k++)
                    {
                        Ql.fetch(m) += edgeData.phi_yl->get(ell,k)*q.get(i,j-1,m,k);
                        Qr.fetch(m) += edgeData.phi_yr->get(ell,k)*q.get(i,j,m,k);
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
                int iside=0;

                if(j==1){
                /*for (int m=1; m<=meqn; m++)
                {
                    Ql.set(m, q.get(i,j,m,1));
                 //   Qr.set(m, q.get(i,j,m,1));
                }
                Ql.set(3,-Ql.get(3));
                */
                iside=1;SetBndFlux(nvec,xedge,Ql,Qr,Auxl,Auxr,iside);
                /*for (int m=1; m<=meqn; m++)
                {
                    Qr.set(m, 0.0 );
                    for (int k=1; k<=kmax; k++)
                    {
                        Qr.fetch(m) += edgeData.phi_yr->get(ell,k)*q.get(i,j,m,k);
                    }
                }*/

                }
                else if (j==jstep && i >=istep){
                //nvec.set(1, 0.0e0 );
                //nvec.set(2, -1.0e0 );
                /*for (int m=1; m<=meqn; m++)
                {
                    Ql.set(m, q.get(i,j,m,1));
                 //   Qr.set(m, q.get(i,j,m,1));
                }
                Ql.set(3,-Ql.get(3));*/
                iside=1;SetBndFlux(nvec,xedge,Ql,Qr,Auxl,Auxr,iside);
                //nvec.set(1, 0.0e0 );
                //nvec.set(2, 1.0e0 );
                }
                else if (j == my+1){
                //nvec.set(1, 0.0e0 );
                //nvec.set(2, -1.0e0 );                
                /*for (int m=1; m<=meqn; m++)
                {
                    Qr.set(m, q.get(i,j-1,m,1));
                 //   Qr.set(m, q.get(i,j,m,1));
                }
                Qr.set(3,-Qr.get(3));
                */
                iside=2;SetBndFlux(nvec,xedge,Ql,Qr,Auxl,Auxr,iside);
                //nvec.set(1, 0.0e0 );
                //nvec.set(2, 1.0e0 );
                }


                // TODO: why is this different than what was used earlier in
                // the code? (-DS)
                const double smax_edge = RiemannSolve(nvec,xedge,Ql,Qr,Auxl,Auxr,Fl,Fr);

                // TODO: this looks like a potential problem for the upper
                // pragma statement.  See above comment. (-DS)
                //
                //smax.set(i, j-1, 2, Max(dx*smax_edge,smax.get(i, j-1,2) ) );
                //smax.set(i, j,   2, Max(dx*smax_edge,smax.get(i, j,  2) ) );


                Smax1.set(i,j,3,Max(dy*smax_edge, smax.get(i,   j, 2)));
                Smax1.set(i,j,4,Max(dy*smax_edge, smax.get(i, j-1, 2)));

                // Construct fluxes
                for (int m=1; m<=meqn; m++)
                {
                    Gm.set(i,j,  m,ell, Fr.get(m) );
                    Gp.set(i,j-1,m,ell, Fl.get(m) );
                }
            }
         }
        }
    }

    // ---------------------------------------------------------------------- //
    // Modify boundary fluxes if necessary
    // ---------------------------------------------------------------------- //
    void SetBndFluxes(const dTensorBC4& q, 
            const dTensorBC4& aux,
            dTensorBC4& Fm,
            dTensorBC4& Fp,
            dTensorBC4& Gm,
            dTensorBC4& Gp);
    SetBndFluxes(q,aux,Fm,Fp,Gm,Gp);

    // ---------------------------------------------------------------------- //
    // Compute ``flux differences'' dF and dG    
    // ---------------------------------------------------------------------- //
    const double half_dx_inv = 0.5/dx;
    const double half_dy_inv = 0.5/dy;
    const int mlength = Lstar.getsize(3);

    //for (int i=(2-mbc); i<=(mx+mbc-1); i++)    
    //for (int j=(2-mbc); j<=(my+mbc-1); j++)
#pragma omp parallel for
    for (int i=1; i<=(mx+1); i++)    
    for (int j=1; j<=(my+1); j++)
    {
    smax.set(i, j,   1, Max(Smax1.get(i,j,1),Smax1.get(i+1,j,2) ) );
    smax.set(i, j,   2, Max(Smax1.get(i,j,3),Smax1.get(i,j+1,4) ) );
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
    }
    // ---------------------------------------------------------------------- //

    // ---------------------------------------------------------------------- //
    // Part III: compute intra-element contributions
    // ---------------------------------------------------------------------- //
    // No need to call this if first-order in space
    if(dogParams.get_space_order()>1)
    {
        L2ProjectGradAdd(1-mbc,mx+mbc,1-mbc,my+mbc,
                space_order,space_order,space_order,space_order,
                &q,&aux,&Lstar,&FluxFunc);
    }

    // ---------------------------------------------------------------------- //
    // Part IV: add extra contributions to Lstar, if any.
    // ---------------------------------------------------------------------- //
    LstarExtra(&q, &aux, &Lstar);

    // ---------------------------------------------------------------------- //
    // Part V: artificial viscosity limiter
    // ---------------------------------------------------------------------- //
    if (dogParams.get_space_order()>1  && dogParams.using_viscosity_limiter())
    {  ArtificialViscosity(&aux, &q, &Lstar);  }

}
