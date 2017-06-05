#include "dogdefs.h"
#include "dog_math.h"
#include "mesh.h"
#include "edge_data_Unst.h"
#include "DogParams.h"
#include "stdlib.h"
#include "LaxWendroff_Unst.h"
#include <iostream>
using namespace std;


// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t = -( f(q,x,y,t)_x + g(q,x,y,t)_y ) + Psi(q,x,y,t)
//
void LaxWendroff_Unst(double dt,
    const mesh& Mesh, const edge_data_Unst& EdgeData,
    dTensor3& aux,                  // SetBndValues modifies ghost cells
    dTensor3& q,                    // SetBndValues modifies ghost cells
    dTensor3& Lstar, dTensor1& smax)
{

    const int NumElems      = Mesh.get_NumElems();
    const int NumPhysElems  = Mesh.get_NumPhysElems();
    const int NumEdges      = Mesh.get_NumEdges();
    const int meqn          = q.getsize(2);
    const int kmax          = q.getsize(3);
    const int maux          = aux.getsize(2);
    const int space_order   = dogParams.get_space_order();
    dTensor3 EdgeFluxIntegral(NumElems,meqn,kmax);
    dTensor3 ElemFluxIntegral(NumElems,meqn,kmax);
    dTensor3              Psi(NumElems,meqn,kmax);

    // ---------------------------------------------------------
    // Boundary Conditions
    SetBndValues_Unst(Mesh, &q, &aux);
    // ---------------------------------------------------------

    // --------------------------------------------------------------------- //
    // Part 0: Compute the Lax-Wendroff "flux" function:
    //
    // Here, we include the extra information about time derivatives.
    // --------------------------------------------------------------------- //
    dTensor3 F(NumElems, meqn, kmax );  F.setall(0.);
    dTensor3 G(NumElems, meqn, kmax );  G.setall(0.);

    dTensor3 FeI(NumElems, meqn, 3 );  FeI.setall(0.);

    //Make this something we compute once and store or something

    /*
    L2ProjectLxW_Unst( dogParams.get_time_order(), 1.0, 0.5*dt, dt*dt/6.0, 1, NumElems,
        space_order, space_order, space_order, space_order, Mesh,
        &q, &aux, &F, &G, &FluxFunc, &DFluxFunc, &D2FluxFunc );*/

    L2ProjectLxW_Unst( 3, 1.0, 0.5*dt, dt*dt/6.0, 1, NumElems,
        3, 3, 3, 3, Mesh,
        &q, &aux, &F, &G, &FluxFunc, &DFluxFunc, &D2FluxFunc );


    // ---------------------------------------------------------
    // Part I: compute source term
    // --------------------------------------------------------- 
    if ( dogParams.get_source_term()>0 )
    {        
        // eprintf("error: have not implemented source term for LxW solver.");
        printf("Source term has not been implemented for LxW solver.  Terminating program.");
        exit(1);
    }
    Lstar.setall(0.);
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part II: compute flux integral on element edges
    // ---------------------------------------------------------

    // Loop over all interior edges
    EdgeFluxIntegral.setall(0.);
    ElemFluxIntegral.setall(0.);

#pragma omp parallel for
    // Loop over all interior edges
    for (int i=1; i<=NumEdges; i++)
    {
        // Edge coordinates
        double x1 = Mesh.get_edge(i,1);
        double y1 = Mesh.get_edge(i,2);
        double x2 = Mesh.get_edge(i,3);
        double y2 = Mesh.get_edge(i,4);

        // Elements on either side of edge
        int ileft  = Mesh.get_eelem(i,1);
        int iright = Mesh.get_eelem(i,2);  
        double Areal = Mesh.get_area_prim(ileft);
        double Arear = Mesh.get_area_prim(iright);

        // Scaled normal to edge
        dTensor1 nhat(2);      
        nhat.set(1, (y2-y1) );
        nhat.set(2, (x1-x2) );

        double nmag = sqrt(pow(nhat.get(1),2) + pow(nhat.get(2),2));

        // Variables to store flux integrals along edge
        dTensor2 Fr_tmp(meqn,dogParams.get_space_order());
        dTensor2 Fl_tmp(meqn,dogParams.get_space_order());

        // Loop over number of quadrature points along each edge
        for (int ell=1; ell<=dogParams.get_space_order(); ell++)
        {
            dTensor1   Ql(meqn),   Qr(meqn);
            dTensor1  ffl(meqn),  ffr(meqn);  // << -- NEW PART -- >>
            dTensor1 Auxl(maux), Auxr(maux);

            // Riemann data - q
            for (int m=1; m<=meqn; m++)
            {
                Ql.set(m, 0.0 );
                Qr.set(m, 0.0 );

                // << -- NEW PART, ffl and ffr -- >> //
                ffl.set(m, 0.0 );
                ffr.set(m, 0.0 );

                for (int k=1; k<=kmax; k++)
                {
                    Ql.set(m, Ql.get(m) + EdgeData.phi_left->get(i,ell,k) 
                            *q.get(ileft, m,k) );
                    Qr.set(m, Qr.get(m) + EdgeData.phi_right->get(i,ell,k)
                            *q.get(iright,m,k) );

                    // << -- NEW PART, ffl and ffr -- >> //
                    // Is this the correct way to use the normal vector?
                    ffl.set(m, ffl.get(m) + EdgeData.phi_left->get (i, ell, k) * ( 
                        nhat.get(1)*F.get( ileft, m, k) + nhat.get(2)*G.get( ileft, m, k) ) );

                    ffr.set(m, ffr.get(m) + EdgeData.phi_right->get(i, ell, k) * (
                        nhat.get(1)*F.get(iright, m, k) + nhat.get(2)*G.get(iright, m, k) ) );

                }
            }

            // Riemann data - aux
            for (int m=1; m<=maux; m++)
            {
                Auxl.set(m, 0.0 );
                Auxr.set(m, 0.0 );

                for (int k=1; k<=kmax; k++)
                {
                    Auxl.set(m, Auxl.get(m) + EdgeData.phi_left->get(i,ell,k)  * aux.get(ileft, m,k) );
                    Auxr.set(m, Auxr.get(m) + EdgeData.phi_right->get(i,ell,k) * aux.get(iright,m,k) );
                }
            }

            // Solve Riemann problem
            dTensor1 xedge(2);
            double s = EdgeData.xpts1d->get(ell);
            xedge.set(1, x1 + 0.5*(s+1.0)*(x2-x1) );
            xedge.set(2, y1 + 0.5*(s+1.0)*(y2-y1) );

            // Solve the Riemann problem for this edge
            dTensor1 Fl(meqn), Fr(meqn);

            // Use the time-averaged fluxes to define left and right values for
            // the Riemann solver.
            const double smax_edge = RiemannSolveLxW(
                    nhat, xedge, Ql, Qr, Auxl, Auxr, ffl, ffr, Fl, Fr);
            smax.set(i, Max(smax_edge,smax.get(i)) );

            // Construct fluxes
            for (int m=1; m<=meqn; m++)
            {
                Fr_tmp.set(m,ell, Fr.get(m) );
                Fl_tmp.set(m,ell, Fl.get(m) );
            }
        }

        // Add edge integral to line integral around the full element
        for (int m=1; m<=meqn; m++)
        for (int k=2; k<=kmax; k++)
        {
            double Fl_sum = 0.0;
            double Fr_sum = 0.0;
            for (int ell=1; ell<=dogParams.get_space_order(); ell++)
            {
                Fl_sum = Fl_sum + 0.5*EdgeData.wgts1d->get(ell)
                    *EdgeData.phi_left->get(i,ell,k) *Fl_tmp.get(m,ell);
                Fr_sum = Fr_sum + 0.5*EdgeData.wgts1d->get(ell)
                    *EdgeData.phi_right->get(i,ell,k)*Fr_tmp.get(m,ell);
            }
            EdgeFluxIntegral.set(ileft, m,k, EdgeFluxIntegral.get(ileft, m,k) + Fl_sum/Areal );
            EdgeFluxIntegral.set(iright,m,k, EdgeFluxIntegral.get(iright,m,k) - Fr_sum/Arear );
        }

     int efound1=0;
     for(int e=1;e<=3;e++)
        {
        if(i==Mesh.get_tedge(ileft,e))
        {efound1=e;}
        }


     int efound2=0;
     for(int e=1;e<=3;e++)
        {
        if(i==Mesh.get_tedge(iright,e))
        {efound2=e;}
        }

       if(efound1==0 ||  efound2==0)
       {cout<<" problem3 "<<i<<" "<<ileft<<" "<<iright<<endl;exit(1);}

        for (int m=1; m<=meqn; m++)
        {int k=1;
            double Fl_sum = 0.0;
            double Fr_sum = 0.0;
            for (int ell=1; ell<=dogParams.get_space_order(); ell++)
            {
                Fl_sum = Fl_sum + 0.5*EdgeData.wgts1d->get(ell)
                    *EdgeData.phi_left->get(i,ell,k) *Fl_tmp.get(m,ell);
                Fr_sum = Fr_sum + 0.5*EdgeData.wgts1d->get(ell)
                    *EdgeData.phi_right->get(i,ell,k)*Fr_tmp.get(m,ell);
            }
            FeI.set(ileft, m,efound1,  Fl_sum );
            FeI.set(iright,m,efound2,  -Fl_sum );
        }
    }

    //Apply Limiter to maintain positive averages.

    if( dogParams.using_positive_limiter() )
    {
        void ApplyPosMPPLimiter_Unst( const double dt, const mesh& Mesh,const edge_data_Unst& EdgeData, 
        const dTensor1& smax,
        const dTensor3& aux, const dTensor3& q, 
        dTensor3& FeI);
        
       ApplyPosMPPLimiter_Unst( dt,Mesh,EdgeData, smax, aux, q, FeI);
    }

    for (int i=1; i<=NumPhysElems; i++)
    {
        double Area = Mesh.get_area_prim(i);
        int k=1;
        for (int m=1; m<=meqn; m++)
        {
            for (int j=1;j<=3;j++)
            {double F_sum = FeI.get(i,m,j);
       if(Mesh.get_tedge(i,j)==0)
       {cout<<" "<<i<<endl;exit(1);}
            EdgeFluxIntegral.set(i, m,k, EdgeFluxIntegral.get(i,m,k) + F_sum/Area );
            }
        }
    }



    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part III: compute intra-element contributions
    // ---------------------------------------------------------
    if( dogParams.get_space_order() > 1 )
    {
        L2ProjectGradAddLegendre_Unst(1, NumPhysElems, space_order, 
            Mesh, &F, &G, &ElemFluxIntegral );
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part IV: construct Lstar
    // ---------------------------------------------------------
    if (dogParams.get_source_term()==0)  // Without Source Term
    { 
#pragma omp parallel for
        for (int i=1; i<=NumPhysElems; i++)	
        for (int m=1; m<=meqn; m++)
        for (int k=1; k<=kmax; k++)
        {
            double tmp = ElemFluxIntegral.get(i,m,k) - EdgeFluxIntegral.get(i,m,k);
            Lstar.set(i,m,k, tmp );	      
        }
    }
    else  // With Source Term
    {
#pragma omp parallel for
        for (int i=1; i<=NumPhysElems; i++)
        for (int m=1; m<=meqn; m++)
        for (int k=1; k<=kmax; k++)
        {
//          double tmp = ElemFluxIntegral.get(i,m,k) 
//              - EdgeFluxIntegral.get(i,m,k)
//              + Psi.get(i,m,k);

//          Lstar.set(i,m,k, tmp );

            printf("Source term has not been implemented for LxW solver.  Terminating program.");
            exit(1);
        }
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part V: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    LstarExtra_Unst(Mesh, &q, &aux, &Lstar);
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part VI: artificial viscosity limiter
    // ---------------------------------------------------------  
//  if (dogParams.get_space_order()>1  &&
//          dogParams.using_viscosity_limiter())
//  {  ArtificialViscosity(&aux,&q,&Lstar);  }
    // ---------------------------------------------------------

}
