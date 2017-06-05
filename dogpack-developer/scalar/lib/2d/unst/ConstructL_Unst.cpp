#include "dogdefs.h"
#include "dog_math.h"
#include "mesh.h"
#include "edge_data_Unst.h"
#include "stdlib.h"
#include "DogParams.h"          // for get_space_order and get_source_term
#include "ConstructL_Unst.h"    // Function definitions

void ConstructL_Unst(
    const double t,
    const dTensor2* vel_vec,
    const mesh& Mesh,
    const edge_data_Unst& EdgeData,
    dTensor3& aux, // SetBndValues_Unst modifies ghost cells
    dTensor3& q,   // SetBndValues_Unst modifies ghost cells
    dTensor3& Lstar, 
    dTensor1& smax)
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
    SetBndValues_Unst(Mesh,&q,&aux);  
    
    // Positivity limiter
    void ApplyPosLimiter_Unst(const mesh& Mesh, const dTensor3& aux, dTensor3& q);
    if( dogParams.using_moment_limiter() )
    { ApplyPosLimiter_Unst(Mesh, aux, q); }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part I: compute flux integral on element edges
    // ---------------------------------------------------------

    // Loop over all interior edges and solve Riemann problems
    // dTensor1 nvec(2);

    // Loop over all interior edges
    EdgeFluxIntegral.setall(0.);
    ElemFluxIntegral.setall(0.);

    // Loop over all interior edges
#pragma omp parallel for
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

        // Variables to store flux integrals along edge
        dTensor2 Fr_tmp(meqn,dogParams.get_space_order());
        dTensor2 Fl_tmp(meqn,dogParams.get_space_order());

        // Loop over number of quadrature points along each edge
        for (int ell=1; ell<=dogParams.get_space_order(); ell++)
        {
            dTensor1 Ql(meqn),Qr(meqn);
            dTensor1 Auxl(maux),Auxr(maux);	  

            // Riemann data - q
            for (int m=1; m<=meqn; m++)
            {
                Ql.set(m, 0.0 );
                Qr.set(m, 0.0 );

                for (int k=1; k<=kmax; k++)
                {
                    Ql.set(m, Ql.get(m) + EdgeData.phi_left->get(i,ell,k) 
                            *q.get(ileft, m,k) );
                    Qr.set(m, Qr.get(m) + EdgeData.phi_right->get(i,ell,k)
                            *q.get(iright,m,k) );
                }

            }


            // Riemann data - aux
            for (int m=1; m<=maux; m++)
            {
                Auxl.set(m, 0.0 );
                Auxr.set(m, 0.0 );

                for (int k=1; k<=kmax; k++)
                {
                    Auxl.set(m, Auxl.get(m) + EdgeData.phi_left->get(i,ell,k)
                            *aux.get(ileft, m,k) );
                    Auxr.set(m, Auxr.get(m) + EdgeData.phi_right->get(i,ell,k)
                            *aux.get(iright,m,k) );
                }
            }

            // Solve Riemann problem
            dTensor1 xedge(2);
            double s = EdgeData.xpts1d->get(ell);
            xedge.set(1, x1 + 0.5*(s+1.0)*(x2-x1) );
            xedge.set(2, y1 + 0.5*(s+1.0)*(y2-y1) );
            dTensor1 Fl(meqn),Fr(meqn);
            const double smax_edge = RiemannSolve(vel_vec, nhat, xedge, Ql, Qr, Auxl, Auxr, Fl, Fr);
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
        for (int k=1; k<=kmax; k++)
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
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part II: compute intra-element contributions
    // ---------------------------------------------------------
    L2ProjectGrad_Unst(vel_vec, 1,NumPhysElems,
            space_order,space_order,space_order,space_order,
            Mesh,&q,&aux,&ElemFluxIntegral,&FluxFunc);
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part III: compute source term
    // --------------------------------------------------------- 
    if ( dogParams.get_source_term()>0 )
    {        
        // Set source term on computational grid
        // Set values and apply L2-projection
        L2Project_Unst(t, vel_vec, 1,NumPhysElems,
                space_order,space_order,space_order,space_order,
                Mesh,&q,&aux,&Psi,&SourceTermFunc);
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
            double tmp = ElemFluxIntegral.get(i,m,k) 
                - EdgeFluxIntegral.get(i,m,k)
                + Psi.get(i,m,k);

            Lstar.set(i,m,k, tmp );
        }
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part V: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    LstarExtra_Unst(Mesh,&q,&aux,&Lstar);
    // ---------------------------------------------------------

}
