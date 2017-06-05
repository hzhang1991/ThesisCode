#include "ConstructL_Unst.h"
#include <iostream>
using namespace std;

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t = -( f(q,x,y,t)_x + g(q,x,y,t)_y ) + Psi(q,x,y,t)
//
void ConstructL_LLF_Unst(const mesh& Mesh,
        const edge_data_Unst& EdgeData,
        const dTensor3& aux,                  // SetBndValues modifies ghost cells
        const dTensor3& q,                    // SetBndValues modifies ghost cells 
        const dTensor1& smax, dTensor3& FeI)
{
    const int NumElems = Mesh.get_NumElems();
    const int NumPhysElems = Mesh.get_NumPhysElems();
    const int NumEdges = Mesh.get_NumEdges();
    const int meqn = q.getsize(2);
    const int kmax = q.getsize(3);
    const int maux = aux.getsize(2);
    const int space_order = dogParams.get_space_order();
    dTensor3 EdgeFluxIntegral(NumElems,meqn,kmax);
    dTensor3 ElemFluxIntegral(NumElems,meqn,kmax);
    dTensor3              Psi(NumElems,meqn,kmax);


    // ---------------------------------------------------------
    // Part I: compute flux integral on element edges
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

        // Variables to store flux integrals along edge
        dTensor2 Fr_tmp(meqn,1);
        dTensor2 Fl_tmp(meqn,1);

        // Loop over number of quadrature points along each edge
        for (int ell=1; ell<=1; ell++)
        {
            dTensor1 Ql(meqn),Qr(meqn);
            dTensor1 Auxl(maux),Auxr(maux);	  

            // Riemann data - q
            for (int m=1; m<=meqn; m++)
            {
                Ql.set(m, 0.0 );
                Qr.set(m, 0.0 );
                //TODO: Fix this hack, figure out how to get edge lengths most efficiently...
               /* {int k=1;
                    Ql.set(m, Ql.get(m) + EdgeData.phi_left->get(i,ell,k) 
                            *q.get(ileft, m,k) );
                    Qr.set(m, Qr.get(m) + EdgeData.phi_right->get(i,ell,k)
                            *q.get(iright,m,k) );
                }*/
            Ql.set(m,q.get(ileft,m,1));
            Qr.set(m,q.get(iright,m,1));
            }
 

            // Riemann data - aux
            for (int m=1; m<=maux; m++)
            {
                Auxl.set(m, 0.0 );
                Auxr.set(m, 0.0 );
                /*
                {int k=1;
                    Auxl.set(m, Auxl.get(m) + EdgeData.phi_left->get(i,ell,k)
                            *aux.get(ileft, m,k) );
                    Auxr.set(m, Auxr.get(m) + EdgeData.phi_right->get(i,ell,k)
                            *aux.get(iright,m,k) );
                }*/
            Auxl.set(m,aux.get(ileft,m,1));
            Auxr.set(m,aux.get(iright,m,1));
            }

            // Solve Riemann problem
            dTensor1 xedge(2);
            //double s = EdgeData.xpts1d->get(ell);
            double s=0.0;
            xedge.set(1, x1 + 0.5*(s+1.0)*(x2-x1) );
            xedge.set(2, y1 + 0.5*(s+1.0)*(y2-y1) );
            dTensor1 Fl(meqn),Fr(meqn);
            const double smax_edge = RiemannSolve(nhat, xedge, Ql, Qr, Auxl, Auxr, Fl, Fr);


            // Construct fluxes
            for (int m=1; m<=meqn; m++)
            {
                Fr_tmp.set(m,1, Fr.get(m) );
                Fl_tmp.set(m,1, Fl.get(m) );
            }
        }

        // Add edge integral to line integral around the full element
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
       {cout<<" "<<i<<" "<<ileft<<" "<<iright<<endl;exit(1);}

        for (int m=1; m<=meqn; m++)
        {int k=1;
            double Fl_sum = 0.0;
            double Fr_sum = 0.0;
            //TODO: More hacking here...fix this 
            for (int ell=1; ell<=1; ell++)
            {
                Fl_sum = Fl_sum + Fl_tmp.get(m,1);
                /*Fr_sum = Fr_sum + 0.5*EdgeData.wgts1d->get(ell)
                    *EdgeData.phi_right->get(i,ell,k)*Fr_tmp.get(m,ell);*/
            }
            FeI.set(ileft, m,efound1,  Fl_sum );
            FeI.set(iright,m,efound2,  -Fl_sum );
        }



    }
    // ---------------------------------------------------------


}
