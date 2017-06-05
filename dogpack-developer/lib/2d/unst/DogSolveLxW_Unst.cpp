#include <stdlib.h>
#include "dog_math.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "DogStateUnst2.h"
#include "mesh.h"
#include "edge_data_Unst.h"
#include "DogSolveLxW_Unst.h"
#include <iostream>
using namespace std;

void SetBndValues_Unst(const mesh& Mesh, dTensor3* q, dTensor3* aux);

// Lax-Wendroff time stepping (for the unstructured code).
//
// See also: DogSolveRK_Unst.
void DogSolveLxW_Unst(const mesh& Mesh, const edge_data_Unst& EdgeData,
    dTensor3& aux, dTensor3& qold, dTensor3& qnew, 
    const double tstart, const double tend, 
    const char* outputdir)
{

    // Quick error check
    if( dogParams.get_time_order() < 1 || dogParams.get_time_order() > 3 )
    {
        printf("Error: get_time_order = %d\n", dogParams.get_time_order() );
        printf("Set to 1, 2 or 3 in order to run Lax-Wendroff method\n");
        exit(1);
    }

    // The usual parameters
    const int meqn = qnew.getsize(2);
    const int kmax = qnew.getsize(3);
    const int maux = aux.getsize(2);

    // Access constant data
    const int time_order   = dogParams.get_time_order();
    const double* cflv     = dogParams.get_cflv();
    const int     nv       = dogParams.get_nv();

    // define local variables
    int n_step          = 0;      // counter for the number of steps
    double t            = tstart;
    double dt           = dogStateUnst2.get_initial_dt();
    double CFL_max      = cflv[1];
    double CFL_target   = cflv[2];
    double cfl          = 0.0;
    double dtmin        = dt;
    double dtmax        = dt;
    const int NumElems  = Mesh.get_NumElems(); // Number of total elements in mesh
    const int NumNodes  = Mesh.get_NumNodes(); // Number of nodes in mesh
    const int NumEdges  = Mesh.get_NumEdges(); // Number of edges in mesh 

    // Storage
    dTensor3  auxold(NumElems,maux,kmax);
    dTensor3    Lold(NumElems,meqn,kmax);

    dTensor3   qstar(NumElems,meqn,kmax);
    dTensor3 auxstar(NumElems,maux,kmax);

    dTensor3   Lstar(NumElems, meqn, kmax);     // Right hand side of PDE
    dTensor1    smax(NumEdges            );     // Maximum wave speed

    qstar.copyfrom( qold );
    auxstar.copyfrom( aux );

    // --------------------------- //
    // MAIN TIME STEPPING LOOP
    // --------------------------- //
    while (t<tend)
    {

        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            eprintf(" Error in DogSolveRK_Unst.cpp: "
                    " Exceeded allowed # of time steps \n"
                    "    n_step = %d\n"
                    "        nv = %d\n\n",
                    n_step,nv);
        }

        // copy qnew into qold (for rejects)
        qold  .copyfrom( qnew );
        auxold.copyfrom( aux  );

        // keep trying until we get a dt that does not violate CFL condition
        while (m_accept==0)
        {   
            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;
            dogStateUnst2.set_time(told);
            dogStateUnst2.set_dt(dt);

            // Set initial maximum wave speed to zero
            for (int i=1; i<=NumEdges; i++)
            {  smax.set(i, 0.0e0 );  }


            // ------------------------------------------------------------- //
            // TAKE A SINGLE TIME STEP HERE!
            // ------------------------------------------------------------- //
            if( dogParams.using_positive_limiter() )
            {
              void ApplyPosLimiter_Unst(const mesh& Mesh, const dTensor3& aux, dTensor3& q);
              ApplyPosLimiter_Unst(Mesh,aux,qnew);
            }


            BeforeStep_Unst(dt, Mesh, aux, qnew);
            LaxWendroff_Unst(dt, Mesh, EdgeData, aux, qnew, Lstar,smax);
            Lold.copyfrom( Lstar );
            StepLxW( dt, qold, Lstar, qnew );
            AfterStep_Unst(dt, Mesh, aux, qnew);
            SetBndValues_Unst(Mesh, &qnew, &aux);
            // ------------------------------------------------------------- //

            // compute cfl number
            cfl = GetCFL_Unst(dt, Mesh, aux, smax);
            //TODO: How to indicate we should be using this limiter?
            if( dogParams.get_use_limiter() )
            {
              void ApplyShockLimiter_Unst(const mesh& Mesh,const edge_data_Unst&, const dTensor3& aux, dTensor3& q);
              ApplyShockLimiter_Unst(Mesh,EdgeData,aux,qnew);
            }

            // output time step information
            if (dogParams.get_verbosity()>0) 
            {
                printf("DogSolve2D ... Step %5d"
                        "   CFL =%6.3f"
                        "   dt =%11.3e"
                        "   t =%11.3e\n",
                        n_step,cfl,dt,t);
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dogParams.get_max_dt(),dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dogParams.get_max_dt();
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)
            { 
                // -- accept-- //
                m_accept = 1; 
                dogStateUnst2.set_time(t);

                // do any extra work
                AfterFullTimeStep_Unst( dogStateUnst2.get_dt(), Mesh,
                        auxold, qold, Lold, aux, qnew);

            }
            else 
            {   
                // -- reject -- //

                t = told;
                dogStateUnst2.set_time(told);
                if (dogParams.get_verbosity()>0)
                {
                    printf("DogSolve2D rejecting step..."
                            "CFL number too large\n");
                }

                // copy qold into qnew
                qnew.copyfrom( qold  );
                aux.copyfrom( auxold );

                // after reject function          
                AfterReject_Unst(Mesh, dt, aux, qnew);
            }      
        }

        // compute conservation and print to file
        ConSoln_Unst(Mesh, aux, qnew, t, outputdir);

    }

    // set initial time step for next call to DogSolveRK
    dogStateUnst2.set_initial_dt(dt);




}

void StepLxW( double dt, const dTensor3& qold,
    const dTensor3& L, dTensor3& qnew )
{

    const int numel = qnew.numel();
    assert_eq(qold.numel(), numel);
    assert_eq(L.numel(),numel);

#pragma omp parallel for
    for(int v=0; v<numel; v++)
    {
        double tmp = qold.vget(v) + dt*L.vget(v);
        qnew.vset(v, tmp );
    }

}
