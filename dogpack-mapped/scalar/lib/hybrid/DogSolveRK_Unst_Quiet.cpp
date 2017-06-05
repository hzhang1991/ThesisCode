#include "dogdefs.h"
#include "dog_math.h"
#include "edge_data_Unst.h"
#include "RKinfo.h"
#include "mesh.h"
#include "DogParams.h"
#include "DogStateUnst2.h"          // For access to 
#include "DogSolveRK_Unst_Quiet.h"

// This file should be identical to DogSolveRK_Unst, with the exception that all
// output printing statements are silenced.
//
// Advance the solution qold to qnew over time interval tstart to tend.
//
// All local information is allocated within this function.  The only part
// that gets shared are time values passed through dogStateUnst2.  This class
// should be modified to accept the state variable, q and aux in place of only
// containing time information as is currently the case.  (-DS)
double DogSolveRK_Unst_Quiet(
    const dTensor2* vel_vec,
    const mesh& Mesh, const edge_data_Unst& EdgeData,
    dTensor3& aux, dTensor3& qold, dTensor3& qnew, 
    const double tstart, const double tend, DogStateUnst2& dogStateUnst2)
{

    const int mx   = qnew.getsize(1);
    const int meqn = qnew.getsize(2);
    const int kmax = qnew.getsize(3);
    const int maux = aux.getsize(2);
    const double* cflv = dogParams.get_cflv();
    const int nv   = dogParams.get_nv();

    RKinfo rk;
    SetRKinfo(dogParams.get_time_order(),rk);

    // define local variables
    int n_step = 0;
    double t  = tstart;
    //double dt = dogStateUnst2.get_initial_dt();

    double dt = dogStateUnst2.get_initial_dt();

    const double CFL_max    = cflv[1];
    const double CFL_target = cflv[2];
    double cfl   = 0.0;
    double dtmin = dt;
    double dtmax = dt;

    const int NumElems = Mesh.get_NumElems(); // Number of total elements in mesh
    const int NumNodes = Mesh.get_NumNodes(); // Number of nodes in mesh
    const int NumEdges = Mesh.get_NumEdges(); // Number of edges in mesh 

    dTensor3   qstar(NumElems,meqn,kmax);
    dTensor3      q1(NumElems,meqn,kmax);
    dTensor3      q2(NumElems,meqn,kmax);
    dTensor3 auxstar(NumElems,maux,kmax);
    dTensor3   Lstar(NumElems,meqn,kmax);
    dTensor3    Lold(NumElems,meqn,kmax);
    dTensor3  auxold(NumElems,maux,kmax);
    dTensor1    smax(NumEdges);

    void L2Project_Unst( 
        const dTensor2* vel_vec, 
        const int istart, const int iend,
        const int QuadOrder, const int BasisOrder_qin, const int BasisOrder_auxin, const
        int BasisOrder_fout, const mesh& Mesh, const dTensor3* qin, const dTensor3*
        auxin, dTensor3* fout, void (*Func)(const dTensor2* vel_vec, const
        dTensor2&,const dTensor2&, const dTensor2&,dTensor2&));

    // JUNK here:
    void AuxFuncWrapper(
        const dTensor2* vel_vec,
        const dTensor2& xpts,
        const dTensor2& NOT_USED_1,
        const dTensor2& NOT_USED_2,
        dTensor2& auxvals);
    const int space_order = dogParams.get_space_order();
    if( maux > 0 )
    { 
        printf("WARNING: maux = %d should be zero for Vlasov routines.", maux);
        printf("    Modify parameters.ini to remove this warning\n" );
        L2Project_Unst(vel_vec,1,NumElems,
                space_order,space_order,space_order,space_order,		       
                Mesh,&qnew,&aux,&aux,&AuxFuncWrapper);  
    }

    // Set initialize qstar and auxstar values
    qstar.copyfrom(qold);
    auxstar.copyfrom(aux);

    // Runge-Kutta time stepping
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

        // copy qnew into qold
        qold.copyfrom(qnew);
        auxold.copyfrom(aux);

        // keep trying until we get a dt that does not violate CFL condition
        while (m_accept==0)
        {

            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // TODO - this needs to be performed at the 'local' level
            dogStateUnst2.set_time ( told );
            dogStateUnst2.set_dt   ( dt   );

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // Take a full time step of size dt
            switch ( dogParams.get_time_order() )
            {

                case 1:  // First order in time (1-stage)


                    // -----------------------------------------------
                    // Stage #1 (the only one in this case)  
                    rk.mstage = 1;
                    BeforeStep_Unst(dt,Mesh,aux,qnew);
                    ConstructL_Unst(told, vel_vec,Mesh,EdgeData,aux,qnew,Lstar,smax);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux, qnew, Lstar, qnew);
                    AfterStep_Unst(dt,Mesh,aux,qnew);
                    // -----------------------------------------------
                    break;
  
                case 2:  // Second order in time (2-stages)

                    // -----------------------------------------------
                    // Stage #1  	        
                    rk.mstage = 1;
                    dogStateUnst2.set_time(told);
                    BeforeStep_Unst(dt,Mesh,aux,qnew);
                    ConstructL_Unst(told,vel_vec,Mesh,EdgeData,aux,qnew,Lstar,smax);
                    UpdateSoln_Unst(
                        rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                        rk.beta->get(rk.mstage), dt, Mesh, aux, qnew, Lstar, qstar);
                    AfterStep_Unst(dt, Mesh, auxstar, qstar);

                    // ------------------------------------------------
                    // Stage #2
                    rk.mstage = 2;
                    dogStateUnst2.set_time(told+dt);
                    BeforeStep_Unst(dt, Mesh, auxstar, qstar);
                    ConstructL_Unst(told+1.0*dt, vel_vec, Mesh, EdgeData, aux, qstar, Lstar, smax);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage), rk.alpha2->get(rk.mstage), 
                            rk.beta->get(rk.mstage), dt, Mesh, auxstar, qstar, Lstar, qnew);
                    AfterStep_Unst(dt, Mesh, aux, qnew);
                    // ------------------------------------------------
                    break;

                case 3:  // Third order in time (3-stages)

//     qnew = alpha1 * qstar + alpha2 * qnew + beta * dt * L( qstar )

// alpha1 = 1.0
// alpha2 = 0.0
// beta   = 1.0

                    // ------------------------------------------------
                    // Stage #1
                    rk.mstage = 1;
                    dogStateUnst2.set_time(told);
                    BeforeStep_Unst(dt,Mesh,aux,qnew);	      
                    ConstructL_Unst(told, vel_vec,Mesh,EdgeData,aux,qnew,Lstar,smax);
                    Lold.copyfrom(Lstar);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux,qnew,Lstar,qstar);
                    AfterStep_Unst(dt,Mesh,aux,qstar);
                    // -------------------------------------------------

// alpha1 = 0.75
// alpha2 = 0.25
// beta   = 0.25

                    // Stage #2
                    rk.mstage = 2;
                    dogStateUnst2.set_time(told+0.5*dt);
                    BeforeStep_Unst(dt,Mesh,aux,qstar);
                    ConstructL_Unst(told+dt,  vel_vec,Mesh,EdgeData,aux,qstar,Lstar,smax);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux,qnew,Lstar,qstar);   
                    AfterStep_Unst(dt,Mesh,aux,qstar);
                    // --------------------------------------------------

// alpha1 = 2/3
// alpha2 = 1/3
// beta   = 2/3

                    // Stage #3
                    rk.mstage = 3;
                    dogStateUnst2.set_time(told+dt);
                    BeforeStep_Unst(dt,Mesh,auxstar,qstar);
                    ConstructL_Unst(told+0.5*dt,vel_vec,Mesh,EdgeData,auxstar,qstar,Lstar,smax);
                    UpdateSoln_Unst(rk.alpha1->get(rk.mstage),rk.alpha2->get(rk.mstage),
                            rk.beta->get(rk.mstage),dt,Mesh,aux,qstar,Lstar,qnew);   
                    AfterStep_Unst(dt,Mesh,aux,qnew);
                    // --------------------------------------------------   
                    break;


                default: unsupported_value_error(dogParams.get_time_order());

            }

            // compute cfl number
            cfl = GetCFL_Unst(dt,Mesh,aux,smax);

            // output time step information
            //if (dogParams.get_verbosity()>0) 
//          {
//              printf("*** DogSolve2D ... Step %5d"
//                      "   CFL =%6.3f"
//                      "   dt =%11.3e"
//                      "   t =%11.3e\n",
//                      n_step,cfl,dt,t);
//          }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dogParams.get_max_dt(), dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dogParams.get_max_dt();
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)
                // accept
            { 
                m_accept = 1; 
                dogStateUnst2.set_time(t);

                // do any extra work
//              AfterFullTimeStep_Unst(dogStateUnst2.get_dt(),Mesh,
//                      auxold,qold,Lold,aux,qnew);
            }
            else 
                //reject
            {   
                t = told;
                dogStateUnst2.set_time(told);
//              if (dogParams.get_verbosity()>0)
//              {
//                  printf("DogSolve2D rejecting step..."
//                          "CFL number too large\n");
//              }

                // copy qold into qnew
                qnew.copyfrom(qold);
                aux.copyfrom(auxold);

                // after reject function	      
//              AfterReject_Unst(Mesh,dt,aux,qnew);
            }      
        }

    }

    // set initial time step for next call to DogSolveRK
    dogStateUnst2.set_initial_dt(dt);

    void DeleteRKInfo(RKinfo& rk);
    DeleteRKInfo(rk);

    return cfl;

}
