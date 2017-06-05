#include <stdlib.h>
#include "dog_math.h"
#include "dogdefs.h"
#include "edge_data.h"
#include "DogParams.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"

// Generic method that user can override
//
void DogSolverCart2::DogSolveUser(double tstart, double tend)
{

    // Data that needs to be updated.  qold -> qnew at time t=tend.
    dTensorBC4& qnew   = fetch_state().fetch_q();
    dTensorBC4& aux    = fetch_state().fetch_aux();
    dTensorBC4& qold   = fetch_state_old().fetch_q();
    dTensorBC4& auxold = fetch_state_old().fetch_aux();
    dTensorBC3& smax   = fetch_smax();

    // access constant data
    const int time_order = dogParams.get_time_order();
    const int nv         = dogParams.get_nv();
    const double* cflv   = dogParams.get_cflv();

    // TODO: can use dogParamsCart2 to get these variables:
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();
    const int maux = aux.getsize(3);

    // ------------------------------------------------------------
    // Function definitions
    void CopyQ(const dTensorBC4& qin,dTensorBC4& qout);
    // ------------------------------------------------------------

    // define local variables
    int n_step = 0;
    double t = tstart;
    double dt = get_dt();
    const double CFL_max = cflv[1];
    const double CFL_target = cflv[2];

    // Example of using extra information and intermediate stages:
    dTensorBC4   qstar(mx,my,meqn,kmax,mbc);
    dTensorBC4      q1(mx,my,meqn,kmax,mbc);
    dTensorBC4      q2(mx,my,meqn,kmax,mbc);
    dTensorBC4 auxstar(mx,my,maux,kmax,mbc);
    dTensorBC4   Lstar(mx,my,meqn,kmax,mbc);

    // Set initialize qstar and auxstar values
    CopyQ(qold,   qstar);
    CopyQ(aux,  auxstar);

    // User-defined time stepping
    while (t<tend)
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            eprintf(" Error in DogSolveUser.cpp: "
                    " Exceeded allowed # of time steps \n"
                    "    n_step = %d\n"
                    "        nv = %d\n\n",
                    n_step,nv);
        }

        // copy qnew into qold
        CopyQ(qnew,qold);
        CopyQ(aux,auxold);

        // keep trying until we get a dt that does not violate CFL condition
        while (m_accept==0)
        {
            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            set_time_hack(told);
            dogParams.set_time(told);
            set_dt(dt);

            // Set initial maximum wave speed to zero
            for (int j=1-mbc; j<=(my+mbc); j++)
            {
                for (int i=1-mbc; i<=(mx+mbc); i++)
                {
                    smax.set(i,j,1, 0.0e0 );
                    smax.set(i,j,2, 0.0e0 );
                }
            }

            // ----------------------------------------------------------------
            //
            //    THIS IS WHERE THE USER-DEFINED TIME-STEPPING SCHEME
            //    SHOULD BE ADDED. IN THE DEFAULT FILE: DogSolveUser.cpp,
            //    THE PROGRAM WILL NOW RETURN AN ERROR MESSAGE.
            // 
            // ----------------------------------------------------------------
            eprintf("\n"
                    " No user-defined time-stepping scheme has been defined yet. \n"
                    " Copy $DOGPACK/lib/2d/cart/DogSolveUser.cpp into the current \n"
                    " directory and modify as needed.\n\n");
            // ----------------------------------------------------------------

            // compute cfl number
            double cfl = GetCFL(dt);

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

                // do any extra work
                ::AfterFullTimeStep(fetch_solver());
            }
            else 
                //reject
            {   
                t = told;
                if (dogParams.get_verbosity()>0)
                {
                    printf("DogSolve2D rejecting step..."
                            "CFL number too large\n");
                }

                // copy qold into qnew
                CopyQ(qold,qnew);
                CopyQ(auxold,aux);
            }      
        }

        // compute conservation and print to file
        ConSoln(aux,qnew,t);
    }

    // set initial time step for next call to DogSolveUser
    set_dt(dt);
}
