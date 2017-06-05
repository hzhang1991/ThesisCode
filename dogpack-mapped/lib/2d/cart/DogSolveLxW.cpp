#include <stdlib.h>
#include "dog_math.h"
#include "dogdefs.h"
#include "edge_data.h"
#include "DogParams.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"

// Lax-Wendroff time stepping
//
// TODO - I guess the idea is to try and incorporate this into DogSolve?  I
// don't see how this can be done with this method for all dimensions ... (-DS)
void DogSolverCart2::DogSolveLxW(double tstart, double tend)
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
    double t   = tstart;
    double dt  = get_dt();
    const double CFL_max    = cflv[1];
    const double CFL_target = cflv[2];

    // Example of using extra information and intermediate stages:
    dTensorBC4   qstar(mx, my, meqn, kmax, mbc);
    dTensorBC4 auxstar(mx, my, maux, kmax, mbc);

    // This memory should have already been allocated:
    // (otherwise, can allocate it here, as in above):
    dTensorBC4& L = fetch_L();

    // Set initialize qstar and auxstar values
    CopyQ(qold,   qstar);
    CopyQ(aux,  auxstar);

    // User-defined time stepping
    while (t<tend)
    {
        // initialize time step
        int m_accept = 0;      
        n_step       = n_step + 1;

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
        CopyQ(qnew,   qold);
        CopyQ(aux,  auxold);

        // keep trying until we get a dt that does not violate CFL condition
        while (m_accept==0)
        {
            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            fetch_state().set_time(told);
            dogParams.set_time(told);      // TODO
            set_dt(dt);

            // Set initial maximum wave speed to zero
            //
            // TODO : CAN CALL reset smax or something ...
            for (int j=1-mbc; j<=(my+mbc); j++)
            {
                for (int i=1-mbc; i<=(mx+mbc); i++)
                {
                    smax.set(i, j, 1, 0.0e0 );
                    smax.set(i, j, 2, 0.0e0 );
                }
            }

            BeforeStep(dt,aux,qnew,*this);
            void SetBndValues(dTensorBC4&,dTensorBC4&);
            SetBndValues(qnew,aux);

            // Construct RHS for LxW formulation:
            void LaxWendroff(double dt, 
                double alpha1,  double beta1,      
                dTensorBC4& aux, dTensorBC4& q,    
                dTensorBC4& Lstar, dTensorBC3& smax);
            LaxWendroff(dt, 1.0, 0.5, aux, qnew, L, smax);

            // Take a single "Euler" time step:
            void StepLxW( double dt, const dTensorBC4& qold, 
                const dTensorBC4& L, dTensorBC4& qnew );
            StepLxW(dt, qnew, L, qnew);

            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,qnew,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt, aux, qnew, *this);

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
                CopyQ(qold,  qnew);
                CopyQ(auxold, aux);
            }      
        }

        // compute conservation and print to file
        ConSoln(aux,qnew,t);
    }

    // set initial time step for next call to DogSolveUser
    set_dt(dt);
}

void StepLxW( double dt, const dTensorBC4& qold,
    const dTensorBC4& L, dTensorBC4& qnew )
{

    const int numel = qnew.numel();
    assert_eq(qold.numel(),numel);
    assert_eq(L.numel(),numel);

#pragma omp parallel for
    for(int v=0;v<numel;v++)
    {
        double tmp = qold.vget(v) + dt*L.vget(v);
        qnew.vset(v, tmp );
    }

}
