#include "dogdefs.h"
#include "dog_math.h"
#include "edge_data.h"
#include "DogParams.h"
#include "DogStateCart2.h"
#include "DogSolverCart2.h"
#include "SDCinfo.h"
#include <string>

// SDC time stepping.
//
void DogSolverCart2::DogSolveSDC(double tstart, double tend)
{
    // access constant data
    const int nv         = dogParams.get_nv();
    const double* cflv   = dogParams.get_cflv();
    const int time_order = dogParams.get_time_order();

    // access variable state data
    dTensorBC4& qnew = fetch_state().fetch_q();
    dTensorBC4& aux  = fetch_state().fetch_aux();
    dTensorBC3& smax = fetch_smax();

    // Define local variables
    int n_step = 0;
    double t = tstart;
    assert_eq(tstart,get_state().get_time());
    double dt = get_dt();
    const double CFL_max    = cflv[1];
    const double CFL_target = cflv[2];
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();
    const int maux = aux.getsize(3);

    // --------------------------------------------------------------
    // create helper arrays
    //
    // create current time step vector
    dTensor1    dtvec(time_order-1);
    dTensor1    tvec(time_order);
    //
    // allocate memory
    //
    dTensorBC4& L0 = fetch_L(); //  L0(mx,my,meqn,kmax,mbc);
    //
    const int arraySize = 5; // why not time_order?
    dTensorBC4** L    = new dTensorBC4* [arraySize];
    dTensorBC4** Lnew = new dTensorBC4* [arraySize];
    dTensorBC4** q    = new dTensorBC4* [arraySize];
    for(int i=0;i<arraySize;i++) L[i] = Lnew[i] = q[i] = 0;
    L[0] = &L0;
    for(int i=1;i<time_order;i++)
    {
        L[i]    = new dTensorBC4(mx,my,meqn,kmax,mbc);
        q[i]    = new dTensorBC4(mx,my,meqn,kmax,mbc);
    }
    for(int i=1;i<time_order-1;i++)
    { Lnew[i]   = new dTensorBC4(mx,my,meqn,kmax,mbc); }

    // Integrated residual
    dTensorBC5 IL(mx,my,meqn,kmax,time_order-1,mbc);
    // --------------------------------------------------------------

    // -----------------------
    // Construct initial L
    // -----------------------  
    if (time_order>2)
    { 
        // Construct initial right-hand side for SDC
        sdc.num_iter  = 1;
        sdc.num_stage = 0;     
        BeforeStep(dt,aux,qnew,*this);
        void SetBndValues(dTensorBC4&,dTensorBC4&);
        SetBndValues(qnew,aux);
        ::ConstructL(aux,qnew,L0,smax);
    }

    // -----------------------
    // Main time-stepping loop
    // -----------------------
    while (t<tend)
    {
        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            eprintf(" Error in DogSolveSDC.cpp: "
                    " Exceeded allowed # of time steps \n"
                    "    n_step = %d\n"
                    "        nv = %d\n\n",
                    n_step,nv);
        }

        // Copy qnew into qold (in order to save data)
        saveState();

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time	  
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            assert_eq(get_state().get_time(), told);
            set_dt(dt);

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            advanceTimeStepSDC(dt,L,Lnew,q,IL,dtvec,tvec);

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
                fetch_state().set_time(t);

                // Copy old L into new L0
                if(time_order>=3)
                    L0.copyfrom(*L[time_order-1]);

                // Call AfterFullTimeStep
                void AfterFullTimeStep(DogSolverCart2& solver);
                AfterFullTimeStep(*this);
            }
            else 
                //reject
            {   
                if (dogParams.get_verbosity()>0)
                {
                    printf("DogSolve2D rejecting step..."
                            "CFL number too large\n");
                }

                // Copy old data back into qnew
                revertToSavedState();
                t = get_state().get_time();
                assert_eq(t, told);
            }
        }

        // compute conservation and print to file
        ConSoln(aux,qnew,t);
    }

    // set initial time step for next call to DogSolve
    set_dt(dt);

    // free memory
    //
    for(int i=1;i<time_order;i++)
        delete L[i];
    for(int i=1;i<time_order-1;i++)
        delete Lnew[i];
    for(int i=1;i<time_order;i++)
        delete q[i];
    delete [] L;
    delete [] Lnew;
    delete [] q;
}

// Take a full time step of size dt
//
// I think this function needs to be changed to maintain
// get_state().get_time() at the current time
// at each time stage (in order to agree with the
// rest of the state data in q and aux).
//
void DogSolverCart2::advanceTimeStepSDC(
        double dt,
        dTensorBC4** L,
        dTensorBC4** Lnew,
        dTensorBC4** q,
        dTensorBC5& IL,
        dTensor1& dtvec,
        dTensor1& tvec)
{
    // declare methods
    void SetBndValues(dTensorBC4&,dTensorBC4&);
    void StepSDC(double,const dTensorBC4&, const dTensorBC4&, 
            const dTensorBC4&, dTensorBC4&);
    void ResInt(double dt, dTensorBC4 const*const*const L, dTensorBC5& ILout);
    void CopyQ(const dTensorBC4& qin,dTensorBC4& qout);

    // create convenience aliases
    //
    dTensorBC4& qnew = fetch_state().fetch_q();
    dTensorBC4& aux  = fetch_state().fetch_aux();
    dTensorBC3& smax = fetch_smax();
    const int mx   = qnew.getsize(1);
    const int my   = qnew.getsize(2);
    const int meqn = qnew.getsize(3);
    const int kmax = qnew.getsize(4);
    const int mbc  = qnew.getmbc();
    const int maux = aux.getsize(3);

    dTensorBC4& q1 = *q[1];
    dTensorBC4& q2 = *q[2];
    dTensorBC4& q3 = *q[3];
    dTensorBC4& q4 = *q[4];
    //
    dTensorBC4& L0 = *L[0];
    dTensorBC4& L1 = *L[1];
    dTensorBC4& L2 = *L[2];
    dTensorBC4& L3 = *L[3];
    dTensorBC4& L4 = *L[4];
    //
    dTensorBC4& L1new = *Lnew[1];
    dTensorBC4& L2new = *Lnew[2];
    dTensorBC4& L3new = *Lnew[3];

    const int time_order = dogParams.get_time_order();
    void SetSDCtimePoints(int morder, double t, double dt, 
            dTensor1& dtvec, dTensor1& tvec);
    SetSDCtimePoints(time_order,get_state().get_time(),dt,dtvec,tvec);

    for (int i=0; i<time_order; i++)
    {  sdc.timevec[i] = tvec.get(i+1);  }
    sdc.max_iter = time_order;

    // Choose the SDC order of accuracy
    switch ( time_order )
    {
        case 2:  // 2nd order in time

            // --------------------------------------------------------
            // Construct initial right-hand side for SDC
            sdc.num_iter = 1;
            sdc.num_stage = 0;                        
            BeforeStep(dt,aux,qnew,*this);
            SetBndValues(qnew,aux);
            ::ConstructL(aux,qnew,L0,smax);

            // Take an Euler time step 
            sdc.num_stage = 1;
            StepSDC(dtvec.get(1),aux,qnew,L0,q1);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q1,*this);

            // Construct new right-hand side            
            BeforeStep(dt,aux,q1,*this);
            SetBndValues(q1,aux);
            ::ConstructL(aux,q1,L1,smax);

            // Integrate residual over time step
            ResInt(dt,L,IL);

            // Correct Q
            sdc.num_iter = 2;
#pragma omp parallel for
            for (int i=(1-mbc); i<=(mx+mbc); i++)
                for (int j=(1-mbc); j<=(my+mbc); j++)
                    for (int m=1; m<=meqn; m++)
                        for (int k=1; k<=kmax; k++)
                        {
                            double tmp = -(q1.get(i,j,m,k)-qnew.get(i,j,m,k));
                            double err1= q1.get(i,j,m,k)+tmp+ IL.get(i,j,m,k,1);

                            qnew.set(i,j,m,k, err1 );
                        }   
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,qnew,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,qnew,*this);
            // --------------------------------------------------------

            break;

        case 3:  // 3rd order in time

            // --------------------------------------------------------
            sdc.num_iter = 1;
            // Take 1st Euler time step             
            sdc.num_stage = 1;
            StepSDC(dtvec.get(1),aux,qnew,L0,q1);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q1,*this);

            // Take 2nd Euler time step                  
            BeforeStep(dt,aux,q1,*this);
            SetBndValues(q1,aux);
            ::ConstructL(aux,q1,L1,smax);
            sdc.num_stage = 2;
            StepSDC(dtvec.get(2),aux,q1,L1,q2);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q2,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q2,*this);

            // Construct new right-hand side
            BeforeStep(dt,aux,q2,*this);
            SetBndValues(q2,aux);
            ::ConstructL(aux,q2,L2,smax);

            // Iterate to construct 
            for (int N=1; N<=(time_order-1); N++)
            {
                sdc.num_iter = N+1;

                // Integrate residual over time step
                ResInt(dt,L,IL);                   

                // First sub-interval
                sdc.num_stage = 1;
#pragma omp parallel for
                for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int j=(1-mbc); j<=(my+mbc); j++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {       
                                // Get error in first sub-interval
                                double err1 = -(q1.get(i,j,m,k)
                                        - qnew.get(i,j,m,k)) + IL.get(i,j,m,k,1);

                                // Correct solution
                                q1.set(i,j,m,k, q1.get(i,j,m,k) + err1 );
                            }
                if(dogParams.using_moment_limiter())
                { ::ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                AfterStep(dt,aux,q1,*this);

                // Construct new right-hand side
                BeforeStep(dt,aux,q1,*this);
                SetBndValues(q1,aux);
                ::ConstructL(aux,q1,L1new,smax);

                // Second sub-interval   
                sdc.num_stage = 2;
#pragma omp parallel for
                for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int j=(1-mbc); j<=(my+mbc); j++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {
                                // Get error in second sub-interval
                                double err2 = dtvec.get(2)*
                                    (L1new.get(i,j,m,k) - L1.get(i,j,m,k)) 
                                    - (q2.get(i,j,m,k)-q1.get(i,j,m,k)) 
                                    + IL.get(i,j,m,k,2);

                                // Correct solution
                                q2.set(i,j,m,k, q2.get(i,j,m,k) + err2 );
                                L1.set(i,j,m,k, L1new.get(i,j,m,k) );
                            }               
                if(dogParams.using_moment_limiter())
                { ::ApplyLimiter(aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                AfterStep(dt,aux,q2,*this);

                // Construct new right-hand side
                BeforeStep(dt,aux,q2,*this);
                SetBndValues(q2,aux);
                ::ConstructL(aux,q2,L2,smax);
            }

            // Update solution
            CopyQ(q2,qnew);
            // --------------------------------------------------------

            break;

        case 4:  // 4th order in time

            // --------------------------------------------------------
            // Take 1st Euler time step   
            sdc.num_iter = 1;
            sdc.num_stage = 1;
            StepSDC(dtvec.get(1),aux,qnew,L0,q1);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q1,*this);

            // Take 2nd Euler time step
            BeforeStep(dt,aux,q1,*this);
            SetBndValues(q1,aux);
            ::ConstructL(aux,q1,L1,smax);
            sdc.num_stage = 2;
            StepSDC(dtvec.get(2),aux,q1,L1,q2);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q2,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q2,*this);

            // Take 3rd Euler time step
            BeforeStep(dt,aux,q2,*this);
            SetBndValues(q2,aux);
            ::ConstructL(aux,q2,L2,smax);
            sdc.num_stage = 3;
            StepSDC(dtvec.get(3),aux,q2,L2,q3);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q3,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q3,*this);

            // Construct new right-hand side
            BeforeStep(dt,aux,q3,*this);
            SetBndValues(q3,aux);
            ::ConstructL(aux,q3,L3,smax);

            // Iterate to construct 
            for (int N=1; N<=(time_order-1); N++)
            {
                sdc.num_iter = N+1;

                // Integrate residual over time step
                ResInt(dt,L,IL);                   

                // First sub-interval     
                sdc.num_stage = 1;
#pragma omp parallel for
                for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int j=(1-mbc); j<=(my+mbc); j++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {       
                                // Get error in first sub-interval
                                double err1 = -(q1.get(i,j,m,k)
                                        - qnew.get(i,j,m,k)) + IL.get(i,j,m,k,1);

                                // Correct solution
                                q1.set(i,j,m,k, q1.get(i,j,m,k) + err1 );
                            }
                if(dogParams.using_moment_limiter())
                { ::ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                AfterStep(dt,aux,q1,*this);

                // Construct new right-hand side
                sdc.num_stage = 2;
                BeforeStep(dt,aux,q1,*this);
                SetBndValues(q1,aux);
                ::ConstructL(aux,q1,L1new,smax);

                // Second sub-interval                  
#pragma omp parallel for
                for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int j=(1-mbc); j<=(my+mbc); j++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {
                                // Get error in second sub-interval
                                double err2 = dtvec.get(2)*
                                    (L1new.get(i,j,m,k) - L1.get(i,j,m,k)) 
                                    - (q2.get(i,j,m,k)-q1.get(i,j,m,k)) 
                                    + IL.get(i,j,m,k,2);

                                // Correct solution
                                q2.set(i,j,m,k, q2.get(i,j,m,k) + err2 );
                                L1.set(i,j,m,k, L1new.get(i,j,m,k) );
                            }               
                if(dogParams.using_moment_limiter())
                { ::ApplyLimiter(aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                AfterStep(dt,aux,q2,*this);

                // Construct new right-hand side
                BeforeStep(dt,aux,q2,*this);
                SetBndValues(q2,aux);
                ::ConstructL(aux,q2,L2new,smax);

                // Third sub-interval
                sdc.num_stage = 3;
#pragma omp parallel for
                for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int j=(1-mbc); j<=(my+mbc); j++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {
                                // Get error in third sub-interval
                                double err3 = dtvec.get(3)*
                                    (L2new.get(i,j,m,k) - L2.get(i,j,m,k)) 
                                    - (q3.get(i,j,m,k)-q2.get(i,j,m,k)) 
                                    + IL.get(i,j,m,k,3);

                                // Correct solution
                                q3.set(i,j,m,k, q3.get(i,j,m,k) + err3 );
                                L2.set(i,j,m,k, L2new.get(i,j,m,k) );
                            }         
                if(dogParams.using_moment_limiter())
                { ::ApplyLimiter(aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                AfterStep(dt,aux,q3,*this);

                // Construct new right-hand side
                BeforeStep(dt,aux,q3,*this);
                SetBndValues(q3,aux);
                ::ConstructL(aux,q3,L3,smax);
            }  

            // Update solution
            CopyQ(q3,qnew);
            // --------------------------------------------------------

            break;

        case 5:  // 5th order in time

            // --------------------------------------------------------
            sdc.num_iter = 1;
            // Take 1st Euler time step             
            sdc.num_stage = 1;
            StepSDC(dtvec.get(1),aux,qnew,L0,q1);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q1,*this);

            // Take 2nd Euler time step
            BeforeStep(dt,aux,q1,*this);
            SetBndValues(q1,aux);
            ::ConstructL(aux,q1,L1,smax);
            sdc.num_stage = 2;
            StepSDC(dtvec.get(2),aux,q1,L1,q2);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q2,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q2,*this);

            // Take 3rd Euler time step
            BeforeStep(dt,aux,q2,*this);
            SetBndValues(q2,aux);
            ::ConstructL(aux,q2,L2,smax);                    
            sdc.num_stage = 3;
            StepSDC(dtvec.get(3),aux,q2,L2,q3);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q3,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q3,*this);

            // Take 4th Euler time step
            BeforeStep(dt,aux,q3,*this);
            SetBndValues(q3,aux);
            ::ConstructL(aux,q3,L3,smax);
            sdc.num_stage = 4;
            StepSDC(dtvec.get(4),aux,q3,L3,q4);
            if(dogParams.using_moment_limiter())
            { ::ApplyLimiter(aux,q4,&ProjectRightEig,&ProjectLeftEig); }
            AfterStep(dt,aux,q4,*this);

            // Construct new right-hand side
            BeforeStep(dt,aux,q4,*this);
            SetBndValues(q4,aux);
            ::ConstructL(aux,q4,L4,smax);

            // Iterate to construct 
            for (int N=1; N<=(time_order-1); N++)
            {
                sdc.num_iter = N+1;

                // Integrate residual over time step
                ResInt(dt,L,IL);                   

                // First sub-interval                   
                sdc.num_stage = 1;
#pragma omp parallel for
                for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int j=(1-mbc); j<=(my+mbc); j++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {       
                                // Get error in first sub-interval
                                double err1 = -(q1.get(i,j,m,k)
                                        - qnew.get(i,j,m,k)) + IL.get(i,j,m,k,1);

                                // Correct solution
                                q1.set(i,j,m,k, q1.get(i,j,m,k) + err1 );
                            }
                if(dogParams.using_moment_limiter())
                { ::ApplyLimiter(aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                AfterStep(dt,aux,q1,*this);

                // Construct new right-hand side
                BeforeStep(dt,aux,q1,*this);
                SetBndValues(q1,aux);
                ::ConstructL(aux,q1,L1new,smax);

                // Second sub-interval                  
                sdc.num_stage = 2;
#pragma omp parallel for
                for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int j=(1-mbc); j<=(my+mbc); j++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {
                                // Get error in second sub-interval
                                double err2 = dtvec.get(2)*
                                    (L1new.get(i,j,m,k) - L1.get(i,j,m,k)) 
                                    - (q2.get(i,j,m,k)-q1.get(i,j,m,k)) 
                                    + IL.get(i,j,m,k,2);

                                // Correct solution
                                q2.set(i,j,m,k, q2.get(i,j,m,k) + err2 );
                                L1.set(i,j,m,k, L1new.get(i,j,m,k) );
                            }
                if(dogParams.using_moment_limiter())
                { ::ApplyLimiter(aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                AfterStep(dt,aux,q2,*this);

                // Construct new right-hand side
                BeforeStep(dt,aux,q2,*this);
                SetBndValues(q2,aux);
                ::ConstructL(aux,q2,L2new,smax);

                // Third sub-interval                   
                sdc.num_stage = 3;
#pragma omp parallel for
                for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int j=(1-mbc); j<=(my+mbc); j++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {
                                // Get error in third sub-interval
                                double err3 = dtvec.get(3)*
                                    (L2new.get(i,j,m,k) - L2.get(i,j,m,k)) 
                                    - (q3.get(i,j,m,k)-q2.get(i,j,m,k)) 
                                    + IL.get(i,j,m,k,3);

                                // Correct solution
                                q3.set(i,j,m,k, q3.get(i,j,m,k) + err3 );
                                L2.set(i,j,m,k, L2new.get(i,j,m,k) );
                            } 
                if(dogParams.using_moment_limiter())
                { ::ApplyLimiter(aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                AfterStep(dt,aux,q3,*this);

                // Construct new right-hand side
                BeforeStep(dt,aux,q3,*this);
                SetBndValues(q3,aux);
                ::ConstructL(aux,q3,L3new,smax);

                // Fourth sub-interval
                sdc.num_stage = 4;
#pragma omp parallel for
                for (int i=(1-mbc); i<=(mx+mbc); i++)
                    for (int j=(1-mbc); j<=(my+mbc); j++)
                        for (int m=1; m<=meqn; m++)
                            for (int k=1; k<=kmax; k++)
                            {
                                // Get error in fourth sub-interval
                                double err4 = dtvec.get(4)*
                                    (L3new.get(i,j,m,k) - L3.get(i,j,m,k)) 
                                    - (q4.get(i,j,m,k)-q3.get(i,j,m,k)) 
                                    + IL.get(i,j,m,k,4);

                                // Correct solution
                                q4.set(i,j,m,k, q4.get(i,j,m,k) + err4 );
                                L3.set(i,j,m,k, L3new.get(i,j,m,k) );
                            }
                if(dogParams.using_moment_limiter())
                { ::ApplyLimiter(aux,q4,&ProjectRightEig,&ProjectLeftEig); }
                AfterStep(dt,aux,q4,*this);

                // Construct new right-hand side
                //
                BeforeStep(dt,aux,q4,*this);
                SetBndValues(q4,aux);
                ::ConstructL(aux,q4,L4,smax);
            }  

            // Update solution
            CopyQ(q4,qnew);
            // --------------------------------------------------------

            break;                  
    }
}
